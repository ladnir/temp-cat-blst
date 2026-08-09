#pragma once

#include "blst/common/defines.h"
#include "blst/curve/GT.h"
#include "blst/fields/Fr.h"
#include "blst/avx/fields/conversion_inversion.h"
#include "blst/avx/fields/matrixnok_fp12.h"
#include "blst/avx/gt/glv4.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>
#include <vector>

namespace blst::avx {

struct AvxGtElement
{
    using Ring = matrixnok::FpRing;
    using Fp12Ring = matrixnok::Fp12Ring;
    using Element = typename Fp12Ring::StandardElement;

    Element data = {};

    AvxGtElement() = default;
    explicit AvxGtElement(const Element& value) : data(value) {}

    inline static constexpr Ring kRing{};
    inline static constexpr Fp12Ring kFp12Ring{};

    static const Ring& ring()
    {
        return kRing;
    }

    static const Fp12Ring& fp12ring()
    {
        return kFp12Ring;
    }

    static AvxGtElement identity()
    {
        return AvxGtElement(fp12ring().one(ring()));
    }

    static AvxGtElement from_scalar(const blst::scalar::GT& value)
    {
        return AvxGtElement(avx_blst_direct_bridge::scalar_fp12_to_ring(ring(), value.data));
    }

    static AvxGtElement from_default(const blst::GT& value)
    {
#if BLST_DEFAULT_AVX
        return AvxGtElement(value.data);
#else
        return from_scalar(value);
#endif
    }

    static blst::GT to_default(const AvxGtElement& value)
    {
#if BLST_DEFAULT_AVX
        blst::avx::GT out = {};
        out.data = value.data;
        return out;
#else
        blst::scalar::GT out = {};
        out.data = avx_blst_direct_bridge::ring_fp12_to_scalar(ring(), value.data);
        return out;
#endif
    }
};

inline std::array<uint64_t, 4> fr_to_limbs_le(const Fr& value)
{
    Fr canonical = {};
    Fr::from_mont(canonical, value);
    return {
        canonical.l[0],
        canonical.l[1],
        canonical.l[2],
        canonical.l[3],
    };
}

struct AvxGtPreparedRoots
{
    std::vector<Fr> scalars;
    std::vector<std::array<uint64_t, 4>> limbs;
    std::vector<AvxGLV4Exp> compiled;
    Fr scale_scalar = {};
    std::array<uint64_t, 4> scale_limbs = {};
    AvxGLV4Exp scale_compiled = {};
    bool has_scale = false;

    explicit AvxGtPreparedRoots(std::span<const Fr> input_scalars = {})
    {
        scalars.assign(input_scalars.begin(), input_scalars.end());
        limbs.reserve(scalars.size());
        compiled.reserve(scalars.size());
        for (const Fr& scalar : scalars) {
            limbs.push_back(fr_to_limbs_le(scalar));
            compiled.push_back(avx_glv4_compile(limbs.back().data()));
        }
    }

    std::span<const Fr> scalar_span() const
    {
        return scalars;
    }

    std::size_t size() const
    {
        return scalars.size();
    }

    void set_scale_from_size(std::size_t n)
    {
        BLST_ASSERT(n != 0);
        Fr::from_uint64(scale_scalar, static_cast<uint64_t>(n));
        Fr::reciprocal(scale_scalar, scale_scalar);
        scale_limbs = fr_to_limbs_le(scale_scalar);
        scale_compiled = avx_glv4_compile(scale_limbs.data());
        has_scale = true;
    }
};

struct AvxGtCoeffCtx
{
    class ScopedRoots
    {
    public:
        explicit ScopedRoots(const AvxGtPreparedRoots& roots) noexcept
            : previous_(prepared_roots_)
        {
            prepared_roots_ = &roots;
        }

        ScopedRoots(const AvxGtPreparedRoots&&) = delete;

        ~ScopedRoots() noexcept
        {
            prepared_roots_ = previous_;
        }

        ScopedRoots(const ScopedRoots&) = delete;
        ScopedRoots& operator=(const ScopedRoots&) = delete;
        ScopedRoots(ScopedRoots&&) = delete;
        ScopedRoots& operator=(ScopedRoots&&) = delete;

    private:
        const AvxGtPreparedRoots* previous_;
    };

    static auto plus(const AvxGtElement& a, const AvxGtElement& b)
    {
        return AvxGtElement(AvxGtElement::fp12ring().fp12_flat_mul(a.data, b.data, AvxGtElement::ring()));
    }

    static auto minus(const AvxGtElement& a, const AvxGtElement& b)
    {
        const auto b_inv = AvxGtElement::fp12ring().fp12_flat_conjugate(b.data, AvxGtElement::ring());
        return AvxGtElement(AvxGtElement::fp12ring().fp12_flat_mul(a.data, b_inv, AvxGtElement::ring()));
    }

    static auto mult(const AvxGtElement& a, const Fr& b)
    {
        const auto& ring = AvxGtElement::ring();
        const auto& fp12ring = AvxGtElement::fp12ring();

        if (prepared_roots_ != nullptr && !prepared_roots_->scalars.empty()) {
            const Fr* base = prepared_roots_->scalars.data();
            const Fr* const end = base + prepared_roots_->scalars.size();
            const Fr* const candidate = &b;
            const std::less<const Fr*> less;
            if (!less(candidate, base) && less(candidate, end)) {
                const std::uintptr_t base_address =
                    reinterpret_cast<std::uintptr_t>(base);
                const std::uintptr_t candidate_address =
                    reinterpret_cast<std::uintptr_t>(candidate);
                if (candidate_address >= base_address &&
                    (candidate_address - base_address) % sizeof(Fr) == 0) {
                    const std::uintptr_t byte_offset =
                        candidate_address - base_address;
                    const std::size_t idx =
                        static_cast<std::size_t>(byte_offset / sizeof(Fr));
                    if (idx < prepared_roots_->compiled.size() &&
                        candidate == &prepared_roots_->scalars[idx]) {
                        return AvxGtElement(avx_cyclotomic_exp_glv4(
                            a.data,
                            prepared_roots_->compiled[idx],
                            fp12ring,
                            ring));
                    }
                }
            }
        }

        if (prepared_roots_ != nullptr &&
            prepared_roots_->has_scale &&
            Fr::equal(b, prepared_roots_->scale_scalar)) {
            return AvxGtElement(avx_cyclotomic_exp_glv4(
                a.data,
                prepared_roots_->scale_compiled,
                fp12ring,
                ring));
        }

        const std::array<uint64_t, 4> scalar_limbs = fr_to_limbs_le(b);
        return AvxGtElement(avx_cyclotomic_exp_glv4(a.data, scalar_limbs.data(), fp12ring, ring));
    }

    static auto inverse(const Fr& a)
    {
        Fr r = {};
        Fr::reciprocal(r, a);
        return r;
    }

    template<typename F>
    static auto from_u64(uint64_t a)
    {
        Fr r = {};
        Fr::from_uint64(r, a);
        return r;
    }

private:
    inline static thread_local const AvxGtPreparedRoots* prepared_roots_ = nullptr;
};

inline std::vector<AvxGtElement> scalar_gt_vector_to_avx(std::span<const blst::GT> values)
{
    std::vector<AvxGtElement> out;
    out.reserve(values.size());
    for (const blst::GT& value : values) {
        out.push_back(AvxGtElement::from_default(value));
    }
    return out;
}

inline std::vector<blst::GT> avx_vector_to_scalar_gt(std::span<const AvxGtElement> values)
{
    std::vector<blst::GT> out;
    out.reserve(values.size());
    for (const AvxGtElement& value : values) {
        out.push_back(AvxGtElement::to_default(value));
    }
    return out;
}

}  // namespace blst::avx
