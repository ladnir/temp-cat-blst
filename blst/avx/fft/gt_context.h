#pragma once

#include "blst/common/defines.h"
#include "blst/curve/GT.h"
#include "blst/fields/Fr.h"
#include "blst/avx/rns/bounded_ring.h"
#include "blst/avx/fields/conversion_inversion.h"
#include "blst/avx/fields/fp12.h"
#include "blst/avx/gt/glv4.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace blst::avx {

inline constexpr char kAvxGt52BaseModulus[] =
    "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";

struct AvxGt52
{
    using Ring = BoundedRing<>;
    using Fp12Ring = FP12Ring<Ring>;
    using Element = typename Fp12Ring::StandardElement;

    Element data = {};

    AvxGt52() = default;
    explicit AvxGt52(const Element& value) : data(value) {}

    static const Ring& ring()
    {
        static const Ring instance(BigInt(kAvxGt52BaseModulus, 10));
        return instance;
    }

    static const Fp12Ring& fp12ring()
    {
        static const Fp12Ring instance(ring());
        return instance;
    }

    static AvxGt52 identity()
    {
        return AvxGt52(fp12ring().one(ring()));
    }

    static AvxGt52 from_scalar(const blst::scalar::GT& value)
    {
        return AvxGt52(avx_blst_direct_bridge::scalar_fp12_to_ring(ring(), value.data));
    }

    static AvxGt52 from_default(const blst::GT& value)
    {
#if BLST_DEFAULT_AVX
        return AvxGt52(value.data);
#else
        return from_scalar(value);
#endif
    }

    static blst::GT to_default(const AvxGt52& value)
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
        canonical.data[0],
        canonical.data[1],
        canonical.data[2],
        canonical.data[3],
    };
}

struct AvxGt52PreparedRoots
{
    std::vector<Fr> scalars;
    std::vector<std::array<uint64_t, 4>> limbs;
    std::vector<AvxGLV4Exp> compiled;
    Fr scale_scalar = {};
    std::array<uint64_t, 4> scale_limbs = {};
    AvxGLV4Exp scale_compiled = {};
    bool has_scale = false;

    explicit AvxGt52PreparedRoots(std::span<const Fr> input_scalars = {})
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

struct AvxGt52CoeffCtx
{
    class ScopedRoots
    {
    public:
        explicit ScopedRoots(const AvxGt52PreparedRoots& roots)
            : previous_(prepared_roots_)
        {
            prepared_roots_ = &roots;
        }

        ~ScopedRoots()
        {
            prepared_roots_ = previous_;
        }

        ScopedRoots(const ScopedRoots&) = delete;
        ScopedRoots& operator=(const ScopedRoots&) = delete;

    private:
        const AvxGt52PreparedRoots* previous_;
    };

    static auto plus(const AvxGt52& a, const AvxGt52& b)
    {
        return AvxGt52(AvxGt52::fp12ring().fp12_flat_mul(a.data, b.data, AvxGt52::ring()));
    }

    static auto minus(const AvxGt52& a, const AvxGt52& b)
    {
        const auto b_inv = AvxGt52::fp12ring().fp12_flat_conjugate(b.data, AvxGt52::ring());
        return AvxGt52(AvxGt52::fp12ring().fp12_flat_mul(a.data, b_inv, AvxGt52::ring()));
    }

    static auto mult(const AvxGt52& a, const Fr& b)
    {
        const auto& ring = AvxGt52::ring();
        const auto& fp12ring = AvxGt52::fp12ring();

        if (prepared_roots_ != nullptr && !prepared_roots_->scalars.empty()) {
            const Fr* base = prepared_roots_->scalars.data();
            if (&b >= base && &b < base + prepared_roots_->scalars.size()) {
                const std::size_t idx = static_cast<std::size_t>(&b - base);
                return AvxGt52(avx_cyclotomic_exp_glv4(
                    a.data,
                    prepared_roots_->compiled[idx],
                    fp12ring,
                    ring));
            }
        }

        if (prepared_roots_ != nullptr &&
            prepared_roots_->has_scale &&
            Fr::equal(b, prepared_roots_->scale_scalar)) {
            return AvxGt52(avx_cyclotomic_exp_glv4(
                a.data,
                prepared_roots_->scale_compiled,
                fp12ring,
                ring));
        }

        const std::array<uint64_t, 4> scalar_limbs = fr_to_limbs_le(b);
        return AvxGt52(avx_cyclotomic_exp_glv4(a.data, scalar_limbs.data(), fp12ring, ring));
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
    inline static thread_local const AvxGt52PreparedRoots* prepared_roots_ = nullptr;
};

inline std::vector<AvxGt52> scalar_gt_vector_to_avx52(std::span<const blst::GT> values)
{
    std::vector<AvxGt52> out;
    out.reserve(values.size());
    for (const blst::GT& value : values) {
        out.push_back(AvxGt52::from_default(value));
    }
    return out;
}

inline std::vector<blst::GT> avx52_vector_to_scalar_gt(std::span<const AvxGt52> values)
{
    std::vector<blst::GT> out;
    out.reserve(values.size());
    for (const AvxGt52& value : values) {
        out.push_back(AvxGt52::to_default(value));
    }
    return out;
}

}  // namespace blst::avx
