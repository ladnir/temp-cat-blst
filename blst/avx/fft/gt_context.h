// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// CatBLST AVX GT FFT adapter, locally composed around the MatrixNoK field
// engine, the CatBLST FFT coefficient interface, and the proved GLV4 executor.

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

/**
 * MatrixNoK storage for one element of the prime-order target group GT.
 *
 * `data` is a standard MatrixNoK Fp12 representation, which may be redundant
 * but must represent a cyclotomic GT element. The default all-zero state is an
 * overwrite-before-read placeholder, not the GT identity. No operation here
 * validates GT membership. In particular, `AvxGtCoeffCtx::minus` uses
 * cyclotomic conjugation as inversion and is not an arbitrary-Fp12 operation.
 */
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

    static const Ring& ring() noexcept
    {
        return kRing;
    }

    static const Fp12Ring& fp12ring() noexcept
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

/** Convert a reduced Montgomery Fr value to four canonical little-endian limbs. */
inline std::array<uint64_t, 4> fr_to_limbs_le(const Fr& value) noexcept
{
    BLST_ASSERT(Fr::is_redc(value));
    Fr canonical = {};
    Fr::from_mont(canonical, value);
    BLST_ASSERT(Fr::is_redc(canonical));
    return {
        canonical.l[0],
        canonical.l[1],
        canonical.l[2],
        canonical.l[3],
    };
}

/**
 * Immutable prepared public Fr exponents for AVX GT FFT multiplication.
 *
 * Construction allocates, copies the reduced Montgomery scalars, and compiles
 * one variable-time GLV4 schedule per scalar. Once a `ScopedRoots` guard can
 * observe this object, its public fields must not be mutated, moved, or
 * destroyed until every guard has left. Immutable objects may be shared by
 * concurrent threads because selection state is thread-local.
 */
struct AvxGtPreparedRoots
{
    std::vector<Fr> scalars;
    std::vector<AvxGLV4Exp> compiled;
    Fr scale_scalar = {};
    AvxGLV4Exp scale_compiled = {};
    bool has_scale = false;

    explicit AvxGtPreparedRoots(std::span<const Fr> input_scalars = {})
        : scalars(input_scalars.begin(), input_scalars.end()),
          compiled(scalars.size())
    {
        for (std::size_t i = 0; i < scalars.size(); ++i) {
            const std::array<uint64_t, 4> limbs = fr_to_limbs_le(scalars[i]);
            compiled[i] = avx_glv4_compile(limbs.data());
        }
    }

    std::span<const Fr> scalar_span() const noexcept
    {
        return scalars;
    }

    std::size_t size() const noexcept
    {
        return scalars.size();
    }

    /** Compile the inverse-FFT scale `1 / size()` for this nonempty table. */
    void prepare_inverse_scale() noexcept
    {
        BLST_ASSERT(!scalars.empty());
        BLST_ASSERT(scalars.size() == compiled.size());
        const std::size_t n = scalars.size();
        Fr::from_uint64(scale_scalar, static_cast<uint64_t>(n));
        Fr::reciprocal(scale_scalar, scale_scalar);
        const std::array<uint64_t, 4> limbs = fr_to_limbs_le(scale_scalar);
        scale_compiled = avx_glv4_compile(limbs.data());
        has_scale = true;
    }
};

/**
 * Variable-time GT coefficient policy for public Fr FFT roots and scales.
 *
 * `plus`, `minus`, and `mult` require GT operands and perform no membership
 * validation or allocation. `mult` selects a precompiled schedule only when
 * its Fr reference is exactly an element of the active prepared table, or its
 * value equals the prepared inverse scale; otherwise it compiles a temporary
 * schedule on the stack. All scalar-dependent branches and work therefore
 * depend only on public FFT coefficients.
 *
 * A `ScopedRoots` object must remain on one thread during synchronous use, and
 * its referenced table must outlive it without mutation or movement. Guards
 * may nest; destruction restores the previous table even during unwinding.
 */
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

    static AvxGtElement plus(const AvxGtElement& a, const AvxGtElement& b)
    {
        return AvxGtElement(AvxGtElement::fp12ring().fp12_flat_mul(a.data, b.data, AvxGtElement::ring()));
    }

    static AvxGtElement minus(const AvxGtElement& a, const AvxGtElement& b)
    {
        const auto b_inv = AvxGtElement::fp12ring().fp12_flat_conjugate(b.data, AvxGtElement::ring());
        return AvxGtElement(AvxGtElement::fp12ring().fp12_flat_mul(a.data, b_inv, AvxGtElement::ring()));
    }

    static AvxGtElement mult(const AvxGtElement& a, const Fr& b)
    {
        const auto& ring = AvxGtElement::ring();
        const auto& fp12ring = AvxGtElement::fp12ring();

        if (prepared_roots_ != nullptr && !prepared_roots_->scalars.empty()) {
            BLST_ASSERT(
                prepared_roots_->scalars.size() ==
                prepared_roots_->compiled.size());
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

    static Fr inverse(const Fr& a)
    {
        Fr r = {};
        Fr::reciprocal(r, a);
        return r;
    }

    template<typename F>
    static Fr from_u64(uint64_t a)
    {
        Fr r = {};
        Fr::from_uint64(r, a);
        return r;
    }

private:
    inline static thread_local const AvxGtPreparedRoots* prepared_roots_ = nullptr;
};

}  // namespace blst::avx
