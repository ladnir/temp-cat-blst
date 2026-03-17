#ifndef BLST_FIELDS_FP12_AVX_H
#define BLST_FIELDS_FP12_AVX_H

/*
 * WARNING: This file is intentionally excluded from the active build during the
 * VROOM backend reset. It is retained only as reference and must not be wired
 * back into active code until the post-reset AVX API is defined.
 */

#include "blst/fields/Fp6_avx.h"
#include "blst/fields/Fp12.h"
#include <array>

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    struct UFp12Pack8x1;

    /**
     * AVX512 Fp12 helpers ported from avxbls.
     *
     * Field model:
     *   Fp12 = Fp6[w] / (w^2 - v)
     *   Fp6  = Fp2[v] / (v^3 - (u + 1))
     *   Fp2  = Fp[u] / (u^2 + 1)
     *
     * All packed values are Montgomery-form field elements.
     */

    /**
     * Canonical W=8 Fp12 pack used by the batch tower path.
     *
     * Component order:
     *   c[0] + w*c[1]
     * where each component is a canonical `Fp6Pack8x1`.
     */
    struct Fp12Pack8x1
    {
        std::array<Fp6Pack8x1, 2> c;

        BLST_FORCEINLINE Fp6Pack8x1& operator[](std::size_t i)
        {
            BLST_ASSERT(i < c.size());
            return c[i];
        }

        BLST_FORCEINLINE const Fp6Pack8x1& operator[](std::size_t i) const
        {
            BLST_ASSERT(i < c.size());
            return c[i];
        }

        BLST_FORCEINLINE Fp6Pack8x1* data()
        {
            return c.data();
        }

        BLST_FORCEINLINE const Fp6Pack8x1* data() const
        {
            return c.data();
        }

        /// Pack 8 scalar `Fp12` values (AoS) into one canonical `Fp12Pack8x1`.
        static void from_aos(Fp12Pack8x1& r, const Fp12* in8);

        /// Unpack one canonical `Fp12Pack8x1` into 8 scalar `Fp12` values (AoS).
        static void to_aos(Fp12* out8, const Fp12Pack8x1& a);

        /// Partial AoS pack for tail handling.
        static void from_aos_partial(Fp12Pack8x1& r, const Fp12* in, std::size_t n_active, const Fp12& fill);

        /// Partial AoS unpack for tail handling.
        static void to_aos_partial(Fp12* out, std::size_t n_active, const Fp12Pack8x1& a);

        /// Component-wise modular add: r = a + b.
        static void add(Fp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b);

        /// Component-wise modular sub: r = a - b.
        static void sub(Fp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b);

        /// Reduced Fp12 multiply.
        static void mul(Fp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b);
        /// Widened Fp12 multiply.
        static void mul(UFp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b);

        /// Reduced Fp12 square.
        static void sqr(Fp12Pack8x1& r, const Fp12Pack8x1& a);
        /// Widened Fp12 square.
        static void sqr(UFp12Pack8x1& r, const Fp12Pack8x1& a);

        /// Cyclotomic square API entry; functionally equal to square.
        static void cyclotomic_sqr(Fp12Pack8x1& r, const Fp12Pack8x1& a);

        /// Frobenius map: r = a^(p^n), n reduced mod 4.
        static void frobenius_map(Fp12Pack8x1& r, const Fp12Pack8x1& a, std::size_t n);

        /// Conjugation: (c0, c1) -> (c0, -c1).
        static void conjugate(Fp12Pack8x1& r, const Fp12Pack8x1& a);

        /// Lane-wise select between two packed values.
        /// Bit `i` of `lane_mask` selects lane `i`: `0 -> a`, `1 -> b`.
        static void select(Fp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b, u8 lane_mask);
    };

    /**
     * Widened companion for canonical `Fp12Pack8x1`.
     */
    struct UFp12Pack8x1
    {
        std::array<UFp6Pack8x1, 2> c;

        BLST_FORCEINLINE UFp6Pack8x1& operator[](std::size_t i)
        {
            BLST_ASSERT(i < c.size());
            return c[i];
        }

        BLST_FORCEINLINE const UFp6Pack8x1& operator[](std::size_t i) const
        {
            BLST_ASSERT(i < c.size());
            return c[i];
        }

        BLST_FORCEINLINE UFp6Pack8x1* data()
        {
            return c.data();
        }

        BLST_FORCEINLINE const UFp6Pack8x1* data() const
        {
            return c.data();
        }

        /// Reduce widened result into reduced `Fp12Pack8x1`.
        static void redc(Fp12Pack8x1& r, const UFp12Pack8x1& a);

        /// Lane-wise select between two widened packed values.
        /// Bit `i` of `lane_mask` selects lane `i`: `0 -> a`, `1 -> b`.
        static void select(UFp12Pack8x1& r, const UFp12Pack8x1& a, const UFp12Pack8x1& b, u8 lane_mask);
    };

    /**
     * Packed schedule of three Fp2 blocks for Fp12 multiplication inputs in `8x1x1w`.
     *
     * Algorithmic role:
     *   packs[0] = ab0
     *   packs[1] = ab1
     *   packs[2] = ab2
     * where `abk` is the pre-arranged product term lane bundle consumed by Fp12 kernels.
     *
     * Visual topology:
     *   [ packs[2] | packs[1] | packs[0] ]
     *   [   ab2    |   ab1    |   ab0    ]
     */
    struct Fp12MulInput8x1x1 : public AvxPackTuple<Fp2Pack8x1x1, 3>
    {
        /// Schoolbook Fp12 multiply variant 1, output as `r01` + `r2` in `4x2x1w`.
        static void mul_v1(
            Fp2Pack4x2x1& r01,
            Fp2Pack4x2x1& r2,
            const Fp12MulInput8x1x1& ab);

        /// Schoolbook Fp12 multiply variant 2, output as `r01` in `4x2x1w` + `r2` in `2x2x2w`.
        static void mul_v2(
            Fp2Pack4x2x1& r01,
            Fp2Pack2x2x2& r2,
            const Fp12MulInput8x1x1& ab);
    };

    /**
     * Packed schedule of three Fp2 blocks for Fp12 multiplication inputs in `4x2x1w`.
     *
     * Algorithmic role:
     *   packs[0] = ab0
     *   packs[1] = ab1
     *   packs[2] = ab2
     */
    struct Fp12MulInput4x2x1 : public AvxPackTuple<Fp2Pack4x2x1, 3>
    {
        /// Karatsuba Fp12 multiply variant 3 (`2x2x2w`-based split output).
        static void mul_v3(
            Fp2Pack2x2x2& r001,
            Fp2Pack2x2x2& r101,
            Fp2Pack2x2x2& r2,
            const Fp12MulInput4x2x1& ab);

        /// Karatsuba Fp12 multiply variant 4 (`4x2x1w` + `2x2x2w` split output).
        static void mul_v4(
            Fp2Pack4x2x1& r0,
            Fp2Pack2x2x2& r101,
            Fp2Pack2x2x2& r12,
            const Fp12MulInput4x2x1& ab);
    };

    /**
     * Sparse Fp12 multiply input in `4x2x1w` schedule for `(x,y,0,0,z,0)`.
     *
     * Field component mapping:
     *   a01 = (a1 | a0)
     *   a2  = (a2 | a2)
     *   b01 = (b1 | b0)
     *   b4  = (b4 | b4)
     */
    struct Fp12SparseInput4x2x1
    {
        Fp2Pack4x2x1 a01;
        Fp2Pack4x2x1 a2;
        Fp2Pack4x2x1 b01;
        Fp2Pack4x2x1 b4;

        /// Sparse Fp12 multiply variant 1.
        static void mul_v1(
            Fp2Pack4x2x1& r0,
            Fp2Pack2x2x2& r1,
            const Fp12SparseInput4x2x1& in);
    };

    /**
     * Sparse Fp12 multiply input in `8x1x1w` schedule for `(x,y,0,0,z,0)`.
     *
     * Field component mapping:
     *   a = packed source block for (a2,a1,a0) schedule lanes
     *   b = packed source block for sparse multiplier lanes (b0,b1,b4)
     */
    struct Fp12SparseInput8x1x1
    {
        Fp2Pack8x1x1 a;
        Fp2Pack8x1x1 b;

        /// Sparse Fp12 multiply variant 2.
        static void mul_v2(
            Fp2Pack4x2x1& r0,
            Fp2Pack2x2x2& r1,
            const Fp12SparseInput8x1x1& in);
    };

    /**
     * Fp12 square input schedule in `4x2x1w`.
     *
     * Layout model:
     *   a0 = low half schedule lanes
     *   a1 = high half schedule lanes
     */
    struct Fp12SquareInput4x2x1
    {
        Fp2Pack4x2x1 a0;
        Fp2Pack4x2x1 a1;

        /// Karatsuba-style Fp12 square variant 1.
        static void sqr_v1(
            Fp2Pack4x2x1& r0,
            Fp2Pack4x2x1& r1,
            const Fp12SquareInput4x2x1& in);
    };

    /**
     * Cyclotomic B/C schedule block for Fp12.
     *
     * Visual mapping (Fp2 slots):
     *   [z5 | z4 | z3 | z2] in AVXBLS `2x2x2x1w` block order.
     */
    struct Fp12CyclotomicBc : public Fp4Pack2x2x2x1
    {
        /// Compressed cyclotomic square update of B/C block only.
        static void compressed_sqr_v1(Fp12CyclotomicBc& rbc, const Fp12CyclotomicBc& bc);
    };

    /**
     * Cyclotomic Fp12 state schedule.
     *
     * Component mapping:
     *   a  = A block (`z0,z1`) in `1x2x2x2w`
     *   bc = B/C block (`z2..z5`) in `2x2x2x1w`
     */
    struct Fp12CyclotomicState
    {
        Fp4Pack1x2x2x2 a;
        Fp12CyclotomicBc bc;

        /// Cyclotomic square variant 1 (double-length Fp4 sub-kernels).
        static void sqr_v1(Fp12CyclotomicState& r, const Fp12CyclotomicState& in);

        /// Cyclotomic square variant 2 (single-length Fp4 sub-kernels).
        static void sqr_v2(Fp12CyclotomicState& r, const Fp12CyclotomicState& in);
    };

    static_assert(std::is_trivially_default_constructible_v<Fp12MulInput8x1x1>);
    static_assert(std::is_trivially_copyable_v<Fp12MulInput8x1x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp12Pack8x1>);
    static_assert(std::is_trivially_copyable_v<Fp12Pack8x1>);
    static_assert(std::is_trivially_default_constructible_v<UFp12Pack8x1>);
    static_assert(std::is_trivially_copyable_v<UFp12Pack8x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp12MulInput4x2x1>);
    static_assert(std::is_trivially_copyable_v<Fp12MulInput4x2x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp12SparseInput4x2x1>);
    static_assert(std::is_trivially_copyable_v<Fp12SparseInput4x2x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp12SparseInput8x1x1>);
    static_assert(std::is_trivially_copyable_v<Fp12SparseInput8x1x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp12SquareInput4x2x1>);
    static_assert(std::is_trivially_copyable_v<Fp12SquareInput4x2x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp12CyclotomicBc>);
    static_assert(std::is_trivially_copyable_v<Fp12CyclotomicBc>);
    static_assert(std::is_trivially_default_constructible_v<Fp12CyclotomicState>);
    static_assert(std::is_trivially_copyable_v<Fp12CyclotomicState>);

    /**
     * AVXBLS wrapper: full Fp12 multiply with scalar-layout inputs/outputs.
     */
    void mul_fp12_vector(Fp12& ret, const Fp12& a, const Fp12& b);

    /**
     * AVXBLS wrapper: sparse Fp12 multiply by `(x,y,0,0,z,0)`.
     */
    void mul_by_xy00z0_fp12_vector(Fp12& ret, const Fp12& a, const Fp6& xy00z0);

    /**
     * AVXBLS wrapper: full Fp12 squaring.
     */
    void sqr_fp12_vector(Fp12& ret, const Fp12& a);

    /**
     * AVXBLS wrapper: cyclotomic square with scalar-layout inputs/outputs.
     */
    void cyclotomic_sqr_fp12_vector(Fp12& ret, const Fp12& a);
}

#endif

#endif // BLST_FIELDS_FP12_AVX_H
