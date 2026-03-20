#ifndef BLST_FIELDS_FP2_AVX_H
#define BLST_FIELDS_FP2_AVX_H

#include "blst/fields/Fp2.h"
#include "blst/fields/Fp_avx.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    /**
     * Lightweight aggregate for fixed-size bundles of homogeneous AVX packs.
     *
     * Visual model:
     *   packs[0], packs[1], ..., packs[N-1]
     *
     * Invariants:
     * - trivial type: no constructors, no ownership, no zero-init side effects
     * - indexing is semantic slot selection (not a layout transform)
     * - debug-only bounds checks via BLST_ASSERT
     */
    template <class PackT, std::size_t N>
    struct AvxPackTuple
    {
        PackT packs[N];

        BLST_FORCEINLINE PackT& operator[](std::size_t i)
        {
            BLST_ASSERT(i < N);
            return packs[i];
        }

        BLST_FORCEINLINE const PackT& operator[](std::size_t i) const
        {
            BLST_ASSERT(i < N);
            return packs[i];
        }

        BLST_FORCEINLINE PackT* data()
        {
            return packs;
        }

        BLST_FORCEINLINE const PackT* data() const
        {
            return packs;
        }
    };

    /**
     * Fp2 schedule 8x1x1: structure-of-arrays over 8 Fp2 elements.
     *
     * Field model used by BLST:
     *   Fp2 = Fp[u] / (u^2 + 1), i.e. u^2 = -1.
     *
     * Representation model (tuple view):
     *   x_i := (x_i0, x_i1)  (equivalently x_i0 + u*x_i1), for i in [0..7].
     *   packs[0] stores all x_i0 terms.
     *   packs[1] stores all x_i1 terms.
     *
     * Visual (per semantic lane i):
     *   packs[0]: [x_70 x_60 x_50 x_40 x_30 x_20 x_10 x_00]
     *   packs[1]: [x_71 x_61 x_51 x_41 x_31 x_21 x_11 x_01]
     */
    struct Fp2Pack8x1x1 : public AvxPackTuple<FpPack8x1, 2>
    {
        /// Pack 8 scalar `Fp2` values (AoS) into one `Fp2Pack8x1x1`.
        static void from_aos(Fp2Pack8x1x1& r, const Fp2* in8);

        /// Unpack one `Fp2Pack8x1x1` into 8 scalar `Fp2` values (AoS).
        static void to_aos(Fp2* out8, const Fp2Pack8x1x1& a);

        /// Partial AoS pack for tail handling.
        static void from_aos_partial(Fp2Pack8x1x1& r, const Fp2* in, std::size_t n_active, const Fp2& fill);

        /// Partial AoS unpack for tail handling.
        static void to_aos_partial(Fp2* out, std::size_t n_active, const Fp2Pack8x1x1& a);

        /// Component-wise modular add: r_i = a_i + b_i.
        static void add(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b);

        /// Component-wise modular sub: r_i = a_i - b_i.
        static void sub(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b);

        /// Reduced Fp2 multiply.
        static void mul(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b);

        /// Reduced Fp2 square.
        static void sqr(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a);

        /// Conjugation: (a0, a1) -> (a0, -a1).
        static void conjugate(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a);

        /// Multiply by `(u+1)`: (a0, a1) -> (a0-a1, a0+a1).
        static void mul_by_u_plus_1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a);

        /// Lane-wise select between two packed values.
        /// Bit `i` of `lane_mask` selects lane `i`: `0 -> a`, `1 -> b`.
        static void select(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b, u8 lane_mask);
    };

    /**
     * Widened companion for 8x1x1 (same lane topology, unreduced coefficients).
     *
     * Each component is UFpPack8x1 and holds double-width partials between mul/redc.
     * Kernel note for mul():
     *   t0 = a0*b0
     *   t1 = a1*b1
     *   t2 = (a0+a1)*(b0+b1)
     *   r0 = t0 - t1
     *   r1 = t2 - t0 - t1
     */
    struct UFp2Pack8x1x1 : public AvxPackTuple<UFpPack8x1, 2>
    {
        /// Unreduced component-wise add: r_i = a_i + b_i.
        static void add(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b);
        /// Unreduced component-wise sub: r_i = a_i - b_i.
        static void sub(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b);
        /// Mixed op on each component vector: r0=a0+b0, r1=a1-b1.
        static void asx4(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b);
        /// Fp2 multiply by (u+1): (x0 + u*x1) -> (x0-x1) + u*(x0+x1).
        static void mul_by_u_plus_1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a);
        /// Karatsuba Fp2 multiply of reduced inputs, unreduced output.
        static void mul(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b);
        /// Fp2 square of reduced input, unreduced output.
        static void sqr(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a);

        /// Reduce widened result into reduced `Fp2Pack8x1x1`.
        static void redc(Fp2Pack8x1x1& r, const UFp2Pack8x1x1& a);

        /// Lane-wise select between two widened packed values.
        /// Bit `i` of `lane_mask` selects lane `i`: `0 -> a`, `1 -> b`.
        static void select(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b, u8 lane_mask);
    };

    /**
     * Fp2 schedule 4x2x1 packed into one FpPack8x1.
     *
     * Let:
     *   A = A0 + u*A1
     *   B = B0 + u*B1
     *   C = C0 + u*C1
     *   D = D0 + u*D1
     *
     * Per-word lane layout (high lane -> low lane):
     *   [D1 | D0 | C1 | C0 | B1 | B0 | A1 | A0]
     *
     * Kernel model:
     * - pair lanes form one Fp2 value
     * - mul/sqr use lane shuffles to build schoolbook intermediates, then redc
     */
    struct Fp2Pack4x2x1 : public FpPack8x1
    {
        /// Per packed Fp2 pair: (x0,x1) + (y0,y1).
        static void add(Fp2Pack4x2x1& r, const Fp2Pack4x2x1& a, const Fp2Pack4x2x1& b);
        /// Per packed Fp2 pair: (x0,x1) - (y0,y1).
        static void sub(Fp2Pack4x2x1& r, const Fp2Pack4x2x1& a, const Fp2Pack4x2x1& b);
        /// Pattern op by packed element [D,C,B,A]: add, add, sub, sub.
        static void aass(Fp2Pack4x2x1& r, const Fp2Pack4x2x1& a, const Fp2Pack4x2x1& b);
        /// Pattern op by packed element [D,C,B,A]: add, sub, sub, add.
        static void assa(Fp2Pack4x2x1& r, const Fp2Pack4x2x1& a, const Fp2Pack4x2x1& b);
        /// Pattern op by packed element [D,C,B,A]: sub, add, sub, sub.
        static void sass(Fp2Pack4x2x1& r, const Fp2Pack4x2x1& a, const Fp2Pack4x2x1& b);
        /// Pattern op by packed element [D,C,B,A]: sub, sub, add, add.
        static void ssaa(Fp2Pack4x2x1& r, const Fp2Pack4x2x1& a, const Fp2Pack4x2x1& b);
        /// Per packed Fp2, multiply by (u+1): (x0,x1)->(x0-x1, x0+x1).
        static void mul_by_u_plus_1(Fp2Pack4x2x1& r, const Fp2Pack4x2x1& a);
        /// Schoolbook Fp2 mul per pair: (x0*y0-x1*y1, x0*y1+x1*y0), then redc.
        static void mul(Fp2Pack4x2x1& r, const Fp2Pack4x2x1& a, const Fp2Pack4x2x1& b);
        /// Fp2 square per pair: ((x0+x1)(x0-x1), 2*x0*x1), then redc.
        static void sqr(Fp2Pack4x2x1& r, const Fp2Pack4x2x1& a);
    };

    /**
     * Widened companion for 4x2x1 (identical lane map, UFp coefficients).
     *
     * Used when postponing Montgomery reduction across fused kernels.
     */
    struct UFp2Pack4x2x1 : public UFpPack8x1
    {
        /// Unreduced per-pair add.
        static void add(UFp2Pack4x2x1& r, const UFp2Pack4x2x1& a, const UFp2Pack4x2x1& b);
        /// Unreduced per-pair sub.
        static void sub(UFp2Pack4x2x1& r, const UFp2Pack4x2x1& a, const UFp2Pack4x2x1& b);
        /// Sparse pattern helper over [D,C,B,A] equivalent to ssaax2 at UFp layer.
        static void sax2(UFp2Pack4x2x1& r, const UFp2Pack4x2x1& a, const UFp2Pack4x2x1& b);
        /// Per packed Fp2, multiply by (u+1) in widened domain.
        static void mul_by_u_plus_1(UFp2Pack4x2x1& r, const UFp2Pack4x2x1& a);
        /// Schoolbook widened Fp2 multiply from reduced inputs.
        static void mul(UFp2Pack4x2x1& r, const Fp2Pack4x2x1& a, const Fp2Pack4x2x1& b);
        /// Specialized widened Fp2 square from reduced input.
        static void sqr(UFp2Pack4x2x1& r, const Fp2Pack4x2x1& a);
        /// Montgomery-reduce widened lanes back to reduced 4x2x1.
        static void redc(Fp2Pack4x2x1& r, const UFp2Pack4x2x1& a);
    };

    /**
     * Fp2 schedule 2x4x1 packed into one FpPack8x1.
     *
     * External wire format matches 4x2x1 lane order:
     *   [D1 | D0 | C1 | C0 | B1 | B0 | A1 | A0]
     *
     * Distinction from 4x2x1 is algorithmic scheduling:
     * - this layout is consumed by permutation patterns used in Fp4/Fp6 kernels.
     */
    struct Fp2Pack2x4x1 : public FpPack8x1
    {
        /// Schoolbook Fp2 multiply + reduction in 2x4x1 permutation schedule.
        /// Accepts any layout-compatible `FpPack8x1`-derived arguments.
        static void mul(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b);
    };

    /** Widened companion for 2x4x1 schedule (same lane order, UFp coefficients). */
    struct UFp2Pack2x4x1 : public UFpPack8x1
    {
        /// Schoolbook widened multiply in 2x4x1 permutation schedule.
        /// Accepts any layout-compatible `FpPack8x1`/`UFpPack8x1`-derived arguments.
        static void mul(UFpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b);
    };

    /**
     * Fp2 schedule 2x2x2 packed into one FpPack4x2.
     *
     * Let:
     *   P = P0 + u*P1
     *   Q = Q0 + u*Q1
     *
     * Per-word lane layout (high lane -> low lane):
     *   [Q1.hi | Q1.lo | Q0.hi | Q0.lo | P1.hi | P1.lo | P0.hi | P0.lo]
     * where `.lo/.hi` denote low/high limb groups in the 4x2 decomposition.
     *
     * Kernel model matches 4x2x1 formulas but operates on FpPack4x2 halves.
     */
    struct Fp2Pack2x2x2 : public FpPack4x2
    {
        /// Per packed Fp2 add.
        static void add(Fp2Pack2x2x2& r, const Fp2Pack2x2x2& a, const Fp2Pack2x2x2& b);
        /// Per packed Fp2 sub.
        static void sub(Fp2Pack2x2x2& r, const Fp2Pack2x2x2& a, const Fp2Pack2x2x2& b);
        /// Mixed pattern over [Q,P]: add on Q, sub on P.
        static void as(Fp2Pack2x2x2& r, const Fp2Pack2x2x2& a, const Fp2Pack2x2x2& b);
        /// Mixed pattern over [Q,P]: sub on Q, add on P.
        static void sa(Fp2Pack2x2x2& r, const Fp2Pack2x2x2& a, const Fp2Pack2x2x2& b);
        /// Per packed Fp2, multiply by (u+1): (x0,x1)->(x0-x1, x0+x1).
        static void mul_by_u_plus_1(Fp2Pack2x2x2& r, const Fp2Pack2x2x2& a);
        /// Schoolbook Fp2 mul per pair, then redc.
        static void mul(Fp2Pack2x2x2& r, const Fp2Pack2x2x2& a, const Fp2Pack2x2x2& b);
        /// Specialized Fp2 square per pair, then redc.
        static void sqr(Fp2Pack2x2x2& r, const Fp2Pack2x2x2& a);
    };

    /** Widened companion for 2x2x2 schedule (same lane map, UFp coefficients). */
    struct UFp2Pack2x2x2 : public UFpPack4x2
    {
        /// Montgomery-reduce widened 2x2x2 value into reduced form.
        static void redc(Fp2Pack2x2x2& r, const UFp2Pack2x2x2& a);
        /// Unreduced per-pair add.
        static void add(UFp2Pack2x2x2& r, const UFp2Pack2x2x2& a, const UFp2Pack2x2x2& b);
        /// Unreduced per-pair sub.
        static void sub(UFp2Pack2x2x2& r, const UFp2Pack2x2x2& a, const UFp2Pack2x2x2& b);
        /// Per packed Fp2, multiply by (u+1) in widened domain.
        static void mul_by_u_plus_1(UFp2Pack2x2x2& r, const UFp2Pack2x2x2& a);
        /// Schoolbook widened Fp2 multiply from reduced inputs.
        static void mul(UFp2Pack2x2x2& r, const Fp2Pack2x2x2& a, const Fp2Pack2x2x2& b);
        /// Specialized widened Fp2 square from reduced input.
        static void sqr(UFp2Pack2x2x2& r, const Fp2Pack2x2x2& a);
    };

    /**
     * Fp2 schedule 1x4x2 packed into one FpPack4x2.
     *
     * External wire format matches 2x2x2:
     *   [Q1.hi | Q1.lo | Q0.hi | Q0.lo | P1.hi | P1.lo | P0.hi | P0.lo]
     *
     * Distinction is algorithmic scheduling only: permutation plan optimized for
     * 1x4x2 consumers.
     */
    struct Fp2Pack1x4x2 : public FpPack4x2
    {
        /// Schoolbook Fp2 multiply + reduction in 1x4x2 permutation schedule.
        /// Accepts any layout-compatible `FpPack4x2`-derived arguments.
        static void mul(FpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b);
    };

    /** Widened companion for 1x4x2 schedule (same lane order, UFp coefficients). */
    struct UFp2Pack1x4x2 : public UFpPack4x2
    {
        /// Schoolbook widened multiply in 1x4x2 permutation schedule.
        /// Accepts any layout-compatible `FpPack4x2`/`UFpPack4x2`-derived arguments.
        static void mul(UFpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b);
    };








    BLST_FORCEINLINE void Fp2Pack8x1x1::add(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
    {
        FpPack8x1::add(r[0], a[0], b[0]);
        FpPack8x1::add(r[1], a[1], b[1]);
    }

    BLST_FORCEINLINE void Fp2Pack8x1x1::sub(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
    {
        FpPack8x1::sub(r[0], a[0], b[0]);
        FpPack8x1::sub(r[1], a[1], b[1]);
    }

    // r0 = a0*b0 - a1*b1
    // r1 = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1 = a0*b1 + a1*b0
    // Karatsuba
    inline void UFp2Pack8x1x1::mul(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
    {
        FpPack8x1 t0, t1;
        UFpPack8x1 tt0, tt1, tt2;

        FpPack8x1::mul(tt0, a[0], b[0]);       // a0*b0
        FpPack8x1::mul(tt1, a[1], b[1]);       // a1*b1
        FpPack8x1::add(t0, a[0], a[1]);          // a0+a1
        FpPack8x1::add(t1, b[0], b[1]);          // b0+b1
        FpPack8x1::mul(tt2, t0, t1);           // (a0+a1)*(b0+b1)
        UFpPack8x1::sub(r[0], tt0, tt1);        // a0*b0-a1*b1
        UFpPack8x1::sub(tt2, tt2, tt0);         // (a0+a1)*(b0+b1)-a0*b0
        UFpPack8x1::sub(r[1], tt2, tt1);        // (a0+a1)*(b0+b1)-a0*b0-a1*b1
    }



    static_assert(std::is_trivially_default_constructible_v<Fp2Pack8x1x1>);
    static_assert(std::is_trivially_copyable_v<Fp2Pack8x1x1>);
    static_assert(std::is_trivially_default_constructible_v<UFp2Pack8x1x1>);
    static_assert(std::is_trivially_copyable_v<UFp2Pack8x1x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp2Pack4x2x1>);
    static_assert(std::is_trivially_copyable_v<Fp2Pack4x2x1>);
    static_assert(std::is_trivially_default_constructible_v<UFp2Pack4x2x1>);
    static_assert(std::is_trivially_copyable_v<UFp2Pack4x2x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp2Pack2x4x1>);
    static_assert(std::is_trivially_copyable_v<Fp2Pack2x4x1>);
    static_assert(std::is_trivially_default_constructible_v<UFp2Pack2x4x1>);
    static_assert(std::is_trivially_copyable_v<UFp2Pack2x4x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp2Pack2x2x2>);
    static_assert(std::is_trivially_copyable_v<Fp2Pack2x2x2>);
    static_assert(std::is_trivially_default_constructible_v<UFp2Pack2x2x2>);
    static_assert(std::is_trivially_copyable_v<UFp2Pack2x2x2>);
    static_assert(std::is_trivially_default_constructible_v<Fp2Pack1x4x2>);
    static_assert(std::is_trivially_copyable_v<Fp2Pack1x4x2>);
    static_assert(std::is_trivially_default_constructible_v<UFp2Pack1x4x2>);
    static_assert(std::is_trivially_copyable_v<UFp2Pack1x4x2>);
}

#endif

#endif // BLST_FIELDS_FP2_AVX_H
