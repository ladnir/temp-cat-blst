#ifndef BLST_FIELDS_FP6_AVX_H
#define BLST_FIELDS_FP6_AVX_H

#include "blst/fields/Fp4_avx.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    /**
     * AVX512 Fp6 helpers ported from avxbls.
     *
     * Field model:
     * - Fp6 element is `a0 + a1*v + a2*v^2` over Fp2.
     * - `v^3 = (u + 1)`.
     * - all components are Montgomery-form Fp2 values.
     *
     * Naming follows avxbls schedules so kernels can be compared line-for-line.
     */

    /**
     * Fp6 schedule backed by `Fp2Pack8x1x1` components.
     *
     * Element model:
     *   X = X.c0 + v*X.c1 + v^2*X.c2
     *   Y = Y.c0 + v*Y.c1 + v^2*Y.c2
     *
     * Component topology:
     *   packs[0] = c0 schedule block
     *   packs[1] = c1 schedule block
     *   packs[2] = c2 schedule block
     *
     * Visual component matrix (rows = Fp6 components, columns = Fp2 components):
     *   packs[0] : [re(c0) | im(c0)]
     *   packs[1] : [re(c1) | im(c1)]
     *   packs[2] : [re(c2) | im(c2)]
     *
     * Each `re/im` column uses inherited `Fp2Pack8x1x1` lane order.
     */
    struct Fp6Pack4x2x1x1 : public AvxPackTuple<Fp2Pack8x1x1, 3>
    {
        /**
         * Karatsuba Fp6 multiply (4-way schedule).
         *
         * Input contract:
         * - `ab0 = b0 | a0` at Fp2 layer
         * - `ab1 = b1 | a1`
         * - `ab2 = b2 | a2`
         *
         * Mathematical outputs:
         * - `r0 = (a1*b2 + a2*b1)*(u+1) + a0*b0`
         * - `r1 = a0*b1 + a1*b0 + a2*b2*(u+1)`
         * - `r2 = a0*b2 + a2*b0 + a1*b1`
         *
         * Output schedule:
         * - `r01` stores widened `r1 | r0`
         * - `r2` stores widened `... | r2` (high half is schedule scratch)
         */
        static void mul(
            UFp2Pack8x1x1& r01,
            UFp2Pack8x1x1& r2,
            const Fp2Pack8x1x1& ab0,
            const Fp2Pack8x1x1& ab1,
            const Fp2Pack8x1x1& ab2);

        /**
         * Sparse Fp6 multiply for `(x,y,0,0,z,0)` in 4-way schedule.
         *
         * Input contract:
         * - `a = a2 | a2 | a1 | a0` at Fp2 layer
         * - `b = b0 | b1 | b0 | b4`
         *
         * Output schedule blocks:
         * - `r01` carries widened schedule for `r0`/`r1`
         * - `r23` carries widened schedule for `r2`/`r3`
         * - `r45` carries widened schedule for `r4`/`r5`
         */
        static void mul_by_xy00z0(
            UFp2Pack8x1x1& r01,
            UFp2Pack8x1x1& r23,
            UFp2Pack8x1x1& r45,
            const Fp2Pack8x1x1& a,
            const Fp2Pack8x1x1& b);
    };

    /** Widened companion storage for `Fp6Pack4x2x1x1`. */
    struct UFp6Pack4x2x1x1 : public AvxPackTuple<UFp2Pack8x1x1, 3> {};

    /**
     * Fp6 schedule backed by `Fp2Pack4x2x1` components.
     *
     * Element model:
     *   X = X.c0 + v*X.c1 + v^2*X.c2
     *   Y = Y.c0 + v*Y.c1 + v^2*Y.c2
     *
     * Component topology:
     *   packs[0] = c0 schedule block
     *   packs[1] = c1 schedule block
     *   packs[2] = c2 schedule block
     *
     * Visual component matrix (rows = Fp6 components, columns = Fp2 components):
     *   packs[0] : [re(c0) | im(c0)]
     *   packs[1] : [re(c1) | im(c1)]
     *   packs[2] : [re(c2) | im(c2)]
     *
     * Each `re/im` column uses inherited `Fp2Pack4x2x1` lane order.
     */
    struct Fp6Pack2x2x2x1 : public AvxPackTuple<Fp2Pack4x2x1, 3>
    {
        /// Karatsuba Fp6 multiply in `4x2x1w` packing; output is widened `r01|r2` blocks.
        static void mul(
            UFp2Pack4x2x1& r01,
            UFp2Pack4x2x1& r2,
            const Fp2Pack4x2x1& ab0,
            const Fp2Pack4x2x1& ab1,
            const Fp2Pack4x2x1& ab2);

        /**
         * Sparse Fp6 multiply for `(x,y,0,0,z,0)` in `4x2x1w` packing.
         *
         * Input contract:
         * - `a01 = a1 | a0`
         * - `a2  = a2 | a2`
         * - `b01 = b1 | b0`
         * - `b4  = b4 | b4`
         *
         * Output schedule:
         * - `r01` stores widened `r1 | r0`
         * - `r2` stores widened `... | r2`
         * - `r3` stores widened `r3 | a2*b1*(u+1)` schedule pair
         * - `r45` stores widened `r5 | r4`
         */
        static void mul_by_xy00z0(
            UFp2Pack4x2x1& r01,
            UFp2Pack4x2x1& r2,
            UFp2Pack4x2x1& r3,
            UFp2Pack4x2x1& r45,
            const Fp2Pack4x2x1& a01,
            const Fp2Pack4x2x1& a2,
            const Fp2Pack4x2x1& b01,
            const Fp2Pack4x2x1& b4);
    };

    /** Widened companion storage for `Fp6Pack2x2x2x1`. */
    struct UFp6Pack2x2x2x1 : public AvxPackTuple<UFp2Pack4x2x1, 3> {};

    /**
     * Fp6 schedule backed by `Fp2Pack2x2x2` components.
     *
     * Element model:
     *   X = X.c0 + v*X.c1 + v^2*X.c2
     *
     * Component topology:
     *   packs[0] = c0 schedule block
     *   packs[1] = c1 schedule block
     *   packs[2] = c2 schedule block
     *
     * Visual component matrix (rows = Fp6 components, columns = Fp2 components):
     *   packs[0] : [re(c0) | im(c0)]
     *   packs[1] : [re(c1) | im(c1)]
     *   packs[2] : [re(c2) | im(c2)]
     *
     * Each `re/im` column uses inherited `Fp2Pack2x2x2` lane order.
     */
    struct Fp6Pack1x2x2x2 : public AvxPackTuple<Fp2Pack2x2x2, 3>
    {
        /// Karatsuba Fp6 multiply in `2x2x2w` packing; output is widened `r01|r2` blocks.
        static void mul(
            UFp2Pack2x2x2& r01,
            UFp2Pack2x2x2& r2,
            const Fp2Pack2x2x2& ab0,
            const Fp2Pack2x2x2& ab1,
            const Fp2Pack2x2x2& ab2);
    };

    /** Widened companion storage for `Fp6Pack1x2x2x2`. */
    struct UFp6Pack1x2x2x2 : public AvxPackTuple<UFp2Pack2x2x2, 3> {};

    static_assert(std::is_trivially_default_constructible_v<Fp6Pack4x2x1x1>);
    static_assert(std::is_trivially_copyable_v<Fp6Pack4x2x1x1>);
    static_assert(std::is_trivially_default_constructible_v<UFp6Pack4x2x1x1>);
    static_assert(std::is_trivially_copyable_v<UFp6Pack4x2x1x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp6Pack2x2x2x1>);
    static_assert(std::is_trivially_copyable_v<Fp6Pack2x2x2x1>);
    static_assert(std::is_trivially_default_constructible_v<UFp6Pack2x2x2x1>);
    static_assert(std::is_trivially_copyable_v<UFp6Pack2x2x2x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp6Pack1x2x2x2>);
    static_assert(std::is_trivially_copyable_v<Fp6Pack1x2x2x2>);
    static_assert(std::is_trivially_default_constructible_v<UFp6Pack1x2x2x2>);
    static_assert(std::is_trivially_copyable_v<UFp6Pack1x2x2x2>);
}

#endif

#endif // BLST_FIELDS_FP6_AVX_H
