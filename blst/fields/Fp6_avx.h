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

    /// Three Fp2 components (`c0,c1,c2`) packed in the `8x1x1w` schedule.
    using Fp6Pack4x2x1x1 = AvxPackTuple<Fp2Pack8x1x1, 3>;
    using fp6_4x2x1x1w = Fp6Pack4x2x1x1;
    /// Three Fp2 components (`c0,c1,c2`) packed in the `4x2x1w` schedule.
    using Fp6Pack2x2x2x1 = AvxPackTuple<Fp2Pack4x2x1, 3>;
    using fp6_2x2x2x1w = Fp6Pack2x2x2x1;
    /// Three Fp2 components (`c0,c1,c2`) packed in the `2x2x2w` schedule.
    using Fp6Pack1x2x2x2 = AvxPackTuple<Fp2Pack2x2x2, 3>;
    using fp6_1x2x2x2w = Fp6Pack1x2x2x2;

    static_assert(std::is_trivially_default_constructible_v<Fp6Pack4x2x1x1>);
    static_assert(std::is_trivially_copyable_v<Fp6Pack4x2x1x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp6Pack2x2x2x1>);
    static_assert(std::is_trivially_copyable_v<Fp6Pack2x2x2x1>);
    static_assert(std::is_trivially_default_constructible_v<Fp6Pack1x2x2x2>);
    static_assert(std::is_trivially_copyable_v<Fp6Pack1x2x2x2>);

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
     * Output layout:
     * - `r01` stores widened `r1 | r0`
     * - `r2` stores widened `... | r2` (high half is scratch schedule data).
     *
     * @pre Inputs are reduced Montgomery Fp2 values in avxbls `8x1x1w` packing.
     * @post Outputs are widened (unreduced) Montgomery intermediates.
     */
    void mul_fp6x2_4x2x1x1w(
        UFp2Pack8x1x1& r01,
        UFp2Pack8x1x1& r2,
        const Fp2Pack8x1x1& ab0,
        const Fp2Pack8x1x1& ab1,
        const Fp2Pack8x1x1& ab2);

    /**
     * Karatsuba Fp6 multiply (2-way schedule).
     *
     * Input/output meanings match `mul_fp6x2_4x2x1x1w`, but in `4x2x1w` packing.
     *
     * @pre Inputs are reduced Montgomery Fp2 values in avxbls `4x2x1w` packing.
     * @post Outputs are widened (unreduced) Montgomery intermediates.
     */
    void mul_fp6x2_2x2x2x1w(
        UFp2Pack4x2x1& r01,
        UFp2Pack4x2x1& r2,
        const Fp2Pack4x2x1& ab0,
        const Fp2Pack4x2x1& ab1,
        const Fp2Pack4x2x1& ab2);

    /**
     * Karatsuba Fp6 multiply (1-way schedule).
     *
     * Input/output meanings match `mul_fp6x2_4x2x1x1w`, but in `2x2x2w` packing.
     *
     * @pre Inputs are reduced Montgomery Fp2 values in avxbls `2x2x2w` packing.
     * @post Outputs are widened (unreduced) Montgomery intermediates.
     */
    void mul_fp6x2_1x2x2x2w(
        UFp2Pack2x2x2& r01,
        UFp2Pack2x2x2& r2,
        const Fp2Pack2x2x2& ab0,
        const Fp2Pack2x2x2& ab1,
        const Fp2Pack2x2x2& ab2);

    /**
     * Sparse Fp6 multiply for `(x,y,0,0,z,0)` pattern (2-way schedule).
     *
     * Input contract:
     * - `a01 = a1 | a0`
     * - `a2  = a2 | a2`
     * - `b01 = b1 | b0`
     * - `b4  = b4 | b4`
     *
     * Mathematical outputs:
     * - `r0 = a2*b1*(u+1) + a0*b0`
     * - `r1 = a0*b1 + a1*b0`
     * - `r2 = a2*b0 + a1*b1`
     * - `r3 = a2*b4*(u+1)`
     * - `r4 = a0*b4`
     * - `r5 = a1*b4`
     *
     * Output layout:
     * - `r01` stores widened `r1 | r0`
     * - `r2` stores widened `... | r2`
     * - `r3` stores widened `r3 | a2*b1*(u+1)` schedule pair
     * - `r45` stores widened `r5 | r4`
     *
     * @pre Inputs are reduced Montgomery Fp2 values in avxbls `4x2x1w` packing.
     * @post Outputs are widened (unreduced) Montgomery intermediates.
     */
    void mul_by_xy00z0_fp6x2_2x2x2x1w(
        UFp2Pack4x2x1& r01,
        UFp2Pack4x2x1& r2,
        UFp2Pack4x2x1& r3,
        UFp2Pack4x2x1& r45,
        const Fp2Pack4x2x1& a01,
        const Fp2Pack4x2x1& a2,
        const Fp2Pack4x2x1& b01,
        const Fp2Pack4x2x1& b4);

    /**
     * Sparse Fp6 multiply for `(x,y,0,0,z,0)` pattern (4-way schedule).
     *
     * Input contract:
     * - `a = a2 | a2 | a1 | a0` at Fp2 layer
     * - `b = b0 | b1 | b0 | b4`
     *
     * This routine emits widened schedule blocks used by higher-level Fp12 kernels:
     * - `r01` carries widened schedule for `r0`/`r1`
     * - `r23` carries widened schedule for `r2`/`r3`
     * - `r45` carries widened schedule for `r4`/`r5`
     *
     * @pre Inputs are reduced Montgomery Fp2 values in avxbls `8x1x1w` packing.
     * @post Outputs are widened (unreduced) Montgomery intermediates.
     */
    void mul_by_xy00z0_fp6x2_2x4x1x1w(
        UFp2Pack8x1x1& r01,
        UFp2Pack8x1x1& r23,
        UFp2Pack8x1x1& r45,
        const Fp2Pack8x1x1& a,
        const Fp2Pack8x1x1& b);
}

#endif

#endif // BLST_FIELDS_FP6_AVX_H
