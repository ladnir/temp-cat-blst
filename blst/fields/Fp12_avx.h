#ifndef BLST_FIELDS_FP12_AVX_H
#define BLST_FIELDS_FP12_AVX_H

#include "blst/fields/Fp6_avx.h"
#include "blst/fields/Fp12.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    /**
     * AVX512 Fp12 helpers ported from avxbls.
     *
     * Fp12 model:
     * - element = a0 + a1 * w, where a0,a1 are Fp6
     * - each Fp6 = c0 + c1 * v + c2 * v^2 over Fp2 with v^3=(u+1)
     * - all component values are Montgomery-form Fp2 values.
     *
     * These kernels operate on AVXBLS schedules. Some outputs are schedule blocks
     * and not direct lane-by-lane scalar tuples unless explicitly stated.
     */

    /**
     * Cyclotomic square (double-length Fp4 square path).
     *
     * Computes the Granger-Scott cyclotomic square update:
     * - `ra` stores the updated A block (`z0,z1` side)
     * - `rbc` stores the updated B/C block (`z2..z5` side)
     *
     * @pre `a` and `bc` are valid packed cyclotomic schedule inputs in Montgomery form.
     * @post `ra` and `rbc` contain reduced Montgomery schedule outputs.
     */
    void cyclotomic_sqr_fp12_vec_v1(
        fp4_1x2x2x2w& ra,
        fp4_2x2x2x1w& rbc,
        const fp4_1x2x2x2w& a,
        const fp4_2x2x2x1w& bc);

    /**
     * Cyclotomic square (single-length Fp4 square path).
     *
     * Same math as `cyclotomic_sqr_fp12_vec_v1`, using the v2 Fp4 square kernels.
     *
     * @pre `a` and `bc` are valid packed cyclotomic schedule inputs in Montgomery form.
     * @post `ra` and `rbc` contain reduced Montgomery schedule outputs.
     */
    void cyclotomic_sqr_fp12_vec_v2(
        fp4_1x2x2x2w& ra,
        fp4_2x2x2x1w& rbc,
        const fp4_1x2x2x2w& a,
        const fp4_2x2x2x1w& bc);

    /**
     * Compressed cyclotomic square update of only the B/C block.
     *
     * This is the B/C half of the cyclotomic square formula used by compressed
     * exponentiation schedules.
     *
     * @pre `bc` is a valid packed B/C block in Montgomery form.
     * @post `rbc` contains reduced Montgomery B/C schedule output.
     */
    void compressed_cyclotomic_sqr_fp12_vec_v1(
        fp4_2x2x2x1w& rbc,
        const fp4_2x2x2x1w& bc);

    /**
     * Schoolbook Fp12 multiply, variant 1.
     *
     * Produces reduced schedule blocks:
     * - `r01` carries `(r1[1], r1[0], r0[2], r0[1])`
     * - `r2` carries `(r1[2], r0[0])` in the low schedule region.
     *
     * @pre `ab0/ab1/ab2` are valid `8x1x1w` AVXBLS packed products schedule inputs.
     * @post `r01` and `r2` contain reduced Montgomery schedule outputs.
     */
    void mul_fp12_vec_v1(
        fp2_4x2x1w& r01,
        fp2_4x2x1w& r2,
        const fp2_8x1x1w& ab0,
        const fp2_8x1x1w& ab1,
        const fp2_8x1x1w& ab2);

    /**
     * Schoolbook Fp12 multiply, variant 2 (double-length sub-routine path).
     *
     * Output meaning matches variant 1, with `r2` emitted in `2x2x2w`.
     *
     * @pre `ab0/ab1/ab2` are valid `8x1x1w` AVXBLS packed products schedule inputs.
     * @post `r01` and `r2` contain reduced Montgomery schedule outputs.
     */
    void mul_fp12_vec_v2(
        fp2_4x2x1w& r01,
        fp2_2x2x2w& r2,
        const fp2_8x1x1w& ab0,
        const fp2_8x1x1w& ab1,
        const fp2_8x1x1w& ab2);

    /**
     * Karatsuba Fp12 multiply, variant 3 (`2x2x2w` sub-routines).
     *
     * Outputs:
     * - `r001`: `(r0[1], r0[0])`
     * - `r101`: `(r1[1], r1[0])`
     * - `r2`: `(r1[2], r0[2])`
     *
     * @pre `ab0/ab1/ab2` are valid `4x2x1w` AVXBLS packed products schedule inputs.
     * @post outputs contain reduced Montgomery schedule blocks.
     */
    void mul_fp12_vec_v3(
        fp2_2x2x2w& r001,
        fp2_2x2x2w& r101,
        fp2_2x2x2w& r2,
        const fp2_4x2x1w& ab0,
        const fp2_4x2x1w& ab1,
        const fp2_4x2x1w& ab2);

    /**
     * Karatsuba Fp12 multiply, variant 4 (`4x2x1w` + `2x2x2w` split).
     *
     * Outputs:
     * - `r0`: low Fp6 block (`r0[0..2]`) schedule block
     * - `r101`: `(r1[1], r1[0])`
     * - `r12`: contains `r1[2]` in low schedule region
     *
     * @pre `ab0/ab1/ab2` are valid `4x2x1w` AVXBLS packed products schedule inputs.
     * @post outputs contain reduced Montgomery schedule blocks.
     */
    void mul_fp12_vec_v4(
        fp2_4x2x1w& r0,
        fp2_2x2x2w& r101,
        fp2_2x2x2w& r12,
        const fp2_4x2x1w& ab0,
        const fp2_4x2x1w& ab1,
        const fp2_4x2x1w& ab2);

    /**
     * Sparse Fp12 multiply by `(x,y,0,0,z,0)`, variant 1 (2-way sparse Fp6 path).
     *
     * @pre `a01/a2/b01/b4` follow AVXBLS sparse schedule packing.
     * @post `r0` and `r1` contain reduced Montgomery schedule blocks.
     */
    void mul_by_xy00z0_fp12_vec_v1(
        fp2_4x2x1w& r0,
        fp2_2x2x2w& r1,
        const fp2_4x2x1w& a01,
        const fp2_4x2x1w& a2,
        const fp2_4x2x1w& b01,
        const fp2_4x2x1w& b4);

    /**
     * Sparse Fp12 multiply by `(x,y,0,0,z,0)`, variant 2 (4-way sparse Fp6 path).
     *
     * @pre `a` and `b` follow AVXBLS sparse `8x1x1w` schedule packing.
     * @post `r0` and `r1` contain reduced Montgomery schedule blocks.
     */
    void mul_by_xy00z0_fp12_vec_v2(
        fp2_4x2x1w& r0,
        fp2_2x2x2w& r1,
        const fp2_8x1x1w& a,
        const fp2_8x1x1w& b);

    /**
     * Karatsuba-style Fp12 square, variant 1.
     *
     * Computes `r = a^2` in packed schedule form:
     * - `r0`: low Fp6 block schedule
     * - `r1`: high Fp6 block schedule
     *
     * @pre `a0/a1` are valid AVXBLS packed square-schedule inputs.
     * @post `r0` and `r1` contain reduced Montgomery schedule outputs.
     */
    void sqr_fp12_vec_v1(
        fp2_4x2x1w& r0,
        fp2_4x2x1w& r1,
        const fp2_4x2x1w& a0,
        const fp2_4x2x1w& a1);

    /**
     * AVXBLS wrapper: full Fp12 multiply with scalar-layout inputs/outputs.
     *
     * This function packs scalar `Fp12` operands into the AVXBLS `4x2x1w`
     * multiply schedule, runs `mul_fp12_vec_v3`, and unpacks back into scalar
     * Montgomery `Fp12` limbs.
     *
     * @pre `a` and `b` are Montgomery-form `Fp12` values.
     * @post `ret` equals `a * b` in Montgomery form.
     */
    void mul_fp12_vector(Fp12& ret, const Fp12& a, const Fp12& b);

    /**
     * AVXBLS wrapper: sparse Fp12 multiply by `(x,y,0,0,z,0)`.
     *
     * Packs scalar `Fp12/Fp6` operands into AVXBLS sparse schedules, runs
     * `mul_by_xy00z0_fp12_vec_v1`, then unpacks to scalar Montgomery layout.
     *
     * @pre `a` and `xy00z0` are Montgomery-form values.
     * @post `ret` equals `a * xy00z0` in Montgomery form.
     */
    void mul_by_xy00z0_fp12_vector(Fp12& ret, const Fp12& a, const Fp6& xy00z0);

    /**
     * AVXBLS wrapper: full Fp12 squaring.
     *
     * Packs scalar `a` into the AVXBLS square schedule, runs `sqr_fp12_vec_v1`,
     * and unpacks to scalar Montgomery layout.
     *
     * @pre `a` is a Montgomery-form `Fp12` value.
     * @post `ret` equals `a^2` in Montgomery form.
     */
    void sqr_fp12_vector(Fp12& ret, const Fp12& a);

    /**
     * AVXBLS wrapper: cyclotomic square with scalar-layout inputs/outputs.
     *
     * Packs scalar `a` into AVXBLS cyclotomic schedules, runs
     * `cyclotomic_sqr_fp12_vec_v1`, and unpacks to scalar Montgomery layout.
     *
     * @pre `a` is in the cyclotomic subgroup and in Montgomery form.
     * @post `ret` equals the cyclotomic square of `a`.
     */
    void cyclotomic_sqr_fp12_vector(Fp12& ret, const Fp12& a);
}

#endif

#endif // BLST_FIELDS_FP12_AVX_H
