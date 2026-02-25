#ifndef BLST_CURVE_G2_AVX_H
#define BLST_CURVE_G2_AVX_H

#include "blst/fields/Fp2_avx.h"
#include "blst/fields/Fp12.h"
#include "blst/curve/point.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    /**
     * AVX512 G2 line-kernel helpers ported from avxbls.
     *
     * These functions implement packed Miller-loop line step kernels. Inputs and outputs
     * are schedule-packed Fp2 values in Montgomery form.
     */

    /**
     * Doubling line kernel (`dbl-2009-alnr`) in `2x2x2w` schedule.
     *
     * Input schedules:
     * - `X1Y1 = X1 | Y1` (Fp2 layer)
     * - `Z1   = Z1 | Z1`
     *
     * Output schedules:
     * - `l0Y3 = Y3 | l0`
     * - `l1` holds `l1` in channels 0/2
     * - `l2` holds `l2` in channels 4/6
     * - `X3` holds `X3` in channels 4/6
     * - `Z3` holds `Z3` in channels 4/6
     *
     * Scalar extraction used by avxbls:
     * - `line0 = (l0Y3 ch0, l0Y3 ch2)`
     * - `line1 = (l1 ch0, l1 ch2)`
     * - `line2 = (l2 ch4, l2 ch6)`
     * - `Y3    = (l0Y3 ch4, l0Y3 ch6)`
     * - `X3    = (X3 ch4, X3 ch6)`
     * - `Z3    = (Z3 ch4, Z3 ch6)`
     */
    void line_dbl_vec_v1(
        fp2_2x2x2w& l0Y3,
        fp2_2x2x2w& l1,
        fp2_2x2x2w& l2,
        fp2_2x2x2w& X3,
        fp2_2x2x2w& Z3,
        const fp2_2x2x2w& X1Y1,
        const fp2_2x2x2w& Z1);

    /**
     * Doubling line kernel (`dbl-2009-alnr`) in `4x2x1w` schedule.
     *
     * Input schedule:
     * - `X1Y1Z1 = Y1 | X1 | Z1 | Y1` (Fp2 layer)
     *
     * Output schedule:
     * - `l0` stores `l0` in channels 2/3
     * - `l12` stores `l1` in channels 4/5 and `l2` in channels 2/3
     * - `X3` stores `X3` in channels 6/7
     * - `Y3` stores `Y3` in channels 6/7
     * - `Z3` stores `Z3` in channels 0/1
     */
    void line_dbl_vec_v2(
        fp2_4x2x1w& l0,
        fp2_4x2x1w& l12,
        fp2_4x2x1w& X3,
        fp2_4x2x1w& Y3,
        fp2_4x2x1w& Z3,
        const fp2_4x2x1w& X1Y1Z1);

    /**
     * Addition line kernel (`madd-2007-bl`) in `2x2x2w` schedule.
     *
     * Function behavior:
     * - evaluates one mixed-addition Miller step for `R=(X1,Y1,Z1)` and affine `Q=(X2,Y2)`,
     * - computes updated Jacobian state `(X3,Y3,Z3)`,
     * - emits line coefficients `(l0,l1,l2)` where `l2` aliases the packed `Z3` output lanes.
     *
     * Input schedules:
     * - `X1Y1 = X1 | Y1`
     * - `Z1Y2 = Z1 | Y2`
     * - `X2   = X2 | X2`
     *
     * Preconditions:
     * - `R` and `Q` are valid on-curve Montgomery points encoded in the schedules above.
     * - This is the non-exceptional `madd-2007-bl` path; require
     *   `H = X2*Z1^2 - X1 != 0` in Fp2.
     * - Inputs that violate the mixed-add preconditions (same point, opposite point, or
     *   malformed non-curve tuples) are outside this kernel contract.
     *
     * Output schedules:
     * - `l0Y3 = Y3 | l0`
     * - `l1` holds `l1` in channels 0/2
     * - `X3` holds `X3` in channels 4/6
     * - `Z3` holds `Z3` in channels 4/6 and also line `l2`
     *
     * Scalar extraction used by avxbls:
     * - `line0 = (l0Y3 ch0, l0Y3 ch2)`
     * - `line1 = (l1 ch0, l1 ch2)`
     * - `line2 = (Z3 ch4, Z3 ch6)`
     * - `Y3    = (l0Y3 ch4, l0Y3 ch6)`
     * - `X3    = (X3 ch4, X3 ch6)`
     * - `Z3    = (Z3 ch4, Z3 ch6)`
     */
    void line_add_vec_v1(
        fp2_2x2x2w& l0Y3,
        fp2_2x2x2w& l1,
        fp2_2x2x2w& X3,
        fp2_2x2x2w& Z3,
        const fp2_2x2x2w& X1Y1,
        const fp2_2x2x2w& Z1Y2,
        const fp2_2x2x2w& X2);

    /**
     * Scale packed line coefficients by packed `Px2` constants in `4x2x1w`.
     *
     * Input contract:
     * - `a` carries `... | l1 | l2 | ...` at Fp2 layer.
     * - `Px2` carries `... | X | Y | ...` at Fp2 layer.
     *
     * Postcondition:
     * - channels containing `l1` are multiplied by `X`
     * - channels containing `l2` are multiplied by `Y`
     */
    void line_by_Px2_4x2x1w(
        fp2_4x2x1w& r,
        const fp2_4x2x1w& a,
        const fp2_4x2x1w& Px2);

    /**
     * Scale packed line coefficients by packed `Px2` constants in `2x2x2w`.
     *
     * Input contract:
     * - `a = l2 | l1`
     * - `Px2 = Y | X`
     *
     * Postcondition:
     * - output keeps `l2*Y | l1*X` in the same schedule.
     */
    void line_by_Px2_2x2x2w(
        fp2_2x2x2w& r,
        const fp2_2x2x2w& a,
        const fp2_2x2x2w& Px2);

    /**
     * AVXBLS wrapper: mixed-addition line step in scalar point layout.
     *
     * Packs scalar Jacobian `R` and affine `Q` into AVXBLS `2x2x2w` schedules,
     * runs `line_add_vec_v1`, and unpacks line coefficients plus updated Jacobian
     * state to scalar Montgomery `Fp2` coordinates.
     *
     * @pre `R` and `Q` satisfy the non-exceptional `madd-2007-bl` preconditions.
     * @post `line` holds `(l0,l1,l2)` and `T` holds the updated Jacobian point.
     */
    void line_add_vector(Fp6& line, G2Point& T, const G2Point& R, const G2AffinePoint& Q);

    /**
     * AVXBLS wrapper: doubling line step in scalar point layout.
     *
     * Packs scalar Jacobian `Q` into AVXBLS `4x2x1w` schedules, runs
     * `line_dbl_vec_v2`, and unpacks line coefficients plus updated Jacobian
     * state to scalar Montgomery `Fp2` coordinates.
     *
     * @pre `Q` is a Jacobian point in Montgomery coordinates.
     * @post `line` holds `(l0,l1,l2)` and `T` holds `2*Q`.
     */
    void line_dbl_vector(Fp6& line, G2Point& T, const G2Point& Q);
}

#endif

#endif // BLST_CURVE_G2_AVX_H
