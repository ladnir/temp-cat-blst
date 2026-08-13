/*
 * Copyright 2026 Simon Langowski
 * Copyright 2026 Category Labs
 * SPDX-License-Identifier: Apache-2.0 AND MIT
 *
 * The complete-add schedule is adapted from VROOM commit d6ecd534. Category
 * Labs integrated it with CatBLST's bounded MatrixNoK representation and
 * proof correspondence.
 */
#pragma once

#include "ec.h"

namespace blst::avx::matrixnok::detail
{
    /**
     * Complete homogeneous addition with signed small-basis preservation.
     *
     * MatrixNoK keeps the first multiplication round in m2. Signed cross
     * terms must be lifted before expansion so reconstructing m1 preserves the
     * represented integer rather than only its residue. The selected Fp and
     * Fp2 policies use this exact schedule.
     *
     * Both inputs must be valid homogeneous points on the same selected a=0
     * curve. Subgroup membership is unnecessary. The fixed schedule performs
     * no validation, allocation, or throwing operation and has no shared
     * mutable state. Return-by-value construction permits exact outer
     * output/input aliasing.
     */
    template<class Ring>
    static BLST_FORCEINLINE ProjectivePoint<typename Ring::StandardElement> add_point(
        const ProjectivePoint<typename Ring::StandardElement>& left,
        const ProjectivePoint<typename Ring::StandardElement>& right,
        const Ring& ring)
    {
        const auto left_xy = left.x + left.y;
        const auto left_xz = left.x + left.z;
        const auto left_yz = left.y + left.z;
        const auto right_xy = right.x + right.y;
        const auto right_yz = right.y + right.z;
        const auto right_xz = right.x + right.z;

        const auto xx_wide = ring.prep_left(left.x) * right.x;
        const auto yy_wide = ring.prep_left(left.y) * right.y;
        const auto zz_wide = ring.prep_left(left.z) * right.z;
        const auto xy_wide =
            ring.prep_left(left_xy) * ring.prep(right_xy);
        const auto yz_wide =
            ring.prep_left(left_yz) * ring.prep(right_yz);
        const auto xz_wide =
            ring.prep_left(left_xz) * ring.prep(right_xz);

        const auto [xx, yy, zz, xy, yz, xz] = ring.batch_reduce(
            xx_wide, yy_wide, zz_wide, xy_wide, yz_wide, xz_wide);
        const auto xx_plus_yy = xx + yy;
        const auto yy_plus_zz = yy + zz;
        const auto xx_plus_zz = xx + zz;
        const auto xz_cross = xz - xx_plus_zz;
        const auto two_xx = xx + xx;
        const auto b3_zz = ring.mul_3b(zz);

        const auto [yz_cross, xy_cross, b3_xz_cross,
            yy_minus_b3zz, three_xx, yy_plus_b3zz] = ring.batch_expand(
                ring.prep_expand(yz - yy_plus_zz),
                ring.prep_expand(xy - xx_plus_yy),
                ring.prep_expand(ring.mul_3b(xz_cross)),
                ring.prep_expand(yy - b3_zz),
                ring.prep_expand(xx + two_xx),
                ring.prep_expand(yy + b3_zz));

        const auto x_negative = ring.prep_left(yz_cross) *
            ring.negate(b3_xz_cross);
        const auto x_positive = ring.prep_left(yy_minus_b3zz) * xy_cross;
        const auto y_left = ring.prep_left(three_xx) * b3_xz_cross;
        const auto y_right = ring.prep_left(yy_minus_b3zz) * yy_plus_b3zz;
        const auto z_left = ring.prep_left(yz_cross) * yy_plus_b3zz;
        const auto z_right = ring.prep_left(three_xx) * xy_cross;

        const auto [x, y, z] = ring.batch_reduce_expand(
            x_positive + x_negative,
            y_right + y_left,
            z_left + z_right);
        return {x, y, z};
    }
} // namespace blst::avx::matrixnok::detail
