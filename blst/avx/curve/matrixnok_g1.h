// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "ec.h"
#include "scalar_mult.h"
#include "../rns/matrixnok_fp.h"

#include <type_traits>

namespace blst::avx::matrixnok
{
    // Isolated production-shaped G1 candidate. It deliberately has no public
    // dispatch entry: callers must opt into the fixed MatrixNoK ring while its
    // correctness, code generation, and performance gates are being closed.
    class G1CurveCandidate
    {
    public:
        using ProjPoint = ProjectivePoint<FpRingCandidate::StandardElement>;
        using AffPoint = AffinePoint<FpRingCandidate::StandardElement>;

        // Keep the first multiplication round in the small m2 representation.
        // Expanding only the normalized cross terms avoids redundant m1 work
        // while preserving the RNS multiple needed by the second round.
        static INLINE ProjPoint add_point(
            const ProjPoint& left, const ProjPoint& right,
            const FpRingCandidate& ring)
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
            return ProjPoint(x, y, z);
        }

        static INLINE ProjPoint add_mixed_point(
            const ProjPoint& left, const AffPoint& right,
            const FpRingCandidate& ring)
        {
            return add_point(left, ProjPoint(right.x, right.y, ring.one()), ring)
                .ctime_replace_if(PointIsZero(right), left);
        }

        static INLINE ProjPoint double_point(
            const ProjPoint& point, const FpRingCandidate& ring)
        {
            return PointDouble<FpRingCandidate>(point, ring);
        }

        static INLINE ProjPoint dadd(
            const ProjPoint& left, const ProjPoint& right,
            const FpRingCandidate& ring,
            const FpRingCandidate::StandardElement* a4)
        {
            ProjPoint result = a4 == nullptr
                ? add_point(left, right, ring)
                : PointDaddA<FpRingCandidate>(left, right, *a4, ring);
            result = result.ctime_replace_if(PointCoordIsZero(right.z), left);
            return result.ctime_replace_if(PointCoordIsZero(left.z), right);
        }

        static INLINE ProjPoint negate(
            const ProjPoint& point, const FpRingCandidate& ring)
        {
            return ProjPoint(point.x, ring.standard_negate(point.y), point.z);
        }

        static INLINE ProjPoint zero(const FpRingCandidate& ring)
        {
            return ProjPoint(ring.zero(), ring.one(), ring.zero());
        }
    };
    using G1ProjectiveCandidate = G1CurveCandidate::ProjPoint;
    using G1AffineCandidate = G1CurveCandidate::AffPoint;
    using G1GLVCandidate = GLVMult<FpRingCandidate, G1CurveCandidate>;

    static_assert(std::is_empty_v<FpRingCandidate>);
    static_assert(std::is_empty_v<G1CurveCandidate>);
} // namespace blst::avx::matrixnok
