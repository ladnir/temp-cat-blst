/*
 * Copyright 2026 Simon Langowski
 * Copyright 2026 Category Labs
 * SPDX-License-Identifier: Apache-2.0 AND MIT
 *
 * Adapted from the VROOM G1 construction at commit d6ecd534.
 */
#pragma once

#include "ec.h"
#include "matrixnok_curve.h"
#include "scalar_mult.h"
#include "../rns/matrixnok_fp.h"

#include <type_traits>

namespace blst::avx::matrixnok
{
    /**
     * Selected fixed G1 construction.
     *
     * Arithmetic accepts valid on-curve homogeneous points; subgroup
     * membership is unnecessary. No operation performs Release validation,
     * allocates, throws, or uses shared mutable state. Point-data control uses
     * fixed schedules and conditional selection. Return values permit exact
     * outer output/input aliasing.
     *
     * `add_mixed_point` recognizes affine infinity only as exact `(0,0)`.
     * `dadd` recognizes projective infinity through an exact zero RNS `z`.
     * A null `a4` selects the a=0 complete-add schedule. Otherwise, `*a4`
     * must be the coefficient of the curve containing both inputs.
     */
    class G1Curve
    {
    public:
        using ProjPoint = ProjectivePoint<FpRingCandidate::StandardElement>;
        using AffPoint = AffinePoint<FpRingCandidate::StandardElement>;

        static BLST_FORCEINLINE ProjPoint add_point(
            const ProjPoint& left, const ProjPoint& right,
            const FpRingCandidate& ring)
        {
            return detail::add_point(left, right, ring);
        }

        static BLST_FORCEINLINE ProjPoint add_mixed_point(
            const ProjPoint& left, const AffPoint& right,
            const FpRingCandidate& ring)
        {
            return add_point(left, ProjPoint(right.x, right.y, ring.one()), ring)
                .ctime_replace_if(AffinePointIsExactIdentity(right), left);
        }

        static BLST_FORCEINLINE ProjPoint double_point(
            const ProjPoint& point, const FpRingCandidate& ring)
        {
            return PointDouble<FpRingCandidate>(point, ring);
        }

        static BLST_FORCEINLINE ProjPoint dadd(
            const ProjPoint& left, const ProjPoint& right,
            const FpRingCandidate& ring,
            const FpRingCandidate::StandardElement* a4)
        {
            ProjPoint result = a4 == nullptr
                ? add_point(left, right, ring)
                : PointDaddA<FpRingCandidate>(left, right, *a4, ring);
            result = result.ctime_replace_if(PointCoordIsExactZero(right.z), left);
            return result.ctime_replace_if(PointCoordIsExactZero(left.z), right);
        }

        static BLST_FORCEINLINE ProjPoint negate(
            const ProjPoint& point, const FpRingCandidate& ring)
        {
            return ProjPoint(point.x, ring.standard_negate(point.y), point.z);
        }

        static BLST_FORCEINLINE ProjPoint zero(const FpRingCandidate& ring)
        {
            return ProjPoint(ring.zero(), ring.one(), ring.zero());
        }
    };
    using G1Projective = G1Curve::ProjPoint;
    using G1Affine = G1Curve::AffPoint;
    using G1GLV = GLVMult<FpRingCandidate, G1Curve>;

    static_assert(std::is_empty_v<FpRingCandidate>);
    static_assert(std::is_empty_v<G1Curve>);
} // namespace blst::avx::matrixnok
