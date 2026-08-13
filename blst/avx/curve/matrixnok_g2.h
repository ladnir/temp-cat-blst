/*
 * Copyright 2026 Simon Langowski
 * Copyright 2026 Category Labs
 * SPDX-License-Identifier: Apache-2.0 AND MIT
 *
 * Adapted from the VROOM G2 construction at commit d6ecd534.
 */
#pragma once

#include "ec.h"
#include "matrixnok_curve.h"
#include "scalar_mult.h"
#include "../fields/fp2.h"
#include "../rns/matrixnok_fp.h"

#include <type_traits>

namespace blst::avx::matrixnok
{
    /**
     * Fixed G2 construction selected by production dispatch.
     *
     * The selected construction reuses the shared complete homogeneous
     * schedule and generic `FP2Ring`; this facade fixes their parameter set.
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
    using G2Fp2Ring = FP2Ring<FpRingCandidate>;

    class G2Curve
    {
    public:
        using ProjPoint = ProjectivePoint<
            G2Fp2Ring::StandardElement>;
        using AffPoint = AffinePoint<
            G2Fp2Ring::StandardElement>;

        static BLST_FORCEINLINE ProjPoint add_point(
            const ProjPoint& left, const ProjPoint& right,
            const G2Fp2Ring& ring)
        {
            return detail::add_point(left, right, ring);
        }

        static BLST_FORCEINLINE ProjPoint add_mixed_point(
            const ProjPoint& left, const AffPoint& right,
            const G2Fp2Ring& ring)
        {
            return add_point(left, ProjPoint(right.x, right.y, ring.one()), ring)
                .ctime_replace_if(AffinePointIsExactIdentity(right), left);
        }

        static BLST_FORCEINLINE ProjPoint double_point(
            const ProjPoint& point, const G2Fp2Ring& ring)
        {
            return PointDouble<G2Fp2Ring>(point, ring);
        }

        static BLST_FORCEINLINE ProjPoint dadd(
            const ProjPoint& left, const ProjPoint& right,
            const G2Fp2Ring& ring,
            const G2Fp2Ring::StandardElement* a4)
        {
            ProjPoint result = a4 == nullptr
                ? add_point(left, right, ring)
                : PointDaddA<G2Fp2Ring>(left, right, *a4, ring);
            result = result.ctime_replace_if(PointCoordIsExactZero(right.z), left);
            return result.ctime_replace_if(PointCoordIsExactZero(left.z), right);
        }

        static BLST_FORCEINLINE ProjPoint negate(
            const ProjPoint& point, const G2Fp2Ring& ring)
        {
            return {point.x, ring.standard_negate(point.y), point.z};
        }

        static BLST_FORCEINLINE ProjPoint zero(const G2Fp2Ring& ring)
        {
            return {ring.zero(), ring.one(), ring.zero()};
        }
    };
    using G2Projective = G2Curve::ProjPoint;
    using G2Affine = G2Curve::AffPoint;
    using G2GLS = GLSMult<FpRingCandidate, G2Curve>;

    static_assert(std::is_empty_v<G2Curve>);
} // namespace blst::avx::matrixnok
