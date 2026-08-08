// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
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
     * Selected proof candidate for G2 over the MatrixNoK base ring.
     *
     * The candidate reuses the shared complete homogeneous schedule and the
     * generic `FP2Ring`; this facade only fixes their parameter set. It
     * deliberately has no public dispatch entry until the correspondence proof
     * and release gates are complete.
     */
    using G2Fp2RingCandidate = FP2Ring<FpRingCandidate>;

    class G2CurveCandidate
    {
    public:
        using ProjPoint = ProjectivePoint<
            G2Fp2RingCandidate::StandardElement>;
        using AffPoint = AffinePoint<
            G2Fp2RingCandidate::StandardElement>;

        static INLINE ProjPoint add_point(
            const ProjPoint& left, const ProjPoint& right,
            const G2Fp2RingCandidate& ring)
        {
            return detail::add_point(left, right, ring);
        }

        static INLINE ProjPoint add_mixed_point(
            const ProjPoint& left, const AffPoint& right,
            const G2Fp2RingCandidate& ring)
        {
            return add_point(left, ProjPoint(right.x, right.y, ring.one()), ring)
                .ctime_replace_if(PointIsZero(right), left);
        }

        static INLINE ProjPoint double_point(
            const ProjPoint& point, const G2Fp2RingCandidate& ring)
        {
            return PointDouble<G2Fp2RingCandidate>(point, ring);
        }

        static INLINE ProjPoint dadd(
            const ProjPoint& left, const ProjPoint& right,
            const G2Fp2RingCandidate& ring,
            const G2Fp2RingCandidate::StandardElement* a4)
        {
            ProjPoint result = a4 == nullptr
                ? add_point(left, right, ring)
                : PointDaddA<G2Fp2RingCandidate>(left, right, *a4, ring);
            result = result.ctime_replace_if(PointCoordIsZero(right.z), left);
            return result.ctime_replace_if(PointCoordIsZero(left.z), right);
        }

        static INLINE ProjPoint negate(
            const ProjPoint& point, const G2Fp2RingCandidate& ring)
        {
            return {point.x, ring.standard_negate(point.y), point.z};
        }

        static INLINE ProjPoint zero(const G2Fp2RingCandidate& ring)
        {
            return {ring.zero(), ring.one(), ring.zero()};
        }
    };
    using G2ProjectiveCandidate = G2CurveCandidate::ProjPoint;
    using G2AffineCandidate = G2CurveCandidate::AffPoint;
    using G2GLSCandidate = GLSMult<FpRingCandidate, G2CurveCandidate>;

    static_assert(std::is_empty_v<G2CurveCandidate>);
} // namespace blst::avx::matrixnok
