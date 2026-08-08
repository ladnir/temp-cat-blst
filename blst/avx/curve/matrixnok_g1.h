// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "ec.h"
#include "matrixnok_curve.h"
#include "scalar_mult.h"
#include "../rns/matrixnok_fp.h"

#include <type_traits>

namespace blst::avx::matrixnok
{
    // Selected proof candidate for G1 over the fixed MatrixNoK ring. It remains
    // deliberately absent from public dispatch until the correspondence proof
    // and release gates are complete.
    class G1CurveCandidate
    {
    public:
        using ProjPoint = ProjectivePoint<FpRingCandidate::StandardElement>;
        using AffPoint = AffinePoint<FpRingCandidate::StandardElement>;

        static INLINE ProjPoint add_point(
            const ProjPoint& left, const ProjPoint& right,
            const FpRingCandidate& ring)
        {
            return detail::add_point(left, right, ring);
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
