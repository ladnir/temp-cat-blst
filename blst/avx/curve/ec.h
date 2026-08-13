/*
 * Copyright 2026 Simon Langowski
 * Copyright 2026 Category Labs
 * SPDX-License-Identifier: MIT
 *
 * The point carriers and PointDouble descend from the VROOM AVX backend at
 * commit d6ecd534. Category Labs added PointDaddA and retains only the
 * schedules selected by the fixed MatrixNoK BLS12-381 construction.
 */
#pragma once
#include "../../common/defines.h"
#include "blst/common/vect.h"
#include <type_traits>

namespace blst::avx {

/** Owning homogeneous projective coordinates `(X/Z, Y/Z)`.
 *
 * Every semantic `Z = 0` represents infinity. Selected construction uses the
 * canonical `(0, 1, 0)` representative. The type records no on-curve or
 * subgroup-validation state. It owns fixed-size storage and performs no
 * allocation, synchronization, or runtime dispatch.
 */
template<class RingElement>
class ProjectivePoint {
    public:
    RingElement x;
    RingElement y;
    RingElement z;
    
    ProjectivePoint() = default;
    ProjectivePoint(const RingElement &x_val, const RingElement &y_val, const RingElement &z_val)
        : x(x_val), y(y_val), z(z_val) {}

    // Return *this when selection is true and other when it is false. Selected
    // RingElement implementations lower this operation to mask selection.
    BLST_FORCEINLINE ProjectivePoint ctime_select(bool selection, const ProjectivePoint &other) const {
        return ProjectivePoint(x.ctime_select(selection, other.x), y.ctime_select(selection, other.y), z.ctime_select(selection, other.z));
    }

    // Clearer API: replace self with replacement if condition is true
    // ctime_select(selection, b) returns *this when selection is true, b when false
    // So we need to negate: if condition is true, we want replacement (not *this)
    BLST_FORCEINLINE ProjectivePoint ctime_replace_if(bool condition, const ProjectivePoint &replacement) const {
        return ctime_select(!condition, replacement);
    }

};

/** Owning affine coordinates `(X, Y)`.
 *
 * The exact all-zero object is the affine infinity sentinel. This is an exact
 * representation contract, not a semantic field-zero test over redundant RNS
 * representatives. The type records no validation state.
 */
template<class RingElement>
class AffinePoint {
    public:
    RingElement x;
    RingElement y;
    
    BLST_FORCEINLINE AffinePoint ctime_select(bool selection, const AffinePoint &other) const {
        return AffinePoint(x.ctime_select(selection, other.x), y.ctime_select(selection, other.y));
    }

    // Clearer API: replace self with replacement if condition is true
    // ctime_select(selection, b) returns *this when selection is true, b when false
    // So we need to negate: if condition is true, we want replacement (not *this)
    BLST_FORCEINLINE AffinePoint ctime_replace_if(bool condition, const AffinePoint &replacement) const {
        return ctime_select(!condition, replacement);
    }

};

template<class Element>
BLST_FORCEINLINE bool PointCoordIsExactZero(const Element& value);

template<class RingElement>
BLST_FORCEINLINE bool AffinePointIsExactIdentity(const AffinePoint<RingElement>& point);

/** Fixed-schedule doubling for `Y^2 Z = X^3 + b Z^3`.
 *
 * `P` must be on the curve represented by `ring`; subgroup membership is not
 * required. `ring.mul_3b` must multiply by that curve's `3b`. The selected
 * fixed MatrixNoK rings make every operation allocation-free and nonthrowing.
 * The returned point is on the same curve, including when `P` is infinity.
 *
 * This is Algorithm 9 of https://eprint.iacr.org/2015/1060.pdf. Preserve its
 * fixed, batched schedule: it is a selected hot-path and formal-proof owner.
 */
template<class Ring>
ProjectivePoint<typename Ring::StandardElement> PointDouble(const ProjectivePoint<typename Ring::StandardElement> &P,const Ring &ring) {
    auto py = ring.prep_left(P.y);
    auto pz = ring.prep_left(P.z);
    constexpr auto _2 = std::integral_constant<int, 2>{};
    constexpr auto _3 = std::integral_constant<int, 3>{};
    constexpr auto _8 = std::integral_constant<int, 8>{};

    auto yy = py * P.y;
    auto zz = pz * P.z;
    auto xy = py * P.x;
    auto yz = py * P.z;

    auto [yy_reduced, zz_reduced, xy_reduced, yz_reduced] = ring.batch_reduce_expand(yy, zz, xy, yz);

    // Add/sub operations
    auto _8yy = ring.prep_left(yy_reduced * _8);
    auto b3_zz = ring.prep(ring.mul_3b(zz_reduced));
    auto yy_plus_b3_zz = ring.prep(yy_reduced + b3_zz);
    auto b9_zz = b3_zz * _3;
    auto yy_minus_b9_zz = ring.prep_left(yy_reduced - b9_zz);
    auto two_xy = ring.prep(xy_reduced * _2);

    // Multiplication round 2
    auto y3_14 = yy_minus_b9_zz * yy_plus_b3_zz;
    auto x3_18 = yy_minus_b9_zz * two_xy;
    auto z3_10 = _8yy * yz_reduced;
    auto y3_15 = y3_14 + (_8yy * b3_zz);

    auto [x3, y3, z3] = ring.batch_reduce_expand(x3_18, y3_15, z3_10);
    return ProjectivePoint<typename Ring::StandardElement>(x3, y3, z3);
}

template<class Element>
BLST_FORCEINLINE bool PointCoordIsExactZero(const Element& value)
{
    return value.is_zero() != 0;
}

template<class RingElement>
BLST_FORCEINLINE bool AffinePointIsExactIdentity(const AffinePoint<RingElement>& point)
{
    return (static_cast<unsigned>(PointCoordIsExactZero(point.x)) &
        static_cast<unsigned>(PointCoordIsExactZero(point.y))) != 0;
}

// Normalize an internal expression before a semantic zero observation. This
// multiplication by one is intentional for MatrixNoK's redundant RNS values.
template<class Ring, class Element>
BLST_FORCEINLINE typename Ring::StandardElement NormalizePointCoord(const Element& value, const Ring& ring)
{
    return ring.batch_reduce_expand(ring.prep_left(value) * ring.one())[0];
}

template<class Ring, class IsZero>
BLST_FORCEINLINE bool ProjectiveEqual(
    const ProjectivePoint<typename Ring::StandardElement>& P,
    const ProjectivePoint<typename Ring::StandardElement>& Q,
    const Ring& ring,
    IsZero is_zero)
{
    const bool p_is_infinity = is_zero(P.z);
    const bool q_is_infinity = is_zero(Q.z);
    auto [u1, u2, s1, s2] = ring.batch_reduce_expand(
        ring.prep_left(P.x) * Q.z,
        ring.prep_left(Q.x) * P.z,
        ring.prep_left(P.y) * Q.z,
        ring.prep_left(Q.y) * P.z);
    const auto h = NormalizePointCoord(u2 - u1, ring);
    const auto r = NormalizePointCoord(s2 - s1, ring);
    return is_zero(h) & is_zero(r) &
        (p_is_infinity ^ q_is_infinity ^ true);
}

/** Fixed-schedule complete addition on a curve with nonzero `a4` coefficient.
 *
 * `P` and `Q` must be on the same curve, and `a4` must be that curve's
 * coefficient in the formula consumed by this schedule. Subgroup membership
 * is not required. Infinity selection is owned by the calling curve facade.
 * Selected rings make the operation allocation-free and nonthrowing.
 *
 * Preserve this fixed, batched schedule. Point coordinates may be secret even
 * though the caller's choice between the a4 and a=0 schedules is public.
 */
template<class Ring>
ProjectivePoint<typename Ring::StandardElement> PointDaddA(
    const ProjectivePoint<typename Ring::StandardElement>& P,
    const ProjectivePoint<typename Ring::StandardElement>& Q,
    const typename Ring::StandardElement& a4,
    const Ring& ring)
{
    constexpr auto _2 = std::integral_constant<int, 2>{};

    auto u1_wide = ring.prep_left(P.x) * Q.z;
    auto u2_wide = ring.prep_left(Q.x) * P.z;
    auto s1_wide = ring.prep_left(P.y) * Q.z;
    auto s2_wide = ring.prep_left(Q.y) * P.z;
    auto z_wide = ring.prep_left(P.z) * Q.z;

    auto [u1, u2, s1, s2, z] = ring.batch_reduce_expand(
        u1_wide,
        u2_wide,
        s1_wide,
        s2_wide,
        z_wide);

    auto t = NormalizePointCoord(u1 + u2, ring);
    auto m = NormalizePointCoord(s1 + s2, ring);

    auto tt_wide = ring.prep_left(t) * t;
    auto u1u2_wide = ring.prep_left(u1) * u2;
    auto zz_wide = ring.prep_left(z) * z;
    auto f_wide = ring.prep_left(z) * m;

    auto [tt, u1u2, zz, f] = ring.batch_reduce_expand(
        tt_wide,
        u1u2_wide,
        zz_wide,
        f_wide);

    auto a_zz = ring.batch_reduce_expand(ring.prep_left(a4) * zz)[0];
    auto r = NormalizePointCoord(tt - u1u2 + a_zz, ring);
    auto l = ring.batch_reduce_expand(ring.prep_left(m) * f)[0];

    auto g_wide = ring.prep_left(t) * l;
    auto r2_wide = ring.prep_left(r) * r;
    auto f2_wide = ring.prep_left(f) * f;
    auto ll_wide = ring.prep_left(l) * l;

    auto [g, r2, f2, ll] = ring.batch_reduce_expand(
        g_wide,
        r2_wide,
        f2_wide,
        ll_wide);

    auto w = NormalizePointCoord(r2 - g, ring);
    auto g_minus_2w = NormalizePointCoord(g - w * _2, ring);
    auto two_f = NormalizePointCoord(f * _2, ring);

    auto [x3, y3_term, z3] = ring.batch_reduce_expand(
        ring.prep_left(two_f) * w,
        ring.prep_left(r) * g_minus_2w,
        ring.prep_left(two_f) * f2);

    auto y3 = NormalizePointCoord(y3_term - ll, ring);
    return ProjectivePoint<typename Ring::StandardElement>(x3, y3, z3);
}

}  // namespace blst::avx
