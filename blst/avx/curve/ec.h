#pragma once
#include "../common/preprocess.h"
#include "blst/common/vect.h"
#include <type_traits>

namespace blst::avx {

template<class RingElement>
class ProjectivePoint {
    public:
    RingElement x;
    RingElement y;
    RingElement z;
    
    ProjectivePoint() = default;
    ProjectivePoint(const RingElement &x_val, const RingElement &y_val, const RingElement &z_val)
        : x(x_val), y(y_val), z(z_val) {}

    INLINE ProjectivePoint ctime_select(bool selection, const ProjectivePoint &other) const {
        return ProjectivePoint(x.ctime_select(selection, other.x), y.ctime_select(selection, other.y), z.ctime_select(selection, other.z));
    }

    // Clearer API: replace self with replacement if condition is true
    // ctime_select(selection, b) returns *this when selection is true, b when false
    // So we need to negate: if condition is true, we want replacement (not *this)
    INLINE ProjectivePoint ctime_replace_if(bool condition, const ProjectivePoint &replacement) const {
        return ctime_select(!condition, replacement);
    }

    // Clearer API: replace self with replacement if condition is false
    INLINE ProjectivePoint ctime_replace_unless(bool condition, const ProjectivePoint &replacement) const {
        return ctime_select(condition, replacement);
    }
};

template<class RingElement>
class AffinePoint {
    public:
    RingElement x;
    RingElement y;
    
    INLINE AffinePoint ctime_select(bool selection, const AffinePoint &other) const {
        return AffinePoint(x.ctime_select(selection, other.x), y.ctime_select(selection, other.y));
    }

    // Clearer API: replace self with replacement if condition is true
    // ctime_select(selection, b) returns *this when selection is true, b when false
    // So we need to negate: if condition is true, we want replacement (not *this)
    INLINE AffinePoint ctime_replace_if(bool condition, const AffinePoint &replacement) const {
        return ctime_select(!condition, replacement);
    }

    // Clearer API: replace self with replacement if condition is false
    INLINE AffinePoint ctime_replace_unless(bool condition, const AffinePoint &replacement) const {
        return ctime_select(condition, replacement);
    }
};

template<class Ring>
class FP2Ring;

template<class Element>
INLINE bool PointCoordIsZero(const Element& value);

template<class RingElement>
INLINE bool PointIsZero(const AffinePoint<RingElement>& point);

// ~5% faster, but doesn't work with new type system yet because it needs to preserve rns bounds across expansion.
/*
template<class Ring, class OffsetElement>
ProjectivePoint<Ring::StandardElement> PointAdd(const ProjectivePoint<Ring::StandardElement> &P, const ProjectivePoint<Ring::StandardElement> &Q, const Ring &ring, const OffsetElement &offset) {
    auto t3_4 = P.x + P.y;
    auto t3_14 = P.x + P.z;
    auto t4_9 = P.y + P.z;
    auto t4_5 = Q.x + Q.y;
    auto x3_10 = Q.y + Q.z;
    auto y3_15 = Q.x + Q.z;

    auto t0_1 = ring.prep_left(P.x) * Q.x;
    auto t1_2 = ring.prep_left(P.y) * Q.y;
    auto t2_3 = ring.prep_left(P.z) * Q.z;
    auto t3_6 = ring.prep_left(t3_4) * ring.prep(t4_5);
    auto t5_11 = ring.prep_left(t4_9) * ring.prep(x3_10);
    auto x3_16 = ring.prep_left(t3_14) * ring.prep(y3_15);

    auto [t0_1_reduced, t1_2_reduced, t2_3_reduced, t3_6_reduced, t5_11_reduced, x3_16_reduced] = ring.batch_reduce(t0_1, t1_2, t2_3, t3_6, t5_11, x3_16);
    auto t4_7 = t0_1_reduced + t1_2_reduced;
    auto x3_12 = t1_2_reduced + t2_3_reduced;
    auto y3_17 = t0_1_reduced + t2_3_reduced;
    auto y3_18 = x3_16_reduced - y3_17;
    auto x3_19 = t0_1_reduced + t0_1_reduced;
    auto t2_21 = ring.mul_3b(t2_3_reduced);

    auto t3_8_s = t3_6_reduced - t4_7;
    auto t4_13_s = t5_11_reduced - x3_12;
    auto t0_20_s = t0_1_reduced + x3_19;
    auto z3_22_s = t1_2_reduced + t2_21;
    auto t1_23_s = t1_2_reduced - t2_21;
    auto y3_24_s = ring.mul_3b(y3_18);

    auto [t4_13, t3_8, y3_24, t1_23, t0_20, z3_22] = ring.batch_expand(t4_13_s, t3_8_s, y3_24_s, t1_23_s, t0_20_s, z3_22_s);
    
    auto t4_13l = ring.prep_left(t4_13);
    auto t1_23l = ring.prep_left(t1_23);
    auto t0_20l = ring.prep_left(t0_20);
    auto x3_25 = t4_13l * ring.negate(y3_24);
    auto t2_26 = t1_23l * t3_8;
    auto y3_28 = t0_20l * y3_24;
    auto t1_29 = t1_23l * z3_22;
    auto t0_31 = t0_20l * t3_8;
    auto z3_32 = t4_13l * z3_22;
    auto x3_27_wide = offset + t2_26 + x3_25;
    auto y3_30_wide = offset + t1_29 + y3_28;
    auto z3_33_wide = offset + z3_32 + t0_31;
    
    auto [x3_27, y3_30, z3_33] = ring.batch_reduce_expand(x3_27_wide, y3_30_wide, z3_33_wide);
    return ProjectivePoint<Ring::StandardElement>(x3_27, y3_30, z3_33);
}
*/

// https://eprint.iacr.org/2015/1060.pdf Algorithm 7
template<class Ring>
ProjectivePoint<typename Ring::StandardElement> PointAdd(const ProjectivePoint<typename Ring::StandardElement> &P, const ProjectivePoint<typename Ring::StandardElement> &Q, const Ring &ring) {
    auto t3_4 = P.x + P.y;
    auto t3_14 = P.x + P.z;
    auto t4_9 = P.y + P.z;
    auto t4_5 = Q.x + Q.y;
    auto x3_10 = Q.y + Q.z;
    auto y3_15 = Q.x + Q.z;

    auto t0_1_wide = ring.prep_left(P.x) * Q.x;
    auto t1_2_wide = ring.prep_left(P.y) * Q.y;
    auto t2_3_wide = ring.prep_left(P.z) * Q.z;
    auto t3_6_wide = ring.prep_left(t3_4) * ring.prep(t4_5);
    auto t5_11_wide = ring.prep_left(t4_9) * ring.prep(x3_10);
    auto x3_16_wide = ring.prep_left(t3_14) * ring.prep(y3_15);

    auto [t0_1, t1_2, t2_3, t3_6, t5_11, x3_16] = ring.batch_reduce(t0_1_wide, t1_2_wide, t2_3_wide, t3_6_wide, t5_11_wide, x3_16_wide);
    auto t4_7 = t0_1 + t1_2;
    auto x3_12 = t1_2 + t2_3;
    auto y3_17 = t0_1 + t2_3;
    auto y3_18 = x3_16 - y3_17;
    auto x3_19 = t0_1 + t0_1;
    auto t2_21 = ring.mul_3b(t2_3);

    auto t3_8_s = ring.prep(t3_6 - t4_7);
    auto y3_24_s = ring.prep(ring.mul_3b(y3_18));
    auto z3_22_s = ring.prep(t1_2 + t2_21);
    auto t1_23_s = ring.prep(t1_2 - t2_21);
    auto t4_13_s = ring.prep(t5_11 - x3_12);
    auto t0_20_s = ring.prep(t0_1 + x3_19);

    auto [t3_8, y3_24, z3_22, t1_23, t4_13, t0_20] = ring.batch_expand(t3_8_s, y3_24_s, z3_22_s, t1_23_s, t4_13_s, t0_20_s);

    auto t4_13l = ring.prep_left(t4_13);
    auto t1_23l = ring.prep_left(t1_23);
    auto t0_20l = ring.prep_left(t0_20);
    
    auto x3_25 = t4_13l * ring.negate(y3_24);
    auto t2_26 = t1_23l * t3_8;
    auto y3_28 = t0_20l * y3_24;
    auto t1_29 = t1_23l * z3_22;
    auto t0_31 = t0_20l * t3_8;
    auto z3_32 = t4_13l * z3_22;
    auto x3_27_wide =  t2_26 + x3_25;
    auto y3_30_wide =  t1_29 + y3_28;
    auto z3_33_wide =  z3_32 + t0_31;
    
    auto [x3_27, y3_30, z3_33] = ring.batch_reduce_expand(x3_27_wide, y3_30_wide, z3_33_wide);
    return ProjectivePoint<typename Ring::StandardElement>(x3_27, y3_30, z3_33);
}

// Not enough time to figure out Algorithm 8; this saves a mult and you could save more additions (by also using small form)
template<class Ring>
ProjectivePoint<typename Ring::StandardElement> PointMixedAdd(const ProjectivePoint<typename Ring::StandardElement> &P, const AffinePoint<typename Ring::StandardElement> &Q, const Ring &ring) {
    auto t3_4 = P.x + P.y;
    auto t3_14 = P.x + P.z;
    auto t4_9 = P.y + P.z;
    auto t4_5 = Q.x + Q.y;
    auto one = ring.one();
    auto x3_10 = Q.y + one;
    auto y3_15 = Q.x + one;

    auto t0_1_wide = ring.prep_left(P.x) * Q.x;
    auto t1_2_wide = ring.prep_left(P.y) * Q.y;
    // auto t2_3_wide = ring.prep_left(P.z) * Q.z;
    auto t3_6_wide = ring.prep_left(t3_4) * ring.prep(t4_5);
    auto t5_11_wide = ring.prep_left(t4_9) * ring.prep(x3_10);
    auto x3_16_wide = ring.prep_left(t3_14) * ring.prep(y3_15);

    auto [t0_1, t1_2, t3_6, t5_11, x3_16] = ring.batch_reduce_expand(t0_1_wide, t1_2_wide, t3_6_wide, t5_11_wide, x3_16_wide);
    auto t2_3 = P.z;
    auto t4_7 = t0_1 + t1_2;
    auto x3_12 = t1_2 + t2_3;
    auto y3_17 = t0_1 + t2_3;
    auto y3_18 = x3_16 - y3_17;
    auto x3_19 = t0_1 + t0_1;
    auto t2_21 = ring.prep(ring.mul_3b(t2_3));

    auto t3_8 = ring.prep(t3_6 - t4_7);
    auto y3_24 = ring.prep(ring.mul_3b(y3_18));
    auto z3_22 = ring.prep(t1_2 + t2_21);
    auto t1_23 = ring.prep_left(t1_2 - t2_21);
    auto t4_13 = ring.prep_left(t5_11 - x3_12);
    auto t0_20 = ring.prep_left(t0_1 + x3_19);
    
    auto x3_25 = t4_13 * ring.negate(y3_24);
    auto t2_26 = t1_23 * t3_8;
    auto y3_28 = t0_20 * y3_24;
    auto t1_29 = t1_23 * z3_22;
    auto t0_31 = t0_20 * t3_8;
    auto z3_32 = t4_13 * z3_22;
    auto x3_27_wide =  t2_26 + x3_25;
    auto y3_30_wide =  t1_29 + y3_28;
    auto z3_33_wide =  z3_32 + t0_31;
    
    auto [x3_27, y3_30, z3_33] = ring.batch_reduce_expand(x3_27_wide, y3_30_wide, z3_33_wide);
    return ProjectivePoint<typename Ring::StandardElement>(x3_27, y3_30, z3_33)
        .ctime_replace_if(PointIsZero(Q), P);
}

template<class Ring>
ProjectivePoint<typename Ring::StandardElement> PointMixedAddUnexpanded(const ProjectivePoint<typename Ring::StandardElement> &P, const AffinePoint<typename Ring::StandardElement> &Q, const Ring &ring) {
    auto t3_4 = P.x + P.y;
    auto t3_14 = P.x + P.z;
    auto t4_9 = P.y + P.z;
    auto t4_5 = Q.x + Q.y;
    auto one = ring.one();
    auto x3_10 = Q.y + one;
    auto y3_15 = Q.x + one;

    auto t0_1_wide = ring.prep_left(P.x) * Q.x;
    auto t1_2_wide = ring.prep_left(P.y) * Q.y;
    auto t3_6_wide = ring.prep_left(t3_4) * ring.prep(t4_5);
    auto t5_11_wide = ring.prep_left(t4_9) * ring.prep(x3_10);
    auto x3_16_wide = ring.prep_left(t3_14) * ring.prep(y3_15);

    auto [t0_1, t1_2, t3_6, t5_11, x3_16] =
        ring.batch_reduce(t0_1_wide, t1_2_wide, t3_6_wide, t5_11_wide, x3_16_wide);
    auto t2_3 = ring.small_m2(P.z);
    auto t4_7 = t0_1 + t1_2;
    auto x3_12 = t1_2 + t2_3;
    auto y3_17 = t0_1 + t2_3;
    auto y3_18 = x3_16 - y3_17;
    auto x3_19 = t0_1 + t0_1;
    auto t2_21 = ring.mul_3b(t2_3);

    auto t3_8_s = ring.prep(t3_6 - t4_7);
    auto y3_24_s = ring.prep(ring.mul_3b(y3_18));
    auto z3_22_s = ring.prep(t1_2 + t2_21);
    auto t1_23_s = ring.prep(t1_2 - t2_21);
    auto t4_13_s = ring.prep(t5_11 - x3_12);
    auto t0_20_s = ring.prep(t0_1 + x3_19);

    auto [t3_8, y3_24, z3_22, t1_23, t4_13, t0_20] =
        ring.batch_expand(t3_8_s, y3_24_s, z3_22_s, t1_23_s, t4_13_s, t0_20_s);

    auto t4_13l = ring.prep_left(t4_13);
    auto t1_23l = ring.prep_left(t1_23);
    auto t0_20l = ring.prep_left(t0_20);
    auto x3_25 = t4_13l * ring.negate(y3_24);
    auto t2_26 = t1_23l * t3_8;
    auto y3_28 = t0_20l * y3_24;
    auto t1_29 = t1_23l * z3_22;
    auto t0_31 = t0_20l * t3_8;
    auto z3_32 = t4_13l * z3_22;
    auto x3_27_wide = t2_26 + x3_25;
    auto y3_30_wide = t1_29 + y3_28;
    auto z3_33_wide = z3_32 + t0_31;

    auto [x3_27, y3_30, z3_33] = ring.batch_reduce_expand(x3_27_wide, y3_30_wide, z3_33_wide);
    return ProjectivePoint<typename Ring::StandardElement>(x3_27, y3_30, z3_33)
        .ctime_replace_if(PointIsZero(Q), P);
}

// https://eprint.iacr.org/2015/1060.pdf Algorithm 9
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
INLINE bool PointCoordIsZero(const Element& value)
{
    return value.is_zero() != 0;
}

template<class RingElement>
INLINE bool PointIsZero(const AffinePoint<RingElement>& point)
{
    return PointCoordIsZero(point.x) & PointCoordIsZero(point.y);
}

template<class Ring, class Element>
INLINE typename Ring::StandardElement NormalizePointCoord(const Element& value, const Ring& ring)
{
    return ring.batch_reduce_expand(ring.prep_left(value) * ring.one())[0];
}

template<class Ring>
INLINE bool ProjectiveEqual(
    const ProjectivePoint<typename Ring::StandardElement>& P,
    const ProjectivePoint<typename Ring::StandardElement>& Q,
    const Ring& ring)
{
    auto [u1, u2, s1, s2] = ring.batch_reduce_expand(
        ring.prep_left(P.x) * Q.z,
        ring.prep_left(Q.x) * P.z,
        ring.prep_left(P.y) * Q.z,
        ring.prep_left(Q.y) * P.z);
    const auto h = u2 - u1;
    const auto r = s2 - s1;
    return PointCoordIsZero(h) && PointCoordIsZero(r);
}

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


template<class Ring>
class G1 {
    public:
    using ProjPoint = ProjectivePoint<typename Ring::StandardElement>;
    using AffPoint = AffinePoint<typename Ring::StandardElement>;

    ProjPoint add_point(const ProjPoint &P, const ProjPoint &Q, const Ring &ring) const {
        return PointAdd<Ring>(P, Q, ring);
    }

    ProjPoint add_mixed_point(const ProjPoint &P, const AffPoint &Q, const Ring &ring) const {
        return PointMixedAddUnexpanded<Ring>(P, Q, ring);
    }

    ProjPoint double_point(const ProjPoint &P, const Ring &ring) const {
        return PointDouble<Ring>(P, ring);
    }

    ProjPoint dadd(const ProjPoint& P, const ProjPoint& Q, const Ring& ring, const typename Ring::StandardElement* a4) const {
        ProjPoint result = {};
        if (a4 != nullptr) {
            result = PointDaddA<Ring>(P, Q, *a4, ring);
        } else {
            result = PointAdd<Ring>(P, Q, ring);
        }
        result = result.ctime_replace_if(PointCoordIsZero(Q.z), P);
        return result.ctime_replace_if(PointCoordIsZero(P.z), Q);
    }

    INLINE ProjPoint negate(const ProjPoint &P, const Ring &ring) const {
        return ProjPoint(P.x, ring.standard_negate(P.y), P.z);
    }

    INLINE ProjPoint zero(const Ring &ring) const {
        return ProjPoint(ring.zero(), ring.one(), ring.zero());
    }

};

template<class Ring>
class G2 {
    public:
    using FP2RingType = FP2Ring<Ring>;
    using ProjPoint = ProjectivePoint<typename FP2RingType::StandardElement>;
    using AffPoint = AffinePoint<typename FP2RingType::StandardElement>;

    ProjPoint add_point(const ProjPoint &P, const ProjPoint &Q, const FP2RingType &ring) const {
        return PointAdd<FP2RingType>(P, Q, ring);
    }

    ProjPoint add_mixed_point(const ProjPoint &P, const AffPoint &Q, const FP2RingType &ring) const {
        return PointMixedAdd<FP2RingType>(P, Q, ring);
    }

    ProjPoint double_point(const ProjPoint &P, const FP2RingType &ring) const {
        return PointDouble<FP2RingType>(P, ring);
    }

    ProjPoint dadd(const ProjPoint& P, const ProjPoint& Q, const FP2RingType& ring, const typename FP2RingType::StandardElement* a4) const {
        ProjPoint result = {};
        if (a4 != nullptr) {
            result = PointDaddA<FP2RingType>(P, Q, *a4, ring);
        } else {
            result = PointAdd<FP2RingType>(P, Q, ring);
        }
        result = result.ctime_replace_if(PointCoordIsZero(Q.z), P);
        return result.ctime_replace_if(PointCoordIsZero(P.z), Q);
    }

    INLINE ProjPoint negate(const ProjPoint &P, const FP2RingType &ring) const {
        return ProjPoint(P.x, ring.standard_negate(P.y), P.z);
    }

    INLINE ProjPoint zero(const FP2RingType &ring) const {
        return ProjPoint(ring.zero(), ring.one(), ring.zero());
    }
};

}  // namespace blst::avx
