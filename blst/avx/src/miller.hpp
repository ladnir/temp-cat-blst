#pragma once
#include "ring_element.hpp"
#include "fp2.hpp"
#include "fp12.hpp"
#include "ec.hpp"

template<class Ring>
class Miller {
    public:
    // Reduce the amount of typename nonsense
    using FPStandardElement = typename Ring::StandardElement;
    using FP2RingStandardElement = typename FP2Ring<Ring>::StandardElement;
    using FP12RingStandardElement = typename FP12Ring<Ring>::StandardElement;

    // An element that is 52 bits, but the RNS bound is not [0, 1]
    template<int lower_bound, int upper_bound>
    using ReducedElement = ExpandedElement<typename Ring::StandardBoundedElement, Bounds<lower_bound, upper_bound>>;
    private:
    // E1 affine point
    const FPStandardElement E1X;
    const FPStandardElement E1Y;
    // Precompute scalar multiples
    const ReducedElement<0, 3> three_E1X;
    const ReducedElement<0, 2> minus2_E1Y;

    using FP2NStandardElement = FP2N<FPStandardElement, FPStandardElement, FPStandardElement>;

    // E2 affine point
    // Precompute negation for FP2 products
    const FP2NStandardElement E2X;
    const FP2NStandardElement E2Y;

    // Shorthand for integral constants
	static constexpr auto _2 = std::integral_constant<int, 2>{};
	static constexpr auto _3 = std::integral_constant<int, 3>{};
    static constexpr auto _8 = std::integral_constant<int, 8>{};

    FP2NStandardElement prep_left_helper(const FP2<FPStandardElement, FPStandardElement> &point, const Ring &ring) {
        auto r_x = point.x();
        auto r_y = point.y();
        auto minus_y = ring.standard_negate(r_y);
        return FP2NStandardElement(r_x, r_y, minus_y);
    }

    // This function does NOT check for identity or infinity points
    public:
    Miller(const AffinePoint<FPStandardElement> &E1_point, const AffinePoint<FP2RingStandardElement> &E2_point, const Ring &ring)
        : 
          E1X(E1_point.x), 
          E1Y(E1_point.y),
          three_E1X(ring.prep_standard(E1_point.x * _3)),
          minus2_E1Y(ring.prep_standard(ring.standard_negate(E1_point.y) * _2)),
          E2X(prep_left_helper(E2_point.x, ring)),
          E2Y(prep_left_helper(E2_point.y, ring)) {
    }

    // returns point, mutates f
    INLINE ProjectivePoint<FP2RingStandardElement> double_step(const ProjectivePoint<FP2RingStandardElement> &P, FP12RingStandardElement &f, const Ring &ring, const FP2Ring<Ring> &fp2ring, const FP12Ring<Ring> &fp12ring) const {
        // https://eprint.iacr.org/2015/1060.pdf Algorithm 9
        // Merged with line computation and update:
        // line_x = x^2 * (3Qx)
        // line_y = yz * (-2Qy),
        // line_c = 3b*z^2 - y^2

        // Multiplication round 1
        auto px = P.x.prep_left(ring);
        auto py = P.y.prep_left(ring);
        auto pz = P.z.prep_left(ring);

        auto yy = py.square();
        auto zz = pz.square();
        auto xx = px.square();
        auto xy = px * py.to_fp2();
        auto yz = py * pz.to_fp2();

        auto [xx_reduced, yy_reduced, zz_reduced, xy_reduced, yz_reduced] = fp2ring.batch_reduce_expand(xx, yy, zz, xy, yz);

        // Add/sub operations
        auto _8yy = yy_reduced.mul_scalar(_8).prep_left(ring);
        auto b3_zz = fp2ring.mul_3b(zz_reduced).prep(ring);
        auto line_c = (b3_zz - yy_reduced).prep(ring);
        auto yy_plus_b3_zz = (yy_reduced + b3_zz).prep(ring);
        auto b9_zz = b3_zz * _3;
        auto yy_minus_b9_zz = (yy_reduced - b9_zz).prep_left(ring);
        auto two_xy = (xy_reduced * _2).prep(ring);

        // Multiplication round 2
        auto xq_term = xx_reduced.mul_fp(three_E1X);
        auto yq_term = yz_reduced.mul_fp(minus2_E1Y);
        auto y3_14 = yy_minus_b9_zz * yy_plus_b3_zz;
        auto x3_18 = yy_minus_b9_zz * two_xy;
        auto z3_10 = _8yy * yz_reduced;
        auto y3_15 = y3_14.mul_accumulate(_8yy, b3_zz);

        auto [line_x, line_y, x3, y3, z3] = fp2ring.batch_reduce_expand(xq_term, yq_term, x3_18, y3_15, z3_10);
        // if constexpr (first_loop) {
        //     f = FP12RingStandardElement {line_c.x(), line_c.y(), line_x.x(), line_x.y(), ring.zero(), ring.zero(), ring.zero(), ring.zero(), line_y.x(), line_y.y(), ring.zero(), ring.zero()};
        // } else {
            fp12ring.flat_mulxy00z0(f, f, line_c, line_x, line_y, ring);
       // }
        return ProjectivePoint<typename FP2Ring<Ring>::StandardElement>(x3, y3, z3);
    }

    // Returns point, mutates f.
    INLINE ProjectivePoint<FP2RingStandardElement> add_step(const ProjectivePoint<FP2RingStandardElement> &P, FP12RingStandardElement &f, const Ring &ring, const FP2Ring<Ring> &fp2ring, const FP12Ring<Ring> &fp12ring) const {
        // https://eprint.iacr.org/2015/1060.pdf Algorithm 8
        // Merged with line computation and update:
        // Line:
        // line_x = (-y1 + y2z1) * Qx
        // line_y = (x1 - x2z1) * Qy,
        // line_c = -y2x1 + x2y1
        // First multiplications
        auto x2x1 = E2X * P.x;
        auto x2y1 = E2X * P.y;
        auto x2z1 = E2X * P.z;
        auto y2x1 = E2Y * P.x;
        auto y2y1 = E2Y * P.y;
        auto y2z1 = E2Y * P.z;

        auto [x2x1_reduced, x2y1_reduced, x2z1_reduced, y2x1_reduced, y2y1_reduced, y2z1_reduced] = fp2ring.batch_reduce_expand(x2x1, x2y1, x2z1, y2x1, y2y1, y2z1);

        // Add/sub operations
        auto line_c = (x2y1_reduced - y2x1_reduced).prep(ring);
        auto xq_term_coef = (y2z1_reduced - P.y).prep(ring);
        auto yq_term_coef = (P.x - x2z1_reduced).prep(ring);
        auto x2y1_plus_y2x1 = (x2y1_reduced + y2x1_reduced).prep_left(ring);
        auto y2z1_plus_y1 = (y2z1_reduced + P.y).prep(ring);
        auto x2z1_plus_x1 = x2z1_reduced + P.x;
        auto _3x2x1 = (x2x1_reduced * _3).prep(ring);
        auto b3z1 = fp2ring.mul_3b(P.z).prep(ring); // # reduce mask -> 2
        auto y2y1_plus_b3z1 = (y2y1_reduced + b3z1).prep_left(ring); // # 3, negate -> 4
        auto y2y1_minus_b3z1 = (y2y1_reduced - b3z1).prep(ring); // # negate -> 4
        auto b3_x2z1_plus_x1 = fp2ring.mul_3b(x2z1_plus_x1).prep_left(ring); // # reduce
  
        // Second round of multiplications
        auto t2_19 = x2y1_plus_y2x1 * y2y1_minus_b3z1; // (x2y1 + y2x1) * (y2y1 - b3z1)
        auto x3_20 = t2_19.mul_sub_accumulate(b3_x2z1_plus_x1, y2z1_plus_y1, ring);
        auto y3_21 = b3_x2z1_plus_x1 * _3x2x1;
        auto y3_23 = y3_21.mul_accumulate(y2y1_plus_b3z1, y2y1_minus_b3z1);
        auto z3_25 = y2y1_plus_b3z1 * y2z1_plus_y1;
        auto z3_26 = z3_25.mul_accumulate(x2y1_plus_y2x1, _3x2x1);
        
        auto xq_term = xq_term_coef.mul_fp(E1X);
        auto yq_term = yq_term_coef.mul_fp(E1Y);

        auto [x3, y3, z3, line_x, line_y] = fp2ring.batch_reduce_expand(x3_20, y3_23, z3_26, xq_term, yq_term);

        fp12ring.flat_mulxy00z0(f, f, line_c, line_x, line_y, ring);
        return ProjectivePoint<typename FP2Ring<Ring>::StandardElement>(x3, y3, z3);
    }

    // Returns point, mutates f.
    INLINE ProjectivePoint<FP2RingStandardElement> add_dbl(const ProjectivePoint<FP2RingStandardElement> &P, FP12RingStandardElement &f, const Ring &ring, const FP2Ring<Ring> &fp2ring, const FP12Ring<Ring> &fp12ring, int iterations) const {
        auto result = add_step(P, f, ring, fp2ring, fp12ring);
#if defined(__clang__)
        #pragma clang loop unroll(disable)
#endif
        for (int i = 0; i < iterations; i++) {
            f = fp12ring.fp12_flat_sqr(f, ring);
            result = double_step(result, f, ring, fp2ring, fp12ring);
        }
        return result;
    }
    
    inline FP12RingStandardElement miller_loop(const Ring &ring, const FP2Ring<Ring> &fp2ring, const FP12Ring<Ring> &fp12ring) const {
        ProjectivePoint<FP2RingStandardElement> P(
            E2X.to_fp2(), 
            E2Y.to_fp2(), 
            FP2Ring<Ring>::one(ring)
        );
        auto f = fp12ring.one(ring);
        // Could optimize first iteration like blst
        // or use goofy goto to avoid code duplication if inlining all
        P = double_step(P, f, ring, fp2ring, fp12ring);
        // P = double_step<true>(P, f, ring, fp2ring, fp12ring);
        std::array<int, 5> iterations = {2, 3, 9, 32, 16};
#if defined(__clang__)
        #pragma clang loop unroll(disable)
#endif
        for (size_t i = 0; i < iterations.size(); i++) {
            P = add_dbl(P, f, ring, fp2ring, fp12ring, iterations[i]);
        }
        f = fp12ring.fp12_flat_conjugate(f, ring);
        return f;
    }
};