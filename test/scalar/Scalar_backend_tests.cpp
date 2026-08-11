#include <gtest/gtest.h>

#include "blst/fields/Fp12.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/map_to_G.h"
#include "blst/scalar/curve/point.h"

namespace {

using blst::scalar::G1AffinePoint;
using blst::scalar::G1Point;
using blst::scalar::G2AffinePoint;
using blst::scalar::G2Point;
using blst::scalar::GT;
using blst::scalar::MillerResult;

TEST(ScalarBackend, G1IsogenyPathNeedsCofactorClear)
{
    G1Point isogenous = {};
    blst::scalar::map_to_isogenous<G1Point>(isogenous, G1AffinePoint::generator().x);

    G1Point mapped = {};
    blst::scalar::isogeny_map_to<G1Point>(mapped, isogenous);

    ASSERT_TRUE(G1Point::on_curve(mapped));

    G1AffinePoint affine = {};
    G1Point::to_affine(affine, mapped);
    EXPECT_FALSE(G1AffinePoint::in_subgroup(affine));

    G1Point cleared = {};
    blst::scalar::clear_cofactor<G1Point>(cleared, mapped);
    G1Point::to_affine(affine, cleared);
    EXPECT_TRUE(G1AffinePoint::in_subgroup(affine));
}

TEST(ScalarBackend, G2PrecomputedLinesMatchDirectMiller)
{
    const G1AffinePoint p = G1AffinePoint::generator();
    const G2AffinePoint q = G2AffinePoint::generator();

    MillerResult direct = {};
    GT::miller_loop(direct, &q, &p);

    blst::Fp6 lines[68] = {};
    GT::precompute_lines(lines, q);

    MillerResult via_lines = {};
    GT::miller_loop_lines(via_lines, lines, p);

    EXPECT_TRUE(blst::Fp12::equal(direct.data, via_lines.data));
}

TEST(ScalarBackend, PrecomputedLinesPreservePairingIdentity)
{
    const G1AffinePoint p = G1AffinePoint::generator();
    const G2AffinePoint q = G2AffinePoint::generator();
    const G1AffinePoint p_infinity = {};
    const G2AffinePoint q_infinity = {};

    blst::Fp6 finite_lines[68] = {};
    GT::precompute_lines(finite_lines, q);
    MillerResult p_infinity_result = {};
    GT::miller_loop_lines(p_infinity_result, finite_lines, p_infinity);
    EXPECT_TRUE(blst::Fp12::is_one(p_infinity_result.data));

    blst::Fp6 infinity_lines[68] = {};
    GT::precompute_lines(infinity_lines, q_infinity);
    MillerResult q_infinity_result = {};
    GT::miller_loop_lines(q_infinity_result, infinity_lines, p);
    EXPECT_TRUE(blst::Fp12::is_one(q_infinity_result.data));
}

}  // namespace
