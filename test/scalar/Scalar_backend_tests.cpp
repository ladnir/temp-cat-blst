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

TEST(ScalarBackend, G1IsogenyPathNeedsCofactorClear)
{
    G1Point isogenous = {};
    blst::scalar::map_to_isogenous<G1Point>(isogenous, G1AffinePoint::generator().X);

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

    GT direct = {};
    GT::miller_loop(direct, &q, &p);

    blst::Fp6 lines[68] = {};
    GT::precompute_lines(lines, q);

    GT via_lines = {};
    GT::miller_loop_lines(via_lines, lines, p);

    EXPECT_TRUE(GT::equal(direct, via_lines));
}

}  // namespace
