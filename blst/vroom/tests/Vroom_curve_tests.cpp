#include <gtest/gtest.h>

#include "blst/vroom/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::VroomCurveFp2Ring;
using blst_vroom_test::VroomCurveG1;
using blst_vroom_test::VroomCurveG2;
using blst_vroom_test::VroomCurveRing;

TEST(VroomCurve, G1AddMixedAddAndDoubleMatchScalarBlst)
{
    const VroomCurveRing ring(blst_vroom_test::modulus_bigint());
    const VroomCurveG1 g1 = {};

    const blst::G1AffinePoint scalar_generator_affine = blst_vroom_test::make_scalar_g1_generator_affine();
    const blst::G1AffinePoint scalar_doubled_affine = blst_vroom_test::make_scalar_g1_doubled_affine();

    const auto vroom_generator_affine = blst_vroom_test::make_vroom_g1_affine(ring, scalar_generator_affine);
    const auto vroom_doubled_affine = blst_vroom_test::make_vroom_g1_affine(ring, scalar_doubled_affine);
    const auto vroom_generator_projective = blst_vroom_test::make_vroom_curve_g1_projective(ring, vroom_generator_affine);
    const auto vroom_doubled_projective = blst_vroom_test::make_vroom_curve_g1_projective(ring, vroom_doubled_affine);

    blst::G1Point scalar_generator_projective = {};
    blst::G1Point::from_affine(scalar_generator_projective, scalar_generator_affine);
    blst::G1Point scalar_doubled_projective = {};
    blst::G1Point::from_affine(scalar_doubled_projective, scalar_doubled_affine);

    blst::G1Point scalar_add_projective = {};
    blst::G1Point::add(scalar_add_projective, scalar_generator_projective, scalar_doubled_projective);
    blst::G1AffinePoint scalar_add_affine = {};
    blst::G1Point::to_affine(scalar_add_affine, scalar_add_projective);

    blst::G1Point scalar_mixed_add_projective = {};
    blst::G1Point::add(scalar_mixed_add_projective, scalar_generator_projective, scalar_doubled_affine);
    blst::G1AffinePoint scalar_mixed_add_affine = {};
    blst::G1Point::to_affine(scalar_mixed_add_affine, scalar_mixed_add_projective);

    blst::G1Point scalar_double_projective = {};
    blst::G1Point::doubling(scalar_double_projective, scalar_generator_projective);
    blst::G1AffinePoint scalar_double_affine = {};
    blst::G1Point::to_affine(scalar_double_affine, scalar_double_projective);

    const auto vroom_add = g1.add_point(vroom_generator_projective, vroom_doubled_projective, ring);
    const auto vroom_mixed_add = g1.add_mixed_point(vroom_generator_projective, vroom_doubled_affine, ring);
    const auto vroom_double = g1.double_point(vroom_generator_projective, ring);

    const auto add_affine = blst_vroom_test::vroom_g1_projective_to_affine_bigints(ring, vroom_add);
    const auto mixed_add_affine = blst_vroom_test::vroom_g1_projective_to_affine_bigints(ring, vroom_mixed_add);
    const auto double_affine = blst_vroom_test::vroom_g1_projective_to_affine_bigints(ring, vroom_double);

    EXPECT_EQ(add_affine.first.to_string(), blst_vroom_test::scalar_fp_to_bigint(scalar_add_affine.X).to_string());
    EXPECT_EQ(add_affine.second.to_string(), blst_vroom_test::scalar_fp_to_bigint(scalar_add_affine.Y).to_string());
    EXPECT_EQ(mixed_add_affine.first.to_string(), blst_vroom_test::scalar_fp_to_bigint(scalar_mixed_add_affine.X).to_string());
    EXPECT_EQ(mixed_add_affine.second.to_string(), blst_vroom_test::scalar_fp_to_bigint(scalar_mixed_add_affine.Y).to_string());
    EXPECT_EQ(double_affine.first.to_string(), blst_vroom_test::scalar_fp_to_bigint(scalar_double_affine.X).to_string());
    EXPECT_EQ(double_affine.second.to_string(), blst_vroom_test::scalar_fp_to_bigint(scalar_double_affine.Y).to_string());
}

TEST(VroomCurve, G2AddMixedAddAndDoubleMatchScalarBlst)
{
    const VroomCurveRing ring(blst_vroom_test::modulus_bigint());
    const VroomCurveFp2Ring fp2ring(ring);
    const VroomCurveG2 g2 = {};

    const blst::G2AffinePoint scalar_generator_affine = blst_vroom_test::make_scalar_g2_generator_affine();
    const blst::G2AffinePoint scalar_doubled_affine = blst_vroom_test::make_scalar_g2_doubled_affine();

    const auto vroom_generator_affine = blst_vroom_test::make_vroom_g2_affine(ring, scalar_generator_affine);
    const auto vroom_doubled_affine = blst_vroom_test::make_vroom_g2_affine(ring, scalar_doubled_affine);
    const auto vroom_generator_projective = blst_vroom_test::make_vroom_curve_g2_projective(fp2ring, vroom_generator_affine);
    const auto vroom_doubled_projective = blst_vroom_test::make_vroom_curve_g2_projective(fp2ring, vroom_doubled_affine);

    blst::G2Point scalar_generator_projective = {};
    blst::G2Point::from_affine(scalar_generator_projective, scalar_generator_affine);
    blst::G2Point scalar_doubled_projective = {};
    blst::G2Point::from_affine(scalar_doubled_projective, scalar_doubled_affine);

    blst::G2Point scalar_add_projective = {};
    blst::G2Point::add(scalar_add_projective, scalar_generator_projective, scalar_doubled_projective);
    blst::G2AffinePoint scalar_add_affine = {};
    blst::G2Point::to_affine(scalar_add_affine, scalar_add_projective);

    blst::G2Point scalar_mixed_add_projective = {};
    blst::G2Point::add(scalar_mixed_add_projective, scalar_generator_projective, scalar_doubled_affine);
    blst::G2AffinePoint scalar_mixed_add_affine = {};
    blst::G2Point::to_affine(scalar_mixed_add_affine, scalar_mixed_add_projective);

    blst::G2Point scalar_double_projective = {};
    blst::G2Point::doubling(scalar_double_projective, scalar_generator_projective);
    blst::G2AffinePoint scalar_double_affine = {};
    blst::G2Point::to_affine(scalar_double_affine, scalar_double_projective);

    const auto vroom_add = g2.add_point(vroom_generator_projective, vroom_doubled_projective, fp2ring);
    const auto vroom_mixed_add = g2.add_mixed_point(vroom_generator_projective, vroom_doubled_affine, fp2ring);
    const auto vroom_double = g2.double_point(vroom_generator_projective, fp2ring);

    const auto add_affine = blst_vroom_test::vroom_g2_projective_to_affine_bigints(ring, vroom_add);
    const auto mixed_add_affine = blst_vroom_test::vroom_g2_projective_to_affine_bigints(ring, vroom_mixed_add);
    const auto double_affine = blst_vroom_test::vroom_g2_projective_to_affine_bigints(ring, vroom_double);

    const auto scalar_add_x = blst_vroom_test::scalar_fp2_to_bigints(scalar_add_affine.X);
    const auto scalar_add_y = blst_vroom_test::scalar_fp2_to_bigints(scalar_add_affine.Y);
    const auto scalar_mixed_add_x = blst_vroom_test::scalar_fp2_to_bigints(scalar_mixed_add_affine.X);
    const auto scalar_mixed_add_y = blst_vroom_test::scalar_fp2_to_bigints(scalar_mixed_add_affine.Y);
    const auto scalar_double_x = blst_vroom_test::scalar_fp2_to_bigints(scalar_double_affine.X);
    const auto scalar_double_y = blst_vroom_test::scalar_fp2_to_bigints(scalar_double_affine.Y);

    EXPECT_EQ(add_affine.first.first.to_string(), scalar_add_x.first.to_string());
    EXPECT_EQ(add_affine.first.second.to_string(), scalar_add_x.second.to_string());
    EXPECT_EQ(add_affine.second.first.to_string(), scalar_add_y.first.to_string());
    EXPECT_EQ(add_affine.second.second.to_string(), scalar_add_y.second.to_string());
    EXPECT_EQ(mixed_add_affine.first.first.to_string(), scalar_mixed_add_x.first.to_string());
    EXPECT_EQ(mixed_add_affine.first.second.to_string(), scalar_mixed_add_x.second.to_string());
    EXPECT_EQ(mixed_add_affine.second.first.to_string(), scalar_mixed_add_y.first.to_string());
    EXPECT_EQ(mixed_add_affine.second.second.to_string(), scalar_mixed_add_y.second.to_string());
    EXPECT_EQ(double_affine.first.first.to_string(), scalar_double_x.first.to_string());
    EXPECT_EQ(double_affine.first.second.to_string(), scalar_double_x.second.to_string());
    EXPECT_EQ(double_affine.second.first.to_string(), scalar_double_y.first.to_string());
    EXPECT_EQ(double_affine.second.second.to_string(), scalar_double_y.second.to_string());
}

}  // namespace
