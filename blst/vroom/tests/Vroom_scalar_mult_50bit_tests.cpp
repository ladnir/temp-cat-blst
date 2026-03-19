#include <gtest/gtest.h>

#define BLST_VROOM_SKIP_TEST_VALUES
#include "blst/vroom/src/scalar_mult.hpp"
#include "blst/vroom/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::VroomCurve50Fp2Ring;
using blst_vroom_test::VroomCurve50G1;
using blst_vroom_test::VroomCurve50G2;
using blst_vroom_test::VroomCurve50Ring;

TEST(VroomScalarMult50Bit, G1ScalarMultGlvMatchesScalarBlst)
{
    const VroomCurve50Ring ring(blst_vroom_test::modulus_bigint());
    const VroomCurve50G1 g1 = {};
    const GLVMult<VroomCurve50Ring> glv(g1, ring);

    for (const auto& test_case : std::array<std::pair<uint64_t, uint64_t>, 2>{{{42, 123}, {777, 321}}}) {
        const blst::G1Point scalar_base = blst_vroom_test::make_scalar_g1_projective_from_seed(test_case.first);
        blst::G1AffinePoint scalar_base_affine = {};
        blst::G1Point::to_affine(scalar_base_affine, scalar_base);

        const auto vroom_base_affine = blst_vroom_test::make_vroom_g1_affine(ring, scalar_base_affine);
        const auto vroom_base = blst_vroom_test::make_vroom_curve50_g1_projective(ring, vroom_base_affine);
        const BigInt scalar_bigint = blst_vroom_test::seeded_scalar_bigint(test_case.second);
        const auto scalar_bytes = blst_vroom_test::seeded_scalar_bytes_le(test_case.second);

        blst::G1Point scalar_result = {};
        blst::G1Point::mult(scalar_result, scalar_base, scalar_bytes.data(), 256);
        blst::G1AffinePoint scalar_result_affine = {};
        blst::G1Point::to_affine(scalar_result_affine, scalar_result);

        const auto vroom_result_bytes = glv.multiply(vroom_base, scalar_bytes, ring);
        const auto vroom_result_bigint = glv.multiply(vroom_base, scalar_bigint, ring);
        const auto vroom_affine_bytes = blst_vroom_test::vroom_g1_projective_to_affine_bigints(ring, vroom_result_bytes);
        const auto vroom_affine_bigint = blst_vroom_test::vroom_g1_projective_to_affine_bigints(ring, vroom_result_bigint);

        const auto expected_x = blst_vroom_test::scalar_fp_to_bigint(scalar_result_affine.X).to_string();
        const auto expected_y = blst_vroom_test::scalar_fp_to_bigint(scalar_result_affine.Y).to_string();

        EXPECT_EQ(vroom_affine_bytes.first.to_string(), expected_x);
        EXPECT_EQ(vroom_affine_bytes.second.to_string(), expected_y);
        EXPECT_EQ(vroom_affine_bigint.first.to_string(), expected_x);
        EXPECT_EQ(vroom_affine_bigint.second.to_string(), expected_y);
    }
}

TEST(VroomScalarMult50Bit, G2ScalarMultGlsMatchesScalarBlst)
{
    const VroomCurve50Ring ring(blst_vroom_test::modulus_bigint());
    const VroomCurve50Fp2Ring fp2ring(ring);
    const VroomCurve50G2 g2 = {};
    const GLSMult<VroomCurve50Ring> gls(ring);

    for (const auto& test_case : std::array<std::pair<uint64_t, uint64_t>, 2>{{{42, 321}, {777, 123}}}) {
        const blst::G2Point scalar_base = blst_vroom_test::make_scalar_g2_projective_from_seed(test_case.first);
        blst::G2AffinePoint scalar_base_affine = {};
        blst::G2Point::to_affine(scalar_base_affine, scalar_base);

        const auto vroom_base_affine = blst_vroom_test::make_vroom_g2_affine(ring, scalar_base_affine);
        const auto vroom_base = blst_vroom_test::make_vroom_curve50_g2_projective(fp2ring, vroom_base_affine);
        const BigInt scalar_bigint = blst_vroom_test::seeded_scalar_bigint(test_case.second);
        const auto scalar_bytes = blst_vroom_test::seeded_scalar_bytes_le(test_case.second);

        blst::G2Point scalar_result = {};
        blst::G2Point::mult(scalar_result, scalar_base, scalar_bytes.data(), 256);
        blst::G2AffinePoint scalar_result_affine = {};
        blst::G2Point::to_affine(scalar_result_affine, scalar_result);

        const auto vroom_result_bytes = gls.multiply(vroom_base, scalar_bytes, fp2ring, g2, ring);
        const auto vroom_result_bigint = gls.multiply(vroom_base, scalar_bigint, fp2ring, g2, ring);
        const auto vroom_affine_bytes = blst_vroom_test::vroom_g2_projective_to_affine_bigints(ring, vroom_result_bytes);
        const auto vroom_affine_bigint = blst_vroom_test::vroom_g2_projective_to_affine_bigints(ring, vroom_result_bigint);

        const auto expected_x = blst_vroom_test::scalar_fp2_to_bigints(scalar_result_affine.X);
        const auto expected_y = blst_vroom_test::scalar_fp2_to_bigints(scalar_result_affine.Y);

        EXPECT_EQ(vroom_affine_bytes.first.first.to_string(), expected_x.first.to_string());
        EXPECT_EQ(vroom_affine_bytes.first.second.to_string(), expected_x.second.to_string());
        EXPECT_EQ(vroom_affine_bytes.second.first.to_string(), expected_y.first.to_string());
        EXPECT_EQ(vroom_affine_bytes.second.second.to_string(), expected_y.second.to_string());
        EXPECT_EQ(vroom_affine_bigint.first.first.to_string(), expected_x.first.to_string());
        EXPECT_EQ(vroom_affine_bigint.first.second.to_string(), expected_x.second.to_string());
        EXPECT_EQ(vroom_affine_bigint.second.first.to_string(), expected_y.first.to_string());
        EXPECT_EQ(vroom_affine_bigint.second.second.to_string(), expected_y.second.to_string());
    }
}

}  // namespace
