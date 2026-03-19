#include <gtest/gtest.h>

#include "blst/avx/src/pairing.hpp"
#include "blst/avx/tests/vroom_basefield_test_utils.hpp"

#include <array>
#include <string>

namespace {

using blst_vroom_test::VroomAffineFp;
using blst_vroom_test::VroomAffineFp2;
using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomFp2Ring;
using blst_vroom_test::VroomRing;

void expect_fp12_bigints_equal(
    const std::array<BigInt, 12>& actual,
    const std::array<BigInt, 12>& expected,
    const std::string& label)
{
    for (size_t i = 0; i < actual.size(); ++i) {
        SCOPED_TRACE(label + " coeff=" + std::to_string(i));
        EXPECT_EQ(actual[i].to_string(), expected[i].to_string());
    }
}

template<class InverterFactory>
void expect_pairing_matches_scalar(const std::string& label, InverterFactory make_inverter)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp2Ring fp2ring(ring);
    const VroomFp12Ring fp12ring(ring);
    const std::array<std::pair<blst::G1AffinePoint, blst::G2AffinePoint>, 3> cases = {{
        {blst_vroom_test::make_scalar_g1_generator_affine(), blst_vroom_test::make_scalar_g2_generator_affine()},
        {blst_vroom_test::make_scalar_g1_doubled_affine(), blst_vroom_test::make_scalar_g2_generator_affine()},
        {blst_vroom_test::make_scalar_g1_generator_affine(), blst_vroom_test::make_scalar_g2_doubled_affine()},
    }};

    for (size_t case_idx = 0; case_idx < cases.size(); ++case_idx) {
        const blst::G1AffinePoint& scalar_g1 = cases[case_idx].first;
        const blst::G2AffinePoint& scalar_g2 = cases[case_idx].second;
        const VroomAffineFp vroom_g1 = blst_vroom_test::make_vroom_g1_affine(ring, scalar_g1);
        const VroomAffineFp2 vroom_g2 = blst_vroom_test::make_vroom_g2_affine(ring, scalar_g2);
        const auto inverter = make_inverter();

        const auto actual = blst_vroom_test::vroom_fp12_to_bigints(
            ring,
            pairing(vroom_g1, vroom_g2, inverter, fp2ring, fp12ring, ring));
        const auto expected = blst_vroom_test::scalar_fp12_to_bigints(
            blst_vroom_test::make_scalar_final_pairing(scalar_g1, scalar_g2).data);

        expect_fp12_bigints_equal(actual, expected, label + " case=" + std::to_string(case_idx));
    }
}

TEST(VroomPairing, PairingBlstInverterMatchesScalarBlst)
{
    expect_pairing_matches_scalar(
        "pairing blst inverter",
        []() { return FP12BLSTInverter<VroomRing>{}; });
}

TEST(VroomPairing, PairingRnsBlstInverterMatchesScalarBlst)
{
    expect_pairing_matches_scalar(
        "pairing rns blst inverter",
        []() { return FP12RNS_BLSTInverter<VroomRing>{}; });
}

}  // namespace
