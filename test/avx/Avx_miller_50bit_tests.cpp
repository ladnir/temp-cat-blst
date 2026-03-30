#include <gtest/gtest.h>

#include "blst/avx/pairing/miller.h"
#include "test/avx/avx_test_utils.h"

#include <array>
#include <string>

using namespace blst::avx;

namespace {

using blst_avx_test::Avx50Fp12Ring;
using blst_avx_test::Avx50Fp2Ring;
using blst_avx_test::Avx50Ring;

void expect_fp12_bigints_equal_50bit_miller(
    const std::array<BigInt, 12>& actual,
    const std::array<BigInt, 12>& expected,
    const std::string& label)
{
    for (size_t i = 0; i < actual.size(); ++i) {
        SCOPED_TRACE(label + " coeff=" + std::to_string(i));
        EXPECT_EQ(actual[i].to_string(), expected[i].to_string());
    }
}

TEST(AvxPairing50Bit, MillerLoopThenScalarFinalExpMatchesScalarBlstOnReferencePairs)
{
    const Avx50Ring ring(blst_avx_test::modulus_bigint());
    const Avx50Fp2Ring fp2ring(ring);
    const Avx50Fp12Ring fp12ring(ring);
    const std::array<std::pair<blst::G1AffinePoint, blst::G2AffinePoint>, 3> cases = {{
        {blst_avx_test::make_scalar_g1_generator_affine(), blst_avx_test::make_scalar_g2_generator_affine()},
        {blst_avx_test::make_scalar_g1_doubled_affine(), blst_avx_test::make_scalar_g2_generator_affine()},
        {blst_avx_test::make_scalar_g1_generator_affine(), blst_avx_test::make_scalar_g2_doubled_affine()},
    }};

    for (size_t case_idx = 0; case_idx < cases.size(); ++case_idx) {
        const blst::G1AffinePoint& scalar_g1 = cases[case_idx].first;
        const blst::G2AffinePoint& scalar_g2 = cases[case_idx].second;
        const auto avx_g1 = blst_avx_test::make_avx_g1_affine(ring, scalar_g1);
        const auto avx_g2 = blst_avx_test::make_avx_g2_affine(ring, scalar_g2);
        Miller<Avx50Ring> miller(avx_g1, avx_g2, ring);

        blst::Fp12 actual_final = blst_avx_test::bigint_array_to_scalar_fp12(
            blst_avx_test::avx_fp12_to_bigints(ring, miller.miller_loop(ring, fp2ring, fp12ring)));
        blst::Fp12::final_exp(actual_final, actual_final);
        const auto actual = blst_avx_test::scalar_fp12_to_bigints(actual_final);
        const auto expected = blst_avx_test::scalar_fp12_to_bigints(
            blst_avx_test::make_scalar_final_pairing(scalar_g1, scalar_g2).data);

        expect_fp12_bigints_equal_50bit_miller(actual, expected, "case=" + std::to_string(case_idx));
    }
}

}  // namespace
