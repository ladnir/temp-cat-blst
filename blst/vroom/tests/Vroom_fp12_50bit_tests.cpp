#include <gtest/gtest.h>

#include "blst/vroom/src/conversion_inversion.hpp"
#include "blst/vroom/src/final_exponentiation.hpp"
#include "blst/vroom/tests/vroom_basefield_test_utils.hpp"

#include <array>
#include <string>

namespace {

using blst_vroom_test::Vroom50Fp12;
using blst_vroom_test::Vroom50Fp12Ring;
using blst_vroom_test::Vroom50Ring;

void expect_fp12_bigints_equal_50bit(
    const std::array<BigInt, 12>& actual,
    const std::array<BigInt, 12>& expected,
    const std::string& label)
{
    for (size_t i = 0; i < actual.size(); ++i) {
        SCOPED_TRACE(label + " coeff=" + std::to_string(i));
        EXPECT_EQ(actual[i].to_string(), expected[i].to_string());
    }
}

TEST(VroomFp1250Bit, MulMatchesScalarBlstFp12)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(456);

    for (int case_idx = 0; case_idx < 4; ++case_idx) {
        std::array<BigInt, 12> left_coeffs;
        std::array<BigInt, 12> right_coeffs;
        for (size_t i = 0; i < left_coeffs.size(); ++i) {
            left_coeffs[i] = BigInt(rng.get_z_range(modulus.get_mpz()));
            right_coeffs[i] = BigInt(rng.get_z_range(modulus.get_mpz()));
        }

        const auto vroom_left = blst_vroom_test::make_vroom_fp12(ring, left_coeffs);
        const auto vroom_right = blst_vroom_test::make_vroom_fp12(ring, right_coeffs);
        const auto vroom_result = fp12ring.fp12_flat_mul(vroom_left, vroom_right, ring);

        const blst::Fp12 scalar_left = blst_vroom_test::bigint_array_to_scalar_fp12(left_coeffs);
        const blst::Fp12 scalar_right = blst_vroom_test::bigint_array_to_scalar_fp12(right_coeffs);
        blst::Fp12 scalar_result = {};
        blst::Fp12::mul(scalar_result, scalar_left, scalar_right);

        expect_fp12_bigints_equal_50bit(
            blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_result),
            blst_vroom_test::scalar_fp12_to_bigints(scalar_result),
            "case=" + std::to_string(case_idx));
    }
}

TEST(VroomFp1250Bit, CyclotomicSqrMatchesScalarBlstFp12)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    const std::array<uint64_t, 4> seeds = {{42, 123, 456, 789}};

    for (size_t case_idx = 0; case_idx < seeds.size(); ++case_idx) {
        const uint64_t seed = seeds[case_idx];
        const auto vroom_value = blst_vroom_test::make_vroom_cyclotomic_fp12_from_seed(ring, seed);
        const auto vroom_result = fp12ring.fp12_flat_cyclotomic_sqr(vroom_value, ring);

        const blst::Fp12 scalar_value = blst_vroom_test::make_scalar_cyclotomic_fp12_from_seed(seed);
        blst::Fp12 scalar_result = {};
        blst::Fp12::cyclotomic_sqr(scalar_result, scalar_value);

        expect_fp12_bigints_equal_50bit(
            blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_result),
            blst_vroom_test::scalar_fp12_to_bigints(scalar_result),
            "case=" + std::to_string(case_idx));
    }
}

TEST(VroomFp1250Bit, MulByXy00z0MatchesScalarBlstFp12)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    const std::array<std::pair<uint64_t, uint64_t>, 4> seeds = {{
        {42, 123},
        {456, 789},
        {17, 33},
        {99, 7},
    }};

    for (size_t case_idx = 0; case_idx < seeds.size(); ++case_idx) {
        const auto [x_seed, sparse_seed] = seeds[case_idx];
        const auto vroom_value = blst_vroom_test::make_vroom_fp12_from_seed(ring, x_seed);
        const auto sparse_coeffs = blst_vroom_test::make_vroom_sparse_fp6_from_seed(ring, sparse_seed);
        const auto c_term = blst_vroom_test::make_vroom_fp2_from_sparse6<Vroom50Ring>(sparse_coeffs, 0);
        const auto x_term = blst_vroom_test::make_vroom_fp2_from_sparse6<Vroom50Ring>(sparse_coeffs, 1);
        const auto y_term = blst_vroom_test::make_vroom_fp2_from_sparse6<Vroom50Ring>(sparse_coeffs, 2);

        Vroom50Fp12 result = {};
        fp12ring.flat_mulxy00z0(result, vroom_value, c_term, x_term, y_term, ring);

        const blst::Fp12 scalar_value = blst_vroom_test::make_scalar_fp12_from_seed(x_seed);
        const blst::Fp6 scalar_sparse = blst_vroom_test::make_scalar_sparse_fp6_from_seed(sparse_seed);
        blst::Fp12 scalar_result = {};
        blst::Fp12::mul_by_xy00z0(scalar_result, scalar_value, scalar_sparse);

        expect_fp12_bigints_equal_50bit(
            blst_vroom_test::vroom_fp12_to_bigints(ring, result),
            blst_vroom_test::scalar_fp12_to_bigints(scalar_result),
            "case=" + std::to_string(case_idx));
    }
}

TEST(VroomFp1250Bit, BlstInverterMatchesScalarBlstInverse)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, 42);
    const auto actual = invert_fp12_via_blst(ring, value);

    blst::Fp12 scalar_value = blst_vroom_test::make_scalar_fp12_from_seed(42);
    blst::Fp12 scalar_inverse = {};
    blst::Fp12::inverse(scalar_inverse, scalar_value);

    expect_fp12_bigints_equal_50bit(
        blst_vroom_test::vroom_fp12_to_bigints(ring, actual),
        blst_vroom_test::scalar_fp12_to_bigints(scalar_inverse),
        "blst inverter");
}

TEST(VroomFp1250Bit, RnsBlstInverterMatchesScalarBlstInverse)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    const FP12RNS_BLSTInverter<Vroom50Ring> inverter = {};
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, 42);
    const auto actual = inverter.invert(value, fp12ring, ring);

    blst::Fp12 scalar_value = blst_vroom_test::make_scalar_fp12_from_seed(42);
    blst::Fp12 scalar_inverse = {};
    blst::Fp12::inverse(scalar_inverse, scalar_value);

    expect_fp12_bigints_equal_50bit(
        blst_vroom_test::vroom_fp12_to_bigints(ring, actual),
        blst_vroom_test::scalar_fp12_to_bigints(scalar_inverse),
        "rns inverter");
}

TEST(VroomFp1250Bit, FinalExpBlstInverterMatchesScalarBlst)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    const FP12BLSTInverter<Vroom50Ring> inverter = {};
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, 42);
    const auto actual = final_exp(value, inverter, fp12ring, ring);

    blst::Fp12 scalar_value = blst_vroom_test::make_scalar_fp12_from_seed(42);
    blst::Fp12 scalar_final = {};
    blst::Fp12::final_exp(scalar_final, scalar_value);

    expect_fp12_bigints_equal_50bit(
        blst_vroom_test::vroom_fp12_to_bigints(ring, actual),
        blst_vroom_test::scalar_fp12_to_bigints(scalar_final),
        "final exp blst");
}

TEST(VroomFp1250Bit, FinalExpRnsBlstInverterMatchesScalarBlst)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    const FP12RNS_BLSTInverter<Vroom50Ring> inverter = {};
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, 42);
    const auto actual = final_exp(value, inverter, fp12ring, ring);

    blst::Fp12 scalar_value = blst_vroom_test::make_scalar_fp12_from_seed(42);
    blst::Fp12 scalar_final = {};
    blst::Fp12::final_exp(scalar_final, scalar_value);

    expect_fp12_bigints_equal_50bit(
        blst_vroom_test::vroom_fp12_to_bigints(ring, actual),
        blst_vroom_test::scalar_fp12_to_bigints(scalar_final),
        "final exp rns");
}

}  // namespace
