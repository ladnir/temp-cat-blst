#include <gtest/gtest.h>

#include "blst/vroom/tests/vroom_basefield_test_utils.hpp"

#include <array>
#include <string>

namespace {

using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomRing;

std::array<BigInt, 12> make_linear_coeffs(const BigInt& modulus, uint64_t base, uint64_t stride)
{
    std::array<BigInt, 12> coeffs;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        coeffs[i] = BigInt((base + stride * i) % 97);
        if (coeffs[i] >= modulus) {
            coeffs[i] %= modulus;
        }
    }
    return coeffs;
}

std::array<BigInt, 12> make_boundary_coeffs(const BigInt& modulus)
{
    const BigInt zero(0);
    const BigInt one(1);
    const BigInt minus_one = modulus - one;
    const BigInt max_reduced = modulus - one;

    return {
        zero, one, minus_one, max_reduced,
        BigInt(2), BigInt(3), BigInt(5), BigInt(8),
        BigInt(13), BigInt(21), BigInt(34), BigInt(55),
    };
}

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

std::array<BigInt, 12> scalar_fp12_mul_expected(
    const std::array<BigInt, 12>& left,
    const std::array<BigInt, 12>& right)
{
    const blst::Fp12 scalar_left = blst_vroom_test::bigint_array_to_scalar_fp12(left);
    const blst::Fp12 scalar_right = blst_vroom_test::bigint_array_to_scalar_fp12(right);
    blst::Fp12 scalar_result = {};
    blst::Fp12::mul(scalar_result, scalar_left, scalar_right);
    return blst_vroom_test::scalar_fp12_to_bigints(scalar_result);
}

TEST(VroomFp12, MulMatchesScalarBlstFp12)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
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

        const auto vroom_bigints = blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_result);
        const auto scalar_bigints = blst_vroom_test::scalar_fp12_to_bigints(scalar_result);

        for (size_t i = 0; i < vroom_bigints.size(); ++i) {
            SCOPED_TRACE("case=" + std::to_string(case_idx) + " coeff=" + std::to_string(i));
            EXPECT_EQ(vroom_bigints[i].to_string(), scalar_bigints[i].to_string());
        }
    }
}

TEST(VroomFp12, MulFixedVectorsMatchScalarBlstFp12)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const std::array<std::pair<std::array<BigInt, 12>, std::array<BigInt, 12>>, 3> cases = {{
        {make_linear_coeffs(modulus, 1, 1), make_linear_coeffs(modulus, 17, 3)},
        {make_boundary_coeffs(modulus), make_linear_coeffs(modulus, 7, 5)},
        {make_linear_coeffs(modulus, 0, 0), make_boundary_coeffs(modulus)},
    }};

    for (size_t case_idx = 0; case_idx < cases.size(); ++case_idx) {
        const auto& left_coeffs = cases[case_idx].first;
        const auto& right_coeffs = cases[case_idx].second;
        const auto vroom_left = blst_vroom_test::make_vroom_fp12(ring, left_coeffs);
        const auto vroom_right = blst_vroom_test::make_vroom_fp12(ring, right_coeffs);
        const auto vroom_result = fp12ring.fp12_flat_mul(vroom_left, vroom_right, ring);
        const auto expected = scalar_fp12_mul_expected(left_coeffs, right_coeffs);

        expect_fp12_bigints_equal(
            blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_result),
            expected,
            "case=" + std::to_string(case_idx));
    }
}

TEST(VroomFp12, SqrMatchesScalarBlstFp12)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(789);

    for (int case_idx = 0; case_idx < 4; ++case_idx) {
        std::array<BigInt, 12> coeffs;
        for (BigInt& coeff : coeffs) {
            coeff = BigInt(rng.get_z_range(modulus.get_mpz()));
        }

        const auto vroom_value = blst_vroom_test::make_vroom_fp12(ring, coeffs);
        const auto vroom_result = fp12ring.fp12_flat_sqr(vroom_value, ring);

        const blst::Fp12 scalar_value = blst_vroom_test::bigint_array_to_scalar_fp12(coeffs);
        blst::Fp12 scalar_result = {};
        blst::Fp12::sqr(scalar_result, scalar_value);

        const auto vroom_bigints = blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_result);
        const auto scalar_bigints = blst_vroom_test::scalar_fp12_to_bigints(scalar_result);

        for (size_t i = 0; i < vroom_bigints.size(); ++i) {
            SCOPED_TRACE("case=" + std::to_string(case_idx) + " coeff=" + std::to_string(i));
            EXPECT_EQ(vroom_bigints[i].to_string(), scalar_bigints[i].to_string());
        }
    }
}

TEST(VroomFp12, RoundTripBoundaryCoefficientsPreserveCanonicalValues)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const std::array<std::array<BigInt, 12>, 3> cases = {{
        make_linear_coeffs(modulus, 1, 1),
        make_boundary_coeffs(modulus),
        make_linear_coeffs(modulus, 23, 11),
    }};

    for (size_t case_idx = 0; case_idx < cases.size(); ++case_idx) {
        const auto& coeffs = cases[case_idx];
        const auto vroom_value = blst_vroom_test::make_vroom_fp12(ring, coeffs);
        const auto round_trip = blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_value);

        expect_fp12_bigints_equal(round_trip, coeffs, "round-trip case=" + std::to_string(case_idx));
    }
}

}  // namespace
