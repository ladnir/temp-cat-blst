#include <gtest/gtest.h>

#include "blst/avx/src/conversion_inversion.hpp"
#include "blst/avx/src/final_exponentiation.hpp"
#include "blst/avx/tests/vroom_basefield_test_utils.hpp"

#include <array>
#include <string>

namespace {

using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomFp12;
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

TEST(VroomFp12, CyclotomicSqrMatchesScalarBlstFp12)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const std::array<uint64_t, 4> seeds = {{42, 123, 456, 789}};

    for (size_t case_idx = 0; case_idx < seeds.size(); ++case_idx) {
        const uint64_t seed = seeds[case_idx];
        const auto vroom_value = blst_vroom_test::make_vroom_cyclotomic_fp12_from_seed(ring, seed);
        const auto vroom_result = fp12ring.fp12_flat_cyclotomic_sqr(vroom_value, ring);

        const blst::Fp12 scalar_value = blst_vroom_test::make_scalar_cyclotomic_fp12_from_seed(seed);
        ASSERT_TRUE(blst::Fp12::is_cyclotomic(scalar_value));

        blst::Fp12 scalar_result = {};
        blst::Fp12::cyclotomic_sqr(scalar_result, scalar_value);
        ASSERT_TRUE(blst::Fp12::is_cyclotomic(scalar_result));

        const auto vroom_bigints = blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_result);
        const auto scalar_bigints = blst_vroom_test::scalar_fp12_to_bigints(scalar_result);

        for (size_t i = 0; i < vroom_bigints.size(); ++i) {
            SCOPED_TRACE("case=" + std::to_string(case_idx) + " coeff=" + std::to_string(i));
            EXPECT_EQ(vroom_bigints[i].to_string(), scalar_bigints[i].to_string());
        }
    }
}

TEST(VroomFp12, MulByXy00z0MatchesScalarBlstFp12)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
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
        const auto c_term = blst_vroom_test::make_vroom_fp2_from_sparse6<VroomRing>(sparse_coeffs, 0);
        const auto x_term = blst_vroom_test::make_vroom_fp2_from_sparse6<VroomRing>(sparse_coeffs, 1);
        const auto y_term = blst_vroom_test::make_vroom_fp2_from_sparse6<VroomRing>(sparse_coeffs, 2);

        VroomFp12 result = {};
        fp12ring.flat_mulxy00z0(result, vroom_value, c_term, x_term, y_term, ring);

        const blst::Fp12 scalar_value = blst_vroom_test::make_scalar_fp12_from_seed(x_seed);
        const blst::Fp6 scalar_sparse = blst_vroom_test::make_scalar_sparse_fp6_from_seed(sparse_seed);
        blst::Fp12 scalar_result = {};
        blst::Fp12::mul_by_xy00z0(scalar_result, scalar_value, scalar_sparse);

        const auto vroom_bigints = blst_vroom_test::vroom_fp12_to_bigints(ring, result);
        const auto scalar_bigints = blst_vroom_test::scalar_fp12_to_bigints(scalar_result);

        for (size_t i = 0; i < vroom_bigints.size(); ++i) {
            SCOPED_TRACE("case=" + std::to_string(case_idx) + " coeff=" + std::to_string(i));
            EXPECT_EQ(vroom_bigints[i].to_string(), scalar_bigints[i].to_string());
        }
    }
}

TEST(VroomFp12, MulByXy00z0FixedVectorsMatchScalarBlstFp12)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const auto fp12_coeffs = make_boundary_coeffs(modulus);
    const std::array<BigInt, 6> sparse_coeffs = {
        BigInt(1), BigInt(2), BigInt(3), BigInt(5), BigInt(8), BigInt(13),
    };

    const auto vroom_value = blst_vroom_test::make_vroom_fp12(ring, fp12_coeffs);
    const auto sparse_vroom = blst_vroom_test::make_vroom_sparse_fp6_coeffs(ring, sparse_coeffs);
    const auto c_term = blst_vroom_test::make_vroom_fp2_from_sparse6<VroomRing>(sparse_vroom, 0);
    const auto x_term = blst_vroom_test::make_vroom_fp2_from_sparse6<VroomRing>(sparse_vroom, 1);
    const auto y_term = blst_vroom_test::make_vroom_fp2_from_sparse6<VroomRing>(sparse_vroom, 2);

    VroomFp12 result = {};
    fp12ring.flat_mulxy00z0(result, vroom_value, c_term, x_term, y_term, ring);

    const blst::Fp12 scalar_value = blst_vroom_test::bigint_array_to_scalar_fp12(fp12_coeffs);
    const blst::Fp6 scalar_sparse = blst_vroom_test::bigint_array_to_scalar_fp6(sparse_coeffs);
    blst::Fp12 scalar_result = {};
    blst::Fp12::mul_by_xy00z0(scalar_result, scalar_value, scalar_sparse);

    expect_fp12_bigints_equal(
        blst_vroom_test::vroom_fp12_to_bigints(ring, result),
        blst_vroom_test::scalar_fp12_to_bigints(scalar_result),
        "fixed sparse mul");
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

TEST(VroomFp12, BlstInverterMatchesScalarBlstInverse)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const std::array<uint64_t, 4> seeds = {{42, 123, 456, 789}};

    for (size_t case_idx = 0; case_idx < seeds.size(); ++case_idx) {
        const uint64_t seed = seeds[case_idx];
        const auto coeffs = blst_vroom_test::seeded_fp12_coeffs(seed);
        const auto vroom_value = blst_vroom_test::make_vroom_fp12(ring, coeffs);
        const auto vroom_inverse = invert_fp12_via_blst(ring, vroom_value);

        const blst::Fp12 scalar_value = blst_vroom_test::bigint_array_to_scalar_fp12(coeffs);
        blst::Fp12 scalar_inverse = {};
        blst::Fp12::inverse(scalar_inverse, scalar_value);

        expect_fp12_bigints_equal(
            blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_inverse),
            blst_vroom_test::scalar_fp12_to_bigints(scalar_inverse),
            "inverse case=" + std::to_string(case_idx));
    }
}

TEST(VroomFp12, ConjugateMatchesScalarBlstFp12)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const std::array<uint64_t, 4> seeds = {{42, 123, 456, 789}};

    for (size_t case_idx = 0; case_idx < seeds.size(); ++case_idx) {
        const uint64_t seed = seeds[case_idx];
        const auto vroom_value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);
        const auto vroom_result = fp12ring.fp12_flat_conjugate(vroom_value, ring);

        blst::Fp12 scalar_result = blst_vroom_test::make_scalar_fp12_from_seed(seed);
        blst::Fp12::conjugate(scalar_result);

        expect_fp12_bigints_equal(
            blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_result),
            blst_vroom_test::scalar_fp12_to_bigints(scalar_result),
            "conjugate case=" + std::to_string(case_idx));
    }
}

TEST(VroomFp12, FrobeniusMapsMatchScalarBlstFp12)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const std::array<uint64_t, 4> seeds = {{42, 123, 456, 789}};

    for (size_t case_idx = 0; case_idx < seeds.size(); ++case_idx) {
        const uint64_t seed = seeds[case_idx];
        const auto vroom_value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);
        const blst::Fp12 scalar_value = blst_vroom_test::make_scalar_fp12_from_seed(seed);

        const auto vroom_frob1 = fp12ring.fp12_flat_frobenius_map1(vroom_value, ring);
        const auto vroom_frob2 = fp12ring.fp12_flat_frobenius_map2(vroom_value, ring);
        const auto vroom_frob3 = fp12ring.fp12_flat_frobenius_map3(vroom_value, ring);

        blst::Fp12 scalar_frob1 = {};
        blst::Fp12 scalar_frob2 = {};
        blst::Fp12 scalar_frob3 = {};
        blst::Fp12::frobenius_map(scalar_frob1, scalar_value, 1);
        blst::Fp12::frobenius_map(scalar_frob2, scalar_value, 2);
        blst::Fp12::frobenius_map(scalar_frob3, scalar_value, 3);

        expect_fp12_bigints_equal(
            blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_frob1),
            blst_vroom_test::scalar_fp12_to_bigints(scalar_frob1),
            "frob1 case=" + std::to_string(case_idx));
        expect_fp12_bigints_equal(
            blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_frob2),
            blst_vroom_test::scalar_fp12_to_bigints(scalar_frob2),
            "frob2 case=" + std::to_string(case_idx));
        expect_fp12_bigints_equal(
            blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_frob3),
            blst_vroom_test::scalar_fp12_to_bigints(scalar_frob3),
            "frob3 case=" + std::to_string(case_idx));
    }
}

TEST(VroomFp12, FinalExpAfterInverseMatchesScalarBlstFinalExp)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const std::array<std::pair<uint64_t, uint64_t>, 4> seeds = {{
        {42, 123},
        {456, 789},
        {17, 33},
        {99, 7},
    }};

    for (size_t case_idx = 0; case_idx < seeds.size(); ++case_idx) {
        const uint64_t value_seed = seeds[case_idx].first;
        const auto coeffs = blst_vroom_test::seeded_fp12_coeffs(value_seed);
        const auto vroom_value = blst_vroom_test::make_vroom_fp12(ring, coeffs);

        const blst::Fp12 scalar_value = blst_vroom_test::bigint_array_to_scalar_fp12(coeffs);
        blst::Fp12 scalar_inverse = {};
        blst::Fp12 scalar_result = {};
        blst::Fp12::inverse(scalar_inverse, scalar_value);
        blst::Fp12::final_exp(scalar_result, scalar_value);

        const auto vroom_inverse =
            blst_vroom_test::make_vroom_fp12(ring, blst_vroom_test::scalar_fp12_to_bigints(scalar_inverse));
        const auto vroom_result = final_exp_after_inverse(vroom_value, vroom_inverse, fp12ring, ring);

        expect_fp12_bigints_equal(
            blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_result),
            blst_vroom_test::scalar_fp12_to_bigints(scalar_result),
            "final-exp-after-inverse case=" + std::to_string(case_idx));
    }
}

TEST(VroomFp12, RnsBlstInverterMatchesScalarBlstInverse)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const FP12RNS_BLSTInverter<VroomRing> inverter = {};
    const std::array<uint64_t, 4> seeds = {{42, 123, 456, 789}};

    for (size_t case_idx = 0; case_idx < seeds.size(); ++case_idx) {
        const uint64_t seed = seeds[case_idx];
        const auto coeffs = blst_vroom_test::seeded_fp12_coeffs(seed);
        const auto vroom_value = blst_vroom_test::make_vroom_fp12(ring, coeffs);
        const auto vroom_inverse = inverter.invert(vroom_value, fp12ring, ring);

        const blst::Fp12 scalar_value = blst_vroom_test::bigint_array_to_scalar_fp12(coeffs);
        blst::Fp12 scalar_inverse = {};
        blst::Fp12::inverse(scalar_inverse, scalar_value);

        expect_fp12_bigints_equal(
            blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_inverse),
            blst_vroom_test::scalar_fp12_to_bigints(scalar_inverse),
            "rns inverse case=" + std::to_string(case_idx));
    }
}

template<class InverterFactory>
void expect_final_exp_matches_scalar(const std::string& label, InverterFactory make_inverter)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const std::array<uint64_t, 4> seeds = {{42, 123, 456, 789}};

    for (size_t case_idx = 0; case_idx < seeds.size(); ++case_idx) {
        const uint64_t seed = seeds[case_idx];
        const auto coeffs = blst_vroom_test::seeded_fp12_coeffs(seed);
        const auto vroom_value = blst_vroom_test::make_vroom_fp12(ring, coeffs);
        const auto inverter = make_inverter();
        const auto vroom_result = final_exp(vroom_value, inverter, fp12ring, ring);

        const blst::Fp12 scalar_value = blst_vroom_test::bigint_array_to_scalar_fp12(coeffs);
        blst::Fp12 scalar_result = {};
        blst::Fp12::final_exp(scalar_result, scalar_value);

        expect_fp12_bigints_equal(
            blst_vroom_test::vroom_fp12_to_bigints(ring, vroom_result),
            blst_vroom_test::scalar_fp12_to_bigints(scalar_result),
            label + " case=" + std::to_string(case_idx));
    }
}

TEST(VroomFp12, FinalExpBlstInverterMatchesScalarBlst)
{
    expect_final_exp_matches_scalar(
        "final exp blst inverter",
        []() { return FP12BLSTInverter<VroomRing>{}; });
}

TEST(VroomFp12, FinalExpRnsBlstInverterMatchesScalarBlst)
{
    expect_final_exp_matches_scalar(
        "final exp rns blst inverter",
        []() { return FP12RNS_BLSTInverter<VroomRing>{}; });
}

}  // namespace
