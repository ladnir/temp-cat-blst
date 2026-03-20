#include <gtest/gtest.h>

#include "blst/avx/src/inversion.hpp"
#include "blst/avx/tests/avx_test_utils.h"
#include "blst/avx/test_data/test_values.hpp"

#include <array>
#include <string>

namespace {

using blst_avx_test::AvxFp2Ring;
using blst_avx_test::AvxRing;

template<int Batch>
void expect_batch_reduce_expand_matches_scalar(const AvxRing& ring, uint64_t seed)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);

    std::array<BigInt, Batch> left_bigints = {};
    std::array<BigInt, Batch> right_bigints = {};
    std::array<AvxRing::ReadyElement<AvxRing::MAX_ADD>, Batch> ready_elements = {};

    for (int i = 0; i < Batch; ++i) {
        left_bigints[i] = BigInt(rng.get_z_range(modulus.get_mpz()));
        right_bigints[i] = BigInt(rng.get_z_range(modulus.get_mpz()));
        const auto left = ring.from_bigint(left_bigints[i]);
        const auto right = ring.from_bigint(right_bigints[i]);
        ready_elements[i] = ring.template ready<AvxRing::MAX_ADD>(left * right);
    }

    const auto result = ring.batch_reduce_expand<Batch>(ready_elements);

    for (int i = 0; i < Batch; ++i) {
        SCOPED_TRACE("batch=" + std::to_string(Batch) + " lane=" + std::to_string(i));

        const blst::Fp left_scalar = blst_avx_test::bigint_to_scalar_fp(left_bigints[i]);
        const blst::Fp right_scalar = blst_avx_test::bigint_to_scalar_fp(right_bigints[i]);
        blst::Fp scalar_result = {};
        blst::Fp::mul(scalar_result, left_scalar, right_scalar);

        EXPECT_EQ(
            ring.to_bigint(result[i]).to_string(),
            blst_avx_test::scalar_fp_to_bigint(scalar_result).to_string());
    }
}

TEST(AvxBaseField, FromBigIntMatchesKnownRnsVectors)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());

    const auto a_element = ring.from_bigint(BigInt(::a, 10));
    const auto b_element = ring.from_bigint(BigInt(::b, 10));
    const auto c_element = ring.from_bigint(BigInt(::c, 10));

    EXPECT_EQ(a_element.m1.data.to_array(), a_m1);
    EXPECT_EQ(a_element.m2.data.to_array(), a_m2);
    EXPECT_EQ(b_element.m1.data.to_array(), b_m1);
    EXPECT_EQ(b_element.m2.data.to_array(), b_m2);
    EXPECT_EQ(c_element.m1.data.to_array(), c_m1);
    EXPECT_EQ(c_element.m2.data.to_array(), c_m2);
}

TEST(AvxBaseField, ModMulMatchesKnownAnswer)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());

    const auto result = ring.modmul(
        ring.from_bigint(BigInt(::a, 10)),
        ring.from_bigint(BigInt(::b, 10)));

    EXPECT_EQ(result.m1.data.to_array(), c_m1);
    EXPECT_EQ(result.m2.data.to_array(), c_m2);
    EXPECT_EQ(ring.to_bigint(result).to_string(), BigInt(::c_correct, 10).to_string());
}

TEST(AvxBaseField, BatchModMul8MatchesScalarBlstFp)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());
    const BigInt modulus = blst_avx_test::modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);

    for (int round = 0; round < 8; ++round) {
        std::array<AvxRing::StandardElement, 8> left_vroom;
        std::array<AvxRing::StandardElement, 8> right_vroom;
        std::array<BigInt, 8> left_bigints;
        std::array<BigInt, 8> right_bigints;

        for (int i = 0; i < 8; ++i) {
            left_bigints[i] = BigInt(rng.get_z_range(modulus.get_mpz()));
            right_bigints[i] = BigInt(rng.get_z_range(modulus.get_mpz()));
            left_vroom[i] = ring.from_bigint(left_bigints[i]);
            right_vroom[i] = ring.from_bigint(right_bigints[i]);
        }

        const auto vroom_result = ring.batch_modmul<8>(left_vroom, right_vroom);

        for (int i = 0; i < 8; ++i) {
            SCOPED_TRACE("round=" + std::to_string(round) + " lane=" + std::to_string(i));

            const blst::Fp left_scalar = blst_avx_test::bigint_to_scalar_fp(left_bigints[i]);
            const blst::Fp right_scalar = blst_avx_test::bigint_to_scalar_fp(right_bigints[i]);
            blst::Fp scalar_result = {};
            blst::Fp::mul(scalar_result, left_scalar, right_scalar);

            EXPECT_EQ(ring.to_bigint(vroom_result[i]).to_string(),
                blst_avx_test::scalar_fp_to_bigint(scalar_result).to_string());
        }
    }
}

TEST(AvxBaseField, AddChainInverseMatchesScalarBlstFp)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const BLS381AddChainInversion<AvxRing> inverter(ring);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);

    for (int i = 0; i < 8; ++i) {
        BigInt value(rng.get_z_range(modulus.get_mpz()));
        if (value == BigInt(0)) {
            value = BigInt(1);
        }

        const auto element = ring.from_bigint(value);
        const auto result = inverter.invert(element, ring);

        const blst::Fp scalar_value = blst_avx_test::bigint_to_scalar_fp(value);
        blst::Fp scalar_inverse = {};
        blst::Fp::reciprocal(scalar_inverse, scalar_value);

        SCOPED_TRACE("case=" + std::to_string(i));
        EXPECT_EQ(
            ring.to_bigint(result).to_string(),
            blst_avx_test::scalar_fp_to_bigint(scalar_inverse).to_string());
    }
}

TEST(AvxBaseField, BatchReduceExpandMatchesScalarBlstFp)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());

    expect_batch_reduce_expand_matches_scalar<1>(ring, 42);
    expect_batch_reduce_expand_matches_scalar<2>(ring, 42);
    expect_batch_reduce_expand_matches_scalar<3>(ring, 42);
    expect_batch_reduce_expand_matches_scalar<4>(ring, 42);
    expect_batch_reduce_expand_matches_scalar<5>(ring, 42);
    expect_batch_reduce_expand_matches_scalar<6>(ring, 42);
    expect_batch_reduce_expand_matches_scalar<12>(ring, 42);
}

TEST(AvxBaseField, Fp2MulMatchesScalarBlstFp2)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());
    const AvxFp2Ring fp2ring(ring);
    const BigInt modulus = blst_avx_test::modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(123);

    for (int i = 0; i < 8; ++i) {
        const BigInt ax(rng.get_z_range(modulus.get_mpz()));
        const BigInt ay(rng.get_z_range(modulus.get_mpz()));
        const BigInt bx(rng.get_z_range(modulus.get_mpz()));
        const BigInt by(rng.get_z_range(modulus.get_mpz()));

        const auto vroom_left = blst_avx_test::make_avx_fp2(ring, ax, ay);
        const auto vroom_right = blst_avx_test::make_avx_fp2(ring, bx, by);
        const auto vroom_result = fp2ring.mul(vroom_left, vroom_right)[0];

        const blst::Fp2 scalar_left = blst_avx_test::bigint_pair_to_scalar_fp2(ax, ay);
        const blst::Fp2 scalar_right = blst_avx_test::bigint_pair_to_scalar_fp2(bx, by);
        blst::Fp2 scalar_result = {};
        blst::Fp2::mul(scalar_result, scalar_left, scalar_right);
        const auto [scalar_x, scalar_y] = blst_avx_test::scalar_fp2_to_bigints(scalar_result);

        SCOPED_TRACE("case=" + std::to_string(i));
        EXPECT_EQ(ring.to_bigint(vroom_result.x()).to_string(), scalar_x.to_string());
        EXPECT_EQ(ring.to_bigint(vroom_result.y()).to_string(), scalar_y.to_string());
    }
}

TEST(AvxBaseField, Fp2SquareMatchesScalarBlstFp2)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());
    const AvxFp2Ring fp2ring(ring);
    const BigInt modulus = blst_avx_test::modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(321);

    for (int i = 0; i < 8; ++i) {
        const BigInt ax(rng.get_z_range(modulus.get_mpz()));
        const BigInt ay(rng.get_z_range(modulus.get_mpz()));

        const auto vroom_value = blst_avx_test::make_avx_fp2(ring, ax, ay);
        const auto vroom_result = fp2ring.square(vroom_value)[0];

        const blst::Fp2 scalar_value = blst_avx_test::bigint_pair_to_scalar_fp2(ax, ay);
        blst::Fp2 scalar_result = {};
        blst::Fp2::sqr(scalar_result, scalar_value);
        const auto [scalar_x, scalar_y] = blst_avx_test::scalar_fp2_to_bigints(scalar_result);

        SCOPED_TRACE("case=" + std::to_string(i));
        EXPECT_EQ(ring.to_bigint(vroom_result.x()).to_string(), scalar_x.to_string());
        EXPECT_EQ(ring.to_bigint(vroom_result.y()).to_string(), scalar_y.to_string());
    }
}

TEST(AvxBaseField, Fp2EdgeValuesAndIdentitiesMatchScalarBlstFp2)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());
    const AvxFp2Ring fp2ring(ring);
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const BigInt zero(0);
    const BigInt one(1);
    const BigInt minus_one = modulus - one;
    const BigInt max_reduced = modulus - one;

    const std::array<std::pair<BigInt, BigInt>, 6> cases = {{
        {zero, zero},
        {one, zero},
        {zero, one},
        {minus_one, zero},
        {zero, minus_one},
        {max_reduced, max_reduced},
    }};

    const auto vroom_one = fp2ring.one();
    const auto vroom_zero = blst_avx_test::make_avx_fp2(ring, zero, zero);
    const blst::Fp2 scalar_one = blst_avx_test::bigint_pair_to_scalar_fp2(one, zero);
    const blst::Fp2 scalar_zero = blst_avx_test::bigint_pair_to_scalar_fp2(zero, zero);

    for (size_t i = 0; i < cases.size(); ++i) {
        const auto& [x, y] = cases[i];
        const auto vroom_value = blst_avx_test::make_avx_fp2(ring, x, y);
        const blst::Fp2 scalar_value = blst_avx_test::bigint_pair_to_scalar_fp2(x, y);

        const auto vroom_square = fp2ring.square(vroom_value)[0];
        blst::Fp2 scalar_square = {};
        blst::Fp2::sqr(scalar_square, scalar_value);
        const auto [scalar_sqr_x, scalar_sqr_y] = blst_avx_test::scalar_fp2_to_bigints(scalar_square);

        const auto vroom_mul_one = fp2ring.mul(vroom_value, vroom_one)[0];
        blst::Fp2 scalar_mul_one = {};
        blst::Fp2::mul(scalar_mul_one, scalar_value, scalar_one);
        const auto [scalar_mul_one_x, scalar_mul_one_y] = blst_avx_test::scalar_fp2_to_bigints(scalar_mul_one);

        const auto vroom_mul_zero = fp2ring.mul(vroom_value, vroom_zero)[0];
        blst::Fp2 scalar_mul_zero = {};
        blst::Fp2::mul(scalar_mul_zero, scalar_value, scalar_zero);
        const auto [scalar_mul_zero_x, scalar_mul_zero_y] = blst_avx_test::scalar_fp2_to_bigints(scalar_mul_zero);

        SCOPED_TRACE("case=" + std::to_string(i));
        EXPECT_EQ(ring.to_bigint(vroom_square.x()).to_string(), scalar_sqr_x.to_string());
        EXPECT_EQ(ring.to_bigint(vroom_square.y()).to_string(), scalar_sqr_y.to_string());
        EXPECT_EQ(ring.to_bigint(vroom_mul_one.x()).to_string(), scalar_mul_one_x.to_string());
        EXPECT_EQ(ring.to_bigint(vroom_mul_one.y()).to_string(), scalar_mul_one_y.to_string());
        EXPECT_EQ(ring.to_bigint(vroom_mul_zero.x()).to_string(), scalar_mul_zero_x.to_string());
        EXPECT_EQ(ring.to_bigint(vroom_mul_zero.y()).to_string(), scalar_mul_zero_y.to_string());
    }
}

}  // namespace
