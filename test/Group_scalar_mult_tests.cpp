#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "bindings/blst.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/map_to_G.h"
#include "blst/scalar/curve/point.h"
#include "blst/scalar/fields/Fp12.h"

namespace {

std::array<blst_byte, 32> test_scalar(size_t nbits)
{
    std::array<blst_byte, 32> scalar = {};
    for (size_t i = 0; i < scalar.size(); ++i) {
        scalar[i] = static_cast<blst_byte>(0xa5U ^ static_cast<unsigned>(i * 29U));
    }
    if (nbits == 0) {
        scalar.fill(0);
        return scalar;
    }

    const size_t nbytes = (nbits + 7) / 8;
    for (size_t i = nbytes; i < scalar.size(); ++i) {
        scalar[i] = 0;
    }
    if ((nbits & 7U) != 0) {
        scalar[nbytes - 1] &= static_cast<blst_byte>((1U << (nbits & 7U)) - 1U);
    }
    return scalar;
}

std::array<blst_byte, 32> scalar_bytes(const std::array<uint64_t, 4>& limbs)
{
    std::array<blst_byte, 32> scalar = {};
    for (size_t limb = 0; limb < limbs.size(); ++limb) {
        for (size_t byte = 0; byte < sizeof(uint64_t); ++byte) {
            scalar[limb * sizeof(uint64_t) + byte] =
                static_cast<blst_byte>(limbs[limb] >> (byte * 8));
        }
    }
    return scalar;
}

constexpr std::array<uint64_t, 4> kScalarOrder = {
    0xffffffff00000001ULL,
    0x53bda402fffe5bfeULL,
    0x3339d80809a1d805ULL,
    0x73eda753299d7d48ULL,
};

constexpr unsigned hex_nibble(char value)
{
    return value >= '0' && value <= '9' ? static_cast<unsigned>(value - '0')
        : value >= 'a' && value <= 'f'
            ? static_cast<unsigned>(value - 'a' + 10)
            : static_cast<unsigned>(value - 'A' + 10);
}

template<size_t Chars>
constexpr auto decode_hex(const char (&encoded)[Chars])
{
    static_assert((Chars - 1) % 2 == 0);
    std::array<blst_byte, (Chars - 1) / 2> decoded = {};
    for (size_t i = 0; i < decoded.size(); ++i) {
        decoded[i] = static_cast<blst_byte>(
            (hex_nibble(encoded[2 * i]) << 4) |
            hex_nibble(encoded[2 * i + 1]));
    }
    return decoded;
}

constexpr auto kG1Order11Encoded = decode_hex(
    "80803f0d09fec09a95f2ee7495323c15c162270c7cceaffa"
    "8566e941c66bcf206e72955d58b3b32e564de3209d672ca5");
constexpr auto kG2Order13Encoded = decode_hex(
    "808514e6f6c41ddb685cd8e351d7cb2ee891d54d5850d42c"
    "440365b15dd8ab06b53da9ebb22bc31cea2234d30499b09f"
    "01b290cae398ecc1fda3d4beb248af58bd351c3d99a5b6c7"
    "70fdb54bf46b123261b55cf0762aeef2341e35e90608fc31");

TEST(GroupScalarMult, ScalarGtConditionalNegation)
{
    using ScalarGT = blst::scalar::GT;

    const ScalarGT original = ScalarGT::generator();
    ScalarGT expected_negation = {};
    ScalarGT::neg(expected_negation, original);

    ScalarGT unchanged = original;
    ScalarGT::cneg(unchanged, 0);
    EXPECT_TRUE(ScalarGT::equal(unchanged, original));

    ScalarGT negated = original;
    ScalarGT::cneg(negated, 1);
    EXPECT_TRUE(ScalarGT::equal(negated, expected_negation));
}

template<class Group>
Group reference_multiply(
    const Group& point, const blst_byte scalar[32], size_t nbits)
{
	Group accumulator = Group::infinity();
	Group power = point;
	for (size_t bit = 0; bit < nbits; ++bit) {
		if ((scalar[bit / 8] >> (bit % 8)) & 1U) {
			Group sum = {};
			Group::dadd(sum, accumulator, power);
			accumulator = sum;
		}
		Group doubled = {};
		Group::doubling(doubled, power);
		power = doubled;
	}
	return accumulator;
}

blst_fp12 reference_gt_multiply(
    const blst_fp12& point, const blst_byte scalar[32], size_t nbits)
{
    blst_fp12 accumulator = *blst_fp12_one();
    blst_fp12 power = point;
    for (size_t bit = 0; bit < nbits; ++bit) {
        if ((scalar[bit / 8] >> (bit % 8)) & 1U) {
            blst_fp12 sum = {};
            blst_fp12_mul(&sum, &accumulator, &power);
            accumulator = sum;
        }
        blst_fp12 doubled = {};
        blst_fp12_cyclotomic_sqr(&doubled, &power);
        power = doubled;
    }
    return accumulator;
}

void expect_single_group_multiplication_matches_oracle(
    const std::array<blst_byte, 32>& scalar, size_t nbits,
    const char* case_name)
{
    ASSERT_LE(nbits, 256U);

    const blst_p1 g1_expected = reference_multiply(
        *blst_p1_generator(), scalar.data(), nbits);
    blst_p1 g1_ct = {};
    blst_p1 g1_vartime = {};
    blst_p1 g1_unchecked = {};
    blst_p1_mult(&g1_ct, blst_p1_generator(), scalar.data(), nbits);
    blst_p1_mult_vartime(
        &g1_vartime, blst_p1_generator(), scalar.data(), nbits);
    blst_p1_unchecked_mult(
        &g1_unchecked, blst_p1_generator(), scalar.data(), nbits);
    EXPECT_TRUE(blst_p1_is_equal(&g1_ct, &g1_expected)) << case_name;
    EXPECT_TRUE(blst_p1_is_equal(&g1_vartime, &g1_expected)) << case_name;
    EXPECT_TRUE(blst_p1_is_equal(&g1_unchecked, &g1_expected)) << case_name;

    const blst_p2 g2_expected = reference_multiply(
        *blst_p2_generator(), scalar.data(), nbits);
    blst_p2 g2_ct = {};
    blst_p2 g2_vartime = {};
    blst_p2 g2_unchecked = {};
    blst_p2_mult(&g2_ct, blst_p2_generator(), scalar.data(), nbits);
    blst_p2_mult_vartime(
        &g2_vartime, blst_p2_generator(), scalar.data(), nbits);
    blst_p2_unchecked_mult(
        &g2_unchecked, blst_p2_generator(), scalar.data(), nbits);
    EXPECT_TRUE(blst_p2_is_equal(&g2_ct, &g2_expected)) << case_name;
    EXPECT_TRUE(blst_p2_is_equal(&g2_vartime, &g2_expected)) << case_name;
    EXPECT_TRUE(blst_p2_is_equal(&g2_unchecked, &g2_expected)) << case_name;

    const blst_fp12 gt_expected = reference_gt_multiply(
        *blst_gt_generator(), scalar.data(), nbits);
    blst_fp12 gt_ct = {};
    blst_fp12 gt_vartime = {};
    blst_gt_mult(&gt_ct, blst_gt_generator(), scalar.data(), nbits);
    blst_gt_mult_vartime(
        &gt_vartime, blst_gt_generator(), scalar.data(), nbits);
    EXPECT_TRUE(blst_gt_is_equal(&gt_ct, &gt_expected)) << case_name;
    EXPECT_TRUE(blst_gt_is_equal(&gt_vartime, &gt_expected)) << case_name;
}

std::array<std::array<uint64_t, 4>, 8> full_width_boundary_scalars()
{
    auto order_minus_one = kScalarOrder;
    --order_minus_one[0];
    auto order_plus_one = kScalarOrder;
    ++order_plus_one[0];

    std::array<uint64_t, 4> twice_order = {};
    uint64_t carry = 0;
    for (size_t i = 0; i < twice_order.size(); ++i) {
        const uint64_t doubled = kScalarOrder[i] + kScalarOrder[i];
        const uint64_t doubled_carry = doubled < kScalarOrder[i];
        const uint64_t with_carry = doubled + carry;
        const uint64_t add_carry = with_carry < doubled;
        twice_order[i] = with_carry;
        carry = doubled_carry | add_carry;
    }
    EXPECT_EQ(carry, 0U);

    auto twice_order_minus_one = twice_order;
    --twice_order_minus_one[0];

    return {{
        {{0, 0, 0, 0}},
        {{1, 0, 0, 0}},
        order_minus_one,
        kScalarOrder,
        order_plus_one,
        twice_order_minus_one,
        twice_order,
        {{UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX}},
    }};
}

TEST(GroupScalarMult, ConstantTimeAndVariableTimeApisAgree)
{
    const std::array<size_t, 6> bit_lengths = {0, 1, 3, 128, 255, 256};
    for (const size_t nbits : bit_lengths) {
        const auto scalar = test_scalar(nbits);

        blst_p1 g1_ct = {};
        blst_p1 g1_vartime = {};
        blst_p1_mult(&g1_ct, blst_p1_generator(), scalar.data(), nbits);
        blst_p1_mult_vartime(&g1_vartime, blst_p1_generator(), scalar.data(), nbits);
        EXPECT_TRUE(blst_p1_is_equal(&g1_ct, &g1_vartime)) << "G1 nbits=" << nbits;

        blst_p2 g2_ct = {};
        blst_p2 g2_vartime = {};
        blst_p2_mult(&g2_ct, blst_p2_generator(), scalar.data(), nbits);
        blst_p2_mult_vartime(&g2_vartime, blst_p2_generator(), scalar.data(), nbits);
        EXPECT_TRUE(blst_p2_is_equal(&g2_ct, &g2_vartime)) << "G2 nbits=" << nbits;

        blst_fp12 gt_ct = {};
        blst_fp12 gt_vartime = {};
        blst_gt_mult(&gt_ct, blst_gt_generator(), scalar.data(), nbits);
        blst_gt_mult_vartime(&gt_vartime, blst_gt_generator(), scalar.data(), nbits);
        EXPECT_TRUE(blst_gt_is_equal(&gt_ct, &gt_vartime)) << "GT nbits=" << nbits;
    }

    blst_p1 g1_zero = {};
    blst_p2 g2_zero = {};
    blst_fp12 gt_zero = {};
    blst_p1_mult(&g1_zero, blst_p1_generator(), nullptr, 0);
    blst_p2_mult(&g2_zero, blst_p2_generator(), nullptr, 0);
    blst_gt_mult(&gt_zero, blst_gt_generator(), nullptr, 0);
    EXPECT_TRUE(blst_p1_is_inf(&g1_zero));
    EXPECT_TRUE(blst_p2_is_inf(&g2_zero));
    EXPECT_TRUE(blst_fp12_is_one(&gt_zero));
}

TEST(GroupScalarMult, DebugContractsRejectInvalidInputs)
{
#if !defined(NDEBUG) && GTEST_HAS_DEATH_TEST
    const auto scalar = test_scalar(256);
    blst_p1 g1_out = {};
    blst_p2 g2_out = {};
    blst_fp12 gt_out = {};

    EXPECT_DEATH(
        blst_p1_mult(&g1_out, blst_p1_generator(), nullptr, 1), "");
    EXPECT_DEATH(
        blst_p2_mult(&g2_out, blst_p2_generator(), scalar.data(), 257), "");

    blst_p1 off_curve_g1 = *blst_p1_generator();
    off_curve_g1.x[0] ^= 1;
    EXPECT_DEATH(
        blst_p1_mult(&g1_out, &off_curve_g1, scalar.data(), 256), "");
    EXPECT_DEATH(
        blst_p1_unchecked_mult(&g1_out, &off_curve_g1, scalar.data(), 256), "");

    blst_p2 off_curve_g2 = *blst_p2_generator();
    off_curve_g2.x[0][0] ^= 1;
    EXPECT_DEATH(
        blst_p2_mult(&g2_out, &off_curve_g2, scalar.data(), 256), "");
    EXPECT_DEATH(
        blst_p2_unchecked_mult(&g2_out, &off_curve_g2, scalar.data(), 256), "");

    const blst_fp12 invalid_gt = {};
    EXPECT_DEATH(
        blst_gt_mult(&gt_out, &invalid_gt, scalar.data(), 256), "");
#endif
}

TEST(GroupScalarMult, DefaultApisSupportExactOutputAliasing)
{
    const auto scalar = test_scalar(256);

    blst_p1 g1_expected = {};
    blst_p1 g1_alias = *blst_p1_generator();
    blst_p1_mult(&g1_expected, &g1_alias, scalar.data(), 256);
    blst_p1_mult(&g1_alias, &g1_alias, scalar.data(), 256);
    EXPECT_TRUE(blst_p1_is_equal(&g1_expected, &g1_alias));
    g1_alias = *blst_p1_generator();
    blst_p1_mult_vartime(&g1_expected, &g1_alias, scalar.data(), 256);
    blst_p1_mult_vartime(&g1_alias, &g1_alias, scalar.data(), 256);
    EXPECT_TRUE(blst_p1_is_equal(&g1_expected, &g1_alias));
    g1_alias = *blst_p1_generator();
    blst_p1_unchecked_mult(&g1_expected, &g1_alias, scalar.data(), 256);
    blst_p1_unchecked_mult(&g1_alias, &g1_alias, scalar.data(), 256);
    EXPECT_TRUE(blst_p1_is_equal(&g1_expected, &g1_alias));

    blst_p2 g2_expected = {};
    blst_p2 g2_alias = *blst_p2_generator();
    blst_p2_mult(&g2_expected, &g2_alias, scalar.data(), 256);
    blst_p2_mult(&g2_alias, &g2_alias, scalar.data(), 256);
    EXPECT_TRUE(blst_p2_is_equal(&g2_expected, &g2_alias));
    g2_alias = *blst_p2_generator();
    blst_p2_mult_vartime(&g2_expected, &g2_alias, scalar.data(), 256);
    blst_p2_mult_vartime(&g2_alias, &g2_alias, scalar.data(), 256);
    EXPECT_TRUE(blst_p2_is_equal(&g2_expected, &g2_alias));
    g2_alias = *blst_p2_generator();
    blst_p2_unchecked_mult(&g2_expected, &g2_alias, scalar.data(), 256);
    blst_p2_unchecked_mult(&g2_alias, &g2_alias, scalar.data(), 256);
    EXPECT_TRUE(blst_p2_is_equal(&g2_expected, &g2_alias));

    blst_fp12 gt_expected = {};
    blst_fp12 gt_alias = *blst_gt_generator();
    blst_gt_mult(&gt_expected, &gt_alias, scalar.data(), 256);
    blst_gt_mult(&gt_alias, &gt_alias, scalar.data(), 256);
    EXPECT_TRUE(blst_gt_is_equal(&gt_expected, &gt_alias));
    gt_alias = *blst_gt_generator();
    blst_gt_mult_vartime(&gt_expected, &gt_alias, scalar.data(), 256);
    blst_gt_mult_vartime(&gt_alias, &gt_alias, scalar.data(), 256);
    EXPECT_TRUE(blst_gt_is_equal(&gt_expected, &gt_alias));
}

TEST(GroupScalarMult, BitsAboveNbitsAreIgnored)
{
    for (const size_t nbits : std::array<size_t, 4>{{1, 127, 129, 255}}) {
        auto clean = test_scalar(nbits);
        auto dirty = clean;
        const size_t top_byte = nbits / 8;
        const unsigned meaningful_top_bits = static_cast<unsigned>(nbits & 7U);
        ASSERT_NE(meaningful_top_bits, 0U);
        dirty[top_byte] |= static_cast<blst_byte>(
            ~static_cast<unsigned>((1U << meaningful_top_bits) - 1U));
        for (size_t i = top_byte + 1; i < dirty.size(); ++i) {
            dirty[i] = 0xff;
        }

        blst_p1 g1_clean = {};
        blst_p1 g1_dirty = {};
        blst_p1_mult(&g1_clean, blst_p1_generator(), clean.data(), nbits);
        blst_p1_mult(&g1_dirty, blst_p1_generator(), dirty.data(), nbits);
        EXPECT_TRUE(blst_p1_is_equal(&g1_clean, &g1_dirty)) << "G1 nbits=" << nbits;

        blst_p2 g2_clean = {};
        blst_p2 g2_dirty = {};
        blst_p2_mult(&g2_clean, blst_p2_generator(), clean.data(), nbits);
        blst_p2_mult(&g2_dirty, blst_p2_generator(), dirty.data(), nbits);
        EXPECT_TRUE(blst_p2_is_equal(&g2_clean, &g2_dirty)) << "G2 nbits=" << nbits;

        blst_fp12 gt_clean = {};
        blst_fp12 gt_dirty = {};
        blst_gt_mult(&gt_clean, blst_gt_generator(), clean.data(), nbits);
        blst_gt_mult(&gt_dirty, blst_gt_generator(), dirty.data(), nbits);
        EXPECT_TRUE(blst_gt_is_equal(&gt_clean, &gt_dirty)) << "GT nbits=" << nbits;
    }
}

TEST(GroupScalarMult, FullWidthBoundaryScalarsAgree)
{
    for (const auto& limbs : full_width_boundary_scalars()) {
        const auto scalar = scalar_bytes(limbs);
		blst_p1 g1_ct = {};
		blst_p1 g1_vartime = {};
		blst_p1 g1_unchecked = {};
		blst_p1_mult(&g1_ct, blst_p1_generator(), scalar.data(), 256);
		blst_p1_mult_vartime(&g1_vartime, blst_p1_generator(), scalar.data(), 256);
		blst_p1_unchecked_mult(&g1_unchecked, blst_p1_generator(), scalar.data(), 256);
		EXPECT_TRUE(blst_p1_is_equal(&g1_ct, &g1_vartime));
		EXPECT_TRUE(blst_p1_is_equal(&g1_ct, &g1_unchecked));

		blst_p2 g2_ct = {};
		blst_p2 g2_vartime = {};
		blst_p2 g2_unchecked = {};
		blst_p2_mult(&g2_ct, blst_p2_generator(), scalar.data(), 256);
		blst_p2_mult_vartime(&g2_vartime, blst_p2_generator(), scalar.data(), 256);
		blst_p2_unchecked_mult(&g2_unchecked, blst_p2_generator(), scalar.data(), 256);
		EXPECT_TRUE(blst_p2_is_equal(&g2_ct, &g2_vartime));
		EXPECT_TRUE(blst_p2_is_equal(&g2_ct, &g2_unchecked));

        blst_fp12 gt_ct = {};
        blst_fp12 gt_vartime = {};
        blst_gt_mult(&gt_ct, blst_gt_generator(), scalar.data(), 256);
        blst_gt_mult_vartime(&gt_vartime, blst_gt_generator(), scalar.data(), 256);
        EXPECT_TRUE(blst_gt_is_equal(&gt_ct, &gt_vartime));
    }
}

TEST(GroupScalarMult, FixedScheduleGtSplitProducesUnsignedRadixDigits)
{
    auto cases = full_width_boundary_scalars();
    uint64_t state = 0x6a09e667f3bcc909ULL;
    for (size_t iteration = 0; iteration < 32; ++iteration) {
        std::array<uint64_t, 4> scalar = {};
        for (auto& limb : scalar) {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            limb = state;
        }

        uint64_t actual[4] = {};
        blst::scalar::Fp12::gt_split_scalar_4d_ct(actual, scalar.data());
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_LT(actual[i], 0xd201000000010000ULL)
                << "random iteration=" << iteration << " digit=" << i;
        }
    }

    for (size_t case_index = 0; case_index < cases.size(); ++case_index) {
        uint64_t actual[4] = {};
        blst::scalar::Fp12::gt_split_scalar_4d_ct(actual, cases[case_index].data());
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_LT(actual[i], 0xd201000000010000ULL)
                << "boundary case=" << case_index << " digit=" << i;
        }
    }

	uint64_t reduced[4] = {};
	blst::scalar::Fp12::gt_split_scalar_4d_ct(reduced, cases[0].data());
	EXPECT_EQ(reduced[0] | reduced[1] | reduced[2] | reduced[3], 0U);
	blst::scalar::Fp12::gt_split_scalar_4d_ct(reduced, cases[1].data());
	EXPECT_EQ(reduced[0], 1U);
	EXPECT_EQ(reduced[1] | reduced[2] | reduced[3], 0U);
	blst::scalar::Fp12::gt_split_scalar_4d_ct(reduced, cases[3].data());
	EXPECT_EQ(reduced[0] | reduced[1] | reduced[2] | reduced[3], 0U);
	blst::scalar::Fp12::gt_split_scalar_4d_ct(reduced, cases[4].data());
	EXPECT_EQ(reduced[0], 1U);
	EXPECT_EQ(reduced[1] | reduced[2] | reduced[3], 0U);
}

TEST(GroupScalarMult, UncheckedMultiplicationSupportsOnCurveNonSubgroupPoints)
{
	using blst::scalar::G1AffinePoint;
	using blst::scalar::G1Point;
	using blst::scalar::G2AffinePoint;
	using blst::scalar::G2Point;

	G1Point g1_isogenous = {};
	G1Point g1 = {};
	blst::scalar::map_to_isogenous<G1Point>(g1_isogenous, G1AffinePoint::generator().x);
	blst::scalar::isogeny_map_to<G1Point>(g1, g1_isogenous);
	ASSERT_TRUE(G1Point::on_curve(g1));
	ASSERT_FALSE(G1Point::in_subgroup(g1));

	G2Point g2_isogenous = {};
	G2Point g2 = {};
	blst::scalar::map_to_isogenous<G2Point>(g2_isogenous, G2AffinePoint::generator().x);
	blst::scalar::isogeny_map_to<G2Point>(g2, g2_isogenous);
	ASSERT_TRUE(G2Point::on_curve(g2));
	ASSERT_FALSE(G2Point::in_subgroup(g2));

	const auto scalar = scalar_bytes({{19, 0, 0, 0}});
	G1Point g1_expected = reference_multiply(g1, scalar.data(), 256);
	G1Point g1_actual = {};
	blst_p1_unchecked_mult(&g1_actual, &g1, scalar.data(), 256);
	EXPECT_TRUE(G1Point::equal(g1_expected, g1_actual));

	G2Point g2_expected = reference_multiply(g2, scalar.data(), 256);
	G2Point g2_actual = {};
	blst_p2_unchecked_mult(&g2_actual, &g2, scalar.data(), 256);
	EXPECT_TRUE(G2Point::equal(g2_expected, g2_actual));

#if !defined(NDEBUG) && GTEST_HAS_DEATH_TEST
	EXPECT_DEATH({
		blst_p1 rejected = {};
		blst_p1_mult(&rejected, &g1, scalar.data(), 256);
	}, "BLST_ASSERT");
	EXPECT_DEATH({
		blst_p2 rejected = {};
		blst_p2_mult(&rejected, &g2, scalar.data(), 256);
	}, "BLST_ASSERT");
#endif
}

TEST(GroupScalarMult, UncheckedMultiplicationSupportsLowOrderPoints)
{
    blst_p1_affine g1_affine = {};
    ASSERT_EQ(
        blst_p1_uncompress(&g1_affine, kG1Order11Encoded.data()),
        BLST_SUCCESS);
    ASSERT_TRUE(blst_p1_affine_on_curve(&g1_affine));
    ASSERT_FALSE(blst_p1_affine_in_g1(&g1_affine));
    blst_p1 g1 = {};
    blst_p1_from_affine(&g1, &g1_affine);

    std::array<blst_byte, 32> scalar = {};
    scalar[0] = 11;
    blst_p1 g1_order = {};
    blst_p1_unchecked_mult(&g1_order, &g1, scalar.data(), 256);
    EXPECT_TRUE(blst_p1_is_inf(&g1_order));
    scalar[0] = 12;
    blst_p1 g1_order_plus_one = {};
    blst_p1_unchecked_mult(&g1_order_plus_one, &g1, scalar.data(), 256);
    EXPECT_TRUE(blst_p1_is_equal(&g1_order_plus_one, &g1));

    blst_p2_affine g2_affine = {};
    ASSERT_EQ(
        blst_p2_uncompress(&g2_affine, kG2Order13Encoded.data()),
        BLST_SUCCESS);
    ASSERT_TRUE(blst_p2_affine_on_curve(&g2_affine));
    ASSERT_FALSE(blst_p2_affine_in_g2(&g2_affine));
    blst_p2 g2 = {};
    blst_p2_from_affine(&g2, &g2_affine);

    scalar.fill(0);
    scalar[0] = 13;
    blst_p2 g2_order = {};
    blst_p2_unchecked_mult(&g2_order, &g2, scalar.data(), 256);
    EXPECT_TRUE(blst_p2_is_inf(&g2_order));
    scalar[0] = 14;
    blst_p2 g2_order_plus_one = {};
    blst_p2_unchecked_mult(&g2_order_plus_one, &g2, scalar.data(), 256);
    EXPECT_TRUE(blst_p2_is_equal(&g2_order_plus_one, &g2));
}

TEST(GroupScalarMult, BasicCPrimitivesSupportOnCurveNonSubgroupPoints)
{
    std::array<blst_byte, 32> scalar = {};

    blst_p1_affine g1_affine = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_p1_uncompress(&g1_affine, kG1Order11Encoded.data()));
    ASSERT_TRUE(blst_p1_affine_on_curve(&g1_affine));
    ASSERT_FALSE(blst_p1_affine_in_g1(&g1_affine));
    EXPECT_FALSE(blst_p1_affine_is_inf(&g1_affine));

    blst_p1 g1 = {};
    blst_p1_from_affine(&g1, &g1_affine);
    ASSERT_TRUE(blst_p1_on_curve(&g1));
    ASSERT_FALSE(blst_p1_in_g1(&g1));
    EXPECT_FALSE(blst_p1_is_inf(&g1));

    scalar[0] = 2;
    blst_p1 g1_twice_expected = {};
    blst_p1_unchecked_mult(
        &g1_twice_expected, &g1, scalar.data(), scalar.size() * 8);
    blst_p1 g1_twice = {};
    blst_p1_double(&g1_twice, &g1);
    EXPECT_TRUE(blst_p1_is_equal(&g1_twice, &g1_twice_expected));
    blst_p1 g1_complete = {};
    blst_p1_add_or_double(&g1_complete, &g1, &g1);
    EXPECT_TRUE(blst_p1_is_equal(&g1_complete, &g1_twice_expected));
    blst_p1_add_or_double_affine(&g1_complete, &g1, &g1_affine);
    EXPECT_TRUE(blst_p1_is_equal(&g1_complete, &g1_twice_expected));

    scalar[0] = 3;
    blst_p1 g1_three_expected = {};
    blst_p1_unchecked_mult(
        &g1_three_expected, &g1, scalar.data(), scalar.size() * 8);
    blst_p1 g1_three = {};
    blst_p1_add(&g1_three, &g1, &g1_twice);
    EXPECT_TRUE(blst_p1_is_equal(&g1_three, &g1_three_expected));
    blst_p1_affine g1_twice_affine = {};
    blst_p1_to_affine(&g1_twice_affine, &g1_twice);
    blst_p1_add_affine(&g1_three, &g1, &g1_twice_affine);
    EXPECT_TRUE(blst_p1_is_equal(&g1_three, &g1_three_expected));

    blst_p1 g1_difference = {};
    blst_p1_sub(&g1_difference, &g1_twice, &g1);
    EXPECT_TRUE(blst_p1_is_equal(&g1_difference, &g1));
    blst_p1_sub_affine(&g1_difference, &g1_twice, &g1_affine);
    EXPECT_TRUE(blst_p1_is_equal(&g1_difference, &g1));
    blst_p1_sub(&g1_difference, &g1, &g1);
    EXPECT_TRUE(blst_p1_is_inf(&g1_difference));
    blst_p1_sub_affine(&g1_difference, &g1, &g1_affine);
    EXPECT_TRUE(blst_p1_is_inf(&g1_difference));

    blst_p1 g1_negated = g1;
    blst_p1_cneg(&g1_negated, true);
    EXPECT_TRUE(blst_p1_on_curve(&g1_negated));
    EXPECT_FALSE(blst_p1_in_g1(&g1_negated));
    blst_p1_cneg(&g1_negated, true);
    EXPECT_TRUE(blst_p1_is_equal(&g1_negated, &g1));
    blst_p1_affine g1_affine_negated = g1_affine;
    blst_p1_affine_cneg(&g1_affine_negated, true);
    blst_p1_affine_cneg(&g1_affine_negated, true);
    EXPECT_TRUE(blst_p1_affine_is_equal(&g1_affine_negated, &g1_affine));

    blst_p1 g1_normalized = {};
    blst_p1_from_jacobian(&g1_normalized, &g1_twice);
    EXPECT_TRUE(blst_p1_is_equal(&g1_normalized, &g1_twice));
    blst_p1_affine g1_roundtrip_affine = {};
    blst_p1_to_affine(&g1_roundtrip_affine, &g1_normalized);
    EXPECT_TRUE(blst_p1_affine_is_equal(
        &g1_roundtrip_affine, &g1_twice_affine));

    std::array<blst_byte, 96> g1_serialized = {};
    std::array<blst_byte, 96> g1_affine_serialized = {};
    blst_p1_serialize(g1_serialized.data(), &g1);
    blst_p1_affine_serialize(g1_affine_serialized.data(), &g1_affine);
    EXPECT_EQ(g1_serialized, g1_affine_serialized);
    blst_p1_affine g1_deserialized = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_p1_deserialize(&g1_deserialized, g1_serialized.data()));
    EXPECT_TRUE(blst_p1_affine_is_equal(&g1_deserialized, &g1_affine));
    EXPECT_FALSE(blst_p1_affine_in_g1(&g1_deserialized));

    std::array<blst_byte, 48> g1_compressed = {};
    std::array<blst_byte, 48> g1_affine_compressed = {};
    blst_p1_compress(g1_compressed.data(), &g1);
    blst_p1_affine_compress(g1_affine_compressed.data(), &g1_affine);
    EXPECT_EQ(g1_compressed, g1_affine_compressed);
    blst_p1_affine g1_uncompressed = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_p1_uncompress(&g1_uncompressed, g1_compressed.data()));
    EXPECT_TRUE(blst_p1_affine_is_equal(&g1_uncompressed, &g1_affine));
    EXPECT_FALSE(blst_p1_affine_in_g1(&g1_uncompressed));

    scalar.fill(0);
    blst_p2_affine g2_affine = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_p2_uncompress(&g2_affine, kG2Order13Encoded.data()));
    ASSERT_TRUE(blst_p2_affine_on_curve(&g2_affine));
    ASSERT_FALSE(blst_p2_affine_in_g2(&g2_affine));
    EXPECT_FALSE(blst_p2_affine_is_inf(&g2_affine));

    blst_p2 g2 = {};
    blst_p2_from_affine(&g2, &g2_affine);
    ASSERT_TRUE(blst_p2_on_curve(&g2));
    ASSERT_FALSE(blst_p2_in_g2(&g2));
    EXPECT_FALSE(blst_p2_is_inf(&g2));

    scalar[0] = 2;
    blst_p2 g2_twice_expected = {};
    blst_p2_unchecked_mult(
        &g2_twice_expected, &g2, scalar.data(), scalar.size() * 8);
    blst_p2 g2_twice = {};
    blst_p2_double(&g2_twice, &g2);
    EXPECT_TRUE(blst_p2_is_equal(&g2_twice, &g2_twice_expected));
    blst_p2 g2_complete = {};
    blst_p2_add_or_double(&g2_complete, &g2, &g2);
    EXPECT_TRUE(blst_p2_is_equal(&g2_complete, &g2_twice_expected));
    blst_p2_add_or_double_affine(&g2_complete, &g2, &g2_affine);
    EXPECT_TRUE(blst_p2_is_equal(&g2_complete, &g2_twice_expected));

    scalar[0] = 3;
    blst_p2 g2_three_expected = {};
    blst_p2_unchecked_mult(
        &g2_three_expected, &g2, scalar.data(), scalar.size() * 8);
    blst_p2 g2_three = {};
    blst_p2_add(&g2_three, &g2, &g2_twice);
    EXPECT_TRUE(blst_p2_is_equal(&g2_three, &g2_three_expected));
    blst_p2_affine g2_twice_affine = {};
    blst_p2_to_affine(&g2_twice_affine, &g2_twice);
    blst_p2_add_affine(&g2_three, &g2, &g2_twice_affine);
    EXPECT_TRUE(blst_p2_is_equal(&g2_three, &g2_three_expected));

    blst_p2 g2_difference = {};
    blst_p2_sub(&g2_difference, &g2_twice, &g2);
    EXPECT_TRUE(blst_p2_is_equal(&g2_difference, &g2));
    blst_p2_sub_affine(&g2_difference, &g2_twice, &g2_affine);
    EXPECT_TRUE(blst_p2_is_equal(&g2_difference, &g2));
    blst_p2_sub(&g2_difference, &g2, &g2);
    EXPECT_TRUE(blst_p2_is_inf(&g2_difference));
    blst_p2_sub_affine(&g2_difference, &g2, &g2_affine);
    EXPECT_TRUE(blst_p2_is_inf(&g2_difference));

    blst_p2 g2_negated = g2;
    blst_p2_cneg(&g2_negated, true);
    EXPECT_TRUE(blst_p2_on_curve(&g2_negated));
    EXPECT_FALSE(blst_p2_in_g2(&g2_negated));
    blst_p2_cneg(&g2_negated, true);
    EXPECT_TRUE(blst_p2_is_equal(&g2_negated, &g2));
    blst_p2_affine g2_affine_negated = g2_affine;
    blst_p2_affine_cneg(&g2_affine_negated, true);
    blst_p2_affine_cneg(&g2_affine_negated, true);
    EXPECT_TRUE(blst_p2_affine_is_equal(&g2_affine_negated, &g2_affine));

    blst_p2 g2_normalized = {};
    blst_p2_from_jacobian(&g2_normalized, &g2_twice);
    EXPECT_TRUE(blst_p2_is_equal(&g2_normalized, &g2_twice));
    blst_p2_affine g2_roundtrip_affine = {};
    blst_p2_to_affine(&g2_roundtrip_affine, &g2_normalized);
    EXPECT_TRUE(blst_p2_affine_is_equal(
        &g2_roundtrip_affine, &g2_twice_affine));

    std::array<blst_byte, 192> g2_serialized = {};
    std::array<blst_byte, 192> g2_affine_serialized = {};
    blst_p2_serialize(g2_serialized.data(), &g2);
    blst_p2_affine_serialize(g2_affine_serialized.data(), &g2_affine);
    EXPECT_EQ(g2_serialized, g2_affine_serialized);
    blst_p2_affine g2_deserialized = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_p2_deserialize(&g2_deserialized, g2_serialized.data()));
    EXPECT_TRUE(blst_p2_affine_is_equal(&g2_deserialized, &g2_affine));
    EXPECT_FALSE(blst_p2_affine_in_g2(&g2_deserialized));

    std::array<blst_byte, 96> g2_compressed = {};
    std::array<blst_byte, 96> g2_affine_compressed = {};
    blst_p2_compress(g2_compressed.data(), &g2);
    blst_p2_affine_compress(g2_affine_compressed.data(), &g2_affine);
    EXPECT_EQ(g2_compressed, g2_affine_compressed);
    blst_p2_affine g2_uncompressed = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_p2_uncompress(&g2_uncompressed, g2_compressed.data()));
    EXPECT_TRUE(blst_p2_affine_is_equal(&g2_uncompressed, &g2_affine));
    EXPECT_FALSE(blst_p2_affine_in_g2(&g2_uncompressed));
}

TEST(GroupScalarMult, PublicMembershipPredicatesRequireCurveAndSubgroup)
{
    using blst::scalar::G1AffinePoint;
    using blst::scalar::G1Point;
    using blst::scalar::G2AffinePoint;
    using blst::scalar::G2Point;

    EXPECT_TRUE(blst_p1_in_g1(blst_p1_generator()));
    EXPECT_TRUE(blst_p2_in_g2(blst_p2_generator()));

    G1Point g1_isogenous = {};
    G1Point g1_non_subgroup = {};
    blst::scalar::map_to_isogenous<G1Point>(
        g1_isogenous, G1AffinePoint::generator().x);
    blst::scalar::isogeny_map_to<G1Point>(g1_non_subgroup, g1_isogenous);
    ASSERT_TRUE(blst_p1_on_curve(&g1_non_subgroup));
    EXPECT_FALSE(blst_p1_in_g1(&g1_non_subgroup));

    G2Point g2_isogenous = {};
    G2Point g2_non_subgroup = {};
    blst::scalar::map_to_isogenous<G2Point>(
        g2_isogenous, G2AffinePoint::generator().x);
    blst::scalar::isogeny_map_to<G2Point>(g2_non_subgroup, g2_isogenous);
    ASSERT_TRUE(blst_p2_on_curve(&g2_non_subgroup));
    EXPECT_FALSE(blst_p2_in_g2(&g2_non_subgroup));

    G1Point g1_off_curve = *blst_p1_generator();
    g1_off_curve.x[0] ^= 1;
    EXPECT_FALSE(blst_p1_on_curve(&g1_off_curve));
    EXPECT_FALSE(blst_p1_in_g1(&g1_off_curve));

    G2Point g2_off_curve = *blst_p2_generator();
    g2_off_curve.x[0][0] ^= 1;
    EXPECT_FALSE(blst_p2_on_curve(&g2_off_curve));
    EXPECT_FALSE(blst_p2_in_g2(&g2_off_curve));

    G1AffinePoint g1_affine = {};
    G2AffinePoint g2_affine = {};
    G1Point::to_affine(g1_affine, g1_off_curve);
    G2Point::to_affine(g2_affine, g2_off_curve);
    EXPECT_FALSE(blst_p1_affine_in_g1(&g1_affine));
    EXPECT_FALSE(blst_p2_affine_in_g2(&g2_affine));

    const G1Point g1_infinity = G1Point::infinity();
    const G2Point g2_infinity = G2Point::infinity();
    EXPECT_TRUE(blst_p1_in_g1(&g1_infinity));
    EXPECT_TRUE(blst_p2_in_g2(&g2_infinity));
}

TEST(GroupScalarMult, ScalarPippengerUsesSubgroupModuloOrderContract)
{
    const auto scalar = scalar_bytes(
        {{UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX}});
    std::array<blst_byte, 64> scalars = {};
    std::copy(scalar.begin(), scalar.end(), scalars.begin());
    std::copy(scalar.begin(), scalar.end(), scalars.begin() + 32);
    const blst_byte* scalar_heads[] = {scalars.data(), nullptr};

    blst_p1 one_g1 = {};
    blst_p1_mult_vartime(
        &one_g1, blst_p1_generator(), scalar.data(), 256);
    blst_p1 expected_g1 = {};
    blst_p1_double(&expected_g1, &one_g1);
    std::array<blst_p1_affine, 2> g1_affines = {};
    blst_p1_to_affine(&g1_affines[0], blst_p1_generator());
    g1_affines[1] = g1_affines[0];
    const blst_p1_affine* g1_heads[] = {g1_affines.data(), nullptr};
    const size_t g1_scratch_bytes =
        blst_p1s_mult_pippenger_scratch_sizeof(g1_affines.size());
    std::vector<blst_limb_t> g1_scratch(
        (g1_scratch_bytes + sizeof(blst_limb_t) - 1) /
        sizeof(blst_limb_t));
    blst_p1 actual_g1 = {};
    blst_p1s_mult_pippenger(
        &actual_g1, g1_heads, g1_affines.size(), scalar_heads, 256,
        g1_scratch.data());
    EXPECT_TRUE(blst_p1_is_equal(&actual_g1, &expected_g1));

    blst_p2 one_g2 = {};
    blst_p2_mult_vartime(
        &one_g2, blst_p2_generator(), scalar.data(), 256);
    blst_p2 expected_g2 = {};
    blst_p2_double(&expected_g2, &one_g2);
    std::array<blst_p2_affine, 2> g2_affines = {};
    blst_p2_to_affine(&g2_affines[0], blst_p2_generator());
    g2_affines[1] = g2_affines[0];
    const blst_p2_affine* g2_heads[] = {g2_affines.data(), nullptr};
    const size_t g2_scratch_bytes =
        blst_p2s_mult_pippenger_scratch_sizeof(g2_affines.size());
    std::vector<blst_limb_t> g2_scratch(
        (g2_scratch_bytes + sizeof(blst_limb_t) - 1) /
        sizeof(blst_limb_t));
    blst_p2 actual_g2 = {};
    blst_p2s_mult_pippenger(
        &actual_g2, g2_heads, g2_affines.size(), scalar_heads, 256,
        g2_scratch.data());
    EXPECT_TRUE(blst_p2_is_equal(&actual_g2, &expected_g2));
}

TEST(GroupScalarMult, IndependentReferenceCoversBoundariesAndCrossovers)
{
    const auto boundaries = full_width_boundary_scalars();
    for (size_t i = 0; i < boundaries.size(); ++i) {
        const auto scalar = scalar_bytes(boundaries[i]);
        const std::string case_name = "boundary=" + std::to_string(i);
        expect_single_group_multiplication_matches_oracle(
            scalar, 256, case_name.c_str());
    }

    for (const size_t nbits : std::array<size_t, 8>{
             1, 127, 128, 143, 144, 175, 176, 255}) {
        const auto scalar = test_scalar(nbits);
        const std::string case_name = "nbits=" + std::to_string(nbits);
        expect_single_group_multiplication_matches_oracle(
            scalar, nbits, case_name.c_str());
    }

    uint64_t state = 0x243f6a8885a308d3ULL;
    for (size_t iteration = 0; iteration < 4; ++iteration) {
        std::array<uint64_t, 4> limbs = {};
        for (uint64_t& limb : limbs) {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            limb = state;
        }
        const auto scalar = scalar_bytes(limbs);
        const std::string case_name =
            "random=" + std::to_string(iteration);
        expect_single_group_multiplication_matches_oracle(
            scalar, 256, case_name.c_str());
    }
}

}  // namespace
