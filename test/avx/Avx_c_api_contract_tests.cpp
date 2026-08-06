/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "bindings/blst.h"
#include "blst/avx/curve/point.h"
#include "blst/scalar/curve/map_to_G.h"

namespace {

using G1 = blst::scalar::G1Point;
using G1Affine = blst::scalar::G1AffinePoint;
using G2 = blst::scalar::G2Point;
using G2Affine = blst::scalar::G2AffinePoint;
using AvxG1 = blst::avx::G1Point;
using AvxG1Affine = blst::avx::G1AffinePoint;
using AvxG2 = blst::avx::G2Point;
using AvxG2Affine = blst::avx::G2AffinePoint;

AvxG1 avx_from_scalar(const G1Affine& point)
{
    AvxG1 out = {};
    AvxG1::from_affine(out, AvxG1Affine(
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.x),
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.y)));
    return out;
}

AvxG2 avx_from_scalar(const G2Affine& point)
{
    AvxG2 out = {};
    AvxG2::from_affine(out, AvxG2Affine(
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.x),
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.y)));
    return out;
}

G1 scalar_from_avx(const AvxG1& point)
{
    return blst::avx::avx_curve_detail::curve_g1_to_scalar(
        {point.X, point.Y, point.Z});
}

G2 scalar_from_avx(const AvxG2& point)
{
    return blst::avx::avx_curve_detail::curve_g2_to_scalar(
        {point.X, point.Y, point.Z});
}

std::array<blst_byte, 32> maximal_scalar()
{
    std::array<blst_byte, 32> scalar = {};
    scalar.fill(0xff);
    return scalar;
}

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

TEST(AvxCApiContract, G1SubgroupMsmHasOneContractAcrossPointCounts)
{
    ASSERT_TRUE(blst_avx_enabled());
    const auto scalar = maximal_scalar();

    G1 expected = {};
    blst_p1_mult(&expected, blst_p1_generator(), scalar.data(), 256);

    G1Affine generator_affine = {};
    blst_p1_to_affine(&generator_affine, blst_p1_generator());
    blst_avx_p1_affine avx_generator_affine = {};
    ASSERT_TRUE(blst_avx_p1_affine_from_scalar(
        &avx_generator_affine, &generator_affine));

    const size_t affine_scratch_bytes =
        blst_avx_p1_mult_pippenger_scratch_sizeof(1);
    const size_t affine_flag_bytes =
        blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof(1);
    ASSERT_GT(affine_scratch_bytes, 0u);
    ASSERT_GT(affine_flag_bytes, 0u);
    std::vector<blst_avx_p1> affine_scratch(
        (affine_scratch_bytes + sizeof(blst_avx_p1) - 1) /
        sizeof(blst_avx_p1));
    std::vector<blst_byte> affine_flags(affine_flag_bytes);

    blst_avx_p1 avx_result = {};
    ASSERT_TRUE(blst_avx_p1_mult_pippenger(
        &avx_result, &avx_generator_affine, scalar.data(), 1, 256,
        affine_scratch.data(), affine_scratch.size() * sizeof(blst_avx_p1),
        affine_flags.data(), affine_flags.size()));
    G1 result = {};
    ASSERT_TRUE(blst_p1_from_avx(&result, &avx_result));
    EXPECT_TRUE(G1::equal(result, expected));

    ASSERT_TRUE(blst_avx_p1s_mult_pippenger(
        &avx_result, blst_avx_p1_generator(), 1,
        scalar.data(), 256, nullptr));
    ASSERT_TRUE(blst_p1_from_avx(&result, &avx_result));
    EXPECT_TRUE(G1::equal(result, expected));

    const std::array<blst_avx_p1, 2> points = {
        *blst_avx_p1_generator(), *blst_avx_p1_generator()};
    std::array<blst_byte, 64> scalars = {};
    scalars.fill(0xff);
    G1 expected_pair = {};
    blst_p1_double(&expected_pair, &expected);
    const size_t scratch_bytes =
        blst_avx_p1s_mult_pippenger_scratch_sizeof(points.size());
    std::vector<blst_avx_p1> scratch(
        (scratch_bytes + sizeof(blst_avx_p1) - 1) / sizeof(blst_avx_p1));
    ASSERT_TRUE(blst_avx_p1s_mult_pippenger(
        &avx_result, points.data(), points.size(), scalars.data(), 256,
        reinterpret_cast<blst_limb_t*>(scratch.data())));
    ASSERT_TRUE(blst_p1_from_avx(&result, &avx_result));
    EXPECT_TRUE(G1::equal(result, expected_pair));
}

TEST(AvxCApiContract, G2SubgroupMsmHasOneContractAcrossPointCounts)
{
    ASSERT_TRUE(blst_avx_enabled());
    const auto scalar = maximal_scalar();

    G2 expected = {};
    blst_p2_mult(&expected, blst_p2_generator(), scalar.data(), 256);
    G2Affine generator_affine = {};
    blst_p2_to_affine(&generator_affine, blst_p2_generator());
    blst_avx_p2_affine avx_generator_affine = {};
    ASSERT_TRUE(blst_avx_p2_affine_from_scalar(
        &avx_generator_affine, &generator_affine));

    const size_t affine_scratch_bytes =
        blst_avx_p2_mult_pippenger_scratch_sizeof(1);
    const size_t affine_flag_bytes =
        blst_avx_p2_mult_pippenger_bucket_used_scratch_sizeof(1);
    ASSERT_GT(affine_scratch_bytes, 0u);
    ASSERT_GT(affine_flag_bytes, 0u);
    std::vector<blst_avx_p2> affine_scratch(
        (affine_scratch_bytes + sizeof(blst_avx_p2) - 1) /
        sizeof(blst_avx_p2));
    std::vector<blst_byte> affine_flags(affine_flag_bytes);

    blst_avx_p2 avx_result = {};
    ASSERT_TRUE(blst_avx_p2_mult_pippenger(
        &avx_result, &avx_generator_affine, scalar.data(), 1, 256,
        affine_scratch.data(), affine_scratch.size() * sizeof(blst_avx_p2),
        affine_flags.data(), affine_flags.size()));
    G2 result = {};
    ASSERT_TRUE(blst_p2_from_avx(&result, &avx_result));
    EXPECT_TRUE(G2::equal(result, expected));
    ASSERT_TRUE(blst_avx_p2s_mult_pippenger(
        &avx_result, blst_avx_p2_generator(), 1,
        scalar.data(), 256, nullptr));
    ASSERT_TRUE(blst_p2_from_avx(&result, &avx_result));
    EXPECT_TRUE(G2::equal(result, expected));

    const std::array<blst_avx_p2, 2> points = {
        *blst_avx_p2_generator(), *blst_avx_p2_generator()};
    std::array<blst_byte, 64> scalars = {};
    scalars.fill(0xff);
    G2 expected_pair = {};
    blst_p2_double(&expected_pair, &expected);
    const size_t scratch_bytes =
        blst_avx_p2s_mult_pippenger_scratch_sizeof(points.size());
    std::vector<blst_avx_p2> scratch(
        (scratch_bytes + sizeof(blst_avx_p2) - 1) / sizeof(blst_avx_p2));
    ASSERT_TRUE(blst_avx_p2s_mult_pippenger(
        &avx_result, points.data(), points.size(), scalars.data(), 256,
        reinterpret_cast<blst_limb_t*>(scratch.data())));
    ASSERT_TRUE(blst_p2_from_avx(&result, &avx_result));
    EXPECT_TRUE(G2::equal(result, expected_pair));
}

TEST(AvxCApiContract, RejectsOversizedScalarsBeforeArithmetic)
{
    const auto scalar = maximal_scalar();
    blst_avx_p1 p1 = {};
    blst_avx_p2 p2 = {};

    EXPECT_FALSE(blst_avx_p1_mult(&p1, &p1, scalar.data(), 257));
    EXPECT_FALSE(blst_avx_p2_mult(&p2, &p2, scalar.data(), 257));
    EXPECT_FALSE(blst_avx_p1s_mult_pippenger(
        &p1, &p1, 1, scalar.data(), 257, nullptr));
    EXPECT_FALSE(blst_avx_p2s_mult_pippenger(
        &p2, &p2, 1, scalar.data(), 257, nullptr));
}

TEST(AvxCApiContract, MembershipAndBatchChecksRejectOffCurveInputs)
{
    ASSERT_TRUE(blst_avx_enabled());

    G1 off_curve_g1 = *blst_p1_generator();
    off_curve_g1.x[0] ^= 1;
    blst_avx_p1 avx_off_curve_g1 = {};
    ASSERT_TRUE(blst_avx_p1_from_scalar(&avx_off_curve_g1, &off_curve_g1));
    EXPECT_FALSE(blst_avx_p1_in_g1(&avx_off_curve_g1));

    G2 off_curve_g2 = *blst_p2_generator();
    off_curve_g2.x[0][0] ^= 1;
    blst_avx_p2 avx_off_curve_g2 = {};
    ASSERT_TRUE(blst_avx_p2_from_scalar(&avx_off_curve_g2, &off_curve_g2));
    EXPECT_FALSE(blst_avx_p2_in_g2(&avx_off_curve_g2));

    const blst_byte coefficient = 0;
    const blst_avx_p1* g1_points[] = {&avx_off_curve_g1};
    const blst_avx_p2* g2_points[] = {&avx_off_curve_g2};
    EXPECT_FALSE(blst_avx_p1_batch_subgroup_check(
        g1_points, 1, &coefficient, 1, blst_avx_p1_generator()));
    EXPECT_FALSE(blst_avx_p2_batch_subgroup_check(
        g2_points, 1, &coefficient, 1, blst_avx_p2_generator()));
    EXPECT_FALSE(blst_avx_p1_batch_subgroup_check(
        nullptr, 0, nullptr, 0, &avx_off_curve_g1));
    EXPECT_FALSE(blst_avx_p2_batch_subgroup_check(
        nullptr, 0, nullptr, 0, &avx_off_curve_g2));
}

TEST(AvxCApiContract, BasicOperationsSupportOnCurveNonSubgroupPoints)
{
    ASSERT_TRUE(blst_avx_enabled());

    G1Affine g1_affine = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_p1_uncompress(&g1_affine, kG1Order11Encoded.data()));
    ASSERT_TRUE(G1Affine::on_curve(g1_affine));
    ASSERT_FALSE(G1Affine::in_subgroup(g1_affine));
    G1 g1 = {};
    G1::from_affine(g1, g1_affine);

    blst_avx_p1 avx_g1 = {};
    ASSERT_TRUE(blst_avx_p1_from_scalar(&avx_g1, &g1));
    EXPECT_FALSE(blst_avx_p1_in_g1(&avx_g1));
    blst_avx_p1 avx_g1_twice = {};
    ASSERT_TRUE(blst_avx_p1_add(&avx_g1_twice, &avx_g1, &avx_g1));
    G1 g1_twice = {};
    ASSERT_TRUE(blst_p1_from_avx(&g1_twice, &avx_g1_twice));
    G1 g1_expected = {};
    G1::doubling(g1_expected, g1);
    EXPECT_TRUE(G1::equal(g1_twice, g1_expected));
    EXPECT_TRUE(G1::on_curve(g1_twice));
    EXPECT_FALSE(G1::in_subgroup(g1_twice));

    ASSERT_TRUE(blst_avx_p1_cneg(&avx_g1, true));
    ASSERT_TRUE(blst_avx_p1_cneg(&avx_g1, true));
    G1 g1_after_negation = {};
    ASSERT_TRUE(blst_p1_from_avx(&g1_after_negation, &avx_g1));
    EXPECT_TRUE(G1::equal(g1_after_negation, g1));

    blst_avx_p1_affine avx_g1_affine = {};
    ASSERT_TRUE(blst_avx_p1_to_affine(&avx_g1_affine, &avx_g1));
    blst_avx_p1 avx_g1_from_affine = {};
    ASSERT_TRUE(blst_avx_p1_from_affine(
        &avx_g1_from_affine, &avx_g1_affine));
    EXPECT_TRUE(blst_avx_p1_is_equal(&avx_g1_from_affine, &avx_g1));
    EXPECT_FALSE(blst_avx_p1_is_inf(&avx_g1_from_affine));

    std::array<blst_byte, 96> g1_serialized = {};
    ASSERT_TRUE(blst_avx_p1_serialize(g1_serialized.data(), &avx_g1));
    blst_avx_p1 avx_g1_deserialized = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_avx_p1_deserialize(
            &avx_g1_deserialized, g1_serialized.data()));
    EXPECT_TRUE(blst_avx_p1_is_equal(&avx_g1_deserialized, &avx_g1));
    EXPECT_FALSE(blst_avx_p1_in_g1(&avx_g1_deserialized));

    std::array<blst_byte, 48> g1_compressed = {};
    ASSERT_TRUE(blst_avx_p1_compress(g1_compressed.data(), &avx_g1));
    blst_avx_p1 avx_g1_uncompressed = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_avx_p1_uncompress(
            &avx_g1_uncompressed, g1_compressed.data()));
    EXPECT_TRUE(blst_avx_p1_is_equal(&avx_g1_uncompressed, &avx_g1));
    EXPECT_FALSE(blst_avx_p1_in_g1(&avx_g1_uncompressed));

    const blst_byte coefficient = 1;
    const blst_avx_p1* g1_points[] = {&avx_g1};
    EXPECT_FALSE(blst_avx_p1_batch_subgroup_check(
        g1_points, 1, &coefficient, 1, blst_avx_p1_generator()));

    G2Affine g2_affine = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_p2_uncompress(&g2_affine, kG2Order13Encoded.data()));
    ASSERT_TRUE(G2Affine::on_curve(g2_affine));
    ASSERT_FALSE(G2Affine::in_subgroup(g2_affine));
    G2 g2 = {};
    G2::from_affine(g2, g2_affine);

    blst_avx_p2 avx_g2 = {};
    ASSERT_TRUE(blst_avx_p2_from_scalar(&avx_g2, &g2));
    EXPECT_FALSE(blst_avx_p2_in_g2(&avx_g2));
    blst_avx_p2 avx_g2_twice = {};
    ASSERT_TRUE(blst_avx_p2_add(&avx_g2_twice, &avx_g2, &avx_g2));
    G2 g2_twice = {};
    ASSERT_TRUE(blst_p2_from_avx(&g2_twice, &avx_g2_twice));
    G2 g2_expected = {};
    G2::doubling(g2_expected, g2);
    EXPECT_TRUE(G2::equal(g2_twice, g2_expected));
    EXPECT_TRUE(G2::on_curve(g2_twice));
    EXPECT_FALSE(G2::in_subgroup(g2_twice));

    ASSERT_TRUE(blst_avx_p2_cneg(&avx_g2, true));
    ASSERT_TRUE(blst_avx_p2_cneg(&avx_g2, true));
    G2 g2_after_negation = {};
    ASSERT_TRUE(blst_p2_from_avx(&g2_after_negation, &avx_g2));
    EXPECT_TRUE(G2::equal(g2_after_negation, g2));

    blst_avx_p2_affine avx_g2_affine = {};
    ASSERT_TRUE(blst_avx_p2_to_affine(&avx_g2_affine, &avx_g2));
    blst_avx_p2 avx_g2_from_affine = {};
    ASSERT_TRUE(blst_avx_p2_from_affine(
        &avx_g2_from_affine, &avx_g2_affine));
    EXPECT_TRUE(blst_avx_p2_is_equal(&avx_g2_from_affine, &avx_g2));
    EXPECT_FALSE(blst_avx_p2_is_inf(&avx_g2_from_affine));

    std::array<blst_byte, 192> g2_serialized = {};
    ASSERT_TRUE(blst_avx_p2_serialize(g2_serialized.data(), &avx_g2));
    blst_avx_p2 avx_g2_deserialized = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_avx_p2_deserialize(
            &avx_g2_deserialized, g2_serialized.data()));
    EXPECT_TRUE(blst_avx_p2_is_equal(&avx_g2_deserialized, &avx_g2));
    EXPECT_FALSE(blst_avx_p2_in_g2(&avx_g2_deserialized));

    std::array<blst_byte, 96> g2_compressed = {};
    ASSERT_TRUE(blst_avx_p2_compress(g2_compressed.data(), &avx_g2));
    blst_avx_p2 avx_g2_uncompressed = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_avx_p2_uncompress(
            &avx_g2_uncompressed, g2_compressed.data()));
    EXPECT_TRUE(blst_avx_p2_is_equal(&avx_g2_uncompressed, &avx_g2));
    EXPECT_FALSE(blst_avx_p2_in_g2(&avx_g2_uncompressed));

    const blst_avx_p2* g2_points[] = {&avx_g2};
    EXPECT_FALSE(blst_avx_p2_batch_subgroup_check(
        g2_points, 1, &coefficient, 1, blst_avx_p2_generator()));
}

TEST(AvxNativePointContract, UncheckedMultiplicationIsExactAndSupportsLowOrderPoints)
{
    ASSERT_TRUE(blst_avx_enabled());

    G1Affine g1_affine = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_p1_uncompress(&g1_affine, kG1Order11Encoded.data()));
    G1 g1 = {};
    G1::from_affine(g1, g1_affine);
    const AvxG1 avx_g1 = avx_from_scalar(g1_affine);

    std::array<blst_byte, 32> scalar = {};
    scalar[0] = 11;
    AvxG1 avx_g1_result = {};
    AvxG1::unchecked_mult(
        avx_g1_result, avx_g1, scalar.data(), 256);
    EXPECT_TRUE(AvxG1::is_inf(avx_g1_result));

    scalar[0] = 12;
    AvxG1::unchecked_mult(
        avx_g1_result, avx_g1, scalar.data(), 256);
    EXPECT_TRUE(AvxG1::equal(avx_g1_result, avx_g1));

    G2Affine g2_affine = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_p2_uncompress(&g2_affine, kG2Order13Encoded.data()));
    G2 g2 = {};
    G2::from_affine(g2, g2_affine);
    const AvxG2 avx_g2 = avx_from_scalar(g2_affine);

    scalar.fill(0);
    scalar[0] = 13;
    AvxG2 avx_g2_result = {};
    AvxG2::unchecked_mult(
        avx_g2_result, avx_g2, scalar.data(), 256);
    EXPECT_TRUE(AvxG2::is_inf(avx_g2_result));

    scalar[0] = 14;
    AvxG2::unchecked_mult(
        avx_g2_result, avx_g2, scalar.data(), 256);
    EXPECT_TRUE(AvxG2::equal(avx_g2_result, avx_g2));

    scalar.fill(0xff);
    G1 scalar_g1_result = {};
    G1::unchecked_mult(scalar_g1_result, g1, scalar.data(), 256);
    AvxG1 avx_g1_alias = avx_g1;
    AvxG1::unchecked_mult(
        avx_g1_alias, avx_g1_alias, scalar.data(), 256);
    const G1 converted_g1_result = scalar_from_avx(avx_g1_alias);
    EXPECT_TRUE(G1::equal(converted_g1_result, scalar_g1_result));

    G2 scalar_g2_result = {};
    G2::unchecked_mult(scalar_g2_result, g2, scalar.data(), 256);
    AvxG2 avx_g2_alias = avx_g2;
    AvxG2::unchecked_mult(
        avx_g2_alias, avx_g2_alias, scalar.data(), 256);
    const G2 converted_g2_result = scalar_from_avx(avx_g2_alias);
    EXPECT_TRUE(G2::equal(converted_g2_result, scalar_g2_result));

    for (const size_t bits : {size_t{1}, size_t{128}, size_t{175}, size_t{176}}) {
        for (size_t i = 0; i < scalar.size(); ++i)
            scalar[i] = static_cast<blst_byte>(0xa5u ^ (17u * i));

        G1 expected_g1 = {};
        G1::unchecked_mult(expected_g1, g1, scalar.data(), bits);
        AvxG1 actual_g1 = {};
        AvxG1::unchecked_mult(actual_g1, avx_g1, scalar.data(), bits);
        EXPECT_TRUE(G1::equal(scalar_from_avx(actual_g1), expected_g1))
            << "G1 nbits=" << bits;

        G2 expected_g2 = {};
        G2::unchecked_mult(expected_g2, g2, scalar.data(), bits);
        AvxG2 actual_g2 = {};
        AvxG2::unchecked_mult(actual_g2, avx_g2, scalar.data(), bits);
        EXPECT_TRUE(G2::equal(scalar_from_avx(actual_g2), expected_g2))
            << "G2 nbits=" << bits;
    }

    AvxG1 zero_g1 = avx_g1;
    AvxG1::unchecked_mult(zero_g1, zero_g1, nullptr, 0);
    EXPECT_TRUE(AvxG1::is_inf(zero_g1));
    AvxG2 zero_g2 = avx_g2;
    AvxG2::unchecked_mult(zero_g2, zero_g2, nullptr, 0);
    EXPECT_TRUE(AvxG2::is_inf(zero_g2));
}

TEST(AvxNativePointContract, RawUncompressFailureDoesNotPublishOutput)
{
    std::array<blst_byte, 48> invalid_g1 = {};
    invalid_g1.back() = 1;
    AvxG1Affine g1 = AvxG1Affine::generator();
    const AvxG1Affine expected_g1 = g1;
    EXPECT_GT(AvxG1Affine::uncompress_BE(
        g1, invalid_g1.data()), 1u);
    EXPECT_TRUE(AvxG1Affine::equal(g1, expected_g1));

    std::array<blst_byte, 96> invalid_g2 = {};
    invalid_g2.back() = 1;
    AvxG2Affine g2 = AvxG2Affine::generator();
    const AvxG2Affine expected_g2 = g2;
    EXPECT_GT(AvxG2Affine::uncompress_BE(
        g2, invalid_g2.data()), 1u);
    EXPECT_TRUE(AvxG2Affine::equal(g2, expected_g2));
}

TEST(AvxCApiContract, GTNameMatchesScalarBackend)
{
    ASSERT_TRUE(blst_avx_enabled());
    const auto scalar = maximal_scalar();

    blst_avx_gt point = {};
    ASSERT_TRUE(blst_avx_fp12_from_scalar(&point, blst_gt_generator()));
    const size_t scratch_bytes =
        blst_avx_gt_mult_pippenger_scratch_sizeof(1);
    const size_t flag_bytes =
        blst_avx_gt_mult_pippenger_bucket_used_scratch_sizeof(1);
    ASSERT_GT(scratch_bytes, 0u);
    ASSERT_GT(flag_bytes, 0u);
    std::vector<blst_avx_gt> scratch(
        (scratch_bytes + sizeof(blst_avx_gt) - 1) / sizeof(blst_avx_gt));
    std::vector<blst_byte> flags(flag_bytes);

    blst_avx_gt actual_avx = {};
    ASSERT_TRUE(blst_avx_gt_mult_pippenger(
        &actual_avx, &point, scalar.data(), 1,
        scratch.data(), scratch.size() * sizeof(blst_avx_gt),
        flags.data(), flags.size()));

    blst_fp12 expected = {};
    blst_gt_mult(&expected, blst_gt_generator(), scalar.data(), 256);
    blst_fp12 actual = {};
    ASSERT_TRUE(blst_fp12_from_avx(&actual, &actual_avx));
    EXPECT_TRUE(blst_gt_is_equal(&actual, &expected));
}

} // namespace
