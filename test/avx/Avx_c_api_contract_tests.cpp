/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstring>
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

G1Affine scalar_from_avx(const AvxG1Affine& point)
{
    return {
        blst::avx::avx_curve_detail::ring_fp_to_scalar(point.X),
        blst::avx::avx_curve_detail::ring_fp_to_scalar(point.Y),
    };
}

G2Affine scalar_from_avx(const AvxG2Affine& point)
{
    return {
        blst::avx::avx_curve_detail::ring_fp2_to_scalar(point.X),
        blst::avx::avx_curve_detail::ring_fp2_to_scalar(point.Y),
    };
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

struct FixedRandomSource {
    std::array<blst_byte, 32> seed = {};
    size_t calls = 0;
    bool succeeds = true;
};

bool fill_fixed_random(void* context, blst_byte* out, size_t out_len)
{
    auto& source = *static_cast<FixedRandomSource*>(context);
    ++source.calls;
    if (!source.succeeds || out == nullptr || out_len != source.seed.size()) {
        return false;
    }
    std::memcpy(out, source.seed.data(), source.seed.size());
    return true;
}

#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
bool throw_from_random_callback(void*, blst_byte*, size_t)
{
    throw 1;
}
#endif

TEST(AvxCApiContract, FriendlyBatchMembershipOwnsChallengeConstruction)
{
    ASSERT_TRUE(blst_avx_enabled());

    EXPECT_TRUE(blst_avx_p1s_in_g1(nullptr, 0, nullptr, nullptr));
    EXPECT_TRUE(blst_avx_p2s_in_g2(
        nullptr, 0, nullptr, nullptr, nullptr, 0));

    const size_t g2_workspace_bytes =
        blst_avx_p2s_in_g2_workspace_sizeof();
    ASSERT_GT(g2_workspace_bytes, 0u);
    std::vector<blst_byte> g2_workspace(g2_workspace_bytes);

    FixedRandomSource random;
    random.seed[0] = 0xa5;
    random.seed[31] = 0x5a;

    const blst_avx_p1* g1_points[] = {blst_avx_p1_generator()};
    const blst_avx_p2* g2_points[] = {blst_avx_p2_generator()};
    EXPECT_FALSE(blst_avx_p2s_in_g2(
        g2_points, 1, fill_fixed_random, &random, nullptr, 0));
    EXPECT_FALSE(blst_avx_p2s_in_g2(
        g2_points, 1, fill_fixed_random, &random,
        g2_workspace.data(), g2_workspace.size() - 1));
    EXPECT_EQ(random.calls, 0u);
    EXPECT_FALSE(blst_avx_p1s_in_g1(nullptr, 1, fill_fixed_random, &random));
    EXPECT_FALSE(blst_avx_p2s_in_g2(
        nullptr, 1, fill_fixed_random, &random,
        g2_workspace.data(), g2_workspace.size()));
    EXPECT_EQ(random.calls, 0u);

    random.succeeds = false;
    EXPECT_FALSE(blst_avx_p1s_in_g1(
        g1_points, 1, fill_fixed_random, &random));
    EXPECT_FALSE(blst_avx_p2s_in_g2(
        g2_points, 1, fill_fixed_random, &random,
        g2_workspace.data(), g2_workspace.size()));
    EXPECT_EQ(random.calls, 2u);

#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
    EXPECT_FALSE(blst_avx_p1s_in_g1(
        g1_points, 1, throw_from_random_callback, nullptr));
    EXPECT_FALSE(blst_avx_p2s_in_g2(
        g2_points, 1, throw_from_random_callback, nullptr,
        g2_workspace.data(), g2_workspace.size()));
#endif

    random.succeeds = true;
    EXPECT_TRUE(blst_avx_p1s_in_g1(
        g1_points, 1, fill_fixed_random, &random));
    EXPECT_TRUE(blst_avx_p2s_in_g2(
        g2_points, 1, fill_fixed_random, &random,
        g2_workspace.data(), g2_workspace.size()));
    EXPECT_EQ(random.calls, 4u);
}

TEST(AvxCApiContract, FriendlyBatchMembershipRejectsLowOrderPoints)
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

    G2Affine g2_affine = {};
    ASSERT_EQ(BLST_SUCCESS,
        blst_p2_uncompress(&g2_affine, kG2Order13Encoded.data()));
    ASSERT_TRUE(G2Affine::on_curve(g2_affine));
    ASSERT_FALSE(G2Affine::in_subgroup(g2_affine));
    G2 g2 = {};
    G2::from_affine(g2, g2_affine);
    blst_avx_p2 avx_g2 = {};
    ASSERT_TRUE(blst_avx_p2_from_scalar(&avx_g2, &g2));

    FixedRandomSource random;
    random.seed.fill(0);
    const blst_avx_p1* g1_points[] = {&avx_g1};
    const blst_avx_p2* g2_points[] = {&avx_g2};
    const size_t g2_workspace_bytes =
        blst_avx_p2s_in_g2_workspace_sizeof();
    ASSERT_GT(g2_workspace_bytes, 0u);
    std::vector<blst_byte> g2_workspace(g2_workspace_bytes);
    EXPECT_FALSE(blst_avx_p1s_in_g1(
        g1_points, 1, fill_fixed_random, &random));
    EXPECT_FALSE(blst_avx_p2s_in_g2(
        g2_points, 1, fill_fixed_random, &random,
        g2_workspace.data(), g2_workspace.size()));
    EXPECT_EQ(random.calls, 2u);

    const blst_avx_p1* mixed_g1_points[] = {
        blst_avx_p1_generator(), &avx_g1, blst_avx_p1_generator()};
    const blst_avx_p2* mixed_g2_points[] = {
        blst_avx_p2_generator(), &avx_g2, blst_avx_p2_generator()};
    EXPECT_FALSE(blst_avx_p1s_in_g1(
        mixed_g1_points, 3, fill_fixed_random, &random));
    EXPECT_FALSE(blst_avx_p2s_in_g2(
        mixed_g2_points, 3, fill_fixed_random, &random,
        g2_workspace.data(), g2_workspace.size()));
    EXPECT_EQ(random.calls, 4u);
}

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

    blst_avx_p1 aliased_one = *blst_avx_p1_generator();
    ASSERT_TRUE(blst_avx_p1s_mult_pippenger(
        &aliased_one, &aliased_one, 1, scalar.data(), 256, nullptr));
    ASSERT_TRUE(blst_p1_from_avx(&result, &aliased_one));
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
    blst_avx_p2 aliased_one = *blst_avx_p2_generator();
    ASSERT_TRUE(blst_avx_p2s_mult_pippenger(
        &aliased_one, &aliased_one, 1, scalar.data(), 256, nullptr));
    ASSERT_TRUE(blst_p2_from_avx(&result, &aliased_one));
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

TEST(AvxCApiContract, RejectsOverlappingMsmLiveRanges)
{
    ASSERT_TRUE(blst_avx_enabled());
    constexpr size_t kCount = 2;
    std::array<blst_byte, 32 * kCount> scalars = {};
    scalars[0] = 1;
    scalars[32] = 2;

    std::array<blst_avx_p1_affine, kCount> g1_points = {};
    std::array<blst_avx_p2_affine, kCount> g2_points = {};
    for (size_t i = 0; i < kCount; ++i) {
        ASSERT_TRUE(blst_avx_p1_to_affine(
            &g1_points[i], blst_avx_p1_generator()));
        ASSERT_TRUE(blst_avx_p2_to_affine(
            &g2_points[i], blst_avx_p2_generator()));
    }

    const size_t g1_scratch_bytes =
        blst_avx_p1_mult_pippenger_scratch_sizeof(kCount);
    const size_t g1_flags_bytes =
        blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof(kCount);
    std::vector<blst_avx_p1> g1_scratch(
        (g1_scratch_bytes + g1_flags_bytes + sizeof(blst_avx_p1) - 1) /
        sizeof(blst_avx_p1));
    std::vector<blst_byte> g1_flags(g1_flags_bytes);
    blst_avx_p1 g1_out = *blst_avx_p1_generator();
    const blst_avx_p1 g1_before = g1_out;
    EXPECT_FALSE(blst_avx_p1_mult_pippenger(
        &g1_out, g1_points.data(), scalars.data(), kCount, 256,
        g1_scratch.data(), g1_scratch_bytes,
        scalars.data(), g1_flags_bytes));
    EXPECT_EQ(std::memcmp(&g1_out, &g1_before, sizeof(g1_out)), 0);

    EXPECT_FALSE(blst_avx_p1_mult_pippenger(
        &g1_out, g1_points.data(), scalars.data(), kCount, 256,
        g1_scratch.data(), g1_scratch_bytes,
        reinterpret_cast<blst_byte*>(g1_scratch.data()) + 1,
        g1_flags_bytes));
    EXPECT_EQ(std::memcmp(&g1_out, &g1_before, sizeof(g1_out)), 0);

    const size_t g2_scratch_bytes =
        blst_avx_p2_mult_pippenger_scratch_sizeof(kCount);
    const size_t g2_flags_bytes =
        blst_avx_p2_mult_pippenger_bucket_used_scratch_sizeof(kCount);
    std::vector<blst_avx_p2> g2_scratch(
        (g2_scratch_bytes + sizeof(blst_avx_p2) - 1) / sizeof(blst_avx_p2));
    std::vector<blst_byte> g2_flags(g2_flags_bytes);
    blst_avx_p2 g2_out = *blst_avx_p2_generator();
    const blst_avx_p2 g2_before = g2_out;
    EXPECT_FALSE(blst_avx_p2_mult_pippenger(
        &g2_out, g2_points.data(), scalars.data(), kCount, 256,
        g2_scratch.data(), g2_scratch_bytes,
        scalars.data(), g2_flags_bytes));
    EXPECT_EQ(std::memcmp(&g2_out, &g2_before, sizeof(g2_out)), 0);

    std::array<blst_avx_gt, kCount> gt_points = {};
    for (auto& point : gt_points) {
        ASSERT_TRUE(blst_avx_fp12_from_scalar(&point, blst_gt_generator()));
    }
    const size_t gt_scratch_bytes =
        blst_avx_gt_mult_pippenger_scratch_sizeof(kCount);
    const size_t gt_flags_bytes =
        blst_avx_gt_mult_pippenger_bucket_used_scratch_sizeof(kCount);
    std::vector<blst_avx_gt> gt_scratch(
        (gt_scratch_bytes + sizeof(blst_avx_gt) - 1) / sizeof(blst_avx_gt));
    std::vector<blst_byte> gt_flags(gt_flags_bytes);
    blst_avx_gt gt_out = {};
    ASSERT_TRUE(blst_avx_fp12_from_scalar(&gt_out, blst_gt_generator()));
    const blst_avx_gt gt_before = gt_out;
    EXPECT_FALSE(blst_avx_gt_mult_pippenger(
        &gt_out, gt_points.data(), scalars.data(), kCount,
        gt_scratch.data(), gt_scratch_bytes,
        scalars.data(), gt_flags_bytes));
    EXPECT_EQ(std::memcmp(&gt_out, &gt_before, sizeof(gt_out)), 0);

    const size_t split_point_bytes =
        blst_avx_p1_mult_pippenger_glv_points_scratch_sizeof(kCount);
    const size_t split_scalar_bytes =
        blst_avx_p1_mult_pippenger_glv_scalars_scratch_sizeof(kCount);
    std::vector<blst_avx_p1_affine> split_points(
        split_point_bytes / sizeof(blst_avx_p1_affine));
    std::vector<blst_byte> split_scalars(split_scalar_bytes);
    EXPECT_FALSE(blst_avx_p1_mult_pippenger_glv_with_scratch(
        &g1_out, g1_points.data(), scalars.data(), kCount,
        g1_scratch.data(), g1_scratch_bytes,
        g1_flags.data(), g1_flags.size(),
        split_points.data(), split_point_bytes,
        scalars.data(), split_scalar_bytes));
    EXPECT_EQ(std::memcmp(&g1_out, &g1_before, sizeof(g1_out)), 0);

    const size_t window_scratch_bytes =
        blst_avx_p1_mult_pippenger_glv_window_scratch_sizeof(kCount, 4);
    const size_t window_flags_bytes =
        blst_avx_p1_mult_pippenger_glv_window_bucket_used_scratch_sizeof(kCount, 4);
    std::vector<blst_avx_p1> window_scratch(
        (window_scratch_bytes + sizeof(blst_avx_p1) - 1) /
        sizeof(blst_avx_p1));
    std::vector<blst_byte> window_flags(window_flags_bytes);
    EXPECT_FALSE(blst_avx_p1_mult_pippenger_glv_window_with_scratch(
        &g1_out, g1_points.data(), scalars.data(), kCount, 4,
        window_scratch.data(), window_scratch_bytes,
        window_flags.data(), window_flags.size(),
        split_points.data(), split_point_bytes,
        scalars.data(), split_scalar_bytes));
    EXPECT_EQ(std::memcmp(&g1_out, &g1_before, sizeof(g1_out)), 0);

    const size_t projective_scratch_bytes =
        blst_avx_p1s_mult_pippenger_scratch_sizeof(kCount);
    std::vector<blst_avx_p1> projective_points(
        std::max(kCount,
            (projective_scratch_bytes + sizeof(blst_avx_p1) - 1) /
                sizeof(blst_avx_p1)),
        *blst_avx_p1_generator());
    EXPECT_FALSE(blst_avx_p1s_mult_pippenger(
        &g1_out, projective_points.data(), kCount, scalars.data(), 256,
        reinterpret_cast<blst_limb_t*>(projective_points.data())));
    EXPECT_EQ(std::memcmp(&g1_out, &g1_before, sizeof(g1_out)), 0);

    const size_t g2_projective_scratch_bytes =
        blst_avx_p2s_mult_pippenger_scratch_sizeof(kCount);
    std::vector<blst_avx_p2> g2_projective_points(
        std::max(kCount,
            (g2_projective_scratch_bytes + sizeof(blst_avx_p2) - 1) /
                sizeof(blst_avx_p2)),
        *blst_avx_p2_generator());
    EXPECT_FALSE(blst_avx_p2s_mult_pippenger(
        &g2_out, g2_projective_points.data(), kCount, scalars.data(), 256,
        reinterpret_cast<blst_limb_t*>(g2_projective_points.data())));
    EXPECT_EQ(std::memcmp(&g2_out, &g2_before, sizeof(g2_out)), 0);
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
    std::vector<blst_byte> g2_workspace(
        blst_avx_p2s_in_g2_workspace_sizeof());
    FixedRandomSource random;
    EXPECT_FALSE(blst_avx_p1s_in_g1(
        g1_points, 1, fill_fixed_random, &random));
    EXPECT_FALSE(blst_avx_p2s_in_g2(
        g2_points, 1, fill_fixed_random, &random,
        g2_workspace.data(), g2_workspace.size()));
    EXPECT_EQ(random.calls, 0u);
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
    EXPECT_TRUE(G1::is_inf(scalar_from_avx(avx_g1_result)));
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
    EXPECT_TRUE(G2::is_inf(scalar_from_avx(avx_g2_result)));
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

TEST(AvxNativePointContract, GeneratedCurveLiteralsMatchScalarBackend)
{
    EXPECT_TRUE(G1::equal(
        scalar_from_avx(AvxG1::generator()), G1::generator()));
    EXPECT_TRUE(G1::equal(
        scalar_from_avx(blst::avx::get_BLST_NEG_G1()),
        blst::scalar::get_BLST_NEG_G1()));
    EXPECT_TRUE(G2::equal(
        scalar_from_avx(AvxG2::generator()), G2::generator()));
    EXPECT_TRUE(G2::equal(
        scalar_from_avx(blst::avx::get_BLST_NEG_G2()),
        blst::scalar::get_BLST_NEG_G2()));

    const blst::Fp aprime_e1 =
        blst::avx::avx_curve_detail::ring_fp_to_scalar(
            blst::avx::get_BLST_Aprime_E1());
    EXPECT_TRUE(blst::Fp::equal(
        aprime_e1, blst::scalar::get_BLST_Aprime_E1()));

    const blst::Fp2 aprime_e2 =
        blst::avx::avx_curve_detail::ring_fp2_to_scalar(
            blst::avx::get_BLST_Aprime_E2());
    const blst::Fp2& scalar_aprime_e2 =
        blst::scalar::get_BLST_Aprime_E2();
    EXPECT_TRUE(blst::Fp::equal(aprime_e2[0], scalar_aprime_e2[0]));
    EXPECT_TRUE(blst::Fp::equal(aprime_e2[1], scalar_aprime_e2[1]));

    blst_fp12 scalar_one = {};
    ASSERT_TRUE(blst_fp12_from_avx(
        &scalar_one, blst_avx_fp12_one()));
    EXPECT_TRUE(blst_fp12_is_equal(&scalar_one, blst_fp12_one()));
}

TEST(AvxNativePointContract, EqualityHandlesEveryInfinityRepresentation)
{
    const auto& ring = blst::avx::avx_curve_detail::curve_ring();
    const auto& fp2ring = blst::avx::avx_curve_detail::curve_fp2_ring();
    const blst::avx::avx_curve_detail::CurveRing::StandardElement modulus{
        blst::avx::avx_curve_detail::CurveRing::StandardResidue(
            blst::avx::AVXVector<8>(
                blst::avx::matrixnok::generated::modulus_in_rns_m1)),
        blst::avx::avx_curve_detail::CurveRing::StandardResidue(
            blst::avx::AVXVector<8>(
                blst::avx::matrixnok::generated::modulus_in_rns_m2))};
    blst::Fp canonical_two = {};
    canonical_two.l[0] = 2;
    const auto scale = ring.from_canonical_limbs(canonical_two);
    const AvxG1 g1 = AvxG1::generator();
    const AvxG1 scaled_g1 = {
        ring.modmul(g1.X, scale),
        ring.modmul(g1.Y, scale),
        ring.modmul(g1.Z, scale),
    };
    AvxG1 doubled_g1 = {};
    AvxG1::doubling(doubled_g1, g1);
    EXPECT_TRUE(AvxG1::equal(g1, scaled_g1));
    EXPECT_FALSE(AvxG1::equal(g1, doubled_g1));

    const AvxG1 canonical_g1_infinity = AvxG1::infinity();
    const AvxG1 zero_g1_infinity = {};
    EXPECT_TRUE(AvxG1::equal(canonical_g1_infinity, zero_g1_infinity));
    EXPECT_FALSE(AvxG1::equal(zero_g1_infinity, g1));
    EXPECT_FALSE(AvxG1::equal(g1, zero_g1_infinity));

    const auto fp2_scale = blst::avx::avx_curve_detail::CurveFp2Ring::StandardElement(
        scale, ring.zero());
    const AvxG2 g2 = AvxG2::generator();
    const AvxG2 scaled_g2 = {
        fp2ring.mul(g2.X, fp2_scale),
        fp2ring.mul(g2.Y, fp2_scale),
        fp2ring.mul(g2.Z, fp2_scale),
    };
    AvxG2 doubled_g2 = {};
    AvxG2::doubling(doubled_g2, g2);
    EXPECT_TRUE(AvxG2::equal(g2, scaled_g2));
    EXPECT_FALSE(AvxG2::equal(g2, doubled_g2));

    const AvxG2 canonical_g2_infinity = AvxG2::infinity();
    const AvxG2 zero_g2_infinity = {};
    EXPECT_TRUE(AvxG2::equal(canonical_g2_infinity, zero_g2_infinity));
    EXPECT_FALSE(AvxG2::equal(zero_g2_infinity, g2));
    EXPECT_FALSE(AvxG2::equal(g2, zero_g2_infinity));

    const AvxG1Affine g1_affine = AvxG1Affine::generator();
    const AvxG1Affine lifted_g1_affine(
        blst::avx::NormalizePointCoord(g1_affine.X + modulus, ring),
        blst::avx::NormalizePointCoord(g1_affine.Y + modulus, ring));
    EXPECT_FALSE(g1_affine.X.equal(lifted_g1_affine.X) &
        g1_affine.Y.equal(lifted_g1_affine.Y));
    AvxG1Affine doubled_g1_affine = {};
    AvxG1::to_affine(doubled_g1_affine, doubled_g1);
    EXPECT_TRUE(AvxG1Affine::equal(g1_affine, g1_affine));
    EXPECT_TRUE(AvxG1Affine::equal(g1_affine, lifted_g1_affine));
    EXPECT_FALSE(AvxG1Affine::equal(g1_affine, doubled_g1_affine));

    const AvxG2Affine g2_affine = AvxG2Affine::generator();
    const auto lifted_fp2 = [&modulus](const auto& value) {
        return blst::avx::avx_curve_detail::CurveFp2Ring::StandardElement(
            blst::avx::NormalizePointCoord(
                value.x() + modulus,
                blst::avx::avx_curve_detail::curve_ring()),
            blst::avx::NormalizePointCoord(
                value.y() + modulus,
                blst::avx::avx_curve_detail::curve_ring()));
    };
    const AvxG2Affine lifted_g2_affine(
        lifted_fp2(g2_affine.X),
        lifted_fp2(g2_affine.Y));
    EXPECT_FALSE(g2_affine.X.equal(lifted_g2_affine.X) &
        g2_affine.Y.equal(lifted_g2_affine.Y));
    AvxG2Affine doubled_g2_affine = {};
    AvxG2::to_affine(doubled_g2_affine, doubled_g2);
    EXPECT_TRUE(AvxG2Affine::equal(g2_affine, g2_affine));
    EXPECT_TRUE(AvxG2Affine::equal(g2_affine, lifted_g2_affine));
    EXPECT_FALSE(AvxG2Affine::equal(g2_affine, doubled_g2_affine));
}

TEST(AvxNativePointContract, SigningAndKeyDerivationStayNativeAndNormalizeOnce)
{
    std::array<blst_byte, 32> key = {};
    key[0] = 3;

    G1 scalar_g1_hash = G1::generator();
    G1::doubling(scalar_g1_hash, scalar_g1_hash);
    G1 scalar_g1_signature = {};
    G1::sign(scalar_g1_signature, scalar_g1_hash, key.data());
    AvxG1 avx_g1_hash = avx_from_scalar(G1Affine::generator());
    AvxG1::doubling(avx_g1_hash, avx_g1_hash);
    AvxG1 avx_g1_signature = {};
    AvxG1::sign(avx_g1_signature, avx_g1_hash, key.data());
    EXPECT_TRUE(G1::equal(scalar_g1_signature, scalar_from_avx(avx_g1_signature)));
    EXPECT_TRUE(avx_g1_signature.Z.equal(
        blst::avx::avx_curve_detail::curve_ring().one()));

    AvxG1 avx_g1_alias = avx_g1_hash;
    AvxG1::sign(avx_g1_alias, avx_g1_alias, key.data());
    EXPECT_TRUE(AvxG1::equal(avx_g1_alias, avx_g1_signature));
    AvxG1 avx_g1_identity_signature = {};
    AvxG1::sign(
        avx_g1_identity_signature, AvxG1::infinity(), key.data());
    EXPECT_TRUE(G1::is_inf(scalar_from_avx(avx_g1_identity_signature)));
    EXPECT_TRUE(AvxG1::is_identity(avx_g1_identity_signature));

    G2 scalar_g2_hash = G2::generator();
    G2::doubling(scalar_g2_hash, scalar_g2_hash);
    G2 scalar_g2_signature = {};
    G2::sign(scalar_g2_signature, scalar_g2_hash, key.data());
    AvxG2 avx_g2_hash = avx_from_scalar(G2Affine::generator());
    AvxG2::doubling(avx_g2_hash, avx_g2_hash);
    AvxG2 avx_g2_signature = {};
    AvxG2::sign(avx_g2_signature, avx_g2_hash, key.data());
    EXPECT_TRUE(G2::equal(scalar_g2_signature, scalar_from_avx(avx_g2_signature)));
    EXPECT_TRUE(avx_g2_signature.Z.equal(
        blst::avx::avx_curve_detail::curve_fp2_ring().one()));

    AvxG2 avx_g2_alias = avx_g2_hash;
    AvxG2::sign(avx_g2_alias, avx_g2_alias, key.data());
    EXPECT_TRUE(AvxG2::equal(avx_g2_alias, avx_g2_signature));
    AvxG2 avx_g2_identity_signature = {};
    AvxG2::sign(
        avx_g2_identity_signature, AvxG2::infinity(), key.data());
    EXPECT_TRUE(G2::is_inf(scalar_from_avx(avx_g2_identity_signature)));
    EXPECT_TRUE(AvxG2::is_identity(avx_g2_identity_signature));

    std::array<blst_byte, 96> scalar_g1_bytes = {};
    std::array<blst_byte, 96> avx_g1_bytes = {};
    G1Affine scalar_g1_affine = {};
    AvxG1Affine avx_g1_affine = {};
    G1::sign2(scalar_g1_bytes.data(), &scalar_g1_affine,
        scalar_g1_hash, key.data());
    AvxG1::sign2(avx_g1_bytes.data(), &avx_g1_affine,
        avx_g1_hash, key.data());
    EXPECT_EQ(avx_g1_bytes, scalar_g1_bytes);
    EXPECT_TRUE(G1Affine::equal(scalar_from_avx(avx_g1_affine), scalar_g1_affine));

    std::array<blst_byte, 192> scalar_g2_bytes = {};
    std::array<blst_byte, 192> avx_g2_bytes = {};
    G2Affine scalar_g2_affine = {};
    AvxG2Affine avx_g2_affine = {};
    G2::sign2(scalar_g2_bytes.data(), &scalar_g2_affine,
        scalar_g2_hash, key.data());
    AvxG2::sign2(avx_g2_bytes.data(), &avx_g2_affine,
        avx_g2_hash, key.data());
    EXPECT_EQ(avx_g2_bytes, scalar_g2_bytes);
    EXPECT_TRUE(G2Affine::equal(scalar_from_avx(avx_g2_affine), scalar_g2_affine));

    blst_scalar c_key = {};
    std::memcpy(c_key.b, key.data(), key.size());
    std::array<blst_byte, 96> c_g1_bytes = {};
    blst_p1_affine c_g1_affine = {};
    blst_sign_pk2_in_g2(
        c_g1_bytes.data(), &c_g1_affine, &scalar_g1_hash, &c_key);
    EXPECT_EQ(c_g1_bytes, scalar_g1_bytes);
    EXPECT_TRUE(G1Affine::equal(c_g1_affine, scalar_g1_affine));
    std::array<blst_byte, 192> c_g2_bytes = {};
    blst_p2_affine c_g2_affine = {};
    blst_sign_pk2_in_g1(
        c_g2_bytes.data(), &c_g2_affine, &scalar_g2_hash, &c_key);
    EXPECT_EQ(c_g2_bytes, scalar_g2_bytes);
    EXPECT_TRUE(G2Affine::equal(c_g2_affine, scalar_g2_affine));

    AvxG1Affine avx_g1_optional = {};
    AvxG1::sign2(nullptr, &avx_g1_optional, avx_g1_hash, key.data());
    EXPECT_TRUE(AvxG1Affine::equal(avx_g1_optional, avx_g1_affine));
    std::array<blst_byte, 96> avx_g1_bytes_only = {};
    AvxG1::sign2(avx_g1_bytes_only.data(), nullptr, avx_g1_hash, key.data());
    EXPECT_EQ(avx_g1_bytes_only, avx_g1_bytes);

    std::array<blst_byte, 96> scalar_g1_pk_bytes = {};
    std::array<blst_byte, 96> avx_g1_pk_bytes = {};
    G1Affine scalar_g1_pk = {};
    AvxG1Affine avx_g1_pk = {};
    G1::to_pk2(scalar_g1_pk_bytes.data(), &scalar_g1_pk, key.data());
    AvxG1::to_pk2(avx_g1_pk_bytes.data(), &avx_g1_pk, key.data());
    EXPECT_EQ(avx_g1_pk_bytes, scalar_g1_pk_bytes);
    EXPECT_TRUE(G1Affine::equal(scalar_from_avx(avx_g1_pk), scalar_g1_pk));

    std::array<blst_byte, 192> scalar_g2_pk_bytes = {};
    std::array<blst_byte, 192> avx_g2_pk_bytes = {};
    G2Affine scalar_g2_pk = {};
    AvxG2Affine avx_g2_pk = {};
    G2::to_pk2(scalar_g2_pk_bytes.data(), &scalar_g2_pk, key.data());
    AvxG2::to_pk2(avx_g2_pk_bytes.data(), &avx_g2_pk, key.data());
    EXPECT_EQ(avx_g2_pk_bytes, scalar_g2_pk_bytes);
    EXPECT_TRUE(G2Affine::equal(scalar_from_avx(avx_g2_pk), scalar_g2_pk));
}

TEST(AvxNativePointContract, PredecoratedOutputsFinalizeToStandardEncodings)
{
    std::array<blst_byte, 32> key = {};
    key[0] = 5;

    std::array<blst_byte, 96> decorated_g1 = {};
    AvxG1Affine avx_g1 = {};
    AvxG1::to_pk2(decorated_g1.data(), &avx_g1, key.data());
    const G1Affine scalar_g1 = scalar_from_avx(avx_g1);
    std::array<blst_byte, 48> compressed_g1 = {};
    std::memcpy(compressed_g1.data(), decorated_g1.data(), compressed_g1.size());
    compressed_g1[0] |= 0x80;
    std::array<blst_byte, 48> expected_compressed_g1 = {};
    G1Affine::compress(expected_compressed_g1.data(), scalar_g1);
    EXPECT_EQ(compressed_g1, expected_compressed_g1);
    decorated_g1[0] &= static_cast<blst_byte>(~0x20u);
    std::array<blst_byte, 96> expected_uncompressed_g1 = {};
    G1Affine::serialize(expected_uncompressed_g1.data(), scalar_g1);
    EXPECT_EQ(decorated_g1, expected_uncompressed_g1);

    std::array<blst_byte, 192> decorated_g2 = {};
    AvxG2Affine avx_g2 = {};
    AvxG2::to_pk2(decorated_g2.data(), &avx_g2, key.data());
    const G2Affine scalar_g2 = scalar_from_avx(avx_g2);
    std::array<blst_byte, 96> compressed_g2 = {};
    std::memcpy(compressed_g2.data(), decorated_g2.data(), compressed_g2.size());
    compressed_g2[0] |= 0x80;
    std::array<blst_byte, 96> expected_compressed_g2 = {};
    G2Affine::compress(expected_compressed_g2.data(), scalar_g2);
    EXPECT_EQ(compressed_g2, expected_compressed_g2);
    decorated_g2[0] &= static_cast<blst_byte>(~0x20u);
    std::array<blst_byte, 192> expected_uncompressed_g2 = {};
    G2Affine::serialize(expected_uncompressed_g2.data(), scalar_g2);
    EXPECT_EQ(decorated_g2, expected_uncompressed_g2);
}

TEST(AvxNativePointContract, TimesMinusZUsesDistinctObjects)
{
    G1 expected_g1 = {};
    blst::scalar::times_minus_z(expected_g1, G1::generator());
    AvxG1 actual_g1 = {};
    blst::avx::times_minus_z(actual_g1, AvxG1::generator());
    EXPECT_TRUE(G1::equal(expected_g1, scalar_from_avx(actual_g1)));

    G2 expected_g2 = {};
    blst::scalar::times_minus_z(expected_g2, G2::generator());
    AvxG2 actual_g2 = {};
    blst::avx::times_minus_z(actual_g2, AvxG2::generator());
    EXPECT_TRUE(G2::equal(expected_g2, scalar_from_avx(actual_g2)));
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

TEST(AvxCApiContract, RetainedG1GlvScratchAndWindowHooksMatchScalar)
{
    ASSERT_TRUE(blst_avx_enabled());

    constexpr size_t kCount = 2;
    std::array<blst_avx_p1_affine, kCount> points = {};
    for (auto& point : points) {
        ASSERT_TRUE(blst_avx_p1_to_affine(
            &point, blst_avx_p1_generator()));
    }

    std::array<blst_byte, 32 * kCount> scalars = {};
    scalars[0] = 3;
    scalars[32] = 5;

    const size_t bucket_bytes =
        blst_avx_p1_mult_pippenger_scratch_sizeof(kCount);
    const size_t bucket_used_bytes =
        blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof(kCount);
    const size_t split_points_bytes =
        blst_avx_p1_mult_pippenger_glv_points_scratch_sizeof(kCount);
    const size_t split_scalars_bytes =
        blst_avx_p1_mult_pippenger_glv_scalars_scratch_sizeof(kCount);
    ASSERT_GT(bucket_bytes, 0u);
    ASSERT_GT(bucket_used_bytes, 0u);
    ASSERT_EQ(split_points_bytes,
        2 * kCount * sizeof(blst_avx_p1_affine));
    ASSERT_EQ(split_scalars_bytes, kCount * 32);

    std::vector<blst_avx_p1> buckets(
        (bucket_bytes + sizeof(blst_avx_p1) - 1) /
        sizeof(blst_avx_p1));
    std::vector<blst_byte> bucket_used(bucket_used_bytes);
    std::vector<blst_avx_p1_affine> split_points(
        split_points_bytes / sizeof(blst_avx_p1_affine));
    std::vector<blst_byte> split_scalars(split_scalars_bytes);

    blst_avx_p1 automatic = {};
    ASSERT_TRUE(blst_avx_p1_mult_pippenger_glv_with_scratch(
        &automatic, points.data(), scalars.data(), kCount,
        buckets.data(), buckets.size() * sizeof(blst_avx_p1),
        bucket_used.data(), bucket_used.size(),
        split_points.data(), split_points_bytes,
        split_scalars.data(), split_scalars.size()));

    constexpr size_t kWindow = 4;
    const size_t window_bucket_bytes =
        blst_avx_p1_mult_pippenger_glv_window_scratch_sizeof(
            kCount, kWindow);
    const size_t window_bucket_used_bytes =
        blst_avx_p1_mult_pippenger_glv_window_bucket_used_scratch_sizeof(
            kCount, kWindow);
    ASSERT_GT(window_bucket_bytes, 0u);
    ASSERT_GT(window_bucket_used_bytes, 0u);

    std::vector<blst_avx_p1> window_buckets(
        (window_bucket_bytes + sizeof(blst_avx_p1) - 1) /
        sizeof(blst_avx_p1));
    std::vector<blst_byte> window_bucket_used(window_bucket_used_bytes);
    blst_avx_p1 explicit_window = {};
    ASSERT_TRUE(blst_avx_p1_mult_pippenger_glv_window_with_scratch(
        &explicit_window, points.data(), scalars.data(), kCount, kWindow,
        window_buckets.data(),
        window_buckets.size() * sizeof(blst_avx_p1),
        window_bucket_used.data(), window_bucket_used.size(),
        split_points.data(), split_points_bytes,
        split_scalars.data(), split_scalars.size()));

    std::array<blst_byte, 32> eight = {};
    eight[0] = 8;
    G1 expected = {};
    G1::mult_vartime(expected, G1::generator(), eight.data(), 256);
    blst_p1 actual_automatic = {};
    blst_p1 actual_window = {};
    ASSERT_TRUE(blst_p1_from_avx(&actual_automatic, &automatic));
    ASSERT_TRUE(blst_p1_from_avx(&actual_window, &explicit_window));
    EXPECT_TRUE(G1::equal(expected, actual_automatic));
    EXPECT_TRUE(G1::equal(expected, actual_window));
}

TEST(AvxCApiContract, PreparedG2LinesMatchDirectPairing)
{
    ASSERT_TRUE(blst_avx_enabled());

    const size_t storage_bytes = blst_avx_prepared_g2_lines_sizeof();
    const size_t storage_alignment = blst_avx_prepared_g2_lines_alignof();
    ASSERT_GT(storage_bytes, 0u);
    ASSERT_LE(storage_alignment, alignof(blst_avx_fp12));

    std::vector<blst_avx_fp12> storage(
        (storage_bytes + sizeof(blst_avx_fp12) - 1) /
        sizeof(blst_avx_fp12));
    ASSERT_TRUE(blst_avx_prepared_g2_lines_init(
        storage.data(),
        storage.size() * sizeof(blst_avx_fp12),
        blst_avx_p2_generator()));

    blst_avx_fp12 prepared = {};
    ASSERT_TRUE(blst_avx_pairing_prepared_g2(
        &prepared,
        storage.data(),
        storage.size() * sizeof(blst_avx_fp12),
        blst_avx_p1_generator()));

    blst_avx_fp12 direct = {};
    ASSERT_TRUE(blst_avx_pairing(
        &direct, blst_avx_p1_generator(), blst_avx_p2_generator()));
    EXPECT_TRUE(blst_avx_fp12_is_equal(&prepared, &direct));
}

} // namespace
