/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "bindings/blst.h"
#include "blst/scalar/curve/map_to_G.h"

namespace {

using G1 = blst::scalar::G1Point;
using G1Affine = blst::scalar::G1AffinePoint;
using G2 = blst::scalar::G2Point;
using G2Affine = blst::scalar::G2AffinePoint;

std::array<blst_byte, 32> maximal_scalar()
{
    std::array<blst_byte, 32> scalar = {};
    scalar.fill(0xff);
    return scalar;
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
