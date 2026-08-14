/*
 * Copyright Supranational LLC and Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "bindings/blst_scalar.h"
#include "blst/scalar/curve/map_to_G.h"
#include "blst/sign/aggregate.h"

#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstring>
#include <type_traits>

namespace {

static_assert(blst::kPairingBatchCapacity == 8);
static_assert(std::is_same_v<
    decltype(&::blst_pairing_get_miller_result),
    void (*)(blst_fp12*, blst_pairing*)>);

struct alignas(blst::PAIRING) PairingStorage
{
    std::array<unsigned char, sizeof(blst::PAIRING) + 64> bytes = {};

    blst_pairing* context()
    {
        return reinterpret_cast<blst_pairing*>(bytes.data());
    }

    blst::PAIRING* native_context()
    {
        return reinterpret_cast<blst::PAIRING*>(bytes.data());
    }
};

struct InvalidPairingPoints
{
    blst_p1_affine g1_non_subgroup = {};
    blst_p2_affine g2_non_subgroup = {};
    blst_p1_affine g1_off_curve = {};
    blst_p2_affine g2_off_curve = {};
};

InvalidPairingPoints invalid_pairing_points()
{
    using blst::scalar::G1AffinePoint;
    using blst::scalar::G1Point;
    using blst::scalar::G2AffinePoint;
    using blst::scalar::G2Point;

    InvalidPairingPoints points;
    G1Point g1_isogenous = {};
    G1Point g1_non_subgroup = {};
    blst::scalar::map_to_isogenous<G1Point>(
        g1_isogenous, G1AffinePoint::generator().x);
    blst::scalar::isogeny_map_to<G1Point>(g1_non_subgroup, g1_isogenous);
    G1Point::to_affine(points.g1_non_subgroup, g1_non_subgroup);

    G2Point g2_isogenous = {};
    G2Point g2_non_subgroup = {};
    blst::scalar::map_to_isogenous<G2Point>(
        g2_isogenous, G2AffinePoint::generator().x);
    blst::scalar::isogeny_map_to<G2Point>(g2_non_subgroup, g2_isogenous);
    G2Point::to_affine(points.g2_non_subgroup, g2_non_subgroup);

    points.g1_off_curve = *blst_p1_affine_generator();
    points.g1_off_curve.x[0] ^= 1;
    points.g2_off_curve = *blst_p2_affine_generator();
    points.g2_off_curve.x[0][0] ^= 1;

    return points;
}

TEST(AggregateAudit, EmptyContextExportsMillerIdentity)
{
    PairingStorage storage;
    blst_pairing_init(storage.context(), true, nullptr, 0);

    blst_fp12 result = {};
    blst_pairing_get_miller_result(&result, storage.context());

    EXPECT_TRUE(blst_fp12_is_equal(&result, blst_fp12_one()));
}

TEST(AggregateAudit, TrailingDstUsesRelocatableContextState)
{
    PairingStorage storage;
    auto* trailing_dst = storage.bytes.data() + blst_pairing_sizeof();
    trailing_dst[0] = 0x42;
    trailing_dst[1] = 0x17;

    blst_pairing_init(storage.context(), true, trailing_dst, 2);

    EXPECT_EQ(storage.native_context()->DST, nullptr);
    EXPECT_EQ(blst_pairing_get_dst(storage.context()), trailing_dst);
}

TEST(AggregateAudit, G1InfinityContributionPreservesAccumulator)
{
    blst_p1 identity = {};
    std::array<blst_byte, 48> encoded_identity = {};
    blst_p1_compress(encoded_identity.data(), &identity);

    const blst_p1 input = *blst_p1_generator();
    blst_p1 output = {};
    ASSERT_EQ(
        blst_aggregate_in_g1(&output, &input, encoded_identity.data()),
        BLST_SUCCESS);
    EXPECT_TRUE(blst_p1_is_equal(&output, &input));

    output = input;
    ASSERT_EQ(
        blst_aggregate_in_g1(&output, &output, encoded_identity.data()),
        BLST_SUCCESS);
    EXPECT_TRUE(blst_p1_is_equal(&output, &input));
}

TEST(AggregateAudit, G2InfinityContributionPreservesAccumulator)
{
    blst_p2 identity = {};
    std::array<blst_byte, 96> encoded_identity = {};
    blst_p2_compress(encoded_identity.data(), &identity);

    const blst_p2 input = *blst_p2_generator();
    blst_p2 output = {};
    ASSERT_EQ(
        blst_aggregate_in_g2(&output, &input, encoded_identity.data()),
        BLST_SUCCESS);
    EXPECT_TRUE(blst_p2_is_equal(&output, &input));

    output = input;
    ASSERT_EQ(
        blst_aggregate_in_g2(&output, &output, encoded_identity.data()),
        BLST_SUCCESS);
    EXPECT_TRUE(blst_p2_is_equal(&output, &input));
}

TEST(AggregateAudit, RejectedPkInG1ItemDoesNotMutateContext)
{
    PairingStorage storage;
    constexpr blst_byte dst[] = "aggregate failure atomicity";
    blst_pairing_init(storage.context(), true, dst, sizeof(dst) - 1);
    const auto before = storage.bytes;

    const blst_p1_affine infinite_pk = {};
    const blst_p2_affine signature = *blst_p2_affine_generator();
    constexpr blst_byte msg[] = "message";
    EXPECT_EQ(
        blst_pairing_chk_n_aggr_pk_in_g1(
            storage.context(),
            &infinite_pk,
            true,
            &signature,
            true,
            msg,
            sizeof(msg) - 1,
            nullptr,
            0),
        BLST_PK_IS_INFINITY);
    EXPECT_EQ(storage.bytes, before);
}

TEST(AggregateAudit, RejectedPkInG2ItemDoesNotMutateContext)
{
    PairingStorage storage;
    constexpr blst_byte dst[] = "aggregate failure atomicity";
    blst_pairing_init(storage.context(), true, dst, sizeof(dst) - 1);
    const auto before = storage.bytes;

    const blst_p2_affine infinite_pk = {};
    const blst_p1_affine signature = *blst_p1_affine_generator();
    constexpr blst_byte msg[] = "message";
    EXPECT_EQ(
        blst_pairing_chk_n_aggr_pk_in_g2(
            storage.context(),
            &infinite_pk,
            true,
            &signature,
            true,
            msg,
            sizeof(msg) - 1,
            nullptr,
            0),
        BLST_PK_IS_INFINITY);
    EXPECT_EQ(storage.bytes, before);
}

TEST(AggregateAudit, RequestedChecksRejectInvalidG1OrientationOperands)
{
    const InvalidPairingPoints invalid = invalid_pairing_points();
    ASSERT_TRUE(blst_p1_affine_on_curve(&invalid.g1_non_subgroup));
    ASSERT_FALSE(blst_p1_affine_in_g1(&invalid.g1_non_subgroup));
    ASSERT_TRUE(blst_p2_affine_on_curve(&invalid.g2_non_subgroup));
    ASSERT_FALSE(blst_p2_affine_in_g2(&invalid.g2_non_subgroup));
    ASSERT_FALSE(blst_p1_affine_on_curve(&invalid.g1_off_curve));
    ASSERT_FALSE(blst_p2_affine_on_curve(&invalid.g2_off_curve));

    PairingStorage storage;
    constexpr blst_byte dst[] = "pairing membership checks";
    constexpr blst_byte msg[] = "message";
    const blst_p1_affine& valid_pk = *blst_p1_affine_generator();
    const blst_p2_affine& valid_sig = *blst_p2_affine_generator();

    const auto expect_rejected_without_mutation =
        [&](const blst_p1_affine& pk, const blst_p2_affine& sig) {
            blst_pairing_init(
                storage.context(), true, dst, sizeof(dst) - 1);
            const auto before = storage.bytes;
            EXPECT_EQ(
                blst_pairing_chk_n_aggr_pk_in_g1(
                    storage.context(), &pk, true, &sig, true,
                    msg, sizeof(msg) - 1, nullptr, 0),
                BLST_POINT_NOT_IN_GROUP);
            EXPECT_EQ(storage.bytes, before);
        };

    expect_rejected_without_mutation(invalid.g1_non_subgroup, valid_sig);
    expect_rejected_without_mutation(invalid.g1_off_curve, valid_sig);
    expect_rejected_without_mutation(valid_pk, invalid.g2_non_subgroup);
    expect_rejected_without_mutation(valid_pk, invalid.g2_off_curve);

    const auto expect_default_signature_rejected =
        [&](const blst_p2_affine& sig) {
            blst_pairing_init(
                storage.context(), true, dst, sizeof(dst) - 1);
            const auto before = storage.bytes;
            EXPECT_EQ(
                blst_pairing_aggregate_pk_in_g1(
                    storage.context(), &valid_pk, &sig,
                    msg, sizeof(msg) - 1, nullptr, 0),
                BLST_POINT_NOT_IN_GROUP);
            EXPECT_EQ(storage.bytes, before);
        };
    expect_default_signature_rejected(invalid.g2_non_subgroup);
    expect_default_signature_rejected(invalid.g2_off_curve);
}

TEST(AggregateAudit, RequestedChecksRejectInvalidG2OrientationOperands)
{
    const InvalidPairingPoints invalid = invalid_pairing_points();
    PairingStorage storage;
    constexpr blst_byte dst[] = "pairing membership checks";
    constexpr blst_byte msg[] = "message";
    const blst_p2_affine& valid_pk = *blst_p2_affine_generator();
    const blst_p1_affine& valid_sig = *blst_p1_affine_generator();

    const auto expect_rejected_without_mutation =
        [&](const blst_p2_affine& pk, const blst_p1_affine& sig) {
            blst_pairing_init(
                storage.context(), true, dst, sizeof(dst) - 1);
            const auto before = storage.bytes;
            EXPECT_EQ(
                blst_pairing_chk_n_aggr_pk_in_g2(
                    storage.context(), &pk, true, &sig, true,
                    msg, sizeof(msg) - 1, nullptr, 0),
                BLST_POINT_NOT_IN_GROUP);
            EXPECT_EQ(storage.bytes, before);
        };

    expect_rejected_without_mutation(invalid.g2_non_subgroup, valid_sig);
    expect_rejected_without_mutation(invalid.g2_off_curve, valid_sig);
    expect_rejected_without_mutation(valid_pk, invalid.g1_non_subgroup);
    expect_rejected_without_mutation(valid_pk, invalid.g1_off_curve);

    const auto expect_default_signature_rejected =
        [&](const blst_p1_affine& sig) {
            blst_pairing_init(
                storage.context(), true, dst, sizeof(dst) - 1);
            const auto before = storage.bytes;
            EXPECT_EQ(
                blst_pairing_aggregate_pk_in_g2(
                    storage.context(), &valid_pk, &sig,
                    msg, sizeof(msg) - 1, nullptr, 0),
                BLST_POINT_NOT_IN_GROUP);
            EXPECT_EQ(storage.bytes, before);
        };
    expect_default_signature_rejected(invalid.g1_non_subgroup);
    expect_default_signature_rejected(invalid.g1_off_curve);
}

}  // namespace
