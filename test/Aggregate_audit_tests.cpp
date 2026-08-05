/*
 * Copyright Supranational LLC and Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "bindings/blst_scalar.h"
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

}  // namespace
