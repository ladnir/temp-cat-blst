/*
 * Copyright Supranational LLC and Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/curve/GT.h"
#include "blst/common/bytes.h"

#include <gtest/gtest.h>

#include <array>
#include <type_traits>
#include <utility>

namespace {

using ScalarFp12 = blst::scalar::Fp12;
using ScalarGt = blst::scalar::GT;

constexpr std::size_t kEncodedGtSize = 12 * 48;

template<class Group>
std::array<unsigned char, kEncodedGtSize> serialize(const Group& value)
{
    std::array<unsigned char, kEncodedGtSize> bytes = {};
    Group::serialize(bytes.data(), value);
    return bytes;
}

template<class Group>
void expect_checked_roundtrip_and_noncanonical_rejection()
{
    static_assert(std::is_same_v<
        decltype(Group::deserialize(
            std::declval<Group&>(),
            std::declval<const unsigned char*>())),
        bool>);

    const Group generator = Group::generator();
    const auto canonical = serialize(generator);

    Group decoded = Group::identity();
    ASSERT_TRUE(Group::deserialize(decoded, canonical.data()));
    EXPECT_TRUE(Group::equal(decoded, generator));

    const Group identity = Group::identity();
    const auto identity_bytes = serialize(identity);
    decoded = generator;
    ASSERT_TRUE(Group::deserialize(decoded, identity_bytes.data()));
    EXPECT_TRUE(Group::is_identity(decoded));

    for (std::size_t coefficient = 0; coefficient < 12; ++coefficient) {
        auto noncanonical = canonical;
        blst::be_bytes_from_limbs(
            noncanonical.data() + coefficient * 48,
            blst::scalar::Fp::modulus().l,
            sizeof(blst::scalar::Fp));

        Group checked_out = generator;
        EXPECT_FALSE(Group::deserialize(checked_out, noncanonical.data()));
        EXPECT_TRUE(Group::equal(checked_out, generator));

        Group unchecked_out = generator;
        EXPECT_FALSE(
            Group::deserialize_unchecked(unchecked_out, noncanonical.data()));
        EXPECT_TRUE(Group::equal(unchecked_out, generator));
    }
}

TEST(GtDeserializeAudit, ModulusMinusOneIsCanonical)
{
    std::array<unsigned char, kEncodedGtSize> bytes = {};
    blst::scalar::Fp modulus_minus_one = blst::scalar::Fp::modulus();
    --modulus_minus_one.l[0];
    blst::be_bytes_from_limbs(
        bytes.data(), modulus_minus_one.l, sizeof(modulus_minus_one));

    ScalarGt candidate = ScalarGt::identity();
    ASSERT_TRUE(ScalarGt::deserialize_unchecked(candidate, bytes.data()));
}

ScalarFp12 make_general_fp12_value()
{
    ScalarFp12 value = {};
    blst::scalar::Fp::add(
        value[0][0][0],
        blst::scalar::Fp::one(),
        blst::scalar::Fp::one());
    return value;
}

TEST(GtDeserializeAudit, ScalarAndSelectedCheckedRoundtrip)
{
    expect_checked_roundtrip_and_noncanonical_rejection<ScalarGt>();
    expect_checked_roundtrip_and_noncanonical_rejection<blst::GT>();
}

TEST(GtDeserializeAudit, CheckedRejectsCanonicalNonGtWithoutMutation)
{
    const ScalarFp12 general_value = make_general_fp12_value();
    ASSERT_FALSE(ScalarFp12::in_group(general_value));

    std::array<unsigned char, kEncodedGtSize> bytes = {};
    ScalarFp12::to_bendian(bytes.data(), general_value);

    const ScalarGt original = ScalarGt::generator();
    ScalarGt out = original;
    EXPECT_FALSE(ScalarGt::deserialize(out, bytes.data()));
    EXPECT_TRUE(ScalarGt::equal(out, original));
}

TEST(GtDeserializeAudit, UncheckedAcceptsCanonicalCandidateOnly)
{
    const ScalarFp12 general_value = make_general_fp12_value();
    std::array<unsigned char, kEncodedGtSize> bytes = {};
    ScalarFp12::to_bendian(bytes.data(), general_value);

    ScalarGt candidate = ScalarGt::identity();
    ASSERT_TRUE(ScalarGt::deserialize_unchecked(candidate, bytes.data()));
    EXPECT_TRUE(ScalarFp12::equal(candidate.data, general_value));
    EXPECT_FALSE(ScalarGt::in_group(candidate));
}

TEST(GtDeserializeAudit, DeferredAggregateMembershipCheck)
{
    const ScalarFp12 value = make_general_fp12_value();
    ScalarFp12 inverse = {};
    ScalarFp12::inverse(inverse, value);

    std::array<unsigned char, kEncodedGtSize> value_bytes = {};
    std::array<unsigned char, kEncodedGtSize> inverse_bytes = {};
    ScalarFp12::to_bendian(value_bytes.data(), value);
    ScalarFp12::to_bendian(inverse_bytes.data(), inverse);

    ScalarGt left;
    ScalarGt right;
    ASSERT_TRUE(ScalarGt::deserialize_unchecked(left, value_bytes.data()));
    ASSERT_TRUE(ScalarGt::deserialize_unchecked(right, inverse_bytes.data()));
    ASSERT_FALSE(ScalarGt::in_group(left));
    ASSERT_FALSE(ScalarGt::in_group(right));

    ScalarFp12 product = {};
    ScalarFp12::mul(product, left.data, right.data);
    ASSERT_TRUE(ScalarFp12::in_group(product));

    const ScalarGt validated(product);
    EXPECT_TRUE(ScalarGt::is_identity(validated));
}

}  // namespace
