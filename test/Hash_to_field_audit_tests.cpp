/*
 * Copyright Supranational LLC and Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "bindings/blst.h"
#include "blst/scalar/fields/hash_to_field.h"

namespace {

using blst::scalar::Fp;
using blst::scalar::Fp2;

unsigned char hex_nibble(char c)
{
    if (c >= '0' && c <= '9')
        return static_cast<unsigned char>(c - '0');
    if (c >= 'a' && c <= 'f')
        return static_cast<unsigned char>(c - 'a' + 10);
    if (c >= 'A' && c <= 'F')
        return static_cast<unsigned char>(c - 'A' + 10);
    ADD_FAILURE() << "invalid hexadecimal test vector";
    return 0;
}

std::vector<unsigned char> from_hex(std::string_view hex)
{
    EXPECT_EQ(hex.size() % 2, 0u);
    std::vector<unsigned char> out(hex.size() / 2);
    for (size_t i = 0; i < out.size(); ++i)
        out[i] = static_cast<unsigned char>(
            (hex_nibble(hex[2 * i]) << 4) | hex_nibble(hex[2 * i + 1]));
    return out;
}

Fp fp_from_hex(const char* hex)
{
    Fp out = {};
    EXPECT_TRUE(blst_fp_from_hexascii(
        &out, reinterpret_cast<const blst_byte*>(hex)));
    return out;
}

std::string fp_to_hex(const Fp& fp)
{
    static constexpr char digits[] = "0123456789abcdef";
    std::array<unsigned char, 48> bytes = {};
    blst_bendian_from_fp(bytes.data(), &fp);
    std::string out = "0x";
    out.reserve(98);
    for (unsigned char byte : bytes) {
        out.push_back(digits[byte >> 4]);
        out.push_back(digits[byte & 15]);
    }
    return out;
}

void expect_xmd(std::string_view expected_hex,
    const unsigned char* msg, size_t msg_len,
    const unsigned char* dst, size_t dst_len)
{
    const std::vector<unsigned char> expected = from_hex(expected_hex);
    std::vector<unsigned char> actual(expected.size());
    blst_expand_message_xmd(actual.data(), actual.size(),
        msg, msg_len, dst, dst_len);
    EXPECT_EQ(actual, expected);
}

TEST(HashToFieldAudit, Rfc9380XmdSha256Vectors)
{
    static constexpr unsigned char dst[] =
        "QUUX-V01-CS02-with-expander-SHA256-128";
    static constexpr unsigned char abc[] = "abc";

    expect_xmd(
        "68a985b87eb6b46952128911f2a4412bbc302a9d759667f87f7a21d803f07235",
        nullptr, 0, dst, sizeof(dst) - 1);
    expect_xmd(
        "abba86a6129e366fc877aab32fc4ffc70120d8996c88aee2fe4b32d6c7b6437a"
        "647e6c3163d40b76a73cf6a5674ef1d890f95b664ee0afa5359a5c4e07985635"
        "bbecbac65d747d3d2da7ec2b8221b17b0ca9dc8a1ac1c07ea6a1e60583e2cb00"
        "058e77b7b72a298425cd1b941ad4ec65e8afc50303a22c0f99b0509b4c895f40",
        abc, sizeof(abc) - 1, dst, sizeof(dst) - 1);
}

TEST(HashToFieldAudit, Rfc9380OversizeDstVector)
{
    static constexpr unsigned char dst[] =
        "QUUX-V01-CS02-with-expander-SHA256-128-long-DST-"
        "1111111111111111111111111111111111111111111111111111111111111111"
        "1111111111111111111111111111111111111111111111111111111111111111"
        "1111111111111111111111111111111111111111111111111111111111111111"
        "1111111111111111";
    static_assert(sizeof(dst) - 1 == 256);

    expect_xmd(
        "e8dc0c8b686b7ef2074086fbdd2f30e3f8bfbd3bdf177f73f04b97ce618a3ed3",
        nullptr, 0, dst, sizeof(dst) - 1);
}

TEST(HashToFieldAudit, ExactPartialOutputSupportsUnalignedByteBuffers)
{
    static constexpr unsigned char msg[] = "abc";
    static constexpr unsigned char dst[] = "BLS-XMD";
    const std::vector<unsigned char> expected = from_hex(
        "e7a0efff4c32ff1dd44e310a93616a8dacd80ef0c4be509debedde68aec3cd61"
        "7d795127f40c6fb09aa35be8ee1f2938aa5b5bd4646e9797ce4f8796b11eeea6"
        "6f");
    std::array<unsigned char, 67> storage;
    storage.fill(0xa5);

    blst_expand_message_xmd(storage.data() + 1, expected.size(),
        msg, sizeof(msg) - 1, dst, sizeof(dst) - 1);

    EXPECT_EQ(storage.front(), 0xa5);
    EXPECT_EQ(storage.back(), 0xa5);
    EXPECT_TRUE(std::equal(expected.begin(), expected.end(), storage.begin() + 1));
}

TEST(HashToFieldAudit, OutputLengthBoundariesFailClosed)
{
    static constexpr unsigned char msg[] = "abc";
    static constexpr unsigned char dst[] = "BLS-XMD";
    std::array<unsigned char, 1> untouched = {0xa5};
    blst_expand_message_xmd(nullptr, 0, nullptr, 0, nullptr, 0);
    blst_expand_message_xmd(untouched.data(), 0, msg, sizeof(msg) - 1,
        dst, sizeof(dst) - 1);
    EXPECT_EQ(untouched[0], 0xa5);

    std::vector<unsigned char> maximum(blst::scalar::xmd_sha256_max_bytes + 2, 0xa5);
    blst_expand_message_xmd(maximum.data() + 1,
        blst::scalar::xmd_sha256_max_bytes,
        msg, sizeof(msg) - 1, dst, sizeof(dst) - 1);
    EXPECT_EQ(maximum.front(), 0xa5);
    EXPECT_EQ(maximum.back(), 0xa5);
    EXPECT_TRUE(std::any_of(maximum.begin() + 1, maximum.end() - 1,
        [](unsigned char byte) { return byte != 0; }));

    std::vector<unsigned char> excessive(
        blst::scalar::xmd_sha256_max_bytes + 1, 0xa5);
    blst_expand_message_xmd(excessive.data(), excessive.size(),
        msg, sizeof(msg) - 1, dst, sizeof(dst) - 1);
    EXPECT_TRUE(std::all_of(excessive.begin(), excessive.end(),
        [](unsigned char byte) { return byte == 0; }));
}

TEST(HashToFieldAudit, Sha256InputDomainBoundariesFailClosed)
{
    constexpr size_t max_bytes =
        static_cast<size_t>(blst::sha256_max_message_bytes);
    constexpr size_t fixed_bytes =
        blst::scalar::hash_to_field_detail::xmd_b0_fixed_bytes;
    constexpr size_t dst_prefix_bytes =
        blst::scalar::hash_to_field_detail::xmd_oversize_dst_prefix_bytes;
    using blst::scalar::hash_to_field_detail::xmd_input_lengths_valid;

    EXPECT_TRUE(xmd_input_lengths_valid(0, max_bytes - fixed_bytes, 0));
    EXPECT_FALSE(xmd_input_lengths_valid(0, max_bytes - fixed_bytes + 1, 0));

    EXPECT_TRUE(xmd_input_lengths_valid(
        max_bytes - fixed_bytes - 255, 0, 255));
    EXPECT_FALSE(xmd_input_lengths_valid(
        max_bytes - fixed_bytes - 255 + 1, 0, 255));

    EXPECT_TRUE(xmd_input_lengths_valid(0, 0, max_bytes - dst_prefix_bytes));
    EXPECT_FALSE(xmd_input_lengths_valid(
        0, 0, max_bytes - dst_prefix_bytes + 1));

    EXPECT_TRUE(xmd_input_lengths_valid(max_bytes - fixed_bytes, 0, 0));
    EXPECT_FALSE(xmd_input_lengths_valid(max_bytes - fixed_bytes + 1, 0, 0));

    const unsigned char dummy = 0;
    std::array<unsigned char, 32> output;

    output.fill(0xa5);
    blst_expand_message_xmd(output.data(), output.size(),
        &dummy, max_bytes - fixed_bytes + 1, nullptr, 0);
    EXPECT_TRUE(std::all_of(output.begin(), output.end(),
        [](unsigned char byte) { return byte == 0; }));

    output.fill(0xa5);
    blst_expand_message_xmd(output.data(), output.size(),
        nullptr, 0, &dummy, max_bytes - dst_prefix_bytes + 1);
    EXPECT_TRUE(std::all_of(output.begin(), output.end(),
        [](unsigned char byte) { return byte == 0; }));
}

#ifdef NDEBUG
TEST(HashToFieldAudit, InvalidInputPointerLengthPairsClearOutputInRelease)
{
    static constexpr unsigned char msg[] = "abc";
    static constexpr unsigned char dst[] = "BLS-XMD";
    std::array<unsigned char, 33> output;

    output.fill(0xa5);
    blst_expand_message_xmd(output.data(), output.size(),
        nullptr, 1, dst, sizeof(dst) - 1);
    EXPECT_TRUE(std::all_of(output.begin(), output.end(),
        [](unsigned char byte) { return byte == 0; }));

    output.fill(0xa5);
    blst_expand_message_xmd(output.data(), output.size(),
        msg, sizeof(msg) - 1, nullptr, 1);
    EXPECT_TRUE(std::all_of(output.begin(), output.end(),
        [](unsigned char byte) { return byte == 0; }));
}
#endif

TEST(HashToFieldAudit, Rfc9380G1HashToFieldIntermediate)
{
    static constexpr unsigned char dst[] =
        "QUUX-V01-CS02-with-BLS12381G1_XMD:SHA-256_SSWU_RO_";
    Fp actual[2] = {};
    blst::scalar::hash_to_field<2>(actual,
        nullptr, 0, nullptr, 0, dst, sizeof(dst) - 1);

    EXPECT_TRUE(Fp::equal(actual[0], fp_from_hex(
        "0x0ba14bd907ad64a016293ee7c2d276b8eae71f25a4b941eece7b0d89f17f75cb3ae5438a614fb61d6835ad59f29c564f")));
    EXPECT_TRUE(Fp::equal(actual[1], fp_from_hex(
        "0x019b9bd7979f12657976de2884c7cce192b82c177c80e0ec604436a7f538d231552f0d96d9f7babe5fa3b19b3ff25ac9")));
}

TEST(HashToFieldAudit, Rfc9380G2HashToFieldIntermediate)
{
    static constexpr unsigned char dst[] =
        "QUUX-V01-CS02-with-BLS12381G2_XMD:SHA-256_SSWU_RO_";
    Fp2 actual[2] = {};
    blst::scalar::hash_to_field<4>(actual,
        nullptr, 0, nullptr, 0, dst, sizeof(dst) - 1);

    const Fp2 expected[2] = {
        {
            fp_from_hex("0x03dbc2cce174e91ba93cbb08f26b917f98194a2ea08d1cce75b2b9cc9f21689d80bd79b594a613d0a68eb807dfdc1cf8"),
            fp_from_hex("0x05a2acec64114845711a54199ea339abd125ba38253b70a92c876df10598bd1986b739cad67961eb94f7076511b3b39a"),
        },
        {
            fp_from_hex("0x02f99798e8a5acdeed60d7e18e9120521ba1f47ec090984662846bc825de191b5b7641148c0dbc237726a334473eee94"),
            fp_from_hex("0x145a81e418d4010cc027a68f14391b30074e89e60ee7a22f87217b2f6eb0c4b94c9115b436e6fa4607e95a98de30a435"),
        },
    };
    EXPECT_EQ(fp_to_hex(actual[0][0]), fp_to_hex(expected[0][0]));
    EXPECT_EQ(fp_to_hex(actual[0][1]), fp_to_hex(expected[0][1]));
    EXPECT_EQ(fp_to_hex(actual[1][0]), fp_to_hex(expected[1][0]));
    EXPECT_EQ(fp_to_hex(actual[1][1]), fp_to_hex(expected[1][1]));
}

} // namespace
