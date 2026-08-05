/*
 * Copyright Category Labs and Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include <gtest/gtest.h>

#include "bindings/blst_scalar.h"

/* Compile the internal HMAC/HKDF layer in its inherited raw-HKDF test mode. */
#define BLST_HKDF_TESTMODE
#include "blst/sign/keygen.cpp"

namespace
{
    unsigned char hex_nibble(char c)
    {
        if (c >= '0' && c <= '9')
            return (unsigned char)(c - '0');
        if (c >= 'a' && c <= 'f')
            return (unsigned char)(c - 'a' + 10);
        return (unsigned char)(c - 'A' + 10);
    }

    std::vector<unsigned char> bytes_from_hex(std::string_view hex)
    {
        EXPECT_EQ(hex.size() % 2, 0u);
        std::vector<unsigned char> bytes(hex.size() / 2);
        for (size_t i = 0; i < bytes.size(); ++i) {
            bytes[i] = (unsigned char)((hex_nibble(hex[2 * i]) << 4) |
                hex_nibble(hex[2 * i + 1]));
        }
        return bytes;
    }

    std::array<unsigned char, 32> hmac_sha256(
        const std::vector<unsigned char>& key,
        const std::vector<unsigned char>& message,
        size_t split)
    {
        blst::HmacSha256Context ctx;
        std::array<unsigned char, 32> digest{};
        blst::hmac_init(&ctx, key.data(), key.size());
        split = (std::min)(split, message.size());
        blst::hmac_update(&ctx, message.data(), split);
        blst::hmac_update(&ctx, message.data() + split,
            message.size() - split);
        blst::hmac_final(digest.data(), &ctx);
        return digest;
    }

    std::array<unsigned char, 32> scalar_bytes(const blst_scalar& scalar)
    {
        std::array<unsigned char, 32> bytes{};
        std::copy_n(scalar.b, bytes.size(), bytes.begin());
        return bytes;
    }

    std::array<unsigned char, 32> scalar_from_hex(std::string_view hex)
    {
        const auto decoded = bytes_from_hex(hex);
        EXPECT_EQ(decoded.size(), 32u);
        std::array<unsigned char, 32> bytes{};
        std::copy_n(decoded.begin(), bytes.size(), bytes.begin());
        return bytes;
    }
}

TEST(KeygenAudit, Rfc4231HmacSha256KnownAnswers)
{
    const std::vector<unsigned char> short_key(20, 0x0b);
    const std::vector<unsigned char> short_message = {
        'H', 'i', ' ', 'T', 'h', 'e', 'r', 'e'
    };
    const auto short_expected = bytes_from_hex(
        "b0344c61d8db38535ca8afceaf0bf12b"
        "881dc200c9833da726e9376c2e32cff7");

    const auto short_digest = hmac_sha256(short_key, short_message, 3);
    EXPECT_TRUE(std::equal(short_digest.begin(), short_digest.end(),
        short_expected.begin(), short_expected.end()));

    const std::vector<unsigned char> long_key(131, 0xaa);
    constexpr std::string_view long_text =
        "Test Using Larger Than Block-Size Key - Hash Key First";
    const std::vector<unsigned char> long_message(
        long_text.begin(), long_text.end());
    const auto long_expected = bytes_from_hex(
        "60e431591ee0b67f0d8a26aacbf5b77f"
        "8e0bc6213728c5140546040f0ee37f54");

    const auto long_digest = hmac_sha256(long_key, long_message, 29);
    EXPECT_TRUE(std::equal(long_digest.begin(), long_digest.end(),
        long_expected.begin(), long_expected.end()));
}

TEST(KeygenAudit, Rfc5869HkdfSha256KnownAnswer)
{
    const std::vector<unsigned char> ikm(22, 0x0b);
    const auto salt = bytes_from_hex("000102030405060708090a0b0c");
    const auto info = bytes_from_hex("f0f1f2f3f4f5f6f7f8f9");
    const auto expected_prk = bytes_from_hex(
        "077709362c2e32df0ddc3f0dc47bba63"
        "90b6c73bb50f9c3122ec844ad7c2b3e5");
    const auto expected_okm = bytes_from_hex(
        "3cb25f25faacd57a90434f64d0362f2a"
        "2d2d0a90cf1a5a4c5db02d56ecc4c5bf"
        "34007208d5b887185865");

    blst::HmacSha256Context ctx;
    std::array<unsigned char, 32> prk{};
    std::vector<unsigned char> okm(expected_okm.size());
    blst::hkdf_extract(prk.data(), salt.data(), salt.size(),
        ikm.data(), ikm.size(), &ctx);
    ASSERT_TRUE(blst::hkdf_expand(okm.data(), okm.size(), prk.data(),
        info.data(), info.size(), &ctx));

    EXPECT_TRUE(std::equal(prk.begin(), prk.end(),
        expected_prk.begin(), expected_prk.end()));
    EXPECT_EQ(okm, expected_okm);
}

TEST(KeygenAudit, HkdfExpandBoundsAndFailureAtomicity)
{
    std::array<unsigned char, 32> prk{};
    for (size_t i = 0; i < prk.size(); ++i)
        prk[i] = (unsigned char)(i * 7 + 3);

    blst::HmacSha256Context ctx;
    std::vector<unsigned char> maximum(255 * 32, 0xa5);
    ASSERT_TRUE(blst::hkdf_expand(maximum.data(), maximum.size(), prk.data(),
        nullptr, 0, &ctx));
    EXPECT_NE(maximum.front(), 0xa5);
    EXPECT_NE(maximum.back(), 0xa5);

    std::array<unsigned char, 4> unchanged = { 0x11, 0x22, 0x33, 0x44 };
    EXPECT_FALSE(blst::hkdf_expand(unchanged.data(), 255 * 32 + 1,
        prk.data(), nullptr, 0, &ctx));
    EXPECT_EQ(unchanged,
        (std::array<unsigned char, 4>{ 0x11, 0x22, 0x33, 0x44 }));

    EXPECT_TRUE(blst::hkdf_expand(nullptr, 0, nullptr, nullptr, 0, &ctx));
    EXPECT_FALSE(blst::hkdf_expand(nullptr, 1, prk.data(), nullptr, 0, &ctx));
    EXPECT_FALSE(blst::hkdf_expand(unchanged.data(), 1, prk.data(),
        nullptr, 1, &ctx));
}

#ifndef NDEBUG
TEST(KeygenAudit, DebugRejectsInvalidInternalPointerShapes)
{
    blst::HmacSha256Context ctx;
    const unsigned char byte = 1;
    std::array<unsigned char, 32> prk{};

    EXPECT_DEATH(blst::hmac_init(nullptr, &byte, 1), "ctx");
    blst::hmac_init(&ctx, &byte, 1);
    EXPECT_DEATH(blst::hmac_update(&ctx, nullptr, 1), "inp");
    EXPECT_DEATH(blst::hkdf_extract(prk.data(), nullptr, 1, &byte, 1, &ctx),
        "salt");
}
#endif

TEST(KeygenAudit, UpstreamVersionCompatibilityVectors)
{
    std::array<unsigned char, 32> ikm{};
    for (size_t i = 0; i < ikm.size(); ++i)
        ikm[i] = (unsigned char)i;
    const std::array<unsigned char, 9> info = {
        0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8
    };
    std::array<unsigned char, 17> salt{};
    for (size_t i = 0; i < salt.size(); ++i)
        salt[i] = (unsigned char)(0x10 + i);

    blst_scalar output{};
    blst_keygen_v3(&output, ikm.data(), ikm.size(), info.data(), info.size());
    EXPECT_EQ(scalar_bytes(output), scalar_from_hex(
        "ad4d7ac9e0d2d2cb001133103862f34a"
        "19ef56ba25ffb83d5805f0a2f908bc13"));

    blst_keygen(&output, ikm.data(), ikm.size(), info.data(), info.size());
    EXPECT_EQ(scalar_bytes(output), scalar_from_hex(
        "796099329e3b4ab11fd9fc80f167ecdb"
        "d6f745247e4372a1a980034e16bd6646"));

    blst_keygen_v4_5(&output, ikm.data(), ikm.size(), salt.data(), salt.size(),
        info.data(), info.size());
    EXPECT_EQ(scalar_bytes(output), scalar_from_hex(
        "e4a09699067f42a127c63ca2164e4cb6"
        "7456b042ffa1037a44a596affbae3109"));

    blst_keygen_v5(&output, ikm.data(), ikm.size(), salt.data(), salt.size(),
        info.data(), info.size());
    EXPECT_EQ(scalar_bytes(output), scalar_from_hex(
        "f6305cdb1bfe7b6d2ef381f88f95ff29"
        "ce01e18fe6c3b6b0468ecd9bbd6f556d"));

    const unsigned char nonnull_empty_salt = 0;
    blst_keygen_v5(&output, ikm.data(), ikm.size(), &nonnull_empty_salt, 0,
        info.data(), info.size());
    EXPECT_EQ(scalar_bytes(output), scalar_from_hex(
        "dbc8994340beee9a1427882ceef88a0a"
        "f000576fe93d33f1f863271c94672705"));
}

TEST(KeygenAudit, UpstreamNullAndAliasSemantics)
{
    std::array<unsigned char, 32> ikm{};
    for (size_t i = 0; i < ikm.size(); ++i)
        ikm[i] = (unsigned char)(0x80 + i);

    blst_scalar default_key{};
    blst_scalar compatibility_key{};
    blst_keygen(&default_key, ikm.data(), ikm.size(), nullptr, 999);
    blst_keygen_v4_5(&compatibility_key, ikm.data(), ikm.size(), nullptr, 999,
        nullptr, 123);
    EXPECT_EQ(scalar_bytes(default_key), scalar_bytes(compatibility_key));

    blst_scalar separate{};
    blst_derive_child_eip2333(&separate, &default_key, 42);
    blst_scalar in_place = default_key;
    blst_derive_child_eip2333(&in_place, &in_place, 42);
    EXPECT_EQ(scalar_bytes(in_place), scalar_bytes(separate));
}

TEST(KeygenAudit, PublishedEip2333MasterAndChildVectors)
{
    struct Vector
    {
        const char* seed;
        const char* master;
        uint32_t child_index;
        const char* child;
    };
    const Vector vectors[] = {
        {
            "c55257c360c07c72029aebc1b53c05ed0362ada38ead3e3e9efa3708e5349553"
            "1f09a6987599d18264c1e1c92f2cf141630c7a3c4ab7c81b2f001698e7463b04",
            "0d7359d57963ab8fbbde1852dcf553fedbc31f464d80ee7d40ae683122b45070",
            0,
            "2d18bd6c14e6d15bf8b5085c9b74f3daae3b03cc2014770a599d8c1539e50f8e",
        },
        {
            "3141592653589793238462643383279502884197169399375105820974944592",
            "41c9e07822b092a93fd6797396338c3ada4170cc81829fdfce6b5d34bd5e7ec7",
            3141592653u,
            "384843fad5f3d777ea39de3e47a8f999ae91f89e42bffa993d91d9782d152a0f",
        },
        {
            "0099ff991111002299dd7744ee3355bbdd8844115566cc55663355668888cc00",
            "3cfa341ab3910a7d00d933d8f7c4fe87c91798a0397421d6b19fd5b815132e80",
            0xffffffffu,
            "40e86285582f35b28821340f6a53b448588efa575bc4d88c32ef8567b8d9479b",
        },
        {
            "d4e56740f876aef8c010b86a40d5f56745a118d0906a34e69aec8c0db1cb8fa3",
            "2a0e28ffa5fbbe2f8e7aad4ed94f745d6bf755c51182e119bb1694fe61d3afca",
            42,
            "455c0dc9fccb3395825d92a60d2672d69416be1c2578a87a7a3d3ced11ebb88d",
        },
    };

    for (const Vector& vector : vectors) {
        const auto seed = bytes_from_hex(vector.seed);
        const auto expected_master = bytes_from_hex(vector.master);
        const auto expected_child = bytes_from_hex(vector.child);
        blst_scalar master{};
        blst_scalar child{};
        std::array<unsigned char, 32> encoded{};

        blst_derive_master_eip2333(&master, seed.data(), seed.size());
        blst_bendian_from_scalar(encoded.data(), &master);
        EXPECT_TRUE(std::equal(encoded.begin(), encoded.end(),
            expected_master.begin(), expected_master.end()));

        blst_derive_child_eip2333(&child, &master, vector.child_index);
        blst_bendian_from_scalar(encoded.data(), &child);
        EXPECT_TRUE(std::equal(encoded.begin(), encoded.end(),
            expected_child.begin(), expected_child.end()));
    }
}
