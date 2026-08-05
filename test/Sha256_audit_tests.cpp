/*
 * Copyright Category Labs and Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include <gtest/gtest.h>

#include "bindings/blst_scalar.h"
#include "blst/sign/sha256.h"

namespace
{
    using Digest = std::array<unsigned char, 32>;

    unsigned char hex_nibble(char c)
    {
        if (c >= '0' && c <= '9')
            return (unsigned char)(c - '0');
        if (c >= 'a' && c <= 'f')
            return (unsigned char)(c - 'a' + 10);
        return (unsigned char)(c - 'A' + 10);
    }

    Digest digest_from_hex(std::string_view hex)
    {
        EXPECT_EQ(hex.size(), 64u);
        Digest digest{};
        for (size_t i = 0; i < digest.size(); ++i)
            digest[i] = (unsigned char)((hex_nibble(hex[2 * i]) << 4) |
                hex_nibble(hex[2 * i + 1]));
        return digest;
    }

    Digest hash_one_shot(const void* data, size_t size)
    {
        Digest digest{};
        blst_sha256(digest.data(), data, size);
        return digest;
    }

    Digest hash_streaming(const unsigned char* data, size_t size,
        const std::vector<size_t>& schedule)
    {
        blst::Sha256Context ctx;
        Digest digest{};
        blst::sha256_init(&ctx);
        blst::sha256_update(&ctx, nullptr, 0);

        size_t offset = 0;
        size_t turn = 0;
        while (offset != size) {
            const size_t chunk =
                (std::min)(schedule[turn % schedule.size()], size - offset);
            blst::sha256_update(&ctx, data + offset, chunk);
            offset += chunk;
            ++turn;
        }

        blst::sha256_final(digest.data(), &ctx);
        blst::sha256_clear(&ctx);
        return digest;
    }

    std::vector<unsigned char> patterned_message(size_t size)
    {
        std::vector<unsigned char> message(size);
        for (size_t i = 0; i < size; ++i)
            message[i] = (unsigned char)((i * 37 + 11) & 0xff);
        return message;
    }
}

TEST(Sha256Audit, StandardKnownAnswers)
{
    struct Vector {
        std::string message;
        const char* digest;
    };
    const Vector vectors[] = {
        { "", "e3b0c44298fc1c149afbf4c8996fb924"
              "27ae41e4649b934ca495991b7852b855" },
        { "abc", "ba7816bf8f01cfea414140de5dae2223"
                 "b00361a396177a9cb410ff61f20015ad" },
        { "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
          "248d6a61d20638b8e5c026930c3e6039"
          "a33ce45964ff2167f6ecedd419db06c1" },
    };

    for (const auto& vector : vectors) {
        const auto expected = digest_from_hex(vector.digest);
        EXPECT_EQ(hash_one_shot(vector.message.data(), vector.message.size()),
            expected);
        EXPECT_EQ(hash_streaming(
            (const unsigned char*)vector.message.data(), vector.message.size(),
            { 1 }), expected);
        EXPECT_EQ(hash_streaming(
            (const unsigned char*)vector.message.data(), vector.message.size(),
            { 7, 64, 3, 19 }), expected);
    }

    const std::string million_a(1000000, 'a');
    const auto million_expected = digest_from_hex(
        "cdc76e5c9914fb9281a1c7e284d73e67"
        "f1809a48a497200e046d39ccc7112cd0");
    EXPECT_EQ(hash_one_shot(million_a.data(), million_a.size()),
        million_expected);
    EXPECT_EQ(hash_streaming((const unsigned char*)million_a.data(),
        million_a.size(), { 1, 31, 64, 257, 4096 }), million_expected);
}

TEST(Sha256Audit, PaddingAndBlockBoundaries)
{
    struct Vector {
        size_t size;
        const char* digest;
    };
    const Vector vectors[] = {
        { 55, "2900465fcb533e05a158fd2b3be0e5e3b03740d83060aa3580e0d98a96bf2384" },
        { 56, "31454ff48ef36af2f08fd511bdc37d9d5855ac23e992e5ff5445cb6b7674a674" },
        { 63, "5f6401b96532c36de4e65beec0409b69b1d181864c8009b7a04f43e5d56350d1" },
        { 64, "94eb5de4943613fd048dc93393ab06877405faa39c11f53e9386083339833e7e" },
        { 65, "fc518669b6eb4b4dd91827ecacef86689c725bd5bab888fd3b26dbb196eec954" },
        { 119, "b0dc41b1a384e2f1203f0351b38fbeaafceef577ce1191d5bfc25da39f721eae" },
        { 120, "5df24dd802ac26132ce608dcb5f09841eef039ee0f152acf98d26d17fe4e88e6" },
        { 127, "0fe729ff19257bd6fec853acc2ea355f6b34b58e6c0f684c3e188fcdfcd9baae" },
        { 128, "0aedd4856f8eba0963627336ad5144a9a7dbe12498e6066f0165fc97d8ddee4c" },
    };

    for (const auto& vector : vectors) {
        const auto message = patterned_message(vector.size);
        const auto expected = digest_from_hex(vector.digest);
        EXPECT_EQ(hash_one_shot(message.data(), message.size()), expected);
        EXPECT_EQ(hash_streaming(message.data(), message.size(), { 1 }), expected);
        EXPECT_EQ(hash_streaming(message.data(), message.size(),
            { 13, 51, 2, 64 }), expected);
    }
}

TEST(Sha256Audit, NullEmptyAndOutputOverlap)
{
    const auto empty = digest_from_hex(
        "e3b0c44298fc1c149afbf4c8996fb924"
        "27ae41e4649b934ca495991b7852b855");
    EXPECT_EQ(hash_one_shot(nullptr, 0), empty);

    std::array<unsigned char, 64> storage{};
    storage[0] = 'a';
    storage[1] = 'b';
    storage[2] = 'c';
    const auto expected = digest_from_hex(
        "ba7816bf8f01cfea414140de5dae2223"
        "b00361a396177a9cb410ff61f20015ad");
    blst_sha256(storage.data(), storage.data(), 3);
    EXPECT_TRUE(std::equal(expected.begin(), expected.end(), storage.begin()));

    storage.fill(0);
    storage[0] = 'a';
    storage[1] = 'b';
    storage[2] = 'c';
    blst_sha256(storage.data() + 1, storage.data(), 3);
    EXPECT_TRUE(std::equal(expected.begin(), expected.end(), storage.begin() + 1));
}

TEST(Sha256Audit, ExplicitClearErasesCompleteContext)
{
    blst::Sha256Context ctx;
    const unsigned char message[] = { 1, 2, 3, 4, 5 };
    blst::sha256_init(&ctx);
    blst::sha256_update(&ctx, message, sizeof(message));
    blst::sha256_clear(&ctx);

    const unsigned char* bytes = (const unsigned char*)&ctx;
    EXPECT_TRUE(std::all_of(bytes, bytes + sizeof(ctx),
        [](unsigned char value) { return value == 0; }));
}

#ifndef NDEBUG
TEST(Sha256Audit, DebugRejectsInvalidContextStateAndLengthDomain)
{
    blst::Sha256Context ctx;
    Digest digest{};
    blst::sha256_init(&ctx);
    ctx.buffered_bytes = sizeof(ctx.buffer);
    EXPECT_DEATH(blst::sha256_final(digest.data(), &ctx), "buffered_bytes");

    blst::sha256_init(&ctx);
    ctx.total_bytes = blst::sha256_max_message_bytes;
    ctx.buffered_bytes =
        (size_t)(ctx.total_bytes & (sizeof(ctx.buffer) - 1));
    const unsigned char byte = 0;
    EXPECT_DEATH(blst::sha256_update(&ctx, &byte, 1), "sha256_max_message_bytes");
}
#endif
