#include "BTE/building_blocks/hkdf_sha256.h"

#include <algorithm>
#include <stdexcept>

#include "blst/sign/sha256.h"

namespace bte::building_blocks
{
    namespace
    {
        constexpr std::size_t kSha256BlockSize = 64;
        constexpr std::size_t kSha256OutputSize = 32;

        std::array<std::uint8_t, kSha256OutputSize> hmac_sha256(
            const std::span<const std::uint8_t> key,
            const std::span<const std::uint8_t> message)
        {
            std::array<std::uint8_t, kSha256BlockSize> key_block{};
            if (key.size() > kSha256BlockSize) {
                const std::array<std::uint8_t, kSha256OutputSize> hashed = sha256(key);
                std::copy(hashed.begin(), hashed.end(), key_block.begin());
            } else {
                std::copy(key.begin(), key.end(), key_block.begin());
            }

            std::array<std::uint8_t, kSha256BlockSize> ipad = key_block;
            std::array<std::uint8_t, kSha256BlockSize> opad = key_block;
            for (std::size_t i = 0; i < kSha256BlockSize; ++i) {
                ipad[i] ^= 0x36U;
                opad[i] ^= 0x5cU;
            }

            blst::SHA256_CTX inner{};
            blst::sha256_init(&inner);
            blst::sha256_update(&inner, ipad.data(), ipad.size());
            if (!message.empty()) {
                blst::sha256_update(&inner, message.data(), message.size());
            }
            std::array<std::uint8_t, kSha256OutputSize> inner_hash{};
            blst::sha256_final(inner_hash.data(), &inner);

            blst::SHA256_CTX outer{};
            blst::sha256_init(&outer);
            blst::sha256_update(&outer, opad.data(), opad.size());
            blst::sha256_update(&outer, inner_hash.data(), inner_hash.size());
            std::array<std::uint8_t, kSha256OutputSize> mac{};
            blst::sha256_final(mac.data(), &outer);
            return mac;
        }
    }

    std::array<std::uint8_t, 32> sha256(const std::span<const std::uint8_t> message)
    {
        blst::SHA256_CTX ctx{};
        blst::sha256_init(&ctx);
        if (!message.empty()) {
            blst::sha256_update(&ctx, message.data(), message.size());
        }
        std::array<std::uint8_t, 32> digest{};
        blst::sha256_final(digest.data(), &ctx);
        return digest;
    }

    std::vector<std::uint8_t> hkdf_sha256(
        const std::span<const std::uint8_t> salt,
        const std::span<const std::uint8_t> ikm,
        const std::span<const std::uint8_t> info,
        const std::size_t output_len)
    {
        if (output_len == 0) {
            return {};
        }
        if (output_len > 255 * kSha256OutputSize) {
            throw std::invalid_argument("HKDF output too long");
        }

        static constexpr std::array<std::uint8_t, kSha256OutputSize> kZeroSalt{};
        const std::span<const std::uint8_t> extract_salt =
            salt.empty() ? std::span<const std::uint8_t>(kZeroSalt.data(), kZeroSalt.size()) : salt;
        const std::array<std::uint8_t, kSha256OutputSize> prk = hmac_sha256(extract_salt, ikm);

        std::vector<std::uint8_t> okm;
        okm.reserve(output_len);

        std::array<std::uint8_t, kSha256OutputSize> previous{};
        std::size_t previous_len = 0;
        std::uint8_t counter = 1;
        while (okm.size() < output_len) {
            std::vector<std::uint8_t> block_input;
            block_input.reserve(previous_len + info.size() + 1);
            if (previous_len != 0) {
                block_input.insert(block_input.end(), previous.begin(), previous.begin() + previous_len);
            }
            block_input.insert(block_input.end(), info.begin(), info.end());
            block_input.push_back(counter);

            previous = hmac_sha256(
                std::span<const std::uint8_t>(prk.data(), prk.size()),
                std::span<const std::uint8_t>(block_input.data(), block_input.size()));
            previous_len = previous.size();

            const std::size_t remaining = output_len - okm.size();
            const std::size_t take = std::min(remaining, previous.size());
            okm.insert(okm.end(), previous.begin(), previous.begin() + take);
            ++counter;
        }

        return okm;
    }
}
