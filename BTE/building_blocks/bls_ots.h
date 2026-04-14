#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <string>

namespace bte::building_blocks
{
    struct BlsMinPkOneTimeSignature
    {
        std::array<std::uint8_t, 48> verification_key{};
        std::array<std::uint8_t, 96> signature{};
    };

    BlsMinPkOneTimeSignature sign_bls_min_pk(
        std::span<const std::uint8_t> ikm,
        std::span<const std::uint8_t> message,
        const std::string& dst);

    bool verify_bls_min_pk(
        std::span<const std::uint8_t> verification_key,
        std::span<const std::uint8_t> signature,
        std::span<const std::uint8_t> message,
        const std::string& dst);
}
