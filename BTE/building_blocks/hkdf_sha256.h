#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace bte::building_blocks
{
    std::array<std::uint8_t, 32> sha256(std::span<const std::uint8_t> message);

    std::vector<std::uint8_t> hkdf_sha256(
        std::span<const std::uint8_t> salt,
        std::span<const std::uint8_t> ikm,
        std::span<const std::uint8_t> info,
        std::size_t output_len);
}
