#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <string>

#include "blst/curve/point.h"

namespace bte::building_blocks
{
    std::array<std::uint8_t, 32> hash_to_scalar_le(
        std::span<const std::uint8_t> message,
        const std::string& dst);

    blst::G1Point hash_to_g1(
        std::span<const std::uint8_t> message,
        const std::string& dst);

    blst::G2Point hash_to_g2(
        std::span<const std::uint8_t> message,
        const std::string& dst);
}
