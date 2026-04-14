#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "blst/fields/Fr.h"

namespace bte::building_blocks
{
    std::vector<blst::Fr> lagrange_at_zero(std::span<const std::uint32_t> indices);
}
