#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

namespace bte::building_blocks
{
    using Bytes = std::vector<std::uint8_t>;

    void append_u8(Bytes& out, std::uint8_t value);
    void append_u32_be(Bytes& out, std::uint32_t value);
    void append_string(Bytes& out, std::string_view value);
    void append_bytes(Bytes& out, std::span<const std::uint8_t> value);

    template<std::size_t N>
    void append_array(Bytes& out, const std::array<std::uint8_t, N>& value)
    {
        append_bytes(out, std::span<const std::uint8_t>(value.data(), value.size()));
    }
}
