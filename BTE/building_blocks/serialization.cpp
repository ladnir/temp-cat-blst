#include "BTE/building_blocks/serialization.h"

namespace bte::building_blocks
{
    void append_u8(Bytes& out, const std::uint8_t value)
    {
        out.push_back(value);
    }

    void append_u32_be(Bytes& out, const std::uint32_t value)
    {
        out.push_back(static_cast<std::uint8_t>((value >> 24) & 0xffU));
        out.push_back(static_cast<std::uint8_t>((value >> 16) & 0xffU));
        out.push_back(static_cast<std::uint8_t>((value >> 8) & 0xffU));
        out.push_back(static_cast<std::uint8_t>(value & 0xffU));
    }

    void append_string(Bytes& out, const std::string_view value)
    {
        append_bytes(out, std::span<const std::uint8_t>(
            reinterpret_cast<const std::uint8_t*>(value.data()),
            value.size()));
    }

    void append_bytes(Bytes& out, const std::span<const std::uint8_t> value)
    {
        out.insert(out.end(), value.begin(), value.end());
    }
}
