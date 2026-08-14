// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>

namespace catblst::fuzz
{
    template<std::size_t Capacity>
    struct DecodedSeed final
    {
        std::array<std::byte, Capacity> bytes{};
        std::size_t size = 0;

        [[nodiscard]] std::span<const std::byte> span() const noexcept
        {
            return std::span(bytes).first(size);
        }
    };

    [[nodiscard]] inline bool has_prefix(
        std::span<const std::byte> input,
        std::string_view prefix) noexcept
    {
        if (input.size() < prefix.size())
            return false;
        for (std::size_t index = 0; index < prefix.size(); ++index)
        {
            if (std::to_integer<unsigned char>(input[index]) !=
                static_cast<unsigned char>(prefix[index]))
                return false;
        }
        return true;
    }

    [[nodiscard]] inline std::span<const std::byte> text_seed_payload(
        std::span<const std::byte> input) noexcept
    {
        constexpr std::string_view prefix = "text:";
        if (!has_prefix(input, prefix))
            return {};
        input = input.subspan(prefix.size());
        while (!input.empty() &&
               (input.back() == std::byte{'\n'} || input.back() == std::byte{'\r'}))
            input = input.first(input.size() - 1);
        return input;
    }

    template<std::size_t Capacity>
    [[nodiscard]] std::optional<DecodedSeed<Capacity>> decode_hex_seed(
        std::span<const std::byte> input) noexcept
    {
        constexpr std::string_view prefix = "hex:";
        if (!has_prefix(input, prefix))
            return std::nullopt;
        input = input.subspan(prefix.size());
        while (!input.empty() &&
               (input.back() == std::byte{'\n'} || input.back() == std::byte{'\r'}))
            input = input.first(input.size() - 1);
        if (input.empty() || (input.size() & 1U) != 0 ||
            input.size() / 2 > Capacity)
            return std::nullopt;

        DecodedSeed<Capacity> result;
        result.size = input.size() / 2;
        for (std::size_t index = 0; index < result.size; ++index)
        {
            unsigned char value = 0;
            for (std::size_t half = 0; half < 2; ++half)
            {
                const unsigned char character =
                    std::to_integer<unsigned char>(input[index * 2 + half]);
                unsigned char nibble = 0;
                if (character >= '0' && character <= '9')
                    nibble = static_cast<unsigned char>(character - '0');
                else if (character >= 'a' && character <= 'f')
                    nibble = static_cast<unsigned char>(character - 'a' + 10);
                else if (character >= 'A' && character <= 'F')
                    nibble = static_cast<unsigned char>(character - 'A' + 10);
                else
                    return std::nullopt;
                value = static_cast<unsigned char>((value << 4) | nibble);
            }
            result.bytes[index] = std::byte{value};
        }
        return result;
    }

    template<std::size_t Capacity>
    [[nodiscard]] std::optional<DecodedSeed<Capacity>> decode_base64_seed(
        std::span<const std::byte> input) noexcept
    {
        constexpr std::string_view prefix = "b64:";
        if (!has_prefix(input, prefix))
            return std::nullopt;
        input = input.subspan(prefix.size());
        while (!input.empty() &&
               (input.back() == std::byte{'\n'} || input.back() == std::byte{'\r'}))
            input = input.first(input.size() - 1);
        if (input.empty() || (input.size() & 3U) != 0)
            return std::nullopt;

        const auto sextet = [](std::byte encoded) -> int {
            const unsigned char value = std::to_integer<unsigned char>(encoded);
            if (value >= 'A' && value <= 'Z')
                return value - 'A';
            if (value >= 'a' && value <= 'z')
                return value - 'a' + 26;
            if (value >= '0' && value <= '9')
                return value - '0' + 52;
            if (value == '+')
                return 62;
            if (value == '/')
                return 63;
            return -1;
        };

        DecodedSeed<Capacity> result;
        for (std::size_t offset = 0; offset < input.size(); offset += 4)
        {
            const bool last = offset + 4 == input.size();
            const bool pad2 = input[offset + 2] == std::byte{'='};
            const bool pad3 = input[offset + 3] == std::byte{'='};
            if ((!last && (pad2 || pad3)) || (pad2 && !pad3))
                return std::nullopt;
            const int a = sextet(input[offset]);
            const int b = sextet(input[offset + 1]);
            const int c = pad2 ? 0 : sextet(input[offset + 2]);
            const int d = pad3 ? 0 : sextet(input[offset + 3]);
            if (a < 0 || b < 0 || c < 0 || d < 0)
                return std::nullopt;
            const std::uint32_t block =
                (static_cast<std::uint32_t>(a) << 18) |
                (static_cast<std::uint32_t>(b) << 12) |
                (static_cast<std::uint32_t>(c) << 6) |
                static_cast<std::uint32_t>(d);
            const std::size_t emitted = pad2 ? 1 : (pad3 ? 2 : 3);
            if (emitted > Capacity || result.size > Capacity - emitted)
                return std::nullopt;
            result.bytes[result.size++] = std::byte{
                static_cast<unsigned char>(block >> 16)};
            if (emitted >= 2)
                result.bytes[result.size++] = std::byte{
                    static_cast<unsigned char>(block >> 8)};
            if (emitted == 3)
                result.bytes[result.size++] = std::byte{
                    static_cast<unsigned char>(block)};
        }
        return result;
    }
}
