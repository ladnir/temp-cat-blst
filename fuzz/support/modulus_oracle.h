// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

namespace catblst::fuzz
{
    enum class ByteOrder : std::uint8_t
    {
        little,
        big,
    };

    template<std::size_t Limbs>
    class ModulusOracle final
    {
      public:
        using Value = std::array<std::uint32_t, Limbs>;

        explicit constexpr ModulusOracle(Value modulus) noexcept
            : modulus_(modulus)
        {}

        [[nodiscard]] Value reduce(
            std::span<const std::byte> input,
            ByteOrder order) const noexcept
        {
            Value remainder{};
            if (order == ByteOrder::big)
            {
                for (const std::byte byte : input)
                    consume_byte(remainder, byte);
            }
            else
            {
                for (std::size_t index = input.size(); index-- != 0;)
                    consume_byte(remainder, input[index]);
            }
            return remainder;
        }

        [[nodiscard]] static bool is_zero(const Value& value) noexcept
        {
            std::uint32_t aggregate = 0;
            for (const std::uint32_t limb : value)
                aggregate |= limb;
            return aggregate == 0;
        }

        [[nodiscard]] static std::array<std::byte, Limbs * 4> to_le_bytes(
            const Value& value) noexcept
        {
            std::array<std::byte, Limbs * 4> output{};
            for (std::size_t limb = 0; limb < Limbs; ++limb)
            {
                for (std::size_t byte = 0; byte < 4; ++byte)
                {
                    output[limb * 4 + byte] = std::byte(
                        static_cast<unsigned char>(value[limb] >> (byte * 8)));
                }
            }
            return output;
        }

      private:
        void consume_byte(Value& value, std::byte byte) const noexcept
        {
            const unsigned char input = std::to_integer<unsigned char>(byte);
            for (unsigned bit = 8; bit-- != 0;)
            {
                std::uint64_t carry = (input >> bit) & 1U;
                for (std::size_t limb = 0; limb < Limbs; ++limb)
                {
                    const std::uint64_t shifted =
                        (static_cast<std::uint64_t>(value[limb]) << 1) | carry;
                    value[limb] = static_cast<std::uint32_t>(shifted);
                    carry = shifted >> 32;
                }
                if (carry != 0 || !less_than(value, modulus_))
                    subtract(value, modulus_);
            }
        }

        [[nodiscard]] static bool less_than(
            const Value& left,
            const Value& right) noexcept
        {
            for (std::size_t index = Limbs; index-- != 0;)
            {
                if (left[index] != right[index])
                    return left[index] < right[index];
            }
            return false;
        }

        static void subtract(Value& left, const Value& right) noexcept
        {
            std::uint64_t borrow = 0;
            for (std::size_t index = 0; index < Limbs; ++index)
            {
                const std::uint64_t subtrahend =
                    static_cast<std::uint64_t>(right[index]) + borrow;
                const std::uint64_t minuend = left[index];
                left[index] = static_cast<std::uint32_t>(minuend - subtrahend);
                borrow = minuend < subtrahend;
            }
        }

        Value modulus_;
    };

    inline constexpr ModulusOracle<8>::Value scalar_modulus32 = {
        0x00000001U,
        0xffffffffU,
        0xfffe5bfeU,
        0x53bda402U,
        0x09a1d805U,
        0x3339d808U,
        0x299d7d48U,
        0x73eda753U,
    };

    inline constexpr ModulusOracle<12>::Value fp_modulus32 = {
        0xffffaaabU,
        0xb9feffffU,
        0xb153ffffU,
        0x1eabfffeU,
        0xf6b0f624U,
        0x6730d2a0U,
        0xf38512bfU,
        0x64774b84U,
        0x434bacd7U,
        0x4b1ba7b6U,
        0x397fe69aU,
        0x1a0111eaU,
    };
}

