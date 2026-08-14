// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "fuzz/internal/factory_access.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>

namespace catblst::fuzz::internal
{
    struct ScalarStorage final
    {
        std::array<std::byte, 32> bytes{};
        std::size_t nbits = 0;
    };

    using ScalarBitsValue = ScalarBits<ScalarStorage>;
    using CanonicalScalarValue = CanonicalScalar<ScalarStorage>;
    using SecretScalarValue = SecretScalar<ScalarStorage>;

    class ScalarFactory final
    {
      public:
        [[nodiscard]] static ScalarBitsValue scalar_bits(
            std::span<const std::byte> input,
            std::size_t requested_bits) noexcept
        {
            ScalarStorage storage;
            storage.nbits = std::min<std::size_t>(requested_bits, 256);
            const std::size_t byte_count = (storage.nbits + 7) / 8;
            const std::size_t copied = std::min(input.size(), byte_count);
            std::copy_n(input.begin(), copied, storage.bytes.begin());

            const std::size_t partial_bits = storage.nbits & 7U;
            if (partial_bits != 0 && byte_count != 0)
            {
                const auto mask = static_cast<unsigned char>(
                    (1U << partial_bits) - 1U);
                storage.bytes[byte_count - 1] &= std::byte{mask};
            }
            return FactoryAccess::make<ScalarStorage, ScalarBitsTag>(
                std::move(storage));
        }

        [[nodiscard]] static std::optional<CanonicalScalarValue>
        canonical_scalar(const ScalarBitsValue& input) noexcept
        {
            const ScalarStorage& storage = FactoryAccess::storage(input);
            if (!less_than_modulus(storage.bytes))
                return std::nullopt;
            return FactoryAccess::make<ScalarStorage, CanonicalScalarTag>(storage);
        }

        [[nodiscard]] static std::optional<SecretScalarValue> secret_scalar(
            const CanonicalScalarValue& input) noexcept
        {
            const ScalarStorage& storage = FactoryAccess::storage(input);
            unsigned char aggregate = 0;
            for (const std::byte value : storage.bytes)
                aggregate |= std::to_integer<unsigned char>(value);
            if (aggregate == 0)
                return std::nullopt;
            return FactoryAccess::make<ScalarStorage, SecretScalarTag>(storage);
        }

      private:
        static constexpr std::array<std::uint64_t, 4> scalar_modulus = {
            0xffffffff00000001ULL,
            0x53bda402fffe5bfeULL,
            0x3339d80809a1d805ULL,
            0x73eda753299d7d48ULL,
        };

        [[nodiscard]] static bool less_than_modulus(
            const std::array<std::byte, 32>& bytes) noexcept
        {
            for (std::size_t limb_index = 4; limb_index-- != 0;)
            {
                std::uint64_t limb = 0;
                for (std::size_t byte_index = 0; byte_index < 8; ++byte_index)
                {
                    limb |= static_cast<std::uint64_t>(
                                std::to_integer<unsigned char>(
                                    bytes[limb_index * 8 + byte_index]))
                            << (byte_index * 8);
                }
                if (limb < scalar_modulus[limb_index])
                    return true;
                if (limb > scalar_modulus[limb_index])
                    return false;
            }
            return false;
        }
    };
}
