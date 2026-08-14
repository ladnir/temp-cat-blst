// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/scalar_factory.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace
{
    std::array<std::byte, 32> little_endian(
        const std::array<std::uint64_t, 4>& limbs)
    {
        std::array<std::byte, 32> result{};
        for (std::size_t i = 0; i < limbs.size(); ++i)
            for (std::size_t j = 0; j < 8; ++j)
                result[i * 8 + j] = std::byte(
                    static_cast<unsigned char>(limbs[i] >> (j * 8)));
        return result;
    }
}

int main()
{
    using namespace catblst::fuzz::internal;
    constexpr std::array<std::uint64_t, 4> modulus = {
        0xffffffff00000001ULL,
        0x53bda402fffe5bfeULL,
        0x3339d80809a1d805ULL,
        0x73eda753299d7d48ULL,
    };

    const std::array<std::byte, 32> zero{};
    const auto zero_bits = ScalarFactory::scalar_bits(zero, 256);
    const auto canonical_zero = ScalarFactory::canonical_scalar(zero_bits);
    if (!canonical_zero || ScalarFactory::secret_scalar(*canonical_zero))
        return 1;

    auto one = zero;
    one[0] = std::byte{1};
    const auto one_bits = ScalarFactory::scalar_bits(one, 256);
    const auto canonical_one = ScalarFactory::canonical_scalar(one_bits);
    if (!canonical_one || !ScalarFactory::secret_scalar(*canonical_one))
        return 2;

    auto modulus_minus_one = modulus;
    --modulus_minus_one[0];
    if (!ScalarFactory::canonical_scalar(
            ScalarFactory::scalar_bits(little_endian(modulus_minus_one), 256)))
        return 3;
    if (ScalarFactory::canonical_scalar(
            ScalarFactory::scalar_bits(little_endian(modulus), 256)))
        return 4;

    std::array<std::byte, 32> all_ones{};
    all_ones.fill(std::byte{0xff});
    const auto masked = ScalarFactory::scalar_bits(all_ones, 13);
    const ScalarStorage& storage = FactoryAccess::storage(masked);
    if (storage.nbits != 13 || storage.bytes[0] != std::byte{0xff} ||
        storage.bytes[1] != std::byte{0x1f} || storage.bytes[2] != std::byte{0})
        return 5;
    return 0;
}
