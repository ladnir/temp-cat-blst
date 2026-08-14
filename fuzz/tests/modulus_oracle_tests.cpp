// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/support/modulus_oracle.h"

#include <array>
#include <cstddef>

int main()
{
    using namespace catblst::fuzz;
    const ModulusOracle<8> scalar(scalar_modulus32);

    const std::array<std::byte, 1> zero{};
    const std::array<std::byte, 1> one{std::byte{1}};
    if (!ModulusOracle<8>::is_zero(scalar.reduce(zero, ByteOrder::little)))
        return 1;
    const auto reduced_one = scalar.reduce(one, ByteOrder::little);
    if (reduced_one[0] != 1 || ModulusOracle<8>::is_zero(reduced_one))
        return 2;

    const auto modulus_bytes = ModulusOracle<8>::to_le_bytes(scalar_modulus32);
    if (!ModulusOracle<8>::is_zero(
            scalar.reduce(modulus_bytes, ByteOrder::little)))
        return 3;

    auto modulus_minus_one = modulus_bytes;
    modulus_minus_one[0] = std::byte{0};
    const auto reduced_minus_one = scalar.reduce(
        modulus_minus_one,
        ByteOrder::little);
    if (ModulusOracle<8>::is_zero(reduced_minus_one) ||
        reduced_minus_one[0] != 0)
        return 4;

    std::array<std::byte, 32> all_ones{};
    all_ones.fill(std::byte{0xff});
    const auto little = scalar.reduce(all_ones, ByteOrder::little);
    const auto big = scalar.reduce(all_ones, ByteOrder::big);
    if (little != big)
        return 5;

    const ModulusOracle<12> fp(fp_modulus32);
    const auto fp_modulus_bytes = ModulusOracle<12>::to_le_bytes(fp_modulus32);
    if (!ModulusOracle<12>::is_zero(
            fp.reduce(fp_modulus_bytes, ByteOrder::little)))
        return 6;
    return 0;
}
