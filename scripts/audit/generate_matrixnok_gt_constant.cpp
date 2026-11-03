// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

#include "blst/avx/fields/conversion_inversion.h"
#include "blst/avx/fields/matrixnok_fp12.h"
#include "blst/scalar/curve/GT.h"

#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>

int main()
{
    using Ring = blst::avx::matrixnok::FpRingCandidate;
    using Fp12 = blst::avx::matrixnok::Fp12RingCandidate;

    const Ring ring;
    const Fp12 fp12ring(ring);
    const auto value = blst::avx::avx_blst_direct_bridge::scalar_fp12_to_ring(
        ring, blst::scalar::GT::generator().data);

    std::array<unsigned char, 12 * 48> serialized = {};
    blst::scalar::GT::serialize(
        serialized.data(), blst::scalar::GT::generator());
    std::cout << "// canonical big-endian Fp12: ";
    std::cout << std::hex << std::setfill('0');
    for (const unsigned char byte : serialized) {
        std::cout << std::setw(2) << static_cast<unsigned>(byte);
    }
    std::cout << std::dec << "\n";

    std::cout << "inline constexpr Fp12Constant gt_generator = {{\n";
    for (const auto& coefficient : value) {
        const auto m1 = coefficient.m1.to_unsigned_array();
        const auto m2 = coefficient.m2.to_unsigned_array();
        std::cout << "    {{{";
        for (std::size_t i = 0; i < m1.size(); ++i) {
            std::cout << (i == 0 ? "" : ", ") << m1[i] << "ULL";
        }
        std::cout << "}, {";
        for (std::size_t i = 0; i < m2.size(); ++i) {
            std::cout << (i == 0 ? "" : ", ") << m2[i] << "ULL";
        }
        std::cout << "}}},\n";
    }
    std::cout << "}};\n";

    (void)fp12ring;
    return 0;
}
