/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "bindings/blst_scalar.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fr.h"

#include <array>
#include <cstring>
#include <limits>

namespace {

blst::Fp fp_from_u64(uint64_t value)
{
    blst::Fp raw = {};
    blst::Fp out;
    raw[0] = value;
    blst::Fp::to_mont(out, raw);
    return out;
}

} // namespace

int main()
{
    blst::Fr fr = blst::Fr::one();
    blst::Fr fr_out;
    blst::Fr::exp(fr_out, fr, nullptr, 0);
    if (!blst::Fr::equal(fr_out, blst::Fr::one()))
        return 1;
    blst::Fr::exp(fr, fr, nullptr, 0);
    if (!blst::Fr::equal(fr, blst::Fr::one()))
        return 2;

    const blst::Fp fp_a = fp_from_u64(2);
    const blst::Fp fp_b = fp_from_u64(3);
    blst::Fp fp_expected;
    blst::Fp fp_out;
    blst::Fp::mul(fp_expected, fp_a, fp_b);
    blst::Fp::sqr_n_mul(fp_out, fp_a, 0, fp_b);
    if (!blst::Fp::equal(fp_out, fp_expected))
        return 3;

    blst::Fp zero_count_alias = fp_b;
    blst::Fp::sqr_n_mul(zero_count_alias, fp_a, 0, zero_count_alias);
    if (!blst::Fp::equal(zero_count_alias, fp_expected))
        return 4;

    blst::Fp aliased = fp_b;
    fp_expected = fp_a;
    for (size_t i = 0; i < 3; ++i)
        blst::Fp::sqr(fp_expected, fp_expected);
    blst::Fp::mul(fp_expected, fp_expected, fp_b);
    blst::Fp::sqr_n_mul(aliased, fp_a, 3, aliased);
    if (!blst::Fp::equal(aliased, fp_expected))
        return 5;

    const std::array<uint32_t, 8> words32 = {1, 2, 3, 4, 5, 6, 7, 8};
    std::array<uint32_t, 8> roundtrip32 = {};
    blst_scalar scalar = {};
    blst_scalar_from_uint32(&scalar, words32.data());
    blst_uint32_from_scalar(roundtrip32.data(), &scalar);
    if (roundtrip32 != words32)
        return 6;

    const std::array<uint64_t, 4> words64 = {1, 2, 3, 4};
    std::array<uint64_t, 4> roundtrip64 = {};
    blst_scalar_from_uint64(&scalar, words64.data());
    blst_uint64_from_scalar(roundtrip64.data(), &scalar);
    if (roundtrip64 != words64)
        return 7;

    std::array<byte, 32> little = {};
    std::array<byte, 32> big = {};
    std::array<byte, 32> encoded = {};
    little.front() = 7;
    big.back() = 7;
    blst_scalar_from_lendian(&scalar, little.data());
    blst_lendian_from_scalar(encoded.data(), &scalar);
    if (encoded != little)
        return 8;
    blst_scalar_from_bendian(&scalar, big.data());
    blst_bendian_from_scalar(encoded.data(), &scalar);
    if (encoded != big)
        return 9;

    blst_fp nonresidue = {};
    bool found_nonresidue = false;
    for (uint64_t value = 2; value < 32; ++value) {
        const std::array<uint64_t, 6> words = {value, 0, 0, 0, 0, 0};
        blst_fp_from_uint64(&nonresidue, words.data());
        if (!blst_fp_is_square(&nonresidue)) {
            found_nonresidue = true;
            break;
        }
    }
    if (!found_nonresidue || blst_fp_sqrt(&nonresidue, &nonresidue))
        return 10;
    const blst_fp zero = {};
    if (std::memcmp(&nonresidue, &zero, sizeof(zero)) != 0)
        return 11;

    blst_scalar scalar_out = {};
    blst_fr field_in = {};
    const size_t invalid_count =
        std::numeric_limits<size_t>::max() / sizeof(blst_scalar) + 1;
    if (blst_scalars_from_fr(&scalar_out, &field_in, invalid_count))
        return 12;

    return 0;
}
