/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "bindings/blst_scalar.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp12.h"
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

    alignas(uint64_t) blst_scalar aliased_scalar = {};
    std::memcpy(&aliased_scalar, words64.data(), sizeof(aliased_scalar));
    blst_scalar_from_uint64(
        &aliased_scalar,
        reinterpret_cast<const uint64_t*>(&aliased_scalar));
    if (std::memcmp(&aliased_scalar, words64.data(), sizeof(aliased_scalar)) != 0)
        return 13;

    std::array<uint64_t, 4> aliased_words = words64;
    blst_uint64_from_scalar(
        aliased_words.data(),
        reinterpret_cast<const blst_scalar*>(aliased_words.data()));
    if (aliased_words != words64)
        return 14;

    blst_fr expected_fr;
    blst_fr_from_scalar(&expected_fr, &aliased_scalar);
    blst_fr aliased_fr;
    std::memcpy(&aliased_fr, &aliased_scalar, sizeof(aliased_fr));
    blst_fr_from_scalar(
        &aliased_fr,
        reinterpret_cast<const blst_scalar*>(&aliased_fr));
    if (std::memcmp(&aliased_fr, &expected_fr, sizeof(aliased_fr)) != 0)
        return 15;

    blst_scalar expected_scalar;
    blst_scalar_from_fr(&expected_scalar, &expected_fr);
    blst_scalar scalar_from_aliased_fr;
    std::memcpy(&scalar_from_aliased_fr, &expected_fr, sizeof(scalar_from_aliased_fr));
    blst_scalar_from_fr(
        &scalar_from_aliased_fr,
        reinterpret_cast<const blst_fr*>(&scalar_from_aliased_fr));
    if (std::memcmp(&scalar_from_aliased_fr, &expected_scalar,
            sizeof(scalar_from_aliased_fr)) != 0)
        return 16;

    const blst::Fp2& fp2_one = blst::Fp2::one();
    const blst::Fp12& fp12_one = blst::Fp12::one();
    if (!blst::Fp::equal(blst::Fp::one(), fp2_one[0]) ||
        !blst::Fp::equal(blst::Fp::zero(), fp2_one[1]) ||
        !blst::Fp::equal(blst::Fp::one(), fp12_one[0][0][0]) ||
        !blst::Fp::equal(blst::Fp::zero(), fp12_one[0][0][1]))
        return 17;

    blst::Fp fp_inverse;
    blst::Fp fp_inverse_check;
    blst::Fp::reciprocal(fp_inverse, fp_a);
    blst::Fp::mul(fp_inverse_check, fp_inverse, fp_a);
    if (!blst::Fp::equal(fp_inverse_check, blst::Fp::one()))
        return 18;

    std::array<blst_byte, 32> little = {};
    std::array<blst_byte, 32> big = {};
    std::array<blst_byte, 32> encoded = {};
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
