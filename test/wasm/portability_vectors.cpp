// Copyright Category Labs
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
#include "bindings/blst.h"
#include "bindings/blst_fft.h"
#include <array>
#include <cstdio>
#include <vector>

namespace {
void emit(const unsigned char* bytes, size_t length)
{
    for (size_t i = 0; i < length; ++i)
        std::printf("%02x", static_cast<unsigned>(bytes[i]));
    std::putchar('\n');
}
void emit_fr(const blst_fr& value)
{
    blst_scalar scalar;
    blst_scalar_from_fr(&scalar, &value);
    emit(scalar.b, sizeof(scalar.b));
}
}

// One byte-for-byte output shared by native and Wasm runs. The ordinary unit
// suite supplies independent known-answer vectors and algorithmic oracles.
int main()
{
    std::array<unsigned char, 257> message;
    for (size_t i = 0; i < message.size(); ++i)
        message[i] = static_cast<unsigned char>(i);
    for (size_t length : {size_t{0}, size_t{1}, size_t{55}, size_t{56},
                          size_t{63}, size_t{64}, size_t{65}, size_t{257}}) {
        unsigned char digest[32];
        blst_sha256(digest, message.data(), length);
        emit(digest, sizeof(digest));
    }
    const unsigned char dst[] = "CatBLST-Wasm-portability-v1";
    blst_p1 p1;
    blst_p2 p2;
    blst_hash_to_g1(&p1, message.data(), message.size(), dst, sizeof(dst)-1, nullptr, 0);
    blst_hash_to_g2(&p2, message.data(), message.size(), dst, sizeof(dst)-1, nullptr, 0);
    unsigned char bytes[576];
    blst_p1_compress(bytes, &p1); emit(bytes, 48);
    blst_p2_compress(bytes, &p2); emit(bytes, 96);
    blst_p1_affine a1;
    blst_p2_affine a2;
    blst_p1_to_affine(&a1, &p1);
    blst_p2_to_affine(&a2, &p2);
    if (!blst_p1_affine_in_g1(&a1) || !blst_p2_affine_in_g2(&a2)) return 1;
    blst_fp12 miller, pairing;
    blst_miller_loop(&miller, &a2, &a1);
    blst_final_exp(&pairing, &miller);
    blst_bendian_from_fp12(bytes, &pairing); emit(bytes, sizeof(bytes));
    std::array<blst_fr, 32> values;
    std::array<blst_p1_affine, 32> points;
    std::uint64_t seed = 0x123456789abcdef0ULL;
    for (size_t i = 0; i < values.size(); ++i) {
        seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
        const std::uint64_t limbs[4] = {seed, 0, 0, 0};
        blst_fr_from_uint64(&values[i], limbs);
        blst_fr product, inverse;
        blst_fr_mul(&product, &values[i], &values[i]);
        blst_fr_inverse(&inverse, &values[i]);
        emit_fr(product); emit_fr(inverse);
        blst_scalar scalar;
        blst_scalar_from_fr(&scalar, &values[i]);
        blst_p1 point;
        blst_p1_mult(&point, blst_p1_generator(), scalar.b, 256);
        blst_p1_to_affine(&points[i], &point);
    }
    const size_t scratch_size = blst_p1s_mult_pippenger_fr_scratch_sizeof(points.size());
    if (scratch_size == 0) return 2;
    std::vector<blst_limb_t> scratch((scratch_size + sizeof(blst_limb_t)-1)/sizeof(blst_limb_t));
    blst_p1 msm;
    if (!blst_p1s_mult_pippenger_fr(&msm, points.data(), values.data(), points.size(), scratch.data(), scratch_size)) return 3;
    blst_p1_compress(bytes, &msm); emit(bytes, 48);
    const auto original = values;
    const size_t domain_length = blst_fr_fft_domain_len(values.size());
    std::vector<blst_fr> domain(domain_length);
    if (!blst_fr_fft_domain_init(domain.data(), domain.size(), values.size())) return 4;
    if (!blst_fr_fft_with_domain(values.data(), values.size(), domain.data(), domain.size(), false)) return 5;
    for (const auto& value : values) emit_fr(value);
    if (!blst_fr_fft_with_domain(values.data(), values.size(), domain.data(), domain.size(), true)) return 6;
    for (size_t i = 0; i < values.size(); ++i) {
        blst_scalar actual, expected;
        blst_scalar_from_fr(&actual, &values[i]); blst_scalar_from_fr(&expected, &original[i]);
        for (size_t j = 0; j < sizeof(actual.b); ++j)
            if (actual.b[j] != expected.b[j]) return 7;
    }
    return 0;
}
