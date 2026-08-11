/*
 * Copyright 2026 Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 *
 * Locally composed around the fixed RNS conversion used by the VROOM AVX
 * backend and CatBLST's scalar field operations.
 */
#pragma once

#include "fp12.h"
#include "../cpu/vector/vector_impl.hpp"

#include "blst/fields/Fp.h"
#include "blst/fields/Fp12.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace blst::avx {

namespace avx_blst_direct_bridge {

BLST_FORCEINLINE void store_u64_unaligned(unsigned char* dst, uint64_t value)
{
    std::memcpy(dst, &value, sizeof(value));
}

// Export the raw radix-2^52 accumulators produced by
// MatrixNoKFp::convert_from_raw. The selected conversion constants guarantee
// that the adjacent digit sums and the final 576-bit addition do not wrap,
// and that the reconstructed integer satisfies Fp::redc's input bound.
// The result is returned in CatBLST's ordinary Montgomery Fp representation.
BLST_FORCEINLINE void from_raw_accumulators(
    const AVXVector<8>& vec_lo,
    const AVXVector<8>& vec_hi,
    blst::Fp& ret)
{
    std::array<uint64_t, 8> digits_lo = {};
    std::array<uint64_t, 8> digits_hi = {};
    vec_hi.store(digits_hi.data());
    vec_lo.store(digits_lo.data());

    std::array<uint64_t, 9> digits = {};
    digits[0] = digits_lo[0];
    for (int i = 1; i < 8; ++i) {
#ifndef NDEBUG
        BLST_ASSERT(digits_hi[i - 1] <=
            std::numeric_limits<uint64_t>::max() - digits_lo[i]);
#endif
        digits[i] = digits_hi[i - 1] + digits_lo[i];
    }
    digits[8] = digits_hi[7];

    const blst::Fp& kRR = blst::get_BLST_RR();
    blst::UFp n_data_even_64 = {};
    blst::UFp n_data_odd_64 = {};
    unsigned char* n_data_even = reinterpret_cast<unsigned char*>(n_data_even_64.data);
    unsigned char* n_data_odd = reinterpret_cast<unsigned char*>(n_data_odd_64.data);

    for (int i = 0; i < 9; i += 2) {
        store_u64_unaligned(n_data_even + (i / 2) * 13, digits[i]);
    }
    for (int i = 1; i < 8; i += 2) {
        store_u64_unaligned(n_data_odd + (i / 2) * 13 + 6, digits[i] << 4);
    }

    blst::UFp sum = {};
    uint64_t carry = 0;
    for (int i = 0; i < 9; ++i) {
        blst::wide_limb_t wide = {};
        blst::cc_add_carry_wide_limb_limb(
            &wide,
            n_data_even_64.data[i],
            n_data_odd_64.data[i],
            carry);
        sum.data[i] = static_cast<uint64_t>(wide.lo);
        carry = static_cast<uint64_t>(wide.hi);
    }
    BLST_ASSERT(carry == 0);

    blst::Fp temp = {};
    blst::Fp::redc(temp, sum);
    blst::Fp::mul(ret, temp, kRR);

    blst::Fp mont = {};
    blst::Fp::to_mont(mont, ret);
    ret = mont;
}

template<class Ring>
BLST_FORCEINLINE void convert_to_scalar_fp(
    const Ring&,
    const typename Ring::StandardElement& value,
    blst::Fp& ret)
{
    AVXVector<8> out_hi = AVXVector<8>(0);
    AVXVector<8> out_lo = AVXVector<8>(0);
    Ring::convert_from_raw(value, out_hi, out_lo);
    from_raw_accumulators(out_lo, out_hi, ret);
}

// Import a valid internal Montgomery Fp value into the selected MatrixNoK
// Ring's standard representation. This performs no allocation or validation.
template<class Ring>
BLST_FORCEINLINE auto convert_to_rns(const Ring&, const blst::Fp& value)
{
    blst::Fp canonical = {};
    blst::Fp::from_mont(canonical, value);
    return Ring::from_canonical_limbs(canonical);
}

template<class Ring>
BLST_FORCEINLINE blst::Fp ring_to_scalar_fp(
    const Ring& ring,
    const typename Ring::StandardElement& value)
{
    blst::Fp mont = {};
    convert_to_scalar_fp(ring, value, mont);
    return mont;
}

template<class Ring>
BLST_FORCEINLINE typename Ring::StandardElement scalar_fp_to_ring(
    const Ring& ring,
    const blst::Fp& value)
{
    return convert_to_rns<Ring>(ring, value);
}

template<class Ring>
blst::Fp12 ring_fp12_to_scalar(
    const Ring& ring,
    const std::array<typename Ring::StandardElement, 12>& value)
{
    blst::Fp12 out = {};
    out[0][0][0] = ring_to_scalar_fp(ring, value[0]);
    out[0][0][1] = ring_to_scalar_fp(ring, value[1]);
    out[0][1][0] = ring_to_scalar_fp(ring, value[2]);
    out[0][1][1] = ring_to_scalar_fp(ring, value[3]);
    out[0][2][0] = ring_to_scalar_fp(ring, value[4]);
    out[0][2][1] = ring_to_scalar_fp(ring, value[5]);
    out[1][0][0] = ring_to_scalar_fp(ring, value[6]);
    out[1][0][1] = ring_to_scalar_fp(ring, value[7]);
    out[1][1][0] = ring_to_scalar_fp(ring, value[8]);
    out[1][1][1] = ring_to_scalar_fp(ring, value[9]);
    out[1][2][0] = ring_to_scalar_fp(ring, value[10]);
    out[1][2][1] = ring_to_scalar_fp(ring, value[11]);
    return out;
}

template<class Ring>
std::array<typename Ring::StandardElement, 12> scalar_fp12_to_ring(
    const Ring& ring,
    const blst::Fp12& value)
{
    return {
        scalar_fp_to_ring(ring, value[0][0][0]),
        scalar_fp_to_ring(ring, value[0][0][1]),
        scalar_fp_to_ring(ring, value[0][1][0]),
        scalar_fp_to_ring(ring, value[0][1][1]),
        scalar_fp_to_ring(ring, value[0][2][0]),
        scalar_fp_to_ring(ring, value[0][2][1]),
        scalar_fp_to_ring(ring, value[1][0][0]),
        scalar_fp_to_ring(ring, value[1][0][1]),
        scalar_fp_to_ring(ring, value[1][1][0]),
        scalar_fp_to_ring(ring, value[1][1][1]),
        scalar_fp_to_ring(ring, value[1][2][0]),
        scalar_fp_to_ring(ring, value[1][2][1]),
    };
}

}  // namespace avx_blst_direct_bridge

/**
 * Invert one nonzero selected-RNS Fp element through scalar CatBLST.
 *
 * The operation is allocation-free and nonthrowing, but variable-time because
 * Fp::reciprocal may take its fixed-chain fallback. The result remains in the
 * selected Ring's standard representation. Zero is outside the contract.
 */
template<class Ring>
auto invert_fp_vartime_via_scalar(
    const Ring& ring,
    const typename Ring::StandardElement& value)
{
    blst::Fp scalar = avx_blst_direct_bridge::ring_to_scalar_fp(ring, value);
    BLST_ASSERT(!blst::Fp::equal(scalar, blst::Fp::zero()));
    blst::Fp inverse = {};
    blst::Fp::reciprocal(inverse, scalar);
    return avx_blst_direct_bridge::scalar_fp_to_ring(ring, inverse);
}

template<class Ring>
class ScalarFpVartimeInverter {
public:
    using StandardElement = typename Ring::StandardElement;

    /** Invert a nonzero element; allocation-free, nonthrowing, variable-time. */
    StandardElement invert(const StandardElement& a, const Ring& ring) const
    {
        return invert_fp_vartime_via_scalar(ring, a);
    }
};

/**
 * Decomposed inversion provider for a nonzero selected-RNS Fp12 element.
 *
 * This is allocation-free and nonthrowing. Its single base-field inversion is
 * variable-time; callers must not use it for secret-dependent Fp12 values.
 */
template<class Ring, class FP12R>
class FP12RNSVartimeInverter {
public:
    using StandardElement = typename FP12R::StandardElement;
    using Inverter = ScalarFpVartimeInverter<Ring>;

    Inverter inverter;

    StandardElement invert(const StandardElement& f, const FP12R& fp12ring, const Ring& ring) const
    {
        return fp12ring.template fp12_flat_inverse<Inverter>(f, inverter, ring);
    }
};

}  // namespace blst::avx
