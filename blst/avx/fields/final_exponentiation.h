// SPDX-License-Identifier: MIT
// Copyright 2026 Simon Langowski
// Copyright 2026 Category Labs
//
// Adapted from the VROOM final-exponentiation schedule at commit 26bd9aee and
// the corresponding BLST scalar construction. CatBLST retains BLST's cubed
// pairing exponent while making ownership, domain, and timing explicit.
#pragma once

#include "fp12.h"
#include <cstddef>

#if defined(_MSC_VER)
#define BLST_AVX_FINAL_EXP_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define BLST_AVX_FINAL_EXP_NOINLINE __attribute__((noinline))
#else
#define BLST_AVX_FINAL_EXP_NOINLINE
#endif

namespace blst::avx {
namespace final_exp_detail {

// The split square reads every input coordinate in both halves. This local
// boundary therefore requires a distinct output owner until both halves have
// completed; exact in-place evaluation would make the second half observe
// partially overwritten input.
template<class Output, class Input, class Fp12RingType, class Ring>
BLST_FORCEINLINE void cyclotomic_square_to(
    Output& out,
    const Input& in,
    const Fp12RingType& fp12ring,
    const Ring& ring)
{
    BLST_ASSERT(static_cast<const void*>(&out) !=
        static_cast<const void*>(&in));
    fp12ring.fp12_flat_cyclotomic_sqr_first_split(out, in, ring);
    fp12ring.fp12_flat_cyclotomic_sqr_second_split(out, in, ring);
}

template<class Element, class Fp12RingType, class Ring>
BLST_FORCEINLINE void conjugate_in_place(
    Element& value,
    const Fp12RingType&,
    const Ring& ring)
{
    value[6] = ring.standard_negate(value[6]);
    value[7] = ring.standard_negate(value[7]);
    value[8] = ring.standard_negate(value[8]);
    value[9] = ring.standard_negate(value[9]);
    value[10] = ring.standard_negate(value[10]);
    value[11] = ring.standard_negate(value[11]);
}

// Multiply current by the fixed addition-chain base, then perform exactly
// Squarings cyclotomic squares. The product is placed according to parity so
// the final value always lands in out without a full-element finishing copy.
template<std::size_t Squarings, class Current, class Base,
    class Fp12RingType, class Ring>
BLST_AVX_FINAL_EXP_NOINLINE void mul_then_squares_to(
    typename Fp12RingType::template CommonElement<6>& out,
    typename Fp12RingType::template CommonElement<6>& scratch,
    const Current& current,
    const Base& base,
    const Fp12RingType& fp12ring,
    const Ring& ring)
{
    static_assert(Squarings > 0);
    BLST_ASSERT(&out != &scratch);

    if constexpr ((Squarings & 1U) == 0U) {
        out = fp12ring.recast_common(
            fp12ring.fp12_flat_mul(current, base, ring));
        for (std::size_t i = 0; i < Squarings; i += 2) {
            cyclotomic_square_to(scratch, out, fp12ring, ring);
            cyclotomic_square_to(out, scratch, fp12ring, ring);
        }
    } else {
        scratch = fp12ring.recast_common(
            fp12ring.fp12_flat_mul(current, base, ring));
        cyclotomic_square_to(out, scratch, fp12ring, ring);
        for (std::size_t i = 1; i < Squarings; i += 2) {
            cyclotomic_square_to(scratch, out, fp12ring, ring);
            cyclotomic_square_to(out, scratch, fp12ring, ring);
        }
    }
}

// Computes a^(-|z|/2). The input must be cyclotomic and must not overlap either
// work owner. Release performs no membership or overlap validation.
template<class Element, class Fp12RingType, class Ring>
BLST_AVX_FINAL_EXP_NOINLINE void raise_to_z_div_by_2_to(
    typename Fp12RingType::template CommonElement<6>& out,
    typename Fp12RingType::template CommonElement<6>& scratch,
    const Element& a,
    const Fp12RingType& fp12ring,
    const Ring& ring)
{
    BLST_ASSERT(&out != &scratch);
    BLST_ASSERT(static_cast<const void*>(&out) !=
        static_cast<const void*>(&a));
    BLST_ASSERT(static_cast<const void*>(&scratch) !=
        static_cast<const void*>(&a));

    cyclotomic_square_to(out, a, fp12ring, ring);                   // 0x2
    mul_then_squares_to<2>(out, scratch, out, a, fp12ring, ring);   // 0xc
    mul_then_squares_to<3>(out, scratch, out, a, fp12ring, ring);   // 0x68
    mul_then_squares_to<9>(out, scratch, out, a, fp12ring, ring);   // 0xd200
    mul_then_squares_to<32>(out, scratch, out, a, fp12ring, ring);  // 0xd20100000000
    mul_then_squares_to<15>(out, scratch, out, a, fp12ring, ring);  // 0x6900800000008000
    conjugate_in_place(out, fp12ring, ring); // z is negative
}

// Computes a^(-|z|). Reversing the two work owners for the half exponent makes
// the final square land directly in out.
template<class Element, class Fp12RingType, class Ring>
BLST_AVX_FINAL_EXP_NOINLINE void raise_to_z_to(
    typename Fp12RingType::template CommonElement<6>& out,
    typename Fp12RingType::template CommonElement<6>& scratch,
    const Element& a,
    const Fp12RingType& fp12ring,
    const Ring& ring)
{
    raise_to_z_div_by_2_to(scratch, out, a, fp12ring, ring);
    cyclotomic_square_to(out, scratch, fp12ring, ring);
}

} // namespace final_exp_detail

// These returning wrappers are retained for the native GT facade. Production
// final exponentiation uses the distinct-output helpers above directly.
template<class Element, class Fp12RingType, class Ring>
inline auto raise_to_z_div_by_2(
    const Element& a, const Fp12RingType& fp12ring, const Ring& ring)
{
    typename Fp12RingType::template CommonElement<6> out;
    typename Fp12RingType::template CommonElement<6> scratch;
    final_exp_detail::raise_to_z_div_by_2_to(
        out, scratch, a, fp12ring, ring);
    return out;
}

template<class Element, class Fp12RingType, class Ring>
inline auto raise_to_z(
    const Element& a, const Fp12RingType& fp12ring, const Ring& ring)
{
    typename Fp12RingType::template CommonElement<6> out;
    typename Fp12RingType::template CommonElement<6> scratch;
    final_exp_detail::raise_to_z_to(out, scratch, a, fp12ring, ring);
    return out;
}

// Requires nonzero f and its exact multiplicative inverse. Both are valid Fp12
// values. The result is f^(3*((p^12-1)/r)) and therefore lies in GT. The fixed
// hard part is allocation-free and has no data-dependent control flow.
template<class Fp12RingType, class Ring>
inline auto final_exp_after_inverse(
    const typename Fp12RingType::StandardElement& f,
    const typename Fp12RingType::StandardElement& f_inverse,
    const Fp12RingType& fp12ring,
    const Ring& ring)
{
    using Standard = typename Fp12RingType::StandardElement;
    using Wide = typename Fp12RingType::template CommonElement<6>;

    // Easy part. Keep one standard temporary and the persistent ret owner.
    Standard ret = fp12ring.fp12_flat_conjugate(f, ring);
    ret = fp12ring.fp12_flat_mul(ret, f_inverse, ring);
    Standard y1 = fp12ring.fp12_flat_frobenius_map2(ret, ring);
    ret = fp12ring.fp12_flat_mul(ret, y1, ring);

    // Hard part. y0 is retained for its final product. y1w/y2w/y3w are the
    // bounded work owners needed by two consecutive fixed-z raises.
    Wide y0;
    Wide y1w;
    Wide y2w;
    Wide y3w;
    final_exp_detail::cyclotomic_square_to(y0, ret, fp12ring, ring);
    final_exp_detail::raise_to_z_to(y1w, y2w, y0, fp12ring, ring);
    final_exp_detail::raise_to_z_div_by_2_to(
        y2w, y3w, y1w, fp12ring, ring);

    y1 = fp12ring.fp12_flat_conjugate(ret, ring);
    y1 = fp12ring.fp12_flat_mul(y1w, y1, ring);
    final_exp_detail::conjugate_in_place(y1, fp12ring, ring);
    y1 = fp12ring.fp12_flat_mul(y1, y2w, ring);

    final_exp_detail::raise_to_z_to(y1w, y2w, y1, fp12ring, ring);
    final_exp_detail::raise_to_z_to(y2w, y3w, y1w, fp12ring, ring);

    Standard y3 = fp12ring.fp12_flat_conjugate(y1, ring);
    y3 = fp12ring.fp12_flat_mul(y2w, y3, ring);

    y1 = fp12ring.fp12_flat_frobenius_map3(y1, ring);
    y1w = fp12ring.fp12_flat_frobenius_map2(y1w, ring);
    y1 = fp12ring.fp12_flat_mul(y1, y1w, ring);

    final_exp_detail::raise_to_z_to(y1w, y2w, y3, fp12ring, ring);
    Standard tail = fp12ring.fp12_flat_mul(y1w, y0, ring);
    tail = fp12ring.fp12_flat_mul(tail, ret, ring);
    y1 = fp12ring.fp12_flat_mul(y1, tail, ring);
    tail = fp12ring.fp12_flat_frobenius_map1(y3, ring);
    ret = fp12ring.fp12_flat_mul(y1, tail, ring);
    return ret;
}

// Requires a valid nonzero Fp12 input. The selected inverter is variable-time;
// callers must not treat this boundary as secret-independent. No exception,
// allocation, or semantic validation is introduced by this wrapper.
template<class Fp12RingType, class Ring, class FP12Inverter>
inline auto final_exp(
    const typename Fp12RingType::StandardElement& f,
    const FP12Inverter& inverter,
    const Fp12RingType& fp12ring,
    const Ring& ring)
{
    const auto f_inverse = inverter.invert(f, fp12ring, ring);
    return final_exp_after_inverse(f, f_inverse, fp12ring, ring);
}

} // namespace blst::avx

#undef BLST_AVX_FINAL_EXP_NOINLINE
