#pragma once
#include "fp12.hpp"
#include <cstddef>

template<class Element1, class Element2, class Ring>
inline auto mul_n_sqr(const Element1 &a, const Element2 &b, size_t n, const FP12Ring<Ring> &fp12ring, const Ring &ring) {
    auto result = fp12ring.recast_common(fp12ring.fp12_flat_mul(a, b, ring));
    for (size_t i = 0; i < n; i++) {
        result = fp12ring.fp12_flat_cyclotomic_sqr(result, ring);
    }
    return result;
}

template<class Element, class Ring>
inline auto raise_to_z_div_by_2(const Element &a, const FP12Ring<Ring> &fp12ring, const Ring &ring) {
    auto ret = fp12ring.fp12_flat_cyclotomic_sqr(a, ring);
    ret = mul_n_sqr(ret, a, 2, fp12ring, ring);
    ret = mul_n_sqr(ret, a, 3, fp12ring, ring);
    ret = mul_n_sqr(ret, a, 9, fp12ring, ring);
    ret = mul_n_sqr(ret, a, 32, fp12ring, ring);
    ret = mul_n_sqr(ret, a, 16-1, fp12ring, ring);
    return fp12ring.fp12_flat_conjugate(ret, ring);
}

template<class Element, class Ring>
inline auto raise_to_z(const Element &b, const FP12Ring<Ring> &fp12ring, const Ring &ring) {
    auto a = raise_to_z_div_by_2(b, fp12ring, ring);
    return fp12ring.fp12_flat_cyclotomic_sqr(a, ring);
}

template<class Ring>
inline auto final_exp_after_inverse(const typename FP12Ring<Ring>::StandardElement &f, const typename FP12Ring<Ring>::StandardElement &f_inverse, const FP12Ring<Ring> &fp12ring, const Ring &ring) {
    // blst pairing.c: 375
    
    auto f_conj = fp12ring.fp12_flat_conjugate(f, ring);
    
    // This is the only use of f_conj and f_inverse and cancels scalar constants...
    // But you could mess with it into f_inverse or something.
    auto ret1 = fp12ring.fp12_flat_mul(f_conj, f_inverse, ring);
    // below is a function of only ret1.
    auto frob2_ret1 = fp12ring.fp12_flat_frobenius_map2(ret1, ring);
    auto ret2 = fp12ring.fp12_flat_mul(frob2_ret1, ret1, ring);
    
    auto ret2_sqr = fp12ring.fp12_flat_cyclotomic_sqr(ret2, ring);
    
    auto y0_to_z = raise_to_z(ret2_sqr, fp12ring, ring);
    
    auto y0_to_z_div2 = raise_to_z_div_by_2(y0_to_z, fp12ring, ring);
    
    auto ret2_conj = fp12ring.fp12_flat_conjugate(ret2, ring);
    auto y1_mul_y3 = fp12ring.fp12_flat_mul(ret2_conj, y0_to_z, ring);
    auto y1_mul_y3_conj = fp12ring.fp12_flat_conjugate(y1_mul_y3, ring);
    auto y1_final = fp12ring.fp12_flat_mul(y1_mul_y3_conj, y0_to_z_div2, ring);
    
    auto y1_to_z = raise_to_z(y1_final, fp12ring, ring);
    
    auto y2_to_z = raise_to_z(y1_to_z, fp12ring, ring);

    auto y1_final_conj = fp12ring.fp12_flat_conjugate(y1_final, ring);
    auto y3_mul_y1 = fp12ring.fp12_flat_mul(y2_to_z, y1_final_conj, ring);
    // auto y1_final_conj2 = fp12ring.fp12_flat_conjugate(y1_final_conj, ring);

    auto y1_frob3 = fp12ring.fp12_flat_frobenius_map3(y1_final, ring);
    auto y1_to_z_frob2 = fp12ring.fp12_flat_frobenius_map2(y1_to_z, ring);
    auto y1_mul_y2 = fp12ring.fp12_flat_mul(y1_frob3, y1_to_z_frob2, ring);

    auto y3_to_z = raise_to_z(y3_mul_y1, fp12ring, ring);
    auto y2_mul_y0 = fp12ring.fp12_flat_mul(y3_to_z, ret2_sqr, ring);
    auto y2_mul_ret = fp12ring.fp12_flat_mul(y2_mul_y0, ret2, ring);
    auto y1_mul_y2_final = fp12ring.fp12_flat_mul(y1_mul_y2, y2_mul_ret, ring);
    auto y3_frob1 = fp12ring.fp12_flat_frobenius_map1(y3_mul_y1, ring);
    return fp12ring.fp12_flat_mul(y1_mul_y2_final, y3_frob1, ring);
}

// FP12 inversion provider interface
// Any type that provides an invert method with signature:
//   std::array<typename Ring::StandardElement, 12> invert(
//       const std::array<typename Ring::StandardElement, 12>& f,
//       const Ring& ring
//   ) const;
template<class Ring, class FP12Inverter>
inline auto final_exp(
    const typename FP12Ring<Ring>::StandardElement &f,
    const FP12Inverter& inverter,
    const FP12Ring<Ring> &fp12ring,
    const Ring &ring
) {
    // Use the inverter provider to compute f inverse
    auto f_inverse = inverter.invert(f, fp12ring, ring);
    
    // Call final_exp_after_inverse with f and f_inverse
    return final_exp_after_inverse(f, f_inverse, fp12ring, ring);
}
