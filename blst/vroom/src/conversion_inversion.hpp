#pragma once
#include "../cpu/vector/vector_impl.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include <array>
#include <stdint.h>
#include <cstring>
#include <string>
#include "constexpr_utils.hpp"
#include "../blst/vect.h"
#include "../blst/fields.h"

// Forward declaration - function is in blst/recip.c
extern "C" {
    void ct_inverse_mod_384_wrapper(vec384 out, const vec384 inp);
    void blst_fp12_inverse(vec384fp12 ret, const vec384fp12 a);
}
#include "../blst/consts.h"
#include "fp12.hpp"

template<int limbs>
AVXVector<limbs> to_vec(uint8_t *le_buf) {
    // Convert from vec384 format (6×64-bit limbs) to RNS format (8×52-bit limbs)
    // Input: 384 bits in little-endian 64-bit limbs
    // Output: 8×52-bit limbs (416 bits, but top bits will be zero)
    
    static_assert(limbs == 8, "to_vec currently only supports 8 limbs");
    
    uint64_t *limbs_64 = (uint64_t *)le_buf;
    std::array<uint64_t, 8> digits;
    
    // Extract 8×52-bit limbs from 6×64-bit limbs
    // Limb 0: bits [0:52) from limbs_64[0]
    digits[0] = limbs_64[0] & ((1ULL << 52) - 1);
    
    // Limb 1: bits [52:104) = bits [52:64) from limbs_64[0] + bits [0:40) from limbs_64[1]
    digits[1] = ((limbs_64[0] >> 52) | (limbs_64[1] << 12)) & ((1ULL << 52) - 1);
    
    // Limb 2: bits [104:156) = bits [40:64) from limbs_64[1] + bits [0:28) from limbs_64[2]
    digits[2] = ((limbs_64[1] >> 40) | (limbs_64[2] << 24)) & ((1ULL << 52) - 1);
    
    // Limb 3: bits [156:208) = bits [28:64) from limbs_64[2] + bits [0:16) from limbs_64[3]
    digits[3] = ((limbs_64[2] >> 28) | (limbs_64[3] << 36)) & ((1ULL << 52) - 1);
    
    // Limb 4: bits [208:260) = bits [16:64) from limbs_64[3] + bits [0:4) from limbs_64[4]
    digits[4] = ((limbs_64[3] >> 16) | (limbs_64[4] << 48)) & ((1ULL << 52) - 1);
    
    // Limb 5: bits [260:312) = bits [4:56) from limbs_64[4]
    digits[5] = (limbs_64[4] >> 4) & ((1ULL << 52) - 1);
    
    // Limb 6: bits [312:364) = bits [56:64) from limbs_64[4] + bits [0:44) from limbs_64[5]
    digits[6] = ((limbs_64[4] >> 56) | (limbs_64[5] << 8)) & ((1ULL << 52) - 1);
    
    // Limb 7: bits [364:416) = bits [44:64) from limbs_64[5] (only 20 bits, rest is zero)
    digits[7] = (limbs_64[5] >> 44) & ((1ULL << 52) - 1);
    
    return AVXVector<limbs>(digits);
}

template<bool montgomery=false>
void from_vec(AVXVector<8> &vec_lo, AVXVector<8> &vec_hi, vec384 ret) {
    std::array<uint64_t, 8> digits_lo;
    std::array<uint64_t, 8> digits_hi;
    vec_hi.store(digits_hi.data());
    vec_lo.store(digits_lo.data());
    std::array<uint64_t, 8 + 1> digits;
    digits[0] = digits_lo[0];
    for (int i = 1; i < 8; i++) {
        digits[i] = digits_hi[i - 1] + digits_lo[i];
    }
    digits[8] = digits_hi[7];
    vec768 n_data_even_64;
    vec768 n_data_odd_64;
    uint8_t *n_data_even = (uint8_t *) n_data_even_64;
    uint8_t *n_data_odd = (uint8_t *) n_data_odd_64;
    memset(n_data_even, 0, sizeof(vec768));
    memset(n_data_odd, 0, sizeof(vec768));
    for (int i = 0; i < 9; i += 2) {
        *(uint64_t *)(n_data_even + (i / 2) * 13) = digits[i];
    }
    for (int i = 1; i < 8; i += 2) {
        *(uint64_t *)(n_data_odd + (i / 2) * 13 + 6) = digits[i] << 4;
    }
    
    // Add n_data_even and n_data_odd together with carry propagation
    vec768 sum;
    memset(sum, 0, sizeof(vec768));
    uint64_t carry = 0;
    for (int i = 0; i < 9; i++) {
        uint64_t a = n_data_even_64[i];
        uint64_t b = n_data_odd_64[i];
        uint64_t s = a + b + carry;
        sum[i] = s;
        carry = (s < a) || (carry && s == a) ? 1 : 0;
    }
    // Always: redc_fp -> multiply by RR -> redc_fp (gives normal form)
    vec384 temp;
    redc_fp(temp, sum);  // temp = sum * R^-1 mod p
    vec768 prod;
    mul_384(prod, temp, BLS12_381_RR);  // prod = sum * R
    redc_fp(ret, prod);  // ret = sum (normal form)
    
    if (montgomery) {
        // For Montgomery form: multiply by RR again
        vec768 prod2;
        mul_384(prod2, ret, BLS12_381_RR);  // prod2 = sum * R^2
        redc_fp(ret, prod2);  // ret = sum * R (Montgomery form)
    }
}


template<class Ring, bool montgomery=false>
void convert_to_int(const Ring & ring, typename Ring::StandardElement & a, vec384 ret) {
    AVXVector<8> out_hi = AVXVector<8>(0);
    AVXVector<8> out_lo = AVXVector<8>(0);
    ring.sys.convert_from.conversion_matrix.rns_reduce_raw(a.m2.data, out_hi, out_lo);
    from_vec<montgomery>(out_lo, out_hi, ret);
}

template<class Ring, bool montgomery=false>
auto convert_to_rns(const Ring & ring, const vec384 ret) {
    vec384 temp_ret;
    if constexpr (montgomery) {
        from_fp(temp_ret, ret);
    } else {
        vec_copy(temp_ret, ret, sizeof(vec384));
    }
    AVXVector<8> l = to_vec<8>((uint8_t *) temp_ret);
    return ring.template convert_to_rns_batch<1>(std::array<AVXVector<8>, 1>{l})[0];
}

template<class Ring>
auto invert_via_blst(const Ring & ring, typename Ring::StandardElement & a) {
    vec384 ret;
    convert_to_int<Ring, false>(ring, a, ret);
    ct_inverse_mod_384_wrapper(ret, ret);
    return convert_to_rns<Ring, false>(ring, ret);
}

template<class Ring>
class BLSTInverter {
    public:
    using StandardElement = typename Ring::StandardElement;

    StandardElement invert(const StandardElement &a, const Ring &ring) const {
        return invert_via_blst(ring, const_cast<StandardElement&>(a));
    }
};

// Convert FP12 (12 ring elements) from ring form to BLST Montgomery form,
// call blst_fp12_inverse, convert back to ring form
template<class Ring>
std::array<typename Ring::StandardElement, 12> invert_fp12_via_blst(
    const Ring& ring,
    const std::array<typename Ring::StandardElement, 12>& fp12_ring
) {
    // Step 1: Convert ring form to BLST Montgomery form using convert_to_int
    vec384fp12 fp12_mont;
    // vec384fp12 is vec384fp6[2], vec384fp6 is vec384x[3], vec384x is vec384[2]
    // For index i (0-11): fp12_mont[i/6][(i%6)/2][i%2]
    for (int i = 0; i < 12; i++) {
        convert_to_int<Ring, true>(ring, const_cast<typename Ring::StandardElement&>(fp12_ring[i]), 
                       fp12_mont[i/6][(i%6)/2][i%2]);
    }
    
    // Step 2: Call blst_fp12_inverse
    vec384fp12 fp12_inv_mont;
    blst_fp12_inverse(fp12_inv_mont, fp12_mont);
    
    // Step 3: Convert BLST Montgomery form back to ring form using convert_to_rns
    std::array<typename Ring::StandardElement, 12> fp12_inv_ring;
    for (int i = 0; i < 12; i++) {
        fp12_inv_ring[i] = convert_to_rns<Ring, true>(ring, fp12_inv_mont[i/6][(i%6)/2][i%2]);
    }
    
    return fp12_inv_ring;
}

// --- BigInt-based FP12 inversion via BLST ---
// This path converts ring -> BigInt (normal form) -> vec384 Montgomery (BLST),
// calls blst_fp12_inverse, then converts back to BigInt (normal form) -> ring.
static inline void bigint_to_vec384_normal(vec384 out, const BigInt& a) {
    // Assumes `a` is already in normal form and non-negative (e.g. output of ring.to_bigint()).
    std::string hex_str = a.to_string(16);
    while (hex_str.length() < 96) {
        hex_str = "0" + hex_str;
    }

    memset(out, 0, sizeof(vec384));
    for (int i = 0; i < 6; i++) {
        std::string limb_hex = hex_str.substr((5 - i) * 16, 16);
        out[i] = std::stoull(limb_hex, nullptr, 16);
    }
}

static inline BigInt vec384_normal_to_bigint(const vec384 a) {
    BigInt result(0);
    BigInt two_to_64 = BigInt(1) << 64;
    for (int i = 5; i >= 0; i--) {
        result = result * two_to_64 + BigInt(static_cast<unsigned long>(a[i]));
    }
    return result;
}

static inline void bigint_to_vec384_montgomery(vec384 out, const BigInt& a) {
    vec384 normal;
    bigint_to_vec384_normal(normal, a);

    // out = normal * R^2 mod p
    vec768 prod;
    mul_384(prod, normal, BLS12_381_RR);
    redc_fp(out, prod);
}

static inline BigInt vec384_montgomery_to_bigint(const vec384 a) {
    vec384 normal;
    from_fp(normal, a);
    return vec384_normal_to_bigint(normal);
}

template<class Ring>
std::array<typename Ring::StandardElement, 12> invert_fp12_via_blst_bigint(
    const Ring& ring,
    const std::array<typename Ring::StandardElement, 12>& fp12_ring
) {
    vec384fp12 fp12_mont;
    for (int i = 0; i < 12; i++) {
        BigInt v = ring.to_bigint(fp12_ring[i]);
        bigint_to_vec384_montgomery(fp12_mont[i/6][(i%6)/2][i%2], v);
    }

    vec384fp12 fp12_inv_mont;
    blst_fp12_inverse(fp12_inv_mont, fp12_mont);

    std::array<BigInt, 12> inv_bigints;
    for (int i = 0; i < 12; i++) {
        inv_bigints[i] = vec384_montgomery_to_bigint(fp12_inv_mont[i/6][(i%6)/2][i%2]);
    }

    return ring.template from_bigint_batch<12>(inv_bigints);
}

// FP12 inversion provider using invert_fp12_via_blst_bigint
template<class Ring>
class FP12BLSTInverter {
public:
    using StandardElement = std::array<typename Ring::StandardElement, 12>;
    using FP12R = FP12Ring<Ring>;
    StandardElement invert(const StandardElement& f, const FP12R& /* fp12ring */, const Ring& ring) const {
        return invert_fp12_via_blst(ring, f);
    }
};

template<class Ring>
class FP12RNS_BLSTInverter {
public:
    using FP12R = FP12Ring<Ring>;
    using StandardElement = typename FP12R::StandardElement;

    using Inverter = BLSTInverter<Ring>;
    Inverter inverter;

    StandardElement invert(const StandardElement& f, const FP12R& fp12ring, const Ring& ring) const {
        return fp12ring.template fp12_flat_inverse<Inverter>(f, inverter, ring);
    }
};