#pragma once
#include "changebase.hpp"
#include "vector_impl.hpp"
#include "../precompute/precompute.hpp"
#include "../reduction/montgomery.hpp"
#include "conversion.hpp"
#include <tuple>

template <template <int> class Reduction, int limbs>
inline AVXVector<limbs> ele_mult(const AVXVector<limbs> &a, const AVXVector<limbs> & b, const Reduction<limbs> &r) {
    auto ab_lo = AVXVector<limbs>(0).mullo(a, b);
    auto ab_hi = AVXVector<limbs>(0).mulhi(a, b);
    return r.reduce_small(ab_hi, ab_lo);
}

template<template <int> class Reduction, int limbs1, int limbs2>
class IntRNS2 {
    
    public:

    const RNSMatrix<limbs1, limbs2> r1;
    const RNSMatrix<limbs2, limbs1> r2;
    const Reduction<limbs1> m1;
    const Reduction<limbs2> m2;

    // Constructor building precompute with separate premult/postmult for each matrix
    inline IntRNS2(
        const std::array<uint64_t, limbs1> &moduli_in_m1,
        const std::array<uint64_t, limbs2> &moduli_in_m2,
        const BigInt &target_modulus,
        const BigInt &premult1,
        const BigInt &postmult1,
        const BigInt &premult2,
        const BigInt &postmult2
    )
    : r1(moduli_in_m1, moduli_in_m2, BigInt(0), premult1, postmult1)  // target_modulus=0 means use modulus_out (matching Python None)
    , r2(moduli_in_m2, moduli_in_m1, BigInt(0), premult2, postmult2)  // target_modulus=0 means use modulus_out (matching Python None)
    , m1(moduli_in_m1)
    , m2(moduli_in_m2)
    {
        (void)target_modulus;
    }

    inline void mul_reduce(const AVXVector<limbs1> &a1, const AVXVector<limbs2> &a2, const AVXVector<limbs1> &b1, const AVXVector<limbs2> &b2, AVXVector<limbs1> &out1, AVXVector<limbs2> &out2) const {
        auto ab_m1 = ele_mult(a1, b1, m1);
        out2 = r1.rns_reduce_with_acc(ab_m1, a2, b2, m2);
        out1 = r2.rns_reduce(out2, m1);
    }

    inline void mul_wide(const AVXVector<limbs1> &a1, const AVXVector<limbs2> &a2, const AVXVector<limbs1> &b1, const AVXVector<limbs2> &b2, AVXVector<limbs1> &out1_lo, AVXVector<limbs1> &out1_hi, AVXVector<limbs2> &out2_lo, AVXVector<limbs2> &out2_hi) const {
        out1_lo = AVXVector<limbs1>(0).mullo(a1, b1);
        out1_hi = AVXVector<limbs1>(0).mulhi(a1, b1);
        out2_lo = AVXVector<limbs2>(0).mullo(a2, b2);
        out2_hi = AVXVector<limbs2>(0).mulhi(a2, b2);
    }

    inline void mul_acc_wide(const AVXVector<limbs1> &a1, const AVXVector<limbs2> &a2, const AVXVector<limbs1> &b1, const AVXVector<limbs2> &b2, AVXVector<limbs1> &out1_lo, AVXVector<limbs1> &out1_hi, AVXVector<limbs2> &out2_lo, AVXVector<limbs2> &out2_hi) const {
        out1_lo = out1_lo.mullo(a1, b1);
        out1_hi = out1_hi.mulhi(a1, b1);
        out2_lo = out2_lo.mullo(a2, b2);
        out2_hi = out2_hi.mulhi(a2, b2);
    }

    inline void reduce_wide(const AVXVector<limbs1> &a1_lo, const AVXVector<limbs1> &a1_hi, const AVXVector<limbs2> &a2_lo, const AVXVector<limbs2> &a2_hi, AVXVector<limbs1> &out1, AVXVector<limbs2> &out2) const {
        auto a1 = m1.reduce_small(a1_hi, a1_lo);
        // need to make copy of a2_hi and a2_lo to avoid modifying the input
        auto a2_hi_copy = a2_hi;
        auto a2_lo_copy = a2_lo;
        r1.rns_reduce_raw(a1, a2_hi_copy, a2_lo_copy);
        out2 = m2.reduce_full(a2_hi_copy, a2_lo_copy);
        out1 = r2.rns_reduce(out2, m1);
    }

    template<int batch>
    inline void reduce_wide_batch(const std::array<AVXVector<limbs1>, batch> &a1_lo, const std::array<AVXVector<limbs1>, batch> &a1_hi, const std::array<AVXVector<limbs2>, batch> &a2_lo, const std::array<AVXVector<limbs2>, batch> &a2_hi, std::array<AVXVector<limbs2>, batch> &out2) const {
        std::array<AVXVector<limbs1>, batch> a1;
        std::array<AVXVector<limbs2>, batch> a2_hi_copy;
        std::array<AVXVector<limbs2>, batch> a2_lo_copy;
        a1 = m1.template reduce_small_batch<batch>(a1_hi, a1_lo);
        for (int i = 0; i < batch; i++) {
            a2_hi_copy[i] = a2_hi[i];
            a2_lo_copy[i] = a2_lo[i];
        }
        r1.template rns_reduce_raw_batch<batch>(a1, a2_hi_copy, a2_lo_copy);
        out2 = m2.template reduce_full_batch<batch>(a2_hi_copy, a2_lo_copy);
    }

    template<int batch>
    inline void batch_expand_wide(const std::array<AVXVector<limbs2>, batch> &out2, std::array<AVXVector<limbs1>, batch> &out1) const {
        std::array<AVXVector<limbs1>, batch> z_hi;
        std::array<AVXVector<limbs1>, batch> z_lo;
        for (int i = 0; i < batch; i++) {
            z_hi[i] = AVXVector<limbs1>(0);
            z_lo[i] = AVXVector<limbs1>(0);
        }
        r2.template rns_reduce_raw_batch<batch>(out2, z_hi, z_lo);
        // Do I also need reduce_full_batch to amortize loading of Montgomery factors?
        out1 = m1.template reduce_full_batch<batch>(z_hi, z_lo);
    }
};

// Factory for Montgomery variant (single template parameter: limbs). Assumes limbs1==limbs2==limbs.
// Generates moduli using GenRNS; aborts if target cannot be reached with the given limb count.
template<int limbs>
inline IntRNS2<MontgomeryReduce, limbs, limbs> make_IntRNS2Montgomery(
    const BigInt &target,
    const int modbits = 52,
    const int max_bound = 9
) {
    // Generate first set of moduli covering target * max_bound
    GenRNS<limbs> g1(target * BigInt(max_bound), -1, BigInt(1), modbits);
    // Generate second set covering ~6*target, chaining state and previous modulus product
    GenRNS<limbs> g2(target * BigInt(6), g1.state, g1.modulus, modbits);

    const BigInt m1 = g1.modulus;
    const BigInt m2 = g2.modulus;
    if (!(m2 > BigInt(6) * target)) {
        throw std::runtime_error("Insufficient limbs for gen RNS: m2 <= 6*target");
    }

    std::array<uint64_t, limbs> moduli1 = g1.moduli;
    std::array<uint64_t, limbs> moduli2 = g2.moduli;

    // Montgomery parameters following python_reduction.IntRNS2
    constexpr int mulbits = 52;
    const BigInt r = BigInt(1) << mulbits;

    // r1 parameters
    const BigInt mont_factor = (BigInt(0) - target).mod_inverse(m1); // pow(-target, -1, m1)
    const BigInt r_inv_m1 = r.mod_inverse(m1);
    const BigInt premult1 = (mont_factor * r_inv_m1) % m1;

    const BigInt inv_factor = m1.mod_inverse(m2); // m1^{-1} mod m2
    const BigInt postmult1 = (target * inv_factor % m2) * (inv_factor % m2) % m2 * (r % m2) % m2 * (r % m2) % m2;

    // r2 parameters
    const BigInt r_inv_m2 = r.mod_inverse(m2);
    const BigInt premult2 = (m1 % m2) * r_inv_m2 % m2;
    const BigInt postmult2 = (r % m1) * (r % m1) % m1; // but postmult is used under modulus_out (m1) inside precompute

    return IntRNS2<MontgomeryReduce, limbs, limbs>(moduli1, moduli2, target, premult1, postmult1, premult2, postmult2);
}

template<class T, std::size_t l1, std::size_t l2>
std::array<T, l1 + l2> concat(const std::array<T, l1> &a, const std::array<T, l2> &b) {
    std::array<T, l1 + l2> result;
    for (std::size_t i = 0; i < l1; i++) {
        result[i] = a[i];
    }
    for (std::size_t i = 0; i < l2; i++) {
        result[l1 + i] = b[i];
    }
    return result;
}

template<int num_moduli>
BigInt rns_reconstruct_full(const std::array<int64_t, num_moduli> &residues, const std::array<uint64_t, num_moduli> &moduli, const BigInt &M1M2) {
    BigInt result = BigInt(0);
    for (int i = 0; i < num_moduli; i++) {
        BigInt m_i = BigInt(moduli[i]);
        BigInt rest = M1M2 / m_i;
        BigInt inv = (rest % m_i).mod_inverse(m_i);
        BigInt icrt = (rest * inv) % M1M2;
        result += icrt * BigInt(residues[i]);
    }
    BigInt result_mod = result % M1M2;
    if (result_mod < 0) {
        result_mod += M1M2;
    }
    return result_mod;
}

template<int bits, int limbs, template <int> class Reduction=MontgomeryReduce>
struct IntRNS2System {
    IntRNS2<Reduction, limbs, limbs> intrns2;
    ConvertToRNS<bits, limbs> convert_to;
    ConvertFromRNS<bits, limbs> convert_from;
    std::array<uint64_t, limbs> moduli1;
    std::array<uint64_t, limbs> moduli2;
    BigInt M1;
    BigInt M2;
    BigInt M1M2;
    BigInt target;
    BigInt total_rotation;
    BigInt total_inv_rotation;
    BigInt inv_factor_target;

    inline IntRNS2System(
        IntRNS2<Reduction, limbs, limbs> intrns2,
        ConvertToRNS<bits, limbs> convert_to,
        ConvertFromRNS<bits, limbs> convert_from,
        const std::array<uint64_t, limbs> &moduli1,
        const std::array<uint64_t, limbs> &moduli2,
        const BigInt &target,
        const BigInt &M1,
        const BigInt &M2,
        const BigInt &M1M2,
        const BigInt &total_rotation,
        const BigInt &total_inv_rotation,
        const BigInt &inv_factor_target
    
    ) : intrns2(intrns2), convert_to(convert_to), convert_from(convert_from), moduli1(moduli1), moduli2(moduli2), M1(M1), M2(M2), M1M2(M1M2), target(target), total_rotation(total_rotation), total_inv_rotation(total_inv_rotation), inv_factor_target(inv_factor_target) {
    }

    // Convert a BigInt to Montgomery form in RNS representation
    // Returns (a_m1, a_m2) where a_m1 is in moduli1 and a_m2 is in moduli2
    inline std::pair<AVXVector<limbs>, AVXVector<limbs>> to_mont_avx(const BigInt &a) const {
        AVXVector<limbs> a_m2 = convert_to.convert_to_rns(a, intrns2.m2);
        AVXVector<limbs> a_m1 = intrns2.r2.rns_reduce(a_m2, intrns2.m1);
        return std::make_pair(a_m1, a_m2);
    }

    // For debugging or for the precomputation of constants, returns the RNS representation of the Montgomery form of the input, with redundance 0.
    inline std::pair<std::array<uint64_t, limbs>, std::array<uint64_t, limbs>> to_mont_exact(const BigInt &a, bool wide=false) const {
        // Multiply by the Montgomery factor (or the Montgomery factor^2 for wide type)
        BigInt redundance = a / target;
        BigInt a_mont = (a * M1) % target;
        if (wide) {
            a_mont = (a_mont * M1) % target;
        }
        a_mont += redundance * target;
        a_mont = (a_mont * total_rotation) % M1M2;
        if (wide) {
            a_mont = (a_mont * total_rotation) % M1M2;
        }
        std::array<uint64_t, limbs> a_m1;
        std::array<uint64_t, limbs> a_m2;
        for (int i = 0; i < limbs; i++) {
            a_m1[i] = (a_mont % BigInt(moduli1[i])).to_ulong();
            a_m2[i] = (a_mont % BigInt(moduli2[i])).to_ulong();
        }
        return std::make_pair(a_m1, a_m2);
    }

    // Mainly for debugging, returns the value and the redundance encoded in the RNS representation
    // inline std::pair<BigInt, BigInt>
    // from_mont_exact(AVXVector<limbs> &a_m1, AVXVector<limbs> &a_m2, bool wide=false) const {
    //     // First reconstruct the integer from the RNS representation

    // }

    // Convert from Montgomery form in RNS representation back to BigInt
    // Takes a_m2 (representation in moduli2) and returns the BigInt value
    inline BigInt from_mont_avx(const AVXVector<limbs> &a_m2) const {
        return convert_from.convert_from_rns(a_m2);
    }

    std::tuple<BigInt, BigInt, BigInt> from_mont_exact(const std::array<int64_t, limbs> &a_m1, const std::array<int64_t, limbs> &a_m2, bool wide=false, bool signed_value=true) const {
        std::array<int64_t, 2*limbs> a = concat<int64_t, limbs, limbs>(a_m1, a_m2);
        std::array<uint64_t, 2*limbs> moduli = concat<uint64_t, limbs, limbs>(moduli1, moduli2);
        BigInt a_mont = rns_reconstruct_full<2*limbs>(a, moduli, M1M2);
        a_mont = (a_mont * total_inv_rotation) % M1M2;
        if (wide) {
            a_mont = (a_mont * total_inv_rotation) % M1M2;
        }
        if (signed_value) {
            if (a_mont > M1M2 / 2) {
                a_mont -= M1M2;
            }
        }
        BigInt value = (a_mont * inv_factor_target) % target;
        if (wide) {
            value = (value * inv_factor_target) % target;
        }
        BigInt redundance = a_mont / target;
        
        if (value < 0) {
            value += target;
            redundance -= BigInt(1);
        }
        return std::make_tuple(value, redundance, a_mont);
    }

    // Multiply and reduce in RNS representation
    // Takes a_m1, a_m2, b_m1, b_m2 and computes (a*b)_m1, (a*b)_m2
    inline void mul_reduce(
        const AVXVector<limbs> &a_m1,
        const AVXVector<limbs> &a_m2,
        const AVXVector<limbs> &b_m1,
        const AVXVector<limbs> &b_m2,
        AVXVector<limbs> &out_m1,
        AVXVector<limbs> &out_m2
    ) const {
        intrns2.mul_reduce(a_m1, a_m2, b_m1, b_m2, out_m1, out_m2);
    }

    inline void raw_mul(
        const AVXVector<limbs> &a_m1,
        const AVXVector<limbs> &a_m2,
        const AVXVector<limbs> &b_m1,
        const AVXVector<limbs> &b_m2,
        AVXVector<limbs> &out_m1,
        AVXVector<limbs> &out_m2
    ) const {
        intrns2.raw_mul(a_m1, a_m2, b_m1, b_m2, out_m1, out_m2);
    }
};

// Factory producing IntRNS2 plus compatible converters following python_reduction.IntRNS2
template<int bits, int limbs, template <int> class Reduction=MontgomeryReduce>
inline IntRNS2System<bits, limbs> make_IntRNS2SystemMontgomery(
    const BigInt &target,
    const int modbits = 52,
    const int max_bound = 9
) {
    GenRNS<limbs> g1(target * BigInt(max_bound), -1, BigInt(1), modbits);
    GenRNS<limbs> g2(target * BigInt(6), g1.state, g1.modulus, modbits);
    const BigInt m1 = g1.modulus;
    const BigInt m2 = g2.modulus;
    if (!(m2 > BigInt(6) * target)) {
        throw std::runtime_error("Insufficient limbs for gen RNS: m2 <= 6*target");
    }
    std::array<uint64_t, limbs> moduli1;
    std::array<uint64_t, limbs> moduli2;
    for (int i = 0; i < limbs; i++) {
        moduli1[i] = static_cast<uint64_t>(g1.moduli[i].to_ulong());
        moduli2[i] = static_cast<uint64_t>(g2.moduli[i].to_ulong());
    }

    constexpr int mulbits = 52;
    const BigInt r = BigInt(1) << mulbits;

    const BigInt mont_factor = (BigInt(0) - target).mod_inverse(m1);
    const BigInt r_inv_m1 = r.mod_inverse(m1);
    const BigInt inv_factor = m1.mod_inverse(m2);

    const BigInt premult1 = (mont_factor * r_inv_m1) % m1;
    const BigInt postmult1 = ((target * inv_factor) % m2 * inv_factor % m2 * r % m2 * r) % m2;

    const BigInt r_inv_m2 = r.mod_inverse(m2);
    const BigInt premult2 = (m1 % m2) * r_inv_m2 % m2;
    const BigInt postmult2 = (r % m1) * (r % m1) % m1;

    const BigInt m1_inv_target = m1.mod_inverse(target);

    BigInt M1 = m1;
    BigInt M2 = m2;
    BigInt M1M2 = M1 * M2;
    BigInt icrt1 = m2 * m2.mod_inverse(m1);
    BigInt icrt2 = m1 * inv_factor;
    BigInt rotation_factor = ((1 * icrt1 + inv_factor * icrt2)) % M1M2;
    BigInt total_rotation = (r * rotation_factor) % M1M2;
    BigInt total_inv_rotation = total_rotation.mod_inverse(M1M2);
    BigInt inv_factor_target = m1.mod_inverse(target);
    IntRNS2System<bits, limbs, Reduction> sys(
        IntRNS2<Reduction, limbs, limbs>(moduli1, moduli2, target, premult1, postmult1, premult2, postmult2),
        ConvertToRNS<bits, limbs>(m1, (inv_factor % m2) * (r % m2) % m2 * (r % m2) % m2, target, moduli2),
        ConvertFromRNS<bits, limbs>((m1 % m2) * r_inv_m2 % m2, m1_inv_target, target, moduli2),
        moduli1,
        moduli2,
        target,
        M1,
        M2,
        M1M2,
        total_rotation,
        total_inv_rotation,
        inv_factor_target
    );
    return sys;
}