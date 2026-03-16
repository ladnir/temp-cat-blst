#pragma once
#include "vector_impl.hpp"
#include "changebase.hpp"
#include "../reduction/montgomery.hpp"
#include "../precompute/gmp_wrapper.hpp"
#include <array>

// Helper function to reinterpret a BigInt as words of specified bit length
template<int limbs, int word_length>
inline void word_reinterpret(const BigInt &integer, std::array<uint64_t, limbs> &digits) {
    BigInt mask = (BigInt(1) << word_length) - 1;
    for (int i = 0; i < limbs; i++) {
        BigInt shifted = integer >> (word_length * i);
        digits[i] = static_cast<uint64_t>((shifted & mask).to_ulong());
    }
}

// The convenience constructor in RNSMatrix now uses IN_DIGITS and OUT_DIGITS from class template

template<int bits, int limbs_out>
class ConvertToRNS {
    public:
    constexpr static int limbs_in = 1;  // Single modulus (target)
    constexpr static int indigits = ceil_div(bits, 52);  // Number of digits needed to represent target
    const RNSMatrix<limbs_in, limbs_out, indigits, 1> conversion_matrix;  // OUT_DIGITS=1 (not bit width!)


    inline ConvertToRNS(const BigInt &premult, const BigInt &postmult, const BigInt &target, const std::array<uint64_t, limbs_out> &moduli_out) 
        : conversion_matrix(std::array<BigInt,1>{target}, moduli_out, BigInt(0), premult, postmult) {
        // Now uses IN_DIGITS=indigits (8) and OUT_DIGITS=1 from class template
    }

    template<template <int> class Reduction>
    inline AVXVector<limbs_out> convert_to_rns(const BigInt &integer, const Reduction<limbs_out> &reducer) const {
        std::array<uint64_t, indigits> digits;
        word_reinterpret<indigits, 52>(integer, digits);
        AVXVector<indigits> input;
        input.load(digits.data());
        return conversion_matrix.rns_reduce(input, reducer);
    }

    template<int batch>
    inline std::pair<std::array<AVXVector<limbs_out>, batch>, std::array<AVXVector<limbs_out>, batch>> batch_convert_to_rns(const std::array<BigInt, batch> &integers) const {
        std::array<AVXVector<indigits>, batch> inputs;
        for (int i = 0; i < batch; i++) {
            // Definitely could be optimized...
            std::array<uint64_t, indigits> digits;
            word_reinterpret<indigits, 52>(integers[i], digits);
            inputs[i].load(digits.data());
        }
        std::array<AVXVector<limbs_out>, batch> out_hi{};
        std::array<AVXVector<limbs_out>, batch> out_lo{};
        conversion_matrix.template rns_reduce_raw_batch<batch>(&inputs, out_hi, out_lo);
        return std::make_pair(out_hi, out_lo);
    }

    // Test accessor
    const RNSMatrix<limbs_in, limbs_out, indigits, 1>& get_conversion_matrix() const {
        return conversion_matrix;
    }
};

// Helper function to reconstruct a BigInt from words of specified bit length
// Could probably be done more efficiently with bit operations.
template<int limbs, int in_limbs, int word_length>
inline BigInt int_reconstruct(const std::array<uint64_t, in_limbs> &digits) {
    static_assert(limbs <= in_limbs, "limbs must be less than or equal to in_limbs");
    BigInt result(0);
    for (int i = 0; i < limbs; i++) {
        result += BigInt(digits[i]) << (word_length * i);
    }
    return result;
}

template<int bits,int limbs_in>
class ConvertFromRNS {
    
    public:
    // Python uses word_len(2*target, outbase) = (bit_length(2*target) + outbase - 1) // outbase
    // For 381-bit target, if applied to the output of a multiply, could be up to 3*target, so we need to add extra bits to handle the redundant output in binary.
    // The conversion process introduces more redundance, so a total of 4 bits are needed.
    static constexpr int outdigits = ceil_div(bits + 4, 52);
    static constexpr int outdigits_rounded = AVXVector<outdigits>::VEC_LIMBS * AVXVector<outdigits>::LIMBS_PER_VEC; // Rounds up to the nearest LIMBS_PER_VEC. Needed to prevent buffer overflow in store
    const RNSMatrix<limbs_in, 1, 1, outdigits> conversion_matrix;  // OUT_DIGITS=outdigits
    const BigInt target;
    inline ConvertFromRNS(const BigInt &premult, const BigInt &postmult, const BigInt &target, const std::array<uint64_t, limbs_in> &moduli_in)
        : conversion_matrix(moduli_in, std::array<BigInt,1>{target}, BigInt(0), premult, postmult), target(target) {}

    inline BigInt convert_from_rns(const AVXVector<limbs_in> &residues) const {
        AVXVector<outdigits> out_hi = AVXVector<outdigits>(0);
        AVXVector<outdigits> out_lo = AVXVector<outdigits>(0);
        conversion_matrix.rns_reduce_raw(residues, out_hi, out_lo);
        // Combine hi and lo into a single BigInt
        std::array<uint64_t, outdigits_rounded> digits_lo;
        std::array<uint64_t, outdigits_rounded> digits_hi;
        // std::cout << "outdigits=" << outdigits << std::endl;
        // std::cout << "outdigits_rounded=" << outdigits_rounded << std::endl;
        out_lo.store(digits_lo.data());
        out_hi.store(digits_hi.data());
        // Eventually we could rotate vectors and use avx512 addition (creating 53 bit numbers), and then do the extraction.
        // We could also get the quotient directly from the matrix instead of % target, to within a factor of 2 (therefore requiring only a conditional correction).
        BigInt i = int_reconstruct<outdigits, outdigits_rounded, 52>(digits_lo);
        i += int_reconstruct<outdigits, outdigits_rounded, 52>(digits_hi) << 52;
        return i % target;
    }
    
    // Test accessor
    const RNSMatrix<limbs_in, 1, 1, outdigits>& get_conversion_matrix() const {
        return conversion_matrix;
    }
};
