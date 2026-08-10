#pragma once
#include "vector_impl.hpp"
#include "changebase.hpp"
#include "../reduction/montgomery.hpp"
#include "../../common/constexpr_utils.h"
#include <array>


namespace blst::avx {


// The convenience constructor in RNSMatrix now uses IN_DIGITS and OUT_DIGITS from class template

template<int bits, int limbs_out, int in_word_bits = 52>
class ConvertToRNS {
    public:
    constexpr static int limbs_in = 1;  // Single modulus (target)
    constexpr static int indigits = ceil_div(bits, in_word_bits);  // Number of digits needed to represent target
    const RNSMatrix<limbs_in, limbs_out, indigits, 1, in_word_bits> conversion_matrix;  // OUT_DIGITS=1 (not bit width!)



    inline ConvertToRNS(
        const std::array<std::array<uint64_t, limbs_out>, indigits>& rns_values_tf,
        const std::array<uint64_t, limbs_out>& correction_vec,
        const std::array<uint64_t, indigits>& sqe_tf,
        const std::array<uint64_t, limbs_out>& moduli_out)
        : conversion_matrix(rns_values_tf, correction_vec, sqe_tf, moduli_out) {}


    // Test accessor
    const RNSMatrix<limbs_in, limbs_out, indigits, 1, in_word_bits>& get_conversion_matrix() const {
        return conversion_matrix;
    }
};

template<int bits,int limbs_in>
class ConvertFromRNS {
    
    public:
    // Python uses word_len(2*target, outbase) = (bit_length(2*target) + outbase - 1) // outbase
    // For 381-bit target, if applied to the output of a multiply, could be up to 3*target, so we need to add extra bits to handle the redundant output in binary.
    // The conversion process introduces more redundance, so a total of 4 bits are needed.
    static constexpr int outdigits = ceil_div(bits + 4, 52);
    static constexpr int outdigits_rounded = AVXVector<outdigits>::VEC_LIMBS * AVXVector<outdigits>::LIMBS_PER_VEC; // Rounds up to the nearest LIMBS_PER_VEC. Needed to prevent buffer overflow in store
    const RNSMatrix<limbs_in, 1, 1, outdigits> conversion_matrix;  // OUT_DIGITS=outdigits

    inline ConvertFromRNS(
        const std::array<std::array<uint64_t, outdigits>, limbs_in>& rns_values_tf,
        const std::array<uint64_t, outdigits>& correction_vec,
        const std::array<uint64_t, limbs_in>& sqe_tf,
        const std::array<uint64_t, 1>& moduli_out)
        : conversion_matrix(rns_values_tf, correction_vec, sqe_tf, moduli_out) {}

    
    // Test accessor
    const RNSMatrix<limbs_in, 1, 1, outdigits>& get_conversion_matrix() const {
        return conversion_matrix;
    }
};

}  // namespace blst::avx
