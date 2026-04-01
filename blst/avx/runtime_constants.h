#pragma once

#include "generated_runtime_constants.hpp"

#include <array>
#include <cstdint>

namespace blst::avx {

struct RuntimeConstantsView {
    const std::array<uint64_t, 8>& moduli1;
    const std::array<uint64_t, 8>& moduli2;
    const std::array<uint64_t, 8>& reducer1_moduli;
    const std::array<uint64_t, 8>& reducer1_mont_factor;
    const std::array<uint64_t, 8>& reducer1_t_squared;
    const std::array<uint64_t, 8>& reducer2_moduli;
    const std::array<uint64_t, 8>& reducer2_mont_factor;
    const std::array<uint64_t, 8>& reducer2_t_squared;
    const std::array<std::array<uint64_t, 8>, 8>& r1_rns_mat;
    const std::array<uint64_t, 8>& r1_correction;
    const std::array<uint64_t, 8>& r1_shifted_quotient_estimations;
    const std::array<std::array<uint64_t, 8>, 8>& r2_rns_mat;
    const std::array<uint64_t, 8>& r2_correction;
    const std::array<uint64_t, 8>& r2_shifted_quotient_estimations;
    const std::array<std::array<uint64_t, 8>, 8>& convert_to_rns_mat;
    const std::array<uint64_t, 8>& convert_to_correction;
    const std::array<uint64_t, 8>& convert_to_shifted_quotient_estimations;
    const std::array<std::array<uint64_t, 8>, 8>& convert_from_rns_mat;
    const std::array<uint64_t, 8>& convert_from_correction;
    const std::array<uint64_t, 8>& convert_from_shifted_quotient_estimations;
    const std::array<uint64_t, 8>& one_m1;
    const std::array<uint64_t, 8>& one_m2;
    const std::array<uint64_t, 8>& modulus_in_rns_m1;
    const std::array<uint64_t, 8>& modulus_in_rns_m2;
    const std::array<uint64_t, 8>& six_modulus_in_rns_m1;
    const std::array<uint64_t, 8>& six_modulus_in_rns_m2;
    const std::array<uint64_t, 8>& max_neg_offset_m1;
    const std::array<uint64_t, 8>& max_neg_offset_m2;
    const std::array<uint64_t, 8>& fp12_offset_m1_high;
    const std::array<uint64_t, 8>& fp12_offset_m1_low;
    const std::array<uint64_t, 8>& fp12_offset_m2_high;
    const std::array<uint64_t, 8>& fp12_offset_m2_low;
};

template<int element_bits>
inline const RuntimeConstantsView& runtime_constants()
{
    static_assert(element_bits == 52 || element_bits == 50, "Only the baked 50-bit and 52-bit AVX configurations are supported.");

    if constexpr (element_bits == 52) {
        static const RuntimeConstantsView constants = {
            generated::fp52::moduli1,
            generated::fp52::moduli2,
            generated::fp52::reducer1_moduli,
            generated::fp52::reducer1_mont_factor,
            generated::fp52::reducer1_t_squared,
            generated::fp52::reducer2_moduli,
            generated::fp52::reducer2_mont_factor,
            generated::fp52::reducer2_t_squared,
            generated::fp52::r1_rns_mat,
            generated::fp52::r1_correction,
            generated::fp52::r1_shifted_quotient_estimations,
            generated::fp52::r2_rns_mat,
            generated::fp52::r2_correction,
            generated::fp52::r2_shifted_quotient_estimations,
            generated::fp52::convert_to_rns_mat,
            generated::fp52::convert_to_correction,
            generated::fp52::convert_to_shifted_quotient_estimations,
            generated::fp52::convert_from_rns_mat,
            generated::fp52::convert_from_correction,
            generated::fp52::convert_from_shifted_quotient_estimations,
            generated::fp52::one_m1,
            generated::fp52::one_m2,
            generated::fp52::modulus_in_rns_m1,
            generated::fp52::modulus_in_rns_m2,
            generated::fp52::six_modulus_in_rns_m1,
            generated::fp52::six_modulus_in_rns_m2,
            generated::fp52::max_neg_offset_m1,
            generated::fp52::max_neg_offset_m2,
            generated::fp52::fp12_offset_m1_high,
            generated::fp52::fp12_offset_m1_low,
            generated::fp52::fp12_offset_m2_high,
            generated::fp52::fp12_offset_m2_low,
        };
        return constants;
    } else {
        static const RuntimeConstantsView constants = {
            generated::fp50::moduli1,
            generated::fp50::moduli2,
            generated::fp50::reducer1_moduli,
            generated::fp50::reducer1_mont_factor,
            generated::fp50::reducer1_t_squared,
            generated::fp50::reducer2_moduli,
            generated::fp50::reducer2_mont_factor,
            generated::fp50::reducer2_t_squared,
            generated::fp50::r1_rns_mat,
            generated::fp50::r1_correction,
            generated::fp50::r1_shifted_quotient_estimations,
            generated::fp50::r2_rns_mat,
            generated::fp50::r2_correction,
            generated::fp50::r2_shifted_quotient_estimations,
            generated::fp50::convert_to_rns_mat,
            generated::fp50::convert_to_correction,
            generated::fp50::convert_to_shifted_quotient_estimations,
            generated::fp50::convert_from_rns_mat,
            generated::fp50::convert_from_correction,
            generated::fp50::convert_from_shifted_quotient_estimations,
            generated::fp50::one_m1,
            generated::fp50::one_m2,
            generated::fp50::modulus_in_rns_m1,
            generated::fp50::modulus_in_rns_m2,
            generated::fp50::six_modulus_in_rns_m1,
            generated::fp50::six_modulus_in_rns_m2,
            generated::fp50::max_neg_offset_m1,
            generated::fp50::max_neg_offset_m2,
            generated::fp50::fp12_offset_m1_high,
            generated::fp50::fp12_offset_m1_low,
            generated::fp50::fp12_offset_m2_high,
            generated::fp50::fp12_offset_m2_low,
        };
        return constants;
    }
}

}  // namespace blst::avx
