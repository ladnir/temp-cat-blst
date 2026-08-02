#pragma once
#include "../cpu/vector/vector_impl.hpp"
#include "../common/constexpr_utils.h"
#include "ring_element.h"
#include "bounded_element.h"
#include "bounds.h"
#include <cstdint>
#include <cstddef>
#include <array>
#include "../common/preprocess.h"

namespace blst::avx {

template<class WideM2, class ReducedM1>
class ReadyToReduce {
    public:
    WideM2 acc;
    ReducedM1 input;
    
    // Constructor for aggregate initialization compatibility
    ReadyToReduce() = default;
    ReadyToReduce(const WideM2 &a, const ReducedM1 &i) : acc(a), input(i) {}
};

// TODO: add chunk size parameter and chunk accordingly?
template<std::size_t batch, class RNSMatrix, class WideM2, class ReducedM1>
INLINE auto batch_change_base(const std::array<ReadyToReduce<WideM2, ReducedM1>, batch> &inputs, const RNSMatrix &rns_matrix) {
    constexpr int limbs = ReducedM1::LIMBS;
    constexpr int element_bits = ReducedM1::ELEMENT_BITS;
    ReducedM1::check_mult_bounds();
    // Note: RNSMatrix doesn't expose template parameters as static members, so we can't static_assert here
    // The template parameters are checked at the call site
    std::array<AVXVector<limbs>, batch> out_hi;
    std::array<AVXVector<limbs>, batch> out_lo;
    std::array<AVXVector<limbs>, batch> residues1;
    for (size_t i = 0; i < batch; i++) {
        residues1[i] = inputs[i].input.data;
        out_hi[i] = inputs[i].acc.high.data;
        out_lo[i] = inputs[i].acc.low.data;
    }
    rns_matrix.template rns_reduce_raw_batch_shifted<batch>(&residues1, out_hi, out_lo);
    // Adds m1_bound limbs time in accumulate loop
    // correction is bounded by m1_bound * limbs, and add_correction adds another m1_bound * limbs
    auto change_base_bounds = Bounds<0, 2*limbs*ReducedM1::bounds.upper>{};
    auto total_bounds = inputs[0].acc.bounds + change_base_bounds;
    std::array<WideElement<decltype(total_bounds), limbs, element_bits>, batch> result;
    for (size_t i = 0; i < batch; i++) {
        result[i].low.data = out_lo[i];
        result[i].high.data = out_hi[i];
    }
    return result;

}

}  // namespace blst::avx
