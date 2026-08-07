#pragma once
#include "miller.h"
#include "blst/avx/fields/final_exponentiation.h"
#include "blst/avx/fields/conversion_inversion.h"

// Pairing function that combines miller loop and final exponentiation
// Similar to blst pairing.c: blst_miller_loop followed by blst_final_exp
namespace blst::avx {

template<class Ring, class FP12RingType, class FP12Inverter>
inline auto pairing(
    const AffinePoint<typename Ring::StandardElement> &P,  // G1 point
    const AffinePoint<typename FP2Ring<Ring>::StandardElement> &Q,  // G2 point
    const FP12Inverter& inverter,
    const FP2Ring<Ring>& fp2ring,
    const FP12RingType& fp12ring,
    const Ring &ring
) {
    // Create Miller instance with G1 and G2 points
    Miller<Ring, FP12RingType> miller(P, Q, ring);
    
    // Step 1: Compute miller loop (returns f)
    auto f = miller.miller_loop(ring, fp2ring, fp12ring);
    
    // Step 2: Compute final exponentiation using the inverter provider
    return final_exp(f, inverter, fp12ring, ring);
}

// Convenience wrapper using FP12BLSTInverter
template<class Ring, class FP12RingType>
inline auto pairing_with_blst_inverter(
    const AffinePoint<typename Ring::StandardElement> &P,  // G1 point
    const AffinePoint<typename FP2Ring<Ring>::StandardElement> &Q,  // G2 point
    const FP2Ring<Ring>& fp2ring,
    const FP12RingType& fp12ring,
    const Ring &ring
) {
    FP12RNS_BLSTInverter<Ring, FP12RingType> inverter;
    return pairing(P, Q, inverter, fp2ring, fp12ring, ring);
}

}  // namespace blst::avx
