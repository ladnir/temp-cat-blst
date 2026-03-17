#pragma once
#include "miller.hpp"
#include "final_exponentiation.hpp"
#include "conversion_inversion.hpp"

// Pairing function that combines miller loop and final exponentiation
// Similar to blst pairing.c: blst_miller_loop followed by blst_final_exp
template<class Ring, class FP12Inverter>
inline auto pairing(
    const AffinePoint<typename Ring::StandardElement> &P,  // G1 point
    const AffinePoint<typename FP2Ring<Ring>::StandardElement> &Q,  // G2 point
    const FP12Inverter& inverter,
    const FP2Ring<Ring>& fp2ring,
    const FP12Ring<Ring>& fp12ring,
    const Ring &ring
) {
    // Create Miller instance with G1 and G2 points
    Miller<Ring> miller(P, Q, ring);
    
    // Step 1: Compute miller loop (returns f)
    auto f = miller.miller_loop(ring, fp2ring, fp12ring);
    
    // Step 2: Compute final exponentiation using the inverter provider
    return final_exp(f, inverter, fp12ring, ring);
}

// Convenience wrapper using FP12BLSTInverter
template<class Ring>
inline auto pairing_with_blst_inverter(
    const AffinePoint<typename Ring::StandardElement> &P,  // G1 point
    const AffinePoint<typename FP2Ring<Ring>::StandardElement> &Q,  // G2 point
    const FP2Ring<Ring>& fp2ring,
    const FP12Ring<Ring>& fp12ring,
    const Ring &ring
) {
    FP12RNS_BLSTInverter<Ring> inverter;
    return pairing(P, Q, inverter, fp2ring, fp12ring, ring);
}

