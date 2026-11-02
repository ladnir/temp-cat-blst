#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fields/Fr.h"

namespace bte::transparent_bte
{
    /**
     * Public material for the construction in ePrint 2026/2205.
     *
     * This experimental implementation deliberately starts after DKG. The
     * caller supplies a Shamir sharing of `sk`; this object contains only the
     * resulting public values used by the online BTE algorithms.
     */
    struct PublicParameters
    {
        std::size_t threshold = 0;
        blst::G1Point public_key = blst::G1Point::infinity();
        blst::G1Point independent_generator = blst::G1Point::infinity();
        blst::GT public_key_in_gt = blst::GT::identity();
        std::vector<blst::G1Point> verification_keys;
    };

    struct SecretShare
    {
        std::uint32_t index = 0;
        blst::Fr value{};
    };

    struct SetupResult
    {
        PublicParameters public_parameters;
        std::vector<SecretShare> secret_shares;
    };

    /**
     * The four algebraic components of Fig. 1.
     *
     * ePrint 2026/2205 specifies the proof relation but does not instantiate
     * its simulation-extractable NIZK. This experiment therefore excludes the
     * proof from both the implementation and the performance comparison.
     */
    struct Ciphertext
    {
        blst::G1Point c1 = blst::G1Point::infinity();
        blst::G2Point c2 = blst::G2Point::infinity();
        blst::G1Point c3 = blst::G1Point::infinity();
        blst::GT c4 = blst::GT::identity();
    };

    struct PartialDecryptionShare
    {
        std::uint32_t index = 0;
        blst::G1Point value = blst::G1Point::infinity();
    };

    /**
     * Ciphertext-only work that may run while decryption shares are in flight.
     * `cross_terms[i]` is the last two pairing terms in Eq. (1), including
     * their signs. Final decryption therefore needs only one pairing per slot.
     */
    struct BatchPrecomputation
    {
        std::vector<blst::GT> cross_terms;
        std::vector<blst::G2AffinePoint> c2_affine;
    };
}
