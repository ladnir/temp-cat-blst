#include "../cpu/precompute/gmp_wrapper.hpp"
extern "C" {
#include "../blst/vect.h"
#include "../blst/fields.h"
#include "../blst/consts.h"
#include "../blst/point.h"
}
#include "pairing.hpp"
#include "bounded_ring.hpp"
#include "fp2.hpp"
#include <iostream>
#include <iomanip>
#include <random>
#include <cstring>

// BLS12-381 modulus q
const char* bls12_381_modulus_hex = "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab";

// External BLST functions
extern "C" {
    extern const POINTonE1 BLS12_381_G1;
    extern const POINTonE2 BLS12_381_G2;
    void blst_p1_mult(POINTonE1 *out, const POINTonE1 *a, const byte *scalar, size_t nbits);
    void blst_p1_to_affine(POINTonE1_affine *out, const POINTonE1 *a);
    void blst_p2_mult(POINTonE2 *out, const POINTonE2 *a, const byte *scalar, size_t nbits);
    void blst_p2_to_affine(POINTonE2_affine *out, const POINTonE2 *a);
    void blst_miller_loop(vec384fp12 ret, const POINTonE2_affine *Q, const POINTonE1_affine *P);
    void blst_final_exp(vec384fp12 ret, const vec384fp12 f);
}

// Helper function to convert vec384 (Montgomery form) to BigInt
// Note: This is a local version for the test file
static BigInt vec384_montgomery_to_bigint_local(const vec384 a) {
    vec384 normal;
    from_fp(normal, a);
    
    BigInt result(0);
    BigInt two_to_64 = BigInt(1);
    two_to_64 = two_to_64 << 64;
    for (int i = 5; i >= 0; i--) {
        result = result * two_to_64 + BigInt(static_cast<unsigned long>(normal[i]));
    }
    return result;
}

// Helper function to convert vec384x (Fp2, Montgomery form) to BigInt pair
std::pair<BigInt, BigInt> vec384x_montgomery_to_bigint(const vec384x a) {
    return std::make_pair(
        vec384_montgomery_to_bigint_local(a[0]),
        vec384_montgomery_to_bigint_local(a[1])
    );
}

// Helper function to convert vec384fp12 (Montgomery form) to array of BigInt
std::array<BigInt, 12> vec384fp12_montgomery_to_bigint(const vec384fp12 a) {
    std::array<BigInt, 12> result;
    // vec384fp12 is vec384fp6[2], vec384fp6 is vec384x[3], vec384x is vec384[2]
    // For index i (0-11): a[i/6][(i%6)/2][i%2]
    for (int i = 0; i < 12; i++) {
        result[i] = vec384_montgomery_to_bigint_local(a[i/6][(i%6)/2][i%2]);
    }
    return result;
}

// Helper function to convert POINTonE1_affine to our AffinePoint
template<class Ring>
AffinePoint<typename Ring::StandardElement> blst_g1_to_affine_point(
    const POINTonE1_affine &p,
    const Ring &ring
) {
    BigInt x = vec384_montgomery_to_bigint_local(p.X);
    BigInt y = vec384_montgomery_to_bigint_local(p.Y);
    
    return AffinePoint<typename Ring::StandardElement>(
        ring.from_bigint(x),
        ring.from_bigint(y)
    );
}

// Helper function to convert POINTonE2_affine to our AffinePoint
template<class Ring>
AffinePoint<typename FP2Ring<Ring>::StandardElement> blst_g2_to_affine_point(
    const POINTonE2_affine &p,
    const Ring &ring
) {
    FP2Ring<Ring> fp2ring(ring);
    auto x_pair = vec384x_montgomery_to_bigint(p.X);
    auto y_pair = vec384x_montgomery_to_bigint(p.Y);
    
    return AffinePoint<typename FP2Ring<Ring>::StandardElement>(
        typename FP2Ring<Ring>::StandardElement(
            ring.from_bigint(x_pair.first),
            ring.from_bigint(x_pair.second)
        ),
        typename FP2Ring<Ring>::StandardElement(
            ring.from_bigint(y_pair.first),
            ring.from_bigint(y_pair.second)
        )
    );
}

// Helper function to compare FP12 results
template<class Ring>
bool compare_fp12_results(
    const std::array<BigInt, 12> &blst_result,
    const typename FP12Ring<Ring>::StandardElement &our_result,
    const Ring &ring,
    const std::string &test_name
) {
    bool all_match = true;
    BigInt q(bls12_381_modulus_hex, 16);
    
    for (int i = 0; i < 12; i++) {
        BigInt our_val = ring.to_bigint(our_result[i]);
        BigInt blst_val = blst_result[i];
        
        BigInt our_mod = our_val % q;
        if (our_mod < 0) our_mod = our_mod + q;
        BigInt blst_mod = blst_val % q;
        if (blst_mod < 0) blst_mod = blst_mod + q;
        
        if (our_mod != blst_mod) {
            all_match = false;
            std::cout << "✗ Mismatch at index " << i << " in " << test_name << std::endl;
            std::cout << "  BLST:  " << blst_mod.to_string(16) << std::endl;
            std::cout << "  Ours:  " << our_mod.to_string(16) << std::endl;
        }
    }
    
    return all_match;
}

template<int modulus_bits=381, int limbs=8, int element_bits=52, int RNS_MAX_NEG=-1932, int RNS_MAX_POS=2377>
void test_pairing_config(const std::string &config_name) {
    std::cout << "=== Testing Pairing (" << config_name << "): BLST vs Our Implementation ===" << std::endl;
    
    BigInt q(bls12_381_modulus_hex, 16);
    BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS> ring(q);
    FP2Ring<BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS>> fp2ring(ring);
    FP12Ring<BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS>> fp12ring(ring);
    
    using RingType = BoundedRing<modulus_bits, limbs, element_bits, RNS_MAX_NEG, RNS_MAX_POS>;
    using StandardElement = typename RingType::StandardElement;
    using FP2StandardElement = typename FP2Ring<RingType>::StandardElement;
    
    int passed = 0;
    int failed = 0;
    const int num_tests = 5;
    
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    std::cout << "Testing " << num_tests << " random pairings..." << std::endl << std::endl;
    
    for (int t = 0; t < num_tests; t++) {
        // Generate random 256-bit scalars
        byte scalar1[32], scalar2[32];
        for (int i = 0; i < 32; i++) {
            scalar1[i] = static_cast<byte>(gen() & 0xff);
            scalar2[i] = static_cast<byte>(gen() & 0xff);
        }
        
        // Compute G1 point: scalar1 * G1
        POINTonE1 p1_proj;
        blst_p1_mult(&p1_proj, &BLS12_381_G1, scalar1, 256);
        POINTonE1_affine p1_affine;
        blst_p1_to_affine(&p1_affine, &p1_proj);
        
        // Compute G2 point: scalar2 * G2
        POINTonE2 p2_proj;
        blst_p2_mult(&p2_proj, &BLS12_381_G2, scalar2, 256);
        POINTonE2_affine p2_affine;
        blst_p2_to_affine(&p2_affine, &p2_proj);
        
        // Compute pairing using BLST
        vec384fp12 blst_miller_result, blst_final_result;
        blst_miller_loop(blst_miller_result, &p2_affine, &p1_affine);
        blst_final_exp(blst_final_result, blst_miller_result);
        
        // Convert BLST result to BigInt array
        std::array<BigInt, 12> blst_result = vec384fp12_montgomery_to_bigint(blst_final_result);
        
        // Convert BLST points to our format
        auto P = blst_g1_to_affine_point(p1_affine, ring);
        auto Q = blst_g2_to_affine_point(p2_affine, ring);
        
        // Compute pairing using our implementation
        auto our_result = pairing_with_blst_inverter(P, Q, fp2ring, fp12ring, ring);
        
        // Compare results
        bool match = compare_fp12_results(blst_result, our_result, ring, 
                                         "Test " + std::to_string(t + 1));
        
        if (match) {
            passed++;
            std::cout << "✓ Test " << (t + 1) << " passed" << std::endl;
        } else {
            failed++;
            std::cout << "✗ Test " << (t + 1) << " failed" << std::endl;
        }
        std::cout << std::endl;
    }
    
    std::cout << "=== Pairing Test Summary (" << config_name << ") ===" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    
    if (failed == 0) {
        std::cout << "✓ All pairing tests passed!" << std::endl;
    } else {
        std::cout << "✗ Some pairing tests failed!" << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    try {
        // Test with 52-bit configuration
        test_pairing_config<381, 8, 52, -1932, 2377>("52-bit");
        
        // Test with 50-bit configuration
        test_pairing_config<381, 8, 50>("50-bit");
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

