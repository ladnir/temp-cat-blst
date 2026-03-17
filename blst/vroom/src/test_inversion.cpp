#include "../test_data/test_values.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include "bounded_ring.hpp"
#include "inversion.hpp"
#include <iostream>
#include <iomanip>
#include <random>
#include <cassert>

// BLS12-381 modulus q
const char* bls12_381_modulus_hex = "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab";

void test_inversion() {
    std::cout << "=== Testing BLS381AddChainInversion ===" << std::endl;
    
    // Create BLS12-381 modulus
    BigInt q(bls12_381_modulus_hex, 16);
    std::cout << "BLS12-381 modulus q: " << q.to_string() << std::endl;
    std::cout << "Modulus bit length: " << q.bit_length() << " bits" << std::endl << std::endl;
    
    // Create BoundedRing
    BoundedRing<381, 8, 52, -1932, 2377> ring(q);
    std::cout << "Created BoundedRing" << std::endl << std::endl;
    
    // Create inversion object
    BLS381AddChainInversion<BoundedRing<381, 8, 52, -1932, 2377>> inverter(ring);
    
    // Test with random field elements
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(1, UINT64_MAX);
    
    const int num_tests = 10;
    int passed = 0;
    int failed = 0;
    
    std::cout << "Testing " << num_tests << " random field elements..." << std::endl << std::endl;
    
    for (int i = 0; i < num_tests; i++) {
        // Generate random field element
        BigInt random_val = BigInt::random(q);
        
        // Skip zero (no inverse)
        if (random_val == BigInt(0)) {
            i--;
            continue;
        }
        
        // Convert to ring element
        auto element = ring.from_bigint(random_val);
        
        // Compute inverse using addition chain
        auto inverted_element = inverter.invert(element, ring);
        
        // Convert back to BigInt
        BigInt inverted_val = ring.to_bigint(inverted_element);
        
        // Compute expected inverse using BigInt mod_inverse
        BigInt expected_inverse = random_val.mod_inverse(q);
        
        // Verify: (inverted_val * random_val) % q == 1
        BigInt product = (inverted_val * random_val) % q;
        if (product < 0) {
            product = product + q;
        }
        
        // Also check if inverted_val matches expected_inverse (mod q)
        BigInt diff = (inverted_val - expected_inverse) % q;
        if (diff < 0) {
            diff = diff + q;
        }
        bool matches_expected = (diff == BigInt(0) || diff == q);
        
        bool test_passed = (product == BigInt(1)) && matches_expected;
        
        if (test_passed) {
            passed++;
            std::cout << "✓ Test " << (i + 1) << " passed" << std::endl;
        } else {
            failed++;
            std::cout << "✗ Test " << (i + 1) << " failed" << std::endl;
            std::cout << "  Input:        " << random_val.to_string() << std::endl;
            std::cout << "  Computed inv: " << inverted_val.to_string() << std::endl;
            std::cout << "  Expected inv: " << expected_inverse.to_string() << std::endl;
            std::cout << "  Product:      " << product.to_string() << " (should be 1)" << std::endl;
            std::cout << "  Difference:   " << diff.to_string() << std::endl;
        }
    }
    
    std::cout << std::endl;
    std::cout << "=== Test Summary ===" << std::endl;
    std::cout << "Passed: " << passed << "/" << num_tests << std::endl;
    std::cout << "Failed: " << failed << "/" << num_tests << std::endl;
    
    if (failed == 0) {
        std::cout << "✓ All tests passed!" << std::endl;
    } else {
        std::cout << "✗ Some tests failed!" << std::endl;
    }
}

int main() {
    try {
        test_inversion();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

