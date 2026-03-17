#include "elementwise_reduce.hpp"
#include "ring_element.hpp"
#include "bounded_element.hpp"
#include "../cpu/vector/vector_impl.hpp"
#include <iostream>
#include <random>
#include <cassert>
#include <cstdint>
#include <array>

typedef __int128_t int128_t;

template<int limbs, int element_bits>
void test_elementwise_reduce() {
    std::cout << "\n=== Testing Elementwise Reduce with limbs=" << limbs << ", element_bits=" << element_bits << " ===\n" << std::endl;
    
    // Use moduli t = 1, 3, 5, 7, 11, 13, 17, 19
    std::array<uint64_t, limbs> t_values = {1, 3, 5, 7, 11, 13, 17, 19};
    std::array<uint64_t, limbs> moduli;
    for (int i = 0; i < limbs; i++) {
        moduli[i] = (1ull << element_bits) - t_values[i];
    }
    
    MontgomeryReduce2<limbs, element_bits> reducer(moduli);
    
    // Test 1: to_mont and from_mont are inverses
    std::cout << "Test 1: to_mont and from_mont are inverses" << std::endl;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(0, (1ULL << element_bits) - 1);
    
    for (int test = 0; test < 10; test++) {
        // Create an array with the same value for all limbs
        std::array<uint64_t, limbs> x_arr;
        uint64_t x = dis(gen);
        for (int i = 0; i < limbs; i++) {
            x_arr[i] = x % moduli[i];
        }
        
        // Convert to AVXVector
        AVXVector<limbs> x_vec;
        x_vec.load(x_arr.data());
        
        // Convert to Montgomery form
        auto mont_x = reducer.to_mont(x_vec);
        
        // Convert back from Montgomery form
        auto recovered_vec = reducer.from_mont(mont_x);
        auto recovered_arr = recovered_vec.to_array();
        
        // Check if they match
        for (int i = 0; i < limbs; i++) {
            uint64_t expected = x_arr[i];
            if (recovered_arr[i] != expected) {
                std::cerr << "ERROR: to_mont/from_mont mismatch at limb " << i 
                          << ": input=" << x_arr[i] << ", recovered=" << recovered_arr[i] 
                          << ", expected=" << expected << std::endl;
                assert(false);
            }
        }
    }
    std::cout << "  ✓ to_mont and from_mont are inverses" << std::endl;
    
    // Test 2: min_sub_reduce and mask_reduce for elementwise arithmetic (x + y + z)
    std::cout << "\nTest 2: min_sub_reduce and mask_reduce for elementwise arithmetic" << std::endl;
    
    // Create random values
    std::array<uint64_t, limbs> arr_x, arr_y, arr_z;
    for (int i = 0; i < limbs; i++) {
        arr_x[i] = dis(gen);
        arr_y[i] = dis(gen);
        arr_z[i] = dis(gen);
    }
    
    auto x = from_constant<limbs, element_bits>(arr_x, arr_x);
    auto y = from_constant<limbs, element_bits>(arr_y, arr_y);
    auto z = from_constant<limbs, element_bits>(arr_z, arr_z);
    
    // x + y + z (already complete, no DelayedProduct)
    auto sum = x + y + z;
    
    // Test mask_reduce on the low part
    auto mask_reduced = reducer.mask_reduce(sum.m1);
    auto mask_reduced_array = mask_reduced.data.to_array();
    
    // Test min_sub_reduce (reduce to bound 1)
    auto min_sub_reduced = reducer.template min_sub_reduce<2, 1>(mask_reduced);
    auto min_sub_reduced_array = min_sub_reduced.data.to_array();
    
    // Compute expected: (x + y + z) mod p for each limb
    for (int i = 0; i < limbs; i++) {
        uint64_t expected_sum = (arr_x[i] + arr_y[i] + arr_z[i]) % moduli[i];
        uint64_t expected_mask = expected_sum % (1ULL << element_bits);
        uint64_t expected_min_sub = expected_mask % moduli[i];
        
        // mask_reduce: takes high bits, multiplies by t, adds to low bits
        // This is a reduction step, so we check it's in the right range
        if (mask_reduced_array[i] >= (1ULL << (element_bits + 1))) {
            std::cerr << "ERROR: mask_reduce[" << i << "] = " << mask_reduced_array[i] 
                      << " exceeds expected range" << std::endl;
            assert(false);
        }
        
        // min_sub_reduce should give us the final reduced value
        if (min_sub_reduced_array[i] != expected_min_sub) {
            std::cerr << "ERROR: min_sub_reduce[" << i << "] = " << min_sub_reduced_array[i] 
                      << ", expected " << expected_min_sub << std::endl;
            assert(false);
        }
    }
    std::cout << "  ✓ min_sub_reduce and mask_reduce work correctly" << std::endl;
    
    // Test 3: reduce_small for a product a*b (after (a*b).complete())
    std::cout << "\nTest 3: reduce_small for product a*b" << std::endl;
    
    std::array<uint64_t, limbs> arr_a, arr_b;
    for (int i = 0; i < limbs; i++) {
        arr_a[i] = dis(gen);
        arr_b[i] = dis(gen);
    }
    
    // Convert to Montgomery form
    // Create arrays with Montgomery form values
    std::array<uint64_t, limbs> mont_arr_a, mont_arr_b;
    for (int i = 0; i < limbs; i++) {
        AVXVector<limbs> a_vec, b_vec;
        std::array<uint64_t, limbs> a_single, b_single;
        for (int j = 0; j < limbs; j++) {
            a_single[j] = (j == i) ? arr_a[i] : 0;
            b_single[j] = (j == i) ? arr_b[i] : 0;
        }
        a_vec.load(a_single.data());
        b_vec.load(b_single.data());
        auto mont_a = reducer.to_mont(a_vec);
        auto mont_b = reducer.to_mont(b_vec);
        auto mont_a_arr = mont_a.data.to_array();
        auto mont_b_arr = mont_b.data.to_array();
        mont_arr_a[i] = mont_a_arr[i];
        mont_arr_b[i] = mont_b_arr[i];
    }
    
    auto a = from_constant<limbs, element_bits>(mont_arr_a, mont_arr_a);
    auto b = from_constant<limbs, element_bits>(mont_arr_b, mont_arr_b);
    auto ab = a * b;
    auto ab_complete = ab.complete();
    
    // reduce_wide expects WideElement with specific bounds
    auto reduced = reducer.reduce_wide(ab_complete.m1);
    auto reduced_array = reduced.data.to_array();
    
    // The result is in Montgomery form, convert back to standard form to verify
    // Reduce bounds to <0, 1> before calling from_mont
    auto reduced_normalized = reducer.template min_sub_reduce<2, 1>(reduced);
    auto recovered_small = reducer.from_mont(reduced_normalized);
    auto recovered_small_arr = recovered_small.to_array();
    
    // Compute expected: (a * b) mod p in standard form
    for (int i = 0; i < limbs; i++) {
        uint64_t expected = ((uint128_t)arr_a[i] * (uint128_t)arr_b[i]) % moduli[i];
        
        if (recovered_small_arr[i] != expected) {
            std::cerr << "ERROR: reduce_small[" << i << "] recovered=" << recovered_small_arr[i] 
                      << ", expected " << expected << " (montgomery=" << reduced_array[i] << ")" << std::endl;
            assert(false);
        }
    }
    std::cout << "  ✓ reduce_small works correctly for product" << std::endl;
    
    // Test 4: reduce_full for sum of products ((a*b).complete()+(c*d))
    std::cout << "\nTest 4: reduce_full for sum of products" << std::endl;
    
    std::array<uint64_t, limbs> arr_c, arr_d;
    for (int i = 0; i < limbs; i++) {
        arr_c[i] = dis(gen);
        arr_d[i] = dis(gen);
    }
    
    // Convert to Montgomery form
    std::array<uint64_t, limbs> mont_arr_c, mont_arr_d;
    for (int i = 0; i < limbs; i++) {
        AVXVector<limbs> c_vec, d_vec;
        std::array<uint64_t, limbs> c_single, d_single;
        for (int j = 0; j < limbs; j++) {
            c_single[j] = (j == i) ? arr_c[i] : 0;
            d_single[j] = (j == i) ? arr_d[i] : 0;
        }
        c_vec.load(c_single.data());
        d_vec.load(d_single.data());
        auto mont_c = reducer.to_mont(c_vec);
        auto mont_d = reducer.to_mont(d_vec);
        auto mont_c_arr = mont_c.data.to_array();
        auto mont_d_arr = mont_d.data.to_array();
        mont_arr_c[i] = mont_c_arr[i];
        mont_arr_d[i] = mont_d_arr[i];
    }
    
    auto c = from_constant<limbs, element_bits>(mont_arr_c, mont_arr_c);
    auto d = from_constant<limbs, element_bits>(mont_arr_d, mont_arr_d);
    auto cd = c * d;
    auto cd_complete = cd.complete();
    
    // (a*b).complete() + (c*d).complete()
    auto sum_products = ab_complete + cd_complete;
    
    // reduce_auto reduces the wide element
    auto reduced_full = reducer.template reduce_auto<1>(sum_products.m1);
    auto reduced_full_array = reduced_full.data.to_array();
    
    // The wide value is already in Montgomery form (since a, b, c, d are in Montgomery form)
    // reduce_full reduces it and keeps it in Montgomery form
    // To verify, we convert back from Montgomery form and check it matches the expected standard form
    // Reduce bounds to <0, 1> before calling from_mont
    auto reduced_full_normalized = reducer.min_sub_reduce(reduced_full);
    auto recovered_full_elem = reducer.from_mont(reduced_full_normalized);
    auto recovered_full = recovered_full_elem.to_array();
    
    // Compute expected: (a*b + c*d) mod p in standard form
    for (int i = 0; i < limbs; i++) {
        uint64_t p = moduli[i];
        uint64_t expected = (((uint128_t)arr_a[i] * (uint128_t)arr_b[i]) + 
                            ((uint128_t)arr_c[i] * (uint128_t)arr_d[i])) % p;
        
        // For 50-bit, the to_array() representation might not directly correspond to the mathematical value
        // because it uses 52-bit shifts for hardware representation
        // But reduce_full should still produce the correct Montgomery form of the reduced value
        if (recovered_full[i] != expected) {
            std::cerr << "ERROR: reduce_full[" << i << "] recovered=" << recovered_full[i] 
                      << ", expected " << expected << " (montgomery=" << reduced_full_array[i] << ")" << std::endl;
            // Check if it's a multiple of p off (might be a representation issue)
            int64_t diff = (int64_t)recovered_full[i] - (int64_t)expected;
            if (diff % (int64_t)p == 0) {
                std::cerr << "  Note: Difference is a multiple of p: " << diff / p << " * p" << std::endl;
            }
            assert(false);
        }
    }
    std::cout << "  ✓ reduce_full works correctly for sum of products" << std::endl;
}

int main() {
    std::cout << "Testing Elementwise Reduce operations" << std::endl;
    std::cout << "Using -DUSE_AVX512_FALLBACK" << std::endl;
    std::cout << "Moduli: t = 1, 3, 5, 7, 11, 13, 17, 19" << std::endl;
    
    // Test with limbs=8, element_bits=52
    test_elementwise_reduce<8, 52>();
    
    // Test with limbs=8, element_bits=50
    test_elementwise_reduce<8, 50>();
    
    std::cout << "\n=== All tests passed! ===" << std::endl;
    return 0;
}

