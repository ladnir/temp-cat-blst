#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "blst/fields/Fr.h"
using namespace blst;

#define FR_LIMBS NLIMBS(256)

// Helper: Set a Fr& to a small integer value
void set_fr(Fr& v, uint64_t x) {
    for (int i = 0; i < FR_LIMBS; ++i)
        v[i] = 0;
    memcpy(v.l, &x, sizeof(x));
}

// Helper: Compare two Fr& for equality
bool fr_equal(const Fr& a, const Fr& b) {
    for (int i = 0; i < FR_LIMBS; ++i)
        if (a[i] != b[i]) return false;
    return true;
}

// Helper: Print Fr& for debugging
void print_fr(const Fr& v) {
    for (int i = FR_LIMBS-1; i >= 0; --i)
        printf("%016llx%s", (unsigned long long)v[i], i ? " " : "\n");
}

// Helper: Convert standard integer to Montgomery form
void to_mont_fr(Fr& r, const Fr& a) {
    // Use *get_get_BLST_r()()RR for conversion to Montgomery form
    Fr::to_mont(r, a);
}

// Helper: Copy one Fr& to another
void fr_copy(Fr& dst, const Fr& src) {
    for (int i = 0; i < FR_LIMBS; ++i)
        dst[i] = src[i];
}

TEST(FrOps, Add) {
    Fr a, b, result, expected;
    
    // Basic addition: 5 + 7 = 12
    set_fr(a, 5);
    set_fr(b, 7);
    Fr::add(result, a, b);
    set_fr(expected, 12);
    EXPECT_TRUE(fr_equal(result, expected)) << "5 + 7 should equal 12";
    
    // Addition with zero
    set_fr(a, 42);
    set_fr(b, 0);
    Fr::add(result, a, b);
    EXPECT_TRUE(fr_equal(result, a)) << "42 + 0 should equal 42";
    
    // Test commutativity: a + b = b + a
    set_fr(a, 123);
    set_fr(b, 456);
    Fr::add(result, a, b);
    Fr result2;
    Fr::add(result2, b, a);
    EXPECT_TRUE(fr_equal(result, result2)) << "Addition should be commutative";
}

TEST(FrOps, AddOverflow) {
    Fr a, b, result, expected;
    
    // Set a = r - 1 (order minus 1)
    fr_copy(a, get_BLST_r());
    Fr one;
    set_fr(one, 1);
    Fr::sub(a, a, one);
    
    // Add 2 to get 1 (modulo r)
    set_fr(b, 2);
    Fr::add(result, a, b);
    set_fr(expected, 1);
    
    EXPECT_TRUE(fr_equal(result, expected)) << "Overflow in addition should wrap modulo r";
}

TEST(FrOps, Sub) {
    Fr a, b, result, expected;
    
    // Basic subtraction: 10 - 3 = 7
    set_fr(a, 10);
    set_fr(b, 3);
    Fr::sub(result, a, b);
    set_fr(expected, 7);
    EXPECT_TRUE(fr_equal(result, expected)) << "10 - 3 should equal 7";
    
    // Subtraction with zero
    set_fr(a, 42);
    set_fr(b, 0);
    Fr::sub(result, a, b);
    EXPECT_TRUE(fr_equal(result, a)) << "42 - 0 should equal 42";
}

TEST(FrOps, SubUnderflow) {
    Fr a, b, result;
    
    // Underflow: 3 - 10 = r - 7
    set_fr(a, 3);
    set_fr(b, 10);
    Fr::sub(result, a, b);
    
    // Expected: r - 7
    Fr expected;
    fr_copy(expected, get_BLST_r());
    Fr seven;
    set_fr(seven, 7);
    Fr::sub(expected, expected, seven);
    
    EXPECT_TRUE(fr_equal(result, expected)) << "Underflow should wrap modulo r";
}

TEST(FrOps, MulBy3) {
    Fr a, result, expected;
    
    // 2 * 3 = 6
    set_fr(a, 2);
    Fr::mul_by_3(result, a);
    set_fr(expected, 6);
    EXPECT_TRUE(fr_equal(result, expected)) << "2 * 3 should equal 6";
    
    // 0 * 3 = 0
    set_fr(a, 0);
    Fr::mul_by_3(result, a);
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result, expected)) << "0 * 3 should equal 0";
    
    // Test with larger value
    set_fr(a, 0x123456789ABCDEF0ULL);
    Fr::mul_by_3(result, a);
    
    // Verify by multiplying by 3 using addition
    Fr temp;
    Fr::add(temp, a, a);    // temp = 2*a
    Fr::add(expected, temp, a);  // expected = 3*a
    
    EXPECT_TRUE(fr_equal(result, expected)) << "MulBy3 should match manual calculation";
}

TEST(FrOps, LShift) {
    Fr a, result, expected;

    set_fr(a, 0x123456789ABCDEF0ULL);
    Fr::lshift(result, a, 0);
    EXPECT_TRUE(fr_equal(result, a)) << "zero-bit left shift should be identity";
    Fr::lshift(a, a, 0);
    EXPECT_TRUE(fr_equal(result, a)) << "in-place zero-bit left shift should be identity";
    
    // 1 << 2 = 4
    set_fr(a, 1);
    Fr::lshift(result, a, 2);
    set_fr(expected, 4);
    EXPECT_TRUE(fr_equal(result, expected)) << "1 << 2 should equal 4";
    
    // 5 << 1 = 10
    set_fr(a, 5);
    Fr::lshift(result, a, 1);
    set_fr(expected, 10);
    EXPECT_TRUE(fr_equal(result, expected)) << "5 << 1 should equal 10";
    
    // 0 << any = 0
    set_fr(a, 0);
    Fr::lshift(result, a, 10);
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result, expected)) << "0 << 10 should equal 0";
}

TEST(FrOps, RShift) {
    Fr a, result, expected;

    set_fr(a, 0x123456789ABCDEF0ULL);
    Fr::rshift(result, a, 0);
    EXPECT_TRUE(fr_equal(result, a)) << "zero-bit right shift should be identity";
    Fr::rshift(a, a, 0);
    EXPECT_TRUE(fr_equal(result, a)) << "in-place zero-bit right shift should be identity";
    
    // 8 >> 2 = 2
    set_fr(a, 8);
    Fr::rshift(result, a, 2);
    set_fr(expected, 2);
    EXPECT_TRUE(fr_equal(result, expected)) << "8 >> 2 should equal 2";
    
    // 10 >> 1 = 5
    set_fr(a, 10);
    Fr::rshift(result, a, 1);
    set_fr(expected, 5);
    EXPECT_TRUE(fr_equal(result, expected)) << "10 >> 1 should equal 5";
    
    // 0 >> any = 0
    set_fr(a, 0);
    Fr::rshift(result, a, 10);
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result, expected)) << "0 >> 10 should equal 0";
}

TEST(FrOps, Mul) {
    Fr a, b, result, result_std, expected;
    
    // Convert to Montgomery form, multiply, and convert back
    set_fr(a, 3);
    set_fr(b, 4);
    
    to_mont_fr(a, a);
    to_mont_fr(b, b);
    
    Fr::mul(result, a, b);
    Fr::from_mont(result_std, result);
    
    set_fr(expected, 12);
    EXPECT_TRUE(fr_equal(result_std, expected)) << "3 * 4 should equal 12";
    
    // Test multiplication by zero
    set_fr(a, 42);
    set_fr(b, 0);
    to_mont_fr(a, a);
    to_mont_fr(b, b);
    
    Fr::mul(result, a, b);
    Fr::from_mont(result_std, result);
    
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result_std, expected)) << "42 * 0 should equal 0";
    
    // Test commutativity: a * b = b * a
    set_fr(a, 123);
    set_fr(b, 456);
    to_mont_fr(a, a);
    to_mont_fr(b, b);
    
    Fr::mul(result, a, b);
    Fr result2;
    Fr::mul(result2, b, a);
    
    EXPECT_TRUE(fr_equal(result, result2)) << "Multiplication should be commutative";
}

TEST(FrOps, Sqr) {
    Fr a, result, result_std, expected;
    
    // 5² = 25
    set_fr(a, 5);
    to_mont_fr(a, a);
    
    Fr::sqr(result, a);
    Fr::from_mont(result_std, result);
    
    set_fr(expected, 25);
    EXPECT_TRUE(fr_equal(result_std, expected)) << "5² should equal 25";
    
    // 0² = 0
    set_fr(a, 0);
    to_mont_fr(a, a);
    
    Fr::sqr(result, a);
    Fr::from_mont(result_std, result);
    
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result_std, expected)) << "0² should equal 0";
    
    // Compare sqr with mul by self
    set_fr(a, 123);
    to_mont_fr(a, a);
    
    Fr::sqr(result, a);
    Fr result2;
    Fr::mul(result2, a, a);
    
    EXPECT_TRUE(fr_equal(result, result2)) << "a² should equal a*a";
}

TEST(FrOps, ExpZeroIsOne) {
    Fr base;
    Fr result;

    set_fr(base, 5);
    Fr::to_mont(base, base);
    Fr::exp(result, base, nullptr, 0);
    EXPECT_TRUE(fr_equal(result, Fr::one()));

    Fr::exp(base, base, nullptr, 0);
    EXPECT_TRUE(fr_equal(base, Fr::one()));
}

TEST(FrOps, CNeg) {
    Fr a, result, expected;
    
    // Negate 5 when flag is true
    set_fr(a, 5);
    Fr::cneg(result, a, 1);
    
    // Expected: r - 5
    fr_copy(expected, get_BLST_r());
    Fr five;
    set_fr(five, 5);
    Fr::sub(expected, expected, five);
    
    EXPECT_TRUE(fr_equal(result, expected)) << "Negation of 5 should equal r-5";
    
    // Don't negate when flag is false
    set_fr(a, 5);
    Fr::cneg(result, a, 0);
    EXPECT_TRUE(fr_equal(result, a)) << "cneg with flag=0 should not change value";
    
    // Negate 0 (should remain 0)
    set_fr(a, 0);
    Fr::cneg(result, a, 1);
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result, expected)) << "Negation of 0 should equal 0";
}

TEST(FrOps, MontgomeryConversion) {
    Fr a, mont_a, back_a;
    
    // Convert 1 to Montgomery form and back
    set_fr(a, 1);
    Fr::to_mont(mont_a, a);
    Fr::from_mont(back_a, mont_a);
    
    EXPECT_TRUE(fr_equal(back_a, a)) << "Montgomery conversion of 1 failed";
    
    // Test with various values
    const uint64_t TEST_VALUES[] = {0, 2, 42, 12345, 0xFFFFFFFFULL, 0xDEADBEEFULL};
    const int NUM_VALUES = sizeof(TEST_VALUES) / sizeof(TEST_VALUES[0]);
    
    for (int i = 0; i < NUM_VALUES; i++) {
        set_fr(a, TEST_VALUES[i]);
        
        // To Montgomery form and back
        Fr::to_mont(mont_a, a);
        Fr::from_mont(back_a, mont_a);
        
        if (!fr_equal(back_a, a)) {
            printf("Failed for value %llu (0x%llx)\n", 
                  (unsigned long long)TEST_VALUES[i], 
                  (unsigned long long)TEST_VALUES[i]);
            printf("Original: "); print_fr(a);
            printf("Mont form: "); print_fr(mont_a);
            printf("Back: "); print_fr(back_a);
        }
        
        EXPECT_TRUE(fr_equal(back_a, a)) 
            << "Montgomery conversion failed for " << TEST_VALUES[i];
    }
}

TEST(FrOps, ArithmeticProperties) {
    Fr a, b, c, result1, result2, mont_a, mont_b, mont_c;
    
    // Test values
    set_fr(a, 17);
    set_fr(b, 23);
    set_fr(c, 42);
    
    // Convert to Montgomery form
    Fr::to_mont(mont_a, a);
    Fr::to_mont(mont_b, b);
    Fr::to_mont(mont_c, c);
    
    // Associativity: (a + b) + c = a + (b + c)
    Fr temp1, temp2;// , std_result1, std_result2;
    
    Fr::add(temp1, a, b);
    Fr::add(result1, temp1, c);
    
    Fr::add(temp2, b, c);
    Fr::add(result2, a, temp2);
    
    EXPECT_TRUE(fr_equal(result1, result2)) << "Addition should be associative";
    
    // Distributivity: a * (b + c) = a*b + a*c
    Fr::add(temp1, mont_b, mont_c);
    Fr::mul(result1, mont_a, temp1);
    
    Fr::mul(temp1, mont_a, mont_b);
    Fr::mul(temp2, mont_a, mont_c);
    Fr::add(result2, temp1, temp2);
    
    EXPECT_TRUE(fr_equal(result1, result2)) << "Multiplication should distribute over addition";
}

TEST(FrOps, EdgeCases) {
    Fr zero, one, a, result;// , expected;
    
    // Test with zero and one
    set_fr(zero, 0);
    set_fr(one, 1);
    
    // 0 + 0 = 0
    Fr::add(result, zero, zero);
    EXPECT_TRUE(fr_equal(result, zero)) << "0 + 0 should equal 0";
    
    // 0 * any = 0
    set_fr(a, 42);
    to_mont_fr(a, a);
    to_mont_fr(zero, zero);
    Fr::mul(result, a, zero);
    Fr::from_mont(result, result);
    EXPECT_TRUE(fr_equal(result, zero)) << "42 * 0 should equal 0";
    
    // 1 * any = any
    set_fr(a, 42);
    to_mont_fr(a, a);
    to_mont_fr(one, one);
    Fr::mul(result, a, one);
    Fr std_a;
    Fr::from_mont(std_a, a);
    Fr::from_mont(result, result);
    EXPECT_TRUE(fr_equal(result, std_a)) << "42 * 1 should equal 42";
}

// Test butterfly operations (used in FFT)
TEST(FrOps, Butterfly) {
    Fr x0, x1, twiddle, expected_x0, expected_x1;
    Fr original_x0, original_x1;
    
    // Setup simple test values
    set_fr(x0, 10);
    set_fr(x1, 20);
    set_fr(twiddle, 2);
    
    // Convert to Montgomery form
    to_mont_fr(x0, x0);
    to_mont_fr(x1, x1);
    to_mont_fr(twiddle, twiddle);
    
    // Save original values
    fr_copy(original_x0, x0);
    fr_copy(original_x1, x1);
    
    // CT butterfly: x0' = x0 + x1*twiddle, x1' = x0 - x1*twiddle
    Fr::ct_bfly(x0, x1, twiddle);
    
    // Calculate expected values manually
    Fr temp;
    Fr::mul(temp, original_x1, twiddle);  // temp = x1 * twiddle
    Fr::add(expected_x0, original_x0, temp);  // x0' = x0 + temp
    Fr::sub(expected_x1, original_x0, temp);  // x1' = x0 - temp
    
    EXPECT_TRUE(fr_equal(x0, expected_x0)) << "CT butterfly x0 calculation incorrect";
    EXPECT_TRUE(fr_equal(x1, expected_x1)) << "CT butterfly x1 calculation incorrect";
    
    // Reset for GS butterfly
    fr_copy(x0, original_x0);
    fr_copy(x1, original_x1);
    
    // GS butterfly: x0' = x0 + x1, x1' = (x0 - x1) * twiddle
    Fr::gs_bfly(x0, x1, twiddle);
    
    // Calculate expected values manually
    Fr::add(expected_x0, original_x0, original_x1);  // x0' = x0 + x1
    Fr::sub(temp, original_x0, original_x1);         // temp = x0 - x1
    Fr::mul(expected_x1, temp, twiddle);             // x1' = temp * twiddle
    
    EXPECT_TRUE(fr_equal(x0, expected_x0)) << "GS butterfly x0 calculation incorrect";
    EXPECT_TRUE(fr_equal(x1, expected_x1)) << "GS butterfly x1 calculation incorrect";
}
