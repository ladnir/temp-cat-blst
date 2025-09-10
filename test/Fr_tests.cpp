#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

extern "C" {
#include "../src/fields.h"
#include "../src/exports.h"
}

#define FR_LIMBS NLIMBS(256)

// Helper: Set a vec256 to a small integer value
void set_fr(vec256 v, uint64_t x) {
    for (int i = 0; i < FR_LIMBS; ++i)
        v[i] = 0;
    memcpy(v, &x, sizeof(x));
}

// Helper: Compare two vec256 for equality
bool fr_equal(const vec256 a, const vec256 b) {
    for (int i = 0; i < FR_LIMBS; ++i)
        if (a[i] != b[i]) return false;
    return true;
}

// Helper: Print vec256 for debugging
void print_fr(const vec256 v) {
    for (int i = FR_LIMBS-1; i >= 0; --i)
        printf("%016llx%s", (unsigned long long)v[i], i ? " " : "\n");
}

// Helper: Convert standard integer to Montgomery form
void to_mont_fr(vec256 r, const vec256 a) {
    // Use BLS12_381_rRR for conversion to Montgomery form
    blst_fr_to(r, a);
}

// Helper: Copy one vec256 to another
void fr_copy(vec256 dst, const vec256 src) {
    for (int i = 0; i < FR_LIMBS; ++i)
        dst[i] = src[i];
}

TEST(FrOps, Add) {
    vec256 a, b, result, expected;
    
    // Basic addition: 5 + 7 = 12
    set_fr(a, 5);
    set_fr(b, 7);
    blst_fr_add(result, a, b);
    set_fr(expected, 12);
    EXPECT_TRUE(fr_equal(result, expected)) << "5 + 7 should equal 12";
    
    // Addition with zero
    set_fr(a, 42);
    set_fr(b, 0);
    blst_fr_add(result, a, b);
    EXPECT_TRUE(fr_equal(result, a)) << "42 + 0 should equal 42";
    
    // Test commutativity: a + b = b + a
    set_fr(a, 123);
    set_fr(b, 456);
    blst_fr_add(result, a, b);
    vec256 result2;
    blst_fr_add(result2, b, a);
    EXPECT_TRUE(fr_equal(result, result2)) << "Addition should be commutative";
}

TEST(FrOps, AddOverflow) {
    vec256 a, b, result, expected;
    
    // Set a = r - 1 (order minus 1)
    fr_copy(a, BLS12_381_r);
    vec256 one;
    set_fr(one, 1);
    blst_fr_sub(a, a, one);
    
    // Add 2 to get 1 (modulo r)
    set_fr(b, 2);
    blst_fr_add(result, a, b);
    set_fr(expected, 1);
    
    EXPECT_TRUE(fr_equal(result, expected)) << "Overflow in addition should wrap modulo r";
}

TEST(FrOps, Sub) {
    vec256 a, b, result, expected;
    
    // Basic subtraction: 10 - 3 = 7
    set_fr(a, 10);
    set_fr(b, 3);
    blst_fr_sub(result, a, b);
    set_fr(expected, 7);
    EXPECT_TRUE(fr_equal(result, expected)) << "10 - 3 should equal 7";
    
    // Subtraction with zero
    set_fr(a, 42);
    set_fr(b, 0);
    blst_fr_sub(result, a, b);
    EXPECT_TRUE(fr_equal(result, a)) << "42 - 0 should equal 42";
}

TEST(FrOps, SubUnderflow) {
    vec256 a, b, result;
    
    // Underflow: 3 - 10 = r - 7
    set_fr(a, 3);
    set_fr(b, 10);
    blst_fr_sub(result, a, b);
    
    // Expected: r - 7
    vec256 expected;
    fr_copy(expected, BLS12_381_r);
    vec256 seven;
    set_fr(seven, 7);
    blst_fr_sub(expected, expected, seven);
    
    EXPECT_TRUE(fr_equal(result, expected)) << "Underflow should wrap modulo r";
}

TEST(FrOps, MulBy3) {
    vec256 a, result, expected;
    
    // 2 * 3 = 6
    set_fr(a, 2);
    blst_fr_mul_by_3(result, a);
    set_fr(expected, 6);
    EXPECT_TRUE(fr_equal(result, expected)) << "2 * 3 should equal 6";
    
    // 0 * 3 = 0
    set_fr(a, 0);
    blst_fr_mul_by_3(result, a);
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result, expected)) << "0 * 3 should equal 0";
    
    // Test with larger value
    set_fr(a, 0x123456789ABCDEF0ULL);
    blst_fr_mul_by_3(result, a);
    
    // Verify by multiplying by 3 using addition
    vec256 temp;
    blst_fr_add(temp, a, a);    // temp = 2*a
    blst_fr_add(expected, temp, a);  // expected = 3*a
    
    EXPECT_TRUE(fr_equal(result, expected)) << "MulBy3 should match manual calculation";
}

TEST(FrOps, LShift) {
    vec256 a, result, expected;
    
    // 1 << 2 = 4
    set_fr(a, 1);
    blst_fr_lshift(result, a, 2);
    set_fr(expected, 4);
    EXPECT_TRUE(fr_equal(result, expected)) << "1 << 2 should equal 4";
    
    // 5 << 1 = 10
    set_fr(a, 5);
    blst_fr_lshift(result, a, 1);
    set_fr(expected, 10);
    EXPECT_TRUE(fr_equal(result, expected)) << "5 << 1 should equal 10";
    
    // 0 << any = 0
    set_fr(a, 0);
    blst_fr_lshift(result, a, 10);
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result, expected)) << "0 << 10 should equal 0";
}

TEST(FrOps, RShift) {
    vec256 a, result, expected;
    
    // 8 >> 2 = 2
    set_fr(a, 8);
    blst_fr_rshift(result, a, 2);
    set_fr(expected, 2);
    EXPECT_TRUE(fr_equal(result, expected)) << "8 >> 2 should equal 2";
    
    // 10 >> 1 = 5
    set_fr(a, 10);
    blst_fr_rshift(result, a, 1);
    set_fr(expected, 5);
    EXPECT_TRUE(fr_equal(result, expected)) << "10 >> 1 should equal 5";
    
    // 0 >> any = 0
    set_fr(a, 0);
    blst_fr_rshift(result, a, 10);
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result, expected)) << "0 >> 10 should equal 0";
}

TEST(FrOps, Mul) {
    vec256 a, b, result, result_std, expected;
    
    // Convert to Montgomery form, multiply, and convert back
    set_fr(a, 3);
    set_fr(b, 4);
    
    to_mont_fr(a, a);
    to_mont_fr(b, b);
    
    blst_fr_mul(result, a, b);
    blst_fr_from(result_std, result);
    
    set_fr(expected, 12);
    EXPECT_TRUE(fr_equal(result_std, expected)) << "3 * 4 should equal 12";
    
    // Test multiplication by zero
    set_fr(a, 42);
    set_fr(b, 0);
    to_mont_fr(a, a);
    to_mont_fr(b, b);
    
    blst_fr_mul(result, a, b);
    blst_fr_from(result_std, result);
    
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result_std, expected)) << "42 * 0 should equal 0";
    
    // Test commutativity: a * b = b * a
    set_fr(a, 123);
    set_fr(b, 456);
    to_mont_fr(a, a);
    to_mont_fr(b, b);
    
    blst_fr_mul(result, a, b);
    vec256 result2;
    blst_fr_mul(result2, b, a);
    
    EXPECT_TRUE(fr_equal(result, result2)) << "Multiplication should be commutative";
}

TEST(FrOps, Sqr) {
    vec256 a, result, result_std, expected;
    
    // 5² = 25
    set_fr(a, 5);
    to_mont_fr(a, a);
    
    blst_fr_sqr(result, a);
    blst_fr_from(result_std, result);
    
    set_fr(expected, 25);
    EXPECT_TRUE(fr_equal(result_std, expected)) << "5² should equal 25";
    
    // 0² = 0
    set_fr(a, 0);
    to_mont_fr(a, a);
    
    blst_fr_sqr(result, a);
    blst_fr_from(result_std, result);
    
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result_std, expected)) << "0² should equal 0";
    
    // Compare sqr with mul by self
    set_fr(a, 123);
    to_mont_fr(a, a);
    
    blst_fr_sqr(result, a);
    vec256 result2;
    blst_fr_mul(result2, a, a);
    
    EXPECT_TRUE(fr_equal(result, result2)) << "a² should equal a*a";
}

TEST(FrOps, CNeg) {
    vec256 a, result, expected;
    
    // Negate 5 when flag is true
    set_fr(a, 5);
    blst_fr_cneg(result, a, 1);
    
    // Expected: r - 5
    fr_copy(expected, BLS12_381_r);
    vec256 five;
    set_fr(five, 5);
    blst_fr_sub(expected, expected, five);
    
    EXPECT_TRUE(fr_equal(result, expected)) << "Negation of 5 should equal r-5";
    
    // Don't negate when flag is false
    set_fr(a, 5);
    blst_fr_cneg(result, a, 0);
    EXPECT_TRUE(fr_equal(result, a)) << "cneg with flag=0 should not change value";
    
    // Negate 0 (should remain 0)
    set_fr(a, 0);
    blst_fr_cneg(result, a, 1);
    set_fr(expected, 0);
    EXPECT_TRUE(fr_equal(result, expected)) << "Negation of 0 should equal 0";
}

TEST(FrOps, MontgomeryConversion) {
    vec256 a, mont_a, back_a;
    
    // Convert 1 to Montgomery form and back
    set_fr(a, 1);
    blst_fr_to(mont_a, a);
    blst_fr_from(back_a, mont_a);
    
    EXPECT_TRUE(fr_equal(back_a, a)) << "Montgomery conversion of 1 failed";
    
    // Test with various values
    const uint64_t TEST_VALUES[] = {0, 2, 42, 12345, 0xFFFFFFFFULL, 0xDEADBEEFULL};
    const int NUM_VALUES = sizeof(TEST_VALUES) / sizeof(TEST_VALUES[0]);
    
    for (int i = 0; i < NUM_VALUES; i++) {
        set_fr(a, TEST_VALUES[i]);
        
        // To Montgomery form and back
        blst_fr_to(mont_a, a);
        blst_fr_from(back_a, mont_a);
        
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
    vec256 a, b, c, result1, result2, mont_a, mont_b, mont_c;
    
    // Test values
    set_fr(a, 17);
    set_fr(b, 23);
    set_fr(c, 42);
    
    // Convert to Montgomery form
    blst_fr_to(mont_a, a);
    blst_fr_to(mont_b, b);
    blst_fr_to(mont_c, c);
    
    // Associativity: (a + b) + c = a + (b + c)
    vec256 temp1, temp2, std_result1, std_result2;
    
    blst_fr_add(temp1, a, b);
    blst_fr_add(result1, temp1, c);
    
    blst_fr_add(temp2, b, c);
    blst_fr_add(result2, a, temp2);
    
    EXPECT_TRUE(fr_equal(result1, result2)) << "Addition should be associative";
    
    // Distributivity: a * (b + c) = a*b + a*c
    blst_fr_add(temp1, mont_b, mont_c);
    blst_fr_mul(result1, mont_a, temp1);
    
    blst_fr_mul(temp1, mont_a, mont_b);
    blst_fr_mul(temp2, mont_a, mont_c);
    blst_fr_add(result2, temp1, temp2);
    
    EXPECT_TRUE(fr_equal(result1, result2)) << "Multiplication should distribute over addition";
}

TEST(FrOps, EdgeCases) {
    vec256 zero, one, a, result, expected;
    
    // Test with zero and one
    set_fr(zero, 0);
    set_fr(one, 1);
    
    // 0 + 0 = 0
    blst_fr_add(result, zero, zero);
    EXPECT_TRUE(fr_equal(result, zero)) << "0 + 0 should equal 0";
    
    // 0 * any = 0
    set_fr(a, 42);
    to_mont_fr(a, a);
    to_mont_fr(zero, zero);
    blst_fr_mul(result, a, zero);
    blst_fr_from(result, result);
    EXPECT_TRUE(fr_equal(result, zero)) << "42 * 0 should equal 0";
    
    // 1 * any = any
    set_fr(a, 42);
    to_mont_fr(a, a);
    to_mont_fr(one, one);
    blst_fr_mul(result, a, one);
    vec256 std_a;
    blst_fr_from(std_a, a);
    blst_fr_from(result, result);
    EXPECT_TRUE(fr_equal(result, std_a)) << "42 * 1 should equal 42";
}

// Test butterfly operations (used in FFT)
TEST(FrOps, Butterfly) {
    vec256 x0, x1, twiddle, expected_x0, expected_x1;
    vec256 original_x0, original_x1;
    
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
    blst_fr_ct_bfly(x0, x1, twiddle);
    
    // Calculate expected values manually
    vec256 temp;
    blst_fr_mul(temp, original_x1, twiddle);  // temp = x1 * twiddle
    blst_fr_add(expected_x0, original_x0, temp);  // x0' = x0 + temp
    blst_fr_sub(expected_x1, original_x0, temp);  // x1' = x0 - temp
    
    EXPECT_TRUE(fr_equal(x0, expected_x0)) << "CT butterfly x0 calculation incorrect";
    EXPECT_TRUE(fr_equal(x1, expected_x1)) << "CT butterfly x1 calculation incorrect";
    
    // Reset for GS butterfly
    fr_copy(x0, original_x0);
    fr_copy(x1, original_x1);
    
    // GS butterfly: x0' = x0 + x1, x1' = (x0 - x1) * twiddle
    blst_fr_gs_bfly(x0, x1, twiddle);
    
    // Calculate expected values manually
    blst_fr_add(expected_x0, original_x0, original_x1);  // x0' = x0 + x1
    blst_fr_sub(temp, original_x0, original_x1);         // temp = x0 - x1
    blst_fr_mul(expected_x1, temp, twiddle);             // x1' = temp * twiddle
    
    EXPECT_TRUE(fr_equal(x0, expected_x0)) << "GS butterfly x0 calculation incorrect";
    EXPECT_TRUE(fr_equal(x1, expected_x1)) << "GS butterfly x1 calculation incorrect";
}