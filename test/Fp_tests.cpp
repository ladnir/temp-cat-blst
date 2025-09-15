#include <gtest/gtest.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>

extern "C" {
#include "../src/fields.h"
}

#define LIMBS NLIMBS(384)


void set_vec384(vec384 v, uint64_t x) {
	for (int i = 0; i < LIMBS; ++i)
		v[i] = 0;

	memcpy(v, &x, sizeof(x));
}

void set_vec384x(vec384x v, uint64_t x0, uint64_t x1) {
	set_vec384(v[0], x0);
	set_vec384(v[1], x1);
}

// Helper: compare two vec384 for equality
bool vec384_equal(const vec384 a, const vec384 b) {
	for (int i = 0; i < LIMBS; ++i)
		if (a[i] != b[i]) return false;
	return true;
}

// Helper: print vec384 (for debugging)
void print_vec384(const vec384 v) {
	for (int i = LIMBS - 1; i >= 0; --i)
		printf("%016llx%s", (unsigned long long)v[i], i ? " " : "\n");
}

//// Helper: Convert standard integer to Montgomery form 
//void to_mont_fp(vec384 r, const vec384 a) {
//	// Use the library's precomputed R² mod p constant for proper conversion to Montgomery form
//	mul_mont_384(r, a, *get_BLST_RR(), *get_BLST_P(), BLST_p0);
//}


TEST(FpOps, Add) {
	vec384 a, b, r, expected;
	set_vec384(a, 5);
	set_vec384(b, 7);
	add_fp(r, a, b);

	// expected = (5 + 7) % p = 12
	set_vec384(expected, 12);
	EXPECT_TRUE(vec384_equal(r, expected));
}

TEST(FpOps, AddOverflow) {
	vec384 a, b, r, expected;
	// a = p - 1, b = 2, expect 1
	memcpy(a, *get_BLST_P(), sizeof(vec384));
	for (int i = 0; i < LIMBS; ++i) a[i] -= (i == 0 ? 1 : 0);
	set_vec384(b, 2);
	add_fp(r, a, b);
	set_vec384(expected, 1);
	EXPECT_TRUE(vec384_equal(r, expected));
}

TEST(FpOps, Sub) {
	vec384 a, b, r, expected;
	set_vec384(a, 10);
	set_vec384(b, 3);
	sub_fp(r, a, b);
	set_vec384(expected, 7);
	EXPECT_TRUE(vec384_equal(r, expected));
}

TEST(FpOps, SubUnderflow) {
	vec384 a, b, r, expected;
	set_vec384(a, 3);
	set_vec384(b, 10);
	sub_fp(r, a, b);
	// expected = (3 - 10) mod p = p - 7
	memcpy(expected, *get_BLST_P(), sizeof(vec384));
	expected[0] -= 7;
	EXPECT_TRUE(vec384_equal(r, expected));
}

TEST(FpOps, MulBy3) {
	vec384 a, r, expected;
	set_vec384(a, 2);
	mul_by_3_fp(r, a);
	set_vec384(expected, 6);
	EXPECT_TRUE(vec384_equal(r, expected));
}

TEST(FpOps, MulBy8) {
	vec384 a, r, expected;
	set_vec384(a, 2);
	mul_by_8_fp(r, a);
	set_vec384(expected, 16);
	EXPECT_TRUE(vec384_equal(r, expected));
}

TEST(FpOps, LShift) {
	vec384 a, r, expected;
	set_vec384(a, 1);
	lshift_fp(r, a, 2);
	set_vec384(expected, 4);
	EXPECT_TRUE(vec384_equal(r, expected));
}

TEST(FpOps, RShift) {
	vec384 a, r, expected;
	set_vec384(a, 8);
	rshift_fp(r, a, 2);
	set_vec384(expected, 2);
	EXPECT_TRUE(vec384_equal(r, expected));
}

TEST(FpOps, DivBy2) {
	vec384 a, r, expected;
	set_vec384(a, 8);
	div_by_2_fp(r, a);
	set_vec384(expected, 4);
	EXPECT_TRUE(vec384_equal(r, expected));
}

TEST(FpOps, Mul) {
	vec384 a, b, r, r_std, expected;

	// Set raw values
	set_vec384(a, 3);
	set_vec384(b, 4);

	// Convert to Montgomery form
	to_fp(a, a);
	to_fp(b, b);

	// Multiply in Montgomery form
	mul_fp(r, a, b);

	// Convert result back to standard form
	from_fp(r_std, r);

	// Expect 12
	set_vec384(expected, 12);

	if (!vec384_equal(r_std, expected)) {
		printf("Actual:   "); print_vec384(r_std);
		printf("Expected: "); print_vec384(expected);
	}

	EXPECT_TRUE(vec384_equal(r_std, expected));
}

TEST(FpOps, Sqr) {
	vec384 a, r, r_std, expected;

	// Set raw value
	set_vec384(a, 5);

	// Convert to Montgomery form
	to_fp(a, a);

	// Square in Montgomery form
	sqr_fp(r, a);

	// Convert result back to standard form
	from_fp(r_std, r);

	set_vec384(expected, 25);
	EXPECT_TRUE(vec384_equal(r_std, expected));
}

TEST(FpOps, FromFp) {
	vec384 a, r, expected;

	// Set raw value and convert to Montgomery form
	set_vec384(a, 9);
	to_fp(a, a);

	// Test converting from Montgomery form
	from_fp(r, a);

	set_vec384(expected, 9);
	EXPECT_TRUE(vec384_equal(r, expected));
}


// Test converting to/from Montgomery representation
TEST(FpOps, MontgomeryConversion) {
	vec384 a, mont_a, back_a;

	// Test with value 1
	set_vec384(a, 1);

	// Convert to Montgomery form using the helper function
	to_fp(mont_a, a);

	// Convert back from Montgomery form
	from_fp(back_a, mont_a);

	// The result should be the original value
	if (!vec384_equal(back_a, a)) {
		printf("Original: "); print_vec384(a);
		printf("Mont form: "); print_vec384(mont_a);
		printf("Back: "); print_vec384(back_a);
	}
	EXPECT_TRUE(vec384_equal(back_a, a)) << "Montgomery conversion failed for value 1";

	// Now test multiplication in Montgomery form
	vec384 b, mont_b, result, std_result;
	set_vec384(a, 3);
	set_vec384(b, 4);

	// Convert to Montgomery form using the helper function
	to_fp(mont_a, a);
	to_fp(mont_b, b);

	// Multiply in Montgomery form
	mul_fp(result, mont_a, mont_b);

	// Convert back to standard form
	from_fp(std_result, result);

	// Check result is 12
	vec384 expected;
	set_vec384(expected, 12);

	if (!vec384_equal(std_result, expected)) {
		printf("Expected: "); print_vec384(expected);
		printf("Got:      "); print_vec384(std_result);
	}

	EXPECT_TRUE(vec384_equal(std_result, expected)) << "Montgomery multiplication failed";
}

// Test Montgomery arithmetic properties
TEST(FpOps, MontgomeryArithmetic) {
	vec384 a, b, c, mont_a, mont_b, mont_c, result1, result2, std_result1, std_result2;

	// Set test values
	set_vec384(a, 17);
	set_vec384(b, 23);

	// Convert to Montgomery form
	to_fp(mont_a, a);
	to_fp(mont_b, b);

	// Test associativity: a * (b * c) = (a * b) * c
	set_vec384(c, 19);
	to_fp(mont_c, c);

	// Calculate a * (b * c) in Montgomery form
	mul_fp(result1, mont_b, mont_c); // result1 = b * c
	mul_fp(result1, mont_a, result1); // result1 = a * (b * c)

	// Calculate (a * b) * c in Montgomery form
	mul_fp(result2, mont_a, mont_b); // result2 = a * b
	mul_fp(result2, result2, mont_c); // result2 = (a * b) * c

	// Convert results back to standard form for comparison
	from_fp(std_result1, result1);
	from_fp(std_result2, result2);

	// Results should be equal
	EXPECT_TRUE(vec384_equal(std_result1, std_result2))
		<< "Montgomery multiplication failed associativity test";

	// Test identity: a * 1 = a
	set_vec384(b, 1);
	to_fp(mont_b, b);

	mul_fp(result1, mont_a, mont_b);
	from_fp(std_result1, result1);
	from_fp(std_result2, mont_a);

	EXPECT_TRUE(vec384_equal(std_result1, std_result2))
		<< "Montgomery identity property failed";

	// Test that squaring works: a² = a * a
	sqr_fp(result1, mont_a);
	mul_fp(result2, mont_a, mont_a);

	EXPECT_TRUE(vec384_equal(result1, result2))
		<< "Montgomery squaring failed";
}

// Add these comprehensive field operation tests

// Test with edge cases
TEST(FpOps, EdgeCases) {
	vec384 zero, one, r, expected;

	// Test operations with zero
	set_vec384(zero, 0);
	set_vec384(one, 1);

	// Zero + Zero = Zero
	add_fp(r, zero, zero);
	EXPECT_TRUE(vec384_equal(r, zero)) << "0 + 0 should equal 0";

	// Zero + One = One
	add_fp(r, zero, one);
	EXPECT_TRUE(vec384_equal(r, one)) << "0 + 1 should equal 1";

	// One - One = Zero
	sub_fp(r, one, one);
	EXPECT_TRUE(vec384_equal(r, zero)) << "1 - 1 should equal 0";

	// Zero - Zero = Zero
	sub_fp(r, zero, zero);
	EXPECT_TRUE(vec384_equal(r, zero)) << "0 - 0 should equal 0";

	// Test with values near modulus
	vec384 near_p, temp;
	memcpy(near_p, *get_BLST_P(), sizeof(vec384));
	for (int i = 0; i < LIMBS; ++i) near_p[i] -= (i == 0 ? 2 : 0);

	// near_p + 1 should be p-1
	add_fp(r, near_p, one);
	memcpy(expected, *get_BLST_P(), sizeof(vec384));
	for (int i = 0; i < LIMBS; ++i) expected[i] -= (i == 0 ? 1 : 0);
	EXPECT_TRUE(vec384_equal(r, expected)) << "(p-2) + 1 should equal (p-1)";

	// Test modular subtraction with underflow
	// (p-1) - p should equal p-1
	memcpy(temp, *get_BLST_P(), sizeof(vec384));
	for (int i = 0; i < LIMBS; ++i) temp[i] -= (i == 0 ? 1 : 0);
	sub_fp(r, temp, *get_BLST_P());
	memcpy(expected, *get_BLST_P(), sizeof(vec384));
	for (int i = 0; i < LIMBS; ++i) expected[i] -= (i == 0 ? 1 : 0);
	EXPECT_TRUE(vec384_equal(r, expected)) << "(p-1) - p should equal (p-1)";
}

// Test Montgomery operations with edge cases and special values
TEST(FpOps, MontgomeryEdgeCases) {
	vec384 zero, one, mont_zero, mont_one, r, expected;

	set_vec384(zero, 0);
	set_vec384(one, 1);

	// Convert 0 and 1 to Montgomery form
	to_fp(mont_zero, zero);
	to_fp(mont_one, one);

	// Test multiplication with 0
	mul_fp(r, mont_zero, mont_one);
	from_fp(expected, r);
	set_vec384(zero, 0);
	EXPECT_TRUE(vec384_equal(expected, zero)) << "0 * 1 in Montgomery form should equal 0";

	// Test multiplication with large values
	vec384 large, mont_large, result;

	// Create a large value (2^128 - 1)
	memset(large, 0, sizeof(vec384));
	large[0] = 0xFFFFFFFFFFFFFFFFULL;
	large[1] = 0xFFFFFFFFFFFFFFFFULL;

	to_fp(mont_large, large);

	// Multiply large * 1
	mul_fp(result, mont_large, mont_one);
	from_fp(r, result);

	// Should get the same large value back (mod p)
	if (!vec384_equal(r, large)) {
		printf("Expected: "); print_vec384(large);
		printf("Got:      "); print_vec384(r);
	}
	EXPECT_TRUE(vec384_equal(r, large)) << "large * 1 should equal large";

	// Test squaring with special value
	vec384 two, mont_two, squared, squared_result;
	set_vec384(two, 2);
	to_fp(mont_two, two);

	// Square 2 (should be 4)
	sqr_fp(result, mont_two);
	from_fp(squared, result);

	set_vec384(squared_result, 4);
	EXPECT_TRUE(vec384_equal(squared, squared_result)) << "2² should equal 4";
}

// Test Montgomery conversion with various values
TEST(FpOps, MontgomeryRoundTrip) {
	const int TEST_VALUES[] = { 0, 1, 2, 3, 7, 10, 100, 12345, 0x7FFFFFFF };
	const int NUM_VALUES = sizeof(TEST_VALUES) / sizeof(TEST_VALUES[0]);

	for (int i = 0; i < NUM_VALUES; i++) {
		vec384 original, mont_form, result;

		set_vec384(original, TEST_VALUES[i]);

		// Convert to Montgomery form and back
		to_fp(mont_form, original);
		from_fp(result, mont_form);

		if (!vec384_equal(result, original)) {
			printf("Failed for value %d\n", TEST_VALUES[i]);
			printf("Original: "); print_vec384(original);
			printf("Mont form: "); print_vec384(mont_form);
			printf("Result: "); print_vec384(result);
		}

		EXPECT_TRUE(vec384_equal(result, original))
			<< "Round-trip Montgomery conversion failed for " << TEST_VALUES[i];
	}
}

// Test distributive property in Montgomery form
TEST(FpOps, MontgomeryDistributive) {
	vec384 a, b, c, mont_a, mont_b, mont_c;
	vec384 left_side, right_side1, right_side2, right_side;
	vec384 std_left, std_right;

	// Set test values
	set_vec384(a, 11);
	set_vec384(b, 13);
	set_vec384(c, 17);

	// Convert to Montgomery form
	to_fp(mont_a, a);
	to_fp(mont_b, b);
	to_fp(mont_c, c);

	// Left side: a * (b + c)
	vec384 sum;
	add_fp(sum, mont_b, mont_c);  // b + c
	mul_fp(left_side, mont_a, sum);  // a * (b + c)

	// Right side: (a * b) + (a * c)
	mul_fp(right_side1, mont_a, mont_b);  // a * b
	mul_fp(right_side2, mont_a, mont_c);  // a * c
	add_fp(right_side, right_side1, right_side2);  // (a * b) + (a * c)

	// Convert back to standard form for comparison
	from_fp(std_left, left_side);
	from_fp(std_right, right_side);

	if (!vec384_equal(std_left, std_right)) {
		printf("Distributive property failed:\n");
		printf("Left side (a * (b + c)): "); print_vec384(std_left);
		printf("Right side ((a * b) + (a * c)): "); print_vec384(std_right);
	}

	EXPECT_TRUE(vec384_equal(std_left, std_right))
		<< "Distributive property failed in Montgomery arithmetic";
}

// Test multiplication with the modulus
TEST(FpOps, MultiplyWithModulus) {
	vec384 a, p, mont_a, mont_p, result, expected;

	// Set a value
	set_vec384(a, 42);

	// Copy the modulus
	memcpy(p, *get_BLST_P(), sizeof(vec384));

	// Convert to Montgomery form
	to_fp(mont_a, a);
	to_fp(mont_p, p);

	// Multiply a * p
	mul_fp(result, mont_a, mont_p);

	// Result should be 0 (since a * p ≡ 0 mod p)
	from_fp(expected, result);

	set_vec384(a, 0);
	EXPECT_TRUE(vec384_equal(expected, a))
		<< "a * p mod p should equal 0";
}

// Helper: copy one vec384 to another
void vec384_copy(vec384 dst, const vec384 src) {
	for (int i = 0; i < LIMBS; ++i)
		dst[i] = src[i];
}

// Test the security properties of constant-time operations
TEST(FpOps, ConstantTimeOperations) {
	// This test doesn't actually verify constant-time execution
	// but ensures the operations work correctly for security-sensitive values

	vec384 a, b, mont_a, mont_b, result, expected;

	// Values that might be used in cryptographic operations
	set_vec384(a, 0xDEADBEEF);  // Some arbitrary "sensitive" value
	set_vec384(b, 0x12345678);  // Another arbitrary value

	// Perform Montgomery conversion and multiplication
	to_fp(mont_a, a);
	to_fp(mont_b, b);
	mul_fp(result, mont_a, mont_b);
	from_fp(expected, result);

	// We're testing that the operation works, but in a real scenario
	// we would want to ensure it's executed in constant time
	// to prevent timing attacks

	// Verify multiplication result by comparing with a direct computation
	// (not practical for large numbers but serves as a test concept)
	uint64_t simple_result = (0xDEADBEEF * 0x12345678) % *get_BLST_P()[0];
	vec384 simple;
	set_vec384(simple, simple_result);

	// This comparison is only meaningful for small values that don't overflow
	// It's just a conceptual test for constant-time operations
	if (simple_result < *get_BLST_P()[0]) {
		EXPECT_EQ(expected[0] & 0xFFFFFFFF, simple_result & 0xFFFFFFFF)
			<< "Lower bits should match for small values";
	}
}

//extern "C" {
//	void sqr_n_mul_fp(vec384 out, const vec384 a, size_t count,
//		const vec384 b);
//}
//
//// Test basic operations to verify they work correctly
//TEST(SqrtDebug, BasicOperations) {
//
//	// Set a = 1 (in Montgomery form)
//	//memset(a, 0, sizeof(vec384));
//
//	vec384 b, c, expected;
//	vec384 a{ ONE_MONT_P };
//
//	// Test squaring: 1² = 1
//	sqr_fp(b, a);
//	//printVec384(a, "a (1 in Mont)");
//	//printVec384(b, "b (1² in Mont)");
//	EXPECT_TRUE(vec_is_equal(a, b, sizeof(vec384)));
//
//	// Test multiplication: 1 × 1 = 1
//	mul_fp(c, a, a);
//	//printVec384(c, "c (1×1 in Mont)");
//	EXPECT_TRUE(vec_is_equal(a, c, sizeof(vec384)));
//}
//
//// Directly test sqr_n_mul_fp with simple inputs
//TEST(SqrtDebug, SquareNMul) {
//
//	vec384 c, expected, result;
//	vec384 a{ ONE_MONT_P };
//	vec384 b{ ONE_MONT_P };
//
//
//	printf("--- Testing sqr_n_mul_fp with count=1 ---\n");
//	// Compute 1² × 1 = 1
//	sqr_n_mul_fp(result, a, 1, b);
//	//printVec384(a, "Input a (1)");
//   // printVec384(b, "Input b (1)");
//	//printVec384(result, "Result (should be 1)");
//
//	// Test intermediate square result
//	vec384 squared;
//	sqr_fp(squared, a);
//	//printVec384(squared, "Squared (should be 1)");
//
//	// Test final multiplication
//	vec384 multiplied;
//	mul_fp(multiplied, squared, b);
//	// printVec384(multiplied, "Multiplied (should be 1)");
//
//	EXPECT_TRUE(vec_is_equal(a, result, sizeof(vec384)));
//
//	// Now try with a larger count
//	printf("--- Testing sqr_n_mul_fp with count=2 ---\n");
//	// Compute 1⁴ × 1 = 1
//	sqr_n_mul_fp(result, a, 2, b);
//	//printVec384(result, "Result (should be 1)");
//	EXPECT_TRUE(vec_is_equal(a, result, sizeof(vec384)));
//
//	// Now try with a different input (a = 2)
//	//memset(a, 0, sizeof(vec384));
//	//a[0] = ONE_MONT_P[0];
//	lshift_fp(a, a, 1);  // a = 2 in Montgomery form
//
//	printf("--- Testing sqr_n_mul_fp with a=2, count=1, b=1 ---\n");
//	// Compute 2² × 1 = 4
//	sqr_n_mul_fp(result, a, 1, b);
//	//printVec384(a, "Input a (2)");
//	//printVec384(result, "Result (should be 4)");
//
//	// Verify by computing 2² = 4 separately
//	sqr_fp(expected, a);
//	//printVec384(expected, "Expected (4)");
//	EXPECT_TRUE(vec_is_equal(expected, result, sizeof(vec384)));
//}
//
//// Test the direct call to the montgomery multiplication with minimal arguments
//TEST(SqrtDebug, DirectMontgomeryCall) {
//	vec384 result;
//
//	vec384 a{ ONE_MONT_P };
//	vec384 b{ ONE_MONT_P };
//
//	printf("--- Testing direct call to sqr_n_mul_mont_383 ---\n");
//	// This tests if the assembly/low-level implementation has issues
//	sqr_n_mul_mont_383(result, a, 1, *get_BLST_P(), BLST_p0, b);
//
//	//printVec384(result, "Direct result");
//	EXPECT_TRUE(vec_is_equal(a, result, sizeof(vec384)));
//}