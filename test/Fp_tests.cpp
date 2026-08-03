#include <gtest/gtest.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "InsecureRng.h"
using namespace blst;

#define LIMBS NLIMBS(384)
using namespace blst;

void set_Fp(Fp& v, uint64_t x) {
	for (int i = 0; i < LIMBS; ++i)
		v[i] = 0;

	memcpy(v.data, &x, sizeof(x));
}

void set_Fp2(Fp2& v, uint64_t x0, uint64_t x1) {
	set_Fp(v[0], x0);
	set_Fp(v[1], x1);
}

void rand_Fp(Fp& v, InsecureRng& rng)
{
	for (u64 i = 0; i < v.num_limb(); ++i)
	{
		v[i] = rng.next_u64();
	}

	while (Fp::is_redc(v) == 0)
	{
		Fp::rshift(v, v, 1);
	}

	BLST_ASSERT( cc_check_mod_n<NLIMBS(384)>((byte*)v.data, get_BLST_P().data));
}

// Helper: compare two Fp for equality
bool Fp_equal(const Fp& a, const Fp& b) {
	return Fp::equal(a, b);
}

// Helper: print Fp (for debugging)
void print_Fp(const Fp& v) {
	for (int i = LIMBS - 1; i >= 0; --i)
		printf("%016llx%s", (unsigned long long)v[i], i ? " " : "\n");
}

//// Helper: Convert standard integer to Montgomery form 
//void to_mont_fp(Fp& r, const Fp& a) {
//	// Use the library's precomputed R² mod p constant for proper conversion to Montgomery form
//	mul_mont_384(r, a, get_BLST_RR(), get_BLST_P(), BLST_p0);
//}


TEST(FpOps, Add) {
	Fp a, b, r, expected;
	set_Fp(a, 5);
	set_Fp(b, 7);
	Fp::add(r, a, b);

	// expected = (5 + 7) % p = 12
	set_Fp(expected, 12);
	EXPECT_TRUE(Fp_equal(r, expected));
}

TEST(FpOps, AddOverflow) {
	Fp a, b, r, expected;
	// a = p - 1, b = 2, expect 1
	memcpy(a.data, get_BLST_P().data, sizeof(Fp));
	for (int i = 0; i < LIMBS; ++i) a[i] -= (i == 0 ? 1 : 0);
	set_Fp(b, 2);
	Fp::add(r, a, b);
	set_Fp(expected, 1);
	EXPECT_TRUE(Fp_equal(r, expected));
}

TEST(FpOps, AddRandom)
{

	// test that (a + b) + c == a + (b + c)
	Fp a, b, c, r1, r2;
	InsecureRng rng;

	for (u64 i = 0; i < 100; ++i)
	{
		rand_Fp(a, rng);
		rand_Fp(b, rng);
		rand_Fp(c, rng);

		Fp::add(r1, a, b);
		Fp::add(r1, r1, c);

		Fp::add(r2, c, b);
		Fp::add(r2, r2, a);


		auto eq = Fp::equal(r1, r2);
		if (!eq)
		{
			printf("a ");
			print_Fp(a);
			printf("\nb ");
			print_Fp(b);
			printf("\nc ");
			print_Fp(c);
			printf("\nr1 ");
			print_Fp(r1);
			printf("\nr2 ");
			print_Fp(r2);



		}
		EXPECT_TRUE(eq);
	}

}

TEST(FpOps, Sub) {
	Fp a, b, r, expected;
	set_Fp(a, 10);
	set_Fp(b, 3);
	Fp::sub(r, a, b);
	set_Fp(expected, 7);
	EXPECT_TRUE(Fp_equal(r, expected));
}

TEST(FpOps, SubUnderflow) {
	Fp a, b, r, expected;
	set_Fp(a, 3);
	set_Fp(b, 10);
	Fp::sub(r, a, b);
	// expected = (3 - 10) mod p = p - 7
	memcpy(expected.data, get_BLST_P().data, sizeof(Fp));
	expected[0] -= 7;
	EXPECT_TRUE(Fp_equal(r, expected));
}


TEST(FpOps, SubRandom)
{
	Fp a, b, c, r1, r2;
	InsecureRng rng;

	for (u64 i = 0; i < 100; ++i)
	{
		rand_Fp(a, rng);
		rand_Fp(b, rng);
		Fp::add(c, a, b);

		Fp::sub(r1, c, a);
		EXPECT_TRUE(Fp::equal(r1, b));

		Fp::sub(r2, c, b);
		EXPECT_TRUE(Fp::equal(r2, a));
	}

}


TEST(FpOps, MulBy3) {
	Fp a, r, expected;
	InsecureRng rng;
	for (u64 i = 0; i < 100; ++i)
	{

		rand_Fp(a, rng);
		Fp::mul_by_3(r, a);

		Fp::add(expected, a, a);
		Fp::add(expected, expected, a);

		EXPECT_TRUE(Fp_equal(r, expected));
	}
}

TEST(FpOps, MulBy8) {
	Fp a, r, expected;
	InsecureRng rng;
	for (u64 i = 0; i < 100; ++i)
	{

		rand_Fp(a, rng);
		Fp::mul_by_8(r, a);

		Fp::add(expected, a, a);
		Fp::add(expected, expected, expected);
		Fp::add(expected, expected, expected);

		EXPECT_TRUE(Fp_equal(r, expected));
	}
}

TEST(FpOps, LShift) {
	Fp a, r, expected;
	InsecureRng rng;
	rand_Fp(a, rng);
	Fp::lshift(r, a, 0);
	EXPECT_TRUE(Fp_equal(r, a));
	Fp::lshift(a, a, 0);
	EXPECT_TRUE(Fp_equal(r, a));
	expected = a;
	for (u64 i = 1; i < 100; ++i)
	{
		Fp::add(expected, expected, expected);
		Fp::lshift(r, a, i);
		EXPECT_TRUE(Fp_equal(r, expected));

	}
}


TEST(FpOps, RShift) {
	Fp a, r, expected, la;
	set_Fp(a, 0x123456789abcdef0ULL);
	Fp::rshift(r, a, 0);
	EXPECT_TRUE(Fp_equal(r, a));
	Fp::rshift(a, a, 0);
	EXPECT_TRUE(Fp_equal(r, a));
	set_Fp(a, 8);
	Fp::rshift(r, a, 2);
	set_Fp(expected, 2);
	EXPECT_TRUE(Fp_equal(r, expected));

	InsecureRng rng;
	for (u64 i = 1; i < 100; ++i)
	{
		memset(&a, 0, sizeof(a));
		a[0] = rng.next_u64();
		Fp::lshift(la, a, i);
		Fp::rshift(r, la, i);
		EXPECT_TRUE(Fp_equal(r, a));
	}
}

TEST(FpOps, DivBy2) {
	Fp a, r, expected;
	set_Fp(a, 8);
	Fp::div_by_2(r, a);
	set_Fp(expected, 4);
	EXPECT_TRUE(Fp_equal(r, expected));
}

TEST(FpOps, Mul) {
	Fp a, b, r, r_std, expected;

	// Set raw values
	set_Fp(a, 3);
	set_Fp(b, 4);

	// Convert to Montgomery form
	Fp::to_mont(a, a);
	Fp::to_mont(b, b);

	// Multiply in Montgomery form
	Fp::mul(r, a, b);

	// Convert result back to standard form
	Fp::from_mont(r_std, r);

	// Expect 12
	set_Fp(expected, 12);

	if (!Fp_equal(r_std, expected)) {
		printf("Actual:   "); print_Fp(r_std);
		printf("Expected: "); print_Fp(expected);
	}

	EXPECT_TRUE(Fp_equal(r_std, expected));
}


TEST(FpOps, MulRand) {
	Fp a, b, c, r0, r1;
	InsecureRng rng;
	for (u64 i = 0; i < 100; ++i)
	{
		// Set raw values
		rand_Fp(a, rng);
		rand_Fp(b, rng);
		rand_Fp(c, rng);

		Fp::mul(r0, a, b);
		Fp::mul(r0, r0, c);

		Fp::mul(r1, a, c);
		Fp::mul(r1, r1, b);

		EXPECT_TRUE(Fp_equal(r0, r1));
	}
}

TEST(FpOps, Sqr) {
	for(u64 i =0; i < 100; ++i)
	{
		Fp a, a2, r;
		InsecureRng rng;
		rand_Fp(a, rng);
		Fp::mul(a2, a, a);
		Fp::sqr(r, a);
		EXPECT_TRUE(Fp_equal(r, a2));
	}
}


TEST(FpOps, SqrNMul) {
	for(u64 i =0; i < 20; ++i)
	{
		for (u64 j = 1; j < 10; ++j)
		{

			Fp a, b, a2, r;
			InsecureRng rng;
			rand_Fp(a, rng);
			rand_Fp(b, rng);

			a2 = a;
			for(u64 jj = 0; jj <j; ++jj)
				Fp::sqr(a2, a2);

			Fp::mul(a2, a2, b);
			Fp::sqr_n_mul(r, a, j, b);
			EXPECT_TRUE(Fp_equal(r, a2));
		}
	}
}



TEST(FpOps, Reciprocal) {
	Fp out, inp, r;
	InsecureRng rng;
	for (u64 i = 0; i < 100; ++i)
	{
		rand_Fp(out, rng);
		Fp::reciprocal(inp, out);
		Fp::mul(r, inp, out);
		EXPECT_TRUE(Fp::equal(r, Fp::one()));
	}
}

TEST(FpOps, RecipSqrt) {
	Fp out; Fp out2; Fp inp; Fp r; Fp t1; Fp t2; 
	InsecureRng rng;
	for (u64 i = 0; i < 100; ++i) {
		do { rand_Fp(out, rng); } while (Fp::equal(out, Fp::zero()));
		Fp::sqr(out2, out);
		Fp::reciprocal(inp, out2);                  // inp = 1/out^2
		bool ok = Fp::recip_sqrt(r, inp);
		EXPECT_TRUE(ok);
		Fp::sqr(t1, r);                             // t1 = r^2
		Fp::mul(t2, t1, inp);                       // t2 = r^2 * inp
		EXPECT_TRUE(Fp::equal(t2, Fp::one()));                 // definition check
		Fp t_r2;
		Fp t_o2;
		Fp::sqr(t_r2, r);          // r^2
		Fp::sqr(t_o2, out);        // out^2
		EXPECT_TRUE(Fp::equal(t_r2, t_o2));
	}
}


TEST(FpOps, exp) {
	Fp a;
	Fp r;
	Fp expected;
	Fp base;
	InsecureRng rng;
	std::vector<u8> e(8);

	auto bit = [&](size_t i) -> bool {
		return (e[i >> 3] >> (i & 7)) & 1;
		};

	for (u64 t = 0; t < 100; ++t) {
		rand_Fp(a, rng);
		for (auto& b : e) b = static_cast<u8>(rng.next_u64());

		size_t bits = e.size() * 8;
		while (bits > 0 && !bit(bits - 1)) --bits;       // trim leading zeros
		if (bits == 0) { e[0] = 1; bits = 1; }           // avoid zero exponent

		Fp::exp(r, a, e.data(), bits);

		expected = Fp::one();
		base = a;
		for (size_t i = 0; i < bits; ++i) {
			if (bit(i)) {
				Fp tmp;
				Fp::mul(tmp, expected, base);
				expected = tmp;
			}
			Fp tmp2;
			Fp::sqr(tmp2, base);
			base = tmp2;
		}

		EXPECT_TRUE(Fp::equal(r, expected));
	}
}

TEST(FpOps, FromFp) {
	Fp a, r, expected;

	// Set raw value and convert to Montgomery form
	set_Fp(a, 9);
	Fp::to_mont(a, a);

	// Test converting from Montgomery form
	Fp::from_mont(r, a);

	set_Fp(expected, 9);
	EXPECT_TRUE(Fp_equal(r, expected));
}


// Test converting to/from Montgomery representation
TEST(FpOps, MontgomeryConversion) {
	Fp a, mont_a, back_a;

	// Test with value 1
	set_Fp(a, 1);

	// Convert to Montgomery form using the helper function
	Fp::to_mont(mont_a, a);

	// Convert back from Montgomery form
	Fp::from_mont(back_a, mont_a);

	// The result should be the original value
	if (!Fp_equal(back_a, a)) {
		printf("Original: "); print_Fp(a);
		printf("Mont form: "); print_Fp(mont_a);
		printf("Back: "); print_Fp(back_a);
	}
	EXPECT_TRUE(Fp_equal(back_a, a)) << "Montgomery conversion failed for value 1";

	// Now test multiplication in Montgomery form
	Fp b, mont_b, result, std_result;
	set_Fp(a, 3);
	set_Fp(b, 4);

	// Convert to Montgomery form using the helper function
	Fp::to_mont(mont_a, a);
	Fp::to_mont(mont_b, b);

	// Multiply in Montgomery form
	Fp::mul(result, mont_a, mont_b);

	// Convert back to standard form
	Fp::from_mont(std_result, result);

	// Check result is 12
	Fp expected;
	set_Fp(expected, 12);

	if (!Fp_equal(std_result, expected)) {
		printf("Expected: "); print_Fp(expected);
		printf("Got:      "); print_Fp(std_result);
	}

	EXPECT_TRUE(Fp_equal(std_result, expected)) << "Montgomery multiplication failed";
}

// Test Montgomery arithmetic properties
TEST(FpOps, MontgomeryArithmetic) {
	Fp a, b, c, mont_a, mont_b, mont_c, result1, result2, std_result1, std_result2;

	// Set test values
	set_Fp(a, 17);
	set_Fp(b, 23);

	// Convert to Montgomery form
	Fp::to_mont(mont_a, a);
	Fp::to_mont(mont_b, b);

	// Test associativity: a * (b * c) = (a * b) * c
	set_Fp(c, 19);
	Fp::to_mont(mont_c, c);

	// Calculate a * (b * c) in Montgomery form
	Fp::mul(result1, mont_b, mont_c); // result1 = b * c
	Fp::mul(result1, mont_a, result1); // result1 = a * (b * c)

	// Calculate (a * b) * c in Montgomery form
	Fp::mul(result2, mont_a, mont_b); // result2 = a * b
	Fp::mul(result2, result2, mont_c); // result2 = (a * b) * c

	// Convert results back to standard form for comparison
	Fp::from_mont(std_result1, result1);
	Fp::from_mont(std_result2, result2);

	// Results should be equal
	EXPECT_TRUE(Fp_equal(std_result1, std_result2))
		<< "Montgomery multiplication failed associativity test";

	// Test identity: a * 1 = a
	set_Fp(b, 1);
	Fp::to_mont(mont_b, b);

	Fp::mul(result1, mont_a, mont_b);
	Fp::from_mont(std_result1, result1);
	Fp::from_mont(std_result2, mont_a);

	EXPECT_TRUE(Fp_equal(std_result1, std_result2))
		<< "Montgomery identity property failed";

	// Test that squaring works: a² = a * a
	Fp::sqr(result1, mont_a);
	Fp::mul(result2, mont_a, mont_a);

	EXPECT_TRUE(Fp_equal(result1, result2))
		<< "Montgomery squaring failed";
}

// Add these comprehensive field operation tests

// Test with edge cases
TEST(FpOps, EdgeCases) {
	Fp zero, one, r, expected;

	// Test operations with zero
	set_Fp(zero, 0);
	set_Fp(one, 1);

	// Zero + Zero = Zero
	Fp::add(r, zero, zero);
	EXPECT_TRUE(Fp_equal(r, zero)) << "0 + 0 should equal 0";

	// Zero + One = One
	Fp::add(r, zero, one);
	EXPECT_TRUE(Fp_equal(r, one)) << "0 + 1 should equal 1";

	// One - One = Zero
	Fp::sub(r, one, one);
	EXPECT_TRUE(Fp_equal(r, zero)) << "1 - 1 should equal 0";

	// Zero - Zero = Zero
	Fp::sub(r, zero, zero);
	EXPECT_TRUE(Fp_equal(r, zero)) << "0 - 0 should equal 0";

	// Test with values near modulus
	Fp near_p;
	memcpy(near_p.data, get_BLST_P().data, sizeof(Fp));
	for (int i = 0; i < LIMBS; ++i) near_p[i] -= (i == 0 ? 2 : 0);

	// near_p + 1 should be p-1
	Fp::add(r, near_p, one);
	memcpy(expected.data, get_BLST_P().data, sizeof(Fp));
	for (int i = 0; i < LIMBS; ++i) expected[i] -= (i == 0 ? 1 : 0);
	EXPECT_TRUE(Fp_equal(r, expected)) << "(p-2) + 1 should equal (p-1)";

}

// Test Montgomery operations with edge cases and special values
TEST(FpOps, MontgomeryEdgeCases) {
	Fp zero, one, mont_zero, mont_one, r, expected;

	set_Fp(zero, 0);
	set_Fp(one, 1);

	// Convert 0 and 1 to Montgomery form
	Fp::to_mont(mont_zero, zero);
	Fp::to_mont(mont_one, one);

	// Test multiplication with 0
	Fp::mul(r, mont_zero, mont_one);
	Fp::from_mont(expected, r);
	set_Fp(zero, 0);
	EXPECT_TRUE(Fp_equal(expected, zero)) << "0 * 1 in Montgomery form should equal 0";

	// Test multiplication with large values
	Fp large, mont_large, result;

	// Create a large value (2^128 - 1)
	memset(large.data, 0, sizeof(Fp));
	large[0] = 0xFFFFFFFFFFFFFFFFULL;
	large[1] = 0xFFFFFFFFFFFFFFFFULL;

	Fp::to_mont(mont_large, large);

	// Multiply large * 1
	Fp::mul(result, mont_large, mont_one);
	Fp::from_mont(r, result);

	// Should get the same large value back (mod p)
	if (!Fp_equal(r, large)) {
		printf("Expected: "); print_Fp(large);
		printf("Got:      "); print_Fp(r);
	}
	EXPECT_TRUE(Fp_equal(r, large)) << "large * 1 should equal large";

	// Test squaring with special value
	Fp two, mont_two, squared, squared_result;
	set_Fp(two, 2);
	Fp::to_mont(mont_two, two);

	// Square 2 (should be 4)
	Fp::sqr(result, mont_two);
	Fp::from_mont(squared, result);

	set_Fp(squared_result, 4);
	EXPECT_TRUE(Fp_equal(squared, squared_result)) << "2² should equal 4";
}

// Test Montgomery conversion with various values
TEST(FpOps, MontgomeryRoundTrip) {
	const int TEST_VALUES[] = { 0, 1, 2, 3, 7, 10, 100, 12345, 0x7FFFFFFF };
	const int NUM_VALUES = sizeof(TEST_VALUES) / sizeof(TEST_VALUES[0]);

	for (int i = 0; i < NUM_VALUES; i++) {
		Fp original, mont_form, result;

		set_Fp(original, TEST_VALUES[i]);

		// Convert to Montgomery form and back
		Fp::to_mont(mont_form, original);
		Fp::from_mont(result, mont_form);

		if (!Fp_equal(result, original)) {
			printf("Failed for value %d\n", TEST_VALUES[i]);
			printf("Original: "); print_Fp(original);
			printf("Mont form: "); print_Fp(mont_form);
			printf("Result: "); print_Fp(result);
		}

		EXPECT_TRUE(Fp_equal(result, original))
			<< "Round-trip Montgomery conversion failed for " << TEST_VALUES[i];
	}
}

// Test distributive property in Montgomery form
TEST(FpOps, MontgomeryDistributive) {
	Fp a, b, c, mont_a, mont_b, mont_c;
	Fp left_side, right_side1, right_side2, right_side;
	Fp std_left, std_right;

	// Set test values
	set_Fp(a, 11);
	set_Fp(b, 13);
	set_Fp(c, 17);

	// Convert to Montgomery form
	Fp::to_mont(mont_a, a);
	Fp::to_mont(mont_b, b);
	Fp::to_mont(mont_c, c);

	// Left side: a * (b + c)
	Fp sum;
	Fp::add(sum, mont_b, mont_c);  // b + c
	Fp::mul(left_side, mont_a, sum);  // a * (b + c)

	// Right side: (a * b) + (a * c)
	Fp::mul(right_side1, mont_a, mont_b);  // a * b
	Fp::mul(right_side2, mont_a, mont_c);  // a * c
	Fp::add(right_side, right_side1, right_side2);  // (a * b) + (a * c)

	// Convert back to standard form for comparison
	Fp::from_mont(std_left, left_side);
	Fp::from_mont(std_right, right_side);

	if (!Fp_equal(std_left, std_right)) {
		printf("Distributive property failed:\n");
		printf("Left side (a * (b + c)): "); print_Fp(std_left);
		printf("Right side ((a * b) + (a * c)): "); print_Fp(std_right);
	}

	EXPECT_TRUE(Fp_equal(std_left, std_right))
		<< "Distributive property failed in Montgomery arithmetic";
}

// Test multiplication with the modulus
TEST(FpOps, MultiplyWithModulus) {
	Fp a, p, mont_a, mont_p, result, expected;

	// Set a value
	set_Fp(a, 42);

	// Copy the modulus
	memcpy(p.data, get_BLST_P().data, sizeof(Fp));

	// Convert to Montgomery form
	Fp::to_mont(mont_a, a);
	Fp::to_mont(mont_p, p);

	// Multiply a * p
	Fp::mul(result, mont_a, mont_p);

	// Result should be 0 (since a * p ≡ 0 mod p)
	Fp::from_mont(expected, result);

	set_Fp(a, 0);
	EXPECT_TRUE(Fp_equal(expected, a))
		<< "a * p mod p should equal 0";
}

// Helper: copy one Fp to another
void Fp_copy(Fp& dst, const Fp& src) {
	for (int i = 0; i < LIMBS; ++i)
		dst[i] = src[i];
}

// Test the security properties of constant-time operations
TEST(FpOps, ConstantTimeOperations) {
	// This test doesn't actually verify constant-time execution
	// but ensures the operations work correctly for security-sensitive values

	Fp a, b, mont_a, mont_b, result, expected;

	// Values that might be used in cryptographic operations
	set_Fp(a, 0xDEADBEEF);  // Some arbitrary "sensitive" value
	set_Fp(b, 0x12345678);  // Another arbitrary value

	// Perform Montgomery conversion and multiplication
	Fp::to_mont(mont_a, a);
	Fp::to_mont(mont_b, b);
	Fp::mul(result, mont_a, mont_b);
	Fp::from_mont(expected, result);

	// We're testing that the operation works, but in a real scenario
	// we would want to ensure it's executed in constant time
	// to prevent timing attacks

	// Verify multiplication result by comparing with a direct computation
	// (not practical for large numbers but serves as a test concept)
	uint64_t simple_result = (0xDEADBEEF * 0x12345678) % get_BLST_P()[0];
	Fp simple;
	set_Fp(simple, simple_result);

	// This comparison is only meaningful for small values that don't overflow
	// It's just a conceptual test for constant-time operations
	if (simple_result < get_BLST_P()[0]) {
		EXPECT_EQ(expected[0] & 0xFFFFFFFF, simple_result & 0xFFFFFFFF)
			<< "Lower bits should match for small values";
	}
}

//extern "C" {
//	void sqr_n_Fp::mul(Fp out, const Fp a, size_t count,
//		const Fp b);
//}
//
//// Test basic operations to verify they work correctly
//TEST(SqrtDebug, BasicOperations) {
//
//	// Set a = 1 (in Montgomery form)
//	//memset(a, 0, sizeof(Fp));
//
//	Fp b, c, expected;
//	Fp a{ ONE_MONT_P };
//
//	// Test squaring: 1² = 1
//	Fp::sqr(b, a);
//	//printVec384(a, "a (1 in Mont)");
//	//printVec384(b, "b (1² in Mont)");
//	EXPECT_TRUE(vec_is_equal(a, b, sizeof(Fp)));
//
//	// Test multiplication: 1 × 1 = 1
//	Fp::mul(c, a, a);
//	//printVec384(c, "c (1×1 in Mont)");
//	EXPECT_TRUE(vec_is_equal(a, c, sizeof(Fp)));
//}
//
//// Directly test sqr_n_Fp::mul with simple inputs
//TEST(SqrtDebug, SquareNMul) {
//
//	Fp c, expected, result;
//	Fp a{ ONE_MONT_P };
//	Fp b{ ONE_MONT_P };
//
//
//	printf("--- Testing sqr_n_Fp::mul with count=1 ---\n");
//	// Compute 1² × 1 = 1
//	sqr_n_Fp::mul(result, a, 1, b);
//	//printVec384(a, "Input a (1)");
//   // printVec384(b, "Input b (1)");
//	//printVec384(result, "Result (should be 1)");
//
//	// Test intermediate square result
//	Fp squared;
//	Fp::sqr(squared, a);
//	//printVec384(squared, "Squared (should be 1)");
//
//	// Test final multiplication
//	Fp multiplied;
//	Fp::mul(multiplied, squared, b);
//	// printVec384(multiplied, "Multiplied (should be 1)");
//
//	EXPECT_TRUE(vec_is_equal(a, result, sizeof(Fp)));
//
//	// Now try with a larger count
//	printf("--- Testing sqr_n_Fp::mul with count=2 ---\n");
//	// Compute 1⁴ × 1 = 1
//	sqr_n_Fp::mul(result, a, 2, b);
//	//printVec384(result, "Result (should be 1)");
//	EXPECT_TRUE(vec_is_equal(a, result, sizeof(Fp)));
//
//	// Now try with a different input (a = 2)
//	//memset(a, 0, sizeof(Fp));
//	//a[0] = ONE_MONT_P[0];
//	lshift_fp(a, a, 1);  // a = 2 in Montgomery form
//
//	printf("--- Testing sqr_n_Fp::mul with a=2, count=1, b=1 ---\n");
//	// Compute 2² × 1 = 4
//	sqr_n_Fp::mul(result, a, 1, b);
//	//printVec384(a, "Input a (2)");
//	//printVec384(result, "Result (should be 4)");
//
//	// Verify by computing 2² = 4 separately
//	Fp::sqr(expected, a);
//	//printVec384(expected, "Expected (4)");
//	EXPECT_TRUE(vec_is_equal(expected, result, sizeof(Fp)));
//}
//
//// Test the direct call to the montgomery multiplication with minimal arguments
//TEST(SqrtDebug, DirectMontgomeryCall) {
//	Fp result;
//
//	Fp a{ ONE_MONT_P };
//	Fp b{ ONE_MONT_P };
//
//	printf("--- Testing direct call to sqr_n_mul_mont_383 ---\n");
//	// This tests if the assembly/low-level implementation has issues
//	sqr_n_mul_mont_383(result, a, 1, get_BLST_P(), BLST_p0, b);
//
//	//printVec384(result, "Direct result");
//	EXPECT_TRUE(vec_is_equal(a, result, sizeof(Fp)));
//}
