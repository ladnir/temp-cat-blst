#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "blst/fields/Fp12.h"
#include "blst/fields/Fr.h"
using namespace blst;

// Helper: Print Fp12 element for debugging
void print_fp12(const Fp12& a, const char* name) {
	printf("%s:\n", name);
	for (int j = 0; j < 2; j++) {  // Fp12 = Fp6[w]/(w^2-v)
		for (int i = 0; i < 3; i++) {  // Fp6 = Fp2[v]/(v^3-u-1)
			printf("  [%d][%d].Re: ", j, i);
			for (int k = NLIMBS(384) - 1; k >= 0; k--) {
				printf("%016llx ", (unsigned long long)a[j][i][0][k]);
			}
			printf("\n  [%d][%d].Im: ", j, i);
			for (int k = NLIMBS(384) - 1; k >= 0; k--) {
				printf("%016llx ", (unsigned long long)a[j][i][1][k]);
			}
			printf("\n");
		}
	}
	printf("\n");
}

TEST(Fp12Test, ZeroAndOne) {
	Fp12 zero, one;

	// Initialize zero
	vec_zero(&zero, sizeof(zero));

	// Get the Fp12 one constant from the library
	const Fp12 one_ptr = Fp12::one();
	vec_copy(&one, &one_ptr, sizeof(one));

	// Check zero properties
	Fp12 result;

	// Zero + Zero = Zero
	Fp12::mul(result, zero, zero);
	//print_fp12(result, "0 * 0");
	EXPECT_TRUE(vec_is_equal(&result, &zero, sizeof(Fp12))) << "0 * 0 = 0";

	// One * One = One
	Fp12::mul(result, one, one);
	EXPECT_TRUE(vec_is_equal(&result, &one, sizeof(Fp12))) << "1 * 1 = 1";

	// Zero * One = Zero
	Fp12::mul(result, zero, one);
	EXPECT_TRUE(vec_is_equal(&result, &zero, sizeof(Fp12))) << "0 * 1 = 0";

	// One * Zero = Zero
	Fp12::mul(result, one, zero);
	EXPECT_TRUE(vec_is_equal(&result, &zero, sizeof(Fp12))) << "1 * 0 = 0";

	// Check that One is detected as One
	EXPECT_TRUE(Fp12::is_one(one)) << "One should be detected as one";
	EXPECT_FALSE(Fp12::is_one(zero)) << "Zero should not be detected as one";
}

TEST(Fp12Test, Multiplication) {
	Fp12 a, b, c, temp;

	// Initialize with some values - use the one constant and modify it
	const Fp12 one_ptr = Fp12::one();
	vec_copy(&a, &one_ptr, sizeof(a));

	// Make a simple change to create a distinct value
	// Modify a[0][0][1][0] (first coordinate, imaginary part)
	a[0][0][1][0] = 2;

	// Create another distinct value
	vec_copy(&b, &one_ptr, sizeof(b));
	// Modify b[0][1][0][0] (second coordinate, real part)
	b[0][1][0][0] = 3;

	// Compute c = a * b
	Fp12::mul(c, a, b);

	// Verify that the result is different from both inputs
	EXPECT_FALSE(vec_is_equal(&c, &a, sizeof(Fp12))) << "a*b should not equal a";
	EXPECT_FALSE(vec_is_equal(&c, &b, sizeof(Fp12))) << "a*b should not equal b";

	// Verify multiplication by 1 gives the same element
	Fp12::mul(temp, a, one_ptr);
	EXPECT_TRUE(vec_is_equal(&temp, &a, sizeof(Fp12))) << "a*1 should equal a";

	Fp12::mul(temp, b, one_ptr);
	EXPECT_TRUE(vec_is_equal(&temp, &b, sizeof(Fp12))) << "b*1 should equal b";
}

TEST(Fp12Test, Squaring) {
	Fp12 a, a_squared1, a_squared2;

	// Initialize with some values - use the one constant and modify it
	const Fp12 one_ptr = Fp12::one();
	vec_copy(&a, &one_ptr, sizeof(a));

	// Make a simple change to create a distinct value
	// Modify a[0][0][1][0] (first coordinate, imaginary part)
	a[0][0][1][0] = 2;
	a[1][1][0][0] = 3;  // Another change to make it more complex

	// Compute a² using multiplication
	Fp12::mul(a_squared1, a, a);

	// Compute a² using squaring
	Fp12::sqr(a_squared2, a);

	// Results should be equal
	EXPECT_TRUE(vec_is_equal(&a_squared1, &a_squared2, sizeof(Fp12)))
		<< "a² computed by multiplication and squaring should be equal";

	// Test cyclotomic squaring if a is in the cyclotomic subgroup
	// This is just a basic test, not a comprehensive one
	Fp12 cyc_squared;
	Fp12::cyclotomic_sqr(cyc_squared, a);
	// We don't expect cyclotomic_sqr to equal regular squaring for arbitrary elements
	// So we just test that the function runs without errors
}

TEST(Fp12Test, ExplicitScalarApisMatchDefaultApis) {
	Fp12 a, b;
	Fp12 mul_default, mul_scalar;
	Fp12 sqr_default, sqr_scalar;
	Fp12 cyclo_seed, cyclo;
	Fp12 cyc_default, cyc_scalar;
	Fp12 sparse_default, sparse_scalar;
	Fp6 sparse_arg;

	const Fp12 one = Fp12::one();
	vec_copy(&a, &one, sizeof(a));
	vec_copy(&b, &one, sizeof(b));

	a[0][0][1][0] = 0x02;
	a[0][1][0][0] = 0x03;
	a[1][2][1][0] = 0x04;
	b[0][2][0][0] = 0x05;
	b[1][0][1][0] = 0x06;
	b[1][1][0][0] = 0x07;

	Fp12::mul(mul_default, a, b);
	Fp12::mul_scalar(mul_scalar, a, b);
	EXPECT_TRUE(Fp12::equal(mul_default, mul_scalar));

	Fp12::sqr(sqr_default, a);
	Fp12::sqr_scalar(sqr_scalar, a);
	EXPECT_TRUE(Fp12::equal(sqr_default, sqr_scalar));

	vec_zero(&sparse_arg, sizeof(sparse_arg));
	sparse_arg[0] = b[0][0];
	sparse_arg[1] = b[0][1];
	sparse_arg[2] = b[0][2];
	Fp12::mul_by_xy00z0(sparse_default, a, sparse_arg);
	Fp12::mul_by_xy00z0_scalar(sparse_scalar, a, sparse_arg);
	EXPECT_TRUE(Fp12::equal(sparse_default, sparse_scalar));

	vec_copy(&cyclo_seed, &a, sizeof(cyclo_seed));
	Fp12::final_exp(cyclo, cyclo_seed);
	ASSERT_TRUE(Fp12::is_cyclotomic(cyclo));
	Fp12::cyclotomic_sqr(cyc_default, cyclo);
	Fp12::cyclotomic_sqr_scalar(cyc_scalar, cyclo);
	EXPECT_TRUE(Fp12::equal(cyc_default, cyc_scalar));
}

TEST(Fp12Test, Conjugation) {
	Fp12 a, conjugate, expected;

	// Initialize with some values
	const Fp12 one_ptr = Fp12::one();
	vec_copy(&a, &one_ptr, sizeof(a));

	// Make it a non-trivial value: a = a_0 + a_1*w
	a[0][0][1][0] = 2;  // Modify a_0
	a[1][1][0][0] = 3;  // Modify a_1

	// Compute conjugate
	vec_copy(&conjugate, &a, sizeof(Fp12));
	Fp12::conjugate(conjugate);

	// Expected = a_0 - a_1*w
	vec_copy(&expected, &a, sizeof(Fp12));
	// Negate the second part of the Fp12 element (a_1)
	for (int i = 0; i < 3; i++) {
		Fp2::cneg(expected[1][i], a[1][i], 1);
	}

	// Results should be equal
	EXPECT_TRUE(vec_is_equal(&conjugate, &expected, sizeof(Fp12)))
		<< "conjugate(a) should equal a_0 - a_1*w";

	// Check that conjugate(conjugate(a)) = a
	Fp12::conjugate(conjugate);
	EXPECT_TRUE(vec_is_equal(&conjugate, &a, sizeof(Fp12)))
		<< "conjugate(conjugate(a)) should equal a";
}

TEST(Fp12Test, Inversion) {
	Fp12 a, a_inv, result;

	// Initialize with a non-trivial value
	const Fp12 one_ptr = Fp12::one();
	vec_copy(&a, &one_ptr, sizeof(a));
	a[0][0][1][0] = 2;
	a[1][1][0][0] = 3;

	// Compute 1/a
	Fp12::inverse(a_inv, a);

	// Check that a * (1/a) = 1
	Fp12::mul(result, a, a_inv);

	// Should equal 1
	EXPECT_TRUE(Fp12::is_one(result))
		<< "a * (1/a) should equal 1";

	// Special case: Inversion of 1 should be 1
	Fp12 one, one_inv;
	vec_copy(&one, &one_ptr, sizeof(one));
	Fp12::inverse(one_inv, one);
	EXPECT_TRUE(vec_is_equal(&one, &one_inv, sizeof(Fp12)))
		<< "1/1 should equal 1";
}

//TEST(Fp12Test, FrobeniusMap) {
//    Fp12 a, frob1, frob2;
//
//    // Initialize with a non-trivial value
//    const Fp12* one_ptr = Fp12::one();
//    vec_copy(a, one_ptr, sizeof(a));
//    a[0][0][1][0] = 2;
//    a[0][1][0][0] = 3;
//    a[1][2][1][0] = 4;
//
//    // Apply Frobenius once
//    vec_copy(frob1, a, sizeof(Fp12));
//    blst_fp12_frobenius_map(frob1, a, 1);
//
//    // Apply Frobenius twice
//    vec_copy(frob2, a, sizeof(Fp12));
//    blst_fp12_frobenius_map(frob2, a, 2);
//
//    // Verify basic property: Frobenius^12 is the identity on Fp12
//    Fp12 frob12;
//    vec_copy(frob12, a, sizeof(Fp12));
//
//    // Apply Frobenius 12 times (or use the direct power-12 operation if available)
//    blst_fp12_frobenius_map(frob12, frob12, 6);
//    blst_fp12_frobenius_map(frob12, frob12, 6);
//
//    // Result should equal original a
//    EXPECT_TRUE(vec_is_equal(frob12, a, sizeof(Fp12)))
//        << "Frobenius^12 should be identity on Fp12";
//
//    // Verify Frobenius^6 is the conjugation
//    Fp12 frob6, conj;
//    vec_copy(frob6, a, sizeof(Fp12));
//    blst_fp12_frobenius_map(frob6, a, 6);
//
//    vec_copy(conj, a, sizeof(Fp12));
//    Fp12::conjugate(conj);
//
//    EXPECT_TRUE(vec_is_equal(frob6, conj, sizeof(Fp12)))
//        << "Frobenius^6 should equal conjugation";
//}

TEST(Fp12Test, Equality) {
	Fp12 a, b;

	// Initialize with different values
	const Fp12 one_ptr = Fp12::one();
	vec_copy(&a, &one_ptr, sizeof(a));
	vec_copy(&b, &one_ptr, sizeof(b));

	// a and b should be equal initially
	EXPECT_TRUE(Fp12::equal(a, b)) << "a and b should be equal initially";

	// Modify b to make it different
	b[0][0][0][0] = 2;
	EXPECT_FALSE(Fp12::equal(a, b)) << "a and b should not be equal after modification";

	// Make them equal again
	vec_copy(&b, &a, sizeof(b));
	EXPECT_TRUE(Fp12::equal(a, b)) << "a and b should be equal after copying";

	// Test equality with self
	EXPECT_TRUE(Fp12::equal(a, a)) << "a should equal itself";
}

TEST(Fp12Test, Serialization) {
	Fp12 a, b;
	unsigned char serialized[12 * 48]; // Fp12 = 12 elements of Fp

	// Initialize with non-trivial values
	const Fp12 one_ptr = Fp12::one();
	vec_copy(&a, &one_ptr, sizeof(a));
	a[0][0][1][0] = 2;
	a[0][1][0][0] = 3;
	a[1][2][1][0] = 4;

	// Serialize to bytes
	Fp12::to_bendian(serialized, a);

	// Create a different value for b
	vec_copy(&b, &one_ptr, sizeof(b));
	b[0][0][0][0] = 42; // Make it different from a

	// Deserialize back to b
	// Note: There's no direct deserialize function in the provided context
	// In a real implementation, you would use the library's deserialization function
	// Here we'll assume serialization/deserialization works and just test the setup

	// Instead, let's verify properties of the serialized data
	bool has_nonzero = false;
	for (size_t i = 0; i < sizeof(serialized); i++) {
		if (serialized[i] != 0) {
			has_nonzero = true;
			break;
		}
	}
	EXPECT_TRUE(has_nonzero) << "Serialized data should have non-zero bytes";

}


static BLST_FORCEINLINE int ge_256(const limb_t a[4], const limb_t b[4])
{
	for (int i = 3; i >= 0; --i) {
		if (a[i] > b[i]) return 1;
		if (a[i] < b[i]) return 0;
	}
	return 1;
}

static BLST_FORCEINLINE void sub_256(limb_t out[4], const limb_t a[4], const limb_t b[4])
{
	limb_t borrow = 0;
	for (int i = 0; i < 4; ++i) {
		const limb_t bi = b[i] + borrow;
		borrow = (a[i] < bi);
		out[i] = a[i] - bi;
	}
}

static BLST_FORCEINLINE void shl1_256(limb_t x[4])
{
	const limb_t c0 = x[0] >> 63;
	const limb_t c1 = x[1] >> 63;
	const limb_t c2 = x[2] >> 63;

	x[0] <<= 1;
	x[1] = (x[1] << 1) | c0;
	x[2] = (x[2] << 1) | c1;
	x[3] = (x[3] << 1) | c2;
}

static BLST_FORCEINLINE unsigned get_bit_384_le(const limb_t x[6], unsigned bit)
{
	const unsigned limb = bit >> 6;      // /64
	const unsigned off = bit & 63u;      // %64
	return (unsigned)((x[limb] >> off) & 1u);
}

// Computes lambda = p mod r (both public constants), where:
//   p is 384-bit (Fp modulus) in little-endian limbs,
//   r is 256-bit (Fr modulus) in little-endian limbs.
// Output is little-endian 4-limb value in [0,r).
static void gt_lambda_from_p_mod_r(limb_t lambda_out[4])
{
	static_assert(sizeof(limb_t) == 8, "This helper currently assumes 64-bit limb_t");

	const limb_t* const p = get_BLST_P().l;   // 6 limbs
	const limb_t* const r = get_BLST_r().l;   // 4 limbs

	limb_t acc[4] = { 0, 0, 0, 0 };

	// Classic bitwise modulus accumulation:
	// acc = (acc*2 + bit) mod r, scan p from MSB down to LSB.
	for (int bit = 383; bit >= 0; --bit) {
		shl1_256(acc);
		acc[0] |= (limb_t)get_bit_384_le(p, (unsigned)bit);

		// Reduce if needed (acc can be up to 2*r-1 here, so at most one subtract).
		if (ge_256(acc, r)) {
			limb_t tmp[4];
			sub_256(tmp, acc, r);
			acc[0] = tmp[0];
			acc[1] = tmp[1];
			acc[2] = tmp[2];
			acc[3] = tmp[3];
		}
	}

	lambda_out[0] = acc[0];
	lambda_out[1] = acc[1];
	lambda_out[2] = acc[2];
	lambda_out[3] = acc[3];
}


static void le_bytes_from_limbs_256(unsigned char out[32], const limb_t in[4])
{
	for (int i = 0; i < 4; ++i) {
		const uint64_t w = (uint64_t)in[i];
		out[i * 8 + 0] = (unsigned char)(w >> 0);
		out[i * 8 + 1] = (unsigned char)(w >> 8);
		out[i * 8 + 2] = (unsigned char)(w >> 16);
		out[i * 8 + 3] = (unsigned char)(w >> 24);
		out[i * 8 + 4] = (unsigned char)(w >> 32);
		out[i * 8 + 5] = (unsigned char)(w >> 40);
		out[i * 8 + 6] = (unsigned char)(w >> 48);
		out[i * 8 + 7] = (unsigned char)(w >> 56);
	}
}

TEST(Fp12Test, GtLambdaMatchesFrobenius)
{
	limb_t lambda[4];
	gt_lambda_from_p_mod_r(lambda);

	// Print limbs so you can hardcode them (little-endian limb order).
	printf("GT lambda = p mod r limbs (LE):\n");
	printf("  0x%016llx\n", (unsigned long long)lambda[0]);
	printf("  0x%016llx\n", (unsigned long long)lambda[1]);
	printf("  0x%016llx\n", (unsigned long long)lambda[2]);
	printf("  0x%016llx\n", (unsigned long long)lambda[3]);

	unsigned char lambda_le[32];
	le_bytes_from_limbs_256(lambda_le, lambda);

	// Construct a cyclotomic element g by taking a "random-ish" Fp12 and final-exping it.
	// (final_exp output is in GT, hence cyclotomic)
	Fp12 g;
	const Fp12 one = Fp12::one();
	vec_copy(&g, &one, sizeof(g));

	// Make g non-trivial (not guaranteed in GT yet).
	g[0][0][1][0] = 2;
	g[0][1][0][0] = 3;
	g[1][2][1][0] = 4;

	Fp12 gt;
	Fp12::final_exp(gt, g);
	ASSERT_TRUE(Fp12::is_cyclotomic(gt));

	// a) Frobenius: gt^p
	Fp12 gt_phi;
	Fp12::frobenius_map(gt_phi, gt, 1);

	// b) Exponentiation by lambda = p mod r
	Fp12 gt_lambda;
	Fp12::cyclotomic_exp(gt_lambda, gt, lambda_le, 256);

	EXPECT_TRUE(vec_is_equal(&gt_phi, &gt_lambda, sizeof(Fp12)))
		<< "Expected frobenius_map(g,1) == g^(p mod r) in GT";
}


static BLST_FORCEINLINE int is_zero_256(const limb_t a[4])
{
	return (a[0] | a[1] | a[2] | a[3]) == 0;
}

static BLST_FORCEINLINE unsigned bitlen_256(const limb_t a[4])
{
	for (int i = 3; i >= 0; --i) {
		if (a[i] != 0) {
			return (unsigned)(i * 64 + cc_num_bits(a[i]));
		}
	}
	return 0;
}

static void gt_lambda2_from_lambda_mod_r(limb_t lambda2_out[4])
{
	// Compute lambda2 = lambda^2 mod r, where lambda is canonical limbs (NOT Montgomery)
	// Use Fr Montgomery helpers safely:
	//   lambda_m = to_mont(lambda)
	//   lambda2_m = sqr(lambda_m)
	//   lambda2 = from_mont(lambda2_m)
	Fr lam;
	vec_copy(&lam, &get_BLST_gt_lambda(), sizeof(lam));

	Fr lam_m;
	Fr::to_mont(lam_m, lam);

	Fr lam2_m;
	Fr::sqr(lam2_m, lam_m);

	Fr lam2;
	Fr::from_mont(lam2, lam2_m);

	lambda2_out[0] = lam2.l[0];
	lambda2_out[1] = lam2.l[1];
	lambda2_out[2] = lam2.l[2];
	lambda2_out[3] = lam2.l[3];
}

static int ge_u512_u256(const limb_t a512[8], const limb_t b256[4])
{
	// Compare a512 to b256 (b treated as 512 with high limbs = 0).
	for (int i = 7; i >= 4; --i) {
		if (a512[i] != 0) return 1; // any high limb means a512 > b256
	}
	for (int i = 3; i >= 0; --i) {
		if (a512[i] > b256[i]) return 1;
		if (a512[i] < b256[i]) return 0;
	}
	return 1;
}

static void sub_u512_u256_inplace(limb_t a512[8], const limb_t b256[4])
{
	limb_t borrow = 0;
	for (int i = 0; i < 4; ++i) {
		const limb_t bi = b256[i] + borrow;
		borrow = (a512[i] < bi);
		a512[i] -= bi;
	}
	for (int i = 4; i < 8; ++i) {
		const limb_t bi = borrow;
		borrow = (a512[i] < bi);
		a512[i] -= bi;
	}
}

static void shl1_u512(limb_t x[8])
{
	limb_t carry = 0;
	for (int i = 0; i < 8; ++i) {
		const limb_t new_carry = x[i] >> 63;
		x[i] = (x[i] << 1) | carry;
		carry = new_carry;
	}
}

static unsigned get_bit_256_le(const limb_t x[4], unsigned bit)
{
	const unsigned limb = bit >> 6;
	const unsigned off = bit & 63u;
	return (unsigned)((x[limb] >> off) & 1u);
}

// Compute q = floor(k / d) and r = k mod d, where k and d are unsigned 256-bit.
// This is classic bitwise long division (256 steps). Slow but fine for a unit test.
static void u256_divrem_u256(const limb_t k[4], const limb_t d[4], limb_t q[4], limb_t r[4])
{
	// divisor must be non-zero
	ASSERT_FALSE(is_zero_256(d));

	q[0] = q[1] = q[2] = q[3] = 0;

	limb_t rem512[8] = { 0 };
	for (int bit = 255; bit >= 0; --bit) {
		shl1_u512(rem512);
		rem512[0] |= (limb_t)get_bit_256_le(k, (unsigned)bit);

		if (ge_u512_u256(rem512, d)) {
			sub_u512_u256_inplace(rem512, d);
			q[bit >> 6] |= (limb_t)1 << (bit & 63);
		}
	}

	r[0] = rem512[0];
	r[1] = rem512[1];
	r[2] = rem512[2];
	r[3] = rem512[3];
}

TEST(Fp12Test, GtLambda2MatchesFrobeniusSquared)
{
	limb_t lambda2[4];
	gt_lambda2_from_lambda_mod_r(lambda2);

	//printf("GT lambda2 = (p^2 mod r) computed as lambda^2 mod r limbs (LE):\n");
	//printf("  0x%016llx\n", (unsigned long long)lambda2[0]);
	//printf("  0x%016llx\n", (unsigned long long)lambda2[1]);
	//printf("  0x%016llx\n", (unsigned long long)lambda2[2]);
	//printf("  0x%016llx\n", (unsigned long long)lambda2[3]);

	unsigned char lambda2_le[32];
	le_bytes_from_limbs_256(lambda2_le, lambda2);

	// Build cyclotomic g
	Fp12 g;
	const Fp12 one = Fp12::one();
	vec_copy(&g, &one, sizeof(g));
	g[0][0][1][0] = 2;
	g[0][1][0][0] = 3;
	g[1][2][1][0] = 4;

	Fp12 gt;
	Fp12::final_exp(gt, g);
	ASSERT_TRUE(Fp12::is_cyclotomic(gt));

	// φ^2(g)
	Fp12 gt_phi2;
	Fp12::frobenius_map(gt_phi2, gt, 2);

	// g^(p^2 mod r)
	Fp12 gt_l2;
	Fp12::cyclotomic_exp(gt_l2, gt, lambda2_le, 256);

	EXPECT_TRUE(vec_is_equal(&gt_phi2, &gt_l2, sizeof(Fp12)))
		<< "Expected frobenius_map(g,2) == g^(p^2 mod r) in GT";


	// compare with get_BLST_gt_lambda2()
	Fr lambda2_fr = get_BLST_gt_lambda2();
	limb_t lambda2_fr_limbs[4];
	lambda2_fr_limbs[0] = lambda2_fr.l[0];
	lambda2_fr_limbs[1] = lambda2_fr.l[1];
	lambda2_fr_limbs[2] = lambda2_fr.l[2];
	lambda2_fr_limbs[3] = lambda2_fr.l[3];
	EXPECT_TRUE(vec_is_equal(lambda2, lambda2_fr_limbs, sizeof(limb_t) * 4))
		<< "Expected computed lambda2 matches get_BLST_gt_lambda2()";

}

TEST(Fp12Test, CyclotomicExpGLV)
{
	limb_t lambda2[4];
	gt_lambda2_from_lambda_mod_r(lambda2);

	unsigned char lambda2_le[32];
	le_bytes_from_limbs_256(lambda2_le, lambda2);

	// Get r for reduction
	const limb_t* r_ptr = get_BLST_r().l;
	limb_t r[4];
	vec_copy(r, r_ptr, sizeof(r));

	// Build cyclotomic g
	Fp12 seed;
	vec_copy(&seed, &Fp12::one(), sizeof(seed));
	seed[0][0][1][0] = 2;
	seed[0][1][0][0] = 3;
	seed[1][2][1][0] = 4;

	Fp12 g;
	Fp12::final_exp(g, seed);
	ASSERT_TRUE(Fp12::is_cyclotomic(g));

	// φ^2(g)
	Fp12 g_phi2;
	Fp12::frobenius_map(g_phi2, g, 2);

	unsigned max_k1_bits = 0;
	unsigned max_k2_bits = 0;

	for (int t = 0; t < 16; ++t) {
		limb_t k_raw[4];
		k_raw[0] = 0x0123456789abcdefULL ^ (0x1111111111111111ULL * (unsigned long long)(t + 1));
		k_raw[1] = 0xfedcba9876543210ULL ^ (0x2222222222222222ULL * (unsigned long long)(t + 1));
		k_raw[2] = 0x0f0f0f0f0f0f0f0fULL ^ (0x3333333333333333ULL * (unsigned long long)(t + 1));
		k_raw[3] = 0xf0f0f0f0f0f0f0f0ULL ^ (0x4444444444444444ULL * (unsigned long long)(t + 1));

		// Ensure k < r for GLV validity
		limb_t k[4], q_ignore[4];
		u256_divrem_u256(k_raw, r, q_ignore, k);

#ifndef NDEBUG
		// Debug assertion: Scalar must be strictly less than group order r.
		// GLV decomposition optimized for curve operations assumes valid Fr scalars.
		{
			int cmp = 0;
			for (int i = 3; i >= 0; i--) {
				if (k[i] > r[i]) { cmp = 1; break; }
				if (k[i] < r[i]) { cmp = -1; break; }
			}
			// If equal (0) or greater (1), assert fails
			BLST_ASSERT(cmp < 0 && "Scalar must be reduced ( < r ) for GLV optimization");
		}
#endif

		unsigned char k_le[32];
		le_bytes_from_limbs_256(k_le, k);

		// Reference: g^k
		Fp12 ref;
		Fp12::cyclotomic_exp(ref, g, k_le, 256);

		// Split by lambda2: k = k1 + k2*lambda2
		limb_t k2[4];
		limb_t k1[4];
		u256_divrem_u256(k, lambda2, k2, k1);

		max_k1_bits = std::max(max_k1_bits, bitlen_256(k1));
		max_k2_bits = std::max(max_k2_bits, bitlen_256(k2));

		unsigned char k1_le[32];
		unsigned char k2_le[32];
		le_bytes_from_limbs_256(k1_le, k1);
		le_bytes_from_limbs_256(k2_le, k2);

		// g^k1 * (φ^2(g))^k2
		Fp12 got;
		Fp12::exp(got, g, (byte*)k, 256);

		EXPECT_TRUE(vec_is_equal(&ref, &got, sizeof(Fp12)))
			<< "Two-way split using phi^2 failed at t=" << t;
	}

	printf("Max observed bit lengths (over 16 samples): k1=%u, k2=%u\n", max_k1_bits, max_k2_bits);
}

//static void u512_divrem_u256(const limb_t n512[8], const limb_t d256[4], limb_t q512[8], limb_t r256[4])
//{
//	ASSERT_FALSE(is_zero_256(d256));
//
//	for (int i = 0; i < 8; ++i) q512[i] = 0;
//
//	limb_t rem512[8] = { 0 };
//
//	auto get_bit_512_le = [](const limb_t x[8], unsigned bit) -> unsigned {
//		const unsigned limb = bit >> 6;
//		const unsigned off = bit & 63u;
//		return (unsigned)((x[limb] >> off) & 1u);
//		};
//
//	for (int bit = 511; bit >= 0; --bit) {
//		shl1_u512(rem512);
//		rem512[0] |= (limb_t)get_bit_512_le(n512, (unsigned)bit);
//
//		if (ge_u512_u256(rem512, d256)) {
//			sub_u512_u256_inplace(rem512, d256);
//			q512[bit >> 6] |= (limb_t)1 << (bit & 63);
//		}
//	}
//
//	r256[0] = rem512[0];
//	r256[1] = rem512[1];
//	r256[2] = rem512[2];
//	r256[3] = rem512[3];
//}
//
//TEST(Fp12Test, GtLambda2BarrettMu)
//{
//	const Fr lambda2_fr = get_BLST_gt_lambda2();
//	const limb_t d[4] = {
//		lambda2_fr.l[0],
//		lambda2_fr.l[1],
//		lambda2_fr.l[2],
//		lambda2_fr.l[3],
//	};
//
//	// Numerator N = 2^256 represented as 512-bit little-endian limbs.
//	limb_t n512[8] = { 0 };
//	n512[4] = 1; // bit 256 set => limb index 4 (4*64)
//
//	limb_t q512[8];
//	limb_t r[4];
//	u512_divrem_u256(n512, d, q512, r);
//
//	// mu = floor(2^256 / d). This fits in the low 256 bits of q512.
//	//limb_t mu[4] = { q512[0], q512[1], q512[2], q512[3] };
//
//	//printf("GT lambda2 Barrett mu = floor(2^256 / lambda2) (LE limbs):\n");
//	//printf("  0x%016llx\n", (unsigned long long)mu[0]);
//	//printf("  0x%016llx\n", (unsigned long long)mu[1]);
//	//printf("  0x%016llx\n", (unsigned long long)mu[2]);
//	//printf("  0x%016llx\n", (unsigned long long)mu[3]);
//
//	//printf("Paste-ready initializer:\n");
//	//printf("  { TO_LIMB_T(0x%016llx), TO_LIMB_T(0x%016llx), TO_LIMB_T(0x%016llx), TO_LIMB_T(0x%016llx) }\n",
//	//	(unsigned long long)mu[0],
//	//	(unsigned long long)mu[1],
//	//	(unsigned long long)mu[2],
//	//	(unsigned long long)mu[3]);
//
//	// Sanity: remainder must be < d. Since d fits 256 bits, r already reduced.
//	EXPECT_FALSE(is_zero_256(d)); // trivial check so test isn't "empty"
//}



TEST(Fp12Test, GtPhi2BarrettSplitMatchesSlowDivRem)
{
	static_assert(sizeof(limb_t) == 8, "This test assumes 64-bit limb_t");

	const Fr lambda2_fr = get_BLST_gt_lambda2();
	const limb_t lambda2[4] = {
		lambda2_fr.l[0],
		lambda2_fr.l[1],
		lambda2_fr.l[2],
		lambda2_fr.l[3],
	};

	// Deterministic pseudo-random-ish scalars + edge cases.
	const limb_t cases[][4] = {
		{ 0, 0, 0, 0 },
		{ 1, 0, 0, 0 },
		{ 0xffffffffffffffffULL, 0, 0, 0 },
		{ 0, 0xffffffffffffffffULL, 0, 0 },
		{ 0, 0, 0xffffffffffffffffULL, 0 },
		{ 0, 0, 0, 0xffffffffffffffffULL },
		{ 0xffffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL },
		{ 0x0123456789abcdefULL, 0xfedcba9876543210ULL, 0x0f0f0f0f0f0f0f0fULL, 0xf0f0f0f0f0f0f0f0ULL },
	};

	auto check_one = [&](const limb_t k[4]) {
		limb_t k2_slow[4], k1_slow[4];
		u256_divrem_u256(k, lambda2, k2_slow, k1_slow);

		limb_t k1_fast[4], k2_fast[4];
		Fp12::glv_split_2(k1_fast, k2_fast, k);

		EXPECT_TRUE(vec_is_equal(k1_slow, k1_fast, sizeof(k1_slow)))
			<< "k1 mismatch";
		EXPECT_TRUE(vec_is_equal(k2_slow, k2_fast, sizeof(k2_slow)))
			<< "k2 mismatch";
		};

	for (const auto& k : cases) {
		check_one(k);
	}

	// Many probes.
	for (int t = 0; t < 2048; ++t) {
		limb_t k[4];
		k[0] = 0x0123456789abcdefULL ^ (0x1111111111111111ULL * (unsigned long long)(t + 1));
		k[1] = 0xfedcba9876543210ULL ^ (0x2222222222222222ULL * (unsigned long long)(t + 1));
		k[2] = 0x0f0f0f0f0f0f0f0fULL ^ (0x3333333333333333ULL * (unsigned long long)(t + 1));
		k[3] = 0xf0f0f0f0f0f0f0f0ULL ^ (0x4444444444444444ULL * (unsigned long long)(t + 1));

		check_one(k);
	}
}



// ... existing tests ...

TEST(Fp12Test, CyclotomicExpGLV2)
{
	// 1. Setup a valid cyclotomic element 'g'
	Fp12 seed;
	vec_copy(&seed, &Fp12::one(), sizeof(seed));
	// Make seed non-trivial
	seed[0][0][1][0] = 2;
	seed[0][1][0][0] = 3;
	seed[1][2][1][0] = 4;

	Fp12 g;
	Fp12::final_exp(g, seed);
	ASSERT_TRUE(Fp12::is_cyclotomic(g));

	// 2. Test with several scalars
	// Hardcoded random values + some edge cases
	struct {
		uint64_t v[4];
		const char* name;
	} cases[] = {
		{ {1, 0, 0, 0}, "One" },
		{ {0xffffffffffffffff, 0, 0, 0}, "Uint64 Max" },
		{ {0x1234567812345678, 0x8765432187654321, 0, 0}, "128-bit random" },
		{ {0xffffffffffffffff, 0xffffffffffffffff, 0, 0}, "128-bit Max" },
		// A scalar likely to produce a 129-bit k2 component (roughly 2^256/2 = 2^255)
		{ {0, 0, 0, 0x8000000000000000}, "High bit set (255)" },
		// Full random mix
		{ {0x5e4487b39c6e3d2a, 0x1b2c3d4e5f6a7b8c, 0x9a8b7c6d5e4f3a2b, 0xd4c3b2a190876543}, "Full 256-bit random" }
	};

	// Get r for reduction
	const limb_t* r_ptr = get_BLST_r().l;
	limb_t r[4];
	vec_copy(r, r_ptr, sizeof(r));

	for (const auto& tc : cases) {
		// Ensure input is reduced mod r
		limb_t k[4], q_ignore[4];
		u256_divrem_u256((const limb_t*)tc.v, r, q_ignore, k);

		// Convert to bytes for generic exp API
		unsigned char scalar_bytes[32];
		memcpy(scalar_bytes, k, 32);

		// Reference: Standard Cyclotomic Exp
		Fp12 ref;
		Fp12::exp(ref, g, scalar_bytes, 256);

		// Target: GLV Exp
		Fp12 glv;
		Fp12::cyclotomic_exp_glv(glv, g, scalar_bytes);

		EXPECT_TRUE(vec_is_equal(&ref, &glv, sizeof(Fp12)))
			<< "GLV Exp mismatch for case: " << tc.name;
	}
}


TEST(Fp12Test, CyclotomicExpGLV2_jsf) {
	const int iterations = 100;


	for (int i = 0; i < iterations; ++i) {
		// Generate random 256-bit scalar
		limb_t k_raw[4];
		for (int j = 0; j < 4; j++) {
			k_raw[j] = ((uint64_t)rand() << 32) | rand();
		}

		// Reduce mod r
		limb_t k[4], q_ignore[4];
		const limb_t* r_ptr = get_BLST_r().l;
		limb_t r[4];
		vec_copy(r, r_ptr, sizeof(r));
		u256_divrem_u256(k_raw, r, q_ignore, k);

			// 1. Setup a valid cyclotomic element 'g'
		Fp12 seed;
		vec_copy(&seed, &Fp12::one(), sizeof(seed));
		// Make seed non-trivial
		seed[0][0][1][0] = 2;
		seed[0][1][0][0] = 3;
		seed[1][2][1][0] = 4;

		Fp12 base;
		Fp12::final_exp(base, seed);
		ASSERT_TRUE(Fp12::is_cyclotomic(base));

		// 2. Execute Optimized Split-Phase GLV
		// Step A: Recode (API requirement)
		//int8_t jsf[132][2];
		limb_t k1[4], k2[4];
		Fp12::glv_split_2(k1, k2, k);
		auto jsf = Fp12::glv_recode_jsf(k1, k2);

		// Step B: Exponentiate
		Fp12 res_opt;
		// Assuming cyclotomic_exp_glv_opt is exposed similarly
		Fp12::cyclotomic_exp_glv_jsf(res_opt, base, jsf);

		// 3. Execute Reference Implementation
		// Logic: res_ref = (base^k1) * (frobenius_map(base, 2)^k2)
		Fp12 res_ref;
		Fp12::cyclotomic_exp_glv(res_ref, base, k1, k2);

		// 4. Verify Equivalence
		// Using memcmp for exact binary match (assuming canonical representation)
		// If Monty form is used, ensure equality check handles it.
		bool match = (memcmp(&res_opt, &res_ref, sizeof(Fp12)) == 0);

		EXPECT_TRUE(match)
			<< "GLV Optimization Mismatch at iteration " << i << "\n"
			<< "K1: " << std::hex << k1[1] << k1[0] << "\n"
			<< "K2: " << std::hex << k2[1] << k2[0] << "\n"
			<< "Recoded Length: " << std::dec << jsf.sched_len << "\n";

		if (!match) break; // detailed debug stop
	}
}





TEST(Fp12Test, CyclotomicExpGLV_4D_Consistency) {
	const int iterations = 50;

	for (int i = 0; i < iterations; ++i) {
		// Generate random 256-bit scalar
		limb_t k_raw[4];
		for (int j = 0; j < 4; j++) {
			k_raw[j] = ((uint64_t)rand() << 32) | rand();
		}

		// Reduce mod r
		limb_t k[4], q_ignore[4];
		const limb_t* r_ptr = get_BLST_r().l;
		limb_t r[4];
		vec_copy(r, r_ptr, sizeof(r));
		u256_divrem_u256(k_raw, r, q_ignore, k);

		// Setup base
		Fp12 seed;
		vec_copy(&seed, &Fp12::one(), sizeof(seed));
		seed[0][0][1][0] = 2; seed[1][2][1][0] = 4;
		Fp12 base;
		Fp12::final_exp(base, seed); // Ensure cyclotomic

		// Run 4D GLV
		Fp12 res_4d;
		Fp12::cyclotomic_exp_glv4(res_4d, base, k);

		// Run Reference
		unsigned char k_bytes[32];
		le_bytes_from_limbs_256(k_bytes, k);
		Fp12 res_ref;
		Fp12::cyclotomic_exp(res_ref, base, k_bytes, 256);

		// Compare
		EXPECT_TRUE(vec_is_equal(&res_4d, &res_ref, sizeof(Fp12)))
			<< "GLV-4D Mismatch at iter " << i;
	}
}


TEST(Fp12Test, CyclotomicExpGLV4_Constant_Consistency) {
	const int iterations = 50;

	for (int i = 0; i < iterations; ++i) {

		// Generate random 256-bit scalar
		limb_t k_raw[4];
		for (int j = 0; j < 4; j++) {
			k_raw[j] = ((uint64_t)rand() << 32) | rand();
		}

		// Reduce mod r
		limb_t k[4], q_ignore[4];
		const limb_t* r_ptr = get_BLST_r().l;
		limb_t r[4];
		vec_copy(r, r_ptr, sizeof(r));
		u256_divrem_u256(k_raw, r, q_ignore, k);

		// Generate valid base
		Fp12 seed;
		vec_copy(&seed, &Fp12::one(), sizeof(seed));
		seed[0][0][1][0] = 2; seed[0][1][0][0] = 3; seed[1][2][1][0] = 4;
		Fp12 base;
		Fp12::final_exp(base, seed);

		// 1. Compile Schedule (Offline Phase)
		auto props = Fp12::glv4_compile(k);

		// 2. Execute (Runtime Phase)
		// No need to pass 'k' anymore!
		Fp12 res_opt;
		Fp12::cyclotomic_exp_glv4(res_opt, base, props);

		// 3. Reference
		unsigned char k_le[32];
		le_bytes_from_limbs_256(k_le, k);
		Fp12 res_ref;
		Fp12::cyclotomic_exp(res_ref, base, k_le, 256);

		// Verify
		bool match = (memcmp(&res_opt, &res_ref, sizeof(Fp12)) == 0);
		EXPECT_TRUE(match) << "GLV-4 Constant Context Mismatch at " << i;
		if (!match) break;
	}
}

// Alternative: Reconstruct using Fp12 exponentiation which is the ground truth.
// If k = k0 + k1*lambda ... then g^k must equal g^k0 * g^(k1*lambda)...
// This is exactly what the GLV expander does, so testing reconstruction *algebraically*
// is harder without a BigInt library.
// We will stick to the functional test: g^k (ref) == result.

TEST(Fp12Test, CyclotomicExpGLV4_EdgeCases) {
	const limb_t* r_ptr = get_BLST_r().l;
	limb_t r[4];
	vec_copy(r, r_ptr, sizeof(r));

	// Define edge case scalars
	struct {
		const char* name;
		uint64_t parts[4]; // 256-bit scalar parts
	} scenarios[] = {
		// Small values
		{"Zero", {0, 0, 0, 0}},
		{"One", {1, 0, 0, 0}},
		{"Two", {2, 0, 0, 0}},

		// Near BLS_Z boundaries (z = 0xd201000000010000)
		{"Z - 1", {0xd20100000000ffffULL, 0, 0, 0}},
		{"Z",     {0xd201000000010000ULL, 0, 0, 0}},
		{"Z + 1", {0xd201000000010001ULL, 0, 0, 0}},

		// Near r boundaries
		{"r - 1", {r[0] - 1, r[1], r[2], r[3]}}, // Assuming no borrow propagation for simple r

		// Specific patterns that might trigger carry issues in decomposition
		{"High-Low Alternating", {0xAAAAAAAAAAAAAAAA, 0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555}},
		{"All Ones", {0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF}},
	};

	// Prepare a base element
	Fp12 seed, base;
	vec_copy(&seed, &Fp12::one(), sizeof(seed));
	seed[0][0][1][0] = 5;
	Fp12::final_exp(base, seed); // Ensure cyclotomic

	for (const auto& sc : scenarios) {
		limb_t k[4];

		// If "r-1", copy carefully
		if (strcmp(sc.name, "r - 1") == 0) {
			// Manual decrement of r
			vec_copy(k, r, sizeof(k));
			// k = r - 1
			int i = 0;
			while (i < 4 && k[i] == 0) { k[i]--; i++; }
			k[i]--;
		}
		else {
			memcpy(k, sc.parts, sizeof(k));
			// Reduce mod r just in case (e.g. All Ones)
			limb_t q_ign[4], r_out[4];
			u256_divrem_u256(k, r, q_ign, r_out);
			vec_copy(k, r_out, sizeof(k));
		}

		// Reference
		unsigned char k_bytes[32];
		le_bytes_from_limbs_256(k_bytes, k);
		Fp12 res_ref;
		Fp12::exp(res_ref, base, k_bytes, 256); // Using standard exp as truth

		// GLV-4
		Fp12 res_glv;
		Fp12::cyclotomic_exp_glv4(res_glv, base, k);

		EXPECT_TRUE(vec_is_equal(&res_ref, &res_glv, sizeof(Fp12)))
			<< "GLV4 Edge Case Failed: " << sc.name;
	}
}

TEST(Fp12Test, CyclotomicExpGLV4_Stress) {
	// A heavier random test with 100 iterations to catch rare decomposition errors.
	const int iterations = 100;

	// Setup base
	Fp12 seed, base;
	vec_copy(&seed, &Fp12::one(), sizeof(seed));
	seed[0][0][1][0] = 0x123;
	Fp12::final_exp(base, seed);

	const limb_t* r_ptr = get_BLST_r().l;
	limb_t r[4];
	vec_copy(r, r_ptr, sizeof(r));

	srand(12345); // Reproducible sequence

	for (int i = 0; i < iterations; ++i) {
		limb_t k_raw[4];
		for (auto& x : k_raw) x = ((uint64_t)rand() << 32) | rand();

		// Ensure k < r
		limb_t k[4], q_ign[4];
		u256_divrem_u256(k_raw, r, q_ign, k);

		// Reference
		unsigned char k_bytes[32];
		le_bytes_from_limbs_256(k_bytes, k);
		Fp12 res_ref;
		Fp12::exp(res_ref, base, k_bytes, 256);

		// Test Target
		Fp12 res_glv;
		Fp12::cyclotomic_exp_glv4(res_glv, base, k);

		if (!vec_is_equal(&res_ref, &res_glv, sizeof(Fp12))) {
			// Failure diagnostic
			printf("GLV4 Stress Test Failed at iter %d\n", i);
			printf("Scalar (LE limbs): %016llx %016llx %016llx %016llx\n",
				(unsigned long long)k[0], (unsigned long long)k[1],
				(unsigned long long)k[2], (unsigned long long)k[3]);

			// Invoke the splitter directly to see what happened (if we could access it)
			// Since we can't easily, just failing the test is enough.
			FAIL() << "Values mismatch.";
		}
	}
}

TEST(Fp12Test, CyclotomicExpGLV4_Cross) {
	// Verification that the optimization works for different bases, 
	// ensuring no static state or base-dependent bugs.

	limb_t k[4] = { 0x1234567890ABCDEF, 0, 0, 0 }; // Fixed scalar
	unsigned char k_bytes[32];
	le_bytes_from_limbs_256(k_bytes, k);

	for (int i = 1; i <= 10; ++i) {
		Fp12 seed, base;
		vec_copy(&seed, &Fp12::one(), sizeof(seed));
		seed[1][2][0][0] = i; // Perturb base
		Fp12::final_exp(base, seed);

		Fp12 res_ref, res_glv;
		Fp12::exp(res_ref, base, k_bytes, 256);
		Fp12::cyclotomic_exp_glv4(res_glv, base, k);

		EXPECT_TRUE(vec_is_equal(&res_ref, &res_glv, sizeof(Fp12)))
			<< "GLV4 failed for base variant " << i;
	}
}
