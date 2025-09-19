#include <gtest/gtest.h>
#include <random>
#include <cstring>

#include <stdint.h>


#include "blst/backend/backend.h"
#include "blst/fields/consts.h"

namespace
{

	// Helper to generate random vec256
	void generate_random_vec256(vec256 v) {
		std::random_device rd;
		std::mt19937_64 gen(rd());
		for (int i = 0; i < NLIMBS(256); ++i) {
			v[i] = gen();
		}
	}


	// Helper to generate random vec256
	void generate_random_vec256(vec256 v, const vec256 mod) {
		std::random_device rd;
		std::mt19937_64 gen(rd());
		for (int i = 0; i < NLIMBS(256)-1; ++i) {
			v[i] = gen();
		}
		v[NLIMBS(256) - 1] = gen() % mod[NLIMBS(256) - 1];
	}




	// Helper to generate random vec384
	void generate_random_vec384(vec384 v, const vec384 p) {
		std::random_device rd;
		std::mt19937_64 gen(rd());
		for (int i = 0; i < NLIMBS(384) - 1; ++i) v[i] = gen();
		// ensure MS limb < p[MS], which guarantees v < p lexicographically
		v[NLIMBS(384) - 1] = gen() % p[NLIMBS(384) - 1];

		assert(cc_check_mod_n<NLIMBS(384)>((byte*)v, p) == 1);
	}

	// Helper to generate random vec384
	void generate_random_vec384(vec384 v) {
		generate_random_vec384(v, *get_BLST_P());
	}

	// Platform-independent random number generator state
	static uint64_t rng_state = 1;

	// Initialize the RNG with a seed
	static void init_rng(uint32_t seed) {
		rng_state = seed;
		if (rng_state == 0) rng_state = 1; // Avoid zero state
	}

	// Platform-independent 64-bit random number generator using xorshift64*
	static uint64_t rand64(void) {
		rng_state ^= rng_state >> 12;
		rng_state ^= rng_state << 25;
		rng_state ^= rng_state >> 27;
		return rng_state * 0x2545F4914F6CDD1DULL;
	}

	// Platform-independent 32-bit random number generator
	//static uint32_t rand32(void) {
	//	return (uint32_t)(rand64() >> 32);
	//}

	// Helper to generate random vec384
	void generate_random_vec384(vec384 v, const vec384 p, int seed) {
		// Initialize our platform-independent RNG
		init_rng(seed);

		// Generate random values for each limb
		for (int i = 0; i < NLIMBS(384) - 1; ++i) {
			v[i] = rand64();
		}

		// Ensure MS limb < p[MS], which guarantees v < p lexicographically
		if (p[NLIMBS(384) - 1] > 0) {
			uint64_t full_rand = rand64();
			v[NLIMBS(384) - 1] = full_rand % p[NLIMBS(384) - 1];
		}
		else {
			v[NLIMBS(384) - 1] = 0;
		}

		// Verify the result is valid
		//assert(check_mod_n((byte*)v, p, NLIMBS(384)) == 1);
	}
	// Helper to generate random vec384x
	void generate_random_vec384x(vec384x v, const vec384 p) {
		generate_random_vec384(v[0], p);
		generate_random_vec384(v[1], p);
	}

	 //Helper to generate random vec384
	void generate_random_vec384x(vec384x v) {
		generate_random_vec384(v[0], *get_BLST_P());
		generate_random_vec384(v[1], *get_BLST_P());
	}

	// Helper to compare two vec256
	bool vec256_equal(const vec256 a, const vec256 b) {
		return std::memcmp(a, b, sizeof(vec256)) == 0;
	}

	// Helper to compare two vec384
	bool vec384_equal(const vec384 a, const vec384 b) {
		return std::memcmp(a, b, sizeof(vec384)) == 0;
	}

	// Add this helper function for printing vectors in hexadecimal format
	void print_vec256(const char* label, const vec256 v) {
		std::cout << label << ": ";
		for (int i = NLIMBS(256) - 1; i >= 0; --i) {
			std::cout << std::hex << std::setfill('0') << std::setw(16) << v[i];
			if (i > 0) std::cout << " ";
		}
		std::cout << std::dec << std::endl;
	}

	void print_vec384(const char* label, const vec384 v) {
		std::cout << label << ": ";
		for (int i = NLIMBS(384) - 1; i >= 0; --i) {
			std::cout << std::hex << std::setfill('0') << std::setw(16) << v[i];
			if (i > 0) std::cout << " ";
		}
		std::cout << std::dec << std::endl;
	}

	void print_vec768(const char* label, const vec768 v) {
		printf("%s: ", label);
		for (int i = NLIMBS(768) - 1; i >= 0; --i) {
			printf("%016llx", (unsigned long long)v[i]);
			if (i > 0) printf(" ");
		}
		printf("\n");
	}


}

// Test suite for BLST_PORTABLE and BLST_ASM implementations
#if defined(BLST_ASM) && defined(BLST_PORTABLE)

TEST(InterfaceTests, MulMontSparseV256) {
	vec256 a, b, p, ret_asm, ret_portable;
	limb_t n0 = BLST_r0;
	int trials = 10;
	memcpy(p, *get_BLST_r(), sizeof(vec256));

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec256(a, p);
		generate_random_vec256(b, p);

		// Compute using both implementations
		mulx_mont_sparse_256(ret_asm, a, b, p, n0); // ASM implementation
		cc_mul_mont_sparse_256(ret_portable, a, b, p, n0); // Portable implementation

		bool equal = vec256_equal(ret_asm, ret_portable);
		EXPECT_TRUE(equal) << "Mismatch in mul_mont_sparse_v256";
	}
}

TEST(InterfaceTests, AddModV384) {
	vec384 a, b, p, ret_asm, ret_portable;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p, *get_BLST_P(), sizeof(vec384));
		generate_random_vec384(a, p);
		generate_random_vec384(b, p);

		// Compute using both implementations
		add_mod_384(ret_asm, a, b, p); // ASM implementation
		cc_add_mod_384(ret_portable, a, b, p); // Portable implementation
		auto eq = vec384_equal(ret_asm, ret_portable);
		if (!eq)
		{
			print_vec384("a", a);
			print_vec384("b", b);
			print_vec384("p", p);
			print_vec384("ret_asm     ", ret_asm);
			print_vec384("ret_portable", ret_portable);
		}

		EXPECT_TRUE(eq) << "Mismatch in add_mod_v384";
	}
}

TEST(InterfaceTests, SqrMontSparseV256) {
	vec256 a, p, ret_asm, ret_portable;
	limb_t n0 = BLST_r0;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec256(a);
		memcpy(p, *get_BLST_r(), sizeof(vec256));

		// Compute using both implementations
		sqrx_mont_sparse_256(ret_asm, a, p, n0); // ASM implementation
		cc_sqr_mont_sparse_256(ret_portable, a, p, n0); // Portable implementation

		EXPECT_TRUE(vec256_equal(ret_asm, ret_portable)) << "Mismatch in sqr_mont_sparse_v256";
	}
}

TEST(InterfaceTests, RedcMontV256) {
	vec512 a;
	vec256 p, ret_asm, ret_portable;
	limb_t n0 = BLST_r0;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec256(reinterpret_cast<vec256&>(a)); // Fill lower 256 bits
		generate_random_vec256(reinterpret_cast<vec256&>(a[4])); // Fill upper 256 bits
		memcpy(p, *get_BLST_r(), sizeof(vec256));
		// Compute using both implementations
		redcx_mont_256(ret_asm, a, p, n0); // ASM implementation
		cc_redc_mont_256(ret_portable, a, p, n0); // Portable implementation
		EXPECT_TRUE(vec256_equal(ret_asm, ret_portable)) << "Mismatch in redc_mont_v256";
	}
}

TEST(InterfaceTests, MulMontV384) {
	vec384 a, b, p, ret_asm, ret_portable;
	limb_t n0 = BLST_p0;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384(a);
		generate_random_vec384(b);
		memcpy(p, *get_BLST_P(), sizeof(vec384));
		// Compute using both implementations
		mulx_mont_384(ret_asm, a, b, p, n0); // ASM implementation
		cc_mul_mont_384(ret_portable, a, b, p, n0); // Portable implementation
		EXPECT_TRUE(vec384_equal(ret_asm, ret_portable)) << "Mismatch in mul_mont_v384";
	}
}

TEST(InterfaceTests, SqrMontV384) {
	vec384 a, p, ret_asm, ret_portable;
	limb_t n0 = BLST_p0;
	int trials = 10;
	for (int t = 0; t < trials; t++)
	{

		generate_random_vec384(a);
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		// Compute using both implementations
		sqrx_mont_384(ret_asm, a, p, n0); // ASM implementation
		cc_sqr_mont_384(ret_portable, a, p, n0); // Portable implementation

		EXPECT_TRUE(vec384_equal(ret_asm, ret_portable)) << "Mismatch in sqr_mont_v384";
	}
}

TEST(InterfaceTests, AddModV256) {
	vec256 a, b, p, ret_asm, ret_portable;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p, *get_BLST_r(), sizeof(vec256));
		generate_random_vec256(a, p);
		generate_random_vec256(b, p);
		// Compute using both implementations
		add_mod_256(ret_asm, a, b, p); // ASM implementation
		cc_add_mod_256(ret_portable, a, b, p); // Portable implementation
		EXPECT_TRUE(vec256_equal(ret_asm, ret_portable)) << "Mismatch in add_mod_v256";
	}
}

TEST(InterfaceTests, SubModV256) {
	vec256 a, b, p, ret_asm, ret_portable;
	int trials = 10;
	for (int t = 0; t < trials; t++)
	{
		memcpy(p, *get_BLST_r(), sizeof(vec256));
		generate_random_vec256(a,p);
		generate_random_vec256(b,p);

		// Compute using both implementations
		sub_mod_256(ret_asm, a, b, p); // ASM implementation
		cc_sub_mod_256(ret_portable, a, b, p); // Portable implementation

		EXPECT_TRUE(vec256_equal(ret_asm, ret_portable)) << "Mismatch in sub_mod_v256";
	}
}

TEST(InterfaceTests, MulBy3ModV256) {
	vec256 a, p, ret_asm, ret_portable;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p, *get_BLST_r(), sizeof(vec256));
		generate_random_vec256(a, p);

		// Compute using both implementations
		mul_by_3_mod_256(ret_asm, a, p); // ASM implementation
		cc_mul_by_3_mod_256(ret_portable, a, p); // Portable implementation

		EXPECT_TRUE(vec256_equal(ret_asm, ret_portable)) << "Mismatch in mul_by_3_mod_v256";
	}
}

TEST(InterfaceTests, CNegModV256) {
	vec256 a, p, ret_asm, ret_portable;
	bool_t flag = true;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p, *get_BLST_r(), sizeof(vec256));
		generate_random_vec256(a,p);
		// Compute using both implementations
		cneg_mod_256(ret_asm, a, flag, p); // ASM implementation
		cc_cneg_mod_256(ret_portable, a, flag, p); // Portable implementation
		EXPECT_TRUE(vec256_equal(ret_asm, ret_portable)) << "Mismatch in cneg_mod_v256";
	}
}

// Test for lshift_mod_v256
TEST(InterfaceTests, LShiftModV256) {
	vec256 a, p, ret_asm, ret_portable;
	size_t count = 3;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p, *get_BLST_r(), sizeof(vec256));
		generate_random_vec256(a,p);

		lshift_mod_256(ret_asm, a, count, p);
		cc_lshift_mod_256(ret_portable, a, count, p);
		auto eq = vec256_equal(ret_asm, ret_portable);
		if (!eq)
		{
			print_vec256("a", a);
			print_vec256("p", p);
			print_vec256("ret_asm", ret_asm);
			print_vec256("ret_portable", ret_portable);
		}
		EXPECT_TRUE(eq) << "Mismatch in lshift_mod_v256";
	}
}

// Test for check_mod_256
TEST(InterfaceTests, CheckMod256) {
	pow256 a;
	vec256 p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		memcpy(p, *get_BLST_r(), sizeof(vec256));
		generate_random_vec256(reinterpret_cast<vec256&>(a));

		// Test both implementations
		limb_t result_asm = check_mod_256(a, p); // ASM implementation
		limb_t result_portable = cc_check_mod_256(a, p); // Portable implementation

		// Compare results
		EXPECT_EQ(result_asm, result_portable) << "Mismatch in check_mod_256";
	}
}

// Test for add_n_check_mod_256
TEST(InterfaceTests, AddNCheckMod256) {
	pow256 a, b, ret_asm, ret_portable;
	vec256 p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		memcpy(p, *get_BLST_r(), sizeof(vec256));
		generate_random_vec256(reinterpret_cast<vec256&>(a),p);
		generate_random_vec256(reinterpret_cast<vec256&>(b),p);

		// Make copies for both implementations
		memcpy(ret_asm, a, sizeof(pow256));
		memcpy(ret_portable, a, sizeof(pow256));

		// Test both implementations
		limb_t result_asm = add_n_check_mod_256(ret_asm, a, b, p); // ASM implementation
		limb_t result_portable = cc_add_n_check_mod_256(ret_portable, a, b, p); // Portable implementation

		// Compare return values
		EXPECT_EQ(result_asm, result_portable) << "Mismatch in add_n_check_mod_256 return value";

		// Compare modified results
		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(pow256)) == 0)
			<< "Mismatch in add_n_check_mod_256 output";
	}
}

// Test for sub_n_check_mod_256
TEST(InterfaceTests, SubNCheckMod256) {
	pow256 a, b, ret_asm, ret_portable;
	vec256 p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		generate_random_vec256(reinterpret_cast<vec256&>(a));
		generate_random_vec256(reinterpret_cast<vec256&>(b));
		memcpy(p, *get_BLST_r(), sizeof(vec256));

		// Make copies for both implementations
		memcpy(ret_asm, a, sizeof(pow256));
		memcpy(ret_portable, a, sizeof(pow256));

		// Test both implementations
		limb_t result_asm = sub_n_check_mod_256(ret_asm, a, b, p); // ASM implementation
		limb_t result_portable = cc_sub_n_check_mod_256(ret_portable, a, b, p); // Portable implementation

		// Compare return values
		EXPECT_EQ(result_asm, result_portable) << "Mismatch in sub_n_check_mod_256 return value";

		// Compare modified results
		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(pow256)) == 0)
			<< "Mismatch in sub_n_check_mod_256 output";
	}
}

// Test for rshift_mod_v256
TEST(InterfaceTests, RShiftModV256) {
	vec256 a, p, ret_asm, ret_portable;
	size_t count = 3;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec256(a);
		memcpy(p, *get_BLST_r(), sizeof(vec256));

		rshift_mod_256(ret_asm, a, count, p);
		cc_rshift_mod_256(ret_portable, a, count, p);

		EXPECT_TRUE(vec256_equal(ret_asm, ret_portable)) << "Mismatch in rshift_mod_v256";
	}
}

// Test for from_mont_v256
TEST(InterfaceTests, FromMontV256) {
	vec256 a, p, ret_asm, ret_portable;
	limb_t n0 = BLST_r0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec256(a);
		memcpy(p, *get_BLST_r(), sizeof(vec256));

		// Compute using both implementations
		fromx_mont_256(ret_asm, a, p, n0); // ASM implementation
		cc_from_mont_256(ret_portable, a, p, n0); // Portable implementation

		EXPECT_TRUE(vec256_equal(ret_asm, ret_portable)) << "Mismatch in from_mont_v256";
	}
}

// Test for from_mont_v384
TEST(InterfaceTests, FromMontV384) {
	vec384 a, p, ret_asm, ret_portable;
	limb_t n0 = BLST_p0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384(a);
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		// Compute using both implementations
		fromx_mont_384(ret_asm, a, p, n0); // ASM implementation
		cc_from_mont_384(ret_portable, a, p, n0); // Portable implementation

		EXPECT_TRUE(vec384_equal(ret_asm, ret_portable)) << "Mismatch in from_mont_v384";
	}
}
//// Test for sqr_n_mul_mont_v384
//TEST(InterfaceTests, SqrNMulMontV384) {
//	vec384 a, b, p, ret_asm, ret_portable;
//	limb_t n0 = BLST_p0;
//	size_t count = 1; // Number of squaring operations
//
//	int trials = 20;
//
//	for (int t = 0; t < trials; t++)
//	{
//		memcpy(p, *get_BLST_P(), sizeof(vec384));
//		generate_random_vec384(a, p);
//		generate_random_vec384(b, p);
//		//memcpy(p, *get_BLST_P(), sizeof(vec384));
//
//		// Compute using both implementations
//		sqrx_n_mul_mont_384(ret_asm, a, count, p, n0, b); // ASM implementation
//		cc_sqr_n_mul_mont_384(ret_portable, a, count, p, n0, b); // Portable implementation
//
//		auto eq = vec384_equal(ret_asm, ret_portable);
//		if (eq == false)
//		{
//			// If they are not equal, print the vectors for debugging
//			std::cout << "Mismatch in sqr_n_mul_mont_v384" << std::endl;
//			print_vec384("a", a);
//			print_vec384("b", b);
//			print_vec384("p", p);
//
//
//			print_vec384("ret_asm     ", ret_asm);
//			print_vec384("ret_portable", ret_portable);
//
//		}
//		EXPECT_TRUE(eq) << "Mismatch in sqr_n_mul_mont_v384";
//	}
//}

//static inline uint64_t mont_n0_from_p0(uint64_t _p0) {
//	// p0 must be odd (true for prime moduli)
//	assert(_p0 & 1);
//
//	// Newton iteration for inverse mod 2^64
//	uint64_t x = 1;                  // initial approx to p0^{-1} mod 2
//	x *= 2 - _p0 * x;                 // 2^2
//	x *= 2 - _p0 * x;                 // 2^4
//	x *= 2 - _p0 * x;                 // 2^8
//	x *= 2 - _p0 * x;                 // 2^16
//	x *= 2 - _p0 * x;                 // 2^32
//	x *= 2 - _p0 * x;                 // 2^64
//
//	// x == p0^{-1} mod 2^64; we want n0 = -(p0^{-1}) mod 2^64
//	return (uint64_t)(0 - x);
//}
//
//static inline uint64_t mont_n0_from_p(vec384 p_le) {
//	return mont_n0_from_p0(p_le[0]); // little-endian limbs: p[0] is LSW
//}

// Test for sqr_n_mul_mont_v383
TEST(InterfaceTests, SqrNMulMontV383) {
	vec384 a, b, aM, bM, p, ret_asm, ret_portable, out_asm, out_c;
	memcpy(p, *get_BLST_P(), sizeof(vec384));

	limb_t n0 = BLST_p0; 
	//auto nn0 = mont_n0_from_p(*get_BLST_P());
	//assert(n0 == nn0);
	size_t count = 5; // Number of squaring operations

	int trials = 40;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384(a,p);
		generate_random_vec384(b,p);

		EXPECT_TRUE(cc_check_mod_n<NLIMBS(384)>((byte*)b, p) == 1);
		//to_fp(a, a);
		//to_fp(b, b);
		// ---- to Montgomery correctly: aM = a * R  and  bM = b * R
		// via MontMul with RR
		cc_mul_mont_n<NLIMBS(384)>(aM, a, *get_BLST_RR(), p, n0);
		cc_mul_mont_n<NLIMBS(384)>(bM, b, *get_BLST_RR(), p, n0);


		// Compute using both implementations
		sqrx_n_mul_mont_383(ret_asm, a, count, p, n0, b); // ASM implementation
		cc_sqr_n_mul_mont_384(ret_portable, a, count, p, n0, b); // Portable implementation

		// ---- back from Montgomery for comparison
		cc_from_mont_384(out_asm, ret_asm, p, n0);
		cc_from_mont_384(out_c, ret_portable, p, n0);


		EXPECT_TRUE(vec384_equal(out_asm, out_c)) << "Mismatch in sqr_n_mul_mont_v383";
	}
}

TEST(InterfaceTests, SubModV384) {
	vec384 a, b, p, ret_asm, ret_portable;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p, *get_BLST_P(), sizeof(vec384));
		generate_random_vec384(a,p);
		generate_random_vec384(b,p);

		sub_mod_384(ret_asm, a, b, p);
		cc_sub_mod_n<NLIMBS(384)>(ret_portable, a, b, p);

		EXPECT_TRUE(vec384_equal(ret_asm, ret_portable)) << "Mismatch in sub_mod_v384";
	}
}

TEST(InterfaceTests, MulBy8ModV384) {
	vec384 a, p, ret_asm, ret_portable;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p, *get_BLST_P(), sizeof(vec384));
		generate_random_vec384(a,p);

		mul_by_8_mod_384(ret_asm, a, p);
		cc_mul_by_8_mod_384(ret_portable, a, p);

		EXPECT_TRUE(vec384_equal(ret_asm, ret_portable)) << "Mismatch in mul_by_8_mod_v384";
	}
}

TEST(InterfaceTests, MulBy3ModV384) {
	vec384 a, p, ret_asm, ret_portable;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p, *get_BLST_P(), sizeof(vec384));
		generate_random_vec384(a, p);

		mul_by_3_mod_384(ret_asm, a, p);
		cc_mul_by_3_mod_384(ret_portable, a, p);

		EXPECT_TRUE(vec384_equal(ret_asm, ret_portable)) << "Mismatch in mul_by_3_mod_v384";
	}
}

TEST(InterfaceTests, CNegModV384) {
	vec384 a, p, ret_asm, ret_portable;
	bool_t flag = true;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p, *get_BLST_P(), sizeof(vec384));
		generate_random_vec384(a,p);

		cneg_mod_384(ret_asm, a, flag, p);
		cc_cneg_mod_384(ret_portable, a, flag, p);

		EXPECT_TRUE(vec384_equal(ret_asm, ret_portable)) << "Mismatch in cneg_mod_v384";
	}
}


TEST(InterfaceTests, RShiftMod384) {
	vec384 a, p, ret_asm, ret_portable;
	size_t count = 3;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p, *get_BLST_P(), sizeof(vec384));
		generate_random_vec384(a,p);

		rshift_mod_384(ret_asm, a, count, p); // ASM implementation
		cc_rshift_mod_384(ret_portable, a, count, p); // Portable implementation

		EXPECT_TRUE(vec384_equal(ret_asm, ret_portable)) << "Mismatch in rshift_mod_384";
	}
}

// Test for redcx_mont_384
TEST(InterfaceTests, RedcxMont384) {
	vec768 a;
	vec384 p, ret_asm, ret_portable;
	limb_t n0 = BLST_p0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		for (int i = 0; i < NLIMBS(768); ++i) {
			a[i] = rand();
		}
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		// Compute using both implementations
		redcx_mont_384(ret_asm, a, p, n0); // ASM implementation
		cc_redc_mont_384(ret_portable, a, p, n0); // Portable implementation

		EXPECT_TRUE(vec384_equal(ret_asm, ret_portable)) << "Mismatch in redcx_mont_384";
	}
}

// Test for redc_mont_384
TEST(InterfaceTests, RedcMont384) {
#if !defined(BLST_ADX) || defined(__MSC_VER)
	vec768 a;
	vec384 p, ret_asm, ret_portable;
	limb_t n0 = BLST_p0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		for (int i = 0; i < NLIMBS(768); ++i) {
			a[i] = rand();
		}
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		// Compute using both implementations
		redc_mont_384(ret_asm, a, p, n0); // ASM implementation
		cc_redc_mont_384(ret_portable, a, p, n0); // Portable implementation

		EXPECT_TRUE(vec384_equal(ret_asm, ret_portable)) << "Mismatch in redc_mont_384";
	}
#endif
}

TEST(InterfaceTests, DivBy2Mod384) {
	vec384 a, p, ret_asm, ret_portable;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384(a);
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		div_by_2_mod_384(ret_asm, a, p); // ASM implementation
		cc_div_by_2_mod_384(ret_portable, a, p); // Portable implementation

		EXPECT_TRUE(vec384_equal(ret_asm, ret_portable)) << "Mismatch in div_by_2_mod_384";
	}
}


static inline unsigned clz64(uint64_t x) {
	if (x == 0) return 64;
	unsigned n = 0;
	if ((x >> 32) == 0) { n += 32; x <<= 32; }
	if ((x >> 48) == 0) { n += 16; x <<= 16; }
	if ((x >> 56) == 0) { n += 8;  x <<= 8; }
	if ((x >> 60) == 0) { n += 4;  x <<= 4; }
	if ((x >> 62) == 0) { n += 2;  x <<= 2; }
	if ((x >> 63) == 0) { n += 1; }
	return n;
}

static inline unsigned bitlen_be(const limb_t* a, size_t n) {
	// find highest nonzero limb
	for (size_t i = n; i-- > 0; ) if (a[i]) {
		return (unsigned)(64 * i + (unsigned)(64 - clz64(a[i])));
	}
	return 0;
}

static inline void left_align_mod(limb_t* modx, const limb_t* mod, size_t n) {
	unsigned bl = bitlen_be(mod, n);
	unsigned sh = (64 - (bl & 63)) & 63; // BLST_LIMB_SIZE == 64
	// modx = mod << sh
	limb_t carry = 0;
	for (size_t i = 0; i < n; i++) {
		limb_t x = mod[i];
		modx[i] = (sh ? (x << sh) | (carry >> (64 - sh)) : x);
		carry = x;
	}
}

const vec384 modx = {    /* left-aligned value of the modulus */
	TO_LIMB_T(0xcff7fffffffd5558), TO_LIMB_T(0xf55ffff58a9ffffd),
	TO_LIMB_T(0x39869507b587b120), TO_LIMB_T(0x23ba5c279c2895fb),
	TO_LIMB_T(0x58dd3db21a5d66bb), TO_LIMB_T(0xd0088f51cbff34d2)
};

TEST(InterfaceTests, CtxInverseMod384) {


	vec384 inp;
	vec768 ret_asm, ret_portable;

	{
		// Verify modx is correctly left-aligned
		vec384 modx_check;
		left_align_mod(modx_check, *get_BLST_P(), NLIMBS(384));
		assert(std::memcmp(modx, modx_check, sizeof(vec384)) == 0);
	}

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384(inp, *get_BLST_P());

		ctx_inverse_mod_384(ret_asm, inp, *get_BLST_P(), modx); // ASM implementation
		cc_ct_inverse_mod_384(ret_portable, inp, *get_BLST_P(), modx); // Portable implementation
		
		
		const size_t n = NLIMBS(384);

		// handy views
		auto* asm_hi = ret_asm + n;
		auto* por_hi = ret_portable + n;

		// 0) High-half equality (raw 2n output)
		int eq0 = std::memcmp(asm_hi, por_hi, sizeof(vec384));

		// 1) REDC stage equality
		union { vec768 T; vec384 r[2]; } A, B;
		std::memcpy(A.T, ret_asm, sizeof(vec768));
		std::memcpy(B.T, ret_portable, sizeof(vec768));
		vec384 redc_asm, redc_por;
		redc_mont_v384(redc_asm, A.T, *get_BLST_P(), BLST_p0);
		redc_mont_v384(redc_por, B.T, *get_BLST_P(), BLST_p0);
		int eq1 = std::memcmp(redc_asm, redc_por, sizeof(vec384));

		// 2) Convert-to-Montgomery equality
		vec384 inv_asm, inv_por;
		mul_mont_v384(inv_asm, redc_asm, *get_BLST_RR(), *get_BLST_P(), BLST_p0);
		mul_mont_v384(inv_por, redc_por, *get_BLST_RR(), *get_BLST_P(), BLST_p0);
		int eq2 = std::memcmp(inv_asm, inv_por, sizeof(vec384));

		// 3) Functional correctness (both sides!)
		vec384 oneA, oneP;
		if (!vec_is_zero(inp, sizeof(vec384))) {
			mul_mont_v384(oneA, inv_asm, inp, *get_BLST_P(), BLST_p0);
			mul_mont_v384(oneP, inv_por, inp, *get_BLST_P(), BLST_p0);
			EXPECT_TRUE(vec_is_equal(oneA, get_BLST_Rx()->p, sizeof(vec384))) << "ASM inverse incorrect";
			EXPECT_TRUE(vec_is_equal(oneP, get_BLST_Rx()->p, sizeof(vec384))) << "Portable inverse incorrect";
		}

		// Assert stage-by-stage (these help pinpoint the break)
		EXPECT_EQ(eq0, 0) << "Mismatch in raw high halves (ret+n)";
		EXPECT_EQ(eq1, 0) << "Mismatch after REDC";
		EXPECT_EQ(eq2, 0) << "Mismatch after convert-to-Mont";
	}
}

TEST(InterfaceTests, CtInverseMod256) {
	vec256 inp;
	vec512 ret_asm, ret_portable;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec256(inp);

		ct_inverse_mod_256(ret_asm, inp, *get_BLST_P(), modx); // ASM implementation
		cc_ct_inverse_mod_256(ret_portable, inp, *get_BLST_P(), modx); // Portable implementation

		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(vec512)) == 0) << "Mismatch in ct_inverse_mod_256";
	}
}

TEST(InterfaceTests, CtIsSquareMod384) {
	vec384 inp;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384(inp, *get_BLST_P(), 24323421 * t + 2142727192);

		bool_t is_square_asm = ct_is_square_mod_384(inp, *get_BLST_P()); // ASM implementation
		bool_t is_square_portable = cc_ct_is_square_mod_384(inp, *get_BLST_P()); // Portable implementation

		//printf("passed? %d\n", is_square_asm == is_square_portable);
		//print_vec384("inp", inp);
		//printf("is_square_asm     = %lld\n", is_square_asm);
		//printf("is_square_portable= %lld\n", is_square_portable);

		EXPECT_EQ(is_square_asm, is_square_portable) << "Mismatch in ct_is_square_mod_384";
	}
}

TEST(InterfaceTests, MulxMont384x) {
	vec384x a, b, ret_asm, ret_portable;
	vec384 p;
	limb_t n0 = BLST_p0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384x(a, *get_BLST_P());
		generate_random_vec384x(b, *get_BLST_P());
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		mulx_mont_384x(ret_asm, a, b, p, n0); // ASM implementation
		cc_mul_mont_384x(ret_portable, a, b, p, n0); // Portable implementation

		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(vec384x)) == 0) << "Mismatch in mulx_mont_384x";
	}
}

TEST(InterfaceTests, SqrxMont384x) {
	vec384x a, ret_asm, ret_portable;
	vec384 p;
	limb_t n0 = BLST_p0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384x(a, *get_BLST_P());
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		sqrx_mont_384x(ret_asm, a, p, n0); // ASM implementation
		cc_sqr_mont_384x(ret_portable, a, p, n0); // Portable implementation

		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(vec384x)) == 0) << "Mismatch in sqrx_mont_384x";
	}
}

TEST(InterfaceTests, SqrxMont382x) {
	vec384x a, ret_asm, ret_portable;
	vec384 p;
	limb_t n0 = BLST_p0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384x(a, *get_BLST_P());
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		sqrx_mont_382x(ret_asm, a, p, n0); // ASM implementation
		cc_sqr_mont_382x(ret_portable, a, p, n0); // Portable implementation

		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(vec384x)) == 0) << "Mismatch in sqrx_mont_382x";
	}
}

TEST(InterfaceTests, AddMod384x) {
	vec384x a, b, ret_asm, ret_portable;
	vec384 p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384x(a, *get_BLST_P());
		generate_random_vec384x(b, *get_BLST_P());
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		add_mod_384x(ret_asm, a, b, p); // ASM implementation
		cc_add_mod_384x(ret_portable, a, b, p); // Portable implementation

		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(vec384x)) == 0) << "Mismatch in add_mod_384x";
	}
}

TEST(InterfaceTests, SubMod384x) {
	vec384x a, b, ret_asm, ret_portable;
	vec384 p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384x(a, *get_BLST_P());
		generate_random_vec384x(b, *get_BLST_P());
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		sub_mod_384x(ret_asm, a, b, p); // ASM implementation
		cc_sub_mod_384x(ret_portable, a, b, p); // Portable implementation

		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(vec384x)) == 0) << "Mismatch in sub_mod_384x";
	}
}

// Test for mul_by_8_mod_384x
TEST(InterfaceTests, MulBy8Mod384x) {
	vec384x a, ret_asm, ret_portable;
	vec384 p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		generate_random_vec384(reinterpret_cast<vec384&>(a[0]));
		generate_random_vec384(reinterpret_cast<vec384&>(a[1]));
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		// Test both implementations
		mul_by_8_mod_384x(ret_asm, a, p); // ASM implementation
		cc_mul_by_8_mod_384x(ret_portable, a, p); // Portable implementation

		// Compare results
		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(vec384x)) == 0) << "Mismatch in mul_by_8_mod_384x";
	}
}

// Test for mul_by_3_mod_384x
TEST(InterfaceTests, MulBy3Mod384x) {
	vec384x a, ret_asm, ret_portable;
	vec384 p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		generate_random_vec384(reinterpret_cast<vec384&>(a[0]));
		generate_random_vec384(reinterpret_cast<vec384&>(a[1]));
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		// Test both implementations
		mul_by_3_mod_384x(ret_asm, a, p); // ASM implementation
		cc_mul_by_3_mod_384x(ret_portable, a, p); // Portable implementation

		// Compare results
		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(vec384x)) == 0) << "Mismatch in mul_by_3_mod_384x";
	}
}

// Test for mul_by_1_plus_i_mod_384x
TEST(InterfaceTests, MulBy1PlusIMod384x) {
	vec384x a, ret_asm, ret_portable;
	vec384 p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		generate_random_vec384(reinterpret_cast<vec384&>(a[0]));
		generate_random_vec384(reinterpret_cast<vec384&>(a[1]));
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		// Test both implementations
		mul_by_1_plus_i_mod_384x(ret_asm, a, p); // ASM implementation
		cc_mul_by_1_plus_i_mod_384x(ret_portable, a, p); // Portable implementation

		// Compare results
		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(vec384x)) == 0) << "Mismatch in mul_by_1_plus_i_mod_384x";
	}
}

// Test for sgn0x_pty_mont_384
TEST(InterfaceTests, Sgn0PtyMont384) {
	vec384 a, p;
	limb_t n0 = BLST_p0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384(a);
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		// Compute using both implementations
		limb_t result_asm = sgn0x_pty_mont_384(a, p, n0); // ASM implementation
		limb_t result_portable = cc_sgn0_pty_mont_384(a, p, n0); // Portable implementation

		EXPECT_EQ(result_asm, result_portable) << "Mismatch in sgn0x_pty_mont_384";
	}
}

// Test for sgn0x_pty_mont_384x (complex number version)
TEST(InterfaceTests, Sgn0PtyMont384x) {
	vec384x a;
	vec384 p;
	limb_t n0 = BLST_p0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384x(a);
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		// Compute using both implementations
		limb_t result_asm = sgn0x_pty_mont_384x(a, p, n0); // ASM implementation
		limb_t result_portable = cc_sgn0_pty_mont_384x(a, p, n0); // Portable implementation

		EXPECT_EQ(result_asm, result_portable) << "Mismatch in sgn0x_pty_mont_384x";
	}
}

// Test for sgn0_pty_mod_384
TEST(InterfaceTests, Sgn0PtyMod384) {
	vec384 a, p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384(a);
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		// Compute using both implementations
		limb_t result_asm = sgn0_pty_mod_384(a, p); // ASM implementation
		limb_t result_portable = cc_sgn0_pty_mod_384(a, p); // Portable implementation

		EXPECT_EQ(result_asm, result_portable) << "Mismatch in sgn0_pty_mod_384";
	}
}

// Test for sgn0_pty_mod_384x (complex number version)
TEST(InterfaceTests, Sgn0PtyMod384x) {
	vec384x a;
	vec384 p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_vec384x(a);
		memcpy(p, *get_BLST_P(), sizeof(vec384));

		// Compute using both implementations
		limb_t result_asm = sgn0_pty_mod_384x(a, p); // ASM implementation
		limb_t result_portable = cc_sgn0_pty_mod_384x(a, p); // Portable implementation

		EXPECT_EQ(result_asm, result_portable) << "Mismatch in sgn0_pty_mod_384x";
	}
}

// Test for SHA256 functions
TEST(InterfaceTests, SHA256Functions) {
	// Test blst_sha256_block_data_order
	const size_t blocks = 2;
	const size_t block_size = 64; // SHA256 block size in bytes
	unsigned int v_asm[8], v_portable[8]; // SHA256 state
	unsigned char input[blocks * block_size];

	// Initialize state and input data
	for (size_t i = 0; i < 8; ++i) {
		v_asm[i] = v_portable[i] = 0x12345678 + (unsigned int)i;
	}

	// Fill input with random data
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, 255);
	for (size_t i = 0; i < blocks * block_size; ++i) {
		input[i] = static_cast<unsigned char>(dis(gen));
	}

	// Process blocks with both implementations
	blst_sha256_block_data_order(v_asm, input, blocks); // ASM implementation
	cc_blst_sha256_block_data_order(v_portable, input, blocks); // Portable implementation

	// Compare results
	EXPECT_TRUE(std::memcmp(v_asm, v_portable, sizeof(v_asm)) == 0) << "Mismatch in blst_sha256_block_data_order";

	// Test blst_sha256_hcopy
	unsigned int dst_asm[8], dst_portable[8];

	blst_sha256_hcopy(dst_asm, v_asm); // ASM implementation
	cc_blst_sha256_hcopy(dst_portable, v_portable); // Portable implementation

	EXPECT_TRUE(std::memcmp(dst_asm, dst_portable, sizeof(dst_asm)) == 0) << "Mismatch in blst_sha256_hcopy";

	// Test blst_sha256_emit
	unsigned char md_asm[32], md_portable[32];

	blst_sha256_emit(md_asm, dst_asm); // ASM implementation
	cc_blst_sha256_emit(md_portable, dst_portable); // Portable implementation

	EXPECT_TRUE(std::memcmp(md_asm, md_portable, sizeof(md_asm)) == 0) << "Mismatch in blst_sha256_emit";

	// Test blst_sha256_bcopy
	unsigned char dst_copy_asm[32], dst_copy_portable[32];

	blst_sha256_bcopy(dst_copy_asm, md_asm, sizeof(md_asm)); // ASM implementation
	cc_blst_sha256_bcopy(dst_copy_portable, md_portable, sizeof(md_portable)); // Portable implementation

	EXPECT_TRUE(std::memcmp(dst_copy_asm, dst_copy_portable, sizeof(dst_copy_asm)) == 0) << "Mismatch in blst_sha256_bcopy";
}

//// Test for division functions
//TEST(InterfaceTests, DivisionFunctions) {
//	// Test div_3_limbs
//	limb_t dividend_top[2] = { 0xABCDEF0123456789, 0x0123456789ABCDEF };
//	limb_t d_lo = 0x1111111111111111;
//	limb_t d_hi = 0x2222222222222222;
//
//	limb_t result_asm = div_3_limbs(dividend_top, d_lo, d_hi); // ASM implementation
//	limb_t result_portable = cc_div_3_limbs(dividend_top, d_lo, d_hi); // Portable implementation
//
//	EXPECT_EQ(result_asm, result_portable) << "Mismatch in div_3_limbs";
//
//	// Test quot_rem_128
//	limb_t quot_rem_asm[4] = { 0xFEDCBA9876543210, 0x0123456789ABCDEF };
//	limb_t quot_rem_portable[4] = { 0xFEDCBA9876543210, 0x0123456789ABCDEF };
//	limb_t divisor[2] = { 0x1111111111111111, 0x2222222222222222 };
//	limb_t quotient = 0x3333333333333333;
//
//	limb_t rem_asm = quot_rem_128(quot_rem_asm, divisor, quotient); // ASM implementation
//	limb_t rem_portable = cc_quot_rem_128(quot_rem_portable, divisor, quotient); // Portable implementation
//
//	EXPECT_EQ(rem_asm, rem_portable) << "Mismatch in quot_rem_128 remainder";
//	EXPECT_TRUE(std::memcmp(quot_rem_asm, quot_rem_portable, sizeof(quot_rem_asm)) == 0) << "Mismatch in quot_rem_128 quotient";
//
//	// Test quot_rem_64
//	limb_t quot_rem_64_asm[2] = { 0xFEDCBA9876543210 };
//	limb_t quot_rem_64_portable[2] = { 0xFEDCBA9876543210 };
//	limb_t divisor_64[1] = { 0x1111111111111111 };
//	limb_t quotient_64 = 0x2222222222222222;
//
//	limb_t rem_64_asm = quot_rem_64(quot_rem_64_asm, divisor_64, quotient_64); // ASM implementation
//	limb_t rem_64_portable = cc_quot_rem_64(quot_rem_64_portable, divisor_64, quotient_64); // Portable implementation
//
//	EXPECT_EQ(rem_64_asm, rem_64_portable) << "Mismatch in quot_rem_64 remainder";
//	EXPECT_TRUE(std::memcmp(quot_rem_64_asm, quot_rem_64_portable, sizeof(quot_rem_64_asm)) == 0) << "Mismatch in quot_rem_64 quotient";
//}
//

// ---- Helpers ----

// Ensure divisor is left-aligned: MSB(d_hi) == 1
static inline std::pair<limb_t, limb_t> make_left_aligned(limb_t hi, limb_t lo) {
	if ((hi >> 63) == 0) hi |= (1ull << 63);
	return { hi, lo };
}

static inline std::pair<limb_t, limb_t> random_left_aligned_divisor(std::mt19937_64& rng) {
	std::uniform_int_distribution<limb_t> dist;
	limb_t d_hi = dist(rng) | (1ull << 63);
	limb_t d_lo = dist(rng);
	return { d_hi, d_lo };
}

// Sample R in [0, D] to satisfy the stated precondition “dividend top limbs not larger than divisor”
static inline std::pair<limb_t, limb_t> random_R_leq_D(std::mt19937_64& rng, limb_t d_hi, limb_t d_lo) {
	std::uniform_int_distribution<limb_t> dist;
	while (true) {
		limb_t r_hi = dist(rng);
		limb_t r_lo = dist(rng);
		if (r_hi < d_hi || (r_hi == d_hi && r_lo <= d_lo)) return { r_hi, r_lo };
	}
}

// Thin wrappers so we don’t duplicate the div_top packing or change call order by mistake.
static inline limb_t uut_call(limb_t r_hi, limb_t r_lo, limb_t d_lo, limb_t d_hi) {
	limb_t div_top[2] = { r_lo, r_hi };
	return cc_div_3_limbs(div_top, d_lo, d_hi);
}
static inline limb_t ref_call(limb_t r_hi, limb_t r_lo, limb_t d_lo, limb_t d_hi) {
	limb_t div_top[2] = { r_lo, r_hi };
	return div_3_limbs(div_top, d_lo, d_hi);
}

// --- Tests ---

TEST(InterfaceTests, Div3LimbsEdgeCases) {
	// D = minimal left-aligned (1<<127), R in {0, D-1, D}
	auto D = make_left_aligned(0, 0);
	limb_t Dhi = D.first, Dlo = D.second;

	struct C { limb_t r_hi, r_lo; const char* name; };
	std::vector<C> cases;

	cases.push_back({ 0, 0, "R=0" });
	// R = D - 1 (handle borrow if Dlo==0)
	{
		limb_t Rhi = Dhi, Rlo = Dlo - 1;
		if (Rlo > Dlo) Rhi -= 1;
		cases.push_back({ Rhi, Rlo, "R=D-1" });
	}
	cases.push_back({ Dhi, Dlo, "R=D" });

	for (auto& tc : cases) {
		limb_t q_uut = uut_call(tc.r_hi, tc.r_lo, Dlo, Dhi);
		limb_t q_ref = ref_call(tc.r_hi, tc.r_lo, Dlo, Dhi);
		ASSERT_EQ(q_uut, q_ref) << "Edge mismatch: " << tc.name
			<< " D=(" << std::hex << Dhi << ":" << Dlo << ")"
			<< " R=(" << tc.r_hi << ":" << tc.r_lo << ")"
			<< " q_uut=" << q_uut << " q_ref=" << q_ref << std::dec;
	}
}

TEST(InterfaceTests, Div3LimbsRandomFuzz) {
	std::mt19937_64 rng(0xC0FFEE);
	for (int t = 0; t < 50000; ++t) {
		auto D = random_left_aligned_divisor(rng);
		limb_t d_hi = D.first, d_lo = D.second;
		auto R = random_R_leq_D(rng, d_hi, d_lo);
		limb_t r_hi = R.first, r_lo = R.second;

		limb_t q_uut = uut_call(r_hi, r_lo, d_lo, d_hi);
		limb_t q_ref = ref_call(r_hi, r_lo, d_lo, d_hi);

		if (q_uut != q_ref) {
			ADD_FAILURE() << "Trial " << t << " mismatch\n"
				<< "D=(" << std::hex << d_hi << ":" << d_lo << ") "
				<< "R=(" << r_hi << ":" << r_lo << ")\n"
				<< "q_uut=" << q_uut << " q_ref=" << q_ref << std::dec;
			break;
		}
	}
}


#endif // defined(BLST_ASM) && defined(BLST_PORTABLE)