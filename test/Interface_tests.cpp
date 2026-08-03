#include <gtest/gtest.h>
#include <random>
#include <cstring>

#include <stdint.h>

#include "blst/fields/Fr.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/backend/backend.h"
#include "InsecureRng.h"

using namespace blst;

#if defined(__x86_64__) || defined(__x86_64) || \
	(defined(_M_X64) && !defined(_M_ARM64EC))
# define CATBLST_TEST_X86_64 1
#else
# define CATBLST_TEST_X86_64 0
#endif

TEST(InterfaceTests, NumBitsHandlesZeroAndBoundaries) {
	EXPECT_EQ(cc_num_bits(0), 0u);
	EXPECT_EQ(cc_num_bits(1), 1u);
	EXPECT_EQ(cc_num_bits(limb_t(1) << 63), 64u);
	EXPECT_EQ(cc_num_bits((limb_t(1) << 63) - 1), 63u);
}

TEST(InterfaceTests, QuotRemPropagatesBorrowPastMaxLimb) {
	const limb_t max = ~limb_t(0);
	const limb_t divisor[2] = { 1, 1 };
	limb_t dividend_and_remainder[3] = { 0, max, 0 };

	const limb_t corrected =
		cc_quot_rem_n<2>(dividend_and_remainder, divisor, max);

	EXPECT_EQ(corrected, max - 1);
	EXPECT_EQ(dividend_and_remainder[0], 2u);
	EXPECT_EQ(dividend_and_remainder[1], 0u);
	EXPECT_EQ(dividend_and_remainder[2], max - 1);
}

TEST(InterfaceTests, EmptyVectorPredicates) {
	EXPECT_TRUE(vec_is_zero(nullptr, 0));
	EXPECT_TRUE(vec_is_equal(nullptr, nullptr, 0));
}

TEST(InterfaceTests, VectorCopyAcceptsExactSelfAlias) {
	unsigned char value[17] = {
		0, 1, 2, 3, 4, 5, 6, 7, 8,
		9, 10, 11, 12, 13, 14, 15, 16
	};
	unsigned char expected[sizeof(value)];
	std::memcpy(expected, value, sizeof(value));

	vec_copy(value, value, sizeof(value));

	EXPECT_EQ(std::memcmp(value, expected, sizeof(value)), 0);
}

TEST(InterfaceTests, ZeroLengthSelectedShaOperationsAreNoOps) {
	unsigned int state[8] = {
		0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
		0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U
	};
	unsigned int expected[8];
	std::memcpy(expected, state, sizeof(state));

	blst_sha256_block_data_order_v(state, nullptr, 0);
	blst_sha256_bcopy_v(nullptr, nullptr, 0);

	EXPECT_EQ(std::memcmp(state, expected, sizeof(state)), 0);
}

namespace
{

	// Helper to generate random Fr
	void generate_random_Fr(Fr& v) {
		std::random_device rd;
		std::mt19937_64 gen(rd());
		for (int i = 0; i < NLIMBS(256); ++i) {
			v.l[i] = gen();
		}
	}


	// Helper to generate random Fr&
	void generate_random_Fr(Fr& v, const Fr& mod) {
		std::random_device rd;
		std::mt19937_64 gen(rd());
		for (int i = 0; i < NLIMBS(256) - 1; ++i) {
			v[i] = gen();
		}
		v[NLIMBS(256) - 1] = gen() % mod[NLIMBS(256) - 1];
	}

	// Helper to generate random Fp
	void generate_random_Fp(Fp& v, InsecureRng& rng) {

		if (0)
		{

			// Generate random values for each limb_t
			for (int i = 0; i < NLIMBS(384) - 1; ++i) {
				v[i] = rng.next_u64();
			}

			// Ensure MS limb_t < p[MS], which guarantees v < p lexicographically
			if (Fp::modulus()[NLIMBS(384) - 1] > 0) {
				uint64_t full_rand = rng.next_u64();
				v[NLIMBS(384) - 1] = full_rand % Fp::modulus()[NLIMBS(384) - 1];
			}
			else {
				v[NLIMBS(384) - 1] = 0;
			}
		}
		else
		{
			for (u64 i = 0; i < v.num_limb(); ++i)
			{
				v[i] = rng.next_u64();
			}

			while (Fp::is_redc(v) == 0)
			{
				Fp::rshift(v, v, 1);
			}
		}
		// Verify the result is valid
		//assert(check_mod_n((byte*)v, p, NLIMBS(384)) == 1);
	}
	// Helper to generate random Fp2
	//void generate_random_Fp2(Fp2& v, const Fp& p) {
	//	generate_random_Fp(v[0], p);
	//	generate_random_Fp(v[1], p);
	//}

	void generate_random_Fp2(Fp2& v, InsecureRng& rng) {
		generate_random_Fp(v[0], rng);
		generate_random_Fp(v[1], rng);
	}

	// Helper to compare two Fr&
	bool Fr_equal(const Fr& a, const Fr& b) {
		return std::memcmp(a.l, b.l, sizeof(Fr)) == 0;
	}

	// Helper to compare two Fp
	bool Fp_equal(const Fp& a, const Fp& b) {
		return std::memcmp(a.l, b.l, sizeof(Fp)) == 0;
	}

	// Add this helper function for printing vectors in hexadecimal format
	void print_Fr(const char* label, const Fr& v) {
		std::cout << label << ": ";
		for (int i = NLIMBS(256) - 1; i >= 0; --i) {
			std::cout << std::hex << std::setfill('0') << std::setw(16) << v[i];
			if (i > 0) std::cout << " ";
		}
		std::cout << std::dec << std::endl;
	}

	void print_Fp(const char* label, const Fp& v) {
		std::cout << label << ": ";
		for (int i = NLIMBS(384) - 1; i >= 0; --i) {
			std::cout << std::hex << std::setfill('0') << std::setw(16) << v[i];
			if (i > 0) std::cout << " ";
		}
		std::cout << std::dec << std::endl;
	}

	//void print_UFp(const char* label, const UFp& v) {
	//	printf("%s: ", label);
	//	for (int i = NLIMBS(768) - 1; i >= 0; --i) {
	//		printf("%016llx", (unsigned long long)v[i]);
	//		if (i > 0) printf(" ");
	//	}
	//	printf("\n");
	//}


}

// Test suite for BLST_PORTABLE and BLST_ASM implementations
#if defined(BLST_ASM) && defined(BLST_PORTABLE)

namespace
{
	void platform_mul_mont_sparse_256(void* ret, const void* a, const void* b,
		const void* p, limb_t n0)
	{
#if CATBLST_TEST_X86_64
		mulx_mont_sparse_256(ret, a, b, p, n0);
#else
		mul_mont_sparse_256(ret, a, b, p, n0);
#endif
	}

	void platform_sqr_mont_sparse_256(void* ret, const void* a,
		const void* p, limb_t n0)
	{
#if CATBLST_TEST_X86_64
		sqrx_mont_sparse_256(ret, a, p, n0);
#else
		sqr_mont_sparse_256(ret, a, p, n0);
#endif
	}

	void platform_redc_mont_256(void* ret, const void* a, const void* p,
		limb_t n0)
	{
#if CATBLST_TEST_X86_64
		redcx_mont_256(ret, a, p, n0);
#else
		redc_mont_256(ret, a, p, n0);
#endif
	}

	void platform_mul_mont_384(void* ret, const void* a, const void* b,
		const void* p, limb_t n0)
	{
#if CATBLST_TEST_X86_64
		mulx_mont_384(ret, a, b, p, n0);
#else
		mul_mont_384(ret, a, b, p, n0);
#endif
	}

	void platform_sqr_mont_384(void* ret, const void* a, const void* p,
		limb_t n0)
	{
#if CATBLST_TEST_X86_64
		sqrx_mont_384(ret, a, p, n0);
#else
		sqr_mont_384(ret, a, p, n0);
#endif
	}

	void platform_from_mont_256(void* ret, const void* a, const void* p,
		limb_t n0)
	{
#if CATBLST_TEST_X86_64
		fromx_mont_256(ret, a, p, n0);
#else
		from_mont_256(ret, a, p, n0);
#endif
	}

	void platform_from_mont_384(void* ret, const void* a, const void* p,
		limb_t n0)
	{
#if CATBLST_TEST_X86_64
		fromx_mont_384(ret, a, p, n0);
#else
		from_mont_384(ret, a, p, n0);
#endif
	}

	void platform_sqr_n_mul_mont_383(void* ret, const void* a, size_t count,
		const void* p, limb_t n0, const void* b)
	{
#if CATBLST_TEST_X86_64
		sqrx_n_mul_mont_383(ret, a, count, p, n0, b);
#else
		sqr_n_mul_mont_383(ret, a, count, p, n0, b);
#endif
	}

	void platform_inverse_mod_384(void* ret, const void* inp, const void* mod,
		const void* modx)
	{
#if CATBLST_TEST_X86_64
		ctx_inverse_mod_384(ret, inp, mod, modx);
#else
		ct_inverse_mod_384(ret, inp, mod, modx);
#endif
	}

	void platform_mul_mont_384x(void* ret, const void* a, const void* b,
		const void* p, limb_t n0)
	{
#if CATBLST_TEST_X86_64
		mulx_mont_384x(ret, a, b, p, n0);
#else
		mul_mont_384x(ret, a, b, p, n0);
#endif
	}

	void platform_sqr_mont_384x(void* ret, const void* a, const void* p,
		limb_t n0)
	{
#if CATBLST_TEST_X86_64
		sqrx_mont_384x(ret, a, p, n0);
#else
		sqr_mont_384x(ret, a, p, n0);
#endif
	}

	limb_t platform_sgn0_pty_mont_384(const void* a, const void* p,
		limb_t n0)
	{
#if CATBLST_TEST_X86_64
		return sgn0x_pty_mont_384(a, p, n0);
#else
		return sgn0_pty_mont_384(a, p, n0);
#endif
	}

	limb_t platform_sgn0_pty_mont_384x(const void* a, const void* p,
		limb_t n0)
	{
#if CATBLST_TEST_X86_64
		return sgn0x_pty_mont_384x(a, p, n0);
#else
		return sgn0_pty_mont_384x(a, p, n0);
#endif
	}
}

TEST(InterfaceTests, MulMontSparseV256) {
	Fr a, b, p, ret_asm, ret_portable;
	limb_t n0 = BLST_r0;
	int trials = 10;
	memcpy(p.l, get_BLST_r().l, sizeof(Fr));

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fr(a, p);
		generate_random_Fr(b, p);

		// Compute using both implementations
		platform_mul_mont_sparse_256(ret_asm.l, a.l, b.l, p.l, n0);
		//cc_mul_mont_sparse_256(ret_portable, a, b, p, n0); // Portable implementation
		cc_mul_mont_n<NLIMBS(256)>(ret_portable.l, a.l, b.l, p.l, n0);


		bool equal = Fr_equal(ret_asm, ret_portable);
		EXPECT_TRUE(equal) << "Mismatch in mul_mont_sparse_v256";
	}
}

TEST(InterfaceTests, AddModV384) {
	Fp a, b, p, ret_asm, ret_portable;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));
		generate_random_Fp(a, rng);
		generate_random_Fp(b, rng);

		// Compute using both implementations
		add_mod_384(ret_asm.l, a.l, b.l, p.l); // ASM implementation
		//cc_add_mod_384(ret_portable, a, b, p); // Portable implementation
		cc_add_mod_n<NLIMBS(384)>(ret_portable.l, a.l, b.l, get_BLST_P().l);

		auto eq = Fp_equal(ret_asm, ret_portable);
		if (!eq)
		{
			print_Fp("a", a);
			print_Fp("b", b);
			print_Fp("p", p);
			print_Fp("ret_asm     ", ret_asm);
			print_Fp("ret_portable", ret_portable);
		}

		EXPECT_TRUE(eq) << "Mismatch in add_mod_v384";
	}
}

TEST(InterfaceTests, SqrMontSparseV256) {
	Fr a, p, ret_asm, ret_portable;
	limb_t n0 = BLST_r0;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fr(a);
		memcpy(p.l, get_BLST_r().l, sizeof(Fr));

		// Compute using both implementations
		platform_sqr_mont_sparse_256(ret_asm.l, a.l, p.l, n0);
		//cc_sqr_mont_sparse_256(ret_portable, a, p, n0); // Portable implementation
		cc_mul_mont_n<NLIMBS(256)>(ret_portable.l, a.l, a.l, p.l, n0);

		EXPECT_TRUE(Fr_equal(ret_asm, ret_portable)) << "Mismatch in sqr_mont_sparse_v256";
	}
}

TEST(InterfaceTests, RedcMontV256) {
	UFr a;
	Fr p, ret_asm, ret_portable;
	limb_t n0 = BLST_r0;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		Fr r;
		generate_random_Fr(r); // Fill lower 256 bits
        memcpy(a.data, r.l, sizeof(Fr));
		generate_random_Fr(r); // Fill upper 256 bits
        memcpy(a.data + NLIMBS(256), r.l, sizeof(Fr));

		memcpy(p.l, get_BLST_r().l, sizeof(Fr));
		// Compute using both implementations
		platform_redc_mont_256(ret_asm.l, a, p.l, n0);
		//cc_redc_mont_256(ret_portable, a, p, n0); // Portable implementation
		cc_redc_mont_n<NLIMBS(256)>(ret_portable.l, a.data, p.l, n0);

		EXPECT_TRUE(Fr_equal(ret_asm, ret_portable)) << "Mismatch in redc_mont_v256";
	}
}

TEST(InterfaceTests, MulMontV384) {
	Fp a, b, p, ret_asm, ret_portable;
	limb_t n0 = BLST_p0;
	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp(a, rng);
		generate_random_Fp(b, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));
		// Compute using both implementations
		platform_mul_mont_384(ret_asm.l, a.l, b.l, p.l, n0);
		//cc_mul_mont_384(ret_portable, a, b, p, n0); // Portable implementation
		cc_mul_mont_n<NLIMBS(384)>(ret_portable.l, a.l, b.l, get_BLST_P().l, BLST_p0);

		EXPECT_TRUE(Fp_equal(ret_asm, ret_portable)) << "Mismatch in mul_mont_v384";
	}
}

TEST(InterfaceTests, SqrMontV384) {
	Fp a, p, ret_asm, ret_portable;
	limb_t n0 = BLST_p0;
	int trials = 10;
	InsecureRng rng;
	for (int t = 0; t < trials; t++)
	{

		generate_random_Fp(a, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		// Compute using both implementations
		platform_sqr_mont_384(ret_asm.l, a.l, p.l, n0);
		//cc_sqr_mont_384(ret_portable, a, p, n0); // Portable implementation
		cc_mul_mont_n<NLIMBS(384)>(ret_portable.l, a.l, a.l, get_BLST_P().l, BLST_p0);

		EXPECT_TRUE(Fp_equal(ret_asm, ret_portable)) << "Mismatch in sqr_mont_v384";
	}
}

TEST(InterfaceTests, AddModV256) {
	Fr a, b, p, ret_asm, ret_portable;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p.l, get_BLST_r().l, sizeof(Fr));
		generate_random_Fr(a, p);
		generate_random_Fr(b, p);
		// Compute using both implementations
		add_mod_256(ret_asm.l, a.l, b.l, p.l); // ASM implementation
		//cc_add_mod_256(ret_portable, a, b, p); // Portable implementation
		cc_add_mod_n<NLIMBS(256)>(ret_portable.l, a.l, b.l, p.l);

		EXPECT_TRUE(Fr_equal(ret_asm, ret_portable)) << "Mismatch in add_mod_v256";
	}
}

TEST(InterfaceTests, SubModV256) {
	Fr a, b, p, ret_asm, ret_portable;
	int trials = 10;
	for (int t = 0; t < trials; t++)
	{
		memcpy(p.l, get_BLST_r().l, sizeof(Fr));
		generate_random_Fr(a, p);
		generate_random_Fr(b, p);

		// Compute using both implementations
		sub_mod_256(ret_asm.l, a.l, b.l, p.l); // ASM implementation
		//cc_sub_mod_256(ret_portable, a, b, p); // Portable implementation
		cc_sub_mod_n<NLIMBS(256)>(ret_portable.l, a.l, b.l, p.l);


		EXPECT_TRUE(Fr_equal(ret_asm, ret_portable)) << "Mismatch in sub_mod_v256";
	}
}

TEST(InterfaceTests, MulBy3ModV256) {
	Fr a, p, ret_asm, ret_portable;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p.l, get_BLST_r().l, sizeof(Fr));
		generate_random_Fr(a, p);

		// Compute using both implementations
		mul_by_3_mod_256(ret_asm.l, a.l, p.l); // ASM implementation
		//cc_mul_by_3_mod_256(ret_portable, a, p); // Portable implementation
		cc_mul_by_3_mod_n<NLIMBS(256)>(ret_portable.l, a.l, p.l);


		EXPECT_TRUE(Fr_equal(ret_asm, ret_portable)) << "Mismatch in mul_by_3_mod_v256";
	}
}

TEST(InterfaceTests, CNegModV256) {
	Fr a, p, ret_asm, ret_portable;
	bool_t flag = true;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p.l, get_BLST_r().l, sizeof(Fr));
		generate_random_Fr(a, p);
		// Compute using both implementations
		cneg_mod_256(ret_asm.l, a.l, flag, p.l); // ASM implementation
		//cc_cneg_mod_256(ret_portable, a, flag, p); // Portable implementation
		cc_cneg_mod_n<NLIMBS(256)>(ret_portable.l, a.l, flag, p.l);

		EXPECT_TRUE(Fr_equal(ret_asm, ret_portable)) << "Mismatch in cneg_mod_v256";
	}
}

// Test for lshift_mod_v256
TEST(InterfaceTests, LShiftModV256) {
	Fr a, p, ret_asm, ret_portable;
	size_t count = 3;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p.l, get_BLST_r().l, sizeof(Fr));
		generate_random_Fr(a, p);

		lshift_mod_256(ret_asm.l, a.l, count, p.l);
		//cc_lshift_mod_256(ret_portable, a, count, p);
		cc_lshift_mod_n<NLIMBS(256)>(ret_portable.l, a.l, count, p.l);

		auto eq = Fr_equal(ret_asm, ret_portable);
		if (!eq)
		{
			print_Fr("a", a);
			print_Fr("p", p);
			print_Fr("ret_asm", ret_asm);
			print_Fr("ret_portable", ret_portable);
		}
		EXPECT_TRUE(eq) << "Mismatch in lshift_mod_v256";
	}
}

// Test for check_mod_256
TEST(InterfaceTests, CheckMod256) {
	pow256 a;
	Fr p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		memcpy(p.l, get_BLST_r().l, sizeof(Fr));
		generate_random_Fr(reinterpret_cast<Fr&>(a));

		// Test both implementations
		limb_t result_asm = check_mod_256(a, p.l); // ASM implementation
		//limb_t result_portable = cc_check_mod_256(a, p); // Portable implementation
		limb_t result_portable = cc_check_mod_n<NLIMBS(256)>(a, p.l);


		// Compare results
		EXPECT_EQ(result_asm, result_portable) << "Mismatch in check_mod_256";
	}
}

// Test for add_n_check_mod_256
TEST(InterfaceTests, AddNCheckMod256) {
	pow256 a, b, ret_asm, ret_portable;
	Fr p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		memcpy(p.l, get_BLST_r().l, sizeof(Fr));
		generate_random_Fr(reinterpret_cast<Fr&>(a), p);
		generate_random_Fr(reinterpret_cast<Fr&>(b), p);

		// Make copies for both implementations
		memcpy(ret_asm, a, sizeof(pow256));
		memcpy(ret_portable, a, sizeof(pow256));

		// Test both implementations
		limb_t result_asm = add_n_check_mod_256(ret_asm, a, b, p.l); // ASM implementation
		//limb_t result_portable = cc_add_n_check_mod_256(ret_portable, a, b, p); // Portable implementation
		limb_t result_portable = cc_add_n_check_mod_n<NLIMBS(256)>(ret_portable, a, b, p.l);

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
	Fr p;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		generate_random_Fr(reinterpret_cast<Fr&>(a));
		generate_random_Fr(reinterpret_cast<Fr&>(b));
		memcpy(p.l, get_BLST_r().l, sizeof(Fr));

		// Make copies for both implementations
		memcpy(ret_asm, a, sizeof(pow256));
		memcpy(ret_portable, a, sizeof(pow256));

		// Test both implementations
		limb_t result_asm = sub_n_check_mod_256(ret_asm, a, b, p.l); // ASM implementation
		//limb_t result_portable = cc_sub_n_check_mod_256(ret_portable, a, b, p); // Portable implementation
		limb_t result_portable = cc_sub_n_check_mod_n<NLIMBS(256)>(ret_portable, a, b, p.l);

		// Compare return values
		EXPECT_EQ(result_asm, result_portable) << "Mismatch in sub_n_check_mod_256 return value";

		// Compare modified results
		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(pow256)) == 0)
			<< "Mismatch in sub_n_check_mod_256 output";
	}
}

// Test for rshift_mod_v256
TEST(InterfaceTests, RShiftModV256) {
	Fr a, p, ret_asm, ret_portable;
	size_t count = 3;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fr(a);
		memcpy(p.l, get_BLST_r().l, sizeof(Fr));

		rshift_mod_256(ret_asm.l, a.l, count, p.l);
		//cc_rshift_mod_256(ret_portable, a, count, p);
		cc_rshift_mod_n<NLIMBS(256)>(ret_portable.l, a.l, count, p.l);


		EXPECT_TRUE(Fr_equal(ret_asm, ret_portable)) << "Mismatch in rshift_mod_v256";
	}
}

// Test for from_mont_v256
TEST(InterfaceTests, FromMontV256) {
	Fr a, p, ret_asm, ret_portable;
	limb_t n0 = BLST_r0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fr(a);
		memcpy(p.l, get_BLST_r().l, sizeof(Fr));

		// Compute using both implementations
		platform_from_mont_256(ret_asm.l, a.l, p.l, n0);
		//cc_from_mont_256(ret_portable, a, p, n0); // Portable implementation
		cc_from_mont_n<NLIMBS(256)>(ret_portable.l, a.l, p.l, n0);


		EXPECT_TRUE(Fr_equal(ret_asm, ret_portable)) << "Mismatch in from_mont_v256";
	}
}

// Test for from_mont_v384
TEST(InterfaceTests, FromMontV384) {
	Fp a, p, ret_asm, ret_portable;
	limb_t n0 = BLST_p0;

	int trials = 10;

	InsecureRng rng;
	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp(a, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		// Compute using both implementations
		platform_from_mont_384(ret_asm.l, a.l, p.l, n0);
		//cc_from_mont_384(ret_portable, a, p, n0); // Portable implementation
		cc_from_mont_n<NLIMBS(384)>(ret_portable.l, a.l, get_BLST_P().l, BLST_p0);


		EXPECT_TRUE(Fp_equal(ret_asm, ret_portable)) << "Mismatch in from_mont_v384";
	}
}
//// Test for sqr_n_mul_mont_v384
//TEST(InterfaceTests, SqrNMulMontV384) {
//	Fp a, b, p, ret_asm, ret_portable;
//	limb_t n0 = BLST_p0;
//	size_t count = 1; // Number of squaring operations
//
//	int trials = 20;
//
//	for (int t = 0; t < trials; t++)
//	{
//		memcpy(p, get_BLST_P(), sizeof(Fp));
//		generate_random_Fp(a, p);
//		generate_random_Fp(b, p);
//		//memcpy(p, get_BLST_P(), sizeof(Fp));
//
//		// Compute using both implementations
//		sqrx_n_mul_mont_384(ret_asm, a, count, p, n0, b); // ASM implementation
//		cc_sqr_n_mul_mont_384(ret_portable, a, count, p, n0, b); // Portable implementation
//
//		auto eq = Fp_equal(ret_asm, ret_portable);
//		if (eq == false)
//		{
//			// If they are not equal, print the vectors for debugging
//			std::cout << "Mismatch in sqr_n_mul_mont_v384" << std::endl;
//			print_Fp("a", a);
//			print_Fp("b", b);
//			print_Fp("p", p);
//
//
//			print_Fp("ret_asm     ", ret_asm);
//			print_Fp("ret_portable", ret_portable);
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
//static inline uint64_t mont_n0_from_p(Fp p_le) {
//	return mont_n0_from_p0(p_le[0]); // little-endian limb_ts: p[0] is LSW
//}

// Test for sqr_n_mul_mont_v383
TEST(InterfaceTests, SqrNMulMontV383) {
	Fp a, b, aM, bM, p, ret_asm, ret_portable, out_asm, out_c;
	memcpy(p.l, get_BLST_P().l, sizeof(Fp));

	limb_t n0 = BLST_p0;
	//auto nn0 = mont_n0_from_p(get_BLST_P());
	//assert(n0 == nn0);
	size_t count = 5; // Number of squaring operations

	int trials = 40;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp(a, rng);
		generate_random_Fp(b, rng);

		EXPECT_TRUE(cc_check_mod_n<NLIMBS(384)>((byte*)b.l, p.l) == 1);
		//to_fp(a, a);
		//to_fp(b, b);
		// ---- to Montgomery correctly: aM = a * R  and  bM = b * R
		// via MontMul with RR
		cc_mul_mont_n<NLIMBS(384)>(aM.l, a.l, get_BLST_RR().l, p.l, n0);
		cc_mul_mont_n<NLIMBS(384)>(bM.l, b.l, get_BLST_RR().l, p.l, n0);


		// Compute using both implementations
		platform_sqr_n_mul_mont_383(ret_asm.l, a.l, count, p.l, n0,
			b.l);
		cc_sqr_n_mul_mont_384(ret_portable.l, a.l, count, p.l, n0, b.l); // Portable implementation

		// ---- back from Montgomery for comparison
		Fp::from_mont(out_asm, ret_asm);
		Fp::from_mont(out_c, ret_portable);


		EXPECT_TRUE(Fp_equal(out_asm, out_c)) << "Mismatch in sqr_n_mul_mont_v383";
	}
}

TEST(InterfaceTests, SubModV384) {
	Fp a, b, p, ret_asm, ret_portable;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));
		generate_random_Fp(a, rng);
		generate_random_Fp(b, rng);

		sub_mod_384(ret_asm.l, a.l, b.l, p.l);
		cc_sub_mod_n<NLIMBS(384)>(ret_portable.l, a.l, b.l, p.l);

		EXPECT_TRUE(Fp_equal(ret_asm, ret_portable)) << "Mismatch in sub_mod_v384";
	}
}

TEST(InterfaceTests, MulBy8ModV384) {
	Fp a, p, ret_asm, ret_portable;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));
		generate_random_Fp(a, rng);

		mul_by_8_mod_384(ret_asm.l, a.l, p.l);
		cc_lshift_mod_n<NLIMBS(384)>(ret_portable.l, a.l, 3, p.l);

		EXPECT_TRUE(Fp_equal(ret_asm, ret_portable)) << "Mismatch in mul_by_8_mod_v384";
	}
}

TEST(InterfaceTests, MulBy3ModV384) {
	Fp a, p, ret_asm, ret_portable;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));
		generate_random_Fp(a, rng);

		mul_by_3_mod_384(ret_asm.l, a.l, p.l);
		//cc_mul_by_3_mod_384(ret_portable, a, p);
		cc_mul_by_3_mod_n<NLIMBS(384)>(ret_portable.l, a.l, get_BLST_P().l);


		EXPECT_TRUE(Fp_equal(ret_asm, ret_portable)) << "Mismatch in mul_by_3_mod_v384";
	}
}

TEST(InterfaceTests, CNegModV384) {
	Fp a, p, ret_asm, ret_portable;
	bool_t flag = true;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));
		generate_random_Fp(a, rng);

		cneg_mod_384(ret_asm.l, a.l, flag, p.l);
		//cc_cneg_mod_384(ret_portable, a, flag, p);
		cc_cneg_mod_n<NLIMBS(384)>(ret_portable.l, a.l, flag, p.l);


		EXPECT_TRUE(Fp_equal(ret_asm, ret_portable)) << "Mismatch in cneg_mod_v384";
	}
}


TEST(InterfaceTests, RShiftMod384) {
	Fp a, p, ret_asm, ret_portable;
	size_t count = 3;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));
		generate_random_Fp(a, rng);

		rshift_mod_384(ret_asm.l, a.l, count, p.l); // ASM implementation
		//cc_rshift_mod_384(ret_portable, a, count, p); // Portable implementation
		cc_rshift_mod_n<NLIMBS(384)>(ret_portable.l, a.l, count, get_BLST_P().l);


		EXPECT_TRUE(Fp_equal(ret_asm, ret_portable)) << "Mismatch in rshift_mod_384";
	}
}

// Test for redcx_mont_384
#if CATBLST_TEST_X86_64
TEST(InterfaceTests, RedcxMont384) {
	UFp a;
	Fp p, ret_asm, ret_portable;
	limb_t n0 = BLST_p0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		for (int i = 0; i < NLIMBS(768); ++i) {
			a[i] = rand();
		}
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		// Compute using both implementations
		redcx_mont_384(ret_asm.l, a, p.l, n0); // ASM implementation
		//cc_redc_mont_384(ret_portable, a, p, n0); // Portable implementation
		cc_redc_mont_n<NLIMBS(384)>(ret_portable.l, a.data, get_BLST_P().l, BLST_p0);


		EXPECT_TRUE(Fp_equal(ret_asm, ret_portable)) << "Mismatch in redcx_mont_384";
	}
}
#endif

// Test for redc_mont_384
TEST(InterfaceTests, RedcMont384) {
	UFp a;
	Fp p, ret_asm, ret_portable;
	limb_t n0 = BLST_p0;

	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		for (int i = 0; i < NLIMBS(768); ++i) {
			a[i] = rand();
		}
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		// Compute using both implementations
		redc_mont_384(ret_asm.l, a.data, p.l, n0); // ASM implementation
		cc_redc_mont_n<NLIMBS(384)>(ret_portable.l, a.data, p.l, n0);

		EXPECT_TRUE(Fp_equal(ret_asm, ret_portable)) << "Mismatch in redc_mont_384";
	}
}


TEST(InterfaceTests, DivBy2Mod384) {
	Fp a, p, ret_asm, ret_portable;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp(a, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		div_by_2_mod_384(ret_asm.l, a.l, p.l); // ASM implementation
		//cc_div_by_2_mod_384(ret_portable, a, p); // Portable implementation
		cc_rshift_mod_n<NLIMBS(384)>(ret_portable.l, a.l, 1, get_BLST_P().l);


		EXPECT_TRUE(Fp_equal(ret_asm, ret_portable)) << "Mismatch in div_by_2_mod_384";
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
	// find highest nonzero limb_t
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

const Fp modx = {    /* left-aligned value of the modulus */
	TO_LIMB_T(0xcff7fffffffd5558), TO_LIMB_T(0xf55ffff58a9ffffd),
	TO_LIMB_T(0x39869507b587b120), TO_LIMB_T(0x23ba5c279c2895fb),
	TO_LIMB_T(0x58dd3db21a5d66bb), TO_LIMB_T(0xd0088f51cbff34d2)
};

TEST(InterfaceTests, CtxInverseMod384) {


	Fp inp;
	UFp ret_asm, ret_portable;

	{
		// Verify modx is correctly left-aligned
		Fp modx_check;
		left_align_mod(modx_check.l, get_BLST_P().l, NLIMBS(384));
		assert(std::memcmp(modx.l, modx_check.l, sizeof(Fp)) == 0);
	}

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp(inp, rng);

		platform_inverse_mod_384(ret_asm, inp.l, get_BLST_P().l,
			modx.l);
		//cc_ct_inverse_mod_384(ret_portable, inp, get_BLST_P(), modx); // Portable implementation
		cc_ct_inverse_mod_n<NLIMBS(384)>(ret_portable.data, inp.l, get_BLST_P().l, modx.l);


		const size_t n = NLIMBS(384);

		// handy views
		auto* asm_hi = ret_asm + n;
		auto* por_hi = ret_portable + n;

		// 0) High-half equality (raw 2n output)
		int eq0 = std::memcmp(asm_hi, por_hi, sizeof(Fp));

		// 1) REDC stage equality
		union { UFp T; Fp r[2]; } A, B;
		std::memcpy(A.T, ret_asm, sizeof(UFp));
		std::memcpy(B.T, ret_portable, sizeof(UFp));
		Fp redc_asm, redc_por;
		Fp::redc(redc_asm, A.T);
		Fp::redc(redc_por, B.T);
		int eq1 = std::memcmp(redc_asm.l, redc_por.l, sizeof(Fp));

		// 2) Convert-to-Montgomery equality
		Fp inv_asm, inv_por;
		Fp::mul(inv_asm, redc_asm, get_BLST_RR());
		Fp::mul(inv_por, redc_por, get_BLST_RR());
		int eq2 = std::memcmp(inv_asm.l, inv_por.l, sizeof(Fp));

		// 3) Functional correctness (both sides!)
		Fp oneA, oneP;
		if (!vec_is_zero(inp.l, sizeof(Fp))) {
			Fp::mul(oneA, inv_asm, inp);
			Fp::mul(oneP, inv_por, inp);
			EXPECT_TRUE(vec_is_equal(oneA.l, get_BLST_Rx().p.l, sizeof(Fp))) << "ASM inverse incorrect";
			EXPECT_TRUE(vec_is_equal(oneP.l, get_BLST_Rx().p.l, sizeof(Fp))) << "Portable inverse incorrect";
		}

		// Assert stage-by-stage (these help pinpoint the break)
		EXPECT_EQ(eq0, 0) << "Mismatch in raw high halves (ret+n)";
		EXPECT_EQ(eq1, 0) << "Mismatch after REDC";
		EXPECT_EQ(eq2, 0) << "Mismatch after convert-to-Mont";
	}
}

TEST(InterfaceTests, CtInverseMod256) {
	Fr modx;
	left_align_mod(modx.l, get_BLST_r().l, NLIMBS(256));

	auto check_inverse = [&modx](const Fr& inp) {
		UFr ret_asm, ret_portable;
		ct_inverse_mod_256(ret_asm, inp.l, get_BLST_r().l, modx.l); // ASM implementation
		cc_ct_inverse_mod_256(ret_portable, inp, get_BLST_r(), modx); // Portable implementation

		EXPECT_TRUE(std::memcmp(ret_asm, ret_portable, sizeof(UFr)) == 0) << "Mismatch in ct_inverse_mod_256";
		if (vec_is_zero(inp.l, sizeof(inp)))
			return;

		Fr inv_asm, inv_portable, product;
		Fr::redc(inv_asm, ret_asm);
		Fr::redc(inv_portable, ret_portable);
		Fr::to_mont(inv_asm, inv_asm);
		Fr::to_mont(inv_portable, inv_portable);

		Fr::mul(product, inp, inv_asm);
		EXPECT_TRUE(Fr_equal(product, Fr::one())) << "ASM inverse is not multiplicative";
		Fr::mul(product, inp, inv_portable);
		EXPECT_TRUE(Fr_equal(product, Fr::one())) << "Portable inverse is not multiplicative";
	};

	Fr zero = {};
	check_inverse(zero);
	check_inverse(Fr::one());

	Fr modulus_minus_one = get_BLST_r();
	modulus_minus_one[0] -= 1;
	check_inverse(modulus_minus_one);

	std::mt19937_64 rng(0x4354424c5354494eULL);
	for (int t = 0; t < 10; t++)
	{
		Fr inp;
		for (size_t i = 0; i < NLIMBS(256) - 1; i++)
			inp[i] = rng();
		inp[NLIMBS(256) - 1] = rng() % get_BLST_r()[NLIMBS(256) - 1];
		check_inverse(inp);
	}
}

TEST(InterfaceTests, CtIsSquareMod384) {
	Fp inp;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp(inp, rng);

		bool_t is_square_asm = ct_is_square_mod_384(inp.l, get_BLST_P().l); // ASM implementation
		//bool_t is_square_portable = cc_ct_is_square_mod_384(inp, get_BLST_P()); // Portable implementation
		bool_t is_square_portable = cc_ct_is_sqr_mod_n<NLIMBS(384)>(inp.l, get_BLST_P().l);

		//printf("passed? %d\n", is_square_asm == is_square_portable);
		//print_Fp("inp", inp);
		//printf("is_square_asm     = %lld\n", is_square_asm);
		//printf("is_square_portable= %lld\n", is_square_portable);

		EXPECT_EQ(is_square_asm, is_square_portable) << "Mismatch in ct_is_square_mod_384";
	}
}

TEST(InterfaceTests, MulxMont384x) {
	Fp2 a, b, ret_asm, ret_portable;
	Fp p;
	limb_t n0 = BLST_p0;

	InsecureRng rng;
	int trials = 10;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp2(a, rng);
		generate_random_Fp2(b, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		platform_mul_mont_384x(&ret_asm, &a, &b, p.l, n0);
		cc_mul_mont_384x(ret_portable, a, b); // Portable implementation

		EXPECT_TRUE(std::memcmp(&ret_asm, &ret_portable, sizeof(Fp2)) == 0) << "Mismatch in mulx_mont_384x";
	}
}

TEST(InterfaceTests, SqrxMont384x) {
	Fp2 a, ret_asm, ret_portable;
	Fp p;
	limb_t n0 = BLST_p0;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp2(a, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		platform_sqr_mont_384x(&ret_asm, &a, p.l, n0);
		cc_sqr_mont_384x(ret_portable, a); // Portable implementation

		EXPECT_TRUE(std::memcmp(&ret_asm, &ret_portable, sizeof(Fp2)) == 0) << "Mismatch in sqrx_mont_384x";
	}
}

TEST(InterfaceTests, Wide382ProvidersMatchPortableFieldArithmetic) {
	InsecureRng rng;

	for (int t = 0; t < 10; ++t) {
		Fp2 a, b, expected_mul, expected_sqr;
		UFp2 wide_mul, wide_sqr;
		generate_random_Fp2(a, rng);
		generate_random_Fp2(b, rng);

#if CATBLST_TEST_X86_64
		UFp2 wide_mul_q, wide_sqr_q;
		mulx_382x(&wide_mul, &a, &b, &get_BLST_P());
		sqrx_382x(&wide_sqr, &a, &get_BLST_P());
		mul_382x(&wide_mul_q, &a, &b, &get_BLST_P());
		sqr_382x(&wide_sqr_q, &a, &get_BLST_P());
		EXPECT_EQ(std::memcmp(&wide_mul, &wide_mul_q, sizeof(wide_mul)), 0);
		EXPECT_EQ(std::memcmp(&wide_sqr, &wide_sqr_q, sizeof(wide_sqr)), 0);
#else
		mul_382x(&wide_mul, &a, &b, &get_BLST_P());
		sqr_382x(&wide_sqr, &a, &get_BLST_P());
#endif

		Fp2 reduced_mul, reduced_sqr;
		Fp::redc(reduced_mul[0], wide_mul[0]);
		Fp::redc(reduced_mul[1], wide_mul[1]);
		Fp::redc(reduced_sqr[0], wide_sqr[0]);
		Fp::redc(reduced_sqr[1], wide_sqr[1]);
		cc_mul_mont_384x(expected_mul, a, b);
		cc_sqr_mont_384x(expected_sqr, a);

		EXPECT_EQ(std::memcmp(&reduced_mul, &expected_mul, sizeof(Fp2)), 0);
		EXPECT_EQ(std::memcmp(&reduced_sqr, &expected_sqr, sizeof(Fp2)), 0);
	}
}

TEST(InterfaceTests, AddMod384x) {
	Fp2 a, b, ret_asm, ret_portable;
	Fp p;

	int trials = 10;

	InsecureRng rng;
	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp2(a, rng);
		generate_random_Fp2(b, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		add_mod_384x(&ret_asm, &a, &b, p.l); // ASM implementation
		cc_add_mod_384x(ret_portable, a, b); // Portable implementation

		EXPECT_TRUE(std::memcmp(&ret_asm, &ret_portable, sizeof(Fp2)) == 0) << "Mismatch in add_mod_384x";
	}
}

TEST(InterfaceTests, SubMod384x) {
	Fp2 a, b, ret_asm, ret_portable;
	Fp p;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp2(a, rng);
		generate_random_Fp2(b, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		sub_mod_384x(&ret_asm, &a, &b, p.l); // ASM implementation
		cc_sub_mod_384x(ret_portable, a, b); // Portable implementation

		EXPECT_TRUE(std::memcmp(&ret_asm, &ret_portable, sizeof(Fp2)) == 0) << "Mismatch in sub_mod_384x";
	}
}

// Test for mul_by_8_mod_384x
TEST(InterfaceTests, MulBy8Mod384x) {
	Fp2 a, ret_asm, ret_portable;
	Fp p;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		generate_random_Fp(a[0], rng);
		generate_random_Fp(a[1], rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		// Test both implementations
		mul_by_8_mod_384x(&ret_asm, &a, p.l); // ASM implementation
		cc_mul_by_8_mod_384x(ret_portable, a); // Portable implementation

		// Compare results
		EXPECT_TRUE(std::memcmp(&ret_asm, &ret_portable, sizeof(Fp2)) == 0) << "Mismatch in mul_by_8_mod_384x";
	}
}

// Test for mul_by_3_mod_384x
TEST(InterfaceTests, MulBy3Mod384x) {
	Fp2 a, ret_asm, ret_portable;
	Fp p;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		generate_random_Fp(a[0], rng);
		generate_random_Fp(a[1], rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		// Test both implementations
		mul_by_3_mod_384x(&ret_asm, &a, p.l); // ASM implementation
		cc_mul_by_3_mod_384x(ret_portable, a); // Portable implementation

		// Compare results
		EXPECT_TRUE(std::memcmp(&ret_asm, &ret_portable, sizeof(Fp2)) == 0) << "Mismatch in mul_by_3_mod_384x";
	}
}

// Test for mul_by_1_plus_i_mod_384x
TEST(InterfaceTests, MulBy1PlusIMod384x) {
	Fp2 a, ret_asm, ret_portable;
	Fp p;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		// Generate random inputs
		generate_random_Fp(a[0], rng);
		generate_random_Fp(a[1], rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		// Test both implementations
		mul_by_1_plus_i_mod_384x(&ret_asm, &a, p.l); // ASM implementation
		cc_mul_by_1_plus_i_mod_384x(ret_portable, a); // Portable implementation

		// Compare results
		EXPECT_TRUE(std::memcmp(&ret_asm, &ret_portable, sizeof(Fp2)) == 0) << "Mismatch in mul_by_1_plus_i_mod_384x";
	}
}

// Test for sgn0x_pty_mont_384
TEST(InterfaceTests, Sgn0PtyMont384) {
	Fp a, p;
	limb_t n0 = BLST_p0;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp(a, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		// Compute using both implementations
		limb_t result_asm = platform_sgn0_pty_mont_384(a.l, p.l, n0);
		//limb_t result_portable = cc_sgn0_pty_mont_384(a, p, n0); // Portable implementation

		limb_t result_portable;
		{
			Fp tmp;
			cc_from_mont_n<NLIMBS(384)>(tmp.l, a.l, p.l, n0);
			result_portable = cc_sgn0_pty_mod_n<NLIMBS(384)>(tmp.l, p.l);
		}

		EXPECT_EQ(result_asm, result_portable) << "Mismatch in sgn0x_pty_mont_384";
	}
}

// Test for sgn0x_pty_mont_384x (complex number version)
TEST(InterfaceTests, Sgn0PtyMont384x) {
	Fp2 a;
	Fp p;
	limb_t n0 = BLST_p0;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp2(a, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		// Compute using both implementations
		limb_t result_asm = platform_sgn0_pty_mont_384x(&a, p.l, n0);
		limb_t result_portable = cc_sgn0_pty_mont_384x(a); // Portable implementation

		EXPECT_EQ(result_asm, result_portable) << "Mismatch in sgn0x_pty_mont_384x";
	}
}

// Test for sgn0_pty_mod_384
TEST(InterfaceTests, Sgn0PtyMod384) {
	Fp a, p;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp(a, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		// Compute using both implementations
		limb_t result_asm = sgn0_pty_mod_384(a.l, p.l); // ASM implementation
		limb_t result_portable = cc_sgn0_pty_mod_n<NLIMBS(384)>(a.l, p.l); // Portable implementation

		EXPECT_EQ(result_asm, result_portable) << "Mismatch in sgn0_pty_mod_384";
	}
}

// Test for sgn0_pty_mod_384x (complex number version)
TEST(InterfaceTests, Sgn0PtyMod384x) {
	Fp2 a;
	Fp p;

	int trials = 10;
	InsecureRng rng;

	for (int t = 0; t < trials; t++)
	{
		generate_random_Fp2(a, rng);
		memcpy(p.l, get_BLST_P().l, sizeof(Fp));

		// Compute using both implementations
		limb_t result_asm = sgn0_pty_mod_384x(&a, p.l); // ASM implementation
		limb_t result_portable = cc_sgn0_pty_mod_384x(a); // Portable implementation

		EXPECT_EQ(result_asm, result_portable) << "Mismatch in sgn0_pty_mod_384x";
	}
}

// Test for SHA256 functions
TEST(InterfaceTests, SHA256Functions) {
	// Test blst_sha256_block_data_order
	const size_t blocks = 2;
	const size_t block_size = 64; // SHA256 block size in bytes
	unsigned int v_asm[8], v_portable[8], v_selected[8]; // SHA256 state
	unsigned char input[blocks * block_size];

	// Initialize state and input data
	for (size_t i = 0; i < 8; ++i) {
		v_asm[i] = v_portable[i] = v_selected[i] =
			0x12345678 + static_cast<unsigned int>(i);
	}

	for (size_t i = 0; i < blocks * block_size; ++i)
		input[i] = static_cast<unsigned char>((i * 29 + 17) & 0xff);

	// Process blocks with both implementations and the selected wrapper.
	blst_sha256_block_data_order(v_asm, input, blocks); // ASM implementation
	cc_blst_sha256_block_data_order(v_portable, input, blocks); // Portable implementation
	blst_sha256_block_data_order_v(v_selected, input, blocks);

	// Compare results
	EXPECT_TRUE(std::memcmp(v_asm, v_portable, sizeof(v_asm)) == 0) << "Mismatch in blst_sha256_block_data_order";
	EXPECT_TRUE(std::memcmp(v_selected, v_portable, sizeof(v_selected)) == 0) << "Mismatch in selected blst_sha256_block_data_order";

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

TEST(InterfaceTests, QuotRemAssemblyMatchesPortable) {
	const limb_t divisor_128[2] = { 5, 1 };
	const limb_t dividends_128[][3] = {
		{ 38, 7, 0 }, // exact estimate: 7 * divisor + 3
		{ 33, 6, 0 }  // estimate 7 is one too high
	};
	const limb_t estimates_128[] = { 7, 7 };

	for (size_t i = 0; i < 2; ++i) {
		limb_t asm_result[3], portable_result[3];
		std::memcpy(asm_result, dividends_128[i], sizeof(asm_result));
		std::memcpy(portable_result, dividends_128[i], sizeof(portable_result));
		const limb_t asm_q = quot_rem_128(asm_result, divisor_128,
			estimates_128[i]);
		const limb_t portable_q = cc_quot_rem_n<NLIMBS(128)>(portable_result,
			divisor_128, estimates_128[i]);

		EXPECT_EQ(asm_q, portable_q);
		EXPECT_EQ(std::memcmp(asm_result, portable_result, sizeof(asm_result)), 0);
	}

	const limb_t divisor_64[1] = { 5 };
	limb_t asm_result_64[2] = { 38, 0 };
	limb_t portable_result_64[2] = { 38, 0 };
	const limb_t asm_q_64 = quot_rem_64(asm_result_64, divisor_64, 7);
	const limb_t portable_q_64 = cc_quot_rem_n<NLIMBS(64)>(portable_result_64,
		divisor_64, 7);

	EXPECT_EQ(asm_q_64, portable_q_64);
	EXPECT_EQ(std::memcmp(asm_result_64, portable_result_64,
		sizeof(asm_result_64)), 0);
}

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
// Sample R in [0, D] to satisfy the stated precondition "dividend top limb_ts not larger than divisor"
static inline std::pair<limb_t, limb_t> random_R_leq_D(std::mt19937_64& rng, limb_t d_hi, limb_t d_lo) {
	std::uniform_int_distribution<limb_t> dist;
	while (true) {
		limb_t r_hi = dist(rng);
		limb_t r_lo = dist(rng);
		if (r_hi < d_hi || (r_hi == d_hi && r_lo <= d_lo)) return { r_hi, r_lo };
	}
}
// Thin wrappers so we don't duplicate the div_top packing or change call order by mistake.
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

#undef CATBLST_TEST_X86_64
