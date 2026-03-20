#include <gtest/gtest.h>


#include <cstdint>
#include <vector>
#include <algorithm>

#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "InsecureRng.h"
#include "blst/cuda/fields/parallel_prefix.h"
#include "blst/backend/cuda/CSA_cuda.h"
#include "blst/backend/cuda/mul_mont_cuda.h"


using namespace blst;


namespace
{

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

		BLST_ASSERT(cc_check_mod_n<NLIMBS(384)>((byte*)v.data, get_BLST_P().data));
	}

	//void set_Fp(Fp& v, uint64_t x) {
	//	for (int i = 0; i < NLIMBS(384); ++i)
	//		v[i] = 0;

	//	memcpy(v.data, &x, sizeof(x));
	//}

	// Helper: Print Fp for debugging
	void print_Fp(const Fp& v) {
		for (int i = NLIMBS(384) - 1; i >= 0; --i)
			printf("%016llx%s", (unsigned long long)v[i], i ? " " : "\n");
	}
}
//
//TEST(ParallelPrefix, AddPart) {
//	InsecureRng rng;
//
//	for (u64 t = 0; t < 100; ++t)
//	{
//
//		Fp a, b;
//		Fp ab, abp;
//		u32 abCarry, abpCarry;
//		constexpr auto n = Fp::num_limb();
//
//		rand_Fp(a, rng);
//		rand_Fp(b, rng);
//		Fp np = Fp::modulus();
//
//		for (u64 i = 0; i < n; ++i)
//			np[i] = ~np[i] + (i == 0 ? 1 : 0); // np = -p mod 2^(64*NLIMBS)
//
//
//		packed_dual_prefix_add_shiftedC<6>(
//			ab.data, abp.data,
//			abCarry, abpCarry,
//			a.data, b.data,
//			Fp::modulus().data, np.data);
//
//		// verify ab = a + b
//		//        abp = a + b - p
//
//		//Fp expected_ab, expected_abp;
//		limb_t mask, carry, borrow;
//		//BLST_VLA_NEW(limb_t, tmp, n);
//		limb_t expected_ab[n], expected_abp[n];
//		size_t i;
//
//		for (carry = 0, i = 0; i < n; i++) {
//			wide_limb_t wide;
//			cc_add_carry_wide_limb_limb(&wide, a[i], b[i], carry);
//			expected_ab[i] = wide.lo;
//			carry = wide.hi;
//
//			if (expected_ab[i] != ab[i]) {
//				printf("Mismatch at limb %zu in a+b: expected %llu, got %llu\n", i, (unsigned long long)expected_ab[i], (unsigned long long)ab[i]);
//				EXPECT_TRUE(false);
//			}
//		}
//
//		for (borrow = 0, i = 0; i < n; i++) {
//			wide_limb_t wide;
//			cc_sub_borrow_wide_limb_limb(&wide, expected_ab[i], Fp::modulus()[i], borrow);
//			expected_abp[i] = wide.lo;
//			borrow = wide.hi & 1;
//
//
//			if (expected_abp[i] != abp[i]) {
//				printf("Mismatch at limb %zu in a+b-p: expected %llu, got %llu\n", i, (unsigned long long)expected_abp[i], (unsigned long long)abp[i]);
//				EXPECT_TRUE(false);
//			}
//		}
//
//	}
//
//}


TEST(ParallelPrefix, Add) {
	InsecureRng rng;

	for (u64 t = 0; t < 100; ++t)
	{

		Fp a, b, c;
		constexpr auto n = Fp::num_limb();

		rand_Fp(a, rng);
		rand_Fp(b, rng);
		Fp np = Fp::neg_modulus();
		Fp p = Fp::modulus();

		for (u64 i = 0; i < n; ++i)
		{
			ASSERT_EQ(np[i], ~p[i] + (i == 0 ? 1 : 0)) << " at " << i;
		}

		parallelPrefix_add<Fp>(c.data, a.data, b.data);

		auto expected = Fp();
		Fp::add(expected, a, b);
		if (!Fp::equal(c, expected)) {
			printf("Mismatch in pp_add:\n");
			printf("a = "); print_Fp(a);
			printf("b = "); print_Fp(b);
			printf("expected = "); print_Fp(expected);
			printf("got = "); print_Fp(c);
			EXPECT_TRUE(false);
		}

	}

}

TEST(ParallelPrefix, LaunchPpAddModAos) {
	InsecureRng rng;

	constexpr size_t n_elems = 100;
	constexpr int NLIMBS = 6; // For 384-bit field elements

	// Host arrays
	std::vector<limb_t> h_A(n_elems * NLIMBS);
	std::vector<limb_t> h_B(n_elems * NLIMBS);
	std::vector<limb_t> h_C(n_elems * NLIMBS);
	std::vector<limb_t> h_expected(n_elems * NLIMBS);

	// Generate test data

	// Fill arrays with random field elements in AoS layout
	for (size_t elem = 0; elem < n_elems; ++elem) {
		Fp a, b;
		rand_Fp(a, rng);
		rand_Fp(b, rng);

		// Copy to AoS layout: base + elem*NLIMBS + limb
		for (int limb = 0; limb < NLIMBS; ++limb) {
			size_t idx = elem * NLIMBS + limb;
			h_A[idx] = a[limb];
			h_B[idx] = b[limb];
		}

		// Compute expected result using CPU implementation
		Fp expected;
		Fp::add(expected, a, b);
		for (int limb = 0; limb < NLIMBS; ++limb) {
			h_expected[elem * NLIMBS + limb] = expected[limb];
		}
	}

	// Device arrays
	limb_t* d_A, * d_B, * d_C;
	size_t bytes = n_elems * NLIMBS * sizeof(limb_t);

	ASSERT_EQ(cudaMalloc(&d_A, bytes), cudaSuccess);
	ASSERT_EQ(cudaMalloc(&d_B, bytes), cudaSuccess);
	ASSERT_EQ(cudaMalloc(&d_C, bytes), cudaSuccess);

	// Copy data to device
	ASSERT_EQ(cudaMemcpy(d_A, h_A.data(), bytes, cudaMemcpyHostToDevice), cudaSuccess);
	ASSERT_EQ(cudaMemcpy(d_B, h_B.data(), bytes, cudaMemcpyHostToDevice), cudaSuccess);

	// Launch the kernel
	cudaError_t result = launch_pp_add_mod_aos<Fp>(d_A, d_B, d_C, n_elems);
	ASSERT_EQ(result, cudaSuccess) << "Kernel launch failed: " << cudaGetErrorString(result);

	// Wait for completion and check for runtime errors
	ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

	// Copy result back to host
	ASSERT_EQ(cudaMemcpy(h_C.data(), d_C, bytes, cudaMemcpyDeviceToHost), cudaSuccess);

	// Verify results
	for (size_t elem = 0; elem < n_elems; ++elem) {
		for (int limb = 0; limb < NLIMBS; ++limb) {
			size_t idx = elem * NLIMBS + limb;
			if (h_C[idx] != h_expected[idx]) {
				printf("Mismatch at element %zu, limb %d: expected %llu, got %llu\n",
					elem, limb, (unsigned long long)h_expected[idx], (unsigned long long)h_C[idx]);
				EXPECT_EQ(h_C[idx], h_expected[idx]);
			}
		}
	}

	// Clean up device memory
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);
}



// ---- paste your compress5to2_true implementation here ----
// (from the previous message)
// ----------------------------------------------------------

//// Helpers
//static inline int bit_at(limb_t* w, int b)
//{
//	auto l = b / (sizeof(limb_t) * 8);
//	b = b % (sizeof(limb_t) * 8);
//
//	return static_cast<int>((w[l] >> b) & limb_t(1));
//}

//// Fetch carry_out bit at column (i,b) from carry_shl1 (which is already <<1 across limbs).
//// carry_out(i,b) equals carry_shl1 bit at (i,b+1), except when b==W-1 it spills to (i+1,0).
//static inline int carry_out_bit(
//	const limb_t* carry_shl1,
//	std::size_t m,
//	int i,
//	int b,
//	limb_t top_carry_bit
//)
//{
//	constexpr int W = sizeof(limb_t) * 8;
//	if (b + 1 < W)
//	{
//		return bit_at(carry_shl1[i], b + 1);
//	}
//	if (i + 1 < static_cast<int>(m))
//	{
//		return bit_at(carry_shl1[i + 1], 0);
//	}
//	return static_cast<int>(top_carry_bit & 1);
//}

template<int m>
static void run_once_and_check(InsecureRng& rng)
{

	//limb_t cin = rng.next_u64() % 3;
	std::vector<limb_t> x1(m);
	std::vector<limb_t> x2(m);
	std::vector<limb_t> x3(m);
	std::vector<limb_t> x4(m);
	std::vector<limb_t> x5(m);
	std::vector<limb_t> exp(m);
	std::vector<limb_t> act(m);


	for (std::size_t i = 0; i < m; i++)
	{
		x1[i] = rng.next_u64();
		x2[i] = rng.next_u64();
		x3[i] = rng.next_u64();
		x4[i] = rng.next_u64();
		x5[i] = rng.next_u64();
	}

	for (u64 cin = 0; cin < 4; cin++)
	{

		auto ci = cin;
		auto add_carry = [&](auto x, auto y) {
			auto r = x + y;
			ci += (r < x) ? 1 : 0;
			return r;
			};


		for (std::size_t i = 0; i < m; i++)
		{
			exp[i] = std::exchange(ci, 0);
			//exp[i] = add_carry(exp[i], cin);
			exp[i] = add_carry(exp[i], x1[i]);
			exp[i] = add_carry(exp[i], x2[i]);
			exp[i] = add_carry(exp[i], x3[i]);
			exp[i] = add_carry(exp[i], x4[i]);
			exp[i] = add_carry(exp[i], x5[i]);
		}

		std::vector<limb_t> sum(m);
		std::vector<limb_t> carry(m);

		compress5to2<m, limb_t>(
			cin,
			x1.data(),
			x2.data(),
			x3.data(),
			x4.data(),
			x5.data(),
			sum.data(),
			carry.data()
		);

		ci = 0;
		for (std::size_t i = 0; i < m; i++)
		{
			act[i] = std::exchange(ci, 0);
			act[i] = add_carry(act[i], sum[i]);
			act[i] = add_carry(act[i], carry[i]);
		}


		for (int i = 0; i < static_cast<int>(m); i++)
		{
			if (exp[i] != act[i])
			{
				std::fprintf(
					stderr,
					"Mismatch at limb %d: expected %016llx, got %016llx\n",
					i,
					(unsigned long long)exp[i],
					(unsigned long long)act[i]
				);
				ASSERT_TRUE(false);
			}
		}
	}
}
//
//static void run_edge_cases()
//{
//	const std::size_t m = 3;
//
//	// All zeros
//	std::vector<limb_t> z(m, 0);
//	std::vector<limb_t> sum(m, 0);
//	std::vector<limb_t> cs(m, 0);
//	std::vector<limb_t> csh(m, 0);
//	limb_t top = 0;
//
//	compress5to2<limb_t>(
//		z.data(), z.data(), z.data(), z.data(), z.data(),
//		sum.data(), cs.data(), m);
//
//	for (std::size_t i = 0; i < m; i++)
//	{
//		ASSERT_TRUE(sum[i] == 0);
//		ASSERT_TRUE(cs[i] == 0);
//		ASSERT_TRUE(csh[i] == 0);
//	}
//	ASSERT_TRUE(top == 0);
//
//	// All ones (stress max local counts)
//	std::vector<limb_t> o(m, ~limb_t(0));
//	std::fill(sum.begin(), sum.end(), 0);
//	std::fill(cs.begin(), cs.end(), 0);
//	std::fill(csh.begin(), csh.end(), 0);
//	top = 0;
//
//	compress5to2<limb_t>(
//		o.data(), o.data(), o.data(), o.data(), o.data(),
//		sum.data(), cs.data(), m);
//
//	// For all-ones inputs, every column has count = 6.
//	// We expect s=0, cs=0, and carry_out=1 everywhere, so csh should be all 1s << 1 across the number.
//	// We validate with the per-bit identity below.
//	constexpr int W = sizeof(limb_t) * 8;
//	for (int i = 0; i < static_cast<int>(m); i++)
//	{
//		for (int b = 0; b < W; b++)
//		{
//			int cnt = 6;
//			int s = bit_at(sum[i], b);
//			int css = bit_at(cs[i], b);
//			int co = carry_out_bit(csh.data(), m, i, b, top);
//			int rhs = s + 2 * css + 2 * co;
//			if (cnt != rhs)
//			{
//				std::fprintf(
//					stderr,
//					"All-ones mismatch at limb=%d bit=%d: rhs=%d=%d+2 %d +2 %d"
//					"\n",
//					i,
//					b,
//					rhs,
//					s,
//					css,
//					co
//				);
//				ASSERT_TRUE(false);
//			}
//		}
//	}
//}
//
//TEST(ParallelPrefix, MultEdgeCase) {
//	// Edge cases first
//	run_edge_cases();
//}
TEST(ParallelPrefix, Mult)
{
	InsecureRng rng;
	// Randomized testing across various sizes
	//for (std::size_t m : {1ul, 2ul, 3ul, 4ul, 8ul})
	{
		for (int iter = 0; iter < 200; iter++)
		{
			run_once_and_check<1>(rng);
			run_once_and_check<2>(rng);
			run_once_and_check<3>(rng);
			run_once_and_check<4>(rng);
			run_once_and_check<5>(rng);
		}
	}

}

// ==== Configure what to test here ====

//// Forward decls for the two functions under test:
//template<typename Field, int n, typename Carry = std::monostate>
//uint64_t montgomery_fios_step(
//	uint64_t* t,
//	const uint64_t* a,
//	uint64_t bj,
//	Carry carry_in);
//
//template<typename Field, int n, typename Carry = std::monostate>
//uint64_t montgomery_fios_step_lo_down_emulated(
//	uint64_t* t,
//	const uint64_t* a,
//	uint64_t bj,
//	Carry carry_in = {});

// ===== deterministic RNG (splitmix64) =====
struct SplitMix64 {
	uint64_t x;
	explicit SplitMix64(uint64_t seed) : x(seed) {}
	uint64_t next() {
		uint64_t z = (x += 0x9E3779B97F4A7C15ull);
		z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
		z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
		return z ^ (z >> 31);
	}
};

// ===== helpers =====
template<int n>
static void random_vec(uint64_t* out, SplitMix64& rng) {
	for (int i = 0; i < n; i++)
		out[i] = rng.next();
}

template<int n>
static void set_all(uint64_t* out, uint64_t v) {
	for (int i = 0; i < n; i++)
		out[i] = v;
}

// Compare first n limbs of t and the returned carry.
// (The reference’s top limb effect should be visible in t[n-1] and in the returned carry.)
template<int n>
static void expect_equal_state(
	const uint64_t* t_ref,
	const uint64_t* t_emu,
	uint64_t carry_ref,
	uint64_t carry_emu)
{
	for (int i = 0; i < n; i++)
		EXPECT_EQ(t_ref[i], t_emu[i]) << "mismatch at limb " << i;
	EXPECT_EQ(carry_ref, carry_emu) << "final carry mismatch";
}

// ==== Test Fixture ====
class FiosStepTest : public ::testing::Test {
protected:
	static constexpr int n = 6;
};


// a single Montgomery fused limb step (FIOS-style) with optional carry in/out.
// Input: digit b_j, multiplicand a[0..n-1], accumulator t[0..n], optional carry_in
// Output: updated t[0..n], returns carry_out
// If no carry_in/out is needed, use std::monostate as Carry.
// 
// compute t = t + a * b_j
// compute m = t[0] * n0
// compute t = (t + m * p) >> 64
// compute update carry_out with carry_in if any
//
template<typename Field, int n, typename Carry = std::monostate>
uint64_t montgomery_fios_test(
	bool print, 
	uint64_t* t,
	const uint64_t* a,
	uint64_t b_j,
	Carry carry_in = {})
{
	using namespace blst;

	constexpr bool has_carry = !std::is_same<Carry, std::monostate>::value;
	constexpr auto p = Field::modulus().data;
	constexpr uint64_t n0 = Field::montgomery_constant();

	// tmp = b[0] * a
	uint64_t hi_ab;
	limb_t mx, hi_mp;
	{
		auto a0 = a[0];
		auto t0 = t[0];
		auto p0 = p[0];

		if constexpr (has_carry)
		{
			// tmp = ab+t
			wide_limb_t  ab_0;
			cc_mul_add_wide_limb_limb(&ab_0, b_j, a0, t0);

			mx = n0 * ab_0.lo;
			// tmp = mx * p + tmp
			wide_limb_t  mp_0;
			cc_mul_add_wide_limb_limb(&mp_0, mx, p[0], ab_0.lo);

			if (print)
				std::cout << "mp = (" << mp_0.hi << "," << mp_0.lo << ") = " << mx << "*" << p[0] 
					<< "+" << ab_0.lo <<"+" << t0 << std::endl;
			hi_mp = mp_0.hi;
			t[0] = mp_0.lo;
			hi_ab = ab_0.hi;
		}
		else
		{
			wide_limb_t  ab_0;
			wide_limb_t  mp_0;
			// tmp = ab
			cc_mul_wide_limb_limb(&ab_0, b_j, a0);
			mx = n0 * ab_0.lo;
			// tmp = mx * p + tmp
			cc_mul_add_wide_limb_limb(&mp_0, mx, p[0], ab_0.lo);
			hi_mp = mp_0.hi;
			if (print)
				std::cout << "mp = (" << mp_0.hi << ","<<mp_0.lo<<") = " << mx <<"*" << p[0] << "+" << ab_0.lo << std::endl;
			t[0] = mp_0.lo;
			hi_ab = ab_0.hi;
		}
	}

	if (print)
		std::cout << "mx = " << mx << "\n";

	for (uint64_t i = 1; i < n; i++) {
		wide_limb_t  ab_i;
		wide_limb_t  mp_i;
		auto ai = a[i];
		auto pi = p[i];
		auto ti = t[i];

		if constexpr (has_carry)
		{
			ab_i.lo = b_j * ai;
			ab_i.hi = umul64hi(b_j, ai);
			mp_i.lo = mx * pi;
			mp_i.hi = umul64hi(mx, pi);

			if (print)
				std::cout << "i " << i - 1 << ":"
				<< " pm.lo " << mp_i.lo
				<< " + ti " << ti
				<< " + pm.hi " << hi_mp
				<< " + ab.lo " << ab_i.lo
				<< " + ab.hi " << hi_ab
				<< "\n";
			auto ahi = ab_i.hi;
			auto phi = mp_i.hi;
			cc_add_wide_wide_limb(&ab_i, ab_i, hi_ab);
			cc_add_wide_wide_limb(&ab_i, ab_i, ti);
			cc_add_wide_wide_limb(&mp_i, mp_i, hi_mp);
			cc_add_wide_wide_limb(&mp_i, mp_i, ab_i.lo);

			auto ab_carry = ab_i.hi - ahi;
			auto mp_carry = mp_i.hi - phi;
			if (print)
				std::cout << " ab_carry " << ab_carry << " mp_carry " << mp_carry << std::endl;
		}
		else
		{
			ab_i.lo = b_j * ai;
			ab_i.hi = umul64hi(b_j, ai);
			mp_i.lo = mx * pi;
			mp_i.hi = umul64hi(mx, pi);
			if (print)
				std::cout << "i " << i - 1 << ":"
				<< " pm.lo " << mp_i.lo
				<< " + pm.hi " << hi_mp
				<< " + ab.lo " << ab_i.lo
				<< " + ab.hi " << hi_ab
				<< "\n";

			cc_add_wide_wide_limb(&ab_i, ab_i, hi_ab);
			cc_add_wide_wide_limb(&mp_i, mp_i, hi_mp);
			cc_add_wide_wide_limb(&mp_i, mp_i, ab_i.lo);
		}
		t[i - 1] = mp_i.lo;
		if (print)
			std::cout << "sum[" << i - 1 << "] = " << t[i - 1] << std::endl;
		hi_ab = ab_i.hi;
		hi_mp = mp_i.hi;
	}
	//t[n] = hi_ab;

	{
		wide_limb_t wide1;

		if constexpr (has_carry) {
			// wide1 = hi_ab + hi_mp + carry_in 
			//wide_limb_t wide;
			//cc_add_wide_limb_limb(&wide, hi_ab, carry_in);
			//cc_add_carry_wide_limb_limb(&wide1, hi_mp, wide.lo, wide.hi);
			// wide1 = hi_ab + hi_mp + carry_in 
			wide_limb_t wide_ab;
			cc_add_wide_limb_limb(&wide_ab, hi_ab, carry_in);

			// Add hi_mp to the low part
			cc_add_wide_limb_limb(&wide1, wide_ab.lo, hi_mp);

			// Add the carry from the first addition to the high part
			wide1.hi += wide_ab.hi;


			if (print)
				std::cout << "i " << n-1 << ":"
				<< " hi_mp " << hi_mp
				<< " + hi_ab " << hi_ab
				<< " + carry_in " << carry_in
				<< "\n";
		}
		else {
			// wide1 = hi_ab + hi_mp
			cc_add_wide_limb_limb(&wide1, hi_mp, hi_ab);
		}

		t[n - 1] = wide1.lo;
		if (print)
			std::cout << "sum[" << n - 1 << "] = " << t[n - 1] << std::endl;
		if (print)
			std::cout << "carry " << wide1.hi << std::endl;
		return wide1.hi;
	}
}

// ============ Tests ============

TEST_F(FiosStepTest, AllZeros_NoCarry) {
	alignas(64) uint64_t a[n];
	alignas(64) uint64_t t_ref[n];
	alignas(64) uint64_t t_emu[n];

	set_all<n>(a, 0ull);
	set_all<n>(t_ref, 0ull);
	std::copy(t_ref, t_ref + n, t_emu);

	uint64_t bj = 0ull;
	bool print = false;
	uint64_t carry_ref = montgomery_fios_test<Fp, n>(print, t_ref, a, bj, std::monostate{});
	uint64_t carry_emu = montgomery_fios_step_lo_down_emulated<Fp, n>(t_emu, a, bj, std::monostate{});

	expect_equal_state<n>(t_ref, t_emu, carry_ref, carry_emu);
}

TEST_F(FiosStepTest, AllOnes_NoCarry) {
	alignas(64) uint64_t a[n];
	alignas(64) uint64_t t_ref[n];
	alignas(64) uint64_t t_emu[n];

	set_all<n>(a, ~0ull);
	set_all<n>(t_ref, ~0ull);
	std::copy(t_ref, t_ref + n, t_emu);

	uint64_t bj = ~0ull;
	bool print = false;


	uint64_t carry_ref = montgomery_fios_test<Fp, n>(print, t_ref, a, bj, std::monostate{});
	uint64_t carry_emu = montgomery_fios_step_lo_down_emulated<Fp, n>(t_emu, a, bj, std::monostate{});

	expect_equal_state<n>(t_ref, t_emu, carry_ref, carry_emu);
}

TEST_F(FiosStepTest, Random_NoCarry) {
	SplitMix64 rng(0x12345678abcdef00ull);

	for (int trial = 0; trial < 10; trial++) {
		alignas(64) uint64_t a[n];
		alignas(64) uint64_t t_ref[n];
		alignas(64) uint64_t t_emu[n];

		random_vec<n>(a, rng);
		random_vec<n>(t_ref, rng);
		std::copy(t_ref, t_ref + n, t_emu);

		uint64_t bj = rng.next();
		bool print = false;

		uint64_t carry_ref = montgomery_fios_test<Fp, n>(print, t_ref, a, bj, std::monostate{});
		uint64_t carry_emu = montgomery_fios_step_lo_down_emulated<Fp, n>(t_emu, a, bj, std::monostate{});

		expect_equal_state<n>(t_ref, t_emu, carry_ref, carry_emu);
	}
}

TEST_F(FiosStepTest, Random_WithCarryBit01) {
	SplitMix64 rng(0xCAFEBABEDEADBEEFull);
	bool print = false;

	for (int trial = 0; trial < 10; trial++) {
		alignas(64) uint64_t a[n];
		alignas(64) uint64_t t_ref[n];
		alignas(64) uint64_t t_emu[n];

		random_vec<n>(a, rng);
		random_vec<n>(t_ref, rng);
		std::copy(t_ref, t_ref + n, t_emu);

		uint64_t bj = rng.next();
		uint64_t carry_in =  rng.next() & 1ull;   // emulate carry bit

		uint64_t carry_ref = montgomery_fios_test<Fp, n>(print, t_ref, a, bj, carry_in);
		uint64_t carry_emu = montgomery_fios_step_lo_down_emulated<Fp, n>(t_emu, a, bj, carry_in);

		expect_equal_state<n>(t_ref, t_emu, carry_ref, carry_emu);
	}
}

TEST_F(FiosStepTest, Random_WithWideCarry) {
	SplitMix64 rng(0x0ddc0ffeebad5eedull);
	bool print = false;

	for (int trial = 0; trial < 10; trial++) {
		alignas(64) uint64_t a[n];
		alignas(64) uint64_t t_ref[n];
		alignas(64) uint64_t t_emu[n];

		random_vec<n>(a, rng);
		random_vec<n>(t_ref, rng);
		std::copy(t_ref, t_ref + n, t_emu);

		uint64_t bj = rng.next();
		uint64_t carry_in = ~0ull;// rng.next();    // allow any 64-bit (your serial tail supports this)

		uint64_t carry_ref = montgomery_fios_test<Fp, n, uint64_t>(print, t_ref, a, bj, carry_in);
		uint64_t carry_emu = montgomery_fios_step_lo_down_emulated<Fp, n, uint64_t>(t_emu, a, bj, carry_in);

		expect_equal_state<n>(t_ref, t_emu, carry_ref, carry_emu);
	}
}

// Stress patterns engineered to cascade carries
TEST_F(FiosStepTest, CarryCascadePatterns) {
	alignas(64) uint64_t a[n];
	alignas(64) uint64_t t_ref[n];
	alignas(64) uint64_t t_emu[n];
	bool print = false;

	// Pattern 1: alternating all-ones / zeros
	for (int i = 0; i < n; i++)
		a[i] = (i & 1) ? ~0ull : 0ull;

	for (int i = 0; i < n; i++)
		t_ref[i] = (i & 1) ? 0ull : ~0ull;

	std::copy(t_ref, t_ref + n, t_emu);
	uint64_t bj = ~0ull;
	uint64_t carry_in = 1ull;

	uint64_t carry_ref = montgomery_fios_test<Fp, n, uint64_t>(print, t_ref, a, bj, carry_in);
	uint64_t carry_emu = montgomery_fios_step_lo_down_emulated<Fp, n, uint64_t>(t_emu, a, bj, carry_in);

	expect_equal_state<n>(t_ref, t_emu, carry_ref, carry_emu);

	// Pattern 2: high words near overflow
	for (int i = 0; i < n; i++)
		a[i] = 0xFFFFFFFFFFFFFFFEull;

	for (int i = 0; i < n; i++)
		t_ref[i] = 0xFFFFFFFFFFFFFFFEull;

	std::copy(t_ref, t_ref + n, t_emu);
	bj = 0xFFFFFFFFFFFFFFFFull;
	carry_in = 0xFFFFFFFFFFFFFFFFull;

	carry_ref = montgomery_fios_test<Fp, n, uint64_t>(print, t_ref, a, bj, carry_in);
	carry_emu = montgomery_fios_step_lo_down_emulated<Fp, n, uint64_t>(t_emu, a, bj, carry_in);

	expect_equal_state<n>(t_ref, t_emu, carry_ref, carry_emu);
}

// --- Helpers ---
void fill_random(std::vector<uint64_t>& v, InsecureRng& rng) {
	for (auto& x : v) x = rng.next_u64();
}

TEST(ParallelPrefix, Compress5to2_Tile_Check) {
	InsecureRng rng;
	constexpr int NLIMBS = 6;
	constexpr int NUM_OPS = 128; // Test 128 instances in parallel

	// 1. Host Data
	std::vector<uint64_t> h_x0(NUM_OPS * NLIMBS), h_x1(NUM_OPS * NLIMBS),
		h_x2(NUM_OPS * NLIMBS), h_x3(NUM_OPS * NLIMBS),
		h_x4(NUM_OPS * NLIMBS);
	std::vector<uint64_t> h_cin(NUM_OPS);

	// Outputs
	std::vector<uint64_t> cpu_s(NUM_OPS * NLIMBS), cpu_c(NUM_OPS * NLIMBS), cpu_spill(NUM_OPS);
	std::vector<uint64_t> gpu_s(NUM_OPS * NLIMBS), gpu_c(NUM_OPS * NLIMBS), gpu_spill(NUM_OPS);

	// Initialize
	fill_random(h_x0, rng); fill_random(h_x1, rng); fill_random(h_x2, rng);
	fill_random(h_x3, rng); fill_random(h_x4, rng);
	for (auto& x : h_cin) x = rng.next_u64() % 4; // cin must be < 4

	// 2. Run CPU Emulation
	for (int i = 0; i < NUM_OPS; ++i) {
		int off = i * NLIMBS;
		cpu_spill[i] = compress5to2<NLIMBS>(
			h_cin[i],
			&h_x0[off], &h_x1[off], &h_x2[off], &h_x3[off], &h_x4[off],
			&cpu_s[off], &cpu_c[off]
		);
	}

	// 3. Run GPU Kernel
	uint64_t* d_x0, * d_x1, * d_x2, * d_x3, * d_x4, * d_cin, * d_s, * d_c, * d_spill;
	size_t vec_bytes = NUM_OPS * NLIMBS * sizeof(uint64_t);
	size_t scal_bytes = NUM_OPS * sizeof(uint64_t);

	cudaMalloc(&d_x0, vec_bytes); cudaMemcpy(d_x0, h_x0.data(), vec_bytes, cudaMemcpyHostToDevice);
	cudaMalloc(&d_x1, vec_bytes); cudaMemcpy(d_x1, h_x1.data(), vec_bytes, cudaMemcpyHostToDevice);
	cudaMalloc(&d_x2, vec_bytes); cudaMemcpy(d_x2, h_x2.data(), vec_bytes, cudaMemcpyHostToDevice);
	cudaMalloc(&d_x3, vec_bytes); cudaMemcpy(d_x3, h_x3.data(), vec_bytes, cudaMemcpyHostToDevice);
	cudaMalloc(&d_x4, vec_bytes); cudaMemcpy(d_x4, h_x4.data(), vec_bytes, cudaMemcpyHostToDevice);
	cudaMalloc(&d_cin, scal_bytes); cudaMemcpy(d_cin, h_cin.data(), scal_bytes, cudaMemcpyHostToDevice);

	cudaMalloc(&d_s, vec_bytes);
	cudaMalloc(&d_c, vec_bytes);
	cudaMalloc(&d_spill, scal_bytes);

	// Launch
	int threads = 256;
	int warps = threads / 32;
	int tiles_per_warp = 32 / NLIMBS;
	int ops_per_block = warps * tiles_per_warp;
	int grids = (NUM_OPS + ops_per_block - 1) / ops_per_block;

	compress5to2_test_kernel<NLIMBS> << <grids, threads >> > (
		d_s, d_c, d_spill, d_cin,
		d_x0, d_x1, d_x2, d_x3, d_x4, NUM_OPS
		);
	cudaDeviceSynchronize();

	// Readback
	cudaMemcpy(gpu_s.data(), d_s, vec_bytes, cudaMemcpyDeviceToHost);
	cudaMemcpy(gpu_c.data(), d_c, vec_bytes, cudaMemcpyDeviceToHost);
	cudaMemcpy(gpu_spill.data(), d_spill, scal_bytes, cudaMemcpyDeviceToHost);

	// 4. Validate
	for (int i = 0; i < NUM_OPS; ++i) {
		// Check Spill
		ASSERT_EQ(cpu_spill[i], gpu_spill[i]) << "Spill mismatch at op " << i;

		// Check S and C vectors
		for (int j = 0; j < NLIMBS; ++j) {
			int idx = i * NLIMBS + j;
			ASSERT_EQ(cpu_s[idx], gpu_s[idx]) << "S mismatch at op " << i << " limb " << j;
			ASSERT_EQ(cpu_c[idx], gpu_c[idx]) << "C mismatch at op " << i << " limb " << j;
		}
	}

	// Cleanup
	cudaFree(d_x0); cudaFree(d_x1); cudaFree(d_x2); cudaFree(d_x3); cudaFree(d_x4);
	cudaFree(d_cin); cudaFree(d_s); cudaFree(d_c); cudaFree(d_spill);
}


TEST(ParallelPrefix, FiosStep_Tile_Check) {
	InsecureRng rng;
	constexpr int NLIMBS = 6;
	constexpr int NUM_OPS = 128;

	// We need a Field definition for P and n0. Using blst::Fp logic.
	using Field = blst::Fp;
	uint64_t n0 = Field::montgomery_constant();

	// 1. Host Data
	std::vector<uint64_t> h_t(NUM_OPS * NLIMBS);
	std::vector<uint64_t> h_a(NUM_OPS * NLIMBS);
	std::vector<uint64_t> h_p(NUM_OPS * NLIMBS);

	std::vector<uint64_t> h_b(NUM_OPS);       // Scalar b_j
	std::vector<uint64_t> h_cin(NUM_OPS);     // Scalar carry_in

	// Initialize
	fill_random(h_t, rng);
	fill_random(h_a, rng);
	for (auto& x : h_b) x = rng.next_u64();
	for (auto& x : h_cin) x = rng.next_u64(); // Full 64-bit carry in allowed

	// Fill P (constant for all ops, but AoS layout)
	for (int i = 0; i < NUM_OPS; ++i) {
		for (int j = 0; j < NLIMBS; ++j) {
			h_p[i * NLIMBS + j] = Field::modulus().data[j];
		}
	}

	// 2. Run CPU Emulation
	std::vector<uint64_t> cpu_t = h_t; // Copy initial state
	std::vector<uint64_t> cpu_cout(NUM_OPS);

	for (int i = 0; i < NUM_OPS; ++i) {
		int off = i * NLIMBS;
		// The CPU emulation updates 't' in place
		cpu_cout[i] = montgomery_fios_step_lo_down_emulated<Field, NLIMBS, uint64_t>(
			&cpu_t[off],
			&h_a[off],
			h_b[i],
			h_cin[i]
		);
	}

	// 3. Run GPU Kernel
	uint64_t* d_t, * d_a, * d_p, * d_b, * d_cin, * d_cout;
	size_t vec_bytes = NUM_OPS * NLIMBS * sizeof(uint64_t);
	size_t scal_bytes = NUM_OPS * sizeof(uint64_t);

	cudaMalloc(&d_t, vec_bytes); cudaMemcpy(d_t, h_t.data(), vec_bytes, cudaMemcpyHostToDevice);
	cudaMalloc(&d_a, vec_bytes); cudaMemcpy(d_a, h_a.data(), vec_bytes, cudaMemcpyHostToDevice);
	cudaMalloc(&d_p, vec_bytes); cudaMemcpy(d_p, h_p.data(), vec_bytes, cudaMemcpyHostToDevice);
	cudaMalloc(&d_b, scal_bytes); cudaMemcpy(d_b, h_b.data(), scal_bytes, cudaMemcpyHostToDevice);
	cudaMalloc(&d_cin, scal_bytes); cudaMemcpy(d_cin, h_cin.data(), scal_bytes, cudaMemcpyHostToDevice);
	cudaMalloc(&d_cout, scal_bytes);

	int threads = 256;
	int warps = threads / 32;
	int tiles_per_warp = 32 / NLIMBS;
	int ops_per_block = warps * tiles_per_warp;
	int grids = (NUM_OPS + ops_per_block - 1) / ops_per_block;

	fios_step_test_kernel<NLIMBS> << <grids, threads >> > (
		d_t, d_cout,
		d_a, d_p, d_b, d_cin, n0, NUM_OPS
		);
	cudaDeviceSynchronize();

	// Readback
	std::vector<uint64_t> gpu_t(NUM_OPS * NLIMBS);
	std::vector<uint64_t> gpu_cout(NUM_OPS);

	cudaMemcpy(gpu_t.data(), d_t, vec_bytes, cudaMemcpyDeviceToHost);
	cudaMemcpy(gpu_cout.data(), d_cout, scal_bytes, cudaMemcpyDeviceToHost);

	// 4. Validate
	for (int i = 0; i < NUM_OPS; ++i) {
		ASSERT_EQ(cpu_cout[i], gpu_cout[i]) << "Carry Out mismatch at op " << i;

		for (int j = 0; j < NLIMBS; ++j) {
			int idx = i * NLIMBS + j;
			ASSERT_EQ(cpu_t[idx], gpu_t[idx]) << "T mismatch at op " << i << " limb " << j;
		}
	}

	cudaFree(d_t); cudaFree(d_a); cudaFree(d_p);
	cudaFree(d_b); cudaFree(d_cin); cudaFree(d_cout);
}

TEST(ParallelPrefix, MulMontCoop) {
	InsecureRng rng;
	constexpr int NLIMBS = 6; // Testing 384-bit
	constexpr int NUM_OPS = 1024;

	using Field = blst::Fp;

	// 1. Host Setup
	std::vector<limb_t> h_a(NUM_OPS * NLIMBS);
	std::vector<limb_t> h_b(NUM_OPS * NLIMBS);
	std::vector<limb_t> h_r_expected(NUM_OPS * NLIMBS);
	std::vector<limb_t> h_r_gpu(NUM_OPS * NLIMBS);

	for (int i = 0; i < NUM_OPS; ++i) {
		Fp a, b, r;
		rand_Fp(a, rng);
		rand_Fp(b, rng);

		// Golden result
		Fp::mul(r, a, b);

		// Pack into AoS buffers
		for (int j = 0; j < NLIMBS; ++j) {
			h_a[i * NLIMBS + j] = a[j];
			h_b[i * NLIMBS + j] = b[j];
			h_r_expected[i * NLIMBS + j] = r[j];
		}
	}

	// 2. Device Setup
	limb_t* d_a, * d_b, * d_r;
	size_t size = NUM_OPS * NLIMBS * sizeof(limb_t);
	cudaMalloc(&d_a, size);
	cudaMalloc(&d_b, size);
	cudaMalloc(&d_r, size);

	cudaMemcpy(d_a, h_a.data(), size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, h_b.data(), size, cudaMemcpyHostToDevice);
	cudaMemset(d_r, 0, size);

	// 3. Launch
	int threads_per_block = 256;
	int warps_per_block = threads_per_block / 32;
	int tiles_per_warp = 32 / NLIMBS;
	int ops_per_block = warps_per_block * tiles_per_warp;
	int grid_size = (NUM_OPS + ops_per_block - 1) / ops_per_block;
	
	blst::mul_mont_coop_kernel<Field> <<<grid_size, threads_per_block >>> (d_r, d_a, d_b, NUM_OPS);

	cudaError_t err = cudaGetLastError();
	ASSERT_EQ(err, cudaSuccess) << "Kernel launch failed: " << cudaGetErrorString(err);
	cudaDeviceSynchronize();

	Fp rr;
	par_mul_mont_u64<Fp>(rr.data, h_a.data(), h_b.data());

	// 4. Verify
	cudaMemcpy(h_r_gpu.data(), d_r, size, cudaMemcpyDeviceToHost);

	int errors = 0;
	for (int i = 0; i < NUM_OPS; ++i) {
		bool mismatch = false;
		for (int j = 0; j < NLIMBS; ++j) {
			if (h_r_gpu[i * NLIMBS + j] != h_r_expected[i * NLIMBS + j]) {
				mismatch = true;
			}
		}
		if (mismatch) {
			if (errors < 5) { // Print first 5 errors
				printf("Mismatch at op %d\n", i);
				printf("Expected: ");
				for (int k = NLIMBS - 1; k >= 0; k--) printf("%016llx ", (unsigned long long)h_r_expected[i * NLIMBS + k]);
				printf("\nGot:      ");
				for (int k = NLIMBS - 1; k >= 0; k--) printf("%016llx ", (unsigned long long)h_r_gpu[i * NLIMBS + k]);
				printf("\n");
			}
			errors++;
		}
	}
	EXPECT_EQ(errors, 0) << "Found " << errors << " mismatches out of " << NUM_OPS;

	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_r);
}