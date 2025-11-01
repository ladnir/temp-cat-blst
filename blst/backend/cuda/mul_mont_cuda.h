
#include <stdint.h>
#include "blst/common/defines.h"
#include "blst/backend/backend_portable.h"
#include <variant>

namespace blst
{
	//// ===== Helper functions for Montgomery multiply =====
	//// r[0..n] = a[0..n-1] * b + carry
	//// 
	//// if there is a carry in, its added to the last limb as part of 
	//// the montgomery update step. The new carry out is returned.
	//template<int n, typename Carry = std::monostate>
	//BLST_HD BLST_FORCEINLINE
	//Carry mul_add_chain(
	//	uint64_t* BLST_restrict r,
	//	const uint64_t* BLST_restrict a,
	//	uint64_t b,
	//	Carry carry = {})
	//{
	//	using namespace blst;
	//
	//	// the first add chain does not have carry input
	//	constexpr bool has_carry = !std::is_same<Carry, std::monostate>::value;
	//
	//	// tmp = b[0] * a
	//	uint64_t hi = 0;
	//	{
	//		wide_limb_t  wide;
	//		if constexpr (has_carry)
	//		{
	//			cc_mul_add_wide_limb_limb(&wide, b, a[0], r[0]);
	//		}
	//		else
	//		{
	//			cc_mul_wide_limb_limb(&wide, b, a[0]);
	//		}
	//		r[0] = wide.lo;
	//		hi = wide.hi;
	//	}
	//#pragma unroll
	//	for (uint64_t i = 1; i < n; i++) {
	//		wide_limb_t  wide;
	//
	//		if constexpr (has_carry)
	//		{
	//			cc_mul_add_add_wide_limb_limb(&wide, b, a[i], hi, r[i]);
	//		}
	//		else
	//		{
	//			cc_mul_add_wide_limb_limb(&wide, b, a[i], hi);
	//		}
	//		r[i] = wide.lo;
	//		hi = wide.hi;
	//	}
	//	r[n] = hi;
	//
	//	if constexpr (has_carry) {
	//		wide_limb_t wide;
	//
	//		cc_add_wide_limb_limb(&wide, r[n], carry);
	//		r[n] = wide.lo;
	//		return wide.hi;
	//	}
	//	else {
	//		return {};
	//	}
	//}
	//
	//template<typename Field, int n, typename Carry = std::monostate>
	//BLST_HD BLST_FORCEINLINE
	//uint64_t mont_update_chain(
	//	uint64_t* BLST_restrict r,
	//	Carry carry = {})
	//{
	//	using namespace blst;
	//
	//	constexpr auto p = Field::modulus().data;
	//	constexpr uint64_t n0 = Field::montgomery_constant();
	//	limb_t mx = n0 * r[0];
	//
	//	// tmp = mx * p + tmp
	//	wide_limb_t  wide;
	//	cc_mul_add_wide_limb_limb(&wide, mx, p[0], r[0]);
	//	limb_t hi = wide.hi;
	//	r[0] = wide.lo;
	//#pragma unroll
	//	for (uint64_t i = 1; i < n; i++) {
	//		cc_mul_add_add_wide_limb_limb(&wide, mx, p[i], hi, r[i]);
	//		r[i - 1] = wide.lo;
	//		hi = wide.hi;
	//	}
	//	if constexpr (!std::is_same<Carry, std::monostate>::value) {
	//		cc_add_carry_wide_limb_limb(&wide, hi, r[n], carry);
	//	}
	//	else {
	//		cc_add_wide_limb_limb(&wide, hi, r[n]);
	//	}
	//	// final limb
	//	//cc_add_carry_wide_limb_limb(&wide, hi, r[n], carry);
	//	r[n - 1] = wide.lo;
	//	return wide.hi;
	//}
	//
	//
	//// ===== Montgomery multiply (n x 64 CIOS, no __int128) =====
	//// r = a * b * R^{-1} mod p, base = 2^64
	//template<typename Field>
	//BLST_HD BLST_FORCEINLINE
	//void ptx_mul_mont_u64_simplified(
	//	uint64_t* BLST_restrict r,
	//	const uint64_t* BLST_restrict a,
	//	const uint64_t* BLST_restrict b)
	//{
	//
	//	constexpr auto p = Field::modulus().data;
	//	constexpr int n = sizeof(Field) / sizeof(uint64_t);
	//
	//	static_assert(n > 0 && n % 2 == 0, "n must be positive and even");
	//	using namespace blst;
	//
	//	//limb_t mask, borrow, mx, hi, carry;
	//	limb_t  tmp[n + 1];
	//
	//	// tmp[0...n] = a[0..n-1] * b[0]
	//	mul_add_chain<n>(tmp, a, b[0]);
	//
	//	// montgomery update: tmp = (tmp + m * p) >> 64
	//	limb_t carry = mont_update_chain<Field, n>(tmp);
	//
	//	for (size_t j = 1; j < n; ++j) {
	//		//wide_limb_t  wide;
	//
	//		// tmp += a * b[j] + carry
	//		carry = mul_add_chain<n>(tmp, a, b[j], carry);
	//
	//		//cc_add_wide_limb_limb(&wide, tmp[n], carry);
	//		//tmp[n] = wide.lo;
	//		//carry = wide.hi;
	//
	//		carry = mont_update_chain<Field, n>(tmp, carry);
	//	}
	//
	//	// conditional subtract: tmp - p and shift
	//	{
	//		limb_t borrow = 0;
	//		for (int i = 0; i < n; i++) {
	//			// wide = tmp - (p + borrow)
	//			wide_limb_t wide;
	//			cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
	//			r[i] = wide.lo;
	//			borrow = wide.hi & 1;
	//		}
	//
	//		limb_t mask = carry - borrow;
	//
	//		// select between tmp and ret
	//		for (int i = 0; i < n; i++)
	//			r[i] = (r[i] & ~mask) | (tmp[i] & mask);
	//	}
	//
	//}
#ifdef __CUDA_ARCH__
#define PTX_MUL_IMPL_CUDA
#endif

	BLST_HD BLST_FORCEINLINE
		uint64_t umul64hi(uint64_t a, uint64_t b) {
#ifdef __CUDA_ARCH__
		return __umul64hi(a, b);
#else
		uint64_t hi;
		static_assert(sizeof(unsigned long long) == sizeof(uint64_t), "uint64_t size mismatch");
		_mulx_u64(a, b, (unsigned long long*) & hi);
		return hi;
#endif
	}


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
	BLST_HD BLST_FORCEINLINE
		uint64_t montgomery_fios_step(
			uint64_t*  t,
			const uint64_t*  a,
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
#if defined(PTX_MUL_IMPL_CUDA) && 1
				wide_limb_t  ab_0;
				wide_limb_t  mp_0;
				asm volatile(
					"mul.lo.u64  %0, %6, %7;\n\t" // ab_lo = b_j * a0
					"mul.hi.u64  %1, %6, %7;\n\t" // ab_hi = umul64hi(b_j, a0)
					"add.cc.u64  %0, %0, %5;\n\t" // ab_lo += t0
					"addc.u64    %1, %1, 0;\n\t"  // ab_hi += carry
					"mul.lo.u64  %2, %0, %8;\n\t"  // mx = ab_0.lo * n0 
					"mul.lo.u64  %3, %2, %9;\n\t" // lo = mx * p0
					"mul.hi.u64  %4, %2, %9;\n\t" // hi = umul64hi(mx, p0)
					"add.cc.u64 %3, %3, %0;\n\t"  // lo += ab_0.lo
					"addc.u64   %4, %4, 0;\n\t"   // hi += carry
					: "=&l"(ab_0.lo), // 0
					"=&l"(ab_0.hi),  // 1
					"=&l"(mx),   // 2
					"=&l"(mp_0.lo), // 3
					"=&l"(mp_0.hi) // 4
					: "l"(t0),  // 5
					"l"(b_j), // 6
					"l"(a0),   // 7
					"l"(n0),    // 8
					"l"(p0)    // 9  
					);

				hi_mp = mp_0.hi;
				t[0] = mp_0.lo;
				hi_ab = ab_0.hi;
#else
				wide_limb_t  ab_0;
				cc_mul_add_wide_limb_limb(&ab_0, b_j, a0, t0);

				mx = n0 * ab_0.lo;
				// tmp = mx * p + tmp
				wide_limb_t  mp_0;
				cc_mul_add_wide_limb_limb(&mp_0, mx, p[0], ab_0.lo);
				hi_mp = mp_0.hi;
				t[0] = mp_0.lo;
				hi_ab = ab_0.hi;
#endif
			}
			else
			{

				// this block is problematic
#if defined(PTX_MUL_IMPL_CUDA) && 1
				wide_limb_t  ab_0;
				wide_limb_t  mp_0;
				asm volatile(
					"mul.lo.u64  %0, %5, %6;\n\t" // ab_lo = b_j * a0
					"mul.hi.u64  %1, %5, %6;\n\t" // ab_hi = umul64hi(b_j, a0)
					"mul.lo.u64  %2, %0, %7;\n\t"  // mx = ab_0.lo * n0 
					"mul.lo.u64  %3, %2, %8;\n\t" // lo = mx * p0
					"mul.hi.u64  %4, %2, %8;\n\t" // hi = umul64hi(mx, p0)
					"add.cc.u64 %3, %3, %0;\n\t"  // lo += ab_0.lo
					"addc.u64   %4, %4, 0;\n\t"   // hi += carry
					: "=&l"(ab_0.lo), // 0
					"=&l"(ab_0.hi),  // 1
					"=&l"(mx),   // 2
					"=&l"(mp_0.lo), // 3
					"=&l"(mp_0.hi) // 4
					:
					"l"(b_j), // 5
					"l"(a0),   // 6
					"l"(n0),    // 7
					"l"(p0)    // 8 
					);

				hi_mp = mp_0.hi;
				t[0] = mp_0.lo;
				hi_ab = ab_0.hi;
#else
				wide_limb_t  ab_0;
				wide_limb_t  mp_0;
				cc_mul_wide_limb_limb(&ab_0, b_j, a0);
				mx = n0 * ab_0.lo;
				// tmp = mx * p + tmp
				cc_mul_add_wide_limb_limb(&mp_0, mx, p[0], ab_0.lo);
				hi_mp = mp_0.hi;
				t[0] = mp_0.lo;
				hi_ab = ab_0.hi;
#endif
			}
		}

#pragma unroll
		for (uint64_t i = 1; i < n; i++) {
			wide_limb_t  ab_i;
			wide_limb_t  mp_i;
			auto ai = a[i];
			auto pi = p[i];
			auto ti = t[i];

			if constexpr (has_carry)
			{
#if defined(PTX_MUL_IMPL_CUDA) && 1
				asm volatile(
					"mul.lo.u64  %0, %7, %8;\n\t"  // ab_lo = b_j * ai
					"mul.hi.u64  %1, %7, %8;\n\t"  // ab_hi = umul64hi(b_j, ai)
					"mul.lo.u64  %2, %9, %10;\n\t" // mp_lo = mx * pi
					"mul.hi.u64  %3, %9, %10;\n\t" // mp_hi = umul64hi(mx, pi)
					"add.cc.u64  %0, %0, %4;\n\t"  // ab_lo += hi_ab
					"addc.u64    %1, %1, 0;\n\t"   // ab_hi += carry
					"add.cc.u64  %0, %0, %5;\n\t"  // ab_lo += ti
					"addc.u64    %1, %1, 0;\n\t"   // ab_hi += carry
					"add.cc.u64  %2, %2, %6;\n\t"  // mp_lo += hi_mp
					"addc.u64    %3, %3, 0;\n\t"   // mp_hi += carry
					"add.cc.u64  %2, %2, %0;\n\t"  // mp_lo += ab_lo (updated)
					"addc.u64    %3, %3, 0;\n\t"   // mp_hi += carry
					: "=l"(ab_i.lo), // 0
					"=l"(ab_i.hi), // 1
					"=l"(mp_i.lo), // 2
					"=l"(mp_i.hi)  // 3
					: "l"(hi_ab),    // 4
					"l"(ti),       // 5
					"l"(hi_mp),    // 6
					"l"(b_j),      // 7
					"l"(ai),       // 8
					"l"(mx),       // 9
					"l"(pi)        // 10
					);
#else


				ab_i.lo = b_j * ai;
				ab_i.hi = umul64hi(b_j, ai);
				mp_i.lo = mx * pi;
				mp_i.hi = umul64hi(mx, pi);

				cc_add_wide_wide_limb(&ab_i, ab_i, hi_ab);
				cc_add_wide_wide_limb(&ab_i, ab_i, ti);
				cc_add_wide_wide_limb(&mp_i, mp_i, hi_mp);
				cc_add_wide_wide_limb(&mp_i, mp_i, ab_i.lo);
#endif

			}
			else
			{
#if defined(PTX_MUL_IMPL_CUDA) && 1
				asm volatile(
					"mul.lo.u64  %0, %6, %7;\n\t"  // ab_lo = b_j * ai
					"mul.hi.u64  %1, %6, %7;\n\t"  // ab_hi = umul64hi(b_j, ai)
					"mul.lo.u64  %2, %8, %9;\n\t"  // mp_lo = mx * pi
					"mul.hi.u64  %3, %8, %9;\n\t"  // mp_hi = umul64hi(mx, pi)
					"add.cc.u64  %0, %0, %4;\n\t"  // ab_i.lo += hi_ab
					"addc.u64    %1, %1, 0;\n\t"   // ab_i.hi += carry
					"add.cc.u64  %2, %2, %5;\n\t"  // mp_i.lo += hi_mp
					"addc.u64    %3, %3, 0;\n\t"   // mp_i.hi += carry
					"add.cc.u64  %2, %2, %0;\n\t"  // mp_i.lo += ab_i.lo (now updated)
					"addc.u64    %3, %3, 0;\n\t"   // mp_i.hi += carry
					: "=l"(ab_i.lo), // 0
					"=l"(ab_i.hi), // 1
					"=l"(mp_i.lo), // 2
					"=l"(mp_i.hi)	 // 3
					: "l"(hi_ab),	 // 4
					"l"(hi_mp)	 // 5
					"l"(b_j),      // 6
					"l"(ai),       // 7
					"l"(mx),       // 8
					"l"(pi)        // 9
					);
#else							


				ab_i.lo = b_j * ai;
				ab_i.hi = umul64hi(b_j, ai);
				mp_i.lo = mx * pi;
				mp_i.hi = umul64hi(mx, pi);
				cc_add_wide_wide_limb(&ab_i, ab_i, hi_ab);
				cc_add_wide_wide_limb(&mp_i, mp_i, hi_mp);
				cc_add_wide_wide_limb(&mp_i, mp_i, ab_i.lo);
#endif
			}
			hi_ab = ab_i.hi;
			t[i - 1] = mp_i.lo;
			hi_mp = mp_i.hi;
		}
		t[n] = hi_ab;

		{
			wide_limb_t wide1;

			if constexpr (has_carry) {
#if defined(PTX_MUL_IMPL_CUDA) && 1
				auto tn = t[n];

				uint64_t s0;       // tn + carry_in
				uint64_t c0;       // carry from s0
				uint64_t sum_lo;   // wide1.lo
				uint64_t sum_hi;   // wide1.hi

				asm volatile(
					// s0 = tn + carry_in; c0 = carry
					"add.cc.u64  %0, %4, %5;\n\t"
					"addc.u64    %1, 0, 0;\n\t"

					// (sum_lo,sum_hi) = hi_mp + s0 + c0
					"add.cc.u64  %2, %6, %0;\n\t"
					"addc.cc.u64 %2, %2, %1;\n\t"
					"addc.u64    %3, 0, 0;\n\t"
					: "=&l"(s0),         // %0  output
					"=&l"(c0),         // %1  output
					"=&l"(sum_lo),     // %2  output
					"=&l"(sum_hi)      // %3  output
					: "l"(tn),           // %4  input
					"l"(carry_in),     // %5  input
					"l"(hi_mp)         // %6  input
					);

				t[n] = s0;
				wide1.lo = sum_lo;
				wide1.hi = sum_hi;
#else
				wide_limb_t wide;
				cc_add_wide_limb_limb(&wide, t[n], carry_in);
				t[n] = wide.lo;
				cc_add_carry_wide_limb_limb(&wide1, hi_mp, t[n], wide.hi);
#endif
			}
			else {
				cc_add_wide_limb_limb(&wide1, hi_mp, t[n]);
			}

			t[n - 1] = wide1.lo;
			return wide1.hi;
		}
	}




	template<typename Field>
	BLST_HD BLST_FORCEINLINE
		void ptx_mul_mont_u64(
			uint64_t*  r,
			const uint64_t*  a,
			const uint64_t*  b)
	{
		using namespace blst;

		constexpr auto p = Field::modulus().data;
		constexpr int n = sizeof(Field) / sizeof(uint64_t);
		static_assert(n > 0 && (n % 2 == 0), "n must be positive and even");

		uint64_t aa[n];
#pragma unroll
		for (int i = 0; i < n; i++) {
			aa[i] = a[i];
		}
		uint64_t bb[n];
#pragma unroll
		for (int i = 0; i < n; i++) {
			bb[i] = b[i];
		}

		uint64_t t[n + 1] = { 0 };

		uint64_t carry = montgomery_fios_step<Field, n>(t, aa, bb[0]);

		// Columns j = 1..n-1 (fused)
		for (int j = 1; j < n; ++j) {
			carry = montgomery_fios_step<Field, n>(t, aa, bb[j], carry);
		}

		// Final conditional subtract (same as yours)
		uint64_t tmp[n];
		{
			uint64_t borrow = 0;
#pragma unroll
			for (int i = 0; i < n; i++) {
				wide_limb_t w;
				cc_sub_borrow_wide_limb_limb(&w, t[i], p[i], borrow);
				tmp[i] = w.lo;
				borrow = w.hi & 1;
			}
			uint64_t mask = carry - borrow;
#pragma unroll
			for (int i = 0; i < n; i++) {
				r[i] = (tmp[i] & ~mask) | (t[i] & mask);
			}
		}

//#ifndef NDEBUG
		uint64_t tt[n];
		constexpr uint64_t n0 = Field::montgomery_constant();
		blst::cc_mul_mont_n<n>(tt, aa, bb, Field::modulus().data, n0);
		for (size_t i = 0; i < n; i++) {
			if (r[i] != tt[i])
			{
				printf("Mismatch at limb %zu: r=%llu, expected=%llu\n", i, r[i], tt[i]);
				assert(false && "ptx_mul_mont_u64 result mismatch");
			}

		}
//#endif

	}

	template<typename Field>
	BLST_HD BLST_FORCEINLINE
		void ptx_sqr_mont_u64(
			uint64_t*  r,
			const uint64_t*  a)
	{
		using namespace blst;

		constexpr auto p = Field::modulus().data;
		constexpr int n = sizeof(Field) / sizeof(uint64_t);
		static_assert(n > 0 && (n % 2 == 0), "n must be positive and even");

		uint64_t aa[n];
#pragma unroll
		for (int i = 0; i < n; i++) {
			aa[i] = a[i];
		}
		uint64_t t[n + 1] = { 0 };

		// 
		uint64_t carry = montgomery_fios_step<Field, n>(t, aa, aa[0]);

		// Columns j = 1..n-1 (fused)
		for (int j = 1; j < n; ++j) {
			carry = montgomery_fios_step<Field, n>(t, aa, aa[j], carry);
		}

		// Final conditional subtract (same as yours)
		uint64_t tmp[n];
		{
			uint64_t borrow = 0;
#pragma unroll
			for (int i = 0; i < n; i++) {
				wide_limb_t w;
				cc_sub_borrow_wide_limb_limb(&w, t[i], p[i], borrow);
				tmp[i] = w.lo;
				borrow = w.hi & 1;
			}
			uint64_t mask = carry - borrow;
#pragma unroll
			for (int i = 0; i < n; i++) {
				r[i] = (tmp[i] & ~mask) | (t[i] & mask);
			}
		}
	}
}