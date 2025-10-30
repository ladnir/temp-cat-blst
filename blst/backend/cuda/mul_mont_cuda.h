
#include <stdint.h>
#include "blst/common/defines.h"
#include "blst/backend/backend_portable.h"
#include <variant>
//
// ===== Montgomery multiply (n x 64 CIOS, no __int128) =====
// r = a * b * R^{-1} mod p, base = 2^64
template<typename Field>
BLST_DEVICE BLST_FORCEINLINE
void ptx_mul_mont_u64_old(
	uint64_t* BLST_restrict r,
	const uint64_t* BLST_restrict a,
	const uint64_t* BLST_restrict b
)
{
	constexpr auto BLST_restrict p = Field::modulus().data;
	constexpr uint64_t n0 = Field::montgomery_constant();
	constexpr int n = sizeof(Field) / sizeof(uint64_t);
	static_assert(sizeof(*p) == sizeof(uint64_t));

	uint64_t carry;
	uint64_t t[n + 1];
	int i;
	for (i = 0; i < n + 1; i++)
		t[i] = 0ull;

	int j;
	for (j = 0; j < n; j++)
	{
		uint64_t bj;
		bj = b[j];

		// ---- t += a * bj  (schoolbook column) ----
		carry = 0ull;
		int k;
#pragma unroll n
		for (k = 0; k < n; k++)
		{
			uint64_t lo;
			uint64_t plo;
			uint64_t phi;
			uint64_t sum;
			uint64_t c1;

			uint64_t ak;
			ak = a[k];

			// current acc = (carry | t[k])
			// p = ak * bj
			// 

			lo = t[k];

			// p = ak * bj
			plo = ak * bj;
			phi = __umul64hi(ak, bj);
			auto oldcarry = carry;

			sum = lo + plo;
			c1 = (sum < lo) ? 1ull : 0ull;

			t[k] = sum;
			carry = carry + phi + c1;

			{
				using namespace blst;
				wide_limb_t wide;
				//wide.lo = _mulx_u64(bj, ak, &wide.hi); 
				wide.lo = plo + oldcarry;
				auto loCarry = (wide.lo < plo) ? 1ull : 0ull;
				wide.hi = phi + loCarry;


				if (wide.lo != sum || wide.hi != carry) {
					printf("Mismatch at j=%d, k=%d: wide=(%llu,%llu), computed=(%llu,%llu)\n",
						j, k, wide.lo, wide.hi, sum, carry);
				}
			}

			if (j == 0 && k < 2)
				printf("a[%d]=%llu, b[%d]=%llu, c=%llu,\n\t ak*bj+c=plo=%llu, phi=%llu\n", k, ak, j, bj, oldcarry, sum, carry);
		}

		// propagate carry into t[n]
		{
			uint64_t old;
			uint64_t sum;
			uint64_t c;

			old = t[n];
			sum = old + carry;
			c = (sum < old) ? 1ull : 0ull;

			t[n] = sum;
			// c would be an extra limb; CIOS keeps t as 5 limbs, so we drop c.
		}

		//printf("After processing b[%d]: carry=%llu, t=", j, carry);
  //      for (int idx = 0; idx < n + 1; idx++) {
  //          printf("%llu ", t[idx]);
  //      }
  //      printf("\n");

		// ---- Montgomery step: m = t[0] * n0 (mod 2^64) ----
		uint64_t m;
		m = t[0] * n0;

		// ---- t += m * p, with the implicit right-shift (drop t[0]) ----
		carry = 0ull;
#pragma unroll n
		for (k = 0; k < n; k++)
		{
			uint64_t lo;
			uint64_t plo;
			uint64_t phi;
			uint64_t sum;
			uint64_t c1;

			lo = t[k];
			plo = p[k] * m;
			phi = __umul64hi(p[k], m);

			sum = lo + plo;
			c1 = (sum < lo) ? 1ull : 0ull;

			if (k != 0)
				t[k - 1] = sum;

			carry = carry + phi + c1;
		}

		// write the most significant kept limb and update t[n]
		{
			uint64_t old;
			uint64_t sum;
			uint64_t c;

			old = t[n];
			sum = old + carry;
			c = (sum < old) ? 1ull : 0ull;

			t[n - 1] = sum;
			t[n] = c;
		}

		//printf("after Montgomery step: carry=%llu, t=", carry);
  //      for (int idx = 0; idx < n + 1; idx++) {
  //          printf("%llu ", t[idx]);
  //      }
  //      printf("\n");
	}

	// ---- Final conditional subtract: diff = t - p (branchless per-limb) ----
	uint64_t diff[n];
	uint64_t borrow;

	borrow = 0ull;
#pragma unroll n
	for (i = 0; i < n; i++)
	{
		uint64_t ti;
		uint64_t pi;
		uint64_t spi;
		uint64_t ri;

		ti = t[i];
		pi = p[i];
		spi = pi + borrow;
		ri = ti - spi;

		borrow = (ti < spi) ? 1ull : 0ull;
		diff[i] = ri;
	}

	// Choose t only if (carry==0 && borrow==1), else choose diff
	uint64_t choose_t = (borrow & (carry ^ 1ull));   // 1 iff (borrow==1 && carry==0)
	uint64_t mask = 0ull - choose_t;             // ~0 or 0

#pragma unroll
	for (int i = 0; i < n; i++) {
		r[i] = (diff[i] & ~mask) | (t[i] & mask);
	}
}


// ===== Helper functions for Montgomery multiply =====
// r[0..n] = a[0..n-1] * b + carry
// 
// if there is a carry in, its added to the last limb as part of 
// the montgomery update step. The new carry out is returned.
template<int n, typename Carry = std::monostate>
BLST_HD BLST_FORCEINLINE
Carry mul_add_chain(
	uint64_t* BLST_restrict r,
	const uint64_t* BLST_restrict a,
	uint64_t b,
	Carry carry = {})
{
	using namespace blst;

	// the first add chain does not have carry input
	constexpr bool has_carry = !std::is_same<Carry, std::monostate>::value;

	// tmp = b[0] * a
	uint64_t hi = 0;
	{
		wide_limb_t  wide;
		if constexpr (has_carry)
		{
			cc_mul_add_wide_limb_limb(&wide, b, a[0], r[0]);
		}
		else
		{
			cc_mul_wide_limb_limb(&wide, b, a[0]);
		}
		r[0] = wide.lo;
		hi = wide.hi;
	}
#pragma unroll
	for (uint64_t i = 1; i < n; i++) {
		wide_limb_t  wide;

		if constexpr (has_carry)
		{
			cc_mul_add_add_wide_limb_limb(&wide, b, a[i], hi, r[i]);
		}
		else
		{
			cc_mul_add_wide_limb_limb(&wide, b, a[i], hi);
		}
		r[i] = wide.lo;
		hi = wide.hi;
	}
	r[n] = hi;

	if constexpr (has_carry) {
		wide_limb_t wide;

		cc_add_wide_limb_limb(&wide, r[n], carry);
		r[n] = wide.lo;
		return wide.hi;
	}
	else {
		return {};
	}
}

template<typename Field, int n, typename Carry = std::monostate>
BLST_HD BLST_FORCEINLINE
uint64_t mont_update_chain(
	uint64_t* BLST_restrict r,
	Carry carry = {})
{
	using namespace blst;

	constexpr auto p = Field::modulus().data;
	constexpr uint64_t n0 = Field::montgomery_constant();
	limb_t mx = n0 * r[0];

	// tmp = mx * p + tmp
	wide_limb_t  wide;
	cc_mul_add_wide_limb_limb(&wide, mx, p[0], r[0]);
	limb_t hi = wide.hi;
	r[0] = wide.lo;
#pragma unroll
	for (uint64_t i = 1; i < n; i++) {
		cc_mul_add_add_wide_limb_limb(&wide, mx, p[i], hi, r[i]);
		r[i - 1] = wide.lo;
		hi = wide.hi;
	}
	if constexpr (!std::is_same<Carry, std::monostate>::value) {
		cc_add_carry_wide_limb_limb(&wide, hi, r[n], carry);
	}
	else {
		cc_add_wide_limb_limb(&wide, hi, r[n]);
	}
	// final limb
	//cc_add_carry_wide_limb_limb(&wide, hi, r[n], carry);
	r[n - 1] = wide.lo;
	return wide.hi;
}


// ===== Montgomery multiply (n x 64 CIOS, no __int128) =====
// r = a * b * R^{-1} mod p, base = 2^64
template<typename Field>
BLST_HD BLST_FORCEINLINE
void ptx_mul_mont_u64_(
	uint64_t* BLST_restrict r,
	const uint64_t* BLST_restrict a,
	const uint64_t* BLST_restrict b)
{

	constexpr auto p = Field::modulus().data;
	constexpr int n = sizeof(Field) / sizeof(uint64_t);

	static_assert(n > 0 && n % 2 == 0, "n must be positive and even");
	using namespace blst;

	//limb_t mask, borrow, mx, hi, carry;
	limb_t  tmp[n + 1];

	// tmp[0...n] = a[0..n-1] * b[0]
	mul_add_chain<n>(tmp, a, b[0]);

	limb_t carry = mont_update_chain<Field, n>(tmp);

	for (size_t j = 1; j < n; ++j) {
		//wide_limb_t  wide;

		// tmp += a * b[j] + carry
		carry = mul_add_chain<n>(tmp, a, b[j], carry);

		//cc_add_wide_limb_limb(&wide, tmp[n], carry);
		//tmp[n] = wide.lo;
		//carry = wide.hi;

		carry = mont_update_chain<Field, n>(tmp, carry);
	}

	// conditional subtract: tmp - p and shift
	{
		limb_t borrow = 0;
		for (int i = 0; i < n; i++) {
			// wide = tmp - (p + borrow)
			wide_limb_t wide;
			cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
			r[i] = wide.lo;
			borrow = wide.hi & 1;
		}

		limb_t mask = carry - borrow;

		// select between tmp and ret
		for (int i = 0; i < n; i++)
			r[i] = (r[i] & ~mask) | (tmp[i] & mask);
	}

}












// Computes one fused limb step (FIOS-style).
// t_im1_out  <- v.lo
// acc_hi_io  <- v.hi   (updated in place)
// Inputs are single 64-bit limbs in registers.
//
// Requirements:
//   - Compile for device (sm_50+ for native 64-bit int ops).
//   - Call from __device__/__host__ __device__ code paths guarded as needed.
//
__device__ __forceinline__
void mont_hot_step_ptx(
	uint64_t& t_im1_out,
	uint64_t& acc_hi_io,
	const uint64_t a_i,
	const uint64_t b_j,
	const uint64_t t_i,
	const uint64_t m,
	const uint64_t p_i)
{
	uint64_t u_lo;
	uint64_t u_hi;
	uint64_t v_lo;
	uint64_t v_hi;

	asm volatile
		(
			"{\n\t"
			// ---- Step A: u = a_i * b_j + t_i + acc_hi_io ----
			"mul.lo.u64 %0, %6, %7;\n\t"
			"mul.hi.u64 %1, %6, %7;\n\t"
			"add.cc.u64 %0, %0, %8;\n\t"
			"addc.u64   %1, %1, %9;\n\t"

			// ---- Step B: v = m * p_i + u_lo + u_hi ----
			"mul.lo.u64 %2, %10, %11;\n\t"
			"mul.hi.u64 %3, %10, %11;\n\t"
			"add.cc.u64 %2, %2, %0;\n\t"
			"addc.u64   %3, %3, %1;\n\t"
			"}\n\t"
			: "=&l"(u_lo),
			"=&l"(u_hi),
			"=&l"(v_lo),
			"=&l"(v_hi)
			: "l"(t_im1_out),   // (dummy to keep matching slots; not used)
			"l"(acc_hi_io),   // (dummy)
			"l"(a_i),         // %6
			"l"(b_j),         // %7
			"l"(t_i),         // %8
			"l"(acc_hi_io),   // %9  (old acc_hi)
			"l"(m),           // %10
			"l"(p_i)          // %11
			);

	t_im1_out = v_lo;
	acc_hi_io = v_hi;
}



// Process one Montgomery column (FIOS-style):
// Input:
//   t[0..n]   : running accumulator (n+1 limbs)
//   a[0..n-1] : multiplicand
//   b_j       : current limb of b
//   carry_in  : carry that you were adding into t[n] in your CIOS version
// Output:
//   t[0..n] updated in place (with t[n-1] finalized and shifted)
// Returns:
//   carry_out to feed into next column j+1
//template<typename Field, int n>
//BLST_HD BLST_FORCEINLINE
//uint64_t montgomery_fios_step(
//	uint64_t* BLST_restrict t,
//	const uint64_t* BLST_restrict a,
//	uint64_t b_j,
//	uint64_t carry_in)
//{
//	using namespace blst;
//
//	constexpr auto p = Field::modulus().data;
//	constexpr uint64_t n0 = Field::montgomery_constant();
//
//	// ---- i = 0: compute new t0, then m = t0 * n0 ----
//	wide_limb_t w0;
//	// u = a0*bj + t0
//	cc_mul_add_wide_limb_limb(&w0, b_j, a[0], t[0]);
//	uint64_t acc_lo = w0.lo;   // u_lo
//	uint64_t acc_hi = w0.hi;   // u_hi
//
//	// m = u_lo * n0  (low limb only)
//	uint64_t m = acc_lo * n0;
//
//	// v0 = m*p0 + u_lo + u_hi
//	cc_mul_add_add_wide_limb_limb(&w0, m, p[0], acc_hi, acc_lo);
//	// carry into i=1:
//	acc_hi = w0.hi;
//
//	// ---- i = 1..n-1: fused steps ----
//#pragma unroll
//	for (int i = 1; i < n; i++) {
//#ifdef __CUDA_ARCH__
//		mont_hot_step_ptx(
//			t[i - 1],
//			acc_hi,
//			a[i],
//			b_j,
//			t[i],
//			m,
//			p[i]);
//#else
//		// Step A
//		cc_mul_add_add_wide_limb_limb(&w0, b_j, a[i], acc_hi, t[i]);
//		uint64_t sum_lo = w0.lo;
//		acc_hi = w0.hi;
//		// Step B
//		wide_limb_t w1;
//		cc_mul_add_add_wide_limb_limb(&w1, m, p[i], acc_hi, sum_lo);
//		t[i - 1] = w1.lo;
//		acc_hi = w1.hi;
//#endif
//	}
//
//	// finalize with t[n] and carry_in
//	cc_add_carry_wide_limb_limb(&w0, acc_hi, t[n], carry_in);
//	t[n - 1] = w0.lo;
//	return w0.hi;
//}

BLST_HD BLST_FORCEINLINE
uint64_t umul64hi(uint64_t a, uint64_t b) {
#ifdef __CUDA_ARCH__
	return __umul64hi(a, b);
#else
	uint64_t hi;
	_mulx_u64(a, b, &hi);
	return hi;
#endif
}

template<typename Field, int n, typename Carry = std::monostate>
BLST_HD BLST_FORCEINLINE
uint64_t montgomery_fios_step(
	uint64_t* BLST_restrict t,
	const uint64_t* BLST_restrict a,
	uint64_t b_j,
	Carry carry_in = {})
{
	using namespace blst;

	constexpr bool has_carry = !std::is_same<Carry, std::monostate>::value;
	constexpr auto p = Field::modulus().data;
	constexpr uint64_t n0 = Field::montgomery_constant();

	//Carry c = mul_add_chain<n>(t, a, b_j, carry_in);
	Carry c = {};

	// tmp = b[0] * a
	uint64_t hi_ab;
	limb_t mx, hi_mp;
	{

		wide_limb_t  wide;
		if constexpr (has_carry)
		{
			cc_mul_add_wide_limb_limb(&wide, b_j, a[0], t[0]);
		}
		else
		{
			cc_mul_wide_limb_limb(&wide, b_j, a[0]);
		}
		//t[0] = wide.lo;
		hi_ab = wide.hi;


		mx = n0 * wide.lo;
		// tmp = mx * p + tmp
		wide_limb_t  wide1;
		cc_mul_add_wide_limb_limb(&wide1, mx, p[0], wide.lo);
		hi_mp = wide1.hi;
		t[0] = wide1.lo;
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
#ifdef __CUDA_ARCH__
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
#ifdef __CUDA_ARCH__
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
			wide_limb_t wide;

			cc_add_wide_limb_limb(&wide, t[n], carry_in);
			t[n] = wide.lo;
			c = wide.hi;

			cc_add_carry_wide_limb_limb(&wide1, hi_mp, t[n], c);
		}
		else {
			cc_add_wide_limb_limb(&wide1, hi_mp, t[n]);
		}
		// final limb
		//cc_add_carry_wide_limb_limb(&wide, hi, r[n], carry);
		t[n - 1] = wide1.lo;
		return wide1.hi;
	}
}




template<typename Field>
BLST_HD BLST_FORCEINLINE
void ptx_mul_mont_u64(
	uint64_t* BLST_restrict r,
	const uint64_t* BLST_restrict a,
	const uint64_t* BLST_restrict b)
{
	using namespace blst;

	constexpr auto p = Field::modulus().data;
	constexpr int n = sizeof(Field) / sizeof(uint64_t);
	static_assert(n > 0 && (n % 2 == 0), "n must be positive and even");

	uint64_t t[n + 1] = { 0 };

	//// Seed t with the first column (pure a*b0) so we can switch into fused steps
	//{
	//	wide_limb_t w;
	//	auto b0 = b[0];
	//	cc_mul_wide_limb_limb(&w, b0, a[0]);   
	//	t[0] = w.lo; 
	//	uint64_t hi = w.hi;
	//	for (int i = 1; i < n; i++) {
	//		cc_mul_add_wide_limb_limb(&w, b0, a[i], hi);
	//		t[i] = w.lo;
	//		hi = w.hi;
	//	}
	//	t[n] = hi;
	//}


	uint64_t carry = montgomery_fios_step<Field, n>(t, a, b[0]);

	// Columns j = 1..n-1 (fused)
	for (int j = 1; j < n; ++j) {
		carry = montgomery_fios_step<Field, n>(t, a, b[j], carry);
	}

	// Final conditional subtract (same as yours)
	uint64_t tmp[n];
	{
		uint64_t borrow = 0;
		for (int i = 0; i < n; i++) {
			wide_limb_t w;
			cc_sub_borrow_wide_limb_limb(&w, t[i], p[i], borrow);
			tmp[i] = w.lo;
			borrow = w.hi & 1;
		}
		uint64_t mask = carry - borrow;
		for (int i = 0; i < n; i++) {
			r[i] = (tmp[i] & ~mask) | (t[i] & mask);
		}
	}
	//// conditional subtract: tmp - p and shift
	//{
	//	limb_t borrow = 0;
	//	for (int i = 0; i < n; i++) {
	//		// wide = tmp - (p + borrow)
	//		wide_limb_t wide;
	//		cc_sub_borrow_wide_limb_limb(&wide, t[i], p[i], borrow);
	//		r[i] = wide.lo;
	//		borrow = wide.hi & 1;
	//	}

	//	limb_t mask = carry - borrow;

	//	// select between tmp and ret
	//	for (int i = 0; i < n; i++)
	//		r[i] = (r[i] & ~mask) | (t[i] & mask);
	//}
}