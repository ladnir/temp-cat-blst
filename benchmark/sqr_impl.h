#pragma once
#include <benchmark/benchmark.h>
#include "blst/fields/Fp12.h"
//#include "blst/backend/cuda/mul_mont_cuda.h"
#include <thrust/device_vector.h>
#include "blst/curve/point.h"
#include <span>
#include "parallel_prefix_bench.h"
#include <chrono>

//BENCHMARK_MAIN();
//void FpMain();
//void Fp12Main();
using namespace blst;



struct SoAFp {
	//std::span<limb_t> data[NLIMBS(384)];
	limb_t* data[NLIMBS(384)];
	BLST_HD BLST_FORCEINLINE auto&& operator[](size_t idx) {
		BLST_ASSERT(idx < NLIMBS(384));
		return data[idx];
	}
	BLST_HD BLST_FORCEINLINE auto&& operator[](size_t idx) const {
		BLST_ASSERT(idx < NLIMBS(384));
		return data[idx];
	}

};

struct SoAFp2 {
	SoAFp data[2];
	BLST_HD BLST_FORCEINLINE auto&& operator[](size_t idx) {
		BLST_ASSERT(idx < 2);
		return data[idx];
	}
	BLST_HD BLST_FORCEINLINE auto&& operator[](size_t idx) const {
		BLST_ASSERT(idx < 2);
		return data[idx];
	}
};

struct SoAFp4 {
	SoAFp2 data[2];
	BLST_HD BLST_FORCEINLINE auto&& operator[](size_t idx) {
		BLST_ASSERT(idx < 2);
		return data[idx];
	}
	BLST_HD BLST_FORCEINLINE auto&& operator[](size_t idx) const {
		BLST_ASSERT(idx < 2);
		return data[idx];
	}
};
struct SoAFp6 {
	SoAFp2 data[3];
	BLST_HD BLST_FORCEINLINE auto&& operator[](size_t idx) {
		BLST_ASSERT(idx < 3);
		return data[idx];
	}
	BLST_HD BLST_FORCEINLINE auto&& operator[](size_t idx) const {
		BLST_ASSERT(idx < 3);
		return data[idx];
	}
};

struct SoAFp12 {
	SoAFp6 data[2];

	BLST_HD BLST_FORCEINLINE auto&& operator[](size_t idx) {
		BLST_ASSERT(idx < 2);
		return data[idx];
	}
	BLST_HD BLST_FORCEINLINE auto&& operator[](size_t idx) const {
		BLST_ASSERT(idx < 2);
		return data[idx];
	}


};

struct NoopIdx
{
	BLST_HD BLST_FORCEINLINE auto&& get(auto&& x)
	{
		static_assert(!std::is_pointer_v<decltype(x)>, "NoopIdx used with pointer");
		return x;
	}
};

// is_span_v
template<typename T>
struct is_span : std::false_type {};
template<typename U>
struct is_span<std::span<U>> : std::true_type {};
template<typename T>
inline constexpr bool is_span_v = is_span<std::remove_cv_t<T>>::value;

struct BasicIdx
{
	size_t idx;
	BLST_HD BLST_FORCEINLINE limb_t& get(limb_t& x) { return x; }
	BLST_HD BLST_FORCEINLINE const limb_t& get(const limb_t& x) { return x; }
	BLST_HD BLST_FORCEINLINE limb_t& get(std::span<limb_t> x) {
		BLST_ASSERT(idx < x.size());
		return x[idx];
	}
	BLST_HD BLST_FORCEINLINE const limb_t& get(std::span<const limb_t> x) {
		BLST_ASSERT(idx < x.size());
		return x[idx];
	}
	BLST_HD BLST_FORCEINLINE limb_t& get(limb_t* x) { return x[idx]; }
	BLST_HD BLST_FORCEINLINE const limb_t& get(const limb_t* x) const { return x[idx]; }

};


	// parallel prefix add + conditional reduction for SoA field elements
	template<typename Field, typename SoAField>
	BLST_HD __forceinline__
		void Fp__ppAdd(
			SoAField& r_,
			const SoAField& a_,
			const SoAField& b_,
			size_t index,
			size_t lane)
	{
#ifdef __CUDA_ARCH__
		// n = 6 limbs (384-bit). If you also need n=4, mirror this structure with group size 4.
		constexpr size_t n = 6;

		// Compute per-warp lane and 6-wide subgroup base and mask.
		const unsigned warp_lane = threadIdx.x & 31;         // my lane in the warp [0..31]
		const unsigned group_base = (warp_lane / n) * n;     // base lane of my 6-lane group
		const unsigned last_in_group = group_base + (n - 1); // last lane in my 6-lane group

		// Build a mask that covers exactly our 6-lane subgroup inside the warp.
		// For lanes [group_base, group_base+5], we set bits to 1.
		const unsigned group_mask = ((unsigned)0x3Fu) << group_base;

		// Sanity: only proceed for the 5 groups in a warp (0..29). If you ever map threads 30/31 here,
		// they’ll be masked out by group_mask and just compute harmlessly.
		// Load inputs from SoA
		const uint64_t ai = a_.data[lane][index];
		const uint64_t bi = b_.data[lane][index];
		const uint64_t pi = Field::modulus().data[lane];

		// ---------------------------
		// 1) Multi-limb add: x = a + b using limb-level Kogge–Stone
		// ---------------------------
		const uint64_t sum0 = ai + bi;

		// Limb-level generate/propagate for carry:
		// carry_{i+1} = g_i | (p_i & carry_i) with carry_0 = 0
		// g_i = (a_i + b_i) overflows -> (sum0 < a_i)
		// p_i = adding 1 to sum0 causes overflow -> (sum0 == UINT64_MAX)
		uint32_t g = (sum0 < ai) ? 1u : 0u;
		uint32_t p = (sum0 == ~0ull) ? 1u : 0u;

		// Kogge–Stone prefix across 6 lanes
		// d = 1,2,4 steps are enough for 6 lanes
		uint32_t g_up;
		uint32_t p_up;

		g_up = __shfl_up_sync(group_mask, g, 1);
		p_up = __shfl_up_sync(group_mask, p, 1);
		if (lane >= 1)
			g = g | (p & g_up);
		if (lane >= 1)
			p = p & p_up;

		g_up = __shfl_up_sync(group_mask, g, 2);
		p_up = __shfl_up_sync(group_mask, p, 2);
		if (lane >= 2)
			g = g | (p & g_up);
		if (lane >= 2)
			p = p & p_up;

		g_up = __shfl_up_sync(group_mask, g, 4);
		p_up = __shfl_up_sync(group_mask, p, 4);
		if (lane >= 4)
			g = g | (p & g_up);
		if (lane >= 4)
			p = p & p_up;

		// Carry-in for this limb is the combined g from the previous limb (or 0 for limb 0)
		uint32_t c_in = 0u;
		uint32_t g_prev = __shfl_up_sync(group_mask, g, 1);
		if (lane >= 1)
			c_in = g_prev;

		const uint64_t x = sum0 + static_cast<uint64_t>(c_in);

		// Final carry out of the 6-limb add (same for all lanes in the group)
		const uint32_t carry1 = __shfl_sync(group_mask, g, last_in_group);

		// ---------------------------
		// 2) Multi-limb subtract: s = x - p to prepare conditional reduction
		// We also need final borrow to decide if x >= p.
		// ---------------------------
		const uint64_t diff0 = x - pi;

		// Limb-level generate/propagate for borrow:
		// borrow_{i+1} = gb_i | (pb_i & borrow_i) with borrow_0 = 0
		// gb_i = (x_i < p_i)
		// pb_i = (diff0 == 0)   // subtracting 1 will borrow iff result is 0
		uint32_t gb = (x < pi) ? 1u : 0u;
		uint32_t pb = (diff0 == 0ull) ? 1u : 0u;

		// Kogge–Stone prefix for borrow
		uint32_t gb_up;
		uint32_t pb_up;

		gb_up = __shfl_up_sync(group_mask, gb, 1);
		pb_up = __shfl_up_sync(group_mask, pb, 1);
		if (lane >= 1)
			gb = gb | (pb & gb_up);
		if (lane >= 1)
			pb = pb & pb_up;

		gb_up = __shfl_up_sync(group_mask, gb, 2);
		pb_up = __shfl_up_sync(group_mask, pb, 2);
		if (lane >= 2)
			gb = gb | (pb & gb_up);
		if (lane >= 2)
			pb = pb & pb_up;

		gb_up = __shfl_up_sync(group_mask, gb, 4);
		pb_up = __shfl_up_sync(group_mask, pb, 4);
		if (lane >= 4)
			gb = gb | (pb & gb_up);
		if (lane >= 4)
			pb = pb & pb_up;

		// Borrow-in for this limb is the combined gb from the previous limb (or 0 for limb 0)
		uint32_t b_in = 0u;
		uint32_t gb_prev = __shfl_up_sync(group_mask, gb, 1);
		if (lane >= 1)
			b_in = gb_prev;

		const uint64_t s = diff0 - static_cast<uint64_t>(b_in);

		// Final borrow out (same for all lanes in group)
		const uint32_t borrow_out = __shfl_sync(group_mask, gb, last_in_group);

		// ---------------------------
		// 3) Decide reduction and blend
		// need_reduce = (carry1 != 0) || (x >= p)
		// Note: x >= p  <=>  no final borrow when subtracting p (borrow_out == 0)
		const uint32_t need_reduce = (carry1 != 0u) | ((borrow_out == 0u) ? 1u : 0u);

		// Branchless select
		const uint64_t ri = (need_reduce != 0u) ? s : x;

		// Store back to SoA
		r_.data[lane][index] = ri;
#else
		(void)r_;
		(void)a_;
		(void)b_;
		(void)index;
		(void)lane;
#endif
	}


	template<typename Field, typename Idx = NoopIdx>
	BLST_HD __forceinline__
		void Fp__add(
			auto& r,
			const auto& a,
			const auto& b,
			Idx idx)
	{
		constexpr auto n = sizeof(Field) / sizeof(limb_t);
		constexpr auto* BLST_restrict p = Field::modulus().data;

#ifdef __CUDA_ARCH__


		if constexpr (n == 6)
		{
			const limb_t& a0 = idx.get(a[0]);
			const limb_t& a1 = idx.get(a[1]);
			const limb_t& a2 = idx.get(a[2]);
			const limb_t& a3 = idx.get(a[3]);
			const limb_t& a4 = idx.get(a[4]);
			const limb_t& a5 = idx.get(a[5]);

			const limb_t& b0 = idx.get(b[0]);
			const limb_t& b1 = idx.get(b[1]);
			const limb_t& b2 = idx.get(b[2]);
			const limb_t& b3 = idx.get(b[3]);
			const limb_t& b4 = idx.get(b[4]);
			const limb_t& b5 = idx.get(b[5]);


			limb_t x0;
			limb_t x1;
			limb_t x2;
			limb_t x3;
			limb_t x4;
			limb_t x5;

			limb_t s0;
			limb_t s1;
			limb_t s2;
			limb_t s3;
			limb_t s4;
			limb_t s5;

			limb_t t0;
			limb_t t1;
			limb_t t2;
			limb_t t3;
			limb_t t4;
			limb_t t5;


			// ===== a + b with carry -> x, carry1 in a GPR ===================================
			uint32_t carry1 = 0;
			asm volatile(
				"add.cc.u64   %0,  %7,  %13;\n\t"
				"addc.cc.u64  %1,  %8,  %14;\n\t"
				"addc.cc.u64  %2,  %9,  %15;\n\t"
				"addc.cc.u64  %3,  %10, %16;\n\t"
				"addc.cc.u64  %4,  %11, %17;\n\t"
				"addc.cc.u64  %5,  %12, %18;\n\t"   // final must be .cc
				"mov.u32      %6,  0;\n\t"
				"addc.u32     %6,  %6,  %6;\n\t"    // carry1 = CC
				: "=&l"(x0),
				"=&l"(x1),
				"=&l"(x2),
				"=&l"(x3),
				"=&l"(x4),
				"=&l"(x5),
				"+r"(carry1)
				: "l"(a0), "l"(a1), "l"(a2), "l"(a3), "l"(a4), "l"(a5),
				"l"(b0), "l"(b1), "l"(b2), "l"(b3), "l"(b4), "l"(b5)
				);

			t0 = x0;
			t1 = x1;
			t2 = x2;
			t3 = x3;
			t4 = x4;
			t5 = x5;

			// ===== x - p -> s (we’ll compute ge(x,p) separately) ============================
			uint32_t noborrow;
			asm volatile(
				// %0..%5  -> s0..s5 (64-bit outputs)
				// %6      -> noborrow (32-bit output)  **NOT** used by the sub lines
				// %7..%12 -> x0..x5 (64-bit inputs)
				// %13..%18-> p0..p5 (64-bit inputs)

				"sub.cc.u64   %0,  %7,  %13;\n\t"
				"subc.cc.u64  %1,  %8,  %14;\n\t"
				"subc.cc.u64  %2,  %9,  %15;\n\t"
				"subc.cc.u64  %3,  %10, %16;\n\t"
				"subc.cc.u64  %4,  %11, %17;\n\t"
				"subc.cc.u64  %5,  %12, %18;\n\t"   // CC now = !borrow
				"mov.u32      %6,  0;\n\t"
				"addc.u32     %6,  %6,  %6;\n\t"    // noborrow = CC (1=no borrow)
				: "=&l"(s0),
				"=&l"(s1),
				"=&l"(s2),
				"=&l"(s3),
				"=&l"(s4),
				"=&l"(s5),
				"=&r"(noborrow)                   // %6 — not referenced above
				: "l"(x0), "l"(x1), "l"(x2), "l"(x3), "l"(x4), "l"(x5),
				"l"(p[0]), "l"(p[1]), "l"(p[2]), "l"(p[3]), "l"(p[4]), "l"(p[5])
				);

			// ===== need_reduce = (carry1 != 0) || (x >= p) ==================================
			uint32_t need_reduce;
			asm volatile("or.b32 %0, %1, %2;" : "=&r"(need_reduce) : "r"(carry1), "r"(noborrow));


			//// ===== blend: r = need_reduce ? s : t ==========================================
			limb_t o0;
			limb_t o1;
			limb_t o2;
			limb_t o3;
			limb_t o4;
			limb_t o5;

			asm volatile(
				"{\n\t"
				".reg .pred psel6;\n\t"
				"setp.ne.u32 psel6, %12, 0;\n\t"
				"selp.u64 %0, %6,  %13, psel6;\n\t"
				"selp.u64 %1, %7,  %14, psel6;\n\t"
				"selp.u64 %2, %8,  %15, psel6;\n\t"
				"selp.u64 %3, %9,  %16, psel6;\n\t"
				"selp.u64 %4, %10, %17, psel6;\n\t"
				"selp.u64 %5, %11, %18, psel6;\n\t"
				"}\n\t"
				: "=&l"(o0),
				"=&l"(o1),
				"=&l"(o2),
				"=&l"(o3),
				"=&l"(o4),
				"=&l"(o5)
				: "l"(s0), "l"(s1), "l"(s2), "l"(s3), "l"(s4), "l"(s5),
				"r"(need_reduce),
				"l"(t0), "l"(t1), "l"(t2), "l"(t3), "l"(t4), "l"(t5)
				);

			idx.get(r[0]) = o0;
			idx.get(r[1]) = o1;
			idx.get(r[2]) = o2;
			idx.get(r[3]) = o3;
			idx.get(r[4]) = o4;
			idx.get(r[5]) = o5;

		}
		else if constexpr (n == 4)
		{

			auto& a0 = idx.get(a[0]);
			auto& a1 = idx.get(a[1]);
			auto& a2 = idx.get(a[2]);
			auto& a3 = idx.get(a[3]);

			auto& b0 = idx.get(b[0]);
			auto& b1 = idx.get(b[1]);
			auto& b2 = idx.get(b[2]);
			auto& b3 = idx.get(b[3]);

			limb_t x0;
			limb_t x1;
			limb_t x2;
			limb_t x3;

			limb_t s0;
			limb_t s1;
			limb_t s2;
			limb_t s3;

			limb_t t0;
			limb_t t1;
			limb_t t2;
			limb_t t3;


			// ===== a + b with carry -> x, carry1 in a GPR ===================================
			uint32_t carry1 = 0;
			asm volatile(
				"add.cc.u64   %0,  %5,  %9;\n\t"
				"addc.cc.u64  %1,  %6,  %10;\n\t"
				"addc.cc.u64  %2,  %7,  %11;\n\t"
				"addc.cc.u64  %3,  %8,  %12;\n\t"  // last limb must be .cc to update CC
				"mov.u32      %4,  0;\n\t"
				"addc.u32     %4,  %4,  %4;\n\t"   // carry1 = CC
				: "=&l"(x0),
				"=&l"(x1),
				"=&l"(x2),
				"=&l"(x3),
				"+r"(carry1)
				: "l"(a0), "l"(a1), "l"(a2), "l"(a3),
				"l"(b0), "l"(b1), "l"(b2), "l"(b3)
				);

			t0 = x0;
			t1 = x1;
			t2 = x2;
			t3 = x3;

			// ===== x - p -> s (we’ll compute ge(x,p) separately) ============================
			uint32_t noborrow;
			asm volatile(
				// %0..%3  -> s0..s3 (64-bit outputs)
				// %4      -> noborrow (32-bit output)  **not used** in sub lines
				// %5..%8  -> x0..x3 (64-bit inputs)
				// %9..%12 -> p0..p3 (64-bit inputs)

				"sub.cc.u64   %0, %5,  %9;\n\t"
				"subc.cc.u64  %1, %6,  %10;\n\t"
				"subc.cc.u64  %2, %7,  %11;\n\t"
				"subc.cc.u64  %3, %8,  %12;\n\t"   // final must be .cc
				"mov.u32      %4, 0;\n\t"
				"addc.u32     %4, %4, %4;\n\t"     // noborrow = CC (1 = no borrow)
				: "=&l"(s0),
				"=&l"(s1),
				"=&l"(s2),
				"=&l"(s3),
				"=&r"(noborrow)
				: "l"(x0), "l"(x1), "l"(x2), "l"(x3),
				"l"(p[0]), "l"(p[1]), "l"(p[2]), "l"(p[3])
				);


			// ===== need_reduce = (carry1 != 0) || (x >= p) ==================================
			uint32_t need_reduce;
			asm volatile("or.b32 %0, %1, %2;" : "=&r"(need_reduce) : "r"(carry1), "r"(noborrow));


			// ===== blend: r = need_reduce ? s : t ==========================================
			limb_t o0;
			limb_t o1;
			limb_t o2;
			limb_t o3;

			asm volatile(
				".reg .pred psel4;\n\t"
				"setp.ne.u32 psel4, %8, 0;\n\t"
				"selp.u64 %0, %4,  %9,  psel4;\n\t"
				"selp.u64 %1, %5,  %10, psel4;\n\t"
				"selp.u64 %2, %6,  %11, psel4;\n\t"
				"selp.u64 %3, %7,  %12, psel4;\n\t"
				: "=&l"(o0),
				"=&l"(o1),
				"=&l"(o2),
				"=&l"(o3)
				: "l"(s0), "l"(s1), "l"(s2), "l"(s3),
				"r"(need_reduce),
				"l"(t0), "l"(t1), "l"(t2), "l"(t3)
				);

			idx.get(r[0]) = o0;
			idx.get(r[1]) = o1;
			idx.get(r[2]) = o2;
			idx.get(r[3]) = o3;
		}
		else
		{
			// static assert to trigger a compile-time error for unsupported Field sizes
			static_assert(n == 6 || n == 4, "ptx_add_mod_u64 only supports 256 or 384-bit fields");
		}
#else
		Fp rr, aa, bb;
		for (u64 i = 0; i < n; ++i)
		{
			aa[i] = idx.get(a[i]);
			bb[i] = idx.get(b[i]);
		}

		cc_add_mod_n<n>(rr.data, aa.data, bb.data, p);
		for (u64 i = 0; i < n; ++i)
		{
			idx.get(r[i]) = rr[i];
		}
#endif // __CUDA_ARCH__

	}
	template<typename Field, typename Idx = NoopIdx >
	BLST_HD BLST_FORCEINLINE
		void Fp__sub(
			auto& r,
			const auto& a,
			const auto& b,
			Idx idx)
	{
		constexpr int n = sizeof(Field) / sizeof(limb_t);
		constexpr auto* BLST_restrict p = Field::modulus().data;


#ifdef __CUDA_ARCH__
		if constexpr (n == 6)
		{
			auto& a0 = idx.get(a[0]);
			auto& a1 = idx.get(a[1]);
			auto& a2 = idx.get(a[2]);
			auto& a3 = idx.get(a[3]);
			auto& a4 = idx.get(a[4]);
			auto& a5 = idx.get(a[5]);

			auto& b0 = idx.get(b[0]);
			auto& b1 = idx.get(b[1]);
			auto& b2 = idx.get(b[2]);
			auto& b3 = idx.get(b[3]);
			auto& b4 = idx.get(b[4]);
			auto& b5 = idx.get(b[5]);

			limb_t r0;
			limb_t r1;
			limb_t r2;
			limb_t r3;
			limb_t r4;
			limb_t r5;

			// ---- Pass 1: r = a - b ; capture final borrow (via CC -> NOT-borrow) ----
			// After the last subc, CC = !borrow. We extract CC with addc.u32 tmp,0,0 and then XOR 1.
			uint32_t borrow;

			// outputs: r0=%0 r1=%1 r2=%2 r3=%3 r4=%4 r5=%5 borrow=%6
			// inputs:  a0=%7 a1=%8 a2=%9 a3=%10 a4=%11 a5=%12  b0=%13 b1=%14 b2=%15 b3=%16 b4=%17 b5=%18
			asm volatile(
				"sub.cc.u64      %0,  %7,  %13;\n\t"
				"subc.cc.u64     %1,  %8,  %14;\n\t"
				"subc.cc.u64     %2,  %9,  %15;\n\t"
				"subc.cc.u64     %3,  %10, %16;\n\t"
				"subc.cc.u64     %4,  %11, %17;\n\t"
				"subc.cc.u64     %5,  %12, %18;\n\t"
				"addc.u32        %6,  0,   0;\n\t"   // %6 = !borrow
				"xor.b32         %6,  %6,  1;\n\t"   // %6 = borrow
				: "=&l"(r0),
				"=&l"(r1),
				"=&l"(r2),
				"=&l"(r3),
				"=&l"(r4),
				"=&l"(r5),
				"=&r"(borrow)
				: "l"(a0),
				"l"(a1),
				"l"(a2),
				"l"(a3),
				"l"(a4),
				"l"(a5),
				"l"(b0),
				"l"(b1),
				"l"(b2),
				"l"(b3),
				"l"(b4),
				"l"(b5)
				);

			//limb_t mask_, carry_, borrow_;
			//limb_t ret[n];
			//{
			//	size_t i;
			//	for (borrow_ = 0, i = 0; i < n; i++) {
			//		wide_limb_t wide;
			//		cc_sub_borrow_wide_limb_limb(&wide, a[i], b[i], borrow_);
			//		ret[i] = wide.lo;
			//		borrow_ = wide.hi & 1;
			//	}
			//	BLST_ASSERT(borrow_ == borrow);
			//	BLST_ASSERT(ret[0] == r0);
			//	BLST_ASSERT(ret[1] == r1);
			//	BLST_ASSERT(ret[2] == r2);
			//	BLST_ASSERT(ret[3] == r3);
			//	BLST_ASSERT(ret[4] == r4);
			//	BLST_ASSERT(ret[5] == r5);
			//}

			// ---- Pass 2: if (borrow) r += p ; predicated, carry-chained; final add need not set CC ----
			asm volatile(
				"{\n\t"
				".reg .pred pb;\n\t"
				"setp.ne.u32     pb, %12, 0;\n\t"

				"@pb add.cc.u64  %0,  %0,  %6;\n\t"
				"@pb addc.cc.u64 %1,  %1,  %7;\n\t"
				"@pb addc.cc.u64 %2,  %2,  %8;\n\t"
				"@pb addc.cc.u64 %3,  %3,  %9;\n\t"
				"@pb addc.cc.u64 %4,  %4,  %10;\n\t"
				"@pb addc.u64     %5,  %5,  %11;\n\t"
				"}\n\t"
				: "+l"(r0),
				"+l"(r1),
				"+l"(r2),
				"+l"(r3),
				"+l"(r4),
				"+l"(r5)
				: "l"(p[0]),
				"l"(p[1]),
				"l"(p[2]),
				"l"(p[3]),
				"l"(p[4]),
				"l"(p[5]),
				"r"(borrow)
				);

			//{
			//	mask_ = 0 - borrow_;
			//	BLST_launder(mask_);
			//	limb_t i;
			//	for (carry_ = 0, i = 0; i < n; i++) {
			//		wide_limb_t wide;
			//		cc_add_carry_wide_limb_limb(&wide, ret[i], p[i] & mask_, carry_);
			//		ret[i] = wide.lo;
			//		carry_ = wide.hi;
			//	}

			//	BLST_ASSERT(borrow_ == borrow);
			//	BLST_ASSERT(ret[0] == r0);
			//	BLST_ASSERT(ret[1] == r1);
			//	BLST_ASSERT(ret[2] == r2);
			//	BLST_ASSERT(ret[3] == r3);
			//	BLST_ASSERT(ret[4] == r4);
			//	BLST_ASSERT(ret[5] == r5);
			//}

			idx.get(r[0]) = r0;
			idx.get(r[1]) = r1;
			idx.get(r[2]) = r2;
			idx.get(r[3]) = r3;
			idx.get(r[4]) = r4;
			idx.get(r[5]) = r5;
		}
		else if constexpr (n == 4)
		{
			auto& a0 = idx.get(a[0]);
			auto& a1 = idx.get(a[1]);
			auto& a2 = idx.get(a[2]);
			auto& a3 = idx.get(a[3]);

			auto& b0 = idx.get(b[0]);
			auto& b1 = idx.get(b[1]);
			auto& b2 = idx.get(b[2]);
			auto& b3 = idx.get(b[3]);

			limb_t r0;
			limb_t r1;
			limb_t r2;
			limb_t r3;

			uint32_t borrow;

			// ---- Pass 1: r = a - b ; get final borrow from CC ----
			asm volatile(
				"sub.cc.u64      %0,  %5,  %9;\n\t"
				"subc.cc.u64     %1,  %6,  %10;\n\t"
				"subc.cc.u64     %2,  %7,  %11;\n\t"
				"subc.cc.u64     %3,  %8,  %12;\n\t"
				"addc.u32        %4,  0,   0;\n\t"   // %4 = !borrow
				"xor.b32         %4,  %4,  1;\n\t"   // %4 = borrow
				: "=&l"(r0),
				"=&l"(r1),
				"=&l"(r2),
				"=&l"(r3),
				"=&r"(borrow)
				: "l"(a0),
				"l"(a1),
				"l"(a2),
				"l"(a3),
				"l"(b0),
				"l"(b1),
				"l"(b2),
				"l"(b3)
				);

			// ---- Pass 2: if (borrow) r += p ; predicated, carry-chained ----
			asm volatile(
				".reg .pred pb;\n\t"
				"setp.ne.u32     pb, %8, 0;\n\t"

				"@pb add.cc.u64  %0,  %0,  %4;\n\t"
				"@pb addc.cc.u64 %1,  %1,  %5;\n\t"
				"@pb addc.cc.u64 %2,  %2,  %6;\n\t"
				"@pb addc.u64     %3,  %3,  %7;\n\t"
				: "+l"(r0),
				"+l"(r1),
				"+l"(r2),
				"+l"(r3)
				: "l"(p[0]),
				"l"(p[1]),
				"l"(p[2]),
				"l"(p[3]),
				"r"(borrow)
				);

			idx.get(r[0]) = r0;
			idx.get(r[1]) = r1;
			idx.get(r[2]) = r2;
			idx.get(r[3]) = r3;
		}
		else
		{
			// static assert to trigger a compile-time error for unsupported Field sizes
			static_assert(n == 6 || n == 4, "ptx_sub_mod_u64 only supports 256 or 384-bit fields");
		}
#else
		Fp rr, aa, bb;
		for (u64 i = 0; i < n; ++i)
		{
			aa[i] = idx.get(a[i]);
			bb[i] = idx.get(b[i]);
		}

		cc_sub_mod_n<n>(rr.data, aa.data, bb.data, p);
		for (u64 i = 0; i < n; ++i)
		{
			idx.get(r[i]) = rr[i];
		}
#endif
	}

#ifdef __CUDA_ARCH__
#define PTX_MUL_IMPL_CUDA
#endif

//	BLST_HD BLST_FORCEINLINE
//		limb_t umul64hi(limb_t a, limb_t b) {
//#ifdef __CUDA_ARCH__
//		return __umul64hi(a, b);
//#else
//		limb_t hi;
//		static_assert(sizeof(unsigned long long) == sizeof(limb_t), "limb_t size mismatch");
//		_mulx_u64(a, b, (unsigned long long*) & hi);
//		return hi;
//#endif
//	}


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
	template<typename Field, typename Carry = std::monostate>
	BLST_HD BLST_FORCEINLINE
		limb_t montgomery_fios_step(
			limb_t* t,
			limb_t* a,
			limb_t b_j,
			Carry carry_in = {})
	{
		using namespace blst;

		constexpr bool has_carry = !std::is_same<Carry, std::monostate>::value;
		constexpr auto p = Field::modulus().data;
		constexpr limb_t n0 = Field::montgomery_constant();
		constexpr size_t n = sizeof(Field) / sizeof(limb_t);

		// tmp = b[0] * a
		limb_t hi_ab;
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
		for (limb_t i = 1; i < n; i++) {
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

				limb_t s0;       // tn + carry_in
				limb_t c0;       // carry from s0
				limb_t sum_lo;   // wide1.lo
				limb_t sum_hi;   // wide1.hi

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




	template<typename Field, typename Idx = NoopIdx>
	BLST_HD BLST_FORCEINLINE
		void ptx_mul_mont_u64(
			auto& r,
			const auto& a,
			const auto& b,
			Idx idx)
	{
		using namespace blst;

		static_assert(
			std::is_same_v<decltype(a), const SoAFp&> ||
			std::is_same_v<decltype(a), const Field&>, "a must be of type Field");

		constexpr auto p = Field::modulus().data;
		constexpr int n = sizeof(Field) / sizeof(limb_t);
		static_assert(n > 0 && (n % 2 == 0), "n must be positive and even");

		limb_t aa[n];
#pragma unroll
		for (int i = 0; i < n; i++) {
			aa[i] = idx.get(a[i]);
		}
		limb_t bb[n];
#pragma unroll
		for (int i = 0; i < n; i++) {
			bb[i] = idx.get(b[i]);
		}

		limb_t t[n + 1] = { 0 };

		limb_t carry = montgomery_fios_step<Field>(t, aa, bb[0]);

		// Columns j = 1..n-1 (fused)
		for (int j = 1; j < n; ++j) {
			carry = montgomery_fios_step<Field>(t, aa, bb[j], carry);
		}

		// Final conditional subtract (same as yours)
		limb_t tmp[n];
		{
			limb_t borrow = 0;
#pragma unroll
			for (int i = 0; i < n; i++) {
				wide_limb_t w;
				cc_sub_borrow_wide_limb_limb(&w, t[i], p[i], borrow);
				tmp[i] = w.lo;
				borrow = w.hi & 1;
			}
			limb_t mask = carry - borrow;
#pragma unroll
			for (int i = 0; i < n; i++) {
				idx.get(r[i]) = (tmp[i] & ~mask) | (t[i] & mask);
			}
		}

		//#ifndef NDEBUG
				//limb_t tt[n];
				//constexpr limb_t n0 = Field::montgomery_constant();
				//blst::cc_mul_mont_n<n>(tt, aa, bb, Field::modulus().data, n0);
				//for (size_t i = 0; i < n; i++) {
				//	if (r[i] != tt[i])
				//	{
				//		printf("Mismatch at limb %zu: r=%llu, expected=%llu\n", i, r[i], tt[i]);
				//		assert(false && "ptx_mul_mont_u64 result mismatch");
				//	}

				//}
		//#endif

	}

	struct Fp4 {
		Fp2 data[2];
		BLST_HD Fp2& operator[](size_t idx) {
			return data[idx];
		}
		BLST_HD const Fp2& operator[](size_t idx)const {
			return data[idx];
		}
	};


	BLST_HD __forceinline__ void Fp2__sqr(Fp2& ret, const auto& a, auto&& idx)
	{
		Fp t0, t1;

		Fp__add<Fp>(t0, a[0], a[1], idx);
		Fp__sub<Fp>(t1, a[0], a[1], idx);

		ptx_mul_mont_u64<Fp>(ret[1], a[0], a[1], idx);
		Fp__add<Fp>(ret[1], ret[1], ret[1], idx);

		ptx_mul_mont_u64<Fp>(ret[0], t0, t1, idx);
	}

	BLST_HD __forceinline__ void Fp2__add(auto& ret, const auto& a, const auto& b, auto&& idx)
	{
		Fp__add<Fp>(ret[0], a[0], b[0], idx);
		Fp__add<Fp>(ret[1], a[1], b[1], idx);
	}


	BLST_HD __forceinline__ void Fp2__sub(auto& ret, const auto& a, const auto& b, auto&& idx)
	{
		Fp__sub<Fp>(ret[0], a[0], b[0], idx);
		Fp__sub<Fp>(ret[1], a[1], b[1], idx);
	}

	BLST_HD __forceinline__  void Fp2__mul_by_u_plus_1(Fp2& ret, const Fp2& a)
	{
		Fp t;

		Fp__add<Fp>(t, a[0], a[1], NoopIdx{});
		Fp__sub<Fp>(ret[0], a[0], a[1], NoopIdx{});
		vec_copy(ret[1].data, t.data, sizeof(t));
	}

#if 1

	BLST_HD __forceinline__ void Fp2__sqr(Fp4& ret, const auto /*Fp2*/& a0, const auto /*Fp2*/& a1, auto&& idx)
	{
		Fp2__add(ret[1], a0, a1, idx);
		Fp2__sqr(ret[1], ret[1], idx);

		{
			Fp2 t1;
			Fp2__sqr(t1, a1, idx);
			Fp2__mul_by_u_plus_1(ret[0], t1);
			Fp2__sub(ret[1], ret[1], t1, idx);
		}
		{
			Fp2 t0;
			Fp2__sqr(t0, a0, idx);
			Fp2__add(ret[0], ret[0], t0, idx);
			Fp2__sub(ret[1], ret[1], t0, idx);
		}
	}

	BLST_HD __forceinline__ void sqr_fp4(Fp4& ret, const auto/*Fp2*/& a0, const auto/*Fp2*/& a1, auto&& idx)
	{
		Fp2__add(ret[1], a0, a1, idx);
		Fp2__sqr(ret[1], ret[1], idx);

		{
			//Fp2 t1;
			Fp2__sqr(ret[0], a1, idx);
			Fp2__sub(ret[1], ret[1], ret[0], idx);
			Fp2__mul_by_u_plus_1(ret[0], ret[0]);

		}
		{

			Fp2 t0;
			Fp2__sqr(t0, a0, idx);
			Fp2__add(ret[0], ret[0], t0, idx);
			Fp2__sub(ret[1], ret[1], t0, idx);
		}
	}


	BLST_HD __forceinline__ void cyclotomic_sqr(auto/*Fp12*/& ret, const auto/*Fp12*/& a, auto&& idx)
	{

		{
			Fp4 t0;
			sqr_fp4(t0, a[0][0], a[1][1], idx);

			// r00 = 4 (a00 a11)[0] - 2 a00
			Fp2__sub(ret[0][0], t0[0], a[0][0], idx);      // r00 = t0_0 - a00
			Fp2__add(ret[0][0], ret[0][0], t0[0], idx);    // r00 += t0_0
			Fp2__add(ret[0][0], ret[0][0], ret[0][0], idx);// r00 *= 2 

			// r11 = 4 (a00 a11)[1] + 2 a11
			Fp2__add(ret[1][1], t0[1], a[1][1], idx);
			Fp2__add(ret[1][1], ret[1][1], t0[1], idx);
			Fp2__add(ret[1][1], ret[1][1], ret[1][1], idx);
		}
		{
			Fp4 t1;
			sqr_fp4(t1, a[1][0], a[0][2], idx);

			Fp2__sub(ret[0][1], t1[0], a[0][1], idx);
			Fp2__add(ret[0][1], ret[0][1], t1[0], idx);
			Fp2__add(ret[0][1], ret[0][1], ret[0][1], idx);

			Fp2__add(ret[1][2], t1[1], a[1][2], idx);
			Fp2__add(ret[1][2], ret[1][2], t1[1], idx);
			Fp2__add(ret[1][2], ret[1][2], ret[1][2], idx);

		}
		{
			Fp4 t2;
			sqr_fp4(t2, a[0][1], a[1][2], idx);

			Fp2__sub(ret[0][2], t2[0], a[0][2], idx);
			Fp2__add(ret[0][2], ret[0][2], ret[0][2], idx);
			Fp2__add(ret[0][2], ret[0][2], t2[0], idx);

			Fp2__mul_by_u_plus_1(t2[1], t2[1]);
			Fp2__add(ret[1][0], t2[1], a[1][0], idx);
			Fp2__add(ret[1][0], ret[1][0], ret[1][0], idx);
			Fp2__add(ret[1][0], ret[1][0], t2[1], idx);
		}
	}
#else

	BLST_HD __forceinline__ void Fp2__sqr(Fp4& ret, const Fp2& a0, const Fp2& a1)
	{
		Fp2 t0, t1;

		Fp2__sqr(t0, a0);
		Fp2__sqr(t1, a1);
		Fp2__add(ret[1], a0, a1);

		Fp2__mul_by_u_plus_1(ret[0], t1);
		Fp2__add(ret[0], ret[0], t0);

		Fp2__sqr(ret[1], ret[1]);
		Fp2__sub(ret[1], ret[1], t0);
		Fp2__sub(ret[1], ret[1], t1);
	}

	BLST_HD __forceinline__ void sqr_fp4(Fp4& ret, const Fp2& a0, const Fp2& a1)
	{
		Fp2 t0, t1;

		Fp2__sqr(t0, a0);
		Fp2__sqr(t1, a1);
		Fp2__add(ret[1], a0, a1);

		Fp2__mul_by_u_plus_1(ret[0], t1);
		Fp2__add(ret[0], ret[0], t0);

		Fp2__sqr(ret[1], ret[1]);
		Fp2__sub(ret[1], ret[1], t0);
		Fp2__sub(ret[1], ret[1], t1);
	}

	BLST_HD __forceinline__ void cyclotomic_sqr(Fp12& ret, const Fp12& a)
	{
		Fp4 t0, t1, t2;

		sqr_fp4(t0, a[0][0], a[1][1]);
		sqr_fp4(t1, a[1][0], a[0][2]);
		sqr_fp4(t2, a[0][1], a[1][2]);

		Fp2__sub(ret[0][0], t0[0], a[0][0]);
		Fp2__add(ret[0][0], ret[0][0], ret[0][0]);
		Fp2__add(ret[0][0], ret[0][0], t0[0]);

		Fp2__sub(ret[0][1], t1[0], a[0][1]);
		Fp2__add(ret[0][1], ret[0][1], ret[0][1]);
		Fp2__add(ret[0][1], ret[0][1], t1[0]);

		Fp2__sub(ret[0][2], t2[0], a[0][2]);
		Fp2__add(ret[0][2], ret[0][2], ret[0][2]);
		Fp2__add(ret[0][2], ret[0][2], t2[0]);

		Fp2__mul_by_u_plus_1(t2[1], t2[1]);
		Fp2__add(ret[1][0], t2[1], a[1][0]);
		Fp2__add(ret[1][0], ret[1][0], ret[1][0]);
		Fp2__add(ret[1][0], ret[1][0], t2[1]);

		Fp2__add(ret[1][1], t0[1], a[1][1]);
		Fp2__add(ret[1][1], ret[1][1], ret[1][1]);
		Fp2__add(ret[1][1], ret[1][1], t0[1]);

		Fp2__add(ret[1][2], t1[1], a[1][2]);
		Fp2__add(ret[1][2], ret[1][2], ret[1][2]);
		Fp2__add(ret[1][2], ret[1][2], t1[1]);
	}
#endif
