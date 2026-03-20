#pragma once

#include "blst/common/defines.h"
#include "blst/fields/Fp.h"



#include <array>
#include <cstdint>
#include <cstdio>
#include <type_traits>
#include <stdint.h>

#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include "blst/backend/cuda/CSA_cuda.h"
namespace cg = cooperative_groups;


namespace blst
{


	__device__ __forceinline__ auto& get_coop(limb_t* in, u64 rank)
	{
		return in[rank];
	}
	__device__ __forceinline__ auto& get_coop(const limb_t* in, u64 rank)
	{
		return in[rank];
	}

	__device__ __forceinline__ auto& get_coop(limb_t& in, u64 rank)
	{
		return in;
	}
	__device__ __forceinline__ auto& get_coop(const limb_t& in, u64 rank)
	{
		return in;
	}

	// ===== helpers =====

	// ===== helpers (same semantics/layout as your non-coop) =====
	BLST_HD BLST_FORCEINLINE
		void add_props(limb_t x, limb_t y, uint32_t& g, uint32_t& p, limb_t& s0)
	{
		s0 = x + y;
		g = (s0 < x) ? 1u : 0u;
		p = (s0 == ~0ull) ? 1u : 0u;
	}

	BLST_HD BLST_FORCEINLINE
		uint32_t pack_gp(uint32_t gA, uint32_t pA, uint32_t gB, uint32_t pB)
	{
		BLST_ASSERT((gA | pA | gB | pB) <= 1u);
		uint32_t g = gA | (gB << 1);
		uint32_t p = pA | (pB << 1);
		return (p << 8) | g;
	}


	//inline void add_props(limb_t x, limb_t y, uint32_t& g, uint32_t& p, limb_t& s0)
	//{
	//	s0 = x + y;
	//	g = (s0 < x) ? 1u : 0u;
	//	p = (s0 == ~0ull) ? 1u : 0u;
	//}

	//inline uint32_t pack_gp(uint32_t gA, uint32_t pA, uint32_t gB, uint32_t pB)
	//{
	//	BLST_ASSERT((gA & ~1u) == 0);
	//	BLST_ASSERT((pA & ~1u) == 0);
	//	BLST_ASSERT((gB & ~1u) == 0);
	//	BLST_ASSERT((pB & ~1u) == 0);

	//	uint32_t w = 0;
	//	w |= gA << 0;
	//	w |= gB << 1;
	//	w |= pA << 8;
	//	w |= pB << 9;
	//	return w;
	//}

	//inline void unpack_gp(uint32_t w, uint32_t& gA, uint32_t& pA, uint32_t& gB, uint32_t& pB)
	//{
	//	gA = (w >> 0) & 1u;
	//	gB = (w >> 1) & 1u;
	//	pA = (w >> 8) & 1u;
	//	pB = (w >> 9) & 1u;
	//}

	inline void shuffle_down(
		uint32_t out[],
		const uint32_t in[],
		int n,
		int shift)
	{
		for (int i = 0; i < n; i++)
		{
			if (i + shift < n)
				out[i] = in[i + shift];
			else
				out[i] = 0u;
		}
	}

	inline void shuffle_up(
		uint32_t out[],
		const uint32_t in[],
		int n,
		int shift)
	{
		for (int i = 0; i < n; i++)
		{
			if (i >= shift)
				out[i] = in[i - shift];
			else
				out[i] = 0u;
		}
	}

	// Kogge–Stone over NLIMBS limbs, packed (updates both adds together)
// Uses shuffle_up(...) to mimic the coop version's shfl_up.
	template<int NLIMBS>
	BLST_INLINE BLST_HD void ks_packed_prefix(uint32_t gp[NLIMBS])
	{
		uint32_t up[NLIMBS];

		int d = 1;
		while (d < NLIMBS)
		{
			shuffle_up(up, gp, NLIMBS, d);

			for (int i = 0; i < NLIMBS; i++)
			{
				auto g = gp[i];
				auto p = gp[i] >> 8;
				auto gu = up[i];
				auto pu = up[i] >> 8;

				g = g | (p & gu);
				p = p & pu;

				gp[i] = (g & 255) | (p << 8);
			}

			d = d << 1;
		}
	}

	// ===== main: single-thread packed dual-prefix add with correct (C<<1) =====
	template<typename Field>
	inline void parallelPrefix_add(
		limb_t* c,          // out: a + b
		const limb_t* a,
		const limb_t* b
	)
	{
		constexpr auto NLIMBS = Field::num_limb();
		constexpr auto pneg = Field::neg_modulus().data;

		limb_t C[NLIMBS];

		uint32_t gp[NLIMBS];

		for (int i = 0; i < NLIMBS; i++)
		{
			C[i] = (a[i] & b[i]) | (a[i] & pneg[i]) | (b[i] & pneg[i]);
			gp[i] = C[i] >> 63;
		}

		uint32_t CPrev[NLIMBS];
		shuffle_up(CPrev, gp, NLIMBS, 1);


		limb_t s0A[NLIMBS];
		limb_t s0B[NLIMBS];
		for (int i = 0; i < NLIMBS; i++)
		{

			limb_t S = a[i] ^ b[i] ^ pneg[i];
			//BLST_ASSERT(CPrev[i] == (i ? C[i - 1] >> 63 : 0));
			limb_t Csh = (C[i] << 1) | CPrev[i];

			uint32_t gA0;
			uint32_t pA0;
			uint32_t gB0;
			uint32_t pB0;

			add_props(a[i], b[i], gA0, pA0, s0A[i]);
			add_props(S, Csh, gB0, pB0, s0B[i]);
			gp[i] = pack_gp(gA0, pA0, gB0, pB0);
		}

		ks_packed_prefix<NLIMBS>(gp);
		u32 prev[NLIMBS];
		shuffle_up(prev, gp, NLIMBS, 1);

		limb_t ab[NLIMBS];
		limb_t* abp = c;

		for (int i = 0; i < NLIMBS; i++)
		{
			uint32_t p = prev[i];
			auto cinA = p & 1u;
			auto cinB = p >> 1;

			ab[i] = s0A[i] + (limb_t)cinA;
			abp[i] = s0B[i] + (limb_t)cinB;
		}

		limb_t carryA = (gp[NLIMBS - 1] >> 0) & 1u;
		limb_t carryB = (gp[NLIMBS - 1] >> 1) & 1u;

		auto mask = carryA - carryB;

		for (int i = 0; i < NLIMBS; i++)
		{
			c[i] = (abp[i] & mask) | (ab[i] & ~mask);
		}
	}

	// ===== Standard Parallel Prefix Add (Sum = S + C) =====
	// Uses the provided ks_packed_prefix helper.
	// No modular reduction, just pure 64-bit limb addition with carry propagation.
	template<int n>
	static inline void parallelPrefix_add(
		uint64_t* sum,
		const uint64_t* s,
		const uint64_t* c,
		uint64_t* carry_out_columns)
	{
		uint32_t gp[n];

		// 1. Generate Local G/P
		// The helper expects: g in low byte (mask 255), p in high byte (<< 8).
		for (int i = 0; i < n; i++)
		{
			uint64_t raw_sum = s[i] + c[i];

			// Generate: Carry generated if sum wraps (overflows 64-bit)
			uint32_t g = (raw_sum < s[i]) ? 1u : 0u;

			// Propagate: Carry propagates if sum is all 1s
			uint32_t p = (raw_sum == 0xFFFFFFFFFFFFFFFFull) ? 1u : 0u;

			gp[i] = g | (p << 8);
		}

		// 2. Run Warp-Level Prefix Scan
		// After this, gp[i] contains the cumulative carry-out state for lane i.
		ks_packed_prefix<n>(gp);

		// 3. Shift to get Carry-In
		// The carry *into* lane i is the carry *out* (the 'g' bit) of lane i-1.
		uint32_t prev_gp[n];
		// Note: Assuming shuffle_up handles uint32_t based on your helper usage
		shuffle_up(prev_gp, gp, n, 1);

		// 4. Final Summation
		for (int i = 0; i < n; i++)
		{
			// Lane 0 has no carry-in from "below" (it is 0).
			// For other lanes, extract the 'g' bit (bit 0) from the previous lane.
			uint64_t cin = (i == 0) ? 0ull : (uint64_t)(prev_gp[i] & 1u);

			sum[i] = s[i] + c[i] + cin;
		}

		// 5. Output Final Carry
		// The total carry out is the carry out of the MSB lane
		*carry_out_columns = (uint64_t)(gp[n - 1] & 1u);
	}

	// ===== Parallel Prefix Subtract (Diff = A - B) =====
	// Computes diff = a - b.
	// Returns the final borrow-out from the MSB lane.
	template<int n>
	static BLST_INLINE BLST_HD
		uint64_t parallelPrefix_sub(
			uint64_t* diff,
			const uint64_t* a,
			const uint64_t* b)
	{
		uint32_t gp[n];

		// 1. Generate Local Borrow/Propagate
		// A borrow is generated if a < b.
		// A borrow is propagated if a == b (incoming borrow causes outgoing borrow).
		for (int i = 0; i < n; i++)
		{
			uint32_t g = (a[i] < b[i]) ? 1u : 0u;
			uint32_t p = (a[i] == b[i]) ? 1u : 0u;
			gp[i] = g | (p << 8);
		}

		// 2. Run Warp-Level Prefix Scan (Kogge-Stone)
		// gp[i] now holds the cumulative borrow-out for lane i
		ks_packed_prefix<n>(gp);

		// 3. Shift to get Borrow-In
		uint32_t prev_gp[n];
		shuffle_up(prev_gp, gp, n, 1);

		// 4. Final Subtraction
		for (int i = 0; i < n; i++)
		{
			// Lane 0 has no borrow-in (0). 
			uint64_t borrow_in = (i == 0) ? 0ull : (uint64_t)(prev_gp[i] & 1u);

			// diff = a - b - borrow_in
			diff[i] = a[i] - b[i] - borrow_in;
		}

		// Return the borrow out of the MSB (the final global borrow)
		return (uint64_t)(gp[n - 1] & 1u);
	}

	template<int n>
	struct Tile
	{
		// bit mask for this tile's group, (n 1s) << base
		unsigned mMask;
		// base index within the group
		int mBase;
		// rank within the group, rank = threadIdx.x % n
		int mRank;
		// group id within the warp. gid = threadIdx.x / n
		int mGid;
		// the number of full groups in the warp
		static constexpr int mNumGroups = 32 / n;

		__device__ __forceinline__
			Tile()
		{
			int lane = threadIdx.x & 31;
			mGid = lane / n;
			mBase = mGid * n;
			mRank = lane - mBase;
			mMask = ((1 << n) - 1) << mBase;
		}

		__device__ __forceinline__ int size() const { return n; }
		__device__ __forceinline__ int thread_rank() const { return mRank; }
		__device__ __forceinline__ bool is_full() const { return mGid < mNumGroups; }

		static constexpr __device__ __forceinline__ int tiles_per_warp() { return mNumGroups; }

		template<typename T>
		__device__ __forceinline__ T shfl_up(T v, int delta) const
		{
			return __shfl_up_sync(mMask, v, delta);
		}

		template<typename T>
		__device__ __forceinline__ T shfl_down(T v, int delta) const {
			return __shfl_down_sync(mMask, v, delta);
		}

		template<typename T>
		__device__ __forceinline__ T shfl(T v, int src_rank) const
		{
			int src_abs = mBase + src_rank;
			return __shfl_sync(mMask, v, src_abs);
		}
	};



	// One-shuffle-per-stage packed Kogge–Stone over a tile
	template<int NLIMBS>
	BLST_DEVICE BLST_FORCEINLINE
		uint32_t ks_packed_prefix_tile(Tile<NLIMBS> tile, uint32_t pack)
	{
		static_assert(NLIMBS == 6 || NLIMBS == 4, "supported: 256/384 bits");

		uint32_t up;

		up = tile.shfl_up(pack, 1);
		if (tile.thread_rank() >= 1)
		{
			uint32_t G = pack & 0x3u;
			uint32_t P = (pack >> 8) & 0x3u;
			uint32_t Gu = up & 0x3u;
			uint32_t Pu = (up >> 8) & 0x3u;
			G = (G | (P & Gu)) & 0x3u;
			P = (P & Pu) & 0x3u;
			pack = (P << 8) | G;
		}

		up = tile.shfl_up(pack, 2);
		if (tile.thread_rank() >= 2)
		{
			uint32_t G = pack & 0x3u;
			uint32_t P = (pack >> 8) & 0x3u;
			uint32_t Gu = up & 0x3u;
			uint32_t Pu = (up >> 8) & 0x3u;
			G = (G | (P & Gu)) & 0x3u;
			P = (P & Pu) & 0x3u;
			pack = (P << 8) | G;
		}

		if constexpr (NLIMBS == 6)
		{
			up = tile.shfl_up(pack, 4);
			if (tile.thread_rank() >= 4)
			{
				uint32_t G = pack & 0x3u;
				uint32_t P = (pack >> 8) & 0x3u;
				uint32_t Gu = up & 0x3u;
				uint32_t Pu = (up >> 8) & 0x3u;
				G = (G | (P & Gu)) & 0x3u;
				P = (P & Pu) & 0x3u;
				pack = (P << 8) | G;
			}
		}

		return pack;
	}

	template<typename Field, int NLIMBS>
	BLST_DEVICE BLST_FORCEINLINE void pp_sub_coop_tile_aos(
		Tile<NLIMBS> tile,
		auto&& A, auto&& B, auto&& R)
	{
		const int rank = tile.thread_rank();
		// Resolve inputs (Register or Memory via get_coop)
		limb_t a = get_coop(A, rank);
		limb_t b = get_coop(B, rank);
		limb_t p_val = Field::modulus().data[rank];

		// 1. Sub: t = a - b
		limb_t t = a - b;
		uint32_t g = (a < b) ? 1u : 0u;     // Borrow generated
		uint32_t p = (a == b) ? 1u : 0u;    // Borrow propagates

		uint32_t pack = g | (p << 8);
		pack = ks_packed_prefix_tile<NLIMBS>(tile, pack);

		uint32_t prev = tile.shfl_up(pack, 1);
		limb_t bin = (rank == 0) ? 0 : (prev & 1u);
		t -= bin;

		limb_t total_borrow = tile.shfl(pack, NLIMBS - 1) & 1u;

		// 2. Add P (if borrow occurred) -> t_plus_p
		limb_t t_plus_p = t + p_val;
		g = (t_plus_p < t) ? 1u : 0u;
		p = (t_plus_p == ~0ull) ? 1u : 0u;

		pack = g | (p << 8);
		pack = ks_packed_prefix_tile<NLIMBS>(tile, pack);
		prev = tile.shfl_up(pack, 1);
		limb_t cin = (rank == 0) ? 0 : (prev & 1u);
		t_plus_p += cin;

		// 3. Select
		// If total_borrow is 1, result was negative, we need t_plus_p.
		// mask is all 1s if borrow=1.
		limb_t mask = 0ull - total_borrow;

		get_coop(R, rank) = (t_plus_p & mask) | (t & ~mask);
	}

	// One tile computes c = (a + b) mod p using p' = -p
	// AoS layout: arrays A,B,PNEG,C are contiguous by element: base + elem*NLIMBS + limb
	template<typename Field, int NLIMBS>
	BLST_DEVICE BLST_FORCEINLINE
		void pp_add_coop_tile_aos(
			Tile<NLIMBS> tile,
			auto&& A,
			auto&& B,
			auto&& C)
	{
		static_assert(NLIMBS == 6 || NLIMBS == 4, "supported: 256/384 bits");

		const int limb = tile.thread_rank();

		//const size_t off = elem_idx * (size_t)NLIMBS + (size_t)limb;

		limb_t a_i = get_coop(A, limb);
		limb_t b_i = get_coop(B, limb);
		limb_t pn = Field::neg_modulus().data[limb];

		limb_t Cmaj = (a_i & b_i) | (a_i & pn) | (b_i & pn);

		uint32_t C_msb = (uint32_t)(Cmaj >> 63);
		uint32_t Cprev = tile.shfl_up(C_msb, 1);
		if (limb == 0)
			Cprev = 0u;

		uint32_t gA;
		uint32_t pA;
		limb_t   s0A;
		add_props(a_i, b_i, gA, pA, s0A);

		limb_t S = a_i ^ b_i ^ pn;
		limb_t Csh = (Cmaj << 1) | (limb_t)Cprev;

		uint32_t gB;
		uint32_t pB;
		limb_t   s0B;
		add_props(S, Csh, gB, pB, s0B);

		uint32_t pack = pack_gp(gA, pA, gB, pB);

		pack = ks_packed_prefix_tile<NLIMBS>(tile, pack);

		uint32_t prevPack = tile.shfl_up(pack, 1);
		if (limb == 0)
			prevPack = 0u;

		uint32_t cinA = prevPack & 1u;
		uint32_t cinB = (prevPack >> 1) & 1u;

		limb_t x_i = s0A + (limb_t)cinA;
		limb_t y_i = s0B + (limb_t)cinB;

		uint32_t last_pack = tile.shfl(pack, NLIMBS - 1);

		uint32_t carryA = last_pack & 1u;
		uint32_t carryB = (last_pack >> 1) & 1u;

		limb_t mask = (limb_t)carryA - (limb_t)carryB;

		limb_t c_i = (y_i & mask) | (x_i & ~mask);

		get_coop(C, limb) = c_i;
	}

	// Kernel: AoS inputs, handles tail tiles safely.
	// For NLIMBS=6, each warp has 5 full tiles + 1 tail of size 2; we skip the tail.
	template<typename Field>
	BLST_GLOBAL void pp_add_coop_aos_kernel(
		const limb_t* BLST_restrict A,
		const limb_t* BLST_restrict B,
		limb_t* BLST_restrict C,
		size_t n_elems)
	{
		constexpr int NLIMBS = Field::num_limb();

		static_assert(NLIMBS == 6 || NLIMBS == 4, "supported: 256/384 bits");

		// thread's tile
		auto tile = Tile<NLIMBS>();

		// Skip tail tile (e.g., size 2 when NLIMBS=6)
		if (tile.is_full() == false)
			return;

		const int tiles_per_warp_full = tile.tiles_per_warp(); // 5 for NLIMBS=6, 8 for NLIMBS=4

		BLST_ASSERT(blockDim.x % 32 == 0);

		int warp_id_in_block = (threadIdx.x / 32);
		int warps_per_block = blockDim.x / 32;
		size_t warp_global_id = (size_t)blockIdx.x * (size_t)warps_per_block + (size_t)warp_id_in_block;
		size_t tile_idx_in_warp = tile.mGid; // 0..tiles_per_warp_full-1
		size_t elem_idx = warp_global_id * (size_t)tiles_per_warp_full + (size_t)tile_idx_in_warp;
		if (elem_idx >= n_elems)
			return;

		A += elem_idx * (size_t)NLIMBS;
		B += elem_idx * (size_t)NLIMBS;
		C += elem_idx * (size_t)NLIMBS;

		pp_add_coop_tile_aos<Field, NLIMBS>(tile, A, B, C);
	}

	template<typename Field>
	inline cudaError_t launch_pp_add_mod_aos(
		const limb_t* A,
		const limb_t* B,
		limb_t* C,
		size_t n_elems,
		int warps_per_block = 4,
		cudaStream_t stream = 0)
	{
		constexpr int NLIMBS = Field::num_limb();

		static_assert(NLIMBS == 6 || NLIMBS == 4, "NLIMBS must be 4 or 6");

		if (n_elems == 0)
			return cudaSuccess;

		if (A == nullptr)
			return cudaErrorInvalidDevicePointer;
		if (B == nullptr)
			return cudaErrorInvalidDevicePointer;
		if (C == nullptr)
			return cudaErrorInvalidDevicePointer;

		if (warps_per_block <= 0)
			warps_per_block = 4;

		const int threads_per_block = 32 * warps_per_block;

		const int tiles_per_warp_full = 32 / NLIMBS;  // 5 for NLIMBS=6, 8 for NLIMBS=4

		size_t warps_needed = (n_elems + tiles_per_warp_full - 1) / tiles_per_warp_full;

		dim3 grid;
		grid.x = static_cast<unsigned>((warps_needed + warps_per_block - 1) / warps_per_block);

		dim3 block;
		block.x = static_cast<unsigned>(threads_per_block);

		pp_add_coop_aos_kernel<Field>
			<< <grid, block, 0, stream >> >
			(A, B, C, n_elems);

		return cudaGetLastError();
	}



	// ===== toy shuffle emulation (u64) =====

	static inline
		void shuffle_down_u64(
			uint64_t* out,
			const uint64_t* in,
			int n,
			int shift)
	{
		if (shift <= 0)
			memcpy(out, in, sizeof(uint64_t) * n);
		else
		{
			int i;
			for (i = 0; i < n; i++)
			{
				int j = i + shift;
				out[i] = (j < n) ? in[j] : 0ull;
			}
		}
	}

	static inline
		void shuffle_up_u64(
			uint64_t* out,
			const uint64_t* in,
			int n,
			int shift)
	{
		if (shift <= 0)
			memcpy(out, in, sizeof(uint64_t) * n);
		else
		{
			int i;
			for (i = 0; i < n; i++)
			{
				int j = i - shift;
				out[i] = (j >= 0) ? in[j] : 0ull;
			}
		}
	}

	// ===== serial stand-in for the parallelPrefix_add over (s,c) =====
	// Computes sum[i] = s[i] + c[i] with carry propagation from i=0..n-1.
	// Returns final carry-out via *carry_out_columns (this is what the coop path would
	// read in the MS lane after the prefix adder).
	template<int n>
	static inline
		void ripple_add_sc(
			const uint64_t* s,
			const uint64_t* c,
			uint64_t* sum,
			uint64_t* carry_out_columns)
	{
		uint64_t carry = 0ull;
		int i;
		for (i = 0; i < n; i++)
		{
			wide_limb_t w;
			cc_add_wide_limb_limb(&w, s[i], c[i]);
			cc_add_wide_wide_limb(&w, w, carry);
			sum[i] = w.lo;
			carry = w.hi;
		}
		*carry_out_columns = carry;
	}

	BLST_HD BLST_FORCEINLINE
		uint64_t umul64hi(uint64_t a, uint64_t b);

	inline uint64_t carry(uint64_t x, uint64_t y, uint64_t z)
	{
		wide_limb_t sum;
		cc_add_wide_limb_limb(&sum, x, y);
		cc_add_wide_wide_limb(&sum, sum, z);
		return sum.hi;
	}

	static bool pp_print = false;

	// 5->2 compressor over m limbs of type Limb (default 64-bit).
	// Inputs : cin  (a one limb carry in. can have value [0,1,2,3])
	// Inputs : x0, x1, x2, x3, x4  (each m-limb n-bit vector)
	// Outputs: s, c  
	// 
	// s + c = x0 + x1 + x2 + x3 + x4 + cin
	//
	// Internals:
	//   Round 1 (5->3, column-local):  X0+...+X4 = S0 + 2*n0 + 4*n1
	//   Round 2 (3->2 CSA):            (S0), (n0<<1), (n1<<2)  ->  s + c
	//
	// No '+' is used inside a limb; everything is XOR/AND/OR/<< with cross-limb “spill” stitching.
	// The two possible overflow bits (at 2^n and 2^(n+1)) are dropped by design in this signature.
	template<int m, typename Limb = uint64_t>
	Limb compress5to2(
		Limb cin,
		const Limb* x0,
		const Limb* x1,
		const Limb* x2,
		const Limb* x3,
		const Limb* x4,
		Limb* s,
		Limb* c
	) {
		static_assert(std::is_unsigned<Limb>::value, "Limb must be unsigned");
		constexpr int W = int(sizeof(Limb) * 8);

		BLST_ASSERT(cin < 4);  // cin must be in {0,1,2,3}
		// -------- Round 1: 5->3 local; compute S0, n0, n1 and pack outgoing spills --------
		Limb S0[m];
		Limb n0[m];
		Limb n1[m];

		// pack3[i] holds 3 bits: [bit2=b2 | bit1=b1 | bit0=b0]
		// b0 = (n0 >> (W-1))       (spill for <<1)
		// b1 = (n1 >> (W-2)) & 1   (LSB spill for <<2)
		// b2 = (n1 >> (W-1)) & 1   (MSB spill for <<2)
		uint32_t pack3[m];

		for (int i = 0; i < m; i++) {
			Limb a = x0[i];
			Limb b = x1[i];
			Limb c0 = x2[i];
			Limb d = x3[i];
			Limb e = x4[i];

			// (a,b,c0) -> sA + 2*cA
			Limb sA = a ^ b ^ c0;
			Limb cA = (a & b) | (b & c0) | (a & c0);

			// (d,e,0) -> sB + 2*cB
			Limb sB = d ^ e;
			Limb cB = d & e;

			// sA + sB -> S0 + 2*cAB
			Limb s0 = sA ^ sB;
			Limb cAB = sA & sB;

			// N = cA + cB + cAB in {0,1,2,3}
			Limb n0_i = (cA ^ cB) ^ cAB;                       // bit0(N)
			Limb n1_i = (cA & cB) ^ (cA & cAB) ^ (cB & cAB);   // bit1(N)

			S0[i] = s0;
			n0[i] = n0_i;
			n1[i] = n1_i;

			uint32_t b0 = uint32_t(n0_i >> (W - 1));
			uint32_t b1 = uint32_t(n1_i >> (W - 2));
			pack3[i] = (b1 << 1) | b0;
		}

		uint32_t top_pack1 = pack3[m - 1];
		Limb spill_carry = (top_pack1 & 1) + (top_pack1 >> 1);

		// Communicate the 3 spill bits up by one limb (shift=1)
		uint32_t in3[m];
		shuffle_up(in3, pack3, m, /*shift=*/1);

		// -------- Form CSA inputs with received spills and do 3:2 CSA --------
		// A = S0
		// B = (n0<<1) | in_b0
		// C = (n1<<2) | (in_b1) | (in_b2<<1)
		// Then S = A^B^C; K = (A&B)|(A&C)|(B&C)
		Limb K_row[m];
		for (int i = 0; i < m; i++) {
			uint32_t in = i ? in3[i] : (cin << 1);
			uint32_t b0 = in & 1;
			uint32_t b1 = in >> 1;

			if (pp_print)
				std::cout << "i " << i << " cin " << b0 << b1 << std::endl;

			Limb A = S0[i];
			Limb B = (n0[i] << 1) | Limb(b0);
			Limb C = (n1[i] << 2) | Limb(b1);

			Limb S = (A ^ B) ^ C;
			Limb K = (A & B) | (A & C) | (B & C);

			s[i] = S;      // final sum row
			K_row[i] = K;      // unshifted carry row (needs <<1 across limbs to become 'c')

			pack3[i] = K_row[i] >> (W - 1);  // single spill bit
		}

		uint32_t top_pack2 = pack3[m - 1]; // Capture before shuffle drops it
		spill_carry += top_pack2;

		shuffle_up(in3, pack3, m, /*shift=*/1);

		for (int i = 0; i < m; i++) {
			auto in = i ? in3[i] : 0;
			c[i] = (K_row[i] << 1) | in;
		}

		if (0)
		{
			Limb exp[m + 1];
			exp[0] = cin;
			Limb cc = 0;
			for (u64 i = 0; i < m; ++i)
			{
				wide_limb_t w{ 0,0 };
				cc_add_wide_wide_limb(&w, w, exp[i]);
				cc_add_wide_wide_limb(&w, w, x0[i]);
				cc_add_wide_wide_limb(&w, w, x1[i]);
				cc_add_wide_wide_limb(&w, w, x2[i]);
				cc_add_wide_wide_limb(&w, w, x3[i]);
				cc_add_wide_wide_limb(&w, w, x4[i]);

				exp[i] = w.lo;
				exp[i + 1] = w.hi;


				//auto act = s[i] + c[i];
				cc_add_wide_limb_limb(&w, s[i], c[i]);
				cc_add_wide_wide_limb(&w, w, cc);
				cc = w.hi;
				if (exp[i] != w.lo)
				{
					std::cout << "bad at " << i << std::endl;
				}
			}
		}

		return spill_carry;
	}

	// ===== the toy serial "coop-emulated" FIOS step (lo-down plan) =====
	//
	// Field must provide:
	//   static constexpr const auto& modulus().data  (p[0..n-1])
	//   static constexpr uint64_t montgomery_constant()  (n0)
	//
	// t is size n+1 (t[0..n]), updated in place like in your serial step.
	// Returns final carry_out of the step (what you'd return from the coop kernel).
	//
	template<typename Field, int n, typename Carry = std::monostate>
	uint64_t montgomery_fios_step_lo_down_emulated(
		uint64_t* t,                    // t[0..n], updated in place
		const uint64_t* a,              // a[0..n-1]
		uint64_t bj,                    // b_j
		Carry carry_in = {})            // optional, like your serial version
	{
		static_assert(n > 0, "n must be > 0");

		constexpr bool has_carry = !std::is_same<Carry, std::monostate>::value;
		const uint64_t* p = Field::modulus().data;
		const uint64_t n0 = Field::montgomery_constant();

		// --- 1) Compute mx from column 0 and broadcast (scalar here) ---
		// mx = (ab_lo[0] + t[0]) * n0  mod 2^64    where ab_lo[0] = a[0] * bj

		// todo, compute all of the a[i]*bj in parallel first. Then we will have 
		// lane 0 compute mx and broadcast.
		uint64_t ab0_lo = a[0] * bj;
		uint64_t mx = [&] {
			if constexpr (has_carry)
				return (ab0_lo + t[0]) * n0;
			else
				return ab0_lo * n0;
			}();

		if (pp_print)
			std::cout << "mx = " << mx << "\n";


		// --- 2) Per-column products (full, unshifted) ---
		uint64_t ab_lo[n];
		uint64_t ab_hi[n];
		uint64_t mp_lo[n];
		uint64_t mp_hi[n];

		for (int i = 0; i < n; i++)
		{
			uint64_t ai = a[i];
			uint64_t pi = p[i];

			wide_limb_t w;

			cc_mul_wide_limb_limb(&w, ai, bj);
			ab_lo[i] = w.lo;
			ab_hi[i] = w.hi;

			cc_mul_wide_limb_limb(&w, pi, mx);
			mp_lo[i] = w.lo;
			mp_hi[i] = w.hi;
			if (i == 0 && pp_print)
			{
				std::cout << "mp = (" << w.hi << "," << w.lo << ") = " << mx << "*" << pi << std::endl;
			}
		}

		uint64_t carry0 = carry(ab0_lo, mp_lo[0], has_carry ? t[0] : 0);

		// --- 3) Shift LO vectors down by 1, and t down by 1 (lane n-1 gets t[n]) ---
		uint64_t ab_lo_d[n];
		uint64_t mp_lo_d[n];
		uint64_t t_d[n];

		shuffle_down_u64(ab_lo_d, ab_lo, n, 1);
		shuffle_down_u64(mp_lo_d, mp_lo, n, 1);
		shuffle_down_u64(t_d, t, n, 1);
		t_d[n - 1] = 0;// t[n];
		mp_lo[n - 1] = 0ull;
		ab_lo_d[n - 1] = 0ull;

		if constexpr (has_carry == false)
		{
			for (int i = 0; i < n; i++)
			{
				t_d[i] = 0;
			}
		}
		else
		{
			t_d[n - 1] = carry_in;
		}

		// --- 4) Column-local 5->2 compression ---
		uint64_t s[n];
		uint64_t c[n];

		if (pp_print)
		{
			std::cout
				<< "c0 " << carry0 << std::endl;
			for (int i = 0; i < n; i++)
			{
				std::cout << "i " << i << ":"
					<< " pm.lo " << mp_lo_d[i]
					<< " + pm.hi " << mp_hi[i]
					<< " + ab.lo " << ab_lo_d[i]
					<< " + ab.hi " << ab_hi[i]
					<< " + t " << t_d[i]
					<< "\n";
			}
		}

		// s + c = ab_lo_d + ab_hi + mp_lo_d + mp_hi + t_d + carry0
		uint64_t compressor_spills = compress5to2<n, uint64_t>(
			carry0,        // extra carry in at the bottom
			ab_lo_d,        // x0
			ab_hi,          // x1
			mp_lo_d,        // x2
			mp_hi,          // x3
			t_d,            // x4
			s,
			c);

		// --- 5) "Prefix" add of (s, c). Serial ripple here for correctness. ---
		uint64_t sum[n];
		uint64_t carry_out_columns = 0ull;

		// (sum, carry_out_columns) = s + c
		parallelPrefix_add<n>(sum, s, c, &carry_out_columns);
		carry_out_columns += compressor_spills;

		if (pp_print)
		{

			for (int i = 0; i < n; i++)
			{
				std::cout << "sum[" << i << "] = " << sum[i] << "\n";
			}
		}

		// --- 8) Write back: t[i] = sum[i] for i=0..n-2; t[n-1] = t_last ---
		for (int i = 0; i < n; i++)
		{
			t[i] = sum[i];
		}


		if (pp_print)
		{
			//t[n - 1] = t_last;
			std::cout << "carry " << carry_out_columns << std::endl;
		}

		return carry_out_columns;
	}


	template<typename Field>
	BLST_HD BLST_FORCEINLINE
		void par_mul_mont_u64(
			uint64_t* r,
			const uint64_t* a,
			const uint64_t* b)
	{
		using namespace blst;

		constexpr auto p = Field::modulus().data;
		constexpr int n = sizeof(Field) / sizeof(uint64_t);
		static_assert(n > 0 && (n % 2 == 0), "n must be positive and even");

		uint64_t aa[n];
		for (int i = 0; i < n; i++) {
			aa[i] = a[i];
		}
		uint64_t bb[n];
		for (int i = 0; i < n; i++) {
			bb[i] = b[i];
		}

		uint64_t t[n];// = { 0 };

		uint64_t carry = montgomery_fios_step_lo_down_emulated<Field, n>(t, aa, bb[0]);

		// Columns j = 1..n-1 (fused)
		for (int j = 1; j < n; ++j) {
			carry = montgomery_fios_step_lo_down_emulated<Field, n>(t, aa, bb[j], carry);
		}

		// Final conditional subtract, (todo: parallel)
		uint64_t tmp[n];
		// 1. Compute tmp = t - p and get the borrow out
		uint64_t borrow = parallelPrefix_sub<n>(tmp, t, p);

		uint64_t mask = carry - borrow;
		for (int i = 0; i < n; i++) {
			r[i] = (tmp[i] & ~mask) | (t[i] & mask);
		}
	}


	//template<typename Field>
	//BLST_HD BLST_FORCEINLINE
	//	void par_mul_mont_u64(
	//		uint64_t* r,
	//		const uint64_t* a,
	//		const uint64_t* b)
	//{
	//	using namespace blst;

	//	constexpr auto p = Field::modulus().data;
	//	constexpr int n = sizeof(Field) / sizeof(uint64_t);
	//	static_assert(n > 0 && (n % 2 == 0), "n must be positive and even");

	//	uint64_t aa[n];
	//	for (int i = 0; i < n; i++) {
	//		aa[i] = a[i];
	//	}
	//	uint64_t bb[n];
	//	for (int i = 0; i < n; i++) {
	//		bb[i] = b[i];
	//	}

	//	// --- DEBUG START ---
	//	bool debug = true;
	//	if (debug) {
	//		printf("\n=== CPU DEBUG ===\n");
	//		printf("n0: %llx\n", (unsigned long long)Field::montgomery_constant());
	//		for (int i = 0; i < n; ++i) {
	//			printf("Rank %d: A=%llx P=%llx B_own=(n/a)\n", i, (unsigned long long)aa[i], (unsigned long long)p[i]);
	//		}
	//	}
	//	// --- DEBUG END ---

	//	uint64_t t[n];// = { 0 };
	//	// Note: For Step 0, t is treated as 0 by the template logic, 
	//	// but for printing purposes we conceptually treat t_in as 0.

	//	uint64_t carry = montgomery_fios_step_lo_down_emulated<Field, n>(t, aa, bb[0]);

	//	// --- DEBUG STEP 0 ---
	//	if (debug) {
	//		printf("--- Step 0 (bj=%llx) ---\n", (unsigned long long)bb[0]);
	//		for (int i = 0; i < n; ++i) {
	//			// For step 0, t_in was implicitly 0
	//			printf("Rank %d: t_in=0 -> t_out=%llx (carry_state=%llx)\n",
	//				i, (unsigned long long)t[i], (unsigned long long)carry);
	//		}
	//	}
	//	// --------------------

	//	// Columns j = 1..n-1 (fused)
	//	for (int j = 1; j < n; ++j) {

	//		// Capture previous state for printing
	//		uint64_t t_prev[n];
	//		if (debug) { for (int k = 0; k < n; ++k) t_prev[k] = t[k]; }

	//		carry = montgomery_fios_step_lo_down_emulated<Field, n>(t, aa, bb[j], carry);

	//		// --- DEBUG LOOP ---
	//		if (debug) {
	//			printf("--- Step %d (bj=%llx) ---\n", j, (unsigned long long)bb[j]);
	//			for (int i = 0; i < n; ++i) {
	//				printf("Rank %d: t_in=%llx -> t_out=%llx (carry_state=%llx)\n",
	//					i, (unsigned long long)t_prev[i], (unsigned long long)t[i], (unsigned long long)carry);
	//			}
	//		}
	//		// ------------------
	//	}

	//	if (debug) {
	//		printf("\n--- Pre-Subtraction ---\n");
	//		for (int i = 0; i < n; ++i) printf("Rank %d: Final T = %llx\n", i, (unsigned long long)t[i]);
	//	}

	//	// Final conditional subtract, (todo: parallel)
	//	uint64_t tmp[n];
	//	// 1. Compute tmp = t - p and get the borrow out
	//	uint64_t borrow = parallelPrefix_sub<n>(tmp, t, p);

	//	uint64_t mask = carry - borrow;

	//	// --- DEBUG FINAL ---
	//	if (debug) {
	//		printf("\n--- Final Logic ---\n");
	//		printf("Carry Top: %llx\n", (unsigned long long)carry);
	//		printf("Total Borrow: %llx\n", (unsigned long long)borrow);
	//		printf("Mask: %llx\n", (unsigned long long)mask);

	//		for (int i = 0; i < n; ++i) {
	//			// Replicate the masking logic for display
	//			uint64_t res = (tmp[i] & ~mask) | (t[i] & mask);
	//			printf("Rank %d: T-P=%llx | Masked T-P=%llx | Masked T=%llx\n",
	//				i, (unsigned long long)tmp[i], (unsigned long long)(tmp[i] & ~mask), (unsigned long long)(t[i] & mask));
	//		}
	//	}
	//	// -------------------

	//	for (int i = 0; i < n; i++) {
	//		r[i] = (tmp[i] & ~mask) | (t[i] & mask);
	//	}

	//}


	// --- Device Helper: 5->2 Compressor (Tile Distributed) ---
	// Returns the spill carry (only valid in the top thread of the tile)
	template<int NLIMBS>
	__device__ __forceinline__ uint64_t compress5to2_tile(
		Tile<NLIMBS> tile,
		uint64_t cin,       // Only relevant for rank 0
		uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4,
		uint64_t& s, uint64_t& c
	) {
		// Round 1: 5->3 local (Bitwise ops)
		uint64_t sA = x0 ^ x1 ^ x2;
		uint64_t cA = (x0 & x1) | (x1 & x2) | (x0 & x2);

		uint64_t sB = x3 ^ x4;
		uint64_t cB = x3 & x4;

		uint64_t s0 = sA ^ sB;
		uint64_t cAB = sA & sB;

		uint64_t n0_i = (cA ^ cB) ^ cAB;
		uint64_t n1_i = (cA & cB) ^ (cA & cAB) ^ (cB & cAB);

		// Pack spills: b0 (from n0 >> 63), b1 (n1 >> 62), b2 (n1 >> 63)
		uint32_t b0 = (uint32_t)(n0_i >> 63);
		uint32_t b1 = (uint32_t)(n1_i >> 62);
		uint32_t pack3 = (b1 << 1) | b0;

		// Shuffle up spills
		uint32_t in3 = tile.shfl_up(pack3, 1);

		// Handle bottom Carry-In for Rank 0
		uint32_t bottom_in = (tile.thread_rank() == 0) ? (uint32_t)(cin << 1) : 0;
		uint32_t in_bits = (tile.thread_rank() == 0) ? bottom_in : in3;

		uint32_t in_b0 = in_bits & 1;
		uint32_t in_b1 = in_bits >> 1;

		// CSA Round 2
		uint64_t A = s0;
		uint64_t B = (n0_i << 1) | (uint64_t)in_b0;
		uint64_t C = (n1_i << 2) | (uint64_t)in_b1;

		s = A ^ B ^ C;
		uint64_t K = (A & B) | (A & C) | (B & C);

		// Final C formation
		uint32_t k_spill = (uint32_t)(K >> 63);
		uint32_t k_in = tile.shfl_up(k_spill, 1);
		if (tile.thread_rank() == 0) k_in = 0;

		c = (K << 1) | k_in;

		// Calculate top spill (only relevant for top thread)
		// From emulation: spill_carry = (top_pack1 & 1) + (top_pack1 >> 1) + top_pack2
		uint64_t spill = 0;
		if (tile.thread_rank() == NLIMBS - 1) {
			spill = (pack3 & 1) + (pack3 >> 1) + k_spill;
		}
		return spill;
	}

	// --- Device Helper: Parallel Prefix Add (S + C) ---
	// Returns carry out from the top limb
	template<int NLIMBS>
	__device__ __forceinline__ uint64_t pp_add_sc_tile(
		Tile<NLIMBS> tile,
		uint64_t& sum_out,
		uint64_t s,
		uint64_t c
	) {
		// 1. Generate Local G/P
		uint64_t raw_sum = s + c;
		uint32_t g = (raw_sum < s) ? 1u : 0u;
		uint32_t p = (raw_sum == ~0ull) ? 1u : 0u;
		uint32_t pack = g | (p << 8);

		// 2. Scan
		pack = ks_packed_prefix_tile<NLIMBS>(tile, pack);

		// 3. Shift to get Carry-In
		uint32_t prev_pack = tile.shfl_up(pack, 1);
		uint64_t cin = (tile.thread_rank() == 0) ? 0ull : (uint64_t)(prev_pack & 1u);

		// 4. Final Sum
		sum_out = s + c + cin;

		// 5. Return global carry out (only valid in top lane)
		return (tile.thread_rank() == NLIMBS - 1) ? (uint64_t)(pack & 1u) : 0;
	}

	// --- Core Device Function: FIOS Step ---
	// Updates local 't' register, returns carry for next iteration
	template<int NLIMBS>
	__device__ __forceinline__ uint64_t mont_fios_step_tile(
		Tile<NLIMBS> tile,
		uint64_t& t,
		uint64_t a_i,       // A[rank]
		uint64_t p_i,       // P[rank]
		uint64_t bj,        // Broadcasted B[j]
		uint64_t n0,
		uint64_t carry_in   // From previous step (fed into top)
	) {
		// 2. Full a[i]*b[j] Product (128-bit)
		uint64_t ab_hi = __umul64hi(a_i, bj);
		uint64_t ab_lo = a_i * bj;

		// 1. Compute m (Lane 0 only) and Broadcast
		// m = (t[0] + a[0]*bj) * n0 mod 2^64
		uint64_t m = 0;
		if (tile.thread_rank() == 0) {
			m = (t + ab_lo) * n0;
		}
		m = tile.shfl(m, 0);


		// 2. Full p[i]*m Product (128-bit)
		uint64_t mp_hi = __umul64hi(p_i, m);
		uint64_t mp_lo = p_i * m;

		// 3. Calculate Carry0 (Lane 0 only) - The "Extra" carry
		// carry0 = carry(ab_lo, mp_lo, t)
		uint64_t carry0 = 0;
		if (tile.thread_rank() == 0) {
			// Emulate: carry(x,y,z) -> (x+y+z) >> 64
			// Using PTX for efficiency or simple C++
			// sum = ab_lo + mp_lo + t;
			// carry is generated if sum wraps or if components wrapped
			uint64_t r1 = ab_lo + mp_lo;
			uint64_t c1 = (r1 < ab_lo);
			uint64_t r2 = r1 + t;
			uint64_t c2 = (r2 < r1);
			carry0 = c1 + c2;
		}

		// 4. Shuffle Down (LO parts and T)
		// rank i gets value from rank i+1. Top gets 0 (or carry_in for T)
		uint64_t ab_lo_d = tile.shfl_down(ab_lo, 1);
		uint64_t mp_lo_d = tile.shfl_down(mp_lo, 1);
		uint64_t t_d = tile.shfl_down(t, 1);

		if (tile.thread_rank() == NLIMBS - 1) {
			ab_lo_d = 0;
			mp_lo_d = 0;
			t_d = carry_in; // Inject carry from previous iteration into top
		}

		// 5. Compressor (5 -> 2)
		// Inputs: carry0 (rank 0), ab_lo_d, ab_hi, mp_lo_d, mp_hi, t_d
		uint64_t s, c;
		uint64_t spill = compress5to2_tile<NLIMBS>(
			tile, carry0,
			ab_lo_d, ab_hi, mp_lo_d, mp_hi, t_d,
			s, c
		);

		// 6. Parallel Prefix Add (S + C)
		uint64_t sum;
		uint64_t pp_carry = pp_add_sc_tile<NLIMBS>(tile, sum, s, c);

		// 7. Update T
		t = sum;

		// 8. Return total carry (Spill + PP Carry) - Only valid in top lane
		return spill + pp_carry;
	}


	// --- Main Kernel Function ---
	template<typename Field>
	__device__ __forceinline__ void mul_mont_coop_tile(
		Tile<Field::num_limb()> tile,
		auto&& R,
		auto&& A,
		auto&& B
	) {
		constexpr int NLIMBS = Field::num_limb();
		const int rank = tile.thread_rank();

		// Load A and P into registers
		uint64_t a_i = get_coop(A, rank);
		uint64_t b_i = get_coop(B, rank);
		uint64_t p_i = Field::modulus().data[rank];
		uint64_t n0 = Field::montgomery_constant();

		// Initialize Accumulator T
		uint64_t t = 0;
		uint64_t carry = 0; // Carry passed between iterations

		// Iterate over limbs of B
#pragma unroll
		for (int j = 0; j < NLIMBS; ++j) {
			// Load B[j] via shuffle (Lane j reads B[rank], then broadcast)
			// Or simpler: Load B[j] from global memory (L1 cache handles it well)
			// Ideally: One thread loads, broadcast.
			// Let's rely on cache broadcast for now or simple load:
			//uint64_t bj = B[j]; // Every thread loads B[j]
			uint64_t bj = tile.shfl(b_i, j);
			carry = mont_fios_step_tile<NLIMBS>(tile, t, a_i, p_i, bj, n0, carry);
		}

		// Final Subtraction (Conditional)
		// 1. Calculate T - P
		// Using parallelPrefix_sub logic locally
		uint64_t t_sub_p;
		uint32_t g = (t < p_i) ? 1u : 0u;
		uint32_t p = (t == p_i) ? 1u : 0u;
		uint32_t pack = g | (p << 8);

		pack = ks_packed_prefix_tile<NLIMBS>(tile, pack);
		uint32_t prev_pack = tile.shfl_up(pack, 1);
		uint64_t borrow_in = (rank == 0) ? 0 : (prev_pack & 1u);

		t_sub_p = t - p_i - borrow_in;
		uint64_t total_borrow = (tile.shfl(pack, NLIMBS - 1) & 1u);

		// Logic: if Carry=1, result > P (definitely), do subtraction (implicitly done by wrapping, result is t_sub_p).
		// If Carry=0, check borrow. If Borrow=1, T < P, use T. If Borrow=0, T >= P, use T-P.
		// Actually: Real result is T (extended by carry).
		// If (Carry - Borrow) == 0, it means T < P, keep T.
		// If (Carry - Borrow) == 1, it means T >= P, use T - P.

		// Broadcast carry from top lane
		uint64_t carry_top = tile.shfl(carry, NLIMBS - 1);

		uint64_t mask = carry_top - total_borrow; // 1 if subtract needed, 0 otherwise
		BLST_ASSERT(mask != 1ull); // should be 0 or all 1s

		// Select
		get_coop(R,rank) = (t_sub_p & ~mask) | (t & mask);
	}
	//
	//	template<typename Field>
	//	__device__ __forceinline__ void mul_mont_coop_tile(
	//		Tile<Field::num_limb()> tile,
	//		limb_t* R,
	//		const limb_t* A,
	//		const limb_t* B
	//	) {
	//		constexpr int NLIMBS = Field::num_limb();
	//		const int rank = tile.thread_rank();
	//
	//		// Debug: Only print for the first instance (Global Tile ID 0)
	//		// We need to calculate global ID to filter
	//		int warp_id = blockIdx.x * (blockDim.x / 32) + (threadIdx.x / 32);
	//		int inst_idx = warp_id * tile.tiles_per_warp() + tile.mGid;
	//		bool debug = (inst_idx == 0);
	//
	//		uint64_t a_i = A[rank];
	//		// Note: Using the original loading logic to minimize variables changed
	//		// (If you kept the optimization, verify b_own here)
	//		uint64_t b_own = B[rank];
	//
	//		uint64_t p_i = Field::modulus().data[rank];
	//		uint64_t n0 = Field::montgomery_constant();
	//
	//		if (debug && rank == 0) printf("\n=== DEBUG INSTANCE 0 ===\nn0: %llx\n", (unsigned long long)n0);
	//		if (debug) printf("Rank %d: A=%llx P=%llx B_own=%llx\n", rank, (unsigned long long)a_i, (unsigned long long)p_i, (unsigned long long)b_own);
	//
	//		uint64_t t = 0;
	//		uint64_t carry = 0;
	//
	//#pragma unroll
	//		for (int j = 0; j < NLIMBS; ++j) {
	//			uint64_t bj = tile.shfl(b_own, j);
	//
	//			// Snapshot t before step
	//			uint64_t t_prev = t;
	//
	//			carry = mont_fios_step_tile<NLIMBS>(tile, t, a_i, p_i, bj, n0, carry);
	//
	//			if (debug) {
	//				// Sync to keep prints orderly
	//				__syncwarp();
	//				if (rank == 0) printf("--- Step %d (bj=%llx) ---\n", j, (unsigned long long)bj);
	//				// Print T state after this step
	//				printf("Rank %d: t_in=%llx -> t_out=%llx (carry_state=%llx)\n",
	//					rank, (unsigned long long)t_prev, (unsigned long long)t, (unsigned long long)carry);
	//			}
	//		}
	//
	//		if (debug) {
	//			__syncwarp();
	//			if (rank == 0) printf("\n--- Pre-Subtraction ---\n");
	//			printf("Rank %d: Final T = %llx\n", rank, (unsigned long long)t);
	//		}
	//
	//		// --- Final Subtraction Logic ---
	//		uint64_t t_sub_p;
	//		uint32_t g = (t < p_i) ? 1u : 0u;
	//		uint32_t p = (t == p_i) ? 1u : 0u;
	//		uint32_t pack = g | (p << 8);
	//
	//		pack = ks_packed_prefix_tile<NLIMBS>(tile, pack);
	//
	//		uint32_t prev_pack = tile.shfl_up(pack, 1);
	//		uint64_t borrow_in = (rank == 0) ? 0 : (prev_pack & 1u);
	//
	//		t_sub_p = t - p_i - borrow_in;
	//
	//		uint64_t total_borrow = (tile.shfl(pack, NLIMBS - 1) & 1u);
	//		uint64_t carry_top = tile.shfl(carry, NLIMBS - 1);
	//
	//		uint64_t mask = carry_top - total_borrow;
	//
	//		if (debug) {
	//			__syncwarp();
	//			if (rank == 0) {
	//				printf("\n--- Final Logic ---\n");
	//				printf("Carry Top: %llx\n", (unsigned long long)carry_top);
	//				printf("Total Borrow: %llx\n", (unsigned long long)total_borrow);
	//				printf("Mask: %llx\n", (unsigned long long)mask);
	//			}
	//			printf("Rank %d: T-P=%llx | Masked T-P=%llx | Masked T=%llx\n",
	//				rank, (unsigned long long)t_sub_p, (unsigned long long)(t_sub_p & mask), (unsigned long long)(t & ~mask));
	//		}
	//
	//		R[rank] = (t_sub_p & ~mask) | (t & mask);
	//	}

		// --- Kernel Wrapper for Testing ---
	template<typename Field>
	__global__ void mul_mont_coop_kernel(
		limb_t* R,
		const limb_t* A,
		const limb_t* B,
		int num_ops
	) {
		constexpr int NLIMBS = Field::num_limb();
		Tile<NLIMBS> tile;

		if (!tile.is_full()) return;

		int warp_id = blockIdx.x * (blockDim.x / 32) + (threadIdx.x / 32);
		int tile_idx = warp_id * tile.tiles_per_warp() + tile.mGid;

		if (tile_idx >= num_ops) return;

		// Pointers for this instance (AoS layout)
		const limb_t* a_ptr = A + tile_idx * NLIMBS;
		const limb_t* b_ptr = B + tile_idx * NLIMBS;
		limb_t* r_ptr = R + tile_idx * NLIMBS;

		mul_mont_coop_tile<Field>(tile, r_ptr, a_ptr, b_ptr);
	}



	// --- Wrapper for compress5to2_tile ---
	template<int NLIMBS>
	__global__ void compress5to2_test_kernel(
		uint64_t* d_s,
		uint64_t* d_c,
		uint64_t* d_spill, // Output array for spills (1 per instance)
		const uint64_t* d_cin,   // Input array of carries (1 per instance)
		const uint64_t* x0, const uint64_t* x1, const uint64_t* x2,
		const uint64_t* x3, const uint64_t* x4,
		int num_ops
	) {
		Tile<NLIMBS> tile;
		if (!tile.is_full()) return;

		int warp_id = blockIdx.x * (blockDim.x / 32) + (threadIdx.x / 32);
		int inst_idx = warp_id * tile.tiles_per_warp() + tile.mGid;

		if (inst_idx >= num_ops) return;

		// AoS Indexing
		int rank = tile.thread_rank();
		int offset = inst_idx * NLIMBS + rank;

		// Load Inputs
		uint64_t v0 = x0[offset];
		uint64_t v1 = x1[offset];
		uint64_t v2 = x2[offset];
		uint64_t v3 = x3[offset];
		uint64_t v4 = x4[offset];

		// Cin is only relevant for rank 0, but we load it to pass it.
		// We assume d_cin has 1 entry per instance.
		uint64_t cin = (rank == 0) ? d_cin[inst_idx] : 0;

		uint64_t s, c;
		uint64_t spill = compress5to2_tile<NLIMBS>(tile, cin, v0, v1, v2, v3, v4, s, c);

		// Store outputs
		d_s[offset] = s;
		d_c[offset] = c;

		// Store spill (only top thread writes)
		if (rank == NLIMBS - 1) {
			d_spill[inst_idx] = spill;
		}
	}

	// --- Wrapper for Wmont_fios_step_tile ---
	template<int NLIMBS>
	__global__ void fios_step_test_kernel(
		uint64_t* d_t,       // In/Out
		uint64_t* d_carry_out, // One per instance
		const uint64_t* d_a,
		const uint64_t* d_p,
		const uint64_t* d_b, // Array of b inputs (one per instance)
		const uint64_t* d_carry_in, // One per instance
		uint64_t n0,
		int num_ops
	) {
		Tile<NLIMBS> tile;
		if (!tile.is_full()) return;

		int warp_id = blockIdx.x * (blockDim.x / 32) + (threadIdx.x / 32);
		int inst_idx = warp_id * tile.tiles_per_warp() + tile.mGid;

		if (inst_idx >= num_ops) return;

		int rank = tile.thread_rank();
		int offset = inst_idx * NLIMBS + rank;

		// Load Inputs
		uint64_t t = d_t[offset];
		uint64_t a = d_a[offset];
		uint64_t p = d_p[offset];

		// Scalar inputs (per instance)
		uint64_t bj = d_b[inst_idx];
		uint64_t cin = d_carry_in[inst_idx];

		uint64_t carry_ret = mont_fios_step_tile<NLIMBS>(tile, t, a, p, bj, n0, cin);

		// Store result
		d_t[offset] = t;

		// Store carry (only top thread writes)
		if (rank == NLIMBS - 1) {
			d_carry_out[inst_idx] = carry_ret;
		}
	}


} // namespace blst