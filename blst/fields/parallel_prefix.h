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

	void shuffle_down(
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

	void shuffle_up(
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
	inline void ks_packed_prefix(uint32_t gp[NLIMBS])
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
	template<int NLIMBS>
	inline void pp_add(
		limb_t c[NLIMBS],          // out: a + b
		const limb_t a[NLIMBS],
		const limb_t b[NLIMBS],
		const limb_t pneg[NLIMBS]  // p' = -p mod 2^(64*NLIMBS)
	)
	{
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

	// One tile computes c = (a + b) mod p using p' = -p
	// AoS layout: arrays A,B,PNEG,C are contiguous by element: base + elem*NLIMBS + limb
	template<int NLIMBS>
	BLST_DEVICE BLST_FORCEINLINE
		void pp_add_coop_tile_aos(Tile<NLIMBS> tile,
			const limb_t* BLST_restrict A,
			const limb_t* BLST_restrict B,
			const limb_t* BLST_restrict PNEG,
			limb_t* BLST_restrict C,
			size_t elem_idx)
	{
		static_assert(NLIMBS == 6 || NLIMBS == 4, "supported: 256/384 bits");

		const int limb = tile.thread_rank();

		const size_t off = elem_idx * (size_t)NLIMBS + (size_t)limb;

		limb_t a_i = A[off];
		limb_t b_i = B[off];
		limb_t pn = PNEG[limb];

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

		C[off] = c_i;
	}

	// Kernel: AoS inputs, handles tail tiles safely.
	// For NLIMBS=6, each warp has 5 full tiles + 1 tail of size 2; we skip the tail.
	template<int NLIMBS>
	BLST_GLOBAL void pp_add_coop_aos_kernel(
		const limb_t* BLST_restrict A,
		const limb_t* BLST_restrict B,
		const limb_t* BLST_restrict PNEG,
		limb_t* BLST_restrict C,
		size_t n_elems)
	{
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

		pp_add_coop_tile_aos<NLIMBS>(tile, A, B, PNEG, C, elem_idx);
	}

	template<int NLIMBS>
	inline cudaError_t launch_pp_add_mod_aos(const limb_t* A,
		const limb_t* B,
		const limb_t* PNEG,
		limb_t* C,
		size_t n_elems,
		int warps_per_block = 4,
		cudaStream_t stream = 0)
	{
		static_assert(NLIMBS == 6 || NLIMBS == 4, "NLIMBS must be 4 or 6");

		if (n_elems == 0)
			return cudaSuccess;

		if (A == nullptr)
			return cudaErrorInvalidDevicePointer;
		if (B == nullptr)
			return cudaErrorInvalidDevicePointer;
		if (PNEG == nullptr)
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

		pp_add_coop_aos_kernel<NLIMBS> << <grid, block, 0, stream >> > (A, B, PNEG, C, n_elems);

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

	uint64_t umul64hi(uint64_t a, uint64_t b);

	uint64_t carry(uint64_t x, uint64_t y, uint64_t z)
	{
		wide_limb_t sum;
		cc_add_wide_limb_limb(&sum, x, y);
		cc_add_wide_wide_limb(&sum, sum, z);
		return sum.hi;
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
		uint64_t ab0_lo = a[0] * bj;
		uint64_t mx = [&] {
			if constexpr (has_carry)
				return (ab0_lo + t[0]) * n0;
			else
				return ab0_lo * n0;
			}();

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
			if (i == 0)
			{
				std::cout << "mp = (" << w.hi << ","<<w.lo<<") = " << mx << "*" << pi << std::endl;

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

		// NOTE: we do NOT shift the HI vectors. We use ab_hi[i], mp_hi[i] as-is.

		// --- 4) Column-local 5->2 compression ---
		uint64_t s[n];
		uint64_t c[n];

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

		// s + c = ab_lo_d + ab_hi + mp_lo_d + mp_hi + t_d + carry0
		compress5to2<n, uint64_t>(
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
		ripple_add_sc<n>(s, c, sum, &carry_out_columns);
		//pp_add<n>(sum, s, c, Field::neg_modulus().data);

		//if constexpr (has_carry)
		//{
		//	// Mirror: wide1 = hi_mp + hi_ab + carry_in
		//	// At this point, sum[n-1] is the top limb as if carry_in == 0,
		//	// and carry_out_columns is the "base" carry.
		//	wide_limb_t w;
		//	cc_add_wide_limb_limb(&w, sum[n - 1], (uint64_t)carry_in);
		//	sum[n - 1] = w.lo;
		//	carry_out_columns += w.hi;   // w.hi is 0 or 1
		//	std::cout << " + " << carry_in << std::endl;
		//}

		for (int i = 0; i < n; i++)
		{
			std::cout << "sum[" << i << "] = " << sum[i] << "\n";
		}

		uint64_t carry_from_tn = 0ull;

		// --- 8) Write back: t[i] = sum[i] for i=0..n-2; t[n-1] = t_last ---
		for (int i = 0; i < n; i++)
		{
			t[i] = sum[i];
		}
		//t[n - 1] = t_last;
		std::cout << "carry " << carry_out_columns << std::endl;

		return carry_out_columns;
	}


} // namespace blst
