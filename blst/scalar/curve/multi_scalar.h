/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/scalar/curve/point.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/bulk_addition.h"

#ifdef __CUDACC__
#include <cuda/std/span>
#else
#include  <span>
#endif
namespace blst::scalar
{
#ifdef __CUDACC__
	template<typename T>
	using span = cuda::std::span<T>;
#else
	template<typename T>
	using span = std::span<T>;
#endif

	namespace detail
	{
		BLST_HD BLST_INLINE constexpr bool valid_wbits(size_t wbits)
		{
			return wbits >= 2 && wbits <= 14;
		}

		BLST_HD BLST_INLINE constexpr bool checked_mul_size(size_t a, size_t b, size_t& out)
		{
			if (a != 0 && b > static_cast<size_t>(-1) / a)
				return false;
			out = a * b;
			return true;
		}

		BLST_HD BLST_INLINE constexpr bool checked_shl_size(size_t value, size_t shift, size_t& out)
		{
			if (shift >= sizeof(size_t) * 8 || value > (static_cast<size_t>(-1) >> shift))
				return false;
			out = value << shift;
			return true;
		}
	}


	//
	BLST_HD BLST_INLINE auto msm_to_jacobian(Point<Fp>& out, const Pointxyzz<Fp>& in) { Pointxyzz<Fp>::to_Jacobian(out, in); }
	BLST_HD BLST_INLINE auto msm_to_jacobian(Point<Fp2>& out, const Pointxyzz<Fp2>& in) { Pointxyzz<Fp2>::to_Jacobian(out, in); }

	BLST_HD BLST_INLINE auto msm_to_jacobian(GT& out, const GT& in) { out = in; }

	BLST_HD BLST_INLINE auto msm_add_or_sub(Pointxyzz<Fp>& out,
		const Pointxyzz<Fp>& p1,
		const AffinePoint<Fp>& p2,
		bool_t subtract)
	{
		Pointxyzz<Fp>::dadd(out, p1, p2, subtract);
	}

	BLST_HD BLST_INLINE auto msm_add_or_sub(Pointxyzz<Fp2>& out,
		const Pointxyzz<Fp2>& p1,
		const AffinePoint<Fp2>& p2,
		bool_t subtract)
	{
		Pointxyzz<Fp2>::dadd(out, p1, p2, subtract);
	}

	BLST_HD BLST_INLINE auto msm_add_or_sub(GT& out,
		const GT& p1,
		const GT& p2,
		bool_t subtract)
	{
		GT neg = p2;
		GT::cneg(neg, subtract);
		GT::add(out, p1, neg);
	}



	BLST_HD inline void print_bits(char* buffer, const void* data, int bits)
	{
		int bytes = (bits + 7) / 8;
		unsigned char* p = (unsigned char*)data;
		for (int i = bytes - 1; i >= 0; i--)
		{
			for (int j = 7; j >= 0; j--)
			{
				int bit_idx = i * 8 + j;
				if (bit_idx < bits)
				{
					*buffer++ = ((p[i] >> j) & 1) ? '1' : '0';
				}
			}
		}
		*buffer = 0;
	}

	template<typename Group>
	BLST_HD BLST_INLINE u64 hash_point(const Group& p)
	{
		u64 hash = 0;
		const byte* data = reinterpret_cast<const byte*>(&p);
		for (size_t i = 0; i < sizeof(Group); i++)
		{
			hash = (hash * 31) + data[i];
			hash ^= (hash >> 4);
		}
		return hash;
	}

	/*
	* This is two-step multi-scalar multiplication procedure. First, given
	* a set of points you pre-compute a table for chosen windowing factor
	* [expressed in bits with value between 2 and 14], and then you pass
	* this table to the actual multiplication procedure along with scalars.
	* Idea is that the pre-computed table will be reused multiple times. In
	* which case multiplication runs faster than below Pippenger algorithm
	* implementation for up to ~16K points for wbits=8, naturally at the
	* expense of multi-megabyte table. One can trade even more memory for
	* performance, but each wbits increment doubles the memory requirement,
	* so at some point it gets prohibively large... For reference, without
	* reusing the table it's faster than Pippenger algorithm for up ~32
	* points [with wbits=5]...
	*/

#define SCRATCH_SZ(ptype) (sizeof(ptype)==sizeof(G1Point) ? 8192 : 4096)

		/* The intermediate infinity points are encoded as [0, 0, 1]. */

	// For one affine point, builds a row of length 2 ^ (wbits - 1) of 
	// consecutive multiples used by Booth / gather(e.g., P·1, P·2, P·3, …) 
	// in Jacobian.
	template<typename Field>
	BLST_HD BLST_CPP_INLINE void point_precompute_row_wbits(
		Point<Field> row[],
		size_t wbits,
		const AffinePoint<Field>& point);

	// Orchestrates precompute for all input points in chunks (to bound memory), 
	// calling the two routines above and producing a flat affine table of size 
	// npoints << (wbits-1).
	// flat |points[n]| can be placed at the end of |table[n<<(wbits-1)]| 
	template<typename Field>
	BLST_HD BLST_CPP_INLINE void msm_precompute_wbits(
		AffinePoint<Field> table[],
		size_t wbits,
		const AffinePoint<Field>* const points[],
		size_t npoints);

	template<typename Field>
	size_t points_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints)
	{
		if (!detail::valid_wbits(wbits))
			return 0;
		size_t bytes = 0;
		if (!detail::checked_mul_size(sizeof(AffinePoint<Field>), npoints, bytes))
			return 0;
		if (!detail::checked_shl_size(bytes, wbits - 1, bytes))
			return 0;
		return bytes;
	}



	// Given a Booth-encoded digit (with sign bit), selects the appropriate
	// affine entry from a row, applies sign, or infinity if zero.
	template<typename Field>
	BLST_HD BLST_INLINE void point_gather_booth_wbits(
		Point<Field>& p,
		const AffinePoint<Field> row[],
		size_t wbits, limb_t booth_encoding)
	{
		bool_t booth_sign = (booth_encoding >> wbits) & 1;
		bool_t idx_is_zero;
		const AffinePoint<Field> infinity = { };

		auto booth_idx = booth_encoding & (((limb_t)1 << wbits) - 1);
		idx_is_zero = is_zero(booth_idx);
		booth_idx -= 1 ^ idx_is_zero;
		vec_select(&p, &infinity, &row[booth_idx], sizeof(row[0]), idx_is_zero);
		Point<Field>::cneg(p, booth_sign);
	}


	// Fixed-table MSM : scans scalar windows (top partial then full windows), 
	// gathers one entry per point(into scratch), periodically accumulates 
	// into ret, and does wbits doublings between windows. (Uses Booth 
	// encoding; avoids per-point additions into shared accumulators.)
	// 
	// best when npoints is in the few hundreds or less and table
	// is reused multiple times; also faster than Pippenger for up to ~32
	//
	// @ret the resulting point
	// @table precomputed matrix of npoints by (1 << (wbits-1))  points. 
	//   The ith contigous row of (1 << (wbits-1) points are the precomputed
	//   for the ith point. Obtained via msm_precompute_wbits() above.
	// @wbits window_size size in bits (between 2 and 14)
	// @npoints number of points (and scalars)
	// @scalars array of npoints pointers to scalars, each scalar is
	//   represented in byte array of length (nbits+7)/8.
	// @nbits number of bits in each scalar
	// @scratch scratch space for staging gathers before accumulation;
	//   suggested size is min(npoints, SCRATCH_SZ(Point<Field>)).
	template<typename Field>
	BLST_HD BLST_CPP_INLINE void msm_mult_table(
		Point<Field>& ret,
		span<const AffinePoint<Field>> table,
		size_t wbits,
		size_t npoints,
		const byte* const scalars[],
		const size_t nbits,
		span<Point<Field>> scratch);

	// Scratch size needed for the gather/accumulate staging.
	template<typename Field>
	BLST_HD BLST_INLINE size_t points_mult_wbits_scratch_sizeof(size_t npoints)
	{
		const size_t scratch_sz = SCRATCH_SZ(Point<Field>);
		return sizeof(Point<Field>) * (npoints < scratch_sz ? npoints : scratch_sz);
	}





	/*
	* Pippenger algorithm implementation, fastest option for larger amount
	* of points...
	*
	* Heuristic: chooses wbits from npoints (smaller for few points, larger up to ~npoints-dependent caps).
	*/
	BLST_HD BLST_INLINE size_t pippenger_window_size(size_t npoints)
	{
		size_t wbits;

		// wbits = log2(npoints)
		for (wbits = 0; npoints >>= 1; wbits++);

		// Heuristic: decrement wbits for small npoints, to avoid overhead.
		if (wbits > 12)
			return wbits - 3;
		else if (wbits > 8)
			return wbits - 2;
		else if (wbits > 4)
			return wbits - 1;

		return wbits ? 2 : 1;
	}

	struct WindowExtract
	{
		size_t window_begin;
		size_t missing_overlap;
		size_t extract_begin;
		size_t window_end;
		size_t extract_size;
		limb_t extracted_mask;
		size_t encode_size;

		WindowExtract() = default;
		BLST_HD BLST_INLINE WindowExtract(size_t window_idx, size_t window_size, size_t nbits)
		{
			window_begin = window_idx * window_size;
			missing_overlap = window_begin == 0 ? 1 : 0;
			extract_begin = window_begin > 0 ? window_begin - 1 : 0;
			window_end = std::min<u64>(window_begin + window_size, nbits);
			extract_size = window_end - extract_begin;
			extracted_mask = ((limb_t)1 << extract_size) - 1;
			encode_size = std::min<u64>(nbits - window_begin + 1, window_size);
		}

		BLST_HD BLST_INLINE limb_t get_wval(const byte* scalar) const
		{
			return (get_wval_limb(scalar, extract_begin, extract_size) & extracted_mask) << missing_overlap;
		}

		BLST_HD BLST_INLINE limb_t extract(const byte* scalar) const
		{
			auto ex = get_wval(scalar);
			return booth_encode(ex, encode_size);
		}
	};


	// MSM with Pippenger algorithm.
	// Uses |buckets|=2^{window_size - 1} scratch space 
	// 
	// the basic idea is that we process windows of bits 
	// from the scalars. We use the booth encoding to 
	// pick the window_size values. See mult_window for an 
	// explainer.
	// 
	// We process all points for window_size i, accumulating
	// the results in the buckets. Then we integrate
	// the buckets to get the result for that window_size.
	// 
	// Let d1,...,dn in [0,...,2^w-1] be the window_size values 
	// for each point within the current window_size. We sort 
	// these into buckets b0,...,b_{2^w-1}. 
	// 
	// acc = 0
	// total = 0
	// for digit = 2^w-1 down to 1:
	//     s_digit = sum_{p in b_digit} p // sum of points in bucket
	//     acc = acc + s_digit
	//     total = total + acc
	// 
	template<typename Group, typename AffineGroup, typename GroupZZ>
	BLST_HD BLST_CPP_INLINE void msm_mult_pippenger(
		Group& ret,
		const AffineGroup* const points[],
		size_t npoints,
		const byte* const scalars[],
		size_t nbits,
		GroupZZ buckets[],
		size_t window_size,
		bool print = false);

	template<typename Group>
	BLST_HD BLST_INLINE size_t msm_mult_pippenger_scratch_sizeof_(size_t npoints)
	{
		return sizeof(Group) << (pippenger_window_size(npoints) - 1);
	}

	/*
	 * Batch-affine Pippenger accumulation for short-Weierstrass groups over
	 * Fp or Fp2. The implementation batches distinct affine bucket updates
	 * behind one variable-time field reciprocal and retains projective
	 * collision buckets. Scratch is caller-owned, uninitialized storage.
	 *
	 * BatchSize is deliberately a compile-time parameter. The production
	 * planner can select among explicitly instantiated kernels without
	 * introducing variable-sized inner loops or heap allocation.
	 */
	template<typename Field, size_t BatchSize>
	size_t msm_mult_batch_affine_scratch_sizeof(size_t window_size);

	template<typename Field, size_t BatchSize>
	bool msm_mult_batch_affine(
		Point<Field>& ret,
		const AffinePoint<Field> points[],
		const byte scalars[],
		size_t npoints,
		size_t nbits,
		limb_t scratch[],
		size_t scratch_size,
		size_t window_size);

	//template<typename Group, typename AffineGroup, typename GroupZZ>
	//BLST_HD BLST_INLINE void msm_tile_pippenger(
	//	Group& ret,
	//	const AffineGroup* const points[],
	//	size_t npoints,
	//	const byte* const scalars[],
	//	size_t nbits,
	//	GroupZZ scratch[],
	//	size_t bit0, size_t window)
	//{
	//	size_t wbits, cbits;

	//	if (bit0 + window > nbits)  wbits = nbits - bit0, cbits = wbits + 1;
	//	else                        wbits = cbits = window;
	//	msm_tile_pippenger<Group, AffineGroup, GroupZZ>(ret, points, npoints, scalars, nbits, scratch,
	//		bit0, wbits, cbits);
	//}

	// dispatches to either fixed-table MSM or Pippenger based on npoints
	template<typename Field>
	BLST_HD BLST_INLINE void msm_mult(
		Point<Field>& ret,
		const AffinePoint<Field>* const points[],
		size_t npoints,
		const byte* const scalars[], size_t nbits,
		Pointxyzz<Field> scratch[])
	{

		if (npoints == 1) {
			Point<Field>::from_affine(ret, *points[0]);
			mult_window<Point<Field>, 5>(ret, ret, scalars[0], nbits);
		}
		else if ((npoints * sizeof(AffinePoint<Field>) * 8 * 3) <= SCRATCH_LIMIT &&
			npoints < 32) {
			AffinePoint<Field>* table = (AffinePoint<Field>*)malloc(npoints * sizeof(AffinePoint<Field>) * 8);
			msm_precompute_wbits<Field>(table, 4, points, npoints);
			msm_mult_table<Field>(ret, table, 4, npoints, scalars, nbits, NULL);
			free(table);
		}
		else
		{
			msm_mult_pippenger<Point<Field>, AffinePoint<Field>, Pointxyzz<Field>>(
				ret, points, npoints, scalars, nbits, scratch, 0);
		}
	}


}
