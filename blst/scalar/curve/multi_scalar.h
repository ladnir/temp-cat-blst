/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/scalar/curve/point.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/bulk_addition.h"

#include <algorithm>
#include <new>
#include <type_traits>

#ifdef __CUDACC__
#include <cuda/std/span>
#else
#include  <span>
#endif
namespace blst::scalar
{
	/** Batch-normalize secret-derived points with a fixed inversion chain.
	 * Null descriptors after the first entry denote contiguous input. */
	template<typename Field>
	void points_to_affine_ct(AffinePoint<Field> dst[],
		const Point<Field>* const points[], size_t npoints);

#ifndef BLST_SCALAR_MSM_SCRATCH_LIMIT
# ifdef __wasm__
#  define BLST_SCALAR_MSM_SCRATCH_LIMIT (45 * 1024)
# else
#  define BLST_SCALAR_MSM_SCRATCH_LIMIT (144 * 1024)
# endif
#endif

#ifdef __CUDACC__
	template<typename T>
	using span = cuda::std::span<T>;
#else
	template<typename T>
	using span = std::span<T>;
#endif

	namespace detail
	{
		inline constexpr size_t kPippengerMaxWindow = 24;

		/*
		 * Pippenger's C boundary receives caller-owned limb storage rather than
		 * a live C++ array. This pointer-only view indexes each bucket by byte
		 * offset, so separately started objects never rely on array pointer
		 * arithmetic. For native C++ callers, from_objects() views an existing
		 * array without restarting its elements. The no-initializer placement
		 * new in from_uninitialized() starts trivial object lifetimes without
		 * writing the scratch buffer.
		 */
		template<typename T>
		class PippengerBucketWorkspace
		{
		public:
			BLST_HD constexpr PippengerBucketWorkspace() : storage_(nullptr)
		#ifndef NDEBUG
				, size_(0)
		#endif
			{}

			BLST_HD static PippengerBucketWorkspace from_objects(
				T* storage, size_t size)
			{
				return PippengerBucketWorkspace(storage, size);
			}

			BLST_HD static PippengerBucketWorkspace from_uninitialized(
				void* storage, size_t size)
			{
				static_assert(std::is_trivially_default_constructible_v<T>);
				static_assert(std::is_trivially_destructible_v<T>);
				PippengerBucketWorkspace out(storage, size);
				for (size_t i = 0; i < size; ++i)
					::new (out.address(i)) T;
				return out;
			}

			BLST_HD BLST_INLINE T& operator[](size_t index) const
			{
				return *std::launder(reinterpret_cast<T*>(address(index)));
			}

		private:
			BLST_HD PippengerBucketWorkspace(void* storage, size_t size)
				: storage_(static_cast<byte*>(storage))
		#ifndef NDEBUG
				, size_(size)
		#endif
			{
				(void)size;
				BLST_ASSERT(storage_ != nullptr);
				BLST_ASSERT(
					reinterpret_cast<uintptr_t>(storage_) % alignof(T) == 0);
			}

			BLST_HD BLST_INLINE void* address(size_t index) const
			{
		#ifndef NDEBUG
				BLST_ASSERT(index < size_);
		#endif
				return storage_ + index * sizeof(T);
			}

			byte* storage_;
		#ifndef NDEBUG
			size_t size_;
		#endif
		};

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

		BLST_HD BLST_INLINE constexpr bool valid_pippenger_window(size_t window)
		{
			/* get_wval_limb supports at most 25 extracted bits. */
			return window >= 1 && window <= kPippengerMaxWindow;
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
		// The signed-window digit is public in this explicitly variable-time path.
		// Keep the branch here rather than charging the constant-time GT API's
		// full conditional-copy cost to every MSM addition.
		if (subtract) {
			Fp12::conjugate(neg.data);
		}
		GT::add(out, p1, neg);
	}



#ifdef BLST_CUDA
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
#endif

	/*
	 * Reusable fixed-table MSM. Ordinary operands are prevalidated subgroup
	 * points and scalar bits are public. A row contains 2^(wbits-1)
	 * consecutive affine multiples for wbits in [2, 14]. Compact pointer
	 * descriptors require a non-null first entry; each later null entry means
	 * contiguous continuation. Exact n-entry arrays need no sentinel.
	 *
	 * Precomputation may overlap a packed input array only when that input is
	 * stored in the final npoints table entries. It reuses free output storage
	 * first, then attempts one bounded nonthrowing allocation, and otherwise
	 * uses a fixed small stack workspace. Multiplication requires matching
	 * caller-owned scratch and accepts nbits in [1, 256]; zero points or bits
	 * produce the identity without reading input storage.
	 */
	template<typename Field>
	BLST_HD BLST_INLINE constexpr size_t fixed_table_scratch_capacity()
	{
		static_assert(std::is_same_v<Field, Fp> || std::is_same_v<Field, Fp2>);
		return std::is_same_v<Field, Fp> ? 8192 : 4096;
	}

		/* The intermediate infinity points are encoded as [0, 0, 1]. */

	// For one affine point, builds a row of length 2 ^ (wbits - 1) of 
	// consecutive multiples used by Booth / gather(e.g., P·1, P·2, P·3, …) 
	// in Jacobian.
	template<typename Field>
	BLST_HD BLST_CPP_INLINE void point_precompute_row_wbits(
		Point<Field> row[],
		size_t wbits,
		const AffinePoint<Field>& point);

	// Produces a flat affine table with npoints * 2^(wbits-1) entries.
	// A packed points array may occupy the final npoints entries of table.
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
		detail::BatchAffineNode<Field>& node,
		const AffinePoint<Field> row[],
		size_t wbits, limb_t booth_encoding)
	{
		bool_t booth_sign = (booth_encoding >> wbits) & 1;
		bool_t idx_is_zero;
		const AffinePoint<Field> infinity = { };

		auto booth_idx = booth_encoding & (((limb_t)1 << wbits) - 1);
		idx_is_zero = is_zero(booth_idx);
		booth_idx -= 1 ^ idx_is_zero;
		vec_select(&node.point, &infinity, &row[booth_idx], sizeof(row[0]), idx_is_zero);
		AffinePoint<Field>::cneg(node.point, booth_sign);
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
	// @nbits number of meaningful bits in each scalar, at most 256
	// @scratch scratch space for staging gathers before accumulation;
	//   suggested size is min(npoints, fixed_table_scratch_capacity<Field>()).
	template<typename Field>
	BLST_HD BLST_CPP_INLINE void msm_mult_table(
		Point<Field>& ret,
		span<const AffinePoint<Field>> table,
		size_t wbits,
		size_t npoints,
		const byte* const scalars[],
		const size_t nbits,
		detail::BatchAffineWorkspace<Field> scratch);

	// Scratch size needed for the gather/accumulate staging.
	template<typename Field>
	BLST_HD BLST_INLINE size_t points_mult_wbits_scratch_sizeof(size_t npoints)
	{
		using Node = detail::BatchAffineNode<Field>;
		const size_t scratch_nodes = fixed_table_scratch_capacity<Field>();
		return sizeof(Node) * (npoints < scratch_nodes ? npoints : scratch_nodes);
	}





	/*
	 * Variable-time Pippenger for public scalar data and prevalidated subgroup
	 * points. Scalars contain at most 256 meaningful little-endian bits.
	 * Explicit windows are in [1, 24]; the default heuristic is capped to that
	 * range. Scratch contains 2^(window-1) buckets. Compact point and scalar
	 * descriptors follow the same exact-array/contiguous-null convention as
	 * fixed-table MSM.
	 */
	BLST_HD BLST_INLINE size_t pippenger_window_size(size_t npoints)
	{
		size_t wbits;

		// wbits = log2(npoints)
		for (wbits = 0; npoints >>= 1; wbits++);

		// Heuristic: decrement wbits for small npoints, to avoid overhead.
		if (wbits > 12)
			wbits -= 3;
		else if (wbits > 8)
			wbits -= 2;
		else if (wbits > 4)
			wbits -= 1;
		else
			return wbits ? 2 : 1;

		return wbits < detail::kPippengerMaxWindow
			? wbits : detail::kPippengerMaxWindow;
	}

	struct WindowExtract
	{
		size_t window_begin = 0;
		size_t missing_overlap = 0;
		size_t extract_begin = 0;
		size_t window_end = 0;
		size_t extract_size = 0;
		limb_t extracted_mask = 0;
		size_t encode_size = 0;

		WindowExtract() = default;
		BLST_HD BLST_INLINE WindowExtract(size_t window_idx, size_t window_size, size_t nbits)
		{
			BLST_ASSERT(detail::valid_pippenger_window(window_size));
			BLST_ASSERT(nbits != 0 && nbits <= 256);
			BLST_ASSERT(window_idx <= nbits / window_size);
			initialize(window_idx * window_size, window_size, nbits);
		}

		BLST_HD BLST_INLINE static WindowExtract from_bit_offset(
			size_t bit0, size_t window_size, size_t nbits)
		{
			BLST_ASSERT(detail::valid_pippenger_window(window_size));
			BLST_ASSERT(nbits != 0 && nbits <= 256);
			BLST_ASSERT(bit0 < nbits);
			WindowExtract out;
			out.initialize(bit0, window_size, nbits);
			return out;
		}

	private:
		BLST_HD BLST_INLINE void initialize(size_t bit0, size_t window_size, size_t nbits)
		{
			window_begin = bit0;
			missing_overlap = window_begin == 0 ? 1 : 0;
			extract_begin = window_begin > 0 ? window_begin - 1 : 0;
			window_end = std::min(window_begin + window_size, nbits);
			extract_size = window_end - extract_begin;
			extracted_mask = ((limb_t)1 << extract_size) - 1;
			encode_size = std::min(nbits - window_begin + 1, window_size);
		}

	public:

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
		detail::PippengerBucketWorkspace<GroupZZ> buckets,
		size_t window_size);

	/* Native C++ convenience overload for an already-live bucket array. */
	template<typename Group, typename AffineGroup, typename GroupZZ>
	BLST_HD BLST_INLINE void msm_mult_pippenger(
		Group& ret,
		const AffineGroup* const points[],
		size_t npoints,
		const byte* const scalars[],
		size_t nbits,
		GroupZZ buckets[],
		size_t window_size)
	{
		const size_t window = window_size
			? window_size : pippenger_window_size(npoints);
		auto workspace = detail::PippengerBucketWorkspace<GroupZZ>::from_objects(
			buckets, static_cast<size_t>(1) << (window - 1));
		msm_mult_pippenger<Group, AffineGroup, GroupZZ>(
			ret, points, npoints, scalars, nbits, workspace, window_size);
	}

	template<typename Group>
	BLST_HD BLST_INLINE size_t msm_mult_pippenger_scratch_sizeof_(size_t npoints)
	{
		if (npoints == 0)
			return 0;
		size_t bytes = 0;
		const size_t window = pippenger_window_size(npoints);
		return detail::checked_shl_size(sizeof(Group), window - 1, bytes)
			? bytes : 0;
	}

	template<typename Group>
	BLST_HD BLST_INLINE size_t msm_mult_pippenger_window_scratch_sizeof_(size_t window)
	{
		if (!detail::valid_pippenger_window(window))
			return 0;
		size_t bytes = 0;
		return detail::checked_shl_size(sizeof(Group), window - 1, bytes)
			? bytes : 0;
	}

	/**
	 * Batch-affine Pippenger accumulation for short-Weierstrass groups over
	 * Fp or Fp2. The implementation batches distinct affine bucket updates
	 * behind one variable-time field reciprocal and retains projective
	 * collision buckets. Scratch is caller-owned, uninitialized storage.
	 *
	 * BatchSize is deliberately a compile-time parameter. The production
	 * planner can select among explicitly instantiated kernels without
	 * introducing variable-sized inner loops or heap allocation.
	 * Windows are in [1, 24], matching WindowExtract's 25-bit maximum including
	 * overlap. Inputs are contiguous, nonidentity subgroup points and public
	 * little-endian scalars with nbits in [1, 256]. Scratch is disjoint,
	 * limb-aligned raw storage of at least the queried byte count. On failure,
	 * the output is unchanged; empty input succeeds and writes infinity.
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

}
