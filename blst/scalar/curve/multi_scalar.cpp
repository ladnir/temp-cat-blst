/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/scalar/curve/multi_scalar.h"

#include <cstdlib>
#include <new>
#include <type_traits>

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
		inline constexpr size_t kPrecomputeFallbackBytes = 8 * 1024;

		//Returns per - block batch size used in batched affine conversion(1536 for Fp, 768 for Fp2).
		template<typename Field>
		constexpr size_t block_stride()
		{
			if constexpr (std::is_same_v<Field, Fp>)
				return 1536;
			else
				return 768;
		}
		//Returns scratch footprint per batch(8192 for Fp, 4096 for Fp2).
		template<typename Field>
		constexpr size_t scratch_stride()
		{
			if constexpr (std::is_same_v<Field, Fp>)
				return 8192;
			else
				return 4096;
		}

		template<typename Field>
		BLST_HD BLST_CPP_INLINE const Field* mont_one()
		{
			return &Field::one();
		}

		template<typename Group>
		BLST_HD BLST_CPP_INLINE void set_msm_identity(Group& out)
		{
			if constexpr (std::is_same_v<Group, GT>) {
				out = GT::identity();
			} else {
				vec_zero(&out, sizeof(out));
			}
		}

		template<typename Group>
		BLST_HD BLST_CPP_INLINE void clear_msm_bucket(Group& out)
		{
			set_msm_identity(out);
		}

		template<typename Field>
		class PrecomputeWorkspace
		{
			using PointT = Point<Field>;

			byte* storage_;

			BLST_HD BLST_FORCEINLINE byte* point_address(size_t index) const
			{
				return storage_ + index * sizeof(PointT);
			}

			BLST_HD BLST_FORCEINLINE byte* prefix_address(
				size_t point_count, size_t index) const
			{
				return storage_ + point_count * sizeof(PointT) + index * sizeof(Field);
			}

		public:
			explicit PrecomputeWorkspace(void* storage)
				: storage_(static_cast<byte*>(storage))
			{
				static_assert(std::is_trivially_default_constructible_v<PointT>);
				static_assert(std::is_trivially_default_constructible_v<Field>);
				static_assert(std::is_trivially_destructible_v<PointT>);
				static_assert(std::is_trivially_destructible_v<Field>);
				static_assert(alignof(AffinePoint<Field>) >= alignof(PointT));
				static_assert(alignof(PointT) >= alignof(Field));
				static_assert(sizeof(PointT) % alignof(Field) == 0);
			}

			BLST_HD BLST_FORCEINLINE PointT& start_point(size_t index) const
			{
				return *::new (point_address(index)) PointT;
			}

			BLST_HD BLST_FORCEINLINE PointT& point(size_t index) const
			{
				return *std::launder(reinterpret_cast<PointT*>(point_address(index)));
			}

			BLST_HD BLST_FORCEINLINE Field& start_prefix(
				size_t point_count, size_t index) const
			{
				return *::new (prefix_address(point_count, index)) Field;
			}

			BLST_HD BLST_FORCEINLINE Field& prefix(
				size_t point_count, size_t index) const
			{
				return *std::launder(
					reinterpret_cast<Field*>(prefix_address(point_count, index)));
			}
		};

		template<typename Field>
		BLST_HD BLST_INLINE bool precompute_workspace_size(
			size_t point_count, size_t& bytes)
		{
			size_t point_bytes = 0;
			size_t prefix_count = 0;
			size_t prefix_bytes = 0;
			if (!checked_mul_size(point_count, sizeof(Point<Field>), point_bytes) ||
				point_count == static_cast<size_t>(-1)) {
				return false;
			}
			prefix_count = point_count + 1;
			if (!checked_mul_size(prefix_count, sizeof(Field), prefix_bytes) ||
				point_bytes > static_cast<size_t>(-1) - prefix_bytes) {
				return false;
			}
			bytes = point_bytes + prefix_bytes;
			return true;
		}

		template<typename Field>
		BLST_HD BLST_FORCEINLINE AffinePoint<Field>& start_affine(
			byte* storage, size_t index)
		{
			static_assert(
				std::is_trivially_default_constructible_v<AffinePoint<Field>>);
			static_assert(std::is_trivially_destructible_v<AffinePoint<Field>>);
			return *::new (storage + index * sizeof(AffinePoint<Field>))
				AffinePoint<Field>;
		}

	} // namespace detail

	// Converts an array of (possibly sparse) Jacobian points 
	// to affine using batched inversion of Z over sized blocks; 
	// preserves infinities.
	template<typename Field>
	BLST_HD BLST_CPP_INLINE void points_to_affine(AffinePoint<Field> dst[],
		const Point<Field>* const points[],
		size_t npoints)
	{
		using PointT = Point<Field>;
		constexpr size_t stride = detail::block_stride<Field>();

		Field ZZ;
		Field ZZZ;
		const PointT* point = nullptr;
		AffinePoint<Field>* out = dst;
		const PointT* const* point_iter = points;

		while (npoints) {
			size_t delta = stride < npoints ? stride : npoints;
			AffinePoint<Field>* block_dst = out;

			point = *point_iter ? *point_iter++ : point + 1;
			bool_t is_inf = PointT::is_identity(*point);
			vec_select(&block_dst[0].x, detail::mont_one<Field>(),
				&point->z, sizeof(Field), is_inf);

			for (size_t i = 1; i < delta; ++i) {
				point = *point_iter ? *point_iter++ : point + 1;
				bool_t next_inf = PointT::is_identity(*point);
				vec_select(&block_dst[i].x, detail::mont_one<Field>(),
					&point->z, sizeof(Field), next_inf);
				Field::mul(block_dst[i].x, block_dst[i].x, block_dst[i - 1].x);
			}

			Field::reciprocal(block_dst[delta - 1].x, block_dst[delta - 1].x);

			const PointT* p = point;
			const PointT* const* walkback = point_iter - 1;
			size_t remaining = delta ? delta - 1 : 0;
			AffinePoint<Field>* dst_iter = block_dst + remaining;

			for (size_t i = remaining; i != 0; --i, --dst_iter) {
				bool_t current_inf = PointT::is_identity(*p);
				Field& current_acc = block_dst[i].x;
				Field& previous_acc = block_dst[i - 1].x;
				Field::mul(previous_acc, previous_acc, current_acc);
				Field::sqr(ZZ, previous_acc);
				Field::mul(ZZZ, ZZ, previous_acc);
				vec_select(&previous_acc, detail::mont_one<Field>(),
					&p->z, sizeof(Field), current_inf);
				Field::mul(previous_acc, previous_acc, current_acc);
				Field::mul(dst_iter->x, p->x, ZZ);
				Field::mul(dst_iter->y, p->y, ZZZ);
				vec_czero(dst_iter, sizeof(*dst_iter), current_inf);
				p = (p == *walkback) ? *--walkback : p - 1;
			}

			Field::sqr(ZZ, block_dst[0].x);
			Field::mul(ZZZ, ZZ, block_dst[0].x);
			Field::mul(dst_iter->x, p->x, ZZ);
			Field::mul(dst_iter->y, p->y, ZZZ);
			vec_czero(dst_iter, sizeof(*dst_iter), PointT::is_identity(*p));

			out = block_dst + delta;
			npoints -= delta;
		}

	}
	template
		BLST_HD BLST_CPP_INLINE void points_to_affine(AffinePoint<Fp> dst[],
			const Point<Fp>* const points[],
			size_t npoints);
	template
		BLST_HD BLST_CPP_INLINE void points_to_affine(AffinePoint<Fp2> dst[],
			const Point<Fp2>* const points[],
			size_t npoints);


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

		/* The intermediate infinity points are encoded as [0, 0, 1]. */

	// For one affine point, builds a row of n consecutive multiples used by
	// Booth / gather (e.g., P·1, P·2, P·3, …) in Jacobian. n is even and at
	// least two for every precompute caller.
	template<typename Field>
	BLST_HD BLST_CPP_INLINE void point_precompute_row(
		Point<Field> row[],
		size_t n,
		const AffinePoint<Field>& point)
	{
		size_t i, j;
		bool_t inf = vec_is_zero(&point, sizeof(point));
		/* row[-1] is implicit infinity */
		vec_copy(&row[0], &point, sizeof(point));           /* row[0]=p*1     */
		vec_copy(&row[0].z, &Field::one(), sizeof(row[0].z));
		Point<Field>::doubling(row[1], row[0]);                  /* row[1]=p*(1+1) */
		vec_select(&row[1].z, &Field::one(), &row[1].z, sizeof(row[1].z), inf);
		for (i = 2, j = 1; i < n; i += 2, j++)
			Point<Field>::add(row[i], row[i - 1], point),  /* row[2]=p*(2+1) */
			Point<Field>::doubling(row[i + 1], row[j]),             /* row[3]=p*(2+2) */
			vec_select(&row[i + 1].z, &Field::one(), &row[i + 1].z, sizeof(row[i + 1].z), inf);
	}

	template<typename Field>
	BLST_HD BLST_CPP_INLINE void point_precompute_row_wbits(
		Point<Field> row[],
		size_t wbits,
		const AffinePoint<Field>& point)
	{
		point_precompute_row<Field>(row, (size_t)1 << (wbits - 1), point);
	}

	template<typename Field>
	BLST_HD BLST_CPP_INLINE void point_precompute_row_workspace(
		const detail::PrecomputeWorkspace<Field>& workspace,
		size_t row_begin,
		size_t n,
		const AffinePoint<Field>& point)
	{
		size_t i, j;
		bool_t inf = vec_is_zero(&point, sizeof(point));
		auto& row0 = workspace.start_point(row_begin);
		vec_copy(&row0, &point, sizeof(point));
		vec_copy(&row0.z, &Field::one(), sizeof(row0.z));

		auto& row1 = workspace.start_point(row_begin + 1);
		Point<Field>::doubling(row1, row0);
		vec_select(&row1.z, &Field::one(), &row1.z, sizeof(row1.z), inf);
		for (i = 2, j = 1; i < n; i += 2, ++j) {
			auto& even = workspace.start_point(row_begin + i);
			Point<Field>::add(
				even, workspace.point(row_begin + i - 1), point);
			auto& odd = workspace.start_point(row_begin + i + 1);
			Point<Field>::doubling(odd, workspace.point(row_begin + j));
			vec_select(&odd.z, &Field::one(), &odd.z, sizeof(odd.z), inf);
		}
	}
	template
	BLST_HD BLST_CPP_INLINE void point_precompute_row_wbits(
		Point<Fp> row[],
		size_t wbits,
		const AffinePoint<Fp>& point);

	template
	BLST_HD BLST_CPP_INLINE void point_precompute_row_wbits(
		Point<Fp2> row[],
		size_t wbits,
		const AffinePoint<Fp2>& point);

	// Converts a block of those rows from Jacobian to affine with one batched 
	// inversion chain; writes dst table laid out as rows.
	template<typename Field>
	BLST_HD BLST_CPP_INLINE void points_to_affine_row_workspace(
		byte* dst_storage,
		const detail::PrecomputeWorkspace<Field>& workspace,
		size_t nwin,
		size_t npoints)
	{
		if (npoints == 0) {
			return;
		}

		const size_t total = npoints * nwin;
		size_t src_index = 0;
		size_t acc_index = 0;
		Field ZZ, ZZZ;

		vec_copy(&workspace.start_prefix(total, acc_index++),
			&Field::one(), sizeof(Field));
		for (size_t row = npoints; row != 0; --row) {
			for (size_t j = nwin; j > 1; --j) {
				const size_t index = (row - 1) * nwin + (j - 1);
				const auto& src = workspace.point(index);
				auto& acc = workspace.start_prefix(total, acc_index);
				Field::mul(acc, workspace.prefix(total, acc_index - 1), src.z);
				++acc_index;
			}
		}

		BLST_ASSERT(acc_index == npoints * (nwin - 1) + 1);
		--acc_index;
		Field::reciprocal(
			workspace.prefix(total, acc_index),
			workspace.prefix(total, acc_index));

		src_index = 0;
		for (size_t i = 0; i < npoints; ++i) {
			const auto& first = workspace.point(src_index);
			if (src_index < 2) {
				const Field x = first.x;
				const Field y = first.y;
				auto& dst = detail::start_affine<Field>(dst_storage, src_index);
				dst.x = x;
				dst.y = y;
			} else {
				auto& dst = detail::start_affine<Field>(dst_storage, src_index);
				vec_copy(&dst, &first, sizeof(dst));
			}
			++src_index;

			for (size_t j = 1; j < nwin; ++j, ++src_index, --acc_index) {
				const auto& src = workspace.point(src_index);
				auto& previous = workspace.prefix(total, acc_index - 1);
				auto& current = workspace.prefix(total, acc_index);
				Field::mul(previous, previous, current); /* 1/Z */
				Field::sqr(ZZ, previous);
				Field::mul(ZZZ, ZZ, previous);
				Field::mul(previous, src.z, current);

				if (src_index < 2) {
					Field x, y;
					Field::mul(x, src.x, ZZ);
					Field::mul(y, src.y, ZZZ);
					auto& dst = detail::start_affine<Field>(dst_storage, src_index);
					dst.x = x;
					dst.y = y;
				} else {
					auto& dst = detail::start_affine<Field>(dst_storage, src_index);
					Field::mul(dst.x, src.x, ZZ);
					Field::mul(dst.y, src.y, ZZZ);
				}
			}
		}
	}

#if defined(_MSC_VER)
# define BLST_MSM_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define BLST_MSM_NOINLINE __attribute__((noinline))
#else
# define BLST_MSM_NOINLINE
#endif

	template<typename Field>
	BLST_HD static BLST_MSM_NOINLINE void msm_precompute_wbits_fallback(
		AffinePoint<Field> table[],
		size_t top,
		size_t nwin,
		const AffinePoint<Field>* point,
		const AffinePoint<Field>* const points[],
		size_t npoints)
	{
		alignas(Point<Field>)
		byte fallback_storage[detail::kPrecomputeFallbackBytes];
		detail::PrecomputeWorkspace<Field> workspace(fallback_storage);
		static_assert(
			detail::kPrecomputeFallbackBytes >=
				2 * (sizeof(Point<Field>) + sizeof(Field)) + sizeof(Field),
			"precompute fallback must hold at least two projective points");
		size_t stride = (detail::kPrecomputeFallbackBytes - sizeof(Field)) /
			(sizeof(Point<Field>) + sizeof(Field));
		stride -= stride % 2;
		if (stride > nwin)
			stride = nwin;

		for (size_t i = 0; i < npoints; ++i, top += nwin) {
			point = *points ? *points++ : point + 1;
			point_precompute_row_workspace<Field>(workspace, 0, stride, *point);
			const Point<Field>* packed_rows[2] = {
				&workspace.point(0), nullptr
			};
			points_to_affine<Field>(&table[top], packed_rows, stride);

			for (size_t j = stride; j < nwin; j += stride) {
				const size_t n = (j + stride) <= nwin ? stride : nwin - j;
				const Point<Field>& previous = workspace.point(stride - 1);
				size_t k = 0;
				for (; k + 1 < n; ++k) {
					auto& out = workspace.start_point(k);
					Point<Field>::add(out, previous, table[top + k]);
				}
				auto& last = k == stride - 1
					? workspace.point(k) : workspace.start_point(k);
				if (j == stride) {
					Point<Field>::doubling(last, previous);
				} else {
					Point<Field>::add(last, previous, table[top + k]);
				}
				points_to_affine<Field>(&table[top + j], packed_rows, n);
			}
		}
	}

#undef BLST_MSM_NOINLINE

	// Orchestrates bounded table reuse, optional heap workspace, and the cold
	// fixed-stack fallback. A packed input array may occupy the table tail.
	template<typename Field>
	BLST_HD BLST_CPP_INLINE void msm_precompute_wbits(
		AffinePoint<Field> table[],
		size_t wbits,
		const AffinePoint<Field>* const points[],
		size_t npoints)
	{
		if (!detail::valid_wbits(wbits) || npoints == 0 || table == nullptr ||
			points == nullptr || points[0] == nullptr)
			return;

		size_t total = 0;
		if (!detail::checked_shl_size(npoints, wbits - 1, total)) {
			return;
		}
		if (total > static_cast<size_t>(-1) / sizeof(AffinePoint<Field>)) {
			return;
		}
		size_t nwin = (size_t)1 << (wbits - 1);
		size_t nmin = wbits > 9 ? (size_t)1 : (size_t)1 << (9 - wbits);
		size_t top = 0;
		const AffinePoint<Field>* point = nullptr;
		byte* const table_storage = reinterpret_cast<byte*>(table);
		size_t stride = ((512 * 1024) / sizeof(AffinePoint<Field>)) >> wbits;
		if (stride == 0) stride = 1;

		while (npoints >= nmin) {
			size_t limit = total - npoints;
			const size_t available_bytes =
				(limit - top) * sizeof(AffinePoint<Field>);
			const size_t per_row_bytes =
				nwin * (sizeof(Point<Field>) + sizeof(Field));
			const size_t max_stride = available_bytes > sizeof(Field)
				? (available_bytes - sizeof(Field)) / per_row_bytes : 0;
			if (stride > max_stride) {
				stride = max_stride;
				if (stride == 0)
					break;
			}

			detail::PrecomputeWorkspace<Field> workspace(
				table_storage + top * sizeof(AffinePoint<Field>));
			for (size_t i = 0; i < stride; ++i) {
				point = *points ? *points++ : point + 1;
				point_precompute_row_workspace<Field>(
					workspace, i * nwin, nwin, *point);
			}
			points_to_affine_row_workspace<Field>(
				table_storage + top * sizeof(AffinePoint<Field>),
				workspace, nwin, stride);
			top += stride * nwin;
			npoints -= stride;
		}
		if (npoints == 0) {
			return;
		}

		size_t row_count = 0;
		size_t workspace_bytes = 0;
		const bool whole_rows_fit =
			detail::checked_mul_size(npoints, nwin, row_count) &&
			detail::precompute_workspace_size<Field>(row_count, workspace_bytes) &&
			workspace_bytes <= BLST_SCALAR_MSM_SCRATCH_LIMIT;

		if (whole_rows_fit) {
#ifndef __CUDA_ARCH__
			void* allocation = std::malloc(workspace_bytes);
			if (allocation != nullptr) {
				detail::PrecomputeWorkspace<Field> workspace(allocation);
				for (size_t i = 0; i < npoints; ++i) {
					point = *points ? *points++ : point + 1;
					point_precompute_row_workspace<Field>(
						workspace, i * nwin, nwin, *point);
				}
				points_to_affine_row_workspace<Field>(
					table_storage + top * sizeof(AffinePoint<Field>),
					workspace, nwin, npoints);
				std::free(allocation);
				return;
			}
#endif
		}

		msm_precompute_wbits_fallback<Field>(
			table, top, nwin, point, points, npoints);
	}
	template 
		BLST_HD BLST_CPP_INLINE void msm_precompute_wbits(
		AffinePoint<Fp> table[],
		size_t wbits,
		const AffinePoint<Fp>* const points[],
			size_t npoints);
	template
		BLST_HD BLST_CPP_INLINE void msm_precompute_wbits(
			AffinePoint<Fp2> table[],
			size_t wbits,
			const AffinePoint<Fp2>* const points[],
			size_t npoints);


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
	//   suggested size is min(npoints, fixed_table_scratch_capacity<Field>()).
	template<typename Field>
	BLST_HD BLST_CPP_INLINE void msm_mult_table(
		Point<Field>& ret,
		span<const AffinePoint<Field>> table,
		size_t wbits,
		size_t npoints,
		const byte* const scalars[],
		const size_t nbits,
		detail::BatchAffineWorkspace<Field> scratch)
	{
		if (npoints == 0 || nbits == 0) {
			ret = Point<Field>::infinity();
			return;
		}
		if (nbits > 256 || !detail::valid_wbits(wbits) || scalars == nullptr ||
			scalars[0] == nullptr || scratch.empty()) {
			ret = Point<Field>::infinity();
			return;
		}

		// scalar_s will walk through the array of scalar pointers.
		// if they passed in one array (scalars[1]==0), we will just keep using that pointer.
		// if they passed in multiple scalars, we will advance the pointer
		const byte* const* scalar_s = scalars;


		auto table_size = (size_t)1 << (wbits - 1);
		size_t expected_table_size = 0;
		if (!detail::checked_mul_size(npoints, table_size, expected_table_size) ||
			table.size() != expected_table_size)
		{
			ret = Point<Field>::infinity();
			return; // or error
		}

		size_t scratch_cap = scratch.size();


		auto nbytes = (nbits + 7) / 8; /* convert |nbits| to bytes */
		auto scalar = *scalar_s++;

		/* top excess bits modulo target window_size size */
		auto cur_wbit = nbits % wbits; /* yes, it may be zero */
		auto window_mask = ((limb_t)1 << (cur_wbit + 1)) - 1;

		auto rem_bits = nbits - cur_wbit;

		// the MSB window_size is allowed to be partial, so we handle it first.
		// In addition, the LSB window_size is special in that it does not have
		// the overlap bit (bit -1). Moreover, we need to load the first 
		// element into the scratch before entering the main loop. Because
		// of all that, we unroll the first iteration here. 
		auto no_overlap = is_zero(rem_bits);
		auto overlap = no_overlap ^ 1;

		// if there is an overlap, we load the overlap bit from the next window_size.
		// otherwise we push in a 0 for the non-existent overlap bit.
		auto wval = (get_wval_limb(scalar, rem_bits - overlap, cur_wbit + overlap) << no_overlap) & window_mask;
		wval = booth_encode(wval, wbits);

		const AffinePoint<Field>* cur_table = table.data();
		// scratch[0] = sign(wval) * table[0][index(wval)]
		point_gather_booth_wbits<Field>(scratch.start(0), cur_table, wbits, wval);
		auto scratch_idx = 1ull;
		cur_table += table_size;

		auto point_idx = 1ull;
		vec_zero(&ret, sizeof(ret));
		while (rem_bits > 0) {
			for (; point_idx < npoints; point_idx++, scratch_idx++, cur_table += table_size) {
				if (scratch_idx == scratch_cap)
				{
					detail::points_accumulate<Field>(ret, scratch, scratch_idx);
					scratch_idx = 0;
				}
				scalar = *scalar_s ? *scalar_s++ : scalar + nbytes;
				wval = get_wval_limb(scalar, rem_bits - 1, cur_wbit + 1) & window_mask;
				wval = booth_encode(wval, wbits);
				point_gather_booth_wbits<Field>(scratch.start(scratch_idx), cur_table, wbits, wval);
			}
			detail::points_accumulate<Field>(ret, scratch, scratch_idx);

			for (size_t j = 0; j < wbits; j++)
				Point<Field>::doubling(ret, ret);

			// reset the window_size (always full now), and restart point/scalar indices
			cur_wbit = wbits;
			window_mask = ((limb_t)1 << (cur_wbit + 1)) - 1;
			rem_bits -= cur_wbit;
			point_idx = 0;
			scratch_idx = 0;
			cur_table = table.data();
			scalar_s = scalars;
		}

		// handle the LSB window_size, which does not have the overlap bit.
		for (; point_idx < npoints; point_idx++, scratch_idx++, cur_table += table_size) {
			if (scratch_idx == scratch_cap)
			{
				detail::points_accumulate<Field>(ret, scratch, scratch_idx);
				scratch_idx = 0;
			}
			scalar = *scalar_s ? *scalar_s++ : scalar + nbytes;
			wval = (get_wval_limb(scalar, 0, cur_wbit) << 1) & window_mask;
			wval = booth_encode(wval, wbits);
			point_gather_booth_wbits<Field>(scratch.start(scratch_idx), cur_table, wbits, wval);
		}
		detail::points_accumulate<Field>(ret, scratch, scratch_idx);
	}
	template 		BLST_HD BLST_CPP_INLINE void msm_mult_table(
		Point<Fp>& ret,
		span<const AffinePoint<Fp>> table,
		size_t wbits,
		size_t npoints,
		const byte* const scalars[],
		const size_t nbits,
		detail::BatchAffineWorkspace<Fp> scratch);
	template 		BLST_HD BLST_CPP_INLINE void msm_mult_table(
		Point<Fp2>& ret,
		span<const AffinePoint<Fp2>> table,
		size_t wbits,
		size_t npoints,
		const byte* const scalars[],
		const size_t nbits,
		detail::BatchAffineWorkspace<Fp2> scratch);


	// High->low running sum over buckets(xyzz form) : 
	// computes Σ d·S_d from pre-summed buckets and writes 
	// the window’s contribution as Jacobian.
	template<typename Group, typename GroupZZ>
	BLST_HD BLST_CPP_INLINE void point_integrate_buckets(
		Group& out,
		GroupZZ buckets[],
		size_t wbits)
	{
		GroupZZ ret, acc;
		size_t n = (size_t)1 << wbits;

		/* Calculate sum of x[i-1]*i for i=1 through 1<<|wbits|. */
		vec_copy(&acc, &buckets[--n], sizeof(acc));
		vec_copy(&ret, &buckets[n], sizeof(ret));
		detail::clear_msm_bucket(buckets[n]);


		// computes Σ d·S_d
		while (n--) {
			GroupZZ::dadd(acc, acc, buckets[n]);
			GroupZZ::dadd(ret, ret, acc);
			detail::clear_msm_bucket(buckets[n]);
		}
		msm_to_jacobian(out, ret);
	}

	// Adds an affine point into the bucket addressed by Booth-encoded digit (handles sign).
	template<typename AffineGroup, typename GroupZZ>
	BLST_HD BLST_CPP_INLINE void point_bucket(GroupZZ buckets[], limb_t booth_idx,
		size_t wbits, const AffineGroup& p)
	{
		bool_t booth_sign = (booth_idx >> wbits) & 1;

		booth_idx &= (static_cast<limb_t>(1) << wbits) - 1;
		if (booth_idx--)
		{
			msm_add_or_sub(buckets[booth_idx], buckets[booth_idx],
				p, booth_sign);
		}
	}

	template<typename GroupZZ>
	BLST_HD BLST_CPP_INLINE void point_prefetch(const GroupZZ buckets[], limb_t booth_idx,
		size_t wbits)
	{
		booth_idx &= (static_cast<limb_t>(1) << wbits) - 1;
		if (booth_idx--)
		{
#ifndef __CUDA_ARCH__
			vec_prefetch(&buckets[booth_idx], sizeof(buckets[booth_idx]));
#endif
		}
	}

	// Processes one window_size across all points: extracts Booth digits for that 
	// bit range, fills buckets, then integrates buckets to produce the window_size
	// partial ret.
	// 
	// acc = 0
	// total = 0
	// for digit = 2^w-1 down to 1:
	//     s_digit = sum_{p in b_digit} p // sum of points in bucket
	//     acc = acc + s_digit
	//     total = total + acc
	//
	// @ret : output window_size partial (Jacobian)
	// @points : input points array
	// @npoints : number of points
	// @scalars : input scalars array, can either be an array of points or 
	// an array of pointers to points. For the former, the array must be
	// the null terminated, e.g. scalars = { ptr, null }, ptr = { scalar0, scalar1, ... }
	// @nbits : number of bits in each scalar
	// @buckets : workspace of size 2^{wbits-1} Pointxyzz. This will hold the sum of the 
	// points with each digit.
	// @window_idx : the index of the current window_size (0-based). The LSB window_size is window_size 0.
	// @window_size : the size of a full sized window_size in bits. 
	template<typename Group, typename AffineGroup, typename GroupZZ>
	BLST_HD BLST_CPP_INLINE void points_tile_pippenger(
		Group& ret,
		const AffineGroup* const points[],
		size_t npoints,
		const byte* const scalars[],
		size_t nbits,
		GroupZZ buckets[],
		size_t window_idx,
		size_t window_size)
	{
		BLST_ASSERT(npoints > 1);

		// we have exactly ceil((nbits+1)/window_size) windows. When 
		// bits is a multiple of window_size, the MSB window_size is "empty"
		// except for the overlap bit from the next window_size.
		BLST_ASSERT(window_idx < 1 + nbits / window_size);

		WindowExtract we(window_idx, window_size, nbits);

		size_t i = 0;
		const byte* scalar = *scalars++;
		const AffineGroup* point = *points++;

		auto nbytes = (nbits + 7) / 8; /* convert |nbits| to bytes */


		// Fetch the current digit. Advance only when another input remains so
		// an exact n-entry descriptor array never requires a sentinel slot.
		auto extract_window = [&]() -> std::pair<const AffineGroup*, limb_t> {
			auto encoding = we.extract(scalar);
			return { point, encoding };
			};
		auto advance_input = [&]() {
			scalar = *scalars ? *scalars++ : scalar + nbytes;
			point = *points ? *points++ : point + 1;
			};

		// current window_size value
		auto nxt_wval = extract_window();
		++i;

		for (; i < npoints; i++) {

			auto cur_wval = nxt_wval;

			// prefetch next window_size value
			advance_input();
			nxt_wval = extract_window();
			point_prefetch(buckets, nxt_wval.second, we.encode_size);

			// buckets[cur_wval] += *point
			point_bucket(buckets, cur_wval.second, we.encode_size, *cur_wval.first);
		}

		// buckets[nxt_wval] += *point
		point_bucket(buckets, nxt_wval.second, we.encode_size, *nxt_wval.first);

		point_integrate_buckets(ret, buckets, we.encode_size - 1);

	}


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
		size_t window_size)
	{
		Group tile;

		window_size = window_size ? window_size : pippenger_window_size(npoints);
		BLST_ASSERT(detail::valid_pippenger_window(window_size));
		BLST_ASSERT(npoints > 1 && nbits != 0 && nbits <= 256);
		for (size_t i = 0, n = (size_t)1 << (window_size - 1); i < n; ++i) {
			detail::clear_msm_bucket(buckets[i]);
		}
		detail::set_msm_identity(ret);

		// we process divCeil(nbits+1, window_size) windows.
		// The plus one is for the overlap bit from the booth encoding.
		u64 num_windows = 1 + nbits / window_size;
		for (u64 window_index = num_windows - 1; window_index < num_windows; --window_index)
		{
			// tile = points_tile_pippenger(...)
			points_tile_pippenger(tile, points, npoints, scalars, nbits,
				buckets, window_index, window_size);
			// ret = (ret + tile) << (window_index ? window_size : 0)
			Group::dadd(ret, ret, tile);
			for (u64 i = 0; i < (window_index ? window_size : 0); i++)
				Group::doubling(ret, ret);
		}

	}
	template
		BLST_HD BLST_CPP_INLINE void msm_mult_pippenger(
			Point<Fp>& ret,
			const AffinePoint<Fp>* const points[],
			size_t npoints,
			const byte* const scalars[],
			size_t nbits,
			Pointxyzz<Fp> buckets[],
			size_t window_size);
	template
		BLST_HD BLST_CPP_INLINE void msm_mult_pippenger(
			Point<Fp2>& ret,
			const AffinePoint<Fp2>* const points[],
			size_t npoints,
			const byte* const scalars[],
			size_t nbits,
			Pointxyzz<Fp2> buckets[],
			size_t window_size);
	template
		BLST_HD BLST_CPP_INLINE void msm_mult_pippenger(
			GT& ret,
			const GT* const points[],
			size_t npoints,
			const byte* const scalars[],
			size_t nbits,
			GT buckets[],
			size_t window_size);

}
