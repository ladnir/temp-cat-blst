/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/scalar/curve/multi_scalar.h"

#include <cstdio>
#include <type_traits>
#if defined(_MSC_VER)
# include <malloc.h>
# define BLST_ALLOCA _alloca
#else
# include <alloca.h>
# define BLST_ALLOCA alloca
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

#ifndef SCRATCH_LIMIT
# ifdef __wasm__
#  define SCRATCH_LIMIT (45 * 1024)
# else
		/* Performance with 144K scratch is within 1-2-3% from optimal */
#  define SCRATCH_LIMIT (144 * 1024)
# endif
#endif
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
			Field* acc = reinterpret_cast<Field*>(block_dst);

			point = *point_iter ? *point_iter++ : point + 1;
			bool_t is_inf = PointT::is_identity(*point);
			vec_select(acc, detail::mont_one<Field>(), &point->Z, sizeof(Field), is_inf);
			++acc;

			for (size_t i = 1; i < delta; ++i, ++acc) {
				point = *point_iter ? *point_iter++ : point + 1;
				bool_t next_inf = PointT::is_identity(*point);
				vec_select(acc, detail::mont_one<Field>(), &point->Z, sizeof(Field), next_inf);
				Field::mul(acc[0], acc[0], acc[-1]);
			}

			--acc;
			Field::reciprocal(acc[0], acc[0]);

			const PointT* p = point;
			const PointT* const* walkback = point_iter - 1;
			size_t remaining = delta ? delta - 1 : 0;
			AffinePoint<Field>* dst_iter = block_dst + remaining;

			for (size_t i = 0; i < remaining; ++i, --acc, --dst_iter) {
				bool_t current_inf = PointT::is_identity(*p);
				Field::mul(acc[-1], acc[-1], acc[0]);
				Field::sqr(ZZ, acc[-1]);
				Field::mul(ZZZ, ZZ, acc[-1]);
				vec_select(&acc[-1], detail::mont_one<Field>(), &p->Z, sizeof(Field), current_inf);
				Field::mul(acc[-1], acc[-1], acc[0]);
				Field::mul(dst_iter->X, p->X, ZZ);
				Field::mul(dst_iter->Y, p->Y, ZZZ);
				vec_czero(dst_iter, sizeof(*dst_iter), current_inf);
				p = (p == *walkback) ? *--walkback : p - 1;
			}

			Field::sqr(ZZ, acc[0]);
			Field::mul(ZZZ, ZZ, acc[0]);
			Field::mul(dst_iter->X, p->X, ZZ);
			Field::mul(dst_iter->Y, p->Y, ZZZ);
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
		vec_copy(&row[0].Z, &Field::one(), sizeof(row[0].Z));
		Point<Field>::doubling(row[1], row[0]);                  /* row[1]=p*(1+1) */
		vec_select(&row[1].Z, &Field::one(), &row[1].Z, sizeof(row[1].Z), inf);
		for (i = 2, j = 1; i < n; i += 2, j++)
			Point<Field>::add(row[i], row[i - 1], point),  /* row[2]=p*(2+1) */
			Point<Field>::doubling(row[i + 1], row[j]),             /* row[3]=p*(2+2) */
			vec_select(&row[i + 1].Z, &Field::one(), &row[i + 1].Z, sizeof(row[i + 1].Z), inf);
	}

	template<typename Field>
	BLST_HD BLST_CPP_INLINE void point_precompute_row_wbits(
		Point<Field> row[],
		size_t wbits,
		const AffinePoint<Field>& point)
	{
		point_precompute_row<Field>(row, (size_t)1 << (wbits - 1), point);
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
	BLST_HD BLST_CPP_INLINE void points_to_affine_row_wbits(
		AffinePoint<Field> dst[],
		Point<Field> src[],
		size_t wbits, size_t npoints)
	{
		if (npoints == 0) {
			return;
		}

		size_t total = npoints << (wbits - 1);
		size_t nwin = (size_t)1 << (wbits - 1);
		size_t i, j;
		Field* acc, ZZ, ZZZ;

		src += total;
		acc = (Field*)src;
		vec_copy(acc++, &Field::one(), sizeof(Field));
		for (i = 0; i < npoints; i++)
			for (j = nwin; --src, --j; acc++)
				Field::mul(acc[0], acc[-1], src->Z);

		--acc; Field::reciprocal(acc[0], acc[0]);

		for (i = 0; i < npoints; i++) {
			vec_copy(dst++, src++, sizeof(AffinePoint<Field>));
			for (j = 1; j < nwin; j++, acc--, src++, dst++) {
				Field::mul(acc[-1], acc[-1], acc[0]);  /* 1/Z        */
				Field::sqr(ZZ, acc[-1]);               /* 1/Z^2      */
				Field::mul(ZZZ, ZZ, acc[-1]);          /* 1/Z^3      */
				Field::mul(acc[-1], src->Z, acc[0]);
				Field::mul(dst->X, src->X, ZZ);        /* X = X'/Z^2 */
				Field::mul(dst->Y, src->Y, ZZZ);       /* Y = Y'/Z^3 */
			}
		}
	}
	template
		BLST_HD BLST_CPP_INLINE void points_to_affine_row_wbits(
		AffinePoint<Fp> dst[],
		Point<Fp> src[],
		size_t wbits, size_t npoints);
	template
		BLST_HD BLST_CPP_INLINE void points_to_affine_row_wbits(
		AffinePoint<Fp2> dst[],
		Point<Fp2> src[],
		size_t wbits, size_t npoints);

	// Orchestrates precompute for all input points in chunks (to bound memory), 
	// calling the two routines above and producing a flat affine table of size 
	// npoints << (wbits-1).
	// flat |points[n]| can be placed at the end of |table[n<<(wbits-1)]| 
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
		size_t table_bytes = 0;
		if (!detail::checked_mul_size(sizeof(AffinePoint<Field>), total, table_bytes)) {
			return;
		}
		size_t nwin = (size_t)1 << (wbits - 1);
		size_t nmin = wbits > 9 ? (size_t)1 : (size_t)1 << (9 - wbits);
		size_t i, top = 0;
		Point<Field>* rows, * row;
		const AffinePoint<Field>* point = NULL;
		size_t stride = ((512 * 1024) / sizeof(AffinePoint<Field>)) >> wbits;
		if (stride == 0) stride = 1;

		while (npoints >= nmin) {
			size_t limit = total - npoints;

			if (stride > ((limit - top) >> wbits)) {
				stride = (limit - top) >> wbits;
				if (stride == 0) break;
			}
			rows = row = (Point<Field>*)(&table[top]);
			for (i = 0; i < stride; i++, row += nwin)
				point = *points ? *points++ : point + 1,
				point_precompute_row<Field>(row, nwin, *point);
			points_to_affine_row_wbits<Field>(&table[top], rows, wbits, stride);
			top += stride << (wbits - 1);
			npoints -= stride;
		}
		if (npoints == 0) {
			return;
		}

		size_t row_count = 0;
		size_t rows_bytes = 0;
		bool whole_rows_fit =
			detail::checked_mul_size(npoints, nwin, row_count) &&
			detail::checked_mul_size(2 * sizeof(AffinePoint<Field>), row_count, rows_bytes) &&
			rows_bytes <= static_cast<size_t>(-1) - sizeof(Field);
		if (whole_rows_fit) {
			rows_bytes += sizeof(Field);
			whole_rows_fit = rows_bytes <= SCRATCH_LIMIT;
		}

		if (whole_rows_fit) {
			rows = row = (Point<Field>*)BLST_ALLOCA(rows_bytes);
			for (i = 0; i < npoints; i++, row += nwin) {
				point = *points ? *points++ : point + 1;
				point_precompute_row<Field>(row, nwin, *point);
			}
			points_to_affine_row_wbits<Field>(&table[top], rows, wbits, npoints);
		} else {
			const Point<Field>* packed_rows[2];

			static_assert(
				SCRATCH_LIMIT >= 2 * sizeof(Point<Field>),
				"SCRATCH_LIMIT must hold at least two projective points");
			stride = SCRATCH_LIMIT / sizeof(Point<Field>);
			stride -= stride % 2;
			if (stride > nwin) {
				stride = nwin;
			}

			packed_rows[0] = row = (Point<Field>*)BLST_ALLOCA(stride * sizeof(Point<Field>));
			packed_rows[1] = nullptr;
			for (i = 0; i < npoints; i++, top += nwin) {
				size_t j, k, n;

				point = *points ? *points++ : point + 1;
				point_precompute_row<Field>(row, stride, *point);
				points_to_affine<Field>(&table[top], packed_rows, stride);
				for (j = stride; j < nwin; j += stride) {
					n = (j + stride) <= nwin ? stride : nwin - j;
					for (k = 0; k < n - 1; k++) {
						Point<Field>::add(row[k], row[stride - 1], table[top + k]);
					}
					if (j == stride) {
						Point<Field>::doubling(row[k], row[stride - 1]);
					} else {
						Point<Field>::add(row[k], row[stride - 1], table[top + k]);
					}
					points_to_affine<Field>(&table[top + j], packed_rows, n);
				}
			}
		}
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
	//   suggested size is min(npoints, SCRATCH_SZ(Point<Field>)).
	template<typename Field>
	BLST_HD BLST_CPP_INLINE void msm_mult_table(
		Point<Field>& ret,
		span<const AffinePoint<Field>> table,
		size_t wbits,
		size_t npoints,
		const byte* const scalars[],
		const size_t nbits,
		span<Point<Field>> scratch)
	{
		if (!detail::valid_wbits(wbits) || npoints == 0 || scalars == nullptr ||
			scalars[0] == nullptr || scratch.empty()) {
			ret = Point<Field>::infinity();
			return;
		}

		// scalar_s will walk through the array of scalar pointers.
		// if they passed in one array (scalars[1]==0), we will just keep using that pointer.
		// if they passed in multiple scalars, we will advance the pointer
		const byte* const* scalar_s = scalars;


		auto table_size = (size_t)1 << (wbits - 1);
		if (table.size() != npoints * table_size)
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
		point_gather_booth_wbits<Field>(scratch[0], cur_table, wbits, wval);
		auto scratch_idx = 1ull;
		cur_table += table_size;

		auto point_idx = 1ull;
		vec_zero(&ret, sizeof(ret));
		while (rem_bits > 0) {
			for (; point_idx < npoints; point_idx++, scratch_idx++, cur_table += table_size) {
				if (scratch_idx == scratch_cap)
				{
					points_accumulate<Field>(ret, scratch.data(), scratch_idx);
					scratch_idx = 0;
				}
				scalar = *scalar_s ? *scalar_s++ : scalar + nbytes;
				wval = get_wval_limb(scalar, rem_bits - 1, cur_wbit + 1) & window_mask;
				wval = booth_encode(wval, wbits);
				point_gather_booth_wbits<Field>(scratch[scratch_idx], cur_table, wbits, wval);
			}
			points_accumulate<Field>(ret, scratch.data(), scratch_idx);

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
				points_accumulate<Field>(ret, scratch.data(), scratch_idx);
				scratch_idx = 0;
			}
			scalar = *scalar_s ? *scalar_s++ : scalar + nbytes;
			wval = (get_wval_limb(scalar, 0, cur_wbit) << 1) & window_mask;
			wval = booth_encode(wval, wbits);
			point_gather_booth_wbits<Field>(scratch[scratch_idx], cur_table, wbits, wval);
		}
		points_accumulate<Field>(ret, scratch.data(), scratch_idx);
	}
	template 		BLST_HD BLST_CPP_INLINE void msm_mult_table(
		Point<Fp>& ret,
		span<const AffinePoint<Fp>> table,
		size_t wbits,
		size_t npoints,
		const byte* const scalars[],
		const size_t nbits,
		span<Point<Fp>> scratch);
	template 		BLST_HD BLST_CPP_INLINE void msm_mult_table(
		Point<Fp2>& ret,
		span<const AffinePoint<Fp2>> table,
		size_t wbits,
		size_t npoints,
		const byte* const scalars[],
		const size_t nbits,
		span<Point<Fp2>> scratch);


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

		booth_idx &= (1 << wbits) - 1;
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
		booth_idx &= (1 << wbits) - 1;
		if (booth_idx--)
		{
#ifndef __CUDA_ARCH__
			vec_prefetch(&buckets[booth_idx], sizeof(buckets[booth_idx]));
#endif
		}
	}

	BLST_HD inline void print_digit(size_t point_idx, limb_t encoding, int wbits)
	{
		size_t digit = encoding & (((limb_t)1 << wbits) - 1);
		size_t sign = (encoding >> wbits) & 1;

		printf("  point[%zu]: digit %zu (%llu), idx=%zu, sign=%zu\n",
			point_idx, digit, (long long unsigned int)encoding, point_idx, sign);
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
		size_t window_size,
		bool print)
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


		// fetch the next window_size value and advance the pointers.
		auto extract_window = [&]() -> std::pair< const AffineGroup*, limb_t> {
			auto encoding = we.extract(scalar);

			if (print)
			{
				auto ex = we.get_wval(scalar);
				size_t digit = encoding & (((limb_t)1 << we.encode_size) - 1);
				size_t sign = (encoding >> we.encode_size) & 1;

				char buffer[257];
				print_bits(buffer, (void*)scalar, std::min<u64>(256, nbits));
				printf("  point[%zu]: %p %s -> digit %zu (%llu), sign=%zu\n",
					i, scalar, buffer, digit, (long long unsigned int)ex, sign);
			}
			auto p = point;
			scalar = *scalars ? *scalars++ : scalar + nbytes;
			point = *points ? *points++ : point + 1;
			return { p, encoding };
			};

		// current window_size value
		auto nxt_wval = extract_window();
		++i;

		for (; i < npoints; i++) {

			auto cur_wval = nxt_wval;

			// prefetch next window_size value
			nxt_wval = extract_window();
			point_prefetch(buckets, nxt_wval.second, we.encode_size);

			// buckets[cur_wval] += *point
			point_bucket(buckets, cur_wval.second, we.encode_size, *cur_wval.first);
		}

		// buckets[nxt_wval] += *point
		point_bucket(buckets, nxt_wval.second, we.encode_size, *nxt_wval.first);

		if (print)
		{

			printf("Pippenger window %llu buckets:\n", (long long unsigned int)window_idx);
			for (size_t b = 0; b < ((size_t)1 << (we.encode_size - 1)); b++)
			{
				Group p;
				msm_to_jacobian(p, buckets[b]);
				u64 hash = hash_point(p);
				printf("  bucket %zu: %llu\n", b, (long long unsigned int)hash);
			}
		}

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
		size_t window_size,
		bool print)
	{
		Group tile;

		window_size = window_size ? window_size : pippenger_window_size(npoints);
		for (size_t i = 0, n = (size_t)1 << (window_size - 1); i < n; ++i) {
			detail::clear_msm_bucket(buckets[i]);
		}
		detail::set_msm_identity(ret);

		if (print)
		{
			auto scalars_ = scalars;
			const byte* scalar = *scalars_++;
			auto nbytes = (nbits + 7) / 8; /* convert |nbits| to bytes */
			for (u64 i = 0; i < npoints; ++i)
			{
				char buffer[257];
				print_bits(buffer, (void*)scalar, std::min<u64>(256, nbits));
				printf("Scalar[%llu]:%p %s\n", (long long unsigned int)i, scalar, buffer);
				scalar = *scalars_ ? *scalars_++ : scalar + nbytes;
			}
		}

		// we process divCeil(nbits+1, window_size) windows.
		// The plus one is for the overlap bit from the booth encoding.
		u64 num_windows = 1 + nbits / window_size;
		for (u64 window_index = num_windows - 1; window_index < num_windows; --window_index)
		{
			// tile = points_tile_pippenger(...)
			points_tile_pippenger(tile, points, npoints, scalars, nbits,
				buckets, window_index, window_size, print);

			if (print)
			{
				u64 hash = hash_point(tile);
				printf("Window %llu result: %llu\n", 
					(long long unsigned int)window_index, 
					(long long unsigned int)hash);
			}
			// ret = (ret + tile) << (window_index ? window_size : 0)
			Group::dadd(ret, ret, tile);
			for (u64 i = 0; i < (window_index ? window_size : 0); i++)
				Group::doubling(ret, ret);
		}

		if (print)
		{
			u64 hash = hash_point(ret);
			printf("Pippenger final result:  %llu\n", (long long unsigned int)hash);
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
			size_t window_size,
			bool print);
	template
		BLST_HD BLST_CPP_INLINE void msm_mult_pippenger(
			Point<Fp2>& ret,
			const AffinePoint<Fp2>* const points[],
			size_t npoints,
			const byte* const scalars[],
			size_t nbits,
			Pointxyzz<Fp2> buckets[],
			size_t window_size,
			bool print);
	template
		BLST_HD BLST_CPP_INLINE void msm_mult_pippenger(
			GT& ret,
			const GT* const points[],
			size_t npoints,
			const byte* const scalars[],
			size_t nbits,
			GT buckets[],
			size_t window_size,
			bool print);

}
