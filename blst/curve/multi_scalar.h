/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/curve/point.h"
#include "blst/curve/bulk_addition.h"

#ifdef __CUDACC__
#include <cuda/std/span>
#else
#include  <span>
#endif
namespace blst
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
		BLST_HD BLST_INLINE const Field* mont_one()
		{
			return &Field::one();
		}

		template<typename Field>
		BLST_HD BLST_INLINE bool_t is_infinity(const Point<Field>& p)
		{
			return vec_is_zero(&p.Z, sizeof(p.Z));
		}

		template<typename Field>
		BLST_HD BLST_INLINE bool_t is_infinity(const AffinePoint<Field>& p)
		{
			return vec_is_zero(&p, sizeof(AffinePoint<Field>));
		}

	} // namespace detail

	// Converts an array of (possibly sparse) Jacobian points 
	// to affine using batched inversion of Z over sized blocks; 
	// preserves infinities.
	template<typename Field>
	BLST_HD BLST_INLINE void points_to_affine(AffinePoint<Field> dst[],
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
			bool_t is_inf = detail::is_infinity(*point);
			vec_select(acc, detail::mont_one<Field>(), &point->Z, sizeof(Field), is_inf);
			++acc;

			for (size_t i = 1; i < delta; ++i, ++acc) {
				point = *point_iter ? *point_iter++ : point + 1;
				bool_t next_inf = detail::is_infinity(*point);
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
				bool_t current_inf = detail::is_infinity(*p);
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
			vec_czero(dst_iter, sizeof(*dst_iter), detail::is_infinity(*p));

			out = block_dst + delta;
			npoints -= delta;
		}
	}

	//void blst_p1s_to_affine(G1AffinePoint dst[],
	//	const G1Point* const points[],
	//	size_t npoints)
	//{
	//	points_to_affine<Fp>(dst, points, npoints);
	//}

	//void blst_p2s_to_affine(G2AffinePoint dst[],
	//	const G2Point* const points[],
	//	size_t npoints)
	//{
	//	points_to_affine<Fp2>(dst, points, npoints);
	//}
	/* X = X'/Z^2 */\





















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

#define PRECOMPUTE_WBITS_IMPL(prefix, ptype, bits, field, one) \
void ptype##_precompute_row_wbits(ptype row[], size_t wbits, \
                                         const ptype##_affine *point) \
{ \
    size_t i, j, n = (size_t)1 << (wbits-1); \
    bool_t inf = vec_is_zero(point, sizeof(*point)); \
                                          /* row[-1] is implicit infinity */\
    vec_copy(&row[0], point, sizeof(*point));           /* row[0]=p*1     */\
    vec_copy(&row[0].Z, one, sizeof(row[0].Z));                             \
    ptype##_double(&row[1],  &row[0]);                  /* row[1]=p*(1+1) */\
    vec_select(&row[1].Z, one, &row[1].Z, sizeof(row[1].Z), inf);           \
    for (i = 2, j = 1; i < n; i += 2, j++) \
        ptype##_add_affine(&row[i], &row[i-1], point),  /* row[2]=p*(2+1) */\
        ptype##_double(&row[i+1], &row[j]),             /* row[3]=p*(2+2) */\
        vec_select(&row[i+1].Z, one, &row[i+1].Z, sizeof(row[i+1].Z), inf); \
}                                                       /* row[4] ...     */\
\
void ptype##s_to_affine_row_wbits(ptype##_affine dst[], ptype src[], \
                                         size_t wbits, size_t npoints) \
{ \
    size_t total = npoints << (wbits-1); \
    size_t nwin = (size_t)1 << (wbits-1); \
    size_t i, j; \
    vec##bits *acc, ZZ, ZZZ; \
\
    src += total; \
    acc = (vec##bits *)src; \
    vec_copy(acc++, one, sizeof(vec##bits)); \
    for (i = 0; i < npoints; i++) \
        for (j = nwin; --src, --j; acc++)    \
            mul_##field(acc[0], acc[-1], src->Z); \
\
    --acc; reciprocal_##field(acc[0], acc[0]); \
\
    for (i = 0; i < npoints; i++) { \
        vec_copy(dst++, src++, sizeof(ptype##_affine)); \
        for (j = 1; j < nwin; j++, acc--, src++, dst++) { \
            mul_##field(acc[-1], acc[-1], acc[0]);  /* 1/Z        */\
            sqr_##field(ZZ, acc[-1]);               /* 1/Z^2      */\
            mul_##field(ZZZ, ZZ, acc[-1]);          /* 1/Z^3      */\
            mul_##field(acc[-1], src->Z, acc[0]);                   \
            mul_##field(dst->X, src->X, ZZ);        /* X = X'/Z^2 */\
            mul_##field(dst->Y, src->Y, ZZZ);       /* Y = Y'/Z^3 */\
        } \
    } \
} \
\
/* flat |points[n]| can be placed at the end of |table[n<<(wbits-1)]| */\
void ptype##s_precompute_wbits(ptype##_affine table[], size_t wbits, \
                                      const ptype##_affine *const points[], \
                                      size_t npoints) \
{ \
    size_t total = npoints << (wbits-1); \
    size_t nwin = (size_t)1 << (wbits-1); \
    size_t nmin = wbits>9 ? (size_t)1: (size_t)1 << (9-wbits); \
    size_t i, top = 0; \
    ptype *rows, *row; \
    const ptype##_affine *point = NULL; \
    size_t stride = ((512*1024)/sizeof(ptype##_affine)) >> wbits; \
    if (stride == 0) stride = 1; \
\
    while (npoints >= nmin) { \
        size_t limit = total - npoints; \
\
        if (top + (stride << wbits) > limit) { \
            stride = (limit - top) >> wbits;   \
            if (stride == 0) break;            \
        } \
        rows = row = (ptype *)(&table[top]); \
        for (i = 0; i < stride; i++, row += nwin) \
            point = *points ? *points++ : point+1, \
            ptype##_precompute_row_wbits(row, wbits, point); \
        ptype##s_to_affine_row_wbits(&table[top], rows, wbits, stride); \
        top += stride << (wbits-1); \
        npoints -= stride; \
    } \
    rows = row = (ptype*)alloca(2*sizeof(ptype##_affine) * npoints * nwin); \
    for (i = 0; i < npoints; i++, row += nwin) \
        point = *points ? *points++ : point+1, \
        ptype##_precompute_row_wbits(row, wbits, point); \
    ptype##s_to_affine_row_wbits(&table[top], rows, wbits, npoints); \
} \
\
size_t prefix##s_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints) \
{ return (sizeof(ptype##_affine)*npoints) << (wbits-1); } \
void prefix##s_mult_wbits_precompute(ptype##_affine table[], size_t wbits, \
                                     const ptype##_affine *const points[], \
                                     size_t npoints) \
{ ptype##s_precompute_wbits(table, wbits, points, npoints); }

#define POINTS_MULT_WBITS_IMPL(prefix, ptype, bits, field, one) \
void ptype##_gather_booth_wbits(ptype *p, const ptype##_affine row[], \
                                       size_t wbits, limb_t booth_idx) \
{ \
    bool_t booth_sign = (booth_idx >> wbits) & 1; \
    bool_t idx_is_zero; \
    const ptype##_affine infinity = { 0 }; \
\
    booth_idx &= ((limb_t)1 << wbits) - 1; \
    idx_is_zero = is_zero(booth_idx); \
    booth_idx -= 1 ^ idx_is_zero; \
    vec_select(p, &infinity, &row[booth_idx], sizeof(row[0]), idx_is_zero); \
    ptype##_cneg(p, booth_sign); \
} \
\
void ptype##s_mult_wbits(ptype *ret, const ptype##_affine table[], \
                                size_t wbits, size_t npoints, \
                                const byte *const scalars[], size_t nbits, \
                                ptype scratch[]) \
{ \
    limb_t wmask, wval; \
    size_t i, j, z, nbytes, window, nwin = (size_t)1 << (wbits-1); \
    const byte *scalar, *const *scalar_s = scalars; \
    const ptype##_affine *row = table; \
\
    size_t scratch_sz = SCRATCH_SZ(ptype); \
    if (scratch == NULL) { \
        scratch_sz /= 4; /* limit to 288K */ \
        scratch_sz = scratch_sz < npoints ? scratch_sz : npoints; \
        scratch = (ptype*)alloca(sizeof(ptype) * scratch_sz); \
    } \
\
    nbytes = (nbits + 7)/8; /* convert |nbits| to bytes */ \
    scalar = *scalar_s++; \
\
    /* top excess bits modulo target window size */ \
    window = nbits % wbits; /* yes, it may be zero */ \
    wmask = ((limb_t)1 << (window + 1)) - 1; \
\
    nbits -= window; \
    z = is_zero(nbits); \
    wval = (get_wval_limb(scalar, nbits - (z^1), window + (z^1)) << z) & wmask; \
    wval = booth_encode(wval, wbits); \
    ptype##_gather_booth_wbits(&scratch[0], row, wbits, wval); \
    row += nwin; \
\
    i = 1; vec_zero(ret, sizeof(*ret)); \
    while (nbits > 0) { \
        for (j = i; i < npoints; i++, j++, row += nwin) { \
            if (j == scratch_sz) \
                ptype##s_accumulate(ret, scratch, j), j = 0; \
            scalar = *scalar_s ? *scalar_s++ : scalar+nbytes; \
            wval = get_wval_limb(scalar, nbits - 1, window + 1) & wmask; \
            wval = booth_encode(wval, wbits); \
            ptype##_gather_booth_wbits(&scratch[j], row, wbits, wval); \
        } \
        ptype##s_accumulate(ret, scratch, j); \
\
        for (j = 0; j < wbits; j++) \
            ptype##_double(ret, ret); \
\
        window = wbits; \
        wmask = ((limb_t)1 << (window + 1)) - 1; \
        nbits -= window; \
        i = 0; row = table; scalar_s = scalars; \
    } \
\
    for (j = i; i < npoints; i++, j++, row += nwin) { \
        if (j == scratch_sz) \
            ptype##s_accumulate(ret, scratch, j), j = 0; \
        scalar = *scalar_s ? *scalar_s++ : scalar+nbytes; \
        wval = (get_wval_limb(scalar, 0, window) << 1) & wmask; \
        wval = booth_encode(wval, wbits); \
        ptype##_gather_booth_wbits(&scratch[j], row, wbits, wval); \
    } \
    ptype##s_accumulate(ret, scratch, j); \
} \
\
size_t prefix##s_mult_wbits_scratch_sizeof(size_t npoints) \
{ \
    const size_t scratch_sz = SCRATCH_SZ(ptype); \
    return sizeof(ptype) * (npoints < scratch_sz ? npoints : scratch_sz); \
} \
void prefix##s_mult_wbits(ptype *ret, const ptype##_affine table[], \
                          size_t wbits, size_t npoints, \
                          const byte *const scalars[], size_t nbits, \
                          ptype scratch[]) \
{ ptype##s_mult_wbits(ret, table, wbits, npoints, scalars, nbits, scratch); }

		//PRECOMPUTE_WBITS_IMPL(blst_p1, G1Point, 384, fp, get_BLST_Rx().p)
		//POINTS_MULT_WBITS_IMPL(blst_p1, G1Point, 384, fp, get_BLST_Rx().p)



		//PRECOMPUTE_WBITS_IMPL(blst_p2, G2Point, 384x, fp2, get_BLST_Rx().p2)
		//POINTS_MULT_WBITS_IMPL(blst_p2, G2Point, 384x, fp2, get_BLST_Rx().p2)



















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
	void point_precompute_row_wbits(
		Point<Field> row[], 
		size_t wbits,
		const AffinePoint<Field>& point)
	{
		size_t i, j, n = (size_t)1 << (wbits - 1);
		bool_t inf = vec_is_zero(&point, sizeof(point));
		/* row[-1] is implicit infinity */
		vec_copy(&row[0], &point, sizeof(point));           /* row[0]=p*1     */
		vec_copy(&row[0].Z, &Field::one(), sizeof(row[0].Z));
		Point<Field>::doubling(row[1], row[0]);                  /* row[1]=p*(1+1) */
		vec_select(&row[1].Z, &Field::one(), &row[1].Z, sizeof(row[1].Z), inf);
		for (i = 2, j = 1; i < n; i += 2, j++)
			Point<Field>::add_affine(row[i], row[i - 1], point),  /* row[2]=p*(2+1) */
			Point<Field>::doubling(row[i + 1], row[j]),             /* row[3]=p*(2+2) */
			vec_select(&row[i + 1].Z, &Field::one(), &row[i + 1].Z, sizeof(row[i + 1].Z), inf);
	}                                                       /* row[4] ...     */

	// Converts a block of those rows from Jacobian to affine with one batched 
	// inversion chain; writes dst table laid out as rows.
	template<typename Field>
	void points_to_affine_row_wbits(
		AffinePoint<Field> dst[], 
		Point<Field> src[],
		size_t wbits, size_t npoints)
	{
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

	// Orchestrates precompute for all input points in chunks (to bound memory), 
	// calling the two routines above and producing a flat affine table of size 
	// npoints << (wbits-1).
	// flat |points[n]| can be placed at the end of |table[n<<(wbits-1)]| 
	template<typename Field>
	void msm_precompute_wbits(
		AffinePoint<Field> table[], 
		size_t wbits,
		const AffinePoint<Field>* const points[],
		size_t npoints)
	{
		size_t total = npoints << (wbits - 1);
		size_t nwin = (size_t)1 << (wbits - 1);
		size_t nmin = wbits > 9 ? (size_t)1 : (size_t)1 << (9 - wbits);
		size_t i, top = 0;
		Point<Field>* rows, * row;
		const AffinePoint<Field>* point = NULL;
		size_t stride = ((512 * 1024) / sizeof(AffinePoint<Field>)) >> wbits;
		if (stride == 0) stride = 1;

		while (npoints >= nmin) {
			size_t limit = total - npoints;

			if (top + (stride << wbits) > limit) {
				stride = (limit - top) >> wbits;
				if (stride == 0) break;
			}
			rows = row = (Point<Field>*)(&table[top]);
			for (i = 0; i < stride; i++, row += nwin)
				point = *points ? *points++ : point + 1,
				point_precompute_row_wbits<Field>(row, wbits, *point);
			points_to_affine_row_wbits<Field>(&table[top], rows, wbits, stride);
			top += stride << (wbits - 1);
			npoints -= stride;
		}
		rows = row = (Point<Field>*)malloc(2 * sizeof(AffinePoint<Field>) * npoints * nwin);
		for (i = 0; i < npoints; i++, row += nwin)
		{
			point = *points ? *points++ : point + 1;
			point_precompute_row_wbits<Field>(row, wbits, *point);
		}
		points_to_affine_row_wbits<Field>(&table[top], rows, wbits, npoints);
		free(rows);
	}

	template<typename Field>
	size_t points_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints)
	{
		return (sizeof(AffinePoint<Field>) * npoints) << (wbits - 1);
	}

	//template<typename Field>
	//void points_mult_wbits_precompute(AffinePoint<Field> table[], size_t wbits,
	//	const AffinePoint<Field>* const points[],
	//	size_t npoints)
	//{
	//	points_precompute_wbits<Field>(table, wbits, points, npoints);
	//}





	// Given a Booth-encoded digit (with sign bit), selects the appropriate
	// affine entry from a row, applies sign, or infinity if zero.
	template<typename Field>
	void point_gather_booth_wbits(
		Point<Field>& p,
		const AffinePoint<Field> row[],
		size_t wbits, limb_t booth_encoding)
	{
		bool_t booth_sign = (booth_encoding >> wbits) & 1;
		bool_t idx_is_zero;
		const AffinePoint<Field> infinity = { 0 };

		auto booth_idx = booth_encoding & ((limb_t)1 << wbits) - 1;
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
	// @wbits window size in bits (between 2 and 14)
	// @npoints number of points (and scalars)
	// @scalars array of npoints pointers to scalars, each scalar is
	//   represented in byte array of length (nbits+7)/8.
	// @nbits number of bits in each scalar
	// @scratch scratch space for staging gathers before accumulation;
	//   suggested size is min(npoints, SCRATCH_SZ(Point<Field>)).
	template<typename Field>
	void msm_mult_table(
		Point<Field>& ret, 
		span<const AffinePoint<Field>> table,
		size_t wbits,
		size_t npoints,
		const byte* const scalars[], 
		size_t nbits,
		span<Point<Field>> scratch)
	{
		limb_t wmask, wval;
		size_t i, j, z, nbytes, window, nwin = (size_t)1 << (wbits - 1);
		const byte* scalar, * const* scalar_s = scalars;

		if (table.size() != npoints * nwin)
		{
			ret = Point<Field>::infinity();
			return; // or error
		}
		const AffinePoint<Field>* row = table.data();

		size_t scratch_sz = scratch.size();
		//size_t scratch_sz = SCRATCH_SZ(Point<Field>);
		//if (scratch == NULL) {
		//	scratch_sz /= 4; /* limit to 288K */
		//	scratch_sz = scratch_sz < npoints ? scratch_sz : npoints;
		//	scratch = (Point<Field>*)alloca(sizeof(Point<Field>) * scratch_sz);
		//}

		nbytes = (nbits + 7) / 8; /* convert |nbits| to bytes */
		scalar = *scalar_s++;

		/* top excess bits modulo target window size */
		window = nbits % wbits; /* yes, it may be zero */
		wmask = ((limb_t)1 << (window + 1)) - 1;

		nbits -= window;
		z = is_zero(nbits);
		wval = (get_wval_limb(scalar, nbits - (z ^ 1), window + (z ^ 1)) << z) & wmask;
		wval = booth_encode(wval, wbits);
		point_gather_booth_wbits<Field>(scratch[0], row, wbits, wval);
		row += nwin;

		i = 1; vec_zero(&ret, sizeof(ret));
		while (nbits > 0) {
			for (j = i; i < npoints; i++, j++, row += nwin) {
				if (j == scratch_sz)
					points_accumulate<Field>(ret, scratch.data(), j), j = 0;
				scalar = *scalar_s ? *scalar_s++ : scalar + nbytes;
				wval = get_wval_limb(scalar, nbits - 1, window + 1) & wmask;
				wval = booth_encode(wval, wbits);
				point_gather_booth_wbits<Field>(scratch[j], row, wbits, wval);
			}
			points_accumulate<Field>(ret, scratch.data(), j);

			for (j = 0; j < wbits; j++)
				Point<Field>::doubling(ret, ret);

			window = wbits;
			wmask = ((limb_t)1 << (window + 1)) - 1;
			nbits -= window;
			i = 0; row = table.data(); scalar_s = scalars;
		}

		for (j = i; i < npoints; i++, j++, row += nwin) {
			if (j == scratch_sz)
			{
				points_accumulate<Field>(ret, scratch.data(), j), j = 0;
			}
			scalar = *scalar_s ? *scalar_s++ : scalar + nbytes;
			wval = (get_wval_limb(scalar, 0, window) << 1) & wmask;
			wval = booth_encode(wval, wbits);
			point_gather_booth_wbits<Field>(scratch[j], row, wbits, wval);
		}
		points_accumulate<Field>(ret, scratch.data(), j);
	}

	// Scratch size needed for the gather/accumulate staging.
	template<typename Field>
	size_t points_mult_wbits_scratch_sizeof(size_t npoints)
	{
		const size_t scratch_sz = SCRATCH_SZ(Point<Field>);
		return sizeof(Point<Field>) * (npoints < scratch_sz ? npoints : scratch_sz);
	}
	//void prefix##s_mult_wbits(ptype *ret, const ptype##_affine table[], 
	//                          size_t wbits, size_t npoints, 
	//                          const byte *const scalars[], size_t nbits, 
	//                          ptype scratch[]) 
	//{ ptype##s_mult_wbits(ret, table, wbits, npoints, scalars, nbits, scratch); }
	//













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


	// High->low running sum over buckets(xyzz form) : 
	// computes Σ d·S_d from pre-summed buckets and writes 
	// the window’s contribution as Jacobian.
	template<typename Field>
	BLST_HD BLST_INLINE void point_integrate_buckets(
		Point<Field>& out, 
		Pointxyzz<Field> buckets[],
		size_t wbits)
	{
		Pointxyzz<Field> ret, acc;
		size_t n = (size_t)1 << wbits;

		/* Calculate sum of x[i-1]*i for i=1 through 1<<|wbits|. */
		vec_copy(&acc, &buckets[--n], sizeof(acc));
		vec_copy(&ret, &buckets[n], sizeof(ret));
		vec_zero(&buckets[n], sizeof(buckets[n]));


		// computes Σ d·S_d
		while (n--) {
			Pointxyzz<Field>::dadd(acc, acc, buckets[n]);
			Pointxyzz<Field>::dadd(ret, ret, acc);
			vec_zero(&buckets[n], sizeof(buckets[n]));
		}
		Pointxyzz<Field>::to_Jacobian(out, ret);
	}


	// Adds an affine point into the bucket addressed by Booth-encoded digit (handles sign).
	template<typename Field>
	BLST_HD BLST_INLINE void point_bucket(Pointxyzz<Field> buckets[], limb_t booth_idx,
		size_t wbits, const AffinePoint<Field>& p)
	{
		bool_t booth_sign = (booth_idx >> wbits) & 1;

		booth_idx &= (1 << wbits) - 1;
		if (booth_idx--)
			Pointxyzz<Field>::dadd_affine(buckets[booth_idx], buckets[booth_idx],
				p, booth_sign);
	}

	template<typename Field>
	BLST_HD BLST_INLINE void point_prefetch(const Pointxyzz<Field> buckets[], limb_t booth_idx,
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

	// Processes one window across all points: extracts Booth digits for that 
	// bit range, fills buckets, then integrates buckets to produce the window
	// partial ret.
	// 
	// acc = 0
	// total = 0
	// for digit = 2^w-1 down to 1:
	//     s_digit = sum_{p in b_digit} p // sum of points in bucket
	//     acc = acc + s_digit
	//     total = total + acc
	//
	// @ret : output window partial (Jacobian)
	// @points : input points array
	// @npoints : number of points
	// @scalars : input scalars array, can either be an array of points or 
	// an array of pointers to points. For the former, the array must be
	// the null terminated, e.g. scalars = { ptr, null }, ptr = { scalar0, scalar1, ... }
	// @nbits : number of bits in each scalar
	// @buckets : workspace of size 2^{wbits-1} Pointxyzz. This will hold the sum of the 
	// points with each digit.
	// @bit0 : index of the least significant bit in the window.
	// @wbits : window size in bits
	// @cbits : wbits + 1 or wbits depending on the call (Booth encoding)
	template<typename Field>
	BLST_HD BLST_INLINE void points_tile_pippenger(
		Point<Field>& ret,
		const AffinePoint<Field>* const points[],
		size_t npoints,
		const byte* const scalars[], 
		size_t nbits,
		Pointxyzz<Field> buckets[],
		size_t bit0, 
		size_t wbits, 
		size_t cbits)
	{
		limb_t wmask, wval, wnxt;
		size_t i, z, nbytes;
		const byte* scalar = *scalars++;
		const AffinePoint<Field>* point = *points++;

		nbytes = (nbits + 7) / 8; /* convert |nbits| to bytes */
		wmask = ((limb_t)1 << (wbits + 1)) - 1;

		// if this is the bottom window, we can't index into -1, so we
		// special case it. Otherwise we decrement the start bit by one
		// to get an overlap for Booth encoding.
		z = is_zero(bit0);
		bit0 -= z ^ 1; wbits += z ^ 1;
		wval = (get_wval_limb(scalar, bit0, wbits) << z) & wmask;
		wval = booth_encode(wval, cbits);
		scalar = *scalars ? *scalars++ : scalar + nbytes;
		wnxt = (get_wval_limb(scalar, bit0, wbits) << z) & wmask;
		wnxt = booth_encode(wnxt, cbits);
		npoints--;  /* account for prefetch */

		// buckets[wval] += *point
		point_bucket(buckets, wval, cbits, *point);

		for (i = 1; i < npoints; i++) {
			wval = wnxt;
			scalar = *scalars ? *scalars++ : scalar + nbytes;
			wnxt = (get_wval_limb(scalar, bit0, wbits) << z) & wmask;
			wnxt = booth_encode(wnxt, cbits);
			point_prefetch(buckets, wnxt, cbits);
			point = *points ? *points++ : point + 1;

			// buckets[wval] += *point
			point_bucket(buckets, wval, cbits, *point);
		}
		point = *points ? *points++ : point + 1;

		// buckets[wnxt] += *point
		point_bucket(buckets, wnxt, cbits, *point);
		point_integrate_buckets(ret, buckets, cbits - 1);
	}

	// MSM with Pippenger algorithm.
	// Uses |buckets|=2^{window - 1} scratch space 
	// 
	// the basic idea is that we process windows of bits 
	// from the scalars. We use the booth encoding to 
	// pick the window values. See mult_window for an 
	// explainer.
	// 
	// We process all points for window i, accumulating
	// the results in the buckets. Then we integrate
	// the buckets to get the result for that window.
	// 
	// Let d1,...,dn in [0,...,2^w-1] be the window values 
	// for each point within the current window. We sort 
	// these into buckets b0,...,b_{2^w-1}. 
	// 
	// acc = 0
	// total = 0
	// for digit = 2^w-1 down to 1:
	//     s_digit = sum_{p in b_digit} p // sum of points in bucket
	//     acc = acc + s_digit
	//     total = total + acc
	// 
	template<typename Field>
	BLST_HD BLST_INLINE void msm_mult_pippenger(
		Point<Field>& ret,
		const AffinePoint<Field>* const points[],
		size_t npoints,
		const byte* const scalars[], 
		size_t nbits,
		Pointxyzz<Field> buckets[], 
		size_t window)
	{
		size_t i, wbits, cbits, bit0 = nbits;
		Point<Field> tile;

		window = window ? window : pippenger_window_size(npoints);
		vec_zero(buckets, sizeof(buckets[0]) << (window - 1));
		vec_zero(&ret, sizeof(ret));

		/* top excess bits modulo target window size */
		wbits = nbits % window; /* yes, it may be zero */
		cbits = wbits + 1;
		while (bit0 -= wbits) {
			points_tile_pippenger(tile, points, npoints, scalars, nbits,
				buckets, bit0, wbits, cbits);
			Point<Field>::dadd(ret, ret, tile, NULL);
			for (i = 0; i < window; i++)
				Point<Field>::doubling(ret, ret);
			cbits = wbits = window;
		}
		points_tile_pippenger(tile, points, npoints, scalars, nbits,
			buckets, 0, wbits, cbits);
		Point<Field>::dadd(ret, ret, tile, NULL);
	}

	template<typename Field>
	BLST_HD BLST_INLINE size_t msm_mult_pippenger_scratch_sizeof(size_t npoints)
	{
		return sizeof(Pointxyzz<Field>) << (pippenger_window_size(npoints) - 1);
	}
	
	template<typename Field>
	BLST_HD BLST_INLINE void msm_tile_pippenger(
		Point<Field>& ret,
		const AffinePoint<Field>* const points[],
		size_t npoints,
		const byte* const scalars[], 
		size_t nbits,
		Pointxyzz<Field> scratch[],
		size_t bit0, size_t window)
	{
		size_t wbits, cbits;

		if (bit0 + window > nbits)  wbits = nbits - bit0, cbits = wbits + 1;
		else                        wbits = cbits = window;
		msm_tile_pippenger<Field>(ret, points, npoints, scalars, nbits, scratch,
			bit0, wbits, cbits);
	}

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
			msm_mult_pippenger<Field>(ret, points, npoints, scalars, nbits, scratch, 0);
		}
	}


}
