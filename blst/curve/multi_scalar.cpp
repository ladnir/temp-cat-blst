/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/curve/point.h"


namespace blst
{


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

		template<typename Field>
		constexpr size_t block_stride()
		{
			if constexpr (std::is_same_v<Field, Fp>)
				return 1536;
			else
				return 768;
		}

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
		BLST_HD BLST_INLINE bool_t is_infinity(const Point_affine<Field>& p)
		{
			return vec_is_zero(&p, sizeof(Point_affine<Field>));
		}

	} // namespace detail

	template<typename Field>
	void points_to_affine(Point_affine<Field> dst[],
		const Point<Field>* const points[],
		size_t npoints)
	{
		using PointT = Point<Field>;
		constexpr size_t stride = detail::block_stride<Field>();

		Field ZZ;
		Field ZZZ;
		const PointT* point = nullptr;
		Point_affine<Field>* out = dst;
		const PointT* const* point_iter = points;

		while (npoints) {
			size_t delta = stride < npoints ? stride : npoints;
			Point_affine<Field>* block_dst = out;
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
			Point_affine<Field>* dst_iter = block_dst + remaining;

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

	void blst_p1s_to_affine(POINTonE1_affine dst[],
		const POINTonE1* const points[],
		size_t npoints)
	{
		points_to_affine<Fp>(dst, points, npoints);
	}

	void blst_p2s_to_affine(POINTonE2_affine dst[],
		const POINTonE2* const points[],
		size_t npoints)
	{
		points_to_affine<Fp2>(dst, points, npoints);
	}
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

#define SCRATCH_SZ(ptype) (sizeof(ptype)==sizeof(POINTonE1) ? 8192 : 4096)

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
                                       size_t wbits, limb booth_idx) \
{ \
    bool_t booth_sign = (booth_idx >> wbits) & 1; \
    bool_t idx_is_zero; \
    const ptype##_affine infinity = { 0 }; \
\
    booth_idx &= ((limb)1 << wbits) - 1; \
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
    limb wmask, wval; \
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
    wmask = ((limb)1 << (window + 1)) - 1; \
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
        wmask = ((limb)1 << (window + 1)) - 1; \
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

		//PRECOMPUTE_WBITS_IMPL(blst_p1, POINTonE1, 384, fp, get_BLST_Rx().p)
		//POINTS_MULT_WBITS_IMPL(blst_p1, POINTonE1, 384, fp, get_BLST_Rx().p)



		//PRECOMPUTE_WBITS_IMPL(blst_p2, POINTonE2, 384x, fp2, get_BLST_Rx().p2)
		//POINTS_MULT_WBITS_IMPL(blst_p2, POINTonE2, 384x, fp2, get_BLST_Rx().p2)



















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

#define SCRATCH_SZ(ptype) (sizeof(ptype)==sizeof(POINTonE1) ? 8192 : 4096)

			 /* The intermediate infinity points are encoded as [0, 0, 1]. */

//#define PRECOMPUTE_WBITS_IMPL(prefix, ptype, bits, field, one) 
template<typename Field>
	void point_precompute_row_wbits(Point<Field> row[], size_t wbits,
		const Point_affine<Field>* point)
	{
		size_t i, j, n = (size_t)1 << (wbits - 1);
		bool_t inf = vec_is_zero(point, sizeof(*point));
		/* row[-1] is implicit infinity */
		vec_copy(&row[0], point, sizeof(*point));           /* row[0]=p*1     */
		vec_copy(&row[0].Z, &Field::one(), sizeof(row[0].Z));
		Point<Field>::doubling(&row[1], &row[0]);                  /* row[1]=p*(1+1) */
		vec_select(&row[1].Z, &Field::one(), &row[1].Z, sizeof(row[1].Z), inf);
		for (i = 2, j = 1; i < n; i += 2, j++)
			Point<Field>::add_affine(&row[i], &row[i - 1], point),  /* row[2]=p*(2+1) */
			Point<Field>::doubling(&row[i + 1], &row[j]),             /* row[3]=p*(2+2) */
			vec_select(&row[i + 1].Z, &Field::one(), &row[i + 1].Z, sizeof(row[i + 1].Z), inf);
	}                                                       /* row[4] ...     */

	template<typename Field>
	void points_to_affine_row_wbits(Point_affine<Field> dst[], Point<Field> src[],
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
			vec_copy(dst++, src++, sizeof(Point_affine<Field>));
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

	/* flat |points[n]| can be placed at the end of |table[n<<(wbits-1)]| */
	template<typename Field>
	void points_precompute_wbits(Point_affine<Field> table[], size_t wbits,
		const Point_affine<Field>* const points[],
		size_t npoints)
	{
		size_t total = npoints << (wbits - 1);
		size_t nwin = (size_t)1 << (wbits - 1);
		size_t nmin = wbits > 9 ? (size_t)1 : (size_t)1 << (9 - wbits);
		size_t i, top = 0;
		Point<Field>* rows, * row;
		const Point_affine<Field>* point = NULL;
		size_t stride = ((512 * 1024) / sizeof(Point_affine<Field>)) >> wbits;
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
				point_precompute_row_wbits<Field>(row, wbits, point);
			points_to_affine_row_wbits<Field>(&table[top], rows, wbits, stride);
			top += stride << (wbits - 1);
			npoints -= stride;
		}
		rows = row = (Point<Field>*)alloca(2 * sizeof(Point_affine<Field>) * npoints * nwin);
		for (i = 0; i < npoints; i++, row += nwin)
		{
			point = *points ? *points++ : point + 1;
			point_precompute_row_wbits<Field>(row, wbits, point);
		}
		points_to_affine_row_wbits<Field>(&table[top], rows, wbits, npoints);
	}

	template<typename Field>
	size_t points_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints)
	{
		return (sizeof(Point_affine<Field>) * npoints) << (wbits - 1);
	}

	template<typename Field>
	void points_mult_wbits_precompute(Point_affine<Field> table[], size_t wbits,
		const Point_affine<Field>* const points[],
		size_t npoints)
	{
		points_precompute_wbits<Field>(table, wbits, points, npoints);
	}





	//#define POINTS_MULT_WBITS_IMPL(prefix, ptype, bits, field, one) 
	template<typename Field>
	void point_gather_booth_wbits(
		Point<Field>* p,
		const Point_affine<Field> row[],
		size_t wbits, limb booth_idx)
	{
		bool_t booth_sign = (booth_idx >> wbits) & 1;
		bool_t idx_is_zero;
		const Point_affine<Field> infinity = { 0 };

		booth_idx &= ((limb)1 << wbits) - 1;
		idx_is_zero = is_zero(booth_idx);
		booth_idx -= 1 ^ idx_is_zero;
		vec_select(p, &infinity, &row[booth_idx], sizeof(row[0]), idx_is_zero);
		Point<Field>::cneg(p, booth_sign);
	}

	template<typename Field>
	void points_mult_wbits(Point<Field>* ret, const Point_affine<Field> table[],
		size_t wbits, size_t npoints,
		const byte* const scalars[], size_t nbits,
		Point<Field> scratch[])
	{
		limb wmask, wval;
		size_t i, j, z, nbytes, window, nwin = (size_t)1 << (wbits - 1);
		const byte* scalar, * const* scalar_s = scalars;
		const Point_affine<Field>* row = table;

		size_t scratch_sz = SCRATCH_SZ(Point<Field>);
		if (scratch == NULL) {
			scratch_sz /= 4; /* limit to 288K */
			scratch_sz = scratch_sz < npoints ? scratch_sz : npoints;
			scratch = (Point<Field>*)alloca(sizeof(Point<Field>) * scratch_sz);
		}

		nbytes = (nbits + 7) / 8; /* convert |nbits| to bytes */
		scalar = *scalar_s++;

		/* top excess bits modulo target window size */
		window = nbits % wbits; /* yes, it may be zero */
		wmask = ((limb)1 << (window + 1)) - 1;

		nbits -= window;
		z = is_zero(nbits);
		wval = (get_wval_limb(scalar, nbits - (z ^ 1), window + (z ^ 1)) << z) & wmask;
		wval = booth_encode(wval, wbits);
		point_gather_booth_wbits<Field>(&scratch[0], row, wbits, wval);
		row += nwin;

		i = 1; vec_zero(ret, sizeof(*ret));
		while (nbits > 0) {
			for (j = i; i < npoints; i++, j++, row += nwin) {
				if (j == scratch_sz)
					points_accumulate<Field>(ret, scratch, j), j = 0;
				scalar = *scalar_s ? *scalar_s++ : scalar + nbytes;
				wval = get_wval_limb(scalar, nbits - 1, window + 1) & wmask;
				wval = booth_encode(wval, wbits);
				point_gather_booth_wbits<Field>(&scratch[j], row, wbits, wval);
			}
			points_accumulate<Field>(ret, scratch, j);

			for (j = 0; j < wbits; j++)
				Point<Field>::doubling(ret, ret);

			window = wbits;
			wmask = ((limb)1 << (window + 1)) - 1;
			nbits -= window;
			i = 0; row = table; scalar_s = scalars;
		}

		for (j = i; i < npoints; i++, j++, row += nwin) {
			if (j == scratch_sz)
				points_accumulate<Field>(ret, scratch, j), j = 0;
			scalar = *scalar_s ? *scalar_s++ : scalar + nbytes;
			wval = (get_wval_limb(scalar, 0, window) << 1) & wmask;
			wval = booth_encode(wval, wbits);
			point_gather_booth_wbits<Field>(&scratch[j], row, wbits, wval);
		}
		pionts_accumulate<Field>(ret, scratch, j);
	}

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
			 */

	size_t pippenger_window_size(size_t npoints)
	{
		size_t wbits;

		for (wbits = 0; npoints >>= 1; wbits++);

		if (wbits > 12)
			return wbits - 3;
		else if (wbits > 8)
			return wbits - 2;
		else if (wbits > 4)
			return wbits - 1;

		return wbits ? 2 : 1;
	}

	//typedef struct { vec##bits X,Y,ZZZ,ZZ; } ptype##xyzz;
	//#define DECLARE_PRIVATE_POINTXYZZ(ptype, bits) 


#define POINTS_MULT_PIPPENGER_IMPL(prefix, ptype) \
void ptype##_integrate_buckets(ptype *out, ptype##xyzz buckets[], \
                                                  size_t wbits) \
{ \
    ptype##xyzz ret[1], acc[1]; \
    size_t n = (size_t)1 << wbits; \
\
    /* Calculate sum of x[i-1]*i for i=1 through 1<<|wbits|. */\
    vec_copy(acc, &buckets[--n], sizeof(acc)); \
    vec_copy(ret, &buckets[n], sizeof(ret)); \
    vec_zero(&buckets[n], sizeof(buckets[n])); \
    while (n--) { \
        ptype##xyzz_dadd(acc, acc, &buckets[n]); \
        ptype##xyzz_dadd(ret, ret, acc); \
        vec_zero(&buckets[n], sizeof(buckets[n])); \
    } \
    ptype##xyzz_to_Jacobian(out, ret); \
} \
\
void ptype##_bucket(ptype##xyzz buckets[], limb booth_idx, \
                           size_t wbits, const ptype##_affine *p) \
{ \
    bool_t booth_sign = (booth_idx >> wbits) & 1; \
\
    booth_idx &= (1<<wbits) - 1; \
    if (booth_idx--) \
        ptype##xyzz_dadd_affine(&buckets[booth_idx], &buckets[booth_idx], \
                                                     p, booth_sign); \
} \
\
void ptype##_prefetch(const ptype##xyzz buckets[], limb booth_idx, \
                             size_t wbits) \
{ \
    booth_idx &= (1<<wbits) - 1; \
    if (booth_idx--) \
        vec_prefetch(&buckets[booth_idx], sizeof(buckets[booth_idx])); \
} \
\
void ptype##s_tile_pippenger(ptype *ret, \
                                    const ptype##_affine *const points[], \
                                    size_t npoints, \
                                    const byte *const scalars[], size_t nbits, \
                                    ptype##xyzz buckets[], \
                                    size_t bit0, size_t wbits, size_t cbits) \
{ \
    limb wmask, wval, wnxt; \
    size_t i, z, nbytes; \
    const byte *scalar = *scalars++; \
    const ptype##_affine *point = *points++; \
\
    nbytes = (nbits + 7)/8; /* convert |nbits| to bytes */ \
    wmask = ((limb)1 << (wbits+1)) - 1; \
    z = is_zero(bit0); \
    bit0 -= z^1; wbits += z^1; \
    wval = (get_wval_limb(scalar, bit0, wbits) << z) & wmask; \
    wval = booth_encode(wval, cbits); \
    scalar = *scalars ? *scalars++ : scalar+nbytes; \
    wnxt = (get_wval_limb(scalar, bit0, wbits) << z) & wmask; \
    wnxt = booth_encode(wnxt, cbits); \
    npoints--;  /* account for prefetch */ \
\
    ptype##_bucket(buckets, wval, cbits, point); \
    for (i = 1; i < npoints; i++) { \
        wval = wnxt; \
        scalar = *scalars ? *scalars++ : scalar+nbytes; \
        wnxt = (get_wval_limb(scalar, bit0, wbits) << z) & wmask; \
        wnxt = booth_encode(wnxt, cbits); \
        ptype##_prefetch(buckets, wnxt, cbits); \
        point = *points ? *points++ : point+1; \
        ptype##_bucket(buckets, wval, cbits, point); \
    } \
    point = *points ? *points++ : point+1; \
    ptype##_bucket(buckets, wnxt, cbits, point); \
    ptype##_integrate_buckets(ret, buckets, cbits - 1); \
} \
\
void ptype##s_mult_pippenger(ptype *ret, \
                                    const ptype##_affine *const points[], \
                                    size_t npoints, \
                                    const byte *const scalars[], size_t nbits, \
                                    ptype##xyzz buckets[], size_t window) \
{ \
    size_t i, wbits, cbits, bit0 = nbits; \
    ptype tile[1]; \
\
    window = window ? window : pippenger_window_size(npoints); \
    vec_zero(buckets, sizeof(buckets[0]) << (window-1)); \
    vec_zero(ret, sizeof(*ret)); \
\
    /* top excess bits modulo target window size */ \
    wbits = nbits % window; /* yes, it may be zero */ \
    cbits = wbits + 1; \
    while (bit0 -= wbits) { \
        ptype##s_tile_pippenger(tile, points, npoints, scalars, nbits, \
                                      buckets, bit0, wbits, cbits); \
        ptype##_dadd(ret, ret, tile, NULL); \
        for (i = 0; i < window; i++) \
            ptype##_double(ret, ret); \
        cbits = wbits = window; \
    } \
    ptype##s_tile_pippenger(tile, points, npoints, scalars, nbits, \
                                  buckets, 0, wbits, cbits); \
    ptype##_dadd(ret, ret, tile, NULL); \
} \
\
size_t prefix##s_mult_pippenger_scratch_sizeof(size_t npoints) \
{   return sizeof(ptype##xyzz) << (pippenger_window_size(npoints)-1);   } \
void prefix##s_tile_pippenger(ptype *ret, \
                              const ptype##_affine *const points[], \
                              size_t npoints, \
                              const byte *const scalars[], size_t nbits, \
                              ptype##xyzz scratch[], \
                              size_t bit0, size_t window) \
{ \
    size_t wbits, cbits; \
\
    if (bit0 + window > nbits)  wbits = nbits - bit0, cbits = wbits + 1; \
    else                        wbits = cbits = window; \
    ptype##s_tile_pippenger(ret, points, npoints, scalars, nbits, scratch, \
                                 bit0, wbits, cbits); \
} \
void prefix##s_mult_pippenger(ptype *ret, \
                              const ptype##_affine *const points[], \
                              size_t npoints, \
                              const byte *const scalars[], size_t nbits, \
                              ptype##xyzz scratch[]) \
{ \
    if (npoints == 1) { \
        prefix##_from_affine(ret, points[0]); \
        ptype##_mult_w5(ret, ret, scalars[0], nbits); \
        return; \
    } \
    if ((npoints * sizeof(ptype##_affine) * 8 * 3) <= SCRATCH_LIMIT && \
        npoints < 32) { \
        ptype##_affine *table = (ptype##_affine *)alloca(npoints * sizeof(ptype##_affine) * 8); \
        ptype##s_precompute_wbits(table, 4, points, npoints); \
        ptype##s_mult_wbits(ret, table, 4, npoints, scalars, nbits, NULL); \
        return; \
    } \
    ptype##s_mult_pippenger(ret, points, npoints, scalars, nbits, scratch, 0); \
}

	//DECLARE_PRIVATE_POINTXYZZ(POINTonE1, 384)
	//POINTXYZZ_TO_JACOBIAN_IMPL(POINTonE1, 384, fp)
	//	POINTXYZZ_DADD_IMPL(POINTonE1, 384, fp)
	//	POINTXYZZ_DADD_AFFINE_IMPL(POINTonE1, 384, fp, get_BLST_Rx().p)
	//	POINTS_MULT_PIPPENGER_IMPL(blst_p1, POINTonE1)

	//	//DECLARE_PRIVATE_POINTXYZZ(POINTonE2, 384x)
	//	POINTXYZZ_TO_JACOBIAN_IMPL(POINTonE2, 384x, fp2)
	//	POINTXYZZ_DADD_IMPL(POINTonE2, 384x, fp2)
	//	POINTXYZZ_DADD_AFFINE_IMPL(POINTonE2, 384x, fp2, get_BLST_Rx().p2)
	//	POINTS_MULT_PIPPENGER_IMPL(blst_p2, POINTonE2)















		//#define POINTS_MULT_PIPPENGER_IMPL(prefix, Point<Field>) 
		template<typename Field>
	void point_integrate_buckets(Point<Field>* out, Pointxyzz<Field> buckets[],
		size_t wbits)
	{
		Pointxyzz<Field> ret[1], acc[1];
		size_t n = (size_t)1 << wbits;

		/* Calculate sum of x[i-1]*i for i=1 through 1<<|wbits|. */
		vec_copy(acc, &buckets[--n], sizeof(acc));
		vec_copy(ret, &buckets[n], sizeof(ret));
		vec_zero(&buckets[n], sizeof(buckets[n]));
		while (n--) {
			Pointxyzz<Field>::dadd(acc, acc, &buckets[n]);
			Pointxyzz<Field>::dadd(ret, ret, acc);
			vec_zero(&buckets[n], sizeof(buckets[n]));
		}
		Pointxyzz<Field>::to_Jacobian(out, ret);
	}

	template<typename Field>
	void point_bucket(Pointxyzz<Field> buckets[], limb booth_idx,
		size_t wbits, const Point_affine<Field>* p)
	{
		bool_t booth_sign = (booth_idx >> wbits) & 1;

		booth_idx &= (1 << wbits) - 1;
		if (booth_idx--)
			pointxyzz_dadd_affine<Field>(&buckets[booth_idx], &buckets[booth_idx],
				p, booth_sign);
	}

	template<typename Field>
	void point_prefetch(const Pointxyzz<Field> buckets[], limb booth_idx,
		size_t wbits)
	{
		booth_idx &= (1 << wbits) - 1;
		if (booth_idx--)
			vec_prefetch(&buckets[booth_idx], sizeof(buckets[booth_idx]));
	}

	template<typename Field>
	void points_tile_pippenger(Point<Field>* ret,
		const Point_affine<Field>* const points[],
		size_t npoints,
		const byte* const scalars[], size_t nbits,
		Pointxyzz<Field> buckets[],
		size_t bit0, size_t wbits, size_t cbits)
	{
		limb wmask, wval, wnxt;
		size_t i, z, nbytes;
		const byte* scalar = *scalars++;
		const Point_affine<Field>* point = *points++;

		nbytes = (nbits + 7) / 8; /* convert |nbits| to bytes */
		wmask = ((limb)1 << (wbits + 1)) - 1;
		z = is_zero(bit0);
		bit0 -= z ^ 1; wbits += z ^ 1;
		wval = (get_wval_limb(scalar, bit0, wbits) << z) & wmask;
		wval = booth_encode(wval, cbits);
		scalar = *scalars ? *scalars++ : scalar + nbytes;
		wnxt = (get_wval_limb(scalar, bit0, wbits) << z) & wmask;
		wnxt = booth_encode(wnxt, cbits);
		npoints--;  /* account for prefetch */

		for (i = 1; i < npoints; i++) {
			wval = wnxt;
			scalar = *scalars ? *scalars++ : scalar + nbytes;
			wnxt = (get_wval_limb(scalar, bit0, wbits) << z) & wmask;
			wnxt = booth_encode(wnxt, cbits);
			point_prefetch(buckets, wnxt, cbits);
			point = *points ? *points++ : point + 1;
			point_bucket(buckets, wval, cbits, point);
		}
		point = *points ? *points++ : point + 1;
		point_bucket(buckets, wnxt, cbits, point);
		point_integrate_buckets(ret, buckets, cbits - 1);
	}

	template<typename Field>
	void points_mult_pippenger(Point<Field>* ret,
		const Point_affine<Field>* const points[],
		size_t npoints,
		const byte* const scalars[], size_t nbits,
		Pointxyzz<Field> buckets[], size_t window)
	{
		size_t i, wbits, cbits, bit0 = nbits;
		Point<Field> tile[1];

		window = window ? window : pippenger_window_size(npoints);
		vec_zero(buckets, sizeof(buckets[0]) << (window - 1));
		vec_zero(ret, sizeof(*ret));

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
	size_t pointss_mult_pippenger_scratch_sizeof(size_t npoints)
	{
		return sizeof(Pointxyzz<Field>) << (pippenger_window_size(npoints) - 1);
	}
	template<typename Field>
	void points_tile_pippenger(Point<Field>* ret,
		const Point_affine<Field>* const points[],
		size_t npoints,
		const byte* const scalars[], size_t nbits,
		Pointxyzz<Field> scratch[],
		size_t bit0, size_t window)
	{
		size_t wbits, cbits;

		if (bit0 + window > nbits)  wbits = nbits - bit0, cbits = wbits + 1;
		else                        wbits = cbits = window;
		points_tile_pippenger<Field>(ret, points, npoints, scalars, nbits, scratch,
			bit0, wbits, cbits);
	}

	template<typename Field>
	void points_mult_pippenger(Point<Field>* ret,
		const Point_affine<Field>* const points[],
		size_t npoints,
		const byte* const scalars[], size_t nbits,
		Pointxyzz<Field> scratch[])
	{
		if (npoints == 1) {
			points_from_affine(ret, points[0]);
			point_mult_w5(ret, ret, scalars[0], nbits);
			return;
		}
		if ((npoints * sizeof(Point_affine<Field>) * 8 * 3) <= SCRATCH_LIMIT &&
			npoints < 32) {
			Point_affine<Field>* table = (Point_affine<Field>*)alloca(npoints * sizeof(Point_affine<Field>) * 8);
			points_precompute_wbits<Field>(table, 4, points, npoints);
			points_mult_wbits<Field>(ret, table, 4, npoints, scalars, nbits, NULL);
			return;
		}
		points_mult_pippenger<Field>(ret, points, npoints, scalars, nbits, scratch, 0);
	}


}
