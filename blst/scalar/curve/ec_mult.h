/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_SCALAR_CURVE_EC_MULT_H
#define BLST_SCALAR_CURVE_EC_MULT_H

#include <cstdlib>
#include <limits>
#include <type_traits>

namespace blst::scalar
{
	namespace detail
	{
		template<int WindowSize>
		struct fixed_window_traits
		{
			static_assert(WindowSize >= 2 && WindowSize <= 8,
				"fixed-window width must lie in [2, 8]");
			static constexpr size_t table_size =
				size_t{1} << (WindowSize - 1);
		};

		struct window_table_allocator
		{
			BLST_HD static void* allocate(size_t size)
			{
				return malloc(size);
			}

			BLST_HD static void deallocate(void* ptr)
			{
				free(ptr);
			}
		};
	}

	/*
	 * Extract window value of given bit size from bit offset in byte array.
	 * The bits are numbered from 0 (least significant) to n-1 (most significant).
	 * The window value is returned in the least significant bits of limb_t.
	 * The caller must mask off any loaded bits above the requested window and
	 * ensure that the source contains every addressed byte.
	 *
	 * Works up to 9 bits 
	 */
	BLST_HD BLST_FORCEINLINE limb_t get_wval(const byte* d, size_t off, size_t bits)
	{
		BLST_ASSERT(d != nullptr);
		BLST_ASSERT(bits >= 1 && bits <= 9);
		BLST_ASSERT(off <= std::numeric_limits<size_t>::max() - (bits - 1));
		size_t top = off + bits - 1;
		limb_t ret;

		ret = ((limb_t)d[top / 8] << 8) | d[off / 8];

		return ret >> (off % 8);
	}

	/* Works up to 25 bits and, like get_wval, leaves final masking to the caller. */
	BLST_HD BLST_INLINE limb_t get_wval_limb(const byte* d, size_t off, size_t bits)
	{
		BLST_ASSERT(d != nullptr);
		BLST_ASSERT(bits >= 1 && bits <= 25);
		BLST_ASSERT(off <= std::numeric_limits<size_t>::max() - (bits - 1));
		size_t i, top = (off + bits - 1) / 8;
		limb_t ret, mask = (limb_t)0 - 1;

		d += off / 8;
		top -= off / 8 - 1;

		/* this is not about constant-time-ness, but branch optimization */
		for (ret = 0, i = 0; i < 4;) {
			ret |= (*d & mask) << (8 * i);
			mask = (limb_t)0 - ((++i - top) >> (8 * sizeof(top) - 1));
			d += 1 & mask;
		}

		return ret >> (off % 8);
	}

	/* Booth encoding, given a window value wval of size w+1 bits
	 * where the LSB is the MSB of the next smaller window. return a  
	 * tuple (sign, idx) where sign is 0 or -1 and idx is in [0, 2^w].
	 * sign will be the w'th bit and up. E.g. for w=3, 
	 * 
	 *  111111 101 = sub index 5
	 *  000000 010 = add index 2
	 * 
	 * Window value encoding that utilizes the fact that -P is trivially
	 * calculated, which allows to halve the size of pre-computed table,
	 * is attributed to A. D. Booth, hence the name of the subroutines...
	 * 
	 * on input wval in [0, 2^(w+1)), returns in [-2^w, 2^w]
	 * 
	 * First we interpret wval as signed (w+1)=4 bit integer 
	 * 
	 * input  0 1 2 3 4 5 6 7  8  9 10 11 12 13 14 15
	 * interp 0 1 2 3 4 5 6 7 -8 -7 -6 -5 -4 -3 -2 -1
	 * sign   0 0 0 0 0 0 0 0  1  1  1  1  1  1  1  1
	 * 
	 * We then divide input by 2 and round up. Then apply sign
	 * 
	 * val    0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8
	 * signEx 0 1 1 2 2 3 3 4 -4 -5 -5 -6 -6 -7 -7 -8
	 * 
	 * the index that is actual used is then the w LSB of this
	 * 
	 * idx    0 1 1 2 2 3 3 4  4  3  3  2  2  1  1  0
	 * 
	 * We will use the sign bit to determine if we should add or subtract.
	 * 
	 * The divide by two is the LSB is the overlap bit for the next smaller
	 * window. If its 1, then the next window going to subtract its value. But this
	 * means that it want to subtract its value from the constant of this 
	 * current window. e.g. 7= 8-1. The next window is going to subtract 1.
	 * So this window we need to add the 8. In this window the `8` is the 
	 * LSB of wval.
	 * 
	 * Hopefully that makes sense. Took me while to understand...
	 * 
	 */
	BLST_HD BLST_FORCEINLINE limb_t booth_encode(limb_t wval, size_t w)
	{
		BLST_ASSERT(w > 0 && w < BLST_LIMB_SIZE - 1);
		BLST_ASSERT(wval < (limb_t{1} << (w + 1)));

		auto sign = wval >> w;
		limb_t mask = 0 - sign;  // "sign", 0 if wval < 2^w, else ~0
		BLST_launder(mask);             

		// val = divCeil(wval, 2)
		limb_t val = (wval + 1) >> 1;

		// ret = sign ? -val : val
		limb_t ret = (val ^ mask) - mask;	 

		/* &0x1f, but <=0x10, is index in table, rest is extended "sign" bit */
		return ret;
	}

	/*
		* Key feature of these constant-time subroutines is that they tolerate
		* zeros in most significant bit positions of the scalar[s], or in other
		* words, zero-padded scalar values. This means that one can and should
		* pass order's bit-length, which is customarily publicly known, instead
		* of the factual scalars' bit-lengths. This is facilitated by point
		* addition subroutines implemented to handle points at infinity, which
		* are encoded as Z==0. [Doubling algorithms handle such points at
		* infinity "naturally," since resulting Z is product of original Z.]
		*/
	template<typename ptype, int window_size>
	BLST_HD BLST_INLINE bool_t gather_booth_w(ptype& BLST_restrict p,
		const ptype table[detail::fixed_window_traits<window_size>::table_size],
		limb_t booth_encoding)
	{
		constexpr size_t table_size =
			detail::fixed_window_traits<window_size>::table_size;
		static_assert(std::is_trivially_default_constructible_v<ptype>);
		static_assert(std::is_trivially_destructible_v<ptype>);

		size_t i;
		bool_t booth_sign = (booth_encoding >> window_size) & 1;

		// this index takes into account the possible `borrow` from the 
		// next window and index the value that we need to add or subtract
		auto booth_idx = booth_encoding & ((limb_t{1} << window_size) - 1);

		// we will default p to the first element.
		// we then select the active row from the rest.
		vec_copy(&p, table, sizeof(ptype));
		/* ~6% with -Os, ~2% with -O3 ... */
		for (i = 2; i <= table_size; i++)
			ptype::ccopy(p, table[i - 1], byte_is_zero((byte)(i ^ booth_idx)));

		ptype::cneg(p, booth_sign);

		return byte_is_zero((byte)booth_idx);
	}


	// Computes the fixed-window table of multiples of `point`.
	// `row` contains fixed_window_traits<window_size>::table_size live objects.
	// Exact `row[0]`/`point` aliasing is supported; partial overlap is not.
	// 
	// row[0] = point * 1
	// row[1] = point * 2
	// ...
	// row[2^(window_size-1)-1] = point * (2^(window_size-1))
	// 
	template<typename ptype, int window_size>
	BLST_HD BLST_INLINE void precompute_w(ptype row[], const ptype& point)
	{
		constexpr size_t table_size =
			detail::fixed_window_traits<window_size>::table_size;
		size_t i, j;
		/* row[-1] is implicit infinity */
		vec_copy(&row[0], &point, sizeof(ptype));        /* row[0]=p*1     */
		ptype::doubling(row[1], point);                  /* row[1]=p*(1+1) */

		for (i = 2, j = 1; i < table_size; i += 2, j++)
		{
			ptype::add(row[i], row[j], row[j - 1]);      /* row[2]=p*(2+1) */
			ptype::doubling(row[i + 1], row[j]);         /* row[3]=p*(2+2) */
		}
	}                                                    /* row[4] ...     */

	// ret = scalar * point, (scalar is `bits` bits long)
	//
	// The operation is constant-time with respect to scalar contents for fixed
	// public `bits` and `window_size`. `bits` is nonzero and `scalar` addresses
	// ceil(bits/8) readable bytes. `ret` may exactly alias `point`; the scalar
	// storage is disjoint from both objects.
	//
	// we will compute this using windowing. The basic idea is 
	// to preprocess multiples of `point`:
	// 
	//   table = { point, 2*point, 3*point, ..., (2^(w-1))*point} 
	// 
	// We will then decompose the scaler into "windows" of w bits each,
	//
	//   s0, s1, s2, ..., sk
	//
	// For each i = k...0, we will select the precomputed multiple
	// corresponding to si, double the current accumulator w times,
	// and add the selected multiple to the accumulator. i.e.
	//
	// for i = k..0: 
	//    acc = acc * (2^w)
	//    acc = acc + (si ? table[si-1] : 0)
	//
	// To improve efficiency, we use Booth encoding. The basic idea
	// is that we can efficiently add and subtract points. However, 
	// the algorithm above only adds. By using Booth encoding, we can
	// essentially halve the table size for free. If a value is in 
	// the upper half of the window, then we will scale the accumulator
	// to the top of the window and then subtract the corresponding 
	// point instead of adding. This allows us to halve the size of the 
	// table.
	// 
	// For example, 7 = 0111 in binary. If we have a window size of 1. 
	// we could represent 7 = 1 + 2 + 4. However, if we use both we can
	// represent 7 = 8 - 1. Halving the table size and reducing the number of
	// operations from 2 to 1.
	// 
	// let 
	// 
	//  scaler = sn || ... || si || ... || s1 || s0 
	// 
	// be the bit decompositon. For a window size of w, 
	// we will ``recode'' each window as 
	// 
	//   d0 = s{ w-1}    |...|s1     |s0   |s{-1}
	//   d1 = s{2w-1}    |...|s{w+1} |sw   |s{w-1}
	// ...					 
	//   dj = s{(j+1)w-1}|...|s{jw+1}|s{jw}|s{jw-1}
	//            							  ^
	//            							  |
	//         overlap bit for correcting ----+
	//         the prior window's subtraction
	//						 
	// Note that each dj is w+1 bits long and the LSB is overlapped
	// with the MSB of the previous window. This overlap allows for a 
	// cleaver method for deciding when we need to scale the accumulator
	// to the top of the window. Note that the first window d0 is special
	// and we will add a zero bit as the LSB.
	// 
	// We then interpret them as signed w+1 bit integers in the 
	// range [-2^w, 2^w - 1]. If s does not have enough bits, we
	// make the MSB windows smaller.
	//
	// We then have this cleaver booth encoding step. Where we use the 
	// overlapping LSB of each window to determin if the next smaller window 
	// wants us to add in an extra term from which they will subtract from.
	// The amount we will add or subtract is 
	// 
	//    di' = (di+1) / 2
	// 
	// If the next smaller window is going to subtract, we will add an extra
	// "+1" to our current window. e.g. 7 = 8 - 1. The next window is going to
	// subtract 1 while we will add 1 which is 8 after scaling.
	// 
	// for i = j..0:
	//    acc = acc * (2^w)
	//    acc = acc + (dj ?  sign(dj') * table[ abs(dj) - 1 ] : 0)
	//
	// 
	// For example, let w=3. We create w+1 = 4 bit overlapping windows.
	// 
	// scalar: 0111011111 | wval    d, sign, index
	// window  01 		  |   1     1   0    1      add 1
	//          1110 	  |  14    -2   1    1      sub 1   our prior window will have added 8
	//             0111   |   7     7   0    4      add 4
	//                1110   14    -2   1    1	    sub 1   our prior window will have added 8      
	//                    |        
	//
	// a bit magical but it works. Tha advantage is that the table is smaller.
	//
	template<typename ptype, int window_size>
	BLST_HD BLST_INLINE void mult_window(ptype& ret, const ptype& point,
		const byte* scalar, const size_t bits)
	{
		constexpr size_t table_size =
			detail::fixed_window_traits<window_size>::table_size;
		BLST_ASSERT(bits != 0);
		BLST_ASSERT(scalar != nullptr);
		limb_t wmask, wval;
		size_t j, window;
		ptype sum, row;
		bool_t sum_is_inf, row_is_inf, ret_is_inf;

		// pre-compute table of multiples
		// point, 2*point, 3*point, ..., (2^(w-1))*point
		ptype table[table_size];
		precompute_w<ptype, window_size>(table, point);

		/* top excess bits modulo target window size */
		window = bits % window_size;  /* yes, it may be zero */
		wmask = ((limb_t)1 << (window + 1)) - 1;

		auto bitIndex = bits - window;
		wval = bitIndex ? get_wval(scalar, bitIndex - 1, window + 1)
			: (limb_t)scalar[0] << 1;

		wval &= wmask;
		auto booth_encoding = booth_encode(wval, window_size);
		ret_is_inf = gather_booth_w<ptype, window_size>(ret, table, booth_encoding);

		while (bitIndex > 0) {
			for (j = 0; j < window_size; j++)
				ptype::doubling(ret, ret);

			window = window_size;
			wmask = ((limb_t)1 << (window + 1)) - 1;
			bitIndex -= window;

			wval = bitIndex ? get_wval(scalar, bitIndex - 1, window + 1)
				: (limb_t)scalar[0] << 1;

			wval &= wmask;
			booth_encoding = booth_encode(wval, window_size);
			row_is_inf = gather_booth_w<ptype, window_size>(row, table, booth_encoding);
			ptype::dadd(sum, ret, row);
			ptype::ccopy(ret, sum, (ret_is_inf | row_is_inf) ^ 1);
			sum_is_inf = vec_is_zero(&ret.z, sizeof(ret.z));
			ret_is_inf |= sum_is_inf;
			row_is_inf |= sum_is_inf;
			ptype::ccopy(ret, row, ret_is_inf);
			ret_is_inf &= row_is_inf;
		}

		vec_czero(&ret.z, sizeof(ret.z), ret_is_inf);
	}

	// same as mult_w but its an MSM for a small number of points
	// 
	//   ret = sum_{i =0..npoints-1} scalars[i] * points[i]
	//
	// Either points or a pre-computed table must be provided. If points is
	// provided, table is populated for the caller or allocated on the heap when
	// it is NULL. If points is NULL, table must already be populated.
	//
	// The first point and scalar pointers must be non-NULL. Later NULL entries
	// select the next contiguous object, preserving the compact BLST convention.
	// On invalid input, size overflow, or allocation failure, ret is unchanged
	// and the function returns zero. `ret` is disjoint from the table and scalar
	// storage. Point/scalar descriptor arrays either contain one explicit head
	// per item or a null entry after a head to select contiguous continuation.
	// For fixed public sizes, table selection and group accumulation are
	// constant-time with respect to scalar contents.
	// 
	// this is better than mult_w because the doubling of the 
	// accumulator is shared between all inputs. This can be a big 
	// savings.
	//
	template<typename ptype, int window_size,
		typename allocator_type = detail::window_table_allocator>
	[[nodiscard]] BLST_HD BLST_FORCEINLINE bool_t mult_window_n(ptype& ret,
		const ptype* points[], size_t npoints,
		const byte* scalars[], size_t bits,
		ptype table[][detail::fixed_window_traits<window_size>::table_size])
	{
		constexpr size_t table_size =
			detail::fixed_window_traits<window_size>::table_size;
		static_assert(std::is_trivially_default_constructible_v<ptype>);
		static_assert(std::is_trivially_destructible_v<ptype>);

		limb_t wmask, wval;
		size_t i, j, window, nbytes;
		const byte* scalar, ** scalar_s = scalars;
		ptype sum, row;
		bool_t sum_is_inf, row_is_inf, ret_is_inf;
		bool_t call_free = 0;

		if (npoints == 0 || bits == 0 || bits > (size_t)-1 - 7 ||
			scalars == NULL || scalars[0] == NULL ||
			(points == NULL && table == NULL) ||
			(points != NULL && points[0] == NULL))
			return 0;

		if (table == NULL)
		{
			const size_t table_row_size = sizeof(*table);
			if (npoints > (size_t)-1 / table_row_size)
				return 0;

			table = static_cast<ptype(*)[table_size]>(
				allocator_type::allocate(table_row_size * npoints));

			if (table == NULL)
				return 0;
			call_free = 1;
		}

		if (points != NULL) {
			const ptype* point = NULL;
			for (i = 0; i < npoints; i++)
				point = *points ? *points++ : point + 1,
				precompute_w<ptype, window_size>(table[i], *point);
		}

		nbytes = (bits + 7) / 8; /* convert |bits| to bytes */
		scalar = *scalar_s++;

		/* top excess bits modulo target window size */
		window = bits % window_size; /* yes, it may be zero */
		wmask = ((limb_t)1 << (window + 1)) - 1;

		bits -= window;
		if (bits > 0)
			wval = get_wval(scalar, bits - 1, window + 1) & wmask;
		else
			wval = (scalar[0] << 1) & wmask;

		wval = booth_encode(wval, window_size);
		ret_is_inf = gather_booth_w<ptype, window_size>(ret, table[0], wval);

		// main loop, we will process `window_size` bits at a time and
		// for each window, we will loop over all points
		i = 1;
		while (bits > 0) {
			for (; i < npoints; i++) {
				scalar = *scalar_s ? *scalar_s++ : scalar + nbytes;
				wval = get_wval(scalar, bits - 1, window + 1) & wmask;
				wval = booth_encode(wval, window_size);
				row_is_inf = gather_booth_w<ptype, window_size>(row, table[i], wval);
				ptype::dadd(sum, ret, row);
				ptype::ccopy(ret, sum, (ret_is_inf | row_is_inf) ^ 1);
				sum_is_inf = vec_is_zero(&ret.z, sizeof(ret.z));
				ret_is_inf |= sum_is_inf;
				row_is_inf |= sum_is_inf;
				ptype::ccopy(ret, row, ret_is_inf);
				ret_is_inf &= row_is_inf;
			}

			// scale the accumulator to the top of the window
			for (j = 0; j < window_size; j++)
				ptype::doubling(ret, ret);

			window = window_size;
			wmask = ((limb_t)1 << (window + 1)) - 1;
			bits -= window;
			i = 0; scalar_s = scalars;
		}

		// final window, if any
		for (; i < npoints; i++) {
			scalar = *scalar_s ? *scalar_s++ : scalar + nbytes;
			wval = (scalar[0] << 1) & wmask;
			wval = booth_encode(wval, window_size);
			row_is_inf = gather_booth_w<ptype, window_size>(row, table[i], wval);
			ptype::dadd(sum, ret, row);
			ptype::ccopy(ret, sum, (ret_is_inf | row_is_inf) ^ 1);
			sum_is_inf = vec_is_zero(&ret.z, sizeof(ret.z));
			ret_is_inf |= sum_is_inf;
			row_is_inf |= sum_is_inf;
			ptype::ccopy(ret, row, ret_is_inf);
			ret_is_inf &= row_is_inf;
		}

		vec_czero(&ret.z, sizeof(ret.z), ret_is_inf);
		if (call_free)
			allocator_type::deallocate(table);

		return 1;
	}

	// A fixed-base MSM with affine lookup rows. The caller prepares all rows
	// from public bases; scalar digits only select through a full linear scan.
	// Mixed dadd handles zero digits, equality, opposites, and infinity without
	// secret-dependent control flow. Scalars are contiguous little-endian bytes.
	template<typename ptype, typename affine_type, int window_size>
	[[nodiscard]] BLST_HD BLST_FORCEINLINE bool_t mult_window_n_affine(
		ptype& ret, const byte* scalars, size_t npoints, size_t bits,
		const affine_type table[][detail::fixed_window_traits<window_size>::table_size])
	{
		constexpr size_t table_size =
			detail::fixed_window_traits<window_size>::table_size;
		static_assert(window_size >= 2 && window_size <= 7);
		static_assert(std::is_trivially_default_constructible_v<affine_type>);
		static_assert(std::is_trivially_destructible_v<affine_type>);
		if (npoints == 0 || bits == 0 || bits > (size_t)-1 - 7 ||
			scalars == NULL || table == NULL)
			return 0;
		const size_t nbytes = (bits + 7) / 8;
		if (npoints > (size_t)-1 / nbytes)
			return 0;

		ret = ptype{};
		size_t window = bits % window_size;
		limb_t wmask = ((limb_t)1 << (window + 1)) - 1;
		bits -= window;
		for (;;) {
			for (size_t i = 0; i < npoints; ++i) {
				const byte* scalar = scalars + i * nbytes;
				const limb_t extracted = bits > 0
					? get_wval(scalar, bits - 1, window + 1) & wmask
					: ((limb_t)scalar[0] << 1) & wmask;
				const limb_t encoded = booth_encode(extracted, window_size);
				const limb_t index = encoded & (((limb_t)1 << window_size) - 1);
				affine_type selected;
				vec_copy(&selected, table[i], sizeof(selected));
				for (size_t j = 2; j <= table_size; ++j)
					vec_select(&selected, &table[i][j - 1], &selected,
						sizeof(selected), byte_is_zero((byte)(j ^ index)));
				affine_type::cneg(selected, (encoded >> window_size) & 1);
				vec_czero(&selected, sizeof(selected),
					byte_is_zero((byte)index));
				ptype sum;
				ptype::dadd(sum, ret, selected);
				ret = sum;
			}
			if (bits == 0) break;
			for (size_t j = 0; j < window_size; ++j)
				ptype::doubling(ret, ret);
			window = window_size;
			wmask = ((limb_t)1 << (window + 1)) - 1;
			bits -= window;
		}
		return 1;
	}

	/* Variable-time wNAF recoding for public little-endian scalar data. Bits in
	 * the final source byte above ScalarBits are ignored. The output contains
	 * ScalarBits+1 entries to retain a possible final positive carry. */
	template<size_t ScalarBits, int Width>
	BLST_HD BLST_INLINE size_t recode_wnaf_vartime(
		int8_t digits[ScalarBits + 1], const byte* scalar)
	{
		static_assert(ScalarBits > 0);
		static_assert(Width >= 2 && Width <= 7);
		constexpr size_t scalar_limbs = (ScalarBits + BLST_LIMB_SIZE - 1) /
			BLST_LIMB_SIZE;
		// A negative top digit can carry one bit beyond ScalarBits. Retaining
		// that carry is required for the final positive wNAF digit.
		constexpr size_t limbs_count = scalar_limbs + 1;
		limb_t value[limbs_count] = {};
		limbs_from_le_bytes(value, scalar, (ScalarBits + 7) / 8);
		if constexpr ((ScalarBits % 8) != 0) {
			constexpr size_t top_bits = ScalarBits % BLST_LIMB_SIZE;
			value[scalar_limbs - 1] &= (limb_t{1} << top_bits) - 1;
		}
		for (size_t i = 0; i <= ScalarBits; ++i)
			digits[i] = 0;

		size_t used = 0;
		while (used <= ScalarBits) {
			limb_t nonzero = 0;
			for (size_t i = 0; i < limbs_count; ++i)
				nonzero |= value[i];
			if (nonzero == 0)
				break;

			if (value[0] & 1) {
				const unsigned raw = static_cast<unsigned>(
					value[0] & (((limb_t)1 << Width) - 1));
				const int digit = raw >= (1U << (Width - 1))
					? static_cast<int>(raw) - (1 << Width)
					: static_cast<int>(raw);
				digits[used] = static_cast<int8_t>(digit);
				if (digit > 0) {
					value[0] -= static_cast<limb_t>(digit);
				}
				else {
					limb_t add = static_cast<limb_t>(-digit);
					for (size_t i = 0; i < limbs_count; ++i) {
						wide_limb_t sum;
						cc_add_wide_limb_limb(&sum, value[i], add);
						value[i] = (limb_t)sum.lo;
						add = (limb_t)sum.hi;
					}
				}
			}

			for (size_t i = 0; i + 1 < limbs_count; ++i)
				value[i] = (value[i] >> 1) |
					(value[i + 1] << (BLST_LIMB_SIZE - 1));
			value[limbs_count - 1] >>= 1;
			++used;
		}

		return used;
	}

}

#endif
