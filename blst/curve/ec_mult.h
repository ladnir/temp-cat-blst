/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_EC_MULT_H__
#define __BLS12_381_ASM_EC_MULT_H__

namespace blst
{

	/*
	 * Extract window value of given bit size from bit offset in byte array.
	 * The bits are numbered from 0 (least significant) to n-1 (most significant).
	 * The window value is returned in the least significant bits of limb_t.
	 * The caller must ensure that 'bits' is not greater than what the
	 * subroutine can handle, see below.
	 *
	 * Works up to 9 bits 
	 */
	BLST_HD BLST_FORCEINLINE limb_t get_wval(const byte* d, size_t off, size_t bits)
	{
		size_t top = off + bits - 1;
		limb_t ret;

		ret = ((limb_t)d[top / 8] << 8) | d[off / 8];

		return ret >> (off % 8);
	}

	/* Works up to 25 bits. */
	BLST_HD BLST_INLINE limb_t get_wval_limb(const byte* d, size_t off, size_t bits)
	{
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
		BLST_ASSERT(wval < (1 << (w + 1)));

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
		const ptype table[1 << (window_size - 1)],
		limb_t booth_encoding)
	{

		size_t i;
		bool_t booth_sign = (booth_encoding >> window_size) & 1;

		// this index takes into account the possible `borrow` from the 
		// next window and index the value that we need to add or subtract
		auto booth_idx = booth_encoding & (1 << window_size) - 1;

		// we will default p to the first element.
		// we then select the active row from the rest.
		vec_copy(&p, table, sizeof(ptype));
		/* ~6% with -Os, ~2% with -O3 ... */
		for (i = 2; i <= 1 << (window_size - 1); i++)
			ptype::ccopy(p, table[i - 1], byte_is_zero((byte)(i ^ booth_idx)));

		ptype::cneg(p, booth_sign);

		return byte_is_zero((byte)booth_idx);
	}


	// computes pre-computation table of multiplies of `point`
	// 
	// row[0] = point * 1
	// row[1] = point * 2
	// ...
	// row[2^(window_size-1)-1] = point * (2^(window_size-1))
	// 
	template<typename ptype, int window_size>
	BLST_HD BLST_INLINE void precompute_w(ptype row[], const ptype& point)
	{
		size_t i, j;
		/* row[-1] is implicit infinity */
		vec_copy(&row[0], &point, sizeof(ptype));        /* row[0]=p*1     */
		ptype::doubling(row[1], point);                  /* row[1]=p*(1+1) */

		for (i = 2, j = 1; i < 1 << (window_size - 1); i += 2, j++)
		{
			ptype::add(row[i], row[j], row[j - 1]);      /* row[2]=p*(2+1) */
			ptype::doubling(row[i + 1], row[j]);         /* row[3]=p*(2+2) */
		}
	}                                                    /* row[4] ...     */

	// ret = scalar * point, (scalar is `bits` bits long)
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
		limb_t wmask, wval;
		size_t j, window;
		ptype sum, row;
		bool_t sum_is_inf, row_is_inf, ret_is_inf;

		// pre-compute table of multiples
		// point, 2*point, 3*point, ..., (2^(w-1))*point
		ptype table[1 << (window_size - 1)];
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
			ptype::dadd(sum, ret, row, NULL);
			ptype::ccopy(ret, sum, (ret_is_inf | row_is_inf) ^ 1);
			sum_is_inf = vec_is_zero(&ret.Z, sizeof(ret.Z));
			ret_is_inf |= sum_is_inf;
			row_is_inf |= sum_is_inf;
			ptype::ccopy(ret, row, ret_is_inf);
			ret_is_inf &= row_is_inf;
		}

		vec_czero(&ret.Z, sizeof(ret.Z), ret_is_inf);
	}

	// same as mult_w but its an MSM for a small number of points
	// 
	//   ret = sum_{i =0..npoints-1} scalars[i] * points[i]
	//
	// Either points or pre-computed table must be provided.
	// If points is provided, the table will be computed on the stack.
	// If table is provided, points is ignored.
	// 
	// this is better than mult_w because the doubling of the 
	// accumulator is shared between all inputs. This can be a big 
	// savings.
	//
	template<typename ptype, int window_size>
	BLST_HD BLST_INLINE void mult_window_n(ptype& ret,
		const ptype* points[], size_t npoints,
		const byte* scalars[], size_t bits,
		ptype table[][1 << (window_size - 1)])
	{
		limb_t wmask, wval;
		size_t i, j, window, nbytes;
		const byte* scalar, ** scalar_s = scalars;
		ptype sum, row;
		bool_t sum_is_inf, row_is_inf, ret_is_inf;

		if (table == NULL)
			table = (ptype(*)[1 << (window_size - 1)])alloca((1 << (window_size - 1)) * sizeof(ptype) *
				npoints);

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
				ptype::dadd(sum, ret, row, NULL);
				ptype::ccopy(ret, sum, (ret_is_inf | row_is_inf) ^ 1);
				sum_is_inf = vec_is_zero(&ret.Z, sizeof(ret.Z));
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
			ptype::dadd(sum, ret, row, NULL);
			ptype::ccopy(ret, sum, (ret_is_inf | row_is_inf) ^ 1);
			sum_is_inf = vec_is_zero(&ret.Z, sizeof(ret.Z));
			ret_is_inf |= sum_is_inf;
			row_is_inf |= sum_is_inf;
			ptype::ccopy(ret, row, ret_is_inf);
			ret_is_inf &= row_is_inf;
		}

		vec_czero(&ret.Z, sizeof(ret.Z), ret_is_inf);
	}

}

#endif
