/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#include "backend_portable.h"



namespace blst
{
	/**
	 * @brief Performs Montgomery multiplication of two large integers.
	 *
	 * Computes `(a * b * R^(-1)) mod p` using the Montgomery multiplication algorithm,
	 * where `R = 2^(BLST_LIMB_SIZE * n)`. The result is stored in `ret`. Assumes `n > 0`
	 * and `n` is even.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_mul_mont_n(limb_t* ret, const limb_t* a, const limb_t* b,
		const limb_t* p, limb_t n0)
	{
		static_assert(n > 0 && n % 2 == 0, "n must be positive and even");

		limb_t mask, borrow, mx, hi, carry;
		limb_t  tmp[n + 1];
		size_t i, j;

		for (mx = b[0], hi = 0, i = 0; i < n; i++) {
			wide_limb_t  wide;

			cc_mul_add_wide_limb_limb(&wide, mx, a[i], hi);
			tmp[i] = (limb_t)wide.lo;
			hi = (limb_t)(wide.hi);
		}
		mx = n0 * tmp[0];
		tmp[i] = hi;

		for (carry = 0, j = 0; ; ) {
			wide_limb_t  wide;
			// wide = mx * p + tmp
			cc_mul_add_wide_limb_limb(&wide, mx, p[0], tmp[0]);
			hi = wide.hi;

			for (i = 1; i < n; i++) {
				// Accumulate the product limb, previous high limb, and current temporary limb.
				cc_mul_add_add_wide_limb_limb(&wide, mx, p[i], hi, tmp[i]);
				tmp[i - 1] = wide.lo;
				hi = wide.hi;
			}
			// wide = hi + carry + tmp
			cc_add_carry_wide_limb_limb(&wide, hi, tmp[i], carry);
			tmp[i - 1] = wide.lo;
			carry = wide.hi;


			if (++j == n)
				break;

			for (mx = b[j], hi = 0, i = 0; i < n; i++) {
				// wide = mx * a + hi + tmp
				cc_mul_add_add_wide_limb_limb(&wide, mx, a[i], hi, tmp[i]);
				tmp[i] = wide.lo;
				hi = wide.hi;
			}
			mx = n0 * tmp[0];
			cc_add_wide_limb_limb(&wide, hi, carry);
			tmp[i] = wide.lo;
			carry = wide.hi;
		}

		for (borrow = 0, i = 0; i < n; i++) {
			// wide = tmp - (p + borrow)
			wide_limb_t wide;
			cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
			ret[i] = wide.lo;
			borrow = wide.hi & 1;
		}

		mask = carry - borrow;
		BLST_launder(mask);

		for (i = 0; i < n; i++)
			ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);

	}

	template
	BLST_HD BLST_CPP_INLINE void cc_mul_mont_n<NLIMBS(256)>(limb_t* ret, const limb_t* a, const limb_t* b,
		const limb_t* p, limb_t n0);
	template
	BLST_HD BLST_CPP_INLINE void cc_mul_mont_n<NLIMBS(384)>(limb_t* ret, const limb_t* a, const limb_t* b,
		const limb_t* p, limb_t n0);


	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_add_mod_n(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t mask, carry, borrow;
		limb_t tmp[n];
		size_t i;

		for (carry = 0, i = 0; i < n; i++) {
			wide_limb_t wide;
			cc_add_carry_wide_limb_limb(&wide, a[i], b[i], carry);
			tmp[i] = wide.lo;
			carry = wide.hi;
		}

		for (borrow = 0, i = 0; i < n; i++) {
			wide_limb_t wide;
			cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
			ret[i] = wide.lo;
			borrow = wide.hi & 1;
		}

		mask = carry - borrow;
		BLST_launder(mask);

		for (i = 0; i < n; i++)
			ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);

	}

	template
	BLST_HD BLST_CPP_INLINE void cc_add_mod_n<NLIMBS(256)>(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[]);
	template
	BLST_HD BLST_CPP_INLINE void cc_add_mod_n<NLIMBS(384)>(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[]);



	/**
	 * @brief Subtracts two large integers modulo a given modulus.
	 *
	 * Requires `n > 0` and reduced operands `0 <= a,b < p`.
	 * Exact output aliasing with either input is supported; partial overlap is not.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_sub_mod_n(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t mask, carry, borrow;
		size_t i;

		for (borrow = 0, i = 0; i < n; i++) {
			wide_limb_t wide;
			cc_sub_borrow_wide_limb_limb(&wide, a[i], b[i], borrow);
			ret[i] = wide.lo;
			borrow = wide.hi & 1;
		}

		mask = 0 - borrow;
		BLST_launder(mask);

		for (carry = 0, i = 0; i < n; i++) {
			wide_limb_t wide;
			cc_add_carry_wide_limb_limb(&wide, ret[i], p[i] & mask, carry);
			ret[i] = wide.lo;
			carry = wide.hi;
		}
	}
	template
	BLST_HD BLST_CPP_INLINE void cc_sub_mod_n<NLIMBS(256)>(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[]);
	template
	BLST_HD BLST_CPP_INLINE void cc_sub_mod_n<NLIMBS(384)>(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[]);


	/**
	 * @brief Multiplies a large integer by 3 modulo a given modulus.
	 *
	 * Requires `n > 0` and `0 <= a < p`.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_mul_by_3_mod_n(limb_t ret[], const limb_t a[], const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t mask, carry, borrow;
		limb_t tmp[n];
		limb_t two_a[n];

		size_t i;

		for (carry = 0, i = 0; i < n; i++) {
			limb_t a_i = a[i];
			tmp[i] = a_i << 1 | carry;
			carry = a_i >> (BLST_LIMB_SIZE - 1);
		}

		for (borrow = 0, i = 0; i < n; i++) {
			wide_limb_t wide;
			cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
			two_a[i] = wide.lo;
			borrow = wide.hi & 1;
		}

		mask = carry - borrow;
		BLST_launder(mask);

		for (i = 0; i < n; i++)
			two_a[i] = (two_a[i] & ~mask) | (tmp[i] & mask);

		for (carry = 0, i = 0; i < n; i++) {
			wide_limb_t wide;
			cc_add_carry_wide_limb_limb(&wide, a[i], two_a[i], carry);
			tmp[i] = wide.lo;
			carry = wide.hi;
		}

		for (borrow = 0, i = 0; i < n; i++) {
			wide_limb_t wide;
			cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
			ret[i] = wide.lo;
			borrow = wide.hi & 1;
		}

		mask = carry - borrow;
		BLST_launder(mask);

		for (i = 0; i < n; i++)
			ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);
	}

	template
	BLST_HD BLST_CPP_INLINE void cc_mul_by_3_mod_n<NLIMBS(256)>(limb_t ret[], const limb_t a[], const limb_t p[]);
	template
	BLST_HD BLST_CPP_INLINE void cc_mul_by_3_mod_n<NLIMBS(384)>(limb_t ret[], const limb_t a[], const limb_t p[]);


	/**
	 * @brief Performs a left shift modulo operation on a large integer.
	 *
	 * Requires `count > 0`, `n > 0`, and `0 <= a < p`.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_lshift_mod_n(limb_t ret[], const limb_t a[], size_t count,
		const limb_t p[])
	{
		BLST_ASSERT(count != 0);
		BLST_ASSUME(count != 0);
		BLST_ASSUME(n != 0);
		limb_t mask, carry, borrow;
		limb_t tmp[n];
		size_t i;

		while (count--) {
			for (carry = 0, i = 0; i < n; i++) {
				limb_t a_i = a[i];
				tmp[i] = a_i << 1 | carry;
				carry = a_i >> (BLST_LIMB_SIZE - 1);
			}

			for (borrow = 0, i = 0; i < n; i++) {
				wide_limb_t wide;
				cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
				ret[i] = wide.lo;
				borrow = wide.hi & 1;
			}

			mask = carry - borrow;
			BLST_launder(mask);

			for (i = 0; i < n; i++)
				ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);

			a = ret;
		}
	}

	template
	BLST_HD BLST_CPP_INLINE void cc_lshift_mod_n<NLIMBS(256)>(limb_t ret[], const limb_t a[], size_t count,
		const limb_t p[]);
	template
	BLST_HD BLST_CPP_INLINE void cc_lshift_mod_n<NLIMBS(384)>(limb_t ret[], const limb_t a[], size_t count,
		const limb_t p[]);


	/**
	 * @brief Conditionally negates a large integer modulo a given modulus.
	 *
	 * If `flag` is nonzero, computes `(-a) mod p`; otherwise, returns `a mod p`.
	 * The result is stored in `ret`. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_cneg_mod_n(limb_t ret[], const limb_t a[], bool_t flag,
		const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t borrow, mask;
		limb_t tmp[n];
		size_t i;

		for (borrow = 0, i = 0; i < n; i++) {
			wide_limb_t wide;
			cc_sub_borrow_wide_limb_limb(&wide, p[i], a[i], borrow);
			tmp[i] = wide.lo;
			borrow = wide.hi & 1;
		}

		/* Match the assembly contract: zero is false and every nonzero value is true. */
		flag = (flag | ((limb_t)0 - flag)) >> (BLST_LIMB_SIZE - 1);
		flag &= vec_is_zero(a, n * sizeof(limb_t)) ^ 1;
		mask = (limb_t)0 - flag;

		for (i = 0; i < n; i++)
			ret[i] = (a[i] & ~mask) | (tmp[i] & mask);
	}

	template
	BLST_HD BLST_CPP_INLINE void cc_cneg_mod_n<NLIMBS(256)>(limb_t ret[], const limb_t a[], bool_t flag,
		const limb_t p[]);
	template
	BLST_HD BLST_CPP_INLINE void cc_cneg_mod_n<NLIMBS(384)>(limb_t ret[], const limb_t a[], bool_t flag,
		const limb_t p[]);


	/**
	 * @brief Checks if a large integer `a` is reduced modulo a given modulus
	 * and not zero.
	 *
	 * Returns one exactly when `0 < a < p`. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE limb_t cc_check_mod_n(const byte a[], const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t borrow, ai, acc;
		size_t i, j;

		for (acc = borrow = 0, i = 0; i < n; i++) {
			for (ai = 0, j = 0; j < 8 * sizeof(limb_t); j += 8)
				ai |= (limb_t)(*a++) << j;
			acc |= ai;
			wide_limb_t wide;
			cc_sub_borrow_wide_limb_limb(&wide, ai, p[i], borrow);
			borrow = wide.hi & 1;
		}

		return borrow & (is_zero(acc) ^ 1);
	}

	template
	BLST_HD BLST_CPP_INLINE limb_t cc_check_mod_n<NLIMBS(256)>(const byte a[], const limb_t p[]);
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_check_mod_n<NLIMBS(384)>(const byte a[], const limb_t p[]);

	/**
	 * @brief Checks if a large integer `a` is reduced modulo a given modulus
	 *
	 * Returns 0 if `a` is greater than or equal to `p`, otherwise returns 1.
	 * Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE limb_t cc_is_redc_mod_n(const byte a[], const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t borrow, ai;
		size_t i, j;

		for (borrow = 0, i = 0; i < n; i++) {
			for (ai = 0, j = 0; j < 8 * sizeof(limb_t); j += 8)
				ai |= (limb_t)(*a++) << j;
			wide_limb_t wide;
			cc_sub_borrow_wide_limb_limb(&wide, ai, p[i], borrow);
			borrow = wide.hi & 1;
		}

		return borrow;
	}
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_is_redc_mod_n<NLIMBS(256)>(const byte a[], const limb_t p[]);
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_is_redc_mod_n<NLIMBS(384)>(const byte a[], const limb_t p[]);

	/**
	 * @brief Adds two large integers modulo a given modulus and checks if the result is zero.
	 *
	 * Inputs are canonical little-endian byte strings in `[0,p)`. Stores
	 * `(a + b) mod p` and returns one exactly when the result is nonzero.
	 * Exact output aliasing with either input is supported.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE limb_t cc_add_n_check_mod_n(byte ret[], const byte a[], const byte b[],
		const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t zero;
		limb_t ret_[n];
		limb_t a_[n];
		limb_t b_[n];

		limbs_from_le_bytes(a_, a, n * sizeof(limb_t));
		limbs_from_le_bytes(b_, b, n * sizeof(limb_t));

		cc_add_mod_n<n>(ret_, a_, b_, p);
		zero = vec_is_zero(ret_, n * sizeof(limb_t));

		le_bytes_from_limbs(ret, ret_, n * sizeof(limb_t));

		return zero ^ 1;
	}
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_add_n_check_mod_n<NLIMBS(256)>(byte ret[], const byte a[], const byte b[],
		const limb_t p[]);
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_add_n_check_mod_n<NLIMBS(384)>(byte ret[], const byte a[], const byte b[],
		const limb_t p[]);

	/**
	 * @brief Subtracts two large integers modulo a given modulus and checks if the result is zero.
	 *
	 * Inputs are canonical little-endian byte strings in `[0,p)`. Stores
	 * `(a - b) mod p` and returns one exactly when the result is nonzero.
	 * Exact output aliasing with either input is supported.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE limb_t cc_sub_n_check_mod_n(byte ret[], const byte a[], const byte b[],
		const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t zero;
		limb_t ret_[n];
		limb_t a_[n];
		limb_t b_[n];

		limbs_from_le_bytes(a_, a, n * sizeof(limb_t));
		limbs_from_le_bytes(b_, b, n * sizeof(limb_t));

		cc_sub_mod_n<n>(ret_, a_, b_, p);
		zero = vec_is_zero(ret_, n * sizeof(limb_t));

		le_bytes_from_limbs(ret, ret_, n * sizeof(limb_t));

		return zero ^ 1;
	}
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_sub_n_check_mod_n<NLIMBS(256)>(byte ret[], const byte a[], const byte b[],
		const limb_t p[]);
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_sub_n_check_mod_n<NLIMBS(384)>(byte ret[], const byte a[], const byte b[],
		const limb_t p[]);

	/**
	 * @brief Converts a large integer from Montgomery form.
	 *
	 * Computes `a * R^(-1) mod p` to convert `a` from Montgomery form to standard form.
	 * The result is stored in `ret`. Assumes `n > 0` and `n` is even.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_from_mont_n(limb_t ret[], const limb_t a[],
		const limb_t p[], limb_t n0)
	{
		BLST_ASSUME(n != 0 && n % 2 == 0);
		wide_limb_t wide;
		limb_t mask, borrow, mx, hi;
		size_t i, j;
		limb_t tmp[n];

		for (j = 0; j < n; j++) {
			mx = n0 * a[0];
			cc_mul_add_wide_limb_limb(&wide, mx, p[0], a[0]);
			hi = wide.hi;

			for (i = 1; i < n; i++) {
				cc_mul_add_add_wide_limb_limb(&wide, mx, p[i], hi, a[i]);
				tmp[i - 1] = wide.lo;
				hi = wide.hi;
			}
			tmp[i - 1] = hi;
			a = tmp;
		}

		/* this is needed only if input can be non-fully-reduced */
		for (borrow = 0, i = 0; i < n; i++) {
			cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
			ret[i] = wide.lo;
			borrow = wide.hi & 1;
		}

		mask = 0 - borrow;
		BLST_launder(mask);

		for (i = 0; i < n; i++)
			ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);

	}
	template
	BLST_HD BLST_CPP_INLINE void cc_from_mont_n<NLIMBS(256)>(limb_t ret[], const limb_t a[],
		const limb_t p[], limb_t n0);
	template
	BLST_HD BLST_CPP_INLINE void cc_from_mont_n<NLIMBS(384)>(limb_t ret[], const limb_t a[],
		const limb_t p[], limb_t n0);

	/**
	 * @brief Performs Montgomery reduction on a large integer.
	 *
	 * Computes `ret = a / R mod p`, where `R = 2^(BLST_LIMB_SIZE * n)`,
	 * using the Montgomery reduction algorithm. This avoids explicit division
	 * by leveraging the precomputed Montgomery constant `n0`. The function
	 * iteratively reduces the input `a` (2n limbs) modulo `p` (n limbs) by
	 * calculating intermediate products and discarding the least significant
	 * limbs. After the reduction loop, it ensures the result is fully reduced
	 * by conditionally subtracting `p`. The result is stored in `ret` (n limbs).
	 * Assumes `n` is non-zero and even.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_redc_mont_n(limb_t ret[], const limb_t a[],
		const limb_t p[], limb_t n0)
	{
		BLST_ASSUME(n != 0 && n % 2 == 0);
		wide_limb_t wide;
		limb_t mask, carry, borrow, mx, hi;
		limb_t tmp[n];

		const limb_t* b = a;
		size_t i, j;

		for (j = 0; j < n; j++) {
			mx = n0 * b[0];
			cc_mul_add_wide_limb_limb(&wide, mx, p[0], b[0]);
			hi = wide.hi;
			for (i = 1; i < n; i++) {
				cc_mul_add_add_wide_limb_limb(&wide, mx, p[i], hi, b[i]);
				tmp[i - 1] = wide.lo;
				hi = wide.hi;
			}
			tmp[i - 1] = hi;
			b = tmp;
		}

		for (carry = 0, i = 0; i < n; i++) {
			cc_add_carry_wide_limb_limb(&wide, a[n + i], tmp[i], carry);
			tmp[i] = wide.lo;
			carry = wide.hi;
		}

		for (borrow = 0, i = 0; i < n; i++) {
			cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
			ret[i] = wide.lo;
			borrow = wide.hi & 1;
		}

		mask = carry - borrow;
		BLST_launder(mask);

		for (i = 0; i < n; i++)
			ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);
	}

	template
	BLST_HD BLST_CPP_INLINE void cc_redc_mont_n<NLIMBS(256)>(limb_t ret[], const limb_t a[],
		const limb_t p[], limb_t n0);
	template
	BLST_HD BLST_CPP_INLINE void cc_redc_mont_n<NLIMBS(384)>(limb_t ret[], const limb_t a[],
		const limb_t p[], limb_t n0);

	/**
	 * @brief Performs a right shift modulo operation on a large integer.
	 *
	 * Computes `a * 2^(-count) mod p`. The operand must be in `[0,p)`, and
	 * @p p must be odd. The result is in `[0,p)`. Exact output aliasing with
	 * @p a is supported; partial overlap is not.
	 *
	 * @param[out] ret The output array to store the result.
	 * @param[in]  a   The input array representing the large integer.
	 * @param[in]  count The number of bits to shift.
	 * @param[in]  p   The modulus array.
	 * @note Assumes `count > 0`, `n > 0`, and `n` is even. The selected
	 * boundary in `backend.h` defines a zero count as an exact copy.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_rshift_mod_n(limb_t ret[], const limb_t a[], size_t count,
		const limb_t p[])
	{

		BLST_ASSUME(count != 0);
		static_assert(n != 0 && n % 2 == 0);
		limb_t mask, carry, limb_, next;
		size_t i;

		while (count--) {
			mask = 0 - (a[0] & 1);
			BLST_launder(mask);
			for (carry = 0, i = 0; i < n; i++) {
				wide_limb_t wide;
				cc_add_carry_wide_limb_limb(&wide, a[i], p[i] & mask, carry);
				ret[i] = wide.lo;
				carry = wide.hi;
			}

			for (next = ret[0], i = 0; i < n - 1; i++) {
				limb_ = next >> 1;
				next = ret[i + 1];
				ret[i] = limb_ | next << (BLST_LIMB_SIZE - 1);
			}
			ret[i] = next >> 1 | carry << (BLST_LIMB_SIZE - 1);

			a = ret;
		}
	}

	template
	BLST_HD BLST_CPP_INLINE void cc_rshift_mod_n<NLIMBS(256)>(limb_t ret[], const limb_t a[], size_t count,
		const limb_t p[]);
	template
	BLST_HD BLST_CPP_INLINE void cc_rshift_mod_n<NLIMBS(384)>(limb_t ret[], const limb_t a[], size_t count,
		const limb_t p[]);


	/* @brief Compute sign and parity of a mod p
	 *
	 * Returns:
	 *  - bit 1: sign, 0 if a is zero or even, 1 if a is odd and positive,
	 *           2 if a is odd and negative
	 *  - bit 0: parity, 0 if a is even, 1 if a is odd
	 *
	 * Note that zero is considered even and non-negative.
	 *
	 * @param a   The input array representing the large integer.
	 * @param p   The modulus array.
	 * @param n   The number of limbs in `a` and `p`.
	 *
	 * @note Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE limb_t cc_sgn0_pty_mod_n(const limb_t a[], const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t carry, borrow, ret;
		limb_t tmp[n];

		size_t i;

		ret = a[0] & 1; /* parity */

		for (carry = 0, i = 0; i < n; i++) {
			limb_t a_i = a[i];
			tmp[i] = a_i << 1 | carry;
			carry = a_i >> (BLST_LIMB_SIZE - 1);
		}

		for (borrow = 0, i = 0; i < n; i++) {
			wide_limb_t wide;
			cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
			borrow = wide.hi & 1;
		}

		ret |= ((carry - borrow) & 2) ^ 2;

		return ret;
	}
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_sgn0_pty_mod_n<NLIMBS(256)>(const limb_t a[], const limb_t p[]);
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_sgn0_pty_mod_n<NLIMBS(384)>(const limb_t a[], const limb_t p[]);


	/*
	 * mul_mont_n without final conditional subtraction, which implies
	 * that modulus is one bit short, which in turn means that there are
	 * no carries to handle between iterations...
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_mul_mont_nonred_n(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[], limb_t n0)
	{
		BLST_ASSUME(n != 0 && n % 2 == 0);
		wide_limb_t wide;
		limb_t mx, hi;
		limb_t tmp[n + 1];

		size_t i, j;

		for (mx = b[0], hi = 0, i = 0; i < n; i++) {
			cc_mul_add_wide_limb_limb(&wide, mx, a[i], hi);
			tmp[i] = wide.lo;
			hi = wide.hi;
		}
		mx = n0 * tmp[0];
		tmp[i] = hi;

		for (j = 0; ; ) {
			cc_mul_add_wide_limb_limb(&wide, mx, p[0], tmp[0]);
			hi = wide.hi;

			for (i = 1; i < n; i++) {
				cc_mul_add_add_wide_limb_limb(&wide, mx, p[i], hi, tmp[i]);
				tmp[i - 1] = wide.lo;
				hi = wide.hi;
			}
			tmp[i - 1] = tmp[i] + hi;

			if (++j == n)
				break;

			for (mx = b[j], hi = 0, i = 0; i < n; i++) {
				cc_mul_add_add_wide_limb_limb(&wide, mx, a[i], hi, tmp[i]);
				tmp[i] = wide.lo;
				hi = wide.hi;
			}
			mx = n0 * tmp[0];
			tmp[i] = hi;

		}

		vec_copy(ret, tmp, n * sizeof(limb_t));
	}

	template
	BLST_HD BLST_CPP_INLINE void cc_mul_mont_nonred_n<NLIMBS(256)>(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[], limb_t n0);
	template
	BLST_HD BLST_CPP_INLINE void cc_mul_mont_nonred_n<NLIMBS(384)>(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[], limb_t n0);



	/**
	 * @brief Approximates two large integers for the extended Euclidean algorithm.
	 * https://eprint.iacr.org/2020/972 with 'k' being BLST_LIMB_SIZE-1.
	 *
	 * This function computes a 2-limb approximation of two large integers `a` and `b`.
	 * The approximation is performed by extracting the most significant limbs of `a` and `b`,
	 * aligning them, and storing the results in `a_` and `b_`. This is used to optimize
	 * the extended Euclidean algorithm by reducing the size of the operands.
	 *
	 * @param[out] a_ The 2-limb approximation of the first large integer.
	 * @param[in]  a  The first large integer.
	 * @param[out] b_ The 2-limb approximation of the second large integer.
	 * @param[in]  b  The second large integer.
	 * @param[in]  n  The number of limbs in `a` and `b`.
	 *
	 * @note Assumes `n > 0` and `n` is even.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_ab_approximation_n(limb_t a_[2], const limb_t a[],
		limb_t b_[2], const limb_t b[])
	{
		BLST_ASSUME(n != 0 && n % 2 == 0);
		limb_t a_hi, a_lo, b_hi, b_lo, mask;
		size_t i;

		i = n - 1;
		a_hi = a[i], a_lo = a[i - 1];
		b_hi = b[i], b_lo = b[i - 1];
		for (i--; --i;) {
			mask = 0 - is_zero(a_hi | b_hi);
			a_hi = ((a_lo ^ a_hi) & mask) ^ a_hi;
			b_hi = ((b_lo ^ b_hi) & mask) ^ b_hi;
			a_lo = ((a[i] ^ a_lo) & mask) ^ a_lo;
			b_lo = ((b[i] ^ b_lo) & mask) ^ b_lo;
		}
		i = BLST_LIMB_SIZE - cc_num_bits(a_hi | b_hi);
		/* |i| can be BLST_LIMB_SIZE if all a[2..]|b[2..] were zeros */

		a_[0] = a[0], a_[1] = cc_lshift_2(a_hi, a_lo, i);
		b_[0] = b[0], b_[1] = cc_lshift_2(b_hi, b_lo, i);
	}

	template
	BLST_HD BLST_CPP_INLINE void cc_ab_approximation_n<NLIMBS(256)>(limb_t a_[2], const limb_t a[],
		limb_t b_[2], const limb_t b[]);
	template
	BLST_HD BLST_CPP_INLINE void cc_ab_approximation_n<NLIMBS(384)>(limb_t a_[2], const limb_t a[],
		limb_t b_[2], const limb_t b[]);


	/**
	 * @brief Performs an inner loop operation for the extended Euclidean algorithm.
	 *
	 * This function performs a series of operations on two 2-limb integers `a_` and `b_`
	 * to compute updated coefficients `f0`, `g0`, `f1`, and `g1` for the extended Euclidean
	 * algorithm. The results are stored in the `factors` structure.
	 *
	 * @param[out] fg  The structure to store the updated coefficients.
	 * @param[in]  a_  The first 2-limb integer.
	 * @param[in]  b_  The second 2-limb integer.
	 * @param[in]  n   The number of iterations to perform.
	 *
	 * @note Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_inner_loop_n(factors* fg, const limb_t a_[2], const limb_t b_[2])
	{
		BLST_ASSUME(n != 0);
		wide_limb_t wide;
		limb_t f0 = 1, g0 = 0, f1 = 0, g1 = 1;
		limb_t a_lo, a_hi, b_lo, b_hi, t_lo, t_hi, odd, borrow, xorm;

		a_lo = a_[0], a_hi = a_[1];
		b_lo = b_[0], b_hi = b_[1];

		auto nn = n;
		while (nn--) {
			odd = 0 - (a_lo & 1);

			/* a_ -= b_ if a_ is odd */
			t_lo = a_lo, t_hi = a_hi;
			cc_sub_borrow_wide_limb_limb(&wide, a_lo, b_lo & odd, 0);
			a_lo = wide.lo;
			borrow = wide.hi & 1;

			cc_sub_borrow_wide_limb_limb(&wide, a_hi, (b_hi & odd), borrow);
			a_hi = wide.lo;
			borrow = wide.hi;

			/* negate a_-b_ if it borrowed */
			a_lo ^= borrow;
			a_hi ^= borrow;
			cc_add_wide_limb_limb(&wide, a_lo, borrow & 1);
			a_lo = wide.lo;
			a_hi += wide.hi & 1;

			/* b_=a_ if a_-b_ borrowed */
			b_lo = ((t_lo ^ b_lo) & borrow) ^ b_lo;
			b_hi = ((t_hi ^ b_hi) & borrow) ^ b_hi;

			/* exchange f0 and f1 if a_-b_ borrowed */
			xorm = (f0 ^ f1) & borrow;
			f0 ^= xorm;
			f1 ^= xorm;

			/* exchange g0 and g1 if a_-b_ borrowed */
			xorm = (g0 ^ g1) & borrow;
			g0 ^= xorm;
			g1 ^= xorm;

			/* subtract if a_ was odd */
			f0 -= f1 & odd;
			g0 -= g1 & odd;

			f1 <<= 1;
			g1 <<= 1;
			a_lo >>= 1; a_lo |= a_hi << (BLST_LIMB_SIZE - 1);
			a_hi >>= 1;
		}

		fg->f0 = f0, fg->g0 = g0, fg->f1 = f1, fg->g1 = g1;
	}
	template
	BLST_HD BLST_CPP_INLINE void cc_inner_loop_n<NLIMBS(256)>(factors* fg, const limb_t a_[2], const limb_t b_[2]);
	template
	BLST_HD BLST_CPP_INLINE void cc_inner_loop_n<NLIMBS(384)>(factors* fg, const limb_t a_[2], const limb_t b_[2]);


	/**
	 * @brief Conditionally negates a large integer.
	 *
	 * This function conditionally negates the large integer `a` based on the `neg` flag.
	 * If `neg` is non-zero, the result is `-a`; otherwise, the result is `a`. The result
	 * is stored in `ret`.
	 *
	 * @param[out] ret The output array to store the result.
	 * @param[in]  a   The large integer to be conditionally negated.
	 * @param[in]  neg The flag indicating whether to negate (`1`) or not (`0`).
	 * @param[in]  n   The number of limbs in `a`.
	 *
	 * @return An all-zero or all-one limb mask representing the result's sign.
	 *
	 * @note Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE limb_t cc_cneg_n(limb_t ret[], const limb_t a[], limb_t neg)
	{
		BLST_ASSUME(n != 0);
		wide_limb_t wide;
		limb_t carry;
		size_t i;
		wide.lo = 0;

		for (carry = neg & 1, i = 0; i < n; i++) {
			cc_add_wide_limb_limb(&wide, a[i] ^ neg, carry);
			ret[i] = wide.lo;
			carry = wide.hi;
		}

		return 0 - cc_msb(wide.lo);
	}
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_cneg_n<NLIMBS(256)>(limb_t ret[], const limb_t a[], limb_t neg);
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_cneg_n<NLIMBS(384)>(limb_t ret[], const limb_t a[], limb_t neg);


	/**
	 * @brief Adds two large integers and returns the carry.
	 *
	 * This function computes the sum of two large integers `a` and `b`, storing the result
	 * in `ret`. It also returns the carry resulting from the addition.
	 *
	 * @param[out] ret The output array to store the result.
	 * @param[in]  a   The first large integer.
	 * @param[in]  b   The second large integer.
	 * @param[in]  n   The number of limbs in `a` and `b`.
	 *
	 * @return The carry resulting from the addition.
	 *
	 * @note Assumes `n > 0`.
	 */
	template<int n>
	BLST_HD BLST_CPP_INLINE limb_t cc_add_n(limb_t ret[], const limb_t a[], limb_t b[])
	{
		BLST_ASSUME(n != 0);
		wide_limb_t wide;
		limb_t carry;
		size_t i;

		for (carry = 0, i = 0; i < n; i++) {
			cc_add_carry_wide_limb_limb(&wide, a[i], b[i], carry);
			ret[i] = wide.lo;
			carry = wide.hi;
		}

		return carry;
	}
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_add_n<NLIMBS(256)>(limb_t ret[], const limb_t a[], limb_t b[]);
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_add_n<NLIMBS(384)>(limb_t ret[], const limb_t a[], limb_t b[]);

	/**
	 * @brief Multiplies a large integer by a single limb and accumulates the result.
	 *
	 * This function computes the product of a large integer `a` and a single limb `b`,
	 * while adding any carry from the previous operation. The result is stored in `ret`,
	 * and the final carry is returned.
	 *
	 * @param[out] ret The output array to store the result.
	 * @param[in]  a   The large integer to be multiplied.
	 * @param[in]  b   The single limb multiplier.
	 * @param[in]  n   The number of limbs in `a`.
	 *
	 * @return The final carry after the multiplication and addition.
	 *
	 * @note Assumes `n > 0`.
	 */
	template<int n>
	BLST_HD BLST_CPP_INLINE limb_t cc_umul_n(limb_t ret[], const limb_t a[], limb_t b)
	{
		BLST_ASSUME(n != 0);
		wide_limb_t wide;
		limb_t hi;
		size_t i;

		for (hi = 0, i = 0; i < n; i++) {
			cc_mul_add_wide_limb_limb(&wide, b, a[i], hi);
			ret[i] = wide.lo;
			hi = wide.hi;
		}

		return hi;
	}
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_umul_n<NLIMBS(256)>(limb_t ret[], const limb_t a[], limb_t b);
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_umul_n<NLIMBS(384)>(limb_t ret[], const limb_t a[], limb_t b);


	/**
	 * @brief Computes a signed modular linear combination of two large integers, followed by a right shift.
	 *
	 * This function calculates `(a * |f| + b * |g|) >> k`, where `k` is `BLST_LIMB_SIZE - 2`.
	 * The coefficients `f` and `g` are adjusted to be positive, and the corresponding integers
	 * `a` and `b` are conditionally negated if their coefficients are negative. The result is
	 * stored in `ret`, and the signs of `f` and `g` are updated to reflect the adjustments.
	 *
	 * @param[out] ret The output array to store the result.
	 * @param[in]  a   The first large integer.
	 * @param[in,out] f_ Pointer to the signed coefficient for `a`, updated to reflect its sign adjustment.
	 * @param[in]  b   The second large integer.
	 * @param[in,out] g_ Pointer to the signed coefficient for `b`, updated to reflect its sign adjustment.
	 * @param[in]  n   The number of limbs in `a` and `b`.
	 *
	 * @return Returns the sign of the result (0 for non-negative, 1 for negative).
	 *
	 * @note Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE limb_t cc_smul_n_shift_n(limb_t ret[], const limb_t a[], limb_t* f_,
		const limb_t b[], limb_t* g_)
	{
		BLST_ASSUME(n != 0);
		limb_t f, g, neg, carry, hi;
		limb_t a_[n + 1];
		limb_t b_[n + 1];
		size_t i;

		/* |a|*|f_| */
		f = *f_;
		neg = 0 - cc_msb(f);
		f = (f ^ neg) - neg;            /* ensure |f| is positive */
		(void)cc_cneg_n<n>(a_, a, neg);
		hi = cc_umul_n<n>(a_, a_, f);
		a_[n] = hi - (f & neg);

		/* |b|*|g_| */
		g = *g_;
		neg = 0 - cc_msb(g);
		g = (g ^ neg) - neg;            /* ensure |g| is positive */
		(void)cc_cneg_n<n>(b_, b, neg);
		hi = cc_umul_n<n>(b_, b_, g);
		b_[n] = hi - (g & neg);

		/* |a|*|f_| + |b|*|g_| */
		(void)cc_add_n<n + 1>(a_, a_, b_);

		/* (|a|*|f_| + |b|*|g_|) >> k */
		for (carry = a_[0], i = 0; i < n; i++) {
			hi = carry >> (BLST_LIMB_SIZE - 2);
			carry = a_[i + 1];
			ret[i] = hi | (carry << 2);
		}

		/* ensure result is non-negative, fix up |f_| and |g_| accordingly */
		neg = 0 - cc_msb(carry);
		*f_ = (*f_ ^ neg) - neg;
		*g_ = (*g_ ^ neg) - neg;
		(void)cc_cneg_n<n>(ret, ret, neg);

		return neg;
	}
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_smul_n_shift_n<NLIMBS(256)>(limb_t ret[], const limb_t a[], limb_t* f_,
		const limb_t b[], limb_t* g_);
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_smul_n_shift_n<NLIMBS(384)>(limb_t ret[], const limb_t a[], limb_t* f_,
		const limb_t b[], limb_t* g_);

	/* @brief Computes |u|*|f_| + |v|*|g_|
	 *
	 * Assumes n>0.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE limb_t cc_smul_2n(limb_t ret[], const limb_t u[], limb_t f,
		const limb_t v[], limb_t g)
	{
		BLST_ASSUME(n != 0);
		limb_t neg, hi;

		limb_t u_[n];
		limb_t v_[n];

		/* |u|*|f_| */
		neg = 0 - cc_msb(f);
		f = (f ^ neg) - neg;            /* ensure |f| is positive */
		neg = cc_cneg_n<n>(u_, u, neg);
		hi = cc_umul_n<n>(u_, u_, f) - (f & neg);

		/* |v|*|g_| */
		neg = 0 - cc_msb(g);
		g = (g ^ neg) - neg;            /* ensure |g| is positive */
		neg = cc_cneg_n<n>(v_, v, neg);
		hi += cc_umul_n<n>(v_, v_, g) - (g & neg);

		/* |u|*|f_| + |v|*|g_| */
		hi += cc_add_n<n>(ret, u_, v_);

		return hi;
	}
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_smul_2n<NLIMBS(256)>(limb_t ret[], const limb_t u[], limb_t f,
		const limb_t v[], limb_t g);
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_smul_2n<NLIMBS(384)>(limb_t ret[], const limb_t u[], limb_t f,
		const limb_t v[], limb_t g);


	/*
	 * Fixed-schedule double-width almost-inverse. The complete contract is at
	 * the declaration in backend_portable.h.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE void cc_ct_inverse_mod_n(limb_t ret[], const limb_t inp[],
		const limb_t mod[], const limb_t modx[])
	{
		BLST_ASSUME(n != 0 && n % 2 == 0);
		wide_limb_t wide;
		limb_t a[n];
		limb_t b[n];
		limb_t u[2 * n];
		limb_t v[2 * n];
		limb_t t[2 * n];

		limb_t a_[2], b_[2], sign, carry, top;
		factors fg;
		size_t i;

		vec_copy(a, inp, n * sizeof(limb_t));
		vec_copy(b, mod, n * sizeof(limb_t));
		vec_zero(u, n * 2 * sizeof(limb_t)); u[0] = 1;
		vec_zero(v, n * 2 * sizeof(limb_t));

		for (i = 0; i < (2 * n * BLST_LIMB_SIZE) / (BLST_LIMB_SIZE - 2); i++) {
			cc_ab_approximation_n<n>(a_, a, b_, b);

			cc_inner_loop_n<BLST_LIMB_SIZE - 2>(&fg, a_, b_);

			(void)cc_smul_n_shift_n<n>(t, a, &fg.f0, b, &fg.g0);
			(void)cc_smul_n_shift_n<n>(b, a, &fg.f1, b, &fg.g1);
			vec_copy(a, t, n * sizeof(limb_t));
			cc_smul_2n<2 * n>(t, u, fg.f0, v, fg.g0);
			cc_smul_2n<2 * n>(v, u, fg.f1, v, fg.g1);
			vec_copy(u, t, 2 * n * sizeof(limb_t));
		}

		cc_inner_loop_n<(2 * n * BLST_LIMB_SIZE) % (BLST_LIMB_SIZE - 2)>(&fg, a, b);
		top = cc_smul_2n<2 * n>(ret, u, fg.f1, v, fg.g1);

		sign = 0 - cc_msb(top);    /* top is 1, 0 or -1 */
		for (carry = 0, i = 0; i < n; i++) {
			cc_add_carry_wide_limb_limb(&wide, ret[n + i], modx[i] & sign, carry);
			ret[n + i] = wide.lo;
			carry = wide.hi;
		}
		top += carry;
		sign = 0 - top;         /* top is 1, 0 or -1 */
		top |= sign;
		for (i = 0; i < n; i++)
			a[i] = modx[i] & top;
		(void)cc_cneg_n<n>(a, a, 0 - cc_msb(sign));
		cc_add_n<n>(ret + n, ret + n, a);
	}
	template
	BLST_HD BLST_CPP_INLINE void cc_ct_inverse_mod_n<NLIMBS(256)>(limb_t ret[], const limb_t inp[],
		const limb_t mod[], const limb_t modx[]);
	template
	BLST_HD BLST_CPP_INLINE void cc_ct_inverse_mod_n<NLIMBS(384)>(limb_t ret[], const limb_t inp[],
		const limb_t mod[], const limb_t modx[]);

	/*
	 * Copy of inner_loop_n above, but with |L| updates.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE limb_t cc_legendre_loop_n(limb_t L, factors* fg, const limb_t a_[2],
		const limb_t b_[2])
	{
		BLST_ASSUME(n != 0);
		wide_limb_t wide;
		limb_t f0 = 1, g0 = 0, f1 = 0, g1 = 1;
		limb_t a_lo, a_hi, b_lo, b_hi, t_lo, t_hi, odd, borrow, xorm;

		a_lo = a_[0], a_hi = a_[1];
		b_lo = b_[0], b_hi = b_[1];

		auto nn = n;
		while (nn--) {
			odd = 0 - (a_lo & 1);

			/* a_ -= b_ if a_ is odd */
			t_lo = a_lo, t_hi = a_hi;
			cc_sub_borrow_wide_limb_limb(&wide, a_lo, b_lo & odd, 0);
			a_lo = wide.lo;
			borrow = wide.hi & 1;
			cc_sub_borrow_wide_limb_limb(&wide, a_hi, (b_hi & odd), borrow);

			a_hi = wide.lo;
			borrow = wide.hi;

			L += ((t_lo & b_lo) >> 1) & borrow;

			/* negate a_-b_ if it borrowed */
			a_lo ^= borrow;
			a_hi ^= borrow;
			cc_add_wide_limb_limb(&wide, a_lo, borrow & 1);
			a_lo = wide.lo;
			a_hi += wide.hi & 1;

			/* b_=a_ if a_-b_ borrowed */
			b_lo = ((t_lo ^ b_lo) & borrow) ^ b_lo;
			b_hi = ((t_hi ^ b_hi) & borrow) ^ b_hi;

			/* exchange f0 and f1 if a_-b_ borrowed */
			xorm = (f0 ^ f1) & borrow;
			f0 ^= xorm;
			f1 ^= xorm;

			/* exchange g0 and g1 if a_-b_ borrowed */
			xorm = (g0 ^ g1) & borrow;
			g0 ^= xorm;
			g1 ^= xorm;

			/* subtract if a_ was odd */
			f0 -= f1 & odd;
			g0 -= g1 & odd;

			f1 <<= 1;
			g1 <<= 1;
			a_lo >>= 1; a_lo |= a_hi << (BLST_LIMB_SIZE - 1);
			a_hi >>= 1;

			L += (b_lo + 2) >> 2;
		}

		fg->f0 = f0, fg->g0 = g0, fg->f1 = f1, fg->g1 = g1;

		return L;
	}
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_legendre_loop_n<NLIMBS(256)>(limb_t L, factors* fg, const limb_t a_[2],
		const limb_t b_[2]);
	template
	BLST_HD BLST_CPP_INLINE limb_t cc_legendre_loop_n<NLIMBS(384)>(limb_t L, factors* fg, const limb_t a_[2],
		const limb_t b_[2]);

	/* Fixed-schedule residue test; see backend_portable.h for the contract. */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE bool_t cc_ct_is_sqr_mod_n(const limb_t inp[], const limb_t mod[])
	{
		BLST_ASSUME(n != 0 && n % 2 == 0);
		limb_t a[n];
		limb_t b[n];
		limb_t t[n];
		limb_t a_[2], b_[2], neg, L = 0;
		factors fg;
		size_t i;

		vec_copy(a, inp, n * sizeof(limb_t));
		vec_copy(b, mod, n * sizeof(limb_t));

		for (i = 0; i < (2 * n * BLST_LIMB_SIZE) / (BLST_LIMB_SIZE - 2); i++) {
			cc_ab_approximation_n<n>(a_, a, b_, b);

			L = cc_legendre_loop_n<BLST_LIMB_SIZE - 2>(L, &fg, a_, b_);

			neg = cc_smul_n_shift_n<n>(t, a, &fg.f0, b, &fg.g0);
			(void)cc_smul_n_shift_n<n>(b, a, &fg.f1, b, &fg.g1);
			vec_copy(a, t, n * sizeof(limb_t));
			L += (b[0] >> 1) & neg;
		}

		L = cc_legendre_loop_n<(2 * n * BLST_LIMB_SIZE) % (BLST_LIMB_SIZE - 2)>(L, &fg, a, b);

		return (L & 1) ^ 1;
	}
	template
	BLST_HD BLST_CPP_INLINE bool_t cc_ct_is_sqr_mod_n<NLIMBS(256)>(const limb_t inp[], const limb_t mod[]);
	template
	BLST_HD BLST_CPP_INLINE bool_t cc_ct_is_sqr_mod_n<NLIMBS(384)>(const limb_t inp[], const limb_t mod[]);


	/**
	 * @brief Squares a 384-bit integer multiple times and multiplies the result in Montgomery form.
	 *
	 * This function computes `(a^(2^count) * b * R^(-1)) mod p` in Montgomery form, where `R` is the Montgomery constant.
	 * The input `a` is squared `count` times, and the result is multiplied by `b` modulo `p`.
	 *
	 * @param[out] ret   The output 384-bit integer in Montgomery form.
	 * @param[in]  a     The input 384-bit integer in Montgomery form.
	 * @param[in]  count The number of squaring operations to perform.
	 * @param[in]  p     The modulus.
	 * @param[in]  n0    The Montgomery constant.
	 * @param[in]  b     The multiplier in Montgomery form.
	 *
	 * @note Assumes `count > 0`; `ret` may alias `a` but must not alias `b`.
	 */
	BLST_HD BLST_CPP_INLINE void cc_sqr_n_mul_mont_384(
		limb_t* ret, const limb_t* a, size_t count,
		const limb_t* p, limb_t n0, const limb_t* b)
	{
		BLST_ASSUME(count != 0);
		BLST_ASSERT(ret != b);
		auto ptr = a;

		while (count--) {
			cc_mul_mont_nonred_n<NLIMBS(384)>(ret, ptr, ptr, p, n0);
			ptr = ret;
		}
		cc_mul_mont_n<NLIMBS(384)>(ret, ret, b, p, n0);
	}


	/*
	 * @brief Divides a 3-limb_t dividend by a 2-limb divisor.
	 *
	 * Estimates one quotient limb for `(div_top * 2^64) / divisor`.
	 * The divisor is normalized and `div_top <= divisor`. Equality saturates
	 * an overflowing quotient to `UINT64_MAX`. Long division supplies the
	 * stricter `div_top < divisor` condition.
	 *
	 * |div_top| points at two most significant limbs of the dividend, |d_hi|
	 * and |d_lo| are two most significant limbs of the divisor. If divisor
	 * is only one limb, it is to be passed in |d_hi| with zero in |d_lo|.
	 * The divisor is required to be "bitwise left-aligned," and dividend's
	 * top limbs to be not larger than the divisor's. The latter limitation
	 * can be problematic in the first iteration of multi-precision division,
	 * where in most general case the condition would have to be "smaller."
	 * The subroutine considers four limbs, two of which are "overlapping,"
	 * hence the name... Another way to look at it is to think of the pair
	 * of the dividend's limbs being suffixed with a zero:
	 *   +-------+-------+-------+
	 * R |       |       |   0   |
	 *   +-------+-------+-------+
	 *           +-------+-------+
	 * D         |       |       |
	 *           +-------+-------+
	 */
	BLST_HD BLST_CPP_INLINE limb_t cc_div_3_limbs(const limb_t div_top[2], limb_t d_lo, limb_t d_hi)
	{
		wide_limb_t Rx;
		limb_t r_lo = div_top[0], r_hi = div_top[1];
		limb_t Q = 0, mask, borrow, rx;
		size_t i;

		for (i = 0; i < BLST_LIMB_SIZE; i++) {
			/* "borrow, Rx = R - D" */
			cc_sub_borrow_wide_limb_limb(&Rx, r_lo, d_lo, 0);
			rx = Rx.lo;
			borrow = Rx.hi & 1;

			cc_sub_borrow_wide_limb_limb(&Rx, r_hi, d_hi, borrow);
			borrow = Rx.hi;// not & 1. want 0000 or 1111


			/* "if (R >= D) R -= D" */
			r_lo = ((r_lo ^ rx) & borrow) ^ rx;
			rx = Rx.lo;
			r_hi = ((r_hi ^ rx) & borrow) ^ rx;

			Q <<= 1;
			Q |= ~borrow & 1;

			/* "D >>= 1" */
			d_lo >>= 1; d_lo |= d_hi << (BLST_LIMB_SIZE - 1);
			d_hi >>= 1;
		}

		mask = 0 - cc_msb(Q);  /* does it overflow? */

		/* "borrow, Rx = R - D" */
		cc_sub_borrow_wide_limb_limb(&Rx, r_lo, d_lo, 0);
		borrow = Rx.hi & 1;

		cc_sub_borrow_wide_limb_limb(&Rx, r_hi, d_hi, borrow);
		borrow = Rx.hi & 1;

		Q <<= 1;
		Q |= borrow ^ 1;

		return (Q | mask);
	}

	/**
	 * @brief Computes the quotient and remainder of a large integer division.
	 *
	 * The caller supplies an exact or one-high quotient estimate for the
	 * `(n+1)`-limb dividend in `div_rem`. The low `n` limbs receive the
	 * remainder; `div_rem[n]` and the return value receive the corrected quotient.
	 */
	template<size_t n>
	BLST_HD BLST_CPP_INLINE limb_t cc_quot_rem_n(limb_t* div_rem, const limb_t* divisor,
		limb_t quotient)
	{
		wide_limb_t wide;
		limb_t carry, mask, borrow;
		limb_t tmp[n + 1];
		size_t i;

		/* divisor*quotient */
		for (carry = 0, i = 0; i < n; i++) {
			cc_mul_add_wide_limb_limb(&wide, quotient, divisor[i], carry);
			tmp[i] = wide.lo;
			carry = wide.hi;
		}
		tmp[i] = carry;

		/* remainder = dividend - divisor*quotient */
		for (borrow = 0, i = 0; i <= n; i++) {
			cc_sub_borrow_wide_limb_limb(&wide, div_rem[i], tmp[i], borrow);
			tmp[i] = wide.lo;
			borrow = wide.hi & 1;
		}

		mask = 0 - borrow;
		BLST_launder(mask);

		/* if quotient was off by one, add divisor to the remainder */
		for (carry = 0, i = 0; i < n; i++) {
			cc_add_carry_wide_limb_limb(&wide, tmp[i], divisor[i] & mask, carry);
			div_rem[i] = wide.lo;
			carry = wide.hi & 1;
		}

		return (div_rem[i] = quotient + mask);
	}
	template
		BLST_HD BLST_CPP_INLINE limb_t cc_quot_rem_n<NLIMBS(128)>(limb_t* div_rem, const limb_t* divisor,
			limb_t quotient);
	template
		BLST_HD BLST_CPP_INLINE limb_t cc_quot_rem_n<NLIMBS(64)>(limb_t* div_rem, const limb_t* divisor,
			limb_t quotient);

}
