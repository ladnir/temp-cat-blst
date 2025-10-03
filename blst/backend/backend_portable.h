/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_NO_ASM_H__
#define __BLS12_NO_ASM_H__

#include "blst/common/bytes.h"
#include "blst/config.h"
#include "blst/common/defines.h"
#include <assert.h>
#include <stddef.h>   // for size_t
 //printf
#include <stdio.h>


#ifdef BLST_ADX
#include <immintrin.h>
#endif


#if defined(__clang__)
# pragma GCC diagnostic ignored "-Wstatic-in-inline"
#endif

//#ifdef __CUDACC__
//#include <cuda/std/array>
//template <typename T, size_t N>
//using blst_array = cuda::std::array<T, N>;;
//#else
//#include <array>
//template <typename T, size_t N>
//using blst_array = std::array<T, N>;
//#endif

namespace blst
{


	//BLST_EXTERN_C_BEGIN
	typedef struct { limb_t lo, hi; } wide_limb_t;

#if BLST_LIMB_SIZE == 64 && defined(BLST_ADX) && !defined(__CUDA_ARCH__)
//	// 64 bit with wide operations
//#define cc_mul_wide_limb_limb(ret, a, b)\
//    (ret)->lo = _mulx_u64(a, b, &(ret)->hi)
//
//#define cc_add_wide_wide_limb(ret, a, b)\
//    (ret)->hi = a.hi + _addcarry_u64(0, (a).lo, b, (limb_t*) &(ret)->lo)
//
//	BLST_HD BLST_FORCEINLINE void  cc_add_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b)
//	{
//		ret->hi = _addcarry_u64(0, a, b, &ret->lo);
//	}
//
//	BLST_HD BLST_FORCEINLINE void cc_add_carry_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b, const limb_t& carry)
//	{
//		assert(carry <= 1);
//		BLST_ASSUME(carry <= 1);
//		(ret)->hi = _addcarry_u64((unsigned char)carry, a, b, (limb_t*)&(ret)->lo);
//	}
//
//#define cc_sub_borrow_wide_limb_limb(ret, a, b, borrow) \
//    do{ assert(borrow <= 1); BLST_ASSUME(borrow <=1); \
//    (ret)->hi = 0-(limb_t)_subborrow_u64((unsigned char)borrow, a, b, (limb_t*) &(ret)->lo); } while(0)

		// 64 bit with wide operations
#define cc_mul_wide_limb_limb(ret, a, b)\
    (ret)->lo = _mulx_u64(a, b, &(ret)->hi)

#define cc_add_wide_wide_limb(ret, a, b)\
    (ret)->hi = a.hi + _addcarry_u64(0, (a).lo, b, (limb_t*) &(ret)->lo)

#define cc_add_wide_limb_limb(ret, a, b)\
    (ret)->hi = _addcarry_u64(0, a, b, (limb_t*) &(ret)->lo)

#define cc_add_carry_wide_limb_limb(ret, a, b, carry) \
    do{ BLST_ASSUME(carry <=1); \
    (ret)->hi = _addcarry_u64((unsigned char)carry, a, b, (limb_t*) &(ret)->lo); } while(0)

#define cc_sub_borrow_wide_limb_limb(ret, a, b, borrow) \
    do{ BLST_ASSUME(borrow <=1); \
    (ret)->hi = 0-(limb_t)_subborrow_u64((unsigned char)borrow, a, b, (limb_t*) &(ret)->lo); } while(0)


	//#define cc_sub_wide_limb_wide(ret, a, b)\
	//    do {												\
	//    limb lo_out, hi_out;								\
	//	unsigned char bor = _subborrow_u64(0, a,			\
	//		(b).lo, &lo_out);								\
	//	(void)_subborrow_u64(bor, 0ULL,						\
	//		(b).hi, &hi_out);								\
	//	(ret)->lo = lo_out;									\
	//	(ret)->hi = hi_out;									\
	//	} while (0)

#elif BLST_LIMB_SIZE == 32 && (defined(_M_X64) || defined(__x86_64__))
	// 32 bit with 64 bit support

	static_assert(&((wide_limb_t*)0)->lo < &((wide_limb_t*)0)->hi,
		"this implementation assumes we can cast to 64 bit");
#define cc_mul_wide_limb_limb(ret, a, b)\
    do{*(unsigned long long*)(ret) = (unsigned long long)(a) * (unsigned long long)(b);}while(0)

#define cc_add_wide_wide_limb(ret, a, b)\
    do{*(unsigned long long*)(ret) = *(unsigned long long*)&(a) + (unsigned long long)(b);}while(0)

#define cc_add_wide_limb_limb(ret, a, b)\
    do{*(unsigned long long*)(ret) = (a) + (unsigned long long)(b);}while(0)

#define cc_sub_wide_limb_wide(ret, a, b)\
    do{*(unsigned long long*)(ret) = (a) - *(unsigned long long*)(&(b));}while(0)

#define cc_add_carry_wide_limb_limb(ret, a, b, carry) \
    do{ assert(carry <= 1); BLST_ASSUME(carry <=1); \
    *(unsigned long long*)(ret) = (a) + (unsigned long long)(b) + carry; } while(0)

#define cc_sub_borrow_wide_limb_limb(ret, a, b, borrow) \
    do{ assert(borrow <= 1); BLST_ASSUME(borrow <=1); \
    *(unsigned long long*)(ret) = (a) - (unsigned long long)(b) - borrow; } while(0)

#else
	// generic

#define cc_mul_wide_limb_limb(ret, a, b)\
    do{int half = sizeof(limb_t) * 4;                                  \
    limb_t low_mask = ((limb_t)1 << half) - 1;                         \
    /* Portable C implementation for 64-bit limbs*/                    \
    /* Break operands into 32-bit halves        */                     \
    limb_t a_lo = (a) & low_mask;                                      \
    limb_t a_hi = (a) >> half;                                         \
    limb_t b_lo = (b) & low_mask;                                      \
    limb_t b_hi = (b) >> half;                                         \
                                                                       \
    /* Four partial products: */                                       \
    limb_t p0 = a_lo * b_lo;                                           \
    limb_t p1 = a_lo * b_hi;                                           \
    limb_t p2 = a_hi * b_lo;                                           \
    limb_t p3 = a_hi * b_hi;                                           \
                                                                       \
    /* Combine middle products with carry propagation */               \
    limb_t middle = (p0 >> half) + (p1 & low_mask) + (p2 & low_mask);  \
    limb_t carry = middle >> half;                                     \
                                                                       \
    /* Assemble final result    */                                     \
    (ret)->lo = (p0 & low_mask) + (middle << half);                    \
    (ret)->hi = p3 + (p1 >> half) + (p2 >> half) + carry;}while(0)


	BLST_HD BLST_FORCEINLINE void  cc_add_wide_wide_limb(wide_limb_t* ret, const  wide_limb_t& a, const limb_t& b)
	{
		wide_limb_t aa = a;
		(ret)->lo = aa.lo + (b);
		(ret)->hi = aa.hi + (((ret)->lo < aa.lo) ? 1 : 0);
	}

	BLST_HD BLST_FORCEINLINE void cc_add_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b)
	{
		limb_t aa = a;
		(ret)->lo = (a)+(b);
		(ret)->hi = ((ret)->lo < aa) ? 1 : 0;
	}

	// wide subtraction *ret = a - b
	// with wide a and limb_t b
	BLST_HD BLST_FORCEINLINE void  cc_sub_wide_wide_limb(wide_limb_t* ret, const  wide_limb_t& a, const  limb_t& b)
	{
		wide_limb_t _aa = (a);
		limb_t      _bb = (limb_t)(b);
		limb_t      _old = _aa.lo;
		(ret)->lo = (limb_t)(_old - _bb);
		limb_t _borrow = (limb_t)(_old < _bb);    /* correct borrow test */
		(ret)->hi = (limb_t)(_aa.hi - _borrow);
	}

	// wide subtraction *ret = a - b
	// with limb_t a and b
	BLST_HD BLST_FORCEINLINE void  cc_sub_wide_limb_limb(wide_limb_t* ret, const  limb_t& a, const  limb_t& b)
	{
		limb_t _aa = (limb_t)(a);
		limb_t _bb = (limb_t)(b);
		(ret)->lo = (limb_t)(_aa - _bb);
		limb_t _borrow = (limb_t)(_aa < _bb);
		(ret)->hi = (limb_t)(-(limb_t)_borrow); /* 0 or 0xFFFFFFFF */
	}

	// ret = a + b + c
	BLST_HD BLST_FORCEINLINE void cc_add_carry_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b, const limb_t& c)
	{
		cc_add_wide_limb_limb(ret, a, b);
		cc_add_wide_wide_limb(ret, *ret, c);
	}

	// ret = a - b - c
	BLST_HD BLST_FORCEINLINE void cc_sub_borrow_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b, const limb_t& c)
	{
		cc_sub_wide_limb_limb(ret, a, b);
		cc_sub_wide_wide_limb(ret, *ret, c);
	}

	//#define cc_sub_wide_limb_wide(ret, a, b)\
	//    do {												\
	//    limb lo = (a) - (b).lo;							\
	//	limb bor = ((a) < (b).lo);						\
	//	limb hi = (limb(0) - (b).hi) - bor;				\
	//	(ret)->lo = lo;										\
	//	(ret)->hi = hi;										\
	//	} while (0)

#endif

// ret = a * b + c
#define cc_mul_add_wide_limb_limb(ret, a, b, c)\
    cc_mul_wide_limb_limb(ret, a, b);\
    cc_add_wide_wide_limb(ret, *ret, c)

// ret = a * b + c + d
#define cc_mul_add_add_wide_limb_limb(ret, a, b, c, d)\
    cc_mul_wide_limb_limb(ret, a, b);\
    cc_add_wide_wide_limb(ret, *ret, c);\
    cc_add_wide_wide_limb(ret, *ret, d)

// ret = a + b + c
#define cc_add_add_wide_limb_limb(ret, a, b, c)\
    cc_add_wide_limb_limb(ret, a, b);\
    cc_add_wide_wide_limb(ret, *ret, c)



	BLST_INLINE void print_vec_n(const char* label, const limb_t* v, int n) {
		printf("%s: ", label);
		for (int i = n - 1; i >= 0; --i) {
			printf("%016llx", v[i]);
			if (i > 0) printf(" ");
		}
		printf("\n");
	}


	/**
	 * @brief Performs Montgomery multiplication of two large integers.
	 *
	 * Computes `(a * b * R^(-1)) mod p` using the Montgomery multiplication algorithm,
	 * where `R = 2^(BLST_LIMB_SIZE * n)`. The result is stored in `ret`. Assumes `n > 0`
	 * and `n` is even.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE void cc_mul_mont_n(limb_t* ret, const limb_t* a, const limb_t* b,
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
				// wide = mx * p + hi + tmp;
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

		//BLST_VLA_FREE(tmp);
	}
	//
	//#define MUL_MONT_IMPL(fn1, fn2, bits) \
	//BLST_HD BLST_INLINE void fn1(vec##bits ret, const vec##bits a, \
	//                            const vec##bits b, const vec##bits p, limb n0) \
	//{   cc_mul_mont_n(ret, a, b, p, n0, NLIMBS(bits));   } \
	//\
	//BLST_HD BLST_INLINE void fn2(vec##bits ret, const vec##bits a, \
	//                            const vec##bits p, limb n0) \
	//{   cc_mul_mont_n(ret, a, a, p, n0, NLIMBS(bits));   }

	/*
	 * 256-bit subroutines can handle arbitrary modulus, even non-"sparse",
	 * but we have to harmonize the naming with assembly.
	 */
	 //MUL_MONT_IMPL(cc_mul_mont_sparse_256, cc_sqr_mont_sparse_256, 256)

	//MUL_MONT_IMPL(cc_mul_mont_384, cc_sqr_mont_384, 384)


	template<size_t n>
	BLST_HD BLST_INLINE void cc_add_mod_n(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t mask, carry, borrow;
		//BLST_VLA_NEW(limb_t, tmp, n);
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

		//BLST_VLA_FREE(tmp);
	}

	//#define ADD_MOD_IMPL(bits) \
	//BLST_HD BLST_INLINE void cc_add_mod_##bits(vec##bits ret, const vec##bits a, \
	//                           const vec##bits b, const vec##bits p) \
	//{   cc_add_mod_n(ret, a, b, p, NLIMBS(bits));   }
	//ADD_MOD_IMPL(256)
	//ADD_MOD_IMPL(384)



	/**
	 * @brief Subtracts two large integers modulo a given modulus.
	 *
	 * Computes `(a - b) mod p` and stores the result in `ret`. Handles borrow propagation
	 * and ensures the result is reduced modulo `p`. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE void cc_sub_mod_n(limb_t ret[], const limb_t a[], const limb_t b[],
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

	//#define SUB_MOD_IMPL(bits) \
	//BLST_HD BLST_INLINE void cc_sub_mod_##bits(vec##bits ret, const vec##bits a, \
	//                           const vec##bits b, const vec##bits p) \
	//{   cc_sub_mod_n(ret, a, b, p, NLIMBS(bits));   }
	//
	//SUB_MOD_IMPL(256)
	//SUB_MOD_IMPL(384)


	/**
	 * @brief Multiplies a large integer by 3 modulo a given modulus.
	 *
	 * Computes `(3 * a) mod p` and stores the result in `ret`. Handles carry propagation
	 * and ensures the result is reduced modulo `p`. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE void cc_mul_by_3_mod_n(limb_t ret[], const limb_t a[], const limb_t p[])
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


		//BLST_VLA_FREE(tmp);
		//BLST_VLA_FREE(two_a);
	}

	//#define MUL_BY_3_MOD_IMPL(bits) \
	//BLST_HD BLST_INLINE void cc_mul_by_3_mod_##bits(vec##bits ret, const vec##bits a, \
	//                                const vec##bits p) \
	//{   cc_mul_by_3_mod_n(ret, a, p, NLIMBS(bits));   }
	//
	//MUL_BY_3_MOD_IMPL(256)
	//MUL_BY_3_MOD_IMPL(384)


	/**
	 * @brief Performs a left shift modulo operation on a large integer.
	 *
	 * Computes `(a << count) mod p` iteratively, ensuring the result remains within
	 * the range `[0, p)`. Assumes `count > 0` and `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE void cc_lshift_mod_n(limb_t ret[], const limb_t a[], size_t count,
		const limb_t p[])
	{
#if 0
		auto ptr = a;
		while (count--)
		{
			cc_add_mod_n<NLIMBS(384)>(ret, ptr, ptr, p);
			ptr = ret;
		}

#else
		BLST_ASSERT(count != 0);
		BLST_ASSUME(count != 0);
		BLST_ASSUME(n != 0);
		limb_t mask, carry, borrow;
		//BLST_VLA_NEW(limb_t, tmp, n);
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
#endif
		//BLST_VLA_FREE(tmp);
	}

	//#define LSHIFT_MOD_IMPL(bits) \
	//BLST_HD BLST_INLINE void cc_lshift_mod_##bits(vec##bits ret, const vec##bits a, size_t count, \
	//                              const vec##bits p) \
	//{   cc_lshift_mod_n(ret, a, count, p, NLIMBS(bits));   }
	//
	//LSHIFT_MOD_IMPL(256)
	//LSHIFT_MOD_IMPL(384)


	/**
	 * @brief Conditionally negates a large integer modulo a given modulus.
	 *
	 * If `flag` is true, computes `(-a) mod p`; otherwise, returns `a mod p`.
	 * The result is stored in `ret`. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE void cc_cneg_mod_n(limb_t ret[], const limb_t a[], bool_t flag,
		const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t borrow, mask;
		//BLST_VLA_NEW(limb_t, tmp, n);
		limb_t tmp[n];
		size_t i;

		for (borrow = 0, i = 0; i < n; i++) {
			wide_limb_t wide;
			cc_sub_borrow_wide_limb_limb(&wide, p[i], a[i], borrow);
			tmp[i] = wide.lo;
			borrow = wide.hi & 1;
		}

		flag &= vec_is_zero(a, n * sizeof(limb_t)) ^ 1;
		mask = (limb_t)0 - flag;

		for (i = 0; i < n; i++)
			ret[i] = (a[i] & ~mask) | (tmp[i] & mask);
		//BLST_VLA_FREE(tmp);
	}
	//
	//#define CNEG_MOD_IMPL(bits) \
	//BLST_HD BLST_INLINE void cc_cneg_mod_##bits(vec##bits ret, const vec##bits a, bool_t flag, \
	//                            const vec##bits p) \
	//{   cc_cneg_mod_n(ret, a, flag, p, NLIMBS(bits));   }
	//
	//CNEG_MOD_IMPL(256)
	//CNEG_MOD_IMPL(384)


	/**
	 * @brief Checks if a large integer `a` is reduced modulo a given modulus
	 * and not zero.
	 *
	 * Returns 0 if `a` is greater than or equal to `p`, otherwise returns 1.
	 * Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE limb_t cc_check_mod_n(const byte a[], const limb_t p[])
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


	/**
	 * @brief Checks if a large integer `a` is reduced modulo a given modulus
	 *
	 * Returns 0 if `a` is greater than or equal to `p`, otherwise returns 1.
	 * Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE limb_t cc_is_redc_mod_n(const byte a[], const limb_t p[])
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


	//#define CHECK_MOD_IMPL(bits) \
	//BLST_HD BLST_INLINE limb_t cc_check_mod_##bits(const pow##bits a, const vec##bits p) \
	//{   return cc_check_mod_n(a, p, NLIMBS(bits));   }
	//
	//CHECK_MOD_IMPL(256)

	/**
	 * @brief Adds two large integers modulo a given modulus and checks if the result is zero.
	 *
	 * Computes `(a + b) mod p` and stores the result in `ret`. Returns 1 if the result is zero,
	 * otherwise returns 0. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE limb_t cc_add_n_check_mod_n(byte ret[], const byte a[], const byte b[],
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

	//#define ADD_N_CHECK_MOD_IMPL(bits) \
	//BLST_HD BLST_INLINE limb_t cc_add_n_check_mod_##bits(pow##bits ret, const pow##bits a, \
	//                                     const pow##bits b, const vec##bits p) \
	//{   return cc_add_n_check_mod_n(ret, a, b, p, NLIMBS(bits));   }
	//
	//ADD_N_CHECK_MOD_IMPL(256)

	/**
	 * @brief Subtracts two large integers modulo a given modulus and checks if the result is zero.
	 *
	 * Computes `(a - b) mod p` and stores the result in `ret`. Returns 1 if the result is zero,
	 * otherwise returns 0. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE limb_t cc_sub_n_check_mod_n(byte ret[], const byte a[], const byte b[],
		const limb_t p[])
	{
		BLST_ASSUME(n != 0);
		limb_t zero;
		//BLST_VLA_NEW(limb_t, ret_, n);
		//BLST_VLA_NEW(limb_t, a_, n);
		//BLST_VLA_NEW(limb_t, b_, n);
		limb_t ret_[n];
		limb_t a_[n];
		limb_t b_[n];

		limbs_from_le_bytes(a_, a, n * sizeof(limb_t));
		limbs_from_le_bytes(b_, b, n * sizeof(limb_t));

		cc_sub_mod_n<n>(ret_, a_, b_, p);
		zero = vec_is_zero(ret_, n * sizeof(limb_t));

		le_bytes_from_limbs(ret, ret_, n * sizeof(limb_t));

		//BLST_VLA_FREE(ret_);
		//BLST_VLA_FREE(a_);
		//BLST_VLA_FREE(b_);

		return zero ^ 1;
	}
	//
	//#define SUB_N_CHECK_MOD_IMPL(bits) \
	//BLST_HD BLST_INLINE limb_t cc_sub_n_check_mod_##bits(pow##bits ret, const pow##bits a, \
	//                                     const pow##bits b, const vec##bits p) \
	//{   return cc_sub_n_check_mod_n(ret, a, b, p, NLIMBS(bits));   }
	//
	//SUB_N_CHECK_MOD_IMPL(256)


	/**
	 * @brief Converts a large integer from Montgomery form.
	 *
	 * Computes `a * R^(-1) mod p` to convert `a` from Montgomery form to standard form.
	 * The result is stored in `ret`. Assumes `n > 0` and `n` is even.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE void cc_from_mont_n(limb_t ret[], const limb_t a[],
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
	BLST_HD BLST_INLINE void cc_redc_mont_n(limb_t ret[], const limb_t a[],
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


	/**
	 * @brief Performs a right shift modulo operation on a large integer.
	 *
	 * This function computes `(a / 2^count) mod p` iteratively, ensuring the result
	 * remains within the range `[0, p)`. For each bit shifted, the modulus `p` is
	 * conditionally added back if the least significant bit of `a` is set before
	 * the shift. This ensures the result is non-negative and properly reduced.
	 *
	 * @param[out] ret The output array to store the result.
	 * @param[in]  a   The input array representing the large integer.
	 * @param[in]  count The number of bits to shift.
	 * @param[in]  p   The modulus array.
	 * @param[in]  n   The number of limb_ts in `a` and `p`.
	 *
	 * @note Assumes `count > 0`, `n > 0`, and `n` is even.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE void cc_rshift_mod_n(limb_t ret[], const limb_t a[], size_t count,
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
	BLST_HD BLST_INLINE limb_t cc_sgn0_pty_mod_n(const limb_t a[], const limb_t p[])
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


	/*
	 * mul_mont_n without final conditional subtraction, which implies
	 * that modulus is one bit short, which in turn means that there are
	 * no carries to handle between iterations...
	 */
	template<size_t n>
	BLST_HD BLST_INLINE void cc_mul_mont_nonred_n(limb_t ret[], const limb_t a[], const limb_t b[],
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



#if defined(__GNUC__) || defined(__clang__)
# define MSB(x) ({ limb_t ret = (x) >> (BLST_LIMB_SIZE-1); BLST_launder(ret); ret; })
#else
# define MSB(x) ((x) >> (BLST_LIMB_SIZE-1))
#endif

	/**
	 * @brief Computes the number of significant bits in a limb.
	 *
	 * This function calculates the number of bits required to represent the given limb `l`.
	 * It efficiently determines the position of the most significant bit set to `1`.
	 *
	 * @param[in] l The input limb whose significant bits are to be counted.
	 *
	 * @return The number of significant bits in `l`. Returns `0` if `l` is `0`.
	 *
	 * @note This function assumes that `limb_t` is either 32-bit or 64-bit.
	 */
	BLST_HD BLST_INLINE size_t cc_num_bits(limb_t l)
	{

#if (defined(__GNUC__) || defined(__clang__)) && !defined(__CUDA_ARCH__) 
		// Use GCC/Clang built-in functions for counting leading zeros
#if BLST_LIMB_SIZE == 64
		return (l != 0) * (64 - __builtin_clzll(l ? l : 1));
#else
		return (l != 0) * (32 - __builtin_clz(l ? l : 1));
#endif

#elif defined(_MSC_VER) && !defined(__CUDA_ARCH__) 
		// Use MSVC intrinsics for counting leading zeros
		unsigned long index;
#if BLST_LIMB_SIZE == 64
		_BitScanReverse64(&index, l);
#else
		_BitScanReverse(&index, l);
#endif
		return (l != 0) * (index + 1);

#else
		// Fallback to portable implementation
		limb_t x, mask;
		size_t bits = is_zero(l) ^ 1;

		if (sizeof(limb_t) == 8) {
			x = l >> (32 & (8 * sizeof(limb_t) - 1));
			mask = 0 - MSB(0 - x);
			bits += 32 & mask;
			l ^= (x ^ l) & mask;
		}

		x = l >> 16;
		mask = 0 - MSB(0 - x);
		bits += 16 & mask;
		l ^= (x ^ l) & mask;

		x = l >> 8;
		mask = 0 - MSB(0 - x);
		bits += 8 & mask;
		l ^= (x ^ l) & mask;

		x = l >> 4;
		mask = 0 - MSB(0 - x);
		bits += 4 & mask;
		l ^= (x ^ l) & mask;

		x = l >> 2;
		mask = 0 - MSB(0 - x);
		bits += 2 & mask;
		l ^= (x ^ l) & mask;

		bits += l >> 1;

		return bits;
#endif
	}

#if defined(__clang_major__) && __clang_major__>7
	__attribute__((optnone))
#endif
		BLST_HD BLST_INLINE  limb_t cc_lshift_2(limb_t hi, limb_t lo, size_t l)
	{
		size_t r = BLST_LIMB_SIZE - l;
		limb_t mask = 0 - (is_zero(l) ^ 1);
		return (hi << (l & (BLST_LIMB_SIZE - 1))) | ((lo & mask) >> (r & (BLST_LIMB_SIZE - 1)));
	}

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
	BLST_HD BLST_INLINE void cc_ab_approximation_n(limb_t a_[2], const limb_t a[],
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

	typedef struct { limb_t f0, g0, f1, g1; } factors;

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
	BLST_HD BLST_INLINE void cc_inner_loop_n(factors* fg, const limb_t a_[2], const limb_t b_[2])
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
	 * @return The sign of the result (0 for non-negative, 1 for negative).
	 *
	 * @note Assumes `n > 0`.
	 */
	template<rsize_t n>
	BLST_HD BLST_INLINE limb_t cc_cneg_n(limb_t ret[], const limb_t a[], limb_t neg)
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

		return 0 - MSB(wide.lo);
	}


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
	BLST_HD BLST_INLINE limb_t cc_add_n(limb_t ret[], const limb_t a[], limb_t b[])
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
	BLST_HD BLST_INLINE limb_t cc_umul_n(limb_t ret[], const limb_t a[], limb_t b)
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
	BLST_HD BLST_INLINE limb_t cc_smul_n_shift_n(limb_t ret[], const limb_t a[], limb_t* f_,
		const limb_t b[], limb_t* g_)
	{
		BLST_ASSUME(n != 0);
		limb_t f, g, neg, carry, hi;
		limb_t a_[n + 1];
		limb_t b_[n + 1];
		size_t i;

		/* |a|*|f_| */
		f = *f_;
		neg = 0 - MSB(f);
		f = (f ^ neg) - neg;            /* ensure |f| is positive */
		(void)cc_cneg_n<n>(a_, a, neg);
		hi = cc_umul_n<n>(a_, a_, f);
		a_[n] = hi - (f & neg);

		/* |b|*|g_| */
		g = *g_;
		neg = 0 - MSB(g);
		g = (g ^ neg) - neg;            /* ensure |g| is positive */
		(void)cc_cneg_n<n>(b_, b, neg);
		hi = cc_umul_n<n>(b_, b_, g);
		b_[n] = hi - (g & neg);

		/* |a|*|f_| + |b|*|g_| */
		(void)cc_add_n<n+1>(a_, a_, b_);

		/* (|a|*|f_| + |b|*|g_|) >> k */
		for (carry = a_[0], i = 0; i < n; i++) {
			hi = carry >> (BLST_LIMB_SIZE - 2);
			carry = a_[i + 1];
			ret[i] = hi | (carry << 2);
		}

		/* ensure result is non-negative, fix up |f_| and |g_| accordingly */
		neg = 0 - MSB(carry);
		*f_ = (*f_ ^ neg) - neg;
		*g_ = (*g_ ^ neg) - neg;
		(void)cc_cneg_n<n>(ret, ret, neg);

		return neg;
	}

	/* @brief Computes |u|*|f_| + |v|*|g_|
	 *
	 * Assumes n>0.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE limb_t cc_smul_2n(limb_t ret[], const limb_t u[], limb_t f,
		const limb_t v[], limb_t g)
	{
		BLST_ASSUME(n != 0);
		limb_t neg, hi;

		limb_t u_[n];
		limb_t v_[n];

		/* |u|*|f_| */
		neg = 0 - MSB(f);
		f = (f ^ neg) - neg;            /* ensure |f| is positive */
		neg = cc_cneg_n<n>(u_, u, neg);
		hi = cc_umul_n<n>(u_, u_, f) - (f & neg);

		/* |v|*|g_| */
		neg = 0 - MSB(g);
		g = (g ^ neg) - neg;            /* ensure |g| is positive */
		neg = cc_cneg_n<n>(v_, v, neg);
		hi += cc_umul_n<n>(v_, v_, g) - (g & neg);

		/* |u|*|f_| + |v|*|g_| */
		hi += cc_add_n<n>(ret, u_, v_);

		return hi;
	}


	/**
	 * @brief Computes the modular inverse of a large integer.
	 *
	 * Computes `ret = a^(-1) mod p` using the constant-time extended Euclidean algorithm.
	 * Assumes `n > 0` and `n` is even. modx should be a left aligned copy of mod.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE void cc_ct_inverse_mod_n(limb_t ret[], const limb_t inp[],
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

		sign = 0 - MSB(top);    /* top is 1, 0 or -1 */
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
		(void)cc_cneg_n<n>(a, a, 0 - MSB(sign));
		cc_add_n<n>(ret + n, ret + n, a);
	}

	/*
	 * Copy of inner_loop_n above, but with |L| updates.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE limb_t cc_legendre_loop_n(limb_t L, factors* fg, const limb_t a_[2],
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

	/**
	 * @brief Checks if a large integer is a quadratic residue modulo a given modulus.
	 *
	 * Returns 1 if `a` is a quadratic residue modulo `p`, otherwise returns 0.
	 * Assumes `n > 0` and `n` is even.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE bool_t cc_ct_is_sqr_mod_n(const limb_t inp[], const limb_t mod[])
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
	 * @note Assumes `count > 0` and `a`, `b`, and `p` are valid 384-bit integers.
	 */
	BLST_HD BLST_INLINE void cc_sqr_n_mul_mont_384(
		limb_t* ret, const limb_t* a, size_t count,
		const limb_t* p, limb_t n0, const limb_t* b)
	{
		BLST_ASSUME(count != 0);
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
	 * Computes the quotient of a 3-limb dividend divided by a 2-limb divisor.
	 * Assumes the divisor is left-aligned and the dividend is not larger than the divisor.
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
	BLST_HD BLST_INLINE limb_t cc_div_3_limbs(const limb_t div_top[2], limb_t d_lo, limb_t d_hi)
	{
		wide_limb_t Rx;
		//llimb_t Rx;
		limb_t r_lo = div_top[0], r_hi = div_top[1];
		limb_t Q = 0, mask, borrow, rx;
		size_t i;

		for (i = 0; i < BLST_LIMB_SIZE; i++) {
			/* "borrow, Rx = R - D" */
			cc_sub_borrow_wide_limb_limb(&Rx, r_lo, d_lo, 0);
			rx = Rx.lo;
			borrow = Rx.hi & 1;

			cc_sub_borrow_wide_limb_limb(&Rx, r_hi, d_hi + (limb_t)borrow, 0);
			borrow = Rx.hi;// not & 1. want 0000 or 1111


			/* "if (R >= D) R -= D" */
			r_lo = ((r_lo ^ rx) & borrow) ^ rx;
			//rx = (limb)Rx;
			rx = Rx.lo;
			r_hi = ((r_hi ^ rx) & borrow) ^ rx;

			Q <<= 1;
			Q |= ~borrow & 1;

			/* "D >>= 1" */
			d_lo >>= 1; d_lo |= d_hi << (BLST_LIMB_SIZE - 1);
			d_hi >>= 1;
		}

		mask = 0 - MSB(Q);  /* does it overflow? */

		/* "borrow, Rx = R - D" */
		cc_sub_borrow_wide_limb_limb(&Rx, r_lo, d_lo, 0);
		rx = Rx.lo;
		borrow = Rx.hi & 1;

		cc_sub_borrow_wide_limb_limb(&Rx, r_hi, d_hi + (limb_t)borrow, 0);
		borrow = Rx.hi & 1;

		Q <<= 1;
		Q |= borrow ^ 1;

		return (Q | mask);
	}

	/**
	 * @brief Computes the quotient and remainder of a large integer division.
	 *
	 * Computes `quotient = dividend / divisor` and updates `div_rem` with the remainder.
	 * Assumes `n > 0` and `n` is even.
	 */
	template<size_t n>
	BLST_HD BLST_INLINE limb_t cc_quot_rem_n(limb_t* div_rem, const limb_t* divisor,
		limb_t quotient)
	{
		wide_limb_t wide;
		limb_t carry, mask, borrow;
		//BLST_VLA_NEW(limb_t, tmp, n + 1);
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
			cc_sub_borrow_wide_limb_limb(&wide, div_rem[i], tmp[i] + borrow, 0);
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

	/**
	 * @brief Computes the quotient and remainder of a large integer division.
	 *
	 * Computes `quotient = dividend / divisor` and updates `div_rem` with the remainder.
	 */
	BLST_HD BLST_INLINE limb_t cc_quot_rem_128(limb_t* div_rem, const limb_t* divisor,
		limb_t quotient)
	{
		return cc_quot_rem_n<NLIMBS(128)>(div_rem, divisor, quotient);
	}


	/**
	 * @brief Computes the quotient and remainder of a large integer division.
	 *
	 * Computes `quotient = dividend / divisor` and updates `div_rem` with the remainder.
	 */
	BLST_HD BLST_INLINE limb_t cc_quot_rem_64(limb_t* div_rem, const limb_t* divisor,
		limb_t quotient)
	{
		return cc_quot_rem_n<NLIMBS(64)>(div_rem, divisor, quotient);
	}


	/*
	 * Prefetching
	 */
	BLST_HD BLST_INLINE void cc_vec_prefetch(const void* ptr, size_t len)
	{
		// prefetch if we have the intrinsic
#if defined(__GNUC__) || defined(__clang__)
		__builtin_prefetch(ptr, 0, 3);
#elif __CUDA_ARCH__
		// todo
		(void)ptr; (void)len;
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
		_mm_prefetch((const char*)ptr, _MM_HINT_T0);
		(void)ptr; (void)len;
#else
		(void)ptr; (void)len;
#endif
	}

	/*
	 * SHA-256
	 */
#define ROTR(x,n)	((x)>>n | (x)<<(32-n))
#define Sigma0(x)	(ROTR((x),2) ^ ROTR((x),13) ^ ROTR((x),22))
#define Sigma1(x)	(ROTR((x),6) ^ ROTR((x),11) ^ ROTR((x),25))
#define sigma0(x)	(ROTR((x),7) ^ ROTR((x),18) ^ ((x)>>3))
#define sigma1(x)	(ROTR((x),17) ^ ROTR((x),19) ^ ((x)>>10))
#define Ch(x,y,z)	(((x) & (y)) ^ ((~(x)) & (z)))
#define Maj(x,y,z)	(((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))

	 /* @brief SHA-256 block compression function.
	  *
	  * Processes `blocks` number of 64-byte blocks from `inp`, updating the state `v`.
	  * The implementation is based on the public domain code from
	  *
	 */
	BLST_HD BLST_INLINE void cc_blst_sha256_block_data_order(unsigned int* v, const void* inp,
		size_t blocks)
	{
		const unsigned int K256[64] = {
			0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
			0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
			0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
			0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
			0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
			0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
			0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
			0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
			0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
			0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
			0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
			0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
			0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
			0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
			0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
			0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
		};
		unsigned int X[16], l, a, b, c, d, e, f, g, h, s0, s1, T1, T2;
		const unsigned char* data = (const unsigned char*)inp;
		size_t round;

		a = v[0];
		b = v[1];
		c = v[2];
		d = v[3];
		e = v[4];
		f = v[5];
		g = v[6];
		h = v[7];

		while (blocks--) {
			for (round = 0; round < 16; round++) {
				l = (unsigned int)data[0] << 24;
				l |= (unsigned int)data[1] << 16;
				l |= (unsigned int)data[2] << 8;
				l |= (unsigned int)data[3];
				data += 4;
				T1 = X[round] = l;
				T1 += h + Sigma1(e) + Ch(e, f, g) + K256[round];
				T2 = Sigma0(a) + Maj(a, b, c);
				h = g;
				g = f;
				f = e;
				e = d + T1;
				d = c;
				c = b;
				b = a;
				a = T1 + T2;
			}

			for (; round < 64; round++) {
				s0 = X[(round + 1) & 0x0f];
				s0 = sigma0(s0);
				s1 = X[(round + 14) & 0x0f];
				s1 = sigma1(s1);

				T1 = X[round & 0xf] += s0 + s1 + X[(round + 9) & 0xf];
				T1 += h + Sigma1(e) + Ch(e, f, g) + K256[round];
				T2 = Sigma0(a) + Maj(a, b, c);
				h = g;
				g = f;
				f = e;
				e = d + T1;
				d = c;
				c = b;
				b = a;
				a = T1 + T2;
			}

			a += v[0]; v[0] = a;
			b += v[1]; v[1] = b;
			c += v[2]; v[2] = c;
			d += v[3]; v[3] = d;
			e += v[4]; v[4] = e;
			f += v[5]; v[5] = f;
			g += v[6]; v[6] = g;
			h += v[7]; v[7] = h;
		}
	}
#undef ROTR
#undef Sigma0
#undef Sigma1
#undef sigma0
#undef sigma1
#undef Ch
#undef Maj

	// copies a 64 byte length array into dst_. i.e memcpy(dst, src, 64)
	BLST_HD BLST_INLINE void cc_blst_sha256_hcopy(unsigned int dst[8], const unsigned int src[8])
	{
		size_t i;

		for (i = 0; i < 8; i++)
			dst[i] = src[i];
	}


	// performs the final compression function for sha256
	BLST_HD BLST_INLINE void cc_blst_sha256_emit(unsigned char md[32], const unsigned int h[8])
	{
		size_t i;

		for (i = 0; i < 8; i++, md += 4) {
			unsigned int h_i = h[i];
			md[0] = (unsigned char)(h_i >> 24);
			md[1] = (unsigned char)(h_i >> 16);
			md[2] = (unsigned char)(h_i >> 8);
			md[3] = (unsigned char)h_i;
		}
	}

	// copies a variable length array into dst_. i.e memcpy(dst, src, len)
	BLST_HD BLST_INLINE void cc_blst_sha256_bcopy(void* dst_, const void* src_, size_t len)
	{
		unsigned char* dst = (unsigned char*)dst_;
		const unsigned char* src = (const unsigned char*)src_;
		size_t i;

		for (i = 0; i < len; i++)
			dst[i] = src[i];
	}


#undef BLST_VLA_NEW
#undef BLST_VLA_FREE
	//BLST_EXTERN_C_END
}

#endif