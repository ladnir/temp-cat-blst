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

namespace blst
{



	//BLST_EXTERN_C_BEGIN
	typedef struct { limb_t lo, hi; } wide_limb_t;


	// generic
	namespace portable
	{
		BLST_HD BLST_FORCEINLINE void cc_mul_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b)
		{
			int half = sizeof(limb_t) * 4;
			limb_t low_mask = ((limb_t)1 << half) - 1;
			/* Portable C implementation for 64-bit limbs*/
			/* Break operands into 32-bit halves        */
			limb_t a_lo = (a)&low_mask;
			limb_t a_hi = (a) >> half;
			limb_t b_lo = (b)&low_mask;
			limb_t b_hi = (b) >> half;

			/* Four partial products: */
			limb_t p0 = a_lo * b_lo;
			limb_t p1 = a_lo * b_hi;
			limb_t p2 = a_hi * b_lo;
			limb_t p3 = a_hi * b_hi;

			/* Combine middle products with carry propagation */
			limb_t middle = (p0 >> half) + (p1 & low_mask) + (p2 & low_mask);
			limb_t carry = middle >> half;

			/* Assemble final result    */
			(ret)->lo = (p0 & low_mask) + (middle << half);
			(ret)->hi = p3 + (p1 >> half) + (p2 >> half) + carry;
		}


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
	}


#if BLST_LIMB_SIZE == 64 && defined(BLST_ADX) && !defined(__CUDA_ARCH__)

	// 64 bit with wide operations with adx support
	static_assert(sizeof(unsigned long long*) == sizeof(limb_t), "limb_t must be 64 bits");

	namespace adx64
	{
		BLST_FORCEINLINE void cc_mul_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b) {
			ret->lo = _mulx_u64(a, b, (unsigned long long*) & ret->hi);
		}

		BLST_FORCEINLINE void cc_add_wide_wide_limb(wide_limb_t* ret, const wide_limb_t& a, const limb_t& b) {
			ret->hi = a.hi + _addcarry_u64(0, a.lo, b, (unsigned long long*) &ret->lo);
		}

		BLST_FORCEINLINE void cc_add_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b) {
			ret->hi = _addcarry_u64(0, a, b, (unsigned long long*) &ret->lo);
		}

		BLST_FORCEINLINE void cc_add_carry_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b, const limb_t& carry) {
			BLST_ASSUME(carry <= 1);
			ret->hi = _addcarry_u64((unsigned char)carry, a, b, (unsigned long long*) &ret->lo);
		}

		BLST_FORCEINLINE void cc_sub_borrow_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b, const limb_t& borrow)
		{
			BLST_ASSUME(borrow <= 1);
			ret->hi = 0 - (limb_t)_subborrow_u64((unsigned char)borrow, a, b, (unsigned long long*) &ret->lo);
		}

	}
#endif

#if BLST_LIMB_SIZE == 32 && (defined(_M_X64) || defined(__x86_64__))
	// 32 bit with 64 bit support
	namespace limb32on64
	{
		BLST_FORCEINLINE void cc_mul_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b) {
			*(unsigned long long*)(ret) = (unsigned long long)(a) * (unsigned long long)(b);
		}
		BLST_FORCEINLINE void cc_add_wide_wide_limb(wide_limb_t* ret, const wide_limb_t& a, const limb_t& b) {
			*(unsigned long long*)(ret) = *(unsigned long long*) & (a)+(unsigned long long)(b);
		}
		BLST_FORCEINLINE void cc_add_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b) {
			*(unsigned long long*)(ret) = (a)+(unsigned long long)(b);
		}
		BLST_FORCEINLINE void cc_add_carry_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b, const limb_t& carry) {
			BLST_ASSUME(carry <= 1);
			*(unsigned long long*)(ret) = (a)+(unsigned long long)(b) + carry;
		}
		BLST_FORCEINLINE void cc_sub_borrow_wide_limb_limb(wide_limb_t* ret, const limb_t& a, const limb_t& b, const limb_t& borrow)
		{
			BLST_ASSUME(borrow <= 1);
			*(unsigned long long*)(ret) = (a)-(unsigned long long)(b) - borrow;
		}
	}
	static_assert(&((wide_limb_t*)0)->lo < &((wide_limb_t*)0)->hi,
		"this implementation assumes we can cast to 64 bit");
	//#define cc_mul_wide_limb_limb(ret, a, b)\
	//    do{*(unsigned long long*)(ret) = (unsigned long long)(a) * (unsigned long long)(b);}while(0)
	//
	//#define cc_add_wide_wide_limb(ret, a, b)\
	//    do{*(unsigned long long*)(ret) = *(unsigned long long*)&(a) + (unsigned long long)(b);}while(0)
	//
	//#define cc_add_wide_limb_limb(ret, a, b)\
	//    do{*(unsigned long long*)(ret) = (a) + (unsigned long long)(b);}while(0)
	//
	//#define cc_sub_wide_limb_wide(ret, a, b)\
	//    do{*(unsigned long long*)(ret) = (a) - *(unsigned long long*)(&(b));}while(0)
	//
	//#define cc_add_carry_wide_limb_limb(ret, a, b, carry) \
	//    do{ assert(carry <= 1); BLST_ASSUME(carry <=1); \
	//    *(unsigned long long*)(ret) = (a) + (unsigned long long)(b) + carry; } while(0)
	//
	//#define cc_sub_borrow_wide_limb_limb(ret, a, b, borrow) \
	//    do{ assert(borrow <= 1); BLST_ASSUME(borrow <=1); \
	//    *(unsigned long long*)(ret) = (a) - (unsigned long long)(b) - borrow; } while(0)
#endif

#if defined(__CUDA_ARCH__)
	static_assert(std::is_same<limb_t, uint64_t>::value,
		"this implementation assumes limb_t is uint64_t");

	namespace cuda64
	{

		// 64x64 -> 128
		BLST_DEVICE BLST_FORCEINLINE
			void cc_mul_wide_limb_limb(
				wide_limb_t* ret,
				const limb_t& a,
				const limb_t& b
			)
		{
			uint64_t lo = a * b;
			uint64_t hi = __umul64hi(a, b);
#ifndef NDEBUG
			// verify against portable
			wide_limb_t portable_res;
			portable::cc_mul_wide_limb_limb(&portable_res, a, b);
			BLST_ASSERT(portable_res.lo == lo);
			BLST_ASSERT(portable_res.hi == hi);
#endif
			ret->lo = lo;
			ret->hi = hi;
		}

		// (a.lo, a.hi) + b -> 128
		BLST_DEVICE BLST_FORCEINLINE
			void cc_add_wide_wide_limb(
				wide_limb_t* ret,
				const wide_limb_t& a,
				const limb_t& b
			)
		{
			uint64_t lo = a.lo;
			uint64_t hi = a.hi;
			asm volatile(
				"add.cc.u64 %0, %0, %2;\n\t"
				"addc.u64   %1, %1, 0;\n\t"
				: "+l"(lo), "+l"(hi)
				: "l"(b)
				);
#ifndef NDEBUG
			// verify against portable
			wide_limb_t portable_res;
			portable::cc_add_wide_wide_limb(&portable_res, a, b);
			BLST_ASSERT(portable_res.lo == lo);
			BLST_ASSERT(portable_res.hi == hi);
#endif // !NDEBUG

			ret->lo = lo;
			ret->hi = hi;
		}

		// a + b -> (lo, carry)
		BLST_DEVICE BLST_FORCEINLINE
			void cc_add_wide_limb_limb(
				wide_limb_t* ret,
				const limb_t& a,
				const limb_t& b
			)
		{
			uint64_t out_lo;
			uint64_t out_hi;
			asm volatile(
				"add.cc.u64 %0, %2, %3;\n\t"
				"addc.u64   %1, 0, 0;\n\t"   // out_hi = carry (0/1)
				: "=&l"(out_lo), "=&l"(out_hi)
				: "l"(a), "l"(b)
				);
#ifndef NDEBUG
			// verify against portable
			wide_limb_t portable_res;
			portable::cc_add_wide_limb_limb(&portable_res, a, b);
			BLST_ASSERT(portable_res.lo == out_lo);
			BLST_ASSERT(portable_res.hi == out_hi);
#endif // !NDEBUG

			ret->lo = out_lo;
			ret->hi = out_hi;

		}

		// (a.lo, a.hi) - b -> 128
		BLST_DEVICE BLST_FORCEINLINE
			void cc_sub_wide_wide_limb(
				wide_limb_t* ret,
				const wide_limb_t& a,
				const limb_t& b
			)
		{
			uint64_t lo = a.lo;
			uint64_t hi = a.hi;
			asm volatile(
				"sub.cc.u64 %0, %0, %2;\n\t"
				"subc.u64   %1, %1, 0;\n\t"
				: "+l"(lo), "+l"(hi)
				: "l"(b)
				);
#ifndef NDEBUG
			// verify against portable
			wide_limb_t portable_res;
			portable::cc_sub_wide_wide_limb(&portable_res, a, b);
			BLST_ASSERT(portable_res.lo == lo);
			BLST_ASSERT(portable_res.hi == hi);
#endif
			ret->lo = lo;
			ret->hi = hi;
		}

		// a - b -> (lo, hi = -(borrow))
		BLST_DEVICE BLST_FORCEINLINE
			void cc_sub_wide_limb_limb(
				wide_limb_t* ret,
				const limb_t& a,
				const limb_t& b
			)
		{
			uint64_t out_lo, out_hi;
			asm volatile(
				"sub.cc.u64 %0, %2, %3;\n\t"  // out_lo = a - b; CC = borrow
				"subc.u64   %1, 0, 0;\n\t"    // out_hi = 0 - borrow  => 0 or ~0
				: "=&l"(out_lo), "=&l"(out_hi)
				: "l"(a), "l"(b)
				);
#ifndef NDEBUG
			// verify against portable
			wide_limb_t portable_res;
			portable::cc_sub_wide_limb_limb(&portable_res, a, b);
			BLST_ASSERT(portable_res.lo == out_lo);
			BLST_ASSERT(portable_res.hi == out_hi);
#endif
			ret->lo = out_lo;
			ret->hi = out_hi;
		}

		//// a - b -> (lo, borrow)
		//BLST_DEVICE BLST_FORCEINLINE
		//	void cc_sub_wide_limb_limb(
		//		wide_limb_t* ret,
		//		const limb_t& a,
		//		const limb_t& b
		//	)
		//{
		//	uint64_t out_lo;
		//	uint64_t out_hi;
		//	asm volatile(
		//		"sub.cc.u64 %0, %2, %3;\n\t"
		//		"addc.u64   %1, 0, 0;\n\t"   // out_hi = borrow (0/1)
		//		: "=&l"(out_lo), "=&l"(out_hi)
		//		: "l"(a), "l"(b)
		//		);
		//	ret->lo = out_lo;
		//	ret->hi = out_hi;
		//}

		// a + b + c -> (lo, carry)
		BLST_DEVICE BLST_FORCEINLINE
			void cc_add_carry_wide_limb_limb(
				wide_limb_t* ret,
				const limb_t& a,
				const limb_t& b,
				const limb_t& c   // 0 or 1
			)
		{
			BLST_ASSERT(c <= 1);
			uint64_t out_lo;
			uint64_t out_hi;
			asm volatile(
				"add.cc.u64  %0, %2, %3;\n\t" // lo = a + b ; CC = carry_ab
				"addc.u64    %1, 0, 0;\n\t"   // hi = carry_ab
				"add.cc.u64  %0, %0, %4;\n\t" // lo += c ; CC = carry_c
				"addc.u64    %1, %1, 0;\n\t"  // hi += carry_c
				: "=&l"(out_lo),
				"=&l"(out_hi)
				: "l"(a),
				"l"(b),
				"l"(c)
				);
#ifndef NDEBUG
			// verify against portable
			wide_limb_t portable_res;
			portable::cc_add_carry_wide_limb_limb(&portable_res, a, b, c);
			
			if (portable_res.lo != out_lo)
			{
				printf("cc_add_carry_wide_limb_limb mismatch lo: portable %llx vs cuda %llx\n a=%llx b=%llx c=%llx\n",
					(unsigned long long)portable_res.lo,
					(unsigned long long)out_lo,
					(unsigned long long)a,
					(unsigned long long)b,
					(unsigned long long)c);
			}
			BLST_ASSERT(portable_res.lo == out_lo);
			BLST_ASSERT(portable_res.hi == out_hi);
#endif
			ret->lo = out_lo;
			ret->hi = out_hi;
		}


		// a - b - c  (match portable: do it in two steps)
		BLST_DEVICE BLST_FORCEINLINE
			void cc_sub_borrow_wide_limb_limb(
				wide_limb_t* ret,
				const limb_t& a,
				const limb_t& b,
				const limb_t& c
			)
		{
			// replicate the portable two-call sequence exactly
			cc_sub_wide_limb_limb(ret, a, b);     // hi = -(borrow_ab)
			cc_sub_wide_wide_limb(ret, *ret, c);  // subtract c from the 128-bit (hi,lo)
		}
	}
#endif



#ifdef __CUDA_ARCH__
	using namespace cuda64;
#elif BLST_LIMB_SIZE == 64 && defined(BLST_ADX) 
	using namespace adx64;
#elif BLST_LIMB_SIZE == 32 && (defined(_M_X64) || defined(__x86_64__))
	using namespace limb32on64;
#else
	using namespace portable;
#endif


	// ret = a * b + c
#define cc_mul_add_wide_limb_limb(ret, a, b, c)\
    do{ cc_mul_wide_limb_limb(ret, a, b);\
    cc_add_wide_wide_limb(ret, *ret, c); } while(0)

// ret = a * b + c + d
#define cc_mul_add_add_wide_limb_limb(ret, a, b, c, d)\
    do {cc_mul_wide_limb_limb(ret, a, b);\
    cc_add_wide_wide_limb(ret, *ret, c);\
    cc_add_wide_wide_limb(ret, *ret, d); } while(0)

// ret = a + b + c
#define cc_add_add_wide_limb_limb(ret, a, b, c)\
    do { cc_add_wide_limb_limb(ret, a, b);\
    cc_add_wide_wide_limb(ret, *ret, c) ;} while(0)








	// Basic limb operations


	/** * @brief Multiplies two large integers using Montgomery multiplication.
	 *
	 * Computes `(a * b * R^-1) mod p`, where `R` is the Montgomery radix, and stores
	 * the result in `ret`. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD void cc_mul_mont_n(limb_t* ret, const limb_t* a, const limb_t* b,
		const limb_t* p, limb_t n0);

	/** * @brief Adds two large integers modulo a given modulus.
	 *
	 * Computes `(a + b) mod p` and stores the result in `ret`. Handles carry propagation
	 * and ensures the result is reduced modulo `p`. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD void cc_add_mod_n(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[]);

	/**
	 * @brief Subtracts two large integers modulo a given modulus.
	 *
	 * Computes `(a - b) mod p` and stores the result in `ret`. Handles borrow propagation
	 * and ensures the result is reduced modulo `p`. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD  void cc_sub_mod_n(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[]);

	/**
	 * @brief Multiplies a large integer by 3 modulo a given modulus.
	 *
	 * Computes `(3 * a) mod p` and stores the result in `ret`. Handles carry propagation
	 * and ensures the result is reduced modulo `p`. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD  void cc_mul_by_3_mod_n(limb_t ret[], const limb_t a[], const limb_t p[]);

	/**
	 * @brief Performs a left shift modulo operation on a large integer.
	 *
	 * Computes `(a << count) mod p` iteratively, ensuring the result remains within
	 * the range `[0, p)`. Assumes `count > 0` and `n > 0`.
	 */
	template<size_t n>
	BLST_HD  void cc_lshift_mod_n(limb_t ret[], const limb_t a[], size_t count,
		const limb_t p[]);

	/**
	 * @brief Conditionally negates a large integer modulo a given modulus.
	 *
	 * If `flag` is true, computes `(-a) mod p`; otherwise, returns `a mod p`.
	 * The result is stored in `ret`. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD  void cc_cneg_mod_n(limb_t ret[], const limb_t a[], bool_t flag,
		const limb_t p[]);


	/**
	 * @brief Checks if a large integer `a` is reduced modulo a given modulus
	 * and not zero.
	 *
	 * Returns 0 if `a` is greater than or equal to `p`, otherwise returns 1.
	 * Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD  limb_t cc_check_mod_n(const byte a[], const limb_t p[]);


	/**
	 * @brief Checks if a large integer `a` is reduced modulo a given modulus
	 *
	 * Returns 0 if `a` is greater than or equal to `p`, otherwise returns 1.
	 * Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD  limb_t cc_is_redc_mod_n(const byte a[], const limb_t p[]);

	/**
	 * @brief Adds two large integers modulo a given modulus and checks if the result is zero.
	 *
	 * Computes `(a + b) mod p` and stores the result in `ret`. Returns 1 if the result is zero,
	 * otherwise returns 0. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD  limb_t cc_add_n_check_mod_n(byte ret[], const byte a[], const byte b[],
		const limb_t p[]);


	/**
	 * @brief Subtracts two large integers modulo a given modulus and checks if the result is zero.
	 *
	 * Computes `(a - b) mod p` and stores the result in `ret`. Returns 1 if the result is zero,
	 * otherwise returns 0. Assumes `n > 0`.
	 */
	template<size_t n>
	BLST_HD  limb_t cc_sub_n_check_mod_n(byte ret[], const byte a[], const byte b[],
		const limb_t p[]);

	/**
	 * @brief Converts a large integer from Montgomery form.
	 *
	 * Computes `a * R^(-1) mod p` to convert `a` from Montgomery form to standard form.
	 * The result is stored in `ret`. Assumes `n > 0` and `n` is even.
	 */
	template<size_t n>
	BLST_HD  void cc_from_mont_n(limb_t ret[], const limb_t a[],
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
	BLST_HD  void cc_redc_mont_n(limb_t ret[], const limb_t a[],
		const limb_t p[], limb_t n0);

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
	BLST_HD  void cc_rshift_mod_n(limb_t ret[], const limb_t a[], size_t count,
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
	BLST_HD  limb_t cc_sgn0_pty_mod_n(const limb_t a[], const limb_t p[]);


	/*
	 * mul_mont_n without final conditional subtraction, which implies
	 * that modulus is one bit short, which in turn means that there are
	 * no carries to handle between iterations...
	 */
	template<size_t n>
	BLST_HD  void cc_mul_mont_nonred_n(limb_t ret[], const limb_t a[], const limb_t b[],
		const limb_t p[], limb_t n0);


#if defined(__GNUC__) || defined(__clang__)
# define BLST_MSB(x) ({ limb_t ret = (x) >> (BLST_LIMB_SIZE-1); BLST_launder(ret); ret; })
#else
# define BLST_MSB(x) ((x) >> (BLST_LIMB_SIZE-1))
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
	BLST_HD BLST_FORCEINLINE size_t cc_num_bits(limb_t l)
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
			mask = 0 - BLST_MSB(0 - x);
			bits += 32 & mask;
			l ^= (x ^ l) & mask;
		}

		x = l >> 16;
		mask = 0 - BLST_MSB(0 - x);
		bits += 16 & mask;
		l ^= (x ^ l) & mask;

		x = l >> 8;
		mask = 0 - BLST_MSB(0 - x);
		bits += 8 & mask;
		l ^= (x ^ l) & mask;

		x = l >> 4;
		mask = 0 - BLST_MSB(0 - x);
		bits += 4 & mask;
		l ^= (x ^ l) & mask;

		x = l >> 2;
		mask = 0 - BLST_MSB(0 - x);
		bits += 2 & mask;
		l ^= (x ^ l) & mask;

		bits += l >> 1;

		return bits;
#endif
	}


#if defined(__clang_major__) && __clang_major__>7
	__attribute__((optnone))
#endif
		BLST_HD BLST_FORCEINLINE limb_t cc_lshift_2(limb_t hi, limb_t lo, size_t l)
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
	BLST_HD  void cc_ab_approximation_n(limb_t a_[2], const limb_t a[],
		limb_t b_[2], const limb_t b[]);

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
	BLST_HD  void cc_inner_loop_n(factors* fg, const limb_t a_[2], const limb_t b_[2]);


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
	template<size_t n>
	BLST_HD  limb_t cc_cneg_n(limb_t ret[], const limb_t a[], limb_t neg);


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
	BLST_HD  limb_t cc_add_n(limb_t ret[], const limb_t a[], limb_t b[]);

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
	BLST_HD  limb_t cc_umul_n(limb_t ret[], const limb_t a[], limb_t b);

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
	BLST_HD  limb_t cc_smul_n_shift_n(limb_t ret[], const limb_t a[], limb_t* f_,
		const limb_t b[], limb_t* g_);

	/* @brief Computes |u|*|f_| + |v|*|g_|
	 *
	 * Assumes n>0.
	 */
	template<size_t n>
	BLST_HD  limb_t cc_smul_2n(limb_t ret[], const limb_t u[], limb_t f,
		const limb_t v[], limb_t g);

	/**
	 * @brief Computes the modular inverse of a large integer.
	 *
	 * Computes `ret = a^(-1) mod p` using the constant-time extended Euclidean algorithm.
	 * Assumes `n > 0` and `n` is even. modx should be a left aligned copy of mod.
	 */
	template<size_t n>
	BLST_HD  void cc_ct_inverse_mod_n(limb_t ret[], const limb_t inp[],
		const limb_t mod[], const limb_t modx[]);

	/*
	 * Copy of inner_loop_n above, but with |L| updates.
	 */
	template<size_t n>
	BLST_HD  limb_t cc_legendre_loop_n(limb_t L, factors* fg, const limb_t a_[2],
		const limb_t b_[2]);

	/**
	 * @brief Checks if a large integer is a quadratic residue modulo a given modulus.
	 *
	 * Returns 1 if `a` is a quadratic residue modulo `p`, otherwise returns 0.
	 * Assumes `n > 0` and `n` is even.
	 */
	template<size_t n>
	BLST_HD  bool_t cc_ct_is_sqr_mod_n(const limb_t inp[], const limb_t mod[]);


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
	BLST_HD  void cc_sqr_n_mul_mont_384(
		limb_t* ret, const limb_t* a, size_t count,
		const limb_t* p, limb_t n0, const limb_t* b);

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
	BLST_HD  limb_t cc_div_3_limbs(const limb_t div_top[2], limb_t d_lo, limb_t d_hi);

	/**
	 * @brief Computes the quotient and remainder of a large integer division.
	 *
	 * Computes `quotient = dividend / divisor` and updates `div_rem` with the remainder.
	 * Assumes `n > 0` and `n` is even.
	 */
	template<size_t n>
	BLST_HD  limb_t cc_quot_rem_n(limb_t* div_rem, const limb_t* divisor,
		limb_t quotient);

	///**
	// * @brief Computes the quotient and remainder of a large integer division.
	// *
	// * Computes `quotient = dividend / divisor` and updates `div_rem` with the remainder.
	// */
	//BLST_HD  limb_t cc_quot_rem_128(limb_t* div_rem, const limb_t* divisor,
	//	limb_t quotient)
	//{
	//	return cc_quot_rem_n<NLIMBS(128)>(div_rem, divisor, quotient);
	//}


	///**
	// * @brief Computes the quotient and remainder of a large integer division.
	// *
	// * Computes `quotient = dividend / divisor` and updates `div_rem` with the remainder.
	// */
	//BLST_HD  limb_t cc_quot_rem_64(limb_t* div_rem, const limb_t* divisor,
	//	limb_t quotient)
	//{
	//	return cc_quot_rem_n<NLIMBS(64)>(div_rem, divisor, quotient);
	//}


	/*
	 * Prefetching
	 */
	BLST_HD BLST_FORCEINLINE void cc_vec_prefetch(const void* ptr, size_t len)
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
	BLST_HD BLST_FORCEINLINE void cc_blst_sha256_hcopy(unsigned int dst[8], const unsigned int src[8])
	{
		size_t i;

		for (i = 0; i < 8; i++)
			dst[i] = src[i];
	}


	// performs the final compression function for sha256
	BLST_HD BLST_FORCEINLINE void cc_blst_sha256_emit(unsigned char md[32], const unsigned int h[8])
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
	BLST_HD BLST_FORCEINLINE void cc_blst_sha256_bcopy(void* dst_, const void* src_, size_t len)
	{
		unsigned char* dst = (unsigned char*)dst_;
		const unsigned char* src = (const unsigned char*)src_;
		size_t i;

		for (i = 0; i < len; i++)
			dst[i] = src[i];
	}


}

#endif