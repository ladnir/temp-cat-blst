/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_VECT_H__
#define __BLS12_381_ASM_VECT_H__

#include "blst/config.h"
#include <stddef.h>
#include "blst/common/defines.h"

//BLST_EXTERN_C_BEGIN

#if defined(BLST_ASM) && !defined(__CUDA_ARCH__)

extern "C"
{
	void vec_select_32(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm
	void vec_select_48(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm
	void vec_select_96(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm
	void vec_select_144(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm
	void vec_select_192(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm
	void vec_select_288(void* ret, const void* a, const void* b, blst::bool_t sel_a);	// asm

	blst::bool_t vec_is_zero_16x(const void* a, size_t num);	// asm
	blst::bool_t vec_is_equal_16x(const void* a, const void* b, size_t num);	// asm
}
#endif

namespace blst
{
	BLST_HD BLST_FORCEINLINE void vec_copy(void* BLST_restrict ret, const void* a, size_t num)
	{
		unsigned char* rp = (unsigned char*)ret;
		const unsigned char* ap = (const unsigned char*)a;
		size_t i;

		for (i = 0; i < num; i++)
			rp[i] = ap[i];
	}


	BLST_HD BLST_INLINE void vec_zero(void* ret, size_t num)
	{
		volatile unsigned char* rp = (volatile unsigned char*)ret;
		size_t i;

		for (i = 0; i < num; i++)
			rp[i] = 0;

		// do not optimize this away
#if defined(_MSC_VER) || defined(__CUDACC__)
	// Compiler barrier for MSVC
		BLST_launder(i);
#elif defined(__GNUC__) || defined(__clang__) 
		__asm__ __volatile__("" : : "r"(ret) : "memory");
#endif
	}



	BLST_HD BLST_FORCEINLINE void vec_select(void* ret, const void* a, const void* b,
		size_t num, bool_t sel_a)
	{
		BLST_launder(sel_a);
#if defined(BLST_ASM) && !defined(__CUDA_ARCH__)
		if (num == 32)          vec_select_32(ret, a, b, sel_a);
		else if (num == 48)     vec_select_48(ret, a, b, sel_a);
		else if (num == 96)     vec_select_96(ret, a, b, sel_a);
		else if (num == 144)    vec_select_144(ret, a, b, sel_a);
		else if (num == 192)    vec_select_192(ret, a, b, sel_a);
		else if (num == 288)    vec_select_288(ret, a, b, sel_a);
#else
		if (0);
#endif
		else {
			unsigned char bi;
			volatile unsigned char* rp = (volatile unsigned char*)ret;
			const unsigned char* ap = (const unsigned char*)a;
			const unsigned char* bp = (const unsigned char*)b;
			unsigned char xorm, mask = (unsigned char)0 - (unsigned char)sel_a;
			size_t i;

			for (i = 0; i < num; i++) {
				xorm = (ap[i] ^ (bi = bp[i])) & mask;
				rp[i] = bi ^ xorm;
			}
		}
	}

	BLST_HD BLST_FORCEINLINE bool_t is_zero(limb_t l)
	{
		limb_t ret = (~l & (l - 1)) >> (BLST_LIMB_SIZE - 1);
		BLST_launder(ret);
		return ret;
	}

	BLST_HD BLST_FORCEINLINE bool_t vec_is_zero(const void* a, size_t num)
	{
		const unsigned char* ap = (const unsigned char*)a;
		unsigned char acc;
		size_t i;

#if defined(BLST_ASM) && !defined(__CUDA_ARCH__)
		if ((num & 15) == 0)
			return vec_is_zero_16x(a, num);
#endif

		for (acc = 0, i = 0; i < num; i++)
			acc |= ap[i];

		limb_t ret = ((limb_t)(acc)-1) >> (BLST_LIMB_SIZE - 1);
		BLST_launder(ret);
		return ret;
	}

	BLST_HD BLST_FORCEINLINE bool_t vec_is_equal(const void* a, const void* b, size_t num)
	{
		const unsigned char* ap = (const unsigned char*)a;
		const unsigned char* bp = (const unsigned char*)b;
		unsigned char acc;
		size_t i;

#if defined(BLST_ASM) && !defined(__CUDA_ARCH__)
		if ((num & 15) == 0)
			return vec_is_equal_16x(a, b, num);
#endif

		for (acc = 0, i = 0; i < num; i++)
			acc |= ap[i] ^ bp[i];

		limb_t ret = ((limb_t)(acc)-1) >> (BLST_LIMB_SIZE - 1);
		BLST_launder(ret);
		return ret;
	}

	BLST_HD BLST_FORCEINLINE void vec_czero(void* ret, size_t num, bool_t cbit)
	{
		unsigned char* rp = (unsigned char*)ret;
		size_t i;
		unsigned char mask;

		BLST_launder(cbit);
		mask = (unsigned char)0 - (unsigned char)(cbit ^ 1);

		for (i = 0; i < num; i++)
			rp[i] &= mask;
	}


	BLST_HD BLST_FORCEINLINE bool_t is_bit_set(const byte* v, size_t i)
	{
		bool_t ret = (v[i / 8] >> (i % 8)) & 1;
		BLST_launder(ret);
		return ret;
	}

	BLST_HD BLST_FORCEINLINE bool_t byte_is_zero(unsigned char c)
	{
		limb_t ret = ((limb_t)(c)-1) >> (BLST_LIMB_SIZE - 1);
		BLST_launder(ret);
		return ret;
	}

	BLST_HD BLST_FORCEINLINE bool_t bytes_are_zero(const unsigned char* a, size_t num)
	{
		unsigned char acc;
		size_t i;

		for (acc = 0, i = 0; i < num; i++)
			acc |= a[i];

		return byte_is_zero(acc);
	}

	BLST_HD BLST_FORCEINLINE void vec_cswap(void* BLST_restrict a, void* BLST_restrict b, size_t num,
		bool_t cbit)
	{
		unsigned char ai, bi;
		unsigned char* ap = (unsigned char*)a;
		unsigned char* bp = (unsigned char*)b;
		unsigned char xorm, mask;
		size_t i;

		BLST_launder(cbit);
		mask = (unsigned char)0 - (unsigned char)cbit;

		for (i = 0; i < num; i++) {
			xorm = ((ai = ap[i]) ^ (bi = bp[i])) & mask;
			ap[i] = ai ^ xorm;
			bp[i] = bi ^ xorm;
		}
	}

	//BLST_EXTERN_C_END


}

#endif /* __BLS12_381_ASM_VECT_H__ */
