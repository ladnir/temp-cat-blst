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
#include "blst/backend/backend_asm.h"

//BLST_EXTERN_C_BEGIN

namespace blst
{

	BLST_HD BLST_FORCEINLINE void vec_copy(void* restrict ret, const void* a, size_t num)
	{
		limb_t* rp = (limb_t*)ret;
		const limb_t* ap = (const limb_t*)a;
		size_t i;

		num /= sizeof(limb_t);

		for (i = 0; i < num; i++)
			rp[i] = ap[i];
	}

	BLST_HD BLST_INLINE void vec_zero(void* ret, size_t num)
	{
		volatile limb_t* rp = (volatile limb_t*)ret;
		size_t i;

		num /= sizeof(limb_t);

		for (i = 0; i < num; i++)
			rp[i] = 0;

		// do not optimize this away
#if defined(_MSC_VER)
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
			limb_t bi;
			volatile limb_t* rp = (limb_t*)ret;
			const limb_t* ap = (const limb_t*)a;
			const limb_t* bp = (const limb_t*)b;
			limb_t xorm, mask = (limb_t)0 - sel_a;
			size_t i;

			num /= sizeof(limb_t);

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
		const limb_t* ap = (const limb_t*)a;
		limb_t acc;
		size_t i;

#if defined(BLST_ASM) && !defined(__CUDA_ARCH__)
		if ((num & 15) == 0)
			return vec_is_zero_16x(a, num);
#endif

		num /= sizeof(limb_t);

		for (acc = 0, i = 0; i < num; i++)
			acc |= ap[i];

		return is_zero(acc);
	}

	BLST_HD BLST_FORCEINLINE bool_t vec_is_equal(const void* a, const void* b, size_t num)
	{
		const limb_t* ap = (const limb_t*)a;
		const limb_t* bp = (const limb_t*)b;
		limb_t acc;
		size_t i;

#if defined(BLST_ASM) && !defined(__CUDA_ARCH__)
		if ((num & 15) == 0)
			return vec_is_equal_16x(a, b, num);
#endif

		num /= sizeof(limb_t);

		for (acc = 0, i = 0; i < num; i++)
			acc |= ap[i] ^ bp[i];

		return is_zero(acc);
	}

	BLST_HD BLST_FORCEINLINE void vec_czero(void* ret, size_t num, bool_t cbit)
	{
		limb_t* rp = (limb_t*)ret;
		size_t i;
		limb_t mask;

		BLST_launder(cbit);
		mask = (limb_t)0 - (cbit ^ 1);

		num /= sizeof(limb_t);

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

	BLST_HD BLST_FORCEINLINE void vec_cswap(void* restrict a, void* restrict b, size_t num,
		bool_t cbit)
	{
		limb_t ai, * ap = (limb_t*)a;
		limb_t bi, * bp = (limb_t*)b;
		limb_t xorm, mask;
		size_t i;

		BLST_launder(cbit);
		mask = (limb_t)0 - cbit;

		num /= sizeof(limb_t);

		for (i = 0; i < num; i++) {
			xorm = ((ai = ap[i]) ^ (bi = bp[i])) & mask;
			ap[i] = ai ^ xorm;
			bp[i] = bi ^ xorm;
		}
	}

	//BLST_EXTERN_C_END


}

#endif /* __BLS12_381_ASM_VECT_H__ */