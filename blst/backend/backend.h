/*
 * Copyright Supranational LLC
 * Copyright 2026 Peter Rindal
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 *
 * This CatBLST header is locally composed. Its fixed-width arithmetic and
 * division lineage includes the upstream BLST src/vect.h and src/vect.c.
 */

#ifndef BLST_BACKEND_H
#define BLST_BACKEND_H

#include "backend_portable.h"
#include "backend_asm.h"
#ifdef BLST_CUDA
#include "backend_cuda.h"
#endif
namespace blst
{

    /*
     * Compile-time-selected scalar operations. The cc_* domain remains the
     * pure C/C++ implementation. These wrappers add no runtime dispatch:
     * optimized host builds call the selected assembly symbol directly,
     * while no-assembly, non-default-assembly, and device builds call cc_*.
     */

    BLST_HD BLST_FORCEINLINE void add_mod_v256(limb_t* ret, const limb_t* a,
        const limb_t* b, const limb_t* p)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ::add_mod_256(ret, a, b, p);
#else
        cc_add_mod_n<NLIMBS(256)>(ret, a, b, p);
#endif
    }

    BLST_HD BLST_FORCEINLINE void sub_mod_v256(limb_t* ret, const limb_t* a,
        const limb_t* b, const limb_t* p)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ::sub_mod_256(ret, a, b, p);
#else
        cc_sub_mod_n<NLIMBS(256)>(ret, a, b, p);
#endif
    }

    BLST_HD BLST_FORCEINLINE void mul_mont_sparse_v256(limb_t* ret,
        const limb_t* a, const limb_t* b, const limb_t* p, limb_t n0)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
# if defined(BLST_ADX)
        ::mulx_mont_sparse_256(ret, a, b, p, n0);
# else
        ::mul_mont_sparse_256(ret, a, b, p, n0);
# endif
#else
        cc_mul_mont_n<NLIMBS(256)>(ret, a, b, p, n0);
#endif
    }

    BLST_HD BLST_FORCEINLINE void sqr_mont_sparse_v256(limb_t* ret,
        const limb_t* a, const limb_t* p, limb_t n0)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
# if defined(BLST_ADX)
        ::sqrx_mont_sparse_256(ret, a, p, n0);
# else
        ::sqr_mont_sparse_256(ret, a, p, n0);
# endif
#else
        cc_mul_mont_n<NLIMBS(256)>(ret, a, a, p, n0);
#endif
    }

    BLST_HD BLST_FORCEINLINE void redc_mont_v256(limb_t* ret,
        const limb_t* a, const limb_t* p, limb_t n0)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
# if defined(BLST_ADX)
        ::redcx_mont_256(ret, a, p, n0);
# else
        ::redc_mont_256(ret, a, p, n0);
# endif
#else
        cc_redc_mont_n<NLIMBS(256)>(ret, a, p, n0);
#endif
    }

    BLST_HD BLST_FORCEINLINE void from_mont_v256(limb_t* ret,
        const limb_t* a, const limb_t* p, limb_t n0)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
# if defined(BLST_ADX)
        ::fromx_mont_256(ret, a, p, n0);
# else
        ::from_mont_256(ret, a, p, n0);
# endif
#else
        cc_from_mont_n<NLIMBS(256)>(ret, a, p, n0);
#endif
    }

    BLST_HD BLST_FORCEINLINE void cneg_mod_v256(limb_t* ret, const limb_t* a,
        bool_t flag, const limb_t* p)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ::cneg_mod_256(ret, a, flag, p);
#else
        cc_cneg_mod_n<NLIMBS(256)>(ret, a, flag, p);
#endif
    }

    BLST_HD BLST_FORCEINLINE void lshift_mod_v256(limb_t* ret,
        const limb_t* a, size_t count, const limb_t* p)
    {
        if (count == 0) {
            for (size_t i = 0; i < NLIMBS(256); ++i)
                ret[i] = a[i];
            return;
        }
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ::lshift_mod_256(ret, a, count, p);
#else
        cc_lshift_mod_n<NLIMBS(256)>(ret, a, count, p);
#endif
    }

    BLST_HD BLST_FORCEINLINE void rshift_mod_v256(limb_t* ret,
        const limb_t* a, size_t count, const limb_t* p)
    {
        if (count == 0) {
            for (size_t i = 0; i < NLIMBS(256); ++i)
                ret[i] = a[i];
            return;
        }
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ::rshift_mod_256(ret, a, count, p);
#else
        cc_rshift_mod_n<NLIMBS(256)>(ret, a, count, p);
#endif
    }

    BLST_HD BLST_FORCEINLINE void ct_inverse_mod_v256(limb_t* ret,
        const limb_t* inp, const limb_t* mod, const limb_t* modx)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ::ct_inverse_mod_256(ret, inp, mod, modx);
#else
        cc_ct_inverse_mod_n<NLIMBS(256)>(ret, inp, mod, modx);
#endif
    }

    BLST_HD BLST_FORCEINLINE void add_mod_v384(limb_t* ret, const limb_t* a,
        const limb_t* b, const limb_t* p)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ::add_mod_384(ret, a, b, p);
#else
        cc_add_mod_n<NLIMBS(384)>(ret, a, b, p);
#endif
    }

    BLST_HD BLST_FORCEINLINE void sub_mod_v384(limb_t* ret, const limb_t* a,
        const limb_t* b, const limb_t* p)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ::sub_mod_384(ret, a, b, p);
#else
        cc_sub_mod_n<NLIMBS(384)>(ret, a, b, p);
#endif
    }

    BLST_HD BLST_FORCEINLINE void mul_mont_v384(limb_t* ret, const limb_t* a,
        const limb_t* b, const limb_t* p, limb_t n0)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
# if defined(BLST_ADX)
        ::mulx_mont_384(ret, a, b, p, n0);
# else
        ::mul_mont_384(ret, a, b, p, n0);
# endif
#else
        cc_mul_mont_n<NLIMBS(384)>(ret, a, b, p, n0);
#endif
    }

    BLST_HD BLST_FORCEINLINE void sqr_mont_v384(limb_t* ret, const limb_t* a,
        const limb_t* p, limb_t n0)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
# if defined(BLST_ADX)
        ::sqrx_mont_384(ret, a, p, n0);
# else
        ::sqr_mont_384(ret, a, p, n0);
# endif
#else
        cc_mul_mont_n<NLIMBS(384)>(ret, a, a, p, n0);
#endif
    }

    BLST_HD BLST_FORCEINLINE void redc_mont_v384(limb_t* ret,
        const limb_t* a, const limb_t* p, limb_t n0)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
# if defined(BLST_ADX)
        ::redcx_mont_384(ret, a, p, n0);
# else
        ::redc_mont_384(ret, a, p, n0);
# endif
#else
        cc_redc_mont_n<NLIMBS(384)>(ret, a, p, n0);
#endif
    }

    BLST_HD BLST_FORCEINLINE void from_mont_v384(limb_t* ret,
        const limb_t* a, const limb_t* p, limb_t n0)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
# if defined(BLST_ADX)
        ::fromx_mont_384(ret, a, p, n0);
# else
        ::from_mont_384(ret, a, p, n0);
# endif
#else
        cc_from_mont_n<NLIMBS(384)>(ret, a, p, n0);
#endif
    }

    BLST_HD BLST_FORCEINLINE void cneg_mod_v384(limb_t* ret, const limb_t* a,
        bool_t flag, const limb_t* p)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ::cneg_mod_384(ret, a, flag, p);
#else
        cc_cneg_mod_n<NLIMBS(384)>(ret, a, flag, p);
#endif
    }

    BLST_HD BLST_FORCEINLINE void lshift_mod_v384(limb_t* ret,
        const limb_t* a, size_t count, const limb_t* p)
    {
        if (count == 0) {
            for (size_t i = 0; i < NLIMBS(384); ++i)
                ret[i] = a[i];
            return;
        }
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ::lshift_mod_384(ret, a, count, p);
#else
        cc_lshift_mod_n<NLIMBS(384)>(ret, a, count, p);
#endif
    }

    BLST_HD BLST_FORCEINLINE void rshift_mod_v384(limb_t* ret,
        const limb_t* a, size_t count, const limb_t* p)
    {
        if (count == 0) {
            for (size_t i = 0; i < NLIMBS(384); ++i)
                ret[i] = a[i];
            return;
        }
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ::rshift_mod_384(ret, a, count, p);
#else
        cc_rshift_mod_n<NLIMBS(384)>(ret, a, count, p);
#endif
    }

    BLST_HD BLST_FORCEINLINE void ct_inverse_mod_v384(limb_t* ret,
        const limb_t* inp, const limb_t* mod, const limb_t* modx)
    {
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
# if defined(BLST_ADX)
        ::ctx_inverse_mod_384(ret, inp, mod, modx);
# else
        ::ct_inverse_mod_384(ret, inp, mod, modx);
# endif
#else
        cc_ct_inverse_mod_n<NLIMBS(384)>(ret, inp, mod, modx);
#endif
    }






#ifndef BLST_ASM
	BLST_HD BLST_FORCEINLINE void vec_prefetch(const void* ptr, size_t len)
	{
		cc_vec_prefetch(ptr, len);
	}
#else
	extern "C" void vec_prefetch(const void* ptr, size_t len); // linked by asm
#endif



	BLST_HD BLST_FORCEINLINE void blst_sha256_block_data_order_v(unsigned int* v, const void* inp,
		size_t blocks)
	{
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		blst_sha256_block_data_order(v, inp, blocks);
#else
		cc_blst_sha256_block_data_order(v, inp, blocks);
#endif
	}


	// Copies the eight-word (32-byte) SHA-256 state.
	BLST_HD BLST_FORCEINLINE void blst_sha256_hcopy_v(unsigned int dst[8], const unsigned int src[8])
	{
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		blst_sha256_hcopy(dst, src);
#else
		cc_blst_sha256_hcopy(dst, src);
#endif
	}

	// performs the final compression function for sha256
	BLST_HD BLST_FORCEINLINE void blst_sha256_emit_v(unsigned char md[32], const unsigned int h[8])
	{
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		blst_sha256_emit(md, h);
#else
		cc_blst_sha256_emit(md, h);
#endif
	}

	// copies a variable length array into dst_. i.e memcpy(dst, src, len)
	BLST_HD BLST_FORCEINLINE void blst_sha256_bcopy_v(void* dst_, const void* src_, size_t len)
	{
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		blst_sha256_bcopy(dst_, src_, len);
#else
		cc_blst_sha256_bcopy(dst_, src_, len);
#endif
	}


	BLST_HD BLST_FORCEINLINE limb_t div_3_limbs_v(const limb_t dividend_top[2], limb_t d_lo, limb_t d_hi)
	{
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return div_3_limbs(dividend_top, d_lo, d_hi);
#else
		return cc_div_3_limbs(dividend_top, d_lo, d_hi);
#endif
	}

	BLST_HD BLST_FORCEINLINE limb_t quot_rem_128_v(limb_t* quot_rem, const limb_t* divisor, limb_t quotient)
	{
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return quot_rem_128(quot_rem, divisor, quotient);
#else
		return cc_quot_rem_n<NLIMBS(128)>(quot_rem, divisor, quotient);
#endif
	}

	BLST_HD BLST_FORCEINLINE limb_t quot_rem_64_v(limb_t* quot_rem, const limb_t* divisor, limb_t quotient)
	{
#if defined(BLST_ASM) && defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return quot_rem_64(quot_rem, divisor, quotient);
#else
		return cc_quot_rem_n<NLIMBS(64)>(quot_rem, divisor, quotient);
#endif
	}



	/*
	 * C subroutines
	 */

	/*
	 * Divide 255-bit |val| by z^2 yielding 128-bit quotient and remainder in place.
	 */
	BLST_HD BLST_INLINE void div_by_zz(limb_t val[])
	{
		const limb_t zz[] = { TO_LIMB_T(0x0000000100000000),
									 TO_LIMB_T(0xac45a4010001a402) };
		size_t loop, zz_len = sizeof(zz) / sizeof(zz[0]);
		limb_t d_lo, d_hi;
		d_lo = zz[zz_len - 2];
		d_hi = zz[zz_len - 1];
		for (loop = zz_len, zz_len--; loop--;) {
			limb_t q = div_3_limbs_v(val + loop + zz_len, d_lo, d_hi);
			(void)quot_rem_128_v(val + loop, zz, q);
		}

		/* remainder is in low half of val[], quotient is in high */
	}

	/*
	 * Divide 128-bit |val| by z yielding 64-bit quotient and remainder in place.
	 */
	BLST_HD BLST_INLINE void div_by_z(limb_t* val)
	{
		const limb_t z[] = { TO_LIMB_T(0xd201000000010000) };
		size_t loop, z_len = sizeof(z) / sizeof(z[0]);
		limb_t d_lo, d_hi;

		d_lo = 0;
		d_hi = z[z_len - 1];

		for (loop = z_len, z_len--; loop--;) {
			limb_t q = div_3_limbs_v(val + loop + z_len, d_lo, d_hi);
			(void)quot_rem_64_v(val + loop, z, q);
		}

		/* remainder is in low half of val[], quotient is in high */
	}

}

#endif // BLST_BACKEND_H
