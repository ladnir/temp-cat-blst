#ifndef BLST_INTERFACE_H
#define BLST_INTERFACE_H

#include "backend_portable.h"
#include "backend_asm.h"
//BLST_EXTERN_C_BEGIN

namespace blst
{






#ifndef BLST_ASM
	BLST_HD BLST_FORCEINLINE void vec_prefetch(const void* ptr, size_t len)
	{
		cc_vec_prefetch(ptr, len);
	}
#else
	extern "C" void vec_prefetch(const void* ptr, size_t len); // linked by asm 
#endif



	BLST_HD BLST_HD BLST_FORCEINLINE void blst_sha256_block_data_order_v(unsigned int* v, const void* inp,
		size_t blocks)
	{
#if defined(__CUDA_ARCH__) || !defined(BLST_ASM) || !defined(BLST_DEFAULT_ASM)
		cc_blst_sha256_block_data_order(v, inp, blocks);
#else
		// asm
#   if (defined(__x86_64__) || defined(__x86_64) || defined(_M_X64)) &&  defined(__SHA__) /* -msha */ 
		blst_sha256_block_data_order_shaext(v, inp, blocks);
#	elif defined(__aarch64__) && defined(__ARM_FEATURE_CRYPTO) 
	#error "not set tested. need to call blst_sha256_block_armv8"
		blst_sha256_block_armv8(v, inp, blocks);
#	else
		#error "unknown asm target for SHA"
#	endif
#endif
	}


	// copies a 64 byte length array into dst_. i.e memcpy(dst, src, 64)
	BLST_HD BLST_HD BLST_FORCEINLINE void blst_sha256_hcopy_v(unsigned int dst[8], const unsigned int src[8])
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		blst_sha256_hcopy(dst, src);
#else
		cc_blst_sha256_hcopy(dst, src);
#endif
	}

	// performs the final compression function for sha256
	BLST_HD BLST_HD BLST_FORCEINLINE void blst_sha256_emit_v(unsigned char md[32], const unsigned int h[8])
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		blst_sha256_emit(md, h);
#else
		cc_blst_sha256_emit(md, h);
#endif
	}

	// copies a variable length array into dst_. i.e memcpy(dst, src, len)
	BLST_HD BLST_HD BLST_FORCEINLINE void blst_sha256_bcopy_v(void* dst_, const void* src_, size_t len)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		blst_sha256_bcopy(dst_, src_, len);
#else
		cc_blst_sha256_bcopy(dst_, src_, len);
#endif
	}


	BLST_HD BLST_FORCEINLINE limb_t div_3_limbs_v(const limb_t dividend_top[2], limb_t d_lo, limb_t d_hi)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return div_3_limbs(dividend_top, d_lo, d_hi);
#else
		return cc_div_3_limbs(dividend_top, d_lo, d_hi);
#endif
	}

	BLST_HD BLST_FORCEINLINE limb_t quot_rem_128_v(limb_t* quot_rem, const limb_t* divisor, limb_t quotient)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return quot_rem_128(quot_rem, divisor, quotient);
#else
		return cc_quot_rem_n<NLIMBS(128)>(quot_rem, divisor, quotient);
#endif
	}

	BLST_HD BLST_FORCEINLINE limb_t quot_rem_64_v(limb_t* quot_rem, const limb_t* divisor, limb_t quotient)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
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

		d_lo = (sizeof(z) == sizeof(limb_t)) ? 0 : z[z_len - 2];
		d_hi = z[z_len - 1];

		for (loop = z_len, z_len--; loop--;) {
			limb_t q = div_3_limbs_v(val + loop + z_len, d_lo, d_hi);
			(void)quot_rem_64_v(val + loop, z, q);
		}

		/* remainder is in low half of val[], quotient is in high */
	}

}

#endif // BLST_INTERFACE_H