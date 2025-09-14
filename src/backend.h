#ifndef BLST_INTERFACE_H
#define BLST_INTERFACE_H

#include "backend_portable.h"
#include "backend_asm.h"
BLST_EXTERN_C_BEGIN

// montgomery multiply: ret = a * b * R^-1 mod p, where R montgomery radix
// optimized for "spare" modulus p, i.e. with many zero limbs
HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_mont_sparse_v256(vec256 ret, const vec256 a, const vec256 b,
	const vec256 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	mulx_mont_sparse_256(ret, a, b, p, n0);// linked by asm
#else
	cc_mul_mont_sparse_256(ret, a, b, p, n0);
#endif
}

// montgomery squaring: ret = a^2 * R^-1 mod p, where R montgomery radix
// optimized for "spare" modulus p, i.e. with many zero limbs
HOST_DEVICE_FUNC BLST_FORCEINLINE void sqr_mont_sparse_v256(vec256 ret, const vec256 a, const vec256 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	sqrx_mont_sparse_256(ret, a, p, n0);// linked by asm
#else
	cc_sqr_mont_sparse_256(ret, a, p, n0);
#endif
}

// montgomery multiply: ret = a * b * R^-1 mod p, where R montgomery radix
HOST_DEVICE_FUNC BLST_FORCEINLINE void redc_mont_v256(vec256 ret, const vec512 a, const vec256 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	redcx_mont_256(ret, a, p, n0);// linked by asm
#else
	cc_redc_mont_256(ret, a, p, n0);
#endif
}

// from montgomery domain: ret = a * R^-1 mod p, where R montgomery radix
HOST_DEVICE_FUNC BLST_FORCEINLINE void from_mont_v256(vec256 ret, const vec256 a, const vec256 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	fromx_mont_256(ret, a, p, n0);// linked by asm
#else
	cc_from_mont_256(ret, a, p, n0);
#endif
}


HOST_DEVICE_FUNC BLST_FORCEINLINE void add_mod_v256(vec256 ret, const vec256 a, const vec256 b, const vec256 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	add_mod_256(ret, a, b, p);// linked by asm, 
#else
	cc_add_mod_256(ret, a, b, p);
#endif

}

HOST_DEVICE_FUNC BLST_FORCEINLINE void sub_mod_v256(vec256 ret, const vec256 a, const vec256 b, const vec256 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	sub_mod_256(ret, a, b, p);// linked by asm
#else
	cc_sub_mod_256(ret, a, b, p);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_by_3_mod_v256(vec256 ret, const vec256 a, const vec256 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	mul_by_3_mod_256(ret, a, p);
#else
	cc_mul_by_3_mod_256(ret, a, p);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void cneg_mod_v256(vec256 ret, const vec256 a, bool_t flag, const vec256 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	cneg_mod_256(ret, a, flag, p);
#else
	cc_cneg_mod_256(ret, a, flag, p);
#endif // BLST_DEFAULT_ASM

}
HOST_DEVICE_FUNC BLST_FORCEINLINE void lshift_mod_v256(vec256 ret, const vec256 a, size_t count, const vec256 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	lshift_mod_256(ret, a, count, p);
#else
	cc_lshift_mod_256(ret, a, count, p);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void rshift_mod_v256(vec256 ret, const vec256 a, size_t count, const vec256 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	rshift_mod_256(ret, a, count, p);
#else
	cc_rshift_mod_256(ret, a, count, p);
#endif
}

//
//bool_t eucl_inverse_mod_v256(vec256 ret, const vec256 a, const vec256 p, const vec256 one)
//{
//#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
//    bool_t eucl_inverse_mod_256(vec256 ret, const vec256 a, const vec256 p, const vec256 one)
//#else
//	return cc_eucl_inverse_mod_256(ret, a, p, one);
//#endif
//
//}



HOST_DEVICE_FUNC BLST_FORCEINLINE limb_t check_mod_v256(const pow256 a, const vec256 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	return check_mod_256(a, p);
#else
	return cc_check_mod_256(a, p);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE limb_t add_n_check_mod_v256(pow256 ret, const pow256 a, const pow256 b,
	const vec256 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	return add_n_check_mod_256(ret, a, b, p);
#else
	return cc_add_n_check_mod_256(ret, a, b, p);
#endif
}

// performs modular subtraction with zero-checking - 
// it computes (a - b) mod p and returns whether the result is non-zero.
HOST_DEVICE_FUNC BLST_FORCEINLINE limb_t sub_n_check_mod_v256(pow256 ret, const pow256 a, const pow256 b,
	const vec256 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	return sub_n_check_mod_256(ret, a, b, p);
#else
	return cc_sub_n_check_mod_256(ret, a, b, p);
#endif
}


#ifndef BLST_ASM
HOST_DEVICE_FUNC BLST_FORCEINLINE void vec_prefetch(const void* ptr, size_t len)
{
	cc_vec_prefetch(ptr, len);
}
#else
void vec_prefetch(const void* ptr, size_t len); // linked by asm
#endif

// 
HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_mont_v384(vec384 ret, const vec384 a, const vec384 b,
	const vec384 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	mulx_mont_384(ret, a, b, p, n0);
#else
	cc_mul_mont_384(ret, a, b, p, n0);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void sqr_mont_v384(vec384 ret, const vec384 a, const vec384 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	sqrx_mont_384(ret, a, p, n0);
#else
	cc_sqr_mont_384(ret, a, p, n0);
#endif
}

//HOST_DEVICE_FUNC BLST_FORCEINLINE void sqr_n_mul_mont_v384(vec384 ret, const vec384 a, size_t count,
//	const vec384 p, limb_t n0, const vec384 b)
//{
//#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
//	sqrx_n_mul_mont_384(ret, a, count, p, n0, b);
//#else
//	cc_sqr_n_mul_mont_384(ret, a, count, p, n0, b);
//#endif
//}

/*
 * "383"-bit variant omits full reductions at the ends of squarings,
 * which results in up to ~15% improvement. [One can improve further
 * by omitting full reductions even after multiplications and
 * performing final reduction at the very end of the chain.]
 * 
 * see sqr_n_mul_fp
 */
HOST_DEVICE_FUNC BLST_FORCEINLINE void sqr_n_mul_mont_v383(vec384 ret, const vec384 a, size_t count,
	const vec384 p, limb_t n0, const vec384 b)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	sqrx_n_mul_mont_383(ret, a, count, p, n0, b);
#else
	cc_sqr_n_mul_mont_384(ret, a, count, p, n0, b);
#endif
}

//void mul_384(vec768 ret, const vec384 a, const vec384 b)
//{
//#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
//	mulx_384(ret, a, b);
//#else
//	cc_mul_384(ret, a, b);
//#endif
//}
//void sqr_384(vec768 ret, const vec384 a);

HOST_DEVICE_FUNC BLST_FORCEINLINE void redc_mont_v384(vec384 ret, const vec768 a, const vec384 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
#ifdef BLST_ADX
	redcx_mont_384(ret, a, p, n0);
#else
	redc_mont_384(ret, a, p, n0);
#endif
#else
	cc_redc_mont_384(ret, a, p, n0);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void from_mont_v384(vec384 ret, const vec384 a, const vec384 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	fromx_mont_384(ret, a, p, n0);
#else
	cc_from_mont_384(ret, a, p, n0);
#endif
}


HOST_DEVICE_FUNC BLST_FORCEINLINE limb_t sgn0_pty_mont_v384(const vec384 a, const vec384 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	return sgn0x_pty_mont_384(a, p, n0);
#else
	return cc_sgn0_pty_mont_384(a, p, n0);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE limb_t sgn0_pty_mont_v384x(const vec384x a, const vec384 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	return sgn0x_pty_mont_384x(a, p, n0);
#else
	return cc_sgn0_pty_mont_384x(a, p, n0);
#endif
}
HOST_DEVICE_FUNC BLST_FORCEINLINE limb_t sgn0_pty_mod_v384(const vec384 a, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	return sgn0_pty_mod_384(a, p);
#else
	return cc_sgn0_pty_mod_384(a, p);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE limb_t sgn0_pty_mod_v384x(const vec384x a, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	return sgn0_pty_mod_384x(a, p);
#else
	return cc_sgn0_pty_mod_384x(a, p);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void add_mod_v384(vec384 ret, const vec384 a, const vec384 b, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	add_mod_384(ret, a, b, p);
#else
	cc_add_mod_384(ret, a, b, p);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void sub_mod_v384(vec384 ret, const vec384 a, const vec384 b, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	sub_mod_384(ret, a, b, p);
#else
	cc_sub_mod_n(ret, a, b, p, NLIMBS(384));
#endif
}
HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_by_8_mod_v384(vec384 ret, const vec384 a, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	mul_by_8_mod_384(ret, a, p);
#else
	cc_mul_by_8_mod_384(ret, a, p);
#endif

}
HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_by_3_mod_v384(vec384 ret, const vec384 a, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	mul_by_3_mod_384(ret, a, p);
#else
	cc_mul_by_3_mod_384(ret, a, p);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void cneg_mod_v384(vec384 ret, const vec384 a, bool_t flag, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	cneg_mod_384(ret, a, flag, p);
#else
	cc_cneg_mod_384(ret, a, flag, p);
#endif
}

void lshift_mod_v384(vec384 ret, const vec384 a, size_t count, const vec384 p);

HOST_DEVICE_FUNC BLST_FORCEINLINE void rshift_mod_v384(vec384 ret, const vec384 a, size_t count, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	rshift_mod_384(ret, a, count, p);
#else
	cc_rshift_mod_384(ret, a, count, p);
#endif
}
HOST_DEVICE_FUNC BLST_FORCEINLINE void div_by_2_mod_v384(vec384 ret, const vec384 a, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	div_by_2_mod_384(ret, a, p);
#else
	cc_div_by_2_mod_384(ret, a, p);
#endif
}


HOST_DEVICE_FUNC BLST_FORCEINLINE void ct_inverse_mod_v384(vec768 ret, const vec384 inp, const vec384 mod,
	const vec384 modx)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	ctx_inverse_mod_384(ret, inp, mod, modx);
#else
	cc_ct_inverse_mod_384(ret, inp, mod, modx);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void ct_inverse_mod_v256(vec512 ret, const vec256 inp, const vec256 mod,
	const vec256 modx)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	ct_inverse_mod_256(ret, inp, mod, modx);
#else	
	cc_ct_inverse_mod_256(ret, inp, mod, modx);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE bool_t ct_is_square_mod_v384(const vec384 inp, const vec384 mod)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	return ct_is_square_mod_384(inp, mod);
#else
	return cc_ct_is_square_mod_384(inp, mod);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_mont_v384x(vec384x ret, const vec384x a, const vec384x b,
	const vec384 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	mulx_mont_384x(ret, a, b, p, n0);
#else
	cc_mul_mont_384x(ret, a, b, p, n0);
#endif
}
HOST_DEVICE_FUNC BLST_FORCEINLINE void sqr_mont_v384x(vec384x ret, const vec384x a, const vec384 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	sqrx_mont_384x(ret, a, p, n0);
#else
	cc_sqr_mont_384x(ret, a, p, n0);
#endif

}
HOST_DEVICE_FUNC BLST_FORCEINLINE void sqr_mont_v382x(vec384x ret, const vec384x a, const vec384 p, limb_t n0)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	sqrx_mont_382x(ret, a, p, n0);
#else
	cc_sqr_mont_382x(ret, a, p, n0);
#endif
}

#if defined(BLST_ASM) && !defined(__CUDA_ARCH__)

HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_v382x(vec768 ret[2], const vec384x a, const vec384x b, const vec384 p)
{
#ifdef BLST_ADX
	mulx_382x(ret, a, b, p);
#else
	mul_382x(ret, a, b, p);
#endif
}
HOST_DEVICE_FUNC BLST_FORCEINLINE void sqr_v382x(vec768 ret[2], const vec384x a, const vec384 p)
{
#ifdef BLST_ADX
	sqrx_382x(ret, a, p);
#else
	sqr_382x(ret, a, p);
#endif
}

#endif

HOST_DEVICE_FUNC BLST_FORCEINLINE void add_mod_v384x(vec384x ret, const vec384x a, const vec384x b,
	const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	add_mod_384x(ret, a, b, p);
#else
	cc_add_mod_384x(ret, a, b, p);
#endif
}
HOST_DEVICE_FUNC BLST_FORCEINLINE void sub_mod_v384x(vec384x ret, const vec384x a, const vec384x b,
	const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	sub_mod_384x(ret, a, b, p);
#else
	cc_sub_mod_384x(ret, a, b, p);
#endif
}
HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_by_8_mod_v384x(vec384x ret, const vec384x a, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	mul_by_8_mod_384x(ret, a, p);
#else
	cc_mul_by_8_mod_384x(ret, a, p);
#endif
}
HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_by_3_mod_v384x(vec384x ret, const vec384x a, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	mul_by_3_mod_384x(ret, a, p);
#else
	cc_mul_by_3_mod_384x(ret, a, p);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_by_1_plus_i_mod_v384x(vec384x ret, const vec384x a, const vec384 p)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	mul_by_1_plus_i_mod_384x(ret, a, p);
#else
	cc_mul_by_1_plus_i_mod_384x(ret, a, p);
#endif
}




HOST_DEVICE_FUNC HOST_DEVICE_FUNC BLST_FORCEINLINE void blst_sha256_block_data_order_v(unsigned int* v, const void* inp,
	size_t blocks)
{
#if (defined(__x86_64__) || defined(__x86_64) || defined(_M_X64)) && \
     defined(__SHA__) /* -msha */ 
#error "not set tested. need to call blst_sha256_block_data_order_shaext"
//# define sha256_block_data_order blst_sha256_block_data_order_shaext
#elif defined(__aarch64__) && \
      defined(__ARM_FEATURE_CRYPTO) 
#error "not set tested. need to call blst_sha256_block_armv8"
//# define sha256_block_data_order blst_sha256_block_armv8
#else
//# define sha256_block_data_order blst_sha256_block_data_order_v
#endif

#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	blst_sha256_block_data_order(v, inp, blocks);
#else
	cc_blst_sha256_block_data_order(v, inp, blocks);
#endif
}


// copies a 64 byte length array into dst_. i.e memcpy(dst, src, 64)
HOST_DEVICE_FUNC HOST_DEVICE_FUNC BLST_FORCEINLINE void blst_sha256_hcopy_v(unsigned int dst[8], const unsigned int src[8])
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	blst_sha256_hcopy(dst, src);
#else
	cc_blst_sha256_hcopy(dst, src);
#endif
}

// performs the final compression function for sha256
HOST_DEVICE_FUNC HOST_DEVICE_FUNC BLST_FORCEINLINE void blst_sha256_emit_v(unsigned char md[32], const unsigned int h[8])
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	blst_sha256_emit(md, h);
#else
	cc_blst_sha256_emit(md, h);
#endif
}

// copies a variable length array into dst_. i.e memcpy(dst, src, len)
HOST_DEVICE_FUNC HOST_DEVICE_FUNC BLST_FORCEINLINE void blst_sha256_bcopy_v(void* dst_, const void* src_, size_t len)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	blst_sha256_bcopy(dst_, src_, len);
#else
	cc_blst_sha256_bcopy(dst_, src_, len);
#endif
}



/*
 * C subroutines
 */
HOST_DEVICE_FUNC void exp_mont_384(vec384 out, const vec384 inp, const byte* pow,
	size_t pow_bits, const vec384 p, limb_t n0);
HOST_DEVICE_FUNC void exp_mont_384x(vec384x out, const vec384x inp, const byte* pow,
	size_t pow_bits, const vec384 p, limb_t n0);



HOST_DEVICE_FUNC BLST_FORCEINLINE void cneg_mod_v384x(vec384x ret, const vec384x a, bool_t flag,
	const vec384 p)
{
	cneg_mod_v384(ret[0], a[0], flag, p);
	cneg_mod_v384(ret[1], a[1], flag, p);
}


HOST_DEVICE_FUNC BLST_FORCEINLINE limb_t div_3_limbs_v(const limb_t dividend_top[2], limb_t d_lo, limb_t d_hi)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	return div_3_limbs(dividend_top, d_lo, d_hi);
#else
	return cc_div_3_limbs(dividend_top, d_lo, d_hi);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE limb_t quot_rem_128_v(limb_t* quot_rem, const limb_t* divisor, limb_t quotient)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	return quot_rem_128(quot_rem, divisor, quotient);
#else
	return cc_quot_rem_128(quot_rem, divisor, quotient);
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE limb_t quot_rem_64_v(limb_t* quot_rem, const limb_t* divisor, limb_t quotient)
{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	return quot_rem_64(quot_rem, divisor, quotient);
#else
	return cc_quot_rem_64(quot_rem, divisor, quotient);
#endif
}

/*
 * Divide 255-bit |val| by z^2 yielding 128-bit quotient and remainder in place.
 */
HOST_DEVICE_FUNC BLST_FORCEINLINE void div_by_zz(limb_t val[])
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
HOST_DEVICE_FUNC BLST_FORCEINLINE void div_by_z(limb_t val[])
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

BLST_EXTERN_C_END
#endif // BLST_INTERFACE_H