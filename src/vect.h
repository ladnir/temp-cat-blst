/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_VECT_H__
#define __BLS12_381_ASM_VECT_H__

#include "config.h"
#include <stddef.h>
#include "defines.h"



static inline void vec_copy(void* restrict ret, const void* a, size_t num)
{
	limb_t* rp = (limb_t*)ret;
	const limb_t* ap = (const limb_t*)a;
	size_t i;

	num /= sizeof(limb_t);

	for (i = 0; i < num; i++)
		rp[i] = ap[i];
}

inline void vec_zero(void* ret, size_t num)
{
	volatile limb_t* rp = (volatile limb_t*)ret;
	size_t i;

	num /= sizeof(limb_t);

	for (i = 0; i < num; i++)
		rp[i] = 0;

#if defined(__GNUC__) || defined(__clang__)
	__asm__ __volatile__("" : : "r"(ret) : "memory");
#endif
}

#ifdef BLST_ASM

/* ret = bit ? a : b */
void vec_select_32(void* ret, const void* a, const void* b, bool_t sel_a);	// asm
void vec_select_48(void* ret, const void* a, const void* b, bool_t sel_a);	// asm
void vec_select_96(void* ret, const void* a, const void* b, bool_t sel_a);	// asm
void vec_select_144(void* ret, const void* a, const void* b, bool_t sel_a);	// asm
void vec_select_192(void* ret, const void* a, const void* b, bool_t sel_a);	// asm
void vec_select_288(void* ret, const void* a, const void* b, bool_t sel_a);	// asm

#endif // BLST_ASM


inline void vec_select(void* ret, const void* a, const void* b,
	size_t num, bool_t sel_a)
{
	BLST_launder(sel_a);
#ifdef BLST_ASM
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

inline bool_t is_zero(limb_t l)
{
	limb_t ret = (~l & (l - 1)) >> (BLST_LIMB_SIZE - 1);
	BLST_launder(ret);
	return ret;
}

inline bool_t vec_is_zero(const void* a, size_t num)
{
	const limb_t* ap = (const limb_t*)a;
	limb_t acc;
	size_t i;

#ifdef BLST_ASM
	bool_t vec_is_zero_16x(const void* a, size_t num);
	if ((num & 15) == 0)
		return vec_is_zero_16x(a, num);
#endif

	num /= sizeof(limb_t);

	for (acc = 0, i = 0; i < num; i++)
		acc |= ap[i];

	return is_zero(acc);
}

inline bool_t vec_is_equal(const void* a, const void* b, size_t num)
{
	const limb_t* ap = (const limb_t*)a;
	const limb_t* bp = (const limb_t*)b;
	limb_t acc;
	size_t i;

#ifdef BLST_ASM
	bool_t vec_is_equal_16x(const void* a, const void* b, size_t num);
	if ((num & 15) == 0)
		return vec_is_equal_16x(a, b, num);
#endif

	num /= sizeof(limb_t);

	for (acc = 0, i = 0; i < num; i++)
		acc |= ap[i] ^ bp[i];

	return is_zero(acc);
}

inline void vec_czero(void* ret, size_t num, bool_t cbit)
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


#include "no_asm.h"

///*
// * Assembly subroutines...
// */
//#if defined(__ADX__) /* e.g. -march=broadwell */ && !defined(__BLST_PORTABLE__)\
//                                                 && !defined(__BLST_NO_ASM__)
//# define mul_mont_sparse_256 mulx_mont_sparse_256
//# define sqr_mont_sparse_256 sqrx_mont_sparse_256
//# define from_mont_256 fromx_mont_256
//# define redc_mont_256 redcx_mont_256
//# define mul_mont_384 mulx_mont_384
//# define sqr_mont_384 sqrx_mont_384
//# define sqr_n_mul_mont_384 sqrx_n_mul_mont_384
//# define sqr_n_mul_mont_383 sqrx_n_mul_mont_383
//# define mul_384 mulx_384
//# define sqr_384 sqrx_384
//# define redc_mont_384 redcx_mont_384
//# define from_mont_384 fromx_mont_384
//# define sgn0_pty_mont_384 sgn0x_pty_mont_384
//# define sgn0_pty_mont_384x sgn0x_pty_mont_384x
//# define ct_inverse_mod_384 ctx_inverse_mod_384
//#endif

inline void mul_mont_sparse_v256(vec256 ret, const vec256 a, const vec256 b,
	const vec256 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	mulx_mont_sparse_256(ret, a, b, p, n0);// linked by asm
#else
	cc_mul_mont_sparse_256(ret, a, b, p, n0);
#endif

}
inline void sqr_mont_sparse_v256(vec256 ret, const vec256 a, const vec256 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	sqrx_mont_sparse_256(ret, a, p, n0);// linked by asm
#else
	cc_sqr_mont_sparse_256(ret, a, p, n0);
#endif
}
inline void redc_mont_v256(vec256 ret, const vec512 a, const vec256 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	redcx_mont_256(ret, a, p, n0);// linked by asm
#else
	cc_redc_mont_256(ret, a, p, n0);
#endif
}
inline void from_mont_v256(vec256 ret, const vec256 a, const vec256 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	fromx_mont_256(ret, a, p, n0);// linked by asm
#else
	cc_from_mont_256(ret, a, p, n0);
#endif
}


inline void add_mod_v256(vec256 ret, const vec256 a, const vec256 b, const vec256 p)
{
#ifdef BLST_DEFAULT_ASM
	addx_mod_256(ret, a, b, p);// linked by asm
#else
	cc_add_mod_256(ret, a, b, p);
#endif

}

inline void sub_mod_v256(vec256 ret, const vec256 a, const vec256 b, const vec256 p)
{
#ifdef BLST_DEFAULT_ASM
	void sub_mod_256(ret, a, b, p);// linked by asm
#else
	cc_sub_mod_256(ret, a, b, p);
#endif
}

inline void mul_by_3_mod_v256(vec256 ret, const vec256 a, const vec256 p)
{
#ifdef BLST_DEFAULT_ASM
	mulx_by_3_mod_256(ret, a, p);
#else
	cc_mul_by_3_mod_256(ret, a, p);
#endif
}

inline void cneg_mod_v256(vec256 ret, const vec256 a, bool_t flag, const vec256 p)
{
#ifdef BLST_DEFAULT_ASM
	cneg_mod_256(ret, a, flag, p);
#else
	cc_cneg_mod_256(ret, a, flag, p);
#endif // BLST_DEFAULT_ASM

}
inline void lshift_mod_v256(vec256 ret, const vec256 a, size_t count, const vec256 p)
{
#ifdef BLST_DEFAULT_ASM
	lshift_mod_256(ret, a, count, p);
#else
	cc_lshift_mod_256(ret, a, count, p);
#endif
}

inline void rshift_mod_v256(vec256 ret, const vec256 a, size_t count, const vec256 p)
{
#ifdef BLST_DEFAULT_ASM
	rshift_mod_256(ret, a, count, p);
#else
	cc_rshift_mod_256(ret, a, count, p);
#endif
}

//
//bool_t eucl_inverse_mod_v256(vec256 ret, const vec256 a, const vec256 p, const vec256 one)
//{
//#ifdef BLST_DEFAULT_ASM
//    bool_t eucl_inverse_mod_256(vec256 ret, const vec256 a, const vec256 p, const vec256 one)
//#else
//	return cc_eucl_inverse_mod_256(ret, a, p, one);
//#endif
//
//}



inline limb_t check_mod_v256(const pow256 a, const vec256 p)
{
#ifdef BLST_DEFAULT_ASM
	return check_mod_256(a, p);
#else
	return cc_check_mod_256(a, p);
#endif
}

inline limb_t add_n_check_mod_v256(pow256 ret, const pow256 a, const pow256 b,
	const vec256 p)
{
#ifdef BLST_DEFAULT_ASM
	return add_n_check_mod_256(ret, a, b, p);
#else
	return cc_add_n_check_mod_256(ret, a, b, p);
#endif
}

// performs modular subtraction with zero-checking - 
// it computes (a - b) mod p and returns whether the result is non-zero.
inline limb_t sub_n_check_mod_v256(pow256 ret, const pow256 a, const pow256 b,
	const vec256 p)
{
#ifdef BLST_DEFAULT_ASM
	return sub_n_check_mod_256(ret, a, b, p);
#else
	return cc_sub_n_check_mod_256(ret, a, b, p);
#endif
}


#ifndef BLST_ASM
inline void vec_prefetch(const void* ptr, size_t len)
{
	cc_vec_prefetch(ptr, len);
}
#else
void vec_prefetch(const void* ptr, size_t len); // linked by asm
#endif

// 
inline void mul_mont_v384(vec384 ret, const vec384 a, const vec384 b,
	const vec384 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	mulx_mont_384(ret, a, b, p, n0);
#else
	cc_mul_mont_384(ret, a, b, p, n0);
#endif
}

inline void sqr_mont_v384(vec384 ret, const vec384 a, const vec384 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	sqrx_mont_384(ret, a, p, n0);
#else
	cc_sqr_mont_384(ret, a, p, n0);
#endif
}

inline void sqr_n_mul_mont_v384(vec384 ret, const vec384 a, size_t count,
	const vec384 p, limb_t n0, const vec384 b)
{
#ifdef BLST_DEFAULT_ASM
	sqrx_n_mul_mont_384(ret, a, count, p, n0, b);
#else
	cc_sqr_n_mul_mont_384(ret, a, count, p, n0, b);
#endif
}

inline void sqr_n_mul_mont_v383(vec384 ret, const vec384 a, size_t count,
	const vec384 p, limb_t n0, const vec384 b)
{
#ifdef BLST_DEFAULT_ASM
	sqrx_n_mul_mont_383(ret, a, count, p, n0, b);
#else
	cc_sqr_n_mul_mont_384(ret, a, count, p, n0, b);
#endif
}

//void mul_384(vec768 ret, const vec384 a, const vec384 b)
//{
//#ifdef BLST_DEFAULT_ASM
//	mulx_384(ret, a, b);
//#else
//	cc_mul_384(ret, a, b);
//#endif
//}
//void sqr_384(vec768 ret, const vec384 a);

inline void redc_mont_v384(vec384 ret, const vec768 a, const vec384 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	redcx_mont_384(ret, a, p, n0);
#else
	cc_redc_mont_384(ret, a, p, n0);
#endif
}

inline void from_mont_v384(vec384 ret, const vec384 a, const vec384 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	fromx_mont_384(ret, a, p, n0);
#else
	cc_from_mont_384(ret, a, p, n0);
#endif
}


inline limb_t sgn0_pty_mont_v384(const vec384 a, const vec384 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	return sgn0x_pty_mont_384(a, p, n0);
#else
	return cc_sgn0_pty_mont_384(a, p, n0);
#endif
}

inline limb_t sgn0_pty_mont_v384x(const vec384x a, const vec384 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	return sgn0x_pty_mont_384x(a, p, n0);
#else
	return cc_sgn0_pty_mont_384x(a, p, n0);
#endif
}
inline limb_t sgn0_pty_mod_v384(const vec384 a, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	return sgn0x_pty_mod_384(a, p);
#else
	return cc_sgn0_pty_mod_384(a, p);
#endif

}
inline limb_t sgn0_pty_mod_v384x(const vec384x a, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	return sgn0x_pty_mod_384x(a, p);
#else
	return cc_sgn0_pty_mod_384x(a, p);
#endif
}

inline void add_mod_v384(vec384 ret, const vec384 a, const vec384 b, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	addx_mod_384(ret, a, b, p);
#else
	cc_add_mod_384(ret, a, b, p);
#endif
}

inline void sub_mod_v384(vec384 ret, const vec384 a, const vec384 b, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	sub_mod_384(ret, a, b, p);
#else
	cc_sub_mod_n(ret, a, b, p, NLIMBS(384));
#endif
}
inline void mul_by_8_mod_v384(vec384 ret, const vec384 a, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	mul_by_8_mod_384(ret, a, p);	
#else
	cc_mul_by_8_mod_384(ret, a, p);
#endif

}
inline void mul_by_3_mod_v384(vec384 ret, const vec384 a, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	mulx_by_3_mod_384(ret, a, p);
#else
	cc_mul_by_3_mod_384(ret, a, p);
#endif
}

inline void cneg_mod_v384(vec384 ret, const vec384 a, bool_t flag, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	cneg_mod_384(ret, a, flag, p);
#else
	cc_cneg_mod_384(ret, a, flag, p);
#endif
}

void lshift_mod_v384(vec384 ret, const vec384 a, size_t count, const vec384 p);

inline void rshift_mod_v384(vec384 ret, const vec384 a, size_t count, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	rshift_mod_384(ret, a, count, p);
#else
	cc_rshift_mod_384(ret, a, count, p);
#endif
}
inline void div_by_2_mod_v384(vec384 ret, const vec384 a, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	div_by_2_mod_384(ret, a, p);
#else
	cc_div_by_2_mod_384(ret, a, p);
#endif
}
inline void ct_inverse_mod_v384(vec768 ret, const vec384 inp, const vec384 mod,
	const vec384 modx)
{
#ifdef BLST_DEFAULT_ASM
	ctx_inverse_mod_384(ret, inp, mod, modx);
#else
	cc_ct_inverse_mod_384(ret, inp, mod, modx);
#endif
}

inline void ct_inverse_mod_v256(vec512 ret, const vec256 inp, const vec256 mod,
	const vec256 modx)
{
#ifdef BLST_DEFAULT_ASM
	ctx_inverse_mod_256(ret, inp, mod, modx);
#else	
	cc_ct_inverse_mod_256(ret, inp, mod, modx);
#endif
}

inline bool_t ct_is_square_mod_v384(const vec384 inp, const vec384 mod)
{
#ifdef BLST_DEFAULT_ASM
	return ct_is_square_mod_384(inp, mod);
#else
	return cc_ct_is_square_mod_384(inp, mod);
#endif
}

//#if defined(__ADX__) /* e.g. -march=broadwell */ && !defined(__BLST_PORTABLE__)
//# define mul_mont_384x mulx_mont_384x
//# define sqr_mont_384x sqrx_mont_384x
//# define sqr_mont_382x sqrx_mont_382x
//# define mul_382x mulx_382x
//# define sqr_382x sqrx_382x
//#endif

inline void mul_mont_v384x(vec384x ret, const vec384x a, const vec384x b,
	const vec384 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	mulx_mont_384x(ret, a, b, p, n0);
#else
	cc_mul_mont_384x(ret, a, b, p, n0);
#endif
}
inline void sqr_mont_v384x(vec384x ret, const vec384x a, const vec384 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	sqrx_mont_384x(ret, a, p, n0);
#else
	cc_sqr_mont_384x(ret, a, p, n0);
#endif

}
inline void sqr_mont_v382x(vec384x ret, const vec384x a, const vec384 p, limb_t n0)
{
#ifdef BLST_DEFAULT_ASM
	sqrx_mont_382x(ret, a, p, n0);
#else
	cc_sqr_mont_382x(ret, a, p, n0);
#endif
}
//void mul_v382x(vec768 ret[2], const vec384x a, const vec384x b, const vec384 p)
//{
//#ifdef BLST_DEFAULT_ASM
//	mulx_382x(ret, a, b, p);
//#else
//	exit(1); // not implemented
//	//cc_mul_382x(ret, a, b, p);
//#endif
//}
//void sqr_v382x(vec768 ret[2], const vec384x a, const vec384 p)
//{
//	#ifdef BLST_DEFAULT_ASM
//	sqrx_382x(ret, a, p);
//#else
//	exit(1); // not implemented
//	//cc_sqr_382x(ret, a, p);
//}

inline void add_mod_v384x(vec384x ret, const vec384x a, const vec384x b,
	const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	add_mod_384x(ret, a, b, p);
#else
	cc_add_mod_384x(ret, a, b, p);
#endif
}
inline void sub_mod_v384x(vec384x ret, const vec384x a, const vec384x b,
	const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	sub_mod_384x(ret, a, b, p);
#else
	cc_sub_mod_384x(ret, a, b, p);
#endif
}
inline void mul_by_8_mod_v384x(vec384x ret, const vec384x a, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	mul_by_8_mod_384x(ret, a, p);
#else
	cc_mul_by_8_mod_384x(ret, a, p);
#endif
}
inline void mul_by_3_mod_v384x(vec384x ret, const vec384x a, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	mul_by_3_mod_384x(ret, a, p);
#else
	cc_mul_by_3_mod_384x(ret, a, p);
#endif
}
inline void mul_by_1_plus_i_mod_384x(vec384x ret, const vec384x a, const vec384 p)
{
#ifdef BLST_DEFAULT_ASM
	mul_by_1_plus_i_mod_384x(ret, a, p);
#else
	cc_mul_by_1_plus_i_mod_384x(ret, a, p);
#endif
}

HOST_DEVICE_FUNC inline void blst_sha256_block_data_order_v(unsigned int* v, const void* inp,
	size_t blocks)
{
#ifdef BLST_DEFAULT_ASM
	blst_sha256_block_data_order(v, inp, blocks);
#else
	cc_blst_sha256_block_data_order(v, inp, blocks);
#endif
}


HOST_DEVICE_FUNC inline void blst_sha256_hcopy_v(unsigned int dst[8], const unsigned int src[8])
{
#ifdef BLST_DEFAULT_ASM
	blst_sha256_hcopy(dst, src);
#else
	cc_blst_sha256_hcopy(dst, src);
#endif
}

HOST_DEVICE_FUNC inline void blst_sha256_emit_v(unsigned char md[32], const unsigned int h[8])
{
#ifdef BLST_DEFAULT_ASM
	blst_sha256_emit(md, h);
#else
	cc_blst_sha256_emit(md, h);
#endif
}

HOST_DEVICE_FUNC inline void blst_sha256_bcopy_v(void* dst_, const void* src_, size_t len)
{
	#ifdef BLST_DEFAULT_ASM
	blst_sha256_bcopy(dst_, src_, len);	
#else
	cc_blst_sha256_bcopy(dst_, src_, len);
#endif
}


//void add_mod_v384x384(vec768 ret, const vec768 a, const vec768 b,
//	const vec384 p)
//{
//#ifdef BLST_DEFAULT_ASM
//	add_mod_384x384(ret, a, b, p);
//#else
//	exit(1); // not implemented
//	//cc_add_mod_384x384(ret, a, b, p);
//#endif
//}
void sub_mod_384x384(vec768 ret, const vec768 a, const vec768 b,
	const vec384 p);

/*
 * C subroutines
 */
void exp_mont_384(vec384 out, const vec384 inp, const byte* pow,
	size_t pow_bits, const vec384 p, limb_t n0);
void exp_mont_384x(vec384x out, const vec384x inp, const byte* pow,
	size_t pow_bits, const vec384 p, limb_t n0);
void div_by_zz(limb_t val[]);
void div_by_z(limb_t val[]);


inline bool_t is_bit_set(const byte* v, size_t i)
{
	bool_t ret = (v[i / 8] >> (i % 8)) & 1;
	BLST_launder(ret);
	return ret;
}

inline bool_t byte_is_zero(unsigned char c)
{
	limb_t ret = ((limb_t)(c)-1) >> (BLST_LIMB_SIZE - 1);
	BLST_launder(ret);
	return ret;
}

inline bool_t bytes_are_zero(const unsigned char* a, size_t num)
{
	unsigned char acc;
	size_t i;

	for (acc = 0, i = 0; i < num; i++)
		acc |= a[i];

	return byte_is_zero(acc);
}

inline void vec_cswap(void* restrict a, void* restrict b, size_t num,
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



/* ret = bit ? a : b */
//void vec_select_32(void* ret, const void* a, const void* b, bool_t sel_a);
//void vec_select_48(void* ret, const void* a, const void* b, bool_t sel_a);
//void vec_select_96(void* ret, const void* a, const void* b, bool_t sel_a);
//void vec_select_144(void* ret, const void* a, const void* b, bool_t sel_a);
//void vec_select_192(void* ret, const void* a, const void* b, bool_t sel_a);
//void vec_select_288(void* ret, const void* a, const void* b, bool_t sel_a);
//
//

inline void cneg_mod_v384x(vec384x ret, const vec384x a, bool_t flag,
	const vec384 p)
{
	cneg_mod_v384(ret[0], a[0], flag, p);
	cneg_mod_v384(ret[1], a[1], flag, p);
}


/*
 * Some compilers get arguably overzealous(*) when passing pointer to
 * multi-dimensional array [such as vec384x] as 'const' argument.
 * General direction seems to be to legitimize such constification,
 * so it's argued that suppressing the warning is appropriate.
 *
 * (*)  http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1923.htm
 */
#if defined(__INTEL_COMPILER)
# pragma warning(disable:167)
# pragma warning(disable:556)
#elif defined(__GNUC__) && !defined(__clang__) && (__STDC_VERSION__-0) < 202311
# pragma GCC diagnostic ignored "-Wpedantic"
#elif defined(_MSC_VER)
# pragma warning(disable: 4127 4189)
#endif

#if !defined(__wasm__) && __STDC_HOSTED__-0 != 0
# include <stdlib.h>
#endif

#if defined(__GNUC__)
# ifndef alloca
#  define alloca(s) __builtin_alloca(s)
# endif
#elif defined(__sun)
# include <alloca.h>
#elif defined(_WIN32)
# include <malloc.h>
# ifndef alloca
#  define alloca(s) _alloca(s)
# endif
#endif


#endif /* __BLS12_381_ASM_VECT_H__ */