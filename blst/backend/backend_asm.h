/*
 * Copyright Supranational LLC
 * Copyright 2026 Peter Rindal
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 *
 * This CatBLST header is locally composed. Its assembly ABI declarations
 * retain upstream BLST src/vect.h lineage.
 */

#ifndef BLST_BACKEND_ASM_H
#define BLST_BACKEND_ASM_H

#include "blst/config.h"
#include "blst/common/defines.h"
#include <stddef.h>

#ifdef BLST_ASM

/*
 * GNU-style portable x86 assembly includes both providers. Non-portable x86
 * includes only the configured provider. The Win64 build compiles both public
 * provider objects. AArch64 has only the generic provider names.
 */
#if defined(__x86_64__) || defined(__x86_64) || \
    (defined(_M_X64) && !defined(_M_ARM64EC))
# if defined(_MSC_VER) || defined(BLST_ADX) || defined(BLST_PORTABLE)
#  define BLST_ASM_HAS_X_PROVIDER 1
# endif
# if defined(_MSC_VER) || !defined(BLST_ADX) || defined(BLST_PORTABLE)
#  define BLST_ASM_HAS_Q_PROVIDER 1
# endif
#elif defined(__aarch64__) || defined(__aarch64) || defined(_M_ARM64) || \
      defined(_M_ARM64EC)
# define BLST_ASM_HAS_Q_PROVIDER 1
#endif

BLST_EXTERN_C_BEGIN_

/* These annotations document raw assembly object shapes; the ABI is pointers. */
#define asm_pow256  /* 32 bytes */ void*
#define asm_Fr      /* 4 limbs  */ void*
#define asm_UFr     /* 8 limbs  */ void*
#define asm_Fp      /* 6 limbs  */ void*
#define asm_UFp     /* 12 limbs */ void*
#define asm_Fp2     /* 12 limbs */ void*
#define asm_UFp2    /* 24 limbs */ void*

/* ret = sel_a ? a : b */
void vec_select_32(void* ret, const void* a, const void* b, blst::bool_t sel_a);
void vec_select_48(void* ret, const void* a, const void* b, blst::bool_t sel_a);
void vec_select_96(void* ret, const void* a, const void* b, blst::bool_t sel_a);
void vec_select_144(void* ret, const void* a, const void* b, blst::bool_t sel_a);
void vec_select_192(void* ret, const void* a, const void* b, blst::bool_t sel_a);
void vec_select_288(void* ret, const void* a, const void* b, blst::bool_t sel_a);

blst::limb_t add_n_check_mod_256(asm_pow256 ret, const asm_pow256 a,
    const asm_pow256 b, const asm_Fr p);
blst::limb_t sub_n_check_mod_256(asm_pow256 ret, const asm_pow256 a,
    const asm_pow256 b, const asm_Fr p);
void add_mod_256(asm_Fr ret, const asm_Fr a, const asm_Fr b, const asm_Fr p);
void sub_mod_256(asm_Fr ret, const asm_Fr a, const asm_Fr b, const asm_Fr p);
void mul_by_3_mod_256(asm_Fr ret, const asm_Fr a, const asm_Fr p);
void cneg_mod_256(asm_Fr ret, const asm_Fr a, blst::bool_t flag, const asm_Fr p);
void lshift_mod_256(asm_Fr ret, const asm_Fr a, size_t count, const asm_Fr p);
void rshift_mod_256(asm_Fr ret, const asm_Fr a, size_t count, const asm_Fr p);
blst::limb_t check_mod_256(const asm_pow256 a, const asm_Fr p);

/* Raw assembly requires a positive byte count divisible by 16. */
blst::bool_t vec_is_zero_16x(const void* a, size_t num);
blst::bool_t vec_is_equal_16x(const void* a, const void* b, size_t num);

/*
 * Single-width modular operations require reduced operands in [0,p). The
 * modulus is odd. Exact output aliasing with either input is supported;
 * partial overlap is not. Raw shift counts must be positive and are public.
 * Every operation uses data-independent control flow for fixed public sizes.
 */

void add_mod_384(asm_Fp ret, const asm_Fp a, const asm_Fp b, const asm_Fp p);
void sub_mod_384(asm_Fp ret, const asm_Fp a, const asm_Fp b, const asm_Fp p);
void mul_by_8_mod_384(asm_Fp ret, const asm_Fp a, const asm_Fp p);
void mul_by_3_mod_384(asm_Fp ret, const asm_Fp a, const asm_Fp p);
void cneg_mod_384(asm_Fp ret, const asm_Fp a, blst::bool_t flag, const asm_Fp p);
void rshift_mod_384(asm_Fp ret, const asm_Fp a, size_t count, const asm_Fp p);
void lshift_mod_384(asm_Fp ret, const asm_Fp a, size_t count, const asm_Fp p);
void div_by_2_mod_384(asm_Fp ret, const asm_Fp a, const asm_Fp p);
blst::limb_t sgn0_pty_mod_384(const asm_Fp a, const asm_Fp p);
blst::limb_t sgn0_pty_mod_384x(const asm_Fp2 a, const asm_Fp p);

/* Raw prefetch requires len > 0. The selected C++ boundary handles zero. */
void vec_prefetch_asm(const void* ptr, size_t len);

void add_mod_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp2 b,
    const asm_Fp p);
void sub_mod_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp2 b,
    const asm_Fp p);
void mul_by_8_mod_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp p);
void mul_by_3_mod_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp p);
void mul_by_1_plus_i_mod_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp p);

/*
 * Fixed-schedule almost-inverse contract. Let w be the input width and
 * R=2^w. mod is an odd positive w-bit integer, 0<=inp<mod, and modx is mod
 * shifted left until its most-significant bit is bit w-1. All four arrays are
 * disjoint. ret has width 2w and satisfies:
 *
 *   inp == 0: ret == 0;
 *   gcd(inp,mod) == 1: ret*inp == R^2 (mod mod);
 *   0 <= ret < modx*R.
 *
 * The routines do not publish a conventional-width inverse. Their field
 * callers normalize this bounded result. Control flow and memory addresses
 * are independent of inp; provider selection and loop counts are public.
 */
void ct_inverse_mod_256(asm_UFr ret, const asm_Fr inp, const asm_Fr mod,
    const asm_Fr modx);

/*
 * Return one exactly for zero or a quadratic residue inp modulo the odd prime
 * mod. The input is reduced: 0<=inp<mod. The fixed schedule has no
 * input-dependent control flow or memory address.
 */
blst::bool_t ct_is_square_mod_384(const asm_Fp inp, const asm_Fp mod);

/* Raw SHA assembly requires positive block/byte counts. */
void blst_sha256_block_data_order(unsigned int* v, const void* inp,
    size_t blocks);
void blst_sha256_hcopy(unsigned int dst[8], const unsigned int src[8]);
void blst_sha256_emit(unsigned char md[32], const unsigned int h[8]);
void blst_sha256_bcopy(void* dst, const void* src, size_t len);

/*
 * Estimate one quotient limb for (dividend_top * 2^64) / divisor.
 * The two-limb divisor is normalized: bit 127 is set. dividend_top is at
 * most divisor. Equality saturates an overflowing quotient to UINT64_MAX.
 * For the long-division correction contract, dividend_top is strictly less
 * than divisor and the returned estimate is exact or one too high.
 */
blst::limb_t div_3_limbs(const blst::limb_t dividend_top[2],
    blst::limb_t d_lo, blst::limb_t d_hi);
/*
 * div_rem initially contains a 192-bit dividend. quotient is exact or one
 * too high. The routine writes the 128-bit remainder followed by the corrected
 * quotient. Exact div_rem/dividend storage reuse is required by the contract.
 */
blst::limb_t quot_rem_128(blst::limb_t* quot_rem,
    const blst::limb_t* divisor, blst::limb_t quotient);
/*
 * div_rem initially contains a 128-bit dividend. quotient must be exact. The
 * routine writes the 64-bit remainder followed by the unchanged quotient.
 */
blst::limb_t quot_rem_64(blst::limb_t* quot_rem,
    const blst::limb_t* divisor, blst::limb_t quotient);

/*
 * Double-width operations reduce modulo p * 2^384. Each operand must be below
 * that boundary. Exact output aliasing with either input is supported;
 * partial overlap is not.
 */
void add_mod_384x384(asm_UFp ret, const asm_UFp a, const asm_UFp b,
    const asm_Fp p);
void sub_mod_384x384(asm_UFp ret, const asm_UFp a, const asm_UFp b,
    const asm_Fp p);

/*
 * Montgomery arithmetic uses R = 2^256 or R = 2^384 and
 * n0 = -p^-1 mod 2^64. The modulus is odd and occupies the declared width.
 * Results from mul, sqr, redc, and from-Montgomery operations are in [0,p).
 * Exact output aliasing is supported where stated; partial overlap is not.
 * Control flow is independent of field values. Provider dispatch and repeated
 * square counts are public control inputs.
 *
 * Sparse 256-bit multiplication additionally requires
 * p <= 2^256 - 2^192 - 1. Its first operand is in [0,p), while its second
 * operand may occupy all of [0,R). Sparse squaring therefore requires its
 * single operand in [0,p). Multiplication supports ret == a or ret == b;
 * squaring supports ret == a.
 * The 256-bit REDC input T satisfies T < pR; ret may equal T's low-width
 * prefix. From-Montgomery accepts a in [0,R) and supports ret == a.
 *
 * The exposed 384-bit mul and sqr boundaries accept reduced operands in
 * [0,p). They support exact output aliasing with either input. The REDC input
 * T satisfies T < pR and permits reuse of its low-width prefix. From-
 * Montgomery accepts a in [0,R) and supports ret == a. The specialized
 * sqr_n_mul_mont_383 boundary requires reduced a and b plus public count > 0;
 * ret may equal a but must not equal b.
 *
 * Fp2 inputs are two contiguous reduced Fp components in real-imaginary
 * order. Fp2 Montgomery mul and sqr support exact output aliasing. The 382x
 * operations return two contiguous double-width components, each below pR;
 * their output must not overlap either input.
 */
#ifdef BLST_ASM_HAS_X_PROVIDER
void mulx_mont_sparse_256(asm_Fr ret, const asm_Fr a, const asm_Fr b,
    const asm_Fr p, blst::limb_t n0);
void sqrx_mont_sparse_256(asm_Fr ret, const asm_Fr a, const asm_Fr p,
    blst::limb_t n0);
void redcx_mont_256(asm_Fr ret, const asm_UFr a, const asm_Fr p,
    blst::limb_t n0);
void fromx_mont_256(asm_Fr ret, const asm_Fr a, const asm_Fr p,
    blst::limb_t n0);

void mulx_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp b,
    const asm_Fp p, blst::limb_t n0);
void sqrx_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp p,
    blst::limb_t n0);
void sqrx_n_mul_mont_383(asm_Fp ret, const asm_Fp a, size_t count,
    const asm_Fp p, blst::limb_t n0, const asm_Fp b);
void redcx_mont_384(asm_Fp ret, const asm_UFp a, const asm_Fp p,
    blst::limb_t n0);
void fromx_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp p,
    blst::limb_t n0);
blst::limb_t sgn0x_pty_mont_384(const asm_Fp a, const asm_Fp p,
    blst::limb_t n0);
blst::limb_t sgn0x_pty_mont_384x(const asm_Fp2 a, const asm_Fp p,
    blst::limb_t n0);
/* Same fixed-schedule almost-inverse contract as ct_inverse_mod_256. */
void ctx_inverse_mod_384(asm_UFp ret, const asm_Fp inp, const asm_Fp mod,
    const asm_Fp modx);
void mulx_mont_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp2 b,
    const asm_Fp p, blst::limb_t n0);
void sqrx_mont_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp p,
    blst::limb_t n0);
void mulx_382x(asm_UFp2 ret, const asm_Fp2 a, const asm_Fp2 b,
    const asm_Fp p);
void sqrx_382x(asm_UFp2 ret, const asm_Fp2 a, const asm_Fp p);
#endif

#ifdef BLST_ASM_HAS_Q_PROVIDER
void mul_mont_sparse_256(asm_Fr ret, const asm_Fr a, const asm_Fr b,
    const asm_Fr p, blst::limb_t n0);
void sqr_mont_sparse_256(asm_Fr ret, const asm_Fr a, const asm_Fr p,
    blst::limb_t n0);
void redc_mont_256(asm_Fr ret, const asm_UFr a, const asm_Fr p,
    blst::limb_t n0);
void from_mont_256(asm_Fr ret, const asm_Fr a, const asm_Fr p,
    blst::limb_t n0);

void mul_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp b,
    const asm_Fp p, blst::limb_t n0);
void sqr_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp p,
    blst::limb_t n0);
void redc_mont_384(asm_Fp ret, const asm_UFp a, const asm_Fp p,
    blst::limb_t n0);
void from_mont_384(asm_Fp ret, const asm_Fp a, const asm_Fp p,
    blst::limb_t n0);
void sqr_n_mul_mont_383(asm_Fp ret, const asm_Fp a, size_t count,
    const asm_Fp p, blst::limb_t n0, const asm_Fp b);
blst::limb_t sgn0_pty_mont_384(const asm_Fp a, const asm_Fp p,
    blst::limb_t n0);
blst::limb_t sgn0_pty_mont_384x(const asm_Fp2 a, const asm_Fp p,
    blst::limb_t n0);
/* Same fixed-schedule almost-inverse contract as ct_inverse_mod_256. */
void ct_inverse_mod_384(asm_UFp ret, const asm_Fp inp, const asm_Fp mod,
    const asm_Fp modx);
void mul_mont_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp2 b,
    const asm_Fp p, blst::limb_t n0);
void sqr_mont_384x(asm_Fp2 ret, const asm_Fp2 a, const asm_Fp p,
    blst::limb_t n0);
void mul_382x(asm_UFp2 ret, const asm_Fp2 a, const asm_Fp2 b,
    const asm_Fp p);
void sqr_382x(asm_UFp2 ret, const asm_Fp2 a, const asm_Fp p);
#endif

#if (defined(__x86_64__) || defined(__x86_64) || defined(_M_X64)) && \
    defined(__SHA__)
void blst_sha256_block_data_order_shaext(unsigned int* v, const void* inp,
    size_t blocks);
#endif

#undef asm_pow256
#undef asm_Fr
#undef asm_UFr
#undef asm_Fp
#undef asm_UFp
#undef asm_Fp2
#undef asm_UFp2

BLST_EXTERN_C_END_

#undef BLST_ASM_HAS_X_PROVIDER
#undef BLST_ASM_HAS_Q_PROVIDER

#endif /* BLST_ASM */
#endif /* BLST_BACKEND_ASM_H */
