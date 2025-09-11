/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_NO_ASM_H__
#define __BLS12_NO_ASM_H__

#include "bytes.h"
#include "config.h"
#include <assert.h>

typedef struct { limb_t lo, hi; } wide_limb_t;
//static_assert(((wide_limb_t*)0));

#if defined(_MSC_VER)
#include <malloc.h>
#define BLST_VLA_NEW(T, name, n) T* name = (T*)_alloca((n) * sizeof(T))
#define BLST_VLA_FREE(name)      ((void)0)


#if defined(_M_X64) || defined(_M_IX86)
#include <immintrin.h>
#endif

#elif defined(__GNUC__) || defined(__clang__)
#include <alloca.h>
#define BLST_VLA_NEW(T, name, n) T name[n]
#define BLST_VLA_FREE(name)      ((void)0)
#else
#include <stdlib.h>
#define BLST_VLA_NEW(T, name, n) T* name = (T*)malloc((n) * sizeof(T))
#define BLST_VLA_FREE(name)      free(name)
#endif

#if defined(__clang__)
# pragma GCC diagnostic ignored "-Wstatic-in-inline"
#endif

#if !defined(__clang__) && !defined(__builtin_assume)
# if defined(__GNUC__) && __GNUC__>=5
#  define __builtin_assume(condition) if (!(condition)) __builtin_unreachable()
# elif defined(_MSC_VER)
#  define __builtin_assume(condition) __assume(condition)
# else
#  define __builtin_assume(condition) (void)(condition)
# endif
#endif


#if BLST_LIMB_SIZE == 64 && defined(BLST_ADX)
// 64 bit with wide operations
#define cc_mul_wide_limb_limb(ret, a, b)\
    (ret)->lo = _mulx_u64(a, b, &(ret)->hi)

#define cc_add_wide_wide_limb(ret, a, b)\
    (ret)->hi = a.hi + _addcarry_u64(0, (a).lo, b, (unsigned long long*) &(ret)->lo)

#define cc_add_wide_limb_limb(ret, a, b)\
    (ret)->hi = _addcarry_u64(0, a, b, (unsigned long long*) &(ret)->lo)

#define cc_add_carry_wide_limb_limb(ret, a, b, carry) \
    do{ assert(carry <= 1); __assume(carry <=1); \
    (ret)->hi = _addcarry_u64(carry, a, b, (unsigned long long*) &(ret)->lo); } while(0)

#define cc_sub_borrow_wide_limb_limb(ret, a, b, borrow) \
    do{ assert(borrow <= 1); __assume(borrow <=1); \
    (ret)->hi = _subborrow_u64(borrow, a, b, (unsigned long long*) &(ret)->lo); } while(0)


#elif BLST_LIMB_SIZE == 32 && (defined(_M_X64) || defined(__x86_64__)) && 0
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
    do{ assert(carry <= 1); __assume(carry <=1); \
    *(unsigned long long*)(ret) = (a) + (unsigned long long)(b) + carry; } while(0)

#define cc_sub_borrow_wide_limb_limb(ret, a, b, borrow) \
    do{ assert(borrow <= 1); __assume(borrow <=1); \
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

#define cc_add_wide_wide_limb(ret, a, b)\
    do{ \
    wide_limb_t aa = a;\
    (ret)->lo = aa.lo + (b);\
    (ret)->hi = aa.hi + (((ret)->lo < aa.lo) ? 1 : 0);}while(0)

#define cc_add_wide_limb_limb(ret, a, b)\
    do{ \
    limb_t aa = a;\
    (ret)->lo = (a) + (b);\
    (ret)->hi = ((ret)->lo < aa) ? 1 : 0;\
    }while(0)

// wide subtraction *ret = a - b
// with wide a and limb b
#define cc_sub_wide_wide_limb(ret, a, b)\
do {                                                                           \
    wide_limb_t _aa = (a);                                                     \
    limb_t      _bb = (limb_t)(b);                                            \
    limb_t      _old = _aa.lo;                                                 \
    (ret)->lo = (limb_t)(_old - _bb);                                          \
    limb_t _borrow = (limb_t)(_old < _bb);    /* correct borrow test */        \
    (ret)->hi = (limb_t)(_aa.hi - _borrow);                                    \
} while (0)

// wide subtraction *ret = a - b
// with limb a and b
#define cc_sub_wide_limb_limb(ret, a, b)\
do {                                                                           \
    limb_t _aa = (limb_t)(a);                                                  \
    limb_t _bb = (limb_t)(b);                                                  \
    (ret)->lo = (limb_t)(_aa - _bb);                                           \
    limb_t _borrow = (limb_t)(_aa < _bb);                                      \
    (ret)->hi = (limb_t)(-(limb_t)_borrow); /* 0 or 0xFFFFFFFF */              \
} while (0)

// ret = a + b + c
#define cc_add_carry_wide_limb_limb(ret, a, b, c)\
    cc_add_wide_limb_limb(ret, a, b);\
    cc_add_wide_wide_limb(ret, *ret, c)

// ret = a - b - c
#define cc_sub_borrow_wide_limb_limb(ret, a, b, c)\
    cc_sub_wide_limb_limb(ret, a, b);\
    cc_sub_wide_wide_limb(ret, *ret, c)
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


HOST_DEVICE_FUNC inline void cc_mul_mont_n(limb_t ret[], const limb_t a[], const limb_t b[],
                       const limb_t p[], limb_t n0, size_t n)
{
    __builtin_assume(n != 0 && n%2 == 0);
    limb_t mask, borrow, mx, hi, carry;
	BLST_VLA_NEW(limb_t, tmp, n + 1);
    size_t i, j;

    for (mx=b[0], hi=0, i=0; i<n; i++) {
        wide_limb_t  wide;

		cc_mul_add_wide_limb_limb(&wide, mx, a[i], hi);
		tmp[i] = (limb_t)wide.lo;
		hi = (limb_t)(wide.hi);
    }
    mx = n0*tmp[0];
    tmp[i] = hi;

    for (carry=0, j=0; ; ) {
		wide_limb_t  wide;
        // wide = mx * p + tmp
		cc_mul_add_wide_limb_limb(&wide, mx, p[0], tmp[0]);
		hi = wide.hi;

        for (i=1; i<n; i++) {
            // wide = mx * p + hi + tmp;
			cc_mul_add_add_wide_limb_limb(&wide, mx, p[i], hi, tmp[i]); 
			tmp[i - 1] = wide.lo;
			hi = wide.hi;
        }
        // wide = hi + carry + tmp
        cc_add_carry_wide_limb_limb(&wide, hi, tmp[i], carry);
		tmp[i - 1] = wide.lo;
		carry = wide.hi;


        if (++j==n)
            break;

        for (mx=b[j], hi=0, i=0; i<n; i++) {
            // wide = mx * a + hi + tmp
            cc_mul_add_add_wide_limb_limb(&wide, mx, a[i], hi, tmp[i]);
			tmp[i] = wide.lo;
			hi = wide.hi;
        }
        mx = n0*tmp[0];
		cc_add_wide_limb_limb(&wide, hi, carry);
		tmp[i] = wide.lo;
		carry = wide.hi;
    }

    for (borrow=0, i=0; i<n; i++) {
        // wide = tmp - (p + borrow)
        wide_limb_t wide;
        cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
        ret[i] = wide.lo;
        borrow = wide.hi & 1;
    }

    mask = carry - borrow;
    BLST_launder(mask);

    for(i=0; i<n; i++)
        ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);

    BLST_VLA_FREE(tmp);
}

#define MUL_MONT_IMPL(fn1, fn2, bits) \
HOST_DEVICE_FUNC inline void fn1(vec##bits ret, const vec##bits a, \
                            const vec##bits b, const vec##bits p, limb_t n0) \
{   cc_mul_mont_n(ret, a, b, p, n0, NLIMBS(bits));   } \
\
HOST_DEVICE_FUNC inline void fn2(vec##bits ret, const vec##bits a, \
                            const vec##bits p, limb_t n0) \
{   cc_mul_mont_n(ret, a, a, p, n0, NLIMBS(bits));   }

/*
 * 256-bit subroutines can handle arbitrary modulus, even non-"sparse",
 * but we have to harmonize the naming with assembly.
 */
MUL_MONT_IMPL(cc_mul_mont_sparse_256, cc_sqr_mont_sparse_256, 256)
MUL_MONT_IMPL(cc_mul_mont_384, cc_sqr_mont_384, 384)

HOST_DEVICE_FUNC inline void cc_add_mod_n(limb_t ret[], const limb_t a[], const limb_t b[],
                      const limb_t p[], size_t n)
{
    __builtin_assume(n != 0);
    limb_t mask, carry, borrow;
    BLST_VLA_NEW(limb_t, tmp, n);
    size_t i;

    for (carry=0, i=0; i<n; i++) {
        wide_limb_t wide;
        cc_add_carry_wide_limb_limb(&wide, a[i], b[i], carry);
        tmp[i] = wide.lo;
        carry = wide.hi;
    }

    for (borrow=0, i=0; i<n; i++) {
        wide_limb_t wide;
        cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
        ret[i] = wide.lo;
        borrow = wide.hi & 1;
    }

    mask = carry - borrow;
    BLST_launder(mask);

    for(i=0; i<n; i++)
        ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);

    BLST_VLA_FREE(tmp);
}

#define ADD_MOD_IMPL(bits) \
HOST_DEVICE_FUNC inline void cc_add_mod_##bits(vec##bits ret, const vec##bits a, \
                           const vec##bits b, const vec##bits p) \
{   cc_add_mod_n(ret, a, b, p, NLIMBS(bits));   }

ADD_MOD_IMPL(256)
ADD_MOD_IMPL(384)

inline void cc_sub_mod_n(limb_t ret[], const limb_t a[], const limb_t b[],
                      const limb_t p[], size_t n)
{
    __builtin_assume(n != 0);
    limb_t mask, carry, borrow;
    size_t i;

    for (borrow=0, i=0; i<n; i++) {
        wide_limb_t wide;
        cc_sub_borrow_wide_limb_limb(&wide, a[i], b[i], borrow);
        ret[i] = wide.lo;
        borrow = wide.hi & 1;
    }

    mask = 0 - borrow;
    BLST_launder(mask);

    for (carry=0, i=0; i<n; i++) {
        wide_limb_t wide;
        cc_add_carry_wide_limb_limb(&wide, ret[i], p[i] & mask, carry);
        ret[i] = wide.lo;
        carry = wide.hi;
    }
}

#define SUB_MOD_IMPL(bits) \
HOST_DEVICE_FUNC inline void cc_sub_mod_##bits(vec##bits ret, const vec##bits a, \
                           const vec##bits b, const vec##bits p) \
{   cc_sub_mod_n(ret, a, b, p, NLIMBS(bits));   }

SUB_MOD_IMPL(256)
SUB_MOD_IMPL(384)

inline void cc_mul_by_3_mod_n(limb_t ret[], const limb_t a[], const limb_t p[],
                           size_t n)
{
    __builtin_assume(n != 0);
    limb_t mask, carry, borrow;
    BLST_VLA_NEW(limb_t, tmp, n);
    BLST_VLA_NEW(limb_t, two_a, n);
    
    size_t i;

    for (carry=0, i=0; i<n; i++) {
        limb_t a_i = a[i];
        tmp[i] = a_i<<1 | carry;
        carry = a_i>>(BLST_LIMB_SIZE-1);
    }

    for (borrow=0, i=0; i<n; i++) {
        wide_limb_t wide;
        cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
        two_a[i] = wide.lo;
        borrow = wide.hi & 1;
    }

    mask = carry - borrow;
    BLST_launder(mask);

    for(i=0; i<n; i++)
        two_a[i] = (two_a[i] & ~mask) | (tmp[i] & mask);

    for (carry=0, i=0; i<n; i++) {
        wide_limb_t wide;
        cc_add_carry_wide_limb_limb(&wide, a[i], two_a[i], carry);
        tmp[i] = wide.lo;
        carry = wide.hi;
    }

    for (borrow=0, i=0; i<n; i++) {
        wide_limb_t wide;
        cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
        ret[i] = wide.lo;
        borrow = wide.hi & 1;
    }

    mask = carry - borrow;
    BLST_launder(mask);

    for(i=0; i<n; i++)
        ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);


    BLST_VLA_FREE(tmp);
    BLST_VLA_FREE(two_a);
}

#define MUL_BY_3_MOD_IMPL(bits) \
HOST_DEVICE_FUNC inline void cc_mul_by_3_mod_##bits(vec##bits ret, const vec##bits a, \
                                const vec##bits p) \
{   cc_mul_by_3_mod_n(ret, a, p, NLIMBS(bits));   }

MUL_BY_3_MOD_IMPL(256)
MUL_BY_3_MOD_IMPL(384)

inline void cc_lshift_mod_n(limb_t ret[], const limb_t a[], size_t count,
                         const limb_t p[], size_t n)
{
    __builtin_assume(count != 0);
    __builtin_assume(n != 0);
    limb_t mask, carry, borrow;
    BLST_VLA_NEW(limb_t, tmp, n);
    size_t i;

    while (count--) {
        for (carry=0, i=0; i<n; i++) {
            limb_t a_i = a[i];
            tmp[i] = a_i<<1 | carry;
            carry = a_i>>(BLST_LIMB_SIZE-1);
        }

        for (borrow=0, i=0; i<n; i++) {
            wide_limb_t wide;
            cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
            ret[i] = wide.lo;
            borrow = wide.hi & 1;
        }

        mask = carry - borrow;
        BLST_launder(mask);

        for(i=0; i<n; i++)
            ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);

        a = ret;
    }

    BLST_VLA_FREE(tmp);
}

#define LSHIFT_MOD_IMPL(bits) \
HOST_DEVICE_FUNC inline void cc_lshift_mod_##bits(vec##bits ret, const vec##bits a, size_t count, \
                              const vec##bits p) \
{   cc_lshift_mod_n(ret, a, count, p, NLIMBS(bits));   }

LSHIFT_MOD_IMPL(256)
LSHIFT_MOD_IMPL(384)

inline void cc_cneg_mod_n(limb_t ret[], const limb_t a[], bool_t flag,
                       const limb_t p[], size_t n)
{
    __builtin_assume(n != 0);
    limb_t borrow, mask;
    BLST_VLA_NEW(limb_t, tmp, n);
    size_t i;

    for (borrow=0, i=0; i<n; i++) {
        wide_limb_t wide;
        cc_sub_borrow_wide_limb_limb(&wide, p[i], a[i], borrow);
        tmp[i] = wide.lo;
        borrow = wide.hi & 1;
    }

    flag &= vec_is_zero(a, n * sizeof(limb_t)) ^ 1;
    mask = (limb_t)0 - flag;

    for(i=0; i<n; i++)
        ret[i] = (a[i] & ~mask) | (tmp[i] & mask);
    BLST_VLA_FREE(tmp);
}

#define CNEG_MOD_IMPL(bits) \
HOST_DEVICE_FUNC inline void cc_cneg_mod_##bits(vec##bits ret, const vec##bits a, bool_t flag, \
                            const vec##bits p) \
{   cc_cneg_mod_n(ret, a, flag, p, NLIMBS(bits));   }

CNEG_MOD_IMPL(256)
CNEG_MOD_IMPL(384)

inline limb_t cc_check_mod_n(const byte a[], const limb_t p[], size_t n)
{
    __builtin_assume(n != 0);
    limb_t borrow, ai, acc;
    size_t i, j;

    for (acc=borrow=0, i=0; i<n; i++) {
        for (ai=0, j=0; j<8*sizeof(limb_t); j+=8)
            ai |= (limb_t)(*a++) << j;
        acc |= ai;
        wide_limb_t wide;
        cc_sub_borrow_wide_limb_limb(&wide, ai, p[i], borrow);
        borrow = wide.hi & 1;
    }

    return borrow & (is_zero(acc) ^ 1);
}

#define CHECK_MOD_IMPL(bits) \
HOST_DEVICE_FUNC inline limb_t cc_check_mod_##bits(const pow##bits a, const vec##bits p) \
{   return cc_check_mod_n(a, p, NLIMBS(bits));   }

CHECK_MOD_IMPL(256)

inline limb_t cc_add_n_check_mod_n(byte ret[], const byte a[], const byte b[],
                                            const limb_t p[], size_t n)
{
    __builtin_assume(n != 0);
    limb_t zero;
    BLST_VLA_NEW(limb_t, ret_, n);
    BLST_VLA_NEW(limb_t, a_, n);
    BLST_VLA_NEW(limb_t, b_, n);

    limbs_from_le_bytes(a_, a, n * sizeof(limb_t));
    limbs_from_le_bytes(b_, b, n * sizeof(limb_t));

    cc_add_mod_n(ret_, a_, b_, p, n);
    zero = vec_is_zero(ret_, n * sizeof(limb_t));

    le_bytes_from_limbs(ret, ret_, n * sizeof(limb_t));

	BLST_VLA_FREE(ret_);
	BLST_VLA_FREE(a_);
	BLST_VLA_FREE(b_);

    return zero^1;
}

#define ADD_N_CHECK_MOD_IMPL(bits) \
HOST_DEVICE_FUNC inline limb_t cc_add_n_check_mod_##bits(pow##bits ret, const pow##bits a, \
                                     const pow##bits b, const vec##bits p) \
{   return cc_add_n_check_mod_n(ret, a, b, p, NLIMBS(bits));   }

ADD_N_CHECK_MOD_IMPL(256)

inline limb_t cc_sub_n_check_mod_n(byte ret[], const byte a[], const byte b[],
                                            const limb_t p[], size_t n)
{
    __builtin_assume(n != 0);
    limb_t zero;
    BLST_VLA_NEW(limb_t, ret_, n);
    BLST_VLA_NEW(limb_t, a_, n);
    BLST_VLA_NEW(limb_t, b_, n);

    limbs_from_le_bytes(a_, a, n * sizeof(limb_t));
    limbs_from_le_bytes(b_, b, n * sizeof(limb_t));

    cc_sub_mod_n(ret_, a_, b_, p, n);
    zero = vec_is_zero(ret_, n * sizeof(limb_t));

    le_bytes_from_limbs(ret, ret_, n * sizeof(limb_t));

    BLST_VLA_FREE(ret_);
    BLST_VLA_FREE(a_);
    BLST_VLA_FREE(b_);

    return zero^1;
}

#define SUB_N_CHECK_MOD_IMPL(bits) \
HOST_DEVICE_FUNC inline limb_t cc_sub_n_check_mod_##bits(pow##bits ret, const pow##bits a, \
                                     const pow##bits b, const vec##bits p) \
{   return cc_sub_n_check_mod_n(ret, a, b, p, NLIMBS(bits));   }

SUB_N_CHECK_MOD_IMPL(256)

inline void cc_from_mont_n(limb_t ret[], const limb_t a[],
                        const limb_t p[], limb_t n0, size_t n)
{
    __builtin_assume(n != 0 && n%2 == 0);
    wide_limb_t wide;
    limb_t mask, borrow, mx, hi;
    size_t i, j;
    BLST_VLA_NEW(limb_t, tmp, n);

    for (j=0; j<n; j++) {
        mx = n0*a[0];
        cc_mul_add_wide_limb_limb(&wide, mx, p[0], a[0]);
        hi = wide.hi;

        for (i=1; i<n; i++) {
            cc_mul_add_add_wide_limb_limb(&wide, mx, p[i], hi, a[i]);
            tmp[i - 1] = wide.lo;
            hi = wide.hi;
        }
        tmp[i-1] = hi;
        a = tmp;
    }

    /* this is needed only if input can be non-fully-reduced */
    for (borrow=0, i=0; i<n; i++) {
        cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
        ret[i] = wide.lo;
        borrow = wide.hi & 1;
    }

    mask = 0 - borrow;
    BLST_launder(mask);

    for(i=0; i<n; i++)
        ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);

    BLST_VLA_FREE(tmp);
}

#define FROM_MONT_IMPL(bits) \
HOST_DEVICE_FUNC inline void cc_from_mont_##bits(vec##bits ret, const vec##bits a, \
                             const vec##bits p, limb_t n0) \
{   cc_from_mont_n(ret, a, p, n0, NLIMBS(bits));   }

FROM_MONT_IMPL(256)
FROM_MONT_IMPL(384)

inline void cc_redc_mont_n(limb_t ret[], const limb_t a[],
                        const limb_t p[], limb_t n0, size_t n)
{
    __builtin_assume(n != 0 && n%2 == 0);
    wide_limb_t wide;
    limb_t mask, carry, borrow, mx, hi;
    BLST_VLA_NEW(limb_t, tmp, n);

    const limb_t *b = a;
    size_t i, j;

    for (j=0; j<n; j++) {
        mx = n0*b[0];
        cc_mul_add_wide_limb_limb(&wide, mx, p[0], b[0]);
        hi = wide.hi;
        for (i=1; i<n; i++) {
			cc_mul_add_add_wide_limb_limb(&wide, mx, p[i], hi, b[i]);
			tmp[i - 1] = wide.lo;
			hi = wide.hi;
        }
        tmp[i-1] = hi;
        b = tmp;
    }

    for (carry=0, i=0; i<n; i++) {
		cc_add_carry_wide_limb_limb(&wide, a[n + i], tmp[i], carry);
		tmp[i] = wide.lo;
		carry = wide.hi;
    }

    for (borrow=0, i=0; i<n; i++) {
		cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
		ret[i] = wide.lo;
		borrow = wide.hi & 1;
    }

    mask = carry - borrow;
    BLST_launder(mask);

    for(i=0; i<n; i++)
        ret[i] = (ret[i] & ~mask) | (tmp[i] & mask);

    BLST_VLA_FREE(tmp);
}

#define REDC_MONT_IMPL(bits, bits2) \
HOST_DEVICE_FUNC inline void cc_redc_mont_##bits(vec##bits ret, const vec##bits2 a, \
                             const vec##bits p, limb_t n0) \
{   cc_redc_mont_n(ret, a, p, n0, NLIMBS(bits));   }

REDC_MONT_IMPL(256, 512)
REDC_MONT_IMPL(384, 768)

inline void cc_rshift_mod_n(limb_t ret[], const limb_t a[], size_t count,
                         const limb_t p[], size_t n)
{
    __builtin_assume(count != 0);
    __builtin_assume(n != 0 && n%2 == 0);
    limb_t mask, carry, limb, next;
    size_t i;

    while (count--) {
        mask = 0 - (a[0] & 1);
        BLST_launder(mask);
        for (carry=0, i=0; i<n; i++) {
			wide_limb_t wide;
			cc_add_carry_wide_limb_limb(&wide, a[i], p[i] & mask, carry);
			ret[i] = wide.lo;
			carry = wide.hi;
        }

        for (next=ret[0], i=0; i<n-1; i++) {
            limb = next >> 1;
            next = ret[i+1];
            ret[i] = limb | next << (BLST_LIMB_SIZE-1);
        }
        ret[i] = next >> 1 | carry << (BLST_LIMB_SIZE-1);

        a = ret;
    }
}

#define RSHIFT_MOD_IMPL(bits) \
HOST_DEVICE_FUNC inline void cc_rshift_mod_##bits(vec##bits ret, const vec##bits a, size_t count, \
                              const vec##bits p) \
{   cc_rshift_mod_n(ret, a, count, p, NLIMBS(bits));   }

RSHIFT_MOD_IMPL(256)
RSHIFT_MOD_IMPL(384)

#define DIV_BY_2_MOD_IMPL(bits) \
HOST_DEVICE_FUNC inline void cc_div_by_2_mod_##bits(vec##bits ret, const vec##bits a, \
                                const vec##bits p) \
{   cc_rshift_mod_n(ret, a, 1, p, NLIMBS(bits));   }

DIV_BY_2_MOD_IMPL(384)

inline limb_t cc_sgn0_pty_mod_n(const limb_t a[], const limb_t p[], size_t n)
{
    __builtin_assume(n != 0);
    limb_t carry, borrow, ret;
    BLST_VLA_NEW(limb_t, tmp, n);

    size_t i;

    ret = a[0] & 1; /* parity */

    for (carry=0, i=0; i<n; i++) {
        limb_t a_i = a[i];
        tmp[i] = a_i<<1 | carry;
        carry = a_i>>(BLST_LIMB_SIZE-1);
    }

    for (borrow=0, i=0; i<n; i++) {
		wide_limb_t wide;
		cc_sub_borrow_wide_limb_limb(&wide, tmp[i], p[i], borrow);
		borrow = wide.hi & 1;
    }

    ret |= ((carry - borrow) & 2) ^ 2;

    BLST_VLA_FREE(tmp);
    return ret;
}

HOST_DEVICE_FUNC inline limb_t cc_sgn0_pty_mod_384(const vec384 a, const vec384 p)
{   return cc_sgn0_pty_mod_n(a, p, NLIMBS(384));   }

HOST_DEVICE_FUNC inline limb_t cc_sgn0_pty_mont_384(const vec384 a, const vec384 p, limb_t n0)
{
    vec384 tmp;

    cc_from_mont_n(tmp, a, p, n0, NLIMBS(384));

    return cc_sgn0_pty_mod_n(tmp, p, NLIMBS(384));
}

HOST_DEVICE_FUNC inline limb_t cc_sgn0_pty_mod_384x(const vec384x a, const vec384 p)
{
    limb_t re, im, sign, prty;

    re = cc_sgn0_pty_mod_n(a[0], p, NLIMBS(384));
    im = cc_sgn0_pty_mod_n(a[1], p, NLIMBS(384));

    /* a->im!=0 ? sgn0(a->im) : sgn0(a->re) */
    sign = (limb_t)0 - vec_is_zero(a[1], sizeof(vec384));
    sign = (re & sign) | (im & ~sign);

    /* a->re==0 ? prty(a->im) : prty(a->re) */
    prty = (limb_t)0 - vec_is_zero(a[0], sizeof(vec384));
    prty = (im & prty) | (re & ~prty);

    return (sign & 2) | (prty & 1);
}

HOST_DEVICE_FUNC inline limb_t cc_sgn0_pty_mont_384x(const vec384x a, const vec384 p, limb_t n0)
{
    vec384x tmp;

    cc_from_mont_n(tmp[0], a[0], p, n0, NLIMBS(384));
    cc_from_mont_n(tmp[1], a[1], p, n0, NLIMBS(384));

    return cc_sgn0_pty_mod_384x(tmp, p);
}

HOST_DEVICE_FUNC inline void cc_mul_mont_384x(vec384x ret, const vec384x a, const vec384x b,
                          const vec384 p, limb_t n0)
{
    vec384 aa, bb, cc;

    cc_add_mod_n(aa, a[0], a[1], p, NLIMBS(384));
    cc_add_mod_n(bb, b[0], b[1], p, NLIMBS(384));
    cc_mul_mont_n(bb, bb, aa, p, n0, NLIMBS(384));
    cc_mul_mont_n(aa, a[0], b[0], p, n0, NLIMBS(384));
    cc_mul_mont_n(cc, a[1], b[1], p, n0, NLIMBS(384));
    cc_sub_mod_n(ret[0], aa, cc, p, NLIMBS(384));
    cc_sub_mod_n(ret[1], bb, aa, p, NLIMBS(384));
    cc_sub_mod_n(ret[1], ret[1], cc, p, NLIMBS(384));
}


/*
 * mul_mont_n without final conditional subtraction, which implies
 * that modulus is one bit short, which in turn means that there are
 * no carries to handle between iterations...
 */
HOST_DEVICE_FUNC inline void cc_mul_mont_nonred_n(limb_t ret[], const limb_t a[], const limb_t b[],
                              const limb_t p[], limb_t n0, size_t n)
{
    __builtin_assume(n != 0 && n%2 == 0);
	wide_limb_t wide;
    limb_t mx, hi;
    BLST_VLA_NEW(limb_t, tmp, n+1);

    size_t i, j;

    for (mx=b[0], hi=0, i=0; i<n; i++) {
		cc_mul_add_wide_limb_limb(&wide, mx, a[i], hi);
		tmp[i] = wide.lo;
		hi = wide.hi;
    }
    mx = n0*tmp[0];
    tmp[i] = hi;

    for (j=0; ; ) {
		cc_mul_add_wide_limb_limb(&wide, mx, p[0], tmp[0]);
		hi = wide.hi;

        for (i=1; i<n; i++) {
			cc_mul_add_add_wide_limb_limb(&wide, mx, p[i], hi, tmp[i]);
			tmp[i - 1] = wide.lo;
			hi = wide.hi;
        }
        tmp[i-1] = tmp[i] + hi;

        if (++j==n)
            break;

        for (mx=b[j], hi=0, i=0; i<n; i++) {
			cc_mul_add_add_wide_limb_limb(&wide, mx, a[i], hi, tmp[i]);
			tmp[i] = wide.lo;
			hi = wide.hi;
        }
        mx = n0*tmp[0];
        tmp[i] = hi;

    }

    vec_copy(ret, tmp,  n*sizeof(limb_t));

    BLST_VLA_FREE(tmp);

}


HOST_DEVICE_FUNC inline void cc_sqr_n_mul_mont_384(vec384 ret, const vec384 a, size_t count,
                        const vec384 p, limb_t n0, const vec384 b)
{
    __builtin_assume(count != 0);
    while(count--) {
        cc_mul_mont_nonred_n(ret, a, a, p, n0, NLIMBS(384));
        a = ret;
    }
    cc_mul_mont_n(ret, ret, b, p, n0, NLIMBS(384));
}

HOST_DEVICE_FUNC inline void cc_sqr_mont_382x(vec384x ret, const vec384x a,
                          const vec384 p, limb_t n0)
{
	wide_limb_t wide;
    limb_t mask, carry, borrow;
    size_t i;
    vec384 t0, t1;

    /* "add_mod_n(t0, a[0], a[1], p, NLIMBS(384));" */
    for (carry=0, i=0; i<NLIMBS(384); i++) {
		cc_add_carry_wide_limb_limb(&wide, a[0][i], a[1][i], carry);
		t0[i] = wide.lo;
		carry = wide.hi;
    }

    /* "sub_mod_n(t1, a[0], a[1], p, NLIMBS(384));" */
    for (borrow=0, i=0; i<NLIMBS(384); i++) {
		cc_sub_borrow_wide_limb_limb(&wide, a[0][i], a[1][i], borrow);
		t1[i] = wide.lo;
		borrow = wide.hi & 1;
    }
    mask = 0 - borrow;
    BLST_launder(mask);

    /* "mul_mont_n(ret[1], a[0], a[1], p, n0, NLIMBS(384));" */
    cc_mul_mont_nonred_n(ret[1], a[0], a[1], p, n0, NLIMBS(384));

    /* "add_mod_n(ret[1], ret[1], ret[1], p, NLIMBS(384));" */
    for (carry=0, i=0; i<NLIMBS(384); i++) {
        limb_t a_i = ret[1][i];
        ret[1][i] = a_i<<1 | carry;
        carry = a_i>>(BLST_LIMB_SIZE-1);
    }

    /* "mul_mont_n(ret[0], t0, t1, p, n0, NLIMBS(384));" */
    cc_mul_mont_nonred_n(ret[0], t0, t1, p, n0, NLIMBS(384));

    /* account for t1's sign... */
    for (borrow=0, i=0; i<NLIMBS(384); i++) {
		cc_sub_borrow_wide_limb_limb(&wide, ret[0][i], t0[i] & mask, borrow);
		ret[0][i] = wide.lo;
		borrow = wide.hi & 1;
    }
    mask = 0 - borrow;
    BLST_launder(mask);
    for (carry=0, i=0; i<NLIMBS(384); i++) {
		cc_add_carry_wide_limb_limb(&wide, ret[0][i], p[i] & mask, carry);
		ret[0][i] = wide.lo;
		carry = wide.hi;
    }
}

HOST_DEVICE_FUNC inline void lshift_mod_v384(vec384 ret, const vec384 a, size_t n,
    const vec384 mod)
{
    while (n--)
        cc_add_mod_384(ret, a, a, mod), a = ret;
}

HOST_DEVICE_FUNC inline void cc_mul_by_8_mod_384(vec384 ret, const vec384 a, const vec384 mod)
{
    lshift_mod_v384(ret, a, 3, mod);
}

HOST_DEVICE_FUNC inline void cc_mul_by_3_mod_384x(vec384x ret, const vec384x a, const vec384 mod)
{
    cc_mul_by_3_mod_384(ret[0], a[0], mod);
    cc_mul_by_3_mod_384(ret[1], a[1], mod);
}
HOST_DEVICE_FUNC inline  void cc_mul_by_8_mod_384x(vec384x ret, const vec384x a, const vec384 mod)
{
    cc_mul_by_8_mod_384(ret[0], a[0], mod);
    cc_mul_by_8_mod_384(ret[1], a[1], mod);
}
HOST_DEVICE_FUNC inline  void cc_mul_by_1_plus_i_mod_384x(vec384x ret, const vec384x a,
    const vec384 mod)
{
    vec384 t;

    cc_add_mod_384(t, a[0], a[1], mod);
    cc_sub_mod_384(ret[0], a[0], a[1], mod);
    vec_copy(ret[1], t, sizeof(t));
}
HOST_DEVICE_FUNC inline void cc_add_mod_384x(vec384x ret, const vec384x a, const vec384x b,
    const vec384 mod)
{
    cc_add_mod_384(ret[0], a[0], b[0], mod);
    cc_add_mod_384(ret[1], a[1], b[1], mod);
}
HOST_DEVICE_FUNC inline void cc_sub_mod_384x(vec384x ret, const vec384x a, const vec384x b,
    const vec384 mod)
{
    cc_sub_mod_384(ret[0], a[0], b[0], mod);
    cc_sub_mod_384(ret[1], a[1], b[1], mod);
}
HOST_DEVICE_FUNC inline void cc_lshift_mod_384x(vec384x ret, const vec384x a, size_t n,
    const vec384 mod)
{
    cc_lshift_mod_384(ret[0], a[0], n, mod);
    cc_lshift_mod_384(ret[1], a[1], n, mod);
}

//#if defined(__BLST_NO_ASM__) && !(defined(__ADX__) && !defined(__BLST_PORTABLE__))
HOST_DEVICE_FUNC inline void cc_sqr_mont_384x(vec384x ret, const vec384x a, const vec384 mod, limb_t n0)
{
    vec384 t0, t1;

    cc_add_mod_384(t0, a[0], a[1], mod);
    cc_sub_mod_384(t1, a[0], a[1], mod);

    cc_mul_mont_384(ret[1], a[0], a[1], mod, n0);
    cc_add_mod_384(ret[1], ret[1], ret[1], mod);

    cc_mul_mont_384(ret[0], t0, t1, mod, n0);
}
//#endif


#if defined(__GNUC__) || defined(__clang__)
# define MSB(x) ({ limb_t ret = (x) >> (BLST_LIMB_SIZE-1); BLST_launder(ret); ret; })
#else
# define MSB(x) ((x) >> (BLST_LIMB_SIZE-1))
#endif

HOST_DEVICE_FUNC inline size_t cc_num_bits(limb_t l)
{
    limb_t x, mask;
    size_t bits = is_zero(l) ^ 1;

    if (sizeof(limb_t) == 8) {
        x = l >> (32 & (8*sizeof(limb_t)-1));
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
}

#if defined(__clang_major__) && __clang_major__>7
__attribute__((optnone))
#endif
HOST_DEVICE_FUNC inline  limb_t cc_lshift_2(limb_t hi, limb_t lo, size_t l)
{
    size_t r = BLST_LIMB_SIZE - l;
    limb_t mask = 0 - (is_zero(l)^1);
    return (hi << (l&(BLST_LIMB_SIZE-1))) | ((lo & mask) >> (r&(BLST_LIMB_SIZE-1)));
}

/*
 * https://eprint.iacr.org/2020/972 with 'k' being BLST_LIMB_SIZE-1.
 */
HOST_DEVICE_FUNC inline void cc_ab_approximation_n(limb_t a_[2], const limb_t a[],
                               limb_t b_[2], const limb_t b[], size_t n)
{
    __builtin_assume(n != 0 && n%2 == 0);
    limb_t a_hi, a_lo, b_hi, b_lo, mask;
    size_t i;

    i = n-1;
    a_hi = a[i],    a_lo = a[i-1];
    b_hi = b[i],    b_lo = b[i-1];
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

HOST_DEVICE_FUNC inline void cc_inner_loop_n(factors *fg, const limb_t a_[2], const limb_t b_[2],
                         size_t n)
{
    __builtin_assume(n != 0);
	wide_limb_t wide;
    limb_t f0 = 1, g0 = 0, f1 = 0, g1 = 1;
    limb_t a_lo, a_hi, b_lo, b_hi, t_lo, t_hi, odd, borrow, xorm;

    a_lo = a_[0], a_hi = a_[1];
    b_lo = b_[0], b_hi = b_[1];

    while(n--) {
        odd = 0 - (a_lo&1);

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
        a_lo >>= 1; a_lo |= a_hi << (BLST_LIMB_SIZE-1);
        a_hi >>= 1;
    }

    fg->f0 = f0, fg->g0 = g0, fg->f1 = f1, fg->g1= g1;
}

HOST_DEVICE_FUNC inline limb_t cc_cneg_n(limb_t ret[], const limb_t a[], limb_t neg, size_t n)
{
    __builtin_assume(n != 0);
	wide_limb_t wide;
    limb_t carry;
    size_t i;

    for (carry=neg&1, i=0; i<n; i++) {
		cc_add_wide_limb_limb(&wide, a[i] ^ neg, carry);
		ret[i] = wide.lo;
		carry = wide.hi;
    }

    return 0 - MSB(wide.lo);
}

HOST_DEVICE_FUNC inline limb_t cc_add_n(limb_t ret[], const limb_t a[], limb_t b[], size_t n)
{
    __builtin_assume(n != 0);
	wide_limb_t wide;
    limb_t carry;
    size_t i;

    for (carry=0, i=0; i<n; i++) {
		cc_add_carry_wide_limb_limb(&wide, a[i], b[i], carry);
		ret[i] = wide.lo;
		carry = wide.hi;
    }

    return carry;
}

HOST_DEVICE_FUNC inline limb_t cc_umul_n(limb_t ret[], const limb_t a[], limb_t b, size_t n)
{
    __builtin_assume(n != 0);
	wide_limb_t wide;
    limb_t hi;
    size_t i;

    for (hi=0, i=0; i<n; i++) {
		cc_mul_add_wide_limb_limb(&wide, b, a[i], hi);
		ret[i] = wide.lo;
		hi = wide.hi;
    }

    return hi;
}

HOST_DEVICE_FUNC inline limb_t cc_smul_n_shift_n(limb_t ret[], const limb_t a[], limb_t *f_,
                                           const limb_t b[], limb_t *g_,
                                           size_t n)
{
    __builtin_assume(n != 0);
    limb_t f, g, neg, carry, hi;
    BLST_VLA_NEW(limb_t, a_, n + 1);
    BLST_VLA_NEW(limb_t, b_, n + 1);
    size_t i;

    /* |a|*|f_| */
    f = *f_;
    neg = 0 - MSB(f);
    f = (f ^ neg) - neg;            /* ensure |f| is positive */
    (void)cc_cneg_n(a_, a, neg, n);
    hi = cc_umul_n(a_, a_, f, n);
    a_[n] = hi - (f & neg);

    /* |b|*|g_| */
    g = *g_;
    neg = 0 - MSB(g);
    g = (g ^ neg) - neg;            /* ensure |g| is positive */
    (void)cc_cneg_n(b_, b, neg, n);
    hi = cc_umul_n(b_, b_, g, n);
    b_[n] = hi - (g & neg);

    /* |a|*|f_| + |b|*|g_| */
    (void)cc_add_n(a_, a_, b_, n+1);

    /* (|a|*|f_| + |b|*|g_|) >> k */
    for (carry=a_[0], i=0; i<n; i++) {
        hi = carry >> (BLST_LIMB_SIZE-2);
        carry = a_[i+1];
        ret[i] = hi | (carry << 2);
    }

    /* ensure result is non-negative, fix up |f_| and |g_| accordingly */
    neg = 0 - MSB(carry);
    *f_ = (*f_ ^ neg) - neg;
    *g_ = (*g_ ^ neg) - neg;
    (void)cc_cneg_n(ret, ret, neg, n);

	BLST_VLA_FREE(a_);
	BLST_VLA_FREE(b_);

    return neg;
}

HOST_DEVICE_FUNC inline limb_t cc_smul_2n(limb_t ret[], const limb_t u[], limb_t f,
                                    const limb_t v[], limb_t g, size_t n)
{
    __builtin_assume(n != 0);
    limb_t neg, hi;

    BLST_VLA_NEW(limb_t, u_, n);
    BLST_VLA_NEW(limb_t, v_, n);

    /* |u|*|f_| */
    neg = 0 - MSB(f);
    f = (f ^ neg) - neg;            /* ensure |f| is positive */
    neg = cc_cneg_n(u_, u, neg, n);
    hi = cc_umul_n(u_, u_, f, n) - (f&neg);

    /* |v|*|g_| */
    neg = 0 - MSB(g);
    g = (g ^ neg) - neg;            /* ensure |g| is positive */
    neg = cc_cneg_n(v_, v, neg, n);
    hi += cc_umul_n(v_, v_, g, n) - (g&neg);

    /* |u|*|f_| + |v|*|g_| */
    hi += cc_add_n(ret, u_, v_, n);

    BLST_VLA_FREE(u_);
    BLST_VLA_FREE(v_);

    return hi;
}

HOST_DEVICE_FUNC inline void cc_ct_inverse_mod_n(limb_t ret[], const limb_t inp[],
                             const limb_t mod[], const limb_t modx[], size_t n)
{
    __builtin_assume(n != 0 && n%2 == 0);
	wide_limb_t wide;
	BLST_VLA_NEW(limb_t, a, n);
	BLST_VLA_NEW(limb_t, b, n);
	BLST_VLA_NEW(limb_t, u, 2 * n);
	BLST_VLA_NEW(limb_t, v, 2 * n);
	BLST_VLA_NEW(limb_t, t, 2 * n);

    limb_t a_[2], b_[2], sign, carry, top;
    factors fg;
    size_t i;

    vec_copy(a, inp, n * sizeof(limb_t));
    vec_copy(b, mod, n * sizeof(limb_t));
    vec_zero(u, n*2); u[0] = 1;
    vec_zero(v, n*2 * sizeof(limb_t));

    for (i=0; i<(2*n*BLST_LIMB_SIZE)/(BLST_LIMB_SIZE-2); i++) {
        cc_ab_approximation_n(a_, a, b_, b, n);
        cc_inner_loop_n(&fg, a_, b_, BLST_LIMB_SIZE-2);
        (void)cc_smul_n_shift_n(t, a, &fg.f0, b, &fg.g0, n);
        (void)cc_smul_n_shift_n(b, a, &fg.f1, b, &fg.g1, n);
        vec_copy(a, t, n * sizeof(limb_t));
        cc_smul_2n(t, u, fg.f0, v, fg.g0, 2*n);
        cc_smul_2n(v, u, fg.f1, v, fg.g1, 2*n);
        vec_copy(u, t, 2*n * sizeof(limb_t));
    }

    cc_inner_loop_n(&fg, a, b, (2*n*BLST_LIMB_SIZE)%(BLST_LIMB_SIZE-2));
    top = cc_smul_2n(ret, u, fg.f1, v, fg.g1, 2*n);

    sign = 0 - MSB(top);    /* top is 1, 0 or -1 */
    for (carry=0, i=0; i<n; i++) {
		cc_add_carry_wide_limb_limb(&wide, ret[n + i], modx[i] & sign, carry);
		ret[n + i] = wide.lo;
		carry = wide.hi;
    }
    top += carry;
    sign = 0 - top;         /* top is 1, 0 or -1 */
    top |= sign;
    for (i=0; i<n; i++)
        a[i] = modx[i] & top;
    (void)cc_cneg_n(a, a, 0 - MSB(sign), n);
    cc_add_n(ret+n, ret+n, a, n);

    BLST_VLA_FREE(a);
	BLST_VLA_FREE(b);
	BLST_VLA_FREE(u);
	BLST_VLA_FREE(v);
	BLST_VLA_FREE(t);
}

#define CT_INVERSE_MOD_IMPL(bits, bits2) \
HOST_DEVICE_FUNC inline void cc_ct_inverse_mod_##bits(vec##bits2 ret, const vec##bits inp, \
                                  const vec##bits mod, const vec##bits modx) \
{   cc_ct_inverse_mod_n(ret, inp, mod, modx, NLIMBS(bits));   }

CT_INVERSE_MOD_IMPL(256, 512)
CT_INVERSE_MOD_IMPL(384, 768)

/*
 * Copy of inner_loop_n above, but with |L| updates.
 */
HOST_DEVICE_FUNC inline limb_t cc_legendre_loop_n(limb_t L, factors *fg, const limb_t a_[2],
                              const limb_t b_[2], size_t n)
{
    __builtin_assume(n != 0);
	wide_limb_t wide;
    limb_t f0 = 1, g0 = 0, f1 = 0, g1 = 1;
    limb_t a_lo, a_hi, b_lo, b_hi, t_lo, t_hi, odd, borrow, xorm;

    a_lo = a_[0], a_hi = a_[1];
    b_lo = b_[0], b_hi = b_[1];

    while(n--) {
        odd = 0 - (a_lo&1);

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
        a_lo >>= 1; a_lo |= a_hi << (BLST_LIMB_SIZE-1);
        a_hi >>= 1;

        L += (b_lo + 2) >> 2;
    }

    fg->f0 = f0, fg->g0 = g0, fg->f1 = f1, fg->g1 = g1;

    return L;
}

HOST_DEVICE_FUNC inline bool_t cc_ct_is_sqr_mod_n(const limb_t inp[], const limb_t mod[], size_t n)
{
    __builtin_assume(n != 0 && n%2 == 0);
	BLST_VLA_NEW(limb_t, a, n);
	BLST_VLA_NEW(limb_t, b, n);
	BLST_VLA_NEW(limb_t, t, n);
    limb_t a_[2], b_[2], neg, L = 0;
    factors fg;
    size_t i;

    vec_copy(a, inp, n * sizeof(limb_t));
    vec_copy(b, mod, n * sizeof(limb_t));

    for (i=0; i<(2*n*BLST_LIMB_SIZE)/(BLST_LIMB_SIZE-2); i++) {
        cc_ab_approximation_n(a_, a, b_, b, n);
        L = cc_legendre_loop_n(L, &fg, a_, b_, BLST_LIMB_SIZE-2);
        neg = cc_smul_n_shift_n(t, a, &fg.f0, b, &fg.g0, n);
        (void)cc_smul_n_shift_n(b, a, &fg.f1, b, &fg.g1, n);
        vec_copy(a, t, n * sizeof(limb_t));
        L += (b[0] >> 1) & neg;
    }

    L = cc_legendre_loop_n(L, &fg, a, b, (2*n*BLST_LIMB_SIZE)%(BLST_LIMB_SIZE-2));

    BLST_VLA_FREE(a);
    BLST_VLA_FREE(b);
    BLST_VLA_FREE(t);

    return (L & 1) ^ 1;
}

#define CT_IS_SQR_MOD_IMPL(bits) \
HOST_DEVICE_FUNC inline bool_t cc_ct_is_square_mod_##bits(const vec##bits inp, \
                                      const vec##bits mod) \
{   return cc_ct_is_sqr_mod_n(inp, mod, NLIMBS(bits));   }

CT_IS_SQR_MOD_IMPL(384)

/*
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
HOST_DEVICE_FUNC inline limb_t cc_div_3_limbs(const limb_t div_top[2], limb_t d_lo, limb_t d_hi)
{
	wide_limb_t wide;
    limb_t r_lo = div_top[0], r_hi = div_top[1];
    limb_t Q = 0, mask, borrow, rx;
    size_t i;

    for (i = 0; i < BLST_LIMB_SIZE; i++) {
        /* "borrow, Rx = R - D" */
		cc_sub_borrow_wide_limb_limb(&wide, r_lo, d_lo, 0);
		rx = wide.lo;
		borrow = wide.hi & 1;
		cc_sub_borrow_wide_limb_limb(&wide, r_hi, d_hi, borrow);
		r_hi = wide.lo;
		borrow = wide.hi & 1;

        /* "if (R >= D) R -= D" */
        r_lo = ((r_lo ^ rx) & borrow) ^ rx;
        r_hi = ((r_hi ^ rx) & borrow) ^ rx;

        Q <<= 1;
        Q |= ~borrow & 1;

        /* "D >>= 1" */
        d_lo >>= 1; d_lo |= d_hi << (BLST_LIMB_SIZE - 1);
        d_hi >>= 1;
    }

    mask = 0 - MSB(Q);  /* does it overflow? */

    /* "borrow, Rx = R - D" */
	cc_sub_borrow_wide_limb_limb(&wide, r_lo, d_lo, 0);
	rx = wide.lo;
	borrow = wide.hi & 1;
	cc_sub_borrow_wide_limb_limb(&wide, r_hi, d_hi, borrow);
	rx = wide.lo;
	borrow = wide.hi & 1;

    Q <<= 1;
    Q |= borrow ^ 1;

    return (Q | mask);
}

HOST_DEVICE_FUNC inline limb_t cc_quot_rem_n(limb_t *div_rem, const limb_t *divisor,
                                          limb_t quotient, size_t n)
{
    __builtin_assume(n != 0 && n%2 == 0);
	wide_limb_t wide;
    limb_t carry, mask, borrow;
	BLST_VLA_NEW(limb_t, tmp, n + 1);
    size_t i;

    /* divisor*quotient */
    for (carry=0, i=0; i<n; i++) {
		cc_mul_add_wide_limb_limb(&wide, quotient, divisor[i], carry);
		tmp[i] = wide.lo;
		carry = wide.hi;
    }
    tmp[i] = carry;

    /* remainder = dividend - divisor*quotient */
    for (borrow=0, i=0; i<=n; i++) {
		cc_sub_borrow_wide_limb_limb(&wide, div_rem[i], tmp[i] + borrow, 0);
		tmp[i] = wide.lo;
		borrow = wide.hi & 1;
    }

    mask = 0 - borrow;
    BLST_launder(mask);

    /* if quotient was off by one, add divisor to the remainder */
    for (carry=0, i=0; i<n; i++) {
		cc_add_carry_wide_limb_limb(&wide, tmp[i], divisor[i] & mask, carry);
		div_rem[i] = wide.lo;
		carry = wide.hi & 1;
    }

	BLST_VLA_FREE(tmp);

    return (div_rem[i] = quotient + mask);
}

HOST_DEVICE_FUNC inline limb_t cc_quot_rem_128(limb_t *div_rem, const limb_t *divisor,
                                            limb_t quotient)
{   return cc_quot_rem_n(div_rem, divisor, quotient, NLIMBS(128));   }

HOST_DEVICE_FUNC inline limb_t cc_quot_rem_64(limb_t *div_rem, const limb_t *divisor,
                                           limb_t quotient)
{   return cc_quot_rem_n(div_rem, divisor, quotient, NLIMBS(64));   }

/*
 * Unlock reference implementations in vect.c
 */
//#define mul_by_8_mod_384 mul_by_8_mod_384
//#define mul_by_8_mod_384x mul_by_8_mod_384x
//#define mul_by_3_mod_384x mul_by_3_mod_384x
//#define mul_by_1_plus_i_mod_384x mul_by_1_plus_i_mod_384x
//#define add_mod_384x add_mod_384x
//#define sub_mod_384x sub_mod_384x
//#define lshift_mod_384x lshift_mod_384x
//#define sqr_mont_384x sqr_mont_384x

HOST_DEVICE_FUNC inline void cc_vec_prefetch(const void *ptr, size_t len)
{  
	// prefetch if we have the intrinsic
#if defined(__GNUC__) || defined(__clang__)
	__builtin_prefetch(ptr, 0, 3);
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
	_mm_prefetch((const char*)ptr, _MM_HINT_T0);
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

HOST_DEVICE_FUNC inline void cc_blst_sha256_block_data_order(unsigned int *v, const void *inp,
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
    const unsigned char *data = (const unsigned char*)inp;
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
            l  = (unsigned int)data[0] << 24;
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

HOST_DEVICE_FUNC inline void cc_blst_sha256_hcopy(unsigned int dst[8], const unsigned int src[8])
{
    size_t i;

    for (i=0; i<8; i++)
        dst[i] = src[i];
}

HOST_DEVICE_FUNC inline void cc_blst_sha256_emit(unsigned char md[32], const unsigned int h[8])
{
    size_t i;

    for (i=0; i<8; i++, md+=4) {
        unsigned int h_i = h[i];
        md[0] = (unsigned char)(h_i >> 24);
        md[1] = (unsigned char)(h_i >> 16);
        md[2] = (unsigned char)(h_i >> 8);
        md[3] = (unsigned char)h_i;
    }
}

HOST_DEVICE_FUNC inline void cc_blst_sha256_bcopy(void *dst_, const void *src_, size_t len)
{
    unsigned char *dst = (unsigned char*)dst_;
    const unsigned char *src = (const unsigned char*)src_;
    size_t i;

    for (i=0; i<len; i++)
        dst[i] = src[i];
}


#undef BLST_VLA_NEW
#undef BLST_VLA_FREE

#endif