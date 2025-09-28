/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
/*
 * Why this file? Overall goal is to ensure that all internal calls
 * remain internal after linking application. This is to both
 *
 * a) minimize possibility of external name conflicts (since all
 *    non-blst-prefixed and [assembly subroutines] remain static);
 * b) preclude possibility of unintentional internal reference
 *    overload in shared library context (one can achieve same
 *    effect with -Bsymbolic, but we don't want to rely on end-user
 *    to remember to use it);
 */

#include "blst/fields/Fr.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/common/bytes.h"
BLST_EXTERN_C_BEGIN_

using namespace blst;

/*
 * BLS12-381-specific Fr shortcuts to assembly.
 */
void blst_fr_add(Fr& ret, const Fr& a, const Fr& b) { fr_add(ret, a, b); }
void blst_fr_sub(Fr& ret, const Fr& a, const Fr& b) { fr_sub(ret, a, b); }
void blst_fr_mul_by_3(Fr& ret, const Fr& a) { fr_mul_by_3(ret, a); }
void blst_fr_lshift(Fr& ret, const Fr& a, size_t count) { fr_lshift(ret, a, count); }
void blst_fr_rshift(Fr& ret, const Fr& a, size_t count) { fr_rshift(ret, a, count); }
void blst_fr_mul(Fr& ret, const Fr& a, const Fr& b) { fr_mul(ret, a, b); }
void blst_fr_ct_bfly(Fr& x0, Fr& x1, const Fr& twiddle) { fr_ct_bfly(x0, x1, twiddle); }
void blst_fr_gs_bfly(Fr& x0, Fr& x1, const Fr& twiddle) { fr_gs_bfly(x0, x1, twiddle); }
void blst_fr_sqr(Fr& ret, const Fr& a) { fr_sqr(ret, a); }
void blst_fr_cneg(Fr& ret, const Fr& a, int flag) { fr_cneg(ret, a, flag); }
void blst_fr_to(Fr& ret, const Fr& a) { fr_to(ret, a); }
void blst_fr_from(Fr& ret, const Fr& a) { fr_from(ret, a); }
void blst_fr_from_scalar(Fr& ret, const pow256 a) { fr_from_scalar(ret, a); }
void blst_scalar_from_fr(pow256 ret, const Fr& a) { scalar_from_fr(ret, a); }
int blst_scalar_fr_check(const pow256 a) { return scalar_fr_check(a); }
int blst_sk_check(const pow256 a) { return sk_check(a); }
int blst_sk_add_n_check(pow256 ret, const pow256 a, const pow256 b) { return sk_add_n_check(ret, a, b); }
int blst_sk_sub_n_check(pow256 ret, const pow256 a, const pow256 b) { return sk_sub_n_check(ret, a, b); }
int blst_sk_mul_n_check(pow256 ret, const pow256 a, const pow256 b) { return sk_mul_n_check(ret, a, b); }
void blst_sk_inverse(pow256 ret, const pow256 a) { sk_inverse(ret, a); }




/*
 * BLS12-381-specific Fp shortcuts to assembly.
 */
void blst_fp_add(Fp& ret, const Fp& a, const Fp& b)
{   Fp::add(ret, a, b);   }

void blst_fp_sub(Fp& ret, const Fp& a, const Fp& b)
{   Fp::sub(ret, a, b);   }

void blst_fp_mul_by_3(Fp& ret, const Fp& a)
{ Fp::mul_by_3(ret, a);   }

void blst_fp_mul_by_8(Fp& ret, const Fp& a)
{  Fp::mul_by_8(ret, a);   }

void blst_fp_lshift(Fp& ret, const Fp& a, size_t count)
{   Fp::lshift(ret, a, count);   }

void blst_fp_mul(Fp& ret, const Fp& a, const Fp& b)
{   Fp::mul(ret, a, b);   }

void blst_fp_sqr(Fp& ret, const Fp& a)
{   Fp::sqr(ret, a);   }

void blst_fp_cneg(Fp& ret, const Fp& a, int flag)
{   Fp::cneg(ret, a, is_zero(flag) ^ 1);   }

void blst_fp_to(Fp& ret, const Fp& a)
{   Fp::mul(ret, a, get_BLST_RR());   }

void blst_fp_from(Fp& ret, const Fp& a)
{   Fp::from_mont(ret, a);   }

/*
 * Fp serialization/deserialization.
 */
void blst_fp_from_uint32(Fp& ret, const unsigned int a[12])
{
    if (sizeof(limb_t) == 8) {
        int i;
        for (i = 0; i < 6; i++)
            ret[i] = a[2*i] | ((limb_t)a[2*i+1] << (32 & (8*sizeof(limb_t)-1)));
        a = (const unsigned int *)&ret;
    }
    Fp::mul(ret, *(const Fp *)a, get_BLST_RR());
}

void blst_uint32_from_fp(unsigned int ret[12], const Fp& a)
{
    if (sizeof(limb_t) == 4) {
        Fp::from_mont(*(Fp *)ret, a);
    } else {
        Fp out;
        int i;

        Fp::from_mont(out, a);
        for (i = 0; i < 6; i++) {
            limb_t limb = out[i];
            ret[2*i]   = (unsigned int)limb;
            ret[2*i+1] = (unsigned int)(limb >> (32 & (8*sizeof(limb)-1)));
        }
    }
}

void blst_fp_from_uint64(Fp& ret, const unsigned long long a[6])
{
    const union {
        long one;
        char little;
    } is_endian = { 1 };

    if (sizeof(limb_t) == 4 && !is_endian.little) {
        int i;
        for (i = 0; i < 6; i++) {
            unsigned long long ll= a[i];
            ret[2*i]   = (limb_t)ll;
            ret[2*i+1] = (limb_t)(ll >> 32);
        }
        a = (const unsigned long long *)ret;
    }
    Fp::mul(ret, *(const Fp*)a, get_BLST_RR());
}

void blst_uint64_from_fp(unsigned long long ret[6], const Fp& a)
{
    const union {
        long one;
        char little;
    } is_endian = { 1 };

    if (sizeof(limb_t) == 8 || is_endian.little) {
        Fp::from_mont(*(Fp *)ret, a);
    } else {
        Fp out;
        int i;

        Fp::from_mont(out, a);
        for (i = 0; i < 6; i++)
            ret[i] = out[2*i] | ((unsigned long long)out[2*i+1] << 32);
    }
}

void blst_fp_from_bendian(Fp& ret, const unsigned char a[48])
{
    Fp out;

    limbs_from_be_bytes(out, a, sizeof(Fp));
    Fp::mul(ret, out, get_BLST_RR());
}

void blst_bendian_from_fp(unsigned char ret[48], const Fp& a)
{
    Fp out;

    Fp::from_mont(out, a);
    be_bytes_from_limbs(ret, out, sizeof(Fp));
}

void blst_fp_from_lendian(Fp& ret, const unsigned char a[48])
{
    Fp out;

    limbs_from_le_bytes(out, a, sizeof(Fp));
    Fp::mul(ret, out, get_BLST_RR());
}

void blst_lendian_from_fp(unsigned char ret[48], const Fp& a)
{
    Fp out;

    Fp::from_mont(out, a);
    le_bytes_from_limbs(ret, out, sizeof(Fp));
}

/*
 * BLS12-381-specific Fp2 shortcuts to assembly.
 */
void blst_fp2_add(Fp2& ret, const Fp2& a, const Fp2& b)
{   Fp2::add(ret, a, b);   }

void blst_fp2_sub(Fp2& ret, const Fp2& a, const Fp2& b)
{   Fp2::sub(ret, a, b);   }

void blst_fp2_mul_by_3(Fp2& ret, const Fp2& a)
{   Fp2::mul_by_3(ret, a);   }

void blst_fp2_mul_by_8(Fp2& ret, const Fp2& a)
{   Fp2::mul_by_8(ret, a);   }

void blst_fp2_lshift(Fp2& ret, const Fp2& a, size_t count)
{   lshift_fp2(ret, a, count);    }

void blst_fp2_mul(Fp2& ret, const Fp2& a, const Fp2& b)
{   Fp2::mul(ret, a, b);   }

void blst_fp2_sqr(Fp2& ret, const Fp2& a)
{   Fp2::sqr(ret, a);   }

void blst_fp2_cneg(Fp2& ret, const Fp2& a, int flag)
{   Fp2::cneg(ret, a, is_zero(flag) ^ 1);   }

/*
 * Scalar serialization/deserialization.
 */
void blst_scalar_from_uint32(pow256 ret, const unsigned int a[8])
{
    const union {
        long one;
        char little;
    } is_endian = { 1 };
    size_t i;

    if ((uptr_t)ret==(uptr_t)a && is_endian.little)
        return;

    for(i = 0; i < 8; i++) {
        unsigned int w = a[i];
        *ret++ = (byte)w;
        *ret++ = (byte)(w >> 8);
        *ret++ = (byte)(w >> 16);
        *ret++ = (byte)(w >> 24);
    }
}

void blst_uint32_from_scalar(unsigned int ret[8], const pow256 a)
{
    const union {
        long one;
        char little;
    } is_endian = { 1 };
    size_t i;

    if ((uptr_t)ret==(uptr_t)a && is_endian.little)
        return;

    for(i = 0; i < 8; i++) {
        unsigned int w = (unsigned int)(*a++);
        w |= (unsigned int)(*a++) << 8;
        w |= (unsigned int)(*a++) << 16;
        w |= (unsigned int)(*a++) << 24;
        ret[i] = w;
    }
}

void blst_scalar_from_uint64(pow256 ret, const unsigned long long a[4])
{
    const union {
        long one;
        char little;
    } is_endian = { 1 };
    size_t i;

    if ((uptr_t)ret==(uptr_t)a && is_endian.little)
        return;

    for(i = 0; i < 4; i++) {
        unsigned long long w = a[i];
        *ret++ = (byte)w;
        *ret++ = (byte)(w >> 8);
        *ret++ = (byte)(w >> 16);
        *ret++ = (byte)(w >> 24);
        *ret++ = (byte)(w >> 32);
        *ret++ = (byte)(w >> 40);
        *ret++ = (byte)(w >> 48);
        *ret++ = (byte)(w >> 56);
    }
}

void blst_uint64_from_scalar(unsigned long long ret[4], const pow256 a)
{
    const union {
        long one;
        char little;
    } is_endian = { 1 };
    size_t i;

    if ((uptr_t)ret==(uptr_t)a && is_endian.little)
        return;

    for(i = 0; i < 4; i++) {
        unsigned long long w = (unsigned long long)(*a++);
        w |= (unsigned long long)(*a++) << 8;
        w |= (unsigned long long)(*a++) << 16;
        w |= (unsigned long long)(*a++) << 24;
        w |= (unsigned long long)(*a++) << 32;
        w |= (unsigned long long)(*a++) << 40;
        w |= (unsigned long long)(*a++) << 48;
        w |= (unsigned long long)(*a++) << 56;
        ret[i] = w;
    }
}

void blst_scalar_from_bendian(pow256 ret, const unsigned char a[32])
{
    Fr out;
    limbs_from_be_bytes(out, a, sizeof(out));
    le_bytes_from_limbs(ret, out, sizeof(out));
    vec_zero(out, sizeof(out));
}

void blst_bendian_from_scalar(unsigned char ret[32], const pow256 a)
{
    Fr out;
    limbs_from_le_bytes(out, a, sizeof(out));
    be_bytes_from_limbs(ret, out, sizeof(out));
    vec_zero(out, sizeof(out));
}

void blst_scalar_from_lendian(pow256 ret, const unsigned char a[32])
{
    size_t i;

    if ((uptr_t)ret==(uptr_t)a)
        return;

    for (i = 0; i < 32; i++)
        ret[i] = a[i];
}

void blst_lendian_from_scalar(unsigned char ret[32], const pow256 a)
{
    size_t i;

    if ((uptr_t)ret==(uptr_t)a)
        return;

    for (i = 0; i < 32; i++)
        ret[i] = a[i];
}


void blst_uint64_from_fr(unsigned long long ret[4], const Fr& a)
{
    const union {
        long one;
        char little;
    } is_endian = { 1 };

    if (sizeof(limb_t) == 8 || is_endian.little) {
        from_mont_v256(*(Fr *)ret, a, get_BLST_r(), BLST_r0);
    } else {
        Fr out;
        int i;

        from_mont_v256(out, a, get_BLST_r(), BLST_r0);
        for (i = 0; i < 4; i++)
            ret[i] = out[2*i] | ((unsigned long long)out[2*i+1] << 32);
        vec_zero(out, sizeof(out));
    }
}

int blst_scalar_from_le_bytes(pow256 out, const unsigned char *bytes, size_t n)
{
    size_t rem = (n - 1) % 32 + 1;
    struct { Fr out, digit; } t;
    limb_t ret;

    vec_zero(t.out, sizeof(t.out));

    n -= rem;
    limbs_from_le_bytes(t.out, bytes += n, rem);
    mul_mont_sparse_v256(t.out, get_BLST_rRR(), t.out, get_BLST_r(), BLST_r0);

    while (n) {
        limbs_from_le_bytes(t.digit, bytes -= 32, 32);
        add_mod_v256(t.out, t.out, t.digit, get_BLST_r());
        mul_mont_sparse_v256(t.out, get_BLST_rRR(), t.out, get_BLST_r(), BLST_r0);
        n -= 32;
    }

    from_mont_v256(t.out, t.out, get_BLST_r(), BLST_r0);

    ret = vec_is_zero(t.out, sizeof(t.out));
    le_bytes_from_limbs(out, t.out, 32);
    vec_zero(&t, sizeof(t));

    return (int)(ret^1);
}

int blst_scalar_from_be_bytes(pow256 out, const unsigned char *bytes, size_t n)
{
    size_t rem = (n - 1) % 32 + 1;
    struct { Fr out, digit; } t;
    limb_t ret;

    vec_zero(t.out, sizeof(t.out));

    limbs_from_be_bytes(t.out, bytes, rem);
    mul_mont_sparse_v256(t.out, get_BLST_rRR(), t.out, get_BLST_r(), BLST_r0);

    while (n -= rem) {
        limbs_from_be_bytes(t.digit, bytes += rem, 32);
        add_mod_v256(t.out, t.out, t.digit, get_BLST_r());
        mul_mont_sparse_v256(t.out, get_BLST_rRR(), t.out, get_BLST_r(), BLST_r0);
        rem = 32;
    }

    from_mont_v256(t.out, t.out, get_BLST_r(), BLST_r0);

    ret = vec_is_zero(t.out, sizeof(t.out));
    le_bytes_from_limbs(out, t.out, 32);
    vec_zero(&t, sizeof(t));

    return (int)(ret^1);
}

/*
 * Single-short SHA-256 hash function.
 */
#include "blst/sign/sha256.h"

void blst_sha256(unsigned char md[32], const void *msg, size_t len)
{
    SHA256_CTX ctx;

    sha256_init(&ctx);
    sha256_update(&ctx, msg, len);
    sha256_final(md, &ctx);
}

/*
 * Test facilitator.
 */
void blst_scalar_from_hexascii(pow256 ret, const char *hex)
{   bytes_from_hexascii(ret, sizeof(pow256), hex);   }

void blst_fp_from_hexascii(Fp& ret, const char *hex)
{
    limbs_from_hexascii(ret, sizeof(Fp), hex);
    Fp::mul(ret, ret, get_BLST_RR());
}
BLST_EXTERN_C_END_
