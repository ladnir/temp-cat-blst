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
#include "bindings/blst_scalar.h"
BLST_EXTERN_C_BEGIN_

using namespace blst;

/*
 * BLS12-381-specific Fr shortcuts to assembly.
 */
void blst_fr_add(blst_fr* ret, const blst_fr* a, const blst_fr* b) { Fr::add(*ret, *a, *b); }
void blst_fr_sub(blst_fr* ret, const blst_fr* a, const blst_fr* b) { Fr::sub(*ret, *a, *b); }
void blst_fr_mul_by_3(blst_fr* ret, const blst_fr* a) { Fr::mul_by_3(*ret, *a); }
void blst_fr_lshift(blst_fr* ret, const blst_fr* a, size_t count) { Fr::lshift(*ret, *a, count); }
void blst_fr_rshift(blst_fr* ret, const blst_fr* a, size_t count) { Fr::rshift(*ret, *a, count); }
void blst_fr_mul(blst_fr* ret, const blst_fr* a, const blst_fr* b) { Fr::mul(*ret, *a, *b); }
void blst_fr_pow(blst_fr* ret, const blst_fr* a, const byte* exponent, size_t nbits)
{
    if (nbits == 0) {
        *ret = Fr::one();
        return;
    }
    Fr::exp(*ret, *a, exponent, nbits);
}
void blst_fr_ct_bfly(blst_fr* x0, blst_fr* x1, const blst_fr* twiddle) { Fr::ct_bfly(*x0, *x1, *twiddle); }
void blst_fr_gs_bfly(blst_fr* x0, blst_fr* x1, const blst_fr* twiddle) { Fr::gs_bfly(*x0, *x1, *twiddle); }
void blst_fr_sqr(blst_fr* ret, const blst_fr* a) { Fr::sqr(*ret, *a); }
void blst_fr_cneg(blst_fr* ret, const blst_fr* a, bool flag) { Fr::cneg(*ret, *a, flag ? 1 : 0); }
void blst_fr_eucl_inverse(blst_fr* ret, const blst_fr* a) { Fr::reciprocal(*ret, *a); }
void blst_fr_inverse(blst_fr* ret, const blst_fr* a) { Fr::reciprocal(*ret, *a); }
void blst_fr_to(blst_fr* ret, const blst_fr* a) { Fr::to_mont(*ret, *a); }
void blst_fr_from(blst_fr* ret, const blst_fr* a) { Fr::from_mont(*ret, *a); }
void blst_fr_from_uint64(blst_fr* ret, const uint64_t a[4])
{
    Fr::from_uint64(*ret, a);
}
void blst_fr_from_scalar(blst_fr* ret, const blst_scalar* a)
{
    if (static_cast<const void*>(ret) == static_cast<const void*>(a)) {
        Fr::to_mont(*ret, *ret);
        return;
    }
    Fr::from_scalar(*ret, a->b);
}

void blst_scalar_from_fr(blst_scalar* ret, const blst_fr* a)
{
    if (static_cast<const void*>(ret) == static_cast<const void*>(a)) {
        Fr in;
        memcpy(&in, ret, sizeof(in));
        Fr::to_scaler(ret->b, in);
        vec_zero(&in, sizeof(in));
        return;
    }
    Fr::to_scaler(ret->b, *a);
}
bool blst_scalar_fr_check(const blst_scalar* a) {
    Fr aa;
    limbs_from_le_bytes(aa.l, a->b, sizeof(a->b));
    return Fr::is_redc(aa) != 0;
}
bool blst_sk_check(const blst_scalar* a) { return sk_check(a->b) != 0; }
bool blst_sk_add_n_check(blst_scalar* ret, const blst_scalar* a, const blst_scalar* b) { return sk_add_n_check(ret->b, a->b, b->b) != 0; }
bool blst_sk_sub_n_check(blst_scalar* ret, const blst_scalar* a, const blst_scalar* b) { return sk_sub_n_check(ret->b, a->b, b->b) != 0; }
bool blst_sk_mul_n_check(blst_scalar* ret, const blst_scalar* a, const blst_scalar* b) { return sk_mul_n_check(ret->b, a->b, b->b) != 0; }
void blst_sk_inverse(blst_scalar* ret, const blst_scalar* a) { sk_inverse(ret->b, a->b); }




/*
 * BLS12-381-specific Fp shortcuts to assembly.
 */
void blst_fp_add(blst_fp* ret, const blst_fp* a, const blst_fp* b)
{   Fp::add(*ret, *a, *b);   }

void blst_fp_sub(blst_fp* ret, const blst_fp* a, const blst_fp* b)
{   Fp::sub(*ret, *a, *b);   }

void blst_fp_mul_by_3(blst_fp* ret, const blst_fp* a)
{ Fp::mul_by_3(*ret, *a);   }

void blst_fp_mul_by_8(blst_fp* ret, const blst_fp* a)
{  Fp::mul_by_8(*ret, *a);   }

void blst_fp_lshift(blst_fp* ret, const blst_fp* a, size_t count)
{   Fp::lshift(*ret, *a, count);   }

void blst_fp_mul(blst_fp* ret, const blst_fp* a, const blst_fp* b)
{   Fp::mul(*ret, *a, *b);   }

void blst_fp_sqr(blst_fp* ret, const blst_fp* a)
{   Fp::sqr(*ret, *a);   }

void blst_fp_cneg(blst_fp* ret, const blst_fp* a, bool flag)
{   Fp::cneg(*ret, *a, flag ? 1 : 0);   }

void blst_fp_eucl_inverse(blst_fp* ret, const blst_fp* a)
{   Fp::reciprocal(*ret, *a);   }

void blst_fp_inverse(blst_fp* ret, const blst_fp* a)
{   Fp::flt_reciprocal(*ret, *a);   }

bool blst_fp_sqrt(blst_fp* ret, const blst_fp* a)
{
    Fp candidate;
    const bool_t is_root = Fp::sqrt(candidate, *a);
    vec_select(ret->l, candidate.l, Fp::zero().l, sizeof(*ret), is_root);
    vec_zero(&candidate, sizeof(candidate));
    return is_root != 0;
}

bool blst_fp_is_square(const blst_fp* a)
{   return Fp::is_sqrt(*a) != 0;   }

void blst_fp_to(blst_fp* ret, const blst_fp* a)
{   Fp::mul(*ret, *a, get_BLST_RR());   }

void blst_fp_from(blst_fp* ret, const blst_fp* a)
{   Fp::from_mont(*ret, *a);   }

/*
 * Fp serialization/deserialization.
 */
void blst_fp_from_uint32(blst_fp* ret, const uint32_t a[12])
{
    Fp in = {};

    for (int i = 0; i < 6; i++) {
        in[i] = static_cast<limb_t>(a[2 * i]) |
                (static_cast<limb_t>(a[2 * i + 1]) << 32);
    }
    Fp::mul(*ret, in, get_BLST_RR());
}

void blst_uint32_from_fp(uint32_t ret[12], const blst_fp* a)
{
    Fp out = {};

    Fp::from_mont(out, *a);
    for (int i = 0; i < 6; i++) {
        limb_t limb = out[i];
        ret[2 * i] = static_cast<uint32_t>(limb);
        ret[2 * i + 1] = static_cast<uint32_t>(limb >> 32);
    }
}

void blst_fp_from_uint64(blst_fp* ret, const uint64_t a[6])
{
    Fp in = {};

    for (int i = 0; i < 6; ++i) {
        in[i] = static_cast<limb_t>(a[i]);
    }
    Fp::mul(*ret, in, get_BLST_RR());
}

void blst_uint64_from_fp(uint64_t ret[6], const blst_fp* a)
{
    Fp out = {};

    Fp::from_mont(out, *a);
    for (int i = 0; i < 6; ++i) {
        ret[i] = static_cast<uint64_t>(out[i]);
    }
}

void blst_fp_from_bendian(blst_fp* ret, const byte a[48])
{
    Fp out;

    limbs_from_be_bytes(out.l, a, sizeof(Fp));
    Fp::mul(*ret, out, get_BLST_RR());
}

void blst_bendian_from_fp(byte ret[48], const blst_fp* a)
{
    Fp out;

    Fp::from_mont(out, *a);
    be_bytes_from_limbs(ret, out.l, sizeof(Fp));
}

void blst_fp_from_lendian(blst_fp* ret, const byte a[48])
{
    Fp out;

    limbs_from_le_bytes(out.l, a, sizeof(Fp));
    Fp::mul(*ret, out, get_BLST_RR());
}

static void fp_add_unreduced(Fp& ret, const Fp& a, const Fp& b)
{
    add_mod_v384(ret.l, a.l, b.l, get_BLST_P().l);
}

void blst_fp_from_le_bytes(blst_fp* out, const byte* bytes, size_t n)
{
    constexpr size_t chunk_size = sizeof(Fp);
    static_assert(chunk_size == 48);

    vec_zero(out->l, sizeof(*out));
    if (n == 0) {
        return;
    }

    const size_t rem = (n - 1) % chunk_size + 1;
    Fp digit;

    n -= rem;
    limbs_from_le_bytes(out->l, bytes + n, rem);
    Fp::mul(*out, get_BLST_RR(), *out);

    while (n != 0) {
        n -= chunk_size;
        limbs_from_le_bytes(digit.l, bytes + n, chunk_size);
        fp_add_unreduced(*out, *out, digit);
        Fp::mul(*out, get_BLST_RR(), *out);
    }
}

void blst_fp_from_be_bytes(blst_fp* out, const byte* bytes, size_t n)
{
    constexpr size_t chunk_size = sizeof(Fp);
    static_assert(chunk_size == 48);

    vec_zero(out->l, sizeof(*out));
    if (n == 0) {
        return;
    }

    size_t rem = (n - 1) % chunk_size + 1;
    Fp digit;

    limbs_from_be_bytes(out->l, bytes, rem);
    Fp::mul(*out, get_BLST_RR(), *out);

    while (n -= rem) {
        bytes += rem;
        limbs_from_be_bytes(digit.l, bytes, chunk_size);
        fp_add_unreduced(*out, *out, digit);
        Fp::mul(*out, get_BLST_RR(), *out);
        rem = chunk_size;
    }
}

void blst_lendian_from_fp(byte ret[48], const blst_fp* a)
{
    Fp out;

    Fp::from_mont(out, *a);
    le_bytes_from_limbs(ret, out.l, sizeof(Fp));
}

/*
 * BLS12-381-specific Fp2 shortcuts to assembly.
 */
void blst_fp2_add(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b)
{   Fp2::add(*ret, *a, *b);   }

void blst_fp2_sub(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b)
{   Fp2::sub(*ret, *a, *b);   }

void blst_fp2_mul_by_3(blst_fp2* ret, const blst_fp2* a)
{   Fp2::mul_by_3(*ret, *a);   }

void blst_fp2_mul_by_8(blst_fp2* ret, const blst_fp2* a)
{   Fp2::mul_by_8(*ret, *a);   }

void blst_fp2_lshift(blst_fp2* ret, const blst_fp2* a, size_t count)
{   Fp2::lshift(*ret, *a, count);    }

void blst_fp2_mul(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b)
{   Fp2::mul(*ret, *a, *b);   }

void blst_fp2_sqr(blst_fp2* ret, const blst_fp2* a)
{   Fp2::sqr(*ret, *a);   }

void blst_fp2_cneg(blst_fp2* ret, const blst_fp2* a, bool flag)
{   Fp2::cneg(*ret, *a, flag ? 1 : 0);   }

void blst_fp2_eucl_inverse(blst_fp2* ret, const blst_fp2* a)
{   Fp2::reciprocal(*ret, *a);   }

void blst_fp2_inverse(blst_fp2* ret, const blst_fp2* a)
{   Fp2::flt_reciprocal(*ret, *a);   }

bool blst_fp2_sqrt(blst_fp2* ret, const blst_fp2* a)
{   return Fp2::sqrt(*ret, *a) != 0;   }

bool blst_fp2_is_square(const blst_fp2* a)
{   return Fp2::is_square(*a) != 0;   }

/*
 * Scalar serialization/deserialization.
 */
void blst_scalar_from_uint32(blst_scalar* ret, const uint32_t a[8])
{
    if (static_cast<const void*>(ret) == static_cast<const void*>(a))
        return;
    memcpy(ret, a, sizeof(*ret));
}

void blst_uint32_from_scalar(uint32_t ret[8], const blst_scalar* a)
{
    if (static_cast<const void*>(ret) == static_cast<const void*>(a))
        return;
    memcpy(ret, a, sizeof(*a));
}

void blst_scalar_from_uint64(blst_scalar* ret, const uint64_t a[4])
{
    if (static_cast<const void*>(ret) == static_cast<const void*>(a))
        return;
    memcpy(ret, a, sizeof(*ret));
}

void blst_uint64_from_scalar(uint64_t ret[4], const blst_scalar* a)
{
    if (static_cast<const void*>(ret) == static_cast<const void*>(a))
        return;
    memcpy(ret, a, sizeof(*a));
}

void blst_scalar_from_bendian(blst_scalar* ret, const byte a[32])
{
    Fr out;
    limbs_from_be_bytes(out.l, a, sizeof(out));
    le_bytes_from_limbs(ret->b, out.l, sizeof(out));
    vec_zero(out.l, sizeof(out));
}

void blst_bendian_from_scalar(byte ret[32], const blst_scalar* a)
{
    Fr out;
    limbs_from_le_bytes(out.l, a->b, sizeof(out));
    be_bytes_from_limbs(ret, out.l, sizeof(out));
    vec_zero(out.l, sizeof(out));
}

void blst_scalar_from_lendian(blst_scalar* ret, const byte a[32])
{
    if (static_cast<const void*>(ret) == static_cast<const void*>(a))
        return;
    memcpy(ret, a, sizeof(*ret));
}

void blst_lendian_from_scalar(byte ret[32], const blst_scalar* a)
{
    if (static_cast<const void*>(ret) == static_cast<const void*>(a))
        return;
    memcpy(ret, a, sizeof(*a));
}


void blst_uint64_from_fr(uint64_t ret[4], const blst_fr* a)
{
    Fr out = {};

    Fr::from_mont(out, *a);
    for (int i = 0; i < 4; ++i) {
        ret[i] = static_cast<uint64_t>(out[i]);
    }
    vec_zero(out.l, sizeof(out));
}

bool blst_scalar_from_le_bytes(blst_scalar* out, const byte* bytes, size_t n)
{
    if (n == 0) {
        vec_zero(out->b, sizeof(out->b));
        return false;
    }

    size_t rem = (n - 1) % 32 + 1;
    struct { Fr out, digit; } t;
    limb_t ret;

    vec_zero(t.out.l, sizeof(t.out));

    n -= rem;
    limbs_from_le_bytes(t.out.l, bytes += n, rem);
    Fr::mul(t.out, get_BLST_rRR(), t.out);

    while (n) {
        limbs_from_le_bytes(t.digit.l, bytes -= 32, 32);
        Fr::add(t.out, t.out, t.digit);
        Fr::mul(t.out, get_BLST_rRR(), t.out);
        n -= 32;
    }

    Fr::from_mont(t.out, t.out);

    ret = vec_is_zero(t.out.l, sizeof(t.out));
    le_bytes_from_limbs(out->b, t.out.l, 32);
    vec_zero(&t, sizeof(t));

    return (ret ^ 1) != 0;
}

bool blst_scalar_from_be_bytes(blst_scalar* out, const byte* bytes, size_t n)
{
    if (n == 0) {
        vec_zero(out->b, sizeof(out->b));
        return false;
    }

    size_t rem = (n - 1) % 32 + 1;
    struct { Fr out, digit; } t;
    limb_t ret;

    vec_zero(t.out.l, sizeof(t.out));

    limbs_from_be_bytes(t.out.l, bytes, rem);
    Fr::mul(t.out, get_BLST_rRR(), t.out);

    while (n -= rem) {
        limbs_from_be_bytes(t.digit.l, bytes += rem, 32);
        Fr::add(t.out, t.out, t.digit);
        Fr::mul(t.out, get_BLST_rRR(), t.out);
        rem = 32;
    }

    Fr::from_mont(t.out, t.out);

    ret = vec_is_zero(t.out.l, sizeof(t.out));
    le_bytes_from_limbs(out->b, t.out.l, 32);
    vec_zero(&t, sizeof(t));

    return (ret ^ 1) != 0;
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
bool blst_scalar_from_hexascii(blst_scalar* ret, const byte* hex)
{
    return bytes_from_hexascii(ret->b, sizeof(ret->b),
                               reinterpret_cast<const char *>(hex));
}

bool blst_fr_from_hexascii(blst_fr* ret, const byte* hex)
{
    if (ret == nullptr)
        return false;
    return Fr::from_hexascii(*ret, reinterpret_cast<const char *>(hex));
}

bool blst_fp_from_hexascii(blst_fp* ret, const byte* hex)
{
    if (ret == nullptr)
        return false;

    Fp parsed = {};
    if (!limbs_from_hexascii(parsed.l, sizeof(Fp),
                             reinterpret_cast<const char *>(hex)))
        return false;
    Fp::mul(parsed, parsed, get_BLST_RR());
    *ret = parsed;
    return true;
}
BLST_EXTERN_C_END_
