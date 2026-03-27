///*
// * Copyright Supranational LLC
// * Licensed under the Apache License, Version 2.0, see LICENSE for details.
// * SPDX-License-Identifier: Apache-2.0
// */
//#ifndef __BLST_AUX_H__
//#define __BLST_AUX_H__
///*
// * This file lists interfaces that might be promoted to blst.h or removed,
// * depending on their proven/unproven worthiness.
// */
//
//void blst_fr_ct_bfly(blst_fr *x0, blst_fr *x1, const blst_fr *twiddle);
//void blst_fr_gs_bfly(blst_fr *x0, blst_fr *x1, const blst_fr *twiddle);
//void blst_fr_to(blst_fr *ret, const blst_fr *a);
//void blst_fr_from(blst_fr *ret, const blst_fr *a);
//#ifdef BLST_FR_PENTAROOT
//void blst_fr_pentaroot(blst_fr *ret, const blst_fr *a);
//void blst_fr_pentapow(blst_fr *ret, const blst_fr *a);
//#endif
//
//void blst_fp_to(blst_fp *ret, const blst_fp *a);
//void blst_fp_from(blst_fp *ret, const blst_fp *a);
//
//bool blst_fp_is_square(const blst_fp *a);
//bool blst_fp2_is_square(const blst_fp2 *a);
//
//void blst_p1_from_jacobian(blst_p1 *out, const blst_p1 *in);
//void blst_p2_from_jacobian(blst_p2 *out, const blst_p2 *in);
//
///*
// * Below functions produce both point and deserialized outcome of
// * SkToPk and Sign. However, deserialized outputs are pre-decorated
// * with sign and infinity bits. This means that you have to bring the
// * output into compliance prior returning to application. If you want
// * compressed point value, then do [equivalent of]
// *
// *  byte temp[96];
// *  blst_sk_to_pk2_in_g1(temp, out_pk, SK);
// *  temp[0] |= 0x80;
// *  memcpy(out, temp, 48);
// *
// * Otherwise do
// *
// *  blst_sk_to_pk2_in_g1(out, out_pk, SK);
// *  out[0] &= ~0x20;
// *
// * Either |out| or |out_<point>| can be NULL.
// */
//void blst_sk_to_pk2_in_g1(byte out[96], blst_p1_affine *out_pk,
//                          const blst_scalar *SK);
//void blst_sign_pk2_in_g1(byte out[192], blst_p2_affine *out_sig,
//                         const blst_p2 *hash, const blst_scalar *SK);
//void blst_sk_to_pk2_in_g2(byte out[192], blst_p2_affine *out_pk,
//                          const blst_scalar *SK);
//void blst_sign_pk2_in_g2(byte out[96], blst_p1_affine *out_sig,
//                         const blst_p1 *hash, const blst_scalar *SK);
//
//#ifdef __BLST_RUST_BINDGEN__
//typedef struct {} blst_uniq;
//#else
//typedef struct blst_opaque blst_uniq;
//#endif
//
//size_t blst_uniq_sizeof(size_t n_nodes);
//void blst_uniq_init(blst_uniq *tree);
//bool blst_uniq_test(blst_uniq *tree, const byte *msg, size_t len);
//
//#ifdef expand_message_xmd
//void expand_message_xmd(unsigned char *bytes, size_t len_in_bytes,
//                        const unsigned char *aug, size_t aug_len,
//                        const unsigned char *msg, size_t msg_len,
//                        const unsigned char *DST, size_t DST_len);
//#else
//void blst_expand_message_xmd(byte *out, size_t out_len,
//                             const byte *msg, size_t msg_len,
//                             const byte *DST, size_t DST_len);
//#endif
//
//void blst_p1_unchecked_mult(blst_p1 *out, const blst_p1 *p, const byte *scalar,
//                                                            size_t nbits);
//void blst_p2_unchecked_mult(blst_p2 *out, const blst_p2 *p, const byte *scalar,
//                                                            size_t nbits);
//
//void blst_pairing_raw_aggregate(blst_pairing *ctx, const blst_p2_affine *q,
//                                                   const blst_p1_affine *p);
//blst_fp12 blst_pairing_as_fp12(blst_pairing *ctx);
//void blst_pairing_as_fp12_to(blst_fp12 *out, blst_pairing *ctx);
//void blst_bendian_from_fp12(byte out[48*12], const blst_fp12 *a);
//
//void blst_keygen_v3(blst_scalar *out_SK, const byte *IKM, size_t IKM_len,
//                    const byte *info DEFNULL, size_t info_len DEFNULL);
//void blst_keygen_v4_5(blst_scalar *out_SK, const byte *IKM, size_t IKM_len,
//                      const byte *salt, size_t salt_len,
//                      const byte *info DEFNULL, size_t info_len DEFNULL);
//void blst_keygen_v5(blst_scalar *out_SK, const byte *IKM, size_t IKM_len,
//                    const byte *salt, size_t salt_len,
//                    const byte *info DEFNULL, size_t info_len DEFNULL);
//void blst_derive_master_eip2333(blst_scalar *out_SK,
//                                const byte *IKM, size_t IKM_len);
//void blst_derive_child_eip2333(blst_scalar *out_SK, const blst_scalar *SK,
//                               uint32_t child_index);
//
//void blst_scalar_from_hexascii(blst_scalar *out, const byte *hex);
//void blst_fr_from_hexascii(blst_fr *ret, const byte *hex);
//void blst_fp_from_hexascii(blst_fp *ret, const byte *hex);
//
//size_t blst_p1_sizeof(void);
//size_t blst_p1_affine_sizeof(void);
//size_t blst_p2_sizeof(void);
//size_t blst_p2_affine_sizeof(void);
//size_t blst_fp12_sizeof(void);
//
///*
// * Single-shot SHA-256 hash function.
// */
//void blst_sha256(byte out[32], const byte *msg, size_t msg_len);
//#endif

#ifndef __BLST_AUX_LIVE_H__
#define __BLST_AUX_LIVE_H__

#ifdef __SIZE_TYPE__
typedef __SIZE_TYPE__ size_t;
#else
#include <stddef.h>
#endif

#if defined(__UINT8_TYPE__) && defined(__UINT32_TYPE__) && defined(__UINT64_TYPE__)
typedef __UINT8_TYPE__ uint8_t;
typedef __UINT32_TYPE__ uint32_t;
typedef __UINT64_TYPE__ uint64_t;
#else
#include <stdint.h>
#endif

#ifdef __cplusplus
extern "C" {
#elif !defined(__STDC_VERSION__) || __STDC_VERSION__ < 202311L
# if !defined(bool)
#  if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#   define bool _Bool
#  else
#   define bool int
#  endif
#  define __blst_aux_h_bool__
# endif
#endif

typedef enum {
    BLST_SUCCESS = 0,
    BLST_BAD_ENCODING,
    BLST_POINT_NOT_ON_CURVE,
    BLST_POINT_NOT_IN_GROUP,
    BLST_AGGR_TYPE_MISMATCH,
    BLST_VERIFY_FAIL,
    BLST_PK_IS_INFINITY,
} BLST_ERROR;

typedef uint8_t byte;
typedef uint64_t limb_t;

typedef struct { byte b[256 / 8]; } blst_scalar;
typedef struct { limb_t l[256 / 8 / sizeof(limb_t)]; } blst_fr;
typedef struct { limb_t l[384 / 8 / sizeof(limb_t)]; } blst_fp;
typedef struct { blst_fp fp[2]; } blst_fp2;
typedef struct { blst_fp2 fp2[3]; } blst_fp6;
typedef struct { blst_fp6 fp6[2]; } blst_fp12;

typedef struct { blst_fp x, y, z; } blst_p1;
typedef struct { blst_fp x, y; } blst_p1_affine;
typedef struct { blst_fp2 x, y, z; } blst_p2;
typedef struct { blst_fp2 x, y; } blst_p2_affine;

typedef struct blst_middle_product_precomputed_b blst_middle_product_precomputed_b;

void blst_fr_add(blst_fr* ret, const blst_fr* a, const blst_fr* b);
void blst_fr_sub(blst_fr* ret, const blst_fr* a, const blst_fr* b);
void blst_fr_mul_by_3(blst_fr* ret, const blst_fr* a);
void blst_fr_lshift(blst_fr* ret, const blst_fr* a, size_t count);
void blst_fr_rshift(blst_fr* ret, const blst_fr* a, size_t count);
void blst_fr_mul(blst_fr* ret, const blst_fr* a, const blst_fr* b);
void blst_fr_sqr(blst_fr* ret, const blst_fr* a);
void blst_fr_cneg(blst_fr* ret, const blst_fr* a, bool flag);
void blst_fr_eucl_inverse(blst_fr* ret, const blst_fr* a);
void blst_fr_inverse(blst_fr* ret, const blst_fr* a);
void blst_fr_from_uint64(blst_fr* ret, const uint64_t a[4]);
void blst_uint64_from_fr(uint64_t ret[4], const blst_fr* a);
void blst_fr_from_scalar(blst_fr* ret, const blst_scalar* a);
void blst_scalar_from_fr(blst_scalar* ret, const blst_fr* a);
void blst_fr_ct_bfly(blst_fr* x0, blst_fr* x1, const blst_fr* twiddle);
void blst_fr_gs_bfly(blst_fr* x0, blst_fr* x1, const blst_fr* twiddle);
void blst_fr_to(blst_fr* ret, const blst_fr* a);
void blst_fr_from(blst_fr* ret, const blst_fr* a);
void blst_fr_from_hexascii(blst_fr* ret, const byte* hex);

void blst_fp_add(blst_fp* ret, const blst_fp* a, const blst_fp* b);
void blst_fp_sub(blst_fp* ret, const blst_fp* a, const blst_fp* b);
void blst_fp_mul_by_3(blst_fp* ret, const blst_fp* a);
void blst_fp_mul_by_8(blst_fp* ret, const blst_fp* a);
void blst_fp_lshift(blst_fp* ret, const blst_fp* a, size_t count);
void blst_fp_mul(blst_fp* ret, const blst_fp* a, const blst_fp* b);
void blst_fp_sqr(blst_fp* ret, const blst_fp* a);
void blst_fp_cneg(blst_fp* ret, const blst_fp* a, bool flag);
void blst_fp_eucl_inverse(blst_fp* ret, const blst_fp* a);
void blst_fp_inverse(blst_fp* ret, const blst_fp* a);
bool blst_fp_sqrt(blst_fp* ret, const blst_fp* a);
void blst_fp_from_uint32(blst_fp* ret, const uint32_t a[12]);
void blst_uint32_from_fp(uint32_t ret[12], const blst_fp* a);
void blst_fp_from_uint64(blst_fp* ret, const uint64_t a[6]);
void blst_uint64_from_fp(uint64_t ret[6], const blst_fp* a);
void blst_fp_from_bendian(blst_fp* ret, const byte a[48]);
void blst_bendian_from_fp(byte ret[48], const blst_fp* a);
void blst_fp_from_lendian(blst_fp* ret, const byte a[48]);
void blst_lendian_from_fp(byte ret[48], const blst_fp* a);
void blst_fp_to(blst_fp* ret, const blst_fp* a);
void blst_fp_from(blst_fp* ret, const blst_fp* a);
bool blst_fp_is_square(const blst_fp* a);
void blst_fp_from_hexascii(blst_fp* ret, const byte* hex);

void blst_fp2_add(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b);
void blst_fp2_sub(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b);
void blst_fp2_mul_by_3(blst_fp2* ret, const blst_fp2* a);
void blst_fp2_mul_by_8(blst_fp2* ret, const blst_fp2* a);
void blst_fp2_lshift(blst_fp2* ret, const blst_fp2* a, size_t count);
void blst_fp2_mul(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b);
void blst_fp2_sqr(blst_fp2* ret, const blst_fp2* a);
void blst_fp2_cneg(blst_fp2* ret, const blst_fp2* a, bool flag);
void blst_fp2_eucl_inverse(blst_fp2* ret, const blst_fp2* a);
void blst_fp2_inverse(blst_fp2* ret, const blst_fp2* a);
bool blst_fp2_sqrt(blst_fp2* ret, const blst_fp2* a);
bool blst_fp2_is_square(const blst_fp2* a);

void blst_fp12_mul(blst_fp12* ret, const blst_fp12* a, const blst_fp12* b);
bool blst_fp12_is_equal(const blst_fp12* a, const blst_fp12* b);
bool blst_fp12_is_one(const blst_fp12* a);
bool blst_fp12_in_group(const blst_fp12* a);
const blst_fp12* blst_fp12_one(void);
void blst_final_exp(blst_fp12* ret, const blst_fp12* f);
bool blst_fp12_finalverify(const blst_fp12* gt1, const blst_fp12* gt2);
void blst_bendian_from_fp12(byte out[48 * 12], const blst_fp12* a);

void blst_p1_add(blst_p1* out, const blst_p1* a, const blst_p1* b);
void blst_p1_add_or_double(blst_p1* out, const blst_p1* a, const blst_p1* b);
void blst_p1_add_affine(blst_p1* out, const blst_p1* a, const blst_p1_affine* b);
void blst_p1_add_or_double_affine(blst_p1* out, const blst_p1* a, const blst_p1_affine* b);
void blst_p1_double(blst_p1* out, const blst_p1* a);
void blst_p1_mult(blst_p1* out, const blst_p1* p, const byte* scalar, size_t nbits);
void blst_p1_cneg(blst_p1* p, bool cbit);
void blst_p1_to_affine(blst_p1_affine* out, const blst_p1* in);
void blst_p1_from_affine(blst_p1* out, const blst_p1_affine* in);
bool blst_p1_on_curve(const blst_p1* p);
bool blst_p1_in_g1(const blst_p1* p);
bool blst_p1_is_equal(const blst_p1* a, const blst_p1* b);
bool blst_p1_is_inf(const blst_p1* a);
const blst_p1* blst_p1_generator(void);
bool blst_p1_affine_on_curve(const blst_p1_affine* p);
bool blst_p1_affine_in_g1(const blst_p1_affine* p);
bool blst_p1_affine_is_equal(const blst_p1_affine* a, const blst_p1_affine* b);
bool blst_p1_affine_is_inf(const blst_p1_affine* a);
const blst_p1_affine* blst_p1_affine_generator(void);
void blst_p1_serialize(byte out[96], const blst_p1* in);
void blst_p1_compress(byte out[48], const blst_p1* in);
void blst_p1_affine_serialize(byte out[96], const blst_p1_affine* in);
void blst_p1_affine_compress(byte out[48], const blst_p1_affine* in);
BLST_ERROR blst_p1_uncompress(blst_p1_affine* out, const byte in[48]);
BLST_ERROR blst_p1_deserialize(blst_p1_affine* out, const byte in[96]);

void blst_p2_add(blst_p2* out, const blst_p2* a, const blst_p2* b);
void blst_p2_add_or_double(blst_p2* out, const blst_p2* a, const blst_p2* b);
void blst_p2_add_affine(blst_p2* out, const blst_p2* a, const blst_p2_affine* b);
void blst_p2_add_or_double_affine(blst_p2* out, const blst_p2* a, const blst_p2_affine* b);
void blst_p2_double(blst_p2* out, const blst_p2* a);
void blst_p2_mult(blst_p2* out, const blst_p2* p, const byte* scalar, size_t nbits);
void blst_p2_cneg(blst_p2* p, bool cbit);
void blst_p2_to_affine(blst_p2_affine* out, const blst_p2* in);
void blst_p2_from_affine(blst_p2* out, const blst_p2_affine* in);
bool blst_p2_on_curve(const blst_p2* p);
bool blst_p2_in_g2(const blst_p2* p);
bool blst_p2_is_equal(const blst_p2* a, const blst_p2* b);
bool blst_p2_is_inf(const blst_p2* a);
const blst_p2* blst_p2_generator(void);
bool blst_p2_affine_on_curve(const blst_p2_affine* p);
bool blst_p2_affine_in_g2(const blst_p2_affine* p);
bool blst_p2_affine_is_equal(const blst_p2_affine* a, const blst_p2_affine* b);
bool blst_p2_affine_is_inf(const blst_p2_affine* a);
const blst_p2_affine* blst_p2_affine_generator(void);
void blst_p2_serialize(byte out[192], const blst_p2* in);
void blst_p2_compress(byte out[96], const blst_p2* in);
void blst_p2_affine_serialize(byte out[192], const blst_p2_affine* in);
void blst_p2_affine_compress(byte out[96], const blst_p2_affine* in);
BLST_ERROR blst_p2_uncompress(blst_p2_affine* out, const byte in[96]);
BLST_ERROR blst_p2_deserialize(blst_p2_affine* out, const byte in[192]);

void blst_miller_loop(blst_fp12* ret, const blst_p2_affine* Q, const blst_p1_affine* P);
void blst_miller_loop_n(blst_fp12* ret, const blst_p2_affine* const Qs[], const blst_p1_affine* const Ps[], size_t n);
void blst_precompute_lines(blst_fp6 Qlines[68], const blst_p2_affine* Q);
void blst_miller_loop_lines(blst_fp12* ret, const blst_fp6 Qlines[68], const blst_p1_affine* P);

blst_middle_product_precomputed_b* blst_middle_product_precompute_b(const blst_p2* b, size_t n);
bool blst_middle_product(const blst_p1* a, size_t n, const blst_middle_product_precomputed_b* precomputed_b, blst_fp12* out);
void blst_middle_product_precomputed_b_free(blst_middle_product_precomputed_b* precomputed_b);

#ifdef __cplusplus
}
#elif defined(__blst_aux_h_bool__)
# undef bool
# undef __blst_aux_h_bool__
#endif

#endif
