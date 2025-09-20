//#include "blst/common/defines.h"
//BLST_EXTERN_C_BEGIN_
//
///*
// * BLS12-381-specific Fr shortcuts to assembly.
// */
//void blst_fr_add(vec256 ret, const vec256 a, const vec256 b);
//void blst_fr_sub(vec256 ret, const vec256 a, const vec256 b);
//void blst_fr_mul_by_3(vec256 ret, const vec256 a);
//void blst_fr_lshift(vec256 ret, const vec256 a, size_t count);
//void blst_fr_rshift(vec256 ret, const vec256 a, size_t count);
//void blst_fr_mul(vec256 ret, const vec256 a, const vec256 b);
//void blst_fr_ct_bfly(vec256 x0, vec256 x1, const vec256 twiddle);
//void blst_fr_gs_bfly(vec256 x0, vec256 x1, const vec256 twiddle);
//void blst_fr_sqr(vec256 ret, const vec256 a);
//void blst_fr_cneg(vec256 ret, const vec256 a, int flag);
//void blst_fr_to(vec256 ret, const vec256 a);
//void blst_fr_from(vec256 ret, const vec256 a);
//void blst_fr_from_scalar(vec256 ret, const pow256 a);
//void blst_scalar_from_fr(pow256 ret, const vec256 a);
//int blst_scalar_fr_check(const pow256 a);
//int blst_sk_check(const pow256 a);
//int blst_sk_add_n_check(pow256 ret, const pow256 a, const pow256 b);
//int blst_sk_sub_n_check(pow256 ret, const pow256 a, const pow256 b);
//int blst_sk_mul_n_check(pow256 ret, const pow256 a, const pow256 b);
//void blst_sk_inverse(pow256 ret, const pow256 a);
//
///*
// * BLS12-381-specific Fp shortcuts to assembly.
// */
//void blst_fp_add(vec384 ret, const vec384 a, const vec384 b);
//void blst_fp_sub(vec384 ret, const vec384 a, const vec384 b);
//void blst_fp_mul_by_3(vec384 ret, const vec384 a);
//void blst_fp_mul_by_8(vec384 ret, const vec384 a);
//void blst_fp_lshift(vec384 ret, const vec384 a, size_t count);
//void blst_fp_mul(vec384 ret, const vec384 a, const vec384 b);
//void blst_fp_sqr(vec384 ret, const vec384 a);
//void blst_fp_cneg(vec384 ret, const vec384 a, int flag);
//void blst_fp_to(vec384 ret, const vec384 a);
//void blst_fp_from(vec384 ret, const vec384 a);
//
///*
// * Fp serialization/deserialization.
// */
//void blst_fp_from_uint32(vec384 ret, const unsigned int a[12]);
//void blst_uint32_from_fp(unsigned int ret[12], const vec384 a);
//void blst_fp_from_uint64(vec384 ret, const unsigned long long a[6]);
//void blst_uint64_from_fp(unsigned long long ret[6], const vec384 a);
//void blst_fp_from_bendian(vec384 ret, const unsigned char a[48]);
//void blst_bendian_from_fp(unsigned char ret[48], const vec384 a);
//void blst_fp_from_lendian(vec384 ret, const unsigned char a[48]);
//void blst_lendian_from_fp(unsigned char ret[48], const vec384 a);
//
///*
// * BLS12-381-specific Fp2 shortcuts to assembly.
// */
//void blst_fp2_add(vec384x ret, const vec384x a, const vec384x b);
//void blst_fp2_sub(vec384x ret, const vec384x a, const vec384x b);
//void blst_fp2_mul_by_3(vec384x ret, const vec384x a);
//void blst_fp2_mul_by_8(vec384x ret, const vec384x a);
//void blst_fp2_lshift(vec384x ret, const vec384x a, size_t count);
//void blst_fp2_mul(vec384x ret, const vec384x a, const vec384x b);
//void blst_fp2_sqr(vec384x ret, const vec384x a);
//void blst_fp2_cneg(vec384x ret, const vec384x a, int flag);
//
///*
// * Scalar serialization/deserialization.
// */
//void blst_scalar_from_uint32(pow256 ret, const unsigned int a[8]);
//void blst_uint32_from_scalar(unsigned int ret[8], const pow256 a);
//void blst_scalar_from_uint64(pow256 ret, const unsigned long long a[4]);
//void blst_uint64_from_scalar(unsigned long long ret[4], const pow256 a);
//void blst_scalar_from_bendian(pow256 ret, const unsigned char a[32]);
//void blst_bendian_from_scalar(unsigned char ret[32], const pow256 a);
//void blst_scalar_from_lendian(pow256 ret, const unsigned char a[32]);
//void blst_lendian_from_scalar(unsigned char ret[32], const pow256 a);
//void blst_fr_from_uint64(vec256 ret, const unsigned long long a[4]);
//void blst_uint64_from_fr(unsigned long long ret[4], const vec256 a);
//int blst_scalar_from_le_bytes(pow256 out, const unsigned char* bytes, size_t n);
//int blst_scalar_from_be_bytes(pow256 out, const unsigned char* bytes, size_t n);
//
///*
// * Single-short SHA-256 hash function.
// */
//void blst_sha256(unsigned char md[32], const void* msg, size_t len);
//
///*
// * Test facilitator.
// */
//void blst_scalar_from_hexascii(pow256 ret, const char* hex);
//void blst_fr_from_hexascii(vec256 ret, const char* hex);
//void blst_fp_from_hexascii(vec384 ret, const char* hex);
//BLST_EXTERN_C_END_
