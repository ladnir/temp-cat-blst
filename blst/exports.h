//#include "blst/common/defines.h"
//BLST_EXTERN_C_BEGIN_
//
///*
// * BLS12-381-specific Fr shortcuts to assembly.
// */
//void blst_fr_add(Fr& ret, const Fr& a, const Fr& b);
//void blst_fr_sub(Fr& ret, const Fr& a, const Fr& b);
//void blst_fr_mul_by_3(Fr& ret, const Fr& a);
//void blst_fr_lshift(Fr& ret, const Fr& a, size_t count);
//void blst_fr_rshift(Fr& ret, const Fr& a, size_t count);
//void blst_fr_mul(Fr& ret, const Fr& a, const Fr& b);
//void blst_fr_ct_bfly(Fr& x0, Fr& x1, const Fr& twiddle);
//void blst_fr_gs_bfly(Fr& x0, Fr& x1, const Fr& twiddle);
//void blst_fr_sqr(Fr& ret, const Fr& a);
//void blst_fr_cneg(Fr& ret, const Fr& a, int flag);
//void blst_fr_to(Fr& ret, const Fr& a);
//void blst_fr_from(Fr& ret, const Fr& a);
//void blst_fr_from_scalar(Fr& ret, const pow256 a);
//void blst_scalar_from_fr(pow256 ret, const Fr& a);
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
//void blst_fp_add(Fp& ret, const Fp& a, const Fp& b);
//void blst_fp_sub(Fp& ret, const Fp& a, const Fp& b);
//void blst_fp_mul_by_3(Fp& ret, const Fp& a);
//void blst_fp_mul_by_8(Fp& ret, const Fp& a);
//void blst_fp_lshift(Fp& ret, const Fp& a, size_t count);
//void blst_fp_mul(Fp& ret, const Fp& a, const Fp& b);
//void blst_fp_sqr(Fp& ret, const Fp& a);
//void blst_fp_cneg(Fp& ret, const Fp& a, int flag);
//void blst_fp_to(Fp& ret, const Fp& a);
//void blst_fp_from(Fp& ret, const Fp& a);
//
///*
// * Fp serialization/deserialization.
// */
//void blst_fp_from_uint32(Fp& ret, const unsigned int a[12]);
//void blst_uint32_from_fp(unsigned int ret[12], const Fp& a);
//void blst_fp_from_uint64(Fp& ret, const unsigned long long a[6]);
//void blst_uint64_from_fp(unsigned long long ret[6], const Fp& a);
//void blst_fp_from_bendian(Fp& ret, const unsigned char a[48]);
//void blst_bendian_from_fp(unsigned char ret[48], const Fp& a);
//void blst_fp_from_lendian(Fp& ret, const unsigned char a[48]);
//void blst_lendian_from_fp(unsigned char ret[48], const Fp& a);
//
///*
// * BLS12-381-specific Fp2 shortcuts to assembly.
// */
//void blst_fp2_add(Fp2& ret, const Fp2& a, const Fp2& b);
//void blst_fp2_sub(Fp2& ret, const Fp2& a, const Fp2& b);
//void blst_fp2_mul_by_3(Fp2& ret, const Fp2& a);
//void blst_fp2_mul_by_8(Fp2& ret, const Fp2& a);
//void blst_fp2_lshift(Fp2& ret, const Fp2& a, size_t count);
//void blst_fp2_mul(Fp2& ret, const Fp2& a, const Fp2& b);
//void blst_fp2_sqr(Fp2& ret, const Fp2& a);
//void blst_fp2_cneg(Fp2& ret, const Fp2& a, int flag);
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
//void blst_fr_from_uint64(Fr& ret, const unsigned long long a[4]);
//void blst_uint64_from_fr(unsigned long long ret[4], const Fr& a);
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
//void blst_fr_from_hexascii(Fr& ret, const char* hex);
//void blst_fp_from_hexascii(Fp& ret, const char* hex);
//BLST_EXTERN_C_END_
