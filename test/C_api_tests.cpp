#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

#include "bindings/blst.h"
#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fr.h"
#include "blst/scalar/fields/hash_to_field.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/map_to_G.h"
#include "blst/scalar/curve/multi_scalar.h"
#include "blst/scalar/curve/point.h"
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/map_to_G.h"
#endif

static_assert(std::is_same_v<blst::BLST_ERROR, ::BLST_ERROR>);
static_assert(static_cast<int>(BLST_SUCCESS) == 0);
static_assert(static_cast<int>(BLST_BAD_SCALAR) == 7);

using BlstAvxP1 = blst_avx_p1;
using BlstAvxP1Affine = blst_avx_p1_affine;
using BlstAvxP2 = blst_avx_p2;
using BlstAvxP2Affine = blst_avx_p2_affine;
using BlstAvxFp12 = blst_avx_fp12;

// Intentionally restate the C surface with internal aliases. Because the
// public umbrella header is included above, the compiler rejects any function
// type that is not exactly compatible with the supported public declaration.
extern "C" {
bool blst_scalar_from_le_bytes(blst_scalar* out, const blst_byte* in, size_t len);
bool blst_scalar_from_be_bytes(blst_scalar* out, const blst_byte* in, size_t len);
bool blst_sk_add_n_check(blst_scalar* out, const blst_scalar* a, const blst_scalar* b);
bool blst_sk_sub_n_check(blst_scalar* out, const blst_scalar* a, const blst_scalar* b);
bool blst_sk_mul_n_check(blst_scalar* out, const blst_scalar* a, const blst_scalar* b);
void blst_sk_inverse(blst_scalar* out, const blst_scalar* a);

void blst_fr_add(blst::Fr* ret, const blst::Fr* a, const blst::Fr* b);
void blst_fr_lshift(blst::Fr* ret, const blst::Fr* a, size_t count);
void blst_fr_rshift(blst::Fr* ret, const blst::Fr* a, size_t count);
void blst_fr_mul(blst::Fr* ret, const blst::Fr* a, const blst::Fr* b);
void blst_fr_pow(blst::Fr* ret, const blst::Fr* a, const unsigned char* exponent, size_t nbits);
void blst_fr_sqr(blst::Fr* ret, const blst::Fr* a);
void blst_fr_eucl_inverse(blst::Fr* ret, const blst::Fr* a);
void blst_fr_inverse(blst::Fr* ret, const blst::Fr* a);
void blst_fr_from_uint64(blst::Fr* ret, const uint64_t a[4]);
void blst_uint64_from_fr(uint64_t ret[4], const blst::Fr* a);
void blst_fr_from_scalar(blst::Fr* ret, const blst_scalar* a);
void blst_scalar_from_fr(blst_scalar* ret, const blst::Fr* a);
bool blst_fr_from_hexascii(blst::Fr* ret, const unsigned char* hex);

void blst_fp_add(blst::Fp* ret, const blst::Fp* a, const blst::Fp* b);
void blst_fp_lshift(blst::Fp* ret, const blst::Fp* a, size_t count);
void blst_fp_mul(blst::Fp* ret, const blst::Fp* a, const blst::Fp* b);
void blst_fp_sqr(blst::Fp* ret, const blst::Fp* a);
void blst_fp_eucl_inverse(blst::Fp* ret, const blst::Fp* a);
void blst_fp_inverse(blst::Fp* ret, const blst::Fp* a);
bool blst_fp_sqrt(blst::Fp* ret, const blst::Fp* a);
void blst_fp_from_uint32(blst::Fp* ret, const uint32_t a[12]);
void blst_uint32_from_fp(uint32_t ret[12], const blst::Fp* a);
void blst_fp_from_uint64(blst::Fp* ret, const uint64_t a[6]);
void blst_uint64_from_fp(uint64_t ret[6], const blst::Fp* a);
void blst_fp_from_bendian(blst::Fp* ret, const unsigned char a[48]);
void blst_bendian_from_fp(unsigned char ret[48], const blst::Fp* a);
bool blst_fp_is_square(const blst::Fp* a);
bool blst_fp_from_hexascii(blst::Fp* ret, const unsigned char* hex);

void blst_fp2_mul(blst::Fp2* ret, const blst::Fp2* a, const blst::Fp2* b);
void blst_fp2_lshift(blst::Fp2* ret, const blst::Fp2* a, size_t count);
void blst_fp2_sqr(blst::Fp2* ret, const blst::Fp2* a);
void blst_fp2_eucl_inverse(blst::Fp2* ret, const blst::Fp2* a);
void blst_fp2_inverse(blst::Fp2* ret, const blst::Fp2* a);
bool blst_fp2_sqrt(blst::Fp2* ret, const blst::Fp2* a);
bool blst_fp2_is_square(const blst::Fp2* a);

void blst_fp12_mul(blst::Fp12* ret, const blst::Fp12* a, const blst::Fp12* b);
void blst_fp12_sqr(blst::Fp12* ret, const blst::Fp12* a);
void blst_fp12_cyclotomic_sqr(blst::Fp12* ret, const blst::Fp12* a);
void blst_fp12_mul_by_xy00z0(blst::Fp12* ret, const blst::Fp12* a, const blst::Fp6* xy00z0);
void blst_fp12_conjugate(blst::Fp12* a);
void blst_fp12_inverse(blst::Fp12* ret, const blst::Fp12* a);
bool blst_fp12_frobenius_map(blst::Fp12* ret, const blst::Fp12* a, size_t n);
bool blst_fp12_is_equal(const blst::Fp12* a, const blst::Fp12* b);
bool blst_fp12_is_one(const blst::Fp12* a);
bool blst_fp12_in_group(const blst::Fp12* a);
const blst::Fp12* blst_fp12_one(void);
void blst_final_exp(blst::Fp12* ret, const blst::Fp12* f);
const blst::Fp12* blst_gt_identity(void);
const blst::Fp12* blst_gt_generator(void);
bool blst_gt_is_equal(const blst::Fp12* a, const blst::Fp12* b);
bool blst_gt_is_identity(const blst::Fp12* p);
void blst_gt_add(blst::Fp12* out, const blst::Fp12* a, const blst::Fp12* b);
void blst_gt_sub(blst::Fp12* out, const blst::Fp12* a, const blst::Fp12* b);
void blst_gt_neg(blst::Fp12* out, const blst::Fp12* p);
void blst_gt_double(blst::Fp12* out, const blst::Fp12* p);
void blst_gt_mult(blst::Fp12* out, const blst::Fp12* p, const unsigned char* scalar, size_t nbits);
bool blst_fp12_finalverify(const blst::Fp12* gt1, const blst::Fp12* gt2);

void blst_p1_add(blst::scalar::G1Point* out, const blst::scalar::G1Point* a, const blst::scalar::G1Point* b);
void blst_p1_sub(blst::scalar::G1Point* out, const blst::scalar::G1Point* a, const blst::scalar::G1Point* b);
void blst_p1_add_or_double(blst::scalar::G1Point* out, const blst::scalar::G1Point* a, const blst::scalar::G1Point* b);
void blst_p1_sub_affine(blst::scalar::G1Point* out, const blst::scalar::G1Point* a, const blst::scalar::G1AffinePoint* b);
void blst_p1_double(blst::scalar::G1Point* out, const blst::scalar::G1Point* a);
void blst_p1_mult(blst::scalar::G1Point* out, const blst::scalar::G1Point* p, const unsigned char* scalar, size_t nbits);
void blst_p1_unchecked_mult(blst::scalar::G1Point* out, const blst::scalar::G1Point* p, const unsigned char* scalar, size_t nbits);
void blst_p1_affine_cneg(blst::scalar::G1AffinePoint* p, bool cbit);
void blst_p1_to_affine(blst::scalar::G1AffinePoint* out, const blst::scalar::G1Point* in);
void blst_p1_from_affine(blst::scalar::G1Point* out, const blst::scalar::G1AffinePoint* in);
void blst_p1_from_jacobian(blst::scalar::G1Point* out, const blst::scalar::G1Point* in);
bool blst_p1_on_curve(const blst::scalar::G1Point* p);
bool blst_p1_in_g1(const blst::scalar::G1Point* p);
bool blst_p1_is_equal(const blst::scalar::G1Point* a, const blst::scalar::G1Point* b);
bool blst_p1_is_inf(const blst::scalar::G1Point* a);
const blst::scalar::G1Point* blst_p1_generator(void);
bool blst_p1_affine_in_g1(const blst::scalar::G1AffinePoint* p);
bool blst_p1_affine_is_equal(const blst::scalar::G1AffinePoint* a, const blst::scalar::G1AffinePoint* b);
bool blst_p1_affine_is_inf(const blst::scalar::G1AffinePoint* a);
const blst::scalar::G1AffinePoint* blst_p1_affine_generator(void);
void blst_p1s_to_affine(blst::scalar::G1AffinePoint dst[], const blst::scalar::G1Point* const points[], size_t npoints);
void blst_p1s_add(blst::scalar::G1Point* ret, const blst::scalar::G1AffinePoint* const points[], size_t npoints);
size_t blst_p1s_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints);
void blst_p1s_mult_wbits_precompute(blst::scalar::G1AffinePoint table[], size_t wbits, const blst::scalar::G1AffinePoint* const points[], size_t npoints);
size_t blst_p1s_mult_wbits_scratch_sizeof(size_t npoints);
void blst_p1s_mult_wbits(blst::scalar::G1Point* ret, const blst::scalar::G1AffinePoint table[], size_t wbits, size_t npoints, const unsigned char* const scalars[], size_t nbits, blst::limb_t* scratch);
size_t blst_p1s_mult_pippenger_scratch_sizeof(size_t npoints);
void blst_p1s_mult_pippenger(blst::scalar::G1Point* ret, const blst::scalar::G1AffinePoint* const points[], size_t npoints, const unsigned char* const scalars[], size_t nbits, blst::limb_t* scratch);
void blst_p1s_tile_pippenger(blst::scalar::G1Point* ret, const blst::scalar::G1AffinePoint* const points[], size_t npoints, const unsigned char* const scalars[], size_t nbits, blst::limb_t* scratch, size_t bit0, size_t window);
void blst_p1_affine_serialize(unsigned char out[96], const blst::scalar::G1AffinePoint* in);
void blst_p1_affine_compress(unsigned char out[48], const blst::scalar::G1AffinePoint* in);
blst::BLST_ERROR blst_p1_uncompress(blst::scalar::G1AffinePoint* out, const unsigned char in[48]);
blst::BLST_ERROR blst_p1_deserialize(blst::scalar::G1AffinePoint* out, const unsigned char in[96]);

void blst_p2_add(blst::scalar::G2Point* out, const blst::scalar::G2Point* a, const blst::scalar::G2Point* b);
void blst_p2_sub(blst::scalar::G2Point* out, const blst::scalar::G2Point* a, const blst::scalar::G2Point* b);
void blst_p2_add_or_double(blst::scalar::G2Point* out, const blst::scalar::G2Point* a, const blst::scalar::G2Point* b);
void blst_p2_sub_affine(blst::scalar::G2Point* out, const blst::scalar::G2Point* a, const blst::scalar::G2AffinePoint* b);
void blst_p2_double(blst::scalar::G2Point* out, const blst::scalar::G2Point* a);
void blst_p2_mult(blst::scalar::G2Point* out, const blst::scalar::G2Point* p, const unsigned char* scalar, size_t nbits);
void blst_p2_unchecked_mult(blst::scalar::G2Point* out, const blst::scalar::G2Point* p, const unsigned char* scalar, size_t nbits);
void blst_p2_affine_cneg(blst::scalar::G2AffinePoint* p, bool cbit);
void blst_p2_to_affine(blst::scalar::G2AffinePoint* out, const blst::scalar::G2Point* in);
void blst_p2_from_affine(blst::scalar::G2Point* out, const blst::scalar::G2AffinePoint* in);
void blst_p2_from_jacobian(blst::scalar::G2Point* out, const blst::scalar::G2Point* in);
bool blst_p2_on_curve(const blst::scalar::G2Point* p);
bool blst_p2_in_g2(const blst::scalar::G2Point* p);
bool blst_p2_is_equal(const blst::scalar::G2Point* a, const blst::scalar::G2Point* b);
bool blst_p2_is_inf(const blst::scalar::G2Point* a);
const blst::scalar::G2Point* blst_p2_generator(void);
bool blst_p2_affine_in_g2(const blst::scalar::G2AffinePoint* p);
bool blst_p2_affine_is_equal(const blst::scalar::G2AffinePoint* a, const blst::scalar::G2AffinePoint* b);
const blst::scalar::G2AffinePoint* blst_p2_affine_generator(void);
void blst_p2s_to_affine(blst::scalar::G2AffinePoint dst[], const blst::scalar::G2Point* const points[], size_t npoints);
void blst_p2s_add(blst::scalar::G2Point* ret, const blst::scalar::G2AffinePoint* const points[], size_t npoints);
size_t blst_p2s_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints);
void blst_p2s_mult_wbits_precompute(blst::scalar::G2AffinePoint table[], size_t wbits, const blst::scalar::G2AffinePoint* const points[], size_t npoints);
size_t blst_p2s_mult_wbits_scratch_sizeof(size_t npoints);
void blst_p2s_mult_wbits(blst::scalar::G2Point* ret, const blst::scalar::G2AffinePoint table[], size_t wbits, size_t npoints, const unsigned char* const scalars[], size_t nbits, blst::limb_t* scratch);
size_t blst_p2s_mult_pippenger_scratch_sizeof(size_t npoints);
void blst_p2s_mult_pippenger(blst::scalar::G2Point* ret, const blst::scalar::G2AffinePoint* const points[], size_t npoints, const unsigned char* const scalars[], size_t nbits, blst::limb_t* scratch);
void blst_p2s_tile_pippenger(blst::scalar::G2Point* ret, const blst::scalar::G2AffinePoint* const points[], size_t npoints, const unsigned char* const scalars[], size_t nbits, blst::limb_t* scratch, size_t bit0, size_t window);
void blst_p2_affine_serialize(unsigned char out[192], const blst::scalar::G2AffinePoint* in);
void blst_p2_affine_compress(unsigned char out[96], const blst::scalar::G2AffinePoint* in);
blst::BLST_ERROR blst_p2_uncompress(blst::scalar::G2AffinePoint* out, const unsigned char in[96]);
blst::BLST_ERROR blst_p2_deserialize(blst::scalar::G2AffinePoint* out, const unsigned char in[192]);

void blst_miller_loop(blst::Fp12* ret, const blst::scalar::G2AffinePoint* Q, const blst::scalar::G1AffinePoint* P);
void blst_miller_loop_n(blst::Fp12* ret, const blst::scalar::G2AffinePoint* const Qs[], const blst::scalar::G1AffinePoint* const Ps[], size_t n);
void blst_precompute_lines(blst::Fp6 Qlines[68], const blst::scalar::G2AffinePoint* Q);
void blst_miller_loop_lines(blst::Fp12* ret, const blst::Fp6 Qlines[68], const blst::scalar::G1AffinePoint* P);
void blst_hash_to_g2(blst::scalar::G2Point* out, const unsigned char* msg, size_t msg_len, const unsigned char* DST, size_t DST_len, const unsigned char* aug, size_t aug_len);
void blst_map_to_g1(blst::scalar::G1Point* out, const blst::Fp* u, const blst::Fp* v);
void blst_map_to_g2(blst::scalar::G2Point* out, const blst::Fp2* u, const blst::Fp2* v);
size_t blst_pairing_sizeof(void);
void blst_pairing_init(blst_pairing* new_ctx, bool hash_or_encode, const blst_byte* DST, size_t DST_len);
blst::BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g1(blst_pairing* ctx, const blst::scalar::G1AffinePoint* PK, bool pk_grpchk, const blst::scalar::G2AffinePoint* signature, bool sig_grpchk, const blst_byte* msg, size_t msg_len, const blst_byte* aug, size_t aug_len);
void blst_pairing_commit(blst_pairing* ctx);
bool blst_pairing_finalverify(const blst_pairing* ctx, const blst::Fp12* gtsig);
void blst_aggregated_in_g2(blst::Fp12* out, const blst::scalar::G2AffinePoint* signature);
size_t blst_uniq_sizeof(size_t capacity);
size_t blst_uniq_alignof(void);
bool blst_uniq_init(blst_uniq* tree, size_t capacity);
BLST_UNIQ_RESULT blst_uniq_insert(blst_uniq* tree, const blst_byte* msg, size_t len);
void blst_keygen(blst_scalar* out_SK, const blst_byte* IKM, size_t IKM_len, const blst_byte* info, size_t info_len);
void blst_keygen_v5(blst_scalar* out_SK, const blst_byte* IKM, size_t IKM_len, const blst_byte* salt, size_t salt_len, const blst_byte* info, size_t info_len);
void blst_derive_master_eip2333(blst_scalar* out_SK, const blst_byte* IKM, size_t IKM_len);
void blst_derive_child_eip2333(blst_scalar* out_SK, const blst_scalar* parent_SK, uint32_t child_index);
void blst_bendian_from_scalar(blst_byte out[32], const blst_scalar* scalar);
void blst_expand_message_xmd(unsigned char* out, size_t out_len, const unsigned char* msg, size_t msg_len, const unsigned char* DST, size_t DST_len);
size_t blst_p1_sizeof(void);
size_t blst_p1_affine_sizeof(void);
size_t blst_p2_sizeof(void);
size_t blst_p2_affine_sizeof(void);
size_t blst_fp12_sizeof(void);
bool blst_avx_enabled(void);
size_t blst_avx_p1_sizeof(void);
size_t blst_avx_p1_alignof(void);
size_t blst_avx_p2_sizeof(void);
size_t blst_avx_p2_alignof(void);
size_t blst_avx_fp12_sizeof(void);
size_t blst_avx_fp12_alignof(void);
bool blst_avx_p1_from_scalar(BlstAvxP1* out, const blst::scalar::G1Point* in);
bool blst_p1_from_avx(blst::scalar::G1Point* out, const BlstAvxP1* in);
bool blst_avx_p1_affine_from_scalar(BlstAvxP1Affine* out, const blst::scalar::G1AffinePoint* in);
bool blst_p1_affine_from_avx(blst::scalar::G1AffinePoint* out, const BlstAvxP1Affine* in);
bool blst_avx_p1_add(BlstAvxP1* out, const BlstAvxP1* a, const BlstAvxP1* b);
bool blst_avx_p1_mult(BlstAvxP1* out, const BlstAvxP1* p, const unsigned char* scalar, size_t nbits);
size_t blst_avx_p1_mult_pippenger_scratch_sizeof(size_t npoints);
size_t blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints);
bool blst_avx_p1_mult_pippenger(
    BlstAvxP1* out,
    const BlstAvxP1Affine* points,
    const unsigned char* scalars,
    size_t npoints,
    size_t nbits,
    void* scratch,
    size_t scratch_bytes,
    unsigned char* bucket_used,
    size_t bucket_used_bytes);
bool blst_avx_p1_to_affine(BlstAvxP1Affine* out, const BlstAvxP1* in);
bool blst_avx_p1_from_affine(BlstAvxP1* out, const BlstAvxP1Affine* in);
bool blst_avx_p1_in_g1(const BlstAvxP1* p);
bool blst_avx_p1_is_equal(const BlstAvxP1* a, const BlstAvxP1* b);
const BlstAvxP1* blst_avx_p1_generator(void);
bool blst_avx_p2_from_scalar(BlstAvxP2* out, const blst::scalar::G2Point* in);
bool blst_p2_from_avx(blst::scalar::G2Point* out, const BlstAvxP2* in);
bool blst_avx_p2_affine_from_scalar(BlstAvxP2Affine* out, const blst::scalar::G2AffinePoint* in);
bool blst_p2_affine_from_avx(blst::scalar::G2AffinePoint* out, const BlstAvxP2Affine* in);
bool blst_avx_p2_add(BlstAvxP2* out, const BlstAvxP2* a, const BlstAvxP2* b);
bool blst_avx_p2_mult(BlstAvxP2* out, const BlstAvxP2* p, const unsigned char* scalar, size_t nbits);
size_t blst_avx_p2_mult_pippenger_scratch_sizeof(size_t npoints);
size_t blst_avx_p2_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints);
bool blst_avx_p2_mult_pippenger(
    BlstAvxP2* out,
    const BlstAvxP2Affine* points,
    const unsigned char* scalars,
    size_t npoints,
    size_t nbits,
    void* scratch,
    size_t scratch_bytes,
    unsigned char* bucket_used,
    size_t bucket_used_bytes);
bool blst_avx_p2_to_affine(BlstAvxP2Affine* out, const BlstAvxP2* in);
bool blst_avx_p2_from_affine(BlstAvxP2* out, const BlstAvxP2Affine* in);
bool blst_avx_p2_in_g2(const BlstAvxP2* p);
bool blst_avx_p2_is_equal(const BlstAvxP2* a, const BlstAvxP2* b);
const BlstAvxP2* blst_avx_p2_generator(void);
bool blst_avx_fp12_from_scalar(BlstAvxFp12* out, const blst::Fp12* in);
bool blst_fp12_from_avx(blst::Fp12* out, const BlstAvxFp12* in);
bool blst_avx_fp12_mul(BlstAvxFp12* ret, const BlstAvxFp12* a, const BlstAvxFp12* b);
bool blst_avx_fp12_sqr(BlstAvxFp12* ret, const BlstAvxFp12* a);
bool blst_avx_fp12_final_exp(BlstAvxFp12* ret, const BlstAvxFp12* a);
bool blst_avx_fp12_is_equal(const BlstAvxFp12* a, const BlstAvxFp12* b);
bool blst_avx_fp12_in_group(const BlstAvxFp12* a);
const BlstAvxFp12* blst_avx_fp12_one(void);
}

namespace {

using BoundaryG1 = blst::scalar::G1Point;
using BoundaryG1Affine = blst::scalar::G1AffinePoint;
using BoundaryG2 = blst::scalar::G2Point;
using BoundaryG2Affine = blst::scalar::G2AffinePoint;
using BoundaryGt = blst::Fp12;
using BoundaryFp6 = blst::Fp6;
using BoundaryFp = blst::Fp;
using BoundaryFp2 = blst::Fp2;
using BoundaryFr = blst::Fr;
using PublicG1Affine = blst::G1AffinePoint;
using PublicG2Affine = blst::G2AffinePoint;
using PublicGt = blst::GT;
using PublicMillerResult = blst::MillerResult;

constexpr std::array<unsigned char, 48> kFpModulusBe = {
    0x1a, 0x01, 0x11, 0xea, 0x39, 0x7f, 0xe6, 0x9a,
    0x4b, 0x1b, 0xa7, 0xb6, 0x43, 0x4b, 0xac, 0xd7,
    0x64, 0x77, 0x4b, 0x84, 0xf3, 0x85, 0x12, 0xbf,
    0x67, 0x30, 0xd2, 0xa0, 0xf6, 0xb0, 0xf6, 0x24,
    0x1e, 0xab, 0xff, 0xfe, 0xb1, 0x53, 0xff, 0xff,
    0xb9, 0xfe, 0xff, 0xff, 0xff, 0xff, 0xaa, 0xab,
};

std::array<unsigned char, 32> scalar_bytes(std::uint64_t value)
{
    std::array<unsigned char, 32> out = {};
    for (size_t i = 0; i < sizeof(value); ++i) {
        out[i] = static_cast<unsigned char>(value >> (8 * i));
    }
    return out;
}

bool scalar_is_zero(const blst_scalar& value)
{
    for (blst_byte byte : value.b) {
        if (byte != 0) {
            return false;
        }
    }
    return true;
}

BoundaryGt gt_to_boundary(const PublicGt& value)
{
    unsigned char serialized[12 * 48] = {};
    PublicGt::serialize(serialized, value);
    blst::scalar::GT scalar_gt = {};
    EXPECT_TRUE(blst::scalar::GT::deserialize_unchecked(scalar_gt, serialized));
    return scalar_gt.data;
}

PublicG1Affine to_public(const BoundaryG1Affine& point)
{
#if BLST_DEFAULT_AVX
    return PublicG1Affine(
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.x),
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.y));
#else
    return point;
#endif
}

PublicG2Affine to_public(const BoundaryG2Affine& point)
{
#if BLST_DEFAULT_AVX
    return PublicG2Affine(
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.x),
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.y));
#else
    return point;
#endif
}

BoundaryG1 p1_mul(std::uint64_t scalar)
{
    BoundaryG1 out = {};
    const auto bytes = scalar_bytes(scalar);
    blst_p1_mult(&out, blst_p1_generator(), bytes.data(), 256);
    return out;
}

BoundaryG2 p2_mul(std::uint64_t scalar)
{
    BoundaryG2 out = {};
    const auto bytes = scalar_bytes(scalar);
    blst_p2_mult(&out, blst_p2_generator(), bytes.data(), 256);
    return out;
}

bool fp2_equal(const BoundaryFp2& a, const BoundaryFp2& b)
{
    return BoundaryFp::equal(a[0], b[0]) != 0 && BoundaryFp::equal(a[1], b[1]) != 0;
}

std::vector<unsigned char> packed_scalars(std::initializer_list<unsigned char> scalars)
{
    return std::vector<unsigned char>(scalars);
}

std::vector<unsigned char> scalar_batch(std::initializer_list<std::uint64_t> scalars)
{
    std::vector<unsigned char> out;
    out.reserve(scalars.size() * 32);
    for (std::uint64_t scalar : scalars) {
        const auto bytes = scalar_bytes(scalar);
        out.insert(out.end(), bytes.begin(), bytes.end());
    }
    return out;
}

std::vector<unsigned char> scalar_batch(std::initializer_list<std::uint64_t> scalars, size_t nbytes)
{
    std::vector<unsigned char> out;
    out.reserve(scalars.size() * nbytes);
    for (std::uint64_t scalar : scalars) {
        const auto bytes = scalar_bytes(scalar);
        out.insert(out.end(), bytes.begin(), bytes.begin() + static_cast<std::ptrdiff_t>(nbytes));
    }
    return out;
}

enum class PrecomputePointLayout
{
    packed,
    pointer_array,
    aliased_packed,
};

template<typename Field>
struct ScalarWbitsApi;

template<>
struct ScalarWbitsApi<blst::Fp>
{
    using Point = BoundaryG1;
    using Affine = BoundaryG1Affine;

    static size_t precompute_sizeof(size_t wbits, size_t npoints)
    {
        return blst_p1s_mult_wbits_precompute_sizeof(wbits, npoints);
    }

    static void to_affine(Affine dst[], const Point* const points[], size_t npoints)
    {
        blst_p1s_to_affine(dst, points, npoints);
    }

    static void precompute(
        Affine table[], size_t wbits, const Affine* const points[], size_t npoints)
    {
        blst_p1s_mult_wbits_precompute(table, wbits, points, npoints);
    }

    static size_t scratch_sizeof(size_t npoints)
    {
        return blst_p1s_mult_wbits_scratch_sizeof(npoints);
    }

    static void mult(
        Point* out,
        const Affine table[],
        size_t wbits,
        size_t npoints,
        const unsigned char* const scalars[],
        size_t nbits,
        blst::limb_t* scratch)
    {
        blst_p1s_mult_wbits(
            out, table, wbits, npoints, scalars, nbits, scratch);
    }
};

template<>
struct ScalarWbitsApi<blst::Fp2>
{
    using Point = BoundaryG2;
    using Affine = BoundaryG2Affine;

    static size_t precompute_sizeof(size_t wbits, size_t npoints)
    {
        return blst_p2s_mult_wbits_precompute_sizeof(wbits, npoints);
    }

    static void to_affine(Affine dst[], const Point* const points[], size_t npoints)
    {
        blst_p2s_to_affine(dst, points, npoints);
    }

    static void precompute(
        Affine table[], size_t wbits, const Affine* const points[], size_t npoints)
    {
        blst_p2s_mult_wbits_precompute(table, wbits, points, npoints);
    }

    static size_t scratch_sizeof(size_t npoints)
    {
        return blst_p2s_mult_wbits_scratch_sizeof(npoints);
    }

    static void mult(
        Point* out,
        const Affine table[],
        size_t wbits,
        size_t npoints,
        const unsigned char* const scalars[],
        size_t nbits,
        blst::limb_t* scratch)
    {
        blst_p2s_mult_wbits(
            out, table, wbits, npoints, scalars, nbits, scratch);
    }
};

template<typename Field>
void run_scalar_wbits_precompute_case(
    size_t wbits,
    PrecomputePointLayout layout)
{
    using Api = ScalarWbitsApi<Field>;
    using Point = typename Api::Point;
    using Affine = typename Api::Affine;
    constexpr size_t npoints = 3;

    std::array<Point, npoints> projective = {};
    std::array<Affine, npoints> affine = {};
    const std::array<unsigned char, npoints> point_scalars = {3, 0, 11};
    for (size_t i = 0; i < npoints; ++i) {
        if (point_scalars[i] != 0) {
            Point::mult(
                projective[i],
                Point::generator(),
                &point_scalars[i],
                8);
        }
        Point::to_affine(affine[i], projective[i]);
    }

    const size_t nwin = (size_t)1 << (wbits - 1);
    const size_t total = npoints * nwin;
    ASSERT_EQ(Api::precompute_sizeof(wbits, npoints), total * sizeof(Affine));

    std::vector<Point> reference_rows(total);
    for (size_t i = 0; i < npoints; ++i) {
        blst::scalar::point_precompute_row_wbits<Field>(
            reference_rows.data() + i * nwin, wbits, affine[i]);
    }
    std::vector<Affine> expected(total);
    const Point* reference_ptrs[2] = {reference_rows.data(), nullptr};
    Api::to_affine(expected.data(), reference_ptrs, total);

    std::vector<Affine> actual(total);
    std::array<const Affine*, npoints> pointer_array = {
        &affine[2], &affine[0], &affine[1]};
    const Affine* packed_ptrs[2] = {affine.data(), nullptr};
    if (layout == PrecomputePointLayout::aliased_packed) {
        std::copy(affine.begin(), affine.end(), actual.end() - npoints);
        packed_ptrs[0] = actual.data() + total - npoints;
    }

    const Affine* const* inputs =
        layout == PrecomputePointLayout::pointer_array
        ? pointer_array.data()
        : packed_ptrs;
    if (layout == PrecomputePointLayout::pointer_array) {
        std::array<Affine, npoints> reordered = {
            affine[2], affine[0], affine[1]};
        for (size_t i = 0; i < npoints; ++i) {
            blst::scalar::point_precompute_row_wbits<Field>(
                reference_rows.data() + i * nwin, wbits, reordered[i]);
        }
        Api::to_affine(expected.data(), reference_ptrs, total);
    }

    Api::precompute(actual.data(), wbits, inputs, npoints);

    for (size_t i = 0; i < total; ++i) {
        SCOPED_TRACE(i);
        EXPECT_EQ(std::memcmp(&actual[i], &expected[i], sizeof(Affine)), 0);
    }

    const std::array<unsigned char, npoints> scalar_blob = {2, 5, 7};
    const unsigned char* scalar_ptrs[2] = {scalar_blob.data(), nullptr};
    std::vector<blst::limb_t> scratch(
        Api::scratch_sizeof(npoints) / sizeof(blst::limb_t));
    Point actual_msm = {};
    Api::mult(
        &actual_msm,
        actual.data(),
        wbits,
        npoints,
        scalar_ptrs,
        8,
        scratch.data());

    Point expected_msm = {};
    bool first = true;
    for (size_t i = 0; i < npoints; ++i) {
        const size_t source =
            layout == PrecomputePointLayout::pointer_array
            ? (i == 0 ? 2 : i - 1)
            : i;
        Point term = {};
        Point::mult(term, projective[source], &scalar_blob[i], 8);
        if (first) {
            expected_msm = term;
            first = false;
        } else {
            Point::dadd(expected_msm, expected_msm, term);
        }
    }
    EXPECT_TRUE(Point::equal(actual_msm, expected_msm));
}

} // namespace

TEST(CApi, FrPrimitivesRoundTripAndInverse)
{
    const uint64_t five_words[4] = {5, 0, 0, 0};
    const uint64_t one_words[4] = {1, 0, 0, 0};
    BoundaryFr five = {};
    BoundaryFr one = {};
    BoundaryFr five_from_hex = {};
    blst_scalar scalar = {};
    BoundaryFr from_scalar = {};

    blst_fr_from_uint64(&five, five_words);
    blst_fr_from_uint64(&one, one_words);
    ASSERT_TRUE(blst_fr_from_hexascii(
        &five_from_hex,
        reinterpret_cast<const unsigned char *>(
            "0000000000000000000000000000000000000000000000000000000000000005")));

    EXPECT_TRUE(BoundaryFr::equal(five, five_from_hex));

    blst_scalar_from_fr(&scalar, &five);
    blst_fr_from_scalar(&from_scalar, &scalar);
    EXPECT_TRUE(BoundaryFr::equal(five, from_scalar));

    uint64_t roundtrip[4] = {};
    blst_uint64_from_fr(roundtrip, &five);
    EXPECT_EQ(roundtrip[0], 5ULL);
    EXPECT_EQ(roundtrip[1], 0ULL);

    BoundaryFr inverse_eucl = {};
    BoundaryFr inverse = {};
    BoundaryFr product = {};
    BoundaryFr power = {};
    BoundaryFr power_zero = {};
    blst_fr_eucl_inverse(&inverse_eucl, &five);
    blst_fr_inverse(&inverse, &five);
    blst_fr_mul(&product, &five, &inverse_eucl);
    const std::array<unsigned char, 1> exponent_three = {3};
    blst_fr_pow(&power, &five, exponent_three.data(), 2);
    blst_fr_pow(&power_zero, &five, nullptr, 0);
    EXPECT_TRUE(BoundaryFr::equal(product, one));
    EXPECT_TRUE(BoundaryFr::equal(inverse_eucl, inverse));
    EXPECT_TRUE(BoundaryFr::equal(power, BoundaryFr::from_uint64(125)));
    EXPECT_TRUE(BoundaryFr::equal(power_zero, one));
}

TEST(CApi, ZeroBitScalarShiftsAreIdentity)
{
    const uint64_t fr_words[4] = {0x123456789abcdef0ULL, 0, 0, 0};
    const uint64_t fp_words[6] = {0x0fedcba987654321ULL, 0, 0, 0, 0, 0};
    BoundaryFr fr = {}, shifted_fr = {};
    BoundaryFp fp = {}, shifted_fp = {};
    BoundaryFp2 fp2 = {}, shifted_fp2 = {};

    blst_fr_from_uint64(&fr, fr_words);
    blst_fp_from_uint64(&fp, fp_words);
    fp2[0] = fp;
    fp2[1] = BoundaryFp::one();

    blst_fr_lshift(&shifted_fr, &fr, 0);
    EXPECT_TRUE(BoundaryFr::equal(shifted_fr, fr));
    blst_fr_rshift(&shifted_fr, &shifted_fr, 0);
    EXPECT_TRUE(BoundaryFr::equal(shifted_fr, fr));

    blst_fp_lshift(&shifted_fp, &fp, 0);
    EXPECT_TRUE(BoundaryFp::equal(shifted_fp, fp));
    blst_fp_lshift(&shifted_fp, &shifted_fp, 0);
    EXPECT_TRUE(BoundaryFp::equal(shifted_fp, fp));

    blst_fp2_lshift(&shifted_fp2, &fp2, 0);
    EXPECT_TRUE(blst::vec_is_equal(&shifted_fp2, &fp2, sizeof(fp2)));
    blst_fp2_lshift(&shifted_fp2, &shifted_fp2, 0);
    EXPECT_TRUE(blst::vec_is_equal(&shifted_fp2, &fp2, sizeof(fp2)));
}

TEST(CApi, PublicFieldShiftsBoundExtremeCounts)
{
    const uint64_t fr_words[4] = {0x123456789abcdef0ULL, 0, 0, 0};
    const uint64_t fp_words[6] = {0x0fedcba987654321ULL, 0, 0, 0, 0, 0};
    BoundaryFr fr = {};
    BoundaryFp fp = {};
    blst_fr_from_uint64(&fr, fr_words);
    blst_fp_from_uint64(&fp, fp_words);

    const size_t counts[] = {257, std::numeric_limits<size_t>::max()};
    for (size_t count : counts) {
        std::array<unsigned char, sizeof(size_t)> exponent = {};
        size_t value = count;
        size_t bits = 0;
        for (size_t i = 0; i < exponent.size(); ++i) {
            exponent[i] = static_cast<unsigned char>(value);
            value >>= 8;
        }
        for (value = count; value != 0; value >>= 1)
            ++bits;

        BoundaryFr fr_two;
        BoundaryFr::add(fr_two, BoundaryFr::one(), BoundaryFr::one());
        BoundaryFr fr_factor;
        BoundaryFr::exp(fr_factor, fr_two, exponent.data(), bits);
        BoundaryFr fr_expected;
        BoundaryFr::mul(fr_expected, fr, fr_factor);
        BoundaryFr fr_actual;
        blst_fr_lshift(&fr_actual, &fr, count);
        EXPECT_TRUE(BoundaryFr::equal(fr_actual, fr_expected));

        BoundaryFr fr_half;
        BoundaryFr::rshift(fr_half, BoundaryFr::one(), 1);
        BoundaryFr::exp(fr_factor, fr_half, exponent.data(), bits);
        BoundaryFr::mul(fr_expected, fr, fr_factor);
        blst_fr_rshift(&fr_actual, &fr, count);
        EXPECT_TRUE(BoundaryFr::equal(fr_actual, fr_expected));

        BoundaryFp fp_two;
        BoundaryFp::add(fp_two, BoundaryFp::one(), BoundaryFp::one());
        BoundaryFp fp_factor;
        BoundaryFp::exp(fp_factor, fp_two, exponent.data(), bits);
        BoundaryFp fp_expected;
        BoundaryFp::mul(fp_expected, fp, fp_factor);
        BoundaryFp fp_actual;
        blst_fp_lshift(&fp_actual, &fp, count);
        EXPECT_TRUE(BoundaryFp::equal(fp_actual, fp_expected));

        BoundaryFp2 fp2 = {};
        fp2[0] = fp;
        fp2[1] = BoundaryFp::one();
        BoundaryFp2 fp2_actual;
        blst_fp2_lshift(&fp2_actual, &fp2, count);
        BoundaryFp fp2_expected_1;
        BoundaryFp::mul(fp2_expected_1, fp2[1], fp_factor);
        EXPECT_TRUE(BoundaryFp::equal(fp2_actual[0], fp_expected));
        EXPECT_TRUE(BoundaryFp::equal(fp2_actual[1], fp2_expected_1));
    }
}

TEST(CApi, FpUint32RoundTripAcceptsOnlyUint32Alignment)
{
    alignas(8) std::array<unsigned char, sizeof(unsigned int) * 12 + 4> input_storage = {};
    alignas(8) std::array<unsigned char, sizeof(unsigned int) * 12 + 4> output_storage = {};
    auto* input_words = reinterpret_cast<unsigned int*>(input_storage.data() + 4);
    auto* output_words = reinterpret_cast<unsigned int*>(output_storage.data() + 4);

    ASSERT_EQ(reinterpret_cast<std::uintptr_t>(input_words) % alignof(unsigned int), 0u);
    ASSERT_EQ(reinterpret_cast<std::uintptr_t>(output_words) % alignof(unsigned int), 0u);
    if (alignof(BoundaryFp) > alignof(unsigned int)) {
        EXPECT_NE(reinterpret_cast<std::uintptr_t>(input_words) % alignof(BoundaryFp), 0u);
        EXPECT_NE(reinterpret_cast<std::uintptr_t>(output_words) % alignof(BoundaryFp), 0u);
    }

    for (size_t i = 0; i < 12; ++i) {
        input_words[i] = 0;
        output_words[i] = 0xA5A5A5A5u;
    }
    input_words[0] = 3;

    BoundaryFp from_uint32 = {};
    BoundaryFp from_uint64 = {};
    const uint64_t three_words[6] = {3, 0, 0, 0, 0, 0};
    blst_fp_from_uint32(&from_uint32, input_words);
    blst_fp_from_uint64(&from_uint64, three_words);
    EXPECT_TRUE(BoundaryFp::equal(from_uint32, from_uint64));

    blst_uint32_from_fp(output_words, &from_uint32);
    EXPECT_EQ(output_words[0], 3u);
    for (size_t i = 1; i < 12; ++i) {
        EXPECT_EQ(output_words[i], 0u);
    }
}

TEST(CApi, AvxExpertP1RoundTripAndOps)
{
    if (!blst_avx_enabled()) {
        GTEST_SKIP() << "AVX expert C API is unavailable in this build.";
    }

    EXPECT_EQ(blst_avx_p1_sizeof(), sizeof(BlstAvxP1));
    EXPECT_EQ(blst_avx_p1_alignof(), alignof(BlstAvxP1));

    const BoundaryG1 p = p1_mul(3);
    const BoundaryG1 q = p1_mul(5);
    BlstAvxP1 avx_p = {};
    BlstAvxP1 avx_q = {};
    BlstAvxP1 avx_sum = {};
    BoundaryG1 sum = {};
    BoundaryG1 expected = {};

    ASSERT_TRUE(blst_avx_p1_from_scalar(&avx_p, &p));
    ASSERT_TRUE(blst_avx_p1_from_scalar(&avx_q, &q));
    EXPECT_TRUE(blst_avx_p1_in_g1(&avx_p));
    ASSERT_TRUE(blst_avx_p1_add(&avx_sum, &avx_p, &avx_q));
    ASSERT_TRUE(blst_p1_from_avx(&sum, &avx_sum));
    blst_p1_add(&expected, &p, &q);
    EXPECT_TRUE(BoundaryG1::equal(sum, expected));

    BoundaryG1Affine p_affine = {};
    blst_p1_to_affine(&p_affine, &p);
    BlstAvxP1Affine avx_p_affine = {};
    BoundaryG1Affine roundtrip_affine = {};
    ASSERT_TRUE(blst_avx_p1_affine_from_scalar(&avx_p_affine, &p_affine));
    ASSERT_TRUE(blst_p1_affine_from_avx(&roundtrip_affine, &avx_p_affine));
    EXPECT_TRUE(BoundaryG1Affine::equal(p_affine, roundtrip_affine));

    const auto* avx_generator = blst_avx_p1_generator();
    ASSERT_NE(avx_generator, nullptr);
    const auto scalar = scalar_bytes(7);
    BlstAvxP1 avx_product = {};
    BoundaryG1 product = {};
    ASSERT_TRUE(blst_avx_p1_mult(&avx_product, avx_generator, scalar.data(), 256));
    ASSERT_TRUE(blst_p1_from_avx(&product, &avx_product));
    EXPECT_TRUE(BoundaryG1::equal(product, p1_mul(7)));
    EXPECT_TRUE(blst_avx_p1_is_equal(avx_generator, avx_generator));
}

TEST(CApi, AvxExpertP2RoundTripAndOps)
{
    if (!blst_avx_enabled()) {
        GTEST_SKIP() << "AVX expert C API is unavailable in this build.";
    }

    EXPECT_EQ(blst_avx_p2_sizeof(), sizeof(BlstAvxP2));
    EXPECT_EQ(blst_avx_p2_alignof(), alignof(BlstAvxP2));

    const BoundaryG2 p = p2_mul(11);
    const BoundaryG2 q = p2_mul(13);
    BlstAvxP2 avx_p = {};
    BlstAvxP2 avx_q = {};
    BlstAvxP2 avx_sum = {};
    BoundaryG2 sum = {};
    BoundaryG2 expected = {};

    ASSERT_TRUE(blst_avx_p2_from_scalar(&avx_p, &p));
    ASSERT_TRUE(blst_avx_p2_from_scalar(&avx_q, &q));
    EXPECT_TRUE(blst_avx_p2_in_g2(&avx_p));
    ASSERT_TRUE(blst_avx_p2_add(&avx_sum, &avx_p, &avx_q));
    ASSERT_TRUE(blst_p2_from_avx(&sum, &avx_sum));
    blst_p2_add(&expected, &p, &q);
    EXPECT_TRUE(BoundaryG2::equal(sum, expected));

    BoundaryG2Affine p_affine = {};
    blst_p2_to_affine(&p_affine, &p);
    BlstAvxP2Affine avx_p_affine = {};
    BoundaryG2Affine roundtrip_affine = {};
    ASSERT_TRUE(blst_avx_p2_affine_from_scalar(&avx_p_affine, &p_affine));
    ASSERT_TRUE(blst_p2_affine_from_avx(&roundtrip_affine, &avx_p_affine));
    EXPECT_TRUE(BoundaryG2Affine::equal(p_affine, roundtrip_affine));

    const auto* avx_generator = blst_avx_p2_generator();
    ASSERT_NE(avx_generator, nullptr);
    const auto scalar = scalar_bytes(17);
    BlstAvxP2 avx_product = {};
    BoundaryG2 product = {};
    ASSERT_TRUE(blst_avx_p2_mult(&avx_product, avx_generator, scalar.data(), 256));
    ASSERT_TRUE(blst_p2_from_avx(&product, &avx_product));
    EXPECT_TRUE(BoundaryG2::equal(product, p2_mul(17)));
    EXPECT_TRUE(blst_avx_p2_is_equal(avx_generator, avx_generator));
}

TEST(CApi, AvxExpertP1PippengerMatchesScalar)
{
    if (!blst_avx_enabled()) {
        GTEST_SKIP() << "AVX expert C API is unavailable in this build.";
    }

    constexpr size_t kCount = 4;
    const std::array<std::uint64_t, kCount> point_scalars = {3, 5, 7, 11};
    const std::array<std::uint64_t, kCount> msm_scalars = {13, 17, 19, 23};
    std::vector<unsigned char> scalars =
        scalar_batch({msm_scalars[0], msm_scalars[1], msm_scalars[2], msm_scalars[3]});
    std::fill(scalars.end() - 32, scalars.end(), 0xff);
    std::array<BlstAvxP1Affine, kCount> avx_points = {};
    BoundaryG1 expected = {};

    for (size_t i = 0; i < kCount; ++i) {
        BoundaryG1 point = p1_mul(point_scalars[i]);
        BoundaryG1Affine affine = {};
        blst_p1_to_affine(&affine, &point);
        ASSERT_TRUE(blst_avx_p1_affine_from_scalar(&avx_points[i], &affine));

        BoundaryG1 term = {};
        blst_p1_mult(&term, &point, scalars.data() + i * 32, 256);
        blst_p1_add(&expected, &expected, &term);
    }

    const size_t scratch_bytes = blst_avx_p1_mult_pippenger_scratch_sizeof(kCount);
    const size_t bucket_used_bytes = blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof(kCount);
    ASSERT_GT(scratch_bytes, 0u);
    ASSERT_GT(bucket_used_bytes, 0u);
    std::vector<BlstAvxP1> scratch((scratch_bytes + sizeof(BlstAvxP1) - 1) / sizeof(BlstAvxP1));
    std::vector<unsigned char> bucket_used(bucket_used_bytes);

    BlstAvxP1 msm = {};
    BoundaryG1 roundtrip = {};
    ASSERT_TRUE(blst_avx_p1_mult_pippenger(
        &msm,
        avx_points.data(),
        scalars.data(),
        kCount,
        256,
        scratch.data(),
        scratch.size() * sizeof(BlstAvxP1),
        bucket_used.data(),
        bucket_used.size()));
    ASSERT_TRUE(blst_p1_from_avx(&roundtrip, &msm));
    EXPECT_TRUE(BoundaryG1::equal(roundtrip, expected));
}

TEST(CApi, AvxExpertP2PippengerMatchesScalar)
{
    if (!blst_avx_enabled()) {
        GTEST_SKIP() << "AVX expert C API is unavailable in this build.";
    }

    constexpr size_t kCount = 4;
    const std::array<std::uint64_t, kCount> point_scalars = {3, 5, 7, 11};
    const std::array<std::uint64_t, kCount> msm_scalars = {13, 17, 19, 23};
    std::vector<unsigned char> scalars =
        scalar_batch({msm_scalars[0], msm_scalars[1], msm_scalars[2], msm_scalars[3]});
    std::fill(scalars.end() - 32, scalars.end(), 0xff);
    std::array<BlstAvxP2Affine, kCount> avx_points = {};
    BoundaryG2 expected = {};

    for (size_t i = 0; i < kCount; ++i) {
        BoundaryG2 point = p2_mul(point_scalars[i]);
        BoundaryG2Affine affine = {};
        blst_p2_to_affine(&affine, &point);
        ASSERT_TRUE(blst_avx_p2_affine_from_scalar(&avx_points[i], &affine));

        BoundaryG2 term = {};
        blst_p2_mult(&term, &point, scalars.data() + i * 32, 256);
        blst_p2_add(&expected, &expected, &term);
    }

    const size_t scratch_bytes = blst_avx_p2_mult_pippenger_scratch_sizeof(kCount);
    const size_t bucket_used_bytes = blst_avx_p2_mult_pippenger_bucket_used_scratch_sizeof(kCount);
    ASSERT_GT(scratch_bytes, 0u);
    ASSERT_GT(bucket_used_bytes, 0u);
    std::vector<BlstAvxP2> scratch((scratch_bytes + sizeof(BlstAvxP2) - 1) / sizeof(BlstAvxP2));
    std::vector<unsigned char> bucket_used(bucket_used_bytes);

    BlstAvxP2 msm = {};
    BoundaryG2 roundtrip = {};
    ASSERT_TRUE(blst_avx_p2_mult_pippenger(
        &msm,
        avx_points.data(),
        scalars.data(),
        kCount,
        256,
        scratch.data(),
        scratch.size() * sizeof(BlstAvxP2),
        bucket_used.data(),
        bucket_used.size()));
    ASSERT_TRUE(blst_p2_from_avx(&roundtrip, &msm));
    EXPECT_TRUE(BoundaryG2::equal(roundtrip, expected));
}

TEST(CApi, AvxExpertP1PippengerShortScalarsMatchScalar)
{
    if (!blst_avx_enabled()) {
        GTEST_SKIP() << "AVX expert C API is unavailable in this build.";
    }

    constexpr size_t kCount = 4;
    constexpr size_t kBits = 128;
    constexpr size_t kBytes = kBits / 8;
    const std::array<std::uint64_t, kCount> point_scalars = {3, 5, 7, 11};
    const std::array<std::uint64_t, kCount> msm_scalars = {13, 17, 19, 23};
    const std::vector<unsigned char> scalars =
        scalar_batch({msm_scalars[0], msm_scalars[1], msm_scalars[2], msm_scalars[3]}, kBytes);
    std::array<BlstAvxP1Affine, kCount> avx_points = {};
    BoundaryG1 expected = {};

    for (size_t i = 0; i < kCount; ++i) {
        BoundaryG1 point = p1_mul(point_scalars[i]);
        BoundaryG1Affine affine = {};
        blst_p1_to_affine(&affine, &point);
        ASSERT_TRUE(blst_avx_p1_affine_from_scalar(&avx_points[i], &affine));

        BoundaryG1 term = {};
        blst_p1_mult(&term, &point, scalars.data() + i * kBytes, kBits);
        blst_p1_add(&expected, &expected, &term);
    }

    const size_t scratch_bytes = blst_avx_p1_mult_pippenger_scratch_sizeof(kCount);
    const size_t bucket_used_bytes = blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof(kCount);
    ASSERT_GT(scratch_bytes, 0u);
    ASSERT_GT(bucket_used_bytes, 0u);
    std::vector<BlstAvxP1> scratch((scratch_bytes + sizeof(BlstAvxP1) - 1) / sizeof(BlstAvxP1));
    std::vector<unsigned char> bucket_used(bucket_used_bytes);

    BlstAvxP1 msm = {};
    BoundaryG1 roundtrip = {};
    ASSERT_TRUE(blst_avx_p1_mult_pippenger(
        &msm,
        avx_points.data(),
        scalars.data(),
        kCount,
        kBits,
        scratch.data(),
        scratch.size() * sizeof(BlstAvxP1),
        bucket_used.data(),
        bucket_used.size()));
    ASSERT_TRUE(blst_p1_from_avx(&roundtrip, &msm));
    EXPECT_TRUE(BoundaryG1::equal(roundtrip, expected));
}

TEST(CApi, AvxExpertP2PippengerShortScalarsMatchScalar)
{
    if (!blst_avx_enabled()) {
        GTEST_SKIP() << "AVX expert C API is unavailable in this build.";
    }

    constexpr size_t kCount = 4;
    constexpr size_t kBits = 128;
    constexpr size_t kBytes = kBits / 8;
    const std::array<std::uint64_t, kCount> point_scalars = {3, 5, 7, 11};
    const std::array<std::uint64_t, kCount> msm_scalars = {13, 17, 19, 23};
    const std::vector<unsigned char> scalars =
        scalar_batch({msm_scalars[0], msm_scalars[1], msm_scalars[2], msm_scalars[3]}, kBytes);
    std::array<BlstAvxP2Affine, kCount> avx_points = {};
    BoundaryG2 expected = {};

    for (size_t i = 0; i < kCount; ++i) {
        BoundaryG2 point = p2_mul(point_scalars[i]);
        BoundaryG2Affine affine = {};
        blst_p2_to_affine(&affine, &point);
        ASSERT_TRUE(blst_avx_p2_affine_from_scalar(&avx_points[i], &affine));

        BoundaryG2 term = {};
        blst_p2_mult(&term, &point, scalars.data() + i * kBytes, kBits);
        blst_p2_add(&expected, &expected, &term);
    }

    const size_t scratch_bytes = blst_avx_p2_mult_pippenger_scratch_sizeof(kCount);
    const size_t bucket_used_bytes = blst_avx_p2_mult_pippenger_bucket_used_scratch_sizeof(kCount);
    ASSERT_GT(scratch_bytes, 0u);
    ASSERT_GT(bucket_used_bytes, 0u);
    std::vector<BlstAvxP2> scratch((scratch_bytes + sizeof(BlstAvxP2) - 1) / sizeof(BlstAvxP2));
    std::vector<unsigned char> bucket_used(bucket_used_bytes);

    BlstAvxP2 msm = {};
    BoundaryG2 roundtrip = {};
    ASSERT_TRUE(blst_avx_p2_mult_pippenger(
        &msm,
        avx_points.data(),
        scalars.data(),
        kCount,
        kBits,
        scratch.data(),
        scratch.size() * sizeof(BlstAvxP2),
        bucket_used.data(),
        bucket_used.size()));
    ASSERT_TRUE(blst_p2_from_avx(&roundtrip, &msm));
    EXPECT_TRUE(BoundaryG2::equal(roundtrip, expected));
}

TEST(CApi, AvxExpertPippengerRejectsOverflowSizedInputs)
{
    if (!blst_avx_enabled()) {
        GTEST_SKIP() << "AVX expert C API is unavailable in this build.";
    }

    constexpr size_t kHuge = std::numeric_limits<size_t>::max();
    EXPECT_EQ(blst_avx_p1_mult_pippenger_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst_avx_p2_mult_pippenger_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst_avx_p2_mult_pippenger_bucket_used_scratch_sizeof(kHuge), 0u);

    BlstAvxP1 p1_out = {};
    BlstAvxP1Affine p1_point = {};
    BlstAvxP1 p1_scratch = {};
    BlstAvxP2 p2_out = {};
    BlstAvxP2Affine p2_point = {};
    BlstAvxP2 p2_scratch = {};
    std::array<unsigned char, 32> scalar = {};
    std::array<unsigned char, 1> bucket_used = {};

    EXPECT_FALSE(blst_avx_p1_mult_pippenger(
        &p1_out,
        &p1_point,
        scalar.data(),
        kHuge,
        256,
        &p1_scratch,
        sizeof(p1_scratch),
        bucket_used.data(),
        bucket_used.size()));
    EXPECT_FALSE(blst_avx_p2_mult_pippenger(
        &p2_out,
        &p2_point,
        scalar.data(),
        kHuge,
        256,
        &p2_scratch,
        sizeof(p2_scratch),
        bucket_used.data(),
        bucket_used.size()));

    const BlstAvxP1* p1_generator = blst_avx_p1_generator();
    const BlstAvxP2* p2_generator = blst_avx_p2_generator();
    ASSERT_NE(p1_generator, nullptr);
    ASSERT_NE(p2_generator, nullptr);
    EXPECT_FALSE(blst_avx_p1_mult(&p1_out, p1_generator, scalar.data(), 257));
    EXPECT_FALSE(blst_avx_p2_mult(&p2_out, p2_generator, scalar.data(), 257));
    EXPECT_FALSE(blst_avx_p1s_mult_pippenger(
        &p1_out, p1_generator, 1, scalar.data(), 257, nullptr));
    EXPECT_FALSE(blst_avx_p2s_mult_pippenger(
        &p2_out, p2_generator, 1, scalar.data(), 257, nullptr));
}

TEST(CApi, AvxExpertFp12RoundTripAndOps)
{
    if (!blst_avx_enabled()) {
        GTEST_SKIP() << "AVX expert C API is unavailable in this build.";
    }

    EXPECT_EQ(blst_avx_fp12_sizeof(), sizeof(BlstAvxFp12));
    EXPECT_EQ(blst_avx_fp12_alignof(), alignof(BlstAvxFp12));

    const BoundaryG1 p = p1_mul(3);
    const BoundaryG2 q = p2_mul(11);
    BoundaryG1Affine p_affine = {};
    BoundaryG2Affine q_affine = {};
    blst_p1_to_affine(&p_affine, &p);
    blst_p2_to_affine(&q_affine, &q);

    BoundaryGt miller = {};
    BoundaryGt expected = {};
    blst_miller_loop(&miller, &q_affine, &p_affine);
    blst_final_exp(&expected, &miller);

    BlstAvxFp12 avx_miller = {};
    BlstAvxFp12 avx_expected = {};
    BlstAvxFp12 avx_product = {};
    BlstAvxFp12 avx_square = {};
    BoundaryGt roundtrip = {};
    BoundaryGt scalar_product = {};
    BoundaryGt scalar_square = {};
    ASSERT_TRUE(blst_avx_fp12_from_scalar(&avx_miller, &miller));
    ASSERT_FALSE(blst_fp12_in_group(&miller));

    blst_fp12_sqr(&scalar_square, &miller);
    ASSERT_TRUE(blst_avx_fp12_sqr(&avx_square, &avx_miller));
    ASSERT_TRUE(blst_fp12_from_avx(&roundtrip, &avx_square));
    EXPECT_TRUE(BoundaryGt::equal(roundtrip, scalar_square));

    avx_product = avx_miller;
    ASSERT_TRUE(blst_avx_fp12_mul(&avx_product, &avx_product, &avx_miller));
    ASSERT_TRUE(blst_fp12_from_avx(&roundtrip, &avx_product));
    EXPECT_TRUE(BoundaryGt::equal(roundtrip, scalar_square));

    avx_square = avx_miller;
    ASSERT_TRUE(blst_avx_fp12_sqr(&avx_square, &avx_square));
    ASSERT_TRUE(blst_fp12_from_avx(&roundtrip, &avx_square));
    EXPECT_TRUE(BoundaryGt::equal(roundtrip, scalar_square));

    ASSERT_TRUE(blst_avx_fp12_final_exp(&avx_expected, &avx_miller));
    ASSERT_TRUE(blst_fp12_from_avx(&roundtrip, &avx_expected));
    EXPECT_TRUE(BoundaryGt::equal(roundtrip, expected));
    EXPECT_TRUE(blst_avx_fp12_in_group(&avx_expected));

    avx_miller = {};
    ASSERT_TRUE(blst_avx_fp12_from_scalar(&avx_miller, &miller));
    ASSERT_TRUE(blst_avx_fp12_final_exp(&avx_miller, &avx_miller));
    ASSERT_TRUE(blst_fp12_from_avx(&roundtrip, &avx_miller));
    EXPECT_TRUE(BoundaryGt::equal(roundtrip, expected));

    const auto* avx_one = blst_avx_fp12_one();
    ASSERT_NE(avx_one, nullptr);
    ASSERT_TRUE(blst_avx_fp12_mul(&avx_product, &avx_expected, avx_one));
    ASSERT_TRUE(blst_fp12_from_avx(&scalar_product, &avx_product));
    EXPECT_TRUE(BoundaryGt::equal(scalar_product, expected));
    EXPECT_TRUE(blst_avx_fp12_is_equal(&avx_product, &avx_expected));
}

TEST(CApi, FpPrimitivesRoundTripInverseAndSqrt)
{
    const uint64_t three_words[6] = {3, 0, 0, 0, 0, 0};
    const uint64_t one_words[6] = {1, 0, 0, 0, 0, 0};
    BoundaryFp three = {};
    BoundaryFp one = {};
    BoundaryFp from_hex = {};
    blst_fp_from_uint64(&three, three_words);
    blst_fp_from_uint64(&one, one_words);
    ASSERT_TRUE(blst_fp_from_hexascii(
        &from_hex,
        reinterpret_cast<const unsigned char *>(
            "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000003")));
    EXPECT_TRUE(BoundaryFp::equal(three, from_hex));

    unsigned char bendian[48] = {};
    BoundaryFp roundtrip = {};
    blst_bendian_from_fp(bendian, &three);
    blst_fp_from_bendian(&roundtrip, bendian);
    EXPECT_TRUE(BoundaryFp::equal(three, roundtrip));

    uint64_t roundtrip_words[6] = {};
    blst_uint64_from_fp(roundtrip_words, &three);
    EXPECT_EQ(roundtrip_words[0], 3ULL);
    EXPECT_EQ(roundtrip_words[1], 0ULL);

    BoundaryFp inverse_eucl = {};
    BoundaryFp inverse = {};
    BoundaryFp product = {};
    blst_fp_eucl_inverse(&inverse_eucl, &three);
    blst_fp_inverse(&inverse, &three);
    blst_fp_mul(&product, &three, &inverse_eucl);
    EXPECT_TRUE(BoundaryFp::equal(product, one));
    blst_fp_mul(&product, &three, &inverse);
    EXPECT_TRUE(BoundaryFp::equal(product, one));

    BoundaryFp square = {};
    BoundaryFp sqrt = {};
    blst_fp_sqr(&square, &three);
    EXPECT_TRUE(blst_fp_is_square(&square));
    ASSERT_TRUE(blst_fp_sqrt(&sqrt, &square));
    blst_fp_sqr(&product, &sqrt);
    EXPECT_TRUE(BoundaryFp::equal(product, square));

    BoundaryFp nonresidue = {};
    bool found_nonresidue = false;
    for (uint64_t value = 2; value < 32; ++value) {
        const uint64_t words[6] = {value, 0, 0, 0, 0, 0};
        blst_fp_from_uint64(&nonresidue, words);
        if (!blst_fp_is_square(&nonresidue)) {
            found_nonresidue = true;
            break;
        }
    }
    ASSERT_TRUE(found_nonresidue);
    ASSERT_FALSE(blst_fp_sqrt(&nonresidue, &nonresidue));
    EXPECT_TRUE(BoundaryFp::equal(nonresidue, BoundaryFp{}));
}

TEST(CApi, Fp2PrimitivesInverseAndSqrt)
{
    const uint64_t two_words[6] = {2, 0, 0, 0, 0, 0};
    const uint64_t one_words[6] = {1, 0, 0, 0, 0, 0};
    BoundaryFp two = {};
    BoundaryFp one = {};
    blst_fp_from_uint64(&two, two_words);
    blst_fp_from_uint64(&one, one_words);

    BoundaryFp2 value = {};
    value[0] = two;
    value[1] = one;

    BoundaryFp2 inverse_eucl = {};
    BoundaryFp2 inverse = {};
    BoundaryFp2 product = {};
    blst_fp2_eucl_inverse(&inverse_eucl, &value);
    blst_fp2_inverse(&inverse, &value);
    blst_fp2_mul(&product, &value, &inverse_eucl);
    EXPECT_TRUE(BoundaryFp::equal(product[0], one));
    EXPECT_TRUE(BoundaryFp::equal(product[1], BoundaryFp{}));
    blst_fp2_mul(&product, &value, &inverse);
    EXPECT_TRUE(BoundaryFp::equal(product[0], one));
    EXPECT_TRUE(BoundaryFp::equal(product[1], BoundaryFp{}));

    BoundaryFp2 square = {};
    BoundaryFp2 sqrt = {};
    blst_fp2_sqr(&square, &value);
    EXPECT_TRUE(blst_fp2_is_square(&square));
    ASSERT_TRUE(blst_fp2_sqrt(&sqrt, &square));
    blst_fp2_sqr(&product, &sqrt);
    EXPECT_TRUE(fp2_equal(product, square));
}

TEST(CApi, G1ProjectiveAndAffineRoundTrip)
{
    const auto scalar = scalar_bytes(5);
    BoundaryG1 point = {};
    blst_p1_mult(&point, blst_p1_generator(), scalar.data(), 256);

    ASSERT_TRUE(blst_p1_on_curve(&point));
    ASSERT_TRUE(blst_p1_in_g1(&point));
    ASSERT_FALSE(blst_p1_is_inf(&point));

    BoundaryG1 expected = {};
    BoundaryG1::mult(expected, BoundaryG1::generator(), scalar.data(), 256);
    EXPECT_TRUE(blst_p1_is_equal(&point, &expected));

    BoundaryG1 doubled = {};
    BoundaryG1 added = {};
    blst_p1_double(&doubled, &point);
    blst_p1_add_or_double(&added, &point, &point);
    EXPECT_TRUE(blst_p1_is_equal(&doubled, &added));

    BoundaryG1Affine affine = {};
    blst_p1_to_affine(&affine, &point);
    ASSERT_TRUE(blst_p1_affine_in_g1(&affine));

    unsigned char compressed[48] = {};
    unsigned char serialized[96] = {};
    blst_p1_affine_compress(compressed, &affine);
    blst_p1_affine_serialize(serialized, &affine);

    BoundaryG1Affine from_compressed = {};
    BoundaryG1Affine from_serialized = {};
    ASSERT_EQ(blst_p1_uncompress(&from_compressed, compressed), blst::BLST_SUCCESS);
    ASSERT_EQ(blst_p1_deserialize(&from_serialized, serialized), blst::BLST_SUCCESS);
    EXPECT_TRUE(blst_p1_affine_is_equal(&affine, &from_compressed));
    EXPECT_TRUE(blst_p1_affine_is_equal(&affine, &from_serialized));

    BoundaryG1 restored = {};
    blst_p1_from_affine(&restored, &from_compressed);
    EXPECT_TRUE(blst_p1_is_equal(&point, &restored));
}

TEST(CApi, G2ProjectiveAndAffineRoundTrip)
{
    const auto scalar = scalar_bytes(7);
    BoundaryG2 point = {};
    blst_p2_mult(&point, blst_p2_generator(), scalar.data(), 256);

    ASSERT_TRUE(blst_p2_on_curve(&point));
    ASSERT_TRUE(blst_p2_in_g2(&point));
    ASSERT_FALSE(blst_p2_is_inf(&point));

    BoundaryG2 expected = {};
    BoundaryG2::mult(expected, BoundaryG2::generator(), scalar.data(), 256);
    EXPECT_TRUE(blst_p2_is_equal(&point, &expected));

    BoundaryG2 doubled = {};
    BoundaryG2 added = {};
    blst_p2_double(&doubled, &point);
    blst_p2_add_or_double(&added, &point, &point);
    EXPECT_TRUE(blst_p2_is_equal(&doubled, &added));

    BoundaryG2Affine affine = {};
    blst_p2_to_affine(&affine, &point);
    ASSERT_TRUE(blst_p2_affine_in_g2(&affine));

    unsigned char compressed[96] = {};
    unsigned char serialized[192] = {};
    blst_p2_affine_compress(compressed, &affine);
    blst_p2_affine_serialize(serialized, &affine);

    BoundaryG2Affine from_compressed = {};
    BoundaryG2Affine from_serialized = {};
    ASSERT_EQ(blst_p2_uncompress(&from_compressed, compressed), blst::BLST_SUCCESS);
    ASSERT_EQ(blst_p2_deserialize(&from_serialized, serialized), blst::BLST_SUCCESS);
    EXPECT_TRUE(blst_p2_affine_is_equal(&affine, &from_compressed));
    EXPECT_TRUE(blst_p2_affine_is_equal(&affine, &from_serialized));

    BoundaryG2 restored = {};
    blst_p2_from_affine(&restored, &from_compressed);
    EXPECT_TRUE(blst_p2_is_equal(&point, &restored));
}

TEST(CApi, PointDecodersRejectCoordinatesAtTheFieldModulus)
{
    BoundaryG1Affine g1 = {};
    std::array<unsigned char, 48> g1_compressed = kFpModulusBe;
    g1_compressed[0] |= 0x80;
    EXPECT_EQ(
        blst_p1_uncompress(&g1, g1_compressed.data()),
        blst::BLST_BAD_ENCODING);

    for (size_t offset : {size_t{0}, size_t{48}}) {
        std::array<unsigned char, 96> encoded = {};
        std::copy(kFpModulusBe.begin(), kFpModulusBe.end(), encoded.begin() + offset);
        EXPECT_EQ(
            blst_p1_deserialize(&g1, encoded.data()),
            blst::BLST_BAD_ENCODING)
            << "coordinate offset " << offset;
    }

    BoundaryG2Affine g2 = {};
    for (size_t offset : {size_t{0}, size_t{48}}) {
        std::array<unsigned char, 96> encoded = {};
        std::copy(kFpModulusBe.begin(), kFpModulusBe.end(), encoded.begin() + offset);
        encoded[0] |= 0x80;
        EXPECT_EQ(
            blst_p2_uncompress(&g2, encoded.data()),
            blst::BLST_BAD_ENCODING)
            << "coordinate offset " << offset;
    }

    for (size_t offset : {size_t{0}, size_t{48}, size_t{96}, size_t{144}}) {
        std::array<unsigned char, 192> encoded = {};
        std::copy(kFpModulusBe.begin(), kFpModulusBe.end(), encoded.begin() + offset);
        EXPECT_EQ(
            blst_p2_deserialize(&g2, encoded.data()),
            blst::BLST_BAD_ENCODING)
            << "coordinate offset " << offset;
    }
}

TEST(CApi, G1MultipointAffineAndAddMatchSinglePointOps)
{
    std::array<BoundaryG1, 3> points = {p1_mul(3), p1_mul(5), p1_mul(7)};
    std::array<BoundaryG1Affine, 3> affines = {};
    const BoundaryG1* point_ptrs[2] = {points.data(), nullptr};
    const BoundaryG1Affine* affine_ptrs[2] = {affines.data(), nullptr};

    blst_p1s_to_affine(affines.data(), point_ptrs, points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        BoundaryG1Affine expected = {};
        blst_p1_to_affine(&expected, &points[i]);
        EXPECT_TRUE(blst_p1_affine_is_equal(&affines[i], &expected));
    }

    BoundaryG1 sum = {};
    blst_p1s_add(&sum, affine_ptrs, affines.size());

    BoundaryG1 expected_sum = points[0];
    blst_p1_add_or_double(&expected_sum, &expected_sum, &points[1]);
    blst_p1_add_or_double(&expected_sum, &expected_sum, &points[2]);
    EXPECT_TRUE(blst_p1_is_equal(&sum, &expected_sum));
}

TEST(CApi, G1MultipointPippengerMatchesSinglePointOps)
{
    std::array<BoundaryG1, 3> points = {p1_mul(3), p1_mul(5), p1_mul(7)};
    std::array<BoundaryG1Affine, 3> affines = {};
    const BoundaryG1* point_ptrs[2] = {points.data(), nullptr};
    blst_p1s_to_affine(affines.data(), point_ptrs, points.size());
    const BoundaryG1Affine* affine_ptrs[2] = {affines.data(), nullptr};

    const std::vector<unsigned char> scalar_blob = packed_scalars({2, 3, 5});
    const unsigned char* scalar_ptrs[2] = {scalar_blob.data(), nullptr};

    std::vector<blst::limb_t> pippenger_scratch(
        blst_p1s_mult_pippenger_scratch_sizeof(affines.size()) / sizeof(blst::limb_t));
    BoundaryG1 pippenger = {};
    blst_p1s_mult_pippenger(&pippenger, affine_ptrs, affines.size(), scalar_ptrs, 8, pippenger_scratch.data());

    BoundaryG1 expected_msm = {};
    BoundaryG1 term = {};
    const auto s0 = scalar_bytes(2);
    blst_p1_mult(&expected_msm, &points[0], s0.data(), 8);
    const auto s1 = scalar_bytes(3);
    blst_p1_mult(&term, &points[1], s1.data(), 8);
    blst_p1_add_or_double(&expected_msm, &expected_msm, &term);
    const auto s2 = scalar_bytes(5);
    blst_p1_mult(&term, &points[2], s2.data(), 8);
    blst_p1_add_or_double(&expected_msm, &expected_msm, &term);
    EXPECT_TRUE(blst_p1_is_equal(&pippenger, &expected_msm));
}

TEST(CApi, G1PippengerAcceptsExactDescriptorArrays)
{
    std::array<BoundaryG1, 3> points = {p1_mul(3), p1_mul(5), p1_mul(7)};
    std::array<BoundaryG1Affine, 3> affines = {};
    for (size_t i = 0; i < points.size(); ++i) {
        blst_p1_to_affine(&affines[i], &points[i]);
    }
    const BoundaryG1Affine* point_ptrs[3] = {
        &affines[0], &affines[1], &affines[2]
    };
    const auto s0 = scalar_bytes(2);
    const auto s1 = scalar_bytes(3);
    const auto s2 = scalar_bytes(5);
    const unsigned char* scalar_ptrs[3] = {s0.data(), s1.data(), s2.data()};
    std::vector<blst::limb_t> scratch(
        blst_p1s_mult_pippenger_scratch_sizeof(points.size()) / sizeof(blst::limb_t));

    BoundaryG1 out = {};
    blst_p1s_mult_pippenger(
        &out, point_ptrs, points.size(), scalar_ptrs, 8, scratch.data());

    BoundaryG1 expected = {};
    BoundaryG1 term = {};
    blst_p1_mult(&expected, &points[0], s0.data(), 8);
    blst_p1_mult(&term, &points[1], s1.data(), 8);
    blst_p1_add_or_double(&expected, &expected, &term);
    blst_p1_mult(&term, &points[2], s2.data(), 8);
    blst_p1_add_or_double(&expected, &expected, &term);
    EXPECT_TRUE(blst_p1_is_equal(&out, &expected));
}

TEST(CApi, G2PippengerAcceptsExactDescriptorArrays)
{
    std::array<BoundaryG2, 3> points = {p2_mul(3), p2_mul(5), p2_mul(7)};
    std::array<BoundaryG2Affine, 3> affines = {};
    for (size_t i = 0; i < points.size(); ++i) {
        blst_p2_to_affine(&affines[i], &points[i]);
    }
    const BoundaryG2Affine* point_ptrs[3] = {
        &affines[0], &affines[1], &affines[2]
    };
    const auto s0 = scalar_bytes(2);
    const auto s1 = scalar_bytes(3);
    const auto s2 = scalar_bytes(5);
    const unsigned char* scalar_ptrs[3] = {s0.data(), s1.data(), s2.data()};
    std::vector<blst::limb_t> scratch(
        blst_p2s_mult_pippenger_scratch_sizeof(points.size()) / sizeof(blst::limb_t));

    BoundaryG2 out = {};
    blst_p2s_mult_pippenger(
        &out, point_ptrs, points.size(), scalar_ptrs, 8, scratch.data());

    BoundaryG2 expected = {};
    BoundaryG2 term = {};
    blst_p2_mult(&expected, &points[0], s0.data(), 8);
    blst_p2_mult(&term, &points[1], s1.data(), 8);
    blst_p2_add_or_double(&expected, &expected, &term);
    blst_p2_mult(&term, &points[2], s2.data(), 8);
    blst_p2_add_or_double(&expected, &expected, &term);
    EXPECT_TRUE(blst_p2_is_equal(&out, &expected));
}

TEST(CApi, G1PippengerWrappersHandleSinglePointAndInvalidScratch)
{
    const BoundaryG1 point = p1_mul(9);
    BoundaryG1Affine affine = {};
    blst_p1_to_affine(&affine, &point);
    const BoundaryG1Affine* single_point[1] = {&affine};
    const std::vector<unsigned char> single_scalar_blob = packed_scalars({7});
    const unsigned char* single_scalar[1] = {single_scalar_blob.data()};

    BoundaryG1 single = {};
    BoundaryG1 expected_single = {};
    blst_p1s_mult_pippenger(&single, single_point, 1, single_scalar, 8, nullptr);
    const auto scalar = scalar_bytes(7);
    blst_p1_mult(&expected_single, &point, scalar.data(), 8);
    EXPECT_TRUE(blst_p1_is_equal(&single, &expected_single));

    std::array<BoundaryG1, 2> points = {p1_mul(9), p1_mul(11)};
    std::array<BoundaryG1Affine, 2> affines = {};
    const BoundaryG1* point_ptrs[2] = {points.data(), nullptr};
    blst_p1s_to_affine(affines.data(), point_ptrs, affines.size());
    const BoundaryG1Affine* affine_ptrs[2] = {affines.data(), nullptr};
    const std::vector<unsigned char> scalar_blob = packed_scalars({7, 3});
    const unsigned char* scalar_ptrs[2] = {scalar_blob.data(), nullptr};

    BoundaryG1 rejected = p1_mul(1);
    blst_p1s_mult_pippenger(&rejected, affine_ptrs, affines.size(), scalar_ptrs, 8, nullptr);
    EXPECT_TRUE(blst_p1_is_inf(&rejected));

    blst::limb_t scratch[1] = {};
    BoundaryG1 tiled = p1_mul(1);
    blst_p1s_tile_pippenger(&tiled, affine_ptrs, affines.size(), scalar_ptrs, 8, scratch, 0, 0);
    EXPECT_TRUE(blst_p1_is_inf(&tiled));
}

TEST(CApi, G1MultipointWbitsMatchesSinglePointOps)
{
    std::array<BoundaryG1, 3> points = {p1_mul(3), p1_mul(5), p1_mul(7)};
    std::array<BoundaryG1Affine, 3> affines = {};
    const BoundaryG1* point_ptrs[2] = {points.data(), nullptr};
    blst_p1s_to_affine(affines.data(), point_ptrs, points.size());
    const BoundaryG1Affine* affine_ptrs[2] = {affines.data(), nullptr};
    const std::vector<unsigned char> scalar_blob = packed_scalars({2, 3, 5});
    const unsigned char* scalar_ptrs[2] = {scalar_blob.data(), nullptr};

    BoundaryG1 expected_msm = {};
    BoundaryG1 term = {};
    const auto s0 = scalar_bytes(2);
    blst_p1_mult(&expected_msm, &points[0], s0.data(), 8);
    const auto s1 = scalar_bytes(3);
    blst_p1_mult(&term, &points[1], s1.data(), 8);
    blst_p1_add_or_double(&expected_msm, &expected_msm, &term);
    const auto s2 = scalar_bytes(5);
    blst_p1_mult(&term, &points[2], s2.data(), 8);
    blst_p1_add_or_double(&expected_msm, &expected_msm, &term);

    const size_t wbits = 4;
    std::vector<BoundaryG1Affine> table(
        blst_p1s_mult_wbits_precompute_sizeof(wbits, affines.size()) / sizeof(BoundaryG1Affine));
    std::vector<blst::limb_t> wbits_scratch(
        blst_p1s_mult_wbits_scratch_sizeof(affines.size()) / sizeof(blst::limb_t));
    blst_p1s_mult_wbits_precompute(table.data(), wbits, affine_ptrs, affines.size());

    BoundaryG1 wbits_result = {};
    blst_p1s_mult_wbits(
        &wbits_result,
        table.data(),
        wbits,
        affines.size(),
        scalar_ptrs,
        8,
        wbits_scratch.data());
    EXPECT_TRUE(blst_p1_is_equal(&wbits_result, &expected_msm));
}

TEST(CApi, G1MultipointWbitsHighWindowMatchesSinglePointOps)
{
    std::array<BoundaryG1, 5> points = {p1_mul(3), p1_mul(5), p1_mul(7), p1_mul(9), p1_mul(11)};
    std::array<BoundaryG1Affine, 5> affines = {};
    const BoundaryG1* point_ptrs[2] = {points.data(), nullptr};
    blst_p1s_to_affine(affines.data(), point_ptrs, points.size());
    const BoundaryG1Affine* affine_ptrs[2] = {affines.data(), nullptr};
    const std::vector<unsigned char> scalar_blob = packed_scalars({2, 3, 5, 7, 11});
    const unsigned char* scalar_ptrs[2] = {scalar_blob.data(), nullptr};

    BoundaryG1 expected_msm = {};
    bool first = true;
    for (size_t i = 0; i < points.size(); ++i) {
        BoundaryG1 term = {};
        const auto scalar = scalar_bytes(static_cast<std::uint64_t>(scalar_blob[i]));
        blst_p1_mult(&term, &points[i], scalar.data(), 8);
        if (first) {
            expected_msm = term;
            first = false;
        } else {
            blst_p1_add_or_double(&expected_msm, &expected_msm, &term);
        }
    }

    const size_t wbits = 10;
    std::vector<BoundaryG1Affine> table(
        blst_p1s_mult_wbits_precompute_sizeof(wbits, affines.size()) / sizeof(BoundaryG1Affine));
    std::vector<blst::limb_t> scratch(
        blst_p1s_mult_wbits_scratch_sizeof(affines.size()) / sizeof(blst::limb_t));
    blst_p1s_mult_wbits_precompute(table.data(), wbits, affine_ptrs, affines.size());

    BoundaryG1 out = {};
    blst_p1s_mult_wbits(&out, table.data(), wbits, affines.size(), scalar_ptrs, 8, scratch.data());
    EXPECT_TRUE(blst_p1_is_equal(&out, &expected_msm));
}

TEST(CApi, ScalarG1WbitsPrecomputeUsesBoundedWorkspace)
{
    run_scalar_wbits_precompute_case<blst::Fp>(
        9, PrecomputePointLayout::packed);
    run_scalar_wbits_precompute_case<blst::Fp>(
        14, PrecomputePointLayout::pointer_array);
    run_scalar_wbits_precompute_case<blst::Fp>(
        12, PrecomputePointLayout::aliased_packed);
}

TEST(CApi, ScalarG2WbitsPrecomputeUsesBoundedWorkspace)
{
    run_scalar_wbits_precompute_case<blst::Fp2>(
        9, PrecomputePointLayout::packed);
    run_scalar_wbits_precompute_case<blst::Fp2>(
        14, PrecomputePointLayout::pointer_array);
    run_scalar_wbits_precompute_case<blst::Fp2>(
        11, PrecomputePointLayout::aliased_packed);
}

TEST(CApi, G1MultipointTileMatchesSinglePointOps)
{
    auto run_case = [](const std::array<BoundaryG1, 3>& points,
                       const std::vector<unsigned char>& scalar_blob,
                       size_t nbits,
                       const BoundaryG1& expected_msm) {
        std::array<BoundaryG1Affine, 3> affines = {};
        const BoundaryG1* point_ptrs[2] = {points.data(), nullptr};
        blst_p1s_to_affine(affines.data(), point_ptrs, points.size());
        const BoundaryG1Affine* affine_ptrs[2] = {affines.data(), nullptr};
        const unsigned char* scalar_ptrs[2] = {scalar_blob.data(), nullptr};

        const size_t window = 4;
        const size_t bucket_words =
            blst_p1s_mult_pippenger_window_scratch_sizeof(window) /
            sizeof(blst::limb_t);
        std::vector<blst::limb_t> scratch(bucket_words);

        BoundaryG1 tiled = {};
        bool first = true;
        for (size_t bit0 = ((nbits - 1) / window) * window;;) {
            BoundaryG1 tile = {};
            std::fill(scratch.begin(), scratch.end(), 0);
            blst_p1s_tile_pippenger(
                &tile,
                affine_ptrs,
                affines.size(),
                scalar_ptrs,
                nbits,
                scratch.data(),
                bit0,
                window);
            if (first) {
                tiled = tile;
                first = false;
            } else {
                for (size_t i = 0; i < window; ++i) {
                    blst_p1_double(&tiled, &tiled);
                }
                blst_p1_add_or_double(&tiled, &tiled, &tile);
            }
            if (bit0 == 0) {
                break;
            }
            bit0 -= window;
        }

        EXPECT_TRUE(blst_p1_is_equal(&tiled, &expected_msm));
    };

    {
        std::array<BoundaryG1, 3> points = {p1_mul(3), p1_mul(5), p1_mul(7)};
        const std::vector<unsigned char> scalar_blob = packed_scalars({2, 3, 5});
        BoundaryG1 expected_msm = {};
        BoundaryG1 term = {};
        const auto s0 = scalar_bytes(2);
        blst_p1_mult(&expected_msm, &points[0], s0.data(), 8);
        const auto s1 = scalar_bytes(3);
        blst_p1_mult(&term, &points[1], s1.data(), 8);
        blst_p1_add_or_double(&expected_msm, &expected_msm, &term);
        const auto s2 = scalar_bytes(5);
        blst_p1_mult(&term, &points[2], s2.data(), 8);
        blst_p1_add_or_double(&expected_msm, &expected_msm, &term);
        run_case(points, scalar_blob, 8, expected_msm);
    }

    {
        std::array<BoundaryG1, 3> points = {p1_mul(3), p1_mul(5), p1_mul(7)};
        const std::array<std::uint16_t, 3> scalar_values = {0x101, 0x155, 0x0ab};
        std::vector<unsigned char> scalar_blob = {
            static_cast<unsigned char>(scalar_values[0] & 0xff),
            static_cast<unsigned char>(scalar_values[0] >> 8),
            static_cast<unsigned char>(scalar_values[1] & 0xff),
            static_cast<unsigned char>(scalar_values[1] >> 8),
            static_cast<unsigned char>(scalar_values[2] & 0xff),
            static_cast<unsigned char>(scalar_values[2] >> 8),
        };

        BoundaryG1 expected_msm = {};
        BoundaryG1 term = {};
        const auto s0 = scalar_bytes(scalar_values[0]);
        blst_p1_mult(&expected_msm, &points[0], s0.data(), 10);
        const auto s1 = scalar_bytes(scalar_values[1]);
        blst_p1_mult(&term, &points[1], s1.data(), 10);
        blst_p1_add_or_double(&expected_msm, &expected_msm, &term);
        const auto s2 = scalar_bytes(scalar_values[2]);
        blst_p1_mult(&term, &points[2], s2.data(), 10);
        blst_p1_add_or_double(&expected_msm, &expected_msm, &term);
        run_case(points, scalar_blob, 10, expected_msm);
    }
}

TEST(CApi, G2MultipointWrappersMatchSinglePointOps)
{
    std::array<BoundaryG2, 3> points = {p2_mul(11), p2_mul(13), p2_mul(17)};
    std::array<BoundaryG2Affine, 3> affines = {};
    const BoundaryG2* point_ptrs[2] = {points.data(), nullptr};
    const BoundaryG2Affine* affine_ptrs[2] = {affines.data(), nullptr};

    blst_p2s_to_affine(affines.data(), point_ptrs, points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        BoundaryG2Affine expected = {};
        blst_p2_to_affine(&expected, &points[i]);
        EXPECT_TRUE(blst_p2_affine_is_equal(&affines[i], &expected));
    }

    BoundaryG2 sum = {};
    blst_p2s_add(&sum, affine_ptrs, affines.size());

    BoundaryG2 expected_sum = points[0];
    blst_p2_add_or_double(&expected_sum, &expected_sum, &points[1]);
    blst_p2_add_or_double(&expected_sum, &expected_sum, &points[2]);
    EXPECT_TRUE(blst_p2_is_equal(&sum, &expected_sum));

    const std::vector<unsigned char> scalar_blob = packed_scalars({2, 3, 5});
    const unsigned char* scalar_ptrs[2] = {scalar_blob.data(), nullptr};
    std::vector<blst::limb_t> scratch(
        blst_p2s_mult_pippenger_scratch_sizeof(affines.size()) / sizeof(blst::limb_t));

    BoundaryG2 pippenger = {};
    blst_p2s_mult_pippenger(&pippenger, affine_ptrs, affines.size(), scalar_ptrs, 8, scratch.data());

    BoundaryG2 expected_msm = {};
    BoundaryG2 term = {};
    const auto s0 = scalar_bytes(2);
    blst_p2_mult(&expected_msm, &points[0], s0.data(), 8);
    const auto s1 = scalar_bytes(3);
    blst_p2_mult(&term, &points[1], s1.data(), 8);
    blst_p2_add_or_double(&expected_msm, &expected_msm, &term);
    const auto s2 = scalar_bytes(5);
    blst_p2_mult(&term, &points[2], s2.data(), 8);
    blst_p2_add_or_double(&expected_msm, &expected_msm, &term);
    EXPECT_TRUE(blst_p2_is_equal(&pippenger, &expected_msm));
}

TEST(CApi, G2PippengerWrappersHandleSinglePointAndInvalidScratch)
{
    const BoundaryG2 point = p2_mul(9);
    BoundaryG2Affine affine = {};
    blst_p2_to_affine(&affine, &point);
    const BoundaryG2Affine* single_point[1] = {&affine};
    const std::vector<unsigned char> single_scalar_blob = packed_scalars({7});
    const unsigned char* single_scalar[1] = {single_scalar_blob.data()};

    BoundaryG2 single = {};
    BoundaryG2 expected_single = {};
    blst_p2s_mult_pippenger(&single, single_point, 1, single_scalar, 8, nullptr);
    const auto scalar = scalar_bytes(7);
    blst_p2_mult(&expected_single, &point, scalar.data(), 8);
    EXPECT_TRUE(blst_p2_is_equal(&single, &expected_single));

    std::array<BoundaryG2, 2> points = {p2_mul(9), p2_mul(11)};
    std::array<BoundaryG2Affine, 2> affines = {};
    const BoundaryG2* point_ptrs[2] = {points.data(), nullptr};
    blst_p2s_to_affine(affines.data(), point_ptrs, affines.size());
    const BoundaryG2Affine* affine_ptrs[2] = {affines.data(), nullptr};
    const std::vector<unsigned char> scalar_blob = packed_scalars({7, 3});
    const unsigned char* scalar_ptrs[2] = {scalar_blob.data(), nullptr};

    BoundaryG2 rejected = p2_mul(1);
    blst_p2s_mult_pippenger(&rejected, affine_ptrs, affines.size(), scalar_ptrs, 8, nullptr);
    EXPECT_TRUE(blst_p2_is_inf(&rejected));

    blst::limb_t scratch[1] = {};
    BoundaryG2 tiled = p2_mul(1);
    blst_p2s_tile_pippenger(&tiled, affine_ptrs, affines.size(), scalar_ptrs, 8, scratch, 0, 0);
    EXPECT_TRUE(blst_p2_is_inf(&tiled));
}

TEST(CApi, PairingPrimitivesMatchBackend)
{
    const BoundaryG1 p = p1_mul(3);
    const BoundaryG2 q = p2_mul(11);

    BoundaryG1Affine p_affine = {};
    BoundaryG2Affine q_affine = {};
    blst_p1_to_affine(&p_affine, &p);
    blst_p2_to_affine(&q_affine, &q);

    BoundaryGt miller = {};
    BoundaryGt final_exp = {};
    BoundaryGt one = *blst_fp12_one();
    BoundaryGt product = {};

    blst_miller_loop(&miller, &q_affine, &p_affine);
    blst_final_exp(&final_exp, &miller);
    blst_fp12_mul(&product, &final_exp, &one);

    EXPECT_TRUE(blst_fp12_is_equal(&final_exp, &product));
    EXPECT_TRUE(blst_fp12_is_one(&one));
    EXPECT_TRUE(blst_fp12_in_group(&final_exp));
    EXPECT_TRUE(blst_fp12_finalverify(&final_exp, &final_exp));

    PublicMillerResult expected_miller = {};
    PublicGt expected_public = {};
    const PublicG2Affine public_q = to_public(q_affine);
    const PublicG1Affine public_p = to_public(p_affine);
    PublicGt::miller_loop(expected_miller, &public_q, &public_p);
    PublicGt::final_exp(expected_public, expected_miller);
    const BoundaryGt expected = gt_to_boundary(expected_public);

    EXPECT_TRUE(blst::Fp12::equal(final_exp, expected));
}

TEST(CApi, PointSubtractionAndAffineNegationMatchExpectedMultiples)
{
    const BoundaryG1 p1_eight = p1_mul(8);
    const BoundaryG1 p1_three = p1_mul(3);
    const BoundaryG1 p1_five = p1_mul(5);
    BoundaryG1Affine p1_three_affine = {};
    BoundaryG1 p1_difference = {};
    BoundaryG1 p1_mixed_difference = {};
    blst_p1_to_affine(&p1_three_affine, &p1_three);
    blst_p1_sub(&p1_difference, &p1_eight, &p1_three);
    blst_p1_sub_affine(&p1_mixed_difference, &p1_eight, &p1_three_affine);
    EXPECT_TRUE(blst_p1_is_equal(&p1_difference, &p1_five));
    EXPECT_TRUE(blst_p1_is_equal(&p1_mixed_difference, &p1_five));

    BoundaryG1 p1_alias_left = p1_eight;
    BoundaryG1 p1_alias_right = p1_three;
    blst_p1_sub(&p1_alias_left, &p1_alias_left, &p1_three);
    blst_p1_sub(&p1_alias_right, &p1_eight, &p1_alias_right);
    EXPECT_TRUE(blst_p1_is_equal(&p1_alias_left, &p1_five));
    EXPECT_TRUE(blst_p1_is_equal(&p1_alias_right, &p1_five));

    BoundaryG1 p1_equal_difference = {};
    blst_p1_sub(&p1_equal_difference, &p1_three, &p1_three);
    EXPECT_TRUE(blst_p1_is_inf(&p1_equal_difference));

    const BoundaryG1Affine p1_original_affine = p1_three_affine;
    blst_p1_affine_cneg(&p1_three_affine, false);
    EXPECT_TRUE(blst_p1_affine_is_equal(&p1_three_affine, &p1_original_affine));
    blst_p1_affine_cneg(&p1_three_affine, true);
    BoundaryG1 p1_negated = {};
    BoundaryG1 p1_identity = {};
    blst_p1_from_affine(&p1_negated, &p1_three_affine);
    blst_p1_add_or_double(&p1_identity, &p1_three, &p1_negated);
    EXPECT_TRUE(blst_p1_is_inf(&p1_identity));

    const BoundaryG2 p2_eight = p2_mul(8);
    const BoundaryG2 p2_three = p2_mul(3);
    const BoundaryG2 p2_five = p2_mul(5);
    BoundaryG2Affine p2_three_affine = {};
    BoundaryG2 p2_difference = {};
    BoundaryG2 p2_mixed_difference = {};
    blst_p2_to_affine(&p2_three_affine, &p2_three);
    blst_p2_sub(&p2_difference, &p2_eight, &p2_three);
    blst_p2_sub_affine(&p2_mixed_difference, &p2_eight, &p2_three_affine);
    EXPECT_TRUE(blst_p2_is_equal(&p2_difference, &p2_five));
    EXPECT_TRUE(blst_p2_is_equal(&p2_mixed_difference, &p2_five));

    BoundaryG2 p2_alias_left = p2_eight;
    BoundaryG2 p2_alias_right = p2_three;
    blst_p2_sub(&p2_alias_left, &p2_alias_left, &p2_three);
    blst_p2_sub(&p2_alias_right, &p2_eight, &p2_alias_right);
    EXPECT_TRUE(blst_p2_is_equal(&p2_alias_left, &p2_five));
    EXPECT_TRUE(blst_p2_is_equal(&p2_alias_right, &p2_five));

    BoundaryG2 p2_equal_difference = {};
    blst_p2_sub(&p2_equal_difference, &p2_three, &p2_three);
    EXPECT_TRUE(blst_p2_is_inf(&p2_equal_difference));

    const BoundaryG2Affine p2_original_affine = p2_three_affine;
    blst_p2_affine_cneg(&p2_three_affine, false);
    EXPECT_TRUE(blst_p2_affine_is_equal(&p2_three_affine, &p2_original_affine));
    blst_p2_affine_cneg(&p2_three_affine, true);
    BoundaryG2 p2_negated = {};
    BoundaryG2 p2_identity = {};
    blst_p2_from_affine(&p2_negated, &p2_three_affine);
    blst_p2_add_or_double(&p2_identity, &p2_three, &p2_negated);
    EXPECT_TRUE(blst_p2_is_inf(&p2_identity));
}

TEST(CApi, TargetGroupScalarMultiplicationMatchesPairingBilinearity)
{
    const BoundaryG1 p = p1_mul(1);
    const BoundaryG1 p_times_five = p1_mul(5);
    const BoundaryG2 q = p2_mul(1);
    BoundaryG1Affine p_affine = {};
    BoundaryG1Affine p_times_five_affine = {};
    BoundaryG2Affine q_affine = {};
    blst_p1_to_affine(&p_affine, &p);
    blst_p1_to_affine(&p_times_five_affine, &p_times_five);
    blst_p2_to_affine(&q_affine, &q);

    BoundaryGt miller = {};
    BoundaryGt base = {};
    BoundaryGt expected_miller = {};
    BoundaryGt expected = {};
    blst_miller_loop(&miller, &q_affine, &p_affine);
    blst_final_exp(&base, &miller);
    blst_miller_loop(&expected_miller, &q_affine, &p_times_five_affine);
    blst_final_exp(&expected, &expected_miller);

    EXPECT_TRUE(blst_gt_is_equal(blst_gt_generator(), &base));
    EXPECT_TRUE(blst_gt_is_identity(blst_gt_identity()));
    EXPECT_FALSE(blst_gt_is_identity(&base));

    BoundaryGt doubled = {};
    BoundaryGt added = {};
    BoundaryGt negated = {};
    BoundaryGt cancelled = {};
    BoundaryGt difference = {};
    blst_gt_double(&doubled, &base);
    blst_gt_add(&added, &base, &base);
    blst_gt_neg(&negated, &base);
    blst_gt_add(&cancelled, &base, &negated);
    blst_gt_sub(&difference, &doubled, &base);
    EXPECT_TRUE(blst_gt_is_equal(&doubled, &added));
    EXPECT_TRUE(blst_gt_is_identity(&cancelled));
    EXPECT_TRUE(blst_gt_is_equal(&difference, &base));

    BoundaryGt alias = base;
    blst_gt_add(&alias, &alias, &base);
    EXPECT_TRUE(blst_gt_is_equal(&alias, &doubled));
    blst_gt_sub(&alias, &alias, &base);
    EXPECT_TRUE(blst_gt_is_equal(&alias, &base));
    blst_gt_neg(&alias, &alias);
    EXPECT_TRUE(blst_gt_is_equal(&alias, &negated));
    blst_gt_double(&alias, &base);
    EXPECT_TRUE(blst_gt_is_equal(&alias, &doubled));

    const auto five = scalar_bytes(5);
    std::array<unsigned char, 33> unaligned_five = {};
    std::copy(five.begin(), five.end(), unaligned_five.begin() + 1);
    BoundaryGt short_scalar = {};
    BoundaryGt full_scalar = {};
    BoundaryGt unaligned_full_scalar = {};
    blst_gt_mult(&short_scalar, &base, five.data(), 3);
    blst_gt_mult(&full_scalar, &base, five.data(), 256);
    blst_gt_mult(&unaligned_full_scalar, &base, unaligned_five.data() + 1, 256);
    EXPECT_TRUE(blst_fp12_is_equal(&short_scalar, &expected));
    EXPECT_TRUE(blst_fp12_is_equal(&full_scalar, &expected));
    EXPECT_TRUE(blst_fp12_is_equal(&unaligned_full_scalar, &expected));

    const auto zero = scalar_bytes(0);
    const auto one = scalar_bytes(1);
    BoundaryGt zero_product = {};
    BoundaryGt alias_product = base;
    blst_gt_mult(&zero_product, &base, zero.data(), 256);
    blst_gt_mult(&alias_product, &alias_product, one.data(), 256);
    EXPECT_TRUE(blst_fp12_is_one(&zero_product));
    EXPECT_TRUE(blst_fp12_is_equal(&alias_product, &base));

    constexpr std::array<unsigned char, 32> scalar_modulus_minus_one = {
        0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff,
        0xfe, 0x5b, 0xfe, 0xff, 0x02, 0xa4, 0xbd, 0x53,
        0x05, 0xd8, 0xa1, 0x09, 0x08, 0xd8, 0x39, 0x33,
        0x48, 0x7d, 0x9d, 0x29, 0x53, 0xa7, 0xed, 0x73,
    };
    BoundaryGt modulus_minus_one_product = {};
    BoundaryGt inverse = base;
    blst_gt_mult(
        &modulus_minus_one_product,
        &base,
        scalar_modulus_minus_one.data(),
        256);
    blst_fp12_conjugate(&inverse);
    EXPECT_TRUE(blst_fp12_is_equal(&modulus_minus_one_product, &inverse));

    constexpr std::uint64_t deterministic_scalar = 0xdeadbeef;
    const BoundaryG1 p_deterministic = p1_mul(deterministic_scalar);
    BoundaryG1Affine p_deterministic_affine = {};
    BoundaryGt deterministic_miller = {};
    BoundaryGt deterministic_expected = {};
    BoundaryGt deterministic_product = {};
    blst_p1_to_affine(&p_deterministic_affine, &p_deterministic);
    blst_miller_loop(&deterministic_miller, &q_affine, &p_deterministic_affine);
    blst_final_exp(&deterministic_expected, &deterministic_miller);
    const auto deterministic_bytes = scalar_bytes(deterministic_scalar);
    blst_gt_mult(&deterministic_product, &base, deterministic_bytes.data(), 256);
    EXPECT_TRUE(blst_fp12_is_equal(&deterministic_product, &deterministic_expected));
}

TEST(CApi, GroupConstantTimeAndVariableTimeMultiplicationAgree)
{
    const std::array<size_t, 6> bit_lengths = {0, 1, 3, 128, 255, 256};
    for (const size_t nbits : bit_lengths) {
        std::array<unsigned char, 32> scalar = {};
        for (size_t i = 0; i < scalar.size(); ++i) {
            scalar[i] = static_cast<unsigned char>(0xa5U ^ static_cast<unsigned>(i * 29U));
        }
        if (nbits == 0) {
            scalar.fill(0);
        } else {
            const size_t nbytes = (nbits + 7) / 8;
            std::fill(scalar.begin() + static_cast<std::ptrdiff_t>(nbytes), scalar.end(), 0);
            if ((nbits & 7U) != 0) {
                scalar[nbytes - 1] &= static_cast<unsigned char>((1U << (nbits & 7U)) - 1U);
            }
        }

        BoundaryG1 g1_ct = {};
        BoundaryG1 g1_vartime = {};
        blst_p1_mult(&g1_ct, blst_p1_generator(), scalar.data(), nbits);
        blst_p1_mult_vartime(&g1_vartime, blst_p1_generator(), scalar.data(), nbits);
        EXPECT_TRUE(blst_p1_is_equal(&g1_ct, &g1_vartime)) << "G1 nbits=" << nbits;

        BoundaryG2 g2_ct = {};
        BoundaryG2 g2_vartime = {};
        blst_p2_mult(&g2_ct, blst_p2_generator(), scalar.data(), nbits);
        blst_p2_mult_vartime(&g2_vartime, blst_p2_generator(), scalar.data(), nbits);
        EXPECT_TRUE(blst_p2_is_equal(&g2_ct, &g2_vartime)) << "G2 nbits=" << nbits;

        BoundaryGt gt_ct = {};
        BoundaryGt gt_vartime = {};
        blst_gt_mult(&gt_ct, blst_gt_generator(), scalar.data(), nbits);
        blst_gt_mult_vartime(&gt_vartime, blst_gt_generator(), scalar.data(), nbits);
        EXPECT_TRUE(blst_gt_is_equal(&gt_ct, &gt_vartime)) << "GT nbits=" << nbits;
    }
}

TEST(CApi, G1TileUsesTheExactBitOffset)
{
    std::array<BoundaryG1, 2> points = {p1_mul(3), p1_mul(5)};
    std::array<BoundaryG1Affine, 2> affines = {};
    blst_p1_to_affine(&affines[0], &points[0]);
    blst_p1_to_affine(&affines[1], &points[1]);
    const BoundaryG1Affine* point_ptrs[2] = {&affines[0], &affines[1]};
    const std::array<unsigned char, 2> scalar_values = {64, 0};
    const unsigned char* scalar_ptrs[2] = {&scalar_values[0], &scalar_values[1]};
    const size_t window = 4;
    std::vector<blst::limb_t> scratch(
        blst_p1s_mult_pippenger_window_scratch_sizeof(window) / sizeof(blst::limb_t));

    BoundaryG1 out = {};
    blst_p1s_tile_pippenger(
        &out, point_ptrs, points.size(), scalar_ptrs, 8,
        scratch.data(), 3, window);

    const auto eight = scalar_bytes(8);
    BoundaryG1 expected = {};
    blst_p1_mult(&expected, &points[0], eight.data(), 8);
    blst_p1_cneg(&expected, true);
    EXPECT_TRUE(blst_p1_is_equal(&out, &expected));
}

TEST(CApi, MillerLoopNMatchesNaiveProduct)
{
    std::array<BoundaryG1, 3> ps = {p1_mul(3), p1_mul(5), p1_mul(7)};
    std::array<BoundaryG2, 3> qs = {p2_mul(11), p2_mul(13), p2_mul(17)};
    std::array<BoundaryG1Affine, 3> ps_affine = {};
    std::array<BoundaryG2Affine, 3> qs_affine = {};
    const BoundaryG1* p_ptrs[2] = {ps.data(), nullptr};
    const BoundaryG2* q_ptrs[2] = {qs.data(), nullptr};
    const BoundaryG1Affine* ps_affine_ptrs[2] = {ps_affine.data(), nullptr};
    const BoundaryG2Affine* qs_affine_ptrs[2] = {qs_affine.data(), nullptr};

    blst_p1s_to_affine(ps_affine.data(), p_ptrs, ps.size());
    blst_p2s_to_affine(qs_affine.data(), q_ptrs, qs.size());

    BoundaryGt combined = {};
    blst_miller_loop_n(&combined, qs_affine_ptrs, ps_affine_ptrs, ps.size());

    BoundaryGt naive = *blst_fp12_one();
    for (size_t i = 0; i < ps.size(); ++i) {
        BoundaryGt term = {};
        blst_miller_loop(&term, &qs_affine[i], &ps_affine[i]);
        blst_fp12_mul(&naive, &naive, &term);
    }

    EXPECT_TRUE(blst_fp12_is_equal(&combined, &naive));
}

TEST(CApi, MillerLoopNMatchesNaiveProductAtBatchBoundaries)
{
    constexpr size_t npoints = 32;
    std::array<BoundaryG1, npoints> ps = {};
    std::array<BoundaryG2, npoints> qs = {};
    std::array<BoundaryG1Affine, npoints> ps_affine = {};
    std::array<BoundaryG2Affine, npoints> qs_affine = {};
    std::array<const BoundaryG1Affine*, npoints> ps_explicit = {};
    std::array<const BoundaryG2Affine*, npoints> qs_explicit = {};

    for (size_t i = 0; i < npoints; ++i) {
        ps[i] = p1_mul(3 + i);
        qs[i] = p2_mul(11 + i);
        ps_explicit[i] = &ps_affine[i];
        qs_explicit[i] = &qs_affine[i];
    }

    const BoundaryG1* p_ptrs[2] = {ps.data(), nullptr};
    const BoundaryG2* q_ptrs[2] = {qs.data(), nullptr};
    const BoundaryG1Affine* ps_affine_ptrs[2] = {ps_affine.data(), nullptr};
    const BoundaryG2Affine* qs_affine_ptrs[2] = {qs_affine.data(), nullptr};

    blst_p1s_to_affine(ps_affine.data(), p_ptrs, ps.size());
    blst_p2s_to_affine(qs_affine.data(), q_ptrs, qs.size());

    BoundaryGt naive = *blst_fp12_one();
    for (size_t i = 0; i < npoints; ++i) {
        BoundaryGt term = {};
        blst_miller_loop(&term, &qs_affine[i], &ps_affine[i]);
        blst_fp12_mul(&naive, &naive, &term);

        const size_t count = i + 1;
        if (count == 1 || count == 15 || count == 16 || count == 17 || count == 32) {
            BoundaryGt compact = {};
            BoundaryGt explicit_pointers = {};
            blst_miller_loop_n(&compact, qs_affine_ptrs, ps_affine_ptrs, count);
            blst_miller_loop_n(
                &explicit_pointers,
                qs_explicit.data(),
                ps_explicit.data(),
                count);
            EXPECT_TRUE(blst_fp12_is_equal(&compact, &naive)) << count;
            EXPECT_TRUE(blst_fp12_is_equal(&explicit_pointers, &naive)) << count;
        }
    }
}

TEST(CApi, MillerLoopNSkipsInfinityPairs)
{
    constexpr size_t npoints = 17;
    std::array<BoundaryG1, npoints> ps = {};
    std::array<BoundaryG2, npoints> qs = {};
    std::array<BoundaryG1Affine, npoints> ps_affine = {};
    std::array<BoundaryG2Affine, npoints> qs_affine = {};
    for (size_t i = 0; i < npoints; ++i) {
        ps[i] = p1_mul(3 + i);
        qs[i] = p2_mul(11 + i);
    }

    const BoundaryG1* p_ptrs[2] = {ps.data(), nullptr};
    const BoundaryG2* q_ptrs[2] = {qs.data(), nullptr};
    const BoundaryG1Affine* ps_affine_ptrs[2] = {ps_affine.data(), nullptr};
    const BoundaryG2Affine* qs_affine_ptrs[2] = {qs_affine.data(), nullptr};
    blst_p1s_to_affine(ps_affine.data(), p_ptrs, ps.size());
    blst_p2s_to_affine(qs_affine.data(), q_ptrs, qs.size());
    ps_affine[7] = {};
    ASSERT_TRUE(blst_p1_affine_is_inf(&ps_affine[7]));

    BoundaryGt combined = {};
    blst_miller_loop_n(&combined, qs_affine_ptrs, ps_affine_ptrs, npoints);

    BoundaryGt naive = *blst_fp12_one();
    for (size_t i = 0; i < npoints; ++i) {
        BoundaryGt term = {};
        blst_miller_loop(&term, &qs_affine[i], &ps_affine[i]);
        blst_fp12_mul(&naive, &naive, &term);
    }
    EXPECT_TRUE(blst_fp12_is_equal(&combined, &naive));
}

TEST(CApi, MillerLoopNEmptyReturnsIdentityAndMalformedInputFailsClosed)
{
    BoundaryGt out = {};
    blst_miller_loop_n(&out, nullptr, nullptr, 0);
    EXPECT_TRUE(blst_fp12_is_one(&out));

    const BoundaryG1Affine* ps[1] = {nullptr};
    const BoundaryG2Affine* qs[1] = {nullptr};
    blst_miller_loop_n(&out, qs, ps, 1);
    const BoundaryGt zero = {};
    EXPECT_TRUE(blst_fp12_is_equal(&out, &zero));
    EXPECT_FALSE(blst_fp12_is_one(&out));
}

TEST(CApi, LinePrecomputeHelpersUseGeneratorsForNullInputs)
{
    BoundaryFp6 lines_from_null[68] = {};
    BoundaryFp6 lines_from_generator[68] = {};
    blst_precompute_lines(lines_from_null, nullptr);
    blst_precompute_lines(lines_from_generator, blst_p2_affine_generator());
    EXPECT_EQ(std::memcmp(lines_from_null, lines_from_generator, sizeof(lines_from_null)), 0);

    BoundaryGt null_loop = {};
    BoundaryGt explicit_loop = {};
    blst_miller_loop_lines(&null_loop, lines_from_null, nullptr);
    blst_miller_loop_lines(&explicit_loop, lines_from_generator, blst_p1_affine_generator());
    EXPECT_TRUE(blst_fp12_is_equal(&null_loop, &explicit_loop));
}

TEST(CApi, LinePrecomputeHelpersPreservePairingIdentity)
{
    const BoundaryG1Affine p_infinity = {};
    const BoundaryG2Affine q_infinity = {};

    BoundaryFp6 finite_lines[68] = {};
    blst_precompute_lines(finite_lines, blst_p2_affine_generator());
    BoundaryGt p_infinity_result = {};
    blst_miller_loop_lines(&p_infinity_result, finite_lines, &p_infinity);
    EXPECT_TRUE(blst_fp12_is_one(&p_infinity_result));

    BoundaryFp6 infinity_lines[68] = {};
    blst_precompute_lines(infinity_lines, &q_infinity);
    BoundaryGt q_infinity_result = {};
    blst_miller_loop_lines(
        &q_infinity_result, infinity_lines, blst_p1_affine_generator());
    EXPECT_TRUE(blst_fp12_is_one(&q_infinity_result));
}

TEST(CApi, ScalarPreparedG2BatchComposesAndRejectsInvalidRangesAtomically)
{
    constexpr size_t npoints = 37;
    std::array<BoundaryG1, npoints> ps = {};
    std::array<BoundaryG2, npoints> qs = {};
    std::array<BoundaryG1Affine, npoints> ps_affine = {};
    std::array<BoundaryG2Affine, npoints> qs_affine = {};
    std::array<size_t, npoints> indices = {};
    for (size_t i = 0; i < npoints; ++i) {
        ps[i] = p1_mul(3 + i);
        qs[i] = p2_mul(11 + i);
        indices[i] = (7 * i + 3) % npoints;
    }
    ps[9] = {};
    qs[13] = {};

    const BoundaryG1* p_ptrs[2] = {ps.data(), nullptr};
    const BoundaryG2* q_ptrs[2] = {qs.data(), nullptr};
    blst_p1s_to_affine(ps_affine.data(), p_ptrs, ps.size());
    blst_p2s_to_affine(qs_affine.data(), q_ptrs, qs.size());

    blst_scalar_prepared_g2_batch* batch =
        blst_scalar_prepared_g2_batch_new(qs.data(), qs.size());
    ASSERT_NE(batch, nullptr);
    EXPECT_EQ(blst_scalar_prepared_g2_batch_len(batch), npoints);
    EXPECT_GT(blst_scalar_prepared_g2_batch_storage_bytes(batch), 0u);

    blst::scalar::MillerResult scalar_expected =
        blst::scalar::MillerResult::identity();
    for (size_t i = 0; i < npoints; ++i) {
        blst::scalar::MillerResult term = {};
        blst::scalar::GT::miller_loop(
            term, &qs_affine[indices[i]], &ps_affine[i]);
        blst::scalar::MillerResult::mul(
            scalar_expected, scalar_expected, term);
    }
    const BoundaryGt expected = scalar_expected.data;

    BoundaryGt complete = {};
    ASSERT_TRUE(blst_scalar_prepared_g2_batch_miller_loop_product_indices(
        &complete, ps.data(), indices.data(), npoints, batch));
    EXPECT_TRUE(blst_fp12_is_equal(&complete, &expected));

    BoundaryGt left = {};
    BoundaryGt right = {};
    BoundaryGt split = {};
    ASSERT_TRUE(blst_scalar_prepared_g2_batch_miller_loop_product_indices(
        &left, ps.data(), indices.data(), 16, batch));
    ASSERT_TRUE(blst_scalar_prepared_g2_batch_miller_loop_product_indices(
        &right, ps.data() + 16, indices.data() + 16, npoints - 16, batch));
    blst_fp12_mul(&split, &left, &right);
    EXPECT_TRUE(blst_fp12_is_equal(&split, &expected));

    BoundaryGt empty = {};
    ASSERT_TRUE(blst_scalar_prepared_g2_batch_miller_loop_product_indices(
        &empty, nullptr, nullptr, 0, batch));
    EXPECT_TRUE(blst_fp12_is_one(&empty));

    constexpr size_t range_begin = 7;
    constexpr size_t range_count = 10;
    std::array<BoundaryGt, range_count> range = {};
    ASSERT_TRUE(blst_scalar_prepared_g2_batch_pairing_range(
        range.data(), range_begin, range.size(), &ps[0], batch));
    for (size_t i = 0; i < range.size(); ++i) {
        BoundaryGt miller = {};
        BoundaryGt expected_pairing = {};
        blst_miller_loop(
            &miller, &qs_affine[range_begin + i], &ps_affine[0]);
        blst_final_exp(&expected_pairing, &miller);
        EXPECT_TRUE(blst_fp12_is_equal(&range[i], &expected_pairing));
    }
    ASSERT_TRUE(blst_scalar_prepared_g2_batch_pairing_range(
        nullptr, npoints, 0, &ps[0], batch));

    const size_t bad_index = npoints;
    BoundaryGt unchanged = expected;
    EXPECT_FALSE(blst_scalar_prepared_g2_batch_miller_loop_product_indices(
        &unchanged, ps.data(), &bad_index, 1, batch));
    EXPECT_TRUE(blst_fp12_is_equal(&unchanged, &expected));
    EXPECT_FALSE(blst_scalar_prepared_g2_batch_miller_loop_product_indices(
        nullptr, ps.data(), indices.data(), npoints, batch));
    EXPECT_FALSE(blst_scalar_prepared_g2_batch_miller_loop_product_indices(
        &unchanged, ps.data(), indices.data(), npoints, nullptr));
    EXPECT_FALSE(blst_scalar_prepared_g2_batch_pairing_range(
        &unchanged, npoints, 1, &ps[0], batch));
    EXPECT_TRUE(blst_fp12_is_equal(&unchanged, &expected));

    blst_scalar_prepared_g2_batch_free(batch);
    blst_scalar_prepared_g2_batch_free(nullptr);
    EXPECT_EQ(blst_scalar_prepared_g2_batch_len(nullptr), 0u);
    EXPECT_EQ(blst_scalar_prepared_g2_batch_storage_bytes(nullptr), 0u);
    EXPECT_EQ(blst_scalar_prepared_g2_batch_new(nullptr, 1), nullptr);
    EXPECT_EQ(blst_scalar_prepared_g2_batch_new(qs.data(), 0), nullptr);
}

TEST(CApi, PairingWeightedAggregateRejectsMalformedScalarWithoutMutation)
{
    const std::array<unsigned char, 1> dst = {'D'};
    const std::array<unsigned char, 33> oversized_scalar = {};
    std::vector<uint64_t> raw((blst_pairing_sizeof() + sizeof(uint64_t) - 1) /
        sizeof(uint64_t));
    auto* ctx = reinterpret_cast<blst_pairing*>(raw.data());

    const auto expect_rejected_without_mutation = [&](const auto& aggregate) {
        blst_pairing_init(ctx, true, dst.data(), dst.size());
        const std::vector<uint64_t> initialized = raw;
        EXPECT_EQ(aggregate(), BLST_BAD_SCALAR);
        EXPECT_EQ(raw, initialized);
    };

    expect_rejected_without_mutation([&] {
        return blst_pairing_mul_n_aggregate_pk_in_g1(
            ctx, nullptr, nullptr, oversized_scalar.data(), 257,
            nullptr, 0, nullptr, 0);
    });
    expect_rejected_without_mutation([&] {
        return blst_pairing_chk_n_mul_n_aggr_pk_in_g1(
            ctx, nullptr, false, nullptr, false, nullptr, 1,
            nullptr, 0, nullptr, 0);
    });
    expect_rejected_without_mutation([&] {
        return blst_pairing_mul_n_aggregate_pk_in_g2(
            ctx, nullptr, nullptr, oversized_scalar.data(), 257,
            nullptr, 0, nullptr, 0);
    });
    expect_rejected_without_mutation([&] {
        return blst_pairing_chk_n_mul_n_aggr_pk_in_g2(
            ctx, nullptr, false, nullptr, false, nullptr, 1,
            nullptr, 0, nullptr, 0);
    });
}

#ifdef NDEBUG
TEST(CApi, FixedScalarCanonicalizerFailsClosedInRelease)
{
    blst::pow256 out = {};
    const std::array<unsigned char, 33> oversized_scalar = {};

    std::memset(out, 0xff, sizeof(out));
    blst::scalar::canonical_scalar_256(out, oversized_scalar.data(), 257);
    EXPECT_TRUE(std::all_of(std::begin(out), std::end(out),
        [](unsigned char value) { return value == 0; }));

    std::memset(out, 0xff, sizeof(out));
    blst::scalar::canonical_scalar_256(out, nullptr, 1);
    EXPECT_TRUE(std::all_of(std::begin(out), std::end(out),
        [](unsigned char value) { return value == 0; }));
}
#endif

TEST(CApi, PairingAggregateAllowsNullSignatureAfterFirstElement)
{
    const auto dst = std::to_array("BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_");
    const auto sk0 = scalar_bytes(3);
    const auto sk1 = scalar_bytes(5);
    const BoundaryG1 pk0 = p1_mul(3);
    const BoundaryG1 pk1 = p1_mul(5);
    const std::array<unsigned char, 3> msg0 = {'a', 'b', 'c'};
    const std::array<unsigned char, 3> msg1 = {'d', 'e', 'f'};
    BoundaryG2 hash0 = {};
    BoundaryG2 hash1 = {};
    BoundaryG2 sig0 = {};
    BoundaryG2 sig1 = {};
    BoundaryG2 agg_sig = {};
    blst_hash_to_g2(&hash0, msg0.data(), msg0.size(), reinterpret_cast<const unsigned char*>(dst.data()), dst.size() - 1, nullptr, 0);
    blst_hash_to_g2(&hash1, msg1.data(), msg1.size(), reinterpret_cast<const unsigned char*>(dst.data()), dst.size() - 1, nullptr, 0);
    blst_p2_mult(&sig0, &hash0, sk0.data(), 255);
    blst_p2_mult(&sig1, &hash1, sk1.data(), 255);
    blst_p2_add_or_double(&agg_sig, &sig0, &sig1);
    BoundaryG1Affine pk0_affine = {};
    BoundaryG1Affine pk1_affine = {};
    BoundaryG2Affine sig_affine = {};
    blst_p1_to_affine(&pk0_affine, &pk0);
    blst_p1_to_affine(&pk1_affine, &pk1);
    blst_p2_to_affine(&sig_affine, &agg_sig);

    std::vector<uint64_t> raw(blst_pairing_sizeof() / sizeof(uint64_t));
    auto* ctx = reinterpret_cast<blst_pairing*>(raw.data());
    blst_pairing_init(ctx, true, reinterpret_cast<const unsigned char*>(dst.data()), dst.size() - 1);
    ASSERT_EQ(
        blst_pairing_chk_n_aggr_pk_in_g1(
            ctx,
            &pk0_affine,
            true,
            &sig_affine,
            true,
            msg0.data(),
            msg0.size(),
            nullptr,
            0),
        blst::BLST_SUCCESS);
    ASSERT_EQ(
        blst_pairing_chk_n_aggr_pk_in_g1(
            ctx,
            &pk1_affine,
            true,
            nullptr,
            false,
            msg1.data(),
            msg1.size(),
            nullptr,
            0),
        blst::BLST_SUCCESS);

    blst_pairing_commit(ctx);
    EXPECT_TRUE(blst_pairing_finalverify(ctx, nullptr));
}

TEST(CApi, Fp12UtilityWrappersMatchScalarBackend)
{
    const BoundaryGt a = *blst_gt_generator();
    BoundaryGt mul = {};
    BoundaryGt sqr = {};
    BoundaryGt cyclo = {};
    BoundaryGt inverse = {};
    BoundaryGt frob = {};
    BoundaryGt conjugated = a;

    blst_fp12_mul(&mul, &a, &a);
    blst_fp12_sqr(&sqr, &a);
    blst_fp12_cyclotomic_sqr(&cyclo, &a);
    blst_fp12_inverse(&inverse, &a);
    EXPECT_TRUE(blst_fp12_frobenius_map(&frob, &a, 2));
    blst_fp12_conjugate(&conjugated);

    BoundaryGt expected_mul = {};
    BoundaryGt expected_inverse = {};
    BoundaryGt expected_frob = {};
    BoundaryGt expected_conjugated = a;
    BoundaryFp6 line = {};
    BoundaryGt xy_mul = {};
    BoundaryGt expected_xy_mul = {};

    const uint64_t five_words[6] = {5, 0, 0, 0, 0, 0};
    const uint64_t seven_words[6] = {7, 0, 0, 0, 0, 0};
    const uint64_t eleven_words[6] = {11, 0, 0, 0, 0, 0};
    BoundaryFp five = {};
    BoundaryFp seven = {};
    BoundaryFp eleven = {};

    BoundaryGt::mul(expected_mul, a, a);
    BoundaryGt::inverse(expected_inverse, a);
    BoundaryGt::frobenius_map(expected_frob, a, 2);
    BoundaryGt::conjugate(expected_conjugated);
    blst_fp_from_uint64(&five, five_words);
    blst_fp_from_uint64(&seven, seven_words);
    blst_fp_from_uint64(&eleven, eleven_words);
    line[0][0] = five;
    line[0][1] = seven;
    line[2][0] = eleven;
    BoundaryGt::mul_by_xy00z0(expected_xy_mul, a, line);
    blst_fp12_mul_by_xy00z0(&xy_mul, &a, &line);

    EXPECT_TRUE(BoundaryGt::equal(mul, expected_mul));
    EXPECT_TRUE(BoundaryGt::equal(sqr, expected_mul));
    EXPECT_TRUE(BoundaryGt::equal(cyclo, expected_mul));
    EXPECT_TRUE(BoundaryGt::equal(inverse, expected_inverse));
    EXPECT_TRUE(BoundaryGt::equal(frob, expected_frob));
    EXPECT_TRUE(BoundaryGt::equal(conjugated, expected_conjugated));
    EXPECT_TRUE(BoundaryGt::equal(xy_mul, expected_xy_mul));

    BoundaryGt rejected = expected_mul;
    EXPECT_FALSE(blst_fp12_frobenius_map(&rejected, &a, 0));
    EXPECT_TRUE(BoundaryGt::equal(rejected, expected_mul));
    EXPECT_FALSE(blst_fp12_frobenius_map(&rejected, &a, 4));
    EXPECT_TRUE(BoundaryGt::equal(rejected, expected_mul));
}

TEST(CApi, MappingAndUncheckedWrappersMatchScalarBackend)
{
    BoundaryFp u = {};
    BoundaryFp v = {};
    BoundaryFp2 u2 = {};
    BoundaryFp2 v2 = {};
    BoundaryG1 mapped_g1 = {};
    BoundaryG2 mapped_g2 = {};
    BoundaryG1 expected_g1 = {};
    BoundaryG2 expected_g2 = {};
    BoundaryG1 checked_g1 = {};
    BoundaryG1 unchecked_g1 = {};
    BoundaryG2 checked_g2 = {};
    BoundaryG2 unchecked_g2 = {};
    BoundaryG1 normalized_g1 = {};
    BoundaryG2 normalized_g2 = {};
    const uint64_t three_words[6] = {3, 0, 0, 0, 0, 0};
    const uint64_t five_words[6] = {5, 0, 0, 0, 0, 0};
    const uint64_t seven_words[6] = {7, 0, 0, 0, 0, 0};
    const uint64_t eleven_words[6] = {11, 0, 0, 0, 0, 0};
    const uint64_t thirteen_words[6] = {13, 0, 0, 0, 0, 0};
    const uint64_t seventeen_words[6] = {17, 0, 0, 0, 0, 0};

    blst_fp_from_uint64(&u, three_words);
    blst_fp_from_uint64(&v, five_words);
    blst_fp_from_uint64(&u2[0], seven_words);
    blst_fp_from_uint64(&u2[1], eleven_words);
    blst_fp_from_uint64(&v2[0], thirteen_words);
    blst_fp_from_uint64(&v2[1], seventeen_words);

    blst_map_to_g1(&mapped_g1, &u, &v);
    blst_map_to_g2(&mapped_g2, &u2, &v2);
    blst::scalar::map_to_curve(expected_g1, u, &v);
    blst::scalar::map_to_curve(expected_g2, u2, &v2);

    const auto scalar = scalar_bytes(19);
    const BoundaryG1 point_g1 = p1_mul(23);
    const BoundaryG2 point_g2 = p2_mul(29);
    blst_p1_mult(&checked_g1, &point_g1, scalar.data(), 256);
    blst_p1_unchecked_mult(&unchecked_g1, &point_g1, scalar.data(), 256);
    blst_p2_mult(&checked_g2, &point_g2, scalar.data(), 256);
    blst_p2_unchecked_mult(&unchecked_g2, &point_g2, scalar.data(), 256);
    blst_p1_from_jacobian(&normalized_g1, &point_g1);
    blst_p2_from_jacobian(&normalized_g2, &point_g2);

    EXPECT_TRUE(BoundaryG1::equal(mapped_g1, expected_g1));
    EXPECT_TRUE(BoundaryG2::equal(mapped_g2, expected_g2));
    EXPECT_TRUE(BoundaryG1::equal(checked_g1, unchecked_g1));
    EXPECT_TRUE(BoundaryG2::equal(checked_g2, unchecked_g2));
    EXPECT_TRUE(BoundaryG1::equal(normalized_g1, point_g1));
    EXPECT_TRUE(BoundaryG2::equal(normalized_g2, point_g2));
    EXPECT_TRUE(BoundaryFp::equal(normalized_g1.z, BoundaryG1::generator().z));
    EXPECT_TRUE(fp2_equal(normalized_g2.z, BoundaryG2::generator().z));
}

TEST(CApi, ExpandMessageAndSizeWrappersWork)
{
    std::array<unsigned char, 64> out = {};
    std::array<unsigned char, 64> expected = {};
    std::array<unsigned char, 33> short_out = {};
    std::array<unsigned char, 33> short_expected = {};
    constexpr unsigned char msg[] = {'a', 'b', 'c'};
    constexpr unsigned char dst[] = {'B', 'L', 'S', '-', 'X', 'M', 'D'};

    blst_expand_message_xmd(out.data(), out.size(), msg, sizeof(msg), dst, sizeof(dst));
    blst::scalar::blst_expand_message_xmd(expected.data(), expected.size(), msg, sizeof(msg), dst, sizeof(dst));
    EXPECT_EQ(out, expected);
    blst_expand_message_xmd(short_out.data(), short_out.size(), msg, sizeof(msg), dst, sizeof(dst));
    blst::scalar::blst_expand_message_xmd(short_expected.data(), short_expected.size(), msg, sizeof(msg), dst, sizeof(dst));
    EXPECT_EQ(short_out, short_expected);

    EXPECT_EQ(blst_p1_sizeof(), sizeof(BoundaryG1));
    EXPECT_EQ(blst_p1_affine_sizeof(), sizeof(BoundaryG1Affine));
    EXPECT_EQ(blst_p2_sizeof(), sizeof(BoundaryG2));
    EXPECT_EQ(blst_p2_affine_sizeof(), sizeof(BoundaryG2Affine));
    EXPECT_EQ(blst_fp12_sizeof(), sizeof(BoundaryGt));
}

TEST(CApi, UniqTreeEnforcesCapacityAndReportsExactResults)
{
    constexpr size_t capacity = 2;
    constexpr unsigned char first[] = {'a'};
    constexpr unsigned char second[] = {'b'};
    constexpr unsigned char third[] = {'c'};
    const size_t bytes = blst_uniq_sizeof(capacity);
    const size_t words = (bytes + sizeof(std::uint64_t) - 1) / sizeof(std::uint64_t);
    std::vector<std::uint64_t> storage(words);
    auto* uniq = reinterpret_cast<blst_uniq*>(storage.data());

    ASSERT_NE(bytes, 0u);
    EXPECT_GE(alignof(std::uint64_t), blst_uniq_alignof());
    ASSERT_TRUE(blst_uniq_init(uniq, capacity));
    EXPECT_EQ(blst_uniq_insert(uniq, first, sizeof(first)), BLST_UNIQ_INSERTED);
    EXPECT_EQ(blst_uniq_insert(uniq, second, sizeof(second)), BLST_UNIQ_INSERTED);
    EXPECT_EQ(blst_uniq_insert(uniq, first, sizeof(first)), BLST_UNIQ_DUPLICATE);
    EXPECT_EQ(blst_uniq_insert(uniq, third, sizeof(third)),
              BLST_UNIQ_CAPACITY_EXCEEDED);
}

TEST(CApi, UniqTreeRejectsInvalidSizesAndInputs)
{
    unsigned char byte = 0;
    std::vector<std::uint64_t> storage(
        (blst_uniq_sizeof(1) + sizeof(std::uint64_t) - 1) /
        sizeof(std::uint64_t));
    auto* uniq = reinterpret_cast<blst_uniq*>(storage.data());

    EXPECT_EQ(blst_uniq_sizeof(0), 0u);
    EXPECT_EQ(blst_uniq_sizeof(std::numeric_limits<size_t>::max()), 0u);
    size_t low = 1;
    size_t high = std::numeric_limits<size_t>::max();
    while (low < high) {
        const size_t distance = high - low;
        const size_t middle = low + distance / 2 + distance % 2;
        if (blst_uniq_sizeof(middle) != 0)
            low = middle;
        else
            high = middle - 1;
    }
    ASSERT_LT(low, std::numeric_limits<size_t>::max());
    EXPECT_NE(blst_uniq_sizeof(low), 0u);
    EXPECT_EQ(blst_uniq_sizeof(low + 1), 0u);

    EXPECT_FALSE(blst_uniq_init(nullptr, 1));
    EXPECT_FALSE(blst_uniq_init(reinterpret_cast<blst_uniq*>(
        reinterpret_cast<unsigned char*>(uniq) + 1), 1));
    EXPECT_FALSE(blst_uniq_init(uniq, 0));
    ASSERT_TRUE(blst_uniq_init(uniq, 1));
    EXPECT_EQ(blst_uniq_insert(nullptr, &byte, 1), BLST_UNIQ_INVALID_INPUT);
    EXPECT_EQ(blst_uniq_insert(reinterpret_cast<blst_uniq*>(
                  reinterpret_cast<unsigned char*>(uniq) + 1), &byte, 1),
              BLST_UNIQ_INVALID_INPUT);
    EXPECT_EQ(blst_uniq_insert(uniq, nullptr, 1), BLST_UNIQ_INVALID_INPUT);
    EXPECT_EQ(blst_uniq_insert(
                  uniq, &byte, (std::numeric_limits<size_t>::max() >> 1) + 1),
              BLST_UNIQ_INVALID_INPUT);
    EXPECT_EQ(blst_uniq_insert(uniq, nullptr, 0), BLST_UNIQ_INSERTED);
    EXPECT_EQ(blst_uniq_insert(uniq, nullptr, 0), BLST_UNIQ_DUPLICATE);
}

TEST(CApi, UniqTreeHandlesLargeAdversarialInsertionOrder)
{
    constexpr size_t capacity = 32768;
    std::vector<std::array<unsigned char, 8>> messages(capacity);
    const size_t bytes = blst_uniq_sizeof(capacity);
    std::vector<std::uint64_t> storage(
        (bytes + sizeof(std::uint64_t) - 1) / sizeof(std::uint64_t));
    auto* uniq = reinterpret_cast<blst_uniq*>(storage.data());

    ASSERT_TRUE(blst_uniq_init(uniq, capacity));
    for (size_t i = 0; i < capacity; i++) {
        size_t value = i;
        for (size_t j = 0; j < messages[i].size(); j++) {
            messages[i][messages[i].size() - 1 - j] =
                static_cast<unsigned char>(value);
            value >>= 8;
        }
        ASSERT_EQ(blst_uniq_insert(uniq, messages[i].data(), messages[i].size()),
                  BLST_UNIQ_INSERTED);
    }

    for (size_t i = 0; i < capacity; i += 257)
        EXPECT_EQ(blst_uniq_insert(uniq, messages[i].data(), messages[i].size()),
                  BLST_UNIQ_DUPLICATE);
}

TEST(CApi, ExpandMessageOversizeOutputClearsBuffer)
{
    constexpr unsigned char msg[] = {'a', 'b', 'c'};
    constexpr unsigned char dst[] = {'B', 'L', 'S', '-', 'X', 'M', 'D'};
    std::vector<unsigned char> out(255 * 32 + 1, 0xA5);

    blst_expand_message_xmd(out.data(), out.size(), msg, sizeof(msg), dst, sizeof(dst));
    EXPECT_TRUE(std::all_of(out.begin(), out.end(), [](unsigned char byte) { return byte == 0; }));
}

TEST(CApi, ScalarByteDecodersRejectEmptyInputAndClearOutput)
{
    blst_scalar out = {};
    std::fill_n(out.b, sizeof(out.b), static_cast<blst_byte>(0xA5));
    EXPECT_FALSE(blst_scalar_from_le_bytes(&out, nullptr, 0));
    EXPECT_TRUE(scalar_is_zero(out));

    std::fill_n(out.b, sizeof(out.b), static_cast<blst_byte>(0x5A));
    EXPECT_FALSE(blst_scalar_from_be_bytes(&out, nullptr, 0));
    EXPECT_TRUE(scalar_is_zero(out));
}

TEST(CApi, KeygenInvalidInputsZeroizeOutput)
{
    blst_scalar out = {};
    std::fill_n(out.b, sizeof(out.b), static_cast<blst_byte>(0xA5));

    const std::array<unsigned char, 31> short_ikm = {};
    blst_keygen(&out, short_ikm.data(), short_ikm.size(), nullptr, 0);
    EXPECT_TRUE(scalar_is_zero(out));

    std::fill_n(out.b, sizeof(out.b), static_cast<blst_byte>(0x5A));
    const auto valid_ikm = scalar_bytes(9);
    blst_keygen_v5(&out, valid_ikm.data(), valid_ikm.size(), nullptr, 0, nullptr, 0);
    EXPECT_TRUE(scalar_is_zero(out));
}

TEST(CApi, KeygenV5AcceptsNonNullEmptySalt)
{
    constexpr unsigned char expected[32] = {
        0x10, 0xc1, 0x86, 0xc2, 0xce, 0xf3, 0x66, 0x40,
        0xd9, 0x97, 0xe4, 0xec, 0x01, 0x4f, 0xe5, 0xbe,
        0xc9, 0xa1, 0x91, 0x5f, 0x8a, 0x1b, 0x2c, 0x43,
        0x27, 0x4c, 0x17, 0x3f, 0xca, 0x3d, 0xba, 0x19
    };
    const auto valid_ikm = scalar_bytes(9);
    const unsigned char non_null_empty_salt = 0;
    blst_scalar first = {};
    blst_scalar second = {};

    blst_keygen_v5(
        &first, valid_ikm.data(), valid_ikm.size(),
        &non_null_empty_salt, 0, nullptr, 0);
    blst_keygen_v5(
        &second, valid_ikm.data(), valid_ikm.size(),
        &non_null_empty_salt, 0, nullptr, 0);

    EXPECT_FALSE(scalar_is_zero(first));
    EXPECT_EQ(0, std::memcmp(first.b, expected, sizeof(expected)));
    EXPECT_EQ(0, std::memcmp(first.b, second.b, sizeof(first.b)));
}

TEST(CApi, Eip2333PublishedVectors)
{
    struct Vector {
        const char* seed;
        const char* master;
        unsigned int child_index;
        const char* child;
    };
    static const Vector vectors[] = {
        {
            "c55257c360c07c72029aebc1b53c05ed0362ada38ead3e3e9efa3708e5349553"
            "1f09a6987599d18264c1e1c92f2cf141630c7a3c4ab7c81b2f001698e7463b04",
            "0d7359d57963ab8fbbde1852dcf553fedbc31f464d80ee7d40ae683122b45070",
            0,
            "2d18bd6c14e6d15bf8b5085c9b74f3daae3b03cc2014770a599d8c1539e50f8e",
        },
        {
            "3141592653589793238462643383279502884197169399375105820974944592",
            "41c9e07822b092a93fd6797396338c3ada4170cc81829fdfce6b5d34bd5e7ec7",
            3141592653u,
            "384843fad5f3d777ea39de3e47a8f999ae91f89e42bffa993d91d9782d152a0f",
        },
        {
            "0099ff991111002299dd7744ee3355bbdd8844115566cc55663355668888cc00",
            "3cfa341ab3910a7d00d933d8f7c4fe87c91798a0397421d6b19fd5b815132e80",
            0xffffffffu,
            "40e86285582f35b28821340f6a53b448588efa575bc4d88c32ef8567b8d9479b",
        },
        {
            "d4e56740f876aef8c010b86a40d5f56745a118d0906a34e69aec8c0db1cb8fa3",
            "2a0e28ffa5fbbe2f8e7aad4ed94f745d6bf755c51182e119bb1694fe61d3afca",
            42,
            "455c0dc9fccb3395825d92a60d2672d69416be1c2578a87a7a3d3ced11ebb88d",
        },
    };

    const auto decode = [](const char* hex) {
        const size_t chars = std::strlen(hex);
        EXPECT_EQ(chars % 2, 0u);
        std::vector<unsigned char> out(chars / 2);
        const auto nibble = [](char c) -> unsigned char {
            if (c >= '0' && c <= '9') return static_cast<unsigned char>(c - '0');
            if (c >= 'a' && c <= 'f') return static_cast<unsigned char>(c - 'a' + 10);
            return 0xff;
        };
        for (size_t i = 0; i < out.size(); ++i) {
            const unsigned char hi = nibble(hex[2 * i]);
            const unsigned char lo = nibble(hex[2 * i + 1]);
            EXPECT_LT(hi, 16u);
            EXPECT_LT(lo, 16u);
            out[i] = static_cast<unsigned char>((hi << 4) | lo);
        }
        return out;
    };

    for (const Vector& vector : vectors) {
        const std::vector<unsigned char> seed = decode(vector.seed);
        const std::vector<unsigned char> expected_master = decode(vector.master);
        const std::vector<unsigned char> expected_child = decode(vector.child);
        blst_scalar master = {}, child = {};
        std::array<unsigned char, 32> encoded = {};

        blst_derive_master_eip2333(&master, seed.data(), seed.size());
        blst_bendian_from_scalar(encoded.data(), &master);
        EXPECT_TRUE(std::equal(encoded.begin(), encoded.end(), expected_master.begin()));

        blst_derive_child_eip2333(&child, &master, vector.child_index);
        blst_bendian_from_scalar(encoded.data(), &child);
        EXPECT_TRUE(std::equal(encoded.begin(), encoded.end(), expected_child.begin()));
    }
}

TEST(CApi, WbitsApisRejectOutOfRangeValues)
{
    EXPECT_EQ(blst_p1s_mult_wbits_precompute_sizeof(0, 4), 0u);
    EXPECT_EQ(blst_p1s_mult_wbits_precompute_sizeof(15, 4), 0u);
    EXPECT_EQ(blst_p1s_mult_wbits_precompute_sizeof(14, std::numeric_limits<size_t>::max()), 0u);

    BoundaryG1Affine affine = {};
    const BoundaryG1Affine* affine_ptrs[1] = {&affine};
    blst_p1s_mult_wbits_precompute(nullptr, 4, affine_ptrs, 1);

    BoundaryG1 out = p1_mul(3);
    const auto scalar = scalar_bytes(7);
    const unsigned char* scalar_ptrs[2] = {scalar.data(), nullptr};
    blst_p1s_mult_wbits(&out, nullptr, 0, 1, scalar_ptrs, 256, nullptr);
    EXPECT_TRUE(blst_p1_is_inf(&out));
}

TEST(CApi, MsmBoundariesRejectInvalidWidthsAndHandleZeroBits)
{
    EXPECT_EQ(blst_p1s_mult_pippenger_scratch_sizeof(0), 0u);
    EXPECT_EQ(blst_p2s_mult_pippenger_scratch_sizeof(0), 0u);
    EXPECT_NE(blst_p1s_mult_pippenger_window_scratch_sizeof(24), 0u);
    EXPECT_NE(blst_p2s_mult_pippenger_window_scratch_sizeof(24), 0u);
    EXPECT_EQ(blst_p1s_mult_pippenger_window_scratch_sizeof(0), 0u);
    EXPECT_EQ(blst_p1s_mult_pippenger_window_scratch_sizeof(25), 0u);
    EXPECT_EQ(blst_p2s_mult_pippenger_window_scratch_sizeof(25), 0u);
    EXPECT_EQ(
        blst_p1s_mult_pippenger_window_scratch_sizeof(
            std::numeric_limits<size_t>::max()),
        0u);

    BoundaryG1 g1 = p1_mul(3);
    blst_p1s_mult_pippenger(&g1, nullptr, 2, nullptr, 0, nullptr);
    EXPECT_TRUE(blst_p1_is_inf(&g1));

    BoundaryG2 g2 = p2_mul(3);
    blst_p2s_mult_pippenger(&g2, nullptr, 2, nullptr, 0, nullptr);
    EXPECT_TRUE(blst_p2_is_inf(&g2));

    g1 = p1_mul(3);
    blst_p1s_mult_wbits(&g1, nullptr, 4, 2, nullptr, 0, nullptr);
    EXPECT_TRUE(blst_p1_is_inf(&g1));
}

TEST(CApi, BatchAffineInternalRejectsOversizedScalarWidths)
{
    constexpr size_t window = 1;
    const size_t scratch_bytes =
        blst::scalar::msm_mult_batch_affine_scratch_sizeof<blst::Fp, 128>(
            window);
    ASSERT_NE(scratch_bytes, 0);
    ASSERT_EQ(scratch_bytes % sizeof(blst::limb_t), 0);
    std::vector<blst::limb_t> scratch(
        scratch_bytes / sizeof(blst::limb_t));
    const blst::scalar::G1AffinePoint point =
        *blst_p1_affine_generator();
    std::array<blst::byte, 33> scalar = {};
    scalar[0] = 1;

    blst::scalar::G1Point actual = *blst_p1_generator();
    const blst::scalar::G1Point unchanged = actual;
    EXPECT_FALSE((blst::scalar::msm_mult_batch_affine<blst::Fp, 128>(
        actual, &point, scalar.data(), 1, 257,
        scratch.data(), scratch_bytes, window)));
    EXPECT_TRUE(blst_p1_is_equal(&actual, &unchanged));

    EXPECT_FALSE((blst::scalar::msm_mult_batch_affine<blst::Fp, 128>(
        actual, &point, scalar.data(), 1, 0,
        scratch.data(), scratch_bytes, window)));
    EXPECT_TRUE(blst_p1_is_equal(&actual, &unchanged));
    EXPECT_EQ(
        (blst::scalar::msm_mult_batch_affine_scratch_sizeof<blst::Fp, 128>(
            25)),
        0);
}

TEST(CApi, LeadingNullPointArraysAreRejectedWithoutPointerArithmetic)
{
    BoundaryG1 sum = p1_mul(3);
    const BoundaryG1Affine* points[1] = {nullptr};
    blst_p1s_add(&sum, points, 1);
    EXPECT_TRUE(blst_p1_is_inf(&sum));

    BoundaryG1Affine affine_out[1];
    affine_out[0] = {};
    const BoundaryG1* proj_points[1] = {nullptr};
    blst_p1s_to_affine(affine_out, proj_points, 1);
    EXPECT_TRUE(blst_p1_affine_is_inf(&affine_out[0]));
}
