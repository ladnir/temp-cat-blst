#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fr.h"
#include "blst/scalar/fields/hash_to_field.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/map_to_G.h"
#include "blst/scalar/curve/point.h"
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/map_to_G.h"
#endif

struct alignas(64) BlstAvxP1
{
    unsigned char b[384];
};

struct alignas(64) BlstAvxP1Affine
{
    unsigned char b[256];
};

struct alignas(64) BlstAvxP2
{
    unsigned char b[768];
};

struct alignas(64) BlstAvxP2Affine
{
    unsigned char b[512];
};

struct alignas(64) BlstAvxFp12
{
    unsigned char b[1536];
};

extern "C" {
bool blst_scalar_from_le_bytes(blst::pow256 out, const unsigned char* in, size_t len);
bool blst_scalar_from_be_bytes(blst::pow256 out, const unsigned char* in, size_t len);
bool blst_sk_add_n_check(blst::pow256 out, const blst::pow256 a, const blst::pow256 b);
bool blst_sk_sub_n_check(blst::pow256 out, const blst::pow256 a, const blst::pow256 b);
bool blst_sk_mul_n_check(blst::pow256 out, const blst::pow256 a, const blst::pow256 b);
void blst_sk_inverse(blst::pow256 out, const blst::pow256 a);

void blst_fr_add(blst::Fr* ret, const blst::Fr* a, const blst::Fr* b);
void blst_fr_mul(blst::Fr* ret, const blst::Fr* a, const blst::Fr* b);
void blst_fr_sqr(blst::Fr* ret, const blst::Fr* a);
void blst_fr_eucl_inverse(blst::Fr* ret, const blst::Fr* a);
void blst_fr_inverse(blst::Fr* ret, const blst::Fr* a);
void blst_fr_from_uint64(blst::Fr* ret, const unsigned long long a[4]);
void blst_uint64_from_fr(unsigned long long ret[4], const blst::Fr* a);
void blst_fr_from_scalar(blst::Fr* ret, const blst::pow256 a);
void blst_scalar_from_fr(blst::pow256 ret, const blst::Fr* a);
void blst_fr_from_hexascii(blst::Fr* ret, const char* hex);

void blst_fp_add(blst::Fp* ret, const blst::Fp* a, const blst::Fp* b);
void blst_fp_mul(blst::Fp* ret, const blst::Fp* a, const blst::Fp* b);
void blst_fp_sqr(blst::Fp* ret, const blst::Fp* a);
void blst_fp_eucl_inverse(blst::Fp* ret, const blst::Fp* a);
void blst_fp_inverse(blst::Fp* ret, const blst::Fp* a);
int blst_fp_sqrt(blst::Fp* ret, const blst::Fp* a);
void blst_fp_from_uint32(blst::Fp* ret, const unsigned int a[12]);
void blst_uint32_from_fp(unsigned int ret[12], const blst::Fp* a);
void blst_fp_from_uint64(blst::Fp* ret, const unsigned long long a[6]);
void blst_uint64_from_fp(unsigned long long ret[6], const blst::Fp* a);
void blst_fp_from_bendian(blst::Fp* ret, const unsigned char a[48]);
void blst_bendian_from_fp(unsigned char ret[48], const blst::Fp* a);
int blst_fp_is_square(const blst::Fp* a);
void blst_fp_from_hexascii(blst::Fp* ret, const char* hex);

void blst_fp2_mul(blst::Fp2* ret, const blst::Fp2* a, const blst::Fp2* b);
void blst_fp2_sqr(blst::Fp2* ret, const blst::Fp2* a);
void blst_fp2_eucl_inverse(blst::Fp2* ret, const blst::Fp2* a);
void blst_fp2_inverse(blst::Fp2* ret, const blst::Fp2* a);
int blst_fp2_sqrt(blst::Fp2* ret, const blst::Fp2* a);
int blst_fp2_is_square(const blst::Fp2* a);

void blst_fp12_mul(blst::Fp12* ret, const blst::Fp12* a, const blst::Fp12* b);
void blst_fp12_sqr(blst::Fp12* ret, const blst::Fp12* a);
void blst_fp12_cyclotomic_sqr(blst::Fp12* ret, const blst::Fp12* a);
void blst_fp12_mul_by_xy00z0(blst::Fp12* ret, const blst::Fp12* a, const blst::Fp6* xy00z0);
void blst_fp12_conjugate(blst::Fp12* a);
void blst_fp12_inverse(blst::Fp12* ret, const blst::Fp12* a);
void blst_fp12_frobenius_map(blst::Fp12* ret, const blst::Fp12* a, size_t n);
bool blst_fp12_is_equal(const blst::Fp12* a, const blst::Fp12* b);
bool blst_fp12_is_one(const blst::Fp12* a);
bool blst_fp12_in_group(const blst::Fp12* a);
const blst::Fp12* blst_fp12_one(void);
void blst_final_exp(blst::Fp12* ret, const blst::Fp12* f);
bool blst_fp12_finalverify(const blst::Fp12* gt1, const blst::Fp12* gt2);

void blst_p1_add(blst::scalar::G1Point* out, const blst::scalar::G1Point* a, const blst::scalar::G1Point* b);
void blst_p1_add_or_double(blst::scalar::G1Point* out, const blst::scalar::G1Point* a, const blst::scalar::G1Point* b);
void blst_p1_double(blst::scalar::G1Point* out, const blst::scalar::G1Point* a);
void blst_p1_mult(blst::scalar::G1Point* out, const blst::scalar::G1Point* p, const unsigned char* scalar, size_t nbits);
void blst_p1_unchecked_mult(blst::scalar::G1Point* out, const blst::scalar::G1Point* p, const unsigned char* scalar, size_t nbits);
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
void blst_p2_add_or_double(blst::scalar::G2Point* out, const blst::scalar::G2Point* a, const blst::scalar::G2Point* b);
void blst_p2_double(blst::scalar::G2Point* out, const blst::scalar::G2Point* a);
void blst_p2_mult(blst::scalar::G2Point* out, const blst::scalar::G2Point* p, const unsigned char* scalar, size_t nbits);
void blst_p2_unchecked_mult(blst::scalar::G2Point* out, const blst::scalar::G2Point* p, const unsigned char* scalar, size_t nbits);
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
void blst_pairing_init(void* new_ctx, bool hash_or_encode, const unsigned char* DST, size_t DST_len);
blst::BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g1(void* ctx, const blst::scalar::G1AffinePoint* PK, bool pk_grpchk, const blst::scalar::G2AffinePoint* signature, bool sig_grpchk, const unsigned char* msg, size_t msg_len, const unsigned char* aug, size_t aug_len);
void blst_pairing_commit(void* ctx);
bool blst_pairing_finalverify(const void* ctx, const blst::Fp12* gtsig);
void blst_aggregated_in_g2(blst::Fp12* out, const blst::scalar::G2AffinePoint* signature);
size_t blst_uniq_sizeof(size_t n_nodes);
void blst_uniq_init(void* tree);
bool blst_uniq_test(void* tree, const unsigned char* msg, size_t len);
void blst_keygen(blst::pow256 out_SK, const unsigned char* IKM, size_t IKM_len, const unsigned char* info, size_t info_len);
void blst_keygen_v5(blst::pow256 out_SK, const unsigned char* IKM, size_t IKM_len, const unsigned char* salt, size_t salt_len, const unsigned char* info, size_t info_len);
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

std::array<unsigned char, 32> scalar_bytes(std::uint64_t value)
{
    std::array<unsigned char, 32> out = {};
    for (size_t i = 0; i < sizeof(value); ++i) {
        out[i] = static_cast<unsigned char>(value >> (8 * i));
    }
    return out;
}

bool scalar_is_zero(const blst::pow256& value)
{
    for (unsigned char byte : value) {
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
    blst::scalar::GT::deserialize(scalar_gt, serialized);
    return scalar_gt.data;
}

PublicG1Affine to_public(const BoundaryG1Affine& point)
{
#if BLST_DEFAULT_AVX
    return PublicG1Affine(
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.X),
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.Y));
#else
    return point;
#endif
}

PublicG2Affine to_public(const BoundaryG2Affine& point)
{
#if BLST_DEFAULT_AVX
    return PublicG2Affine(
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.X),
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.Y));
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

} // namespace

TEST(CApi, FrPrimitivesRoundTripAndInverse)
{
    const unsigned long long five_words[4] = {5, 0, 0, 0};
    const unsigned long long one_words[4] = {1, 0, 0, 0};
    BoundaryFr five = {};
    BoundaryFr one = {};
    BoundaryFr five_from_hex = {};
    blst::pow256 scalar = {};
    BoundaryFr from_scalar = {};

    blst_fr_from_uint64(&five, five_words);
    blst_fr_from_uint64(&one, one_words);
    blst_fr_from_hexascii(&five_from_hex, "0000000000000000000000000000000000000000000000000000000000000005");

    EXPECT_TRUE(BoundaryFr::equal(five, five_from_hex));

    blst_scalar_from_fr(scalar, &five);
    blst_fr_from_scalar(&from_scalar, scalar);
    EXPECT_TRUE(BoundaryFr::equal(five, from_scalar));

    unsigned long long roundtrip[4] = {};
    blst_uint64_from_fr(roundtrip, &five);
    EXPECT_EQ(roundtrip[0], 5ULL);
    EXPECT_EQ(roundtrip[1], 0ULL);

    BoundaryFr inverse_eucl = {};
    BoundaryFr inverse = {};
    BoundaryFr product = {};
    blst_fr_eucl_inverse(&inverse_eucl, &five);
    blst_fr_inverse(&inverse, &five);
    blst_fr_mul(&product, &five, &inverse_eucl);
    EXPECT_TRUE(BoundaryFr::equal(product, one));
    EXPECT_TRUE(BoundaryFr::equal(inverse_eucl, inverse));
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
    const unsigned long long three_words[6] = {3, 0, 0, 0, 0, 0};
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
    const std::vector<unsigned char> scalars =
        scalar_batch({msm_scalars[0], msm_scalars[1], msm_scalars[2], msm_scalars[3]});
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
    const std::vector<unsigned char> scalars =
        scalar_batch({msm_scalars[0], msm_scalars[1], msm_scalars[2], msm_scalars[3]});
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
    BoundaryGt roundtrip = {};
    BoundaryGt scalar_product = {};
    ASSERT_TRUE(blst_avx_fp12_from_scalar(&avx_miller, &miller));
    ASSERT_TRUE(blst_avx_fp12_final_exp(&avx_expected, &avx_miller));
    ASSERT_TRUE(blst_fp12_from_avx(&roundtrip, &avx_expected));
    EXPECT_TRUE(BoundaryGt::equal(roundtrip, expected));
    EXPECT_TRUE(blst_avx_fp12_in_group(&avx_expected));

    const auto* avx_one = blst_avx_fp12_one();
    ASSERT_NE(avx_one, nullptr);
    ASSERT_TRUE(blst_avx_fp12_mul(&avx_product, &avx_expected, avx_one));
    ASSERT_TRUE(blst_fp12_from_avx(&scalar_product, &avx_product));
    EXPECT_TRUE(BoundaryGt::equal(scalar_product, expected));
    EXPECT_TRUE(blst_avx_fp12_is_equal(&avx_product, &avx_expected));
}

TEST(CApi, FpPrimitivesRoundTripInverseAndSqrt)
{
    const unsigned long long three_words[6] = {3, 0, 0, 0, 0, 0};
    const unsigned long long one_words[6] = {1, 0, 0, 0, 0, 0};
    BoundaryFp three = {};
    BoundaryFp one = {};
    BoundaryFp from_hex = {};
    blst_fp_from_uint64(&three, three_words);
    blst_fp_from_uint64(&one, one_words);
    blst_fp_from_hexascii(&from_hex, "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000003");
    EXPECT_TRUE(BoundaryFp::equal(three, from_hex));

    unsigned char bendian[48] = {};
    BoundaryFp roundtrip = {};
    blst_bendian_from_fp(bendian, &three);
    blst_fp_from_bendian(&roundtrip, bendian);
    EXPECT_TRUE(BoundaryFp::equal(three, roundtrip));

    unsigned long long roundtrip_words[6] = {};
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
}

TEST(CApi, Fp2PrimitivesInverseAndSqrt)
{
    const unsigned long long two_words[6] = {2, 0, 0, 0, 0, 0};
    const unsigned long long one_words[6] = {1, 0, 0, 0, 0, 0};
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
        const size_t bucket_words = blst_p1s_mult_pippenger_scratch_sizeof(0) / sizeof(blst::limb_t);
        std::vector<blst::limb_t> scratch(bucket_words << (window - 1));

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

    PublicGt expected_public = {};
    const PublicG2Affine public_q = to_public(q_affine);
    const PublicG1Affine public_p = to_public(p_affine);
    PublicGt::miller_loop(expected_public, &public_q, &public_p);
    PublicGt::final_mult(expected_public, expected_public);
    const BoundaryGt expected = gt_to_boundary(expected_public);

    EXPECT_TRUE(blst::Fp12::equal(final_exp, expected));
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

TEST(CApi, MillerLoopNMatchesNaiveProductBeyondScalarBatchLimit)
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

    BoundaryGt combined = {};
    blst_miller_loop_n(&combined, qs_affine_ptrs, ps_affine_ptrs, ps.size());

    BoundaryGt naive = *blst_fp12_one();
    for (size_t i = 0; i < npoints; ++i) {
        BoundaryGt term = {};
        blst_miller_loop(&term, &qs_affine[i], &ps_affine[i]);
        blst_fp12_mul(&naive, &naive, &term);
    }

    EXPECT_TRUE(blst_fp12_is_equal(&combined, &naive));
}

TEST(CApi, MillerLoopNEmptyAndMalformedLeadingPointersReturnIdentity)
{
    BoundaryGt out = {};
    blst_miller_loop_n(&out, nullptr, nullptr, 0);
    EXPECT_TRUE(blst_fp12_is_one(&out));

    const BoundaryG1Affine* ps[1] = {nullptr};
    const BoundaryG2Affine* qs[1] = {nullptr};
    blst_miller_loop_n(&out, qs, ps, 1);
    EXPECT_TRUE(blst_fp12_is_one(&out));
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

    std::vector<unsigned long long> raw(blst_pairing_sizeof() / sizeof(unsigned long long));
    void* ctx = raw.data();
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
    const BoundaryGt a = blst::scalar::GT::generator().data;
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
    blst_fp12_frobenius_map(&frob, &a, 2);
    blst_fp12_conjugate(&conjugated);

    BoundaryGt expected_mul = {};
    BoundaryGt expected_inverse = {};
    BoundaryGt expected_frob = {};
    BoundaryGt expected_conjugated = a;
    BoundaryFp6 line = {};
    BoundaryGt xy_mul = {};
    BoundaryGt expected_xy_mul = {};

    const unsigned long long five_words[6] = {5, 0, 0, 0, 0, 0};
    const unsigned long long seven_words[6] = {7, 0, 0, 0, 0, 0};
    const unsigned long long eleven_words[6] = {11, 0, 0, 0, 0, 0};
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
    const unsigned long long three_words[6] = {3, 0, 0, 0, 0, 0};
    const unsigned long long five_words[6] = {5, 0, 0, 0, 0, 0};
    const unsigned long long seven_words[6] = {7, 0, 0, 0, 0, 0};
    const unsigned long long eleven_words[6] = {11, 0, 0, 0, 0, 0};
    const unsigned long long thirteen_words[6] = {13, 0, 0, 0, 0, 0};
    const unsigned long long seventeen_words[6] = {17, 0, 0, 0, 0, 0};

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
    EXPECT_TRUE(BoundaryFp::equal(normalized_g1.Z, BoundaryG1::generator().Z));
    EXPECT_TRUE(fp2_equal(normalized_g2.Z, BoundaryG2::generator().Z));
}

TEST(CApi, ExpandMessageUniqAndSizeWrappersWork)
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

    EXPECT_EQ(blst_uniq_sizeof(0), 0u);
    std::vector<std::uint64_t> uniq_bytes(blst_uniq_sizeof(4) / sizeof(std::uint64_t));
    void* uniq = uniq_bytes.data();
    blst_uniq_init(uniq);
    EXPECT_TRUE(blst_uniq_test(uniq, msg, sizeof(msg)));
    EXPECT_FALSE(blst_uniq_test(uniq, msg, sizeof(msg)));

    EXPECT_EQ(blst_p1_sizeof(), sizeof(BoundaryG1));
    EXPECT_EQ(blst_p1_affine_sizeof(), sizeof(BoundaryG1Affine));
    EXPECT_EQ(blst_p2_sizeof(), sizeof(BoundaryG2));
    EXPECT_EQ(blst_p2_affine_sizeof(), sizeof(BoundaryG2Affine));
    EXPECT_EQ(blst_fp12_sizeof(), sizeof(BoundaryGt));
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
    blst::pow256 out = {};
    std::fill_n(out, sizeof(out), static_cast<unsigned char>(0xA5));
    EXPECT_FALSE(blst_scalar_from_le_bytes(out, nullptr, 0));
    EXPECT_TRUE(scalar_is_zero(out));

    std::fill_n(out, sizeof(out), static_cast<unsigned char>(0x5A));
    EXPECT_FALSE(blst_scalar_from_be_bytes(out, nullptr, 0));
    EXPECT_TRUE(scalar_is_zero(out));
}

TEST(CApi, KeygenInvalidInputsZeroizeOutput)
{
    blst::pow256 out = {};
    std::fill_n(out, sizeof(out), static_cast<unsigned char>(0xA5));

    const std::array<unsigned char, 31> short_ikm = {};
    blst_keygen(out, short_ikm.data(), short_ikm.size(), nullptr, 0);
    EXPECT_TRUE(scalar_is_zero(out));

    std::fill_n(out, sizeof(out), static_cast<unsigned char>(0x5A));
    const auto valid_ikm = scalar_bytes(9);
    blst_keygen_v5(out, valid_ikm.data(), valid_ikm.size(), nullptr, 0, nullptr, 0);
    EXPECT_TRUE(scalar_is_zero(out));
}

TEST(CApi, WbitsApisRejectOutOfRangeValues)
{
    EXPECT_EQ(blst_p1s_mult_wbits_precompute_sizeof(0, 4), 0u);
    EXPECT_EQ(blst_p1s_mult_wbits_precompute_sizeof(15, 4), 0u);
    EXPECT_EQ(blst_p1s_mult_wbits_precompute_sizeof(14, std::numeric_limits<size_t>::max()), 0u);

    BoundaryG1 out = p1_mul(3);
    const auto scalar = scalar_bytes(7);
    const unsigned char* scalar_ptrs[2] = {scalar.data(), nullptr};
    blst_p1s_mult_wbits(&out, nullptr, 0, 1, scalar_ptrs, 256, nullptr);
    EXPECT_TRUE(blst_p1_is_inf(&out));
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
