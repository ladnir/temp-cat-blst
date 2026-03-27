#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fr.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/point.h"
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/map_to_G.h"
#endif

extern "C" {
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
void blst_p1_to_affine(blst::scalar::G1AffinePoint* out, const blst::scalar::G1Point* in);
void blst_p1_from_affine(blst::scalar::G1Point* out, const blst::scalar::G1AffinePoint* in);
bool blst_p1_on_curve(const blst::scalar::G1Point* p);
bool blst_p1_in_g1(const blst::scalar::G1Point* p);
bool blst_p1_is_equal(const blst::scalar::G1Point* a, const blst::scalar::G1Point* b);
bool blst_p1_is_inf(const blst::scalar::G1Point* a);
const blst::scalar::G1Point* blst_p1_generator(void);
bool blst_p1_affine_in_g1(const blst::scalar::G1AffinePoint* p);
bool blst_p1_affine_is_equal(const blst::scalar::G1AffinePoint* a, const blst::scalar::G1AffinePoint* b);
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
void blst_p2_to_affine(blst::scalar::G2AffinePoint* out, const blst::scalar::G2Point* in);
void blst_p2_from_affine(blst::scalar::G2Point* out, const blst::scalar::G2AffinePoint* in);
bool blst_p2_on_curve(const blst::scalar::G2Point* p);
bool blst_p2_in_g2(const blst::scalar::G2Point* p);
bool blst_p2_is_equal(const blst::scalar::G2Point* a, const blst::scalar::G2Point* b);
bool blst_p2_is_inf(const blst::scalar::G2Point* a);
const blst::scalar::G2Point* blst_p2_generator(void);
bool blst_p2_affine_in_g2(const blst::scalar::G2AffinePoint* p);
bool blst_p2_affine_is_equal(const blst::scalar::G2AffinePoint* a, const blst::scalar::G2AffinePoint* b);
void blst_p2s_to_affine(blst::scalar::G2AffinePoint dst[], const blst::scalar::G2Point* const points[], size_t npoints);
void blst_p2s_add(blst::scalar::G2Point* ret, const blst::scalar::G2AffinePoint* const points[], size_t npoints);
size_t blst_p2s_mult_pippenger_scratch_sizeof(size_t npoints);
void blst_p2s_mult_pippenger(blst::scalar::G2Point* ret, const blst::scalar::G2AffinePoint* const points[], size_t npoints, const unsigned char* const scalars[], size_t nbits, blst::limb_t* scratch);
void blst_p2_affine_serialize(unsigned char out[192], const blst::scalar::G2AffinePoint* in);
void blst_p2_affine_compress(unsigned char out[96], const blst::scalar::G2AffinePoint* in);
blst::BLST_ERROR blst_p2_uncompress(blst::scalar::G2AffinePoint* out, const unsigned char in[96]);
blst::BLST_ERROR blst_p2_deserialize(blst::scalar::G2AffinePoint* out, const unsigned char in[192]);

void blst_miller_loop(blst::Fp12* ret, const blst::scalar::G2AffinePoint* Q, const blst::scalar::G1AffinePoint* P);
void blst_miller_loop_n(blst::Fp12* ret, const blst::scalar::G2AffinePoint* const Qs[], const blst::scalar::G1AffinePoint* const Ps[], size_t n);
void blst_precompute_lines(blst::Fp6 Qlines[68], const blst::scalar::G2AffinePoint* Q);
void blst_miller_loop_lines(blst::Fp12* ret, const blst::Fp6 Qlines[68], const blst::scalar::G1AffinePoint* P);
void blst_hash_to_g2(blst::scalar::G2Point* out, const unsigned char* msg, size_t msg_len, const unsigned char* DST, size_t DST_len, const unsigned char* aug, size_t aug_len);
size_t blst_pairing_sizeof(void);
void blst_pairing_init(void* new_ctx, bool hash_or_encode, const unsigned char* DST, size_t DST_len);
blst::BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g1(void* ctx, const blst::scalar::G1AffinePoint* PK, bool pk_grpchk, const blst::scalar::G2AffinePoint* signature, bool sig_grpchk, const unsigned char* msg, size_t msg_len, const unsigned char* aug, size_t aug_len);
void blst_pairing_commit(void* ctx);
bool blst_pairing_finalverify(const void* ctx, const blst::Fp12* gtsig);
void blst_aggregated_in_g2(blst::Fp12* out, const blst::scalar::G2AffinePoint* signature);
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

} // namespace

TEST(LegacyCApi, FrPrimitivesRoundTripAndInverse)
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

TEST(LegacyCApi, FpPrimitivesRoundTripInverseAndSqrt)
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

TEST(LegacyCApi, Fp2PrimitivesInverseAndSqrt)
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

TEST(LegacyCApi, G1ProjectiveAndAffineRoundTrip)
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

TEST(LegacyCApi, G2ProjectiveAndAffineRoundTrip)
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

TEST(LegacyCApi, G1MultipointAffineAndAddMatchSinglePointOps)
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

TEST(LegacyCApi, G1MultipointPippengerMatchesSinglePointOps)
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

TEST(LegacyCApi, G1MultipointWbitsMatchesSinglePointOps)
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

TEST(LegacyCApi, G1MultipointTileMatchesSinglePointOps)
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

TEST(LegacyCApi, G2MultipointWrappersMatchSinglePointOps)
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

TEST(LegacyCApi, PairingPrimitivesMatchBackend)
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

TEST(LegacyCApi, MillerLoopNMatchesNaiveProduct)
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

TEST(LegacyCApi, PairingAggregateAllowsNullSignatureAfterFirstElement)
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
