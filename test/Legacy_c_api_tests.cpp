#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <type_traits>

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
void blst_p2_affine_serialize(unsigned char out[192], const blst::scalar::G2AffinePoint* in);
void blst_p2_affine_compress(unsigned char out[96], const blst::scalar::G2AffinePoint* in);
blst::BLST_ERROR blst_p2_uncompress(blst::scalar::G2AffinePoint* out, const unsigned char in[96]);
blst::BLST_ERROR blst_p2_deserialize(blst::scalar::G2AffinePoint* out, const unsigned char in[192]);

void blst_miller_loop(blst::Fp12* ret, const blst::scalar::G2AffinePoint* Q, const blst::scalar::G1AffinePoint* P);
void blst_precompute_lines(blst::Fp6 Qlines[68], const blst::scalar::G2AffinePoint* Q);
void blst_miller_loop_lines(blst::Fp12* ret, const blst::Fp6 Qlines[68], const blst::scalar::G1AffinePoint* P);
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
