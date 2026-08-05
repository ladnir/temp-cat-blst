/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <gtest/gtest.h>

#include "bindings/blst.h"
#include "blst/curve/map_to_G.h"
#include "blst/scalar/curve/map_to_G.h"

namespace {

using blst::scalar::Fp;
using blst::scalar::Fp2;
using blst::scalar::G1AffinePoint;
using blst::scalar::G1Point;
using blst::scalar::G2AffinePoint;
using blst::scalar::G2Point;

Fp fp_from_hex(const char* hex)
{
    Fp out = {};
    if (!blst_fp_from_hexascii(&out,
            reinterpret_cast<const blst_byte*>(hex))) {
        ADD_FAILURE() << "invalid test-vector field element";
    }
    return out;
}

void expect_g1_affine(const G1Point& point, const char* x, const char* y)
{
    G1AffinePoint actual = {};
    G1Point::to_affine(actual, point);
    const G1AffinePoint expected = {fp_from_hex(x), fp_from_hex(y)};
    EXPECT_TRUE(G1AffinePoint::equal(actual, expected));
    EXPECT_TRUE(G1Point::in_subgroup(point));
}

void expect_g2_affine(const G2Point& point,
    const char* x0, const char* x1, const char* y0, const char* y1)
{
    G2AffinePoint actual = {};
    G2Point::to_affine(actual, point);
    const G2AffinePoint expected = {
        Fp2{fp_from_hex(x0), fp_from_hex(x1)},
        Fp2{fp_from_hex(y0), fp_from_hex(y1)},
    };
    EXPECT_TRUE(G2AffinePoint::equal(actual, expected));
    EXPECT_TRUE(G2Point::in_subgroup(point));
}

TEST(MapToGroupAudit, Rfc9380G1EncodeEmptyMessage)
{
    static constexpr unsigned char dst[] =
        "QUUX-V01-CS02-with-BLS12381G1_XMD:SHA-256_SSWU_NU_";
    G1Point point = {};
    blst_encode_to_g1(&point, nullptr, 0, dst, sizeof(dst) - 1, nullptr, 0);

    expect_g1_affine(point,
        "0x184bb665c37ff561a89ec2122dd343f20e0f4cbcaec84e3c3052ea81d1834e192c426074b02ed3dca4e7676ce4ce48ba",
        "0x04407b8d35af4dacc809927071fc0405218f1401a6d15af775810e4e460064bcc9468beeba82fdc751be70476c888bf3");
}

TEST(MapToGroupAudit, Rfc9380G1HashEmptyMessage)
{
    static constexpr unsigned char dst[] =
        "QUUX-V01-CS02-with-BLS12381G1_XMD:SHA-256_SSWU_RO_";
    G1Point point = {};
    blst_hash_to_g1(&point, nullptr, 0, dst, sizeof(dst) - 1, nullptr, 0);

    expect_g1_affine(point,
        "0x052926add2207b76ca4fa57a8734416c8dc95e24501772c814278700eed6d1e4e8cf62d9c09db0fac349612b759e79a1",
        "0x08ba738453bfed09cb546dbb0783dbb3a5f1f566ed67bb6be0e8c67e2e81a4cc68ee29813bb7994998f3eae0c9c6a265");
}

TEST(MapToGroupAudit, Rfc9380G2EncodeEmptyMessage)
{
    static constexpr unsigned char dst[] =
        "QUUX-V01-CS02-with-BLS12381G2_XMD:SHA-256_SSWU_NU_";
    G2Point point = {};
    blst_encode_to_g2(&point, nullptr, 0, dst, sizeof(dst) - 1, nullptr, 0);

    expect_g2_affine(point,
        "0x00e7f4568a82b4b7dc1f14c6aaa055edf51502319c723c4dc2688c7fe5944c213f510328082396515734b6612c4e7bb7",
        "0x126b855e9e69b1f691f816e48ac6977664d24d99f8724868a184186469ddfd4617367e94527d4b74fc86413483afb35b",
        "0x0caead0fd7b6176c01436833c79d305c78be307da5f6af6c133c47311def6ff1e0babf57a0fb5539fce7ee12407b0a42",
        "0x1498aadcf7ae2b345243e281ae076df6de84455d766ab6fcdaad71fab60abb2e8b980a440043cd305db09d283c895e3d");
}

TEST(MapToGroupAudit, Rfc9380G2HashEmptyMessage)
{
    static constexpr unsigned char dst[] =
        "QUUX-V01-CS02-with-BLS12381G2_XMD:SHA-256_SSWU_RO_";
    G2Point point = {};
    blst_hash_to_g2(&point, nullptr, 0, dst, sizeof(dst) - 1, nullptr, 0);

    expect_g2_affine(point,
        "0x0141ebfbdca40eb85b87142e130ab689c673cf60f1a3e98d69335266f30d9b8d4ac44c1038e9dcdd5393faf5c41fb78a",
        "0x05cb8437535e20ecffaef7752baddf98034139c38452458baeefab379ba13dff5bf5dd71b72418717047f5b0f37da03d",
        "0x0503921d7f6a12805e72940b963c0cf3471c7b2a524950ca195d11062ee75ec076daf2d4bc358c4b190c0c98064fdd92",
        "0x12424ac32561493f3fe3c260708a12b7c620e7be00099a974e259ddc7d1f6395c3c811cdd19f1e8dbf3e9ecfdcbab8d6");
}

TEST(MapToGroupAudit, ZeroFieldExercisesExceptionalSswuPath)
{
    G1Point g1 = {};
    blst::scalar::map_to_curve(g1, Fp::zero(), nullptr);
    EXPECT_TRUE(G1Point::on_curve(g1));
    EXPECT_TRUE(G1Point::in_subgroup(g1));

    G2Point g2 = {};
    blst::scalar::map_to_curve(g2, Fp2::zero(), nullptr);
    EXPECT_TRUE(G2Point::on_curve(g2));
    EXPECT_TRUE(G2Point::in_subgroup(g2));
}

TEST(MapToGroupAudit, IsogenyMapSupportsExactAliasing)
{
    G1Point g1_input = {};
    blst::scalar::map_to_isogenous(g1_input, Fp::one());
    G1Point g1_expected = {};
    blst::scalar::isogeny_map_to(g1_expected, g1_input);
    blst::scalar::isogeny_map_to(g1_input, g1_input);
    EXPECT_TRUE(G1Point::equal(g1_input, g1_expected));

    G2Point g2_input = {};
    blst::scalar::map_to_isogenous(g2_input, Fp2::one());
    G2Point g2_expected = {};
    blst::scalar::isogeny_map_to(g2_expected, g2_input);
    blst::scalar::isogeny_map_to(g2_input, g2_input);
    EXPECT_TRUE(G2Point::equal(g2_input, g2_expected));
}

#if BLST_DEFAULT_AVX
TEST(MapToGroupAudit, SelectedAvxBridgePreservesProjectivePoints)
{
    static constexpr unsigned char msg[] = "inversion-free bridge";
    static constexpr unsigned char g1_dst[] =
        "BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_TEST";

    G1Point scalar_g1 = {};
    blst::scalar::Hash_to_curve(
        scalar_g1, msg, sizeof(msg) - 1, g1_dst, sizeof(g1_dst) - 1, nullptr, 0);
    const blst::avx::G1Point avx_g1 =
        blst::avx::avx_map_to_G_detail::from_scalar(scalar_g1);
    const G1Point round_trip_g1 =
        blst::avx::avx_map_to_G_detail::to_scalar(avx_g1);
    EXPECT_TRUE(G1Point::equal(scalar_g1, round_trip_g1));

    const blst::avx::G1Point avx_identity =
        blst::avx::avx_map_to_G_detail::from_scalar(G1Point::identity());
    EXPECT_TRUE(G1Point::is_inf(
        blst::avx::avx_map_to_G_detail::to_scalar(avx_identity)));

    static constexpr unsigned char g2_dst[] =
        "BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_TEST";
    G2Point scalar_g2 = {};
    blst::scalar::Hash_to_curve(
        scalar_g2, msg, sizeof(msg) - 1, g2_dst, sizeof(g2_dst) - 1, nullptr, 0);
    const blst::avx::G2Point avx_g2 =
        blst::avx::avx_map_to_G_detail::from_scalar(scalar_g2);
    const G2Point round_trip_g2 =
        blst::avx::avx_map_to_G_detail::to_scalar(avx_g2);
    EXPECT_TRUE(G2Point::equal(scalar_g2, round_trip_g2));

    const blst::avx::G2Point avx_g2_identity =
        blst::avx::avx_map_to_G_detail::from_scalar(G2Point::identity());
    EXPECT_TRUE(G2Point::is_inf(
        blst::avx::avx_map_to_G_detail::to_scalar(avx_g2_identity)));
}
#endif

} // namespace
