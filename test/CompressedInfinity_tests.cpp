#include <array>

#include <gtest/gtest.h>

#include "bindings/blst_avx.h"
#include "bindings/blst_scalar.h"
#include "blst/scalar/curve/point.h"

namespace
{
template <std::size_t N>
std::array<byte, N> compressed_infinity()
{
    std::array<byte, N> encoded{};
    encoded[0] = 0xc0; // compressed + infinity, with sign clear
    return encoded;
}
}

TEST(CompressedInfinity, ScalarCApiAcceptsExactG1AndG2Encodings)
{
    const auto g1_encoded = compressed_infinity<48>();
    blst_p1_affine g1 = {};
    ASSERT_EQ(BLST_SUCCESS, blst_p1_uncompress(&g1, g1_encoded.data()));
    EXPECT_TRUE(blst_p1_affine_is_inf(&g1));

    const auto g2_encoded = compressed_infinity<96>();
    blst_p2_affine g2 = {};
    ASSERT_EQ(BLST_SUCCESS, blst_p2_uncompress(&g2, g2_encoded.data()));
    EXPECT_TRUE(blst_p2_affine_is_inf(&g2));
}

TEST(CompressedInfinity, NativeCppAcceptsExactG1AndG2Encodings)
{
    const auto g1_encoded = compressed_infinity<48>();
    blst::scalar::G1AffinePoint g1 = {};
    ASSERT_EQ(
        blst::BLST_SUCCESS,
        blst::scalar::G1AffinePoint::uncompress(g1, g1_encoded.data()));
    EXPECT_TRUE(blst::scalar::G1AffinePoint::is_inf(g1));

    const auto g2_encoded = compressed_infinity<96>();
    blst::scalar::G2AffinePoint g2 = {};
    ASSERT_EQ(
        blst::BLST_SUCCESS,
        blst::scalar::G2AffinePoint::uncompress(g2, g2_encoded.data()));
    EXPECT_TRUE(blst::scalar::G2AffinePoint::is_inf(g2));
}

TEST(CompressedInfinity, AvxCApiAcceptsExactG1AndG2EncodingsWhenAvailable)
{
    if (!blst_avx_enabled()) {
        GTEST_SKIP() << "AVX backend is not enabled in this build";
    }

    const auto g1_encoded = compressed_infinity<48>();
    blst_avx_p1 g1 = {};
    ASSERT_EQ(
        BLST_SUCCESS,
        blst_avx_p1_uncompress(&g1, g1_encoded.data()));
    EXPECT_TRUE(blst_avx_p1_is_inf(&g1));

    const auto g2_encoded = compressed_infinity<96>();
    blst_avx_p2 g2 = {};
    ASSERT_EQ(
        BLST_SUCCESS,
        blst_avx_p2_uncompress(&g2, g2_encoded.data()));
    EXPECT_TRUE(blst_avx_p2_is_inf(&g2));
}
