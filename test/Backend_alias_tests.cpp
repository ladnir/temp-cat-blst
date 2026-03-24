#include <gtest/gtest.h>

#include <type_traits>

#include "blst/curve/GT.h"
#include "blst/curve/map_to_G.h"
#include "blst/curve/point.h"

using namespace blst;

TEST(BackendAlias, TopLevelAliasesMatchConfig)
{
#if BLST_DEFAULT_AVX
    EXPECT_TRUE((std::is_same_v<G1Point, avx::G1Point>));
    EXPECT_TRUE((std::is_same_v<G1AffinePoint, avx::G1AffinePoint>));
    EXPECT_TRUE((std::is_same_v<G2Point, avx::G2Point>));
    EXPECT_TRUE((std::is_same_v<G2AffinePoint, avx::G2AffinePoint>));
    EXPECT_TRUE((std::is_same_v<GT, avx::GT>));
#else
    EXPECT_TRUE((std::is_same_v<G1Point, scalar::G1Point>));
    EXPECT_TRUE((std::is_same_v<G1AffinePoint, scalar::G1AffinePoint>));
    EXPECT_TRUE((std::is_same_v<G2Point, scalar::G2Point>));
    EXPECT_TRUE((std::is_same_v<G2AffinePoint, scalar::G2AffinePoint>));
    EXPECT_TRUE((std::is_same_v<GT, scalar::GT>));
#endif
}

TEST(BackendAlias, TopLevelHashToCurveCompilesAndProducesSubgroupPoints)
{
    static constexpr unsigned char msg[] = "top-level backend alias";
    static constexpr unsigned char dst_g1[] = "BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_TEST";
    static constexpr unsigned char dst_g2[] = "BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_TEST";

    G1Point g1 = {};
    Hash_to_curve(g1, msg, sizeof(msg) - 1, dst_g1, sizeof(dst_g1) - 1, nullptr, 0);
    EXPECT_TRUE(G1Point::on_curve(g1));
    EXPECT_TRUE(G1Point::in_subgroup(g1));

    G2Point g2 = {};
    Hash_to_curve(g2, msg, sizeof(msg) - 1, dst_g2, sizeof(dst_g2) - 1, nullptr, 0);
    EXPECT_TRUE(G2Point::on_curve(g2));
    EXPECT_TRUE(G2Point::in_subgroup(g2));
}

TEST(BackendAlias, TopLevelHashToCurveRoundTripsThroughAffine)
{
    static constexpr unsigned char msg_g1[] = "public-api-g1";
    static constexpr unsigned char dst_g1[] = "BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_TEST";

    G1Point g1 = {};
    Hash_to_curve(g1, msg_g1, sizeof(msg_g1) - 1, dst_g1, sizeof(dst_g1) - 1, nullptr, 0);

    G1AffinePoint g1_affine = {};
    G1Point::to_affine(g1_affine, g1);
    EXPECT_TRUE(G1AffinePoint::in_subgroup(g1_affine));

    unsigned char g1_compressed[48] = {};
    G1Point::compress(g1_compressed, g1);

    G1AffinePoint g1_decoded = {};
    ASSERT_EQ(BLST_SUCCESS, G1AffinePoint::uncompress(g1_decoded, g1_compressed));

    G1Point g1_roundtrip = {};
    G1Point::from_affine(g1_roundtrip, g1_decoded);
    EXPECT_TRUE(G1Point::equal(g1, g1_roundtrip));

    static constexpr unsigned char msg_g2[] = "public-api-g2";
    static constexpr unsigned char dst_g2[] = "BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_TEST";

    G2Point g2 = {};
    Hash_to_curve(g2, msg_g2, sizeof(msg_g2) - 1, dst_g2, sizeof(dst_g2) - 1, nullptr, 0);

    G2AffinePoint g2_affine = {};
    G2Point::to_affine(g2_affine, g2);
    EXPECT_TRUE(G2AffinePoint::in_subgroup(g2_affine));

    unsigned char g2_compressed[96] = {};
    G2Point::compress(g2_compressed, g2);

    G2AffinePoint g2_decoded = {};
    ASSERT_EQ(BLST_SUCCESS, G2AffinePoint::uncompress(g2_decoded, g2_compressed));

    G2Point g2_roundtrip = {};
    G2Point::from_affine(g2_roundtrip, g2_decoded);
    EXPECT_TRUE(G2Point::equal(g2, g2_roundtrip));
}

TEST(BackendAlias, TopLevelPairingSurfaceMatchesBilinearDoubling)
{
    const G1Point g1 = G1Point::generator();
    const G2Point g2 = G2Point::generator();

    G1Point g1_double = {};
    G1Point::doubling(g1_double, g1);

    G2Point g2_double = {};
    G2Point::doubling(g2_double, g2);

    G1AffinePoint g1_affine = {};
    G1AffinePoint g1_double_affine = {};
    G1Point::to_affine(g1_affine, g1);
    G1Point::to_affine(g1_double_affine, g1_double);

    G2AffinePoint g2_affine = {};
    G2AffinePoint g2_double_affine = {};
    G2Point::to_affine(g2_affine, g2);
    G2Point::to_affine(g2_double_affine, g2_double);

    GT base = {};
    GT::miller_loop(base, &g2_affine, &g1_affine);
    GT::final_mult(base, base);

    GT lhs = {};
    GT::miller_loop(lhs, &g2_affine, &g1_double_affine);
    GT::final_mult(lhs, lhs);

    GT rhs = {};
    GT::miller_loop(rhs, &g2_double_affine, &g1_affine);
    GT::final_mult(rhs, rhs);

    GT doubled = {};
    GT::doubling(doubled, base);

    EXPECT_TRUE(GT::equal(lhs, rhs));
    EXPECT_TRUE(GT::equal(lhs, doubled));
}
