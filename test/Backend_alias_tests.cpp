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
