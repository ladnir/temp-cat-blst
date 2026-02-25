#include <array>
#include <span>

#include "gtest/gtest.h"

#include "blst/curve/pairing_avx.h"
#include "blst/curve/point.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

using namespace blst;

namespace
{
void make_pairing_inputs(
    std::array<G2AffinePoint, 2>& q,
    std::array<G1AffinePoint, 2>& p)
{
    G1Point g1 = G1Point::generator();
    G2Point g2 = G2Point::generator();
    G1Point g1_2;
    G2Point g2_2;
    G1Point::doubling(g1_2, g1);
    G2Point::doubling(g2_2, g2);

    G1Point::to_affine(p[0], g1);
    G1Point::to_affine(p[1], g1_2);
    G2Point::to_affine(q[0], g2);
    G2Point::to_affine(q[1], g2_2);
}
} // namespace

TEST(PairingAvx, MillerLoopSpanMatchesPointer)
{
    std::array<G2AffinePoint, 2> q{};
    std::array<G1AffinePoint, 2> p{};
    make_pairing_inputs(q, p);

    GT ptr_f{};
    GT span_f{};

    avx::miller_loop_n(ptr_f, q.data(), p.data(), q.size());
    avx::miller_loop_n(span_f, std::span<const G2AffinePoint>(q), std::span<const G1AffinePoint>(p));
    EXPECT_TRUE(GT::equal(ptr_f, span_f));

    GT ptr_exp{};
    GT span_exp{};
    avx::final_exp(ptr_exp, ptr_f);
    avx::final_exp(span_exp, span_f);
    EXPECT_TRUE(GT::equal(ptr_exp, span_exp));
}

TEST(PairingAvx, OptimalAteSpanMatchesPointer)
{
    std::array<G2AffinePoint, 2> q{};
    std::array<G1AffinePoint, 2> p{};
    make_pairing_inputs(q, p);

    GT ptr{};
    GT span{};
    avx::optimal_ate_pairing(ptr, q.data(), p.data(), q.size());
    avx::optimal_ate_pairing(span, std::span<const G2AffinePoint>(q), std::span<const G1AffinePoint>(p));

    EXPECT_TRUE(GT::equal(ptr, span));
}

TEST(PairingAvx, EmptyInputReturnsIdentity)
{
    GT ptr{};
    GT span{};

    avx::miller_loop_n(ptr, static_cast<const G2AffinePoint*>(nullptr), static_cast<const G1AffinePoint*>(nullptr), 0);
    avx::miller_loop_n(span, std::span<const G2AffinePoint>(), std::span<const G1AffinePoint>());

    EXPECT_TRUE(GT::is_identity(ptr));
    EXPECT_TRUE(GT::is_identity(span));
}

#endif
