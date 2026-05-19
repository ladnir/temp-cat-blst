#include <gtest/gtest.h>

#include "blst/avx/curve/point.h"

#include <array>
#include <cstring>

namespace {

using ScalarBytes = std::array<blst::byte, 32>;
using G1JjOp = void (*)(blst::avx::G1Point&, const blst::avx::G1Point&, const blst::avx::G1Point&);
using G2JjOp = void (*)(blst::avx::G2Point&, const blst::avx::G2Point&, const blst::avx::G2Point&);
using G1JaOp = void (*)(blst::avx::G1Point&, const blst::avx::G1Point&, const blst::avx::G1AffinePoint&);
using G2JaOp = void (*)(blst::avx::G2Point&, const blst::avx::G2Point&, const blst::avx::G2AffinePoint&);

ScalarBytes scalar_bytes(uint64_t seed)
{
    ScalarBytes out = {};
    uint64_t x = seed * 0x9e3779b97f4a7c15ULL + 0xbf58476d1ce4e5b9ULL;
    for (size_t i = 0; i < out.size(); ++i) {
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        out[i] = static_cast<blst::byte>(x >> ((i & 7u) * 8u));
    }
    out[31] &= 0x0f;
    return out;
}

blst::scalar::G1Point scalar_g1(uint64_t seed)
{
    const ScalarBytes scalar = scalar_bytes(seed);
    blst::scalar::G1Point out = {};
    blst::scalar::G1Point::mult(out, blst::scalar::G1Point::generator(), scalar.data(), 256);
    return out;
}

blst::scalar::G2Point scalar_g2(uint64_t seed)
{
    const ScalarBytes scalar = scalar_bytes(seed);
    blst::scalar::G2Point out = {};
    blst::scalar::G2Point::mult(out, blst::scalar::G2Point::generator(), scalar.data(), 256);
    return out;
}

blst::avx::G1Point avx_from_scalar(const blst::scalar::G1Point& point)
{
    blst::scalar::G1AffinePoint affine = {};
    blst::scalar::G1Point::to_affine(affine, point);
    const blst::avx::G1AffinePoint avx_affine = {
        blst::avx::avx_curve_detail::scalar_fp_to_ring(affine.X),
        blst::avx::avx_curve_detail::scalar_fp_to_ring(affine.Y),
    };
    blst::avx::G1Point out = {};
    blst::avx::G1Point::from_affine(out, avx_affine);
    return out;
}

blst::avx::G2Point avx_from_scalar(const blst::scalar::G2Point& point)
{
    blst::scalar::G2AffinePoint affine = {};
    blst::scalar::G2Point::to_affine(affine, point);
    const blst::avx::G2AffinePoint avx_affine = {
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(affine.X),
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(affine.Y),
    };
    blst::avx::G2Point out = {};
    blst::avx::G2Point::from_affine(out, avx_affine);
    return out;
}

blst::scalar::G1Point scalar_from_avx(const blst::avx::G1Point& point)
{
    blst::avx::G1AffinePoint affine = {};
    blst::avx::G1Point::to_affine(affine, point);
    const blst::scalar::G1AffinePoint scalar_affine = {
        blst::avx::avx_curve_detail::ring_fp_to_scalar(affine.X),
        blst::avx::avx_curve_detail::ring_fp_to_scalar(affine.Y),
    };
    blst::scalar::G1Point out = {};
    blst::scalar::G1Point::from_affine(out, scalar_affine);
    return out;
}

blst::scalar::G2Point scalar_from_avx(const blst::avx::G2Point& point)
{
    blst::avx::G2AffinePoint affine = {};
    blst::avx::G2Point::to_affine(affine, point);
    const blst::scalar::G2AffinePoint scalar_affine = {
        blst::avx::avx_curve_detail::ring_fp2_to_scalar(affine.X),
        blst::avx::avx_curve_detail::ring_fp2_to_scalar(affine.Y),
    };
    blst::scalar::G2Point out = {};
    blst::scalar::G2Point::from_affine(out, scalar_affine);
    return out;
}

blst::scalar::G1AffinePoint scalar_from_avx(const blst::avx::G1AffinePoint& point)
{
    return {
        blst::avx::avx_curve_detail::ring_fp_to_scalar(point.X),
        blst::avx::avx_curve_detail::ring_fp_to_scalar(point.Y),
    };
}

blst::scalar::G2AffinePoint scalar_from_avx(const blst::avx::G2AffinePoint& point)
{
    return {
        blst::avx::avx_curve_detail::ring_fp2_to_scalar(point.X),
        blst::avx::avx_curve_detail::ring_fp2_to_scalar(point.Y),
    };
}

void expect_same(const blst::scalar::G1Point& expected, const blst::avx::G1Point& actual)
{
    unsigned char expected_bytes[96] = {};
    unsigned char actual_bytes[96] = {};
    blst::scalar::G1Point::serialize(expected_bytes, expected);
    const blst::scalar::G1Point actual_scalar = scalar_from_avx(actual);
    blst::scalar::G1Point::serialize(actual_bytes, actual_scalar);
    EXPECT_EQ(0, std::memcmp(expected_bytes, actual_bytes, sizeof(expected_bytes)));
}

void expect_same(const blst::scalar::G2Point& expected, const blst::avx::G2Point& actual)
{
    unsigned char expected_bytes[192] = {};
    unsigned char actual_bytes[192] = {};
    blst::scalar::G2Point::serialize(expected_bytes, expected);
    const blst::scalar::G2Point actual_scalar = scalar_from_avx(actual);
    blst::scalar::G2Point::serialize(actual_bytes, actual_scalar);
    EXPECT_EQ(0, std::memcmp(expected_bytes, actual_bytes, sizeof(expected_bytes)));
}

void expect_g1_jj_matches_scalar(
    G1JjOp op,
    const blst::avx::G1Point& a,
    const blst::avx::G1Point& b)
{
    blst::avx::G1Point actual = {};
    op(actual, a, b);

    blst::scalar::G1Point expected = {};
    const blst::scalar::G1Point scalar_a = scalar_from_avx(a);
    const blst::scalar::G1Point scalar_b = scalar_from_avx(b);
    blst::scalar::G1Point::dadd(expected, scalar_a, scalar_b, nullptr);
    expect_same(expected, actual);
}

void expect_g2_jj_matches_scalar(
    G2JjOp op,
    const blst::avx::G2Point& a,
    const blst::avx::G2Point& b)
{
    blst::avx::G2Point actual = {};
    op(actual, a, b);

    blst::scalar::G2Point expected = {};
    const blst::scalar::G2Point scalar_a = scalar_from_avx(a);
    const blst::scalar::G2Point scalar_b = scalar_from_avx(b);
    blst::scalar::G2Point::dadd(expected, scalar_a, scalar_b, nullptr);
    expect_same(expected, actual);
}

void expect_g1_ja_matches_scalar(
    G1JaOp op,
    const blst::avx::G1Point& a,
    const blst::avx::G1AffinePoint& b)
{
    blst::avx::G1Point actual = {};
    op(actual, a, b);

    blst::scalar::G1Point expected = {};
    const blst::scalar::G1Point scalar_a = scalar_from_avx(a);
    const blst::scalar::G1AffinePoint scalar_b = scalar_from_avx(b);
    blst::scalar::G1Point::dadd(expected, scalar_a, scalar_b);
    expect_same(expected, actual);
}

void expect_g2_ja_matches_scalar(
    G2JaOp op,
    const blst::avx::G2Point& a,
    const blst::avx::G2AffinePoint& b)
{
    blst::avx::G2Point actual = {};
    op(actual, a, b);

    blst::scalar::G2Point expected = {};
    const blst::scalar::G2Point scalar_a = scalar_from_avx(a);
    const blst::scalar::G2AffinePoint scalar_b = scalar_from_avx(b);
    blst::scalar::G2Point::dadd(expected, scalar_a, scalar_b);
    expect_same(expected, actual);
}

}  // namespace

TEST(AvxAddCornerCases, G1JjAddAndDaddMatchScalar)
{
    const blst::avx::G1Point a = avx_from_scalar(scalar_g1(11));
    const blst::avx::G1Point b = avx_from_scalar(scalar_g1(12));
    const blst::avx::G1Point inf = blst::avx::G1Point::infinity();
    blst::avx::G1Point neg_a = {};
    blst::avx::G1Point::neg(neg_a, a);
    blst::avx::G1Point nontrivial_z = {};
    blst::avx::G1Point::add(nontrivial_z, a, b);

    const std::array<G1JjOp, 2> ops = {
        static_cast<G1JjOp>(&blst::avx::G1Point::add),
        static_cast<G1JjOp>(&blst::avx::G1Point::dadd),
    };
    for (const G1JjOp op : ops) {
        expect_g1_jj_matches_scalar(op, a, b);
        expect_g1_jj_matches_scalar(op, a, a);
        expect_g1_jj_matches_scalar(op, a, neg_a);
        expect_g1_jj_matches_scalar(op, inf, a);
        expect_g1_jj_matches_scalar(op, a, inf);
        expect_g1_jj_matches_scalar(op, inf, inf);
        expect_g1_jj_matches_scalar(op, nontrivial_z, b);
    }

    blst::avx::G1Point alias = a;
    blst::avx::G1Point::add(alias, alias, b);
    blst::scalar::G1Point expected = {};
    blst::scalar::G1Point::dadd(expected, scalar_from_avx(a), scalar_from_avx(b), nullptr);
    expect_same(expected, alias);
}

TEST(AvxAddCornerCases, G2JjAddAndDaddMatchScalar)
{
    const blst::avx::G2Point a = avx_from_scalar(scalar_g2(21));
    const blst::avx::G2Point b = avx_from_scalar(scalar_g2(22));
    const blst::avx::G2Point inf = blst::avx::G2Point::infinity();
    blst::avx::G2Point neg_a = {};
    blst::avx::G2Point::neg(neg_a, a);
    blst::avx::G2Point nontrivial_z = {};
    blst::avx::G2Point::add(nontrivial_z, a, b);

    const std::array<G2JjOp, 2> ops = {
        static_cast<G2JjOp>(&blst::avx::G2Point::add),
        static_cast<G2JjOp>(&blst::avx::G2Point::dadd),
    };
    for (const G2JjOp op : ops) {
        expect_g2_jj_matches_scalar(op, a, b);
        expect_g2_jj_matches_scalar(op, a, a);
        expect_g2_jj_matches_scalar(op, a, neg_a);
        expect_g2_jj_matches_scalar(op, inf, a);
        expect_g2_jj_matches_scalar(op, a, inf);
        expect_g2_jj_matches_scalar(op, inf, inf);
        expect_g2_jj_matches_scalar(op, nontrivial_z, b);
    }

    blst::avx::G2Point alias = a;
    blst::avx::G2Point::add(alias, alias, b);
    blst::scalar::G2Point expected = {};
    blst::scalar::G2Point::dadd(expected, scalar_from_avx(a), scalar_from_avx(b), nullptr);
    expect_same(expected, alias);
}

TEST(AvxAddCornerCases, G1JaAddAndDaddHonorAffineIdentity)
{
    const blst::avx::G1Point a = avx_from_scalar(scalar_g1(31));
    const blst::avx::G1Point inf = blst::avx::G1Point::infinity();
    blst::avx::G1AffinePoint a_affine = {};
    blst::avx::G1Point::to_affine(a_affine, a);
    blst::avx::G1AffinePoint neg_affine = a_affine;
    blst::avx::G1AffinePoint::cneg(neg_affine, 1);
    const blst::avx::G1AffinePoint affine_inf = {};

    const std::array<G1JaOp, 2> ops = {
        static_cast<G1JaOp>(&blst::avx::G1Point::add),
        static_cast<G1JaOp>(&blst::avx::G1Point::dadd),
    };
    for (const G1JaOp op : ops) {
        expect_g1_ja_matches_scalar(op, a, a_affine);
        expect_g1_ja_matches_scalar(op, a, neg_affine);
        expect_g1_ja_matches_scalar(op, inf, a_affine);
        expect_g1_ja_matches_scalar(op, a, affine_inf);
    }
}

TEST(AvxAddCornerCases, G2JaAddAndDaddHonorAffineIdentity)
{
    const blst::avx::G2Point a = avx_from_scalar(scalar_g2(41));
    const blst::avx::G2Point inf = blst::avx::G2Point::infinity();
    blst::avx::G2AffinePoint a_affine = {};
    blst::avx::G2Point::to_affine(a_affine, a);
    blst::avx::G2AffinePoint neg_affine = a_affine;
    blst::avx::G2AffinePoint::cneg(neg_affine, 1);
    const blst::avx::G2AffinePoint affine_inf = {};

    const std::array<G2JaOp, 2> ops = {
        static_cast<G2JaOp>(&blst::avx::G2Point::add),
        static_cast<G2JaOp>(&blst::avx::G2Point::dadd),
    };
    for (const G2JaOp op : ops) {
        expect_g2_ja_matches_scalar(op, a, a_affine);
        expect_g2_ja_matches_scalar(op, a, neg_affine);
        expect_g2_ja_matches_scalar(op, inf, a_affine);
        expect_g2_ja_matches_scalar(op, a, affine_inf);
    }
}

TEST(AvxAddCornerCases, FromAffineZeroProducesIdentity)
{
    blst::avx::G1Point g1 = {};
    blst::avx::G1Point::from_affine(g1, blst::avx::G1AffinePoint{});
    EXPECT_TRUE(blst::avx::G1Point::is_identity(g1));

    blst::avx::G2Point g2 = {};
    blst::avx::G2Point::from_affine(g2, blst::avx::G2AffinePoint{});
    EXPECT_TRUE(blst::avx::G2Point::is_identity(g2));
}

TEST(AvxAddCornerCases, CnegMatchesScalarSemantics)
{
    blst::avx::G1AffinePoint g1_affine = {};
    blst::avx::G1Point::to_affine(g1_affine, avx_from_scalar(scalar_g1(51)));
    const blst::avx::G1AffinePoint g1_original = g1_affine;
    blst::avx::G1AffinePoint::cneg(g1_affine, 0);
    EXPECT_TRUE(blst::avx::G1AffinePoint::equal(g1_affine, g1_original));
    blst::avx::G1AffinePoint::cneg(g1_affine, 1);
    blst::avx::G1AffinePoint g1_expected = g1_original;
    g1_expected.Y = blst::avx::avx_curve_detail::curve_ring().standard_negate(g1_expected.Y);
    EXPECT_TRUE(blst::avx::G1AffinePoint::equal(g1_affine, g1_expected));

    blst::avx::G2AffinePoint g2_affine = {};
    blst::avx::G2Point::to_affine(g2_affine, avx_from_scalar(scalar_g2(52)));
    const blst::avx::G2AffinePoint g2_original = g2_affine;
    blst::avx::G2AffinePoint::cneg(g2_affine, 0);
    EXPECT_TRUE(blst::avx::G2AffinePoint::equal(g2_affine, g2_original));
    blst::avx::G2AffinePoint::cneg(g2_affine, 1);
    blst::avx::G2AffinePoint g2_expected = g2_original;
    g2_expected.Y = blst::avx::avx_curve_detail::curve_fp2_ring().standard_negate(g2_expected.Y);
    EXPECT_TRUE(blst::avx::G2AffinePoint::equal(g2_affine, g2_expected));

    blst::avx::G1AffinePoint g1_zero = {};
    blst::avx::G1AffinePoint::cneg(g1_zero, 1);
    EXPECT_TRUE(blst::avx::G1AffinePoint::is_inf(g1_zero));

    blst::avx::G2AffinePoint g2_zero = {};
    blst::avx::G2AffinePoint::cneg(g2_zero, 1);
    EXPECT_TRUE(blst::avx::G2AffinePoint::is_inf(g2_zero));
}

TEST(AvxAddCornerCases, DaddWithA4HandlesIdentity)
{
    const blst::avx::G1Point g1 = avx_from_scalar(scalar_g1(61));
    const blst::avx::G1Point g1_inf = blst::avx::G1Point::infinity();
    blst::avx::G1Point g1_actual = {};
    blst::avx::G1Point::dadd(g1_actual, g1, g1_inf, &blst::avx::get_BLST_Aprime_E1());
    expect_same(scalar_from_avx(g1), g1_actual);

    const blst::avx::G2Point g2 = avx_from_scalar(scalar_g2(62));
    const blst::avx::G2Point g2_inf = blst::avx::G2Point::infinity();
    blst::avx::G2Point g2_actual = {};
    blst::avx::G2Point::dadd(g2_actual, g2_inf, g2, &blst::avx::get_BLST_Aprime_E2());
    expect_same(scalar_from_avx(g2), g2_actual);
}
