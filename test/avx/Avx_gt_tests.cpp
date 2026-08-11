#include <gtest/gtest.h>

#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/point.h"
#include "blst/avx/curve/GT.h"

#include <array>
#include <cstring>

namespace {

blst::scalar::GT scalar_gt_base()
{
    const blst::scalar::G1Point g1 = blst::scalar::G1Point::generator();
    const blst::scalar::G2Point g2 = blst::scalar::G2Point::generator();
    blst::scalar::G1AffinePoint g1_affine = {};
    blst::scalar::G2AffinePoint g2_affine = {};
    blst::scalar::G1Point::to_affine(g1_affine, g1);
    blst::scalar::G2Point::to_affine(g2_affine, g2);

    blst::scalar::MillerResult miller = {};
    blst::scalar::GT out = {};
    blst::scalar::GT::miller_loop(miller, &g2_affine, &g1_affine);
    blst::scalar::GT::final_exp(out, miller);
    return out;
}

blst::avx::GT avx_from_scalar(const blst::scalar::GT& value)
{
    unsigned char bytes[12 * 48] = {};
    blst::scalar::GT::serialize(bytes, value);

    blst::avx::GT out = {};
    EXPECT_TRUE(blst::avx::GT::deserialize(out, bytes));
    return out;
}

void expect_same(const blst::scalar::GT& scalar, const blst::avx::GT& avx)
{
    unsigned char scalar_bytes[12 * 48] = {};
    unsigned char avx_bytes[12 * 48] = {};
    blst::scalar::GT::serialize(scalar_bytes, scalar);
    blst::avx::GT::serialize(avx_bytes, avx);
    EXPECT_EQ(0, std::memcmp(scalar_bytes, avx_bytes, sizeof(scalar_bytes)));
}

void expect_same(const blst::avx::GT& left, const blst::avx::GT& right)
{
    unsigned char left_bytes[12 * 48] = {};
    unsigned char right_bytes[12 * 48] = {};
    blst::avx::GT::serialize(left_bytes, left);
    blst::avx::GT::serialize(right_bytes, right);
    EXPECT_EQ(0, std::memcmp(left_bytes, right_bytes, sizeof(left_bytes)));
}

void expect_same_after_bound_sensitive_ops(
    const blst::scalar::GT& scalar,
    const blst::avx::GT& avx)
{
    expect_same(scalar, avx);

    blst::scalar::GT scalar_conj = scalar;
    blst::avx::GT avx_conj = avx;
    blst::scalar::GT::conjugate(scalar_conj);
    blst::avx::GT::conjugate(avx_conj);
    expect_same(scalar_conj, avx_conj);

    blst::scalar::GT scalar_neg = {};
    blst::avx::GT avx_neg = {};
    blst::scalar::GT::neg(scalar_neg, scalar);
    blst::avx::GT::neg(avx_neg, avx);
    expect_same(scalar_neg, avx_neg);

    const blst::scalar::GT scalar_base = scalar_gt_base();
    const blst::avx::GT avx_base = avx_from_scalar(scalar_base);
    blst::scalar::GT scalar_add = {};
    blst::avx::GT avx_add = {};
    blst::scalar::GT::add(scalar_add, scalar, scalar_base);
    blst::avx::GT::add(avx_add, avx, avx_base);
    expect_same(scalar_add, avx_add);
}

bool avx_miller_equal(
    const blst::avx::MillerResult& left,
    const blst::avx::MillerResult& right)
{
    return blst::Fp12::equal(
        blst::avx::avx_gt_detail::to_scalar(left.data),
        blst::avx::avx_gt_detail::to_scalar(right.data));
}

}  // namespace

TEST(AvxGt, NarrowToStandardMatchesMulByOneReference)
{
    const blst::avx::GT avx_base = avx_from_scalar(scalar_gt_base());
    const auto wide = blst::avx::avx_gt_detail::fp12ring().fp12_flat_cyclotomic_sqr(
        avx_base.data, blst::avx::avx_gt_detail::ring());

    blst::avx::GT reference = {};
    reference.data = blst::avx::avx_gt_detail::to_standard(wide);

    blst::avx::GT narrowed = {};
    narrowed.data = blst::avx::avx_gt_detail::narrow_to_standard(wide);

    expect_same(reference, narrowed);
}

TEST(AvxGt, DoublingMatchesScalar)
{
    const blst::scalar::GT scalar_base = scalar_gt_base();
    const blst::avx::GT avx_base = avx_from_scalar(scalar_base);

    blst::scalar::GT scalar_out = {};
    blst::scalar::GT::doubling(scalar_out, scalar_base);

    blst::avx::GT avx_out = {};
    blst::avx::GT::doubling(avx_out, avx_base);

    expect_same_after_bound_sensitive_ops(scalar_out, avx_out);
}

TEST(AvxGt, PreparedLinesPreservePairingIdentityAndCompactBatching)
{
    using blst::avx::G1AffinePoint;
    using blst::avx::G2AffinePoint;
    using blst::avx::GT;
    using blst::avx::MillerResult;
    using blst::avx::PreparedLine;

    const G1AffinePoint p = G1AffinePoint::generator();
    const G2AffinePoint q = G2AffinePoint::generator();
    const G1AffinePoint p_infinity = {};
    const G2AffinePoint q_infinity = {};

    std::array<PreparedLine, 3 * 68> lines = {};
    GT::precompute_lines(lines.data(), q);
    GT::precompute_lines(lines.data() + 68, q_infinity);
    GT::precompute_lines(lines.data() + 2 * 68, q);

    MillerResult p_infinity_result = {};
    GT::miller_loop_lines(p_infinity_result, lines.data(), p_infinity);
    EXPECT_TRUE(avx_miller_equal(p_infinity_result, MillerResult::identity()));

    MillerResult q_infinity_result = {};
    GT::miller_loop_lines(q_infinity_result, lines.data() + 68, p);
    EXPECT_TRUE(avx_miller_equal(q_infinity_result, MillerResult::identity()));

    const std::array<G1AffinePoint, 3> ps = {p_infinity, p, p};
    const std::array<const PreparedLine*, 3> line_ptrs = {
        lines.data(), nullptr, nullptr
    };
    const std::array<const G1AffinePoint*, 3> p_ptrs = {
        ps.data(), nullptr, nullptr
    };
    MillerResult batched = {};
    GT::miller_loop_lines_n(
        batched, line_ptrs.data(), p_ptrs.data(), line_ptrs.size());

    MillerResult expected = {};
    GT::miller_loop_lines(expected, lines.data() + 2 * 68, p);
    EXPECT_TRUE(avx_miller_equal(batched, expected));
}

TEST(AvxGt, PreparedLineValueInitializationUsesZeroCarriers)
{
    const blst::avx::PreparedLine line = {};
    EXPECT_EQ(line.c.is_zero(), 1U);
    EXPECT_EQ(line.x.to_fp2().is_zero(), 1U);
    EXPECT_EQ(line.y.to_fp2().is_zero(), 1U);
}

#if !defined(NDEBUG) && GTEST_HAS_DEATH_TEST
TEST(AvxGt, RejectsInvalidConditionalBitsAndFrobeniusPowersInDebug)
{
    const blst::avx::GT identity = blst::avx::GT::identity();
    EXPECT_DEATH(
        {
            blst::avx::GT left = identity;
            blst::avx::GT right = identity;
            blst::avx::GT::cswap(left, right, 2);
        },
        "");
    EXPECT_DEATH(
        {
            blst::avx::GT out = identity;
            blst::avx::GT::ccopy(out, identity, 2);
        },
        "");
    EXPECT_DEATH(
        {
            blst::avx::GT out = {};
            blst::avx::GT::frobenius_map(out, identity, 0);
        },
        "");
}
#endif

TEST(AvxGt, MultByZMatchesScalar)
{
    const blst::scalar::GT scalar_base = scalar_gt_base();
    const blst::avx::GT avx_base = avx_from_scalar(scalar_base);

    blst::scalar::GT scalar_out = {};
    blst::scalar::GT::mult_by_z(scalar_out, scalar_base);

    blst::avx::GT avx_out = {};
    blst::avx::GT::mult_by_z(avx_out, avx_base);

    expect_same_after_bound_sensitive_ops(scalar_out, avx_out);
}

TEST(AvxGt, MultByZDivBy2MatchesScalar)
{
    const blst::scalar::GT scalar_base = scalar_gt_base();
    const blst::avx::GT avx_base = avx_from_scalar(scalar_base);

    blst::scalar::GT scalar_out = {};
    blst::scalar::GT::mult_by_z_div_by_2(scalar_out, scalar_base);

    blst::avx::GT avx_out = {};
    blst::avx::GT::mult_by_z_div_by_2(avx_out, avx_base);

    expect_same_after_bound_sensitive_ops(scalar_out, avx_out);
}

TEST(AvxGt, FullWidthScalarsMatchScalarBasicExponentiation)
{
    const uint64_t* order = blst::get_BLST_r().l;
    std::array<uint64_t, 4> twice_order = {};
    uint64_t carry = 0;
    for (size_t i = 0; i < twice_order.size(); ++i) {
        const uint64_t doubled = order[i] + order[i];
        const uint64_t next = doubled + carry;
        carry = (doubled < order[i]) | (next < doubled);
        twice_order[i] = next;
    }
    ASSERT_EQ(carry, 0U);

    const std::array<uint64_t, 4> order_value = {
        order[0], order[1], order[2], order[3]
    };
    std::array<uint64_t, 4> order_plus_one = order_value;
    std::array<uint64_t, 4> twice_order_minus_one = twice_order;
    ++order_plus_one[0];
    --twice_order_minus_one[0];

    const std::array<std::array<uint64_t, 4>, 5> scalars = {{
        order_value,
        order_plus_one,
        twice_order_minus_one,
        twice_order,
        {UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX},
    }};
    const std::array<const char*, 5> names = {
        "r", "r + 1", "2r - 1", "2r", "2^256 - 1"
    };

    const blst::scalar::GT scalar_base = scalar_gt_base();
    const blst::avx::GT avx_base = avx_from_scalar(scalar_base);
    for (size_t i = 0; i < scalars.size(); ++i) {
        SCOPED_TRACE(names[i]);
        unsigned char scalar_bytes[32] = {};
        std::memcpy(scalar_bytes, scalars[i].data(), sizeof(scalar_bytes));

        blst::scalar::GT expected = {};
        blst::scalar::Fp12::cyclotomic_exp_basic(
            expected.data, scalar_base.data, scalar_bytes, 256);

        blst::avx::GT actual = {};
        blst::avx::GT::mult(actual, avx_base, scalar_bytes, 256);
        expect_same(expected, actual);

        blst::avx::GT actual_vartime = {};
        blst::avx::GT::mult_vartime(actual_vartime, avx_base, scalar_bytes, 256);
        expect_same(expected, actual_vartime);
    }

    const std::array<unsigned char, 16> scalar_128 = {
        0x6d, 0xb4, 0x23, 0x91, 0xe5, 0x77, 0xa2, 0x0c,
        0x3f, 0x18, 0xc8, 0x44, 0x9a, 0x52, 0x01, 0xd7,
    };
    blst::scalar::GT expected_128 = {};
    blst::scalar::Fp12::cyclotomic_exp_basic(
        expected_128.data, scalar_base.data, scalar_128.data(), 128);
    blst::avx::GT actual_128 = {};
    blst::avx::GT::mult_vartime(
        actual_128, avx_base, scalar_128.data(), 128);
    expect_same(expected_128, actual_128);
}
