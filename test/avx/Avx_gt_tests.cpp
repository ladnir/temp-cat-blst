#include <gtest/gtest.h>

#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/point.h"
#include "blst/avx/curve/GT.h"

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

    blst::scalar::GT out = {};
    blst::scalar::GT::miller_loop(out, &g2_affine, &g1_affine);
    blst::scalar::GT::final_mult(out, out);
    return out;
}

blst::avx::GT avx_from_scalar(const blst::scalar::GT& value)
{
    unsigned char bytes[12 * 48] = {};
    blst::scalar::GT::serialize(bytes, value);

    blst::avx::GT out = {};
    blst::avx::GT::deserialize(out, bytes);
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
