#include <gtest/gtest.h>

#include "blst/avx/curve/multi_scalar.h"
#include "blst/scalar/curve/multi_scalar.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "bindings/blst_avx.h"

namespace {

using ScalarBytes = std::array<blst::byte, 32>;

static_assert(blst::avx::valid_precomputed_table_log_size<1>);
static_assert(blst::avx::valid_precomputed_table_log_size<4>);
static_assert(blst::avx::valid_precomputed_table_log_size<7>);
static_assert(!blst::avx::valid_precomputed_table_log_size<0>);
static_assert(!blst::avx::valid_precomputed_table_log_size<8>);

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

std::vector<blst::byte> scalar_batch(size_t count)
{
    std::vector<blst::byte> scalars(count * 32);
    for (size_t i = 0; i < count; ++i) {
        const ScalarBytes scalar = scalar_bytes(1000 + i);
        std::copy(scalar.begin(), scalar.end(), scalars.begin() + static_cast<std::ptrdiff_t>(i * 32));
    }
    return scalars;
}

std::vector<blst::byte> scalar_batch(size_t count, size_t nbytes)
{
    std::vector<blst::byte> scalars(count * nbytes);
    for (size_t i = 0; i < count; ++i) {
        const ScalarBytes scalar = scalar_bytes(1000 + i);
        std::copy_n(scalar.begin(), nbytes, scalars.begin() + static_cast<std::ptrdiff_t>(i * nbytes));
    }
    return scalars;
}

std::vector<blst::avx::G1AffinePoint> g1_points(size_t count)
{
    std::vector<blst::avx::G1AffinePoint> out(count);
    for (size_t i = 0; i < count; ++i) {
        const ScalarBytes scalar = scalar_bytes(42 + i);
        blst::avx::G1Point point = {};
        blst::avx::G1Point::mult(point, blst::avx::G1Point::generator(), scalar.data(), 256);
        blst::avx::G1Point::to_affine(out[i], point);
    }
    return out;
}

std::vector<blst::avx::G2AffinePoint> g2_points(size_t count)
{
    std::vector<blst::avx::G2AffinePoint> out(count);
    for (size_t i = 0; i < count; ++i) {
        const ScalarBytes scalar = scalar_bytes(42 + i);
        blst::avx::G2Point point = {};
        blst::avx::G2Point::mult(point, blst::avx::G2Point::generator(), scalar.data(), 256);
        blst::avx::G2Point::to_affine(out[i], point);
    }
    return out;
}

std::vector<blst::avx::GT> gt_points(size_t count)
{
    const std::vector<blst::avx::G1AffinePoint> g1 = g1_points(count);
    const std::vector<blst::avx::G2AffinePoint> g2 = g2_points(count);
    std::vector<blst::avx::GT> out(count);
    for (size_t i = 0; i < count; ++i) {
        blst::avx::GT miller = {};
        blst::avx::GT::miller_loop(miller, &g2[i], &g1[i]);
        blst::avx::GT::final_mult(out[i], miller);
    }
    return out;
}

blst::avx::GT gt_naive(
    const std::vector<blst::avx::GT>& points,
    const std::vector<blst::byte>& scalars,
    size_t nbits)
{
    const size_t nbytes = (nbits + 7) / 8;
    blst::avx::GT expected = blst::avx::GT::identity();
    for (size_t i = 0; i < points.size(); ++i) {
        blst::avx::GT term = {};
        blst::avx::GT::mult(term, points[i], scalars.data() + i * nbytes, nbits);
        blst::avx::GT::add(expected, expected, term);
    }
    return expected;
}

blst::scalar::GT scalar_from_avx(const blst::avx::GT& point)
{
    unsigned char bytes[12 * 48] = {};
    blst::avx::GT::serialize(bytes, point);
    blst::scalar::GT out = {};
    blst::scalar::GT::deserialize(out, bytes);
    return out;
}

blst_fp12 public_fp12_from_scalar_gt(const blst::scalar::GT& point)
{
    static_assert(sizeof(blst_fp12) == sizeof(blst::scalar::Fp12));
    blst_fp12 out = {};
    std::memcpy(&out, &point.data, sizeof(out));
    return out;
}

blst::scalar::GT scalar_gt_from_public_fp12(const blst_fp12& point)
{
    static_assert(sizeof(blst_fp12) == sizeof(blst::scalar::Fp12));
    blst::scalar::Fp12 fp12 = {};
    std::memcpy(&fp12, &point, sizeof(fp12));
    return blst::scalar::GT(fp12);
}

}  // namespace

TEST(AvxScalarMult, G1DefaultWindowMatchesScalar)
{
    for (const uint64_t seed : {1ULL, 42ULL, 777ULL}) {
        const ScalarBytes scalar = scalar_bytes(seed);
        blst::avx::G1Point avx = {};
        blst::scalar::G1Point scalar_point = {};
        blst::avx::G1Point::mult(avx, blst::avx::G1Point::generator(), scalar.data(), 256);
        blst::scalar::G1Point::mult(
            scalar_point, blst::scalar::G1Point::generator(), scalar.data(), 256);

        std::array<unsigned char, 96> avx_bytes = {};
        std::array<unsigned char, 96> scalar_serialized = {};
        blst::avx::G1Point::serialize(avx_bytes.data(), avx);
        blst::scalar::G1Point::serialize(scalar_serialized.data(), scalar_point);
        EXPECT_EQ(avx_bytes, scalar_serialized);
    }
}

TEST(AvxScalarMult, G2DefaultWindowMatchesScalar)
{
    for (const uint64_t seed : {1ULL, 42ULL, 777ULL}) {
        const ScalarBytes scalar = scalar_bytes(seed);
        blst::avx::G2Point avx = {};
        blst::scalar::G2Point scalar_point = {};
        blst::avx::G2Point::mult(avx, blst::avx::G2Point::generator(), scalar.data(), 256);
        blst::scalar::G2Point::mult(
            scalar_point, blst::scalar::G2Point::generator(), scalar.data(), 256);

        std::array<unsigned char, 192> avx_bytes = {};
        std::array<unsigned char, 192> scalar_serialized = {};
        blst::avx::G2Point::serialize(avx_bytes.data(), avx);
        blst::scalar::G2Point::serialize(scalar_serialized.data(), scalar_point);
        EXPECT_EQ(avx_bytes, scalar_serialized);
    }
}

TEST(AvxMsm, G1PippengerMatchesNaive)
{
    constexpr size_t kCount = 17;
    const std::vector<blst::avx::G1AffinePoint> points = g1_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount);

    blst::avx::G1Point msm = {};
    blst::avx::avx_g1_mult_pippenger(msm, points.data(), scalars.data(), kCount, 256);

    blst::avx::G1Point expected = blst::avx::G1Point::infinity();
    for (size_t i = 0; i < kCount; ++i) {
        if (blst::avx::G1AffinePoint::is_inf(points[i])) {
            continue;
        }
        blst::avx::G1Point point = {};
        blst::avx::G1Point::from_affine(point, points[i]);
        blst::avx::G1Point::mult(point, point, scalars.data() + i * 32, 256);
        blst::avx::G1Point::dadd(expected, expected, point);
    }

    EXPECT_TRUE(blst::avx::G1Point::equal(msm, expected));
}

TEST(AvxMsm, G1PippengerWithCallerScratchMatchesNaive)
{
    constexpr size_t kCount = 33;
    const std::vector<blst::avx::G1AffinePoint> points = g1_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount);
    const size_t window_size = blst::avx::avx_g1_pippenger_window_size(kCount);
    std::vector<blst::avx::AvxG1PippengerBucket> buckets(
        blst::avx::avx_g1_mult_pippenger_scratch_sizeof(kCount) / sizeof(blst::avx::AvxG1PippengerBucket));
    std::vector<uint8_t> bucket_used(blst::avx::avx_g1_mult_pippenger_bucket_used_scratch_sizeof(kCount));

    blst::avx::G1Point msm = {};
    blst::avx::avx_g1_mult_pippenger(
        msm, points.data(), scalars.data(), kCount, 256, buckets.data(), window_size, bucket_used.data());

    blst::avx::G1Point expected = blst::avx::G1Point::infinity();
    for (size_t i = 0; i < kCount; ++i) {
        if (blst::avx::G1AffinePoint::is_inf(points[i])) {
            continue;
        }
        blst::avx::G1Point point = {};
        blst::avx::G1Point::from_affine(point, points[i]);
        blst::avx::G1Point::mult(point, point, scalars.data() + i * 32, 256);
        blst::avx::G1Point::dadd(expected, expected, point);
    }

    EXPECT_TRUE(blst::avx::G1Point::equal(msm, expected));
}

TEST(AvxMsm, G1PippengerFullWidthDefaultMatchesPlainWindow)
{
    constexpr size_t kCount = 33;
    const std::vector<blst::avx::G1AffinePoint> points = g1_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount);
    std::vector<blst::avx::AvxG1PippengerBucket> buckets(
        blst::avx::avx_g1_mult_pippenger_scratch_sizeof(kCount) / sizeof(blst::avx::AvxG1PippengerBucket));
    std::vector<uint8_t> bucket_used(blst::avx::avx_g1_mult_pippenger_bucket_used_scratch_sizeof(kCount));

    blst::avx::G1Point glv_default = {};
    blst::avx::avx_g1_mult_pippenger(
        glv_default, points.data(), scalars.data(), kCount, 256, buckets.data(), 0, bucket_used.data());

    blst::avx::G1Point plain_explicit = {};
    blst::avx::avx_g1_mult_pippenger(
        plain_explicit,
        points.data(),
        scalars.data(),
        kCount,
        256,
        buckets.data(),
        blst::avx::avx_g1_pippenger_window_size(kCount),
        bucket_used.data());

    EXPECT_TRUE(blst::avx::G1Point::equal(glv_default, plain_explicit));
}

TEST(AvxMsm, G1PippengerShortScalarsUseCompactStride)
{
    constexpr size_t kCount = 17;
    constexpr size_t kBits = 128;
    constexpr size_t kBytes = kBits / 8;
    const std::vector<blst::avx::G1AffinePoint> points = g1_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount, kBytes);

    blst::avx::G1Point msm = {};
    blst::avx::avx_g1_mult_pippenger(msm, points.data(), scalars.data(), kCount, kBits);

    blst::avx::G1Point expected = blst::avx::G1Point::infinity();
    for (size_t i = 0; i < kCount; ++i) {
        blst::avx::G1Point point = {};
        blst::avx::G1Point::from_affine(point, points[i]);
        blst::avx::G1Point::mult(point, point, scalars.data() + i * kBytes, kBits);
        blst::avx::G1Point::dadd(expected, expected, point);
    }

    EXPECT_TRUE(blst::avx::G1Point::equal(msm, expected));
}

TEST(AvxMsm, G1PippengerHonorsAffineIdentitySentinel)
{
    constexpr size_t kCount = 17;
    std::vector<blst::avx::G1AffinePoint> points = g1_points(kCount);
    points[3] = {};
    const std::vector<blst::byte> scalars = scalar_batch(kCount);

    blst::avx::G1Point msm = {};
    blst::avx::avx_g1_mult_pippenger(msm, points.data(), scalars.data(), kCount, 256);

    blst::avx::G1Point expected = blst::avx::G1Point::infinity();
    for (size_t i = 0; i < kCount; ++i) {
        if (blst::avx::G1AffinePoint::is_inf(points[i])) {
            continue;
        }
        blst::avx::G1Point point = {};
        blst::avx::G1Point::from_affine(point, points[i]);
        blst::avx::G1Point::mult(point, point, scalars.data() + i * 32, 256);
        blst::avx::G1Point::dadd(expected, expected, point);
    }

    EXPECT_TRUE(blst::avx::G1Point::equal(msm, expected));
}

TEST(AvxMsm, G2PippengerMatchesNaive)
{
    constexpr size_t kCount = 17;
    const std::vector<blst::avx::G2AffinePoint> points = g2_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount);

    blst::avx::G2Point msm = {};
    blst::avx::avx_g2_mult_pippenger(msm, points.data(), scalars.data(), kCount, 256);

    blst::avx::G2Point expected = blst::avx::G2Point::infinity();
    for (size_t i = 0; i < kCount; ++i) {
        blst::avx::G2Point point = {};
        blst::avx::G2Point::from_affine(point, points[i]);
        blst::avx::G2Point::mult(point, point, scalars.data() + i * 32, 256);
        blst::avx::G2Point::dadd(expected, expected, point);
    }

    EXPECT_TRUE(blst::avx::G2Point::equal(msm, expected));
}

TEST(AvxMsm, G2PippengerHonorsAffineIdentitySentinel)
{
    constexpr size_t kCount = 17;
    std::vector<blst::avx::G2AffinePoint> points = g2_points(kCount);
    points[3] = {};
    const std::vector<blst::byte> scalars = scalar_batch(kCount);

    blst::avx::G2Point msm = {};
    blst::avx::avx_g2_mult_pippenger(msm, points.data(), scalars.data(), kCount, 256);

    blst::avx::G2Point expected = blst::avx::G2Point::infinity();
    for (size_t i = 0; i < kCount; ++i) {
        if (blst::avx::G2AffinePoint::is_inf(points[i])) {
            continue;
        }
        blst::avx::G2Point point = {};
        blst::avx::G2Point::from_affine(point, points[i]);
        blst::avx::G2Point::mult(point, point, scalars.data() + i * 32, 256);
        blst::avx::G2Point::dadd(expected, expected, point);
    }

    EXPECT_TRUE(blst::avx::G2Point::equal(msm, expected));
}

TEST(AvxMsm, G2PippengerFullWidthDefaultMatchesPlainWindow)
{
    constexpr size_t kCount = 33;
    const std::vector<blst::avx::G2AffinePoint> points = g2_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount);
    std::vector<blst::avx::AvxG2PippengerBucket> buckets(
        blst::avx::avx_g2_mult_pippenger_scratch_sizeof(kCount) / sizeof(blst::avx::AvxG2PippengerBucket));
    std::vector<uint8_t> bucket_used(blst::avx::avx_g2_mult_pippenger_bucket_used_scratch_sizeof(kCount));

    blst::avx::G2Point gls_default = {};
    blst::avx::avx_g2_mult_pippenger(
        gls_default, points.data(), scalars.data(), kCount, 256, buckets.data(), 0, bucket_used.data());

    blst::avx::G2Point plain_explicit = {};
    blst::avx::avx_g2_mult_pippenger(
        plain_explicit,
        points.data(),
        scalars.data(),
        kCount,
        256,
        buckets.data(),
        blst::avx::avx_g2_pippenger_window_size(kCount),
        bucket_used.data());

    EXPECT_TRUE(blst::avx::G2Point::equal(gls_default, plain_explicit));
}

TEST(AvxMsm, G2PippengerShortScalarsUseCompactStride)
{
    constexpr size_t kCount = 17;
    constexpr size_t kBits = 128;
    constexpr size_t kBytes = kBits / 8;
    const std::vector<blst::avx::G2AffinePoint> points = g2_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount, kBytes);

    blst::avx::G2Point msm = {};
    blst::avx::avx_g2_mult_pippenger(msm, points.data(), scalars.data(), kCount, kBits);

    blst::avx::G2Point expected = blst::avx::G2Point::infinity();
    for (size_t i = 0; i < kCount; ++i) {
        blst::avx::G2Point point = {};
        blst::avx::G2Point::from_affine(point, points[i]);
        blst::avx::G2Point::mult(point, point, scalars.data() + i * kBytes, kBits);
        blst::avx::G2Point::dadd(expected, expected, point);
    }

    EXPECT_TRUE(blst::avx::G2Point::equal(msm, expected));
}

TEST(AvxMsm, G2PippengerWithCallerScratchMatchesNaive)
{
    constexpr size_t kCount = 33;
    const std::vector<blst::avx::G2AffinePoint> points = g2_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount);
    const size_t window_size = blst::avx::avx_g2_pippenger_window_size(kCount);
    std::vector<blst::avx::AvxG2PippengerBucket> buckets(
        blst::avx::avx_g2_mult_pippenger_scratch_sizeof(kCount) / sizeof(blst::avx::AvxG2PippengerBucket));
    std::vector<uint8_t> bucket_used(blst::avx::avx_g2_mult_pippenger_bucket_used_scratch_sizeof(kCount));

    blst::avx::G2Point msm = {};
    blst::avx::avx_g2_mult_pippenger(
        msm, points.data(), scalars.data(), kCount, 256, buckets.data(), window_size, bucket_used.data());

    blst::avx::G2Point expected = blst::avx::G2Point::infinity();
    for (size_t i = 0; i < kCount; ++i) {
        blst::avx::G2Point point = {};
        blst::avx::G2Point::from_affine(point, points[i]);
        blst::avx::G2Point::mult(point, point, scalars.data() + i * 32, 256);
        blst::avx::G2Point::dadd(expected, expected, point);
    }

    EXPECT_TRUE(blst::avx::G2Point::equal(msm, expected));
}

TEST(AvxMsm, ScratchSizeHelpersRejectOverflowSizedInputs)
{
    constexpr size_t kHuge = std::numeric_limits<size_t>::max();

    EXPECT_LE(blst::avx::avx_g1_pippenger_window_size(kHuge), 16u);
    EXPECT_LE(blst::avx::avx_g2_pippenger_window_size(kHuge), 16u);
    EXPECT_LE(blst::avx::avx_g1_glv_pippenger_window_size(kHuge), 16u);
    EXPECT_LE(blst::avx::avx_g2_gls_pippenger_window_size(kHuge), 16u);

    EXPECT_GE(
        blst::avx::avx_g1_mult_pippenger_scratch_sizeof(256),
        blst::avx::avx_g1_mult_pippenger_glv_scratch_sizeof(256));
    EXPECT_GE(
        blst::avx::avx_g1_mult_pippenger_bucket_used_scratch_sizeof(256),
        blst::avx::avx_g1_mult_pippenger_glv_bucket_used_scratch_sizeof(256));
    EXPECT_GE(
        blst::avx::avx_g2_mult_pippenger_scratch_sizeof(256),
        blst::avx::avx_g2_mult_pippenger_gls_scratch_sizeof(256));
    EXPECT_GE(
        blst::avx::avx_g2_mult_pippenger_bucket_used_scratch_sizeof(256),
        blst::avx::avx_g2_mult_pippenger_gls_bucket_used_scratch_sizeof(256));

    EXPECT_EQ(blst::avx::avx_g1_mult_pippenger_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_g1_mult_pippenger_bucket_used_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_g1_mult_pippenger_glv_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_g1_mult_pippenger_glv_bucket_used_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_g1_mult_pippenger_glv_points_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_g1_mult_pippenger_glv_scalars_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_g2_mult_pippenger_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_g2_mult_pippenger_bucket_used_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_g2_mult_pippenger_gls_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_g2_mult_pippenger_gls_bucket_used_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_g2_mult_pippenger_gls_points_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_g2_mult_pippenger_gls_scalars_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_gt_mult_pippenger_glv4_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_gt_mult_pippenger_glv4_bucket_used_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_gt_mult_pippenger_glv4_points_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_gt_mult_pippenger_glv4_scalars_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_gt_mult_pippenger_scratch_sizeof(kHuge), 0u);
    EXPECT_EQ(blst::avx::avx_gt_mult_pippenger_bucket_used_scratch_sizeof(kHuge), 0u);
}

TEST(AvxMsm, GTDefaultFullWidthGLV4MatchesNaive)
{
    constexpr size_t kCount = 17;
    const std::vector<blst::avx::GT> points = gt_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount);

    blst::avx::GT msm = {};
    blst::avx::avx_gt_mult_pippenger(msm, points.data(), scalars.data(), kCount, 256);

    const blst::avx::GT expected = gt_naive(points, scalars, 256);
    EXPECT_TRUE(blst::avx::GT::equal(msm, expected));
}

TEST(AvxMsm, GTGLV4WithCallerScratchMatchesNaive)
{
    constexpr size_t kCount = 33;
    const std::vector<blst::avx::GT> points = gt_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount);
    std::vector<blst::avx::AvxGtPippengerBucket> buckets(
        blst::avx::avx_gt_mult_pippenger_scratch_sizeof(kCount) / sizeof(blst::avx::AvxGtPippengerBucket));
    std::vector<uint8_t> bucket_used(blst::avx::avx_gt_mult_pippenger_bucket_used_scratch_sizeof(kCount));

    blst::avx::GT msm = {};
    blst::avx::avx_gt_mult_pippenger(
        msm, points.data(), scalars.data(), kCount, 256, buckets.data(), 0, bucket_used.data());

    const blst::avx::GT expected = gt_naive(points, scalars, 256);
    EXPECT_TRUE(blst::avx::GT::equal(msm, expected));
}

TEST(AvxMsm, GTGLV4CallerSplitScratchMatchesAllocated)
{
    constexpr size_t kCount = 33;
    const std::vector<blst::avx::GT> points = gt_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount);
    std::vector<blst::avx::AvxGtPippengerBucket> buckets(
        blst::avx::avx_gt_mult_pippenger_scratch_sizeof(kCount) / sizeof(blst::avx::AvxGtPippengerBucket));
    std::vector<uint8_t> bucket_used(blst::avx::avx_gt_mult_pippenger_bucket_used_scratch_sizeof(kCount));
    std::vector<blst::avx::GT> split_points(
        blst::avx::avx_gt_mult_pippenger_glv4_points_scratch_sizeof(kCount) / sizeof(blst::avx::GT));
    std::vector<blst::byte> split_scalars(blst::avx::avx_gt_mult_pippenger_glv4_scalars_scratch_sizeof(kCount));

    blst::avx::GT allocated = {};
    blst::avx::avx_gt_mult_pippenger_glv4(allocated, points.data(), scalars.data(), kCount);

    blst::avx::GT caller_scratch = {};
    blst::avx::avx_gt_mult_pippenger_glv4(
        caller_scratch,
        points.data(),
        scalars.data(),
        kCount,
        buckets.data(),
        0,
        bucket_used.data(),
        split_points.data(),
        split_scalars.data());

    EXPECT_TRUE(blst::avx::GT::equal(caller_scratch, allocated));
}

TEST(AvxMsm, GTIdentityInputsAndZeroScalarsDoNotPoisonBuckets)
{
    constexpr size_t kCount = 17;
    std::vector<blst::avx::GT> points = gt_points(kCount);
    std::vector<blst::byte> scalars = scalar_batch(kCount);
    points[3] = blst::avx::GT::identity();
    std::fill_n(scalars.begin() + static_cast<std::ptrdiff_t>(5 * 32), 32, 0);

    blst::avx::GT msm = {};
    blst::avx::avx_gt_mult_pippenger(msm, points.data(), scalars.data(), kCount, 256);

    const blst::avx::GT expected = gt_naive(points, scalars, 256);
    EXPECT_TRUE(blst::avx::GT::equal(msm, expected));
}

TEST(AvxMsm, GTCApiPippengerMatchesCppGLV4)
{
    constexpr size_t kCount = 17;
    const std::vector<blst::avx::GT> points = gt_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount);
    std::vector<blst_avx_fp12> c_points(kCount);
    for (size_t i = 0; i < kCount; ++i) {
        const blst_fp12 public_point = public_fp12_from_scalar_gt(scalar_from_avx(points[i]));
        ASSERT_TRUE(blst_avx_fp12_from_scalar(&c_points[i], &public_point));
    }

    const size_t scratch_bytes = blst_avx_fp12_mult_pippenger_scratch_sizeof(kCount);
    const size_t bucket_used_bytes = blst_avx_fp12_mult_pippenger_bucket_used_scratch_sizeof(kCount);
    ASSERT_GT(scratch_bytes, 0u);
    ASSERT_GT(bucket_used_bytes, 0u);
    std::vector<blst_avx_fp12> scratch((scratch_bytes + sizeof(blst_avx_fp12) - 1) / sizeof(blst_avx_fp12));
    std::vector<blst::byte> bucket_used(bucket_used_bytes);

    blst_avx_fp12 c_msm = {};
    ASSERT_TRUE(blst_avx_fp12_mult_pippenger(
        &c_msm,
        c_points.data(),
        scalars.data(),
        kCount,
        scratch.data(),
        scratch.size() * sizeof(blst_avx_fp12),
        bucket_used.data(),
        bucket_used.size()));

    blst_fp12 c_public = {};
    ASSERT_TRUE(blst_fp12_from_avx(&c_public, &c_msm));
    blst::avx::GT expected_avx = {};
    blst::avx::avx_gt_mult_pippenger(expected_avx, points.data(), scalars.data(), kCount, 256);
    const blst::scalar::GT expected = scalar_from_avx(expected_avx);
    const blst::scalar::GT actual = scalar_gt_from_public_fp12(c_public);
    EXPECT_TRUE(blst::scalar::GT::equal(actual, expected));
}

TEST(AvxMsm, ScalarGenericGTPippengerMatchesNaive)
{
    constexpr size_t kCount = 17;
    const std::vector<blst::avx::GT> avx_points = gt_points(kCount);
    const std::vector<blst::byte> scalars = scalar_batch(kCount);
    std::vector<blst::scalar::GT> points(kCount);
    for (size_t i = 0; i < kCount; ++i) {
        points[i] = scalar_from_avx(avx_points[i]);
    }

    const blst::scalar::GT* point_ptrs[] = {points.data(), nullptr};
    const blst::byte* scalar_ptrs[] = {scalars.data(), nullptr};
    std::vector<blst::scalar::GT> buckets(
        blst::scalar::msm_mult_pippenger_scratch_sizeof_<blst::scalar::GT>(kCount) / sizeof(blst::scalar::GT));

    blst::scalar::GT msm = {};
    blst::scalar::msm_mult_pippenger<blst::scalar::GT, blst::scalar::GT, blst::scalar::GT>(
        msm, point_ptrs, kCount, scalar_ptrs, 256, buckets.data(), blst::scalar::pippenger_window_size(kCount));

    blst::scalar::GT expected = blst::scalar::GT::identity();
    for (size_t i = 0; i < kCount; ++i) {
        blst::scalar::GT term = {};
        blst::scalar::GT::mult(term, points[i], scalars.data() + i * 32, 256);
        blst::scalar::GT::add(expected, expected, term);
    }

    EXPECT_TRUE(blst::scalar::GT::equal(msm, expected));
}
