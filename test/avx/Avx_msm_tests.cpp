#include <gtest/gtest.h>

#include "blst/avx/curve/multi_scalar.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace {

using ScalarBytes = std::array<blst::byte, 32>;

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

}  // namespace

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
}
