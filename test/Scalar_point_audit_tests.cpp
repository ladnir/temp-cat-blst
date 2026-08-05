/*
 * Copyright Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <array>
#include <cstring>
#include <vector>

#include <gtest/gtest.h>

#include "blst/scalar/curve/point.h"
#include "blst/scalar/curve/bulk_addition.h"
#include "bindings/blst_scalar.h"

namespace
{
using blst::BLST_POINT_NOT_IN_GROUP;
using blst::BLST_SUCCESS;
using blst::byte;
using blst::scalar::G1AffinePoint;
using blst::scalar::G1Point;
using blst::scalar::G2AffinePoint;
using blst::scalar::G2Point;

struct FailingBulkAllocator
{
    static void* allocate(size_t) noexcept { return nullptr; }
    static void deallocate(void*) noexcept {}
};

template<typename Point>
using AffineFor = blst::scalar::AffinePoint<typename Point::Field_type>;

template<typename Point>
std::vector<AffineFor<Point>> repeated_generator(size_t count)
{
    std::vector<AffineFor<Point>> points(count);
    for (auto& point : points)
        Point::to_affine(point, Point::generator());
    return points;
}

template<typename Point>
Point serial_sum(const std::vector<AffineFor<Point>>& points)
{
    Point sum = Point::infinity();
    for (const auto& point : points)
        Point::dadd(sum, sum, point);
    return sum;
}

template<typename Point>
void expect_compact_batch_sum(size_t count)
{
    using Field = typename Point::Field_type;
    auto affine = repeated_generator<Point>(count);
    const AffineFor<Point>* compact[] = {affine.data(), nullptr};
    Point actual = {};
    blst::scalar::points_add<Field>(actual, compact, affine.size());
    EXPECT_TRUE(Point::equal(actual, serial_sum<Point>(affine)));
}

template<typename Point>
void expect_bulk_boundaries()
{
    using Field = typename Point::Field_type;
    using Node = blst::scalar::detail::BatchAffineNode<Field>;
    constexpr size_t capacity = BLST_SCALAR_BULK_SCRATCH_LIMIT / sizeof(Node);

    expect_compact_batch_sum<Point>(0);
    expect_compact_batch_sum<Point>(1);
    expect_compact_batch_sum<Point>(15);
    expect_compact_batch_sum<Point>(16);
    expect_compact_batch_sum<Point>(17);
    expect_compact_batch_sum<Point>(capacity + 1);
}

template<typename Point>
void expect_special_and_allocation_failure()
{
    using Field = typename Point::Field_type;
    using Affine = AffineFor<Point>;

    Affine generator;
    Point::to_affine(generator, Point::generator());
    Affine opposite = generator;
    Affine::cneg(opposite, 1);
    const Affine infinity = {};

    std::vector<Affine> affine(17, generator);
    affine[1] = generator; // equal-point case
    affine[2] = opposite;  // opposite-point case
    affine[3] = infinity;

    std::vector<const Affine*> explicit_points;
    explicit_points.reserve(affine.size());
    for (const auto& point : affine)
        explicit_points.push_back(&point);

    const Point expected = serial_sum<Point>(affine);
    Point batched = {};
    blst::scalar::points_add<Field>(
        batched, explicit_points.data(), explicit_points.size());
    EXPECT_TRUE(Point::equal(batched, expected));

    const Affine* compact[] = {affine.data(), nullptr};
    Point fallback = {};
    blst::scalar::points_add<Field, FailingBulkAllocator>(
        fallback, compact, affine.size());
    EXPECT_TRUE(Point::equal(fallback, expected));
}

bool fp2_equal(const blst::scalar::Fp2& a, const blst::scalar::Fp2& b)
{
    return blst::vec_is_equal(&a, &b, sizeof(a)) != 0;
}

TEST(ScalarPointAudit, TypedAffineGeneratorsMatchProjectiveGenerators)
{
    const G1Point& g1 = G1Point::generator();
    const G1AffinePoint& g1_affine = G1AffinePoint::generator();
    EXPECT_TRUE(blst::scalar::Fp::equal(g1.x, g1_affine.x));
    EXPECT_TRUE(blst::scalar::Fp::equal(g1.y, g1_affine.y));
    EXPECT_TRUE(blst::scalar::Fp::equal(g1.z, blst::scalar::Fp::one()));

    const G2Point& g2 = G2Point::generator();
    const G2AffinePoint& g2_affine = G2AffinePoint::generator();
    EXPECT_TRUE(fp2_equal(g2.x, g2_affine.x));
    EXPECT_TRUE(fp2_equal(g2.y, g2_affine.y));
    EXPECT_TRUE(fp2_equal(g2.z, blst::scalar::Fp2::one()));
}

TEST(ScalarPointAudit, ProjectiveSerializationUsesTypedCoordinates)
{
    std::array<byte, 192> projective = {};
    std::array<byte, 192> affine = {};

    G1Point::serialize(projective.data(), G1Point::generator());
    G1AffinePoint::serialize(affine.data(), G1AffinePoint::generator());
    EXPECT_EQ(0, std::memcmp(projective.data(), affine.data(), 96));

    G2Point::serialize(projective.data(), G2Point::generator());
    G2AffinePoint::serialize(affine.data(), G2AffinePoint::generator());
    EXPECT_EQ(0, std::memcmp(projective.data(), affine.data(), 192));

    G1Point::compress(projective.data(), G1Point::generator());
    G1AffinePoint::compress(affine.data(), G1AffinePoint::generator());
    EXPECT_EQ(0, std::memcmp(projective.data(), affine.data(), 48));

    G2Point::compress(projective.data(), G2Point::generator());
    G2AffinePoint::compress(affine.data(), G2AffinePoint::generator());
    EXPECT_EQ(0, std::memcmp(projective.data(), affine.data(), 96));
}

TEST(ScalarPointAudit, G2MixedAdditionSelectsTheFp2Identity)
{
    const G2AffinePoint& generator_affine = G2AffinePoint::generator();
    const G2Point& generator = G2Point::generator();

    G2Point result = {};
    G2Point::add(result, G2Point::infinity(), generator_affine);
    EXPECT_TRUE(G2Point::equal(result, generator));

    G2Point aliased = G2Point::infinity();
    G2Point::add(aliased, aliased, generator_affine);
    EXPECT_TRUE(G2Point::equal(aliased, generator));

    const G2AffinePoint affine_infinity = {};
    G2Point::add(result, generator, affine_infinity);
    EXPECT_TRUE(G2Point::equal(result, generator));
}

TEST(ScalarPointAudit, CompleteAdditionHandlesEqualAndOppositePoints)
{
    G1Point doubled = {};
    G1Point complete = {};
    G1Point::doubling(doubled, G1Point::generator());
    G1Point::dadd(complete, G1Point::generator(), G1Point::generator());
    EXPECT_TRUE(G1Point::equal(complete, doubled));

    G1Point opposite = G1Point::generator();
    G1Point::cneg(opposite, 1);
    G1Point::dadd(complete, G1Point::generator(), opposite);
    EXPECT_TRUE(G1Point::is_inf(complete));
}

TEST(ScalarPointAudit, G2CompleteAdditionHandlesEqualAndOppositePoints)
{
    G2Point doubled = {};
    G2Point complete = {};
    G2Point::doubling(doubled, G2Point::generator());
    G2Point::dadd(complete, G2Point::generator(), G2Point::generator());
    EXPECT_TRUE(G2Point::equal(complete, doubled));

    G2Point opposite = G2Point::generator();
    G2Point::cneg(opposite, 1);
    G2Point::dadd(complete, G2Point::generator(), opposite);
    EXPECT_TRUE(G2Point::is_inf(complete));
}

TEST(ScalarPointAudit, BatchAffineScratchPreservesTheProjectiveFootprint)
{
    using G1Node = blst::scalar::detail::BatchAffineNode<blst::scalar::Fp>;
    using G2Node = blst::scalar::detail::BatchAffineNode<blst::scalar::Fp2>;
    static_assert(sizeof(G1Node) == sizeof(G1Point));
    static_assert(alignof(G1Node) == alignof(G1Point));
    static_assert(sizeof(G2Node) == sizeof(G2Point));
    static_assert(alignof(G2Node) == alignof(G2Point));
    SUCCEED();
}

TEST(ScalarPointAudit, G1BulkAdditionCoversBatchAndChunkBoundaries)
{
    expect_bulk_boundaries<G1Point>();
    expect_special_and_allocation_failure<G1Point>();
}

TEST(ScalarPointAudit, G2BulkAdditionCoversBatchAndChunkBoundaries)
{
    expect_bulk_boundaries<G2Point>();
    expect_special_and_allocation_failure<G2Point>();
}

TEST(ScalarPointAudit, CBatchAdditionUsesCompactFormAndRejectsNullHead)
{
    auto g1 = repeated_generator<G1Point>(17);
    const blst_p1_affine* g1_compact[] = {g1.data(), nullptr};
    blst_p1 g1_sum = {};
    blst_p1s_add(&g1_sum, g1_compact, g1.size());
    EXPECT_TRUE(G1Point::equal(g1_sum, serial_sum<G1Point>(g1)));

    const blst_p1_affine* invalid[] = {nullptr};
    blst_p1s_add(&g1_sum, invalid, 1);
    EXPECT_TRUE(G1Point::is_inf(g1_sum));

    auto g2 = repeated_generator<G2Point>(17);
    const blst_p2_affine* g2_compact[] = {g2.data(), nullptr};
    blst_p2 g2_sum = {};
    blst_p2s_add(&g2_sum, g2_compact, g2.size());
    EXPECT_TRUE(G2Point::equal(g2_sum, serial_sum<G2Point>(g2)));

    const blst_p2_affine* invalid_g2[] = {nullptr};
    blst_p2s_add(&g2_sum, invalid_g2, 1);
    EXPECT_TRUE(G2Point::is_inf(g2_sum));
}

TEST(ScalarPointAudit, CFixedTableMsmConstructsTypedScratchInRawStorage)
{
    constexpr size_t count = 17;
    constexpr size_t wbits = 4;
    auto affine = repeated_generator<G1Point>(count);
    const blst_p1_affine* compact_points[] = {affine.data(), nullptr};

    std::vector<blst_p1_affine> table(
        blst_p1s_mult_wbits_precompute_sizeof(wbits, count) /
        sizeof(blst_p1_affine));
    blst_p1s_mult_wbits_precompute(
        table.data(), wbits, compact_points, count);

    std::array<byte, count> scalar_bytes = {};
    G1Point expected = G1Point::infinity();
    for (size_t i = 0; i < count; ++i) {
        scalar_bytes[i] = static_cast<byte>(i + 1);
        G1Point product = {};
        G1Point::mult(product, G1Point::generator(), &scalar_bytes[i], 7);
        G1Point::dadd(expected, expected, product);
    }
    const byte* compact_scalars[] = {scalar_bytes.data(), nullptr};

    const size_t scratch_bytes = blst_p1s_mult_wbits_scratch_sizeof(count);
    ASSERT_EQ(0u, scratch_bytes % sizeof(blst_limb_t));
    std::vector<blst_limb_t> raw_scratch(
        scratch_bytes / sizeof(blst_limb_t), ~blst_limb_t{0});
    blst_p1 actual = {};
    blst_p1s_mult_wbits(
        &actual, table.data(), wbits, count, compact_scalars, 7,
        raw_scratch.data());
    EXPECT_TRUE(G1Point::equal(actual, expected));
}

TEST(ScalarPointAudit, FailedNativeG1DecodingLeavesOutputUnchanged)
{
    const G1AffinePoint expected = G1AffinePoint::generator();
    G1AffinePoint out = expected;

    std::array<byte, 48> compressed_x_zero = {};
    compressed_x_zero[0] = 0x80;
    EXPECT_EQ(
        BLST_POINT_NOT_IN_GROUP,
        G1AffinePoint::uncompress(out, compressed_x_zero.data()));
    EXPECT_TRUE(G1AffinePoint::equal(out, expected));

    std::array<byte, 96> uncompressed_x_zero_y_two = {};
    uncompressed_x_zero_y_two.back() = 2;
    EXPECT_EQ(
        BLST_POINT_NOT_IN_GROUP,
        G1AffinePoint::deserialize(out, uncompressed_x_zero_y_two.data()));
    EXPECT_TRUE(G1AffinePoint::equal(out, expected));

    std::array<byte, 48> valid = {};
    G1AffinePoint::compress(valid.data(), expected);
    out = {};
    ASSERT_EQ(BLST_SUCCESS, G1AffinePoint::uncompress(out, valid.data()));
    EXPECT_TRUE(G1AffinePoint::equal(out, expected));
}
}
