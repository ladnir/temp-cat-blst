/*
 * Copyright Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <array>
#include <cstring>

#include <gtest/gtest.h>

#include "blst/scalar/curve/point.h"

namespace
{
using blst::BLST_POINT_NOT_IN_GROUP;
using blst::BLST_SUCCESS;
using blst::byte;
using blst::scalar::G1AffinePoint;
using blst::scalar::G1Point;
using blst::scalar::G2AffinePoint;
using blst::scalar::G2Point;

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
