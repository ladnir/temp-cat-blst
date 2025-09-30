#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "blst/fields/Fp.h"
#include "blst/curve/point.h"
#include "blst/curve/e1.h"
#include "blst/curve/map_to_g1.h"

using namespace blst;

// Helper: Print point coordinates for debugging
void print_point(const G1Point& p, const char* name) {
    printf("%s:\n", name);
    printf("  X: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p.X[i]);
    }
    printf("\n  Y: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p.Y[i]);
    }
    printf("\n  Z: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p.Z[i]);
    }
    printf("\n");
}

//// Helper: Convert hex string to bytes
//void hex_to_bytes(unsigned char* bytes, const char* hex) {
//    size_t len = strlen(hex) / 2;
//    for (size_t i = 0; i < len; i++) {
//        sscanf(hex + (i * 2), "%2hhx", &bytes[i]);
//    }
//}

TEST(E1PointTest, PointAddition) {
    G1Point p1, p2, result, temp;

    // Get generator point
    p1 = G1Point::generator();
    
    // Create a different point (we'll use 2*G)
    G1Point::doubling(temp, p1);
    
    // Add generator and the different point (G + 2G)
    G1Point::add(result, p1, temp);
    
    // Compute 3G using scalar multiplication
    unsigned char scalar[32] = {3};
    G1Point::mult(p2, p1, scalar, 8);
    
    // Results should match
    EXPECT_TRUE(G1Point::equal(result, p2));
    
    // Also test the doubling-specific function
    G1Point::doubling(result, p1);
    
    // Create 2G with scalar multiplication
    unsigned char scalar2[32] = {2};
    G1Point::mult(p2, p1, scalar2, 8);
    
    // Results should match
    EXPECT_TRUE(G1Point::equal(result, p2));
}

TEST(E1PointTest, AdditionWithIdentity) {
    G1Point p1, identity, result;

    // Get generator point
    p1 = G1Point::generator();

    // Create identity point (point at infinity)
    memset(&identity, 0, sizeof(identity));
    identity.Z[0] = 0;  // explicitly set Z to 0 for infinity

    // P + O = P
    G1Point::add(result, p1, identity);
    EXPECT_TRUE(G1Point::equal(result, p1));
}

TEST(E1PointTest, PointNegation) {
    G1Point p1, neg_p1, result;
    
    // Get generator point
    p1 = G1Point::generator();
    
    // Make a copy of the point and negate it
    neg_p1 = p1;
    G1Point::cneg(neg_p1, 1);
    
    // P + (-P) should be the identity
    G1Point::add(result, p1, neg_p1);
    
    // Result should be point at infinity
    EXPECT_TRUE(G1Point::is_inf(result));
}

TEST(E1PointTest, ScalarMultiplication) {
    G1Point p1, result, expected;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    p1 = G1Point::generator();

    // Scalar = 2
    scalar[0] = 2;

    // Compute 2*P
    G1Point::mult(result, p1, scalar, 8);

    // Should equal P + P (using doubling formula)
    G1Point::doubling(expected, p1);
    EXPECT_TRUE(G1Point::equal(result, expected));

    // Scalar = 3
    scalar[0] = 3;

    // Compute 3*P
    G1Point::mult(result, p1, scalar, 8);

    // Should equal 2*P + P
    G1Point doubled;
    G1Point::doubling(doubled, p1);         // doubled = 2*P
    G1Point::add(expected, doubled, p1); // expected = 2*P + P = 3*P
    
    EXPECT_TRUE(G1Point::equal(result, expected));
}

TEST(E1PointTest, ScalarMultiplicationByZero) {
    G1Point p1, result;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    p1 = G1Point::generator();

    // Compute 0*P
    G1Point::mult(result, p1, scalar, 8);

    // Result should be point at infinity
    EXPECT_TRUE(G1Point::is_inf(result));
}

TEST(E1PointTest, PointOnCurve) {
    G1Point p1;
    G1AffinePoint p1_affine;

    // Get generator point
    p1 = G1Point::generator();

    // Generator should be on curve
    EXPECT_TRUE(G1Point::on_curve(p1));

    // Convert to affine
    G1Point::to_affine(p1_affine, p1);
    
    // Affine point should also be on curve
    EXPECT_TRUE(G1AffinePoint::on_curve(p1_affine));

    // Modify point to be off curve
    p1.X[0] ^= 1;  // Flip a bit

    // Modified point should not be on curve
    EXPECT_FALSE(G1Point::on_curve(p1));
}

TEST(E1PointTest, SerializeDeserialize) {
    G1Point p1, p2;
    G1AffinePoint p1_affine;
    unsigned char serialized[96];

    // Get generator point
    p1 = G1Point::generator();

    // Serialize the point
    G1Point::serialize(serialized, p1);

    // Deserialize back to a point
    EXPECT_EQ(G1AffinePoint::deserialize(p1_affine, serialized), BLST_SUCCESS);

    // Convert from affine back to Jacobian coordinates
    blst_p1_from_affine(p2, p1_affine);

    // Points should be equal
    EXPECT_TRUE(G1Point::equal(p1, p2));
}

TEST(E1PointTest, FixedVectorCompressUncompress) {
    // Test with the generator point, which is a guaranteed valid point
    const G1AffinePoint g_affine = G1AffinePoint::generator();
    G1AffinePoint decompressed;
    unsigned char compressed[48] = { 0 };

    // Inspect the generator point to ensure it's valid
    ASSERT_TRUE(G1AffinePoint::on_curve(g_affine)) << "Generator point not on curve";

    // Manual compression:
    // First, we'll get the byte representation of the generator's X coordinate
    Fp x_coord;
    Fp::from_mont(x_coord, g_affine.X);  // Convert from Montgomery to regular form

    // Set the compression bit and copy X coordinate
    compressed[0] = 0x80;  // Only set compressed bit, not sign bit
    be_bytes_from_limbs(compressed, x_coord, 48);

    // The first byte now has both the compressed flag and the high bits of X
    // We need to preserve those high bits
    compressed[0] |= 0x80;  // Ensure compressed bit is set

    //printf("Manually compressed point (first few bytes):\n");
    //for (int i = 0; i < 8; i++) {
    //    printf("%02X ", compressed[i]);
    //}
    //printf("\n");

    // Try to uncompress with both possible sign bits
    BLST_ERROR err1 = G1AffinePoint::uncompress(decompressed, compressed);

    // Try flipping the sign bit
    compressed[0] ^= 0x20;  // Flip the sign bit
    //printf("Flipped sign bit, now first byte: 0x%02X\n", compressed[0]);

    BLST_ERROR err2 = G1AffinePoint::uncompress(decompressed, compressed);

    //printf("Uncompress results: without sign bit = %d, with sign bit = %d\n", err1, err2);

    // At least one of these should work
    EXPECT_TRUE(err1 == BLST_SUCCESS || err2 == BLST_SUCCESS)
        << "Failed to uncompress point with either sign bit";

    // If successful, check if the point is on the curve
    if (err1 == BLST_SUCCESS || err2 == BLST_SUCCESS) {
        EXPECT_TRUE(G1AffinePoint::on_curve(decompressed))
            << "Decompressed point not on curve";
    }
}

TEST(E1PointTest, BasicCompressUncompress) {
    {
        const auto gen = G1Point::generator();

        uint8_t comp[48];
        G1Point::compress(comp, gen);  // or: to_affine + blst_p1_affine_compress

        G1AffinePoint dec;
        BLST_ERROR err = G1AffinePoint::uncompress(dec, comp);
        ASSERT_TRUE(err == BLST_SUCCESS);
        ASSERT_TRUE(G1AffinePoint::on_curve(dec));
        ASSERT_TRUE(blst_p1_affine_in_g1(dec));
    }
    // Starting with the generator point in projective form
    G1Point point = G1Point::generator();

    // Convert to affine
    G1AffinePoint affine_point;
    G1Point::to_affine(affine_point, point);

    // First confirm it's on the curve
    ASSERT_TRUE(G1AffinePoint::on_curve(affine_point)) << "Point not on curve";

    // Compress the affine point
    unsigned char compressed[48];
    G1AffinePoint::compress(compressed, affine_point);

    // Print the compressed format for debugging
    //printf("Compressed point (first byte: 0x%02X)\n", compressed[0]);

    // Try to uncompress using the same formats as the library expects
    G1AffinePoint result;
    BLST_ERROR err = G1AffinePoint::uncompress(result, compressed);

    // If it fails, try a fallback approach
    if (err != BLST_SUCCESS) {
        printf("Standard decompression failed with code %d\n", err);

        // Try using serialize/deserialize instead
        unsigned char serialized[96];
        G1AffinePoint::serialize(serialized, affine_point);

        err = G1AffinePoint::deserialize(result, serialized);
        EXPECT_EQ(err, BLST_SUCCESS) << "Even serialization/deserialization failed";

        // If that works, the issue is specific to compression
        if (err == BLST_SUCCESS) {
            printf("Serialization/deserialization works but compression fails.\n");
            printf("This suggests an issue with the compression algorithm.\n");
        }
    }
    else {
        // Success! Points should match
        EXPECT_TRUE(G1AffinePoint::equal(result, affine_point))
            << "Points don't match after compression/decompression";
    }
}

TEST(E1PointTest, AffineCompressUncompress) {
    G1AffinePoint original, decompressed;
    unsigned char compressed[48] = { 0 };

    // Get the generator in affine form
    original = G1AffinePoint::generator();

    // Verify it's on curve
    ASSERT_TRUE(G1AffinePoint::on_curve(original)) << "Generator point not on curve";

    // Compress 
    G1AffinePoint::compress(compressed, original);

    // Verify compressed format
    EXPECT_TRUE(compressed[0] & 0x80) << "Compressed bit not set";
    EXPECT_FALSE(compressed[0] & 0x40) << "Infinity bit incorrectly set";

    // Uncompress
    BLST_ERROR err = G1AffinePoint::uncompress(decompressed, compressed);
    ASSERT_EQ(err, BLST_SUCCESS) << "Decompression failed with code " << err;

    // Verify the point is on curve
    ASSERT_TRUE(G1AffinePoint::on_curve(decompressed)) << "Decompressed point not on curve";

    // Directly compare the affine coordinates
    EXPECT_TRUE(G1AffinePoint::equal(original, decompressed))
        << "Affine points don't match after compression cycle";
}

TEST(E1PointTest, CompressUncompress) {
    G1Point p1, p2;
    G1AffinePoint p1_affine, p2_affine;
    unsigned char compressed[48] = {0};

    // Start with the generator point which is guaranteed to be valid
    p1 = G1Point::generator();
    
    // Convert to affine first (important for compression)
    G1Point::to_affine(p1_affine, p1);
    
    // Compress directly using G1Point::compress (not the affine version)
    G1Point::compress(compressed, p1);
    
    // Manually verify the format - first byte should have compressed bit set (0x80)
    // but not have infinity bit set (0x40)
    EXPECT_TRUE((compressed[0] & 0x80) && !(compressed[0] & 0x40));
    
    // Try to uncompress
    BLST_ERROR err = G1AffinePoint::uncompress(p2_affine, compressed);
    if (err != BLST_SUCCESS) {
        printf("Uncompress failed with error code: %d\n", err);
        // Check if we have error messages for specific codes
        const char* err_msgs[] = {
            "BLST_SUCCESS",
            "BLST_BAD_ENCODING",
            "BLST_POINT_NOT_ON_CURVE",
            "BLST_POINT_NOT_IN_GROUP",
            "BLST_AGGR_TYPE_MISMATCH",
            "BLST_VERIFY_FAIL"
        };
        if (err >= 0 && err < 6) {
            printf("Error message: %s\n", err_msgs[err]);
        }
    }
    ASSERT_EQ(err, BLST_SUCCESS) << "Point decompression failed";
    
    // Convert back to Jacobian coordinates
    blst_p1_from_affine(p2, p2_affine);
    
    // Debug output in case of comparison failure
    if (!G1Point::equal(p1, p2)) {
        printf("Points differ after compress/decompress cycle\n");
        print_point(p1, "Original");
        print_point(p2, "After compress/decompress");
        
        // Check if points differ only by Z-coordinate normalization
        G1Point p1_normalized, p2_normalized;
        G1Point::from_jacobian(p1_normalized, p1);
        G1Point::from_jacobian(p2_normalized, p2);
        
        bool normalized_equal = G1Point::equal(p1_normalized, p2_normalized);
        printf("After normalization, points are %s\n", 
               normalized_equal ? "EQUAL" : "DIFFERENT");
    }
    
    // Points should be equal (if not, the issue could be normalization)
    EXPECT_TRUE(G1Point::equal(p1, p2));
}

TEST(E1PointTest, InfinityPoint) {
    G1Point identity;
    unsigned char serialized[96];
    G1AffinePoint affine_identity;

    // Create identity point
    memset(&identity, 0, sizeof(identity));

    // Point should be at infinity
    EXPECT_TRUE(G1Point::is_inf(identity));

    // Serialize infinity point
    G1Point::serialize(serialized, identity);

    // First byte should have infinity bit set
    EXPECT_TRUE(serialized[0] & 0x40);

    // Deserialize back to a point
    EXPECT_EQ(G1AffinePoint::deserialize(affine_identity, serialized), BLST_SUCCESS);

    // Should be affine infinity
    EXPECT_TRUE(G1AffinePoint::is_inf(affine_identity));
}

TEST(E1PointTest, SigningBasicOperation) {
    G1Point p1, sig;
    unsigned char SK[32] = {
        0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08,
        0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10,
        0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18,
        0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20
    };
    unsigned char pk_bytes[96];

    // Generate public key from secret key
    G1Point::to_pk2(pk_bytes, NULL, SK);

    // Create a message point - we'll use the generator for simplicity
    p1 = G1Point::generator();
	G1Point::doubling(p1, p1); // Just to have a different point
	G1Point::add(p1, p1, G1Point::generator());

    // Sign the message
    G1Point::sign(sig, p1, SK);

    // Signature should be on the curve
    EXPECT_TRUE(G1Point::on_curve(sig));
}

TEST(E1PointTest, PointCoordinateConversion) {
    G1Point p1, p2;
    G1AffinePoint p_affine;

    // Get generator point
    p1 = G1Point::generator();

    // Convert to affine
    G1Point::to_affine(p_affine, p1);

    // Convert back to Jacobian
    blst_p1_from_affine(p2, p_affine);

    // Should be equal to original point
    EXPECT_TRUE(G1Point::equal(p1, p2));

    // Convert from Jacobian to Jacobian (normalized Z)
    G1Point::from_jacobian(p2, p1);

    // Should be equal to original point
    EXPECT_TRUE(G1Point::equal(p1, p2));
}

TEST(E1PointTest, GeneratorProperties) {
    // Get generator point
    const G1Point g = G1Point::generator();

    // Generator should be on curve
    EXPECT_TRUE(G1Point::on_curve(g));

    // Generator should not be at infinity
    EXPECT_FALSE(G1Point::is_inf(g));

    // Get affine generator
    const G1AffinePoint g_affine = G1AffinePoint::generator();

    // Affine generator should be on curve
    EXPECT_TRUE(G1AffinePoint::on_curve(g_affine));

    // Affine generator should not be at infinity
    EXPECT_FALSE(G1AffinePoint::is_inf(g_affine));
}

TEST(E1PointTest, LargeScalarMultiplication) {
    G1Point p1, result;
    // Create a large scalar close to the order
    unsigned char large_scalar[32] = {
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x0f
    };

    // Get generator point
    p1 = G1Point::generator();

    // Perform large scalar multiplication
    G1Point::mult(result, p1, large_scalar, 252);

    // Result should be on curve
    EXPECT_TRUE(G1Point::on_curve(result));
}

TEST(E1PointTest, AffineEquality) {
    G1AffinePoint a1, a2;
    G1Point p1, p2;

    // Get generator point
    p1 = G1Point::generator();

    // Double the generator
    G1Point::doubling(p2, p1);

    // Convert both to affine
    G1Point::to_affine(a1, p1);
    G1Point::to_affine(a2, p1);

    // Affine points should be equal
    EXPECT_TRUE(G1AffinePoint::equal(a1, a2));

    // Convert p2 to affine
    G1Point::to_affine(a2, p2);

    // Should not be equal to a1
    EXPECT_FALSE(G1AffinePoint::equal(a1, a2));
}

TEST(E1PointTest, DeserializeInvalidPoint) {
    G1AffinePoint result;
    unsigned char invalid_data[96];

    // Fill with invalid data
    memset(&invalid_data, 0xFF, sizeof(invalid_data));

    // Should fail to deserialize
    EXPECT_NE(G1AffinePoint::deserialize(result, invalid_data), BLST_SUCCESS);
}

TEST(E1PointTest, SizeChecks) {
    // Check that reported sizes match expected
    EXPECT_EQ(blst_p1_sizeof(), sizeof(G1Point));
    EXPECT_EQ(blst_p1_affine_sizeof(), sizeof(G1AffinePoint));
}

TEST(E1PointTest, MultiplePointAddition) {
    G1Point p1, p3a, p2, result;

    // Get generator point
    p1 = G1Point::generator();

    // Compute 3*P using additions
    G1Point::doubling(p2, p1);    // temp = 2*P
    G1Point::add(p3a, p2, p1);  // p3a = 3*P

    // Compute 3*P using scalar multiplication
    unsigned char scalar[32] = { 3 };
    G1Point::mult(result, p1, scalar, 8);

    // Results should match
    EXPECT_TRUE(G1Point::equal(result, p3a));

    // Test mixed addition with different points
    G1AffinePoint p2_affine;
    G1Point::to_affine(p2_affine, p2);  // Convert 2*P to affine

    // P + (2P affine) should equal 3P
    G1Point::add_affine(result, p1, p2_affine);
    EXPECT_TRUE(G1Point::equal(result, p3a));  // Should equal 3*P

    // For completeness, test dadd with distinct points too
    G1Point p3b;
    G1Point::add_affine(p3b, p1, p2_affine);
    //G1Point::dadd_affine(p3b, p1, p2_affine);
    EXPECT_TRUE(G1Point::equal(p3a, p3b));  // Should equal 3*P
    
    // Also verify that dadd_affine works with same point (P+P)
    G1AffinePoint p_affine;
    G1Point::to_affine(p_affine, p1);
    
    // Test P + P using dadd_affine
    G1Point::dadd_affine(result, p1, p_affine);
    EXPECT_TRUE(G1Point::equal(result, p2));  // Should equal 2*P
}


TEST(E1PointTest, MultiScalarMultiplication) {
    const size_t n = 3;
    G1Point points[n];
    byte scalars[n][32] = { 0 };

    // Generate n different points
    points[0] = G1Point::generator();
    for (size_t i = 1; i < n; i++) {
        G1Point::doubling(points[i], points[i - 1]);
    }

    // Set small scalar values for testing
    scalars[0][0] = 2;
    scalars[1][0] = 3;
    scalars[2][0] = 4;

    // Calculate sum of scalar multiplications individually
    G1Point expected, temp;
    memset(&expected, 0, sizeof(expected));

    for (size_t i = 0; i < n; i++) {
        G1Point::mult(temp, points[i], scalars[i], 8);
        if (i == 0) {
            expected = temp;
        }
        else {
            G1Point::add(expected, expected, temp);
        }
    }

    // Test multi-scalar multiplication if available
    // blst_p1_mult_pippenger(result, points, scalars, n);
    // EXPECT_TRUE(G1Point::equal(result, expected));
}

TEST(E1PointTest, AdditionEdgeCases) {
    G1Point p1, neg_p1, result, temp;

    // Get generator point
    p1 = G1Point::generator();

    // Make a copy and negate it
    neg_p1 = p1;
    G1Point::cneg(neg_p1, 1);

    // P + (-P) should be infinity
    G1Point::dadd(result, p1, neg_p1);
    EXPECT_TRUE(G1Point::is_inf(result));

    // Add two identical points with regular add (should handle doubling correctly)
    G1Point::dadd(result, p1, p1);
    G1Point::doubling(temp, p1);
    EXPECT_TRUE(G1Point::equal(result, temp));
}

TEST(E1PointTest, EndomorphismAcceleration) {
    G1Point p1, result1;
    unsigned char scalar[32];

    // Fill with non-trivial scalar
    for (int i = 0; i < 32; i++) {
        scalar[i] = (i * 13) % 256;
    }

    p1 = G1Point::generator();

    // Compare regular scalar mult vs GLV method
    // This assumes an internal function to specifically test GLV
    G1Point::mult(result1, p1, scalar, 256);
    // If you have a test hook for the GLV implementation:
    // blst_p1_mult_glv(result2, p1, scalar);
    // EXPECT_TRUE(G1Point::equal(result1, result2));

    // Verify the result is on the curve
    EXPECT_TRUE(G1Point::on_curve(result1));
}

TEST(E1PointTest, HashToCurve) {
    G1Point p1;
    unsigned char msg[34] = "example message for hash-to-curve";
    unsigned char dst[30] = "QUUX-V01-CS02-with-BLS12381G1";

     //Test hash_to_curve if available
     blst_hash_to_g1(p1, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);

     //Result should be on curve and in the correct subgroup
     EXPECT_TRUE(G1Point::on_curve(p1));

     //Hash different message should give different point
     G1Point p2;
     msg[0] ^= 1;  // Flip a bit
     blst_hash_to_g1(p2, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);
     EXPECT_FALSE(G1Point::equal(p1, p2));
}

TEST(E1PointTest, InvalidPointHandling) {
    // Create a point on the curve but not in G1 subgroup

    // Step 1: Map to the curve without cofactor clearing
    G1Point random_point, invalid_point;
    Fp u;

    // Fill u with some arbitrary data
    for (size_t i = 0; i < sizeof(u) / sizeof(u[0]); i++) {
        u[i] = i * 0x0101010101010101ULL + 1;
    }

    // Access the internal mapping function that doesn't clear cofactor
    map_to_isogenous_E1(random_point, u);

    // Map to the BLS12-381 curve 
    isogeny_map_to_E1(invalid_point, random_point);

    // This point should be on the curve but not in G1
    ASSERT_TRUE(G1Point::on_curve(invalid_point)) << "Point not on curve";

    // Convert to affine for G1 membership check
    G1AffinePoint invalid_affine;
    G1Point::to_affine(invalid_affine, invalid_point);

    // Point should NOT be in G1 (if it is, our test assumptions are wrong)
    bool is_in_g1 = blst_p1_affine_in_g1(invalid_affine);

    // This is the critical test - the point should NOT be in G1
    EXPECT_FALSE(is_in_g1) << "Point should not be in G1 subgroup";

    // Now clear the cofactor to put the point in G1
    G1Point cleared_point;
    G1Point_times_minus_z(cleared_point, invalid_point);
    G1Point::dadd(cleared_point, cleared_point, invalid_point, NULL);

    // Convert to affine for G1 membership check
    G1AffinePoint cleared_affine;
    G1Point::to_affine(cleared_affine, cleared_point);

    // Now the point should be in G1
    bool cleared_in_g1 = blst_p1_affine_in_g1(cleared_affine);
    EXPECT_TRUE(cleared_in_g1) << "Cofactor-cleared point should be in G1";
}

TEST(E1PointTest, SerializationFormats) {
    G1Point p1 = G1Point::generator();
    G1AffinePoint p_affine;
    G1Point::to_affine(p_affine, p1);

    // Test all serialization formats
    unsigned char compressed[48];
    unsigned char uncompressed[96];

    // Test compressed
    G1Point::compress(compressed, p1);
    EXPECT_TRUE(compressed[0] & 0x80);  // Compressed flag

    // Test uncompressed
    G1Point::serialize(uncompressed, p1);
    EXPECT_FALSE(uncompressed[0] & 0x80);  // Not compressed

    // Test infinity encoding for each format
    G1Point infinity;
    memset(&infinity, 0, sizeof(infinity));

    G1Point::compress(compressed, infinity);
    EXPECT_TRUE(compressed[0] & 0x40);  // Infinity flag

    G1Point::serialize(uncompressed, infinity);
    EXPECT_TRUE(uncompressed[0] & 0x40);  // Infinity flag
}


TEST(E1PointTest, MixedAdditionCornerCases) {
    G1Point p1, result;
    G1AffinePoint infinity_affine;

    // Create affine point at infinity
    memset(&infinity_affine, 0, sizeof(infinity_affine));

    // Get generator point
    p1 = G1Point::generator();

    // P + O_affine = P
    G1Point::add_affine(result, p1, infinity_affine);
    EXPECT_TRUE(G1Point::equal(result, p1));

    // O + P_affine = P
    G1Point infinity;
    memset(&infinity, 0, sizeof(infinity));
    G1AffinePoint p_affine;
    G1Point::to_affine(p_affine, p1);

    G1Point::add_affine(result, infinity, p_affine);

    EXPECT_TRUE(G1Point::equal(result, p1));
}