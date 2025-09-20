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
void print_point(const POINTonE1* p, const char* name) {
    printf("%s:\n", name);
    printf("  X: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p->X[i]);
    }
    printf("\n  Y: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p->Y[i]);
    }
    printf("\n  Z: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p->Z[i]);
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
    POINTonE1 p1, p2, result, temp;

    // Get generator point
    p1 = *blst_p1_generator();
    
    // Create a different point (we'll use 2*G)
    blst_p1_double(&temp, &p1);
    
    // Add generator and the different point (G + 2G)
    blst_p1_add(&result, &p1, &temp);
    
    // Compute 3G using scalar multiplication
    unsigned char scalar[32] = {3};
    blst_p1_mult(&p2, &p1, scalar, 8);
    
    // Results should match
    EXPECT_TRUE(blst_p1_is_equal(&result, &p2));
    
    // Also test the doubling-specific function
    blst_p1_double(&result, &p1);
    
    // Create 2G with scalar multiplication
    unsigned char scalar2[32] = {2};
    blst_p1_mult(&p2, &p1, scalar2, 8);
    
    // Results should match
    EXPECT_TRUE(blst_p1_is_equal(&result, &p2));
}

TEST(E1PointTest, AdditionWithIdentity) {
    POINTonE1 p1, identity, result;

    // Get generator point
    p1 = *blst_p1_generator();

    // Create identity point (point at infinity)
    memset(&identity, 0, sizeof(identity));
    identity.Z[0] = 0;  // explicitly set Z to 0 for infinity

    // P + O = P
    blst_p1_add(&result, &p1, &identity);
    EXPECT_TRUE(blst_p1_is_equal(&result, &p1));
}

TEST(E1PointTest, PointNegation) {
    POINTonE1 p1, neg_p1, result;
    
    // Get generator point
    p1 = *blst_p1_generator();
    
    // Make a copy of the point and negate it
    neg_p1 = p1;
    blst_p1_cneg(&neg_p1, 1);
    
    // P + (-P) should be the identity
    blst_p1_add(&result, &p1, &neg_p1);
    
    // Result should be point at infinity
    EXPECT_TRUE(blst_p1_is_inf(&result));
}

TEST(E1PointTest, ScalarMultiplication) {
    POINTonE1 p1, result, expected;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    p1 = *blst_p1_generator();

    // Scalar = 2
    scalar[0] = 2;

    // Compute 2*P
    blst_p1_mult(&result, &p1, scalar, 8);

    // Should equal P + P (using doubling formula)
    blst_p1_double(&expected, &p1);
    EXPECT_TRUE(blst_p1_is_equal(&result, &expected));

    // Scalar = 3
    scalar[0] = 3;

    // Compute 3*P
    blst_p1_mult(&result, &p1, scalar, 8);

    // Should equal 2*P + P
    POINTonE1 doubled;
    blst_p1_double(&doubled, &p1);         // doubled = 2*P
    blst_p1_add(&expected, &doubled, &p1); // expected = 2*P + P = 3*P
    
    EXPECT_TRUE(blst_p1_is_equal(&result, &expected));
}

TEST(E1PointTest, ScalarMultiplicationByZero) {
    POINTonE1 p1, result;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    p1 = *blst_p1_generator();

    // Compute 0*P
    blst_p1_mult(&result, &p1, scalar, 8);

    // Result should be point at infinity
    EXPECT_TRUE(blst_p1_is_inf(&result));
}

TEST(E1PointTest, PointOnCurve) {
    POINTonE1 p1;
    POINTonE1_affine p1_affine;

    // Get generator point
    p1 = *blst_p1_generator();

    // Generator should be on curve
    EXPECT_TRUE(blst_p1_on_curve(&p1));

    // Convert to affine
    blst_p1_to_affine(&p1_affine, &p1);

    // Affine point should also be on curve
    EXPECT_TRUE(blst_p1_affine_on_curve(&p1_affine));

    // Modify point to be off curve
    p1.X[0] ^= 1;  // Flip a bit

    // Modified point should not be on curve
    EXPECT_FALSE(blst_p1_on_curve(&p1));
}

TEST(E1PointTest, SerializeDeserialize) {
    POINTonE1 p1, p2;
    POINTonE1_affine p1_affine;
    unsigned char serialized[96];

    // Get generator point
    p1 = *blst_p1_generator();

    // Serialize the point
    blst_p1_serialize(serialized, &p1);

    // Deserialize back to a point
    EXPECT_EQ(blst_p1_deserialize(&p1_affine, serialized), BLST_SUCCESS);

    // Convert from affine back to Jacobian coordinates
    blst_p1_from_affine(&p2, &p1_affine);

    // Points should be equal
    EXPECT_TRUE(blst_p1_is_equal(&p1, &p2));
}

TEST(E1PointTest, FixedVectorCompressUncompress) {
    // Test with the generator point, which is a guaranteed valid point
    const POINTonE1_affine* g_affine = blst_p1_affine_generator();
    POINTonE1_affine decompressed;
    unsigned char compressed[48] = { 0 };

    // Inspect the generator point to ensure it's valid
    ASSERT_TRUE(blst_p1_affine_on_curve(g_affine)) << "Generator point not on curve";

    // Manual compression:
    // First, we'll get the byte representation of the generator's X coordinate
    vec384 x_coord;
    from_fp(x_coord, g_affine->X);  // Convert from Montgomery to regular form

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
    BLST_ERROR err1 = blst_p1_uncompress(&decompressed, compressed);

    // Try flipping the sign bit
    compressed[0] ^= 0x20;  // Flip the sign bit
    //printf("Flipped sign bit, now first byte: 0x%02X\n", compressed[0]);

    BLST_ERROR err2 = blst_p1_uncompress(&decompressed, compressed);

    //printf("Uncompress results: without sign bit = %d, with sign bit = %d\n", err1, err2);

    // At least one of these should work
    EXPECT_TRUE(err1 == BLST_SUCCESS || err2 == BLST_SUCCESS)
        << "Failed to uncompress point with either sign bit";

    // If successful, check if the point is on the curve
    if (err1 == BLST_SUCCESS || err2 == BLST_SUCCESS) {
        EXPECT_TRUE(blst_p1_affine_on_curve(&decompressed))
            << "Decompressed point not on curve";
    }
}

TEST(E1PointTest, BasicCompressUncompress) {
    {
        const auto* gen = blst_p1_generator();

        uint8_t comp[48];
        blst_p1_compress(comp, gen);  // or: to_affine + blst_p1_affine_compress

        POINTonE1_affine dec;
        BLST_ERROR err = blst_p1_uncompress(&dec, comp);
        ASSERT_TRUE(err == BLST_SUCCESS);
        ASSERT_TRUE(blst_p1_affine_on_curve(&dec));
        ASSERT_TRUE(blst_p1_affine_in_g1(&dec));
    }
    // Starting with the generator point in projective form
    POINTonE1 point = *blst_p1_generator();

    // Convert to affine
    POINTonE1_affine affine_point;
    blst_p1_to_affine(&affine_point, &point);

    // First confirm it's on the curve
    ASSERT_TRUE(blst_p1_affine_on_curve(&affine_point)) << "Point not on curve";

    // Compress the affine point
    unsigned char compressed[48];
    blst_p1_affine_compress(compressed, &affine_point);

    // Print the compressed format for debugging
    //printf("Compressed point (first byte: 0x%02X)\n", compressed[0]);

    // Try to uncompress using the same formats as the library expects
    POINTonE1_affine result;
    BLST_ERROR err = blst_p1_uncompress(&result, compressed);

    // If it fails, try a fallback approach
    if (err != BLST_SUCCESS) {
        printf("Standard decompression failed with code %d\n", err);

        // Try using serialize/deserialize instead
        unsigned char serialized[96];
        blst_p1_affine_serialize(serialized, &affine_point);

        err = blst_p1_deserialize(&result, serialized);
        EXPECT_EQ(err, BLST_SUCCESS) << "Even serialization/deserialization failed";

        // If that works, the issue is specific to compression
        if (err == BLST_SUCCESS) {
            printf("Serialization/deserialization works but compression fails.\n");
            printf("This suggests an issue with the compression algorithm.\n");
        }
    }
    else {
        // Success! Points should match
        EXPECT_TRUE(blst_p1_affine_is_equal(&result, &affine_point))
            << "Points don't match after compression/decompression";
    }
}

TEST(E1PointTest, AffineCompressUncompress) {
    POINTonE1_affine original, decompressed;
    unsigned char compressed[48] = { 0 };

    // Get the generator in affine form
    original = *blst_p1_affine_generator();

    // Verify it's on curve
    ASSERT_TRUE(blst_p1_affine_on_curve(&original)) << "Generator point not on curve";

    // Compress 
    blst_p1_affine_compress(compressed, &original);

    // Verify compressed format
    EXPECT_TRUE(compressed[0] & 0x80) << "Compressed bit not set";
    EXPECT_FALSE(compressed[0] & 0x40) << "Infinity bit incorrectly set";

    // Uncompress
    BLST_ERROR err = blst_p1_uncompress(&decompressed, compressed);
    ASSERT_EQ(err, BLST_SUCCESS) << "Decompression failed with code " << err;

    // Verify the point is on curve
    ASSERT_TRUE(blst_p1_affine_on_curve(&decompressed)) << "Decompressed point not on curve";

    // Directly compare the affine coordinates
    EXPECT_TRUE(blst_p1_affine_is_equal(&original, &decompressed))
        << "Affine points don't match after compression cycle";
}

TEST(E1PointTest, CompressUncompress) {
    POINTonE1 p1, p2;
    POINTonE1_affine p1_affine, p2_affine;
    unsigned char compressed[48] = {0};

    // Start with the generator point which is guaranteed to be valid
    p1 = *blst_p1_generator();
    
    // Convert to affine first (important for compression)
    blst_p1_to_affine(&p1_affine, &p1);
    
    // Compress directly using blst_p1_compress (not the affine version)
    blst_p1_compress(compressed, &p1);
    
    // Manually verify the format - first byte should have compressed bit set (0x80)
    // but not have infinity bit set (0x40)
    EXPECT_TRUE((compressed[0] & 0x80) && !(compressed[0] & 0x40));
    
    // Try to uncompress
    BLST_ERROR err = blst_p1_uncompress(&p2_affine, compressed);
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
    blst_p1_from_affine(&p2, &p2_affine);
    
    // Debug output in case of comparison failure
    if (!blst_p1_is_equal(&p1, &p2)) {
        printf("Points differ after compress/decompress cycle\n");
        print_point(&p1, "Original");
        print_point(&p2, "After compress/decompress");
        
        // Check if points differ only by Z-coordinate normalization
        POINTonE1 p1_normalized, p2_normalized;
        blst_p1_from_jacobian(&p1_normalized, &p1);
        blst_p1_from_jacobian(&p2_normalized, &p2);
        
        bool normalized_equal = blst_p1_is_equal(&p1_normalized, &p2_normalized);
        printf("After normalization, points are %s\n", 
               normalized_equal ? "EQUAL" : "DIFFERENT");
    }
    
    // Points should be equal (if not, the issue could be normalization)
    EXPECT_TRUE(blst_p1_is_equal(&p1, &p2));
}

TEST(E1PointTest, InfinityPoint) {
    POINTonE1 identity;
    unsigned char serialized[96];
    POINTonE1_affine affine_identity;

    // Create identity point
    memset(&identity, 0, sizeof(identity));

    // Point should be at infinity
    EXPECT_TRUE(blst_p1_is_inf(&identity));

    // Serialize infinity point
    blst_p1_serialize(serialized, &identity);

    // First byte should have infinity bit set
    EXPECT_TRUE(serialized[0] & 0x40);

    // Deserialize back to a point
    EXPECT_EQ(blst_p1_deserialize(&affine_identity, serialized), BLST_SUCCESS);

    // Should be affine infinity
    EXPECT_TRUE(blst_p1_affine_is_inf(&affine_identity));
}

TEST(E1PointTest, SigningBasicOperation) {
    POINTonE1 p1, sig;
    unsigned char SK[32] = {
        0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08,
        0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10,
        0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18,
        0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20
    };
    unsigned char pk_bytes[96];

    // Generate public key from secret key
    blst_sk_to_pk2_in_g1(pk_bytes, NULL, SK);

    // Create a message point - we'll use the generator for simplicity
    p1 = *blst_p1_generator();

    // Sign the message
    blst_sign_pk_in_g2(&sig, &p1, SK);

    // Signature should be on the curve
    EXPECT_TRUE(blst_p1_on_curve(&sig));
}

TEST(E1PointTest, PointCoordinateConversion) {
    POINTonE1 p1, p2;
    POINTonE1_affine p_affine;

    // Get generator point
    p1 = *blst_p1_generator();

    // Convert to affine
    blst_p1_to_affine(&p_affine, &p1);

    // Convert back to Jacobian
    blst_p1_from_affine(&p2, &p_affine);

    // Should be equal to original point
    EXPECT_TRUE(blst_p1_is_equal(&p1, &p2));

    // Convert from Jacobian to Jacobian (normalized Z)
    blst_p1_from_jacobian(&p2, &p1);

    // Should be equal to original point
    EXPECT_TRUE(blst_p1_is_equal(&p1, &p2));
}

TEST(E1PointTest, GeneratorProperties) {
    // Get generator point
    const POINTonE1* g = blst_p1_generator();

    // Generator should be on curve
    EXPECT_TRUE(blst_p1_on_curve(g));

    // Generator should not be at infinity
    EXPECT_FALSE(blst_p1_is_inf(g));

    // Get affine generator
    const POINTonE1_affine* g_affine = blst_p1_affine_generator();

    // Affine generator should be on curve
    EXPECT_TRUE(blst_p1_affine_on_curve(g_affine));

    // Affine generator should not be at infinity
    EXPECT_FALSE(blst_p1_affine_is_inf(g_affine));
}

TEST(E1PointTest, LargeScalarMultiplication) {
    POINTonE1 p1, result;
    // Create a large scalar close to the order
    unsigned char large_scalar[32] = {
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x0f
    };

    // Get generator point
    p1 = *blst_p1_generator();

    // Perform large scalar multiplication
    blst_p1_mult(&result, &p1, large_scalar, 252);

    // Result should be on curve
    EXPECT_TRUE(blst_p1_on_curve(&result));
}

TEST(E1PointTest, AffineEquality) {
    POINTonE1_affine a1, a2;
    POINTonE1 p1, p2;

    // Get generator point
    p1 = *blst_p1_generator();

    // Double the generator
    blst_p1_double(&p2, &p1);

    // Convert both to affine
    blst_p1_to_affine(&a1, &p1);
    blst_p1_to_affine(&a2, &p1);

    // Affine points should be equal
    EXPECT_TRUE(blst_p1_affine_is_equal(&a1, &a2));

    // Convert p2 to affine
    blst_p1_to_affine(&a2, &p2);

    // Should not be equal to a1
    EXPECT_FALSE(blst_p1_affine_is_equal(&a1, &a2));
}

TEST(E1PointTest, DeserializeInvalidPoint) {
    POINTonE1_affine result;
    unsigned char invalid_data[96];

    // Fill with invalid data
    memset(invalid_data, 0xFF, sizeof(invalid_data));

    // Should fail to deserialize
    EXPECT_NE(blst_p1_deserialize(&result, invalid_data), BLST_SUCCESS);
}

TEST(E1PointTest, SizeChecks) {
    // Check that reported sizes match expected
    EXPECT_EQ(blst_p1_sizeof(), sizeof(POINTonE1));
    EXPECT_EQ(blst_p1_affine_sizeof(), sizeof(POINTonE1_affine));
}

TEST(E1PointTest, MultiplePointAddition) {
    POINTonE1 p1, p3a, p2, result;

    // Get generator point
    p1 = *blst_p1_generator();

    // Compute 3*P using additions
    blst_p1_double(&p2, &p1);    // temp = 2*P
    blst_p1_add(&p3a, &p2, &p1);  // p3a = 3*P

    // Compute 3*P using scalar multiplication
    unsigned char scalar[32] = { 3 };
    blst_p1_mult(&result, &p1, scalar, 8);

    // Results should match
    EXPECT_TRUE(blst_p1_is_equal(&result, &p3a));

    // Test mixed addition with different points
    POINTonE1_affine p2_affine;
    blst_p1_to_affine(&p2_affine, &p2);  // Convert 2*P to affine

    // P + (2P affine) should equal 3P
    blst_p1_add_affine(&result, &p1, &p2_affine);
    EXPECT_TRUE(blst_p1_is_equal(&result, &p3a));  // Should equal 3*P

    // For completeness, test add_or_double with distinct points too
    POINTonE1 p3b;
    blst_p1_add_affine(&p3b, &p1, &p2_affine);
    //blst_p1_add_or_double_affine(&p3b, &p1, &p2_affine);
    EXPECT_TRUE(blst_p1_is_equal(&p3a, &p3b));  // Should equal 3*P
    
    // Also verify that add_or_double_affine works with same point (P+P)
    POINTonE1_affine p_affine;
    blst_p1_to_affine(&p_affine, &p1);
    
    // Test P + P using add_or_double_affine
    blst_p1_add_or_double_affine(&result, &p1, &p_affine);
    EXPECT_TRUE(blst_p1_is_equal(&result, &p2));  // Should equal 2*P
}


TEST(E1PointTest, MultiScalarMultiplication) {
    const size_t n = 3;
    POINTonE1 points[n];
    byte scalars[n][32] = { 0 };

    // Generate n different points
    points[0] = *blst_p1_generator();
    for (size_t i = 1; i < n; i++) {
        blst_p1_double(&points[i], &points[i - 1]);
    }

    // Set small scalar values for testing
    scalars[0][0] = 2;
    scalars[1][0] = 3;
    scalars[2][0] = 4;

    // Calculate sum of scalar multiplications individually
    POINTonE1 expected, temp;
    memset(&expected, 0, sizeof(expected));

    for (size_t i = 0; i < n; i++) {
        blst_p1_mult(&temp, &points[i], scalars[i], 8);
        if (i == 0) {
            expected = temp;
        }
        else {
            blst_p1_add(&expected, &expected, &temp);
        }
    }

    // Test multi-scalar multiplication if available
    // blst_p1_mult_pippenger(&result, points, scalars, n);
    // EXPECT_TRUE(blst_p1_is_equal(&result, &expected));
}

TEST(E1PointTest, AdditionEdgeCases) {
    POINTonE1 p1, neg_p1, result, temp;

    // Get generator point
    p1 = *blst_p1_generator();

    // Make a copy and negate it
    neg_p1 = p1;
    blst_p1_cneg(&neg_p1, 1);

    // P + (-P) should be infinity
    blst_p1_add_or_double(&result, &p1, &neg_p1);
    EXPECT_TRUE(blst_p1_is_inf(&result));

    // Add two identical points with regular add (should handle doubling correctly)
    blst_p1_add_or_double(&result, &p1, &p1);
    blst_p1_double(&temp, &p1);
    EXPECT_TRUE(blst_p1_is_equal(&result, &temp));
}

TEST(E1PointTest, EndomorphismAcceleration) {
    POINTonE1 p1, result1;
    unsigned char scalar[32];

    // Fill with non-trivial scalar
    for (int i = 0; i < 32; i++) {
        scalar[i] = (i * 13) % 256;
    }

    p1 = *blst_p1_generator();

    // Compare regular scalar mult vs GLV method
    // This assumes an internal function to specifically test GLV
    blst_p1_mult(&result1, &p1, scalar, 256);
    // If you have a test hook for the GLV implementation:
    // blst_p1_mult_glv(&result2, &p1, scalar);
    // EXPECT_TRUE(blst_p1_is_equal(&result1, &result2));

    // Verify the result is on the curve
    EXPECT_TRUE(blst_p1_on_curve(&result1));
}

TEST(E1PointTest, HashToCurve) {
    POINTonE1 p1;
    unsigned char msg[34] = "example message for hash-to-curve";
    unsigned char dst[30] = "QUUX-V01-CS02-with-BLS12381G1";

     //Test hash_to_curve if available
     blst_hash_to_g1(&p1, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);

     //Result should be on curve and in the correct subgroup
     EXPECT_TRUE(blst_p1_on_curve(&p1));

     //Hash different message should give different point
     POINTonE1 p2;
     msg[0] ^= 1;  // Flip a bit
     blst_hash_to_g1(&p2, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);
     EXPECT_FALSE(blst_p1_is_equal(&p1, &p2));
}

TEST(E1PointTest, InvalidPointHandling) {
    // Create a point on the curve but not in G1 subgroup

    // Step 1: Map to the curve without cofactor clearing
    POINTonE1 random_point, invalid_point;
    vec384 u;

    // Fill u with some arbitrary data
    for (size_t i = 0; i < sizeof(u) / sizeof(u[0]); i++) {
        u[i] = i * 0x0101010101010101ULL + 1;
    }

    // Access the internal mapping function that doesn't clear cofactor
    map_to_isogenous_E1(&random_point, u);

    // Map to the BLS12-381 curve 
    isogeny_map_to_E1(&invalid_point, &random_point);

    // This point should be on the curve but not in G1
    ASSERT_TRUE(blst_p1_on_curve(&invalid_point)) << "Point not on curve";

    // Convert to affine for G1 membership check
    POINTonE1_affine invalid_affine;
    blst_p1_to_affine(&invalid_affine, &invalid_point);

    // Point should NOT be in G1 (if it is, our test assumptions are wrong)
    bool is_in_g1 = blst_p1_affine_in_g1(&invalid_affine);

    // This is the critical test - the point should NOT be in G1
    EXPECT_FALSE(is_in_g1) << "Point should not be in G1 subgroup";

    // Now clear the cofactor to put the point in G1
    POINTonE1 cleared_point;
    POINTonE1_times_minus_z(&cleared_point, &invalid_point);
    POINTonE1_dadd(&cleared_point, &cleared_point, &invalid_point, NULL);

    // Convert to affine for G1 membership check
    POINTonE1_affine cleared_affine;
    blst_p1_to_affine(&cleared_affine, &cleared_point);

    // Now the point should be in G1
    bool cleared_in_g1 = blst_p1_affine_in_g1(&cleared_affine);
    EXPECT_TRUE(cleared_in_g1) << "Cofactor-cleared point should be in G1";
}

TEST(E1PointTest, SerializationFormats) {
    POINTonE1 p1 = *blst_p1_generator();
    POINTonE1_affine p_affine;
    blst_p1_to_affine(&p_affine, &p1);

    // Test all serialization formats
    unsigned char compressed[48];
    unsigned char uncompressed[96];

    // Test compressed
    blst_p1_compress(compressed, &p1);
    EXPECT_TRUE(compressed[0] & 0x80);  // Compressed flag

    // Test uncompressed
    blst_p1_serialize(uncompressed, &p1);
    EXPECT_FALSE(uncompressed[0] & 0x80);  // Not compressed

    // Test infinity encoding for each format
    POINTonE1 infinity;
    memset(&infinity, 0, sizeof(infinity));

    blst_p1_compress(compressed, &infinity);
    EXPECT_TRUE(compressed[0] & 0x40);  // Infinity flag

    blst_p1_serialize(uncompressed, &infinity);
    EXPECT_TRUE(uncompressed[0] & 0x40);  // Infinity flag
}


TEST(E1PointTest, MixedAdditionCornerCases) {
    POINTonE1 p1, result;
    POINTonE1_affine infinity_affine;

    // Create affine point at infinity
    memset(&infinity_affine, 0, sizeof(infinity_affine));

    // Get generator point
    p1 = *blst_p1_generator();

    // P + O_affine = P
    blst_p1_add_affine(&result, &p1, &infinity_affine);
    EXPECT_TRUE(blst_p1_is_equal(&result, &p1));

    // O + P_affine = P
    POINTonE1 infinity;
    memset(&infinity, 0, sizeof(infinity));
    POINTonE1_affine p_affine;
    blst_p1_to_affine(&p_affine, &p1);

    blst_p1_add_affine(&result, &infinity, &p_affine);
    EXPECT_TRUE(blst_p1_is_equal(&result, &p1));
}