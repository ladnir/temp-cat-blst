#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

extern "C" {
#include "../src/point.h"
#include "../src/fields.h"
#include "../src/e1.h"
#include "../src/e2.h"
#include "../src/map_to_g2.h"
}

// Helper: Print G2 point coordinates for debugging
void print_point_e2(const POINTonE2* p, const char* name) {
    printf("%s:\n", name);
    printf("  X.Re: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p->X[0][i]);
    }
    printf("\n  X.Im: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p->X[1][i]);
    }
    printf("\n  Y.Re: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p->Y[0][i]);
    }
    printf("\n  Y.Im: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p->Y[1][i]);
    }
    printf("\n  Z.Re: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p->Z[0][i]);
    }
    printf("\n  Z.Im: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)p->Z[1][i]);
    }
    printf("\n");
}

// Helper: Convert hex string to bytes
void hex_to_bytes_e2(unsigned char* bytes, const char* hex) {
    size_t len = strlen(hex) / 2;
    for (size_t i = 0; i < len; i++) {
        sscanf(hex + (i * 2), "%2hhx", &bytes[i]);
    }
}

TEST(E2PointTest, PointAddition) {
    POINTonE2 p1, p2, result, temp;

    // Get generator point
    p1 = *blst_p2_generator();

    // Create a different point (we'll use 2*G)
    blst_p2_double(&temp, &p1);

    // Add generator and the different point (G + 2G)
    blst_p2_add(&result, &p1, &temp);

    // Compute 3G using scalar multiplication
    unsigned char scalar[32] = { 3 };
    blst_p2_mult(&p2, &p1, scalar, 8);

    // Results should match
    EXPECT_TRUE(blst_p2_is_equal(&result, &p2));

    // Also test the doubling-specific function
    blst_p2_double(&result, &p1);

    // Create 2G with scalar multiplication
    unsigned char scalar2[32] = { 2 };
    blst_p2_mult(&p2, &p1, scalar2, 8);

    // Results should match
    EXPECT_TRUE(blst_p2_is_equal(&result, &p2));
}

TEST(E2PointTest, AdditionWithIdentity) {
    POINTonE2 p1, identity, result;

    // Get generator point
    p1 = *blst_p2_generator();

    // Create identity point (point at infinity)
    memset(&identity, 0, sizeof(identity));
    identity.Z[0][0] = 0;  // explicitly set Z to 0 for infinity

    // P + O = P
    blst_p2_add(&result, &p1, &identity);
    EXPECT_TRUE(blst_p2_is_equal(&result, &p1));
}

TEST(E2PointTest, PointNegation) {
    POINTonE2 p1, neg_p1, result;

    // Get generator point
    p1 = *blst_p2_generator();

    // Make a copy of the point and negate it
    neg_p1 = p1;
    blst_p2_cneg(&neg_p1, 1);

    // P + (-P) should be the identity
    blst_p2_add(&result, &p1, &neg_p1);

    // Result should be point at infinity
    EXPECT_TRUE(blst_p2_is_inf(&result));
}

TEST(E2PointTest, ScalarMultiplication) {
    POINTonE2 p1, result, expected;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    p1 = *blst_p2_generator();

    // Scalar = 2
    scalar[0] = 2;

    // Compute 2*P
    blst_p2_mult(&result, &p1, scalar, 8);

    // Should equal P + P (using doubling formula)
    blst_p2_double(&expected, &p1);
    EXPECT_TRUE(blst_p2_is_equal(&result, &expected));

    // Scalar = 3
    scalar[0] = 3;

    // Compute 3*P
    blst_p2_mult(&result, &p1, scalar, 8);

    // Should equal 2*P + P
    POINTonE2 doubled;
    blst_p2_double(&doubled, &p1);         // doubled = 2*P
    blst_p2_add(&expected, &doubled, &p1); // expected = 2*P + P = 3*P

    EXPECT_TRUE(blst_p2_is_equal(&result, &expected));
}

TEST(E2PointTest, ScalarMultiplicationByZero) {
    POINTonE2 p1, result;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    p1 = *blst_p2_generator();

    // Compute 0*P
    blst_p2_mult(&result, &p1, scalar, 8);

    // Result should be point at infinity
    EXPECT_TRUE(blst_p2_is_inf(&result));
}

TEST(E2PointTest, PointOnCurve) {
    POINTonE2 p1;
    POINTonE2_affine p1_affine;

    // Get generator point
    p1 = *blst_p2_generator();

    // Generator should be on curve
    EXPECT_TRUE(blst_p2_on_curve(&p1));

    // Convert to affine
    blst_p2_to_affine(&p1_affine, &p1);

    // Affine point should also be on curve
    EXPECT_TRUE(blst_p2_affine_on_curve(&p1_affine));

    // Modify point to be off curve
    p1.X[0][0] ^= 1;  // Flip a bit

    // Modified point should not be on curve
    EXPECT_FALSE(blst_p2_on_curve(&p1));
}

TEST(E2PointTest, SerializeDeserialize) {
    POINTonE2 p1, p2;
    POINTonE2_affine p1_affine;
    unsigned char serialized[192];  // 2 * 96 bytes for G2 point

    // Get generator point
    p1 = *blst_p2_generator();

    // Serialize the point
    blst_p2_serialize(serialized, &p1);

    // Deserialize back to a point
    EXPECT_EQ(blst_p2_deserialize(&p1_affine, serialized), BLST_SUCCESS);

    // Convert from affine back to Jacobian coordinates
    blst_p2_from_affine(&p2, &p1_affine);

    // Points should be equal
    EXPECT_TRUE(blst_p2_is_equal(&p1, &p2));
}

extern "C" {
    int blst_p2_affine_in_g2(const POINTonE2_affine* p);
}

TEST(E2PointTest, BasicCompressUncompress) {
    // Starting with the generator point in projective form
    POINTonE2 point = *blst_p2_generator();

    // Convert to affine
    POINTonE2_affine affine_point;
    blst_p2_to_affine(&affine_point, &point);

    // First confirm it's on the curve
    ASSERT_TRUE(blst_p2_affine_on_curve(&affine_point)) << "Point not on curve";

    // Compress the affine point
    unsigned char compressed[96];  // 2 * 48 bytes for compressed G2
    blst_p2_affine_compress(compressed, &affine_point);

    // Uncompress
    POINTonE2_affine result;
    BLST_ERROR err = blst_p2_uncompress(&result, compressed);
    ASSERT_EQ(err, BLST_SUCCESS) << "Decompression failed with code " << err;

    // Points should match
    EXPECT_TRUE(blst_p2_affine_is_equal(&result, &affine_point))
        << "Points don't match after compression/decompression";
}

TEST(E2PointTest, AffineCompressUncompress) {
    POINTonE2_affine original, decompressed;
    unsigned char compressed[96] = { 0 };

    // Get the generator in affine form
    original = *blst_p2_affine_generator();

    // Verify it's on curve
    ASSERT_TRUE(blst_p2_affine_on_curve(&original)) << "Generator point not on curve";

    // Compress 
    blst_p2_affine_compress(compressed, &original);

    // Verify compressed format
    EXPECT_TRUE(compressed[0] & 0x80) << "Compressed bit not set";
    EXPECT_FALSE(compressed[0] & 0x40) << "Infinity bit incorrectly set";

    // Uncompress
    BLST_ERROR err = blst_p2_uncompress(&decompressed, compressed);
    ASSERT_EQ(err, BLST_SUCCESS) << "Decompression failed with code " << err;

    // Verify the point is on curve
    ASSERT_TRUE(blst_p2_affine_on_curve(&decompressed)) << "Decompressed point not on curve";

    // Directly compare the affine coordinates
    EXPECT_TRUE(blst_p2_affine_is_equal(&original, &decompressed))
        << "Affine points don't match after compression cycle";
}

TEST(E2PointTest, CompressUncompress) {
    POINTonE2 p1, p2;
    POINTonE2_affine p1_affine, p2_affine;
    unsigned char compressed[96] = { 0 };

    // Start with the generator point which is guaranteed to be valid
    p1 = *blst_p2_generator();

    // Convert to affine first (important for compression)
    blst_p2_to_affine(&p1_affine, &p1);

    // Compress directly
    blst_p2_compress(compressed, &p1);

    // Manually verify the format - first byte should have compressed bit set (0x80)
    // but not have infinity bit set (0x40)
    EXPECT_TRUE((compressed[0] & 0x80) && !(compressed[0] & 0x40));

    // Try to uncompress
    BLST_ERROR err = blst_p2_uncompress(&p2_affine, compressed);
    ASSERT_EQ(err, BLST_SUCCESS) << "Point decompression failed";

    // Convert back to Jacobian coordinates
    blst_p2_from_affine(&p2, &p2_affine);

    // Points should be equal
    EXPECT_TRUE(blst_p2_is_equal(&p1, &p2));
}

TEST(E2PointTest, InfinityPoint) {
    POINTonE2 identity;
    unsigned char serialized[192];
    POINTonE2_affine affine_identity;

    // Create identity point
    memset(&identity, 0, sizeof(identity));

    // Point should be at infinity
    EXPECT_TRUE(blst_p2_is_inf(&identity));

    // Serialize infinity point
    blst_p2_serialize(serialized, &identity);

    // First byte should have infinity bit set
    EXPECT_TRUE(serialized[0] & 0x40);

    // Deserialize back to a point
    EXPECT_EQ(blst_p2_deserialize(&affine_identity, serialized), BLST_SUCCESS);

    // Should be affine infinity
    EXPECT_TRUE(blst_p2_affine_is_inf(&affine_identity));
}

TEST(E2PointTest, SigningBasicOperation) {
    POINTonE2 p1, sig;
    unsigned char SK[32] = {
        0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08,
        0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10,
        0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18,
        0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20
    };
    unsigned char pk_bytes[192];

    // Generate public key from secret key
    blst_sk_to_pk2_in_g2(pk_bytes, NULL, SK);

    // Create a message point - we'll use the generator for simplicity
    p1 = *blst_p2_generator();

    // Sign the message
    blst_sign_pk_in_g1(&sig, &p1, SK);

    // Signature should be on the curve
    EXPECT_TRUE(blst_p2_on_curve(&sig));
}

TEST(E2PointTest, PointCoordinateConversion) {
    POINTonE2 p1, p2;
    POINTonE2_affine p_affine;

    // Get generator point
    p1 = *blst_p2_generator();

    // Convert to affine
    blst_p2_to_affine(&p_affine, &p1);

    // Convert back to Jacobian
    blst_p2_from_affine(&p2, &p_affine);

    // Should be equal to original point
    EXPECT_TRUE(blst_p2_is_equal(&p1, &p2));

    // Convert from Jacobian to Jacobian (normalized Z)
    blst_p2_from_jacobian(&p2, &p1);

    // Should be equal to original point
    EXPECT_TRUE(blst_p2_is_equal(&p1, &p2));
}

TEST(E2PointTest, GeneratorProperties) {
    // Get generator point
    const POINTonE2* g = blst_p2_generator();

    // Generator should be on curve
    EXPECT_TRUE(blst_p2_on_curve(g));

    // Generator should not be at infinity
    EXPECT_FALSE(blst_p2_is_inf(g));

    // Get affine generator
    const POINTonE2_affine* g_affine = blst_p2_affine_generator();

    // Affine generator should be on curve
    EXPECT_TRUE(blst_p2_affine_on_curve(g_affine));

    // Affine generator should not be at infinity
    EXPECT_FALSE(blst_p2_affine_is_inf(g_affine));
}

TEST(E2PointTest, LargeScalarMultiplication) {
    POINTonE2 p1, result;
    // Create a large scalar close to the order
    unsigned char large_scalar[32] = {
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x0f
    };

    // Get generator point
    p1 = *blst_p2_generator();

    // Perform large scalar multiplication
    blst_p2_mult(&result, &p1, large_scalar, 252);

    // Result should be on curve
    EXPECT_TRUE(blst_p2_on_curve(&result));
}

TEST(E2PointTest, AffineEquality) {
    POINTonE2_affine a1, a2;
    POINTonE2 p1, p2;

    // Get generator point
    p1 = *blst_p2_generator();

    // Double the generator
    blst_p2_double(&p2, &p1);

    // Convert both to affine
    blst_p2_to_affine(&a1, &p1);
    blst_p2_to_affine(&a2, &p1);

    // Affine points should be equal
    EXPECT_TRUE(blst_p2_affine_is_equal(&a1, &a2));

    // Convert p2 to affine
    blst_p2_to_affine(&a2, &p2);

    // Should not be equal to a1
    EXPECT_FALSE(blst_p2_affine_is_equal(&a1, &a2));
}

TEST(E2PointTest, DeserializeInvalidPoint) {
    POINTonE2_affine result;
    unsigned char invalid_data[192];

    // Fill with invalid data
    memset(invalid_data, 0xFF, sizeof(invalid_data));

    // Should fail to deserialize
    EXPECT_NE(blst_p2_deserialize(&result, invalid_data), BLST_SUCCESS);
}

TEST(E2PointTest, SizeChecks) {
    // Check that reported sizes match expected
    EXPECT_EQ(blst_p2_sizeof(), sizeof(POINTonE2));
    EXPECT_EQ(blst_p2_affine_sizeof(), sizeof(POINTonE2_affine));
}

TEST(E2PointTest, MultiplePointAddition) {
    POINTonE2 p1, p3a, p2, result;

    // Get generator point
    p1 = *blst_p2_generator();

    // Compute 3*P using additions
    blst_p2_double(&p2, &p1);    // temp = 2*P
    blst_p2_add(&p3a, &p2, &p1);  // p3a = 3*P

    // Compute 3*P using scalar multiplication
    unsigned char scalar[32] = { 3 };
    blst_p2_mult(&result, &p1, scalar, 8);

    // Results should match
    EXPECT_TRUE(blst_p2_is_equal(&result, &p3a));

    // Test mixed addition with different points
    POINTonE2_affine p2_affine;
    blst_p2_to_affine(&p2_affine, &p2);  // Convert 2*P to affine

    // P + (2P affine) should equal 3P
    blst_p2_add_affine(&result, &p1, &p2_affine);
    EXPECT_TRUE(blst_p2_is_equal(&result, &p3a));  // Should equal 3*P

    // Also verify that add_or_double_affine works with same point (P+P)
    POINTonE2_affine p_affine;
    blst_p2_to_affine(&p_affine, &p1);

    // Test P + P using add_or_double_affine
    blst_p2_add_or_double_affine(&result, &p1, &p_affine);
    EXPECT_TRUE(blst_p2_is_equal(&result, &p2));  // Should equal 2*P
}

TEST(E2PointTest, AdditionEdgeCases) {
    POINTonE2 p1, neg_p1, result, temp;

    // Get generator point
    p1 = *blst_p2_generator();

    // Make a copy and negate it
    neg_p1 = p1;
    blst_p2_cneg(&neg_p1, 1);

    // P + (-P) should be infinity
    blst_p2_add_or_double(&result, &p1, &neg_p1);
    EXPECT_TRUE(blst_p2_is_inf(&result));

    // Add two identical points with regular add (should handle doubling correctly)
    blst_p2_add_or_double(&result, &p1, &p1);
    blst_p2_double(&temp, &p1);
    EXPECT_TRUE(blst_p2_is_equal(&result, &temp));
}

TEST(E2PointTest, HashToCurve) {
    POINTonE2 p1;
    unsigned char msg[34] = "example message for hash-to-curve";
    unsigned char dst[30] = "QUUX-V01-CS02-with-BLS12381G2";

    // Test hash_to_curve
    blst_hash_to_g2(&p1, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);

    // Result should be on curve
    EXPECT_TRUE(blst_p2_on_curve(&p1));

    // Hash different message should give different point
    POINTonE2 p2;
    msg[0] ^= 1;  // Flip a bit
    blst_hash_to_g2(&p2, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);
    EXPECT_FALSE(blst_p2_is_equal(&p1, &p2));
}

TEST(E2PointTest, InvalidPointHandling) {
    // This test requires internal functions to create a point on E2 but not in G2
    // We can simulate by hashing to curve without cofactor clearing

    // For now, just verify that the generator is in G2
    POINTonE2 point = *blst_p2_generator();
    POINTonE2_affine point_affine;
    blst_p2_to_affine(&point_affine, &point);

    // Generator should be in G2
    EXPECT_TRUE(blst_p2_affine_in_g2(&point_affine));
}

TEST(E2PointTest, SerializationFormats) {
    POINTonE2 p1 = *blst_p2_generator();
    POINTonE2_affine p_affine;
    blst_p2_to_affine(&p_affine, &p1);

    // Test all serialization formats
    unsigned char compressed[96];
    unsigned char uncompressed[192];

    // Test compressed
    blst_p2_compress(compressed, &p1);
    EXPECT_TRUE(compressed[0] & 0x80);  // Compressed flag

    // Test uncompressed
    blst_p2_serialize(uncompressed, &p1);
    EXPECT_FALSE(uncompressed[0] & 0x80);  // Not compressed

    // Test infinity encoding for each format
    POINTonE2 infinity;
    memset(&infinity, 0, sizeof(infinity));

    blst_p2_compress(compressed, &infinity);
    EXPECT_TRUE(compressed[0] & 0x40);  // Infinity flag

    blst_p2_serialize(uncompressed, &infinity);
    EXPECT_TRUE(uncompressed[0] & 0x40);  // Infinity flag
}

TEST(E2PointTest, MixedAdditionCornerCases) {
    POINTonE2 p1, result;
    POINTonE2_affine infinity_affine;

    // Create affine point at infinity
    memset(&infinity_affine, 0, sizeof(infinity_affine));

    // Get generator point
    p1 = *blst_p2_generator();

    // P + O_affine = P
    blst_p2_add_affine(&result, &p1, &infinity_affine);
    EXPECT_TRUE(blst_p2_is_equal(&result, &p1));

    // O + P_affine = P
    POINTonE2 infinity;
    memset(&infinity, 0, sizeof(infinity));
    POINTonE2_affine p_affine;
    blst_p2_to_affine(&p_affine, &p1);

    blst_p2_add_affine(&result, &infinity, &p_affine);
    EXPECT_TRUE(blst_p2_is_equal(&result, &p1));
}
//
//extern "C" void frobenius_map_fp2(vec384x ret, const vec384x a, size_t n);
//
//
//TEST(E2PointTest, FrobeniusMap) {
//    // Test Frobenius map properties on Fp2 field elements, not directly on points
//
//    // Start with a G2 generator point (just to get sample Fp2 values)
//    POINTonE2 p = *blst_p2_generator();
//
//    // Create sample Fp2 elements from point coordinates
//    vec384x sample_fp2;
//    vec_copy(sample_fp2[0], p.X[0], sizeof(vec384));  // Real part
//    vec_copy(sample_fp2[1], p.X[1], sizeof(vec384));  // Imaginary part
//
//    // Verify basic property of Frobenius in Fp2: π(a+bi) = a-bi
//    vec384x result1;
//    frobenius_map_fp2(result1, sample_fp2, 1);
//
//    // Check that real part is preserved and imaginary part is negated
//    bool real_preserved = vec_is_equal(result1[0], sample_fp2[0], sizeof(vec384));
//
//    // To check negation, negate the original and compare
//    vec384 negated_imag;
//    cneg_fp(negated_imag, sample_fp2[1], 1);
//    bool imag_negated = vec_is_equal(result1[1], negated_imag, sizeof(vec384));
//
//    EXPECT_TRUE(real_preserved) << "Frobenius should preserve the real part in Fp2";
//    EXPECT_TRUE(imag_negated) << "Frobenius should negate the imaginary part in Fp2";
//
//    // Apply Frobenius twice (π²) - should get back to original element for real part
//    // and negated imaginary part
//    vec384x result2;
//    frobenius_map_fp2(result2, result1, 1);
//
//    bool real_same_after_double = vec_is_equal(result2[0], sample_fp2[0], sizeof(vec384));
//    bool imag_negated_after_double = vec_is_equal(result2[1], negated_imag, sizeof(vec384));
//
//    EXPECT_TRUE(real_same_after_double)
//        << "Double Frobenius should preserve real part in Fp2";
//    EXPECT_FALSE(vec_is_equal(result2[1], sample_fp2[1], sizeof(vec384)))
//        << "Double Frobenius should not preserve imaginary part in Fp2";
//    EXPECT_TRUE(imag_negated_after_double)
//        << "Double Frobenius should negate imaginary part in Fp2";
//
//    // Subgroup membership test - G2 generators should remain in G2
//    POINTonE2_affine gen_affine;
//    blst_p2_to_affine(&gen_affine, &p);
//    EXPECT_TRUE(blst_p2_affine_in_g2(&gen_affine))
//        << "Generator should be in G2";
//
//    // Note: Testing π(P) directly on points requires proper twist constants
//    // We've verified the field-level Frobenius works correctly
//}

TEST(E2PointTest, Fp2ArithmeticWithPoints) {
    // Test that points with Fp2 coordinates behave correctly
    POINTonE2 p1, p2;
    p1 = *blst_p2_generator();
    
    // Make a copy of the generator
    vec_copy(&p2, &p1, sizeof(POINTonE2));
    
    // Negate the imaginary parts of all coordinates using proper function
    cneg_fp(p2.X[1], p2.X[1], 1);  // Negate X.Im
    cneg_fp(p2.Y[1], p2.Y[1], 1);  // Negate Y.Im
    cneg_fp(p2.Z[1], p2.Z[1], 1);  // Negate Z.Im

    // The resulting point should have conjugated Fp2 coordinates
    // It's likely not on the curve anymore, but we can check
    EXPECT_FALSE(blst_p2_is_equal(&p1, &p2)) 
        << "Conjugating Fp2 components didn't change the point";

    // Check if it's still on the curve (likely not)
    bool still_on_curve = blst_p2_on_curve(&p2);
    //printf("After Fp2 component conjugation, point is %s the curve\n",
    //    still_on_curve ? "ON" : "NOT ON");
        
    // Let's try another approach: create a point with real parts only
    POINTonE2 p3;
    vec_copy(&p3, &p1, sizeof(POINTonE2));
    
    // Zero out all imaginary parts
    vec_zero(p3.X[1], sizeof(vec384));
    vec_zero(p3.Y[1], sizeof(vec384));
    vec_zero(p3.Z[1], sizeof(vec384));
    
    // This should definitely not be on the curve or equal to original
    EXPECT_FALSE(blst_p2_is_equal(&p1, &p3)) 
        << "Zeroing imaginary parts didn't change the point";
        
    // Check if it's on the curve (unlikely)
    bool p3_on_curve = blst_p2_on_curve(&p3);
    //printf("Point with zeroed imaginary parts is %s the curve\n",
    //    p3_on_curve ? "ON" : "NOT ON");
}