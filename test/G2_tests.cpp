#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "blst/curve/point.h"

#include "blst/curve/map_to_G.h"

using namespace blst;

// Helper: Print G2 point coordinates for debugging
void print_point_e2(const G2Point& p, const char* name) {
    print(name, p);
}

// Helper: Convert hex string to bytes
//void hex_to_bytes_e2(unsigned char* bytes, const char* hex) {
//    size_t len = strlen(hex) / 2;
//    for (size_t i = 0; i < len; i++) {
//        // unsafe
//        //sscanf(hex + (i * 2), "%2hhx", &bytes[i]);
//		
//
//    }
//}

TEST(E2PointTest, PointAddition) {
    G2Point p1, p2, result, temp;

    // Get generator point
    p1 = G2Point::generator();

    // Create a different point (we'll use 2*G)
    G2Point::doubling(temp, p1);

    // Add generator and the different point (G + 2G)
    G2Point::add(result, p1, temp);

    // Compute 3G using scalar multiplication
    unsigned char scalar[32] = { 3 };
    G2Point::mult(p2, p1, scalar, 8);

    // Results should match
    EXPECT_TRUE(G2Point::equal(result, p2));

    // Also test the doubling-specific function
    G2Point::doubling(result, p1);

    // Create 2G with scalar multiplication
    unsigned char scalar2[32] = { 2 };
    G2Point::mult(p2, p1, scalar2, 8);

    // Results should match
    EXPECT_TRUE(G2Point::equal(result, p2));
}

TEST(E2PointTest, AdditionWithIdentity) {
    G2Point p1, identity, result;

    // Get generator point
    p1 = G2Point::generator();

    // Create identity point (point at infinity)
    identity = G2Point::identity();

    // P + O = P
    G2Point::add(result, p1, identity);
    EXPECT_TRUE(G2Point::equal(result, p1));
}

TEST(E2PointTest, PointNegation) {
    G2Point p1, neg_p1, result;

    // Get generator point
    p1 = G2Point::generator();

    // Make a copy of the point and negate it
    neg_p1 = p1;
    G2Point::cneg(neg_p1, 1);

    // P + (-P) should be the identity
    G2Point::add(result, p1, neg_p1);

    // Result should be point at infinity
    EXPECT_TRUE(G2Point::is_inf(result));
}

TEST(E2PointTest, ScalarMultiplication) {
    G2Point p1, result, expected;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    p1 = G2Point::generator();

    // Scalar = 2
    scalar[0] = 2;

    // Compute 2*P
    G2Point::mult(result, p1, scalar, 8);

    // Should equal P + P (using doubling formula)
    G2Point::doubling(expected, p1);
    EXPECT_TRUE(G2Point::equal(result, expected));

    // Scalar = 3
    scalar[0] = 3;

    // Compute 3*P
    G2Point::mult(result, p1, scalar, 8);

    // Should equal 2*P + P
    G2Point doubled;
    G2Point::doubling(doubled, p1);         // doubled = 2*P
    G2Point::add(expected, doubled, p1); // expected = 2*P + P = 3*P

    EXPECT_TRUE(G2Point::equal(result, expected));
}

TEST(E2PointTest, ScalarMultiplicationByZero) {
    G2Point p1, result;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    p1 = G2Point::generator();

    // Compute 0*P
    G2Point::mult(result, p1, scalar, 8);

    // Result should be point at infinity
    EXPECT_TRUE(G2Point::is_inf(result));
}

TEST(E2PointTest, PointOnCurve) {
    G2Point p1;
    G2AffinePoint p1_affine;

    // Get generator point
    p1 = G2Point::generator();

    // Generator should be on curve
    EXPECT_TRUE(G2Point::on_curve(p1));

    // Convert to affine
    G2Point::to_affine(p1_affine, p1);
    
    // Affine point should also be on curve
    EXPECT_TRUE(G2AffinePoint::on_curve(p1_affine));

#if !BLST_DEFAULT_AVX
    // Modify point to be off curve
    p1.X[0][0] ^= 1;  // Flip a bit

    // Modified point should not be on curve
    EXPECT_FALSE(G2Point::on_curve(p1));
#endif
}

TEST(E2PointTest, SerializeDeserialize) {
    G2Point p1, p2;
    G2AffinePoint p1_affine;
    unsigned char serialized[192];  // 2 * 96 bytes for G2 point

    // Get generator point
    p1 = G2Point::generator();

    // Serialize the point
    G2Point::serialize(serialized, p1);

    // Deserialize back to a point
    EXPECT_EQ(G2AffinePoint::deserialize(p1_affine, serialized), BLST_SUCCESS);

    // Convert from affine back to Jacobian coordinates
    G2Point::from_affine(p2, p1_affine);

    // Points should be equal
    EXPECT_TRUE(G2Point::equal(p1, p2));
}

TEST(E2PointTest, BasicCompressUncompress) {
    // Starting with the generator point in projective form
    G2Point point = G2Point::generator();

    // Convert to affine
    G2AffinePoint affine_point;
    G2Point::to_affine(affine_point, point);

    // First confirm it's on the curve
    ASSERT_TRUE(G2AffinePoint::on_curve(affine_point)) << "Point not on curve";

    // Compress the affine point
    unsigned char compressed[96];  // 2 * 48 bytes for compressed G2
    G2AffinePoint::compress(compressed, affine_point);

    // Uncompress
    G2AffinePoint result;
    BLST_ERROR err = G2AffinePoint::uncompress(result, compressed);
    ASSERT_EQ(err, BLST_SUCCESS) << "Decompression failed with code " << err;

    // Points should match
    EXPECT_TRUE(G2AffinePoint::equal(result, affine_point))
        << "Points don't match after compression/decompression";
}

TEST(E2PointTest, AffineCompressUncompress) {
    G2AffinePoint original, decompressed;
    unsigned char compressed[96] = { 0 };

    // Get the generator in affine form
    original = G2AffinePoint::generator();

    // Verify it's on curve
    ASSERT_TRUE(G2AffinePoint::on_curve(original)) << "Generator point not on curve";

    // Compress 
    G2AffinePoint::compress(compressed, original);

    // Verify compressed format
    EXPECT_TRUE(compressed[0] & 0x80) << "Compressed bit not set";
    EXPECT_FALSE(compressed[0] & 0x40) << "Infinity bit incorrectly set";

    // Uncompress
    BLST_ERROR err = G2AffinePoint::uncompress(decompressed, compressed);
    ASSERT_EQ(err, BLST_SUCCESS) << "Decompression failed with code " << err;

    // Verify the point is on curve
    ASSERT_TRUE(G2AffinePoint::on_curve(decompressed)) << "Decompressed point not on curve";

    // Directly compare the affine coordinates
    EXPECT_TRUE(G2AffinePoint::equal(original, decompressed))
        << "Affine points don't match after compression cycle";
}

TEST(E2PointTest, CompressUncompress) {
    G2Point p1, p2;
    G2AffinePoint p1_affine, p2_affine;
    unsigned char compressed[96] = { 0 };

    // Start with the generator point which is guaranteed to be valid
    p1 = G2Point::generator();

    // Convert to affine first (important for compression)
    G2Point::to_affine(p1_affine, p1);

    // Compress directly
    G2Point::compress(compressed, p1);

    // Manually verify the format - first byte should have compressed bit set (0x80)
    // but not have infinity bit set (0x40)
    EXPECT_TRUE((compressed[0] & 0x80) && !(compressed[0] & 0x40));

    // Try to uncompress
    BLST_ERROR err = G2AffinePoint::uncompress(p2_affine, compressed);
    ASSERT_EQ(err, BLST_SUCCESS) << "Point decompression failed";

    // Convert back to Jacobian coordinates
    G2Point::from_affine(p2, p2_affine);

    // Points should be equal
    EXPECT_TRUE(G2Point::equal(p1, p2));
}

TEST(E2PointTest, InfinityPoint) {
    G2Point identity;
    unsigned char serialized[192];
    G2AffinePoint affine_identity;

    // Create identity point
    memset(&identity, 0, sizeof(identity));

    // Point should be at infinity
    EXPECT_TRUE(G2Point::is_inf(identity));

    // Serialize infinity point
    G2Point::serialize(serialized, identity);

    // First byte should have infinity bit set
    EXPECT_TRUE(serialized[0] & 0x40);

    // Deserialize back to a point
    EXPECT_EQ(G2AffinePoint::deserialize(affine_identity, serialized), BLST_SUCCESS);

    // Should be affine infinity
    EXPECT_TRUE(G2AffinePoint::is_inf(affine_identity));
}

TEST(E2PointTest, SigningBasicOperation) {
    G2Point p1, sig;
    unsigned char SK[32] = {
        0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08,
        0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10,
        0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18,
        0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20
    };
    unsigned char pk_bytes[192];

    // Generate public key from secret key
    G2Point::to_pk2(pk_bytes, NULL, SK);

    // Create a message point - we'll use the generator for simplicity
    p1 = G2Point::generator();

    // Sign the message
    G2Point::sign(sig, p1, SK);

    // Signature should be on the curve
    EXPECT_TRUE(G2Point::on_curve(sig));
}

TEST(E2PointTest, PointCoordinateConversion) {
    G2Point p1, p2;
    G2AffinePoint p_affine;

    // Get generator point
    p1 = G2Point::generator();

    // Convert to affine
    G2Point::to_affine(p_affine, p1);

    // Convert back to Jacobian
    G2Point::from_affine(p2, p_affine);

    // Should be equal to original point
    EXPECT_TRUE(G2Point::equal(p1, p2));

    // Convert from Jacobian to Jacobian (normalized Z)
    G2Point::from_jacobian(p2, p1);

    // Should be equal to original point
    EXPECT_TRUE(G2Point::equal(p1, p2));
}

TEST(E2PointTest, GeneratorProperties) {
    // Get generator point
    const G2Point g = G2Point::generator();

    // Generator should be on curve
    EXPECT_TRUE(G2Point::on_curve(g));

    // Generator should not be at infinity
    EXPECT_FALSE(G2Point::is_inf(g));

    // Get affine generator
    const G2AffinePoint g_affine = G2AffinePoint::generator();

    // Affine generator should be on curve
    EXPECT_TRUE(G2AffinePoint::on_curve(g_affine));

    // Affine generator should not be at infinity
    EXPECT_FALSE(G2AffinePoint::is_inf(g_affine));
}

TEST(E2PointTest, LargeScalarMultiplication) {
    G2Point p1, result;
    // Create a large scalar close to the order
    unsigned char large_scalar[32] = {
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x0f
    };

    // Get generator point
    p1 = G2Point::generator();

    // Perform large scalar multiplication
    G2Point::mult(result, p1, large_scalar, 252);

    // Result should be on curve
    EXPECT_TRUE(G2Point::on_curve(result));
}

TEST(E2PointTest, AffineEquality) {
    G2AffinePoint a1, a2;
    G2Point p1, p2;

    // Get generator point
    p1 = G2Point::generator();

    // Double the generator
    G2Point::doubling(p2, p1);

    // Convert both to affine
    G2Point::to_affine(a1, p1);
    G2Point::to_affine(a2, p1);

    // Affine points should be equal
    EXPECT_TRUE(G2AffinePoint::equal(a1, a2));

    // Convert p2 to affine
    G2Point::to_affine(a2, p2);

    // Should not be equal to a1
    EXPECT_FALSE(G2AffinePoint::equal(a1, a2));
}

TEST(E2PointTest, DeserializeInvalidPoint) {
    G2AffinePoint result;
    unsigned char invalid_data[192];

    // Fill with invalid data
    memset(invalid_data, 0xFF, sizeof(invalid_data));

    // Should fail to deserialize
    EXPECT_NE(G2AffinePoint::uncompress(result, invalid_data), BLST_SUCCESS);
}

TEST(E2PointTest, MultiplePointAddition) {
    G2Point p1, p3a, p2, result;

    // Get generator point
    p1 = G2Point::generator();

    // Compute 3*P using additions
    G2Point::doubling(p2, p1);    // temp = 2*P
    G2Point::add(p3a, p2, p1);  // p3a = 3*P

    // Compute 3*P using scalar multiplication
    unsigned char scalar[32] = { 3 };
    G2Point::mult(result, p1, scalar, 8);

    // Results should match
    EXPECT_TRUE(G2Point::equal(result, p3a));

    // Test mixed addition with different points
    G2AffinePoint p2_affine;
    G2Point::to_affine(p2_affine, p2);  // Convert 2*P to affine

    // P + (2P affine) should equal 3P
    G2Point::add(result, p1, p2_affine);
    EXPECT_TRUE(G2Point::equal(result, p3a));  // Should equal 3*P

    // Also verify that dadd works with same point (P+P)
    G2AffinePoint p_affine;
    G2Point::to_affine(p_affine, p1);

    // Test P + P using dadd
    G2Point::dadd(result, p1, p_affine);
    EXPECT_TRUE(G2Point::equal(result, p2));  // Should equal 2*P
}

TEST(E2PointTest, AdditionEdgeCases) {
    G2Point p1, neg_p1, result, temp;

    // Get generator point
    p1 = G2Point::generator();

    // Make a copy and negate it
    neg_p1 = p1;
    G2Point::cneg(neg_p1, 1);

    // P + (-P) should be infinity
    G2Point::dadd(result, p1, neg_p1);
    EXPECT_TRUE(G2Point::is_inf(result));

    // Add two identical points with regular add (should handle doubling correctly)
    G2Point::dadd(result, p1, p1);
    G2Point::doubling(temp, p1);
    EXPECT_TRUE(G2Point::equal(result, temp));
}

TEST(E2PointTest, HashToCurve) {
    G2Point p1;
    unsigned char msg[34] = "example message for hash-to-curve";
    unsigned char dst[30] = "QUUX-V01-CS02-with-BLS12381G2";

    // Test hash_to_curve
    Hash_to_curve<G2Point>(p1, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);

    // Result should be on curve
    EXPECT_TRUE(G2Point::on_curve(p1));

    // Hash different message should give different point
    G2Point p2;
    msg[0] ^= 1;  // Flip a bit
    Hash_to_curve<G2Point>(p2, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);
    EXPECT_FALSE(G2Point::equal(p1, p2));
}

TEST(E2PointTest, InvalidPointHandling) {
    // This test requires internal functions to create a point on E2 but not in G2
    // We can simulate by hashing to curve without cofactor clearing

    // For now, just verify that the generator is in G2
    G2Point point = G2Point::generator();
    G2AffinePoint point_affine;
    G2Point::to_affine(point_affine, point);

    // Generator should be in G2
    EXPECT_TRUE(G2AffinePoint::in_subgroup(point_affine));
}

TEST(E2PointTest, SerializationFormats) {
    G2Point p1 = G2Point::generator();
    G2AffinePoint p_affine;
    G2Point::to_affine(p_affine, p1);

    // Test all serialization formats
    unsigned char compressed[96];
    unsigned char uncompressed[192];

    // Test compressed
    G2Point::compress(compressed, p1);
    EXPECT_TRUE(compressed[0] & 0x80);  // Compressed flag

    // Test uncompressed
    G2Point::serialize(uncompressed, p1);
    EXPECT_FALSE(uncompressed[0] & 0x80);  // Not compressed

    // Test infinity encoding for each format
    G2Point infinity;
    memset(&infinity, 0, sizeof(infinity));

    G2Point::compress(compressed, infinity);
    EXPECT_TRUE(compressed[0] & 0x40);  // Infinity flag

    G2Point::serialize(uncompressed, infinity);
    EXPECT_TRUE(uncompressed[0] & 0x40);  // Infinity flag
}

TEST(E2PointTest, MixedAdditionCornerCases) {
    G2Point p1, result;
    G2AffinePoint infinity_affine;

    // Create affine point at infinity
    memset(&infinity_affine, 0, sizeof(infinity_affine));

    // Get generator point
    p1 = G2Point::generator();

    // P + O_affine = P
    G2Point::add(result, p1, infinity_affine);
    EXPECT_TRUE(G2Point::equal(result, p1));

    // O + P_affine = P
    G2Point infinity;
    memset(&infinity, 0, sizeof(infinity));
    G2AffinePoint p_affine;
    G2Point::to_affine(p_affine, p1);

    G2Point::add(result, infinity, p_affine);
    EXPECT_TRUE(G2Point::equal(result, p1));
}
//
//extern "C" void frobenius_map_fp2(Fp2 ret, const Fp2 a, size_t n);
//
//
//TEST(E2PointTest, FrobeniusMap) {
//    // Test Frobenius map properties on Fp2 field elements, not directly on points
//
//    // Start with a G2 generator point (just to get sample Fp2 values)
//    G2Point p = G2Point::generator();
//
//    // Create sample Fp2 elements from point coordinates
//    Fp2 sample_fp2;
//    vec_copy(sample_fp2[0], p.X[0], sizeof(Fp));  // Real part
//    vec_copy(sample_fp2[1], p.X[1], sizeof(Fp));  // Imaginary part
//
//    // Verify basic property of Frobenius in Fp2: π(a+bi) = a-bi
//    Fp2 result1;
//    frobenius_map_fp2(result1, sample_fp2, 1);
//
//    // Check that real part is preserved and imaginary part is negated
//    bool real_preserved = vec_is_equal(result1[0], sample_fp2[0], sizeof(Fp));
//
//    // To check negation, negate the original and compare
//    Fp negated_imag;
//    cFp::neg(negated_imag, sample_fp2[1], 1);
//    bool imag_negated = vec_is_equal(result1[1], negated_imag, sizeof(Fp));
//
//    EXPECT_TRUE(real_preserved) << "Frobenius should preserve the real part in Fp2";
//    EXPECT_TRUE(imag_negated) << "Frobenius should negate the imaginary part in Fp2";
//
//    // Apply Frobenius twice (π²) - should get back to original element for real part
//    // and negated imaginary part
//    Fp2 result2;
//    frobenius_map_fp2(result2, result1, 1);
//
//    bool real_same_after_double = vec_is_equal(result2[0], sample_fp2[0], sizeof(Fp));
//    bool imag_negated_after_double = vec_is_equal(result2[1], negated_imag, sizeof(Fp));
//
//    EXPECT_TRUE(real_same_after_double)
//        << "Double Frobenius should preserve real part in Fp2";
//    EXPECT_FALSE(vec_is_equal(result2[1], sample_fp2[1], sizeof(Fp)))
//        << "Double Frobenius should not preserve imaginary part in Fp2";
//    EXPECT_TRUE(imag_negated_after_double)
//        << "Double Frobenius should negate imaginary part in Fp2";
//
//    // Subgroup membership test - G2 generators should remain in G2
//    G2AffinePoint gen_affine;
//    blst_p2_to_affine(gen_affine, p);
//    EXPECT_TRUE(blst_p2_affine_in_g2(gen_affine))
//        << "Generator should be in G2";
//
//    // Note: Testing π(P) directly on points requires proper twist constants
//    // We've verified the field-level Frobenius works correctly
//}

#if !BLST_DEFAULT_AVX
TEST(E2PointTest, Fp2ArithmeticWithPoints) {
    // Test that points with Fp2 coordinates behave correctly
    G2Point p1, p2;
    p1 = G2Point::generator();
    
    // Make a copy of the generator
    vec_copy(&p2, &p1, sizeof(G2Point));
    
    // Negate the imaginary parts of all coordinates using proper function
    Fp::cneg(p2.X[1], p2.X[1], 1);  // Negate X.Im
    Fp::cneg(p2.Y[1], p2.Y[1], 1);  // Negate Y.Im
    Fp::cneg(p2.Z[1], p2.Z[1], 1);  // Negate Z.Im

    // The resulting point should have conjugated Fp2 coordinates
    // It's likely not on the curve anymore, but we can check
    EXPECT_FALSE(G2Point::equal(p1, p2)) 
        << "Conjugating Fp2 components didn't change the point";

    // Check if it's still on the curve (likely not)
    EXPECT_FALSE(G2Point::on_curve(p2));
        
    // Let's try another approach: create a point with real parts only
    G2Point p3;
    vec_copy(&p3, &p1, sizeof(G2Point));
    
    // Zero out all imaginary parts
    vec_zero(p3.X[1].data, sizeof(Fp));
    vec_zero(p3.Y[1].data, sizeof(Fp));
    vec_zero(p3.Z[1].data, sizeof(Fp));
    
    // This should definitely not be on the curve or equal to original
    EXPECT_FALSE(G2Point::equal(p1, p3)) 
        << "Zeroing imaginary parts didn't change the point";
        
    // Check if it's on the curve (unlikely)
    EXPECT_FALSE(G2Point::on_curve(p3));
}
#endif
