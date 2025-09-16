#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

extern "C" {
#include "blst/curve/e2.h"
#include "blst/curve/map_to_g2.h"
}

#ifndef DEVICE_FUNC
#define DEVICE_FUNC
#endif

#ifndef  GLOBAL_FUNC
#define GLOBAL_FUNC
#endif // ! GLOBAL_FUNC

DEVICE_FUNC bool cuda_point_equal_g2(const POINTonE2* a, const POINTonE2* b) {
    return blst_p2_is_equal(a, b);
}

DEVICE_FUNC void cuda_set_point_generator_g2(POINTonE2* p) {
    *p = *blst_p2_generator();
}

DEVICE_FUNC void cuda_set_point_infinity_g2(POINTonE2* p) {
    memset(p, 0, sizeof(POINTonE2));
}

DEVICE_FUNC bool cuda_point_is_infinity_g2(const POINTonE2* p) {
    return blst_p2_is_inf(p);
}

DEVICE_FUNC bool cuda_point_on_curve_g2(const POINTonE2* p) {
    return blst_p2_on_curve(p);
}

// CUDA kernels for G2 tests

GLOBAL_FUNC void cuda_test_g2_point_addition_kernel(bool* results) {
    POINTonE2 p1, p2, result, temp;

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Create a different point (we'll use 2*G)
    blst_p2_double(&temp, &p1);

    // Add generator and the different point (G + 2G)
    blst_p2_add(&result, &p1, &temp);

    // Compute 3G using scalar multiplication
    unsigned char scalar[32] = { 3 };
    blst_p2_mult(&p2, &p1, scalar, 8);

    // Results should match
    results[0] = cuda_point_equal_g2(&result, &p2);

    // Also test the doubling-specific function
    blst_p2_double(&result, &p1);

    // Create 2G with scalar multiplication
    unsigned char scalar2[32] = { 2 };
    blst_p2_mult(&p2, &p1, scalar2, 8);

    // Results should match
    results[1] = cuda_point_equal_g2(&result, &p2);
}

GLOBAL_FUNC void cuda_test_g2_identity_operations_kernel(bool* results) {
    POINTonE2 p1, identity, result;

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Create identity point (point at infinity)
    cuda_set_point_infinity_g2(&identity);
    identity.Z[0][0] = 0;  // explicitly set Z to 0 for infinity

    // P + O = P
    blst_p2_add(&result, &p1, &identity);
    results[0] = cuda_point_equal_g2(&result, &p1);

    // O + P = P
    blst_p2_add(&result, &identity, &p1);
    results[1] = cuda_point_equal_g2(&result, &p1);

    // Check identity properties
    results[2] = cuda_point_is_infinity_g2(&identity);
    results[3] = !cuda_point_is_infinity_g2(&p1);
}

GLOBAL_FUNC void cuda_test_g2_point_negation_kernel(bool* results) {
    POINTonE2 p1, neg_p1, result;

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Make a copy of the point and negate it
    neg_p1 = p1;
    blst_p2_cneg(&neg_p1, 1);

    // P + (-P) should be the identity
    blst_p2_add(&result, &p1, &neg_p1);

    // Result should be point at infinity
    results[0] = cuda_point_is_infinity_g2(&result);

    // Check that negation doesn't affect curve membership
    results[1] = cuda_point_on_curve_g2(&neg_p1);
}

GLOBAL_FUNC void cuda_test_g2_scalar_multiplication_kernel(bool* results) {
    POINTonE2 p1, result, expected;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Scalar = 2
    scalar[0] = 2;

    // Compute 2*P
    blst_p2_mult(&result, &p1, scalar, 8);

    // Should equal P + P (using doubling formula)
    blst_p2_double(&expected, &p1);
    results[0] = cuda_point_equal_g2(&result, &expected);

    // Scalar = 3
    scalar[0] = 3;

    // Compute 3*P
    blst_p2_mult(&result, &p1, scalar, 8);

    // Should equal 2*P + P
    POINTonE2 doubled;
    blst_p2_double(&doubled, &p1);         // doubled = 2*P
    blst_p2_add(&expected, &doubled, &p1); // expected = 2*P + P = 3*P

    results[1] = cuda_point_equal_g2(&result, &expected);

    // Scalar = 0
    scalar[0] = 0;
    blst_p2_mult(&result, &p1, scalar, 8);
    results[2] = cuda_point_is_infinity_g2(&result);
}

GLOBAL_FUNC void cuda_test_g2_curve_membership_kernel(bool* results) {
    POINTonE2 p1, modified;
    POINTonE2_affine p1_affine;

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Generator should be on curve
    results[0] = cuda_point_on_curve_g2(&p1);

    // Convert to affine
    blst_p2_to_affine(&p1_affine, &p1);

    // Affine point should also be on curve
    results[1] = blst_p2_affine_on_curve(&p1_affine);

    // Modify point to be off curve
    modified = p1;
    modified.X[0][0] ^= 1;  // Flip a bit

    // Modified point should not be on curve
    results[2] = !cuda_point_on_curve_g2(&modified);
}

GLOBAL_FUNC void cuda_test_g2_serialization_kernel(bool* results) {
    POINTonE2 p1, p2;
    POINTonE2_affine p1_affine;
    unsigned char serialized[192];  // 2 * 96 bytes for G2 point

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Serialize the point
    blst_p2_serialize(serialized, &p1);

    // Deserialize back to a point
    BLST_ERROR err = blst_p2_deserialize(&p1_affine, serialized);
    results[0] = (err == BLST_SUCCESS);

    if (err == BLST_SUCCESS) {
        // Convert from affine back to Jacobian coordinates
        blst_p2_from_affine(&p2, &p1_affine);

        // Points should be equal
        results[1] = cuda_point_equal_g2(&p1, &p2);
    }
    else {
        results[1] = false;
    }
}

GLOBAL_FUNC void cuda_test_g2_compression_kernel(bool* results) {
    POINTonE2 p1;
    POINTonE2_affine p1_affine, decompressed;
    unsigned char compressed[96];  // 2 * 48 bytes for compressed G2

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Convert to affine first
    blst_p2_to_affine(&p1_affine, &p1);

    // Compress the point
    blst_p2_affine_compress(compressed, &p1_affine);

    // Verify compressed format
    results[0] = (compressed[0] & 0x80) && !(compressed[0] & 0x40);

    // Try to uncompress
    BLST_ERROR err = blst_p2_uncompress(&decompressed, compressed);
    results[1] = (err == BLST_SUCCESS);

    if (err == BLST_SUCCESS) {
        // Points should be equal
        results[2] = blst_p2_affine_is_equal(&p1_affine, &decompressed);
    }
    else {
        results[2] = false;
    }
}

GLOBAL_FUNC void cuda_test_g2_coordinate_conversion_kernel(bool* results) {
    POINTonE2 p1, p2;
    POINTonE2_affine p_affine;

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Convert to affine
    blst_p2_to_affine(&p_affine, &p1);

    // Convert back to Jacobian
    blst_p2_from_affine(&p2, &p_affine);

    // Should be equal to original point
    results[0] = cuda_point_equal_g2(&p1, &p2);

    // Convert from Jacobian to Jacobian (normalized Z)
    blst_p2_from_jacobian(&p2, &p1);

    // Should be equal to original point
    results[1] = cuda_point_equal_g2(&p1, &p2);
}

GLOBAL_FUNC void cuda_test_g2_hash_to_curve_kernel(bool* results) {
    POINTonE2 p1, p2;
    unsigned char msg[34] = "example message for hash-to-curve";
    unsigned char dst[30] = "QUUX-V01-CS02-with-BLS12381G2";

    // Test hash_to_curve
    blst_hash_to_g2(&p1, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);

    // Result should be on curve
    results[0] = cuda_point_on_curve_g2(&p1);

    // Hash different message should give different point
    msg[0] ^= 1;  // Flip a bit
    blst_hash_to_g2(&p2, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);
    results[1] = !cuda_point_equal_g2(&p1, &p2);
}

GLOBAL_FUNC void cuda_test_g2_mixed_addition_kernel(bool* results) {
    POINTonE2 p1, p2, result;
    POINTonE2_affine p2_affine, p_affine;

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Compute 2*P and 3*P
    blst_p2_double(&p2, &p1);
    POINTonE2 p3;
    blst_p2_add(&p3, &p2, &p1);

    // Convert 2*P to affine
    blst_p2_to_affine(&p2_affine, &p2);

    // P + (2P affine) should equal 3P
    blst_p2_add_affine(&result, &p1, &p2_affine);
    results[0] = cuda_point_equal_g2(&result, &p3);

    // Test P + P using add_or_double_affine
    blst_p2_to_affine(&p_affine, &p1);
    blst_p2_add_or_double_affine(&result, &p1, &p_affine);
    results[1] = cuda_point_equal_g2(&result, &p2);

    // Test mixed addition with infinity
    POINTonE2_affine infinity_affine;
    memset(&infinity_affine, 0, sizeof(infinity_affine));
    blst_p2_add_affine(&result, &p1, &infinity_affine);
    results[2] = cuda_point_equal_g2(&result, &p1);
}

GLOBAL_FUNC void cuda_test_g2_advanced_operations_kernel(bool* results) {
    POINTonE2 p1, neg_p1, result, temp;

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Test add_or_double with negation
    neg_p1 = p1;
    blst_p2_cneg(&neg_p1, 1);

    // P + (-P) should be infinity
    blst_p2_add_or_double(&result, &p1, &neg_p1);
    results[0] = cuda_point_is_infinity_g2(&result);

    // Test add_or_double with same point (should handle doubling)
    blst_p2_add_or_double(&result, &p1, &p1);
    blst_p2_double(&temp, &p1);
    results[1] = cuda_point_equal_g2(&result, &temp);

    // Test generator properties
    results[2] = !cuda_point_is_infinity_g2(&p1);
    results[3] = cuda_point_on_curve_g2(&p1);
}

extern "C" {
    int blst_p2_affine_in_g2(const POINTonE2_affine* p);
}

GLOBAL_FUNC void cuda_test_g2_subgroup_membership_kernel(bool* results) {
    POINTonE2 point;
    POINTonE2_affine point_affine;

    // Get generator point
    cuda_set_point_generator_g2(&point);
    blst_p2_to_affine(&point_affine, &point);

    // Generator should be in G2
    results[0] = blst_p2_affine_in_g2(&point_affine);

    // Get affine generator directly
    const POINTonE2_affine* g_affine = blst_p2_affine_generator();
    results[1] = blst_p2_affine_on_curve(g_affine);
    results[2] = !blst_p2_affine_is_inf(g_affine);
    results[3] = blst_p2_affine_in_g2(g_affine);
}

GLOBAL_FUNC void cuda_test_g2_serialization_formats_kernel(bool* results) {
    POINTonE2 p1 = *blst_p2_generator();
    POINTonE2_affine p_affine;
    blst_p2_to_affine(&p_affine, &p1);

    unsigned char compressed[96];
    unsigned char uncompressed[192];

    // Test compressed
    blst_p2_compress(compressed, &p1);
    results[0] = (compressed[0] & 0x80);  // Compressed flag

    // Test uncompressed
    blst_p2_serialize(uncompressed, &p1);
    results[1] = !(uncompressed[0] & 0x80);  // Not compressed

    // Test infinity encoding for each format
    POINTonE2 infinity;
    memset(&infinity, 0, sizeof(infinity));

    blst_p2_compress(compressed, &infinity);
    results[2] = (compressed[0] & 0x40);  // Infinity flag

    blst_p2_serialize(uncompressed, &infinity);
    results[3] = (uncompressed[0] & 0x40);  // Infinity flag
}

GLOBAL_FUNC void cuda_test_g2_large_scalar_multiplication_kernel(bool* results) {
    POINTonE2 p1, result;
    // Create a large scalar close to the order
    unsigned char large_scalar[32] = {
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x0f
    };

    p1 = *blst_p2_generator();
    blst_p2_mult(&result, &p1, large_scalar, 252);

    // Result should be on curve
    results[0] = blst_p2_on_curve(&result);
}

// Batch processing tests for G2
GLOBAL_FUNC void cuda_test_g2_batch_operations_kernel(bool* results, int batch_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= batch_size) return;

    POINTonE2 p1, result;

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Create different scalar values for each thread
    unsigned char scalar1[32] = { 0 };
    unsigned char scalar2[32] = { 0 };
    scalar1[0] = (idx % 7) + 1;
    scalar2[0] = ((idx + 3) % 5) + 1;

    // Perform scalar multiplications
    POINTonE2 temp1, temp2;
    blst_p2_mult(&temp1, &p1, scalar1, 8);     // scalar1 * G
    blst_p2_mult(&temp2, &p1, scalar2, 8);     // scalar2 * G

    // Add the results
    blst_p2_add(&result, &temp1, &temp2);

    // Verify result is on curve and not infinity
    bool on_curve = cuda_point_on_curve_g2(&result);
    bool not_infinity = !cuda_point_is_infinity_g2(&result);

    results[idx] = on_curve && not_infinity;
}

GLOBAL_FUNC void cuda_test_g2_performance_kernel(int* operations_completed) {
    POINTonE2 p1, result;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    cuda_set_point_generator_g2(&p1);

    // Set a non-trivial scalar
    scalar[0] = 5;
    scalar[1] = 7;

    int ops = 0;
    // Perform many operations to test performance (fewer than G1 since G2 is more complex)
    for (int i = 0; i < 500; i++) {
        blst_p2_mult(&result, &p1, scalar, 16);
        blst_p2_double(&p1, &result);
        blst_p2_add(&result, &p1, &result);
        ops += 3;

        // Update scalar to vary operations
        scalar[0] = (scalar[0] + 1) % 256;
    }

    *operations_completed = ops;
}