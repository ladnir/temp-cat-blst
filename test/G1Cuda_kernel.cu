#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

extern "C" {
#include "blst/curve/e1.h"
#include "blst/curve/map_to_g1.h"
}

#ifndef DEVICE_FUNC
#define DEVICE_FUNC
#endif

#ifndef  GLOBAL_FUNC
#define GLOBAL_FUNC
#endif // ! GLOBAL_FUNC

DEVICE_FUNC bool cuda_point_equal(const POINTonE1* a, const POINTonE1* b) {
    return blst_p1_is_equal(a, b);
}

DEVICE_FUNC void cuda_set_point_generator(POINTonE1* p) {
    *p = *blst_p1_generator();
}

DEVICE_FUNC void cuda_set_point_infinity(POINTonE1* p) {
    memset(p, 0, sizeof(POINTonE1));
}

DEVICE_FUNC bool cuda_point_is_infinity(const POINTonE1* p) {
    return blst_p1_is_inf(p);
}

DEVICE_FUNC bool cuda_point_on_curve(const POINTonE1* p) {
    return blst_p1_on_curve(p);
}

// CUDA kernels for G1 tests

GLOBAL_FUNC void cuda_test_g1_point_addition_kernel(bool* results) {
    POINTonE1 p1, p2, result, temp;

    // Get generator point
    cuda_set_point_generator(&p1);

    // Create a different point (we'll use 2*G)
    blst_p1_double(&temp, &p1);

    // Add generator and the different point (G + 2G)
    blst_p1_add(&result, &p1, &temp);

    // Compute 3G using scalar multiplication
    unsigned char scalar[32] = { 3 };
    blst_p1_mult(&p2, &p1, scalar, 8);

    // Results should match
    results[0] = cuda_point_equal(&result, &p2);

    // Also test the doubling-specific function
    blst_p1_double(&result, &p1);

    // Create 2G with scalar multiplication
    unsigned char scalar2[32] = { 2 };
    blst_p1_mult(&p2, &p1, scalar2, 8);

    // Results should match
    results[1] = cuda_point_equal(&result, &p2);
}

GLOBAL_FUNC void cuda_test_g1_identity_operations_kernel(bool* results) {
    POINTonE1 p1, identity, result;

    // Get generator point
    cuda_set_point_generator(&p1);

    // Create identity point (point at infinity)
    cuda_set_point_infinity(&identity);

    // P + O = P
    blst_p1_add(&result, &p1, &identity);
    results[0] = cuda_point_equal(&result, &p1);

    // O + P = P
    blst_p1_add(&result, &identity, &p1);
    results[1] = cuda_point_equal(&result, &p1);

    // Check identity properties
    results[2] = cuda_point_is_infinity(&identity);
    results[3] = !cuda_point_is_infinity(&p1);
}

GLOBAL_FUNC void cuda_test_g1_point_negation_kernel(bool* results) {
    POINTonE1 p1, neg_p1, result;

    // Get generator point
    cuda_set_point_generator(&p1);

    // Make a copy of the point and negate it
    neg_p1 = p1;
    blst_p1_cneg(&neg_p1, 1);

    // P + (-P) should be the identity
    blst_p1_add(&result, &p1, &neg_p1);

    // Result should be point at infinity
    results[0] = cuda_point_is_infinity(&result);

    // Check that negation doesn't affect curve membership
    results[1] = cuda_point_on_curve(&neg_p1);
}

GLOBAL_FUNC void cuda_test_g1_scalar_multiplication_kernel(bool* results) {
    POINTonE1 p1, result, expected;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    cuda_set_point_generator(&p1);

    // Scalar = 2
    scalar[0] = 2;

    // Compute 2*P
    blst_p1_mult(&result, &p1, scalar, 8);

    // Should equal P + P (using doubling formula)
    blst_p1_double(&expected, &p1);
    results[0] = cuda_point_equal(&result, &expected);

    // Scalar = 3
    scalar[0] = 3;

    // Compute 3*P
    blst_p1_mult(&result, &p1, scalar, 8);

    // Should equal 2*P + P
    POINTonE1 doubled;
    blst_p1_double(&doubled, &p1);         // doubled = 2*P
    blst_p1_add(&expected, &doubled, &p1); // expected = 2*P + P = 3*P

    results[1] = cuda_point_equal(&result, &expected);

    // Scalar = 0
    scalar[0] = 0;
    blst_p1_mult(&result, &p1, scalar, 8);
    results[2] = cuda_point_is_infinity(&result);
}

GLOBAL_FUNC void cuda_test_g1_curve_membership_kernel(bool* results) {
    POINTonE1 p1, modified;
    POINTonE1_affine p1_affine;

    // Get generator point
    cuda_set_point_generator(&p1);

    // Generator should be on curve
    results[0] = cuda_point_on_curve(&p1);

    // Convert to affine
    blst_p1_to_affine(&p1_affine, &p1);

    // Affine point should also be on curve
    results[1] = blst_p1_affine_on_curve(&p1_affine);

    // Modify point to be off curve
    modified = p1;
    modified.X[0] ^= 1;  // Flip a bit

    // Modified point should not be on curve
    results[2] = !cuda_point_on_curve(&modified);
}

GLOBAL_FUNC void cuda_test_g1_serialization_kernel(bool* results) {
    POINTonE1 p1, p2;
    POINTonE1_affine p1_affine;
    unsigned char serialized[96];

    // Get generator point
    cuda_set_point_generator(&p1);

    // Serialize the point
    blst_p1_serialize(serialized, &p1);

    // Deserialize back to a point
    BLST_ERROR err = blst_p1_deserialize(&p1_affine, serialized);
    results[0] = (err == BLST_SUCCESS);

    if (err == BLST_SUCCESS) {
        // Convert from affine back to Jacobian coordinates
        blst_p1_from_affine(&p2, &p1_affine);

        // Points should be equal
        results[1] = cuda_point_equal(&p1, &p2);
    }
    else {
        results[1] = false;
    }
}

GLOBAL_FUNC void cuda_test_g1_compression_kernel(bool* results) {
    POINTonE1 p1;
    POINTonE1_affine p1_affine, decompressed;
    unsigned char compressed[48];

    // Get generator point
    cuda_set_point_generator(&p1);

    // Convert to affine first
    blst_p1_to_affine(&p1_affine, &p1);

    // Compress the point
    blst_p1_affine_compress(compressed, &p1_affine);

    // Verify compressed format
    results[0] = (compressed[0] & 0x80) && !(compressed[0] & 0x40);

    // Try to uncompress
    BLST_ERROR err = blst_p1_uncompress(&decompressed, compressed);
    results[1] = (err == BLST_SUCCESS);

    if (err == BLST_SUCCESS) {
        // Points should be equal
        results[2] = blst_p1_affine_is_equal(&p1_affine, &decompressed);
    }
    else {
        results[2] = false;
    }
}

GLOBAL_FUNC void cuda_test_g1_coordinate_conversion_kernel(bool* results) {
    POINTonE1 p1, p2;
    POINTonE1_affine p_affine;

    // Get generator point
    cuda_set_point_generator(&p1);

    // Convert to affine
    blst_p1_to_affine(&p_affine, &p1);

    // Convert back to Jacobian
    blst_p1_from_affine(&p2, &p_affine);

    // Should be equal to original point
    results[0] = cuda_point_equal(&p1, &p2);

    // Convert from Jacobian to Jacobian (normalized Z)
    blst_p1_from_jacobian(&p2, &p1);

    // Should be equal to original point
    results[1] = cuda_point_equal(&p1, &p2);
}

GLOBAL_FUNC void cuda_test_g1_hash_to_curve_kernel(bool* results) {
    POINTonE1 p1, p2;
    unsigned char msg[34] = "example message for hash-to-curve";
    unsigned char dst[30] = "QUUX-V01-CS02-with-BLS12381G1";

    // Test hash_to_curve
    blst_hash_to_g1(&p1, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);

    // Result should be on curve
    results[0] = cuda_point_on_curve(&p1);

    // Hash different message should give different point
    msg[0] ^= 1;  // Flip a bit
    blst_hash_to_g1(&p2, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);
    results[1] = !cuda_point_equal(&p1, &p2);
}

// Batch processing tests for G1
GLOBAL_FUNC void cuda_test_g1_batch_operations_kernel(bool* results, int batch_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= batch_size) return;

    POINTonE1 p1, result;

    // Get generator point
    cuda_set_point_generator(&p1);

    // Create different scalar values for each thread
    unsigned char scalar1[32] = { 0 };
    unsigned char scalar2[32] = { 0 };
    scalar1[0] = (idx % 7) + 1;
    scalar2[0] = ((idx + 3) % 5) + 1;

    // Perform scalar multiplications
    POINTonE1 temp1, temp2;
    blst_p1_mult(&temp1, &p1, scalar1, 8);     // scalar1 * G
    blst_p1_mult(&temp2, &p1, scalar2, 8);     // scalar2 * G

    // Add the results
    blst_p1_add(&result, &temp1, &temp2);

    // Verify result is on curve and not infinity
    bool on_curve = cuda_point_on_curve(&result);
    bool not_infinity = !cuda_point_is_infinity(&result);

    results[idx] = on_curve && not_infinity;
}

GLOBAL_FUNC void cuda_test_g1_performance_kernel(int* operations_completed) {
    POINTonE1 p1, result;
    unsigned char scalar[32] = { 0 };

    // Get generator point
    cuda_set_point_generator(&p1);

    // Set a non-trivial scalar
    scalar[0] = 5;
    scalar[1] = 7;

    int ops = 0;
    // Perform many operations to test performance
    for (int i = 0; i < 1000; i++) {  // More operations for G1 (simpler than Fp12)
        blst_p1_mult(&result, &p1, scalar, 16);
        blst_p1_double(&p1, &result);
        blst_p1_add(&result, &p1, &result);
        ops += 3;

        // Update scalar to vary operations
        scalar[0] = (scalar[0] + 1) % 256;
    }

    *operations_completed = ops;
}