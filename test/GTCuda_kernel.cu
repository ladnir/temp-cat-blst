#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

#include "blst/point.h"
#include "blst/fields/fields.h"
#include "blst/e1.h"
#include "blst/e2.h"
#include "blst/pairing.h"
#include "blst/fp12_tower.h"
#include "blst/aggregate.h"
using namespace blst;

#ifndef DEVICE_FUNC
#define DEVICE_FUNC
#endif

#ifndef  GLOBAL_FUNC
#define GLOBAL_FUNC
#endif // ! GLOBAL_FUNC

DEVICE_FUNC bool cuda_gt_equal(const vec384fp12 a, const vec384fp12 b) {
    return blst_fp12_is_equal(a, b);
}

DEVICE_FUNC void cuda_set_gt_identity(vec384fp12 a) {
    blst_fp12_one_copy(a);
}

DEVICE_FUNC bool cuda_gt_is_identity(const vec384fp12 a) {
    return blst_fp12_is_one(a);
}

DEVICE_FUNC bool cuda_gt_in_group(const vec384fp12 a) {
    return blst_fp12_in_group(a);
}

DEVICE_FUNC void cuda_compute_gt_element(vec384fp12 gt_element) {
    // Compute basic GT element using generators
    POINTonE1 g1_point = *blst_p1_generator();
    POINTonE2 g2_point = *blst_p2_generator();

    POINTonE1_affine g1_affine;
    POINTonE2_affine g2_affine;

    blst_p1_to_affine(&g1_affine, &g1_point);
    blst_p2_to_affine(&g2_affine, &g2_point);

    blst_miller_loop(gt_element, &g2_affine, &g1_affine);
    blst_final_exp(gt_element, gt_element);
}

// CUDA kernels for GT tests

GLOBAL_FUNC void cuda_test_gt_identity_kernel(bool* results) {
    vec384fp12 identity, gt_element, temp;

    // Test identity element
    cuda_set_gt_identity(identity);
    results[0] = cuda_gt_is_identity(identity);

    // Test multiplying by identity
    cuda_compute_gt_element(gt_element);
    blst_fp12_mul(temp, gt_element, identity);
    results[1] = cuda_gt_equal(temp, gt_element);

    // Test that gt_element is not identity
    results[2] = !cuda_gt_is_identity(gt_element);

    // Test identity is in group
    results[3] = cuda_gt_in_group(identity);
}

GLOBAL_FUNC void cuda_test_gt_multiplication_kernel(bool* results) {
    vec384fp12 gt1, gt2, product1, product2;
    POINTonE1 g1_point, g1_double;
    POINTonE2 g2_point;
    POINTonE1_affine g1_affine, g1_double_affine;
    POINTonE2_affine g2_affine;

    // Create two different GT elements
    g1_point = *blst_p1_generator();
    g2_point = *blst_p2_generator();

    blst_p1_double(&g1_double, &g1_point);

    blst_p1_to_affine(&g1_affine, &g1_point);
    blst_p1_to_affine(&g1_double_affine, &g1_double);
    blst_p2_to_affine(&g2_affine, &g2_point);

    // gt1 = e(P, Q), gt2 = e(2P, Q)
    blst_miller_loop(gt1, &g2_affine, &g1_affine);
    blst_final_exp(gt1, gt1);

    blst_miller_loop(gt2, &g2_affine, &g1_double_affine);
    blst_final_exp(gt2, gt2);

    // Test multiplication
    blst_fp12_mul(product1, gt1, gt2);
    results[0] = !cuda_gt_equal(product1, gt1);
    results[1] = !cuda_gt_equal(product1, gt2);

    // Test commutativity: a * b = b * a
    blst_fp12_mul(product2, gt2, gt1);
    results[2] = cuda_gt_equal(product1, product2);

    // Test that products are in group
    results[3] = cuda_gt_in_group(product1);
}

GLOBAL_FUNC void cuda_test_gt_inversion_kernel(bool* results) {
    vec384fp12 gt_element, inverse, product, identity;

    // Create GT element
    cuda_compute_gt_element(gt_element);

    // Compute inverse
    blst_fp12_inverse(inverse, gt_element);

    // Test that element * inverse = identity
    blst_fp12_mul(product, gt_element, inverse);
    cuda_set_gt_identity(identity);
    results[0] = cuda_gt_equal(product, identity);

    // Test that inverse is in group
    results[1] = cuda_gt_in_group(inverse);

    // Test that inverse of identity is identity
    vec384fp12 identity_inv;
    blst_fp12_inverse(identity_inv, identity);
    results[2] = cuda_gt_equal(identity, identity_inv);
}

GLOBAL_FUNC void cuda_test_gt_squaring_kernel(bool* results) {
    vec384fp12 gt_element, square1, square2;

    // Create GT element
    cuda_compute_gt_element(gt_element);

    // Test that squaring equals multiplying by self
    blst_fp12_sqr(square1, gt_element);
    blst_fp12_mul(square2, gt_element, gt_element);
    results[0] = cuda_gt_equal(square1, square2);

    // Test that square is in group
    results[1] = cuda_gt_in_group(square1);

    // Test cyclotomic squaring (for GT elements)
    vec384fp12 cyclo_square;
    cyclotomic_sqr_fp12(cyclo_square, gt_element);
    results[2] = cuda_gt_equal(cyclo_square, square1);
}

GLOBAL_FUNC void cuda_test_gt_exponentiation_kernel(bool* results) {
    vec384fp12 gt_element, result, expected;

    cuda_compute_gt_element(gt_element);

    // Test e^0 = 1
    limb_t zero[4] = { 0, 0, 0, 0 };
    blst_fp12_pow(result, gt_element, zero, sizeof(zero) * 8);
    vec384fp12 identity;
    cuda_set_gt_identity(identity);
    results[0] = cuda_gt_equal(result, identity);

    // Test e^1 = e
    limb_t one[4] = { 1, 0, 0, 0 };
    blst_fp12_pow(result, gt_element, one, sizeof(one) * 8);
    results[1] = cuda_gt_equal(result, gt_element);

    // Test e^2 = e * e
    limb_t two[4] = { 2, 0, 0, 0 };
    blst_fp12_pow(result, gt_element, two, sizeof(two) * 8);
    blst_fp12_sqr(expected, gt_element);
    results[2] = cuda_gt_equal(result, expected);

    // Test that exponentiation result is in group
    results[3] = cuda_gt_in_group(result);
}

GLOBAL_FUNC void cuda_test_gt_bilinearity_kernel(bool* results) {
    POINTonE1 g1_point, g1_double;
    POINTonE2 g2_point, g2_double;
    POINTonE1_affine g1_affine, g1_double_affine;
    POINTonE2_affine g2_affine, g2_double_affine;
    vec384fp12 e_2P_Q, e_P_2Q, e_PQ_squared, base_pairing;

    // Setup points
    g1_point = *blst_p1_generator();
    g2_point = *blst_p2_generator();

    blst_p1_double(&g1_double, &g1_point);
    blst_p2_double(&g2_double, &g2_point);

    blst_p1_to_affine(&g1_affine, &g1_point);
    blst_p1_to_affine(&g1_double_affine, &g1_double);
    blst_p2_to_affine(&g2_affine, &g2_point);
    blst_p2_to_affine(&g2_double_affine, &g2_double);

    // Base pairing e(P, Q)
    blst_miller_loop(base_pairing, &g2_affine, &g1_affine);
    blst_final_exp(base_pairing, base_pairing);

    // e(2P, Q)
    blst_miller_loop(e_2P_Q, &g2_affine, &g1_double_affine);
    blst_final_exp(e_2P_Q, e_2P_Q);

    // e(P, 2Q)
    blst_miller_loop(e_P_2Q, &g2_double_affine, &g1_affine);
    blst_final_exp(e_P_2Q, e_P_2Q);

    // e(P,Q)^2
    blst_fp12_sqr(e_PQ_squared, base_pairing);

    // Test bilinearity: e(2P, Q) = e(P, 2Q) = e(P,Q)^2
    results[0] = cuda_gt_equal(e_2P_Q, e_P_2Q);
    results[1] = cuda_gt_equal(e_2P_Q, e_PQ_squared);
    results[2] = cuda_gt_equal(e_P_2Q, e_PQ_squared);
}

GLOBAL_FUNC void cuda_test_gt_group_membership_kernel(bool* results) {
    vec384fp12 gt_element, identity, not_finalized;
    POINTonE1_affine g1_affine;
    POINTonE2_affine g2_affine;

    // Create GT element
    cuda_compute_gt_element(gt_element);
    results[0] = cuda_gt_in_group(gt_element);

    // Identity should be in group
    cuda_set_gt_identity(identity);
    results[1] = cuda_gt_in_group(identity);

    // Element before final exponentiation should not be in group
    POINTonE1 g1_point = *blst_p1_generator();
    POINTonE2 g2_point = *blst_p2_generator();
    blst_p1_to_affine(&g1_affine, &g1_point);
    blst_p2_to_affine(&g2_affine, &g2_point);

    blst_miller_loop(not_finalized, &g2_affine, &g1_affine);
    results[2] = !cuda_gt_in_group(not_finalized);

    // Operations should preserve group membership
    vec384fp12 squared;
    blst_fp12_sqr(squared, gt_element);
    results[3] = cuda_gt_in_group(squared);
}

GLOBAL_FUNC void cuda_test_gt_conjugate_kernel(bool* results) {
    vec384fp12 gt_element, conjugate, double_conjugate;

    cuda_compute_gt_element(gt_element);

    // Test conjugate
    memcpy(conjugate, gt_element, sizeof(vec384fp12));
    blst_fp12_conjugate(conjugate);

    // Conjugate of conjugate should give original
    memcpy(double_conjugate, conjugate, sizeof(vec384fp12));
    blst_fp12_conjugate(double_conjugate);
    results[0] = cuda_gt_equal(gt_element, double_conjugate);

    // Conjugate should be in GT group
    results[1] = cuda_gt_in_group(conjugate);

    // Test that conjugate is different from original (unless special case)
    results[2] = !cuda_gt_equal(gt_element, conjugate) || cuda_gt_is_identity(gt_element);
}

GLOBAL_FUNC void cuda_test_gt_serialization_kernel(bool* results) {
    vec384fp12 original, deserialized;
    unsigned char serialized[12 * 48]; // Fp12 consists of 12 base field elements

    // Create GT element
    cuda_compute_gt_element(original);

    // Serialize to bytes
    blst_bendian_from_fp12(serialized, original);

    // Deserialize back
    blst_fp12_from_bendian(deserialized, serialized);

    // Should be equal
    results[0] = cuda_gt_equal(original, deserialized);

    // Verify serialized data is non-trivial
    bool has_nonzero = false;
    for (size_t i = 0; i < sizeof(serialized); i++) {
        if (serialized[i] != 0) {
            has_nonzero = true;
            break;
        }
    }
    results[1] = has_nonzero;
}

GLOBAL_FUNC void cuda_test_gt_multi_pairing_kernel(bool* results) {
    POINTonE1 g1_point, g1_double;
    POINTonE2 g2_point;
    POINTonE1_affine g1_points[2];
    POINTonE2_affine g2_points[2];
    vec384fp12 individual_product, multi_result, e1, e2;

    // Setup points
    g1_point = *blst_p1_generator();
    g2_point = *blst_p2_generator();
    blst_p1_double(&g1_double, &g1_point);

    blst_p1_to_affine(&g1_points[0], &g1_point);
    blst_p1_to_affine(&g1_points[1], &g1_double);
    blst_p2_to_affine(&g2_points[0], &g2_point);
    blst_p2_to_affine(&g2_points[1], &g2_point);

    // Calculate individual pairings and their product
    blst_miller_loop(e1, &g2_points[0], &g1_points[0]);
    blst_miller_loop(e2, &g2_points[1], &g1_points[1]);
    blst_fp12_mul(individual_product, e1, e2);
    blst_final_exp(individual_product, individual_product);

    // Calculate using multi-pairing
    const POINTonE1_affine* g1_ptrs[3] = { &g1_points[0], &g1_points[1], nullptr };
    const POINTonE2_affine* g2_ptrs[3] = { &g2_points[0], &g2_points[1], nullptr };

    blst_miller_loop_n(multi_result, g2_ptrs, g1_ptrs, 2);
    blst_final_exp(multi_result, multi_result);

    // Should match
    results[0] = cuda_gt_equal(individual_product, multi_result);
    results[1] = cuda_gt_in_group(multi_result);
}

GLOBAL_FUNC void cuda_test_gt_final_verify_kernel(bool* results) {
    vec384fp12 gt1, gt2, identity;
    POINTonE1 g1_point, g1_double;
    POINTonE2 g2_point;
    POINTonE1_affine g1_affine, g1_double_affine;
    POINTonE2_affine g2_affine;

    // Create two different GT elements
    g1_point = *blst_p1_generator();
    g2_point = *blst_p2_generator();
    blst_p1_double(&g1_double, &g1_point);

    blst_p1_to_affine(&g1_affine, &g1_point);
    blst_p1_to_affine(&g1_double_affine, &g1_double);
    blst_p2_to_affine(&g2_affine, &g2_point);

    blst_miller_loop(gt1, &g2_affine, &g1_affine);
    blst_final_exp(gt1, gt1);

    blst_miller_loop(gt2, &g2_affine, &g1_double_affine);
    blst_final_exp(gt2, gt2);

    // Test finalverify with different elements
    results[0] = !blst_fp12_finalverify(gt1, gt2);

    // Test finalverify with same element
    results[1] = blst_fp12_finalverify(gt1, gt1);

    // Test with identity
    cuda_set_gt_identity(identity);
    results[2] = blst_fp12_finalverify(identity, identity);
}

// Batch processing tests for GT
GLOBAL_FUNC void cuda_test_gt_batch_operations_kernel(bool* results, int batch_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= batch_size) return;

    vec384fp12 gt_base, result;
    POINTonE1 g1_point;
    POINTonE2 g2_point;
    POINTonE1_affine g1_affine;
    POINTonE2_affine g2_affine;

    // Create base GT element
    g1_point = *blst_p1_generator();
    g2_point = *blst_p2_generator();

    // Scale the G1 point based on thread index
    for (int i = 0; i < (idx % 5) + 1; i++) {
        POINTonE1 temp;
        blst_p1_add(&temp, &g1_point, blst_p1_generator());
        g1_point = temp;
    }

    blst_p1_to_affine(&g1_affine, &g1_point);
    blst_p2_to_affine(&g2_affine, &g2_point);

    // Compute pairing
    blst_miller_loop(gt_base, &g2_affine, &g1_affine);
    blst_final_exp(gt_base, gt_base);

    // Perform operations
    vec384fp12 temp1, temp2;
    blst_fp12_sqr(temp1, gt_base);              // gt_base^2
    blst_fp12_mul(temp2, gt_base, temp1);       // gt_base^3
    blst_fp12_inverse(result, temp2);           // 1/gt_base^3

    // Verify result is in group and not identity
    bool in_group = cuda_gt_in_group(result);
    bool not_identity = !cuda_gt_is_identity(result);

    results[idx] = in_group && not_identity;
}

GLOBAL_FUNC void cuda_test_gt_performance_kernel(int* operations_completed) {
    vec384fp12 gt_element, result;

    // Create GT element
    cuda_compute_gt_element(gt_element);

    int ops = 0;
    // Perform many operations to test performance (fewer than basic field ops due to complexity)
    for (int i = 0; i < 50; i++) {
        blst_fp12_sqr(result, gt_element);
        blst_fp12_mul(gt_element, result, gt_element);
        blst_fp12_conjugate(result);
        memcpy(result, gt_element, sizeof(vec384fp12));
        blst_fp12_conjugate(result);
        ops += 4;
    }

    *operations_completed = ops;
}

GLOBAL_FUNC void cuda_test_gt_powers_of_two_kernel(bool* results) {
    vec384fp12 gt_element, manual, result;

    cuda_compute_gt_element(gt_element);
    memcpy(manual, gt_element, sizeof(vec384fp12));

    // Test powers of 2: g^(2^i) for small i
    const int max_power = 10; // Test 2^0, 2^1, ..., 2^9

    for (int i = 0; i < max_power; i++) {
        // Create exponent 2^i
        limb_t exp[2];
        memset(exp, 0, sizeof(exp));

        if (i < sizeof(limb_t) * 8) {
            exp[0] = (limb_t)1 << i;
        }
        else {
            exp[1] = (limb_t)1 << (i - sizeof(limb_t) * 8);
        }

        // Calculate using blst_fp12_pow
        blst_fp12_pow(result, gt_element, exp, sizeof(exp) * 8);

        // They should be equal
        if (!cuda_gt_equal(result, manual)) {
            results[0] = false;
            return;
        }

        // Prepare for next iteration by squaring
        vec384fp12 temp;
        blst_fp12_sqr(temp, manual);
        memcpy(manual, temp, sizeof(vec384fp12));
    }

    results[0] = true;
}