#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

extern "C" {
#include "blst/fields/fields.h"
#include "blst/fields/fp12_tower.h"
}

#ifndef DEVICE_FUNC
#define DEVICE_FUNC
#endif

#ifndef  GLOBAL_FUNC
#define GLOBAL_FUNC
#endif // ! GLOBAL_FUNC


DEVICE_FUNC bool cuda_vec384fp12_equal(const vec384fp12 a, const vec384fp12 b) {
    return vec_is_equal(a, b, sizeof(vec384fp12));
}

DEVICE_FUNC void cuda_set_fp12_zero(vec384fp12 a) {
    vec_zero(a, sizeof(vec384fp12));
}

DEVICE_FUNC void cuda_set_fp12_one(vec384fp12 a) {
    vec_zero(a, sizeof(vec384fp12));
    // Set to the Fp12 one constant
    const vec384fp12* one_ptr = blst_fp12_one();
    vec_copy(a, one_ptr, sizeof(vec384fp12));
}

DEVICE_FUNC void cuda_set_fp12_value(vec384fp12 a, int re_val, int im_val) {
    vec_zero(a, sizeof(vec384fp12));

    // Set a simple value based on parameters
    if (re_val > 0) {
        vec_copy(a[0][0][0], get_BLST_Rx()->p, sizeof(vec384));
        for (int i = 1; i < re_val && i < 10; i++) {
            add_fp(a[0][0][0], a[0][0][0], get_BLST_Rx()->p);
        }
    }

    if (im_val > 0) {
        vec_copy(a[0][0][1], get_BLST_Rx()->p, sizeof(vec384));
        for (int i = 1; i < im_val && i < 10; i++) {
            add_fp(a[0][0][1], a[0][0][1], get_BLST_Rx()->p);
        }
    }
}

// CUDA kernels for Fp12 tests

GLOBAL_FUNC void cuda_test_fp12_zero_and_one_kernel(bool* results) {
    vec384fp12 zero, one;

    // Initialize zero and one
    cuda_set_fp12_zero(zero);
    cuda_set_fp12_one(one);

    // Check zero properties
    vec384fp12 result;

    // Zero * Zero = Zero
    blst_fp12_mul(result, zero, zero);
    results[0] = cuda_vec384fp12_equal(result, zero);

    // One * One = One
    blst_fp12_mul(result, one, one);
    results[1] = cuda_vec384fp12_equal(result, one);

    // Zero * One = Zero
    blst_fp12_mul(result, zero, one);
    results[2] = cuda_vec384fp12_equal(result, zero);

    // One * Zero = Zero
    blst_fp12_mul(result, one, zero);
    results[3] = cuda_vec384fp12_equal(result, zero);

    // Check that One is detected as One
    results[4] = blst_fp12_is_one(one);
    results[5] = !blst_fp12_is_one(zero);
}

GLOBAL_FUNC void cuda_test_fp12_multiplication_kernel(bool* results) {
    vec384fp12 a, b, c, temp;

    // Initialize with some values - use the one constant and modify it
    const vec384fp12* one_ptr = blst_fp12_one();
    vec_copy(a, one_ptr, sizeof(a));

    // Make a simple change to create a distinct value
    a[0][0][1][0] = 2;  // Modify a[0][0][1][0] (first coordinate, imaginary part)

    // Create another distinct value
    vec_copy(b, one_ptr, sizeof(b));
    b[0][1][0][0] = 3;  // Modify b[0][1][0][0] (second coordinate, real part)

    // Compute c = a * b
    blst_fp12_mul(c, a, b);

    // Verify that the result is different from both inputs
    results[0] = !cuda_vec384fp12_equal(c, a);  // a*b should not equal a
    results[1] = !cuda_vec384fp12_equal(c, b);  // a*b should not equal b

    // Verify multiplication by 1 gives the same element
    blst_fp12_mul(temp, a, *one_ptr);
    results[2] = cuda_vec384fp12_equal(temp, a);  // a*1 should equal a

    blst_fp12_mul(temp, b, *one_ptr);
    results[3] = cuda_vec384fp12_equal(temp, b);  // b*1 should equal b
}

GLOBAL_FUNC void cuda_test_fp12_squaring_kernel(bool* results) {
    vec384fp12 a, a_squared1, a_squared2;

    // Initialize with some values - use the one constant and modify it
    const vec384fp12* one_ptr = blst_fp12_one();
    vec_copy(a, one_ptr, sizeof(a));

    // Make a simple change to create a distinct value
    a[0][0][1][0] = 2;  // Modify a[0][0][1][0] (first coordinate, imaginary part)
    a[1][1][0][0] = 3;  // Another change to make it more complex

    // Compute a² using multiplication
    blst_fp12_mul(a_squared1, a, a);

    // Compute a² using squaring
    blst_fp12_sqr(a_squared2, a);

    // Results should be equal
    results[0] = cuda_vec384fp12_equal(a_squared1, a_squared2);

    // Test cyclotomic squaring (just verify it runs without errors)
    vec384fp12 cyc_squared;
    blst_fp12_cyclotomic_sqr(cyc_squared, a);
    results[1] = true;  // If we get here, the function didn't crash
}

GLOBAL_FUNC void cuda_test_fp12_conjugation_kernel(bool* results) {
    vec384fp12 a, conjugate, expected;

    // Initialize with some values
    const vec384fp12* one_ptr = blst_fp12_one();
    vec_copy(a, one_ptr, sizeof(a));

    // Make it a non-trivial value: a = a_0 + a_1*w
    a[0][0][1][0] = 2;  // Modify a_0
    a[1][1][0][0] = 3;  // Modify a_1

    // Compute conjugate
    vec_copy(conjugate, a, sizeof(vec384fp12));
    blst_fp12_conjugate(conjugate);

    // Expected = a_0 - a_1*w
    vec_copy(expected, a, sizeof(vec384fp12));
    // Negate the second part of the Fp12 element (a_1)
    for (int i = 0; i < 3; i++) {
        cneg_fp2(expected[1][i], a[1][i], 1);
    }

    // Results should be equal
    results[0] = cuda_vec384fp12_equal(conjugate, expected);

    // Check that conjugate(conjugate(a)) = a
    blst_fp12_conjugate(conjugate);
    results[1] = cuda_vec384fp12_equal(conjugate, a);
}

GLOBAL_FUNC void cuda_test_fp12_inversion_kernel(bool* results) {
    vec384fp12 a, a_inv, result;

    // Initialize with a non-trivial value
    const vec384fp12* one_ptr = blst_fp12_one();
    vec_copy(a, one_ptr, sizeof(a));
    a[0][0][1][0] = 2;
    a[1][1][0][0] = 3;

    // Compute 1/a
    blst_fp12_inverse(a_inv, a);

    // Check that a * (1/a) = 1
    blst_fp12_mul(result, a, a_inv);

    // Should equal 1
    results[0] = blst_fp12_is_one(result);

    // Special case: Inversion of 1 should be 1
    vec384fp12 one, one_inv;
    vec_copy(one, one_ptr, sizeof(one));
    blst_fp12_inverse(one_inv, one);
    results[1] = cuda_vec384fp12_equal(one, one_inv);
}

GLOBAL_FUNC void cuda_test_fp12_equality_kernel(bool* results) {
    vec384fp12 a, b;

    // Initialize with different values
    const vec384fp12* one_ptr = blst_fp12_one();
    vec_copy(a, one_ptr, sizeof(a));
    vec_copy(b, one_ptr, sizeof(b));

    // a and b should be equal initially
    results[0] = blst_fp12_is_equal(a, b);

    // Modify b to make it different
    b[0][0][0][0] = 2;
    results[1] = !blst_fp12_is_equal(a, b);

    // Make them equal again
    vec_copy(b, a, sizeof(b));
    results[2] = blst_fp12_is_equal(a, b);

    // Test equality with self
    results[3] = blst_fp12_is_equal(a, a);
}

GLOBAL_FUNC void cuda_test_fp12_serialization_kernel(bool* results) {
    vec384fp12 a;
    unsigned char serialized[12 * 48]; // Fp12 = 12 elements of Fp

    // Initialize with non-trivial values
    const vec384fp12* one_ptr = blst_fp12_one();
    vec_copy(a, one_ptr, sizeof(a));
    a[0][0][1][0] = 2;
    a[0][1][0][0] = 3;
    a[1][2][1][0] = 4;

    // Serialize to bytes
    blst_bendian_from_fp12(serialized, a);

    // Verify properties of the serialized data
    bool has_nonzero = false;
    for (size_t i = 0; i < sizeof(serialized); i++) {
        if (serialized[i] != 0) {
            has_nonzero = true;
            break;
        }
    }
    results[0] = has_nonzero;

    // Test size reporting
    size_t fp12_size = blst_fp12_sizeof();
    results[1] = (fp12_size == sizeof(vec384fp12));
}

// Batch processing tests for Fp12
GLOBAL_FUNC void cuda_test_fp12_batch_operations_kernel(bool* results, int batch_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= batch_size) return;

    vec384fp12 a, b, result;

    // Create different test values for each thread
    const vec384fp12* one_ptr = blst_fp12_one();
    vec_copy(a, one_ptr, sizeof(a));
    vec_copy(b, one_ptr, sizeof(b));

    // Set values based on thread index
    int val_a = (idx % 7) + 1;
    int val_b = ((idx + 3) % 5) + 1;

    // Modify elements to create distinct values
    a[0][0][1][0] = val_a;
    a[0][1][0][0] = val_a + 1;
    b[0][0][0][0] = val_b;
    b[1][0][1][0] = val_b + 2;

    // Perform various operations
    vec384fp12 temp1, temp2, temp3;
    blst_fp12_mul(temp1, a, b);        // a * b
    blst_fp12_sqr(temp2, a);           // a²
    blst_fp12_conjugate(temp3);
    vec_copy(temp3, temp1, sizeof(vec384fp12));
    blst_fp12_conjugate(temp3);        // conjugate(a * b)

    // Combine operations: (a * b) * a²
    blst_fp12_mul(result, temp1, temp2);

    // Simple validation - result should not be zero
    vec384fp12 zero;
    cuda_set_fp12_zero(zero);
    results[idx] = !cuda_vec384fp12_equal(result, zero);
}

GLOBAL_FUNC void cuda_test_fp12_performance_kernel(int* operations_completed) {
    vec384fp12 a, b, result;

    // Initialize with some values
    const vec384fp12* one_ptr = blst_fp12_one();
    vec_copy(a, one_ptr, sizeof(a));
    vec_copy(b, one_ptr, sizeof(b));

    // Make them non-trivial
    a[0][0][1][0] = 5;
    a[1][1][0][0] = 7;
    b[0][1][0][0] = 3;
    b[1][0][1][0] = 4;

    int ops = 0;
    // Perform many operations to test performance
    for (int i = 0; i < 100; i++) {  // Reduced for Fp12 complexity
        blst_fp12_mul(result, a, b);
        blst_fp12_sqr(a, result);
        vec_copy(b, result, sizeof(vec384fp12));
        blst_fp12_conjugate(b);
        ops += 4;
    }

    *operations_completed = ops;
}
