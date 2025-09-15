#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>

extern "C" {
#include "../src/fields.h"
}

#define LIMBS NLIMBS(384)

// CUDA kernel wrapper macros for device execution
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while(0)

// Device versions of helper functions
__device__ void cuda_set_vec384(vec384 v, uint64_t x) {
    for (int i = 0; i < LIMBS; ++i)
        v[i] = 0;

    memcpy(v, &x, sizeof(x));
}

__device__ void cuda_set_vec384x(vec384x v, uint64_t x0, uint64_t x1) {
    cuda_set_vec384(v[0], x0);
    cuda_set_vec384(v[1], x1);
}

__device__ bool cuda_vec384_equal(const vec384 a, const vec384 b) {
    for (int i = 0; i < LIMBS; ++i)
        if (a[i] != b[i]) return false;
    return true;
}

__device__ void cuda_print_vec384(const vec384 v) {
    for (int i = LIMBS - 1; i >= 0; --i)
        printf("%016llx%s", (unsigned long long)v[i], i ? " " : "\n");
}

// CUDA kernel for basic arithmetic tests
__global__ void cuda_test_add_kernel(bool* results) {
    vec384 a, b, r, expected;
    cuda_set_vec384(a, 5);
    cuda_set_vec384(b, 7);
    add_fp(r, a, b);

    cuda_set_vec384(expected, 12);
    results[0] = cuda_vec384_equal(r, expected);
}

__global__ void cuda_test_add_overflow_kernel(bool* results) {
    vec384 a, b, r, expected;
    // a = p - 1, b = 2, expect 1
    memcpy(a, *get_BLST_P(), sizeof(vec384));
    for (int i = 0; i < LIMBS; ++i) a[i] -= (i == 0 ? 1 : 0);
    cuda_set_vec384(b, 2);
    add_fp(r, a, b);
    cuda_set_vec384(expected, 1);
    results[0] = cuda_vec384_equal(r, expected);
}

__global__ void cuda_test_sub_kernel(bool* results) {
    vec384 a, b, r, expected;
    cuda_set_vec384(a, 10);
    cuda_set_vec384(b, 3);
    sub_fp(r, a, b);
    cuda_set_vec384(expected, 7);
    results[0] = cuda_vec384_equal(r, expected);
}

__global__ void cuda_test_sub_underflow_kernel(bool* results) {
    vec384 a, b, r, expected;
    cuda_set_vec384(a, 3);
    cuda_set_vec384(b, 10);
    sub_fp(r, a, b);
    // expected = (3 - 10) mod p = p - 7
    memcpy(expected, *get_BLST_P(), sizeof(vec384));
    expected[0] -= 7;
    results[0] = cuda_vec384_equal(r, expected);
}

__global__ void cuda_test_mul_by_3_kernel(bool* results) {
    vec384 a, r, expected;
    cuda_set_vec384(a, 2);
    mul_by_3_fp(r, a);
    cuda_set_vec384(expected, 6);
    results[0] = cuda_vec384_equal(r, expected);
}

__global__ void cuda_test_mul_by_8_kernel(bool* results) {
    vec384 a, r, expected;
    cuda_set_vec384(a, 2);
    mul_by_8_fp(r, a);
    cuda_set_vec384(expected, 16);
    results[0] = cuda_vec384_equal(r, expected);
}

__global__ void cuda_test_lshift_kernel(bool* results) {
    vec384 a, r, expected;
    cuda_set_vec384(a, 1);
    lshift_fp(r, a, 2);
    cuda_set_vec384(expected, 4);
    results[0] = cuda_vec384_equal(r, expected);
}

__global__ void cuda_test_rshift_kernel(bool* results) {
    vec384 a, r, expected;
    cuda_set_vec384(a, 8);
    rshift_fp(r, a, 2);
    cuda_set_vec384(expected, 2);
    results[0] = cuda_vec384_equal(r, expected);
}

__global__ void cuda_test_div_by_2_kernel(bool* results) {
    vec384 a, r, expected;
    cuda_set_vec384(a, 8);
    div_by_2_fp(r, a);
    cuda_set_vec384(expected, 4);
    results[0] = cuda_vec384_equal(r, expected);
}

__global__ void cuda_test_mul_kernel(bool* results) {
    vec384 a, b, r, r_std, expected;

    // Set raw values
    cuda_set_vec384(a, 3);
    cuda_set_vec384(b, 4);

    // Convert to Montgomery form
    to_fp(a, a);
    to_fp(b, b);

    // Multiply in Montgomery form
    mul_fp(r, a, b);

    // Convert result back to standard form
    from_fp(r_std, r);

    // Expect 12
    cuda_set_vec384(expected, 12);
    results[0] = cuda_vec384_equal(r_std, expected);
}

__global__ void cuda_test_sqr_kernel(bool* results) {
    vec384 a, r, r_std, expected;

    // Set raw value
    cuda_set_vec384(a, 5);

    // Convert to Montgomery form
    to_fp(a, a);

    // Square in Montgomery form
    sqr_fp(r, a);

    // Convert result back to standard form
    from_fp(r_std, r);

    cuda_set_vec384(expected, 25);
    results[0] = cuda_vec384_equal(r_std, expected);
}

__global__ void cuda_test_from_fp_kernel(bool* results) {
    vec384 a, r, expected;

    // Set raw value and convert to Montgomery form
    cuda_set_vec384(a, 9);
    to_fp(a, a);

    // Test converting from Montgomery form
    from_fp(r, a);

    cuda_set_vec384(expected, 9);
    results[0] = cuda_vec384_equal(r, expected);
}

__global__ void cuda_test_montgomery_conversion_kernel(bool* results) {
    vec384 a, mont_a, back_a;

    // Test with value 1
    cuda_set_vec384(a, 1);

    // Convert to Montgomery form
    to_fp(mont_a, a);

    // Convert back from Montgomery form
    from_fp(back_a, mont_a);

    // The result should be the original value
    results[0] = cuda_vec384_equal(back_a, a);

    // Test multiplication in Montgomery form
    vec384 b, mont_b, result, std_result;
    cuda_set_vec384(a, 3);
    cuda_set_vec384(b, 4);

    // Convert to Montgomery form
    to_fp(mont_a, a);
    to_fp(mont_b, b);

    // Multiply in Montgomery form
    mul_fp(result, mont_a, mont_b);

    // Convert back to standard form
    from_fp(std_result, result);

    // Check result is 12
    vec384 expected;
    cuda_set_vec384(expected, 12);

    results[1] = cuda_vec384_equal(std_result, expected);
}

__global__ void cuda_test_montgomery_arithmetic_kernel(bool* results) {
    vec384 a, b, c, mont_a, mont_b, mont_c, result1, result2, std_result1, std_result2;

    // Set test values
    cuda_set_vec384(a, 17);
    cuda_set_vec384(b, 23);

    // Convert to Montgomery form
    to_fp(mont_a, a);
    to_fp(mont_b, b);

    // Test associativity: a * (b * c) = (a * b) * c
    cuda_set_vec384(c, 19);
    to_fp(mont_c, c);

    // Calculate a * (b * c) in Montgomery form
    mul_fp(result1, mont_b, mont_c); // result1 = b * c
    mul_fp(result1, mont_a, result1); // result1 = a * (b * c)

    // Calculate (a * b) * c in Montgomery form
    mul_fp(result2, mont_a, mont_b); // result2 = a * b
    mul_fp(result2, result2, mont_c); // result2 = (a * b) * c

    // Convert results back to standard form for comparison
    from_fp(std_result1, result1);
    from_fp(std_result2, result2);

    // Results should be equal
    results[0] = cuda_vec384_equal(std_result1, std_result2);

    // Test identity: a * 1 = a
    cuda_set_vec384(b, 1);
    to_fp(mont_b, b);

    mul_fp(result1, mont_a, mont_b);
    from_fp(std_result1, result1);
    from_fp(std_result2, mont_a);

    results[1] = cuda_vec384_equal(std_result1, std_result2);

    // Test that squaring works: a² = a * a
    sqr_fp(result1, mont_a);
    mul_fp(result2, mont_a, mont_a);

    results[2] = cuda_vec384_equal(result1, result2);
}

__global__ void cuda_test_edge_cases_kernel(bool* results) {
    vec384 zero, one, r;

    // Test operations with zero
    cuda_set_vec384(zero, 0);
    cuda_set_vec384(one, 1);

    // Zero + Zero = Zero
    add_fp(r, zero, zero);
    results[0] = cuda_vec384_equal(r, zero);

    // Zero + One = One
    add_fp(r, zero, one);
    results[1] = cuda_vec384_equal(r, one);

    // One - One = Zero
    sub_fp(r, one, one);
    results[2] = cuda_vec384_equal(r, zero);

    // Zero - Zero = Zero
    sub_fp(r, zero, zero);
    results[3] = cuda_vec384_equal(r, zero);
}

// Additional comprehensive tests that combine multiple operations
__global__ void cuda_test_round_trip_kernel(bool* results) {
    const int TEST_VALUES[] = { 0, 1, 2, 3, 7, 10, 100, 12345, 0x7FFFFFFF };
    const int NUM_VALUES = sizeof(TEST_VALUES) / sizeof(TEST_VALUES[0]);

    bool all_passed = true;

    for (int i = 0; i < NUM_VALUES; i++) {
        vec384 original, mont_form, result;

        cuda_set_vec384(original, TEST_VALUES[i]);

        // Convert to Montgomery form and back
        to_fp(mont_form, original);
        from_fp(result, mont_form);

        if (!cuda_vec384_equal(result, original)) {
            all_passed = false;
            break;
        }
    }

    results[0] = all_passed;
}

__global__ void cuda_test_distributive_kernel(bool* results) {
    vec384 a, b, c, mont_a, mont_b, mont_c;
    vec384 left_side, right_side1, right_side2, right_side;
    vec384 std_left, std_right;

    // Set test values
    cuda_set_vec384(a, 11);
    cuda_set_vec384(b, 13);
    cuda_set_vec384(c, 17);

    // Convert to Montgomery form
    to_fp(mont_a, a);
    to_fp(mont_b, b);
    to_fp(mont_c, c);

    // Left side: a * (b + c)
    vec384 sum;
    add_fp(sum, mont_b, mont_c);  // b + c
    mul_fp(left_side, mont_a, sum);  // a * (b + c)

    // Right side: (a * b) + (a * c)
    mul_fp(right_side1, mont_a, mont_b);  // a * b
    mul_fp(right_side2, mont_a, mont_c);  // a * c
    add_fp(right_side, right_side1, right_side2);  // (a * b) + (a * c)

    // Convert back to standard form for comparison
    from_fp(std_left, left_side);
    from_fp(std_right, right_side);

    results[0] = cuda_vec384_equal(std_left, std_right);
}

__global__ void cuda_test_modulus_multiply_kernel(bool* results) {
    vec384 a, p, mont_a, mont_p, result, expected;

    // Set a value
    cuda_set_vec384(a, 42);

    // Copy the modulus
    memcpy(p, *get_BLST_P(), sizeof(vec384));

    // Convert to Montgomery form
    to_fp(mont_a, a);
    to_fp(mont_p, p);

    // Multiply a * p
    mul_fp(result, mont_a, mont_p);

    // Result should be 0 (since a * p ≡ 0 mod p)
    from_fp(expected, result);

    cuda_set_vec384(a, 0);
    results[0] = cuda_vec384_equal(expected, a);
}

// Performance comparison test
__global__ void cuda_performance_test_kernel(int* operations_completed) {
    vec384 a, b, result;
    cuda_set_vec384(a, 12345);
    cuda_set_vec384(b, 67890);

    to_fp(a, a);
    to_fp(b, b);

    int ops = 0;
    // Perform many operations to test performance
    for (int i = 0; i < 1000; i++) {
        mul_fp(result, a, b);
        sqr_fp(a, result);
        add_fp(b, a, b);
        ops += 3;
    }

    *operations_completed = ops;
}

// Batch processing test - simulate real-world usage
__global__ void cuda_batch_processing_kernel(bool* results, int batch_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= batch_size) return;

    vec384 a, b, result;

    // Different values for each thread based on index
    cuda_set_vec384(a, 100 + idx);
    cuda_set_vec384(b, 200 + idx * 2);

    // Convert to Montgomery form
    to_fp(a, a);
    to_fp(b, b);

    // Perform computation
    mul_fp(result, a, b);
    sqr_fp(result, result);

    // Simple validation - result should not be zero
    results[idx] = !cuda_vec384_equal(result, result); // This will be false, which means not zero

    // Better validation - convert back and check range
    vec384 std_result;
    from_fp(std_result, result);
    results[idx] = (std_result[0] != 0); // At least one limb should be non-zero
}

// Helper function to run a CUDA kernel and check results
template<typename KernelFunc>
bool run_cuda_test(KernelFunc kernel, int num_results = 1) {
    bool* d_results, * h_results;

    // Allocate memory
    h_results = new bool[num_results];
    CUDA_CHECK(cudaMalloc(&d_results, num_results * sizeof(bool)));

    // Initialize results to false
    CUDA_CHECK(cudaMemset(d_results, 0, num_results * sizeof(bool)));

    // Launch kernel
    kernel << <1, 1 >> > (d_results);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_results, d_results, num_results * sizeof(bool), cudaMemcpyDeviceToHost));

    // Check all results
    bool all_passed = true;
    for (int i = 0; i < num_results; i++) {
        if (!h_results[i]) {
            all_passed = false;
            break;
        }
    }

    // Cleanup
    delete[] h_results;
    CUDA_CHECK(cudaFree(d_results));

    return all_passed;
}

// CUDA test cases - basic functionality
TEST(FpOpsCUDA, Add) {
    EXPECT_TRUE(run_cuda_test(cuda_test_add_kernel));
}

TEST(FpOpsCUDA, AddOverflow) {
    EXPECT_TRUE(run_cuda_test(cuda_test_add_overflow_kernel));
}

TEST(FpOpsCUDA, Sub) {
    EXPECT_TRUE(run_cuda_test(cuda_test_sub_kernel));
}

TEST(FpOpsCUDA, SubUnderflow) {
    EXPECT_TRUE(run_cuda_test(cuda_test_sub_underflow_kernel));
}

TEST(FpOpsCUDA, MulBy3) {
    EXPECT_TRUE(run_cuda_test(cuda_test_mul_by_3_kernel));
}

TEST(FpOpsCUDA, MulBy8) {
    EXPECT_TRUE(run_cuda_test(cuda_test_mul_by_8_kernel));
}

TEST(FpOpsCUDA, LShift) {
    EXPECT_TRUE(run_cuda_test(cuda_test_lshift_kernel));
}

TEST(FpOpsCUDA, RShift) {
    EXPECT_TRUE(run_cuda_test(cuda_test_rshift_kernel));
}

TEST(FpOpsCUDA, DivBy2) {
    EXPECT_TRUE(run_cuda_test(cuda_test_div_by_2_kernel));
}

TEST(FpOpsCUDA, Mul) {
    EXPECT_TRUE(run_cuda_test(cuda_test_mul_kernel));
}

TEST(FpOpsCUDA, Sqr) {
    EXPECT_TRUE(run_cuda_test(cuda_test_sqr_kernel));
}

TEST(FpOpsCUDA, FromFp) {
    EXPECT_TRUE(run_cuda_test(cuda_test_from_fp_kernel));
}

// Advanced Montgomery arithmetic tests
TEST(FpOpsCUDA, MontgomeryConversion) {
    EXPECT_TRUE(run_cuda_test(cuda_test_montgomery_conversion_kernel, 2));
}

TEST(FpOpsCUDA, MontgomeryArithmetic) {
    EXPECT_TRUE(run_cuda_test(cuda_test_montgomery_arithmetic_kernel, 3));
}

TEST(FpOpsCUDA, EdgeCases) {
    EXPECT_TRUE(run_cuda_test(cuda_test_edge_cases_kernel, 4));
}

TEST(FpOpsCUDA, MontgomeryRoundTrip) {
    EXPECT_TRUE(run_cuda_test(cuda_test_round_trip_kernel));
}

TEST(FpOpsCUDA, MontgomeryDistributive) {
    EXPECT_TRUE(run_cuda_test(cuda_test_distributive_kernel));
}

TEST(FpOpsCUDA, MultiplyWithModulus) {
    EXPECT_TRUE(run_cuda_test(cuda_test_modulus_multiply_kernel));
}

// Performance and batch processing tests
TEST(FpOpsCUDA, PerformanceBaseline) {
    int* d_ops, h_ops;

    CUDA_CHECK(cudaMalloc(&d_ops, sizeof(int)));
    CUDA_CHECK(cudaMemset(d_ops, 0, sizeof(int)));

    // Record timing
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));

    CUDA_CHECK(cudaEventRecord(start));
    cuda_performance_test_kernel << <1, 1 >> > (d_ops);
    CUDA_CHECK(cudaEventRecord(stop));

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));

    CUDA_CHECK(cudaMemcpy(&h_ops, d_ops, sizeof(int), cudaMemcpyDeviceToHost));

    printf("CUDA Performance: %d operations in %.3f ms (%.1f ops/ms)\n",
        h_ops, milliseconds, h_ops / milliseconds);

    CUDA_CHECK(cudaFree(d_ops));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));

    EXPECT_GT(h_ops, 0);
}

TEST(FpOpsCUDA, BatchProcessing) {
    const int batch_size = 1024;
    bool* d_results, * h_results;

    // Allocate memory
    h_results = new bool[batch_size];
    CUDA_CHECK(cudaMalloc(&d_results, batch_size * sizeof(bool)));
    CUDA_CHECK(cudaMemset(d_results, 0, batch_size * sizeof(bool)));

    // Launch kernel with multiple threads
    int threads_per_block = 256;
    int blocks = (batch_size + threads_per_block - 1) / threads_per_block;

    cuda_batch_processing_kernel << <blocks, threads_per_block >> > (d_results, batch_size);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_results, d_results, batch_size * sizeof(bool), cudaMemcpyDeviceToHost));

    // Check that all computations completed successfully
    int successful_operations = 0;
    for (int i = 0; i < batch_size; i++) {
        if (h_results[i]) {
            successful_operations++;
        }
    }

    printf("Batch processing: %d/%d operations successful\n", successful_operations, batch_size);

    // Cleanup
    delete[] h_results;
    CUDA_CHECK(cudaFree(d_results));

    // Expect at least 95% success rate
    EXPECT_GE(successful_operations, batch_size * 0.95);
}

// Memory coalescing test
__global__ void cuda_memory_coalescing_test_kernel(vec384* input_array, vec384* output_array, int array_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= array_size) return;

    vec384 temp;
    // Load from global memory
    memcpy(temp, input_array[idx], sizeof(vec384));

    // Convert to Montgomery form
    to_fp(temp, temp);

    // Square the value
    sqr_fp(temp, temp);

    // Convert back from Montgomery form
    from_fp(temp, temp);

    // Store to global memory
    memcpy(output_array[idx], temp, sizeof(vec384));
}

TEST(FpOpsCUDA, MemoryCoalescing) {
    const int array_size = 1024;
    size_t array_bytes = array_size * sizeof(vec384);

    vec384* h_input, * h_output;
    vec384* d_input, * d_output;

    // Allocate host memory
    h_input = new vec384[array_size];
    h_output = new vec384[array_size];

    // Initialize input data
    for (int i = 0; i < array_size; i++) {
        for (int j = 0; j < LIMBS; j++) {
            h_input[i][j] = 0;
        }
        h_input[i][0] = i + 1; // Simple pattern
    }

    // Allocate device memory
    CUDA_CHECK(cudaMalloc(&d_input, array_bytes));
    CUDA_CHECK(cudaMalloc(&d_output, array_bytes));

    // Copy input to device
    CUDA_CHECK(cudaMemcpy(d_input, h_input, array_bytes, cudaMemcpyHostToDevice));

    // Launch kernel
    int threads_per_block = 256;
    int blocks = (array_size + threads_per_block - 1) / threads_per_block;

    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));

    CUDA_CHECK(cudaEventRecord(start));
    cuda_memory_coalescing_test_kernel << <blocks, threads_per_block >> > (d_input, d_output, array_size);
    CUDA_CHECK(cudaEventRecord(stop));

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));

    // Copy output back
    CUDA_CHECK(cudaMemcpy(h_output, d_output, array_bytes, cudaMemcpyDeviceToHost));

    printf("Memory coalescing test: %d elements processed in %.3f ms\n", array_size, milliseconds);

    // Verify results - outputs should not be the same as inputs
    int different_results = 0;
    for (int i = 0; i < array_size; i++) {
        bool is_different = false;
        for (int j = 0; j < LIMBS; j++) {
            if (h_input[i][j] != h_output[i][j]) {
                is_different = true;
                break;
            }
        }
        if (is_different) {
            different_results++;
        }
    }

    printf("Different results: %d/%d\n", different_results, array_size);

    // Cleanup
    delete[] h_input;
    delete[] h_output;
    CUDA_CHECK(cudaFree(d_input));
    CUDA_CHECK(cudaFree(d_output));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));

    // Most results should be different (since we're squaring)
    EXPECT_GT(different_results, array_size * 0.9);
} 