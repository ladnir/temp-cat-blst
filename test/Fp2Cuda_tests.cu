#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

extern "C" {
#include "../src/fields.h"
}

// CUDA kernel wrapper macros for device execution
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while(0)

// Device helper functions for Fp2 operations
DEVICE_FUNC void cuda_print_fp2(const vec384x a, const char* name) {
    printf("%s:\n", name);
    printf("  Re: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)a[0][i]);
    }
    printf("\n  Im: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)a[1][i]);
    }
    printf("\n");
}

DEVICE_FUNC bool cuda_vec384x_equal(const vec384x a, const vec384x b) {
    return vec_is_equal(a, b, sizeof(vec384x));
}

DEVICE_FUNC void cuda_set_fp2_zero(vec384x a) {
    vec_zero(a, sizeof(vec384x));
}

DEVICE_FUNC void cuda_set_fp2_one(vec384x a) {
    vec_zero(a, sizeof(vec384x));
    vec_copy(a[0], get_BLST_Rx()->p, sizeof(vec384));  // Re(1) = 1
}

DEVICE_FUNC void cuda_set_fp2_value(vec384x a, int re_val, int im_val) {
    vec_zero(a, sizeof(vec384x));

    // Set real part
    if (re_val > 0) {
        vec_copy(a[0], get_BLST_Rx()->p, sizeof(vec384));
        for (int i = 1; i < re_val; i++) {
            add_fp(a[0], a[0], get_BLST_Rx()->p);
        }
    }

    // Set imaginary part
    if (im_val > 0) {
        vec_copy(a[1], get_BLST_Rx()->p, sizeof(vec384));
        for (int i = 1; i < im_val; i++) {
            add_fp(a[1], a[1], get_BLST_Rx()->p);
        }
    }
}

// CUDA kernels for Fp2 tests

GLOBAL_FUNC void cuda_test_fp2_zero_and_one_kernel(bool* results) {
    vec384x zero, one;

    // Initialize zero and one
    cuda_set_fp2_zero(zero);
    cuda_set_fp2_one(one);

    // Check zero properties
    vec384x result;
    add_fp2(result, zero, zero);
    results[0] = cuda_vec384x_equal(result, zero); // 0 + 0 = 0

    mul_fp2(result, zero, one);
    results[1] = cuda_vec384x_equal(result, zero); // 0 * 1 = 0

    // Check one properties
    mul_fp2(result, one, one);
    results[2] = cuda_vec384x_equal(result, one); // 1 * 1 = 1

    vec384x temp;
    vec_copy(temp, one, sizeof(vec384x));
    add_fp2(result, temp, temp);

    // 1 + 1 should be 2 (in Montgomery representation)
    vec384x expected_two;
    cuda_set_fp2_zero(expected_two);
    lshift_fp(expected_two[0], one[0], 1);  // Re(2) = 2

    results[3] = cuda_vec384x_equal(result, expected_two); // 1 + 1 = 2
}

GLOBAL_FUNC void cuda_test_fp2_addition_kernel(bool* results) {
    vec384x a, b, c, expected;

    // Initialize with some values
    // a = 2 + 3i (in Montgomery representation)
    cuda_set_fp2_zero(a);
    lshift_fp(a[0], get_BLST_Rx()->p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], get_BLST_Rx()->p);      // Im(a) = 3

    // b = 4 + 5i (in Montgomery representation)
    cuda_set_fp2_zero(b);
    lshift_fp(b[0], get_BLST_Rx()->p, 2);     // Re(b) = 4
    add_fp(b[1], get_BLST_Rx()->p, b[0]);     // Im(b) = 1 + 4 = 5

    // c = a + b should be 6 + 8i
    add_fp2(c, a, b);

    cuda_set_fp2_zero(expected);
    // Calculate Re(expected) = 6
    mul_by_3_fp(expected[0], get_BLST_Rx()->p);       // Re(expected) = 3
    lshift_fp(expected[0], expected[0], 1);         // Re(expected) = 6
    // Calculate Im(expected) = 8
    lshift_fp(expected[1], get_BLST_Rx()->p, 3);      // Im(expected) = 8

    results[0] = cuda_vec384x_equal(c, expected); // (2+3i) + (4+5i) = (6+8i)

    // Test commutativity
    add_fp2(c, b, a);
    results[1] = cuda_vec384x_equal(c, expected); // (4+5i) + (2+3i) = (6+8i)
}

GLOBAL_FUNC void cuda_test_fp2_subtraction_kernel(bool* results) {
    vec384x a, b, c, expected;

    // Initialize with some values
    // a = 7 + 4i (in Montgomery representation)
    cuda_set_fp2_zero(a);
    mul_by_3_fp(a[0], get_BLST_Rx()->p);      // Re(a) = 3
    add_fp(a[0], a[0], get_BLST_Rx()->p);     // Re(a) = 3+1 = 4
    add_fp(a[0], a[0], get_BLST_Rx()->p);     // Re(a) = 4+1 = 5
    add_fp(a[0], a[0], get_BLST_Rx()->p);     // Re(a) = 5+1 = 6
    add_fp(a[0], a[0], get_BLST_Rx()->p);     // Re(a) = 6+1 = 7
    lshift_fp(a[1], get_BLST_Rx()->p, 2);     // Im(a) = 4

    // b = 2 + 3i (in Montgomery representation)
    cuda_set_fp2_zero(b);
    lshift_fp(b[0], get_BLST_Rx()->p, 1);     // Re(b) = 2
    mul_by_3_fp(b[1], get_BLST_Rx()->p);      // Im(b) = 3

    // c = a - b should be 5 + 1i
    sub_fp2(c, a, b);

    cuda_set_fp2_zero(expected);
    add_fp(expected[0], get_BLST_Rx()->p, get_BLST_Rx()->p);  // Re(expected) = 2
    add_fp(expected[0], expected[0], get_BLST_Rx()->p);     // Re(expected) = 2+1 = 3
    add_fp(expected[0], expected[0], get_BLST_Rx()->p);     // Re(expected) = 3+1 = 4
    add_fp(expected[0], expected[0], get_BLST_Rx()->p);     // Re(expected) = 4+1 = 5
    vec_copy(expected[1], get_BLST_Rx()->p, sizeof(vec384));  // Im(expected) = 1

    results[0] = cuda_vec384x_equal(c, expected); // (7+4i) - (2+3i) = (5+1i)

    // c = b - a should be -5 - 1i
    sub_fp2(c, b, a);

    // Negate expected to get -5 - 1i
    cneg_fp2(expected, expected, 1);

    results[1] = cuda_vec384x_equal(c, expected); // (2+3i) - (7+4i) = (-5-1i)
}

GLOBAL_FUNC void cuda_test_fp2_multiplication_kernel(bool* results) {
    vec384x a, b, c, expected;

    // a = 2 + 3i
    cuda_set_fp2_zero(a);
    lshift_fp(a[0], get_BLST_Rx()->p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], get_BLST_Rx()->p);      // Im(a) = 3

    // b = 4 + 5i
    cuda_set_fp2_zero(b);
    lshift_fp(b[0], get_BLST_Rx()->p, 2);     // Re(b) = 4
    add_fp(b[1], get_BLST_Rx()->p, b[0]);     // Im(b) = 1 + 4 = 5

    // c = a * b should be (2*4 - 3*5) + (2*5 + 3*4)i = (8 - 15) + (10 + 12)i = -7 + 22i
    mul_fp2(c, a, b);

    cuda_set_fp2_zero(expected);
    // Re(expected) = -7
    mul_by_3_fp(expected[0], get_BLST_Rx()->p);      // 3
    lshift_fp(expected[0], expected[0], 1);        // 6
    add_fp(expected[0], expected[0], get_BLST_Rx()->p); // 7
    cneg_fp(expected[0], expected[0], 1);           // -7

    // Im(expected) = 22
    vec384 temp;
    lshift_fp(temp, get_BLST_Rx()->p, 1);             // 2
    mul_by_3_fp(expected[1], get_BLST_Rx()->p);       // 3
    lshift_fp(expected[1], expected[1], 3);         // 24
    sub_fp(expected[1], expected[1], temp);         // 22

    results[0] = cuda_vec384x_equal(c, expected); // (2+3i) * (4+5i) = (-7+22i)

    // Test commutativity
    mul_fp2(c, b, a);
    results[1] = cuda_vec384x_equal(c, expected); // (4+5i) * (2+3i) = (-7+22i)
}

GLOBAL_FUNC void cuda_test_fp2_square_and_multiply_kernel(bool* results) {
    vec384x a, a_squared1, a_squared2;

    // a = 2 + 3i
    cuda_set_fp2_zero(a);
    lshift_fp(a[0], get_BLST_Rx()->p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], get_BLST_Rx()->p);      // Im(a) = 3

    // Compute a² using multiplication
    mul_fp2(a_squared1, a, a);

    // Compute a² using squaring
    sqr_fp2(a_squared2, a);

    // Results should be equal: (2+3i)² = 4-9 + 12i = -5 + 12i
    results[0] = cuda_vec384x_equal(a_squared1, a_squared2);
}


//HOST_DEVICE_FUNC BLST_FORCEINLINE void reciprocal_fp2_my(vec384x out, const vec384x inp)
//{
//    vec384 t0, t1;
//
//    /*
//     * |out| = 1/(a + b*i) = a/(a^2+b^2) - b/(a^2+b^2)*i
//     */
//    sqr_fp(t0, inp[0]);
//    sqr_fp(t1, inp[1]);
//    add_fp(t0, t0, t1);
//    reciprocal_fp(t1, t0);
//    mul_fp(out[0], inp[0], t1);
//    mul_fp(out[1], inp[1], t1);
//    neg_fp(out[1], out[1]);
//}


GLOBAL_FUNC void cuda_test_fp2_inversion_kernel(bool* results) {
    vec384x a, a_inv, result;

    // a = 2 + 3i
    cuda_set_fp2_zero(a);
    lshift_fp(a[0], get_BLST_Rx()->p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], get_BLST_Rx()->p);      // Im(a) = 3

    // Compute 1/a
    reciprocal_fp2(a_inv, a);

    // Check that a * (1/a) = 1
    mul_fp2(result, a, a_inv);

    // Should equal 1 (real part 1, imaginary part 0)
    vec384x one;
    cuda_set_fp2_one(one);

    results[0] = cuda_vec384x_equal(result, one);
}

GLOBAL_FUNC void cuda_test_fp2_frobenius_kernel(bool* results) {
    vec384x a, frob_a, expected;

    // a = 2 + 3i
    cuda_set_fp2_zero(a);
    lshift_fp(a[0], get_BLST_Rx()->p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], get_BLST_Rx()->p);      // Im(a) = 3

    // Frobenius map in Fp2 is complex conjugation: (a+bi) -> (a-bi)
    // expected = 2 - 3i
    vec_copy(expected, a, sizeof(vec384x));
    cneg_fp(expected[1], expected[1], 1);  // Negate imaginary part

    // Apply Frobenius
    conjugate_fp2(frob_a, a);

    results[0] = cuda_vec384x_equal(frob_a, expected);
}

GLOBAL_FUNC void cuda_test_fp2_square_root_kernel(bool* results) {
    vec384x a, root, root_squared;

    // Let's try with a perfect square: a = (2+3i)²
    vec384x temp;
    cuda_set_fp2_zero(temp);
    lshift_fp(temp[0], get_BLST_Rx()->p, 1);    // Re(temp) = 2
    mul_by_3_fp(temp[1], get_BLST_Rx()->p);     // Im(temp) = 3

    sqr_fp2(a, temp);  // a = (2+3i)²

    // Find square root
    bool sqrt_success = sqrt_fp2(root, a);
    results[0] = sqrt_success;

    if (sqrt_success) {
        // Verify by squaring the root
        sqr_fp2(root_squared, root);

        // Should get back our original value
        results[1] = cuda_vec384x_equal(root_squared, a);

        // Verify that either root or -root equals our original value
        bool matches_orig = cuda_vec384x_equal(root, temp);
        cneg_fp2(root, root, 1);
        bool matches_neg = cuda_vec384x_equal(root, temp);

        results[2] = matches_orig || matches_neg;
    }
    else {
        results[1] = false;
        results[2] = false;
    }
}

GLOBAL_FUNC void cuda_test_fp2_serialization_kernel(bool* results) {
    vec384x a, b;
    unsigned char serialized[96]; // 2 * 48 bytes

    // a = 2 + 3i
    cuda_set_fp2_zero(a);
    lshift_fp(a[0], get_BLST_Rx()->p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], get_BLST_Rx()->p);      // Im(a) = 3

    // Serialize to bytes: convert from Montgomery to standard form
    vec384 temp_re, temp_im;
    from_fp(temp_re, a[0]);  // Convert real part from Montgomery
    from_fp(temp_im, a[1]);  // Convert imaginary part from Montgomery

    // Copy to serialized buffer
    memcpy(serialized, temp_re, 48);
    memcpy(serialized + 48, temp_im, 48);

    // Deserialize back: copy from buffer then convert to Montgomery form
    vec384 raw_re, raw_im;
    memcpy(raw_re, serialized, 48);
    memcpy(raw_im, serialized + 48, 48);

    // Convert to Montgomery form by multiplying to R^2
    mul_fp(b[0], raw_re, *get_BLST_RR());
    mul_fp(b[1], raw_im, *get_BLST_RR());

    // Should get back our original value
    results[0] = cuda_vec384x_equal(a, b);
}

// Batch processing tests for Fp2
GLOBAL_FUNC void cuda_test_fp2_batch_operations_kernel(bool* results, int batch_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= batch_size) return;

    vec384x a, b, result;

    // Create different test values for each thread
    cuda_set_fp2_zero(a);
    cuda_set_fp2_zero(b);

    // Set values based on thread index
    int re_a = (idx % 10) + 1;
    int im_a = ((idx + 3) % 7) + 1;
    int re_b = ((idx + 5) % 8) + 1;
    int im_b = ((idx + 7) % 6) + 1;

    cuda_set_fp2_value(a, re_a, im_a);
    cuda_set_fp2_value(b, re_b, im_b);

    // Perform various operations
    vec384x temp1, temp2, temp3;
    add_fp2(temp1, a, b);      // a + b
    mul_fp2(temp2, a, b);      // a * b
    sqr_fp2(temp3, a);         // a²

    // Combine operations: (a + b) * a²
    mul_fp2(result, temp1, temp3);

    // Simple validation - result should not be zero
    vec384x zero;
    cuda_set_fp2_zero(zero);
    results[idx] = !cuda_vec384x_equal(result, zero);
}

GLOBAL_FUNC void cuda_test_fp2_performance_kernel(int* operations_completed) {
    vec384x a, b, result;

    // Initialize with some values
    cuda_set_fp2_value(a, 5, 7);
    cuda_set_fp2_value(b, 3, 4);

    int ops = 0;
    // Perform many operations to test performance
    for (int i = 0; i < 500; i++) {  // Reduced from 1000 for Fp2 complexity
        mul_fp2(result, a, b);
        sqr_fp2(a, result);
        add_fp2(b, a, b);
        conjugate_fp2(result, result);
        ops += 4;
    }

    *operations_completed = ops;
}

// Helper function to run a CUDA kernel and check results
template<typename KernelFunc>
bool run_cuda_fp2_test(KernelFunc kernel, int num_results = 1) {
    bool* d_results, * h_results;

    // Allocate memory
    h_results = new bool[num_results];
    CUDA_CHECK(cudaMalloc(&d_results, num_results * sizeof(bool)));

    // Initialize results to false
    CUDA_CHECK(cudaMemset(d_results, 0, num_results * sizeof(bool)));

    // Launch kernel
    kernel<<<1, 1>>>(d_results);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_results, d_results, num_results * sizeof(bool), cudaMemcpyDeviceToHost));

    // Check all results
    bool all_passed = true;
    for (int i = 0; i < num_results; i++) {
        if (!h_results[i]) {
            all_passed = false;
            printf("Test failed at index %d\n", i);
            break;
        }
    }

    // Cleanup
    delete[] h_results;
    CUDA_CHECK(cudaFree(d_results));

    return all_passed;
}

// Global flag for CUDA availability
static bool cuda_device_available = false;

// CUDA Environment for runtime detection
class CUDAEnvironment : public ::testing::Environment {
public:
    void SetUp() override {
        int device_count = 0;
        cudaError_t error = cudaGetDeviceCount(&device_count);

        if (error != cudaSuccess || device_count == 0) {
            std::cout << "No CUDA devices available. CUDA Fp2 tests will be skipped.\n";
            cuda_device_available = false;
        }
        else {
            std::cout << "Found " << device_count << " CUDA device(s) for Fp2 tests\n";
            cuda_device_available = true;
        }
    }

    static bool IsCudaAvailable() { return cuda_device_available; }
};

// Initialize the test environment
struct Fp2CudaTestRegistration {
    Fp2CudaTestRegistration() {
        ::testing::AddGlobalTestEnvironment(new CUDAEnvironment());
    }
};

static Fp2CudaTestRegistration fp2_cuda_test_registration;

// CUDA test cases for Fp2 operations
TEST(Fp2OpsCUDA, ZeroAndOne) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp2_test(cuda_test_fp2_zero_and_one_kernel, 4));
}

TEST(Fp2OpsCUDA, Addition) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp2_test(cuda_test_fp2_addition_kernel, 2));
}

TEST(Fp2OpsCUDA, Subtraction) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp2_test(cuda_test_fp2_subtraction_kernel, 2));
}

TEST(Fp2OpsCUDA, Multiplication) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp2_test(cuda_test_fp2_multiplication_kernel, 2));
}

TEST(Fp2OpsCUDA, SquareAndMultiply) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp2_test(cuda_test_fp2_square_and_multiply_kernel, 1));
}

TEST(Fp2OpsCUDA, Inversion) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp2_test(cuda_test_fp2_inversion_kernel, 1));
}

TEST(Fp2OpsCUDA, Frobenius) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp2_test(cuda_test_fp2_frobenius_kernel, 1));
}

TEST(Fp2OpsCUDA, SquareRoot) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp2_test(cuda_test_fp2_square_root_kernel, 3));
}

TEST(Fp2OpsCUDA, Serialization) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp2_test(cuda_test_fp2_serialization_kernel, 1));
}

// Performance and batch processing tests
TEST(Fp2OpsCUDA, PerformanceBaseline) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    int* d_ops, h_ops;

    CUDA_CHECK(cudaMalloc(&d_ops, sizeof(int)));
    CUDA_CHECK(cudaMemset(d_ops, 0, sizeof(int)));

    // Record timing
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));

    CUDA_CHECK(cudaEventRecord(start));
    cuda_test_fp2_performance_kernel<<<1, 1>>>(d_ops);
    CUDA_CHECK(cudaEventRecord(stop));

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));

    CUDA_CHECK(cudaMemcpy(&h_ops, d_ops, sizeof(int), cudaMemcpyDeviceToHost));

    printf("CUDA Fp2 Performance: %d operations in %.3f ms (%.1f ops/ms)\n",
        h_ops, milliseconds, h_ops / milliseconds);

    CUDA_CHECK(cudaFree(d_ops));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));

    EXPECT_GT(h_ops, 0);
}

TEST(Fp2OpsCUDA, BatchProcessing) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const int batch_size = 512;  // Reduced for Fp2 complexity
    bool* d_results, * h_results;

    // Allocate memory
    h_results = new bool[batch_size];
    CUDA_CHECK(cudaMalloc(&d_results, batch_size * sizeof(bool)));
    CUDA_CHECK(cudaMemset(d_results, 0, batch_size * sizeof(bool)));

    // Launch kernel with multiple threads
    int threads_per_block = 128;  // Reduced for Fp2 complexity
    int blocks = (batch_size + threads_per_block - 1) / threads_per_block;

    cuda_test_fp2_batch_operations_kernel<<<blocks, threads_per_block>>>(d_results, batch_size);
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

    printf("Fp2 Batch processing: %d/%d operations successful\n", successful_operations, batch_size);

    // Cleanup
    delete[] h_results;
    CUDA_CHECK(cudaFree(d_results));

    // Expect at least 95% success rate
    EXPECT_GE(successful_operations, batch_size * 0.95);
}