#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <memory>
#include <random>
#include <cuda_runtime.h>

extern "C" {
#include "blst/fields/fields.h"
#include "blst/fields/Fp2_cuda.h"
#include "blst/exports.h"
}

#define LIMBS NLIMBS(384)

class Fp2CudaTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize random number generator
        gen.seed(42); // Fixed seed for reproducible tests

        // Initialize CUDA device
        cudaError_t err = cudaSetDevice(0);
        if (err != cudaSuccess) {
            GTEST_SKIP() << "CUDA device not available: " << cudaGetErrorString(err);
        }
    }

    void TearDown() override {
        // Reset CUDA device
        cudaDeviceReset();
    }

    // Helper: Set a vec384x to a small integer value (real + i*imaginary)
    void set_vec384x(vec384x v, uint64_t real, uint64_t imaginary) {
        for (int i = 0; i < LIMBS; ++i) {
            v[0][i] = 0;
            v[1][i] = 0;
        }
        memcpy(v[0], &real, sizeof(real));
        memcpy(v[1], &imaginary, sizeof(imaginary));
    }

    // Helper: Compare two vec384x for equality
    bool vec384x_equal(const vec384x a, const vec384x b) {
        for (int i = 0; i < LIMBS; ++i) {
            if (a[0][i] != b[0][i] || a[1][i] != b[1][i]) return false;
        }
        return true;
    }

    // Helper: Print vec384x for debugging
    void print_vec384x(const vec384x v) {
        printf("Real: ");
        for (int i = LIMBS - 1; i >= 0; --i)
            printf("%016llx%s", (unsigned long long)v[0][i], i ? " " : "\n");
        printf("Imag: ");
        for (int i = LIMBS - 1; i >= 0; --i)
            printf("%016llx%s", (unsigned long long)v[1][i], i ? " " : "\n");
    }

    // Helper: Copy one vec384x to another
    void fp2_copy(vec384x dst, const vec384x src) {
        for (int i = 0; i < LIMBS; ++i) {
            dst[0][i] = src[0][i];
            dst[1][i] = src[1][i];
        }
    }

    // Helper: Generate random Fp2 element
    void random_fp2(vec384x v) {
        // Generate random limbs for both real and imaginary parts
        for (int i = 0; i < LIMBS; ++i) {
            v[0][i] = gen();
            v[1][i] = gen();
        }
        // Reduce modulo to ensure valid elements
        if (LIMBS > 0) {
            v[0][LIMBS - 1] &= 0x7FFFFFFFFFFFFFFFULL;
            v[1][LIMBS - 1] &= 0x7FFFFFFFFFFFFFFFULL;
        }
    }

    // Helper: Convert standard integer to Montgomery form for Fp2
    void to_mont_fp2(vec384x r, const vec384x a) {
        to_fp(r[0], a[0]);
        to_fp(r[1], a[1]);
    }

    // Helper: Allocate device memory for vec384x arrays
    vec384x* allocate_device_fp2_array(size_t num_elements) {
        vec384x* d_ptr = nullptr;
        cudaError_t err = cudaMalloc((void**)&d_ptr, num_elements * sizeof(vec384x));
        if (err != cudaSuccess) {
            return nullptr;
        }
        return d_ptr;
    }

    // Helper: Free device memory
    void free_device_fp2_array(vec384x* d_ptr) {
        if (d_ptr) {
            cudaFree(d_ptr);
        }
    }

    // Helper: Copy host vec384x array to device
    batched_error_code copy_host_to_device_fp2(vec384x* d_dst, const vec384x* h_src, size_t num_elements) {
        if (!d_dst || !h_src) {
            return BATCHED_ERROR_NULL_POINTER;
        }

        cudaError_t err = cudaMemcpy(d_dst, h_src, num_elements * sizeof(vec384x), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
        }

        return BATCHED_SUCCESS;
    }

    // Helper: Copy device vec384x array to host
    batched_error_code copy_device_to_host_fp2(vec384x* h_dst, const vec384x* d_src, size_t num_elements) {
        if (!h_dst || !d_src) {
            return BATCHED_ERROR_NULL_POINTER;
        }

        cudaError_t err = cudaMemcpy(h_dst, d_src, num_elements * sizeof(vec384x), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
        }

        return BATCHED_SUCCESS;
    }

    // Helper: Compare vector results against CPU reference
    bool compare_fp2_vectors(const vec384x* vec_result, const vec384x* cpu_result, size_t num, const char* operation) {
        bool all_equal = true;
        for (size_t i = 0; i < num; i++) {
            if (!vec384x_equal(vec_result[i], cpu_result[i])) {
                printf("Mismatch in %s at index %zu:\n", operation, i);
                printf("Vector result: "); print_vec384x(vec_result[i]);
                printf("CPU result:    "); print_vec384x(cpu_result[i]);
                all_equal = false;
                // Only show first few mismatches to avoid flooding output
                if (i > 5) break;
            }
        }
        return all_equal;
    }

    // Helper: Check batched error code and provide meaningful error message
    void expect_success(batched_error_code result, const char* operation) {
        EXPECT_EQ(result, BATCHED_SUCCESS)
            << "Operation " << operation << " failed with error: "
            << batched_error_code_string(result);
    }

    std::mt19937_64 gen;
};

TEST_F(Fp2CudaTest, VectorAdd) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384x(h_a[i], 5 + i, 3 + i);
        set_vec384x(h_b[i], 7 + i * 2, 4 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        add_fp2(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_b = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp2(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = add_fp2_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "add_fp2_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "add_fp2_vec"));
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_b);
    free_device_fp2_array(d_result);
}

TEST_F(Fp2CudaTest, VectorSub) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384x(h_a[i], 10 + i * 2, 8 + i);
        set_vec384x(h_b[i], 3 + i, 2 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        sub_fp2(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_b = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp2(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = sub_fp2_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "sub_fp2_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "sub_fp2_vec"));
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_b);
    free_device_fp2_array(d_result);
}

TEST_F(Fp2CudaTest, VectorMul) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data with Montgomery form conversion
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        vec384x temp_a, temp_b;
        set_vec384x(temp_a, 3 + i, 2 + i);
        set_vec384x(temp_b, 4 + i, 5 + i);

        // Convert to Montgomery form
        to_mont_fp2(h_a[i], temp_a);
        to_mont_fp2(h_b[i], temp_b);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        mul_fp2(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_b = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp2(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = mul_fp2_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "mul_fp2_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "mul_fp2_vec"));
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_b);
    free_device_fp2_array(d_result);
}

TEST_F(Fp2CudaTest, VectorSqr) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data with Montgomery form conversion
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        vec384x temp_a;
        set_vec384x(temp_a, 5 + i, 3 + i);

        // Convert to Montgomery form
        to_mont_fp2(h_a[i], temp_a);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        sqr_fp2(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = sqr_fp2_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "sqr_fp2_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "sqr_fp2_vec"));
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_result);
}

TEST_F(Fp2CudaTest, VectorCNeg) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result_true = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result_true = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result_false = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result_false = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384x(h_a[i], 10 + i, 5 + i);
    }

    // CPU reference computation with flag = 1
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        cneg_fp2(cpu_result_true[i], h_a[i], 1);
    }

    // CPU reference computation with flag = 0
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        cneg_fp2(cpu_result_false[i], h_a[i], 0);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result_true = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result_false = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result_true, nullptr) << "Failed to allocate device memory for result_true";
    ASSERT_NE(d_result_false, nullptr) << "Failed to allocate device memory for result_false";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with flag = 1 and error checking
    batched_error_code result_true = cneg_fp2_vec(d_result_true, d_a, 1, NUM_ELEMENTS);
    expect_success(result_true, "cneg_fp2_vec (flag=1)");

    // Vector computation with flag = 0 and error checking
    batched_error_code result_false = cneg_fp2_vec(d_result_false, d_a, 0, NUM_ELEMENTS);
    expect_success(result_false, "cneg_fp2_vec (flag=0)");

    // Copy results back to host and compare only if operations succeeded
    if (result_true == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result_true.get(), d_result_true, NUM_ELEMENTS), "copy result_true to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result_true.get(), cpu_result_true.get(), NUM_ELEMENTS, "cneg_fp2_vec (flag=1)"));
    }
    if (result_false == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result_false.get(), d_result_false, NUM_ELEMENTS), "copy result_false to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result_false.get(), cpu_result_false.get(), NUM_ELEMENTS, "cneg_fp2_vec (flag=0)"));
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_result_true);
    free_device_fp2_array(d_result_false);
}

TEST_F(Fp2CudaTest, VectorMul3) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384x(h_a[i], 2 + i, 3 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        mul_by_3_fp(cpu_result[i][0], h_a[i][0]);
        mul_by_3_fp(cpu_result[i][1], h_a[i][1]);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = mul3_fp2_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "mul3_fp2_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "mul3_fp2_vec"));
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_result);
}

TEST_F(Fp2CudaTest, VectorMul8) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384x(h_a[i], 2 + i, 3 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        mul_by_8_fp(cpu_result[i][0], h_a[i][0]);
        mul_by_8_fp(cpu_result[i][1], h_a[i][1]);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = mul8_fp2_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "mul8_fp2_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "mul8_fp2_vec"));
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_result);
}

TEST_F(Fp2CudaTest, VectorLshift) {
    const size_t NUM_ELEMENTS = 1000;
    const size_t shift_count = 2;
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384x(h_a[i], 1 + i, 2 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        lshift_fp(cpu_result[i][0], h_a[i][0], shift_count);
        lshift_fp(cpu_result[i][1], h_a[i][1], shift_count);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = lshift_fp2_vec(d_result, d_a, shift_count, NUM_ELEMENTS);
    expect_success(result, "lshift_fp2_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "lshift_fp2_vec"));
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_result);
}

TEST_F(Fp2CudaTest, VectorConjugate) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384x(h_a[i], 5 + i, 7 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        conjugate_fp2(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = conjugate_fp2_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "conjugate_fp2_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "conjugate_fp2_vec"));
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_result);
}

TEST_F(Fp2CudaTest, VectorInverse) {
    const size_t NUM_ELEMENTS = 100; // Smaller for inverse operations
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data with Montgomery form conversion
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        vec384x temp_a;
        set_vec384x(temp_a, 5 + i, 3 + i); // Non-zero values

        // Convert to Montgomery form
        to_mont_fp2(h_a[i], temp_a);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        reciprocal_fp2(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = inverse_fp2_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "inverse_fp2_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "inverse_fp2_vec"));
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_result);
}

TEST_F(Fp2CudaTest, VectorSqrt) {
    const size_t NUM_ELEMENTS = 100; // Smaller for sqrt operations
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data with perfect squares in Montgomery form
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        vec384x temp_base, temp_squared;
        set_vec384x(temp_base, 2 + i, 1 + i);

        // Convert to Montgomery form and square to ensure it has a square root
        to_mont_fp2(temp_base, temp_base);
        sqr_fp2(h_a[i], temp_base);
    }

    // CPU reference computation
    bool all_sqrt_success = true;
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        bool sqrt_success = sqrt_fp2(cpu_result[i], h_a[i]);
        if (!sqrt_success) {
            all_sqrt_success = false;
            break;
        }
    }

    if (!all_sqrt_success) {
        GTEST_SKIP() << "Some CPU square root operations failed";
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = sqrt_fp2_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "sqrt_fp2_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");

        // For square root, we need to verify that the result squared gives the original
        // rather than exact equality (since there are two square roots)
        bool all_correct = true;
        for (size_t i = 0; i < NUM_ELEMENTS; i++) {
            vec384x squared_result;
            sqr_fp2(squared_result, h_result[i]);
            if (!vec384x_equal(squared_result, h_a[i])) {
                all_correct = false;
                printf("Square root verification failed at index %zu\n", i);
                break;
            }
        }
        EXPECT_TRUE(all_correct) << "Square root verification failed";
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_result);
}

TEST_F(Fp2CudaTest, VectorErrorHandling) {
    const size_t NUM_ELEMENTS = 100;

    // Test null pointer handling - these should fail immediately without device allocation
    batched_error_code null_out_result = add_fp2_vec(nullptr, nullptr, nullptr, NUM_ELEMENTS);
    EXPECT_EQ(null_out_result, BATCHED_ERROR_NULL_POINTER)
        << "Expected NULL_POINTER error for null output array";

    // Test zero size handling - zero size should be OK
    vec384x* d_dummy = allocate_device_fp2_array(NUM_ELEMENTS);
    ASSERT_NE(d_dummy, nullptr) << "Failed to allocate device memory for dummy array";

    batched_error_code zero_size_result = add_fp2_vec(d_dummy, d_dummy, d_dummy, 0);
    EXPECT_EQ(zero_size_result, BATCHED_SUCCESS)
        << "Expected SUCCESS for zero element count (zero-size operations should be OK)";

    free_device_fp2_array(d_dummy);
}

TEST_F(Fp2CudaTest, VectorMixedOperations) {
    const size_t MIXED_SIZE = 256;
    auto h_a = std::make_unique<vec384x[]>(MIXED_SIZE);
    auto h_b = std::make_unique<vec384x[]>(MIXED_SIZE);
    auto h_result = std::make_unique<vec384x[]>(MIXED_SIZE);
    auto cpu_result = std::make_unique<vec384x[]>(MIXED_SIZE);

    // Initialize test data
    for (size_t i = 0; i < MIXED_SIZE; i++) {
        set_vec384x(h_a[i], 3 + i, 2 + i);
        set_vec384x(h_b[i], 5 + i * 2, 4 + i);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(MIXED_SIZE);
    vec384x* d_b = allocate_device_fp2_array(MIXED_SIZE);
    vec384x* d_result = allocate_device_fp2_array(MIXED_SIZE);
    vec384x* d_temp = allocate_device_fp2_array(MIXED_SIZE);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";
    ASSERT_NE(d_temp, nullptr) << "Failed to allocate device memory for temp";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), MIXED_SIZE), "copy h_a to device");
    expect_success(copy_host_to_device_fp2(d_b, h_b.get(), MIXED_SIZE), "copy h_b to device");

    // Perform a sequence of operations: (a + b) * 3 then conjugate
    // Step 1: Add
    batched_error_code add_result = add_fp2_vec(d_temp, d_a, d_b, MIXED_SIZE);
    expect_success(add_result, "mixed operations add");

    // Step 2: Multiply by 3
    if (add_result == BATCHED_SUCCESS) {
        batched_error_code mul3_result = mul3_fp2_vec(d_result, d_temp, MIXED_SIZE);
        expect_success(mul3_result, "mixed operations mul3");

        // Step 3: Conjugate
        if (mul3_result == BATCHED_SUCCESS) {
            batched_error_code conj_result = conjugate_fp2_vec(d_result, d_result, MIXED_SIZE);
            expect_success(conj_result, "mixed operations conjugate");

            if (conj_result == BATCHED_SUCCESS) {
                expect_success(copy_device_to_host_fp2(h_result.get(), d_result, MIXED_SIZE), "copy mixed result to host");

                // Calculate expected results on host
                for (size_t i = 0; i < MIXED_SIZE; i++) {
                    vec384x temp_sum, temp_mul3;
                    add_fp2(temp_sum, h_a[i], h_b[i]);
                    mul_by_3_fp(temp_mul3[0], temp_sum[0]);
                    mul_by_3_fp(temp_mul3[1], temp_sum[1]);
                    conjugate_fp2(cpu_result[i], temp_mul3);
                }

                EXPECT_TRUE(compare_fp2_vectors(h_result.get(), cpu_result.get(), MIXED_SIZE, "mixed operations"));
            }
        }
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_b);
    free_device_fp2_array(d_result);
    free_device_fp2_array(d_temp);
}

TEST_F(Fp2CudaTest, VectorPerformanceBaseline) {
    const size_t NUM_ELEMENTS = 5000; // Smaller due to Fp2 complexity
    auto h_a = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384x[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        random_fp2(h_a[i]);
        random_fp2(h_b[i]);
    }

    // This test primarily validates correctness at scale
    // Performance measurements would require timing infrastructure

    // CPU computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        add_fp2(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    vec384x* d_a = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_b = allocate_device_fp2_array(NUM_ELEMENTS);
    vec384x* d_result = allocate_device_fp2_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp2(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp2(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = add_fp2_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "performance baseline");

    // Verify correctness only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp2(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp2_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "performance baseline"));
    }

    // Cleanup
    free_device_fp2_array(d_a);
    free_device_fp2_array(d_b);
    free_device_fp2_array(d_result);
}