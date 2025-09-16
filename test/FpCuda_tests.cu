#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <memory>
#include <random>
#include <cuda_runtime.h>

extern "C" {
#include "blst/fields/fields.h"
#include "blst/fields/fields_cuda.h"
#include "blst/exports.h"
}

#define LIMBS NLIMBS(384)

class FpCudaTest : public ::testing::Test {
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

    // Helper: Set a vec384 to a small integer value
    void set_vec384(vec384 v, uint64_t x) {
        for (int i = 0; i < LIMBS; ++i)
            v[i] = 0;
        memcpy(v, &x, sizeof(x));
    }

    // Helper: Compare two vec384 for equality
    bool vec384_equal(const vec384 a, const vec384 b) {
        for (int i = 0; i < LIMBS; ++i)
            if (a[i] != b[i]) return false;
        return true;
    }

    // Helper: Print vec384 for debugging
    void print_vec384(const vec384 v) {
        for (int i = LIMBS - 1; i >= 0; --i)
            printf("%016llx%s", (unsigned long long)v[i], i ? " " : "\n");
    }

    // Helper: Copy one vec384 to another
    void fp_copy(vec384 dst, const vec384 src) {
        for (int i = 0; i < LIMBS; ++i)
            dst[i] = src[i];
    }

    // Helper: Generate random Fp element
    void random_fp(vec384 v) {
        // Generate random limbs
        for (int i = 0; i < LIMBS; ++i) {
            v[i] = gen();
        }
        // Reduce modulo Fp to ensure valid element
        // Note: This is a simplified approach; real implementation would use proper modular reduction
        if (LIMBS > 0) {
            v[LIMBS - 1] &= 0x7FFFFFFFFFFFFFFFULL; // Ensure MSB is clear to avoid overflow
        }
    }

    // Helper: Convert standard integer to Montgomery form
    void to_mont_fp(vec384 r, const vec384 a) {
        to_fp(r, a);
    }

    // Helper: Allocate device memory for vec384 arrays
    vec384* allocate_device_fp_array(size_t num_elements) {
        vec384* d_ptr = nullptr;
        cudaError_t err = cudaMalloc((void**)&d_ptr, num_elements * sizeof(vec384));
        if (err != cudaSuccess) {
            return nullptr;
        }
        return d_ptr;
    }

    // Helper: Free device memory
    void free_device_fp_array(vec384* d_ptr) {
        if (d_ptr) {
            cudaFree(d_ptr);
        }
    }

    // Helper: Copy host vec384 array to device
    batched_error_code copy_host_to_device_fp(vec384* d_dst, const vec384* h_src, size_t num_elements) {
        if (!d_dst || !h_src) {
            return BATCHED_ERROR_NULL_POINTER;
        }

        cudaError_t err = cudaMemcpy(d_dst, h_src, num_elements * sizeof(vec384), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
        }

        return BATCHED_SUCCESS;
    }

    // Helper: Copy device vec384 array to host
    batched_error_code copy_device_to_host_fp(vec384* h_dst, const vec384* d_src, size_t num_elements) {
        if (!h_dst || !d_src) {
            return BATCHED_ERROR_NULL_POINTER;
        }

        cudaError_t err = cudaMemcpy(h_dst, d_src, num_elements * sizeof(vec384), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
        }

        return BATCHED_SUCCESS;
    }

    // Helper: Compare vector results against CPU reference
    bool compare_fp_vectors(const vec384* vec_result, const vec384* cpu_result, size_t num, const char* operation) {
        bool all_equal = true;
        for (size_t i = 0; i < num; i++) {
            if (!vec384_equal(vec_result[i], cpu_result[i])) {
                printf("Mismatch in %s at index %zu:\n", operation, i);
                printf("Vector result: "); print_vec384(vec_result[i]);
                printf("CPU result:    "); print_vec384(cpu_result[i]);
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

TEST_F(FpCudaTest, VectorAdd) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_a[i], 5 + i);
        set_vec384(h_b[i], 7 + i * 2);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        add_fp(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_b = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = add_fp_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "add_fp_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "add_fp_vec"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_b);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorSub) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_a[i], 10 + i * 2);
        set_vec384(h_b[i], 3 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        sub_fp(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_b = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = sub_fp_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "sub_fp_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "sub_fp_vec"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_b);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorMul) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data with Montgomery form conversion
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        vec384 temp_a, temp_b;
        set_vec384(temp_a, 3 + i);
        set_vec384(temp_b, 4 + i);

        // Convert to Montgomery form
        to_mont_fp(h_a[i], temp_a);
        to_mont_fp(h_b[i], temp_b);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        mul_fp(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_b = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = mul_fp_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "mul_fp_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "mul_fp_vec"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_b);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorSqr) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data with Montgomery form conversion
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        vec384 temp_a;
        set_vec384(temp_a, 5 + i);

        // Convert to Montgomery form
        to_mont_fp(h_a[i], temp_a);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        sqr_fp(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = sqr_fp_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "sqr_fp_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "sqr_fp_vec"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorCNeg) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result_true = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result_true = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result_false = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result_false = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_a[i], 10 + i);
    }

    // CPU reference computation with flag = 1
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        cneg_fp(cpu_result_true[i], h_a[i], 1);
    }

    // CPU reference computation with flag = 0
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        cneg_fp(cpu_result_false[i], h_a[i], 0);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result_true = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result_false = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result_true, nullptr) << "Failed to allocate device memory for result_true";
    ASSERT_NE(d_result_false, nullptr) << "Failed to allocate device memory for result_false";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with flag = 1 and error checking
    batched_error_code result_true = cneg_fp_vec(d_result_true, d_a, 1, NUM_ELEMENTS);
    expect_success(result_true, "cneg_fp_vec (flag=1)");

    // Vector computation with flag = 0 and error checking
    batched_error_code result_false = cneg_fp_vec(d_result_false, d_a, 0, NUM_ELEMENTS);
    expect_success(result_false, "cneg_fp_vec (flag=0)");

    // Copy results back to host and compare only if operations succeeded
    if (result_true == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result_true.get(), d_result_true, NUM_ELEMENTS), "copy result_true to host");
        EXPECT_TRUE(compare_fp_vectors(h_result_true.get(), cpu_result_true.get(), NUM_ELEMENTS, "cneg_fp_vec (flag=1)"));
    }
    if (result_false == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result_false.get(), d_result_false, NUM_ELEMENTS), "copy result_false to host");
        EXPECT_TRUE(compare_fp_vectors(h_result_false.get(), cpu_result_false.get(), NUM_ELEMENTS, "cneg_fp_vec (flag=0)"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_result_true);
    free_device_fp_array(d_result_false);
}

TEST_F(FpCudaTest, VectorMul3) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_a[i], 2 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        mul_by_3_fp(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = mul3_fp_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "mul3_fp_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "mul3_fp_vec"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorMul8) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_a[i], 2 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        mul_by_8_fp(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = mul8_fp_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "mul8_fp_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "mul8_fp_vec"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorDiv2) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_a[i], 8 + i * 2);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        div_by_2_fp(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = div2_fp_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "div2_fp_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "div2_fp_vec"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorLshift) {
    const size_t NUM_ELEMENTS = 1000;
    const size_t shift_count = 2;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_a[i], 1 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        lshift_fp(cpu_result[i], h_a[i], shift_count);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = lshift_fp_vec(d_result, d_a, shift_count, NUM_ELEMENTS);
    expect_success(result, "lshift_fp_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "lshift_fp_vec"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorToFp) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_a[i], 9 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        to_fp(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = to_fp_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "to_fp_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "to_fp_vec"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorFromFp) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data in Montgomery form
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        vec384 temp;
        set_vec384(temp, 9 + i);
        to_fp(h_a[i], temp);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        from_fp(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = from_fp_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "from_fp_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "from_fp_vec"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorMontgomeryRoundTrip) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_a[i], 42 + i);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Convert to Montgomery form
    batched_error_code to_result = to_fp_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(to_result, "to_fp_vec");

    // Convert back from Montgomery form
    if (to_result == BATCHED_SUCCESS) {
        batched_error_code from_result = from_fp_vec(d_result, d_result, NUM_ELEMENTS);
        expect_success(from_result, "from_fp_vec");

        if (from_result == BATCHED_SUCCESS) {
            expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
            EXPECT_TRUE(compare_fp_vectors(h_result.get(), h_a.get(), NUM_ELEMENTS, "round-trip conversion"));
        }
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorMixedOperations) {
    const size_t MIXED_SIZE = 256;
    auto h_a = std::make_unique<vec384[]>(MIXED_SIZE);
    auto h_b = std::make_unique<vec384[]>(MIXED_SIZE);
    auto h_result = std::make_unique<vec384[]>(MIXED_SIZE);
    auto cpu_result = std::make_unique<vec384[]>(MIXED_SIZE);

    // Initialize test data
    for (size_t i = 0; i < MIXED_SIZE; i++) {
        set_vec384(h_a[i], 3 + i);
        set_vec384(h_b[i], 5 + i * 2);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(MIXED_SIZE);
    vec384* d_b = allocate_device_fp_array(MIXED_SIZE);
    vec384* d_result = allocate_device_fp_array(MIXED_SIZE);
    vec384* d_temp = allocate_device_fp_array(MIXED_SIZE);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";
    ASSERT_NE(d_temp, nullptr) << "Failed to allocate device memory for temp";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), MIXED_SIZE), "copy h_a to device");
    expect_success(copy_host_to_device_fp(d_b, h_b.get(), MIXED_SIZE), "copy h_b to device");

    // Perform a sequence of operations: (a + b) * 3 / 2
    // Step 1: Add
    batched_error_code add_result = add_fp_vec(d_temp, d_a, d_b, MIXED_SIZE);
    expect_success(add_result, "mixed operations add");

    // Step 2: Multiply by 3
    if (add_result == BATCHED_SUCCESS) {
        batched_error_code mul3_result = mul3_fp_vec(d_result, d_temp, MIXED_SIZE);
        expect_success(mul3_result, "mixed operations mul3");

        // Step 3: Divide by 2
        if (mul3_result == BATCHED_SUCCESS) {
            batched_error_code div2_result = div2_fp_vec(d_result, d_result, MIXED_SIZE);
            expect_success(div2_result, "mixed operations div2");

            if (div2_result == BATCHED_SUCCESS) {
                expect_success(copy_device_to_host_fp(h_result.get(), d_result, MIXED_SIZE), "copy mixed result to host");

                // Calculate expected results on host
                for (size_t i = 0; i < MIXED_SIZE; i++) {
                    vec384 temp_sum, temp_mul3;
                    add_fp(temp_sum, h_a[i], h_b[i]);
                    mul_by_3_fp(temp_mul3, temp_sum);
                    div_by_2_fp(cpu_result[i], temp_mul3);
                }

                EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), MIXED_SIZE, "mixed operations"));
            }
        }
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_b);
    free_device_fp_array(d_result);
    free_device_fp_array(d_temp);
}

TEST_F(FpCudaTest, VectorEdgeCases) {
    const size_t NUM_ELEMENTS = 10;

    // Test with zero elements
    auto h_zeros = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_zeros[i], 0);
    }

    // Allocate device memory
    vec384* d_zeros = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_zeros, nullptr) << "Failed to allocate device memory for zeros";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_zeros, h_zeros.get(), NUM_ELEMENTS), "copy zeros to device");

    // Test addition with zeros
    batched_error_code add_result = add_fp_vec(d_result, d_zeros, d_zeros, NUM_ELEMENTS);
    expect_success(add_result, "add_fp_vec with zeros");

    if (add_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy add result to host");
        for (size_t i = 0; i < NUM_ELEMENTS; i++) {
            EXPECT_TRUE(vec384_equal(h_result[i], h_zeros[i])) << "0 + 0 should equal 0 at index " << i;
        }
    }

    // Test with one
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_zeros[i], 1);
    }

    expect_success(copy_host_to_device_fp(d_zeros, h_zeros.get(), NUM_ELEMENTS), "copy ones to device");

    batched_error_code add_one_result = add_fp_vec(d_result, d_zeros, d_zeros, NUM_ELEMENTS);
    expect_success(add_one_result, "add_fp_vec with ones");

    if (add_one_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy add ones result to host");
        for (size_t i = 0; i < NUM_ELEMENTS; i++) {
            vec384 expected;
            set_vec384(expected, 2);
            EXPECT_TRUE(vec384_equal(h_result[i], expected)) << "1 + 1 should equal 2 at index " << i;
        }
    }

    // Cleanup
    free_device_fp_array(d_zeros);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorErrorHandling) {
    const size_t NUM_ELEMENTS = 100;

    // Test null pointer handling - these should fail immediately without device allocation
    batched_error_code null_out_result = add_fp_vec(nullptr, nullptr, nullptr, NUM_ELEMENTS);
    EXPECT_EQ(null_out_result, BATCHED_ERROR_NULL_POINTER)
        << "Expected NULL_POINTER error for null output array";

    // Test zero size handling - zero size should be OK
    vec384* d_dummy = allocate_device_fp_array(NUM_ELEMENTS);
    ASSERT_NE(d_dummy, nullptr) << "Failed to allocate device memory for dummy array";

    batched_error_code zero_size_result = add_fp_vec(d_dummy, d_dummy, d_dummy, 0);
    EXPECT_EQ(zero_size_result, BATCHED_SUCCESS)
        << "Expected SUCCESS for zero element count (zero-size operations should be OK)";

    free_device_fp_array(d_dummy);
}

TEST_F(FpCudaTest, VectorSmallValues) {
    const size_t NUM_ELEMENTS = 5;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Test with small known values
    const uint64_t TEST_VALUES[] = { 1, 2, 5, 10, 42 };

    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_vec384(h_a[i], TEST_VALUES[i]);
        set_vec384(h_b[i], TEST_VALUES[(i + 1) % NUM_ELEMENTS]);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_b = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Test addition
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        add_fp(cpu_result[i], h_a[i], h_b[i]);
    }
    batched_error_code add_result = add_fp_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(add_result, "small values add");
    if (add_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy add result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "small values add"));
    }

    // Test subtraction
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        sub_fp(cpu_result[i], h_a[i], h_b[i]);
    }
    batched_error_code sub_result = sub_fp_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(sub_result, "small values sub");
    if (sub_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy sub result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "small values sub"));
    }

    // Test multiplication by 3
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        mul_by_3_fp(cpu_result[i], h_a[i]);
    }
    batched_error_code mul3_result = mul3_fp_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(mul3_result, "small values mul3");
    if (mul3_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy mul3 result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "small values mul3"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_b);
    free_device_fp_array(d_result);
}

TEST_F(FpCudaTest, VectorPerformanceBaseline) {
    const size_t NUM_ELEMENTS = 10000;
    auto h_a = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        random_fp(h_a[i]);
        random_fp(h_b[i]);
    }

    // This test primarily validates correctness at scale
    // Performance measurements would require timing infrastructure

    // CPU computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        add_fp(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    vec384* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_b = allocate_device_fp_array(NUM_ELEMENTS);
    vec384* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = add_fp_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "performance baseline");

    // Verify correctness only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "performance baseline"));
    }

    // Cleanup
    free_device_fp_array(d_a);
    free_device_fp_array(d_b);
    free_device_fp_array(d_result);
}