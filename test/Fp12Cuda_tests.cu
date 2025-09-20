#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <memory>
#include <random>
#include <cuda_runtime.h>

#include "blst/fields/Fp12.h"
#include "blst/fields/Fp12_cuda.h"
#include "blst/exports.h"

using namespace blst;

class Fp12CudaTest : public ::testing::Test {
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

    // Helper: Set a vec384fp12 to a small integer value (only c0[0][0] component)
    void set_fp12_simple(vec384fp12 v, uint64_t x) {
        // Clear entire Fp12 element
        memset(v, 0, sizeof(vec384fp12));

        // Set only the first component to the integer value
        for (int i = 0; i < NLIMBS(384); ++i)
            v[0][0][0][i] = 0;
        memcpy(v[0][0][0], &x, sizeof(x));

        // Convert to Montgomery form
        to_fp(v[0][0][0], v[0][0][0]);
    }

    // Helper: Set vec384fp12 to one
    void set_fp12_one(vec384fp12 v) {
        blst_fp12_one_copy(v);
    }

    // Helper: Set vec384fp12 to zero  
    void set_fp12_zero(vec384fp12 v) {
        memset(v, 0, sizeof(vec384fp12));
    }

    // Helper: Compare two vec384fp12 for equality
    bool fp12_equal(const vec384fp12 a, const vec384fp12 b) {
        return blst_fp12_is_equal(a, b);
    }

    // Helper: Print vec384fp12 for debugging (simplified)
    void print_fp12(const vec384fp12 v) {
        printf("Fp12(c0[0][0]: ");
        for (int i = NLIMBS(384) - 1; i >= 0; --i)
            printf("%016llx%s", (unsigned long long)v[0][0][0][i], i ? " " : "");
        printf(", is_one: %s)\n", blst_fp12_is_one(v) ? "true" : "false");
    }

    // Helper: Copy one vec384fp12 to another
    void fp12_copy(vec384fp12 dst, const vec384fp12 src) {
        memcpy(dst, src, sizeof(vec384fp12));
    }

    // Helper: Generate random Fp12 element
    void random_fp12(vec384fp12 v) {
        // Generate random components for all 12 Fp elements
        for (int i = 0; i < 2; i++) {          // Fp12 components
            for (int j = 0; j < 3; j++) {      // Fp6 components
                for (int k = 0; k < 2; k++) {  // Fp2 components
                    for (int l = 0; l < NLIMBS(384); l++) {
                        v[i][j][k][l] = gen();
                    }
                    // Ensure MSB is clear to avoid overflow
                    if (NLIMBS(384) > 0) {
                        v[i][j][k][NLIMBS(384) - 1] &= 0x7FFFFFFFFFFFFFFFULL;
                    }
                    // Convert to Montgomery form
                    to_fp(v[i][j][k], v[i][j][k]);
                }
            }
        }
    }

    // Helper: Allocate device memory for vec384fp12 arrays
    vec384fp12* allocate_device_fp12_array(size_t num_elements) {
        vec384fp12* d_ptr = nullptr;
        cudaError_t err = cudaMalloc((void**)&d_ptr, num_elements * sizeof(vec384fp12));
        if (err != cudaSuccess) {
            return nullptr;
        }
        return d_ptr;
    }

    // Helper: Free device memory
    void free_device_fp12_array(vec384fp12* d_ptr) {
        if (d_ptr) {
            cudaFree(d_ptr);
        }
    }

    // Helper: Copy host vec384fp12 array to device
    batched_error_code copy_host_to_device_fp12(vec384fp12* d_dst, const vec384fp12* h_src, size_t num_elements) {
        if (!d_dst || !h_src) {
            return BATCHED_ERROR_NULL_POINTER;
        }

        cudaError_t err = cudaMemcpy(d_dst, h_src, num_elements * sizeof(vec384fp12), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
        }

        return BATCHED_SUCCESS;
    }

    // Helper: Copy device vec384fp12 array to host
    batched_error_code copy_device_to_host_fp12(vec384fp12* h_dst, const vec384fp12* d_src, size_t num_elements) {
        if (!h_dst || !d_src) {
            return BATCHED_ERROR_NULL_POINTER;
        }

        cudaError_t err = cudaMemcpy(h_dst, d_src, num_elements * sizeof(vec384fp12), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
        }

        return BATCHED_SUCCESS;
    }

    // Helper: Compare vector results against CPU reference
    bool compare_fp12_vectors(const vec384fp12* vec_result, const vec384fp12* cpu_result, size_t num, const char* operation) {
        bool all_equal = true;
        for (size_t i = 0; i < num; i++) {
            if (!fp12_equal(vec_result[i], cpu_result[i])) {
                printf("Mismatch in %s at index %zu:\n", operation, i);
                printf("Vector result: "); print_fp12(vec_result[i]);
                printf("CPU result:    "); print_fp12(cpu_result[i]);
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

TEST_F(Fp12CudaTest, VectorMul) {
    const size_t NUM_ELEMENTS = 100;  // Reduced for Fp12 complexity
    auto h_a = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        if (i < 5) {
            // Include some simple cases
            set_fp12_simple(h_a[i], 3 + i);
            set_fp12_simple(h_b[i], 4 + i);
        }
        else {
            random_fp12(h_a[i]);
            random_fp12(h_b[i]);
        }
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        blst_fp12_mul(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_b = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_result = allocate_device_fp12_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp12(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = mul_fp12_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "mul_fp12_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp12(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp12_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "mul_fp12_vec"));
    }

    // Cleanup
    free_device_fp12_array(d_a);
    free_device_fp12_array(d_b);
    free_device_fp12_array(d_result);
}

TEST_F(Fp12CudaTest, VectorSqr) {
    const size_t NUM_ELEMENTS = 100;
    auto h_a = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        if (i < 5) {
            set_fp12_simple(h_a[i], 5 + i);
        }
        else {
            random_fp12(h_a[i]);
        }
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        blst_fp12_sqr(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_result = allocate_device_fp12_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = sqr_fp12_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "sqr_fp12_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp12(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp12_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "sqr_fp12_vec"));
    }

    // Cleanup
    free_device_fp12_array(d_a);
    free_device_fp12_array(d_result);
}

TEST_F(Fp12CudaTest, VectorCyclotomicSqr) {
    const size_t NUM_ELEMENTS = 100;
    auto h_a = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);

    // Generate test data - for cyclotomic squaring, elements should be in cyclotomic subgroup
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        if (i < 5) {
            set_fp12_simple(h_a[i], 2 + i);
        }
        else {
            random_fp12(h_a[i]);
        }
        // Ensure element is in cyclotomic subgroup by conjugating and multiplying
        vec384fp12 temp;
        fp12_copy(temp, h_a[i]);
        blst_fp12_conjugate(temp);
        blst_fp12_mul(h_a[i], h_a[i], temp);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        blst_fp12_cyclotomic_sqr(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_result = allocate_device_fp12_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = cyclotomic_sqr_fp12_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "cyclotomic_sqr_fp12_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp12(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp12_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "cyclotomic_sqr_fp12_vec"));
    }

    // Cleanup
    free_device_fp12_array(d_a);
    free_device_fp12_array(d_result);
}

TEST_F(Fp12CudaTest, VectorInverse) {
    const size_t NUM_ELEMENTS = 50; // Further reduced for expensive inverse operation
    auto h_a = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);

    // Generate test data - avoid zero elements for inversion
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        if (i < 3) {
            set_fp12_simple(h_a[i], 7 + i); // Small non-zero values
        }
        else {
            random_fp12(h_a[i]);
        }
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        blst_fp12_inverse(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_result = allocate_device_fp12_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = inverse_fp12_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "inverse_fp12_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp12(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp12_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "inverse_fp12_vec"));
    }

    // Cleanup
    free_device_fp12_array(d_a);
    free_device_fp12_array(d_result);
}

TEST_F(Fp12CudaTest, VectorConjugate) {
    const size_t NUM_ELEMENTS = 100;
    auto h_a = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        if (i < 5) {
            set_fp12_simple(h_a[i], 3 + i);
        }
        else {
            random_fp12(h_a[i]);
        }
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        fp12_copy(cpu_result[i], h_a[i]);
        blst_fp12_conjugate(cpu_result[i]);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_result = allocate_device_fp12_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = conjugate_fp12_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "conjugate_fp12_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp12(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp12_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "conjugate_fp12_vec"));
    }

    // Cleanup
    free_device_fp12_array(d_a);
    free_device_fp12_array(d_result);
}

TEST_F(Fp12CudaTest, VectorFrobeniusMap) {
    const size_t NUM_ELEMENTS = 100;
    const size_t power = 1; // Test with power = 1
    auto h_a = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        if (i < 5) {
            set_fp12_simple(h_a[i], 4 + i);
        }
        else {
            random_fp12(h_a[i]);
        }
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        blst_fp12_frobenius_map(cpu_result[i], h_a[i], power);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_result = allocate_device_fp12_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = frobenius_map_fp12_vec(d_result, d_a, power, NUM_ELEMENTS);
    expect_success(result, "frobenius_map_fp12_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp12(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp12_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "frobenius_map_fp12_vec"));
    }

    // Cleanup
    free_device_fp12_array(d_a);
    free_device_fp12_array(d_result);
}

TEST_F(Fp12CudaTest, VectorEquality) {
    const size_t NUM_ELEMENTS = 10;
    auto h_a = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<bool[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<bool[]>(NUM_ELEMENTS);

    // Generate test data - some equal, some different
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_fp12_simple(h_a[i], 5 + i);
        if (i % 2 == 0) {
            fp12_copy(h_b[i], h_a[i]); // Equal
        }
        else {
            set_fp12_simple(h_b[i], 10 + i); // Different
        }
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        cpu_result[i] = blst_fp12_is_equal(h_a[i], h_b[i]);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_b = allocate_device_fp12_array(NUM_ELEMENTS);
    bool* d_result = nullptr;
    cudaMalloc(&d_result, NUM_ELEMENTS * sizeof(bool));

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp12(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = is_equal_fp12_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "is_equal_fp12_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        cudaMemcpy(h_result.get(), d_result, NUM_ELEMENTS * sizeof(bool), cudaMemcpyDeviceToHost);

        bool all_correct = true;
        for (size_t i = 0; i < NUM_ELEMENTS; i++) {
            if (h_result[i] != cpu_result[i]) {
                printf("Equality test failed at index %zu: got %s, expected %s\n",
                    i, h_result[i] ? "true" : "false", cpu_result[i] ? "true" : "false");
                all_correct = false;
            }
        }
        EXPECT_TRUE(all_correct);
    }

    // Cleanup
    free_device_fp12_array(d_a);
    free_device_fp12_array(d_b);
    cudaFree(d_result);
}

TEST_F(Fp12CudaTest, VectorIsOne) {
    const size_t NUM_ELEMENTS = 10;
    auto h_a = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<bool[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<bool[]>(NUM_ELEMENTS);

    // Generate test data - some ones, some non-ones
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        if (i % 3 == 0) {
            set_fp12_one(h_a[i]); // One
        }
        else {
            set_fp12_simple(h_a[i], 2 + i); // Non-one
        }
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        cpu_result[i] = blst_fp12_is_one(h_a[i]);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(NUM_ELEMENTS);
    bool* d_result = nullptr;
    cudaMalloc(&d_result, NUM_ELEMENTS * sizeof(bool));

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = is_one_fp12_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "is_one_fp12_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        cudaMemcpy(h_result.get(), d_result, NUM_ELEMENTS * sizeof(bool), cudaMemcpyDeviceToHost);

        bool all_correct = true;
        for (size_t i = 0; i < NUM_ELEMENTS; i++) {
            if (h_result[i] != cpu_result[i]) {
                printf("IsOne test failed at index %zu: got %s, expected %s\n",
                    i, h_result[i] ? "true" : "false", cpu_result[i] ? "true" : "false");
                all_correct = false;
            }
        }
        EXPECT_TRUE(all_correct);
    }

    // Cleanup
    free_device_fp12_array(d_a);
    cudaFree(d_result);
}

TEST_F(Fp12CudaTest, VectorSerialization) {
    const size_t NUM_ELEMENTS = 20;
    auto h_a = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_serialized = std::make_unique<unsigned char[]>(NUM_ELEMENTS * 12 * 48);
    auto h_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto cpu_serialized = std::make_unique<unsigned char[]>(NUM_ELEMENTS * 12 * 48);
    auto cpu_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        if (i < 3) {
            set_fp12_simple(h_a[i], 6 + i);
        }
        else {
            random_fp12(h_a[i]);
        }
    }

    // CPU reference computation - serialize and deserialize
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        blst_bendian_from_fp12(&cpu_serialized[i * 12 * 48], h_a[i]);
        blst_fp12_from_bendian(cpu_result[i], &cpu_serialized[i * 12 * 48]);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(NUM_ELEMENTS);
    unsigned char* d_serialized = nullptr;
    vec384fp12* d_result = allocate_device_fp12_array(NUM_ELEMENTS);
    cudaMalloc(&d_serialized, NUM_ELEMENTS * 12 * 48);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_serialized, nullptr) << "Failed to allocate device memory for serialized";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector serialization
    batched_error_code serialize_result = bendian_from_fp12_vec(d_serialized, d_a, NUM_ELEMENTS);
    expect_success(serialize_result, "bendian_from_fp12_vec");

    if (serialize_result == BATCHED_SUCCESS) {
        // Vector deserialization
        batched_error_code deserialize_result = fp12_from_bendian_vec(d_result, d_serialized, NUM_ELEMENTS);
        expect_success(deserialize_result, "fp12_from_bendian_vec");

        if (deserialize_result == BATCHED_SUCCESS) {
            expect_success(copy_device_to_host_fp12(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
            EXPECT_TRUE(compare_fp12_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "serialization round-trip"));
        }
    }

    // Cleanup
    free_device_fp12_array(d_a);
    cudaFree(d_serialized);
    free_device_fp12_array(d_result);
}

TEST_F(Fp12CudaTest, VectorMixedOperations) {
    const size_t MIXED_SIZE = 32; // Reduced for Fp12 complexity
    auto h_a = std::make_unique<vec384fp12[]>(MIXED_SIZE);
    auto h_b = std::make_unique<vec384fp12[]>(MIXED_SIZE);
    auto h_result = std::make_unique<vec384fp12[]>(MIXED_SIZE);
    auto cpu_result = std::make_unique<vec384fp12[]>(MIXED_SIZE);

    // Initialize test data
    for (size_t i = 0; i < MIXED_SIZE; i++) {
        set_fp12_simple(h_a[i], 3 + i);
        set_fp12_simple(h_b[i], 5 + i * 2);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(MIXED_SIZE);
    vec384fp12* d_b = allocate_device_fp12_array(MIXED_SIZE);
    vec384fp12* d_result = allocate_device_fp12_array(MIXED_SIZE);
    vec384fp12* d_temp = allocate_device_fp12_array(MIXED_SIZE);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";
    ASSERT_NE(d_temp, nullptr) << "Failed to allocate device memory for temp";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), MIXED_SIZE), "copy h_a to device");
    expect_success(copy_host_to_device_fp12(d_b, h_b.get(), MIXED_SIZE), "copy h_b to device");

    // Perform a sequence of operations: (a * b)^2
    // Step 1: Multiply
    batched_error_code mul_result = mul_fp12_vec(d_temp, d_a, d_b, MIXED_SIZE);
    expect_success(mul_result, "mixed operations mul");

    // Step 2: Square
    if (mul_result == BATCHED_SUCCESS) {
        batched_error_code sqr_result = sqr_fp12_vec(d_result, d_temp, MIXED_SIZE);
        expect_success(sqr_result, "mixed operations sqr");

        if (sqr_result == BATCHED_SUCCESS) {
            expect_success(copy_device_to_host_fp12(h_result.get(), d_result, MIXED_SIZE), "copy mixed result to host");

            // Calculate expected results on host
            for (size_t i = 0; i < MIXED_SIZE; i++) {
                vec384fp12 temp_mul;
                blst_fp12_mul(temp_mul, h_a[i], h_b[i]);
                blst_fp12_sqr(cpu_result[i], temp_mul);
            }

            EXPECT_TRUE(compare_fp12_vectors(h_result.get(), cpu_result.get(), MIXED_SIZE, "mixed operations"));
        }
    }

    // Cleanup
    free_device_fp12_array(d_a);
    free_device_fp12_array(d_b);
    free_device_fp12_array(d_result);
    free_device_fp12_array(d_temp);
}

TEST_F(Fp12CudaTest, VectorErrorHandling) {
    const size_t NUM_ELEMENTS = 100;

    // Test null pointer handling - these should fail immediately without device allocation
    batched_error_code null_out_result = mul_fp12_vec(nullptr, nullptr, nullptr, NUM_ELEMENTS);
    EXPECT_EQ(null_out_result, BATCHED_ERROR_NULL_POINTER)
        << "Expected NULL_POINTER error for null output array";

    // Test zero size handling - zero size should be OK
    vec384fp12* d_dummy = allocate_device_fp12_array(NUM_ELEMENTS);
    ASSERT_NE(d_dummy, nullptr) << "Failed to allocate device memory for dummy array";

    batched_error_code zero_size_result = mul_fp12_vec(d_dummy, d_dummy, d_dummy, 0);
    EXPECT_EQ(zero_size_result, BATCHED_SUCCESS)
        << "Expected SUCCESS for zero element count (zero-size operations should be OK)";

    free_device_fp12_array(d_dummy);
}

TEST_F(Fp12CudaTest, VectorIdentityProperties) {
    const size_t NUM_ELEMENTS = 20;
    auto h_a = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_one = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_fp12_simple(h_a[i], 3 + i);
        set_fp12_one(h_one[i]);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_one = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_result = allocate_device_fp12_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_one, nullptr) << "Failed to allocate device memory for array one";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp12(d_one, h_one.get(), NUM_ELEMENTS), "copy h_one to device");

    // Test a * 1 = a
    batched_error_code mul_result = mul_fp12_vec(d_result, d_a, d_one, NUM_ELEMENTS);
    expect_success(mul_result, "identity multiplication");

    if (mul_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp12(h_result.get(), d_result, NUM_ELEMENTS), "copy mul result to host");
        EXPECT_TRUE(compare_fp12_vectors(h_result.get(), h_a.get(), NUM_ELEMENTS, "a * 1 = a"));
    }

    // Cleanup
    free_device_fp12_array(d_a);
    free_device_fp12_array(d_one);
    free_device_fp12_array(d_result);
}

TEST_F(Fp12CudaTest, VectorPerformanceBaseline) {
    const size_t NUM_ELEMENTS = 1000; // Reasonable size for Fp12
    auto h_a = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<vec384fp12[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        random_fp12(h_a[i]);
        random_fp12(h_b[i]);
    }

    // This test primarily validates correctness at scale
    // Performance measurements would require timing infrastructure

    // CPU computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        blst_fp12_mul(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    vec384fp12* d_a = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_b = allocate_device_fp12_array(NUM_ELEMENTS);
    vec384fp12* d_result = allocate_device_fp12_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp12(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp12(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = mul_fp12_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "performance baseline");

    // Verify correctness only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp12(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fp12_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "performance baseline"));
    }

    // Cleanup
    free_device_fp12_array(d_a);
    free_device_fp12_array(d_b);
    free_device_fp12_array(d_result);
}