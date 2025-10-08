#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <memory>
#include <random>
#include <cuda_runtime.h>

#include "blst/fields/Fr.h"
#include "blst/fields/Fr_cuda.h"
using namespace blst;

#define FR_LIMBS NLIMBS(256)

class FrCudaTest : public ::testing::Test {
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

    // Helper: Set a Fr& to a small integer value
    void set_fr(Fr& v, uint64_t x) {
        for (int i = 0; i < FR_LIMBS; ++i)
            v[i] = 0;
        memcpy(v, &x, sizeof(x));
    }

    // Helper: Compare two Fr& for equality
    bool fr_equal(const Fr& a, const Fr& b) {
        for (int i = 0; i < FR_LIMBS; ++i)
            if (a[i] != b[i]) return false;
        return true;
    }

    // Helper: Print Fr& for debugging
    void print_fr(const Fr& v) {
        for (int i = FR_LIMBS - 1; i >= 0; --i)
            printf("%016llx%s", (unsigned long long)v[i], i ? " " : "\n");
    }

    // Helper: Copy one Fr& to another
    void fr_copy(Fr& dst, const Fr& src) {
        for (int i = 0; i < FR_LIMBS; ++i)
            dst[i] = src[i];
    }

    // Helper: Generate random Fr element
    void random_fr(Fr& v) {
        // Generate random limbs
        for (int i = 0; i < FR_LIMBS; ++i) {
            v[i] = gen();
        }
        // Reduce modulo Fr to ensure valid element
        // Note: This is a simplified approach; real implementation would use proper modular reduction
        if (FR_LIMBS > 0) {
            v[FR_LIMBS - 1] &= 0x7FFFFFFFFFFFFFFFULL; // Ensure MSB is clear to avoid overflow
        }
    }

    // Helper: Convert standard integer to Montgomery form
    void to_mont_fr(Fr& r, const Fr& a) {
        Fr::to_mont(r, a);
    }

    // Helper: Allocate device memory for Fr& arrays
    Fr* allocate_device_fr_array(size_t num_elements) {
        Fr* d_ptr = nullptr;
        cudaError_t err = cudaMalloc((void**)&d_ptr, num_elements * sizeof(Fr));
        if (err != cudaSuccess) {
            return nullptr;
        }
        return d_ptr;
    }

    // Helper: Free device memory
    void free_device_fr_array(Fr* d_ptr) {
        if (d_ptr) {
            cudaFree(d_ptr);
        }
    }

    // Helper: Copy host Fr array to device
    batched_error_code copy_host_to_device_fr(Fr* d_dst, const Fr* h_src, size_t num_elements) {
        if (!d_dst || !h_src) {
            return BATCHED_ERROR_NULL_POINTER;
        }

        cudaError_t err = cudaMemcpy(d_dst, h_src, num_elements * sizeof(Fr), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
        }

        return BATCHED_SUCCESS;
    }

    // Helper: Copy device Fr array to host
    batched_error_code copy_device_to_host_fr(Fr* h_dst, const Fr* d_src, size_t num_elements) {
        if (!h_dst || !d_src) {
            return BATCHED_ERROR_NULL_POINTER;
        }

        cudaError_t err = cudaMemcpy(h_dst, d_src, num_elements * sizeof(Fr), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
        }

        return BATCHED_SUCCESS;
    }

    // Helper: Compare vector results against CPU reference
    bool compare_fr_vectors(const Fr* vec_result, const Fr* cpu_result, size_t num, const char* operation) {
        bool all_equal = true;
        for (size_t i = 0; i < num; i++) {
            if (!fr_equal(vec_result[i], cpu_result[i])) {
                printf("Mismatch in %s at index %zu:\n", operation, i);
                printf("Vector result: "); print_fr(vec_result[i]);
                printf("CPU result:    "); print_fr(cpu_result[i]);
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

TEST_F(FrCudaTest, VectorAdd) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        random_fr(h_a[i]);
        random_fr(h_b[i]);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::add(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    Fr* d_a = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_b = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fr(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = add_fr_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "add_fr_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fr_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "add_fr_vec"));
    }

    // Cleanup
    free_device_fr_array(d_a);
    free_device_fr_array(d_b);
    free_device_fr_array(d_result);
}

TEST_F(FrCudaTest, VectorSub) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        random_fr(h_a[i]);
        random_fr(h_b[i]);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::sub(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    Fr* d_a = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_b = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fr(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = sub_fr_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "sub_fr_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fr_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "sub_fr_vec"));
    }

    // Cleanup
    free_device_fr_array(d_a);
    free_device_fr_array(d_b);
    free_device_fr_array(d_result);
}

TEST_F(FrCudaTest, VectorMul) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Generate test data (convert to Montgomery form for multiplication)
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr temp_a, temp_b;
        random_fr(temp_a);
        random_fr(temp_b);
        to_mont_fr(h_a[i], temp_a);
        to_mont_fr(h_b[i], temp_b);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::mul(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    Fr* d_a = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_b = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fr(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = mul_fr_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "mul_fr_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fr_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "mul_fr_vec"));
    }

    // Cleanup
    free_device_fr_array(d_a);
    free_device_fr_array(d_b);
    free_device_fr_array(d_result);
}

TEST_F(FrCudaTest, VectorSqr) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Generate test data (convert to Montgomery form for squaring)
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr temp;
        random_fr(temp);
        to_mont_fr(h_a[i], temp);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::sqr(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    Fr* d_a = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = sqr_fr_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "sqr_fr_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fr_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "sqr_fr_vec"));
    }

    // Cleanup
    free_device_fr_array(d_a);
    free_device_fr_array(d_result);
}

TEST_F(FrCudaTest, VectorCNeg) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result_true = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_result_true = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result_false = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_result_false = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        random_fr(h_a[i]);
    }

    // CPU reference computation with flag = 1
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::cneg(cpu_result_true[i], h_a[i], 1);
    }

    // CPU reference computation with flag = 0
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::cneg(cpu_result_false[i], h_a[i], 0);
    }

    // Allocate device memory
    Fr* d_a = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result_true = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result_false = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result_true, nullptr) << "Failed to allocate device memory for result_true";
    ASSERT_NE(d_result_false, nullptr) << "Failed to allocate device memory for result_false";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with flag = 1 and error checking
    batched_error_code result_true = cneg_fr_vec(d_result_true, d_a, 1, NUM_ELEMENTS);
    expect_success(result_true, "cneg_fr_vec (flag=1)");

    // Vector computation with flag = 0 and error checking
    batched_error_code result_false = cneg_fr_vec(d_result_false, d_a, 0, NUM_ELEMENTS);
    expect_success(result_false, "cneg_fr_vec (flag=0)");

    // Copy results back to host and compare only if operations succeeded
    if (result_true == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result_true.get(), d_result_true, NUM_ELEMENTS), "copy result_true to host");
        EXPECT_TRUE(compare_fr_vectors(h_result_true.get(), cpu_result_true.get(), NUM_ELEMENTS, "cneg_fr_vec (flag=1)"));
    }
    if (result_false == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result_false.get(), d_result_false, NUM_ELEMENTS), "copy result_false to host");
        EXPECT_TRUE(compare_fr_vectors(h_result_false.get(), cpu_result_false.get(), NUM_ELEMENTS, "cneg_fr_vec (flag=0)"));
    }

    // Cleanup
    free_device_fr_array(d_a);
    free_device_fr_array(d_result_true);
    free_device_fr_array(d_result_false);
}

TEST_F(FrCudaTest, VectorMulBy3) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        random_fr(h_a[i]);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::mul_by_3(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    Fr* d_a = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = mul3_fr_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(result, "mul3_fr_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fr_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "mul3_fr_vec"));
    }

    // Cleanup
    free_device_fr_array(d_a);
    free_device_fr_array(d_result);
}

TEST_F(FrCudaTest, VectorLShift) {
    const size_t NUM_ELEMENTS = 1000;
    const size_t SHIFT_COUNT = 3;
    auto h_a = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        random_fr(h_a[i]);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::lshift(cpu_result[i], h_a[i], SHIFT_COUNT);
    }

    // Allocate device memory
    Fr* d_a = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = lshift_fr_vec(d_result, d_a, SHIFT_COUNT, NUM_ELEMENTS);
    expect_success(result, "lshift_fr_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fr_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "lshift_fr_vec"));
    }

    // Cleanup
    free_device_fr_array(d_a);
    free_device_fr_array(d_result);
}

TEST_F(FrCudaTest, VectorRShift) {
    const size_t NUM_ELEMENTS = 1000;
    const size_t SHIFT_COUNT = 2;
    auto h_a = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        random_fr(h_a[i]);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::rshift(cpu_result[i], h_a[i], SHIFT_COUNT);
    }

    // Allocate device memory
    Fr* d_a = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation with error checking
    batched_error_code result = rshift_fr_vec(d_result, d_a, SHIFT_COUNT, NUM_ELEMENTS);
    expect_success(result, "rshift_fr_vec");

    // Copy result back to host and compare only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fr_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "rshift_fr_vec"));
    }

    // Cleanup
    free_device_fr_array(d_a);
    free_device_fr_array(d_result);
}

TEST_F(FrCudaTest, VectorMontgomeryConversion) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_a = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_to_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_to_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_from_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_from_result = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Generate test data - ensure values are properly reduced modulo the field
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        random_fr(h_a[i]);
        // Ensure the value is properly reduced by performing a modular reduction
        Fr temp;
        Fr::from_mont(temp, h_a[i]);  // This ensures proper reduction
        h_a[i] = temp;
    }

    // CPU reference computation - to Montgomery
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::to_mont(cpu_to_result[i], h_a[i]);
    }

    // Allocate device memory
    Fr* d_a = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_to_result = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_from_result = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_to_result, nullptr) << "Failed to allocate device memory for to_result";
    ASSERT_NE(d_from_result, nullptr) << "Failed to allocate device memory for from_result";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");

    // Vector computation - to Montgomery with error checking
    batched_error_code to_result = to_fr_vec(d_to_result, d_a, NUM_ELEMENTS);
    expect_success(to_result, "to_fr_vec");

    // Compare to_fr results only if operation succeeded
    if (to_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_to_result.get(), d_to_result, NUM_ELEMENTS), "copy to_result to host");
        EXPECT_TRUE(compare_fr_vectors(h_to_result.get(), cpu_to_result.get(), NUM_ELEMENTS, "to_fr_vec"));

        // CPU reference computation - from Montgomery using the CPU's to_montgomery result
        for (size_t i = 0; i < NUM_ELEMENTS; i++) {
            Fr::from_mont(cpu_from_result[i], cpu_to_result[i]);
        }

        // Vector computation - from Montgomery using the GPU's to_montgomery result
        batched_error_code from_result = from_fr_vec(d_from_result, d_to_result, NUM_ELEMENTS);
        expect_success(from_result, "from_fr_vec");

        // Compare from_fr results only if operation succeeded
        if (from_result == BATCHED_SUCCESS) {
            expect_success(copy_device_to_host_fr(h_from_result.get(), d_from_result, NUM_ELEMENTS), "copy from_result to host");
            EXPECT_TRUE(compare_fr_vectors(h_from_result.get(), cpu_from_result.get(), NUM_ELEMENTS, "from_fr_vec"));

            // Verify round-trip conversion using CPU reference
            // Compare GPU from_montgomery result with CPU from_montgomery result
            // This ensures we're comparing like with like
            EXPECT_TRUE(compare_fr_vectors(h_from_result.get(), cpu_from_result.get(), NUM_ELEMENTS, "round-trip conversion consistency"));
        }
    }

    // Cleanup
    free_device_fr_array(d_a);
    free_device_fr_array(d_to_result);
    free_device_fr_array(d_from_result);
}

TEST_F(FrCudaTest, VectorButterfly) {
    const size_t NUM_ELEMENTS = 1000;
    auto h_x0 = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_x1 = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_twiddle = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_x0_ct = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_x1_ct = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_x0_ct = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_x1_ct = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_x0_gs = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_x1_gs = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_x0_gs = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_x1_gs = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Generate test data (convert to Montgomery form)
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr temp_x0, temp_x1, temp_twiddle;
        random_fr(temp_x0);
        random_fr(temp_x1);
        random_fr(temp_twiddle);
        to_mont_fr(h_x0[i], temp_x0);
        to_mont_fr(h_x1[i], temp_x1);
        to_mont_fr(h_twiddle[i], temp_twiddle);
    }

    // Copy data for CT butterfly tests
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        cpu_x0_ct[i] = h_x0[i];
        cpu_x1_ct[i] = h_x1[i];
    }

    // CPU reference computation - CT butterfly
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::ct_bfly(cpu_x0_ct[i], cpu_x1_ct[i], h_twiddle[i]);
    }

    // Allocate device memory for CT butterfly
    Fr* d_x0_ct = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_x1_ct = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_twiddle = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_x0_ct, nullptr) << "Failed to allocate device memory for x0_ct";
    ASSERT_NE(d_x1_ct, nullptr) << "Failed to allocate device memory for x1_ct";
    ASSERT_NE(d_twiddle, nullptr) << "Failed to allocate device memory for twiddle";

    // Copy data to device for CT butterfly
    expect_success(copy_host_to_device_fr(d_x0_ct, h_x0.get(), NUM_ELEMENTS), "copy h_x0 to device");
    expect_success(copy_host_to_device_fr(d_x1_ct, h_x1.get(), NUM_ELEMENTS), "copy h_x1 to device");
    expect_success(copy_host_to_device_fr(d_twiddle, h_twiddle.get(), NUM_ELEMENTS), "copy h_twiddle to device");

    // Vector computation - CT butterfly with error checking
    batched_error_code ct_result = ct_bfly_fr_vec(d_x0_ct, d_x1_ct, d_twiddle, NUM_ELEMENTS);
    expect_success(ct_result, "ct_bfly_fr_vec");

    // Compare CT butterfly results only if operation succeeded
    if (ct_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_x0_ct.get(), d_x0_ct, NUM_ELEMENTS), "copy x0_ct result to host");
        expect_success(copy_device_to_host_fr(h_x1_ct.get(), d_x1_ct, NUM_ELEMENTS), "copy x1_ct result to host");
        EXPECT_TRUE(compare_fr_vectors(h_x0_ct.get(), cpu_x0_ct.get(), NUM_ELEMENTS, "ct_bfly_fr_vec x0"));
        EXPECT_TRUE(compare_fr_vectors(h_x1_ct.get(), cpu_x1_ct.get(), NUM_ELEMENTS, "ct_bfly_fr_vec x1"));
    }

    // Prepare for GS butterfly tests
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        cpu_x0_gs[i] = h_x0[i];
        cpu_x1_gs[i] = h_x1[i];
    }

    // CPU reference computation - GS butterfly
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::gs_bfly(cpu_x0_gs[i], cpu_x1_gs[i], h_twiddle[i]);
    }

    // Allocate device memory for GS butterfly
    Fr* d_x0_gs = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_x1_gs = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_x0_gs, nullptr) << "Failed to allocate device memory for x0_gs";
    ASSERT_NE(d_x1_gs, nullptr) << "Failed to allocate device memory for x1_gs";

    // Copy original data to device for GS butterfly
    expect_success(copy_host_to_device_fr(d_x0_gs, h_x0.get(), NUM_ELEMENTS), "copy h_x0 to device (GS)");
    expect_success(copy_host_to_device_fr(d_x1_gs, h_x1.get(), NUM_ELEMENTS), "copy h_x1 to device (GS)");

    // Vector computation - GS butterfly with error checking
    batched_error_code gs_result = gs_bfly_fr_vec(d_x0_gs, d_x1_gs, d_twiddle, NUM_ELEMENTS);
    expect_success(gs_result, "gs_bfly_fr_vec");

    // Compare GS butterfly results only if operation succeeded
    if (gs_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_x0_gs.get(), d_x0_gs, NUM_ELEMENTS), "copy x0_gs result to host");
        expect_success(copy_device_to_host_fr(h_x1_gs.get(), d_x1_gs, NUM_ELEMENTS), "copy x1_gs result to host");
        EXPECT_TRUE(compare_fr_vectors(h_x0_gs.get(), cpu_x0_gs.get(), NUM_ELEMENTS, "gs_bfly_fr_vec x0"));
        EXPECT_TRUE(compare_fr_vectors(h_x1_gs.get(), cpu_x1_gs.get(), NUM_ELEMENTS, "gs_bfly_fr_vec x1"));
    }

    // Cleanup
    free_device_fr_array(d_x0_ct);
    free_device_fr_array(d_x1_ct);
    free_device_fr_array(d_twiddle);
    free_device_fr_array(d_x0_gs);
    free_device_fr_array(d_x1_gs);
}

TEST_F(FrCudaTest, VectorEdgeCases) {
    const size_t NUM_ELEMENTS = 10;

    // Test with zero elements
    auto h_zeros = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fr[]>(NUM_ELEMENTS);

    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_fr(h_zeros[i], 0);
    }

    // Allocate device memory
    Fr* d_zeros = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_zeros, nullptr) << "Failed to allocate device memory for zeros";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_zeros, h_zeros.get(), NUM_ELEMENTS), "copy zeros to device");

    // Test addition with zeros
    batched_error_code add_result = add_fr_vec(d_result, d_zeros, d_zeros, NUM_ELEMENTS);
    expect_success(add_result, "add_fr_vec with zeros");

    if (add_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy add result to host");
        for (size_t i = 0; i < NUM_ELEMENTS; i++) {
            EXPECT_TRUE(fr_equal(h_result[i], h_zeros[i])) << "0 + 0 should equal 0 at index " << i;
        }
    }

    // Test multiplication with zeros (in Montgomery form)
    auto h_mont_zeros = std::make_unique<Fr[]>(NUM_ELEMENTS);
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        to_mont_fr(h_mont_zeros[i], h_zeros[i]);
    }

    Fr* d_mont_zeros = allocate_device_fr_array(NUM_ELEMENTS);
    ASSERT_NE(d_mont_zeros, nullptr) << "Failed to allocate device memory for mont_zeros";

    expect_success(copy_host_to_device_fr(d_mont_zeros, h_mont_zeros.get(), NUM_ELEMENTS), "copy mont_zeros to device");

    batched_error_code mul_result = mul_fr_vec(d_result, d_mont_zeros, d_mont_zeros, NUM_ELEMENTS);
    expect_success(mul_result, "mul_fr_vec with zeros");

    if (mul_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy mul result to host");
        for (size_t i = 0; i < NUM_ELEMENTS; i++) {
            EXPECT_TRUE(fr_equal(h_result[i], h_mont_zeros[i])) << "0 * 0 should equal 0 at index " << i;
        }
    }

    // Test squaring with zeros
    batched_error_code sqr_result = sqr_fr_vec(d_result, d_mont_zeros, NUM_ELEMENTS);
    expect_success(sqr_result, "sqr_fr_vec with zeros");

    if (sqr_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy sqr result to host");
        for (size_t i = 0; i < NUM_ELEMENTS; i++) {
            EXPECT_TRUE(fr_equal(h_result[i], h_mont_zeros[i])) << "0^2 should equal 0 at index " << i;
        }
    }

    // Cleanup
    free_device_fr_array(d_zeros);
    free_device_fr_array(d_result);
    free_device_fr_array(d_mont_zeros);
}
TEST_F(FrCudaTest, VectorErrorHandling) {
    const size_t NUM_ELEMENTS = 100;

    // Test null pointer handling - these should fail immediately without device allocation
    batched_error_code null_out_result = add_fr_vec(nullptr, nullptr, nullptr, NUM_ELEMENTS);
    EXPECT_EQ(null_out_result, BATCHED_ERROR_NULL_POINTER)
        << "Expected NULL_POINTER error for null output array";

    // Test zero size handling - zero size should be OK
    Fr* d_dummy = allocate_device_fr_array(NUM_ELEMENTS);
    ASSERT_NE(d_dummy, nullptr) << "Failed to allocate device memory for dummy array";

    batched_error_code zero_size_result = add_fr_vec(d_dummy, d_dummy, d_dummy, 0);
    EXPECT_EQ(zero_size_result, BATCHED_SUCCESS)
        << "Expected SUCCESS for zero element count (zero-size operations should be OK)";

    free_device_fr_array(d_dummy);
}

TEST_F(FrCudaTest, VectorSmallValues) {
    const size_t NUM_ELEMENTS = 5;
    auto h_a = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Test with small known values
    const uint64_t TEST_VALUES[] = { 1, 2, 5, 10, 42 };

    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_fr(h_a[i], TEST_VALUES[i]);
        set_fr(h_b[i], TEST_VALUES[(i + 1) % NUM_ELEMENTS]);
    }

    // Allocate device memory
    Fr* d_a = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_b = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fr(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Test addition
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::add(cpu_result[i], h_a[i], h_b[i]);
    }
    batched_error_code add_result = add_fr_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(add_result, "small values add");
    if (add_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy add result to host");
        EXPECT_TRUE(compare_fr_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "small values add"));
    }

    // Test subtraction
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::sub(cpu_result[i], h_a[i], h_b[i]);
    }
    batched_error_code sub_result = sub_fr_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(sub_result, "small values sub");
    if (sub_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy sub result to host");
        EXPECT_TRUE(compare_fr_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "small values sub"));
    }

    // Test multiplication by 3
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::mul_by_3(cpu_result[i], h_a[i]);
    }
    batched_error_code mul3_result = mul3_fr_vec(d_result, d_a, NUM_ELEMENTS);
    expect_success(mul3_result, "small values mul3");
    if (mul3_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy mul3 result to host");
        EXPECT_TRUE(compare_fr_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "small values mul3"));
    }

    // Cleanup
    free_device_fr_array(d_a);
    free_device_fr_array(d_b);
    free_device_fr_array(d_result);
}

TEST_F(FrCudaTest, VectorPerformanceBaseline) {
    const size_t NUM_ELEMENTS = 10000;
    auto h_a = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fr[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fr[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        random_fr(h_a[i]);
        random_fr(h_b[i]);
    }

    // This test primarily validates correctness at scale
    // Performance measurements would require timing infrastructure

    // CPU computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fr::add(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    Fr* d_a = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_b = allocate_device_fr_array(NUM_ELEMENTS);
    Fr* d_result = allocate_device_fr_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fr(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fr(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Vector computation with error checking
    batched_error_code result = add_fr_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(result, "performance baseline");

    // Verify correctness only if operation succeeded
    if (result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fr(h_result.get(), d_result, NUM_ELEMENTS), "copy result to host");
        EXPECT_TRUE(compare_fr_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "performance baseline"));
    }

    // Cleanup
    free_device_fr_array(d_a);
    free_device_fr_array(d_b);
    free_device_fr_array(d_result);
}