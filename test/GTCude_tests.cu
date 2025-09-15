#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include "config.h"

// CUDA kernel wrapper macros for device execution
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while(0)

// CUDA kernels for GT tests
GLOBAL_FUNC void cuda_test_gt_identity_kernel(bool* results);
GLOBAL_FUNC void cuda_test_gt_multiplication_kernel(bool* results);
GLOBAL_FUNC void cuda_test_gt_inversion_kernel(bool* results);
GLOBAL_FUNC void cuda_test_gt_squaring_kernel(bool* results);
GLOBAL_FUNC void cuda_test_gt_exponentiation_kernel(bool* results);
GLOBAL_FUNC void cuda_test_gt_bilinearity_kernel(bool* results);
GLOBAL_FUNC void cuda_test_gt_group_membership_kernel(bool* results);
GLOBAL_FUNC void cuda_test_gt_conjugate_kernel(bool* results);
GLOBAL_FUNC void cuda_test_gt_serialization_kernel(bool* results);
GLOBAL_FUNC void cuda_test_gt_multi_pairing_kernel(bool* results);
GLOBAL_FUNC void cuda_test_gt_final_verify_kernel(bool* results);
GLOBAL_FUNC void cuda_test_gt_powers_of_two_kernel(bool* results);

// Batch processing tests for GT
GLOBAL_FUNC void cuda_test_gt_batch_operations_kernel(bool* results, int batch_size);
GLOBAL_FUNC void cuda_test_gt_performance_kernel(int* operations_completed);

// Helper function to run a CUDA kernel and check results
template<typename KernelFunc>
bool run_cuda_gt_test(KernelFunc kernel, int num_results = 1) {
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
            printf("GT Test failed at index %d\n", i);
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
class GTCUDAEnvironment : public ::testing::Environment {
public:
    void SetUp() override {
        int device_count = 0;
        cudaError_t error = cudaGetDeviceCount(&device_count);

        if (error != cudaSuccess || device_count == 0) {
            std::cout << "No CUDA devices available. CUDA GT tests will be skipped.\n";
            cuda_device_available = false;
        }
        else {
            std::cout << "Found " << device_count << " CUDA device(s) for GT tests\n";
            cuda_device_available = true;
        }
    }

    static bool IsCudaAvailable() { return cuda_device_available; }
};

// Initialize the test environment
struct GTCudaTestRegistration {
    GTCudaTestRegistration() {
        ::testing::AddGlobalTestEnvironment(new GTCUDAEnvironment());
    }
};

static GTCudaTestRegistration gt_cuda_test_registration;

// CUDA test cases for GT operations
TEST(GTOpsCUDA, Identity) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_identity_kernel, 4));
}

TEST(GTOpsCUDA, Multiplication) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_multiplication_kernel, 4));
}

TEST(GTOpsCUDA, Inversion) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_inversion_kernel, 3));
}

TEST(GTOpsCUDA, Squaring) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_squaring_kernel, 3));
}

TEST(GTOpsCUDA, Exponentiation) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_exponentiation_kernel, 4));
}

TEST(GTOpsCUDA, Bilinearity) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_bilinearity_kernel, 3));
}

TEST(GTOpsCUDA, GroupMembership) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_group_membership_kernel, 4));
}

TEST(GTOpsCUDA, Conjugate) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_conjugate_kernel, 3));
}

TEST(GTOpsCUDA, Serialization) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_serialization_kernel, 2));
}

TEST(GTOpsCUDA, MultiPairing) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_multi_pairing_kernel, 2));
}

TEST(GTOpsCUDA, FinalVerify) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_final_verify_kernel, 3));
}

TEST(GTOpsCUDA, PowersOfTwo) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_gt_test(cuda_test_gt_powers_of_two_kernel, 1));
}

// Performance and batch processing tests
TEST(GTOpsCUDA, PerformanceBaseline) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
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
    cuda_test_gt_performance_kernel << <1, 1 >> > (d_ops);
    CUDA_CHECK(cudaEventRecord(stop));

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));

    CUDA_CHECK(cudaMemcpy(&h_ops, d_ops, sizeof(int), cudaMemcpyDeviceToHost));

    printf("CUDA GT Performance: %d operations in %.3f ms (%.1f ops/ms)\n",
        h_ops, milliseconds, h_ops / milliseconds);

    CUDA_CHECK(cudaFree(d_ops));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));

    EXPECT_GT(h_ops, 0);
}

TEST(GTOpsCUDA, BatchProcessing) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const int batch_size = 128;  // Reduced for GT complexity (most expensive operations)
    bool* d_results, * h_results;

    // Allocate memory
    h_results = new bool[batch_size];
    CUDA_CHECK(cudaMalloc(&d_results, batch_size * sizeof(bool)));
    CUDA_CHECK(cudaMemset(d_results, 0, batch_size * sizeof(bool)));

    // Launch kernel with multiple threads
    int threads_per_block = 32;  // Reduced for GT complexity
    int blocks = (batch_size + threads_per_block - 1) / threads_per_block;

    cuda_test_gt_batch_operations_kernel << <blocks, threads_per_block >> > (d_results, batch_size);
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

    printf("GT Batch processing: %d/%d operations successful\n", successful_operations, batch_size);

    // Cleanup
    delete[] h_results;
    CUDA_CHECK(cudaFree(d_results));

    // Expect at least 95% success rate
    EXPECT_GE(successful_operations, batch_size * 0.95);
}

// Additional GT-specific tests that require more complex setups
TEST(GTOpsCUDA, AdvancedExponentiation) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    bool* d_results, * h_results;
    const int num_results = 3;

    h_results = new bool[num_results];
    CUDA_CHECK(cudaMalloc(&d_results, num_results * sizeof(bool)));
    CUDA_CHECK(cudaMemset(d_results, 0, num_results * sizeof(bool)));

    // Custom kernel for advanced exponentiation tests
    auto advanced_exp_test = [d_results] __global__() {
        vec384fp12 gt_element, result, expected;

        // Create GT element
        POINTonE1 g1_point = *blst_p1_generator();
        POINTonE2 g2_point = *blst_p2_generator();
        POINTonE1_affine g1_affine;
        POINTonE2_affine g2_affine;
        blst_p1_to_affine(&g1_affine, &g1_point);
        blst_p2_to_affine(&g2_affine, &g2_point);
        blst_miller_loop(gt_element, &g2_affine, &g1_affine);
        blst_final_exp(gt_element, gt_element);

        // Test power of 10
        limb_t ten[4] = { 10, 0, 0, 0 };
        blst_fp12_pow(result, gt_element, ten, sizeof(ten) * 8);

        // Calculate e^10 by repeated multiplication
        blst_fp12_one_copy(expected);
        for (int i = 0; i < 10; i++) {
            blst_fp12_mul(expected, expected, gt_element);
        }
        d_results[0] = blst_fp12_is_equal(result, expected);

        // Test with scattered bits
        limb_t scattered_bits[4] = { 0, 0, 0, 0 };
        scattered_bits[0] = 1 | (1 << 4) | (1 << 8);  // Bits 0, 4, and 8 set
        blst_fp12_pow(result, gt_element, scattered_bits, sizeof(scattered_bits) * 8);

        // Calculate manually
        blst_fp12_one_copy(expected);
        blst_fp12_mul(expected, expected, gt_element);  // Bit 0
        for (int i = 0; i < 4; i++) {
            blst_fp12_sqr(expected, expected);
        }
        blst_fp12_mul(expected, expected, gt_element);  // Bit 4
        for (int i = 0; i < 4; i++) {
            blst_fp12_sqr(expected, expected);
        }
        blst_fp12_mul(expected, expected, gt_element);  // Bit 8

        d_results[1] = blst_fp12_is_equal(result, expected);

        // Test that result is in group
        d_results[2] = blst_fp12_in_group(result);
    };

    advanced_exp_test << <1, 1 >> > ();
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_results, d_results, num_results * sizeof(bool), cudaMemcpyDeviceToHost));

    for (int i = 0; i < num_results; i++) {
        EXPECT_TRUE(h_results[i]) << "Advanced exponentiation test failed at index " << i;
    }

    delete[] h_results;
    CUDA_CHECK(cudaFree(d_results));
}

TEST(GTOpsCUDA, PairingWithPrecomputation) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    bool* d_results, * h_results;
    const int num_results = 2;

    h_results = new bool[num_results];
    CUDA_CHECK(cudaMalloc(&d_results, num_results * sizeof(bool)));
    CUDA_CHECK(cudaMemset(d_results, 0, num_results * sizeof(bool)));

    // Custom kernel for precomputation tests
    auto precomp_test = [d_results] __global__() {
        vec384fp6 lines[68];  // Size from the API
        vec384fp12 pre_result, direct_result;
        POINTonE1 g1_point = *blst_p1_generator();
        POINTonE2 g2_point = *blst_p2_generator();
        POINTonE1_affine g1_affine;
        POINTonE2_affine g2_affine;

        blst_p1_to_affine(&g1_affine, &g1_point);
        blst_p2_to_affine(&g2_affine, &g2_point);

        // Precompute lines for pairing with g2_point
        blst_precompute_lines(lines, &g2_affine);

        // Compute pairing using precomputed lines
        blst_miller_loop_lines(pre_result, lines, &g1_affine);
        blst_final_exp(pre_result, pre_result);

        // Compute direct pairing
        blst_miller_loop(direct_result, &g2_affine, &g1_affine);
        blst_final_exp(direct_result, direct_result);

        d_results[0] = blst_fp12_is_equal(pre_result, direct_result);

        // Test that both results are in group
        d_results[1] = blst_fp12_in_group(pre_result) && blst_fp12_in_group(direct_result);
    };

    precomp_test << <1, 1 >> > ();
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_results, d_results, num_results * sizeof(bool), cudaMemcpyDeviceToHost));

    for (int i = 0; i < num_results; i++) {
        EXPECT_TRUE(h_results[i]) << "Precomputation test failed at index " << i;
    }

    delete[] h_results;
    CUDA_CHECK(cudaFree(d_results));
}

TEST(GTOpsCUDA, AssociativityProperty) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    bool* d_results, * h_results;
    h_results = new bool[1];
    CUDA_CHECK(cudaMalloc(&d_results, sizeof(bool)));
    CUDA_CHECK(cudaMemset(d_results, 0, sizeof(bool)));

    // Custom kernel for associativity test
    auto associativity_test = [d_results] __global__() {
        vec384fp12 a, b, c, temp1, temp2;
        POINTonE1 g1_point = *blst_p1_generator();
        POINTonE2 g2_point = *blst_p2_generator();
        POINTonE1_affine g1_affine;
        POINTonE2_affine g2_affine;

        blst_p1_to_affine(&g1_affine, &g1_point);
        blst_p2_to_affine(&g2_affine, &g2_point);

        // Create three GT elements
        blst_miller_loop(a, &g2_affine, &g1_affine);
        blst_final_exp(a, a);

        // b = a^2
        blst_fp12_sqr(b, a);

        // c = a^3
        blst_fp12_mul(c, a, b);

        // Test associativity: (a * b) * c = a * (b * c)
        blst_fp12_mul(temp1, a, b);      // a * b
        blst_fp12_mul(temp1, temp1, c);  // (a * b) * c

        blst_fp12_mul(temp2, b, c);      // b * c
        blst_fp12_mul(temp2, a, temp2);  // a * (b * c)

        d_results[0] = blst_fp12_is_equal(temp1, temp2);
    };

    associativity_test << <1, 1 >> > ();
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_results, d_results, sizeof(bool), cudaMemcpyDeviceToHost));
    EXPECT_TRUE(h_results[0]) << "Associativity property failed";

    delete[] h_results;
    CUDA_CHECK(cudaFree(d_results));
}