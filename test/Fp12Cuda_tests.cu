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

// CUDA kernels for Fp12 tests

GLOBAL_FUNC void cuda_test_fp12_zero_and_one_kernel(bool* results);
GLOBAL_FUNC void cuda_test_fp12_multiplication_kernel(bool* results);
GLOBAL_FUNC void cuda_test_fp12_squaring_kernel(bool* results);
GLOBAL_FUNC void cuda_test_fp12_conjugation_kernel(bool* results);
GLOBAL_FUNC void cuda_test_fp12_inversion_kernel(bool* results);
GLOBAL_FUNC void cuda_test_fp12_equality_kernel(bool* results);
GLOBAL_FUNC void cuda_test_fp12_serialization_kernel(bool* results);

// Batch processing tests for Fp12
GLOBAL_FUNC void cuda_test_fp12_batch_operations_kernel(bool* results, int batch_size);
GLOBAL_FUNC void cuda_test_fp12_performance_kernel(int* operations_completed);



// Helper function to run a CUDA kernel and check results
template<typename KernelFunc>
bool run_cuda_fp12_test(KernelFunc kernel, int num_results = 1) {
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
            std::cout << "No CUDA devices available. CUDA Fp12 tests will be skipped.\n";
            cuda_device_available = false;
        }
        else {
            std::cout << "Found " << device_count << " CUDA device(s) for Fp12 tests\n";
            cuda_device_available = true;
        }
    }

    static bool IsCudaAvailable() { return cuda_device_available; }
};

// Initialize the test environment
struct Fp12CudaTestRegistration {
    Fp12CudaTestRegistration() {
        ::testing::AddGlobalTestEnvironment(new CUDAEnvironment());
    }
};

static Fp12CudaTestRegistration fp12_cuda_test_registration;

// CUDA test cases for Fp12 operations
TEST(Fp12OpsCUDA, ZeroAndOne) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp12_test(cuda_test_fp12_zero_and_one_kernel, 6));
}

TEST(Fp12OpsCUDA, Multiplication) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp12_test(cuda_test_fp12_multiplication_kernel, 4));
}

TEST(Fp12OpsCUDA, Squaring) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp12_test(cuda_test_fp12_squaring_kernel, 2));
}

TEST(Fp12OpsCUDA, Conjugation) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp12_test(cuda_test_fp12_conjugation_kernel, 2));
}

TEST(Fp12OpsCUDA, Inversion) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp12_test(cuda_test_fp12_inversion_kernel, 2));
}

TEST(Fp12OpsCUDA, Equality) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp12_test(cuda_test_fp12_equality_kernel, 4));
}

TEST(Fp12OpsCUDA, Serialization) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }
    EXPECT_TRUE(run_cuda_fp12_test(cuda_test_fp12_serialization_kernel, 2));
}

// Performance and batch processing tests
TEST(Fp12OpsCUDA, PerformanceBaseline) {
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
    cuda_test_fp12_performance_kernel << <1, 1 >> > (d_ops);
    CUDA_CHECK(cudaEventRecord(stop));

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));

    CUDA_CHECK(cudaMemcpy(&h_ops, d_ops, sizeof(int), cudaMemcpyDeviceToHost));

    printf("CUDA Fp12 Performance: %d operations in %.3f ms (%.1f ops/ms)\n",
        h_ops, milliseconds, h_ops / milliseconds);

    CUDA_CHECK(cudaFree(d_ops));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));

    EXPECT_GT(h_ops, 0);
}

TEST(Fp12OpsCUDA, BatchProcessing) {
    if (!CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const int batch_size = 256;  // Reduced for Fp12 complexity
    bool* d_results, * h_results;

    // Allocate memory
    h_results = new bool[batch_size];
    CUDA_CHECK(cudaMalloc(&d_results, batch_size * sizeof(bool)));
    CUDA_CHECK(cudaMemset(d_results, 0, batch_size * sizeof(bool)));

    // Launch kernel with multiple threads
    int threads_per_block = 64;  // Reduced for Fp12 complexity
    int blocks = (batch_size + threads_per_block - 1) / threads_per_block;

    cuda_test_fp12_batch_operations_kernel << <blocks, threads_per_block >> > (d_results, batch_size);
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

    printf("Fp12 Batch processing: %d/%d operations successful\n", successful_operations, batch_size);

    // Cleanup
    delete[] h_results;
    CUDA_CHECK(cudaFree(d_results));

    // Expect at least 95% success rate
    EXPECT_GE(successful_operations, batch_size * 0.95);
}