#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>

extern "C" {
#include "src/fields.h"
}

#define LIMBS NLIMBS(384)
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, \
                    cudaGetErrorString(err)); \
            exit(1); \
        } \
    } while(0)

// CUDA kernel for testing add_fp on GPU
__global__ void test_add_fp_kernel(vec384* results, const vec384* a_vals,
    const vec384* b_vals, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        add_fp(results[idx], a_vals[idx], b_vals[idx]);
    }
}

// CUDA kernel for batch testing with different values
__global__ void test_add_fp_batch_kernel(vec384* results, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        vec384 a, b;

        // Initialize different test values for each thread
        for (int i = 0; i < LIMBS; ++i) {
            a[i] = 0;
            b[i] = 0;
        }

        // Set unique values based on thread index
        a[0] = idx + 1;      // a = 1, 2, 3, 4, ...
        b[0] = (idx + 1) * 2; // b = 2, 4, 6, 8, ...

        add_fp(results[idx], a, b);
    }
}

// Helper functions (same as CPU tests)
void set_vec384_cuda(vec384 v, uint64_t x) {
    for (int i = 0; i < LIMBS; ++i)
        v[i] = 0;
    memcpy(v, &x, sizeof(x));
}

bool vec384_equal_cuda(const vec384 a, const vec384 b) {
    for (int i = 0; i < LIMBS; ++i)
        if (a[i] != b[i]) return false;
    return true;
}

void print_vec384_cuda(const vec384 v) {
    for (int i = LIMBS - 1; i >= 0; --i)
        printf("%016llx%s", (unsigned long long)v[i], i ? " " : "\n");
}

class FpCudaOps : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize CUDA device
        int device;
        CUDA_CHECK(cudaGetDevice(&device));

        // Print device info for debugging
        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, device));
        printf("Running CUDA tests on: %s\n", prop.name);
    }

    void TearDown() override {
        CUDA_CHECK(cudaDeviceReset());
    }
};

TEST_F(FpCudaOps, AddSingle) {
    // Host data
    vec384 h_a, h_b, h_result, h_expected;
    set_vec384_cuda(h_a, 5);
    set_vec384_cuda(h_b, 7);
    set_vec384_cuda(h_expected, 12);

    // Device data
    vec384* d_a, * d_b, * d_result;
    CUDA_CHECK(cudaMalloc(&d_a, sizeof(vec384)));
    CUDA_CHECK(cudaMalloc(&d_b, sizeof(vec384)));
    CUDA_CHECK(cudaMalloc(&d_result, sizeof(vec384)));

    // Copy to device
    CUDA_CHECK(cudaMemcpy(d_a, h_a, sizeof(vec384), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b, sizeof(vec384), cudaMemcpyHostToDevice));

    // Launch kernel
    test_add_fp_kernel << <1, 1 >> > (d_result, d_a, d_b, 1);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy result back
    CUDA_CHECK(cudaMemcpy(h_result, d_result, sizeof(vec384), cudaMemcpyDeviceToHost));

    // Verify result
    EXPECT_TRUE(vec384_equal_cuda(h_result, h_expected))
        << "CUDA: 5 + 7 should equal 12";

    // Cleanup
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_result));
}

TEST_F(FpCudaOps, AddOverflow) {
    // Host data
    vec384 h_a, h_b, h_result, h_expected;

    // a = p - 1, b = 2, expect 1
    memcpy(h_a, BLS12_381_P, sizeof(vec384));
    for (int i = 0; i < LIMBS; ++i) h_a[i] -= (i == 0 ? 1 : 0);
    set_vec384_cuda(h_b, 2);
    set_vec384_cuda(h_expected, 1);

    // Device data
    vec384* d_a, * d_b, * d_result;
    CUDA_CHECK(cudaMalloc(&d_a, sizeof(vec384)));
    CUDA_CHECK(cudaMalloc(&d_b, sizeof(vec384)));
    CUDA_CHECK(cudaMalloc(&d_result, sizeof(vec384)));

    // Copy to device
    CUDA_CHECK(cudaMemcpy(d_a, h_a, sizeof(vec384), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b, sizeof(vec384), cudaMemcpyHostToDevice));

    // Launch kernel
    test_add_fp_kernel << <1, 1 >> > (d_result, d_a, d_b, 1);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy result back
    CUDA_CHECK(cudaMemcpy(h_result, d_result, sizeof(vec384), cudaMemcpyDeviceToHost));

    // Verify result
    EXPECT_TRUE(vec384_equal_cuda(h_result, h_expected))
        << "CUDA: (p-1) + 2 should equal 1";

    // Cleanup
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_result));
}

TEST_F(FpCudaOps, AddBatch) {
    const int N = 1024;  // Test with 1024 parallel additions

    // Host data
    vec384* h_a = new vec384[N];
    vec384* h_b = new vec384[N];
    vec384* h_results = new vec384[N];
    vec384* h_expected = new vec384[N];

    // Initialize test data
    for (int i = 0; i < N; i++) {
        set_vec384_cuda(h_a[i], i + 1);       // a[i] = i+1
        set_vec384_cuda(h_b[i], (i + 1) * 2); // b[i] = 2*(i+1)
        set_vec384_cuda(h_expected[i], (i + 1) * 3); // expected = 3*(i+1)
    }

    // Device data
    vec384* d_a, * d_b, * d_results;
    CUDA_CHECK(cudaMalloc(&d_a, N * sizeof(vec384)));
    CUDA_CHECK(cudaMalloc(&d_b, N * sizeof(vec384)));
    CUDA_CHECK(cudaMalloc(&d_results, N * sizeof(vec384)));

    // Copy to device
    CUDA_CHECK(cudaMemcpy(d_a, h_a, N * sizeof(vec384), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b, N * sizeof(vec384), cudaMemcpyHostToDevice));

    // Launch kernel
    int blockSize = 256;
    int gridSize = (N + blockSize - 1) / blockSize;
    test_add_fp_kernel << <gridSize, blockSize >> > (d_results, d_a, d_b, N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_results, d_results, N * sizeof(vec384), cudaMemcpyDeviceToHost));

    // Verify all results
    int errors = 0;
    for (int i = 0; i < N; i++) {
        if (!vec384_equal_cuda(h_results[i], h_expected[i])) {
            if (errors < 5) {  // Only print first 5 errors
                printf("Error at index %d:\n", i);
                printf("  Expected: "); print_vec384_cuda(h_expected[i]);
                printf("  Got:      "); print_vec384_cuda(h_results[i]);
            }
            errors++;
        }
    }

    EXPECT_EQ(errors, 0) << "Found " << errors << " errors in batch addition test";

    // Cleanup
    delete[] h_a;
    delete[] h_b;
    delete[] h_results;
    delete[] h_expected;
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_results));
}

TEST_F(FpCudaOps, AddBatchGenerated) {
    const int N = 512;

    // Host data for results
    vec384* h_results = new vec384[N];

    // Device data
    vec384* d_results;
    CUDA_CHECK(cudaMalloc(&d_results, N * sizeof(vec384)));

    // Launch kernel that generates test data internally
    int blockSize = 256;
    int gridSize = (N + blockSize - 1) / blockSize;
    test_add_fp_batch_kernel << <gridSize, blockSize >> > (d_results, N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_results, d_results, N * sizeof(vec384), cudaMemcpyDeviceToHost));

    // Verify results on CPU
    int errors = 0;
    for (int i = 0; i < N; i++) {
        vec384 expected_a, expected_b, expected_result;
        set_vec384_cuda(expected_a, i + 1);
        set_vec384_cuda(expected_b, (i + 1) * 2);

        // Compute expected result on CPU
        add_fp(expected_result, expected_a, expected_b);

        if (!vec384_equal_cuda(h_results[i], expected_result)) {
            if (errors < 3) {  // Only print first 3 errors
                printf("Error at index %d:\n", i);
                printf("  a:        "); print_vec384_cuda(expected_a);
                printf("  b:        "); print_vec384_cuda(expected_b);
                printf("  Expected: "); print_vec384_cuda(expected_result);
                printf("  Got:      "); print_vec384_cuda(h_results[i]);
            }
            errors++;
        }
    }

    EXPECT_EQ(errors, 0) << "Found " << errors << " errors in batch generated test";

    // Cleanup
    delete[] h_results;
    CUDA_CHECK(cudaFree(d_results));
}

// Test CPU vs GPU consistency
TEST_F(FpCudaOps, CPUvsGPUConsistency) {
    const int N = 100;

    // Host data
    vec384* h_a = new vec384[N];
    vec384* h_b = new vec384[N];
    vec384* h_cpu_results = new vec384[N];
    vec384* h_gpu_results = new vec384[N];

    // Generate random test data
    for (int i = 0; i < N; i++) {
        set_vec384_cuda(h_a[i], rand() % 10000);
        set_vec384_cuda(h_b[i], rand() % 10000);

        // Compute CPU result
        add_fp(h_cpu_results[i], h_a[i], h_b[i]);
    }

    // GPU computation
    vec384* d_a, * d_b, * d_results;
    CUDA_CHECK(cudaMalloc(&d_a, N * sizeof(vec384)));
    CUDA_CHECK(cudaMalloc(&d_b, N * sizeof(vec384)));
    CUDA_CHECK(cudaMalloc(&d_results, N * sizeof(vec384)));

    CUDA_CHECK(cudaMemcpy(d_a, h_a, N * sizeof(vec384), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b, N * sizeof(vec384), cudaMemcpyHostToDevice));

    int blockSize = 256;
    int gridSize = (N + blockSize - 1) / blockSize;
    test_add_fp_kernel << <gridSize, blockSize >> > (d_results, d_a, d_b, N);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(h_gpu_results, d_results, N * sizeof(vec384), cudaMemcpyDeviceToHost));

    // Compare CPU and GPU results
    int errors = 0;
    for (int i = 0; i < N; i++) {
        if (!vec384_equal_cuda(h_cpu_results[i], h_gpu_results[i])) {
            if (errors < 3) {
                printf("Mismatch at index %d:\n", i);
                printf("  CPU: "); print_vec384_cuda(h_cpu_results[i]);
                printf("  GPU: "); print_vec384_cuda(h_gpu_results[i]);
            }
            errors++;
        }
    }

    EXPECT_EQ(errors, 0) << "CPU and GPU results differ for " << errors << " cases";

    // Cleanup
    delete[] h_a;
    delete[] h_b;
    delete[] h_cpu_results;
    delete[] h_gpu_results;
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_results));
}