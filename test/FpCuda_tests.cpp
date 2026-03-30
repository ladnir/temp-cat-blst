#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <memory>
#include <random>
#include <cuda_runtime.h>

#include "blst/fields/Fp.h"
#include "blst/fields/Fp_cuda.h"
#include "blst/cuda/fields/parallel_prefix.h"
#include "blst/exports.h"
#include "InsecureRng.h"

using namespace blst;

#define LIMBS NLIMBS(384)


class FpCudaTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize random number generator
        //gen.seed(42); // Fixed seed for reproducible tests

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

    // Helper: Compare two Fp& for equality
    bool Fp_equal(const Fp& a, const Fp& b) {
        for (int i = 0; i < LIMBS; ++i)
            if (a[i] != b[i]) return false;
        return true;
    }

    // Helper: Print Fp for debugging
    void print_Fp(const Fp& v) {
        for (int i = LIMBS - 1; i >= 0; --i)
            printf("%016llx%s", (unsigned long long)v[i], i ? " " : "\n");
    }

    // Helper: Copy one Fp to another
    void fp_copy(Fp& dst, const Fp& src) {
        for (int i = 0; i < LIMBS; ++i)
            dst[i] = src[i];
    }
    void set_fp(Fp& v, uint64_t x) {
        for (int i = 0; i < LIMBS; ++i)
            v[i] = 0;
        memcpy(&v, &x, sizeof(x));
    }

    //// Helper: Generate random Fp element
    //void random_fp(Fp& v) {
    //    // Generate random limbs
    //    for (int i = 0; i < LIMBS; ++i) {
    //        v[i] = gen();
    //    }
    //    // Reduce modulo Fp to ensure valid element
    //    // Note: This is a simplified approach; real implementation would use proper modular reduction
    //    if (LIMBS > 0) {
    //        v[LIMBS - 1] &= 0x7FFFFFFFFFFFFFFFULL; // Ensure MSB is clear to avoid overflow
    //    }
    //}
    void print_Fp(const char* label, const Fp& v) {
        printf("%s: ", label);
        for (int i = LIMBS - 1; i >= 0; --i)
            printf("%016llx%s", (unsigned long long)v[i], i ? " " : "\n");
    }
    // Manually set all limbs to a specific pattern (good for triggering overflows)
    void set_fp_pattern(Fp& v, uint64_t pattern) {
        for (int i = 0; i < LIMBS; ++i) v[i] = pattern;
        // We must mask the top limb to ensure it is < P (approx) 
        // or allow it to be fully reduced later. 
        // For raw multiplication tests, we usually want valid field elements.
        // But to test the multiplier internals, unreduced inputs are sometimes useful.
        // Here we adhere to strict Fp < P by reducing.
        while (Fp::is_redc(v) == 0) {
            Fp::rshift(v, v, 1);
        }
    }
    void set_fp_val(Fp& v, uint64_t val) {
        for (int i = 0; i < LIMBS; ++i) v[i] = 0;
        v[0] = val;
    }

    void rand_fp(Fp& v, InsecureRng& rng)
    {
        for (u64 i = 0; i < v.num_limb(); ++i)
        {
            v[i] = rng.next_u64();
        }

        while (Fp::is_redc(v) == 0)
        {
            Fp::rshift(v, v, 1);
        }

        BLST_ASSERT(cc_check_mod_n<NLIMBS(384)>((byte*)v.data, get_BLST_P().data));
    }
    void rand_fp(Fp& v, int seed)
    {
        InsecureRng rng(seed);
        rand_fp(v, rng);
    }
    // Helper: Convert standard integer to Montgomery form
    void to_mont_fp(Fp& r, const Fp& a) {
        Fp::to_mont(r, a);
    }

    // Helper: Allocate device memory for Fp arrays
    Fp* allocate_device_fp_array(size_t num_elements) {
        Fp* d_ptr = nullptr;
        cudaError_t err = cudaMalloc((void**)&d_ptr, num_elements * sizeof(Fp));
        if (err != cudaSuccess) {
            return nullptr;
        }
        return d_ptr;
    }

    // Helper: Free device memory
    void free_device_fp_array(Fp* d_ptr) {
        if (d_ptr) {
            cudaFree(d_ptr);
        }
    }

    // Helper: Copy host Fp array to device
    batched_error_code copy_host_to_device_fp(Fp* d_dst, const Fp* h_src, size_t num_elements) {
        if (!d_dst || !h_src) {
            return BATCHED_ERROR_NULL_POINTER;
        }

        cudaError_t err = cudaMemcpy(d_dst, h_src, num_elements * sizeof(Fp), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
        }

        return BATCHED_SUCCESS;
    }

    // Helper: Copy device Fp array to host
    batched_error_code copy_device_to_host_fp(Fp* h_dst, const Fp* d_src, size_t num_elements) {
        if (!h_dst || !d_src) {
            return BATCHED_ERROR_NULL_POINTER;
        }

        cudaError_t err = cudaMemcpy(h_dst, d_src, num_elements * sizeof(Fp), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
        }

        return BATCHED_SUCCESS;
    }

    // Helper: Compare vector results against CPU reference
    bool compare_fp_vectors(const Fp* vec_result, const Fp* cpu_result, size_t num, const char* operation) {
        bool all_equal = true;
        for (size_t i = 0; i < num; i++) {
            if (!Fp_equal(vec_result[i], cpu_result[i])) {
                printf("Mismatch in %s at index %zu:\n", operation, i);
                printf("Vector result: "); print_Fp(vec_result[i]);
                printf("CPU result:    "); print_Fp(cpu_result[i]);
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
    // CUDA Management Wrappers
    Fp* d_alloc(size_t n) {
        Fp* p = nullptr;
        cudaMalloc((void**)&p, n * sizeof(Fp));
        return p;
    }

    void to_dev(Fp* d, const Fp* h, size_t n) {
        cudaMemcpy(d, h, n * sizeof(Fp), cudaMemcpyHostToDevice);
    }

    void to_host(Fp* h, const Fp* d, size_t n) {
        cudaMemcpy(h, d, n * sizeof(Fp), cudaMemcpyDeviceToHost);
    }
    //std::mt19937_64 gen;
};

// ============================================================================
// Test 1: Host Logic Verification
// Checks the C++ implementation of the FIOS/CIOS logic directly on CPU.
// This runs faster and separates logic bugs from CUDA memory bugs.
// ============================================================================
TEST_F(FpCudaTest, Host_PTX_Logic_Check) {
    // We test specific edge cases that cause carry propagation
    std::vector<std::pair<std::string, uint64_t>> patterns = {
        {"Zero", 0x0000000000000000ULL},
        {"One",  0x0000000000000001ULL},
        {"Max",  0xFFFFFFFFFFFFFFFFULL},
        {"Alt1", 0xAAAAAAAAAAAAAAAAULL},
        {"Alt2", 0x5555555555555555ULL},
        {"High", 0x8000000000000000ULL}
    };

    for (const auto& p1 : patterns) {
        for (const auto& p2 : patterns) {
            Fp a, b, expected, result;

            // Construct inputs
            set_fp_pattern(a, p1.second);
            set_fp_pattern(b, p2.second);

            // 1. Golden Reference (BLST optimized asm)
            Fp::mul(expected, a, b);

            // 2. Your C++ implementation of the CUDA logic
            ptx_mul_mont_u64<Fp>(result.data, a.data, b.data);

            if (!Fp_equal(expected, result)) {
                printf("FAILURE: Logic Check %s * %s\n", p1.first.c_str(), p2.first.c_str());
                print_Fp("A       ", a);
                print_Fp("B       ", b);
                print_Fp("Expected", expected);
                print_Fp("Got     ", result);
                FAIL() << "Host logic mismatch";
            }
        }
    }
}

// ============================================================================
// Test 2: CUDA Vector Multiply with Edge Cases
// Specifically targets the carry bug by using saturated limbs (All 1s).
// ============================================================================
TEST_F(FpCudaTest, VectorMul_EdgeCases) {
    // Define specific test cases
    struct TestCase { Fp a; Fp b; std::string desc; };
    std::vector<TestCase> cases;

    Fp zero, one, max_limbs, alt_limbs;
    set_fp_val(zero, 0);
    set_fp_val(one, 1);
    set_fp_pattern(max_limbs, ~0ull);
    set_fp_pattern(alt_limbs, 0xAAAAAAAAAAAAAAAAull);

    // 1. Identity
    cases.push_back({ one, one, "1 * 1" });
    cases.push_back({ max_limbs, one, "Max * 1" });

    // 2. Zero prop
    cases.push_back({ max_limbs, zero, "Max * 0" });

    // 3. Max Saturation (The Carry Killer)
    // A * B where both are dense with 1s creates maximum partial products
    // and maximum carry propagation chains.
    cases.push_back({ max_limbs, max_limbs, "Max * Max" });

    // 4. Alternating bits (Stress tests internal adders)
    cases.push_back({ alt_limbs, alt_limbs, "Alt * Alt" });
    cases.push_back({ max_limbs, alt_limbs, "Max * Alt" });

    // 5. Random Dense (Simulate high entropy)
    InsecureRng rng(12345);
    for (int i = 0; i < 50; ++i) {
        Fp ra, rb;
        rand_fp(ra, rng);
        rand_fp(rb, rng);
        cases.push_back({ ra, rb, "Random Dense" });
    }

    size_t N = cases.size();
    auto h_a = std::make_unique<Fp[]>(N);
    auto h_b = std::make_unique<Fp[]>(N);
    auto h_res = std::make_unique<Fp[]>(N);
    auto cpu_res = std::make_unique<Fp[]>(N);

    // Setup Host Data
    for (size_t i = 0; i < N; ++i) {
        h_a[i] = cases[i].a;
        h_b[i] = cases[i].b;

        // Compute Gold Standard
        Fp::mul(cpu_res[i], h_a[i], h_b[i]);

        // Double check Host Logic (Debugging the debugger)
        Fp host_logic_res;
        ptx_mul_mont_u64<Fp>(host_logic_res.data, h_a[i].data, h_b[i].data);
        if (!Fp_equal(cpu_res[i], host_logic_res)) {
            printf("CRITICAL: Host logic mismatch at index %zu (%s)\n", i, cases[i].desc.c_str());
            FAIL();
        }
    }

    // Setup Device
    Fp* d_a = d_alloc(N);
    Fp* d_b = d_alloc(N);
    Fp* d_r = d_alloc(N);

    to_dev(d_a, h_a.get(), N);
    to_dev(d_b, h_b.get(), N);

    // Run Kernel
    batched_error_code err = mul_fp_vec(d_r, d_a, d_b, N);
    ASSERT_EQ(err, BATCHED_SUCCESS);

    // Verify
    to_host(h_res.get(), d_r, N);

    for (size_t i = 0; i < N; ++i) {
        if (!Fp_equal(cpu_res[i], h_res[i])) {
            printf("\n[MISMATCH] Index %zu: %s\n", i, cases[i].desc.c_str());
            print_Fp("Input A ", h_a[i]);
            print_Fp("Input B ", h_b[i]);
            print_Fp("GPU Res ", h_res[i]);
            print_Fp("CPU Res ", cpu_res[i]);

            // Check if GPU result matches the buggy logic?
            // (Optional: helps identify if the GPU code is the old version)

            FAIL() << "GPU Multiplication failed";
        }
    }

    cudaFree(d_a); cudaFree(d_b); cudaFree(d_r);
}

// ============================================================================
// Test 3: High Volume Random (Existing test, slightly tuned)
// ============================================================================
TEST_F(FpCudaTest, VectorMul_MassiveRandom) {
    const size_t NUM_ELEMENTS = 10000;
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_res = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_res = std::make_unique<Fp[]>(NUM_ELEMENTS);

    InsecureRng rng(999); // Different seed

    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        rand_fp(h_a[i], rng);
        rand_fp(h_b[i], rng);
        Fp::mul(cpu_res[i], h_a[i], h_b[i]);
    }

    Fp* d_a = d_alloc(NUM_ELEMENTS);
    Fp* d_b = d_alloc(NUM_ELEMENTS);
    Fp* d_r = d_alloc(NUM_ELEMENTS);

    to_dev(d_a, h_a.get(), NUM_ELEMENTS);
    to_dev(d_b, h_b.get(), NUM_ELEMENTS);

    mul_fp_vec(d_r, d_a, d_b, NUM_ELEMENTS);
    to_host(h_res.get(), d_r, NUM_ELEMENTS);

    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        if (!Fp_equal(h_res[i], cpu_res[i])) {
            print_Fp("Expected", cpu_res[i]);
            print_Fp("Got     ", h_res[i]);
            FAIL() << "Random mismatch at index " << i;
        }
    }

    cudaFree(d_a); cudaFree(d_b); cudaFree(d_r);
}

TEST_F(FpCudaTest, par_PTX_Logic_Check) {
    // We test specific edge cases that cause carry propagation
    std::vector<std::pair<std::string, uint64_t>> patterns = {
        {"Zero", 0x0000000000000000ULL},
        {"One",  0x0000000000000001ULL},
        {"Max",  0xFFFFFFFFFFFFFFFFULL},
        {"Alt1", 0xAAAAAAAAAAAAAAAAULL},
        {"Alt2", 0x5555555555555555ULL},
        {"High", 0x8000000000000000ULL}
    };

    for (const auto& p1 : patterns) {
        for (const auto& p2 : patterns) {
            Fp a, b, expected, result;

            // Construct inputs
            set_fp_pattern(a, p1.second);
            set_fp_pattern(b, p2.second);

            // 1. Golden Reference (BLST optimized asm)
            Fp::mul(expected, a, b);

            // 2. Your C++ implementation of the CUDA logic
            par_mul_mont_u64<Fp>(result.data, a.data, b.data);

            if (!Fp_equal(expected, result)) {
                printf("FAILURE: Logic Check %s * %s\n", p1.first.c_str(), p2.first.c_str());
                print_Fp("A       ", a);
                print_Fp("B       ", b);
                print_Fp("Expected", expected);
                print_Fp("Got     ", result);
                FAIL() << "Host logic mismatch";
            }
        }
    }
}


TEST_F(FpCudaTest, VectorAdd) {
    const size_t NUM_ELEMENTS = 55;
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

	InsecureRng rng;
    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        //rand_fp(h_a[i], 5 + i);
        //rand_fp(h_b[i], 7 + i * 2);
		rand_fp(h_a[i], rng);
		rand_fp(h_b[i], rng);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::add(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_b = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        rand_fp(h_a[i], 10 + i * 2);
        rand_fp(h_b[i], 3 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::sub(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_b = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data with Montgomery form conversion
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp temp_a, temp_b;
        rand_fp(temp_a, 3 + i);
        rand_fp(temp_b, 4 + i);

        // Convert to Montgomery form
        to_mont_fp(h_a[i], temp_a);
        to_mont_fp(h_b[i], temp_b);
    }
    //printf("\n\n");

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::mul(cpu_result[i], h_a[i], h_b[i]);
        Fp rr;
		ptx_mul_mont_u64<Fp>(rr.data, h_a[i].data, h_b[i].data);
        if (!Fp_equal(cpu_result[i], rr)) {
            printf("CPU mul mismatch at index %zu:\n", i);
            printf("Expected: "); print_Fp(cpu_result[i]);
            printf("Got:      "); print_Fp(rr);
		}
    }
    //printf("\n\n");
    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_b = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data with Montgomery form conversion
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp temp_a;
        rand_fp(temp_a, 5 + i);

        // Convert to Montgomery form
        to_mont_fp(h_a[i], temp_a);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::sqr(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result_true = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result_true = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result_false = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result_false = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        rand_fp(h_a[i], 10 + i);
    }

    // CPU reference computation with flag = 1
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::cneg(cpu_result_true[i], h_a[i], 1);
    }

    // CPU reference computation with flag = 0
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::cneg(cpu_result_false[i], h_a[i], 0);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result_true = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result_false = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        rand_fp(h_a[i], 2 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::mul_by_3(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        rand_fp(h_a[i], 2 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::mul_by_8(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        rand_fp(h_a[i], 8 + i * 2);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::div_by_2(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        rand_fp(h_a[i], 1 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::lshift(cpu_result[i], h_a[i], shift_count);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        rand_fp(h_a[i], 9 + i);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::to_mont(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data in Montgomery form
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp temp;
        rand_fp(temp, 9 + i);
        Fp::to_mont(h_a[i], temp);
    }

    // CPU reference computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::from_mont(cpu_result[i], h_a[i]);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        rand_fp(h_a[i], 42 + i);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
    auto h_a = std::make_unique<Fp[]>(MIXED_SIZE);
    auto h_b = std::make_unique<Fp[]>(MIXED_SIZE);
    auto h_result = std::make_unique<Fp[]>(MIXED_SIZE);
    auto cpu_result = std::make_unique<Fp[]>(MIXED_SIZE);

    // Initialize test data
    for (size_t i = 0; i < MIXED_SIZE; i++) {
        rand_fp(h_a[i], 3 + i);
        rand_fp(h_b[i], 5 + i * 2);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(MIXED_SIZE);
    Fp* d_b = allocate_device_fp_array(MIXED_SIZE);
    Fp* d_result = allocate_device_fp_array(MIXED_SIZE);
    Fp* d_temp = allocate_device_fp_array(MIXED_SIZE);

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
                    Fp temp_sum, temp_mul3;
                    Fp::add(temp_sum, h_a[i], h_b[i]);
                    Fp::mul_by_3(temp_mul3, temp_sum);
                    Fp::div_by_2(cpu_result[i], temp_mul3);
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
    auto h_zeros = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_fp(h_zeros[i], 0);
    }

    // Allocate device memory
    Fp* d_zeros = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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
            EXPECT_TRUE(Fp_equal(h_result[i], h_zeros[i])) << "0 + 0 should equal 0 at index " << i;
        }
    }

    // Test with one
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        set_fp(h_zeros[i], 1);
    }

    expect_success(copy_host_to_device_fp(d_zeros, h_zeros.get(), NUM_ELEMENTS), "copy ones to device");

    batched_error_code add_one_result = add_fp_vec(d_result, d_zeros, d_zeros, NUM_ELEMENTS);
    expect_success(add_one_result, "add_fp_vec with ones");

    if (add_one_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy add ones result to host");
        for (size_t i = 0; i < NUM_ELEMENTS; i++) {
            Fp expected;
            set_fp(expected, 2);
            EXPECT_TRUE(Fp_equal(h_result[i], expected)) << "1 + 1 should equal 2 at index " << i;
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
    Fp* d_dummy = allocate_device_fp_array(NUM_ELEMENTS);
    ASSERT_NE(d_dummy, nullptr) << "Failed to allocate device memory for dummy array";

    batched_error_code zero_size_result = add_fp_vec(d_dummy, d_dummy, d_dummy, 0);
    EXPECT_EQ(zero_size_result, BATCHED_SUCCESS)
        << "Expected SUCCESS for zero element count (zero-size operations should be OK)";

    free_device_fp_array(d_dummy);
}

TEST_F(FpCudaTest, VectorSmallValues) {
    const size_t NUM_ELEMENTS = 5;
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Test with small known values
    const uint64_t TEST_VALUES[] = { 1, 2, 5, 10, 42 };

    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        rand_fp(h_a[i], TEST_VALUES[i]);
        rand_fp(h_b[i], TEST_VALUES[(i + 1) % NUM_ELEMENTS]);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_b = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

    ASSERT_NE(d_a, nullptr) << "Failed to allocate device memory for array a";
    ASSERT_NE(d_b, nullptr) << "Failed to allocate device memory for array b";
    ASSERT_NE(d_result, nullptr) << "Failed to allocate device memory for result";

    // Copy data to device
    expect_success(copy_host_to_device_fp(d_a, h_a.get(), NUM_ELEMENTS), "copy h_a to device");
    expect_success(copy_host_to_device_fp(d_b, h_b.get(), NUM_ELEMENTS), "copy h_b to device");

    // Test addition
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::add(cpu_result[i], h_a[i], h_b[i]);
    }
    batched_error_code add_result = add_fp_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(add_result, "small values add");
    if (add_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy add result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "small values add"));
    }

    // Test subtraction
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::sub(cpu_result[i], h_a[i], h_b[i]);
    }
    batched_error_code sub_result = sub_fp_vec(d_result, d_a, d_b, NUM_ELEMENTS);
    expect_success(sub_result, "small values sub");
    if (sub_result == BATCHED_SUCCESS) {
        expect_success(copy_device_to_host_fp(h_result.get(), d_result, NUM_ELEMENTS), "copy sub result to host");
        EXPECT_TRUE(compare_fp_vectors(h_result.get(), cpu_result.get(), NUM_ELEMENTS, "small values sub"));
    }

    // Test multiplication by 3
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::mul_by_3(cpu_result[i], h_a[i]);
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
    auto h_a = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_b = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto h_result = std::make_unique<Fp[]>(NUM_ELEMENTS);
    auto cpu_result = std::make_unique<Fp[]>(NUM_ELEMENTS);

    // Generate test data
	InsecureRng rng;
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        rand_fp(h_a[i], rng);
        rand_fp(h_b[i], rng);
    }

    // This test primarily validates correctness at scale
    // Performance measurements would require timing infrastructure

    // CPU computation
    for (size_t i = 0; i < NUM_ELEMENTS; i++) {
        Fp::add(cpu_result[i], h_a[i], h_b[i]);
    }

    // Allocate device memory
    Fp* d_a = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_b = allocate_device_fp_array(NUM_ELEMENTS);
    Fp* d_result = allocate_device_fp_array(NUM_ELEMENTS);

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