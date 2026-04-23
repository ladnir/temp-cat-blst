#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <array>
#include <iostream>
#include "config.h"
#include "blst/curve/GT.h"
#include "blst/curve/GT_cuda.h"
#include "blst/sign/aggregate.h"
#include "blst/curve/point.h"
#include "blst/curve/pairing.h"
#include "InsecureRng.h"

using namespace blst;

// CUDA kernel wrapper macros for device execution
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while(0)

namespace
{
    // Helper function to create test GT elements on host
    void create_test_gt_elements(GT* elements, size_t num, InsecureRng& rng) {
        G1Point g1_gen = G1Point::generator();
        G2Point g2_gen = G2Point::generator();
        G1AffinePoint g1_affine;
        G2AffinePoint g2_affine;

        G1Point::to_affine(g1_affine, g1_gen);
        G2Point::to_affine(g2_affine, g2_gen);

        for (size_t i = 0; i < num; i++) {
            // Create different GT elements by varying the scalar multiplier
            G1Point g1_mult;
            unsigned char scalar[32] = { 0 };
            scalar[0] = (rng.uniform_u64() % 254) + 1; // Avoid 0
            scalar[1] = (rng.uniform_u64() % 255);
            scalar[2] = (i + 1) % 255;
            if (scalar[0] == 0) scalar[0] = 1;

            G1Point::mult(g1_mult, g1_gen, scalar, 24);
            G1Point::to_affine(g1_affine, g1_mult);

            // Compute pairing to get GT element
            GT::miller_loop(elements[i], &g2_affine, &g1_affine);
            GT::final_mult(elements[i], elements[i]);
        }
    }

    // Helper function to check batched_error_code
    bool check_batched_success(batched_error_code code)
    {
        if (code != BATCHED_SUCCESS) {
            fprintf(stderr, "Batched operation failed with error code: %d\n", code);
            return false;
        }
        return true;
    }

    // Helper function to compare CUDA results with CPU results
    template<typename T>
    bool compare_results(const T* cuda_results, const T* cpu_results, size_t num,
        auto equality_func) {
        for (size_t i = 0; i < num; i++) {
            if (!equality_func(cuda_results[i], cpu_results[i])) {
                printf("Mismatch at index %zu\n", i);
                return false;
            }
        }
        return true;
    }

    // Equality functions for different types
    bool gt_element_equal(const GT& a, const GT& b) {
        return GT::equal(a, b);
    }

    bool bool_equal(const bool& a, const bool& b) {
        return a == b;
    }

    std::array<u64, 4> scalar_from_modulus_minus(std::array<u64, 4> delta)
    {
        const Fr& modulus = get_BLST_r();
        std::array<u64, 4> out = { modulus[0], modulus[1], modulus[2], modulus[3] };

        u64 borrow = 0;
        for (size_t i = 0; i < out.size(); ++i) {
            const u64 sub = delta[i] + borrow;
            borrow = (out[i] < sub) ? 1 : 0;
            out[i] -= sub;
        }

        return out;
    }

    // Global flag for CUDA availability
    static bool cuda_device_available = false;
}

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

// CUDA test cases for GT operations using the GT_cuda interface
TEST(GTOpsCUDA, Identity) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 10;
    GT* d_elements;
    bool* d_results;
    GT h_elements[num_elements];
    bool h_cuda_results[num_elements], h_cpu_results[num_elements];

    // Create test data including identity element
    InsecureRng rng;
    create_test_gt_elements(h_elements, num_elements - 1, rng);

    // Add identity element
	h_elements[num_elements - 1] = GT::identity();

    // CPU computation
    for (size_t i = 0; i < num_elements; i++) {
        h_cpu_results[i] = GT::is_identity(h_elements[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_elements, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_results, num_elements * sizeof(bool)));

    CUDA_CHECK(cudaMemcpy(d_elements, h_elements, num_elements * sizeof(GT), cudaMemcpyHostToDevice));

    batched_error_code result = is_one_gt_vec(d_results, d_elements, num_elements);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_results, d_results, num_elements * sizeof(bool), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_results, h_cpu_results, num_elements, bool_equal));

    // Verify that the last element (identity) is correctly identified
    EXPECT_TRUE(h_cuda_results[num_elements - 1]) << "Identity element should be identified as one";

    // Cleanup
    CUDA_CHECK(cudaFree(d_elements));
    CUDA_CHECK(cudaFree(d_results));
}

TEST(GTOpsCUDA, Multiplication) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 10;
    GT* d_a, * d_b, * d_result;
    GT h_a[num_elements], h_b[num_elements];
    GT h_cuda_result[num_elements], h_cpu_result[num_elements];

    // Create test data
    InsecureRng rng;
    create_test_gt_elements(h_a, num_elements, rng);
    create_test_gt_elements(h_b, num_elements, rng);

    // CPU computation
    for (size_t i = 0; i < num_elements; i++) {
        GT::add(h_cpu_result[i], h_a[i], h_b[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_a, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_b, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_result, num_elements * sizeof(GT)));

    CUDA_CHECK(cudaMemcpy(d_a, h_a, num_elements * sizeof(GT), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b, num_elements * sizeof(GT), cudaMemcpyHostToDevice));

    batched_error_code result = mult_gt_vec(d_result, d_a, d_b, num_elements);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_elements, gt_element_equal));

    // Test commutativity on one element
    GT h_commute_test;
    GT::add(h_commute_test, h_b[0], h_a[0]);
    EXPECT_TRUE(GT::equal(h_commute_test, h_cuda_result[0]))
        << "GT multiplication should be commutative";

    // Cleanup
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_result));
}

TEST(GTOpsCUDA, Inversion) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 10;
    GT* d_elements, * d_result;
    GT h_elements[num_elements];
    GT h_cuda_result[num_elements], h_cpu_result[num_elements];

    // Create test data
    InsecureRng rng;
    create_test_gt_elements(h_elements, num_elements, rng);

    // CPU computation
    for (size_t i = 0; i < num_elements; i++) {
        GT::neg(h_cpu_result[i], h_elements[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_elements, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_result, num_elements * sizeof(GT)));

    CUDA_CHECK(cudaMemcpy(d_elements, h_elements, num_elements * sizeof(GT), cudaMemcpyHostToDevice));

    batched_error_code result = inverse_gt_vec(d_result, d_elements, num_elements);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_elements, gt_element_equal));

    // Test that element * inverse = identity
    GT product, identity;
	identity = GT::identity();
    GT::add(product, h_elements[0], h_cuda_result[0]);
    EXPECT_TRUE(GT::equal(product, identity))
        << "Element * inverse should equal identity";

    // Cleanup
    CUDA_CHECK(cudaFree(d_elements));
    CUDA_CHECK(cudaFree(d_result));
}

TEST(GTOpsCUDA, Squaring) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 10;
    GT* d_elements, * d_result, * d_cyclo_result;
    GT h_elements[num_elements];
    GT h_cuda_result[num_elements], h_cpu_result[num_elements];
    GT h_cuda_cyclo_result[num_elements], h_cpu_cyclo_result[num_elements];

    // Create test data
    InsecureRng rng;
    create_test_gt_elements(h_elements, num_elements, rng);

    // CPU computation - regular squaring
    for (size_t i = 0; i < num_elements; i++) {
        GT::doubling(h_cpu_result[i], h_elements[i]);
    }

    // CPU computation - cyclotomic squaring
    for (size_t i = 0; i < num_elements; i++) {
        GT::doubling(h_cpu_cyclo_result[i], h_elements[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_elements, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_result, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_cyclo_result, num_elements * sizeof(GT)));

    CUDA_CHECK(cudaMemcpy(d_elements, h_elements, num_elements * sizeof(GT), cudaMemcpyHostToDevice));

    // Test regular squaring
    batched_error_code result = square_gt_vec(d_result, d_elements, num_elements);
    EXPECT_TRUE(check_batched_success(result));

    // Test cyclotomic squaring
    batched_error_code cyclo_result = cyclotomic_square_gt_vec(d_cyclo_result, d_elements, num_elements);
    EXPECT_TRUE(check_batched_success(cyclo_result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_cuda_cyclo_result, d_cyclo_result, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_elements, gt_element_equal));
    EXPECT_TRUE(compare_results(h_cuda_cyclo_result, h_cpu_cyclo_result, num_elements, gt_element_equal));

    // Test that squaring equals multiplication by self
    GT mult_by_self;
    GT::add(mult_by_self, h_elements[0], h_elements[0]);
    EXPECT_TRUE(GT::equal(mult_by_self, h_cuda_result[0]))
        << "Squaring should equal multiplication by self";

    // For GT elements, cyclotomic and regular squaring should be equal
    EXPECT_TRUE(GT::equal(h_cuda_result[0], h_cuda_cyclo_result[0]))
        << "Cyclotomic and regular squaring should be equal for GT elements";

    // Cleanup
    CUDA_CHECK(cudaFree(d_elements));
    CUDA_CHECK(cudaFree(d_result));
    CUDA_CHECK(cudaFree(d_cyclo_result));
}

TEST(GTOpsCUDA, Exponentiation) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 10;
    const size_t exp_bits = 64;
    const size_t exp_bytes = (exp_bits + 7) / 8;

    GT* d_bases, * d_result;
    unsigned char* d_exponents;
    GT h_bases[num_elements];
    unsigned char h_exponents[num_elements * exp_bytes];
    GT h_cuda_result[num_elements], h_cpu_result[num_elements];

    // Create test data
    InsecureRng rng;
    create_test_gt_elements(h_bases, num_elements, rng);

    // Create small test exponents
    for (size_t i = 0; i < num_elements; i++) {
        memset(&h_exponents[i * exp_bytes], 0, exp_bytes);
        h_exponents[i * exp_bytes] = (rng.uniform_u64() % 254) + 2; // Avoid 0 and 1
        h_exponents[i * exp_bytes + 1] = (rng.uniform_u64() % 255);
    }

    // CPU computation
    for (size_t i = 0; i < num_elements; i++) {
        GT::mult(h_cpu_result[i], h_bases[i], &h_exponents[i * exp_bytes], exp_bits);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_bases, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_result, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_exponents, num_elements * exp_bytes));

    CUDA_CHECK(cudaMemcpy(d_bases, h_bases, num_elements * sizeof(GT), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_exponents, h_exponents, num_elements * exp_bytes, cudaMemcpyHostToDevice));

    batched_error_code result = pow_gt_vec(d_result, d_bases, d_exponents, exp_bits, num_elements);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_elements, gt_element_equal));

    // Test special cases - exponent 0 should give identity
    unsigned char zero_exp[exp_bytes] = { 0 };
    GT zero_result, identity;
	identity = GT::identity();
    GT::mult(zero_result, h_bases[0], zero_exp, exp_bits);
    EXPECT_TRUE(GT::equal(zero_result, identity))
        << "Exponentiation with 0 should give identity";

    // Test exponent 1 should give the base
    unsigned char one_exp[exp_bytes] = { 1 };
    GT one_result;
    GT::mult(one_result, h_bases[0], one_exp, exp_bits);
    EXPECT_TRUE(GT::equal(one_result, h_bases[0]))
        << "Exponentiation with 1 should give the base";

    // Cleanup
    CUDA_CHECK(cudaFree(d_bases));
    CUDA_CHECK(cudaFree(d_result));
    CUDA_CHECK(cudaFree(d_exponents));
}

TEST(GTOpsCUDA, Bilinearity) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_pairs = 5;
    GT* d_result;
    G1AffinePoint* d_g1_points;
    G2AffinePoint* d_g2_points;

    G1AffinePoint h_g1_points[num_pairs];
    G2AffinePoint h_g2_points[num_pairs];
    GT h_cuda_result[num_pairs], h_cpu_result[num_pairs];

    // Create test point pairs
    G1Point g1_gen = G1Point::generator();
    G2Point g2_gen = G2Point::generator();
    InsecureRng rng;

    for (size_t i = 0; i < num_pairs; i++) {
        G1Point g1_mult;
        G2Point g2_mult;
        unsigned char scalar[32] = { 0 };
        scalar[0] = (rng.uniform_u64() % 254) + 1;
        scalar[1] = (rng.uniform_u64() % 255);

        G1Point::mult(g1_mult, g1_gen, scalar, 16);
        G2Point::mult(g2_mult, g2_gen, scalar, 16);

        G1Point::to_affine(h_g1_points[i], g1_mult);
        G2Point::to_affine(h_g2_points[i], g2_mult);
    }

    // CPU computation
    for (size_t i = 0; i < num_pairs; i++) {
        GT::miller_loop(h_cpu_result[i],&h_g2_points[i], &h_g1_points[i]);
        GT::final_mult(h_cpu_result[i], h_cpu_result[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_g1_points, num_pairs * sizeof(G1AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_g2_points, num_pairs * sizeof(G2AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_result, num_pairs * sizeof(GT)));

    CUDA_CHECK(cudaMemcpy(d_g1_points, h_g1_points, num_pairs * sizeof(G1AffinePoint), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_g2_points, h_g2_points, num_pairs * sizeof(G2AffinePoint), cudaMemcpyHostToDevice));

    batched_error_code result = pairing_gt_vec(d_result, d_g2_points, d_g1_points, num_pairs);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_pairs * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_pairs, gt_element_equal));

    // Test bilinearity property: e(aP, Q) = e(P, aQ) = e(P,Q)^a
    G1Point g1_double;
    G2Point g2_double;
    G1Point::doubling(g1_double, g1_gen);
    G2Point::doubling(g2_double, g2_gen);

    G1AffinePoint g1_affine, g1_double_affine;
    G2AffinePoint g2_affine, g2_double_affine;
    G1Point::to_affine(g1_affine, g1_gen);
    G2Point::to_affine(g2_affine, g2_gen);
    G1Point::to_affine(g1_double_affine, g1_double);
    G2Point::to_affine(g2_double_affine, g2_double);

    GT e_2P_Q, e_P_2Q, e_PQ_squared, base_pairing;

    // e(P, Q)
    GT::miller_loop(base_pairing, &g2_affine, &g1_affine);
    GT::final_mult(base_pairing, base_pairing);

    // e(2P, Q)
    GT::miller_loop(e_2P_Q, &g2_affine, &g1_double_affine);
    GT::final_mult(e_2P_Q, e_2P_Q);

    // e(P, 2Q)
    GT::miller_loop(e_P_2Q, &g2_double_affine, &g1_affine);
    GT::final_mult(e_P_2Q, e_P_2Q);

    // e(P,Q)^2
    GT::doubling(e_PQ_squared, base_pairing);

    // Test bilinearity
    EXPECT_TRUE(GT::equal(e_2P_Q, e_P_2Q))
        << "e(2P, Q) should equal e(P, 2Q)";
    EXPECT_TRUE(GT::equal(e_2P_Q, e_PQ_squared))
        << "e(2P, Q) should equal e(P,Q)^2";

    // Cleanup
    CUDA_CHECK(cudaFree(d_g1_points));
    CUDA_CHECK(cudaFree(d_g2_points));
    CUDA_CHECK(cudaFree(d_result));
}

TEST(GTOpsCUDA, GroupMembership) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 10;
    GT* d_elements;
    bool* d_results;
    GT h_elements[num_elements];
    bool h_cuda_results[num_elements], h_cpu_results[num_elements];

    // Create test data - mix of valid GT elements and invalid ones
    InsecureRng rng;
    create_test_gt_elements(h_elements, num_elements - 2, rng);

    // Add some potentially invalid elements
    memcpy(&h_elements[num_elements - 2], &h_elements[0], sizeof(GT));
    // Modify slightly to make it potentially not in group
    *(unsigned char*)&h_elements[num_elements - 2] ^= 1;

    memcpy(&h_elements[num_elements - 1], &h_elements[1], sizeof(GT));
    // Modify more significantly
    *(unsigned char*)&h_elements[num_elements - 1] ^= 0xFF;

    // CPU computation
    for (size_t i = 0; i < num_elements; i++) {
        h_cpu_results[i] = GT::in_group(h_elements[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_elements, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_results, num_elements * sizeof(bool)));

    CUDA_CHECK(cudaMemcpy(d_elements, h_elements, num_elements * sizeof(GT), cudaMemcpyHostToDevice));

    batched_error_code result = in_group_gt_vec(d_results, d_elements, num_elements);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_results, d_results, num_elements * sizeof(bool), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_results, h_cpu_results, num_elements, bool_equal));

    // Valid GT elements should be in the group
    for (size_t i = 0; i < num_elements - 2; i++) {
        EXPECT_TRUE(h_cuda_results[i]) << "Valid GT element " << i << " should be in group";
    }

    // Cleanup
    CUDA_CHECK(cudaFree(d_elements));
    CUDA_CHECK(cudaFree(d_results));
}

TEST(GTOpsCUDA, Conjugate) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 10;
    GT* d_elements, * d_result;
    GT h_elements[num_elements];
    GT h_cuda_result[num_elements], h_cpu_result[num_elements];

    // Create test data
    InsecureRng rng;
    create_test_gt_elements(h_elements, num_elements, rng);

    // CPU computation
    for (size_t i = 0; i < num_elements; i++) {
        memcpy(&h_cpu_result[i], &h_elements[i], sizeof(GT));
        GT::conjugate(h_cpu_result[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_elements, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_result, num_elements * sizeof(GT)));

    CUDA_CHECK(cudaMemcpy(d_elements, h_elements, num_elements * sizeof(GT), cudaMemcpyHostToDevice));

    batched_error_code result = conjugate_gt_vec(d_result, d_elements, num_elements);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_elements, gt_element_equal));

    // Test conjugate properties: conjugate(conjugate(x)) = x
    GT double_conjugate;
    memcpy(&double_conjugate, &h_cuda_result[0], sizeof(GT));
    GT::conjugate(double_conjugate);
    EXPECT_TRUE(GT::equal(double_conjugate, h_elements[0]))
        << "Conjugate of conjugate should equal original";

    // Test that conjugate is also in GT group
    bool conjugate_in_group = GT::in_group(h_cuda_result[0]);
    EXPECT_TRUE(conjugate_in_group) << "Conjugate should be in GT group";

    // Cleanup
    CUDA_CHECK(cudaFree(d_elements));
    CUDA_CHECK(cudaFree(d_result));
}

TEST(GTOpsCUDA, Serialization) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 5;
    GT* d_elements, * d_deserialized;
    unsigned char* d_serialized;
    GT h_elements[num_elements];
    GT h_cuda_deserialized[num_elements], h_cpu_deserialized[num_elements];
    unsigned char h_serialized[num_elements * 576]; // 576 = 12 * 48 bytes

    // Create test data
    InsecureRng rng;
    create_test_gt_elements(h_elements, num_elements, rng);

    // CPU computation
    for (size_t i = 0; i < num_elements; i++) {
        Fp12::to_bendian(&h_serialized[i * 576], h_elements[i].data);
        Fp12::from_bendian(h_cpu_deserialized[i].data, &h_serialized[i * 576]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_elements, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_serialized, num_elements * 576));
    CUDA_CHECK(cudaMalloc(&d_deserialized, num_elements * sizeof(GT)));

    CUDA_CHECK(cudaMemcpy(d_elements, h_elements, num_elements * sizeof(GT), cudaMemcpyHostToDevice));

    batched_error_code serialize_result = serialize_gt_vec(d_serialized, d_elements, num_elements);
    EXPECT_TRUE(check_batched_success(serialize_result));

    batched_error_code deserialize_result = deserialize_gt_vec(d_deserialized, d_serialized, num_elements);
    EXPECT_TRUE(check_batched_success(deserialize_result));

    CUDA_CHECK(cudaMemcpy(h_cuda_deserialized, d_deserialized, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_deserialized, h_cpu_deserialized, num_elements, gt_element_equal));

    // Test round-trip: deserialize(serialize(x)) = x
    EXPECT_TRUE(compare_results(h_cuda_deserialized, h_elements, num_elements, gt_element_equal))
        << "Round-trip serialization should preserve elements";

    // Cleanup
    CUDA_CHECK(cudaFree(d_elements));
    CUDA_CHECK(cudaFree(d_serialized));
    CUDA_CHECK(cudaFree(d_deserialized));
}

TEST(GTOpsCUDA, MultiPairing) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_pairs = 5;
    GT* d_miller_results, * d_final_results;
    G1AffinePoint* d_g1_points;
    G2AffinePoint* d_g2_points;

    G1AffinePoint h_g1_points[num_pairs];
    G2AffinePoint h_g2_points[num_pairs];
    GT h_cuda_miller[num_pairs], h_cpu_miller[num_pairs];
    GT h_cuda_final[num_pairs], h_cpu_final[num_pairs];

    // Create test point pairs
    G1Point g1_gen = G1Point::generator();
    G2Point g2_gen = G2Point::generator();
    InsecureRng rng;

    for (size_t i = 0; i < num_pairs; i++) {
        G1Point g1_mult;
        G2Point g2_mult;
        unsigned char scalar[32] = { 0 };
        scalar[0] = (rng.uniform_u64() % 254) + 1;
        scalar[1] = (rng.uniform_u64() % 255);

        G1Point::mult(g1_mult, g1_gen, scalar, 16);
        G2Point::mult(g2_mult, g2_gen, scalar, 16);

        G1Point::to_affine(h_g1_points[i], g1_mult);
        G2Point::to_affine(h_g2_points[i], g2_mult);
    }

    // CPU computation
    for (size_t i = 0; i < num_pairs; i++) {
        GT::miller_loop(h_cpu_miller[i], &h_g2_points[i], &h_g1_points[i]);
        GT::final_mult(h_cpu_final[i], h_cpu_miller[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_g1_points, num_pairs * sizeof(G1AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_g2_points, num_pairs * sizeof(G2AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_miller_results, num_pairs * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_final_results, num_pairs * sizeof(GT)));

    CUDA_CHECK(cudaMemcpy(d_g1_points, h_g1_points, num_pairs * sizeof(G1AffinePoint), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_g2_points, h_g2_points, num_pairs * sizeof(G2AffinePoint), cudaMemcpyHostToDevice));

    batched_error_code miller_result = miller_loop_gt_vec(d_miller_results, d_g2_points, d_g1_points, num_pairs);
    EXPECT_TRUE(check_batched_success(miller_result));

    batched_error_code final_result = final_exp_gt_vec(d_final_results, d_miller_results, num_pairs);
    EXPECT_TRUE(check_batched_success(final_result));

    CUDA_CHECK(cudaMemcpy(h_cuda_miller, d_miller_results, num_pairs * sizeof(GT), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_cuda_final, d_final_results, num_pairs * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_miller, h_cpu_miller, num_pairs, gt_element_equal));
    EXPECT_TRUE(compare_results(h_cuda_final, h_cpu_final, num_pairs, gt_element_equal));

    // Test that miller results are not in GT group but final results are
    for (size_t i = 0; i < num_pairs; i++) {
        EXPECT_FALSE(GT::in_group(h_cuda_miller[i]))
            << "Miller loop result should not be in GT group";
        EXPECT_TRUE(GT::in_group(h_cuda_final[i]))
            << "Final exponentiation result should be in GT group";
    }

    // Cleanup
    CUDA_CHECK(cudaFree(d_g1_points));
    CUDA_CHECK(cudaFree(d_g2_points));
    CUDA_CHECK(cudaFree(d_miller_results));
    CUDA_CHECK(cudaFree(d_final_results));
}

TEST(GTOpsCUDA, FinalVerify) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_pairs = 5;
    GT* d_gt1, * d_gt2;
    bool* d_results;
    GT h_gt1[num_pairs], h_gt2[num_pairs];
    bool h_cuda_results[num_pairs], h_cpu_results[num_pairs];

    // Create test data - pairs of GT elements
    InsecureRng rng;
    create_test_gt_elements(h_gt1, num_pairs, rng);
    create_test_gt_elements(h_gt2, num_pairs, rng);

    // Make some pairs equal for testing
    memcpy(&h_gt2[0], &h_gt1[0], sizeof(GT)); // Should return true
    memcpy(&h_gt2[2], &h_gt1[2], sizeof(GT)); // Should return true

    // CPU computation
    for (size_t i = 0; i < num_pairs; i++) {
        h_cpu_results[i] = blst_fp12_finalverify(h_gt1[i], h_gt2[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_gt1, num_pairs * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_gt2, num_pairs * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_results, num_pairs * sizeof(bool)));

    CUDA_CHECK(cudaMemcpy(d_gt1, h_gt1, num_pairs * sizeof(GT), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gt2, h_gt2, num_pairs * sizeof(GT), cudaMemcpyHostToDevice));

    batched_error_code result = finalverify_gt_vec(d_results, d_gt1, d_gt2, num_pairs);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_results, d_results, num_pairs * sizeof(bool), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_results, h_cpu_results, num_pairs, bool_equal));

    // Verify expected results
    EXPECT_TRUE(h_cuda_results[0]) << "Equal elements should verify as true";
    EXPECT_TRUE(h_cuda_results[2]) << "Equal elements should verify as true";

    // Cleanup
    CUDA_CHECK(cudaFree(d_gt1));
    CUDA_CHECK(cudaFree(d_gt2));
    CUDA_CHECK(cudaFree(d_results));
}

TEST(GTOpsCUDA, PowersOfTwo) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_bases = 16; // Multiple bases for parallel processing
    const size_t max_power = 8; // Test 2^0, 2^1, ..., 2^7
    const size_t exp_bits = 64;
    const size_t exp_bytes = (exp_bits + 7) / 8;

    GT* d_bases, * d_result;
    unsigned char* d_exponents;
    GT h_bases[num_bases];
    unsigned char h_exponents[num_bases * exp_bytes];
    GT h_cuda_result[num_bases], h_cpu_result[num_bases];

    // Create multiple base GT elements for parallel testing
    InsecureRng rng;
    create_test_gt_elements(h_bases, num_bases, rng);

    // Create exponents that are powers of 2, cycling through different powers
    for (size_t i = 0; i < num_bases; i++) {
        memset(&h_exponents[i * exp_bytes], 0, exp_bytes);
        size_t power = i % max_power; // Cycle through powers 0, 1, 2, ..., 7

        if (power == 0) {
            // Exponent = 0, result should be identity
            // Already zeroed
        }
        else {
            // Exponent = 2^(power-1), so power=1 gives exp=1, power=2 gives exp=2, etc.
            if (power <= sizeof(limb_t) * 8) {
                *((limb_t*)&h_exponents[i * exp_bytes]) = (limb_t)1 << (power - 1);
            }
        }
    }

    // CPU computation for comparison
    for (size_t i = 0; i < num_bases; i++) {
        GT::mult(h_cpu_result[i], h_bases[i], &h_exponents[i * exp_bytes], exp_bits);
    }

    // CUDA parallel computation
    CUDA_CHECK(cudaMalloc(&d_bases, num_bases * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_result, num_bases * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_exponents, num_bases * exp_bytes));

    CUDA_CHECK(cudaMemcpy(d_bases, h_bases, num_bases * sizeof(GT), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_exponents, h_exponents, num_bases * exp_bytes, cudaMemcpyHostToDevice));

    // Time the CUDA operation
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));

    CUDA_CHECK(cudaEventRecord(start));
    batched_error_code result = pow_gt_vec(d_result, d_bases, d_exponents, exp_bits, num_bases);
    CUDA_CHECK(cudaEventRecord(stop));

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float cuda_milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&cuda_milliseconds, start, stop));

    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_bases * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_bases, gt_element_equal));

    // Test specific properties for each power
	GT identity = GT::identity();

    for (size_t i = 0; i < num_bases; i++) {
        size_t power = i % max_power;
        limb_t actual_exponent = *((limb_t*)&h_exponents[i * exp_bytes]);

        switch (power) {
        case 0:
            // Exponent = 0 should give identity
            EXPECT_TRUE(GT::equal(h_cuda_result[i], identity))
                << "Base " << i << ": exponent 0 should give identity";
            break;
        case 1:
            // Exponent = 1 should give the base itself
            EXPECT_TRUE(GT::equal(h_cuda_result[i], h_bases[i]))
                << "Base " << i << ": exponent 1 should give base";
            break;
        case 2:
        {

            // Exponent = 2 should give base^2
            GT base_squared;
            GT::doubling(base_squared, h_bases[i]);
            EXPECT_TRUE(GT::equal(h_cuda_result[i], base_squared))
                << "Base " << i << ": exponent 2 should give base^2";
            break;
        }
        default:
            // All results should be in GT group
            EXPECT_TRUE(GT::in_group(h_cuda_result[i]))
                << "Base " << i << ": result should be in GT group (exponent = " << actual_exponent << ")";
            break;
        }
    }

    // Test mathematical properties across multiple bases
    // Verify that (base^(2^k))^2 = base^(2^(k+1)) for available pairs
    for (size_t i = 0; i < num_bases - 1; i++) {
        size_t power1 = i % max_power;
        size_t power2 = (i + 1) % max_power;

        // Only test when we have consecutive powers and both > 0
        if (power2 == power1 + 1 && power1 > 0 && power1 < 6) { // Limit to avoid overflow
            limb_t exp1 = *((limb_t*)&h_exponents[i * exp_bytes]);
            limb_t exp2 = *((limb_t*)&h_exponents[(i + 1) * exp_bytes]);

            // Verify exp2 = 2 * exp1
            if (exp2 == 2 * exp1) {
                GT squared_result;
                GT::doubling(squared_result, h_cuda_result[i]);

                // Find a base with the same element but the next exponent
                size_t next_idx = i + 1;
                if (next_idx < num_bases && GT::equal(h_bases[i], h_bases[next_idx])) {
                    EXPECT_TRUE(GT::equal(squared_result, h_cuda_result[next_idx]))
                        << "Mathematical property: (base^" << exp1 << ")^2 should equal base^" << exp2;
                }
            }
        }
    }

    // Test parallel consistency - all operations should complete successfully
    for (size_t i = 0; i < num_bases; i++) {
        EXPECT_TRUE(GT::in_group(h_cuda_result[i]))
            << "Parallel result " << i << " should be in GT group";
        EXPECT_FALSE(vec_is_zero(&h_cuda_result[i], sizeof(GT)))
            << "Parallel result " << i << " should not be zero";
    }

    // Enhanced debugging output
    //printf("Powers of Two Test Results:\n");
    //printf("  GPU Performance: %zu operations in %.3f ms (%.1f ops/ms)\n",
    //   num_bases, cuda_milliseconds, num_bases / cuda_milliseconds);

    for (size_t i = 0; i < std::min(num_bases, size_t(8)); i++) {
        size_t power = i % max_power;
        limb_t actual_exponent = *((limb_t*)&h_exponents[i * exp_bytes]);

        const char* status = "PASS";
        switch (power) {
        case 0:
            if (!GT::equal(h_cuda_result[i], identity)) status = "FAIL";
            break;
        case 1:
            if (!GT::equal(h_cuda_result[i], h_bases[i])) status = "FAIL";
            break;
        case 2: {
            GT base_squared;
            GT::doubling(base_squared, h_bases[i]);
            if (!GT::equal(h_cuda_result[i], base_squared)) status = "FAIL";
            break;
        }
        }

        if (status[0] != 'P') {
            printf("  Base %zu (exponent=%lu): %s\n", i, (unsigned long)actual_exponent, status);
        }
        //printf("  Base %zu (exponent=%lu): %s\n", i, (unsigned long)actual_exponent, status);
    }

    // Cleanup
    CUDA_CHECK(cudaFree(d_bases));
    CUDA_CHECK(cudaFree(d_result));
    CUDA_CHECK(cudaFree(d_exponents));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
}

// Performance and batch processing tests
TEST(GTOpsCUDA, PerformanceBaseline) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 100; // Reduced due to GT complexity
    GT* d_elements, * d_result, * h_elements;

    // Allocate host memory
    h_elements = new GT[num_elements];

    // Create test data
    InsecureRng rng;
    create_test_gt_elements(h_elements, num_elements, rng);

    // Allocate device memory
    CUDA_CHECK(cudaMalloc(&d_elements, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_result, num_elements * sizeof(GT)));

    // Copy to device
    CUDA_CHECK(cudaMemcpy(d_elements, h_elements, num_elements * sizeof(GT), cudaMemcpyHostToDevice));

    // Record timing
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));

    CUDA_CHECK(cudaEventRecord(start));
    batched_error_code result = square_gt_vec(d_result, d_elements, num_elements);
    CUDA_CHECK(cudaEventRecord(stop));

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));

    printf("CUDA GT Performance: %zu squaring operations in %.3f ms (%.1f ops/ms)\n",
        num_elements, milliseconds, num_elements / milliseconds);

    EXPECT_TRUE(check_batched_success(result));

    // Cleanup
    delete[] h_elements;
    CUDA_CHECK(cudaFree(d_elements));
    CUDA_CHECK(cudaFree(d_result));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
}

TEST(GTOpsCUDA, BatchProcessing) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t batch_size = 128; // Reduced for GT complexity
    GT* d_a, * d_b, * d_result, * h_a, * h_b;
    GT* h_cuda_result, * h_cpu_result;

    // Allocate host memory
    h_a = new GT[batch_size];
    h_b = new GT[batch_size];
    h_cuda_result = new GT[batch_size];
    h_cpu_result = new GT[batch_size];

    // Create test data
    InsecureRng rng;
    create_test_gt_elements(h_a, batch_size, rng);
    create_test_gt_elements(h_b, batch_size, rng);

    // CPU computation
    for (size_t i = 0; i < batch_size; i++) {
        GT::add(h_cpu_result[i], h_a[i], h_b[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_a, batch_size * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_b, batch_size * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_result, batch_size * sizeof(GT)));

    CUDA_CHECK(cudaMemcpy(d_a, h_a, batch_size * sizeof(GT), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b, batch_size * sizeof(GT), cudaMemcpyHostToDevice));

    batched_error_code mult_result = mult_gt_vec(d_result, d_a, d_b, batch_size);
    EXPECT_TRUE(check_batched_success(mult_result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, batch_size * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    bool all_match = compare_results(h_cuda_result, h_cpu_result, batch_size, gt_element_equal);

    printf("GT Batch processing: %s (CPU vs CUDA comparison)\n",
        all_match ? "All results match" : "Some results differ");

    EXPECT_TRUE(all_match);

    // Cleanup
    delete[] h_a;
    delete[] h_b;
    delete[] h_cuda_result;
    delete[] h_cpu_result;
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_result));
}

//// Advanced GT operation tests
//TEST(GTOpsCUDA, MultiExponentiation) {
//    if (!GTCUDAEnvironment::IsCudaAvailable()) {
//        GTEST_SKIP() << "CUDA not available";
//    }
//
//    const size_t num_bases = 5;
//    const size_t exp_bits = 32;
//    const size_t exp_bytes = (exp_bits + 7) / 8;
//
//    GT* d_bases, * d_result;
//    unsigned char* d_exponents;
//    GT h_bases[num_bases];
//    unsigned char h_exponents[num_bases * exp_bytes];
//    GT h_cuda_result, h_cpu_result;
//
//    // Create test data
//    InsecureRng rng;
//    create_test_gt_elements(h_bases, num_bases, rng);
//
//    // Create test exponents
//    for (size_t i = 0; i < num_bases; i++) {
//        memset(&h_exponents[i * exp_bytes], 0, exp_bytes);
//        h_exponents[i * exp_bytes] = (rng.uniform_u64() % 254) + 1;
//        h_exponents[i * exp_bytes + 1] = (rng.uniform_u64() % 255);
//    }
//
//    // CPU computation (simple multi-exponentiation)
//    GT::identity_copy(h_cpu_result);
//    for (size_t i = 0; i < num_bases; i++) {
//        GT temp;
//        GT::pow(temp, h_bases[i], (limb_t*)&h_exponents[i * exp_bytes], exp_bits);
//        if (i == 0) {
//            memcpy(h_cpu_result, temp, sizeof(GT));
//        }
//        else {
//            GT::add(h_cpu_result, h_cpu_result, temp);
//        }
//    }
//
//    // CUDA computation
//    CUDA_CHECK(cudaMalloc(&d_bases, num_bases * sizeof(GT)));
//    CUDA_CHECK(cudaMalloc(&d_result, sizeof(GT)));
//    CUDA_CHECK(cudaMalloc(&d_exponents, num_bases * exp_bytes));
//
//    CUDA_CHECK(cudaMemcpy(d_bases, h_bases, num_bases * sizeof(GT), cudaMemcpyHostToDevice));
//    CUDA_CHECK(cudaMemcpy(d_exponents, h_exponents, num_bases * exp_bytes, cudaMemcpyHostToDevice));
//
//    batched_error_code result = multi_exp_gt(d_result, d_bases, d_exponents, exp_bits, num_bases);
//    EXPECT_TRUE(check_batched_success(result));
//
//    CUDA_CHECK(cudaMemcpy(&h_cuda_result, d_result, sizeof(GT), cudaMemcpyDeviceToHost));
//
//    // Compare results
//    EXPECT_TRUE(gt_element_equal(&h_cuda_result, &h_cpu_result));
//
//    // Cleanup
//    CUDA_CHECK(cudaFree(d_bases));
//    CUDA_CHECK(cudaFree(d_result));
//    CUDA_CHECK(cudaFree(d_exponents));
//}

TEST(GTOpsCUDA, PairingVerification) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_pairs = 3;
    bool* d_results;
    G1AffinePoint* d_g1_points;
    G2AffinePoint* d_g2_points;
    GT* d_target_values;

    G1AffinePoint h_g1_points[num_pairs];
    G2AffinePoint h_g2_points[num_pairs];
    GT h_target_values[num_pairs];
    bool h_cuda_results[num_pairs], h_cpu_results[num_pairs];

    // Create test point pairs and compute target values
    G1Point g1_gen = G1Point::generator();
    G2Point g2_gen = G2Point::generator();
    InsecureRng rng;

    for (size_t i = 0; i < num_pairs; i++) {
        G1Point g1_mult;
        G2Point g2_mult;
        unsigned char scalar[32] = { 0 };
        scalar[0] = (rng.uniform_u64() % 254) + 1;
        scalar[1] = (rng.uniform_u64() % 255);

        G1Point::mult(g1_mult, g1_gen, scalar, 16);
        G2Point::mult(g2_mult, g2_gen, scalar, 16);

        G1Point::to_affine(h_g1_points[i], g1_mult);
        G2Point::to_affine(h_g2_points[i], g2_mult);

        // Compute target pairing value
        GT::miller_loop(h_target_values[i], &h_g2_points[i], &h_g1_points[i]);
        GT::final_mult(h_target_values[i], h_target_values[i]);
    }

    // CPU computation
    for (size_t i = 0; i < num_pairs; i++) {
        GT computed;
        GT::miller_loop(computed, &h_g2_points[i], &h_g1_points[i]);
        GT::final_mult(computed, computed);
        h_cpu_results[i] = GT::equal(computed, h_target_values[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_g1_points, num_pairs * sizeof(G1AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_g2_points, num_pairs * sizeof(G2AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_target_values, num_pairs * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_results, num_pairs * sizeof(bool)));

    CUDA_CHECK(cudaMemcpy(d_g1_points, h_g1_points, num_pairs * sizeof(G1AffinePoint), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_g2_points, h_g2_points, num_pairs * sizeof(G2AffinePoint), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_target_values, h_target_values, num_pairs * sizeof(GT), cudaMemcpyHostToDevice));

    batched_error_code result = pairing_verify_gt_vec(d_results, d_g2_points, d_g1_points, d_target_values, num_pairs);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_results, d_results, num_pairs * sizeof(bool), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_results, h_cpu_results, num_pairs, bool_equal));

    // All results should be true since we computed the target values correctly
    for (size_t i = 0; i < num_pairs; i++) {
        EXPECT_TRUE(h_cuda_results[i]) << "Pairing verification failed at index " << i;
    }

    // Cleanup
    CUDA_CHECK(cudaFree(d_g1_points));
    CUDA_CHECK(cudaFree(d_g2_points));
    CUDA_CHECK(cudaFree(d_target_values));
    CUDA_CHECK(cudaFree(d_results));
}

TEST(GTOpsCUDA, EqualityCheck) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 10;
    GT* d_a, * d_b;
    bool* d_results;
    GT h_a[num_elements], h_b[num_elements];
    bool h_cuda_results[num_elements], h_cpu_results[num_elements];

    // Create test data
    InsecureRng rng;
    create_test_gt_elements(h_a, num_elements, rng);
    create_test_gt_elements(h_b, num_elements, rng);

    // Make some elements equal for testing
    memcpy(&h_b[0], &h_a[0], sizeof(GT));
    memcpy(&h_b[3], &h_a[3], sizeof(GT));
    memcpy(&h_b[7], &h_a[7], sizeof(GT));

    // CPU computation
    for (size_t i = 0; i < num_elements; i++) {
        h_cpu_results[i] = GT::equal(h_a[i], h_b[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_a, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_b, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_results, num_elements * sizeof(bool)));

    CUDA_CHECK(cudaMemcpy(d_a, h_a, num_elements * sizeof(GT), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b, num_elements * sizeof(GT), cudaMemcpyHostToDevice));

    batched_error_code result = is_equal_gt_vec(d_results, d_a, d_b, num_elements);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_results, d_results, num_elements * sizeof(bool), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_results, h_cpu_results, num_elements, bool_equal));

    // Verify expected results
    EXPECT_TRUE(h_cuda_results[0]) << "Elements at index 0 should be equal";
    EXPECT_TRUE(h_cuda_results[3]) << "Elements at index 3 should be equal";
    EXPECT_TRUE(h_cuda_results[7]) << "Elements at index 7 should be equal";

    // Test reflexivity: element should equal itself
	GT identity = GT::identity();
    EXPECT_TRUE(GT::equal(identity, identity)) << "Identity should equal itself";

    // Cleanup
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_results));
}

TEST(GTOpsCUDA, AdvancedExponentiation) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 5;
    const size_t exp_bits = 128;
    const size_t exp_bytes = (exp_bits + 7) / 8;

    GT* d_bases, * d_result;
    unsigned char* d_exponents;
    GT h_bases[num_elements];
    unsigned char h_exponents[num_elements * exp_bytes];
    GT h_cuda_result[num_elements], h_cpu_result[num_elements];

    // Create test data
    InsecureRng rng;
    create_test_gt_elements(h_bases, num_elements, rng);

    // Create more complex exponents including edge cases
    for (size_t i = 0; i < num_elements; i++) {
        memset(&h_exponents[i * exp_bytes], 0, exp_bytes);

        switch (i) {
        case 0:
            // Test with 1
            h_exponents[i * exp_bytes] = 1;
            break;
        case 1:
            // Test with small prime
            h_exponents[i * exp_bytes] = 7;
            break;
        case 2:
            // Test with power of 2
            h_exponents[i * exp_bytes] = 16;
            break;
        case 3:
            // Test with larger value
            h_exponents[i * exp_bytes] = 255;
            h_exponents[i * exp_bytes + 1] = 1;
            break;
        case 4:
            // Test with high bit set
            h_exponents[i * exp_bytes + exp_bytes - 1] = 0x80;
            break;
        }
    }

    // CPU computation
    for (size_t i = 0; i < num_elements; i++) {
        GT::mult(h_cpu_result[i], h_bases[i], &h_exponents[i * exp_bytes], exp_bits);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_bases, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_result, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_exponents, num_elements * exp_bytes));

    CUDA_CHECK(cudaMemcpy(d_bases, h_bases, num_elements * sizeof(GT), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_exponents, h_exponents, num_elements * exp_bytes, cudaMemcpyHostToDevice));

    batched_error_code result = pow_gt_vec(d_result, d_bases, d_exponents, exp_bits, num_elements);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_elements, gt_element_equal));

    // Verify specific properties
    EXPECT_TRUE(GT::equal(h_cuda_result[0], h_bases[0]))
        << "base^1 should equal base";

    // Test that all results are in GT group
    for (size_t i = 0; i < num_elements; i++) {
        EXPECT_TRUE(GT::in_group(h_cuda_result[i]))
            << "Exponentiation result " << i << " should be in GT group";
    }

    // Cleanup
    CUDA_CHECK(cudaFree(d_bases));
    CUDA_CHECK(cudaFree(d_result));
    CUDA_CHECK(cudaFree(d_exponents));
}

TEST(GTOpsCUDA, Exponentiation256BitReducedScalars) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    constexpr size_t num_elements = 5;
    constexpr size_t exp_bits = 256;
    constexpr size_t exp_bytes = exp_bits / 8;

    GT* d_bases = nullptr;
    GT* d_result = nullptr;
    unsigned char* d_exponents = nullptr;
    GT h_bases[num_elements];
    unsigned char h_exponents[num_elements * exp_bytes] = {};
    GT h_cuda_result[num_elements], h_cpu_result[num_elements];

    InsecureRng rng;
    create_test_gt_elements(h_bases, num_elements, rng);

    const std::array<std::array<u64, 4>, num_elements> reduced_scalars = {
        std::array<u64, 4>{1, 0, 0, 0},
        scalar_from_modulus_minus({1, 0, 0, 0}),
        scalar_from_modulus_minus({0x12345, 0, 0, 0}),
        scalar_from_modulus_minus({0, 0x6789, 0, 0}),
        scalar_from_modulus_minus({0, 0, 0, 0x55})
    };

    for (size_t i = 0; i < num_elements; ++i) {
        memcpy(&h_exponents[i * exp_bytes], reduced_scalars[i].data(), exp_bytes);
        GT::mult(h_cpu_result[i], h_bases[i], &h_exponents[i * exp_bytes], exp_bits);
    }

    CUDA_CHECK(cudaMalloc(&d_bases, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_result, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_exponents, num_elements * exp_bytes));

    CUDA_CHECK(cudaMemcpy(d_bases, h_bases, num_elements * sizeof(GT), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_exponents, h_exponents, num_elements * exp_bytes, cudaMemcpyHostToDevice));

    batched_error_code result = pow_gt_vec(d_result, d_bases, d_exponents, exp_bits, num_elements);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));

    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_elements, gt_element_equal));

    CUDA_CHECK(cudaFree(d_bases));
    CUDA_CHECK(cudaFree(d_result));
    CUDA_CHECK(cudaFree(d_exponents));
}

TEST(GTOpsCUDA, SpecialElements) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 5;
    GT* d_elements;
    bool* d_results;
    GT h_elements[num_elements];
    bool h_cuda_results[num_elements], h_cpu_results[num_elements];

    // Create special test elements
    InsecureRng rng;
    create_test_gt_elements(h_elements, num_elements - 2, rng);

    // Add identity element
	h_elements[num_elements - 2] = GT::identity();

    // Add conjugate of first element
    memcpy(&h_elements[num_elements - 1], &h_elements[0], sizeof(GT));
    GT::conjugate(h_elements[num_elements - 1]);

    // CPU computation - test identity check
    for (size_t i = 0; i < num_elements; i++) {
        h_cpu_results[i] = GT::is_identity(h_elements[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_elements, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_results, num_elements * sizeof(bool)));

    CUDA_CHECK(cudaMemcpy(d_elements, h_elements, num_elements * sizeof(GT), cudaMemcpyHostToDevice));

    batched_error_code result = is_one_gt_vec(d_results, d_elements, num_elements);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_results, d_results, num_elements * sizeof(bool), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_results, h_cpu_results, num_elements, bool_equal));

    // Verify expected results
    EXPECT_FALSE(h_cuda_results[0]) << "Random GT element should not be identity";
    EXPECT_TRUE(h_cuda_results[num_elements - 2]) << "Identity element should be detected";
    EXPECT_FALSE(h_cuda_results[num_elements - 1]) << "Conjugate should not be identity (unless original was)";

    // Test that conjugate is in GT group
    EXPECT_TRUE(GT::in_group(h_elements[num_elements - 1]))
        << "Conjugate should be in GT group";

    // Cleanup
    CUDA_CHECK(cudaFree(d_elements));
    CUDA_CHECK(cudaFree(d_results));
}

TEST(GTOpsCUDA, PairingProperties) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_tests = 3;
    GT* d_results;
    G1AffinePoint* d_g1_points;
    G2AffinePoint* d_g2_points;

    G1AffinePoint h_g1_points[num_tests];
    G2AffinePoint h_g2_points[num_tests];
    GT h_cuda_results[num_tests], h_cpu_results[num_tests];

    // Setup test points for pairing properties
    G1Point g1_gen = G1Point::generator();
    G2Point g2_gen = G2Point::generator();

    // Test 1: e(G1_gen, G2_gen) - basic pairing
    G1Point::to_affine(h_g1_points[0], g1_gen);
    G2Point::to_affine(h_g2_points[0], g2_gen);

    // Test 2: e(2*G1_gen, G2_gen) 
    G1Point g1_double;
    G1Point::doubling(g1_double, g1_gen);
    G1Point::to_affine(h_g1_points[1], g1_double);
    G2Point::to_affine(h_g2_points[1], g2_gen);

    // Test 3: e(G1_gen, 2*G2_gen)
    G2Point g2_double;
    G2Point::doubling(g2_double, g2_gen);
    G1Point::to_affine(h_g1_points[2], g1_gen);
    G2Point::to_affine(h_g2_points[2], g2_double);

    // CPU computation
    for (size_t i = 0; i < num_tests; i++) {
        GT::miller_loop(h_cpu_results[i], &h_g2_points[i], &h_g1_points[i]);
        GT::final_mult(h_cpu_results[i], h_cpu_results[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_g1_points, num_tests * sizeof(G1AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_g2_points, num_tests * sizeof(G2AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_results, num_tests * sizeof(GT)));

    CUDA_CHECK(cudaMemcpy(d_g1_points, h_g1_points, num_tests * sizeof(G1AffinePoint), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_g2_points, h_g2_points, num_tests * sizeof(G2AffinePoint), cudaMemcpyHostToDevice));

    batched_error_code result = pairing_gt_vec(d_results, d_g2_points, d_g1_points, num_tests);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_results, d_results, num_tests * sizeof(GT), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_results, h_cpu_results, num_tests, gt_element_equal));

    // Test bilinearity: e(2P, Q) = e(P, 2Q) = e(P,Q)^2
    GT base_squared;
    GT::doubling(base_squared, h_cuda_results[0]);

    EXPECT_TRUE(GT::equal(h_cuda_results[1], h_cuda_results[2]))
        << "e(2P, Q) should equal e(P, 2Q)";
    EXPECT_TRUE(GT::equal(h_cuda_results[1], base_squared))
        << "e(2P, Q) should equal e(P,Q)^2";

    // All pairing results should be in GT group
    for (size_t i = 0; i < num_tests; i++) {
        EXPECT_TRUE(GT::in_group(h_cuda_results[i]))
            << "Pairing result " << i << " should be in GT group";
        EXPECT_FALSE(GT::is_identity(h_cuda_results[i]))
            << "Pairing result " << i << " should not be identity (non-degenerate)";
    }

    // Cleanup
    CUDA_CHECK(cudaFree(d_g1_points));
    CUDA_CHECK(cudaFree(d_g2_points));
    CUDA_CHECK(cudaFree(d_results));
}

TEST(GTOpsCUDA, ErrorHandling) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    // Test null pointer handling
    batched_error_code result;

    // Test null output pointer
    result = mult_gt_vec(nullptr, nullptr, nullptr, 1);
    EXPECT_EQ(result, BATCHED_ERROR_NULL_POINTER);

    // Test null input pointers
    GT* d_valid;
    CUDA_CHECK(cudaMalloc(&d_valid, sizeof(GT)));

    result = mult_gt_vec(d_valid, nullptr, d_valid, 1);
    EXPECT_EQ(result, BATCHED_ERROR_NULL_POINTER);

    result = mult_gt_vec(d_valid, d_valid, nullptr, 1);
    EXPECT_EQ(result, BATCHED_ERROR_NULL_POINTER);

    // Test zero size operations (should succeed)
    result = mult_gt_vec(d_valid, d_valid, d_valid, 0);
    EXPECT_EQ(result, BATCHED_SUCCESS);

    // Test other operations with null pointers
    result = square_gt_vec(nullptr, d_valid, 1);
    EXPECT_EQ(result, BATCHED_ERROR_NULL_POINTER);

    result = inverse_gt_vec(d_valid, nullptr, 1);
    EXPECT_EQ(result, BATCHED_ERROR_NULL_POINTER);

    bool* d_bool_result;
    CUDA_CHECK(cudaMalloc(&d_bool_result, sizeof(bool)));

    result = is_one_gt_vec(nullptr, d_valid, 1);
    EXPECT_EQ(result, BATCHED_ERROR_NULL_POINTER);

    result = is_one_gt_vec(d_bool_result, nullptr, 1);
    EXPECT_EQ(result, BATCHED_ERROR_NULL_POINTER);

    // Cleanup
    CUDA_CHECK(cudaFree(d_valid));
    CUDA_CHECK(cudaFree(d_bool_result));
}

TEST(GTOpsCUDA, LargeScaleOperations) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t large_batch_size = 500; // Test with larger batch
    GT* d_elements, * d_result, * h_elements, * h_cuda_result, * h_cpu_result;

    // Allocate host memory
    h_elements = new GT[large_batch_size];
    h_cuda_result = new GT[large_batch_size];
    h_cpu_result = new GT[large_batch_size];

    // Create test data
    InsecureRng rng;
    create_test_gt_elements(h_elements, large_batch_size, rng);

    // CPU computation (sample - full computation would be slow)
    const size_t sample_size = std::min(large_batch_size, size_t(20));
    for (size_t i = 0; i < sample_size; i++) {
        GT::doubling(h_cpu_result[i], h_elements[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_elements, large_batch_size * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_result, large_batch_size * sizeof(GT)));

    CUDA_CHECK(cudaMemcpy(d_elements, h_elements, large_batch_size * sizeof(GT), cudaMemcpyHostToDevice));

    // Time the operation
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));

    CUDA_CHECK(cudaEventRecord(start));
    batched_error_code result = square_gt_vec(d_result, d_elements, large_batch_size);
    CUDA_CHECK(cudaEventRecord(stop));

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));

    printf("Large scale GT squaring: %zu operations in %.3f ms (%.1f ops/ms)\n",
        large_batch_size, milliseconds, large_batch_size / milliseconds);

    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, large_batch_size * sizeof(GT), cudaMemcpyDeviceToHost));

    // Verify a sample of results
    bool sample_match = compare_results(h_cuda_result, h_cpu_result, sample_size, gt_element_equal);
    EXPECT_TRUE(sample_match) << "Sample of large batch results should match CPU computation";

    // Verify all results are in GT group
    for (size_t i = 0; i < std::min(large_batch_size, size_t(50)); i += 10) {
        EXPECT_TRUE(GT::in_group(h_cuda_result[i]))
            << "Large batch result " << i << " should be in GT group";
    }

    // Cleanup
    delete[] h_elements;
    delete[] h_cuda_result;
    delete[] h_cpu_result;
    CUDA_CHECK(cudaFree(d_elements));
    CUDA_CHECK(cudaFree(d_result));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
}

TEST(GTOpsCUDA, CyclotomicProperties) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_elements = 10;
    GT* d_elements, * d_regular_sqr, * d_cyclo_sqr;
    GT h_elements[num_elements];
    GT h_regular_result[num_elements], h_cyclo_result[num_elements];

    // Create test data (GT elements are in cyclotomic subgroup)
    InsecureRng rng;
    create_test_gt_elements(h_elements, num_elements, rng);

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_elements, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_regular_sqr, num_elements * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_cyclo_sqr, num_elements * sizeof(GT)));

    CUDA_CHECK(cudaMemcpy(d_elements, h_elements, num_elements * sizeof(GT), cudaMemcpyHostToDevice));

    // Compare regular vs cyclotomic squaring
    batched_error_code regular_result = square_gt_vec(d_regular_sqr, d_elements, num_elements);
    EXPECT_TRUE(check_batched_success(regular_result));

    batched_error_code cyclo_result = cyclotomic_square_gt_vec(d_cyclo_sqr, d_elements, num_elements);
    EXPECT_TRUE(check_batched_success(cyclo_result));

    CUDA_CHECK(cudaMemcpy(h_regular_result, d_regular_sqr, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_cyclo_result, d_cyclo_sqr, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));

    // For GT elements, both methods should give identical results
    EXPECT_TRUE(compare_results(h_regular_result, h_cyclo_result, num_elements, gt_element_equal))
        << "Regular and cyclotomic squaring should be identical for GT elements";

    // Test multiple squarings
    for (int rounds = 0; rounds < 3; rounds++) {
        batched_error_code round_result = cyclotomic_square_gt_vec(d_cyclo_sqr, d_cyclo_sqr, num_elements);
        EXPECT_TRUE(check_batched_success(round_result));
    }

    CUDA_CHECK(cudaMemcpy(h_cyclo_result, d_cyclo_sqr, num_elements * sizeof(GT), cudaMemcpyDeviceToHost));

    // Results should still be in GT group after multiple squarings
    for (size_t i = 0; i < num_elements; i++) {
        EXPECT_TRUE(GT::in_group(h_cyclo_result[i]))
            << "Multiple cyclotomic squaring result " << i << " should be in GT group";
    }

    // Cleanup
    CUDA_CHECK(cudaFree(d_elements));
    CUDA_CHECK(cudaFree(d_regular_sqr));
    CUDA_CHECK(cudaFree(d_cyclo_sqr));
}

TEST(GTOpsCUDA, ComprehensivePairingTest) {
    if (!GTCUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_pairs = 8;
    GT* d_miller_results, * d_final_results, * d_complete_results;
    G1AffinePoint* d_g1_points;
    G2AffinePoint* d_g2_points;
    bool* d_verify_results;

    G1AffinePoint h_g1_points[num_pairs];
    G2AffinePoint h_g2_points[num_pairs];
    GT h_miller_results[num_pairs], h_final_results[num_pairs], h_complete_results[num_pairs];
    bool h_verify_results[num_pairs];

    // Create diverse test point pairs
    G1Point g1_gen = G1Point::generator();
    G2Point g2_gen = G2Point::generator();
    InsecureRng rng;

    for (size_t i = 0; i < num_pairs; i++) {
        G1Point g1_mult;
        G2Point g2_mult;
        unsigned char scalar1[32] = { 0 }, scalar2[32] = { 0 };

        // Create different scalars for G1 and G2
        scalar1[0] = (rng.uniform_u64() % 254) + 1;
        scalar1[1] = (rng.uniform_u64() % 255);
        scalar2[0] = (rng.uniform_u64() % 254) + 1;
        scalar2[1] = (rng.uniform_u64() % 255);

        G1Point::mult(g1_mult, g1_gen, scalar1, 16);
        G2Point::mult(g2_mult, g2_gen, scalar2, 16);

        G1Point::to_affine(h_g1_points[i], g1_mult);
        G2Point::to_affine(h_g2_points[i], g2_mult);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_g1_points, num_pairs * sizeof(G1AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_g2_points, num_pairs * sizeof(G2AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_miller_results, num_pairs * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_final_results, num_pairs * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_complete_results, num_pairs * sizeof(GT)));
    CUDA_CHECK(cudaMalloc(&d_verify_results, num_pairs * sizeof(bool)));

    CUDA_CHECK(cudaMemcpy(d_g1_points, h_g1_points, num_pairs * sizeof(G1AffinePoint), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_g2_points, h_g2_points, num_pairs * sizeof(G2AffinePoint), cudaMemcpyHostToDevice));

    // Test separate Miller loop and final exponentiation
    batched_error_code miller_result = miller_loop_gt_vec(d_miller_results, d_g2_points, d_g1_points, num_pairs);
    EXPECT_TRUE(check_batched_success(miller_result));

    batched_error_code final_result = final_exp_gt_vec(d_final_results, d_miller_results, num_pairs);
    EXPECT_TRUE(check_batched_success(final_result));

    // Test complete pairing
    batched_error_code complete_result = pairing_gt_vec(d_complete_results, d_g2_points, d_g1_points, num_pairs);
    EXPECT_TRUE(check_batched_success(complete_result));

    // Test pairing verification (should all be true since we're verifying against themselves)
    batched_error_code verify_result = pairing_verify_gt_vec(d_verify_results, d_g2_points, d_g1_points, d_complete_results, num_pairs);
    EXPECT_TRUE(check_batched_success(verify_result));

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_miller_results, d_miller_results, num_pairs * sizeof(GT), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_final_results, d_final_results, num_pairs * sizeof(GT), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_complete_results, d_complete_results, num_pairs * sizeof(GT), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_verify_results, d_verify_results, num_pairs * sizeof(bool), cudaMemcpyDeviceToHost));

    // Verify consistency between separate and complete pairing
    EXPECT_TRUE(compare_results(h_final_results, h_complete_results, num_pairs, gt_element_equal))
        << "Separate Miller loop + final exp should equal complete pairing";

    // Verify all pairing verification results are true
    for (size_t i = 0; i < num_pairs; i++) {
        EXPECT_TRUE(h_verify_results[i]) << "Pairing verification " << i << " should be true";
    }

    // Verify group membership properties
    for (size_t i = 0; i < num_pairs; i++) {
        EXPECT_FALSE(GT::in_group(h_miller_results[i]))
            << "Miller loop result " << i << " should not be in GT group";
        EXPECT_TRUE(GT::in_group(h_final_results[i]))
            << "Final exponentiation result " << i << " should be in GT group";
        EXPECT_TRUE(GT::in_group(h_complete_results[i]))
            << "Complete pairing result " << i << " should be in GT group";
    }

    // Cleanup
    CUDA_CHECK(cudaFree(d_g1_points));
    CUDA_CHECK(cudaFree(d_g2_points));
    CUDA_CHECK(cudaFree(d_miller_results));
    CUDA_CHECK(cudaFree(d_final_results));
    CUDA_CHECK(cudaFree(d_complete_results));
    CUDA_CHECK(cudaFree(d_verify_results));
}
