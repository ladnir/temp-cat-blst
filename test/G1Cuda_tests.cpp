#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include "config.h"
#include "blst/curve/G1_cuda.h"
#include "blst/curve/map_to_G.h"
#include "blst/curve/multi_scalar_cuda.h"
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

// Helper function to create test points on host
void create_test_points(G1Point* points, size_t num, InsecureRng& rng) {
    G1Point generator = G1Point::generator();
    unsigned char scalar[32] = { 0 };

    for (size_t i = 0; i < num; i++) {
        scalar[0] = rng.uniform_u64();
        G1Point::mult(points[i], generator, scalar, 8);
    }
}

// Helper function to check batched_error_code
bool check_batched_success(batched_error_code code) {
    return code == BATCHED_SUCCESS;
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
bool point_equal(const G1Point& a, const G1Point& b) {
    return G1Point::equal(a, b);
}

bool affine_equal(const G1AffinePoint& a, const G1AffinePoint& b) {
    return G1AffinePoint::equal(a, b);
}

bool error_equal(const BLST_ERROR& a, const BLST_ERROR& b) {
    return a == b;
}

bool bool_equal(const bool& a, const bool& b) {
    return a == b;
}

// Global flag for CUDA availability
static bool cuda_device_available = false;

// CUDA Environment for runtime detection
class G1CUDAEnvironment : public ::testing::Environment {
public:
    void SetUp() override {
        cudaGetLastError();
        int device_count = 0;
        cudaError_t error = cudaGetDeviceCount(&device_count);

        if (error != cudaSuccess || device_count == 0) {
            std::cout << "No CUDA devices available. CUDA G1 tests will be skipped.\n";
            cuda_device_available = false;
        }
        else {
            std::cout << "Found " << device_count << " CUDA device(s) for G1 tests\n";
            cuda_device_available = true;
        }
    }

    static bool IsCudaAvailable() { return cuda_device_available; }
};

// Initialize the test environment
struct G1CudaTestRegistration {
    G1CudaTestRegistration() {
        ::testing::AddGlobalTestEnvironment(new G1CUDAEnvironment());
    }
};

static G1CudaTestRegistration g1_cuda_test_registration;




// CUDA test cases comparing against CPU versions
TEST(G1OpsCUDA, PointAddition) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_points = 10;
    G1Point* d_a, * d_b, * d_result;
    G1Point h_a[num_points], h_b[num_points];
    G1Point h_cuda_result[num_points], h_cpu_result[num_points];
    InsecureRng rng;

    // Create test data
    create_test_points(h_a, num_points, rng);
    create_test_points(h_b, num_points, rng);

    // CPU computation
    for (size_t i = 0; i < num_points; i++) {
        G1Point::add(h_cpu_result[i], h_a[i], h_b[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_a, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_b, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(G1Point)));

    CUDA_CHECK(cudaMemcpy(d_a, h_a, num_points * sizeof(G1Point), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b, num_points * sizeof(G1Point), cudaMemcpyHostToDevice));

    batched_error_code result = add_g1_vec(d_result, d_a, d_b, num_points);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_points * sizeof(G1Point), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results<G1Point>(h_cuda_result, h_cpu_result, num_points, point_equal));

    // Cleanup
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_result));
}

TEST(G1OpsCUDA, PointDoubling) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_points = 10;
    G1Point* d_points, * d_result;
    G1Point h_points[num_points];
    G1Point h_cuda_result[num_points], h_cpu_result[num_points];

    // Create test data
    InsecureRng rng;
    create_test_points(h_points, num_points, rng);

    // CPU computation
    for (size_t i = 0; i < num_points; i++) {
        G1Point::doubling(h_cpu_result[i], h_points[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(G1Point)));

    CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(G1Point), cudaMemcpyHostToDevice));

    batched_error_code result = double_g1_vec(d_result, d_points, num_points);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_points * sizeof(G1Point), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_points, point_equal));

    // Cleanup
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_result));
}

TEST(G1OpsCUDA, PointNegation) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_points = 10;
    G1Point* d_points, * d_result;
    G1Point h_points[num_points];
    G1Point h_cuda_result[num_points], h_cpu_result[num_points];
	int flags[num_points];
    int* d_flags;
    for (size_t i = 0; i < num_points; i++) {
        flags[i] = i&1; 
	}
    // Create test data
    InsecureRng rng;
    create_test_points(h_points, num_points, rng);

    // CPU computation
    for (size_t i = 0; i < num_points; i++) {
        h_cpu_result[i] = h_points[i];
        G1Point::cneg(h_cpu_result[i], flags[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_flags, num_points * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(G1Point), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_flags, flags, num_points * sizeof(int), cudaMemcpyHostToDevice));

    batched_error_code result = cneg_g1_vec(d_result, d_points, d_flags, num_points);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_points * sizeof(G1Point), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_points, point_equal));

    // Cleanup
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_result));
}

TEST(G1OpsCUDA, ScalarMultiplication) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_points = 10;
    const size_t scalar_bits = 128;
	const size_t scalar_bytes = (scalar_bits + 7) / 8;
    G1Point* d_points, * d_result;
    unsigned char* d_scalars;
    G1Point h_points[num_points];
    unsigned char h_scalars[num_points * scalar_bytes];
    G1Point h_cuda_result[num_points], h_cpu_result[num_points];

    // Create test data
    InsecureRng rng;
    create_test_points(h_points, num_points, rng);
    for (size_t i = 0; i < num_points; i++) {
        memset(&h_scalars[i * scalar_bytes], 0, scalar_bytes);
        for (size_t j = 0; j < scalar_bytes; j++) {
            h_scalars[i * scalar_bytes + j] = rng.uniform_u32(0, 255);
		}


    }

    // CPU computation
    for (size_t i = 0; i < num_points; i++) {
        G1Point::mult(h_cpu_result[i],h_points[i], &h_scalars[i * scalar_bytes], scalar_bits);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_scalars, num_points * (scalar_bits/8)));

    CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(G1Point), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_scalars, h_scalars, num_points * scalar_bytes, cudaMemcpyHostToDevice));

    batched_error_code result = mult_g1_vec(d_result, d_points, d_scalars, scalar_bits, num_points);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_points * sizeof(G1Point), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_points, point_equal));

    // Cleanup
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_result));
    CUDA_CHECK(cudaFree(d_scalars));
}

TEST(G1OpsCUDA, CurveMembership) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_points = 10;
    G1Point* d_points;
    bool* d_results;
    G1Point h_points[num_points];
    bool h_cuda_results[num_points], h_cpu_results[num_points];

    // Create test data (mix of valid and invalid points)
    InsecureRng rng;
    create_test_points(h_points, num_points - 2, rng);

    // Add some invalid points
    h_points[num_points - 2] = h_points[0];
    h_points[num_points - 2].X[0] ^= 1;  // Flip a bit to make invalid

    h_points[num_points - 1] = h_points[1];
    h_points[num_points - 1].Y[0] ^= 1;  // Flip a bit to make invalid

    // CPU computation
    for (size_t i = 0; i < num_points; i++) {
        h_cpu_results[i] = G1Point::on_curve(h_points[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_results, num_points * sizeof(bool)));

    CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(G1Point), cudaMemcpyHostToDevice));

    batched_error_code result = on_curve_g1_vec(d_results, d_points, num_points);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_results, d_results, num_points * sizeof(bool), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_results, h_cpu_results, num_points, bool_equal));

    // Cleanup
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_results));
}

TEST(G1OpsCUDA, Serialization) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_points = 10;
    G1Point* d_points;
    G1AffinePoint* d_affine_result;
    unsigned char* d_serialized;
    BLST_ERROR* d_errors;
    G1Point h_points[num_points];
    G1AffinePoint h_cuda_affine[num_points], h_cpu_affine[num_points];
    unsigned char h_serialized[num_points * 96];
    BLST_ERROR h_errors[num_points], h_errors_exp[num_points];


    // Create test data
    InsecureRng rng;
    create_test_points(h_points, num_points, rng);

    // CPU computation
    for (size_t i = 0; i < num_points; i++) {
        G1Point::serialize(&h_serialized[i * 96], h_points[i]);
        h_errors_exp[i] = G1AffinePoint::deserialize(h_cpu_affine[i], &h_serialized[i * 96]);
        EXPECT_EQ(h_errors_exp[i], BLST_SUCCESS);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_serialized, num_points * 96));
    CUDA_CHECK(cudaMalloc(&d_errors, num_points * 96));
    CUDA_CHECK(cudaMalloc(&d_affine_result, num_points * sizeof(G1AffinePoint)));

    CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(G1Point), cudaMemcpyHostToDevice));

    batched_error_code serialize_result = serialize_g1_vec(d_serialized, d_points, num_points);
    EXPECT_TRUE(check_batched_success(serialize_result));

    batched_error_code deserialize_result = deserialize_g1_vec(d_affine_result, d_serialized, d_errors, num_points);
    EXPECT_TRUE(check_batched_success(deserialize_result));

    CUDA_CHECK(cudaMemcpy(h_cuda_affine, d_affine_result, num_points * sizeof(G1AffinePoint), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_errors, d_errors, num_points * sizeof(BLST_ERROR), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_affine, h_cpu_affine, num_points, affine_equal));
    EXPECT_TRUE(compare_results(h_errors, h_errors_exp, num_points, error_equal));

    // Cleanup
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_serialized));
    CUDA_CHECK(cudaFree(d_affine_result));
}

TEST(G1OpsCUDA, Compression) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_points = 10;
    G1Point* d_points;
    G1AffinePoint* d_affine_result;
    unsigned char* d_compressed;
    BLST_ERROR* d_errors;
    G1Point h_points[num_points];
    G1AffinePoint h_points_affine[num_points];
    G1AffinePoint h_cuda_result[num_points], h_cpu_result[num_points];
    unsigned char h_compressed[num_points * 48];
    BLST_ERROR h_errors[num_points], h_errors_exp[num_points];

    // Create test data and convert to affine for CPU operations
    InsecureRng rng;
    create_test_points(h_points, num_points, rng);
    for (size_t i = 0; i < num_points; i++) {
        G1Point::to_affine(h_points_affine[i], h_points[i]);
    }

    // CPU computation
    for (size_t i = 0; i < num_points; i++) {
        G1AffinePoint::compress(&h_compressed[i * 48], h_points_affine[i]);
        h_errors_exp[i] = G1AffinePoint::uncompress(h_cpu_result[i], &h_compressed[i * 48]);
        EXPECT_EQ(h_errors_exp[i], BLST_SUCCESS);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_compressed, num_points * 48));
    CUDA_CHECK(cudaMalloc(&d_affine_result, num_points * sizeof(G1AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_errors, num_points * sizeof(BLST_ERROR)));

    CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(G1Point), cudaMemcpyHostToDevice));

    batched_error_code compress_result = compress_g1_vec(d_compressed, d_points, num_points);
    EXPECT_TRUE(check_batched_success(compress_result));

    batched_error_code uncompress_result = uncompress_g1_vec(d_affine_result, d_compressed, d_errors, num_points);
    EXPECT_TRUE(check_batched_success(uncompress_result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_affine_result, num_points * sizeof(G1AffinePoint), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_errors, d_errors, num_points * sizeof(BLST_ERROR), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_points, affine_equal));
    EXPECT_TRUE(compare_results(h_errors, h_errors_exp, num_points, error_equal));

    // Cleanup
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_compressed));
    CUDA_CHECK(cudaFree(d_affine_result));
    CUDA_CHECK(cudaFree(d_errors));
}

TEST(G1OpsCUDA, CoordinateConversion) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_points = 10;
    G1Point* d_points, * d_jacobian_result;
    G1AffinePoint* d_affine;
    G1Point h_points[num_points];
    G1AffinePoint h_cuda_affine[num_points], h_cpu_affine[num_points];
    G1Point h_cuda_jacobian[num_points], h_cpu_jacobian[num_points];

    // Create test data
    InsecureRng rng;
    create_test_points(h_points, num_points, rng);

    // CPU computation - Jacobian to affine
    for (size_t i = 0; i < num_points; i++) {
        G1Point::to_affine(h_cpu_affine[i], h_points[i]);
    }

    // CPU computation - affine back to Jacobian
    for (size_t i = 0; i < num_points; i++) {
        G1Point::from_affine(h_cpu_jacobian[i], h_cpu_affine[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_affine, num_points * sizeof(G1AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_jacobian_result, num_points * sizeof(G1Point)));

    CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(G1Point), cudaMemcpyHostToDevice));

    // Test Jacobian to affine conversion
    batched_error_code to_affine_result = to_affine_g1_vec(d_affine, d_points, num_points);
    EXPECT_TRUE(check_batched_success(to_affine_result));

    CUDA_CHECK(cudaMemcpy(h_cuda_affine, d_affine, num_points * sizeof(G1AffinePoint), cudaMemcpyDeviceToHost));

    // Test affine to Jacobian conversion
    batched_error_code from_affine_result = from_affine_g1_vec(d_jacobian_result, d_affine, num_points);
    EXPECT_TRUE(check_batched_success(from_affine_result));

    CUDA_CHECK(cudaMemcpy(h_cuda_jacobian, d_jacobian_result, num_points * sizeof(G1Point), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_affine, h_cpu_affine, num_points, affine_equal));
    EXPECT_TRUE(compare_results(h_cuda_jacobian, h_cpu_jacobian, num_points, point_equal));

    // Cleanup
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_affine));
    CUDA_CHECK(cudaFree(d_jacobian_result));
}

TEST(G1OpsCUDA, HashToCurve) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_msgs = 5;
    G1Point* d_result;
    unsigned char** d_msg_ptrs;
    unsigned char* d_dst;
    size_t* d_msg_lens;
    G1Point h_cuda_result[num_msgs], h_cpu_result[num_msgs];

    // Create test data
    std::vector<std::string> messages = {
        "message1", "message2", "message3", "message4", "message5"
    };
    std::string dst = "QUUX-V01-CS02-with-BLS12381G1";

    size_t h_msg_lens[num_msgs];
    unsigned char* h_msg_ptrs[num_msgs];
    unsigned char* d_msg_data[num_msgs];

    // Prepare message data
    for (size_t i = 0; i < num_msgs; i++) {
        h_msg_lens[i] = messages[i].length();

        // Allocate device memory for each message
        CUDA_CHECK(cudaMalloc(&d_msg_data[i], messages[i].length()));
        CUDA_CHECK(cudaMemcpy(d_msg_data[i], messages[i].c_str(), messages[i].length(), cudaMemcpyHostToDevice));

        h_msg_ptrs[i] = d_msg_data[i];
    }

    // CPU computation
    for (size_t i = 0; i < num_msgs; i++) {
        Hash_to_curve<G1Point>(h_cpu_result[i],
            reinterpret_cast<const unsigned char*>(messages[i].c_str()),
            messages[i].length(),
            reinterpret_cast<const unsigned char*>(dst.c_str()),
            dst.length(),
            NULL, 0);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_result, num_msgs * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_msg_ptrs, num_msgs * sizeof(unsigned char*)));
    CUDA_CHECK(cudaMalloc(&d_msg_lens, num_msgs * sizeof(size_t)));
    CUDA_CHECK(cudaMalloc(&d_dst, dst.length()));

    // Copy message pointer array and lengths to device
    CUDA_CHECK(cudaMemcpy(d_msg_ptrs, h_msg_ptrs, num_msgs * sizeof(unsigned char*), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_msg_lens, h_msg_lens, num_msgs * sizeof(size_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_dst, dst.c_str(), dst.length(), cudaMemcpyHostToDevice));

    batched_error_code result = hash_to_g1_vec(d_result, d_msg_ptrs, d_msg_lens, d_dst, dst.length(), num_msgs);
    EXPECT_TRUE(check_batched_success(result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_msgs * sizeof(G1Point), cudaMemcpyDeviceToHost));

    // Compare results
    EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_msgs, point_equal));

    // Cleanup
    for (size_t i = 0; i < num_msgs; i++) {
        CUDA_CHECK(cudaFree(d_msg_data[i]));
    }
    CUDA_CHECK(cudaFree(d_result));
    CUDA_CHECK(cudaFree(d_msg_ptrs));
    CUDA_CHECK(cudaFree(d_msg_lens));
    CUDA_CHECK(cudaFree(d_dst));
}

// Performance and batch processing tests
TEST(G1OpsCUDA, PerformanceBaseline) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t num_points = 1000;
    G1Point* d_points, * d_result, * h_points;
    unsigned char* d_scalars, * h_scalars;

    // Allocate host memory
    h_points = new G1Point[num_points];
    h_scalars = new unsigned char[num_points * 32];

    // Create test data
    InsecureRng rng;
    create_test_points(h_points, num_points, rng);
    for (size_t i = 0; i < num_points; i++) {
        memset(&h_scalars[i * 32], 0, 32);
        h_scalars[i * 32] = (i % 255) + 1;
    }

    // Allocate device memory
    CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_scalars, num_points * 32));

    // Copy to device
    CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(G1Point), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_scalars, h_scalars, num_points * 32, cudaMemcpyHostToDevice));

    // Record timing
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));

    CUDA_CHECK(cudaEventRecord(start));
    batched_error_code result = mult_g1_vec(d_result, d_points, d_scalars, 8, num_points);
    CUDA_CHECK(cudaEventRecord(stop));

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));

    printf("CUDA G1 Performance: %zu scalar multiplications in %.3f ms (%.1f ops/ms)\n",
        num_points, milliseconds, num_points / milliseconds);

    EXPECT_TRUE(check_batched_success(result));

    // Cleanup
    delete[] h_points;
    delete[] h_scalars;
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_result));
    CUDA_CHECK(cudaFree(d_scalars));
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
}

TEST(G1OpsCUDA, BatchProcessing) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    const size_t batch_size = 512;
    G1Point* d_points, * d_result, * h_points;
    G1Point* h_cuda_result, * h_cpu_result;

    // Allocate host memory
    h_points = new G1Point[batch_size];
    h_cuda_result = new G1Point[batch_size];
    h_cpu_result = new G1Point[batch_size];

    // Create test data
    InsecureRng rng;
    create_test_points(h_points, batch_size, rng);

    // CPU computation
    for (size_t i = 0; i < batch_size; i++) {
        G1Point::doubling(h_cpu_result[i], h_points[i]);
    }

    // CUDA computation
    CUDA_CHECK(cudaMalloc(&d_points, batch_size * sizeof(G1Point)));
    CUDA_CHECK(cudaMalloc(&d_result, batch_size * sizeof(G1Point)));

    CUDA_CHECK(cudaMemcpy(d_points, h_points, batch_size * sizeof(G1Point), cudaMemcpyHostToDevice));

    batched_error_code double_result = double_g1_vec(d_result, d_points, batch_size);
    EXPECT_TRUE(check_batched_success(double_result));

    CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, batch_size * sizeof(G1Point), cudaMemcpyDeviceToHost));

    // Compare results
    bool all_match = compare_results(h_cuda_result, h_cpu_result, batch_size, point_equal);

    printf("G1 Batch processing: %s (CPU vs CUDA comparison)\n",
        all_match ? "All results match" : "Some results differ");

    EXPECT_TRUE(all_match);

    // Cleanup
    delete[] h_points;
    delete[] h_cuda_result;
    delete[] h_cpu_result;
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_result));
}

TEST(G1OpsCUDA, SortDigitsParallel) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    // Test parameters
    const int N = 1000;  // Number of elements per digit
    const int D = 4;     // Number of digits  
    const int wbits = 5; // Bits per digit (max value = 31)
    const int W = 1 << wbits; // Number of possible values per digit

    // Allocate host memory - using u8 for keys and IdxSign for values
    u8* h_keys_input = new u8[N * D];
    u8* h_keys_output = new u8[N * D];
    blst::IdxSign* h_vals_input = new blst::IdxSign[N * D];
    blst::IdxSign* h_vals_output = new blst::IdxSign[N * D];
    u8* h_keys_expected = new u8[N * D];
    blst::IdxSign* h_vals_expected = new blst::IdxSign[N * D];

    // Generate test data - random keys and sequential values
    InsecureRng rng;
    for (int d = 0; d < D; d++) {
        for (int i = 0; i < N; i++) {
            int idx = d * N + i;
            // Generate random keys in valid range [0, 2^wbits)
            h_keys_input[idx] = static_cast<u8>(rng.uniform_u32(0, W - 1));
            // Use sequential values with IdxSign structure: (index << 1) | sign
            bool sign = (i % 3 == 0); // Every 3rd element is negative
            h_vals_input[idx] = blst::IdxSign{ static_cast<uint32_t>((i << 1) | (sign ? 1 : 0)) };
        }
    }

    // Create expected results by sorting each digit slice on CPU
    memcpy(h_keys_expected, h_keys_input, sizeof(u8) * N * D);
    memcpy(h_vals_expected, h_vals_input, sizeof(blst::IdxSign) * N * D);

    // CPU sort each digit slice separately
    for (int d = 0; d < D; d++) {
        u8* keys_slice = h_keys_expected + d * N;
        blst::IdxSign* vals_slice = h_vals_expected + d * N;

        // Simple bubble sort for reference (stable sort)
        for (int i = 0; i < N - 1; i++) {
            for (int j = 0; j < N - 1 - i; j++) {
                if (keys_slice[j] > keys_slice[j + 1]) {
                    // Swap keys
                    u8 temp_key = keys_slice[j];
                    keys_slice[j] = keys_slice[j + 1];
                    keys_slice[j + 1] = temp_key;

                    // Swap corresponding values
                    blst::IdxSign temp_val = vals_slice[j];
                    vals_slice[j] = vals_slice[j + 1];
                    vals_slice[j + 1] = temp_val;
                }
            }
        }
    }

    // Allocate device memory
    u8* d_keys_input;
    u8* d_keys_output;
    blst::IdxSign* d_vals_input;
    blst::IdxSign* d_vals_output;
    CUDA_CHECK(cudaMalloc(&d_keys_input, sizeof(u8) * N * D));
    CUDA_CHECK(cudaMalloc(&d_keys_output, sizeof(u8) * N * D));
    CUDA_CHECK(cudaMalloc(&d_vals_input, sizeof(blst::IdxSign) * N * D));
    CUDA_CHECK(cudaMalloc(&d_vals_output, sizeof(blst::IdxSign) * N * D));

    // Copy test data to device
    CUDA_CHECK(cudaMemcpy(d_keys_input, h_keys_input, sizeof(u8) * N * D, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vals_input, h_vals_input, sizeof(blst::IdxSign) * N * D, cudaMemcpyHostToDevice));

    // Query workspace size
    size_t workspace_bytes = 0;
    cudaError_t result = blst::sort_digits_segmented(
        nullptr,          // d_workspace (nullptr to query size)
        workspace_bytes,  // in/out workspace size
        d_keys_input,
        d_keys_output,
        d_vals_input,
        d_vals_output,
        N,
        D,
        wbits,
        0);               // default stream

    EXPECT_EQ(result, cudaSuccess);
    EXPECT_GT(workspace_bytes, 0) << "Workspace size should be positive";

    // Allocate workspace
    void* d_workspace = nullptr;
    CUDA_CHECK(cudaMalloc(&d_workspace, workspace_bytes));

    // Call the function under test
    result = blst::sort_digits_segmented(
        d_workspace,
        workspace_bytes,
        d_keys_input,
        d_keys_output,
        d_vals_input,
        d_vals_output,
        N,
        D,
        wbits,
        0);

    EXPECT_EQ(result, cudaSuccess);

    // Copy results back to host
    CUDA_CHECK(cudaMemcpy(h_keys_output, d_keys_output, sizeof(u8) * N * D, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vals_output, d_vals_output, sizeof(blst::IdxSign) * N * D, cudaMemcpyDeviceToHost));

    // Verify results - each digit slice should be sorted
    bool all_sorted = true;
    for (int d = 0; d < D; d++) {
        u8* keys_slice = h_keys_output + d * N;
        blst::IdxSign* vals_slice = h_vals_output + d * N;

        // Check that keys are sorted within each digit
        for (int i = 0; i < N - 1; i++) {
            if (keys_slice[i] > keys_slice[i + 1]) {
                printf("Digit %d not sorted at position %d: %u > %u\n",
                    d, i, keys_slice[i], keys_slice[i + 1]);
                all_sorted = false;
                break;
            }
        }

        // Verify that values follow their corresponding keys correctly
        for (int i = 0; i < N - 1; i++) {
            if (keys_slice[i] == keys_slice[i + 1]) {
                // For equal keys, check that indices are valid
                uint32_t idx1 = vals_slice[i].idx();
                uint32_t idx2 = vals_slice[i + 1].idx();
                if (idx1 >= N || idx2 >= N) {
                    printf("Invalid index found in digit %d\n", d);
                    all_sorted = false;
                    break;
                }
            }
        }
    }

    EXPECT_TRUE(all_sorted) << "Sort verification failed";

    // Additional check: verify all original indices are present
    for (int d = 0; d < D; d++) {
        std::vector<bool> found(N, false);
        blst::IdxSign* vals_slice = h_vals_output + d * N;

        for (int i = 0; i < N; i++) {
            uint32_t idx = vals_slice[i].idx();
            if (idx < N) {
                found[idx] = true;
            }
        }

        bool all_found = true;
        for (int i = 0; i < N; i++) {
            if (!found[i]) {
                printf("Index %d missing in digit %d\n", i, d);
                all_found = false;
            }
        }
        EXPECT_TRUE(all_found) << "Some indices were lost during sorting in digit " << d;
    }

    // Cleanup
    delete[] h_keys_input;
    delete[] h_keys_output;
    delete[] h_vals_input;
    delete[] h_vals_output;
    delete[] h_keys_expected;
    delete[] h_vals_expected;

    CUDA_CHECK(cudaFree(d_keys_input));
    CUDA_CHECK(cudaFree(d_keys_output));
    CUDA_CHECK(cudaFree(d_vals_input));
    CUDA_CHECK(cudaFree(d_vals_output));
    CUDA_CHECK(cudaFree(d_workspace));
}

// Additional edge case test
TEST(G1OpsCUDA, SortDigitsParallelEdgeCases) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    // Test with minimal data
    const int N_small = 10;
    const int D_small = 2;
    const int wbits_small = 3; // 8 possible values

    u8* d_keys_input;
    u8* d_keys_output;
    blst::IdxSign* d_vals_input;
    blst::IdxSign* d_vals_output;
    CUDA_CHECK(cudaMalloc(&d_keys_input, sizeof(u8) * N_small * D_small));
    CUDA_CHECK(cudaMalloc(&d_keys_output, sizeof(u8) * N_small * D_small));
    CUDA_CHECK(cudaMalloc(&d_vals_input, sizeof(blst::IdxSign) * N_small * D_small));
    CUDA_CHECK(cudaMalloc(&d_vals_output, sizeof(blst::IdxSign) * N_small * D_small));

    // Test data: already sorted
    u8 h_sorted_keys[] = { 0, 1, 2, 3, 4, 5, 6, 7, 0, 1,  // digit 0
                          0, 0, 1, 1, 2, 2, 3, 3, 4, 4 }; // digit 1
    blst::IdxSign h_sorted_vals[N_small * D_small];
    for (int i = 0; i < N_small * D_small; i++) {
        h_sorted_vals[i] = blst::IdxSign{ static_cast<uint32_t>(i << 1) }; // No sign bit
    }

    CUDA_CHECK(cudaMemcpy(d_keys_input, h_sorted_keys,
        sizeof(u8) * N_small * D_small, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vals_input, h_sorted_vals,
        sizeof(blst::IdxSign) * N_small * D_small, cudaMemcpyHostToDevice));

    // Query workspace size
    size_t workspace_bytes = 0;
    cudaError_t result = blst::sort_digits_segmented(
        nullptr,
        workspace_bytes,
        d_keys_input,
        d_keys_output,
        d_vals_input,
        d_vals_output,
        N_small,
        D_small,
        wbits_small,
        0);

    EXPECT_EQ(result, cudaSuccess);

    // Allocate workspace
    void* d_workspace = nullptr;
    CUDA_CHECK(cudaMalloc(&d_workspace, workspace_bytes));

    // Should handle already sorted data gracefully
    result = blst::sort_digits_segmented(
        d_workspace,
        workspace_bytes,
        d_keys_input,
        d_keys_output,
        d_vals_input,
        d_vals_output,
        N_small,
        D_small,
        wbits_small,
        0);

    EXPECT_EQ(result, cudaSuccess);

    // Test with reverse sorted data
    u8 h_reverse_keys[] = { 7, 6, 5, 4, 3, 2, 1, 0, 1, 0,  // digit 0
                           4, 4, 3, 3, 2, 2, 1, 1, 0, 0 }; // digit 1

    CUDA_CHECK(cudaMemcpy(d_keys_input, h_reverse_keys,
        sizeof(u8) * N_small * D_small, cudaMemcpyHostToDevice));

    result = blst::sort_digits_segmented(
        d_workspace,
        workspace_bytes,
        d_keys_input,
        d_keys_output,
        d_vals_input,
        d_vals_output,
        N_small,
        D_small,
        wbits_small,
        0);

    EXPECT_EQ(result, cudaSuccess);

    // Verify the reverse sorted data is now properly sorted
    u8 h_result_keys[N_small * D_small];
    CUDA_CHECK(cudaMemcpy(h_result_keys, d_keys_output,
        sizeof(u8) * N_small * D_small, cudaMemcpyDeviceToHost));

    // Check first digit is sorted
    for (int i = 0; i < N_small - 1; i++) {
        EXPECT_LE(h_result_keys[i], h_result_keys[i + 1])
            << "First digit not sorted at position " << i;
    }

    // Check second digit is sorted  
    for (int i = N_small; i < 2 * N_small - 1; i++) {
        EXPECT_LE(h_result_keys[i], h_result_keys[i + 1])
            << "Second digit not sorted at position " << (i - N_small);
    }

    CUDA_CHECK(cudaFree(d_keys_input));
    CUDA_CHECK(cudaFree(d_keys_output));
    CUDA_CHECK(cudaFree(d_vals_input));
    CUDA_CHECK(cudaFree(d_vals_output));
    CUDA_CHECK(cudaFree(d_workspace));
}

TEST(G1OpsCUDA, ReduceGroups) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    // Test parameters
    const int M = 4;         // number of MSMs (reduced for easier testing)
    const int N = 100;       // number of elements per MSM
    const int wbits = 4;     // Window size (4 bits = 16 possible values)
    const int W = (1 << (wbits - 1)); // Number of groups = 2^(wbits-1) = 8

    // Host data - now arrays sized for M batches
    std::vector<u8> h_keys(M * N);
    std::vector<blst::IdxSign> h_vals(M * N);
    std::vector<G1Point> h_points(N);  // Shared points across all MSMs
    std::vector<Pointxyzz<Fp>> h_cuda_groups_xyzz(M * W);
    std::vector<Pointxyzz<Fp>> h_cpu_groups_xyzz(M * W);
    std::vector<G1Point> h_cuda_groups(M * W), h_cpu_groups(M * W);
    std::vector<G1AffinePoint> h_cuda_groups_affine(M * W), h_cpu_groups_affine(M * W);

    // Initialize groups to zero
    std::fill(h_cuda_groups_xyzz.begin(), h_cuda_groups_xyzz.end(), Pointxyzz<Fp>{});
    std::fill(h_cpu_groups_xyzz.begin(), h_cpu_groups_xyzz.end(), Pointxyzz<Fp>{});
    std::fill(h_cuda_groups.begin(), h_cuda_groups.end(), G1Point{});
    std::fill(h_cpu_groups.begin(), h_cpu_groups.end(), G1Point{});
    std::fill(h_cuda_groups_affine.begin(), h_cuda_groups_affine.end(), G1AffinePoint{});
    std::fill(h_cpu_groups_affine.begin(), h_cpu_groups_affine.end(), G1AffinePoint{});

    // Create test data
    InsecureRng rng;
    create_test_points(h_points.data(), N, rng);

    // Generate sorted keys and packed values for each MSM
    for (int m = 0; m < M; m++) {
        for (int i = 0; i < N; i++) {
            int global_idx = m * N + i;

            // Keys range from 1 to W (group 0 is skipped)
            // Use different patterns per MSM for variety
            h_keys[global_idx] = ((i + m) % W) + 1;  // Groups 1 to W

            // Pack values: (index << 1) | sign using IdxSign
            // The index refers to the point in the shared points array
            bool sign = ((i + m) % 3 == 0); // Every 3rd element is negative (different per MSM)
            h_vals[global_idx] = blst::IdxSign{ static_cast<uint32_t>((i << 1) | (sign ? 1 : 0)) };
        }

        // Sort by keys within each MSM (bubble sort for test data)
        for (int i = 0; i < N - 1; i++) {
            for (int j = 0; j < N - 1 - i; j++) {
                int idx1 = m * N + j;
                int idx2 = m * N + j + 1;
                if (h_keys[idx1] > h_keys[idx2]) {
                    // Swap keys
                    std::swap(h_keys[idx1], h_keys[idx2]);
                    // Swap values
                    std::swap(h_vals[idx1], h_vals[idx2]);
                }
            }
        }
    }

    // Convert input points to affine for CUDA
    std::vector<G1AffinePoint> h_affine_points(N);
    for (int i = 0; i < N; i++) {
        G1Point::to_affine(h_affine_points[i], h_points[i]);
    }

    // CPU computation - compute expected group sums for each MSM using XYZZ accumulation
    for (int m = 0; m < M; m++) {
        for (int i = 0; i < N; i++) {
            int global_idx = m * N + i;
            u8 group = h_keys[global_idx];
            if (group == 0 || group > W) {
                continue; // Skip invalid groups
            }

            uint32_t idx = h_vals[global_idx].idx();
            bool sign = h_vals[global_idx].sign();

            if (idx >= N) {
                continue; // Skip invalid indices
            }

            G1AffinePoint point_affine = h_affine_points[idx];
            if (sign) {
                G1AffinePoint::cneg(point_affine, sign);
            }

            Pointxyzz<Fp> point_xyzz;
            Pointxyzz<Fp>::from_affine(point_xyzz, point_affine);

            int group_idx = m * W + (group - 1);
            Pointxyzz<Fp> sum;
            Pointxyzz<Fp>::dadd(sum, h_cpu_groups_xyzz[group_idx], point_xyzz);
            h_cpu_groups_xyzz[group_idx] = sum;
        }
    }

    // Convert CPU results to Jacobian and affine for comparison
    for (int m = 0; m < M; m++) {
        for (int g = 0; g < W; g++) {
            int group_idx = m * W + g;
            Pointxyzz<Fp>::to_Jacobian(h_cpu_groups[group_idx], h_cpu_groups_xyzz[group_idx]);
            G1Point::to_affine(h_cpu_groups_affine[group_idx], h_cpu_groups[group_idx]);
        }
    }

    // Device memory allocation
    u8* d_keys;
    G1AffinePoint* d_points;
    Pointxyzz<Fp>* d_groups;
    blst::IdxSign* d_vals;
    void* d_workspace = nullptr;

    CUDA_CHECK(cudaMalloc(&d_keys, M * N * sizeof(u8)));
    CUDA_CHECK(cudaMalloc(&d_points, N * sizeof(G1AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_groups, M * W * sizeof(Pointxyzz<Fp>)));
    CUDA_CHECK(cudaMalloc(&d_vals, M * N * sizeof(blst::IdxSign)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_keys, h_keys.data(), M * N * sizeof(u8), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_points, h_affine_points.data(), N * sizeof(G1AffinePoint), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vals, h_vals.data(), M * N * sizeof(blst::IdxSign), cudaMemcpyHostToDevice));

    // Create transform iterator for gathering points with sign
    auto device_gather_iter = thrust::make_transform_iterator(
        d_vals,
        blst::PackedToPointXYZZ<Fp>{ d_points, N }
    );

    // Query workspace size - using batched version (M MSMs, N elements each, W groups)
    size_t workspace_bytes = 0;
    cudaError_t result = blst::reduce_groups<blst::Fp>(
        nullptr,          // d_workspace (nullptr to query size)
        workspace_bytes,  // in/out workspace size
        d_keys,
        device_gather_iter,
        d_groups,
        M,                // Number of MSMs
        N,                // Number of elements per MSM
        W,                // Number of groups per MSM
        0);               // default stream

    EXPECT_EQ(result, cudaSuccess);
    EXPECT_GT(workspace_bytes, 0) << "Workspace size should be positive";

    // Allocate workspace
    CUDA_CHECK(cudaMalloc(&d_workspace, workspace_bytes));

    // Perform the actual batched reduction
    result = blst::reduce_groups<blst::Fp>(
        d_workspace,
        workspace_bytes,
        d_keys,
        device_gather_iter,
        d_groups,
        M,                // Number of MSMs
        N,                // Number of elements per MSM
        W,                // Number of groups per MSM
        0);

    // Block on stream to catch errors
    CUDA_CHECK(cudaStreamSynchronize(0));

    if (result != cudaSuccess) {
        const char* err_str = cudaGetErrorString(result);
        printf("batched reduce_groups failed: %s\n", err_str);
    }

    EXPECT_EQ(result, cudaSuccess);

    // Copy results back to host
    CUDA_CHECK(cudaMemcpy(h_cuda_groups_xyzz.data(), d_groups, M * W * sizeof(Pointxyzz<Fp>), cudaMemcpyDeviceToHost));

    // Convert CUDA results to Jacobian and affine for comparison
    for (int m = 0; m < M; m++) {
        for (int g = 0; g < W; g++) {
            int group_idx = m * W + g;
            Pointxyzz<Fp>::to_Jacobian(h_cuda_groups[group_idx], h_cuda_groups_xyzz[group_idx]);
            G1Point::to_affine(h_cuda_groups_affine[group_idx], h_cuda_groups[group_idx]);
        }
    }

    // Compare results for each MSM
    bool all_match = true;
    for (int m = 0; m < M; m++) {
        for (int g = 0; g < W; g++) {
            int group_idx = m * W + g;
            if (!G1AffinePoint::equal(h_cuda_groups_affine[group_idx], h_cpu_groups_affine[group_idx])) {
                printf("MSM %d, Group %d mismatch\n", m, g + 1);
                all_match = false;
            }
        }
    }

    EXPECT_TRUE(all_match) << "CUDA and CPU batched group reductions should match";

    // Additional verification: check that some groups have non-zero results
    bool has_nonzero = false;
    for (int m = 0; m < M; m++) {
        for (int g = 0; g < W; g++) {
            int group_idx = m * W + g;
            if (!G1AffinePoint::is_inf(h_cuda_groups_affine[group_idx])) {
                has_nonzero = true;
                break;
            }
        }
        if (has_nonzero) break;
    }
    EXPECT_TRUE(has_nonzero) << "At least some groups should have non-zero results";

    //printf("Batched reduce_groups test: %d MSMs, %d points each, %d groups per MSM - %s\n",
    //    M, N, W, all_match ? "PASSED" : "FAILED");

    // Cleanup
    CUDA_CHECK(cudaFree(d_keys));
    CUDA_CHECK(cudaFree(d_vals));
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_groups));
    CUDA_CHECK(cudaFree(d_workspace));
}
TEST(G1OpsCUDA, ReduceGroupsEdgeCases) {
    if (!G1CUDAEnvironment::IsCudaAvailable()) {
        GTEST_SKIP() << "CUDA not available";
    }

    // Test with minimal data and edge cases
    const int M = 2;         // number of MSMs (small for edge case testing)
    const int N = 10;        // number of elements per MSM
    const int wbits = 3;     // 3 bits = 8 possible values
    const int W = (1 << (wbits - 1)); // 4 groups

    // Properly sorted keys: all instances of same key are contiguous within each MSM
    std::vector<u8> h_keys(M * N);
    std::vector<blst::IdxSign> h_vals(M * N);
    std::vector<G1Point> h_points(N);
    std::vector<G1AffinePoint> h_affine_points(N);

    // Create test points
    InsecureRng rng;
    create_test_points(h_points.data(), N, rng);
    for (int i = 0; i < N; i++) {
        G1Point::to_affine(h_affine_points[i], h_points[i]);
    }

    // Set up keys and values for each MSM
    for (int m = 0; m < M; m++) {
        // Create sorted keys pattern for this MSM
        std::vector<u8> msm_keys = { 1, 1, 1, 2, 2, 2, 3, 3, 4, 4 }; // Actually sorted

        for (int i = 0; i < N; i++) {
            int global_idx = m * N + i;
            h_keys[global_idx] = msm_keys[i];
            h_vals[global_idx] = blst::IdxSign{ static_cast<uint32_t>(i << 1) }; // All positive (sign bit = 0)
        }
    }

    // Device allocation
    u8* d_keys;
    blst::IdxSign* d_vals;
    G1AffinePoint* d_points;
    Pointxyzz<Fp>* d_groups;

    CUDA_CHECK(cudaMalloc(&d_keys, M * N * sizeof(u8)));
    CUDA_CHECK(cudaMalloc(&d_vals, M * N * sizeof(blst::IdxSign)));
    CUDA_CHECK(cudaMalloc(&d_points, N * sizeof(G1AffinePoint)));
    CUDA_CHECK(cudaMalloc(&d_groups, M * W * sizeof(Pointxyzz<Fp>)));

    CUDA_CHECK(cudaMemcpy(d_keys, h_keys.data(), M * N * sizeof(u8), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vals, h_vals.data(), M * N * sizeof(blst::IdxSign), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_points, h_affine_points.data(), N * sizeof(G1AffinePoint), cudaMemcpyHostToDevice));

    // Create transform iterator
    auto gather_iter = thrust::make_transform_iterator(
        d_vals,
        blst::PackedToPointXYZZ<Fp>{ d_points, N }
    );

    // Test workspace size query for batched version
    size_t workspace_bytes = 0;
    cudaError_t result = blst::reduce_groups<blst::Fp>(
        nullptr, workspace_bytes, d_keys, gather_iter, d_groups, M, N, W, 0);

    EXPECT_EQ(result, cudaSuccess);

    // Test with insufficient workspace
    void* d_workspace;
    size_t small_workspace = workspace_bytes / 2;
    CUDA_CHECK(cudaMalloc(&d_workspace, small_workspace));

    result = blst::reduce_groups<blst::Fp>(
        d_workspace, small_workspace, d_keys, gather_iter, d_groups, M, N, W, 0);

    EXPECT_EQ(result, cudaErrorInvalidValue) << "Should fail with insufficient workspace";

    // Test with correct workspace size
    CUDA_CHECK(cudaFree(d_workspace));
    CUDA_CHECK(cudaMalloc(&d_workspace, workspace_bytes));

    result = blst::reduce_groups<blst::Fp>(
        d_workspace, workspace_bytes, d_keys, gather_iter, d_groups, M, N, W, 0);

    EXPECT_EQ(result, cudaSuccess);

    //printf("Batched reduce_groups edge cases: %d MSMs, %d points each, %d groups per MSM - PASSED\n", M, N, W);

    // Cleanup
    CUDA_CHECK(cudaFree(d_keys));
    CUDA_CHECK(cudaFree(d_vals));
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_groups));
    CUDA_CHECK(cudaFree(d_workspace));
}



