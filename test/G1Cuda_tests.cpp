#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include "config.h"
#include "blst/curve/G1_cuda.h"
#include "blst/curve/map_to_G.h"
//#include "blst/curve/multi_scalar_par.h"
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
