#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include "config.h"
#include "blst/curve/G2_cuda.h"
#include "blst/curve/e2.h"
#include "blst/curve/map_to_g2.h"
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

	// Helper function to create test points on host
	void create_test_points_g2(POINTonE2* points, size_t num, InsecureRng& rng) {
		POINTonE2 generator = *blst_p2_generator();
		unsigned char scalar[32] = { 0 };

		for (size_t i = 0; i < num; i++) {
			scalar[0] = rng.uniform_u64() & 0xFF;
			if (scalar[0] == 0) scalar[0] = 1; // Avoid zero scalar
			blst_p2_mult(&points[i], &generator, scalar, 8);
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
		bool (*equality_func)(const T*, const T*)) {
		for (size_t i = 0; i < num; i++) {
			if (!equality_func(&cuda_results[i], &cpu_results[i])) {
				printf("Mismatch at index %zu\n", i);
				return false;
			}
		}
		return true;
	}

	// Equality functions for different types
	bool point_equal_g2(const POINTonE2* a, const POINTonE2* b) {
		return blst_p2_is_equal(a, b);
	}

	bool affine_equal_g2(const POINTonE2_affine* a, const POINTonE2_affine* b) {
		return blst_p2_affine_is_equal(a, b);
	}

	bool error_equal(const BLST_ERROR* a, const BLST_ERROR* b) {
		return *a == *b;
	}

	bool bool_equal(const bool* a, const bool* b) {
		return *a == *b;
	}

	// Global flag for CUDA availability
	static bool cuda_device_available = false;

}

// CUDA Environment for runtime detection
class G2CUDAEnvironment : public ::testing::Environment {
public:
	void SetUp() override {
		int device_count = 0;
		cudaError_t error = cudaGetDeviceCount(&device_count);

		if (error != cudaSuccess || device_count == 0) {
			std::cout << "No CUDA devices available. CUDA G2 tests will be skipped.\n";
			cuda_device_available = false;
		}
		else {
			std::cout << "Found " << device_count << " CUDA device(s) for G2 tests\n";
			cuda_device_available = true;
		}
	}

	static bool IsCudaAvailable() { return cuda_device_available; }
};

// Initialize the test environment
struct G2CudaTestRegistration {
	G2CudaTestRegistration() {
		::testing::AddGlobalTestEnvironment(new G2CUDAEnvironment());
	}
};

static G2CudaTestRegistration g2_cuda_test_registration;

// CUDA test cases comparing against CPU versions
TEST(G2OpsCUDA, PointAddition) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 10;
	POINTonE2* d_a, * d_b, * d_result;
	POINTonE2 h_a[num_points], h_b[num_points];
	POINTonE2 h_cuda_result[num_points], h_cpu_result[num_points];
	InsecureRng rng;

	// Create test data
	create_test_points_g2(h_a, num_points, rng);
	create_test_points_g2(h_b, num_points, rng);

	// CPU computation
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_add(&h_cpu_result[i], &h_a[i], &h_b[i]);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_a, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_b, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(POINTonE2)));

	CUDA_CHECK(cudaMemcpy(d_a, h_a, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_b, h_b, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));

	batched_error_code result = add_g2_vec(d_result, d_a, d_b, num_points);
	EXPECT_TRUE(check_batched_success(result));

	CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_points * sizeof(POINTonE2), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_points, point_equal_g2));

	// Cleanup
	CUDA_CHECK(cudaFree(d_a));
	CUDA_CHECK(cudaFree(d_b));
	CUDA_CHECK(cudaFree(d_result));
}

TEST(G2OpsCUDA, PointDoubling) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 10;
	POINTonE2* d_points, * d_result;
	POINTonE2 h_points[num_points];
	POINTonE2 h_cuda_result[num_points], h_cpu_result[num_points];

	// Create test data
	InsecureRng rng;
	create_test_points_g2(h_points, num_points, rng);

	// CPU computation
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_double(&h_cpu_result[i], &h_points[i]);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(POINTonE2)));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));

	batched_error_code result = double_g2_vec(d_result, d_points, num_points);
	EXPECT_TRUE(check_batched_success(result));

	CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_points * sizeof(POINTonE2), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_points, point_equal_g2));

	// Cleanup
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_result));
}

TEST(G2OpsCUDA, PointNegation) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 10;
	POINTonE2* d_points, * d_result;
	POINTonE2 h_points[num_points];
	POINTonE2 h_cuda_result[num_points], h_cpu_result[num_points];

	// Create test data
	InsecureRng rng;
	create_test_points_g2(h_points, num_points, rng);

	// CPU computation
	for (size_t i = 0; i < num_points; i++) {
		h_cpu_result[i] = h_points[i];
		blst_p2_cneg(&h_cpu_result[i], 1);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(POINTonE2)));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));

	batched_error_code result = cneg_g2_vec(d_result, d_points, 1, num_points);
	EXPECT_TRUE(check_batched_success(result));

	CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_points * sizeof(POINTonE2), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_points, point_equal_g2));

	// Cleanup
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_result));
}

TEST(G2OpsCUDA, ScalarMultiplication) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 10;
	const size_t scalar_bits = 128;
	const size_t scalar_bytes = (scalar_bits + 7) / 8;

	POINTonE2* d_points, * d_result;
	unsigned char* d_scalars;
	POINTonE2 h_points[num_points];
	unsigned char h_scalars[num_points * scalar_bytes];
	POINTonE2 h_cuda_result[num_points], h_cpu_result[num_points];

	// Create test data
	InsecureRng rng;
	create_test_points_g2(h_points, num_points, rng);
	for (size_t i = 0; i < num_points; i++) {
		memset(&h_scalars[i * scalar_bytes], 0, scalar_bytes);
		h_scalars[i * scalar_bytes] = (i + 1) % 255;
		if (h_scalars[i * scalar_bytes] == 0) h_scalars[i * scalar_bytes] = 1;
		h_scalars[i * scalar_bytes + 1] = (i * 3) % 255;
	}

	// CPU computation
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_mult(&h_cpu_result[i], &h_points[i], &h_scalars[i * scalar_bytes], scalar_bits);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_scalars, num_points * scalar_bytes));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_scalars, h_scalars, num_points * scalar_bytes, cudaMemcpyHostToDevice));

	batched_error_code result = mult_g2_vec(d_result, d_points, d_scalars, scalar_bits, num_points);
	EXPECT_TRUE(check_batched_success(result));

	CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_points * sizeof(POINTonE2), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_points, point_equal_g2));

	// Cleanup
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_result));
	CUDA_CHECK(cudaFree(d_scalars));
}

TEST(G2OpsCUDA, CurveMembership) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 10;
	POINTonE2* d_points;
	bool* d_results;
	POINTonE2 h_points[num_points];
	bool h_cuda_results[num_points], h_cpu_results[num_points];

	// Create test data (mix of valid and invalid points)
	InsecureRng rng;
	create_test_points_g2(h_points, num_points - 2, rng);

	// Add some invalid points
	h_points[num_points - 2] = h_points[0];
	h_points[num_points - 2].X[0][0] ^= 1;  // Flip a bit to make invalid

	h_points[num_points - 1] = h_points[1];
	h_points[num_points - 1].Y[0][0] ^= 1;  // Flip a bit to make invalid

	// CPU computation
	for (size_t i = 0; i < num_points; i++) {
		h_cpu_results[i] = blst_p2_on_curve(&h_points[i]);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_results, num_points * sizeof(bool)));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));

	batched_error_code result = on_curve_g2_vec(d_results, d_points, num_points);
	EXPECT_TRUE(check_batched_success(result));

	CUDA_CHECK(cudaMemcpy(h_cuda_results, d_results, num_points * sizeof(bool), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_results, h_cpu_results, num_points, bool_equal));

	// Cleanup
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_results));
}

TEST(G2OpsCUDA, Serialization) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 10;
	POINTonE2* d_points;
	POINTonE2_affine* d_affine_result;
	unsigned char* d_serialized;
	POINTonE2 h_points[num_points];
	POINTonE2_affine h_cuda_affine[num_points], h_cpu_affine[num_points];
	unsigned char h_serialized[num_points * 192];
	BLST_ERROR h_errors[num_points], h_errors_exp[num_points];

	// Create test data
	InsecureRng rng;
	create_test_points_g2(h_points, num_points, rng);

	// CPU computation
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_serialize(&h_serialized[i * 192], &h_points[i]);
		h_errors_exp[i] = blst_p2_deserialize(&h_cpu_affine[i], &h_serialized[i * 192]);
		EXPECT_EQ(h_errors_exp[i], BLST_SUCCESS);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_serialized, num_points * 192));
	CUDA_CHECK(cudaMalloc(&d_affine_result, num_points * sizeof(POINTonE2_affine)));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));

	batched_error_code serialize_result = serialize_g2_vec(d_serialized, d_points, num_points);
	EXPECT_TRUE(check_batched_success(serialize_result));

	batched_error_code deserialize_result = deserialize_g2_vec(d_affine_result, d_serialized, h_errors, num_points);
	EXPECT_TRUE(check_batched_success(deserialize_result));

	CUDA_CHECK(cudaMemcpy(h_cuda_affine, d_affine_result, num_points * sizeof(POINTonE2_affine), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_affine, h_cpu_affine, num_points, affine_equal_g2));
	EXPECT_TRUE(compare_results(h_errors, h_errors_exp, num_points, error_equal));

	// Cleanup
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_serialized));
	CUDA_CHECK(cudaFree(d_affine_result));
}

TEST(G2OpsCUDA, Compression) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 10;
	POINTonE2* d_points;
	POINTonE2_affine* d_affine_result;
	unsigned char* d_compressed;
	POINTonE2 h_points[num_points];
	POINTonE2_affine h_points_affine[num_points];
	POINTonE2_affine h_cuda_result[num_points], h_cpu_result[num_points];
	unsigned char h_compressed[num_points * 96];
	BLST_ERROR h_errors[num_points], h_errors_exp[num_points];

	// Create test data and convert to affine for CPU operations
	InsecureRng rng;
	create_test_points_g2(h_points, num_points, rng);
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_to_affine(&h_points_affine[i], &h_points[i]);
	}

	// CPU computation
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_affine_compress(&h_compressed[i * 96], &h_points_affine[i]);
		h_errors_exp[i] = blst_p2_uncompress(&h_cpu_result[i], &h_compressed[i * 96]);
		EXPECT_EQ(h_errors_exp[i], BLST_SUCCESS);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_compressed, num_points * 96));
	CUDA_CHECK(cudaMalloc(&d_affine_result, num_points * sizeof(POINTonE2_affine)));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));

	batched_error_code compress_result = compress_g2_vec(d_compressed, d_points, num_points);
	EXPECT_TRUE(check_batched_success(compress_result));

	batched_error_code uncompress_result = uncompress_g2_vec(d_affine_result, d_compressed, h_errors, num_points);
	EXPECT_TRUE(check_batched_success(uncompress_result));

	CUDA_CHECK(cudaMemcpy(h_cuda_result, d_affine_result, num_points * sizeof(POINTonE2_affine), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_points, affine_equal_g2));
	EXPECT_TRUE(compare_results(h_errors, h_errors_exp, num_points, error_equal));

	// Cleanup
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_compressed));
	CUDA_CHECK(cudaFree(d_affine_result));
}

TEST(G2OpsCUDA, CoordinateConversion) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 10;
	POINTonE2* d_points, * d_jacobian_result;
	POINTonE2_affine* d_affine;
	POINTonE2 h_points[num_points];
	POINTonE2_affine h_cuda_affine[num_points], h_cpu_affine[num_points];
	POINTonE2 h_cuda_jacobian[num_points], h_cpu_jacobian[num_points];

	// Create test data
	InsecureRng rng;
	create_test_points_g2(h_points, num_points, rng);

	// CPU computation - Jacobian to affine
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_to_affine(&h_cpu_affine[i], &h_points[i]);
	}

	// CPU computation - affine back to Jacobian
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_from_affine(&h_cpu_jacobian[i], &h_cpu_affine[i]);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_affine, num_points * sizeof(POINTonE2_affine)));
	CUDA_CHECK(cudaMalloc(&d_jacobian_result, num_points * sizeof(POINTonE2)));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));

	// Test Jacobian to affine conversion
	batched_error_code to_affine_result = to_affine_g2_vec(d_affine, d_points, num_points);
	EXPECT_TRUE(check_batched_success(to_affine_result));

	CUDA_CHECK(cudaMemcpy(h_cuda_affine, d_affine, num_points * sizeof(POINTonE2_affine), cudaMemcpyDeviceToHost));

	// Test affine to Jacobian conversion
	batched_error_code from_affine_result = from_affine_g2_vec(d_jacobian_result, d_affine, num_points);
	EXPECT_TRUE(check_batched_success(from_affine_result));

	CUDA_CHECK(cudaMemcpy(h_cuda_jacobian, d_jacobian_result, num_points * sizeof(POINTonE2), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_affine, h_cpu_affine, num_points, affine_equal_g2));
	EXPECT_TRUE(compare_results(h_cuda_jacobian, h_cpu_jacobian, num_points, point_equal_g2));

	// Cleanup
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_affine));
	CUDA_CHECK(cudaFree(d_jacobian_result));
}

TEST(G2OpsCUDA, HashToCurve) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_msgs = 5;
	POINTonE2* d_result;
	unsigned char** d_msg_ptrs;
	unsigned char* d_dst;
	size_t* d_msg_lens;
	POINTonE2 h_cuda_result[num_msgs], h_cpu_result[num_msgs];

	// Create test data
	std::vector<std::string> messages = {
		"message1", "message2", "message3", "message4", "message5"
	};
	std::string dst = "QUUX-V01-CS02-with-BLS12381G2";

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
		blst_hash_to_g2(&h_cpu_result[i],
			reinterpret_cast<const unsigned char*>(messages[i].c_str()),
			messages[i].length(),
			reinterpret_cast<const unsigned char*>(dst.c_str()),
			dst.length(),
			NULL, 0);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_result, num_msgs * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_msg_ptrs, num_msgs * sizeof(unsigned char*)));
	CUDA_CHECK(cudaMalloc(&d_msg_lens, num_msgs * sizeof(size_t)));
	CUDA_CHECK(cudaMalloc(&d_dst, dst.length()));

	// Copy message pointer array and lengths to device
	CUDA_CHECK(cudaMemcpy(d_msg_ptrs, h_msg_ptrs, num_msgs * sizeof(unsigned char*), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_msg_lens, h_msg_lens, num_msgs * sizeof(size_t), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_dst, dst.c_str(), dst.length(), cudaMemcpyHostToDevice));

	batched_error_code result = hash_to_g2_vec(d_result, d_msg_ptrs, d_msg_lens, d_dst, dst.length(), num_msgs);
	EXPECT_TRUE(check_batched_success(result));

	CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_msgs * sizeof(POINTonE2), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_msgs, point_equal_g2));

	// Cleanup
	for (size_t i = 0; i < num_msgs; i++) {
		CUDA_CHECK(cudaFree(d_msg_data[i]));
	}
	CUDA_CHECK(cudaFree(d_result));
	CUDA_CHECK(cudaFree(d_msg_ptrs));
	CUDA_CHECK(cudaFree(d_msg_lens));
	CUDA_CHECK(cudaFree(d_dst));
}

TEST(G2OpsCUDA, MixedAddition) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 10;
	POINTonE2* d_jacobian, * d_result;
	POINTonE2_affine* d_affine;
	POINTonE2 h_jacobian[num_points];
	POINTonE2_affine h_affine[num_points];
	POINTonE2 h_cuda_result[num_points], h_cpu_result[num_points];

	// Create test data
	InsecureRng rng;
	create_test_points_g2(h_jacobian, num_points, rng);

	// Create affine points (different from Jacobian points)
	POINTonE2 temp_points[num_points];
	create_test_points_g2(temp_points, num_points, rng);
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_to_affine(&h_affine[i], &temp_points[i]);
	}

	// CPU computation
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_add_affine(&h_cpu_result[i], &h_jacobian[i], &h_affine[i]);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_jacobian, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_affine, num_points * sizeof(POINTonE2_affine)));
	CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(POINTonE2)));

	CUDA_CHECK(cudaMemcpy(d_jacobian, h_jacobian, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_affine, h_affine, num_points * sizeof(POINTonE2_affine), cudaMemcpyHostToDevice));

	batched_error_code result = add_mixed_g2_vec(d_result, d_jacobian, d_affine, num_points);
	EXPECT_TRUE(check_batched_success(result));

	CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_points * sizeof(POINTonE2), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_points, point_equal_g2));

	// Cleanup
	CUDA_CHECK(cudaFree(d_jacobian));
	CUDA_CHECK(cudaFree(d_affine));
	CUDA_CHECK(cudaFree(d_result));
}

TEST(G2OpsCUDA, AdvancedOperations) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 10;
	POINTonE2* d_points, * d_result;
	bool* d_inf_results;
	POINTonE2 h_points[num_points];
	//POINTonE2 h_cuda_result[num_points], h_cpu_result[num_points];
	bool h_cuda_inf_results[num_points], h_cpu_inf_results[num_points];

	// Create test data
	InsecureRng rng;
	create_test_points_g2(h_points, num_points, rng);

	// CPU computation - test infinity check
	for (size_t i = 0; i < num_points; i++) {
		h_cpu_inf_results[i] = blst_p2_is_inf(&h_points[i]);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_inf_results, num_points * sizeof(bool)));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));

	// Test infinity check
	batched_error_code inf_result = is_inf_g2_vec(d_inf_results, d_points, num_points);
	EXPECT_TRUE(check_batched_success(inf_result));

	CUDA_CHECK(cudaMemcpy(h_cuda_inf_results, d_inf_results, num_points * sizeof(bool), cudaMemcpyDeviceToHost));

	// Compare infinity results
	EXPECT_TRUE(compare_results(h_cuda_inf_results, h_cpu_inf_results, num_points, bool_equal));

	// Cleanup
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_result));
	CUDA_CHECK(cudaFree(d_inf_results));
}

TEST(G2OpsCUDA, SubgroupMembership) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 10;
	POINTonE2_affine* d_points;
	bool* d_results;
	POINTonE2 h_points_jac[num_points];
	POINTonE2_affine h_points[num_points];
	bool h_cuda_results[num_points], h_cpu_results[num_points];

	// Create test data and convert to affine
	InsecureRng rng;
	create_test_points_g2(h_points_jac, num_points, rng);
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_to_affine(&h_points[i], &h_points_jac[i]);
	}

	// CPU computation
	for (size_t i = 0; i < num_points; i++) {
		h_cpu_results[i] = blst_p2_affine_in_g2(&h_points[i]);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2_affine)));
	CUDA_CHECK(cudaMalloc(&d_results, num_points * sizeof(bool)));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2_affine), cudaMemcpyHostToDevice));

	batched_error_code result = in_group_g2_vec(d_results, d_points, num_points);
	EXPECT_TRUE(check_batched_success(result));

	CUDA_CHECK(cudaMemcpy(h_cuda_results, d_results, num_points * sizeof(bool), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_results, h_cpu_results, num_points, bool_equal));

	// Cleanup
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_results));
}

// Performance and batch processing tests
TEST(G2OpsCUDA, PerformanceBaseline) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 1000;
	POINTonE2* d_points, * d_result, * h_points;
	unsigned char* d_scalars, * h_scalars;

	// Allocate host memory
	h_points = new POINTonE2[num_points];
	h_scalars = new unsigned char[num_points * 32];

	// Create test data
	InsecureRng rng;
	create_test_points_g2(h_points, num_points, rng);
	for (size_t i = 0; i < num_points; i++) {
		memset(&h_scalars[i * 32], 0, 32);
		h_scalars[i * 32] = (i % 255) + 1;
	}

	// Allocate device memory
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_scalars, num_points * 32));

	// Copy to device
	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_scalars, h_scalars, num_points * 32, cudaMemcpyHostToDevice));

	// Record timing
	cudaEvent_t start, stop;
	CUDA_CHECK(cudaEventCreate(&start));
	CUDA_CHECK(cudaEventCreate(&stop));

	CUDA_CHECK(cudaEventRecord(start));
	batched_error_code result = mult_g2_vec(d_result, d_points, d_scalars, 8, num_points);
	CUDA_CHECK(cudaEventRecord(stop));

	CUDA_CHECK(cudaGetLastError());
	CUDA_CHECK(cudaDeviceSynchronize());

	float milliseconds = 0;
	CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));

	printf("CUDA G2 Performance: %zu scalar multiplications in %.3f ms (%.1f ops/ms)\n",
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

TEST(G2OpsCUDA, BatchProcessing) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t batch_size = 256;  // Lower than G1 since G2 operations are more complex
	POINTonE2* d_points, * d_result, * h_points;
	POINTonE2* h_cuda_result, * h_cpu_result;

	// Allocate host memory
	h_points = new POINTonE2[batch_size];
	h_cuda_result = new POINTonE2[batch_size];
	h_cpu_result = new POINTonE2[batch_size];

	// Create test data
	InsecureRng rng;
	create_test_points_g2(h_points, batch_size, rng);

	// CPU computation
	for (size_t i = 0; i < batch_size; i++) {
		blst_p2_double(&h_cpu_result[i], &h_points[i]);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, batch_size * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_result, batch_size * sizeof(POINTonE2)));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, batch_size * sizeof(POINTonE2), cudaMemcpyHostToDevice));

	batched_error_code double_result = double_g2_vec(d_result, d_points, batch_size);
	EXPECT_TRUE(check_batched_success(double_result));

	CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, batch_size * sizeof(POINTonE2), cudaMemcpyDeviceToHost));

	// Compare results
	bool all_match = compare_results(h_cuda_result, h_cpu_result, batch_size, point_equal_g2);

	printf("G2 Batch processing: %s (CPU vs CUDA comparison)\n",
		all_match ? "All results match" : "Some results differ");

	EXPECT_TRUE(all_match);

	// Cleanup
	delete[] h_points;
	delete[] h_cuda_result;
	delete[] h_cpu_result;
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_result));
}

// Additional G2-specific tests
TEST(G2OpsCUDA, SerializationFormats) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 5;
	POINTonE2* d_points;
	bool* d_results;
	POINTonE2 h_points[num_points];
	bool h_cuda_results[num_points], h_cpu_results[num_points];

	// Create test data including infinity point
	InsecureRng rng;
	create_test_points_g2(h_points, num_points - 1, rng);

	// Add infinity point
	memset(&h_points[num_points - 1], 0, sizeof(POINTonE2));

	// CPU computation - test point validation
	for (size_t i = 0; i < num_points; i++) {
		h_cpu_results[i] = blst_p2_on_curve(&h_points[i]);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_results, num_points * sizeof(bool)));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));

	batched_error_code result = on_curve_g2_vec(d_results, d_points, num_points);
	EXPECT_TRUE(check_batched_success(result));

	CUDA_CHECK(cudaMemcpy(h_cuda_results, d_results, num_points * sizeof(bool), cudaMemcpyDeviceToHost));

	// Compare results
	EXPECT_TRUE(compare_results(h_cuda_results, h_cpu_results, num_points, bool_equal));

	// Cleanup
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_results));
}

TEST(G2OpsCUDA, LargeScalarMultiplication) {
	if (!G2CUDAEnvironment::IsCudaAvailable()) {
		GTEST_SKIP() << "CUDA not available";
	}

	const size_t num_points = 1;
	const size_t scalar_bits = 252;
	const size_t scalar_bytes = (scalar_bits + 7) / 8;
	POINTonE2* d_points, * d_result;
	unsigned char* d_scalars;
	POINTonE2 h_points[num_points];
	POINTonE2 h_cuda_result[num_points], h_cpu_result[num_points];

	// Create large scalars close to the order
	unsigned char h_scalars[num_points * scalar_bytes];
	for (size_t i = 0; i < num_points; i++) {
		// Initialize with large values
		for (size_t j = 0; j < scalar_bytes; j++) {
			h_scalars[i * scalar_bytes + j] = (j < scalar_bytes-1) ? 0xFF : 0x0F;
		}
		// Make each scalar slightly different
		h_scalars[i * scalar_bytes] = (i + 1) & 0xFF;
	}

	// Create test data
	InsecureRng rng;
	create_test_points_g2(h_points, num_points, rng);
	// CPU computation
	for (size_t i = 0; i < num_points; i++) {
		blst_p2_mult(&h_cpu_result[i], &h_points[i], &h_scalars[i * scalar_bytes], scalar_bits);
	}

	// CUDA computation
	CUDA_CHECK(cudaMalloc(&d_points, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_result, num_points * sizeof(POINTonE2)));
	CUDA_CHECK(cudaMalloc(&d_scalars, num_points * scalar_bytes));

	CUDA_CHECK(cudaMemcpy(d_points, h_points, num_points * sizeof(POINTonE2), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d_scalars, h_scalars, num_points * scalar_bytes, cudaMemcpyHostToDevice));

	batched_error_code result = mult_g2_vec(
		d_result, d_points, d_scalars, scalar_bits, num_points);
	EXPECT_TRUE(check_batched_success(result));


	CUDA_CHECK(cudaMemcpy(h_cuda_result, d_result, num_points * sizeof(POINTonE2), cudaMemcpyDeviceToHost));

	// Compare results
	//bool results_match = true;
	for (size_t i = 0; i < num_points; i++) {
		// Verify CUDA result is on curve
		if (!blst_p2_on_curve(&h_cuda_result[i])) {
			printf("Error G2 Large: CUDA result %zu is not on curve\n", i);
			//results_match = false;
		}

		if (!blst_p2_is_equal(&h_cuda_result[i], &h_cpu_result[i])) {
			printf("G2 Large Mismatch at index %zu\n", i);
			printf("CPU  result: X[0][0]=%llx, Y[0][0]=%llx, Z[0][0]=%llx\n",
				h_cpu_result[i].X[0][0], h_cpu_result[i].Y[0][0], h_cpu_result[i].Z[0][0]);
			printf("CUDA result: X[0][0]=%llx, Y[0][0]=%llx, Z[0][0]=%llx\n",
				h_cuda_result[i].X[0][0], h_cuda_result[i].Y[0][0], h_cuda_result[i].Z[0][0]);
			printf("Large scalar used (first 8 bytes): ");
			for (int j = 0; j < 8; j++) {
				printf("%02x ", h_scalars[i * 32 + j]);
			}
			printf("\nLarge scalar used (last 8 bytes): ");
			for (int j = 24; j < 32; j++) {
				printf("%02x ", h_scalars[i * 32 + j]);
			}
			printf("\n");

			// Additional G2-specific debugging - check both components of Fp2
			printf("CPU  X[1][0]=%llx, Y[1][0]=%llx, Z[1][0]=%llx\n",
				h_cpu_result[i].X[1][0], h_cpu_result[i].Y[1][0], h_cpu_result[i].Z[1][0]);
			printf("CUDA X[1][0]=%llx, Y[1][0]=%llx, Z[1][0]=%llx\n",
				h_cuda_result[i].X[1][0], h_cuda_result[i].Y[1][0], h_cuda_result[i].Z[1][0]);

			//results_match = false;
		}
	}
	EXPECT_TRUE(compare_results(h_cuda_result, h_cpu_result, num_points, point_equal_g2));

	// Cleanup
	CUDA_CHECK(cudaFree(d_points));
	CUDA_CHECK(cudaFree(d_result));
	CUDA_CHECK(cudaFree(d_scalars));
}