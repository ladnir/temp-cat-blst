
#include <gtest/gtest.h>

#ifdef BLST_CUDA
#include <cuda_runtime.h>
#endif

#include "blst/curve/GT.h"
#include "blst/curve/multi_scalar_par.h"
#include "InsecureRng.h"
#include "device.h"
using namespace blst;



template<Backend B>
void sort_digits_segmented_test()
{
	//if (!G1CUDAEnvironment::IsCudaAvailable()) {
	//    GTEST_SKIP() << "CUDA not available";
	//}

//#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)
//	const Backend B = Backend::CUDA;
//#else
//	const Backend B = Backend::CPU;
//#endif

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
	CHECK(deviceMalloc<B>(&d_keys_input, sizeof(u8) * N * D));
	CHECK(deviceMalloc<B>(&d_keys_output, sizeof(u8) * N * D));
	CHECK(deviceMalloc<B>(&d_vals_input, sizeof(blst::IdxSign) * N * D));
	CHECK(deviceMalloc<B>(&d_vals_output, sizeof(blst::IdxSign) * N * D));

	// Copy test data to device
	CHECK(deviceMemcpy<B>(d_keys_input, h_keys_input, sizeof(u8) * N * D, MemcpyKind::HostToDevice));
	CHECK(deviceMemcpy<B>(d_vals_input, h_vals_input, sizeof(blst::IdxSign) * N * D, MemcpyKind::HostToDevice));

	// Query workspace size
	size_t workspace_bytes = 0;
	DeviceError result = blst::sort_digits_segmented<B>(
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

	EXPECT_EQ(result, DeviceError::cudaSuccess);
	EXPECT_GT(workspace_bytes, 0) << "Workspace size should be positive";

	// Allocate workspace
	void* d_workspace = nullptr;
	CHECK(deviceMalloc<B>(&d_workspace, workspace_bytes));

	// Call the function under test
	result = blst::sort_digits_segmented<B>(
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

	EXPECT_EQ(result, DeviceError::cudaSuccess);

	// Copy results back to host
	CHECK(deviceMemcpy<B>(h_keys_output, d_keys_output, sizeof(u8) * N * D, MemcpyKind::DeviceToHost));
	CHECK(deviceMemcpy<B>(h_vals_output, d_vals_output, sizeof(blst::IdxSign) * N * D, MemcpyKind::DeviceToHost));

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

	CHECK(deviceFree<B>(d_keys_input));
	CHECK(deviceFree<B>(d_keys_output));
	CHECK(deviceFree<B>(d_vals_input));
	CHECK(deviceFree<B>(d_vals_output));
	CHECK(deviceFree<B>(d_workspace));
}



TEST(MSM, CpuSortDigitsParallel) {
	sort_digits_segmented_test<Backend::CPU>();
}
#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)
TEST(MSM, CudaSortDigitsParallel) {
	sort_digits_segmented_test<Backend::CUDA>();
}
#endif

// Additional edge case test
//TEST(MSM, SortDigitsParallelEdgeCases) {

template<Backend B>
void sort_digits_segmented_edge_cases_test()
{
	//if (!G1CUDAEnvironment::IsCudaAvailable()) {
	//    GTEST_SKIP() << "CUDA not available";
	//}

	// Test with minimal data
	const int N_small = 10;
	const int D_small = 2;
	const int wbits_small = 3; // 8 possible values

	u8* d_keys_input;
	u8* d_keys_output;
	blst::IdxSign* d_vals_input;
	blst::IdxSign* d_vals_output;
	CHECK(deviceMalloc<B>(&d_keys_input, sizeof(u8) * N_small * D_small));
	CHECK(deviceMalloc<B>(&d_keys_output, sizeof(u8) * N_small * D_small));
	CHECK(deviceMalloc<B>(&d_vals_input, sizeof(blst::IdxSign) * N_small * D_small));
	CHECK(deviceMalloc<B>(&d_vals_output, sizeof(blst::IdxSign) * N_small * D_small));

	// Test data: already sorted
	u8 h_sorted_keys[] = { 0, 1, 2, 3, 4, 5, 6, 7, 0, 1,  // digit 0
						  0, 0, 1, 1, 2, 2, 3, 3, 4, 4 }; // digit 1
	blst::IdxSign h_sorted_vals[N_small * D_small];
	for (int i = 0; i < N_small * D_small; i++) {
		h_sorted_vals[i] = blst::IdxSign{ static_cast<uint32_t>(i << 1) }; // No sign bit
	}

	CHECK(deviceMemcpy<B>(d_keys_input, h_sorted_keys,
		sizeof(u8) * N_small * D_small, MemcpyKind::HostToDevice));
	CHECK(deviceMemcpy<B>(d_vals_input, h_sorted_vals,
		sizeof(blst::IdxSign) * N_small * D_small, MemcpyKind::HostToDevice));

	// Query workspace size
	size_t workspace_bytes = 0;
	DeviceError result = blst::sort_digits_segmented<B>(
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

	EXPECT_EQ(result, DeviceError::cudaSuccess);

	// Allocate workspace
	void* d_workspace = nullptr;
	CHECK(deviceMalloc<B>(&d_workspace, workspace_bytes));

	// Should handle already sorted data gracefully
	result = blst::sort_digits_segmented<B>(
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

	EXPECT_EQ(result, DeviceError::cudaSuccess);

	// Test with reverse sorted data
	u8 h_reverse_keys[] = { 7, 6, 5, 4, 3, 2, 1, 0, 1, 0,  // digit 0
						   4, 4, 3, 3, 2, 2, 1, 1, 0, 0 }; // digit 1

	CHECK(deviceMemcpy<B>(d_keys_input, h_reverse_keys,
		sizeof(u8) * N_small * D_small, MemcpyKind::HostToDevice));

	result = blst::sort_digits_segmented<B>(
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

	EXPECT_EQ(result, DeviceError::cudaSuccess);

	// Verify the reverse sorted data is now properly sorted
	u8 h_result_keys[N_small * D_small];
	CHECK(deviceMemcpy<B>(h_result_keys, d_keys_output,
		sizeof(u8) * N_small * D_small, MemcpyKind::DeviceToHost));

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

	CHECK(deviceFree<B>(d_keys_input));
	CHECK(deviceFree<B>(d_keys_output));
	CHECK(deviceFree<B>(d_vals_input));
	CHECK(deviceFree<B>(d_vals_output));
	CHECK(deviceFree<B>(d_workspace));
}

TEST(MSM, CpuSortDigitsParallelEdgeCases) {
	sort_digits_segmented_edge_cases_test<Backend::CPU>();
}
#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)
	
TEST(MSM, CudaSortDigitsParallelEdgeCases) {
	sort_digits_segmented_edge_cases_test<Backend::CUDA>();
}
#endif

//TEST(MSM, ReduceGroups) {
template<Backend B>
void reduce_groups_test()
{

	//if (!G1CUDAEnvironment::IsCudaAvailable()) {
	//    GTEST_SKIP() << "CUDA not available";
	//}
//
//#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)
//	const Backend B = Backend::CUDA;
//#else
//	const Backend B = Backend::CPU;
//#endif
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
	create_test_group(h_points.data(), N, rng);

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

	CHECK(deviceMalloc<B>(&d_keys, M * N * sizeof(u8)));
	CHECK(deviceMalloc<B>(&d_points, M * N * sizeof(G1AffinePoint)));
	CHECK(deviceMalloc<B>(&d_groups, M * W * sizeof(Pointxyzz<Fp>)));
	CHECK(deviceMalloc<B>(&d_vals, M * N * sizeof(blst::IdxSign)));

	// Copy data to device
	CHECK(deviceMemcpy<B>(d_keys, h_keys.data(), M * N * sizeof(u8), MemcpyKind::HostToDevice));
	CHECK(deviceMemcpy<B>(d_points, h_affine_points.data(), M * N * sizeof(G1AffinePoint), MemcpyKind::HostToDevice));
	CHECK(deviceMemcpy<B>(d_vals, h_vals.data(), M * N * sizeof(blst::IdxSign), MemcpyKind::HostToDevice));

	// Create transform iterator for gathering points with sign
	auto device_gather_iter = transform_iterator(
		d_vals,
		blst::PackedToPointXYZZ<AffinePoint<Fp>, Defaulted<Pointxyzz<Fp>>>{ d_points, N }
	);


	// Query workspace size - using batched version (M MSMs, N elements each, W groups)
	size_t workspace_bytes = 0;
	DeviceError result = blst::reduce_groups<B, Pointxyzz<Fp>>(
		nullptr,          // d_workspace (nullptr to query size)
		workspace_bytes,  // in/out workspace size
		d2Span<u8>(d_keys, M, N),
		device_gather_iter,
		d2Span<Pointxyzz<Fp>>(d_groups, M, W),
		M,                // Number of MSMs
		N,                // Number of elements per MSM
		W,                // Number of groups per MSM
		0);               // default stream

	EXPECT_EQ(result, DeviceError::cudaSuccess);
	EXPECT_GT(workspace_bytes, 0) << "Workspace size should be positive";

	// Allocate workspace
	CHECK(deviceMalloc<B>(&d_workspace, workspace_bytes));

	// Perform the actual batched reduction
	result = blst::reduce_groups<B, Pointxyzz<Fp>>(
		d_workspace,
		workspace_bytes,
		d2Span<u8>(d_keys, M, N),
		device_gather_iter,
		d2Span<Pointxyzz<Fp>>(d_groups, M, W),
		M,                // Number of MSMs
		N,                // Number of elements per MSM
		W,                // Number of groups per MSM
		0);

	// Block on stream to catch errors
	CHECK(streamSynchronize<B>(0));

	if (result != DeviceError::cudaSuccess) {
		const char* err_str = getErrorString<B>(result);
		printf("batched reduce_groups failed: %s\n", err_str);
	}

	EXPECT_EQ(result, DeviceError::cudaSuccess);

	// Copy results back to host
	CHECK(deviceMemcpy<B>(h_cuda_groups_xyzz.data(), d_groups, M * W * sizeof(Pointxyzz<Fp>), MemcpyKind::DeviceToHost));

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
	CHECK(deviceFree<B>(d_keys));
	CHECK(deviceFree<B>(d_vals));
	CHECK(deviceFree<B>(d_points));
	CHECK(deviceFree<B>(d_groups));
	CHECK(deviceFree<B>(d_workspace));
}

TEST(MSM, CpuReduceGroups) {
	reduce_groups_test<Backend::CPU>();
}
#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)
TEST(MSM, CudaReduceGroups) {
	reduce_groups_test<Backend::CUDA>();
}
#endif


//TEST(MSM, ReduceGroupsEdgeCases) 
template<Backend B>
void reduce_groups_edge_cases_test()
{
//#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)
//	const Backend B = Backend::CUDA;
//#else
//	const Backend B = Backend::CPU;
//#endif
	//if (!G1CUDAEnvironment::IsCudaAvailable()) {
	//    GTEST_SKIP() << "CUDA not available";
	//}

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
	create_test_group(h_points.data(), N, rng);
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

	CHECK(deviceMalloc<B>(&d_keys, M * N * sizeof(u8)));
	CHECK(deviceMalloc<B>(&d_vals, M * N * sizeof(blst::IdxSign)));
	CHECK(deviceMalloc<B>(&d_points, N * sizeof(G1AffinePoint)));
	CHECK(deviceMalloc<B>(&d_groups, M * W * sizeof(Pointxyzz<Fp>)));

	CHECK(deviceMemcpy<B>(d_keys, h_keys.data(), M * N * sizeof(u8), MemcpyKind::HostToDevice));
	CHECK(deviceMemcpy<B>(d_vals, h_vals.data(), M * N * sizeof(blst::IdxSign), MemcpyKind::HostToDevice));
	CHECK(deviceMemcpy<B>(d_points, h_affine_points.data(), N * sizeof(G1AffinePoint), MemcpyKind::HostToDevice));

	// Create transform iterator
	auto gather_iter = transform_iterator(
		d_vals,
		blst::PackedToPointXYZZ<AffinePoint<Fp>, Defaulted<Pointxyzz<Fp>>>{ d_points, N }
	);


	// Test workspace size query for batched version
	size_t workspace_bytes = 0;
	DeviceError result = blst::reduce_groups<B, Pointxyzz<Fp>>(
		nullptr, workspace_bytes,
		d2Span<u8>(d_keys, M, N), gather_iter,
		d2Span<Pointxyzz<Fp>>(d_groups, M, W),
		M, N, W, 0);

	EXPECT_EQ(result, DeviceError::cudaSuccess);

	// Test with insufficient workspace
	void* d_workspace;
	size_t small_workspace = workspace_bytes / 2;
	CHECK(deviceMalloc<B>(&d_workspace, small_workspace));

	result = blst::reduce_groups<B, Pointxyzz<Fp>>(
		d_workspace, small_workspace,
		d2Span<u8>(d_keys, M, N), gather_iter,
		d2Span<Pointxyzz<Fp>>(d_groups, M, W),
		M, N, W, 0);

	EXPECT_EQ(result, DeviceError::cudaErrorInvalidValue) << "Should fail with insufficient workspace";

	// Test with correct workspace size
	CHECK(deviceFree<B>(d_workspace));
	CHECK(deviceMalloc<B>(&d_workspace, workspace_bytes));

	result = blst::reduce_groups<B, Pointxyzz<Fp>>(
		d_workspace, workspace_bytes,
		d2Span<u8>(d_keys, M, N), gather_iter,
		d2Span<Pointxyzz<Fp>>(d_groups, M, W),
		M, N, W, 0);

	EXPECT_EQ(result, DeviceError::cudaSuccess);

	//printf("Batched reduce_groups edge cases: %d MSMs, %d points each, %d groups per MSM - PASSED\n", M, N, W);

	// Cleanup
	CHECK(deviceFree<B>(d_keys));
	CHECK(deviceFree<B>(d_vals));
	CHECK(deviceFree<B>(d_points));
	CHECK(deviceFree<B>(d_groups));
	CHECK(deviceFree<B>(d_workspace));
}
TEST(MSM, CpuReduceGroupsEdgeCases) {
	reduce_groups_edge_cases_test<Backend::CPU>();
}
#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)
	
TEST(MSM, CudaReduceGroupsEdgeCases) {
	reduce_groups_edge_cases_test<Backend::CUDA>();
}
#endif


template<Backend B, typename Group, typename AffineGroup, typename GroupZZ>
void msm_booth_batched_pipeline(
	int N,
	int M,
	int nbits,
	int wbits,
	bool debug = false,
	bool print = false
)
{
	auto nbytes = (nbits + 7) / 8;

	// Host data - now M*N points total, N points per MSM
	std::vector<AffineGroup> h_points(M * N);  // M*N total points
	std::vector<unsigned char> h_scalars(M * N * nbytes);
	std::vector<Group> h_cuda_outputs(M), h_cpu_outputs(M);

	// Create test data
	InsecureRng rng;

	// Generate M*N unique points (each MSM gets its own set of N points)
	std::vector<Group> temp_points(M * N);
	create_test_group(temp_points.data(), M * N, rng);
	for (int i = 0; i < M * N; i++) {
		msm_to_affine(h_points[i], temp_points[i]);
	}

	// Generate random scalars for each MSM
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			for (int b = 0; b < nbytes; b++) {
				h_scalars[m * N * nbytes + n * nbytes + b] = rng.uniform_u32(0, 255);
			}
		}
	}

	// CPU computation for verification - each MSM uses its own set of N points
	for (int m = 0; m < M; m++) {
		std::vector<const AffineGroup*> points_for_msm(2);
		std::vector<const byte*> scalars_for_msm(2);

		points_for_msm[0] = &h_points[m * N];
		scalars_for_msm[0] = &h_scalars[m * N * nbytes];

		// Each MSM m uses points [m*N ... (m+1)*N-1]
		//for (int n = 0; n < N; n++) {
		//    points_for_msm[n] = &h_points[m * N + n];  // MSM m uses points starting at m*N
		//    scalars_for_msm[n] = &h_scalars[m * N * nbytes + n * nbytes];
		//}

		// Use CPU Pippenger for reference
		size_t scratch_size = blst::msm_mult_pippenger_scratch_sizeof_<GroupZZ>(N);
		GroupZZ* scratch = (GroupZZ*)malloc(scratch_size);

		blst::msm_mult_pippenger<Group, AffineGroup, GroupZZ>(
			h_cpu_outputs[m],
			points_for_msm.data(),
			N,
			scalars_for_msm.data(),
			nbits,
			scratch,
			wbits,
			m == 1 && print
		);

		//     auto w = wbits;
		//     auto D = 1 + nbits / w;
		//     // Test data with known bit patterns
		//     //std::vector<uint8_t> scalars(N * nbytes);
		//     std::vector<u8> digits(D * N);
		//     std::vector<blst::IdxSign> index_signs(D * N);

		//     // Extract digits using the function under test
		//     for (int d = D-1; d >= 0; --d) {
				 //std::cout << "Pippenger window " << d << " digits:\n";
		//         for (int i = 0; i < N; i++) {
		//             blst::extract_digits_flat(
		//                 d, i,
		//                 scalars_for_msm[0],
		//                 N, nbytes, nbits, wbits, D,
		//                 digits.data(),
		//                 index_signs.data()
		//             );
		//         }
		//     }

		free(scratch);
	}

	// Device memory allocation - now for M*N points
	AffineGroup* d_points_ptr;
	uint8_t* d_scalars_ptr;
	Group* d_outputs_ptr;

	CHECK(deviceMalloc<B>(&d_points_ptr, M * N * sizeof(AffineGroup)));  // M*N points
	CHECK(deviceMalloc<B>(&d_scalars_ptr, M * N * nbytes));
	CHECK(deviceMalloc<B>(&d_outputs_ptr, M * sizeof(Group)));

	d2Span<AffineGroup> d_points(d_points_ptr, M, N);
	d3Span<uint8_t> d_scalars(d_scalars_ptr, M, N, nbytes);
	span<Group> d_outputs(d_outputs_ptr, M);

	// Copy data to device - all M*N points
	CHECK(deviceMemcpy<B>(d_points.data_handle(), h_points.data(), M * N * sizeof(AffineGroup), MemcpyKind::HostToDevice));
	CHECK(deviceMemcpy<B>(d_scalars.data_handle(), h_scalars.data(), M * N * nbytes, MemcpyKind::HostToDevice));

	// Query workspace size
	size_t workspace_bytes = 0;
	auto result = blst::msm_booth_batched_pipeline<B, Group, AffineGroup, GroupZZ>(
		nullptr,
		workspace_bytes,
		d_points,
		d_scalars,
		nbytes,
		nbits,
		N,
		M,
		wbits,
		d_outputs,
		0,
		debug);

	EXPECT_EQ(result, DeviceError::cudaSuccess);
	EXPECT_GT(workspace_bytes, 0) << "Workspace size should be positive";

	// Allocate workspace
	void* d_workspace;
	CHECK(deviceMalloc<B>(&d_workspace, workspace_bytes));

	// Run the batched MSM pipeline
	result = blst::msm_booth_batched_pipeline<B, Group, AffineGroup, GroupZZ>(
		d_workspace,
		workspace_bytes,
		d_points,
		d_scalars,
		nbytes,
		nbits,
		N,
		M,
		wbits,
		d_outputs,
		0,
		debug, print);

	EXPECT_EQ(result, DeviceError::cudaSuccess);

	// sync the stream
	CHECK(streamSynchronize<B>(0));

	// Copy results back
	CHECK(deviceMemcpy<B>(h_cuda_outputs.data(), d_outputs.data(), M * sizeof(Group), MemcpyKind::DeviceToHost));

	// Compare CUDA results with CPU results
	bool all_match = true;
	for (int m = 0; m < M; m++) {
		if (!Group::equal(h_cuda_outputs[m], h_cpu_outputs[m])) {
			printf("MSM %d mismatch\n", m);
			all_match = false;
		}
	}

	EXPECT_TRUE(all_match) << "CUDA and CPU MSM results should match";

	// Cleanup
	CHECK(deviceFree<B>(d_points.data_handle()));
	CHECK(deviceFree<B>(d_scalars.data_handle()));
	CHECK(deviceFree<B>(d_outputs.data()));
	CHECK(deviceFree<B>(d_workspace));
}

TEST(MSM, CpuBoothBatchedPipeline) {


	bool debug = false;
	bool print = false;

	// Test parameters
	const int N = 33;        // points per MSM
	const int M = 3;         // number of MSMs
	const int nbits = 7;   // scalar bit length
	const int wbits = 3;     // window size

	   msm_booth_batched_pipeline<Backend::CPU, Point<Fp>, AffinePoint<Fp>, Pointxyzz<Fp>>(
	       N,
	       M,
	       nbits,
	       wbits,
	       debug,
	       print
	   );


	msm_booth_batched_pipeline<Backend::CPU, Point<Fp2>, AffinePoint<Fp2>, Pointxyzz<Fp2>>(
		N,
		M,
		nbits,
		wbits,
		debug,
		print
	);


	msm_booth_batched_pipeline<Backend::CPU, GT, GT, GT>(
	    N,
	    M,
	    nbits,
	    wbits,
	    debug,
	    print
	);


}

#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)

TEST(MSM, CudaBoothBatchedPipeline) {


	bool debug = false;
	bool print = false;

	// Test parameters
	const int N = 33;        // points per MSM
	const int M = 3;         // number of MSMs
	const int nbits = 7;   // scalar bit length
	const int wbits = 3;     // window size

	msm_booth_batched_pipeline<Backend::CUDA, Point<Fp>, AffinePoint<Fp>, Pointxyzz<Fp>>(
		N,
		M,
		nbits,
		wbits,
		debug,
		print
	);


	msm_booth_batched_pipeline<Backend::CUDA, Point<Fp2>, AffinePoint<Fp2>, Pointxyzz<Fp2>>(
		N,
		M,
		nbits,
		wbits,
		debug,
		print
	);


	msm_booth_batched_pipeline<Backend::CUDA, GT, GT, GT>(
		N,
		M,
		nbits,
		wbits,
		debug,
		print
	);


}

#endif

//TEST(MSM, MSMBoothBatchedPipelineLargeScale) {
//    if (!G1CUDAEnvironment::IsCudaAvailable()) {
//        GTEST_SKIP() << "CUDA not available";
//    }
//
//    // Test with larger scale for performance validation
//    const int N = 256;       // points per MSM
//    const int M = 8;         // number of MSMs
//    const int nbits = 256;   // full scalar size
//    const int nbytes = (nbits + 7) / 8;
//    const int wbits = 6;     // larger window
//
//    // Use dynamic allocation for larger arrays - M*N total points
//    std::vector<G1AffinePoint> h_points(M * N);  // M*N total points
//    std::vector<unsigned char> h_scalars(M * N * nbytes);
//    std::vector<G1Point> h_cuda_outputs(M);
//
//    // Create test data
//    InsecureRng rng;
//
//    // Generate M*N unique points
//    std::vector<G1Point> temp_points(M * N);
//    create_test_group(temp_points.data(), M * N, rng);
//    for (int i = 0; i < M * N; i++) {
//        G1Point::to_affine(h_points[i], temp_points[i]);
//    }
//
//    // Generate scalars
//    for (size_t i = 0; i < h_scalars.size(); i++) {
//        h_scalars[i] = rng.uniform_u32(0, 255);
//    }
//
//    // Device allocation - M*N points
//    G1AffinePoint* d_points;
//    uint8_t* d_scalars;
//    G1Point* d_outputs;
//
//    CHECK(deviceMalloc(&d_points, M * N * sizeof(G1AffinePoint)));  // M*N points
//    CHECK(deviceMalloc(&d_scalars, M * N * nbytes));
//    CHECK(deviceMalloc(&d_outputs, M * sizeof(G1Point)));
//
//    CHECK(deviceMemcpy(d_points, h_points.data(), M * N * sizeof(G1AffinePoint), MemcpyKind::HostToDevice));
//    CHECK(deviceMemcpy(d_scalars, h_scalars.data(), M * N * nbytes, MemcpyKind::HostToDevice));
//
//    // Timing test
//    cudaEvent_t start, stop;
//    CHECK(cudaEventCreate(&start));
//    CHECK(cudaEventCreate(&stop));
//
//    size_t workspace_bytes = 0;
//    DeviceError result = blst::msm_booth_batched_pipeline<blst::Fp>(
//        nullptr, workspace_bytes, d_points, d_scalars, nbytes, nbits, N, M, wbits, d_outputs, 0);
//
//    EXPECT_EQ(result,DeviceError::cudaSuccess);
//
//    void* d_workspace;
//    CHECK(deviceMalloc(&d_workspace, workspace_bytes));
//
//    CHECK(cudaEventRecord(start));
//    result = blst::msm_booth_batched_pipeline<blst::Fp>(
//        d_workspace, workspace_bytes, d_points, d_scalars, nbytes, nbits, N, M, wbits, d_outputs, 0);
//    CHECK(cudaEventRecord(stop));
//
//    CHECK(cudaEventSynchronize(stop));
//
//    float milliseconds = 0;
//    CHECK(cudaEventElapsedTime(&milliseconds, start, stop));
//
//    printf("Batched MSM Performance: %d MSMs of %d points each (%d total points) in %.3f ms (%.1f MSMs/sec)\n",
//        M, N, M * N, milliseconds, M * 1000.0f / milliseconds);
//
//    EXPECT_EQ(result,DeviceError::cudaSuccess);
//
//    // Copy results back
//    CHECK(deviceMemcpy(h_cuda_outputs.data(), d_outputs, M * sizeof(G1Point), MemcpyKind::DeviceToHost));
//
//    // Basic sanity check - ensure results are not all infinity
//    bool has_non_infinity = false;
//    for (int m = 0; m < M; m++) {
//        if (!blst::detail::is_infinity(h_cuda_outputs[m])) {
//            has_non_infinity = true;
//            break;
//        }
//    }
//    EXPECT_TRUE(has_non_infinity) << "At least some MSM results should be non-infinity";
//
//    // Cleanup
//    CHECK(deviceFree(d_points));
//    CHECK(deviceFree(d_scalars));
//    CHECK(deviceFree(d_outputs));
//    CHECK(deviceFree(d_workspace));
//    CHECK(deviceEventDestroy(start));
//    CHECK(deviceEventDestroy(stop));
//}

//TEST(MSM, MSMBoothBatchedPipelineWindowSizes) {
//    if (!G1CUDAEnvironment::IsCudaAvailable()) {
//        GTEST_SKIP() << "CUDA not available";
//    }
//
//    // Test different window sizes
//    const int N = 32;
//    const int M = 2;
//    const int nbits = 128;
//    const int nbytes = (nbits + 7) / 8;
//    const std::vector<int> window_sizes = { 3, 4, 5, 6, 7, 8 };
//
//    // Generate consistent test data - M*N total points
//    InsecureRng rng(12345);  // Fixed seed for reproducibility
//
//    std::vector<G1AffinePoint> h_points(M * N);  // M*N total points
//    std::vector<unsigned char> h_scalars(M * N * nbytes);
//
//    std::vector<G1Point> temp_points(M * N);
//    create_test_group(temp_points.data(), M * N, rng);
//    for (int i = 0; i < M * N; i++) {
//        G1Point::to_affine(h_points[i], temp_points[i]);
//    }
//
//    for (size_t i = 0; i < h_scalars.size(); i++) {
//        h_scalars[i] = rng.uniform_u32(0, 255);
//    }
//
//    // Device allocation - M*N points
//    G1AffinePoint* d_points;
//    uint8_t* d_scalars;
//    G1Point* d_outputs;
//
//    CHECK(deviceMalloc(&d_points, M * N * sizeof(G1AffinePoint)));  // M*N points
//    CHECK(deviceMalloc(&d_scalars, M * N * nbytes));
//    CHECK(deviceMalloc(&d_outputs, M * sizeof(G1Point)));
//
//    CHECK(deviceMemcpy(d_points, h_points.data(), M * N * sizeof(G1AffinePoint), MemcpyKind::HostToDevice));
//    CHECK(deviceMemcpy(d_scalars, h_scalars.data(), M * N * nbytes, MemcpyKind::HostToDevice));
//
//    std::vector<G1Point> reference_outputs(M);
//    bool first_run = true;
//
//    for (int wbits : window_sizes) {
//        size_t workspace_bytes = 0;
//        DeviceError result = blst::msm_booth_batched_pipeline<blst::Fp>(
//            nullptr, workspace_bytes, d_points, d_scalars, nbytes, nbits, N, M, wbits, d_outputs, 0);
//
//        EXPECT_EQ(result,DeviceError::cudaSuccess) << "Window size query failed for wbits=" << wbits;
//
//        void* d_workspace;
//        CHECK(deviceMalloc(&d_workspace, workspace_bytes));
//
//        result = blst::msm_booth_batched_pipeline<blst::Fp>(
//            d_workspace, workspace_bytes, d_points, d_scalars, nbytes, nbits, N, M, wbits, d_outputs, 0);
//
//        EXPECT_EQ(result,DeviceError::cudaSuccess) << "MSM failed for wbits=" << wbits;
//
//        std::vector<G1Point> current_outputs(M);
//        CHECK(deviceMemcpy(current_outputs.data(), d_outputs, M * sizeof(G1Point), MemcpyKind::DeviceToHost));
//
//        if (first_run) {
//            reference_outputs = current_outputs;
//            first_run = false;
//        }
//        else {
//            // All window sizes should produce the same result
//            for (int m = 0; m < M; m++) {
//                EXPECT_TRUE(G1Point::equal(current_outputs[m], reference_outputs[m]))
//                    << "MSM " << m << " differs for wbits=" << wbits;
//            }
//        }
//
//        printf("Window size %d: workspace = %zu bytes\n", wbits, workspace_bytes);
//
//        CHECK(deviceFree(d_workspace));
//    }
//
//    // Cleanup
//    CHECK(deviceFree(d_points));
//    CHECK(deviceFree(d_scalars));
//    CHECK(deviceFree(d_outputs));
//}




