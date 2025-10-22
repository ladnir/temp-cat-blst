
#include <gtest/gtest.h>

#include "blst/curve/multi_scalar_par.h"
#include "InsecureRng.h"
using namespace blst;


// CUDA kernel wrapper macros for device execution
#define CHECK(call) do { \
    DeviceError err = call; \
    if (err != DeviceError::cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, getErrorString(err)); \
        exit(1); \
    } \
} while(0)


namespace 
{
    // Helper function to create test points on host
    void create_test_points(G1Point* points, size_t num, InsecureRng& rng) {
        G1Point generator = G1Point::generator();
        unsigned char scalar[32] = { 0 };

        for (size_t i = 0; i < num; i++) {
            scalar[0] = rng.uniform_u64();
            G1Point::mult(points[i], generator, scalar, 8);
        }
    }
}

template<typename T>
auto deviceMalloc(T** ptr, size_t size) -> DeviceError {
#ifdef BLST_CUDA

    cudaError_t err = cudaMalloc(ptr, size);
    return static_cast<DeviceError>(err);
#else
    *ptr = (T*)malloc(size);
    return DeviceError::cudaSuccess;
#endif // BLST_CUDA
}

DeviceError deviceMemcpy(void* dest, const void* src, size_t size, MemcpyKind kind)
{
#ifdef BLST_CUDA
    return cudaMemcpy(dest, src, size, static_cast<cudaMemcpyKind>(kind));
#else
    memcpy(dest, src, size);
    return DeviceError::cudaSuccess;
#endif // BLST_CUDA
}

auto deviceFree(void* ptr) -> DeviceError
{
#ifdef BLST_CUDA
    cudaError_t err = cudaFree(ptr);
    return static_cast<DeviceError>(err);
#else
    free(ptr);
    return DeviceError::cudaSuccess;
#endif // BLST_CUDA
}



TEST(Msm, BoothBatchedPipeline) {
    //if (!G1CUDAEnvironment::IsCudaAvailable()) {
    //    GTEST_SKIP() << "CUDA not available";
    //}

    bool debug = false;
    bool print = false;

    // Test parameters
    const int N = 16;        // points per MSM
    const int M = 1;         // number of MSMs
    const int nbits = 6;   // scalar bit length
    const int nbytes = (nbits + 7) / 8;
    const int wbits = 3;     // window size

    // Host data - now M*N points total, N points per MSM
    std::vector<G1AffinePoint> h_points(M * N);  // M*N total points
    std::vector<unsigned char> h_scalars(M * N * nbytes);
    std::vector<G1Point> h_cuda_outputs(M), h_cpu_outputs(M);

    // Create test data
    InsecureRng rng;

    // Generate M*N unique points (each MSM gets its own set of N points)
    std::vector<G1Point> temp_points(M * N);
    create_test_points(temp_points.data(), M * N, rng);
    for (int i = 0; i < M * N; i++) {
        G1Point::to_affine(h_points[i], temp_points[i]);
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
        std::vector<const G1AffinePoint*> points_for_msm(2);
        std::vector<const byte*> scalars_for_msm(2);

        points_for_msm[0] = &h_points[m * N];
        scalars_for_msm[0] = &h_scalars[m * N * nbytes];

        // Each MSM m uses points [m*N ... (m+1)*N-1]
        //for (int n = 0; n < N; n++) {
        //    points_for_msm[n] = &h_points[m * N + n];  // MSM m uses points starting at m*N
        //    scalars_for_msm[n] = &h_scalars[m * N * nbytes + n * nbytes];
        //}

        // Use CPU Pippenger for reference
        size_t scratch_size = blst::msm_mult_pippenger_scratch_sizeof<blst::Fp>(N);
        blst::Pointxyzz<blst::Fp>* scratch =
            (blst::Pointxyzz<blst::Fp>*)malloc(scratch_size);

        blst::msm_mult_pippenger<blst::Fp>(
            h_cpu_outputs[m],
            points_for_msm.data(),
            N,
            scalars_for_msm.data(),
            nbits,
            scratch,
            wbits
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
    G1AffinePoint* d_points;
    uint8_t* d_scalars;
    G1Point* d_outputs;

    CHECK(deviceMalloc(&d_points, M * N * sizeof(G1AffinePoint)));  // M*N points
    CHECK(deviceMalloc(&d_scalars, M * N * nbytes));
    CHECK(deviceMalloc(&d_outputs, M * sizeof(G1Point)));

    // Copy data to device - all M*N points
    CHECK(deviceMemcpy(d_points, h_points.data(), M * N * sizeof(G1AffinePoint), MemcpyKind::HostToDevice));
    CHECK(deviceMemcpy(d_scalars, h_scalars.data(), M * N * nbytes, MemcpyKind::HostToDevice));

    // Query workspace size
    size_t workspace_bytes = 0;
    auto result = blst::msm_booth_batched_pipeline<blst::Fp>(
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
    CHECK(deviceMalloc(&d_workspace, workspace_bytes));

    // Run the batched MSM pipeline
    result = blst::msm_booth_batched_pipeline<blst::Fp>(
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
    CHECK(streamSynchronize(0));

    // Copy results back
    CHECK(deviceMemcpy(h_cuda_outputs.data(), d_outputs, M * sizeof(G1Point), MemcpyKind::DeviceToHost));

    // Compare CUDA results with CPU results
    bool all_match = true;
    for (int m = 0; m < M; m++) {
        if (!G1Point::equal(h_cuda_outputs[m], h_cpu_outputs[m])) {
            printf("MSM %d mismatch\n", m);
            all_match = false;
        }
    }

    EXPECT_TRUE(all_match) << "CUDA and CPU MSM results should match";

    // Cleanup
    CHECK(deviceFree(d_points));
    CHECK(deviceFree(d_scalars));
    CHECK(deviceFree(d_outputs));
    CHECK(deviceFree(d_workspace));
}

//TEST(G1OpsCUDA, MSMBoothBatchedPipelineLargeScale) {
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
//    create_test_points(temp_points.data(), M * N, rng);
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
//    cudaError_t result = blst::msm_booth_batched_pipeline<blst::Fp>(
//        nullptr, workspace_bytes, d_points, d_scalars, nbytes, nbits, N, M, wbits, d_outputs, 0);
//
//    EXPECT_EQ(result, cudaSuccess);
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
//    EXPECT_EQ(result, cudaSuccess);
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

//TEST(G1OpsCUDA, MSMBoothBatchedPipelineWindowSizes) {
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
//    create_test_points(temp_points.data(), M * N, rng);
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
//        cudaError_t result = blst::msm_booth_batched_pipeline<blst::Fp>(
//            nullptr, workspace_bytes, d_points, d_scalars, nbytes, nbits, N, M, wbits, d_outputs, 0);
//
//        EXPECT_EQ(result, cudaSuccess) << "Window size query failed for wbits=" << wbits;
//
//        void* d_workspace;
//        CHECK(deviceMalloc(&d_workspace, workspace_bytes));
//
//        result = blst::msm_booth_batched_pipeline<blst::Fp>(
//            d_workspace, workspace_bytes, d_points, d_scalars, nbytes, nbits, N, M, wbits, d_outputs, 0);
//
//        EXPECT_EQ(result, cudaSuccess) << "MSM failed for wbits=" << wbits;
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




