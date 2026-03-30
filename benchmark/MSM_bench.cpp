#include <benchmark/benchmark.h>
#include <vector>
#include <chrono>
#include <iostream>
#include <string>


#include "blst/curve/multi_scalar_par.h"
#include "test/InsecureRng.h"
#include "test/device.h"

#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)

using namespace blst;

// Template for different group types
template<typename Group, typename AffineGroup, typename GroupZZ, Backend B>
class MSMBenchmark {
private:
    std::vector<AffineGroup> h_points;
    std::vector<unsigned char> h_scalars;
    std::vector<Group> h_cuda_outputs, h_cpu_outputs;

    // Device memory
    AffineGroup* d_points_ptr = nullptr;
    uint8_t* d_scalars_ptr = nullptr;
    Group* d_outputs_ptr = nullptr;
    void* d_workspace = nullptr;

    int N, M, nbits, wbits, nbytes;
    size_t workspace_bytes = 0;

public:
    MSMBenchmark(int n_points, int n_msms, int num_bits = 256, int window_bits = 3)
        : N(n_points), M(n_msms), nbits(num_bits), wbits(window_bits) {
        nbytes = (nbits + 7) / 8;
        setupData();
        setupDevice();
    }

    ~MSMBenchmark() {
        cleanup();
    }

    void setupData() {
        // Allocate host memory
        h_points.resize(M * N);
        h_scalars.resize(M * N * nbytes);
        h_cuda_outputs.resize(M);
        h_cpu_outputs.resize(M);

        // Generate test data
        InsecureRng rng;

        // Generate M*N unique points
        std::vector<Group> temp_points(M * N);
        create_test_group(temp_points.data(), M * N, rng);
        for (int i = 0; i < M * N; i++) {
            msm_to_affine(h_points[i], temp_points[i]);
        }

        // Generate random scalars
        for (int m = 0; m < M; m++) {
            for (int n = 0; n < N; n++) {
                for (int b = 0; b < nbytes; b++) {
                    h_scalars[m * N * nbytes + n * nbytes + b] = rng.uniform_u32(0, 255);
                }
            }
        }
    }

    void setupDevice() {
        // Allocate device memory
        CHECK(deviceMalloc<B>(&d_points_ptr, M * N * sizeof(AffineGroup)));
        CHECK(deviceMalloc<B>(&d_scalars_ptr, M * N * nbytes));
        CHECK(deviceMalloc<B>(&d_outputs_ptr, M * sizeof(Group)));

        // Copy data to device
        CHECK(deviceMemcpy<B>(d_points_ptr, h_points.data(), M * N * sizeof(AffineGroup), MemcpyKind::HostToDevice));
        CHECK(deviceMemcpy<B>(d_scalars_ptr, h_scalars.data(), M * N * nbytes, MemcpyKind::HostToDevice));

        // Query workspace size
        d2Span<AffineGroup> d_points(d_points_ptr, M, N);
        d3Span<uint8_t> d_scalars(d_scalars_ptr, M, N, nbytes);
        span<Group> d_outputs(d_outputs_ptr, M);

        auto result = blst::msm_booth_batched_pipeline<B, Group, AffineGroup, GroupZZ>(
            nullptr, workspace_bytes, d_points, d_scalars, nbytes, nbits, N, M, wbits, d_outputs, 0, false);

        if (result == DeviceError::cudaSuccess && workspace_bytes > 0) {
            CHECK(deviceMalloc<B>(&d_workspace, workspace_bytes));
        }
    }

    void benchmarkCPU() {
        for (int m = 0; m < M; m++) {
            std::vector<const AffineGroup*> points_for_msm(1);
            std::vector<const byte*> scalars_for_msm(1);

            points_for_msm[0] = &h_points[m * N];
            scalars_for_msm[0] = &h_scalars[m * N * nbytes];

            size_t scratch_size = blst::msm_mult_pippenger_scratch_sizeof_<GroupZZ>(N);
            GroupZZ* scratch = (GroupZZ*)malloc(scratch_size);

            blst::msm_mult_pippenger<Group, AffineGroup, GroupZZ>(
                h_cpu_outputs[m], points_for_msm.data(), N, scalars_for_msm.data(),
                nbits, scratch, wbits, false);

            free(scratch);
        }
    }

    void benchmarkGPU() {
        d2Span<AffineGroup> d_points(d_points_ptr, M, N);
        d3Span<uint8_t> d_scalars(d_scalars_ptr, M, N, nbytes);
        span<Group> d_outputs(d_outputs_ptr, M);

        auto result = blst::msm_booth_batched_pipeline<B, Group, AffineGroup, GroupZZ>(
            d_workspace, workspace_bytes, d_points, d_scalars, nbytes, nbits, N, M, wbits, d_outputs, 0, false, false);

        CHECK(streamSynchronize<B>(0));
    }

    void cleanup() {
        if (d_points_ptr) CHECK(deviceFree<B>(d_points_ptr));
        if (d_scalars_ptr) CHECK(deviceFree<B>(d_scalars_ptr));
        if (d_outputs_ptr) CHECK(deviceFree<B>(d_outputs_ptr));
        if (d_workspace) CHECK(deviceFree<B>(d_workspace));
    }

    int getTotalPoints() const { return M * N; }
    int getBatchCount() const { return M; }
    int getPointsPerMSM() const { return N; }
};

// Benchmark function for CPU MSM
template<typename Group, typename AffineGroup, typename GroupZZ>
static void BM_MSM_CPU(benchmark::State& state) {
    int points_per_msm = state.range(0);
    int num_msms = state.range(1);

    MSMBenchmark<Group, AffineGroup, GroupZZ, Backend::CPU> bench(points_per_msm, num_msms);

    for (auto _ : state) {
        bench.benchmarkCPU();
    }

    // Set custom counters
    state.counters["MSMs"] = num_msms;
    state.counters["Points/MSM"] = points_per_msm;
    state.counters["Total Points"] = bench.getTotalPoints();
    state.counters["Points/sec"] = benchmark::Counter(bench.getTotalPoints(), benchmark::Counter::kIsRate);
}

// Benchmark function for GPU MSM
template<typename Group, typename AffineGroup, typename GroupZZ>
static void BM_MSM_GPU(benchmark::State& state) {
    int points_per_msm = state.range(0);
    int num_msms = state.range(1);

    MSMBenchmark<Group, AffineGroup, GroupZZ, Backend::CUDA> bench(points_per_msm, num_msms);

    for (auto _ : state) {
        bench.benchmarkGPU();
    }

    // Set custom counters
    state.counters["MSMs"] = num_msms;
    state.counters["Points/MSM"] = points_per_msm;
    state.counters["Total Points"] = bench.getTotalPoints();
    state.counters["Points/sec"] = benchmark::Counter(bench.getTotalPoints(), benchmark::Counter::kIsRate);
}

// Register benchmarks for different group types and configurations

#define PARAMS \
->Args({ 64, 1 })->Args({ 128, 1 })->Args({ 256, 1 })->Args({ 512, 1 })->Args({ 1028, 1 }) \
->Args({ 2048, 1 })->Args({ 4096, 1 })->Args({ 8192, 1 })->Args({ 10324, 1 })->Args({ 20648, 1 })

//->Args({ 64, 1 })->Args({ 128, 128 })->Args({ 256, 256 })->Args({ 512, 512 })->Args({ 1028, 1024 })

// G1 benchmarks
BENCHMARK_TEMPLATE(BM_MSM_CPU, G1Point, G1AffinePoint, G1PointZZ)
PARAMS
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK_TEMPLATE(BM_MSM_GPU, G1Point, G1AffinePoint, G1PointZZ)
PARAMS
->Unit(benchmark::kMillisecond)
->UseRealTime();

// G2 benchmarks
BENCHMARK_TEMPLATE(BM_MSM_CPU, G2Point, G2AffinePoint, G2PointZZ)
PARAMS
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK_TEMPLATE(BM_MSM_GPU, G2Point, G2AffinePoint, G2PointZZ)
PARAMS
->Unit(benchmark::kMillisecond)
->UseRealTime();

// GT benchmarks (if applicable)

BENCHMARK_TEMPLATE(BM_MSM_CPU, GT, GT, GT)
PARAMS
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK_TEMPLATE(BM_MSM_GPU, GT, GT, GT)
PARAMS
->Unit(benchmark::kMillisecond)
->UseRealTime();

// Scaling benchmarks - fixed number of MSMs, varying points per MSM
static void BM_MSM_GPU_G1_Scaling(benchmark::State& state) {
    int points_per_msm = state.range(0);
    int num_msms = 4; // Fixed number of MSMs

    MSMBenchmark<G1Point, G1AffinePoint, G1PointZZ, Backend::CUDA> bench(points_per_msm, num_msms);

    for (auto _ : state) {
        bench.benchmarkGPU();
    }

    state.counters["Points/MSM"] = points_per_msm;
    state.counters["Total Points"] = bench.getTotalPoints();
    state.counters["Points/sec"] = benchmark::Counter(bench.getTotalPoints(), benchmark::Counter::kIsRate);
}

BENCHMARK(BM_MSM_GPU_G1_Scaling)
->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)->Arg(8192)->Arg(16384)
->Unit(benchmark::kMillisecond)
->UseRealTime();

// Batch scaling benchmarks - fixed points per MSM, varying number of MSMs
static void BM_MSM_GPU_G1_BatchScaling(benchmark::State& state) {
    int points_per_msm = 2048; // Fixed points per MSM
    int num_msms = state.range(0);

    MSMBenchmark<G1Point, G1AffinePoint, G1PointZZ, Backend::CUDA> bench(points_per_msm, num_msms);

    for (auto _ : state) {
        bench.benchmarkGPU();
    }

    state.counters["MSMs"] = num_msms;
    state.counters["Total Points"] = bench.getTotalPoints();
    state.counters["Points/sec"] = benchmark::Counter(bench.getTotalPoints(), benchmark::Counter::kIsRate);
}

BENCHMARK(BM_MSM_GPU_G1_BatchScaling)
->Arg(1)->Arg(2)->Arg(4)->Arg(8)->Arg(16)->Arg(32)
->Unit(benchmark::kMillisecond)
->UseRealTime();

#endif