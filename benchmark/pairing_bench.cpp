#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>

#include "blst/curve/point.h"
#include "blst/curve/e1.h"
#include "blst/curve/e2.h"
#include "blst/curve/pairing.h"
#include "blst/curve/GT_cuda.h"


using namespace blst;

// Helper to generate a random 256-bit scalar
static void generate_random_scalar(byte* scalar) {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint64_t> dis;
    uint64_t* as_limbs = reinterpret_cast<uint64_t*>(scalar);
    for (size_t i = 0; i < 256 / 64; ++i) {
        as_limbs[i] = dis(gen);
    }
}

// Helper to generate a random G1 point
static void generate_random_g1_point(POINTonE1& p) {
    byte scalar[32];
    generate_random_scalar(scalar);
    blst_p1_mult(p, blst_p1_generator(), scalar, 256);
}

// Helper to generate a random G2 point
static void generate_random_g2_point(POINTonE2& p) {
    byte scalar[32];
    generate_random_scalar(scalar);
    blst_p2_mult(p, blst_p2_generator(), scalar, 256);
}

// ========== PAIRING BENCHMARKS ==========

// Benchmark for Miller loop (single pairing)
static void BM_MillerLoop(benchmark::State& state) {
    POINTonE1 g1;
    POINTonE2 g2;
    generate_random_g1_point(g1);
    generate_random_g2_point(g2);

    POINTonE1_affine g1_affine;
    POINTonE2_affine g2_affine;
    blst_p1_to_affine(g1_affine, g1);
    blst_p2_to_affine(g2_affine, g2);

    Fp12 result;
    for (auto _ : state) {
        blst_miller_loop(result, &g2_affine, &g1_affine);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_MillerLoop);

// Benchmark for final exponentiation
static void BM_FinalExp(benchmark::State& state) {
    POINTonE1 g1;
    POINTonE2 g2;
    generate_random_g1_point(g1);
    generate_random_g2_point(g2);

    POINTonE1_affine g1_affine;
    POINTonE2_affine g2_affine;
    blst_p1_to_affine(g1_affine, g1);
    blst_p2_to_affine(g2_affine, g2);

    Fp12 result;
    blst_miller_loop(result, &g2_affine, &g1_affine);

    for (auto _ : state) {
        blst_final_exp(result, result);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FinalExp);

// Benchmark for multi-pairing (n=2)
static void BM_MultiPairing2(benchmark::State& state) {
    POINTonE1 g1[2];
    POINTonE2 g2[2];
    for (int i = 0; i < 2; ++i) {
        generate_random_g1_point(g1[i]);
        generate_random_g2_point(g2[i]);
    }

    POINTonE1_affine g1_affine[2];
    POINTonE2_affine g2_affine[2];
    for (int i = 0; i < 2; ++i) {
        blst_p1_to_affine(g1_affine[i], g1[i]);
        blst_p2_to_affine(g2_affine[i], g2[i]);
    }

    // Arrays of pointers for API
    const POINTonE1_affine* g1_ptrs[2] = { &g1_affine[0], &g1_affine[1] };
    const POINTonE2_affine* g2_ptrs[2] = { &g2_affine[0], &g2_affine[1] };

    Fp12 result;
    for (auto _ : state) {
        blst_miller_loop_n(result, g2_ptrs, g1_ptrs, 2);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_MultiPairing2);

// Benchmark for multi-pairing (n=4)
static void BM_MultiPairing4(benchmark::State& state) {
    POINTonE1 g1[4];
    POINTonE2 g2[4];
    for (int i = 0; i < 4; ++i) {
        generate_random_g1_point(g1[i]);
        generate_random_g2_point(g2[i]);
    }

    POINTonE1_affine g1_affine[4];
    POINTonE2_affine g2_affine[4];
    for (int i = 0; i < 4; ++i) {
        blst_p1_to_affine(g1_affine[i], g1[i]);
        blst_p2_to_affine(g2_affine[i], g2[i]);
    }

    // Arrays of pointers for API
    const POINTonE1_affine* g1_ptrs[4] = { &g1_affine[0], &g1_affine[1], &g1_affine[2], &g1_affine[3] };
    const POINTonE2_affine* g2_ptrs[4] = { &g2_affine[0], &g2_affine[1], &g2_affine[2], &g2_affine[3] };

    Fp12 result;
    for (auto _ : state) {
        blst_miller_loop_n(result, g2_ptrs, g1_ptrs, 4);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_MultiPairing4);

// Benchmark for pairing with identity points
static void BM_PairingWithIdentity(benchmark::State& state) {
    POINTonE1 identity_e1;
    POINTonE2 identity_e2;
    memset(&identity_e1, 0, sizeof(identity_e1));
    memset(&identity_e2, 0, sizeof(identity_e2));

    POINTonE1_affine identity_e1_affine;
    POINTonE2_affine identity_e2_affine;
    blst_p1_to_affine(identity_e1_affine, identity_e1);
    blst_p2_to_affine(identity_e2_affine, identity_e2);

    Fp12 result;
    for (auto _ : state) {
        blst_miller_loop(result, &identity_e2_affine, &identity_e1_affine);
        blst_final_exp(result, result);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_PairingWithIdentity);


// Benchmark for batched Miller loop + final exponentiation (CPU)
static void BM_CPU_MillerFinal(benchmark::State& state) {
    const size_t num_pairs = static_cast<size_t>(state.range(0));
    auto g1_points = new POINTonE1_affine[num_pairs];
    auto g2_points = new POINTonE2_affine[num_pairs];
    auto miller_results = new Fp12[num_pairs];
    auto final_results = new Fp12[num_pairs];

    // Generate random G1/G2 pairs
    for (size_t i = 0; i < num_pairs; ++i) {
        POINTonE1 g1;
        POINTonE2 g2;
        byte scalar1[32], scalar2[32];
        generate_random_scalar(scalar1);
        generate_random_scalar(scalar2);
        blst_p1_mult(g1, blst_p1_generator(), scalar1, 256);
        blst_p2_mult(g2, blst_p2_generator(), scalar2, 256);
        blst_p1_to_affine(g1_points[i], g1);
        blst_p2_to_affine(g2_points[i], g2);
    }

    for (auto _ : state) {
        // Miller loop (batched, CPU)
        for (size_t i = 0; i < num_pairs; ++i) {
            blst_miller_loop(miller_results[i], &g2_points[i], &g1_points[i]);
        }
        // Final exponentiation (batched, CPU)
        for (size_t i = 0; i < num_pairs; ++i) {
            blst_final_exp(final_results[i], miller_results[i]);
        }
        benchmark::DoNotOptimize(final_results);
    }
    state.SetItemsProcessed(state.iterations() * num_pairs);

    delete[] g1_points;
    delete[] g2_points;
    delete[] miller_results;
    delete[] final_results;
}
BENCHMARK(BM_CPU_MillerFinal)->Arg(64)->Arg(256)->Arg(1024);


// ========== CUDA PAIRING BENCHMARKS ==========

#ifdef BLST_CUDA
static void BM_CUDA_MillerFinal(benchmark::State& state) {
    const size_t num_pairs = static_cast<size_t>(state.range(0));
    auto h_g1_points = new POINTonE1_affine[num_pairs];
    auto h_g2_points = new POINTonE2_affine[num_pairs];
    auto h_miller_results = new Fp12[num_pairs];
    auto h_final_results = new Fp12[num_pairs];

    // Generate random G1/G2 pairs
    for (size_t i = 0; i < num_pairs; ++i) {
        POINTonE1 g1;
        POINTonE2 g2;
        byte scalar1[32], scalar2[32];
        generate_random_scalar(scalar1);
        generate_random_scalar(scalar2);
        blst_p1_mult(g1, blst_p1_generator(), scalar1, 256);
        blst_p2_mult(g2, blst_p2_generator(), scalar2, 256);
        blst_p1_to_affine(h_g1_points[i], g1);
        blst_p2_to_affine(h_g2_points[i], g2);
    }

    // Allocate device memory
    POINTonE1_affine* d_g1_points = nullptr;
    POINTonE2_affine* d_g2_points = nullptr;
    Fp12* d_miller_results = nullptr;
    Fp12* d_final_results = nullptr;

    cudaMalloc(&d_g1_points, num_pairs * sizeof(POINTonE1_affine));
    cudaMalloc(&d_g2_points, num_pairs * sizeof(POINTonE2_affine));
    cudaMalloc(&d_miller_results, num_pairs * sizeof(Fp12));
    cudaMalloc(&d_final_results, num_pairs * sizeof(Fp12));

    cudaMemcpy(d_g1_points, h_g1_points, num_pairs * sizeof(POINTonE1_affine), cudaMemcpyHostToDevice);
    cudaMemcpy(d_g2_points, h_g2_points, num_pairs * sizeof(POINTonE2_affine), cudaMemcpyHostToDevice);

    for (auto _ : state) {
        // Miller loop (batched)
        batched_error_code miller_result = miller_loop_gt_vec(d_miller_results, d_g2_points, d_g1_points, num_pairs);
        benchmark::DoNotOptimize(miller_result);

        // Final exponentiation (batched)
        batched_error_code final_result = final_exp_gt_vec(d_final_results, d_miller_results, num_pairs);
        benchmark::DoNotOptimize(final_result);

        // Optionally, synchronize to ensure all work is done
        //cudaDeviceSynchronize();
    }
    state.SetItemsProcessed(state.iterations() * num_pairs);

    // Cleanup
    cudaFree(d_g1_points);
    cudaFree(d_g2_points);
    cudaFree(d_miller_results);
    cudaFree(d_final_results);
	delete[] h_g1_points;
	delete[] h_g2_points;
	delete[] h_miller_results;
	delete[] h_final_results;
}

BENCHMARK(BM_CUDA_MillerFinal)->Arg(64)->Arg(256)->Arg(1024)->Arg(4096)->Arg(16384)->Arg(65536);
#endif