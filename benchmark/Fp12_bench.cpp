#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>

// Include the actual BLST field headers
#include "blst/fields/Fp12.h"

using namespace blst;

// Helper function to generate random Fp12 elements
void generate_random_fp12(Fp12& fp12) {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint64_t> dis;

    limb* p = (limb*)&fp12;
    for (size_t i = 0; i < sizeof(Fp12) / sizeof(limb); ++i) {
        p[i] = dis(gen);
    }
}

// ========== MULTIPLICATION BENCHMARKS ==========

// Benchmark for Fp12 multiplication
static void BM_Fp12Multiplication(benchmark::State& state) {
    Fp12 a, b, c;
    generate_random_fp12(a);
    generate_random_fp12(b);

    for (auto _ : state) {
        blst_fp12_mul(c, a, b);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Multiplication);

// ========== SQUARING BENCHMARKS ==========

// Benchmark for Fp12 squaring
static void BM_Fp12Squaring(benchmark::State& state) {
    Fp12 a, c;
    generate_random_fp12(a);

    for (auto _ : state) {
        blst_fp12_sqr(c, a);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Squaring);

// Benchmark for Fp12 cyclotomic squaring
static void BM_Fp12CyclotomicSquaring(benchmark::State& state) {
    Fp12 a, c;
    generate_random_fp12(a);

    for (auto _ : state) {
        blst_fp12_cyclotomic_sqr(c, a);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12CyclotomicSquaring);

// ========== INVERSION BENCHMARK ==========

// Benchmark for Fp12 inversion
static void BM_Fp12Inversion(benchmark::State& state) {
    Fp12 a, c;
    generate_random_fp12(a);

    for (auto _ : state) {
        blst_fp12_inverse(c, a);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Inversion);

// ========== CONJUGATION BENCHMARK ==========

// Benchmark for Fp12 conjugation
static void BM_Fp12Conjugate(benchmark::State& state) {
    Fp12 a;
    generate_random_fp12(a);

    for (auto _ : state) {
        // Conjugation is in-place
        blst_fp12_conjugate(a);
        benchmark::DoNotOptimize(a);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Conjugate);

// ========== FROBENIUS MAP BENCHMARK ==========

// Benchmark for Fp12 Frobenius map
static void BM_Fp12Frobenius(benchmark::State& state) {
    Fp12 a, c;
    generate_random_fp12(a);

    for (auto _ : state) {
        blst_fp12_frobenius_map(c, a, 1);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Frobenius);

// ========== EXPONENTIATION BENCHMARK ==========

// Benchmark for Fp12 exponentiation
static void BM_Fp12Exponentiation(benchmark::State& state) {
    Fp12 a, c;
    generate_random_fp12(a);

    // A small, simple exponent for benchmarking purposes
    limb exponent[] = { 0x7, 0x0 };

    for (auto _ : state) {
        blst_fp12_pow(c, a, exponent, 3);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Exponentiation);

