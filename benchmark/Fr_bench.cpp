#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>

// Include the actual BLST field headers
#include "blst/fields/Fr.h"
#include "blst/backend/backend.h"

using namespace blst;

#define LIMBS NLIMBS(256)

// Helper function to set vec256 from uint64_t value
void set_vec256(vec256 v, uint64_t x) {
    for (int i = 0; i < LIMBS; ++i)
        v[i] = 0;
    memcpy(v, &x, sizeof(x));
}

// Helper function to generate random Fr elements
void generate_random_fr(vec256 fr) {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint64_t> dis;

    // Initialize with random values
    for (int i = 0; i < LIMBS; ++i) {
        fr[i] = dis(gen);
    }
}

// ========== ADDITION BENCHMARKS ==========

// Benchmark for Fr addition (ASM version)
static void BM_FrAddition_asm(benchmark::State& state) {
    vec256 a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);
    auto r = get_BLST_r();

    for (auto _ : state) {
        add_mod_256(c, a, b, *r);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrAddition_asm);

// Benchmark for Fr addition (C version)
static void BM_FrAddition_cc(benchmark::State& state) {
    vec256 a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);
    auto r = get_BLST_r();

    for (auto _ : state) {
        cc_add_mod_256(c, a, b, *r);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrAddition_cc);

// Benchmark for Fr addition (High-level API)
static void BM_FrAddition_api(benchmark::State& state) {
    vec256 a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);

    for (auto _ : state) {
        fr_add(c, a, b);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrAddition_api);

// ========== SUBTRACTION BENCHMARKS ==========

// Benchmark for Fr subtraction (ASM version)
static void BM_FrSubtraction_asm(benchmark::State& state) {
    vec256 a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);
    auto r = get_BLST_r();

    for (auto _ : state) {
        sub_mod_256(c, a, b, *r);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSubtraction_asm);

// Benchmark for Fr subtraction (C version)
static void BM_FrSubtraction_cc(benchmark::State& state) {
    vec256 a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);
    auto r = get_BLST_r();

    for (auto _ : state) {
        cc_sub_mod_256(c, a, b, *r);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSubtraction_cc);

// Benchmark for Fr subtraction (High-level API)
static void BM_FrSubtraction_api(benchmark::State& state) {
    vec256 a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);

    for (auto _ : state) {
        fr_sub(c, a, b);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSubtraction_api);

// ========== MULTIPLICATION BENCHMARKS ==========

// Benchmark for Fr multiplication (ASM version)
static void BM_FrMultiplication_asm(benchmark::State& state) {
    vec256 a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);

    fr_to(a, a);
    fr_to(b, b);

    auto r = get_BLST_r();

    for (auto _ : state) {
        mulx_mont_sparse_256(c, a, b, *r, BLST_r0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrMultiplication_asm);

// Benchmark for Fr multiplication (C version)
static void BM_FrMultiplication_cc(benchmark::State& state) {
    vec256 a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);

    fr_to(a, a);
    fr_to(b, b);

    auto r = get_BLST_r();

    for (auto _ : state) {
        cc_mul_mont_sparse_256(c, a, b, *r, BLST_r0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrMultiplication_cc);

// Benchmark for Fr multiplication (High-level API)
static void BM_FrMultiplication_api(benchmark::State& state) {
    vec256 a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);

    fr_to(a, a);
    fr_to(b, b);

    for (auto _ : state) {
        fr_mul(c, a, b);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrMultiplication_api);

// ========== SQUARING BENCHMARKS ==========

// Benchmark for Fr squaring (ASM version)
static void BM_FrSquaring_asm(benchmark::State& state) {
    vec256 a, c;
    generate_random_fr(a);

    fr_to(a, a);

    auto r = get_BLST_r();

    for (auto _ : state) {
        sqrx_mont_sparse_256(c, a, *r, BLST_r0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSquaring_asm);

// Benchmark for Fr squaring (C version)
static void BM_FrSquaring_cc(benchmark::State& state) {
    vec256 a, c;
    generate_random_fr(a);

    fr_to(a, a);

    auto r = get_BLST_r();

    for (auto _ : state) {
        cc_sqr_mont_sparse_256(c, a, *r, BLST_r0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSquaring_cc);

// Benchmark for Fr squaring (High-level API)
static void BM_FrSquaring_api(benchmark::State& state) {
    vec256 a, c;
    generate_random_fr(a);

    fr_to(a, a);

    for (auto _ : state) {
        fr_sqr(c, a);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSquaring_api);

// ========== INVERSION BENCHMARKS ==========

// Benchmark for Fr inversion (High-level API)
static void BM_FrInversion_api(benchmark::State& state) {
    vec256 a, c;
    generate_random_fr(a);
    fr_to(a, a);

    for (auto _ : state) {
        fr_inverse(c, a);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrInversion_api);

// ========== MONTGOMERY CONVERSION BENCHMARKS ==========

// Benchmark for Montgomery conversion to (High-level API)
static void BM_FrToMontgomery_api(benchmark::State& state) {
    vec256 a, c;
    generate_random_fr(a);

    for (auto _ : state) {
        fr_to(c, a);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrToMontgomery_api);

// Benchmark for Montgomery conversion from (ASM version)
static void BM_FrFromMontgomery_asm(benchmark::State& state) {
    vec256 a, c;
    generate_random_fr(a);
    fr_to(a, a);
    auto r = get_BLST_r();

    for (auto _ : state) {
        fromx_mont_256(c, a, *r, BLST_r0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrFromMontgomery_asm);

// Benchmark for Montgomery conversion from (C version)
static void BM_FrFromMontgomery_cc(benchmark::State& state) {
    vec256 a, c;
    generate_random_fr(a);
    fr_to(a, a);
    auto r = get_BLST_r();

    for (auto _ : state) {
        cc_from_mont_256(c, a, *r, BLST_r0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrFromMontgomery_cc);

// Benchmark for Montgomery conversion from (High-level API)
static void BM_FrFromMontgomery_api(benchmark::State& state) {
    vec256 a, c;
    generate_random_fr(a);
    fr_to(a, a);

    for (auto _ : state) {
        fr_from(c, a);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrFromMontgomery_api);