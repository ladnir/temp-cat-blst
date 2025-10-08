#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>

// Include the actual BLST field headers
#include "blst/fields/Fr.h"
#include "blst/backend/backend.h"

using namespace blst;

#define LIMBS NLIMBS(256)

// Helper function to set Fr& from uint64_t value
void set_Fr(Fr& v, uint64_t x) {
    for (int i = 0; i < LIMBS; ++i)
        v[i] = 0;
    memcpy(v, &x, sizeof(x));
}

// Helper function to generate random Fr elements
void generate_random_fr(Fr& fr) {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint64_t> dis;

    // Initialize with random values
    for (int i = 0; i < LIMBS; ++i) {
        fr[i] = dis(gen);
    }

    while (Fr::is_redc(fr) == 0)
    {
        Fr::rshift(fr, fr, 1);
    }
}

// ========== ADDITION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fr addition (ASM version)
static void BM_FrAddition_asm(benchmark::State& state) {
    Fr a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);
    auto r = get_BLST_r();

    for (auto _ : state) {
        add_mod_256(&c, &a, &b, &r);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrAddition_asm);
#endif

// Benchmark for Fr addition (C version)
static void BM_FrAddition_cc(benchmark::State& state) {
    Fr a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);
    auto r = get_BLST_r();

    for (auto _ : state) {
        //cc_add_mod_256(c, a, b, r);
        cc_add_mod_n<NLIMBS(256)>(c.data, a.data, b.data, r.data);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrAddition_cc);

// Benchmark for Fr addition (High-level API)
static void BM_FrAddition_api(benchmark::State& state) {
    Fr a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);

    for (auto _ : state) {
        Fr::add(c, a, b);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrAddition_api);

// ========== SUBTRACTION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fr subtraction (ASM version)
static void BM_FrSubtraction_asm(benchmark::State& state) {
    Fr a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);
    auto r = get_BLST_r();

    for (auto _ : state) {
        sub_mod_256(&c, &a, &b, &r);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSubtraction_asm);
#endif

// Benchmark for Fr subtraction (C version)
static void BM_FrSubtraction_cc(benchmark::State& state) {
    Fr a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);
    auto r = get_BLST_r();

    for (auto _ : state) {
        //cc_sub_mod_256(c, a, b, r);
        cc_sub_mod_n<NLIMBS(256)>(c.data, a.data, b.data, r.data);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSubtraction_cc);

// Benchmark for Fr subtraction (High-level API)
static void BM_FrSubtraction_api(benchmark::State& state) {
    Fr a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);

    for (auto _ : state) {
        Fr::sub(c, a, b);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSubtraction_api);

// ========== MULTIPLICATION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fr multiplication (ASM version)
static void BM_FrMultiplication_asm(benchmark::State& state) {
    Fr a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);

    Fr::to_mont(a, a);
    Fr::to_mont(b, b);

    auto r = get_BLST_r();

    for (auto _ : state) {
        mulx_mont_sparse_256(&c, &a, &b, &r, BLST_r0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrMultiplication_asm);
#endif

// Benchmark for Fr multiplication (C version)
static void BM_FrMultiplication_cc(benchmark::State& state) {
    Fr a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);

    Fr::to_mont(a, a);
    Fr::to_mont(b, b);

    auto r = get_BLST_r();

    for (auto _ : state) {
        //cc_mul_mont_sparse_256(c, a, b, r, BLST_r0);
        cc_mul_mont_n<NLIMBS(256)>(c.data, a.data, b.data, r.data, BLST_r0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrMultiplication_cc);

// Benchmark for Fr multiplication (High-level API)
static void BM_FrMultiplication_api(benchmark::State& state) {
    Fr a, b, c;
    generate_random_fr(a);
    generate_random_fr(b);

    Fr::to_mont(a, a);
    Fr::to_mont(b, b);

    for (auto _ : state) {
        Fr::mul(c, a, b);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrMultiplication_api);

// ========== SQUARING BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fr squaring (ASM version)
static void BM_FrSquaring_asm(benchmark::State& state) {
    Fr a, c;
    generate_random_fr(a);

    Fr::to_mont(a, a);

    auto r = get_BLST_r();

    for (auto _ : state) {
        sqrx_mont_sparse_256(&c, &a, &r, BLST_r0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSquaring_asm);
#endif

// Benchmark for Fr squaring (C version)
static void BM_FrSquaring_cc(benchmark::State& state) {
    Fr a, c;
    generate_random_fr(a);

    Fr::to_mont(a, a);

    auto r = get_BLST_r();

    for (auto _ : state) {
        //cc_sqr_mont_sparse_256(c, a, r, BLST_r0);
        cc_mul_mont_n<NLIMBS(256)>(c.data, a.data, a.data, r.data, BLST_r0);

        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSquaring_cc);

// Benchmark for Fr squaring (High-level API)
static void BM_FrSquaring_api(benchmark::State& state) {
    Fr a, c;
    generate_random_fr(a);

    Fr::to_mont(a, a);

    for (auto _ : state) {
        Fr::sqr(c, a);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrSquaring_api);

// ========== INVERSION BENCHMARKS ==========


// Benchmark for Fr inversion (High-level API)
static void BM_FrInversion_api(benchmark::State& state) {
    Fr a, c;
    generate_random_fr(a);
    Fr::to_mont(a, a);

    for (auto _ : state) {
        Fr::reciprocal(c, a);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrInversion_api);

// ========== MONTGOMERY CONVERSION BENCHMARKS ==========

// Benchmark for Montgomery conversion to (High-level API)
static void BM_FrToMontgomery_api(benchmark::State& state) {
    Fr a, c;
    generate_random_fr(a);

    for (auto _ : state) {
        Fr::to_mont(c, a);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrToMontgomery_api);

#ifdef BLST_ASM

// Benchmark for Montgomery conversion from (ASM version)
static void BM_FrFromMontgomery_asm(benchmark::State& state) {
    Fr a, c;
    generate_random_fr(a);
    Fr::to_mont(a, a);
    auto r = get_BLST_r();

    for (auto _ : state) {
        fromx_mont_256(&c, &a, &r, BLST_r0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrFromMontgomery_asm);
#endif

// Benchmark for Montgomery conversion from (C version)
static void BM_FrFromMontgomery_cc(benchmark::State& state) {
    Fr a, c;
    generate_random_fr(a);
    Fr::to_mont(a, a);
    auto r = get_BLST_r();

    for (auto _ : state) {
        //cc_from_mont_256(c, a, r, BLST_r0);
        cc_from_mont_n<NLIMBS(256)>(c.data, a.data, r.data, BLST_r0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrFromMontgomery_cc);

// Benchmark for Montgomery conversion from (High-level API)
static void BM_FrFromMontgomery_api(benchmark::State& state) {
    Fr a, c;
    generate_random_fr(a);
    Fr::to_mont(a, a);

    for (auto _ : state) {
        Fr::from_mont(c, a);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FrFromMontgomery_api);