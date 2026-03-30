#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>

// Include the actual BLST field headers
#include "blst/fields/Fp2.h"
#include "blst/backend/backend.h"

using namespace blst;

#define LIMBS NLIMBS(384)

// Helper function to generate random Fp2 elements
void generate_random_fp2(Fp2& fp2) {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint64_t> dis;

    // Initialize with random values
    for (int i = 0; i < LIMBS; ++i) {
        fp2[0][i] = dis(gen);
        fp2[1][i] = dis(gen);
    }

    for (auto& v : fp2.data)
    {
        while (Fp::is_redc(v) == 0)
        {
            Fp::rshift(v, v, 1);
        }
    }
}

// ========== ADDITION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp2 addition (ASM version)
static void BM_Fp2Addition_asm(benchmark::State& state) {
    Fp2 a, b, c;
    generate_random_fp2(a);
    generate_random_fp2(b);
    auto p = get_BLST_P();

    for (auto _ : state) {
        add_mod_384x(&c, &a, &b, &p);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Addition_asm);
#endif

// Benchmark for Fp2 addition (C version)
static void BM_Fp2Addition_cc(benchmark::State& state) {
    Fp2 a, b, c;
    generate_random_fp2(a);
    generate_random_fp2(b);
    auto p = get_BLST_P();

    for (auto _ : state) {
        cc_add_mod_384x(c, a, b, p);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Addition_cc);

// Benchmark for Fp2 addition (High-level API)
static void BM_Fp2Addition_api(benchmark::State& state) {
    Fp2 a, b, c;
    generate_random_fp2(a);
    generate_random_fp2(b);

    for (auto _ : state) {
        Fp2::add(c, a, b);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Addition_api);

// ========== SUBTRACTION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp2 subtraction (ASM version)
static void BM_Fp2Subtraction_asm(benchmark::State& state) {
    Fp2 a, b, c;
    generate_random_fp2(a);
    generate_random_fp2(b);
    auto p = get_BLST_P();

    for (auto _ : state) {
        sub_mod_384x(&c, &a, &b, &p);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Subtraction_asm);
#endif

// Benchmark for Fp2 subtraction (C version)
static void BM_Fp2Subtraction_cc(benchmark::State& state) {
    Fp2 a, b, c;
    generate_random_fp2(a);
    generate_random_fp2(b);
    auto p = get_BLST_P();

    for (auto _ : state) {
        cc_sub_mod_384x(c, a, b, p);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Subtraction_cc);

// Benchmark for Fp2 subtraction (High-level API)
static void BM_Fp2Subtraction_api(benchmark::State& state) {
    Fp2 a, b, c;
    generate_random_fp2(a);
    generate_random_fp2(b);

    for (auto _ : state) {
        Fp2::sub(c, a, b);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Subtraction_api);

// ========== MULTIPLICATION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp2 multiplication (ASM version)
static void BM_Fp2Multiplication_asm(benchmark::State& state) {
    Fp2 a, b, c;
    generate_random_fp2(a);
    generate_random_fp2(b);
    auto p = get_BLST_P();

    for (auto _ : state) {
        mulx_mont_384x(&c, &a, &b, &p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Multiplication_asm);
#endif

// Benchmark for Fp2 multiplication (C version)
static void BM_Fp2Multiplication_cc(benchmark::State& state) {
    Fp2 a, b, c;
    generate_random_fp2(a);
    generate_random_fp2(b);
    auto p = get_BLST_P();

    for (auto _ : state) {
        cc_mul_mont_384x(c, a, b, p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Multiplication_cc);

// Benchmark for Fp2 multiplication (High-level API)
static void BM_Fp2Multiplication_api(benchmark::State& state) {
    Fp2 a, b, c;
    generate_random_fp2(a);
    generate_random_fp2(b);

    for (auto _ : state) {
        Fp2::mul(c, a, b);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Multiplication_api);

// ========== SQUARING BENCHMARKS ==========


#ifdef BLST_ASM
// Benchmark for Fp2 squaring (ASM version)
static void BM_Fp2Squaring_asm(benchmark::State& state) {
    Fp2 a, c;
    generate_random_fp2(a);
    auto p = get_BLST_P();

    for (auto _ : state) {
        sqrx_mont_384x(&c, &a, &p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Squaring_asm);
#endif

// Benchmark for Fp2 squaring (C version)
static void BM_Fp2Squaring_cc(benchmark::State& state) {
    Fp2 a, c;
    generate_random_fp2(a);
    auto p = get_BLST_P();

    for (auto _ : state) {
        cc_sqr_mont_384x(c, a, p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Squaring_cc);

// Benchmark for Fp2 squaring (High-level API)
static void BM_Fp2Squaring_api(benchmark::State& state) {
    Fp2 a, c;
    generate_random_fp2(a);

    for (auto _ : state) {
        Fp2::sqr(c, a);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Squaring_api);

// ========== INVERSION BENCHMARKS ==========

// Benchmark for Fp2 inversion (fast reciprocal)
static void BM_Fp2Inversion_fast(benchmark::State& state) {
    Fp2 a, result;
    generate_random_fp2(a);

    for (auto _ : state) {
        Fp2::flt_reciprocal(result, a);
        benchmark::DoNotOptimize(result);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Inversion_fast);

// Benchmark for Fp2 inversion (constant-time reciprocal)
static void BM_Fp2Inversion_ct(benchmark::State& state) {
    Fp2 a, result;
    generate_random_fp2(a);

    for (auto _ : state) {
        Fp2::reciprocal(result, a);
        benchmark::DoNotOptimize(result);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Inversion_ct);

// ========== SQUARE ROOT BENCHMARKS ==========

// Benchmark for Fp2 square root
static void BM_Fp2Sqrt(benchmark::State& state) {
    Fp2 a, b;
    generate_random_fp2(a);
    Fp2::sqr(a, a); // Ensure it has a square root

    for (auto _ : state) {
        Fp2::sqrt(b, a);
        benchmark::DoNotOptimize(b);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Sqrt);

// ========== CONJUGATION BENCHMARKS ==========

// Benchmark for Fp2 conjugation
static void BM_Fp2Conjugate(benchmark::State& state) {
    Fp2 a, c;
    generate_random_fp2(a);

    for (auto _ : state) {
        Fp2::conjugate(c, a);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp2Conjugate);

