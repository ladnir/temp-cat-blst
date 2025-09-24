#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>

// Include the actual BLST field headers
#include "blst/fields/Fp.h"
#include "blst/backend/backend.h"

using namespace blst;

#define LIMBS NLIMBS(384)

// Helper function to set vec384 from uint64_t value
void set_vec384(vec384 v, uint64_t x) {
    for (int i = 0; i < LIMBS; ++i)
        v[i] = 0;
    memcpy(v, &x, sizeof(x));
}

// Helper function to generate random Fp elements
void generate_random_fp(vec384 fp) {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint64_t> dis;
    
    // Initialize with random values
    for (int i = 0; i < LIMBS; ++i) {
        fp[i] = dis(gen);
    }
}

// ========== ADDITION BENCHMARKS ==========

// Benchmark for Fp addition (ASM version)
static void BM_FpAddition_asm(benchmark::State& state) {
    vec384 a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        add_mod_384(c, a, b, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpAddition_asm);

// Benchmark for Fp addition (C version)
static void BM_FpAddition_cc(benchmark::State& state) {
    vec384 a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_add_mod_384(c, a, b, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpAddition_cc);

// Benchmark for Fp addition (High-level API)
static void BM_FpAddition_api(benchmark::State& state) {
    vec384 a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    
    for (auto _ : state) {
        add_fp(c, a, b);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpAddition_api);

// ========== SUBTRACTION BENCHMARKS ==========

// Benchmark for Fp subtraction (ASM version)
static void BM_FpSubtraction_asm(benchmark::State& state) {
    vec384 a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        sub_mod_384(c, a, b, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSubtraction_asm);

// Benchmark for Fp subtraction (C version)
static void BM_FpSubtraction_cc(benchmark::State& state) {
    vec384 a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_sub_mod_n<NLIMBS(384)>(c, a, b, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSubtraction_cc);

// Benchmark for Fp subtraction (High-level API)
static void BM_FpSubtraction_api(benchmark::State& state) {
    vec384 a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    
    for (auto _ : state) {
        sub_fp(c, a, b);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSubtraction_api);

// ========== MULTIPLICATION BENCHMARKS ==========

// Benchmark for Fp multiplication (ASM version)
static void BM_FpMultiplication_asm(benchmark::State& state) {
    vec384 a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    
    // Convert to Montgomery form for multiplication
    to_fp(a, a);
    to_fp(b, b);
    
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        mulx_mont_384(c, a, b, *p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMultiplication_asm);

// Benchmark for Fp multiplication (C version)
static void BM_FpMultiplication_cc(benchmark::State& state) {
    vec384 a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    
    // Convert to Montgomery form for multiplication
    to_fp(a, a);
    to_fp(b, b);
    
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_mul_mont_384(c, a, b, *p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMultiplication_cc);

// Benchmark for Fp multiplication (High-level API)
static void BM_FpMultiplication_api(benchmark::State& state) {
    vec384 a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    
    // Convert to Montgomery form for multiplication
    to_fp(a, a);
    to_fp(b, b);
    
    for (auto _ : state) {
        mul_fp(c, a, b);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMultiplication_api);

// ========== SQUARING BENCHMARKS ==========

// Benchmark for Fp squaring (ASM version)
static void BM_FpSquaring_asm(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    
    // Convert to Montgomery form for squaring
    to_fp(a, a);
    
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        sqrx_mont_384(c, a, *p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSquaring_asm);

// Benchmark for Fp squaring (C version)
static void BM_FpSquaring_cc(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    
    // Convert to Montgomery form for squaring
    to_fp(a, a);
    
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_sqr_mont_384(c, a, *p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSquaring_cc);

// Benchmark for Fp squaring (High-level API)
static void BM_FpSquaring_api(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    
    // Convert to Montgomery form for squaring
    to_fp(a, a);
    
    for (auto _ : state) {
        sqr_fp(c, a);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSquaring_api);

// ========== MULTIPLY BY 3 BENCHMARKS ==========

// Benchmark for Fp multiplication by 3 (ASM version)
static void BM_FpMulBy3_asm(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        mul_by_3_mod_384(c, a, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy3_asm);

// Benchmark for Fp multiplication by 3 (C version)
static void BM_FpMulBy3_cc(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_mul_by_3_mod_384(c, a, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy3_cc);

// Benchmark for Fp multiplication by 3 (High-level API)
static void BM_FpMulBy3_api(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    
    for (auto _ : state) {
        mul_by_3_fp(c, a);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy3_api);

// ========== MULTIPLY BY 8 BENCHMARKS ==========

// Benchmark for Fp multiplication by 8 (ASM version)
static void BM_FpMulBy8_asm(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        mul_by_8_mod_384(c, a, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy8_asm);

// Benchmark for Fp multiplication by 8 (C version)
static void BM_FpMulBy8_cc(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_mul_by_8_mod_384(c, a, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy8_cc);

// Benchmark for Fp multiplication by 8 (High-level API)
static void BM_FpMulBy8_api(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    
    for (auto _ : state) {
        mul_by_8_fp(c, a);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy8_api);

// ========== RIGHT SHIFT BENCHMARKS ==========

// Benchmark for Fp right shift (ASM version)
static void BM_FpRshift_asm(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        rshift_mod_384(c, a, 2, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpRshift_asm);

// Benchmark for Fp right shift (C version)
static void BM_FpRshift_cc(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_rshift_mod_384(c, a, 2, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpRshift_cc);

// Benchmark for Fp right shift (High-level API)
static void BM_FpRshift_api(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    
    for (auto _ : state) {
        rshift_fp(c, a, 2);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpRshift_api);

// ========== DIVISION BY 2 BENCHMARKS ==========

// Benchmark for Fp division by 2 (ASM version)
static void BM_FpDivBy2_asm(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        div_by_2_mod_384(c, a, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpDivBy2_asm);

// Benchmark for Fp division by 2 (C version)
static void BM_FpDivBy2_cc(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_div_by_2_mod_384(c, a, *p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpDivBy2_cc);

// Benchmark for Fp division by 2 (High-level API)
static void BM_FpDivBy2_api(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    
    for (auto _ : state) {
        div_by_2_fp(c, a);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpDivBy2_api);

// ========== MONTGOMERY CONVERSION BENCHMARKS ==========

// Benchmark for Montgomery conversion from (ASM version)
static void BM_FpFromMontgomery_asm(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    to_fp(a, a);  // Convert to Montgomery form first
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        fromx_mont_384(c, a, *p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpFromMontgomery_asm);

// Benchmark for Montgomery conversion from (C version)
static void BM_FpFromMontgomery_cc(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    to_fp(a, a);  // Convert to Montgomery form first
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_from_mont_384(c, a, *p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpFromMontgomery_cc);

// Benchmark for Montgomery conversion from (High-level API)
static void BM_FpFromMontgomery_api(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    to_fp(a, a);  // Convert to Montgomery form first
    
    for (auto _ : state) {
        from_fp(c, a);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpFromMontgomery_api);

// ========== CONDITIONAL NEGATION BENCHMARKS ==========

// Benchmark for Fp conditional negation (ASM version)
static void BM_FpCneg_asm(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    bool flag = true;
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cneg_mod_384(c, a, flag, *p);
        benchmark::DoNotOptimize(c);
        flag = !flag;
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpCneg_asm);

// Benchmark for Fp conditional negation (C version)
static void BM_FpCneg_cc(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    bool flag = true;
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_cneg_mod_384(c, a, flag, *p);
        benchmark::DoNotOptimize(c);
        flag = !flag;
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpCneg_cc);

// Benchmark for Fp conditional negation (High-level API)
static void BM_FpCneg_api(benchmark::State& state) {
    vec384 a, c;
    generate_random_fp(a);
    bool flag = true;
    
    for (auto _ : state) {
        cneg_fp(c, a, flag);
        benchmark::DoNotOptimize(c);
        flag = !flag;
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpCneg_api);

// ========== EXPONENTIATION BENCHMARKS ==========

// Benchmark for exponentiation using sqr_n_mul (ASM version)
static void BM_FpExponentiation_asm(benchmark::State& state) {
    vec384 a, result, one;
    set_vec384(a, 3);  // Base value
    to_fp(a, a);       // Convert to Montgomery form
    set_vec384(one, 1);
    to_fp(one, one);   // One in Montgomery form
    
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        // Use sqr_n_mul_mont_383 for efficient exponentiation
        // This computes a^(2^count) * b
        sqrx_n_mul_mont_383(result, a, 3, *p, BLST_p0, one); // a^8 * 1 = a^8
        benchmark::DoNotOptimize(result);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpExponentiation_asm);

// Benchmark for exponentiation using sqr_n_mul (C version)
static void BM_FpExponentiation_cc(benchmark::State& state) {
    vec384 a, result, one;
    set_vec384(a, 3);  // Base value
    to_fp(a, a);       // Convert to Montgomery form
    set_vec384(one, 1);
    to_fp(one, one);   // One in Montgomery form
    
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        // Use C version of sqr_n_mul_mont
        cc_sqr_n_mul_mont_384(result, a, 3, *p, BLST_p0, one); // a^8 * 1 = a^8
        benchmark::DoNotOptimize(result);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpExponentiation_cc);

// Benchmark for exponentiation using sqr_n_mul (High-level API)
static void BM_FpExponentiation_api(benchmark::State& state) {
    vec384 a, result, one;
    set_vec384(a, 3);  // Base value
    to_fp(a, a);       // Convert to Montgomery form
    set_vec384(one, 1);
    to_fp(one, one);   // One in Montgomery form
    
    for (auto _ : state) {
        // Use high-level API
        sqr_n_mul_fp(result, a, 3, one); // a^8 * 1 = a^8
        benchmark::DoNotOptimize(result);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpExponentiation_api);

// ========== INVERSION BENCHMARKS ==========

// Benchmark for Fp inversion (fast reciprocal)
static void BM_FpInversion_fast(benchmark::State& state) {
    vec384 a, result;
    generate_random_fp(a);
    to_fp(a, a);  // Convert to Montgomery form
    
    for (auto _ : state) {
        flt_reciprocal_fp(result, a);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpInversion_fast);

// Benchmark for Fp inversion (constant-time reciprocal)
static void BM_FpInversion_ct(benchmark::State& state) {
    vec384 a, result;
    generate_random_fp(a);
    to_fp(a, a);  // Convert to Montgomery form
    
    for (auto _ : state) {
        reciprocal_fp(result, a);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpInversion_ct);

// ========== COMPLEX OPERATIONS BENCHMARKS ==========

// Benchmark for complex operation combining multiple operations
static void BM_FpComplexOperation(benchmark::State& state) {
    vec384 a, b, c, d, temp1, temp2;
    generate_random_fp(a);
    generate_random_fp(b);
    
    // Convert to Montgomery form
    to_fp(a, a);
    to_fp(b, b);
    
    for (auto _ : state) {
        // Perform: d = (a * b) + (a^2) - b
        mul_fp(temp1, a, b);      // temp1 = a * b
        sqr_fp(c, a);             // c = a^2  
        add_fp(temp2, temp1, c);  // temp2 = temp1 + c
        sub_fp(d, temp2, b);      // d = temp2 - b
        benchmark::DoNotOptimize(d);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpComplexOperation);

// Benchmark for Montgomery round-trip (to_fp + from_fp)
static void BM_FpMontgomeryRoundTrip(benchmark::State& state) {
    vec384 a, mont, result;
    generate_random_fp(a);
    
    for (auto _ : state) {
        to_fp(mont, a);         // Convert to Montgomery form
        from_fp(result, mont);  // Convert back from Montgomery form
        benchmark::DoNotOptimize(result);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMontgomeryRoundTrip);

