#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>

// Include the actual BLST field headers
#include "blst/fields/Fp.h"
#include "blst/backend/backend.h"

using namespace blst;

#define LIMBS NLIMBS(384)

// Helper function to set Fp from uint64_t value
void set_Fp(Fp& v, uint64_t x) {
    for (int i = 0; i < LIMBS; ++i)
        v[i] = 0;
    memcpy(v, &x, sizeof(x));
}

// Helper function to generate random Fp elements
void generate_random_fp(Fp& fp) {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint64_t> dis;
    
    // Initialize with random values
    for (int i = 0; i < LIMBS; ++i) {
        fp[i] = dis(gen);
    }
}

// ========== ADDITION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp addition (ASM version)
static void BM_FpAddition_asm(benchmark::State& state) {
    Fp a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        add_mod_384(&c, &a, &b, &p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpAddition_asm);
#endif

// Benchmark for Fp addition (C version)
static void BM_FpAddition_cc(benchmark::State& state) {
    Fp a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        //Fp::add(c, a, b);
        cc_add_mod_n<NLIMBS(384)>(c.data, a.data, b.data, get_BLST_P().data);

        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpAddition_cc);

// Benchmark for Fp addition (High-level API)
static void BM_FpAddition_api(benchmark::State& state) {
    Fp a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    
    for (auto _ : state) {
        Fp::add(c, a, b);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpAddition_api);

// ========== SUBTRACTION BENCHMARKS ==========

#ifdef BLST_ASM

// Benchmark for Fp subtraction (ASM version)
static void BM_FpSubtraction_asm(benchmark::State& state) {
    Fp a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        sub_mod_384(&c, &a, &b, &p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSubtraction_asm);

#endif

// Benchmark for Fp subtraction (C version)
static void BM_FpSubtraction_cc(benchmark::State& state) {
    Fp a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_sub_mod_n<NLIMBS(384)>(c, a, b, p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSubtraction_cc);

// Benchmark for Fp subtraction (High-level API)
static void BM_FpSubtraction_api(benchmark::State& state) {
    Fp a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    
    for (auto _ : state) {
        Fp::sub(c, a, b);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSubtraction_api);

// ========== MULTIPLICATION BENCHMARKS ==========
#ifdef BLST_ASM

// Benchmark for Fp multiplication (ASM version)
static void BM_FpMultiplication_asm(benchmark::State& state) {
    Fp a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    
    // Convert to Montgomery form for multiplication
    Fp::to_mont(a, a);
    Fp::to_mont(b, b);
    
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        mulx_mont_384(&c, &a, &b, &p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMultiplication_asm);
#endif

// Benchmark for Fp multiplication (C version)
static void BM_FpMultiplication_cc(benchmark::State& state) {
    Fp a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    
    // Convert to Montgomery form for multiplication
    Fp::to_mont(a, a);
    Fp::to_mont(b, b);
    
    //auto p = get_BLST_P();
    
    for (auto _ : state) {
        //Fp::mul(c, a, b);
        cc_mul_mont_n<NLIMBS(384)>(c.data, a.data, b.data, get_BLST_P().data, BLST_p0);

        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMultiplication_cc);

// Benchmark for Fp multiplication (High-level API)
static void BM_FpMultiplication_api(benchmark::State& state) {
    Fp a, b, c;
    generate_random_fp(a);
    generate_random_fp(b);
    
    // Convert to Montgomery form for multiplication
    Fp::to_mont(a, a);
    Fp::to_mont(b, b);
    
    for (auto _ : state) {
        Fp::mul(c, a, b);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMultiplication_api);

// ========== SQUARING BENCHMARKS ==========
#ifdef BLST_ASM

// Benchmark for Fp squaring (ASM version)
static void BM_FpSquaring_asm(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    
    // Convert to Montgomery form for squaring
    Fp::to_mont(a, a);
    
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        sqrx_mont_384(&c, &a, &p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSquaring_asm);
#endif

// Benchmark for Fp squaring (C version)
static void BM_FpSquaring_cc(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    
    // Convert to Montgomery form for squaring
    Fp::to_mont(a, a);
    
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        //Fp::sqr(c, a);
        cc_mul_mont_n<NLIMBS(384)>(c.data, a.data, a.data, get_BLST_P().data, BLST_p0);

        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSquaring_cc);

// Benchmark for Fp squaring (High-level API)
static void BM_FpSquaring_api(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    
    // Convert to Montgomery form for squaring
    Fp::to_mont(a, a);
    
    for (auto _ : state) {
        Fp::sqr(c, a);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpSquaring_api);

// ========== MULTIPLY BY 3 BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp multiplication by 3 (ASM version)
static void BM_FpMulBy3_asm(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        mul_by_3_mod_384(&c, &a, &p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy3_asm);
#endif

// Benchmark for Fp multiplication by 3 (C version)
static void BM_FpMulBy3_cc(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        //Fp::mul_by_3(c, a);
        cc_mul_by_3_mod_n<NLIMBS(384)>(c.data, a.data, get_BLST_P().data);

        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy3_cc);

// Benchmark for Fp multiplication by 3 (High-level API)
static void BM_FpMulBy3_api(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    
    for (auto _ : state) {
        Fp::mul_by_3(c, a);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy3_api);

// ========== MULTIPLY BY 8 BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp multiplication by 8 (ASM version)
static void BM_FpMulBy8_asm(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        mul_by_8_mod_384(&c,&a, &p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy8_asm);
#endif

// Benchmark for Fp multiplication by 8 (C version)
static void BM_FpMulBy8_cc(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cc_mul_by_8_mod_384(c, a, p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy8_cc);

// Benchmark for Fp multiplication by 8 (High-level API)
static void BM_FpMulBy8_api(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    
    for (auto _ : state) {
        Fp::mul_by_8(c, a);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMulBy8_api);

// ========== RIGHT SHIFT BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp right shift (ASM version)
static void BM_FpRshift_asm(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        rshift_mod_384(&c, &a, 2, &p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpRshift_asm);
#endif

// Benchmark for Fp right shift (C version)
static void BM_FpRshift_cc(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        //Fp::rshift(c, a, 2);
        cc_rshift_mod_n<NLIMBS(384)>(c.data, a.data, 2, get_BLST_P().data);

        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpRshift_cc);

// Benchmark for Fp right shift (High-level API)
static void BM_FpRshift_api(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    
    for (auto _ : state) {
        Fp::rshift(c, a, 2);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpRshift_api);



// ========== LEFT  SHIFT BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp right shift (ASM version)
static void BM_FpLshift_asm(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();

    for (auto _ : state) {
        lshift_mod_384(&c, &a, 2, &p);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpLshift_asm);
#endif

// Benchmark for Fp left  shift (C version)
static void BM_FpLshift_cc(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();

    for (auto _ : state) {
        //Fp::rshift(c, a, 2);
        cc_lshift_mod_n<NLIMBS(384)>(c.data, a.data, 2, get_BLST_P().data);

        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpLshift_cc);

// Benchmark for Fp left  shift (High-level API)
static void BM_FpLshift_api(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);

    for (auto _ : state) {
        Fp::lshift(c, a, 2);
        benchmark::DoNotOptimize(c);
    }

    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpLshift_api);


// ========== DIVISION BY 2 BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp division by 2 (ASM version)
static void BM_FpDivBy2_asm(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        div_by_2_mod_384(&c, &a, &p);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpDivBy2_asm);
#endif

// Benchmark for Fp division by 2 (C version)
static void BM_FpDivBy2_cc(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        //cc_div_by_2_mod_384(c, a, p);
        cc_rshift_mod_n<NLIMBS(384)>(c.data, a.data, 1, get_BLST_P().data);

        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpDivBy2_cc);

// Benchmark for Fp division by 2 (High-level API)
static void BM_FpDivBy2_api(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    
    for (auto _ : state) {
        Fp::div_by_2(c, a);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpDivBy2_api);

// ========== MONTGOMERY CONVERSION BENCHMARKS ==========


#ifdef BLST_ASM
// Benchmark for Montgomery conversion from (ASM version)
static void BM_FpFromMontgomery_asm(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    Fp::to_mont(a, a);  // Convert to Montgomery form first
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        fromx_mont_384(&c, &a, &p, BLST_p0);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpFromMontgomery_asm);
#endif

// Benchmark for Montgomery conversion from (C version)
static void BM_FpFromMontgomery_cc(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    Fp::to_mont(a, a);  // Convert to Montgomery form first
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        //cc_from_mont_384(c, a, p, BLST_p0);
        cc_from_mont_n<NLIMBS(384)>(c.data, a.data, get_BLST_P().data, BLST_p0);

        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpFromMontgomery_cc);

// Benchmark for Montgomery conversion from (High-level API)
static void BM_FpFromMontgomery_api(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    Fp::to_mont(a, a);  // Convert to Montgomery form first
    
    for (auto _ : state) {
        Fp::from_mont(c, a);
        benchmark::DoNotOptimize(c);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpFromMontgomery_api);

// ========== CONDITIONAL NEGATION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for Fp conditional negation (ASM version)
static void BM_FpCneg_asm(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    bool flag = true;
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        cneg_mod_384(&c, &a, flag, &p);
        benchmark::DoNotOptimize(c);
        flag = !flag;
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpCneg_asm);
#endif

// Benchmark for Fp conditional negation (C version)
static void BM_FpCneg_cc(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    bool flag = true;
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        //cc_cneg_mod_384(c, a, flag, p);
        cc_cneg_mod_n<NLIMBS(384)>(c.data, a.data, flag, p.data);

        benchmark::DoNotOptimize(c);
        flag = !flag;
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpCneg_cc);

// Benchmark for Fp conditional negation (High-level API)
static void BM_FpCneg_api(benchmark::State& state) {
    Fp a, c;
    generate_random_fp(a);
    bool flag = true;
    
    for (auto _ : state) {
        Fp::cneg(c, a, flag);
        benchmark::DoNotOptimize(c);
        flag = !flag;
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpCneg_api);

// ========== EXPONENTIATION BENCHMARKS ==========

#ifdef BLST_ASM
// Benchmark for exponentiation using sqr_n_mul (ASM version)
static void BM_FpExponentiation_asm(benchmark::State& state) {
    Fp a, result, one;
    set_Fp(a, 3);  // Base value
    Fp::to_mont(a, a);       // Convert to Montgomery form
    set_Fp(one, 1);
    Fp::to_mont(one, one);   // One in Montgomery form
    
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        // Use sqr_n_mul_mont_383 for efficient exponentiation
        // This computes a^(2^count) * b
        sqrx_n_mul_mont_383(&result, &a, 3, &p, BLST_p0, &one); // a^8 * 1 = a^8
        benchmark::DoNotOptimize(result);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpExponentiation_asm);
#endif

// Benchmark for exponentiation using sqr_n_mul (C version)
static void BM_FpExponentiation_cc(benchmark::State& state) {
    Fp a, result, one;
    set_Fp(a, 3);  // Base value
    Fp::to_mont(a, a);       // Convert to Montgomery form
    set_Fp(one, 1);
    Fp::to_mont(one, one);   // One in Montgomery form
    
    auto p = get_BLST_P();
    
    for (auto _ : state) {
        // Use C version of sqr_n_mul_mont
        cc_sqr_n_mul_mont_384(result, a, 3, p, BLST_p0, one); // a^8 * 1 = a^8
        benchmark::DoNotOptimize(result);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpExponentiation_cc);

// Benchmark for exponentiation using sqr_n_mul (High-level API)
static void BM_FpExponentiation_api(benchmark::State& state) {
    Fp a, result, one;
    set_Fp(a, 3);  // Base value
    Fp::to_mont(a, a);       // Convert to Montgomery form
    set_Fp(one, 1);
    Fp::to_mont(one, one);   // One in Montgomery form
    
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
    Fp a, result;
    generate_random_fp(a);
    Fp::to_mont(a, a);  // Convert to Montgomery form
    
    for (auto _ : state) {
        flt_reciprocal_fp(result, a);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpInversion_fast);

// Benchmark for Fp inversion (constant-time reciprocal)
static void BM_FpInversion_ct(benchmark::State& state) {
    Fp a, result;
    generate_random_fp(a);
    Fp::to_mont(a, a);  // Convert to Montgomery form
    
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
    Fp a, b, c, d, temp1, temp2;
    generate_random_fp(a);
    generate_random_fp(b);
    
    // Convert to Montgomery form
    Fp::to_mont(a, a);
    Fp::to_mont(b, b);
    
    for (auto _ : state) {
        // Perform: d = (a * b) + (a^2) - b
        Fp::mul(temp1, a, b);      // temp1 = a * b
        Fp::sqr(c, a);             // c = a^2  
        Fp::add(temp2, temp1, c);  // temp2 = temp1 + c
        Fp::sub(d, temp2, b);      // d = temp2 - b
        benchmark::DoNotOptimize(d);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpComplexOperation);

// Benchmark for Montgomery round-trip (Fp::to_mont + from_fp)
static void BM_FpMontgomeryRoundTrip(benchmark::State& state) {
    Fp a, mont, result;
    generate_random_fp(a);
    
    for (auto _ : state) {
        Fp::to_mont(mont, a);         // Convert to Montgomery form
        Fp::from_mont(result, mont);  // Convert back from Montgomery form
        benchmark::DoNotOptimize(result);
    }
    
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_FpMontgomeryRoundTrip);

