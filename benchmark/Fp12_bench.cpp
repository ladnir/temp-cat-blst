#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>
#include <thrust/device_vector.h>

// Include the actual BLST field headers
#include "blst/fields/Fp12.h"

using namespace blst;

// Helper function to generate random Fp12 elements
void generate_random_fp12(Fp12& fp12) {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint64_t> dis;

    limb_t* p = (limb_t*)&fp12;
    for (size_t i = 0; i < sizeof(Fp12) / sizeof(limb_t); ++i) {
        p[i] = dis(gen);
    }

	// Ensure the generated element is reduced modulo the field characteristic
	for (auto& comp : fp12.data) {
        for (auto& part : comp.data) {
            for (auto& fp : part.data) {
                while (Fp::is_redc(fp) == 0) {
                    Fp::rshift(fp, fp, 1);
                }
            }
        }
    }
}

// ========== MULTIPLICATION BENCHMARKS ==========

// Benchmark for Fp12 multiplication
static void BM_Fp12Multiplication(benchmark::State& state) {
    Fp12 a, b, c;
    generate_random_fp12(a);
    generate_random_fp12(b);

    for (auto _ : state) {
        Fp12::mul(c, a, b);
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
        Fp12::sqr(c, a);
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
        Fp12::cyclotomic_sqr(c, a);
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
        Fp12::inverse(c, a);
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
        Fp12::conjugate(a);
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
        Fp12::frobenius_map(c, a, 1);
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
    byte exponent[32];

	std::memset(exponent, 128 | 0xcc, 32); // Example exponent

    for (auto _ : state) {
        Fp12::exp(c, a, exponent, 256);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12Exponentiation);


// Helper: Generate a random scalar (256-bit)
void generate_random_scalar(byte scalar[33]) {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint64_t> dis;

    uint64_t* p = (uint64_t*)scalar;
    for (int i = 0; i < 4; ++i) {
        p[i] = dis(gen);
    }
}

// ... existing benchmarks ...

static void BM_Fp12CyclotomicExponentiation(benchmark::State& state) {
    Fp12 a, base, c;
    generate_random_fp12(a);
    Fp12::final_exp(base, a); // Ensure base is in cyclotomic subgroup

    byte scalar[33] = { 0 };
    generate_random_scalar(scalar);

    for (auto _ : state) {
        Fp12::cyclotomic_exp(c, base, scalar, 257);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12CyclotomicExponentiation);

static void BM_Fp12ExponentiationGLV(benchmark::State& state) {
    Fp12 a, base, c;
    generate_random_fp12(a);
    Fp12::final_exp(base, a); // Ensure base is in cyclotomic subgroup

    byte scalar[32];
    generate_random_scalar(scalar);

    for (auto _ : state) {
        Fp12::cyclotomic_exp_glv(c, base, scalar);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12ExponentiationGLV);

// We need a helper to split scalars for the pre-split benchmark.
// Since the library function might be internal, we duplicate the call logic
// or use the Fr helper if available. 
// Assuming Fr.h is available and has the constants/helpers we need.
// If gt_split_scalar_phi2_u256_barrett is not public, we can't call it directly.
// HOWEVER, we can just benchmark exp_glv with "random" 128-bit k1/k2 
// which is effectively what the split produces performance-wise.

static void BM_Fp12ExponentiationGLV_PreSplit(benchmark::State& state) {
    Fp12 a, base, c;
    generate_random_fp12(a);
    Fp12::final_exp(base, a);

    // Generate random 128-bit splits (simulation of pre-computation)
    // k1, k2 are approx 128 bits for 256-bit scalars.
    std::vector<std::pair<std::array<limb_t, 4>, std::array<limb_t, 4>>> splits(128);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;

    for (auto& s : splits) {
        // Fill 128 bits (2 limbs) randomly, zero upper limbs
        s.first[0] = dis(gen); s.first[1] = dis(gen); s.first[2] = 0; s.first[3] = 0;
        s.second[0] = dis(gen); s.second[1] = dis(gen); s.second[2] = 0; s.second[3] = 0;
    }

    size_t idx = 0;
    const auto& s = splits[idx & 127];
    for (auto _ : state) {
        Fp12::cyclotomic_exp_glv(c, base, s.first.data(), s.second.data());
        benchmark::DoNotOptimize(c);
        idx++;
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12ExponentiationGLV_PreSplit);


//// We need a helper to split scalars for the pre-split benchmark.
//// Since the library function might be internal, we duplicate the call logic
//// or use the Fr helper if available. 
//// Assuming Fr.h is available and has the constants/helpers we need.
//// If gt_split_scalar_phi2_u256_barrett is not public, we can't call it directly.
//// HOWEVER, we can just benchmark exp_glv with "random" 128-bit k1/k2 
//// which is effectively what the split produces performance-wise.
//
//static void BM_Fp12ExponentiationGLV_opt(benchmark::State& state) {
//    Fp12 a, base, c;
//    generate_random_fp12(a);
//    Fp12::final_exp(base, a);
//
//    // Generate random 128-bit splits (simulation of pre-computation)
//    // k1, k2 are approx 128 bits for 256-bit scalars.
//    std::vector<std::pair<std::array<limb_t, 4>, std::array<limb_t, 4>>> splits(128);
//
//    std::random_device rd;
//    std::mt19937_64 gen(rd());
//    std::uniform_int_distribution<uint64_t> dis;
//
//    for (auto& s : splits) {
//        // Fill 128 bits (2 limbs) randomly, zero upper limbs
//        s.first[0] = dis(gen); s.first[1] = dis(gen); s.first[2] = 0; s.first[3] = 0;
//        s.second[0] = dis(gen); s.second[1] = dis(gen); s.second[2] = 0; s.second[3] = 0;
//    }
//
//    size_t idx = 0;
//    const auto& s = splits[idx & 127];
//    for (auto _ : state) {
//        Fp12::cyclotomic_exp_glv_optimized(c, base, s.first.data(), s.second.data());
//        benchmark::DoNotOptimize(c);
//        idx++;
//    }
//    state.SetItemsProcessed(state.iterations());
//}
//BENCHMARK(BM_Fp12ExponentiationGLV_opt);
//namespace blst {
//    BLST_HD int recode_glv_jsf(int8_t jsf[][2], const limb_t* k1, const limb_t* k2);
//
//    BLST_HD void cyclotomic_exp_glv_opt(Fp12& out, const Fp12& base, int8_t jsf[132][2], int len);
//}

// We need a helper to split scalars for the pre-split benchmark.
// Since the library function might be internal, we duplicate the call logic
// or use the Fr helper if available. 
// Assuming Fr.h is available and has the constants/helpers we need.
// If gt_split_scalar_phi2_u256_barrett is not public, we can't call it directly.
// HOWEVER, we can just benchmark exp_glv with "random" 128-bit k1/k2 
// which is effectively what the split produces performance-wise.

static void BM_Fp12ExponentiationGLV_opt(benchmark::State& state) {
    Fp12 a, base, c;
    generate_random_fp12(a);
    Fp12::final_exp(base, a);

    // Generate random 128-bit splits (simulation of pre-computation)
    // k1, k2 are approx 128 bits for 256-bit scalars.
    //std::vector<std::pair<std::array<limb_t, 4>, std::array<limb_t, 4>>> splits(128);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;
    std::pair<std::array<limb_t, 4>, std::array<limb_t, 4>> s;

    //for (auto& s : splits) 
    {
        // Fill 128 bits (2 limbs) randomly, zero upper limbs
        s.first[0] = dis(gen); s.first[1] = dis(gen); s.first[2] = 0; s.first[3] = 0;
        s.second[0] = dis(gen); s.second[1] = dis(gen); s.second[2] = 0; s.second[3] = 0;
    }

    int8_t jsf[132][2];
	auto len = Fp12::glv_recode_jsf(jsf, s.first.data(), s.second.data());

    size_t idx = 0;
    //const auto& s = splits[idx & 127];
    for (auto _ : state) {
        //Fp12::cyclotomic_exp_glv_optimized(c, base, s.first.data(), s.second.data());
		Fp12::cyclotomic_exp_glv_jsf(c, base, jsf, len);
        benchmark::DoNotOptimize(c);
        idx++;
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12ExponentiationGLV_opt);



// 1. Basic GLV-4 Benchmark (Runtime Split + Exp)
static void BM_Fp12ExponentiationGLV4(benchmark::State& state) {
    Fp12 a, base, c;
    generate_random_fp12(a);
    Fp12::final_exp(base, a); // Ensure base is in cyclotomic subgroup

    byte scalar[32];
    generate_random_scalar(scalar);

    // Convert scalar bytes to limbs for the API
    limb_t k[4];
    std::memcpy(k, scalar, 32);

    for (auto _ : state) {
        // Calls the full GLV-4 logic including decomposition
        Fp12::cyclotomic_exp_glv4(c, base, k);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12ExponentiationGLV4);

// 2. Constant Scalar GLV-4 Benchmark (Offline Split, Runtime Exp)
static void BM_Fp12ExponentiationGLV4_Constant(benchmark::State& state) {
    Fp12 a, base, c;
    generate_random_fp12(a);
    Fp12::final_exp(base, a);

    byte scalar[32];
    generate_random_scalar(scalar);
    limb_t k[4];
    std::memcpy(k, scalar, 32);

    // --- Offline Compilation Phase (Not Timed) ---
    // Compile the scalar into properties (Schedule + Sign Mask)
    GLV4_Exp props = Fp12::glv4_compile(k);

    // --- Runtime Benchmark Phase ---
    for (auto _ : state) {
        // Only run the execution phase, simulating fixed scalar usage
        Fp12::cyclotomic_exp_glv4(c, base, props);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Fp12ExponentiationGLV4_Constant);