#include <benchmark/benchmark.h>

#include "blst/vroom/src/final_exponentiation.hpp"
#include "blst/vroom/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::Vroom50Fp12;
using blst_vroom_test::Vroom50Fp12Ring;
using blst_vroom_test::Vroom50Ring;

static void BM_FP12_Mul_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    volatile uint64_t left_seed = 42;
    volatile uint64_t right_seed = 123;
    const auto left = blst_vroom_test::make_vroom_fp12_from_seed(ring, left_seed);
    const auto right = blst_vroom_test::make_vroom_fp12_from_seed(ring, right_seed);

    benchmark::DoNotOptimize(left);
    benchmark::DoNotOptimize(right);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_mul(left, right, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_Mul_50bit");
}

static void BM_FP12_Sqr_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_sqr(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_Sqr_50bit");
}

static void BM_FP12_CyclotomicSqr_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_cyclotomic_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_cyclotomic_sqr(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_CyclotomicSqr_50bit");
}

static void BM_FP12_Conjugate_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_conjugate(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_Conjugate_50bit");
}

static void BM_FP12_FrobeniusMap1_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_frobenius_map1(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FrobeniusMap1_50bit");
}

static void BM_FP12_FrobeniusMap2_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_frobenius_map2(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FrobeniusMap2_50bit");
}

static void BM_FP12_FrobeniusMap3_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_frobenius_map3(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FrobeniusMap3_50bit");
}

static void BM_FP12_MulByXy0z00_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    volatile uint64_t value_seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, value_seed);
    const auto sparse = blst_vroom_test::make_vroom_sparse_fp6_from_seed(ring, 123);
    const auto c_term = blst_vroom_test::make_vroom_fp2_from_sparse6<Vroom50Ring>(sparse, 0);
    const auto x_term = blst_vroom_test::make_vroom_fp2_from_sparse6<Vroom50Ring>(sparse, 1);
    const auto y_term = blst_vroom_test::make_vroom_fp2_from_sparse6<Vroom50Ring>(sparse, 2);

    benchmark::DoNotOptimize(value);
    benchmark::DoNotOptimize(sparse);

    for (auto _ : state) {
        Vroom50Fp12 result;
        fp12ring.flat_mulxy00z0(result, value, c_term, x_term, y_term, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_MulByXy0z00_50bit");
}

static void BM_FP12_FinalExpAfterInverse_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    volatile uint64_t seed_f = 42;
    volatile uint64_t seed_f_inverse = 123;
    const auto f = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed_f);
    const auto f_inverse = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed_f_inverse);

    benchmark::DoNotOptimize(f);
    benchmark::DoNotOptimize(f_inverse);

    for (auto _ : state) {
        const auto result = final_exp_after_inverse(f, f_inverse, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FinalExpAfterInverse_50bit");
}

BENCHMARK(BM_FP12_Mul_50bit);
BENCHMARK(BM_FP12_Sqr_50bit);
BENCHMARK(BM_FP12_CyclotomicSqr_50bit);
BENCHMARK(BM_FP12_Conjugate_50bit);
BENCHMARK(BM_FP12_FrobeniusMap1_50bit);
BENCHMARK(BM_FP12_FrobeniusMap2_50bit);
BENCHMARK(BM_FP12_FrobeniusMap3_50bit);
BENCHMARK(BM_FP12_MulByXy0z00_50bit);
BENCHMARK(BM_FP12_FinalExpAfterInverse_50bit);

}  // namespace
