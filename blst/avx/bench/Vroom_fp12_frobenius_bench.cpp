#include <benchmark/benchmark.h>

#include "blst/avx/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomRing;

static void BM_FP12_Conjugate(benchmark::State& state)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_conjugate(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_Conjugate");
}

static void BM_FP12_FrobeniusMap1(benchmark::State& state)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_frobenius_map1(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FrobeniusMap1");
}

static void BM_FP12_FrobeniusMap2(benchmark::State& state)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_frobenius_map2(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FrobeniusMap2");
}

static void BM_FP12_FrobeniusMap3(benchmark::State& state)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_frobenius_map3(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FrobeniusMap3");
}

BENCHMARK(BM_FP12_Conjugate);
BENCHMARK(BM_FP12_FrobeniusMap1);
BENCHMARK(BM_FP12_FrobeniusMap2);
BENCHMARK(BM_FP12_FrobeniusMap3);

}  // namespace
