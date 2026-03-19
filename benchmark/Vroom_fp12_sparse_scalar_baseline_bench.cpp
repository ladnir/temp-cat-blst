#include <benchmark/benchmark.h>

#include "blst/avx/tests/vroom_basefield_test_utils.hpp"

namespace {

static void BM_FP12_ScalarMulByXy0z00Baseline(benchmark::State& state)
{
    const blst::Fp12 value = blst_vroom_test::make_scalar_fp12_from_seed(42);
    const blst::Fp6 sparse = blst_vroom_test::make_scalar_sparse_fp6_from_seed(123);

    benchmark::DoNotOptimize(value);
    benchmark::DoNotOptimize(sparse);

    for (auto _ : state) {
        blst::Fp12 result = {};
        blst::Fp12::mul_by_xy00z0(result, value, sparse);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_ScalarMulByXy0z00Baseline");
}

BENCHMARK(BM_FP12_ScalarMulByXy0z00Baseline);

}  // namespace
