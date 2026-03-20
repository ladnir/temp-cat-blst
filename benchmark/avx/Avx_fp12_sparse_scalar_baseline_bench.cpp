#include <benchmark/benchmark.h>

#include "test/avx/avx_test_utils.h"

using namespace blst::avx;

namespace {

static void BM_FP12_ScalarMulByXy0z00Baseline(benchmark::State& state)
{
    const blst::Fp12 value = blst_avx_test::make_scalar_fp12_from_seed(42);
    const blst::Fp6 sparse = blst_avx_test::make_scalar_sparse_fp6_from_seed(123);

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
