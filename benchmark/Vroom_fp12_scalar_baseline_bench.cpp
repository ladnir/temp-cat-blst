#include <benchmark/benchmark.h>

#include "blst/avx/tests/avx_test_utils.h"

namespace {

static void BM_FP12_ScalarMulBaseline(benchmark::State& state)
{
    const blst::Fp12 left = blst_vroom_test::make_scalar_fp12_from_seed(42);
    const blst::Fp12 right = blst_vroom_test::make_scalar_fp12_from_seed(123);

    benchmark::DoNotOptimize(left);
    benchmark::DoNotOptimize(right);

    for (auto _ : state) {
        blst::Fp12 result = {};
        blst::Fp12::mul(result, left, right);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_ScalarMulBaseline");
}

BENCHMARK(BM_FP12_ScalarMulBaseline);

}  // namespace
