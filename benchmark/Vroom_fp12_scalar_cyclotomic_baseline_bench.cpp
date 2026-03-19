#include <benchmark/benchmark.h>

#include "blst/avx/tests/vroom_basefield_test_utils.hpp"

namespace {

static void BM_FP12_ScalarCyclotomicSqrBaseline(benchmark::State& state)
{
    const blst::Fp12 value = blst_vroom_test::make_scalar_cyclotomic_fp12_from_seed(42);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        blst::Fp12 result = {};
        blst::Fp12::cyclotomic_sqr(result, value);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_ScalarCyclotomicSqrBaseline");
}

BENCHMARK(BM_FP12_ScalarCyclotomicSqrBaseline);

}  // namespace
