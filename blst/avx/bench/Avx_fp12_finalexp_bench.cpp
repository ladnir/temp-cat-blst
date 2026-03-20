#include <benchmark/benchmark.h>

#include "blst/avx/fp12/final_exponentiation.h"
#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_avx_test::AvxFp12Ring;
using blst_avx_test::AvxRing;

static void BM_FP12_FinalExpAfterInverse(benchmark::State& state)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());
    const AvxFp12Ring fp12ring(ring);
    volatile uint64_t seed_f = 42;
    volatile uint64_t seed_inv = 123;
    const auto value = blst_avx_test::make_avx_fp12_from_seed(ring, seed_f);
    const auto inverse = blst_avx_test::make_avx_fp12_from_seed(ring, seed_inv);

    benchmark::DoNotOptimize(value);
    benchmark::DoNotOptimize(inverse);

    for (auto _ : state) {
        const auto result = final_exp_after_inverse(value, inverse, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FinalExpAfterInverse");
}

BENCHMARK(BM_FP12_FinalExpAfterInverse);

}  // namespace
