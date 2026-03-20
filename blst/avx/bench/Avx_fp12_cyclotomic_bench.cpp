#include <benchmark/benchmark.h>

#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_avx_test::AvxFp12Ring;
using blst_avx_test::AvxRing;

static void BM_FP12_CyclotomicSqr(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_avx_test::make_avx_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_cyclotomic_sqr(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_CyclotomicSqr");
}

BENCHMARK(BM_FP12_CyclotomicSqr);

}  // namespace
