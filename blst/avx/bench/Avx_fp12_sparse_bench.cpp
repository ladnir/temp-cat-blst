#include <benchmark/benchmark.h>

#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_avx_test::AvxFp12Ring;
using blst_avx_test::AvxRing;

static void BM_FP12_MulByXy0z00(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp12Ring fp12ring(ring);
    volatile uint64_t x_seed = 42;
    const auto value = blst_avx_test::make_avx_fp12_from_seed(ring, x_seed);

    const auto sparse_coeffs = blst_avx_test::make_avx_sparse_fp6_from_seed(ring, 123);
    const auto c_term = blst_avx_test::make_avx_fp2_from_sparse6<AvxRing>(sparse_coeffs, 0);
    const auto x_term = blst_avx_test::make_avx_fp2_from_sparse6<AvxRing>(sparse_coeffs, 1);
    const auto y_term = blst_avx_test::make_avx_fp2_from_sparse6<AvxRing>(sparse_coeffs, 2);

    benchmark::DoNotOptimize(value);
    benchmark::DoNotOptimize(sparse_coeffs);

    for (auto _ : state) {
        AvxFp12Ring::StandardElement result = {};
        fp12ring.flat_mulxy00z0(result, value, c_term, x_term, y_term, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_MulByXy0z00");
}

BENCHMARK(BM_FP12_MulByXy0z00);

}  // namespace
