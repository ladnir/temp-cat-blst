#include <benchmark/benchmark.h>

#include "blst/avx/fields/conversion_inversion.h"
#include "test/avx/avx_test_utils.h"

using namespace blst::avx;

namespace {

using blst_avx_test::AvxFp12Ring;
using blst_avx_test::AvxRing;

static void BM_FP12_Inverse_BLST(benchmark::State& state)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());
    volatile uint64_t seed = 42;
    const auto value = blst_avx_test::make_avx_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = invert_fp12_via_blst(ring, value);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_Inverse_BLST");
}

static void BM_FP12_Inverse_RNS_BLST(benchmark::State& state)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());
    const AvxFp12Ring fp12ring(ring);
    const FP12RNS_BLSTInverter<AvxRing> inverter = {};
    volatile uint64_t seed = 42;
    const auto value = blst_avx_test::make_avx_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = inverter.invert(value, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_Inverse_RNS_BLST");
}

BENCHMARK(BM_FP12_Inverse_BLST);
BENCHMARK(BM_FP12_Inverse_RNS_BLST);

}  // namespace
