#include <benchmark/benchmark.h>

#include "blst/avx/fp12/conversion_inversion.h"
#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_vroom_test::Vroom50Fp12Ring;
using blst_vroom_test::Vroom50Ring;

static void BM_FP12_Inverse_BLST_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = invert_fp12_via_blst(ring, value);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_Inverse_BLST_50bit");
}

static void BM_FP12_Inverse_RNS_BLST_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    const FP12RNS_BLSTInverter<Vroom50Ring> inverter = {};
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = inverter.invert(value, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_Inverse_RNS_BLST_50bit");
}

BENCHMARK(BM_FP12_Inverse_BLST_50bit);
BENCHMARK(BM_FP12_Inverse_RNS_BLST_50bit);

}  // namespace
