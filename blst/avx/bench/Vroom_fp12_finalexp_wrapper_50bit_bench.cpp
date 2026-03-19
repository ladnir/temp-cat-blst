#include <benchmark/benchmark.h>

#include "blst/avx/fp12/conversion_inversion.h"
#include "blst/avx/fp12/final_exponentiation.h"
#include "blst/avx/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::Vroom50Fp12Ring;
using blst_vroom_test::Vroom50Ring;

static void BM_FP12_FinalExp_BLST_Inverter_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    const FP12BLSTInverter<Vroom50Ring> inverter = {};
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = final_exp(value, inverter, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FinalExp_BLST_Inverter_50bit");
}

static void BM_FP12_FinalExp_RNS_BLST_Inverter_50bit(benchmark::State& state)
{
    const Vroom50Ring ring(blst_vroom_test::modulus_bigint());
    const Vroom50Fp12Ring fp12ring(ring);
    const FP12RNS_BLSTInverter<Vroom50Ring> inverter = {};
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = final_exp(value, inverter, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FinalExp_RNS_BLST_Inverter_50bit");
}

BENCHMARK(BM_FP12_FinalExp_BLST_Inverter_50bit);
BENCHMARK(BM_FP12_FinalExp_RNS_BLST_Inverter_50bit);

}  // namespace
