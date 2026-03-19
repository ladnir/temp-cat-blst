#include <benchmark/benchmark.h>

#include "blst/avx/src/conversion_inversion.hpp"
#include "blst/avx/src/final_exponentiation.hpp"
#include "blst/avx/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomRing;

static void BM_FP12_FinalExp_BLST_Inverter(benchmark::State& state)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const FP12BLSTInverter<VroomRing> inverter = {};
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = final_exp(value, inverter, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FinalExp_BLST_Inverter");
}

static void BM_FP12_FinalExp_RNS_BLST_Inverter(benchmark::State& state)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const FP12RNS_BLSTInverter<VroomRing> inverter = {};
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = final_exp(value, inverter, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FinalExp_RNS_BLST_Inverter");
}

BENCHMARK(BM_FP12_FinalExp_BLST_Inverter);
BENCHMARK(BM_FP12_FinalExp_RNS_BLST_Inverter);

}  // namespace
