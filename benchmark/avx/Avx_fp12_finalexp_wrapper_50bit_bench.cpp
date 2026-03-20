#include <benchmark/benchmark.h>

#include "blst/avx/fp12/conversion_inversion.h"
#include "blst/avx/fp12/final_exponentiation.h"
#include "test/avx/avx_test_utils.h"

using namespace blst::avx;

namespace {

using blst_avx_test::Avx50Fp12Ring;
using blst_avx_test::Avx50Ring;

static void BM_FP12_FinalExp_BLST_Inverter_50bit(benchmark::State& state)
{
    const Avx50Ring ring(blst_avx_test::modulus_bigint());
    const Avx50Fp12Ring fp12ring(ring);
    const FP12BLSTInverter<Avx50Ring> inverter = {};
    volatile uint64_t seed = 42;
    const auto value = blst_avx_test::make_avx_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = final_exp(value, inverter, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FinalExp_BLST_Inverter_50bit");
}

static void BM_FP12_FinalExp_RNS_BLST_Inverter_50bit(benchmark::State& state)
{
    const Avx50Ring ring(blst_avx_test::modulus_bigint());
    const Avx50Fp12Ring fp12ring(ring);
    const FP12RNS_BLSTInverter<Avx50Ring> inverter = {};
    volatile uint64_t seed = 42;
    const auto value = blst_avx_test::make_avx_fp12_from_seed(ring, seed);

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
