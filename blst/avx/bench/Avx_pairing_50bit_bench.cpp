#include <benchmark/benchmark.h>

#include "blst/avx/src/pairing.hpp"
#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_avx_test::Avx50Fp12Ring;
using blst_avx_test::Avx50Fp2Ring;
using blst_avx_test::Avx50Ring;

static void BM_Pairing_BLST_Inverter_50bit(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const Avx50Ring ring(modulus);
    const Avx50Fp2Ring fp2ring(ring);
    const Avx50Fp12Ring fp12ring(ring);
    const auto g1 = blst_avx_test::make_avx_g1_generator_affine(ring);
    const auto g2 = blst_avx_test::make_avx_g2_generator_affine(ring);
    FP12BLSTInverter<Avx50Ring> inverter;

    for (auto _ : state) {
        const auto result = pairing(g1, g2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_BLST_Inverter_50bit");
}

static void BM_Pairing_RNS_BLST_Inverter_50bit(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const Avx50Ring ring(modulus);
    const Avx50Fp2Ring fp2ring(ring);
    const Avx50Fp12Ring fp12ring(ring);
    const auto g1 = blst_avx_test::make_avx_g1_generator_affine(ring);
    const auto g2 = blst_avx_test::make_avx_g2_generator_affine(ring);
    FP12RNS_BLSTInverter<Avx50Ring> inverter;

    for (auto _ : state) {
        const auto result = pairing(g1, g2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_RNS_BLST_Inverter_50bit");
}

BENCHMARK(BM_Pairing_BLST_Inverter_50bit);
BENCHMARK(BM_Pairing_RNS_BLST_Inverter_50bit);

}  // namespace
