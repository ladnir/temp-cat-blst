#include <benchmark/benchmark.h>

#include "blst/avx/src/pairing.hpp"
#include "blst/avx/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::Vroom50Fp12Ring;
using blst_vroom_test::Vroom50Fp2Ring;
using blst_vroom_test::Vroom50Ring;

static void BM_Pairing_BLST_Inverter_50bit(benchmark::State& state)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const Vroom50Ring ring(modulus);
    const Vroom50Fp2Ring fp2ring(ring);
    const Vroom50Fp12Ring fp12ring(ring);
    const auto g1 = blst_vroom_test::make_vroom_g1_generator_affine(ring);
    const auto g2 = blst_vroom_test::make_vroom_g2_generator_affine(ring);
    FP12BLSTInverter<Vroom50Ring> inverter;

    for (auto _ : state) {
        const auto result = pairing(g1, g2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_BLST_Inverter_50bit");
}

static void BM_Pairing_RNS_BLST_Inverter_50bit(benchmark::State& state)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const Vroom50Ring ring(modulus);
    const Vroom50Fp2Ring fp2ring(ring);
    const Vroom50Fp12Ring fp12ring(ring);
    const auto g1 = blst_vroom_test::make_vroom_g1_generator_affine(ring);
    const auto g2 = blst_vroom_test::make_vroom_g2_generator_affine(ring);
    FP12RNS_BLSTInverter<Vroom50Ring> inverter;

    for (auto _ : state) {
        const auto result = pairing(g1, g2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_RNS_BLST_Inverter_50bit");
}

BENCHMARK(BM_Pairing_BLST_Inverter_50bit);
BENCHMARK(BM_Pairing_RNS_BLST_Inverter_50bit);

}  // namespace
