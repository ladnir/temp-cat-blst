#include <benchmark/benchmark.h>

#include "blst/avx/pairing/pairing.h"
#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_avx_test::AvxAffineFp;
using blst_avx_test::AvxAffineFp2;
using blst_avx_test::AvxFp12Ring;
using blst_avx_test::AvxFp2Ring;
using blst_avx_test::AvxRing;

static void BM_Pairing_BLST_Inverter(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp2Ring fp2ring(ring);
    const AvxFp12Ring fp12ring(ring);
    const AvxAffineFp g1 = blst_avx_test::make_avx_g1_generator_affine(ring);
    const AvxAffineFp2 g2 = blst_avx_test::make_avx_g2_generator_affine(ring);
    FP12BLSTInverter<AvxRing> inverter;

    for (auto _ : state) {
        auto result = pairing(g1, g2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_BLST_Inverter");
}

static void BM_Pairing_RNS_BLST_Inverter(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp2Ring fp2ring(ring);
    const AvxFp12Ring fp12ring(ring);
    const AvxAffineFp g1 = blst_avx_test::make_avx_g1_generator_affine(ring);
    const AvxAffineFp2 g2 = blst_avx_test::make_avx_g2_generator_affine(ring);
    FP12RNS_BLSTInverter<AvxRing> inverter;

    for (auto _ : state) {
        auto result = pairing(g1, g2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_RNS_BLST_Inverter");
}

BENCHMARK(BM_Pairing_BLST_Inverter);
BENCHMARK(BM_Pairing_RNS_BLST_Inverter);

}  // namespace
