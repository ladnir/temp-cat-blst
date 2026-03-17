#include <benchmark/benchmark.h>

#include "blst/vroom/src/pairing.hpp"
#include "blst/vroom/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::VroomAffineFp;
using blst_vroom_test::VroomAffineFp2;
using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomFp2Ring;
using blst_vroom_test::VroomRing;

static void BM_Pairing_BLST_Inverter(benchmark::State& state)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const VroomRing ring(modulus);
    const VroomFp2Ring fp2ring(ring);
    const VroomFp12Ring fp12ring(ring);
    const VroomAffineFp g1 = blst_vroom_test::make_vroom_g1_generator_affine(ring);
    const VroomAffineFp2 g2 = blst_vroom_test::make_vroom_g2_generator_affine(ring);
    FP12BLSTInverter<VroomRing> inverter;

    for (auto _ : state) {
        auto result = pairing(g1, g2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_BLST_Inverter");
}

static void BM_Pairing_RNS_BLST_Inverter(benchmark::State& state)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const VroomRing ring(modulus);
    const VroomFp2Ring fp2ring(ring);
    const VroomFp12Ring fp12ring(ring);
    const VroomAffineFp g1 = blst_vroom_test::make_vroom_g1_generator_affine(ring);
    const VroomAffineFp2 g2 = blst_vroom_test::make_vroom_g2_generator_affine(ring);
    FP12RNS_BLSTInverter<VroomRing> inverter;

    for (auto _ : state) {
        auto result = pairing(g1, g2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_RNS_BLST_Inverter");
}

BENCHMARK(BM_Pairing_BLST_Inverter);
BENCHMARK(BM_Pairing_RNS_BLST_Inverter);

}  // namespace
