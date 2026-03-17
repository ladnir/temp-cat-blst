#include <benchmark/benchmark.h>

#include "blst/vroom/src/miller.hpp"
#include "blst/vroom/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::VroomAffineFp;
using blst_vroom_test::VroomAffineFp2;
using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomFp2Ring;
using blst_vroom_test::VroomProjectiveFp2;
using blst_vroom_test::VroomRing;

static void BM_MillerLoop_DoubleStep(benchmark::State& state)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const VroomRing ring(modulus);
    const VroomFp2Ring fp2ring(ring);
    const VroomFp12Ring fp12ring(ring);
    const VroomAffineFp g1 = blst_vroom_test::make_vroom_g1_generator_affine(ring);
    const VroomAffineFp2 g2 = blst_vroom_test::make_vroom_g2_generator_affine(ring);
    Miller<VroomRing> miller(g1, g2, ring);

    VroomProjectiveFp2 point(g2.x, g2.y, VroomFp2Ring::one(ring));
    auto f = fp12ring.one(ring);

    benchmark::DoNotOptimize(point);
    benchmark::DoNotOptimize(f);

    for (auto _ : state) {
        point = miller.double_step(point, f, ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(point);
        benchmark::DoNotOptimize(f);
    }

    state.SetLabel("MillerLoop_DoubleStep");
}

static void BM_MillerLoop_AddStep(benchmark::State& state)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const VroomRing ring(modulus);
    const VroomFp2Ring fp2ring(ring);
    const VroomFp12Ring fp12ring(ring);
    const VroomAffineFp g1 = blst_vroom_test::make_vroom_g1_generator_affine(ring);
    const VroomAffineFp2 g2 = blst_vroom_test::make_vroom_g2_generator_affine(ring);
    Miller<VroomRing> miller(g1, g2, ring);

    VroomProjectiveFp2 point(g2.x, g2.y, VroomFp2Ring::one(ring));
    auto f = fp12ring.one(ring);
    point = miller.double_step(point, f, ring, fp2ring, fp12ring);

    benchmark::DoNotOptimize(point);
    benchmark::DoNotOptimize(f);

    for (auto _ : state) {
        point = miller.add_step(point, f, ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(point);
        benchmark::DoNotOptimize(f);
    }

    state.SetLabel("MillerLoop_AddStep");
}

static void BM_MillerLoop(benchmark::State& state)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const VroomRing ring(modulus);
    const VroomFp2Ring fp2ring(ring);
    const VroomFp12Ring fp12ring(ring);
    const VroomAffineFp g1 = blst_vroom_test::make_vroom_g1_generator_affine(ring);
    const VroomAffineFp2 g2 = blst_vroom_test::make_vroom_g2_generator_affine(ring);
    Miller<VroomRing> miller(g1, g2, ring);

    for (auto _ : state) {
        auto result = miller.miller_loop(ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("MillerLoop");
}

BENCHMARK(BM_MillerLoop_DoubleStep);
BENCHMARK(BM_MillerLoop_AddStep);
BENCHMARK(BM_MillerLoop);

}  // namespace
