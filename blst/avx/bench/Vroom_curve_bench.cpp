#include <benchmark/benchmark.h>

#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_vroom_test::VroomCurveFp2Ring;
using blst_vroom_test::VroomCurveG1;
using blst_vroom_test::VroomCurveG2;
using blst_vroom_test::VroomCurveRing;

static void BM_G1_Add(benchmark::State& state)
{
    const VroomCurveRing ring(blst_vroom_test::modulus_bigint());
    const VroomCurveG1 g1 = {};
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    const auto p = blst_vroom_test::make_random_vroom_curve_g1_projective(ring, seed1);
    const auto q = blst_vroom_test::make_random_vroom_curve_g1_projective(ring, seed2);

    benchmark::DoNotOptimize(p);
    benchmark::DoNotOptimize(q);

    for (auto _ : state) {
        const auto result = g1.add_point(p, q, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G1_Add");
}

static void BM_G1_MixedAdd(benchmark::State& state)
{
    const VroomCurveRing ring(blst_vroom_test::modulus_bigint());
    const VroomCurveG1 g1 = {};
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    const auto p = blst_vroom_test::make_random_vroom_curve_g1_projective(ring, seed1);
    const auto q_affine = blst_vroom_test::make_random_vroom_curve_g1_affine(ring, seed2);

    benchmark::DoNotOptimize(p);
    benchmark::DoNotOptimize(q_affine);

    for (auto _ : state) {
        const auto result = g1.add_mixed_point(p, q_affine, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G1_MixedAdd");
}

static void BM_G1_Double(benchmark::State& state)
{
    const VroomCurveRing ring(blst_vroom_test::modulus_bigint());
    const VroomCurveG1 g1 = {};
    volatile uint64_t seed = 42;
    const auto p = blst_vroom_test::make_random_vroom_curve_g1_projective(ring, seed);

    benchmark::DoNotOptimize(p);

    for (auto _ : state) {
        const auto result = g1.double_point(p, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G1_Double");
}

static void BM_G2_Add(benchmark::State& state)
{
    const VroomCurveRing ring(blst_vroom_test::modulus_bigint());
    const VroomCurveFp2Ring fp2ring(ring);
    const VroomCurveG2 g2 = {};
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    const auto p = blst_vroom_test::make_random_vroom_curve_g2_projective(fp2ring, ring, seed1);
    const auto q = blst_vroom_test::make_random_vroom_curve_g2_projective(fp2ring, ring, seed2);

    benchmark::DoNotOptimize(p);
    benchmark::DoNotOptimize(q);

    for (auto _ : state) {
        const auto result = g2.add_point(p, q, fp2ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G2_Add");
}

static void BM_G2_MixedAdd(benchmark::State& state)
{
    const VroomCurveRing ring(blst_vroom_test::modulus_bigint());
    const VroomCurveFp2Ring fp2ring(ring);
    const VroomCurveG2 g2 = {};
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    const auto p = blst_vroom_test::make_random_vroom_curve_g2_projective(fp2ring, ring, seed1);
    const auto q_affine = blst_vroom_test::make_random_vroom_curve_g2_affine(ring, seed2);

    benchmark::DoNotOptimize(p);
    benchmark::DoNotOptimize(q_affine);

    for (auto _ : state) {
        const auto result = g2.add_mixed_point(p, q_affine, fp2ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G2_MixedAdd");
}

static void BM_G2_Double(benchmark::State& state)
{
    const VroomCurveRing ring(blst_vroom_test::modulus_bigint());
    const VroomCurveFp2Ring fp2ring(ring);
    const VroomCurveG2 g2 = {};
    volatile uint64_t seed = 42;
    const auto p = blst_vroom_test::make_random_vroom_curve_g2_projective(fp2ring, ring, seed);

    benchmark::DoNotOptimize(p);

    for (auto _ : state) {
        const auto result = g2.double_point(p, fp2ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G2_Double");
}

BENCHMARK(BM_G1_Add);
BENCHMARK(BM_G1_MixedAdd);
BENCHMARK(BM_G1_Double);
BENCHMARK(BM_G2_Add);
BENCHMARK(BM_G2_MixedAdd);
BENCHMARK(BM_G2_Double);

}  // namespace
