#include <benchmark/benchmark.h>

#include "blst/vroom/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::VroomCurve50Fp2Ring;
using blst_vroom_test::VroomCurve50G1;
using blst_vroom_test::VroomCurve50G2;
using blst_vroom_test::VroomCurve50Ring;

static void BM_G1_Add_50bit(benchmark::State& state)
{
    const VroomCurve50Ring ring(blst_vroom_test::modulus_bigint());
    const VroomCurve50G1 g1 = {};
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    const auto p = blst_vroom_test::make_random_vroom_curve50_g1_projective(ring, seed1);
    const auto q = blst_vroom_test::make_random_vroom_curve50_g1_projective(ring, seed2);

    benchmark::DoNotOptimize(p);
    benchmark::DoNotOptimize(q);

    for (auto _ : state) {
        const auto result = g1.add_point(p, q, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G1_Add_50bit");
}

static void BM_G1_MixedAdd_50bit(benchmark::State& state)
{
    const VroomCurve50Ring ring(blst_vroom_test::modulus_bigint());
    const VroomCurve50G1 g1 = {};
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    const auto p = blst_vroom_test::make_random_vroom_curve50_g1_projective(ring, seed1);
    const auto q_affine = blst_vroom_test::make_random_vroom_curve50_g1_affine(ring, seed2);

    benchmark::DoNotOptimize(p);
    benchmark::DoNotOptimize(q_affine);

    for (auto _ : state) {
        const auto result = g1.add_mixed_point(p, q_affine, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G1_MixedAdd_50bit");
}

static void BM_G1_Double_50bit(benchmark::State& state)
{
    const VroomCurve50Ring ring(blst_vroom_test::modulus_bigint());
    const VroomCurve50G1 g1 = {};
    volatile uint64_t seed = 42;
    const auto p = blst_vroom_test::make_random_vroom_curve50_g1_projective(ring, seed);

    benchmark::DoNotOptimize(p);

    for (auto _ : state) {
        const auto result = g1.double_point(p, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G1_Double_50bit");
}

static void BM_G2_Add_50bit(benchmark::State& state)
{
    const VroomCurve50Ring ring(blst_vroom_test::modulus_bigint());
    const VroomCurve50Fp2Ring fp2ring(ring);
    const VroomCurve50G2 g2 = {};
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    const auto p = blst_vroom_test::make_random_vroom_curve50_g2_projective(fp2ring, ring, seed1);
    const auto q = blst_vroom_test::make_random_vroom_curve50_g2_projective(fp2ring, ring, seed2);

    benchmark::DoNotOptimize(p);
    benchmark::DoNotOptimize(q);

    for (auto _ : state) {
        const auto result = g2.add_point(p, q, fp2ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G2_Add_50bit");
}

static void BM_G2_MixedAdd_50bit(benchmark::State& state)
{
    const VroomCurve50Ring ring(blst_vroom_test::modulus_bigint());
    const VroomCurve50Fp2Ring fp2ring(ring);
    const VroomCurve50G2 g2 = {};
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    const auto p = blst_vroom_test::make_random_vroom_curve50_g2_projective(fp2ring, ring, seed1);
    const auto q_affine = blst_vroom_test::make_random_vroom_curve50_g2_affine(ring, seed2);

    benchmark::DoNotOptimize(p);
    benchmark::DoNotOptimize(q_affine);

    for (auto _ : state) {
        const auto result = g2.add_mixed_point(p, q_affine, fp2ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G2_MixedAdd_50bit");
}

static void BM_G2_Double_50bit(benchmark::State& state)
{
    const VroomCurve50Ring ring(blst_vroom_test::modulus_bigint());
    const VroomCurve50Fp2Ring fp2ring(ring);
    const VroomCurve50G2 g2 = {};
    volatile uint64_t seed = 42;
    const auto p = blst_vroom_test::make_random_vroom_curve50_g2_projective(fp2ring, ring, seed);

    benchmark::DoNotOptimize(p);

    for (auto _ : state) {
        const auto result = g2.double_point(p, fp2ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G2_Double_50bit");
}

BENCHMARK(BM_G1_Add_50bit);
BENCHMARK(BM_G1_MixedAdd_50bit);
BENCHMARK(BM_G1_Double_50bit);
BENCHMARK(BM_G2_Add_50bit);
BENCHMARK(BM_G2_MixedAdd_50bit);
BENCHMARK(BM_G2_Double_50bit);

}  // namespace
