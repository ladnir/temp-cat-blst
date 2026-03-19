#include <benchmark/benchmark.h>

#include "blst/avx/src/miller.hpp"
#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_avx_test::Avx50Fp12Ring;
using blst_avx_test::Avx50Fp2Ring;
using blst_avx_test::Avx50Ring;

static void BM_MillerLoop_DoubleStep_50bit(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const Avx50Ring ring(modulus);
    const Avx50Fp2Ring fp2ring(ring);
    const Avx50Fp12Ring fp12ring(ring);
    const auto g1 = blst_avx_test::make_vroom_g1_generator_affine(ring);
    const auto g2 = blst_avx_test::make_vroom_g2_generator_affine(ring);
    Miller<Avx50Ring> miller(g1, g2, ring);

    ProjectivePoint<Avx50Fp2Ring::StandardElement> point(g2.x, g2.y, Avx50Fp2Ring::one(ring));
    auto f = fp12ring.one(ring);

    benchmark::DoNotOptimize(point);
    benchmark::DoNotOptimize(f);

    for (auto _ : state) {
        point = miller.double_step(point, f, ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(point);
        benchmark::DoNotOptimize(f);
    }

    state.SetLabel("MillerLoop_DoubleStep_50bit");
}

static void BM_MillerLoop_AddStep_50bit(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const Avx50Ring ring(modulus);
    const Avx50Fp2Ring fp2ring(ring);
    const Avx50Fp12Ring fp12ring(ring);
    const auto g1 = blst_avx_test::make_vroom_g1_generator_affine(ring);
    const auto g2 = blst_avx_test::make_vroom_g2_generator_affine(ring);
    Miller<Avx50Ring> miller(g1, g2, ring);

    ProjectivePoint<Avx50Fp2Ring::StandardElement> point(g2.x, g2.y, Avx50Fp2Ring::one(ring));
    auto f = fp12ring.one(ring);
    point = miller.double_step(point, f, ring, fp2ring, fp12ring);

    benchmark::DoNotOptimize(point);
    benchmark::DoNotOptimize(f);

    for (auto _ : state) {
        point = miller.add_step(point, f, ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(point);
        benchmark::DoNotOptimize(f);
    }

    state.SetLabel("MillerLoop_AddStep_50bit");
}

static void BM_MillerLoop_50bit(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const Avx50Ring ring(modulus);
    const Avx50Fp2Ring fp2ring(ring);
    const Avx50Fp12Ring fp12ring(ring);
    const auto g1 = blst_avx_test::make_vroom_g1_generator_affine(ring);
    const auto g2 = blst_avx_test::make_vroom_g2_generator_affine(ring);
    Miller<Avx50Ring> miller(g1, g2, ring);

    for (auto _ : state) {
        auto result = miller.miller_loop(ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("MillerLoop_50bit");
}

BENCHMARK(BM_MillerLoop_DoubleStep_50bit);
BENCHMARK(BM_MillerLoop_AddStep_50bit);
BENCHMARK(BM_MillerLoop_50bit);

}  // namespace
