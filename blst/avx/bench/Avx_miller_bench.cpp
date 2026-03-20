#include <benchmark/benchmark.h>

#include "blst/avx/pairing/miller.h"
#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_avx_test::AvxAffineFp;
using blst_avx_test::AvxAffineFp2;
using blst_avx_test::AvxFp12Ring;
using blst_avx_test::AvxFp2Ring;
using blst_avx_test::AvxProjectiveFp2;
using blst_avx_test::AvxRing;

static void BM_MillerLoop_DoubleStep(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp2Ring fp2ring(ring);
    const AvxFp12Ring fp12ring(ring);
    const AvxAffineFp g1 = blst_avx_test::make_avx_g1_generator_affine(ring);
    const AvxAffineFp2 g2 = blst_avx_test::make_avx_g2_generator_affine(ring);
    Miller<AvxRing> miller(g1, g2, ring);

    AvxProjectiveFp2 point(g2.x, g2.y, AvxFp2Ring::one(ring));
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
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp2Ring fp2ring(ring);
    const AvxFp12Ring fp12ring(ring);
    const AvxAffineFp g1 = blst_avx_test::make_avx_g1_generator_affine(ring);
    const AvxAffineFp2 g2 = blst_avx_test::make_avx_g2_generator_affine(ring);
    Miller<AvxRing> miller(g1, g2, ring);

    AvxProjectiveFp2 point(g2.x, g2.y, AvxFp2Ring::one(ring));
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
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp2Ring fp2ring(ring);
    const AvxFp12Ring fp12ring(ring);
    const AvxAffineFp g1 = blst_avx_test::make_avx_g1_generator_affine(ring);
    const AvxAffineFp2 g2 = blst_avx_test::make_avx_g2_generator_affine(ring);
    Miller<AvxRing> miller(g1, g2, ring);

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
