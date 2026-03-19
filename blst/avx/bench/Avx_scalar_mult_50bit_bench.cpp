#include <benchmark/benchmark.h>

#define BLST_VROOM_SKIP_TEST_VALUES
#include "blst/avx/src/scalar_mult.hpp"
#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_avx_test::AvxCurve50Fp2Ring;
using blst_avx_test::AvxCurve50G1;
using blst_avx_test::AvxCurve50G2;
using blst_avx_test::AvxCurve50Ring;

static void BM_G1_ScalarMult_GLV_50bit(benchmark::State& state)
{
    const AvxCurve50Ring ring(blst_avx_test::modulus_bigint());
    const AvxCurve50G1 g1 = {};
    const GLVMult<AvxCurve50Ring> glv(g1, ring);

    volatile uint64_t seed_point = 42;
    const auto point = blst_avx_test::make_random_vroom_curve50_g1_projective(ring, seed_point);
    const auto scalar_bytes = blst_avx_test::seeded_scalar_bytes_le(123);

    benchmark::DoNotOptimize(point);
    benchmark::DoNotOptimize(scalar_bytes);

    for (auto _ : state) {
        const auto result = glv.multiply(point, scalar_bytes, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G1_ScalarMult_GLV_50bit");
}

static void BM_G2_ScalarMult_GLS_50bit(benchmark::State& state)
{
    const AvxCurve50Ring ring(blst_avx_test::modulus_bigint());
    const AvxCurve50Fp2Ring fp2ring(ring);
    const AvxCurve50G2 g2 = {};
    const GLSMult<AvxCurve50Ring> gls(ring);

    volatile uint64_t seed_point = 42;
    const auto point = blst_avx_test::make_random_vroom_curve50_g2_projective(fp2ring, ring, seed_point);
    const auto scalar_bytes = blst_avx_test::seeded_scalar_bytes_le(321);

    benchmark::DoNotOptimize(point);
    benchmark::DoNotOptimize(scalar_bytes);

    for (auto _ : state) {
        const auto result = gls.multiply(point, scalar_bytes, fp2ring, g2, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G2_ScalarMult_GLS_50bit");
}

BENCHMARK(BM_G1_ScalarMult_GLV_50bit);
BENCHMARK(BM_G2_ScalarMult_GLS_50bit);

}  // namespace
