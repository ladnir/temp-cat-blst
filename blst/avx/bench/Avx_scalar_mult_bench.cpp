#include <benchmark/benchmark.h>

#define BLST_VROOM_SKIP_TEST_VALUES
#include "blst/avx/src/scalar_mult.hpp"
#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_avx_test::AvxCurveFp2Ring;
using blst_avx_test::AvxCurveG1;
using blst_avx_test::AvxCurveG2;
using blst_avx_test::AvxCurveRing;

static void BM_G1_ScalarMult_GLV(benchmark::State& state)
{
    const AvxCurveRing ring(blst_avx_test::modulus_bigint());
    const AvxCurveG1 g1 = {};
    const GLVMult<AvxCurveRing> glv(g1, ring);

    volatile uint64_t seed_point = 42;
    const auto point = blst_avx_test::make_random_vroom_curve_g1_projective(ring, seed_point);
    const auto scalar_bytes = blst_avx_test::seeded_scalar_bytes_le(123);

    benchmark::DoNotOptimize(point);
    benchmark::DoNotOptimize(scalar_bytes);

    for (auto _ : state) {
        const auto result = glv.multiply(point, scalar_bytes, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G1_ScalarMult_GLV");
}

static void BM_G2_ScalarMult_GLS(benchmark::State& state)
{
    const AvxCurveRing ring(blst_avx_test::modulus_bigint());
    const AvxCurveFp2Ring fp2ring(ring);
    const AvxCurveG2 g2 = {};
    const GLSMult<AvxCurveRing> gls(ring);

    volatile uint64_t seed_point = 42;
    const auto point = blst_avx_test::make_random_vroom_curve_g2_projective(fp2ring, ring, seed_point);
    const auto scalar_bytes = blst_avx_test::seeded_scalar_bytes_le(321);

    benchmark::DoNotOptimize(point);
    benchmark::DoNotOptimize(scalar_bytes);

    for (auto _ : state) {
        const auto result = gls.multiply(point, scalar_bytes, fp2ring, g2, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("G2_ScalarMult_GLS");
}

BENCHMARK(BM_G1_ScalarMult_GLV);
BENCHMARK(BM_G2_ScalarMult_GLS);

}  // namespace
