#include <benchmark/benchmark.h>

#define BLST_VROOM_SKIP_TEST_VALUES
#include "blst/vroom/src/scalar_mult.hpp"
#include "blst/vroom/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::VroomCurveFp2Ring;
using blst_vroom_test::VroomCurveG1;
using blst_vroom_test::VroomCurveG2;
using blst_vroom_test::VroomCurveRing;

static void BM_G1_ScalarMult_GLV(benchmark::State& state)
{
    const VroomCurveRing ring(blst_vroom_test::modulus_bigint());
    const VroomCurveG1 g1 = {};
    const GLVMult<VroomCurveRing> glv(g1, ring);

    volatile uint64_t seed_point = 42;
    const auto point = blst_vroom_test::make_random_vroom_curve_g1_projective(ring, seed_point);
    const auto scalar_bytes = blst_vroom_test::seeded_scalar_bytes_le(123);

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
    const VroomCurveRing ring(blst_vroom_test::modulus_bigint());
    const VroomCurveFp2Ring fp2ring(ring);
    const VroomCurveG2 g2 = {};
    const GLSMult<VroomCurveRing> gls(ring);

    volatile uint64_t seed_point = 42;
    const auto point = blst_vroom_test::make_random_vroom_curve_g2_projective(fp2ring, ring, seed_point);
    const auto scalar_bytes = blst_vroom_test::seeded_scalar_bytes_le(321);

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
