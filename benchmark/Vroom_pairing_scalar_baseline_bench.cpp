#include <benchmark/benchmark.h>

#include "blst/avx/tests/avx_test_utils.h"

namespace {

static void BM_BlstScalarPairingFinalExp(benchmark::State& state)
{
    const blst::G1AffinePoint g1 = blst_vroom_test::make_scalar_g1_generator_affine();
    const blst::G2AffinePoint g2 = blst_vroom_test::make_scalar_g2_generator_affine();

    for (auto _ : state) {
        blst::GT result = {};
        blst::GT::miller_loop(result, &g2, &g1);
        blst::GT::final_mult(result, result);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("BlstScalarPairingFinalExp");
}

BENCHMARK(BM_BlstScalarPairingFinalExp);

}  // namespace
