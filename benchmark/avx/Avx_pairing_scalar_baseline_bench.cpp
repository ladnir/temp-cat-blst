#include <benchmark/benchmark.h>

#include "test/avx/avx_test_utils.h"

using namespace blst::avx;

namespace {

static void BM_BlstScalarPairingFinalExp(benchmark::State& state)
{
    const blst::scalar::G1AffinePoint g1 = blst_avx_test::make_scalar_g1_generator_affine();
    const blst::scalar::G2AffinePoint g2 = blst_avx_test::make_scalar_g2_generator_affine();

    for (auto _ : state) {
        blst::scalar::GT result = {};
        blst::scalar::GT::miller_loop(result, &g2, &g1);
        blst::scalar::GT::final_mult(result, result);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("BlstScalarPairingFinalExp");
}

BENCHMARK(BM_BlstScalarPairingFinalExp);

}  // namespace
