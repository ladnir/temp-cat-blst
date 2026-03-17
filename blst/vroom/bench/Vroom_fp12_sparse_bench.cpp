#include <benchmark/benchmark.h>

#include "blst/vroom/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomRing;

static void BM_FP12_MulByXy0z00(benchmark::State& state)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const VroomRing ring(modulus);
    const VroomFp12Ring fp12ring(ring);
    volatile uint64_t x_seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, x_seed);

    const auto sparse_coeffs = blst_vroom_test::make_vroom_sparse_fp6_from_seed(ring, 123);
    const auto c_term = blst_vroom_test::make_vroom_fp2_from_sparse6(sparse_coeffs, 0);
    const auto x_term = blst_vroom_test::make_vroom_fp2_from_sparse6(sparse_coeffs, 1);
    const auto y_term = blst_vroom_test::make_vroom_fp2_from_sparse6(sparse_coeffs, 2);

    benchmark::DoNotOptimize(value);
    benchmark::DoNotOptimize(sparse_coeffs);

    for (auto _ : state) {
        VroomFp12Ring::StandardElement result = {};
        fp12ring.flat_mulxy00z0(result, value, c_term, x_term, y_term, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_MulByXy0z00");
}

BENCHMARK(BM_FP12_MulByXy0z00);

}  // namespace
