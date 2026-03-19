#include <benchmark/benchmark.h>

#include "blst/avx/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomRing;

static void BM_FP12_CyclotomicSqr(benchmark::State& state)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const VroomRing ring(modulus);
    const VroomFp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = blst_vroom_test::make_vroom_fp12_from_seed(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_cyclotomic_sqr(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_CyclotomicSqr");
}

BENCHMARK(BM_FP12_CyclotomicSqr);

}  // namespace
