#include <benchmark/benchmark.h>

#include "test/avx/avx_test_utils.h"

#include <array>

using namespace blst::avx;

namespace {

using blst_avx_test::AvxFp12Ring;
using blst_avx_test::AvxRing;

AvxFp12Ring::StandardElement generate_random_fp12(const AvxRing& ring, uint64_t seed)
{
    return blst_avx_test::make_avx_fp12_from_seed(ring, seed);
}

static void BM_FP12_Mul(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp12Ring fp12ring(ring);
    volatile uint64_t seed_left = 42;
    volatile uint64_t seed_right = 123;
    const auto left = generate_random_fp12(ring, seed_left);
    const auto right = generate_random_fp12(ring, seed_right);

    benchmark::DoNotOptimize(left);
    benchmark::DoNotOptimize(right);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_mul(left, right, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_Mul");
}

static void BM_FP12_Sqr(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp12Ring fp12ring(ring);
    volatile uint64_t seed = 42;
    const auto value = generate_random_fp12(ring, seed);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp12ring.fp12_flat_sqr(value, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_Sqr");
}

BENCHMARK(BM_FP12_Mul);
BENCHMARK(BM_FP12_Sqr);

}  // namespace
