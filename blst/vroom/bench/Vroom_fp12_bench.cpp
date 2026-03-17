#include <benchmark/benchmark.h>

#include "blst/vroom/tests/vroom_basefield_test_utils.hpp"

#include <array>

namespace {

using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomRing;

VroomFp12Ring::StandardElement generate_random_fp12(const VroomRing& ring, uint64_t seed)
{
    VroomFp12Ring::StandardElement result = {};
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    volatile uint64_t stable_seed = seed;
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(stable_seed);

    for (size_t i = 0; i < result.size(); ++i) {
        result[i] = ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())));
    }
    return result;
}

static void BM_FP12_Mul(benchmark::State& state)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const VroomRing ring(modulus);
    const VroomFp12Ring fp12ring(ring);
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
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const VroomRing ring(modulus);
    const VroomFp12Ring fp12ring(ring);
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
