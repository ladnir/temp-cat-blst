#include <benchmark/benchmark.h>

#include "blst/avx/tests/avx_test_utils.h"

#include <array>

namespace {

using blst_avx_test::AvxFp2Ring;
using blst_avx_test::AvxRing;

AvxFp2Ring::StandardElement random_fp2(gmp_randclass& rng, const BigInt& modulus, const AvxRing& ring)
{
    return blst_avx_test::make_vroom_fp2(
        ring,
        BigInt(rng.get_z_range(modulus.get_mpz())),
        BigInt(rng.get_z_range(modulus.get_mpz())));
}

static void BM_ModMul(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    const BigInt left_value(rng.get_z_range(modulus.get_mpz()));
    const BigInt right_value(rng.get_z_range(modulus.get_mpz()));

    const auto left = ring.from_bigint(left_value);
    const auto right = ring.from_bigint(right_value);

    benchmark::DoNotOptimize(left);
    benchmark::DoNotOptimize(right);

    for (auto _ : state) {
        const auto result = ring.modmul(left, right);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("ModMul");
}

static void BM_BatchModMul_8(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    std::array<AvxRing::StandardElement, 8> left;
    std::array<AvxRing::StandardElement, 8> right;

    for (int i = 0; i < 8; ++i) {
        left[i] = ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())));
        right[i] = ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())));
    }

    benchmark::DoNotOptimize(left);
    benchmark::DoNotOptimize(right);

    for (auto _ : state) {
        const auto result = ring.batch_modmul<8>(left, right);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("BatchModMul_8");
}

static void BM_FP2_Mul(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp2Ring fp2ring(ring);

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(52);
    const auto left = random_fp2(rng, modulus, ring);
    const auto right = random_fp2(rng, modulus, ring);

    benchmark::DoNotOptimize(left);
    benchmark::DoNotOptimize(right);

    for (auto _ : state) {
        const auto result = fp2ring.mul(left, right);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP2_Mul");
}

static void BM_FP2_Sqr(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp2Ring fp2ring(ring);

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(53);
    const auto value = random_fp2(rng, modulus, ring);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp2ring.square(value);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP2_Sqr");
}

BENCHMARK(BM_ModMul);
BENCHMARK(BM_BatchModMul_8);
BENCHMARK(BM_FP2_Mul);
BENCHMARK(BM_FP2_Sqr);

}  // namespace
