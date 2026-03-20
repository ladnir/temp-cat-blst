#include <benchmark/benchmark.h>

#include "blst/avx/tests/avx_test_utils.h"

#include <array>

namespace {

using blst_avx_test::Avx50Fp2Ring;
using blst_avx_test::Avx50Ring;

Avx50Fp2Ring::StandardElement random_fp2_50bit(gmp_randclass& rng, const BigInt& modulus, const Avx50Ring& ring)
{
    return blst_avx_test::make_avx_fp2(
        ring,
        BigInt(rng.get_z_range(modulus.get_mpz())),
        BigInt(rng.get_z_range(modulus.get_mpz())));
}

template<int Batch>
std::array<Avx50Ring::StandardElement, Batch> make_modmul_inputs_50bit(const Avx50Ring& ring, gmp_randclass& rng, const BigInt& modulus)
{
    std::array<Avx50Ring::StandardElement, Batch> values = {};
    for (int i = 0; i < Batch; ++i) {
        values[i] = ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())));
    }
    return values;
}

static void BM_ModMul_50bit(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const Avx50Ring ring(modulus);

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

    state.SetLabel("ModMul_50bit");
}

#define BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50(BATCH)                                \
    static void BM_BatchModMul_##BATCH##_50bit(benchmark::State& state)               \
    {                                                                                  \
        const BigInt modulus = blst_avx_test::modulus_bigint();                     \
        const Avx50Ring ring(modulus);                                               \
        gmp_randclass rng(gmp_randinit_default);                                       \
        rng.seed(42);                                                                  \
        const auto left = make_modmul_inputs_50bit<BATCH>(ring, rng, modulus);        \
        const auto right = make_modmul_inputs_50bit<BATCH>(ring, rng, modulus);       \
                                                                                       \
        benchmark::DoNotOptimize(left);                                                \
        benchmark::DoNotOptimize(right);                                               \
                                                                                       \
        for (auto _ : state) {                                                         \
            const auto result = ring.batch_modmul<BATCH>(left, right);                \
            benchmark::DoNotOptimize(result);                                          \
        }                                                                              \
                                                                                       \
        state.SetLabel("BatchModMul_" #BATCH "_50bit");                                \
    }                                                                                  \
    BENCHMARK(BM_BatchModMul_##BATCH##_50bit)

BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50(1);
BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50(2);
BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50(3);
BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50(4);
BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50(5);
BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50(6);
BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50(7);
BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50(8);
BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50(10);
BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50(12);

#undef BLST_VROOM_DEFINE_BATCH_MODMUL_BENCH_50

static void BM_FP2_Mul_50bit(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const Avx50Ring ring(modulus);
    const Avx50Fp2Ring fp2ring(ring);

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(52);
    const auto left = random_fp2_50bit(rng, modulus, ring);
    const auto right = random_fp2_50bit(rng, modulus, ring);

    benchmark::DoNotOptimize(left);
    benchmark::DoNotOptimize(right);

    for (auto _ : state) {
        const auto result = fp2ring.mul(left, right);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP2_Mul_50bit");
}

static void BM_FP2_Sqr_50bit(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const Avx50Ring ring(modulus);
    const Avx50Fp2Ring fp2ring(ring);

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(53);
    const auto value = random_fp2_50bit(rng, modulus, ring);

    benchmark::DoNotOptimize(value);

    for (auto _ : state) {
        const auto result = fp2ring.square(value);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP2_Sqr_50bit");
}

BENCHMARK(BM_ModMul_50bit);
BENCHMARK(BM_FP2_Mul_50bit);
BENCHMARK(BM_FP2_Sqr_50bit);

}  // namespace
