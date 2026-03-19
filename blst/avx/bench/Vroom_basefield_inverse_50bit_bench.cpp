#include <benchmark/benchmark.h>

#include "blst/avx/src/inversion.hpp"
#include "blst/avx/tests/vroom_basefield_test_utils.hpp"

namespace {

using blst_vroom_test::Vroom50Ring;

static void BM_Inverse_50bit(benchmark::State& state)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    const Vroom50Ring ring(modulus);
    const BLS381AddChainInversion<Vroom50Ring> inverter(ring);

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    BigInt random_value(rng.get_z_range(modulus.get_mpz()));
    if (random_value == BigInt(0)) {
        random_value = BigInt(1);
    }

    const auto element = ring.from_bigint(random_value);
    benchmark::DoNotOptimize(element);

    for (auto _ : state) {
        const auto result = inverter.invert(element, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Inverse_50bit");
}

BENCHMARK(BM_Inverse_50bit);

}  // namespace
