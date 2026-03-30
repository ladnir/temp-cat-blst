#include <benchmark/benchmark.h>

#include "blst/avx/rns/inversion.h"
#include "test/avx/avx_test_utils.h"

using namespace blst::avx;

namespace {

using blst_avx_test::AvxRing;

static void BM_Inverse(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const BLS381AddChainInversion<AvxRing> inverter(ring);

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

    state.SetLabel("Inverse");
}

BENCHMARK(BM_Inverse);

}  // namespace
