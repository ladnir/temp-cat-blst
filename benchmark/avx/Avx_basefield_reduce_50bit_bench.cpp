#include <benchmark/benchmark.h>

#include "test/avx/avx_test_utils.h"

#include <array>

using namespace blst::avx;

namespace {

using blst_avx_test::Avx50Ring;

template<int Batch>
std::array<Avx50Ring::ReadyElement<Avx50Ring::MAX_ADD>, Batch> make_ready_elements_50bit(const Avx50Ring& ring)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);

    std::array<Avx50Ring::ReadyElement<Avx50Ring::MAX_ADD>, Batch> ready_elements = {};
    for (int i = 0; i < Batch; ++i) {
        const auto left = ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())));
        const auto right = ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())));
        ready_elements[i] = ring.template ready<Avx50Ring::MAX_ADD>(left * right);
    }
    return ready_elements;
}

#define BLST_AVX_DEFINE_BATCH_REDUCE_BENCH_50(BATCH)                                  \
    static void BM_BatchReduceExpand_##BATCH##_50bit(benchmark::State& state)           \
    {                                                                                    \
        const Avx50Ring ring(blst_avx_test::modulus_bigint());                      \
        const auto ready_elements = make_ready_elements_50bit<BATCH>(ring);             \
                                                                                         \
        benchmark::DoNotOptimize(ready_elements);                                        \
                                                                                         \
        for (auto _ : state) {                                                           \
            const auto result = ring.batch_reduce_expand<BATCH>(ready_elements);         \
            benchmark::DoNotOptimize(result);                                            \
        }                                                                                \
                                                                                         \
        state.SetLabel("BatchReduceExpand_" #BATCH "_50bit");                            \
    }                                                                                    \
    BENCHMARK(BM_BatchReduceExpand_##BATCH##_50bit)

BLST_AVX_DEFINE_BATCH_REDUCE_BENCH_50(1);
BLST_AVX_DEFINE_BATCH_REDUCE_BENCH_50(2);
BLST_AVX_DEFINE_BATCH_REDUCE_BENCH_50(3);
BLST_AVX_DEFINE_BATCH_REDUCE_BENCH_50(4);
BLST_AVX_DEFINE_BATCH_REDUCE_BENCH_50(5);
BLST_AVX_DEFINE_BATCH_REDUCE_BENCH_50(6);
BLST_AVX_DEFINE_BATCH_REDUCE_BENCH_50(12);

#undef BLST_AVX_DEFINE_BATCH_REDUCE_BENCH_50

}  // namespace
