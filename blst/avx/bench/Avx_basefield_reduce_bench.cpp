#include <benchmark/benchmark.h>

#include "blst/avx/tests/avx_test_utils.h"

#include <array>

namespace {

using blst_vroom_test::VroomRing;

template<int Batch>
std::array<VroomRing::ReadyElement<VroomRing::MAX_ADD>, Batch> make_ready_elements(const VroomRing& ring)
{
    const BigInt modulus = blst_vroom_test::modulus_bigint();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);

    std::array<VroomRing::ReadyElement<VroomRing::MAX_ADD>, Batch> ready_elements = {};
    for (int i = 0; i < Batch; ++i) {
        const auto left = ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())));
        const auto right = ring.from_bigint(BigInt(rng.get_z_range(modulus.get_mpz())));
        ready_elements[i] = ring.template ready<VroomRing::MAX_ADD>(left * right);
    }
    return ready_elements;
}

#define BLST_VROOM_DEFINE_BATCH_REDUCE_BENCH(BATCH)                                \
    static void BM_BatchReduceExpand_##BATCH(benchmark::State& state)              \
    {                                                                              \
        const VroomRing ring(blst_vroom_test::modulus_bigint());                   \
        const auto ready_elements = make_ready_elements<BATCH>(ring);              \
                                                                                   \
        benchmark::DoNotOptimize(ready_elements);                                  \
                                                                                   \
        for (auto _ : state) {                                                     \
            const auto result = ring.batch_reduce_expand<BATCH>(ready_elements);   \
            benchmark::DoNotOptimize(result);                                      \
        }                                                                          \
                                                                                   \
        state.SetLabel("BatchReduceExpand_" #BATCH);                               \
    }                                                                              \
    BENCHMARK(BM_BatchReduceExpand_##BATCH)

BLST_VROOM_DEFINE_BATCH_REDUCE_BENCH(1);
BLST_VROOM_DEFINE_BATCH_REDUCE_BENCH(2);
BLST_VROOM_DEFINE_BATCH_REDUCE_BENCH(3);
BLST_VROOM_DEFINE_BATCH_REDUCE_BENCH(4);
BLST_VROOM_DEFINE_BATCH_REDUCE_BENCH(5);
BLST_VROOM_DEFINE_BATCH_REDUCE_BENCH(6);
BLST_VROOM_DEFINE_BATCH_REDUCE_BENCH(12);

#undef BLST_VROOM_DEFINE_BATCH_REDUCE_BENCH

}  // namespace
