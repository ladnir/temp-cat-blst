// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

// This benchmark-only translation unit compiles the selected implementation
// directly from a pinned VROOM source tree. Keeping it separate from CatBLST's
// implementation avoids header and macro collisions while putting both
// implementations in the same executable and benchmark process.

#include <benchmark/benchmark.h>

#include "bls12_381_ring_params.hpp"

#include <array>
#include <cstddef>

namespace
{
    using VroomRing = bls12_381_params::BlsRingMatrixNoK;

    constexpr const char* field_modulus_decimal =
        "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";
    constexpr const char* left_decimal =
        "10688470598384192112760140533683975028842723510450044479374037626765189060660250327871479868795648192029502413";
    constexpr const char* right_decimal =
        "6303071894749576408239319750618994437764284500129504072810508811036909324432781187590536039296235607926728003";

    void BM_VROOM_MatrixNoK_ModMul(benchmark::State& state)
    {
        const VroomRing ring(BigInt(field_modulus_decimal, 10));
        auto left = ring.from_bigint(BigInt(left_decimal, 10));
        auto right = ring.from_bigint(BigInt(right_decimal, 10));

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = ring.modmul(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

    template<std::size_t batch>
    void BM_VROOM_MatrixNoK_BatchModMul(benchmark::State& state)
    {
        const VroomRing ring(BigInt(field_modulus_decimal, 10));
        const auto left_value = ring.from_bigint(BigInt(left_decimal, 10));
        const auto right_value = ring.from_bigint(BigInt(right_decimal, 10));
        std::array<VroomRing::StandardElement, batch> left;
        std::array<VroomRing::StandardElement, batch> right;
        left.fill(left_value);
        right.fill(right_value);

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = ring.template batch_modmul<batch>(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

#define BLST_VROOM_MATRIXNOK_BATCH_BENCH(BATCH) \
    BENCHMARK_TEMPLATE(BM_VROOM_MatrixNoK_BatchModMul, BATCH)

    BENCHMARK(BM_VROOM_MatrixNoK_ModMul);
    BLST_VROOM_MATRIXNOK_BATCH_BENCH(1);
    BLST_VROOM_MATRIXNOK_BATCH_BENCH(2);
    BLST_VROOM_MATRIXNOK_BATCH_BENCH(4);
    BLST_VROOM_MATRIXNOK_BATCH_BENCH(8);
    BLST_VROOM_MATRIXNOK_BATCH_BENCH(12);

#undef BLST_VROOM_MATRIXNOK_BATCH_BENCH
}
