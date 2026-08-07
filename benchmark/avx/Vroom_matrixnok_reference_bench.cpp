// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

// This benchmark-only translation unit compiles the selected implementation
// directly from a pinned VROOM source tree. Keeping it separate from CatBLST's
// implementation avoids header and macro collisions while putting both
// implementations in the same executable and benchmark process.

#include <benchmark/benchmark.h>

#include "bls12_381_ring_params.hpp"
#include "fp2.hpp"

#include <array>
#include <cstddef>

namespace
{
    using VroomRing = bls12_381_params::BlsRingMatrixNoK;
    using VroomFp2Ring = FP2Ring<VroomRing>;

    constexpr const char* field_modulus_decimal =
        "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";
    constexpr const char* left_decimal =
        "10688470598384192112760140533683975028842723510450044479374037626765189060660250327871479868795648192029502413";
    constexpr const char* right_decimal =
        "6303071894749576408239319750618994437764284500129504072810508811036909324432781187590536039296235607926728003";
    constexpr const char* third_decimal =
        "21376941196768384225520281067367950057685447020900088958748075253530378121320500655742959737591296384059004826";
    constexpr const char* fourth_decimal =
        "18909215684248729224717959251856983313292853500388512218431526433110727973298343562771608117888706823780184009";

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

    void BM_VROOM_MatrixNoK_FP2_Mul(benchmark::State& state)
    {
        const VroomRing ring(BigInt(field_modulus_decimal, 10));
        const VroomFp2Ring fp2(ring);
        auto left = VroomFp2Ring::StandardElement(
            ring.from_bigint(BigInt(left_decimal, 10)),
            ring.from_bigint(BigInt(right_decimal, 10)));
        auto right = VroomFp2Ring::StandardElement(
            ring.from_bigint(BigInt(third_decimal, 10)),
            ring.from_bigint(BigInt(fourth_decimal, 10)));

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = fp2.mul(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_VROOM_MatrixNoK_FP2_Square(benchmark::State& state)
    {
        const VroomRing ring(BigInt(field_modulus_decimal, 10));
        const VroomFp2Ring fp2(ring);
        auto value = VroomFp2Ring::StandardElement(
            ring.from_bigint(BigInt(left_decimal, 10)),
            ring.from_bigint(BigInt(right_decimal, 10)));

        benchmark::DoNotOptimize(value);
        for (auto _ : state) {
            auto result = fp2.square(value);
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
    BENCHMARK(BM_VROOM_MatrixNoK_FP2_Mul);
    BENCHMARK(BM_VROOM_MatrixNoK_FP2_Square);
    BLST_VROOM_MATRIXNOK_BATCH_BENCH(1);
    BLST_VROOM_MATRIXNOK_BATCH_BENCH(2);
    BLST_VROOM_MATRIXNOK_BATCH_BENCH(4);
    BLST_VROOM_MATRIXNOK_BATCH_BENCH(8);
    BLST_VROOM_MATRIXNOK_BATCH_BENCH(12);

#undef BLST_VROOM_MATRIXNOK_BATCH_BENCH
}
