// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

#include <benchmark/benchmark.h>

#include "blst/avx/cpu/precompute/gmp_wrapper.hpp"
#include "blst/avx/fields/fp2.h"
#include "blst/avx/fields/matrixnok_fp2.h"
#include "blst/avx/rns/bounded_ring.h"
#include "blst/avx/rns/matrixnok_fp.h"

namespace
{
    using blst::avx::BigInt;
    using Incumbent52Ring = blst::avx::BoundedRing<>;
    using Incumbent50Ring = blst::avx::BoundedRing<381, 8, 50>;
    using Incumbent50Fp2Ring = blst::avx::FP2Ring<Incumbent50Ring>;
    using MatrixNoKRing = blst::avx::matrixnok::FpRingCandidate;
    using MatrixNoKFp2Ring = blst::avx::matrixnok::Fp2RingCandidate;

    constexpr const char* left_decimal =
        "10688470598384192112760140533683975028842723510450044479374037626765189060660250327871479868795648192029502413";
    constexpr const char* right_decimal =
        "6303071894749576408239319750618994437764284500129504072810508811036909324432781187590536039296235607926728003";
    constexpr const char* third_decimal =
        "21376941196768384225520281067367950057685447020900088958748075253530378121320500655742959737591296384059004826";
    constexpr const char* fourth_decimal =
        "18909215684248729224717959251856983313292853500388512218431526433110727973298343562771608117888706823780184009";

    template<class Ring>
    void incumbent_modmul(benchmark::State& state)
    {
        const Ring ring;
        auto left = ring.from_bigint(BigInt(left_decimal, 10));
        auto right = ring.from_bigint(BigInt(right_decimal, 10));

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = ring.modmul(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_CatBLST_52_ModMul(benchmark::State& state)
    {
        incumbent_modmul<Incumbent52Ring>(state);
    }

    void BM_CatBLST_50_ModMul(benchmark::State& state)
    {
        incumbent_modmul<Incumbent50Ring>(state);
    }

    void BM_CatBLST_MatrixNoK_ModMul(benchmark::State& state)
    {
        const MatrixNoKRing ring;
        auto left = ring.from_canonical_50bit_digits({
            0x123456789abcdULL,
            0x2bcdef0123456ULL,
            0x3456789abcdefULL,
            0x0fedcba987654ULL,
            0x23456789abcdeULL,
            0x3abcdef012345ULL,
            0x156789abcdef0ULL,
            0x0000000001234ULL,
        });
        auto right = ring.from_canonical_50bit_digits({
            0x3edcba9876543ULL,
            0x176543210fedcULL,
            0x2ba9876543210ULL,
            0x3456789abcdefULL,
            0x1cdef01234567ULL,
            0x289abcdef0123ULL,
            0x0f0123456789aULL,
            0x0000000000abcULL,
        });
        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = ring.multiply(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_CatBLST_50_FP2_Mul(benchmark::State& state)
    {
        const Incumbent50Ring ring;
        const Incumbent50Fp2Ring fp2(ring);
        auto left = Incumbent50Fp2Ring::StandardElement(
            ring.from_bigint(BigInt(left_decimal, 10)),
            ring.from_bigint(BigInt(right_decimal, 10)));
        auto right = Incumbent50Fp2Ring::StandardElement(
            ring.from_bigint(BigInt(third_decimal, 10)),
            ring.from_bigint(BigInt(fourth_decimal, 10)));

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = fp2.mul(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_CatBLST_50_FP2_Square(benchmark::State& state)
    {
        const Incumbent50Ring ring;
        const Incumbent50Fp2Ring fp2(ring);
        auto value = Incumbent50Fp2Ring::StandardElement(
            ring.from_bigint(BigInt(left_decimal, 10)),
            ring.from_bigint(BigInt(right_decimal, 10)));

        benchmark::DoNotOptimize(value);
        for (auto _ : state) {
            auto result = fp2.square(value);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_CatBLST_MatrixNoK_FP2_Mul(benchmark::State& state)
    {
        const MatrixNoKRing ring;
        const MatrixNoKFp2Ring fp2(ring);
        auto left = MatrixNoKFp2Ring::Element{
            ring.from_canonical_50bit_digits({
                0x123456789abcdULL, 0x2bcdef0123456ULL,
                0x3456789abcdefULL, 0x0fedcba987654ULL,
                0x23456789abcdeULL, 0x3abcdef012345ULL,
                0x156789abcdef0ULL, 0x0000000001234ULL,
            }),
            ring.from_canonical_50bit_digits({
                0x3edcba9876543ULL, 0x176543210fedcULL,
                0x2ba9876543210ULL, 0x3456789abcdefULL,
                0x1cdef01234567ULL, 0x289abcdef0123ULL,
                0x0f0123456789aULL, 0x0000000000abcULL,
            }),
        };
        auto right = MatrixNoKFp2Ring::Element{
            ring.from_canonical_50bit_digits({
                0x2468acf13579aULL, 0x179bde02468acULL,
                0x28acf13579bdfULL, 0x1fdb97530eca9ULL,
                0x068acf13579bcULL, 0x3579bde02468bULL,
                0x2acf13579bde1ULL, 0x0000000002468ULL,
            }),
            ring.from_canonical_50bit_digits({
                0x3c962fc962fc9ULL, 0x062fc9632fc96ULL,
                0x02fc962fc9631ULL, 0x1d0369d0369cfULL,
                0x169cd0369d037ULL, 0x39d0369cd036aULL,
                0x2d0369d0369cfULL, 0x0000000002034ULL,
            }),
        };

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = fp2.multiply(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_CatBLST_MatrixNoK_FP2_Square(benchmark::State& state)
    {
        const MatrixNoKRing ring;
        const MatrixNoKFp2Ring fp2(ring);
        auto value = MatrixNoKFp2Ring::Element{
            ring.from_canonical_50bit_digits({
                0x123456789abcdULL, 0x2bcdef0123456ULL,
                0x3456789abcdefULL, 0x0fedcba987654ULL,
                0x23456789abcdeULL, 0x3abcdef012345ULL,
                0x156789abcdef0ULL, 0x0000000001234ULL,
            }),
            ring.from_canonical_50bit_digits({
                0x3edcba9876543ULL, 0x176543210fedcULL,
                0x2ba9876543210ULL, 0x3456789abcdefULL,
                0x1cdef01234567ULL, 0x289abcdef0123ULL,
                0x0f0123456789aULL, 0x0000000000abcULL,
            }),
        };

        benchmark::DoNotOptimize(value);
        for (auto _ : state) {
            auto result = fp2.square(value);
            benchmark::DoNotOptimize(result);
        }
    }

    template<std::size_t batch>
    void BM_CatBLST_50_BatchModMul(benchmark::State& state)
    {
        const Incumbent50Ring ring;
        std::array<Incumbent50Ring::StandardElement, batch> left;
        std::array<Incumbent50Ring::StandardElement, batch> right;
        for (std::size_t index = 0; index < batch; ++index) {
            left[index] = ring.from_bigint(BigInt(left_decimal, 10));
            right[index] = ring.from_bigint(BigInt(right_decimal, 10));
        }

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = ring.template batch_modmul<batch>(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

    template<std::size_t batch>
    void BM_CatBLST_MatrixNoK_BatchModMul(benchmark::State& state)
    {
        const MatrixNoKRing ring;
        const auto left_value = ring.from_canonical_50bit_digits({
            0x123456789abcdULL, 0x2bcdef0123456ULL,
            0x3456789abcdefULL, 0x0fedcba987654ULL,
            0x23456789abcdeULL, 0x3abcdef012345ULL,
            0x156789abcdef0ULL, 0x0000000001234ULL,
        });
        const auto right_value = ring.from_canonical_50bit_digits({
            0x3edcba9876543ULL, 0x176543210fedcULL,
            0x2ba9876543210ULL, 0x3456789abcdefULL,
            0x1cdef01234567ULL, 0x289abcdef0123ULL,
            0x0f0123456789aULL, 0x0000000000abcULL,
        });
        std::array<MatrixNoKRing::Element, batch> left;
        std::array<MatrixNoKRing::Element, batch> right;
        left.fill(left_value);
        right.fill(right_value);

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            auto result = ring.template multiply_batch<batch>(left, right);
            benchmark::DoNotOptimize(result);
        }
    }

#define BLST_MATRIXNOK_BATCH_BENCH(BATCH)                                              \
    BENCHMARK_TEMPLATE(BM_CatBLST_50_BatchModMul, BATCH);                             \
    BENCHMARK_TEMPLATE(BM_CatBLST_MatrixNoK_BatchModMul, BATCH)

    BENCHMARK(BM_CatBLST_52_ModMul);
    BENCHMARK(BM_CatBLST_50_ModMul);
    BENCHMARK(BM_CatBLST_MatrixNoK_ModMul);
    BENCHMARK(BM_CatBLST_50_FP2_Mul);
    BENCHMARK(BM_CatBLST_50_FP2_Square);
    BENCHMARK(BM_CatBLST_MatrixNoK_FP2_Mul);
    BENCHMARK(BM_CatBLST_MatrixNoK_FP2_Square);
    BLST_MATRIXNOK_BATCH_BENCH(1);
    BLST_MATRIXNOK_BATCH_BENCH(2);
    BLST_MATRIXNOK_BATCH_BENCH(4);
    BLST_MATRIXNOK_BATCH_BENCH(8);
    BLST_MATRIXNOK_BATCH_BENCH(12);

#undef BLST_MATRIXNOK_BATCH_BENCH
}
