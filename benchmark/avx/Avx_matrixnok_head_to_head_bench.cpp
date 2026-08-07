// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

#include <benchmark/benchmark.h>

#include "blst/avx/cpu/precompute/gmp_wrapper.hpp"
#include "blst/avx/rns/bounded_ring.h"
#include "blst/avx/rns/matrixnok_fp.h"

namespace
{
    using blst::avx::BigInt;
    using Incumbent52Ring = blst::avx::BoundedRing<>;
    using Incumbent50Ring = blst::avx::BoundedRing<381, 8, 50>;
    using MatrixNoKRing = blst::avx::matrixnok::FpRingCandidate;

    constexpr const char* left_decimal =
        "10688470598384192112760140533683975028842723510450044479374037626765189060660250327871479868795648192029502413";
    constexpr const char* right_decimal =
        "6303071894749576408239319750618994437764284500129504072810508811036909324432781187590536039296235607926728003";

    template<class Ring>
    void incumbent_modmul(benchmark::State& state)
    {
        const Ring ring;
        const auto left = ring.from_bigint(BigInt(left_decimal, 10));
        const auto right = ring.from_bigint(BigInt(right_decimal, 10));

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            const auto result = ring.modmul(left, right);
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
        const auto left = ring.from_canonical_50bit_digits({
            0x123456789abcdULL,
            0x2bcdef0123456ULL,
            0x3456789abcdefULL,
            0x0fedcba987654ULL,
            0x23456789abcdeULL,
            0x3abcdef012345ULL,
            0x156789abcdef0ULL,
            0x0000000001234ULL,
        });
        const auto right = ring.from_canonical_50bit_digits({
            0x3edcba9876543ULL,
            0x176543210fedcULL,
            0x2ba9876543210ULL,
            0x3456789abcdefULL,
            0x1cdef01234567ULL,
            0x289abcdef0123ULL,
            0x0f0123456789aULL,
            0x0000000000abcULL,
        });
        MatrixNoKRing::Element result;

        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        for (auto _ : state) {
            ring.multiply(result, left, right);
            benchmark::DoNotOptimize(result);
        }
    }

    BENCHMARK(BM_CatBLST_52_ModMul);
    BENCHMARK(BM_CatBLST_50_ModMul);
    BENCHMARK(BM_CatBLST_MatrixNoK_ModMul);
}
