// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

#include <benchmark/benchmark.h>

#include "blst/avx/curve/matrixnok_g2.h"
#include "blst/avx/curve/point.h"

#include <array>
#include <cstdint>

namespace
{
    namespace matrixnok = blst::avx::matrixnok;
    namespace incumbent = blst::avx::avx_curve_detail;

    constexpr std::array<uint8_t, 32> dense_scalar = {
        0xa5, 0x5a, 0xc3, 0x3c, 0x96, 0x69, 0xf0, 0x0f,
        0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0,
        0x0f, 0xed, 0xcb, 0xa9, 0x87, 0x65, 0x43, 0x21,
        0x10, 0x32, 0x54, 0x76, 0x08, 0x07, 0x06, 0x3f,
    };

    matrixnok::G2Fp2RingCandidate::StandardElement candidate_fp2(
        const matrixnok::FpRingCandidate& ring, const blst::Fp2& value)
    {
        return {ring.from_scalar(value[0]), ring.from_scalar(value[1])};
    }

    matrixnok::G2ProjectiveCandidate candidate_generator(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::G2Fp2RingCandidate& fp2ring)
    {
        const auto affine = blst::scalar::G2AffinePoint::generator();
        return {
            candidate_fp2(ring, affine.x),
            candidate_fp2(ring, affine.y),
            fp2ring.one(),
        };
    }

    incumbent::CurveG2Projective incumbent_generator()
    {
        const auto& generator = blst::avx::G2Point::generator();
        return {generator.X, generator.Y, generator.Z};
    }

    void BM_CatBLST_Incumbent_G2Add(benchmark::State& state)
    {
        const auto& curve = incumbent::curve_g2();
        const auto& fp2ring = incumbent::curve_fp2_ring();
        const auto generator = incumbent_generator();
        auto current = generator;
        for (auto _ : state) {
            current = curve.add_point(current, generator, fp2ring);
            benchmark::DoNotOptimize(current);
        }
    }

    void BM_CatBLST_MatrixNoKCandidate_G2Add(benchmark::State& state)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        const matrixnok::G2CurveCandidate curve;
        const auto generator = candidate_generator(ring, fp2ring);
        auto current = generator;
        for (auto _ : state) {
            current = curve.add_point(current, generator, fp2ring);
            benchmark::DoNotOptimize(current);
        }
    }

    void BM_CatBLST_Incumbent_G2Double(benchmark::State& state)
    {
        const auto& curve = incumbent::curve_g2();
        const auto& fp2ring = incumbent::curve_fp2_ring();
        auto current = incumbent_generator();
        for (auto _ : state) {
            current = curve.double_point(current, fp2ring);
            benchmark::DoNotOptimize(current);
        }
    }

    void BM_CatBLST_MatrixNoKCandidate_G2Double(benchmark::State& state)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        const matrixnok::G2CurveCandidate curve;
        auto current = candidate_generator(ring, fp2ring);
        for (auto _ : state) {
            current = curve.double_point(current, fp2ring);
            benchmark::DoNotOptimize(current);
        }
    }

    void BM_CatBLST_Incumbent_G2GlsCt(benchmark::State& state)
    {
        const auto& ring = incumbent::curve_ring();
        const auto& fp2ring = incumbent::curve_fp2_ring();
        const auto& curve = incumbent::curve_g2();
        const auto& gls = incumbent::curve_g2_gls();
        auto current = incumbent_generator();
        for (auto _ : state) {
            current = gls.multiply(current, dense_scalar, fp2ring, curve, ring);
            benchmark::DoNotOptimize(current);
        }
    }

    void BM_CatBLST_MatrixNoKCandidate_G2GlsCt(benchmark::State& state)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        const matrixnok::G2CurveCandidate curve;
        const matrixnok::G2GLSCandidate gls(ring);
        auto current = candidate_generator(ring, fp2ring);
        for (auto _ : state) {
            current = gls.multiply(current, dense_scalar, fp2ring, curve, ring);
            benchmark::DoNotOptimize(current);
        }
    }

    void BM_CatBLST_Incumbent_G2GlsVt(benchmark::State& state)
    {
        const auto& ring = incumbent::curve_ring();
        const auto& fp2ring = incumbent::curve_fp2_ring();
        const auto& curve = incumbent::curve_g2();
        const auto& gls = incumbent::curve_g2_gls();
        auto current = incumbent_generator();
        for (auto _ : state) {
            current = gls.multiply_vartime(
                current, dense_scalar, fp2ring, curve, ring);
            benchmark::DoNotOptimize(current);
        }
    }

    void BM_CatBLST_MatrixNoKCandidate_G2GlsVt(benchmark::State& state)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        const matrixnok::G2CurveCandidate curve;
        const matrixnok::G2GLSCandidate gls(ring);
        auto current = candidate_generator(ring, fp2ring);
        for (auto _ : state) {
            current = gls.multiply_vartime(
                current, dense_scalar, fp2ring, curve, ring);
            benchmark::DoNotOptimize(current);
        }
    }

    BENCHMARK(BM_CatBLST_Incumbent_G2Add);
    BENCHMARK(BM_CatBLST_MatrixNoKCandidate_G2Add);
    BENCHMARK(BM_CatBLST_Incumbent_G2Double);
    BENCHMARK(BM_CatBLST_MatrixNoKCandidate_G2Double);
    BENCHMARK(BM_CatBLST_Incumbent_G2GlsCt);
    BENCHMARK(BM_CatBLST_MatrixNoKCandidate_G2GlsCt);
    BENCHMARK(BM_CatBLST_Incumbent_G2GlsVt);
    BENCHMARK(BM_CatBLST_MatrixNoKCandidate_G2GlsVt);
} // namespace
