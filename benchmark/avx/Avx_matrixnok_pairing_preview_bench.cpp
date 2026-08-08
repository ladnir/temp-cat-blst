// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

#include <benchmark/benchmark.h>

#include "blst/avx/curve/matrixnok_g1.h"
#include "blst/avx/curve/matrixnok_g2.h"
#include "blst/avx/curve/GT.h"
#include "blst/avx/pairing/matrixnok_pairing.h"
#include "blst/avx/pairing/pairing.h"

namespace
{
    namespace incumbent = blst::avx::avx_gt_detail;
    namespace matrixnok = blst::avx::matrixnok;

    matrixnok::G1AffineCandidate matrixnok_g1_generator(
        const matrixnok::FpRingCandidate& ring)
    {
        const auto& point = blst::scalar::G1AffinePoint::generator();
        return {ring.from_scalar(point.x), ring.from_scalar(point.y)};
    }

    matrixnok::G2AffineCandidate matrixnok_g2_generator(
        const matrixnok::FpRingCandidate& ring)
    {
        const auto& point = blst::scalar::G2AffinePoint::generator();
        return {
            {ring.from_scalar(point.x[0]), ring.from_scalar(point.x[1])},
            {ring.from_scalar(point.y[0]), ring.from_scalar(point.y[1])}};
    }

    void BM_CatBLST_Incumbent_Miller(benchmark::State& state)
    {
        const auto& ring = incumbent::ring();
        const auto& fp2ring = incumbent::fp2ring();
        const auto& fp12ring = incumbent::fp12ring();
        const auto& public_g1 = blst::avx::G1AffinePoint::generator();
        const auto& public_g2 = blst::avx::G2AffinePoint::generator();
        const blst::avx::AffinePoint<incumbent::Ring::StandardElement> g1{
            public_g1.X, public_g1.Y};
        const blst::avx::AffinePoint<incumbent::PairingFp2Element> g2{
            public_g2.X, public_g2.Y};
        const blst::avx::Miller<incumbent::Ring> miller(g1, g2, ring);
        for (auto _ : state) {
            auto result = miller.miller_loop(ring, fp2ring, fp12ring);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_CatBLST_MatrixNoKCandidate_Miller(benchmark::State& state)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        const matrixnok::Fp12RingCandidate fp12ring(ring);
        const auto g1 = matrixnok_g1_generator(ring);
        const auto g2 = matrixnok_g2_generator(ring);
        const matrixnok::MillerCandidate miller(g1, g2, ring);
        for (auto _ : state) {
            auto result = miller.miller_loop(ring, fp2ring, fp12ring);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_CatBLST_Incumbent_Pairing(benchmark::State& state)
    {
        const auto& ring = incumbent::ring();
        const auto& fp2ring = incumbent::fp2ring();
        const auto& fp12ring = incumbent::fp12ring();
        const auto& public_g1 = blst::avx::G1AffinePoint::generator();
        const auto& public_g2 = blst::avx::G2AffinePoint::generator();
        const blst::avx::AffinePoint<incumbent::Ring::StandardElement> g1{
            public_g1.X, public_g1.Y};
        const blst::avx::AffinePoint<incumbent::PairingFp2Element> g2{
            public_g2.X, public_g2.Y};
        const blst::avx::FP12RNS_BLSTInverter<incumbent::Ring> inverter{};
        for (auto _ : state) {
            auto result = blst::avx::pairing(
                g1, g2, inverter, fp2ring, fp12ring, ring);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_CatBLST_MatrixNoKCandidate_Pairing(benchmark::State& state)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        const matrixnok::Fp12RingCandidate fp12ring(ring);
        const auto g1 = matrixnok_g1_generator(ring);
        const auto g2 = matrixnok_g2_generator(ring);
        const matrixnok::PairingInverterCandidate inverter{};
        for (auto _ : state) {
            auto result = blst::avx::pairing(
                g1, g2, inverter, fp2ring, fp12ring, ring);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_CatBLST_MatrixNoKCandidate_Fp12Inverse(
        benchmark::State& state)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        const matrixnok::Fp12RingCandidate fp12ring(ring);
        const auto g1 = matrixnok_g1_generator(ring);
        const auto g2 = matrixnok_g2_generator(ring);
        const matrixnok::MillerCandidate miller(g1, g2, ring);
        const auto value = miller.miller_loop(ring, fp2ring, fp12ring);
        const matrixnok::PairingInverterCandidate inverter{};
        for (auto _ : state) {
            auto result = inverter.invert(value, fp12ring, ring);
            benchmark::DoNotOptimize(result);
        }
    }

    void BM_CatBLST_MatrixNoKCandidate_FinalExponentAfterInverse(
        benchmark::State& state)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        const matrixnok::Fp12RingCandidate fp12ring(ring);
        const auto g1 = matrixnok_g1_generator(ring);
        const auto g2 = matrixnok_g2_generator(ring);
        const matrixnok::MillerCandidate miller(g1, g2, ring);
        const auto value = miller.miller_loop(ring, fp2ring, fp12ring);
        const matrixnok::PairingInverterCandidate inverter{};
        const auto inverse = inverter.invert(value, fp12ring, ring);
        for (auto _ : state) {
            auto result = blst::avx::final_exp_after_inverse(
                value, inverse, fp12ring, ring);
            benchmark::DoNotOptimize(result);
        }
    }

    BENCHMARK(BM_CatBLST_Incumbent_Miller);
    BENCHMARK(BM_CatBLST_MatrixNoKCandidate_Miller);
    BENCHMARK(BM_CatBLST_Incumbent_Pairing);
    BENCHMARK(BM_CatBLST_MatrixNoKCandidate_Pairing);
    BENCHMARK(BM_CatBLST_MatrixNoKCandidate_Fp12Inverse);
    BENCHMARK(BM_CatBLST_MatrixNoKCandidate_FinalExponentAfterInverse);
} // namespace
