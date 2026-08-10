// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

#include <benchmark/benchmark.h>

#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/matrixnok_g1.h"
#include "blst/avx/curve/matrixnok_g2.h"
#include "blst/avx/pairing/matrixnok_pairing.h"
#include "blst/avx/pairing/pairing.h"
#include "blst/scalar/curve/GT.h"

namespace {

static void BM_ScalarPairing_MillerLoop(benchmark::State& state)
{
    const auto& p = blst::scalar::G1AffinePoint::generator();
    const auto& q = blst::scalar::G2AffinePoint::generator();

    for (auto _ : state) {
        blst::scalar::MillerResult out = {};
        blst::scalar::GT::miller_loop(out, &q, &p);
        benchmark::DoNotOptimize(out);
    }
}

static void BM_ScalarPairing_Full(benchmark::State& state)
{
    const auto& p = blst::scalar::G1AffinePoint::generator();
    const auto& q = blst::scalar::G2AffinePoint::generator();

    for (auto _ : state) {
        blst::scalar::MillerResult miller = {};
        blst::scalar::GT out = {};
        blst::scalar::GT::miller_loop(miller, &q, &p);
        blst::scalar::GT::final_exp(out, miller);
        benchmark::DoNotOptimize(out);
    }
}

static void BM_AvxPairing_MillerLoop(benchmark::State& state)
{
    const auto& p = blst::avx::G1AffinePoint::generator();
    const auto& q = blst::avx::G2AffinePoint::generator();

    for (auto _ : state) {
        blst::avx::MillerResult out = {};
        blst::avx::GT::miller_loop(out, &q, &p);
        benchmark::DoNotOptimize(out);
    }
}

static void BM_AvxPairing_Full(benchmark::State& state)
{
    const auto& p = blst::avx::G1AffinePoint::generator();
    const auto& q = blst::avx::G2AffinePoint::generator();

    for (auto _ : state) {
        blst::avx::MillerResult miller = {};
        blst::avx::GT out = {};
        blst::avx::GT::miller_loop(miller, &q, &p);
        blst::avx::GT::final_exp(out, miller);
        benchmark::DoNotOptimize(out);
    }
}

static void BM_DirectMatrixNoKPairing_MillerLoop(benchmark::State& state)
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::G2Fp2RingCandidate fp2ring(ring);
    const blst::avx::matrixnok::Fp12RingCandidate fp12ring(ring);
    const auto& public_p = blst::avx::G1AffinePoint::generator();
    const auto& public_q = blst::avx::G2AffinePoint::generator();
    const blst::avx::matrixnok::G1AffineCandidate p{public_p.X, public_p.Y};
    const blst::avx::matrixnok::G2AffineCandidate q{public_q.X, public_q.Y};
    const blst::avx::matrixnok::MillerCandidate miller(p, q, ring);

    for (auto _ : state) {
        auto out = miller.miller_loop(ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(out);
    }
}

static void BM_DirectMatrixNoKPairing_Full(benchmark::State& state)
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::G2Fp2RingCandidate fp2ring(ring);
    const blst::avx::matrixnok::Fp12RingCandidate fp12ring(ring);
    const auto& public_p = blst::avx::G1AffinePoint::generator();
    const auto& public_q = blst::avx::G2AffinePoint::generator();
    const blst::avx::matrixnok::G1AffineCandidate p{public_p.X, public_p.Y};
    const blst::avx::matrixnok::G2AffineCandidate q{public_q.X, public_q.Y};
    const blst::avx::matrixnok::PairingInverterCandidate inverter;

    for (auto _ : state) {
        auto out = blst::avx::pairing(
            p, q, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(out);
    }
}

BENCHMARK(BM_ScalarPairing_MillerLoop);
BENCHMARK(BM_ScalarPairing_Full);
BENCHMARK(BM_AvxPairing_MillerLoop);
BENCHMARK(BM_AvxPairing_Full);
BENCHMARK(BM_DirectMatrixNoKPairing_MillerLoop);
BENCHMARK(BM_DirectMatrixNoKPairing_Full);

}  // namespace
