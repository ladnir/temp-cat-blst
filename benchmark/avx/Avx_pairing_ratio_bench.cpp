// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

#include <benchmark/benchmark.h>

#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/matrixnok_g1.h"
#include "blst/avx/curve/matrixnok_g2.h"
#include "blst/avx/pairing/matrixnok_pairing.h"
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

// Measures the online Miller work used by BTX after the G2 operand has been
// validated, normalized, and expanded into its immutable 68-line table.
static void BM_AvxPreparedPairing_MillerLoop(benchmark::State& state)
{
    const auto& p = blst::avx::G1AffinePoint::generator();
    const auto& q = blst::avx::G2AffinePoint::generator();
    blst::avx::PreparedLineTable lines = {};
    blst::avx::GT::precompute_lines(lines.lines, q);

    for (auto _ : state) {
        blst::avx::MillerResult out = {};
        blst::avx::GT::miller_loop_lines_finite(out, lines.lines, p);
        benchmark::DoNotOptimize(out);
    }
}

// Measures only the map from one valid, nonzero Miller result into GT. The
// Miller input is fixed outside the timed region, matching repeated independent
// final exponentiations in the FFT middle-product pointwise phase.
static void BM_AvxPairing_FinalExp(benchmark::State& state)
{
    const auto& p = blst::avx::G1AffinePoint::generator();
    const auto& q = blst::avx::G2AffinePoint::generator();
    blst::avx::PreparedLineTable lines = {};
    blst::avx::GT::precompute_lines(lines.lines, q);
    blst::avx::MillerResult miller = {};
    blst::avx::GT::miller_loop_lines_finite(miller, lines.lines, p);

    for (auto _ : state) {
        blst::avx::GT out = {};
        blst::avx::GT::final_exp(out, miller);
        benchmark::DoNotOptimize(out);
    }
}

// Measures the exact prepared-pairing composition used by one BTX frequency
// point. Line preparation and affine conversion remain outside the timed loop.
static void BM_AvxPreparedPairing_Full(benchmark::State& state)
{
    const auto& p = blst::avx::G1AffinePoint::generator();
    const auto& q = blst::avx::G2AffinePoint::generator();
    blst::avx::PreparedLineTable lines = {};
    blst::avx::GT::precompute_lines(lines.lines, q);

    for (auto _ : state) {
        blst::avx::MillerResult miller = {};
        blst::avx::GT out = {};
        blst::avx::GT::miller_loop_lines_finite(miller, lines.lines, p);
        blst::avx::GT::final_exp(out, miller);
        benchmark::DoNotOptimize(out);
    }
}

template<int64_t left_bound, int64_t right_bound>
static void BM_AvxFp12DenseMulImpl(
    benchmark::State& state,
    const blst::avx::avx_gt_detail::FP12::CommonElement<left_bound>& left,
    const blst::avx::avx_gt_detail::FP12::CommonElement<right_bound>& right)
{
    const auto& ring = blst::avx::avx_gt_detail::ring();
    const auto& fp12ring = blst::avx::avx_gt_detail::fp12ring();
    for (auto _ : state) {
        auto out = fp12ring.fp12_flat_mul(left, right, ring);
        benchmark::DoNotOptimize(out);
    }
}

// Measure dense Fp12 multiplication independently of pairing. Bounds one and
// six are the complete set of producer contracts used by the selected AVX
// backend. The operands are distinct nonidentity GT values and are constructed
// outside the timed region. These benchmarks intentionally instantiate all
// four bound combinations because each produces separate template code.
static void BM_AvxFp12DenseMul_Bounds1x1(benchmark::State& state)
{
    const auto& ring = blst::avx::avx_gt_detail::ring();
    const auto& fp12ring = blst::avx::avx_gt_detail::fp12ring();
    const auto left = blst::avx::GT::generator().data;
    const auto right = fp12ring.fp12_flat_frobenius_map1(left, ring);
    BM_AvxFp12DenseMulImpl(state, left, right);
}

static void BM_AvxFp12DenseMul_Bounds1x6(benchmark::State& state)
{
    const auto& ring = blst::avx::avx_gt_detail::ring();
    const auto& fp12ring = blst::avx::avx_gt_detail::fp12ring();
    const auto left = blst::avx::GT::generator().data;
    const auto right = fp12ring.recast_common(
        fp12ring.fp12_flat_frobenius_map1(left, ring));
    BM_AvxFp12DenseMulImpl(state, left, right);
}

static void BM_AvxFp12DenseMul_Bounds6x1(benchmark::State& state)
{
    const auto& ring = blst::avx::avx_gt_detail::ring();
    const auto& fp12ring = blst::avx::avx_gt_detail::fp12ring();
    const auto left = fp12ring.recast_common(blst::avx::GT::generator().data);
    const auto right = fp12ring.fp12_flat_frobenius_map1(
        blst::avx::GT::generator().data, ring);
    BM_AvxFp12DenseMulImpl(state, left, right);
}

static void BM_AvxFp12DenseMul_Bounds6x6(benchmark::State& state)
{
    const auto& ring = blst::avx::avx_gt_detail::ring();
    const auto& fp12ring = blst::avx::avx_gt_detail::fp12ring();
    const auto left = fp12ring.recast_common(blst::avx::GT::generator().data);
    const auto right = fp12ring.recast_common(
        fp12ring.fp12_flat_frobenius_map1(
            blst::avx::GT::generator().data, ring));
    BM_AvxFp12DenseMulImpl(state, left, right);
}

static void BM_DirectMatrixNoKPairing_MillerLoop(benchmark::State& state)
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::G2Fp2Ring fp2ring(ring);
    const blst::avx::matrixnok::Fp12Ring fp12ring{};
    const auto& public_p = blst::avx::G1AffinePoint::generator();
    const auto& public_q = blst::avx::G2AffinePoint::generator();
    const blst::avx::matrixnok::G1Affine p{public_p.X, public_p.Y};
    const blst::avx::matrixnok::G2Affine q{public_q.X, public_q.Y};
    const blst::avx::matrixnok::MillerEngine miller(p, q, ring);

    for (auto _ : state) {
        auto out = miller.miller_loop(ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(out);
    }
}

static void BM_DirectMatrixNoKPairing_Full(benchmark::State& state)
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::G2Fp2Ring fp2ring(ring);
    const blst::avx::matrixnok::Fp12Ring fp12ring{};
    const auto& public_p = blst::avx::G1AffinePoint::generator();
    const auto& public_q = blst::avx::G2AffinePoint::generator();
    const blst::avx::matrixnok::G1Affine p{public_p.X, public_p.Y};
    const blst::avx::matrixnok::G2Affine q{public_q.X, public_q.Y};
    const blst::avx::matrixnok::PairingInverter inverter{};

    for (auto _ : state) {
        const blst::avx::matrixnok::MillerEngine miller(p, q, ring);
        const auto f = miller.miller_loop(ring, fp2ring, fp12ring);
        auto out = blst::avx::final_exp(f, inverter, fp12ring, ring);
        benchmark::DoNotOptimize(out);
    }
}

BENCHMARK(BM_ScalarPairing_MillerLoop);
BENCHMARK(BM_ScalarPairing_Full);
BENCHMARK(BM_AvxPairing_MillerLoop);
BENCHMARK(BM_AvxPairing_Full);
BENCHMARK(BM_AvxPreparedPairing_MillerLoop);
BENCHMARK(BM_AvxPairing_FinalExp);
BENCHMARK(BM_AvxPreparedPairing_Full);
BENCHMARK(BM_AvxFp12DenseMul_Bounds1x1);
BENCHMARK(BM_AvxFp12DenseMul_Bounds1x6);
BENCHMARK(BM_AvxFp12DenseMul_Bounds6x1);
BENCHMARK(BM_AvxFp12DenseMul_Bounds6x6);
BENCHMARK(BM_DirectMatrixNoKPairing_MillerLoop);
BENCHMARK(BM_DirectMatrixNoKPairing_Full);

}  // namespace
