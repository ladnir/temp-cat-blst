#include <benchmark/benchmark.h>

#include <vector>

#include "blst/avx/curve/GT.h"
#include "blst/avx/pairing/pairing.h"
#include "test/avx/avx_test_utils.h"

using namespace blst::avx;

namespace {

using blst_avx_test::AvxAffineFp;
using blst_avx_test::AvxAffineFp2;
using blst_avx_test::AvxFp12Ring;
using blst_avx_test::AvxFp2Ring;
using blst_avx_test::AvxRing;

constexpr int kNumG1Samples = 1024;

std::vector<G1AffinePoint> make_g1_samples(const AvxRing& ring, int count)
{
    (void)ring;
    std::vector<G1AffinePoint> points;
    points.reserve(static_cast<size_t>(count));
    for (int i = 0; i < count; ++i) {
        const std::array<uint8_t, 32> scalar_bytes =
            blst_avx_test::seeded_scalar_bytes_le(static_cast<uint64_t>(i + 1));
        G1Point point = {};
        G1Point::mult(point, G1Point::generator(), scalar_bytes.data(), 256);

        G1AffinePoint affine = {};
        G1Point::to_affine(affine, point);
        points.push_back(affine);
    }
    return points;
}

static void BM_Pairing_BLST_Inverter(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp2Ring fp2ring(ring);
    const AvxFp12Ring fp12ring(ring);
    const AvxAffineFp g1 = blst_avx_test::make_avx_g1_generator_affine(ring);
    const AvxAffineFp2 g2 = blst_avx_test::make_avx_g2_generator_affine(ring);
    FP12BLSTInverter<AvxRing> inverter;

    for (auto _ : state) {
        auto result = pairing(g1, g2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_BLST_Inverter");
}

static void BM_Pairing_RNS_BLST_Inverter(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const AvxFp2Ring fp2ring(ring);
    const AvxFp12Ring fp12ring(ring);
    const AvxAffineFp g1 = blst_avx_test::make_avx_g1_generator_affine(ring);
    const AvxAffineFp2 g2 = blst_avx_test::make_avx_g2_generator_affine(ring);
    FP12RNS_BLSTInverter<AvxRing> inverter;

    for (auto _ : state) {
        auto result = pairing(g1, g2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_RNS_BLST_Inverter");
}

static void BM_MillerLoop_Direct(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const std::vector<G1AffinePoint> g1_points = make_g1_samples(ring, kNumG1Samples);
    const G2AffinePoint g2 = G2AffinePoint::generator();

    size_t point_idx = 0;
    for (auto _ : state) {
        MillerResult result = {};
        GT::miller_loop(result, &g2, &g1_points[point_idx]);
        benchmark::DoNotOptimize(result);

        point_idx += 1;
        if (point_idx == g1_points.size()) {
            point_idx = 0;
        }
    }

    state.SetItemsProcessed(state.iterations());
    state.SetLabel("AVX_MillerLoop_Direct");
}

static void BM_MillerLoop_PrecomputedLines(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const std::vector<G1AffinePoint> g1_points = make_g1_samples(ring, kNumG1Samples);
    const G2AffinePoint g2 = G2AffinePoint::generator();
    PreparedLineTable lines = {};
    GT::precompute_lines(lines.lines, g2);

    size_t point_idx = 0;
    for (auto _ : state) {
        MillerResult result = {};
        GT::miller_loop_lines(result, lines.lines, g1_points[point_idx]);
        benchmark::DoNotOptimize(result);

        point_idx += 1;
        if (point_idx == g1_points.size()) {
            point_idx = 0;
        }
    }

    state.SetItemsProcessed(state.iterations());
    state.SetLabel("AVX_MillerLoop_PrecomputedLines");
}

static void BM_MillerLoop_PrecomputeLines(benchmark::State& state)
{
    const G2AffinePoint g2 = G2AffinePoint::generator();

    for (auto _ : state) {
        PreparedLineTable lines = {};
        GT::precompute_lines(lines.lines, g2);
        benchmark::DoNotOptimize(lines);
    }

    state.SetItemsProcessed(state.iterations());
    state.SetLabel("AVX_MillerLoop_PrecomputeLines");
}

static void BM_MillerLoop_PrecomputedLines_FinalExp(benchmark::State& state)
{
    const BigInt modulus = blst_avx_test::modulus_bigint();
    const AvxRing ring(modulus);
    const std::vector<G1AffinePoint> g1_points = make_g1_samples(ring, kNumG1Samples);
    const G2AffinePoint g2 = G2AffinePoint::generator();
    PreparedLineTable lines = {};
    GT::precompute_lines(lines.lines, g2);

    size_t point_idx = 0;
    for (auto _ : state) {
        MillerResult miller = {};
        GT result = {};
        GT::miller_loop_lines(miller, lines.lines, g1_points[point_idx]);
        GT::final_exp(result, miller);
        benchmark::DoNotOptimize(result);

        point_idx += 1;
        if (point_idx == g1_points.size()) {
            point_idx = 0;
        }
    }

    state.SetItemsProcessed(state.iterations());
    state.SetLabel("AVX_MillerLoop_PrecomputedLines_FinalExp");
}

BENCHMARK(BM_Pairing_BLST_Inverter);
BENCHMARK(BM_Pairing_RNS_BLST_Inverter);
BENCHMARK(BM_MillerLoop_Direct);
BENCHMARK(BM_MillerLoop_PrecomputedLines);
BENCHMARK(BM_MillerLoop_PrecomputeLines);
BENCHMARK(BM_MillerLoop_PrecomputedLines_FinalExp);

}  // namespace
