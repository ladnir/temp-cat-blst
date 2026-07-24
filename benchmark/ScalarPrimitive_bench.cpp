#include <benchmark/benchmark.h>

#include <array>
#include <cstddef>
#include <cstdint>

#include "blst/curve/point.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fr.h"

namespace {

using namespace blst;

constexpr std::size_t LANES = 8;

Fp fp_from_u64(std::uint64_t value)
{
    Fp raw = {};
    raw[0] = value;
    Fp output = {};
    Fp::to_mont(output, raw);
    return output;
}

Fr fr_from_u64(std::uint64_t value)
{
    Fr raw = {};
    raw[0] = value;
    Fr output = {};
    Fr::to_mont(output, raw);
    return output;
}

std::array<G1Point, LANES> projective_points(std::uint8_t offset)
{
    std::array<G1Point, LANES> points = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        std::array<byte, 32> scalar = {};
        scalar[0] = static_cast<byte>(offset + lane);
        G1Point::mult(
            points[lane], G1Point::generator(), scalar.data(), 8);
    }
    return points;
}

std::array<G1AffinePoint, LANES> affine_points(std::uint8_t offset)
{
    const auto projective = projective_points(offset);
    std::array<G1AffinePoint, LANES> affine = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        G1Point::to_affine(affine[lane], projective[lane]);
    }
    return affine;
}

void set_items(benchmark::State& state)
{
    state.SetItemsProcessed(
        static_cast<std::int64_t>(state.iterations() * LANES));
}

void BM_Primitive_FpAdd(benchmark::State& state)
{
    std::array<Fp, LANES> accumulators = {};
    std::array<Fp, LANES> addends = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        accumulators[lane] = fp_from_u64(11 + lane);
        addends[lane] = fp_from_u64(101 + 2 * lane);
    }

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            Fp::add(
                accumulators[lane], accumulators[lane], addends[lane]);
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_FpAdd);

void BM_Primitive_FpMul(benchmark::State& state)
{
    std::array<Fp, LANES> accumulators = {};
    std::array<Fp, LANES> factors = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        accumulators[lane] = fp_from_u64(11 + lane);
        factors[lane] = fp_from_u64(101 + 2 * lane);
    }

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            Fp::mul(
                accumulators[lane], accumulators[lane], factors[lane]);
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_FpMul);

void BM_Primitive_FpSquare(benchmark::State& state)
{
    std::array<Fp, LANES> accumulators = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        accumulators[lane] = fp_from_u64(11 + lane);
    }

    for (auto _ : state) {
        for (auto& value : accumulators) {
            Fp::sqr(value, value);
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_FpSquare);

void BM_Primitive_FpInverse(benchmark::State& state)
{
    std::array<Fp, LANES> accumulators = {};
    std::array<Fp, LANES> output = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        accumulators[lane] = fp_from_u64(11 + lane);
    }

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            Fp::flt_reciprocal(output[lane], accumulators[lane]);
            accumulators[lane] = output[lane];
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_FpInverse);

void BM_Primitive_FpReciprocal(benchmark::State& state)
{
    std::array<Fp, LANES> accumulators = {};
    std::array<Fp, LANES> output = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        accumulators[lane] = fp_from_u64(11 + lane);
    }

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            Fp::reciprocal(output[lane], accumulators[lane]);
            accumulators[lane] = output[lane];
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_FpReciprocal);

void BM_Primitive_FrAdd(benchmark::State& state)
{
    std::array<Fr, LANES> accumulators = {};
    std::array<Fr, LANES> addends = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        accumulators[lane] = fr_from_u64(11 + lane);
        addends[lane] = fr_from_u64(101 + 2 * lane);
    }

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            Fr::add(
                accumulators[lane], accumulators[lane], addends[lane]);
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_FrAdd);

void BM_Primitive_FrMul(benchmark::State& state)
{
    std::array<Fr, LANES> accumulators = {};
    std::array<Fr, LANES> factors = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        accumulators[lane] = fr_from_u64(11 + lane);
        factors[lane] = fr_from_u64(101 + 2 * lane);
    }

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            Fr::mul(
                accumulators[lane], accumulators[lane], factors[lane]);
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_FrMul);

void BM_Primitive_FrSquare(benchmark::State& state)
{
    std::array<Fr, LANES> accumulators = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        accumulators[lane] = fr_from_u64(11 + lane);
    }

    for (auto _ : state) {
        for (auto& value : accumulators) {
            Fr::sqr(value, value);
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_FrSquare);

void BM_Primitive_FrInverse(benchmark::State& state)
{
    std::array<Fr, LANES> accumulators = {};
    std::array<Fr, LANES> output = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        accumulators[lane] = fr_from_u64(11 + lane);
    }

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            Fr::reciprocal(output[lane], accumulators[lane]);
            accumulators[lane] = output[lane];
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_FrInverse);

void BM_Primitive_G1ProjectiveAdd(benchmark::State& state)
{
    auto accumulators = projective_points(3);
    const auto addends = projective_points(41);
    std::array<G1Point, LANES> output = {};

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            G1Point::add(
                output[lane], accumulators[lane], addends[lane]);
            accumulators[lane] = output[lane];
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_G1ProjectiveAdd);

void BM_Primitive_G1MixedAdd(benchmark::State& state)
{
    auto accumulators = projective_points(3);
    const auto addends = affine_points(41);
    std::array<G1Point, LANES> output = {};

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            G1Point::add(
                output[lane], accumulators[lane], addends[lane]);
            accumulators[lane] = output[lane];
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_G1MixedAdd);

void BM_Primitive_G1Double(benchmark::State& state)
{
    auto accumulators = projective_points(3);
    std::array<G1Point, LANES> output = {};

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            G1Point::doubling(output[lane], accumulators[lane]);
            accumulators[lane] = output[lane];
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_G1Double);

void BM_Primitive_G1XyzzMixedAdd(benchmark::State& state)
{
    const auto initial = affine_points(3);
    const auto addends = affine_points(41);
    std::array<G1PointZZ, LANES> accumulators = {};
    std::array<G1PointZZ, LANES> output = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        G1PointZZ::from_affine(accumulators[lane], initial[lane]);
    }

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            G1PointZZ::dadd(
                output[lane], accumulators[lane], addends[lane], false);
            accumulators[lane] = output[lane];
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_G1XyzzMixedAdd);

void BM_Primitive_G1XyzzAdd(benchmark::State& state)
{
    const auto initial = affine_points(3);
    const auto other = affine_points(41);
    std::array<G1PointZZ, LANES> accumulators = {};
    std::array<G1PointZZ, LANES> addends = {};
    std::array<G1PointZZ, LANES> output = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        G1PointZZ::from_affine(accumulators[lane], initial[lane]);
        G1PointZZ::from_affine(addends[lane], other[lane]);
    }

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            G1PointZZ::dadd(
                output[lane], accumulators[lane], addends[lane]);
            accumulators[lane] = output[lane];
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_G1XyzzAdd);

void BM_Primitive_G1ScalarMul(benchmark::State& state)
{
    auto accumulators = projective_points(3);
    std::array<G1Point, LANES> output = {};
    const Fr multiplier = fr_from_u64(6364136223846793005ULL);
    const Fr increment = fr_from_u64(1442695040888963407ULL);
    Fr scalar_state = fr_from_u64(17);
    for (std::size_t round = 0; round < 16; ++round) {
        Fr::mul(scalar_state, scalar_state, multiplier);
        Fr::add(scalar_state, scalar_state, increment);
    }
    std::array<std::array<byte, 32>, LANES> scalars = {};
    for (std::size_t lane = 0; lane < LANES; ++lane) {
        Fr::mul(scalar_state, scalar_state, multiplier);
        Fr::add(scalar_state, scalar_state, increment);
        Fr::to_scaler(scalars[lane].data(), scalar_state);
    }

    for (auto _ : state) {
        for (std::size_t lane = 0; lane < LANES; ++lane) {
            G1Point::mult(
                output[lane], accumulators[lane],
                scalars[lane].data(), 254);
            accumulators[lane] = output[lane];
        }
        benchmark::DoNotOptimize(accumulators);
    }
    set_items(state);
}
BENCHMARK(BM_Primitive_G1ScalarMul);

} // namespace
