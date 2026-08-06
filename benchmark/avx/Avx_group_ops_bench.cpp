#include <benchmark/benchmark.h>

#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/multi_scalar.h"
#include "blst/avx/curve/point.h"
#include "blst/scalar/fields/Fp.h"
#include "blst/scalar/fields/Fp2.h"
#include "test/avx/avx_test_utils.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

using namespace blst::avx;

namespace {

using blst::BLST_ERROR;
using blst::bool_t;
using ScalarBytes = std::array<blst::byte, 32>;
constexpr size_t kBatchSubgroupWidth = 40;

ScalarBytes scalar_bytes(uint64_t seed)
{
    const auto bytes = blst_avx_test::seeded_scalar_bytes_le(seed);
    ScalarBytes out = {};
    for (size_t i = 0; i < out.size(); ++i) {
        out[i] = bytes[i];
    }
    return out;
}

G1Point make_g1(uint64_t seed)
{
    const ScalarBytes scalar = scalar_bytes(seed);
    G1Point out = {};
    G1Point::mult(out, G1Point::generator(), scalar.data(), 256);
    return out;
}

G2Point make_g2(uint64_t seed)
{
    const ScalarBytes scalar = scalar_bytes(seed);
    G2Point out = {};
    G2Point::mult(out, G2Point::generator(), scalar.data(), 256);
    return out;
}

GT make_gt(uint64_t g1_seed)
{
    G1Point g1 = make_g1(g1_seed);
    G1AffinePoint g1_affine = {};
    G1Point::to_affine(g1_affine, g1);

    MillerResult miller = {};
    GT out = {};
    GT::miller_loop(miller, &G2AffinePoint::generator(), &g1_affine);
    GT::final_exp(out, miller);
    return out;
}

blst::scalar::G1AffinePoint make_scalar_g1_affine(uint64_t seed)
{
    blst::scalar::G1AffinePoint out = {};
    blst::scalar::G1Point::to_affine(out, blst_avx_test::make_scalar_g1_projective_from_seed(seed));
    return out;
}

blst::scalar::G2AffinePoint make_scalar_g2_affine(uint64_t seed)
{
    blst::scalar::G2AffinePoint out = {};
    blst::scalar::G2Point::to_affine(out, blst_avx_test::make_scalar_g2_projective_from_seed(seed));
    return out;
}

blst::scalar::GT make_scalar_gt(uint64_t g1_seed)
{
    const blst::scalar::G1AffinePoint p = make_scalar_g1_affine(g1_seed);
    blst::scalar::MillerResult miller = {};
    blst::scalar::GT out = {};
    blst::scalar::GT::miller_loop(miller, &blst::scalar::G2AffinePoint::generator(), &p);
    blst::scalar::GT::final_exp(out, miller);
    return out;
}

blst::scalar::Fp make_scalar_fp_square(uint64_t seed)
{
    const blst::scalar::G1AffinePoint point = make_scalar_g1_affine(seed);
    blst::scalar::Fp out = {};
    blst::scalar::Fp::sqr(out, point.X);
    return out;
}

blst::scalar::Fp2 make_scalar_fp2_square(uint64_t seed)
{
    const blst::scalar::G2AffinePoint point = make_scalar_g2_affine(seed);
    blst::scalar::Fp2 out = {};
    blst::scalar::Fp2::sqr(out, point.X);
    return out;
}

std::vector<G1Point> make_g1_batch(size_t count)
{
    std::vector<G1Point> points;
    points.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        points.push_back(make_g1(42 + i));
    }
    return points;
}

std::vector<G2Point> make_g2_batch(size_t count)
{
    std::vector<G2Point> points;
    points.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        points.push_back(make_g2(42 + i));
    }
    return points;
}

std::vector<blst::byte> make_scalar_batch(size_t count)
{
    std::vector<blst::byte> scalars(count * 32);
    for (size_t i = 0; i < count; ++i) {
        const ScalarBytes scalar = scalar_bytes(1000 + i);
        std::copy(scalar.begin(), scalar.end(), scalars.begin() + static_cast<std::ptrdiff_t>(i * 32));
    }
    return scalars;
}

void batch_to_affine_g1_nonzero(
    G1AffinePoint* out,
    const G1Point* points,
    size_t count,
    std::vector<G1Point::Field_type>& prefix_products);

void batch_to_affine_g2_nonzero(
    G2AffinePoint* out,
    const G2Point* points,
    size_t count,
    std::vector<G2Point::Field_type>& prefix_products);

std::vector<G1AffinePoint> make_g1_affine_batch(size_t count)
{
    std::vector<G1Point> points = make_g1_batch(count);
    std::vector<G1AffinePoint> affine(count);
    std::vector<G1Point::Field_type> prefix_products(count);
    batch_to_affine_g1_nonzero(affine.data(), points.data(), count, prefix_products);
    return affine;
}

std::vector<G2AffinePoint> make_g2_affine_batch(size_t count)
{
    std::vector<G2Point> points = make_g2_batch(count);
    std::vector<G2AffinePoint> affine(count);
    std::vector<G2Point::Field_type> prefix_products(count);
    batch_to_affine_g2_nonzero(affine.data(), points.data(), count, prefix_products);
    return affine;
}

std::vector<G1AffinePoint> make_g1_affine_doubles(const std::vector<G1AffinePoint>& points)
{
    std::vector<G1Point> doubled(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        G1Point projective = {};
        G1Point::from_affine(projective, points[i]);
        G1Point::doubling(doubled[i], projective);
    }

    std::vector<G1AffinePoint> affine(points.size());
    std::vector<G1Point::Field_type> prefix_products(points.size());
    batch_to_affine_g1_nonzero(affine.data(), doubled.data(), doubled.size(), prefix_products);
    return affine;
}

std::vector<G2AffinePoint> make_g2_affine_doubles(const std::vector<G2AffinePoint>& points)
{
    std::vector<G2Point> doubled(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        G2Point projective = {};
        G2Point::from_affine(projective, points[i]);
        G2Point::doubling(doubled[i], projective);
    }

    std::vector<G2AffinePoint> affine(points.size());
    std::vector<G2Point::Field_type> prefix_products(points.size());
    batch_to_affine_g2_nonzero(affine.data(), doubled.data(), doubled.size(), prefix_products);
    return affine;
}

std::vector<uint8_t> make_mod3_coefficients(size_t count, size_t width)
{
    uint64_t rng = 0x9e3779b97f4a7c15ull ^ (count << 32) ^ width;
    std::vector<uint8_t> coefficients(count * width);
    for (size_t i = 0; i < coefficients.size(); ++i) {
        rng ^= rng >> 12;
        rng ^= rng << 25;
        rng ^= rng >> 27;
        coefficients[i] = static_cast<uint8_t>((rng * 0x2545f4914f6cdd1dull) % 3);
    }
    return coefficients;
}

std::vector<uint8_t> make_binary_coefficients(size_t count, size_t width)
{
    uint64_t rng = 0xd1b54a32d192ed03ull ^ (count << 32) ^ width;
    std::vector<uint8_t> coefficients(count * width);
    for (size_t i = 0; i < coefficients.size(); ++i) {
        rng ^= rng >> 12;
        rng ^= rng << 25;
        rng ^= rng >> 27;
        coefficients[i] = static_cast<uint8_t>((rng * 0x2545f4914f6cdd1dull) & 1);
    }
    return coefficients;
}

void batch_to_affine_g1_nonzero(
    G1AffinePoint* out,
    const G1Point* points,
    size_t count,
    std::vector<G1Point::Field_type>& prefix_products)
{
    const auto& ring = avx_curve_detail::curve_ring();

    prefix_products[0] = points[0].Z;
    for (size_t i = 1; i < count; ++i) {
        prefix_products[i] = ring.modmul(prefix_products[i - 1], points[i].Z);
    }

    auto running_inverse = avx_curve_detail::curve_inverse(prefix_products[count - 1]);
    for (size_t i = count; i-- > 0;) {
        const auto z_inverse = i == 0
            ? running_inverse
            : ring.modmul(prefix_products[i - 1], running_inverse);
        running_inverse = ring.modmul(running_inverse, points[i].Z);

        auto [x, y] = ring.batch_reduce_expand(
            ring.prep_left(points[i].X) * z_inverse,
            ring.prep_left(points[i].Y) * z_inverse);
        out[i] = G1AffinePoint(x, y);
    }
}

void batch_to_affine_g2_nonzero(
    G2AffinePoint* out,
    const G2Point* points,
    size_t count,
    std::vector<G2Point::Field_type>& prefix_products)
{
    const auto& fp2ring = avx_curve_detail::curve_fp2_ring();

    prefix_products[0] = points[0].Z;
    for (size_t i = 1; i < count; ++i) {
        prefix_products[i] = fp2ring.mul(prefix_products[i - 1], points[i].Z)[0];
    }

    auto running_inverse = avx_curve_detail::curve_fp2_inverse(prefix_products[count - 1]);
    for (size_t i = count; i-- > 0;) {
        const auto z_inverse = i == 0
            ? running_inverse
            : fp2ring.mul(prefix_products[i - 1], running_inverse)[0];
        running_inverse = fp2ring.mul(running_inverse, points[i].Z)[0];

        out[i] = G2AffinePoint(
            fp2ring.mul(points[i].X, z_inverse)[0],
            fp2ring.mul(points[i].Y, z_inverse)[0]);
    }
}

bool_t batch_in_subgroup_g1_assume_on_curve(
    const G1AffinePoint* points,
    size_t count,
    const uint8_t* coefficients,
    size_t width)
{
    bool_t ok = true;
    for (size_t j = 0; j < width; ++j) {
        G1Point sum = G1Point::infinity();
        const uint8_t* row = coefficients + j * count;
        for (size_t i = 0; i < count; ++i) {
            if (row[i] == 1) {
                G1Point::add(sum, sum, points[i]);
            } else if (row[i] == 2) {
                G1Point::add(sum, sum, points[i]);
                G1Point::add(sum, sum, points[i]);
            }
        }
        ok &= G1Point::in_subgroup(sum);
    }
    return ok;
}

bool_t batch_in_subgroup_g1_precomputed_twos_assume_on_curve(
    const G1AffinePoint* points,
    const G1AffinePoint* doubles,
    size_t count,
    const uint8_t* coefficients,
    size_t width)
{
    bool_t ok = true;
    for (size_t j = 0; j < width; ++j) {
        G1Point sum = G1Point::infinity();
        const uint8_t* row = coefficients + j * count;
        for (size_t i = 0; i < count; ++i) {
            if (row[i] == 1) {
                G1Point::add(sum, sum, points[i]);
            } else if (row[i] == 2) {
                G1Point::add(sum, sum, doubles[i]);
            }
        }
        ok &= G1Point::in_subgroup(sum);
    }
    return ok;
}

bool_t batch_in_subgroup_g1_binary_assume_on_curve(
    const G1AffinePoint* points,
    size_t count,
    const uint8_t* coefficients,
    size_t width)
{
    bool_t ok = true;
    for (size_t j = 0; j < width; ++j) {
        G1Point sum = G1Point::infinity();
        const uint8_t* row = coefficients + j * count;
        for (size_t i = 0; i < count; ++i) {
            if (row[i] != 0) {
                G1Point::add(sum, sum, points[i]);
            }
        }
        ok &= G1Point::in_subgroup(sum);
    }
    return ok;
}

bool_t batch_in_subgroup_g2_assume_on_curve(
    const G2AffinePoint* points,
    size_t count,
    const uint8_t* coefficients,
    size_t width)
{
    bool_t ok = true;
    for (size_t j = 0; j < width; ++j) {
        G2Point sum = G2Point::infinity();
        const uint8_t* row = coefficients + j * count;
        for (size_t i = 0; i < count; ++i) {
            if (row[i] == 1) {
                G2Point::add(sum, sum, points[i]);
            } else if (row[i] == 2) {
                G2Point::add(sum, sum, points[i]);
                G2Point::add(sum, sum, points[i]);
            }
        }
        ok &= G2Point::in_subgroup(sum);
    }
    return ok;
}

bool_t batch_in_subgroup_g2_precomputed_twos_assume_on_curve(
    const G2AffinePoint* points,
    const G2AffinePoint* doubles,
    size_t count,
    const uint8_t* coefficients,
    size_t width)
{
    bool_t ok = true;
    for (size_t j = 0; j < width; ++j) {
        G2Point sum = G2Point::infinity();
        const uint8_t* row = coefficients + j * count;
        for (size_t i = 0; i < count; ++i) {
            if (row[i] == 1) {
                G2Point::add(sum, sum, points[i]);
            } else if (row[i] == 2) {
                G2Point::add(sum, sum, doubles[i]);
            }
        }
        ok &= G2Point::in_subgroup(sum);
    }
    return ok;
}

bool_t batch_in_subgroup_g2_binary_assume_on_curve(
    const G2AffinePoint* points,
    size_t count,
    const uint8_t* coefficients,
    size_t width)
{
    bool_t ok = true;
    for (size_t j = 0; j < width; ++j) {
        G2Point sum = G2Point::infinity();
        const uint8_t* row = coefficients + j * count;
        for (size_t i = 0; i < count; ++i) {
            if (row[i] != 0) {
                G2Point::add(sum, sum, points[i]);
            }
        }
        ok &= G2Point::in_subgroup(sum);
    }
    return ok;
}

static void BM_AvxG1_Add(benchmark::State& state)
{
    const G1Point a = make_g1(42);
    const G1Point b = make_g1(123);
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        G1Point out = {};
        G1Point::add(out, a, b);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1_Add");
}

static void BM_AvxG1_AddMixed(benchmark::State& state)
{
    const G1Point a = make_g1(42);
    G1AffinePoint b = {};
    G1Point::to_affine(b, make_g1(123));
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        G1Point out = {};
        G1Point::add(out, a, b);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1_AddMixed_JA");
}

static void BM_AvxG1_Dadd(benchmark::State& state)
{
    const G1Point a = make_g1(42);
    const G1Point b = make_g1(123);
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        G1Point out = {};
        G1Point::dadd(out, a, b);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1_Dadd");
}

static void BM_AvxG1_DaddMixed(benchmark::State& state)
{
    const G1Point a = make_g1(42);
    G1AffinePoint b = {};
    G1Point::to_affine(b, make_g1(123));
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        G1Point out = {};
        G1Point::dadd(out, a, b);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1_DaddMixed_JA");
}

static void BM_AvxG1_Double(benchmark::State& state)
{
    const G1Point a = make_g1(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        G1Point out = {};
        G1Point::doubling(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1_Double");
}

static void BM_AvxG1_Mult(benchmark::State& state)
{
    const G1Point a = make_g1(42);
    const ScalarBytes scalar = scalar_bytes(123);
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(scalar);

    for (auto _ : state) {
        G1Point out = {};
        G1Point::mult(out, a, scalar.data(), 256);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1_Mult");
}

static void BM_AvxG1_MsmPippenger(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G1AffinePoint> points = make_g1_affine_batch(count);
    const std::vector<blst::byte> scalars = make_scalar_batch(count);
    std::vector<AvxG1PippengerBucket> buckets(
        avx_g1_mult_pippenger_scratch_sizeof(count) / sizeof(AvxG1PippengerBucket));
    benchmark::DoNotOptimize(points.data());
    benchmark::DoNotOptimize(scalars.data());

    for (auto _ : state) {
        G1Point out = {};
        avx_g1_mult_pippenger(out, points.data(), scalars.data(), count, 256, buckets.data());
        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG1_MsmPippenger");
}

static void BM_AvxG1_OnCurve(benchmark::State& state)
{
    const G1Point a = make_g1(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = G1Point::on_curve(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1_OnCurve");
}

static void BM_AvxG1_InSubgroup(benchmark::State& state)
{
    const G1Point a = make_g1(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = G1Point::in_subgroup(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1_InSubgroup");
}

static void BM_AvxG1_Serialize(benchmark::State& state)
{
    const G1Point a = make_g1(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        unsigned char out[96] = {};
        G1Point::serialize(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1_Serialize");
}

static void BM_AvxG1Affine_OnCurve(benchmark::State& state)
{
    G1AffinePoint a = {};
    G1Point::to_affine(a, make_g1(42));
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = G1AffinePoint::on_curve(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine_OnCurve");
}

static void BM_AvxG1Affine_InSubgroup(benchmark::State& state)
{
    G1AffinePoint a = {};
    G1Point::to_affine(a, make_g1(42));
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = G1AffinePoint::in_subgroup(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine_InSubgroup");
}

static void BM_AvxG1Affine_Serialize(benchmark::State& state)
{
    G1AffinePoint a = {};
    G1Point::to_affine(a, make_g1(42));
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        unsigned char out[96] = {};
        G1AffinePoint::serialize(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine_Serialize");
}

static void BM_AvxG1Affine_Deserialize(benchmark::State& state)
{
    G1AffinePoint a = {};
    G1Point::to_affine(a, make_g1(42));
    unsigned char serialized[96] = {};
    G1AffinePoint::serialize(serialized, a);
    benchmark::DoNotOptimize(serialized);

    for (auto _ : state) {
        G1AffinePoint out = {};
        const BLST_ERROR err = G1AffinePoint::deserialize(out, serialized);
        benchmark::DoNotOptimize(err);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine_Deserialize");
}

static void BM_AvxG2_Add(benchmark::State& state)
{
    const G2Point a = make_g2(42);
    const G2Point b = make_g2(123);
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        G2Point out = {};
        G2Point::add(out, a, b);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2_Add");
}

static void BM_AvxG2_AddMixed(benchmark::State& state)
{
    const G2Point a = make_g2(42);
    G2AffinePoint b = {};
    G2Point::to_affine(b, make_g2(123));
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        G2Point out = {};
        G2Point::add(out, a, b);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2_AddMixed_JA");
}

static void BM_AvxG2_Dadd(benchmark::State& state)
{
    const G2Point a = make_g2(42);
    const G2Point b = make_g2(123);
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        G2Point out = {};
        G2Point::dadd(out, a, b);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2_Dadd");
}

static void BM_AvxG2_DaddMixed(benchmark::State& state)
{
    const G2Point a = make_g2(42);
    G2AffinePoint b = {};
    G2Point::to_affine(b, make_g2(123));
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        G2Point out = {};
        G2Point::dadd(out, a, b);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2_DaddMixed_JA");
}

static void BM_AvxG2_Double(benchmark::State& state)
{
    const G2Point a = make_g2(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        G2Point out = {};
        G2Point::doubling(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2_Double");
}

static void BM_AvxG2_Mult(benchmark::State& state)
{
    const G2Point a = make_g2(42);
    const ScalarBytes scalar = scalar_bytes(123);
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(scalar);

    for (auto _ : state) {
        G2Point out = {};
        G2Point::mult(out, a, scalar.data(), 256);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2_Mult");
}

static void BM_AvxG2_MsmPippenger(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G2AffinePoint> points = make_g2_affine_batch(count);
    const std::vector<blst::byte> scalars = make_scalar_batch(count);
    std::vector<AvxG2PippengerBucket> buckets(
        avx_g2_mult_pippenger_scratch_sizeof(count) / sizeof(AvxG2PippengerBucket),
        G2Point::infinity());
    benchmark::DoNotOptimize(points.data());
    benchmark::DoNotOptimize(scalars.data());

    for (auto _ : state) {
        G2Point out = {};
        avx_g2_mult_pippenger(out, points.data(), scalars.data(), count, 256, buckets.data());
        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG2_MsmPippenger");
}

static void BM_AvxG2_OnCurve(benchmark::State& state)
{
    const G2Point a = make_g2(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = G2Point::on_curve(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2_OnCurve");
}

static void BM_AvxG2_InSubgroup(benchmark::State& state)
{
    const G2Point a = make_g2(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = G2Point::in_subgroup(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2_InSubgroup");
}

static void BM_AvxG2_Serialize(benchmark::State& state)
{
    const G2Point a = make_g2(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        unsigned char out[192] = {};
        G2Point::serialize(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2_Serialize");
}

static void BM_AvxG2Affine_OnCurve(benchmark::State& state)
{
    G2AffinePoint a = {};
    G2Point::to_affine(a, make_g2(42));
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = G2AffinePoint::on_curve(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine_OnCurve");
}

static void BM_AvxG2Affine_InSubgroup(benchmark::State& state)
{
    G2AffinePoint a = {};
    G2Point::to_affine(a, make_g2(42));
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = G2AffinePoint::in_subgroup(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine_InSubgroup");
}

static void BM_AvxG2Affine_Serialize(benchmark::State& state)
{
    G2AffinePoint a = {};
    G2Point::to_affine(a, make_g2(42));
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        unsigned char out[192] = {};
        G2AffinePoint::serialize(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine_Serialize");
}

static void BM_AvxG2Affine_Deserialize(benchmark::State& state)
{
    G2AffinePoint a = {};
    G2Point::to_affine(a, make_g2(42));
    unsigned char serialized[192] = {};
    G2AffinePoint::serialize(serialized, a);
    benchmark::DoNotOptimize(serialized);

    for (auto _ : state) {
        G2AffinePoint out = {};
        const BLST_ERROR err = G2AffinePoint::deserialize(out, serialized);
        benchmark::DoNotOptimize(err);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine_Deserialize");
}

static void BM_AvxGT_Add(benchmark::State& state)
{
    const GT a = make_gt(42);
    const GT b = make_gt(123);
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        GT out = {};
        GT::add(out, a, b);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT_Add");
}

static void BM_AvxGT_Dadd(benchmark::State& state)
{
    const GT a = make_gt(42);
    const GT b = make_gt(123);
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        GT out = {};
        GT::dadd(out, a, b);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT_Dadd");
}

static void BM_AvxGT_Double(benchmark::State& state)
{
    const GT a = make_gt(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        GT out = {};
        GT::doubling(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT_Double");
}

static void BM_AvxGT_Mult(benchmark::State& state)
{
    const GT a = make_gt(42);
    const ScalarBytes scalar = scalar_bytes(123);
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(scalar);

    for (auto _ : state) {
        GT out = {};
        GT::mult(out, a, scalar.data(), 256);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT_Mult");
}

static void BM_AvxGT_InGroup(benchmark::State& state)
{
    const GT a = make_gt(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = GT::in_group(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT_InGroup");
}

static void BM_AvxGT_Serialize(benchmark::State& state)
{
    const GT a = make_gt(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        unsigned char out[12 * 48] = {};
        GT::serialize(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT_Serialize");
}

static void BM_AvxGT_Deserialize(benchmark::State& state)
{
    const GT a = make_gt(42);
    unsigned char serialized[12 * 48] = {};
    GT::serialize(serialized, a);
    benchmark::DoNotOptimize(serialized);

    for (auto _ : state) {
        GT out = {};
        benchmark::DoNotOptimize(GT::deserialize(out, serialized));
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT_Deserialize");
}

static void BM_AvxG1_Breakdown_ProjectiveToAffine(benchmark::State& state)
{
    const G1Point a = make_g1(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        G1AffinePoint out = {};
        G1Point::to_affine(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1.Breakdown.ProjectiveToAffine");
}

static void BM_AvxG1_Breakdown_ProjectiveToAffineRepeated(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G1Point> points = make_g1_batch(count);
    std::vector<G1AffinePoint> out(count);
    benchmark::DoNotOptimize(points.data());

    for (auto _ : state) {
        for (size_t i = 0; i < count; ++i) {
            G1Point::to_affine(out[i], points[i]);
        }
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG1.Breakdown.ProjectiveToAffineRepeated");
}

static void BM_AvxG1_Breakdown_BatchProjectiveToAffine(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G1Point> points = make_g1_batch(count);
    std::vector<G1AffinePoint> out(count);
    std::vector<G1Point::Field_type> prefix_products(count);
    benchmark::DoNotOptimize(points.data());

    for (auto _ : state) {
        batch_to_affine_g1_nonzero(out.data(), points.data(), count, prefix_products);
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG1.Breakdown.BatchProjectiveToAffine");
}

static void BM_AvxG1Affine_Breakdown_AvxToScalar(benchmark::State& state)
{
    G1AffinePoint a = {};
    G1Point::to_affine(a, make_g1(42));
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const auto out = avx_curve_detail::curve_g1_affine_to_scalar({a.X, a.Y});
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine.Breakdown.AvxToScalar");
}

static void BM_AvxG1Affine_Breakdown_ScalarSerialize(benchmark::State& state)
{
    const blst::scalar::G1AffinePoint a = make_scalar_g1_affine(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        unsigned char out[96] = {};
        blst::scalar::G1AffinePoint::serialize(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine.Breakdown.ScalarSerialize");
}

static void BM_AvxG1Affine_Breakdown_ScalarCompress(benchmark::State& state)
{
    const blst::scalar::G1AffinePoint a = make_scalar_g1_affine(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        unsigned char out[48] = {};
        blst::scalar::G1AffinePoint::compress(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine.Breakdown.ScalarCompress");
}

static void BM_AvxG1Affine_Breakdown_BytesToScalar(benchmark::State& state)
{
    const blst::scalar::G1AffinePoint a = make_scalar_g1_affine(42);
    unsigned char serialized[96] = {};
    blst::scalar::G1AffinePoint::serialize(serialized, a);
    benchmark::DoNotOptimize(serialized);

    for (auto _ : state) {
        blst::scalar::G1AffinePoint out = {};
        const BLST_ERROR err = blst::scalar::G1AffinePoint::deserialize(out, serialized);
        benchmark::DoNotOptimize(err);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine.Breakdown.BytesToScalar");
}

static void BM_AvxG1Affine_Breakdown_CompressedBytesToScalar(benchmark::State& state)
{
    const blst::scalar::G1AffinePoint a = make_scalar_g1_affine(42);
    unsigned char compressed[48] = {};
    blst::scalar::G1AffinePoint::compress(compressed, a);
    benchmark::DoNotOptimize(compressed);

    for (auto _ : state) {
        blst::scalar::G1AffinePoint out = {};
        const BLST_ERROR err = blst::scalar::G1AffinePoint::uncompress(out, compressed);
        benchmark::DoNotOptimize(err);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine.Breakdown.CompressedBytesToScalar");
}

static void BM_AvxG1Affine_Breakdown_ScalarToAvx(benchmark::State& state)
{
    const blst::scalar::G1AffinePoint a = make_scalar_g1_affine(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const G1AffinePoint out(
            avx_curve_detail::scalar_fp_to_ring(a.X),
            avx_curve_detail::scalar_fp_to_ring(a.Y));
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine.Breakdown.ScalarToAvx");
}

static void BM_AvxG1Affine_Breakdown_ScalarOnCurve(benchmark::State& state)
{
    const blst::scalar::G1AffinePoint a = make_scalar_g1_affine(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = blst::scalar::G1AffinePoint::on_curve(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine.Breakdown.ScalarOnCurve");
}

static void BM_AvxG1Affine_Breakdown_ScalarInSubgroup(benchmark::State& state)
{
    const blst::scalar::G1AffinePoint a = make_scalar_g1_affine(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = blst::scalar::G1AffinePoint::in_subgroup(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG1Affine.Breakdown.ScalarInSubgroup");
}

static void BM_AvxG1Affine_Breakdown_InSubgroupRepeated(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G1AffinePoint> points = make_g1_affine_batch(count);
    benchmark::DoNotOptimize(points.data());

    for (auto _ : state) {
        bool_t ok = true;
        for (size_t i = 0; i < count; ++i) {
            ok &= G1AffinePoint::in_subgroup(points[i]);
        }
        benchmark::DoNotOptimize(ok);
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG1Affine.Breakdown.InSubgroupRepeated");
}

static void BM_AvxG1Affine_Breakdown_BatchInSubgroupMod3W40(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G1AffinePoint> points = make_g1_affine_batch(count);
    const std::vector<uint8_t> coefficients = make_mod3_coefficients(count, kBatchSubgroupWidth);
    benchmark::DoNotOptimize(points.data());
    benchmark::DoNotOptimize(coefficients.data());

    for (auto _ : state) {
        const bool_t ok = batch_in_subgroup_g1_assume_on_curve(
            points.data(),
            count,
            coefficients.data(),
            kBatchSubgroupWidth);
        benchmark::DoNotOptimize(ok);
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG1Affine.Breakdown.BatchInSubgroupMod3W40");
}

static void BM_AvxG1Affine_Breakdown_BatchInSubgroupMod3Precompute2W40(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G1AffinePoint> points = make_g1_affine_batch(count);
    const std::vector<G1AffinePoint> doubles = make_g1_affine_doubles(points);
    const std::vector<uint8_t> coefficients = make_mod3_coefficients(count, kBatchSubgroupWidth);
    benchmark::DoNotOptimize(points.data());
    benchmark::DoNotOptimize(doubles.data());
    benchmark::DoNotOptimize(coefficients.data());

    for (auto _ : state) {
        const bool_t ok = batch_in_subgroup_g1_precomputed_twos_assume_on_curve(
            points.data(),
            doubles.data(),
            count,
            coefficients.data(),
            kBatchSubgroupWidth);
        benchmark::DoNotOptimize(ok);
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG1Affine.Breakdown.BatchInSubgroupMod3Precompute2W40");
}

static void BM_AvxG1Affine_Breakdown_BatchInSubgroupBinaryW40(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G1AffinePoint> points = make_g1_affine_batch(count);
    const std::vector<uint8_t> coefficients = make_binary_coefficients(count, kBatchSubgroupWidth);
    benchmark::DoNotOptimize(points.data());
    benchmark::DoNotOptimize(coefficients.data());

    for (auto _ : state) {
        const bool_t ok = batch_in_subgroup_g1_binary_assume_on_curve(
            points.data(),
            count,
            coefficients.data(),
            kBatchSubgroupWidth);
        benchmark::DoNotOptimize(ok);
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG1Affine.Breakdown.BatchInSubgroupBinaryW40");
}

static void BM_AvxG2_Breakdown_ProjectiveToAffine(benchmark::State& state)
{
    const G2Point a = make_g2(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        G2AffinePoint out = {};
        G2Point::to_affine(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2.Breakdown.ProjectiveToAffine");
}

static void BM_AvxG2_Breakdown_ProjectiveToAffineRepeated(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G2Point> points = make_g2_batch(count);
    std::vector<G2AffinePoint> out(count);
    benchmark::DoNotOptimize(points.data());

    for (auto _ : state) {
        for (size_t i = 0; i < count; ++i) {
            G2Point::to_affine(out[i], points[i]);
        }
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG2.Breakdown.ProjectiveToAffineRepeated");
}

static void BM_AvxG2_Breakdown_BatchProjectiveToAffine(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G2Point> points = make_g2_batch(count);
    std::vector<G2AffinePoint> out(count);
    std::vector<G2Point::Field_type> prefix_products(count);
    benchmark::DoNotOptimize(points.data());

    for (auto _ : state) {
        batch_to_affine_g2_nonzero(out.data(), points.data(), count, prefix_products);
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG2.Breakdown.BatchProjectiveToAffine");
}

static void BM_AvxG2Affine_Breakdown_AvxToScalar(benchmark::State& state)
{
    G2AffinePoint a = {};
    G2Point::to_affine(a, make_g2(42));
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const auto out = avx_curve_detail::curve_g2_affine_to_scalar({a.X, a.Y});
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine.Breakdown.AvxToScalar");
}

static void BM_AvxG2Affine_Breakdown_ScalarSerialize(benchmark::State& state)
{
    const blst::scalar::G2AffinePoint a = make_scalar_g2_affine(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        unsigned char out[192] = {};
        blst::scalar::G2AffinePoint::serialize(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine.Breakdown.ScalarSerialize");
}

static void BM_AvxG2Affine_Breakdown_ScalarCompress(benchmark::State& state)
{
    const blst::scalar::G2AffinePoint a = make_scalar_g2_affine(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        unsigned char out[96] = {};
        blst::scalar::G2AffinePoint::compress(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine.Breakdown.ScalarCompress");
}

static void BM_AvxG2Affine_Breakdown_BytesToScalar(benchmark::State& state)
{
    const blst::scalar::G2AffinePoint a = make_scalar_g2_affine(42);
    unsigned char serialized[192] = {};
    blst::scalar::G2AffinePoint::serialize(serialized, a);
    benchmark::DoNotOptimize(serialized);

    for (auto _ : state) {
        blst::scalar::G2AffinePoint out = {};
        const BLST_ERROR err = blst::scalar::G2AffinePoint::deserialize(out, serialized);
        benchmark::DoNotOptimize(err);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine.Breakdown.BytesToScalar");
}

static void BM_AvxG2Affine_Breakdown_CompressedBytesToScalar(benchmark::State& state)
{
    const blst::scalar::G2AffinePoint a = make_scalar_g2_affine(42);
    unsigned char compressed[96] = {};
    blst::scalar::G2AffinePoint::compress(compressed, a);
    benchmark::DoNotOptimize(compressed);

    for (auto _ : state) {
        blst::scalar::G2AffinePoint out = {};
        const BLST_ERROR err = blst::scalar::G2AffinePoint::uncompress(out, compressed);
        benchmark::DoNotOptimize(err);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine.Breakdown.CompressedBytesToScalar");
}

static void BM_AvxG2Affine_Breakdown_ScalarToAvx(benchmark::State& state)
{
    const blst::scalar::G2AffinePoint a = make_scalar_g2_affine(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const G2AffinePoint out(
            avx_curve_detail::scalar_fp2_to_ring(a.X),
            avx_curve_detail::scalar_fp2_to_ring(a.Y));
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine.Breakdown.ScalarToAvx");
}

static void BM_AvxG2Affine_Breakdown_ScalarOnCurve(benchmark::State& state)
{
    const blst::scalar::G2AffinePoint a = make_scalar_g2_affine(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = blst::scalar::G2AffinePoint::on_curve(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine.Breakdown.ScalarOnCurve");
}

static void BM_AvxG2Affine_Breakdown_ScalarInSubgroup(benchmark::State& state)
{
    const blst::scalar::G2AffinePoint a = make_scalar_g2_affine(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = blst::scalar::G2AffinePoint::in_subgroup(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxG2Affine.Breakdown.ScalarInSubgroup");
}

static void BM_AvxG2Affine_Breakdown_InSubgroupRepeated(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G2AffinePoint> points = make_g2_affine_batch(count);
    benchmark::DoNotOptimize(points.data());

    for (auto _ : state) {
        bool_t ok = true;
        for (size_t i = 0; i < count; ++i) {
            ok &= G2AffinePoint::in_subgroup(points[i]);
        }
        benchmark::DoNotOptimize(ok);
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG2Affine.Breakdown.InSubgroupRepeated");
}

static void BM_AvxG2Affine_Breakdown_BatchInSubgroupMod3W40(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G2AffinePoint> points = make_g2_affine_batch(count);
    const std::vector<uint8_t> coefficients = make_mod3_coefficients(count, kBatchSubgroupWidth);
    benchmark::DoNotOptimize(points.data());
    benchmark::DoNotOptimize(coefficients.data());

    for (auto _ : state) {
        const bool_t ok = batch_in_subgroup_g2_assume_on_curve(
            points.data(),
            count,
            coefficients.data(),
            kBatchSubgroupWidth);
        benchmark::DoNotOptimize(ok);
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG2Affine.Breakdown.BatchInSubgroupMod3W40");
}

static void BM_AvxG2Affine_Breakdown_BatchInSubgroupMod3Precompute2W40(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G2AffinePoint> points = make_g2_affine_batch(count);
    const std::vector<G2AffinePoint> doubles = make_g2_affine_doubles(points);
    const std::vector<uint8_t> coefficients = make_mod3_coefficients(count, kBatchSubgroupWidth);
    benchmark::DoNotOptimize(points.data());
    benchmark::DoNotOptimize(doubles.data());
    benchmark::DoNotOptimize(coefficients.data());

    for (auto _ : state) {
        const bool_t ok = batch_in_subgroup_g2_precomputed_twos_assume_on_curve(
            points.data(),
            doubles.data(),
            count,
            coefficients.data(),
            kBatchSubgroupWidth);
        benchmark::DoNotOptimize(ok);
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG2Affine.Breakdown.BatchInSubgroupMod3Precompute2W40");
}

static void BM_AvxG2Affine_Breakdown_BatchInSubgroupBinaryW40(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<G2AffinePoint> points = make_g2_affine_batch(count);
    const std::vector<uint8_t> coefficients = make_binary_coefficients(count, kBatchSubgroupWidth);
    benchmark::DoNotOptimize(points.data());
    benchmark::DoNotOptimize(coefficients.data());

    for (auto _ : state) {
        const bool_t ok = batch_in_subgroup_g2_binary_assume_on_curve(
            points.data(),
            count,
            coefficients.data(),
            kBatchSubgroupWidth);
        benchmark::DoNotOptimize(ok);
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxG2Affine.Breakdown.BatchInSubgroupBinaryW40");
}

static void BM_AvxGT_Breakdown_AvxToScalar(benchmark::State& state)
{
    const GT a = make_gt(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const blst::Fp12 out = avx_gt_detail::to_scalar(a.data);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT.Breakdown.AvxToScalar");
}

static void BM_AvxGT_Breakdown_ScalarSerialize(benchmark::State& state)
{
    const blst::scalar::GT a = make_scalar_gt(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        unsigned char out[12 * 48] = {};
        blst::scalar::GT::serialize(out, a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT.Breakdown.ScalarSerialize");
}

static void BM_AvxGT_Breakdown_BytesToScalar(benchmark::State& state)
{
    const blst::scalar::GT a = make_scalar_gt(42);
    unsigned char serialized[12 * 48] = {};
    blst::scalar::GT::serialize(serialized, a);
    benchmark::DoNotOptimize(serialized);

    for (auto _ : state) {
        blst::scalar::GT out = {};
        benchmark::DoNotOptimize(
            blst::scalar::GT::deserialize(out, serialized));
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT.Breakdown.BytesToScalar");
}

static void BM_AvxGT_Breakdown_ScalarToAvx(benchmark::State& state)
{
    const blst::scalar::GT a = make_scalar_gt(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const GT out(avx_gt_detail::from_scalar(a.data));
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT.Breakdown.ScalarToAvx");
}

static void BM_AvxGT_Breakdown_ScalarInGroup(benchmark::State& state)
{
    const blst::scalar::GT a = make_scalar_gt(42);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        const bool_t out = blst::scalar::GT::in_group(a);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxGT.Breakdown.ScalarInGroup");
}

static void BM_BlsFp_Sqrt(benchmark::State& state)
{
    const blst::scalar::Fp square = make_scalar_fp_square(42);
    benchmark::DoNotOptimize(square);

    for (auto _ : state) {
        blst::scalar::Fp out = {};
        const bool_t ok = blst::scalar::Fp::sqrt(out, square);
        benchmark::DoNotOptimize(ok);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("BlsFp_Sqrt");
}

static void BM_BlsFp2_Sqrt(benchmark::State& state)
{
    const blst::scalar::Fp2 square = make_scalar_fp2_square(42);
    benchmark::DoNotOptimize(square);

    for (auto _ : state) {
        blst::scalar::Fp2 out = {};
        const bool_t ok = blst::scalar::Fp2::sqrt(out, square);
        benchmark::DoNotOptimize(ok);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("BlsFp2_Sqrt");
}

static void BM_AvxPairing_MillerLoop(benchmark::State& state)
{
    G1Point g1 = make_g1(42);
    G1AffinePoint p = {};
    G1Point::to_affine(p, g1);
    const G2AffinePoint& q = G2AffinePoint::generator();
    benchmark::DoNotOptimize(p);
    benchmark::DoNotOptimize(q);

    for (auto _ : state) {
        MillerResult out = {};
        GT::miller_loop(out, &q, &p);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxPairing_MillerLoop");
}

static void BM_AvxPairing_FinalExp(benchmark::State& state)
{
    G1Point g1 = make_g1(42);
    G1AffinePoint p = {};
    G1Point::to_affine(p, g1);
    const G2AffinePoint& q = G2AffinePoint::generator();
    MillerResult miller = {};
    GT::miller_loop(miller, &q, &p);
    benchmark::DoNotOptimize(miller);

    for (auto _ : state) {
        GT out = {};
        GT::final_exp(out, miller);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxPairing_FinalExp");
}

static void BM_AvxPairing_Full(benchmark::State& state)
{
    G1Point g1 = make_g1(42);
    G1AffinePoint p = {};
    G1Point::to_affine(p, g1);
    const G2AffinePoint& q = G2AffinePoint::generator();
    benchmark::DoNotOptimize(p);
    benchmark::DoNotOptimize(q);

    for (auto _ : state) {
        MillerResult miller = {};
        GT out = {};
        GT::miller_loop(miller, &q, &p);
        GT::final_exp(out, miller);
        benchmark::DoNotOptimize(out);
    }

    state.SetLabel("AvxPairing_Full");
}

BENCHMARK(BM_AvxG1_Add);
BENCHMARK(BM_AvxG1_AddMixed);
BENCHMARK(BM_AvxG1_Dadd);
BENCHMARK(BM_AvxG1_DaddMixed);
BENCHMARK(BM_AvxG1_Double);
BENCHMARK(BM_AvxG1_Mult);
BENCHMARK(BM_AvxG1_MsmPippenger)->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)->Arg(8192)->Arg(16384)->Arg(32768);
BENCHMARK(BM_AvxG1_OnCurve);
BENCHMARK(BM_AvxG1_InSubgroup);
BENCHMARK(BM_AvxG1_Serialize);
BENCHMARK(BM_AvxG1Affine_OnCurve);
BENCHMARK(BM_AvxG1Affine_InSubgroup);
BENCHMARK(BM_AvxG1Affine_Serialize);
BENCHMARK(BM_AvxG1Affine_Deserialize);
BENCHMARK(BM_AvxG2_Add);
BENCHMARK(BM_AvxG2_AddMixed);
BENCHMARK(BM_AvxG2_Dadd);
BENCHMARK(BM_AvxG2_DaddMixed);
BENCHMARK(BM_AvxG2_Double);
BENCHMARK(BM_AvxG2_Mult);
BENCHMARK(BM_AvxG2_MsmPippenger)->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)->Arg(8192)->Arg(16384)->Arg(32768);
BENCHMARK(BM_AvxG2_OnCurve);
BENCHMARK(BM_AvxG2_InSubgroup);
BENCHMARK(BM_AvxG2_Serialize);
BENCHMARK(BM_AvxG2Affine_OnCurve);
BENCHMARK(BM_AvxG2Affine_InSubgroup);
BENCHMARK(BM_AvxG2Affine_Serialize);
BENCHMARK(BM_AvxG2Affine_Deserialize);
BENCHMARK(BM_AvxGT_Add);
BENCHMARK(BM_AvxGT_Dadd);
BENCHMARK(BM_AvxGT_Double);
BENCHMARK(BM_AvxGT_Mult);
BENCHMARK(BM_AvxGT_InGroup);
BENCHMARK(BM_AvxGT_Serialize);
BENCHMARK(BM_AvxGT_Deserialize);
BENCHMARK(BM_AvxG1_Breakdown_ProjectiveToAffine);
BENCHMARK(BM_AvxG1_Breakdown_ProjectiveToAffineRepeated)->Arg(1)->Arg(4)->Arg(16)->Arg(64)->Arg(256);
BENCHMARK(BM_AvxG1_Breakdown_BatchProjectiveToAffine)->Arg(1)->Arg(4)->Arg(16)->Arg(64)->Arg(256);
BENCHMARK(BM_AvxG1Affine_Breakdown_AvxToScalar);
BENCHMARK(BM_AvxG1Affine_Breakdown_ScalarSerialize);
BENCHMARK(BM_AvxG1Affine_Breakdown_ScalarCompress);
BENCHMARK(BM_AvxG1Affine_Breakdown_BytesToScalar);
BENCHMARK(BM_AvxG1Affine_Breakdown_CompressedBytesToScalar);
BENCHMARK(BM_AvxG1Affine_Breakdown_ScalarToAvx);
BENCHMARK(BM_AvxG1Affine_Breakdown_ScalarOnCurve);
BENCHMARK(BM_AvxG1Affine_Breakdown_ScalarInSubgroup);
BENCHMARK(BM_AvxG1Affine_Breakdown_InSubgroupRepeated)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_AvxG1Affine_Breakdown_BatchInSubgroupMod3W40)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_AvxG1Affine_Breakdown_BatchInSubgroupMod3Precompute2W40)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_AvxG1Affine_Breakdown_BatchInSubgroupBinaryW40)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_AvxG2_Breakdown_ProjectiveToAffine);
BENCHMARK(BM_AvxG2_Breakdown_ProjectiveToAffineRepeated)->Arg(1)->Arg(4)->Arg(16)->Arg(64)->Arg(256);
BENCHMARK(BM_AvxG2_Breakdown_BatchProjectiveToAffine)->Arg(1)->Arg(4)->Arg(16)->Arg(64)->Arg(256);
BENCHMARK(BM_AvxG2Affine_Breakdown_AvxToScalar);
BENCHMARK(BM_AvxG2Affine_Breakdown_ScalarSerialize);
BENCHMARK(BM_AvxG2Affine_Breakdown_ScalarCompress);
BENCHMARK(BM_AvxG2Affine_Breakdown_BytesToScalar);
BENCHMARK(BM_AvxG2Affine_Breakdown_CompressedBytesToScalar);
BENCHMARK(BM_AvxG2Affine_Breakdown_ScalarToAvx);
BENCHMARK(BM_AvxG2Affine_Breakdown_ScalarOnCurve);
BENCHMARK(BM_AvxG2Affine_Breakdown_ScalarInSubgroup);
BENCHMARK(BM_AvxG2Affine_Breakdown_InSubgroupRepeated)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_AvxG2Affine_Breakdown_BatchInSubgroupMod3W40)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_AvxG2Affine_Breakdown_BatchInSubgroupMod3Precompute2W40)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_AvxG2Affine_Breakdown_BatchInSubgroupBinaryW40)->Arg(64)->Arg(256)->Arg(1024);
BENCHMARK(BM_AvxGT_Breakdown_AvxToScalar);
BENCHMARK(BM_AvxGT_Breakdown_ScalarSerialize);
BENCHMARK(BM_AvxGT_Breakdown_BytesToScalar);
BENCHMARK(BM_AvxGT_Breakdown_ScalarToAvx);
BENCHMARK(BM_AvxGT_Breakdown_ScalarInGroup);
BENCHMARK(BM_BlsFp_Sqrt);
BENCHMARK(BM_BlsFp2_Sqrt);
BENCHMARK(BM_AvxPairing_MillerLoop);
BENCHMARK(BM_AvxPairing_FinalExp);
BENCHMARK(BM_AvxPairing_Full);

}  // namespace
