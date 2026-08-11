#include <benchmark/benchmark.h>

#include "blst/avx/curve/multi_scalar.h"
#include "blst/scalar/curve/multi_scalar.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace {

using ScalarBytes = std::array<blst::byte, 32>;

ScalarBytes scalar_bytes(uint64_t seed)
{
    ScalarBytes out = {};
    uint64_t x = seed * 0x9e3779b97f4a7c15ULL + 0xbf58476d1ce4e5b9ULL;
    for (size_t i = 0; i < out.size(); ++i) {
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        out[i] = static_cast<blst::byte>(x >> ((i & 7u) * 8u));
    }
    out[31] &= 0x0f;
    return out;
}

std::vector<blst::byte> scalar_batch(size_t count)
{
    std::vector<blst::byte> scalars(count * 32);
    for (size_t i = 0; i < count; ++i) {
        const ScalarBytes scalar = scalar_bytes(1000 + i);
        std::copy(scalar.begin(), scalar.end(), scalars.begin() + static_cast<std::ptrdiff_t>(i * 32));
    }
    return scalars;
}

blst::avx::G1AffinePoint make_g1_affine(uint64_t seed)
{
    const ScalarBytes scalar = scalar_bytes(seed);
    blst::avx::G1Point point = {};
    blst::avx::G1Point::mult(point, blst::avx::G1Point::generator(), scalar.data(), 256);
    blst::avx::G1AffinePoint out = {};
    blst::avx::G1Point::to_affine(out, point);
    return out;
}

blst::avx::G2AffinePoint make_g2_affine(uint64_t seed)
{
    const ScalarBytes scalar = scalar_bytes(seed);
    blst::avx::G2Point point = {};
    blst::avx::G2Point::mult(point, blst::avx::G2Point::generator(), scalar.data(), 256);
    blst::avx::G2AffinePoint out = {};
    blst::avx::G2Point::to_affine(out, point);
    return out;
}

std::vector<blst::avx::GT> make_avx_gt_batch(size_t count)
{
    std::vector<blst::avx::GT> out(count);
    for (size_t i = 0; i < count; ++i) {
        const blst::avx::G1AffinePoint g1 = make_g1_affine(42 + i);
        const blst::avx::G2AffinePoint g2 = make_g2_affine(4242 + i);
        blst::avx::MillerResult miller = {};
        blst::avx::GT::miller_loop(miller, &g2, &g1);
        blst::avx::GT::final_exp(out[i], miller);
    }
    return out;
}

blst::scalar::GT scalar_from_avx(const blst::avx::GT& point)
{
    unsigned char bytes[12 * 48] = {};
    blst::avx::GT::serialize(bytes, point);
    blst::scalar::GT out = {};
    [[maybe_unused]] const bool decoded =
        blst::scalar::GT::deserialize_unchecked(out, bytes);
    BLST_ASSERT(decoded);
    return out;
}

std::vector<blst::scalar::GT> make_scalar_gt_batch(const std::vector<blst::avx::GT>& points)
{
    std::vector<blst::scalar::GT> out(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        out[i] = scalar_from_avx(points[i]);
    }
    return out;
}

static void BM_GTScalarGenericMsm(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<blst::avx::GT> avx_points = make_avx_gt_batch(count);
    const std::vector<blst::scalar::GT> points = make_scalar_gt_batch(avx_points);
    const std::vector<blst::byte> scalars = scalar_batch(count);
    const blst::scalar::GT* point_ptrs[] = {points.data(), nullptr};
    const blst::byte* scalar_ptrs[] = {scalars.data(), nullptr};
    std::vector<blst::scalar::GT> buckets(
        blst::scalar::msm_mult_pippenger_scratch_sizeof_<blst::scalar::GT>(count) / sizeof(blst::scalar::GT));

    benchmark::DoNotOptimize(points.data());
    benchmark::DoNotOptimize(scalars.data());
    for (auto _ : state) {
        blst::scalar::GT out = {};
        blst::scalar::msm_mult_pippenger<blst::scalar::GT, blst::scalar::GT, blst::scalar::GT>(
            out,
            point_ptrs,
            count,
            scalar_ptrs,
            256,
            buckets.data(),
            blst::scalar::pippenger_window_size(count));
        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("GTScalarGenericMsm");
}

static void BM_AvxGTMsmGLV4(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<blst::avx::GT> points = make_avx_gt_batch(count);
    const std::vector<blst::byte> scalars = scalar_batch(count);
    std::vector<blst::avx::AvxGtPippengerBucket> buckets(
        blst::avx::avx_gt_mult_pippenger_scratch_sizeof(count) / sizeof(blst::avx::AvxGtPippengerBucket));
    std::vector<uint8_t> bucket_used(blst::avx::avx_gt_mult_pippenger_bucket_used_scratch_sizeof(count));
    std::vector<blst::avx::GT> split_points(
        blst::avx::avx_gt_mult_pippenger_glv4_points_scratch_sizeof(count) / sizeof(blst::avx::GT));
    std::vector<blst::byte> split_scalars(blst::avx::avx_gt_mult_pippenger_glv4_scalars_scratch_sizeof(count));

    benchmark::DoNotOptimize(points.data());
    benchmark::DoNotOptimize(scalars.data());
    for (auto _ : state) {
        blst::avx::GT out = {};
        blst::avx::avx_gt_mult_pippenger_glv4(
            out,
            points.data(),
            scalars.data(),
            count,
            buckets.data(),
            0,
            bucket_used.data(),
            split_points.data(),
            split_scalars.data());
        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("AvxGTMsmGLV4");
}

static void BM_GTNaiveMsm(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<blst::avx::GT> points = make_avx_gt_batch(count);
    const std::vector<blst::byte> scalars = scalar_batch(count);

    benchmark::DoNotOptimize(points.data());
    benchmark::DoNotOptimize(scalars.data());
    for (auto _ : state) {
        blst::avx::GT out = blst::avx::GT::identity();
        for (size_t i = 0; i < count; ++i) {
            blst::avx::GT term = {};
            blst::avx::GT::mult(term, points[i], scalars.data() + i * 32, 256);
            blst::avx::GT::add(out, out, term);
        }
        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
    state.SetLabel("GTNaiveMsm");
}

BENCHMARK(BM_GTScalarGenericMsm)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)->Arg(8192)->Arg(16384);
BENCHMARK(BM_AvxGTMsmGLV4)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)->Arg(8192)->Arg(16384);
BENCHMARK(BM_GTNaiveMsm)->Arg(16)->Arg(32)->Arg(64);

}  // namespace
