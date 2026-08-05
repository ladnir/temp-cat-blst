#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "bindings/blst.h"

namespace {

std::array<blst_byte, 32> benchmark_scalar(size_t nbits, size_t variant = 0)
{
    std::array<blst_byte, 32> scalar = {};
    for (size_t i = 0; i < scalar.size(); ++i) {
        scalar[i] = static_cast<blst_byte>(
            (0x6dU + static_cast<unsigned>(variant * 0x53U)) ^
            static_cast<unsigned>(i * (37U + 2U * variant)));
    }
    const size_t nbytes = (nbits + 7) / 8;
    for (size_t i = nbytes; i < scalar.size(); ++i) {
        scalar[i] = 0;
    }
    if ((nbits & 7U) != 0) {
        scalar[nbytes - 1] &= static_cast<blst_byte>((1U << (nbits & 7U)) - 1U);
    }
    // Keep full-width benchmark scalars canonical so the optimized G1/G2
    // paths exercise their intended Fr-scalar fast path.
    scalar[31] &= 0x3fU;
    return scalar;
}

template<class Point, const Point* (*Generator)(),
    void (*Multiply)(Point*, const Point*, const blst_byte*, size_t)>
void benchmark_point_mult(benchmark::State& state)
{
    const size_t nbits = static_cast<size_t>(state.range(0));
    const auto scalar = benchmark_scalar(nbits, static_cast<size_t>(state.range(1)));
    const Point base = *Generator();
    Point out = {};

    for (auto _ : state) {
        Multiply(&out, &base, scalar.data(), nbits);
        benchmark::DoNotOptimize(out);
    }
    state.SetItemsProcessed(state.iterations());
}

void BM_G1MultCt(benchmark::State& state)
{
    benchmark_point_mult<blst_p1, blst_p1_generator, blst_p1_mult>(state);
}

void BM_G1MultVartime(benchmark::State& state)
{
    benchmark_point_mult<blst_p1, blst_p1_generator, blst_p1_mult_vartime>(state);
}

void BM_G2MultCt(benchmark::State& state)
{
    benchmark_point_mult<blst_p2, blst_p2_generator, blst_p2_mult>(state);
}

void BM_G2MultVartime(benchmark::State& state)
{
    benchmark_point_mult<blst_p2, blst_p2_generator, blst_p2_mult_vartime>(state);
}

template<void (*Multiply)(blst_fp12*, const blst_fp12*, const blst_byte*, size_t)>
void benchmark_gt_mult(benchmark::State& state)
{
    const size_t nbits = static_cast<size_t>(state.range(0));
    const auto scalar = benchmark_scalar(nbits, static_cast<size_t>(state.range(1)));
    const blst_fp12 base = *blst_gt_generator();
    blst_fp12 out = {};

    for (auto _ : state) {
        Multiply(&out, &base, scalar.data(), nbits);
        benchmark::DoNotOptimize(out);
    }
    state.SetItemsProcessed(state.iterations());
}

void BM_GTMultCt(benchmark::State& state)
{
    benchmark_gt_mult<blst_gt_mult>(state);
}

void BM_GTMultVartime(benchmark::State& state)
{
    benchmark_gt_mult<blst_gt_mult_vartime>(state);
}

template<class Point, class Affine,
    const Point* (*Generator)(),
    void (*ToAffine)(Affine*, const Point*),
    size_t (*ScratchSize)(size_t),
    void (*Pippenger)(Point*, const Affine* const[], size_t,
        const blst_byte* const[], size_t, blst_limb_t*)>
void benchmark_point_msm(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const Point generator = *Generator();
    Affine generator_affine = {};
    ToAffine(&generator_affine, &generator);

    std::vector<Affine> points(count, generator_affine);
    std::vector<blst_byte> scalars(count * 32);
    for (size_t i = 0; i < count; ++i) {
        const auto scalar = benchmark_scalar(256, i);
        std::copy(scalar.begin(), scalar.end(), scalars.begin() + i * 32);
    }

    const Affine* point_ptrs[2] = {points.data(), nullptr};
    const blst_byte* scalar_ptrs[2] = {scalars.data(), nullptr};
    const size_t scratch_bytes = ScratchSize(count);
    std::vector<blst_limb_t> scratch(
        (scratch_bytes + sizeof(blst_limb_t) - 1) / sizeof(blst_limb_t));
    Point out = {};

    for (auto _ : state) {
        Pippenger(&out, point_ptrs, count, scalar_ptrs, 256, scratch.data());
        benchmark::DoNotOptimize(out);
    }
    state.SetItemsProcessed(
        static_cast<int64_t>(state.iterations() * count));
}

void BM_G1Msm(benchmark::State& state)
{
    benchmark_point_msm<blst_p1, blst_p1_affine,
        blst_p1_generator, blst_p1_to_affine,
        blst_p1s_mult_pippenger_scratch_sizeof,
        blst_p1s_mult_pippenger>(state);
}

void BM_G2Msm(benchmark::State& state)
{
    benchmark_point_msm<blst_p2, blst_p2_affine,
        blst_p2_generator, blst_p2_to_affine,
        blst_p2s_mult_pippenger_scratch_sizeof,
        blst_p2s_mult_pippenger>(state);
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
void BM_GTMsm(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    blst_avx_gt generator = {};
    if (!blst_avx_fp12_from_scalar(&generator, blst_gt_generator())) {
        state.SkipWithError("AVX GT conversion is unavailable");
        return;
    }

    std::vector<blst_avx_gt> points(count, generator);
    std::vector<blst_byte> scalars(count * 32);
    for (size_t i = 0; i < count; ++i) {
        const auto scalar = benchmark_scalar(256, i);
        std::copy(scalar.begin(), scalar.end(), scalars.begin() + i * 32);
    }

    const size_t scratch_bytes =
        blst_avx_gt_mult_pippenger_scratch_sizeof(count);
    const size_t bucket_used_bytes =
        blst_avx_gt_mult_pippenger_bucket_used_scratch_sizeof(count);
    std::vector<blst_avx_gt> scratch(
        (scratch_bytes + sizeof(blst_avx_gt) - 1) / sizeof(blst_avx_gt));
    std::vector<blst_byte> bucket_used(bucket_used_bytes);
    blst_avx_gt out = {};

    for (auto _ : state) {
        bool ok = blst_avx_gt_mult_pippenger(
            &out, points.data(), scalars.data(), count,
            scratch.data(), scratch.size() * sizeof(blst_avx_gt),
            bucket_used.data(), bucket_used.size());
        benchmark::DoNotOptimize(ok);
        benchmark::DoNotOptimize(out);
    }
    state.SetItemsProcessed(
        static_cast<int64_t>(state.iterations() * count));
}
#endif

#define CATBLST_GROUP_MULT_BENCHMARK(function_name) \
    BENCHMARK(function_name)->Args({128, 0})->Args({256, 0})->Args({256, 1})->Args({256, 2})

CATBLST_GROUP_MULT_BENCHMARK(BM_G1MultCt);
CATBLST_GROUP_MULT_BENCHMARK(BM_G1MultVartime);
CATBLST_GROUP_MULT_BENCHMARK(BM_G2MultCt);
CATBLST_GROUP_MULT_BENCHMARK(BM_G2MultVartime);
CATBLST_GROUP_MULT_BENCHMARK(BM_GTMultCt);
CATBLST_GROUP_MULT_BENCHMARK(BM_GTMultVartime);
#undef CATBLST_GROUP_MULT_BENCHMARK

#define CATBLST_GROUP_MSM_BENCHMARK(function_name) \
    BENCHMARK(function_name)->Arg(1)->Arg(8)->Arg(64)->Arg(256)

CATBLST_GROUP_MSM_BENCHMARK(BM_G1Msm);
CATBLST_GROUP_MSM_BENCHMARK(BM_G2Msm);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
CATBLST_GROUP_MSM_BENCHMARK(BM_GTMsm);
#endif
#undef CATBLST_GROUP_MSM_BENCHMARK

}  // namespace
