#include <benchmark/benchmark.h>

#include <array>
#include <cstddef>

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

#define CATBLST_GROUP_MULT_BENCHMARK(function_name) \
    BENCHMARK(function_name)->Args({128, 0})->Args({256, 0})->Args({256, 1})->Args({256, 2})

CATBLST_GROUP_MULT_BENCHMARK(BM_G1MultCt);
CATBLST_GROUP_MULT_BENCHMARK(BM_G1MultVartime);
CATBLST_GROUP_MULT_BENCHMARK(BM_G2MultCt);
CATBLST_GROUP_MULT_BENCHMARK(BM_G2MultVartime);
CATBLST_GROUP_MULT_BENCHMARK(BM_GTMultCt);
CATBLST_GROUP_MULT_BENCHMARK(BM_GTMultVartime);
#undef CATBLST_GROUP_MULT_BENCHMARK

}  // namespace
