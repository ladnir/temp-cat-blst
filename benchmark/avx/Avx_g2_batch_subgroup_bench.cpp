#include <benchmark/benchmark.h>

#include "bindings/blst.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace {

struct BenchmarkRandomSource {
    std::uint64_t counter = 1;
};

bool fill_benchmark_random(void* context, blst_byte* out, size_t out_len)
{
    if (context == nullptr || out == nullptr || out_len != 32) {
        return false;
    }
    auto& source = *static_cast<BenchmarkRandomSource*>(context);
    std::uint64_t value = source.counter++;
    for (size_t i = 0; i < out_len; ++i) {
        value += 0x9e3779b97f4a7c15ULL;
        std::uint64_t mixed = value;
        mixed = (mixed ^ (mixed >> 30)) * 0xbf58476d1ce4e5b9ULL;
        mixed = (mixed ^ (mixed >> 27)) * 0x94d049bb133111ebULL;
        mixed ^= mixed >> 31;
        out[i] = static_cast<blst_byte>(mixed);
    }
    return true;
}

std::vector<blst_avx_p2> make_distinct_g2_points(size_t count)
{
    std::vector<blst_avx_p2> points(count);
    if (count == 0) {
        return points;
    }
    const blst_avx_p2* generator = blst_avx_p2_generator();
    BLST_ASSERT(generator != nullptr);
    points[0] = *generator;
    for (size_t i = 1; i < count; ++i) {
        const bool added = blst_avx_p2_add(
            &points[i], &points[i - 1], generator);
        BLST_ASSERT(added);
    }
    return points;
}

std::vector<const blst_avx_p2*> point_pointers(
    const std::vector<blst_avx_p2>& points)
{
    std::vector<const blst_avx_p2*> out(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        out[i] = &points[i];
    }
    return out;
}

using EncodedG2 = std::array<blst_byte, 96>;

std::vector<EncodedG2> encode_points(
    const std::vector<blst_avx_p2>& points)
{
    std::vector<EncodedG2> out(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        const bool encoded = blst_avx_p2_compress(out[i].data(), &points[i]);
        BLST_ASSERT(encoded);
    }
    return out;
}

void BM_G2IndividualSubgroupChecks(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<blst_avx_p2> points = make_distinct_g2_points(count);
    for (auto _ : state) {
        bool accepted = true;
        for (const blst_avx_p2& point : points) {
            accepted &= blst_avx_p2_in_g2(&point);
        }
        benchmark::DoNotOptimize(accepted);
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
}

void BM_G2BatchSubgroupCheck(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<blst_avx_p2> points = make_distinct_g2_points(count);
    const std::vector<const blst_avx_p2*> pointers = point_pointers(points);
    std::vector<blst_byte> workspace(
        blst_avx_p2s_in_g2_workspace_sizeof());
    BenchmarkRandomSource random;
    for (auto _ : state) {
        bool accepted = blst_avx_p2s_in_g2(
            pointers.data(), pointers.size(), fill_benchmark_random, &random,
            workspace.data(), workspace.size());
        benchmark::DoNotOptimize(accepted);
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
}

void BM_G2DecodeAndIndividualSubgroupChecks(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<EncodedG2> encodings =
        encode_points(make_distinct_g2_points(count));
    std::vector<blst_avx_p2> decoded(count);
    for (auto _ : state) {
        bool accepted = true;
        for (size_t i = 0; i < count; ++i) {
            accepted &= blst_avx_p2_uncompress(
                &decoded[i], encodings[i].data()) == BLST_SUCCESS;
            accepted &= blst_avx_p2_in_g2(&decoded[i]);
        }
        benchmark::DoNotOptimize(accepted);
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
}

void BM_G2DecodeAndBatchSubgroupCheck(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<EncodedG2> encodings =
        encode_points(make_distinct_g2_points(count));
    std::vector<blst_avx_p2> decoded(count);
    const std::vector<const blst_avx_p2*> pointers = point_pointers(decoded);
    std::vector<blst_byte> workspace(
        blst_avx_p2s_in_g2_workspace_sizeof());
    BenchmarkRandomSource random;
    for (auto _ : state) {
        bool accepted = true;
        for (size_t i = 0; i < count; ++i) {
            accepted &= blst_avx_p2_uncompress(
                &decoded[i], encodings[i].data()) == BLST_SUCCESS;
        }
        accepted &= blst_avx_p2s_in_g2(
            pointers.data(), pointers.size(), fill_benchmark_random, &random,
            workspace.data(), workspace.size());
        benchmark::DoNotOptimize(accepted);
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * static_cast<int64_t>(count));
}

void apply_g2_batch_sizes(benchmark::internal::Benchmark* registration)
{
    registration->Arg(2048)->Arg(18432)->Unit(benchmark::kMillisecond)->MinTime(1.0);
}

BENCHMARK(BM_G2IndividualSubgroupChecks)->Apply(apply_g2_batch_sizes);
BENCHMARK(BM_G2BatchSubgroupCheck)->Apply(apply_g2_batch_sizes);
BENCHMARK(BM_G2DecodeAndIndividualSubgroupChecks)->Apply(apply_g2_batch_sizes);
BENCHMARK(BM_G2DecodeAndBatchSubgroupCheck)->Apply(apply_g2_batch_sizes);

}  // namespace
