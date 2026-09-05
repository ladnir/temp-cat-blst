#include <benchmark/benchmark.h>

#include "bindings/blst.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <thread>
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

std::vector<blst_avx_p1> make_distinct_g1_points(size_t count)
{
    std::vector<blst_avx_p1> points(count);
    if (count == 0) {
        return points;
    }
    const blst_avx_p1* generator = blst_avx_p1_generator();
    BLST_ASSERT(generator != nullptr);
    points[0] = *generator;
    for (size_t i = 1; i < count; ++i) {
        const bool added = blst_avx_p1_add(
            &points[i], &points[i - 1], generator);
        BLST_ASSERT(added);
    }
    return points;
}

std::vector<const blst_avx_p1*> point_pointers(
    const std::vector<blst_avx_p1>& points)
{
    std::vector<const blst_avx_p1*> out(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        out[i] = &points[i];
    }
    return out;
}

void BM_G1FriendlyBatchSubgroupCheck(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const std::vector<blst_avx_p1> points = make_distinct_g1_points(count);
    const std::vector<const blst_avx_p1*> pointers = point_pointers(points);
    BenchmarkRandomSource random;
    for (auto _ : state) {
        bool accepted = blst_avx_p1s_in_g1(
            pointers.data(), pointers.size(), fill_benchmark_random, &random);
        benchmark::DoNotOptimize(accepted);
        benchmark::ClobberMemory();
    }
}

void BM_G1ChallengeRangeBatchSubgroupCheck(benchmark::State& state)
{
    const size_t count = static_cast<size_t>(state.range(0));
    const size_t worker_count = static_cast<size_t>(state.range(1));
    const std::vector<blst_avx_p1> points = make_distinct_g1_points(count);
    const std::vector<const blst_avx_p1*> pointers = point_pointers(points);
    const size_t challenge_bytes =
        blst_avx_p1_batch_subgroup_challenge_sizeof(count);
    std::vector<blst_byte> challenge(challenge_bytes);
    BenchmarkRandomSource random;

    for (auto _ : state) {
        bool initialized = blst_avx_p1_batch_subgroup_challenge_init(
            challenge.data(), challenge.size(), count,
            fill_benchmark_random, &random);
        std::atomic<bool> accepted{initialized};
        std::vector<std::thread> workers;
        workers.reserve(worker_count);
        const size_t rounds = challenge_bytes / count;
        for (size_t worker = 0; worker < worker_count; ++worker) {
            const size_t begin = rounds * worker / worker_count;
            const size_t end = rounds * (worker + 1) / worker_count;
            workers.emplace_back([&, begin, end] {
                if (!blst_avx_p1_batch_subgroup_challenge_check_range(
                        pointers.data(), count,
                        challenge.data(), challenge.size(), begin, end)) {
                    accepted.store(false, std::memory_order_relaxed);
                }
            });
        }
        for (std::thread& worker : workers) {
            worker.join();
        }
        benchmark::DoNotOptimize(accepted.load(std::memory_order_relaxed));
        benchmark::ClobberMemory();
    }
}

void apply_batch_sizes(benchmark::internal::Benchmark* registration)
{
    registration->Arg(64)->Arg(256)->Arg(1024)
        ->Unit(benchmark::kMillisecond)->MinTime(0.2);
}

void apply_parallel_sizes(benchmark::internal::Benchmark* registration)
{
    for (const int64_t count : {64, 256, 1024}) {
        for (const int64_t workers : {1, 2, 4, 8, 16}) {
            registration->Args({count, workers});
        }
    }
    registration->Unit(benchmark::kMillisecond)->MinTime(0.2)->UseRealTime();
}

BENCHMARK(BM_G1FriendlyBatchSubgroupCheck)->Apply(apply_batch_sizes);
BENCHMARK(BM_G1ChallengeRangeBatchSubgroupCheck)->Apply(apply_parallel_sizes);

}  // namespace
