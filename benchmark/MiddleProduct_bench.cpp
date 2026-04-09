#include <benchmark/benchmark.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

#include "../test/MiddleProduct_reference.h"

namespace {

using blst::GT;
using blst::G1Point;
using blst::G2Point;
using blst::MiddleProductPairingMode;
using blst::MiddleProductScalarContext;
using blst::TimingStats;
using blst::test_support::middle_product_direct;
using blst::test_support::precompute_middle_product_b_direct;
using blst::test_support::update_direct;

G1Point make_g1_point(uint64_t scalar)
{
    G1Point out = {};
    blst::Fr s = {};
    blst::Fr::from_uint64(s, scalar);
    blst::Fr s_non_mont = {};
    blst::Fr::from_mont(s_non_mont, s);
    G1Point::mult(
        out,
        G1Point::generator(),
        reinterpret_cast<const blst::byte*>(&s_non_mont),
        256);
    return out;
}

G2Point make_g2_point(uint64_t scalar)
{
    G2Point out = {};
    blst::Fr s = {};
    blst::Fr::from_uint64(s, scalar);
    blst::Fr s_non_mont = {};
    blst::Fr::from_mont(s_non_mont, s);
    G2Point::mult(
        out,
        G2Point::generator(),
        reinterpret_cast<const blst::byte*>(&s_non_mont),
        256);
    return out;
}

std::vector<G1Point> make_g1_data(std::size_t n, uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist(1, 1024);

    std::vector<G1Point> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        if ((i % 17) == 4) {
            out.push_back(G1Point::infinity());
            continue;
        }
        out.push_back(make_g1_point(dist(rng)));
    }
    return out;
}

std::vector<G2Point> make_g2_data(std::size_t n, uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist(1, 1024);

    std::vector<G2Point> out;
    out.reserve(2 * n);
    for (std::size_t i = 0; i < (2 * n); ++i) {
        if ((i % 19) == 7) {
            out.push_back(G2Point::infinity());
            continue;
        }
        out.push_back(make_g2_point(dist(rng)));
    }
    return out;
}

template<class Context>
void bench_middle_product_direct_impl(
    benchmark::State& state,
    MiddleProductPairingMode pairing_mode,
    const char* offline_counter_name)
{
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    const std::vector<G1Point> a = make_g1_data(n, 0xD1EC7ULL);
    const std::vector<G2Point> b = make_g2_data(n, 0xBEEFULL);
    TimingStats pre_stats;
    const auto precomputed =
        precompute_middle_product_b_direct<Context>(b, n, &pre_stats, pairing_mode);
    std::vector<GT> out(n, GT::identity());

    double fill_ns = 0.0;
    double pointwise_ns = 0.0;
    double total_ns = 0.0;
    int64_t iterations = 0;

    for (auto _ : state) {
        TimingStats iter_stats;
        middle_product_direct<Context>(a, precomputed, std::span<GT>(out), &iter_stats);
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();

        fill_ns += static_cast<double>(iter_stats.get_ns("middle_product_direct.fill_ns"));
        pointwise_ns += static_cast<double>(iter_stats.get_ns("middle_product_direct.pointwise_ns"));
        total_ns += static_cast<double>(iter_stats.get_ns("middle_product_direct.total_ns"));
        iterations++;
    }

    state.counters["n"] = static_cast<double>(n);
    state.counters["pairs"] = benchmark::Counter(static_cast<double>(n) * static_cast<double>(n));
    state.counters["fill_ms"] = benchmark::Counter(fill_ns / 1e6 / iterations);
    state.counters["pointwise_ms"] = benchmark::Counter(pointwise_ns / 1e6 / iterations);
    state.counters["online_total_ms"] = benchmark::Counter(total_ns / 1e6 / iterations);
    state.counters["pairings/sec"] = benchmark::Counter(
        static_cast<double>(n) * static_cast<double>(n),
        benchmark::Counter::kIsIterationInvariantRate);
    state.counters[offline_counter_name] = benchmark::Counter(
        static_cast<double>(pre_stats.get_ns("middle_product_direct_precompute.total_ns")) / 1e6);
}

template<class Context>
void bench_middle_product_update_direct_impl(
    benchmark::State& state,
    MiddleProductPairingMode pairing_mode,
    const char* offline_counter_name)
{
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    const std::vector<G1Point> a = make_g1_data(n, 0xA77EULL);
    const std::vector<G2Point> b = make_g2_data(n, 0xA11CEULL);
    TimingStats pre_stats;
    const auto precomputed =
        precompute_middle_product_b_direct<Context>(b, n, &pre_stats, pairing_mode);
    std::vector<GT> out(n, GT::identity());

    double pointwise_ns = 0.0;
    double total_ns = 0.0;
    int64_t iterations = 0;

    for (auto _ : state) {
        TimingStats iter_stats;
        update_direct<Context>(a[0], precomputed, std::span<GT>(out), &iter_stats);
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();

        pointwise_ns +=
            static_cast<double>(iter_stats.get_ns("middle_product_direct_update.pointwise_ns"));
        total_ns += static_cast<double>(iter_stats.get_ns("middle_product_direct_update.total_ns"));
        iterations++;
    }

    state.counters["n"] = static_cast<double>(n);
    state.counters["pairs"] = benchmark::Counter(static_cast<double>(n));
    state.counters["pointwise_ms"] = benchmark::Counter(pointwise_ns / 1e6 / iterations);
    state.counters["online_total_ms"] = benchmark::Counter(total_ns / 1e6 / iterations);
    state.counters["pairings/sec"] = benchmark::Counter(
        static_cast<double>(n),
        benchmark::Counter::kIsIterationInvariantRate);
    state.counters[offline_counter_name] = benchmark::Counter(
        static_cast<double>(pre_stats.get_ns("middle_product_direct_precompute.total_ns")) / 1e6);
}

static void bench_middle_product_direct_scalar(benchmark::State& state)
{
    bench_middle_product_direct_impl<MiddleProductScalarContext>(
        state,
        MiddleProductPairingMode::Direct,
        "offline_scalar_pre_ms");
}

static void bench_middle_product_direct_scalar_lines(benchmark::State& state)
{
    bench_middle_product_direct_impl<MiddleProductScalarContext>(
        state,
        MiddleProductPairingMode::PrecomputedLines,
        "offline_scalar_lines_pre_ms");
}

static void bench_middle_product_update_direct_scalar(benchmark::State& state)
{
    bench_middle_product_update_direct_impl<MiddleProductScalarContext>(
        state,
        MiddleProductPairingMode::Direct,
        "offline_scalar_pre_ms");
}

static void bench_middle_product_update_direct_scalar_lines(benchmark::State& state)
{
    bench_middle_product_update_direct_impl<MiddleProductScalarContext>(
        state,
        MiddleProductPairingMode::PrecomputedLines,
        "offline_scalar_lines_pre_ms");
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
static void bench_middle_product_direct_avx(benchmark::State& state)
{
    bench_middle_product_direct_impl<blst::MiddleProductAvxContext>(
        state,
        MiddleProductPairingMode::Direct,
        "offline_avx_pre_ms");
}

static void bench_middle_product_direct_avx_lines(benchmark::State& state)
{
    bench_middle_product_direct_impl<blst::MiddleProductAvxContext>(
        state,
        MiddleProductPairingMode::PrecomputedLines,
        "offline_avx_lines_pre_ms");
}

static void bench_middle_product_update_direct_avx(benchmark::State& state)
{
    bench_middle_product_update_direct_impl<blst::MiddleProductAvxContext>(
        state,
        MiddleProductPairingMode::Direct,
        "offline_avx_pre_ms");
}

static void bench_middle_product_update_direct_avx_lines(benchmark::State& state)
{
    bench_middle_product_update_direct_impl<blst::MiddleProductAvxContext>(
        state,
        MiddleProductPairingMode::PrecomputedLines,
        "offline_avx_lines_pre_ms");
}
#endif

} // namespace

#define DIRECT_MIDDLE_PROD_SIZES \
->DenseRange(64, 256, 32)

BENCHMARK(bench_middle_product_direct_scalar)
DIRECT_MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_direct_scalar_lines)
DIRECT_MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_update_direct_scalar)
DIRECT_MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_update_direct_scalar_lines)
DIRECT_MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
BENCHMARK(bench_middle_product_direct_avx)
DIRECT_MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_direct_avx_lines)
DIRECT_MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_update_direct_avx)
DIRECT_MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_update_direct_avx_lines)
DIRECT_MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
#endif
