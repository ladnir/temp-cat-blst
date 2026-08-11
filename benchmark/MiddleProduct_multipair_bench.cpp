#include <benchmark/benchmark.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <type_traits>
#include <vector>

#include "blst/fft/middle_product.h"

namespace {

using blst::GT;
using blst::G1Point;
using blst::G2Point;
using blst::MiddleProductScalarContext;
using blst::TimingStats;

template<class Context>
auto precompute_middle_product_b_multipair_for_context(
    std::span<const G2Point> b,
    std::size_t n,
    TimingStats* timings)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        return blst::precompute_middle_product_b_multipair_scalar(b, n, timings);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else {
        return blst::precompute_middle_product_b_multipair_avx(b, n, timings);
#endif
    }
}

template<class Context, class PrecomputedB>
void middle_product_multipair_for_context(
    std::span<const G1Point> a,
    const PrecomputedB& precomputed_b,
    std::span<GT> out,
    TimingStats* timings)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        blst::middle_product_multipair_scalar(a, precomputed_b, out, timings);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else {
        blst::middle_product_multipair_avx(a, precomputed_b, out, timings);
#endif
    }
}

template<class Context, class PrecomputedB>
void update_multipair_for_context(
    const G1Point& a0,
    const PrecomputedB& precomputed_b,
    std::span<GT> out,
    TimingStats* timings)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        blst::update_multipair_scalar(a0, precomputed_b, out, timings);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else {
        blst::update_multipair_avx(a0, precomputed_b, out, timings);
#endif
    }
}

template<class Context>
void middle_product_one_shot_for_context(
    std::span<const G1Point> a,
    std::span<const G2Point> b,
    std::span<GT> out,
    blst::MiddleProductAlgorithm algorithm)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        blst::middle_product_scalar(a, b, out, algorithm);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else {
        blst::middle_product_avx(a, b, out, algorithm);
#endif
    }
}

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
void bench_middle_product_multipair_impl(
    benchmark::State& state,
    const char* offline_counter_name)
{
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    const std::vector<G1Point> a = make_g1_data(n, 0xD1EC7ULL);
    const std::vector<G2Point> b = make_g2_data(n, 0xBEEFULL);
    TimingStats pre_stats;
    const auto precomputed =
        precompute_middle_product_b_multipair_for_context<Context>(b, n, &pre_stats);
    std::vector<GT> out(n, GT::identity());

    double pointwise_ns = 0.0;
    double total_ns = 0.0;
    int64_t iterations = 0;

    for (auto _ : state) {
        TimingStats iter_stats;
        middle_product_multipair_for_context<Context>(
            a,
            precomputed,
            std::span<GT>(out),
            &iter_stats);
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();

        pointwise_ns += static_cast<double>(iter_stats.get_ns("middle_product_multipair.pointwise_ns"));
        total_ns += static_cast<double>(iter_stats.get_ns("middle_product_multipair.total_ns"));
        iterations++;
    }

    state.counters["n"] = static_cast<double>(n);
    state.counters["pairs"] = benchmark::Counter(static_cast<double>(n) * static_cast<double>(n));
    state.counters["pointwise_ms"] = benchmark::Counter(pointwise_ns / 1e6 / iterations);
    state.counters["online_total_ms"] = benchmark::Counter(total_ns / 1e6 / iterations);
    state.counters["pairings/sec"] = benchmark::Counter(
        static_cast<double>(n) * static_cast<double>(n),
        benchmark::Counter::kIsIterationInvariantRate);
    state.counters[offline_counter_name] = benchmark::Counter(
        static_cast<double>(pre_stats.get_ns("middle_product_multipair_precompute.total_ns")) / 1e6);
}

template<class Context>
void bench_middle_product_update_multipair_impl(
    benchmark::State& state,
    const char* offline_counter_name)
{
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    const std::vector<G1Point> a = make_g1_data(n, 0xA77EULL);
    const std::vector<G2Point> b = make_g2_data(n, 0xA11CEULL);
    TimingStats pre_stats;
    const auto precomputed =
        precompute_middle_product_b_multipair_for_context<Context>(b, n, &pre_stats);
    std::vector<GT> out(n, GT::identity());

    double pointwise_ns = 0.0;
    double total_ns = 0.0;
    int64_t iterations = 0;

    for (auto _ : state) {
        TimingStats iter_stats;
        update_multipair_for_context<Context>(a[0], precomputed, std::span<GT>(out), &iter_stats);
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();

        pointwise_ns +=
            static_cast<double>(iter_stats.get_ns("middle_product_multipair_update.pointwise_ns"));
        total_ns +=
            static_cast<double>(iter_stats.get_ns("middle_product_multipair_update.total_ns"));
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
        static_cast<double>(pre_stats.get_ns("middle_product_multipair_precompute.total_ns")) / 1e6);
}

template<class Context, blst::MiddleProductAlgorithm Algorithm>
void bench_middle_product_one_shot_impl(benchmark::State& state)
{
    const std::size_t n = static_cast<std::size_t>(state.range(0));
    const std::vector<G1Point> a = make_g1_data(n, 0xA0170ULL);
    const std::vector<G2Point> b = make_g2_data(n, 0xA1170ULL);
    std::vector<GT> out(n, GT::identity());

    for (auto _ : state) {
        middle_product_one_shot_for_context<Context>(a, b, out, Algorithm);
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }

    state.counters["n"] = static_cast<double>(n);
}

static void bench_middle_product_multipair_scalar(benchmark::State& state)
{
    bench_middle_product_multipair_impl<MiddleProductScalarContext>(
        state,
        "offline_scalar_multipair_pre_ms");
}

static void bench_middle_product_update_multipair_scalar(benchmark::State& state)
{
    bench_middle_product_update_multipair_impl<MiddleProductScalarContext>(
        state,
        "offline_scalar_multipair_pre_ms");
}

static void bench_middle_product_one_shot_scalar_fft(benchmark::State& state)
{
    bench_middle_product_one_shot_impl<
        MiddleProductScalarContext,
        blst::MiddleProductAlgorithm::Fft>(state);
}

static void bench_middle_product_one_shot_scalar_multipair(benchmark::State& state)
{
    bench_middle_product_one_shot_impl<
        MiddleProductScalarContext,
        blst::MiddleProductAlgorithm::Multipair>(state);
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
static void bench_middle_product_multipair_avx(benchmark::State& state)
{
    bench_middle_product_multipair_impl<blst::MiddleProductAvxContext>(
        state,
        "offline_avx_multipair_pre_ms");
}

static void bench_middle_product_update_multipair_avx(benchmark::State& state)
{
    bench_middle_product_update_multipair_impl<blst::MiddleProductAvxContext>(
        state,
        "offline_avx_multipair_pre_ms");
}

static void bench_middle_product_one_shot_avx_fft(benchmark::State& state)
{
    bench_middle_product_one_shot_impl<
        blst::MiddleProductAvxContext,
        blst::MiddleProductAlgorithm::Fft>(state);
}

static void bench_middle_product_one_shot_avx_multipair(benchmark::State& state)
{
    bench_middle_product_one_shot_impl<
        blst::MiddleProductAvxContext,
        blst::MiddleProductAlgorithm::Multipair>(state);
}
#endif

} // namespace

#define MIDDLE_PRODUCT_MULTIPAIR_SIZES \
->Arg(2)->Arg(4)->Arg(8)->Arg(12)->Arg(16)->Arg(24)->Arg(32)->Arg(48)->Arg(64) \
->DenseRange(96, 256, 32)

#define MIDDLE_PRODUCT_AUTO_CROSSOVER_SIZES \
->Arg(2)->Arg(4)->Arg(8)->Arg(12)->Arg(16)->Arg(20)->Arg(24)->Arg(28)->Arg(30) \
->Arg(31)->Arg(32)->Arg(36)->Arg(40)->Arg(44)->Arg(48)->Arg(49)->Arg(50)->Arg(51) \
->Arg(52)->Arg(56)->Arg(64)->Arg(65)->Arg(68)->Arg(70)->Arg(72)->Arg(74)->Arg(76) \
->Arg(78)->Arg(80)->Arg(88)->Arg(96)->Arg(112)->Arg(128)

BENCHMARK(bench_middle_product_multipair_scalar)
MIDDLE_PRODUCT_MULTIPAIR_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_update_multipair_scalar)
MIDDLE_PRODUCT_MULTIPAIR_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_one_shot_scalar_fft)
MIDDLE_PRODUCT_AUTO_CROSSOVER_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_one_shot_scalar_multipair)
MIDDLE_PRODUCT_AUTO_CROSSOVER_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
BENCHMARK(bench_middle_product_multipair_avx)
MIDDLE_PRODUCT_MULTIPAIR_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_update_multipair_avx)
MIDDLE_PRODUCT_MULTIPAIR_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_one_shot_avx_fft)
MIDDLE_PRODUCT_AUTO_CROSSOVER_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_one_shot_avx_multipair)
MIDDLE_PRODUCT_AUTO_CROSSOVER_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
#endif
