#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "blst/fft/experimental_g1_dft_msm.h"

namespace {

template<std::size_t InputSize>
std::vector<blst::avx::G1Point> make_input()
{
    std::vector<blst::avx::G1Point> input;
    input.reserve(InputSize);
    for (std::size_t index = 0; index < InputSize; ++index) {
        blst::Fr scalar = {};
        blst::Fr::from_uint64(
            scalar,
            static_cast<std::uint64_t>(0x9e3779b9ULL * (index + 1)));
        blst::Fr canonical = {};
        blst::Fr::from_mont(canonical, scalar);
        blst::avx::G1Point point = {};
        blst::avx::G1Point::mult_vartime(
            point,
            blst::avx::G1Point::generator(),
            reinterpret_cast<const blst::byte*>(&canonical),
            256);
        input.push_back(point);
    }
    return input;
}

template<std::size_t InputSize>
std::vector<blst::scalar::G1Point> make_scalar_input()
{
    std::vector<blst::scalar::G1Point> input;
    input.reserve(InputSize);
    for (std::size_t index = 0; index < InputSize; ++index) {
        blst::Fr scalar = {};
        blst::Fr::from_uint64(
            scalar,
            static_cast<std::uint64_t>(0x9e3779b9ULL * (index + 1)));
        blst::Fr canonical = {};
        blst::Fr::from_mont(canonical, scalar);
        blst::scalar::G1Point point = {};
        blst::scalar::G1Point::mult_vartime(
            point,
            blst::scalar::G1Point::generator(),
            reinterpret_cast<const blst::byte*>(&canonical),
            256);
        input.push_back(point);
    }
    return input;
}

template<std::size_t InputSize>
void bench_middle_forward_fft(benchmark::State& state)
{
    constexpr std::size_t transform_size = 2 * InputSize;
    const std::vector<blst::Fr> roots = blst::make_fft_roots(transform_size);
    const std::vector<blst::avx::G1Point> input = make_input<InputSize>();
    const std::size_t xn =
        blst::fft_admissible_size(InputSize, transform_size);
    const FftTruncatedForwardPlan plan =
        make_fft_truncated_forward_plan(transform_size, xn);
    std::vector<blst::avx::G1Point> output(
        transform_size, blst::avx::G1Point::infinity());

    for (auto _ : state) {
        std::copy(input.begin(), input.end(), output.begin());
        std::fill(
            output.begin() + InputSize,
            output.end(),
            blst::avx::G1Point::infinity());
        benchmark::DoNotOptimize(output.data());
        fft_forward_truncated_prefix<
            blst::avx::G1Point,
            blst::Fr,
            blst::detail::GroupFftCoeffCtx<blst::avx::G1Point>>(
                std::span<blst::avx::G1Point>(output),
                std::span<const blst::Fr>(roots),
                xn,
                plan);
        benchmark::ClobberMemory();
    }
    state.counters["batch"] = static_cast<double>(InputSize);
    state.counters["transform"] = static_cast<double>(transform_size);
}

template<std::size_t InputSize>
void bench_middle_forward_msm(benchmark::State& state)
{
    constexpr std::size_t transform_size = 2 * InputSize;
    const std::vector<blst::Fr> roots = blst::make_fft_roots(transform_size);
    const std::vector<blst::avx::G1Point> input = make_input<InputSize>();
    blst::experimental::AvxG1MiddleDftMsm<InputSize> transform(roots);
    std::vector<blst::avx::G1Point> output(transform_size);

    for (auto _ : state) {
        benchmark::DoNotOptimize(input.data());
        transform.transform(output, input);
        benchmark::DoNotOptimize(output.data());
        benchmark::ClobberMemory();
    }
    state.counters["batch"] = static_cast<double>(InputSize);
    state.counters["transform"] = static_cast<double>(transform_size);
}

template<std::size_t InputSize, int Width>
void bench_middle_forward_msm_width(benchmark::State& state)
{
    constexpr std::size_t transform_size = 2 * InputSize;
    const std::vector<blst::Fr> roots = blst::make_fft_roots(transform_size);
    const std::vector<blst::avx::G1Point> input = make_input<InputSize>();
    blst::experimental::AvxG1MiddleDftMsm<InputSize, Width> transform(roots);
    std::vector<blst::avx::G1Point> output(transform_size);

    for (auto _ : state) {
        benchmark::DoNotOptimize(input.data());
        transform.transform(output, input);
        benchmark::DoNotOptimize(output.data());
        benchmark::ClobberMemory();
    }
    state.counters["batch"] = static_cast<double>(InputSize);
    state.counters["transform"] = static_cast<double>(transform_size);
    state.counters["width"] = static_cast<double>(Width);
}

template<std::size_t InputSize, int Width>
void bench_middle_forward_msm_prepared(benchmark::State& state)
{
    constexpr std::size_t transform_size = 2 * InputSize;
    const std::vector<blst::Fr> roots = blst::make_fft_roots(transform_size);
    const std::vector<blst::avx::G1Point> input = make_input<InputSize>();
    blst::experimental::AvxG1MiddleDftMsm<InputSize, Width> transform(roots);
    transform.prepare_bases(input);
    std::vector<blst::avx::G1Point> output(transform_size);

    for (auto _ : state) {
        transform.transform_prepared(output);
        benchmark::DoNotOptimize(output.data());
        benchmark::ClobberMemory();
    }
    state.counters["batch"] = static_cast<double>(InputSize);
    state.counters["transform"] = static_cast<double>(transform_size);
    state.counters["width"] = static_cast<double>(Width);
}

template<std::size_t InputSize>
void bench_middle_forward_pippenger(benchmark::State& state)
{
    constexpr std::size_t transform_size = 2 * InputSize;
    const std::vector<blst::Fr> roots = blst::make_fft_roots(transform_size);
    const std::vector<blst::avx::G1Point> input = make_input<InputSize>();
    blst::experimental::AvxG1MiddleDftPippenger<InputSize> transform(roots);
    std::vector<blst::avx::G1Point> output(transform_size);

    for (auto _ : state) {
        benchmark::DoNotOptimize(input.data());
        transform.transform(output, input);
        benchmark::DoNotOptimize(output.data());
        benchmark::ClobberMemory();
    }
    state.counters["batch"] = static_cast<double>(InputSize);
    state.counters["transform"] = static_cast<double>(transform_size);
}

template<std::size_t InputSize>
void bench_middle_forward_scalar_fft(benchmark::State& state)
{
    constexpr std::size_t transform_size = 2 * InputSize;
    const std::vector<blst::Fr> roots = blst::make_fft_roots(transform_size);
    const std::vector<blst::scalar::G1Point> input =
        make_scalar_input<InputSize>();
    const std::size_t xn =
        blst::fft_admissible_size(InputSize, transform_size);
    const FftTruncatedForwardPlan plan =
        make_fft_truncated_forward_plan(transform_size, xn);
    std::vector<blst::scalar::G1Point> output(
        transform_size, blst::scalar::G1Point::infinity());

    for (auto _ : state) {
        std::copy(input.begin(), input.end(), output.begin());
        std::fill(
            output.begin() + InputSize,
            output.end(),
            blst::scalar::G1Point::infinity());
        benchmark::DoNotOptimize(output.data());
        fft_forward_truncated_prefix<
            blst::scalar::G1Point,
            blst::Fr,
            blst::detail::GroupFftCoeffCtx<blst::scalar::G1Point>>(
                std::span<blst::scalar::G1Point>(output),
                std::span<const blst::Fr>(roots),
                xn,
                plan);
        benchmark::ClobberMemory();
    }
    state.counters["batch"] = static_cast<double>(InputSize);
    state.counters["transform"] = static_cast<double>(transform_size);
}

template<std::size_t InputSize, std::size_t Wbits>
void bench_middle_forward_scalar_fixed(benchmark::State& state)
{
    constexpr std::size_t transform_size = 2 * InputSize;
    const std::vector<blst::Fr> roots = blst::make_fft_roots(transform_size);
    const std::vector<blst::scalar::G1Point> input =
        make_scalar_input<InputSize>();
    blst::experimental::ScalarG1MiddleDftFixedTable<InputSize, Wbits>
        transform(roots);
    std::vector<blst::scalar::G1Point> output(transform_size);

    for (auto _ : state) {
        benchmark::DoNotOptimize(input.data());
        transform.transform(output, input);
        benchmark::DoNotOptimize(output.data());
        benchmark::ClobberMemory();
    }
    state.counters["batch"] = static_cast<double>(InputSize);
    state.counters["transform"] = static_cast<double>(transform_size);
    state.counters["wbits"] = static_cast<double>(Wbits);
}

#define REGISTER_SIZE(size) \
    BENCHMARK_TEMPLATE(bench_middle_forward_fft, size) \
        ->Name("middle_forward/fft/" #size); \
    BENCHMARK_TEMPLATE(bench_middle_forward_msm, size) \
        ->Name("middle_forward/msm_wnaf/" #size); \
    BENCHMARK_TEMPLATE(bench_middle_forward_pippenger, size) \
        ->Name("middle_forward/msm_pippenger/" #size)

REGISTER_SIZE(1);
REGISTER_SIZE(2);
REGISTER_SIZE(4);
REGISTER_SIZE(8);
REGISTER_SIZE(16);
REGISTER_SIZE(32);
REGISTER_SIZE(64);
REGISTER_SIZE(128);
REGISTER_SIZE(256);

#undef REGISTER_SIZE

#define REGISTER_WIDTH(size, width) \
    BENCHMARK_TEMPLATE(bench_middle_forward_msm_width, size, width) \
        ->Name("middle_forward/msm_w" #width "/" #size)

#define REGISTER_WIDTH_SWEEP(size) \
    REGISTER_WIDTH(size, 4); \
    REGISTER_WIDTH(size, 6); \
    REGISTER_WIDTH(size, 7); \
    REGISTER_WIDTH(size, 8)

REGISTER_WIDTH_SWEEP(2);
REGISTER_WIDTH_SWEEP(4);
REGISTER_WIDTH_SWEEP(8);
REGISTER_WIDTH_SWEEP(16);
REGISTER_WIDTH_SWEEP(32);

#undef REGISTER_WIDTH_SWEEP
#undef REGISTER_WIDTH

#define REGISTER_PREPARED(size, width) \
    BENCHMARK_TEMPLATE(bench_middle_forward_msm_prepared, size, width) \
        ->Name("middle_forward/msm_prepared_w" #width "/" #size)

REGISTER_PREPARED(2, 4);
REGISTER_PREPARED(4, 6);
REGISTER_PREPARED(8, 7);
REGISTER_PREPARED(16, 8);
REGISTER_PREPARED(32, 8);

#undef REGISTER_PREPARED

#define REGISTER_SCALAR_FFT(size) \
    BENCHMARK_TEMPLATE(bench_middle_forward_scalar_fft, size) \
        ->Name("middle_forward/scalar_fft/" #size)

#define REGISTER_SCALAR_FIXED(size, width) \
    BENCHMARK_TEMPLATE(bench_middle_forward_scalar_fixed, size, width) \
        ->Name("middle_forward/scalar_fixed_w" #width "/" #size)

#define REGISTER_SCALAR_FIXED_SWEEP(size) \
    REGISTER_SCALAR_FIXED(size, 2); \
    REGISTER_SCALAR_FIXED(size, 3); \
    REGISTER_SCALAR_FIXED(size, 4); \
    REGISTER_SCALAR_FIXED(size, 5); \
    REGISTER_SCALAR_FIXED(size, 6); \
    REGISTER_SCALAR_FIXED(size, 7); \
    REGISTER_SCALAR_FIXED(size, 8)

REGISTER_SCALAR_FFT(1);
REGISTER_SCALAR_FIXED(1, 2);
REGISTER_SCALAR_FFT(2);
REGISTER_SCALAR_FIXED_SWEEP(2);
REGISTER_SCALAR_FFT(4);
REGISTER_SCALAR_FIXED_SWEEP(4);
REGISTER_SCALAR_FFT(8);
REGISTER_SCALAR_FIXED_SWEEP(8);
REGISTER_SCALAR_FFT(16);
REGISTER_SCALAR_FIXED_SWEEP(16);
REGISTER_SCALAR_FFT(32);
REGISTER_SCALAR_FIXED_SWEEP(32);

#undef REGISTER_SCALAR_FIXED_SWEEP
#undef REGISTER_SCALAR_FIXED
#undef REGISTER_SCALAR_FFT

} // namespace
