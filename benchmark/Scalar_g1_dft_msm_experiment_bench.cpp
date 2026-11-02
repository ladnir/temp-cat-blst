#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "blst/fft/experimental_scalar_g1_dft_msm.h"

namespace {

template<std::size_t InputSize>
std::vector<blst::scalar::G1Point> make_scalar_input()
{
    std::vector<blst::scalar::G1Point> input;
    input.reserve(InputSize);
    for (std::size_t index = 0; index < InputSize; ++index) {
        blst::Fr scalar{};
        blst::Fr::from_uint64(
            scalar,
            static_cast<std::uint64_t>(0x9e3779b9ULL * (index + 1)));
        blst::Fr canonical{};
        blst::Fr::from_mont(canonical, scalar);
        blst::scalar::G1Point point{};
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
void bench_scalar_fft(benchmark::State& state)
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
void bench_scalar_msm(benchmark::State& state)
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

#define REGISTER_SCALAR_FFT(size) \
    BENCHMARK_TEMPLATE(bench_scalar_fft, size) \
        ->Name("middle_forward/scalar_fft/" #size)

#define REGISTER_SCALAR_MSM(size, width) \
    BENCHMARK_TEMPLATE(bench_scalar_msm, size, width) \
        ->Name("middle_forward/scalar_msm_w" #width "/" #size)

#define REGISTER_SCALAR_MSM_SWEEP(size) \
    REGISTER_SCALAR_MSM(size, 2); \
    REGISTER_SCALAR_MSM(size, 3); \
    REGISTER_SCALAR_MSM(size, 4); \
    REGISTER_SCALAR_MSM(size, 5); \
    REGISTER_SCALAR_MSM(size, 6); \
    REGISTER_SCALAR_MSM(size, 7); \
    REGISTER_SCALAR_MSM(size, 8); \
    REGISTER_SCALAR_MSM(size, 9); \
    REGISTER_SCALAR_MSM(size, 10)

REGISTER_SCALAR_FFT(1);
REGISTER_SCALAR_MSM(1, 2);
REGISTER_SCALAR_FFT(2);
REGISTER_SCALAR_MSM_SWEEP(2);
REGISTER_SCALAR_FFT(4);
REGISTER_SCALAR_MSM_SWEEP(4);
REGISTER_SCALAR_FFT(8);
REGISTER_SCALAR_MSM_SWEEP(8);
REGISTER_SCALAR_FFT(16);
REGISTER_SCALAR_MSM_SWEEP(16);
REGISTER_SCALAR_FFT(32);
REGISTER_SCALAR_MSM_SWEEP(32);
REGISTER_SCALAR_FFT(64);
REGISTER_SCALAR_MSM_SWEEP(64);
REGISTER_SCALAR_FFT(128);
REGISTER_SCALAR_MSM_SWEEP(128);
REGISTER_SCALAR_FFT(256);
REGISTER_SCALAR_MSM_SWEEP(256);

#undef REGISTER_SCALAR_MSM_SWEEP
#undef REGISTER_SCALAR_MSM
#undef REGISTER_SCALAR_FFT

} // namespace
