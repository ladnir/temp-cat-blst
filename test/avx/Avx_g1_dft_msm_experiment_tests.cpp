#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "blst/fft/experimental_g1_dft_msm.h"

namespace {

template<std::size_t InputSize, int Width = 5>
void expect_msm_dft_matches_middle_product_fft()
{
    constexpr std::size_t transform_size = 2 * InputSize;
    const std::vector<blst::Fr> roots = blst::make_fft_roots(transform_size);
    std::vector<blst::avx::G1Point> input;
    input.reserve(InputSize);
    for (std::size_t index = 0; index < InputSize; ++index) {
        blst::Fr scalar = {};
        blst::Fr::from_uint64(
            scalar,
            static_cast<std::uint64_t>(17 * index + 3));
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

    std::vector<blst::avx::G1Point> expected(
        transform_size, blst::avx::G1Point::infinity());
    std::copy(input.begin(), input.end(), expected.begin());
    const std::size_t xn =
        blst::fft_admissible_size(InputSize, transform_size);
    const FftTruncatedForwardPlan plan =
        make_fft_truncated_forward_plan(transform_size, xn);
    fft_forward_truncated_prefix<
        blst::avx::G1Point,
        blst::Fr,
        blst::detail::GroupFftCoeffCtx<blst::avx::G1Point>>(
            std::span<blst::avx::G1Point>(expected),
            std::span<const blst::Fr>(roots),
            xn,
            plan);

    blst::experimental::AvxG1MiddleDftMsm<InputSize, Width> transform(roots);
    std::vector<blst::avx::G1Point> actual(transform_size);
    transform.transform(actual, input);

    for (std::size_t index = 0; index < transform_size; ++index) {
        EXPECT_NE(
            blst::avx::G1Point::equal(expected[index], actual[index]),
            0u)
            << "transform output " << index;
    }

    blst::experimental::AvxG1MiddleDftPippenger<InputSize>
        pippenger_transform(roots);
    pippenger_transform.transform(actual, input);
    for (std::size_t index = 0; index < transform_size; ++index) {
        EXPECT_NE(
            blst::avx::G1Point::equal(expected[index], actual[index]),
            0u)
            << "Pippenger transform output " << index;
    }
}

template<std::size_t InputSize, std::size_t Wbits>
void expect_scalar_fixed_dft_matches_middle_product_fft()
{
    constexpr std::size_t transform_size = 2 * InputSize;
    const std::vector<blst::Fr> roots = blst::make_fft_roots(transform_size);
    std::vector<blst::scalar::G1Point> input;
    input.reserve(InputSize);
    for (std::size_t index = 0; index < InputSize; ++index) {
        blst::Fr scalar = {};
        blst::Fr::from_uint64(
            scalar,
            static_cast<std::uint64_t>(17 * index + 3));
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

    std::vector<blst::scalar::G1Point> expected(
        transform_size, blst::scalar::G1Point::infinity());
    std::copy(input.begin(), input.end(), expected.begin());
    const std::size_t xn =
        blst::fft_admissible_size(InputSize, transform_size);
    const FftTruncatedForwardPlan plan =
        make_fft_truncated_forward_plan(transform_size, xn);
    fft_forward_truncated_prefix<
        blst::scalar::G1Point,
        blst::Fr,
        blst::detail::GroupFftCoeffCtx<blst::scalar::G1Point>>(
            std::span<blst::scalar::G1Point>(expected),
            std::span<const blst::Fr>(roots),
            xn,
            plan);

    blst::experimental::ScalarG1MiddleDftFixedTable<InputSize, Wbits>
        transform(roots);
    std::vector<blst::scalar::G1Point> actual(transform_size);
    transform.transform(actual, input);
    for (std::size_t index = 0; index < transform_size; ++index) {
        EXPECT_NE(
            blst::scalar::G1Point::equal(expected[index], actual[index]),
            0u)
            << "scalar fixed-table transform output " << index;
    }
}

TEST(AvxG1DftMsmExperiment, MatchesMiddleProductForwardTransform)
{
    expect_msm_dft_matches_middle_product_fft<1>();
    expect_msm_dft_matches_middle_product_fft<2>();
    expect_msm_dft_matches_middle_product_fft<4>();
    expect_msm_dft_matches_middle_product_fft<8>();
    expect_msm_dft_matches_middle_product_fft<16>();
    expect_msm_dft_matches_middle_product_fft<32>();
    expect_msm_dft_matches_middle_product_fft<64>();
    expect_msm_dft_matches_middle_product_fft<128>();
    expect_msm_dft_matches_middle_product_fft<256>();
    expect_msm_dft_matches_middle_product_fft<8, 4>();
    expect_msm_dft_matches_middle_product_fft<8, 6>();
    expect_msm_dft_matches_middle_product_fft<8, 7>();
    expect_msm_dft_matches_middle_product_fft<8, 8>();
}

TEST(AvxG1DftMsmExperiment, RejectsAliasingAndWrongShapes)
{
    constexpr std::size_t input_size = 8;
    constexpr std::size_t transform_size = 2 * input_size;
    const std::vector<blst::Fr> roots = blst::make_fft_roots(transform_size);
    blst::experimental::AvxG1MiddleDftMsm<input_size> transform(roots);
    std::vector<blst::avx::G1Point> storage(transform_size);

    EXPECT_THROW(
        transform.transform(
            std::span<blst::avx::G1Point>(storage),
            std::span<const blst::avx::G1Point>(storage.data(), input_size)),
        std::invalid_argument);
    EXPECT_THROW(
        transform.transform(
            std::span<blst::avx::G1Point>(storage),
            std::span<const blst::avx::G1Point>(storage.data(), input_size - 1)),
        std::invalid_argument);
}

TEST(AvxG1DftMsmExperiment, ScalarFixedTableMatchesScalarFft)
{
    expect_scalar_fixed_dft_matches_middle_product_fft<2, 2>();
    expect_scalar_fixed_dft_matches_middle_product_fft<4, 4>();
    expect_scalar_fixed_dft_matches_middle_product_fft<8, 6>();
    expect_scalar_fixed_dft_matches_middle_product_fft<16, 7>();
    expect_scalar_fixed_dft_matches_middle_product_fft<32, 8>();
}

} // namespace
