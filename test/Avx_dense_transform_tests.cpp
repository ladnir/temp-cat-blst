#include <array>
#include <gtest/gtest.h>

#include <random>
#include <span>
#include <vector>

#include "blst/fields/Fr.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/fft/g1_dense_transform.h"
#endif

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

using namespace blst;

namespace {

template<std::size_t M>
std::array<blst::avx::G1Point, M> make_avx_g1_data(uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist(1, 1024);

    std::array<blst::avx::G1Point, M> out = {};
    for (std::size_t i = 0; i < M; ++i) {
        blst::Fr scalar = {};
        blst::Fr::from_uint64(scalar, dist(rng));
        blst::Fr scalar_non_mont = {};
        blst::Fr::from_mont(scalar_non_mont, scalar);
        blst::avx::G1Point::mult(
            out[i],
            blst::avx::G1Point::generator(),
            reinterpret_cast<const blst::byte*>(&scalar_non_mont),
            256);
    }
    return out;
}

template<std::size_t M>
void expect_dense_transform_matches_fft()
{
    using AvxG1 = blst::avx::G1Point;
    using AvxCoeffCtx = blst::avx::DenseTransformCoeffCtx<AvxG1>;

    const auto plan = blst::avx::make_g1_dense_transform_plan<M>();
    const auto input = make_avx_g1_data<M>(0xD35EULL + M);

    std::array<AvxG1, M> dense_out = {};
    std::array<AvxG1, M> fft_out = input;

    blst::avx::g1_dense_transform<M>(
        std::span<const AvxG1>(input.data(), input.size()),
        plan,
        std::span<AvxG1>(dense_out.data(), dense_out.size()));

    fft<AvxG1, Fr, AvxCoeffCtx>(
        std::span<AvxG1>(fft_out.data(), fft_out.size()),
        std::span<const Fr>(plan.roots.data(), plan.roots.size()),
        false);

    for (std::size_t i = 0; i < M; ++i) {
        EXPECT_TRUE(AvxG1::equal(dense_out[i], fft_out[i]) != 0)
            << "Mismatch at size " << M << ", output index " << i;
    }
}

} // namespace

TEST(Fft, AvxDenseTransformMatchesFft)
{
    expect_dense_transform_matches_fft<2>();
    expect_dense_transform_matches_fft<4>();
    expect_dense_transform_matches_fft<8>();
    expect_dense_transform_matches_fft<16>();
    expect_dense_transform_matches_fft<32>();
}

#endif
