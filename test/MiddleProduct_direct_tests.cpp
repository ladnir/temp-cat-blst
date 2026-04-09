#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

#include "MiddleProduct_reference.h"

namespace {

using blst::GT;
using blst::G1Point;
using blst::G2Point;
using blst::MiddleProductPairingMode;
using blst::MiddleProductScalarContext;
using blst::test_support::combine_direct;
using blst::test_support::middle_product_direct;
using blst::test_support::precompute_middle_product_b_direct;
using blst::test_support::update_direct;

G1Point g1(uint64_t scalar)
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

G2Point g2(uint64_t scalar)
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

std::vector<G1Point> make_g1_data(std::size_t n)
{
    std::vector<G1Point> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        if ((i % 7) == 3) {
            out.push_back(G1Point::infinity());
            continue;
        }
        out.push_back(g1(static_cast<uint64_t>((37 * i + 11) % 257 + 1)));
    }
    return out;
}

std::vector<G2Point> make_g2_data(std::size_t n)
{
    std::vector<G2Point> out;
    out.reserve(2 * n);
    for (std::size_t i = 0; i < (2 * n); ++i) {
        if ((i % 11) == 5) {
            out.push_back(G2Point::infinity());
            continue;
        }
        out.push_back(g2(static_cast<uint64_t>((43 * i + 17) % 263 + 1)));
    }
    return out;
}

template<class Context>
auto fft_precompute_for_context(
    std::span<const G2Point> b,
    std::size_t n,
    blst::TimingStats* timings)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        return blst::precompute_middle_product_b_scalar(b, n, timings);
    } else {
        return blst::precompute_middle_product_b(b, n, timings);
    }
}

template<class Context, class PrecomputedB>
void fft_middle_product_for_context(
    std::span<const G1Point> a,
    const PrecomputedB& precomputed_b,
    std::span<GT> out)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        blst::middle_product_scalar(a, precomputed_b, out);
    } else {
        blst::middle_product(a, precomputed_b, out);
    }
}

template<class Context>
void expect_direct_matches_fft(
    std::span<const G1Point> a,
    std::span<const G2Point> b,
    MiddleProductPairingMode pairing_mode)
{
    const std::size_t n = a.size();
    const auto direct_precomputed =
        precompute_middle_product_b_direct<Context>(b, n, nullptr, pairing_mode);
    const auto fft_precomputed =
        fft_precompute_for_context<Context>(b, n, nullptr);

    std::vector<GT> direct_out(n, GT::identity());
    std::vector<GT> fft_out(n, GT::identity());

    middle_product_direct<Context>(a, direct_precomputed, std::span<GT>(direct_out));
    fft_middle_product_for_context<Context>(a, fft_precomputed, std::span<GT>(fft_out));

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(direct_out[i], fft_out[i]))
            << "direct mismatch at n=" << n << ", i=" << i;
    }
}

template<class Context>
void expect_update_and_combine_matches_fft(
    std::span<const G1Point> a,
    std::span<const G2Point> b,
    MiddleProductPairingMode pairing_mode)
{
    const std::size_t n = a.size();
    const auto direct_precomputed =
        precompute_middle_product_b_direct<Context>(b, n, nullptr, pairing_mode);

    std::vector<GT> expected(n, GT::identity());
    std::vector<GT> mid(n, GT::identity());
    std::vector<GT> cnst(n, GT::identity());
    std::vector<GT> combined(n, GT::identity());

    std::vector<G1Point> a_without_constant(a.begin(), a.end());
    a_without_constant[0] = G1Point::infinity();

    middle_product_direct<Context>(a, direct_precomputed, std::span<GT>(expected));
    middle_product_direct<Context>(
        std::span<const G1Point>(a_without_constant),
        direct_precomputed,
        std::span<GT>(mid));
    update_direct<Context>(a[0], direct_precomputed, std::span<GT>(cnst));
    combine_direct<Context>(
        std::span<const GT>(mid),
        std::span<const GT>(cnst),
        std::span<GT>(combined));

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(expected[i], combined[i]))
            << "update/combine mismatch at n=" << n << ", i=" << i;
    }
}

} // namespace

TEST(MiddleProductDirect, ScalarMatchesScalarFft)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_direct_matches_fft<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductPairingMode::Direct);
    }
}

TEST(MiddleProductDirect, ScalarPrecomputedLinesMatchesScalarFft)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_direct_matches_fft<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductPairingMode::PrecomputedLines);
    }
}

TEST(MiddleProductDirect, ScalarUpdateAndCombineReconstructsFullResult)
{
    const std::vector<std::size_t> sizes = { 4, 6, 10, 18 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_update_and_combine_matches_fft<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductPairingMode::PrecomputedLines);
    }
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
TEST(MiddleProductDirect, AvxPortableMatchesDefaultFft)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_direct_matches_fft<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductPairingMode::Direct);
    }
}

TEST(MiddleProductDirect, AvxPortablePrecomputedLinesMatchesDefaultFft)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_direct_matches_fft<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductPairingMode::PrecomputedLines);
    }
}
#endif
