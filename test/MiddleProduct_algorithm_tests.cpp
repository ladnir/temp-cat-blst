#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "blst/fft/middle_product.h"

namespace {

using blst::GT;
using blst::G1Point;
using blst::G2Point;
using blst::MiddleProductAlgorithm;
using blst::MiddleProductScalarContext;

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
        return blst::precompute_middle_product_b_avx(b, n, timings);
    }
}

template<class Context>
auto multipair_precompute_for_context(
    std::span<const G2Point> b,
    std::size_t n,
    blst::TimingStats* timings)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        return blst::precompute_middle_product_b_multipair_scalar(b, n, timings);
    } else {
        return blst::precompute_middle_product_b_multipair_avx(b, n, timings);
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
        blst::middle_product_avx(a, precomputed_b, out);
    }
}

template<class Context, class PrecomputedB>
void multipair_middle_product_for_context(
    std::span<const G1Point> a,
    const PrecomputedB& precomputed_b,
    std::span<GT> out)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        blst::middle_product_multipair_scalar(a, precomputed_b, out);
    } else {
        blst::middle_product_multipair_avx(a, precomputed_b, out);
    }
}

template<class Context>
void algorithm_middle_product_for_context(
    std::span<const G1Point> a,
    std::span<const G2Point> b,
    std::span<GT> out,
    MiddleProductAlgorithm algorithm,
    blst::TimingStats* timings = nullptr)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        blst::middle_product_scalar(a, b, out, algorithm, timings);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    } else if constexpr (std::is_same_v<Context, blst::MiddleProductAvxContext>) {
        blst::middle_product_avx(a, b, out, algorithm, timings);
#endif
    } else {
        static_assert(std::is_same_v<Context, blst::MiddleProductScalarContext>);
    }
}

template<class Context, class PrecomputedB>
void multipair_update_for_context(
    const G1Point& a0,
    const PrecomputedB& precomputed_b,
    std::span<GT> out)
{
    if constexpr (std::is_same_v<Context, blst::MiddleProductScalarContext>) {
        blst::update_multipair_scalar(a0, precomputed_b, out);
    } else {
        blst::update_multipair_avx(a0, precomputed_b, out);
    }
}

template<class Context>
void expect_algorithm_dispatch_matches_explicit(
    std::span<const G1Point> a,
    std::span<const G2Point> b,
    MiddleProductAlgorithm algorithm)
{
    const std::size_t n = a.size();
    std::vector<GT> dispatch_out(n, GT::identity());
    std::vector<GT> explicit_out(n, GT::identity());

    algorithm_middle_product_for_context<Context>(a, b, std::span<GT>(dispatch_out), algorithm);

    switch (algorithm) {
    case MiddleProductAlgorithm::Fft: {
        const auto fft_precomputed = fft_precompute_for_context<Context>(b, n, nullptr);
        fft_middle_product_for_context<Context>(a, fft_precomputed, std::span<GT>(explicit_out));
        break;
    }
    case MiddleProductAlgorithm::Multipair: {
        const auto multipair_precomputed = multipair_precompute_for_context<Context>(b, n, nullptr);
        multipair_middle_product_for_context<Context>(
            a,
            multipair_precomputed,
            std::span<GT>(explicit_out));
        break;
    }
    case MiddleProductAlgorithm::Auto:
    default:
        BLST_ASSERT(false && "Explicit dispatch helper expects Fft or Multipair");
    }

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(dispatch_out[i], explicit_out[i]))
            << "algorithm dispatch mismatch at n=" << n << ", i=" << i;
    }
}

template<class Context>
void expect_auto_dispatch_uses_expected_algorithm(
    std::span<const G1Point> a,
    std::span<const G2Point> b)
{
    const std::size_t n = a.size();
    const MiddleProductAlgorithm expected_algorithm =
        (n <= blst::middle_product_auto_multipair_max_n)
        ? MiddleProductAlgorithm::Multipair
        : MiddleProductAlgorithm::Fft;

    blst::TimingStats timings;
    std::vector<GT> out(n, GT::identity());

    algorithm_middle_product_for_context<Context>(
        a,
        b,
        std::span<GT>(out),
        MiddleProductAlgorithm::Auto,
        &timings);

    if (expected_algorithm == MiddleProductAlgorithm::Multipair) {
        EXPECT_GT(timings.get_ns("middle_product_multipair.total_ns"), 0);
        EXPECT_EQ(timings.get_ns("middle_product.total_ns"), 0);
    } else {
        EXPECT_GT(timings.get_ns("middle_product.total_ns"), 0);
        EXPECT_EQ(timings.get_ns("middle_product_multipair.total_ns"), 0);
    }
}

template<class Context>
void expect_multipair_matches_fft(
    std::span<const G1Point> a,
    std::span<const G2Point> b)
{
    const std::size_t n = a.size();
    const auto multipair_precomputed =
        multipair_precompute_for_context<Context>(b, n, nullptr);
    const auto fft_precomputed =
        fft_precompute_for_context<Context>(b, n, nullptr);

    std::vector<GT> multipair_out(n, GT::identity());
    std::vector<GT> fft_out(n, GT::identity());

    multipair_middle_product_for_context<Context>(
        a,
        multipair_precomputed,
        std::span<GT>(multipair_out));
    fft_middle_product_for_context<Context>(
        a,
        fft_precomputed,
        std::span<GT>(fft_out));

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(multipair_out[i], fft_out[i]))
            << "multipair mismatch at n=" << n << ", i=" << i;
    }
}

template<class Context>
void expect_multipair_update_and_combine_matches_full(
    std::span<const G1Point> a,
    std::span<const G2Point> b)
{
    const std::size_t n = a.size();
    const auto multipair_precomputed =
        multipair_precompute_for_context<Context>(b, n, nullptr);

    std::vector<GT> expected(n, GT::identity());
    std::vector<GT> mid(n, GT::identity());
    std::vector<GT> cnst(n, GT::identity());
    std::vector<GT> combined(n, GT::identity());

    std::vector<G1Point> a_without_constant(a.begin(), a.end());
    a_without_constant[0] = G1Point::infinity();

    multipair_middle_product_for_context<Context>(
        a,
        multipair_precomputed,
        std::span<GT>(expected));
    multipair_middle_product_for_context<Context>(
        std::span<const G1Point>(a_without_constant),
        multipair_precomputed,
        std::span<GT>(mid));
    multipair_update_for_context<Context>(a[0], multipair_precomputed, std::span<GT>(cnst));
    blst::combine(
        std::span<const GT>(mid),
        std::span<const GT>(cnst),
        std::span<GT>(combined));

    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(expected[i], combined[i]))
            << "multipair update/combine mismatch at n=" << n << ", i=" << i;
    }
}

} // namespace

TEST(MiddleProductSizing, RejectsUnsupportedLengthsWithoutOverflow)
{
    EXPECT_THROW(blst::middle_product_fft_size(0), std::invalid_argument);
    EXPECT_EQ(
        blst::middle_product_fft_size(blst::middle_product_max_n),
        blst::middle_product_max_fft_size);
    EXPECT_THROW(
        blst::middle_product_fft_size(blst::middle_product_max_n + 1),
        std::length_error);
    EXPECT_THROW(
        blst::middle_product_fft_size(std::numeric_limits<std::size_t>::max()),
        std::length_error);
}

TEST(MiddleProductSizing, BasicOracleAndAdmissibleSizesAreChecked)
{
    constexpr std::size_t basic_max_n =
        (blst::middle_product_max_fft_size + 1) / 3;

    EXPECT_EQ(
        blst::middle_product_basic_fft_size(basic_max_n),
        blst::middle_product_max_fft_size);
    EXPECT_THROW(
        blst::middle_product_basic_fft_size(basic_max_n + 1),
        std::length_error);
    EXPECT_EQ(
        blst::fft_admissible_size(std::numeric_limits<std::size_t>::max(), 1024),
        1024);
    EXPECT_THROW(blst::fft_admissible_size(1, 0), std::invalid_argument);
    EXPECT_THROW(blst::fft_admissible_size(1, 3), std::invalid_argument);
}

TEST(MiddleProductPrecompute, InverseRootsReverseForwardPowers)
{
    EXPECT_TRUE(blst::make_ifft_roots(std::span<const blst::Fr>{}).empty());

    for (const std::size_t size : {std::size_t{1}, std::size_t{2}, std::size_t{8}}) {
        const std::vector<blst::Fr> roots = blst::make_fft_roots(size);
        const std::vector<blst::Fr> inv_roots = blst::make_ifft_roots(roots);

        ASSERT_EQ(inv_roots.size(), roots.size());
        EXPECT_TRUE(blst::Fr::equal(inv_roots[0], roots[0]));
        for (std::size_t i = 1; i < size; ++i) {
            EXPECT_TRUE(blst::Fr::equal(inv_roots[i], roots[size - i]));
            blst::Fr product = {};
            blst::Fr::mul(product, roots[i], inv_roots[i]);
            EXPECT_TRUE(blst::Fr::equal(product, blst::Fr::one()));
        }
    }
}

TEST(MiddleProductPrecompute, ScalarPlansRemainStableAcrossEvaluations)
{
    constexpr std::size_t n = 3;
    const std::vector<G1Point> a = make_g1_data(n);
    const std::vector<G2Point> b = make_g2_data(n);
    const auto precomputed = blst::precompute_middle_product_b_scalar(b, n);

    EXPECT_EQ(precomputed.forward_plan.n, precomputed.fft_size);
    EXPECT_EQ(
        precomputed.forward_plan.xn,
        blst::fft_admissible_size(n, precomputed.fft_size));
    EXPECT_EQ(precomputed.inverse_plan.n, precomputed.fft_size);
    EXPECT_EQ(precomputed.inverse_plan.start, n);
    EXPECT_EQ(precomputed.inverse_plan.length, n);

    const auto* const forward_storage = precomputed.forward_plan.stage_work.data();
    const auto* const inverse_storage = precomputed.inverse_plan.stage_left.data();
    std::vector<GT> first(n, GT::identity());
    std::vector<GT> second(n, GT::identity());
    blst::MiddleProductWorkspaceScalar first_workspace;
    blst::MiddleProductWorkspaceScalar second_workspace;

    blst::middle_product_scalar(a, precomputed, first, first_workspace);
    blst::middle_product_scalar(a, precomputed, second, second_workspace);

    EXPECT_EQ(precomputed.forward_plan.stage_work.data(), forward_storage);
    EXPECT_EQ(precomputed.inverse_plan.stage_left.data(), inverse_storage);
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(GT::equal(first[i], second[i]));
    }
}

TEST(MiddleProductMultipair, ScalarMatchesFft)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_multipair_matches_fft<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}

TEST(MiddleProductMultipair, ScalarUpdateAndCombineReconstructsFullResult)
{
    const std::vector<std::size_t> sizes = { 4, 6, 10, 18 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_multipair_update_and_combine_matches_full<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}

TEST(MiddleProductAlgorithm, ScalarFftMatchesExplicitFft)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16, 33 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_algorithm_dispatch_matches_explicit<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductAlgorithm::Fft);
    }
}

TEST(MiddleProductAlgorithm, ScalarMultipairMatchesExplicitMultipair)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16, 24 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_algorithm_dispatch_matches_explicit<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductAlgorithm::Multipair);
    }
}

TEST(MiddleProductAlgorithm, ScalarAutoUsesExpectedAlgorithm)
{
    const std::vector<std::size_t> sizes = {
        blst::middle_product_auto_multipair_max_n,
        blst::middle_product_auto_multipair_max_n + 1,
    };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_auto_dispatch_uses_expected_algorithm<MiddleProductScalarContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
TEST(MiddleProductMultipair, AvxPortableMatchesFft)
{
    const std::vector<std::size_t> sizes = { 3, 5, 9, 16 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_multipair_matches_fft<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}

TEST(MiddleProductMultipair, AvxPortableUpdateAndCombineReconstructsFullResult)
{
    const std::vector<std::size_t> sizes = { 4, 6, 10, 18 };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_multipair_update_and_combine_matches_full<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}

TEST(MiddleProductAlgorithm, AvxPortableFftMatchesExplicitFft)
{
    const std::vector<std::size_t> sizes = {
        3,
        16,
        blst::middle_product_auto_multipair_max_n + 1,
    };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_algorithm_dispatch_matches_explicit<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductAlgorithm::Fft);
    }
}

TEST(MiddleProductAlgorithm, AvxPortableMultipairMatchesExplicitMultipair)
{
    const std::vector<std::size_t> sizes = {
        3,
        16,
        blst::middle_product_auto_multipair_max_n,
    };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_algorithm_dispatch_matches_explicit<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b),
            MiddleProductAlgorithm::Multipair);
    }
}

TEST(MiddleProductAlgorithm, AvxPortableAutoUsesExpectedAlgorithm)
{
    const std::vector<std::size_t> sizes = {
        blst::middle_product_auto_multipair_max_n,
        blst::middle_product_auto_multipair_max_n + 1,
    };

    for (std::size_t n : sizes) {
        const std::vector<G1Point> a = make_g1_data(n);
        const std::vector<G2Point> b = make_g2_data(n);
        expect_auto_dispatch_uses_expected_algorithm<blst::MiddleProductAvxContext>(
            std::span<const G1Point>(a),
            std::span<const G2Point>(b));
    }
}
#endif
