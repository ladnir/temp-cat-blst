#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/domain.h"
#include "blst/fft/fr_interpolation.h"

#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{
    using blst::Fr;
    using InterpolationStatus = blst::fr_interpolation::Status;

    struct RootTestField
    {
        std::uint64_t value = 0;

        constexpr RootTestField() = default;
        constexpr RootTestField(const std::uint64_t input) : value(input % order()) {}

        static constexpr std::uint64_t order() { return 17; }
        constexpr std::uint64_t integer() const { return value; }

        constexpr RootTestField pow(std::uint64_t exponent) const
        {
            RootTestField result(1);
            RootTestField base = *this;
            while (exponent != 0) {
                if ((exponent & 1) != 0) {
                    result *= base;
                }
                base *= base;
                exponent >>= 1;
            }
            return result;
        }

        constexpr RootTestField& operator*=(const RootTestField rhs)
        {
            value = value * rhs.value % order();
            return *this;
        }

        friend constexpr bool operator==(
            const RootTestField&, const RootTestField&) = default;
        friend constexpr bool operator==(const RootTestField lhs, const std::uint64_t rhs)
        {
            return lhs.value == rhs % order();
        }
    };

    static_assert(blst::fft_detail::valid_size(1));
    static_assert(blst::fft_detail::valid_size(8));
    static_assert(!blst::fft_detail::valid_size(0));
    static_assert(!blst::fft_detail::valid_size(3));
    static_assert(blst::fft_detail::domain_storage_len(8) == 10);
    static_assert(
        blst::fr_interpolation::domain_suffix_scratch_len(5, 8) == 33);
    static_assert(blst::log2ceil(1) == 0);
    static_assert(blst::log2ceil(2) == 1);
    static_assert(blst::log2ceil(3) == 2);
    static_assert(blst::log2ceil(std::numeric_limits<std::uint64_t>::max()) == 64);
    static_assert(blst::ceil_div(0, 1) == 0);
    static_assert(blst::ceil_div(1, 1) == 1);
    static_assert(blst::ceil_div(std::numeric_limits<std::size_t>::max(), 2) ==
        std::numeric_limits<std::size_t>::max() / 2 + 1);
#ifdef NDEBUG
    static_assert(blst::log2ceil(0) == 0);
    static_assert(blst::ceil_div(1, 0) == 0);
#endif

    Fr from_u64(const std::uint64_t value)
    {
        return Fr::from_uint64(value);
    }

    Fr add(const Fr& lhs, const Fr& rhs)
    {
        Fr out;
        Fr::add(out, lhs, rhs);
        return out;
    }

    Fr multiply(const Fr& lhs, const Fr& rhs)
    {
        Fr out;
        Fr::mul(out, lhs, rhs);
        return out;
    }

    Fr evaluate_polynomial(
        const std::span<const Fr> coefficients,
        const Fr& point)
    {
        Fr value = {};
        for (auto coefficient = coefficients.rbegin();
             coefficient != coefficients.rend();
             ++coefficient) {
            value = add(multiply(value, point), *coefficient);
        }
        return value;
    }

    template<std::size_t N>
    bool equal(
        const std::array<Fr, N>& lhs,
        const std::array<Fr, N>& rhs)
    {
        for (std::size_t i = 0; i < N; ++i) {
            if (!Fr::equal(lhs[i], rhs[i])) {
                return false;
            }
        }
        return true;
    }

    bool equal(
        const std::span<const Fr> lhs,
        const std::span<const Fr> rhs)
    {
        if (lhs.size() != rhs.size()) {
            return false;
        }
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            if (!Fr::equal(lhs[i], rhs[i])) {
                return false;
            }
        }
        return true;
    }
}

TEST(FftApi, RootPreparationRejectsEmptyDomain)
{
    const RootTestField one(1);
    EXPECT_THROW(
        (blst::precomputeRootsOfUnity<RootTestField>(
            one,
            std::span<RootTestField>{})),
        std::invalid_argument);

    std::array<RootTestField, 1> roots = {};
    EXPECT_NO_THROW(blst::precomputeRootsOfUnity<RootTestField>(one, roots));
    EXPECT_EQ(roots[0], one);
}

TEST(FftApi, RootPreparationRejectsNonprimitiveRootWithoutWriting)
{
    const RootTestField one(1);
    std::array<RootTestField, 4> roots = {
        RootTestField(2), RootTestField(3), RootTestField(4), RootTestField(5)
    };
    const auto original = roots;

    EXPECT_THROW(
        blst::precomputeRootsOfUnity<RootTestField>(one, roots),
        std::invalid_argument);
    EXPECT_EQ(roots, original);
}

TEST(FftApi, ExplicitGeneratorProducesAnExactOrderRoot)
{
    const auto root = blst::primRootOfUnity<RootTestField>(8, RootTestField(3));
    EXPECT_TRUE(blst::isPrimRootOfUnity<RootTestField>(8, root));
    EXPECT_THROW(
        (void)blst::primRootOfUnity<RootTestField>(4, RootTestField(1)),
        std::invalid_argument);
}

TEST(FftApi, RootBuilderRejectsUnsupportedTransformSizes)
{
    EXPECT_THROW(
        (static_cast<void>(
            blst::fft_detail::make_roots<Fr, blst::FrCoeffCtx>(0, false))),
        std::invalid_argument);
    EXPECT_THROW(
        (static_cast<void>(
            blst::fft_detail::make_roots<Fr, blst::FrCoeffCtx>(3, false))),
        std::invalid_argument);
}

TEST(CommonDefines, InvalidSizingIsFailClosed)
{
#ifdef NDEBUG
    EXPECT_EQ(blst::log2ceil(0), 0);
    EXPECT_EQ(blst::ceil_div(1, 0), 0);
#else
    EXPECT_DEATH((void)blst::log2ceil(0), "");
    EXPECT_DEATH((void)blst::ceil_div(1, 0), "");
#endif
}

TEST(FftApi, GenericPreparedDomainMatchesGenericTransform)
{
    constexpr std::size_t n = 8;
    std::array<Fr, n> input = {};
    for (std::size_t i = 0; i < n; ++i) {
        input[i] = from_u64(i * 17 + 3);
    }

    auto expected = input;
    const std::vector<Fr> roots =
        blst::fft_detail::make_roots<Fr, blst::FrCoeffCtx>(n, false);
    fft<Fr, Fr, blst::FrCoeffCtx>(expected, roots, false);

    std::array<Fr, blst::fft_detail::domain_storage_len(n)> domain = {};
    blst::fft_detail::initialize_domain<Fr, blst::FrCoeffCtx>(domain, n);

    auto prepared = input;
    blst::fft_detail::transform_with_domain<Fr, Fr, blst::FrCoeffCtx>(
        prepared,
        domain,
        false);
    EXPECT_TRUE(equal(prepared, expected));

    blst::fft_detail::transform_with_domain<Fr, Fr, blst::FrCoeffCtx>(
        prepared,
        domain,
        true);
    EXPECT_TRUE(equal(prepared, input));
}

TEST(FftApi, CoefficientContextPreservesSelectedGroupOperations)
{
    const Fr scalar = from_u64(3);
    Fr scalar_non_montgomery = {};
    Fr::from_mont(scalar_non_montgomery, scalar);
    const auto* scalar_bytes =
        reinterpret_cast<const blst::byte*>(&scalar_non_montgomery);

    const blst::G1Point& g1 = blst::G1Point::generator();
    blst::G1Point expected_g1 = {};
    blst::G1Point::mult_vartime(expected_g1, g1, scalar_bytes, 256);
    const blst::G1Point actual_g1 = blst::G1CoeffCtx::mult(g1, scalar);
    EXPECT_TRUE(blst::G1Point::equal(actual_g1, expected_g1));

    const blst::G2Point& g2 = blst::G2Point::generator();
    blst::G2Point expected_g2 = {};
    blst::G2Point::mult_vartime(expected_g2, g2, scalar_bytes, 256);
    const blst::G2Point actual_g2 = blst::G2CoeffCtx::mult(g2, scalar);
    EXPECT_TRUE(blst::G2Point::equal(actual_g2, expected_g2));

    const blst::GT& gt = blst::GT::generator();
    blst::GT expected_gt = {};
    blst::GT::mult_vartime(expected_gt, gt, scalar_bytes, 256);
    const blst::GT actual_gt = blst::GTCoeffCtx::mult(gt, scalar);
    EXPECT_TRUE(blst::GT::equal(actual_gt, expected_gt));

    const blst::GT pairing = blst::CoeffCtx::mult(g1, g2);
    EXPECT_TRUE(blst::GT::equal(pairing, blst::GT::generator()));
}

TEST(FftCore, PlanConstructionRejectsInvalidAndOverflowedShapes)
{
    EXPECT_THROW(make_fft_truncated_forward_plan(0, 1), std::invalid_argument);
    EXPECT_THROW(make_fft_truncated_forward_plan(3, 1), std::invalid_argument);
    EXPECT_THROW(make_fft_truncated_forward_plan(8, 0), std::invalid_argument);
    EXPECT_THROW(make_fft_truncated_forward_plan(8, 9), std::invalid_argument);

    EXPECT_THROW(
        make_fft_inverse_truncated_range_plan(0, 0, 1),
        std::invalid_argument);
    EXPECT_THROW(
        make_fft_inverse_truncated_range_plan(8, 8, 1),
        std::invalid_argument);
    EXPECT_THROW(
        make_fft_inverse_truncated_range_plan(8, 7, 2),
        std::invalid_argument);
    EXPECT_THROW(
        make_fft_inverse_truncated_range_plan(
            8,
            7,
            std::numeric_limits<std::size_t>::max()),
        std::invalid_argument);
}

TEST(FftCore, ReusableTruncatedPlansMatchFullTransforms)
{
    constexpr std::size_t n = 16;
    const std::vector<Fr> roots =
        blst::fft_detail::make_roots<Fr, blst::FrCoeffCtx>(n, false);
    const std::vector<Fr> inverse_roots =
        blst::fft_detail::make_roots<Fr, blst::FrCoeffCtx>(n, true);

    for (const std::size_t xn : { std::size_t{1}, std::size_t{5}, n }) {
        std::array<Fr, n> expected = {};
        for (std::size_t i = 0; i < xn; ++i) {
            expected[i] = from_u64(17 * i + 3);
        }
        auto planned = expected;
        fft<Fr, Fr, blst::FrCoeffCtx>(expected, roots, false);

        const FftTruncatedForwardPlan plan =
            make_fft_truncated_forward_plan(n, xn);
        EXPECT_EQ(
            plan.butterflies_executed + plan.butterflies_skipped,
            plan.butterflies_total);
        blst::TimingStats timings;
        fft_forward_truncated_prefix<Fr, Fr, blst::FrCoeffCtx>(
            planned,
            roots,
            xn,
            plan,
            &timings,
            "fft_test");
        EXPECT_TRUE(equal(expected, planned));
        EXPECT_EQ(timings.get_ns("fft_test.xn"), xn);
    }

    std::array<Fr, n> frequency = {};
    for (std::size_t i = 0; i < n; ++i) {
        frequency[i] = from_u64(29 * i + 11);
    }
    auto full_inverse = frequency;
    fft<Fr, Fr, blst::FrCoeffCtx>(full_inverse, inverse_roots, true);

    for (const auto [start, length] : {
             std::pair<std::size_t, std::size_t>{0, 1},
             {3, 7},
             {15, 1},
             {0, n} }) {
        auto planned = frequency;
        const FftInverseTruncatedRangePlan plan =
            make_fft_inverse_truncated_range_plan(n, start, length);
        EXPECT_EQ(
            plan.butterflies_executed + plan.butterflies_skipped,
            plan.butterflies_total);
        fft_inverse_truncated_range<Fr, Fr, blst::FrCoeffCtx>(
            planned,
            inverse_roots,
            start,
            length,
            plan);
        EXPECT_TRUE(equal(
            std::span<const Fr>(planned).subspan(start, length),
            std::span<const Fr>(full_inverse).subspan(start, length)));
    }
}

#ifndef NDEBUG
TEST(FftCore, DebugBuildRejectsInvalidHotTransformShapes)
{
    std::array<Fr, 3> non_power_of_two = {};
    std::array<Fr, 3> roots = {};
    EXPECT_DEATH(
        (fft<Fr, Fr, blst::FrCoeffCtx>(non_power_of_two, roots, false)),
        "");

    std::array<Fr, 2> values = {};
    EXPECT_DEATH(
        (fft<Fr, Fr, blst::FrCoeffCtx>(
            values,
            std::span<const Fr>{},
            false)),
        "");

    std::array<Fr, 3> invalid_roots = {};
    EXPECT_DEATH(
        (blst::fft_detail::fill_roots<Fr, blst::FrCoeffCtx>(
            invalid_roots,
            false)),
        "");

    std::array<Fr, 2> short_domain = {};
    EXPECT_DEATH(
        (blst::fft_detail::initialize_domain<Fr, blst::FrCoeffCtx>(
            short_domain,
            1)),
        "");

    std::array<Fr, 2> prepared_values = {};
    std::array<Fr, 3> short_prepared_domain = {};
    EXPECT_DEATH(
        (blst::fft_detail::transform_with_domain<
            Fr,
            Fr,
            blst::FrCoeffCtx>(
                prepared_values,
                short_prepared_domain,
                false)),
        "");
}
#endif

TEST(FftApi, GeneralInterpolationReportsDuplicatePoints)
{
    const std::array<Fr, 3> coefficients = {
        from_u64(9), from_u64(2), from_u64(7)
    };
    const std::array<Fr, 3> points = {
        from_u64(1), from_u64(2), from_u64(3)
    };
    std::array<Fr, 3> values = {};
    for (std::size_t i = 0; i < points.size(); ++i) {
        Fr acc = {};
        for (auto coefficient = coefficients.rbegin();
             coefficient != coefficients.rend();
             ++coefficient) {
            acc = add(multiply(acc, points[i]), *coefficient);
        }
        values[i] = acc;
    }

    std::array<Fr, 3> recovered = {};
    ASSERT_EQ(
        blst::fr_interpolation::interpolate(recovered, points, values),
        InterpolationStatus::success);
    EXPECT_TRUE(equal(recovered, coefficients));

    auto duplicate_points = points;
    duplicate_points[1] = duplicate_points[0];
    const std::array<Fr, 3> sentinel = {
        from_u64(101), from_u64(102), from_u64(103)
    };
    recovered = sentinel;
    EXPECT_EQ(
        blst::fr_interpolation::interpolate(
            recovered,
            duplicate_points,
            values),
        InterpolationStatus::duplicate_point);
    EXPECT_TRUE(equal(recovered, sentinel));
}

TEST(FftApi, GeneralInterpolationCrossesFastPathThresholds)
{
    for (const std::size_t n : { 63U, 64U, 65U }) {
        SCOPED_TRACE(n);
        std::vector<Fr> coefficients(n);
        std::vector<Fr> points(n);
        std::vector<Fr> values(n);
        for (std::size_t i = 0; i < n; ++i) {
            coefficients[i] = (i % 7 == 3) ? Fr{} : from_u64(17 * i + 5);
            points[i] = from_u64(i + 1);
        }
        for (std::size_t i = 0; i < n; ++i) {
            values[i] = evaluate_polynomial(coefficients, points[i]);
        }

        std::vector<Fr> recovered(n, from_u64(0xBAD));
        ASSERT_EQ(
            blst::fr_interpolation::interpolate(recovered, points, values),
            InterpolationStatus::success);
        EXPECT_TRUE(equal(recovered, coefficients));

        std::fill(values.begin(), values.end(), Fr{});
        std::fill(recovered.begin(), recovered.end(), from_u64(0xBAD));
        ASSERT_EQ(
            blst::fr_interpolation::interpolate(recovered, points, values),
            InterpolationStatus::success);
        EXPECT_TRUE(std::all_of(
            recovered.begin(),
            recovered.end(),
            [](const Fr& value) { return Fr::equal(value, Fr{}); }));
    }

    constexpr std::size_t n = 65;
    std::vector<Fr> points(n);
    std::vector<Fr> values(n, Fr::one());
    std::vector<Fr> output(n, from_u64(0xBAD));
    for (std::size_t i = 0; i < n; ++i) {
        points[i] = from_u64(i + 1);
    }
    points.back() = points.front();
    const std::vector<Fr> sentinel = output;
    EXPECT_EQ(
        blst::fr_interpolation::interpolate(output, points, values),
        InterpolationStatus::duplicate_point);
    EXPECT_TRUE(equal(output, sentinel));
}

TEST(FftApi, DomainSuffixUsesCallerScratchAndReportsErrors)
{
    constexpr std::size_t fft_len = 8;
    constexpr std::size_t n = 5;
    constexpr std::size_t ell = 2;
    const std::array<Fr, n> coefficients = {
        from_u64(13),
        from_u64(21),
        from_u64(34),
        from_u64(55),
        from_u64(89),
    };
    const std::array<std::size_t, n> indices = { 0, 1, 3, 4, 6 };

    const std::vector<Fr> roots =
        blst::fft_detail::make_roots<Fr, blst::FrCoeffCtx>(fft_len, false);
    std::array<Fr, n> values = {};
    for (std::size_t i = 0; i < n; ++i) {
        Fr acc = {};
        for (auto coefficient = coefficients.rbegin();
             coefficient != coefficients.rend();
             ++coefficient) {
            acc = add(multiply(acc, roots[indices[i]]), *coefficient);
        }
        values[i] = acc;
    }

    const std::size_t required =
        blst::fr_interpolation::domain_suffix_scratch_len(n, fft_len);
    std::vector<Fr> scratch(required);
    std::array<Fr, n - ell> suffix = {};
    ASSERT_EQ(
        blst::fr_interpolation::interpolate_domain_suffix_with_scratch(
            suffix,
            indices,
            values,
            ell,
            fft_len,
            scratch),
        InterpolationStatus::success);
    for (std::size_t i = 0; i < suffix.size(); ++i) {
        EXPECT_TRUE(Fr::equal(suffix[i], coefficients[ell + i]));
    }

    EXPECT_EQ(
        blst::fr_interpolation::interpolate_domain_suffix_with_scratch(
            suffix,
            indices,
            values,
            ell,
            fft_len,
            std::span<Fr>(scratch).first(required - 1)),
        InterpolationStatus::insufficient_scratch);
    EXPECT_TRUE(equal(
        std::span<const Fr>(suffix),
        std::span<const Fr>(coefficients).subspan(ell)));

    auto values_alias = values;
    ASSERT_EQ(
        blst::fr_interpolation::interpolate_domain_suffix_with_scratch(
            std::span<Fr>(values_alias).first(n - ell),
            indices,
            values_alias,
            ell,
            fft_len,
            scratch),
        InterpolationStatus::success);
    EXPECT_TRUE(equal(
        std::span<const Fr>(values_alias).first(n - ell),
        std::span<const Fr>(coefficients).subspan(ell)));

    auto invalid_indices = indices;
    invalid_indices[0] = fft_len;
    const auto sentinel = suffix;
    EXPECT_EQ(
        blst::fr_interpolation::interpolate_domain_suffix_with_scratch(
            suffix,
            invalid_indices,
            values,
            ell,
            fft_len,
            scratch),
        InterpolationStatus::invalid_index);
    EXPECT_TRUE(equal(suffix, sentinel));

    auto duplicate_indices = indices;
    duplicate_indices[1] = duplicate_indices[0];
    EXPECT_EQ(
        blst::fr_interpolation::interpolate_domain_suffix_with_scratch(
            suffix,
            duplicate_indices,
            values,
            ell,
            fft_len,
            scratch),
        InterpolationStatus::duplicate_point);
    EXPECT_TRUE(equal(suffix, sentinel));
}
