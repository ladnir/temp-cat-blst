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

        friend constexpr bool operator==(const RootTestField, const RootTestField) = default;
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
        blst::fr_interpolation::domain_suffix_scratch_len(5, 8) == 41);
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
}

TEST(FftApi, RootPreparationRejectsEmptyDomain)
{
    const RootTestField one(1);
    EXPECT_THROW(
        (blst::precomputeRootsOfUnity<RootTestField>(
            one,
            std::span<RootTestField>{})),
        std::runtime_error);

    std::array<RootTestField, 1> roots = {};
    EXPECT_NO_THROW(blst::precomputeRootsOfUnity<RootTestField>(one, roots));
    EXPECT_EQ(roots[0], one);
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
}
