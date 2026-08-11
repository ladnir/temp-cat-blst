#include "blst/avx/generated_matrixnok_constants.hpp"
#include "blst/avx/rns/matrixnok_change_base.h"
#include "blst/avx/rns/matrixnok_fp.h"
#include "blst/avx/fields/matrixnok_fp12.h"
#include "blst/avx/fields/final_exponentiation.h"
#include "blst/avx/gt/glv4.h"
#include "blst/avx/curve/matrixnok_g1.h"
#include "blst/avx/curve/matrixnok_g2.h"
#include "blst/avx/curve/point.h"
#include "blst/avx/curve/GT.h"
#include "blst/avx/pairing/matrixnok_pairing.h"
#include "blst/common/bytes.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/primitveRoot.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fr.h"
#include "blst/scalar/curve/point.h"
#include "blst/scalar/curve/GT.h"

#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

namespace blst::avx::test
{
    namespace matrixnok_constants = matrixnok::generated;

    template<std::size_t lanes>
    struct ScalarRnsAccumulators
    {
        std::array<uint64_t, lanes> high;
        std::array<uint64_t, lanes> low;
    };

    static constexpr uint64_t kIfmaDigitMask =
        (uint64_t{1} << 52) - 1;

    static constexpr uint64_t inverse_power_of_two(
        uint64_t modulus, unsigned int exponent)
    {
        uint64_t result = 1;
        for (unsigned int bit = 0; bit < exponent; ++bit) {
            result = (result & 1) != 0
                ? (result + modulus) >> 1
                : result >> 1;
        }
        return result;
    }

    // Scalar model of VPMADD52: multiply the low 52-bit digits and add the
    // selected half with ordinary uint64_t accumulator wrap.
    static void scalar_ifma_accumulate(
        uint64_t& high, uint64_t& low, uint64_t left, uint64_t right)
    {
        const uint128_t product =
            static_cast<uint128_t>(left & kIfmaDigitMask) *
            static_cast<uint128_t>(right & kIfmaDigitMask);
        low += static_cast<uint64_t>(product) & kIfmaDigitMask;
        high += static_cast<uint64_t>(product >> 52) & kIfmaDigitMask;
    }

    template<bool shifted, std::size_t input_count, std::size_t output_count>
    static ScalarRnsAccumulators<output_count> scalar_rns_accumulate(
        const std::array<std::array<uint64_t, output_count>, input_count>& matrix,
        const std::array<uint64_t, output_count>& correction,
        const std::array<uint64_t, input_count>& quotient_estimates,
        const std::array<uint64_t, input_count>& residues,
        const std::array<uint64_t, output_count>& initial_high,
        const std::array<uint64_t, output_count>& initial_low,
        const std::array<uint64_t, output_count>& output_moduli)
    {
        // This is deliberately independent of FixedRNSMatrix and the vector
        // wrappers. It reconstructs the matrix/correction schedule directly
        // from scalar products and the generated constants.
        ScalarRnsAccumulators<output_count> result{initial_high, initial_low};
        uint128_t quotient_numerator = 0;
        for (std::size_t input = 0; input < input_count; ++input) {
            quotient_numerator += static_cast<uint128_t>(residues[input]) *
                static_cast<uint128_t>(quotient_estimates[input]);
            for (std::size_t output = 0; output < output_count; ++output) {
                scalar_ifma_accumulate(
                    result.high[output], result.low[output],
                    matrix[input][output], residues[input]);
            }
        }

        const uint64_t quotient =
            static_cast<uint64_t>(quotient_numerator >> 64);
        const uint64_t quotient_high = quotient >> 52;
        for (std::size_t output = 0; output < output_count; ++output) {
            scalar_ifma_accumulate(
                result.high[output], result.low[output],
                correction[output], quotient);
            if constexpr (shifted) {
                const uint64_t shifted_correction = mod_mul(
                    correction[output], uint64_t{1} << 52,
                    output_moduli[output]);
                scalar_ifma_accumulate(
                    result.high[output], result.low[output],
                    shifted_correction, quotient_high);
            } else {
                const uint128_t high_product =
                    static_cast<uint128_t>(correction[output]) *
                    static_cast<uint128_t>(quotient_high);
                result.high[output] +=
                    static_cast<uint64_t>(high_product) & kIfmaDigitMask;
                result.high[output] +=
                    (static_cast<uint64_t>(high_product >> 52) &
                        kIfmaDigitMask) << 52;
            }
        }
        return result;
    }

    static_assert(matrixnok_constants::BLS12_381_RNS_DIM == 8);
    static_assert(matrixnok_constants::BLS12_381_CONVERT_TO_RADIX_BITS == 50);
    static_assert(matrixnok_constants::moduli1.size() == 8);
    static_assert(matrixnok_constants::moduli2.size() == 8);
    static_assert(matrixnok_constants::r1_perm_mat_fixed.size() == 8);
    static_assert(sizeof(uint128_t) == 16);
    static_assert(sizeof(int128_t) == 16);
    static constexpr int128_t kSignedAboveInt64 =
        int128_t(INT64_MAX) + int128_t(int64_t{1});
    static constexpr int128_t kSignedBelowInt64 =
        int128_t(INT64_MIN) - int128_t(int64_t{1});
    static constexpr int128_t kSignedInt64Square =
        int128_t(INT64_MAX) * int128_t(INT64_MAX);
    static constexpr uint128_t kUnsignedInt64Square =
        uint128_t(UINT64_MAX) * uint128_t(UINT64_MAX);
    static constexpr uint128_t kUnsignedWideProduct =
        (uint128_t(uint64_t{1}) << 64 | uint128_t(uint64_t{3})) *
        (uint128_t(uint64_t{1}) << 64 | uint128_t(uint64_t{5}));
    static_assert(kSignedAboveInt64 > int128_t(INT64_MAX));
    static_assert(kSignedBelowInt64 < int128_t(INT64_MIN));
    static_assert(kSignedInt64Square > int128_t(INT64_MAX));
    static_assert(static_cast<uint64_t>(kUnsignedInt64Square) == 1);
    static_assert(
        static_cast<uint64_t>(kUnsignedInt64Square >> 64) ==
        UINT64_MAX - 1);
    static_assert(static_cast<uint64_t>(kUnsignedWideProduct) == 15);
    static_assert(static_cast<uint64_t>(kUnsignedWideProduct >> 64) == 8);
    static_assert(
        int128_t(int64_t{-7}) * int128_t(int64_t{9}) ==
        int128_t(int64_t{-63}));
    using ShiftedSignedBounds = decltype(
        Bounds<-2, 3>{}.template shift_left<2>());
    static_assert(ShiftedSignedBounds::lower == -8);
    static_assert(ShiftedSignedBounds::upper == 12);
    using LargestSelectedWide = WideElement<Bounds<0, 832>, 8, 50>;
    static_assert(LargestSelectedWide::max_wide_bound == 2048);
    using BorrowFreeWideLeft = WideElement<Bounds<80, 81>, 8, 50>;
    using BorrowFreeWideRight = WideElement<Bounds<0, 40>, 8, 50>;
    using BorrowFreeWideDifference = decltype(
        std::declval<BorrowFreeWideLeft>() -
        std::declval<BorrowFreeWideRight>());
    static_assert(BorrowFreeWideDifference::bounds.lower == 40);
    static_assert(BorrowFreeWideDifference::bounds.upper == 81);
    using RingCarrierResidue = BoundedElement<Bounds<0, 1>, 8, 50>;
    using RingSmall = SmallElement<RingCarrierResidue, Bounds<-1, 1>>;
    using RingSmallWideRecast = decltype(
        std::declval<RingSmall>().template recast_rns_bounds<
            INT64_C(-2147483649), INT64_C(2147483648)>());
    static_assert(RingSmallWideRecast::rns_bounds.lower ==
        INT64_C(-2147483649));
    static_assert(RingSmallWideRecast::rns_bounds.upper ==
        INT64_C(2147483648));
    using RingExpanded = ExpandedElement<RingCarrierResidue, Bounds<0, 1>>;
    using RingCompleted = decltype(std::declval<RingExpanded>().complete());
    static_assert(std::is_same_v<
        typename RingCompleted::rns_bounds_type, Bounds<0, 1>>);
    static_assert(!std::is_const_v<typename RingCompleted::rns_bounds_type>);
    using RingDelayedProduct = decltype(
        std::declval<RingExpanded>() * std::declval<RingExpanded>());
    using RingCompletedProduct = decltype(
        std::declval<RingDelayedProduct>().complete());
    static_assert(RingCompletedProduct::rns_bounds.lower == 0);
    static_assert(RingCompletedProduct::rns_bounds.upper == 1);
    using RingAccumulatedProduct = decltype(
        std::declval<RingCompletedProduct>() +
        std::declval<RingDelayedProduct>());
    static_assert(RingAccumulatedProduct::rns_bounds.lower == 0);
    static_assert(RingAccumulatedProduct::rns_bounds.upper == 2);
    static_assert(std::is_empty_v<matrixnok::FpRingCandidate>);
    static_assert(std::is_empty_v<matrixnok::Fp12Ring>);
    static_assert(std::is_empty_v<matrixnok::G1CurveCandidate>);
    static_assert(std::is_empty_v<matrixnok::G2CurveCandidate>);
    static_assert(std::is_trivially_default_constructible_v<
        matrixnok::PairingInverter>);
    static_assert(std::is_same_v<avx_curve_detail::CurveRing,
        matrixnok::FpRingCandidate>);
    static_assert(std::is_same_v<avx_curve_detail::CurveFp2Ring,
        matrixnok::G2Fp2RingCandidate>);
    static_assert(std::is_same_v<avx_curve_detail::CurveG1,
        matrixnok::G1CurveCandidate>);
    static_assert(std::is_same_v<avx_curve_detail::CurveG2,
        matrixnok::G2CurveCandidate>);

    using MatrixNoKFp12 = matrixnok::Fp12Ring;
    static_assert(std::is_default_constructible_v<
        typename MatrixNoKFp12::StandardElement>);

    TEST(AvxMatrixNoKPrimitives, RingCarrierContractsAreExact)
    {
        std::array<uint64_t, 8> left_lanes{};
        std::array<uint64_t, 8> right_lanes{};
        for (std::size_t lane = 0; lane < left_lanes.size(); ++lane) {
            left_lanes[lane] = lane + 1;
            right_lanes[lane] = 3 * lane + 2;
        }

        const RingExpanded left{
            RingCarrierResidue(AVXVector<8>(left_lanes)),
            RingCarrierResidue(AVXVector<8>(left_lanes))};
        const RingExpanded right{
            RingCarrierResidue(AVXVector<8>(right_lanes)),
            RingCarrierResidue(AVXVector<8>(right_lanes))};

        EXPECT_EQ(left.ctime_select(true, right).m1.data.to_array(), left_lanes);
        EXPECT_EQ(left.ctime_select(false, right).m1.data.to_array(), right_lanes);
        EXPECT_EQ(left.ctime_replace_if(true, right).m2.data.to_array(),
            right_lanes);
        EXPECT_EQ(left.ctime_replace_if(false, right).m2.data.to_array(),
            left_lanes);

        const RingExpanded value_initialized{};
        EXPECT_EQ(value_initialized.m1.data.to_array(),
            (std::array<uint64_t, 8>{}));
        EXPECT_EQ(value_initialized.m2.data.to_array(),
            (std::array<uint64_t, 8>{}));
        EXPECT_EQ(left.is_zero(), 0);
        EXPECT_EQ(value_initialized.is_zero(), 1);
        EXPECT_EQ(left.equal(right), 0);
        EXPECT_EQ(left.equal(left), 1);

        const auto product = (left * right).complete();
        const auto accumulated = product + (left * right);
        const auto product_high = product.m1.high.data.to_array();
        const auto product_low = product.m1.low.data.to_array();
        const auto accumulated_high = accumulated.m1.high.data.to_array();
        const auto accumulated_low = accumulated.m1.low.data.to_array();
        for (std::size_t lane = 0; lane < left_lanes.size(); ++lane) {
            const uint64_t expected = left_lanes[lane] * right_lanes[lane];
            EXPECT_EQ(product_high[lane], expected >> 52);
            EXPECT_EQ(product_low[lane], expected & kIfmaDigitMask);
            EXPECT_EQ(accumulated_high[lane], (2 * expected) >> 52);
            EXPECT_EQ(accumulated_low[lane], (2 * expected) & kIfmaDigitMask);
        }
    }

    TEST(AvxMatrixNoKPrimitives, WideIntegerCompatibilityIsExact)
    {
        const auto explicit_zero =
            BoundedElement<Bounds<0, 1>, 8, 50>::zero().data.to_array();
        EXPECT_EQ(explicit_zero, (std::array<uint64_t, 8>{}));

        const uint128_t full_product =
            uint128_t(UINT64_MAX) * uint128_t(UINT64_MAX);
        EXPECT_EQ(static_cast<uint64_t>(full_product), uint64_t{1});
        EXPECT_EQ(
            static_cast<uint64_t>(full_product >> 64), UINT64_MAX - 1);

        const uint128_t shifted =
            (uint128_t(uint64_t{1}) << 64) + uint128_t(uint64_t{5});
        EXPECT_EQ(shifted % uint64_t{17}, uint64_t{6});
#ifdef NDEBUG
        EXPECT_EQ(mod_mul(1, 1, 0), 0);
#else
        EXPECT_DEATH((void)mod_mul(1, 1, 0), "");
#endif

        const int128_t positive_wide =
            int128_t(INT64_MAX) * int128_t(int64_t{2});
        const int128_t negative_wide =
            int128_t(INT64_MIN) * int128_t(int64_t{2});
        EXPECT_TRUE(positive_wide > int128_t(INT64_MAX));
        EXPECT_TRUE(negative_wide < int128_t(INT64_MIN));
        EXPECT_EQ(
            static_cast<int64_t>(
                int128_t(int64_t{-1234567}) * int128_t(int64_t{17})),
            int64_t{-20987639});

#if defined(_MSC_VER) && !defined(__clang__)
        constexpr std::array<int64_t, 8> signed_values = {
            INT64_MIN, INT64_MIN + 1, -123456789, -1,
            0, 1, 123456789, INT64_MAX};
        for (const int64_t left : signed_values) {
            for (const int64_t right : signed_values) {
                int64_t expected_high = 0;
                const int64_t expected_low = _mul128(
                    left, right, &expected_high);
                const int128_t actual = int128_t(left) * int128_t(right);
                EXPECT_EQ(actual.lo, std::bit_cast<uint64_t>(expected_low));
                EXPECT_EQ(actual.hi, std::bit_cast<uint64_t>(expected_high));
            }
        }
#endif
    }

    TEST(AvxMatrixNoKPrimitives, PermutationRowsAreCacheLineAligned)
    {
        const auto expect_aligned = [](const auto& values) {
            const auto address = reinterpret_cast<std::uintptr_t>(values.data());
            EXPECT_EQ(address % 64, 0U);
        };
        expect_aligned(matrixnok_constants::moduli1);
        expect_aligned(matrixnok_constants::moduli2);
        expect_aligned(matrixnok_constants::r1_perm_mat_fixed);
        expect_aligned(matrixnok_constants::r2_rns_mat);
        expect_aligned(matrixnok_constants::r2_correction);
        expect_aligned(matrixnok_constants::convert_to_rns_mat);
        expect_aligned(matrixnok_constants::convert_to_correction);
        expect_aligned(matrixnok_constants::convert_from_rns_mat);
        expect_aligned(matrixnok_constants::convert_from_correction);
        expect_aligned(matrixnok::detail::reducer_m1_constants.modulus_multiples);
        expect_aligned(matrixnok::detail::reducer_m1_constants.t);
        expect_aligned(matrixnok::detail::reducer_m1_constants.negative_moduli);
        expect_aligned(matrixnok::detail::reducer_m1_constants.montgomery_factor);
        expect_aligned(matrixnok::detail::reducer_m2_constants.modulus_multiples);
        expect_aligned(matrixnok::detail::reducer_m2_constants.t);
        expect_aligned(matrixnok::detail::reducer_m2_constants.negative_moduli);
        expect_aligned(matrixnok::detail::reducer_m2_constants.montgomery_factor);
    }

    TEST(AvxMatrixNoKPrimitives, FixedContextMatchesScalarAccumulatorOracle)
    {
        EXPECT_EQ(sizeof(matrixnok::FpRingCandidate), 1U);
        EXPECT_EQ(sizeof(matrixnok::G2Fp2RingCandidate), 1U);

        std::array<uint64_t, 8> r2_input = {};
        std::array<uint64_t, 8> digit_input = {};
        std::array<uint64_t, 8> initial_high = {};
        std::array<uint64_t, 8> initial_low = {};
        for (std::size_t lane = 0; lane < 8; ++lane) {
            r2_input[lane] = matrixnok_constants::moduli2[lane] - lane - 1;
            digit_input[lane] = (0x123456789abULL * (lane + 1)) &
                ((uint64_t{1} << 50) - 1);
            initial_high[lane] = 17 * lane + 3;
            initial_low[lane] = 29 * lane + 5;
        }

        const auto expected_r2 = scalar_rns_accumulate<true>(
            matrixnok_constants::r2_rns_mat,
            matrixnok_constants::r2_correction,
            matrixnok_constants::r2_shifted_quotient_estimations,
            r2_input, initial_high, initial_low,
            matrixnok_constants::moduli1);
        AVXVector<8> fixed_r2_high(initial_high);
        AVXVector<8> fixed_r2_low(initial_low);
        matrixnok::detail::R2Matrix::reduce_raw_shifted(
            AVXVector<8>(r2_input), fixed_r2_high, fixed_r2_low);
        EXPECT_EQ(fixed_r2_high.to_array(), expected_r2.high);
        EXPECT_EQ(fixed_r2_low.to_array(), expected_r2.low);

        const auto expected_convert_to = scalar_rns_accumulate<true>(
            matrixnok_constants::convert_to_rns_mat,
            matrixnok_constants::convert_to_correction,
            matrixnok_constants::convert_to_shifted_quotient_estimations,
            digit_input, initial_high, initial_low,
            matrixnok_constants::moduli2);
        AVXVector<8> fixed_to_high(initial_high);
        AVXVector<8> fixed_to_low(initial_low);
        matrixnok::detail::ConvertToMatrix::reduce_raw_shifted(
            AVXVector<8>(digit_input), fixed_to_high, fixed_to_low);
        EXPECT_EQ(fixed_to_high.to_array(), expected_convert_to.high);
        EXPECT_EQ(fixed_to_low.to_array(), expected_convert_to.low);

        const auto expected_convert_from = scalar_rns_accumulate<false>(
            matrixnok_constants::convert_from_rns_mat,
            matrixnok_constants::convert_from_correction,
            matrixnok_constants::convert_from_shifted_quotient_estimations,
            r2_input, initial_high, initial_low,
            std::array<uint64_t, 8>{});
        AVXVector<8> fixed_from_high(initial_high);
        AVXVector<8> fixed_from_low(initial_low);
        matrixnok::detail::ConvertFromMatrix::reduce_raw(
            AVXVector<8>(r2_input), fixed_from_high, fixed_from_low);
        EXPECT_EQ(fixed_from_high.to_array(), expected_convert_from.high);
        EXPECT_EQ(fixed_from_low.to_array(), expected_convert_from.low);

        std::array<AVXVector<8>, 2> batch_input = {
            AVXVector<8>(r2_input), AVXVector<8>(digit_input)};
        std::array<AVXVector<8>, 2> initial_batch_high = {
            AVXVector<8>(initial_high), AVXVector<8>(initial_low)};
        std::array<AVXVector<8>, 2> initial_batch_low = {
            AVXVector<8>(initial_low), AVXVector<8>(initial_high)};
        auto fixed_batch_high = initial_batch_high;
        auto fixed_batch_low = initial_batch_low;
        matrixnok::detail::R2Matrix::reduce_raw_batch_shifted<2>(
            batch_input, fixed_batch_high, fixed_batch_low);
        for (std::size_t item = 0; item < 2; ++item) {
            const auto expected = scalar_rns_accumulate<true>(
                matrixnok_constants::r2_rns_mat,
                matrixnok_constants::r2_correction,
                matrixnok_constants::r2_shifted_quotient_estimations,
                batch_input[item].to_array(),
                initial_batch_high[item].to_array(),
                initial_batch_low[item].to_array(),
                matrixnok_constants::moduli1);
            EXPECT_EQ(fixed_batch_high[item].to_array(), expected.high);
            EXPECT_EQ(fixed_batch_low[item].to_array(), expected.low);
        }

        using WideBound = Bounds<0, 24>;
        using WideResidue = BoundedElement<WideBound, 8, 50>;
        std::array<uint64_t, 8> wide_high = {};
        std::array<uint64_t, 8> wide_low = {};
        for (std::size_t lane = 0; lane < 8; ++lane) {
            wide_high[lane] = 7 * matrixnok_constants::moduli1[lane] + lane;
            wide_low[lane] = ((uint64_t{1} << 52) - 1) - lane;
        }
        const WideElement<WideBound, 8, 50> wide{
            WideResidue(AVXVector<8>(wide_high)),
            WideResidue(AVXVector<8>(wide_low))};
        const auto fixed_reduced = matrixnok::detail::ReducerM1::template reduce_auto<2>(
            matrixnok::detail::ReducerM1::reduce_wide(wide));
        const auto reduced_lanes = fixed_reduced.data.to_array();
        const auto doubled_moduli =
            matrixnok::detail::ReducerM1::template moduli<1>().data.to_array();
        for (std::size_t lane = 0; lane < 8; ++lane) {
            const uint64_t modulus = matrixnok_constants::moduli1[lane];
            const uint64_t radix_inverse = inverse_power_of_two(modulus, 52);
            const uint64_t expected =
                (wide_high[lane] % modulus +
                    mod_mul(wide_low[lane], radix_inverse, modulus)) % modulus;
            EXPECT_LT(reduced_lanes[lane], 2 * modulus);
            EXPECT_EQ(reduced_lanes[lane] % modulus, expected);
            EXPECT_EQ(doubled_moduli[lane], 2 * modulus);
        }
    }

    TEST(AvxMatrixNoKPrimitives, LowLaneExtractionAndBroadcast)
    {
        const std::array<uint64_t, 8> low = {
            (uint64_t{7} << 52) | 5, 11, 13, 17, 19, 23, 29, 31};
        const std::array<uint64_t, 8> high = {
            3, 37, 41, 43, 47, 53, 59, 61};

        const ShortVector low_short = ShortVector::from_avx_lo(AVXVector<8>(low));
        const ShortVector high_short = ShortVector::from_avx_lo(AVXVector<8>(high));
        const auto broadcast = high_short.add(low_short.srli(52)).broadcast_to_avx<8>().to_array();

        for (const uint64_t lane : broadcast) {
            EXPECT_EQ(lane, 10);
        }
        const AVXVector<8> vector(low);
        EXPECT_EQ(vector.extract<0>(), low[0]);
        EXPECT_EQ(vector.extract<1>(), low[1]);
        EXPECT_EQ(vector.extract<2>(), low[2]);
        EXPECT_EQ(vector.extract<3>(), low[3]);
        EXPECT_EQ(vector.extract<4>(), low[4]);
        EXPECT_EQ(vector.extract<5>(), low[5]);
        EXPECT_EQ(vector.extract<6>(), low[6]);
        EXPECT_EQ(vector.extract<7>(), low[7]);
    }

    template<int64_t upper_bound, int64_t stop>
    void check_min_sub_schedule()
    {
        constexpr int limbs = 8;
        constexpr int element_bits = 50;
        const auto& moduli = matrixnok_constants::moduli1;

        constexpr std::array<int64_t, 16> quotient_seeds = {
            0, 1, 2, 3, 4, 5, 7, 8, 15, 16, 31, 32,
            upper_bound / 2, upper_bound - 3, upper_bound - 2, upper_bound - 1,
        };
        for (const int64_t seed : quotient_seeds) {
            const int64_t quotient = seed < 0 ? 0 : (seed >= upper_bound ? upper_bound - 1 : seed);
            std::array<uint64_t, limbs> values{};
            for (int lane = 0; lane < limbs; ++lane) {
                const uint64_t remainder = (lane & 1) == 0 ? 0 : moduli[lane] - 1;
                values[lane] = static_cast<uint64_t>(quotient) * moduli[lane] + remainder;
            }

            using Input = BoundedElement<Bounds<0, upper_bound>, limbs, element_bits>;
            const auto reduced = matrixnok::detail::ReducerM1::template min_sub_reduce<
                upper_bound, stop>(Input(AVXVector<limbs>(values)));
            const auto actual = reduced.data.to_array();
            for (int lane = 0; lane < limbs; ++lane) {
                EXPECT_LT(actual[lane], static_cast<uint64_t>(stop) * moduli[lane]);
                EXPECT_EQ(actual[lane] % moduli[lane], values[lane] % moduli[lane]);
            }
        }
    }

    TEST(AvxMatrixNoKPrimitives, ConditionalSubtractionCoversConfiguredGap)
    {
        check_min_sub_schedule<2, 1>();
        check_min_sub_schedule<4, 1>();
        check_min_sub_schedule<40, 1>();
        check_min_sub_schedule<100, 95>();
    }

    struct Product52
    {
        uint64_t low;
        uint64_t high;
    };

    constexpr Product52 split_product_52(uint64_t left, uint64_t right)
    {
        constexpr uint64_t half_mask = (uint64_t{1} << 26) - 1;
        constexpr uint64_t limb_mask = (uint64_t{1} << 52) - 1;
        const uint64_t left_low = left & half_mask;
        const uint64_t left_high = left >> 26;
        const uint64_t right_low = right & half_mask;
        const uint64_t right_high = right >> 26;
        const uint64_t low_product = left_low * right_low;
        const uint64_t cross = left_low * right_high + left_high * right_low;
        const uint64_t lower = low_product + ((cross & half_mask) << 26);
        return {lower & limb_mask, left_high * right_high + (cross >> 26) + (lower >> 52)};
    }

    template<std::size_t batch>
    void check_matrixnok_change_base(
        const std::array<std::array<uint64_t, 8>, batch>& input,
        const std::array<std::array<uint64_t, 8>, batch>& initial_hi,
        const std::array<std::array<uint64_t, 8>, batch>& initial_lo)
    {
        std::array<AVXVector<8>, batch> residues;
        std::array<AVXVector<8>, batch> actual_hi;
        std::array<AVXVector<8>, batch> actual_lo;
        for (std::size_t item = 0; item < batch; ++item) {
            residues[item] = AVXVector<8>(input[item]);
            actual_hi[item] = AVXVector<8>(initial_hi[item]);
            actual_lo[item] = AVXVector<8>(initial_lo[item]);
        }

        matrixnok::r1_change_base_raw_batch(residues, actual_hi, actual_lo);
        for (std::size_t item = 0; item < batch; ++item) {
            auto expected_hi = initial_hi[item];
            auto expected_lo = initial_lo[item];
            for (std::size_t input_lane = 0; input_lane < 8; ++input_lane) {
                for (std::size_t output_lane = 0; output_lane < 8; ++output_lane) {
                    const Product52 product = split_product_52(
                        input[item][input_lane],
                        matrixnok_constants::r1_rns_mat[input_lane][output_lane]);
                    expected_hi[output_lane] += product.high;
                    expected_lo[output_lane] += product.low;
                }
            }
            EXPECT_EQ(actual_hi[item].to_array(), expected_hi);
            EXPECT_EQ(actual_lo[item].to_array(), expected_lo);
        }
    }

    TEST(AvxMatrixNoKPrimitives, FixedR1PermutationMatchesScalarMatrix)
    {
        constexpr std::size_t batch = 3;
        std::array<std::array<uint64_t, 8>, batch> input{};
        std::array<std::array<uint64_t, 8>, batch> initial_hi{};
        std::array<std::array<uint64_t, 8>, batch> initial_lo{};

        for (std::size_t lane = 0; lane < 8; ++lane) {
            input[0][lane] = matrixnok_constants::moduli1[lane] - 1;
            input[1][lane] = lane;
            initial_hi[1][lane] = UINT64_MAX - lane;
            initial_lo[1][lane] = UINT64_MAX - 2 * lane;
        }

        uint64_t state = 0x9e3779b97f4a7c15ULL;
        for (std::size_t lane = 0; lane < 8; ++lane) {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            input[2][lane] = state % matrixnok_constants::moduli1[lane];
            initial_hi[2][lane] = state ^ 0xa5a5a5a5a5a5a5a5ULL;
            initial_lo[2][lane] = state ^ 0x5a5a5a5a5a5a5a5aULL;
        }

        check_matrixnok_change_base(input, initial_hi, initial_lo);
    }

    std::array<uint64_t, 8> canonical_50bit_digits(const blst::Fp& montgomery)
    {
        blst::Fp canonical = {};
        blst::cc_from_mont_n<blst::Fp::num_limb()>(
            canonical.l, montgomery.l, blst::get_BLST_P().l, blst::Fp::montgomery_constant());
        constexpr uint64_t mask = (uint64_t{1} << 50) - 1;
        std::array<uint64_t, 8> digits{};
        for (std::size_t digit = 0; digit < digits.size(); ++digit) {
            const std::size_t bit = digit * 50;
            const std::size_t limb = bit / 64;
            const unsigned shift = static_cast<unsigned>(bit % 64);
            uint64_t value = canonical.l[limb] >> shift;
            if (shift != 0 && limb + 1 < blst::Fp::num_limb()) {
                value |= canonical.l[limb + 1] << (64 - shift);
            }
            digits[digit] = value & mask;
        }
        return digits;
    }

    blst::Fp canonical_to_montgomery(std::array<uint64_t, blst::Fp::num_limb()> limbs)
    {
        blst::Fp canonical = {};
        for (std::size_t index = 0; index < limbs.size(); ++index) {
            canonical.l[index] = limbs[index];
        }
        blst::Fp montgomery = {};
        blst::cc_mul_mont_n<blst::Fp::num_limb()>(
            montgomery.l,
            canonical.l,
            blst::get_BLST_RR().l,
            blst::get_BLST_P().l,
            blst::Fp::montgomery_constant());
        return montgomery;
    }

    blst::Fp matrixnok_to_scalar(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::FpRingCandidate::Element& value)
    {
        AVXVector<8> high;
        AVXVector<8> low;
        ring.convert_from_raw(value, high, low);

        const auto low_digits = low.to_array();
        const auto high_digits = high.to_array();
        std::array<uint64_t, 9> digits = {};
        digits[0] = low_digits[0];
        for (std::size_t index = 1; index < 8; ++index) {
            digits[index] = high_digits[index - 1] + low_digits[index];
        }
        digits[8] = high_digits[7];

        // Independently decode the base-2^52 accumulator with the portable
        // scalar backend. Keeping this oracle out of the configured backend
        // dispatch lets the same small test run in portable-only and
        // native-plus-portable configurations without linking the ASM product.
        blst::UFp even = {};
        blst::UFp odd = {};
        auto* even_bytes = reinterpret_cast<unsigned char*>(even.data);
        auto* odd_bytes = reinterpret_cast<unsigned char*>(odd.data);
        for (std::size_t index = 0; index < digits.size(); index += 2) {
            std::memcpy(even_bytes + (index / 2) * 13, &digits[index], sizeof(uint64_t));
        }
        for (std::size_t index = 1; index < 8; index += 2) {
            const uint64_t shifted = digits[index] << 4;
            std::memcpy(odd_bytes + (index / 2) * 13 + 6, &shifted, sizeof(uint64_t));
        }

        blst::UFp wide = {};
        uint64_t carry = 0;
        for (std::size_t index = 0; index < std::size(wide.data); ++index) {
            const uint64_t left = even.data[index];
            const uint64_t right = odd.data[index];
            const uint64_t sum = left + right + carry;
            wide.data[index] = sum;
            carry = (sum < left) || (carry != 0 && sum == left);
        }

        blst::Fp reduced = {};
        blst::cc_redc_mont_n<blst::Fp::num_limb()>(
            reduced.l,
            wide.data,
            blst::get_BLST_P().l,
            blst::Fp::montgomery_constant());
        blst::Fp canonical = {};
        blst::cc_mul_mont_n<blst::Fp::num_limb()>(
            canonical.l,
            reduced.l,
            blst::get_BLST_RR().l,
            blst::get_BLST_P().l,
            blst::Fp::montgomery_constant());
        blst::Fp montgomery = {};
        blst::cc_mul_mont_n<blst::Fp::num_limb()>(
            montgomery.l,
            canonical.l,
            blst::get_BLST_RR().l,
            blst::get_BLST_P().l,
            blst::Fp::montgomery_constant());
        return montgomery;
    }

    void expect_same_scalar_fp(const blst::Fp& actual, const blst::Fp& expected)
    {
        for (std::size_t limb = 0; limb < blst::Fp::num_limb(); ++limb) {
            EXPECT_EQ(actual.l[limb], expected.l[limb]) << "limb " << limb;
        }
    }

    void expect_standard_residue_bounds(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::FpRingCandidate::Element& value)
    {
        const auto m1 = value.m1.data.to_array();
        const auto m2 = value.m2.data.to_array();
        for (std::size_t lane = 0; lane < m1.size(); ++lane) {
            EXPECT_LT(m1[lane], matrixnok::FpRingCandidate::standard_bound * ring.moduli1()[lane]);
            EXPECT_LT(m2[lane], matrixnok::FpRingCandidate::standard_bound * ring.moduli2()[lane]);
        }
    }

    blst::Fp scalar_product(const blst::Fp& left, const blst::Fp& right)
    {
        blst::Fp product = {};
        blst::cc_mul_mont_n<blst::Fp::num_limb()>(
            product.l,
            left.l,
            right.l,
            blst::get_BLST_P().l,
            blst::Fp::montgomery_constant());
        return product;
    }

    std::array<blst::Fp, 8> matrixnok_scalar_values()
    {
        auto p_minus_one = std::to_array<uint64_t>({
            blst::get_BLST_P().l[0], blst::get_BLST_P().l[1], blst::get_BLST_P().l[2],
            blst::get_BLST_P().l[3], blst::get_BLST_P().l[4], blst::get_BLST_P().l[5]});
        --p_minus_one[0];
        return {
            canonical_to_montgomery({0, 0, 0, 0, 0, 0}),
            canonical_to_montgomery({1, 0, 0, 0, 0, 0}),
            canonical_to_montgomery({2, 0, 0, 0, 0, 0}),
            canonical_to_montgomery(p_minus_one),
            canonical_to_montgomery({0x0123456789abcdefULL, 0xfedcba9876543210ULL, 0, 0, 0, 0}),
            canonical_to_montgomery({UINT64_MAX, UINT64_MAX, UINT64_MAX, 0, 0, 0}),
            canonical_to_montgomery({0xaaaaaaaaaaaaaaaaULL, 0x5555555555555555ULL,
                0x3333333333333333ULL, 0x1111111111111111ULL, 0, 0}),
            canonical_to_montgomery({UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX,
                UINT64_MAX, 0x1000000000000000ULL}),
        };
    }

    template<std::size_t batch>
    void expect_batch_products(
        const matrixnok::FpRingCandidate& ring,
        const std::array<blst::Fp, 8>& values)
    {
        std::array<matrixnok::FpRingCandidate::Element, batch> left;
        std::array<matrixnok::FpRingCandidate::Element, batch> right;
        std::array<blst::Fp, batch> expected;
        for (std::size_t index = 0; index < batch; ++index) {
            const blst::Fp& left_value = values[index % values.size()];
            const blst::Fp& right_value = values[(3 * index + 5) % values.size()];
            left[index] = ring.from_canonical_50bit_digits(canonical_50bit_digits(left_value));
            right[index] = ring.from_canonical_50bit_digits(canonical_50bit_digits(right_value));
            expected[index] = scalar_product(left_value, right_value);
        }

        const auto products = ring.multiply_batch(left, right);
        const auto production_products = ring.batch_modmul(left, right);
        for (std::size_t index = 0; index < batch; ++index) {
            expect_standard_residue_bounds(ring, products[index]);
            expect_standard_residue_bounds(ring, production_products[index]);
            expect_same_scalar_fp(matrixnok_to_scalar(ring, products[index]), expected[index]);
            expect_same_scalar_fp(
                matrixnok_to_scalar(ring, production_products[index]), expected[index]);
            EXPECT_TRUE(products[index].equal(production_products[index]));
        }
    }

#ifndef NDEBUG
    TEST(AvxMatrixNoKPrimitives, FixedRawContractsRejectInvalidInputs)
    {
        const matrixnok::FpRingCandidate ring;

        std::array<uint64_t, matrixnok::FpRingCandidate::limbs> oversized_digits = {};
        oversized_digits[0] = uint64_t{1} << matrixnok::FpRingCandidate::element_bits;
        EXPECT_DEATH((void)ring.from_canonical_50bit_digits(oversized_digits),
            "BLST_ASSERT");

        blst::Fp noncanonical = blst::get_BLST_P();
        EXPECT_DEATH((void)ring.from_canonical_limbs(noncanonical), "BLST_ASSERT");

        AVXVector<8> residues(0);
        AVXVector<8> low(0);
        EXPECT_DEATH(
            matrixnok::r1_change_base_raw(residues, residues, low),
            "BLST_ASSERT");

        std::array<uint64_t, 8> invalid_r1 = {};
        invalid_r1[0] = 2 * matrixnok::generated::moduli1[0];
        AVXVector<8> invalid_residues(invalid_r1);
        AVXVector<8> high(0);
        EXPECT_DEATH(
            matrixnok::r1_change_base_raw(invalid_residues, high, low),
            "BLST_ASSERT");
    }
#endif

    TEST(AvxMatrixNoKPrimitives, CanonicalExportCollapsesRedundantZeroLifts)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::FpRingCandidate::StandardElement modulus{
            matrixnok::FpRingCandidate::StandardResidue(
                AVXVector<8>(matrixnok_constants::modulus_in_rns_m1)),
            matrixnok::FpRingCandidate::StandardResidue(
                AVXVector<8>(matrixnok_constants::modulus_in_rns_m2))};
        const auto normalized = ring.batch_reduce_expand(
            ring.prep_left(modulus) * ring.one())[0];

        EXPECT_TRUE(blst::Fp::equal(
            matrixnok_to_scalar(ring, ring.zero()), blst::Fp::zero()));
        EXPECT_TRUE(blst::Fp::equal(
            matrixnok_to_scalar(ring, modulus), blst::Fp::zero()));
        EXPECT_TRUE(blst::Fp::equal(
            matrixnok_to_scalar(ring, normalized), blst::Fp::zero()));
    }

    TEST(AvxMatrixNoKPrimitives, ProductionPointEqualityUsesFieldSemantics)
    {
        const auto& ring = avx_curve_detail::curve_ring();
        const auto& fp2ring = avx_curve_detail::curve_fp2_ring();
        const matrixnok::FpRingCandidate::StandardElement modulus{
            matrixnok::FpRingCandidate::StandardResidue(
                AVXVector<8>(matrixnok_constants::modulus_in_rns_m1)),
            matrixnok::FpRingCandidate::StandardResidue(
                AVXVector<8>(matrixnok_constants::modulus_in_rns_m2))};

        const G1Point g1 = G1Point::generator();
        blst::Fp canonical_two = {};
        canonical_two.l[0] = 2;
        const auto scale = ring.from_canonical_limbs(canonical_two);
        const G1Point scaled_g1{
            ring.modmul(g1.X, scale),
            ring.modmul(g1.Y, scale),
            ring.modmul(g1.Z, scale)};
        EXPECT_TRUE(G1Point::equal(g1, scaled_g1));

        const G1AffinePoint g1_affine = G1AffinePoint::generator();
        const G1AffinePoint lifted_g1_affine(
            NormalizePointCoord(g1_affine.X + modulus, ring),
            NormalizePointCoord(g1_affine.Y + modulus, ring));
        EXPECT_FALSE(g1_affine.X.equal(lifted_g1_affine.X) &
            g1_affine.Y.equal(lifted_g1_affine.Y));
        EXPECT_TRUE(G1AffinePoint::equal(g1_affine, lifted_g1_affine));

        const G2Point g2 = G2Point::generator();
        const auto fp2_scale = avx_curve_detail::CurveFp2Ring::StandardElement(
            scale, ring.zero());
        const G2Point scaled_g2{
            fp2ring.mul(g2.X, fp2_scale),
            fp2ring.mul(g2.Y, fp2_scale),
            fp2ring.mul(g2.Z, fp2_scale)};
        EXPECT_TRUE(G2Point::equal(g2, scaled_g2));

        const auto lift_fp2 = [&modulus, &ring](const auto& value) {
            return avx_curve_detail::CurveFp2Ring::StandardElement(
                NormalizePointCoord(value.x() + modulus, ring),
                NormalizePointCoord(value.y() + modulus, ring));
        };
        const G2AffinePoint g2_affine = G2AffinePoint::generator();
        const G2AffinePoint lifted_g2_affine(
            lift_fp2(g2_affine.X),
            lift_fp2(g2_affine.Y));
        EXPECT_FALSE(g2_affine.X.equal(lifted_g2_affine.X) &
            g2_affine.Y.equal(lifted_g2_affine.Y));
        EXPECT_TRUE(G2AffinePoint::equal(g2_affine, lifted_g2_affine));
    }

    TEST(AvxMatrixNoKPrimitives, BaseFieldMultiplyMatchesScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const auto values = matrixnok_scalar_values();

        expect_same_scalar_fp(matrixnok_to_scalar(ring, ring.zero()), values[0]);
        expect_same_scalar_fp(matrixnok_to_scalar(ring, ring.one()), values[1]);

        for (const blst::Fp& value : values) {
            const auto encoded = ring.from_canonical_50bit_digits(canonical_50bit_digits(value));
            expect_standard_residue_bounds(ring, encoded);
            expect_same_scalar_fp(matrixnok_to_scalar(ring, encoded), value);
        }

        for (const blst::Fp& left : values) {
            for (const blst::Fp& right : values) {
                const auto candidate_left = ring.from_canonical_50bit_digits(canonical_50bit_digits(left));
                const auto candidate_right = ring.from_canonical_50bit_digits(canonical_50bit_digits(right));
                matrixnok::FpRingCandidate::Element candidate_product;
                ring.multiply(candidate_product, candidate_left, candidate_right);
                expect_standard_residue_bounds(ring, candidate_product);
                expect_same_scalar_fp(
                    matrixnok_to_scalar(ring, candidate_product),
                    scalar_product(left, right));
                expect_same_scalar_fp(
                    matrixnok_to_scalar(
                        ring, ring.modmul(candidate_left, candidate_right)),
                    scalar_product(left, right));
            }
        }

        expect_batch_products<1>(ring, values);
        expect_batch_products<2>(ring, values);
        expect_batch_products<4>(ring, values);
        expect_batch_products<8>(ring, values);
        expect_batch_products<12>(ring, values);

        uint64_t state = 0xd1b54a32d192ed03ULL;
        auto next_random = [&state]() {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            return state;
        };
        for (int trial = 0; trial < 64; ++trial) {
            std::array<uint64_t, blst::Fp::num_limb()> left_limbs = {};
            std::array<uint64_t, blst::Fp::num_limb()> right_limbs = {};
            for (std::size_t limb = 0; limb + 1 < blst::Fp::num_limb(); ++limb) {
                left_limbs[limb] = next_random();
                right_limbs[limb] = next_random();
            }
            // A top limb below p's top limb makes these dense values canonical
            // without sharing a reduction implementation with the oracle.
            left_limbs.back() = next_random() & 0x0fffffffffffffffULL;
            right_limbs.back() = next_random() & 0x0fffffffffffffffULL;
            const blst::Fp left = canonical_to_montgomery(left_limbs);
            const blst::Fp right = canonical_to_montgomery(right_limbs);
            const auto candidate_left = ring.from_canonical_50bit_digits(canonical_50bit_digits(left));
            const auto candidate_right = ring.from_canonical_50bit_digits(canonical_50bit_digits(right));
            matrixnok::FpRingCandidate::Element candidate_product;
            ring.multiply(candidate_product, candidate_left, candidate_right);
            expect_standard_residue_bounds(ring, candidate_product);
            expect_same_scalar_fp(
                matrixnok_to_scalar(ring, candidate_product),
                scalar_product(left, right));
        }

        auto candidate_accumulator = ring.from_canonical_50bit_digits(canonical_50bit_digits(values[4]));
        blst::Fp scalar_accumulator = values[4];
        for (int iteration = 0; iteration < 64; ++iteration) {
            ring.multiply(
                candidate_accumulator,
                candidate_accumulator,
                ring.from_canonical_50bit_digits(canonical_50bit_digits(values[6])));
            scalar_accumulator = scalar_product(scalar_accumulator, values[6]);
            expect_standard_residue_bounds(ring, candidate_accumulator);
            expect_same_scalar_fp(matrixnok_to_scalar(ring, candidate_accumulator), scalar_accumulator);
        }

        auto right_alias = ring.from_canonical_50bit_digits(canonical_50bit_digits(values[5]));
        const auto left_for_alias = ring.from_canonical_50bit_digits(canonical_50bit_digits(values[7]));
        ring.multiply(right_alias, left_for_alias, right_alias);
        expect_standard_residue_bounds(ring, right_alias);
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, right_alias),
            scalar_product(values[7], values[5]));
    }

    TEST(AvxMatrixNoKPrimitives, CanonicalImportKeepsBothRnsBasesConsistent)
    {
        const matrixnok::FpRingCandidate ring;
        const auto canonical_limbs = std::to_array<uint64_t>({
            11721118646859857034ULL,
            14515276896522647911ULL,
            16519252582114724462ULL,
            14272158255930264238ULL,
            5025824079061901589ULL,
            593286955712011240ULL,
        });
        const auto digits = std::to_array<uint64_t>({
            500616628141194ULL,
            1029601075259562ULL,
            686812189624822ULL,
            769198589660543ULL,
            581522751939018ULL,
            26173626664305ULL,
            89497429892085ULL,
            552541535ULL,
        });
        const blst::Fp expected = canonical_to_montgomery(canonical_limbs);
        const auto imported = ring.from_canonical_50bit_digits(digits);

        // Conversion back to the scalar field must preserve every canonical
        // input. The inherited ceiling estimate selects a negative lift for
        // this witness, outside ConvertFrom's nonnegative contract.
        expect_same_scalar_fp(matrixnok_to_scalar(ring, imported), expected);

        // Multiplication consumes both bases. Multiplying by the field identity
        // must therefore preserve the imported value and exposes any mismatch.
        const auto multiplied_by_one = ring.multiply(imported, ring.one());
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, multiplied_by_one), expected);
    }

    TEST(AvxMatrixNoKPrimitives, SignedDelayedProductsMatchScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const auto values = matrixnok_scalar_values();
        const auto left = ring.from_scalar(values[4]);
        const auto right = ring.from_scalar(values[6]);

        blst::Fp scalar_difference = {};
        blst::Fp::sub(scalar_difference, values[4], values[6]);
        const auto normalized_difference = ring.batch_reduce_expand(
            ring.prep_left(left - right) * ring.one())[0];
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, normalized_difference), scalar_difference);

        blst::Fp scalar_sum = {};
        blst::Fp::add(scalar_sum, values[4], values[6]);
        const auto signed_product = ring.batch_reduce_expand(
            ring.prep_left(left - right) * ring.prep(left + right))[0];
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, signed_product),
            scalar_product(scalar_difference, scalar_sum));

        const auto [ll, rr, lr, rl] = ring.batch_reduce_expand(
            ring.prep_left(left) * left,
            ring.prep_left(right) * right,
            ring.prep_left(left) * right,
            ring.prep_left(right) * left);
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, ll), scalar_product(values[4], values[4]));
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, rr), scalar_product(values[6], values[6]));
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, lr), scalar_product(values[4], values[6]));
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, rl), scalar_product(values[6], values[4]));
    }

    matrixnok::G1ProjectiveCandidate matrixnok_g1_from_scalar(
        const matrixnok::FpRingCandidate& ring,
        const blst::scalar::G1Point& point)
    {
        if (blst::scalar::G1Point::is_inf(point)) {
            return matrixnok::G1CurveCandidate{}.zero(ring);
        }
        blst::scalar::G1AffinePoint affine = {};
        blst::scalar::G1Point::to_affine(affine, point);
        return matrixnok::G1ProjectiveCandidate(
            ring.from_canonical_50bit_digits(canonical_50bit_digits(affine.x)),
            ring.from_canonical_50bit_digits(canonical_50bit_digits(affine.y)),
            ring.one());
    }

    blst::scalar::G1Point matrixnok_g1_to_scalar(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::G1ProjectiveCandidate& point)
    {
        // The AVX formulas use homogeneous coordinates. Convert to the scalar
        // Jacobian representation without an inversion:
        // (Xh*Zh, Yh*Zh^2, Zh).
        const auto z_squared = ring.modmul(point.z, point.z);
        blst::scalar::G1Point result = {};
        result.x = matrixnok_to_scalar(ring, ring.modmul(point.x, point.z));
        result.y = matrixnok_to_scalar(ring, ring.modmul(point.y, z_squared));
        result.z = matrixnok_to_scalar(ring, point.z);
        return result;
    }

    void expect_same_scalar_g1(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::G1ProjectiveCandidate& actual,
        const blst::scalar::G1Point& expected)
    {
        const auto decoded = matrixnok_g1_to_scalar(ring, actual);
        EXPECT_TRUE(blst::scalar::G1Point::on_curve(decoded));
        EXPECT_TRUE(blst::scalar::G1Point::equal(decoded, expected));
    }

    TEST(AvxMatrixNoKPrimitives, G1CompleteOperationsMatchScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G1CurveCandidate curve;
        const auto generator = blst::scalar::G1Point::generator();
        const auto candidate_generator = matrixnok_g1_from_scalar(ring, generator);
        const auto candidate_infinity = curve.zero(ring);
        ASSERT_TRUE(blst::scalar::G1Point::on_curve(generator));
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, candidate_generator.x), generator.x);
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, candidate_generator.y), generator.y);
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, candidate_generator.z), blst::Fp::one());
        const auto decoded_generator = matrixnok_g1_to_scalar(ring, candidate_generator);
        expect_same_scalar_fp(decoded_generator.x, generator.x);
        expect_same_scalar_fp(decoded_generator.y, generator.y);
        expect_same_scalar_fp(decoded_generator.z, generator.z);
        expect_same_scalar_g1(ring, candidate_generator, generator);
        expect_same_scalar_g1(
            ring, candidate_infinity, blst::scalar::G1Point::infinity());

        blst::scalar::G1Point twice = {};
        blst::scalar::G1Point::doubling(twice, generator);
        const auto candidate_doubled = curve.double_point(candidate_generator, ring);
        const auto fp_mul = [](const blst::Fp& a, const blst::Fp& b) {
            blst::Fp out = {};
            blst::Fp::mul(out, a, b);
            return out;
        };
        const auto fp_add = [](const blst::Fp& a, const blst::Fp& b) {
            blst::Fp out = {};
            blst::Fp::add(out, a, b);
            return out;
        };
        const auto fp_sub = [](const blst::Fp& a, const blst::Fp& b) {
            blst::Fp out = {};
            blst::Fp::sub(out, a, b);
            return out;
        };
        const auto fp_scale = [&](const blst::Fp& a, unsigned multiplier) {
            blst::Fp out = {};
            for (unsigned index = 0; index < multiplier; ++index) {
                out = fp_add(out, a);
            }
            return out;
        };
        const auto yy = fp_mul(generator.y, generator.y);
        const auto xy = fp_mul(generator.x, generator.y);
        const auto yz = generator.y;
        const auto eight_yy = fp_scale(yy, 8);
        const auto b3_zz = fp_scale(blst::Fp::one(), 12);
        const auto yy_plus_b3_zz = fp_add(yy, b3_zz);
        const auto yy_minus_b9_zz = fp_sub(
            yy, fp_scale(b3_zz, 3));
        const auto expected_x = fp_mul(
            yy_minus_b9_zz, fp_scale(xy, 2));
        const auto expected_z = fp_mul(eight_yy, yz);
        const auto expected_y = fp_add(
            fp_mul(yy_minus_b9_zz, yy_plus_b3_zz),
            fp_mul(eight_yy, b3_zz));
        {
            SCOPED_TRACE("MatrixNoK doubled homogeneous x");
            expect_same_scalar_fp(
                matrixnok_to_scalar(ring, candidate_doubled.x), expected_x);
        }
        {
            SCOPED_TRACE("MatrixNoK doubled homogeneous y");
            expect_same_scalar_fp(
                matrixnok_to_scalar(ring, candidate_doubled.y), expected_y);
        }
        {
            SCOPED_TRACE("MatrixNoK doubled homogeneous z");
            expect_same_scalar_fp(
                matrixnok_to_scalar(ring, candidate_doubled.z), expected_z);
        }
        {
            SCOPED_TRACE("MatrixNoK doubled Jacobian x product");
            expect_same_scalar_fp(
                matrixnok_to_scalar(
                    ring, ring.modmul(candidate_doubled.x, candidate_doubled.z)),
                fp_mul(expected_x, expected_z));
        }
        {
            SCOPED_TRACE("MatrixNoK doubled Jacobian y product");
            const auto z2 = ring.modmul(candidate_doubled.z, candidate_doubled.z);
            expect_same_scalar_fp(
                matrixnok_to_scalar(ring, ring.modmul(candidate_doubled.y, z2)),
                fp_mul(expected_y, fp_mul(expected_z, expected_z)));
        }
        expect_same_scalar_g1(
            ring, candidate_doubled, twice);
        expect_same_scalar_g1(
            ring, curve.add_point(candidate_generator, candidate_generator, ring), twice);

        blst::scalar::G1Point thrice = {};
        blst::scalar::G1Point::dadd(thrice, twice, generator);
        const auto candidate_twice = matrixnok_g1_from_scalar(ring, twice);
        expect_same_scalar_g1(
            ring, curve.add_point(candidate_twice, candidate_generator, ring), thrice);

        blst::scalar::G1Point negative = generator;
        blst::scalar::G1Point::cneg(negative, 1);
        const auto candidate_negative = matrixnok_g1_from_scalar(ring, negative);
        expect_same_scalar_g1(
            ring,
            curve.add_point(candidate_generator, candidate_negative, ring),
            blst::scalar::G1Point::infinity());
        expect_same_scalar_g1(
            ring, curve.add_point(candidate_infinity, candidate_generator, ring),
            generator);
        expect_same_scalar_g1(
            ring, curve.add_point(candidate_generator, candidate_infinity, ring),
            generator);

        uint64_t state = 0x94d049bb133111ebULL;
        blst::scalar::G1Point running = generator;
        auto candidate_running = candidate_generator;
        for (int trial = 0; trial < 32; ++trial) {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            const unsigned steps = static_cast<unsigned>(state & 7U) + 1U;
            blst::scalar::G1Point addend = generator;
            for (unsigned step = 1; step < steps; ++step) {
                blst::scalar::G1Point next = {};
                blst::scalar::G1Point::dadd(next, addend, generator);
                addend = next;
            }
            const auto candidate_addend = matrixnok_g1_from_scalar(ring, addend);
            blst::scalar::G1Point expected = {};
            blst::scalar::G1Point::dadd(expected, running, addend);
            candidate_running = curve.add_point(
                candidate_running, candidate_addend, ring);
            expect_same_scalar_g1(ring, candidate_running, expected);
            running = expected;
        }
    }

    TEST(AvxMatrixNoKPrimitives, G1GlvMultiplicationMatchesScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G1CurveCandidate curve;
        const auto& glv = avx_curve_detail::curve_g1_glv();
        const auto generator = blst::scalar::G1Point::generator();
        blst::scalar::G1Point doubled = {};
        blst::scalar::G1Point::doubling(doubled, generator);
        const std::array<blst::scalar::G1Point, 2> bases = {
            generator, doubled};
        std::vector<std::array<uint8_t, 32>> scalars = {{
            {0},
            {1},
            {2},
            {0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88,
             0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x01},
            {0xa5, 0x5a, 0xc3, 0x3c, 0x96, 0x69, 0xf0, 0x0f,
             0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0,
             0x0f, 0xed, 0xcb, 0xa9, 0x87, 0x65, 0x43, 0x21,
             0x10, 0x32, 0x54, 0x76, 0x08, 0x00, 0x00, 0x00},
        }};

        constexpr std::array<uint64_t, 4> scalar_order_minus_one = {
            0xffffffff00000000ULL,
            0x53bda402fffe5bfeULL,
            0x3339d80809a1d805ULL,
            0x73eda753299d7d48ULL,
        };
        std::array<uint8_t, 32> boundary = {};
        for (std::size_t limb = 0; limb < scalar_order_minus_one.size(); ++limb) {
            for (std::size_t byte = 0; byte < sizeof(uint64_t); ++byte) {
                boundary[limb * sizeof(uint64_t) + byte] =
                    static_cast<uint8_t>(
                        scalar_order_minus_one[limb] >> (byte * 8));
            }
        }
        scalars.push_back(boundary);

        uint64_t state = 0x243f6a8885a308d3ULL;
        for (int trial = 0; trial < 16; ++trial) {
            std::array<uint8_t, 32> scalar = {};
            for (uint8_t& byte : scalar) {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                byte = static_cast<uint8_t>(state);
            }
            scalar.back() &= 0x3f;
            scalars.push_back(scalar);
        }

        for (const auto& base : bases) {
            const auto candidate_base = matrixnok_g1_from_scalar(ring, base);
            for (const auto& scalar : scalars) {
                blst::scalar::G1Point expected = {};
                blst::scalar::G1Point::mult_vartime(
                    expected, base, scalar.data(), 256);
                expect_same_scalar_g1(
                    ring,
                    glv.multiply_vartime(candidate_base, scalar, ring),
                    expected);
                expect_same_scalar_g1(
                    ring,
                    glv.multiply(candidate_base, scalar, ring),
                    expected);
            }
        }
    }

    FP2Ring<matrixnok::FpRingCandidate>::StandardElement
    matrixnok_generic_fp2_from_scalar(
        const matrixnok::FpRingCandidate& ring,
        const blst::Fp2& value)
    {
        return {ring.from_scalar(value[0]), ring.from_scalar(value[1])};
    }

    TEST(AvxMatrixNoKPrimitives, SelectedFp2CompletionBoundaryMatchesScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const FP2Ring<matrixnok::FpRingCandidate> fp2ring(ring);
        // These canonical operands make the abstract signed real-coordinate
        // expression land on a deliberately weak R2 rounding-margin witness.
        // The selected physical path must apply its negative-range completion
        // before reduction, retain scalar semantics, and return one Fp2 value
        // rather than leaking the one-element batch used internally.
        const blst::Fp2 left = {{
            canonical_to_montgomery({1, 0, 0, 0, 0, 0}),
            canonical_to_montgomery({
                13402431016077863594ULL,
                2210141511517208575ULL,
                7435674573564081700ULL,
                7239337960414712511ULL,
                5412103778470702295ULL,
                1873798617647539866ULL,
            }),
        }};
        const blst::Fp2 right = {{
            canonical_to_montgomery({
                1618306854713067172ULL,
                6995949995587336674ULL,
                6898866514068117168ULL,
                3807786787093763224ULL,
                7056063295788903324ULL,
                203523531925434294ULL,
            }),
            canonical_to_montgomery({
                9396560402482380242ULL,
                6808606215268122910ULL,
                13634768106530980159ULL,
                11464363280807193251ULL,
                9602730034778851917ULL,
                8073234691677ULL,
            }),
        }};

        blst::Fp2 expected = {};
        blst::Fp2::mul(expected, left, right);
        const auto actual = fp2ring.mul(
            matrixnok_generic_fp2_from_scalar(ring, left),
            matrixnok_generic_fp2_from_scalar(ring, right));
        expect_standard_residue_bounds(ring, actual.x());
        expect_standard_residue_bounds(ring, actual.y());
        expect_same_scalar_fp(matrixnok_to_scalar(ring, actual.x()), expected[0]);
        expect_same_scalar_fp(matrixnok_to_scalar(ring, actual.y()), expected[1]);

        blst::Fp2 expected_square = {};
        blst::Fp2::sqr(expected_square, left);
        const auto actual_square = fp2ring.square(
            matrixnok_generic_fp2_from_scalar(ring, left));
        expect_standard_residue_bounds(ring, actual_square.x());
        expect_standard_residue_bounds(ring, actual_square.y());
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, actual_square.x()), expected_square[0]);
        expect_same_scalar_fp(
            matrixnok_to_scalar(ring, actual_square.y()), expected_square[1]);
    }

    matrixnok::G2ProjectiveCandidate matrixnok_g2_from_scalar(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::G2Fp2RingCandidate& fp2ring,
        const blst::scalar::G2Point& point)
    {
        const matrixnok::G2CurveCandidate curve;
        if (blst::scalar::G2Point::is_inf(point))
            return curve.zero(fp2ring);
        blst::scalar::G2AffinePoint affine = {};
        blst::scalar::G2Point::to_affine(affine, point);
        return matrixnok::G2ProjectiveCandidate(
            matrixnok_generic_fp2_from_scalar(ring, affine.x),
            matrixnok_generic_fp2_from_scalar(ring, affine.y),
            fp2ring.one());
    }

    blst::scalar::G2Point matrixnok_g2_to_scalar(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::G2ProjectiveCandidate& point)
    {
        const auto decode = [&](const auto& value) {
            return blst::Fp2{{
                matrixnok_to_scalar(ring, value.x()),
                matrixnok_to_scalar(ring, value.y())}};
        };
        const auto multiply = [](const blst::Fp2& left, const blst::Fp2& right) {
            blst::Fp2 result = {};
            blst::Fp2::mul(result, left, right);
            return result;
        };
        const blst::Fp2 x = decode(point.x);
        const blst::Fp2 y = decode(point.y);
        const blst::Fp2 z = decode(point.z);
        const blst::Fp2 z_squared = multiply(z, z);
        blst::scalar::G2Point result = {};
        result.x = multiply(x, z);
        result.y = multiply(y, z_squared);
        result.z = z;
        return result;
    }

    void expect_same_scalar_g2(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::G2ProjectiveCandidate& actual,
        const blst::scalar::G2Point& expected)
    {
        const auto decoded = matrixnok_g2_to_scalar(ring, actual);
        EXPECT_TRUE(blst::scalar::G2Point::on_curve(decoded));
        EXPECT_TRUE(blst::scalar::G2Point::equal(decoded, expected));
    }

    TEST(AvxMatrixNoKPrimitives, G2CompleteOperationsMatchScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        const matrixnok::G2CurveCandidate curve;
        const auto generator = blst::scalar::G2Point::generator();
        const auto candidate_generator =
            matrixnok_g2_from_scalar(ring, fp2ring, generator);
        const auto candidate_infinity = curve.zero(fp2ring);
        ASSERT_TRUE(blst::scalar::G2Point::on_curve(generator));
        {
            SCOPED_TRACE("generator conversion");
            expect_same_scalar_g2(ring, candidate_generator, generator);
        }
        {
            SCOPED_TRACE("infinity conversion");
            expect_same_scalar_g2(
                ring, candidate_infinity, blst::scalar::G2Point::infinity());
        }

        blst::scalar::G2Point twice = {};
        blst::scalar::G2Point::doubling(twice, generator);
        const auto candidate_doubled =
            curve.double_point(candidate_generator, fp2ring);
        {
            SCOPED_TRACE("double");
            expect_same_scalar_g2(ring, candidate_doubled, twice);
        }
        {
            SCOPED_TRACE("self add");
            expect_same_scalar_g2(
                ring,
                curve.add_point(candidate_generator, candidate_generator, fp2ring),
                twice);
        }

        blst::scalar::G2Point thrice = {};
        blst::scalar::G2Point::dadd(thrice, twice, generator);
        {
            SCOPED_TRACE("distinct add");
            expect_same_scalar_g2(
                ring,
                curve.add_point(candidate_doubled, candidate_generator, fp2ring),
                thrice);
        }

        blst::scalar::G2Point negative = generator;
        blst::scalar::G2Point::cneg(negative, 1);
        const auto candidate_negative =
            matrixnok_g2_from_scalar(ring, fp2ring, negative);
        {
            SCOPED_TRACE("opposites");
            expect_same_scalar_g2(
                ring,
                curve.add_point(candidate_generator, candidate_negative, fp2ring),
                blst::scalar::G2Point::infinity());
        }
        expect_same_scalar_g2(
            ring,
            curve.add_point(candidate_infinity, candidate_generator, fp2ring),
            generator);
        expect_same_scalar_g2(
            ring,
            curve.add_point(candidate_generator, candidate_infinity, fp2ring),
            generator);

        auto candidate_running = candidate_generator;
        blst::scalar::G2Point running = generator;
        for (unsigned step = 1; step <= 16; ++step) {
            blst::scalar::G2Point expected = {};
            blst::scalar::G2Point::dadd(expected, running, generator);
            candidate_running = curve.add_point(
                candidate_running, candidate_generator, fp2ring);
            SCOPED_TRACE(step);
            expect_same_scalar_g2(ring, candidate_running, expected);
            running = expected;
        }
    }

    TEST(AvxMatrixNoKPrimitives, G2GlsMultiplicationMatchesScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        const matrixnok::G2CurveCandidate curve;
        const auto& gls = avx_curve_detail::curve_g2_gls();
        const auto generator = blst::scalar::G2Point::generator();
        blst::scalar::G2Point doubled = {};
        blst::scalar::G2Point::doubling(doubled, generator);
        const std::array<blst::scalar::G2Point, 2> bases = {
            generator, doubled};
        std::vector<std::array<uint8_t, 32>> scalars = {{
            {0},
            {1},
            {2},
            {0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88,
             0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11, 0x01},
            {0xa5, 0x5a, 0xc3, 0x3c, 0x96, 0x69, 0xf0, 0x0f,
             0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0,
             0x0f, 0xed, 0xcb, 0xa9, 0x87, 0x65, 0x43, 0x21,
             0x10, 0x32, 0x54, 0x76, 0x08, 0x00, 0x00, 0x00},
        }};
        constexpr std::array<uint64_t, 4> scalar_order_minus_one = {
            0xffffffff00000000ULL,
            0x53bda402fffe5bfeULL,
            0x3339d80809a1d805ULL,
            0x73eda753299d7d48ULL,
        };
        std::array<uint8_t, 32> boundary = {};
        for (std::size_t limb = 0; limb < scalar_order_minus_one.size(); ++limb) {
            for (std::size_t byte = 0; byte < sizeof(uint64_t); ++byte) {
                boundary[limb * sizeof(uint64_t) + byte] =
                    static_cast<uint8_t>(
                        scalar_order_minus_one[limb] >> (byte * 8));
            }
        }
        scalars.push_back(boundary);

        uint64_t state = 0x13198a2e03707344ULL;
        for (int trial = 0; trial < 16; ++trial) {
            std::array<uint8_t, 32> scalar = {};
            for (uint8_t& byte : scalar) {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                byte = static_cast<uint8_t>(state);
            }
            scalar.back() &= 0x3f;
            scalars.push_back(scalar);
        }

        for (const auto& base : bases) {
            const auto candidate_base =
                matrixnok_g2_from_scalar(ring, fp2ring, base);
            for (const auto& scalar : scalars) {
                blst::scalar::G2Point expected = {};
                blst::scalar::G2Point::mult_vartime(
                    expected, base, scalar.data(), 256);
                expect_same_scalar_g2(
                    ring,
                    gls.multiply_vartime(
                        candidate_base, scalar, fp2ring, curve, ring),
                    expected);
                expect_same_scalar_g2(
                    ring,
                    gls.multiply(candidate_base, scalar, fp2ring, curve, ring),
                    expected);
            }
        }
    }

    struct MatrixNoKG1FftCtx
    {
        inline static const matrixnok::FpRingCandidate ring = {};
        inline static const matrixnok::G1CurveCandidate curve = {};
        inline static const matrixnok::G1GLVCandidate& glv =
            avx_curve_detail::kCurveG1Glv;

        static matrixnok::G1ProjectiveCandidate plus(
            const matrixnok::G1ProjectiveCandidate& left,
            const matrixnok::G1ProjectiveCandidate& right)
        {
            return curve.add_point(left, right, ring);
        }

        static matrixnok::G1ProjectiveCandidate minus(
            const matrixnok::G1ProjectiveCandidate& left,
            const matrixnok::G1ProjectiveCandidate& right)
        {
            return plus(left, curve.negate(right, ring));
        }

        static matrixnok::G1ProjectiveCandidate mult(
            const matrixnok::G1ProjectiveCandidate& point,
            const blst::Fr& scalar)
        {
            blst::Fr canonical = {};
            blst::Fr::from_mont(canonical, scalar);
            std::array<uint8_t, 32> bytes = {};
            blst::le_bytes_from_limbs(bytes.data(), canonical.l, bytes.size());
            return glv.multiply_vartime(point, bytes, ring);
        }

        static blst::Fr inverse(const blst::Fr& value)
        {
            blst::Fr result = {};
            blst::Fr::reciprocal(result, value);
            return result;
        }

        template<class>
        static blst::Fr from_u64(const uint64_t value)
        {
            return blst::Fr::from_uint64(value);
        }
    };

    TEST(AvxMatrixNoKPrimitives, G1FftRoundTripMatchesScalar)
    {
        constexpr std::size_t size = 16;
        const matrixnok::FpRingCandidate ring;
        std::vector<blst::Fr> roots(size);
        std::vector<blst::Fr> inverse_roots(size);
        blst::Fr root = blst::primRootOfUnity<blst::Fr>(size);
        blst::Fr inverse_root = {};
        blst::Fr::reciprocal(inverse_root, root);
        roots[0] = blst::Fr::one();
        inverse_roots[0] = roots[0];
        for (std::size_t index = 1; index < size; ++index) {
            blst::Fr::mul(roots[index], roots[index - 1], root);
            blst::Fr::mul(
                inverse_roots[index], inverse_roots[index - 1], inverse_root);
        }

        std::vector<blst::scalar::G1Point> expected(size);
        std::vector<matrixnok::G1ProjectiveCandidate> points(size);
        const auto generator = blst::scalar::G1Point::generator();
        expected[0] = blst::scalar::G1Point::infinity();
        points[0] = matrixnok_g1_from_scalar(ring, expected[0]);
        for (std::size_t index = 1; index < size; ++index) {
            blst::scalar::G1Point::dadd(expected[index], expected[index - 1], generator);
            points[index] = matrixnok_g1_from_scalar(ring, expected[index]);
        }

        fft<matrixnok::G1ProjectiveCandidate, blst::Fr, MatrixNoKG1FftCtx>(
            points, roots, false);
        fft<matrixnok::G1ProjectiveCandidate, blst::Fr, MatrixNoKG1FftCtx>(
            points, inverse_roots, true);

        for (std::size_t index = 0; index < size; ++index) {
            SCOPED_TRACE(index);
            expect_same_scalar_g1(ring, points[index], expected[index]);
        }
    }

    blst::Fp2 matrixnok_fp2_to_scalar(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::G2Fp2RingCandidate::StandardElement& value)
    {
        return {{matrixnok_to_scalar(ring, value.x()),
            matrixnok_to_scalar(ring, value.y())}};
    }

    void expect_same_scalar_fp2(const blst::Fp2& actual, const blst::Fp2& expected)
    {
        expect_same_scalar_fp(actual[0], expected[0]);
        expect_same_scalar_fp(actual[1], expected[1]);
    }

    void expect_standard_fp2_bounds(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::G2Fp2RingCandidate::StandardElement& value)
    {
        expect_standard_residue_bounds(ring, value.x());
        expect_standard_residue_bounds(ring, value.y());
    }

    TEST(AvxMatrixNoKPrimitives, SelectedFp2MultiplyAndSquareMatchScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        testing::Test::RecordProperty("fp_ring_candidate_bytes", sizeof(ring));
        testing::Test::RecordProperty("fp2_ring_candidate_bytes", sizeof(fp2ring));

        const auto base_values = matrixnok_scalar_values();
        std::array<blst::Fp2, 8> values = {};
        for (std::size_t index = 0; index < values.size(); ++index) {
            values[index][0] = base_values[index];
            values[index][1] = base_values[(3 * index + 1) % base_values.size()];
        }

        for (const blst::Fp2& value : values) {
            const auto encoded = matrixnok_generic_fp2_from_scalar(ring, value);
            const auto candidate_square = fp2ring.square(encoded);
            blst::Fp2 expected_square = {};
            blst::Fp2::sqr(expected_square, value);
            expect_standard_fp2_bounds(ring, candidate_square);
            expect_same_scalar_fp2(
                matrixnok_fp2_to_scalar(ring, candidate_square), expected_square);
        }

        for (const blst::Fp2& left : values) {
            for (const blst::Fp2& right : values) {
                const auto candidate_left = matrixnok_generic_fp2_from_scalar(ring, left);
                const auto candidate_right = matrixnok_generic_fp2_from_scalar(ring, right);
                const auto candidate_product = fp2ring.mul(candidate_left, candidate_right);
                blst::Fp2 expected_product = {};
                blst::Fp2::mul(expected_product, left, right);
                expect_standard_fp2_bounds(ring, candidate_product);
                expect_same_scalar_fp2(
                    matrixnok_fp2_to_scalar(ring, candidate_product), expected_product);
            }
        }

        auto right_alias = matrixnok_generic_fp2_from_scalar(ring, values[5]);
        const auto left_for_alias = matrixnok_generic_fp2_from_scalar(ring, values[7]);
        right_alias = fp2ring.mul(left_for_alias, right_alias);
        blst::Fp2 expected_alias = {};
        blst::Fp2::mul(expected_alias, values[7], values[5]);
        expect_standard_fp2_bounds(ring, right_alias);
        expect_same_scalar_fp2(matrixnok_fp2_to_scalar(ring, right_alias), expected_alias);

        auto left_alias = matrixnok_generic_fp2_from_scalar(ring, values[4]);
        const auto right_for_alias = matrixnok_generic_fp2_from_scalar(ring, values[6]);
        left_alias = fp2ring.mul(left_alias, right_for_alias);
        blst::Fp2 expected_left_alias = {};
        blst::Fp2::mul(expected_left_alias, values[4], values[6]);
        expect_standard_fp2_bounds(ring, left_alias);
        expect_same_scalar_fp2(
            matrixnok_fp2_to_scalar(ring, left_alias), expected_left_alias);

        auto candidate_product_chain = matrixnok_generic_fp2_from_scalar(ring, values[4]);
        blst::Fp2 scalar_product_chain = values[4];
        const auto candidate_factor = matrixnok_generic_fp2_from_scalar(ring, values[6]);
        for (int iteration = 0; iteration < 64; ++iteration) {
            candidate_product_chain = fp2ring.mul(candidate_product_chain, candidate_factor);
            blst::Fp2 next = {};
            blst::Fp2::mul(next, scalar_product_chain, values[6]);
            scalar_product_chain = next;
            expect_standard_fp2_bounds(ring, candidate_product_chain);
            expect_same_scalar_fp2(
                matrixnok_fp2_to_scalar(ring, candidate_product_chain), scalar_product_chain);
        }

        auto candidate_square_chain = matrixnok_generic_fp2_from_scalar(ring, values[5]);
        blst::Fp2 scalar_square_chain = values[5];
        for (int iteration = 0; iteration < 16; ++iteration) {
            candidate_square_chain = fp2ring.square(candidate_square_chain);
            blst::Fp2 next = {};
            blst::Fp2::sqr(next, scalar_square_chain);
            scalar_square_chain = next;
            expect_standard_fp2_bounds(ring, candidate_square_chain);
            expect_same_scalar_fp2(
                matrixnok_fp2_to_scalar(ring, candidate_square_chain), scalar_square_chain);
        }
    }

    MatrixNoKFp12::StandardElement matrixnok_fp12_from_scalar(
        const matrixnok::FpRingCandidate& ring, const blst::Fp12& value)
    {
        MatrixNoKFp12::StandardElement out = {};
        for (std::size_t coefficient = 0; coefficient < out.size(); ++coefficient) {
            out[coefficient] = ring.from_scalar(
                value[coefficient / 6][(coefficient % 6) / 2][coefficient % 2]);
        }
        return out;
    }

    blst::Fp12 matrixnok_fp12_to_scalar(
        const matrixnok::FpRingCandidate& ring,
        const MatrixNoKFp12::StandardElement& value)
    {
        blst::Fp12 out = {};
        for (std::size_t coefficient = 0; coefficient < value.size(); ++coefficient) {
            out[coefficient / 6][(coefficient % 6) / 2][coefficient % 2] =
                matrixnok_to_scalar(ring, value[coefficient]);
        }
        return out;
    }

    void expect_same_scalar_fp12(
        const blst::Fp12& actual, const blst::Fp12& expected)
    {
        for (std::size_t coefficient = 0; coefficient < 12; ++coefficient) {
            SCOPED_TRACE(coefficient);
            expect_same_scalar_fp(
                actual[coefficient / 6][(coefficient % 6) / 2][coefficient % 2],
                expected[coefficient / 6][(coefficient % 6) / 2][coefficient % 2]);
        }
    }

    TEST(AvxMatrixNoKPrimitives, GtGeneratorLiteralMatchesScalar)
    {
        const matrixnok::FpRingCandidate ring;
        MatrixNoKFp12::StandardElement literal = {};
        for (std::size_t coefficient = 0; coefficient < literal.size(); ++coefficient) {
            literal[coefficient] = matrixnok::FpRingCandidate::StandardElement(
                matrixnok::FpRingCandidate::StandardResidue(
                    AVXVector<8>(matrixnok_constants::gt_generator[coefficient][0])),
                matrixnok::FpRingCandidate::StandardResidue(
                    AVXVector<8>(matrixnok_constants::gt_generator[coefficient][1])));
        }

        const auto converted = matrixnok_fp12_from_scalar(
            ring, blst::scalar::GT::generator().data);
        for (std::size_t coefficient = 0; coefficient < literal.size(); ++coefficient) {
            EXPECT_TRUE(literal[coefficient].equal(converted[coefficient]));
        }
        expect_same_scalar_fp12(
            matrixnok_fp12_to_scalar(ring, literal),
            blst::scalar::GT::generator().data);
    }

    TEST(AvxMatrixNoKPrimitives, Fp12MultiplyAndSquareMatchScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const MatrixNoKFp12 fp12ring{};
        const auto values = matrixnok_scalar_values();
        blst::Fp12 left = {};
        blst::Fp12 right = {};
        for (std::size_t coefficient = 0; coefficient < 12; ++coefficient) {
            left[coefficient / 6][(coefficient % 6) / 2][coefficient % 2] =
                values[coefficient % values.size()];
            right[coefficient / 6][(coefficient % 6) / 2][coefficient % 2] =
                values[(3 * coefficient + 1) % values.size()];
        }

        const auto candidate_left = matrixnok_fp12_from_scalar(ring, left);
        const auto candidate_right = matrixnok_fp12_from_scalar(ring, right);

        {
            SCOPED_TRACE("left times one");
            expect_same_scalar_fp12(
                matrixnok_fp12_to_scalar(
                    ring,
                    fp12ring.fp12_flat_mul(
                        candidate_left, fp12ring.one(ring), ring)),
                left);
        }

        blst::Fp12 expected_product = {};
        blst::Fp12::mul(expected_product, left, right);
        {
            SCOPED_TRACE("general product");
            expect_same_scalar_fp12(
                matrixnok_fp12_to_scalar(
                    ring, fp12ring.fp12_flat_mul(
                        candidate_left, candidate_right, ring)),
                expected_product);
        }

        blst::Fp12 expected_square = {};
        blst::Fp12::sqr(expected_square, left);
        {
            SCOPED_TRACE("square");
            expect_same_scalar_fp12(
                matrixnok_fp12_to_scalar(
                    ring, fp12ring.fp12_flat_sqr(candidate_left, ring)),
                expected_square);
        }

        {
            SCOPED_TRACE("one");
            expect_same_scalar_fp12(
                matrixnok_fp12_to_scalar(ring, fp12ring.one(ring)),
                blst::Fp12::one());
        }
    }

    TEST(AvxMatrixNoKPrimitives, Fp12InverseDecompositionMatchesScalar)
    {
        const matrixnok::FpRing ring;
        const MatrixNoKFp12 fp12ring{};
        const matrixnok::PairingInverter inverter{};
        const auto values = matrixnok_scalar_values();

        std::array<blst::Fp12, 4> inputs{};
        inputs[0] = blst::Fp12::one();
        for (std::size_t input = 1; input < inputs.size(); ++input) {
            for (std::size_t coefficient = 0; coefficient < 12; ++coefficient) {
                inputs[input][coefficient / 6]
                    [(coefficient % 6) / 2][coefficient % 2] =
                    values[(5 * coefficient + 3 * input) % values.size()];
            }
            // State the nonzero precondition independently of the fixture mix.
            inputs[input][0][0][0] = blst::Fp::one();
        }

        for (std::size_t input = 0; input < inputs.size(); ++input) {
            SCOPED_TRACE("input " + std::to_string(input));
            const auto candidate = matrixnok_fp12_from_scalar(ring, inputs[input]);
            const auto candidate_inverse =
                inverter.invert(candidate, fp12ring, ring);

            blst::Fp12 expected_inverse = {};
            blst::Fp12::inverse(expected_inverse, inputs[input]);
            expect_same_scalar_fp12(
                matrixnok_fp12_to_scalar(ring, candidate_inverse),
                expected_inverse);

            const auto candidate_product = fp12ring.fp12_flat_mul(
                candidate, candidate_inverse, ring);
            expect_same_scalar_fp12(
                matrixnok_fp12_to_scalar(ring, candidate_product),
                blst::Fp12::one());

            const auto candidate_round_trip =
                inverter.invert(candidate_inverse, fp12ring, ring);
            expect_same_scalar_fp12(
                matrixnok_fp12_to_scalar(ring, candidate_round_trip),
                inputs[input]);
        }
    }

    TEST(AvxMatrixNoKPrimitives, Fp12MapsAndGtExponentiationMatchScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const MatrixNoKFp12 fp12ring{};
        const auto values = matrixnok_scalar_values();
        blst::Fp12 general = {};
        for (std::size_t coefficient = 0; coefficient < 12; ++coefficient) {
            general[coefficient / 6][(coefficient % 6) / 2][coefficient % 2] =
                values[(5 * coefficient + 2) % values.size()];
        }
        const auto candidate_general = matrixnok_fp12_from_scalar(ring, general);

        for (std::size_t power = 1; power <= 3; ++power) {
            blst::Fp12 expected = {};
            blst::Fp12::frobenius_map(expected, general, power);
            MatrixNoKFp12::StandardElement actual = {};
            if (power == 1) {
                actual = fp12ring.fp12_flat_frobenius_map1(candidate_general, ring);
            } else if (power == 2) {
                actual = fp12ring.narrow_to_standard(
                    fp12ring.fp12_flat_frobenius_map2(candidate_general, ring), ring);
            } else {
                actual = fp12ring.fp12_flat_frobenius_map3(candidate_general, ring);
            }
            SCOPED_TRACE("Frobenius power " + std::to_string(power));
            expect_same_scalar_fp12(matrixnok_fp12_to_scalar(ring, actual), expected);
        }

        blst::Fp12 expected_conjugate = general;
        blst::Fp12::conjugate(expected_conjugate);
        expect_same_scalar_fp12(
            matrixnok_fp12_to_scalar(
                ring, fp12ring.fp12_flat_conjugate(candidate_general, ring)),
            expected_conjugate);

        blst::Fp12 gt = {};
        blst::Fp12::final_exp(gt, general);
        const auto candidate_gt = matrixnok_fp12_from_scalar(ring, gt);

        blst::Fp12 expected_square = {};
        blst::Fp12::cyclotomic_sqr(expected_square, gt);
        const auto candidate_square = fp12ring.narrow_to_standard(
            fp12ring.fp12_flat_cyclotomic_sqr(candidate_gt, ring), ring);
        expect_same_scalar_fp12(
            matrixnok_fp12_to_scalar(ring, candidate_square), expected_square);

        blst::Fp12 expected_z = {};
        blst::Fp12 expected_z_half = {};
        blst::Fp12::raise_to_z(expected_z, gt);
        blst::Fp12::raise_to_z_div_by_2(expected_z_half, gt);
        const auto candidate_z = fp12ring.narrow_to_standard(
            raise_to_z(candidate_gt, fp12ring, ring), ring);
        const auto candidate_z_half = fp12ring.narrow_to_standard(
            raise_to_z_div_by_2(candidate_gt, fp12ring, ring), ring);
        expect_same_scalar_fp12(
            matrixnok_fp12_to_scalar(ring, candidate_z), expected_z);
        expect_same_scalar_fp12(
            matrixnok_fp12_to_scalar(ring, candidate_z_half), expected_z_half);

        const std::array<uint64_t, 4> scalar = {
            0xfedcba9876543210ULL,
            0x0123456789abcdefULL,
            0x1020304050607080ULL,
            0x42eda753299d7d48ULL,
        };
        blst::Fp12 expected_ct = {};
        blst::Fp12 expected_vt = {};
        blst::Fp12::cyclotomic_exp_glv4_ct(expected_ct, gt, scalar.data());
        blst::Fp12::cyclotomic_exp_glv4(expected_vt, gt, scalar.data());
        const auto candidate_ct = avx_cyclotomic_exp_glv4_ct(
            candidate_gt, scalar.data(), fp12ring, ring);
        const auto candidate_vt = avx_cyclotomic_exp_glv4_unsigned_vartime(
            candidate_gt, scalar.data(), fp12ring, ring);
        expect_same_scalar_fp12(
            matrixnok_fp12_to_scalar(ring, candidate_ct), expected_ct);
        expect_same_scalar_fp12(
            matrixnok_fp12_to_scalar(ring, candidate_vt), expected_vt);
    }

    TEST(AvxMatrixNoKPrimitives, Fp12SparseMultiplyAndFinalExpMatchScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const MatrixNoKFp12 fp12ring{};
        const auto values = matrixnok_scalar_values();

        blst::Fp12 general = {};
        for (std::size_t coefficient = 0; coefficient < 12; ++coefficient) {
            general[coefficient / 6][(coefficient % 6) / 2][coefficient % 2] =
                values[(7 * coefficient + 3) % values.size()];
        }
        blst::Fp6 sparse = {};
        for (std::size_t coefficient = 0; coefficient < 6; ++coefficient) {
            sparse[coefficient / 2][coefficient % 2] =
                values[(5 * coefficient + 1) % values.size()];
        }

        const auto candidate_general = matrixnok_fp12_from_scalar(ring, general);
        const auto constant_term = matrixnok_generic_fp2_from_scalar(ring, sparse[0]);
        const auto x_term = matrixnok_generic_fp2_from_scalar(ring, sparse[1]);
        const auto y_term = matrixnok_generic_fp2_from_scalar(ring, sparse[2]);

        blst::Fp12 expected_sparse = {};
        blst::Fp12::mul_by_xy00z0(expected_sparse, general, sparse);
        MatrixNoKFp12::StandardElement candidate_sparse = {};
        fp12ring.flat_mulxy00z0(candidate_sparse, candidate_general,
            constant_term, x_term, y_term, ring);
        expect_same_scalar_fp12(
            matrixnok_fp12_to_scalar(ring, candidate_sparse), expected_sparse);

        auto candidate_alias = candidate_general;
        fp12ring.flat_mulxy00z0(candidate_alias, candidate_alias,
            constant_term, x_term, y_term, ring);
        expect_same_scalar_fp12(
            matrixnok_fp12_to_scalar(ring, candidate_alias), expected_sparse);

        for (const std::size_t shift : {0U, 1U, 3U}) {
            blst::Fp12 final_input = {};
            for (std::size_t coefficient = 0; coefficient < 12; ++coefficient) {
                final_input[coefficient / 6]
                    [(coefficient % 6) / 2][coefficient % 2] =
                    values[(7 * coefficient + 3 + shift) % values.size()];
            }
            blst::Fp12 scalar_inverse = {};
            blst::Fp12::inverse(scalar_inverse, final_input);
            const auto candidate_input =
                matrixnok_fp12_from_scalar(ring, final_input);
            const auto candidate_inverse =
                matrixnok_fp12_from_scalar(ring, scalar_inverse);
            const auto candidate_final = final_exp_after_inverse(
                candidate_input, candidate_inverse, fp12ring, ring);
            blst::Fp12 expected_final = {};
            blst::Fp12::final_exp(expected_final, final_input);
            SCOPED_TRACE("final exponentiation shift " + std::to_string(shift));
            expect_same_scalar_fp12(
                matrixnok_fp12_to_scalar(ring, candidate_final), expected_final);
        }
    }

    matrixnok::G1AffineCandidate matrixnok_g1_affine_from_scalar(
        const matrixnok::FpRingCandidate& ring,
        const blst::scalar::G1AffinePoint& point)
    {
        return {ring.from_scalar(point.x), ring.from_scalar(point.y)};
    }

    matrixnok::G2AffineCandidate matrixnok_g2_affine_from_scalar(
        const matrixnok::FpRingCandidate& ring,
        const blst::scalar::G2AffinePoint& point)
    {
        return {
            matrixnok_generic_fp2_from_scalar(ring, point.x),
            matrixnok_generic_fp2_from_scalar(ring, point.y)};
    }

    TEST(AvxMatrixNoKPrimitives, MillerAndPairingMatchScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::G2Fp2RingCandidate fp2ring(ring);
        const matrixnok::Fp12Ring fp12ring{};
        const matrixnok::PairingInverter inverter{};

        const auto g1_generator = blst::scalar::G1Point::generator();
        blst::scalar::G1Point g1_double = {};
        blst::scalar::G1Point g1_triple = {};
        blst::scalar::G1Point::doubling(g1_double, g1_generator);
        blst::scalar::G1Point::dadd(g1_triple, g1_double, g1_generator);
        const auto g2_generator = blst::scalar::G2Point::generator();
        blst::scalar::G2Point g2_double = {};
        blst::scalar::G2Point g2_triple = {};
        blst::scalar::G2Point::doubling(g2_double, g2_generator);
        blst::scalar::G2Point::dadd(g2_triple, g2_double, g2_generator);

        std::array<blst::scalar::G1AffinePoint, 3> g1 = {};
        std::array<blst::scalar::G2AffinePoint, 3> g2 = {};
        blst::scalar::G1Point::to_affine(g1[0], g1_generator);
        blst::scalar::G1Point::to_affine(g1[1], g1_double);
        blst::scalar::G1Point::to_affine(g1[2], g1_triple);
        blst::scalar::G2Point::to_affine(g2[0], g2_generator);
        blst::scalar::G2Point::to_affine(g2[1], g2_double);
        blst::scalar::G2Point::to_affine(g2[2], g2_triple);

        const std::array<std::pair<std::size_t, std::size_t>, 4> cases = {{
            {0, 0}, {1, 0}, {0, 1}, {2, 2},
        }};
        std::array<blst::avx::MillerResult, cases.size()> batch_inputs{};
        std::array<blst::avx::GT, cases.size()> individual_outputs{};
        for (std::size_t case_index = 0; case_index < cases.size(); ++case_index) {
            const auto [g1_index, g2_index] = cases[case_index];
            SCOPED_TRACE(
                "g1=" + std::to_string(g1_index) +
                " g2=" + std::to_string(g2_index));
            const auto candidate_g1 =
                matrixnok_g1_affine_from_scalar(ring, g1[g1_index]);
            const auto candidate_g2 =
                matrixnok_g2_affine_from_scalar(ring, g2[g2_index]);

            const matrixnok::MillerEngine miller(
                candidate_g1, candidate_g2, ring);
            const auto candidate_miller =
                miller.miller_loop(ring, fp2ring, fp12ring);
            batch_inputs[case_index] = blst::avx::MillerResult(candidate_miller);
            blst::avx::GT::final_exp(
                individual_outputs[case_index], batch_inputs[case_index]);
            blst::Fp12 candidate_miller_final = {};
            blst::Fp12::final_exp(
                candidate_miller_final,
                matrixnok_fp12_to_scalar(ring, candidate_miller));

            blst::scalar::MillerResult expected_miller = {};
            blst::scalar::GT::miller_loop(
                expected_miller, &g2[g2_index], &g1[g1_index]);
            blst::scalar::GT expected = {};
            blst::scalar::GT::final_exp(expected, expected_miller);
            expect_same_scalar_fp12(candidate_miller_final, expected.data);

            const auto candidate_pairing = blst::avx::pairing(
                candidate_g1, candidate_g2, inverter,
                fp2ring, fp12ring, ring);
            expect_same_scalar_fp12(
                matrixnok_fp12_to_scalar(ring, candidate_pairing),
                expected.data);
        }

        std::array<blst::avx::GT, cases.size()> batch_outputs{};
        blst::avx::GT::FinalExpScratch final_exp_scratch;
        EXPECT_FALSE(final_exp_scratch.has_capacity(batch_inputs.size()));
        final_exp_scratch.reserve(batch_inputs.size());
        EXPECT_TRUE(final_exp_scratch.has_capacity(batch_inputs.size()));
        blst::avx::GT::final_exp_n(
            batch_outputs.data(),
            batch_inputs.data(),
            batch_inputs.size(),
            final_exp_scratch);

        EXPECT_TRUE(final_exp_scratch.has_capacity(batch_inputs.size()));
        std::array<blst::avx::GT, cases.size()> repeated_batch_outputs{};
        std::array<blst::avx::GT, 2> prefix_batch_outputs{};
        blst::avx::GT::final_exp_n(
            prefix_batch_outputs.data(),
            batch_inputs.data(),
            prefix_batch_outputs.size(),
            final_exp_scratch);
        blst::avx::GT::final_exp_n(
            repeated_batch_outputs.data(),
            batch_inputs.data(),
            batch_inputs.size(),
            final_exp_scratch);
        blst::avx::GT::final_exp_n(nullptr, nullptr, 0, final_exp_scratch);

        EXPECT_TRUE(final_exp_scratch.has_capacity(batch_inputs.size()));
        for (std::size_t case_index = 0; case_index < cases.size(); ++case_index) {
            SCOPED_TRACE("batch case " + std::to_string(case_index));
            expect_same_scalar_fp12(
                matrixnok_fp12_to_scalar(ring, batch_outputs[case_index].data),
                matrixnok_fp12_to_scalar(ring, individual_outputs[case_index].data));
            expect_same_scalar_fp12(
                matrixnok_fp12_to_scalar(ring, repeated_batch_outputs[case_index].data),
                matrixnok_fp12_to_scalar(ring, batch_outputs[case_index].data));
            if (case_index < prefix_batch_outputs.size()) {
                expect_same_scalar_fp12(
                    matrixnok_fp12_to_scalar(ring, prefix_batch_outputs[case_index].data),
                    matrixnok_fp12_to_scalar(ring, batch_outputs[case_index].data));
            }
        }
    }
} // namespace blst::avx::test
