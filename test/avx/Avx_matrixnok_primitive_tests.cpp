#include "blst/avx/generated_matrixnok_constants.hpp"
#include "blst/avx/cpu/vector/conversion.hpp"
#include "blst/avx/rns/elementwise_reduce.h"
#include "blst/avx/rns/matrixnok_change_base.h"
#include "blst/avx/rns/matrixnok_fp.h"
#include "blst/avx/fields/matrixnok_fp2.h"
#include "blst/avx/fields/matrixnok_fp12.h"
#include "blst/avx/gt/glv4.h"
#include "blst/avx/curve/matrixnok_g1.h"
#include "blst/common/bytes.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/primitveRoot.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fr.h"
#include "blst/scalar/curve/point.h"

#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

namespace blst::avx::test
{
    namespace matrixnok_constants = matrixnok::generated;

    static_assert(ConvertToRNS<381, 8>::indigits == 8);
    static_assert(ConvertToRNS<381, 8, 50>::indigits == 8);
    static_assert(ConvertToRNS<401, 8>::indigits == 8);
    static_assert(ConvertToRNS<401, 8, 50>::indigits == 9);
    static_assert(matrixnok_constants::BLS12_381_RNS_DIM == 8);
    static_assert(matrixnok_constants::BLS12_381_CONVERT_TO_RADIX_BITS == 50);
    static_assert(matrixnok_constants::moduli1.size() == 8);
    static_assert(matrixnok_constants::moduli2.size() == 8);
    static_assert(matrixnok_constants::r1_perm_mat_fixed.size() == 8);
    static_assert(std::is_empty_v<matrixnok::FpRingCandidate>);
    static_assert(std::is_empty_v<matrixnok::Fp2RingCandidate>);
    static_assert(std::is_empty_v<matrixnok::Fp12RingCandidate>);
    static_assert(std::is_empty_v<matrixnok::G1CurveCandidate>);

    using MatrixNoKFp12Scaffold = matrixnok::Fp12RingCandidate;
    static_assert(std::is_default_constructible_v<
        typename MatrixNoKFp12Scaffold::StandardElement>);

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

    TEST(AvxMatrixNoKPrimitives, StatelessFixedContextMatchesGenericHelpers)
    {
        EXPECT_EQ(sizeof(matrixnok::FpRingCandidate), 1U);
        EXPECT_EQ(sizeof(matrixnok::Fp2RingCandidate), 1U);

        const RNSMatrix<8, 8> generic_r2(
            matrixnok_constants::r2_rns_mat,
            matrixnok_constants::r2_correction,
            matrixnok_constants::r2_shifted_quotient_estimations,
            matrixnok_constants::moduli1);
        const RNSMatrix<1, 8, 8, 1, 50> generic_convert_to(
            matrixnok_constants::convert_to_rns_mat,
            matrixnok_constants::convert_to_correction,
            matrixnok_constants::convert_to_shifted_quotient_estimations,
            matrixnok_constants::moduli2);
        const RNSMatrix<8, 1, 1, 8> generic_convert_from(
            matrixnok_constants::convert_from_rns_mat,
            matrixnok_constants::convert_from_correction,
            matrixnok_constants::convert_from_shifted_quotient_estimations,
            std::array<uint64_t, 1>{1});

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

        const auto compare_shifted = [&](const auto& generic, auto fixed_reduce, const auto& input) {
            AVXVector<8> generic_high(initial_high);
            AVXVector<8> generic_low(initial_low);
            AVXVector<8> fixed_high(initial_high);
            AVXVector<8> fixed_low(initial_low);
            generic.rns_reduce_raw_shifted(
                AVXVector<8>(input), generic_high, generic_low);
            fixed_reduce(AVXVector<8>(input), fixed_high, fixed_low);
            EXPECT_EQ(fixed_high.to_array(), generic_high.to_array());
            EXPECT_EQ(fixed_low.to_array(), generic_low.to_array());
        };
        compare_shifted(
            generic_r2,
            matrixnok::detail::R2Matrix::reduce_raw_shifted,
            r2_input);
        compare_shifted(
            generic_convert_to,
            matrixnok::detail::ConvertToMatrix::reduce_raw_shifted,
            digit_input);

        AVXVector<8> generic_from_high(initial_high);
        AVXVector<8> generic_from_low(initial_low);
        AVXVector<8> fixed_from_high(initial_high);
        AVXVector<8> fixed_from_low(initial_low);
        generic_convert_from.rns_reduce_raw(
            AVXVector<8>(r2_input), generic_from_high, generic_from_low);
        matrixnok::detail::ConvertFromMatrix::reduce_raw(
            AVXVector<8>(r2_input), fixed_from_high, fixed_from_low);
        EXPECT_EQ(fixed_from_high.to_array(), generic_from_high.to_array());
        EXPECT_EQ(fixed_from_low.to_array(), generic_from_low.to_array());

        std::array<AVXVector<8>, 2> batch_input = {
            AVXVector<8>(r2_input), AVXVector<8>(digit_input)};
        std::array<AVXVector<8>, 2> generic_batch_high = {
            AVXVector<8>(initial_high), AVXVector<8>(initial_low)};
        std::array<AVXVector<8>, 2> generic_batch_low = {
            AVXVector<8>(initial_low), AVXVector<8>(initial_high)};
        auto fixed_batch_high = generic_batch_high;
        auto fixed_batch_low = generic_batch_low;
        generic_r2.rns_reduce_raw_batch_shifted<2>(
            &batch_input, generic_batch_high, generic_batch_low);
        matrixnok::detail::R2Matrix::reduce_raw_batch_shifted<2>(
            batch_input, fixed_batch_high, fixed_batch_low);
        for (std::size_t item = 0; item < 2; ++item) {
            EXPECT_EQ(fixed_batch_high[item].to_array(), generic_batch_high[item].to_array());
            EXPECT_EQ(fixed_batch_low[item].to_array(), generic_batch_low[item].to_array());
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
        const MontgomeryReduce2<8, 50, 6> generic_reducer(matrixnok_constants::moduli1);
        const auto generic_reduced = generic_reducer.template reduce_auto<2>(
            generic_reducer.reduce_wide(wide));
        const auto fixed_reduced = matrixnok::detail::ReducerM1::template reduce_auto<2>(
            matrixnok::detail::ReducerM1::reduce_wide(wide));
        EXPECT_EQ(fixed_reduced.data.to_array(), generic_reduced.data.to_array());
        EXPECT_EQ(
            matrixnok::detail::ReducerM1::template moduli<1>().data.to_array(),
            generic_reducer.template moduli<1>().data.to_array());
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

    template<int64_t upper_bound, int64_t stop, int log_multiples>
    void check_min_sub_schedule()
    {
        constexpr int limbs = 8;
        constexpr int element_bits = 50;
        const std::array<uint64_t, limbs> moduli = {
            (uint64_t{1} << element_bits) - 1,
            (uint64_t{1} << element_bits) - 3,
            (uint64_t{1} << element_bits) - 5,
            (uint64_t{1} << element_bits) - 7,
            (uint64_t{1} << element_bits) - 11,
            (uint64_t{1} << element_bits) - 13,
            (uint64_t{1} << element_bits) - 17,
            (uint64_t{1} << element_bits) - 19,
        };
        const MontgomeryReduce2<limbs, element_bits, log_multiples> reducer(moduli);

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
            const auto reduced = reducer.template min_sub_reduce<upper_bound, stop>(Input(AVXVector<limbs>(values)));
            const auto actual = reduced.data.to_array();
            for (int lane = 0; lane < limbs; ++lane) {
                EXPECT_LT(actual[lane], static_cast<uint64_t>(stop) * moduli[lane]);
                EXPECT_EQ(actual[lane] % moduli[lane], values[lane] % moduli[lane]);
            }
        }
    }

    TEST(AvxMatrixNoKPrimitives, ConditionalSubtractionCoversConfiguredGap)
    {
        check_min_sub_schedule<2, 1, 1>();
        check_min_sub_schedule<4, 1, 2>();
        check_min_sub_schedule<40, 1, 6>();
        check_min_sub_schedule<100, 95, 3>();
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
        for (std::size_t index = 0; index < batch; ++index) {
            expect_standard_residue_bounds(ring, products[index]);
            expect_same_scalar_fp(matrixnok_to_scalar(ring, products[index]), expected[index]);
        }
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
        const matrixnok::G1GLVCandidate glv(curve, ring);
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

    struct MatrixNoKG1FftCtx
    {
        inline static const matrixnok::FpRingCandidate ring = {};
        inline static const matrixnok::G1CurveCandidate curve = {};
        inline static const matrixnok::G1GLVCandidate glv = {curve, ring};

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

    matrixnok::Fp2RingCandidate::Element matrixnok_fp2_from_scalar(
        const matrixnok::FpRingCandidate& ring,
        const blst::Fp2& value)
    {
        return {
            ring.from_canonical_50bit_digits(canonical_50bit_digits(value[0])),
            ring.from_canonical_50bit_digits(canonical_50bit_digits(value[1]))};
    }

    blst::Fp2 matrixnok_fp2_to_scalar(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::Fp2RingCandidate::Element& value)
    {
        return {{matrixnok_to_scalar(ring, value.c0), matrixnok_to_scalar(ring, value.c1)}};
    }

    void expect_same_scalar_fp2(const blst::Fp2& actual, const blst::Fp2& expected)
    {
        expect_same_scalar_fp(actual[0], expected[0]);
        expect_same_scalar_fp(actual[1], expected[1]);
    }

    void expect_standard_fp2_bounds(
        const matrixnok::FpRingCandidate& ring,
        const matrixnok::Fp2RingCandidate::Element& value)
    {
        expect_standard_residue_bounds(ring, value.c0);
        expect_standard_residue_bounds(ring, value.c1);
    }

    TEST(AvxMatrixNoKPrimitives, Fp2MultiplyAndSquareMatchScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const matrixnok::Fp2RingCandidate fp2ring;
        testing::Test::RecordProperty("fp_ring_candidate_bytes", sizeof(ring));
        testing::Test::RecordProperty("fp2_ring_candidate_bytes", sizeof(fp2ring));

        const auto base_values = matrixnok_scalar_values();
        std::array<blst::Fp2, 8> values = {};
        for (std::size_t index = 0; index < values.size(); ++index) {
            values[index][0] = base_values[index];
            values[index][1] = base_values[(3 * index + 1) % base_values.size()];
        }

        for (const blst::Fp2& value : values) {
            const auto encoded = matrixnok_fp2_from_scalar(ring, value);
            const auto candidate_square = fp2ring.square(encoded);
            blst::Fp2 expected_square = {};
            blst::Fp2::sqr(expected_square, value);
            expect_standard_fp2_bounds(ring, candidate_square);
            expect_same_scalar_fp2(
                matrixnok_fp2_to_scalar(ring, candidate_square), expected_square);

            auto square_alias = encoded;
            fp2ring.square(square_alias, square_alias);
            expect_same_scalar_fp2(
                matrixnok_fp2_to_scalar(ring, square_alias), expected_square);
        }

        for (const blst::Fp2& left : values) {
            for (const blst::Fp2& right : values) {
                const auto candidate_left = matrixnok_fp2_from_scalar(ring, left);
                const auto candidate_right = matrixnok_fp2_from_scalar(ring, right);
                const auto candidate_product = fp2ring.multiply(candidate_left, candidate_right);
                blst::Fp2 expected_product = {};
                blst::Fp2::mul(expected_product, left, right);
                expect_standard_fp2_bounds(ring, candidate_product);
                expect_same_scalar_fp2(
                    matrixnok_fp2_to_scalar(ring, candidate_product), expected_product);
            }
        }

        auto right_alias = matrixnok_fp2_from_scalar(ring, values[5]);
        const auto left_for_alias = matrixnok_fp2_from_scalar(ring, values[7]);
        fp2ring.multiply(right_alias, left_for_alias, right_alias);
        blst::Fp2 expected_alias = {};
        blst::Fp2::mul(expected_alias, values[7], values[5]);
        expect_standard_fp2_bounds(ring, right_alias);
        expect_same_scalar_fp2(matrixnok_fp2_to_scalar(ring, right_alias), expected_alias);

        auto left_alias = matrixnok_fp2_from_scalar(ring, values[4]);
        const auto right_for_alias = matrixnok_fp2_from_scalar(ring, values[6]);
        fp2ring.multiply(left_alias, left_alias, right_for_alias);
        blst::Fp2 expected_left_alias = {};
        blst::Fp2::mul(expected_left_alias, values[4], values[6]);
        expect_standard_fp2_bounds(ring, left_alias);
        expect_same_scalar_fp2(
            matrixnok_fp2_to_scalar(ring, left_alias), expected_left_alias);

        auto candidate_product_chain = matrixnok_fp2_from_scalar(ring, values[4]);
        blst::Fp2 scalar_product_chain = values[4];
        const auto candidate_factor = matrixnok_fp2_from_scalar(ring, values[6]);
        for (int iteration = 0; iteration < 64; ++iteration) {
            candidate_product_chain = fp2ring.multiply(candidate_product_chain, candidate_factor);
            blst::Fp2 next = {};
            blst::Fp2::mul(next, scalar_product_chain, values[6]);
            scalar_product_chain = next;
            expect_standard_fp2_bounds(ring, candidate_product_chain);
            expect_same_scalar_fp2(
                matrixnok_fp2_to_scalar(ring, candidate_product_chain), scalar_product_chain);
        }

        auto candidate_square_chain = matrixnok_fp2_from_scalar(ring, values[5]);
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

    MatrixNoKFp12Scaffold::StandardElement matrixnok_fp12_from_scalar(
        const matrixnok::FpRingCandidate& ring, const blst::Fp12& value)
    {
        MatrixNoKFp12Scaffold::StandardElement out = {};
        for (std::size_t coefficient = 0; coefficient < out.size(); ++coefficient) {
            out[coefficient] = ring.from_scalar(
                value[coefficient / 6][(coefficient % 6) / 2][coefficient % 2]);
        }
        return out;
    }

    blst::Fp12 matrixnok_fp12_to_scalar(
        const matrixnok::FpRingCandidate& ring,
        const MatrixNoKFp12Scaffold::StandardElement& value)
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

    TEST(AvxMatrixNoKPrimitives, Fp12ScaffoldMultiplyAndSquareMatchScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const MatrixNoKFp12Scaffold fp12ring(ring);
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

    TEST(AvxMatrixNoKPrimitives, Fp12MapsAndGtExponentiationMatchScalar)
    {
        const matrixnok::FpRingCandidate ring;
        const MatrixNoKFp12Scaffold fp12ring(ring);
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
            MatrixNoKFp12Scaffold::StandardElement actual = {};
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
} // namespace blst::avx::test
