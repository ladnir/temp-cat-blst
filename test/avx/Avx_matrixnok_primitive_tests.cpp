#include "blst/avx/generated_matrixnok_constants.hpp"
#include "blst/avx/cpu/vector/conversion.hpp"
#include "blst/avx/rns/elementwise_reduce.h"
#include "blst/avx/rns/matrixnok_change_base.h"
#include "blst/avx/rns/matrixnok_fp.h"
#include "blst/fields/Fp.h"

#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <cstring>

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
        auto p_minus_one = std::to_array<uint64_t>({
            blst::get_BLST_P().l[0], blst::get_BLST_P().l[1], blst::get_BLST_P().l[2],
            blst::get_BLST_P().l[3], blst::get_BLST_P().l[4], blst::get_BLST_P().l[5]});
        --p_minus_one[0];
        const std::array<blst::Fp, 8> values = {
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
} // namespace blst::avx::test
