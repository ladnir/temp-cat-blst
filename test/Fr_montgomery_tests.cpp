#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <limits>

#include "blst/fields/Fr.h"
#include "blst/blst_t.hpp"

using blst::Fr;
using blst::limb_t;

namespace {

using FrLimbs = std::array<limb_t, Fr::num_limbs()>;

extern const blst::vec_256 NATIVE_FR_MODULUS = { BLST_r_LIMBS };
extern const blst::vec_256 NATIVE_FR_RR = { BLST_rRR_LIMBS };
extern const blst::vec_256 NATIVE_FR_ONE = {
    TO_LIMB_T(0x00000001fffffffe),
    TO_LIMB_T(0x5884b7fa00034802),
    TO_LIMB_T(0x998c4fefecbc4ff5),
    TO_LIMB_T(0x1824b159acc5056f),
};

using NativeFr = blst_256_t<
    255, NATIVE_FR_MODULUS, BLST_r0, NATIVE_FR_RR, NATIVE_FR_ONE>;

Fr make_fr(const FrLimbs& limbs)
{
    Fr value = {};
    for (size_t i = 0; i < limbs.size(); ++i)
        value[i] = limbs[i];
    return value;
}

FrLimbs limbs_of(const Fr& value)
{
    FrLimbs limbs = {};
    for (size_t i = 0; i < limbs.size(); ++i)
        limbs[i] = value[i];
    return limbs;
}

FrLimbs limbs_of(const NativeFr& value)
{
    FrLimbs limbs = {};
    value.store(limbs.data());
    return limbs;
}

bool greater_than_or_equal(const Fr& lhs, const Fr& rhs)
{
    for (size_t i = Fr::num_limbs(); i-- > 0;) {
        if (lhs[i] != rhs[i])
            return lhs[i] > rhs[i];
    }
    return true;
}

TEST(FrMontgomery, NativeConditionalSelectionNormalizesTruthyIntegers)
{
    NativeFr selector{0};
    const NativeFr a{5};
    const NativeFr b{7};
    const NativeFr zero{0};

    EXPECT_EQ(limbs_of(selector.csel(a, b, 0)), limbs_of(b));
    EXPECT_EQ(limbs_of(selector.csel(a, b, 1)), limbs_of(a));
    EXPECT_EQ(limbs_of(selector.csel(a, b, 2)), limbs_of(a));
    EXPECT_EQ(limbs_of(selector.csel(a, b, -1)), limbs_of(a));
    EXPECT_EQ(limbs_of(czero(a, 0)), limbs_of(a));
    EXPECT_EQ(limbs_of(czero(a, 2)), limbs_of(zero));
}

Fr subtract_no_underflow(const Fr& lhs, const Fr& rhs)
{
    Fr result = {};
    limb_t borrow = 0;

    for (size_t i = 0; i < Fr::num_limbs(); ++i) {
        const limb_t lhs_limb = lhs[i];
        const limb_t rhs_limb = rhs[i];
        result[i] = lhs_limb - rhs_limb - borrow;
        borrow = static_cast<limb_t>(
            lhs_limb < rhs_limb || (borrow != 0 && lhs_limb == rhs_limb));
    }

    EXPECT_EQ(borrow, 0);
    return result;
}

Fr add_without_overflow(const Fr& lhs, const Fr& rhs)
{
    Fr result = {};
    limb_t carry = 0;

    for (size_t i = 0; i < Fr::num_limbs(); ++i) {
        const limb_t partial = lhs[i] + rhs[i];
        const limb_t carry_from_operands = static_cast<limb_t>(partial < lhs[i]);
        const limb_t sum = partial + carry;
        const limb_t carry_from_carry = static_cast<limb_t>(sum < partial);
        result[i] = sum;
        carry = carry_from_operands | carry_from_carry;
    }

    EXPECT_EQ(carry, 0);
    return result;
}

Fr add_small_without_overflow(const Fr& value, limb_t increment)
{
    Fr small = {};
    small[0] = increment;
    return add_without_overflow(value, small);
}

Fr subtract_small(const Fr& value, limb_t decrement)
{
    Fr small = {};
    small[0] = decrement;
    return subtract_no_underflow(value, small);
}

Fr reduce_reference(Fr value)
{
    const Fr& modulus = blst::get_BLST_r();
    while (greater_than_or_equal(value, modulus))
        value = subtract_no_underflow(value, modulus);
    return value;
}

Fr to_mont_reference(const Fr& value)
{
    const Fr reduced = reduce_reference(value);
    Fr result = {};

    // The sparse assembly contract accepts a non-reduced value only as the
    // second input (the third function argument). Keep fixed R^2 first.
    Fr::mul(result, blst::get_BLST_rRR(), reduced);
    return result;
}

Fr next_wide_input(uint64_t& state)
{
    Fr value = {};
    for (size_t i = 0; i < Fr::num_limbs(); ++i) {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        const uint64_t word = state * UINT64_C(0x2545f4914f6cdd1d);
        value[i] = static_cast<limb_t>(word);
    }
    return value;
}

std::array<Fr, 10> boundary_inputs()
{
    const Fr& modulus = blst::get_BLST_r();
    const Fr twice_modulus = add_without_overflow(modulus, modulus);

    Fr alternating = {};
    for (size_t i = 0; i < Fr::num_limbs(); ++i) {
        alternating[i] = (i & 1) == 0
            ? std::numeric_limits<limb_t>::max()
            : static_cast<limb_t>(0x5555555555555555ULL);
    }

    Fr maximum = {};
    for (size_t i = 0; i < Fr::num_limbs(); ++i)
        maximum[i] = std::numeric_limits<limb_t>::max();

    return {
        Fr{},
        add_small_without_overflow(Fr{}, 1),
        subtract_small(modulus, 1),
        modulus,
        add_small_without_overflow(modulus, 1),
        subtract_small(twice_modulus, 1),
        twice_modulus,
        add_small_without_overflow(twice_modulus, 1),
        alternating,
        maximum,
    };
}

} // namespace

TEST(FrMontgomery, ConvertsReducedAndUnreducedBoundaryInputs)
{
    for (const Fr& input : boundary_inputs()) {
        SCOPED_TRACE(::testing::PrintToString(limbs_of(input)));

        Fr actual = {};
        Fr::to_mont(actual, input);

        EXPECT_EQ(limbs_of(actual), limbs_of(to_mont_reference(input)));

        Fr round_trip = {};
        Fr::from_mont(round_trip, actual);
        EXPECT_EQ(limbs_of(round_trip), limbs_of(reduce_reference(input)));
    }
}

TEST(FrMontgomery, SupportsInputOutputAliasing)
{
    for (const Fr& input : boundary_inputs()) {
        SCOPED_TRACE(::testing::PrintToString(limbs_of(input)));

        const Fr expected = to_mont_reference(input);
        Fr in_place = input;
        Fr::to_mont(in_place, in_place);

        EXPECT_EQ(limbs_of(in_place), limbs_of(expected));
    }
}

TEST(FrMontgomery, MatchesSafeOperandOrderForWideInputs)
{
    uint64_t state = UINT64_C(0x6a09e667f3bcc909);

    for (size_t i = 0; i < 4096; ++i) {
        const Fr input = next_wide_input(state);
        SCOPED_TRACE(i);
        SCOPED_TRACE(::testing::PrintToString(limbs_of(input)));

        Fr actual = {};
        Fr::to_mont(actual, input);

        EXPECT_EQ(limbs_of(actual), limbs_of(to_mont_reference(input)));
    }
}

TEST(FrMontgomery, NativeCppToMontHandlesUnreducedInputs)
{
    for (const Fr& input : boundary_inputs()) {
        SCOPED_TRACE(::testing::PrintToString(limbs_of(input)));

        NativeFr actual{input.data};
        actual.to();

        Fr expected = {};
        Fr::to_mont(expected, input);
        EXPECT_EQ(limbs_of(actual), limbs_of(expected));
    }
}

TEST(FrMontgomery, NativeCppReciprocalMatchesScalarField)
{
    const Fr& modulus = blst::get_BLST_r();
    const std::array<Fr, 4> inputs = {
        add_small_without_overflow(Fr{}, 1),
        add_small_without_overflow(Fr{}, 2),
        subtract_small(modulus, 1),
        subtract_small(modulus, 2),
    };

    for (const Fr& input : inputs) {
        SCOPED_TRACE(::testing::PrintToString(limbs_of(input)));

        NativeFr native_input{input.data};
        native_input.to();
        const NativeFr actual = native_input.reciprocal();

        Fr montgomery_input = {};
        Fr expected = {};
        Fr::to_mont(montgomery_input, input);
        Fr::reciprocal(expected, montgomery_input);
        EXPECT_EQ(limbs_of(actual), limbs_of(expected));
    }
}
