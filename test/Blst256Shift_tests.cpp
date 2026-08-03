#include "blst/blst_t.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

#include <gtest/gtest.h>

namespace
{

extern const blst::vec_256 SHIFT_MODULUS = {
    UINT64_C(0xffffffffffffffed),
    UINT64_C(0xffffffffffffffff),
    UINT64_C(0xffffffffffffffff),
    UINT64_C(0x7fffffffffffffff),
};
extern const blst::vec_256 UNUSED_RR = { 0, 0, 0, 0 };
extern const blst::vec_256 UNUSED_ONE = { 1, 0, 0, 0 };
extern const blst::vec_384 UNUSED_384_MODULUS = { 1, 0, 0, 0, 0, 0 };
extern const blst::vec_384 UNUSED_384_RR = { 0, 0, 0, 0, 0, 0 };
extern const blst::vec_384 UNUSED_384_ONE = { 1, 0, 0, 0, 0, 0 };

using ShiftElement = blst_256_t<255, SHIFT_MODULUS, 0, UNUSED_RR, UNUSED_ONE>;
using Header384Element = blst_384_t<
    381, UNUSED_384_MODULUS, 0, UNUSED_384_RR, UNUSED_384_ONE>;
using Limbs = std::array<uint64_t, 4>;
using Limbs384 = std::array<uint64_t, 6>;

static_assert(sizeof(ShiftElement) == sizeof(blst::vec_256));
static_assert(ShiftElement::nbits == 255);
static_assert(ShiftElement::degree == 1);
static_assert(sizeof(Header384Element) == sizeof(blst::vec_384));
static_assert(Header384Element::nbits == 381);
static_assert(Header384Element::degree == 1);

Limbs store(const ShiftElement& value)
{
    Limbs ret{};
    value.store(ret.data());
    return ret;
}

Limbs384 store(const Header384Element& value)
{
    Limbs384 ret{};
    value.store(ret.data());
    return ret;
}

Limbs modular_right_shift(Limbs value, unsigned count)
{
    constexpr Limbs modulus = {
        UINT64_C(0xffffffffffffffed),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x7fffffffffffffff),
    };

    while (count--) {
        if ((value[0] & 1) != 0) {
            uint64_t carry = 0;
            for (size_t i = 0; i < value.size(); ++i) {
                const uint64_t sum = value[i] + modulus[i];
                const uint64_t carry_from_sum = sum < value[i];
                const uint64_t with_carry = sum + carry;
                const uint64_t carry_from_input = with_carry < sum;
                value[i] = with_carry;
                carry = carry_from_sum | carry_from_input;
            }
            EXPECT_EQ(carry, 0u);
        }

        uint64_t high_bits = 0;
        for (size_t i = value.size(); i-- > 0;) {
            const uint64_t next_high_bits = value[i] << 63;
            value[i] = (value[i] >> 1) | high_bits;
            high_bits = next_high_bits;
        }
    }

    return value;
}

ShiftElement make_element(const Limbs& limbs)
{
    blst::vec_256 value = { limbs[0], limbs[1], limbs[2], limbs[3] };
    return ShiftElement{value};
}

Header384Element make_element(const Limbs384& limbs)
{
    blst::vec_384 value = {
        limbs[0], limbs[1], limbs[2], limbs[3], limbs[4], limbs[5]
    };
    return Header384Element{value};
}

TEST(Blst256Shift, LeftShiftByZeroIsIdentity)
{
    const Limbs input = {
        UINT64_C(0x0123456789abcdef),
        UINT64_C(0xfedcba9876543210),
        UINT64_C(0x1111111111111111),
        UINT64_C(0x2222222222222222),
    };
    const ShiftElement value = make_element(input);

    EXPECT_EQ(store(value << 0), input);

    ShiftElement assigned = value;
    assigned <<= 0;
    EXPECT_EQ(store(assigned), input);
}

TEST(Blst256Shift, RightShiftReturnsModularQuotientWithoutMutatingInput)
{
    const Limbs input = { 5, 0, 0, 0 };
    const ShiftElement value = make_element(input);

    const ShiftElement shifted = value >> 1;

    EXPECT_EQ(store(shifted), modular_right_shift(input, 1));
    EXPECT_EQ(store(value), input);
}

TEST(Blst256Shift, RightShiftByZeroIsIdentity)
{
    const Limbs input = {
        UINT64_C(0x0123456789abcdef),
        UINT64_C(0xfedcba9876543210),
        UINT64_C(0x1111111111111111),
        UINT64_C(0x2222222222222222),
    };
    const ShiftElement value = make_element(input);

    EXPECT_EQ(store(value >> 0), input);

    ShiftElement assigned = value;
    assigned >>= 0;
    EXPECT_EQ(store(assigned), input);
}

TEST(Blst256Shift, RightShiftAssignmentHandlesEvenAndRepeatedOddValues)
{
    const Limbs odd_input = {
        UINT64_C(0x0123456789abcdef),
        UINT64_C(0xfedcba9876543210),
        UINT64_C(0x1111111111111111),
        UINT64_C(0x2222222222222222),
    };
    ShiftElement odd_value = make_element(odd_input);
    odd_value >>= 7;
    EXPECT_EQ(store(odd_value), modular_right_shift(odd_input, 7));

    const Limbs even_input = {
        UINT64_C(0xfffffffffffffffe),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0),
        UINT64_C(0x4000000000000000),
    };
    ShiftElement even_value = make_element(even_input);
    even_value >>= 1;
    EXPECT_EQ(store(even_value), modular_right_shift(even_input, 1));
}

TEST(Blst384Shift, ZeroBitShiftsAreIdentity)
{
    const Limbs384 input = {
        UINT64_C(0x0123456789abcdef),
        UINT64_C(0xfedcba9876543210),
        UINT64_C(0x1111111111111111),
        UINT64_C(0x2222222222222222),
        UINT64_C(0x3333333333333333),
        UINT64_C(0x0444444444444444),
    };
    const Header384Element value = make_element(input);

    EXPECT_EQ(store(value << 0), input);
    EXPECT_EQ(store(value >> 0), input);

    Header384Element left_assigned = value;
    left_assigned <<= 0;
    EXPECT_EQ(store(left_assigned), input);

    Header384Element right_assigned = value;
    right_assigned >>= 0;
    EXPECT_EQ(store(right_assigned), input);
}

} // namespace
