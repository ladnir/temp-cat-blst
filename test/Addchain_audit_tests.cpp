/*
 * Copyright Supranational LLC and Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include <gtest/gtest.h>

#include <array>

#include "blst/scalar/fields/Fp.h"
#include "blst/scalar/fields/Fp2.h"

namespace {

using blst::byte;
using blst::limb_t;
using blst::scalar::Fp;
using blst::scalar::Fp2;

Fp fp_from_u64(limb_t value)
{
    Fp raw = {};
    Fp out;
    raw[0] = value;
    Fp::to_mont(out, raw);
    return out;
}

constexpr byte p_minus_2[] = {
    TO_BYTES(0xb9feffffffffaaa9), TO_BYTES(0x1eabfffeb153ffff),
    TO_BYTES(0x6730d2a0f6b0f624), TO_BYTES(0x64774b84f38512bf),
    TO_BYTES(0x4b1ba7b6434bacd7), TO_BYTES(0x1a0111ea397fe69a)
};

constexpr byte p_minus_3_div_4[] = {
    TO_BYTES(0xee7fbfffffffeaaa), TO_BYTES(0x07aaffffac54ffff),
    TO_BYTES(0xd9cc34a83dac3d89), TO_BYTES(0xd91dd2e13ce144af),
    TO_BYTES(0x92c6e9ed90d2eb35), TO_BYTES(0x0680447a8e5ff9a6)
};

std::array<Fp, 5> boundary_values()
{
    Fp minus_one;
    Fp::neg(minus_one, Fp::one());
    return {Fp::zero(), Fp::one(), fp_from_u64(2), fp_from_u64(3), minus_one};
}

TEST(AddchainAudit, ReciprocalMatchesIndependentExponentiation)
{
    for (const Fp& input : boundary_values()) {
        Fp expected;
        Fp actual;
        Fp::exp(expected, input, p_minus_2, 381);
        Fp::flt_reciprocal(actual, input);
        EXPECT_TRUE(Fp::equal(actual, expected));

        Fp aliased = input;
        Fp::flt_reciprocal(aliased, aliased);
        EXPECT_TRUE(Fp::equal(aliased, expected));

        if (!Fp::equal(input, Fp::zero())) {
            Fp product;
            Fp::mul(product, actual, input);
            EXPECT_TRUE(Fp::equal(product, Fp::one()));
        }
    }
}

TEST(AddchainAudit, ReciprocalSquareRootMatchesIndependentExponentiation)
{
    for (const Fp& input : boundary_values()) {
        Fp expected;
        Fp actual;
        Fp::exp(expected, input, p_minus_3_div_4, 379);
        Fp::recip_sqrt_3mod4(actual, input);
        EXPECT_TRUE(Fp::equal(actual, expected));

        Fp aliased = input;
        Fp::recip_sqrt_3mod4(aliased, aliased);
        EXPECT_TRUE(Fp::equal(aliased, expected));
    }
}

TEST(AddchainAudit, ZeroSemanticsMatchBranchlessInfinityUse)
{
    Fp fp_inverse = Fp::one();
    Fp::flt_reciprocal(fp_inverse, Fp::zero());
    EXPECT_TRUE(Fp::equal(fp_inverse, Fp::zero()));

    Fp fp_recip_sqrt = Fp::one();
    EXPECT_TRUE(Fp::recip_sqrt(fp_recip_sqrt, Fp::zero()));
    EXPECT_TRUE(Fp::equal(fp_recip_sqrt, Fp::zero()));

    Fp2 fp2_inverse = Fp2::one();
    Fp2::flt_reciprocal(fp2_inverse, Fp2::zero());
    EXPECT_TRUE(blst::vec_is_equal(
        &fp2_inverse, &Fp2::zero(), sizeof(fp2_inverse)));
}

} // namespace
