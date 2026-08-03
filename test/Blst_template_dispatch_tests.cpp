#include "blst/fields/Fp.h"
#include "blst/blst_t.hpp"

#include <array>
#include <cstdint>

#include <gtest/gtest.h>

namespace
{

using blst::Fp;
using blst::limb_t;
using FpLimbs = std::array<limb_t, Fp::num_limb()>;

extern const blst::vec_384 NATIVE_FP_MODULUS = { BLST_P_LIMBS };
extern const blst::vec_384 NATIVE_FP_RR = { BLST_RR_LIMBS };
extern const blst::vec_384 NATIVE_FP_ONE = { ONE_MONT_P_LIMBS };

using NativeFp = blst_384_t<
    381, NATIVE_FP_MODULUS, BLST_p0, NATIVE_FP_RR, NATIVE_FP_ONE>;

Fp raw_fp(uint64_t value)
{
    Fp ret = {};
    ret[0] = value;
    return ret;
}

FpLimbs limbs_of(const Fp& value)
{
    FpLimbs ret = {};
    for (size_t i = 0; i < ret.size(); ++i)
        ret[i] = value[i];
    return ret;
}

FpLimbs limbs_of(const NativeFp& value)
{
    FpLimbs ret = {};
    value.store(ret.data());
    return ret;
}

NativeFp native_montgomery(const Fp& raw)
{
    NativeFp ret{raw.data};
    ret.to();
    return ret;
}

Fp scalar_montgomery(const Fp& raw)
{
    Fp ret = {};
    Fp::to_mont(ret, raw);
    return ret;
}

TEST(BlstTemplateDispatch, NativeFpMatchesSelectedScalarOperations)
{
    const Fp raw_a = raw_fp(5);
    const Fp raw_b = raw_fp(3);
    const Fp a = scalar_montgomery(raw_a);
    const Fp b = scalar_montgomery(raw_b);
    const NativeFp native_a = native_montgomery(raw_a);
    const NativeFp native_b = native_montgomery(raw_b);

    Fp expected = {};

    Fp::add(expected, a, b);
    EXPECT_EQ(limbs_of(native_a + native_b), limbs_of(expected));

    Fp::sub(expected, a, b);
    EXPECT_EQ(limbs_of(native_a - native_b), limbs_of(expected));

    Fp::mul(expected, a, b);
    EXPECT_EQ(limbs_of(native_a * native_b), limbs_of(expected));

    Fp::lshift(expected, a, 3);
    EXPECT_EQ(limbs_of(native_a << 3), limbs_of(expected));

    Fp::rshift(expected, a, 1);
    EXPECT_EQ(limbs_of(native_a >> 1), limbs_of(expected));

    Fp::cneg(expected, a, true);
    NativeFp native_negated = native_a;
    native_negated.cneg(true);
    EXPECT_EQ(limbs_of(native_negated), limbs_of(expected));

    Fp::reciprocal(expected, a);
    EXPECT_EQ(limbs_of(native_a.reciprocal()), limbs_of(expected));

    NativeFp native_round_trip = native_a;
    native_round_trip.from();
    EXPECT_EQ(limbs_of(native_round_trip), limbs_of(raw_a));
}

} // namespace
