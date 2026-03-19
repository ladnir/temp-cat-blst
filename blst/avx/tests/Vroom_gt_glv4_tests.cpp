#include <gtest/gtest.h>

#include "blst/avx/src/gt_glv4.hpp"
#include "blst/avx/tests/vroom_basefield_test_utils.hpp"

#include <array>
#include <string>

namespace {

using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomRing;

void expect_fp12_bigints_equal(
    const std::array<BigInt, 12>& actual,
    const std::array<BigInt, 12>& expected,
    const std::string& label)
{
    for (size_t i = 0; i < actual.size(); ++i) {
        SCOPED_TRACE(label + " coeff=" + std::to_string(i));
        EXPECT_EQ(actual[i].to_string(), expected[i].to_string());
    }
}

void expect_vroom_glv4_matches_scalar(
    const VroomRing& ring,
    const VroomFp12Ring& fp12ring,
    uint64_t base_seed,
    uint64_t scalar_seed,
    const std::string& label)
{
    const auto vroom_base = blst_vroom_test::make_vroom_cyclotomic_fp12_from_seed(ring, base_seed);
    const auto scalar_limbs = blst_vroom_test::seeded_scalar_limbs_le(scalar_seed);
    const auto scalar_bytes = blst_vroom_test::seeded_scalar_bytes_le(scalar_seed);

    const auto basic = vroom_cyclotomic_exp_basic(vroom_base, scalar_limbs.data(), fp12ring, ring);
    const auto runtime_glv4 = vroom_cyclotomic_exp_glv4(vroom_base, scalar_limbs.data(), fp12ring, ring);
    const VroomGLV4Exp props = vroom_glv4_compile(scalar_limbs.data());
    const auto compiled_glv4 = vroom_cyclotomic_exp_glv4(vroom_base, props, fp12ring, ring);

    const blst::Fp12 scalar_base = blst_vroom_test::make_scalar_cyclotomic_fp12_from_seed(base_seed);
    blst::Fp12 scalar_basic = {};
    blst::Fp12 scalar_glv4 = {};
    blst::Fp12::cyclotomic_exp(scalar_basic, scalar_base, scalar_bytes.data(), 256);
    blst::Fp12::cyclotomic_exp_glv4(scalar_glv4, scalar_base, scalar_limbs.data());

    const auto expected = blst_vroom_test::scalar_fp12_to_bigints(scalar_glv4);
    const auto expected_basic = blst_vroom_test::scalar_fp12_to_bigints(scalar_basic);

    const auto basic_bigints = blst_vroom_test::vroom_fp12_to_bigints(ring, basic);
    const auto runtime_bigints = blst_vroom_test::vroom_fp12_to_bigints(ring, runtime_glv4);
    const auto compiled_bigints = blst_vroom_test::vroom_fp12_to_bigints(ring, compiled_glv4);

    expect_fp12_bigints_equal(basic_bigints, expected_basic, label + " basic");
    expect_fp12_bigints_equal(runtime_bigints, expected, label + " runtime");
    expect_fp12_bigints_equal(compiled_bigints, expected, label + " compiled");
    expect_fp12_bigints_equal(runtime_bigints, compiled_bigints, label + " runtime_vs_compiled");
}

TEST(VroomGtGlv4, RuntimeAndCompiledMatchScalarBlst)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);

    for (const auto& test_case : std::array<std::pair<uint64_t, uint64_t>, 4>{{{42, 123}, {7, 99}, {456, 321}, {999, 17}}}) {
        expect_vroom_glv4_matches_scalar(ring, fp12ring, test_case.first, test_case.second,
                                         "base=" + std::to_string(test_case.first) + ",scalar=" + std::to_string(test_case.second));
    }
}

TEST(VroomGtGlv4, ZeroScalarReturnsIdentity)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const auto base = blst_vroom_test::make_vroom_cyclotomic_fp12_from_seed(ring, 42);
    const std::array<uint64_t, 4> zero_scalar = {};

    const auto basic = vroom_cyclotomic_exp_basic(base, zero_scalar.data(), fp12ring, ring);
    const auto runtime_glv4 = vroom_cyclotomic_exp_glv4(base, zero_scalar.data(), fp12ring, ring);
    const auto compiled_glv4 = vroom_cyclotomic_exp_glv4(base, vroom_glv4_compile(zero_scalar.data()), fp12ring, ring);
    const auto expected = blst_vroom_test::vroom_fp12_to_bigints(ring, fp12ring.one(ring));

    expect_fp12_bigints_equal(blst_vroom_test::vroom_fp12_to_bigints(ring, basic), expected, "zero basic");
    expect_fp12_bigints_equal(blst_vroom_test::vroom_fp12_to_bigints(ring, runtime_glv4), expected, "zero runtime");
    expect_fp12_bigints_equal(blst_vroom_test::vroom_fp12_to_bigints(ring, compiled_glv4), expected, "zero compiled");
}

TEST(VroomGtGlv4, SplitMatchesScalarBlstSplit)
{
    for (const uint64_t scalar_seed : std::array<uint64_t, 4>{{123, 321, 777, 999}}) {
        const auto scalar_limbs = blst_vroom_test::seeded_scalar_limbs_le(scalar_seed);
        int64_t scalar_split[4] = {};
        int64_t vroom_split[4] = {};

        blst::Fp12::gt_split_scalar_4d(scalar_split, scalar_limbs.data());
        vroom_gt_glv4_detail::gt_split_scalar_4d(vroom_split, scalar_limbs.data());

        for (size_t i = 0; i < 4; ++i) {
            EXPECT_EQ(vroom_split[i], scalar_split[i]) << "seed=" << scalar_seed << " idx=" << i;
        }
    }
}

}  // namespace
