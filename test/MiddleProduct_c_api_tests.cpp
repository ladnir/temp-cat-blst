#include <gtest/gtest.h>

#include <limits>
#include <type_traits>
#include <vector>

#include "blst/config.h"
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#define BLST_AVX_USE_CPP_TYPES 1
#endif
#include "bindings/blst_middle_product.h"
#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fft/context.h"
#include "blst/fft/middle_product.h"
#include "blst/fields/Fr.h"
#include "blst/fields/Fp12.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/point.h"
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/map_to_G.h"
#endif

namespace {

using BoundaryG1 = blst::scalar::G1Point;
using BoundaryG2 = blst::scalar::G2Point;
using BoundaryGt = blst::Fp12;
using PublicG1 = blst::G1Point;
using PublicG2 = blst::G2Point;
using PublicGt = blst::GT;

BoundaryG1 g1(uint64_t scalar)
{
    BoundaryG1 out = {};
    blst::Fr s = {};
    blst::Fr::from_uint64(s, scalar);
    blst::Fr s_non_mont = {};
    blst::Fr::from_mont(s_non_mont, s);
    BoundaryG1::mult(
        out,
        BoundaryG1::generator(),
        reinterpret_cast<const blst::byte*>(&s_non_mont),
        256);
    return out;
}

BoundaryG2 g2(uint64_t scalar)
{
    BoundaryG2 out = {};
    blst::Fr s = {};
    blst::Fr::from_uint64(s, scalar);
    blst::Fr s_non_mont = {};
    blst::Fr::from_mont(s_non_mont, s);
    BoundaryG2::mult(
        out,
        BoundaryG2::generator(),
        reinterpret_cast<const blst::byte*>(&s_non_mont),
        256);
    return out;
}

BoundaryGt gt_to_fp12(const PublicGt& value)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    {
        unsigned char serialized[12 * 48] = {};
        PublicGt::serialize(serialized, value);
        blst::scalar::GT scalar_gt = {};
        EXPECT_TRUE(blst::scalar::GT::deserialize_unchecked(scalar_gt, serialized));
        return scalar_gt.data;
    }
#else
    if constexpr (std::is_same_v<PublicGt, blst::scalar::GT>) {
        return value.data;
    } else {
        unsigned char serialized[12 * 48] = {};
        PublicGt::serialize(serialized, value);
        blst::scalar::GT scalar_gt = {};
        EXPECT_TRUE(blst::scalar::GT::deserialize_unchecked(scalar_gt, serialized));
        return scalar_gt.data;
    }
#endif
}

std::vector<PublicG1> to_public_g1(const std::vector<BoundaryG1>& input)
{
    auto convert = []<class PublicPoint, class BoundaryPoint>(const BoundaryPoint& point) -> PublicPoint {
        if constexpr (std::is_same_v<PublicPoint, BoundaryPoint>) {
            return point;
        } else {
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
            return blst::avx::avx_map_to_G_detail::from_scalar(point);
#else
            static_assert(std::is_same_v<PublicPoint, BoundaryPoint>, "Scalar builds should not require AVX point conversion.");
            return point;
#endif
        }
    };

    std::vector<PublicG1> out;
    out.reserve(input.size());
    for (const BoundaryG1& point : input) {
        out.push_back(convert.template operator()<PublicG1>(point));
    }
    return out;
}

std::vector<PublicG2> to_public_g2(const std::vector<BoundaryG2>& input)
{
    auto convert = []<class PublicPoint, class BoundaryPoint>(const BoundaryPoint& point) -> PublicPoint {
        if constexpr (std::is_same_v<PublicPoint, BoundaryPoint>) {
            return point;
        } else {
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
            return blst::avx::avx_map_to_G_detail::from_scalar(point);
#else
            static_assert(std::is_same_v<PublicPoint, BoundaryPoint>, "Scalar builds should not require AVX point conversion.");
            return point;
#endif
        }
    };

    std::vector<PublicG2> out;
    out.reserve(input.size());
    for (const BoundaryG2& point : input) {
        out.push_back(convert.template operator()<PublicG2>(point));
    }
    return out;
}

std::vector<BoundaryGt> reduced_middle_product(
    const std::vector<BoundaryG1>& a,
    const std::vector<BoundaryG2>& b)
{
    const std::vector<PublicG1> public_a = to_public_g1(a);
    const std::vector<PublicG2> public_b = to_public_g2(b);

    const auto precomputed = blst::precompute_middle_product_b_scalar(
        std::span<const PublicG2>(public_b),
        a.size());
    std::vector<PublicGt> raw(a.size(), PublicGt::identity());
    blst::middle_product_scalar(
        std::span<const PublicG1>(public_a),
        precomputed,
        std::span<PublicGt>(raw));

    std::vector<BoundaryGt> out;
    out.reserve(raw.size());
    for (const PublicGt& value : raw) {
        out.push_back(gt_to_fp12(value));
    }
    return out;
}

bool run_scalar_ranges_serial(
    const BoundaryG1* input,
    const std::size_t n,
    const blst_middle_product_precomputed_b_scalar* precomputed,
    blst_middle_product_workspace_scalar* workspace,
    BoundaryGt* output)
{
    const std::size_t fft_size = blst::middle_product_fft_size(n);
    const std::size_t stages = fft_log2_size(fft_size);
    if (!blst_middle_product_scalar_input_range(
            input, n, precomputed, workspace, 0, fft_size)) {
        return false;
    }
    const auto execute_full = [&](
        const blst_middle_product_range_phase phase,
        const std::size_t stage = 0) {
        std::size_t count = 0;
        return blst_middle_product_scalar_range_work_count(
                   precomputed, phase, stage, &count) &&
            blst_middle_product_scalar_execute_range(
                   precomputed, workspace, phase, stage, 0, count);
    };
    if (!execute_full(BLST_MIDDLE_PRODUCT_RANGE_FORWARD_BIT_REVERSE)) {
        return false;
    }
    for (std::size_t stage = 0; stage < stages; ++stage) {
        if (!execute_full(BLST_MIDDLE_PRODUCT_RANGE_FORWARD_J0, stage) ||
            !execute_full(BLST_MIDDLE_PRODUCT_RANGE_FORWARD_NONZERO, stage)) {
            return false;
        }
    }
    if (!execute_full(BLST_MIDDLE_PRODUCT_RANGE_AFFINE_CONVERSION) ||
        !execute_full(BLST_MIDDLE_PRODUCT_RANGE_POINTWISE_PAIRING) ||
        !execute_full(BLST_MIDDLE_PRODUCT_RANGE_INVERSE_BIT_REVERSE)) {
        return false;
    }
    for (std::size_t stage = 0; stage < stages; ++stage) {
        if (!execute_full(
                BLST_MIDDLE_PRODUCT_RANGE_INVERSE_BUTTERFLIES, stage)) {
            return false;
        }
    }
    return execute_full(BLST_MIDDLE_PRODUCT_RANGE_NORMALIZATION) &&
        blst_middle_product_scalar_output_range(
            precomputed, workspace, output, n, 0, n);
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
bool run_avx_native_ranges_serial(
    const blst::avx::G1Point* input,
    const std::size_t n,
    const blst_middle_product_precomputed_b_avx* precomputed,
    blst_middle_product_workspace_avx_native* workspace,
    blst::avx::GT* output)
{
    const std::size_t fft_size = blst::middle_product_fft_size(n);
    const std::size_t stages = fft_log2_size(fft_size);
    if (!blst_middle_product_avx_native_input_range(
            input, n, precomputed, workspace, 0, fft_size)) {
        return false;
    }
    const auto execute_full = [&](
        const blst_middle_product_range_phase phase,
        const std::size_t stage = 0) {
        std::size_t count = 0;
        return blst_middle_product_avx_native_range_work_count(
                   precomputed, phase, stage, &count) &&
            blst_middle_product_avx_native_execute_range(
                   precomputed, workspace, phase, stage, 0, count);
    };
    if (!execute_full(BLST_MIDDLE_PRODUCT_RANGE_FORWARD_BIT_REVERSE)) {
        return false;
    }
    for (std::size_t stage = 0; stage < stages; ++stage) {
        if (!execute_full(BLST_MIDDLE_PRODUCT_RANGE_FORWARD_J0, stage) ||
            !execute_full(BLST_MIDDLE_PRODUCT_RANGE_FORWARD_NONZERO, stage)) {
            return false;
        }
    }
    if (!execute_full(BLST_MIDDLE_PRODUCT_RANGE_AFFINE_CONVERSION) ||
        !execute_full(BLST_MIDDLE_PRODUCT_RANGE_POINTWISE_PAIRING) ||
        !execute_full(BLST_MIDDLE_PRODUCT_RANGE_INVERSE_BIT_REVERSE)) {
        return false;
    }
    for (std::size_t stage = 0; stage < stages; ++stage) {
        if (!execute_full(
                BLST_MIDDLE_PRODUCT_RANGE_INVERSE_BUTTERFLIES, stage)) {
            return false;
        }
    }
    return execute_full(BLST_MIDDLE_PRODUCT_RANGE_NORMALIZATION) &&
        blst_middle_product_avx_native_output_range(
            precomputed, workspace, output, n, 0, n);
}
#endif

} // namespace

TEST(MiddleProductCApi, ComputesExpectedMiddleWindow)
{
    const std::vector<BoundaryG1> a = {g1(3), g1(5), g1(7)};
    const std::vector<BoundaryG2> b = {g2(11), g2(13), g2(17), g2(19), g2(23), g2(29)};
    const std::vector<BoundaryGt> expected = reduced_middle_product(a, b);

    blst_middle_product_precomputed_b* precomputed =
        blst_middle_product_precompute_b(b.data(), a.size());
    ASSERT_NE(precomputed, nullptr);

    std::vector<BoundaryGt> out(a.size());
    ASSERT_TRUE(blst_middle_product(a.data(), a.size(), precomputed, out.data()));

    for (std::size_t k = 0; k < a.size(); ++k) {
        EXPECT_TRUE(blst::Fp12::equal(out[k], expected[k]))
            << "mismatch at output index " << k;
    }

    blst_middle_product_precomputed_b_free(precomputed);
}

TEST(MiddleProductCApi, RejectsLengthMismatch)
{
    const std::vector<BoundaryG1> a = {g1(3), g1(5)};
    const std::vector<BoundaryG2> b = {g2(11), g2(13), g2(17), g2(19)};

    blst_middle_product_precomputed_b* precomputed =
        blst_middle_product_precompute_b(b.data(), 2);
    ASSERT_NE(precomputed, nullptr);

    std::vector<BoundaryGt> out(2, BoundaryGt::one());
    const std::vector<BoundaryGt> before = out;
    EXPECT_FALSE(blst_middle_product(a.data(), 1, precomputed, out.data()));
    for (std::size_t i = 0; i < out.size(); ++i) {
        EXPECT_TRUE(blst::Fp12::equal(out[i], before[i]));
    }

    blst_middle_product_precomputed_b_free(precomputed);
}

TEST(MiddleProductCApi, RejectsUnsupportedLengthBeforeReadingInput)
{
    const BoundaryG2 dummy = BoundaryG2::infinity();

    EXPECT_EQ(
        blst_middle_product_precompute_b(
            &dummy,
            blst::middle_product_max_n + 1),
        nullptr);
    EXPECT_EQ(
        blst_middle_product_precompute_b(
            &dummy,
            std::numeric_limits<std::size_t>::max()),
        nullptr);
}

TEST(MiddleProductCApi, UpdateAndCombineReconstructFullResult)
{
    const std::vector<BoundaryG1> a_full = {g1(3), g1(5), g1(7), g1(11)};
    const std::vector<BoundaryG1> a_without_constant = {
        BoundaryG1::infinity(),
        g1(5),
        g1(7),
        g1(11),
    };
    const std::vector<BoundaryG2> b = {
        g2(13), g2(17), g2(19), g2(23),
        g2(29), g2(31), g2(37), g2(41),
    };
    const std::vector<BoundaryGt> scalar_expected =
        reduced_middle_product(a_full, b);
    const std::vector<BoundaryGt> scalar_mid =
        reduced_middle_product(a_without_constant, b);

    blst_middle_product_precomputed_b* precomputed =
        blst_middle_product_precompute_b(b.data(), a_full.size());
    ASSERT_NE(precomputed, nullptr);

    std::vector<BoundaryGt> expected(a_full.size());
    std::vector<BoundaryGt> mid(a_full.size());
    std::vector<BoundaryGt> cnst(a_full.size());
    std::vector<BoundaryGt> combined(a_full.size());
    std::vector<BoundaryGt> scalar_mid_plus_native_update(a_full.size());
    std::vector<BoundaryGt> exact_alias_cnst;

    ASSERT_TRUE(blst_middle_product(a_full.data(), a_full.size(), precomputed, expected.data()));
    ASSERT_TRUE(blst_middle_product(
        a_without_constant.data(),
        a_without_constant.size(),
        precomputed,
        mid.data()));
    ASSERT_TRUE(blst_middle_product_update(&a_full[0], precomputed, cnst.data()));
    ASSERT_TRUE(blst_middle_product_combine(mid.data(), cnst.data(), combined.size(), combined.data()));
    ASSERT_TRUE(blst_middle_product_combine(
        scalar_mid.data(),
        cnst.data(),
        scalar_mid_plus_native_update.size(),
        scalar_mid_plus_native_update.data()));

    std::vector<BoundaryGt> exact_alias = mid;
    ASSERT_TRUE(blst_middle_product_combine(
        exact_alias.data(),
        cnst.data(),
        exact_alias.size(),
        exact_alias.data()));
    exact_alias_cnst = cnst;
    ASSERT_TRUE(blst_middle_product_combine(
        mid.data(),
        exact_alias_cnst.data(),
        exact_alias_cnst.size(),
        exact_alias_cnst.data()));

    for (std::size_t i = 0; i < expected.size(); ++i) {
        EXPECT_TRUE(blst::Fp12::equal(expected[i], scalar_expected[i]))
            << "native full middle product differs from scalar at output index " << i;
        EXPECT_TRUE(blst::Fp12::equal(mid[i], scalar_mid[i]))
            << "native zero-constant middle product differs from scalar at output index " << i;
        EXPECT_TRUE(blst::Fp12::equal(scalar_expected[i], scalar_mid_plus_native_update[i]))
            << "native update differs from scalar contribution at output index " << i;
        EXPECT_TRUE(blst::Fp12::equal(expected[i], combined[i]))
            << "update/combine mismatch at output index " << i;
        EXPECT_TRUE(blst::Fp12::equal(expected[i], exact_alias[i]))
            << "exact-alias combine mismatch at output index " << i;
        EXPECT_TRUE(blst::Fp12::equal(expected[i], exact_alias_cnst[i]))
            << "second-input exact-alias combine mismatch at output index " << i;
    }

    blst_middle_product_precomputed_b_free(precomputed);
}

TEST(MiddleProductCApi, ScalarRangeBoundaryMatchesOneShotAndRejectsBadJobs)
{
    constexpr std::size_t n = 17;
    std::vector<BoundaryG1> a;
    std::vector<BoundaryG2> b;
    a.reserve(n);
    b.reserve(2 * n);
    for (std::size_t i = 0; i < n; ++i) {
        a.push_back(g1(3 + i));
    }
    for (std::size_t i = 0; i < 2 * n; ++i) {
        b.push_back(g2(11 + i));
    }

    blst_middle_product_precomputed_b_scalar* precomputed =
        blst_middle_product_precompute_b_scalar(b.data(), n);
    ASSERT_NE(precomputed, nullptr);
    blst_middle_product_workspace_scalar* workspace =
        blst_middle_product_workspace_scalar_new(precomputed);
    ASSERT_NE(workspace, nullptr);
    std::vector<BoundaryGt> expected(n);
    std::vector<BoundaryGt> actual(n);
    ASSERT_TRUE(blst_middle_product_scalar(
        a.data(), n, precomputed, expected.data()));
    ASSERT_TRUE(run_scalar_ranges_serial(
        a.data(), n, precomputed, workspace, actual.data()));
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(blst::Fp12::equal(expected[i], actual[i])) << i;
    }

    std::size_t count = 0xfeed;
    EXPECT_FALSE(blst_middle_product_scalar_range_work_count(
        precomputed,
        static_cast<blst_middle_product_range_phase>(255),
        0,
        &count));
    EXPECT_EQ(count, 0xfeed);
    EXPECT_FALSE(blst_middle_product_scalar_range_work_count(
        precomputed,
        BLST_MIDDLE_PRODUCT_RANGE_FORWARD_BIT_REVERSE,
        1,
        &count));
    ASSERT_TRUE(blst_middle_product_scalar_range_work_count(
        precomputed,
        BLST_MIDDLE_PRODUCT_RANGE_FORWARD_BIT_REVERSE,
        0,
        &count));
    EXPECT_FALSE(blst_middle_product_scalar_execute_range(
        precomputed,
        workspace,
        BLST_MIDDLE_PRODUCT_RANGE_FORWARD_BIT_REVERSE,
        0,
        count,
        count + 1));
    EXPECT_FALSE(blst_middle_product_scalar_execute_range(
        precomputed,
        workspace,
        BLST_MIDDLE_PRODUCT_RANGE_INPUT_LOAD,
        0,
        0,
        0));

    const std::vector<BoundaryGt> canary(n, BoundaryGt::one());
    actual = canary;
    EXPECT_FALSE(blst_middle_product_scalar_output_range(
        precomputed, workspace, actual.data(), n, n, n + 1));
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(blst::Fp12::equal(actual[i], canary[i])) << i;
    }

    blst_middle_product_workspace_scalar_free(workspace);
    blst_middle_product_precomputed_b_scalar_free(precomputed);
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
TEST(MiddleProductCApi, NativeAvxRangeBoundaryMatchesOneShot)
{
    constexpr std::size_t n = 17;
    std::vector<blst::avx::G1Point> a;
    std::vector<blst::avx::G2Point> b;
    a.reserve(n);
    b.reserve(2 * n);
    for (std::size_t i = 0; i < n; ++i) {
        a.push_back(blst::avx::avx_map_to_G_detail::from_scalar(g1(3 + i)));
    }
    for (std::size_t i = 0; i < 2 * n; ++i) {
        b.push_back(blst::avx::avx_map_to_G_detail::from_scalar(g2(11 + i)));
    }

    blst_middle_product_precomputed_b_avx* precomputed =
        blst_middle_product_precompute_b_avx_native(b.data(), n);
    ASSERT_NE(precomputed, nullptr);
    blst_middle_product_workspace_avx_native* workspace =
        blst_middle_product_workspace_avx_native_new(precomputed);
    ASSERT_NE(workspace, nullptr);
    std::vector<blst::avx::GT> expected(n);
    std::vector<blst::avx::GT> actual(n);
    ASSERT_TRUE(blst_middle_product_avx_native_with_workspace(
        a.data(), n, precomputed, workspace, expected.data()));
    ASSERT_TRUE(run_avx_native_ranges_serial(
        a.data(), n, precomputed, workspace, actual.data()));
    for (std::size_t i = 0; i < n; ++i) {
        EXPECT_TRUE(blst::avx::GT::equal(expected[i], actual[i])) << i;
    }

    blst_middle_product_workspace_avx_native_free(workspace);
    blst_middle_product_precomputed_b_avx_free(precomputed);
}
#endif
