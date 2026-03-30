#include <gtest/gtest.h>

#include <type_traits>
#include <vector>

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

struct blst_middle_product_precomputed_b;

extern "C" {
blst_middle_product_precomputed_b* blst_middle_product_precompute_b(
    const blst::scalar::G2Point* b,
    size_t n);
bool blst_middle_product(
    const blst::scalar::G1Point* a,
    size_t n,
    const blst_middle_product_precomputed_b* precomputed_b,
    blst::Fp12* out);
void blst_middle_product_precomputed_b_free(
    blst_middle_product_precomputed_b* precomputed_b);
}

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
        blst::scalar::GT::deserialize(scalar_gt, serialized);
        return scalar_gt.data;
    }
#else
    if constexpr (std::is_same_v<PublicGt, blst::scalar::GT>) {
        return value.data;
    } else {
        unsigned char serialized[12 * 48] = {};
        PublicGt::serialize(serialized, value);
        blst::scalar::GT scalar_gt = {};
        blst::scalar::GT::deserialize(scalar_gt, serialized);
        return scalar_gt.data;
    }
#endif
}

std::vector<PublicG1> to_public_g1(const std::vector<BoundaryG1>& input)
{
    std::vector<PublicG1> out;
    out.reserve(input.size());
    for (const BoundaryG1& point : input) {
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
        out.push_back(blst::avx::avx_map_to_G_detail::from_scalar(point));
#else
        out.push_back(point);
#endif
    }
    return out;
}

std::vector<PublicG2> to_public_g2(const std::vector<BoundaryG2>& input)
{
    std::vector<PublicG2> out;
    out.reserve(input.size());
    for (const BoundaryG2& point : input) {
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
        out.push_back(blst::avx::avx_map_to_G_detail::from_scalar(point));
#else
        out.push_back(point);
#endif
    }
    return out;
}

std::vector<BoundaryGt> reduced_middle_product(
    const std::vector<BoundaryG1>& a,
    const std::vector<BoundaryG2>& b)
{
    const std::vector<PublicG1> public_a = to_public_g1(a);
    const std::vector<PublicG2> public_b = to_public_g2(b);

    const auto precomputed = blst::precompute_middle_product_b(
        std::span<const PublicG2>(public_b),
        a.size());
    std::vector<PublicGt> raw(a.size(), PublicGt::identity());
    blst::middle_product(
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

    std::vector<BoundaryGt> out(2);
    EXPECT_FALSE(blst_middle_product(a.data(), 1, precomputed, out.data()));

    blst_middle_product_precomputed_b_free(precomputed);
}
