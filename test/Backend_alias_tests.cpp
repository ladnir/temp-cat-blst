#include <gtest/gtest.h>

#include <span>
#include <type_traits>
#include <vector>

#include "blst/curve/GT.h"
#include "blst/curve/map_to_G.h"
#include "blst/curve/point.h"
#include "blst/scalar/curve/GT.h"
#include "blst/fft/middle_product.h"
#include "blst/fields/Fr.h"

using namespace blst;

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
using AvxScopedRoots = blst::avx::AvxGtCoeffCtx::ScopedRoots;
using AvxPreparedRoots = blst::avx::AvxGtPreparedRoots;

static_assert(std::is_constructible_v<AvxScopedRoots, const AvxPreparedRoots&>);
static_assert(std::is_nothrow_constructible_v<
    AvxScopedRoots,
    const AvxPreparedRoots&>);
static_assert(std::is_nothrow_destructible_v<AvxScopedRoots>);
static_assert(!std::is_constructible_v<AvxScopedRoots, AvxPreparedRoots&&>);
static_assert(!std::is_constructible_v<AvxScopedRoots, const AvxPreparedRoots&&>);
static_assert(!std::is_copy_constructible_v<AvxScopedRoots>);
static_assert(!std::is_copy_assignable_v<AvxScopedRoots>);
static_assert(!std::is_move_constructible_v<AvxScopedRoots>);
static_assert(!std::is_move_assignable_v<AvxScopedRoots>);
#endif

static G1Point point_from_scalar_g1(u64 scalar)
{
    G1Point out = {};
    Fr scalar_mont = {};
    Fr::from_uint64(scalar_mont, scalar);
    Fr scalar_non_mont = {};
    Fr::from_mont(scalar_non_mont, scalar_mont);
    G1Point::mult(out, G1Point::generator(), reinterpret_cast<const byte*>(&scalar_non_mont), 256);
    return out;
}

static G2Point point_from_scalar_g2(u64 scalar)
{
    G2Point out = {};
    Fr scalar_mont = {};
    Fr::from_uint64(scalar_mont, scalar);
    Fr scalar_non_mont = {};
    Fr::from_mont(scalar_non_mont, scalar_mont);
    G2Point::mult(out, G2Point::generator(), reinterpret_cast<const byte*>(&scalar_non_mont), 256);
    return out;
}

TEST(BackendAlias, TopLevelAliasesMatchConfig)
{
#if BLST_DEFAULT_AVX
    EXPECT_TRUE((std::is_same_v<G1Point, avx::G1Point>));
    EXPECT_TRUE((std::is_same_v<G1AffinePoint, avx::G1AffinePoint>));
    EXPECT_TRUE((std::is_same_v<G2Point, avx::G2Point>));
    EXPECT_TRUE((std::is_same_v<G2AffinePoint, avx::G2AffinePoint>));
    EXPECT_TRUE((std::is_same_v<GT, avx::GT>));
#else
    EXPECT_TRUE((std::is_same_v<G1Point, scalar::G1Point>));
    EXPECT_TRUE((std::is_same_v<G1AffinePoint, scalar::G1AffinePoint>));
    EXPECT_TRUE((std::is_same_v<G2Point, scalar::G2Point>));
    EXPECT_TRUE((std::is_same_v<G2AffinePoint, scalar::G2AffinePoint>));
    EXPECT_TRUE((std::is_same_v<GT, scalar::GT>));
#endif
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
TEST(BackendAlias, ScalarInfinityBridgeProducesMultiplicationSafeAvxIdentity)
{
    const avx::G1Point g1 = avx::avx_map_to_G_detail::from_scalar(
        scalar::G1Point::infinity());
    const avx::G2Point g2 = avx::avx_map_to_G_detail::from_scalar(
        scalar::G2Point::infinity());
    const byte scalar_bytes[1] = {7};

    avx::G1Point g1_product = {};
    avx::G1Point::mult_vartime(g1_product, g1, scalar_bytes, 8);
    EXPECT_TRUE(avx::G1Point::is_inf(g1));
    EXPECT_TRUE(avx::G1Point::is_inf(g1_product));

    avx::G2Point g2_product = {};
    avx::G2Point::mult_vartime(g2_product, g2, scalar_bytes, 8);
    EXPECT_TRUE(avx::G2Point::is_inf(g2));
    EXPECT_TRUE(avx::G2Point::is_inf(g2_product));
}
#endif

TEST(BackendAlias, TopLevelHashToCurveCompilesAndProducesSubgroupPoints)
{
    static constexpr unsigned char msg[] = "top-level backend alias";
    static constexpr unsigned char dst_g1[] = "BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_TEST";
    static constexpr unsigned char dst_g2[] = "BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_TEST";

    G1Point g1 = {};
    Hash_to_curve(g1, msg, sizeof(msg) - 1, dst_g1, sizeof(dst_g1) - 1, nullptr, 0);
    EXPECT_TRUE(G1Point::on_curve(g1));
    EXPECT_TRUE(G1Point::in_subgroup(g1));

    G2Point g2 = {};
    Hash_to_curve(g2, msg, sizeof(msg) - 1, dst_g2, sizeof(dst_g2) - 1, nullptr, 0);
    EXPECT_TRUE(G2Point::on_curve(g2));
    EXPECT_TRUE(G2Point::in_subgroup(g2));
}

TEST(BackendAlias, TopLevelHashToCurveRoundTripsThroughAffine)
{
    static constexpr unsigned char msg_g1[] = "public-api-g1";
    static constexpr unsigned char dst_g1[] = "BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_TEST";

    G1Point g1 = {};
    Hash_to_curve(g1, msg_g1, sizeof(msg_g1) - 1, dst_g1, sizeof(dst_g1) - 1, nullptr, 0);

    G1AffinePoint g1_affine = {};
    G1Point::to_affine(g1_affine, g1);
    EXPECT_TRUE(G1AffinePoint::in_subgroup(g1_affine));

    unsigned char g1_compressed[48] = {};
    G1Point::compress(g1_compressed, g1);

    G1AffinePoint g1_decoded = {};
    ASSERT_EQ(BLST_SUCCESS, G1AffinePoint::uncompress(g1_decoded, g1_compressed));

    G1Point g1_roundtrip = {};
    G1Point::from_affine(g1_roundtrip, g1_decoded);
    EXPECT_TRUE(G1Point::equal(g1, g1_roundtrip));

    static constexpr unsigned char msg_g2[] = "public-api-g2";
    static constexpr unsigned char dst_g2[] = "BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_TEST";

    G2Point g2 = {};
    Hash_to_curve(g2, msg_g2, sizeof(msg_g2) - 1, dst_g2, sizeof(dst_g2) - 1, nullptr, 0);

    G2AffinePoint g2_affine = {};
    G2Point::to_affine(g2_affine, g2);
    EXPECT_TRUE(G2AffinePoint::in_subgroup(g2_affine));

    unsigned char g2_compressed[96] = {};
    G2Point::compress(g2_compressed, g2);

    G2AffinePoint g2_decoded = {};
    ASSERT_EQ(BLST_SUCCESS, G2AffinePoint::uncompress(g2_decoded, g2_compressed));

    G2Point g2_roundtrip = {};
    G2Point::from_affine(g2_roundtrip, g2_decoded);
    EXPECT_TRUE(G2Point::equal(g2, g2_roundtrip));
}

TEST(BackendAlias, TopLevelPairingSurfaceMatchesBilinearDoubling)
{
    const G1Point g1 = G1Point::generator();
    const G2Point g2 = G2Point::generator();

    G1Point g1_double = {};
    G1Point::doubling(g1_double, g1);

    G2Point g2_double = {};
    G2Point::doubling(g2_double, g2);

    G1AffinePoint g1_affine = {};
    G1AffinePoint g1_double_affine = {};
    G1Point::to_affine(g1_affine, g1);
    G1Point::to_affine(g1_double_affine, g1_double);

    G2AffinePoint g2_affine = {};
    G2AffinePoint g2_double_affine = {};
    G2Point::to_affine(g2_affine, g2);
    G2Point::to_affine(g2_double_affine, g2_double);

    MillerResult base_miller = {};
    GT base = {};
    GT::miller_loop(base_miller, &g2_affine, &g1_affine);
    GT::final_exp(base, base_miller);

    MillerResult lhs_miller = {};
    GT lhs = {};
    GT::miller_loop(lhs_miller, &g2_affine, &g1_double_affine);
    GT::final_exp(lhs, lhs_miller);

    MillerResult rhs_miller = {};
    GT rhs = {};
    GT::miller_loop(rhs_miller, &g2_double_affine, &g1_affine);
    GT::final_exp(rhs, rhs_miller);

    GT doubled = {};
    GT::doubling(doubled, base);

    EXPECT_TRUE(GT::equal(lhs, rhs));
    EXPECT_TRUE(GT::equal(lhs, doubled));
}

TEST(BackendAlias, ScalarMillerLoopNEmptyReturnsIdentityAndMalformedInputFailsClosed)
{
    blst::scalar::MillerResult scalar_multi = {};
    const blst::scalar::G2AffinePoint* const* scalar_empty_qs = nullptr;
    const blst::scalar::G1AffinePoint* const* scalar_empty_ps = nullptr;
    blst::scalar::GT::miller_loop_n(scalar_multi, scalar_empty_qs, scalar_empty_ps, 0);
    blst::scalar::GT scalar_final = {};
    blst::scalar::GT::final_exp(scalar_final, scalar_multi);
    EXPECT_TRUE(blst::scalar::GT::is_identity(scalar_final));

    const std::array<const blst::scalar::G1AffinePoint*, 1> scalar_ps = {nullptr};
    const std::array<const blst::scalar::G2AffinePoint*, 1> scalar_qs = {nullptr};
    blst::scalar::GT::miller_loop_n(
        scalar_multi,
        scalar_qs.data(),
        scalar_ps.data(),
        scalar_ps.size());
    blst::scalar::GT::final_exp(scalar_final, scalar_multi);
    EXPECT_TRUE(blst::scalar::GT::equal(scalar_final, {}));
    EXPECT_FALSE(blst::scalar::GT::is_identity(scalar_final));
}

#if BLST_DEFAULT_AVX
TEST(BackendAlias, AvxMillerLoopNEmptyReturnsIdentityAndMalformedInputFailsClosed)
{
    blst::avx::MillerResult avx_multi = {};
    const blst::avx::G2AffinePoint* const* avx_empty_qs = nullptr;
    const blst::avx::G1AffinePoint* const* avx_empty_ps = nullptr;
    blst::avx::GT::miller_loop_n(avx_multi, avx_empty_qs, avx_empty_ps, 0);
    blst::avx::GT avx_final = {};
    blst::avx::GT::final_exp(avx_final, avx_multi);
    EXPECT_TRUE(blst::avx::GT::is_identity(avx_final));

    const std::array<const blst::avx::G1AffinePoint*, 1> avx_ps = {nullptr};
    const std::array<const blst::avx::G2AffinePoint*, 1> avx_qs = {nullptr};
    blst::avx::GT::miller_loop_n(
        avx_multi,
        avx_qs.data(),
        avx_ps.data(),
        avx_ps.size());
    blst::avx::GT::final_exp(avx_final, avx_multi);
    EXPECT_TRUE(blst::avx::GT::equal(avx_final, {}));
    EXPECT_FALSE(blst::avx::GT::is_identity(avx_final));
}
#endif

TEST(BackendAlias, MiddleProductReducedMatchesBasicOnSmallVector)
{
    using ScalarCtx = blst::detail::ScalarMiddleProductContext;

    const std::vector<u64> a_scalars = { 3, 5, 7 };
    const std::vector<u64> b_scalars = { 11, 13, 17, 19, 23, 29 };

    std::vector<G1Point> a_points;
    std::vector<G2Point> b_points;
    a_points.reserve(a_scalars.size());
    b_points.reserve(b_scalars.size());

    for (u64 scalar : a_scalars) {
        a_points.push_back(point_from_scalar_g1(scalar));
    }
    for (u64 scalar : b_scalars) {
        b_points.push_back(point_from_scalar_g2(scalar));
    }

    std::vector<GT> basic(a_points.size(), GT::identity());
    std::vector<GT> reduced_scalar(a_points.size(), GT::identity());
    std::vector<GT> reduced_raw(a_points.size(), GT::identity());
    std::vector<GT> reduced_pre(a_points.size(), GT::identity());

    middle_product_basic(
        std::span<const G1Point>(a_points),
        std::span<const G2Point>(b_points),
        std::span<GT>(basic));

    middle_product(
        std::span<const G1Point>(a_points),
        std::span<const G2Point>(b_points),
        std::span<GT>(reduced_raw));

    middle_product_scalar(
        std::span<const G1Point>(a_points),
        std::span<const G2Point>(b_points),
        std::span<GT>(reduced_scalar));

    const auto precomputed_b = precompute_middle_product_b(
        std::span<const G2Point>(b_points),
        a_points.size());
    middle_product(
        std::span<const G1Point>(a_points),
        precomputed_b,
        std::span<GT>(reduced_pre));

    for (std::size_t i = 0; i < a_points.size(); ++i) {
        EXPECT_TRUE(GT::equal(basic[i], reduced_scalar[i]))
            << "Reduced scalar mismatch vs basic at i=" << i;
        EXPECT_TRUE(GT::equal(basic[i], reduced_raw[i]))
            << "Reduced raw mismatch vs basic at i=" << i;
        EXPECT_TRUE(GT::equal(basic[i], reduced_pre[i]))
            << "Reduced precomputed mismatch vs basic at i=" << i;
    }

    const auto scalar_precomputed = blst::detail::precompute_middle_product_b_public<ScalarCtx>(
        std::span<const G2Point>(b_points),
        a_points.size(),
        nullptr);
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    using AvxCtx = blst::detail::AvxMiddleProductContext;
    const auto avx_precomputed = blst::detail::precompute_middle_product_b_public<AvxCtx>(
        std::span<const G2Point>(b_points),
        a_points.size(),
        nullptr);

    const auto scalar_a = blst::detail::convert_g1_inputs<ScalarCtx>(std::span<const G1Point>(a_points));
    const auto avx_a = blst::detail::convert_g1_inputs<AvxCtx>(std::span<const G1Point>(a_points));

    std::vector<ScalarCtx::G1Point> scalar_a_fft(
        scalar_precomputed.fft_size,
        ScalarCtx::g1_infinity());
    std::vector<AvxCtx::G1Point> avx_a_fft(
        avx_precomputed.fft_size,
        AvxCtx::g1_infinity());
    std::copy(scalar_a.begin(), scalar_a.end(), scalar_a_fft.begin());
    std::copy(avx_a.begin(), avx_a.end(), avx_a_fft.begin());

    const std::size_t xn_a = fft_admissible_size(a_points.size(), scalar_precomputed.fft_size);
    fft_forward_truncated_prefix<ScalarCtx::G1Point, Fr, ScalarCtx::G1CoeffCtx>(
        std::span<ScalarCtx::G1Point>(scalar_a_fft),
        std::span<const Fr>(scalar_precomputed.roots),
        xn_a,
        nullptr,
        "backend_alias.scalar_fft");
    fft_forward_truncated_prefix<AvxCtx::G1Point, Fr, AvxCtx::G1CoeffCtx>(
        std::span<AvxCtx::G1Point>(avx_a_fft),
        std::span<const Fr>(avx_precomputed.roots),
        xn_a,
        nullptr,
        "backend_alias.avx_fft");

    std::vector<ScalarCtx::GtFftElement> scalar_c_fft(
        scalar_precomputed.fft_size,
        ScalarCtx::identity());
    std::vector<AvxCtx::GtFftElement> avx_c_fft(
        avx_precomputed.fft_size,
        AvxCtx::identity());
    for (std::size_t i = 0; i < scalar_precomputed.fft_size; ++i) {
        scalar_c_fft[i] = ScalarCtx::pointwise_pair_direct(scalar_a_fft[i], scalar_precomputed.b_fft[i]);
        avx_c_fft[i] = AvxCtx::pointwise_pair_direct(avx_a_fft[i], avx_precomputed.b_fft[i]);
        EXPECT_TRUE(GT::equal(ScalarCtx::to_public_gt(scalar_c_fft[i]), AvxCtx::to_public_gt(avx_c_fft[i])))
            << "Scalar/AVX reduced pointwise mismatch at fft index " << i;
    }

    {
        const AvxCtx::GtFftElement sample = avx_c_fft[1];
        for (std::size_t i = 0; i < avx_precomputed.inv_roots.size(); ++i) {
            const Fr& root = avx_precomputed.inv_roots.scalar_span()[i];

            blst::avx::AvxGtElement via_prepared = {};
            {
                blst::avx::AvxGtCoeffCtx::ScopedRoots scoped_roots(avx_precomputed.inv_roots);
                via_prepared = blst::avx::AvxGtCoeffCtx::mult(sample, root);
            }

            Fr root_copy = root;
            const blst::avx::AvxGtElement via_fallback = blst::avx::AvxGtCoeffCtx::mult(sample, root_copy);
            EXPECT_TRUE(GT::equal(
                blst::avx::AvxGtElement::to_default(via_prepared),
                blst::avx::AvxGtElement::to_default(via_fallback)))
                << "Prepared root mismatch at inverse root index " << i;
        }

        const std::array<Fr, 1> outer_scalars = { Fr::from_uint64(2) };
        const std::array<Fr, 1> replacement_scalars = { Fr::from_uint64(3) };
        const std::array<Fr, 1> inner_scalars = { Fr::from_uint64(5) };
        blst::avx::AvxGtPreparedRoots outer_roots(outer_scalars);
        const blst::avx::AvxGtPreparedRoots replacement_roots(
            replacement_scalars);
        const blst::avx::AvxGtPreparedRoots inner_roots(inner_scalars);

        // Deliberately substitute a different compiled exponent so this test
        // can distinguish prepared-table selection from the scalar fallback.
        outer_roots.compiled[0] = replacement_roots.compiled[0];
        Fr outer_scalar_copy = outer_roots.scalars[0];
        const blst::avx::AvxGtElement fallback =
            blst::avx::AvxGtCoeffCtx::mult(sample, outer_scalar_copy);

        blst::avx::AvxGtElement selected = {};
        blst::avx::AvxGtElement restored_after_exception = {};
        {
            blst::avx::AvxGtCoeffCtx::ScopedRoots outer_scope(outer_roots);
            selected = blst::avx::AvxGtCoeffCtx::mult(
                sample,
                outer_roots.scalars[0]);
            try {
                blst::avx::AvxGtCoeffCtx::ScopedRoots inner_scope(inner_roots);
                throw 1;
            } catch (int) {
            }
            restored_after_exception = blst::avx::AvxGtCoeffCtx::mult(
                sample,
                outer_roots.scalars[0]);
        }
        const blst::avx::AvxGtElement restored_fallback =
            blst::avx::AvxGtCoeffCtx::mult(sample, outer_roots.scalars[0]);

        EXPECT_FALSE(GT::equal(
            blst::avx::AvxGtElement::to_default(selected),
            blst::avx::AvxGtElement::to_default(fallback)));
        EXPECT_TRUE(GT::equal(
            blst::avx::AvxGtElement::to_default(selected),
            blst::avx::AvxGtElement::to_default(restored_after_exception)));
        EXPECT_TRUE(GT::equal(
            blst::avx::AvxGtElement::to_default(fallback),
            blst::avx::AvxGtElement::to_default(restored_fallback)));
    }

    auto scalar_c_fft_full = scalar_c_fft;
    auto avx_c_fft_full_prepared = avx_c_fft;
    auto avx_c_fft_full_fallback = avx_c_fft;
    auto avx_c_fft_trunc_fallback = avx_c_fft;

    fft<ScalarCtx::GtFftElement, Fr, ScalarCtx::GtCoeffCtx>(
        std::span<ScalarCtx::GtFftElement>(scalar_c_fft_full),
        std::span<const Fr>(scalar_precomputed.inv_roots),
        true);

    {
        blst::avx::AvxGtCoeffCtx::ScopedRoots scoped_roots(avx_precomputed.inv_roots);
        fft<AvxCtx::GtFftElement, Fr, blst::avx::AvxGtCoeffCtx>(
            std::span<AvxCtx::GtFftElement>(avx_c_fft_full_prepared),
            avx_precomputed.inv_roots.scalar_span(),
            true);
    }

    {
        std::vector<Fr> inv_roots_copy(
            avx_precomputed.inv_roots.scalar_span().begin(),
            avx_precomputed.inv_roots.scalar_span().end());
        fft<AvxCtx::GtFftElement, Fr, blst::avx::AvxGtCoeffCtx>(
            std::span<AvxCtx::GtFftElement>(avx_c_fft_full_fallback),
            std::span<const Fr>(inv_roots_copy),
            true);
    }

    {
        auto scalar_stage = scalar_c_fft;
        auto avx_stage = avx_c_fft;
        std::vector<Fr> inv_roots_copy(
            avx_precomputed.inv_roots.scalar_span().begin(),
            avx_precomputed.inv_roots.scalar_span().end());

        bit_reverse(std::span<ScalarCtx::GtFftElement>(scalar_stage));
        bit_reverse(std::span<AvxCtx::GtFftElement>(avx_stage));

        for (std::size_t len = 2, stage = 0; len <= scalar_stage.size(); len <<= 1, ++stage) {
            const std::size_t half = len >> 1;
            for (std::size_t i = 0; i < scalar_stage.size(); i += len) {
                {
                    const auto u_scalar = scalar_stage[i];
                    const auto t_scalar = scalar_stage[i + half];
                    scalar_stage[i] = ScalarCtx::GtCoeffCtx::plus(u_scalar, t_scalar);
                    scalar_stage[i + half] = ScalarCtx::GtCoeffCtx::minus(u_scalar, t_scalar);

                    const auto u_avx = avx_stage[i];
                    const auto t_avx = avx_stage[i + half];
                    avx_stage[i] = blst::avx::AvxGtCoeffCtx::plus(u_avx, t_avx);
                    avx_stage[i + half] = blst::avx::AvxGtCoeffCtx::minus(u_avx, t_avx);
                }

                for (std::size_t j = 1; j < half; ++j) {
                    const Fr& w_scalar = scalar_precomputed.inv_roots[j * (scalar_stage.size() / len)];
                    const Fr& w_avx = inv_roots_copy[j * (avx_stage.size() / len)];

                    const auto u_scalar = scalar_stage[i + j];
                    const auto t_scalar =
                        ScalarCtx::GtCoeffCtx::mult(scalar_stage[i + j + half], w_scalar);
                    const auto u_avx = avx_stage[i + j];
                    const auto base_avx = avx_stage[i + j + half];
                    const auto t_avx =
                        blst::avx::AvxGtCoeffCtx::mult(base_avx, w_avx);
                    const std::array<uint64_t, 4> w_limbs = blst::avx::fr_to_limbs_le(w_avx);
                    const blst::avx::AvxGtElement t_avx_basic(
                        blst::avx::avx_cyclotomic_exp_basic(
                            base_avx.data,
                            w_limbs.data(),
                            blst::avx::AvxGtElement::fp12ring(),
                            blst::avx::AvxGtElement::ring()));

                    EXPECT_TRUE(GT::equal(
                        ScalarCtx::to_public_gt(u_scalar),
                        AvxCtx::to_public_gt(u_avx)))
                        << "Scalar/AVX stage input mismatch at stage=" << stage
                        << ", len=" << len
                        << ", block=" << i
                        << ", j=" << j;
                    EXPECT_TRUE(GT::equal(
                        ScalarCtx::to_public_gt(t_scalar),
                        AvxCtx::to_public_gt(t_avx)))
                        << "Scalar/AVX twiddle-product mismatch at stage=" << stage
                        << ", len=" << len
                        << ", block=" << i
                        << ", j=" << j;
                    EXPECT_TRUE(GT::equal(
                        AvxCtx::to_public_gt(t_avx),
                        AvxCtx::to_public_gt(t_avx_basic)))
                        << "AVX GLV/basic twiddle-product mismatch at stage=" << stage
                        << ", len=" << len
                        << ", block=" << i
                        << ", j=" << j;

                    scalar_stage[i + j] = ScalarCtx::GtCoeffCtx::plus(u_scalar, t_scalar);
                    scalar_stage[i + j + half] = ScalarCtx::GtCoeffCtx::minus(u_scalar, t_scalar);
                    avx_stage[i + j] = blst::avx::AvxGtCoeffCtx::plus(u_avx, t_avx);
                    avx_stage[i + j + half] = blst::avx::AvxGtCoeffCtx::minus(u_avx, t_avx);
                }
            }

            for (std::size_t i = 0; i < scalar_stage.size(); ++i) {
                EXPECT_TRUE(GT::equal(
                    ScalarCtx::to_public_gt(scalar_stage[i]),
                    AvxCtx::to_public_gt(avx_stage[i])))
                    << "Scalar/AVX stage mismatch at stage=" << stage
                    << ", len=" << len
                    << ", fft index=" << i;
            }
        }
    }

    ScalarCtx::inverse_fft_middle_window(
        std::span<ScalarCtx::GtFftElement>(scalar_c_fft),
        scalar_precomputed.inv_roots,
        a_points.size(),
        a_points.size());
    AvxCtx::inverse_fft_middle_window(
        std::span<AvxCtx::GtFftElement>(avx_c_fft),
        avx_precomputed.inv_roots,
        a_points.size(),
        a_points.size());
    {
        std::vector<Fr> inv_roots_copy(
            avx_precomputed.inv_roots.scalar_span().begin(),
            avx_precomputed.inv_roots.scalar_span().end());
        fft_inverse_truncated_range<AvxCtx::GtFftElement, Fr, blst::avx::AvxGtCoeffCtx>(
            std::span<AvxCtx::GtFftElement>(avx_c_fft_trunc_fallback),
            std::span<const Fr>(inv_roots_copy),
            a_points.size(),
            a_points.size());
    }

    for (std::size_t i = 0; i < scalar_c_fft_full.size(); ++i) {
        EXPECT_TRUE(GT::equal(
            ScalarCtx::to_public_gt(scalar_c_fft_full[i]),
            AvxCtx::to_public_gt(avx_c_fft_full_prepared[i])))
            << "Scalar/AVX full-iFFT prepared mismatch at fft index " << i;
        EXPECT_TRUE(GT::equal(
            AvxCtx::to_public_gt(avx_c_fft_full_prepared[i]),
            AvxCtx::to_public_gt(avx_c_fft_full_fallback[i])))
            << "AVX full-iFFT prepared/fallback mismatch at fft index " << i;
    }

    for (std::size_t i = 0; i < a_points.size(); ++i) {
        EXPECT_TRUE(GT::equal(
            ScalarCtx::to_public_gt(scalar_c_fft[a_points.size() + i]),
            AvxCtx::to_public_gt(avx_c_fft[a_points.size() + i])))
            << "Scalar/AVX reduced post-iFFT mismatch at output index " << i;
        EXPECT_TRUE(GT::equal(
            AvxCtx::to_public_gt(avx_c_fft[a_points.size() + i]),
            AvxCtx::to_public_gt(avx_c_fft_trunc_fallback[a_points.size() + i])))
            << "AVX truncated-iFFT prepared/fallback mismatch at output index " << i;
        EXPECT_TRUE(GT::equal(
            ScalarCtx::to_public_gt(scalar_c_fft[a_points.size() + i]),
            ScalarCtx::to_public_gt(scalar_c_fft_full[a_points.size() + i])))
            << "Scalar truncated/full mismatch at output index " << i;
        EXPECT_TRUE(GT::equal(
            AvxCtx::to_public_gt(avx_c_fft[a_points.size() + i]),
            AvxCtx::to_public_gt(avx_c_fft_full_prepared[a_points.size() + i])))
            << "AVX truncated/full prepared mismatch at output index " << i;
    }
#endif
}

TEST(BackendAlias, MiddleProductUpdateAndCombineReconstructFullResult)
{
    const std::vector<u64> a_scalars = { 3, 5, 7, 11 };
    const std::vector<u64> b_scalars = { 13, 17, 19, 23, 29, 31, 37, 41 };

    std::vector<G1Point> a_full;
    std::vector<G1Point> a_without_constant;
    std::vector<G2Point> b_points;
    a_full.reserve(a_scalars.size());
    a_without_constant.reserve(a_scalars.size());
    b_points.reserve(b_scalars.size());

    for (std::size_t i = 0; i < a_scalars.size(); ++i) {
        const G1Point point = point_from_scalar_g1(a_scalars[i]);
        a_full.push_back(point);
        a_without_constant.push_back((i == 0) ? G1Point::infinity() : point);
    }
    for (u64 scalar : b_scalars) {
        b_points.push_back(point_from_scalar_g2(scalar));
    }

    const auto precomputed_b = precompute_middle_product_b(
        std::span<const G2Point>(b_points),
        a_full.size());

    std::vector<GT> expected(a_full.size(), GT::identity());
    std::vector<GT> mid(a_full.size(), GT::identity());
    std::vector<GT> cnst(a_full.size(), GT::identity());
    std::vector<GT> combined(a_full.size(), GT::identity());

    middle_product(
        std::span<const G1Point>(a_full),
        precomputed_b,
        std::span<GT>(expected));
    middle_product(
        std::span<const G1Point>(a_without_constant),
        precomputed_b,
        std::span<GT>(mid));
    update(a_full[0], precomputed_b, std::span<GT>(cnst));
    combine(std::span<const GT>(mid), std::span<const GT>(cnst), std::span<GT>(combined));

    for (std::size_t i = 0; i < expected.size(); ++i) {
        EXPECT_TRUE(GT::equal(expected[i], combined[i]))
            << "Update/combine mismatch at i=" << i;
    }
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
TEST(BackendAlias, AvxGtInverseFftStageReplayMatchesScalar)
{
    using ScalarCtx = blst::detail::ScalarMiddleProductContext;
    using AvxCtx = blst::detail::AvxMiddleProductContext;

    const std::vector<u64> a_scalars = { 3, 5, 7 };
    const std::vector<u64> b_scalars = { 11, 13, 17, 19, 23, 29 };

    std::vector<G1Point> a_points;
    std::vector<G2Point> b_points;
    for (u64 scalar : a_scalars) {
        a_points.push_back(point_from_scalar_g1(scalar));
    }
    for (u64 scalar : b_scalars) {
        b_points.push_back(point_from_scalar_g2(scalar));
    }

    const auto scalar_precomputed = blst::detail::precompute_middle_product_b_public<ScalarCtx>(
        std::span<const G2Point>(b_points),
        a_points.size(),
        nullptr);
    const auto avx_precomputed = blst::detail::precompute_middle_product_b_public<AvxCtx>(
        std::span<const G2Point>(b_points),
        a_points.size(),
        nullptr);

    const auto scalar_a = blst::detail::convert_g1_inputs<ScalarCtx>(std::span<const G1Point>(a_points));
    const auto avx_a = blst::detail::convert_g1_inputs<AvxCtx>(std::span<const G1Point>(a_points));

    std::vector<ScalarCtx::G1Point> scalar_a_fft(
        scalar_precomputed.fft_size,
        ScalarCtx::g1_infinity());
    std::vector<AvxCtx::G1Point> avx_a_fft(
        avx_precomputed.fft_size,
        AvxCtx::g1_infinity());
    std::copy(scalar_a.begin(), scalar_a.end(), scalar_a_fft.begin());
    std::copy(avx_a.begin(), avx_a.end(), avx_a_fft.begin());

    const std::size_t xn_a = fft_admissible_size(a_points.size(), scalar_precomputed.fft_size);
    fft_forward_truncated_prefix<ScalarCtx::G1Point, Fr, ScalarCtx::G1CoeffCtx>(
        std::span<ScalarCtx::G1Point>(scalar_a_fft),
        std::span<const Fr>(scalar_precomputed.roots),
        xn_a,
        nullptr,
        "backend_alias.scalar_fft");
    fft_forward_truncated_prefix<AvxCtx::G1Point, Fr, AvxCtx::G1CoeffCtx>(
        std::span<AvxCtx::G1Point>(avx_a_fft),
        std::span<const Fr>(avx_precomputed.roots),
        xn_a,
        nullptr,
        "backend_alias.avx_fft");

    std::vector<ScalarCtx::GtFftElement> scalar_stage(
        scalar_precomputed.fft_size,
        ScalarCtx::identity());
    std::vector<AvxCtx::GtFftElement> avx_stage(
        avx_precomputed.fft_size,
        AvxCtx::identity());
    for (std::size_t i = 0; i < scalar_precomputed.fft_size; ++i) {
        scalar_stage[i] = ScalarCtx::pointwise_pair_direct(scalar_a_fft[i], scalar_precomputed.b_fft[i]);
        avx_stage[i] = AvxCtx::pointwise_pair_direct(avx_a_fft[i], avx_precomputed.b_fft[i]);
    }

    bit_reverse(std::span<ScalarCtx::GtFftElement>(scalar_stage));
    bit_reverse(std::span<AvxCtx::GtFftElement>(avx_stage));

    std::vector<Fr> inv_roots_copy(
        avx_precomputed.inv_roots.scalar_span().begin(),
        avx_precomputed.inv_roots.scalar_span().end());

    for (std::size_t len = 2, stage = 0; len <= scalar_stage.size(); len <<= 1, ++stage) {
        const std::size_t half = len >> 1;
        for (std::size_t i = 0; i < scalar_stage.size(); i += len) {
            {
                const auto u_scalar = scalar_stage[i];
                const auto t_scalar = scalar_stage[i + half];
                scalar_stage[i] = ScalarCtx::GtCoeffCtx::plus(u_scalar, t_scalar);
                scalar_stage[i + half] = ScalarCtx::GtCoeffCtx::minus(u_scalar, t_scalar);

                const auto u_avx = avx_stage[i];
                const auto t_avx = avx_stage[i + half];
                avx_stage[i] = blst::avx::AvxGtCoeffCtx::plus(u_avx, t_avx);
                avx_stage[i + half] = blst::avx::AvxGtCoeffCtx::minus(u_avx, t_avx);
            }

            for (std::size_t j = 1; j < half; ++j) {
                const Fr& w_scalar = scalar_precomputed.inv_roots[j * (scalar_stage.size() / len)];
                const Fr& w_avx = inv_roots_copy[j * (avx_stage.size() / len)];

                const auto u_scalar = scalar_stage[i + j];
                const auto t_scalar =
                    ScalarCtx::GtCoeffCtx::mult(scalar_stage[i + j + half], w_scalar);
                const auto u_avx = avx_stage[i + j];
                const auto t_avx =
                    blst::avx::AvxGtCoeffCtx::mult(avx_stage[i + j + half], w_avx);

                EXPECT_TRUE(GT::equal(
                    ScalarCtx::to_public_gt(u_scalar),
                    AvxCtx::to_public_gt(u_avx)))
                    << "Scalar/AVX stage input mismatch at stage=" << stage
                    << ", len=" << len
                    << ", block=" << i
                    << ", j=" << j;
                EXPECT_TRUE(GT::equal(
                    ScalarCtx::to_public_gt(t_scalar),
                    AvxCtx::to_public_gt(t_avx)))
                    << "Scalar/AVX twiddle-product mismatch at stage=" << stage
                    << ", len=" << len
                    << ", block=" << i
                    << ", j=" << j;

                scalar_stage[i + j] = ScalarCtx::GtCoeffCtx::plus(u_scalar, t_scalar);
                scalar_stage[i + j + half] = ScalarCtx::GtCoeffCtx::minus(u_scalar, t_scalar);
                avx_stage[i + j] = blst::avx::AvxGtCoeffCtx::plus(u_avx, t_avx);
                avx_stage[i + j + half] = blst::avx::AvxGtCoeffCtx::minus(u_avx, t_avx);
            }
        }

        for (std::size_t i = 0; i < scalar_stage.size(); ++i) {
            EXPECT_TRUE(GT::equal(
                ScalarCtx::to_public_gt(scalar_stage[i]),
                AvxCtx::to_public_gt(avx_stage[i])))
                << "Scalar/AVX stage mismatch at stage=" << stage
                << ", len=" << len
                << ", fft index=" << i;
        }
    }
}
#endif
