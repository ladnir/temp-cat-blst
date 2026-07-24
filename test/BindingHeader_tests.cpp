#include <gtest/gtest.h>

#include <array>
#include <cstring>
#include <vector>

#include "bindings/blst.h"
#include "bindings/blst_btx.h"
#include "bindings/blst_fft.h"
#include "bindings/blst_trx.h"

namespace {

blst_fr test_fr_from_u64(uint64_t value)
{
    const std::array<uint64_t, 4> words = {value, 0, 0, 0};
    blst_fr out = {};
    blst_fr_from_uint64(&out, words.data());
    return out;
}

std::vector<blst_fr> test_full_width_scalars(size_t count)
{
    const blst_fr multiplier = test_fr_from_u64(6364136223846793005ULL);
    const blst_fr increment = test_fr_from_u64(1442695040888963407ULL);
    blst_fr state = test_fr_from_u64(17);
    std::vector<blst_fr> scalars(count);
    for (auto& scalar : scalars) {
        scalar = state;
        blst_fr_mul(&state, &state, &multiplier);
        blst_fr_add(&state, &state, &increment);
    }
    return scalars;
}

std::vector<blst_p1_affine> test_distinct_g1_points(size_t count)
{
    std::vector<blst_p1> projective(count);
    std::vector<blst_p1_affine> affine(count);
    if (count == 0) {
        return affine;
    }
    const blst_p1 generator = *blst_p1_generator();
    blst_p1 current = generator;
    for (auto& point : projective) {
        point = current;
        blst_p1 next = {};
        blst_p1_add_or_double(&next, &current, &generator);
        current = next;
    }
    const blst_p1* descriptors[] = {projective.data(), nullptr};
    blst_p1s_to_affine(affine.data(), descriptors, count);
    return affine;
}

std::vector<blst_p2_affine> test_distinct_g2_points(size_t count)
{
    std::vector<blst_p2> projective(count);
    std::vector<blst_p2_affine> affine(count);
    if (count == 0) {
        return affine;
    }
    const blst_p2 generator = *blst_p2_generator();
    blst_p2 current = generator;
    for (auto& point : projective) {
        point = current;
        blst_p2 next = {};
        blst_p2_add_or_double(&next, &current, &generator);
        current = next;
    }
    const blst_p2* descriptors[] = {projective.data(), nullptr};
    blst_p2s_to_affine(affine.data(), descriptors, count);
    return affine;
}

void expect_batch_affine_g1_matches(
    const std::vector<blst_p1_affine>& points,
    const std::vector<blst_fr>& scalars)
{
    ASSERT_EQ(points.size(), scalars.size());
    const size_t pippenger_bytes =
        blst_p1s_mult_pippenger_fr_scratch_sizeof(points.size());
    std::vector<limb_t> pippenger_scratch(
        pippenger_bytes / sizeof(limb_t));
    blst_p1 expected = {};
    ASSERT_TRUE(blst_p1s_mult_pippenger_fr(
        &expected, points.data(), scalars.data(), points.size(),
        pippenger_scratch.data(), pippenger_bytes));

    for (const size_t window : {size_t{1}, size_t{2}, size_t{7},
                                size_t{10}, size_t{14}}) {
        for (const size_t batch_size : {size_t{128}, size_t{384}}) {
            const size_t scratch_bytes =
                blst_p1s_mult_batch_affine_fr_scratch_sizeof(
                    points.size(), window, batch_size);
            ASSERT_NE(scratch_bytes, 0);
            std::vector<limb_t> scratch(
                scratch_bytes / sizeof(limb_t), ~limb_t{0});
            for (size_t reuse = 0; reuse < 2; ++reuse) {
                blst_p1 actual = {};
                ASSERT_TRUE(blst_p1s_mult_batch_affine_fr(
                    &actual, points.data(), scalars.data(), points.size(),
                    scratch.data(), scratch_bytes, window, batch_size));
                EXPECT_TRUE(blst_p1_is_equal(&actual, &expected))
                    << "window=" << window
                    << " batch_size=" << batch_size
                    << " reuse=" << reuse;
            }
        }
    }
}

} // namespace

TEST(BindingHeader, BlstUmbrellaExposesKeygenAndPkApis)
{
    std::array<byte, 32> ikm = {};
    for (size_t i = 0; i < ikm.size(); ++i) {
        ikm[i] = static_cast<byte>(i + 1);
    }

    blst_scalar sk = {};
    blst_p1 pk_g1 = {};
    blst_p2 pk_g2 = {};

    blst_keygen(&sk, ikm.data(), ikm.size(), nullptr, 0);
    blst_sk_to_pk_in_g1(&pk_g1, &sk);
    blst_sk_to_pk_in_g2(&pk_g2, &sk);

    EXPECT_FALSE(blst_p1_is_inf(&pk_g1));
    EXPECT_TRUE(blst_p1_in_g1(&pk_g1));
    EXPECT_FALSE(blst_p2_is_inf(&pk_g2));
    EXPECT_TRUE(blst_p2_in_g2(&pk_g2));
}

TEST(BindingHeader, KeygenIsDeterministicForSameInputs)
{
    std::array<byte, 64> ikm = {};
    for (size_t i = 0; i < ikm.size(); ++i) {
        ikm[i] = static_cast<byte>(0x7f ^ i);
    }

    blst_scalar sk1 = {};
    blst_scalar sk2 = {};
    blst_p1 pk1 = {};
    blst_p1 pk2 = {};

    blst_keygen(&sk1, ikm.data(), ikm.size(), nullptr, 0);
    blst_keygen(&sk2, ikm.data(), ikm.size(), nullptr, 0);

    EXPECT_EQ(std::memcmp(&sk1, &sk2, sizeof(sk1)), 0);

    blst_sk_to_pk_in_g1(&pk1, &sk1);
    blst_sk_to_pk_in_g1(&pk2, &sk2);
    EXPECT_TRUE(blst_p1_is_equal(&pk1, &pk2));
}

TEST(BindingHeader, DeriveMasterIsDeterministicForSameInputs)
{
    std::array<byte, 64> ikm = {};
    for (size_t i = 0; i < ikm.size(); ++i) {
        ikm[i] = static_cast<byte>(0xa5 + i);
    }

    blst_scalar sk1 = {};
    blst_scalar sk2 = {};

    blst_derive_master_eip2333(&sk1, ikm.data(), ikm.size());
    blst_derive_master_eip2333(&sk2, ikm.data(), ikm.size());

    EXPECT_EQ(std::memcmp(&sk1, &sk2, sizeof(sk1)), 0);
}

TEST(BindingHeader, ScalarHelperSurfaceIsExposedAndStable)
{
    EXPECT_EQ(static_cast<int>(BLST_BAD_SCALAR), 7);

    std::array<byte, 1> one_be = {1};
    blst_scalar one = {};
    blst_scalar two = {};
    blst_scalar zero = {};

    EXPECT_TRUE(blst_scalar_from_be_bytes(&one, one_be.data(), one_be.size()));
    EXPECT_TRUE(blst_sk_add_n_check(&two, &one, &one));
    EXPECT_TRUE(blst_sk_sub_n_check(&one, &two, &one));
    EXPECT_TRUE(blst_sk_mul_n_check(&two, &one, &two));
    EXPECT_FALSE(blst_scalar_from_be_bytes(&zero, nullptr, 0));
}

TEST(BindingHeader, PreparedFrFftUsesCallerOwnedDomain)
{
    constexpr size_t n = 8;
    EXPECT_EQ(blst_fr_fft_domain_len(0), 0);
    EXPECT_EQ(blst_fr_fft_domain_len(3), 0);
    EXPECT_EQ(blst_fr_fft_domain_len(n), n + 1);

    std::array<blst_fr, n + 1> domain = {};
    EXPECT_FALSE(blst_fr_fft_domain_init(nullptr, domain.size(), n));
    EXPECT_FALSE(blst_fr_fft_domain_init(domain.data(), n, n));
    ASSERT_TRUE(blst_fr_fft_domain_init(domain.data(), domain.size(), n));

    std::array<blst_fr, n> values = {};
    for (size_t i = 0; i < values.size(); ++i) {
        std::array<uint64_t, 4> scalar = {i + 1, 0, 0, 0};
        blst_fr_from_uint64(&values[i], scalar.data());
    }
    const auto input = values;

    EXPECT_FALSE(blst_fr_fft_with_domain(
        nullptr, n, domain.data(), domain.size(), false));
    EXPECT_FALSE(blst_fr_fft_with_domain(
        values.data(), n, nullptr, domain.size(), false));
    EXPECT_FALSE(blst_fr_fft_with_domain(
        values.data(), n, domain.data(), n, false));

    ASSERT_TRUE(blst_fr_fft_with_domain(
        values.data(), n, domain.data(), domain.size(), false));
    ASSERT_TRUE(blst_fr_fft_with_domain(
        values.data(), n, domain.data(), domain.size(), true));
    EXPECT_EQ(std::memcmp(values.data(), input.data(), sizeof(values)), 0);
}

TEST(BindingHeader, NativeFrPippengerUsesCallerOwnedScratch)
{
    constexpr size_t n = 4;
    std::array<blst_fr, n> scalars = {};
    for (size_t i = 0; i < n; ++i) {
        const std::array<uint64_t, 4> words = {i + 1, 0, 0, 0};
        blst_fr_from_uint64(&scalars[i], words.data());
    }
    const std::array<byte, 32> ten = {10};

    EXPECT_EQ(blst_p1s_mult_pippenger_fr_scratch_sizeof(0), 0);
    EXPECT_EQ(blst_p1s_mult_pippenger_fr_scratch_sizeof(1), 32);
    EXPECT_EQ(blst_p2s_mult_pippenger_fr_scratch_sizeof(0), 0);
    EXPECT_EQ(blst_p2s_mult_pippenger_fr_scratch_sizeof(1), 32);

    blst_p1 empty_g1 = *blst_p1_generator();
    blst_p2 empty_g2 = *blst_p2_generator();
    ASSERT_TRUE(blst_p1s_mult_pippenger_fr(
        &empty_g1, nullptr, nullptr, 0, nullptr, 0));
    ASSERT_TRUE(blst_p2s_mult_pippenger_fr(
        &empty_g2, nullptr, nullptr, 0, nullptr, 0));
    EXPECT_TRUE(blst_p1_is_inf(&empty_g1));
    EXPECT_TRUE(blst_p2_is_inf(&empty_g2));

    std::array<blst_p1_affine, n> g1_points = {};
    std::array<blst_p2_affine, n> g2_points = {};
    g1_points.fill(*blst_p1_affine_generator());
    g2_points.fill(*blst_p2_affine_generator());

    const size_t g1_scratch_size =
        blst_p1s_mult_pippenger_fr_scratch_sizeof(n);
    const size_t g2_scratch_size =
        blst_p2s_mult_pippenger_fr_scratch_sizeof(n);
    ASSERT_NE(g1_scratch_size, 0);
    ASSERT_NE(g2_scratch_size, 0);
    ASSERT_EQ(g1_scratch_size % sizeof(limb_t), 0);
    ASSERT_EQ(g2_scratch_size % sizeof(limb_t), 0);
    std::vector<limb_t> g1_scratch(g1_scratch_size / sizeof(limb_t));
    std::vector<limb_t> g2_scratch(g2_scratch_size / sizeof(limb_t));

    blst_p1 actual_g1 = {};
    blst_p1 expected_g1 = {};
    blst_p1_mult(&expected_g1, blst_p1_generator(), ten.data(), 255);
    EXPECT_FALSE(blst_p1s_mult_pippenger_fr(
        &actual_g1, g1_points.data(), scalars.data(), n,
        g1_scratch.data(), g1_scratch_size - 1));
    ASSERT_TRUE(blst_p1s_mult_pippenger_fr(
        &actual_g1, g1_points.data(), scalars.data(), n,
        g1_scratch.data(), g1_scratch_size));
    EXPECT_TRUE(blst_p1_is_equal(&actual_g1, &expected_g1));

    std::array<blst_scalar, n> canonical = {};
    EXPECT_TRUE(blst_scalars_from_fr(nullptr, nullptr, 0));
    EXPECT_FALSE(blst_scalars_from_fr(nullptr, scalars.data(), n));
    EXPECT_FALSE(blst_scalars_from_fr(canonical.data(), nullptr, n));
    ASSERT_TRUE(blst_scalars_from_fr(canonical.data(), scalars.data(), n));
    const blst_p1_affine* g1_point_descriptors[] = {
        g1_points.data(), nullptr};
    const byte* scalar_descriptors[] = {
        reinterpret_cast<const byte*>(canonical.data()), nullptr};
    const size_t canonical_g1_scratch_size =
        blst_p1s_mult_pippenger_scratch_sizeof(n);
    std::vector<limb_t> canonical_g1_scratch(
        canonical_g1_scratch_size / sizeof(limb_t));
    blst_p1 canonical_g1 = {};
    blst_p1s_mult_pippenger(
        &canonical_g1, g1_point_descriptors, n, scalar_descriptors, 255,
        canonical_g1_scratch.data());
    EXPECT_TRUE(blst_p1_is_equal(&actual_g1, &canonical_g1));

    blst_p2 actual_g2 = {};
    blst_p2 expected_g2 = {};
    blst_p2_mult(&expected_g2, blst_p2_generator(), ten.data(), 255);
    EXPECT_FALSE(blst_p2s_mult_pippenger_fr(
        &actual_g2, g2_points.data(), scalars.data(), n,
        g2_scratch.data(), g2_scratch_size - 1));
    ASSERT_TRUE(blst_p2s_mult_pippenger_fr(
        &actual_g2, g2_points.data(), scalars.data(), n,
        g2_scratch.data(), g2_scratch_size));
    EXPECT_TRUE(blst_p2_is_equal(&actual_g2, &expected_g2));

    const blst_p2_affine* g2_point_descriptors[] = {
        g2_points.data(), nullptr};
    const size_t canonical_g2_scratch_size =
        blst_p2s_mult_pippenger_scratch_sizeof(n);
    std::vector<limb_t> canonical_g2_scratch(
        canonical_g2_scratch_size / sizeof(limb_t));
    blst_p2 canonical_g2 = {};
    blst_p2s_mult_pippenger(
        &canonical_g2, g2_point_descriptors, n, scalar_descriptors, 255,
        canonical_g2_scratch.data());
    EXPECT_TRUE(blst_p2_is_equal(&actual_g2, &canonical_g2));
}

TEST(BindingHeader, NativeFrBatchAffineG1MatchesPippenger)
{
    constexpr size_t n = 513;
    constexpr size_t window = 10;

    const auto fr_from_u64 = [](uint64_t value) {
        const std::array<uint64_t, 4> words = {value, 0, 0, 0};
        blst_fr out = {};
        blst_fr_from_uint64(&out, words.data());
        return out;
    };
    const blst_fr multiplier = fr_from_u64(6364136223846793005ULL);
    const blst_fr increment = fr_from_u64(1442695040888963407ULL);
    blst_fr state = fr_from_u64(17);

    std::array<blst_fr, n> scalars = {};
    for (auto& scalar : scalars) {
        scalar = state;
        blst_fr_mul(&state, &state, &multiplier);
        blst_fr_add(&state, &state, &increment);
    }
    scalars[11] = {};

    std::array<blst_p1_affine, n> points = {};
    points.fill(*blst_p1_affine_generator());

    std::array<blst_scalar, n> canonical = {};
    ASSERT_TRUE(blst_scalars_from_fr(
        canonical.data(), scalars.data(), scalars.size()));
    const blst_p1_affine* point_descriptors[] = {
        points.data(), nullptr};
    const byte* scalar_descriptors[] = {
        reinterpret_cast<const byte*>(canonical.data()), nullptr};
    const size_t pippenger_scratch_size =
        blst_p1s_mult_pippenger_scratch_sizeof(n);
    ASSERT_NE(pippenger_scratch_size, 0);
    std::vector<limb_t> pippenger_scratch(
        pippenger_scratch_size / sizeof(limb_t));
    blst_p1 expected = {};
    blst_p1s_mult_pippenger(
        &expected, point_descriptors, n, scalar_descriptors, 255,
        pippenger_scratch.data());

    EXPECT_EQ(
        blst_p1s_mult_batch_affine_fr_scratch_sizeof(n, 0, 128), 0);
    EXPECT_EQ(
        blst_p1s_mult_batch_affine_fr_scratch_sizeof(n, window, 64), 0);

    for (const size_t batch_size : {size_t{128}, size_t{384}}) {
        const size_t scratch_size =
            blst_p1s_mult_batch_affine_fr_scratch_sizeof(
                n, window, batch_size);
        ASSERT_NE(scratch_size, 0);
        ASSERT_EQ(scratch_size % sizeof(limb_t), 0);
        std::vector<limb_t> scratch(scratch_size / sizeof(limb_t));

        blst_p1 actual = {};
        EXPECT_FALSE(blst_p1s_mult_batch_affine_fr(
            &actual, points.data(), scalars.data(), n,
            scratch.data(), scratch_size - 1, window, batch_size));
        ASSERT_TRUE(blst_p1s_mult_batch_affine_fr(
            &actual, points.data(), scalars.data(), n,
            scratch.data(), scratch_size, window, batch_size));
        EXPECT_TRUE(blst_p1_is_equal(&actual, &expected));
    }

    blst_p1 empty = *blst_p1_generator();
    EXPECT_EQ(
        blst_p1s_mult_batch_affine_fr_scratch_sizeof(0, window, 128), 0);
    ASSERT_TRUE(blst_p1s_mult_batch_affine_fr(
        &empty, nullptr, nullptr, 0, nullptr, 0, window, 128));
    EXPECT_TRUE(blst_p1_is_inf(&empty));

    const size_t single_scratch_size =
        blst_p1s_mult_batch_affine_fr_scratch_sizeof(1, window, 128);
    EXPECT_EQ(single_scratch_size, sizeof(blst_scalar));
    std::vector<limb_t> single_scratch(
        single_scratch_size / sizeof(limb_t));
    blst_p1 single = {};
    ASSERT_TRUE(blst_p1s_mult_batch_affine_fr(
        &single, points.data(), scalars.data(), 1,
        single_scratch.data(), single_scratch_size, window, 128));
    blst_scalar first_scalar = {};
    blst_scalar_from_fr(&first_scalar, &scalars[0]);
    blst_p1 single_expected = {};
    blst_p1_mult(
        &single_expected, blst_p1_generator(),
        first_scalar.b, 255);
    EXPECT_TRUE(blst_p1_is_equal(&single, &single_expected));
}

TEST(BindingHeader, NativeFrBatchAffineG1DifferentialDistributions)
{
    constexpr size_t n = 257;
    const auto distinct_points = test_distinct_g1_points(n);
    auto full_width = test_full_width_scalars(n);
    full_width[0] = {};
    full_width[128] = {};
    expect_batch_affine_g1_matches(distinct_points, full_width);

    std::vector<blst_p1_affine> repeated_points(
        n, *blst_p1_affine_generator());
    expect_batch_affine_g1_matches(repeated_points, full_width);

    std::vector<blst_fr> all_zero(n);
    expect_batch_affine_g1_matches(distinct_points, all_zero);

    std::vector<blst_fr> collision_heavy(
        n, test_fr_from_u64(0x15555));
    expect_batch_affine_g1_matches(distinct_points, collision_heavy);
}

TEST(BindingHeader, NativeFrBatchAffineG2MatchesPippenger)
{
    constexpr size_t n = 257;
    const auto points = test_distinct_g2_points(n);
    auto scalars = test_full_width_scalars(n);
    scalars[0] = {};
    scalars[128] = {};

    const size_t pippenger_bytes =
        blst_p2s_mult_pippenger_fr_scratch_sizeof(n);
    std::vector<limb_t> pippenger_scratch(
        pippenger_bytes / sizeof(limb_t));
    blst_p2 expected = {};
    ASSERT_TRUE(blst_p2s_mult_pippenger_fr(
        &expected, points.data(), scalars.data(), n,
        pippenger_scratch.data(), pippenger_bytes));

    for (const size_t window : {size_t{1}, size_t{7}, size_t{14}}) {
        for (const size_t batch_size : {size_t{128}, size_t{384}}) {
            const size_t scratch_bytes =
                blst_p2s_mult_batch_affine_fr_scratch_sizeof(
                    n, window, batch_size);
            ASSERT_NE(scratch_bytes, 0);
            std::vector<limb_t> scratch(
                scratch_bytes / sizeof(limb_t), ~limb_t{0});
            blst_p2 actual = {};
            ASSERT_TRUE(blst_p2s_mult_batch_affine_fr(
                &actual, points.data(), scalars.data(), n,
                scratch.data(), scratch_bytes, window, batch_size));
            EXPECT_TRUE(blst_p2_is_equal(&actual, &expected))
                << "window=" << window
                << " batch_size=" << batch_size;
        }
    }
}

TEST(BindingHeader, SelectedNativeFrMsmScratchIsCapacityMonotone)
{
    size_t previous_g1 = 0;
    size_t previous_g2 = 0;
    for (size_t n = 1; n <= (size_t{1} << 18) + 1; ++n) {
        const size_t current_g1 = blst_p1s_mult_fr_scratch_sizeof(n);
        const size_t current_g2 = blst_p2s_mult_fr_scratch_sizeof(n);
        ASSERT_NE(current_g1, 0);
        ASSERT_NE(current_g2, 0);
        EXPECT_GE(current_g1, previous_g1) << "n=" << n;
        EXPECT_GE(current_g2, previous_g2) << "n=" << n;
        previous_g1 = current_g1;
        previous_g2 = current_g2;
    }
}

TEST(BindingHeader, SelectedNativeFrG1MsmHandlesCapacityAndIdentity)
{
    constexpr size_t capacity = size_t{1} << 12;
    constexpr size_t n = size_t{1} << 11;
    auto points = test_distinct_g1_points(capacity);
    const auto scalars = test_full_width_scalars(capacity);

    const size_t selected_bytes =
        blst_p1s_mult_fr_scratch_sizeof(capacity);
    std::vector<limb_t> selected_scratch(
        selected_bytes / sizeof(limb_t), ~limb_t{0});
    const size_t pippenger_bytes =
        blst_p1s_mult_pippenger_fr_scratch_sizeof(n);
    std::vector<limb_t> pippenger_scratch(
        pippenger_bytes / sizeof(limb_t));

    blst_p1 expected = {};
    blst_p1 actual = {};
    ASSERT_TRUE(blst_p1s_mult_pippenger_fr(
        &expected, points.data(), scalars.data(), n,
        pippenger_scratch.data(), pippenger_bytes));
    ASSERT_TRUE(blst_p1s_mult_fr(
        &actual, points.data(), scalars.data(), n,
        selected_scratch.data(), selected_bytes));
    EXPECT_TRUE(blst_p1_is_equal(&actual, &expected));

    points[17] = {};
    ASSERT_TRUE(blst_p1s_mult_pippenger_fr(
        &expected, points.data(), scalars.data(), capacity,
        selected_scratch.data(), selected_bytes));
    EXPECT_FALSE(blst_p1s_mult_fr(
        &actual, points.data(), scalars.data(), capacity,
        selected_scratch.data(), selected_bytes - 1));
    ASSERT_TRUE(blst_p1s_mult_fr(
        &actual, points.data(), scalars.data(), capacity,
        selected_scratch.data(), selected_bytes));
    EXPECT_TRUE(blst_p1_is_equal(&actual, &expected));
}

TEST(BindingHeader, SelectedNativeFrG2MsmHandlesBatchRange)
{
    constexpr size_t n = size_t{1} << 12;
    const auto points = test_distinct_g2_points(n);
    const auto scalars = test_full_width_scalars(n);
    const size_t selected_bytes =
        blst_p2s_mult_fr_scratch_sizeof(n);
    std::vector<limb_t> selected_scratch(
        selected_bytes / sizeof(limb_t), ~limb_t{0});
    const size_t pippenger_bytes =
        blst_p2s_mult_pippenger_fr_scratch_sizeof(n);
    std::vector<limb_t> pippenger_scratch(
        pippenger_bytes / sizeof(limb_t));

    blst_p2 expected = {};
    blst_p2 actual = {};
    ASSERT_TRUE(blst_p2s_mult_pippenger_fr(
        &expected, points.data(), scalars.data(), n,
        pippenger_scratch.data(), pippenger_bytes));
    ASSERT_TRUE(blst_p2s_mult_fr(
        &actual, points.data(), scalars.data(), n,
        selected_scratch.data(), selected_bytes));
    EXPECT_TRUE(blst_p2_is_equal(&actual, &expected));
}
