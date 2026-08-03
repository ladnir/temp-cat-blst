#include <gtest/gtest.h>

#include <array>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>
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

blst_fp test_fp_from_u64(uint64_t value)
{
    const std::array<uint64_t, 6> words = {value, 0, 0, 0, 0, 0};
    blst_fp out = {};
    blst_fp_from_uint64(&out, words.data());
    return out;
}

blst_fp test_fp_from_bytes_oracle(
    const std::vector<byte>& input,
    bool little_endian)
{
    const blst_fp radix = test_fp_from_u64(256);
    blst_fp acc = test_fp_from_u64(0);

    for (size_t i = 0; i < input.size(); ++i) {
        const size_t index = little_endian ? input.size() - 1 - i : i;
        const blst_fp digit = test_fp_from_u64(input[index]);
        blst_fp shifted = {};
        blst_fp_mul(&shifted, &acc, &radix);
        blst_fp_add(&acc, &shifted, &digit);
    }
    return acc;
}

std::array<byte, 48> test_fp_to_bendian(const blst_fp& value)
{
    std::array<byte, 48> out = {};
    blst_bendian_from_fp(out.data(), &value);
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

static_assert(std::is_same_v<
    decltype(blst_scalar_from_hexascii(nullptr, nullptr)), bool>);
static_assert(std::is_same_v<
    decltype(blst_fr_from_hexascii(nullptr, nullptr)), bool>);
static_assert(std::is_same_v<
    decltype(blst_fp_from_hexascii(nullptr, nullptr)), bool>);

TEST(BindingHeader, HexParsersRequireACompleteRepresentableValue)
{
    const auto ascii = [](const char* text) {
        return reinterpret_cast<const byte*>(text);
    };

    blst_scalar scalar = {};
    ASSERT_TRUE(blst_scalar_from_hexascii(&scalar, ascii("0x2a")));
    blst_scalar expected_scalar = {};
    expected_scalar.b[0] = 0x2a;
    EXPECT_EQ(std::memcmp(&scalar, &expected_scalar, sizeof(scalar)), 0);

    blst_fr fr = {};
    ASSERT_TRUE(blst_fr_from_hexascii(&fr, ascii("2A")));
    const blst_fr expected_fr = test_fr_from_u64(42);
    EXPECT_EQ(std::memcmp(&fr, &expected_fr, sizeof(fr)), 0);

    blst_fp fp = {};
    ASSERT_TRUE(blst_fp_from_hexascii(&fp, ascii("0x2a")));
    const blst_fp expected_fp = test_fp_from_u64(42);
    EXPECT_EQ(std::memcmp(&fp, &expected_fp, sizeof(fp)), 0);

    for (const char* invalid : {"", "0x", "12xyz", "-1", " 2a"}) {
        SCOPED_TRACE(invalid);
        std::memset(&scalar, 0xa5, sizeof(scalar));
        std::memset(&fr, 0xa5, sizeof(fr));
        std::memset(&fp, 0xa5, sizeof(fp));
        const blst_scalar scalar_before = scalar;
        const blst_fr fr_before = fr;
        const blst_fp fp_before = fp;

        EXPECT_FALSE(blst_scalar_from_hexascii(&scalar, ascii(invalid)));
        EXPECT_FALSE(blst_fr_from_hexascii(&fr, ascii(invalid)));
        EXPECT_FALSE(blst_fp_from_hexascii(&fp, ascii(invalid)));
        EXPECT_EQ(std::memcmp(&scalar, &scalar_before, sizeof(scalar)), 0);
        EXPECT_EQ(std::memcmp(&fr, &fr_before, sizeof(fr)), 0);
        EXPECT_EQ(std::memcmp(&fp, &fp_before, sizeof(fp)), 0);
    }

    const std::string scalar_too_wide(65, '1');
    const std::string fp_too_wide(97, '1');
    EXPECT_FALSE(blst_scalar_from_hexascii(&scalar, ascii(scalar_too_wide.c_str())));
    EXPECT_FALSE(blst_fr_from_hexascii(&fr, ascii(scalar_too_wide.c_str())));
    EXPECT_FALSE(blst_fp_from_hexascii(&fp, ascii(fp_too_wide.c_str())));

    EXPECT_FALSE(blst_scalar_from_hexascii(nullptr, ascii("1")));
    EXPECT_FALSE(blst_fr_from_hexascii(nullptr, ascii("1")));
    EXPECT_FALSE(blst_fp_from_hexascii(nullptr, ascii("1")));
    EXPECT_FALSE(blst_scalar_from_hexascii(&scalar, nullptr));
    EXPECT_FALSE(blst_fr_from_hexascii(&fr, nullptr));
    EXPECT_FALSE(blst_fp_from_hexascii(&fp, nullptr));
}

TEST(BindingHeader, VariableLengthFpDecodersMatchBytewiseOracle)
{
    const blst_fp zero = test_fp_from_u64(0);
    blst_fp actual_le = test_fp_from_u64(7);
    blst_fp actual_be = test_fp_from_u64(9);
    blst_fp_from_le_bytes(&actual_le, nullptr, 0);
    blst_fp_from_be_bytes(&actual_be, nullptr, 0);
    EXPECT_EQ(test_fp_to_bendian(actual_le), test_fp_to_bendian(zero));
    EXPECT_EQ(test_fp_to_bendian(actual_be), test_fp_to_bendian(zero));

    for (const size_t len : {size_t{1}, size_t{2}, size_t{47}, size_t{48},
                             size_t{49}, size_t{95}, size_t{96}, size_t{97},
                             size_t{255}, size_t{511}, size_t{4097}}) {
        SCOPED_TRACE(::testing::Message() << "len=" << len);
        std::vector<byte> input(len);
        for (size_t i = 0; i < len; ++i) {
            input[i] = static_cast<byte>((i * 131 + len * 17 + 0x5b) & 0xff);
        }

        blst_fp_from_le_bytes(&actual_le, input.data(), input.size());
        blst_fp_from_be_bytes(&actual_be, input.data(), input.size());
        const blst_fp expected_le = test_fp_from_bytes_oracle(input, true);
        const blst_fp expected_be = test_fp_from_bytes_oracle(input, false);
        EXPECT_EQ(test_fp_to_bendian(actual_le), test_fp_to_bendian(expected_le));
        EXPECT_EQ(test_fp_to_bendian(actual_be), test_fp_to_bendian(expected_be));

        std::reverse(input.begin(), input.end());
        blst_fp reversed_be = {};
        blst_fp_from_be_bytes(&reversed_be, input.data(), input.size());
        EXPECT_EQ(test_fp_to_bendian(actual_le), test_fp_to_bendian(reversed_be));
    }

    for (const size_t len : {size_t{48}, size_t{49}, size_t{96}}) {
        std::vector<byte> input(len, 0xff);
        blst_fp_from_le_bytes(&actual_le, input.data(), input.size());
        blst_fp_from_be_bytes(&actual_be, input.data(), input.size());
        EXPECT_EQ(test_fp_to_bendian(actual_le),
                  test_fp_to_bendian(test_fp_from_bytes_oracle(input, true)));
        EXPECT_EQ(test_fp_to_bendian(actual_be),
                  test_fp_to_bendian(test_fp_from_bytes_oracle(input, false)));
    }
}

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

    const std::array<uint32_t, 8> words32 = {1, 2, 3, 4, 5, 6, 7, 8};
    std::array<uint32_t, 8> roundtrip32 = {};
    blst_scalar_from_uint32(&one, words32.data());
    blst_uint32_from_scalar(roundtrip32.data(), &one);
    EXPECT_EQ(roundtrip32, words32);

    const std::array<uint64_t, 4> words64 = {1, 2, 3, 4};
    std::array<uint64_t, 4> roundtrip64 = {};
    blst_scalar_from_uint64(&one, words64.data());
    blst_uint64_from_scalar(roundtrip64.data(), &one);
    EXPECT_EQ(roundtrip64, words64);

    std::array<byte, 32> little = {};
    std::array<byte, 32> big = {};
    std::array<byte, 32> encoded = {};
    little.front() = 7;
    big.back() = 7;
    blst_scalar_from_lendian(&one, little.data());
    blst_lendian_from_scalar(encoded.data(), &one);
    EXPECT_EQ(encoded, little);
    blst_scalar_from_bendian(&one, big.data());
    blst_bendian_from_scalar(encoded.data(), &one);
    EXPECT_EQ(encoded, big);

    const blst_fr fr_one = test_fr_from_u64(1);
    const blst_fr fr_zero = {};
    blst_fr selected = {};
    blst_fr sum = {};
    blst_fr_cneg(&selected, &fr_one, false);
    EXPECT_EQ(std::memcmp(&selected, &fr_one, sizeof(fr_one)), 0);
    blst_fr_cneg(&selected, &fr_one, true);
    blst_fr_add(&sum, &selected, &fr_one);
    EXPECT_EQ(std::memcmp(&sum, &fr_zero, sizeof(sum)), 0);
}

TEST(BindingHeader, PreparedFrFftUsesCallerOwnedDomain)
{
    constexpr size_t n = 8;
    EXPECT_FALSE(blst_fr_fft(nullptr, 0, false));
    EXPECT_EQ(blst_fr_fft_domain_len(0), 0);
    EXPECT_EQ(blst_fr_fft_domain_len(3), 0);
    EXPECT_EQ(blst_fr_fft_domain_len(n), n + 2);

    std::array<blst_fr, n + 2> domain = {};
    std::array<blst_fr, n + 3> oversized_domain = {};
    EXPECT_FALSE(blst_fr_fft_domain_init(nullptr, domain.size(), n));
    EXPECT_FALSE(blst_fr_fft_domain_init(domain.data(), domain.size() - 1, n));
    EXPECT_FALSE(blst_fr_fft_domain_init(
        oversized_domain.data(), oversized_domain.size(), n));
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
        values.data(), n, domain.data(), domain.size() - 1, false));
    EXPECT_FALSE(blst_fr_fft_with_domain(
        values.data(), n, oversized_domain.data(), oversized_domain.size(), false));

    std::array<blst_fr, 4 + 2> domain_four = {};
    ASSERT_TRUE(blst_fr_fft_domain_init(domain_four.data(), domain_four.size(), 4));
    EXPECT_FALSE(blst_fr_fft_with_domain(
        values.data(), 4, domain.data(), domain_four.size(), false));

    auto corrupt_tag_domain = domain;
    corrupt_tag_domain[n + 1] = {};
    EXPECT_FALSE(blst_fr_fft_with_domain(
        values.data(), n, corrupt_tag_domain.data(), corrupt_tag_domain.size(), false));
    EXPECT_EQ(std::memcmp(values.data(), input.data(), sizeof(values)), 0);

    ASSERT_TRUE(blst_fr_fft_with_domain(
        values.data(), n, domain.data(), domain.size(), false));
    ASSERT_TRUE(blst_fr_fft_with_domain(
        values.data(), n, domain.data(), domain.size(), true));
    EXPECT_EQ(std::memcmp(values.data(), input.data(), sizeof(values)), 0);
}

TEST(BindingHeader, FrFftInterpolationRejectsOversubscribedDomain)
{
    constexpr size_t fft_len = 8;
    constexpr size_t n = fft_len + 1;
    std::array<size_t, n> indices = {};
    std::array<blst_fr, n> values = {};
    std::array<blst_fr, n> out = {};
    std::array<blst_fr, 1> scratch = {};

    EXPECT_EQ(blst_fr_interpolate_fft_domain_suffix_scratch_len(n, fft_len), 0);
    EXPECT_FALSE(blst_fr_interpolate_fft_domain_suffix_with_scratch(
        out.data(),
        indices.data(),
        values.data(),
        n,
        0,
        fft_len,
        scratch.data(),
        scratch.size()));
}

TEST(BindingHeader, FrFftSuffixScratchLengthIsChecked)
{
    constexpr size_t largest_supported_fft = size_t{1} << 30;
    const size_t required = blst_fr_interpolate_fft_domain_suffix_scratch_len(
        largest_supported_fft, largest_supported_fft);

    if constexpr (sizeof(size_t) == sizeof(uint32_t)) {
        EXPECT_EQ(required, 0);
    } else {
        EXPECT_EQ(required, static_cast<size_t>(6ULL * largest_supported_fft + 2));
    }
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
    EXPECT_FALSE(blst_scalars_from_fr(
        canonical.data(), scalars.data(),
        std::numeric_limits<size_t>::max() / sizeof(blst_scalar) + 1));
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
