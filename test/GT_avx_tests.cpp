#include <gtest/gtest.h>

#include <array>
#include <span>

#include "InsecureRng.h"
#include "blst/curve/GT_avx.h"
#include "blst/fields/Fr.h"

using namespace blst;

namespace
{
void rand_redc_fp(Fp& v, InsecureRng& rng)
{
    for (u64 i = 0; i < v.num_limb(); ++i) {
        v[i] = rng.next_u64();
    }
    while (Fp::is_redc(v) == 0) {
        Fp::rshift(v, v, 1);
    }
}

void rand_redc_fp2(Fp2& v, InsecureRng& rng)
{
    rand_redc_fp(v[0], rng);
    rand_redc_fp(v[1], rng);
}

void rand_redc_fp12(Fp12& v, InsecureRng& rng)
{
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            rand_redc_fp2(v[i][j], rng);
        }
    }
}

void make_random_gt(GT& out, InsecureRng& rng)
{
    static const std::array<GT, 16> pool = []() {
        std::array<GT, 16> values{};
        InsecureRng local_rng;
        for (size_t i = 0; i < values.size(); ++i) {
            Fp12 seed{};
            rand_redc_fp12(seed, local_rng);
            Fp12::final_exp(values[i].data, seed);
        }
        return values;
    }();

    static size_t next_idx = 0;
    (void)rng;
    out = pool[next_idx % pool.size()];
    ++next_idx;
}

void fill_scalar_little_endian(byte* out, const size_t stride, const size_t nbits, InsecureRng& rng)
{
    for (size_t i = 0; i < stride; ++i) {
        out[i] = static_cast<byte>(rng.next_u64());
    }

    if (nbits >= (stride * 8)) {
        return;
    }

    const size_t full_bytes = nbits / 8;
    const size_t rem_bits = nbits % 8;
    if (rem_bits != 0) {
        const byte mask = static_cast<byte>((1u << rem_bits) - 1u);
        out[full_bytes] = static_cast<byte>(out[full_bytes] & mask);
        for (size_t i = full_bytes + 1; i < stride; ++i) {
            out[i] = 0;
        }
        return;
    }

    for (size_t i = full_bytes; i < stride; ++i) {
        out[i] = 0;
    }
}

u64 load_le_u64(const byte* in)
{
    u64 v = 0;
    for (size_t i = 0; i < sizeof(v); ++i) {
        v |= static_cast<u64>(in[i]) << (8 * i);
    }
    return v;
}

bool scalar_ge_group_order(const byte* scalar_32)
{
    const Fr& order = get_BLST_r();
    for (int limb = 3; limb >= 0; --limb) {
        const u64 s = load_le_u64(scalar_32 + limb * 8);
        const u64 r = order[limb];
        if (s > r) {
            return true;
        }
        if (s < r) {
            return false;
        }
    }
    return true;
}

void scalar_rshift_1(byte* scalar_32)
{
    byte carry = 0;
    for (int i = 31; i >= 0; --i) {
        const byte next_carry = static_cast<byte>(scalar_32[i] & 1u);
        scalar_32[i] = static_cast<byte>((scalar_32[i] >> 1) | (carry << 7));
        carry = next_carry;
    }
}

void make_reduced_scalar_little_endian(byte* out, const size_t stride, const size_t nbits, InsecureRng& rng)
{
    fill_scalar_little_endian(out, stride, nbits, rng);
    if (stride < 32) {
        return;
    }
    while (scalar_ge_group_order(out)) {
        scalar_rshift_1(out);
    }
}

} // namespace

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

TEST(GTBatchAvx, PackRoundTripAndSelect)
{
    InsecureRng rng;
    std::array<GT, 8> a{};
    std::array<GT, 8> b{};
    std::array<GT, 8> out{};
    GTPack8x1 pa{};
    GTPack8x1 pb{};
    GTPack8x1 pr{};

    for (size_t lane = 0; lane < 8; ++lane) {
        make_random_gt(a[lane], rng);
        make_random_gt(b[lane], rng);
    }

    GTPack8x1::from_aos(pa, a.data());
    GTPack8x1::to_aos(out.data(), pa);
    for (size_t lane = 0; lane < 8; ++lane) {
        EXPECT_TRUE(GT::equal(out[lane], a[lane]));
    }

    GTPack8x1::from_aos(pb, b.data());
    const std::array<u8, 2> masks = {0x00, 0x93};
    for (const u8 mask : masks) {
        GTPack8x1::select(pr, pa, pb, mask);
        GTPack8x1::to_aos(out.data(), pr);
        for (size_t lane = 0; lane < 8; ++lane) {
            const GT& expected = ((mask >> lane) & 1u) ? b[lane] : a[lane];
            EXPECT_TRUE(GT::equal(out[lane], expected));
        }
    }
}

TEST(GTBatchAvx, PackPartialRoundTrip)
{
    InsecureRng rng;
    const GT fill = GT::identity();

    const std::array<size_t, 2> n_active_cases = {0, 8};
    for (const size_t n_active : n_active_cases) {
        std::array<GT, 8> in{};
        std::array<GT, 8> out{};
        std::array<GT, 8> full{};
        GTPack8x1 pack{};

        for (size_t lane = 0; lane < 8; ++lane) {
            make_random_gt(in[lane], rng);
        }

        GTPack8x1::from_aos_partial(pack, in.data(), n_active, fill);
        GTPack8x1::to_aos_partial(out.data(), n_active, pack);
        GTPack8x1::to_aos(full.data(), pack);

        for (size_t lane = 0; lane < n_active; ++lane) {
            EXPECT_TRUE(GT::equal(out[lane], in[lane]));
            EXPECT_TRUE(GT::equal(full[lane], in[lane]));
        }
        for (size_t lane = n_active; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(full[lane], fill));
        }
    }
}

static void run_add_matches_scalar()
{
    InsecureRng rng;

    for (int iter = 0; iter < 1; ++iter) {
        std::array<GT, 8> a{};
        std::array<GT, 8> b{};
        std::array<GT, 8> got{};
        std::array<GT, 8> expected{};
        GTPack8x1 pa{};
        GTPack8x1 pb{};
        GTPack8x1 pr{};

        for (size_t lane = 0; lane < 8; ++lane) {
            make_random_gt(a[lane], rng);
            make_random_gt(b[lane], rng);
            GT::add(expected[lane], a[lane], b[lane]);
        }

        GTPack8x1::from_aos(pa, a.data());
        GTPack8x1::from_aos(pb, b.data());

        GTPack8x1::add(pr, pa, pb);
        GTPack8x1::to_aos(got.data(), pr);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }

        GTPack8x1::add(pr, pa, pb);
        GTPack8x1::to_aos(got.data(), pr);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }
    }
}

static void run_add_aliasing()
{
    InsecureRng rng;

    for (int iter = 0; iter < 1; ++iter) {
        std::array<GT, 8> a{};
        std::array<GT, 8> b{};
        std::array<GT, 8> expected{};
        std::array<GT, 8> got{};
        GTPack8x1 pa{};
        GTPack8x1 pb{};
        GTPack8x1 pa_orig{};
        GTPack8x1 pb_orig{};

        for (size_t lane = 0; lane < 8; ++lane) {
            make_random_gt(a[lane], rng);
            make_random_gt(b[lane], rng);
            GT::add(expected[lane], a[lane], b[lane]);
        }

        GTPack8x1::from_aos(pa, a.data());
        GTPack8x1::from_aos(pb, b.data());
        pa_orig = pa;
        pb_orig = pb;

        GTPack8x1::add(pa, pa, pb);
        GTPack8x1::to_aos(got.data(), pa);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }

        GTPack8x1::add(pb, pa_orig, pb);
        GTPack8x1::to_aos(got.data(), pb);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }
    }
}

static void run_doubling_matches_scalar_and_aliasing()
{
    InsecureRng rng;

    for (int iter = 0; iter < 1; ++iter) {
        std::array<GT, 8> a{};
        std::array<GT, 8> got{};
        std::array<GT, 8> expected{};
        GTPack8x1 pa{};
        GTPack8x1 pinout{};
        GTPack8x1 pr{};

        for (size_t lane = 0; lane < 8; ++lane) {
            make_random_gt(a[lane], rng);
            ASSERT_TRUE(GT::is_cyclotomic(a[lane]));
            GT::doubling(expected[lane], a[lane]);
        }

        GTPack8x1::from_aos(pa, a.data());

        GTPack8x1::doubling(pr, pa);
        GTPack8x1::to_aos(got.data(), pr);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }

        pinout = pa;
        GTPack8x1::doubling(pinout, pinout);
        GTPack8x1::to_aos(got.data(), pinout);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }
    }
}

static void run_scalar_multiplication_matches_scalar()
{
    InsecureRng rng;
    constexpr size_t kLaneCount = 8;
    constexpr size_t kStride = 32;
    const std::array<size_t, 2> nbits_cases = {0, 5};

    for (const size_t nbits : nbits_cases) {
        for (int iter = 0; iter < 1; ++iter) {
            std::array<GT, kLaneCount> base{};
            std::array<GT, kLaneCount> expected{};
            std::array<GT, kLaneCount> got_pack{};
            std::array<byte, kLaneCount * kStride> scalars{};
            GTPack8x1 base_pack{};
            GTPack8x1 out_pack{};

            for (size_t lane = 0; lane < kLaneCount; ++lane) {
                make_random_gt(base[lane], rng);
                make_reduced_scalar_little_endian(&scalars[lane * kStride], kStride, nbits, rng);
                GT::mult(expected[lane], base[lane], &scalars[lane * kStride], nbits);
            }

            GTPack8x1::from_aos(base_pack, base.data());
            GTPack8x1::mult(out_pack, base_pack, scalars.data(), kStride, nbits);
            GTPack8x1::to_aos(got_pack.data(), out_pack);

            for (size_t lane = 0; lane < kLaneCount; ++lane) {
                EXPECT_TRUE(GT::equal(got_pack[lane], expected[lane]));
            }
        }
    }
}

static void run_scalar_multiplication_aliasing()
{
    InsecureRng rng;
    constexpr size_t kLaneCount = 8;
    constexpr size_t kStride = 32;
    constexpr size_t kNBits = 5;

    for (int iter = 0; iter < 1; ++iter) {
        std::array<GT, kLaneCount> base{};
        std::array<GT, kLaneCount> expected{};
        std::array<GT, kLaneCount> got{};
        std::array<byte, kLaneCount * kStride> scalars{};
        GTPack8x1 inout_pack{};

        for (size_t lane = 0; lane < kLaneCount; ++lane) {
            make_random_gt(base[lane], rng);
            make_reduced_scalar_little_endian(&scalars[lane * kStride], kStride, kNBits, rng);
            GT::mult(expected[lane], base[lane], &scalars[lane * kStride], kNBits);
        }

        GTPack8x1::from_aos(inout_pack, base.data());
        GTPack8x1::mult(inout_pack, inout_pack, scalars.data(), kStride, kNBits);
        GTPack8x1::to_aos(got.data(), inout_pack);

        for (size_t lane = 0; lane < kLaneCount; ++lane) {
            EXPECT_TRUE(GT::equal(got[lane], expected[lane]));
        }
    }
}

static void run_scalar_multiplication_glv4_matches_scalar()
{
    InsecureRng rng;
    constexpr size_t kLaneCount = 8;
    constexpr size_t kStride = 32;
    const std::array<size_t, 3> nbits_cases = {129, 200, 256};

    for (const size_t kNBits : nbits_cases) {
        std::array<GT, kLaneCount> base{};
        std::array<GT, kLaneCount> expected{};
        std::array<GT, kLaneCount> got_pack{};
        std::array<byte, kLaneCount * kStride> scalars{};
        GTPack8x1 base_pack{};
        GTPack8x1 out_pack{};

        for (size_t lane = 0; lane < kLaneCount; ++lane) {
            base[lane] = GT::identity();
            expected[lane] = GT::identity();
            for (size_t i = 0; i < kStride; ++i) {
                scalars[lane * kStride + i] = 0;
            }
        }

        make_random_gt(base[0], rng);
        scalars[0] = 3;
        GT::mult(expected[0], base[0], &scalars[0], kNBits);

        GTPack8x1::from_aos(base_pack, base.data());
        GTPack8x1::mult(out_pack, base_pack, scalars.data(), kStride, kNBits);
        GTPack8x1::to_aos(got_pack.data(), out_pack);

        for (size_t lane = 0; lane < kLaneCount; ++lane) {
            EXPECT_TRUE(GT::equal(got_pack[lane], expected[lane]));
        }
    }
}

static void run_span_add_and_doubling_match_scalar()
{
    InsecureRng rng;
    constexpr size_t kMaxTuples = 16;
    constexpr size_t kMaxPacks = 2;
    const std::array<size_t, 2> n_pack_cases = {0, 2};

    for (const size_t n_packs : n_pack_cases) {
        const size_t n_padded = n_packs * 8;

        std::array<GT, kMaxTuples> a{};
        std::array<GT, kMaxTuples> b{};
        std::array<GT, kMaxTuples> expected_add{};
        std::array<GT, kMaxTuples> expected_dbl{};
        std::array<GT, kMaxTuples> got{};
        std::array<GTPack8x1, kMaxPacks> pa{};
        std::array<GTPack8x1, kMaxPacks> pb{};
        std::array<GTPack8x1, kMaxPacks> pout{};

        for (size_t i = 0; i < kMaxTuples; ++i) {
            if (i < n_padded) {
                make_random_gt(a[i], rng);
                make_random_gt(b[i], rng);
                GT::add(expected_add[i], a[i], b[i]);
                GT::doubling(expected_dbl[i], a[i]);
            } else {
                make_random_gt(a[i], rng);
                make_random_gt(b[i], rng);
                expected_add[i] = GT::identity();
                expected_dbl[i] = GT::identity();
            }
        }

        for (size_t p = 0; p < n_packs; ++p) {
            GTPack8x1::from_aos(pa[p], &a[p * 8]);
            GTPack8x1::from_aos(pb[p], &b[p * 8]);
        }

        GTPack8x1::add(
            std::span<GTPack8x1>(pout.data(), n_packs),
            std::span<const GTPack8x1>(pa.data(), n_packs),
            std::span<const GTPack8x1>(pb.data(), n_packs));
        for (size_t p = 0; p < n_packs; ++p) {
            GTPack8x1::to_aos(&got[p * 8], pout[p]);
        }
        for (size_t i = 0; i < n_padded; ++i) {
            EXPECT_TRUE(GT::equal(got[i], expected_add[i]));
        }

        GTPack8x1::doubling(
            std::span<GTPack8x1>(pout.data(), n_packs),
            std::span<const GTPack8x1>(pa.data(), n_packs));
        for (size_t p = 0; p < n_packs; ++p) {
            GTPack8x1::to_aos(&got[p * 8], pout[p]);
        }
        for (size_t i = 0; i < n_padded; ++i) {
            EXPECT_TRUE(GT::equal(got[i], expected_dbl[i]));
        }
    }
}

static void run_span_scalar_multiplication_matches_scalar()
{
    InsecureRng rng;
    constexpr size_t kMaxTuples = 16;
    constexpr size_t kMaxPacks = 2;
    constexpr size_t kStride = 32;
    constexpr size_t kNBits = 5;
    const std::array<size_t, 2> n_pack_cases = {0, 2};

    for (const size_t n_packs : n_pack_cases) {
        const size_t n_padded = n_packs * 8;

        std::array<GT, kMaxTuples> base{};
        std::array<GT, kMaxTuples> expected{};
        std::array<GT, kMaxTuples> got{};
        std::array<GTPack8x1, kMaxPacks> pbase{};
        std::array<GTPack8x1, kMaxPacks> pout{};
        std::array<byte, kMaxTuples * kStride> scalars{};

        for (size_t i = 0; i < kMaxTuples; ++i) {
            make_random_gt(base[i], rng);
            make_reduced_scalar_little_endian(&scalars[i * kStride], kStride, kNBits, rng);
            if (i < n_padded) {
                GT::mult(expected[i], base[i], &scalars[i * kStride], kNBits);
            } else {
                expected[i] = GT::identity();
            }
        }

        for (size_t p = 0; p < n_packs; ++p) {
            GTPack8x1::from_aos(pbase[p], &base[p * 8]);
        }

        GTPack8x1::mult(
            std::span<GTPack8x1>(pout.data(), n_packs),
            std::span<const GTPack8x1>(pbase.data(), n_packs),
            std::span<const byte>(scalars.data(), n_packs * 8 * kStride),
            kStride,
            kNBits);

        for (size_t p = 0; p < n_packs; ++p) {
            GTPack8x1::to_aos(&got[p * 8], pout[p]);
        }
        for (size_t i = 0; i < n_padded; ++i) {
            EXPECT_TRUE(GT::equal(got[i], expected[i]));
        }
    }
}

static void run_span_aliasing()
{
    InsecureRng rng;
    constexpr size_t kNPacks = 2;
    constexpr size_t kStride = 32;
    constexpr size_t kNBits = 5;

    std::array<GT, 16> a{};
    std::array<GT, 16> b{};
    std::array<GT, 16> base{};
    std::array<GT, 16> expected_add{};
    std::array<GT, 16> expected_dbl{};
    std::array<GT, 16> expected_mul{};
    std::array<GT, 16> got{};
    std::array<GTPack8x1, kNPacks> pa{};
    std::array<GTPack8x1, kNPacks> pb{};
    std::array<GTPack8x1, kNPacks> pbase{};
    std::array<byte, 16 * kStride> scalars{};

    for (size_t i = 0; i < 16; ++i) {
        make_random_gt(a[i], rng);
        make_random_gt(b[i], rng);
        make_random_gt(base[i], rng);
        make_reduced_scalar_little_endian(&scalars[i * kStride], kStride, kNBits, rng);
        GT::add(expected_add[i], a[i], b[i]);
        GT::doubling(expected_dbl[i], a[i]);
        GT::mult(expected_mul[i], base[i], &scalars[i * kStride], kNBits);
    }

    for (size_t p = 0; p < kNPacks; ++p) {
        GTPack8x1::from_aos(pa[p], &a[p * 8]);
        GTPack8x1::from_aos(pb[p], &b[p * 8]);
        GTPack8x1::from_aos(pbase[p], &base[p * 8]);
    }

    auto pa_orig = pa;

    GTPack8x1::add(
        std::span<GTPack8x1>(pa.data(), kNPacks),
        std::span<const GTPack8x1>(pa.data(), kNPacks),
        std::span<const GTPack8x1>(pb.data(), kNPacks));
    for (size_t p = 0; p < kNPacks; ++p) {
        GTPack8x1::to_aos(&got[p * 8], pa[p]);
    }
    for (size_t i = 0; i < 16; ++i) {
        EXPECT_TRUE(GT::equal(got[i], expected_add[i]));
    }

    GTPack8x1::doubling(
        std::span<GTPack8x1>(pa_orig.data(), kNPacks),
        std::span<const GTPack8x1>(pa_orig.data(), kNPacks));
    for (size_t p = 0; p < kNPacks; ++p) {
        GTPack8x1::to_aos(&got[p * 8], pa_orig[p]);
    }
    for (size_t i = 0; i < 16; ++i) {
        EXPECT_TRUE(GT::equal(got[i], expected_dbl[i]));
    }

    GTPack8x1::mult(
        std::span<GTPack8x1>(pbase.data(), kNPacks),
        std::span<const GTPack8x1>(pbase.data(), kNPacks),
        std::span<const byte>(scalars.data(), kNPacks * 8 * kStride),
        kStride,
        kNBits);
    for (size_t p = 0; p < kNPacks; ++p) {
        GTPack8x1::to_aos(&got[p * 8], pbase[p]);
    }
    for (size_t i = 0; i < 16; ++i) {
        EXPECT_TRUE(GT::equal(got[i], expected_mul[i]));
    }
}

TEST(GTBatchAvx, FastCoverage)
{
    run_add_matches_scalar();
    run_add_aliasing();
    run_doubling_matches_scalar_and_aliasing();
    run_scalar_multiplication_matches_scalar();
    run_scalar_multiplication_aliasing();
    run_scalar_multiplication_glv4_matches_scalar();
    run_span_add_and_doubling_match_scalar();
    run_span_scalar_multiplication_matches_scalar();
    run_span_aliasing();
}

#else

TEST(GTBatchAvx, Disabled)
{
    SUCCEED();
}

#endif
