#include <gtest/gtest.h>

#include <array>

#include "InsecureRng.h"
#include "blst/fields/Fp2.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include "blst/fields/Fp2_avx.h"

using namespace blst;

namespace
{
using FpVec8 = std::array<Fp, 8>;
using FpVec4 = std::array<Fp, 4>;
using Fp2Vec8 = std::array<Fp2, 8>;
using Fp2Vec4 = std::array<Fp2, 4>;
using Fp2Vec2 = std::array<Fp2, 2>;

enum class ElemOp
{
    Add,
    Sub
};

constexpr std::array<ElemOp, 4> kAass4 = {
    ElemOp::Sub, ElemOp::Sub, ElemOp::Add, ElemOp::Add,
};

constexpr std::array<ElemOp, 4> kAssa4 = {
    ElemOp::Add, ElemOp::Sub, ElemOp::Sub, ElemOp::Add,
};

constexpr std::array<ElemOp, 4> kSass4 = {
    ElemOp::Sub, ElemOp::Sub, ElemOp::Add, ElemOp::Sub,
};

constexpr std::array<ElemOp, 4> kSsaa4 = {
    ElemOp::Add, ElemOp::Add, ElemOp::Sub, ElemOp::Sub,
};

constexpr std::array<ElemOp, 2> kAs2 = {
    ElemOp::Sub, ElemOp::Add,
};

constexpr std::array<ElemOp, 2> kSa2 = {
    ElemOp::Add, ElemOp::Sub,
};

constexpr std::array<ElemOp, 8> kAsx4_8 = {
    ElemOp::Sub, ElemOp::Add, ElemOp::Sub, ElemOp::Add,
    ElemOp::Sub, ElemOp::Add, ElemOp::Sub, ElemOp::Add,
};

constexpr std::array<ElemOp, 4> kSax2_4 = {
    ElemOp::Add, ElemOp::Sub, ElemOp::Add, ElemOp::Sub,
};

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

void pack_FpPack8x1(FpPack8x1& out, const FpVec8& in)
{
    std::array<std::array<u64, NWORDS>, 8> words{};
    for (size_t lane = 0; lane < 8; ++lane) {
        conv_64to48_fp(words[lane].data(), reinterpret_cast<const u64*>(in[lane].data));
    }

    for (size_t i = 0; i < NWORDS; ++i) {
        out[i] = _cc512_set_epi64(
            words[7][i], words[6][i], words[5][i], words[4][i],
            words[3][i], words[2][i], words[1][i], words[0][i]);
    }
}

void unpack_FpPack8x1(FpVec8& out, const FpPack8x1& in)
{
    std::array<std::array<u64, NWORDS>, 8> words{};
    std::array<u64, 8> lanes{};
    for (size_t i = 0; i < NWORDS; ++i) {
        _cc512_storeu_si512(lanes.data(), in[i]);
        for (size_t lane = 0; lane < 8; ++lane) {
            words[lane][i] = lanes[lane];
        }
    }
    for (size_t lane = 0; lane < 8; ++lane) {
        conv_48to64_fp(reinterpret_cast<u64*>(out[lane].data), words[lane].data());
    }
}

void pack_FpPack4x2(FpPack4x2& out, const FpVec4& in)
{
    std::array<std::array<u64, NWORDS>, 4> words{};
    for (size_t lane = 0; lane < 4; ++lane) {
        conv_64to48_fp(words[lane].data(), reinterpret_cast<const u64*>(in[lane].data));
    }

    for (size_t i = 0; i < VWORDS; ++i) {
        out[i] = _cc512_set_epi64(
            words[3][i + VWORDS], words[3][i], words[2][i + VWORDS], words[2][i],
            words[1][i + VWORDS], words[1][i], words[0][i + VWORDS], words[0][i]);
    }
}

void unpack_FpPack4x2(FpVec4& out, const FpPack4x2& in)
{
    std::array<std::array<u64, NWORDS>, 4> words{};
    std::array<u64, 8> lanes{};
    for (size_t i = 0; i < VWORDS; ++i) {
        _cc512_storeu_si512(lanes.data(), in[i]);
        for (size_t lane = 0; lane < 4; ++lane) {
            words[lane][i] = lanes[2 * lane];
            words[lane][i + VWORDS] = lanes[2 * lane + 1];
        }
    }
    for (size_t lane = 0; lane < 4; ++lane) {
        conv_48to64_fp(reinterpret_cast<u64*>(out[lane].data), words[lane].data());
    }
}

void pack_fp2_8x1x1w(fp2_8x1x1w& out, const Fp2Vec8& in)
{
    FpVec8 c0{};
    FpVec8 c1{};
    for (size_t i = 0; i < 8; ++i) {
        c0[i] = in[i][0];
        c1[i] = in[i][1];
    }
    pack_FpPack8x1(out[0], c0);
    pack_FpPack8x1(out[1], c1);
}

void unpack_fp2_8x1x1w(Fp2Vec8& out, const fp2_8x1x1w& in)
{
    FpVec8 c0{};
    FpVec8 c1{};
    unpack_FpPack8x1(c0, in[0]);
    unpack_FpPack8x1(c1, in[1]);
    for (size_t i = 0; i < 8; ++i) {
        out[i][0] = c0[i];
        out[i][1] = c1[i];
    }
}

void pack_fp2_4x2x1w(fp2_4x2x1w& out, const Fp2Vec4& in)
{
    FpVec8 lanes{};
    for (size_t i = 0; i < 4; ++i) {
        lanes[2 * i] = in[i][0];
        lanes[2 * i + 1] = in[i][1];
    }
    pack_FpPack8x1(out, lanes);
}

void unpack_fp2_4x2x1w(Fp2Vec4& out, const fp2_4x2x1w& in)
{
    FpVec8 lanes{};
    unpack_FpPack8x1(lanes, in);
    for (size_t i = 0; i < 4; ++i) {
        out[i][0] = lanes[2 * i];
        out[i][1] = lanes[2 * i + 1];
    }
}

void pack_fp2_2x2x2w(fp2_2x2x2w& out, const Fp2Vec2& in)
{
    FpVec4 lanes = {in[0][0], in[0][1], in[1][0], in[1][1]};
    pack_FpPack4x2(out, lanes);
}

void unpack_fp2_2x2x2w(Fp2Vec2& out, const fp2_2x2x2w& in)
{
    FpVec4 lanes{};
    unpack_FpPack4x2(lanes, in);
    out[0][0] = lanes[0];
    out[0][1] = lanes[1];
    out[1][0] = lanes[2];
    out[1][1] = lanes[3];
}

void redc_fp2x2_8x1x1w_local(fp2_8x1x1w& r, const fp2x2_8x1x1w& a)
{
    redc_UFpPack8x1(r[0], a[0]);
    redc_UFpPack8x1(r[1], a[1]);
}

void apply_pattern(Fp2Vec4& out, const Fp2Vec4& a, const Fp2Vec4& b, const std::array<ElemOp, 4>& pattern)
{
    for (size_t i = 0; i < 4; ++i) {
        if (pattern[i] == ElemOp::Add) {
            Fp2::add(out[i], a[i], b[i]);
        } else {
            Fp2::sub(out[i], a[i], b[i]);
        }
    }
}

void apply_pattern(Fp2Vec2& out, const Fp2Vec2& a, const Fp2Vec2& b, const std::array<ElemOp, 2>& pattern)
{
    for (size_t i = 0; i < 2; ++i) {
        if (pattern[i] == ElemOp::Add) {
            Fp2::add(out[i], a[i], b[i]);
        } else {
            Fp2::sub(out[i], a[i], b[i]);
        }
    }
}

void apply_pattern(Fp2Vec8& out, const Fp2Vec8& a, const Fp2Vec8& b, const std::array<ElemOp, 8>& pattern)
{
    for (size_t i = 0; i < 8; ++i) {
        if (pattern[i] == ElemOp::Add) {
            Fp2::add(out[i], a[i], b[i]);
        } else {
            Fp2::sub(out[i], a[i], b[i]);
        }
    }
}

void expect_vec_words_eq(const cc512i* got, const cc512i* expected, size_t words)
{
    for (size_t i = 0; i < words; ++i) {
        std::array<u64, 8> got_lanes{};
        std::array<u64, 8> exp_lanes{};
        _cc512_storeu_si512(got_lanes.data(), got[i]);
        _cc512_storeu_si512(exp_lanes.data(), expected[i]);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_EQ(got_lanes[lane], exp_lanes[lane]);
        }
    }
}
} // namespace

TEST(Fp2Avx, Add8x1x1w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 128; ++iter) {
        Fp2Vec8 a{}, b{}, got{}, expected{};
        for (size_t i = 0; i < 8; ++i) {
            rand_redc_fp2(a[i], rng);
            rand_redc_fp2(b[i], rng);
            Fp2::add(expected[i], a[i], b[i]);
        }

        fp2_8x1x1w avx_a{}, avx_b{}, avx_r{};
        pack_fp2_8x1x1w(avx_a, a);
        pack_fp2_8x1x1w(avx_b, b);
        add_fp2_8x1x1w(avx_r, avx_a, avx_b);
        unpack_fp2_8x1x1w(got, avx_r);

        for (size_t i = 0; i < 8; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }
    }
}

TEST(Fp2Avx, Core4x2x1w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 192; ++iter) {
        Fp2Vec4 a{}, b{}, got{}, expected{};
        for (size_t i = 0; i < 4; ++i) {
            rand_redc_fp2(a[i], rng);
            rand_redc_fp2(b[i], rng);
        }

        fp2_4x2x1w avx_a{}, avx_b{}, avx_r{};
        pack_fp2_4x2x1w(avx_a, a);
        pack_fp2_4x2x1w(avx_b, b);

        add_fp2_4x2x1w(avx_r, avx_a, avx_b);
        unpack_fp2_4x2x1w(got, avx_r);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::add(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sub_fp2_4x2x1w(avx_r, avx_a, avx_b);
        unpack_fp2_4x2x1w(got, avx_r);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::sub(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        aass_fp2_4x2x1w(avx_r, avx_a, avx_b);
        unpack_fp2_4x2x1w(got, avx_r);
        apply_pattern(expected, a, b, kAass4);
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        assa_fp2_4x2x1w(avx_r, avx_a, avx_b);
        unpack_fp2_4x2x1w(got, avx_r);
        apply_pattern(expected, a, b, kAssa4);
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sass_fp2_4x2x1w(avx_r, avx_a, avx_b);
        unpack_fp2_4x2x1w(got, avx_r);
        apply_pattern(expected, a, b, kSass4);
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        ssaa_fp2_4x2x1w(avx_r, avx_a, avx_b);
        unpack_fp2_4x2x1w(got, avx_r);
        apply_pattern(expected, a, b, kSsaa4);
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        mul_by_u_plus_1_fp2_4x2x1w(avx_r, avx_a);
        unpack_fp2_4x2x1w(got, avx_r);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::mul_by_u_plus_1(expected[i], a[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        mul_fp2_4x2x1w(avx_r, avx_a, avx_b);
        unpack_fp2_4x2x1w(got, avx_r);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::mul(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sqr_fp2_4x2x1w(avx_r, avx_a);
        unpack_fp2_4x2x1w(got, avx_r);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::sqr(expected[i], a[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }
    }
}

TEST(Fp2Avx, Core2x2x2w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 192; ++iter) {
        Fp2Vec2 a{}, b{}, got{}, expected{};
        rand_redc_fp2(a[0], rng);
        rand_redc_fp2(a[1], rng);
        rand_redc_fp2(b[0], rng);
        rand_redc_fp2(b[1], rng);

        fp2_2x2x2w avx_a{}, avx_b{}, avx_r{};
        pack_fp2_2x2x2w(avx_a, a);
        pack_fp2_2x2x2w(avx_b, b);

        add_fp2_2x2x2w(avx_r, avx_a, avx_b);
        unpack_fp2_2x2x2w(got, avx_r);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::add(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sub_fp2_2x2x2w(avx_r, avx_a, avx_b);
        unpack_fp2_2x2x2w(got, avx_r);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::sub(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        as_fp2_2x2x2w(avx_r, avx_a, avx_b);
        unpack_fp2_2x2x2w(got, avx_r);
        apply_pattern(expected, a, b, kAs2);
        for (size_t i = 0; i < 2; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sa_fp2_2x2x2w(avx_r, avx_a, avx_b);
        unpack_fp2_2x2x2w(got, avx_r);
        apply_pattern(expected, a, b, kSa2);
        for (size_t i = 0; i < 2; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        mul_by_u_plus_1_fp2_2x2x2w(avx_r, avx_a);
        unpack_fp2_2x2x2w(got, avx_r);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::mul_by_u_plus_1(expected[i], a[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        mul_fp2_2x2x2w(avx_r, avx_a, avx_b);
        unpack_fp2_2x2x2w(got, avx_r);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::mul(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sqr_fp2_2x2x2w(avx_r, avx_a);
        unpack_fp2_2x2x2w(got, avx_r);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::sqr(expected[i], a[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }
    }
}

TEST(Fp2Avx, Widened8x1x1w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 96; ++iter) {
        Fp2Vec8 a{}, b{}, c{}, d{}, ab{}, cd{}, got{}, expected{};
        for (size_t i = 0; i < 8; ++i) {
            rand_redc_fp2(a[i], rng);
            rand_redc_fp2(b[i], rng);
            rand_redc_fp2(c[i], rng);
            rand_redc_fp2(d[i], rng);
            Fp2::mul(ab[i], a[i], b[i]);
            Fp2::mul(cd[i], c[i], d[i]);
        }

        fp2_8x1x1w avx_a{}, avx_b{}, avx_c{}, avx_d{}, red{};
        fp2x2_8x1x1w wide_ab{}, wide_cd{}, wide_r{};
        pack_fp2_8x1x1w(avx_a, a);
        pack_fp2_8x1x1w(avx_b, b);
        pack_fp2_8x1x1w(avx_c, c);
        pack_fp2_8x1x1w(avx_d, d);

        mul_fp2x2_8x1x1w(wide_ab, avx_a, avx_b);
        mul_fp2x2_8x1x1w(wide_cd, avx_c, avx_d);

        add_fp2x2_8x1x1w(wide_r, wide_ab, wide_cd);
        redc_fp2x2_8x1x1w_local(red, wide_r);
        unpack_fp2_8x1x1w(got, red);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::add(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sub_fp2x2_8x1x1w(wide_r, wide_ab, wide_cd);
        redc_fp2x2_8x1x1w_local(red, wide_r);
        unpack_fp2_8x1x1w(got, red);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::sub(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        asx4_fp2x2_8x1x1w(wide_r, wide_ab, wide_cd);
        redc_fp2x2_8x1x1w_local(red, wide_r);
        unpack_fp2_8x1x1w(got, red);
        apply_pattern(expected, ab, cd, kAsx4_8);
        for (size_t i = 0; i < 8; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        mul_by_u_plus_1_fp2x2_8x1x1w(wide_r, wide_ab);
        redc_fp2x2_8x1x1w_local(red, wide_r);
        unpack_fp2_8x1x1w(got, red);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::mul_by_u_plus_1(expected[i], ab[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }
    }
}

TEST(Fp2Avx, Widened4x2x1w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 96; ++iter) {
        Fp2Vec4 a{}, b{}, c{}, d{}, ab{}, cd{}, got{}, expected{};
        for (size_t i = 0; i < 4; ++i) {
            rand_redc_fp2(a[i], rng);
            rand_redc_fp2(b[i], rng);
            rand_redc_fp2(c[i], rng);
            rand_redc_fp2(d[i], rng);
            Fp2::mul(ab[i], a[i], b[i]);
            Fp2::mul(cd[i], c[i], d[i]);
        }

        fp2_4x2x1w avx_a{}, avx_b{}, avx_c{}, avx_d{}, red{};
        fp2x2_4x2x1w wide_ab{}, wide_cd{}, wide_r{};
        pack_fp2_4x2x1w(avx_a, a);
        pack_fp2_4x2x1w(avx_b, b);
        pack_fp2_4x2x1w(avx_c, c);
        pack_fp2_4x2x1w(avx_d, d);

        mul_fp2x2_4x2x1w(wide_ab, avx_a, avx_b);
        mul_fp2x2_4x2x1w(wide_cd, avx_c, avx_d);

        add_fp2x2_4x2x1w(wide_r, wide_ab, wide_cd);
        redc_fp2x2_4x2x1w(red, wide_r);
        unpack_fp2_4x2x1w(got, red);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::add(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sub_fp2x2_4x2x1w(wide_r, wide_ab, wide_cd);
        redc_fp2x2_4x2x1w(red, wide_r);
        unpack_fp2_4x2x1w(got, red);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::sub(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sax2_fp2x2_4x2x1w(wide_r, wide_ab, wide_cd);
        redc_fp2x2_4x2x1w(red, wide_r);
        unpack_fp2_4x2x1w(got, red);
        apply_pattern(expected, ab, cd, kSax2_4);
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        mul_by_u_plus_1_fp2x2_4x2x1w(wide_r, wide_ab);
        redc_fp2x2_4x2x1w(red, wide_r);
        unpack_fp2_4x2x1w(got, red);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::mul_by_u_plus_1(expected[i], ab[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sqr_fp2x2_4x2x1w(wide_r, avx_a);
        redc_fp2x2_4x2x1w(red, wide_r);
        unpack_fp2_4x2x1w(got, red);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::sqr(expected[i], a[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }
    }
}

TEST(Fp2Avx, Widened2x2x2w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 96; ++iter) {
        Fp2Vec2 a{}, b{}, c{}, d{}, ab{}, cd{}, got{}, expected{};
        rand_redc_fp2(a[0], rng);
        rand_redc_fp2(a[1], rng);
        rand_redc_fp2(b[0], rng);
        rand_redc_fp2(b[1], rng);
        rand_redc_fp2(c[0], rng);
        rand_redc_fp2(c[1], rng);
        rand_redc_fp2(d[0], rng);
        rand_redc_fp2(d[1], rng);

        Fp2::mul(ab[0], a[0], b[0]);
        Fp2::mul(ab[1], a[1], b[1]);
        Fp2::mul(cd[0], c[0], d[0]);
        Fp2::mul(cd[1], c[1], d[1]);

        fp2_2x2x2w avx_a{}, avx_b{}, avx_c{}, avx_d{}, red{};
        fp2x2_2x2x2w wide_ab{}, wide_cd{}, wide_r{};
        pack_fp2_2x2x2w(avx_a, a);
        pack_fp2_2x2x2w(avx_b, b);
        pack_fp2_2x2x2w(avx_c, c);
        pack_fp2_2x2x2w(avx_d, d);

        mul_fp2x2_2x2x2w(wide_ab, avx_a, avx_b);
        mul_fp2x2_2x2x2w(wide_cd, avx_c, avx_d);

        add_fp2x2_2x2x2w(wide_r, wide_ab, wide_cd);
        redc_fp2x2_2x2x2w(red, wide_r);
        unpack_fp2_2x2x2w(got, red);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::add(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sub_fp2x2_2x2x2w(wide_r, wide_ab, wide_cd);
        redc_fp2x2_2x2x2w(red, wide_r);
        unpack_fp2_2x2x2w(got, red);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::sub(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        mul_by_u_plus_1_fp2x2_2x2x2w(wide_r, wide_ab);
        redc_fp2x2_2x2x2w(red, wide_r);
        unpack_fp2_2x2x2w(got, red);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::mul_by_u_plus_1(expected[i], ab[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        sqr_fp2x2_2x2x2w(wide_r, avx_a);
        redc_fp2x2_2x2x2w(red, wide_r);
        unpack_fp2_2x2x2w(got, red);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::sqr(expected[i], a[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }
    }
}

TEST(Fp2Avx, Schedule2x4x1w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 192; ++iter) {
        Fp2Vec4 a{}, b{};
        for (size_t i = 0; i < 4; ++i) {
            rand_redc_fp2(a[i], rng);
            rand_redc_fp2(b[i], rng);
        }

        fp2_2x4x1w avx_a{}, avx_b{}, red_direct{}, red_from_wide{};
        fp2x2_2x4x1w wide{};
        pack_fp2_4x2x1w(avx_a, a);
        pack_fp2_4x2x1w(avx_b, b);

        mul_fp2_2x4x1w(red_direct, avx_a, avx_b);
        mul_fp2x2_2x4x1w(wide, avx_a, avx_b);
        redc_fp2x2_4x2x1w(red_from_wide, wide);

        expect_vec_words_eq(red_direct.data(), red_from_wide.data(), NWORDS);
    }
}

TEST(Fp2Avx, Schedule2x4x1wScalar)
{
    InsecureRng rng;
    for (int iter = 0; iter < 128; ++iter) {
        Fp2 a{}, b{}, expected{};
        rand_redc_fp2(a, rng);
        rand_redc_fp2(b, rng);
        Fp2::mul(expected, a, b);

        Fp2Vec4 lanes_a = {a, a, a, a};
        Fp2Vec4 lanes_b = {b, b, b, b};
        Fp2Vec4 got{};
        fp2_2x4x1w avx_a{}, avx_b{}, avx_r{};

        pack_fp2_4x2x1w(avx_a, lanes_a);
        pack_fp2_4x2x1w(avx_b, lanes_b);
        mul_fp2_2x4x1w(avx_r, avx_a, avx_b);
        unpack_fp2_4x2x1w(got, avx_r);

        bool found_expected = false;
        for (size_t i = 0; i < 4; ++i) {
            if (Fp::equal(got[i][0], expected[0]) && Fp::equal(got[i][1], expected[1])) {
                found_expected = true;
                break;
            }
        }
        EXPECT_TRUE(found_expected);
    }
}

TEST(Fp2Avx, Schedule1x4x2w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 192; ++iter) {
        Fp2Vec2 a{}, b{};
        rand_redc_fp2(a[0], rng);
        rand_redc_fp2(a[1], rng);
        rand_redc_fp2(b[0], rng);
        rand_redc_fp2(b[1], rng);

        fp2_1x4x2w avx_a{}, avx_b{}, red_direct{}, red_from_wide{};
        fp2x2_1x4x2w wide{};
        pack_fp2_2x2x2w(avx_a, a);
        pack_fp2_2x2x2w(avx_b, b);

        mul_fp2_1x4x2w(red_direct, avx_a, avx_b);
        mul_fp2x2_1x4x2w(wide, avx_a, avx_b);
        redc_fp2x2_2x2x2w(red_from_wide, wide);

        expect_vec_words_eq(red_direct.data(), red_from_wide.data(), VWORDS);
    }
}

TEST(Fp2Avx, Schedule1x4x2wScalar)
{
    InsecureRng rng;
    for (int iter = 0; iter < 128; ++iter) {
        Fp2 a{}, b{}, expected{};
        rand_redc_fp2(a, rng);
        rand_redc_fp2(b, rng);
        Fp2::mul(expected, a, b);

        Fp2Vec2 lanes_a = {a, a};
        Fp2Vec2 lanes_b = {b, b};
        Fp2Vec2 got{};
        fp2_1x4x2w avx_a{}, avx_b{}, avx_r{};

        pack_fp2_2x2x2w(avx_a, lanes_a);
        pack_fp2_2x2x2w(avx_b, lanes_b);
        mul_fp2_1x4x2w(avx_r, avx_a, avx_b);
        unpack_fp2_2x2x2w(got, avx_r);

        bool found_expected = false;
        for (size_t i = 0; i < 2; ++i) {
            if (Fp::equal(got[i][0], expected[0]) && Fp::equal(got[i][1], expected[1])) {
                found_expected = true;
                break;
            }
        }
        EXPECT_TRUE(found_expected);
    }
}

#else

TEST(Fp2Avx, Disabled)
{
    SUCCEED();
}

#endif
