#include <gtest/gtest.h>

#include <array>

#include "InsecureRng.h"
#include "blst/fields/Fp2.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include "blst/fields/Fp4_avx.h"

using namespace blst;

namespace
{
using FpVec8 = std::array<Fp, 8>;
using FpVec4 = std::array<Fp, 4>;
using Fp2Vec4 = std::array<Fp2, 4>;
using Fp2Vec2 = std::array<Fp2, 2>;

struct Fp4Elem
{
    Fp2 c0;
    Fp2 c1;
};

using Fp4Vec2 = std::array<Fp4Elem, 2>;

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

void pack_fp_8x1w(fp_8x1w& out, const FpVec8& in)
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

void unpack_fp_8x1w(FpVec8& out, const fp_8x1w& in)
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

void pack_fp_4x2w(fp_4x2w& out, const FpVec4& in)
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

void unpack_fp_4x2w(FpVec4& out, const fp_4x2w& in)
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

void pack_fp2_4x2x1w(fp2_4x2x1w& out, const Fp2Vec4& in)
{
    FpVec8 lanes{};
    for (size_t i = 0; i < 4; ++i) {
        lanes[2 * i] = in[i][0];
        lanes[2 * i + 1] = in[i][1];
    }
    pack_fp_8x1w(out, lanes);
}

void unpack_fp2_4x2x1w(Fp2Vec4& out, const fp2_4x2x1w& in)
{
    FpVec8 lanes{};
    unpack_fp_8x1w(lanes, in);
    for (size_t i = 0; i < 4; ++i) {
        out[i][0] = lanes[2 * i];
        out[i][1] = lanes[2 * i + 1];
    }
}

void pack_fp2_2x2x2w(fp2_2x2x2w& out, const Fp2Vec2& in)
{
    FpVec4 lanes = {in[0][0], in[0][1], in[1][0], in[1][1]};
    pack_fp_4x2w(out, lanes);
}

void unpack_fp2_2x2x2w(Fp2Vec2& out, const fp2_2x2x2w& in)
{
    FpVec4 lanes{};
    unpack_fp_4x2w(lanes, in);
    out[0][0] = lanes[0];
    out[0][1] = lanes[1];
    out[1][0] = lanes[2];
    out[1][1] = lanes[3];
}

void pack_fp4_2x2x2x1w(fp4_2x2x2x1w& out, const Fp4Vec2& in)
{
    Fp2Vec4 lanes = {
        in[0].c0, in[0].c1,
        in[1].c0, in[1].c1,
    };
    pack_fp2_4x2x1w(out, lanes);
}

void unpack_fp4_2x2x2x1w(Fp4Vec2& out, const fp4_2x2x2x1w& in)
{
    Fp2Vec4 lanes{};
    unpack_fp2_4x2x1w(lanes, in);
    out[0].c0 = lanes[0];
    out[0].c1 = lanes[1];
    out[1].c0 = lanes[2];
    out[1].c1 = lanes[3];
}

void pack_fp4_1x2x2x2w(fp4_1x2x2x2w& out, const Fp4Elem& in)
{
    Fp2Vec2 lanes = {in.c0, in.c1};
    pack_fp2_2x2x2w(out, lanes);
}

void unpack_fp4_1x2x2x2w(Fp4Elem& out, const fp4_1x2x2x2w& in)
{
    Fp2Vec2 lanes{};
    unpack_fp2_2x2x2w(lanes, in);
    out.c0 = lanes[0];
    out.c1 = lanes[1];
}

void sqr_fp4_scalar(Fp4Elem& out, const Fp4Elem& in)
{
    Fp2 t0{}, t1{}, t2{};
    Fp2::sqr(t0, in.c0);
    Fp2::sqr(t1, in.c1);
    Fp2::mul_by_u_plus_1(t1, t1);
    Fp2::add(out.c0, t0, t1);
    Fp2::mul(t2, in.c0, in.c1);
    Fp2::add(out.c1, t2, t2);
}

void expect_fp2_eq(const Fp2& got, const Fp2& expected)
{
    EXPECT_TRUE(Fp::equal(got[0], expected[0]));
    EXPECT_TRUE(Fp::equal(got[1], expected[1]));
}

void fill_edge_value(Fp4Elem& v, int idx)
{
    Fp two;
    Fp p_minus_one;
    Fp::add(two, Fp::one(), Fp::one());
    Fp::sub(p_minus_one, Fp::zero(), Fp::one());

    const Fp* fp_cycle[4] = {&Fp::zero(), &Fp::one(), &two, &p_minus_one};

    v.c0[0] = *fp_cycle[(idx + 0) & 3];
    v.c0[1] = *fp_cycle[(idx + 1) & 3];
    v.c1[0] = *fp_cycle[(idx + 2) & 3];
    v.c1[1] = *fp_cycle[(idx + 3) & 3];
}
} // namespace

TEST(Fp4Avx, Sqr2x2x2x1wVariants)
{
    InsecureRng rng;

    for (int iter = 0; iter < 192; ++iter) {
        Fp4Vec2 a{};
        if (iter < 16) {
            fill_edge_value(a[0], iter);
            fill_edge_value(a[1], iter + 1);
        } else {
            rand_redc_fp2(a[0].c0, rng);
            rand_redc_fp2(a[0].c1, rng);
            rand_redc_fp2(a[1].c0, rng);
            rand_redc_fp2(a[1].c1, rng);
        }

        fp4_2x2x2x1w avx_a{}, r1{}, r2{}, r3{};
        Fp4Vec2 got1{}, got2{}, got3{}, expected{};

        pack_fp4_2x2x2x1w(avx_a, a);
        sqr_fp4_2x2x2x1w_v1(r1, avx_a);
        sqr_fp4_2x2x2x1w_v2(r2, avx_a);
        sqr_fp4_2x2x2x1w_v3(r3, avx_a);

        unpack_fp4_2x2x2x1w(got1, r1);
        unpack_fp4_2x2x2x1w(got2, r2);
        unpack_fp4_2x2x2x1w(got3, r3);

        sqr_fp4_scalar(expected[0], a[0]);
        sqr_fp4_scalar(expected[1], a[1]);

        for (size_t i = 0; i < 2; ++i) {
            expect_fp2_eq(got1[i].c0, expected[i].c0);
            expect_fp2_eq(got1[i].c1, expected[i].c1);
            expect_fp2_eq(got2[i].c0, expected[i].c0);
            expect_fp2_eq(got2[i].c1, expected[i].c1);
            expect_fp2_eq(got3[i].c0, expected[i].c0);
            expect_fp2_eq(got3[i].c1, expected[i].c1);

            expect_fp2_eq(got1[i].c0, got2[i].c0);
            expect_fp2_eq(got1[i].c1, got2[i].c1);
            expect_fp2_eq(got1[i].c0, got3[i].c0);
            expect_fp2_eq(got1[i].c1, got3[i].c1);
        }
    }
}

TEST(Fp4Avx, Sqr1x2x2x2wVariants)
{
    InsecureRng rng;

    for (int iter = 0; iter < 192; ++iter) {
        Fp4Elem a{};
        if (iter < 16) {
            fill_edge_value(a, iter);
        } else {
            rand_redc_fp2(a.c0, rng);
            rand_redc_fp2(a.c1, rng);
        }

        fp4_1x2x2x2w avx_a{}, r1{}, r2{};
        Fp4Elem got1{}, got2{}, expected{};

        pack_fp4_1x2x2x2w(avx_a, a);
        sqr_fp4_1x2x2x2w_v1(r1, avx_a);
        sqr_fp4_1x2x2x2w_v2(r2, avx_a);

        unpack_fp4_1x2x2x2w(got1, r1);
        unpack_fp4_1x2x2x2w(got2, r2);
        sqr_fp4_scalar(expected, a);

        expect_fp2_eq(got1.c0, expected.c0);
        expect_fp2_eq(got1.c1, expected.c1);
        expect_fp2_eq(got2.c0, expected.c0);
        expect_fp2_eq(got2.c1, expected.c1);
        expect_fp2_eq(got1.c0, got2.c0);
        expect_fp2_eq(got1.c1, got2.c1);
    }
}

#else

TEST(Fp4Avx, Disabled)
{
    SUCCEED();
}

#endif
