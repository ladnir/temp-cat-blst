#include <gtest/gtest.h>

#include <array>
#include <cstring>

#include "blst/fields/Fp.h"
#include "InsecureRng.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include "blst/fields/Fp_avx.h"

using namespace blst;

namespace
{
static_assert(BLST_LIMB_SIZE == 64, "AVX Fp tests require 64-bit limbs.");

using FpVec8 = std::array<Fp, 8>;
using FpVec4 = std::array<Fp, 4>;

enum class LaneOp
{
    Add,
    Sub
};

constexpr std::array<LaneOp, 8> kAsx4Pattern = {
    LaneOp::Sub, LaneOp::Add, LaneOp::Sub, LaneOp::Add,
    LaneOp::Sub, LaneOp::Add, LaneOp::Sub, LaneOp::Add,
};

constexpr std::array<LaneOp, 8> kAsax6Pattern = {
    LaneOp::Add, LaneOp::Add, LaneOp::Add, LaneOp::Add,
    LaneOp::Add, LaneOp::Add, LaneOp::Sub, LaneOp::Add,
};

constexpr std::array<LaneOp, 8> kAssx6Pattern = {
    LaneOp::Sub, LaneOp::Sub, LaneOp::Sub, LaneOp::Sub,
    LaneOp::Sub, LaneOp::Sub, LaneOp::Sub, LaneOp::Add,
};

constexpr std::array<LaneOp, 4> kAsx2Pattern = {
    LaneOp::Sub, LaneOp::Add, LaneOp::Sub, LaneOp::Add,
};

auto load_lanes(const cc512i& v) -> std::array<u64, 8>
{
    std::array<u64, 8> out{};
    ::blst::_cc512_storeu_si512(out.data(), v);
    return out;
}

auto make_vec(const std::array<u64, 8>& lanes) -> cc512i
{
    return ::blst::_cc512_loadu_si512(lanes.data());
}

void rand_redc_fp(Fp& v, InsecureRng& rng)
{
    for (u64 i = 0; i < v.num_limb(); ++i) {
        v[i] = rng.next_u64();
    }

    while (Fp::is_redc(v) == 0) {
        Fp::rshift(v, v, 1);
    }
}

void force_edge_lanes(FpVec8& v)
{
    Fp two;
    Fp p_minus_one;
    Fp::add(two, Fp::one(), Fp::one());
    Fp::sub(p_minus_one, Fp::zero(), Fp::one());

    v[0] = Fp::zero();
    v[1] = Fp::one();
    v[2] = p_minus_one;
    v[3] = two;
}

void force_edge_lanes(FpVec4& v)
{
    Fp two;
    Fp p_minus_one;
    Fp::add(two, Fp::one(), Fp::one());
    Fp::sub(p_minus_one, Fp::zero(), Fp::one());

    v[0] = Fp::zero();
    v[1] = Fp::one();
    v[2] = p_minus_one;
    v[3] = two;
}

void pack_fp_8x1w(fp_8x1w out, const FpVec8& in)
{
    std::array<std::array<u64, NWORDS>, 8> words{};
    for (size_t lane = 0; lane < 8; ++lane) {
        conv_64to48_fp(words[lane].data(), reinterpret_cast<const u64*>(in[lane].data));
    }

    for (size_t i = 0; i < NWORDS; ++i) {
        std::array<u64, 8> lanes{};
        for (size_t lane = 0; lane < 8; ++lane) {
            lanes[lane] = words[lane][i];
        }
        out[i] = make_vec(lanes);
    }
}

void unpack_fp_8x1w(FpVec8& out, const fp_8x1w in)
{
    std::array<std::array<u64, NWORDS>, 8> words{};
    for (size_t i = 0; i < NWORDS; ++i) {
        const auto lanes = load_lanes(in[i]);
        for (size_t lane = 0; lane < 8; ++lane) {
            words[lane][i] = lanes[lane];
        }
    }

    for (size_t lane = 0; lane < 8; ++lane) {
        conv_48to64_fp(reinterpret_cast<u64*>(out[lane].data), words[lane].data());
    }
}

void pack_fp_4x2w(fp_4x2w out, const FpVec4& in)
{
    std::array<std::array<u64, NWORDS>, 4> words{};
    for (size_t lane = 0; lane < 4; ++lane) {
        conv_64to48_fp(words[lane].data(), reinterpret_cast<const u64*>(in[lane].data));
    }

    for (size_t i = 0; i < VWORDS; ++i) {
        std::array<u64, 8> lanes{};
        for (size_t lane = 0; lane < 4; ++lane) {
            lanes[2 * lane] = words[lane][i];
            lanes[2 * lane + 1] = words[lane][i + VWORDS];
        }
        out[i] = make_vec(lanes);
    }
}

void unpack_fp_4x2w(FpVec4& out, const fp_4x2w in)
{
    std::array<std::array<u64, NWORDS>, 4> words{};
    for (size_t i = 0; i < VWORDS; ++i) {
        const auto lanes = load_lanes(in[i]);
        for (size_t lane = 0; lane < 4; ++lane) {
            words[lane][i] = lanes[2 * lane];
            words[lane][i + VWORDS] = lanes[2 * lane + 1];
        }
    }

    for (size_t lane = 0; lane < 4; ++lane) {
        conv_48to64_fp(reinterpret_cast<u64*>(out[lane].data), words[lane].data());
    }
}

void apply_pattern(FpVec8& out, const FpVec8& a, const FpVec8& b, const std::array<LaneOp, 8>& pattern)
{
    for (size_t i = 0; i < 8; ++i) {
        if (pattern[i] == LaneOp::Add) {
            Fp::add(out[i], a[i], b[i]);
        } else {
            Fp::sub(out[i], a[i], b[i]);
        }
    }
}

void apply_pattern(FpVec4& out, const FpVec4& a, const FpVec4& b, const std::array<LaneOp, 4>& pattern)
{
    for (size_t i = 0; i < 4; ++i) {
        if (pattern[i] == LaneOp::Add) {
            Fp::add(out[i], a[i], b[i]);
        } else {
            Fp::sub(out[i], a[i], b[i]);
        }
    }
}

} // namespace

TEST(FpAvx, ConversionRoundTrip)
{
    InsecureRng rng;

    for (int iter = 0; iter < 128; ++iter) {
        FpVec8 src8{};
        for (auto& e : src8) {
            rand_redc_fp(e, rng);
        }
        force_edge_lanes(src8);

        cc512i in64_8[6];
        for (size_t limb = 0; limb < 6; ++limb) {
            std::array<u64, 8> lanes{};
            for (size_t lane = 0; lane < 8; ++lane) {
                lanes[lane] = static_cast<u64>(src8[lane][limb]);
            }
            in64_8[limb] = make_vec(lanes);
        }

        fp_8x1w fp8{};
        cc512i out64_8[6];
        conv_64to48_fp_8x1w(fp8, in64_8);
        conv_48to64_fp_8x1w(out64_8, fp8);

        for (size_t limb = 0; limb < 6; ++limb) {
            const auto in_lanes = load_lanes(in64_8[limb]);
            const auto out_lanes = load_lanes(out64_8[limb]);
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_EQ(out_lanes[lane], in_lanes[lane]);
            }
        }

        FpVec4 src4{};
        for (size_t i = 0; i < 4; ++i) {
            src4[i] = src8[i];
        }
        force_edge_lanes(src4);

        cc512i in64_4[3];
        for (size_t limb = 0; limb < 3; ++limb) {
            std::array<u64, 8> lanes{};
            for (size_t lane = 0; lane < 4; ++lane) {
                lanes[2 * lane] = static_cast<u64>(src4[lane][limb]);
                lanes[2 * lane + 1] = static_cast<u64>(src4[lane][limb + 3]);
            }
            in64_4[limb] = make_vec(lanes);
        }

        fp_4x2w fp4{};
        cc512i out64_4[3];
        conv_64to48_fp_4x2w(fp4, in64_4);
        conv_48to64_fp_4x2w(out64_4, fp4);

        for (size_t limb = 0; limb < 3; ++limb) {
            const auto in_lanes = load_lanes(in64_4[limb]);
            const auto out_lanes = load_lanes(out64_4[limb]);
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_EQ(out_lanes[lane], in_lanes[lane]);
            }
        }
    }
}

TEST(FpAvx, Fp8x1wOps)
{
    InsecureRng rng;

    for (int iter = 0; iter < 256; ++iter) {
        FpVec8 a{}, b{};
        for (size_t i = 0; i < 8; ++i) {
            rand_redc_fp(a[i], rng);
            rand_redc_fp(b[i], rng);
        }
        force_edge_lanes(a);
        force_edge_lanes(b);

        fp_8x1w avx_a{}, avx_b{}, avx_r{};
        pack_fp_8x1w(avx_a, a);
        pack_fp_8x1w(avx_b, b);

        FpVec8 got{}, expected{};

        add_fp_8x1w(avx_r, avx_a, avx_b);
        unpack_fp_8x1w(got, avx_r);
        for (size_t lane = 0; lane < 8; ++lane) {
            Fp::add(expected[lane], a[lane], b[lane]);
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        sub_fp_8x1w(avx_r, avx_a, avx_b);
        unpack_fp_8x1w(got, avx_r);
        for (size_t lane = 0; lane < 8; ++lane) {
            Fp::sub(expected[lane], a[lane], b[lane]);
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        asx4_fp_8x1w(avx_r, avx_a, avx_b);
        unpack_fp_8x1w(got, avx_r);
        apply_pattern(expected, a, b, kAsx4Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        asax6_fp_8x1w(avx_r, avx_a, avx_b);
        unpack_fp_8x1w(got, avx_r);
        apply_pattern(expected, a, b, kAsax6Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        assx6_fp_8x1w(avx_r, avx_a, avx_b);
        unpack_fp_8x1w(got, avx_r);
        apply_pattern(expected, a, b, kAssx6Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }
    }
}

TEST(FpAvx, Fp4x2wOps)
{
    InsecureRng rng;

    for (int iter = 0; iter < 256; ++iter) {
        FpVec4 a{}, b{};
        for (size_t i = 0; i < 4; ++i) {
            rand_redc_fp(a[i], rng);
            rand_redc_fp(b[i], rng);
        }
        force_edge_lanes(a);
        force_edge_lanes(b);

        fp_4x2w avx_a{}, avx_b{}, avx_r{};
        pack_fp_4x2w(avx_a, a);
        pack_fp_4x2w(avx_b, b);

        FpVec4 got{}, expected{};

        add_fp_4x2w(avx_r, avx_a, avx_b);
        unpack_fp_4x2w(got, avx_r);
        for (size_t lane = 0; lane < 4; ++lane) {
            Fp::add(expected[lane], a[lane], b[lane]);
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        sub_fp_4x2w(avx_r, avx_a, avx_b);
        unpack_fp_4x2w(got, avx_r);
        for (size_t lane = 0; lane < 4; ++lane) {
            Fp::sub(expected[lane], a[lane], b[lane]);
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        asx2_fp_4x2w(avx_r, avx_a, avx_b);
        unpack_fp_4x2w(got, avx_r);
        apply_pattern(expected, a, b, kAsx2Pattern);
        for (size_t lane = 0; lane < 4; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }
    }
}

TEST(FpAvx, MulRedc8x1w)
{
    InsecureRng rng;

    for (int iter = 0; iter < 192; ++iter) {
        FpVec8 a{}, b{};
        for (size_t i = 0; i < 8; ++i) {
            rand_redc_fp(a[i], rng);
            rand_redc_fp(b[i], rng);
        }
        force_edge_lanes(a);
        force_edge_lanes(b);

        fp_8x1w avx_a{}, avx_b{};
        pack_fp_8x1w(avx_a, a);
        pack_fp_8x1w(avx_b, b);

        fpx2_8x1w wide_v1{}, wide_v2{}, wide_v3{}, wide_v4{};
        fp_8x1w red_v1{}, red_v2{}, red_v3{}, red_v4{};
        FpVec8 out_v1{}, out_v2{}, out_v3{}, out_v4{};

        mul_fpx2_8x1w_v1(wide_v1, avx_a, avx_b);
        redc_fpx2_8x1w(red_v1, wide_v1);
        unpack_fp_8x1w(out_v1, red_v1);

        mul_fpx2_8x1w_v2(wide_v2, avx_a, avx_b);
        redc_fpx2_8x1w(red_v2, wide_v2);
        unpack_fp_8x1w(out_v2, red_v2);

        mul_fpx2_8x1w_v3(wide_v3, avx_a, avx_b);
        redc_fpx2_8x1w(red_v3, wide_v3);
        unpack_fp_8x1w(out_v3, red_v3);

        mul_fpx2_8x1w_v4(wide_v4, avx_a, avx_b);
        redc_fpx2_8x1w(red_v4, wide_v4);
        unpack_fp_8x1w(out_v4, red_v4);

        for (size_t lane = 0; lane < 8; ++lane) {
            Fp expected{};
            Fp::mul(expected, a[lane], b[lane]);

            EXPECT_TRUE(Fp::equal(out_v1[lane], expected));
            EXPECT_TRUE(Fp::equal(out_v2[lane], expected));
            EXPECT_TRUE(Fp::equal(out_v3[lane], expected));
            EXPECT_TRUE(Fp::equal(out_v4[lane], expected));

            EXPECT_TRUE(Fp::equal(out_v1[lane], out_v2[lane]));
            EXPECT_TRUE(Fp::equal(out_v1[lane], out_v3[lane]));
            EXPECT_TRUE(Fp::equal(out_v1[lane], out_v4[lane]));
        }
    }
}

TEST(FpAvx, MulRedc4x2w)
{
    InsecureRng rng;

    for (int iter = 0; iter < 192; ++iter) {
        FpVec4 a4{}, b4{};
        for (size_t i = 0; i < 4; ++i) {
            rand_redc_fp(a4[i], rng);
            rand_redc_fp(b4[i], rng);
        }
        force_edge_lanes(a4);
        force_edge_lanes(b4);

        fp_4x2w avx4_a{}, avx4_b{}, avx4_red{};
        fpx2_4x2w avx4_wide{};
        pack_fp_4x2w(avx4_a, a4);
        pack_fp_4x2w(avx4_b, b4);

        mul_fpx2_4x2w_v1(avx4_wide, avx4_a, avx4_b);
        redc_fpx2_4x2w(avx4_red, avx4_wide);

        FpVec4 out4{};
        unpack_fp_4x2w(out4, avx4_red);

        for (size_t lane = 0; lane < 4; ++lane) {
            Fp expected{};
            Fp::mul(expected, a4[lane], b4[lane]);
            EXPECT_TRUE(Fp::equal(out4[lane], expected));
        }

        FpVec8 a8{}, b8{}, out8{};
        for (size_t lane = 0; lane < 4; ++lane) {
            a8[lane] = a4[lane];
            b8[lane] = b4[lane];
        }
        for (size_t lane = 4; lane < 8; ++lane) {
            rand_redc_fp(a8[lane], rng);
            rand_redc_fp(b8[lane], rng);
        }

        fp_8x1w avx8_a{}, avx8_b{}, avx8_red{};
        fpx2_8x1w avx8_wide{};
        pack_fp_8x1w(avx8_a, a8);
        pack_fp_8x1w(avx8_b, b8);

        mul_fpx2_8x1w_v2(avx8_wide, avx8_a, avx8_b);
        redc_fpx2_8x1w(avx8_red, avx8_wide);
        unpack_fp_8x1w(out8, avx8_red);

        for (size_t lane = 0; lane < 4; ++lane) {
            EXPECT_TRUE(Fp::equal(out4[lane], out8[lane]));
        }
    }
}

#else

TEST(FpAvx, Disabled)
{
    SUCCEED();
}

#endif
