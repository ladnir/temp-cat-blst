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

void pack_Fp2Pack8x1x1(Fp2Pack8x1x1& out, const Fp2Vec8& in)
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

void unpack_Fp2Pack8x1x1(Fp2Vec8& out, const Fp2Pack8x1x1& in)
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

void pack_Fp2Pack4x2x1(Fp2Pack4x2x1& out, const Fp2Vec4& in)
{
    FpVec8 lanes{};
    for (size_t i = 0; i < 4; ++i) {
        lanes[2 * i] = in[i][0];
        lanes[2 * i + 1] = in[i][1];
    }
    pack_FpPack8x1(out, lanes);
}

void unpack_Fp2Pack4x2x1(Fp2Vec4& out, const Fp2Pack4x2x1& in)
{
    FpVec8 lanes{};
    unpack_FpPack8x1(lanes, in);
    for (size_t i = 0; i < 4; ++i) {
        out[i][0] = lanes[2 * i];
        out[i][1] = lanes[2 * i + 1];
    }
}

void pack_Fp2Pack2x4x1(Fp2Pack2x4x1& out, const Fp2Vec4& in)
{
    FpVec8 lanes{};
    for (size_t i = 0; i < 4; ++i) {
        lanes[2 * i] = in[i][0];
        lanes[2 * i + 1] = in[i][1];
    }
    pack_FpPack8x1(out, lanes);
}

void unpack_Fp2Pack2x4x1(Fp2Vec4& out, const Fp2Pack2x4x1& in)
{
    FpVec8 lanes{};
    unpack_FpPack8x1(lanes, in);
    for (size_t i = 0; i < 4; ++i) {
        out[i][0] = lanes[2 * i];
        out[i][1] = lanes[2 * i + 1];
    }
}

void pack_Fp2Pack2x2x2(Fp2Pack2x2x2& out, const Fp2Vec2& in)
{
    FpVec4 lanes = {in[0][0], in[0][1], in[1][0], in[1][1]};
    pack_FpPack4x2(out, lanes);
}

void unpack_Fp2Pack2x2x2(Fp2Vec2& out, const Fp2Pack2x2x2& in)
{
    FpVec4 lanes{};
    unpack_FpPack4x2(lanes, in);
    out[0][0] = lanes[0];
    out[0][1] = lanes[1];
    out[1][0] = lanes[2];
    out[1][1] = lanes[3];
}

void pack_Fp2Pack1x4x2(Fp2Pack1x4x2& out, const Fp2Vec2& in)
{
    FpVec4 lanes = {in[0][0], in[0][1], in[1][0], in[1][1]};
    pack_FpPack4x2(out, lanes);
}

void unpack_Fp2Pack1x4x2(Fp2Vec2& out, const Fp2Pack1x4x2& in)
{
    FpVec4 lanes{};
    unpack_FpPack4x2(lanes, in);
    out[0][0] = lanes[0];
    out[0][1] = lanes[1];
    out[1][0] = lanes[2];
    out[1][1] = lanes[3];
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

void apply_lane_mask_select(Fp2Vec8& out, const Fp2Vec8& a, const Fp2Vec8& b, u8 lane_mask)
{
    for (size_t lane = 0; lane < 8; ++lane) {
        const bool select_b = ((lane_mask >> lane) & 1U) != 0;
        out[lane] = select_b ? b[lane] : a[lane];
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

        Fp2Pack8x1x1 avx_a{}, avx_b{}, avx_r{};
        pack_Fp2Pack8x1x1(avx_a, a);
        pack_Fp2Pack8x1x1(avx_b, b);
        Fp2Pack8x1x1::add(avx_r, avx_a, avx_b);
        unpack_Fp2Pack8x1x1(got, avx_r);

        for (size_t i = 0; i < 8; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }
    }
}

TEST(Fp2Avx, Canonical8x1x1wOps)
{
    InsecureRng rng;
    for (int iter = 0; iter < 128; ++iter) {
        Fp2Vec8 a{}, b{}, got{}, expected{};
        for (size_t i = 0; i < 8; ++i) {
            rand_redc_fp2(a[i], rng);
            rand_redc_fp2(b[i], rng);
        }

        Fp2Pack8x1x1 avx_a{}, avx_b{}, avx_r{};
        Fp2Pack8x1x1::from_aos(avx_a, a.data());
        Fp2Pack8x1x1::from_aos(avx_b, b.data());

        Fp2Pack8x1x1::add(avx_r, avx_a, avx_b);
        Fp2Pack8x1x1::to_aos(got.data(), avx_r);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::add(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack8x1x1::sub(avx_r, avx_a, avx_b);
        Fp2Pack8x1x1::to_aos(got.data(), avx_r);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::sub(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack8x1x1::mul(avx_r, avx_a, avx_b);
        Fp2Pack8x1x1::to_aos(got.data(), avx_r);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::mul(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack8x1x1::sqr(avx_r, avx_a);
        Fp2Pack8x1x1::to_aos(got.data(), avx_r);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::sqr(expected[i], a[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack8x1x1::mul_by_u_plus_1(avx_r, avx_a);
        Fp2Pack8x1x1::to_aos(got.data(), avx_r);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::mul_by_u_plus_1(expected[i], a[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack8x1x1::conjugate(avx_r, avx_a);
        Fp2Pack8x1x1::to_aos(got.data(), avx_r);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::conjugate(expected[i], a[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }
    }
}

TEST(Fp2Avx, CanonicalAosPartial8x1x1w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 96; ++iter) {
        Fp2Vec8 src{};
        for (auto& e : src) {
            rand_redc_fp2(e, rng);
        }

        for (size_t n_active = 0; n_active <= 8; ++n_active) {
            Fp2Pack8x1x1 packed{};
            Fp2Vec8 got{};

            Fp2Pack8x1x1::from_aos_partial(packed, src.data(), n_active, Fp2::zero());
            Fp2Pack8x1x1::to_aos(got.data(), packed);
            for (size_t lane = 0; lane < 8; ++lane) {
                const Fp2& expected = (lane < n_active) ? src[lane] : Fp2::zero();
                EXPECT_TRUE(Fp::equal(got[lane][0], expected[0]));
                EXPECT_TRUE(Fp::equal(got[lane][1], expected[1]));
            }

            Fp2Pack8x1x1::from_aos_partial(packed, src.data(), n_active, Fp2::one());
            Fp2Pack8x1x1::to_aos(got.data(), packed);
            for (size_t lane = 0; lane < 8; ++lane) {
                const Fp2& expected = (lane < n_active) ? src[lane] : Fp2::one();
                EXPECT_TRUE(Fp::equal(got[lane][0], expected[0]));
                EXPECT_TRUE(Fp::equal(got[lane][1], expected[1]));
            }

            Fp2Vec8 out_partial{};
            for (auto& e : out_partial) {
                rand_redc_fp2(e, rng);
            }
            const Fp2Vec8 out_before = out_partial;

            Fp2Pack8x1x1::from_aos(packed, src.data());
            Fp2Pack8x1x1::to_aos_partial(out_partial.data(), n_active, packed);
            for (size_t lane = 0; lane < 8; ++lane) {
                if (lane < n_active) {
                    EXPECT_TRUE(Fp::equal(out_partial[lane][0], src[lane][0]));
                    EXPECT_TRUE(Fp::equal(out_partial[lane][1], src[lane][1]));
                } else {
                    EXPECT_TRUE(Fp::equal(out_partial[lane][0], out_before[lane][0]));
                    EXPECT_TRUE(Fp::equal(out_partial[lane][1], out_before[lane][1]));
                }
            }
        }
    }
}

TEST(Fp2Avx, CanonicalSelect8x1x1w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 48; ++iter) {
        Fp2Vec8 a{}, b{}, got{}, expected{};
        for (size_t i = 0; i < 8; ++i) {
            rand_redc_fp2(a[i], rng);
            rand_redc_fp2(b[i], rng);
        }

        Fp2Pack8x1x1 avx_a{}, avx_b{}, avx_r{};
        Fp2Pack8x1x1::from_aos(avx_a, a.data());
        Fp2Pack8x1x1::from_aos(avx_b, b.data());

        for (u16 mask = 0; mask < 256; ++mask) {
            Fp2Pack8x1x1::select(avx_r, avx_a, avx_b, static_cast<u8>(mask));
            Fp2Pack8x1x1::to_aos(got.data(), avx_r);
            apply_lane_mask_select(expected, a, b, static_cast<u8>(mask));
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane][0], expected[lane][0]));
                EXPECT_TRUE(Fp::equal(got[lane][1], expected[lane][1]));
            }

            Fp2Pack8x1x1 avx_alias_a{};
            Fp2Pack8x1x1::from_aos(avx_alias_a, a.data());
            Fp2Pack8x1x1::select(avx_alias_a, avx_alias_a, avx_b, static_cast<u8>(mask));
            Fp2Pack8x1x1::to_aos(got.data(), avx_alias_a);
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane][0], expected[lane][0]));
                EXPECT_TRUE(Fp::equal(got[lane][1], expected[lane][1]));
            }

            Fp2Pack8x1x1 avx_alias_b{};
            Fp2Pack8x1x1::from_aos(avx_alias_b, b.data());
            Fp2Pack8x1x1::select(avx_alias_b, avx_a, avx_alias_b, static_cast<u8>(mask));
            Fp2Pack8x1x1::to_aos(got.data(), avx_alias_b);
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane][0], expected[lane][0]));
                EXPECT_TRUE(Fp::equal(got[lane][1], expected[lane][1]));
            }
        }
    }
}

TEST(Fp2Avx, CanonicalSelectWide8x1x1w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 32; ++iter) {
        Fp2Vec8 a{}, b{}, c{}, d{}, ab{}, cd{}, got{}, expected{};
        for (size_t i = 0; i < 8; ++i) {
            rand_redc_fp2(a[i], rng);
            rand_redc_fp2(b[i], rng);
            rand_redc_fp2(c[i], rng);
            rand_redc_fp2(d[i], rng);
            Fp2::mul(ab[i], a[i], b[i]);
            Fp2::mul(cd[i], c[i], d[i]);
        }

        Fp2Pack8x1x1 avx_a{}, avx_b{}, avx_c{}, avx_d{}, red{};
        UFp2Pack8x1x1 wide_ab{}, wide_cd{}, wide_sel{};
        Fp2Pack8x1x1::from_aos(avx_a, a.data());
        Fp2Pack8x1x1::from_aos(avx_b, b.data());
        Fp2Pack8x1x1::from_aos(avx_c, c.data());
        Fp2Pack8x1x1::from_aos(avx_d, d.data());
        UFp2Pack8x1x1::mul(wide_ab, avx_a, avx_b);
        UFp2Pack8x1x1::mul(wide_cd, avx_c, avx_d);

        for (u16 mask = 0; mask < 256; ++mask) {
            UFp2Pack8x1x1::select(wide_sel, wide_ab, wide_cd, static_cast<u8>(mask));
            UFp2Pack8x1x1::redc(red, wide_sel);
            Fp2Pack8x1x1::to_aos(got.data(), red);
            apply_lane_mask_select(expected, ab, cd, static_cast<u8>(mask));
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane][0], expected[lane][0]));
                EXPECT_TRUE(Fp::equal(got[lane][1], expected[lane][1]));
            }

            UFp2Pack8x1x1 wide_alias_a = wide_ab;
            UFp2Pack8x1x1::select(wide_alias_a, wide_alias_a, wide_cd, static_cast<u8>(mask));
            UFp2Pack8x1x1::redc(red, wide_alias_a);
            Fp2Pack8x1x1::to_aos(got.data(), red);
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane][0], expected[lane][0]));
                EXPECT_TRUE(Fp::equal(got[lane][1], expected[lane][1]));
            }

            UFp2Pack8x1x1 wide_alias_b = wide_cd;
            UFp2Pack8x1x1::select(wide_alias_b, wide_ab, wide_alias_b, static_cast<u8>(mask));
            UFp2Pack8x1x1::redc(red, wide_alias_b);
            Fp2Pack8x1x1::to_aos(got.data(), red);
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane][0], expected[lane][0]));
                EXPECT_TRUE(Fp::equal(got[lane][1], expected[lane][1]));
            }
        }
    }
}

TEST(Fp2Avx, Redc8x1x1wMatchesSerialComponentReduction)
{
    InsecureRng rng;
    for (int iter = 0; iter < 128; ++iter) {
        Fp2Vec8 a{}, b{}, got{}, expected{};
        for (size_t i = 0; i < 8; ++i) {
            rand_redc_fp2(a[i], rng);
            rand_redc_fp2(b[i], rng);
        }

        Fp2Pack8x1x1 avx_a{}, avx_b{}, red{}, expected_red{};
        UFp2Pack8x1x1 wide{};
        Fp2Pack8x1x1::from_aos(avx_a, a.data());
        Fp2Pack8x1x1::from_aos(avx_b, b.data());

        UFp2Pack8x1x1::mul(wide, avx_a, avx_b);
        UFp2Pack8x1x1::redc(red, wide);
        UFpPack8x1::redc(expected_red[0], wide[0]);
        UFpPack8x1::redc(expected_red[1], wide[1]);

        Fp2Pack8x1x1::to_aos(got.data(), red);
        Fp2Pack8x1x1::to_aos(expected.data(), expected_red);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane][0], expected[lane][0]));
            EXPECT_TRUE(Fp::equal(got[lane][1], expected[lane][1]));
        }
    }
}

TEST(Fp2Avx, CanonicalPartialTailMul8x1x1w)
{
    InsecureRng rng;
    for (int iter = 0; iter < 128; ++iter) {
        Fp2Vec8 a{}, b{}, out{};
        for (size_t i = 0; i < 8; ++i) {
            rand_redc_fp2(a[i], rng);
            rand_redc_fp2(b[i], rng);
        }

        for (size_t n_active = 0; n_active <= 8; ++n_active) {
            Fp2Pack8x1x1 avx_a{}, avx_b{}, avx_r{};
            Fp2Pack8x1x1::from_aos_partial(avx_a, a.data(), n_active, Fp2::one());
            Fp2Pack8x1x1::from_aos_partial(avx_b, b.data(), n_active, Fp2::one());
            Fp2Pack8x1x1::mul(avx_r, avx_a, avx_b);

            for (auto& e : out) {
                rand_redc_fp2(e, rng);
            }
            const Fp2Vec8 out_before = out;

            Fp2Pack8x1x1::to_aos_partial(out.data(), n_active, avx_r);
            for (size_t lane = 0; lane < 8; ++lane) {
                if (lane < n_active) {
                    Fp2 expected{};
                    Fp2::mul(expected, a[lane], b[lane]);
                    EXPECT_TRUE(Fp::equal(out[lane][0], expected[0]));
                    EXPECT_TRUE(Fp::equal(out[lane][1], expected[1]));
                } else {
                    EXPECT_TRUE(Fp::equal(out[lane][0], out_before[lane][0]));
                    EXPECT_TRUE(Fp::equal(out[lane][1], out_before[lane][1]));
                }
            }
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

        Fp2Pack4x2x1 avx_a{}, avx_b{}, avx_r{};
        pack_Fp2Pack4x2x1(avx_a, a);
        pack_Fp2Pack4x2x1(avx_b, b);

        Fp2Pack4x2x1::add(avx_r, avx_a, avx_b);
        unpack_Fp2Pack4x2x1(got, avx_r);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::add(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack4x2x1::sub(avx_r, avx_a, avx_b);
        unpack_Fp2Pack4x2x1(got, avx_r);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::sub(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack4x2x1::aass(avx_r, avx_a, avx_b);
        unpack_Fp2Pack4x2x1(got, avx_r);
        apply_pattern(expected, a, b, kAass4);
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack4x2x1::assa(avx_r, avx_a, avx_b);
        unpack_Fp2Pack4x2x1(got, avx_r);
        apply_pattern(expected, a, b, kAssa4);
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack4x2x1::sass(avx_r, avx_a, avx_b);
        unpack_Fp2Pack4x2x1(got, avx_r);
        apply_pattern(expected, a, b, kSass4);
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack4x2x1::ssaa(avx_r, avx_a, avx_b);
        unpack_Fp2Pack4x2x1(got, avx_r);
        apply_pattern(expected, a, b, kSsaa4);
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack4x2x1::mul_by_u_plus_1(avx_r, avx_a);
        unpack_Fp2Pack4x2x1(got, avx_r);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::mul_by_u_plus_1(expected[i], a[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack4x2x1::mul(avx_r, avx_a, avx_b);
        unpack_Fp2Pack4x2x1(got, avx_r);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::mul(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack4x2x1::sqr(avx_r, avx_a);
        unpack_Fp2Pack4x2x1(got, avx_r);
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

        Fp2Pack2x2x2 avx_a{}, avx_b{}, avx_r{};
        pack_Fp2Pack2x2x2(avx_a, a);
        pack_Fp2Pack2x2x2(avx_b, b);

        Fp2Pack2x2x2::add(avx_r, avx_a, avx_b);
        unpack_Fp2Pack2x2x2(got, avx_r);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::add(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack2x2x2::sub(avx_r, avx_a, avx_b);
        unpack_Fp2Pack2x2x2(got, avx_r);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::sub(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack2x2x2::as(avx_r, avx_a, avx_b);
        unpack_Fp2Pack2x2x2(got, avx_r);
        apply_pattern(expected, a, b, kAs2);
        for (size_t i = 0; i < 2; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack2x2x2::sa(avx_r, avx_a, avx_b);
        unpack_Fp2Pack2x2x2(got, avx_r);
        apply_pattern(expected, a, b, kSa2);
        for (size_t i = 0; i < 2; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack2x2x2::mul_by_u_plus_1(avx_r, avx_a);
        unpack_Fp2Pack2x2x2(got, avx_r);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::mul_by_u_plus_1(expected[i], a[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack2x2x2::mul(avx_r, avx_a, avx_b);
        unpack_Fp2Pack2x2x2(got, avx_r);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::mul(expected[i], a[i], b[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        Fp2Pack2x2x2::sqr(avx_r, avx_a);
        unpack_Fp2Pack2x2x2(got, avx_r);
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

        Fp2Pack8x1x1 avx_a{}, avx_b{}, avx_c{}, avx_d{}, red{};
        UFp2Pack8x1x1 wide_ab{}, wide_cd{}, wide_r{};
        pack_Fp2Pack8x1x1(avx_a, a);
        pack_Fp2Pack8x1x1(avx_b, b);
        pack_Fp2Pack8x1x1(avx_c, c);
        pack_Fp2Pack8x1x1(avx_d, d);

        UFp2Pack8x1x1::mul(wide_ab, avx_a, avx_b);
        UFp2Pack8x1x1::mul(wide_cd, avx_c, avx_d);

        UFp2Pack8x1x1::add(wide_r, wide_ab, wide_cd);
        UFp2Pack8x1x1::redc(red, wide_r);
        unpack_Fp2Pack8x1x1(got, red);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::add(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        UFp2Pack8x1x1::sub(wide_r, wide_ab, wide_cd);
        UFp2Pack8x1x1::redc(red, wide_r);
        unpack_Fp2Pack8x1x1(got, red);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::sub(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        UFp2Pack8x1x1::asx4(wide_r, wide_ab, wide_cd);
        UFp2Pack8x1x1::redc(red, wide_r);
        unpack_Fp2Pack8x1x1(got, red);
        apply_pattern(expected, ab, cd, kAsx4_8);
        for (size_t i = 0; i < 8; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        UFp2Pack8x1x1::mul_by_u_plus_1(wide_r, wide_ab);
        UFp2Pack8x1x1::redc(red, wide_r);
        unpack_Fp2Pack8x1x1(got, red);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::mul_by_u_plus_1(expected[i], ab[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        UFp2Pack8x1x1::sqr(wide_r, avx_a);
        UFp2Pack8x1x1::redc(red, wide_r);
        unpack_Fp2Pack8x1x1(got, red);
        for (size_t i = 0; i < 8; ++i) {
            Fp2::sqr(expected[i], a[i]);
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

        Fp2Pack4x2x1 avx_a{}, avx_b{}, avx_c{}, avx_d{}, red{};
        UFp2Pack4x2x1 wide_ab{}, wide_cd{}, wide_r{};
        pack_Fp2Pack4x2x1(avx_a, a);
        pack_Fp2Pack4x2x1(avx_b, b);
        pack_Fp2Pack4x2x1(avx_c, c);
        pack_Fp2Pack4x2x1(avx_d, d);

        UFp2Pack4x2x1::mul(wide_ab, avx_a, avx_b);
        UFp2Pack4x2x1::mul(wide_cd, avx_c, avx_d);

        UFp2Pack4x2x1::add(wide_r, wide_ab, wide_cd);
        UFp2Pack4x2x1::redc(red, wide_r);
        unpack_Fp2Pack4x2x1(got, red);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::add(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        UFp2Pack4x2x1::sub(wide_r, wide_ab, wide_cd);
        UFp2Pack4x2x1::redc(red, wide_r);
        unpack_Fp2Pack4x2x1(got, red);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::sub(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        UFp2Pack4x2x1::sax2(wide_r, wide_ab, wide_cd);
        UFp2Pack4x2x1::redc(red, wide_r);
        unpack_Fp2Pack4x2x1(got, red);
        apply_pattern(expected, ab, cd, kSax2_4);
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        UFp2Pack4x2x1::mul_by_u_plus_1(wide_r, wide_ab);
        UFp2Pack4x2x1::redc(red, wide_r);
        unpack_Fp2Pack4x2x1(got, red);
        for (size_t i = 0; i < 4; ++i) {
            Fp2::mul_by_u_plus_1(expected[i], ab[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        UFp2Pack4x2x1::sqr(wide_r, avx_a);
        UFp2Pack4x2x1::redc(red, wide_r);
        unpack_Fp2Pack4x2x1(got, red);
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

        Fp2Pack2x2x2 avx_a{}, avx_b{}, avx_c{}, avx_d{}, red{};
        UFp2Pack2x2x2 wide_ab{}, wide_cd{}, wide_r{};
        pack_Fp2Pack2x2x2(avx_a, a);
        pack_Fp2Pack2x2x2(avx_b, b);
        pack_Fp2Pack2x2x2(avx_c, c);
        pack_Fp2Pack2x2x2(avx_d, d);

        UFp2Pack2x2x2::mul(wide_ab, avx_a, avx_b);
        UFp2Pack2x2x2::mul(wide_cd, avx_c, avx_d);

        UFp2Pack2x2x2::add(wide_r, wide_ab, wide_cd);
        UFp2Pack2x2x2::redc(red, wide_r);
        unpack_Fp2Pack2x2x2(got, red);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::add(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        UFp2Pack2x2x2::sub(wide_r, wide_ab, wide_cd);
        UFp2Pack2x2x2::redc(red, wide_r);
        unpack_Fp2Pack2x2x2(got, red);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::sub(expected[i], ab[i], cd[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        UFp2Pack2x2x2::mul_by_u_plus_1(wide_r, wide_ab);
        UFp2Pack2x2x2::redc(red, wide_r);
        unpack_Fp2Pack2x2x2(got, red);
        for (size_t i = 0; i < 2; ++i) {
            Fp2::mul_by_u_plus_1(expected[i], ab[i]);
            EXPECT_TRUE(Fp::equal(got[i][0], expected[i][0]));
            EXPECT_TRUE(Fp::equal(got[i][1], expected[i][1]));
        }

        UFp2Pack2x2x2::sqr(wide_r, avx_a);
        UFp2Pack2x2x2::redc(red, wide_r);
        unpack_Fp2Pack2x2x2(got, red);
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

        Fp2Pack2x4x1 avx_a{}, avx_b{}, red_direct{}, red_from_wide{};
        UFp2Pack2x4x1 wide{};
        pack_Fp2Pack2x4x1(avx_a, a);
        pack_Fp2Pack2x4x1(avx_b, b);

        Fp2Pack2x4x1::mul(red_direct, avx_a, avx_b);
        UFp2Pack2x4x1::mul(wide, avx_a, avx_b);
        UFpPack8x1::redc(red_from_wide, wide);

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
        Fp2Pack2x4x1 avx_a{}, avx_b{}, avx_r{};

        pack_Fp2Pack2x4x1(avx_a, lanes_a);
        pack_Fp2Pack2x4x1(avx_b, lanes_b);
        Fp2Pack2x4x1::mul(avx_r, avx_a, avx_b);
        unpack_Fp2Pack2x4x1(got, avx_r);

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

        Fp2Pack1x4x2 avx_a{}, avx_b{}, red_direct{}, red_from_wide{};
        UFp2Pack1x4x2 wide{};
        pack_Fp2Pack1x4x2(avx_a, a);
        pack_Fp2Pack1x4x2(avx_b, b);

        Fp2Pack1x4x2::mul(red_direct, avx_a, avx_b);
        UFp2Pack1x4x2::mul(wide, avx_a, avx_b);
        UFpPack4x2::redc(red_from_wide, wide);

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
        Fp2Pack1x4x2 avx_a{}, avx_b{}, avx_r{};

        pack_Fp2Pack1x4x2(avx_a, lanes_a);
        pack_Fp2Pack1x4x2(avx_b, lanes_b);
        Fp2Pack1x4x2::mul(avx_r, avx_a, avx_b);
        unpack_Fp2Pack1x4x2(got, avx_r);

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
