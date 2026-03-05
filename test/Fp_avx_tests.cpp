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

constexpr std::array<LaneOp, 8> kAllAdd8Pattern = {
    LaneOp::Add, LaneOp::Add, LaneOp::Add, LaneOp::Add,
    LaneOp::Add, LaneOp::Add, LaneOp::Add, LaneOp::Add,
};

constexpr std::array<LaneOp, 8> kAllSub8Pattern = {
    LaneOp::Sub, LaneOp::Sub, LaneOp::Sub, LaneOp::Sub,
    LaneOp::Sub, LaneOp::Sub, LaneOp::Sub, LaneOp::Sub,
};

constexpr std::array<LaneOp, 8> kAaasx2Pattern = {
    LaneOp::Sub, LaneOp::Add, LaneOp::Add, LaneOp::Add,
    LaneOp::Sub, LaneOp::Add, LaneOp::Add, LaneOp::Add,
};

constexpr std::array<LaneOp, 8> kAassx2Pattern = {
    LaneOp::Sub, LaneOp::Sub, LaneOp::Add, LaneOp::Add,
    LaneOp::Sub, LaneOp::Sub, LaneOp::Add, LaneOp::Add,
};

constexpr std::array<LaneOp, 8> kSsaax2Pattern = {
    LaneOp::Add, LaneOp::Add, LaneOp::Sub, LaneOp::Sub,
    LaneOp::Add, LaneOp::Add, LaneOp::Sub, LaneOp::Sub,
};

constexpr std::array<LaneOp, 8> kSsasx2Pattern = {
    LaneOp::Sub, LaneOp::Add, LaneOp::Sub, LaneOp::Sub,
    LaneOp::Sub, LaneOp::Add, LaneOp::Sub, LaneOp::Sub,
};

constexpr std::array<LaneOp, 4> kAllAdd4Pattern = {
    LaneOp::Add, LaneOp::Add, LaneOp::Add, LaneOp::Add,
};

constexpr std::array<LaneOp, 4> kAllSub4Pattern = {
    LaneOp::Sub, LaneOp::Sub, LaneOp::Sub, LaneOp::Sub,
};

constexpr std::array<LaneOp, 4> kAass4Pattern = {
    LaneOp::Sub, LaneOp::Sub, LaneOp::Add, LaneOp::Add,
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

void pack_FpPack8x1(FpPack8x1& out, const FpVec8& in)
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

void unpack_FpPack8x1(FpVec8& out, const FpPack8x1& in)
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

void pack_FpPack4x2(FpPack4x2& out, const FpVec4& in)
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

void unpack_FpPack4x2(FpVec4& out, const FpPack4x2& in)
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

void apply_lane_mask_select(FpVec8& out, const FpVec8& a, const FpVec8& b, u8 lane_mask)
{
    for (size_t lane = 0; lane < 8; ++lane) {
        const bool select_b = ((lane_mask >> lane) & 1U) != 0;
        out[lane] = select_b ? b[lane] : a[lane];
    }
}

template <int K>
void run_kwise_mul_redc_matches_scalar(InsecureRng& rng)
{
    static_assert((K == 1) || (K == 2) || (K == 4));
    constexpr size_t kBlocks = 3;
    constexpr size_t kPacks = K * kBlocks;

    std::array<FpVec8, kPacks> a{};
    std::array<FpVec8, kPacks> b{};
    std::array<FpVec8, kPacks> expected{};
    std::array<FpPack8x1, kPacks> pa{};
    std::array<FpPack8x1, kPacks> pb{};
    std::array<FpPack8x1, kPacks> got_from_split{};
    std::array<FpPack8x1, kPacks> got_from_fused{};
    std::array<UFpPack8x1, kPacks> wide{};

    for (size_t pack = 0; pack < kPacks; ++pack) {
        for (size_t lane = 0; lane < 8; ++lane) {
            rand_redc_fp(a[pack][lane], rng);
            rand_redc_fp(b[pack][lane], rng);
            Fp::mul(expected[pack][lane], a[pack][lane], b[pack][lane]);
        }
        force_edge_lanes(a[pack]);
        force_edge_lanes(b[pack]);
        for (size_t lane = 0; lane < 8; ++lane) {
            Fp::mul(expected[pack][lane], a[pack][lane], b[pack][lane]);
        }
        pack_FpPack8x1(pa[pack], a[pack]);
        pack_FpPack8x1(pb[pack], b[pack]);
    }

    fp8x1_mul_k<K>(wide.data(), pa.data(), pb.data(), kBlocks);
    fp8x1_redc_k<K>(got_from_split.data(), wide.data(), kBlocks);

    fp8x1_mul_redc_k<K>(got_from_fused.data(), pa.data(), pb.data(), kBlocks);

    for (size_t pack = 0; pack < kPacks; ++pack) {
        FpVec8 split_lanes{};
        FpVec8 fused_lanes{};
        unpack_FpPack8x1(split_lanes, got_from_split[pack]);
        unpack_FpPack8x1(fused_lanes, got_from_fused[pack]);

        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(split_lanes[lane], expected[pack][lane]));
            EXPECT_TRUE(Fp::equal(fused_lanes[lane], expected[pack][lane]));
            EXPECT_TRUE(Fp::equal(split_lanes[lane], fused_lanes[lane]));
        }
    }
}

void run_kwise_dispatch_matches_scalar(InsecureRng& rng)
{
    constexpr size_t kPacks = 7;
    std::array<FpVec8, kPacks> a{};
    std::array<FpVec8, kPacks> b{};
    std::array<FpVec8, kPacks> expected{};
    std::array<FpPack8x1, kPacks> pa{};
    std::array<FpPack8x1, kPacks> pb{};
    std::array<FpPack8x1, kPacks> out_direct{};
    std::array<FpPack8x1, kPacks> out_dispatch{};
    std::array<UFpPack8x1, kPacks> wide_direct{};
    std::array<UFpPack8x1, kPacks> wide_dispatch{};

    for (size_t pack = 0; pack < kPacks; ++pack) {
        for (size_t lane = 0; lane < 8; ++lane) {
            rand_redc_fp(a[pack][lane], rng);
            rand_redc_fp(b[pack][lane], rng);
            Fp::mul(expected[pack][lane], a[pack][lane], b[pack][lane]);
        }
        force_edge_lanes(a[pack]);
        force_edge_lanes(b[pack]);
        for (size_t lane = 0; lane < 8; ++lane) {
            Fp::mul(expected[pack][lane], a[pack][lane], b[pack][lane]);
        }

        pack_FpPack8x1(pa[pack], a[pack]);
        pack_FpPack8x1(pb[pack], b[pack]);
    }

    for (size_t pack = 0; pack < kPacks; ++pack) {
        FpPack8x1::mul(wide_direct[pack], pa[pack], pb[pack]);
        UFpPack8x1::redc(out_direct[pack], wide_direct[pack]);
    }

    fp8x1_mul_dispatch(wide_dispatch.data(), pa.data(), pb.data(), kPacks);
    fp8x1_redc_dispatch(out_dispatch.data(), wide_dispatch.data(), kPacks);
    for (size_t pack = 0; pack < kPacks; ++pack) {
        FpVec8 got{};
        FpVec8 ref{};
        unpack_FpPack8x1(got, out_dispatch[pack]);
        unpack_FpPack8x1(ref, out_direct[pack]);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], ref[lane]));
            EXPECT_TRUE(Fp::equal(got[lane], expected[pack][lane]));
        }
    }

    fp8x1_mul_redc_dispatch(out_dispatch.data(), pa.data(), pb.data(), kPacks);
    for (size_t pack = 0; pack < kPacks; ++pack) {
        FpVec8 got{};
        unpack_FpPack8x1(got, out_dispatch[pack]);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[pack][lane]));
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

        FpPack8x1 fp8{};
        cc512i out64_8[6];
        FpPack8x1::conv_64to48(fp8, in64_8);
        FpPack8x1::conv_48to64(out64_8, fp8);

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

        FpPack4x2 fp4{};
        cc512i out64_4[3];
        FpPack4x2::conv_64to48(fp4, in64_4);
        FpPack4x2::conv_48to64(out64_4, fp4);

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

        FpPack8x1 avx_a{}, avx_b{}, avx_r{};
        pack_FpPack8x1(avx_a, a);
        pack_FpPack8x1(avx_b, b);

        FpVec8 got{}, expected{};

        FpPack8x1::add(avx_r, avx_a, avx_b);
        unpack_FpPack8x1(got, avx_r);
        for (size_t lane = 0; lane < 8; ++lane) {
            Fp::add(expected[lane], a[lane], b[lane]);
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        FpPack8x1::sub(avx_r, avx_a, avx_b);
        unpack_FpPack8x1(got, avx_r);
        for (size_t lane = 0; lane < 8; ++lane) {
            Fp::sub(expected[lane], a[lane], b[lane]);
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        FpPack8x1::asx4(avx_r, avx_a, avx_b);
        unpack_FpPack8x1(got, avx_r);
        apply_pattern(expected, a, b, kAsx4Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        FpPack8x1::asax6(avx_r, avx_a, avx_b);
        unpack_FpPack8x1(got, avx_r);
        apply_pattern(expected, a, b, kAsax6Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        FpPack8x1::assx6(avx_r, avx_a, avx_b);
        unpack_FpPack8x1(got, avx_r);
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

        FpPack4x2 avx_a{}, avx_b{}, avx_r{};
        pack_FpPack4x2(avx_a, a);
        pack_FpPack4x2(avx_b, b);

        FpVec4 got{}, expected{};

        FpPack4x2::add(avx_r, avx_a, avx_b);
        unpack_FpPack4x2(got, avx_r);
        for (size_t lane = 0; lane < 4; ++lane) {
            Fp::add(expected[lane], a[lane], b[lane]);
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        FpPack4x2::sub(avx_r, avx_a, avx_b);
        unpack_FpPack4x2(got, avx_r);
        for (size_t lane = 0; lane < 4; ++lane) {
            Fp::sub(expected[lane], a[lane], b[lane]);
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        FpPack4x2::asx2(avx_r, avx_a, avx_b);
        unpack_FpPack4x2(got, avx_r);
        apply_pattern(expected, a, b, kAsx2Pattern);
        for (size_t lane = 0; lane < 4; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }
    }
}

TEST(FpAvx, Widened8x1wOps)
{
    InsecureRng rng;

    for (int iter = 0; iter < 128; ++iter) {
        FpVec8 a{}, b{}, c{}, d{};
        for (size_t i = 0; i < 8; ++i) {
            rand_redc_fp(a[i], rng);
            rand_redc_fp(b[i], rng);
            rand_redc_fp(c[i], rng);
            rand_redc_fp(d[i], rng);
        }
        force_edge_lanes(a);
        force_edge_lanes(b);
        force_edge_lanes(c);
        force_edge_lanes(d);

        FpPack8x1 avx_a{}, avx_b{}, avx_c{}, avx_d{};
        UFpPack8x1 wide_ab{}, wide_cd{}, wide_r{};
        FpPack8x1 red_r{};
        FpVec8 mul_ab{}, mul_cd{}, got{}, expected{};

        pack_FpPack8x1(avx_a, a);
        pack_FpPack8x1(avx_b, b);
        pack_FpPack8x1(avx_c, c);
        pack_FpPack8x1(avx_d, d);

        FpPack8x1::mul(wide_ab, avx_a, avx_b);
        FpPack8x1::mul(wide_cd, avx_c, avx_d);

        for (size_t lane = 0; lane < 8; ++lane) {
            Fp::mul(mul_ab[lane], a[lane], b[lane]);
            Fp::mul(mul_cd[lane], c[lane], d[lane]);
        }

        UFpPack8x1::add(wide_r, wide_ab, wide_cd);
        UFpPack8x1::redc(red_r, wide_r);
        unpack_FpPack8x1(got, red_r);
        apply_pattern(expected, mul_ab, mul_cd, kAllAdd8Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        UFpPack8x1::sub(wide_r, wide_ab, wide_cd);
        UFpPack8x1::redc(red_r, wide_r);
        unpack_FpPack8x1(got, red_r);
        apply_pattern(expected, mul_ab, mul_cd, kAllSub8Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        UFpPack8x1::aaasx2(wide_r, wide_ab, wide_cd);
        UFpPack8x1::redc(red_r, wide_r);
        unpack_FpPack8x1(got, red_r);
        apply_pattern(expected, mul_ab, mul_cd, kAaasx2Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        UFpPack8x1::aassx2(wide_r, wide_ab, wide_cd);
        UFpPack8x1::redc(red_r, wide_r);
        unpack_FpPack8x1(got, red_r);
        apply_pattern(expected, mul_ab, mul_cd, kAassx2Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        UFpPack8x1::asx4(wide_r, wide_ab, wide_cd);
        UFpPack8x1::redc(red_r, wide_r);
        unpack_FpPack8x1(got, red_r);
        apply_pattern(expected, mul_ab, mul_cd, kAsx4Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        UFpPack8x1::ssaax2(wide_r, wide_ab, wide_cd);
        UFpPack8x1::redc(red_r, wide_r);
        unpack_FpPack8x1(got, red_r);
        apply_pattern(expected, mul_ab, mul_cd, kSsaax2Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        UFpPack8x1::ssasx2(wide_r, wide_ab, wide_cd);
        UFpPack8x1::redc(red_r, wide_r);
        unpack_FpPack8x1(got, red_r);
        apply_pattern(expected, mul_ab, mul_cd, kSsasx2Pattern);
        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }
    }
}

TEST(FpAvx, Widened4x2wOps)
{
    InsecureRng rng;

    for (int iter = 0; iter < 128; ++iter) {
        FpVec4 a{}, b{}, c{}, d{};
        for (size_t i = 0; i < 4; ++i) {
            rand_redc_fp(a[i], rng);
            rand_redc_fp(b[i], rng);
            rand_redc_fp(c[i], rng);
            rand_redc_fp(d[i], rng);
        }
        force_edge_lanes(a);
        force_edge_lanes(b);
        force_edge_lanes(c);
        force_edge_lanes(d);

        FpPack4x2 avx_a{}, avx_b{}, avx_c{}, avx_d{};
        UFpPack4x2 wide_ab{}, wide_cd{}, wide_r{};
        FpPack4x2 red_r{};
        FpVec4 mul_ab{}, mul_cd{}, got{}, expected{};

        pack_FpPack4x2(avx_a, a);
        pack_FpPack4x2(avx_b, b);
        pack_FpPack4x2(avx_c, c);
        pack_FpPack4x2(avx_d, d);

        FpPack4x2::mul_v1(wide_ab, avx_a, avx_b);
        FpPack4x2::mul_v1(wide_cd, avx_c, avx_d);

        for (size_t lane = 0; lane < 4; ++lane) {
            Fp::mul(mul_ab[lane], a[lane], b[lane]);
            Fp::mul(mul_cd[lane], c[lane], d[lane]);
        }

        UFpPack4x2::add(wide_r, wide_ab, wide_cd);
        UFpPack4x2::redc(red_r, wide_r);
        unpack_FpPack4x2(got, red_r);
        apply_pattern(expected, mul_ab, mul_cd, kAllAdd4Pattern);
        for (size_t lane = 0; lane < 4; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        UFpPack4x2::sub(wide_r, wide_ab, wide_cd);
        UFpPack4x2::redc(red_r, wide_r);
        unpack_FpPack4x2(got, red_r);
        apply_pattern(expected, mul_ab, mul_cd, kAllSub4Pattern);
        for (size_t lane = 0; lane < 4; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        UFpPack4x2::asx2(wide_r, wide_ab, wide_cd);
        UFpPack4x2::redc(red_r, wide_r);
        unpack_FpPack4x2(got, red_r);
        apply_pattern(expected, mul_ab, mul_cd, kAsx2Pattern);
        for (size_t lane = 0; lane < 4; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }

        UFpPack4x2::aass(wide_r, wide_ab, wide_cd);
        UFpPack4x2::redc(red_r, wide_r);
        unpack_FpPack4x2(got, red_r);
        apply_pattern(expected, mul_ab, mul_cd, kAass4Pattern);
        for (size_t lane = 0; lane < 4; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
        }
    }
}

TEST(FpAvx, LayoutHelpers)
{
    InsecureRng rng;
    const auto sra48 = [](u64 v) -> u64 {
        return static_cast<u64>(static_cast<i64>(v) >> BRADIX);
    };

    for (int iter = 0; iter < 128; ++iter) {
        std::array<u64, NWORDS> carry_in{};
        for (size_t i = 0; i < NWORDS; ++i) {
            carry_in[i] = rng.next_u64();
        }
        std::array<u64, NWORDS> carry_ref = carry_in;

        for (size_t i = 0; i < NWORDS - 1; ++i) {
            carry_ref[i + 1] += (carry_ref[i] >> BRADIX);
            carry_ref[i] &= BMASK;
        }

        carryp_mpi48(carry_in.data());
        for (size_t i = 0; i < NWORDS; ++i) {
            EXPECT_EQ(carry_in[i], carry_ref[i]);
        }

        FpPack4x2 carry4{};
        std::array<std::array<u64, 8>, VWORDS> carry4_ref{};
        for (size_t i = 0; i < VWORDS; ++i) {
            std::array<u64, 8> lanes{};
            for (size_t lane = 0; lane < 8; ++lane) {
                lanes[lane] = rng.next_u64();
                carry4_ref[i][lane] = lanes[lane];
            }
            carry4[i] = make_vec(lanes);
        }

        for (size_t lane = 0; lane < 8; ++lane) {
            carry4_ref[1][lane] += sra48(carry4_ref[0][lane]);
            carry4_ref[0][lane] &= BMASK;
            carry4_ref[2][lane] += sra48(carry4_ref[1][lane]);
            carry4_ref[1][lane] &= BMASK;
            carry4_ref[3][lane] += sra48(carry4_ref[2][lane]);
            carry4_ref[2][lane] &= BMASK;
        }

        FpPack4x2::carryp(carry4);
        for (size_t i = 0; i < VWORDS; ++i) {
            const auto got = load_lanes(carry4[i]);
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_EQ(got[lane], carry4_ref[i][lane]);
            }
        }

        FpPack8x1 ch8{};
        std::array<std::array<u64, 8>, NWORDS> ch8_ref{};
        for (size_t i = 0; i < NWORDS; ++i) {
            std::array<u64, 8> lanes{};
            for (size_t lane = 0; lane < 8; ++lane) {
                lanes[lane] = rng.next_u64();
                ch8_ref[i][lane] = lanes[lane];
            }
            ch8[i] = make_vec(lanes);
        }

        for (size_t ch = 0; ch < 8; ++ch) {
            std::array<u64, NWORDS> out{};
            FpPack8x1::get_channel(out.data(), ch8, static_cast<int>(ch));
            for (size_t i = 0; i < NWORDS; ++i) {
                EXPECT_EQ(out[i], ch8_ref[i][ch]);
            }
        }

        FpPack4x2 ch4{};
        std::array<std::array<u64, 8>, VWORDS> ch4_ref{};
        for (size_t i = 0; i < VWORDS; ++i) {
            std::array<u64, 8> lanes{};
            for (size_t lane = 0; lane < 8; ++lane) {
                lanes[lane] = rng.next_u64();
                ch4_ref[i][lane] = lanes[lane];
            }
            ch4[i] = make_vec(lanes);
        }

        for (size_t ch = 0; ch <= 6; ch += 2) {
            std::array<u64, NWORDS> out{};
            FpPack4x2::get_channel(out.data(), ch4, static_cast<int>(ch));
            for (size_t i = 0; i < VWORDS; ++i) {
                EXPECT_EQ(out[i], ch4_ref[i][ch]);
                EXPECT_EQ(out[i + VWORDS], ch4_ref[i][ch + 1]);
            }
        }
    }
}

TEST(FpAvx, AosRoundTrip8x1)
{
    InsecureRng rng;

    for (int iter = 0; iter < 256; ++iter) {
        FpVec8 src{};
        for (auto& e : src) {
            rand_redc_fp(e, rng);
        }
        force_edge_lanes(src);

        FpPack8x1 packed{};
        FpVec8 got{};
        FpPack8x1::from_aos(packed, src.data());
        FpPack8x1::to_aos(got.data(), packed);

        for (size_t lane = 0; lane < 8; ++lane) {
            EXPECT_TRUE(Fp::equal(got[lane], src[lane]));
        }
    }
}

TEST(FpAvx, AosPartialPackUnpack8x1)
{
    InsecureRng rng;

    for (int iter = 0; iter < 128; ++iter) {
        FpVec8 src{};
        for (auto& e : src) {
            rand_redc_fp(e, rng);
        }
        force_edge_lanes(src);

        for (size_t n_active = 0; n_active <= 8; ++n_active) {
            FpPack8x1 packed{};
            FpVec8 got{};

            FpPack8x1::from_aos_partial(packed, src.data(), n_active, Fp::zero());
            FpPack8x1::to_aos(got.data(), packed);
            for (size_t lane = 0; lane < 8; ++lane) {
                const Fp& expected = (lane < n_active) ? src[lane] : Fp::zero();
                EXPECT_TRUE(Fp::equal(got[lane], expected));
            }

            FpPack8x1::from_aos_partial(packed, src.data(), n_active, Fp::one());
            FpPack8x1::to_aos(got.data(), packed);
            for (size_t lane = 0; lane < 8; ++lane) {
                const Fp& expected = (lane < n_active) ? src[lane] : Fp::one();
                EXPECT_TRUE(Fp::equal(got[lane], expected));
            }

            FpVec8 out_partial{};
            for (auto& e : out_partial) {
                rand_redc_fp(e, rng);
            }
            const FpVec8 out_before = out_partial;

            FpPack8x1::from_aos(packed, src.data());
            FpPack8x1::to_aos_partial(out_partial.data(), n_active, packed);
            for (size_t lane = 0; lane < 8; ++lane) {
                if (lane < n_active) {
                    EXPECT_TRUE(Fp::equal(out_partial[lane], src[lane]));
                } else {
                    EXPECT_TRUE(Fp::equal(out_partial[lane], out_before[lane]));
                }
            }
        }
    }
}

TEST(FpAvx, Select8x1Masks)
{
    InsecureRng rng;

    for (int iter = 0; iter < 48; ++iter) {
        FpVec8 a{}, b{}, got{}, expected{};
        for (size_t lane = 0; lane < 8; ++lane) {
            rand_redc_fp(a[lane], rng);
            rand_redc_fp(b[lane], rng);
        }
        force_edge_lanes(a);
        force_edge_lanes(b);

        FpPack8x1 avx_a{}, avx_b{}, avx_r{};
        FpPack8x1::from_aos(avx_a, a.data());
        FpPack8x1::from_aos(avx_b, b.data());

        for (u16 mask = 0; mask < 256; ++mask) {
            FpPack8x1::select(avx_r, avx_a, avx_b, static_cast<u8>(mask));
            FpPack8x1::to_aos(got.data(), avx_r);
            apply_lane_mask_select(expected, a, b, static_cast<u8>(mask));
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
            }

            FpPack8x1 avx_alias_a{};
            FpPack8x1::from_aos(avx_alias_a, a.data());
            FpPack8x1::select(avx_alias_a, avx_alias_a, avx_b, static_cast<u8>(mask));
            FpPack8x1::to_aos(got.data(), avx_alias_a);
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
            }

            FpPack8x1 avx_alias_b{};
            FpPack8x1::from_aos(avx_alias_b, b.data());
            FpPack8x1::select(avx_alias_b, avx_a, avx_alias_b, static_cast<u8>(mask));
            FpPack8x1::to_aos(got.data(), avx_alias_b);
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
            }
        }
    }
}

TEST(FpAvx, SelectWide8x1Masks)
{
    InsecureRng rng;

    for (int iter = 0; iter < 40; ++iter) {
        FpVec8 a{}, b{}, c{}, d{};
        FpVec8 mul_ab{}, mul_cd{}, got{}, expected{};
        for (size_t lane = 0; lane < 8; ++lane) {
            rand_redc_fp(a[lane], rng);
            rand_redc_fp(b[lane], rng);
            rand_redc_fp(c[lane], rng);
            rand_redc_fp(d[lane], rng);
        }
        force_edge_lanes(a);
        force_edge_lanes(b);
        force_edge_lanes(c);
        force_edge_lanes(d);

        for (size_t lane = 0; lane < 8; ++lane) {
            Fp::mul(mul_ab[lane], a[lane], b[lane]);
            Fp::mul(mul_cd[lane], c[lane], d[lane]);
        }

        FpPack8x1 avx_a{}, avx_b{}, avx_c{}, avx_d{}, red{};
        UFpPack8x1 wide_ab{}, wide_cd{}, wide_sel{};

        FpPack8x1::from_aos(avx_a, a.data());
        FpPack8x1::from_aos(avx_b, b.data());
        FpPack8x1::from_aos(avx_c, c.data());
        FpPack8x1::from_aos(avx_d, d.data());
        FpPack8x1::mul(wide_ab, avx_a, avx_b);
        FpPack8x1::mul(wide_cd, avx_c, avx_d);

        for (u16 mask = 0; mask < 256; ++mask) {
            UFpPack8x1::select(wide_sel, wide_ab, wide_cd, static_cast<u8>(mask));
            UFpPack8x1::redc(red, wide_sel);
            FpPack8x1::to_aos(got.data(), red);
            apply_lane_mask_select(expected, mul_ab, mul_cd, static_cast<u8>(mask));
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
            }

            UFpPack8x1 wide_alias_a = wide_ab;
            UFpPack8x1::select(wide_alias_a, wide_alias_a, wide_cd, static_cast<u8>(mask));
            UFpPack8x1::redc(red, wide_alias_a);
            FpPack8x1::to_aos(got.data(), red);
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
            }

            UFpPack8x1 wide_alias_b = wide_cd;
            UFpPack8x1::select(wide_alias_b, wide_ab, wide_alias_b, static_cast<u8>(mask));
            UFpPack8x1::redc(red, wide_alias_b);
            FpPack8x1::to_aos(got.data(), red);
            for (size_t lane = 0; lane < 8; ++lane) {
                EXPECT_TRUE(Fp::equal(got[lane], expected[lane]));
            }
        }
    }
}

TEST(FpAvx, PartialTailMulRedc8x1w)
{
    InsecureRng rng;

    for (int iter = 0; iter < 192; ++iter) {
        FpVec8 a{}, b{}, out{};
        for (size_t lane = 0; lane < 8; ++lane) {
            rand_redc_fp(a[lane], rng);
            rand_redc_fp(b[lane], rng);
        }
        force_edge_lanes(a);
        force_edge_lanes(b);

        for (size_t n_active = 0; n_active <= 8; ++n_active) {
            FpPack8x1 avx_a{}, avx_b{}, red{};
            UFpPack8x1 wide{};
            FpPack8x1::from_aos_partial(avx_a, a.data(), n_active, Fp::one());
            FpPack8x1::from_aos_partial(avx_b, b.data(), n_active, Fp::one());

            FpPack8x1::mul(wide, avx_a, avx_b);
            UFpPack8x1::redc(red, wide);

            for (auto& e : out) {
                rand_redc_fp(e, rng);
            }
            const FpVec8 out_before = out;

            FpPack8x1::to_aos_partial(out.data(), n_active, red);
            for (size_t lane = 0; lane < 8; ++lane) {
                if (lane < n_active) {
                    Fp expected{};
                    Fp::mul(expected, a[lane], b[lane]);
                    EXPECT_TRUE(Fp::equal(out[lane], expected));
                } else {
                    EXPECT_TRUE(Fp::equal(out[lane], out_before[lane]));
                }
            }
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

        FpPack8x1 avx_a{}, avx_b{};
        pack_FpPack8x1(avx_a, a);
        pack_FpPack8x1(avx_b, b);

        UFpPack8x1 wide_alias{}, wide_v2{};
        FpPack8x1 red_alias{}, red_v2{};
        FpVec8 out_alias{}, out_v2{};

        FpPack8x1::mul(wide_alias, avx_a, avx_b);
        UFpPack8x1::redc(red_alias, wide_alias);
        unpack_FpPack8x1(out_alias, red_alias);

        FpPack8x1::mul(wide_v2, avx_a, avx_b);
        UFpPack8x1::redc(red_v2, wide_v2);
        unpack_FpPack8x1(out_v2, red_v2);

        for (size_t lane = 0; lane < 8; ++lane) {
            Fp expected{};
            Fp::mul(expected, a[lane], b[lane]);

            EXPECT_TRUE(Fp::equal(out_alias[lane], expected));
            EXPECT_TRUE(Fp::equal(out_v2[lane], expected));
            EXPECT_TRUE(Fp::equal(out_alias[lane], out_v2[lane]));
        }
    }
}

TEST(FpAvx, RedcCoreNormalize8x1w)
{
    InsecureRng rng;
    fp8x1_lfp_bounds_reset();

    for (int iter = 0; iter < 192; ++iter) {
        FpVec8 a{}, b{};
        for (size_t lane = 0; lane < 8; ++lane) {
            rand_redc_fp(a[lane], rng);
            rand_redc_fp(b[lane], rng);
        }
        force_edge_lanes(a);
        force_edge_lanes(b);

        FpPack8x1 avx_a{}, avx_b{};
        UFpPack8x1 wide{};
        FpPack8x1 red_direct{}, red_split{};
        LFpPack8x1 lazy{};
        FpVec8 out_direct{}, out_split{};

        pack_FpPack8x1(avx_a, a);
        pack_FpPack8x1(avx_b, b);

        FpPack8x1::mul(wide, avx_a, avx_b);
        UFpPack8x1::redc(red_direct, wide);
        UFpPack8x1::redc_core(lazy, wide);
        UFpPack8x1::normalize(red_split, lazy);

        unpack_FpPack8x1(out_direct, red_direct);
        unpack_FpPack8x1(out_split, red_split);
        for (size_t lane = 0; lane < 8; ++lane) {
            Fp expected{};
            Fp::mul(expected, a[lane], b[lane]);
            EXPECT_TRUE(Fp::equal(out_direct[lane], expected));
            EXPECT_TRUE(Fp::equal(out_split[lane], expected));
            EXPECT_TRUE(Fp::equal(out_direct[lane], out_split[lane]));
            EXPECT_EQ(Fp::is_redc(out_split[lane]), 1);
        }
    }

    const LFpBoundsStats8x1 stats = fp8x1_lfp_bounds_get();
#if defined(NDEBUG)
    EXPECT_EQ(stats.total_lanes, 0u);
    EXPECT_EQ(stats.ge_p_lanes, 0u);
    EXPECT_EQ(stats.max_limb, 0u);
#else
    EXPECT_GT(stats.total_lanes, 0u);
    EXPECT_LE(stats.ge_p_lanes, stats.total_lanes);
    EXPECT_LE(stats.max_limb, BMASK);
#endif
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

        FpPack4x2 avx4_a{}, avx4_b{}, avx4_red_alias{}, avx4_red_v1{};
        UFpPack4x2 avx4_wide_alias{}, avx4_wide_v1{};
        pack_FpPack4x2(avx4_a, a4);
        pack_FpPack4x2(avx4_b, b4);

        FpPack4x2::mul(avx4_wide_alias, avx4_a, avx4_b);
        UFpPack4x2::redc(avx4_red_alias, avx4_wide_alias);

        FpPack4x2::mul_v1(avx4_wide_v1, avx4_a, avx4_b);
        UFpPack4x2::redc(avx4_red_v1, avx4_wide_v1);

        FpVec4 out4_alias{}, out4_v1{};
        unpack_FpPack4x2(out4_alias, avx4_red_alias);
        unpack_FpPack4x2(out4_v1, avx4_red_v1);

        for (size_t lane = 0; lane < 4; ++lane) {
            Fp expected{};
            Fp::mul(expected, a4[lane], b4[lane]);
            EXPECT_TRUE(Fp::equal(out4_alias[lane], expected));
            EXPECT_TRUE(Fp::equal(out4_v1[lane], expected));
            EXPECT_TRUE(Fp::equal(out4_alias[lane], out4_v1[lane]));
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

        FpPack8x1 avx8_a{}, avx8_b{}, avx8_red{};
        UFpPack8x1 avx8_wide{};
        pack_FpPack8x1(avx8_a, a8);
        pack_FpPack8x1(avx8_b, b8);

        FpPack8x1::mul(avx8_wide, avx8_a, avx8_b);
        UFpPack8x1::redc(avx8_red, avx8_wide);
        unpack_FpPack8x1(out8, avx8_red);

        for (size_t lane = 0; lane < 4; ++lane) {
            EXPECT_TRUE(Fp::equal(out4_alias[lane], out8[lane]));
            EXPECT_TRUE(Fp::equal(out4_v1[lane], out8[lane]));
        }
    }
}

TEST(FpAvx, KWiseMulRedcKernels)
{
    InsecureRng rng;
    for (int iter = 0; iter < 32; ++iter) {
        run_kwise_mul_redc_matches_scalar<1>(rng);
        run_kwise_mul_redc_matches_scalar<2>(rng);
        run_kwise_mul_redc_matches_scalar<4>(rng);
        run_kwise_dispatch_matches_scalar(rng);
    }
}

#else

TEST(FpAvx, Disabled)
{
    SUCCEED();
}

#endif
