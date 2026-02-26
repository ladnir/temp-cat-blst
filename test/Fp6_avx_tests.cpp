#include <gtest/gtest.h>

#include <array>

#include "InsecureRng.h"
#include "blst/fields/Fp2.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include "blst/fields/Fp6_avx.h"

using namespace blst;

namespace
{
using FpVec8 = std::array<Fp, 8>;
using FpVec4 = std::array<Fp, 4>;
using FpVec2 = std::array<Fp, 2>;
using Fp2Vec8 = std::array<Fp2, 8>;
using Fp2Vec4 = std::array<Fp2, 4>;
using Fp2Vec2 = std::array<Fp2, 2>;

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

void fill_edge_fp2(Fp2& v, int idx)
{
    Fp two;
    Fp p_minus_one;
    Fp::add(two, Fp::one(), Fp::one());
    Fp::sub(p_minus_one, Fp::zero(), Fp::one());

    const Fp* cycle[4] = {&Fp::zero(), &Fp::one(), &two, &p_minus_one};
    v[0] = *cycle[(idx + 0) & 3];
    v[1] = *cycle[(idx + 1) & 3];
}

bool fp2_equal(const Fp2& a, const Fp2& b)
{
    return Fp::equal(a[0], b[0]) && Fp::equal(a[1], b[1]);
}

template <size_t N>
std::array<Fp2, N> filled_fp2(const Fp2& v)
{
    std::array<Fp2, N> out{};
    for (auto& x : out) {
        x = v;
    }
    return out;
}

template <size_t N>
bool contains_fp2(const std::array<Fp2, N>& lanes, const Fp2& value)
{
    for (const auto& lane : lanes) {
        if (fp2_equal(lane, value)) {
            return true;
        }
    }
    return false;
}

void expect_fp2_eq(const Fp2& got, const Fp2& expected)
{
    EXPECT_TRUE(fp2_equal(got, expected));
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

void redc_UFp2Pack8x1x1_local(Fp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
{
    UFpPack8x1::redc(r[0], a[0]);
    UFpPack8x1::redc(r[1], a[1]);
}

void mul_fp6_scalar(Fp2& r0, Fp2& r1, Fp2& r2, const Fp2& a0, const Fp2& a1, const Fp2& a2, const Fp2& b0, const Fp2& b1, const Fp2& b2)
{
    Fp2 t00{}, t01{}, t02{}, t10{}, t11{}, t20{};
    Fp2 t21{}, t22{}, tmp{};

    Fp2::mul(t00, a0, b0);
    Fp2::mul(t01, a1, b2);
    Fp2::mul(t02, a2, b1);
    Fp2::add(t10, t01, t02);
    Fp2::mul_by_u_plus_1(t10, t10);
    Fp2::add(r0, t00, t10);

    Fp2::mul(t11, a0, b1);
    Fp2::mul(t20, a1, b0);
    Fp2::mul(t21, a2, b2);
    Fp2::mul_by_u_plus_1(t21, t21);
    Fp2::add(tmp, t11, t20);
    Fp2::add(r1, tmp, t21);

    Fp2::mul(t22, a0, b2);
    Fp2::mul(t01, a1, b1);
    Fp2::mul(t02, a2, b0);
    Fp2::add(tmp, t22, t01);
    Fp2::add(r2, tmp, t02);
}

void mul_sparse_xy00z0_scalar(
    Fp2& r0,
    Fp2& r1,
    Fp2& r2,
    Fp2& r3,
    Fp2& r4,
    Fp2& r5,
    const Fp2& a0,
    const Fp2& a1,
    const Fp2& a2,
    const Fp2& b0,
    const Fp2& b1,
    const Fp2& b4)
{
    Fp2 t0{}, t1{}, t2{};

    Fp2::mul(t0, a2, b1);
    Fp2::mul_by_u_plus_1(t0, t0);
    Fp2::mul(t1, a0, b0);
    Fp2::add(r0, t0, t1);

    Fp2::mul(t0, a0, b1);
    Fp2::mul(t1, a1, b0);
    Fp2::add(r1, t0, t1);

    Fp2::mul(t0, a2, b0);
    Fp2::mul(t1, a1, b1);
    Fp2::add(r2, t0, t1);

    Fp2::mul(t2, a2, b4);
    Fp2::mul_by_u_plus_1(r3, t2);

    Fp2::mul(r4, a0, b4);
    Fp2::mul(r5, a1, b4);
}

} // namespace

TEST(Fp6Avx, Mul2x2x2x1w)
{
    InsecureRng rng;

    for (int iter = 0; iter < 192; ++iter) {
        Fp2 a0{}, a1{}, a2{};
        if (iter < 16) {
            fill_edge_fp2(a0, iter + 0);
            fill_edge_fp2(a1, iter + 1);
            fill_edge_fp2(a2, iter + 2);
        } else {
            rand_redc_fp2(a0, rng);
            rand_redc_fp2(a1, rng);
            rand_redc_fp2(a2, rng);
        }

        Fp2Pack4x2x1 ab0{}, ab1{}, ab2{}, red01{}, red2{};
        UFp2Pack4x2x1 wide01{}, wide2{};
        Fp2Vec4 out01{}, out2{};

        pack_Fp2Pack4x2x1(ab0, filled_fp2<4>(a0));
        pack_Fp2Pack4x2x1(ab1, filled_fp2<4>(a1));
        pack_Fp2Pack4x2x1(ab2, filled_fp2<4>(a2));

        Fp6Pack2x2x2x1::mul(wide01, wide2, ab0, ab1, ab2);
        UFp2Pack4x2x1::redc(red01, wide01);
        UFp2Pack4x2x1::redc(red2, wide2);
        unpack_Fp2Pack4x2x1(out01, red01);
        unpack_Fp2Pack4x2x1(out2, red2);

        Fp2 e0{}, e1{}, e2{};
        mul_fp6_scalar(e0, e1, e2, a0, a1, a2, a0, a1, a2);
        EXPECT_TRUE(contains_fp2(out01, e0));
        EXPECT_TRUE(contains_fp2(out01, e1));
        EXPECT_TRUE(contains_fp2(out2, e2));
    }
}

TEST(Fp6Avx, Mul1x2x2x2w)
{
    InsecureRng rng;

    for (int iter = 0; iter < 192; ++iter) {
        Fp2 a0{}, a1{}, a2{}, b0{}, b1{}, b2{};
        if (iter < 16) {
            fill_edge_fp2(a0, iter + 0);
            fill_edge_fp2(a1, iter + 1);
            fill_edge_fp2(a2, iter + 2);
            fill_edge_fp2(b0, iter + 3);
            fill_edge_fp2(b1, iter + 4);
            fill_edge_fp2(b2, iter + 5);
        } else {
            rand_redc_fp2(a0, rng);
            rand_redc_fp2(a1, rng);
            rand_redc_fp2(a2, rng);
            rand_redc_fp2(b0, rng);
            rand_redc_fp2(b1, rng);
            rand_redc_fp2(b2, rng);
        }

        Fp2Pack2x2x2 ab0{}, ab1{}, ab2{}, red01{}, red2{};
        UFp2Pack2x2x2 wide01{}, wide2{};
        Fp2Vec2 out01{}, out2{};

        pack_Fp2Pack2x2x2(ab0, {a0, b0});
        pack_Fp2Pack2x2x2(ab1, {a1, b1});
        pack_Fp2Pack2x2x2(ab2, {a2, b2});

        Fp6Pack1x2x2x2::mul(wide01, wide2, ab0, ab1, ab2);
        UFp2Pack2x2x2::redc(red01, wide01);
        UFp2Pack2x2x2::redc(red2, wide2);
        unpack_Fp2Pack2x2x2(out01, red01);
        unpack_Fp2Pack2x2x2(out2, red2);

        Fp2 e0{}, e1{}, e2{};
        mul_fp6_scalar(e0, e1, e2, a0, a1, a2, b0, b1, b2);
        expect_fp2_eq(out01[0], e0);
        expect_fp2_eq(out01[1], e1);
        expect_fp2_eq(out2[0], e2);
    }
}

TEST(Fp6Avx, Mul4x2x1x1w)
{
    InsecureRng rng;

    for (int iter = 0; iter < 160; ++iter) {
        Fp2 a0{}, a1{}, a2{};
        if (iter < 16) {
            fill_edge_fp2(a0, iter + 0);
            fill_edge_fp2(a1, iter + 1);
            fill_edge_fp2(a2, iter + 2);
        } else {
            rand_redc_fp2(a0, rng);
            rand_redc_fp2(a1, rng);
            rand_redc_fp2(a2, rng);
        }

        Fp2Pack8x1x1 ab0{}, ab1{}, ab2{}, red01{}, red2{};
        UFp2Pack8x1x1 wide01{}, wide2{};
        Fp2Vec8 out01{}, out2{};

        pack_Fp2Pack8x1x1(ab0, filled_fp2<8>(a0));
        pack_Fp2Pack8x1x1(ab1, filled_fp2<8>(a1));
        pack_Fp2Pack8x1x1(ab2, filled_fp2<8>(a2));

        Fp6Pack4x2x1x1::mul(wide01, wide2, ab0, ab1, ab2);
        redc_UFp2Pack8x1x1_local(red01, wide01);
        redc_UFp2Pack8x1x1_local(red2, wide2);
        unpack_Fp2Pack8x1x1(out01, red01);
        unpack_Fp2Pack8x1x1(out2, red2);

        Fp2 e0{}, e1{}, e2{};
        mul_fp6_scalar(e0, e1, e2, a0, a1, a2, a0, a1, a2);
        EXPECT_TRUE(contains_fp2(out01, e0));
        EXPECT_TRUE(contains_fp2(out01, e1));
        EXPECT_TRUE(contains_fp2(out2, e2));
    }
}

TEST(Fp6Avx, Sparse2x2x2x1w)
{
    InsecureRng rng;

    for (int iter = 0; iter < 192; ++iter) {
        Fp2 a01v{}, a2v{}, b01v{}, b4v{};
        if (iter < 16) {
            fill_edge_fp2(a01v, iter + 0);
            fill_edge_fp2(a2v, iter + 1);
            fill_edge_fp2(b01v, iter + 2);
            fill_edge_fp2(b4v, iter + 3);
        } else {
            rand_redc_fp2(a01v, rng);
            rand_redc_fp2(a2v, rng);
            rand_redc_fp2(b01v, rng);
            rand_redc_fp2(b4v, rng);
        }

        Fp2Pack4x2x1 a01{}, a2p{}, b01{}, b4p{}, red01{}, red2{}, red3{}, red45{};
        UFp2Pack4x2x1 wide01{}, wide2{}, wide3{}, wide45{};
        Fp2Vec4 out01{}, out2{}, out3{}, out45{};

        pack_Fp2Pack4x2x1(a01, filled_fp2<4>(a01v));
        pack_Fp2Pack4x2x1(a2p, filled_fp2<4>(a2v));
        pack_Fp2Pack4x2x1(b01, filled_fp2<4>(b01v));
        pack_Fp2Pack4x2x1(b4p, filled_fp2<4>(b4v));

        Fp6Pack2x2x2x1::mul_by_xy00z0(wide01, wide2, wide3, wide45, a01, a2p, b01, b4p);
        UFp2Pack4x2x1::redc(red01, wide01);
        UFp2Pack4x2x1::redc(red2, wide2);
        UFp2Pack4x2x1::redc(red3, wide3);
        UFp2Pack4x2x1::redc(red45, wide45);
        unpack_Fp2Pack4x2x1(out01, red01);
        unpack_Fp2Pack4x2x1(out2, red2);
        unpack_Fp2Pack4x2x1(out3, red3);
        unpack_Fp2Pack4x2x1(out45, red45);

        Fp2 e0{}, e1{}, e2{}, e3{}, e4{}, e5{}, e3low{};
        mul_sparse_xy00z0_scalar(e0, e1, e2, e3, e4, e5, a01v, a01v, a2v, b01v, b01v, b4v);
        Fp2::mul(e3low, a2v, b01v);
        Fp2::mul_by_u_plus_1(e3low, e3low);

        EXPECT_TRUE(contains_fp2(out01, e0));
        EXPECT_TRUE(contains_fp2(out01, e1));
        EXPECT_TRUE(contains_fp2(out2, e2));
        EXPECT_TRUE(contains_fp2(out3, e3low));
        EXPECT_TRUE(contains_fp2(out3, e3));
        EXPECT_TRUE(contains_fp2(out45, e4));
        EXPECT_TRUE(contains_fp2(out45, e5));
    }
}

TEST(Fp6Avx, Sparse2x4x1x1w)
{
    InsecureRng rng;

    for (int iter = 0; iter < 160; ++iter) {
        Fp2 av{}, bv{};
        if (iter < 16) {
            fill_edge_fp2(av, iter + 0);
            fill_edge_fp2(bv, iter + 1);
        } else {
            rand_redc_fp2(av, rng);
            rand_redc_fp2(bv, rng);
        }

        Fp2Pack8x1x1 a{}, b{}, red01_4{}, red23_4{}, red45_4{};
        UFp2Pack8x1x1 wide01_4{}, wide23_4{}, wide45_4{};
        Fp2Vec8 out01_4{}, out23_4{}, out45_4{};

        pack_Fp2Pack8x1x1(a, filled_fp2<8>(av));
        pack_Fp2Pack8x1x1(b, filled_fp2<8>(bv));

        Fp6Pack4x2x1x1::mul_by_xy00z0(wide01_4, wide23_4, wide45_4, a, b);
        redc_UFp2Pack8x1x1_local(red01_4, wide01_4);
        redc_UFp2Pack8x1x1_local(red23_4, wide23_4);
        redc_UFp2Pack8x1x1_local(red45_4, wide45_4);
        unpack_Fp2Pack8x1x1(out01_4, red01_4);
        unpack_Fp2Pack8x1x1(out23_4, red23_4);
        unpack_Fp2Pack8x1x1(out45_4, red45_4);

        Fp2 r0{}, r1{}, r2{}, r3{}, r4{}, r5{}, r3low{}, raw_a2b1{}, raw_a2b4{};
        mul_sparse_xy00z0_scalar(r0, r1, r2, r3, r4, r5, av, av, av, bv, bv, bv);
        Fp2::mul(raw_a2b1, av, bv);
        Fp2::mul(raw_a2b4, av, bv);
        Fp2::mul_by_u_plus_1(r3low, raw_a2b1);

        EXPECT_TRUE(contains_fp2(out01_4, r0));
        EXPECT_TRUE(contains_fp2(out01_4, r1));
        EXPECT_TRUE(contains_fp2(out23_4, r2));
        EXPECT_TRUE(contains_fp2(out23_4, r3));
        EXPECT_TRUE(contains_fp2(out23_4, r3low));
        EXPECT_TRUE(contains_fp2(out45_4, r4));
        EXPECT_TRUE(contains_fp2(out45_4, r5));
        EXPECT_TRUE(contains_fp2(out45_4, raw_a2b1));
        EXPECT_TRUE(contains_fp2(out45_4, raw_a2b4));
    }
}

#else

TEST(Fp6Avx, Disabled)
{
    SUCCEED();
}

#endif
