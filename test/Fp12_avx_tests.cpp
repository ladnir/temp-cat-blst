#include <gtest/gtest.h>

#include <array>

#include "InsecureRng.h"
#include "blst/fields/Fp12.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include "blst/fields/Fp12_avx.h"

using namespace blst;

namespace
{
using FpVec8 = std::array<Fp, 8>;
using FpVec4 = std::array<Fp, 4>;
using FpVec2 = std::array<Fp, 2>;
using Fp2Vec8 = std::array<Fp2, 8>;
using Fp2Vec4 = std::array<Fp2, 4>;
using Fp2Vec2 = std::array<Fp2, 2>;
using Fp2Vec6 = std::array<Fp2, 6>;

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

void fill_edge_fp12(Fp12& v, int idx)
{
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            fill_edge_fp2(v[i][j], idx + static_cast<int>(3 * i + j));
        }
    }
}

bool fp2_equal(const Fp2& a, const Fp2& b)
{
    return Fp::equal(a[0], b[0]) && Fp::equal(a[1], b[1]);
}

bool fp12_equal(const Fp12& a, const Fp12& b)
{
    return Fp12::equal(a, b) != 0;
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

template <size_t N0, size_t N1>
bool contains_any(const std::array<Fp2, N0>& a, const std::array<Fp2, N1>& b, const Fp2& v)
{
    return contains_fp2(a, v) || contains_fp2(b, v);
}

template <size_t N0, size_t N1, size_t N2>
bool contains_any(const std::array<Fp2, N0>& a, const std::array<Fp2, N1>& b, const std::array<Fp2, N2>& c, const Fp2& v)
{
    return contains_fp2(a, v) || contains_fp2(b, v) || contains_fp2(c, v);
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

void pack_fp4_2x2x2x1w(Fp4Pack2x2x2x1& out, const Fp2Vec4& in)
{
    pack_Fp2Pack4x2x1(out, in);
}

void unpack_fp4_2x2x2x1w(Fp2Vec4& out, const Fp4Pack2x2x2x1& in)
{
    unpack_Fp2Pack4x2x1(out, in);
}

void pack_fp4_1x2x2x2w(Fp4Pack1x2x2x2& out, const Fp2Vec2& in)
{
    pack_Fp2Pack2x2x2(out, in);
}

void unpack_fp4_1x2x2x2w(Fp2Vec2& out, const Fp4Pack1x2x2x2& in)
{
    unpack_Fp2Pack2x2x2(out, in);
}

Fp2Vec6 flatten_fp12(const Fp12& f)
{
    Fp2Vec6 out{};
    out[0] = f[0][0];
    out[1] = f[0][1];
    out[2] = f[0][2];
    out[3] = f[1][0];
    out[4] = f[1][1];
    out[5] = f[1][2];
    return out;
}

} // namespace

TEST(Fp12Avx, MulVariantsBroadcastAgainstScalar)
{
    InsecureRng rng;

    for (int iter = 0; iter < 144; ++iter) {
        Fp2 c0{}, c1{}, c2{};
        if (iter < 16) {
            fill_edge_fp2(c0, iter + 0);
            fill_edge_fp2(c1, iter + 1);
            fill_edge_fp2(c2, iter + 2);
        } else {
            rand_redc_fp2(c0, rng);
            rand_redc_fp2(c1, rng);
            rand_redc_fp2(c2, rng);
        }

        Fp12 a{}, b{}, expected{};
        for (size_t i = 0; i < 2; ++i) {
            a[i][0] = c0;
            a[i][1] = c1;
            a[i][2] = c2;
            b[i][0] = c0;
            b[i][1] = c1;
            b[i][2] = c2;
        }
        Fp12::mul(expected, a, b);
        const auto exp = flatten_fp12(expected);

        Fp12MulInput8x1x1 ab8{};
        Fp12MulInput4x2x1 ab4{};
        pack_Fp2Pack8x1x1(ab8[0], filled_fp2<8>(c0));
        pack_Fp2Pack8x1x1(ab8[1], filled_fp2<8>(c1));
        pack_Fp2Pack8x1x1(ab8[2], filled_fp2<8>(c2));
        pack_Fp2Pack4x2x1(ab4[0], filled_fp2<4>(c0));
        pack_Fp2Pack4x2x1(ab4[1], filled_fp2<4>(c1));
        pack_Fp2Pack4x2x1(ab4[2], filled_fp2<4>(c2));

        Fp2Pack4x2x1 r01_v1{}, r2_v1{};
        Fp2Vec4 out01_v1{}, out2_v1{};
        Fp12MulInput8x1x1::mul_v1(r01_v1, r2_v1, ab8);
        unpack_Fp2Pack4x2x1(out01_v1, r01_v1);
        unpack_Fp2Pack4x2x1(out2_v1, r2_v1);

        Fp2Pack4x2x1 r01_v2{};
        Fp2Pack2x2x2 r2_v2{};
        Fp2Vec4 out01_v2{};
        Fp2Vec2 out2_v2{};
        Fp12MulInput8x1x1::mul_v2(r01_v2, r2_v2, ab8);
        unpack_Fp2Pack4x2x1(out01_v2, r01_v2);
        unpack_Fp2Pack2x2x2(out2_v2, r2_v2);

        Fp2Pack2x2x2 r001_v3{}, r101_v3{}, r2_v3{};
        Fp2Vec2 out001_v3{}, out101_v3{}, out2_v3{};
        Fp12MulInput4x2x1::mul_v3(r001_v3, r101_v3, r2_v3, ab4);
        unpack_Fp2Pack2x2x2(out001_v3, r001_v3);
        unpack_Fp2Pack2x2x2(out101_v3, r101_v3);
        unpack_Fp2Pack2x2x2(out2_v3, r2_v3);

        Fp2Pack4x2x1 r0_v4{};
        Fp2Pack2x2x2 r101_v4{}, r12_v4{};
        Fp2Vec4 out0_v4{};
        Fp2Vec2 out101_v4{}, out12_v4{};
        Fp12MulInput4x2x1::mul_v4(r0_v4, r101_v4, r12_v4, ab4);
        unpack_Fp2Pack4x2x1(out0_v4, r0_v4);
        unpack_Fp2Pack2x2x2(out101_v4, r101_v4);
        unpack_Fp2Pack2x2x2(out12_v4, r12_v4);

        for (const auto& e : exp) {
            EXPECT_TRUE(contains_any(out01_v1, out2_v1, e));
            EXPECT_TRUE(contains_any(out01_v2, out2_v2, e));
            EXPECT_TRUE(contains_any(out001_v3, out101_v3, out2_v3, e));
            EXPECT_TRUE(contains_any(out0_v4, out101_v4, out12_v4, e));
        }
    }
}

TEST(Fp12Avx, SqrBroadcastDeterministic)
{
    InsecureRng rng;

    for (int iter = 0; iter < 128; ++iter) {
        Fp2 v{};
        if (iter < 16) {
            fill_edge_fp2(v, iter);
        } else {
            rand_redc_fp2(v, rng);
        }

        Fp12SquareInput4x2x1 in{};
        Fp2Pack4x2x1 r0{}, r1{};
        Fp2Pack4x2x1 r0b{}, r1b{};
        Fp2Vec4 out0{}, out1{};
        Fp2Vec4 out0b{}, out1b{};
        pack_Fp2Pack4x2x1(in.a0, filled_fp2<4>(v));
        pack_Fp2Pack4x2x1(in.a1, filled_fp2<4>(v));
        Fp12SquareInput4x2x1::sqr_v1(r0, r1, in);
        Fp12SquareInput4x2x1::sqr_v1(r0b, r1b, in);
        unpack_Fp2Pack4x2x1(out0, r0);
        unpack_Fp2Pack4x2x1(out1, r1);
        unpack_Fp2Pack4x2x1(out0b, r0b);
        unpack_Fp2Pack4x2x1(out1b, r1b);

        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(fp2_equal(out0[i], out0b[i]));
            EXPECT_TRUE(fp2_equal(out1[i], out1b[i]));
        }
    }
}

TEST(Fp12Avx, CyclotomicVariantsAndCompressed)
{
    InsecureRng rng;

    for (int iter = 0; iter < 128; ++iter) {
        Fp2 va{}, vbc{};
        if (iter < 16) {
            fill_edge_fp2(va, iter + 0);
            fill_edge_fp2(vbc, iter + 1);
        } else {
            rand_redc_fp2(va, rng);
            rand_redc_fp2(vbc, rng);
        }

        Fp12CyclotomicState in{}, out1{}, out2{};
        Fp12CyclotomicBc rbc_comp{};
        Fp2Vec2 out_ra1{}, out_ra2{};
        Fp2Vec4 out_rbc1{}, out_rbc2{}, out_rbc_comp{};

        pack_fp4_1x2x2x2w(in.a, filled_fp2<2>(va));
        pack_fp4_2x2x2x1w(in.bc, filled_fp2<4>(vbc));

        Fp12CyclotomicState::sqr_v1(out1, in);
        Fp12CyclotomicState::sqr_v2(out2, in);
        Fp12CyclotomicBc::compressed_sqr_v1(rbc_comp, in.bc);

        unpack_fp4_1x2x2x2w(out_ra1, out1.a);
        unpack_fp4_1x2x2x2w(out_ra2, out2.a);
        unpack_fp4_2x2x2x1w(out_rbc1, out1.bc);
        unpack_fp4_2x2x2x1w(out_rbc2, out2.bc);
        unpack_fp4_2x2x2x1w(out_rbc_comp, rbc_comp);

        for (size_t i = 0; i < 2; ++i) {
            EXPECT_TRUE(fp2_equal(out_ra1[i], out_ra2[i]));
        }
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(fp2_equal(out_rbc1[i], out_rbc2[i]));
            EXPECT_TRUE(fp2_equal(out_rbc1[i], out_rbc_comp[i]));
        }
    }
}

TEST(Fp12Avx, SparseMulVariantsBroadcast)
{
    InsecureRng rng;

    for (int iter = 0; iter < 128; ++iter) {
        Fp2 av{}, bv{};
        if (iter < 16) {
            fill_edge_fp2(av, iter + 0);
            fill_edge_fp2(bv, iter + 1);
        } else {
            rand_redc_fp2(av, rng);
            rand_redc_fp2(bv, rng);
        }

        Fp12SparseInput4x2x1 in4{};
        Fp12SparseInput8x1x1 in8{};
        Fp2Pack4x2x1 r0_v1{};
        Fp2Pack2x2x2 r1_v1{}, r1_v2{};
        Fp2Pack4x2x1 r0_v2{};
        Fp2Vec4 out0_v1{}, out0_v2{};
        Fp2Vec2 out1_v1{}, out1_v2{};

        pack_Fp2Pack4x2x1(in4.a01, filled_fp2<4>(av));
        pack_Fp2Pack4x2x1(in4.a2, filled_fp2<4>(av));
        pack_Fp2Pack4x2x1(in4.b01, filled_fp2<4>(bv));
        pack_Fp2Pack4x2x1(in4.b4, filled_fp2<4>(bv));
        pack_Fp2Pack8x1x1(in8.a, filled_fp2<8>(av));
        pack_Fp2Pack8x1x1(in8.b, filled_fp2<8>(bv));

        Fp12SparseInput4x2x1::mul_v1(r0_v1, r1_v1, in4);
        Fp12SparseInput8x1x1::mul_v2(r0_v2, r1_v2, in8);

        unpack_Fp2Pack4x2x1(out0_v1, r0_v1);
        unpack_Fp2Pack2x2x2(out1_v1, r1_v1);
        unpack_Fp2Pack4x2x1(out0_v2, r0_v2);
        unpack_Fp2Pack2x2x2(out1_v2, r1_v2);

        Fp12 a_scalar{}, expected{};
        Fp6 sparse{};
        for (size_t i = 0; i < 2; ++i) {
            a_scalar[i][0] = av;
            a_scalar[i][1] = av;
            a_scalar[i][2] = av;
        }
        sparse[0] = bv;
        sparse[1] = bv;
        sparse[2] = bv;
        Fp12::mul_by_xy00z0(expected, a_scalar, sparse);
        const auto exp = flatten_fp12(expected);

        for (const auto& e : exp) {
            EXPECT_TRUE(contains_any(out0_v1, out1_v1, e));
            EXPECT_TRUE(contains_any(out0_v2, out1_v2, e));
        }
    }
}

TEST(Fp12Avx, MulVectorMatchesScalar)
{
    InsecureRng rng;

    for (int iter = 0; iter < 96; ++iter) {
        Fp12 a{}, b{}, expected{}, got{};
        if (iter < 16) {
            fill_edge_fp12(a, iter + 0);
            fill_edge_fp12(b, iter + 7);
        } else {
            rand_redc_fp12(a, rng);
            rand_redc_fp12(b, rng);
        }

        Fp12::mul(expected, a, b);
        mul_fp12_vector(got, a, b);
        EXPECT_TRUE(fp12_equal(got, expected));
    }
}

TEST(Fp12Avx, SqrVectorMatchesScalar)
{
    InsecureRng rng;

    for (int iter = 0; iter < 96; ++iter) {
        Fp12 a{}, expected{}, got{};
        if (iter < 16) {
            fill_edge_fp12(a, iter);
        } else {
            rand_redc_fp12(a, rng);
        }

        Fp12::sqr(expected, a);
        sqr_fp12_vector(got, a);
        EXPECT_TRUE(fp12_equal(got, expected));
    }
}

TEST(Fp12Avx, SparseMulVectorMatchesScalar)
{
    InsecureRng rng;

    for (int iter = 0; iter < 96; ++iter) {
        Fp12 a{}, expected{}, got{};
        Fp6 sparse{};
        if (iter < 16) {
            fill_edge_fp12(a, iter);
            for (size_t j = 0; j < 3; ++j) {
                fill_edge_fp2(sparse[j], iter + static_cast<int>(j) + 3);
            }
        } else {
            rand_redc_fp12(a, rng);
            for (size_t j = 0; j < 3; ++j) {
                rand_redc_fp2(sparse[j], rng);
            }
        }

        Fp12::mul_by_xy00z0(expected, a, sparse);
        mul_by_xy00z0_fp12_vector(got, a, sparse);
        EXPECT_TRUE(fp12_equal(got, expected));
    }
}

TEST(Fp12Avx, ExplicitAvxApisMatchScalarApis)
{
    InsecureRng rng;

    for (int iter = 0; iter < 96; ++iter) {
        Fp12 a{}, b{};
        Fp12 mul_scalar{}, mul_avx{};
        Fp12 sqr_scalar{}, sqr_avx{};
        Fp12 sparse_scalar{}, sparse_avx{};
        Fp12 cyclo_seed{}, cyclo{};
        Fp12 cyc_scalar{}, cyc_avx{};
        Fp6 sparse{};

        if (iter < 16) {
            fill_edge_fp12(a, iter + 0);
            fill_edge_fp12(b, iter + 7);
            for (size_t j = 0; j < 3; ++j) {
                fill_edge_fp2(sparse[j], iter + static_cast<int>(j) + 11);
            }
        } else {
            rand_redc_fp12(a, rng);
            rand_redc_fp12(b, rng);
            for (size_t j = 0; j < 3; ++j) {
                rand_redc_fp2(sparse[j], rng);
            }
        }

        Fp12::mul_scalar(mul_scalar, a, b);
        Fp12::mul_avx(mul_avx, a, b);
        EXPECT_TRUE(fp12_equal(mul_avx, mul_scalar));

        Fp12::sqr_scalar(sqr_scalar, a);
        Fp12::sqr_avx(sqr_avx, a);
        EXPECT_TRUE(fp12_equal(sqr_avx, sqr_scalar));

        Fp12::mul_by_xy00z0_scalar(sparse_scalar, a, sparse);
        Fp12::mul_by_xy00z0_avx(sparse_avx, a, sparse);
        EXPECT_TRUE(fp12_equal(sparse_avx, sparse_scalar));

        cyclo_seed = a;
        Fp12::final_exp(cyclo, cyclo_seed);
        ASSERT_TRUE(Fp12::is_cyclotomic(cyclo));
        Fp12::cyclotomic_sqr_scalar(cyc_scalar, cyclo);
        Fp12::cyclotomic_sqr_avx(cyc_avx, cyclo);
        EXPECT_TRUE(fp12_equal(cyc_avx, cyc_scalar));
    }
}

TEST(Fp12Avx, CyclotomicSqrVectorMatchesScalar)
{
    InsecureRng rng;

    for (int iter = 0; iter < 48; ++iter) {
        Fp12 seed{}, cyclo{}, expected{}, got{};
        rand_redc_fp12(seed, rng);
        Fp12::final_exp(cyclo, seed);

        Fp12::cyclotomic_sqr(expected, cyclo);
        cyclotomic_sqr_fp12_vector(got, cyclo);
        EXPECT_TRUE(fp12_equal(got, expected));
    }
}

#else

TEST(Fp12Avx, Disabled)
{
    SUCCEED();
}

#endif
