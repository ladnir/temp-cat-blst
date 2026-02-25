#include "blst/curve/G2_avx.h"
#include "blst/backend/avx512/avx_layout.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include <array>

namespace blst
{
namespace
{
    inline void store_u64_lanes(std::array<u64, 8>& lanes, const cc512i v)
    {
        _cc512_storeu_si512(lanes.data(), v);
    }

#if defined(NDEBUG)
#define BLST_ASSERT_NO_OVERLAP_PTR(a, an, b, bn) ((void)0)
#else
    inline bool g2_ptr_ranges_overlap(const void* a, std::size_t a_bytes, const void* b, std::size_t b_bytes)
    {
        const auto* ap = static_cast<const u8*>(a);
        const auto* bp = static_cast<const u8*>(b);
        return (ap < (bp + b_bytes)) && (bp < (ap + a_bytes));
    }
#define BLST_ASSERT_NO_OVERLAP_PTR(a, an, b, bn) \
    BLST_ASSERT(!g2_ptr_ranges_overlap((a), (an), (b), (bn)))
#endif
} // namespace

void line_dbl_vec_v1(
    fp2_2x2x2w& l0Y3,
    fp2_2x2x2w& l1,
    fp2_2x2x2w& l2,
    fp2_2x2x2w& X3,
    fp2_2x2x2w& Z3,
    const fp2_2x2x2w& X1Y1,
    const fp2_2x2x2w& Z1)
{
    fp2_2x2x2w AB, C, D, E, FZZ;
    fp2_2x2x2w t0, t1, t2, t3, t4, t5;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 2, 1, 0, 3, 2, 1, 0);

    sqr_fp2_2x2x2w(AB, X1Y1);
    perm_var_hl(t3.data(), X1Y1.data(), m0);
    blend_0x0F_hl(t0.data(), t3.data(), AB.data());
    add_fp2_2x2x2w(t0, t0, AB);
    perm_var_hl(t1.data(), AB.data(), m0);
    blend_0x0F_hl(t2.data(), t0.data(), t1.data());
    sqr_fp2_2x2x2w(C, t2);
    blend_0x0F_hl(t2.data(), C.data(), t0.data());
    blend_0x0F_hl(t0.data(), t1.data(), AB.data());
    sa_fp2_2x2x2w(E, t2, t0);
    perm_var_hl(t5.data(), C.data(), m0);
    blend_0x0F_hl(t2.data(), t5.data(), X1Y1.data());
    sa_fp2_2x2x2w(t0, E, t2);
    blend_0x0F_hl(t2.data(), Z1.data(), E.data());
    sqr_fp2_2x2x2w(FZZ, t2);
    blend_0x0F_hl(t3.data(), t0.data(), t3.data());
    blend_0x0F_hl(t2.data(), t0.data(), Z1.data());
    add_fp2_2x2x2w(D, t3, t2);
    perm_var_hl(t0.data(), t0.data(), m0);
    blend_0x0F_hl(t0.data(), t0.data(), D.data());
    sqr_fp2_2x2x2w(t0, t0);
    perm_var_hl(t2.data(), FZZ.data(), m0);
    blend_0x0F_hl(t3.data(), t2.data(), t0.data());
    blend_0x0F_hl(t1.data(), D.data(), t1.data());
    sub_fp2_2x2x2w(t3, t3, t1);
    perm_var_hl(t0.data(), t0.data(), m0);
    blend_0x0F_hl(t0.data(), t3.data(), t0.data());
    blend_0x0F_hl(t4.data(), D.data(), AB.data());
    sub_fp2_2x2x2w(X3, t0, t4);
    perm_var_hl(t3.data(), t3.data(), m0);
    blend_0x0F_hl(t3.data(), t3.data(), C.data());
    blend_0x0F_hl(t4.data(), FZZ.data(), C.data());
    sa_fp2_2x2x2w(Z3, t3, t4);
    perm_var_hl(t4.data(), Z3.data(), m0);
    blend_0x0F_hl(t1.data(), t4.data(), t1.data());
    add_fp2_2x2x2w(t1, t1, t1);
    add_fp2_2x2x2w(t1, t1, t1);
    blend_0x0F_hl(t0.data(), D.data(), X3.data());
    blend_0x0F_hl(t3.data(), X3.data(), FZZ.data());
    sub_fp2_2x2x2w(t0, t0, t3);
    blend_0x0F_hl(t2.data(), t0.data(), t2.data());
    perm_var_hl(t3.data(), E.data(), m1);
    mul_fp2_2x2x2w(l1, t2, t3);
    blend_0x0F_hl(t0.data(), l1.data(), t0.data());
    sub_fp2_2x2x2w(l0Y3, t0, t1);
    mul_fp2_2x2x2w(l2, Z3, FZZ);
}

void line_dbl_vec_v2(
    fp2_4x2x1w& l0,
    fp2_4x2x1w& l12,
    fp2_4x2x1w& X3,
    fp2_4x2x1w& Y3,
    fp2_4x2x1w& Z3,
    const fp2_4x2x1w& X1Y1Z1)
{
    fp2_4x2x1w ABZZ, CF, E;
    fp2_4x2x1w t0, t1, t2, t3;
    const cc512i m0 = VSET(7, 6, 5, 4, 7, 6, 5, 4);
    const cc512i m1 = VSET(5, 4, 5, 4, 5, 4, 5, 4);
    const cc512i m2 = VSET(1, 0, 0, 0, 3, 2, 0, 0);
    const cc512i m3 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m4 = VSET(5, 4, 7, 6, 5, 4, 7, 6);
    const cc512i m5 = VSET(1, 0, 7, 6, 7, 6, 7, 6);
    const cc512i m6 = VSET(7, 6, 1, 0, 7, 6, 1, 0);
    const cc512i m7 = VSET(7, 6, 3, 2, 7, 6, 1, 0);
    const cc512i m8 = VSET(5, 4, 3, 2, 3, 2, 1, 0);

    perm_1032(t0.data(), X1Y1Z1.data());
    add_fp2_4x2x1w(t0, X1Y1Z1, t0);
    blend_0x03(t0.data(), X1Y1Z1.data(), t0.data());
    sqr_fp2_4x2x1w(ABZZ, t0);
    perm_var(t0.data(), ABZZ.data(), m0);
    perm_var(t1.data(), X1Y1Z1.data(), m1);
    blend_0xC3(t1.data(), t0.data(), t1.data());
    add_fp2_4x2x1w(t0, t0, t1);
    perm_1032(t1.data(), t0.data());
    blend_0xC3(t1.data(), t0.data(), t1.data());
    perm_var(t2.data(), t0.data(), m2);
    blend_0x33(t2.data(), t2.data(), ABZZ.data());
    add_fp2_4x2x1w(E, t1, t2);
    perm_var(t2.data(), E.data(), m3);
    blend_0x0F(t2.data(), t0.data(), t2.data());
    perm_var(t0.data(), ABZZ.data(), m4);
    blend_0x30(t2.data(), t2.data(), t0.data());
    sqr_fp2_4x2x1w(CF, t2);
    blend_0x03(t2.data(), CF.data(), ABZZ.data());
    blend_0x0F(t0.data(), CF.data(), t0.data());
    aass_fp2_4x2x1w(t0, t2, t0);
    blend_0x03(t2.data(), ABZZ.data(), CF.data());
    perm_1032(t2.data(), t2.data());
    perm_1010(t3.data(), t0.data());
    blend_0x0F(t2.data(), t3.data(), t2.data());
    sass_fp2_4x2x1w(Z3, t0, t2);
    blend_0xC0(t1.data(), Z3.data(), t1.data());
    blend_0x0F(t1.data(), t1.data(), E.data());
    sass_fp2_4x2x1w(l0, Z3, t1);
    blend_0x0F(t0.data(), l0.data(), CF.data());
    perm_var(t1.data(), t0.data(), m5);
    perm_var(t0.data(), t0.data(), m6);
    ssaa_fp2_4x2x1w(t0, t1, t0);
    blend_0x0F(t1.data(), l0.data(), t0.data());
    perm_var(t1.data(), t1.data(), m7);
    sass_fp2_4x2x1w(X3, t0, t1);
    blend_0x0F(t0.data(), X3.data(), Z3.data());
    perm_1032(t0.data(), t0.data());
    blend_0x33(t0.data(), t0.data(), E.data());
    blend_0x0F(t1.data(), E.data(), ABZZ.data());
    perm_var(t1.data(), t1.data(), m8);
    mul_fp2_4x2x1w(l12, t0, t1);
    perm_1010(t0.data(), l0.data());
    sub_fp2_4x2x1w(Y3, l12, t0);
}

void line_add_vec_v1(
    fp2_2x2x2w& l0Y3,
    fp2_2x2x2w& l1,
    fp2_2x2x2w& X3,
    fp2_2x2x2w& Z3,
    const fp2_2x2x2w& X1Y1,
    const fp2_2x2x2w& Z1Y2,
    const fp2_2x2x2w& X2)
{
    fp2_2x2x2w H, J, ZZ, US, V;
    fp2_2x2x2w t0, t1, t2, t3, t4;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(7, 6, 5, 4, 7, 6, 5, 4);

    perm_var_hl(t0.data(), Z1Y2.data(), m1);
    mul_fp2_2x2x2w(ZZ, t0, Z1Y2);
    blend_0x0F_hl(t0.data(), X2.data(), ZZ.data());
    perm_var_hl(t1.data(), ZZ.data(), m1);
    mul_fp2_2x2x2w(US, t0, t1);
    sub_fp2_2x2x2w(H, US, X1Y1);
    blend_0x0F_hl(t0.data(), Z1Y2.data(), H.data());
    add_fp2_2x2x2w(l1, H, t0);
    perm_var_hl(t0.data(), H.data(), m0);
    blend_0x0F_hl(t0.data(), l1.data(), t0.data());
    sqr_fp2_2x2x2w(t0, t0);
    blend_0x0F_hl(t1.data(), ZZ.data(), t0.data());
    sa_fp2_2x2x2w(t0, t0, t1);
    perm_var_hl(t1.data(), t1.data(), m0);
    blend_0x0F_hl(t1.data(), t1.data(), t0.data());
    sa_fp2_2x2x2w(Z3, t0, t1);
    blend_0x0F_hl(t0.data(), H.data(), l1.data());
    perm_var_hl(t3.data(), Z3.data(), m0);
    blend_0x0F_hl(t1.data(), t3.data(), X2.data());
    mul_fp2_2x2x2w(J, t0, t1);
    perm_var_hl(t1.data(), l1.data(), m0);
    blend_0x0F_hl(t2.data(), t1.data(), Z1Y2.data());
    blend_0x0F_hl(t0.data(), t1.data(), t3.data());
    mul_fp2_2x2x2w(t0, t0, t2);
    perm_var_hl(t2.data(), J.data(), m0);
    blend_0x0F_hl(t2.data(), t3.data(), t2.data());
    mul_fp2_2x2x2w(V, X1Y1, t2);
    perm_var_hl(t3.data(), V.data(), m0);
    blend_0x0F_hl(t2.data(), t0.data(), t3.data());
    blend_0x0F_hl(t3.data(), J.data(), t3.data());
    sa_fp2_2x2x2w(t2, t2, t3);
    blend_0x0F_hl(t3.data(), t2.data(), V.data());
    perm_var_hl(t2.data(), t2.data(), m0);
    blend_0x0F_hl(t2.data(), t2.data(), t3.data());
    sa_fp2_2x2x2w(X3, t3, t2);
    blend_0x0F_hl(t3.data(), V.data(), J.data());
    blend_0x0F_hl(t0.data(), X3.data(), t0.data());
    sub_fp2_2x2x2w(t0, t3, t0);
    mul_fp2_2x2x2w(t2, t0, t1);
    blend_0x0F_hl(t0.data(), t2.data(), t0.data());
    perm_var_hl(t1.data(), X3.data(), m0);
    blend_0x0F_hl(t1.data(), t1.data(), t0.data());
    sa_fp2_2x2x2w(l0Y3, t0, t1);
}

void line_by_Px2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& Px2)
{
    fp2x2_4x2x1w tt0;
    FpPack8x1::mul(tt0, a, Px2);
    UFpPack8x1::redc(r, tt0);
}

void line_by_Px2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& Px2)
{
    fp2x2_2x2x2w tt0;
    FpPack4x2::mul(tt0, a, Px2);
    UFpPack4x2::redc(r, tt0);
}

void line_add_vector(Fp6& line, G2Point& T, const G2Point& R, const G2AffinePoint& Q)
{
    BLST_ASSERT_NO_OVERLAP_PTR(&line, sizeof(line), &T, sizeof(T));
    BLST_ASSERT_NO_OVERLAP_PTR(&line, sizeof(line), &R, sizeof(R));
    BLST_ASSERT_NO_OVERLAP_PTR(&line, sizeof(line), &Q, sizeof(Q));

    fp2_2x2x2w X1Y1, Z1Y2, X2, l0Y3, l1, X3, Z3;
    std::array<cc512i, SWORDS / 2> in0, in1, in2;
    std::array<cc512i, SWORDS / 2> out0, out1, out2, out3;

    for (u64 i = 0; i < SWORDS / 2; ++i) {
        in0[i] = VSET(
            static_cast<i64>(R.X[1][i + (SWORDS / 2)]), static_cast<i64>(R.X[1][i]),
            static_cast<i64>(R.X[0][i + (SWORDS / 2)]), static_cast<i64>(R.X[0][i]),
            static_cast<i64>(R.Y[1][i + (SWORDS / 2)]), static_cast<i64>(R.Y[1][i]),
            static_cast<i64>(R.Y[0][i + (SWORDS / 2)]), static_cast<i64>(R.Y[0][i]));
        in1[i] = VSET(
            static_cast<i64>(R.Z[1][i + (SWORDS / 2)]), static_cast<i64>(R.Z[1][i]),
            static_cast<i64>(R.Z[0][i + (SWORDS / 2)]), static_cast<i64>(R.Z[0][i]),
            static_cast<i64>(Q.Y[1][i + (SWORDS / 2)]), static_cast<i64>(Q.Y[1][i]),
            static_cast<i64>(Q.Y[0][i + (SWORDS / 2)]), static_cast<i64>(Q.Y[0][i]));
        in2[i] = VSET(
            static_cast<i64>(Q.X[1][i + (SWORDS / 2)]), static_cast<i64>(Q.X[1][i]),
            static_cast<i64>(Q.X[0][i + (SWORDS / 2)]), static_cast<i64>(Q.X[0][i]),
            static_cast<i64>(Q.X[1][i + (SWORDS / 2)]), static_cast<i64>(Q.X[1][i]),
            static_cast<i64>(Q.X[0][i + (SWORDS / 2)]), static_cast<i64>(Q.X[0][i]));
    }

    FpPack4x2::conv_64to48(X1Y1, in0.data());
    FpPack4x2::conv_64to48(Z1Y2, in1.data());
    FpPack4x2::conv_64to48(X2, in2.data());

    line_add_vec_v1(l0Y3, l1, X3, Z3, X1Y1, Z1Y2, X2);

    FpPack4x2::carryp(l0Y3);
    FpPack4x2::carryp(l1);
    FpPack4x2::carryp(X3);
    FpPack4x2::carryp(Z3);

    FpPack4x2::conv_48to64(out0.data(), l0Y3);
    FpPack4x2::conv_48to64(out1.data(), l1);
    FpPack4x2::conv_48to64(out2.data(), X3);
    FpPack4x2::conv_48to64(out3.data(), Z3);

    for (u64 i = 0; i < SWORDS / 2; ++i) {
        std::array<u64, 8> lanes0, lanes1, lanes2, lanes3;
        store_u64_lanes(lanes0, out0[i]);
        store_u64_lanes(lanes1, out1[i]);
        store_u64_lanes(lanes2, out2[i]);
        store_u64_lanes(lanes3, out3[i]);

        line[0][0][i] = lanes0[0];
        line[0][0][i + (SWORDS / 2)] = lanes0[1];
        line[0][1][i] = lanes0[2];
        line[0][1][i + (SWORDS / 2)] = lanes0[3];

        line[1][0][i] = lanes1[0];
        line[1][0][i + (SWORDS / 2)] = lanes1[1];
        line[1][1][i] = lanes1[2];
        line[1][1][i + (SWORDS / 2)] = lanes1[3];

        T.X[0][i] = lanes2[4];
        T.X[0][i + (SWORDS / 2)] = lanes2[5];
        T.X[1][i] = lanes2[6];
        T.X[1][i + (SWORDS / 2)] = lanes2[7];

        T.Y[0][i] = lanes0[4];
        T.Y[0][i + (SWORDS / 2)] = lanes0[5];
        T.Y[1][i] = lanes0[6];
        T.Y[1][i + (SWORDS / 2)] = lanes0[7];

        T.Z[0][i] = lanes3[4];
        T.Z[0][i + (SWORDS / 2)] = lanes3[5];
        T.Z[1][i] = lanes3[6];
        T.Z[1][i + (SWORDS / 2)] = lanes3[7];

        line[2][0][i] = T.Z[0][i];
        line[2][0][i + (SWORDS / 2)] = T.Z[0][i + (SWORDS / 2)];
        line[2][1][i] = T.Z[1][i];
        line[2][1][i + (SWORDS / 2)] = T.Z[1][i + (SWORDS / 2)];
    }
}

void line_dbl_vector(Fp6& line, G2Point& T, const G2Point& Q)
{
    BLST_ASSERT_NO_OVERLAP_PTR(&line, sizeof(line), &T, sizeof(T));
    BLST_ASSERT_NO_OVERLAP_PTR(&line, sizeof(line), &Q, sizeof(Q));

    fp2_4x2x1w X1Y1Z1, l0, l12, X3, Y3, Z3;
    std::array<cc512i, SWORDS> in0;
    std::array<cc512i, SWORDS> out0, out1, out2, out3, out4;

    for (u64 i = 0; i < SWORDS; ++i) {
        in0[i] = VSET(
            static_cast<i64>(Q.Y[1][i]), static_cast<i64>(Q.Y[0][i]),
            static_cast<i64>(Q.X[1][i]), static_cast<i64>(Q.X[0][i]),
            static_cast<i64>(Q.Z[1][i]), static_cast<i64>(Q.Z[0][i]),
            static_cast<i64>(Q.Y[1][i]), static_cast<i64>(Q.Y[0][i]));
    }

    FpPack8x1::conv_64to48(X1Y1Z1, in0.data());

    line_dbl_vec_v2(l0, l12, X3, Y3, Z3, X1Y1Z1);

    FpPack8x1::conv_48to64(out0.data(), l0);
    FpPack8x1::conv_48to64(out1.data(), l12);
    FpPack8x1::conv_48to64(out2.data(), X3);
    FpPack8x1::conv_48to64(out3.data(), Y3);
    FpPack8x1::conv_48to64(out4.data(), Z3);

    for (u64 i = 0; i < SWORDS; ++i) {
        std::array<u64, 8> lanes0, lanes1, lanes2, lanes3, lanes4;
        store_u64_lanes(lanes0, out0[i]);
        store_u64_lanes(lanes1, out1[i]);
        store_u64_lanes(lanes2, out2[i]);
        store_u64_lanes(lanes3, out3[i]);
        store_u64_lanes(lanes4, out4[i]);

        line[0][0][i] = lanes0[2];
        line[0][1][i] = lanes0[3];
        line[1][0][i] = lanes1[4];
        line[1][1][i] = lanes1[5];
        line[2][0][i] = lanes1[2];
        line[2][1][i] = lanes1[3];

        T.X[0][i] = lanes2[6];
        T.X[1][i] = lanes2[7];
        T.Y[0][i] = lanes3[6];
        T.Y[1][i] = lanes3[7];
        T.Z[0][i] = lanes4[0];
        T.Z[1][i] = lanes4[1];
    }
}

} // namespace blst

#endif
