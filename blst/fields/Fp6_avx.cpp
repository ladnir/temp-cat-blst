#include "blst/fields/Fp6_avx.h"
#include "blst/backend/avx512/avx_layout.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
namespace
{
    inline void shuf_01_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        shuf_01(r[0].data(), a[0].data());
        shuf_01(r[1].data(), a[1].data());
    }

    inline void blend_0x11_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
    {
        blend_0x11(r[0].data(), a[0].data(), b[0].data());
        blend_0x11(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x33_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
    {
        blend_0x33(r[0].data(), a[0].data(), b[0].data());
        blend_0x33(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x55_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
    {
        blend_0x55(r[0].data(), a[0].data(), b[0].data());
        blend_0x55(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x77_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
    {
        blend_0x77(r[0].data(), a[0].data(), b[0].data());
        blend_0x77(r[1].data(), a[1].data(), b[1].data());
    }

    inline void perm_0123_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        perm_0123(r[0].data(), a[0].data());
        perm_0123(r[1].data(), a[1].data());
    }

    inline void perm_0300_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        perm_0300(r[0].data(), a[0].data());
        perm_0300(r[1].data(), a[1].data());
    }

    inline void perm_1032_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        perm_1032(r[0].data(), a[0].data());
        perm_1032(r[1].data(), a[1].data());
    }

    inline void perm_1200_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        perm_1200(r[0].data(), a[0].data());
        perm_1200(r[1].data(), a[1].data());
    }

    inline void perm_2000_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        perm_2000(r[0].data(), a[0].data());
        perm_2000(r[1].data(), a[1].data());
    }

    inline void shuf_00_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        shuf_00_dl(r[0].data(), a[0].data());
        shuf_00_dl(r[1].data(), a[1].data());
    }

    inline void shuf_01_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        shuf_01_dl(r[0].data(), a[0].data());
        shuf_01_dl(r[1].data(), a[1].data());
    }

    inline void shuf_11_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        shuf_11_dl(r[0].data(), a[0].data());
        shuf_11_dl(r[1].data(), a[1].data());
    }

    inline void blend_0x33_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
    {
        blend_0x33_dl(r[0].data(), a[0].data(), b[0].data());
        blend_0x33_dl(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x55_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
    {
        blend_0x55_dl(r[0].data(), a[0].data(), b[0].data());
        blend_0x55_dl(r[1].data(), a[1].data(), b[1].data());
    }

    inline void perm_2323_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        perm_2323_dl(r[0].data(), a[0].data());
        perm_2323_dl(r[1].data(), a[1].data());
    }

    inline void perm_3000_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        perm_3000_dl(r[0].data(), a[0].data());
        perm_3000_dl(r[1].data(), a[1].data());
    }
} // namespace

void mul_fp6x2_4x2x1x1w(
    UFp2Pack8x1x1& r01,
    UFp2Pack8x1x1& r2,
    const Fp2Pack8x1x1& ab0,
    const Fp2Pack8x1x1& ab1,
    const Fp2Pack8x1x1& ab2)
{
    Fp2Pack8x1x1 t0, t1, t2, t3, t4;
    UFp2Pack8x1x1 tt0, tt1, tt2, tt3, tt4, tt5;

    Fp2Pack8x1x1::add(t0, ab1, ab2);
    Fp2Pack8x1x1::add(t1, ab0, ab1);
    Fp2Pack8x1x1::add(t2, ab0, ab2);
    shuf_01_Fp2Pack8x1x1(t3, ab1);
    blend_0x55_Fp2Pack8x1x1(t3, t3, ab0);
    shuf_01_Fp2Pack8x1x1(t4, ab0);
    blend_0x55_Fp2Pack8x1x1(t4, ab1, t4);
    UFp2Pack8x1x1::mul(tt0, t3, t4);
    shuf_01_Fp2Pack8x1x1(t3, ab2);
    blend_0x55_Fp2Pack8x1x1(t3, t3, t0);
    shuf_01_Fp2Pack8x1x1(t4, t0);
    blend_0x55_Fp2Pack8x1x1(t4, ab2, t4);
    UFp2Pack8x1x1::mul(tt1, t3, t4);
    shuf_01_Fp2Pack8x1x1(t3, t1);
    blend_0x55_Fp2Pack8x1x1(t3, t3, t2);
    shuf_01_Fp2Pack8x1x1(t4, t2);
    blend_0x55_Fp2Pack8x1x1(t4, t1, t4);
    UFp2Pack8x1x1::mul(tt2, t3, t4);
    blend_0x55_UFp2Pack8x1x1(tt3, tt2, tt1);
    shuf_11_UFp2Pack8x1x1(tt4, tt0);
    UFp2Pack8x1x1::sub(tt3, tt3, tt4);
    shuf_01_UFp2Pack8x1x1(tt2, tt2);
    blend_0x55_UFp2Pack8x1x1(tt2, tt2, tt3);
    shuf_11_UFp2Pack8x1x1(tt1, tt1);
    UFp2Pack8x1x1::sub(tt2, tt2, tt1);
    shuf_11_UFp2Pack8x1x1(tt5, tt2);
    blend_0x55_UFp2Pack8x1x1(tt5, tt3, tt5);
    shuf_00_UFp2Pack8x1x1(tt0, tt0);
    UFp2Pack8x1x1::sub(tt5, tt5, tt0);
    blend_0x55_UFp2Pack8x1x1(tt2, tt1, tt2);
    UFp2Pack8x1x1::mul_by_u_plus_1(tt1, tt2);
    blend_0x55_UFp2Pack8x1x1(tt0, tt5, tt0);
    UFp2Pack8x1x1::add(r01, tt0, tt1);
    UFp2Pack8x1x1::add(r2, tt5, tt4);
}

void mul_fp6x2_2x2x2x1w(
    UFp2Pack4x2x1& r01,
    UFp2Pack4x2x1& r2,
    const Fp2Pack4x2x1& ab0,
    const Fp2Pack4x2x1& ab1,
    const Fp2Pack4x2x1& ab2)
{
    Fp2Pack4x2x1 t0, t1, t2, t3, t4;
    UFp2Pack4x2x1 tt0, tt1, tt2, tt3, tt4, tt5;

    Fp2Pack4x2x1::add(t0, ab1, ab2);
    Fp2Pack4x2x1::add(t1, ab0, ab1);
    Fp2Pack4x2x1::add(t2, ab0, ab2);
    perm_1032(t3.data(), ab1.data());
    blend_0x33(t3.data(), t3.data(), ab0.data());
    perm_1032(t4.data(), ab0.data());
    blend_0x33(t4.data(), ab1.data(), t4.data());
    UFp2Pack4x2x1::mul(tt0, t3, t4);
    perm_1032(t3.data(), ab2.data());
    blend_0x33(t3.data(), t3.data(), t0.data());
    perm_1032(t4.data(), t0.data());
    blend_0x33(t4.data(), ab2.data(), t4.data());
    UFp2Pack4x2x1::mul(tt1, t3, t4);
    perm_1032(t3.data(), t1.data());
    blend_0x33(t3.data(), t3.data(), t2.data());
    perm_1032(t4.data(), t2.data());
    blend_0x33(t4.data(), t1.data(), t4.data());
    UFp2Pack4x2x1::mul(tt2, t3, t4);
    blend_0x33_dl(tt3.data(), tt2.data(), tt1.data());
    perm_3232_dl(tt4.data(), tt0.data());
    UFp2Pack4x2x1::sub(tt3, tt3, tt4);
    perm_1032_dl(tt2.data(), tt2.data());
    blend_0x33_dl(tt2.data(), tt2.data(), tt3.data());
    perm_3232_dl(tt1.data(), tt1.data());
    UFp2Pack4x2x1::sub(tt2, tt2, tt1);
    perm_3232_dl(tt5.data(), tt2.data());
    blend_0x33_dl(tt5.data(), tt3.data(), tt5.data());
    perm_1010_dl(tt0.data(), tt0.data());
    UFp2Pack4x2x1::sub(tt5, tt5, tt0);
    blend_0x33_dl(tt2.data(), tt1.data(), tt2.data());
    UFp2Pack4x2x1::mul_by_u_plus_1(tt1, tt2);
    blend_0x33_dl(tt0.data(), tt5.data(), tt0.data());
    UFp2Pack4x2x1::add(r01, tt0, tt1);
    UFp2Pack4x2x1::add(r2, tt5, tt4);
}

void mul_fp6x2_1x2x2x2w(
    UFp2Pack2x2x2& r01,
    UFp2Pack2x2x2& r2,
    const Fp2Pack2x2x2& ab0,
    const Fp2Pack2x2x2& ab1,
    const Fp2Pack2x2x2& ab2)
{
    Fp2Pack2x2x2 t0, t1, t2, t3, t4;
    UFp2Pack2x2x2 tt0, tt1, tt2, tt3, tt4, tt5;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(7, 6, 5, 4, 7, 6, 5, 4);
    const cc512i m2 = VSET(3, 2, 1, 0, 3, 2, 1, 0);

    Fp2Pack2x2x2::add(t0, ab1, ab2);
    Fp2Pack2x2x2::add(t1, ab0, ab1);
    Fp2Pack2x2x2::add(t2, ab0, ab2);
    perm_var_hl(t3.data(), ab1.data(), m0);
    blend_0x0F_hl(t3.data(), t3.data(), ab0.data());
    perm_var_hl(t4.data(), ab0.data(), m0);
    blend_0x0F_hl(t4.data(), ab1.data(), t4.data());
    UFp2Pack2x2x2::mul(tt0, t3, t4);
    perm_var_hl(t3.data(), ab2.data(), m0);
    blend_0x0F_hl(t3.data(), t3.data(), t0.data());
    perm_var_hl(t4.data(), t0.data(), m0);
    blend_0x0F_hl(t4.data(), ab2.data(), t4.data());
    UFp2Pack2x2x2::mul(tt1, t3, t4);
    perm_var_hl(t3.data(), t1.data(), m0);
    blend_0x0F_hl(t3.data(), t3.data(), t2.data());
    perm_var_hl(t4.data(), t2.data(), m0);
    blend_0x0F_hl(t4.data(), t1.data(), t4.data());
    UFp2Pack2x2x2::mul(tt2, t3, t4);
    blend_0x0F_vl(tt3.data(), tt2.data(), tt1.data());
    perm_var_vl(tt4.data(), tt0.data(), m1);
    UFp2Pack2x2x2::sub(tt3, tt3, tt4);
    perm_var_vl(tt2.data(), tt2.data(), m0);
    blend_0x0F_vl(tt2.data(), tt2.data(), tt3.data());
    perm_var_vl(tt1.data(), tt1.data(), m1);
    UFp2Pack2x2x2::sub(tt2, tt2, tt1);
    perm_var_vl(tt5.data(), tt2.data(), m1);
    blend_0x0F_vl(tt5.data(), tt3.data(), tt5.data());
    perm_var_vl(tt0.data(), tt0.data(), m2);
    UFp2Pack2x2x2::sub(tt5, tt5, tt0);
    blend_0x0F_vl(tt2.data(), tt1.data(), tt2.data());
    UFp2Pack2x2x2::mul_by_u_plus_1(tt1, tt2);
    blend_0x0F_vl(tt0.data(), tt5.data(), tt0.data());
    UFp2Pack2x2x2::add(r01, tt0, tt1);
    UFp2Pack2x2x2::add(r2, tt5, tt4);
}

void mul_by_xy00z0_fp6x2_2x2x2x1w(
    UFp2Pack4x2x1& r01,
    UFp2Pack4x2x1& r2,
    UFp2Pack4x2x1& r3,
    UFp2Pack4x2x1& r45,
    const Fp2Pack4x2x1& a01,
    const Fp2Pack4x2x1& a2,
    const Fp2Pack4x2x1& b01,
    const Fp2Pack4x2x1& b4)
{
    Fp2Pack4x2x1 t0, t1;
    UFp2Pack4x2x1 tt0, tt1, tt2;

    perm_1032(t0.data(), b01.data());
    blend_0x33(t1.data(), b4.data(), t0.data());
    UFp2Pack4x2x1::mul(tt0, t1, a2);
    UFp2Pack4x2x1::mul_by_u_plus_1(r3, tt0);
    UFp2Pack4x2x1::mul(r45, a01, b4);
    UFp2Pack4x2x1::mul(tt0, a01, b01);
    blend_0x33(t0.data(), t0.data(), a01.data());
    perm_1032(t1.data(), a01.data());
    blend_0x33(t1.data(), b01.data(), t1.data());
    Fp2Pack4x2x1::add(t0, t0, t1);
    blend_0x33(t1.data(), t0.data(), b01.data());
    perm_1032(t0.data(), t0.data());
    blend_0x33(t0.data(), t0.data(), a2.data());
    UFp2Pack4x2x1::mul(tt1, t1, t0);
    perm_1032_dl(tt2.data(), tt0.data());
    UFp2Pack4x2x1::sax2(r2, tt1, tt2);
    blend_0x33_dl(tt1.data(), r2.data(), r3.data());
    UFp2Pack4x2x1::sax2(r01, tt1, tt0);
}

void mul_by_xy00z0_fp6x2_2x4x1x1w(
    UFp2Pack8x1x1& r01,
    UFp2Pack8x1x1& r23,
    UFp2Pack8x1x1& r45,
    const Fp2Pack8x1x1& a,
    const Fp2Pack8x1x1& b)
{
    Fp2Pack8x1x1 t0, t1;
    UFp2Pack8x1x1 tt0, tt1, tt2;

    perm_2000_Fp2Pack8x1x1(t0, b);
    UFp2Pack8x1x1::mul(r45, a, t0);
    blend_0x33_Fp2Pack8x1x1(t0, b, a);
    perm_0300_Fp2Pack8x1x1(t1, t0);
    perm_1200_Fp2Pack8x1x1(t0, t0);
    Fp2Pack8x1x1::add(t0, t0, t1);
    blend_0x77_Fp2Pack8x1x1(t1, t0, a);
    perm_0123_Fp2Pack8x1x1(t1, t1);
    perm_1032_Fp2Pack8x1x1(t0, t0);
    blend_0x11_Fp2Pack8x1x1(t0, b, t0);
    UFp2Pack8x1x1::mul(tt2, t0, t1);
    blend_0x33_UFp2Pack8x1x1(tt1, r45, tt2);
    perm_2323_UFp2Pack8x1x1(tt0, tt2);
    blend_0x33_UFp2Pack8x1x1(tt0, tt1, tt0);
    blend_0x33_dl(tt1[1].data(), tt1[0].data(), tt1[1].data());
    blend_0x33_dl(tt0[0].data(), tt0[1].data(), tt0[0].data());
    UFpPack8x1::ssasx2(r23[0], tt1[0], tt0[0]);
    UFpPack8x1::aaasx2(r23[1], tt1[1], tt0[1]);
    perm_3000_UFp2Pack8x1x1(tt0, r23);
    UFp2Pack8x1x1::asx4(r01, tt0, tt2);
}

} // namespace blst

#endif
