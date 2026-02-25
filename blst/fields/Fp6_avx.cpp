#include "blst/fields/Fp6_avx.h"
#include "blst/backend/avx512/avx_layout.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
namespace
{
    inline void shuf_01_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a)
    {
        shuf_01(r[0], a[0]);
        shuf_01(r[1], a[1]);
    }

    inline void blend_0x11_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a, const fp2_8x1x1w& b)
    {
        blend_0x11(r[0], a[0], b[0]);
        blend_0x11(r[1], a[1], b[1]);
    }

    inline void blend_0x33_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a, const fp2_8x1x1w& b)
    {
        blend_0x33(r[0], a[0], b[0]);
        blend_0x33(r[1], a[1], b[1]);
    }

    inline void blend_0x55_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a, const fp2_8x1x1w& b)
    {
        blend_0x55(r[0], a[0], b[0]);
        blend_0x55(r[1], a[1], b[1]);
    }

    inline void blend_0x77_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a, const fp2_8x1x1w& b)
    {
        blend_0x77(r[0], a[0], b[0]);
        blend_0x77(r[1], a[1], b[1]);
    }

    inline void perm_0123_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a)
    {
        perm_0123(r[0], a[0]);
        perm_0123(r[1], a[1]);
    }

    inline void perm_0300_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a)
    {
        perm_0300(r[0], a[0]);
        perm_0300(r[1], a[1]);
    }

    inline void perm_1032_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a)
    {
        perm_1032(r[0], a[0]);
        perm_1032(r[1], a[1]);
    }

    inline void perm_1200_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a)
    {
        perm_1200(r[0], a[0]);
        perm_1200(r[1], a[1]);
    }

    inline void perm_2000_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a)
    {
        perm_2000(r[0], a[0]);
        perm_2000(r[1], a[1]);
    }

    inline void shuf_00_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a)
    {
        shuf_00_dl(r[0], a[0]);
        shuf_00_dl(r[1], a[1]);
    }

    inline void shuf_01_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a)
    {
        shuf_01_dl(r[0], a[0]);
        shuf_01_dl(r[1], a[1]);
    }

    inline void shuf_11_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a)
    {
        shuf_11_dl(r[0], a[0]);
        shuf_11_dl(r[1], a[1]);
    }

    inline void blend_0x33_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const fp2x2_8x1x1w& b)
    {
        blend_0x33_dl(r[0], a[0], b[0]);
        blend_0x33_dl(r[1], a[1], b[1]);
    }

    inline void blend_0x55_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const fp2x2_8x1x1w& b)
    {
        blend_0x55_dl(r[0], a[0], b[0]);
        blend_0x55_dl(r[1], a[1], b[1]);
    }

    inline void perm_2323_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a)
    {
        perm_2323_dl(r[0], a[0]);
        perm_2323_dl(r[1], a[1]);
    }

    inline void perm_3000_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a)
    {
        perm_3000_dl(r[0], a[0]);
        perm_3000_dl(r[1], a[1]);
    }
} // namespace

void mul_fp6x2_4x2x1x1w(
    fp2x2_8x1x1w& r01,
    fp2x2_8x1x1w& r2,
    const fp2_8x1x1w& ab0,
    const fp2_8x1x1w& ab1,
    const fp2_8x1x1w& ab2)
{
    fp2_8x1x1w t0, t1, t2, t3, t4;
    fp2x2_8x1x1w tt0, tt1, tt2, tt3, tt4, tt5;

    add_fp2_8x1x1w(t0, ab1, ab2);
    add_fp2_8x1x1w(t1, ab0, ab1);
    add_fp2_8x1x1w(t2, ab0, ab2);
    shuf_01_fp2_8x1x1w(t3, ab1);
    blend_0x55_fp2_8x1x1w(t3, t3, ab0);
    shuf_01_fp2_8x1x1w(t4, ab0);
    blend_0x55_fp2_8x1x1w(t4, ab1, t4);
    mul_fp2x2_8x1x1w(tt0, t3, t4);
    shuf_01_fp2_8x1x1w(t3, ab2);
    blend_0x55_fp2_8x1x1w(t3, t3, t0);
    shuf_01_fp2_8x1x1w(t4, t0);
    blend_0x55_fp2_8x1x1w(t4, ab2, t4);
    mul_fp2x2_8x1x1w(tt1, t3, t4);
    shuf_01_fp2_8x1x1w(t3, t1);
    blend_0x55_fp2_8x1x1w(t3, t3, t2);
    shuf_01_fp2_8x1x1w(t4, t2);
    blend_0x55_fp2_8x1x1w(t4, t1, t4);
    mul_fp2x2_8x1x1w(tt2, t3, t4);
    blend_0x55_fp2x2_8x1x1w(tt3, tt2, tt1);
    shuf_11_fp2x2_8x1x1w(tt4, tt0);
    sub_fp2x2_8x1x1w(tt3, tt3, tt4);
    shuf_01_fp2x2_8x1x1w(tt2, tt2);
    blend_0x55_fp2x2_8x1x1w(tt2, tt2, tt3);
    shuf_11_fp2x2_8x1x1w(tt1, tt1);
    sub_fp2x2_8x1x1w(tt2, tt2, tt1);
    shuf_11_fp2x2_8x1x1w(tt5, tt2);
    blend_0x55_fp2x2_8x1x1w(tt5, tt3, tt5);
    shuf_00_fp2x2_8x1x1w(tt0, tt0);
    sub_fp2x2_8x1x1w(tt5, tt5, tt0);
    blend_0x55_fp2x2_8x1x1w(tt2, tt1, tt2);
    mul_by_u_plus_1_fp2x2_8x1x1w(tt1, tt2);
    blend_0x55_fp2x2_8x1x1w(tt0, tt5, tt0);
    add_fp2x2_8x1x1w(r01, tt0, tt1);
    add_fp2x2_8x1x1w(r2, tt5, tt4);
}

void mul_fp6x2_2x2x2x1w(
    fp2x2_4x2x1w& r01,
    fp2x2_4x2x1w& r2,
    const fp2_4x2x1w& ab0,
    const fp2_4x2x1w& ab1,
    const fp2_4x2x1w& ab2)
{
    fp2_4x2x1w t0, t1, t2, t3, t4;
    fp2x2_4x2x1w tt0, tt1, tt2, tt3, tt4, tt5;

    add_fp2_4x2x1w(t0, ab1, ab2);
    add_fp2_4x2x1w(t1, ab0, ab1);
    add_fp2_4x2x1w(t2, ab0, ab2);
    perm_1032(t3, ab1);
    blend_0x33(t3, t3, ab0);
    perm_1032(t4, ab0);
    blend_0x33(t4, ab1, t4);
    mul_fp2x2_4x2x1w(tt0, t3, t4);
    perm_1032(t3, ab2);
    blend_0x33(t3, t3, t0);
    perm_1032(t4, t0);
    blend_0x33(t4, ab2, t4);
    mul_fp2x2_4x2x1w(tt1, t3, t4);
    perm_1032(t3, t1);
    blend_0x33(t3, t3, t2);
    perm_1032(t4, t2);
    blend_0x33(t4, t1, t4);
    mul_fp2x2_4x2x1w(tt2, t3, t4);
    blend_0x33_dl(tt3, tt2, tt1);
    perm_3232_dl(tt4, tt0);
    sub_fp2x2_4x2x1w(tt3, tt3, tt4);
    perm_1032_dl(tt2, tt2);
    blend_0x33_dl(tt2, tt2, tt3);
    perm_3232_dl(tt1, tt1);
    sub_fp2x2_4x2x1w(tt2, tt2, tt1);
    perm_3232_dl(tt5, tt2);
    blend_0x33_dl(tt5, tt3, tt5);
    perm_1010_dl(tt0, tt0);
    sub_fp2x2_4x2x1w(tt5, tt5, tt0);
    blend_0x33_dl(tt2, tt1, tt2);
    mul_by_u_plus_1_fp2x2_4x2x1w(tt1, tt2);
    blend_0x33_dl(tt0, tt5, tt0);
    add_fp2x2_4x2x1w(r01, tt0, tt1);
    add_fp2x2_4x2x1w(r2, tt5, tt4);
}

void mul_fp6x2_1x2x2x2w(
    fp2x2_2x2x2w& r01,
    fp2x2_2x2x2w& r2,
    const fp2_2x2x2w& ab0,
    const fp2_2x2x2w& ab1,
    const fp2_2x2x2w& ab2)
{
    fp2_2x2x2w t0, t1, t2, t3, t4;
    fp2x2_2x2x2w tt0, tt1, tt2, tt3, tt4, tt5;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(7, 6, 5, 4, 7, 6, 5, 4);
    const cc512i m2 = VSET(3, 2, 1, 0, 3, 2, 1, 0);

    add_fp2_2x2x2w(t0, ab1, ab2);
    add_fp2_2x2x2w(t1, ab0, ab1);
    add_fp2_2x2x2w(t2, ab0, ab2);
    perm_var_hl(t3, ab1, m0);
    blend_0x0F_hl(t3, t3, ab0);
    perm_var_hl(t4, ab0, m0);
    blend_0x0F_hl(t4, ab1, t4);
    mul_fp2x2_2x2x2w(tt0, t3, t4);
    perm_var_hl(t3, ab2, m0);
    blend_0x0F_hl(t3, t3, t0);
    perm_var_hl(t4, t0, m0);
    blend_0x0F_hl(t4, ab2, t4);
    mul_fp2x2_2x2x2w(tt1, t3, t4);
    perm_var_hl(t3, t1, m0);
    blend_0x0F_hl(t3, t3, t2);
    perm_var_hl(t4, t2, m0);
    blend_0x0F_hl(t4, t1, t4);
    mul_fp2x2_2x2x2w(tt2, t3, t4);
    blend_0x0F_vl(tt3, tt2, tt1);
    perm_var_vl(tt4, tt0, m1);
    sub_fp2x2_2x2x2w(tt3, tt3, tt4);
    perm_var_vl(tt2, tt2, m0);
    blend_0x0F_vl(tt2, tt2, tt3);
    perm_var_vl(tt1, tt1, m1);
    sub_fp2x2_2x2x2w(tt2, tt2, tt1);
    perm_var_vl(tt5, tt2, m1);
    blend_0x0F_vl(tt5, tt3, tt5);
    perm_var_vl(tt0, tt0, m2);
    sub_fp2x2_2x2x2w(tt5, tt5, tt0);
    blend_0x0F_vl(tt2, tt1, tt2);
    mul_by_u_plus_1_fp2x2_2x2x2w(tt1, tt2);
    blend_0x0F_vl(tt0, tt5, tt0);
    add_fp2x2_2x2x2w(r01, tt0, tt1);
    add_fp2x2_2x2x2w(r2, tt5, tt4);
}

void mul_by_xy00z0_fp6x2_2x2x2x1w(
    fp2x2_4x2x1w& r01,
    fp2x2_4x2x1w& r2,
    fp2x2_4x2x1w& r3,
    fp2x2_4x2x1w& r45,
    const fp2_4x2x1w& a01,
    const fp2_4x2x1w& a2,
    const fp2_4x2x1w& b01,
    const fp2_4x2x1w& b4)
{
    fp2_4x2x1w t0, t1;
    fp2x2_4x2x1w tt0, tt1, tt2;

    perm_1032(t0, b01);
    blend_0x33(t1, b4, t0);
    mul_fp2x2_4x2x1w(tt0, t1, a2);
    mul_by_u_plus_1_fp2x2_4x2x1w(r3, tt0);
    mul_fp2x2_4x2x1w(r45, a01, b4);
    mul_fp2x2_4x2x1w(tt0, a01, b01);
    blend_0x33(t0, t0, a01);
    perm_1032(t1, a01);
    blend_0x33(t1, b01, t1);
    add_fp2_4x2x1w(t0, t0, t1);
    blend_0x33(t1, t0, b01);
    perm_1032(t0, t0);
    blend_0x33(t0, t0, a2);
    mul_fp2x2_4x2x1w(tt1, t1, t0);
    perm_1032_dl(tt2, tt0);
    sax2_fp2x2_4x2x1w(r2, tt1, tt2);
    blend_0x33_dl(tt1, r2, r3);
    sax2_fp2x2_4x2x1w(r01, tt1, tt0);
}

void mul_by_xy00z0_fp6x2_2x4x1x1w(
    fp2x2_8x1x1w& r01,
    fp2x2_8x1x1w& r23,
    fp2x2_8x1x1w& r45,
    const fp2_8x1x1w& a,
    const fp2_8x1x1w& b)
{
    fp2_8x1x1w t0, t1;
    fp2x2_8x1x1w tt0, tt1, tt2;

    perm_2000_fp2_8x1x1w(t0, b);
    mul_fp2x2_8x1x1w(r45, a, t0);
    blend_0x33_fp2_8x1x1w(t0, b, a);
    perm_0300_fp2_8x1x1w(t1, t0);
    perm_1200_fp2_8x1x1w(t0, t0);
    add_fp2_8x1x1w(t0, t0, t1);
    blend_0x77_fp2_8x1x1w(t1, t0, a);
    perm_0123_fp2_8x1x1w(t1, t1);
    perm_1032_fp2_8x1x1w(t0, t0);
    blend_0x11_fp2_8x1x1w(t0, b, t0);
    mul_fp2x2_8x1x1w(tt2, t0, t1);
    blend_0x33_fp2x2_8x1x1w(tt1, r45, tt2);
    perm_2323_fp2x2_8x1x1w(tt0, tt2);
    blend_0x33_fp2x2_8x1x1w(tt0, tt1, tt0);
    blend_0x33_dl(tt1[1], tt1[0], tt1[1]);
    blend_0x33_dl(tt0[0], tt0[1], tt0[0]);
    ssasx2_fpx2_8x1w(r23[0], tt1[0], tt0[0]);
    aaasx2_fpx2_8x1w(r23[1], tt1[1], tt0[1]);
    perm_3000_fp2x2_8x1x1w(tt0, r23);
    asx4_fp2x2_8x1x1w(r01, tt0, tt2);
}

} // namespace blst

#endif
