#include "blst/fields/Fp12_avx.h"
#include "blst/backend/avx512/avx_layout.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include <array>

namespace blst
{
namespace
{
    inline void shuf_01_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a)
    {
        shuf_01_dl(r[0].data(), a[0].data());
        shuf_01_dl(r[1].data(), a[1].data());
    }

    inline void blend_0x01_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const fp2x2_8x1x1w& b)
    {
        blend_0x01_dl(r[0].data(), a[0].data(), b[0].data());
        blend_0x01_dl(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x03_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const fp2x2_8x1x1w& b)
    {
        blend_0x03_dl(r[0].data(), a[0].data(), b[0].data());
        blend_0x03_dl(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x33_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const fp2x2_8x1x1w& b)
    {
        blend_0x33_dl(r[0].data(), a[0].data(), b[0].data());
        blend_0x33_dl(r[1].data(), a[1].data(), b[1].data());
    }

    inline void perm_1032_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a)
    {
        perm_1032_dl(r[0].data(), a[0].data());
        perm_1032_dl(r[1].data(), a[1].data());
    }

    inline void perm_var_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const cc512i idx)
    {
        perm_var_dl(r[0].data(), a[0].data(), idx);
        perm_var_dl(r[1].data(), a[1].data(), idx);
    }

    inline void conv_sltohl(cc512i* r, const cc512i* a)
    {
        blend_0x55_hl(r, &a[VWORDS], a);
    }

    template<class RPackT, class APackT>
    inline void conv_sltohl(RPackT& r, const APackT& a)
    {
        conv_sltohl(r.data(), a.data());
    }

    inline void conv_dltovl(cc512i* r, const cc512i* a)
    {
        const cc512i m0 = VSET(0, FMASK, 0, FMASK, 0, FMASK, 0, FMASK);
        for (u64 i = 0; i < 2 * VWORDS; ++i) {
            r[i] = VAND(a[i], m0);
        }
        blend_0x55_hl(&r[2 * VWORDS], &a[3 * VWORDS], &a[2 * VWORDS]);
    }

    template<class RPackT, class APackT>
    inline void conv_dltovl(RPackT& r, const APackT& a)
    {
        conv_dltovl(r.data(), a.data());
    }

} // namespace

void cyclotomic_sqr_fp12_vec_v1(
    fp4_1x2x2x2w& ra,
    fp4_2x2x2x1w& rbc,
    const fp4_1x2x2x2w& a,
    const fp4_2x2x2x1w& bc)
{
    fp4_1x2x2x2w ta;
    fp4_2x2x2x1w tbc, t0;
    const cc512i m0 = VSET(3, 2, 1, 0, 5, 4, 7, 6);

    sqr_fp4_1x2x2x2w_v1(ta, a);
    as_fp2_2x2x2w(ra, ta, a);
    add_fp2_2x2x2w(ra, ra, ra);
    add_fp2_2x2x2w(ra, ra, ta);

    sqr_fp4_2x2x2x1w_v1(tbc, bc);
    mul_by_u_plus_1_fp2_4x2x1w(t0, tbc);
    blend_0xC0(tbc.data(), tbc.data(), t0.data());
    perm_var(tbc.data(), tbc.data(), m0);
    assa_fp2_4x2x1w(rbc, tbc, bc);
    add_fp2_4x2x1w(rbc, rbc, rbc);
    add_fp2_4x2x1w(rbc, rbc, tbc);
}

void cyclotomic_sqr_fp12_vec_v2(
    fp4_1x2x2x2w& ra,
    fp4_2x2x2x1w& rbc,
    const fp4_1x2x2x2w& a,
    const fp4_2x2x2x1w& bc)
{
    fp4_1x2x2x2w ta;
    fp4_2x2x2x1w tbc, t0;
    const cc512i m0 = VSET(3, 2, 1, 0, 5, 4, 7, 6);

    sqr_fp4_1x2x2x2w_v2(ta, a);
    as_fp2_2x2x2w(ra, ta, a);
    add_fp2_2x2x2w(ra, ra, ra);
    add_fp2_2x2x2w(ra, ra, ta);

    sqr_fp4_2x2x2x1w_v2(tbc, bc);
    mul_by_u_plus_1_fp2_4x2x1w(t0, tbc);
    blend_0xC0(tbc.data(), tbc.data(), t0.data());
    perm_var(tbc.data(), tbc.data(), m0);
    assa_fp2_4x2x1w(rbc, tbc, bc);
    add_fp2_4x2x1w(rbc, rbc, rbc);
    add_fp2_4x2x1w(rbc, rbc, tbc);
}

void compressed_cyclotomic_sqr_fp12_vec_v1(fp4_2x2x2x1w& rbc, const fp4_2x2x2x1w& bc)
{
    fp4_2x2x2x1w tbc, t0;
    const cc512i m0 = VSET(3, 2, 1, 0, 5, 4, 7, 6);

    sqr_fp4_2x2x2x1w_v1(tbc, bc);
    mul_by_u_plus_1_fp2_4x2x1w(t0, tbc);
    blend_0xC0(tbc.data(), tbc.data(), t0.data());
    perm_var(tbc.data(), tbc.data(), m0);
    assa_fp2_4x2x1w(rbc, tbc, bc);
    add_fp2_4x2x1w(rbc, rbc, rbc);
    add_fp2_4x2x1w(rbc, rbc, tbc);
}

void mul_fp12_vec_v1(
    fp2_4x2x1w& r01,
    fp2_4x2x1w& r2,
    const fp2_8x1x1w& ab0,
    const fp2_8x1x1w& ab1,
    const fp2_8x1x1w& ab2)
{
    fp2x2_8x1x1w tt01, tt2, tt3, tt4;
    fp2x2_4x2x1w ss0, ss1, ss2, ss3;
    fp2_4x2x1w s0, s1;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 2, 1, 0, 5, 4, 7, 6);

    mul_fp6x2_4x2x1x1w(tt01, tt2, ab0, ab1, ab2);
    shuf_01_fp2x2_8x1x1w(tt3, tt01);
    shuf_01_fp2x2_8x1x1w(tt4, tt2);
    blend_0x01_fp2x2_8x1x1w(tt4, tt4, tt3);
    blend_0x01_fp2x2_8x1x1w(tt2, tt2, tt01);
    blend_0x03_fp2x2_8x1x1w(tt01, tt01, tt4);
    shuf_01_fp2x2_8x1x1w(tt3, tt01);
    blend_0x55_dl(ss0.data(), tt3[1].data(), tt01[0].data());
    blend_0x55_dl(ss1.data(), tt01[1].data(), tt3[0].data());
    perm_1032_dl(ss2.data(), ss0.data());
    blend_0x33_dl(ss2.data(), ss1.data(), ss2.data());
    perm_var_dl(ss2.data(), ss2.data(), m0);
    perm_1032_dl(ss3.data(), ss1.data());
    blend_0x33_dl(ss3.data(), ss3.data(), ss0.data());
    add_fp2x2_4x2x1w(ss2, ss2, ss3);
    redc_fp2x2_4x2x1w(r01, ss2);

    shuf_01_dl(tt3[0].data(), tt2[1].data());
    blend_0x55_dl(ss0.data(), tt3[0].data(), tt2[0].data());
    redc_fp2x2_4x2x1w(s0, ss0);
    mul_by_u_plus_1_fp2_4x2x1w(s1, s0);
    blend_0xC0(s1.data(), s0.data(), s1.data());
    perm_var(s0.data(), s1.data(), m1);
    add_fp2_4x2x1w(r2, s0, s1);
}

void mul_fp12_vec_v2(
    fp2_4x2x1w& r01,
    fp2_2x2x2w& r2,
    const fp2_8x1x1w& ab0,
    const fp2_8x1x1w& ab1,
    const fp2_8x1x1w& ab2)
{
    fp2x2_8x1x1w tt01, tt2, tt3, tt4;
    fp2x2_4x2x1w ss0, ss1, ss2, ss3;
    fp2_4x2x1w s0, s1;
    fp2x2_2x2x2w hh0, hh1, hh2;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 3, 2, 2, 1, 1, 0, 0);
    const cc512i m2 = VSET(5, 5, 4, 4, 7, 7, 6, 6);

    mul_fp6x2_4x2x1x1w(tt01, tt2, ab0, ab1, ab2);
    shuf_01_fp2x2_8x1x1w(tt3, tt01);
    shuf_01_fp2x2_8x1x1w(tt4, tt2);
    blend_0x01_fp2x2_8x1x1w(tt4, tt4, tt3);
    blend_0x01_fp2x2_8x1x1w(tt2, tt2, tt01);
    blend_0x03_fp2x2_8x1x1w(tt01, tt01, tt4);
    shuf_01_fp2x2_8x1x1w(tt3, tt01);
    blend_0x55_dl(ss0.data(), tt3[1].data(), tt01[0].data());
    blend_0x55_dl(ss1.data(), tt01[1].data(), tt3[0].data());
    perm_1032_dl(ss2.data(), ss0.data());
    blend_0x33_dl(ss2.data(), ss1.data(), ss2.data());
    perm_var_dl(ss2.data(), ss2.data(), m0);
    perm_1032_dl(ss3.data(), ss1.data());
    blend_0x33_dl(ss3.data(), ss3.data(), ss0.data());
    add_fp2x2_4x2x1w(ss2, ss2, ss3);
    redc_fp2x2_4x2x1w(r01, ss2);

    shuf_01_dl(tt3[0].data(), tt2[1].data());
    blend_0x55_dl(ss0.data(), tt3[0].data(), tt2[0].data());
    perm_var_dl(ss1.data(), ss0.data(), m1);
    conv_dltovl(hh0, ss1);
    perm_var_dl(ss1.data(), ss0.data(), m2);
    conv_dltovl(hh1, ss1);
    mul_by_u_plus_1_fp2x2_2x2x2w(hh2, hh1);
    blend_0x0F_vl(hh1.data(), hh1.data(), hh2.data());
    add_fp2x2_2x2x2w(hh0, hh0, hh1);
    redc_fp2x2_2x2x2w(r2, hh0);
}

void mul_fp12_vec_v3(
    fp2_2x2x2w& r001,
    fp2_2x2x2w& r101,
    fp2_2x2x2w& r2,
    const fp2_4x2x1w& ab0,
    const fp2_4x2x1w& ab1,
    const fp2_4x2x1w& ab2)
{
    fp2x2_4x2x1w tt0, tt1, tt2, tt3, tt4, tt5;
    fp2_4x2x1w t0;
    fp2x2_2x2x2w ss0, ss1, ss2, ss3, ss4, ss5;
    fp2_2x2x2w s0, s1, s2, s3;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 3, 2, 2, 1, 1, 0, 0);
    const cc512i m2 = VSET(7, 7, 6, 6, 7, 7, 6, 6);
    const cc512i m3 = VSET(7, 7, 6, 6, 5, 5, 4, 4);

    mul_fp6x2_2x2x2x1w(tt0, tt1, ab0, ab1, ab2);

    perm_1100_dl(tt2.data(), tt0.data());
    perm_1100_dl(tt1.data(), tt1.data());
    perm_var_dl(tt3.data(), tt1.data(), m0);
    blend_0x0F_dl(tt4.data(), tt2.data(), tt3.data());
    perm_var_dl(tt5.data(), tt0.data(), m1);
    conv_dltovl(ss0, tt4);
    mul_by_u_plus_1_fp2x2_2x2x2w(ss1, ss0);
    blend_0x0F_vl(ss1.data(), ss0.data(), ss1.data());
    conv_dltovl(ss4, tt5);
    add_fp2x2_2x2x2w(ss1, ss4, ss1);
    redc_fp2x2_2x2x2w(r001, ss1);
    conv_dltovl(ss2, tt1);
    perm_var_dl(tt5.data(), tt0.data(), m2);
    conv_dltovl(ss1, tt5);

    perm_var(t0.data(), ab0.data(), m1);
    conv_sltohl(s0, t0);
    perm_var(t0.data(), ab0.data(), m3);
    conv_sltohl(s1, t0);
    add_fp2_2x2x2w(s0, s0, s1);
    perm_var(t0.data(), ab1.data(), m1);
    conv_sltohl(s1, t0);
    perm_var(t0.data(), ab1.data(), m3);
    conv_sltohl(s2, t0);
    add_fp2_2x2x2w(s1, s1, s2);
    perm_var(t0.data(), ab2.data(), m1);
    conv_sltohl(s2, t0);
    perm_var(t0.data(), ab2.data(), m3);
    conv_sltohl(s3, t0);
    add_fp2_2x2x2w(s2, s2, s3);
    mul_fp6x2_1x2x2x2w(ss3, ss5, s0, s1, s2);
    sub_fp2x2_2x2x2w(ss3, ss3, ss4);
    perm_var_dl(tt5.data(), tt0.data(), m3);
    conv_dltovl(ss4, tt5);
    sub_fp2x2_2x2x2w(ss3, ss3, ss4);
    redc_fp2x2_2x2x2w(r101, ss3);
    sub_fp2x2_2x2x2w(ss5, ss5, ss2);
    add_fp2x2_2x2x2w(ss1, ss2, ss1);
    sub_fp2x2_2x2x2w(ss5, ss5, ss0);
    perm_var_vl(ss5.data(), ss5.data(), m0);
    blend_0x0F_vl(ss1.data(), ss5.data(), ss1.data());
    redc_fp2x2_2x2x2w(r2, ss1);
}

void mul_fp12_vec_v4(
    fp2_4x2x1w& r0,
    fp2_2x2x2w& r101,
    fp2_2x2x2w& r12,
    const fp2_4x2x1w& ab0,
    const fp2_4x2x1w& ab1,
    const fp2_4x2x1w& ab2)
{
    fp2x2_4x2x1w tt0, tt1, tt2, tt3, tt4, tt5;
    fp2_4x2x1w t0;
    fp2x2_2x2x2w ss0, ss1, ss2, ss3, ss4;
    fp2_2x2x2w s0, s1, s2, s3;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 3, 2, 2, 1, 1, 0, 0);
    const cc512i m2 = VSET(1, 0, 7, 6, 5, 4, 3, 2);
    const cc512i m3 = VSET(7, 7, 6, 6, 5, 5, 4, 4);

    mul_fp6x2_2x2x2x1w(tt0, tt1, ab0, ab1, ab2);

    perm_var_dl(tt2.data(), tt1.data(), m0);
    blend_0x0F_dl(tt3.data(), tt2.data(), tt0.data());
    mul_by_u_plus_1_fp2x2_4x2x1w(tt4, tt2);
    perm_var_dl(tt2.data(), tt0.data(), m2);
    blend_0x03_dl(tt2.data(), tt2.data(), tt4.data());
    add_fp2x2_4x2x1w(tt3, tt3, tt2);
    redc_fp2x2_4x2x1w(r0, tt3);

    perm_1100_dl(tt2.data(), tt0.data());
    perm_1100_dl(tt1.data(), tt1.data());
    perm_var_dl(tt3.data(), tt1.data(), m0);
    blend_0x0F_dl(tt4.data(), tt2.data(), tt3.data());
    perm_var_dl(tt5.data(), tt0.data(), m1);
    conv_dltovl(ss0, tt4);
    conv_dltovl(ss4, tt5);
    conv_dltovl(ss2, tt1);

    perm_var(t0.data(), ab0.data(), m1);
    conv_sltohl(s0, t0);
    perm_var(t0.data(), ab0.data(), m3);
    conv_sltohl(s1, t0);
    add_fp2_2x2x2w(s0, s0, s1);
    perm_var(t0.data(), ab1.data(), m1);
    conv_sltohl(s1, t0);
    perm_var(t0.data(), ab1.data(), m3);
    conv_sltohl(s2, t0);
    add_fp2_2x2x2w(s1, s1, s2);
    perm_var(t0.data(), ab2.data(), m1);
    conv_sltohl(s2, t0);
    perm_var(t0.data(), ab2.data(), m3);
    conv_sltohl(s3, t0);
    add_fp2_2x2x2w(s2, s2, s3);
    mul_fp6x2_1x2x2x2w(ss3, ss1, s0, s1, s2);
    sub_fp2x2_2x2x2w(ss3, ss3, ss4);
    perm_var_dl(tt5.data(), tt0.data(), m3);
    conv_dltovl(ss4, tt5);
    sub_fp2x2_2x2x2w(ss3, ss3, ss4);
    redc_fp2x2_2x2x2w(r101, ss3);
    sub_fp2x2_2x2x2w(ss1, ss1, ss2);
    sub_fp2x2_2x2x2w(ss1, ss1, ss0);
    redc_fp2x2_2x2x2w(r12, ss1);
}

void mul_by_xy00z0_fp12_vec_v1(
    fp2_4x2x1w& r0,
    fp2_2x2x2w& r1,
    const fp2_4x2x1w& a01,
    const fp2_4x2x1w& a2,
    const fp2_4x2x1w& b01,
    const fp2_4x2x1w& b4)
{
    fp2x2_4x2x1w tt0, tt1, tt2, tt3, tt4, tt5;
    fp2x2_2x2x2w ss0, ss1, ss2;
    const cc512i m0 = VSET(1, 0, 3, 2, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 3, 2, 2, 3, 3, 2, 2);
    const cc512i m2 = VSET(7, 7, 6, 6, 7, 7, 6, 6);

    mul_by_xy00z0_fp6x2_2x2x2x1w(tt0, tt1, tt2, tt3, a01, a2, b01, b4);

    blend_0x03_dl(tt5.data(), tt0.data(), tt1.data());
    blend_0x33_dl(tt4.data(), tt2.data(), tt3.data());
    perm_var_dl(tt4.data(), tt4.data(), m0);
    add_fp2x2_4x2x1w(tt5, tt5, tt4);
    redc_fp2x2_4x2x1w(r0, tt5);

    perm_1100_dl(tt0.data(), tt0.data());
    perm_var_dl(tt4.data(), tt3.data(), m1);
    blend_0x0F_dl(tt4.data(), tt4.data(), tt0.data());
    conv_dltovl(ss0, tt4);
    perm_1100_dl(tt1.data(), tt1.data());
    perm_var_dl(tt3.data(), tt3.data(), m2);
    blend_0x0F_dl(tt1.data(), tt1.data(), tt3.data());
    conv_dltovl(ss1, tt1);
    mul_by_u_plus_1_fp2x2_2x2x2w(ss2, ss1);
    blend_0x0F_vl(ss1.data(), ss1.data(), ss2.data());
    add_fp2x2_2x2x2w(ss1, ss1, ss0);
    redc_fp2x2_2x2x2w(r1, ss1);
}

void mul_by_xy00z0_fp12_vec_v2(
    fp2_4x2x1w& r0,
    fp2_2x2x2w& r1,
    const fp2_8x1x1w& a,
    const fp2_8x1x1w& b)
{
    fp2x2_8x1x1w tt0, tt1, tt2, tt3;
    fp2x2_4x2x1w ss0, ss1;
    fp2x2_2x2x2w hh0, hh1, hh2;
    const cc512i m0 = VSET(6, 6, 7, 7, 2, 2, 1, 1);
    const cc512i m1 = VSET(0, 0, 2, 2, 6, 6, 4, 4);
    const cc512i m2 = VSET(1, 1, 1, 1, 3, 3, 3, 3);
    const cc512i m3 = VSET(5, 5, 5, 5, 7, 7, 7, 7);

    mul_by_xy00z0_fp6x2_2x4x1x1w(tt0, tt1, tt2, a, b);

    blend_0x33_fp2x2_8x1x1w(tt3, tt0, tt1);
    perm_var_fp2x2_8x1x1w(tt3, tt3, m0);
    blend_0x55_dl(ss0.data(), tt3[1].data(), tt3[0].data());
    blend_0x33_fp2x2_8x1x1w(tt3, tt1, tt2);
    perm_var_fp2x2_8x1x1w(tt3, tt3, m1);
    blend_0x55_dl(ss1.data(), tt3[1].data(), tt3[0].data());
    add_fp2x2_4x2x1w(ss0, ss0, ss1);
    redc_fp2x2_4x2x1w(r0, ss0);

    blend_0x33_fp2x2_8x1x1w(tt3, tt0, tt2);
    perm_var_fp2x2_8x1x1w(tt3, tt3, m2);
    blend_0x33_dl(ss0.data(), tt3[1].data(), tt3[0].data());
    conv_dltovl(hh0, ss0);
    perm_1032_fp2x2_8x1x1w(tt3, tt2);
    blend_0x33_fp2x2_8x1x1w(tt3, tt3, tt1);
    perm_var_fp2x2_8x1x1w(tt3, tt3, m3);
    blend_0x33_dl(ss1.data(), tt3[1].data(), tt3[0].data());
    conv_dltovl(hh1, ss1);
    mul_by_u_plus_1_fp2x2_2x2x2w(hh2, hh1);
    blend_0x0F_vl(hh1.data(), hh1.data(), hh2.data());
    add_fp2x2_2x2x2w(hh1, hh1, hh0);
    redc_fp2x2_2x2x2w(r1, hh1);
}

void sqr_fp12_vec_v1(
    fp2_4x2x1w& r0,
    fp2_4x2x1w& r1,
    const fp2_4x2x1w& a0,
    const fp2_4x2x1w& a1)
{
    fp2_4x2x1w t0, t1, t2, t3, t4;
    fp2_4x2x1w s0, s1, s2;
    fp2x2_4x2x1w tt0, tt1;
    const cc512i m0 = VSET(5, 4, 3, 2, 1, 0, 7, 6);
    const cc512i m1 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m2 = VSET(5, 4, 1, 0, 7, 6, 3, 2);
    const cc512i m3 = VSET(7, 6, 3, 2, 5, 4, 1, 0);
    const cc512i m4 = VSET(4, 4, 1, 0, 7, 6, 3, 2);
    const cc512i m5 = VSET(5, 5, 5, 4, 7, 6, 3, 2);

    FpPack8x1::asax6(t0, a0, a1);
    blend_0x0F(t1.data(), t0.data(), a1.data());
    perm_var(t1.data(), t1.data(), m0);
    add_fp2_4x2x1w(t1, a0, t1);
    perm_var(t2.data(), a1.data(), m1);
    blend_0x0F(t2.data(), t2.data(), a0.data());
    perm_var(s0.data(), t2.data(), m2);
    perm_1032(t3.data(), t0.data());
    blend_0x03(t3.data(), t3.data(), t1.data());
    blend_0x0F(s0.data(), s0.data(), t3.data());
    perm_var(s1.data(), t2.data(), m3);
    perm_1032(t3.data(), t1.data());
    blend_0x03(t3.data(), t0.data(), t3.data());
    blend_0x0F(s1.data(), s1.data(), t3.data());
    perm_1032(s2.data(), a1.data());
    blend_0xC0(s2.data(), a0.data(), s2.data());
    perm_var(t2.data(), t0.data(), m3);
    perm_var(t3.data(), t1.data(), m1);
    blend_0x03(t2.data(), t2.data(), t3.data());
    blend_0x0F(s2.data(), s2.data(), t2.data());

    mul_fp6x2_2x2x2x1w(tt0, tt1, s0, s1, s2);
    redc_fp2x2_4x2x1w(t0, tt0);
    redc_fp2x2_4x2x1w(t1, tt1);

    perm_var(t2.data(), t0.data(), m1);
    blend_0x0F(t2.data(), t1.data(), t2.data());
    add_fp2_4x2x1w(r1, t2, t2);

    perm_var(t4.data(), t1.data(), m4);
    blend_0x0F(t4.data(), t4.data(), t0.data());
    perm_var(t3.data(), t1.data(), m5);
    blend_0x0F(t3.data(), t3.data(), t2.data());

    FpPack8x1::assx6(t3, t4, t3);
    blend_0x0F(t4.data(), t3.data(), t2.data());
    perm_var(t4.data(), t4.data(), m0);
    sub_fp2_4x2x1w(r0, t3, t4);
}

void mul_fp12_vector(Fp12& ret, const Fp12& a, const Fp12& b)
{
    fp2_4x2x1w ab0, ab1, ab2;
    fp2_2x2x2w r001, r101, r2;
    std::array<cc512i, SWORDS> in0, in1, in2;
    std::array<cc512i, SWORDS / 2> out0, out1, out2;

    for (u64 i = 0; i < SWORDS; ++i) {
        in0[i] = VSET(
            static_cast<i64>(b[1][0][1][i]), static_cast<i64>(b[1][0][0][i]),
            static_cast<i64>(a[1][0][1][i]), static_cast<i64>(a[1][0][0][i]),
            static_cast<i64>(b[0][0][1][i]), static_cast<i64>(b[0][0][0][i]),
            static_cast<i64>(a[0][0][1][i]), static_cast<i64>(a[0][0][0][i]));
        in1[i] = VSET(
            static_cast<i64>(b[1][1][1][i]), static_cast<i64>(b[1][1][0][i]),
            static_cast<i64>(a[1][1][1][i]), static_cast<i64>(a[1][1][0][i]),
            static_cast<i64>(b[0][1][1][i]), static_cast<i64>(b[0][1][0][i]),
            static_cast<i64>(a[0][1][1][i]), static_cast<i64>(a[0][1][0][i]));
        in2[i] = VSET(
            static_cast<i64>(b[1][2][1][i]), static_cast<i64>(b[1][2][0][i]),
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(b[0][2][1][i]), static_cast<i64>(b[0][2][0][i]),
            static_cast<i64>(a[0][2][1][i]), static_cast<i64>(a[0][2][0][i]));
    }

    FpPack8x1::conv_64to48(ab0, in0.data());
    FpPack8x1::conv_64to48(ab1, in1.data());
    FpPack8x1::conv_64to48(ab2, in2.data());

    mul_fp12_vec_v3(r001, r101, r2, ab0, ab1, ab2);

    FpPack4x2::carryp(r001);
    FpPack4x2::carryp(r101);
    FpPack4x2::carryp(r2);

    FpPack4x2::conv_48to64(out0.data(), r001);
    FpPack4x2::conv_48to64(out1.data(), r101);
    FpPack4x2::conv_48to64(out2.data(), r2);

    for (u64 i = 0; i < SWORDS / 2; ++i) {
        const u64* lanes0 = reinterpret_cast<const u64*>(&out0[i]);
        const u64* lanes1 = reinterpret_cast<const u64*>(&out1[i]);
        const u64* lanes2 = reinterpret_cast<const u64*>(&out2[i]);

        ret[0][0][0][i] = lanes0[0];
        ret[0][0][0][i + (SWORDS / 2)] = lanes0[1];
        ret[0][0][1][i] = lanes0[2];
        ret[0][0][1][i + (SWORDS / 2)] = lanes0[3];
        ret[0][1][0][i] = lanes0[4];
        ret[0][1][0][i + (SWORDS / 2)] = lanes0[5];
        ret[0][1][1][i] = lanes0[6];
        ret[0][1][1][i + (SWORDS / 2)] = lanes0[7];

        ret[1][0][0][i] = lanes1[0];
        ret[1][0][0][i + (SWORDS / 2)] = lanes1[1];
        ret[1][0][1][i] = lanes1[2];
        ret[1][0][1][i + (SWORDS / 2)] = lanes1[3];
        ret[1][1][0][i] = lanes1[4];
        ret[1][1][0][i + (SWORDS / 2)] = lanes1[5];
        ret[1][1][1][i] = lanes1[6];
        ret[1][1][1][i + (SWORDS / 2)] = lanes1[7];

        ret[0][2][0][i] = lanes2[0];
        ret[0][2][0][i + (SWORDS / 2)] = lanes2[1];
        ret[0][2][1][i] = lanes2[2];
        ret[0][2][1][i + (SWORDS / 2)] = lanes2[3];
        ret[1][2][0][i] = lanes2[4];
        ret[1][2][0][i + (SWORDS / 2)] = lanes2[5];
        ret[1][2][1][i] = lanes2[6];
        ret[1][2][1][i + (SWORDS / 2)] = lanes2[7];
    }
}

void mul_by_xy00z0_fp12_vector(Fp12& ret, const Fp12& a, const Fp6& xy00z0)
{
    fp2_4x2x1w r0, a01, a2, b01, b4;
    fp2_2x2x2w r1;
    std::array<cc512i, SWORDS> in0, in1, in2, in3;
    std::array<cc512i, SWORDS> out0;
    std::array<cc512i, SWORDS / 2> out1;

    for (u64 i = 0; i < SWORDS; ++i) {
        in0[i] = VSET(
            static_cast<i64>(a[1][1][1][i]), static_cast<i64>(a[1][1][0][i]),
            static_cast<i64>(a[1][0][1][i]), static_cast<i64>(a[1][0][0][i]),
            static_cast<i64>(a[0][1][1][i]), static_cast<i64>(a[0][1][0][i]),
            static_cast<i64>(a[0][0][1][i]), static_cast<i64>(a[0][0][0][i]));
        in1[i] = VSET(
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(a[0][2][1][i]), static_cast<i64>(a[0][2][0][i]),
            static_cast<i64>(a[0][2][1][i]), static_cast<i64>(a[0][2][0][i]));
        in2[i] = VSET(
            static_cast<i64>(xy00z0[1][1][i]), static_cast<i64>(xy00z0[1][0][i]),
            static_cast<i64>(xy00z0[0][1][i]), static_cast<i64>(xy00z0[0][0][i]),
            static_cast<i64>(xy00z0[1][1][i]), static_cast<i64>(xy00z0[1][0][i]),
            static_cast<i64>(xy00z0[0][1][i]), static_cast<i64>(xy00z0[0][0][i]));
        in3[i] = VSET(
            static_cast<i64>(xy00z0[2][1][i]), static_cast<i64>(xy00z0[2][0][i]),
            static_cast<i64>(xy00z0[2][1][i]), static_cast<i64>(xy00z0[2][0][i]),
            static_cast<i64>(xy00z0[2][1][i]), static_cast<i64>(xy00z0[2][0][i]),
            static_cast<i64>(xy00z0[2][1][i]), static_cast<i64>(xy00z0[2][0][i]));
    }

    FpPack8x1::conv_64to48(a01, in0.data());
    FpPack8x1::conv_64to48(a2, in1.data());
    FpPack8x1::conv_64to48(b01, in2.data());
    FpPack8x1::conv_64to48(b4, in3.data());

    mul_by_xy00z0_fp12_vec_v1(r0, r1, a01, a2, b01, b4);

    FpPack4x2::carryp(r1);

    FpPack8x1::conv_48to64(out0.data(), r0);
    FpPack4x2::conv_48to64(out1.data(), r1);

    for (u64 i = 0; i < SWORDS; ++i) {
        const u64* lanes = reinterpret_cast<const u64*>(&out0[i]);
        ret[0][2][0][i] = lanes[0];
        ret[0][2][1][i] = lanes[1];
        ret[0][1][0][i] = lanes[2];
        ret[0][1][1][i] = lanes[3];
        ret[1][0][0][i] = lanes[4];
        ret[1][0][1][i] = lanes[5];
        ret[1][1][0][i] = lanes[6];
        ret[1][1][1][i] = lanes[7];
    }

    for (u64 i = 0; i < SWORDS / 2; ++i) {
        const u64* lanes = reinterpret_cast<const u64*>(&out1[i]);
        ret[0][0][0][i] = lanes[0];
        ret[0][0][0][i + (SWORDS / 2)] = lanes[1];
        ret[0][0][1][i] = lanes[2];
        ret[0][0][1][i + (SWORDS / 2)] = lanes[3];
        ret[1][2][0][i] = lanes[4];
        ret[1][2][0][i + (SWORDS / 2)] = lanes[5];
        ret[1][2][1][i] = lanes[6];
        ret[1][2][1][i + (SWORDS / 2)] = lanes[7];
    }
}

void sqr_fp12_vector(Fp12& ret, const Fp12& a)
{
    fp2_4x2x1w r0, r1, a0, a1;
    std::array<cc512i, SWORDS> in0, in1, out0, out1;

    for (u64 i = 0; i < SWORDS; ++i) {
        in0[i] = VSET(
            static_cast<i64>(a[1][2][0][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(a[0][2][1][i]), static_cast<i64>(a[0][2][0][i]),
            static_cast<i64>(a[0][1][1][i]), static_cast<i64>(a[0][1][0][i]),
            static_cast<i64>(a[0][0][1][i]), static_cast<i64>(a[0][0][0][i]));
        in1[i] = VSET(
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][1][i]),
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(a[1][1][1][i]), static_cast<i64>(a[1][1][0][i]),
            static_cast<i64>(a[1][0][1][i]), static_cast<i64>(a[1][0][0][i]));
    }

    FpPack8x1::conv_64to48(a0, in0.data());
    FpPack8x1::conv_64to48(a1, in1.data());

    sqr_fp12_vec_v1(r0, r1, a0, a1);

    FpPack8x1::conv_48to64(out0.data(), r0);
    FpPack8x1::conv_48to64(out1.data(), r1);

    for (u64 i = 0; i < SWORDS; ++i) {
        const u64* lanes0 = reinterpret_cast<const u64*>(&out0[i]);
        const u64* lanes1 = reinterpret_cast<const u64*>(&out1[i]);
        ret[0][0][0][i] = lanes0[0];
        ret[0][0][1][i] = lanes0[1];
        ret[0][1][0][i] = lanes0[2];
        ret[0][1][1][i] = lanes0[3];
        ret[0][2][0][i] = lanes0[4];
        ret[0][2][1][i] = lanes0[5];
        ret[1][0][0][i] = lanes1[0];
        ret[1][0][1][i] = lanes1[1];
        ret[1][1][0][i] = lanes1[2];
        ret[1][1][1][i] = lanes1[3];
        ret[1][2][0][i] = lanes1[4];
        ret[1][2][1][i] = lanes1[5];
    }
}

void cyclotomic_sqr_fp12_vector(Fp12& ret, const Fp12& a)
{
    fp4_1x2x2x2w packed_a, out_a;
    fp4_2x2x2x1w packed_bc, out_bc;
    std::array<cc512i, SWORDS / 2> lanes_a;
    std::array<cc512i, SWORDS> lanes_bc;

    for (u64 i = 0; i < SWORDS / 2; ++i) {
        lanes_a[i] = VSET(
            static_cast<i64>(a[1][1][1][i + (SWORDS / 2)]), static_cast<i64>(a[1][1][1][i]),
            static_cast<i64>(a[1][1][0][i + (SWORDS / 2)]), static_cast<i64>(a[1][1][0][i]),
            static_cast<i64>(a[0][0][1][i + (SWORDS / 2)]), static_cast<i64>(a[0][0][1][i]),
            static_cast<i64>(a[0][0][0][i + (SWORDS / 2)]), static_cast<i64>(a[0][0][0][i]));
    }
    FpPack4x2::conv_64to48(packed_a, lanes_a.data());

    for (u64 i = 0; i < SWORDS; ++i) {
        lanes_bc[i] = VSET(
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(a[0][1][1][i]), static_cast<i64>(a[0][1][0][i]),
            static_cast<i64>(a[0][2][1][i]), static_cast<i64>(a[0][2][0][i]),
            static_cast<i64>(a[1][0][1][i]), static_cast<i64>(a[1][0][0][i]));
    }
    FpPack8x1::conv_64to48(packed_bc, lanes_bc.data());

    cyclotomic_sqr_fp12_vec_v1(out_a, out_bc, packed_a, packed_bc);

    FpPack4x2::carryp(out_a);
    FpPack4x2::conv_48to64(lanes_a.data(), out_a);
    for (u64 i = 0; i < SWORDS / 2; ++i) {
        const u64* lanes = reinterpret_cast<const u64*>(&lanes_a[i]);
        ret[0][0][0][i] = lanes[0];
        ret[0][0][0][i + (SWORDS / 2)] = lanes[1];
        ret[0][0][1][i] = lanes[2];
        ret[0][0][1][i + (SWORDS / 2)] = lanes[3];
        ret[1][1][0][i] = lanes[4];
        ret[1][1][0][i + (SWORDS / 2)] = lanes[5];
        ret[1][1][1][i] = lanes[6];
        ret[1][1][1][i + (SWORDS / 2)] = lanes[7];
    }

    FpPack8x1::conv_48to64(lanes_bc.data(), out_bc);
    for (u64 i = 0; i < SWORDS; ++i) {
        const u64* lanes = reinterpret_cast<const u64*>(&lanes_bc[i]);
        ret[1][0][0][i] = lanes[0];
        ret[1][0][1][i] = lanes[1];
        ret[0][2][0][i] = lanes[2];
        ret[0][2][1][i] = lanes[3];
        ret[0][1][0][i] = lanes[4];
        ret[0][1][1][i] = lanes[5];
        ret[1][2][0][i] = lanes[6];
        ret[1][2][1][i] = lanes[7];
    }
}

} // namespace blst

#endif
