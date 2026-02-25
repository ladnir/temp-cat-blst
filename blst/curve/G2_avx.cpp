#include "blst/curve/G2_avx.h"
#include "blst/backend/avx512/avx_layout.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{

void line_dbl_vec_v1(
    fp2_2x2x2w l0Y3,
    fp2_2x2x2w l1,
    fp2_2x2x2w l2,
    fp2_2x2x2w X3,
    fp2_2x2x2w Z3,
    const fp2_2x2x2w X1Y1,
    const fp2_2x2x2w Z1)
{
    fp2_2x2x2w AB, C, D, E, FZZ;
    fp2_2x2x2w t0, t1, t2, t3, t4, t5;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 2, 1, 0, 3, 2, 1, 0);

    sqr_fp2_2x2x2w(AB, X1Y1);
    perm_var_hl(t3, X1Y1, m0);
    blend_0x0F_hl(t0, t3, AB);
    add_fp2_2x2x2w(t0, t0, AB);
    perm_var_hl(t1, AB, m0);
    blend_0x0F_hl(t2, t0, t1);
    sqr_fp2_2x2x2w(C, t2);
    blend_0x0F_hl(t2, C, t0);
    blend_0x0F_hl(t0, t1, AB);
    sa_fp2_2x2x2w(E, t2, t0);
    perm_var_hl(t5, C, m0);
    blend_0x0F_hl(t2, t5, X1Y1);
    sa_fp2_2x2x2w(t0, E, t2);
    blend_0x0F_hl(t2, Z1, E);
    sqr_fp2_2x2x2w(FZZ, t2);
    blend_0x0F_hl(t3, t0, t3);
    blend_0x0F_hl(t2, t0, Z1);
    add_fp2_2x2x2w(D, t3, t2);
    perm_var_hl(t0, t0, m0);
    blend_0x0F_hl(t0, t0, D);
    sqr_fp2_2x2x2w(t0, t0);
    perm_var_hl(t2, FZZ, m0);
    blend_0x0F_hl(t3, t2, t0);
    blend_0x0F_hl(t1, D, t1);
    sub_fp2_2x2x2w(t3, t3, t1);
    perm_var_hl(t0, t0, m0);
    blend_0x0F_hl(t0, t3, t0);
    blend_0x0F_hl(t4, D, AB);
    sub_fp2_2x2x2w(X3, t0, t4);
    perm_var_hl(t3, t3, m0);
    blend_0x0F_hl(t3, t3, C);
    blend_0x0F_hl(t4, FZZ, C);
    sa_fp2_2x2x2w(Z3, t3, t4);
    perm_var_hl(t4, Z3, m0);
    blend_0x0F_hl(t1, t4, t1);
    add_fp2_2x2x2w(t1, t1, t1);
    add_fp2_2x2x2w(t1, t1, t1);
    blend_0x0F_hl(t0, D, X3);
    blend_0x0F_hl(t3, X3, FZZ);
    sub_fp2_2x2x2w(t0, t0, t3);
    blend_0x0F_hl(t2, t0, t2);
    perm_var_hl(t3, E, m1);
    mul_fp2_2x2x2w(l1, t2, t3);
    blend_0x0F_hl(t0, l1, t0);
    sub_fp2_2x2x2w(l0Y3, t0, t1);
    mul_fp2_2x2x2w(l2, Z3, FZZ);
}

void line_dbl_vec_v2(
    fp2_4x2x1w l0,
    fp2_4x2x1w l12,
    fp2_4x2x1w X3,
    fp2_4x2x1w Y3,
    fp2_4x2x1w Z3,
    const fp2_4x2x1w X1Y1Z1)
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

    perm_1032(t0, X1Y1Z1);
    add_fp2_4x2x1w(t0, X1Y1Z1, t0);
    blend_0x03(t0, X1Y1Z1, t0);
    sqr_fp2_4x2x1w(ABZZ, t0);
    perm_var(t0, ABZZ, m0);
    perm_var(t1, X1Y1Z1, m1);
    blend_0xC3(t1, t0, t1);
    add_fp2_4x2x1w(t0, t0, t1);
    perm_1032(t1, t0);
    blend_0xC3(t1, t0, t1);
    perm_var(t2, t0, m2);
    blend_0x33(t2, t2, ABZZ);
    add_fp2_4x2x1w(E, t1, t2);
    perm_var(t2, E, m3);
    blend_0x0F(t2, t0, t2);
    perm_var(t0, ABZZ, m4);
    blend_0x30(t2, t2, t0);
    sqr_fp2_4x2x1w(CF, t2);
    blend_0x03(t2, CF, ABZZ);
    blend_0x0F(t0, CF, t0);
    aass_fp2_4x2x1w(t0, t2, t0);
    blend_0x03(t2, ABZZ, CF);
    perm_1032(t2, t2);
    perm_1010(t3, t0);
    blend_0x0F(t2, t3, t2);
    sass_fp2_4x2x1w(Z3, t0, t2);
    blend_0xC0(t1, Z3, t1);
    blend_0x0F(t1, t1, E);
    sass_fp2_4x2x1w(l0, Z3, t1);
    blend_0x0F(t0, l0, CF);
    perm_var(t1, t0, m5);
    perm_var(t0, t0, m6);
    ssaa_fp2_4x2x1w(t0, t1, t0);
    blend_0x0F(t1, l0, t0);
    perm_var(t1, t1, m7);
    sass_fp2_4x2x1w(X3, t0, t1);
    blend_0x0F(t0, X3, Z3);
    perm_1032(t0, t0);
    blend_0x33(t0, t0, E);
    blend_0x0F(t1, E, ABZZ);
    perm_var(t1, t1, m8);
    mul_fp2_4x2x1w(l12, t0, t1);
    perm_1010(t0, l0);
    sub_fp2_4x2x1w(Y3, l12, t0);
}

void line_add_vec_v1(
    fp2_2x2x2w l0Y3,
    fp2_2x2x2w l1,
    fp2_2x2x2w X3,
    fp2_2x2x2w Z3,
    const fp2_2x2x2w X1Y1,
    const fp2_2x2x2w Z1Y2,
    const fp2_2x2x2w X2)
{
    fp2_2x2x2w H, J, ZZ, US, V;
    fp2_2x2x2w t0, t1, t2, t3, t4;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(7, 6, 5, 4, 7, 6, 5, 4);

    perm_var_hl(t0, Z1Y2, m1);
    mul_fp2_2x2x2w(ZZ, t0, Z1Y2);
    blend_0x0F_hl(t0, X2, ZZ);
    perm_var_hl(t1, ZZ, m1);
    mul_fp2_2x2x2w(US, t0, t1);
    sub_fp2_2x2x2w(H, US, X1Y1);
    blend_0x0F_hl(t0, Z1Y2, H);
    add_fp2_2x2x2w(l1, H, t0);
    perm_var_hl(t0, H, m0);
    blend_0x0F_hl(t0, l1, t0);
    sqr_fp2_2x2x2w(t0, t0);
    blend_0x0F_hl(t1, ZZ, t0);
    sa_fp2_2x2x2w(t0, t0, t1);
    perm_var_hl(t1, t1, m0);
    blend_0x0F_hl(t1, t1, t0);
    sa_fp2_2x2x2w(Z3, t0, t1);
    blend_0x0F_hl(t0, H, l1);
    perm_var_hl(t3, Z3, m0);
    blend_0x0F_hl(t1, t3, X2);
    mul_fp2_2x2x2w(J, t0, t1);
    perm_var_hl(t1, l1, m0);
    blend_0x0F_hl(t2, t1, Z1Y2);
    blend_0x0F_hl(t0, t1, t3);
    mul_fp2_2x2x2w(t0, t0, t2);
    perm_var_hl(t2, J, m0);
    blend_0x0F_hl(t2, t3, t2);
    mul_fp2_2x2x2w(V, X1Y1, t2);
    perm_var_hl(t3, V, m0);
    blend_0x0F_hl(t2, t0, t3);
    blend_0x0F_hl(t3, J, t3);
    sa_fp2_2x2x2w(t2, t2, t3);
    blend_0x0F_hl(t3, t2, V);
    perm_var_hl(t2, t2, m0);
    blend_0x0F_hl(t2, t2, t3);
    sa_fp2_2x2x2w(X3, t3, t2);
    blend_0x0F_hl(t3, V, J);
    blend_0x0F_hl(t0, X3, t0);
    sub_fp2_2x2x2w(t0, t3, t0);
    mul_fp2_2x2x2w(t2, t0, t1);
    blend_0x0F_hl(t0, t2, t0);
    perm_var_hl(t1, X3, m0);
    blend_0x0F_hl(t1, t1, t0);
    sa_fp2_2x2x2w(l0Y3, t0, t1);
}

void line_by_Px2_4x2x1w(fp2_4x2x1w r, const fp2_4x2x1w a, const fp2_4x2x1w Px2)
{
    fp2x2_4x2x1w tt0;
    mul_fpx2_8x1w(tt0, a, Px2);
    redc_fpx2_8x1w(r, tt0);
}

void line_by_Px2_2x2x2w(fp2_2x2x2w r, const fp2_2x2x2w a, const fp2_2x2x2w Px2)
{
    fp2x2_2x2x2w tt0;
    mul_fpx2_4x2w(tt0, a, Px2);
    redc_fpx2_4x2w(r, tt0);
}

} // namespace blst

#endif
