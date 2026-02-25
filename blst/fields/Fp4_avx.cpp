#include "blst/fields/Fp4_avx.h"
#include "blst/backend/avx512/avx_layout.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{

void sqr_fp4_2x2x2x1w_v1(fp4_2x2x2x1w& r, const fp4_2x2x2x1w& a)
{
    fp4_2x2x2x1w t0;
    fp4x2_2x2x2x1w tt0, tt1, tt2;

    // a = A1 | A0 at Fp2 layer
    sqr_fp2x2_4x2x1w(tt0, a);               //        A1^2 |              A0^2
    mul_by_u_plus_1_fp2x2_4x2x1w(tt1, tt0); //  (u+1)*A1^2 |               ...
    perm_zz32_dl(tt1.data(), tt1.data());   //           0 |        (u+1)*A1^2
    perm_10zz(t0.data(), a.data());         //          A0 |                 0
    mul_fp2x2_2x4x1w(tt2, a, t0);           //       A0*A1 |               ...
    blend_0x33_dl(tt0.data(), tt2.data(), tt0.data()); //       A0*A1 |              A0^2
    blend_0x33_dl(tt2.data(), tt2.data(), tt1.data()); //       A0*A1 |      (u+1)*(A1^2)
    add_fp2x2_4x2x1w(tt0, tt0, tt2);        //     2*A0*A1 | A0^2+(u+1)*(A1^2)
    redc_fp2x2_4x2x1w(r, tt0);              //     2*A0*A1 | A0^2+(u+1)*(A1^2)
}

void sqr_fp4_2x2x2x1w_v2(fp4_2x2x2x1w& r, const fp4_2x2x2x1w& a)
{
    fp4_2x2x2x1w t0, t1, t2, t3;

    // a = A1 | A0 at Fp2 layer
    sqr_fp2_4x2x1w(t0, a);               //        A1^2 |              A0^2
    mul_by_u_plus_1_fp2_4x2x1w(t1, t0);  //  (u+1)*A1^2 |               ...
    perm_zz32(t1.data(), t1.data());     //           0 |        (u+1)*A1^2
    perm_10zz(t3.data(), a.data());      //          A0 |                 0
    mul_fp2_2x4x1w(t2, a, t3);           //       A0*A1 |               ...
    blend_0x33(t0.data(), t2.data(), t0.data()); //       A0*A1 |              A0^2
    blend_0x33(t2.data(), t2.data(), t1.data()); //       A0*A1 |      (u+1)*(A1^2)
    add_fp2_4x2x1w(r, t2, t0);           //     2*A0*A1 | A0^2+(u+1)*(A1^2)
}

void sqr_fp4_2x2x2x1w_v3(fp4_2x2x2x1w& r, const fp4_2x2x2x1w& a)
{
    fp4_2x2x2x1w t0, t1, t2;
    fp4x2_2x2x2x1w tt0, tt1, tt2;

    // a = A1 | A0 at Fp2 layer
    mul_by_u_plus_1_fp2_4x2x1w(t0, a);      //    A1*(u+1) |                                     ...
    perm_1010(t1.data(), a.data());          //          A0 |                                      A0
    perm_1032(t2.data(), a.data());          //          A0 |                                      A1
    blend_0x33(t2.data(), t0.data(), t2.data()); //    A1*(u+1) |                                      A1
    add_fp2_4x2x1w(t2, t1, t2);              // A0+A1*(u+1) |                                   A0+A1
    blend_0x33(t0.data(), t1.data(), t2.data()); //          A0 |                                   A0+A1
    perm_1032(t2.data(), t2.data());         //       A0+A1 |                             A0+A1*(u+1)
    blend_0x33(t2.data(), a.data(), t2.data()); //          A1 |                             A0+A1*(u+1)
    mul_fp2x2_4x2x1w(tt0, t2, t0);           //       A0*A1 |                   (A0+A1)*(A0+A1*(u+1))
    mul_by_u_plus_1_fp2x2_4x2x1w(tt1, tt0);  // A0*A1*(u+1) |                                     ...
    perm_1032_dl(tt2.data(), tt0.data());    //         ... |                                   A0*A1
    sub_fp2x2_4x2x1w(tt2, tt0, tt2);         //         ... |             (A0+A1)*(A0+A1*(u+1))-A0*A1
    blend_0x33_dl(tt2.data(), tt0.data(), tt2.data()); //       A0*A1 |             (A0+A1)*(A0+A1*(u+1))-A0*A1
    perm_1032_dl(tt1.data(), tt1.data());    //         ... |                             A0*A1*(u+1)
    blend_0x33_dl(tt1.data(), tt0.data(), tt1.data()); //       A0*A1 |                             A0*A1*(u+1)
    UFpPack8x1::aassx2(tt2, tt2, tt1);         //     2*A0*A1 | (A0+A1)*(A0+A1*(u+1))-A0*A1-A0*A1*(u+1)
    redc_fp2x2_4x2x1w(r, tt2);               //     2*A0*A1 | (A0+A1)*(A0+A1*(u+1))-A0*A1-A0*A1*(u+1)
}

void sqr_fp4_1x2x2x2w_v1(fp4_1x2x2x2w& r, const fp4_1x2x2x2w& a)
{
    fp4x2_1x2x2x2w tt0, tt1, tt2;
    fp4_1x2x2x2w t0;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);

    // a = A1 | A0 at Fp2 layer
    sqr_fp2x2_2x2x2w(tt0, a);               //        A1^2 |              A0^2
    mul_by_u_plus_1_fp2x2_2x2x2w(tt1, tt0); //  (u+1)*A1^2 |               ...
    perm_var_vl(tt1.data(), tt1.data(), m0); //         ... |        (u+1)*A1^2
    perm_var_hl(t0.data(), a.data(), m0);   //          A0 |                 0
    mul_fp2x2_1x4x2w(tt2, a, t0);           //       A0*A1 |               ...
    blend_0x0F_vl(tt0.data(), tt2.data(), tt0.data()); //       A0*A1 |              A0^2
    blend_0x0F_vl(tt2.data(), tt2.data(), tt1.data()); //       A0*A1 |      (u+1)*(A1^2)
    add_fp2x2_2x2x2w(tt0, tt0, tt2);        //     2*A0*A1 | A0^2+(u+1)*(A1^2)
    redc_fp2x2_2x2x2w(r, tt0);              //     2*A0*A1 | A0^2+(u+1)*(A1^2)
}

void sqr_fp4_1x2x2x2w_v2(fp4_1x2x2x2w& r, const fp4_1x2x2x2w& a)
{
    fp4_1x2x2x2w t0, t1, t2, t3;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);

    // a = A1 | A0 at Fp2 layer
    sqr_fp2_2x2x2w(t0, a);               //        A1^2 |              A0^2
    mul_by_u_plus_1_fp2_2x2x2w(t1, t0);  //  (u+1)*A1^2 |               ...
    perm_var_hl(t1.data(), t1.data(), m0); //         ... |        (u+1)*A1^2
    perm_var_hl(t3.data(), a.data(), m0);  //          A0 |                 0
    mul_fp2_1x4x2w(t2, a, t3);           //       A0*A1 |               ...
    blend_0x0F_hl(t0.data(), t2.data(), t0.data()); //       A0*A1 |              A0^2
    blend_0x0F_hl(t2.data(), t2.data(), t1.data()); //       A0*A1 |      (u+1)*(A1^2)
    add_fp2_2x2x2w(r, t2, t0);           //     2*A0*A1 | A0^2+(u+1)*(A1^2)
}

} // namespace blst

#endif
