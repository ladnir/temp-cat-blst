

#pragma once


#include "blst/common/defines.h"
#include "blst/backend/backend.h"
#include "blst/fields/Fp2.h"



namespace blst
{


    struct Fp6 {
        Fp2 data[3];

        BLST_HD BLST_FORCEINLINE auto& operator[](size_t idx) {
            return data[idx];
        }
        BLST_HD BLST_FORCEINLINE auto& operator[](size_t idx) const
        {
            return data[idx];
        }
    };

    struct Fp12
    {
        Fp6 data[2];

        BLST_HD BLST_FORCEINLINE auto& operator[](size_t idx) {
            return data[idx];
        }

        BLST_HD BLST_FORCEINLINE auto& operator[](size_t idx) const {
            return data[idx];
        }

        BLST_HD static void add(Fp12& ret, const Fp12& a, const Fp12& b);
        
        BLST_HD static void sub(Fp12& ret, const Fp12& a, const Fp12& b);

        BLST_HD static void sqr(Fp12& ret, const Fp12& a);

        BLST_HD static void mul(Fp12& ret, const Fp12& a, const Fp12& b);

        BLST_HD static void exp(Fp12& ret, const Fp12& a, const byte* exp, size_t bits);


        BLST_HD static const Fp12& one();
        BLST_HD static void cyclotomic_sqr(Fp12& ret, const Fp12& a);
        BLST_HD static void conjugate(Fp12& a);
        BLST_HD static void inverse(Fp12& ret, const Fp12& a);

        /* caveat lector! |n| has to be non-zero and not more than 3! */
        BLST_HD static void frobenius_map(Fp12& ret, const Fp12& a, size_t n);

        BLST_HD static void mul_by_xy00z0(Fp12& ret, const Fp12& a, const Fp6& xy00z0);


        BLST_HD static int equal(const Fp12& a, const Fp12& b);

        BLST_HD static int is_one(const Fp12& a);

        BLST_HD static void from_bendian(Fp12& ret, const unsigned char serialized[12 * 48]);

        BLST_HD static void to_bendian(unsigned char serialized[12 * 48], const Fp12& a);


        BLST_HD static int in_group(const Fp12& f);

        BLST_HD static bool_t is_cyclotomic(const Fp12& f);


        BLST_HD static void mul_n_sqr(Fp12& ret, const Fp12& a, size_t n);

        BLST_HD static void raise_to_z_div_by_2(Fp12& ret, const Fp12& a);

        BLST_HD static void raise_to_z(Fp12& a, const Fp12& b);

        BLST_HD static void final_exp(Fp12& ret, const Fp12& f);

    };



//
//
//
//#if defined(BLST_ASM) && !defined(__CUDA_ARCH__)
//#define __FP2x2__
///*
// * Fp2x2 is a "widened" version of Fp2, which allows to consolidate
// * reductions from several multiplications. In other words instead of
// * "mul_redc-mul_redc-add" we get "mul-mul-add-redc," where latter
// * addition is double-width... To be more specific this gives ~7-10%
// * faster pairing depending on platform...
// */
//
//    BLST_HD BLST_INLINE void add_fp2x2(UFp2& ret, const UFp2& a, const UFp2& b)
//    {
//        add_mod_384x384(&ret[0], &a[0], &b[0], &get_BLST_P());
//        add_mod_384x384(&ret[1], &a[1], &b[1], &get_BLST_P());
//    }
//
//    BLST_HD BLST_INLINE void sub_fp2x2(UFp2& ret, const UFp2& a, const UFp2& b)
//    {
//        sub_mod_384x384(&ret[0], &a[0], &b[0], &get_BLST_P());
//        sub_mod_384x384(&ret[1], &a[1], &b[1], &get_BLST_P());
//    }
//
//    BLST_HD BLST_INLINE void mul_by_u_plus_1_fp2x2(UFp2& ret, const UFp2& a)
//    {
//        /* caveat lector! |ret| may not be same as |a| */
//        sub_mod_384x384(&ret[0], &a[0], &a[1], &get_BLST_P());
//        add_mod_384x384(&ret[1], &a[0], &a[1], &get_BLST_P());
//    }
//
//    BLST_HD BLST_INLINE void redc_fp2x2(Fp2& ret, const UFp2& a)
//    {
//        Fp::redc(ret[0], a[0]);
//        Fp::redc(ret[1], a[1]);
//    }
//
//    BLST_HD BLST_INLINE void mul_fp2x2(UFp2& ret, const Fp2& a, const Fp2& b)
//    {
//#if 1
//        mul_v382x(ret, a, b, get_BLST_P());   /* +~6% in Miller loop */
//#else
//        union { Fp x[2]; UFp& x2; } t;
//
//        add_mod_384(t.x[0], a[0], a[1], get_BLST_P());
//        add_mod_384(t.x[1], b[0], b[1], get_BLST_P());
//        mul_384(ret[1], t.x[0], t.x[1]);
//
//        mul_384(ret[0], a[0], b[0]);
//        mul_384(t.x2, a[1], b[1]);
//
//        sub_mod_384x384(ret[1], ret[1], ret[0], get_BLST_P());
//        sub_mod_384x384(ret[1], ret[1], t.x2, get_BLST_P());
//
//        sub_mod_384x384(ret[0], ret[0], t.x2, get_BLST_P());
//#endif
//    }
//
//    BLST_HD BLST_INLINE void sqr_fp2x2(UFp2& ret, const Fp2& a)
//    {
//#if 1
//        sqr_v382x(ret, a, get_BLST_P());      /* +~5% in final exponentiation */
//#else
//        Fp& t0, t1;
//
//        add_mod_384(t0, a[0], a[1], get_BLST_P());
//        sub_mod_384(t1, a[0], a[1], get_BLST_P());
//
//        mul_384(ret[1], a[0], a[1]);
//        add_mod_384x384(ret[1], ret[1], ret[1], get_BLST_P());
//
//        mul_384(ret[0], t0, t1);
//#endif
//    }
//
//    //typedef UFp2 UFp6&[3];
//    struct UFp6 {
//        UFp2 data[3];
//
//        UFp2& operator[](size_t idx) {
//            return data[idx];
//		} 
//        const UFp2& operator[](size_t idx) const {
//            return data[idx];
//        }
//	};
//
//    BLST_HD BLST_INLINE void sub_fp6x2(UFp6& ret, const UFp6& a,
//        const UFp6& b)
//    {
//        sub_fp2x2(ret[0], a[0], b[0]);
//        sub_fp2x2(ret[1], a[1], b[1]);
//        sub_fp2x2(ret[2], a[2], b[2]);
//    }
//
//    BLST_HD BLST_INLINE void mul_fp6x2(UFp6& ret, const Fp6& a, const Fp6& b)
//    {
//        UFp2 t0, t1, t2;
//        Fp2 aa, bb;
//        const Fp2& a0 = a[0];
//        mul_fp2x2(t0, a0, b[0]);
//        mul_fp2x2(t1, a[1], b[1]);
//        mul_fp2x2(t2, a[2], b[2]);
//
//        /* ret[0] = ((a1 + a2)*(b1 + b2) - a1*b1 - a2*b2)*(u+1) + a0*b0
//                  = (a1*b2 + a2*b1)*(u+1) + a0*b0 */
//        Fp2::add(aa, a[1], a[2]);
//        Fp2::add(bb, b[1], b[2]);
//        mul_fp2x2(ret[0], aa, bb);
//        sub_fp2x2(ret[0], ret[0], t1);
//        sub_fp2x2(ret[0], ret[0], t2);
//        mul_by_u_plus_1_fp2x2(ret[1], ret[0]);  /* borrow ret[1] for a moment */
//        add_fp2x2(ret[0], ret[1], t0);
//
//        /* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1 + a2*b2*(u+1)
//                  = a0*b1 + a1*b0 + a2*b2*(u+1) */
//        Fp2::add(aa, a[0], a[1]);
//        Fp2::add(bb, b[0], b[1]);
//        mul_fp2x2(ret[1], aa, bb);
//        sub_fp2x2(ret[1], ret[1], t0);
//        sub_fp2x2(ret[1], ret[1], t1);
//        mul_by_u_plus_1_fp2x2(ret[2], t2);      /* borrow ret[2] for a moment */
//        add_fp2x2(ret[1], ret[1], ret[2]);
//
//        /* ret[2] = (a0 + a2)*(b0 + b2) - a0*b0 - a2*b2 + a1*b1
//                  = a0*b2 + a2*b0 + a1*b1 */
//        Fp2::add(aa, a[0], a[2]);
//        Fp2::add(bb, b[0], b[2]);
//        mul_fp2x2(ret[2], aa, bb);
//        sub_fp2x2(ret[2], ret[2], t0);
//        sub_fp2x2(ret[2], ret[2], t2);
//        add_fp2x2(ret[2], ret[2], t1);
//    }
//
//    BLST_HD BLST_INLINE void redc_fp6x2(Fp6& ret, const UFp6& a)
//    {
//        redc_fp2x2(ret[0], a[0]);
//        redc_fp2x2(ret[1], a[1]);
//        redc_fp2x2(ret[2], a[2]);
//    }
//
//    BLST_HD BLST_INLINE void mul_fp6(Fp6& ret, const Fp6& a, const Fp6& b)
//    {
//        UFp6 r;
//
//        mul_fp6x2(r, a, b);
//        redc_fp6x2(ret, r); /* narrow to normal width */
//    }
//
//    BLST_HD BLST_INLINE void sqr_fp6(Fp6& ret, const Fp6& a)
//    {
//        UFp2 s0, m01, m12, s2, rx;
//
//        sqr_fp2x2(s0, a[0]);
//
//        mul_fp2x2(m01, a[0], a[1]);
//        add_fp2x2(m01, m01, m01);
//
//        mul_fp2x2(m12, a[1], a[2]);
//        add_fp2x2(m12, m12, m12);
//
//        sqr_fp2x2(s2, a[2]);
//
//        /* ret[2] = (a0 + a1 + a2)^2 - a0^2 - a2^2 - 2*(a0*a1) - 2*(a1*a2)
//                  = a1^2 + 2*(a0*a2) */
//        Fp2::add(ret[2], a[2], a[1]);
//        Fp2::add(ret[2], ret[2], a[0]);
//        sqr_fp2x2(rx, ret[2]);
//        sub_fp2x2(rx, rx, s0);
//        sub_fp2x2(rx, rx, s2);
//        sub_fp2x2(rx, rx, m01);
//        sub_fp2x2(rx, rx, m12);
//        redc_fp2x2(ret[2], rx);
//
//        /* ret[0] = a0^2 + 2*(a1*a2)*(u+1) */
//        mul_by_u_plus_1_fp2x2(rx, m12);
//        add_fp2x2(rx, rx, s0);
//        redc_fp2x2(ret[0], rx);
//
//        /* ret[1] = a2^2*(u+1) + 2*(a0*a1) */
//        mul_by_u_plus_1_fp2x2(rx, s2);
//        add_fp2x2(rx, rx, m01);
//        redc_fp2x2(ret[1], rx);
//    }
//
//
//
//#else // __FP2x2__
//
//
//
//
//
//
//    BLST_HD BLST_INLINE void mul_fp6(Fp6& ret, const Fp6& a, const Fp6& b)
//    {
//        Fp2 t0, t1, t2, t3, t4, t5;
//
//        Fp2::mul(t0, a[0], b[0]);
//        Fp2::mul(t1, a[1], b[1]);
//        Fp2::mul(t2, a[2], b[2]);
//
//        /* ret[0] = ((a1 + a2)*(b1 + b2) - a1*b1 - a2*b2)*(u+1) + a0*b0
//                  = (a1*b2 + a2*b1)*(u+1) + a0*b0 */
//        Fp2::add(t4, a[1], a[2]);
//        Fp2::add(t5, b[1], b[2]);
//        Fp2::mul(t3, t4, t5);
//        Fp2::sub(t3, t3, t1);
//        Fp2::sub(t3, t3, t2);
//        Fp2::mul_by_u_plus_1(t3, t3);
//        /* Fp2::add(ret[0], t3, t0); considering possible aliasing... */
//
//        /* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1 + a2*b2*(u+1)
//                  = a0*b1 + a1*b0 + a2*b2*(u+1) */
//        Fp2::add(t4, a[0], a[1]);
//        Fp2::add(t5, b[0], b[1]);
//        Fp2::mul(ret[1], t4, t5);
//        Fp2::sub(ret[1], ret[1], t0);
//        Fp2::sub(ret[1], ret[1], t1);
//        Fp2::mul_by_u_plus_1(t4, t2);
//        Fp2::add(ret[1], ret[1], t4);
//
//        /* ret[2] = (a0 + a2)*(b0 + b2) - a0*b0 - a2*b2 + a1*b1
//                  = a0*b2 + a2*b0 + a1*b1 */
//        Fp2::add(t4, a[0], a[2]);
//        Fp2::add(t5, b[0], b[2]);
//        Fp2::mul(ret[2], t4, t5);
//        Fp2::sub(ret[2], ret[2], t0);
//        Fp2::sub(ret[2], ret[2], t2);
//        Fp2::add(ret[2], ret[2], t1);
//
//        Fp2::add(ret[0], t3, t0);    /* ... moved from above */
//    }
//
//    BLST_HD BLST_INLINE void sqr_fp6(Fp6& ret, const Fp6& a)
//    {
//        Fp2 s0, m01, m12, s2;
//
//        Fp2::sqr(s0, a[0]);
//
//        Fp2::mul(m01, a[0], a[1]);
//        Fp2::add(m01, m01, m01);
//
//        Fp2::mul(m12, a[1], a[2]);
//        Fp2::add(m12, m12, m12);
//
//        Fp2::sqr(s2, a[2]);
//
//        /* ret[2] = (a0 + a1 + a2)^2 - a0^2 - a2^2 - 2*(a0*a1) - 2*(a1*a2)
//                  = a1^2 + 2*(a0*a2) */
//        Fp2::add(ret[2], a[2], a[1]);
//        Fp2::add(ret[2], ret[2], a[0]);
//        Fp2::sqr(ret[2], ret[2]);
//        Fp2::sub(ret[2], ret[2], s0);
//        Fp2::sub(ret[2], ret[2], s2);
//        Fp2::sub(ret[2], ret[2], m01);
//        Fp2::sub(ret[2], ret[2], m12);
//
//        /* ret[0] = a0^2 + 2*(a1*a2)*(u+1) */
//        Fp2::mul_by_u_plus_1(ret[0], m12);
//        Fp2::add(ret[0], ret[0], s0);
//
//        /* ret[1] = a2^2*(u+1) + 2*(a0*a1) */
//        Fp2::mul_by_u_plus_1(ret[1], s2);
//        Fp2::add(ret[1], ret[1], m01);
//    }
//#endif // __FP2x2__
//
//    BLST_HD BLST_INLINE void add_fp6(Fp6& ret, const Fp6& a, const Fp6& b)
//    {
//        Fp2::add(ret[0], a[0], b[0]);
//        Fp2::add(ret[1], a[1], b[1]);
//        Fp2::add(ret[2], a[2], b[2]);
//    }
//
//    BLST_HD BLST_INLINE void sub_fp6(Fp6& ret, const Fp6& a, const Fp6& b)
//    {
//        Fp2::sub(ret[0], a[0], b[0]);
//        Fp2::sub(ret[1], a[1], b[1]);
//        Fp2::sub(ret[2], a[2], b[2]);
//    }
//
//    BLST_HD BLST_INLINE void neg_fp6(Fp6& ret, const Fp6& a)
//    {
//        Fp2::neg(ret[0], a[0]);
//        Fp2::neg(ret[1], a[1]);
//        Fp2::neg(ret[2], a[2]);
//    }
//
//
//
//
//    BLST_HD BLST_INLINE void Fp12::add(Fp12& ret, const Fp12& a, const Fp12& b)
//    {
//		add_fp6(ret[0], a[0], b[0]);
//		add_fp6(ret[1], a[1], b[1]);
//    }
//    BLST_HD BLST_INLINE void Fp12::sub(Fp12& ret, const Fp12& a, const Fp12& b)
//    {
//		sub_fp6(ret[0], a[0], b[0]);
//		sub_fp6(ret[1], a[1], b[1]);
//    }
//
//
////#if 0
////#define mul_by_v_fp6 mul_by_v_fp6
////    BLST_HD BLST_INLINE void mul_by_v_fp6(Fp6& ret, const Fp6& a)
////    {
////        Fp2 t;
////
////        Fp2::mul_by_u_plus_1(t, a[2]);
////        vec_copy(ret[2], a[1], sizeof(a[1]));
////        vec_copy(ret[1], a[0], sizeof(a[0]));
////        vec_copy(ret[0], t, sizeof(t));
////    }
////#endif
//
//    /*
//     * Fp12 extension
//     */
//#if defined(__FP2x2__)
//    BLST_HD BLST_INLINE void Fp12::mul(Fp12& ret, const Fp12& a, const Fp12& b)
//    {
//        UFp6 t0, t1, rx;
//        Fp6 t2;
//
//        mul_fp6x2(t0, a[0], b[0]);
//        mul_fp6x2(t1, a[1], b[1]);
//
//        /* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1
//                  = a0*b1 + a1*b0 */
//        add_fp6(t2, a[0], a[1]);
//        add_fp6(ret[1], b[0], b[1]);
//        mul_fp6x2(rx, ret[1], t2);
//        sub_fp6x2(rx, rx, t0);
//        sub_fp6x2(rx, rx, t1);
//        redc_fp6x2(ret[1], rx);
//
//        /* ret[0] = a0*b0 + a1*b1*v */
//        mul_by_u_plus_1_fp2x2(rx[0], t1[2]);
//        add_fp2x2(rx[0], t0[0], rx[0]);
//        add_fp2x2(rx[1], t0[1], t1[0]);
//        add_fp2x2(rx[2], t0[2], t1[1]);
//        redc_fp6x2(ret[0], rx);
//    }
//
//    BLST_HD BLST_INLINE void mul_by_0y0_fp6x2(UFp6& ret, const Fp6& a,
//        const Fp2& b)
//    {
//        mul_fp2x2(ret[1], a[2], b);     /* borrow ret[1] for a moment */
//        mul_by_u_plus_1_fp2x2(ret[0], ret[1]);
//        mul_fp2x2(ret[1], a[0], b);
//        mul_fp2x2(ret[2], a[1], b);
//    }
//
//    BLST_HD BLST_INLINE void mul_by_xy0_fp6x2(UFp6& ret, const Fp6& a,
//        const Fp6& b)
//    {
//        UFp2 t0, t1;
//        Fp2 aa, bb;
//
//        mul_fp2x2(t0, a[0], b[0]);
//        mul_fp2x2(t1, a[1], b[1]);
//
//        /* ret[0] = ((a1 + a2)*(b1 + 0) - a1*b1 - a2*0)*(u+1) + a0*b0
//                  = (a1*0 + a2*b1)*(u+1) + a0*b0 */
//        mul_fp2x2(ret[1], a[2], b[1]);  /* borrow ret[1] for a moment */
//        mul_by_u_plus_1_fp2x2(ret[0], ret[1]);
//        add_fp2x2(ret[0], ret[0], t0);
//
//        /* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1 + a2*0*(u+1)
//                  = a0*b1 + a1*b0 + a2*0*(u+1) */
//        Fp2::add(aa, a[0], a[1]);
//        Fp2::add(bb, b[0], b[1]);
//        mul_fp2x2(ret[1], aa, bb);
//        sub_fp2x2(ret[1], ret[1], t0);
//        sub_fp2x2(ret[1], ret[1], t1);
//
//        /* ret[2] = (a0 + a2)*(b0 + 0) - a0*b0 - a2*0 + a1*b1
//                  = a0*0 + a2*b0 + a1*b1 */
//        mul_fp2x2(ret[2], a[2], b[0]);
//        add_fp2x2(ret[2], ret[2], t1);
//    }
//
//    BLST_HD BLST_INLINE void Fp12::mul_by_xy00z0(Fp12& ret, const Fp12& a,
//        const Fp6& xy00z0)
//    {
//        UFp6 t0, t1, rr;
//        Fp6 t2;
//
//        mul_by_xy0_fp6x2(t0, a[0], xy00z0);
//        mul_by_0y0_fp6x2(t1, a[1], xy00z0[2]);
//
//        /* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1
//                  = a0*b1 + a1*b0 */
//        vec_copy(&t2[0], &xy00z0[0], sizeof(t2[0]));
//        Fp2::add(t2[1], xy00z0[1], xy00z0[2]);
//        add_fp6(ret[1], a[0], a[1]);
//        mul_by_xy0_fp6x2(rr, ret[1], t2);
//        sub_fp6x2(rr, rr, t0);
//        sub_fp6x2(rr, rr, t1);
//        redc_fp6x2(ret[1], rr);
//
//        /* ret[0] = a0*b0 + a1*b1*v */
//        mul_by_u_plus_1_fp2x2(rr[0], t1[2]);
//        add_fp2x2(rr[0], t0[0], rr[0]);
//        add_fp2x2(rr[1], t0[1], t1[0]);
//        add_fp2x2(rr[2], t0[2], t1[1]);
//        redc_fp6x2(ret[0], rr);
//    }
//#else // __FP2x2__
//    BLST_HD BLST_INLINE void Fp12::mul(Fp12& ret, const Fp12& a, const Fp12& b)
//    {
//        Fp6 t0, t1, t2;
//
//        mul_fp6(t0, a[0], b[0]);
//        mul_fp6(t1, a[1], b[1]);
//
//        /* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1
//                  = a0*b1 + a1*b0 */
//        add_fp6(t2, a[0], a[1]);
//        add_fp6(ret[1], b[0], b[1]);
//        mul_fp6(ret[1], ret[1], t2);
//        sub_fp6(ret[1], ret[1], t0);
//        sub_fp6(ret[1], ret[1], t1);
//
//        /* ret[0] = a0*b0 + a1*b1*v */
//#ifdef mul_by_v_fp6
//        mul_by_v_fp6(t1, t1);
//        add_fp6(ret[0], t0, t1);
//#else
//        Fp2::mul_by_u_plus_1(t1[2], t1[2]);
//        Fp2::add(ret[0][0], t0[0], t1[2]);
//        Fp2::add(ret[0][1], t0[1], t1[0]);
//        Fp2::add(ret[0][2], t0[2], t1[1]);
//#endif
//    }
//
//    BLST_HD BLST_INLINE void mul_by_0y0_fp6(Fp6& ret, const Fp6& a,
//        const Fp2& b)
//    {
//        Fp2 t;
//
//        Fp2::mul(t, a[2], b);
//        Fp2::mul(ret[2], a[1], b);
//        Fp2::mul(ret[1], a[0], b);
//        Fp2::mul_by_u_plus_1(ret[0], t);
//    }
//
//    BLST_HD BLST_INLINE void mul_by_xy0_fp6(Fp6& ret, const Fp6& a, const Fp6& b)
//    {
//        Fp2 t0, t1, /*t2,*/ t3, t4, t5;
//
//        Fp2::mul(t0, a[0], b[0]);
//        Fp2::mul(t1, a[1], b[1]);
//
//        /* ret[0] = ((a1 + a2)*(b1 + 0) - a1*b1 - a2*0)*(u+1) + a0*b0
//                  = (a1*0 + a2*b1)*(u+1) + a0*b0 */
//        Fp2::mul(t3, a[2], b[1]);
//        Fp2::mul_by_u_plus_1(t3, t3);
//        /* Fp2::add(ret[0], t3, t0); considering possible aliasing... */
//
//        /* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1 + a2*0*(u+1)
//                  = a0*b1 + a1*b0 + a2*0*(u+1) */
//        Fp2::add(t4, a[0], a[1]);
//        Fp2::add(t5, b[0], b[1]);
//        Fp2::mul(ret[1], t4, t5);
//        Fp2::sub(ret[1], ret[1], t0);
//        Fp2::sub(ret[1], ret[1], t1);
//
//        /* ret[2] = (a0 + a2)*(b0 + 0) - a0*b0 - a2*0 + a1*b1
//                  = a0*0 + a2*b0 + a1*b1 */
//        Fp2::mul(ret[2], a[2], b[0]);
//        Fp2::add(ret[2], ret[2], t1);
//
//        Fp2::add(ret[0], t3, t0);    /* ... moved from above */
//    }
//
//    BLST_HD BLST_INLINE void Fp12::mul_by_xy00z0(Fp12& ret, const Fp12& a,
//        const Fp6& xy00z0)
//    {
//        Fp6 t0, t1, t2;
//
//        mul_by_xy0_fp6(t0, a[0], xy00z0);
//        mul_by_0y0_fp6(t1, a[1], xy00z0[2]);
//
//        /* ret[1] = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1
//                  = a0*b1 + a1*b0 */
//        vec_copy(&t2[0], &xy00z0[0], sizeof(t2[0]));
//        Fp2::add(t2[1], xy00z0[1], xy00z0[2]);
//        add_fp6(ret[1], a[0], a[1]);
//        mul_by_xy0_fp6(ret[1], ret[1], t2);
//        sub_fp6(ret[1], ret[1], t0);
//        sub_fp6(ret[1], ret[1], t1);
//
//        /* ret[0] = a0*b0 + a1*b1*v */
//#ifdef mul_by_v_fp6
//        mul_by_v_fp6(t1, t1);
//        add_fp6(ret[0], t0, t1);
//#else
//        Fp2::mul_by_u_plus_1(t1[2], t1[2]);
//        Fp2::add(ret[0][0], t0[0], t1[2]);
//        Fp2::add(ret[0][1], t0[1], t1[0]);
//        Fp2::add(ret[0][2], t0[2], t1[1]);
//#endif
//    }
//#endif // __FP2x2__
//
//    BLST_HD BLST_INLINE void Fp12::sqr(Fp12& ret, const Fp12& a)
//    {
//        Fp6 t0, t1;
//
//        add_fp6(t0, a[0], a[1]);
//#ifdef mul_by_v_fp6
//        mul_by_v_fp6(t1, a[1]);
//        add_fp6(t1, a[0], t1);
//#else
//        Fp2::mul_by_u_plus_1(t1[2], a[1][2]);
//        Fp2::add(t1[0], a[0][0], t1[2]);
//        Fp2::add(t1[1], a[0][1], a[1][0]);
//        Fp2::add(t1[2], a[0][2], a[1][1]);
//#endif
//        mul_fp6(t0, t0, t1);
//        mul_fp6(t1, a[0], a[1]);
//
//        /* ret[1] = 2*(a0*a1) */
//        add_fp6(ret[1], t1, t1);
//
//        /* ret[0] = (a0 + a1)*(a0 + a1*v) - a0*a1 - a0*a1*v
//                  = a0^2 + a1^2*v */
//        sub_fp6(ret[0], t0, t1);
//#ifdef mul_by_v_fp6
//        mul_by_v_fp6(t1, t1);
//        sub_fp6(ret[0], ret[0], t1);
//#else
//        Fp2::mul_by_u_plus_1(t1[2], t1[2]);
//        Fp2::sub(ret[0][0], ret[0][0], t1[2]);
//        Fp2::sub(ret[0][1], ret[0][1], t1[0]);
//        Fp2::sub(ret[0][2], ret[0][2], t1[1]);
//#endif
//    }
//
//    BLST_HD BLST_INLINE void Fp12::conjugate(Fp12& a)
//    {
//        neg_fp6(a[1], a[1]);
//    }
//
//    BLST_HD BLST_INLINE void inverse_fp6(Fp6& ret, const Fp6& a)
//    {
//        Fp2 c0, c1, c2, t0, t1;
//
//        /* c0 = a0^2 - (a1*a2)*(u+1) */
//        Fp2::sqr(c0, a[0]);
//        Fp2::mul(t0, a[1], a[2]);
//        Fp2::mul_by_u_plus_1(t0, t0);
//        Fp2::sub(c0, c0, t0);
//
//        /* c1 = a2^2*(u+1) - (a0*a1) */
//        Fp2::sqr(c1, a[2]);
//        Fp2::mul_by_u_plus_1(c1, c1);
//        Fp2::mul(t0, a[0], a[1]);
//        Fp2::sub(c1, c1, t0);
//
//        /* c2 = a1^2 - a0*a2 */
//        Fp2::sqr(c2, a[1]);
//        Fp2::mul(t0, a[0], a[2]);
//        Fp2::sub(c2, c2, t0);
//
//        /* (a2*c1 + a1*c2)*(u+1) + a0*c0 */
//        Fp2::mul(t0, c1, a[2]);
//        Fp2::mul(t1, c2, a[1]);
//        Fp2::add(t0, t0, t1);
//        Fp2::mul_by_u_plus_1(t0, t0);
//        Fp2::mul(t1, c0, a[0]);
//        Fp2::add(t0, t0, t1);
//
//        Fp2::reciprocal(t1, t0);
//
//        Fp2::mul(ret[0], c0, t1);
//        Fp2::mul(ret[1], c1, t1);
//        Fp2::mul(ret[2], c2, t1);
//    }
//
//    BLST_HD BLST_INLINE void Fp12::inverse(Fp12& ret, const Fp12& a)
//    {
//        Fp6 t0, t1;
//
//        sqr_fp6(t0, a[0]);
//        sqr_fp6(t1, a[1]);
//#ifdef mul_by_v_fp6
//        mul_by_v_fp6(t1, t1);
//        sub_fp6(t0, t0, t1);
//#else
//        Fp2::mul_by_u_plus_1(t1[2], t1[2]);
//        Fp2::sub(t0[0], t0[0], t1[2]);
//        Fp2::sub(t0[1], t0[1], t1[0]);
//        Fp2::sub(t0[2], t0[2], t1[1]);
//#endif
//
//        inverse_fp6(t1, t0);
//
//        mul_fp6(ret[0], a[0], t1);
//        mul_fp6(ret[1], a[1], t1);
//        neg_fp6(ret[1], ret[1]);
//    }
//
//    //typedef Fp2 Fpfp4[2];
//    struct Fp4 {
//        Fp2 data[2];
//        BLST_HD Fp2& operator[](size_t idx) {
//            return data[idx];
//        }
//        BLST_HD const Fp2& operator[](size_t idx)const {
//            return data[idx];
//        }
//    };
//
//#if defined(__FP2x2__)
//    BLST_HD BLST_INLINE void sqr_fp4(Fp4& ret, const Fp2& a0, const Fp2& a1)
//    {
//        UFp2 t0, t1, t2;
//
//        sqr_fp2x2(t0, a0);
//        sqr_fp2x2(t1, a1);
//        Fp2::add(ret[1], a0, a1);
//
//        mul_by_u_plus_1_fp2x2(t2, t1);
//        add_fp2x2(t2, t2, t0);
//        redc_fp2x2(ret[0], t2);
//
//        sqr_fp2x2(t2, ret[1]);
//        sub_fp2x2(t2, t2, t0);
//        sub_fp2x2(t2, t2, t1);
//        redc_fp2x2(ret[1], t2);
//    }
//#else //__FP2x2__
//    BLST_HD BLST_INLINE void sqr_fp4(Fp4& ret, const Fp2& a0, const Fp2& a1)
//    {
//        Fp2 t0, t1;
//
//        Fp2::sqr(t0, a0);
//        Fp2::sqr(t1, a1);
//        Fp2::add(ret[1], a0, a1);
//
//        Fp2::mul_by_u_plus_1(ret[0], t1);
//        Fp2::add(ret[0], ret[0], t0);
//
//        Fp2::sqr(ret[1], ret[1]);
//        Fp2::sub(ret[1], ret[1], t0);
//        Fp2::sub(ret[1], ret[1], t1);
//    }
//#endif // __FP2x2__
//
//    BLST_HD BLST_INLINE void Fp12::cyclotomic_sqr(Fp12& ret, const Fp12& a)
//    {
//        Fp4 t0, t1, t2;
//
//        sqr_fp4(t0, a[0][0], a[1][1]);
//        sqr_fp4(t1, a[1][0], a[0][2]);
//        sqr_fp4(t2, a[0][1], a[1][2]);
//
//        Fp2::sub(ret[0][0], t0[0], a[0][0]);
//        Fp2::add(ret[0][0], ret[0][0], ret[0][0]);
//        Fp2::add(ret[0][0], ret[0][0], t0[0]);
//
//        Fp2::sub(ret[0][1], t1[0], a[0][1]);
//        Fp2::add(ret[0][1], ret[0][1], ret[0][1]);
//        Fp2::add(ret[0][1], ret[0][1], t1[0]);
//
//        Fp2::sub(ret[0][2], t2[0], a[0][2]);
//        Fp2::add(ret[0][2], ret[0][2], ret[0][2]);
//        Fp2::add(ret[0][2], ret[0][2], t2[0]);
//
//        Fp2::mul_by_u_plus_1(t2[1], t2[1]);
//        Fp2::add(ret[1][0], t2[1], a[1][0]);
//        Fp2::add(ret[1][0], ret[1][0], ret[1][0]);
//        Fp2::add(ret[1][0], ret[1][0], t2[1]);
//
//        Fp2::add(ret[1][1], t0[1], a[1][1]);
//        Fp2::add(ret[1][1], ret[1][1], ret[1][1]);
//        Fp2::add(ret[1][1], ret[1][1], t0[1]);
//
//        Fp2::add(ret[1][2], t1[1], a[1][2]);
//        Fp2::add(ret[1][2], ret[1][2], ret[1][2]);
//        Fp2::add(ret[1][2], ret[1][2], t1[1]);
//    }
//
//    BLST_HD BLST_INLINE void frobenius_map_fp6(Fp6& ret, const Fp6& a, size_t n)
//    {
//        const Fp2 coeffs1[] = {  /* (u + 1)^((P^n - 1) / 3) */
//          Fp2{ Fp{ 0 },
//            Fp{ TO_LIMB_T(0xcd03c9e48671f071), TO_LIMB_T(0x5dab22461fcda5d2),
//              TO_LIMB_T(0x587042afd3851b95), TO_LIMB_T(0x8eb60ebe01bacb9e),
//              TO_LIMB_T(0x03f97d6e83d050d2), TO_LIMB_T(0x18f0206554638741) } },
//          Fp2{ 
//            Fp{ TO_LIMB_T(0x30f1361b798a64e8), TO_LIMB_T(0xf3b8ddab7ece5a2a),
//              TO_LIMB_T(0x16a8ca3ac61577f7), TO_LIMB_T(0xc26a2ff874fd029b),
//              TO_LIMB_T(0x3636b76660701c6e), TO_LIMB_T(0x051ba4ab241b6160) } },
//          Fp2{ Fp{ 0 }, Fp{ ONE_MONT_P_LIMBS } }
//        };
//        const Fp coeffs2[] = {  /* (u + 1)^((2P^n - 2) / 3) */
//          {   TO_LIMB_T(0x890dc9e4867545c3), TO_LIMB_T(0x2af322533285a5d5),
//              TO_LIMB_T(0x50880866309b7e2c), TO_LIMB_T(0xa20d1b8c7e881024),
//              TO_LIMB_T(0x14e4f04fe2db9068), TO_LIMB_T(0x14e56d3f1564853a)   },
//          {   TO_LIMB_T(0xcd03c9e48671f071), TO_LIMB_T(0x5dab22461fcda5d2),
//              TO_LIMB_T(0x587042afd3851b95), TO_LIMB_T(0x8eb60ebe01bacb9e),
//              TO_LIMB_T(0x03f97d6e83d050d2), TO_LIMB_T(0x18f0206554638741)   },
//          {   TO_LIMB_T(0x43f5fffffffcaaae), TO_LIMB_T(0x32b7fff2ed47fffd),
//              TO_LIMB_T(0x07e83a49a2e99d69), TO_LIMB_T(0xeca8f3318332bb7a),
//              TO_LIMB_T(0xef148d1ea0f4c069), TO_LIMB_T(0x040ab3263eff0206)   }
//        };
//
//        frobenius_map_fp2(ret[0], a[0], n);
//        frobenius_map_fp2(ret[1], a[1], n);
//        frobenius_map_fp2(ret[2], a[2], n);
//        --n;    /* implied ONE_MONT_P at index 0 */
//        Fp2::mul(ret[1], ret[1], coeffs1[n]);
//        Fp::mul(ret[2][0], ret[2][0], coeffs2[n]);
//        Fp::mul(ret[2][1], ret[2][1], coeffs2[n]);
//    }
//
//    BLST_HD BLST_INLINE void Fp12::frobenius_map(Fp12& ret, const Fp12& a, size_t n)
//    {
//        const Fp2 coeffs[] = {  /* (u + 1)^((P^n - 1) / 6) */
//          Fp2{ 
//            Fp{ TO_LIMB_T(0x07089552b319d465), TO_LIMB_T(0xc6695f92b50a8313),
//              TO_LIMB_T(0x97e83cccd117228f), TO_LIMB_T(0xa35baecab2dc29ee),
//              TO_LIMB_T(0x1ce393ea5daace4d), TO_LIMB_T(0x08f2220fb0fb66eb) },
//            Fp{ TO_LIMB_T(0xb2f66aad4ce5d646), TO_LIMB_T(0x5842a06bfc497cec),
//              TO_LIMB_T(0xcf4895d42599d394), TO_LIMB_T(0xc11b9cba40a8e8d0),
//              TO_LIMB_T(0x2e3813cbe5a0de89), TO_LIMB_T(0x110eefda88847faf) } },
//          Fp2{ 
//            Fp{ TO_LIMB_T(0xecfb361b798dba3a), TO_LIMB_T(0xc100ddb891865a2c),
//              TO_LIMB_T(0x0ec08ff1232bda8e), TO_LIMB_T(0xd5c13cc6f1ca4721),
//              TO_LIMB_T(0x47222a47bf7b5c04), TO_LIMB_T(0x0110f184e51c5f59) } },
//          Fp2{ 
//            Fp{ TO_LIMB_T(0x3e2f585da55c9ad1), TO_LIMB_T(0x4294213d86c18183),
//              TO_LIMB_T(0x382844c88b623732), TO_LIMB_T(0x92ad2afd19103e18),
//              TO_LIMB_T(0x1d794e4fac7cf0b9), TO_LIMB_T(0x0bd592fc7d825ec8) },
//            Fp{ TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),
//              TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),
//              TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2) } },
//        };
//
//        frobenius_map_fp6(ret[0], a[0], n);
//        frobenius_map_fp6(ret[1], a[1], n);
//        --n;    /* implied ONE_MONT_P at index 0 */
//        Fp2::mul(ret[1][0], ret[1][0], coeffs[n]);
//        Fp2::mul(ret[1][1], ret[1][1], coeffs[n]);
//        Fp2::mul(ret[1][2], ret[1][2], coeffs[n]);
//    }
//
//
//    /*
//     * BLS12-381-specific Fp12 shortcuts.
//     */
//    //BLST_HD BLST_INLINE void blst_fp12_sqr(Fp12& ret, const Fp12& a)
//    //{
//    //    sqr_fp12(ret, a);
//    //}
//
//    //BLST_HD BLST_INLINE void blst_fp12_cyclotomic_sqr(Fp12& ret, const Fp12& a)
//    //{
//    //    cyclotomic_sqr_fp12(ret, a);
//    //}
//
//    //BLST_HD BLST_INLINE void blst_fp12_mul(Fp12& ret, const Fp12& a, const Fp12& b)
//    //{
//    //    mul_fp12(ret, a, b);
//    //}
//
//    //BLST_HD BLST_INLINE void blst_fp12_mul_by_xy00z0(Fp12& ret, const Fp12& a,
//    //    const Fp6& xy00z0)
//    //{
//    //    mul_by_xy00z0_fp12(ret, a, xy00z0);
//    //}
//
//    //BLST_HD BLST_INLINE void blst_fp12_conjugate(Fp12& a)
//    //{
//    //    conjugate_fp12(a);
//    //}
//
//    //BLST_HD BLST_INLINE void blst_fp12_inverse(Fp12& ret, const Fp12& a)
//    //{
//    //    inverse_fp12(ret, a);
//    //}
//
//    /* caveat lector! |n| has to be non-zero and not more than 3! */
//    //BLST_HD BLST_INLINE void blst_fp12_frobenius_map(Fp12& ret, const Fp12& a, size_t n)
//    //{
//    //    frobenius_map_fp12(ret, a, n);
//    //}
//
//    BLST_HD BLST_INLINE int Fp12::equal(const Fp12& a, const Fp12& b)
//    {
//        return (int)vec_is_equal(&a, &b, sizeof(Fp12));
//    }
//
//    BLST_HD BLST_INLINE int Fp12::is_one(const Fp12& a)
//    {
//        return 
//            (int)(vec_is_equal(
//                &a[0][0], 
//                &get_BLST_Rx().p2, sizeof(a[0][0])) &
//            vec_is_zero(&a[0][1], sizeof(Fp12) - sizeof(a[0][0])));
//    }
//
//
//    BLST_HD BLST_INLINE const Fp12& Fp12::one()
//    {
//        return *(Fp12*)get_BLST_Rx().p12;
//    }
//
//
//    BLST_HD BLST_INLINE void blst_fp12_one_copy(Fp12& ret)
//    {
//        char* src = (char*)get_BLST_Rx().p12;
//        char* dst = (char*)&ret;
//        int size = sizeof(Fp12);
//        for (int i = 0; i < size; ++i)
//            dst[i] = src[i];
//    }
//
//    BLST_HD BLST_INLINE void Fp12::from_bendian(Fp12& ret, const unsigned char serialized[12 * 48]) {
//        // Clear the output first
//        for (int i = 0; i < sizeof(Fp12); ++i)
//            ((unsigned char*)&ret)[i] = 0;
//
//        // Fp12& is structured as ret[2][3][2] representing:
//        // - 2 components of Fp6 (for Fp12 = Fp6[w])
//        // - 3 components of Fp2 each (for Fp6 = Fp2[v])  
//        // - 2 components of Fp each (for Fp2 = Fp[u])
//
//        size_t byte_offset = 0;
//        for (size_t fp6_comp = 0; fp6_comp < 3; fp6_comp++) {      // Fp6 components  
//            for (size_t fp12_comp = 0; fp12_comp < 2; fp12_comp++) {       // Fp12 components
//                for (size_t fp2_comp = 0; fp2_comp < 2; fp2_comp++) {  // Fp2 components
//                    Fp r;
//                    // Convert from big-endian bytes to Fp (limbs)
//                    limbs_from_be_bytes(r, &serialized[byte_offset], 48);
//                    Fp::to_mont(ret[fp12_comp][fp6_comp][fp2_comp], r);
//                    byte_offset += 48;
//                }
//            }
//        }
//    }
//
//    BLST_HD BLST_INLINE void Fp12::to_bendian(unsigned char ret[48 * 12], const Fp12& a)
//    {
//        size_t i, j;
//        Fp out;
//
//        for (i = 0; i < 3; i++) {
//            for (j = 0; j < 2; j++) {
//                for (int k = 0; k < 2; ++k)
//                {
//                    Fp::from_mont(out, a[j][i][k]);
//                    be_bytes_from_limbs(ret, out, sizeof(Fp));
//                    ret += 48;
//                }
//            }
//        }
//    }
//
//    BLST_HD BLST_INLINE size_t blst_fp12_sizeof(void)
//    {
//        return sizeof(Fp12);
//    }
//
//
//    BLST_HD BLST_INLINE void Fp12::exp(Fp12& ret, const Fp12& a, const byte* exp, size_t bits)
//    {
//        // For bits=0 or if exponent is 0, return 1
//        bool_t is_zero = 1;
//        for (size_t i = 0; i < (bits + 8 * sizeof(byte) - 1) / (8 * sizeof(byte)); i++) {
//            if (exp[i] != 0) {
//                is_zero = 0;
//                break;
//            }
//        }
//
//        if (bits == 0 || is_zero) {
//            ret = Fp12::one();
//            return;
//        }
//
//        // Make a copy of the base to prevent issues if a == ret
//        Fp12 base;
//        //memcpy(base, a, sizeof(Fp12));
//        for (int i = 0; i < sizeof(Fp12); ++i)
//            ((unsigned char*)&base)[i] = ((unsigned char*)&a)[i];
//
//        // Start with 1
//        ret = Fp12::one();
//
//        // Process bits from least to most significant
//        for (size_t i = 0; i < bits; i++) {
//            // If the current bit is set
//            size_t limb_idx = i / (8 * sizeof(byte));
//            size_t bit_idx = i % (8 * sizeof(byte));
//
//            if ((exp[limb_idx] >> bit_idx) & 1) {
//                Fp12::mul(ret, ret, base);
//            }
//
//            // Square the base for the next bit position (if not the last bit)
//            if (i < bits - 1) {
//                Fp12::sqr(base, base);
//            }
//
//
//        }
//    }
//
//
//
//
//    BLST_HD BLST_INLINE bool_t Fp12::is_cyclotomic(const Fp12& f)
//    {
//        Fp12 a, b;
//
//        Fp12::frobenius_map(a, f, 2);
//        Fp12::frobenius_map(b, a, 2);
//        Fp12::mul(b, b, f);
//
//        return vec_is_equal(&a, &b, sizeof(a));
//    }
//
//
//    BLST_HD BLST_INLINE void Fp12::mul_n_sqr(Fp12& ret, const Fp12& a, size_t n)
//    {
//        Fp12::mul(ret, ret, a);
//        while (n--)
//            Fp12::cyclotomic_sqr(ret, ret);
//    }
//
//    BLST_HD BLST_INLINE void Fp12::raise_to_z_div_by_2(Fp12& ret, const Fp12& a)
//    {
//        Fp12::cyclotomic_sqr(ret, a);                /* 0x2                  */
//        Fp12::mul_n_sqr(ret, a, 2);                       /* ..0xc                */
//        Fp12::mul_n_sqr(ret, a, 3);                       /* ..0x68               */
//        Fp12::mul_n_sqr(ret, a, 9);                       /* ..0xd200             */
//        Fp12::mul_n_sqr(ret, a, 32);                      /* ..0xd20100000000     */
//        Fp12::mul_n_sqr(ret, a, 16 - 1);                    /* ..0x6900800000008000 */
//        Fp12::conjugate(ret);                /* account for z being negative */
//    }
//
//    BLST_HD BLST_INLINE void Fp12::raise_to_z(Fp12& a, const Fp12& b)
//    {
//        raise_to_z_div_by_2(a, b);
//        Fp12::cyclotomic_sqr(a, a);
//    }
//
//    BLST_HD BLST_INLINE int Fp12::in_group(const Fp12& f)
//    {
//        Fp12 a, b;
//
//        if (vec_is_zero(&f, sizeof(Fp12)) || !is_cyclotomic(f))
//            return 0;
//
//        Fp12::frobenius_map(a, f, 1);
//        raise_to_z(b, f);
//
//        return (int)vec_is_equal(&a, &b, sizeof(a));
//    }
//
//
//
//    /*
//     * Adaptation from <zkcrypto>/pairing/src/bls12_381/mod.rs
//     */
//    BLST_HD BLST_INLINE void Fp12::final_exp(Fp12& ret, const Fp12& f)
//    {
//        Fp12 y0, y1, y2, y3;
//
//        vec_copy(&y1, &f, sizeof(y1));
//        Fp12::conjugate(y1);
//        Fp12::inverse(y2, f);
//        Fp12::mul(ret, y1, y2);
//        Fp12::frobenius_map(y2, ret, 2);
//        Fp12::mul(ret, ret, y2);
//
//        Fp12::cyclotomic_sqr(y0, ret);
//        Fp12::raise_to_z(y1, y0);
//        Fp12::raise_to_z_div_by_2(y2, y1);
//        vec_copy(&y3, &ret, sizeof(y3));
//        Fp12::conjugate(y3);
//        Fp12::mul(y1, y1, y3);
//        Fp12::conjugate(y1);
//        Fp12::mul(y1, y1, y2);
//        Fp12::raise_to_z(y2, y1);
//        Fp12::raise_to_z(y3, y2);
//        Fp12::conjugate(y1);
//        Fp12::mul(y3, y3, y1);
//        Fp12::conjugate(y1);
//        Fp12::frobenius_map(y1, y1, 3);
//        Fp12::frobenius_map(y2, y2, 2);
//        Fp12::mul(y1, y1, y2);
//        Fp12::raise_to_z(y2, y3);
//        Fp12::mul(y2, y2, y0);
//        Fp12::mul(y2, y2, ret);
//        Fp12::mul(y1, y1, y2);
//        Fp12::frobenius_map(y2, y3, 1);
//        Fp12::mul(ret, y1, y2);
//    }
}



