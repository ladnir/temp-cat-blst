/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/curve/point.h"

#ifndef MILLER_LOOP_N_MAX
# define MILLER_LOOP_N_MAX 16
#endif

namespace blst
{

    /*
    * Line evaluations from  https://eprint.iacr.org/2010/354.pdf
    * with a twist moving common expression to line_by_Px2.
    */
    BLST_HD void line_add(Fp6& line, POINTonE2& T, const POINTonE2& R,
        const POINTonE2_affine& Q)
    {
        Fp2 Z1Z1, U2, S2, H, HH, I, J, V;
#if 1
# define r line[1]
#else
        Fp2 r;
#endif

        /*
         * https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
         * with XYZ3 being |T|, XYZ1 - |R|, XY2 - |Q|, i.e. Q is affine
         */
        Fp2::sqr(Z1Z1, R.Z);                /* Z1Z1 = Z1^2 */
        Fp2::mul(U2, Q.X, Z1Z1);            /* U2 = X2*Z1Z1 */

        Fp2::mul(S2, Q.Y, R.Z);
        Fp2::mul(S2, S2, Z1Z1);              /* S2 = Y2*Z1*Z1Z1 */

        Fp2::sub(H, U2, R.X);               /* H = U2-X1 */

        Fp2::sqr(HH, H);                     /* HH = H^2 */
        Fp2::add(I, HH, HH);
        Fp2::add(I, I, I);                   /* I = 4*HH */

        Fp2::mul(J, H, I);                   /* J = H*I */

        Fp2::sub(r, S2, R.Y);
        Fp2::add(r, r, r);                   /* r = 2*(S2-Y1) */

        Fp2::mul(V, R.X, I);                /* V = X1*I */

        Fp2::sqr(T.X, r);
        Fp2::sub(T.X, T.X, J);
        Fp2::sub(T.X, T.X, V);
        Fp2::sub(T.X, T.X, V);             /* X3 = r^2-J-2*V */

        Fp2::mul(J, J, R.Y);
        Fp2::sub(T.Y, V, T.X);
        Fp2::mul(T.Y, T.Y, r);
        Fp2::sub(T.Y, T.Y, J);
        Fp2::sub(T.Y, T.Y, J);             /* Y3 = r*(V-X3)-2*Y1*J */

        Fp2::add(T.Z, R.Z, H);
        Fp2::sqr(T.Z, T.Z);
        Fp2::sub(T.Z, T.Z, Z1Z1);
        Fp2::sub(T.Z, T.Z, HH);            /* Z3 = (Z1+H)^2-Z1Z1-HH */

        /*
         * line evaluation
         */
        Fp2::mul(I, r, Q.X);
        Fp2::mul(J, Q.Y, T.Z);
        Fp2::sub(I, I, J);
        Fp2::add(line[0], I, I);          /* 2*(r*X2 - Y2*Z3) */
#ifdef r
# undef r
#else
        vec_copy(&line[1], &r, sizeof(r));
#endif
        vec_copy(&line[2], &T.Z, sizeof(T.Z));
    }

    BLST_HD void line_dbl(Fp6& line, POINTonE2& T, const POINTonE2& Q)
    {
        Fp2 ZZ, A, B, C, D, E, F;

        /*
         * https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-alnr
         */
        Fp2::sqr(A, Q.X);                   /* A = X1^2 */
        Fp2::sqr(B, Q.Y);                   /* B = Y1^2 */
        Fp2::sqr(ZZ, Q.Z);                  /* ZZ = Z1^2 */
        Fp2::sqr(C, B);                      /* C = B^2 */

        Fp2::add(D, Q.X, B);                /* X1+B */
        Fp2::sqr(D, D);                      /* (X1+B)^2 */
        Fp2::sub(D, D, A);                   /* (X1+B)^2-A */
        Fp2::sub(D, D, C);                   /* (X1+B)^2-A-C */
        Fp2::add(D, D, D);                   /* D = 2*((X1+B)^2-A-C) */

        mul_by_3_fp2(E, A);                 /* E = 3*A */
        Fp2::sqr(F, E);                      /* F = E^2 */

        Fp2::add(line[0], E, Q.X);          /* 3*A+X1 for line evaluation */

        Fp2::sub(T.X, F, D);
        Fp2::sub(T.X, T.X, D);             /* X3 = F-2*D */

        Fp2::add(T.Z, Q.Y, Q.Z);
        Fp2::sqr(T.Z, T.Z);
        Fp2::sub(T.Z, T.Z, B);
        Fp2::sub(T.Z, T.Z, ZZ);            /* Z3 = (Y1+Z1)^2-B-ZZ */

        mul_by_8_fp2(C, C);                 /* 8*C */
        Fp2::sub(T.Y, D, T.X);             /* D-X3 */
        Fp2::mul(T.Y, T.Y, E);             /* E*(D-X3) */
        Fp2::sub(T.Y, T.Y, C);             /* Y3 = E*(D-X3)-8*C */

        /*
         * line evaluation
         */
        Fp2::sqr(line[0], line[0]);
        Fp2::sub(line[0], line[0], A);
        Fp2::sub(line[0], line[0], F);       /* (3*A+X1)^2 - X1^2 - 9*A^2 */
        lshift_fp2(B, B, 2);
        Fp2::sub(line[0], line[0], B);       /* 6*X1^3 - 4*Y1^2 */

        Fp2::mul(line[1], E, ZZ);            /* 3*X1^2 * Z1^2 */

        Fp2::mul(line[2], T.Z, ZZ);         /* Z3 * Z1^2 */
    }

    BLST_HD void line_by_Px2(Fp6& line, const POINTonE1_affine& Px2)
    {
        Fp::mul(line[1][0], line[1][0], Px2.X);   /* "b01" *= -2*P.X */
        Fp::mul(line[1][1], line[1][1], Px2.X);

        Fp::mul(line[2][0], line[2][0], Px2.Y);   /* "b11" *= 2*P.Y */
        Fp::mul(line[2][1], line[2][1], Px2.Y);
    }

#if 0
    BLST_HD  void add_n_dbl(Fp12* ret, POINTonE2& T, const POINTonE2_affine* Q,
        const POINTonE1_affine* Px2, Fp6& line, size_t n)
    {
        line_add(line, T, T, Q);    line_by_Px2(line, Px2);
        mul_by_xy00z0_fp12(ret, ret, line);
        while (n--) {
            sqr_fp12(ret, ret);
            line_dbl(line, T, T);   line_by_Px2(line, Px2);
            mul_by_xy00z0_fp12(ret, ret, line);
        }
    }

    BLST_HD void miller_loop(Fp12* ret, const POINTonE2& Q, const POINTonE1* P)
    {
#define Q ((const POINTonE2_affine *)Q)
        POINTonE2 T[1];
        POINTonE1_affine Px2[1];
        Fp6& line; /* it's not actual fp6, but 3 packed fp2, "xy00z0"  */

        /* Move common expression from line evaluation to line_by_Px2. */
        Fp::add(Px2.X, P.X, P.X);
        Fp::neg(Px2.X, Px2.X);
        Fp::add(Px2.Y, P.Y, P.Y);

        vec_copy(T.X, Q.X, 2 * sizeof(T.X));
        vec_copy(T.Z, get_BLST_Rx().p2, sizeof(T.Z));

        /* first step is ret = 1^2*line, which is replaced with ret = line  */
        line_dbl(line, T, T);                       /* 0x2                  */
        line_by_Px2(line, Px2);
        vec_zero(ret, sizeof(Fp12));
        vec_copy(ret[0][0], line[0], 2 * sizeof(Fp2));
        vec_copy(ret[1][1], line[2], sizeof(Fp2));
        add_n_dbl(ret, T, Q, Px2, line, 2);         /* ..0xc                */
        add_n_dbl(ret, T, Q, Px2, line, 3);         /* ..0x68               */
        add_n_dbl(ret, T, Q, Px2, line, 9);         /* ..0xd200             */
        add_n_dbl(ret, T, Q, Px2, line, 32);        /* ..0xd20100000000     */
        add_n_dbl(ret, T, Q, Px2, line, 16);        /* ..0xd201000000010000 */
        conjugate_fp12(ret);                /* account for z being negative */
#undef Q
    }
#endif

    BLST_HD void start_dbl_n(Fp12& ret, POINTonE2 T[],
        const POINTonE1_affine Px2[], size_t n)
    {
        size_t i;
        Fp6 line; /* it's not actual fp6, but 3 packed fp2, "xy00z0"  */

        /* first step is ret = 1^2*line, which is replaced with ret = line  */
        line_dbl(line, T[0], T[0]);           
        line_by_Px2(line, Px2[0]);
        vec_zero(&ret, sizeof(Fp12));
        vec_copy(&ret[0][0], &line[0], 2 * sizeof(Fp2));
        vec_copy(&ret[1][1], &line[2], sizeof(Fp2));

        for (i = 1; i < n; i++) {
            line_dbl(line, T[i], T[i]);       
            line_by_Px2(line, Px2[i]);
            mul_by_xy00z0_fp12(ret, ret, line);
        }
    }

    BLST_HD  void add_n_dbl_n(Fp12& ret, 
        POINTonE2 T[],
        const POINTonE2_affine Q[],
        const POINTonE1_affine Px2[],
        size_t n, 
        size_t k)
    {
        size_t i;
        Fp6 line; /* it's not actual fp6, but 3 packed fp2, "xy00z0"  */

        for (i = 0; i < n; i++) {
            line_add(line, T[i], T[i], Q[i]);  line_by_Px2(line, Px2[i]);
            mul_by_xy00z0_fp12(ret, ret, line);
        }
        while (k--) {
            sqr_fp12(ret, ret);
            for (i = 0; i < n; i++) {
                line_dbl(line, T[i], T[i]);   line_by_Px2(line, Px2[i]);
                mul_by_xy00z0_fp12(ret, ret, line);
            }
        }
    }

    BLST_HD void miller_loop_n(
        Fp12& ret, 
        const POINTonE2_affine Q[],
        const POINTonE1_affine P[], 
        size_t n)
    {
//#if !defined(__STDC_VERSION__) || __STDC_VERSION__<199901 \
//                               || defined(__STDC_NO_VLA__)
//        POINTonE2& T = (POINTonE2&)alloca(n * sizeof(POINTonE2));
//        POINTonE1_affine* Px2 = (POINTonE1_affine*)alloca(n * sizeof(POINTonE1_affine));
//#else
        if (n > MILLER_LOOP_N_MAX)
            return;

        POINTonE2 T[MILLER_LOOP_N_MAX];
        POINTonE1_affine Px2[MILLER_LOOP_N_MAX];
//#endif
        size_t i;

        if ((n == 1) && (vec_is_zero(&Q[0], sizeof(Q[0])) |
            vec_is_zero(&P[0], sizeof(P[0])))) {
            /*
             * Special case of infinite aggregated signature, pair the additive
             * group's identity with the multiplicative group's identity.
             */
            vec_copy(&ret, get_BLST_Rx().p12, sizeof(Fp12));
            return;
        }

        for (i = 0; i < n; i++) {
            /* Move common expression from line evaluation to line_by_Px2.  */
            Fp::add(Px2[i].X, P[i].X, P[i].X);
            Fp::neg(Px2[i].X, Px2[i].X);
            Fp::add(Px2[i].Y, P[i].Y, P[i].Y);

            vec_copy(&T[i].X, &Q[i].X, 2 * sizeof(T[i].X));
            vec_copy(&T[i].Z, &get_BLST_Rx().p2, sizeof(T[i].Z));
        }

        /* first step is ret = 1^2*line, which is replaced with ret = line  */
        start_dbl_n(ret, T, Px2, n);                /* 0x2                  */
        add_n_dbl_n(ret, T, Q, Px2, n, 2);          /* ..0xc                */
        add_n_dbl_n(ret, T, Q, Px2, n, 3);          /* ..0x68               */
        add_n_dbl_n(ret, T, Q, Px2, n, 9);          /* ..0xd200             */
        add_n_dbl_n(ret, T, Q, Px2, n, 32);         /* ..0xd20100000000     */
        add_n_dbl_n(ret, T, Q, Px2, n, 16);         /* ..0xd201000000010000 */
        conjugate_fp12(ret);                /* account for z being negative */
    }


    BLST_HD void pre_add_n_dbl(Fp6 lines[], POINTonE2& T,
        const POINTonE2_affine& Q,
        size_t n)
    {
        line_add(lines++[0], T, T, Q);
        while (n--)
            line_dbl(lines++[0], T, T);
    }

    BLST_HD void precompute_lines(Fp6 Qlines[68], const POINTonE2_affine& Q)
    {
        POINTonE2 T;

        vec_copy(&T.X, &Q.X, 2 * sizeof(T.X));
        vec_copy(&T.Z, &get_BLST_Rx().p2, sizeof(T.Z));

        line_dbl(Qlines[0], T, T);                  /* 0x2                  */
        pre_add_n_dbl(&Qlines[1], T, Q, 2);        /* ..0xc                */
        pre_add_n_dbl(&Qlines[4], T, Q, 3);        /* ..0x68               */
        pre_add_n_dbl(&Qlines[8], T, Q, 9);        /* ..0xd200             */
        pre_add_n_dbl(&Qlines[18], T, Q, 32);       /* ..0xd20100000000     */
        pre_add_n_dbl(&Qlines[51], T, Q, 16);       /* ..0xd201000000010000 */
    }

    BLST_HD void post_line_by_Px2(Fp6& out, const Fp6& in,
        const POINTonE1_affine& Px2)
    {
        vec_copy(&out[0], &in[0], sizeof(out[0]));

        Fp::mul(out[1][0], in[1][0], Px2.X);        /* "b01" *= -2*P.X */
        Fp::mul(out[1][1], in[1][1], Px2.X);

        Fp::mul(out[2][0], in[2][0], Px2.Y);        /* "b11" *= 2*P.Y */
        Fp::mul(out[2][1], in[2][1], Px2.Y);
    }

    BLST_HD void post_add_n_dbl(Fp12& ret, const Fp6 lines[],
        const POINTonE1_affine& Px2, size_t n)
    {
        Fp6 line;

        post_line_by_Px2(line, lines++[0], Px2);
        mul_by_xy00z0_fp12(ret, ret, line);
        while (n--) {
            sqr_fp12(ret, ret);
            post_line_by_Px2(line, lines++[0], Px2);
            mul_by_xy00z0_fp12(ret, ret, line);
        }
    }

    BLST_HD void miller_loop_lines(Fp12& ret, const Fp6 Qlines[68],
        const POINTonE1_affine& P)
    {
        POINTonE1_affine Px2;
        Fp6 line; /* it's not actual fp6, but 3 packed fp2, "xy00z0"  */

        /* Move common expression from line evaluation to line_by_Px2. */
        Fp::add(Px2.X, P.X, P.X);
        Fp::neg(Px2.X, Px2.X);
        Fp::add(Px2.Y, P.Y, P.Y);

        /* first step is ret = 1^2*line, which is replaced with ret = line  */
        post_line_by_Px2(line, Qlines[0], Px2);     /* 0x2                  */
        vec_zero(&ret, sizeof(Fp12));
        vec_copy(&ret[0][0], &line[0], 2 * sizeof(Fp2));
        vec_copy(&ret[1][1], &line[2], sizeof(Fp2));
        post_add_n_dbl(ret, &Qlines[1], Px2, 2);   /* ..0xc                */
        post_add_n_dbl(ret, &Qlines[4], Px2, 3);   /* ..0x68               */
        post_add_n_dbl(ret, &Qlines[8], Px2, 9);   /* ..0xd200             */
        post_add_n_dbl(ret, &Qlines[18], Px2, 32);  /* ..0xd20100000000     */
        post_add_n_dbl(ret, &Qlines[51], Px2, 16);  /* ..0xd201000000010000 */
        conjugate_fp12(ret);                /* account for z being negative */
    }

#ifdef INTERNAL_TESTMODE
    void miller_loop_alt(Fp12* ret, const POINTonE2_affine* Q,
        const POINTonE1_affine* P)
    {
        Fp6& lines[68];

        precompute_lines(lines, Q);
        miller_loop_lines(ret, lines, P);
    }
#endif

    BLST_HD void mul_n_sqr(Fp12& ret, const Fp12& a, size_t n)
    {
        mul_fp12(ret, ret, a);
        while (n--)
            cyclotomic_sqr_fp12(ret, ret);
    }

    BLST_HD void raise_to_z_div_by_2(Fp12& ret, const Fp12& a)
    {
        cyclotomic_sqr_fp12(ret, a);                /* 0x2                  */
        mul_n_sqr(ret, a, 2);                       /* ..0xc                */
        mul_n_sqr(ret, a, 3);                       /* ..0x68               */
        mul_n_sqr(ret, a, 9);                       /* ..0xd200             */
        mul_n_sqr(ret, a, 32);                      /* ..0xd20100000000     */
        mul_n_sqr(ret, a, 16 - 1);                    /* ..0x6900800000008000 */
        conjugate_fp12(ret);                /* account for z being negative */
    }

#define raise_to_z(a, b) (raise_to_z_div_by_2(a, b), cyclotomic_sqr_fp12(a, a))

    /*
     * Adaptation from <zkcrypto>/pairing/src/bls12_381/mod.rs
     */
    BLST_HD void final_exp(Fp12& ret, const Fp12& f)
    {
        Fp12 y0, y1, y2, y3;

        vec_copy(&y1, &f, sizeof(y1));
        conjugate_fp12(y1);
        inverse_fp12(y2, f);
        mul_fp12(ret, y1, y2);
        frobenius_map_fp12(y2, ret, 2);
        mul_fp12(ret, ret, y2);

        cyclotomic_sqr_fp12(y0, ret);
        raise_to_z(y1, y0);
        raise_to_z_div_by_2(y2, y1);
        vec_copy(&y3, &ret, sizeof(y3));
        conjugate_fp12(y3);
        mul_fp12(y1, y1, y3);
        conjugate_fp12(y1);
        mul_fp12(y1, y1, y2);
        raise_to_z(y2, y1);
        raise_to_z(y3, y2);
        conjugate_fp12(y1);
        mul_fp12(y3, y3, y1);
        conjugate_fp12(y1);
        frobenius_map_fp12(y1, y1, 3);
        frobenius_map_fp12(y2, y2, 2);
        mul_fp12(y1, y1, y2);
        raise_to_z(y2, y3);
        mul_fp12(y2, y2, y0);
        mul_fp12(y2, y2, ret);
        mul_fp12(y1, y1, y2);
        frobenius_map_fp12(y2, y3, 1);
        mul_fp12(ret, y1, y2);
    }

    BLST_HD void blst_miller_loop(Fp12& ret, const POINTonE2_affine* Q,
        const POINTonE1_affine* P)
    {
        miller_loop_n(ret, 
            Q ? Q : (const POINTonE2_affine*)&get_BLST_G2(),
            P ? P : (const POINTonE1_affine*)&get_BLST_G1(), 
            1);
    }


    BLST_HD void blst_miller_loop_n(
        Fp12& out, 
        const POINTonE2_affine* const Qs[],
        const POINTonE1_affine* const Ps[],
        size_t n)
    {   /* ~10KB of stack storage */
        POINTonE2 T[MILLER_LOOP_N_MAX];
        POINTonE2_affine Q[MILLER_LOOP_N_MAX];
        POINTonE1_affine Px2[MILLER_LOOP_N_MAX];
        const POINTonE2_affine* Qptr = NULL;
        const POINTonE1_affine* Pptr = NULL;
        size_t i, j;

        for (i = 0, j = 0; j < n; j++) {
            Qptr = *Qs ? *Qs++ : Qptr + 1;
            Pptr = *Ps ? *Ps++ : Pptr + 1;

            /* Move common expression from line evaluation to line_by_Px2.  */
            Fp::add(Px2[i].X, Pptr->X, Pptr->X);
            Fp::neg(Px2[i].X, Px2[i].X);
            Fp::add(Px2[i].Y, Pptr->Y, Pptr->Y);

            vec_copy(&Q[i].X, &Qptr->X, 2 * sizeof(Q[i].X));
            vec_copy(&T[i].X, &Qptr->X, 2 * sizeof(T[i].X));
            vec_copy(&T[i].Z, &get_BLST_Rx().p2, sizeof(T[i].Z));

            if (++i == MILLER_LOOP_N_MAX || j == n - 1) {
                Fp12 tmp;
                Fp12& ret = j < MILLER_LOOP_N_MAX ? out : tmp;

                /* first step is ret = 1^2*line, which is just ret = line       */
                start_dbl_n(ret, T, Px2, i);            /* 0x2                  */
                add_n_dbl_n(ret, T, Q, Px2, i, 2);      /* ..0xc                */
                add_n_dbl_n(ret, T, Q, Px2, i, 3);      /* ..0x68               */
                add_n_dbl_n(ret, T, Q, Px2, i, 9);      /* ..0xd200             */
                add_n_dbl_n(ret, T, Q, Px2, i, 32);     /* ..0xd20100000000     */
                add_n_dbl_n(ret, T, Q, Px2, i, 16);     /* ..0xd201000000010000 */
                conjugate_fp12(ret);            /* account for z being negative */

                if (j >= MILLER_LOOP_N_MAX)
                    mul_fp12(out, out, ret);

                i = 0;
            }
        }
    }

    BLST_HD void blst_final_exp(Fp12& ret, const Fp12& f)
    {
        final_exp(ret, f);
    }

    BLST_HD void blst_precompute_lines(Fp6 Qlines[68], const POINTonE2_affine& Q)
    {
        precompute_lines(Qlines, Q);
    }

    BLST_HD  void blst_miller_loop_lines(Fp12& ret, const Fp6 Qlines[68],
        const POINTonE1_affine& P)
    {
        miller_loop_lines(ret, Qlines, P);
    }

    BLST_HD  bool_t is_cyclotomic(const Fp12& f)
    {
        Fp12 a, b;

        frobenius_map_fp12(a, f, 2);
        frobenius_map_fp12(b, a, 2);
        mul_fp12(b, b, f);

        return vec_is_equal(&a, &b, sizeof(a));
    }

    BLST_HD int blst_fp12_in_group(const Fp12& f)
    {
        Fp12 a, b;

        if (vec_is_zero(&f, sizeof(Fp12)) || !is_cyclotomic(f))
            return 0;

        frobenius_map_fp12(a, f, 1);
        raise_to_z(b, f);

        return (int)vec_is_equal(&a, &b, sizeof(a));
    }
}