/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/scalar/curve/point.h"
#include "blst/scalar/curve/GT.h"

#ifndef MILLER_LOOP_N_MAX
# define MILLER_LOOP_N_MAX 16
#endif

namespace blst::scalar
{
    /*
    * Line evaluations from  https://eprint.iacr.org/2010/354.pdf
    * with a twist moving common expression to line_by_Px2.
    */
    BLST_HD void line_add(Fp6& line, G2Point& T, const G2Point& R,
        const G2AffinePoint& Q)
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

    BLST_HD void line_dbl(Fp6& line, G2Point& T, const G2Point& Q)
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

        Fp2::mul_by_3(E, A);                 /* E = 3*A */
        Fp2::sqr(F, E);                      /* F = E^2 */

        Fp2::add(line[0], E, Q.X);          /* 3*A+X1 for line evaluation */

        Fp2::sub(T.X, F, D);
        Fp2::sub(T.X, T.X, D);             /* X3 = F-2*D */

        Fp2::add(T.Z, Q.Y, Q.Z);
        Fp2::sqr(T.Z, T.Z);
        Fp2::sub(T.Z, T.Z, B);
        Fp2::sub(T.Z, T.Z, ZZ);            /* Z3 = (Y1+Z1)^2-B-ZZ */

        Fp2::mul_by_8(C, C);                 /* 8*C */
        Fp2::sub(T.Y, D, T.X);             /* D-X3 */
        Fp2::mul(T.Y, T.Y, E);             /* E*(D-X3) */
        Fp2::sub(T.Y, T.Y, C);             /* Y3 = E*(D-X3)-8*C */

        /*
         * line evaluation
         */
        Fp2::sqr(line[0], line[0]);
        Fp2::sub(line[0], line[0], A);
        Fp2::sub(line[0], line[0], F);       /* (3*A+X1)^2 - X1^2 - 9*A^2 */
        Fp2::lshift(B, B, 2);
        Fp2::sub(line[0], line[0], B);       /* 6*X1^3 - 4*Y1^2 */

        Fp2::mul(line[1], E, ZZ);            /* 3*X1^2 * Z1^2 */

        Fp2::mul(line[2], T.Z, ZZ);         /* Z3 * Z1^2 */
    }

    BLST_HD void line_by_Px2(Fp6& line, const G1AffinePoint& Px2)
    {
        Fp::mul(line[1][0], line[1][0], Px2.X);   /* "b01" *= -2*P.X */
        Fp::mul(line[1][1], line[1][1], Px2.X);

        Fp::mul(line[2][0], line[2][0], Px2.Y);   /* "b11" *= 2*P.Y */
        Fp::mul(line[2][1], line[2][1], Px2.Y);
    }

//#if 0
//    BLST_HD  void add_n_dbl(Fp12* ret, G2Point& T, const G2AffinePoint* Q,
//        const G1AffinePoint* Px2, Fp6& line, size_t n)
//    {
//        line_add(line, T, T, Q);    line_by_Px2(line, Px2);
//        Fp12::mul_by_xy00z0(ret, ret, line);
//        while (n--) {
//            sqr_fp12(ret, ret);
//            line_dbl(line, T, T);   line_by_Px2(line, Px2);
//            Fp12::mul_by_xy00z0(ret, ret, line);
//        }
//    }
//
//    BLST_HD void miller_loop(Fp12* ret, const G2Point& Q, const G1Point* P)
//    {
//#define Q ((const G2AffinePoint *)Q)
//        G2Point T[1];
//        G1AffinePoint Px2[1];
//        Fp6& line; /* it's not actual fp6, but 3 packed fp2, "xy00z0"  */
//
//        /* Move common expression from line evaluation to line_by_Px2. */
//        Fp::add(Px2.X, P.X, P.X);
//        Fp::neg(Px2.X, Px2.X);
//        Fp::add(Px2.Y, P.Y, P.Y);
//
//        vec_copy(T.X, Q.X, 2 * sizeof(T.X));
//        vec_copy(T.Z, get_BLST_Rx().p2, sizeof(T.Z));
//
//        /* first step is ret = 1^2*line, which is replaced with ret = line  */
//        line_dbl(line, T, T);                       /* 0x2                  */
//        line_by_Px2(line, Px2);
//        vec_zero(ret, sizeof(Fp12));
//        vec_copy(ret[0][0], line[0], 2 * sizeof(Fp2));
//        vec_copy(ret[1][1], line[2], sizeof(Fp2));
//        add_n_dbl(ret, T, Q, Px2, line, 2);         /* ..0xc                */
//        add_n_dbl(ret, T, Q, Px2, line, 3);         /* ..0x68               */
//        add_n_dbl(ret, T, Q, Px2, line, 9);         /* ..0xd200             */
//        add_n_dbl(ret, T, Q, Px2, line, 32);        /* ..0xd20100000000     */
//        add_n_dbl(ret, T, Q, Px2, line, 16);        /* ..0xd201000000010000 */
//        Fp12::conjugate(ret);                /* account for z being negative */
//#undef Q
//    }
//#endif

    BLST_HD void start_dbl_n(
        GT& ret, 
        G2Point T[],
        const G1AffinePoint Px2[], size_t n)
    {
        size_t i;
        Fp6 line; /* it's not actual fp6, but 3 packed fp2, "xy00z0"  */

        /* first step is ret = 1^2*line, which is replaced with ret = line  */
        line_dbl(line, T[0], T[0]);           
        line_by_Px2(line, Px2[0]);
        vec_zero(&ret, sizeof(GT));
        vec_copy(&ret.data[0][0], &line[0], 2 * sizeof(Fp2));
        vec_copy(&ret.data[1][1], &line[2], sizeof(Fp2));

        for (i = 1; i < n; i++) {
            line_dbl(line, T[i], T[i]);       
            line_by_Px2(line, Px2[i]);
            Fp12::mul_by_xy00z0(ret.data, ret.data, line);
        }
    }

    BLST_HD  void add_n_dbl_n(GT& ret, 
        G2Point T[],
        const G2AffinePoint Q[],
        const G1AffinePoint Px2[],
        size_t n, 
        size_t k)
    {
        size_t i;
        Fp6 line; /* it's not actual fp6, but 3 packed fp2, "xy00z0"  */

        for (i = 0; i < n; i++) {
            line_add(line, T[i], T[i], Q[i]);  line_by_Px2(line, Px2[i]);
            Fp12::mul_by_xy00z0(ret.data, ret.data, line);
        }
        while (k--) {
            Fp12::sqr(ret.data, ret.data);
            for (i = 0; i < n; i++) {
                line_dbl(line, T[i], T[i]);   line_by_Px2(line, Px2[i]);
                Fp12::mul_by_xy00z0(ret.data, ret.data, line);
            }
        }
    }

    BLST_HD void GT::miller_loop_n(
        GT& ret, 
        const G2AffinePoint Q[],
        const G1AffinePoint P[], 
        size_t n)
    {
//#if !defined(__STDC_VERSION__) || __STDC_VERSION__<199901 \
//                               || defined(__STDC_NO_VLA__)
//        G2Point& T = (G2Point&)alloca(n * sizeof(G2Point));
//        G1AffinePoint* Px2 = (G1AffinePoint*)alloca(n * sizeof(G1AffinePoint));
//#else
        if (n > MILLER_LOOP_N_MAX)
            return;

        G2Point T[MILLER_LOOP_N_MAX];
        G1AffinePoint Px2[MILLER_LOOP_N_MAX];
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
        Fp12::conjugate(ret.data);                /* account for z being negative */
    }


    BLST_HD void pre_add_n_dbl(Fp6 lines[], G2Point& T,
        const G2AffinePoint& Q,
        size_t n)
    {
        line_add(lines++[0], T, T, Q);
        while (n--)
            line_dbl(lines++[0], T, T);
    }

    BLST_HD void GT::precompute_lines(Fp6 Qlines[68], const G2AffinePoint& Q)
    {
        G2Point T;

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
        const G1AffinePoint& Px2)
    {
        vec_copy(&out[0], &in[0], sizeof(out[0]));

        Fp::mul(out[1][0], in[1][0], Px2.X);        /* "b01" *= -2*P.X */
        Fp::mul(out[1][1], in[1][1], Px2.X);

        Fp::mul(out[2][0], in[2][0], Px2.Y);        /* "b11" *= 2*P.Y */
        Fp::mul(out[2][1], in[2][1], Px2.Y);
    }

    BLST_HD void post_add_n_dbl(GT& ret, const Fp6 lines[],
        const G1AffinePoint& Px2, size_t n)
    {
        Fp6 line;

        post_line_by_Px2(line, lines++[0], Px2);
        Fp12::mul_by_xy00z0(ret.data, ret.data, line);
        while (n--) {
            Fp12::sqr(ret.data, ret.data);
            post_line_by_Px2(line, lines++[0], Px2);
            Fp12::mul_by_xy00z0(ret.data, ret.data, line);
        }
    }

    BLST_HD void GT::miller_loop_lines(GT& ret, const Fp6 Qlines[68],
        const G1AffinePoint& P)
    {
        G1AffinePoint Px2;
        Fp6 line; /* it's not actual fp6, but 3 packed fp2, "xy00z0"  */

        /* Move common expression from line evaluation to line_by_Px2. */
        Fp::add(Px2.X, P.X, P.X);
        Fp::neg(Px2.X, Px2.X);
        Fp::add(Px2.Y, P.Y, P.Y);

        /* first step is ret = 1^2*line, which is replaced with ret = line  */
        post_line_by_Px2(line, Qlines[0], Px2);     /* 0x2                  */
        vec_zero(&ret, sizeof(GT));
        vec_copy(&ret.data[0][0], &line[0], 2 * sizeof(Fp2));
        vec_copy(&ret.data[1][1], &line[2], sizeof(Fp2));
        post_add_n_dbl(ret, &Qlines[1], Px2, 2);   /* ..0xc                */
        post_add_n_dbl(ret, &Qlines[4], Px2, 3);   /* ..0x68               */
        post_add_n_dbl(ret, &Qlines[8], Px2, 9);   /* ..0xd200             */
        post_add_n_dbl(ret, &Qlines[18], Px2, 32);  /* ..0xd20100000000     */
        post_add_n_dbl(ret, &Qlines[51], Px2, 16);  /* ..0xd201000000010000 */
        Fp12::conjugate(ret.data);                /* account for z being negative */
    }

#ifdef INTERNAL_TESTMODE
    void miller_loop_alt(GT* ret, const G2AffinePoint* Q,
        const G1AffinePoint* P)
    {
        Fp6& lines[68];

        precompute_lines(lines, Q);
        miller_loop_lines(ret, lines, P);
    }
#endif


    BLST_HD void GT::miller_loop(GT& ret, const G2AffinePoint* Q,
        const G1AffinePoint* P)
    {
        miller_loop_n(ret, 
            Q ? Q : (const G2AffinePoint*)&get_BLST_G2(),
            P ? P : (const G1AffinePoint*)&get_BLST_G1(), 
            1);
    }


    BLST_HD void GT::miller_loop_n(
        GT& out,
        const G2AffinePoint* const Qs[],
        const G1AffinePoint* const Ps[],
        size_t n)
    {   /* ~10KB of stack storage */
        G2Point T[MILLER_LOOP_N_MAX];
        G2AffinePoint Q[MILLER_LOOP_N_MAX];
        G1AffinePoint Px2[MILLER_LOOP_N_MAX];
        const G2AffinePoint* Qptr = NULL;
        const G1AffinePoint* Pptr = NULL;
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
                GT tmp;
                GT& ret = j < MILLER_LOOP_N_MAX ? out : tmp;

                /* first step is ret = 1^2*line, which is just ret = line       */
                start_dbl_n(ret, T, Px2, i);            /* 0x2                  */
                add_n_dbl_n(ret, T, Q, Px2, i, 2);      /* ..0xc                */
                add_n_dbl_n(ret, T, Q, Px2, i, 3);      /* ..0x68               */
                add_n_dbl_n(ret, T, Q, Px2, i, 9);      /* ..0xd200             */
                add_n_dbl_n(ret, T, Q, Px2, i, 32);     /* ..0xd20100000000     */
                add_n_dbl_n(ret, T, Q, Px2, i, 16);     /* ..0xd201000000010000 */
                Fp12::conjugate(ret.data);            /* account for z being negative */

                if (j >= MILLER_LOOP_N_MAX)
                    Fp12::mul(out.data, out.data, ret.data);

                i = 0;
            }
        }
    }

    //BLST_HD void blst_precompute_lines(Fp6 Qlines[68], const G2AffinePoint& Q)
    //{
    //    precompute_lines(Qlines, Q);
    //}

    //BLST_HD  void blst_miller_loop_lines(GT& ret, const Fp6 Qlines[68],
    //    const G1AffinePoint& P)
    //{
    //    miller_loop_lines(ret, Qlines, P);
    //}

}
