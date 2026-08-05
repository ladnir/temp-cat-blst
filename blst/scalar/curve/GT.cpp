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
#define BLST_GT_GENERATOR_LIMBS                                           \
Fp6{                                                                      \
    Fp2{                                                                  \
        Fp{ TO_LIMB_T(0x1972e433a01f85c5), TO_LIMB_T(0x97d32b76fd772538), \
            TO_LIMB_T(0xc8ce546fc96bcdf9), TO_LIMB_T(0xcef63e7366d40614), \
            TO_LIMB_T(0xa611342781843780), TO_LIMB_T(0x13f3448a3fc6d825) },\
        Fp{ TO_LIMB_T(0xd26331b02e9d6995), TO_LIMB_T(0x9d68a482f7797e7d), \
            TO_LIMB_T(0x9c9b29248d39ea92), TO_LIMB_T(0xf4801ca2e13107aa), \
            TO_LIMB_T(0xa16c0732bdbcb066), TO_LIMB_T(0x083ca4afba360478) } \
    },                                                                    \
    Fp2{                                                                  \
        Fp{ TO_LIMB_T(0x59e261db0916b641), TO_LIMB_T(0x2716b6f4b23e960d), \
            TO_LIMB_T(0xc8e55b10a0bd9c45), TO_LIMB_T(0x0bdb0bd99c4deda8), \
            TO_LIMB_T(0x8cf89ebf57fdaac5), TO_LIMB_T(0x12d6b7929e777a5e) },\
        Fp{ TO_LIMB_T(0x5fc85188b0e15f35), TO_LIMB_T(0x34a06e3a8f096365), \
            TO_LIMB_T(0xdb3126a6e02ad62c), TO_LIMB_T(0xfc6f5aa97d9a990b), \
            TO_LIMB_T(0xa12f55f5eb89c210), TO_LIMB_T(0x1723703a926f8889) } \
    },                                                                    \
    Fp2{                                                                  \
        Fp{ TO_LIMB_T(0x93588f2971828778), TO_LIMB_T(0x43f65b8611ab7585), \
            TO_LIMB_T(0x3183aaf5ec279fdf), TO_LIMB_T(0xfa73d7e18ac99df6), \
            TO_LIMB_T(0x64e176a6a64c99b0), TO_LIMB_T(0x179fa78c58388f1f) },\
        Fp{ TO_LIMB_T(0x672a0a11ca2aef12), TO_LIMB_T(0x0d11b9b52aa3f16b), \
            TO_LIMB_T(0xa44412d0699d056e), TO_LIMB_T(0xc01d0177221a5ba5), \
            TO_LIMB_T(0x66e0cede6c735529), TO_LIMB_T(0x05f5a71e9fddc339) } \
    }                                                                     \
},                                                                        \
Fp6{                                                                      \
    Fp2{                                                                  \
        Fp{ TO_LIMB_T(0xd30a88a1b062c679), TO_LIMB_T(0x5ac56a5d35fc8304), \
            TO_LIMB_T(0xd0c834a6a81f290d), TO_LIMB_T(0xcd5430c2da3707c7), \
            TO_LIMB_T(0xf0c27ff780500af0), TO_LIMB_T(0x09245da6e2d72eae) },\
        Fp{ TO_LIMB_T(0x9f2e0676791b5156), TO_LIMB_T(0xe2d1c8234918fe13), \
            TO_LIMB_T(0x4c9e459f3c561bf4), TO_LIMB_T(0xa3e85e53b9d3e3c1), \
            TO_LIMB_T(0x820a121e21a70020), TO_LIMB_T(0x15af618341c59acc) } \
    },                                                                    \
    Fp2{                                                                  \
        Fp{ TO_LIMB_T(0x7c95658c24993ab1), TO_LIMB_T(0x73eb38721ca886b9), \
            TO_LIMB_T(0x5256d749477434bc), TO_LIMB_T(0x8ba41902ea504a8b), \
            TO_LIMB_T(0x04a3d3f80c86ce6d), TO_LIMB_T(0x18a64a87fb686eaa) },\
        Fp{ TO_LIMB_T(0xbb83e71bb920cf26), TO_LIMB_T(0x2a5277ac92a73945), \
            TO_LIMB_T(0xfc0ee59f94f046a0), TO_LIMB_T(0x7158cdf3786058f7), \
            TO_LIMB_T(0x7cc1061b82f945f6), TO_LIMB_T(0x03f847aa9fdbe567) } \
    },                                                                    \
    Fp2{                                                                  \
        Fp{ TO_LIMB_T(0x8078dba56134e657), TO_LIMB_T(0x1cd7ec9a43998a6e), \
            TO_LIMB_T(0xb1aa599a1a993766), TO_LIMB_T(0xc9a0f62f0842ee44), \
            TO_LIMB_T(0x8e159be3b605dffa), TO_LIMB_T(0x0c86ba0d4af13fc2) },\
        Fp{ TO_LIMB_T(0xe80ff2a06a52ffb1), TO_LIMB_T(0x7694ca48721a906c), \
            TO_LIMB_T(0x7583183e03b08514), TO_LIMB_T(0xf567afdd40cee4e2), \
            TO_LIMB_T(0x9a6d96d2e526a5fc), TO_LIMB_T(0x197e9f49861f2242) } \
    }                                                                     \
}

    constexpr GT BLST_GT_GENERATOR_HOST =
        Fp12{ BLST_GT_GENERATOR_LIMBS };

    const GT& GT::generator()
    {
        return BLST_GT_GENERATOR_HOST;
    }

    /*
    * Line evaluations from  https://eprint.iacr.org/2010/354.pdf
    * with a twist moving common expression to line_by_Px2.
    */
    BLST_HD void line_add(Fp6& line, G2Point& T, const G2Point& R,
        const G2AffinePoint& Q)
    {
        Fp2 Z1Z1, U2, S2, H, HH, I, J, V;
# define r line[1]

        /*
         * https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
         * with XYZ3 being |T|, XYZ1 - |R|, XY2 - |Q|, i.e. Q is affine
         */
        Fp2::sqr(Z1Z1, R.z);                /* Z1Z1 = Z1^2 */
        Fp2::mul(U2, Q.x, Z1Z1);            /* U2 = X2*Z1Z1 */

        Fp2::mul(S2, Q.y, R.z);
        Fp2::mul(S2, S2, Z1Z1);              /* S2 = Y2*Z1*Z1Z1 */

        Fp2::sub(H, U2, R.x);               /* H = U2-X1 */

        Fp2::sqr(HH, H);                     /* HH = H^2 */
        Fp2::add(I, HH, HH);
        Fp2::add(I, I, I);                   /* I = 4*HH */

        Fp2::mul(J, H, I);                   /* J = H*I */

        Fp2::sub(r, S2, R.y);
        Fp2::add(r, r, r);                   /* r = 2*(S2-Y1) */

        Fp2::mul(V, R.x, I);                /* V = X1*I */

        Fp2::sqr(T.x, r);
        Fp2::sub(T.x, T.x, J);
        Fp2::sub(T.x, T.x, V);
        Fp2::sub(T.x, T.x, V);             /* X3 = r^2-J-2*V */

        Fp2::mul(J, J, R.y);
        Fp2::sub(T.y, V, T.x);
        Fp2::mul(T.y, T.y, r);
        Fp2::sub(T.y, T.y, J);
        Fp2::sub(T.y, T.y, J);             /* Y3 = r*(V-X3)-2*Y1*J */

        Fp2::add(T.z, R.z, H);
        Fp2::sqr(T.z, T.z);
        Fp2::sub(T.z, T.z, Z1Z1);
        Fp2::sub(T.z, T.z, HH);            /* Z3 = (Z1+H)^2-Z1Z1-HH */

        /*
         * line evaluation
         */
        Fp2::mul(I, r, Q.x);
        Fp2::mul(J, Q.y, T.z);
        Fp2::sub(I, I, J);
        Fp2::add(line[0], I, I);          /* 2*(r*X2 - Y2*Z3) */
#ifdef r
# undef r
#else
        vec_copy(&line[1], &r, sizeof(r));
#endif
        vec_copy(&line[2], &T.z, sizeof(T.z));
    }

    BLST_HD void line_dbl(Fp6& line, G2Point& T, const G2Point& Q)
    {
        Fp2 ZZ, A, B, C, D, E, F;

        /*
         * https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-alnr
         */
        Fp2::sqr(A, Q.x);                   /* A = X1^2 */
        Fp2::sqr(B, Q.y);                   /* B = Y1^2 */
        Fp2::sqr(ZZ, Q.z);                  /* ZZ = Z1^2 */
        Fp2::sqr(C, B);                      /* C = B^2 */

        Fp2::add(D, Q.x, B);                /* X1+B */
        Fp2::sqr(D, D);                      /* (X1+B)^2 */
        Fp2::sub(D, D, A);                   /* (X1+B)^2-A */
        Fp2::sub(D, D, C);                   /* (X1+B)^2-A-C */
        Fp2::add(D, D, D);                   /* D = 2*((X1+B)^2-A-C) */

        Fp2::mul_by_3(E, A);                 /* E = 3*A */
        Fp2::sqr(F, E);                      /* F = E^2 */

        Fp2::add(line[0], E, Q.x);          /* 3*A+X1 for line evaluation */

        Fp2::sub(T.x, F, D);
        Fp2::sub(T.x, T.x, D);             /* X3 = F-2*D */

        Fp2::add(T.z, Q.y, Q.z);
        Fp2::sqr(T.z, T.z);
        Fp2::sub(T.z, T.z, B);
        Fp2::sub(T.z, T.z, ZZ);            /* Z3 = (Y1+Z1)^2-B-ZZ */

        Fp2::mul_by_8(C, C);                 /* 8*C */
        Fp2::sub(T.y, D, T.x);             /* D-X3 */
        Fp2::mul(T.y, T.y, E);             /* E*(D-X3) */
        Fp2::sub(T.y, T.y, C);             /* Y3 = E*(D-X3)-8*C */

        /*
         * line evaluation
         */
        Fp2::sqr(line[0], line[0]);
        Fp2::sub(line[0], line[0], A);
        Fp2::sub(line[0], line[0], F);       /* (3*A+X1)^2 - X1^2 - 9*A^2 */
        Fp2::lshift(B, B, 2);
        Fp2::sub(line[0], line[0], B);       /* 6*X1^3 - 4*Y1^2 */

        Fp2::mul(line[1], E, ZZ);            /* 3*X1^2 * Z1^2 */

        Fp2::mul(line[2], T.z, ZZ);         /* Z3 * Z1^2 */
    }

    BLST_HD void line_by_Px2(Fp6& line, const G1AffinePoint& Px2)
    {
        Fp::mul(line[1][0], line[1][0], Px2.x);   /* "b01" *= -2*P.x */
        Fp::mul(line[1][1], line[1][1], Px2.x);

        Fp::mul(line[2][0], line[2][0], Px2.y);   /* "b11" *= 2*P.y */
        Fp::mul(line[2][1], line[2][1], Px2.y);
    }

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
        vec_copy(&ret, &get_BLST_Rx(), sizeof(Fp12));
        if (n == 0) {
            return;
        }

        if ((n == 1) && (vec_is_zero(&Q[0], sizeof(Q[0])) |
            vec_is_zero(&P[0], sizeof(P[0])))) {
            /*
             * Special case of infinite aggregated signature, pair the additive
             * group's identity with the multiplicative group's identity.
             */
            return;
        }

        size_t offset = 0;
        while (offset < n) {
            G2Point T[MILLER_LOOP_N_MAX];
            G2AffinePoint Qbatch[MILLER_LOOP_N_MAX];
            G1AffinePoint Px2[MILLER_LOOP_N_MAX];
            GT chunk = {};
            size_t batch = 0;

            while (offset < n && batch < MILLER_LOOP_N_MAX) {
                if (G2AffinePoint::is_inf(Q[offset]) || G1AffinePoint::is_inf(P[offset])) {
                    ++offset;
                    continue;
                }

                /* Move common expression from line evaluation to line_by_Px2.  */
                Fp::add(Px2[batch].x, P[offset].x, P[offset].x);
                Fp::neg(Px2[batch].x, Px2[batch].x);
                Fp::add(Px2[batch].y, P[offset].y, P[offset].y);

                vec_copy(&Qbatch[batch].x, &Q[offset].x, 2 * sizeof(Qbatch[batch].x));
                vec_copy(&T[batch].x, &Q[offset].x, 2 * sizeof(T[batch].x));
                vec_copy(&T[batch].z, &Fp2::one(), sizeof(T[batch].z));
                ++batch;
                ++offset;
            }

            if (batch == 0) {
                continue;
            }

            /* first step is ret = 1^2*line, which is replaced with ret = line  */
            start_dbl_n(chunk, T, Px2, batch);                /* 0x2                  */
            add_n_dbl_n(chunk, T, Qbatch, Px2, batch, 2);     /* ..0xc                */
            add_n_dbl_n(chunk, T, Qbatch, Px2, batch, 3);     /* ..0x68               */
            add_n_dbl_n(chunk, T, Qbatch, Px2, batch, 9);     /* ..0xd200             */
            add_n_dbl_n(chunk, T, Qbatch, Px2, batch, 32);    /* ..0xd20100000000     */
            add_n_dbl_n(chunk, T, Qbatch, Px2, batch, 16);    /* ..0xd201000000010000 */
            Fp12::conjugate(chunk.data);                      /* account for z being negative */
            Fp12::mul(ret.data, ret.data, chunk.data);
        }
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

        vec_copy(&T.x, &Q.x, 2 * sizeof(T.x));
        vec_copy(&T.z, &Fp2::one(), sizeof(T.z));

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

        Fp::mul(out[1][0], in[1][0], Px2.x);        /* "b01" *= -2*P.x */
        Fp::mul(out[1][1], in[1][1], Px2.x);

        Fp::mul(out[2][0], in[2][0], Px2.y);        /* "b11" *= 2*P.y */
        Fp::mul(out[2][1], in[2][1], Px2.y);
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

    BLST_HD void start_lines_n(
        GT& ret,
        const Fp6* const Qlines[],
        const G1AffinePoint Px2[],
        size_t n)
    {
        BLST_ASSERT(n != 0);

        Fp6 line;
        post_line_by_Px2(line, Qlines[0][0], Px2[0]);
        vec_zero(&ret, sizeof(GT));
        vec_copy(&ret.data[0][0], &line[0], 2 * sizeof(Fp2));
        vec_copy(&ret.data[1][1], &line[2], sizeof(Fp2));

        for (size_t i = 1; i < n; ++i) {
            post_line_by_Px2(line, Qlines[i][0], Px2[i]);
            Fp12::mul_by_xy00z0(ret.data, ret.data, line);
        }
    }

    BLST_HD void post_add_n_dbl_n(
        GT& ret,
        const Fp6* const Qlines[],
        const G1AffinePoint Px2[],
        size_t n,
        size_t start,
        size_t k)
    {
        Fp6 line;
        size_t line_idx = start;

        for (size_t i = 0; i < n; ++i) {
            post_line_by_Px2(line, Qlines[i][line_idx], Px2[i]);
            Fp12::mul_by_xy00z0(ret.data, ret.data, line);
        }

        while (k--) {
            ++line_idx;
            Fp12::sqr(ret.data, ret.data);
            for (size_t i = 0; i < n; ++i) {
                post_line_by_Px2(line, Qlines[i][line_idx], Px2[i]);
                Fp12::mul_by_xy00z0(ret.data, ret.data, line);
            }
        }
    }

    BLST_HD void GT::miller_loop_lines(GT& ret, const Fp6 Qlines[68],
        const G1AffinePoint& P)
    {
        G1AffinePoint Px2;
        Fp6 line; /* it's not actual fp6, but 3 packed fp2, "xy00z0"  */

        /* Move common expression from line evaluation to line_by_Px2. */
        Fp::add(Px2.x, P.x, P.x);
        Fp::neg(Px2.x, Px2.x);
        Fp::add(Px2.y, P.y, P.y);

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

    BLST_HD void GT::miller_loop_lines_n(
        GT& out,
        const Fp6* const Qlines[],
        const G1AffinePoint* const Ps[],
        size_t n)
    {
        vec_copy(&out, &get_BLST_Rx(), sizeof(Fp12));
        if (n == 0 || Qlines == NULL || Ps == NULL || Qlines[0] == NULL || Ps[0] == NULL) {
            return;
        }

        G1AffinePoint Px2[MILLER_LOOP_N_MAX];
        const Fp6* line_batch[MILLER_LOOP_N_MAX];
        size_t batch = 0;
        const Fp6* lines_ptr = NULL;
        const G1AffinePoint* pptr = NULL;

        // Preserve the old generic/null-skipping entry point, but steer hot call
        // sites to the no-chunk preprocessed variant below.
        for (size_t i = 0; i < n; ++i) {
            lines_ptr = Qlines[i] ? Qlines[i] : lines_ptr + 68;
            pptr = Ps[i] ? Ps[i] : pptr + 1;
            BLST_ASSERT(lines_ptr != NULL);
            BLST_ASSERT(pptr != NULL);

            if (G1AffinePoint::is_inf(*pptr)) {
                continue;
            }

            line_batch[batch] = lines_ptr;
            Fp::add(Px2[batch].x, pptr->x, pptr->x);
            Fp::neg(Px2[batch].x, Px2[batch].x);
            Fp::add(Px2[batch].y, pptr->y, pptr->y);
            ++batch;

            if (batch == MILLER_LOOP_N_MAX) {
                GT chunk = {};
                GT::miller_loop_lines_n_preprocessed(chunk, line_batch, Px2, batch);
                Fp12::mul(out.data, out.data, chunk.data);
                batch = 0;
            }
        }

        if (batch == 0) {
            return;
        }

        GT tail = {};
        GT::miller_loop_lines_n_preprocessed(tail, line_batch, Px2, batch);
        Fp12::mul(out.data, out.data, tail.data);
    }

    BLST_HD void GT::miller_loop_lines_n_preprocessed(
        GT& out,
        const Fp6* const Qlines[],
        const G1AffinePoint Px2[],
        size_t n)
    {
        vec_copy(&out, &get_BLST_Rx(), sizeof(Fp12));
        if (n == 0 || Qlines == NULL || Px2 == NULL || Qlines[0] == NULL) {
            return;
        }

        start_lines_n(out, Qlines, Px2, n);             /* 0x2                  */
        post_add_n_dbl_n(out, Qlines, Px2, n, 1, 2);    /* ..0xc                */
        post_add_n_dbl_n(out, Qlines, Px2, n, 4, 3);    /* ..0x68               */
        post_add_n_dbl_n(out, Qlines, Px2, n, 8, 9);    /* ..0xd200             */
        post_add_n_dbl_n(out, Qlines, Px2, n, 18, 32);  /* ..0xd20100000000     */
        post_add_n_dbl_n(out, Qlines, Px2, n, 51, 16);  /* ..0xd201000000010000 */
        Fp12::conjugate(out.data);
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
            Q ? Q : &G2AffinePoint::generator(),
            P ? P : &G1AffinePoint::generator(),
            1);
    }


    BLST_HD void GT::miller_loop_n(
        GT& out,
        const G2AffinePoint* const Qs[],
        const G1AffinePoint* const Ps[],
        size_t n)
    {   /* ~10KB of stack storage */
        vec_copy(&out, &get_BLST_Rx(), sizeof(Fp12));
        if (n == 0) {
            return;
        }
        if (Qs == NULL || Ps == NULL || Qs[0] == NULL || Ps[0] == NULL) {
            return;
        }

        G2Point T[MILLER_LOOP_N_MAX];
        G2AffinePoint Q[MILLER_LOOP_N_MAX];
        G1AffinePoint Px2[MILLER_LOOP_N_MAX];
        const G2AffinePoint* Qptr = NULL;
        const G1AffinePoint* Pptr = NULL;
        size_t i, j;

        for (i = 0, j = 0; j < n; j++) {
            Qptr = *Qs ? *Qs++ : Qptr + 1;
            Pptr = *Ps ? *Ps++ : Pptr + 1;

            if (!G2AffinePoint::is_inf(*Qptr) && !G1AffinePoint::is_inf(*Pptr)) {
                /* Move common expression from line evaluation to line_by_Px2.  */
                Fp::add(Px2[i].x, Pptr->x, Pptr->x);
                Fp::neg(Px2[i].x, Px2[i].x);
                Fp::add(Px2[i].y, Pptr->y, Pptr->y);

                vec_copy(&Q[i].x, &Qptr->x, 2 * sizeof(Q[i].x));
                vec_copy(&T[i].x, &Qptr->x, 2 * sizeof(T[i].x));
                vec_copy(&T[i].z, &Fp2::one(), sizeof(T[i].z));
                ++i;
            }

            if (i == MILLER_LOOP_N_MAX || (j == n - 1 && i != 0)) {
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

}
