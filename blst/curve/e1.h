/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BLST_E1_H
#define BLST_E1_H

#include "blst/curve/point.h"
#include "blst/fields/Fr.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/common/errors.h"
#include "blst/fields/consts.h"

#include "ec_mult.h"
#include "ec_ops.h"


#if defined(DEBUG) || !defined(NDEBUG)
#include <stdio.h>
#endif

namespace blst
{

    //BLST_EXTERN_C_BEGIN

    BLST_HD BLST_FORCEINLINE void mul_by_b(Fp& out, const Fp& in)
    {
        Fp::lshift(out, in, 2);
    }

    BLST_HD BLST_FORCEINLINE void mul_by_4b(Fp& out, const Fp& in)
    {
        Fp::lshift(out, in, 4);
    }

    BLST_HD BLST_FORCEINLINE void POINTonE1_cneg(POINTonE1& p, bool_t cbit)
    {
        Fp::cneg(p.Y, p.Y, cbit);
    }


    template<>
    BLST_HD BLST_FORCEINLINE void Point<Fp>::cneg(POINTonE1& p, bool_t cbit)
    {
		POINTonE1_cneg(p, cbit);
    }

    BLST_HD BLST_FORCEINLINE void blst_p1_cneg(POINTonE1& a, int cbit)
    {
        POINTonE1_cneg(a, is_zero(cbit) ^ 1);
    }

    BLST_HD BLST_INLINE void POINTonE1_from_Jacobian(POINTonE1& out, const POINTonE1& in)
    {
        Fp Z, ZZ;
        limb_t inf = vec_is_zero(in.Z, sizeof(in.Z));

        Fp::reciprocal(Z, in.Z);                            /* 1/Z   */

        Fp::sqr(ZZ, Z);
        Fp::mul(out.X, in.X, ZZ);                          /* X = X/Z^2 */

        Fp::mul(ZZ, ZZ, Z);
        Fp::mul(out.Y, in.Y, ZZ);                          /* Y = Y/Z^3 */

        vec_select(out.Z, in.Z, get_BLST_G1().Z,
            sizeof(get_BLST_G1().Z), inf);    /* Z = inf ? 0 : 1 */
    }

    BLST_HD BLST_FORCEINLINE void blst_p1_from_jacobian(POINTonE1& out, const POINTonE1& a)
    {
        POINTonE1_from_Jacobian(out, a);
    }

    BLST_HD BLST_INLINE void POINTonE1_to_affine(POINTonE1_affine& out, const POINTonE1& in)
    {
        POINTonE1 p;

		auto ptr = &in;
        if (!vec_is_equal(ptr->Z, get_BLST_Rx().p, sizeof(ptr->Z))) {
            POINTonE1_from_Jacobian(p, *ptr);
            ptr = &p;
        }
        vec_copy(&out, ptr, sizeof(out));
    }

    BLST_HD BLST_FORCEINLINE void blst_p1_to_affine(POINTonE1_affine& out, const POINTonE1& a)
    {
        POINTonE1_to_affine(out, a);
    }

    BLST_HD BLST_INLINE void blst_p1_from_affine(POINTonE1& out, const POINTonE1_affine& a)
    {
        vec_copy(&out, &a, sizeof(a));
        vec_select(&out.Z, &a.X, &get_BLST_Rx().p, sizeof(out.Z),
            vec_is_zero(&a, sizeof(a)));
    }

    BLST_HD BLST_INLINE bool_t POINTonE1_affine_on_curve(const POINTonE1_affine& p)
    {
        Fp XXX, YY;

        Fp::sqr(XXX, p.X);
        Fp::mul(XXX, XXX, p.X);                             /* X^3 */
        Fp::add(XXX, XXX, get_BLST_B_E1());                             /* X^3 + B */

        Fp::sqr(YY, p.Y);                                   /* Y^2 */

        return vec_is_equal(XXX, YY, sizeof(XXX));
    }

    BLST_HD BLST_FORCEINLINE int blst_p1_affine_on_curve(const POINTonE1_affine& p)
    {
        return (int)(POINTonE1_affine_on_curve(p) | vec_is_zero(&p, sizeof(p)));
    }

    BLST_HD BLST_INLINE  bool_t POINTonE1_on_curve(const POINTonE1& p)
    {
        Fp XXX, YY, BZ6;
        limb_t inf = vec_is_zero(&p.Z, sizeof(p.Z));

        Fp::sqr(BZ6, p.Z);
        Fp::mul(BZ6, BZ6, p.Z);
        Fp::sqr(BZ6, BZ6);                                   /* Z^6 */
        mul_by_b(BZ6, BZ6);                            /* B*Z^6 */

        Fp::sqr(XXX, p.X);
        Fp::mul(XXX, XXX, p.X);                             /* X^3 */
        Fp::add(XXX, XXX, BZ6);                              /* X^3 + B*Z^6 */

        Fp::sqr(YY, p.Y);                                   /* Y^2 */

        return vec_is_equal(XXX, YY, sizeof(XXX)) | inf;
    }

    BLST_HD BLST_FORCEINLINE int blst_p1_on_curve(const POINTonE1& p)
    {
        return (int)POINTonE1_on_curve(p);
    }

    BLST_HD BLST_INLINE limb_t POINTonE1_affine_Serialize_BE(unsigned char out[96],
        const POINTonE1_affine& in)
    {
        Fp temp;

        Fp::from_mont(temp, in.X);
        be_bytes_from_limbs(out, temp, sizeof(temp));

        Fp::from_mont(temp, in.Y);
        be_bytes_from_limbs(out + 48, temp, sizeof(temp));

        return Fp::sgn0_pty(temp);
    }

    BLST_HD BLST_INLINE void blst_p1_affine_serialize(unsigned char out[96],
        const POINTonE1_affine& in)
    {
        if (vec_is_zero(in.X, 2 * sizeof(in.X))) {
            bytes_zero(out, 96);
            out[0] = 0x40;    /* infinity bit */
        }
        else {
            (void)POINTonE1_affine_Serialize_BE(out, in);
        }
    }

    BLST_HD BLST_INLINE limb_t POINTonE1_Serialize_BE(unsigned char out[96],
        const POINTonE1& in)
    {
        POINTonE1 p;
		auto ptr = &in;
        if (!vec_is_equal(ptr->Z, get_BLST_Rx().p, sizeof(ptr->Z))) {
            POINTonE1_from_Jacobian(p, *ptr);
            ptr = &p;
        }

        return POINTonE1_affine_Serialize_BE(out, (const POINTonE1_affine&)*ptr);
    }

    BLST_HD BLST_INLINE void POINTonE1_Serialize(unsigned char out[96], const POINTonE1& in)
    {
        if (vec_is_zero(in.Z, sizeof(in.Z))) {
            bytes_zero(out, 96);
            out[0] = 0x40;    /* infinity bit */
        }
        else {
            (void)POINTonE1_Serialize_BE(out, in);
        }
    }

    BLST_HD BLST_FORCEINLINE void blst_p1_serialize(unsigned char out[96], const POINTonE1& in)
    {
        POINTonE1_Serialize(out, in);
    }


    BLST_HD BLST_INLINE limb_t POINTonE1_affine_Compress_BE(unsigned char out[48],
        const POINTonE1_affine& in)
    {
        Fp temp;

        Fp::from_mont(temp, in.X);
        be_bytes_from_limbs(out, temp, sizeof(temp));

        return Fp::sgn0_pty_mont(in.Y);
    }

    BLST_HD BLST_INLINE void blst_p1_affine_compress(unsigned char out[48], const POINTonE1_affine& in)
    {
        if (vec_is_zero(in.X, 2 * sizeof(in.X))) {
            bytes_zero(out, 48);
            out[0] = 0xc0;    /* compressed and infinity bits */
        }
        else {
            limb_t sign = POINTonE1_affine_Compress_BE(out, in);
            out[0] |= (unsigned char)(0x80 | ((sign & 2) << 4));
        }
    }

    BLST_HD BLST_INLINE limb_t POINTonE1_Compress_BE(unsigned char out[48],
        const POINTonE1& in)
    {
        POINTonE1 p;
		auto ptr = &in;

        if (!vec_is_equal(ptr->Z, get_BLST_Rx().p, sizeof(in.Z))) {
            POINTonE1_from_Jacobian(p, *ptr);
            ptr = &p;
        }

        return POINTonE1_affine_Compress_BE(out, (const POINTonE1_affine&)*ptr);
    }

    BLST_HD BLST_INLINE void blst_p1_compress(unsigned char out[48], const POINTonE1& in)
    {
        if (vec_is_zero(in.Z, sizeof(in.Z))) {
            bytes_zero(out, 48);
            out[0] = 0xc0;    /* compressed and infinity bits */
        }
        else {
            limb_t sign = POINTonE1_Compress_BE(out, in);
            out[0] |= (unsigned char)(0x80 | ((sign & 2) << 4));
        }
    }


    BLST_HD BLST_INLINE limb_t POINTonE1_Uncompress_BE(POINTonE1_affine& out,
        const unsigned char in[48])
    {
        POINTonE1_affine ret;
        Fp temp;

        limbs_from_be_bytes(ret.X, in, sizeof(ret.X));
        /* clear top 3 bits in case caller was conveying some information there */
        ret.X[sizeof(ret.X) / sizeof(limb_t) - 1] &= ((limb_t)0 - 1) >> 3;
        Fp::add(temp, ret.X, Fp::zero());  /* less than modulus? */
        if (!vec_is_equal(temp, ret.X, sizeof(temp)))
            return (limb_t)0 - BLST_BAD_ENCODING;
        Fp::mul(ret.X, ret.X, get_BLST_RR());

        Fp::sqr(ret.Y, ret.X);
        Fp::mul(ret.Y, ret.Y, ret.X);
        Fp::add(ret.Y, ret.Y, get_BLST_B_E1());                         /* X^3 + B */
        if (!Fp::sqrt(ret.Y, ret.Y))
            return (limb_t)0 - BLST_POINT_NOT_ON_CURVE;

        vec_copy(&out, &ret, sizeof(ret));

        return Fp::sgn0_pty_mont(out.Y);
    }

    BLST_HD BLST_INLINE BLST_ERROR POINTonE1_Uncompress_Z(POINTonE1_affine& out,
        const unsigned char in[48])
    {
        unsigned char in0 = in[0];
        limb_t sgn0_pty;

        if ((in0 & 0x80) == 0)      /* compressed bit */
            return BLST_BAD_ENCODING;

        if (in0 & 0x40) {           /* infinity bit */
            if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in + 1, 47)) {
                vec_zero(&out, sizeof(out));
                return BLST_SUCCESS;
            }
            else {
                return BLST_BAD_ENCODING;
            }
        }

        sgn0_pty = POINTonE1_Uncompress_BE(out, in);

        if (sgn0_pty > 3)
            return (BLST_ERROR)(0 - sgn0_pty); /* POINT_NOT_ON_CURVE */

        sgn0_pty >>= 1; /* skip over parity bit */
        sgn0_pty ^= (in0 & 0x20) >> 5;
        Fp::cneg(out.Y, out.Y, sgn0_pty);

        /* (0,±2) is not in group, but application might want to ignore? */
        return vec_is_zero(out.X, sizeof(out.X)) ? BLST_POINT_NOT_IN_GROUP
            : BLST_SUCCESS;
    }

    BLST_HD BLST_FORCEINLINE BLST_ERROR blst_p1_uncompress(POINTonE1_affine& out, const unsigned char in[48])
    {
        return POINTonE1_Uncompress_Z(out, in);
    }

    BLST_HD BLST_INLINE BLST_ERROR POINTonE1_Deserialize_BE(POINTonE1_affine& out,
        const unsigned char in[96])
    {
        POINTonE1_affine ret;
        Fp temp;

        limbs_from_be_bytes(ret.X, in, sizeof(ret.X));
        limbs_from_be_bytes(ret.Y, in + 48, sizeof(ret.Y));

        /* clear top 3 bits in case caller was conveying some information there */
        ret.X[sizeof(ret.X) / sizeof(limb_t) - 1] &= ((limb_t)0 - 1) >> 3;
        Fp::add(temp, ret.X, Fp::zero());  /* less than modulus? */
        if (!vec_is_equal(temp, ret.X, sizeof(temp)))
            return BLST_BAD_ENCODING;

        Fp::add(temp, ret.Y, Fp::zero());  /* less than modulus? */
        if (!vec_is_equal(temp, ret.Y, sizeof(temp)))
            return BLST_BAD_ENCODING;

        Fp::mul(ret.X, ret.X, get_BLST_RR());
        Fp::mul(ret.Y, ret.Y, get_BLST_RR());

        if (!POINTonE1_affine_on_curve(ret))
            return BLST_POINT_NOT_ON_CURVE;

        vec_copy(&out, &ret, sizeof(ret));

        /* (0,±2) is not in group, but application might want to ignore? */
        return vec_is_zero(out.X, sizeof(out.X)) ? BLST_POINT_NOT_IN_GROUP
            : BLST_SUCCESS;
    }

    BLST_HD BLST_INLINE BLST_ERROR POINTonE1_Deserialize_Z(POINTonE1_affine& out,
        const unsigned char in[96])
    {
        unsigned char in0 = in[0];

        if ((in0 & 0xe0) == 0)
            return POINTonE1_Deserialize_BE(out, in);

        if (in0 & 0x80)             /* compressed bit */
            return POINTonE1_Uncompress_Z(out, in);

        if (in0 & 0x40) {           /* infinity bit */
            if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in + 1, 95)) {
                vec_zero(&out, sizeof(out));
                return BLST_SUCCESS;
            }
        }

        return BLST_BAD_ENCODING;
    }

    BLST_HD BLST_FORCEINLINE BLST_ERROR blst_p1_deserialize(POINTonE1_affine& out,
        const unsigned char in[96])
    {
        return POINTonE1_Deserialize_Z(out, in);
    }

    //POINT_DADD_IMPL(POINTonE1, 384, fp)
    //    POINT_DADD_AFFINE_IMPL_A0(POINTonE1, 384, fp, get_BLST_Rx().p)
    //    POINT_ADD_IMPL(POINTonE1, 384, fp)
    //    POINT_ADD_AFFINE_IMPL(POINTonE1, 384, fp, get_BLST_Rx().p)
    //    POINT_DOUBLE_IMPL_A0(POINTonE1, 384, fp)
    //    POINT_IS_EQUAL_IMPL(POINTonE1, 384, fp)

        BLST_HD BLST_INLINE void blst_p1_add(POINTonE1& out, const POINTonE1& a, const POINTonE1& b)
    {
#if (defined(DEBUG) || !defined(NDEBUG)) && !defined(__CUDA_ARCH__)
        // Check if we're trying to add a point to itself in debug builds
        if (!vec_is_zero(a.Z, sizeof(a.Z)) && !vec_is_zero(b.Z, sizeof(b.Z))) { // If neither point is infinity
            // Convert both points to affine to check if they're the same
            Fp Z_a, ZZ_a, Z_b, ZZ_b;

            // Convert a to affine
            Fp::reciprocal(Z_a, a.Z);            // 1/Z_a
            Fp::sqr(ZZ_a, Z_a);                   // 1/Z_a^2
            Fp a_X_affine, a_Y_affine;
            Fp::mul(a_X_affine, a.X, ZZ_a);      // X_a/Z_a^2

            // Convert b to affine
            Fp::reciprocal(Z_b, b.Z);            // 1/Z_b
            Fp::sqr(ZZ_b, Z_b);                   // 1/Z_b^2
            Fp b_X_affine, b_Y_affine;
            Fp::mul(b_X_affine, b.X, ZZ_b);      // X_b/Z_b^2

            // Compare X coordinates
            if (vec_is_equal(a_X_affine, b_X_affine, sizeof(Fp))) {
                // Also check Y coordinates
                Fp::mul(ZZ_a, ZZ_a, Z_a);         // 1/Z_a^3
                Fp::mul(a_Y_affine, a.Y, ZZ_a);  // Y_a/Z_a^3

                Fp::mul(ZZ_b, ZZ_b, Z_b);         // 1/Z_b^3
                Fp::mul(b_Y_affine, b.Y, ZZ_b);  // Y_b/Z_b^3

                if (vec_is_equal(a_Y_affine, b_Y_affine, sizeof(Fp))) {
                    printf("ERROR: blst_p1_add called with the same point. "
                        "Use blst_p1_add_or_double instead.\n");
                    exit(1);
                }
            }
        }
#endif

        POINTonE1::add(out, a, b);
    }

    BLST_HD BLST_FORCEINLINE void blst_p1_add_or_double(POINTonE1& out, const POINTonE1& a,
        const POINTonE1& b)
    {
        POINTonE1::dadd(out, a, b, NULL);
    }

    BLST_HD BLST_INLINE void blst_p1_add_affine(POINTonE1& out, const POINTonE1& a,
        const POINTonE1_affine& b)
    {
#if (defined(DEBUG) || !defined(NDEBUG)) && !defined(__CUDA_ARCH__)
        // Check if we're trying to add a point to itself in debug builds
        if (!vec_is_zero(a.Z, sizeof(a.Z))) {  // If a is not infinity
            Fp Z, ZZ;
            Fp::reciprocal(Z, a.Z);              // 1/Z
            Fp::sqr(ZZ, Z);                       // 1/Z^2

            // Compare X coordinates: a.X/Z^2 == b.X?
            Fp a_X_affine;
            Fp::mul(a_X_affine, a.X, ZZ);        // X/Z^2

            if (vec_is_equal(a_X_affine, b.X, sizeof(Fp))) {
                // Check Y coordinates as well
                Fp a_Y_affine;
                Fp::mul(ZZ, ZZ, Z);               // 1/Z^3
                Fp::mul(a_Y_affine, a.Y, ZZ);    // Y/Z^3

                if (vec_is_equal(a_Y_affine, b.Y, sizeof(Fp))) {
                    printf("ERROR: blst_p1_add_affine called with the same point. "
                        "Use blst_p1_add_or_double_affine instead.\n");
                    exit(1);
                    // We could assert here, but let's just warn and proceed
                    // assert(0 && "blst_p1_add_affine called with same point");
                }
            }
        }
#endif

        POINTonE1::add_affine(out, a, b);
    }

    BLST_HD BLST_FORCEINLINE void blst_p1_add_or_double_affine(POINTonE1& out, const POINTonE1& a,
        const POINTonE1_affine& b)
    {
        POINTonE1::dadd_affine(out, a, b);
    }

    BLST_HD BLST_FORCEINLINE void blst_p1_double(POINTonE1& out, const POINTonE1& a)
    {
        POINTonE1::doubling(out, a);
    }

    BLST_HD BLST_FORCEINLINE int blst_p1_is_equal(const POINTonE1& a, const POINTonE1& b)
    {
        return (int)POINTonE1::is_equal(a, b);
    }

    //POINT_MULT_SCALAR_WX_IMPL(POINTonE1, 4)
    //    POINT_MULT_SCALAR_WX_IMPL(POINTonE1, 5)

#ifdef __BLST_PRIVATE_TESTMODE__
        POINT_AFFINE_MULT_SCALAR_IMPL(POINTonE1)

        DECLARE_PRIVATE_POINTXZ(POINTonE1, 384)
    POINT_LADDER_PRE_IMPL(POINTonE1, 384, fp)
        POINT_LADDER_STEP_IMPL_A0(POINTonE1, 384, fp, onE1)
        POINT_LADDER_POST_IMPL_A0(POINTonE1, 384, fp, onE1)
        POINT_MULT_SCALAR_LADDER_IMPL(POINTonE1)
#endif


        BLST_HD BLST_FORCEINLINE void sigma(POINTonE1& out, const POINTonE1& in)
    {
        vec_copy(out.X, in.X, 2 * sizeof(out.X));
        Fp::mul(out.Z, in.Z, get_BLST_beta());
    }

    /* Gallant-Lambert-Vanstone, ~45% faster than mult_w<POINTonE1,5> */
    BLST_HD BLST_INLINE void POINTonE1_mult_glv(POINTonE1& out, const POINTonE1& in,
        const pow256 SK)
    {
        union { Fr l; pow256 s; } val;

        /* SK/z^2 [in constant time] */

        limbs_from_le_bytes(val.l, SK, 32);
        div_by_zz(val.l);
        le_bytes_from_limbs(val.s, val.l, 32);

        {
            const byte* scalars[2] = { val.s + 16, val.s };
            POINTonE1 table[2][1 << (5 - 1)];   /* 4.5KB */
            size_t i;

            precompute_w<POINTonE1, 5>(table[0], in);
            for (i = 0; i < 1 << (5 - 1); i++) {
                Fp::mul(table[1][i].X, table[0][i].X, get_BLST_beta());
                Fp::cneg(table[1][i].Y, table[0][i].Y, 1);
                vec_copy(table[1][i].Z, table[0][i].Z, sizeof(table[1][i].Z));
            }

            s_mult_w<POINTonE1, 5>(out, NULL, 2, scalars, 128, table);
            POINTonE1_cneg(out, 1);
            Fp::mul(out.Z, out.Z, get_BLST_beta());
            Fp::mul(out.Z, out.Z, get_BLST_beta());
        }

        vec_zero(val.l, sizeof(val));   /* scrub the copy of SK */
    }

    BLST_HD BLST_INLINE void POINTonE1_sign(POINTonE1& out, const POINTonE1& in, const pow256 SK)
    {
        Fp Z, ZZ;
        limb_t inf;

        POINTonE1_mult_glv(out, in, SK);

        /* convert to affine to remove possible bias in out.Z */
        inf = vec_is_zero(out.Z, sizeof(out.Z));
#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
        Fp::flt_reciprocal(Z, out.Z);                       /* 1/Z   */
#else
        Fp::reciprocal(Z, out.Z);                           /* 1/Z   */
#endif

        Fp::sqr(ZZ, Z);
        Fp::mul(out.X, out.X, ZZ);                         /* X = X/Z^2 */

        Fp::mul(ZZ, ZZ, Z);
        Fp::mul(out.Y, out.Y, ZZ);                         /* Y = Y/Z^3 */

        vec_select(out.Z, out.Z, get_BLST_G1().Z, sizeof(get_BLST_G1().Z),
            inf);                            /* Z = inf ? 0 : 1 */
    }

    BLST_HD BLST_FORCEINLINE void blst_sk_to_pk_in_g1(POINTonE1& out, const pow256 SK)
    {
        POINTonE1_sign(out, get_BLST_G1(), SK);
    }

    BLST_HD BLST_FORCEINLINE void blst_sign_pk_in_g2(POINTonE1& out, const POINTonE1& msg, const pow256 SK)
    {
        POINTonE1_sign(out, msg, SK);
    }

    BLST_HD BLST_INLINE void blst_sk_to_pk2_in_g1(unsigned char out[96], POINTonE1_affine* PK,
        const pow256 SK)
    {
        POINTonE1 P;

        POINTonE1_sign(P, get_BLST_G1(), SK);
        if (PK != NULL)
            vec_copy(PK, &P, sizeof(*PK));
        if (out != NULL) {
            limb_t sgn0_pty = POINTonE1_Serialize_BE(out, P);
            out[0] |= (sgn0_pty & 2) << 4;      /* pre-decorate */
            out[0] |= vec_is_zero(P.Z, sizeof(P.Z)) << 6;
        }
    }

    BLST_HD BLST_INLINE void blst_sign_pk2_in_g2(unsigned char out[96], POINTonE1_affine* sig,
        const POINTonE1& hash, const pow256 SK)
    {
        POINTonE1 P;

        POINTonE1_sign(P, hash, SK);
        if (sig != NULL)
            vec_copy(sig, &P, sizeof(*sig));
        if (out != NULL) {
            limb_t sgn0_pty = POINTonE1_Serialize_BE(out, P);
            out[0] |= (sgn0_pty & 2) << 4;      /* pre-decorate */
            out[0] |= vec_is_zero(P.Z, sizeof(P.Z)) << 6;
        }
    }

    BLST_HD BLST_INLINE void blst_p1_mult(POINTonE1& out, const POINTonE1& a,
        const byte* scalar, size_t nbits)
    {
        if (nbits < 176) {
            if (nbits)
                mult_w<POINTonE1,4>(out, a, scalar, nbits);
            else
                vec_zero(&out, sizeof(out));
        }
        else if (nbits <= 256) {
            union { Fr l; pow256 s; } val;
            size_t i, j, top, mask = (size_t)0 - 1;

            /* this is not about constant-time-ness, but branch optimization */
            for (top = (nbits + 7) / 8, i = 0, j = 0; i < sizeof(val.s);) {
                val.s[i++] = scalar[j] & mask;
                mask = 0 - ((i - top) >> (8 * sizeof(top) - 1));
                j += 1 & mask;
            }

            if (check_mod_v256(val.s, get_BLST_r()))  /* z^4 is the formal limit */
                POINTonE1_mult_glv(out, a, val.s);
            else    /* should never be the case, added for formal completeness */
                mult_w<POINTonE1,5>(out, a, scalar, nbits);

            vec_zero(val.l, sizeof(val));
        }
        else {    /* should never be the case, added for formal completeness */
            mult_w<POINTonE1,5>(out, a, scalar, nbits);
        }
    }

    BLST_HD BLST_FORCEINLINE void blst_p1_unchecked_mult(POINTonE1& out, const POINTonE1& a,
        const byte* scalar, size_t nbits)
    {
        if (nbits)
            mult_w<POINTonE1,4>(out, a, scalar, nbits);
        else
            vec_zero(&out, sizeof(out));
    }

    BLST_HD BLST_FORCEINLINE int blst_p1_affine_is_equal(const POINTonE1_affine& a,
        const POINTonE1_affine& b)
    {
        return (int)vec_is_equal(&a, &b, sizeof(a));
    }

    BLST_HD BLST_FORCEINLINE int blst_p1_is_inf(const POINTonE1& p)
    {
        return (int)vec_is_zero(p.Z, sizeof(p.Z));
    }

    BLST_HD BLST_FORCEINLINE const POINTonE1& blst_p1_generator(void)
    {
        return get_BLST_G1();
    }

    BLST_HD BLST_FORCEINLINE int blst_p1_affine_is_inf(const POINTonE1_affine& p)
    {
        return (int)vec_is_zero(&p, sizeof(p));
    }

    BLST_HD BLST_FORCEINLINE const POINTonE1_affine& blst_p1_affine_generator(void)
    {
        return (const POINTonE1_affine&)get_BLST_G1();
    }

    BLST_HD BLST_FORCEINLINE size_t blst_p1_sizeof(void)
    {
        return sizeof(POINTonE1);
    }

    BLST_HD BLST_FORCEINLINE size_t blst_p1_affine_sizeof(void)
    {
        return sizeof(POINTonE1_affine);
    }

    //BLST_EXTERN_C_END
}

#endif /* BLST_E1_H */