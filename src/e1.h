/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "point.h"
#include "fields.h"
#include "errors.h"


#if defined(DEBUG) || !defined(NDEBUG)
#include <stdio.h>
#endif

/*
 * y^2 = x^3 + B
 */
const vec384 B_E1 = {        /* (4 << 384) % P */
    TO_LIMB_T(0xaa270000000cfff3), TO_LIMB_T(0x53cc0032fc34000a),
    TO_LIMB_T(0x478fe97a6b0a807f), TO_LIMB_T(0xb1d37ebee6ba24d7),
    TO_LIMB_T(0x8ec9733bbf78ab2f), TO_LIMB_T(0x09d645513d83de7e)
};

const POINTonE1 BLS12_381_G1 = {    /* generator point [in Montgomery] */
  /* (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905
   *    a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb << 384) % P */
  { TO_LIMB_T(0x5cb38790fd530c16), TO_LIMB_T(0x7817fc679976fff5),
    TO_LIMB_T(0x154f95c7143ba1c1), TO_LIMB_T(0xf0ae6acdf3d0e747),
    TO_LIMB_T(0xedce6ecc21dbf440), TO_LIMB_T(0x120177419e0bfb75) },
  /* (0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af6
   *    00db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1 << 384) % P */
  { TO_LIMB_T(0xbaac93d50ce72271), TO_LIMB_T(0x8c22631a7918fd8e),
    TO_LIMB_T(0xdd595f13570725ce), TO_LIMB_T(0x51ac582950405194),
    TO_LIMB_T(0x0e1c8c3fad0059c0), TO_LIMB_T(0x0bbc3efc5008a26a) },
  { ONE_MONT_P }
};

const POINTonE1 BLS12_381_NEG_G1 = { /* negative generator [in Montgomery] */
  /* (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905
   *    a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb << 384) % P */
  { TO_LIMB_T(0x5cb38790fd530c16), TO_LIMB_T(0x7817fc679976fff5),
    TO_LIMB_T(0x154f95c7143ba1c1), TO_LIMB_T(0xf0ae6acdf3d0e747),
    TO_LIMB_T(0xedce6ecc21dbf440), TO_LIMB_T(0x120177419e0bfb75) },
  /* (0x114d1d6855d545a8aa7d76c8cf2e21f267816aef1db507c9
   *    6655b9d5caac42364e6f38ba0ecb751bad54dcd6b939c2ca << 384) % P */
  { TO_LIMB_T(0xff526c2af318883a), TO_LIMB_T(0x92899ce4383b0270),
    TO_LIMB_T(0x89d7738d9fa9d055), TO_LIMB_T(0x12caf35ba344c12a),
    TO_LIMB_T(0x3cff1b76964b5317), TO_LIMB_T(0x0e44d2ede9774430) },
  { ONE_MONT_P }
};

inline void mul_by_b_onE1(vec384 out, const vec384 in)
{   lshift_fp(out, in, 2);   }

inline void mul_by_4b_onE1(vec384 out, const vec384 in)
{   lshift_fp(out, in, 4);   }

inline void POINTonE1_cneg(POINTonE1 *p, bool_t cbit)
{   cneg_fp(p->Y, p->Y, cbit);   }

inline void blst_p1_cneg(POINTonE1 *a, int cbit)
{   POINTonE1_cneg(a, is_zero(cbit) ^ 1);   }

inline void POINTonE1_from_Jacobian(POINTonE1 *out, const POINTonE1 *in)
{
    vec384 Z, ZZ;
    limb_t inf = vec_is_zero(in->Z, sizeof(in->Z));

    reciprocal_fp(Z, in->Z);                            /* 1/Z   */

    sqr_fp(ZZ, Z);
    mul_fp(out->X, in->X, ZZ);                          /* X = X/Z^2 */

    mul_fp(ZZ, ZZ, Z);
    mul_fp(out->Y, in->Y, ZZ);                          /* Y = Y/Z^3 */

    vec_select(out->Z, in->Z, BLS12_381_G1.Z,
                       sizeof(BLS12_381_G1.Z), inf);    /* Z = inf ? 0 : 1 */
}

inline void blst_p1_from_jacobian(POINTonE1 *out, const POINTonE1 *a)
{   POINTonE1_from_Jacobian(out, a);   }

inline void POINTonE1_to_affine(POINTonE1_affine *out, const POINTonE1 *in)
{
    POINTonE1 p;

    if (!vec_is_equal(in->Z, get_BLST_Rx()->p, sizeof(in->Z))) {
        POINTonE1_from_Jacobian(&p, in);
        in = &p;
    }
    vec_copy(out, in, sizeof(*out));
}

inline void blst_p1_to_affine(POINTonE1_affine *out, const POINTonE1 *a)
{   POINTonE1_to_affine(out, a);   }

inline void blst_p1_from_affine(POINTonE1 *out, const POINTonE1_affine *a)
{
    vec_copy(out, a, sizeof(*a));
    vec_select(out->Z, a->X, get_BLST_Rx()->p, sizeof(out->Z),
                       vec_is_zero(a, sizeof(*a)));
}

inline bool_t POINTonE1_affine_on_curve(const POINTonE1_affine *p)
{
    vec384 XXX, YY;

    sqr_fp(XXX, p->X);
    mul_fp(XXX, XXX, p->X);                             /* X^3 */
    add_fp(XXX, XXX, B_E1);                             /* X^3 + B */

    sqr_fp(YY, p->Y);                                   /* Y^2 */

    return vec_is_equal(XXX, YY, sizeof(XXX));
}

inline int blst_p1_affine_on_curve(const POINTonE1_affine *p)
{   return (int)(POINTonE1_affine_on_curve(p) | vec_is_zero(p, sizeof(*p)));   }

inline  bool_t POINTonE1_on_curve(const POINTonE1 *p)
{
    vec384 XXX, YY, BZ6;
    limb_t inf = vec_is_zero(p->Z, sizeof(p->Z));

    sqr_fp(BZ6, p->Z);
    mul_fp(BZ6, BZ6, p->Z);
    sqr_fp(BZ6, BZ6);                                   /* Z^6 */
    mul_by_b_onE1(BZ6, BZ6);                            /* B*Z^6 */

    sqr_fp(XXX, p->X);
    mul_fp(XXX, XXX, p->X);                             /* X^3 */
    add_fp(XXX, XXX, BZ6);                              /* X^3 + B*Z^6 */

    sqr_fp(YY, p->Y);                                   /* Y^2 */

    return vec_is_equal(XXX, YY, sizeof(XXX)) | inf;
}

inline int blst_p1_on_curve(const POINTonE1 *p)
{   return (int)POINTonE1_on_curve(p);   }

inline limb_t POINTonE1_affine_Serialize_BE(unsigned char out[96],
                                            const POINTonE1_affine *in)
{
    vec384 temp;

    from_fp(temp, in->X);
    be_bytes_from_limbs(out, temp, sizeof(temp));

    from_fp(temp, in->Y);
    be_bytes_from_limbs(out + 48, temp, sizeof(temp));

    return sgn0_pty_mod_v384(temp, *get_BLST_P());
}

inline void blst_p1_affine_serialize(unsigned char out[96],
                              const POINTonE1_affine *in)
{
    if (vec_is_zero(in->X, 2*sizeof(in->X))) {
        bytes_zero(out, 96);
        out[0] = 0x40;    /* infinity bit */
    } else {
        (void)POINTonE1_affine_Serialize_BE(out, in);
    }
}

inline limb_t POINTonE1_Serialize_BE(unsigned char out[96],
                                     const POINTonE1 *in)
{
    POINTonE1 p;

    if (!vec_is_equal(in->Z, get_BLST_Rx()->p, sizeof(in->Z))) {
        POINTonE1_from_Jacobian(&p, in);
        in = &p;
    }

    return POINTonE1_affine_Serialize_BE(out, (const POINTonE1_affine *)in);
}

inline void POINTonE1_Serialize(unsigned char out[96], const POINTonE1 *in)
{
    if (vec_is_zero(in->Z, sizeof(in->Z))) {
        bytes_zero(out, 96);
        out[0] = 0x40;    /* infinity bit */
    } else {
        (void)POINTonE1_Serialize_BE(out, in);
    }
}

inline void blst_p1_serialize(unsigned char out[96], const POINTonE1 *in)
{   POINTonE1_Serialize(out, in);   }


inline limb_t POINTonE1_affine_Compress_BE(unsigned char out[48],
                                           const POINTonE1_affine *in)
{
    vec384 temp;

    from_fp(temp, in->X);
    be_bytes_from_limbs(out, temp, sizeof(temp));

    return sgn0_pty_mont_v384(in->Y, *get_BLST_P(), BLST_p0);
}

inline void blst_p1_affine_compress(unsigned char out[48], const POINTonE1_affine *in)
{
    if (vec_is_zero(in->X, 2*sizeof(in->X))) {
        bytes_zero(out, 48);
        out[0] = 0xc0;    /* compressed and infinity bits */
    } else {
        limb_t sign = POINTonE1_affine_Compress_BE(out, in);
        out[0] |= (unsigned char)(0x80 | ((sign & 2) << 4));
    }
}

inline limb_t POINTonE1_Compress_BE(unsigned char out[48],
                                    const POINTonE1 *in)
{
    POINTonE1 p;

    if (!vec_is_equal(in->Z, get_BLST_Rx()->p, sizeof(in->Z))) {
        POINTonE1_from_Jacobian(&p, in);
        in = &p;
    }

    return POINTonE1_affine_Compress_BE(out, (const POINTonE1_affine *)in);
}

inline void blst_p1_compress(unsigned char out[48], const POINTonE1 *in)
{
    if (vec_is_zero(in->Z, sizeof(in->Z))) {
        bytes_zero(out, 48);
        out[0] = 0xc0;    /* compressed and infinity bits */
    } else {
        limb_t sign = POINTonE1_Compress_BE(out, in);
        out[0] |= (unsigned char)(0x80 | ((sign & 2) << 4));
    }
}


inline limb_t POINTonE1_Uncompress_BE(POINTonE1_affine *out,
                                      const unsigned char in[48])
{
    POINTonE1_affine ret;
    vec384 temp;

    limbs_from_be_bytes(ret.X, in, sizeof(ret.X));
    /* clear top 3 bits in case caller was conveying some information there */
    ret.X[sizeof(ret.X)/sizeof(limb_t)-1] &= ((limb_t)0-1) >> 3;
    add_fp(temp, ret.X, ZERO_384);  /* less than modulus? */
    if (!vec_is_equal(temp, ret.X, sizeof(temp)))
        return (limb_t)0 - BLST_BAD_ENCODING;
    mul_fp(ret.X, ret.X, *get_BLST_RR());

    sqr_fp(ret.Y, ret.X);
    mul_fp(ret.Y, ret.Y, ret.X);
    add_fp(ret.Y, ret.Y, B_E1);                         /* X^3 + B */
    if (!sqrt_fp(ret.Y, ret.Y))
        return (limb_t)0 - BLST_POINT_NOT_ON_CURVE;

    vec_copy(out, &ret, sizeof(ret));

    return sgn0_pty_mont_v384(out->Y, *get_BLST_P(), BLST_p0);
}

inline BLST_ERROR POINTonE1_Uncompress_Z(POINTonE1_affine *out,
                                         const unsigned char in[48])
{
    unsigned char in0 = in[0];
    limb_t sgn0_pty;

    if ((in0 & 0x80) == 0)      /* compressed bit */
        return BLST_BAD_ENCODING;

    if (in0 & 0x40) {           /* infinity bit */
        if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in+1, 47)) {
            vec_zero(out, sizeof(*out));
            return BLST_SUCCESS;
        } else {
            return BLST_BAD_ENCODING;
        }
    }

    sgn0_pty = POINTonE1_Uncompress_BE(out, in);

    if (sgn0_pty > 3)
        return (BLST_ERROR)(0 - sgn0_pty); /* POINT_NOT_ON_CURVE */

    sgn0_pty >>= 1; /* skip over parity bit */
    sgn0_pty ^= (in0 & 0x20) >> 5;
    cneg_fp(out->Y, out->Y, sgn0_pty);

    /* (0,±2) is not in group, but application might want to ignore? */
    return vec_is_zero(out->X, sizeof(out->X)) ? BLST_POINT_NOT_IN_GROUP
                                               : BLST_SUCCESS;
}

inline BLST_ERROR blst_p1_uncompress(POINTonE1_affine *out, const unsigned char in[48])
{   return POINTonE1_Uncompress_Z(out, in);   }

inline BLST_ERROR POINTonE1_Deserialize_BE(POINTonE1_affine *out,
                                           const unsigned char in[96])
{
    POINTonE1_affine ret;
    vec384 temp;

    limbs_from_be_bytes(ret.X, in, sizeof(ret.X));
    limbs_from_be_bytes(ret.Y, in + 48, sizeof(ret.Y));

    /* clear top 3 bits in case caller was conveying some information there */
    ret.X[sizeof(ret.X)/sizeof(limb_t)-1] &= ((limb_t)0-1) >> 3;
    add_fp(temp, ret.X, ZERO_384);  /* less than modulus? */
    if (!vec_is_equal(temp, ret.X, sizeof(temp)))
        return BLST_BAD_ENCODING;

    add_fp(temp, ret.Y, ZERO_384);  /* less than modulus? */
    if (!vec_is_equal(temp, ret.Y, sizeof(temp)))
        return BLST_BAD_ENCODING;

    mul_fp(ret.X, ret.X, *get_BLST_RR());
    mul_fp(ret.Y, ret.Y, *get_BLST_RR());

    if (!POINTonE1_affine_on_curve(&ret))
        return BLST_POINT_NOT_ON_CURVE;

    vec_copy(out, &ret, sizeof(ret));

    /* (0,±2) is not in group, but application might want to ignore? */
    return vec_is_zero(out->X, sizeof(out->X)) ? BLST_POINT_NOT_IN_GROUP
                                               : BLST_SUCCESS;
}

inline BLST_ERROR POINTonE1_Deserialize_Z(POINTonE1_affine *out,
                                          const unsigned char in[96])
{
    unsigned char in0 = in[0];

    if ((in0 & 0xe0) == 0)
        return POINTonE1_Deserialize_BE(out, in);

    if (in0 & 0x80)             /* compressed bit */
        return POINTonE1_Uncompress_Z(out, in);

    if (in0 & 0x40) {           /* infinity bit */
        if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in+1, 95)) {
            vec_zero(out, sizeof(*out));
            return BLST_SUCCESS;
        }
    }

    return BLST_BAD_ENCODING;
}

inline BLST_ERROR blst_p1_deserialize(POINTonE1_affine *out,
                               const unsigned char in[96])
{   return POINTonE1_Deserialize_Z(out, in);   }

#include "ec_ops.h"
POINT_DADD_IMPL(POINTonE1, 384, fp)
POINT_DADD_AFFINE_IMPL_A0(POINTonE1, 384, fp, get_BLST_Rx()->p)
POINT_ADD_IMPL(POINTonE1, 384, fp)
POINT_ADD_AFFINE_IMPL(POINTonE1, 384, fp, get_BLST_Rx()->p)
POINT_DOUBLE_IMPL_A0(POINTonE1, 384, fp)
POINT_IS_EQUAL_IMPL(POINTonE1, 384, fp)

inline void blst_p1_add(POINTonE1 *out, const POINTonE1 *a, const POINTonE1 *b)
{  
#if defined(DEBUG) || !defined(NDEBUG)
    // Check if we're trying to add a point to itself in debug builds
    if (!vec_is_zero(a->Z, sizeof(a->Z)) && !vec_is_zero(b->Z, sizeof(b->Z))) { // If neither point is infinity
        // Convert both points to affine to check if they're the same
        vec384 Z_a, ZZ_a, Z_b, ZZ_b;

        // Convert a to affine
        reciprocal_fp(Z_a, a->Z);            // 1/Z_a
        sqr_fp(ZZ_a, Z_a);                   // 1/Z_a^2
        vec384 a_X_affine, a_Y_affine;
        mul_fp(a_X_affine, a->X, ZZ_a);      // X_a/Z_a^2

        // Convert b to affine
        reciprocal_fp(Z_b, b->Z);            // 1/Z_b
        sqr_fp(ZZ_b, Z_b);                   // 1/Z_b^2
        vec384 b_X_affine, b_Y_affine;
        mul_fp(b_X_affine, b->X, ZZ_b);      // X_b/Z_b^2

        // Compare X coordinates
        if (vec_is_equal(a_X_affine, b_X_affine, sizeof(vec384))) {
            // Also check Y coordinates
            mul_fp(ZZ_a, ZZ_a, Z_a);         // 1/Z_a^3
            mul_fp(a_Y_affine, a->Y, ZZ_a);  // Y_a/Z_a^3

            mul_fp(ZZ_b, ZZ_b, Z_b);         // 1/Z_b^3
            mul_fp(b_Y_affine, b->Y, ZZ_b);  // Y_b/Z_b^3

            if (vec_is_equal(a_Y_affine, b_Y_affine, sizeof(vec384))) {
                printf("ERROR: blst_p1_add called with the same point. "
                    "Use blst_p1_add_or_double instead.\n");
                exit(1);
            }
        }
    }
#endif

    POINTonE1_add(out, a, b);   
}

inline void blst_p1_add_or_double(POINTonE1 *out, const POINTonE1 *a,
                                           const POINTonE1 *b)
{   POINTonE1_dadd(out, a, b, NULL);   }

inline void blst_p1_add_affine(POINTonE1 *out, const POINTonE1 *a,
                                        const POINTonE1_affine *b)
{   
#if defined(DEBUG) || !defined(NDEBUG)
    // Check if we're trying to add a point to itself in debug builds
    if (!vec_is_zero(a->Z, sizeof(a->Z))) {  // If a is not infinity
        vec384 Z, ZZ;
        reciprocal_fp(Z, a->Z);              // 1/Z
        sqr_fp(ZZ, Z);                       // 1/Z^2

        // Compare X coordinates: a->X/Z^2 == b->X?
        vec384 a_X_affine;
        mul_fp(a_X_affine, a->X, ZZ);        // X/Z^2

        if (vec_is_equal(a_X_affine, b->X, sizeof(vec384))) {
            // Check Y coordinates as well
            vec384 a_Y_affine;
            mul_fp(ZZ, ZZ, Z);               // 1/Z^3
            mul_fp(a_Y_affine, a->Y, ZZ);    // Y/Z^3

            if (vec_is_equal(a_Y_affine, b->Y, sizeof(vec384))) {
                printf("ERROR: blst_p1_add_affine called with the same point. "
                    "Use blst_p1_add_or_double_affine instead.\n");
                exit(1);
                // We could assert here, but let's just warn and proceed
                // assert(0 && "blst_p1_add_affine called with same point");
            }
        }
    }
#endif

    POINTonE1_add_affine(out, a, b);   
}

inline void blst_p1_add_or_double_affine(POINTonE1 *out, const POINTonE1 *a,
                                                  const POINTonE1_affine *b)
{   POINTonE1_dadd_affine(out, a, b);   }

inline void blst_p1_double(POINTonE1 *out, const POINTonE1 *a)
{   POINTonE1_double(out, a);   }

inline int blst_p1_is_equal(const POINTonE1 *a, const POINTonE1 *b)
{   return (int)POINTonE1_is_equal(a, b);   }

#include "ec_mult.h"
POINT_MULT_SCALAR_WX_IMPL(POINTonE1, 4)
POINT_MULT_SCALAR_WX_IMPL(POINTonE1, 5)

#ifdef __BLST_PRIVATE_TESTMODE__
POINT_AFFINE_MULT_SCALAR_IMPL(POINTonE1)

DECLARE_PRIVATE_POINTXZ(POINTonE1, 384)
POINT_LADDER_PRE_IMPL(POINTonE1, 384, fp)
POINT_LADDER_STEP_IMPL_A0(POINTonE1, 384, fp, onE1)
POINT_LADDER_POST_IMPL_A0(POINTonE1, 384, fp, onE1)
POINT_MULT_SCALAR_LADDER_IMPL(POINTonE1)
#endif

const vec384 beta = {            /* such that beta^3 - 1 = 0  */
    /* -1/2 * (1 + sqrt(-3)) = ((P-2)^(P-2)) * (1 + (P-3)^((P+1)/4)) */
    /* (0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4
          897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac << 384) % P */
    TO_LIMB_T(0xcd03c9e48671f071), TO_LIMB_T(0x5dab22461fcda5d2),
    TO_LIMB_T(0x587042afd3851b95), TO_LIMB_T(0x8eb60ebe01bacb9e),
    TO_LIMB_T(0x03f97d6e83d050d2), TO_LIMB_T(0x18f0206554638741)
};

inline void sigma(POINTonE1 *out, const POINTonE1 *in)
{
    vec_copy(out->X, in->X, 2*sizeof(out->X));
    mul_fp(out->Z, in->Z, beta);
}

/* Gallant-Lambert-Vanstone, ~45% faster than POINTonE1_mult_w5 */
inline void POINTonE1_mult_glv(POINTonE1 *out, const POINTonE1 *in,
                               const pow256 SK)
{
    union { vec256 l; pow256 s; } val;

    /* SK/z^2 [in constant time] */

    limbs_from_le_bytes(val.l, SK, 32);
    div_by_zz(val.l);
    le_bytes_from_limbs(val.s, val.l, 32);

    {
        const byte *scalars[2] = { val.s+16, val.s };
        POINTonE1 table[2][1<<(5-1)];   /* 4.5KB */
        size_t i;

        POINTonE1_precompute_w5(table[0], in);
        for (i = 0; i < 1<<(5-1); i++) {
            mul_fp(table[1][i].X, table[0][i].X, beta);
            cneg_fp(table[1][i].Y, table[0][i].Y, 1);
            vec_copy(table[1][i].Z, table[0][i].Z, sizeof(table[1][i].Z));
        }

        POINTonE1s_mult_w5(out, NULL, 2, scalars, 128, table);
        POINTonE1_cneg(out, 1);
        mul_fp(out->Z, out->Z, beta);
        mul_fp(out->Z, out->Z, beta);
    }

    vec_zero(val.l, sizeof(val));   /* scrub the copy of SK */
}

inline void POINTonE1_sign(POINTonE1 *out, const POINTonE1 *in, const pow256 SK)
{
    vec384 Z, ZZ;
    limb_t inf;

    POINTonE1_mult_glv(out, in, SK);

    /* convert to affine to remove possible bias in out->Z */
    inf = vec_is_zero(out->Z, sizeof(out->Z));
#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    flt_reciprocal_fp(Z, out->Z);                       /* 1/Z   */
#else
    reciprocal_fp(Z, out->Z);                           /* 1/Z   */
#endif

    sqr_fp(ZZ, Z);
    mul_fp(out->X, out->X, ZZ);                         /* X = X/Z^2 */

    mul_fp(ZZ, ZZ, Z);
    mul_fp(out->Y, out->Y, ZZ);                         /* Y = Y/Z^3 */

    vec_select(out->Z, out->Z, BLS12_381_G1.Z, sizeof(BLS12_381_G1.Z),
                       inf);                            /* Z = inf ? 0 : 1 */
}

inline void blst_sk_to_pk_in_g1(POINTonE1 *out, const pow256 SK)
{   POINTonE1_sign(out, &BLS12_381_G1, SK);   }

inline void blst_sign_pk_in_g2(POINTonE1 *out, const POINTonE1 *msg, const pow256 SK)
{   POINTonE1_sign(out, msg, SK);   }

inline void blst_sk_to_pk2_in_g1(unsigned char out[96], POINTonE1_affine *PK,
                          const pow256 SK)
{
    POINTonE1 P[1];

    POINTonE1_sign(P, &BLS12_381_G1, SK);
    if (PK != NULL)
        vec_copy(PK, P, sizeof(*PK));
    if (out != NULL) {
        limb_t sgn0_pty = POINTonE1_Serialize_BE(out, P);
        out[0] |= (sgn0_pty & 2) << 4;      /* pre-decorate */
        out[0] |= vec_is_zero(P->Z, sizeof(P->Z)) << 6;
    }
}

inline void blst_sign_pk2_in_g2(unsigned char out[96], POINTonE1_affine *sig,
                         const POINTonE1 *hash, const pow256 SK)
{
    POINTonE1 P[1];

    POINTonE1_sign(P, hash, SK);
    if (sig != NULL)
        vec_copy(sig, P, sizeof(*sig));
    if (out != NULL) {
        limb_t sgn0_pty = POINTonE1_Serialize_BE(out, P);
        out[0] |= (sgn0_pty & 2) << 4;      /* pre-decorate */
        out[0] |= vec_is_zero(P->Z, sizeof(P->Z)) << 6;
    }
}

inline void blst_p1_mult(POINTonE1 *out, const POINTonE1 *a,
                                  const byte *scalar, size_t nbits)
{
    if (nbits < 176) {
        if (nbits)
            POINTonE1_mult_w4(out, a, scalar, nbits);
        else
            vec_zero(out, sizeof(*out));
    } else if (nbits <= 256) {
        union { vec256 l; pow256 s; } val;
        size_t i, j, top, mask = (size_t)0 - 1;

        /* this is not about constant-time-ness, but branch optimization */
        for (top = (nbits + 7)/8, i=0, j=0; i<sizeof(val.s);) {
            val.s[i++] = scalar[j] & mask;
            mask = 0 - ((i - top) >> (8*sizeof(top)-1));
            j += 1 & mask;
        }

        if (check_mod_v256(val.s, *get_BLST_r()))  /* z^4 is the formal limit */
            POINTonE1_mult_glv(out, a, val.s);
        else    /* should never be the case, added for formal completeness */
            POINTonE1_mult_w5(out, a, scalar, nbits);

        vec_zero(val.l, sizeof(val));
    } else {    /* should never be the case, added for formal completeness */
        POINTonE1_mult_w5(out, a, scalar, nbits);
    }
}

inline void blst_p1_unchecked_mult(POINTonE1 *out, const POINTonE1 *a,
                                            const byte *scalar, size_t nbits)
{
    if (nbits)
        POINTonE1_mult_w4(out, a, scalar, nbits);
    else
        vec_zero(out, sizeof(*out));
}

inline int blst_p1_affine_is_equal(const POINTonE1_affine *a,
                            const POINTonE1_affine *b)
{   return (int)vec_is_equal(a, b, sizeof(*a));   }

inline int blst_p1_is_inf(const POINTonE1 *p)
{   return (int)vec_is_zero(p->Z, sizeof(p->Z));   }

inline const POINTonE1 *blst_p1_generator(void)
{   return &BLS12_381_G1;   }

inline int blst_p1_affine_is_inf(const POINTonE1_affine *p)
{   return (int)vec_is_zero(p, sizeof(*p));   }

inline const POINTonE1_affine *blst_p1_affine_generator(void)
{   return (const POINTonE1_affine *)&BLS12_381_G1;   }

inline size_t blst_p1_sizeof(void)
{   return sizeof(POINTonE1);   }

inline size_t blst_p1_affine_sizeof(void)
{   return sizeof(POINTonE1_affine);   }
