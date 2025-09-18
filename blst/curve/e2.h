/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/curve/point.h"
#include "blst/fields/Fp2.h"
#include "blst/common/errors.h"
//BLST_EXTERN_C_BEGIN


HOST_DEVICE_FUNC BLST_INLINE void mul_by_b_onE2(vec384x out, const vec384x in)
{
    sub_fp(out[0], in[0], in[1]);
    add_fp(out[1], in[0], in[1]);
    lshift_fp(out[0], out[0], 2);
    lshift_fp(out[1], out[1], 2);
}

HOST_DEVICE_FUNC BLST_INLINE void mul_by_4b_onE2(vec384x out, const vec384x in)
{
    sub_fp(out[0], in[0], in[1]);
    add_fp(out[1], in[0], in[1]);
    lshift_fp(out[0], out[0], 4);
    lshift_fp(out[1], out[1], 4);
}

HOST_DEVICE_FUNC BLST_INLINE void POINTonE2_cneg(POINTonE2 *p, bool_t cbit)
{   cneg_fp2(p->Y, p->Y, cbit);   }

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_cneg(POINTonE2 *a, int cbit)
{   POINTonE2_cneg(a, is_zero(cbit) ^ 1);   }

 HOST_DEVICE_FUNC BLST_INLINE void POINTonE2_from_Jacobian(POINTonE2 *out, const POINTonE2 *in)
{
    vec384x Z, ZZ;
    limb_t inf = vec_is_zero(in->Z, sizeof(in->Z));

    reciprocal_fp2(Z, in->Z);                           /* 1/Z */

    sqr_fp2(ZZ, Z);
    mul_fp2(out->X, in->X, ZZ);                         /* X = X/Z^2 */

    mul_fp2(ZZ, ZZ, Z);
    mul_fp2(out->Y, in->Y, ZZ);                         /* Y = Y/Z^3 */

    vec_select(out->Z, in->Z, get_BLST_G2()->Z,
                       sizeof(get_BLST_G2()->Z), inf);    /* Z = inf ? 0 : 1 */
}

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_from_jacobian(POINTonE2 *out, const POINTonE2 *a)
{   POINTonE2_from_Jacobian(out, a);   }

HOST_DEVICE_FUNC BLST_INLINE void POINTonE2_to_affine(POINTonE2_affine *out, const POINTonE2 *in)
{
    POINTonE2 p;

    if (!vec_is_equal(in->Z, get_BLST_Rx()->p2, sizeof(in->Z))) {
        POINTonE2_from_Jacobian(&p, in);
        in = &p;
    }
    vec_copy(out, in, sizeof(*out));
}

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_to_affine(POINTonE2_affine *out, const POINTonE2 *a)
{   POINTonE2_to_affine(out, a);   }

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_from_affine(POINTonE2 *out, const POINTonE2_affine *a)
{
    vec_copy(out, a, sizeof(*a));
    vec_select(out->Z, a->X, get_BLST_Rx()->p2, sizeof(out->Z),
                       vec_is_zero(a, sizeof(*a)));
}

HOST_DEVICE_FUNC BLST_INLINE  bool_t POINTonE2_affine_on_curve(const POINTonE2_affine *p)
{
    vec384x XXX, YY;

    sqr_fp2(XXX, p->X);
    mul_fp2(XXX, XXX, p->X);                            /* X^3 */
    add_fp2(XXX, XXX, *get_BLST_B_E2());                            /* X^3 + B */

    sqr_fp2(YY, p->Y);                                  /* Y^2 */

    return vec_is_equal(XXX, YY, sizeof(XXX));
}

HOST_DEVICE_FUNC BLST_INLINE int blst_p2_affine_on_curve(const POINTonE2_affine *p)
{   return (int)(POINTonE2_affine_on_curve(p) | vec_is_zero(p, sizeof(*p)));   }

HOST_DEVICE_FUNC BLST_INLINE bool_t POINTonE2_on_curve(const POINTonE2 *p)
{
    vec384x XXX, YY, BZ6;
    limb_t inf = vec_is_zero(p->Z, sizeof(p->Z));

    sqr_fp2(BZ6, p->Z);
    mul_fp2(BZ6, BZ6, p->Z);
    sqr_fp2(XXX, BZ6);                                  /* Z^6 */
    mul_by_b_onE2(BZ6, XXX);                            /* B*Z^6 */

    sqr_fp2(XXX, p->X);
    mul_fp2(XXX, XXX, p->X);                            /* X^3 */
    add_fp2(XXX, XXX, BZ6);                             /* X^3 + B*Z^6 */

    sqr_fp2(YY, p->Y);                                  /* Y^2 */

    return vec_is_equal(XXX, YY, sizeof(XXX)) | inf;
}

HOST_DEVICE_FUNC BLST_INLINE int blst_p2_on_curve(const POINTonE2 *p)
{   return (int)POINTonE2_on_curve(p);   }

HOST_DEVICE_FUNC BLST_INLINE limb_t POINTonE2_affine_Serialize_BE(unsigned char out[192],
                                            const POINTonE2_affine *in)
{
    vec384x temp;

    from_fp(temp[1], in->X[1]);
    be_bytes_from_limbs(out, temp[1], sizeof(temp[1]));
    from_fp(temp[0], in->X[0]);
    be_bytes_from_limbs(out + 48, temp[0], sizeof(temp[0]));

    from_fp(temp[1], in->Y[1]);
    be_bytes_from_limbs(out + 96, temp[1], sizeof(temp[1]));
    from_fp(temp[0], in->Y[0]);
    be_bytes_from_limbs(out + 144, temp[0], sizeof(temp[0]));

    return sgn0_pty_mod_v384x(temp, *get_BLST_P());
}

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_affine_serialize(unsigned char out[192],
                              const POINTonE2_affine *in)
{
    if (vec_is_zero(in->X, 2*sizeof(in->X))) {
        bytes_zero(out, 192);
        out[0] = 0x40;    /* infinity bit */
    } else {
        (void)POINTonE2_affine_Serialize_BE(out, in);
    }
}

HOST_DEVICE_FUNC BLST_INLINE  limb_t POINTonE2_Serialize_BE(unsigned char out[192],
                                     const POINTonE2 *in)
{
    POINTonE2 p;

    if (!vec_is_equal(in->Z, get_BLST_Rx()->p2, sizeof(in->Z))) {
        POINTonE2_from_Jacobian(&p, in);
        in = &p;
    }

    return POINTonE2_affine_Serialize_BE(out, (const POINTonE2_affine *)in);
}

HOST_DEVICE_FUNC BLST_INLINE void POINTonE2_Serialize(unsigned char out[192], const POINTonE2 *in)
{
    if (vec_is_zero(in->Z, sizeof(in->Z))) {
        bytes_zero(out, 192);
        out[0] = 0x40;    /* infinity bit */
    } else {
        (void)POINTonE2_Serialize_BE(out, in);
    }
}

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_serialize(unsigned char out[192], const POINTonE2 *in)
{   POINTonE2_Serialize(out, in);   }

HOST_DEVICE_FUNC BLST_INLINE limb_t POINTonE2_affine_Compress_BE(unsigned char out[96],
                                           const POINTonE2_affine *in)
{
    vec384 temp;

    from_fp(temp, in->X[1]);
    be_bytes_from_limbs(out, temp, sizeof(temp));
    from_fp(temp, in->X[0]);
    be_bytes_from_limbs(out + 48, temp, sizeof(temp));

    return sgn0_pty_mont_v384x(in->Y, *get_BLST_P(), BLST_p0);
}

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_affine_compress(unsigned char out[96], const POINTonE2_affine *in)
{
    if (vec_is_zero(in->X, 2*sizeof(in->X))) {
        bytes_zero(out, 96);
        out[0] = 0xc0;    /* compressed and infinity bits */
    } else {
        limb_t sign = POINTonE2_affine_Compress_BE(out, in);
        out[0] |= (unsigned char)(0x80 | ((sign & 2) << 4));
    }
}

HOST_DEVICE_FUNC BLST_INLINE limb_t POINTonE2_Compress_BE(unsigned char out[96],
                                    const POINTonE2 *in)
{
    POINTonE2 p;

    if (!vec_is_equal(in->Z, get_BLST_Rx()->p, sizeof(in->Z))) {
        POINTonE2_from_Jacobian(&p, in);
        in = &p;
    }

    return POINTonE2_affine_Compress_BE(out, (const POINTonE2_affine *)in);
}

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_compress(unsigned char out[96], const POINTonE2 *in)
{
    if (vec_is_zero(in->Z, sizeof(in->Z))) {
        bytes_zero(out, 96);
        out[0] = 0xc0;    /* compressed and infinity bits */
    } else {
        limb_t sign = POINTonE2_Compress_BE(out, in);
        out[0] |= (unsigned char)(0x80 | ((sign & 2) << 4));
    }
}

HOST_DEVICE_FUNC BLST_INLINE limb_t POINTonE2_Uncompress_BE(POINTonE2_affine *out,
                                      const unsigned char in[96])
{
    POINTonE2_affine ret;
    vec384 temp;

    limbs_from_be_bytes(ret.X[1], in, sizeof(ret.X[1]));
    limbs_from_be_bytes(ret.X[0], in + 48, sizeof(ret.X[0]));

    /* clear top 3 bits in case caller was conveying some information there */
    ret.X[1][sizeof(ret.X[1])/sizeof(limb_t)-1] &= ((limb_t)0-1) >> 3;
    add_fp(temp, ret.X[1], ZERO_384);  /* less than modulus? */
    if (!vec_is_equal(temp, ret.X[1], sizeof(temp)))
        return (limb_t)0 - BLST_BAD_ENCODING;

    add_fp(temp, ret.X[0], ZERO_384);  /* less than modulus? */
    if (!vec_is_equal(temp, ret.X[0], sizeof(temp)))
        return (limb_t)0 - BLST_BAD_ENCODING;

    mul_fp(ret.X[0], ret.X[0], *get_BLST_RR());
    mul_fp(ret.X[1], ret.X[1], *get_BLST_RR());

    sqr_fp2(ret.Y, ret.X);
    mul_fp2(ret.Y, ret.Y, ret.X);
    add_fp2(ret.Y, ret.Y, *get_BLST_B_E2());                        /* X^3 + B */
    if (!sqrt_fp2(ret.Y, ret.Y))
        return (limb_t)0 - BLST_POINT_NOT_ON_CURVE;

    vec_copy(out, &ret, sizeof(ret));

    return sgn0_pty_mont_v384x(out->Y, *get_BLST_P(), BLST_p0);
}

HOST_DEVICE_FUNC BLST_INLINE BLST_ERROR POINTonE2_Uncompress_Z(POINTonE2_affine *out,
                                         const unsigned char in[96])
{
    unsigned char in0 = in[0];
    limb_t sgn0_pty;

    if ((in0 & 0x80) == 0)      /* compressed bit */
        return BLST_BAD_ENCODING;

    if (in0 & 0x40) {           /* infinity bit */
        if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in+1, 95)) {
            vec_zero(out, sizeof(*out));
            return BLST_SUCCESS;
        } else {
            return BLST_BAD_ENCODING;
        }
    }

    sgn0_pty = POINTonE2_Uncompress_BE(out, in);

    if (sgn0_pty > 3)
        return (BLST_ERROR)(0 - sgn0_pty); /* POINT_NOT_ON_CURVE */

    sgn0_pty >>= 1; /* skip over parity bit */
    sgn0_pty ^= (in0 & 0x20) >> 5;
    cneg_fp2(out->Y, out->Y, sgn0_pty);

    return BLST_SUCCESS;
}

HOST_DEVICE_FUNC BLST_INLINE BLST_ERROR blst_p2_uncompress(POINTonE2_affine *out, const unsigned char in[96])
{   return POINTonE2_Uncompress_Z(out, in);   }

HOST_DEVICE_FUNC BLST_INLINE BLST_ERROR POINTonE2_Deserialize_BE(POINTonE2_affine *out,
                                           const unsigned char in[192])
{
    POINTonE2_affine ret;
    vec384 temp;

    limbs_from_be_bytes(ret.X[1], in, sizeof(ret.X[1]));
    limbs_from_be_bytes(ret.X[0], in + 48, sizeof(ret.X[0]));
    limbs_from_be_bytes(ret.Y[1], in + 96, sizeof(ret.Y[1]));
    limbs_from_be_bytes(ret.Y[0], in + 144, sizeof(ret.Y[0]));

    /* clear top 3 bits in case caller was conveying some information there */
    ret.X[1][sizeof(ret.X[1])/sizeof(limb_t)-1] &= ((limb_t)0-1) >> 3;
    add_fp(temp, ret.X[1], ZERO_384);  /* less than modulus? */
    if (!vec_is_equal(temp, ret.X[1], sizeof(temp)))
        return BLST_BAD_ENCODING;

    add_fp(temp, ret.X[0], ZERO_384);  /* less than modulus? */
    if (!vec_is_equal(temp, ret.X[0], sizeof(temp)))
        return BLST_BAD_ENCODING;

    add_fp(temp, ret.Y[1], ZERO_384);  /* less than modulus? */
    if (!vec_is_equal(temp, ret.Y[1], sizeof(temp)))
        return BLST_BAD_ENCODING;

    add_fp(temp, ret.Y[0], ZERO_384);  /* less than modulus? */
    if (!vec_is_equal(temp, ret.Y[0], sizeof(temp)))
        return BLST_BAD_ENCODING;

    mul_fp(ret.X[0], ret.X[0], *get_BLST_RR());
    mul_fp(ret.X[1], ret.X[1], *get_BLST_RR());
    mul_fp(ret.Y[0], ret.Y[0], *get_BLST_RR());
    mul_fp(ret.Y[1], ret.Y[1], *get_BLST_RR());

    if (!POINTonE2_affine_on_curve(&ret))
        return BLST_POINT_NOT_ON_CURVE;

    vec_copy(out, &ret, sizeof(ret));

    return BLST_SUCCESS;
}

HOST_DEVICE_FUNC BLST_INLINE BLST_ERROR POINTonE2_Deserialize_Z(POINTonE2_affine *out,
                                          const unsigned char in[192])
{
    unsigned char in0 = in[0];

    if ((in0 & 0xe0) == 0)
        return POINTonE2_Deserialize_BE(out, in);

    if (in0 & 0x80)             /* compressed bit */
        return POINTonE2_Uncompress_Z(out, in);

    if (in0 & 0x40) {           /* infinity bit */
        if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in+1, 191)) {
            vec_zero(out, sizeof(*out));
            return BLST_SUCCESS;
        }
    }

    return BLST_BAD_ENCODING;
}

HOST_DEVICE_FUNC BLST_INLINE BLST_ERROR blst_p2_deserialize(POINTonE2_affine *out,
                               const unsigned char in[192])
{   return POINTonE2_Deserialize_Z(out, in);   }

#include "ec_ops.h"
POINT_DADD_IMPL(POINTonE2, 384x, fp2)
POINT_DADD_AFFINE_IMPL_A0(POINTonE2, 384x, fp2, get_BLST_Rx()->p2)
POINT_ADD_IMPL(POINTonE2, 384x, fp2)
POINT_ADD_AFFINE_IMPL(POINTonE2, 384x, fp2, get_BLST_Rx()->p2)
POINT_DOUBLE_IMPL_A0(POINTonE2, 384x, fp2)
POINT_IS_EQUAL_IMPL(POINTonE2, 384x, fp2)

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_add(POINTonE2 *out, const POINTonE2 *a, const POINTonE2 *b)
{   POINTonE2_add(out, a, b);   }

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_add_or_double(POINTonE2 *out, const POINTonE2 *a,
                                           const POINTonE2 *b)
{   POINTonE2_dadd(out, a, b, NULL);   }

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_add_affine(POINTonE2 *out, const POINTonE2 *a,
                                        const POINTonE2_affine *b)
{   POINTonE2_add_affine(out, a, b);   }

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_add_or_double_affine(POINTonE2 *out, const POINTonE2 *a,
                                                  const POINTonE2_affine *b)
{   POINTonE2_dadd_affine(out, a, b);   }

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_double(POINTonE2 *out, const POINTonE2 *a)
{   POINTonE2_double(out, a);   }

HOST_DEVICE_FUNC BLST_INLINE int blst_p2_is_equal(const POINTonE2 *a, const POINTonE2 *b)
{   return (int)POINTonE2_is_equal(a, b);   }

#include "ec_mult.h"
POINT_MULT_SCALAR_WX_IMPL(POINTonE2, 4)
POINT_MULT_SCALAR_WX_IMPL(POINTonE2, 5)

#ifdef __BLST_PRIVATE_TESTMODE__
POINT_AFFINE_MULT_SCALAR_IMPL(POINTonE2)

DECLARE_PRIVATE_POINTXZ(POINTonE2, 384x)
POINT_LADDER_PRE_IMPL(POINTonE2, 384x, fp2)
POINT_LADDER_STEP_IMPL_A0(POINTonE2, 384x, fp2, onE2)
POINT_LADDER_POST_IMPL_A0(POINTonE2, 384x, fp2, onE2)
POINT_MULT_SCALAR_LADDER_IMPL(POINTonE2)
#endif

HOST_DEVICE_FUNC BLST_INLINE void psi(POINTonE2 *out, const POINTonE2 *in)
{
    const vec384x frobenius_x = { /* 1/(1 + i)^((P-1)/3) */
      { 0 },
      { /* (0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4
              897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad << 384) % P */
        TO_LIMB_T(0x890dc9e4867545c3), TO_LIMB_T(0x2af322533285a5d5),
        TO_LIMB_T(0x50880866309b7e2c), TO_LIMB_T(0xa20d1b8c7e881024),
        TO_LIMB_T(0x14e4f04fe2db9068), TO_LIMB_T(0x14e56d3f1564853a) }
    };
    const vec384x frobenius_y = { /* 1/(1 + i)^((P-1)/2) */
      { /* (0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60
              ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2 << 384) % P */
        TO_LIMB_T(0x3e2f585da55c9ad1), TO_LIMB_T(0x4294213d86c18183),
        TO_LIMB_T(0x382844c88b623732), TO_LIMB_T(0x92ad2afd19103e18),
        TO_LIMB_T(0x1d794e4fac7cf0b9), TO_LIMB_T(0x0bd592fc7d825ec8) },
      { /* (0x06af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e
              77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09 << 384) % P */
        TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),
        TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),
        TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2) },
    };

    vec_copy(out, in, sizeof(*out));
    cneg_fp(out->X[1], out->X[1], 1);   mul_fp2(out->X, out->X, frobenius_x);
    cneg_fp(out->Y[1], out->Y[1], 1);   mul_fp2(out->Y, out->Y, frobenius_y);
    cneg_fp(out->Z[1], out->Z[1], 1);
}

/* Galbraith-Lin-Scott, ~67% faster than POINTonE2_mul_w5 */
HOST_DEVICE_FUNC BLST_INLINE void POINTonE2_mult_gls(POINTonE2 *out, const POINTonE2 *in,
                               const pow256 SK)
{
    union { vec256 l; pow256 s; } val;

    /* break down SK to "digits" with |z| as radix [in constant time] */

    limbs_from_le_bytes(val.l, SK, 32);
    div_by_zz(val.l);
    div_by_z(val.l);
    div_by_z(val.l + NLIMBS(256)/2);
    le_bytes_from_limbs(val.s, val.l, 32);

    {
        const byte *scalars[2] = { val.s, NULL };
        POINTonE2 table[4][1<<(5-1)];   /* 18KB */
        size_t i;

        POINTonE2_precompute_w5(table[0], in);
        for (i = 0; i < 1<<(5-1); i++) {
            psi(&table[1][i], &table[0][i]);
            psi(&table[2][i], &table[1][i]);
            psi(&table[3][i], &table[2][i]);
            POINTonE2_cneg(&table[1][i], 1); /* account for z being negative */
            POINTonE2_cneg(&table[3][i], 1);
        }

        POINTonE2s_mult_w5(out, NULL, 4, scalars, 64, table);
    }

    vec_zero(val.l, sizeof(val));   /* scrub the copy of SK */
}

HOST_DEVICE_FUNC BLST_INLINE void POINTonE2_sign(POINTonE2 *out, const POINTonE2 *in, const pow256 SK)
{
    vec384x Z, ZZ;
    limb_t inf;

    POINTonE2_mult_gls(out, in, SK);

    /* convert to affine to remove possible bias in out->Z */
    inf = vec_is_zero(out->Z, sizeof(out->Z));
#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    flt_reciprocal_fp2(Z, out->Z);                      /* 1/Z   */
#else
    reciprocal_fp2(Z, out->Z);                          /* 1/Z   */
#endif

    sqr_fp2(ZZ, Z);
    mul_fp2(out->X, out->X, ZZ);                        /* X = X/Z^2 */

    mul_fp2(ZZ, ZZ, Z);
    mul_fp2(out->Y, out->Y, ZZ);                        /* Y = Y/Z^3 */

    vec_select(out->Z, out->Z, get_BLST_G2()->Z, sizeof(get_BLST_G2()->Z),
                       inf);                            /* Z = inf ? 0 : 1 */
}

HOST_DEVICE_FUNC BLST_INLINE void blst_sk_to_pk_in_g2(POINTonE2 *out, const pow256 SK)
{   POINTonE2_sign(out, get_BLST_G2(), SK);   }

HOST_DEVICE_FUNC BLST_INLINE void blst_sign_pk_in_g1(POINTonE2 *out, const POINTonE2 *msg, const pow256 SK)
{   POINTonE2_sign(out, msg, SK);   }

HOST_DEVICE_FUNC BLST_INLINE void blst_sk_to_pk2_in_g2(unsigned char out[192], POINTonE2_affine *PK,
                          const pow256 SK)
{
    POINTonE2 P[1];

    POINTonE2_sign(P, get_BLST_G2(), SK);
    if (PK != NULL)
        vec_copy(PK, P, sizeof(*PK));
    if (out != NULL) {
        limb_t sgn0_pty = POINTonE2_Serialize_BE(out, P);
        out[0] |= (sgn0_pty & 2) << 4;      /* pre-decorate */
        out[0] |= vec_is_zero(P->Z, sizeof(P->Z)) << 6;
    }
}

HOST_DEVICE_FUNC BLST_INLINE void blst_sign_pk2_in_g1(unsigned char out[192], POINTonE2_affine *sig,
                         const POINTonE2 *hash, const pow256 SK)
{
    POINTonE2 P[1];

    POINTonE2_sign(P, hash, SK);
    if (sig != NULL)
        vec_copy(sig, P, sizeof(*sig));
    if (out != NULL) {
        limb_t sgn0_pty = POINTonE2_Serialize_BE(out, P);
        out[0] |= (sgn0_pty & 2) << 4;      /* pre-decorate */
        out[0] |= vec_is_zero(P->Z, sizeof(P->Z)) << 6;
    }
}

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_mult(POINTonE2 *out, const POINTonE2 *a,
                                  const byte *scalar, size_t nbits)
{
    if (nbits < 144) {
        if (nbits)
            POINTonE2_mult_w4(out, a, scalar, nbits);
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
            POINTonE2_mult_gls(out, a, val.s);
        else    /* should never be the case, added for formal completeness */
            POINTonE2_mult_w5(out, a, scalar, nbits);

        vec_zero(val.l, sizeof(val));
    } else {    /* should never be the case, added for formal completeness */
        POINTonE2_mult_w5(out, a, scalar, nbits);
    }
}

HOST_DEVICE_FUNC BLST_INLINE void blst_p2_unchecked_mult(POINTonE2 *out, const POINTonE2 *a,
                                            const byte *scalar, size_t nbits)
{
    if (nbits)
        POINTonE2_mult_w4(out, a, scalar, nbits);
    else
        vec_zero(out, sizeof(*out));
}

HOST_DEVICE_FUNC BLST_INLINE int blst_p2_affine_is_equal(const POINTonE2_affine *a,
                            const POINTonE2_affine *b)
{   return (int)vec_is_equal(a, b, sizeof(*a));   }

HOST_DEVICE_FUNC BLST_INLINE int blst_p2_is_inf(const POINTonE2 *p)
{   return (int)vec_is_zero(p->Z, sizeof(p->Z));   }

HOST_DEVICE_FUNC BLST_INLINE const POINTonE2 *blst_p2_generator(void)
{   return get_BLST_G2();   }

HOST_DEVICE_FUNC BLST_INLINE int blst_p2_affine_is_inf(const POINTonE2_affine *p)
{   return (int)vec_is_zero(p, sizeof(*p));   }

HOST_DEVICE_FUNC BLST_INLINE const POINTonE2_affine *blst_p2_affine_generator(void)
{   return (const POINTonE2_affine *)get_BLST_G2();   }

HOST_DEVICE_FUNC BLST_INLINE size_t blst_p2_sizeof(void)
{   return sizeof(POINTonE2);   }

HOST_DEVICE_FUNC BLST_INLINE size_t blst_p2_affine_sizeof(void)
{   return sizeof(POINTonE2_affine);   }

//BLST_EXTERN_C_END