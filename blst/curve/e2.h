/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include "blst/curve/point.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fr.h"
#include "blst/common/errors.h"
#include "blst/common/Logger.h"

#include "ec_mult.h"
#include "ec_ops.h"

#include <string>

namespace blst
{

    BLST_HD BLST_INLINE void mul_by_b_onE2(Fp2& out, const Fp2& in)
    {
        Fp::sub(out[0], in[0], in[1]);
        Fp::add(out[1], in[0], in[1]);
        Fp::lshift(out[0], out[0], 2);
        Fp::lshift(out[1], out[1], 2);
    }

    BLST_HD BLST_INLINE void mul_by_4b_onE2(Fp2& out, const Fp2& in)
    {
        Fp::sub(out[0], in[0], in[1]);
        Fp::add(out[1], in[0], in[1]);
        Fp::lshift(out[0], out[0], 4);
        Fp::lshift(out[1], out[1], 4);
    }

    template<>
    BLST_HD BLST_INLINE void POINTonE2::cneg(POINTonE2& p, bool_t cbit)
    {
        Fp2::cneg(p.Y, p.Y, cbit);
    }

    BLST_HD BLST_INLINE void blst_p2_cneg(POINTonE2& a, int cbit)
    {
        POINTonE2::cneg(a, is_zero(cbit) ^ 1);
    }

    BLST_HD BLST_INLINE void POINTonE2_from_Jacobian(POINTonE2& out, const POINTonE2& in)
    {
        Fp2 Z, ZZ;
        limb_t inf = vec_is_zero(in.Z.data, sizeof(in.Z));

        reciprocal_fp2(Z, in.Z);                           /* 1/Z */

        Fp2::sqr(ZZ, Z);
        Fp2::mul(out.X, in.X, ZZ);                         /* X = X/Z^2 */

        Fp2::mul(ZZ, ZZ, Z);
        Fp2::mul(out.Y, in.Y, ZZ);                         /* Y = Y/Z^3 */

        vec_select(&out.Z, &in.Z, &get_BLST_G2().Z,
            sizeof(get_BLST_G2().Z), inf);    /* Z = inf ? 0 : 1 */
    }

    template<>
    BLST_HD BLST_INLINE void POINTonE2::from_Jacobian(POINTonE2& out, const POINTonE2& in)
    {
		POINTonE2_from_Jacobian(out, in);
    }


    BLST_HD BLST_INLINE void blst_p2_from_jacobian(POINTonE2& out, const POINTonE2& a)
    {
        POINTonE2_from_Jacobian(out, a);
    }

    BLST_HD BLST_INLINE void POINTonE2_to_affine(POINTonE2_affine& out, const POINTonE2& in)
    {
        POINTonE2 p;
        auto ptr = &in;
        if (!vec_is_equal(ptr->Z.data, get_BLST_Rx().p2, sizeof(in.Z))) {
            POINTonE2_from_Jacobian(p, *ptr);
            ptr = &p;
        }
        vec_copy(&out, ptr, sizeof(out));
    }

    template<>
    BLST_HD BLST_INLINE void POINTonE2::to_affine(POINTonE2_affine& out, const POINTonE2& in)
    {
		POINTonE2_to_affine(out, in);
    }

    BLST_HD BLST_INLINE void blst_p2_to_affine(POINTonE2_affine& out, const POINTonE2& a)
    {
        POINTonE2_to_affine(out, a);
    }

    BLST_HD BLST_INLINE void blst_p2_from_affine(POINTonE2& out, const POINTonE2_affine& a)
    {
        vec_copy(&out, &a, sizeof(a));
        vec_select(&out.Z, &a.X, &get_BLST_Rx().p2, sizeof(out.Z),
            vec_is_zero(&a, sizeof(a)));
    }

    BLST_HD BLST_INLINE  bool_t POINTonE2_affine_on_curve(const POINTonE2_affine& p)
    {
        Fp2 XXX, YY;

        Fp2::sqr(XXX, p.X);
        Fp2::mul(XXX, XXX, p.X);                            /* X^3 */
        Fp2::add(XXX, XXX, get_BLST_B_E2());                            /* X^3 + B */

        Fp2::sqr(YY, p.Y);                                  /* Y^2 */

        return vec_is_equal(&XXX, &YY, sizeof(XXX));
    }

    BLST_HD BLST_INLINE int blst_p2_affine_on_curve(const POINTonE2_affine& p)
    {
        return (int)(POINTonE2_affine_on_curve(p) | vec_is_zero(&p, sizeof(p)));
    }

    BLST_HD BLST_INLINE bool_t POINTonE2_on_curve(const POINTonE2& p)
    {
        Fp2 XXX, YY, BZ6;
        limb_t inf = vec_is_zero(&p.Z, sizeof(p.Z));

        Fp2::sqr(BZ6, p.Z);
        Fp2::mul(BZ6, BZ6, p.Z);
        Fp2::sqr(XXX, BZ6);                                  /* Z^6 */
        mul_by_b_onE2(BZ6, XXX);                            /* B*Z^6 */

        Fp2::sqr(XXX, p.X);
        Fp2::mul(XXX, XXX, p.X);                            /* X^3 */
        Fp2::add(XXX, XXX, BZ6);                             /* X^3 + B*Z^6 */

        Fp2::sqr(YY, p.Y);                                  /* Y^2 */

        return vec_is_equal(&XXX, &YY, sizeof(XXX)) | inf;
    }

    BLST_HD BLST_INLINE int blst_p2_on_curve(const POINTonE2& p)
    {
        return (int)POINTonE2_on_curve(p);
    }

    BLST_HD BLST_INLINE limb_t POINTonE2_affine_Serialize_BE(unsigned char out[192],
        const POINTonE2_affine& in)
    {
        Fp2 temp;

        Fp::from_mont(temp[1], in.X[1]);
        be_bytes_from_limbs(out, temp[1], sizeof(temp[1]));
        Fp::from_mont(temp[0], in.X[0]);
        be_bytes_from_limbs(out + 48, temp[0], sizeof(temp[0]));

        Fp::from_mont(temp[1], in.Y[1]);
        be_bytes_from_limbs(out + 96, temp[1], sizeof(temp[1]));
        Fp::from_mont(temp[0], in.Y[0]);
        be_bytes_from_limbs(out + 144, temp[0], sizeof(temp[0]));

        return sgn0_pty_mod_v384x(temp, get_BLST_P());
    }

    BLST_HD BLST_INLINE void blst_p2_affine_serialize(unsigned char out[192],
        const POINTonE2_affine& in)
    {
        if (vec_is_zero(&in.X, 2 * sizeof(in.X))) {
            bytes_zero(out, 192);
            out[0] = 0x40;    /* infinity bit */
        }
        else {
            (void)POINTonE2_affine_Serialize_BE(out, in);
        }
    }

    BLST_HD BLST_INLINE  limb_t POINTonE2_Serialize_BE(unsigned char out[192],
        const POINTonE2& in)
    {
        POINTonE2 p;
        auto ptr = &in;

        if (!vec_is_equal(&ptr->Z, &get_BLST_Rx().p2, sizeof(in.Z))) {
            POINTonE2_from_Jacobian(p, *ptr);
            ptr = &p;
        }

        return POINTonE2_affine_Serialize_BE(out, (const POINTonE2_affine&)*ptr);
    }

    BLST_HD BLST_INLINE void POINTonE2_Serialize(unsigned char out[192], const POINTonE2& in)
    {
        if (vec_is_zero(&in.Z, sizeof(in.Z))) {
            bytes_zero(out, 192);
            out[0] = 0x40;    /* infinity bit */
        }
        else {
            (void)POINTonE2_Serialize_BE(out, in);
        }
    }

    BLST_HD BLST_INLINE void blst_p2_serialize(unsigned char out[192], const POINTonE2& in)
    {
        POINTonE2_Serialize(out, in);
    }

    BLST_HD BLST_INLINE limb_t POINTonE2_affine_Compress_BE(unsigned char out[96],
        const POINTonE2_affine& in)
    {
        Fp temp;

        Fp::from_mont(temp, in.X[1]);
        be_bytes_from_limbs(out, temp, sizeof(temp));
        Fp::from_mont(temp, in.X[0]);
        be_bytes_from_limbs(out + 48, temp, sizeof(temp));

        return sgn0_pty_mont_v384x(in.Y, get_BLST_P(), BLST_p0);
    }

    BLST_HD BLST_INLINE void blst_p2_affine_compress(unsigned char out[96], const POINTonE2_affine& in)
    {
        if (vec_is_zero(&in.X, 2 * sizeof(in.X))) {
            bytes_zero(out, 96);
            out[0] = 0xc0;    /* compressed and infinity bits */
        }
        else {
            limb_t sign = POINTonE2_affine_Compress_BE(out, in);
            out[0] |= (unsigned char)(0x80 | ((sign & 2) << 4));
        }
    }

    BLST_HD BLST_INLINE limb_t POINTonE2_Compress_BE(unsigned char out[96],
        const POINTonE2& in)
    {
        POINTonE2 p;
        auto ptr = &in;

        if (!vec_is_equal(&ptr->Z, &get_BLST_Rx().p, sizeof(in.Z))) {
            POINTonE2_from_Jacobian(p, *ptr);
            ptr = &p;
        }

        return POINTonE2_affine_Compress_BE(out, (const POINTonE2_affine&)*ptr);
    }

    BLST_HD BLST_INLINE void blst_p2_compress(unsigned char out[96], const POINTonE2& in)
    {
        if (vec_is_zero(&in.Z, sizeof(in.Z))) {
            bytes_zero(out, 96);
            out[0] = 0xc0;    /* compressed and infinity bits */
        }
        else {
            limb_t sign = POINTonE2_Compress_BE(out, in);
            out[0] |= (unsigned char)(0x80 | ((sign & 2) << 4));
        }
    }

    BLST_HD BLST_INLINE limb_t POINTonE2_Uncompress_BE(POINTonE2_affine& out,
        const unsigned char in[96])
    {
        POINTonE2_affine ret;
        Fp temp;

        limbs_from_be_bytes(ret.X[1], in, sizeof(ret.X[1]));
        limbs_from_be_bytes(ret.X[0], in + 48, sizeof(ret.X[0]));

        /* clear top 3 bits in case caller was conveying some information there */
        ret.X[1][sizeof(ret.X[1]) / sizeof(limb_t) - 1] &= ((limb_t)0 - 1) >> 3;
        Fp::add(temp, ret.X[1], Fp::zero());  /* less than modulus? */
        if (!vec_is_equal(temp, ret.X[1], sizeof(temp)))
            return (limb_t)0 - BLST_BAD_ENCODING;

        Fp::add(temp, ret.X[0], Fp::zero());  /* less than modulus? */
        if (!vec_is_equal(temp, ret.X[0], sizeof(temp)))
            return (limb_t)0 - BLST_BAD_ENCODING;

        Fp::mul(ret.X[0], ret.X[0], get_BLST_RR());
        Fp::mul(ret.X[1], ret.X[1], get_BLST_RR());

        Fp2::sqr(ret.Y, ret.X);
        Fp2::mul(ret.Y, ret.Y, ret.X);
        Fp2::add(ret.Y, ret.Y, get_BLST_B_E2());                        /* X^3 + B */
        if (!sqrt_fp2(ret.Y, ret.Y))
            return (limb_t)0 - BLST_POINT_NOT_ON_CURVE;

        vec_copy(&out, &ret, sizeof(ret));

        return sgn0_pty_mont_v384x(out.Y, get_BLST_P(), BLST_p0);
    }

    BLST_HD BLST_INLINE BLST_ERROR POINTonE2_Uncompress_Z(POINTonE2_affine& out,
        const unsigned char in[96])
    {
        unsigned char in0 = in[0];
        limb_t sgn0_pty;

        if ((in0 & 0x80) == 0)      /* compressed bit */
            return BLST_BAD_ENCODING;

        if (in0 & 0x40) {           /* infinity bit */
            if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in + 1, 95)) {
                vec_zero(&out, sizeof(out));
                return BLST_SUCCESS;
            }
            else {
                return BLST_BAD_ENCODING;
            }
        }

        sgn0_pty = POINTonE2_Uncompress_BE(out, in);

        if (sgn0_pty > 3)
            return (BLST_ERROR)(0 - sgn0_pty); /* POINT_NOT_ON_CURVE */

        sgn0_pty >>= 1; /* skip over parity bit */
        sgn0_pty ^= (in0 & 0x20) >> 5;
        Fp2::cneg(out.Y, out.Y, sgn0_pty);

        return BLST_SUCCESS;
    }

    BLST_HD BLST_INLINE BLST_ERROR blst_p2_uncompress(POINTonE2_affine& out, const unsigned char in[96])
    {
        return POINTonE2_Uncompress_Z(out, in);
    }

    BLST_HD BLST_INLINE BLST_ERROR POINTonE2_Deserialize_BE(POINTonE2_affine& out,
        const unsigned char in[192])
    {
        POINTonE2_affine ret;
        Fp temp;

        limbs_from_be_bytes(ret.X[1], in, sizeof(ret.X[1]));
        limbs_from_be_bytes(ret.X[0], in + 48, sizeof(ret.X[0]));
        limbs_from_be_bytes(ret.Y[1], in + 96, sizeof(ret.Y[1]));
        limbs_from_be_bytes(ret.Y[0], in + 144, sizeof(ret.Y[0]));

        /* clear top 3 bits in case caller was conveying some information there */
        ret.X[1][sizeof(ret.X[1]) / sizeof(limb_t) - 1] &= ((limb_t)0 - 1) >> 3;
        Fp::add(temp, ret.X[1], Fp::zero());  /* less than modulus? */
        if (!vec_is_equal(temp, ret.X[1], sizeof(temp)))
            return BLST_BAD_ENCODING;

        Fp::add(temp, ret.X[0], Fp::zero());  /* less than modulus? */
        if (!vec_is_equal(temp, ret.X[0], sizeof(temp)))
            return BLST_BAD_ENCODING;

        Fp::add(temp, ret.Y[1], Fp::zero());  /* less than modulus? */
        if (!vec_is_equal(temp, ret.Y[1], sizeof(temp)))
            return BLST_BAD_ENCODING;

        Fp::add(temp, ret.Y[0], Fp::zero());  /* less than modulus? */
        if (!vec_is_equal(temp, ret.Y[0], sizeof(temp)))
            return BLST_BAD_ENCODING;

        Fp::mul(ret.X[0], ret.X[0], get_BLST_RR());
        Fp::mul(ret.X[1], ret.X[1], get_BLST_RR());
        Fp::mul(ret.Y[0], ret.Y[0], get_BLST_RR());
        Fp::mul(ret.Y[1], ret.Y[1], get_BLST_RR());

        if (!POINTonE2_affine_on_curve(ret))
            return BLST_POINT_NOT_ON_CURVE;

        vec_copy(&out, &ret, sizeof(ret));

        return BLST_SUCCESS;
    }

    BLST_HD BLST_INLINE BLST_ERROR POINTonE2_Deserialize_Z(POINTonE2_affine& out,
        const unsigned char in[192])
    {
        unsigned char in0 = in[0];

        if ((in0 & 0xe0) == 0)
            return POINTonE2_Deserialize_BE(out, in);

        if (in0 & 0x80)             /* compressed bit */
            return POINTonE2_Uncompress_Z(out, in);

        if (in0 & 0x40) {           /* infinity bit */
            if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in + 1, 191)) {
                vec_zero(&out, sizeof(out));
                return BLST_SUCCESS;
            }
        }

        return BLST_BAD_ENCODING;
    }

    BLST_HD BLST_INLINE BLST_ERROR blst_p2_deserialize(POINTonE2_affine& out,
        const unsigned char in[192])
    {
        return POINTonE2_Deserialize_Z(out, in);
    }

    //POINT_DADD_IMPL(POINTonE2, 384x, fp2)
        //POINT_DADD_AFFINE_IMPL_A0(POINTonE2, 384x, fp2, get_BLST_Rx().p2)
        //POINT_ADD_IMPL(POINTonE2, 384x, fp2)
        //POINT_ADD_AFFINE_IMPL(POINTonE2, 384x, fp2, get_BLST_Rx().p2)
        //POINT_DOUBLE_IMPL_A0(POINTonE2, 384x, fp2)
        //POINT_IS_EQUAL_IMPL(POINTonE2, 384x, fp2)

        BLST_HD BLST_INLINE void blst_p2_add(POINTonE2& out, const POINTonE2& a, const POINTonE2& b)
    {
        POINTonE2::add(out, a, b);
    }

    BLST_HD BLST_INLINE void blst_p2_add_or_double(POINTonE2& out, const POINTonE2& a,
        const POINTonE2& b)
    {
        POINTonE2::dadd(out, a, b, NULL);
    }

    BLST_HD BLST_INLINE void blst_p2_add_affine(POINTonE2& out, const POINTonE2& a,
        const POINTonE2_affine& b)
    {
        POINTonE2::add_affine(out, a, b);
    }

    BLST_HD BLST_INLINE void blst_p2_add_or_double_affine(POINTonE2& out, const POINTonE2& a,
        const POINTonE2_affine& b)
    {
        POINTonE2::dadd_affine(out, a, b);
    }

    BLST_HD BLST_INLINE void blst_p2_double(POINTonE2& out, const POINTonE2& a)
    {
        POINTonE2::doubling(out, a);
    }

    BLST_HD BLST_INLINE int blst_p2_is_equal(const POINTonE2& a, const POINTonE2& b)
    {
        return (int)POINTonE2::is_equal(a, b);
    }

    //POINT_MULT_SCALAR_WX_IMPL(POINTonE2, 4)
    //    POINT_MULT_SCALAR_WX_IMPL(POINTonE2, 5)

#ifdef __BLST_PRIVATE_TESTMODE__
        POINT_AFFINE_MULT_SCALAR_IMPL(POINTonE2)

        DECLARE_PRIVATE_POINTXZ(POINTonE2, 384x)
    POINT_LADDER_PRE_IMPL(POINTonE2, 384x, fp2)
        POINT_LADDER_STEP_IMPL_A0(POINTonE2, 384x, fp2, onE2)
        POINT_LADDER_POST_IMPL_A0(POINTonE2, 384x, fp2, onE2)
        POINT_MULT_SCALAR_LADDER_IMPL(POINTonE2)
#endif

        BLST_HD BLST_INLINE void psi(POINTonE2& out, const POINTonE2& in)
    {
        const Fp2 frobenius_x = { /* 1/(1 + i)^((P-1)/3) */
          Fp{ 0 },
          Fp{ /* (0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4
                  897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad << 384) % P */
            TO_LIMB_T(0x890dc9e4867545c3), TO_LIMB_T(0x2af322533285a5d5),
            TO_LIMB_T(0x50880866309b7e2c), TO_LIMB_T(0xa20d1b8c7e881024),
            TO_LIMB_T(0x14e4f04fe2db9068), TO_LIMB_T(0x14e56d3f1564853a) }
        };
        const Fp2 frobenius_y = { /* 1/(1 + i)^((P-1)/2) */
          Fp{ /* (0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60
                  ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2 << 384) % P */
            TO_LIMB_T(0x3e2f585da55c9ad1), TO_LIMB_T(0x4294213d86c18183),
            TO_LIMB_T(0x382844c88b623732), TO_LIMB_T(0x92ad2afd19103e18),
            TO_LIMB_T(0x1d794e4fac7cf0b9), TO_LIMB_T(0x0bd592fc7d825ec8) },
          Fp{ /* (0x06af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e
                  77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09 << 384) % P */
            TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),
            TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),
            TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2) },
        };

        vec_copy(&out, &in, sizeof(out));
        Fp::cneg(out.X[1], out.X[1], 1);   Fp2::mul(out.X, out.X, frobenius_x);
        Fp::cneg(out.Y[1], out.Y[1], 1);   Fp2::mul(out.Y, out.Y, frobenius_y);
        Fp::cneg(out.Z[1], out.Z[1], 1);
    }

    BLST_HD BLST_INLINE u64 hash(const void* d, size_t n) {
		u64 s = 234523763451;
		const byte* p = (const byte*)d;
		for (size_t i = 0; i < n; i++) {
			s = (s << 5) + s + p[i];
			s *= 1234567891;
		}

        return s;
	}

    /* Galbraith-Lin-Scott, ~67% faster than POINTonE2_mul_w5 */
    BLST_HD BLST_INLINE void POINTonE2_mult_gls(POINTonE2& out, const POINTonE2& in,
        const pow256 SK)
    {
        Fr l;
        pow256 s; 

        /* break down SK to "digits" with |z| as radix [in constant time] */
        limbs_from_le_bytes(l, SK, 32);
        div_by_zz(l);
        div_by_z(l);
        div_by_z(l + NLIMBS(256) / 2);


        le_bytes_from_limbs(s, l, 32);

        {
            const byte* scalars[2] = { s, NULL };
            POINTonE2 table[4][1 << (5 - 1)];   /* 18KB */
            size_t i;

            //POINTonE2_precompute_w5(table[0], in);
			precompute_w<POINTonE2,5>(table[0], in);

            for (i = 0; i < 1 << (5 - 1); i++) {
                psi(table[1][i], table[0][i]);
                psi(table[2][i], table[1][i]);
                psi(table[3][i], table[2][i]);
                POINTonE2::cneg(table[1][i], 1); /* account for z being negative */
                POINTonE2::cneg(table[3][i], 1);
            }

            s_mult_w<POINTonE2, 5>(out, NULL, 4, scalars, 64, table);
        }

		// copy back the scalar to scrub the stack
        for(u64 i = 0; i < NLIMBS(256); i++) {
            l[i] = s[i];
		}

        vec_zero(l, sizeof(l));   /* scrub the copy of SK */
    }

    BLST_HD BLST_INLINE void POINTonE2_sign(POINTonE2& out, const POINTonE2& in, const pow256 SK)
    {
        Fp2 Z, ZZ;
        limb_t inf;

        POINTonE2_mult_gls(out, in, SK);

        /* convert to affine to remove possible bias in out.Z */
        inf = vec_is_zero(&out.Z, sizeof(out.Z));
#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
        flt_reciprocal_fp2(Z, out.Z);                      /* 1/Z   */
#else
        reciprocal_fp2(Z, out.Z);                          /* 1/Z   */
#endif

        Fp2::sqr(ZZ, Z);
        Fp2::mul(out.X, out.X, ZZ);                        /* X = X/Z^2 */

        Fp2::mul(ZZ, ZZ, Z);
        Fp2::mul(out.Y, out.Y, ZZ);                        /* Y = Y/Z^3 */

        vec_select(&out.Z, &out.Z, &get_BLST_G2().Z, sizeof(get_BLST_G2().Z),
            inf);                            /* Z = inf ? 0 : 1 */
    }

    BLST_HD BLST_INLINE void blst_sk_to_pk_in_g2(POINTonE2& out, const pow256 SK)
    {
        POINTonE2_sign(out, get_BLST_G2(), SK);
    }

    BLST_HD BLST_INLINE void blst_sign_pk_in_g1(POINTonE2& out, const POINTonE2& msg, const pow256 SK)
    {
        POINTonE2_sign(out, msg, SK);
    }

    BLST_HD BLST_INLINE void blst_sk_to_pk2_in_g2(unsigned char out[192], POINTonE2_affine* PK,
        const pow256 SK)
    {
        POINTonE2 P;

        POINTonE2_sign(P, get_BLST_G2(), SK);
        if (PK != NULL)
            vec_copy(PK,& P, sizeof(*PK));
        if (out != NULL) {
            limb_t sgn0_pty = POINTonE2_Serialize_BE(out, P);
            out[0] |= (sgn0_pty & 2) << 4;      /* pre-decorate */
            out[0] |= vec_is_zero(&P.Z, sizeof(P.Z)) << 6;
        }
    }

    BLST_HD BLST_INLINE void blst_sign_pk2_in_g1(unsigned char out[192], POINTonE2_affine* sig,
        const POINTonE2& hash, const pow256 SK)
    {
        POINTonE2 P;

        POINTonE2_sign(P, hash, SK);
        if (sig != NULL)
            vec_copy(sig, &P, sizeof(sig));
        if (out != NULL) {
            limb_t sgn0_pty = POINTonE2_Serialize_BE(out, P);
            out[0] |= (sgn0_pty & 2) << 4;      /* pre-decorate */
            out[0] |= vec_is_zero(&P.Z, sizeof(P.Z)) << 6;
        }
    }

    BLST_HD BLST_INLINE void blst_p2_mult(POINTonE2& out, const POINTonE2& a,
        const byte* scalar, size_t nbits)
    {
        if (nbits < 144) {
            if (nbits)
                mult_w<POINTonE2,4>(out, a, scalar, nbits);
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

            if (Fr::check_mod(val.l))  /* z^4 is the formal limit */
                POINTonE2_mult_gls(out, a, val.s);
            else    /* should never be the case, added for formal completeness */
                mult_w<POINTonE2,5>(out, a, scalar, nbits);

            vec_zero(val.l, sizeof(val));
        }
        else {    /* should never be the case, added for formal completeness */
            mult_w<POINTonE2,5>(out, a, scalar, nbits);
        }
    }

    BLST_HD BLST_INLINE void blst_p2_unchecked_mult(POINTonE2& out, const POINTonE2& a,
        const byte* scalar, size_t nbits)
    {
        if (nbits)
            mult_w< POINTonE2,4>(out, a, scalar, nbits);
        else
            vec_zero(&out, sizeof(out));
    }

    BLST_HD BLST_INLINE int blst_p2_affine_is_equal(const POINTonE2_affine& a,
        const POINTonE2_affine& b)
    {
        return (int)vec_is_equal(&a, &b, sizeof(a));
    }

    BLST_HD BLST_INLINE int blst_p2_is_inf(const POINTonE2& p)
    {
        return (int)vec_is_zero(&p.Z, sizeof(p.Z));
    }

    BLST_HD BLST_INLINE const POINTonE2& blst_p2_generator(void)
    {
        return get_BLST_G2();
    }

    BLST_HD BLST_INLINE int blst_p2_affine_is_inf(const POINTonE2_affine& p)
    {
        return (int)vec_is_zero(&p, sizeof(p));
    }

    BLST_HD BLST_INLINE const POINTonE2_affine& blst_p2_affine_generator(void)
    {
        return (const POINTonE2_affine&)get_BLST_G2();
    }

    BLST_HD BLST_INLINE size_t blst_p2_sizeof(void)
    {
        return sizeof(POINTonE2);
    }

    BLST_HD BLST_INLINE size_t blst_p2_affine_sizeof(void)
    {
        return sizeof(POINTonE2_affine);
    }

}
//BLST_EXTERN_C_END