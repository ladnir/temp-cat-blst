/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_FIELDS_H__
#define __BLS12_381_ASM_FIELDS_H__

#include "blst/common/vect.h"
#include "blst/config.h"
#include "blst/common/defines.h" 
#include "blst/fields/recip-addchain.h"
#include "blst/fields/exp.h"
#include "blst/fields/Fp.h"


namespace blst
{


    struct Fp2
    {
        Fp data[2];

        BLST_HD Fp& operator[](size_t idx) { return data[idx]; }
        BLST_HD const Fp& operator[](size_t idx) const { return data[idx]; }


        BLST_HD static void add(Fp2& ret, const Fp2& a, const Fp2& b);
        BLST_HD static void sub(Fp2& ret, const Fp2& a, const Fp2& b);
        BLST_HD static void mul(Fp2& ret, const Fp2& a, const Fp2& b);
        BLST_HD static void sqr(Fp2& ret, const Fp2& a);
        BLST_HD static void exp(Fp2& out, const Fp2& inp, const byte* pow, size_t pow_bits);


        BLST_HD static void cneg(Fp2& ret, const Fp2& a, bool_t flag);
        BLST_HD static void neg(Fp2& ret, const Fp2& a);
        BLST_HD static void from_mont(Fp2& ret, const Fp2& a);
        BLST_HD static void to_mont(Fp2& ret, const Fp2& a);
        
        BLST_HD static void reciprocal(Fp2& out, const Fp2& inp);
        BLST_HD static void flt_reciprocal(Fp2& out, const Fp2& inp);

        BLST_HD static bool_t sqrt(Fp2& out, const Fp2& inp);
        BLST_HD static void mul_by_3(Fp2& ret, const Fp2& a);
        BLST_HD static void mul_by_8(Fp2& ret, const Fp2& a);


        BLST_HD static void lshift(Fp2& ret, const Fp2& a, size_t count);


        BLST_HD static void conjugate(Fp2& ret, const Fp2& a);


        BLST_HD static void mul_by_u_plus_1(Fp2& ret, const Fp2& a);

        //BLST_HD static int is_square(const Fp2& inp);
        BLST_HD static int is_square(const Fp2& inp);
        

        BLST_HD static bool_t recip_sqrt(Fp2& out, const Fp2& inp);


        BLST_HD static bool_t sqrt_align(Fp2& out, const Fp2& ret,
            const Fp2& sqrt, const Fp2& inp);


        BLST_HD static limb_t sgn0_pty_mont(const Fp2& a);

        BLST_HD static limb_t sgn0_pty(const Fp2& a);


        BLST_HD static const Fp2& zero();
        BLST_HD static const Fp2& one();
    };



#define BLST_SQRT_SQRT_MINUS_1_LIMBS                                \
Fp{ TO_LIMB_T(0x3e2f585da55c9ad1), TO_LIMB_T(0x4294213d86c18183),     \
  TO_LIMB_T(0x382844c88b623732), TO_LIMB_T(0x92ad2afd19103e18),     \
  TO_LIMB_T(0x1d794e4fac7cf0b9), TO_LIMB_T(0x0bd592fc7d825ec8)},    \
Fp{ TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),     \
  TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),     \
  TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2) }

#define BLST_SQRT_MINUS_SQRT_MINUS_1_LIMBS                          \
Fp{ TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),     \
  TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),     \
  TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2) },   \
Fp{ TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),     \
  TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),     \
  TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2) }    



#if defined(BLST_CUDA) && defined(__CUDACC__)

    static __constant__ Fp2 BLST_sqrt_minus_1_DEV = Fp2({ Fp{ { 0 } }, Fp{ ONE_MONT_P_LIMBS } });
    /*
        * "magic" number is ±2^((p-3)/4)%p, which is "1/sqrt(2)",
        * in quotes because 2*"1/sqrt(2)"^2 == -1 mod p, not 1,
        * but it pivots into "complex" plane nevertheless...
        */
    static __constant__ Fp2 BLST_sqrt_sqrt_minus_1_DEV = { BLST_SQRT_SQRT_MINUS_1_LIMBS };
    static __constant__ Fp2 BLST_sqrt_minus_sqrt_minus_1_DEV = { BLST_SQRT_MINUS_SQRT_MINUS_1_LIMBS };

#endif // __CUDA_ARCH__





    extern const Fp2 BLST_sqrt_minus_1_HOST;// = { { 0 }, { ONE_MONT_P_LIMBS } };
    /*
        * "magic" number is ±2^((p-3)/4)%p, which is "1/sqrt(2)",
        * in quotes because 2*"1/sqrt(2)"^2 == -1 mod p, not 1,
        * but it pivots into "complex" plane nevertheless...
        */
    extern const Fp2 BLST_sqrt_sqrt_minus_1_HOST; /*= { BLST_SQRT_SQRT_MINUS_1_LIMBS };*/
    extern const Fp2 BLST_sqrt_minus_sqrt_minus_1_HOST; /*= { BLST_SQRT_MINUS_SQRT_MINUS_1_LIMBS };*/


    BLST_HD BLST_FORCEINLINE const Fp2& get_BLST_sqrt_minus_1(void)
    {
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
        return BLST_sqrt_minus_1_DEV;
#else
        return BLST_sqrt_minus_1_HOST;
#endif
    }

    BLST_HD BLST_FORCEINLINE const Fp2& get_BLST_sqrt_sqrt_minus_1(void)
    {
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
        return BLST_sqrt_sqrt_minus_1_DEV;
#else
        return BLST_sqrt_sqrt_minus_1_HOST;
#endif
    }

    BLST_HD BLST_FORCEINLINE const Fp2& get_BLST_sqrt_minus_sqrt_minus_1(void)
    {
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
        return BLST_sqrt_minus_sqrt_minus_1_DEV;
#else
        return BLST_sqrt_minus_sqrt_minus_1_HOST;
#endif
    }








    /**
     * @brief Computes the sign and parity of a 384-bit complex number modulo a given modulus.
     *
     * This function determines the sign and parity of a complex number `a` modulo `p`.
     * The sign is determined based on the imaginary part (`a[1]`) if it is non-zero;
     * otherwise, the real part (`a[0]`) is used. The parity is determined based on
     * the least significant bit of the selected part.
     *
     * @param[in]  a The input complex number represented as two 384-bit integers.
     * @param[in]  p The modulus.
     *
     * @return A 2-bit value:
     *         - Bit 1: Sign (0 if non-negative, 1 if negative).
     *         - Bit 0: Parity (0 if even, 1 if odd).
     *
     * @note Assumes `a` and `p` are valid 384-bit integers.
     */
    BLST_HD BLST_INLINE limb_t cc_sgn0_pty_mod_384x(const Fp2& a, const Fp& p)
    {
        limb_t re, im, sign, prty;
        const Fp& a0 = a[0];
        const Fp& a1 = a[1];

        re = cc_sgn0_pty_mod_n<NLIMBS(384)>(a0, p.data);
        im = cc_sgn0_pty_mod_n<NLIMBS(384)>(a1, p.data);

        /* a->im!=0 ? sgn0(a->im) : sgn0(a->re) */
        sign = (limb_t)0 - vec_is_zero(a[1].data, sizeof(Fp));
        sign = (re & sign) | (im & ~sign);

        /* a->re==0 ? prty(a->im) : prty(a->re) */
        prty = (limb_t)0 - vec_is_zero(a[0].data, sizeof(Fp));
        prty = (im & prty) | (re & ~prty);

        return (sign & 2) | (prty & 1);
    }


    /**
     * @brief Computes the sign and parity of a 384-bit complex number in Montgomery form modulo a given modulus.
     *
     * This function converts the complex number `a` from Montgomery form to standard form
     * and then determines its sign and parity modulo `p`. The sign and parity are computed
     * in the same way as `cc_sgn0_pty_mod_384x`.
     *
     * @param[in]  a   The input complex number in Montgomery form.
     * @param[in]  p   The modulus.
     * @param[in]  n0  The Montgomery constant.
     *
     * @return A 2-bit value:
     *         - Bit 1: Sign (0 if non-negative, 1 if negative).
     *         - Bit 0: Parity (0 if even, 1 if odd).
     *
     * @note Assumes `a` and `p` are valid 384-bit integers.
     */
    BLST_HD BLST_INLINE limb_t cc_sgn0_pty_mont_384x(const Fp2& a, const Fp& p, limb_t n0)
    {
        Fp2 tmp;

        cc_from_mont_n<NLIMBS(384)>(tmp[0].data, a[0].data, p.data, n0);
        cc_from_mont_n<NLIMBS(384)>(tmp[1].data, a[1].data, p.data, n0);

        return cc_sgn0_pty_mod_384x(tmp, p);
    }

    /**
     * @brief Performs Montgomery multiplication on two 384-bit complex numbers.
     *
     * Computes `(a * b * R^(-1)) mod p` for complex numbers `a` and `b` in Montgomery form.
     * The result is stored in `ret`.
     */
    BLST_HD BLST_INLINE void cc_mul_mont_384x(Fp2& ret, const Fp2& a, const Fp2& b,
        const Fp& p, limb_t n0)
    {
        Fp aa, bb, cc;

        cc_add_mod_n<NLIMBS(384)>(aa, a[0], a[1], p);
        cc_add_mod_n<NLIMBS(384)>(bb, b[0], b[1], p);
        cc_mul_mont_n<NLIMBS(384)>(bb, bb, aa, p, n0);
        cc_mul_mont_n<NLIMBS(384)>(aa, a[0], b[0], p, n0);
        cc_mul_mont_n<NLIMBS(384)>(cc, a[1], b[1], p, n0);
        cc_sub_mod_n<NLIMBS(384)>(ret[0], aa, cc, p);
        cc_sub_mod_n<NLIMBS(384)>(ret[1], bb, aa, p);
        cc_sub_mod_n<NLIMBS(384)>(ret[1], ret[1], cc, p);
    }


    /**
     * @brief Squares a 384-bit complex number in Montgomery form modulo a given modulus.
     *
     * Computes `(a^2 * R^(-1)) mod p` for a complex number `a` in Montgomery form.
     * The result is stored in `ret`.
     *
     * @param[out] ret The output complex number in Montgomery form.
     * @param[in]  a   The input complex number in Montgomery form.
     * @param[in]  p   The modulus.
     * @param[in]  n0  The Montgomery constant.
     *
     * @note Assumes `a` and `p` are valid 384-bit integers.
     */
    BLST_HD BLST_INLINE void cc_sqr_mont_382x(Fp2& ret, const Fp2& a,
        const Fp& p, limb_t n0)
    {
        wide_limb_t wide;
        limb_t mask, carry, borrow;
        size_t i;
        Fp t0, t1;

        /* "add_mod_n(t0, a[0], a[1], p, NLIMBS(384));" */
        for (carry = 0, i = 0; i < NLIMBS(384); i++) {
            cc_add_carry_wide_limb_limb(&wide, a[0][i], a[1][i], carry);
            t0[i] = wide.lo;
            carry = wide.hi;
        }

        /* "sub_mod_n(t1, a[0], a[1], p, NLIMBS(384));" */
        for (borrow = 0, i = 0; i < NLIMBS(384); i++) {
            cc_sub_borrow_wide_limb_limb(&wide, a[0][i], a[1][i], borrow);
            t1[i] = wide.lo;
            borrow = wide.hi & 1;
        }
        mask = 0 - borrow;
        BLST_launder(mask);

        /* "mul_mont_n(ret[1], a[0], a[1], p, n0, NLIMBS(384));" */
        cc_mul_mont_nonred_n<NLIMBS(384)>(ret[1], a[0], a[1], p, n0);

        /* "add_mod_n(ret[1], ret[1], ret[1], p, NLIMBS(384));" */
        for (carry = 0, i = 0; i < NLIMBS(384); i++) {
            limb_t a_i = ret[1][i];
            ret[1][i] = a_i << 1 | carry;
            carry = a_i >> (BLST_LIMB_SIZE - 1);
        }

        /* "mul_mont_n(ret[0], t0, t1, p, n0, NLIMBS(384));" */
        cc_mul_mont_nonred_n<NLIMBS(384)>(ret[0], t0, t1, p, n0);

        /* account for t1's sign... */
        for (borrow = 0, i = 0; i < NLIMBS(384); i++) {
            cc_sub_borrow_wide_limb_limb(&wide, ret[0][i], t0[i] & mask, borrow);
            ret[0][i] = wide.lo;
            borrow = wide.hi & 1;
        }
        mask = 0 - borrow;
        BLST_launder(mask);
        for (carry = 0, i = 0; i < NLIMBS(384); i++) {
            cc_add_carry_wide_limb_limb(&wide, ret[0][i], p[i] & mask, carry);
            ret[0][i] = wide.lo;
            carry = wide.hi;
        }
    }

    BLST_HD BLST_INLINE  void cc_mul_by_8_mod_384x(Fp2& ret, const Fp2& a, const Fp& mod)
    {
        Fp::mul_by_8(ret[0], a[0]);
        Fp::mul_by_8(ret[1], a[1]);
    }

    BLST_HD BLST_INLINE  void cc_mul_by_1_plus_i_mod_384x(Fp2& ret, const Fp2& a,
        const Fp& mod)
    {
        Fp t;

        Fp::add(t, a[0], a[1]);
        Fp::sub(ret[0], a[0], a[1]);
        vec_copy(ret[1], t, sizeof(t));
    }
    BLST_HD BLST_INLINE void cc_add_mod_384x(Fp2& ret, const Fp2& a, const Fp2& b,
        const Fp& mod)
    {
        Fp::add(ret[0], a[0], b[0]);
        Fp::add(ret[1], a[1], b[1]);
    }
    BLST_HD BLST_INLINE void cc_sub_mod_384x(Fp2& ret, const Fp2& a, const Fp2& b,
        const Fp& mod)
    {
        Fp::sub(ret[0], a[0], b[0]);
        Fp::sub(ret[1], a[1], b[1]);
    }
    BLST_HD BLST_INLINE void cc_lshift_mod_384x(Fp2& ret, const Fp2& a, size_t n,
        const Fp& mod)
    {
        Fp::lshift(ret[0], a[0], n);
        Fp::lshift(ret[1], a[1], n);
    }

    /* @brief Squares a 384-bit complex number in Montgomery form modulo a given modulus.
     *
     * Computes `(a^2 * R^(-1)) mod p` for a complex number `a` in Montgomery form.
     * The result is stored in `ret`.
     *
     * @param[out] ret The output complex number in Montgomery form.
     * @param[in]  a   The input complex number in Montgomery form.
     * @param[in]  mod The modulus.
     * @param[in]  n0  The Montgomery constant.
     *
     * @note Assumes `a` and `mod` are valid 384-bit integers.
     */
    BLST_HD BLST_INLINE void cc_sqr_mont_384x(Fp2& ret, const Fp2& a, const Fp& mod, limb_t n0)
    {
        Fp t0, t1;

        Fp::add(t0, a[0], a[1]);
        Fp::sub(t1, a[0], a[1]);

        Fp::mul(ret[1], a[0], a[1]);
        Fp::add(ret[1], ret[1], ret[1]);

        Fp::mul(ret[0], t0, t1);
    }
    //#endif


    //BLST_HD BLST_FORCEINLINE void cneg_mod_v384x(Fp2& ret, const Fp2& a, bool_t flag,
    //    const Fp& p)
    //{
    //    //Fp::cneg(ret[0], a[0], flag/*, p*/);
    //    //Fp::cneg(ret[1], a[1], flag/*, p*/);
    //}


    BLST_HD BLST_INLINE void Fp2::exp(Fp2& out, const Fp2& inp, const byte* pow,
        size_t pow_bits)
    {
        Fp2 ret = inp;

        //vec_copy(ret, inp, sizeof(ret));  /* |ret| = |inp|^1 */
        --pow_bits; /* most significant bit is accounted for, skip over */
        while (pow_bits--) {
            Fp2::sqr(ret, ret);
            if (is_bit_set(pow, pow_bits))
                Fp2::mul(ret, ret, inp/*, p, n0*/);
        }
        out = ret;
        //vec_copy(out, ret, sizeof(ret));  /* |out| = |ret| */
    }

    BLST_HD BLST_INLINE void cc_mul_by_3_mod_384x(Fp2& ret, const Fp2& a, const Fp& mod)
    {
        Fp::mul_by_3(ret[0], a[0]);
        Fp::mul_by_3(ret[1], a[1]);
    }


    BLST_HD BLST_FORCEINLINE limb_t Fp2::sgn0_pty_mont(const Fp2& a)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        return sgn0x_pty_mont_384x(&a, &Fp::modulus(), BLST_p0);
#else
        return cc_sgn0_pty_mont_384x(a, Fp::modulus(), BLST_p0);
#endif
    }
    BLST_HD BLST_FORCEINLINE limb_t Fp2::sgn0_pty(const Fp2& a)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        return sgn0_pty_mod_384x(&a, &Fp::modulus());
#else
        return cc_sgn0_pty_mod_384x(a, Fp::modulus());
#endif
    }


#if defined(BLST_ASM) && !defined(__CUDA_ARCH__)

    BLST_HD BLST_FORCEINLINE void mul_v382x(UFp2& ret, const Fp2& a, const Fp2& b, const Fp& p)
    {
#ifdef BLST_ADX
        mulx_382x(&ret, &a, &b, &p);
#else
        mul_382x(&ret, &a, &b, &p);
#endif
    }
    BLST_HD BLST_FORCEINLINE void sqr_v382x(UFp2& ret, const Fp2& a, const Fp& p)
    {
#ifdef BLST_ADX
        sqrx_382x(&ret, &a, &p);
#else
        sqr_382x(&ret, &a, &p);
#endif
    }

#endif

    BLST_HD BLST_FORCEINLINE void Fp2::lshift(Fp2& ret, const Fp2& a, size_t count)
    {
        Fp::lshift(ret[0], a[0], count);
		Fp::lshift(ret[1], a[1], count);
    }

#define vec_load_global vec_copy


    BLST_HD BLST_FORCEINLINE void Fp2::conjugate(Fp2& ret, const Fp2& a)
    {
        vec_copy(ret[0], a[0], sizeof(Fp));  /* Re(ret) = Re(a) */
        Fp::cneg(ret[1], a[1], 1);               /* Im(ret) = -Im(a) */
    }


    BLST_HD BLST_INLINE bool_t Fp2::sqrt_align(Fp2& out, const Fp2& ret,
        const Fp2& sqrt, const Fp2& inp)
    {
        Fp2 coeff;
        bool_t is_sqrt;

        {
            Fp2 t0;
            bool_t flag;
            /*
             * Instead of multiple trial squarings we can perform just one
             * and see if the result is "rotated by multiple of 90°" in
             * relation to |inp|, and "rotate" |ret| accordingly.
             */
            Fp2::sqr(t0, sqrt);
            /* "sqrt(|inp|)"^2 = (a + b*i)^2 = (a^2-b^2) + 2ab*i */

            /* (a^2-b^2) + 2ab*i == |inp| ? |ret| is spot on */
            Fp2 t1;
            Fp2::sub(t1, t0, inp);
            is_sqrt = vec_is_zero(&t1, sizeof(t1));
            vec_copy(&coeff, &get_BLST_Rx().p2, sizeof(coeff));

            /* -(a^2-b^2) - 2ab*i == |inp| ? "rotate |ret| by 90°" */
            Fp2::add(t1, t0, inp);
            vec_select(&coeff, &get_BLST_sqrt_minus_1(), &coeff, sizeof(coeff),
                flag = vec_is_zero(&t1, sizeof(t1)));
            is_sqrt |= flag;

            /* 2ab - (a^2-b^2)*i == |inp| ? "rotate |ret| by 135°" */
            Fp::sub(t1[0], t0[0], inp[1]);
            Fp::add(t1[1], t0[1], inp[0]);
            vec_select(&coeff, &get_BLST_sqrt_sqrt_minus_1(), &coeff, sizeof(coeff),
                flag = vec_is_zero(&t1, sizeof(t1)));
            is_sqrt |= flag;

            /* -2ab + (a^2-b^2)*i == |inp| ? "rotate |ret| by 45°" */
            Fp::add(t1[0], t0[0], inp[1]);
            Fp::sub(t1[1], t0[1], inp[0]);
            vec_select(&coeff, &get_BLST_sqrt_minus_sqrt_minus_1(), &coeff, sizeof(coeff),
                flag = vec_is_zero(&t1, sizeof(t1)));
            is_sqrt |= flag;
        }

        /* actual "rotation" */
        Fp2::mul(out, ret, coeff);

        return is_sqrt;
    }

    /*
     * |inp| = a + b*i
     */
    BLST_HD BLST_INLINE bool_t Fp2::recip_sqrt(Fp2& out, const Fp2& inp)
    {

        const Fp2 recip_ZZZ = {      /* 1/(Z^3) */
          Fp{ TO_LIMB_T(0x65018f5c28f598eb), TO_LIMB_T(0xe6020417f022d916),
            TO_LIMB_T(0xd6327313288369c7), TO_LIMB_T(0x622ded8eb447156f),
            TO_LIMB_T(0xe52a2aee72c2a01f), TO_LIMB_T(0x089812fb8481ffe4) },
          Fp{ TO_LIMB_T(0x2574eb851eb8619f), TO_LIMB_T(0xdba2e97912925604),
            TO_LIMB_T(0x67e495a909e7a18e), TO_LIMB_T(0xdf2da23b8145b8f7),
            TO_LIMB_T(0xcf5d3728310ebf6d), TO_LIMB_T(0x11be446236f4c116) }
        };
        const Fp2 magic_ZZZ = {      /* 1/Z^3 = a + b*i */
            /* a^2 + b^2 */
    Fp{ TO_LIMB_T(0xaa7eb851eb8508e0), TO_LIMB_T(0x1c54fdf360989374),
      TO_LIMB_T(0xc87f2fc6e716c62e), TO_LIMB_T(0x0124aefb1f9efea7),
      TO_LIMB_T(0xb2f8be63e844865c), TO_LIMB_T(0x08b47f775a7ef35a) },
      /* (a^2 + b^2)^((P-3)/4) */
    Fp{ TO_LIMB_T(0xe4132bbd838cf70a), TO_LIMB_T(0x01d769ac83772c19),
      TO_LIMB_T(0xa83dd6e974c22e45), TO_LIMB_T(0xbc8ec3e777b08dff),
      TO_LIMB_T(0xc035c2042ecf5da3), TO_LIMB_T(0x073929e97f0850bf) }
        };

        Fp aa, bb, cc;
        Fp2 inp_;
        bool_t is_sqrt;

        Fp::sqr(aa, inp[0]);
        Fp::sqr(bb, inp[1]);
        Fp::add(aa, aa, bb);

        is_sqrt = Fp::recip_sqrt(cc, aa);  /* 1/sqrt(a²+b²)                    */

        /* if |inp| doesn't have quadratic residue, multiply by "1/Z³" ...    */
        Fp2::mul(inp_, inp, recip_ZZZ);
        /* ... and adjust |aa| and |cc| accordingly                           */
        {
            Fp za, zc;

            Fp::mul(za, aa, magic_ZZZ[0]); /* aa*(za² + zb²)                   */
            Fp::mul(zc, cc, magic_ZZZ[1]); /* cc*(za² + zb²)^((p-3)/4)         */
            vec_select(aa, aa, za, sizeof(aa), is_sqrt);
            vec_select(cc, cc, zc, sizeof(cc), is_sqrt);
        }
        vec_select(&inp_, &inp, &inp_, sizeof(inp_), is_sqrt);

        Fp::mul(aa, aa, cc);               /* sqrt(a²+b²)                      */

        Fp::sub(bb, inp_[0], aa);
        Fp::add(aa, inp_[0], aa);
        vec_select(aa, bb, aa, sizeof(aa), vec_is_zero(aa, sizeof(aa)));
        Fp::div_by_2(aa, aa);              /* (a ± sqrt(a²+b²))/2              */

        /* if it says "no sqrt," final "align" will find right one...         */
        (void)Fp::recip_sqrt(out[0], aa);  /* 1/sqrt((a ± sqrt(a²+b²))/2)      */

        Fp::div_by_2(out[1], inp_[1]);
        Fp::mul(out[1], out[1], out[0]);   /* b/(2*sqrt((a ± sqrt(a²+b²))/2))  */
        Fp::mul(out[0], out[0], aa);       /* sqrt((a ± sqrt(a²+b²))/2)        */

        /* bound to succeed                                                   */
        (void)Fp2::sqrt_align(out, out, out, inp_);

        Fp::mul(out[0], out[0], cc);       /* inverse the result               */
        Fp::mul(out[1], out[1], cc);
        Fp::neg(out[1], out[1]);

        return is_sqrt;
    }

    BLST_HD BLST_INLINE bool_t Fp2::sqrt(Fp2& out, const Fp2& inp)
    {
        Fp2 ret;

        {
            Fp aa, bb;

            Fp::sqr(aa, inp[0]);
            Fp::sqr(bb, inp[1]);
            Fp::add(aa, aa, bb);

            /* don't pay attention to return value, final "align" will tell...    */
            (void)Fp::sqrt(aa, aa);            /* sqrt(a²+b²)                      */

            Fp::sub(bb, inp[0], aa);
            Fp::add(aa, inp[0], aa);
            vec_select(aa, bb, aa, sizeof(aa), vec_is_zero(aa, sizeof(aa)));
            Fp::div_by_2(aa, aa);              /* (a ± sqrt(a²+b²))/2              */

            /* if it says "no sqrt," final "align" will find right one...         */
            (void)Fp::recip_sqrt(ret[0], aa);  /* 1/sqrt((a ± sqrt(a²+b²))/2)      */

            Fp::div_by_2(ret[1], inp[1]);
            Fp::mul(ret[1], ret[1], ret[0]);   /* b/(2*sqrt((a ± sqrt(a²+b²))/2))  */
            Fp::mul(ret[0], ret[0], aa);       /* sqrt((a ± sqrt(a²+b²))/2)        */
        }

        /*
         * Now see if |ret| is or can be made sqrt(|inp|)...
         */

        return Fp2::sqrt_align(out, ret, ret, inp);
    }

    BLST_HD BLST_FORCEINLINE int Fp2::is_square(const Fp2& inp)
    {
        Fp aa, bb;

        Fp::sqr(aa, inp[0]);
        Fp::sqr(bb, inp[1]);
        Fp::add(aa, aa, bb);

        return (int)Fp::is_sqrt(aa);
    }



    /*
     * Fp2  = Fp[u]  / (u^2 + 1)
     * Fp6  = Fp2[v] / (v^3 - u - 1)
     * Fp12 = Fp6[w] / (w^2 - v)
     */

    BLST_HD BLST_INLINE void Fp2::mul_by_u_plus_1(Fp2& ret, const Fp2& a)
    {
        //mul_by_1_plus_i_mod_v384x(ret, a, get_BLST_P());

#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        mul_by_1_plus_i_mod_384x(&ret, &a, &get_BLST_P());
#else
        cc_mul_by_1_plus_i_mod_384x(ret, a, get_BLST_P());
#endif
    }

    /*
     * caveat lector! |n| has to be non-zero and not more than 3!
     */
    BLST_HD BLST_INLINE void frobenius_map_fp2(Fp2& ret, const Fp2& a, size_t n)
    {
        vec_copy(ret[0], a[0], sizeof(ret[0]));
        Fp::cneg(ret[1], a[1], n & 1);
    }


    BLST_HD BLST_INLINE void Fp2::add(Fp2& ret, const Fp2& a, const Fp2& b)
    {
        //add_mod_v384x(ret, a, b, get_BLST_P());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        add_mod_384x(&ret, &a, &b, &get_BLST_P());
#else
        cc_add_mod_384x(ret, a, b, get_BLST_P());
#endif
    }
    BLST_HD BLST_INLINE void Fp2::sub(Fp2& ret, const Fp2& a, const Fp2& b)
    {
        //sub_mod_v384x(ret, a, b, get_BLST_P());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        sub_mod_384x(&ret, &a, &b, &get_BLST_P());
#else
        cc_sub_mod_384x(ret, a, b, get_BLST_P());
#endif
    }
    BLST_HD BLST_INLINE void Fp2::mul(Fp2& ret, const Fp2& a, const Fp2& b)
    {
        //mul_mont_v384x(ret, a, b, get_BLST_P(), BLST_p0);
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        mulx_mont_384x(&ret, &a, &b, &get_BLST_P(), BLST_p0);
#else
        cc_mul_mont_384x(ret, a, b, get_BLST_P(), BLST_p0);
#endif
    }

    BLST_HD BLST_INLINE void Fp2::sqr(Fp2& ret, const Fp2& a)
    {
        //sqr_mont_v384x(ret, a, get_BLST_P(), BLST_p0);
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        sqrx_mont_384x(&ret, &a, &get_BLST_P(), BLST_p0);
#else
        cc_sqr_mont_384x(ret, a, get_BLST_P(), BLST_p0);
#endif
    }
    BLST_HD BLST_INLINE void Fp2::cneg(Fp2& ret, const Fp2& a, bool_t flag)
    {
        //cneg_mod_v384x(ret, a, flag, get_BLST_P());
        Fp::cneg(ret[0], a[0], flag/*, p*/);
        Fp::cneg(ret[1], a[1], flag/*, p*/);
    }
    BLST_HD BLST_INLINE void Fp2::neg(Fp2& ret, const Fp2& a)
    {
        //cneg_mod_v384x(ret, a, 1, get_BLST_P());
        cneg(ret, a, 1);
    }
    //  BLST_HD void Fp2::from_mont(Fp2& ret, const Fp2& a)
    //  {
          //from_mont_v384x(ret, a, get_BLST_P(), BLST_p0);
    //  }
      //BLST_HD void Fp2::to_mont(Fp2& ret, const Fp2& a);
    BLST_HD BLST_INLINE void Fp2::reciprocal(Fp2& out, const Fp2& inp)
    {
        //reciprocal_fp2(out, inp);
        Fp t0, t1;

        /*
         * |out| = 1/(a + b*i) = a/(a^2+b^2) - b/(a^2+b^2)*i
         */
        Fp::sqr(t0, inp[0]);
        Fp::sqr(t1, inp[1]);
        Fp::add(t0, t0, t1);
        Fp::reciprocal(t1, t0);
        Fp::mul(out[0], inp[0], t1);
        Fp::mul(out[1], inp[1], t1);
        Fp::neg(out[1], out[1]);
    }


    BLST_HD BLST_INLINE void Fp2::flt_reciprocal(Fp2& out, const Fp2& inp)
    {
        Fp t0, t1;

        /*
         * |out| = 1/(a + b*i) = a/(a^2+b^2) - b/(a^2+b^2)*i
         */
        Fp::sqr(t0, inp[0]);
        Fp::sqr(t1, inp[1]);
        Fp::add(t0, t0, t1);
        Fp::flt_reciprocal(t1, t0);
        Fp::mul(out[0], inp[0], t1);
        Fp::mul(out[1], inp[1], t1);
        Fp::neg(out[1], out[1]);
    }



    //BLST_HD BLST_INLINE bool_t Fp2::sqrt(Fp2& out, const Fp2& inp)
    //{
    //    return sqrt_fp2(out, inp);
    //}
    BLST_HD BLST_INLINE void Fp2::mul_by_3(Fp2& ret, const Fp2& a)
    {
        //mul_by_3_mod_v384x(ret, a, get_BLST_P());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        mul_by_3_mod_384x(&ret, &a, &get_BLST_P());
#else
        cc_mul_by_3_mod_384x(ret, a, get_BLST_P());
#endif
    }
    BLST_HD BLST_INLINE void Fp2::mul_by_8(Fp2& ret, const Fp2& a)
    {
        //mul_by_8_mod_v384x(ret, a, get_BLST_P());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        mul_by_8_mod_384x(&ret, &a, &get_BLST_P());
#else
        cc_mul_by_8_mod_384x(ret, a, get_BLST_P());
#endif
    }

    BLST_HD BLST_INLINE const Fp2& Fp2::one()
    {
		return *(const Fp2*)get_BLST_Rx().p2;
    }


    //BLST_EXTERN_C_END
}

#endif /* __BLS12_381_ASM_FIELDS_H__ */

