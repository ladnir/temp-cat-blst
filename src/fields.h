/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_FIELDS_H__
#define __BLS12_381_ASM_FIELDS_H__

#include "vect.h"
#include "consts.h"
#include "config.h"
#include "defines.h"
#include "backend.h"
#include "recip-addchain.h"
#include "exp.h"

BLST_EXTERN_C_BEGIN

/*
 * BLS12-381-specific Fp shortcuts to assembly.
 */
HOST_DEVICE_FUNC BLST_FORCEINLINE void add_fp(vec384 ret, const vec384 a, const vec384 b)
{   add_mod_v384(ret, a, b, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void sub_fp(vec384 ret, const vec384 a, const vec384 b)
{   sub_mod_v384(ret, a, b, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_by_3_fp(vec384 ret, const vec384 a)
{   mul_by_3_mod_v384(ret, a, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_by_8_fp(vec384 ret, const vec384 a)
{   mul_by_8_mod_v384(ret, a, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void lshift_fp(vec384 ret, const vec384 a, size_t count)
{   lshift_mod_v384(ret, a, count, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void rshift_fp(vec384 ret, const vec384 a, size_t count)
{   rshift_mod_v384(ret, a, count, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void div_by_2_fp(vec384 ret, const vec384 a)
{   div_by_2_mod_v384(ret, a, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_fp(vec384 ret, const vec384 a, const vec384 b)
{   mul_mont_v384(ret, a, b, *get_BLST_P(), BLST_p0);   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void sqr_fp(vec384 ret, const vec384 a)
{   sqr_mont_v384(ret, a, *get_BLST_P(), BLST_p0);   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void cneg_fp(vec384 ret, const vec384 a, bool_t flag)
{   cneg_mod_v384(ret, a, flag, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void from_fp(vec384 ret, const vec384 a)
{   from_mont_v384(ret, a, *get_BLST_P(), BLST_p0);   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void to_fp(vec384 ret, const vec384 a)
{ mul_fp(ret, a, *get_BLST_RR()); }

HOST_DEVICE_FUNC BLST_FORCEINLINE void redc_fp(vec384 ret, const vec768 a)
{   redc_mont_v384(ret, a, *get_BLST_P(), BLST_p0);   }

/*
 * BLS12-381-specific Fp2 shortcuts to assembly.
 */
HOST_DEVICE_FUNC BLST_FORCEINLINE void add_fp2(vec384x ret, const vec384x a, const vec384x b)
{   add_mod_v384x(ret, a, b, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void sub_fp2(vec384x ret, const vec384x a, const vec384x b)
{   sub_mod_v384x(ret, a, b, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_by_3_fp2(vec384x ret, const vec384x a)
{   mul_by_3_mod_v384x(ret, a, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_by_8_fp2(vec384x ret, const vec384x a)
{   mul_by_8_mod_v384x(ret, a, *get_BLST_P());   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void lshift_fp2(vec384x ret, const vec384x a, size_t count)
{
    lshift_mod_v384(ret[0], a[0], count, *get_BLST_P());
    lshift_mod_v384(ret[1], a[1], count, *get_BLST_P());
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void mul_fp2(vec384x ret, const vec384x a, const vec384x b)
{   mul_mont_v384x(ret, a, b, *get_BLST_P(), BLST_p0);   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void sqr_fp2(vec384x ret, const vec384x a)
{   sqr_mont_v384x(ret, a, *get_BLST_P(), BLST_p0);   }

HOST_DEVICE_FUNC BLST_FORCEINLINE void cneg_fp2(vec384x ret, const vec384x a, bool_t flag)
{
    cneg_mod_v384(ret[0], a[0], flag, *get_BLST_P());
    cneg_mod_v384(ret[1], a[1], flag, *get_BLST_P());
}

#define vec_load_global vec_copy

HOST_DEVICE_FUNC void reciprocal_fp(vec384 out, const vec384 inp);
HOST_DEVICE_FUNC void flt_reciprocal_fp(vec384 out, const vec384 inp);
HOST_DEVICE_FUNC bool_t recip_sqrt_fp(vec384 out, const vec384 inp);
HOST_DEVICE_FUNC bool_t sqrt_fp(vec384 out, const vec384 inp);

HOST_DEVICE_FUNC void reciprocal_fp2(vec384x out, const vec384x inp);
HOST_DEVICE_FUNC void flt_reciprocal_fp2(vec384x out, const vec384x inp);
HOST_DEVICE_FUNC bool_t recip_sqrt_fp2(vec384x out, const vec384x inp,
                             const vec384x recip_ZZZ, const vec384x magic_ZZZ);
HOST_DEVICE_FUNC bool_t sqrt_fp2(vec384x out, const vec384x inp);
HOST_DEVICE_FUNC bool_t sqrt_align_fp2(vec384x out, const vec384x ret,
                             const vec384x sqrt, const vec384x inp);

typedef vec384x   vec384fp2;
typedef vec384fp2 vec384fp6[3];
typedef vec384fp6 vec384fp12[2];

HOST_DEVICE_FUNC void sqr_fp12(vec384fp12 ret, const vec384fp12 a);
HOST_DEVICE_FUNC void cyclotomic_sqr_fp12(vec384fp12 ret, const vec384fp12 a);
HOST_DEVICE_FUNC void mul_fp12(vec384fp12 ret, const vec384fp12 a, const vec384fp12 b);
HOST_DEVICE_FUNC void mul_by_xy00z0_fp12(vec384fp12 ret, const vec384fp12 a,
                                               const vec384fp6 xy00z0);
HOST_DEVICE_FUNC BLST_FORCEINLINE void conjugate_fp2(vec384x ret, const vec384x a)
{
    vec_copy(ret[0], a[0], sizeof(vec384));  /* Re(ret) = Re(a) */
    cneg_fp(ret[1], a[1], 1);               /* Im(ret) = -Im(a) */
}

HOST_DEVICE_FUNC void conjugate_fp12(vec384fp12 a);
HOST_DEVICE_FUNC void inverse_fp12(vec384fp12 ret, const vec384fp12 a);
/* caveat lector! |n| has to be non-zero and not more than 3! */
HOST_DEVICE_FUNC void frobenius_map_fp12(vec384fp12 ret, const vec384fp12 a, size_t n);

#define neg_fp(r,a) cneg_fp((r),(a),1)
#define neg_fp2(r,a) cneg_fp2((r),(a),1)




//// Vectorized addition of Fp elements (CPU)
//BLST_INLINE void add_fp_vec_cc(vec384 *out, const vec384 *a, const vec384 *b, size_t num)
//{
//    while (num--) {
//        add_fp(*out, *a, *b);
//        out++; a++; b++;
//    }
//}
//
//// kernel for addition of Fp elements (GPU)
//GLOBAL_FUNC BLST_INLINE void add_fp_vec_ker(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//    if (idx < num) {
//        add_fp(out[idx], a[idx], b[idx]);
//	}
//}
//
//// Vectorized addition of Fp elements (GPU)
//BLST_INLINE void add_fp_vec_cu(vec384 *out, const vec384 *a, const vec384 *b, size_t num)
//{
//    size_t threads_per_block = 256;
//    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
//    add_fp_vec_ker<<<num_blocks, threads_per_block>>>(out, a, b, num);
//	cudaDeviceSynchronize();
//}
//
//
//BLST_INLINE void add_fp_vec(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//#ifdef BLST_CUDA
//	add_fp_vec_cu(out, a, b, num);
//#else
//	add_fp_vec_cc(out, a, b, num);
//#endif
//}


































#if defined(__OPTIMIZE_SIZE__)  || defined(__CUDA_ARCH__)
HOST_DEVICE_FUNC BLST_FORCEINLINE void recip_sqrt_fp_3mod4(vec384 out, const vec384 inp)
{
    const byte BLS_12_381_P_minus_3_div_4[] = {
        TO_BYTES(0xee7fbfffffffeaaa), TO_BYTES(0x07aaffffac54ffff),
        TO_BYTES(0xd9cc34a83dac3d89), TO_BYTES(0xd91dd2e13ce144af),
        TO_BYTES(0x92c6e9ed90d2eb35), TO_BYTES(0x0680447a8e5ff9a6)
    };

    exp_mont_384(out, inp, BLS_12_381_P_minus_3_div_4, 379, *get_BLST_P(), BLST_p0);
}
#else
# if 1
/*
 * "383"-bit variant omits full reductions at the ends of squarings,
 * which results in up to ~15% improvement. [One can improve further
 * by omitting full reductions even after multiplications and
 * performing final reduction at the very end of the chain.]
 */
HOST_DEVICE_FUNC BLST_FORCEINLINE void sqr_n_mul_fp(vec384 out, const vec384 a, size_t count,
    const vec384 b)
{
    sqr_n_mul_mont_v383(out, a, count, *get_BLST_P(), BLST_p0, b);
}
# else
HOST_DEVICE_FUNC void sqr_n_mul_fp(vec384 out, const vec384 a, size_t count,
    const vec384 b)
{
    while (count--) {
        sqr_fp(out, a);
        a = out;
    }
    mul_fp(out, out, b);
}
# endif

# define sqr(ret,a)		sqr_fp(ret,a)
# define mul(ret,a,b)		mul_fp(ret,a,b)
# define sqr_n_mul(ret,a,n,b)	sqr_n_mul_fp(ret,a,n,b)

# include "sqrt-addchain.h"
HOST_DEVICE_FUNC BLST_INLINE void recip_sqrt_fp_3mod4(vec384 out, const vec384 inp)
{
    RECIP_SQRT_MOD_BLS12_381_P(out, inp, vec384);
}
# undef RECIP_SQRT_MOD_BLS12_381_P

# undef sqr_n_mul
# undef sqr
# undef mul
#endif

HOST_DEVICE_FUNC BLST_INLINE bool_t recip_sqrt_fp(vec384 out, const vec384 inp)
{
    vec384 t0, t1;
    bool_t ret;

    recip_sqrt_fp_3mod4(t0, inp);

    mul_fp(t1, t0, inp);
    sqr_fp(t1, t1);
    ret = vec_is_equal(t1, inp, sizeof(t1));
    vec_copy(out, t0, sizeof(t0));

    return ret;
}

HOST_DEVICE_FUNC BLST_INLINE bool_t sqrt_fp(vec384 out, const vec384 inp)
{
    vec384 t0, t1;
    bool_t ret;

    recip_sqrt_fp_3mod4(t0, inp);

    mul_fp(t0, t0, inp);
    sqr_fp(t1, t0);
    ret = vec_is_equal(t1, inp, sizeof(t1));
    vec_copy(out, t0, sizeof(t0));

    return ret;
}

HOST_DEVICE_FUNC BLST_FORCEINLINE int blst_fp_sqrt(vec384 out, const vec384 inp)
{
    return (int)sqrt_fp(out, inp);
}

HOST_DEVICE_FUNC BLST_FORCEINLINE int blst_fp_is_square(const vec384 inp)
{
    return (int)ct_is_square_mod_v384(inp, *get_BLST_P());
}

HOST_DEVICE_FUNC BLST_INLINE bool_t sqrt_align_fp2(vec384x out, const vec384x ret,
    const vec384x sqrt, const vec384x inp)
{
    const vec384x sqrt_minus_1 = { { 0 }, { ONE_MONT_P_LIMBS } };
    const vec384x sqrt_sqrt_minus_1 = {
        /*
         * "magic" number is ±2^((p-3)/4)%p, which is "1/sqrt(2)",
         * in quotes because 2*"1/sqrt(2)"^2 == -1 mod p, not 1,
         * but it pivots into "complex" plane nevertheless...
         */
        { TO_LIMB_T(0x3e2f585da55c9ad1), TO_LIMB_T(0x4294213d86c18183),
          TO_LIMB_T(0x382844c88b623732), TO_LIMB_T(0x92ad2afd19103e18),
          TO_LIMB_T(0x1d794e4fac7cf0b9), TO_LIMB_T(0x0bd592fc7d825ec8) },
        { TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),
          TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),
          TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2) }
    };
    const vec384x sqrt_minus_sqrt_minus_1 = {
      { TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),
        TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),
        TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2) },
      { TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),
        TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),
        TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2) }
    };
    vec384x coeff, t0, t1;
    bool_t is_sqrt, flag;

    /*
     * Instead of multiple trial squarings we can perform just one
     * and see if the result is "rotated by multiple of 90°" in
     * relation to |inp|, and "rotate" |ret| accordingly.
     */
    sqr_fp2(t0, sqrt);
    /* "sqrt(|inp|)"^2 = (a + b*i)^2 = (a^2-b^2) + 2ab*i */

    /* (a^2-b^2) + 2ab*i == |inp| ? |ret| is spot on */
    sub_fp2(t1, t0, inp);
    is_sqrt = vec_is_zero(t1, sizeof(t1));
    vec_copy(coeff, get_BLST_Rx()->p2, sizeof(coeff));

    /* -(a^2-b^2) - 2ab*i == |inp| ? "rotate |ret| by 90°" */
    add_fp2(t1, t0, inp);
    vec_select(coeff, sqrt_minus_1, coeff, sizeof(coeff),
        flag = vec_is_zero(t1, sizeof(t1)));
    is_sqrt |= flag;

    /* 2ab - (a^2-b^2)*i == |inp| ? "rotate |ret| by 135°" */
    sub_fp(t1[0], t0[0], inp[1]);
    add_fp(t1[1], t0[1], inp[0]);
    vec_select(coeff, sqrt_sqrt_minus_1, coeff, sizeof(coeff),
        flag = vec_is_zero(t1, sizeof(t1)));
    is_sqrt |= flag;

    /* -2ab + (a^2-b^2)*i == |inp| ? "rotate |ret| by 45°" */
    add_fp(t1[0], t0[0], inp[1]);
    sub_fp(t1[1], t0[1], inp[0]);
    vec_select(coeff, sqrt_minus_sqrt_minus_1, coeff, sizeof(coeff),
        flag = vec_is_zero(t1, sizeof(t1)));
    is_sqrt |= flag;

    /* actual "rotation" */
    mul_fp2(out, ret, coeff);

    return is_sqrt;
}

/*
 * |inp| = a + b*i
 */
HOST_DEVICE_FUNC BLST_INLINE bool_t recip_sqrt_fp2(vec384x out, const vec384x inp,
    const vec384x recip_ZZZ,
    const vec384x magic_ZZZ)
{
    vec384 aa, bb, cc;
    vec384x inp_;
    bool_t is_sqrt;

    sqr_fp(aa, inp[0]);
    sqr_fp(bb, inp[1]);
    add_fp(aa, aa, bb);

    is_sqrt = recip_sqrt_fp(cc, aa);  /* 1/sqrt(a²+b²)                    */

    /* if |inp| doesn't have quadratic residue, multiply by "1/Z³" ...    */
    mul_fp2(inp_, inp, recip_ZZZ);
    /* ... and adjust |aa| and |cc| accordingly                           */
    {
        vec384 za, zc;

        mul_fp(za, aa, magic_ZZZ[0]); /* aa*(za² + zb²)                   */
        mul_fp(zc, cc, magic_ZZZ[1]); /* cc*(za² + zb²)^((p-3)/4)         */
        vec_select(aa, aa, za, sizeof(aa), is_sqrt);
        vec_select(cc, cc, zc, sizeof(cc), is_sqrt);
    }
    vec_select(inp_, inp, inp_, sizeof(inp_), is_sqrt);

    mul_fp(aa, aa, cc);               /* sqrt(a²+b²)                      */

    sub_fp(bb, inp_[0], aa);
    add_fp(aa, inp_[0], aa);
    vec_select(aa, bb, aa, sizeof(aa), vec_is_zero(aa, sizeof(aa)));
    div_by_2_fp(aa, aa);              /* (a ± sqrt(a²+b²))/2              */

    /* if it says "no sqrt," final "align" will find right one...         */
    (void)recip_sqrt_fp(out[0], aa);  /* 1/sqrt((a ± sqrt(a²+b²))/2)      */

    div_by_2_fp(out[1], inp_[1]);
    mul_fp(out[1], out[1], out[0]);   /* b/(2*sqrt((a ± sqrt(a²+b²))/2))  */
    mul_fp(out[0], out[0], aa);       /* sqrt((a ± sqrt(a²+b²))/2)        */

    /* bound to succeed                                                   */
    (void)sqrt_align_fp2(out, out, out, inp_);

    mul_fp(out[0], out[0], cc);       /* inverse the result               */
    mul_fp(out[1], out[1], cc);
    neg_fp(out[1], out[1]);

    return is_sqrt;
}

HOST_DEVICE_FUNC BLST_INLINE bool_t sqrt_fp2(vec384x out, const vec384x inp)
{
    vec384x ret;
    vec384 aa, bb;

    sqr_fp(aa, inp[0]);
    sqr_fp(bb, inp[1]);
    add_fp(aa, aa, bb);

    /* don't pay attention to return value, final "align" will tell...    */
    (void)sqrt_fp(aa, aa);            /* sqrt(a²+b²)                      */

    sub_fp(bb, inp[0], aa);
    add_fp(aa, inp[0], aa);
    vec_select(aa, bb, aa, sizeof(aa), vec_is_zero(aa, sizeof(aa)));
    div_by_2_fp(aa, aa);              /* (a ± sqrt(a²+b²))/2              */

    /* if it says "no sqrt," final "align" will find right one...         */
    (void)recip_sqrt_fp(ret[0], aa);  /* 1/sqrt((a ± sqrt(a²+b²))/2)      */

    div_by_2_fp(ret[1], inp[1]);
    mul_fp(ret[1], ret[1], ret[0]);   /* b/(2*sqrt((a ± sqrt(a²+b²))/2))  */
    mul_fp(ret[0], ret[0], aa);       /* sqrt((a ± sqrt(a²+b²))/2)        */

    /*
     * Now see if |ret| is or can be made sqrt(|inp|)...
     */

    return sqrt_align_fp2(out, ret, ret, inp);
}

HOST_DEVICE_FUNC BLST_FORCEINLINE int blst_fp2_sqrt(vec384x out, const vec384x inp)
{
    return (int)sqrt_fp2(out, inp);
}

HOST_DEVICE_FUNC BLST_FORCEINLINE int blst_fp2_is_square(const vec384x inp)
{
    vec384 aa, bb;

    sqr_fp(aa, inp[0]);
    sqr_fp(bb, inp[1]);
    add_fp(aa, aa, bb);

    return (int)ct_is_square_mod_v384(aa, *get_BLST_P());
}




#if defined(__OPTIMIZE_SIZE__) || defined(__CUDA_ARCH__)
/*
 * 608 multiplications for scalar inversion modulo BLS12-381 prime, 32%
 * more than corresponding optimal addition-chain, plus mispredicted
 * branch penalties on top of that... The addition chain below was
 * measured to be >50% faster.
 */
HOST_DEVICE_FUNC BLST_INLINE void flt_reciprocal_fp(vec384 out, const vec384 inp)
{
    const byte BLS12_381_P_minus_2[] = {
        TO_BYTES(0xb9feffffffffaaa9), TO_BYTES(0x1eabfffeb153ffff),
        TO_BYTES(0x6730d2a0f6b0f624), TO_BYTES(0x64774b84f38512bf),
        TO_BYTES(0x4b1ba7b6434bacd7), TO_BYTES(0x1a0111ea397fe69a)
    };

    exp_mont_384(out, inp, BLS12_381_P_minus_2, 381, *get_BLST_P(), BLST_p0);
}
#else
# define sqr(ret,a)		sqr_fp(ret,a)
# define mul(ret,a,b)		mul_fp(ret,a,b)
# define sqr_n_mul(ret,a,n,b)	sqr_n_mul_fp(ret,a,n,b)

HOST_DEVICE_FUNC BLST_INLINE void flt_reciprocal_fp(vec384 out, const vec384 inp)
{
    RECIPROCAL_MOD_BLS12_381_P(out, inp, vec384);
}
# undef RECIPROCAL_MOD_BLS12_381_P
# undef sqr_n_mul
# undef mul
# undef sqr
#endif

HOST_DEVICE_FUNC BLST_INLINE void flt_reciprocal_fp2(vec384x out, const vec384x inp)
{
    vec384 t0, t1;

    /*
     * |out| = 1/(a + b*i) = a/(a^2+b^2) - b/(a^2+b^2)*i
     */
    sqr_fp(t0, inp[0]);
    sqr_fp(t1, inp[1]);
    add_fp(t0, t0, t1);
    flt_reciprocal_fp(t1, t0);
    mul_fp(out[0], inp[0], t1);
    mul_fp(out[1], inp[1], t1);
    neg_fp(out[1], out[1]);
}

HOST_DEVICE_FUNC BLST_INLINE void reciprocal_fp(vec384 out, const vec384 inp)
{
    const vec384 Px8 = {    /* left-aligned value of the modulus */
        TO_LIMB_T(0xcff7fffffffd5558), TO_LIMB_T(0xf55ffff58a9ffffd),
        TO_LIMB_T(0x39869507b587b120), TO_LIMB_T(0x23ba5c279c2895fb),
        TO_LIMB_T(0x58dd3db21a5d66bb), TO_LIMB_T(0xd0088f51cbff34d2)
    };
    union { vec768 x; vec384 r[2]; } temp;

    ct_inverse_mod_v384(temp.x, inp, *get_BLST_P(), Px8);
    redc_mont_v384(temp.r[0], temp.x, *get_BLST_P(), BLST_p0);
    mul_mont_v384(temp.r[0], temp.r[0], *get_BLST_RR(), *get_BLST_P(), BLST_p0);

#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    /* sign goes straight to flt_reciprocal */
    mul_mont_v384(temp.r[1], temp.r[0], inp, *get_BLST_P(), BLST_p0);
    if (vec_is_equal(temp.r[1], get_BLST_Rx()->p, sizeof(vec384)) |
        vec_is_zero(temp.r[1], sizeof(vec384)))
        vec_copy(out, temp.r[0], sizeof(vec384));
    else
        flt_reciprocal_fp(out, inp);
#else
    vec_copy(out, temp.r[0], sizeof(vec384));
#endif
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void blst_fp_inverse(vec384 out, const vec384 inp)
{
    reciprocal_fp(out, inp);
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void blst_fp_eucl_inverse(vec384 ret, const vec384 a)
{
    reciprocal_fp(ret, a);
}

HOST_DEVICE_FUNC BLST_INLINE void reciprocal_fp2(vec384x out, const vec384x inp)
{
    vec384 t0, t1;

    /*
     * |out| = 1/(a + b*i) = a/(a^2+b^2) - b/(a^2+b^2)*i
     */
    sqr_fp(t0, inp[0]);
    sqr_fp(t1, inp[1]);
    add_fp(t0, t0, t1);
    reciprocal_fp(t1, t0);
    mul_fp(out[0], inp[0], t1);
    mul_fp(out[1], inp[1], t1);
    neg_fp(out[1], out[1]);
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void blst_fp2_inverse(vec384x out, const vec384x inp)
{
    reciprocal_fp2(out, inp);
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void blst_fp2_eucl_inverse(vec384x out, const vec384x inp)
{
    reciprocal_fp2(out, inp);
}

HOST_DEVICE_FUNC BLST_INLINE void reciprocal_fr(vec256 out, const vec256 inp)
{
    const vec256 rx2 = { /* left-aligned value of the modulus */
        TO_LIMB_T(0xfffffffe00000002), TO_LIMB_T(0xa77b4805fffcb7fd),
        TO_LIMB_T(0x6673b0101343b00a), TO_LIMB_T(0xe7db4ea6533afa90),
    };
    vec512 temp;

    ct_inverse_mod_v256(temp, inp, *get_BLST_r(), rx2);
    redc_mont_v256(out, temp, *get_BLST_r(), BLST_r0);
    mul_mont_sparse_v256(out, out, *get_BLST_rRR(), *get_BLST_r(), BLST_r0);
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void blst_fr_inverse(vec256 out, const vec256 inp)
{
    reciprocal_fr(out, inp);
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void blst_fr_eucl_inverse(vec256 out, const vec256 inp)
{
    reciprocal_fr(out, inp);
}

BLST_EXTERN_C_END

#endif /* __BLS12_381_ASM_FIELDS_H__ */
