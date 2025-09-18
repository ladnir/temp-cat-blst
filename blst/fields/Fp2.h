/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_FIELDS_H__
#define __BLS12_381_ASM_FIELDS_H__

#include "blst/common/vect.h"
#include "blst/fields/consts.h"
#include "blst/config.h"
#include "blst/common/defines.h"
#include "blst/backend/backend.h"
#include "blst/fields/recip-addchain.h"
#include "blst/fields/exp.h"
#include "blst/fields/Fp.h"

BLST_EXTERN_C_BEGIN

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


HOST_DEVICE_FUNC void reciprocal_fp2(vec384x out, const vec384x inp);
HOST_DEVICE_FUNC void flt_reciprocal_fp2(vec384x out, const vec384x inp);
HOST_DEVICE_FUNC bool_t recip_sqrt_fp2(vec384x out, const vec384x inp,
                             const vec384x recip_ZZZ, const vec384x magic_ZZZ);
HOST_DEVICE_FUNC bool_t sqrt_fp2(vec384x out, const vec384x inp);
HOST_DEVICE_FUNC bool_t sqrt_align_fp2(vec384x out, const vec384x ret,
                             const vec384x sqrt, const vec384x inp);


HOST_DEVICE_FUNC BLST_FORCEINLINE void conjugate_fp2(vec384x ret, const vec384x a)
{
    vec_copy(ret[0], a[0], sizeof(vec384));  /* Re(ret) = Re(a) */
    cneg_fp(ret[1], a[1], 1);               /* Im(ret) = -Im(a) */
}

//#define neg_fp(r,a) cneg_fp((r),(a),1)
#define neg_fp2(r,a) cneg_fp2((r),(a),1)


































HOST_DEVICE_FUNC BLST_INLINE bool_t sqrt_align_fp2(vec384x out, const vec384x ret,
    const vec384x sqrt, const vec384x inp)
{
    vec384x coeff;
    bool_t is_sqrt;

    {
        vec384x t0;
        bool_t flag;
        /*
         * Instead of multiple trial squarings we can perform just one
         * and see if the result is "rotated by multiple of 90°" in
         * relation to |inp|, and "rotate" |ret| accordingly.
         */
        sqr_fp2(t0, sqrt);
        /* "sqrt(|inp|)"^2 = (a + b*i)^2 = (a^2-b^2) + 2ab*i */

        /* (a^2-b^2) + 2ab*i == |inp| ? |ret| is spot on */
        vec384x t1;
        sub_fp2(t1, t0, inp);
        is_sqrt = vec_is_zero(t1, sizeof(t1));
        vec_copy(coeff, get_BLST_Rx()->p2, sizeof(coeff));

        /* -(a^2-b^2) - 2ab*i == |inp| ? "rotate |ret| by 90°" */
        add_fp2(t1, t0, inp);
        vec_select(coeff, *get_BLST_sqrt_minus_1(), coeff, sizeof(coeff),
            flag = vec_is_zero(t1, sizeof(t1)));
        is_sqrt |= flag;

        /* 2ab - (a^2-b^2)*i == |inp| ? "rotate |ret| by 135°" */
        sub_fp(t1[0], t0[0], inp[1]);
        add_fp(t1[1], t0[1], inp[0]);
        vec_select(coeff, *get_BLST_sqrt_sqrt_minus_1(), coeff, sizeof(coeff),
            flag = vec_is_zero(t1, sizeof(t1)));
        is_sqrt |= flag;

        /* -2ab + (a^2-b^2)*i == |inp| ? "rotate |ret| by 45°" */
        add_fp(t1[0], t0[0], inp[1]);
        sub_fp(t1[1], t0[1], inp[0]);
        vec_select(coeff, *get_BLST_sqrt_minus_sqrt_minus_1(), coeff, sizeof(coeff),
            flag = vec_is_zero(t1, sizeof(t1)));
        is_sqrt |= flag;
    }

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

    {
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
    }

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


/*
 * Fp2  = Fp[u]  / (u^2 + 1)
 * Fp6  = Fp2[v] / (v^3 - u - 1)
 * Fp12 = Fp6[w] / (w^2 - v)
 */

HOST_DEVICE_FUNC BLST_INLINE void mul_by_u_plus_1_fp2(vec384x ret, const vec384x a)
{
    mul_by_1_plus_i_mod_v384x(ret, a, *get_BLST_P());
}

/*
 * caveat lector! |n| has to be non-zero and not more than 3!
 */
HOST_DEVICE_FUNC BLST_INLINE void frobenius_map_fp2(vec384x ret, const vec384x a, size_t n)
{
    vec_copy(ret[0], a[0], sizeof(ret[0]));
    cneg_fp(ret[1], a[1], n & 1);
}


BLST_EXTERN_C_END

#endif /* __BLS12_381_ASM_FIELDS_H__ */
