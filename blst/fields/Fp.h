/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_FP_H
#define BLST_FP_H

#include "blst/common/vect.h"
#include "blst/fields/consts.h"
#include "blst/config.h"
#include "blst/common/defines.h"
#include "blst/backend/backend.h"
#include "blst/fields/recip-addchain.h"
#include "blst/fields/exp.h"

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

#define vec_load_global vec_copy

HOST_DEVICE_FUNC void reciprocal_fp(vec384 out, const vec384 inp);
HOST_DEVICE_FUNC void flt_reciprocal_fp(vec384 out, const vec384 inp);
HOST_DEVICE_FUNC bool_t recip_sqrt_fp(vec384 out, const vec384 inp);
HOST_DEVICE_FUNC bool_t sqrt_fp(vec384 out, const vec384 inp);

#define neg_fp(r,a) cneg_fp((r),(a),1)
//#define neg_fp2(r,a) cneg_fp2((r),(a),1)



























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

BLST_EXTERN_C_END

#endif /* BLST_FP_H */
