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

#endif /* __BLS12_381_ASM_FIELDS_H__ */
