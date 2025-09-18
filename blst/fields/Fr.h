

#ifndef BLST_FR_H
#define BLST_FR_H

#include "blst/common/vect.h"
#include "blst/backend/backend.h"
#include "blst/fields/consts.h"


//BLST_EXTERN_C_BEGIN






/*
 * BLS12-381-specific Fr shortcuts to assembly.
 */
HOST_DEVICE_FUNC void fr_add(vec256 ret, const vec256 a, const vec256 b);
HOST_DEVICE_FUNC void fr_sub(vec256 ret, const vec256 a, const vec256 b);
HOST_DEVICE_FUNC void fr_mul_by_3(vec256 ret, const vec256 a);
HOST_DEVICE_FUNC void fr_lshift(vec256 ret, const vec256 a, size_t count);
HOST_DEVICE_FUNC void fr_rshift(vec256 ret, const vec256 a, size_t count);
HOST_DEVICE_FUNC void fr_mul(vec256 ret, const vec256 a, const vec256 b);
HOST_DEVICE_FUNC void fr_ct_bfly(vec256 x0, vec256 x1, const vec256 twiddle);
HOST_DEVICE_FUNC void fr_gs_bfly(vec256 x0, vec256 x1, const vec256 twiddle);
HOST_DEVICE_FUNC void fr_sqr(vec256 ret, const vec256 a);
HOST_DEVICE_FUNC void fr_cneg(vec256 ret, const vec256 a, int flag);
HOST_DEVICE_FUNC void fr_to(vec256 ret, const vec256 a);
HOST_DEVICE_FUNC void fr_from(vec256 ret, const vec256 a);
HOST_DEVICE_FUNC void fr_from_scalar(vec256 ret, const pow256 a);
HOST_DEVICE_FUNC void scalar_from_fr(pow256 ret, const vec256 a);
HOST_DEVICE_FUNC int scalar_fr_check(const pow256 a);
HOST_DEVICE_FUNC int sk_check(const pow256 a);
HOST_DEVICE_FUNC int sk_add_n_check(pow256 ret, const pow256 a, const pow256 b);
HOST_DEVICE_FUNC int sk_sub_n_check(pow256 ret, const pow256 a, const pow256 b);
HOST_DEVICE_FUNC int sk_mul_n_check(pow256 ret, const pow256 a, const pow256 b);
HOST_DEVICE_FUNC void sk_inverse(pow256 ret, const pow256 a);







/*
 * BLS12-381-specific Fr shortcuts to assembly.
 */
HOST_DEVICE_FUNC BLST_INLINE void fr_add(vec256 ret, const vec256 a, const vec256 b)
{
    add_mod_v256(ret, a, b, *get_BLST_r());
}

HOST_DEVICE_FUNC BLST_INLINE void fr_sub(vec256 ret, const vec256 a, const vec256 b)
{
    sub_mod_v256(ret, a, b, *get_BLST_r());
}

HOST_DEVICE_FUNC BLST_INLINE void fr_mul_by_3(vec256 ret, const vec256 a)
{
    mul_by_3_mod_v256(ret, a, *get_BLST_r());
}

HOST_DEVICE_FUNC BLST_INLINE void fr_lshift(vec256 ret, const vec256 a, size_t count)
{
    lshift_mod_v256(ret, a, count, *get_BLST_r());
}

HOST_DEVICE_FUNC BLST_INLINE void fr_rshift(vec256 ret, const vec256 a, size_t count)
{
    rshift_mod_v256(ret, a, count, *get_BLST_r());
}

HOST_DEVICE_FUNC BLST_INLINE void fr_mul(vec256 ret, const vec256 a, const vec256 b)
{
    mul_mont_sparse_v256(ret, a, b, *get_BLST_r(), BLST_r0);
}

HOST_DEVICE_FUNC BLST_INLINE void fr_ct_bfly(vec256 x0, vec256 x1, const vec256 twiddle)
{
    vec256 x2;

    mul_mont_sparse_v256(x2, x1, twiddle, *get_BLST_r(), BLST_r0);
    sub_mod_v256(x1, x0, x2, *get_BLST_r());
    add_mod_v256(x0, x0, x2, *get_BLST_r());
}

HOST_DEVICE_FUNC BLST_INLINE void fr_gs_bfly(vec256 x0, vec256 x1, const vec256 twiddle)
{
    vec256 x2;

    sub_mod_v256(x2, x0, x1, *get_BLST_r());
    add_mod_v256(x0, x0, x1, *get_BLST_r());
    mul_mont_sparse_v256(x1, x2, twiddle, *get_BLST_r(), BLST_r0);
}

HOST_DEVICE_FUNC BLST_INLINE void fr_sqr(vec256 ret, const vec256 a)
{
    sqr_mont_sparse_v256(ret, a, *get_BLST_r(), BLST_r0);
}

HOST_DEVICE_FUNC BLST_INLINE void fr_cneg(vec256 ret, const vec256 a, int flag)
{
    cneg_mod_v256(ret, a, is_zero(flag) ^ 1, *get_BLST_r());
}

HOST_DEVICE_FUNC BLST_INLINE void fr_to(vec256 ret, const vec256 a)
{
    mul_mont_sparse_v256(ret, a, *get_BLST_rRR(), *get_BLST_r(), BLST_r0);
}

HOST_DEVICE_FUNC BLST_INLINE void fr_from(vec256 ret, const vec256 a)
{
    from_mont_v256(ret, a, *get_BLST_r(), BLST_r0);
}

HOST_DEVICE_FUNC BLST_INLINE void fr_from_scalar(vec256 ret, const pow256 a)
{
    const union {
        long one;
        char little;
    } is_endian = { 1 };

    if ((uptr_t)ret == (uptr_t)a && is_endian.little) {
        mul_mont_sparse_v256(ret, (const limb_t*)a, *get_BLST_rRR(),
            *get_BLST_r(), BLST_r0);
    }
    else {
        vec256 out;
        limbs_from_le_bytes(out, a, 32);
        mul_mont_sparse_v256(ret, out, *get_BLST_rRR(), *get_BLST_r(), BLST_r0);
        vec_zero(out, sizeof(out));
    }
}

HOST_DEVICE_FUNC BLST_INLINE void scalar_from_fr(pow256 ret, const vec256 a)
{
    const union {
        long one;
        char little;
    } is_endian = { 1 };

    if ((uptr_t)ret == (uptr_t)a && is_endian.little) {
        from_mont_v256((limb_t*)ret, a, *get_BLST_r(), BLST_r0);
    }
    else {
        vec256 out;
        from_mont_v256(out, a, *get_BLST_r(), BLST_r0);
        le_bytes_from_limbs(ret, out, 32);
        vec_zero(out, sizeof(out));
    }
}

HOST_DEVICE_FUNC BLST_INLINE int scalar_fr_check(const pow256 a)
{
    return (int)(check_mod_v256(a, *get_BLST_r()) |
        bytes_are_zero(a, sizeof(pow256)));
}



HOST_DEVICE_FUNC BLST_INLINE void fr_from_uint64(vec256 ret, const unsigned long long a[4])
{
    const union {
        long one;
        char little;
    } is_endian = { 1 };

    if (sizeof(limb_t) == 4 && !is_endian.little) {
        int i;
        for (i = 0; i < 4; i++) {
            unsigned long long limb = a[i];
            ret[2 * i] = (limb_t)limb;
            ret[2 * i + 1] = (limb_t)(limb >> 32);
        }
        a = (const unsigned long long*)ret;
    }
    mul_mont_sparse_v256(ret, (const limb_t*)a, *get_BLST_rRR(), *get_BLST_r(), BLST_r0);
}

HOST_DEVICE_FUNC BLST_INLINE void fr_from_hexascii(vec256 ret, const char* hex)
{
    limbs_from_hexascii(ret, sizeof(vec256), hex);
    mul_mont_sparse_v256(ret, ret, *get_BLST_rRR(), *get_BLST_r(), BLST_r0);
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

HOST_DEVICE_FUNC BLST_FORCEINLINE void fr_inverse(vec256 out, const vec256 inp)
{
    reciprocal_fr(out, inp);
}

HOST_DEVICE_FUNC BLST_FORCEINLINE void fr_eucl_inverse(vec256 out, const vec256 inp)
{
    reciprocal_fr(out, inp);
}

HOST_DEVICE_FUNC BLST_FORCEINLINE int sk_check(const pow256 a)
{
    return (int)check_mod_v256(a, *get_BLST_r());
}

HOST_DEVICE_FUNC BLST_FORCEINLINE int sk_add_n_check(pow256 ret, const pow256 a, const pow256 b)
{
    return (int)add_n_check_mod_v256(ret, a, b, *get_BLST_r());
}

HOST_DEVICE_FUNC BLST_FORCEINLINE int sk_sub_n_check(pow256 ret, const pow256 a, const pow256 b)
{
    return (int)sub_n_check_mod_v256(ret, a, b, *get_BLST_r());
}

HOST_DEVICE_FUNC BLST_INLINE int sk_mul_n_check(pow256 ret, const pow256 a, const pow256 b)
{
    vec256 t[2];
    const union {
        long one;
        char little;
    } is_endian = { 1 };
    bool_t is_zero;

    if (((size_t)a | (size_t)b) % sizeof(limb_t) != 0 || !is_endian.little) {
        limbs_from_le_bytes(t[0], a, sizeof(pow256));
        limbs_from_le_bytes(t[1], b, sizeof(pow256));
        a = (const byte*)t[0];
        b = (const byte*)t[1];
    }
    mul_mont_sparse_v256(t[0], *get_BLST_rRR(), (const limb_t*)a, *get_BLST_r(), BLST_r0);
    mul_mont_sparse_v256(t[0], t[0], (const limb_t*)b, *get_BLST_r(), BLST_r0);
    le_bytes_from_limbs(ret, t[0], sizeof(pow256));
    is_zero = vec_is_zero(t[0], sizeof(vec256));
    vec_zero(t, sizeof(t));

    return (int)(is_zero ^ 1);
}

HOST_DEVICE_FUNC BLST_INLINE void sk_inverse(pow256 ret, const pow256 a)
{
    const union {
        long one;
        char little;
    } is_endian = { 1 };

    if (((size_t)a | (size_t)ret) % sizeof(limb_t) == 0 && is_endian.little) {
        limb_t* out = (limb_t*)ret;
        mul_mont_sparse_v256(out, (const limb_t*)a, *get_BLST_rRR(),
            *get_BLST_r(), BLST_r0);
        reciprocal_fr(out, out);
        from_mont_v256(out, out, *get_BLST_r(), BLST_r0);
    }
    else {
        vec256 out;
        limbs_from_le_bytes(out, a, 32);
        mul_mont_sparse_v256(out, out, *get_BLST_rRR(), *get_BLST_r(), BLST_r0);
        reciprocal_fr(out, out);
        from_mont_v256(out, out, *get_BLST_r(), BLST_r0);
        le_bytes_from_limbs(ret, out, 32);
        vec_zero(out, sizeof(out));
    }
}

//BLST_EXTERN_C_END


#endif