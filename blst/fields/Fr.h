

#ifndef BLST_FR_H
#define BLST_FR_H

#include "blst/common/vect.h"

#include "blst/backend/backend_asm.h"
#include "blst/backend/backend_portable.h"

//BLST_EXTERN_C_BEGIN



namespace blst
{

    // -1/r, this is the modular inverse of r modulo 2^256
    // used in montgomery reduction as the "magic" number
#define BLST_r0 (limb_t)0xfffffffeffffffff  /* -1/r */

    struct Fr
    {
        vec_256 data;

        BLST_HD operator limb_t* () { return data; }
        BLST_HD operator const limb_t* () const { return data; }


        BLST_HD static const Fr& zero();
        BLST_HD static const Fr& one();
    };

#define BLST_r_LIMBS \
    TO_LIMB_T(0xffffffff00000001), TO_LIMB_T(0x53bda402fffe5bfe),\
    TO_LIMB_T(0x3339d80809a1d805), TO_LIMB_T(0x73eda753299d7d48)

#define BLST_rRR_LIMBS \
    TO_LIMB_T(0xc999e990f3f29c6d), TO_LIMB_T(0x2b6cedcb87925c23),\
    TO_LIMB_T(0x05d314967254398f), TO_LIMB_T(0x0748d9d99f59ff11)

    /* z^4 - z^2 + 1, group order */
    constexpr Fr BLST_r_HOST = { BLST_r_LIMBS };

    /* (1<<512)%r, "radix"^2, to-Montgomery */
    constexpr Fr BLST_rRR_HOST = { BLST_rRR_LIMBS };


#if defined(BLST_CUDA) && defined(__CUDACC__)

    static __constant__ Fr BLST_r_DEV = { BLST_r_LIMBS };
    static __constant__ Fr BLST_rRR_DEV = { BLST_rRR_LIMBS };


#endif // __CUDA_ARCH__

    // returns the order of the group G1. This prime define the scalar field.
    // Will return host or device memory depending on compilation target.
    BLST_HD BLST_FORCEINLINE const Fr& get_BLST_r(void)
    {
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
        return BLST_r_DEV;
#else
        return BLST_r_HOST;
#endif
    }

    // returns R^2 mod r, the to-Montgomery form value for the group order r.
    // Perform montgomery multiplication with this value to convert a normal
    // integer into montgomery form.
    // Will return host or device memory depending on compilation target.
    BLST_HD BLST_FORCEINLINE const Fr& get_BLST_rRR(void)
    {
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
        return BLST_rRR_DEV;
#else
        return BLST_rRR_HOST;
#endif
    }



    ///*
    // * BLS12-381-specific Fr shortcuts to assembly.
    // */
    //BLST_HD void fr_add(Fr& ret, const Fr& a, const Fr& b);
    //BLST_HD void fr_sub(Fr& ret, const Fr& a, const Fr& b);
    //BLST_HD void fr_mul_by_3(Fr& ret, const Fr& a);
    //BLST_HD void fr_lshift(Fr& ret, const Fr& a, size_t count);
    //BLST_HD void fr_rshift(Fr& ret, const Fr& a, size_t count);
    //BLST_HD void fr_mul(Fr& ret, const Fr& a, const Fr& b);
    //BLST_HD void fr_ct_bfly(Fr& x0, Fr& x1, const Fr& twiddle);
    //BLST_HD void fr_gs_bfly(Fr& x0, Fr& x1, const Fr& twiddle);
    //BLST_HD void fr_sqr(Fr& ret, const Fr& a);
    //BLST_HD void fr_cneg(Fr& ret, const Fr& a, int flag);
    //BLST_HD void fr_to(Fr& ret, const Fr& a);
    //BLST_HD void fr_from(Fr& ret, const Fr& a);
    //BLST_HD void fr_from_scalar(Fr& ret, const pow256 a);
    //BLST_HD void scalar_from_fr(pow256 ret, const Fr& a);
    //BLST_HD int scalar_fr_check(const pow256 a);
    //BLST_HD int sk_check(const pow256 a);
    //BLST_HD int sk_add_n_check(pow256 ret, const pow256 a, const pow256 b);
    //BLST_HD int sk_sub_n_check(pow256 ret, const pow256 a, const pow256 b);
    //BLST_HD int sk_mul_n_check(pow256 ret, const pow256 a, const pow256 b);
    //BLST_HD void sk_inverse(pow256 ret, const pow256 a);



    BLST_HD BLST_FORCEINLINE limb_t check_mod_v256(const pow256 a, const Fr& p)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        return check_mod_256(a, &p);
#else
        return cc_check_mod_n<NLIMBS(256)>(a, p.data);
#endif
    }

    BLST_HD BLST_FORCEINLINE limb_t add_n_check_mod_v256(pow256 ret, const pow256 a, const pow256 b,
        const Fr& p)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        return add_n_check_mod_256(ret, a, b, &p);
#else
        //return cc_add_n_check_mod_256(ret, a, b, p);
        return cc_add_n_check_mod_n<NLIMBS(256)>(ret, a, b, p.data);
#endif
    }

    // performs modular subtraction with zero-checking - 
    // it computes (a - b) mod p and returns whether the result is non-zero.
    BLST_HD BLST_FORCEINLINE limb_t sub_n_check_mod_v256(pow256 ret, const pow256 a, const pow256 b,
        const Fr& p)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        return sub_n_check_mod_256(ret, a, b, &p);
#else
        //return cc_sub_n_check_mod_256(ret, a, b, p);
        return cc_sub_n_check_mod_n<NLIMBS(256)>(ret, a, b, p.data);
#endif
    }


    //BLST_HD BLST_INLINE void cc_mul_mont_sparse_256(
    //    Fr& ret, const Fr& a, const Fr& b, const Fr& p, limb n0)
    //{
    //}
    //BLST_HD BLST_INLINE void cc_sqr_mont_sparse_256(
    //    Fr& ret, const Fr& a, const Fr& p, limb n0)
    //{
    //}


    //BLST_HD BLST_INLINE void cc_add_mod_256(Fr& ret, const Fr& a,
    //    const Fr& b, const Fr& p)
    //{
    //}
    //BLST_HD BLST_INLINE void cc_sub_mod_256(Fr& ret, const Fr& a,
    //    const Fr& b, const Fr& p)
    //{
    //}
    //BLST_HD BLST_INLINE void cc_mul_by_3_mod_256(Fr& ret, const Fr& a,
    //    const Fr& p)
    //{
    //}

    //BLST_HD BLST_INLINE void cc_lshift_mod_256(Fr& ret, const Fr& a, size_t count,
    //    const Fr& p)
    //{
    //}

    //BLST_HD BLST_INLINE void cc_cneg_mod_256(Fr& ret, const Fr& a, bool_t flag,
    //    const Fr& p)
    //{
    //}	
    
    //BLST_HD BLST_INLINE limb cc_check_mod_256(const pow256 a, const Fr& p)
    //{
    //}

    //BLST_HD BLST_INLINE limb cc_add_n_check_mod_256(pow256 ret, const pow256 a,
    //    const pow256 b, const Fr& p)
    //{
    //}
    //BLST_HD BLST_INLINE limb cc_sub_n_check_mod_256(pow256 ret, const pow256 a,
    //    const pow256 b, const Fr& p)
    //{
    //}
    //BLST_HD BLST_INLINE void cc_from_mont_256(Fr& ret, const Fr& a,
    //    const Fr& p, limb n0)
    //{
    //}

    //BLST_HD BLST_INLINE void cc_redc_mont_256(Fr& ret, const UFr& a,
    //    const Fr& p, limb n0)
    //{
    //}


    //BLST_HD BLST_INLINE void cc_rshift_mod_256(Fr& ret, const Fr& a, size_t count,
    //    const Fr& p)
    //{
    //}


    BLST_HD BLST_INLINE void cc_ct_inverse_mod_256(UFr& ret, const Fr& inp,
        const Fr& mod, const Fr& modx)
    {
        cc_ct_inverse_mod_n<NLIMBS(256)>(ret, inp, mod, modx);
    }











    // montgomery multiply: ret = a * b * R^-1 mod p, where R montgomery radix
    // optimized for "spare" modulus p, i.e. with many zero limbs
    BLST_HD BLST_FORCEINLINE void mul_mont_sparse_v256(Fr& ret, const Fr& a, const Fr& b,
        const Fr& p, limb_t n0)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        mulx_mont_sparse_256(&ret, &a, &b, &p, n0);// linked by asm
#else
        //cc_mul_mont_sparse_256(ret, a, b, p, n0);
        cc_mul_mont_n<NLIMBS(256)>(ret.data, a.data, b.data, p.data, n0);
#endif
    }

    // montgomery squaring: ret = a^2 * R^-1 mod p, where R montgomery radix
    // optimized for "spare" modulus p, i.e. with many zero limbs
    BLST_HD BLST_FORCEINLINE void sqr_mont_sparse_v256(Fr& ret, const Fr& a, const Fr& p, limb_t n0)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        sqrx_mont_sparse_256(&ret, &a, &p, n0);// linked by asm
#else
        //cc_sqr_mont_sparse_256(ret, a, p, n0);
        cc_mul_mont_n<NLIMBS(256)>(ret.data, a.data, a.data, p.data, n0);

#endif
    }

    // montgomery multiply: ret = a * b * R^-1 mod p, where R montgomery radix
    BLST_HD BLST_FORCEINLINE void redc_mont_v256(Fr& ret, const UFr& a, const Fr& p, limb_t n0)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        redcx_mont_256(&ret, &a, &p, n0);// linked by asm
#else
        //cc_redc_mont_256(ret, a, p, n0);
        cc_redc_mont_n<NLIMBS(256)>(ret.data, a.data, p.data, n0);
#endif
    }

    // from montgomery domain: ret = a * R^-1 mod p, where R montgomery radix
    BLST_HD BLST_FORCEINLINE void from_mont_v256(Fr& ret, const Fr& a, const Fr& p, limb_t n0)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        fromx_mont_256(&ret, &a, &p, n0);// linked by asm
#else
        //cc_from_mont_256(ret, a, p, n0);
        cc_from_mont_n<NLIMBS(256)>(ret.data, a.data, p.data, n0);
#endif
    }


    BLST_HD BLST_FORCEINLINE void add_mod_v256(Fr& ret, const Fr& a, const Fr& b, const Fr& p)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        add_mod_256(&ret, &a, &b, &p);// linked by asm, 
#else
        //cc_add_mod_256(ret, a, b, p);
        cc_add_mod_n<NLIMBS(256)>(ret.data, a.data, b.data, p.data);
#endif

    }

    BLST_HD BLST_FORCEINLINE void sub_mod_v256(Fr& ret, const Fr& a, const Fr& b, const Fr& p)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        sub_mod_256(&ret, &a, &b, &p);// linked by asm
#else
        //cc_sub_mod_256(ret, a, b, p);
        cc_sub_mod_n<NLIMBS(256)>(ret.data, a.data, b.data, p.data);
#endif
    }

    BLST_HD BLST_FORCEINLINE void mul_by_3_mod_v256(Fr& ret, const Fr& a, const Fr& p)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        mul_by_3_mod_256(&ret, &a, &p);
#else
        //cc_mul_by_3_mod_256(ret, a, p);
        cc_mul_by_3_mod_n<NLIMBS(256)>(ret.data, a.data, p.data);
#endif
    }

    BLST_HD BLST_FORCEINLINE void cneg_mod_v256(Fr& ret, const Fr& a, bool_t flag, const Fr& p)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        cneg_mod_256(&ret, &a, flag, &p);
#else
        //cc_cneg_mod_256(ret, a, flag, p);
        cc_cneg_mod_n<NLIMBS(256)>(ret.data, a.data, flag, p.data);
#endif // BLST_DEFAULT_ASM

    }
    BLST_HD BLST_FORCEINLINE void lshift_mod_v256(Fr& ret, const Fr& a, size_t count, const Fr& p)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        lshift_mod_256(&ret, &a, count, &p);
#else
        //cc_lshift_mod_256(ret, a, count, p);
        cc_lshift_mod_n<NLIMBS(256)>(ret.data, a.data, count, p.data);
#endif
    }

    BLST_HD BLST_FORCEINLINE void rshift_mod_v256(Fr& ret, const Fr& a, size_t count, const Fr& p)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        rshift_mod_256(&ret, &a, count, &p);
#else
        //cc_rshift_mod_256(ret, a, count, p);
        cc_rshift_mod_n<NLIMBS(256)>(ret.data, a.data, count, p.data);
#endif
    }

    BLST_HD BLST_FORCEINLINE void ct_inverse_mod_v256(UFr& ret, const Fr& inp, const Fr& mod,
        const Fr& modx)
    {
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
        ct_inverse_mod_256(&ret, &inp, &mod, &modx);
#else	
        cc_ct_inverse_mod_256(ret, inp, mod, modx);
#endif
    }
    //
    //bool_t eucl_inverse_mod_v256(Fr& ret, const Fr& a, const Fr& p, const Fr& one)
    //{
    //#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
    //    bool_t eucl_inverse_mod_256(Fr& ret, const Fr& a, const Fr& p, const Fr& one)
    //#else
    //	return cc_eucl_inverse_mod_256(ret, a, p, one);
    //#endif
    //
    //}









    /*
     * BLS12-381-specific Fr shortcuts to assembly.
     */
    BLST_HD BLST_INLINE void fr_add(Fr& ret, const Fr& a, const Fr& b)
    {
        add_mod_v256(ret, a, b, get_BLST_r());
    }

    BLST_HD BLST_INLINE void fr_sub(Fr& ret, const Fr& a, const Fr& b)
    {
        sub_mod_v256(ret, a, b, get_BLST_r());
    }

    BLST_HD BLST_INLINE void fr_mul_by_3(Fr& ret, const Fr& a)
    {
        mul_by_3_mod_v256(ret, a, get_BLST_r());
    }

    BLST_HD BLST_INLINE void fr_lshift(Fr& ret, const Fr& a, size_t count)
    {
        lshift_mod_v256(ret, a, count, get_BLST_r());
    }

    BLST_HD BLST_INLINE void fr_rshift(Fr& ret, const Fr& a, size_t count)
    {
        rshift_mod_v256(ret, a, count, get_BLST_r());
    }

    BLST_HD BLST_INLINE void fr_mul(Fr& ret, const Fr& a, const Fr& b)
    {
        mul_mont_sparse_v256(ret, a, b, get_BLST_r(), BLST_r0);
    }

    BLST_HD BLST_INLINE void fr_ct_bfly(Fr& x0, Fr& x1, const Fr& twiddle)
    {
        Fr x2;

        mul_mont_sparse_v256(x2, x1, twiddle, get_BLST_r(), BLST_r0);
        sub_mod_v256(x1, x0, x2, get_BLST_r());
        add_mod_v256(x0, x0, x2, get_BLST_r());
    }

    BLST_HD BLST_INLINE void fr_gs_bfly(Fr& x0, Fr& x1, const Fr& twiddle)
    {
        Fr x2;

        sub_mod_v256(x2, x0, x1, get_BLST_r());
        add_mod_v256(x0, x0, x1, get_BLST_r());
        mul_mont_sparse_v256(x1, x2, twiddle, get_BLST_r(), BLST_r0);
    }

    BLST_HD BLST_INLINE void fr_sqr(Fr& ret, const Fr& a)
    {
        sqr_mont_sparse_v256(ret, a, get_BLST_r(), BLST_r0);
    }

    BLST_HD BLST_INLINE void fr_cneg(Fr& ret, const Fr& a, int flag)
    {
        cneg_mod_v256(ret, a, is_zero(flag) ^ 1, get_BLST_r());
    }

    BLST_HD BLST_INLINE void fr_to(Fr& ret, const Fr& a)
    {
        mul_mont_sparse_v256(ret, a, get_BLST_rRR(), get_BLST_r(), BLST_r0);
    }

    BLST_HD BLST_INLINE void fr_from(Fr& ret, const Fr& a)
    {
        from_mont_v256(ret, a, get_BLST_r(), BLST_r0);
    }

    BLST_HD BLST_INLINE void fr_from_scalar(Fr& ret, const pow256 a)
    {
        const union {
            long one;
            char little;
        } is_endian = { 1 };

        if ((void*)ret == (void*)a && is_endian.little) {
            mul_mont_sparse_v256(ret, (const Fr&)a, get_BLST_rRR(),
                get_BLST_r(), BLST_r0);
        }
        else {
            Fr out;
            limbs_from_le_bytes(out, a, 32);
            mul_mont_sparse_v256(ret, out, get_BLST_rRR(), get_BLST_r(), BLST_r0);
            vec_zero(out, sizeof(out));
        }
    }

    BLST_HD BLST_INLINE void scalar_from_fr(pow256 ret, const Fr& a)
    {
        const union {
            long one;
            char little;
        } is_endian = { 1 };

        if ((void*)ret == (void*)a.data && is_endian.little) {
            from_mont_v256((Fr&)ret, a, get_BLST_r(), BLST_r0);
        }
        else {
            Fr out;
            from_mont_v256(out, a, get_BLST_r(), BLST_r0);
            le_bytes_from_limbs(ret, out, 32);
            vec_zero(&out, sizeof(out));
        }
    }

    BLST_HD BLST_INLINE int scalar_fr_check(const pow256 a)
    {
        return (int)(check_mod_v256(a, get_BLST_r()) |
            bytes_are_zero(a, sizeof(pow256)));
    }



    BLST_HD BLST_INLINE void fr_from_uint64(Fr& ret, const unsigned long long a[4])
    {
        const union {
            long one;
            char little;
        } is_endian = { 1 };

        if (sizeof(limb_t) == 4 && !is_endian.little) {
            int i;
            for (i = 0; i < 4; i++) {
                unsigned long long L = a[i];
                ret[2 * i] = (limb_t)L;
                ret[2 * i + 1] = (limb_t)(L >> 32);
            }
            a = (const unsigned long long*)ret;
        }
        mul_mont_sparse_v256(ret, (const Fr&)a, get_BLST_rRR(), get_BLST_r(), BLST_r0);
    }

    BLST_HD BLST_INLINE void fr_from_hexascii(Fr& ret, const char* hex)
    {
        limbs_from_hexascii(ret, sizeof(Fr), hex);
        mul_mont_sparse_v256(ret, ret, get_BLST_rRR(), get_BLST_r(), BLST_r0);
    }





    BLST_HD BLST_INLINE void reciprocal_fr(Fr& out, const Fr& inp)
    {
        const Fr rx2 = { /* left-aligned value of the modulus */
            TO_LIMB_T(0xfffffffe00000002), TO_LIMB_T(0xa77b4805fffcb7fd),
            TO_LIMB_T(0x6673b0101343b00a), TO_LIMB_T(0xe7db4ea6533afa90),
        };
        UFr temp;

        ct_inverse_mod_v256(temp, inp, get_BLST_r(), rx2);
        redc_mont_v256(out, temp, get_BLST_r(), BLST_r0);
        mul_mont_sparse_v256(out, out, get_BLST_rRR(), get_BLST_r(), BLST_r0);
    }

    BLST_HD BLST_FORCEINLINE void fr_inverse(Fr& out, const Fr& inp)
    {
        reciprocal_fr(out, inp);
    }

    BLST_HD BLST_FORCEINLINE void fr_eucl_inverse(Fr& out, const Fr& inp)
    {
        reciprocal_fr(out, inp);
    }

    BLST_HD BLST_FORCEINLINE int sk_check(const pow256 a)
    {
        return (int)check_mod_v256(a, get_BLST_r());
    }

    BLST_HD BLST_FORCEINLINE int sk_add_n_check(pow256 ret, const pow256 a, const pow256 b)
    {
        return (int)add_n_check_mod_v256(ret, a, b, get_BLST_r());
    }

    BLST_HD BLST_FORCEINLINE int sk_sub_n_check(pow256 ret, const pow256 a, const pow256 b)
    {
        return (int)sub_n_check_mod_v256(ret, a, b, get_BLST_r());
    }

    BLST_HD BLST_INLINE int sk_mul_n_check(pow256 ret, const pow256 a, const pow256 b)
    {
        Fr t[2];
        const union {
            long one;
            char little;
        } is_endian = { 1 };
        bool_t is_zero;

        if (((size_t)a | (size_t)b) % sizeof(limb_t) != 0 || !is_endian.little) {
            limbs_from_le_bytes(t[0], a, sizeof(pow256));
            limbs_from_le_bytes(t[1], b, sizeof(pow256));
            a = (const byte*)&t[0];
            b = (const byte*)&t[1];
        }
        mul_mont_sparse_v256(t[0], get_BLST_rRR(), (const Fr&)a, get_BLST_r(), BLST_r0);
        mul_mont_sparse_v256(t[0], t[0], (const Fr&)b, get_BLST_r(), BLST_r0);
        le_bytes_from_limbs(ret, t[0], sizeof(pow256));
        is_zero = vec_is_zero(t[0], sizeof(Fr&));
        vec_zero(t, sizeof(t));

        return (int)(is_zero ^ 1);
    }

    BLST_HD BLST_INLINE void sk_inverse(pow256 ret, const pow256 a)
    {
        const union {
            long one;
            char little;
        } is_endian = { 1 };

        if (((size_t)a | (size_t)ret) % sizeof(limb_t) == 0 && is_endian.little) {
            auto& out = (Fr&)ret;
            mul_mont_sparse_v256(out, (const Fr&)a, get_BLST_rRR(),
                get_BLST_r(), BLST_r0);
            reciprocal_fr(out, out);
            from_mont_v256(out, out, get_BLST_r(), BLST_r0);
        }
        else {
            Fr out;
            limbs_from_le_bytes(out, a, 32);
            mul_mont_sparse_v256(out, out, get_BLST_rRR(), get_BLST_r(), BLST_r0);
            reciprocal_fr(out, out);
            from_mont_v256(out, out, get_BLST_r(), BLST_r0);
            le_bytes_from_limbs(ret, out, 32);
            vec_zero(out, sizeof(out));
        }
    }


    BLST_HD BLST_FORCEINLINE const Fr& Fr::one() {
        static Fr one_in_montgomery;
        static bool initialized = false;

        if (!initialized) {
            // Compute 1 in Montgomery form: (1 * R) mod r
            Fr one_normal;
            one_normal.data[0] = 1; // Set the first limb_t to 1 (normal integer 1)
            for (size_t i = 1; i < NLIMBS(256); ++i) {
                one_normal.data[i] = 0; // Zero out the rest of the limbs
            }

            // Convert to Montgomery form: (1 * R) mod r
            mul_mont_sparse_v256(one_in_montgomery, one_normal, get_BLST_rRR(), get_BLST_r(), BLST_r0);

            initialized = true;
        }

        return one_in_montgomery;
    }


    //BLST_EXTERN_C_END
}

#endif