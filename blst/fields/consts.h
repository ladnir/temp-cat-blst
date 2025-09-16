/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_CONST_H__
#define __BLS12_381_ASM_CONST_H__
#include "blst/common/defines.h"
#include "blst/config.h"
#include "blst/common/vect.h"

BLST_EXTERN_C_BEGIN

// -1/P, this is the modular inverse of P modulo 2^64
//   used in montgomery reduction as the "magic" number
#define BLST_p0  (limb_t)0x89f3fffcfffcfffd
// -1/r, this is the modular inverse of r modulo 2^256
// used in montgomery reduction as the "magic" number
#define BLST_r0 (limb_t)0xfffffffeffffffff  /* -1/r */

#define ONE_MONT_P_LIMBS TO_LIMB_T(0x760900000002fffd), \
                   TO_LIMB_T(0xebf4000bc40c0002), \
                   TO_LIMB_T(0x5f48985753c758ba), \
                   TO_LIMB_T(0x77ce585370525745), \
                   TO_LIMB_T(0x5c071a97a256ec6d), \
                   TO_LIMB_T(0x15f65ec3fa80e493)


 /* (z-1)^2 * (z^4 - z^2 + 1)/3 + z */
#define BLST_P_LIMBS  \
    TO_LIMB_T(0xb9feffffffffaaab), TO_LIMB_T(0x1eabfffeb153ffff),\
    TO_LIMB_T(0x6730d2a0f6b0f624), TO_LIMB_T(0x64774b84f38512bf),\
    TO_LIMB_T(0x4b1ba7b6434bacd7), TO_LIMB_T(0x1a0111ea397fe69a)

#define BLST_RR_LIMBS \
    TO_LIMB_T(0xf4df1f341c341746), TO_LIMB_T(0x0a76e6a609d104f1),\
    TO_LIMB_T(0x8de5476c4c95b6d5), TO_LIMB_T(0x67eb88a9939d83c0),\
    TO_LIMB_T(0x9a793e85b519952d), TO_LIMB_T(0x11988fe592cae3aa)

#define BLST_r_LIMBS \
    TO_LIMB_T(0xffffffff00000001), TO_LIMB_T(0x53bda402fffe5bfe),\
    TO_LIMB_T(0x3339d80809a1d805), TO_LIMB_T(0x73eda753299d7d48)

#define BLST_rRR_LIMBS \
    TO_LIMB_T(0xc999e990f3f29c6d), TO_LIMB_T(0x2b6cedcb87925c23),\
    TO_LIMB_T(0x05d314967254398f), TO_LIMB_T(0x0748d9d99f59ff11)

#define BLST_RRRR_LIMBS \
    TO_LIMB_T(0xed48ac6bd94ca1e0), TO_LIMB_T(0x315f831e03a7adf8),\
    TO_LIMB_T(0x9a53352a615e29dd), TO_LIMB_T(0x34c04e5e921e1761),\
    TO_LIMB_T(0x2512d43565724728), TO_LIMB_T(0x0aa6346091755d4d)



typedef union { vec384 p12[12]; vec384x p2; vec384 p; } radix384;

extern const vec384 BLST_P_HOST;
extern const vec384 BLST_RR_HOST;   /* (1<<768)%P, "radix"^2, to-Montgomery   */
extern const vec384 BLST_RRRR_HOST; /* (1<<1536)%P, "radix"^4, hash_to_field */
extern const vec256 BLST_r_HOST;    /* order */
extern const vec256 BLST_rRR_HOST;  /* (1<<512)%r, "radix"^2, to-Montgomery   */
extern const radix384 BLST_Rx_HOST; /* (1<<384)%P, "radix", one-in-Montgomery */


#if defined(BLST_CUDA) && defined(__CUDACC__)

//static __constant__ vec384 BLST_P_DEV;
//static __constant__ vec384 BLST_RR_DEV;
//static __constant__ vec256 BLST_r_DEV;
//static __constant__ vec256 BLST_rRR_DEV;
//static __constant__ vec384 BLST_RRRR_DEV;
//static __constant__ radix384 BLST_Rx_DEV;
static __constant__ vec384 BLST_P_DEV = { BLST_P_LIMBS };
static __constant__ vec384 BLST_RR_DEV = { BLST_RR_LIMBS };
static __constant__ vec256 BLST_r_DEV = { BLST_r_LIMBS };
static __constant__ vec256 BLST_rRR_DEV = { BLST_rRR_LIMBS };
static __constant__ vec384 BLST_RRRR_DEV = { BLST_RRRR_LIMBS };
static __constant__ radix384 BLST_Rx_DEV = { {
    { ONE_MONT_P_LIMBS },
    { 0 }
} };


#endif // __CUDA_ARCH__

// returns the prime modulus for BLS12-381. 
// This is the field over which G1 is defined.
// The order of the group G1 is a different prime, which can be obtained
// via get_BLST_r().
// Will return host or device memory depending on compilation target.
HOST_DEVICE_FUNC BLST_FORCEINLINE const vec384* get_BLST_P(void)
{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_P_DEV;
#else
	return &BLST_P_HOST;
#endif
}

// returns R, the one-in-Montgomery form. 
// Will return host or device memory depending on compilation target.
HOST_DEVICE_FUNC BLST_FORCEINLINE const radix384* get_BLST_Rx(void)
{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
    return &BLST_Rx_DEV;
#else
    return &BLST_Rx_HOST;
#endif
}

// returns R^2, the to-Montgomery form value. Perform montgomery multiplication
// with this value to convert a normal integer into montgomery form.
// Will return host or device memory depending on compilation target.
HOST_DEVICE_FUNC BLST_FORCEINLINE const vec384* get_BLST_RR(void)
{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_RR_DEV;
#else
	return &BLST_RR_HOST;
#endif
}

// returns R^4, used in hash_to_field to convert 256-bit input into montgomery form.
HOST_DEVICE_FUNC BLST_FORCEINLINE const vec384* get_BLST_RRRR(void)
{
    #if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_RRRR_DEV;
#else
	return &BLST_RRRR_HOST;
#endif
}

// returns the order of the group G1. This prime define the scalar field.
// Will return host or device memory depending on compilation target.
HOST_DEVICE_FUNC BLST_FORCEINLINE const vec256* get_BLST_r(void)
{
    #if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_r_DEV;
#else
	return &BLST_r_HOST;
#endif
}

// returns R^2 mod r, the to-Montgomery form value for the group order r.
// Perform montgomery multiplication with this value to convert a normal
// integer into montgomery form.
// Will return host or device memory depending on compilation target.
HOST_DEVICE_FUNC BLST_FORCEINLINE const vec256* get_BLST_rRR(void)
{
    #if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_rRR_DEV;
#else
	return &BLST_rRR_HOST;
#endif
}

#define ZERO_384 (get_BLST_Rx()->p2[1])





BLST_EXTERN_C_END
#endif
