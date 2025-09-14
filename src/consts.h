/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_CONST_H__
#define __BLS12_381_ASM_CONST_H__
#include "defines.h"
#include "config.h"

extern const vec384 BLST_P_HOST;

// -1/P, this is the modular inverse of P modulo 2^64
//   used in montgomery reduction as the "magic" number
#define BLST_p0  (limb_t)0x89f3fffcfffcfffd

typedef union { vec384 p12[12]; vec384x p2; vec384 p; } radix384;

extern const radix384 BLST_Rx_HOST; /* (1<<384)%P, "radix", one-in-Montgomery */
extern const vec384 BLST_RR_HOST;   /* (1<<768)%P, "radix"^2, to-Montgomery   */


extern const vec384 BLST_RRRR_HOST; /* (1<<1536)%P, "radix"^4, hash_to_field */

#define ONE_MONT_P TO_LIMB_T(0x760900000002fffd), \
                   TO_LIMB_T(0xebf4000bc40c0002), \
                   TO_LIMB_T(0x5f48985753c758ba), \
                   TO_LIMB_T(0x77ce585370525745), \
                   TO_LIMB_T(0x5c071a97a256ec6d), \
                   TO_LIMB_T(0x15f65ec3fa80e493)


extern const vec256 BLST_r_HOST;    /* order */

// -1/r, this is the modular inverse of r modulo 2^256
// used in montgomery reduction as the "magic" number
#define BLST_r0 (limb_t)0xfffffffeffffffff  /* -1/r */

extern const vec256 BLST_rRR_HOST;  /* (1<<512)%r, "radix"^2, to-Montgomery   */

#ifdef BLST_CUDA
extern const vec384 BLST_P_DEV;
extern const radix384 BLST_Rx_DEV;
extern const vec384 BLST_RR_DEV;
extern const vec256 BLST_r_DEV;
extern const vec256 BLST_rRR_DEV;
extern const vec384 BLST_RRRR_DEV;
#endif

// returns the prime modulus for BLS12-381. 
// This is the field over which G1 is defined.
// The order of the group G1 is a different prime, which can be obtained
// via get_BLST_r().
// Will return host or device memory depending on compilation target.
BLST_FORCEINLINE const vec384* get_BLST_P(void)
{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_P_DEV;
#else
	return &BLST_P_HOST;
#endif
}

// returns R, the one-in-Montgomery form. 
// Will return host or device memory depending on compilation target.
BLST_FORCEINLINE const radix384* get_BLST_Rx(void)
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
BLST_FORCEINLINE const vec384* get_BLST_RR(void)
{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_RR_DEV;
#else
	return &BLST_RR_HOST;
#endif
}

// returns R^4, used in hash_to_field to convert 256-bit input into montgomery form.
BLST_FORCEINLINE const vec384* get_BLST_RRRR(void)
{
    #if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_RRRR_DEV;
#else
	return &BLST_RRRR_HOST;
#endif
}

// returns the order of the group G1. This prime define the scalar field.
// Will return host or device memory depending on compilation target.
BLST_FORCEINLINE const vec256* get_BLST_r(void)
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
BLST_FORCEINLINE const vec256* get_BLST_rRR(void)
{
    #if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
	return &BLST_rRR_DEV;
#else
	return &BLST_rRR_HOST;
#endif
}

#define ZERO_384 (get_BLST_Rx()->p2[1])


#endif
