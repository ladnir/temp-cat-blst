/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/consts.h"
#include "blst/config.h"


namespace blst
{

    /* z = -0xd201000000010000  ,  (z-1)^2 * (z^4 - z^2 + 1)/3 + z */
    const vec384 BLST_P_HOST = { BLST_P_LIMBS };

    const radix384 BLST_Rx_HOST = { /* (1<<384)%P, "radix", one-in-Montgomery */
      { { ONE_MONT_P_LIMBS },
        { 0 } }
    };

    /* (1<<768)%P, "radix"^2, to-Montgomery */
    const vec384 BLST_RR_HOST = { BLST_RR_LIMBS };

    /* z^4 - z^2 + 1, group order */
    const vec256 BLST_r_HOST = { BLST_r_LIMBS };

    /* (1<<512)%r, "radix"^2, to-Montgomery */
    const vec256 BLST_rRR_HOST = { BLST_rRR_LIMBS };


    /* RR^2 */
    const vec384 BLST_RRRR_HOST = { BLST_RRRR_LIMBS };


    const vec384x BLST_sqrt_minus_1_HOST = { { 0 }, { ONE_MONT_P_LIMBS } };
    /*
        * "magic" number is ±2^((p-3)/4)%p, which is "1/sqrt(2)",
        * in quotes because 2*"1/sqrt(2)"^2 == -1 mod p, not 1,
        * but it pivots into "complex" plane nevertheless...
        */
    const vec384x BLST_sqrt_sqrt_minus_1_HOST = { BLST_SQRT_SQRT_MINUS_1_LIMBS };
    const vec384x BLST_sqrt_minus_sqrt_minus_1_HOST = { BLST_SQRT_MINUS_SQRT_MINUS_1_LIMBS };


#if defined(BLST_CUDA) && defined(__CUDACC__)

    //__constant__ vec384 BLST_P_DEV = { BLST_P_LIMBS };
    //__constant__ vec384 BLST_RR_DEV = { BLST_RR_LIMBS };
    //__constant__ vec256 BLST_r_DEV = { BLST_r_LIMBS };
    //__constant__ vec256 BLST_rRR_DEV = { BLST_rRR_LIMBS };
    //__constant__ vec384 BLST_RRRR_DEV = { BLST_RRRR_LIMBS };
    //__constant__ radix384 BLST_Rx_DEV = { {
    //    { ONE_MONT_P_LIMBS },
    //    { 0 }
    //} };
}

#endif // __CUDA_ARCH__


