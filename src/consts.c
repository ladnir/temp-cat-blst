/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "consts.h"
#include "config.h"

/* z = -0xd201000000010000 */
const vec384 BLST_P_HOST = {    /* (z-1)^2 * (z^4 - z^2 + 1)/3 + z */
    TO_LIMB_T(0xb9feffffffffaaab), TO_LIMB_T(0x1eabfffeb153ffff),
    TO_LIMB_T(0x6730d2a0f6b0f624), TO_LIMB_T(0x64774b84f38512bf),
    TO_LIMB_T(0x4b1ba7b6434bacd7), TO_LIMB_T(0x1a0111ea397fe69a)
};

const radix384 BLST_Rx_HOST = { /* (1<<384)%P, "radix", one-in-Montgomery */
  { { ONE_MONT_P },
    { 0 } }
};

const vec384 BLST_RR_HOST = {   /* (1<<768)%P, "radix"^2, to-Montgomery */
    TO_LIMB_T(0xf4df1f341c341746), TO_LIMB_T(0x0a76e6a609d104f1),
    TO_LIMB_T(0x8de5476c4c95b6d5), TO_LIMB_T(0x67eb88a9939d83c0),
    TO_LIMB_T(0x9a793e85b519952d), TO_LIMB_T(0x11988fe592cae3aa)
};

const vec256 BLST_r_HOST = {    /* z^4 - z^2 + 1, group order */
    TO_LIMB_T(0xffffffff00000001), TO_LIMB_T(0x53bda402fffe5bfe),
    TO_LIMB_T(0x3339d80809a1d805), TO_LIMB_T(0x73eda753299d7d48)
};

const vec256 BLST_rRR_HOST = {  /* (1<<512)%r, "radix"^2, to-Montgomery */
    TO_LIMB_T(0xc999e990f3f29c6d), TO_LIMB_T(0x2b6cedcb87925c23),
    TO_LIMB_T(0x05d314967254398f), TO_LIMB_T(0x0748d9d99f59ff11)
};


const vec384 BLST_RRRR_HOST = {  /* RR^2 */
    TO_LIMB_T(0xed48ac6bd94ca1e0), TO_LIMB_T(0x315f831e03a7adf8),
    TO_LIMB_T(0x9a53352a615e29dd), TO_LIMB_T(0x34c04e5e921e1761),
    TO_LIMB_T(0x2512d43565724728), TO_LIMB_T(0x0aa6346091755d4d)
};

#ifdef __CUDA_ARCH__

__device__ const vec384 BLST_P_DEV = {    /* (z-1)^2 * (z^4 - z^2 + 1)/3 + z */
    TO_LIMB_T(0xb9feffffffffaaab), TO_LIMB_T(0x1eabfffeb153ffff),
    TO_LIMB_T(0x6730d2a0f6b0f624), TO_LIMB_T(0x64774b84f38512bf),
    TO_LIMB_T(0x4b1ba7b6434bacd7), TO_LIMB_T(0x1a0111ea397fe69a)
};
__device__ const radix384 BLST_Rx_DEV = { /* (1<<384)%P, "radix", one-in-Montgomery */
  { { ONE_MONT_P },
    { 0 } }
};
__device__ const vec384 BLST_RR_DEV = {   /* (1<<768)%P, "radix"^2, to-Montgomery */
    TO_LIMB_T(0xf4df1f341c341746), TO_LIMB_T(0x0a76e6a609d104f1),
    TO_LIMB_T(0x8de5476c4c95b6d5), TO_LIMB_T(0x67eb88a9939d83c0),
    TO_LIMB_T(0x9a793e85b519952d), TO_LIMB_T(0x11988fe592cae3aa)
};
__device__ const vec256 BLST_r_DEV = {    /* z^4 - z^2 + 1, group order */
    TO_LIMB_T(0xffffffff00000001), TO_LIMB_T(0x53bda402fffe5bfe),
    TO_LIMB_T(0x3339d80809a1d805), TO_LIMB_T(0x73eda753299d7d48)
};
__device__ const vec256 BLST_rRR_DEV = {  /* (1<<512)%r, "radix"^2, to-Montgomery */
    TO_LIMB_T(0xc999e990f3f29c6d), TO_LIMB_T(0x2b6cedcb87925c23),
    TO_LIMB_T(0x05d314967254398f), TO_LIMB_T(0x0748d9d99f59ff11)
};

const vec384 BLST_RRRR_DEV = {  /* RR^2 */
    TO_LIMB_T(0xed48ac6bd94ca1e0), TO_LIMB_T(0x315f831e03a7adf8),
    TO_LIMB_T(0x9a53352a615e29dd), TO_LIMB_T(0x34c04e5e921e1761),
    TO_LIMB_T(0x2512d43565724728), TO_LIMB_T(0x0aa6346091755d4d)
};
#endif // __CUDA_ARCH__
