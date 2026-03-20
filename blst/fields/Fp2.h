/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_FIELDS_H__
#define __BLS12_381_ASM_FIELDS_H__

#include "blst/fields/Fp.h"
#include "blst/scalar/fields/Fp2.h"

namespace blst
{
    using scalar::Fp2;

    using scalar::BLST_sqrt_minus_1_HOST;
    using scalar::BLST_sqrt_minus_sqrt_minus_1_HOST;
    using scalar::BLST_sqrt_sqrt_minus_1_HOST;

    using scalar::cc_add_mod_384x;
    using scalar::cc_lshift_mod_384x;
    using scalar::cc_mul_by_1_plus_i_mod_384x;
    using scalar::cc_mul_by_3_mod_384x;
    using scalar::cc_mul_by_8_mod_384x;
    using scalar::cc_mul_mont_384x;
    using scalar::cc_sgn0_pty_mod_384x;
    using scalar::cc_sgn0_pty_mont_384x;
    using scalar::cc_sqr_mont_382x;
    using scalar::cc_sqr_mont_384x;
    using scalar::cc_sub_mod_384x;

    using scalar::get_BLST_sqrt_minus_1;
    using scalar::get_BLST_sqrt_minus_sqrt_minus_1;
    using scalar::get_BLST_sqrt_sqrt_minus_1;

    using scalar::mul_v382x;
    using scalar::sqr_v382x;
}

#endif /* __BLS12_381_ASM_FIELDS_H__ */
