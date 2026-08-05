/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_FIELDS_FP12_H
#define BLST_FIELDS_FP12_H

#include "blst/fields/Fp2.h"
#include "blst/scalar/fields/Fp12.h"

namespace blst
{
    using scalar::Fp4;
    using scalar::Fp6;
    /**
     * @brief Main native C++ BLS12-381 degree-12 extension-field type.
     *
     * This is an exact alias of `scalar::Fp12`. General Fp12, cyclotomic, and
     * GT-only method domains remain distinct as documented on that type.
     */
    using scalar::Fp12;
    using scalar::GLV4_Exp;
    using scalar::GLV_JSF_Exp;
    using scalar::BLST_Rx_HOST;
    using scalar::get_BLST_Rx;
}

#endif /* BLST_FIELDS_FP12_H */
