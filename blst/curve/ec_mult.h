/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_CURVE_EC_MULT_H
#define BLST_CURVE_EC_MULT_H

#include "blst/scalar/curve/ec_mult.h"

namespace blst
{
    using scalar::booth_encode;
    using scalar::gather_booth_w;
    using scalar::get_wval;
    using scalar::get_wval_limb;
    using scalar::mult_window;
    using scalar::mult_window_n;
    using scalar::precompute_w;
}

#endif
