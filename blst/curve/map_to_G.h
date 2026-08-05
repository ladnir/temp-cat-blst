/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include "blst/common/defines.h"
#if BLST_DEFAULT_AVX
#include "blst/avx/curve/map_to_G.h"
#else
#include "blst/scalar/curve/map_to_G.h"
#endif

namespace blst
{
#if BLST_DEFAULT_AVX
    using avx::Encode_to_curve;
    using avx::Hash_to_curve;
#else
    using scalar::Encode_to_curve;
    using scalar::Hash_to_curve;
#endif
}
