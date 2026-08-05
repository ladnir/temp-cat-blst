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
    /// Selected RFC 9380 nonuniform G1/G2 encoder; see the AVX declaration for
    /// pointer, validation, allocation, exception, and subgroup contracts.
    using avx::Encode_to_curve;
    /// Selected RFC 9380 random-oracle G1/G2 hash; its contract matches the
    /// selected Encode_to_curve declaration.
    using avx::Hash_to_curve;
#else
    /// Selected RFC 9380 nonuniform G1/G2 encoder; see the scalar declaration
    /// for pointer, validation, allocation, exception, and subgroup contracts.
    using scalar::Encode_to_curve;
    /// Selected RFC 9380 random-oracle G1/G2 hash; its contract matches the
    /// selected Encode_to_curve declaration.
    using scalar::Hash_to_curve;
#endif
}
