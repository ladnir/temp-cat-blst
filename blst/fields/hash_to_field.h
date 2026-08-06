/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_HASH_TO_FIELD_H
#define BLST_HASH_TO_FIELD_H

#include "blst/scalar/fields/hash_to_field.h"

namespace blst
{
    /**
     * Selected native C++ access to allocation-free RFC 9380 XMD:SHA-256.
     * See the owning scalar declaration for the exact pointer, length,
     * nonoverlap, maximum-output, and fail-closed contracts.
     */
    using scalar::blst_expand_message_xmd;
}

#endif
