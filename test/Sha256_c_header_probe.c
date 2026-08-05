/*
 * Copyright Category Labs and Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "bindings/blst_scalar.h"

void blst_sha256_c_header_probe(blst_byte out[32], const void* input,
    size_t input_len)
{
    blst_sha256(out, input, input_len);
}
