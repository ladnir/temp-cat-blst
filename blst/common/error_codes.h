/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BLST_ERROR_CODES_H
#define BLST_ERROR_CODES_H

/**
 * Stable error codes shared by the core implementation and the plain C API.
 *
 * Invariants:
 * - Numeric values are part of the public ABI.
 * - The scalar, AVX, and middle-product APIs use the same values.
 */
typedef enum {
    BLST_SUCCESS = 0,
    BLST_BAD_ENCODING = 1,
    BLST_POINT_NOT_ON_CURVE = 2,
    BLST_POINT_NOT_IN_GROUP = 3,
    BLST_AGGR_TYPE_MISMATCH = 4,
    BLST_VERIFY_FAIL = 5,
    BLST_PK_IS_INFINITY = 6,
    BLST_BAD_SCALAR = 7,
} BLST_ERROR;

#endif /* BLST_ERROR_CODES_H */
