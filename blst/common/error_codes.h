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
 * - The scalar and AVX APIs use the same values.
 * - New values may be appended, so FFI consumers must preserve unknown
 *   four-byte results rather than model this as a permanently closed enum.
 */
typedef enum {
    /** Operation completed successfully. */
    BLST_SUCCESS = 0,
    /** Byte input is malformed or noncanonical. */
    BLST_BAD_ENCODING = 1,
    /** Decoded coordinates do not satisfy the curve equation. */
    BLST_POINT_NOT_ON_CURVE = 2,
    /** An on-curve point is outside the required prime-order subgroup. */
    BLST_POINT_NOT_IN_GROUP = 3,
    /** Pairing accumulator orientations are incompatible. */
    BLST_AGGR_TYPE_MISMATCH = 4,
    /** A verification equation did not hold. */
    BLST_VERIFY_FAIL = 5,
    /** A public key that must be finite was infinity. */
    BLST_PK_IS_INFINITY = 6,
    /** A scalar pointer, width, or value violated the entry-point contract. */
    BLST_BAD_SCALAR = 7,
} BLST_ERROR;

/* The Rust FFI uses a transparent 32-bit integer wrapper for forward-safe
 * handling of error values added by future compatible releases. */
typedef char blst_error_code_must_be_4_bytes[
    sizeof(BLST_ERROR) == 4 ? 1 : -1];

#endif /* BLST_ERROR_CODES_H */
