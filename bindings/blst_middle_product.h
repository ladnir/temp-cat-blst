/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef __BLST_MIDDLE_PRODUCT_H__
#define __BLST_MIDDLE_PRODUCT_H__

#include "blst_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct blst_middle_product_precomputed_b blst_middle_product_precomputed_b;

/**
 * Precomputes the `b` operand for repeated middle-product evaluations.
 *
 * Parameter interpretation:
 * - `b` points to `2 * n` consecutive G2 projective points in the scalar/public
 *   ABI layout.
 * - `n` is the logical output length of the middle product.
 *
 * Return value:
 * - Non-null opaque handle on success.
 * - Null when `b` is null, `n == 0`, allocation fails, or precomputation fails.
 *
 * Invariants:
 * - The returned handle owns native resources and must be released with
 *   `blst_middle_product_precomputed_b_free`.
 * - The Rust-facing ABI never exposes the internal layout of the precompute.
 */
blst_middle_product_precomputed_b* blst_middle_product_precompute_b(
    const blst_p2* b,
    size_t n);

/**
 * Evaluates the middle product using a precomputed `b` operand.
 *
 * Parameter interpretation:
 * - `a` points to `n` consecutive G1 projective points in the scalar/public ABI
 *   layout.
 * - `n` must match the length used to create `precomputed_b`.
 * - `precomputed_b` is the handle returned by
 *   `blst_middle_product_precompute_b`.
 * - `out` points to writable storage for `n` consecutive `blst_fp12` outputs.
 *
 * Return value:
 * - `true` on success.
 * - `false` when any pointer is null, `n == 0`, the precompute length does not
 *   match `n`, or the native computation throws/fails.
 *
 * Invariants:
 * - `out[i]` receives the `i`th middle-product result in scalar/public Fp12
 *   layout.
 * - The implementation may internally convert the scalar/public inputs to AVX
 *   values, but the exposed ABI remains scalar/public.
 */
bool blst_middle_product(
    const blst_p1* a,
    size_t n,
    const blst_middle_product_precomputed_b* precomputed_b,
    blst_fp12* out);

/**
 * Releases a precomputed middle-product handle.
 *
 * Parameters:
 * - `precomputed_b` may be null.
 *
 * Return value:
 * - None.
 *
 * Invariants:
 * - Passing null is a no-op.
 * - The handle must not be used after this call returns.
 */
void blst_middle_product_precomputed_b_free(
    blst_middle_product_precomputed_b* precomputed_b);

#ifdef __cplusplus
}
#endif

#endif
