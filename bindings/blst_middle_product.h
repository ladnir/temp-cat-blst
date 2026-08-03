/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BLST_MIDDLE_PRODUCT_H
#define BLST_MIDDLE_PRODUCT_H

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
 * Computes the constant-term contribution against a precomputed `b`.
 *
 * Parameter interpretation:
 * - `a0` points to one G1 projective point in the scalar/public ABI layout.
 * - `precomputed_b` is the handle returned by
 *   `blst_middle_product_precompute_b`.
 * - `out` points to writable storage for `n` consecutive `blst_fp12` outputs,
 *   where `n` is the length captured by `precomputed_b`.
 *
 * Return value:
 * - `true` on success.
 * - `false` when any pointer is null, the handle is invalid, or the native
 *   computation throws/fails.
 *
 * Invariants:
 * - `out[k]` receives the pairing result for `a0` against the stored tail
 *   point `b[n + k]`.
 */
bool blst_middle_product_update(
    const blst_p1* a0,
    const blst_middle_product_precomputed_b* precomputed_b,
    blst_fp12* out);

/**
 * Combines two length-`n` middle-product vectors elementwise in GT.
 *
 * Parameter interpretation:
 * - `mid` and `cnst` point to `n` consecutive `blst_fp12` values in the
 *   scalar/public ABI layout.
 * - `out` points to writable storage for `n` combined outputs.
 *
 * Return value:
 * - `true` on success.
 * - `false` when any pointer is null or `n == 0`.
 *
 * Invariants:
 * - Combination uses the GT group law.
 * - At the `blst_fp12` storage boundary this is represented by elementwise
 *   multiplication, since GT is the cyclotomic subgroup inside Fp12.
 */
bool blst_middle_product_combine(
    const blst_fp12* mid,
    const blst_fp12* cnst,
    size_t n,
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
