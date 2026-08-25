/*
 * Copyright Supranational LLC
 * Copyright 2026 Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 *
 * CatBLST middle-product C API, built from BLST pairing and FFT primitives.
 */

#ifndef BLST_MIDDLE_PRODUCT_H
#define BLST_MIDDLE_PRODUCT_H

#include "blst_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct blst_middle_product_precomputed_b blst_middle_product_precomputed_b;
typedef struct blst_middle_product_precomputed_b_scalar
    blst_middle_product_precomputed_b_scalar;
typedef struct blst_middle_product_precomputed_b_avx
    blst_middle_product_precomputed_b_avx;

/**
 * Explicit scalar and AVX middle-product handles.
 *
 * These APIs accept and produce the same scalar/public ABI values as the
 * default interface, but the opaque precompute fixes the arithmetic backend
 * at construction. Scalar and AVX handles are different C types and are not
 * interchangeable. The AVX constructor returns null, and AVX evaluation
 * returns false, when CatBLST was built without AVX support.
 */
BLST_NODISCARD blst_middle_product_precomputed_b_scalar*
blst_middle_product_precompute_b_scalar(const blst_p2* b, size_t n);
BLST_NODISCARD bool blst_middle_product_scalar(
    const blst_p1* a,
    size_t n,
    const blst_middle_product_precomputed_b_scalar* precomputed_b,
    blst_fp12* out);
BLST_NODISCARD bool blst_middle_product_update_scalar(
    const blst_p1* a0,
    const blst_middle_product_precomputed_b_scalar* precomputed_b,
    blst_fp12* out);
void blst_middle_product_precomputed_b_scalar_free(
    blst_middle_product_precomputed_b_scalar* precomputed_b);

BLST_NODISCARD blst_middle_product_precomputed_b_avx*
blst_middle_product_precompute_b_avx(const blst_p2* b, size_t n);
BLST_NODISCARD bool blst_middle_product_avx(
    const blst_p1* a,
    size_t n,
    const blst_middle_product_precomputed_b_avx* precomputed_b,
    blst_fp12* out);
BLST_NODISCARD bool blst_middle_product_update_avx(
    const blst_p1* a0,
    const blst_middle_product_precomputed_b_avx* precomputed_b,
    blst_fp12* out);
void blst_middle_product_precomputed_b_avx_free(
    blst_middle_product_precomputed_b_avx* precomputed_b);

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
 * - Null when `b` is null, `n == 0`, `n` exceeds the supported native FFT
 *   range, allocation fails, or precomputation fails.
 *
 * Invariants:
 * - Every input point must already be a reduced, on-curve member of the
 *   prime-order G2 subgroup. Infinity is accepted as the identity. This
 *   function performs no Release point or subgroup validation.
 * - The returned handle owns native resources and must be released with
 *   `blst_middle_product_precomputed_b_free`.
 * - A successfully constructed handle is immutable. Concurrent evaluations
 *   may share it, but freeing it must not race with an evaluation.
 * - Precomputation allocates and is variable-time with respect to the public
 *   points and length. Native exceptions, unsupported lengths, and allocation
 *   failures are translated to null.
 * - The Rust-facing ABI never exposes the internal layout of the precompute.
 */
BLST_NODISCARD blst_middle_product_precomputed_b* blst_middle_product_precompute_b(
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
 * - Every input point must already be a reduced, on-curve member of the
 *   prime-order G1 subgroup. Infinity is accepted as the identity. No Release
 *   point or subgroup validation is performed.
 * - `out[i]` receives the `i`th GT middle-product result in scalar/public
 *   `blst_fp12` storage.
 * - Specifically, `out[k] = product_i e(a[i], b[n+k-i])` for
 *   `0 <= i,k < n`, where `b` is the operand captured by `precomputed_b`.
 * - The function is variable-time with respect to the public points and
 *   length and may allocate on every call.
 * - On failure no element of `out` is modified.
 * - `out` must not overlap `a` or storage reachable through `precomputed_b`.
 * - A non-null `precomputed_b` must be a live handle returned by this API.
 *   Arbitrary, freed, or concurrently freed pointers are outside the contract
 *   and cannot be validated by the C boundary.
 * - The implementation may internally convert the scalar/public inputs to AVX
 *   values, but the exposed ABI remains scalar/public.
 */
BLST_NODISCARD bool blst_middle_product(
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
 * - `false` when any pointer is null or the native computation throws/fails.
 *
 * Invariants:
 * - `a0` must already be a reduced, on-curve member of the prime-order G1
 *   subgroup. Infinity is accepted as the identity. No Release validation is
 *   performed.
 * - `out[k]` receives the GT pairing result for `a0` against the stored tail
 *   point `b[n + k]`, represented in `blst_fp12` storage.
 * - The function is variable-time with respect to the public point and may
 *   allocate on every call.
 * - On failure no element of `out` is modified.
 * - `out` must not overlap `a0` or storage reachable through `precomputed_b`.
 * - A non-null `precomputed_b` must be a live handle returned by this API.
 *   Arbitrary, freed, or concurrently freed pointers are outside the contract.
 */
BLST_NODISCARD bool blst_middle_product_update(
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
 * - Every input must be a reduced GT element. No Release membership validation
 *   is performed.
 * - Exact aliasing (`out == mid` or `out == cnst`) is supported. Partial
 *   overlap is forbidden.
 * - The operation allocates nothing and is nonthrowing for valid inputs.
 */
BLST_NODISCARD bool blst_middle_product_combine(
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
 * - The call must not race with another use of the same handle.
 */
void blst_middle_product_precomputed_b_free(
    blst_middle_product_precomputed_b* precomputed_b);

#ifdef __cplusplus
}
#endif

#endif
