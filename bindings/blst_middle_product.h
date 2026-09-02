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

#include "blst_avx.h"
#include "blst_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct blst_middle_product_precomputed_b blst_middle_product_precomputed_b;
typedef struct blst_middle_product_precomputed_b_scalar
    blst_middle_product_precomputed_b_scalar;
typedef struct blst_middle_product_precomputed_b_avx
    blst_middle_product_precomputed_b_avx;
typedef struct blst_middle_product_workspace_scalar
    blst_middle_product_workspace_scalar;
typedef struct blst_middle_product_workspace_avx
    blst_middle_product_workspace_avx;
typedef struct blst_middle_product_workspace_avx_native
    blst_middle_product_workspace_avx_native;

/**
 * Stateless range phases for one middle-product workspace.
 *
 * `INPUT_LOAD` and `EXTRACTION` use their typed input/output operations.
 * `FORWARD_J0`, `FORWARD_NONZERO`, and `INVERSE_BUTTERFLIES` use a zero-based
 * FFT stage. Every other phase requires `stage == 0`. The work-count query
 * defines the phase-local `[0, count)` work-item domain.
 */
typedef enum blst_middle_product_range_phase {
    BLST_MIDDLE_PRODUCT_RANGE_INPUT_LOAD = 0,
    BLST_MIDDLE_PRODUCT_RANGE_FORWARD_BIT_REVERSE = 1,
    BLST_MIDDLE_PRODUCT_RANGE_FORWARD_J0 = 2,
    BLST_MIDDLE_PRODUCT_RANGE_FORWARD_NONZERO = 3,
    BLST_MIDDLE_PRODUCT_RANGE_AFFINE_CONVERSION = 4,
    BLST_MIDDLE_PRODUCT_RANGE_POINTWISE_PAIRING = 5,
    BLST_MIDDLE_PRODUCT_RANGE_INVERSE_BIT_REVERSE = 6,
    BLST_MIDDLE_PRODUCT_RANGE_INVERSE_BUTTERFLIES = 7,
    BLST_MIDDLE_PRODUCT_RANGE_NORMALIZATION = 8,
    BLST_MIDDLE_PRODUCT_RANGE_EXTRACTION = 9,
} blst_middle_product_range_phase;

/**
 * Stateless scalar range execution.
 *
 * The precompute is immutable and may be shared. The workspace owns fixed,
 * 64-byte-aligned scalar buffers and must outlive every call. Concurrent calls
 * may share one workspace only when they execute disjoint work-item ranges of
 * the same phase and stage; the caller supplies every stage and phase barrier.
 * Input, precompute, and output must describe the same logical `n`. A workspace
 * created for a larger precompute may be reused here; its FFT capacity must be
 * at least the required precompute size. Input and output must not alias
 * workspace or precompute storage.
 *
 * `*_range_work_count` returns the constant-time work count for one valid
 * `(phase, stage)`. Input work items cover the full padded FFT buffer: indices
 * below `n` copy an admitted G1 point and the remainder write identity.
 * `*_execute_range` accepts only non-input/non-extraction phases. Output work
 * items cover the logical result `[0, n)`.
 *
 * Empty ranges are successful no-ops. Invalid phases, stages, ranges, lengths,
 * handles, or workspace capacity return false before the addressed operation
 * writes. A false output-range call does not publish external output. Calls
 * allocate nothing and are variable-time only with respect to public points,
 * sizes, and schedules. CatBLST retains no cursor, task, cancellation,
 * priority, callback, thread-pool, or phase-order state. Executing phases out
 * of order, overlapping ranges concurrently, or crossing a stage without a
 * barrier violates the contract and leaves workspace contents unspecified.
 */
BLST_NODISCARD bool blst_middle_product_scalar_range_work_count(
    const blst_middle_product_precomputed_b_scalar* precomputed_b,
    blst_middle_product_range_phase phase,
    size_t stage,
    size_t* count);
BLST_NODISCARD bool blst_middle_product_scalar_input_range(
    const blst_p1* a,
    size_t n,
    const blst_middle_product_precomputed_b_scalar* precomputed_b,
    blst_middle_product_workspace_scalar* workspace,
    size_t begin,
    size_t end);
BLST_NODISCARD bool blst_middle_product_scalar_execute_range(
    const blst_middle_product_precomputed_b_scalar* precomputed_b,
    blst_middle_product_workspace_scalar* workspace,
    blst_middle_product_range_phase phase,
    size_t stage,
    size_t begin,
    size_t end);
BLST_NODISCARD bool blst_middle_product_scalar_output_range(
    const blst_middle_product_precomputed_b_scalar* precomputed_b,
    blst_middle_product_workspace_scalar* workspace,
    blst_fp12* out,
    size_t n,
    size_t begin,
    size_t end);

/**
 * Conversion-free native-AVX range execution.
 *
 * This is the native-AVX analogue of the scalar range contract above. Input
 * and output use `blst_avx_p1` and `blst_avx_gt` directly; no canonical/scalar
 * representation appears inside a range call. The API returns false when AVX
 * support is unavailable. Native values and opaque handles must have been
 * created by the matching CatBLST AVX API; alignment is carried by their types
 * and the opaque 64-byte-aligned workspace allocation.
 */
BLST_NODISCARD bool blst_middle_product_avx_native_range_work_count(
    const blst_middle_product_precomputed_b_avx* precomputed_b,
    blst_middle_product_range_phase phase,
    size_t stage,
    size_t* count);
BLST_NODISCARD bool blst_middle_product_avx_native_input_range(
    const blst_avx_p1* a,
    size_t n,
    const blst_middle_product_precomputed_b_avx* precomputed_b,
    blst_middle_product_workspace_avx_native* workspace,
    size_t begin,
    size_t end);
BLST_NODISCARD bool blst_middle_product_avx_native_execute_range(
    const blst_middle_product_precomputed_b_avx* precomputed_b,
    blst_middle_product_workspace_avx_native* workspace,
    blst_middle_product_range_phase phase,
    size_t stage,
    size_t begin,
    size_t end);
BLST_NODISCARD bool blst_middle_product_avx_native_output_range(
    const blst_middle_product_precomputed_b_avx* precomputed_b,
    blst_middle_product_workspace_avx_native* workspace,
    blst_avx_gt* out,
    size_t n,
    size_t begin,
    size_t end);

/**
 * Explicit scalar and AVX middle-product handles.
 *
 * These APIs accept and produce the same scalar/public ABI values as the
 * default interface, but the opaque precompute fixes the arithmetic backend
 * at construction. Scalar and AVX handles are different C types and are not
 * interchangeable. The AVX constructor returns null, and AVX evaluation
 * returns false, when CatBLST was built without AVX support.
 *
 * The `*_workspace_*_new` constructors allocate all mutable FFT, output
 * staging, and (for AVX) input-conversion storage required by the matching
 * precompute. The `*_with_workspace` evaluations allocate nothing. A
 * workspace is mutable and must not be used concurrently; create one per
 * evaluation thread. Its length and backend must match the precompute.
 */
BLST_NODISCARD blst_middle_product_precomputed_b_scalar*
blst_middle_product_precompute_b_scalar(const blst_p2* b, size_t n);
BLST_NODISCARD bool blst_middle_product_scalar(
    const blst_p1* a,
    size_t n,
    const blst_middle_product_precomputed_b_scalar* precomputed_b,
    blst_fp12* out);
BLST_NODISCARD blst_middle_product_workspace_scalar*
blst_middle_product_workspace_scalar_new(
    const blst_middle_product_precomputed_b_scalar* precomputed_b);
BLST_NODISCARD bool blst_middle_product_scalar_with_workspace(
    const blst_p1* a,
    size_t n,
    const blst_middle_product_precomputed_b_scalar* precomputed_b,
    blst_middle_product_workspace_scalar* workspace,
    blst_fp12* out);
void blst_middle_product_workspace_scalar_free(
    blst_middle_product_workspace_scalar* workspace);
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
BLST_NODISCARD blst_middle_product_workspace_avx*
blst_middle_product_workspace_avx_new(
    const blst_middle_product_precomputed_b_avx* precomputed_b);
BLST_NODISCARD bool blst_middle_product_avx_with_workspace(
    const blst_p1* a,
    size_t n,
    const blst_middle_product_precomputed_b_avx* precomputed_b,
    blst_middle_product_workspace_avx* workspace,
    blst_fp12* out);
void blst_middle_product_workspace_avx_free(
    blst_middle_product_workspace_avx* workspace);
BLST_NODISCARD bool blst_middle_product_update_avx(
    const blst_p1* a0,
    const blst_middle_product_precomputed_b_avx* precomputed_b,
    blst_fp12* out);
void blst_middle_product_precomputed_b_avx_free(
    blst_middle_product_precomputed_b_avx* precomputed_b);

/**
 * Conversion-free AVX middle-product boundary.
 *
 * All G1, G2, and GT values use the native AVX ABI from `blst_avx.h`. The
 * precompute and workspace retain native storage, and evaluation writes
 * directly into caller-owned native GT output. No scalar representation or
 * scalar/native conversion is used by these functions.
 *
 * The precompute owns immutable data and may be shared between threads. The
 * workspace is mutable, allocation-free after construction, and must remain
 * exclusive to one evaluation at a time. The input, precompute, workspace,
 * and output lengths must match exactly.
 */
BLST_NODISCARD blst_middle_product_precomputed_b_avx*
blst_middle_product_precompute_b_avx_native(const blst_avx_p2* b, size_t n);
BLST_NODISCARD blst_middle_product_workspace_avx_native*
blst_middle_product_workspace_avx_native_new(
    const blst_middle_product_precomputed_b_avx* precomputed_b);
BLST_NODISCARD bool blst_middle_product_avx_native_with_workspace(
    const blst_avx_p1* a,
    size_t n,
    const blst_middle_product_precomputed_b_avx* precomputed_b,
    blst_middle_product_workspace_avx_native* workspace,
    blst_avx_gt* out);
void blst_middle_product_workspace_avx_native_free(
    blst_middle_product_workspace_avx_native* workspace);

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
