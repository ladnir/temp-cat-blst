/*
 * Copyright Supranational LLC
 * Copyright Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BLST_FFT_H
#define BLST_FFT_H

#include "blst_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Runs the owned BLST reference FFT in-place over Montgomery-form Fr values.
 *
 * Parameter interpretation:
 * - `values` points to `n` consecutive reduced Montgomery-form `blst_fr`
 *   values.
 * - `n` must be a non-zero power of two no larger than `2^30`, the largest
 *   size supported by the signed-index radix-two engine.
 * - `inverse == false` evaluates with powers of the primitive `n`th root.
 * - `inverse == true` evaluates with reciprocal root powers and applies the
 *   usual `1/n` inverse-FFT scale.
 *
 * Return value:
 * - `true` on success.
 * - `false` when `values` is null, `n` is invalid, or the native computation
 *   throws/fails.
 *
 * Root preparation may allocate and fail before modifying `values`. Once the
 * in-place transform begins, an unexpected internal failure may leave values
 * partially transformed. The arithmetic schedule and addresses depend on the
 * public size and direction, not coefficient values. The function owns no
 * global mutable state; concurrent calls are safe when their mutable ranges
 * are disjoint.
 */
BLST_NODISCARD bool blst_fr_fft(blst_fr* values, size_t n, bool inverse);

/**
 * Returns the caller-owned `blst_fr` element count required for a prepared
 * radix-two FFT domain of size `n`.
 *
 * The prepared storage contains forward root powers, the inverse-size scale,
 * and an encoded size tag. A zero return means `n` is not a non-zero power of
 * two no larger than `2^30`. This query allocates nothing and cannot fail for
 * any other reason.
 */
BLST_NODISCARD size_t blst_fr_fft_domain_len(size_t n);

/**
 * Initializes caller-owned storage for a radix-two FFT domain of size `n`.
 *
 * `domain` must point to exactly `domain_len` writable `blst_fr` elements, and
 * `domain_len` must equal `blst_fr_fft_domain_len(n)`. The prepared values
 * remain in Montgomery form and may be reused by any number of forward or
 * inverse transforms of exactly `n` values. The function allocates nothing.
 * Publish or reuse the domain only after a `true` result; on `false`, its
 * contents are unspecified. Calls with disjoint mutable domain storage are
 * concurrency-safe.
 */
BLST_NODISCARD bool blst_fr_fft_domain_init(
    blst_fr* domain,
    size_t domain_len,
    size_t n);

/**
 * Runs the owned BLST FFT in place using a caller-prepared domain.
 *
 * The function allocates no memory and generates no roots. `values` points to
 * `n` writable Montgomery-form values. `domain` must have been initialized for
 * exactly `n` by `blst_fr_fft_domain_init`, and `domain_len` must exactly equal
 * `blst_fr_fft_domain_len(n)`; the two ranges must not overlap. The stored size
 * tag is checked before transforming. The inverse transform uses reciprocal
 * powers derived from the stored forward roots and applies the prepared `1/n`
 * scale. Every value and domain element must be a reduced Montgomery-form Fr
 * value. Invalid pointers, shapes, or size tags are rejected before modifying
 * `values`; an unexpected internal failure after transformation begins may
 * leave values partially transformed. Execution schedules and addresses depend
 * only on the public size and direction. An immutable prepared domain may be
 * shared concurrently, while callers synchronize or separate mutable value
 * buffers.
 */
BLST_NODISCARD bool blst_fr_fft_with_domain(
    blst_fr* values,
    size_t n,
    const blst_fr* domain,
    size_t domain_len,
    bool inverse);

/**
 * In-place radix-two FFT over projective G1 or G2 points. Inputs must already
 * be valid subgroup points, including the point at infinity. The transform is
 * linear over the group and uses the same roots and natural output order as
 * `blst_fr_fft`. Sizes must be nonzero powers of two through `2^30`.
 *
 * The one-shot calls allocate a root table. The prepared calls allocate
 * nothing and require an immutable domain initialized with
 * `blst_fr_fft_domain_init` for exactly `n`; domain and values must not
 * overlap. These functions use variable-time group multiplication and are
 * intended for public group points and roots, not secret points. Invalid
 * pointers, sizes, domain lengths, or domain tags return false before writing
 * values. A native failure after transformation begins may leave values
 * partially transformed. Concurrent calls need disjoint mutable values.
 */
BLST_NODISCARD bool blst_p1_fft(blst_p1* values, size_t n, bool inverse);
BLST_NODISCARD bool blst_p2_fft(blst_p2* values, size_t n, bool inverse);
BLST_NODISCARD bool blst_p1_fft_with_domain(
    blst_p1* values,
    size_t n,
    const blst_fr* domain,
    size_t domain_len,
    bool inverse);
BLST_NODISCARD bool blst_p2_fft_with_domain(
    blst_p2* values,
    size_t n,
    const blst_fr* domain,
    size_t domain_len,
    bool inverse);

/**
 * Inverse group FFT without the usual `1/n` normalization. The output is
 * `sum_j values[j] * root^(-j*k)` at each index `k`. This is useful when a
 * later linear transform applies one combined scale. Domain, subgroup,
 * aliasing, failure, and concurrency requirements match the prepared group
 * FFTs above. Neither call allocates.
 */
BLST_NODISCARD bool blst_p1_ifft_unscaled_with_domain(
    blst_p1* values, size_t n, const blst_fr* domain, size_t domain_len);
BLST_NODISCARD bool blst_p2_ifft_unscaled_with_domain(
    blst_p2* values, size_t n, const blst_fr* domain, size_t domain_len);

/**
 * Interpolates the unique degree-<n Fr polynomial through n distinct points.
 *
 * Parameter interpretation:
 * - `out_coefficients` points to n writable Fr values.
 * - `points` and `values` each point to n reduced Montgomery-form Fr values.
 * - On success, `out_coefficients[i]` is the monomial coefficient of X^i.
 *
 * Return value:
 * - `true` on success.
 * - `false` when any pointer is null, n is zero or unsupported, points are not
 *   distinct, or the native computation throws/fails.
 *
 * Output is unchanged on failure and may overlap either input. The operation
 * may allocate and is variable-time in the points and public length. For fixed
 * points its reconstruction schedule does not depend on values. Ordinary heap
 * storage contains value-derived intermediates and is not securely erased.
 * Calls own no global mutable state and are concurrency-safe when mutable
 * output ranges are disjoint.
 */
BLST_NODISCARD bool blst_fr_interpolate(
    blst_fr* out_coefficients,
    const blst_fr* points,
    const blst_fr* values,
    size_t n);

/**
 * Returns the caller-owned Fr scratch length required by the FFT-domain
 * coefficient-suffix interpolation operation below.
 *
 * The returned unit is `blst_fr` elements, not bytes. A zero return means the
 * shape is invalid, including when `n` exceeds `fft_len` or the calculation
 * would overflow `size_t`.
 */
BLST_NODISCARD size_t blst_fr_interpolate_fft_domain_suffix_scratch_len(
    size_t n,
    size_t fft_len);

/**
 * Interpolates a high coefficient suffix from FFT-domain samples.
 *
 * Given `n` values at FFT-domain indices and an offset `ell`, writes
 * coefficients `ell..n-1` of the unique degree-<n interpolant to `out_suffix`.
 * Every input value must be a reduced Montgomery-form Fr element.
 * The output must contain `n-ell` writable Fr values. `scratch` must point to at
 * least
 * `blst_fr_interpolate_fft_domain_suffix_scratch_len(n, fft_len)` writable Fr
 * elements. The scratch buffer is fully caller-owned and may be reused across
 * calls with the same or smaller shape. Scratch must not overlap the output,
 * indices, or values ranges. Invalid shapes, indices, duplicate points, and
 * short scratch are rejected before writing the output. Every pointer argument
 * must be non-null, including `out_suffix` when `ell == n` and its logical
 * output length is zero.
 * Output may overlap
 * values, but not indices or scratch. Output is unchanged on failure; scratch
 * may be modified after initial validation succeeds. The operation allocates
 * nothing and its schedule depends on public lengths and indices, not values.
 * Scratch contains value-derived intermediates and must be cleared by the
 * caller when that lifecycle policy is required.
 * Calls own no global mutable state and are concurrency-safe when mutable
 * output and scratch ranges are disjoint between calls.
 */
BLST_NODISCARD bool blst_fr_interpolate_fft_domain_suffix_with_scratch(
    blst_fr* out_suffix,
    const size_t* indices,
    const blst_fr* values,
    size_t n,
    size_t ell,
    size_t fft_len,
    blst_fr* scratch,
    size_t scratch_len);

#ifdef __cplusplus
}
#endif

#endif
