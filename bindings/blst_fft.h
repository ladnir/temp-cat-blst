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
 * - `values` points to `n` consecutive `blst_fr` values.
 * - `n` must be a non-zero power of two and fit within the BLS12-381 scalar
 *   field 2-adicity.
 * - `inverse == false` evaluates with powers of the primitive `n`th root.
 * - `inverse == true` evaluates with reciprocal root powers and applies the
 *   usual `1/n` inverse-FFT scale.
 *
 * Return value:
 * - `true` on success.
 * - `false` when `values` is null, `n` is invalid, or the native computation
 *   throws/fails.
 */
bool blst_fr_fft(blst_fr* values, size_t n, bool inverse);

/**
 * Returns the caller-owned `blst_fr` element count required for a prepared
 * radix-two FFT domain of size `n`.
 *
 * The prepared storage contains forward root powers, the inverse-size scale,
 * and an encoded size tag. A zero return means `n` is not a supported non-zero
 * power of two.
 */
size_t blst_fr_fft_domain_len(size_t n);

/**
 * Initializes caller-owned storage for a radix-two FFT domain of size `n`.
 *
 * `domain` must point to exactly `domain_len` writable `blst_fr` elements, and
 * `domain_len` must equal `blst_fr_fft_domain_len(n)`. The prepared values
 * remain in Montgomery form and may be reused by any number of forward or
 * inverse transforms of exactly `n` values.
 */
bool blst_fr_fft_domain_init(
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
 * scale.
 */
bool blst_fr_fft_with_domain(
    blst_fr* values,
    size_t n,
    const blst_fr* domain,
    size_t domain_len,
    bool inverse);

/**
 * Interpolates the unique degree-<n Fr polynomial through n distinct points.
 *
 * Parameter interpretation:
 * - `out_coefficients` points to n writable Fr values.
 * - `points` and `values` each point to n Fr values.
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
 * The output must contain `n-ell` writable Fr values. `scratch` must point to at
 * least
 * `blst_fr_interpolate_fft_domain_suffix_scratch_len(n, fft_len)` writable Fr
 * elements. The scratch buffer is fully caller-owned and may be reused across
 * calls with the same or smaller shape. Scratch must not overlap the output,
 * indices, or values ranges. Invalid shapes, indices, duplicate points, and
 * short scratch are rejected before writing the output. Output may overlap
 * values, but not indices or scratch. Output is unchanged on failure; scratch
 * may be modified after initial validation succeeds. The operation allocates
 * nothing and its schedule depends on public lengths and indices, not values.
 * Scratch contains value-derived intermediates and must be cleared by the
 * caller when that lifecycle policy is required.
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
