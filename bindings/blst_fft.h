/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef __BLST_FFT_H__
#define __BLST_FFT_H__

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
 * Runs the transpose-equivalent of the forward Fr FFT in-place.
 *
 * This is the verifier-side BVE scalar kernel: callers scatter sparse
 * coefficients into a length-n vector and read the low coefficients after this
 * transform. The current BLST Fr FFT exposes the natural-order DFT matrix,
 * whose forward matrix is symmetric, so the native implementation can share the
 * same butterfly kernel while giving protocol code the correct operation name.
 */
bool blst_fr_fft_transpose_forward(blst_fr* values, size_t n);

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
 * - `false` when any pointer is null, n is zero, points are not distinct, or
 *   the native computation throws/fails.
 */
bool blst_fr_interpolate(
    blst_fr* out_coefficients,
    const blst_fr* points,
    const blst_fr* values,
    size_t n);

/**
 * Interpolates using points selected by index from a size-fft_len Fr FFT domain.
 *
 * `indices[i]` selects omega^indices[i], where omega is the primitive
 * fft_len-th root used by `blst_fr_fft`. This avoids materializing the domain
 * points at the Rust/protocol boundary. `n` must not exceed `fft_len`.
 */
bool blst_fr_interpolate_fft_domain(
    blst_fr* out_coefficients,
    const size_t* indices,
    const blst_fr* values,
    size_t n,
    size_t fft_len);

/**
 * Interpolates only the high coefficient suffix used by BVE receivers.
 *
 * Given n values at FFT-domain indices and an offset ell, writes coefficients
 * ell..n-1 of the unique degree-<n interpolant to out_suffix. The out buffer
 * must contain n-ell writable Fr values. `n` must not exceed `fft_len`.
 */
bool blst_fr_interpolate_fft_domain_suffix(
    blst_fr* out_suffix,
    const size_t* indices,
    const blst_fr* values,
    size_t n,
    size_t ell,
    size_t fft_len);

/**
 * Returns the caller-owned Fr scratch length required by
 * `blst_fr_interpolate_fft_domain_suffix_with_scratch`.
 *
 * The returned unit is `blst_fr` elements, not bytes. A zero return means the
 * shape is invalid, including when `n` exceeds `fft_len` or the calculation
 * would overflow `size_t`.
 */
size_t blst_fr_interpolate_fft_domain_suffix_scratch_len(size_t n, size_t fft_len);

/**
 * Scratch-backed variant of `blst_fr_interpolate_fft_domain_suffix`.
 *
 * `scratch` must point to at least
 * `blst_fr_interpolate_fft_domain_suffix_scratch_len(n, fft_len)` writable Fr
 * elements. The scratch buffer is fully caller-owned and may be reused across
 * calls with the same or smaller shape.
 */
bool blst_fr_interpolate_fft_domain_suffix_with_scratch(
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
