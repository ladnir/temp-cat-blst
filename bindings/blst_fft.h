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
 * points at the Rust/protocol boundary.
 */
bool blst_fr_interpolate_fft_domain(
    blst_fr* out_coefficients,
    const size_t* indices,
    const blst_fr* values,
    size_t n,
    size_t fft_len);

#ifdef __cplusplus
}
#endif

#endif
