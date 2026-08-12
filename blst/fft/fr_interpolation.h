// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// CatBLST arbitrary-point and FFT-domain scalar-field interpolation.

#pragma once

#include "blst/fields/Fr.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>

namespace blst::fr_interpolation
{
    /** Native interpolation result. No exception crosses these entry points. */
    enum class Status : std::uint8_t
    {
        success = 0,
        invalid_size,
        invalid_shape,
        invalid_index,
        duplicate_point,
        insufficient_scratch,
        allocation_failure,
        internal_failure,
    };

    /**
     * Recover the n monomial coefficients of the unique degree-<n polynomial
     * through `(points[i], values[i])`.
     *
     * Preconditions: all three spans have the same nonzero supported length and
     * points are pairwise distinct. Output may overlap either input because it is
     * published only after the complete result is available. On failure output
     * is unchanged. The operation may allocate and is variable-time in the
     * public points and length; reconstruction has a fixed shape with respect
     * to values for a fixed point set.
     */
    [[nodiscard]] Status interpolate(
        std::span<Fr> out_coefficients,
        std::span<const Fr> points,
        std::span<const Fr> values) noexcept;

    /** Return required caller scratch in Fr elements, or zero if invalid. */
    [[nodiscard]] constexpr std::size_t domain_suffix_scratch_len(
        const std::size_t n,
        const std::size_t fft_len) noexcept
    {
        if (n == 0 || n > fft_len || fft_len == 0 ||
            (fft_len & (fft_len - 1)) != 0 ||
            fft_len > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            return 0;
        }

        constexpr std::size_t max = std::numeric_limits<std::size_t>::max();
        if (fft_len > (max - 2) / 2) {
            return 0;
        }
        const std::size_t fft_storage = 2 * fft_len;
        const std::size_t remaining = max - fft_storage - 2;
        if (n > remaining / 3) {
            return 0;
        }
        return fft_storage + 3 * n + 2;
    }

    /**
     * Recover coefficients `[ell,n)` of a degree-<n polynomial from values at
     * distinct positions in the ordered `fft_len` root-of-unity domain.
     *
     * Preconditions: `fft_len` is supported, `1 <= n <= fft_len`, indices and
     * values have length n, every index is below fft_len and unique, output has
     * length `n-ell`, and scratch has at least `domain_suffix_scratch_len`
     * elements. Scratch must not overlap any other argument. Output may overlap
     * values, but not indices, and is unchanged on failure; scratch is
     * unspecified after a call that passes initial shape validation. Execution
     * allocates nothing. Its schedule and addresses depend on public lengths and
     * indices, not on values.
     */
    [[nodiscard]] Status interpolate_domain_suffix_with_scratch(
        std::span<Fr> out_suffix,
        std::span<const std::size_t> indices,
        std::span<const Fr> values,
        std::size_t ell,
        std::size_t fft_len,
        std::span<Fr> scratch) noexcept;
}
