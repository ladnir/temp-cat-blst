#pragma once

#include "blst/fields/Fr.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>

namespace blst::fr_interpolation
{
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

    [[nodiscard]] Status interpolate(
        std::span<Fr> out_coefficients,
        std::span<const Fr> points,
        std::span<const Fr> values) noexcept;

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
        if (fft_len > (max - 2) / 3) {
            return 0;
        }
        const std::size_t fft_storage = 3 * fft_len;
        const std::size_t remaining = max - fft_storage - 2;
        if (n > remaining / 3) {
            return 0;
        }
        return fft_storage + 3 * n + 2;
    }

    [[nodiscard]] Status interpolate_domain_suffix_with_scratch(
        std::span<Fr> out_suffix,
        std::span<const std::size_t> indices,
        std::span<const Fr> values,
        std::size_t ell,
        std::size_t fft_len,
        std::span<Fr> scratch) noexcept;
}
