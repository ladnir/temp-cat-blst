/*
 * Copyright Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "bindings/blst_fft.h"

#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/domain.h"
#include "blst/fft/fr_interpolation.h"

#include <span>
#include <type_traits>
#include <vector>

namespace
{
    using NativeFr = blst::Fr;
    using InterpolationStatus = blst::fr_interpolation::Status;
    using NativeP1 = blst::scalar::G1Point;
    using NativeP2 = blst::scalar::G2Point;

    static_assert(std::is_same_v<blst_fr, NativeFr>);
    static_assert(std::is_standard_layout_v<NativeFr>);
    static_assert(std::is_trivially_copyable_v<NativeFr>);
    static_assert(std::is_same_v<blst_p1, NativeP1>);
    static_assert(std::is_same_v<blst_p2, NativeP2>);

    [[nodiscard]] bool succeeded(const InterpolationStatus status) noexcept
    {
        return status == InterpolationStatus::success;
    }

    std::span<NativeFr> native_span(
        blst_fr* const values,
        const std::size_t size) noexcept
    {
        return { values, size };
    }

    std::span<const NativeFr> native_span(
        const blst_fr* const values,
        const std::size_t size) noexcept
    {
        return { values, size };
    }

    template<class Group>
    bool group_fft(Group* values, const std::size_t n, const bool inverse)
    {
        if (values == nullptr || !blst::fft_detail::valid_size(n)) {
            return false;
        }
        try {
            const auto roots = blst::fft_detail::make_roots<
                NativeFr, blst::FrCoeffCtx>(n, inverse);
            fft<Group, NativeFr, blst::detail::GroupFftCoeffCtx<Group>>(
                { values, n }, roots, inverse);
            return true;
        } catch (...) {
            return false;
        }
    }

    template<class Group, bool UnscaledInverse = false>
    bool group_fft_with_domain(
        Group* values,
        const std::size_t n,
        const blst_fr* domain,
        const std::size_t domain_len,
        const bool inverse)
    {
        const std::size_t required = blst::fft_detail::domain_storage_len(n);
        if (values == nullptr || domain == nullptr || required == 0 ||
            domain_len != required) {
            return false;
        }
        try {
            const auto native_domain = native_span(domain, domain_len);
            if (!NativeFr::equal(
                    native_domain[n + 1], NativeFr::from_uint64(n))) {
                return false;
            }
            if constexpr (UnscaledInverse) {
                (void)inverse;
                fft_butterflies<
                    true, Group, NativeFr,
                    blst::detail::GroupFftCoeffCtx<Group>>(
                    { values, n }, native_domain.first(n));
            } else {
                blst::fft_detail::transform_with_domain<
                    Group, NativeFr, blst::detail::GroupFftCoeffCtx<Group>>(
                    { values, n }, native_domain, inverse);
            }
            return true;
        } catch (...) {
            return false;
        }
    }
}

extern "C" {

// No exception may cross this C ABI. Allocation-bearing FFT functions catch
// locally; interpolation reaches only native noexcept status entry points.

bool blst_fr_fft(blst_fr* values, size_t n, bool inverse)
{
    if (values == nullptr || !blst::fft_detail::valid_size(n)) {
        return false;
    }

    try {
        const std::vector<NativeFr> roots =
            blst::fft_detail::make_roots<NativeFr, blst::FrCoeffCtx>(n, inverse);
        fft<NativeFr, NativeFr, blst::FrCoeffCtx>(
            native_span(values, n),
            roots,
            inverse);
        return true;
    } catch (...) {
        return false;
    }
}

size_t blst_fr_fft_domain_len(size_t n)
{
    return blst::fft_detail::domain_storage_len(n);
}

bool blst_fr_fft_domain_init(blst_fr* domain, size_t domain_len, size_t n)
{
    const std::size_t required = blst::fft_detail::domain_storage_len(n);
    if (domain == nullptr || required == 0 || domain_len != required) {
        return false;
    }

    try {
        blst::fft_detail::initialize_domain<NativeFr, blst::FrCoeffCtx>(
            native_span(domain, domain_len),
            n);
        return true;
    } catch (...) {
        return false;
    }
}

bool blst_fr_fft_with_domain(
    blst_fr* values,
    size_t n,
    const blst_fr* domain,
    size_t domain_len,
    bool inverse)
{
    const std::size_t required = blst::fft_detail::domain_storage_len(n);
    if (values == nullptr || domain == nullptr || required == 0 ||
        domain_len != required) {
        return false;
    }

    try {
        const std::span<const NativeFr> native_domain =
            native_span(domain, domain_len);
        const NativeFr expected_size_tag = NativeFr::from_uint64(n);
        if (!NativeFr::equal(native_domain[n + 1], expected_size_tag)) {
            return false;
        }

        blst::fft_detail::transform_with_domain<
            NativeFr,
            NativeFr,
            blst::FrCoeffCtx>(
                native_span(values, n),
                native_domain,
                inverse);
        return true;
    } catch (...) {
        return false;
    }
}

bool blst_p1_fft(blst_p1* values, size_t n, bool inverse)
{
    return group_fft(values, n, inverse);
}

bool blst_p2_fft(blst_p2* values, size_t n, bool inverse)
{
    return group_fft(values, n, inverse);
}

bool blst_p1_fft_with_domain(
    blst_p1* values, size_t n, const blst_fr* domain,
    size_t domain_len, bool inverse)
{
    return group_fft_with_domain(values, n, domain, domain_len, inverse);
}

bool blst_p2_fft_with_domain(
    blst_p2* values, size_t n, const blst_fr* domain,
    size_t domain_len, bool inverse)
{
    return group_fft_with_domain(values, n, domain, domain_len, inverse);
}

bool blst_p1_ifft_unscaled_with_domain(
    blst_p1* values, size_t n, const blst_fr* domain, size_t domain_len)
{
    return group_fft_with_domain<NativeP1, true>(
        values, n, domain, domain_len, true);
}

bool blst_p2_ifft_unscaled_with_domain(
    blst_p2* values, size_t n, const blst_fr* domain, size_t domain_len)
{
    return group_fft_with_domain<NativeP2, true>(
        values, n, domain, domain_len, true);
}

bool blst_fr_interpolate(
    blst_fr* out_coefficients,
    const blst_fr* points,
    const blst_fr* values,
    size_t n)
{
    return out_coefficients != nullptr && points != nullptr && values != nullptr &&
        succeeded(blst::fr_interpolation::interpolate(
            native_span(out_coefficients, n),
            native_span(points, n),
            native_span(values, n)));
}

size_t blst_fr_interpolate_fft_domain_suffix_scratch_len(size_t n, size_t fft_len)
{
    return blst::fr_interpolation::domain_suffix_scratch_len(n, fft_len);
}

bool blst_fr_interpolate_fft_domain_suffix_with_scratch(
    blst_fr* out_suffix,
    const size_t* indices,
    const blst_fr* values,
    size_t n,
    size_t ell,
    size_t fft_len,
    blst_fr* scratch,
    size_t scratch_len)
{
    if (out_suffix == nullptr || indices == nullptr || values == nullptr ||
        scratch == nullptr || ell > n) {
        return false;
    }

    return succeeded(
        blst::fr_interpolation::interpolate_domain_suffix_with_scratch(
            native_span(out_suffix, n - ell),
            std::span<const std::size_t>(indices, n),
            native_span(values, n),
            ell,
            fft_len,
            native_span(scratch, scratch_len)));
}

}
