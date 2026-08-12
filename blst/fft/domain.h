// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// CatBLST FFT-domain construction helpers.

#pragma once

#include "blst/fft/ctFft.h"
#include "blst/fft/primitveRoot.h"

#include <cstddef>
#include <span>
#include <vector>

namespace blst::fft_detail
{
    [[nodiscard]] constexpr bool valid_size(const std::size_t n) noexcept
    {
        return valid_transform_size(n);
    }

    [[nodiscard]] constexpr std::size_t domain_storage_len(
        const std::size_t n) noexcept
    {
        return valid_size(n) ? n + 2 : 0;
    }

    template<typename F, typename CoeffCtx>
    void fill_roots(const std::span<F> roots, const bool inverse)
    {
        F primitive = primRootOfUnity<F>(roots.size());
        if (inverse) {
            primitive = CoeffCtx::inverse(primitive);
        }

        F current = CoeffCtx::template from_u64<F>(1);
        for (F& root : roots) {
            root = current;
            current = CoeffCtx::mult(current, primitive);
        }
    }

    template<typename F, typename CoeffCtx>
    [[nodiscard]] std::vector<F> make_roots(
        const std::size_t n,
        const bool inverse)
    {
        std::vector<F> roots(n);
        fill_roots<F, CoeffCtx>(roots, inverse);
        return roots;
    }

    // Domain layout: n forward roots, 1/n, and an n tag. The tag is checked at
    // typed ABI boundaries; generic internal callers are responsible for using
    // a domain initialized for the transform length.
    template<typename F, typename CoeffCtx>
    void initialize_domain(const std::span<F> domain, const std::size_t n)
    {
        fill_roots<F, CoeffCtx>(domain.first(n), false);
        domain[n] = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(n));
        domain[n + 1] = CoeffCtx::template from_u64<F>(n);
    }

    template<typename G, typename F, typename CoeffCtx>
    void transform_with_domain(
        const std::span<G> values,
        const std::span<const F> domain,
        const bool inverse)
    {
        fft_with_forward_roots<G, F, CoeffCtx>(
            values,
            domain.first(values.size()),
            inverse,
            domain[values.size()]);
    }
}
