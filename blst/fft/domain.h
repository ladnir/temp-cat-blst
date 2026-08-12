// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// CatBLST FFT-domain construction helpers.

#pragma once

#include "blst/fft/ctFft.h"
#include "blst/fft/primitive_root.h"

#include <cstddef>
#include <span>
#include <stdexcept>
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
        // Narrow native helper: supported ABI adapters validate before entry.
        BLST_ASSERT(valid_size(roots.size()));
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
        if (!valid_size(n)) {
            throw std::invalid_argument(
                "FFT root count must be a supported nonzero power of two");
        }
        std::vector<F> roots(n);
        fill_roots<F, CoeffCtx>(roots, inverse);
        return roots;
    }

    // Initializes exactly n + 2 elements as n ordered forward roots, 1/n, and
    // an n tag. n must be a supported transform size. This narrow native helper
    // performs Debug-only shape checks; supported C/Rust boundaries validate
    // before entry. Failure during construction may leave domain partially
    // written, so callers must not publish it until this function returns.
    template<typename F, typename CoeffCtx>
    void initialize_domain(const std::span<F> domain, const std::size_t n)
    {
        BLST_ASSERT(valid_size(n));
        BLST_ASSERT(domain.size() == domain_storage_len(n));
        fill_roots<F, CoeffCtx>(domain.first(n), false);
        domain[n] = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(n));
        domain[n + 1] = CoeffCtx::template from_u64<F>(n);
    }

    // Transforms values in place using a domain produced by initialize_domain
    // for exactly values.size(). The two spans must not overlap. Domain roots,
    // inverse size, and tag are immutable public preparation data. The helper
    // validates only shapes in Debug and does not authenticate domain contents.
    template<typename G, typename F, typename CoeffCtx>
    void transform_with_domain(
        const std::span<G> values,
        const std::span<const F> domain,
        const bool inverse)
    {
        BLST_ASSERT(valid_size(values.size()));
        BLST_ASSERT(domain.size() == domain_storage_len(values.size()));
        fft_with_forward_roots<G, F, CoeffCtx>(
            values,
            domain.first(values.size()),
            inverse,
            domain[values.size()]);
    }
}
