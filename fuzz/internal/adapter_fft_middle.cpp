// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter_common.h"
#include "bindings/blst_fft.h"
#include "bindings/blst_middle_product.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"

#include <algorithm>
#include <memory>
#include <vector>

namespace catblst::fuzz::adapter
{
    namespace
    {
        [[nodiscard]] blst_fr fr_from_u64(std::uint64_t value) noexcept
        {
            const std::array<std::uint64_t, 4> words{value, 0, 0, 0};
            blst_fr result{};
            ::blst_fr_from_uint64(&result, words.data());
            return result;
        }

        [[nodiscard]] blst_fr evaluate(
            const std::vector<blst_fr>& coefficients,
            const blst_fr& point) noexcept
        {
            blst_fr result{};
            for (std::size_t index = coefficients.size(); index-- != 0;)
            {
                ::blst_fr_mul(&result, &result, &point);
                ::blst_fr_add(&result, &result, &coefficients[index]);
            }
            return result;
        }
    }

    void run_fft_middle_case(std::span<const std::byte> input)
    {
        using detail::object_equal;
        using detail::require;
        const unsigned log_n = input.empty()
            ? 0
            : std::to_integer<unsigned char>(input.front()) % 7U;
        const std::size_t n = std::size_t{1} << log_n;
        const blst_fr seed = detail::native<blst_fr>(blst_scalar_from_le_bytes(input));
        std::vector<blst_fr> coefficients(n);
        for (std::size_t index = 0; index < n; ++index)
        {
            const blst_fr offset = fr_from_u64(index + 1);
            ::blst_fr_add(&coefficients[index], &seed, &offset);
        }

        auto allocating = coefficients;
        require(::blst_fr_fft(allocating.data(), n, false));
        auto roundtrip = allocating;
        require(::blst_fr_fft(roundtrip.data(), n, true));
        require(roundtrip == coefficients); // FUZZ_PATH: fft.roundtrip

        const std::size_t domain_len = ::blst_fr_fft_domain_len(n);
        require(domain_len != 0);
        std::vector<blst_fr> domain(domain_len);
        require(::blst_fr_fft_domain_init(domain.data(), domain.size(), n));
        auto prepared = coefficients;
        require(::blst_fr_fft_with_domain(
            prepared.data(), n, domain.data(), domain.size(), false));
        require(prepared == allocating);
        require(::blst_fr_fft_with_domain(
            prepared.data(), n, domain.data(), domain.size(), true));
        require(prepared == coefficients);

        const std::size_t xn = input.size() < 4
            ? n
            : 1 + std::to_integer<unsigned char>(input[3]) % n;
        auto truncated_input = coefficients;
        std::fill(truncated_input.begin() + xn, truncated_input.end(), blst_fr{});
        auto truncated_reference = truncated_input;
        require(::blst_fr_fft(truncated_reference.data(), n, false));
        auto truncated_native = truncated_input;
        // FUZZ_COVER: make_fft_truncated_forward_plan
        const auto forward_plan = ::make_fft_truncated_forward_plan(n, xn);
        // FUZZ_COVER: fft_forward_truncated_prefix
        ::fft_forward_truncated_prefix<blst_fr, blst_fr, blst::FrCoeffCtx>(
            std::span<blst_fr>(truncated_native),
            std::span<const blst_fr>(domain.data(), n),
            xn,
            forward_plan);
        require(truncated_native == truncated_reference);

        std::vector<blst_fr> inverse_roots(n);
        for (std::size_t index = 0; index < n; ++index)
            inverse_roots[index] = blst::FrCoeffCtx::inverse(domain[index]);
        const std::size_t inverse_start = input.size() < 5
            ? 0
            : std::to_integer<unsigned char>(input[4]) % n;
        const std::size_t inverse_length = input.size() < 6
            ? n - inverse_start
            : 1 + std::to_integer<unsigned char>(input[5]) % (n - inverse_start);
        auto inverse_native = allocating;
        // FUZZ_COVER: make_fft_inverse_truncated_range_plan
        const auto inverse_plan = ::make_fft_inverse_truncated_range_plan(
            n, inverse_start, inverse_length);
        // FUZZ_COVER: fft_inverse_truncated_range
        ::fft_inverse_truncated_range<blst_fr, blst_fr, blst::FrCoeffCtx>(
            std::span<blst_fr>(inverse_native),
            std::span<const blst_fr>(inverse_roots),
            inverse_start,
            inverse_length,
            inverse_plan);
        for (std::size_t index = inverse_start;
             index < inverse_start + inverse_length;
             ++index)
            require(object_equal(inverse_native[index], coefficients[index]));

        std::array<blst_fr, 3> invalid{seed, seed, seed};
        const auto invalid_before = invalid;
        require(!::blst_fr_fft(invalid.data(), invalid.size(), false));
        require(invalid == invalid_before);
        require(::blst_fr_fft_domain_len(3) == 0);

        const std::size_t interpolation_n = std::min<std::size_t>(n, 32);
        std::vector<blst_fr> points(interpolation_n);
        std::vector<blst_fr> values(interpolation_n);
        std::vector<blst_fr> expected_coefficients(
            coefficients.begin(), coefficients.begin() + interpolation_n);
        for (std::size_t index = 0; index < interpolation_n; ++index)
        {
            points[index] = fr_from_u64(index + 1);
            values[index] = evaluate(expected_coefficients, points[index]);
        }
        std::vector<blst_fr> interpolated(interpolation_n);
        require(::blst_fr_interpolate(
            interpolated.data(), points.data(), values.data(), interpolation_n));
        require(interpolated == expected_coefficients);

        std::vector<std::size_t> indices(n);
        for (std::size_t index = 0; index < n; ++index)
            indices[index] = index;
        const std::size_t ell = input.size() < 2
            ? n / 2
            : std::to_integer<unsigned char>(input[1]) % (n + 1);
        const std::size_t scratch_len =
            ::blst_fr_interpolate_fft_domain_suffix_scratch_len(n, n);
        require(scratch_len != 0);
        std::vector<blst_fr> scratch(scratch_len);
        std::vector<blst_fr> suffix(std::max<std::size_t>(n - ell, 1));
        require(::blst_fr_interpolate_fft_domain_suffix_with_scratch(
            suffix.data(), indices.data(), allocating.data(), n, ell, n,
            scratch.data(), scratch.size()));
        for (std::size_t index = ell; index < n; ++index)
            require(object_equal(suffix[index - ell], coefficients[index]));

        // Middle-product final exponentiation is intentionally expensive in
        // sanitizer Debug builds.  A third input byte explicitly selects this
        // subcampaign so ordinary FFT seeds do not repeat the same work.
        if (input.size() < 3)
            return;
        const std::size_t middle_n =
            1 + std::to_integer<unsigned char>(input[2]) % 16;
        std::vector<blst_p1> a(middle_n);
        std::vector<blst_p2> b(2 * middle_n);
        const auto scalar = detail::fixed_bytes<32>(input);
        for (std::size_t index = 0; index < middle_n; ++index)
        {
            std::array<blst_byte, 32> multiplier = scalar;
            multiplier[0] = static_cast<blst_byte>(index + 1);
            ::blst_p1_mult_vartime(
                &a[index], ::blst_p1_generator(), multiplier.data(), 8);
        }
        for (std::size_t index = 0; index < b.size(); ++index)
        {
            std::array<blst_byte, 32> multiplier = scalar;
            multiplier[0] = static_cast<blst_byte>(index + 1);
            ::blst_p2_mult_vartime(
                &b[index], ::blst_p2_generator(), multiplier.data(), 8);
        }
        blst_middle_product_precomputed_b* raw =
            ::blst_middle_product_precompute_b(b.data(), middle_n);
        require(raw != nullptr);
        const std::unique_ptr<
            blst_middle_product_precomputed_b,
            decltype(&::blst_middle_product_precomputed_b_free)>
            handle(raw, &::blst_middle_product_precomputed_b_free);
        std::vector<blst_fp12> full(middle_n);
        require(::blst_middle_product(a.data(), middle_n, handle.get(), full.data()));
        for (const auto& value : full)
            require(::blst_fp12_in_group(&value));

        std::vector<blst_fp12> constant(middle_n);
        require(::blst_middle_product_update(&a[0], handle.get(), constant.data()));
        auto without_constant = a;
        without_constant[0] = {};
        std::vector<blst_fp12> middle(middle_n);
        require(::blst_middle_product(
            without_constant.data(), middle_n, handle.get(), middle.data()));
        std::vector<blst_fp12> combined(middle_n);
        require(::blst_middle_product_combine(
            middle.data(), constant.data(), middle_n, combined.data()));
        for (std::size_t index = 0; index < middle_n; ++index)
            require(::blst_gt_is_equal(&combined[index], &full[index]));

        auto aliased = middle;
        require(::blst_middle_product_combine(
            aliased.data(), constant.data(), middle_n, aliased.data()));
        for (std::size_t index = 0; index < middle_n; ++index)
            require(::blst_gt_is_equal(&aliased[index], &full[index]));

        // Exercise each explicit backend on a one-element instance so the
        // differential campaign covers every ownership mode without
        // multiplying the expensive final-exponentiation work at fuzz sizes.
        constexpr std::size_t explicit_n = 1;
        blst_middle_product_precomputed_b_scalar* scalar_precomputed =
            ::blst_middle_product_precompute_b_scalar(b.data(), explicit_n);
        require(scalar_precomputed != nullptr);
        std::array<blst_fp12, explicit_n> scalar_allocating{};
        require(::blst_middle_product_scalar(
            a.data(), explicit_n, scalar_precomputed,
            scalar_allocating.data()));
        blst_middle_product_workspace_scalar* scalar_workspace =
            ::blst_middle_product_workspace_scalar_new(scalar_precomputed);
        require(scalar_workspace != nullptr);
        std::array<blst_fp12, explicit_n> scalar_reused{};
        require(::blst_middle_product_scalar_with_workspace(
            a.data(), explicit_n, scalar_precomputed, scalar_workspace,
            scalar_reused.data()));
        require(::blst_gt_is_equal(
            &scalar_reused[0], &scalar_allocating[0]));
        std::array<blst_fp12, explicit_n> scalar_update{};
        require(::blst_middle_product_update_scalar(
            &a[0], scalar_precomputed, scalar_update.data()));
        require(::blst_gt_is_equal(
            &scalar_update[0], &scalar_allocating[0]));
        ::blst_middle_product_workspace_scalar_free(scalar_workspace);
        ::blst_middle_product_precomputed_b_scalar_free(scalar_precomputed);

        if (::blst_avx_enabled())
        {
            blst_middle_product_precomputed_b_avx* avx_precomputed =
                ::blst_middle_product_precompute_b_avx(b.data(), explicit_n);
            require(avx_precomputed != nullptr);
            std::array<blst_fp12, explicit_n> avx_allocating{};
            require(::blst_middle_product_avx(
                a.data(), explicit_n, avx_precomputed,
                avx_allocating.data()));
            require(::blst_gt_is_equal(
                &avx_allocating[0], &scalar_allocating[0]));
            blst_middle_product_workspace_avx* avx_workspace =
                ::blst_middle_product_workspace_avx_new(avx_precomputed);
            require(avx_workspace != nullptr);
            std::array<blst_fp12, explicit_n> avx_reused{};
            require(::blst_middle_product_avx_with_workspace(
                a.data(), explicit_n, avx_precomputed, avx_workspace,
                avx_reused.data()));
            require(::blst_gt_is_equal(
                &avx_reused[0], &scalar_allocating[0]));
            std::array<blst_fp12, explicit_n> avx_update{};
            require(::blst_middle_product_update_avx(
                &a[0], avx_precomputed, avx_update.data()));
            require(::blst_gt_is_equal(
                &avx_update[0], &scalar_allocating[0]));
            ::blst_middle_product_workspace_avx_free(avx_workspace);
            ::blst_middle_product_precomputed_b_avx_free(avx_precomputed);

            std::array<blst_avx_p1, explicit_n> native_a{};
            std::array<blst_avx_p2, 2 * explicit_n> native_b{};
            require(::blst_avx_p1_from_scalar(&native_a[0], &a[0]));
            require(::blst_avx_p2_from_scalar(&native_b[0], &b[0]));
            require(::blst_avx_p2_from_scalar(&native_b[1], &b[1]));
            blst_middle_product_precomputed_b_avx* native_precomputed =
                ::blst_middle_product_precompute_b_avx_native(
                    native_b.data(), explicit_n);
            require(native_precomputed != nullptr);
            blst_middle_product_workspace_avx_native* native_workspace =
                ::blst_middle_product_workspace_avx_native_new(
                    native_precomputed);
            require(native_workspace != nullptr);
            std::array<blst_avx_gt, explicit_n> native_output{};
            require(::blst_middle_product_avx_native_with_workspace(
                native_a.data(), explicit_n, native_precomputed,
                native_workspace, native_output.data()));
            blst_fp12 native_scalar_output{};
            require(::blst_fp12_from_avx(
                &native_scalar_output, &native_output[0]));
            require(::blst_gt_is_equal(
                &native_scalar_output, &scalar_allocating[0]));
            ::blst_middle_product_workspace_avx_native_free(native_workspace);
            ::blst_middle_product_precomputed_b_avx_free(native_precomputed);
        }
    }
}
