// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter_common.h"

#include <vector>

namespace catblst::fuzz::adapter
{
    namespace
    {
        [[nodiscard]] std::vector<blst_limb_t> scratch(std::size_t bytes)
        {
            return std::vector<blst_limb_t>(
                (bytes + sizeof(blst_limb_t) - 1) / sizeof(blst_limb_t));
        }

        void exercise_g1_msm(
            std::span<const std::byte> input,
            std::size_t n)
        {
            using detail::require;
            std::vector<blst_p1_affine> points(n);
            std::vector<std::array<blst_byte, 1>> scalar_storage(n);
            std::vector<const blst_byte*> scalar_heads(n);
            blst_p1 expected{};
            for (std::size_t index = 0; index < n; ++index)
            {
                const std::array<blst_byte, 1> point_scalar{
                    static_cast<blst_byte>(index + 1)};
                blst_p1 point{};
                ::blst_p1_mult_vartime(
                    &point, ::blst_p1_generator(), point_scalar.data(), 8);
                ::blst_p1_to_affine(&points[index], &point);
                scalar_storage[index][0] = input.empty()
                    ? static_cast<blst_byte>(index + 1)
                    : static_cast<blst_byte>(
                          std::to_integer<unsigned char>(input[index % input.size()]) +
                          static_cast<unsigned char>(index + 1));
                scalar_heads[index] = scalar_storage[index].data();
                blst_p1 term{};
                ::blst_p1_mult_vartime(
                    &term, &point, scalar_heads[index], 8);
                blst_p1 sum{};
                ::blst_p1_add_or_double(&sum, &expected, &term);
                expected = sum;
            }
            std::vector<const blst_p1_affine*> point_heads(n, nullptr);
            point_heads[0] = points.data();

            const std::size_t auto_bytes =
                ::blst_p1s_mult_pippenger_scratch_sizeof(n);
            auto auto_scratch = scratch(auto_bytes);
            blst_p1 result{};
            ::blst_p1s_mult_pippenger(
                &result, point_heads.data(), n, scalar_heads.data(), 8,
                auto_scratch.data());
            require(::blst_p1_is_equal(&result, &expected));

            constexpr std::size_t window = 8;
            auto window_scratch = scratch(
                ::blst_p1s_mult_pippenger_window_scratch_sizeof(window));
            ::blst_p1s_mult_pippenger_window(
                &result, point_heads.data(), n, scalar_heads.data(), 8,
                window_scratch.data(), window);
            require(::blst_p1_is_equal(&result, &expected));

            // A tile is one signed Booth digit, not a complete MSM.  Keep the
            // top overlap bit clear so the bit-zero, eight-bit digit equals
            // the direct unsigned value used by this focused oracle.  The
            // complete MSM calls above continue to exercise all eight bits.
            std::vector<std::array<blst_byte, 1>> tile_scalar_storage =
                scalar_storage;
            std::vector<const blst_byte*> tile_scalar_heads(n);
            blst_p1 tile_expected{};
            for (std::size_t index = 0; index < n; ++index)
            {
                tile_scalar_storage[index][0] &= 0x7f;
                tile_scalar_heads[index] = tile_scalar_storage[index].data();
                blst_p1 point{};
                ::blst_p1_from_affine(&point, &points[index]);
                blst_p1 term{};
                ::blst_p1_mult_vartime(
                    &term, &point, tile_scalar_heads[index], 8);
                blst_p1 sum{};
                ::blst_p1_add_or_double(&sum, &tile_expected, &term);
                tile_expected = sum;
            }
            ::blst_p1s_tile_pippenger(
                &result, point_heads.data(), n, tile_scalar_heads.data(), 8,
                window_scratch.data(), 0, window);
            require(::blst_p1_is_equal(&result, &tile_expected));

            constexpr std::size_t wbits = 4;
            std::vector<blst_p1_affine> table(
                ::blst_p1s_mult_wbits_precompute_sizeof(wbits, n) /
                sizeof(blst_p1_affine));
            ::blst_p1s_mult_wbits_precompute(
                table.data(), wbits, point_heads.data(), n);
            auto wbits_scratch = scratch(
                ::blst_p1s_mult_wbits_scratch_sizeof(n));
            ::blst_p1s_mult_wbits(
                &result, table.data(), wbits, n, scalar_heads.data(), 8,
                wbits_scratch.data());
            require(::blst_p1_is_equal(&result, &expected));

            std::vector<blst_fr> fr_scalars(n);
            for (std::size_t index = 0; index < n; ++index)
            {
                blst_scalar scalar{};
                ::blst_scalar_from_le_bytes(
                    &scalar, scalar_storage[index].data(), 1);
                ::blst_fr_from_scalar(&fr_scalars[index], &scalar);
            }
            const std::size_t fr_bytes =
                ::blst_p1s_mult_pippenger_fr_scratch_sizeof(n);
            auto fr_scratch = scratch(fr_bytes);
            require(::blst_p1s_mult_pippenger_fr(
                &result, points.data(), fr_scalars.data(), n,
                fr_scratch.data(), fr_bytes));
            require(::blst_p1_is_equal(&result, &expected));
            const std::size_t selected_bytes = ::blst_p1s_mult_fr_scratch_sizeof(n);
            auto selected_scratch = scratch(selected_bytes);
            require(::blst_p1s_mult_fr(
                &result, points.data(), fr_scalars.data(), n,
                selected_scratch.data(), selected_bytes));
            require(::blst_p1_is_equal(&result, &expected));

            const std::size_t batch_bytes =
                ::blst_p1s_mult_batch_affine_fr_scratch_sizeof(n, 4, 128);
            auto batch_scratch = scratch(batch_bytes);
            require(::blst_p1s_mult_batch_affine_fr(
                &result, points.data(), fr_scalars.data(), n,
                batch_scratch.data(), batch_bytes, 4, 128));
            require(::blst_p1_is_equal(&result, &expected));
        }

        void exercise_g2_msm(
            std::span<const std::byte> input,
            std::size_t n)
        {
            using detail::require;
            std::vector<blst_p2_affine> points(n);
            std::vector<std::array<blst_byte, 1>> scalar_storage(n);
            std::vector<const blst_byte*> scalar_heads(n);
            blst_p2 expected{};
            for (std::size_t index = 0; index < n; ++index)
            {
                const std::array<blst_byte, 1> point_scalar{
                    static_cast<blst_byte>(index + 1)};
                blst_p2 point{};
                ::blst_p2_mult_vartime(
                    &point, ::blst_p2_generator(), point_scalar.data(), 8);
                ::blst_p2_to_affine(&points[index], &point);
                scalar_storage[index][0] = input.empty()
                    ? static_cast<blst_byte>(index + 1)
                    : static_cast<blst_byte>(
                          std::to_integer<unsigned char>(input[index % input.size()]) +
                          static_cast<unsigned char>(index + 1));
                scalar_heads[index] = scalar_storage[index].data();
                blst_p2 term{};
                ::blst_p2_mult_vartime(
                    &term, &point, scalar_heads[index], 8);
                blst_p2 sum{};
                ::blst_p2_add_or_double(&sum, &expected, &term);
                expected = sum;
            }
            std::vector<const blst_p2_affine*> point_heads(n, nullptr);
            point_heads[0] = points.data();
            const std::size_t auto_bytes =
                ::blst_p2s_mult_pippenger_scratch_sizeof(n);
            auto auto_scratch = scratch(auto_bytes);
            blst_p2 result{};
            ::blst_p2s_mult_pippenger(
                &result, point_heads.data(), n, scalar_heads.data(), 8,
                auto_scratch.data());
            require(::blst_p2_is_equal(&result, &expected));
            constexpr std::size_t window = 8;
            auto window_scratch = scratch(
                ::blst_p2s_mult_pippenger_window_scratch_sizeof(window));
            ::blst_p2s_mult_pippenger_window(
                &result, point_heads.data(), n, scalar_heads.data(), 8,
                window_scratch.data(), window);
            require(::blst_p2_is_equal(&result, &expected));

            std::vector<std::array<blst_byte, 1>> tile_scalar_storage =
                scalar_storage;
            std::vector<const blst_byte*> tile_scalar_heads(n);
            blst_p2 tile_expected{};
            for (std::size_t index = 0; index < n; ++index)
            {
                tile_scalar_storage[index][0] &= 0x7f;
                tile_scalar_heads[index] = tile_scalar_storage[index].data();
                blst_p2 point{};
                ::blst_p2_from_affine(&point, &points[index]);
                blst_p2 term{};
                ::blst_p2_mult_vartime(
                    &term, &point, tile_scalar_heads[index], 8);
                blst_p2 sum{};
                ::blst_p2_add_or_double(&sum, &tile_expected, &term);
                tile_expected = sum;
            }
            ::blst_p2s_tile_pippenger(
                &result, point_heads.data(), n, tile_scalar_heads.data(), 8,
                window_scratch.data(), 0, window);
            require(::blst_p2_is_equal(&result, &tile_expected));

            constexpr std::size_t wbits = 4;
            std::vector<blst_p2_affine> table(
                ::blst_p2s_mult_wbits_precompute_sizeof(wbits, n) /
                sizeof(blst_p2_affine));
            ::blst_p2s_mult_wbits_precompute(
                table.data(), wbits, point_heads.data(), n);
            auto wbits_scratch = scratch(::blst_p2s_mult_wbits_scratch_sizeof(n));
            ::blst_p2s_mult_wbits(
                &result, table.data(), wbits, n, scalar_heads.data(), 8,
                wbits_scratch.data());
            require(::blst_p2_is_equal(&result, &expected));

            std::vector<blst_fr> fr_scalars(n);
            for (std::size_t index = 0; index < n; ++index)
            {
                blst_scalar scalar{};
                ::blst_scalar_from_le_bytes(
                    &scalar, scalar_storage[index].data(), 1);
                ::blst_fr_from_scalar(&fr_scalars[index], &scalar);
            }
            const std::size_t fr_bytes =
                ::blst_p2s_mult_pippenger_fr_scratch_sizeof(n);
            auto fr_scratch = scratch(fr_bytes);
            require(::blst_p2s_mult_pippenger_fr(
                &result, points.data(), fr_scalars.data(), n,
                fr_scratch.data(), fr_bytes));
            require(::blst_p2_is_equal(&result, &expected));
            const std::size_t selected_bytes = ::blst_p2s_mult_fr_scratch_sizeof(n);
            auto selected_scratch = scratch(selected_bytes);
            require(::blst_p2s_mult_fr(
                &result, points.data(), fr_scalars.data(), n,
                selected_scratch.data(), selected_bytes));
            require(::blst_p2_is_equal(&result, &expected));
            const std::size_t batch_bytes =
                ::blst_p2s_mult_batch_affine_fr_scratch_sizeof(n, 4, 128);
            auto batch_scratch = scratch(batch_bytes);
            require(::blst_p2s_mult_batch_affine_fr(
                &result, points.data(), fr_scalars.data(), n,
                batch_scratch.data(), batch_bytes, 4, 128));
            require(::blst_p2_is_equal(&result, &expected));
        }
    }

    void run_msm_case(std::span<const std::byte> input)
    {
        const std::size_t n = input.empty()
            ? 1
            : (std::to_integer<unsigned char>(input.front()) % 16U) + 1;
        exercise_g1_msm(input, n);
        exercise_g2_msm(input, n);
    }
}
