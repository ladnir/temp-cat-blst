// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter_common.h"
#include "bindings/blst_avx.h"
#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/point.h"

#include <vector>

namespace catblst::fuzz::adapter
{
    namespace
    {
        struct RandomSeed
        {
            std::array<blst_byte, 32> bytes{};
        };

        bool fill_random(
            void* context,
            blst_byte* output,
            std::size_t output_size) noexcept
        {
            if (context == nullptr || output == nullptr || output_size != 32)
                return false;
            const auto& seed = *static_cast<const RandomSeed*>(context);
            std::memcpy(output, seed.bytes.data(), seed.bytes.size());
            return true;
        }

        template<class Storage>
        [[nodiscard]] std::vector<Storage> aligned_scratch(std::size_t bytes)
        {
            return std::vector<Storage>(
                (bytes + sizeof(Storage) - 1) / sizeof(Storage));
        }

        void require_g1_equal(
            const blst_avx_p1& avx,
            const blst_p1& scalar) noexcept
        {
            blst_p1 converted{};
            detail::require(::blst_p1_from_avx(&converted, &avx));
            detail::require(::blst_p1_is_equal(&converted, &scalar));
        }

        void require_g2_equal(
            const blst_avx_p2& avx,
            const blst_p2& scalar) noexcept
        {
            blst_p2 converted{};
            detail::require(::blst_p2_from_avx(&converted, &avx));
            detail::require(::blst_p2_is_equal(&converted, &scalar));
        }

        void exercise_decode(std::span<const std::byte> input)
        {
            using detail::object_equal;
            using detail::require;
            const auto g1_wire = detail::fixed_bytes<96>(input);
            blst_avx_p1 g1_sentinel = *::blst_avx_p1_generator();
            blst_avx_p1 g1 = g1_sentinel;
            const BLST_ERROR avx_g1_status =
                ::blst_avx_p1_deserialize(&g1, g1_wire.data());
            blst_p1_affine scalar_g1_affine{};
            const BLST_ERROR scalar_g1_status =
                ::blst_p1_deserialize(&scalar_g1_affine, g1_wire.data());
            require(avx_g1_status == scalar_g1_status);
            if (avx_g1_status != BLST_SUCCESS)
                require(object_equal(g1, g1_sentinel));
            else
            {
                blst_p1 scalar_g1{};
                ::blst_p1_from_affine(&scalar_g1, &scalar_g1_affine);
                require_g1_equal(g1, scalar_g1);
            }

            g1 = g1_sentinel;
            const BLST_ERROR avx_g1_compressed =
                ::blst_avx_p1_uncompress(&g1, g1_wire.data());
            const BLST_ERROR scalar_g1_compressed =
                ::blst_p1_uncompress(&scalar_g1_affine, g1_wire.data());
            require(avx_g1_compressed == scalar_g1_compressed);
            if (avx_g1_compressed != BLST_SUCCESS)
                require(object_equal(g1, g1_sentinel));

            const auto g2_wire = detail::fixed_bytes<192>(input);
            blst_avx_p2 g2_sentinel = *::blst_avx_p2_generator();
            blst_avx_p2 g2 = g2_sentinel;
            const BLST_ERROR avx_g2_status =
                ::blst_avx_p2_deserialize(&g2, g2_wire.data());
            blst_p2_affine scalar_g2_affine{};
            const BLST_ERROR scalar_g2_status =
                ::blst_p2_deserialize(&scalar_g2_affine, g2_wire.data());
            require(avx_g2_status == scalar_g2_status);
            if (avx_g2_status != BLST_SUCCESS)
                require(object_equal(g2, g2_sentinel));
            else
            {
                blst_p2 scalar_g2{};
                ::blst_p2_from_affine(&scalar_g2, &scalar_g2_affine);
                require_g2_equal(g2, scalar_g2);
            }

            g2 = g2_sentinel;
            const BLST_ERROR avx_g2_compressed =
                ::blst_avx_p2_uncompress(&g2, g2_wire.data());
            const BLST_ERROR scalar_g2_compressed =
                ::blst_p2_uncompress(&scalar_g2_affine, g2_wire.data());
            require(avx_g2_compressed == scalar_g2_compressed);
            if (avx_g2_compressed != BLST_SUCCESS)
                require(object_equal(g2, g2_sentinel));
        }

        void exercise_points(
            std::span<const std::byte> input,
            blst_p1& scalar_g1,
            blst_p2& scalar_g2,
            blst_avx_p1& avx_g1,
            blst_avx_p2& avx_g2)
        {
            using detail::require;
            const auto scalar = detail::fixed_bytes<32>(input);
            ::blst_p1_mult_vartime(
                &scalar_g1, ::blst_p1_generator(), scalar.data(), 256);
            ::blst_p2_mult_vartime(
                &scalar_g2, ::blst_p2_generator(), scalar.data(), 256);
            require(::blst_avx_p1_from_scalar(&avx_g1, &scalar_g1));
            require(::blst_avx_p2_from_scalar(&avx_g2, &scalar_g2));
            require_g1_equal(avx_g1, scalar_g1);
            require_g2_equal(avx_g2, scalar_g2);

            blst_p1_affine scalar_g1_affine{};
            blst_p2_affine scalar_g2_affine{};
            ::blst_p1_to_affine(&scalar_g1_affine, &scalar_g1);
            ::blst_p2_to_affine(&scalar_g2_affine, &scalar_g2);
            blst_avx_p1_affine avx_g1_affine{};
            blst_avx_p2_affine avx_g2_affine{};
            require(::blst_avx_p1_affine_from_scalar(
                &avx_g1_affine, &scalar_g1_affine));
            require(::blst_avx_p2_affine_from_scalar(
                &avx_g2_affine, &scalar_g2_affine));
            blst_p1_affine g1_affine_roundtrip{};
            blst_p2_affine g2_affine_roundtrip{};
            require(::blst_p1_affine_from_avx(
                &g1_affine_roundtrip, &avx_g1_affine));
            require(::blst_p2_affine_from_avx(
                &g2_affine_roundtrip, &avx_g2_affine));
            require(::blst_p1_affine_is_equal(
                &g1_affine_roundtrip, &scalar_g1_affine));
            require(::blst_p2_affine_is_equal(
                &g2_affine_roundtrip, &scalar_g2_affine));

            blst_avx_p1 affine_g1_roundtrip{};
            blst_avx_p2 affine_g2_roundtrip{};
            require(::blst_avx_p1_from_affine(
                &affine_g1_roundtrip, &avx_g1_affine));
            require(::blst_avx_p2_from_affine(
                &affine_g2_roundtrip, &avx_g2_affine));
            require(::blst_avx_p1_is_equal(&affine_g1_roundtrip, &avx_g1));
            require(::blst_avx_p2_is_equal(&affine_g2_roundtrip, &avx_g2));
            require(::blst_avx_p1_to_affine(&avx_g1_affine, &avx_g1));
            require(::blst_avx_p2_to_affine(&avx_g2_affine, &avx_g2));

            blst_avx_p1 avx_g1_sum{};
            blst_avx_p2 avx_g2_sum{};
            require(::blst_avx_p1_add(&avx_g1_sum, &avx_g1, &avx_g1));
            require(::blst_avx_p2_add(&avx_g2_sum, &avx_g2, &avx_g2));
            blst_p1 scalar_g1_sum{};
            blst_p2 scalar_g2_sum{};
            ::blst_p1_double(&scalar_g1_sum, &scalar_g1);
            ::blst_p2_double(&scalar_g2_sum, &scalar_g2);
            require_g1_equal(avx_g1_sum, scalar_g1_sum);
            require_g2_equal(avx_g2_sum, scalar_g2_sum);

            require(::blst_avx_p1_cneg(&avx_g1_sum, true));
            require(::blst_avx_p2_cneg(&avx_g2_sum, true));
            ::blst_p1_cneg(&scalar_g1_sum, true);
            ::blst_p2_cneg(&scalar_g2_sum, true);
            require_g1_equal(avx_g1_sum, scalar_g1_sum);
            require_g2_equal(avx_g2_sum, scalar_g2_sum);

            blst_avx_p1 avx_g1_product{};
            blst_avx_p2 avx_g2_product{};
            require(::blst_avx_p1_mult(
                &avx_g1_product, &avx_g1, scalar.data(), 256));
            require(::blst_avx_p2_mult(
                &avx_g2_product, &avx_g2, scalar.data(), 256));
            blst_p1 scalar_g1_product{};
            blst_p2 scalar_g2_product{};
            ::blst_p1_mult_vartime(
                &scalar_g1_product, &scalar_g1, scalar.data(), 256);
            ::blst_p2_mult_vartime(
                &scalar_g2_product, &scalar_g2, scalar.data(), 256);
            require_g1_equal(avx_g1_product, scalar_g1_product);
            require_g2_equal(avx_g2_product, scalar_g2_product);

            require(::blst_avx_p1_in_g1(&avx_g1));
            require(::blst_avx_p2_in_g2(&avx_g2));
            require(::blst_avx_p1_is_inf(&avx_g1) ==
                    ::blst_p1_is_inf(&scalar_g1));
            require(::blst_avx_p2_is_inf(&avx_g2) ==
                    ::blst_p2_is_inf(&scalar_g2));

            std::array<blst_byte, 96> g1_serialized{};
            std::array<blst_byte, 48> g1_compressed{};
            std::array<blst_byte, 192> g2_serialized{};
            std::array<blst_byte, 96> g2_compressed{};
            require(::blst_avx_p1_serialize(g1_serialized.data(), &avx_g1));
            require(::blst_avx_p1_compress(g1_compressed.data(), &avx_g1));
            require(::blst_avx_p2_serialize(g2_serialized.data(), &avx_g2));
            require(::blst_avx_p2_compress(g2_compressed.data(), &avx_g2));
            std::array<blst_byte, 96> scalar_g1_serialized{};
            std::array<blst_byte, 48> scalar_g1_compressed{};
            std::array<blst_byte, 192> scalar_g2_serialized{};
            std::array<blst_byte, 96> scalar_g2_compressed{};
            ::blst_p1_serialize(scalar_g1_serialized.data(), &scalar_g1);
            ::blst_p1_compress(scalar_g1_compressed.data(), &scalar_g1);
            ::blst_p2_serialize(scalar_g2_serialized.data(), &scalar_g2);
            ::blst_p2_compress(scalar_g2_compressed.data(), &scalar_g2);
            require(g1_serialized == scalar_g1_serialized);
            require(g1_compressed == scalar_g1_compressed);
            require(g2_serialized == scalar_g2_serialized);
            require(g2_compressed == scalar_g2_compressed);

            RandomSeed random{detail::fixed_bytes<32>(input)};
            const blst_avx_p1* g1_points[]{&avx_g1};
            const blst_avx_p2* g2_points[]{&avx_g2};
            require(::blst_avx_p1s_in_g1(
                g1_points, 1, fill_random, &random));
            std::vector<blst_byte> g2_membership_workspace(
                ::blst_avx_p2s_in_g2_workspace_sizeof());
            require(::blst_avx_p2s_in_g2(
                g2_points, 1, fill_random, &random,
                g2_membership_workspace.data(),
                g2_membership_workspace.size()));
            const std::array<blst_byte, 1> coefficients{1};
            require(::blst_avx_p1_batch_subgroup_check(
                g1_points, 1, coefficients.data(), 1,
                ::blst_avx_p1_generator()));
            require(::blst_avx_p2_batch_subgroup_check(
                g2_points, 1, coefficients.data(), 1,
                ::blst_avx_p2_generator()));
        }

        void exercise_msm(std::span<const std::byte> input)
        {
            using detail::require;
            const std::size_t count = input.empty()
                ? 1
                : 1 + std::to_integer<unsigned char>(input.front()) % 4;
            std::vector<blst_p1> scalar_g1(count);
            std::vector<blst_p2> scalar_g2(count);
            std::vector<blst_p1_affine> scalar_g1_affine(count);
            std::vector<blst_p2_affine> scalar_g2_affine(count);
            std::vector<blst_avx_p1> avx_g1(count);
            std::vector<blst_avx_p2> avx_g2(count);
            std::vector<blst_avx_p1_affine> avx_g1_affine(count);
            std::vector<blst_avx_p2_affine> avx_g2_affine(count);
            std::vector<blst_byte> scalars(count * 32);
            blst_p1 expected_g1{};
            blst_p2 expected_g2{};
            for (std::size_t index = 0; index < count; ++index)
            {
                scalars[index * 32] = static_cast<blst_byte>(index + 1);
                const std::array<blst_byte, 1> point_scalar{
                    static_cast<blst_byte>(index + 2)};
                ::blst_p1_mult_vartime(
                    &scalar_g1[index], ::blst_p1_generator(),
                    point_scalar.data(), 8);
                ::blst_p2_mult_vartime(
                    &scalar_g2[index], ::blst_p2_generator(),
                    point_scalar.data(), 8);
                ::blst_p1_to_affine(
                    &scalar_g1_affine[index], &scalar_g1[index]);
                ::blst_p2_to_affine(
                    &scalar_g2_affine[index], &scalar_g2[index]);
                require(::blst_avx_p1_from_scalar(
                    &avx_g1[index], &scalar_g1[index]));
                require(::blst_avx_p2_from_scalar(
                    &avx_g2[index], &scalar_g2[index]));
                require(::blst_avx_p1_affine_from_scalar(
                    &avx_g1_affine[index], &scalar_g1_affine[index]));
                require(::blst_avx_p2_affine_from_scalar(
                    &avx_g2_affine[index], &scalar_g2_affine[index]));
                blst_p1 g1_term{};
                blst_p2 g2_term{};
                ::blst_p1_mult_vartime(
                    &g1_term, &scalar_g1[index], scalars.data() + index * 32, 256);
                ::blst_p2_mult_vartime(
                    &g2_term, &scalar_g2[index], scalars.data() + index * 32, 256);
                ::blst_p1_add_or_double(&expected_g1, &expected_g1, &g1_term);
                ::blst_p2_add_or_double(&expected_g2, &expected_g2, &g2_term);
            }

            auto g1_scratch = aligned_scratch<blst_avx_p1>(
                ::blst_avx_p1_mult_pippenger_scratch_sizeof(count));
            std::vector<blst_byte> g1_flags(
                ::blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof(count));
            blst_avx_p1 g1_result{};
            require(::blst_avx_p1_mult_pippenger(
                &g1_result, avx_g1_affine.data(), scalars.data(), count, 256,
                g1_scratch.data(),
                g1_scratch.size() * sizeof(blst_avx_p1),
                g1_flags.data(), g1_flags.size()));
            require_g1_equal(g1_result, expected_g1);

            auto g2_scratch = aligned_scratch<blst_avx_p2>(
                ::blst_avx_p2_mult_pippenger_scratch_sizeof(count));
            std::vector<blst_byte> g2_flags(
                ::blst_avx_p2_mult_pippenger_bucket_used_scratch_sizeof(count));
            blst_avx_p2 g2_result{};
            require(::blst_avx_p2_mult_pippenger(
                &g2_result, avx_g2_affine.data(), scalars.data(), count, 256,
                g2_scratch.data(),
                g2_scratch.size() * sizeof(blst_avx_p2),
                g2_flags.data(), g2_flags.size()));
            require_g2_equal(g2_result, expected_g2);

            auto g1_projective_scratch = aligned_scratch<blst_limb_t>(
                ::blst_avx_p1s_mult_pippenger_scratch_sizeof(count));
            auto g2_projective_scratch = aligned_scratch<blst_limb_t>(
                ::blst_avx_p2s_mult_pippenger_scratch_sizeof(count));
            require(::blst_avx_p1s_mult_pippenger(
                &g1_result, avx_g1.data(), count, scalars.data(), 256,
                g1_projective_scratch.empty()
                    ? nullptr : g1_projective_scratch.data()));
            require(::blst_avx_p2s_mult_pippenger(
                &g2_result, avx_g2.data(), count, scalars.data(), 256,
                g2_projective_scratch.empty()
                    ? nullptr : g2_projective_scratch.data()));
            require_g1_equal(g1_result, expected_g1);
            require_g2_equal(g2_result, expected_g2);

            std::vector<blst_fr> secret_scalars(count);
            for (std::size_t index = 0; index < count; ++index)
            {
                blst_scalar reduced{};
                ::blst_scalar_from_lendian(
                    &reduced, scalars.data() + index * 32);
                ::blst_fr_from_scalar(&secret_scalars[index], &reduced);
            }
            auto secret_scratch = aligned_scratch<blst_avx_p1>(
                ::blst_avx_p1s_mult_ct_scratch_sizeof());
            require(::blst_avx_p1s_mult_ct(
                &g1_result, avx_g1.data(), secret_scalars.data(), count,
                secret_scratch.data(),
                secret_scratch.size() * sizeof(blst_avx_p1)));
            require_g1_equal(g1_result, expected_g1);

            auto public_scratch = aligned_scratch<blst_avx_p1>(
                ::blst_avx_p1s_mult_fr_vartime_scratch_sizeof(count));
            require(::blst_avx_p1s_mult_fr_vartime(
                &g1_result, avx_g1.data(), secret_scalars.data(), count,
                public_scratch.empty() ? nullptr : public_scratch.data(),
                public_scratch.size() * sizeof(blst_avx_p1)));
            require_g1_equal(g1_result, expected_g1);

            auto split_points = aligned_scratch<blst_avx_p1_affine>(
                ::blst_avx_p1_mult_pippenger_glv_points_scratch_sizeof(count));
            std::vector<blst_byte> split_scalars(
                ::blst_avx_p1_mult_pippenger_glv_scalars_scratch_sizeof(count));
            require(::blst_avx_p1_mult_pippenger_glv_with_scratch(
                &g1_result, avx_g1_affine.data(), scalars.data(), count,
                g1_scratch.data(), g1_scratch.size() * sizeof(blst_avx_p1),
                g1_flags.data(), g1_flags.size(),
                split_points.data(),
                split_points.size() * sizeof(blst_avx_p1_affine),
                split_scalars.data(), split_scalars.size()));
            require_g1_equal(g1_result, expected_g1);

            constexpr std::size_t window = 4;
            auto window_scratch = aligned_scratch<blst_avx_p1>(
                ::blst_avx_p1_mult_pippenger_glv_window_scratch_sizeof(
                    count, window));
            std::vector<blst_byte> window_flags(
                ::blst_avx_p1_mult_pippenger_glv_window_bucket_used_scratch_sizeof(
                    count, window));
            require(::blst_avx_p1_mult_pippenger_glv_window_with_scratch(
                &g1_result, avx_g1_affine.data(), scalars.data(), count, window,
                window_scratch.data(),
                window_scratch.size() * sizeof(blst_avx_p1),
                window_flags.data(), window_flags.size(),
                split_points.data(),
                split_points.size() * sizeof(blst_avx_p1_affine),
                split_scalars.data(), split_scalars.size()));
            require_g1_equal(g1_result, expected_g1);
        }

        void exercise_fp12_and_pairing(
            std::span<const std::byte> input,
            const blst_p1& scalar_g1,
            const blst_p2& scalar_g2,
            const blst_avx_p1& avx_g1,
            const blst_avx_p2& avx_g2)
        {
            using detail::require;
            const blst_fp12 scalar_a =
                detail::native<blst_fp12>(reduce_fp12(input));
            const blst_fp12 scalar_b =
                detail::native<blst_fp12>(reduce_fp12(input.subspan(
                    std::min<std::size_t>(input.size(), 1))));
            blst_avx_fp12 avx_a{};
            blst_avx_fp12 avx_b{};
            require(::blst_avx_fp12_from_scalar(&avx_a, &scalar_a));
            require(::blst_avx_fp12_from_scalar(&avx_b, &scalar_b));
            blst_fp12 roundtrip{};
            require(::blst_fp12_from_avx(&roundtrip, &avx_a));
            require(::blst_fp12_is_equal(&roundtrip, &scalar_a));

            blst_avx_fp12 avx_product{};
            blst_fp12 scalar_product{};
            require(::blst_avx_fp12_mul(&avx_product, &avx_a, &avx_b));
            ::blst_fp12_mul(&scalar_product, &scalar_a, &scalar_b);
            require(::blst_fp12_from_avx(&roundtrip, &avx_product));
            require(::blst_fp12_is_equal(&roundtrip, &scalar_product));
            require(::blst_avx_fp12_sqr(&avx_product, &avx_a));
            ::blst_fp12_sqr(&scalar_product, &scalar_a);
            require(::blst_fp12_from_avx(&roundtrip, &avx_product));
            require(::blst_fp12_is_equal(&roundtrip, &scalar_product));

            blst_p1_affine g1_affine{};
            blst_p2_affine g2_affine{};
            ::blst_p1_to_affine(&g1_affine, &scalar_g1);
            ::blst_p2_to_affine(&g2_affine, &scalar_g2);
            blst_fp12 miller{};
            ::blst_miller_loop(&miller, &g2_affine, &g1_affine);
            blst_fp12 scalar_gt{};
            ::blst_final_exp(&scalar_gt, &miller);
            blst_avx_fp12 avx_gt{};
            require(::blst_avx_pairing(&avx_gt, &avx_g1, &avx_g2));
            require(::blst_fp12_from_avx(&roundtrip, &avx_gt));
            require(::blst_fp12_is_equal(&roundtrip, &scalar_gt));
            require(::blst_avx_fp12_in_group(&avx_gt));
            require(::blst_avx_fp12_is_equal(&avx_gt, &avx_gt));

            blst_avx_fp12 avx_miller{};
            require(::blst_avx_fp12_from_scalar(&avx_miller, &miller));
            require(::blst_avx_fp12_final_exp(&avx_product, &avx_miller));
            require(::blst_avx_fp12_is_equal(&avx_product, &avx_gt));
            blst_fp12 scalar_one = *::blst_fp12_one();
            require(::blst_fp12_from_avx(
                &roundtrip, ::blst_avx_fp12_one()));
            require(::blst_fp12_is_equal(&roundtrip, &scalar_one));

            blst_avx_gt avx_gt_inverse{};
            require(::blst_avx_gt_inverse(&avx_gt_inverse, &avx_gt));
            require(::blst_avx_fp12_mul(
                &avx_product, &avx_gt, &avx_gt_inverse));
            require(::blst_avx_fp12_is_equal(
                &avx_product, ::blst_avx_fp12_one()));

            std::array<blst_byte, 576> encoded_gt{};
            require(::blst_avx_gt_serialize(encoded_gt.data(), &avx_gt));
            blst_avx_gt decoded_gt{};
            require(::blst_avx_gt_deserialize(&decoded_gt, encoded_gt.data()));
            require(::blst_avx_fp12_is_equal(&decoded_gt, &avx_gt));

            const auto scalar = detail::fixed_bytes<32>(input);
            const std::array<blst_avx_gt, 1> points{avx_gt};
            auto gt_scratch = aligned_scratch<blst_avx_fp12>(
                ::blst_avx_gt_mult_pippenger_scratch_sizeof(1));
            std::vector<blst_byte> gt_flags(
                ::blst_avx_gt_mult_pippenger_bucket_used_scratch_sizeof(1));
            blst_avx_gt gt_msm{};
            blst_fp12 scalar_gt_product{};
            ::blst_gt_mult_vartime(
                &scalar_gt_product, &scalar_gt, scalar.data(), 256);
            require(::blst_avx_gt_mult(
                &gt_msm, &avx_gt, scalar.data(), 256));
            require(::blst_fp12_from_avx(&roundtrip, &gt_msm));
            require(::blst_fp12_is_equal(&roundtrip, &scalar_gt_product));

            require(::blst_avx_gt_mult_pippenger(
                &gt_msm, points.data(), scalar.data(), 1,
                gt_scratch.data(),
                gt_scratch.size() * sizeof(blst_avx_fp12),
                gt_flags.data(), gt_flags.size()));
            require(::blst_fp12_from_avx(&roundtrip, &gt_msm));
            require(::blst_fp12_is_equal(&roundtrip, &scalar_gt_product));

            blst_avx_p1 g1_identity{};
            blst_avx_p2 g2_identity{};
            require(::blst_avx_p1_identity(&g1_identity));
            require(::blst_avx_p2_identity(&g2_identity));
            require(::blst_avx_p1_is_inf(&g1_identity));
            require(::blst_avx_p2_is_inf(&g2_identity));

            blst_avx_p1_affine avx_g1_affine{};
            blst_avx_p2_affine avx_g2_affine{};
            require(::blst_avx_p1_to_affine(&avx_g1_affine, &avx_g1));
            require(::blst_avx_p2_to_affine(&avx_g2_affine, &avx_g2));
            blst_avx_fp12 batched_pairing{};
            require(::blst_avx_pairing_n(
                &batched_pairing, &avx_g1_affine, &avx_g2_affine, 1));
            require(::blst_avx_fp12_is_equal(&batched_pairing, &avx_gt));

            auto prepared_lines = aligned_scratch<blst_avx_fp12>(
                ::blst_avx_prepared_g2_lines_sizeof());
            require(
                alignof(blst_avx_fp12) >=
                ::blst_avx_prepared_g2_lines_alignof());
            require(::blst_avx_prepared_g2_lines_init(
                prepared_lines.data(),
                prepared_lines.size() * sizeof(blst_avx_fp12), &avx_g2));
            blst_avx_fp12 prepared_pairing{};
            require(::blst_avx_pairing_prepared_g2(
                &prepared_pairing, prepared_lines.data(),
                prepared_lines.size() * sizeof(blst_avx_fp12), &avx_g1));
            require(::blst_avx_fp12_is_equal(&prepared_pairing, &avx_gt));

            blst_avx_prepared_g2_batch* raw_batch =
                ::blst_avx_prepared_g2_batch_new(&avx_g2, 1);
            require(raw_batch != nullptr);
            require(::blst_avx_prepared_g2_batch_len(raw_batch) == 1);
            require(
                ::blst_avx_prepared_g2_batch_storage_bytes(raw_batch)
                != 0);
            blst_avx_prepared_g2_batch_workspace* raw_workspace =
                ::blst_avx_prepared_g2_batch_workspace_new(
                    raw_batch);
            require(raw_workspace != nullptr);
            std::array<blst_avx_gt, 1> prepared_outputs{};
            require(::blst_avx_prepared_g2_batch_pairing_prefix(
                prepared_outputs.data(), 1, &avx_g1, raw_batch,
                raw_workspace));
            require(::blst_avx_fp12_is_equal(
                &prepared_outputs[0], &avx_gt));
            constexpr std::array<std::size_t, 1> indices{0};
            require(::blst_avx_prepared_g2_batch_pairing_indices(
                prepared_outputs.data(), indices.data(), indices.size(),
                &avx_g1, raw_batch, raw_workspace));
            require(::blst_avx_fp12_is_equal(
                &prepared_outputs[0], &avx_gt));
            blst_avx_fp12 prepared_miller{};
            require(
                ::blst_avx_prepared_g2_batch_miller_loop_product_indices(
                    &prepared_miller, &avx_g1, indices.data(),
                    indices.size(), raw_batch));
            blst_avx_fp12 prepared_product{};
            require(::blst_avx_fp12_final_exp(
                &prepared_product, &prepared_miller));
            require(::blst_avx_fp12_is_equal(
                &prepared_product, &avx_gt));
            constexpr std::array<std::size_t, 1> bad_indices{1};
            blst_avx_fp12 unchanged = avx_gt;
            require(
                !::blst_avx_prepared_g2_batch_miller_loop_product_indices(
                    &unchanged, &avx_g1, bad_indices.data(),
                    bad_indices.size(), raw_batch));
            require(::blst_avx_fp12_is_equal(&unchanged, &avx_gt));
            ::blst_avx_prepared_g2_batch_workspace_free(raw_workspace);
            ::blst_avx_prepared_g2_batch_free(raw_batch);
        }

        void exercise_unique_native_avx_methods()
        {
            using detail::require;
            using blst::avx::G1Point;
            using blst::avx::G2Point;
            using blst::avx::GT;

            G1Point g1_sum{};
            // FUZZ_COVER: blst::avx::G1Point::dadd
            G1Point::dadd(g1_sum, G1Point::generator(), G1Point::generator());
            G1Point g1_double{};
            G1Point::doubling(g1_double, G1Point::generator());
            require(G1Point::equal(g1_sum, g1_double));

            G2Point g2_sum{};
            // FUZZ_COVER: blst::avx::G2Point::dadd
            G2Point::dadd(g2_sum, G2Point::generator(), G2Point::generator());
            G2Point g2_double{};
            G2Point::doubling(g2_double, G2Point::generator());
            require(G2Point::equal(g2_sum, g2_double));

            GT gt_sum{};
            // FUZZ_COVER: blst::avx::GT::dadd
            GT::dadd(gt_sum, GT::identity(), GT::generator());
            require(GT::equal(gt_sum, GT::generator()));

            std::array<unsigned char, 12 * 48> encoded{};
            GT::serialize(encoded.data(), GT::generator());
            GT decoded = GT::identity();
            // FUZZ_COVER: blst::avx::GT::deserialize_unchecked
            require(GT::deserialize_unchecked(decoded, encoded.data()));
            require(GT::equal(decoded, GT::generator()));
        }
    }

    void run_avx_case(std::span<const std::byte> input)
    {
        using detail::require;
        require(::blst_avx_enabled());
        require(::blst_avx_p1_sizeof() == sizeof(blst_avx_p1));
        require(::blst_avx_p1_alignof() == alignof(blst_avx_p1));
        require(::blst_avx_p1_affine_sizeof() == sizeof(blst_avx_p1_affine));
        require(::blst_avx_p1_affine_alignof() == alignof(blst_avx_p1_affine));
        require(::blst_avx_p2_sizeof() == sizeof(blst_avx_p2));
        require(::blst_avx_p2_alignof() == alignof(blst_avx_p2));
        require(::blst_avx_p2_affine_sizeof() == sizeof(blst_avx_p2_affine));
        require(::blst_avx_p2_affine_alignof() == alignof(blst_avx_p2_affine));
        require(::blst_avx_fp12_sizeof() == sizeof(blst_avx_fp12));
        require(::blst_avx_fp12_alignof() == alignof(blst_avx_fp12));

        exercise_decode(input);
        blst_p1 scalar_g1{};
        blst_p2 scalar_g2{};
        blst_avx_p1 avx_g1{};
        blst_avx_p2 avx_g2{};
        exercise_points(input, scalar_g1, scalar_g2, avx_g1, avx_g2);
        exercise_msm(input);
        exercise_fp12_and_pairing(
            input, scalar_g1, scalar_g2, avx_g1, avx_g2);
        exercise_unique_native_avx_methods();
    }
}
