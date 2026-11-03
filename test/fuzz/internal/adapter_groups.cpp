// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "test/fuzz/internal/adapter_common.h"

#include <algorithm>

namespace catblst::fuzz::adapter
{
    namespace
    {
        [[nodiscard]] std::size_t scalar_bits(
            std::span<const std::byte> input) noexcept
        {
            if (input.empty())
                return 0;
            return std::to_integer<unsigned char>(input.front()) % 257U;
        }

        void direct_mult(
            blst_p1& out,
            const blst_p1& point,
            const std::array<blst_byte, 32>& scalar,
            std::size_t nbits) noexcept
        {
            out = {};
            for (std::size_t bit = nbits; bit-- != 0;)
            {
                blst_p1 doubled{};
                ::blst_p1_double(&doubled, &out);
                out = doubled;
                if (((scalar[bit / 8] >> (bit & 7U)) & 1U) != 0)
                {
                    blst_p1 sum{};
                    ::blst_p1_add_or_double(&sum, &out, &point);
                    out = sum;
                }
            }
        }

        void direct_mult(
            blst_p2& out,
            const blst_p2& point,
            const std::array<blst_byte, 32>& scalar,
            std::size_t nbits) noexcept
        {
            out = {};
            for (std::size_t bit = nbits; bit-- != 0;)
            {
                blst_p2 doubled{};
                ::blst_p2_double(&doubled, &out);
                out = doubled;
                if (((scalar[bit / 8] >> (bit & 7U)) & 1U) != 0)
                {
                    blst_p2 sum{};
                    ::blst_p2_add_or_double(&sum, &out, &point);
                    out = sum;
                }
            }
        }

        void exercise_g1(
            std::span<const std::byte> input,
            const std::array<blst_byte, 32>& scalar,
            std::size_t nbits)
        {
            using detail::require;
            const bool use_wide = !input.empty() &&
                (std::to_integer<unsigned char>(input.back()) & 1U) != 0;
            const OnCurveG1AffineValue wide_value = use_wide
                ? retained_non_subgroup_g1()
                : as_on_curve(blst_p1_affine_generator());
            const blst_p1_affine affine =
                detail::native<blst_p1_affine>(wide_value);
            blst_p1 p{};
            ::blst_p1_from_affine(&p, &affine);
            require(::blst_p1_on_curve(&p));
            require(::blst_p1_affine_on_curve(&affine));
            require(::blst_p1_in_g1(&p) == ::blst_p1_affine_in_g1(&affine));

            blst_p1 q{};
            ::blst_p1_double(&q, &p);
            blst_p1 complete{};
            blst_p1 incomplete{};
            ::blst_p1_add_or_double(&complete, &p, &q);
            ::blst_p1_add(&incomplete, &p, &q);
            require(::blst_p1_is_equal(&complete, &incomplete));
            require(::blst_p1_on_curve(&complete));

            blst_p1 equal_sum{};
            ::blst_p1_add_or_double(&equal_sum, &p, &p);
            require(::blst_p1_is_equal(&equal_sum, &q));
            blst_p1 opposite = p;
            ::blst_p1_cneg(&opposite, true);
            ::blst_p1_add_or_double(&equal_sum, &p, &opposite);
            require(::blst_p1_is_inf(&equal_sum));

            blst_p1_affine q_affine{};
            ::blst_p1_to_affine(&q_affine, &q);
            blst_p1 mixed{};
            ::blst_p1_add_affine(&mixed, &p, &q_affine);
            require(::blst_p1_is_equal(&mixed, &complete));
            ::blst_p1_add_or_double_affine(&mixed, &p, &q_affine);
            require(::blst_p1_is_equal(&mixed, &complete));
            ::blst_p1_sub(&mixed, &complete, &q);
            require(::blst_p1_is_equal(&mixed, &p));
            ::blst_p1_sub_affine(&mixed, &complete, &q_affine);
            require(::blst_p1_is_equal(&mixed, &p));

            blst_p1 normalized{};
            ::blst_p1_from_jacobian(&normalized, &complete);
            require(::blst_p1_is_equal(&normalized, &complete));
            blst_p1 roundtrip{};
            blst_p1_affine complete_affine{};
            ::blst_p1_to_affine(&complete_affine, &complete);
            ::blst_p1_from_affine(&roundtrip, &complete_affine);
            require(::blst_p1_is_equal(&roundtrip, &complete));
            require(::blst_p1_affine_is_inf(&complete_affine) ==
                    ::blst_p1_is_inf(&complete));

            blst_p1 affine_neg_projective{};
            blst_p1_affine affine_neg = complete_affine;
            ::blst_p1_affine_cneg(&affine_neg, true);
            ::blst_p1_from_affine(&affine_neg_projective, &affine_neg);
            blst_p1 projective_neg = complete;
            ::blst_p1_cneg(&projective_neg, true);
            require(::blst_p1_is_equal(&affine_neg_projective, &projective_neg));
            blst_p1 unchanged = complete;
            ::blst_p1_cneg(&unchanged, false);
            require(::blst_p1_is_equal(&unchanged, &complete));

            std::array<blst_byte, 96> serialized{};
            std::array<blst_byte, 96> affine_serialized{};
            std::array<blst_byte, 48> compressed{};
            std::array<blst_byte, 48> affine_compressed{};
            ::blst_p1_serialize(serialized.data(), &complete);
            ::blst_p1_affine_serialize(affine_serialized.data(), &complete_affine);
            ::blst_p1_compress(compressed.data(), &complete);
            ::blst_p1_affine_compress(affine_compressed.data(), &complete_affine);
            require(serialized == affine_serialized);
            require(compressed == affine_compressed);

            require(::blst_p1_is_equal(::blst_p1_generator(),
                                       ::blst_p1_generator()));
            require(::blst_p1_affine_is_equal(::blst_p1_affine_generator(),
                                              ::blst_p1_affine_generator()));

            const std::array<blst_p1, 2> batch_points{p, q};
            const std::array<const blst_p1*, 2> batch_heads{
                batch_points.data(), nullptr};
            std::array<blst_p1_affine, 2> batch_affine{};
            ::blst_p1s_to_affine(
                batch_affine.data(), batch_heads.data(), batch_points.size());
            require(::blst_p1_affine_is_equal(&batch_affine[0], &affine));
            require(::blst_p1_affine_is_equal(&batch_affine[1], &q_affine));
            const std::array<const blst_p1_affine*, 2> affine_heads{
                batch_affine.data(), nullptr};
            blst_p1 batch_sum{};
            ::blst_p1s_add(&batch_sum, affine_heads.data(), batch_affine.size());
            require(::blst_p1_is_equal(&batch_sum, &complete));

            const blst_p1 subgroup = *::blst_p1_generator();
            const blst_byte* scalar_pointer = nbits == 0 ? nullptr : scalar.data();
            blst_p1 ct{};
            blst_p1 vt{};
            blst_p1 unchecked{};
            ::blst_p1_mult(&ct, &subgroup, scalar_pointer, nbits);
            ::blst_p1_mult_vartime(&vt, &subgroup, scalar_pointer, nbits);
            ::blst_p1_unchecked_mult(&unchecked, &subgroup, scalar_pointer, nbits);
            require(::blst_p1_is_equal(&ct, &vt));
            require(::blst_p1_is_equal(&ct, &unchecked));
            blst_p1 alias = subgroup;
            ::blst_p1_mult(&alias, &alias, scalar_pointer, nbits);
            require(::blst_p1_is_equal(&alias, &ct)); // FUZZ_PATH: groups.exact_alias

            const std::size_t direct_bits = std::min<std::size_t>(nbits, 16);
            blst_p1 direct{};
            direct_mult(direct, subgroup, scalar, direct_bits);
            blst_p1 bounded{};
            ::blst_p1_mult(&bounded, &subgroup, scalar.data(), direct_bits);
            require(::blst_p1_is_equal(&direct, &bounded));
            direct_mult(direct, p, scalar, direct_bits);
            ::blst_p1_unchecked_mult(&bounded, &p, scalar.data(), direct_bits);
            require(::blst_p1_is_equal(&direct, &bounded));
        }

        void exercise_g2(
            std::span<const std::byte> input,
            const std::array<blst_byte, 32>& scalar,
            std::size_t nbits)
        {
            using detail::require;
            const bool use_wide = !input.empty() &&
                (std::to_integer<unsigned char>(input.back()) & 2U) != 0;
            const OnCurveG2AffineValue wide_value = use_wide
                ? retained_non_subgroup_g2()
                : as_on_curve(blst_p2_affine_generator());
            const blst_p2_affine affine =
                detail::native<blst_p2_affine>(wide_value);
            blst_p2 p{};
            ::blst_p2_from_affine(&p, &affine);
            require(::blst_p2_on_curve(&p));
            require(::blst_p2_affine_on_curve(&affine));
            require(::blst_p2_in_g2(&p) == ::blst_p2_affine_in_g2(&affine));

            blst_p2 q{};
            ::blst_p2_double(&q, &p);
            blst_p2 complete{};
            blst_p2 incomplete{};
            ::blst_p2_add_or_double(&complete, &p, &q);
            ::blst_p2_add(&incomplete, &p, &q);
            require(::blst_p2_is_equal(&complete, &incomplete));
            require(::blst_p2_on_curve(&complete));

            blst_p2 equal_sum{};
            ::blst_p2_add_or_double(&equal_sum, &p, &p);
            require(::blst_p2_is_equal(&equal_sum, &q));
            blst_p2 opposite = p;
            ::blst_p2_cneg(&opposite, true);
            ::blst_p2_add_or_double(&equal_sum, &p, &opposite);
            require(::blst_p2_is_inf(&equal_sum));

            blst_p2_affine q_affine{};
            ::blst_p2_to_affine(&q_affine, &q);
            blst_p2 mixed{};
            ::blst_p2_add_affine(&mixed, &p, &q_affine);
            require(::blst_p2_is_equal(&mixed, &complete));
            ::blst_p2_add_or_double_affine(&mixed, &p, &q_affine);
            require(::blst_p2_is_equal(&mixed, &complete));
            ::blst_p2_sub(&mixed, &complete, &q);
            require(::blst_p2_is_equal(&mixed, &p));
            ::blst_p2_sub_affine(&mixed, &complete, &q_affine);
            require(::blst_p2_is_equal(&mixed, &p));

            blst_p2 normalized{};
            ::blst_p2_from_jacobian(&normalized, &complete);
            require(::blst_p2_is_equal(&normalized, &complete));
            blst_p2 roundtrip{};
            blst_p2_affine complete_affine{};
            ::blst_p2_to_affine(&complete_affine, &complete);
            ::blst_p2_from_affine(&roundtrip, &complete_affine);
            require(::blst_p2_is_equal(&roundtrip, &complete));
            require(::blst_p2_affine_is_inf(&complete_affine) ==
                    ::blst_p2_is_inf(&complete));

            blst_p2 affine_neg_projective{};
            blst_p2_affine affine_neg = complete_affine;
            ::blst_p2_affine_cneg(&affine_neg, true);
            ::blst_p2_from_affine(&affine_neg_projective, &affine_neg);
            blst_p2 projective_neg = complete;
            ::blst_p2_cneg(&projective_neg, true);
            require(::blst_p2_is_equal(&affine_neg_projective, &projective_neg));

            std::array<blst_byte, 192> serialized{};
            std::array<blst_byte, 192> affine_serialized{};
            std::array<blst_byte, 96> compressed{};
            std::array<blst_byte, 96> affine_compressed{};
            ::blst_p2_serialize(serialized.data(), &complete);
            ::blst_p2_affine_serialize(affine_serialized.data(), &complete_affine);
            ::blst_p2_compress(compressed.data(), &complete);
            ::blst_p2_affine_compress(affine_compressed.data(), &complete_affine);
            require(serialized == affine_serialized);
            require(compressed == affine_compressed);
            require(::blst_p2_is_equal(::blst_p2_generator(),
                                       ::blst_p2_generator()));
            require(::blst_p2_affine_is_equal(::blst_p2_affine_generator(),
                                              ::blst_p2_affine_generator()));

            const std::array<blst_p2, 2> batch_points{p, q};
            const std::array<const blst_p2*, 2> batch_heads{
                batch_points.data(), nullptr};
            std::array<blst_p2_affine, 2> batch_affine{};
            ::blst_p2s_to_affine(
                batch_affine.data(), batch_heads.data(), batch_points.size());
            require(::blst_p2_affine_is_equal(&batch_affine[0], &affine));
            require(::blst_p2_affine_is_equal(&batch_affine[1], &q_affine));
            const std::array<const blst_p2_affine*, 2> affine_heads{
                batch_affine.data(), nullptr};
            blst_p2 batch_sum{};
            ::blst_p2s_add(&batch_sum, affine_heads.data(), batch_affine.size());
            require(::blst_p2_is_equal(&batch_sum, &complete));

            const blst_p2 subgroup = *::blst_p2_generator();
            const blst_byte* scalar_pointer = nbits == 0 ? nullptr : scalar.data();
            blst_p2 ct{};
            blst_p2 vt{};
            blst_p2 unchecked{};
            ::blst_p2_mult(&ct, &subgroup, scalar_pointer, nbits);
            ::blst_p2_mult_vartime(&vt, &subgroup, scalar_pointer, nbits);
            ::blst_p2_unchecked_mult(&unchecked, &subgroup, scalar_pointer, nbits);
            require(::blst_p2_is_equal(&ct, &vt));
            require(::blst_p2_is_equal(&ct, &unchecked));
            blst_p2 alias = subgroup;
            ::blst_p2_mult(&alias, &alias, scalar_pointer, nbits);
            require(::blst_p2_is_equal(&alias, &ct));

            const std::size_t direct_bits = std::min<std::size_t>(nbits, 16);
            blst_p2 direct{};
            direct_mult(direct, subgroup, scalar, direct_bits);
            blst_p2 bounded{};
            ::blst_p2_mult(&bounded, &subgroup, scalar.data(), direct_bits);
            require(::blst_p2_is_equal(&direct, &bounded));
            direct_mult(direct, p, scalar, direct_bits);
            ::blst_p2_unchecked_mult(&bounded, &p, scalar.data(), direct_bits);
            require(::blst_p2_is_equal(&direct, &bounded));
        }

        template<class Point, class Affine, class PointXyzz>
        void exercise_unique_native_group_methods()
        {
            using detail::require;
            const Point& generator = Point::generator();
            Point doubled{};
            // FUZZ_COVER: blst::scalar::Point::dadd
            Point::dadd(doubled, generator, generator);
            Point expected{};
            Point::doubling(expected, generator);
            require(Point::equal(doubled, expected));

            PointXyzz left{};
            PointXyzz right{};
            PointXyzz sum{};
            const Affine& affine = Affine::generator();
            PointXyzz::from_affine(left, affine);
            PointXyzz::from_affine(right, affine);
            // FUZZ_COVER: blst::scalar::Pointxyzz::dadd
            PointXyzz::dadd(sum, left, right);
            Point roundtrip{};
            PointXyzz::to_Jacobian(roundtrip, sum);
            require(Point::equal(roundtrip, expected));
        }
    }

    void run_groups_case(std::span<const std::byte> input)
    {
        const auto scalar = detail::fixed_bytes<32>(input);
        const std::size_t nbits = scalar_bits(input);
        exercise_g1(input, scalar, nbits);
        exercise_g2(input, scalar, nbits);
        exercise_unique_native_group_methods<
            blst::scalar::G1Point,
            blst::scalar::G1AffinePoint,
            blst::scalar::G1PointZZ>();
        exercise_unique_native_group_methods<
            blst::scalar::G2Point,
            blst::scalar::G2AffinePoint,
            blst::scalar::G2PointZZ>();
    }
}
