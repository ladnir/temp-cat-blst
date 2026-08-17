// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter_common.h"
#include "blst/scalar/curve/GT.h"

#include <algorithm>

namespace catblst::fuzz::adapter
{
    namespace
    {
        [[nodiscard]] std::pair<
            std::span<const std::byte>,
            std::span<const std::byte>>
        split(std::span<const std::byte> input) noexcept
        {
            const std::size_t middle = input.size() / 2;
            return {input.first(middle), input.subspan(middle)};
        }

        [[nodiscard]] bool zero(const blst_fr& value) noexcept
        {
            const blst_fr z{};
            return detail::object_equal(value, z);
        }

        [[nodiscard]] bool zero(const blst_fp& value) noexcept
        {
            const blst_fp z{};
            return detail::object_equal(value, z);
        }

        [[nodiscard]] bool zero(const blst_fp2& value) noexcept
        {
            const blst_fp2 z{};
            return detail::object_equal(value, z);
        }

        [[nodiscard]] bool zero(const blst_fp12& value) noexcept
        {
            const blst_fp12 z{};
            return detail::object_equal(value, z);
        }

        [[nodiscard]] std::size_t exponent_bits(
            const std::array<blst_byte, 32>& exponent) noexcept
        {
            for (std::size_t byte = exponent.size(); byte-- != 0;)
            {
                unsigned char value = exponent[byte];
                if (value == 0)
                    continue;
                std::size_t bits = byte * 8;
                while (value != 0)
                {
                    ++bits;
                    value >>= 1;
                }
                return bits;
            }
            return 0;
        }

        void exercise_fr(
            const blst_fr& a,
            const blst_fr& b,
            const std::array<blst_byte, 32>& exponent,
            std::size_t nbits,
            std::size_t shift)
        {
            using detail::object_equal;
            using detail::require;
            const std::array<std::uint64_t, 4> one_words{1, 0, 0, 0};
            blst_fr one{};
            ::blst_fr_from_uint64(&one, one_words.data());

            blst_fr result{};
            blst_fr intermediate{};
            ::blst_fr_add(&result, &a, &b);
            ::blst_fr_sub(&result, &result, &b);
            require(object_equal(result, a)); // FUZZ_PATH: fields.fr_add_sub
            ::blst_fr_sub(&result, &a, &b);
            ::blst_fr_add(&result, &result, &b);
            require(object_equal(result, a));
            ::blst_fr_sqr(&result, &a);
            ::blst_fr_mul(&intermediate, &a, &a);
            require(object_equal(result, intermediate));
            ::blst_fr_mul_by_3(&result, &a);
            ::blst_fr_add(&intermediate, &a, &a);
            ::blst_fr_add(&intermediate, &intermediate, &a);
            require(object_equal(result, intermediate));

            ::blst_fr_lshift(&result, &a, shift);
            ::blst_fr_rshift(&result, &result, shift);
            require(object_equal(result, a));
            ::blst_fr_cneg(&result, &a, false);
            require(object_equal(result, a));
            ::blst_fr_cneg(&result, &a, true);
            ::blst_fr_add(&result, &result, &a);
            require(zero(result));

            ::blst_fr_pow(
                &result,
                &a,
                nbits == 0 ? nullptr : exponent.data(),
                nbits);
            blst_fr pow_zero{};
            ::blst_fr_pow(&pow_zero, &a, nullptr, 0);
            require(object_equal(pow_zero, one));
            const std::array<blst_byte, 1> exponent_one{1};
            ::blst_fr_pow(&intermediate, &a, exponent_one.data(), 1);
            require(object_equal(intermediate, a));

            if (!zero(a))
            {
                blst_fr inverse{};
                blst_fr euclidean{};
                ::blst_fr_inverse(&inverse, &a);
                ::blst_fr_eucl_inverse(&euclidean, &a);
                require(object_equal(inverse, euclidean));
                ::blst_fr_mul(&result, &a, &inverse);
                require(object_equal(result, one));
            }

            blst_fr ordinary{};
            blst_fr montgomery{};
            ::blst_fr_from(&ordinary, &a);
            ::blst_fr_to(&montgomery, &ordinary);
            require(object_equal(montgomery, a));
            std::array<std::uint64_t, 4> exported{};
            ::blst_uint64_from_fr(exported.data(), &a);
            blst_fr imported{};
            ::blst_fr_from_uint64(&imported, exported.data());
            require(object_equal(imported, a));

            blst_scalar scalar_a{};
            blst_scalar scalar_b{};
            blst_scalar scalar_roundtrip{};
            ::blst_scalar_from_fr(&scalar_a, &a);
            ::blst_scalar_from_fr(&scalar_b, &b);
            const std::array<blst_fr, 2> batch_input{a, b};
            std::array<blst_scalar, 2> batch_output{};
            require(::blst_scalars_from_fr(
                batch_output.data(), batch_input.data(), batch_input.size()));
            require(object_equal(batch_output[0], scalar_a));
            require(object_equal(batch_output[1], scalar_b));
            if (::blst_sk_check(&scalar_a) && ::blst_sk_check(&scalar_b))
            {
                blst_scalar scalar_result{};
                blst_fr field_result{};
                bool nonzero = ::blst_sk_add_n_check(
                    &scalar_result, &scalar_a, &scalar_b);
                ::blst_fr_add(&field_result, &a, &b);
                ::blst_scalar_from_fr(&scalar_roundtrip, &field_result);
                require(object_equal(scalar_result, scalar_roundtrip));
                require(nonzero == ::blst_sk_check(&scalar_result));

                nonzero = ::blst_sk_sub_n_check(
                    &scalar_result, &scalar_a, &scalar_b);
                ::blst_fr_sub(&field_result, &a, &b);
                ::blst_scalar_from_fr(&scalar_roundtrip, &field_result);
                require(object_equal(scalar_result, scalar_roundtrip));
                require(nonzero == ::blst_sk_check(&scalar_result));

                nonzero = ::blst_sk_mul_n_check(
                    &scalar_result, &scalar_a, &scalar_b);
                ::blst_fr_mul(&field_result, &a, &b);
                ::blst_scalar_from_fr(&scalar_roundtrip, &field_result);
                require(object_equal(scalar_result, scalar_roundtrip));
                require(nonzero == ::blst_sk_check(&scalar_result));

                blst_scalar scalar_inverse{};
                ::blst_sk_inverse(&scalar_inverse, &scalar_a);
                require(::blst_sk_mul_n_check(
                    &scalar_result, &scalar_a, &scalar_inverse));
                std::array<blst_byte, 32> inverse_product{};
                ::blst_lendian_from_scalar(inverse_product.data(), &scalar_result);
                require(inverse_product[0] == 1);
                for (std::size_t index = 1; index < inverse_product.size(); ++index)
                    require(inverse_product[index] == 0);
            }

            blst_fr ct0 = a;
            blst_fr ct1 = b;
            blst_fr expected0{};
            blst_fr expected1{};
            ::blst_fr_mul(&intermediate, &b, &one);
            ::blst_fr_add(&expected0, &a, &intermediate);
            ::blst_fr_sub(&expected1, &a, &intermediate);
            ::blst_fr_ct_bfly(&ct0, &ct1, &one);
            require(object_equal(ct0, expected0));
            require(object_equal(ct1, expected1));

            blst_fr gs0 = a;
            blst_fr gs1 = b;
            ::blst_fr_add(&expected0, &a, &b);
            ::blst_fr_sub(&expected1, &a, &b);
            ::blst_fr_gs_bfly(&gs0, &gs1, &one);
            require(object_equal(gs0, expected0));
            require(object_equal(gs1, expected1));
        }

        void exercise_fp(const blst_fp& a, const blst_fp& b, std::size_t shift)
        {
            using detail::object_equal;
            using detail::require;
            const std::array<std::uint64_t, 6> one_words{1, 0, 0, 0, 0, 0};
            blst_fp one{};
            ::blst_fp_from_uint64(&one, one_words.data());
            blst_fp result{};
            blst_fp intermediate{};
            ::blst_fp_add(&result, &a, &b);
            ::blst_fp_sub(&result, &result, &b);
            require(object_equal(result, a));
            ::blst_fp_sqr(&result, &a);
            ::blst_fp_mul(&intermediate, &a, &a);
            require(object_equal(result, intermediate));
            ::blst_fp_mul_by_3(&result, &a);
            ::blst_fp_add(&intermediate, &a, &a);
            ::blst_fp_add(&intermediate, &intermediate, &a);
            require(object_equal(result, intermediate));
            ::blst_fp_mul_by_8(&result, &a);
            intermediate = a;
            for (unsigned i = 0; i < 3; ++i)
                ::blst_fp_add(&intermediate, &intermediate, &intermediate);
            require(object_equal(result, intermediate));
            ::blst_fp_lshift(&result, &a, shift);
            intermediate = a;
            for (std::size_t i = 0; i < shift; ++i)
                ::blst_fp_add(&intermediate, &intermediate, &intermediate);
            require(object_equal(result, intermediate));
            ::blst_fp_cneg(&result, &a, false);
            require(object_equal(result, a));
            ::blst_fp_cneg(&result, &a, true);
            ::blst_fp_add(&result, &result, &a);
            require(zero(result));

            if (!zero(a))
            {
                blst_fp inverse{};
                blst_fp euclidean{};
                ::blst_fp_inverse(&inverse, &a);
                ::blst_fp_eucl_inverse(&euclidean, &a);
                require(object_equal(inverse, euclidean));
                ::blst_fp_mul(&result, &a, &inverse);
                require(object_equal(result, one));
            }

            const bool square = ::blst_fp_is_square(&a);
            const bool sqrt_ok = ::blst_fp_sqrt(&result, &a);
            require(square == sqrt_ok);
            if (sqrt_ok)
            {
                ::blst_fp_sqr(&intermediate, &result);
                require(object_equal(intermediate, a));
            }

            blst_fp ordinary{};
            blst_fp montgomery{};
            ::blst_fp_from(&ordinary, &a);
            ::blst_fp_to(&montgomery, &ordinary);
            require(object_equal(montgomery, a));
            std::array<std::uint64_t, 6> words64{};
            std::array<std::uint32_t, 12> words32{};
            ::blst_uint64_from_fp(words64.data(), &a);
            ::blst_fp_from_uint64(&montgomery, words64.data());
            require(object_equal(montgomery, a));
            ::blst_uint32_from_fp(words32.data(), &a);
            ::blst_fp_from_uint32(&montgomery, words32.data());
            require(object_equal(montgomery, a));
        }

        void exercise_fp2(const blst_fp2& a, const blst_fp2& b, std::size_t shift)
        {
            using detail::object_equal;
            using detail::require;
            const std::array<std::uint64_t, 6> one_words{1, 0, 0, 0, 0, 0};
            blst_fp2 one{};
            ::blst_fp_from_uint64(&one.fp[0], one_words.data());
            blst_fp2 result{};
            blst_fp2 intermediate{};
            ::blst_fp2_add(&result, &a, &b);
            ::blst_fp2_sub(&result, &result, &b);
            require(object_equal(result, a));
            ::blst_fp2_sqr(&result, &a);
            ::blst_fp2_mul(&intermediate, &a, &a);
            require(object_equal(result, intermediate));
            ::blst_fp2_mul_by_3(&result, &a);
            ::blst_fp2_add(&intermediate, &a, &a);
            ::blst_fp2_add(&intermediate, &intermediate, &a);
            require(object_equal(result, intermediate));
            ::blst_fp2_mul_by_8(&result, &a);
            intermediate = a;
            for (unsigned i = 0; i < 3; ++i)
                ::blst_fp2_add(&intermediate, &intermediate, &intermediate);
            require(object_equal(result, intermediate));
            ::blst_fp2_lshift(&result, &a, shift);
            intermediate = a;
            for (std::size_t i = 0; i < shift; ++i)
                ::blst_fp2_add(&intermediate, &intermediate, &intermediate);
            require(object_equal(result, intermediate));
            ::blst_fp2_cneg(&result, &a, false);
            require(object_equal(result, a));
            ::blst_fp2_cneg(&result, &a, true);
            ::blst_fp2_add(&result, &result, &a);
            require(zero(result));
            if (!zero(a))
            {
                blst_fp2 inverse{};
                blst_fp2 euclidean{};
                ::blst_fp2_inverse(&inverse, &a);
                ::blst_fp2_eucl_inverse(&euclidean, &a);
                require(object_equal(inverse, euclidean));
                ::blst_fp2_mul(&result, &a, &inverse);
                require(object_equal(result, one));
            }
            const bool square = ::blst_fp2_is_square(&a);
            const bool sqrt_ok = ::blst_fp2_sqrt(&result, &a);
            require(square == sqrt_ok);
            if (sqrt_ok)
            {
                ::blst_fp2_sqr(&intermediate, &result);
                require(object_equal(intermediate, a));
            }
        }

        void exercise_fp12_and_gt(
            const blst_fp12& a,
            const blst_fp12& b,
            const std::array<blst_byte, 32>& scalar,
            std::size_t nbits,
            std::size_t frobenius)
        {
            using detail::object_equal;
            using detail::require;
            blst_fp12 result{};
            blst_fp12 intermediate{};
            ::blst_fp12_sqr(&result, &a);
            ::blst_fp12_mul(&intermediate, &a, &a);
            require(::blst_fp12_is_equal(&result, &intermediate));
            ::blst_fp12_mul(&result, &a, ::blst_fp12_one());
            require(::blst_fp12_is_equal(&result, &a));
            intermediate = a;
            ::blst_fp12_conjugate(&intermediate);
            ::blst_fp12_conjugate(&intermediate);
            require(::blst_fp12_is_equal(&intermediate, &a));
            if (!zero(a))
            {
                ::blst_fp12_inverse(&intermediate, &a);
                ::blst_fp12_mul(&result, &a, &intermediate);
                require(::blst_fp12_is_one(&result));
            }

            blst_fp6 sparse{};
            sparse.fp2[0] = b.fp6[0].fp2[0];
            sparse.fp2[1] = b.fp6[0].fp2[1];
            sparse.fp2[2] = b.fp6[1].fp2[1];
            blst_fp12 sparse_full{};
            sparse_full.fp6[0].fp2[0] = sparse.fp2[0];
            sparse_full.fp6[0].fp2[1] = sparse.fp2[1];
            sparse_full.fp6[1].fp2[1] = sparse.fp2[2];
            ::blst_fp12_mul_by_xy00z0(&result, &a, &sparse);
            ::blst_fp12_mul(&intermediate, &a, &sparse_full);
            require(::blst_fp12_is_equal(&result, &intermediate));

            blst_fp12 before{};
            std::memset(&before, 0xa5, sizeof(before));
            result = before;
            const bool frobenius_ok =
                ::blst_fp12_frobenius_map(&result, &a, frobenius);
            require(frobenius_ok == (frobenius >= 1 && frobenius <= 3));
            if (!frobenius_ok)
                require(object_equal(result, before));

            const blst_fp12& source = !zero(a) ? a : *::blst_fp12_one();
            blst_fp12 gt{};
            ::blst_final_exp(&gt, &source);
            require(::blst_fp12_in_group(&gt));
            ::blst_fp12_cyclotomic_sqr(&result, &gt);
            ::blst_fp12_sqr(&intermediate, &gt);
            require(::blst_fp12_is_equal(&result, &intermediate));

            require(::blst_gt_is_identity(::blst_gt_identity()));
            require(!::blst_gt_is_identity(::blst_gt_generator()));
            require(::blst_fp12_in_group(::blst_gt_generator()));
            ::blst_gt_add(&result, &gt, ::blst_gt_generator());
            ::blst_fp12_mul(&intermediate, &gt, ::blst_gt_generator());
            require(::blst_gt_is_equal(&result, &intermediate));
            ::blst_gt_sub(&result, &result, ::blst_gt_generator());
            require(::blst_gt_is_equal(&result, &gt));
            ::blst_gt_neg(&result, &gt);
            ::blst_gt_add(&result, &result, &gt);
            require(::blst_gt_is_identity(&result));
            ::blst_gt_double(&result, &gt);
            ::blst_gt_add(&intermediate, &gt, &gt);
            require(::blst_gt_is_equal(&result, &intermediate));

            const blst_byte* scalar_pointer = nbits == 0 ? nullptr : scalar.data();
            ::blst_gt_mult(&result, &gt, scalar_pointer, nbits);
            ::blst_gt_mult_vartime(&intermediate, &gt, scalar_pointer, nbits);
            require(::blst_gt_is_equal(&result, &intermediate));
            blst_fp12 alias = gt;
            ::blst_gt_mult(&alias, &alias, scalar_pointer, nbits);
            require(::blst_gt_is_equal(&alias, &result));
        }

        void exercise_unique_native_fields()
        {
            using blst::scalar::GT;
            using detail::require;

            const GT identity = GT::identity();
            const GT generator = GT::generator();
            GT result{};
            // FUZZ_COVER: blst::scalar::GT::dadd
            GT::dadd(result, identity, generator);
            require(GT::equal(result, generator));

            std::array<unsigned char, 12 * 48> encoded{};
            GT::serialize(encoded.data(), generator);
            GT decoded = identity;
            // FUZZ_COVER: blst::scalar::GT::deserialize_unchecked
            require(GT::deserialize_unchecked(decoded, encoded.data()));
            require(GT::equal(decoded, generator));
        }
    }

    void run_fields_case(std::span<const std::byte> input)
    {
        input = input.first(std::min<std::size_t>(input.size(), 256));
        const auto [left, right] = split(input);
        const blst_fr fr_a = detail::native<blst_fr>(
            blst_scalar_from_le_bytes(left));
        const blst_fr fr_b = detail::native<blst_fr>(
            blst_scalar_from_le_bytes(right));
        const blst_fp fp_a = detail::native<blst_fp>(
            blst_fp_from_le_bytes(left));
        const blst_fp fp_b = detail::native<blst_fp>(
            blst_fp_from_le_bytes(right));
        const blst_fp2 fp2_a = detail::native<blst_fp2>(reduce_fp2(left));
        const blst_fp2 fp2_b = detail::native<blst_fp2>(reduce_fp2(right));
        const blst_fp12 fp12_a = detail::native<blst_fp12>(reduce_fp12(left));
        const blst_fp12 fp12_b = detail::native<blst_fp12>(reduce_fp12(right));
        const auto scalar = detail::fixed_bytes<32>(input);
        const std::size_t nbits = input.empty()
            ? 0
            : std::to_integer<unsigned char>(input.front()) % 257U;
        const std::size_t shift = input.empty()
            ? 0
            : std::to_integer<unsigned char>(input.back()) % 17U;
        const std::size_t frobenius = input.size() < 2
            ? 0
            : std::to_integer<unsigned char>(input[1]) % 6U;

        exercise_fr(fr_a, fr_b, scalar, exponent_bits(scalar), shift);
        exercise_fp(fp_a, fp_b, shift);
        exercise_fp2(fp2_a, fp2_b, shift);
        exercise_fp12_and_gt(fp12_a, fp12_b, scalar, nbits, frobenius);
        exercise_unique_native_fields();
    }
}
