// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "../generated_matrixnok_constants.hpp"
#include "matrixnok_change_base.h"
#include "matrixnok_fixed.h"
#include "ring_element.h"

#include <array>
#include <type_traits>

namespace blst::avx::matrixnok
{
    class Fp2RingCandidate;

    class FpRingCandidate
    {
    public:
        static constexpr int limbs = 8;
        static constexpr int element_bits = 50;
        static constexpr int standard_bound = 2;
        static constexpr int quotient_error_bound = 16;

        using StandardResidue = BoundedElement<Bounds<0, standard_bound>, limbs, element_bits>;

        struct Element
        {
            StandardResidue m1;
            StandardResidue m2;
        };

    private:
        static INLINE StandardResidue expand_m2(const StandardResidue& m2)
        {
            AVXVector<limbs> high(0);
            AVXVector<limbs> low(0);
            detail::R2Matrix::reduce_raw_shifted(m2.data, high, low);

            using WideBound = Bounds<0, 2 * limbs * standard_bound>;
            using Bounded = BoundedElement<WideBound, limbs, element_bits>;
            const WideElement<WideBound, limbs, element_bits> wide{
                Bounded(high), Bounded(low)};
            return detail::ReducerM1::template reduce_auto<standard_bound>(
                detail::ReducerM1::reduce_wide(wide));
        }

        template<std::size_t batch>
        static INLINE std::array<StandardResidue, batch> expand_m2_batch(
            const std::array<StandardResidue, batch>& m2)
        {
            std::array<AVXVector<limbs>, batch> residues;
            std::array<AVXVector<limbs>, batch> high{};
            std::array<AVXVector<limbs>, batch> low{};
            BLST_UNROLL
            for (std::size_t index = 0; index < batch; ++index) {
                residues[index] = m2[index].data;
            }
            detail::R2Matrix::template reduce_raw_batch_shifted<batch>(
                residues, high, low);

            using WideBound = Bounds<0, 2 * limbs * standard_bound>;
            using Bounded = BoundedElement<WideBound, limbs, element_bits>;
            std::array<StandardResidue, batch> result;
            BLST_UNROLL
            for (std::size_t index = 0; index < batch; ++index) {
                const WideElement<WideBound, limbs, element_bits> wide{
                    Bounded(high[index]), Bounded(low[index])};
                result[index] = detail::ReducerM1::template reduce_auto<standard_bound>(
                    detail::ReducerM1::reduce_wide(wide));
            }
            return result;
        }

        // Finish a fixed batch of already accumulated wide products. The
        // caller supplies the exact compile-time product bound; this is the
        // narrow bridge used by extension-field formulas to combine terms
        // before paying for reduction and expansion.
        template<int product_bound, std::size_t batch>
        static INLINE std::array<Element, batch> finish_wide_products(
            std::array<AVXVector<limbs>, batch>& m1_high,
            std::array<AVXVector<limbs>, batch>& m1_low,
            std::array<AVXVector<limbs>, batch>& m2_high,
            std::array<AVXVector<limbs>, batch>& m2_low)
        {
            static_assert(product_bound > 0);
            static_assert(batch > 0);
            using ProductBound = Bounds<0, product_bound>;
            using ProductElement = BoundedElement<ProductBound, limbs, element_bits>;

            std::array<AVXVector<limbs>, batch> residues;
            BLST_UNROLL
            for (std::size_t index = 0; index < batch; ++index) {
                const WideElement<ProductBound, limbs, element_bits> product{
                    ProductElement(m1_high[index]), ProductElement(m1_low[index])};
                residues[index] = detail::ReducerM1::template reduce_auto<standard_bound>(
                    detail::ReducerM1::reduce_wide(product)).data;
            }

            r1_change_base_raw_batch(residues, m2_high, m2_low);

            constexpr int change_base_bound = limbs * standard_bound;
            constexpr int combined_bound = product_bound + change_base_bound;
            static_assert(combined_bound < 40);
            using CombinedBound = Bounds<0, combined_bound>;
            using CombinedElement = BoundedElement<CombinedBound, limbs, element_bits>;
            std::array<StandardResidue, batch> m2;
            BLST_UNROLL
            for (std::size_t index = 0; index < batch; ++index) {
                const WideElement<CombinedBound, limbs, element_bits> combined{
                    CombinedElement(m2_high[index]), CombinedElement(m2_low[index])};
                m2[index] = detail::ReducerM2::template reduce_auto<standard_bound>(
                    detail::ReducerM2::reduce_wide(combined));
            }

            const auto m1 = expand_m2_batch(m2);
            std::array<Element, batch> result;
            BLST_UNROLL
            for (std::size_t index = 0; index < batch; ++index) {
                result[index] = {m1[index], m2[index]};
            }
            return result;
        }

        static INLINE Element negate_for_product(const Element& value)
        {
            return {
                StandardResidue(
                    detail::ReducerM1::template moduli<1>().data.sub(value.m1.data)),
                StandardResidue(
                    detail::ReducerM2::template moduli<1>().data.sub(value.m2.data))};
        }

    public:
        constexpr FpRingCandidate() noexcept = default;

        static INLINE Element from_canonical_50bit_digits(
            const std::array<uint64_t, limbs>& digits)
        {
            AVXVector<limbs> high(0);
            AVXVector<limbs> low(0);
            detail::ConvertToMatrix::reduce_raw_shifted(
                AVXVector<limbs>(digits), high, low);

            using ConvertBound = Bounds<0, 800>;
            using ConvertElement = BoundedElement<ConvertBound, limbs, element_bits>;
            const WideElement<ConvertBound, limbs, element_bits> converted{
                ConvertElement(high), ConvertElement(low)};
            const StandardResidue m2 = detail::ReducerM2::template reduce_auto<standard_bound>(
                detail::ReducerM2::reduce_wide(converted));
            return {expand_m2(m2), m2};
        }

        static INLINE Element multiply(const Element& left, const Element& right)
        {
            const AVXVector<limbs> zero(0);
            const AVXVector<limbs> m1_high = zero.mulhi(left.m1.data, right.m1.data);
            const AVXVector<limbs> m1_low = zero.mullo(left.m1.data, right.m1.data);
            const AVXVector<limbs> m2_high = zero.mulhi(left.m2.data, right.m2.data);
            const AVXVector<limbs> m2_low = zero.mullo(left.m2.data, right.m2.data);

            using ProductBound = Bounds<0, standard_bound * standard_bound>;
            using ProductElement = BoundedElement<ProductBound, limbs, element_bits>;
            const WideElement<ProductBound, limbs, element_bits> m1_product{
                ProductElement(m1_high), ProductElement(m1_low)};
            const StandardResidue m1_reduced =
                detail::ReducerM1::template reduce_auto<standard_bound>(
                    detail::ReducerM1::reduce_wide(m1_product));

            AVXVector<limbs> residues = m1_reduced.data;
            AVXVector<limbs> high = m2_high;
            AVXVector<limbs> low = m2_low;
            r1_change_base_raw(residues, high, low);

            constexpr int change_base_bound = limbs * standard_bound;
            constexpr int combined_bound = standard_bound * standard_bound + change_base_bound;
            static_assert(combined_bound == 20);
            static_assert(combined_bound < 40);
            using CombinedBound = Bounds<0, combined_bound>;
            using CombinedElement = BoundedElement<CombinedBound, limbs, element_bits>;
            const WideElement<CombinedBound, limbs, element_bits> combined{
                CombinedElement(high), CombinedElement(low)};
            const StandardResidue m2 = detail::ReducerM2::template reduce_auto<standard_bound>(
                detail::ReducerM2::reduce_wide(combined));
            const StandardResidue m1 = expand_m2(m2);
            return {m1, m2};
        }

        // Multiply modulo p. Exact output aliasing with either input is
        // supported; partial overlap is not.
        static INLINE void multiply(
            Element& out, const Element& left, const Element& right)
        {
            out = multiply(left, right);
        }

        template<std::size_t batch>
        static INLINE std::array<Element, batch> multiply_batch(
            const std::array<Element, batch>& left,
            const std::array<Element, batch>& right)
        {
            static_assert(batch > 0, "MatrixNoK field multiplication requires a nonempty fixed batch");
            const AVXVector<limbs> zero(0);
            std::array<AVXVector<limbs>, batch> residues;
            std::array<AVXVector<limbs>, batch> high;
            std::array<AVXVector<limbs>, batch> low;

            BLST_UNROLL
            for (std::size_t index = 0; index < batch; ++index) {
                const AVXVector<limbs> m1_high = zero.mulhi(left[index].m1.data, right[index].m1.data);
                const AVXVector<limbs> m1_low = zero.mullo(left[index].m1.data, right[index].m1.data);
                high[index] = zero.mulhi(left[index].m2.data, right[index].m2.data);
                low[index] = zero.mullo(left[index].m2.data, right[index].m2.data);

                using ProductBound = Bounds<0, standard_bound * standard_bound>;
                using ProductElement = BoundedElement<ProductBound, limbs, element_bits>;
                const WideElement<ProductBound, limbs, element_bits> m1_product{
                    ProductElement(m1_high), ProductElement(m1_low)};
                residues[index] = detail::ReducerM1::template reduce_auto<standard_bound>(
                    detail::ReducerM1::reduce_wide(m1_product)).data;
            }

            r1_change_base_raw_batch(residues, high, low);

            constexpr int change_base_bound = limbs * standard_bound;
            constexpr int combined_bound = standard_bound * standard_bound + change_base_bound;
            static_assert(combined_bound == 20);
            static_assert(combined_bound < 40);
            using CombinedBound = Bounds<0, combined_bound>;
            using CombinedElement = BoundedElement<CombinedBound, limbs, element_bits>;
            std::array<StandardResidue, batch> m2;
            BLST_UNROLL
            for (std::size_t index = 0; index < batch; ++index) {
                const WideElement<CombinedBound, limbs, element_bits> combined{
                    CombinedElement(high[index]), CombinedElement(low[index])};
                m2[index] = detail::ReducerM2::template reduce_auto<standard_bound>(
                    detail::ReducerM2::reduce_wide(combined));
            }

            const auto m1 = expand_m2_batch(m2);
            std::array<Element, batch> result;
            BLST_UNROLL
            for (std::size_t index = 0; index < batch; ++index) {
                result[index] = {m1[index], m2[index]};
            }
            return result;
        }

        static INLINE void convert_from_raw(
            const Element& value,
            AVXVector<8>& out_high,
            AVXVector<8>& out_low)
        {
            out_high = AVXVector<8>(0);
            out_low = AVXVector<8>(0);
            detail::ConvertFromMatrix::reduce_raw(value.m2.data, out_high, out_low);
        }

        static constexpr const std::array<uint64_t, limbs>& moduli1() noexcept
        {
            return generated::moduli1;
        }

        static constexpr const std::array<uint64_t, limbs>& moduli2() noexcept
        {
            return generated::moduli2;
        }

    private:
        friend class Fp2RingCandidate;

    };

    static_assert(std::is_empty_v<FpRingCandidate>);
    static_assert(std::is_trivially_default_constructible_v<FpRingCandidate>);
    static_assert(std::is_trivially_copyable_v<FpRingCandidate>);
} // namespace blst::avx::matrixnok
