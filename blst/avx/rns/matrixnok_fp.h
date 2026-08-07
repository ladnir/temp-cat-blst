// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "../generated_matrixnok_constants.hpp"
#include "change_base.h"
#include "matrixnok_change_base.h"
#include "matrixnok_fixed.h"
#include "ring_element.h"
#include "../../common/bytes.h"
#include "../../fields/Fp.h"

#include <array>
#include <tuple>
#include <type_traits>
#include <utility>

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
        static constexpr int64_t max_add = 800;
        static constexpr int64_t max_rns_multiple = 4309;
        static constexpr int64_t max_negative_rns = 1932;
        static constexpr int64_t base_redundant_multiple = 40;

        static constexpr uint64_t MULT_OK_BOUND = uint64_t{1} << (52 - element_bits);
        static constexpr uint64_t STANDARD_BOUND = standard_bound;
        static constexpr int64_t MAX_ADD = max_add;
        using StandardResidue = BoundedElement<Bounds<0, standard_bound>, limbs, element_bits>;
        using FullyReducedResidue = BoundedElement<Bounds<0, 1>, limbs, element_bits>;
        using MultOkResidue = BoundedElement<Bounds<0, MULT_OK_BOUND>, limbs, element_bits>;
        using StandardBoundedElement = StandardResidue;
        using StandardElement = ExpandedElement<StandardResidue, Bounds<0, 1>>;
        using Element = StandardElement;
        using SmallStandardElement = SmallElement<StandardResidue, Bounds<0, 1>>;
        template<int64_t element_bound, int64_t rns_bound>
        using WideOffsetElement = ExpandedElement<
            WideElement<Bounds<element_bound, element_bound + 1>, limbs, element_bits>,
            Bounds<rns_bound, rns_bound>>;
        template<int64_t upper_bound>
        using ReducedElement = ExpandedElement<
            StandardBoundedElement, Bounds<0, upper_bound>>;
        template<int64_t upper_bound>
        using SmallReducedElement = SmallElement<
            StandardBoundedElement, Bounds<0, upper_bound>>;
        template<int64_t A>
        using ReadyElement = ReadyToReduce<
            WideElement<Bounds<0, A>, limbs, element_bits>, StandardResidue>;
        template<int64_t rns_bound>
        using OffsetElement = ExpandedElement<FullyReducedResidue, Bounds<rns_bound, rns_bound>>;

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
                result[index] = Element(m1[index], m2[index]);
            }
            return result;
        }

        static INLINE Element negate_for_product(const Element& value)
        {
            return Element(
                StandardResidue(
                    detail::ReducerM1::template moduli<1>().data.sub(value.m1.data)),
                StandardResidue(
                    detail::ReducerM2::template moduli<1>().data.sub(value.m2.data)));
        }

    public:
        constexpr FpRingCandidate() noexcept = default;

        static INLINE Element from_canonical_limbs(const blst::Fp& canonical)
        {
            constexpr uint64_t mask = (uint64_t{1} << element_bits) - 1;
            std::array<uint64_t, limbs> digits = {};
            BLST_UNROLL
            for (std::size_t digit = 0; digit < digits.size(); ++digit) {
                const std::size_t bit = digit * element_bits;
                const std::size_t limb = bit / 64;
                const unsigned shift = static_cast<unsigned>(bit % 64);
                uint64_t value = canonical.l[limb] >> shift;
                if (shift != 0 && limb + 1 < blst::Fp::num_limb()) {
                    value |= canonical.l[limb + 1] << (64 - shift);
                }
                digits[digit] = value & mask;
            }
            return from_canonical_50bit_digits(digits);
        }

        static INLINE Element from_scalar(const blst::Fp& montgomery)
        {
            blst::Fp canonical = {};
            blst::Fp::from_mont(canonical, montgomery);
            return from_canonical_limbs(canonical);
        }

        template<int bytes>
        static INLINE Element from_hex(const char (&literal)[bytes])
        {
            const char* hex = literal;
            if constexpr (bytes >= 3) {
                if (literal[0] == '0' &&
                    (literal[1] == 'x' || literal[1] == 'X')) {
                    hex += 2;
                }
            }
            blst::Fp canonical = {};
            const bool parsed = blst::limbs_from_hexascii(
                canonical.l, sizeof(canonical), hex);
            BLST_ASSERT(parsed);
            (void)parsed;
            return from_canonical_limbs(canonical);
        }

        static INLINE StandardElement zero()
        {
            return StandardElement{};
        }

        static INLINE StandardElement one()
        {
            return StandardElement(
                StandardResidue(AVXVector<limbs>(generated::one_m1)),
                StandardResidue(AVXVector<limbs>(generated::one_m2)));
        }

        template<int multiple = 1>
        static INLINE auto elementwise_moduli()
        {
            static_assert(multiple >= 1);
            constexpr int index = ceil_log2(multiple);
            auto m1 = detail::ReducerM1::template moduli<index>();
            auto m2 = detail::ReducerM2::template moduli<index>();
            return ExpandedElement<decltype(m1), Bounds<0, 0>>(m1, m2);
        }

        template<int64_t lower, int64_t upper, class RNSBoundsType>
        static INLINE auto wideify(
            const ExpandedElement<
                BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
                RNSBoundsType>& element)
        {
            auto zero = BoundedElement<Bounds<0, 0>, limbs, element_bits>::zero()
                .template recast<lower, upper>();
            using Wide = WideElement<Bounds<lower, upper>, limbs, element_bits>;
            return ExpandedElement<Wide, RNSBoundsType>(
                Wide(zero, element.m1), Wide(zero, element.m2));
        }

        static INLINE auto negative_range_offset()
        {
            const auto offset = OffsetElement<max_negative_rns>(
                FullyReducedResidue(AVXVector<limbs>(generated::max_neg_offset_m1)),
                FullyReducedResidue(AVXVector<limbs>(generated::max_neg_offset_m2)));
            return wideify(offset);
        }

        template<int64_t lower, int64_t upper, class RNSBoundsType,
            int to = MULT_OK_BOUND>
        static INLINE auto prep(const ExpandedElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            RNSBoundsType>& element)
        {
            auto m1 = detail::ReducerM1::template reduce_auto<to>(element.m1);
            auto m2 = detail::ReducerM2::template reduce_auto<to>(element.m2);
            return ExpandedElement<decltype(m1), RNSBoundsType>(m1, m2);
        }

        template<class ElementType>
        static INLINE auto prep_left(const ElementType& element)
        {
            return prep(element);
        }

        template<int64_t lower, int64_t upper, class RNSBoundsType,
            int to = MULT_OK_BOUND>
        static INLINE auto prep(const SmallElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            RNSBoundsType>& element)
        {
            auto m2 = detail::ReducerM2::template reduce_auto<to>(element.element);
            return SmallElement<decltype(m2), RNSBoundsType>(m2);
        }

        template<int64_t upper, int64_t lower, class RNSBoundsType>
        static INLINE auto negate(const ExpandedElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            RNSBoundsType>& value)
        {
            return elementwise_moduli<upper>() - value;
        }

        template<int64_t lower, int64_t upper, class RNSBoundsType>
        static INLINE auto negate(const SmallElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            RNSBoundsType>& value)
        {
            const auto moduli = elementwise_moduli<upper>();
            const SmallElement<decltype(moduli.m2), Bounds<0, 0>> offset(
                moduli.m2);
            return offset - value;
        }

        template<int64_t lower, int64_t upper>
        static INLINE StandardElement standard_negate(const ExpandedElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            Bounds<0, 1>>& value)
        {
            const auto offset = OffsetElement<1>(
                FullyReducedResidue(AVXVector<limbs>(generated::modulus_in_rns_m1)),
                FullyReducedResidue(AVXVector<limbs>(generated::modulus_in_rns_m2)));
            const auto negated = negate(value) + offset;
            return prep<decltype(negated.m1)::bounds.lower,
                decltype(negated.m1)::bounds.upper, Bounds<0, 1>, STANDARD_BOUND>(negated);
        }

        template<int64_t lower, int64_t upper>
        static INLINE auto standard_negate(const ExpandedElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            Bounds<0, 6>>& value)
        {
            const auto offset = OffsetElement<6>(
                FullyReducedResidue(AVXVector<limbs>(generated::six_modulus_in_rns_m1)),
                FullyReducedResidue(AVXVector<limbs>(generated::six_modulus_in_rns_m2)));
            const auto negated = negate(value) + offset;
            return prep<decltype(negated.m1)::bounds.lower,
                decltype(negated.m1)::bounds.upper, Bounds<0, 6>, STANDARD_BOUND>(negated);
        }

        template<class ElementType, int64_t rns_lower, int64_t rns_upper>
        static INLINE auto complete(const ExpandedElement<
            ElementType, Bounds<rns_lower, rns_upper>>& value)
        {
            if constexpr (rns_lower < 0) {
                return negative_range_offset() + value;
            } else {
                return value.complete();
            }
        }

        template<int64_t A, class ElementType, class RNSBoundsType>
        static INLINE ReadyElement<A> ready(
            const ExpandedElement<ElementType, RNSBoundsType>& value)
        {
            const auto completed = complete(value);
            using Completed = decltype(completed);
            static_assert(Completed::rns_bounds.lower >= 0);
            static_assert(Completed::rns_bounds.upper <= max_rns_multiple);
            static_assert(decltype(completed.m2.low)::bounds.lower >= 0);
            static_assert(decltype(completed.m2.high)::bounds.lower >= 0);
            const auto m2 = WideElement<Bounds<0, A>, limbs, element_bits>(
                completed.m2.high.template recast<0, A>(),
                completed.m2.low.template recast<0, A>());
            const auto m1_reduced =
                detail::ReducerM1::template reduce_auto<STANDARD_BOUND>(completed.m1);
            return ReadyElement<A>(
                m2, m1_reduced.template recast<0, STANDARD_BOUND>());
        }

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
            return Element(expand_m2(m2), m2);
        }

        template<int batch, int64_t rns_bound>
        static INLINE std::array<
            ExpandedElement<StandardResidue, Bounds<0, rns_bound>>, batch>
        batch_expand(const std::array<
            SmallElement<StandardResidue, Bounds<0, rns_bound>>, batch>& values)
        {
            static_assert(batch > 0);
            std::array<StandardResidue, batch> m2;
            BLST_UNROLL
            for (int index = 0; index < batch; ++index) {
                m2[index] = values[index].element;
            }
            const auto m1 = expand_m2_batch(m2);
            std::array<
                ExpandedElement<StandardResidue, Bounds<0, rns_bound>>, batch> result;
            BLST_UNROLL
            for (int index = 0; index < batch; ++index) {
                result[index] = ExpandedElement<
                    StandardResidue, Bounds<0, rns_bound>>(m1[index], m2[index]);
            }
            return result;
        }

        template<int batch>
        static INLINE std::array<StandardResidue, batch> batch_expand(
            const std::array<StandardResidue, batch>& m2)
        {
            static_assert(batch > 0);
            return expand_m2_batch(m2);
        }

        template<class ElementType, class RNSBoundsType>
        static INLINE auto small_m2(
            const ExpandedElement<ElementType, RNSBoundsType>& element)
        {
            return SmallElement<ElementType, RNSBoundsType>(element.m2);
        }

        template<class ElementType, class RNSBoundsType>
        static INLINE auto prep_expand(
            const SmallElement<ElementType, RNSBoundsType>& element)
        {
            if constexpr (RNSBoundsType::lower < 0) {
                constexpr int offset_shift =
                    ceil_log2(-RNSBoundsType::lower);
                const auto modulus = SmallElement<
                    FullyReducedResidue,
                    Bounds<1, 1>>(
                        FullyReducedResidue(
                            AVXVector<limbs>(generated::modulus_in_rns_m2)));
                return prep_expand(
                    element + modulus.template shift_left<offset_shift>());
            }
            else {
                return SmallElement<StandardResidue, RNSBoundsType>(
                    detail::ReducerM2::template reduce_auto<STANDARD_BOUND>(
                        element.element));
            }
        }

    private:
        template<class... Args, std::size_t... indices>
        static INLINE auto make_expand_array(
            const std::tuple<Args...>& args, std::index_sequence<indices...>)
        {
            return std::array<StandardResidue, sizeof...(Args)>{
                std::get<indices>(args).element...};
        }

        template<class Small, class Expanded>
        static INLINE auto reattach_expand_bounds(
            const Small& small, const Expanded& expanded)
        {
            return ExpandedElement<
                StandardResidue,
                std::remove_cv_t<decltype(small.rns_bounds)>>(expanded, small);
        }

        template<class... Args, std::size_t... indices>
        static INLINE auto reattach_expand_tuple(
            const std::tuple<Args...>& args,
            const std::array<StandardResidue, sizeof...(Args)>& expanded,
            std::index_sequence<indices...>)
        {
            return std::tuple<decltype(reattach_expand_bounds(
                std::get<indices>(args), expanded[indices]))...>{
                reattach_expand_bounds(std::get<indices>(args), expanded[indices])...};
        }

    public:
        template<class... Args>
        static INLINE auto batch_expand(const Args&... args)
        {
            constexpr int batch = static_cast<int>(sizeof...(Args));
            static_assert(batch > 0);
            const auto tuple = std::tie(args...);
            const auto m2 = make_expand_array(
                tuple, std::make_index_sequence<sizeof...(Args)>{});
            const auto m1 = batch_expand<batch>(m2);
            return reattach_expand_tuple(
                tuple, m1, std::make_index_sequence<sizeof...(Args)>{});
        }

        template<int batch, int64_t A = MAX_ADD>
        static INLINE auto batch_reduce(
            const std::array<ReadyElement<A>, batch>& values)
        {
            static_assert(batch > 0);
            std::array<AVXVector<limbs>, batch> residues;
            std::array<AVXVector<limbs>, batch> high;
            std::array<AVXVector<limbs>, batch> low;
            BLST_UNROLL
            for (int index = 0; index < batch; ++index) {
                residues[index] = values[index].input.data;
                high[index] = values[index].acc.high.data;
                low[index] = values[index].acc.low.data;
            }
            r1_change_base_raw_batch(residues, high, low);

            constexpr int64_t change_base_bound =
                2 * limbs * STANDARD_BOUND;
            constexpr int64_t combined_bound = A + change_base_bound;
            using CombinedBounds = Bounds<0, combined_bound>;
            using CombinedResidue =
                BoundedElement<CombinedBounds, limbs, element_bits>;
            using ReducedResidue = decltype(
                detail::ReducerM2::reduce_wide(
                    WideElement<CombinedBounds, limbs, element_bits>{}));
            using Result = SmallElement<ReducedResidue, Bounds<0, 1>>;
            std::array<Result, batch> result;
            BLST_UNROLL
            for (int index = 0; index < batch; ++index) {
                const WideElement<CombinedBounds, limbs, element_bits> combined{
                    CombinedResidue(high[index]), CombinedResidue(low[index])};
                result[index] = Result(detail::ReducerM2::reduce_wide(combined));
            }
            return result;
        }

        template<int batch>
        static INLINE std::array<StandardElement, batch> batch_reduce_expand(
            const std::array<ReadyElement<MAX_ADD>, batch>& values)
        {
            const auto reduced = batch_reduce<batch>(values);
            std::array<SmallStandardElement, batch> m2;
            BLST_UNROLL
            for (int index = 0; index < batch; ++index) {
                m2[index] = SmallStandardElement(
                    detail::ReducerM2::template reduce_auto<STANDARD_BOUND>(
                        reduced[index].element));
            }
            return batch_expand<batch>(m2);
        }

    private:
        template<int64_t A, class... Args, std::size_t... indices>
        static INLINE auto make_ready_array(
            const std::tuple<Args...>& args, std::index_sequence<indices...>)
        {
            return std::array<ReadyElement<A>, sizeof...(Args)>{
                ready<A>(std::get<indices>(args))...};
        }

    public:
        template<class... Args>
        static INLINE auto batch_reduce_expand(const Args&... args)
        {
            constexpr int batch = static_cast<int>(sizeof...(Args));
            static_assert(batch > 0);
            const auto tuple = std::tie(args...);
            const auto values = make_ready_array<MAX_ADD>(
                tuple, std::make_index_sequence<sizeof...(Args)>{});
            return batch_reduce_expand<batch>(values);
        }

        template<class... Args>
        static INLINE auto batch_reduce(const Args&... args)
        {
            constexpr int batch = static_cast<int>(sizeof...(Args));
            static_assert(batch > 0);
            const auto tuple = std::tie(args...);
            constexpr int64_t product_bound = MULT_OK_BOUND * MULT_OK_BOUND;
            const auto values = make_ready_array<product_bound>(
                tuple, std::make_index_sequence<sizeof...(Args)>{});
            return batch_reduce<batch>(values);
        }

        static INLINE StandardElement modmul(
            const StandardElement& left, const StandardElement& right)
        {
            return batch_reduce_expand(left * right)[0];
        }

        template<int batch>
        static INLINE std::array<StandardElement, batch> batch_modmul(
            const std::array<StandardElement, batch>& left,
            const std::array<StandardElement, batch>& right)
        {
            std::array<ReadyElement<MAX_ADD>, batch> products;
            BLST_UNROLL
            for (int index = 0; index < batch; ++index) {
                products[index] = ready<MAX_ADD>(left[index] * right[index]);
            }
            return batch_reduce_expand<batch>(products);
        }

        template<class Value>
        static INLINE auto mul_3b(const Value& value)
        {
            return value * std::integral_constant<int, 12>{};
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
            return Element(m1, m2);
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
                result[index] = Element(m1[index], m2[index]);
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

        static constexpr const std::array<uint64_t, limbs>&
        fp12_offset_m1_high() noexcept
        {
            return generated::fp12_offset_m1_high;
        }

        static constexpr const std::array<uint64_t, limbs>&
        fp12_offset_m1_low() noexcept
        {
            return generated::fp12_offset_m1_low;
        }

        static constexpr const std::array<uint64_t, limbs>&
        fp12_offset_m2_high() noexcept
        {
            return generated::fp12_offset_m2_high;
        }

        static constexpr const std::array<uint64_t, limbs>&
        fp12_offset_m2_low() noexcept
        {
            return generated::fp12_offset_m2_low;
        }

    private:
        friend class Fp2RingCandidate;

    };

    static_assert(std::is_empty_v<FpRingCandidate>);
    static_assert(std::is_trivially_default_constructible_v<FpRingCandidate>);
    static_assert(std::is_trivially_copyable_v<FpRingCandidate>);
} // namespace blst::avx::matrixnok
