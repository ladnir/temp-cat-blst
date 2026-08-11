// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// Fixed BLS12-381 MatrixNoK construction adapted from the algorithm and
// parameters in Simon Langowski's VROOM commit d6ecd534. This implementation
// is locally composed around CatBLST's field, bound, and proof interfaces.
#pragma once

#include "../generated_matrixnok_constants.hpp"
#include "change_base.h"
#include "matrixnok_change_base.h"
#include "matrixnok_fixed.h"
#include "ring_element.h"
#include "../../fields/Fp.h"

#include <array>
#include <tuple>
#include <type_traits>
#include <utility>

namespace blst::avx::matrixnok
{
    // A field value carries the same element of Fp in two fixed 8-lane RNS
    // bases. Standard residues are in [0, 2*m_i) lane by lane. The separate
    // rns_bounds metadata tracks the represented multiple of p across the two
    // bases; it is a compile-time proof claim, not runtime validation.
    // `is_zero()` and `equal()` on the carrier compare redundant RNS storage,
    // not general field semantics. Internal formulas may use them only where
    // their exact-representation invariant is established. Public observation
    // boundaries canonicalize through ConvertFrom as required by AUD-228.
    class FpRingCandidate
    {
    public:
        static constexpr int limbs = 8;
        static constexpr int element_bits = 50;
        static constexpr int standard_bound = 2;
        static constexpr int64_t max_add = 800;
        // VROOM's BLS12-381 bound ledger admits represented multiples in
        // [-1932, 2377]. Negative-range completion adds 1932*p, translating
        // that interval to [0, 4309] without changing the field value.
        static constexpr int64_t min_tracked_rns_multiple = -1932;
        static constexpr int64_t max_tracked_rns_multiple = 2377;
        static constexpr int64_t max_negative_rns =
            -min_tracked_rns_multiple;
        static constexpr int64_t max_rns_multiple =
            max_tracked_rns_multiple - min_tracked_rns_multiple;
        static constexpr int64_t base_redundant_multiple = 40;
        static constexpr int modulus_bits = 381;
        static constexpr int rns_bits = element_bits * limbs - 1;
        static constexpr int64_t wide_redundant_multiple =
            base_redundant_multiple * base_redundant_multiple;
        static constexpr int combined_rns_multiple_bits = static_cast<int>(
            blst::log2ceil(static_cast<uint64_t>(
                max_rns_multiple * wide_redundant_multiple)));
        static constexpr int multiply_product_bound =
            standard_bound * standard_bound;
        static constexpr int r1_accumulation_bound =
            limbs * standard_bound;
        // The generic delayed-reduction path deliberately retains its selected
        // conservative type bound. Narrowing it to the exact R1 contribution
        // changes template-selected tower codegen and requires a separate
        // benchmark gate; Lean still proves the stronger bound above.
        static constexpr int delayed_r1_type_bound =
            2 * r1_accumulation_bound;
        static constexpr int multiply_combined_bound =
            multiply_product_bound + r1_accumulation_bound;

        static_assert(multiply_product_bound == 4);
        static_assert(r1_accumulation_bound == 16);
        static_assert(delayed_r1_type_bound == 32);
        static_assert(multiply_combined_bound == 20);
        static_assert(max_negative_rns == 1932);
        static_assert(max_rns_multiple == 4309);
        static_assert(combined_rns_multiple_bits == 23);
        static_assert(
            combined_rns_multiple_bits + 2 * modulus_bits <= 2 * rns_bits,
            "tracked MatrixNoK range must fit in the combined RNS bases");
        static_assert(multiply_combined_bound < base_redundant_multiple);

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
        static BLST_FORCEINLINE StandardResidue expand_m2(const StandardResidue& m2)
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
        static BLST_FORCEINLINE std::array<StandardResidue, batch> expand_m2_batch(
            const std::array<StandardResidue, batch>& m2)
        {
            static_assert(batch > 0);
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

    public:
        constexpr FpRingCandidate() noexcept = default;

        static BLST_FORCEINLINE Element from_canonical_limbs(const blst::Fp& canonical)
        {
            // Narrow contract: canonical is an ordinary representative in
            // [0, p), not Montgomery form. The assertion has no Release cost.
            BLST_ASSERT(blst::Fp::is_redc(canonical));
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

        static BLST_FORCEINLINE Element from_scalar(const blst::Fp& montgomery)
        {
            blst::Fp canonical = {};
            blst::Fp::from_mont(canonical, montgomery);
            return from_canonical_limbs(canonical);
        }

        static BLST_FORCEINLINE StandardElement zero()
        {
            return StandardElement{};
        }

        static BLST_FORCEINLINE StandardElement one()
        {
            return StandardElement(
                StandardResidue(AVXVector<limbs>(generated::one_m1)),
                StandardResidue(AVXVector<limbs>(generated::one_m2)));
        }

        template<int multiple = 1>
        static BLST_FORCEINLINE auto elementwise_moduli()
        {
            static_assert(multiple >= 1);
            constexpr int index = static_cast<int>(blst::log2ceil(multiple));
            auto m1 = detail::ReducerM1::template moduli<index>();
            auto m2 = detail::ReducerM2::template moduli<index>();
            return ExpandedElement<decltype(m1), Bounds<0, 0>>(m1, m2);
        }

        template<int64_t lower, int64_t upper, class RNSBoundsType>
        static BLST_FORCEINLINE auto wideify(
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

        static BLST_FORCEINLINE auto negative_range_offset()
        {
            const auto offset = OffsetElement<max_negative_rns>(
                FullyReducedResidue(AVXVector<limbs>(generated::max_neg_offset_m1)),
                FullyReducedResidue(AVXVector<limbs>(generated::max_neg_offset_m2)));
            return wideify(offset);
        }

        template<int64_t lower, int64_t upper, class RNSBoundsType,
            int to = MULT_OK_BOUND>
        static BLST_FORCEINLINE auto prep(const ExpandedElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            RNSBoundsType>& element)
        {
            auto m1 = detail::ReducerM1::template reduce_auto<to>(element.m1);
            auto m2 = detail::ReducerM2::template reduce_auto<to>(element.m2);
            return ExpandedElement<decltype(m1), RNSBoundsType>(m1, m2);
        }

        template<class ElementType>
        static BLST_FORCEINLINE auto prep_left(const ElementType& element)
        {
            return prep(element);
        }

        template<int64_t lower, int64_t upper, class RNSBoundsType>
        static BLST_FORCEINLINE auto prep_standard(const ExpandedElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            RNSBoundsType>& element)
        {
            return prep<lower, upper, RNSBoundsType, STANDARD_BOUND>(element);
        }

        template<int64_t lower, int64_t upper, class RNSBoundsType,
            int to = MULT_OK_BOUND>
        static BLST_FORCEINLINE auto prep(const SmallElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            RNSBoundsType>& element)
        {
            auto m2 = detail::ReducerM2::template reduce_auto<to>(element.element);
            return SmallElement<decltype(m2), RNSBoundsType>(m2);
        }

        template<int64_t lower, int64_t upper, class RNSBoundsType>
        static BLST_FORCEINLINE auto prep_standard(const SmallElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            RNSBoundsType>& element)
        {
            return prep<lower, upper, RNSBoundsType, STANDARD_BOUND>(element);
        }

        template<int64_t lower, int64_t upper, class RNSBoundsType>
        static BLST_FORCEINLINE auto negate(const ExpandedElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            RNSBoundsType>& value)
        {
            return elementwise_moduli<upper>() - value;
        }

        template<int64_t lower, int64_t upper, class RNSBoundsType>
        static BLST_FORCEINLINE auto negate(const SmallElement<
            BoundedElement<Bounds<lower, upper>, limbs, element_bits>,
            RNSBoundsType>& value)
        {
            const auto moduli = elementwise_moduli<upper>();
            const SmallElement<decltype(moduli.m2), Bounds<0, 0>> offset(
                moduli.m2);
            return offset - value;
        }

        template<int64_t lower, int64_t upper>
        static BLST_FORCEINLINE StandardElement standard_negate(const ExpandedElement<
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
        static BLST_FORCEINLINE auto standard_negate(const ExpandedElement<
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
        static BLST_FORCEINLINE auto complete(const ExpandedElement<
            ElementType, Bounds<rns_lower, rns_upper>>& value)
        {
            if constexpr (rns_lower < 0) {
                return negative_range_offset() + value;
            } else {
                return value.complete();
            }
        }

        template<int64_t A, class ElementType, class RNSBoundsType>
        static BLST_FORCEINLINE ReadyElement<A> ready(
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

        static BLST_FORCEINLINE Element from_canonical_50bit_digits(
            const std::array<uint64_t, limbs>& digits)
        {
            // Narrow contract: each lane is a radix-2^50 digit and the
            // reconstructed ordinary integer is in [0, p). The checks below
            // are Debug evidence only; production constants and conversions
            // establish this precondition before entering the hot path.
#ifndef NDEBUG
            constexpr uint64_t mask = (uint64_t{1} << element_bits) - 1;
            blst::Fp canonical = {};
            for (std::size_t digit = 0; digit < digits.size(); ++digit) {
                BLST_ASSERT((digits[digit] & ~mask) == 0);
                const std::size_t bit = digit * element_bits;
                const std::size_t limb = bit / 64;
                const unsigned shift = static_cast<unsigned>(bit % 64);
                if (limb < blst::Fp::num_limb()) {
                    canonical.l[limb] |= digits[digit] << shift;
                }
                if (shift != 0) {
                    const uint64_t upper = digits[digit] >> (64 - shift);
                    if (limb + 1 < blst::Fp::num_limb()) {
                        canonical.l[limb + 1] |= upper;
                    }
                    else {
                        BLST_ASSERT(upper == 0);
                    }
                }
            }
            BLST_ASSERT(blst::Fp::is_redc(canonical));
#endif
            AVXVector<limbs> high(0);
            AVXVector<limbs> low(0);
            detail::ConvertToMatrix::reduce_raw_shifted(
                AVXVector<limbs>(digits), high, low);

            // ConvertTo maps canonical radix-2^50 digits directly into M2.
            // M1 is then derived from that reduced M2 value so both bases are
            // linked by one explicit conversion path.
            using ConvertBound = Bounds<0, max_add>;
            using ConvertElement = BoundedElement<ConvertBound, limbs, element_bits>;
            const WideElement<ConvertBound, limbs, element_bits> converted{
                ConvertElement(high), ConvertElement(low)};
            const StandardResidue m2 = detail::ReducerM2::template reduce_auto<standard_bound>(
                detail::ReducerM2::reduce_wide(converted));
            return Element(expand_m2(m2), m2);
        }

        template<int batch, int64_t rns_bound>
        static BLST_FORCEINLINE std::array<
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
        static BLST_FORCEINLINE std::array<StandardResidue, batch> batch_expand(
            const std::array<StandardResidue, batch>& m2)
        {
            static_assert(batch > 0);
            return expand_m2_batch(m2);
        }

        template<class ElementType, class RNSBoundsType>
        static BLST_FORCEINLINE auto small_m2(
            const ExpandedElement<ElementType, RNSBoundsType>& element)
        {
            return SmallElement<ElementType, RNSBoundsType>(element.m2);
        }

        template<class ElementType, class RNSBoundsType>
        static BLST_FORCEINLINE auto prep_expand(
            const SmallElement<ElementType, RNSBoundsType>& element)
        {
            if constexpr (RNSBoundsType::lower < 0) {
                static_assert(
                    RNSBoundsType::lower >= min_tracked_rns_multiple,
                    "MatrixNoK negative RNS bound exceeds the reviewed range");
                constexpr int offset_shift = static_cast<int>(
                    blst::log2ceil(static_cast<uint64_t>(-RNSBoundsType::lower)));
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
        static BLST_FORCEINLINE auto make_expand_array(
            const std::tuple<Args...>& args, std::index_sequence<indices...>)
        {
            return std::array<StandardResidue, sizeof...(Args)>{
                std::get<indices>(args).element...};
        }

        template<class Small, class Expanded>
        static BLST_FORCEINLINE auto reattach_expand_bounds(
            const Small& small, const Expanded& expanded)
        {
            return ExpandedElement<
                StandardResidue,
                typename Small::rns_bounds_type>(expanded, small);
        }

        template<class... Args, std::size_t... indices>
        static BLST_FORCEINLINE auto reattach_expand_tuple(
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
        static BLST_FORCEINLINE auto batch_expand(const Args&... args)
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
        static BLST_FORCEINLINE auto batch_reduce(
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

            constexpr int64_t combined_bound = A + delayed_r1_type_bound;
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
        static BLST_FORCEINLINE std::array<StandardElement, batch> batch_reduce_expand(
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
        static BLST_FORCEINLINE auto make_ready_array(
            const std::tuple<Args...>& args, std::index_sequence<indices...>)
        {
            return std::array<ReadyElement<A>, sizeof...(Args)>{
                ready<A>(std::get<indices>(args))...};
        }

    public:
        template<class... Args>
        static BLST_FORCEINLINE auto batch_reduce_expand(const Args&... args)
        {
            constexpr int batch = static_cast<int>(sizeof...(Args));
            static_assert(batch > 0);
            const auto tuple = std::tie(args...);
            const auto values = make_ready_array<MAX_ADD>(
                tuple, std::make_index_sequence<sizeof...(Args)>{});
            return batch_reduce_expand<batch>(values);
        }

        template<class... Args>
        static BLST_FORCEINLINE auto batch_reduce(const Args&... args)
        {
            constexpr int batch = static_cast<int>(sizeof...(Args));
            static_assert(batch > 0);
            const auto tuple = std::tie(args...);
            constexpr int64_t product_bound = MULT_OK_BOUND * MULT_OK_BOUND;
            const auto values = make_ready_array<product_bound>(
                tuple, std::make_index_sequence<sizeof...(Args)>{});
            return batch_reduce<batch>(values);
        }

        static BLST_FORCEINLINE StandardElement modmul(
            const StandardElement& left, const StandardElement& right)
        {
            // Production one-shot multiplication uses the exact-bound,
            // register-latency MatrixNoK kernel. Delayed reduction remains for
            // formulas that intentionally combine multiple products.
            return multiply(left, right);
        }

        template<int batch>
        static BLST_FORCEINLINE std::array<StandardElement, batch> batch_modmul(
            const std::array<StandardElement, batch>& left,
            const std::array<StandardElement, batch>& right)
        {
            return multiply_batch(left, right);
        }

        template<class Value>
        static BLST_FORCEINLINE auto mul_3b(const Value& value)
        {
            return value * std::integral_constant<int, 12>{};
        }

        static BLST_FORCEINLINE Element multiply(const Element& left, const Element& right)
        {
            // Selected MatrixNoK multiplication pipeline:
            //   1. form products independently in M1 and M2;
            //   2. Montgomery-reduce M1;
            //   3. add correction-free R1(M1) to the existing M2 product;
            //   4. Montgomery-reduce the combined M2 value; and
            //   5. expand reduced M2 back to M1.
            // The batch implementation below deliberately spells out the same
            // stages to preserve its independent throughput-oriented codegen.
            const AVXVector<limbs> zero(0);
            const AVXVector<limbs> m1_high = zero.mulhi(left.m1.data, right.m1.data);
            const AVXVector<limbs> m1_low = zero.mullo(left.m1.data, right.m1.data);
            const AVXVector<limbs> m2_high = zero.mulhi(left.m2.data, right.m2.data);
            const AVXVector<limbs> m2_low = zero.mullo(left.m2.data, right.m2.data);

            using ProductBound = Bounds<0, multiply_product_bound>;
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

            using CombinedBound = Bounds<0, multiply_combined_bound>;
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
        static BLST_FORCEINLINE void multiply(
            Element& out, const Element& left, const Element& right)
        {
            out = multiply(left, right);
        }

        template<std::size_t batch>
        static BLST_FORCEINLINE std::array<Element, batch> multiply_batch(
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

                using ProductBound = Bounds<0, multiply_product_bound>;
                using ProductElement = BoundedElement<ProductBound, limbs, element_bits>;
                const WideElement<ProductBound, limbs, element_bits> m1_product{
                    ProductElement(m1_high), ProductElement(m1_low)};
                residues[index] = detail::ReducerM1::template reduce_auto<standard_bound>(
                    detail::ReducerM1::reduce_wide(m1_product)).data;
            }

            r1_change_base_raw_batch(residues, high, low);

            using CombinedBound = Bounds<0, multiply_combined_bound>;
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

        static BLST_FORCEINLINE void convert_from_raw(
            const Element& value,
            AVXVector<limbs>& out_high,
            AVXVector<limbs>& out_low)
        {
            // The three vector objects must be distinct. The outputs are the
            // raw high/low radix-2^52 accumulators consumed by
            // `from_raw_accumulators`, not independently reduced field or RNS
            // values. FixedRNSMatrix
            // asserts the disjointness contract in Debug.
            out_high = AVXVector<limbs>(0);
            out_low = AVXVector<limbs>(0);
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

    };

    static_assert(std::is_empty_v<FpRingCandidate>);
    static_assert(std::is_trivially_default_constructible_v<FpRingCandidate>);
    static_assert(std::is_trivially_copyable_v<FpRingCandidate>);

    // Production spelling. The Candidate name remains as the exact theorem
    // correspondence owner until the post-integration Lean audit is closed.
    using FpRing = FpRingCandidate;
} // namespace blst::avx::matrixnok
