// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "../cpu/vector/multiplication.hpp"
#include "../generated_matrixnok_constants.hpp"
#include "elementwise_reduce.h"
#include "matrixnok_change_base.h"
#include "ring_element.h"

#include <array>

namespace blst::avx::matrixnok
{
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
        using System = IntRNS2System<381, limbs, MontgomeryReduce, 50>;

        static System make_system()
        {
            return System(
                IntRNS2<MontgomeryReduce, limbs, limbs>(
                    generated::moduli1,
                    generated::moduli2,
                    generated::r1_rns_mat,
                    generated::r1_correction,
                    generated::r1_shifted_quotient_estimations,
                    generated::r2_rns_mat,
                    generated::r2_correction,
                    generated::r2_shifted_quotient_estimations),
                ConvertToRNS<381, limbs, 50>(
                    generated::convert_to_rns_mat,
                    generated::convert_to_correction,
                    generated::convert_to_shifted_quotient_estimations,
                    generated::moduli2),
                ConvertFromRNS<381, limbs>(
                    generated::convert_from_rns_mat,
                    generated::convert_from_correction,
                    generated::convert_from_shifted_quotient_estimations,
                    std::array<uint64_t, 1>{1}),
                generated::moduli1,
                generated::moduli2);
        }

        INLINE StandardResidue expand_m2(const StandardResidue& m2) const
        {
            AVXVector<limbs> high(0);
            AVXVector<limbs> low(0);
            system_.intrns2.r2.rns_reduce_raw_shifted(m2.data, high, low);

            using WideBound = Bounds<0, 2 * limbs * standard_bound>;
            using Bounded = BoundedElement<WideBound, limbs, element_bits>;
            const WideElement<WideBound, limbs, element_bits> wide{
                Bounded(high), Bounded(low)};
            return reducer_m1_.template reduce_auto<standard_bound>(
                reducer_m1_.reduce_wide(wide));
        }

    public:
        FpRingCandidate()
            : system_(make_system())
            , reducer_m1_(system_.moduli1)
            , reducer_m2_(system_.moduli2)
        {
        }

        INLINE Element from_canonical_50bit_digits(const std::array<uint64_t, limbs>& digits) const
        {
            AVXVector<limbs> high(0);
            AVXVector<limbs> low(0);
            system_.convert_to.conversion_matrix.rns_reduce_raw_shifted(
                AVXVector<limbs>(digits), high, low);

            using ConvertBound = Bounds<0, 800>;
            using ConvertElement = BoundedElement<ConvertBound, limbs, element_bits>;
            const WideElement<ConvertBound, limbs, element_bits> converted{
                ConvertElement(high), ConvertElement(low)};
            const StandardResidue m2 = reducer_m2_.template reduce_auto<standard_bound>(
                reducer_m2_.reduce_wide(converted));
            return {expand_m2(m2), m2};
        }

        // Multiply modulo p. Exact output aliasing with either input is
        // supported; partial overlap is not.
        INLINE void multiply(Element& out, const Element& left, const Element& right) const
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
            const StandardResidue m1_reduced = reducer_m1_.template reduce_auto<standard_bound>(
                reducer_m1_.reduce_wide(m1_product));

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
            const StandardResidue m2 = reducer_m2_.template reduce_auto<standard_bound>(
                reducer_m2_.reduce_wide(combined));
            const StandardResidue m1 = expand_m2(m2);
            out = {m1, m2};
        }

        INLINE void convert_from_raw(
            const Element& value,
            AVXVector<8>& out_high,
            AVXVector<8>& out_low) const
        {
            out_high = AVXVector<8>(0);
            out_low = AVXVector<8>(0);
            system_.convert_from.conversion_matrix.rns_reduce_raw(
                value.m2.data, out_high, out_low);
        }

        INLINE const std::array<uint64_t, limbs>& moduli1() const
        {
            return system_.moduli1;
        }

        INLINE const std::array<uint64_t, limbs>& moduli2() const
        {
            return system_.moduli2;
        }

    private:
        const System system_;
        const MontgomeryReduce2<limbs, element_bits, 6> reducer_m1_;
        const MontgomeryReduce2<limbs, element_bits, 6> reducer_m2_;
    };
} // namespace blst::avx::matrixnok
