// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "fp12.h"
#include "../rns/matrixnok_fp.h"

#include <algorithm>
#include <type_traits>

namespace blst::avx::matrixnok
{
    // Parameter-bound Fp12 constants for the selected MatrixNoK field tower.
    // This policy owns no state: every accessor materializes directly from
    // aligned, independently generated lane constants. Coefficient ordering is
    // inherited unchanged from FP12Ring and documented in the construction map.
    class Fp12ConstantsCandidate
    {
        using Ring = FpRingCandidate;
        using StandardElement = Ring::StandardElement;
        using StandardResidue = Ring::StandardResidue;

        static constexpr int offset_rns_bound =
            10 * Ring::STANDARD_BOUND * Ring::STANDARD_BOUND * 6 * 6;
        static constexpr int offset_ele_bound =
            10 * Ring::STANDARD_BOUND *
            std::max(Ring::STANDARD_BOUND, Ring::MULT_OK_BOUND);
        using OffsetElement = BoundedElement<
            Bounds<offset_ele_bound, offset_ele_bound + 1>,
            Ring::limbs, Ring::element_bits>;
        using OffsetWide = WideElement<
            Bounds<offset_ele_bound, offset_ele_bound + 1>,
            Ring::limbs, Ring::element_bits>;
        using Offset = Ring::WideOffsetElement<
            offset_ele_bound, offset_rns_bound>;

        static INLINE StandardElement standard(
            const std::array<uint64_t, Ring::limbs>& m1,
            const std::array<uint64_t, Ring::limbs>& m2)
        {
            return StandardElement(
                StandardResidue(AVXVector<Ring::limbs>(m1)),
                StandardResidue(AVXVector<Ring::limbs>(m2)));
        }

    public:
        static constexpr bool split_cyclotomic_v3 = true;

        Fp12ConstantsCandidate() noexcept = default;
        explicit Fp12ConstantsCandidate(const Ring&) noexcept {}

        static INLINE StandardElement one_third()
        {
            return standard(generated::fp12_one_third_m1,
                generated::fp12_one_third_m2);
        }

        static INLINE StandardElement neg_two_thirds()
        {
            return standard(generated::fp12_neg_two_thirds_m1,
                generated::fp12_neg_two_thirds_m2);
        }

        static INLINE StandardElement frob1_2()
        {
            return standard(generated::fp12_frob1_2_m1,
                generated::fp12_frob1_2_m2);
        }

        static INLINE StandardElement frob1_3()
        {
            return standard(generated::fp12_frob1_3_m1,
                generated::fp12_frob1_3_m2);
        }

        static INLINE StandardElement frob1_4()
        {
            return standard(generated::fp12_frob1_4_m1,
                generated::fp12_frob1_4_m2);
        }

        static INLINE StandardElement frob1_5()
        {
            return standard(generated::fp12_frob1_5_m1,
                generated::fp12_frob1_5_m2);
        }

        static INLINE StandardElement frob1_6()
        {
            return standard(generated::fp12_frob1_6_m1,
                generated::fp12_frob1_6_m2);
        }

        static INLINE StandardElement frob1_7()
        {
            return standard(generated::fp12_frob1_7_m1,
                generated::fp12_frob1_7_m2);
        }

        static INLINE StandardElement frob1_8()
        {
            return standard(generated::fp12_frob1_8_m1,
                generated::fp12_frob1_8_m2);
        }

        static INLINE StandardElement frob1_9()
        {
            return standard(generated::fp12_frob1_9_m1,
                generated::fp12_frob1_9_m2);
        }

        static INLINE StandardElement frob1_10()
        {
            return standard(generated::fp12_frob1_10_m1,
                generated::fp12_frob1_10_m2);
        }

        static INLINE StandardElement frob2_1() { return frob1_4(); }
        static INLINE StandardElement frob2_2() { return frob1_2(); }

        static INLINE StandardElement frob2_3()
        {
            return standard(generated::fp12_frob2_3_m1,
                generated::fp12_frob2_3_m2);
        }

        static INLINE StandardElement frob2_5() { return frob1_3(); }
        static INLINE StandardElement frob3_2() { return frob1_8(); }
        static INLINE StandardElement frob3_3() { return frob1_7(); }

        static INLINE Offset offset()
        {
            return Offset(
                OffsetWide(
                    OffsetElement(AVXVector<Ring::limbs>(
                        generated::fp12_offset_m1_high)),
                    OffsetElement(AVXVector<Ring::limbs>(
                        generated::fp12_offset_m1_low))),
                OffsetWide(
                    OffsetElement(AVXVector<Ring::limbs>(
                        generated::fp12_offset_m2_high)),
                    OffsetElement(AVXVector<Ring::limbs>(
                        generated::fp12_offset_m2_low))));
        }
    };

    using Fp12RingCandidate =
        FP12Ring<FpRingCandidate, Fp12ConstantsCandidate, 3>;

    static_assert(
        10 * FpRingCandidate::STANDARD_BOUND *
        FpRingCandidate::STANDARD_BOUND * 6 * 6 == 1440);
    static_assert(
        10 * FpRingCandidate::STANDARD_BOUND *
        std::max(FpRingCandidate::STANDARD_BOUND,
            FpRingCandidate::MULT_OK_BOUND) == 80);
    static_assert(std::is_empty_v<Fp12ConstantsCandidate>);
    static_assert(std::is_trivially_default_constructible_v<
        Fp12ConstantsCandidate>);
    static_assert(std::is_trivially_copyable_v<Fp12ConstantsCandidate>);
    static_assert(std::is_empty_v<Fp12RingCandidate>);
    static_assert(std::is_trivially_default_constructible_v<
        Fp12RingCandidate>);

    using Fp12Ring = Fp12RingCandidate;
} // namespace blst::avx::matrixnok
