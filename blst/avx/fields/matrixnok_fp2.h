// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "../rns/matrixnok_fp.h"

#include <array>
#include <cstddef>
#include <type_traits>

namespace blst::avx::matrixnok
{
    /**
     * Isolated BLS12-381 Fp2 candidate over the MatrixNoK base field.
     *
     * Elements represent `c0 + u*c1` with `u^2 = -1`. Inputs and outputs use
     * the owning base ring's standard two-basis representation. Operations are
     * allocation-free, nonthrowing, and constant-schedule with respect to field
     * values. Exact output aliasing with either input is supported by the
     * output-reference wrappers; partial overlap is not.
     */
    class Fp2RingCandidate
    {
    public:
        using Fp = FpRingCandidate;

        struct Element
        {
            Fp::Element c0;
            Fp::Element c1;
        };

        constexpr explicit Fp2RingCandidate(const Fp& ring) noexcept
            : ring_(&ring)
        {
        }

        INLINE Element multiply(const Element& left, const Element& right) const
        {
            // (a + bu)(c + du) = (ac - bd) + (ad + bc)u.
            const Fp::Element minus_b = ring_->negate_for_product(left.c1);
            const AVXVector<Fp::limbs> zero(0);
            std::array<AVXVector<Fp::limbs>, 2> m1_high;
            std::array<AVXVector<Fp::limbs>, 2> m1_low;
            std::array<AVXVector<Fp::limbs>, 2> m2_high;
            std::array<AVXVector<Fp::limbs>, 2> m2_low;

            m1_high[0] = zero.mulhi(left.c0.m1.data, right.c0.m1.data)
                .mulhi(minus_b.m1.data, right.c1.m1.data);
            m1_low[0] = zero.mullo(left.c0.m1.data, right.c0.m1.data)
                .mullo(minus_b.m1.data, right.c1.m1.data);
            m2_high[0] = zero.mulhi(left.c0.m2.data, right.c0.m2.data)
                .mulhi(minus_b.m2.data, right.c1.m2.data);
            m2_low[0] = zero.mullo(left.c0.m2.data, right.c0.m2.data)
                .mullo(minus_b.m2.data, right.c1.m2.data);

            m1_high[1] = zero.mulhi(left.c0.m1.data, right.c1.m1.data)
                .mulhi(left.c1.m1.data, right.c0.m1.data);
            m1_low[1] = zero.mullo(left.c0.m1.data, right.c1.m1.data)
                .mullo(left.c1.m1.data, right.c0.m1.data);
            m2_high[1] = zero.mulhi(left.c0.m2.data, right.c1.m2.data)
                .mulhi(left.c1.m2.data, right.c0.m2.data);
            m2_low[1] = zero.mullo(left.c0.m2.data, right.c1.m2.data)
                .mullo(left.c1.m2.data, right.c0.m2.data);

            constexpr int product_bound = 2 * Fp::standard_bound * Fp::standard_bound;
            static_assert(product_bound == 8);
            const auto result = ring_->template finish_wide_products<product_bound>(
                m1_high, m1_low, m2_high, m2_low);
            return {result[0], result[1]};
        }

        INLINE void multiply(Element& out, const Element& left, const Element& right) const
        {
            out = multiply(left, right);
        }

        INLINE Element square(const Element& value) const
        {
            // (a + bu)^2 = (a^2 - b^2) + 2ab*u. The imaginary
            // accumulator doubles one product rather than issuing it twice.
            const Fp::Element minus_b = ring_->negate_for_product(value.c1);
            const AVXVector<Fp::limbs> zero(0);
            std::array<AVXVector<Fp::limbs>, 2> m1_high;
            std::array<AVXVector<Fp::limbs>, 2> m1_low;
            std::array<AVXVector<Fp::limbs>, 2> m2_high;
            std::array<AVXVector<Fp::limbs>, 2> m2_low;

            m1_high[0] = zero.mulhi(value.c0.m1.data, value.c0.m1.data)
                .mulhi(minus_b.m1.data, value.c1.m1.data);
            m1_low[0] = zero.mullo(value.c0.m1.data, value.c0.m1.data)
                .mullo(minus_b.m1.data, value.c1.m1.data);
            m2_high[0] = zero.mulhi(value.c0.m2.data, value.c0.m2.data)
                .mulhi(minus_b.m2.data, value.c1.m2.data);
            m2_low[0] = zero.mullo(value.c0.m2.data, value.c0.m2.data)
                .mullo(minus_b.m2.data, value.c1.m2.data);

            const AVXVector<Fp::limbs> m1_ab_high =
                zero.mulhi(value.c0.m1.data, value.c1.m1.data);
            const AVXVector<Fp::limbs> m1_ab_low =
                zero.mullo(value.c0.m1.data, value.c1.m1.data);
            const AVXVector<Fp::limbs> m2_ab_high =
                zero.mulhi(value.c0.m2.data, value.c1.m2.data);
            const AVXVector<Fp::limbs> m2_ab_low =
                zero.mullo(value.c0.m2.data, value.c1.m2.data);
            m1_high[1] = m1_ab_high.add(m1_ab_high);
            m1_low[1] = m1_ab_low.add(m1_ab_low);
            m2_high[1] = m2_ab_high.add(m2_ab_high);
            m2_low[1] = m2_ab_low.add(m2_ab_low);

            constexpr int product_bound = 2 * Fp::standard_bound * Fp::standard_bound;
            static_assert(product_bound == 8);
            const auto result = ring_->template finish_wide_products<product_bound>(
                m1_high, m1_low, m2_high, m2_low);
            return {result[0], result[1]};
        }

        INLINE void square(Element& out, const Element& value) const
        {
            out = square(value);
        }

    private:
        const Fp* ring_;
    };

    static_assert(sizeof(Fp2RingCandidate) == sizeof(const FpRingCandidate*));
    static_assert(std::is_trivially_copyable_v<Fp2RingCandidate>);
    static_assert(std::is_trivially_destructible_v<Fp2RingCandidate>);
}
