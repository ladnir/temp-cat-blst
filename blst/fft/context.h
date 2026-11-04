// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// CatBLST scalar-field and group coefficient contexts for the generic FFT.

#pragma once

#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fields/Fr.h"

#include <cstdint>
#include <type_traits>

namespace blst
{
    /**
     * Arithmetic context for an FFT whose values and roots are both Fr.
     *
     * Inputs and outputs use reduced Montgomery Fr representation. Root and
     * size inversions require a nonzero operand. These operations allocate no
     * storage and do not throw. FFT schedules depend only on public shapes;
     * this context does not inspect coefficient values.
     */
    struct FrCoeffCtx
    {
        [[nodiscard]] static Fr plus(const Fr& a, const Fr& b)
        {
            Fr result;
            Fr::add(result, a, b);
            return result;
        }

        [[nodiscard]] static Fr minus(const Fr& a, const Fr& b)
        {
            Fr result;
            Fr::sub(result, a, b);
            return result;
        }

        [[nodiscard]] static Fr mult(const Fr& a, const Fr& b)
        {
            Fr result;
            Fr::mul(result, a, b);
            return result;
        }

        /** @pre @p a is nonzero. */
        [[nodiscard]] static Fr inverse(const Fr& a)
        {
            Fr result;
            Fr::reciprocal(result, a);
            return result;
        }

        template<typename F>
        [[nodiscard]] static F from_u64(const std::uint64_t value)
        {
            static_assert(
                std::is_same_v<F, Fr>,
                "FrCoeffCtx constructs only Fr coefficients");
            return Fr::from_uint64(value);
        }
    };

    /**
     * Arithmetic context for FFTs over G1, G2, or GT with public Fr roots.
     *
     * G1 and G2 operands must already be on-curve prime-order subgroup points;
     * GT operands must already satisfy GT membership. Release builds perform
     * no membership validation. Addition, subtraction, and negation preserve
     * those domains. Multiplication applies a public Fr root with
     * `mult_vartime`; it must not be used when the Fr operand is secret.
     *
     * The G1-by-G2 overload computes the pairing into GT. Its operands have the
     * same subgroup preconditions. All operations return a new value, allocate
     * no storage, retain no references, and do not throw.
     */
    struct CoeffCtx
    {
    private:
        template<typename Group>
        static constexpr bool is_supported_group =
            std::is_same_v<Group, G1Point> ||
            std::is_same_v<Group, G2Point> ||
            std::is_same_v<Group, GT>;

    public:
        template<typename Group>
        [[nodiscard]] static Group plus(const Group& a, const Group& b)
        {
            static_assert(
                is_supported_group<Group>,
                "CoeffCtx supports only G1, G2, and GT values");
            Group result;
            Group::dadd(result, a, b);
            return result;
        }

        template<typename Group>
        [[nodiscard]] static Group minus(const Group& a, const Group& b)
        {
            static_assert(
                is_supported_group<Group>,
                "CoeffCtx supports only G1, G2, and GT values");
            Group negated;
            Group::neg(negated, b);
            Group result;
            Group::dadd(result, a, negated);
            return result;
        }

        template<typename Group>
        [[nodiscard]] static Group neg(const Group& value)
        {
            static_assert(
                is_supported_group<Group>,
                "CoeffCtx supports only G1, G2, and GT values");
            Group result;
            Group::neg(result, value);
            return result;
        }

        template<typename Group>
        [[nodiscard]] static Group mult(const Group& value, const Fr& scalar)
        {
            static_assert(
                is_supported_group<Group>,
                "CoeffCtx supports only G1, G2, and GT values");
            Fr scalar_non_montgomery;
            Fr::from_mont(scalar_non_montgomery, scalar);
            Group result;
            Group::mult_vartime(
                result,
                value,
                reinterpret_cast<const byte*>(&scalar_non_montgomery),
                256);
            return result;
        }

        [[nodiscard]] static GT mult(const G1Point& p, const G2Point& q)
        {
            G1AffinePoint p_affine;
            G1Point::to_affine(p_affine, p);
            G2AffinePoint q_affine;
            G2Point::to_affine(q_affine, q);

            MillerResult miller;
            GT::miller_loop(miller, &q_affine, &p_affine);
            GT result;
            GT::final_exp(result, miller);
            return result;
        }

        /** @pre @p value is nonzero. */
        [[nodiscard]] static Fr inverse(const Fr& value)
        {
            Fr result;
            Fr::reciprocal(result, value);
            return result;
        }

        template<typename F>
        [[nodiscard]] static F from_u64(const std::uint64_t value)
        {
            static_assert(
                std::is_same_v<F, Fr>,
                "CoeffCtx constructs only Fr coefficients");
            return Fr::from_uint64(value);
        }
    };

    using G1CoeffCtx = CoeffCtx;
    using G2CoeffCtx = CoeffCtx;
    using GTCoeffCtx = CoeffCtx;

    namespace detail
    {
        // Unlike CoeffCtx, this works with the fixed scalar point types used
        // by the C ABI even when the selected C++ backend is AVX.
        template<class Group>
        struct GroupFftCoeffCtx
        {
            static Group plus(const Group& a, const Group& b)
            {
                Group result;
                Group::dadd(result, a, b);
                return result;
            }

            static Group minus(const Group& a, const Group& b)
            {
                Group negated;
                Group::neg(negated, b);
                Group result;
                Group::dadd(result, a, negated);
                return result;
            }

            static Group mult(const Group& a, const Fr& b)
            {
                Fr canonical;
                Fr::from_mont(canonical, b);
                Group result;
                Group::mult_vartime(
                    result, a, reinterpret_cast<const byte*>(&canonical), 256);
                return result;
            }

            static Fr inverse(const Fr& a)
            {
                Fr result;
                Fr::reciprocal(result, a);
                return result;
            }

            template<typename F>
            static Fr from_u64(const std::uint64_t a)
            {
                static_assert(std::is_same_v<F, Fr>);
                return Fr::from_uint64(a);
            }
        };
    }
}
