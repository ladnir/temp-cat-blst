// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// Fixed BLS12-381 MatrixNoK construction adapted from the algorithm and
// parameters in Simon Langowski's VROOM commit d6ecd534. This implementation
// is locally composed around CatBLST's vector, bound, and proof interfaces.
#pragma once

#include "../common/constexpr_utils.h"
#include "../common/preprocess.h"
#include "../cpu/vector/vector_impl.hpp"
#include "../generated_matrixnok_constants.hpp"
#include "bounded_element.h"
#include "ring_element.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace blst::avx::matrixnok::detail
{
    // Fixed BLS12-381 MatrixNoK support. All stored data is constant-initialized
    // read-only storage derived from the independently validated generated
    // parameters. The operation classes are empty and introduce no constructor,
    // allocation, guard variable, or mutable process state.
    inline constexpr int limbs = 8;
    inline constexpr int element_bits = 50;
    // G1's complete formulas transiently create signed bounds well beyond the
    // base-field-only reduction gap. Match VROOM's EC-depth table while
    // retaining the same fixed, compile-time reducer representation.
    inline constexpr int log_modulus_multiples = 12;
    inline constexpr int multiply_bits = 52;
    inline constexpr int quotient_precision = 64;

    static_assert(limbs == generated::BLS12_381_RNS_DIM);
    static_assert(limbs == generated::BLS12_381_PADDED_DIM);
    static_assert(limbs == generated::BLS12_381_CONVERT_TO_DIGITS);
    static_assert(element_bits == generated::BLS12_381_CONVERT_TO_RADIX_BITS);
    static_assert(element_bits < multiply_bits);
    static_assert(multiply_bits < 64);

    struct ReducerConstants
    {
        std::array<std::array<uint64_t, limbs>, log_modulus_multiples> modulus_multiples;
        std::array<uint64_t, limbs> t;
        std::array<uint64_t, limbs> negative_moduli;
        std::array<uint64_t, limbs> montgomery_factor;
    };

    constexpr ReducerConstants make_reducer_constants(
        const std::array<uint64_t, limbs>& moduli,
        const std::array<uint64_t, limbs>& montgomery_factor)
    {
        ReducerConstants result = {};
        for (std::size_t lane = 0; lane < limbs; ++lane) {
            result.t[lane] = (uint64_t{1} << element_bits) - moduli[lane];
            result.negative_moduli[lane] = uint64_t{0} - moduli[lane];
            result.montgomery_factor[lane] = montgomery_factor[lane];
            for (std::size_t index = 0; index < log_modulus_multiples; ++index) {
                result.modulus_multiples[index][lane] = moduli[lane] << index;
            }
        }
        return result;
    }

    constexpr std::array<uint64_t, limbs> make_shifted_correction(
        const std::array<uint64_t, limbs>& correction,
        const std::array<uint64_t, limbs>& output_moduli)
    {
        std::array<uint64_t, limbs> result = {};
        for (std::size_t lane = 0; lane < limbs; ++lane) {
            result[lane] = mod_mul(
                correction[lane], uint64_t{1} << multiply_bits, output_moduli[lane]);
        }
        return result;
    }

    alignas(64) inline constexpr ReducerConstants reducer_m1_constants =
        make_reducer_constants(generated::moduli1, generated::reducer1_mont_factor);
    alignas(64) inline constexpr ReducerConstants reducer_m2_constants =
        make_reducer_constants(generated::moduli2, generated::reducer2_mont_factor);
    alignas(64) inline constexpr std::array<uint64_t, limbs> r2_correction_shift =
        make_shifted_correction(generated::r2_correction, generated::moduli1);
    alignas(64) inline constexpr std::array<uint64_t, limbs> convert_to_correction_shift =
        make_shifted_correction(generated::convert_to_correction, generated::moduli2);

    INLINE void debug_assert_residues_below(
        const AVXVector<limbs>& value,
        const std::array<uint64_t, limbs>& moduli,
        uint64_t multiple)
    {
#ifndef NDEBUG
        const auto lanes = value.to_array();
        for (std::size_t lane = 0; lane < limbs; ++lane) {
            BLST_ASSERT(lanes[lane] < multiple * moduli[lane]);
        }
#else
        (void)value;
        (void)moduli;
        (void)multiple;
#endif
    }

    INLINE void debug_assert_radix_50_digits(const AVXVector<limbs>& value)
    {
#ifndef NDEBUG
        constexpr uint64_t radix = uint64_t{1} << element_bits;
        const auto lanes = value.to_array();
        for (const uint64_t lane : lanes) {
            BLST_ASSERT(lane < radix);
        }
#else
        (void)value;
#endif
    }

    struct ReducerM1Parameters
    {
        static constexpr const auto& moduli() noexcept { return generated::moduli1; }
        static constexpr const auto& constants() noexcept { return reducer_m1_constants; }
    };

    struct ReducerM2Parameters
    {
        static constexpr const auto& moduli() noexcept { return generated::moduli2; }
        static constexpr const auto& constants() noexcept { return reducer_m2_constants; }
    };

    template<class Parameters>
    class FixedMontgomeryReduce
    {
        template<int64_t upper_bound, int64_t stop, int index>
        static INLINE void min_sub_peel(AVXVector<limbs>& data)
        {
            if constexpr (index >= 0) {
                constexpr int64_t multiple = int64_t{1} << index;
                if constexpr ((upper_bound - stop) >= multiple) {
                    const AVXVector<limbs> modulus_multiple(
                        Parameters::constants().modulus_multiples[index]);
                    data = data.min(data.sub(modulus_multiple));
                }
                if constexpr (index > 0) {
                    min_sub_peel<upper_bound, stop, index - 1>(data);
                }
            }
        }

    public:
        template<int64_t lower_bound, int64_t upper_bound>
        static INLINE auto mask_reduce(
            const BoundedElement<Bounds<lower_bound, upper_bound>, limbs, element_bits>& value)
        {
            static_assert(lower_bound >= 0);
            // Per lane, split x at 2^50 and replace each high digit by
            // high * (2^50 - modulus). The result is congruent to x modulo
            // the lane modulus and lies below two moduli. This is an ordinary
            // 64-bit multiply-add, not an RNS-basis change.
            const AVXVector<limbs> t(Parameters::constants().t);
            const auto data = value.data.and_scalar(get_scalar_mask(element_bits))
                .mullo(value.data.srli(element_bits), t);
            return BoundedElement<Bounds<0, 2>, limbs, element_bits>(data);
        }

        template<int64_t upper_bound, int stop = 1>
        static INLINE auto min_sub_reduce(
            const BoundedElement<Bounds<0, upper_bound>, limbs, element_bits>& value)
        {
            static_assert(stop >= 1);
            static_assert(upper_bound >= stop);
            static_assert(
                (upper_bound - stop) < (int64_t{1} << log_modulus_multiples));
            auto data = value.data;
            // Unsigned min implements a conditional subtraction: when data is
            // smaller, data - multiple underflows and min retains data.
            min_sub_peel<upper_bound, stop, log_modulus_multiples - 1>(data);
            return BoundedElement<Bounds<0, stop>, limbs, element_bits>(data);
        }

        template<int64_t upper_bound>
        static INLINE auto reduce_small(
            const BoundedElement<Bounds<0, upper_bound>, limbs, element_bits>& high,
            const AVXVector<limbs>& low)
        {
            // One radix-2^52 Montgomery step in every lane. The low half
            // chooses q = low * modulus^-1 (mod 2^52); subtracting the high
            // half of q*modulus cancels the low radix digit. The +modulus
            // embedded in negative_moduli keeps the unsigned representative
            // nonnegative. P1/P2 must prove both the congruence and this bound.
            const AVXVector<limbs> montgomery_factor(
                Parameters::constants().montgomery_factor);
            const AVXVector<limbs> negative_moduli(
                Parameters::constants().negative_moduli);
            const AVXVector<limbs> moduli(Parameters::moduli());
            const auto quotient = AVXVector<limbs>(0).mullo(low, montgomery_factor);
            const auto negative_high = negative_moduli.mulhi(quotient, moduli);
            const auto bounded_negative_high =
                BoundedElement<Bounds<-1, 0>, limbs, element_bits>(negative_high);
            constexpr int64_t multiple = int64_t{1} << (multiply_bits - element_bits);
            const auto rebounded_high =
                BoundedElement<Bounds<0, ceil_div(upper_bound, multiple)>, limbs, element_bits>(
                    high.data);
            return rebounded_high - bounded_negative_high;
        }

        template<int64_t upper_bound>
        static INLINE auto reduce_wide(
            const WideElement<Bounds<0, upper_bound>, limbs, element_bits>& value)
        {
            if constexpr (upper_bound > 1) {
                const auto lower_high = value.low.data.srli(multiply_bits);
                const auto upper =
                    BoundedElement<Bounds<0, upper_bound>, limbs, element_bits>(
                        value.high.data.add(lower_high));
                return reduce_small(upper, value.low.data);
            } else {
                return reduce_small(value.high, value.low.data);
            }
        }

        template<int index = 0>
        static INLINE auto moduli()
        {
            static_assert(index >= 0 && index < log_modulus_multiples);
            return BoundedElement<Bounds<1 << index, 1 << index>, limbs, element_bits>(
                AVXVector<limbs>(Parameters::constants().modulus_multiples[index]));
        }

        template<int stop = 1, int64_t upper_bound, int64_t lower_bound>
        static INLINE auto reduce_auto(
            const BoundedElement<Bounds<lower_bound, upper_bound>, limbs, element_bits>& value)
        {
            if constexpr (lower_bound < 0) {
                return reduce_auto<stop>(value + moduli<ceil_log2(-lower_bound)>());
            } else if constexpr (upper_bound <= stop) {
                return value;
            } else if constexpr (lower_bound > 0) {
                return reduce_auto<stop>(value.template recast<0, upper_bound>());
            } else if constexpr (stop == 1) {
                if constexpr (upper_bound <= 4) {
                    return min_sub_reduce<upper_bound, stop>(value);
                } else {
                    return min_sub_reduce<2, 1>(mask_reduce(value));
                }
            } else {
                if constexpr (upper_bound <= stop * 2) {
                    return min_sub_reduce<upper_bound, stop>(value);
                } else {
                    return mask_reduce(value);
                }
            }
        }

        template<int stop = 1, int64_t upper_bound, int64_t lower_bound>
        static INLINE auto reduce_auto(
            const WideElement<Bounds<lower_bound, upper_bound>, limbs, element_bits>& value)
        {
            static_assert(lower_bound >= 0);
            return reduce_auto<stop>(reduce_wide(value.template recast<0, upper_bound>()));
        }
    };

    struct R2MatrixParameters
    {
        static constexpr const auto& matrix() noexcept { return generated::r2_rns_mat; }
        static constexpr const auto& correction() noexcept { return generated::r2_correction; }
        static constexpr const auto& correction_shift() noexcept { return r2_correction_shift; }
        static constexpr const auto& quotient_estimates() noexcept
        {
            return generated::r2_shifted_quotient_estimations;
        }
        static INLINE void debug_assert_input(const AVXVector<limbs>& value)
        {
            debug_assert_residues_below(value, generated::moduli2, 2);
        }
    };

    struct ConvertToMatrixParameters
    {
        static constexpr const auto& matrix() noexcept { return generated::convert_to_rns_mat; }
        static constexpr const auto& correction() noexcept
        {
            return generated::convert_to_correction;
        }
        static constexpr const auto& correction_shift() noexcept
        {
            return convert_to_correction_shift;
        }
        static constexpr const auto& quotient_estimates() noexcept
        {
            return generated::convert_to_shifted_quotient_estimations;
        }
        static INLINE void debug_assert_input(const AVXVector<limbs>& value)
        {
            debug_assert_radix_50_digits(value);
        }
    };

    struct ConvertFromMatrixParameters
    {
        static constexpr const auto& matrix() noexcept { return generated::convert_from_rns_mat; }
        static constexpr const auto& correction() noexcept
        {
            return generated::convert_from_correction;
        }
        static constexpr const auto& quotient_estimates() noexcept
        {
            return generated::convert_from_shifted_quotient_estimations;
        }
        static INLINE void debug_assert_input(const AVXVector<limbs>& value)
        {
            debug_assert_residues_below(value, generated::moduli2, 2);
        }
    };

    template<class Parameters>
    class FixedRNSMatrix
    {
        // All three entry points below are in/out accumulators. They add the
        // matrix product and its approximate-quotient correction to the
        // caller-provided high/low pair; callers pass zero for a pure basis
        // conversion and an existing product for fused MatrixNoK reduction.
        // The shifted variants use correction * 2^52 modulo each output
        // modulus so IFMA high/low accumulation stays in the same radix.
        template<std::size_t index>
        static INLINE void accumulate_registers(
            AVXVector<limbs>& out_high,
            AVXVector<limbs>& out_low,
            const AVXVector<limbs>& residues,
            uint128_t& quotient_numerator)
        {
            static_assert(index < limbs);
            const uint64_t scalar = residues.template extract<index>();
            const AVXScalar broadcast = Scalar(scalar);
            const AVXVector<limbs> row(Parameters::matrix()[index]);
            out_high = out_high.mulhi_scalar(row, broadcast);
            out_low = out_low.mullo_scalar(row, broadcast);
            quotient_numerator += static_cast<uint128_t>(scalar) *
                static_cast<uint128_t>(Parameters::quotient_estimates()[index]);
            if constexpr (index + 1 < limbs) {
                accumulate_registers<index + 1>(
                    out_high, out_low, residues, quotient_numerator);
            }
        }

        static INLINE void add_correction(
            AVXVector<limbs>& out_high,
            AVXVector<limbs>& out_low,
            uint128_t quotient_numerator)
        {
            const AVXVector<limbs> correction(Parameters::correction());
            AVXScalar quotient = Scalar(
                static_cast<uint64_t>(quotient_numerator >> quotient_precision));
            out_high = out_high.mulhi_scalar(correction, quotient);
            out_low = out_low.mullo_scalar(correction, quotient);
            quotient = scalar_shift(quotient, multiply_bits);
            out_high = out_high.mullo_scalar(correction, quotient);
            const auto lower_product =
                AVXVector<limbs>(0).mulhi_scalar(correction, quotient);
            out_high = out_high.add(lower_product.slli(multiply_bits));
        }

    public:
        static INLINE void reduce_raw(
            const AVXVector<limbs>& residues,
            AVXVector<limbs>& out_high,
            AVXVector<limbs>& out_low)
        {
            // Raw matrix operations require disjoint input and accumulator
            // objects. Their IFMA additions have uint64_t wrapping semantics;
            // production callers separately prove that their bounds do not wrap.
            BLST_ASSERT(&residues != &out_high);
            BLST_ASSERT(&residues != &out_low);
            BLST_ASSERT(&out_high != &out_low);
            Parameters::debug_assert_input(residues);
            std::array<uint64_t, limbs> scalars;
            residues.store(scalars.data());
            uint128_t quotient_numerator = 0;
            BLST_UNROLL
            for (std::size_t index = 0; index < limbs; ++index) {
                const AVXScalar broadcast = Scalar(scalars[index]);
                const AVXVector<limbs> row(Parameters::matrix()[index]);
                out_high = out_high.mulhi_scalar(row, broadcast);
                out_low = out_low.mullo_scalar(row, broadcast);
                quotient_numerator += static_cast<uint128_t>(scalars[index]) *
                    static_cast<uint128_t>(Parameters::quotient_estimates()[index]);
            }
            add_correction(out_high, out_low, quotient_numerator);
        }

        static INLINE void reduce_raw_shifted(
            const AVXVector<limbs>& residues,
            AVXVector<limbs>& out_high,
            AVXVector<limbs>& out_low)
        {
            BLST_ASSERT(&residues != &out_high);
            BLST_ASSERT(&residues != &out_low);
            BLST_ASSERT(&out_high != &out_low);
            Parameters::debug_assert_input(residues);
            uint128_t quotient_numerator = 0;
            accumulate_registers<0>(
                out_high, out_low, residues, quotient_numerator);
            const AVXVector<limbs> correction(Parameters::correction());
            const AVXVector<limbs> correction_shift(Parameters::correction_shift());
            AVXScalar quotient = Scalar(
                static_cast<uint64_t>(quotient_numerator >> quotient_precision));
            out_high = out_high.mulhi_scalar(correction, quotient);
            out_low = out_low.mullo_scalar(correction, quotient);
            quotient = scalar_shift(quotient, multiply_bits);
            out_low = out_low.mullo_scalar(correction_shift, quotient);
            out_high = out_high.mulhi_scalar(correction_shift, quotient);
        }

        template<std::size_t batch>
        static INLINE void reduce_raw_batch_shifted(
            const std::array<AVXVector<limbs>, batch>& residues,
            std::array<AVXVector<limbs>, batch>& out_high,
            std::array<AVXVector<limbs>, batch>& out_low)
        {
            static_assert(batch > 0);
            BLST_ASSERT(&residues != &out_high);
            BLST_ASSERT(&residues != &out_low);
            BLST_ASSERT(&out_high != &out_low);
            alignas(64) std::array<std::array<uint64_t, limbs>, batch> scalars;
            BLST_UNROLL
            for (std::size_t item = 0; item < batch; ++item) {
                Parameters::debug_assert_input(residues[item]);
                residues[item].store(scalars[item].data());
            }
            std::array<uint128_t, batch> quotient_numerators = {};
            for (std::size_t index = 0; index < limbs; ++index) {
                const AVXVector<limbs> row(Parameters::matrix()[index]);
                BLST_UNROLL
                for (std::size_t item = 0; item < batch; ++item) {
                    const AVXScalar broadcast = Scalar(scalars[item][index]);
                    quotient_numerators[item] +=
                        static_cast<uint128_t>(scalars[item][index]) *
                        static_cast<uint128_t>(Parameters::quotient_estimates()[index]);
                    out_high[item] = out_high[item].mulhi_scalar(row, broadcast);
                    out_low[item] = out_low[item].mullo_scalar(row, broadcast);
                }
            }

            const AVXVector<limbs> correction(Parameters::correction());
            const AVXVector<limbs> correction_shift(Parameters::correction_shift());
            AVXScalarArray<batch> quotients;
            BLST_UNROLL
            for (std::size_t item = 0; item < batch; ++item) {
                quotients[item] = Scalar(
                    static_cast<uint64_t>(
                        quotient_numerators[item] >> quotient_precision));
                out_high[item] = out_high[item].mulhi_scalar(correction, quotients[item]);
                out_low[item] = out_low[item].mullo_scalar(correction, quotients[item]);
                quotients[item] = scalar_shift(quotients[item], multiply_bits);
                out_low[item] = out_low[item].mullo_scalar(
                    correction_shift, quotients[item]);
                out_high[item] = out_high[item].mulhi_scalar(
                    correction_shift, quotients[item]);
            }
        }
    };

    using ReducerM1 = FixedMontgomeryReduce<ReducerM1Parameters>;
    using ReducerM2 = FixedMontgomeryReduce<ReducerM2Parameters>;
    using R2Matrix = FixedRNSMatrix<R2MatrixParameters>;
    using ConvertToMatrix = FixedRNSMatrix<ConvertToMatrixParameters>;
    using ConvertFromMatrix = FixedRNSMatrix<ConvertFromMatrixParameters>;
}
