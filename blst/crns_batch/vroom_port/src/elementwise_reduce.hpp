#pragma once
#include "../cpu/vector/vector_impl.hpp"
#include "bounded_element.hpp"
#include "ring_element.hpp"
#include "constexpr_utils.hpp"
#include "preprocess.hpp"
#include <cassert>

template<int limbs, int element_bits = 52, int log_multiples = 4>
class MontgomeryReduce2 {
    static constexpr int mulbits = 52;
    static constexpr int hi_bit_shift = element_bits + element_bits - mulbits;

    AVXScalar hi_mask;
    AVXScalar lo_mask;

    std::array<AVXVector<limbs>, log_multiples> moduli_multiples;
    AVXVector<limbs> t;
    AVXVector<limbs> mont_arr;
    AVXVector<limbs> neg_moduli;
    AVXVector<limbs> mont_convert_factor;

    public:
    INLINE MontgomeryReduce2(const std::array<uint64_t, limbs> &moduli) 
        : hi_mask(get_scalar_mask(hi_bit_shift))
        , lo_mask(get_scalar_mask(mulbits))
    {
        std::array<std::array<uint64_t,limbs>, log_multiples> moduli_multiples{};
        std::array<uint64_t, limbs> mont_arr{};
        std::array<uint64_t, limbs> neg_moduli{};
        std::array<uint64_t, limbs> t{};
        std::array<uint64_t, limbs> mont_convert_factor{};

        static constexpr uint64_t r_mod = 1ull << mulbits;
        
        for (int i = 0; i < limbs; i++) {
            t[i] = (1ull << element_bits) - moduli[i];
            mont_arr[i] = mod_inverse(moduli[i], r_mod);
            neg_moduli[i] = -moduli[i];
            mont_convert_factor[i] = ((uint128_t)(1) << (2*mulbits)) % moduli[i];
            for (int j = 0; j < log_multiples; j++) {
                moduli_multiples[j][i] = moduli[i] << j;
            }
        }
        this->t.load(t.data());
        this->mont_arr.load(mont_arr.data());
        this->neg_moduli.load(neg_moduli.data());
        this->mont_convert_factor.load(mont_convert_factor.data());
        for (int j = 0; j < log_multiples; j++) {
            this->moduli_multiples[j].load(moduli_multiples[j].data());
        }
    }

    template<int64_t lower_bound, int64_t upper_bound>
    INLINE auto mask_reduce(const BoundedElement<Bounds<lower_bound, upper_bound>, limbs, element_bits> &a) const {
        static_assert(lower_bound >= 0, "Lower bound must be non-negative for mask reduction");
        auto data = a.data.and_scalar(get_scalar_mask(element_bits)).mullo(a.data.srli(element_bits), t);
        return BoundedElement<Bounds<0, 2>, limbs, element_bits>(data);
    }

    template<int64_t upper_bound, int stop = 1>
    INLINE auto min_sub_reduce(const BoundedElement<Bounds<0, upper_bound>, limbs, element_bits> &a) const {
        static_assert(stop >= 1, "Stop must be at least 1");
        static_assert(upper_bound <= (1 << (log_multiples+1)), "Upper bound cannot exceed log multiples");
        auto data = a.data;
        if constexpr ((upper_bound > 8) && (stop < 16)) {
            data = data.min(data.sub(moduli_multiples[3])); // -= 8p
        }
        if constexpr ((upper_bound > 4) && (stop < 8)) {
            data = data.min(data.sub(moduli_multiples[2])); // -= 4p
        }
        if constexpr ((upper_bound > 2) && (stop < 4)) {
            data = data.min(data.sub(moduli_multiples[1])); // -= 2p
        }
        if constexpr ((upper_bound > 1) && (stop < 2)) {
            data = data.min(data.sub(moduli_multiples[0])); // -= p
        }
        return BoundedElement<Bounds<0, stop>, limbs, element_bits>(data);
    }

    template<int64_t upper_bound>
    INLINE auto reduce_small(const BoundedElement<Bounds<0, upper_bound>, limbs, element_bits> &high, const AVXVector<limbs> &low) const {
        auto q = AVXVector<limbs>(0).mullo(low, mont_arr);
        // For b bits: h = q * moduli = [0, 2**b - t] because 2**52 * (2**b - t) >> 52 = 2**b - t
        // h_n = [-(2**b - t), 0] so this is [-1, 0] in the output bounds
        auto h_n = neg_moduli.mulhi(q, moduli_multiples[0]);
        auto h_n_bounded = BoundedElement<Bounds<-1, 0>, limbs, element_bits>(h_n);
        // For bounds multiple k:
        // k*(2**b - t)**2 < k*(2**b - t)2**b 
        // Taking the upper bits by shifting by 52:
        // (k*(2**b - t)**2 >> 52) < (k*(2**b - t) >> (52 - b)) <= ceil_div(k, 1 << (52 - b)) * (2**b - t)
        constexpr int64_t multiple = 1ull << (52 - element_bits);
        auto hi_rebounded = BoundedElement<Bounds<0, ceil_div(upper_bound, multiple)>, limbs, element_bits>(high.data);
        // Note that in inclusive bounds, this can return 2**b - t (Only occurs for the all zero input).
        return hi_rebounded - h_n_bounded;
    }

    template<int64_t upper_bound>
    INLINE auto reduce_wide(const WideElement<Bounds<0, upper_bound>, limbs, element_bits> &a) const {
        if constexpr (upper_bound > 1) {
            // Move the excess lower bits to the upper bits
            // Requires unsigned interpretation (lower_bound >= 0)
            auto lo_hi = a.low.data.srli(mulbits);
            // Adding the lower bits back to the higher bit; does not change the value when interpreted as a uint128_t.
            // No need to mask the lower bits - this happens implicitly in the mullo instruction
            auto new_upper = BoundedElement<Bounds<0, upper_bound>, limbs, element_bits>(a.high.data.add(lo_hi));
            return reduce_small(new_upper, a.low.data);
        } else {
            return reduce_small(a.high, a.low.data);
        }
    }

    template<int index=0>
    INLINE auto moduli() const {
        static_assert(index >= 0, "Index must be non-negative");
        static_assert(index < log_multiples, "Index must be less than log multiples; increase precomputed table size");
        // Returns an exact multiple of the modulus, with tight inclusive bounds.
        return BoundedElement<Bounds<1 << index, 1 << index>, limbs, element_bits>(moduli_multiples[index]);
    }

    template<int stop=1, int64_t upper_bound, int64_t lower_bound>
    INLINE auto reduce_auto(const BoundedElement<Bounds<lower_bound, upper_bound>, limbs, element_bits> & a) const {
        if constexpr (lower_bound < 0) {
            return reduce_auto<stop>(a + moduli<ceil_log2(-lower_bound)>());
        } else if constexpr (upper_bound <= stop) {
            return a;
        } else if constexpr (lower_bound > 0) {
            return reduce_auto<stop>(a.template recast<0, upper_bound>());
        } else if constexpr (stop == 1) {
            // mask reduce: slli, and, mullo = 1 + 0.5 + 0.5 =2 cycles
            // minsub = min + sub = 1 + 0.5 = 1.5 cycles
            if constexpr(upper_bound <= 4) {
                // up to 2 min sub required = 3
                return min_sub_reduce<upper_bound, stop>(a);
            } else {
                // 1 mask, 1 minsub, total 3.5, any more than 2 min sub should be masked reduce.
                return min_sub_reduce<2, 1>(mask_reduce(a));
            }
        } else {
            if constexpr (upper_bound <= stop * 2) {
                // 1 min sub required: 1.5
                return min_sub_reduce<upper_bound, stop>(a);
            } else {
                // if more than 1 min sub would be required, do mask reduce: 2 - stop is at least 2 so mask reduce is full sufficient.
                return mask_reduce(a);
            }
        }
    }

    template<int stop=1, int64_t upper_bound, int64_t lower_bound>
    INLINE auto reduce_auto(const WideElement<Bounds<lower_bound, upper_bound>, limbs, element_bits> &a) const {
        static_assert(lower_bound >= 0, "Lower bound must be non-negative; add offsets to make positive");
        return reduce_auto<stop>(reduce_wide(a.template recast<0, upper_bound>()));
    }

    INLINE auto to_mont(AVXVector<limbs> &a) const {
        // Multiply and reduce
        using BoundedElementType = BoundedElement<Bounds<0, 1>, limbs, element_bits>;
        WideElement<Bounds<0, 1>, limbs, element_bits> bounded(BoundedElementType(AVXVector<limbs>(0).mulhi(a, mont_convert_factor)), BoundedElementType(AVXVector<limbs>(0).mullo(a, mont_convert_factor)));
        return min_sub_reduce<2, 1>(reduce_wide(bounded));
    }

    INLINE auto from_mont(BoundedElement<Bounds<0, 1>, limbs, element_bits> &a) const {
        // Multiply by 1 (omitted since multiplication by 1 does nothing)
        // Then Montgomery reduce: create WideElement with high=0, low=a and reduce
        using BoundedElementType = BoundedElement<Bounds<0, 1>, limbs, element_bits>;
        WideElement<Bounds<0, 1>, limbs, element_bits> bounded(BoundedElementType(AVXVector<limbs>(0)), a);
        // By articially setting the high upper bound to [0, 1] (rather than [0, 0]), the output will be in the standard representation.
        return min_sub_reduce<2, 1>(reduce_wide(bounded)).data;
    }

};