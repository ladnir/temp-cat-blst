#pragma once
#include <algorithm>
#include <cstdint>
#include <climits>
#include "constexpr_utils.hpp"

// All ranges are inclusive of the lower bound and upper bound
// This makes subtraction of bounds simpler, and allows us to express [1, 1] as the modulus itself.
// It does mean that [0, 1] could represent 0 by the modulus, while everything else has exactly one representative.
template<int64_t lower_bound, int64_t upper_bound>
class Bounds {
    static_assert(lower_bound <= upper_bound, "Lower bound must be less than or equal to upper bound");

    public:
    // Expose bounds as static constexpr members for access
    static constexpr int64_t lower = lower_bound;
    static constexpr int64_t upper = upper_bound;
    
    // Addition: [a, b] + [c, d] = [a+c, b+d]
    // Compute sums in int128_t to avoid overflow during calculation, then check they fit in int64_t
    template<int64_t other_lower, int64_t other_upper>
    constexpr auto operator+(const Bounds<other_lower, other_upper> &) const {
        // Use template parameters directly instead of function parameter for clang constexpr compatibility
        constexpr int128_t new_lower = (int128_t)(lower) + (int128_t)(other_lower);
        constexpr int128_t new_upper = (int128_t)(upper) + (int128_t)(other_upper);
        static_assert(new_lower >= INT64_MIN, "Sum underflows 64 bit type");
        static_assert(new_upper <= INT64_MAX, "Sum overflows 64 bit type");
        return Bounds<(int64_t)new_lower, (int64_t)new_upper>{};
    }

    // Subtraction: [a, b] - [c, d] = [a-d, b-c]
    // Compute differences in int128_t to avoid overflow during calculation, then check they fit in int64_t
    template<int64_t other_lower, int64_t other_upper>
    constexpr auto operator-(const Bounds<other_lower, other_upper> &) const {
        // Use template parameters directly instead of function parameter for clang constexpr compatibility
        constexpr int128_t new_lower = (int128_t)(lower) - (int128_t)(other_upper);
        constexpr int128_t new_upper = (int128_t)(upper) - (int128_t)(other_lower);
        static_assert(new_lower >= INT64_MIN, "Difference underflows 64 bit type");
        static_assert(new_upper <= INT64_MAX, "Difference overflows 64 bit type");
        return Bounds<(int64_t)new_lower, (int64_t)new_upper>{};
    }

    // Multiplication: [a, b] * [c, d] = [min(a*c, a*d, b*c, b*d), max(a*c, a*d, b*c, b*d)]
    // Need to consider all four products to handle mixed signs correctly
    // Compute products in int128_t to avoid overflow during calculation, then check they fit in int64_t
    template<int64_t other_lower, int64_t other_upper>
    constexpr auto operator*(const Bounds<other_lower, other_upper> &) const {
        // Use template parameters directly instead of function parameter for clang constexpr compatibility
        constexpr int128_t ac = (int128_t)(lower) * (int128_t)(other_lower);
        constexpr int128_t ad = (int128_t)(lower) * (int128_t)(other_upper);
        constexpr int128_t bc = (int128_t)(upper) * (int128_t)(other_lower);
        constexpr int128_t bd = (int128_t)(upper) * (int128_t)(other_upper);
        // Use nested min/max for constexpr compatibility
        constexpr int128_t min_val = std::min(std::min(ac, ad), std::min(bc, bd));
        constexpr int128_t max_val = std::max(std::max(ac, ad), std::max(bc, bd));
        static_assert(min_val >= INT64_MIN, "Product underflows 64 bit type");
        static_assert(max_val <= INT64_MAX, "Product overflows 64 bit type");
        return Bounds<(int64_t)min_val, (int64_t)max_val>{};
    }

    // Unary negation: -[a, b] = [-b, -a]
    constexpr auto operator-() const {
        return Bounds<(-upper), (-lower)>{};
    }

    template<int64_t left_shift>
    constexpr auto shift_left() const {
        static_assert(left_shift >= 0, "Left shift must be non-negative");
        // Use multiplication by power of 2 instead of bit shift to avoid undefined behavior
        // with negative numbers. Shifting negative numbers left is undefined in C++.
        constexpr int64_t multiplier = static_cast<int64_t>(1) << left_shift;
        constexpr int128_t new_lower = (int128_t)(lower) * (int128_t)(multiplier);
        constexpr int128_t new_upper = (int128_t)(upper) * (int128_t)(multiplier);
        static_assert(new_lower >= INT64_MIN, "Shift left underflows 64 bit type");
        static_assert(new_upper <= INT64_MAX, "Shift left overflows 64 bit type");
        return Bounds<(int64_t)new_lower, (int64_t)new_upper>{};
    }

    template<int64_t new_lower, int64_t new_upper>
    constexpr auto recast() const {
        static_assert(new_lower <= lower, "New lower bound must be less than or equal to lower bound");
        static_assert(new_upper >= upper, "New upper bound must be greater than or equal to upper bound");
        return Bounds<new_lower, new_upper>{};
    }
};
