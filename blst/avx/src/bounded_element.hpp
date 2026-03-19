#pragma once
#include "bounds.hpp"
#include "../cpu/vector/vector_impl.hpp"
#include "preprocess.hpp"
#include <string>
#include <iostream>

// All bounds are inclusive of both the upper and lower endpoint
// Implcit multiples
// Normal: [low * modulus, upper * modulus]
// Wide: [low * modulus**2, upper * modulus**2]
// However, we have an additional requirement for wide that the low 52 bits are in
// the range [2^52 * low, 2^52 * upper].  This prevents underflow when subtracting the high/low parts without borrowing.
// Therefore, except for [0, 0], no wide element can have low == upper.

template<class BoundsType, int limbs, int element_bits=52>
class BoundedElement {
    public:
    static constexpr BoundsType bounds = BoundsType{};
    AVXVector<limbs> data;

    static constexpr int LIMBS = limbs;
    static constexpr int ELEMENT_BITS = element_bits;
    
    static constexpr int64_t max_bound = (1ull << (64 - element_bits - 1));
    // 52(50) bit number stored in 64 bit type
    static_assert(BoundsType::lower > -max_bound, "Value can underflow type");
    static_assert(BoundsType::upper < max_bound, "Value can overflow type");

    // Constructor
    INLINE BoundedElement() = default;
    INLINE explicit BoundedElement(const AVXVector<limbs> &d) : data(d) {}

    template<class OtherBounds>
    INLINE auto operator+(const BoundedElement<OtherBounds, limbs, element_bits> &other) const {
        return BoundedElement<decltype(bounds + other.bounds), limbs, element_bits>(data.add(other.data));
    }

    template<class OtherBounds>
    INLINE auto operator-(const BoundedElement<OtherBounds, limbs, element_bits> &other) const {
        return BoundedElement<decltype(bounds - other.bounds), limbs, element_bits>(data.sub(other.data));
    }
    
    template<int64_t left_shift>
    INLINE auto shift_left() const {
        auto shifted_bounds = bounds.template shift_left<left_shift>();
        return BoundedElement<decltype(shifted_bounds), limbs, element_bits>(data.slli(left_shift));
    }

    INLINE static void check_mult_bounds() {
        static_assert(bounds.lower >= 0, "Lower bound must be non-negative for unsigned multiplication");
        // We use AVX IFMA 52-bit multiplication hardware.
        // For 52-bit numbers: can multiply directly (bound <= 1, representing 1 * 2^52 = full 52-bit range)
        // For 50-bit numbers: can sum up to 4 numbers (bound <= 4, representing 4 * 2^50 = 2^52, still fits in 52-bit mult)
        // Any more would exceed 2^52 and overflow the 52-bit multiplication
        static_assert(bounds.upper <= (1ull << (52 - element_bits)), "Upper bound must be <= 2^(52 - element_bits) for 52-bit multiplication");
    }

    INLINE static auto zero() {
        return BoundedElement<Bounds<0, 0>, limbs, element_bits>{};
    }

    INLINE std::array<int64_t, limbs> to_array() const {
        auto uint_array = data.to_array();
        std::array<int64_t, limbs> result;
        for (int i = 0; i < limbs; i++) {
            result[i] = (int64_t)uint_array[i];
        }
        return result;
    }

    INLINE std::array<uint64_t, limbs> to_unsigned_array() const {
        auto uint_array = data.to_array();
        std::array<uint64_t, limbs> result;
        for (int i = 0; i < limbs; i++) {
            result[i] = uint_array[i];
        }
        return result;
    }

    INLINE bool check_bounds(const std::array<int64_t, limbs> &moduli, const std::string &context = "") const {
        auto data_array = to_array();
        for (int i = 0; i < limbs; i++) {
            // Inclusive lower bound, inclusive upper bound
            int64_t lower_bound = bounds.lower * moduli[i];
            int64_t upper_bound = bounds.upper * moduli[i];
            if (data_array[i] < lower_bound || data_array[i] > upper_bound) {
                if (!context.empty()) {
                    std::cerr << "check_bounds failed for " << context << ":" << std::endl;
                    std::cerr << "  Limb " << i << ": value = " << data_array[i] 
                              << ", expected range = [" << lower_bound << ", " << upper_bound << "]" << std::endl;
                    std::cerr << "  Bounds: [" << bounds.lower << ", " << bounds.upper << "]" << std::endl;
                    std::cerr << "  Modulus[" << i << "] = " << moduli[i] << std::endl;
                }
                return false;
            }
        }
        return true;
    }

    INLINE auto complete() const {
        return *this;
    }
    
    template<int64_t new_lower, int64_t new_upper>
    INLINE auto recast() const {
        auto new_bounds = bounds.template recast<new_lower, new_upper>();
        return BoundedElement<decltype(new_bounds), limbs, element_bits>(data);
    }

    INLINE auto operator*(std::integral_constant<int, 1>) const {
        return *this;
    }

    INLINE auto operator*(std::integral_constant<int, 2>) const {
        return *this + *this;
    }

    INLINE auto operator*(std::integral_constant<int, 3>) const {
        return *this + *this + *this;
    }

    INLINE auto operator*(std::integral_constant<int, 6>) const {
        return (*this * std::integral_constant<int, 3>()) * std::integral_constant<int, 2>();
    }

    INLINE auto ctime_select(bool selection, const BoundedElement &other) const {
        return BoundedElement(data.ctime_select(selection, other.data));
    }

    // Clearer API: replace self with replacement if condition is true
    // ctime_select(selection, b) returns *this when selection is true, b when false
    // So we need to negate: if condition is true, we want replacement (not *this)
    INLINE auto ctime_replace_if(bool condition, const BoundedElement &replacement) const {
        return ctime_select(!condition, replacement);
    }
};

// Forward declaration
template<class ElementType, int limbs, int element_bits>
class WideElement;

template<class Left, class Right>
class DelayedProduct {
    public:
    Left left;
    Right right;

    INLINE DelayedProduct(Left left, Right right) : left(left), right(right) {}

    INLINE auto complete() const {
        auto zero = Left::zero();
        constexpr int limbs = Left::LIMBS;
        constexpr int element_bits = Left::ELEMENT_BITS;
        return WideElement<Bounds<0, 0>, limbs, element_bits>(zero, zero) + *this;
    }
};

template<class BoundsType, int limbs, int element_bits>
class WideElement {
    using ElementType = BoundedElement<BoundsType, limbs, element_bits>;
    public:
    ElementType high;
    ElementType low;
    static constexpr BoundsType bounds = BoundsType{};

    INLINE WideElement() = default;
    INLINE explicit WideElement(const ElementType &high, const ElementType &low) : high(high), low(low) {}
    INLINE explicit WideElement(const AVXVector<limbs> &high_data, const AVXVector<limbs> &low_data) : high(high_data), low(low_data) {}
    
    template<class OtherElementType>
    INLINE auto operator+(const WideElement<OtherElementType, limbs, element_bits> &other) const {
        auto result_low = low + other.low;
        auto result_high = high + other.high;
        return WideElement<decltype(bounds + other.bounds), limbs, element_bits>(result_high, result_low);
    }

    template<class OtherElementType>
    INLINE auto operator-(const WideElement<OtherElementType, limbs, element_bits> &other) const {
        auto result_low = low - other.low;
        auto result_high = high - other.high;
        return WideElement<decltype(bounds - other.bounds), limbs, element_bits>(result_high, result_low);
    }

    template<int64_t left_shift>
    INLINE auto shift_left() const {
        auto high_shifted = high.template shift_left<left_shift>();
        auto low_shifted = low.template shift_left<left_shift>();
        return WideElement<decltype(bounds.template shift_left<left_shift>()), limbs, element_bits>(high_shifted, low_shifted);
    }

    template<class LeftBounds, class RightBounds>
    INLINE auto mulacc(const BoundedElement<LeftBounds, limbs, element_bits> &left, const BoundedElement<RightBounds, limbs, element_bits> &right) const {
        left.check_mult_bounds();
        right.check_mult_bounds();
        auto result_low_data = low.data.mullo(left.data, right.data);
        auto result_high_data = high.data.mulhi(left.data, right.data);
        using product_bounds = decltype(bounds + (left.bounds * right.bounds));
        auto result_low = BoundedElement<product_bounds, limbs, element_bits>(result_low_data);
        auto result_high = BoundedElement<product_bounds, limbs, element_bits>(result_high_data);
        return WideElement<product_bounds, limbs, element_bits>(result_high, result_low);
    }

    template<class LeftBounds, class RightBounds>
    INLINE auto operator+(const DelayedProduct<BoundedElement<LeftBounds, limbs, element_bits>, BoundedElement<RightBounds, limbs, element_bits>> &delayed_product) const {
        return mulacc(delayed_product.left, delayed_product.right);
    }

    // No action required for non delayed products.
    INLINE WideElement complete() const {
        return *this;
    }

    template<int64_t new_lower, int64_t new_upper>
    INLINE auto recast() const {
        auto new_bounds = bounds.template recast<new_lower, new_upper>();
        return WideElement<decltype(new_bounds), limbs, element_bits>(high.template recast<new_lower, new_upper>(), low.template recast<new_lower, new_upper>());
    }

    INLINE std::array<int128_t, limbs> to_array() const {
        auto low_array = low.data.to_array();
        auto high_array = high.data.to_array();
        std::array<int128_t, limbs> result;
        for (int i = 0; i < limbs; i++) {
            result[i] = (int128_t)low_array[i] + ((int128_t)high_array[i] << 52);
        }
        return result;
    }

    INLINE bool check_bounds(const std::array<int64_t, limbs> &moduli, const std::string &context = "") const {
        std::array<int64_t, limbs> low_mask;
        std::array<int128_t, limbs> reconstructed_value = to_array();
        bool high_ok = true;
        for (int i = 0; i < limbs; i++) {
            int128_t lower_bound = (int128_t)high.bounds.lower * (int128_t)moduli[i] * (int128_t)moduli[i];
            int128_t upper_bound = (int128_t)high.bounds.upper * (int128_t)moduli[i] * (int128_t)moduli[i];
            if (reconstructed_value[i] < lower_bound || reconstructed_value[i] > upper_bound) {
                if (!context.empty()) {
                    std::cerr << "check_bounds failed for " << context << ":" << std::endl;
                    std::cerr << "  Limb " << i << ":  hi value = " << high.data.to_array()[i] << ", lo value = " << low.data.to_array()[i] << std::endl;
                    std::cerr << "  Bounds: [" << high.bounds.lower << ", " << high.bounds.upper << "]" << std::endl;
                    std::cerr << "  Modulus[" << i << "] = " << moduli[i] << std::endl;
                }
                high_ok = false;
            }
            low_mask[i] = (1ull << 52) - 1;
        }
        std::string low_context = context.empty() ? "low" : context + ".low";
        bool low_ok = low.check_bounds(low_mask, low_context);
        return low_ok && high_ok;
    }
};

template<int limbs, int element_bits>
INLINE auto WideZero() {
    auto zero = BoundedElement<Bounds<0, 0>, limbs, element_bits>::zero();
    return WideElement<Bounds<0, 0>, limbs, element_bits>(zero, zero);
}