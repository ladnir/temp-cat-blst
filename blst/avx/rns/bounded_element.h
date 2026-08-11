#pragma once
#include "bounds.h"
#include "../cpu/vector/vector_impl.hpp"
#include "../../common/defines.h"

// All bounds are inclusive of both the upper and lower endpoint
// Implicit multiples
// Normal: [low * modulus, upper * modulus]
// Wide: [low * modulus**2, upper * modulus**2]
// A wide value additionally keeps its low accumulator in the interval
// [2^52 * low, 2^52 * upper]. WideElement enforces the corresponding capacity
// and borrow-free subtraction requirements at compile time.

namespace blst::avx {

template<class BoundsType, int limbs, int element_bits=52>
class BoundedElement {
    public:
    static constexpr BoundsType bounds = BoundsType{};
    AVXVector<limbs> data;

    static constexpr int LIMBS = limbs;
    static constexpr int ELEMENT_BITS = element_bits;

    static_assert(limbs > 0, "A bounded element requires at least one lane");
    static_assert(element_bits > 0 && element_bits <= 52,
        "Bounded elements require a radix between 2^1 and the IFMA 2^52 radix");
    static constexpr int64_t max_bound =
        int64_t{1} << (64 - element_bits - 1);
    // An element_bits-wide value is stored in a signed 64-bit lane.
    static_assert(BoundsType::lower > -max_bound, "Value can underflow type");
    static_assert(BoundsType::upper < max_bound, "Value can overflow type");

    // Default initialization intentionally leaves the vector uninitialized for
    // fixed arrays whose callers overwrite every element before use.
    BLST_FORCEINLINE BoundedElement() = default;
    constexpr BLST_FORCEINLINE explicit BoundedElement(const AVXVector<limbs> &d) : data(d) {}

    template<class OtherBounds>
    BLST_FORCEINLINE auto operator+(const BoundedElement<OtherBounds, limbs, element_bits> &other) const {
        return BoundedElement<decltype(bounds + other.bounds), limbs, element_bits>(data.add(other.data));
    }

    template<class OtherBounds>
    BLST_FORCEINLINE auto operator-(const BoundedElement<OtherBounds, limbs, element_bits> &other) const {
        return BoundedElement<decltype(bounds - other.bounds), limbs, element_bits>(data.sub(other.data));
    }
    
    template<int64_t left_shift>
    BLST_FORCEINLINE auto shift_left() const {
        auto shifted_bounds = bounds.template shift_left<left_shift>();
        return BoundedElement<decltype(shifted_bounds), limbs, element_bits>(data.slli(left_shift));
    }

    BLST_FORCEINLINE static void check_mult_bounds() {
        static_assert(bounds.lower >= 0, "Lower bound must be non-negative for unsigned multiplication");
        // We use AVX IFMA 52-bit multiplication hardware.
        // For 52-bit numbers: can multiply directly (bound <= 1, representing 1 * 2^52 = full 52-bit range)
        // For 50-bit numbers: can sum up to 4 numbers (bound <= 4, representing 4 * 2^50 = 2^52, still fits in 52-bit mult)
        // Any more would exceed 2^52 and overflow the 52-bit multiplication
        static_assert(bounds.upper <= (int64_t{1} << (52 - element_bits)), "Upper bound must be <= 2^(52 - element_bits) for 52-bit multiplication");
    }

    BLST_FORCEINLINE static auto zero() {
        using Zero = BoundedElement<Bounds<0, 0>, limbs, element_bits>;
        return Zero(AVXVector<limbs>(0));
    }

    BLST_FORCEINLINE auto complete() const {
        return *this;
    }
    
    template<int64_t new_lower, int64_t new_upper>
    BLST_FORCEINLINE auto recast() const {
        auto new_bounds = bounds.template recast<new_lower, new_upper>();
        return BoundedElement<decltype(new_bounds), limbs, element_bits>(data);
    }

    BLST_FORCEINLINE auto operator*(std::integral_constant<int, 1>) const {
        return *this;
    }

    BLST_FORCEINLINE auto operator*(std::integral_constant<int, 2>) const {
        return *this + *this;
    }

    BLST_FORCEINLINE auto operator*(std::integral_constant<int, 3>) const {
        return *this + *this + *this;
    }

    BLST_FORCEINLINE auto operator*(std::integral_constant<int, 6>) const {
        return (*this * std::integral_constant<int, 3>()) * std::integral_constant<int, 2>();
    }

    BLST_FORCEINLINE auto ctime_select(bool selection, const BoundedElement &other) const {
        return BoundedElement(data.ctime_select(selection, other.data));
    }

    BLST_FORCEINLINE bool_t is_zero() const {
        return data.is_zero();
    }

    BLST_FORCEINLINE bool_t equal(const BoundedElement &other) const {
        return data.equal(other.data);
    }

    // Clearer API: replace self with replacement if condition is true
    // ctime_select(selection, b) returns *this when selection is true, b when false
    // So we need to negate: if condition is true, we want replacement (not *this)
    BLST_FORCEINLINE auto ctime_replace_if(bool condition, const BoundedElement &replacement) const {
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

    BLST_FORCEINLINE DelayedProduct(Left left, Right right) : left(left), right(right) {}

    BLST_FORCEINLINE auto complete() const {
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

    // The low half accumulates 52-bit IFMA contributions regardless of the
    // input RNS radix. Keep it in the signed 64-bit envelope used by the
    // representation and require a nonnegative, nondegenerate interval.
    static constexpr int64_t max_wide_bound = int64_t{1} << (64 - 52 - 1);
    static_assert(BoundsType::lower >= 0,
        "Wide accumulators require nonnegative bounds");
    static_assert(BoundsType::upper < max_wide_bound,
        "Wide accumulator bound exceeds the signed 52-bit contribution capacity");
    static_assert((BoundsType::lower == 0 && BoundsType::upper == 0) ||
        BoundsType::lower < BoundsType::upper,
        "Only the zero wide accumulator may have equal lower and upper bounds");

    // As with BoundedElement, default initialization is overwrite-before-read.
    BLST_FORCEINLINE WideElement() = default;
    BLST_FORCEINLINE explicit WideElement(const ElementType &high, const ElementType &low) : high(high), low(low) {}
    BLST_FORCEINLINE explicit WideElement(const AVXVector<limbs> &high_data, const AVXVector<limbs> &low_data) : high(high_data), low(low_data) {}
    
    template<class OtherBounds>
    BLST_FORCEINLINE auto operator+(const WideElement<OtherBounds, limbs, element_bits> &other) const {
        auto result_low = low + other.low;
        auto result_high = high + other.high;
        return WideElement<decltype(bounds + other.bounds), limbs, element_bits>(result_high, result_low);
    }

    template<class OtherBounds>
    BLST_FORCEINLINE auto operator-(const WideElement<OtherBounds, limbs, element_bits> &other) const {
        // Wide halves are deliberately redundant rather than normalized. The
        // low halves can therefore be subtracted without propagating a borrow
        // only when the left interval dominates the complete right interval.
        static_assert(BoundsType::lower >= OtherBounds::upper,
            "Wide subtraction requires a borrow-free low accumulator");
        auto result_low = low - other.low;
        auto result_high = high - other.high;
        return WideElement<decltype(bounds - other.bounds), limbs, element_bits>(result_high, result_low);
    }

    template<int64_t left_shift>
    BLST_FORCEINLINE auto shift_left() const {
        auto high_shifted = high.template shift_left<left_shift>();
        auto low_shifted = low.template shift_left<left_shift>();
        return WideElement<decltype(bounds.template shift_left<left_shift>()), limbs, element_bits>(high_shifted, low_shifted);
    }

    template<class LeftBounds, class RightBounds>
    BLST_FORCEINLINE auto mulacc(const BoundedElement<LeftBounds, limbs, element_bits> &left, const BoundedElement<RightBounds, limbs, element_bits> &right) const {
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
    BLST_FORCEINLINE auto operator+(const DelayedProduct<BoundedElement<LeftBounds, limbs, element_bits>, BoundedElement<RightBounds, limbs, element_bits>> &delayed_product) const {
        return mulacc(delayed_product.left, delayed_product.right);
    }

    // No action required for non delayed products.
    BLST_FORCEINLINE WideElement complete() const {
        return *this;
    }

    template<int64_t new_lower, int64_t new_upper>
    BLST_FORCEINLINE auto recast() const {
        auto new_bounds = bounds.template recast<new_lower, new_upper>();
        return WideElement<decltype(new_bounds), limbs, element_bits>(high.template recast<new_lower, new_upper>(), low.template recast<new_lower, new_upper>());
    }

};

template<int limbs, int element_bits>
BLST_FORCEINLINE auto WideZero() {
    auto zero = BoundedElement<Bounds<0, 0>, limbs, element_bits>::zero();
    return WideElement<Bounds<0, 0>, limbs, element_bits>(zero, zero);

}

}  // namespace blst::avx
