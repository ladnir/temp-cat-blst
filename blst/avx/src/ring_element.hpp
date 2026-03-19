#pragma once
#include "bounded_element.hpp"
#include "preprocess.hpp"


template<class Element, class RNSBounds>
class SmallElement {
    public:
    Element element;
    static constexpr RNSBounds rns_bounds = RNSBounds{};

    INLINE SmallElement() = default;
    INLINE explicit SmallElement(const Element &e) : element(e) {}

    template<class OtherElement, class OtherRNSBounds>
    INLINE auto operator+(const SmallElement<OtherElement, OtherRNSBounds> &other) const {
        auto result = element + other.element;
        return SmallElement<decltype(result), decltype(rns_bounds + other.rns_bounds)>(result);
    }

    template<class OtherElement, class OtherRNSBounds>
    INLINE auto operator-(const SmallElement<OtherElement, OtherRNSBounds> &other) const {
        auto result = element - other.element;
        return SmallElement<decltype(result), decltype(rns_bounds - other.rns_bounds)>(result);
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

    template<int64_t left_shift>
    INLINE auto shift_left() const {
        auto element_shifted = element.template shift_left<left_shift>();
        auto bounds_shifted = rns_bounds.template shift_left<left_shift>();
        return SmallElement<decltype(element_shifted), decltype(bounds_shifted)>(element_shifted);
    }

    INLINE auto operator*(std::integral_constant<int, 12>) const {
        return (*this * std::integral_constant<int, 3>()).template shift_left<2>();
    }

    template<int new_lower, int new_upper>
    INLINE auto recast_rns_bounds() const {
        // Fails to compile if new bounds are invalid.
        auto new_rns_bounds = rns_bounds.template recast<new_lower, new_upper>();
        return SmallElement<decltype(element), decltype(new_rns_bounds)>(element);
    }
};

template<class Element, class RNSBounds>
class ExpandedElement {
    public:
    Element m1;
    Element m2;

    static constexpr RNSBounds rns_bounds = RNSBounds{};

    INLINE ExpandedElement() = default;
    INLINE explicit ExpandedElement(const Element &m1, const Element &m2) : m1(m1), m2(m2) {}

    INLINE explicit ExpandedElement(const Element &m1, const SmallElement<Element, RNSBounds> &m2) : m1(m1), m2(m2.element) {}

    template<class OtherElement, class OtherRNSBounds>
    INLINE auto operator+(const ExpandedElement<OtherElement, OtherRNSBounds> &other) const {
        auto result_m1 = m1 + other.m1;
        auto result_m2 = m2 + other.m2;
        return ExpandedElement<decltype(result_m1), decltype(rns_bounds + other.rns_bounds)>(result_m1, result_m2);
    }

    template<class OtherElement, class OtherRNSBounds>
    INLINE auto operator-(const ExpandedElement<OtherElement, OtherRNSBounds> &other) const {
        auto result_m1 = m1 - other.m1;
        auto result_m2 = m2 - other.m2;
        return ExpandedElement<decltype(result_m1), decltype(rns_bounds - other.rns_bounds)>(result_m1, result_m2);
    }

    template<class OtherElement, class OtherRNSBounds>
    INLINE auto operator*(const ExpandedElement<OtherElement, OtherRNSBounds> &other) const {
        using DelayedProductType = DelayedProduct<Element, OtherElement>;
        using ResultBoundsType = decltype(rns_bounds * other.rns_bounds);
        return ExpandedElement<DelayedProductType, ResultBoundsType>(
            DelayedProductType(m1, other.m1), 
            DelayedProductType(m2, other.m2)
        );
    }

    template<class LeftElement, class RightElement, class OtherRNSBounds>
    INLINE auto operator+(const ExpandedElement<DelayedProduct<LeftElement, RightElement>, OtherRNSBounds> &other) const {
        auto x = complete();
        auto result_m1 = x.m1.mulacc(other.m1.left, other.m1.right);
        auto result_m2 = x.m2.mulacc(other.m2.left, other.m2.right);
        using ResultBoundsType = decltype(x.rns_bounds + other.rns_bounds);
        return ExpandedElement<decltype(result_m1), ResultBoundsType>(result_m1, result_m2);
    }
    
    template<int64_t left_shift>
    INLINE auto shift_left() const {
        auto m1_shifted = m1.template shift_left<left_shift>();
        auto m2_shifted = m2.template shift_left<left_shift>();
        auto bounds_shifted = rns_bounds.template shift_left<left_shift>();
        return ExpandedElement<decltype(m1_shifted), decltype(bounds_shifted)>(m1_shifted, m2_shifted);
    }

    INLINE auto complete() const {
        auto result_m1 = m1.complete();
        auto result_m2 = m2.complete();
        return ExpandedElement<decltype(result_m1), decltype(rns_bounds)>(result_m1, result_m2);
    }

    // Convenience operators
    // Beyond 3 you may want to do true scalar multiplication (mullo, mulhi)
    INLINE auto operator*(std::integral_constant<int, 2>) const {
        auto x = complete();
        return x + x;
    }

    INLINE auto operator*(std::integral_constant<int, 3>) const {
        auto x = complete();
        return x + x + x;
    }

    INLINE auto operator*(std::integral_constant<int, 6>) const {
        auto x = complete();
        return (x * std::integral_constant<int, 3>{}) * std::integral_constant<int, 2>{};
    }

    INLINE auto operator*(std::integral_constant<int, 8>) const {
        auto x = complete();
        return x.template shift_left<3>();
    }

    INLINE auto operator*(std::integral_constant<int, 12>) const {
        return (*this * std::integral_constant<int, 3>()).template shift_left<2>();
    }

    template<int64_t new_lower, int64_t new_upper>
    INLINE auto recast_rns_bounds() const {
        auto new_rns_bounds = rns_bounds.template recast<new_lower, new_upper>();
        return ExpandedElement<decltype(m1), decltype(new_rns_bounds)>(m1, m2);
    }

    INLINE auto ctime_select(bool selection, const ExpandedElement &other) const {
        return ExpandedElement(m1.ctime_select(selection, other.m1), m2.ctime_select(selection, other.m2));
    }

    // Clearer API: replace self with replacement if condition is true
    // ctime_select(selection, b) returns *this when selection is true, b when false
    // So we need to negate: if condition is true, we want replacement (not *this)
    INLINE auto ctime_replace_if(bool condition, const ExpandedElement &replacement) const {
        return ctime_select(!condition, replacement);
    }
};

template<int limbs, int element_bits=52>
INLINE auto from_constant(const std::array<uint64_t, limbs> &array_m1, const std::array<uint64_t, limbs> &array_m2) {
    auto data_m1 = AVXVector<limbs>(array_m1);
    auto data_m2 = AVXVector<limbs>(array_m2);
    using BoundsType = Bounds<0, 1>;
    using ElementType = BoundedElement<BoundsType, limbs, element_bits>;
    return ExpandedElement<ElementType, BoundsType>(
        ElementType(data_m1), 
        ElementType(data_m2)
    );
}

template<int limbs, int element_bits=52>
INLINE auto ElementWideZero() {
    auto wide_zero = WideZero<limbs, element_bits>();
    return ExpandedElement<decltype(wide_zero), Bounds<0, 0>>(wide_zero, wide_zero);
}