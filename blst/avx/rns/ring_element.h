// SPDX-License-Identifier: MIT
// Copyright 2026 Simon Langowski
// Copyright 2026 Category Labs
//
// Adapted from the VROOM RNS carrier at commit d6ecd534. CatBLST retains the
// compile-time bound model while giving its selected operations explicit
// contracts below.
#pragma once

#include "bounded_element.h"
#include "../../common/defines.h"

namespace blst::avx {

// One-basis carrier. RNSBounds tracks the represented multiple of the field
// modulus; it is compile-time proof metadata, not runtime validation.
template<class Element, class RNSBounds>
class SmallElement {
    public:
    using element_type = Element;
    using rns_bounds_type = RNSBounds;

    Element element;
    static constexpr RNSBounds rns_bounds = RNSBounds{};

    // Default initialization is overwrite-before-read. Value initialization
    // (`SmallElement{}`) zero-initializes the selected trivial element types.
    BLST_FORCEINLINE SmallElement() = default;
    BLST_FORCEINLINE explicit SmallElement(const Element &e) : element(e) {}

    template<class OtherElement, class OtherRNSBounds>
    BLST_FORCEINLINE auto operator+(const SmallElement<OtherElement, OtherRNSBounds> &other) const {
        auto result = element + other.element;
        return SmallElement<decltype(result), decltype(rns_bounds + other.rns_bounds)>(result);
    }

    template<class OtherElement, class OtherRNSBounds>
    BLST_FORCEINLINE auto operator-(const SmallElement<OtherElement, OtherRNSBounds> &other) const {
        auto result = element - other.element;
        return SmallElement<decltype(result), decltype(rns_bounds - other.rns_bounds)>(result);
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

    template<int64_t left_shift>
    BLST_FORCEINLINE auto shift_left() const {
        auto element_shifted = element.template shift_left<left_shift>();
        using ShiftedRNSBounds = decltype(
            RNSBounds{}.template shift_left<left_shift>());
        return SmallElement<decltype(element_shifted), ShiftedRNSBounds>(
            element_shifted);
    }

    BLST_FORCEINLINE auto operator*(std::integral_constant<int, 12>) const {
        return (*this * std::integral_constant<int, 3>()).template shift_left<2>();
    }

    template<int64_t new_lower, int64_t new_upper>
    BLST_FORCEINLINE auto recast_rns_bounds() const {
        using NewRNSBounds = decltype(
            RNSBounds{}.template recast<new_lower, new_upper>());
        return SmallElement<Element, NewRNSBounds>(element);
    }

    // Exact-representation predicates only. They are not semantic field
    // observations when the RNS representation permits redundant zero lifts.
    BLST_FORCEINLINE bool_t is_zero() const {
        return element.is_zero();
    }

    BLST_FORCEINLINE bool_t equal(const SmallElement &other) const {
        return element.equal(other.element);
    }
};

// Two-basis carrier. m1 and m2 must encode the same field value, and
// RNSBounds must describe that shared value. Constructors do not establish or
// validate either invariant; only reviewed producers may attach these claims.
template<class Element, class RNSBounds>
class ExpandedElement {
    public:
    using element_type = Element;
    using rns_bounds_type = RNSBounds;

    Element m1;
    Element m2;

    static constexpr RNSBounds rns_bounds = RNSBounds{};

    // Default initialization is overwrite-before-read. Value initialization
    // (`ExpandedElement{}`) zero-initializes the selected trivial element types.
    BLST_FORCEINLINE ExpandedElement() = default;
    constexpr BLST_FORCEINLINE explicit ExpandedElement(const Element &m1, const Element &m2) : m1(m1), m2(m2) {}

    BLST_FORCEINLINE explicit ExpandedElement(const Element &m1, const SmallElement<Element, RNSBounds> &m2) : m1(m1), m2(m2.element) {}

    template<class OtherElement, class OtherRNSBounds>
    BLST_FORCEINLINE auto operator+(const ExpandedElement<OtherElement, OtherRNSBounds> &other) const {
        auto result_m1 = m1 + other.m1;
        auto result_m2 = m2 + other.m2;
        return ExpandedElement<decltype(result_m1), decltype(rns_bounds + other.rns_bounds)>(result_m1, result_m2);
    }

    template<class OtherElement, class OtherRNSBounds>
    BLST_FORCEINLINE auto operator-(const ExpandedElement<OtherElement, OtherRNSBounds> &other) const {
        auto result_m1 = m1 - other.m1;
        auto result_m2 = m2 - other.m2;
        return ExpandedElement<decltype(result_m1), decltype(rns_bounds - other.rns_bounds)>(result_m1, result_m2);
    }

    template<class OtherElement, class OtherRNSBounds>
    BLST_FORCEINLINE auto operator*(const ExpandedElement<OtherElement, OtherRNSBounds> &other) const {
        using DelayedProductType = DelayedProduct<Element, OtherElement>;
        using ResultBoundsType = decltype(rns_bounds * other.rns_bounds);
        return ExpandedElement<DelayedProductType, ResultBoundsType>(
            DelayedProductType(m1, other.m1), 
            DelayedProductType(m2, other.m2)
        );
    }

    template<class LeftElement, class RightElement, class OtherRNSBounds>
    BLST_FORCEINLINE auto operator+(const ExpandedElement<DelayedProduct<LeftElement, RightElement>, OtherRNSBounds> &other) const {
        // Delayed products are accumulated into a completed wide value. This
        // overload is intentionally ordered: the left operand supplies the
        // existing accumulator and the right operand supplies one product.
        auto x = complete();
        auto result_m1 = x.m1.mulacc(other.m1.left, other.m1.right);
        auto result_m2 = x.m2.mulacc(other.m2.left, other.m2.right);
        using ResultBoundsType = decltype(x.rns_bounds + other.rns_bounds);
        return ExpandedElement<decltype(result_m1), ResultBoundsType>(result_m1, result_m2);
    }
    
    template<int64_t left_shift>
    BLST_FORCEINLINE auto shift_left() const {
        auto m1_shifted = m1.template shift_left<left_shift>();
        auto m2_shifted = m2.template shift_left<left_shift>();
        using ShiftedRNSBounds = decltype(
            RNSBounds{}.template shift_left<left_shift>());
        return ExpandedElement<decltype(m1_shifted), ShiftedRNSBounds>(
            m1_shifted, m2_shifted);
    }

    BLST_FORCEINLINE auto complete() const {
        auto result_m1 = m1.complete();
        auto result_m2 = m2.complete();
        return ExpandedElement<decltype(result_m1), RNSBounds>(
            result_m1, result_m2);
    }

    // Fixed compile-time scalings used by the selected curve and field-tower
    // formulas. These are additions and shifts, not field multiplication.
    BLST_FORCEINLINE auto operator*(std::integral_constant<int, 2>) const {
        auto x = complete();
        return x + x;
    }

    BLST_FORCEINLINE auto operator*(std::integral_constant<int, 3>) const {
        auto x = complete();
        return x + x + x;
    }

    BLST_FORCEINLINE auto operator*(std::integral_constant<int, 6>) const {
        auto x = complete();
        return (x * std::integral_constant<int, 3>{}) * std::integral_constant<int, 2>{};
    }

    BLST_FORCEINLINE auto operator*(std::integral_constant<int, 8>) const {
        auto x = complete();
        return x.template shift_left<3>();
    }

    BLST_FORCEINLINE auto operator*(std::integral_constant<int, 12>) const {
        return (*this * std::integral_constant<int, 3>()).template shift_left<2>();
    }

    template<int64_t new_lower, int64_t new_upper>
    BLST_FORCEINLINE auto recast_rns_bounds() const {
        using NewRNSBounds = decltype(
            RNSBounds{}.template recast<new_lower, new_upper>());
        return ExpandedElement<Element, NewRNSBounds>(m1, m2);
    }

    // Selects between values with identical physical and proof-metadata types.
    BLST_FORCEINLINE auto ctime_select(bool selection, const ExpandedElement &other) const {
        return ExpandedElement(m1.ctime_select(selection, other.m1), m2.ctime_select(selection, other.m2));
    }

    // Exact-representation predicates only. They are not semantic field
    // observations when the RNS representation permits redundant zero lifts.
    BLST_FORCEINLINE bool_t is_zero() const {
        return m1.is_zero() & m2.is_zero();
    }

    BLST_FORCEINLINE bool_t equal(const ExpandedElement &other) const {
        return m1.equal(other.m1) & m2.equal(other.m2);
    }

    // Clearer API: replace self with replacement if condition is true
    // ctime_select(selection, b) returns *this when selection is true, b when false
    // So we need to negate: if condition is true, we want replacement (not *this)
    BLST_FORCEINLINE auto ctime_replace_if(bool condition, const ExpandedElement &replacement) const {
        return ctime_select(!condition, replacement);
    }
};

}  // namespace blst::avx
