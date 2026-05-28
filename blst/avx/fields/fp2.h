#pragma once
#include "../common/preprocess.h"
#include "../rns/ring_element.h"
#include <tuple>
#include <array>
#include <type_traits>
// Forward declaration
namespace blst::avx {

template<class X, class Y, class minusY>
class FP2N;

template<class X, class Y>
class FP2 {
    private:
    X x_;
    Y y_;

    public:
    // Default constructor - default-constructs members
    INLINE FP2() : x_(), y_() {}
    
    INLINE FP2(const X &x, const Y &y) : x_(x), y_(y) {}

    // Const getters to access members
    INLINE const X& x() const { return x_; }
    INLINE const Y& y() const { return y_; }

    template<class OtherX, class OtherY>
    INLINE auto operator+(const FP2<OtherX, OtherY> &other) const {
        auto x_sum = x_ + other.x();
        auto y_sum = y_ + other.y();
        return FP2<decltype(x_sum), decltype(y_sum)>(x_sum, y_sum);
    }

    template<class OtherX, class OtherY>
    INLINE auto operator-(const FP2<OtherX, OtherY> &other) const {
        auto x_diff = x_ - other.x();
        auto y_diff = y_ - other.y();
        return FP2<decltype(x_diff), decltype(y_diff)>(x_diff, y_diff);
    }

    template<class Ring>
    INLINE auto prep(const Ring &ring) const {
        auto x_resolved = ring.prep(x_);
        auto y_resolved = ring.prep(y_);
        return FP2<decltype(x_resolved), decltype(y_resolved)>(x_resolved, y_resolved);
    }

    // And only the negated prepared ones can do products?
    // Basically we can reuse the negation result a lot.
    template<class Ring>
    INLINE auto prep_left(const Ring &ring) const {
        auto r_x = ring.prep(x_);
        auto r_y = ring.prep(y_);
        auto minus_y = ring.prep(ring.negate(r_y));
        return FP2N<decltype(r_x), decltype(r_y), decltype(minus_y)>(r_x, r_y, minus_y);
    }

    template<class OtherXN, class OtherYN, class OtherMinusYN, class OtherX, class OtherY>
    INLINE auto mul_accumulate(const FP2N<OtherXN, OtherYN, OtherMinusYN> &other1, const FP2<OtherX, OtherY> &other2) const {
        auto result_x = x_ + other1.x() * other2.x() + other1.minus_y() * other2.y();
        auto result_y = y_ + other1.x() * other2.y() + other1.y() * other2.x();
        return FP2<decltype(result_x), decltype(result_y)>(result_x, result_y);
    }

    template<class Ring, class OtherXN, class OtherYN, class OtherMinusYN, class OtherX, class OtherY>
    INLINE auto mul_sub_accumulate(FP2N<OtherXN, OtherYN, OtherMinusYN> &other1, const FP2<OtherX, OtherY> &other2, const Ring &ring) const {
        auto minus_x = ring.prep(ring.negate(other1.x()));
        auto result_x = x_ + minus_x * other2.x() + other1.y() * other2.y();
        auto result_y = y_ + minus_x * other2.y() + other1.minus_y() * other2.x();
        return FP2<decltype(result_x), decltype(result_y)>(result_x, result_y);
    }

    template<class FP>
    INLINE auto mul_scalar(const FP &scalar) const {
        auto result_x = x_ * scalar;
        auto result_y = y_ * scalar;
        return FP2<decltype(result_x), decltype(result_y)>(result_x, result_y);
    }

    // Multiply FP2 by a base field element (scalar multiplication in Fp)
    template<class FP>
    INLINE auto mul_fp(const FP &other) const {
       return mul_scalar(other);
    }

    template<int N>
    INLINE auto operator*(std::integral_constant<int, N>) const {
        return mul_scalar(std::integral_constant<int, N>{});
    }

    INLINE auto ctime_select(const bool condition, const FP2<X, Y> &other) const {
        return FP2<X, Y>(x_.ctime_select(condition, other.x()), y_.ctime_select(condition, other.y()));
    }

    // Clearer API: replace self with replacement if condition is true
    // ctime_select(selection, b) returns *this when selection is true, b when false
    // So we need to negate: if condition is true, we want replacement (not *this)
    INLINE auto ctime_replace_if(bool condition, const FP2<X, Y> &replacement) const {
        return ctime_select(!condition, replacement);
    }

    // Clearer API: replace self with replacement if condition is false
    INLINE auto ctime_replace_unless(bool condition, const FP2<X, Y> &replacement) const {
        return ctime_select(condition, replacement);
    }
};

template<class X, class Y, class minusY>
class FP2N {
    private:
    X x_;
    Y y_;
    minusY minus_y_;

    public:
    INLINE FP2N(const X &x, const Y &y, const minusY &minus_y) : x_(x), y_(y), minus_y_(minus_y) {}

    // Converting constructor from another FP2N type
    template<class OtherX, class OtherY, class OtherMinusY>
    INLINE FP2N(const FP2N<OtherX, OtherY, OtherMinusY> &other) : x_(other.x()), y_(other.y()), minus_y_(other.minus_y()) {}

    // Const getters to access members
    INLINE const X& x() const { return x_; }
    INLINE const Y& y() const { return y_; }
    INLINE const minusY& minus_y() const { return minus_y_; }

    template<class OtherX, class OtherY>
    INLINE auto operator*(const FP2<OtherX, OtherY> &other) const {
        auto result_x = x_ * other.x() + minus_y_ * other.y();
        auto result_y = x_ * other.y() + y_ * other.x();
        return FP2<decltype(result_x), decltype(result_y)>(result_x, result_y);
    }

    INLINE FP2<X, Y> to_fp2() const {
        return FP2<X, Y>(x_, y_);
    }

    template<class FP>
    INLINE auto mul_fp(const FP &scalar) const {
        auto fp2 = to_fp2();
        return fp2.mul_fp(scalar);
    }
    
    // Conversion operator
    INLINE operator FP2<X, Y>() const {
        return FP2<X, Y>(x_, y_);
    }

    INLINE auto square() const {
        auto result_x = x_ * x_ + minus_y_ * y_;
        auto result_y = x_ * y_ * std::integral_constant<int, 2>{};
        return FP2<decltype(result_x), decltype(result_y)>(result_x, result_y);
    }
};

template<class Ring>
class FP2Ring {
    const Ring &ring;

    public:
    constexpr static auto _12 = std::integral_constant<int, 12>{};
    using StandardElement = FP2<typename Ring::StandardElement, typename Ring::StandardElement>;

    // May need to remove if it is recalling the constructor for the ring.
    INLINE FP2Ring(const Ring &r) : ring(r) {}


    // Helper to convert variadic arguments to array using index_sequence
    // We extract two ring elements per FP2 element, which makes the total size 2*
    // Layout: [x0, x1, ..., xN, y0, y1, ..., yN] (strided)
    template<typename... Args, std::size_t... I>
    INLINE auto make_ready_array_impl(const std::tuple<Args...> &args, std::index_sequence<I...>) const {
         return std::array<typename Ring::template ReadyElement<Ring::MAX_ADD>, 2*sizeof...(Args)>{ring.template ready<Ring::MAX_ADD>(std::get<I>(args).x())..., ring.template ready<Ring::MAX_ADD>(std::get<I>(args).y())...};
     }

    // Helper to convert variadic arguments to array for batch_reduce (uses ready<1>)
    template<typename... Args, std::size_t... I>
    INLINE auto make_ready_array_impl_reduce(const std::tuple<Args...> &args, std::index_sequence<I...>) const {
        constexpr int64_t reduce_bound = 2 * Ring::MULT_OK_BOUND * Ring::MULT_OK_BOUND + 1;
        return std::array<typename Ring::template ReadyElement<reduce_bound>, 2*sizeof...(Args)>{ring.template ready<reduce_bound>(std::get<I>(args).x())..., ring.template ready<reduce_bound>(std::get<I>(args).y())...};
     }

    // Helper to extract elements from SmallElement FP2 types for batch_expand
    template<typename... Args, std::size_t... I>
    INLINE auto make_element_array_impl(const std::tuple<Args...> &args, std::index_sequence<I...>) const {
        return std::array<typename Ring::StandardBoundedElement, 2*sizeof...(Args)>{std::get<I>(args).x().element..., std::get<I>(args).y().element...};
    }

    // Variadic template version that captures all the different arities
    template<typename... Args>
    INLINE auto batch_reduce_expand(const Args&... args) const {
        constexpr int batch = static_cast<int>(sizeof...(Args));
        static_assert(batch > 0, "batch_reduce_expand requires at least one argument");
        auto tuple = std::tie(args...);
        auto ready_array = make_ready_array_impl(tuple, std::make_index_sequence<sizeof...(Args)>{});
        // We have 2*#FP2 elements FP elements in strided layout: [x0...xN, y0...yN]
        auto reduced = ring.template batch_reduce_expand<2*batch>(ready_array);
        // Reconstruct FP2 elements by pairing x[i] with y[i] from the strided array
        std::array<StandardElement, batch> result;
        for (int i = 0; i < batch; i++) {
            result[i] = StandardElement(reduced[i], reduced[i + batch]);
        }
        return result;
    }

    // Variadic template version for batch_reduce
    template<typename... Args>
    INLINE auto batch_reduce(const Args&... args) const {
        constexpr int batch = static_cast<int>(sizeof...(Args));
        static_assert(batch > 0, "batch_reduce requires at least one argument");
        auto tuple = std::tie(args...);
        auto ready_array = make_ready_array_impl_reduce(tuple, std::make_index_sequence<sizeof...(Args)>{});
        // We have 2*#FP2 elements FP elements in strided layout: [x0...xN, y0...yN]
        auto reduced = ring.template batch_reduce<2*batch>(ready_array);
        // Reconstruct FP2 elements by pairing x[i] with y[i] from the strided array
        using ReducedXType = std::remove_reference_t<decltype(reduced[0])>;
        using ReducedYType = std::remove_reference_t<decltype(reduced[batch])>;
        using ReducedFP2Type = FP2<ReducedXType, ReducedYType>;
        std::array<ReducedFP2Type, batch> result;
        for (int i = 0; i < batch; i++) {
            result[i] = ReducedFP2Type(ReducedXType(reduced[i]), ReducedYType(reduced[i + batch]));
        }
        return result;
    }

    template<typename... Args, std::size_t... I>
    INLINE auto make_expand_array_impl(const std::tuple<Args...>& args, std::index_sequence<I...>) const {
        return std::array<typename Ring::StandardBoundedElement, 2 * sizeof...(Args)>{
            std::get<I>(args).x().element..., std::get<I>(args).y().element...};
    }

    template<class E, class ExpandedX, class ExpandedY>
    INLINE auto reattach_expand_bounds(const E& e, const ExpandedX& expanded_x, const ExpandedY& expanded_y) const {
        using OutX = ExpandedElement<typename Ring::StandardBoundedElement, std::remove_cv_t<decltype(e.x().rns_bounds)>>;
        using OutY = ExpandedElement<typename Ring::StandardBoundedElement, std::remove_cv_t<decltype(e.y().rns_bounds)>>;
        return FP2<OutX, OutY>(OutX(expanded_x, e.x()), OutY(expanded_y, e.y()));
    }

    template<typename... Args, std::size_t... I>
    INLINE auto reattach_expand_tuple_impl(
        const std::tuple<Args...>& args,
        const std::array<typename Ring::StandardBoundedElement, 2 * sizeof...(Args)>& expanded,
        std::index_sequence<I...>) const {
        constexpr size_t batch = sizeof...(Args);
        return std::tuple<decltype(reattach_expand_bounds(std::get<I>(args), expanded[I], expanded[I + batch]))...>{
            reattach_expand_bounds(std::get<I>(args), expanded[I], expanded[I + batch])...};
    }

    template<typename... Args>
    INLINE auto batch_expand(const Args&... args) const {
        constexpr int batch = static_cast<int>(sizeof...(Args));
        static_assert(batch > 0, "batch_expand requires at least one argument");
        const auto tuple = std::tie(args...);
        const auto data = make_expand_array_impl(tuple, std::make_index_sequence<sizeof...(Args)>{});
        const auto expanded = ring.template batch_expand<2 * batch>(data);
        return reattach_expand_tuple_impl(tuple, expanded, std::make_index_sequence<sizeof...(Args)>{});
    }

    template<class X, class Y>
    INLINE auto small_m2(const FP2<X, Y>& fp2) const {
        return FP2<decltype(ring.small_m2(fp2.x())), decltype(ring.small_m2(fp2.y()))>(
            ring.small_m2(fp2.x()),
            ring.small_m2(fp2.y()));
    }

    template<int upper_bound, class X, class Y>
    INLINE auto recast(const FP2<X, Y> &fp2) const {
        auto x_recast = ring.template recast<0, upper_bound>(fp2.x());
        auto y_recast = ring.template recast<0, upper_bound>(fp2.y());
        return FP2<decltype(x_recast), decltype(y_recast)>(x_recast, y_recast);
    }

    static INLINE StandardElement zero(const Ring &ring) {
        auto zero_val = ring.zero();
        return FP2<decltype(zero_val), decltype(zero_val)>(zero_val, zero_val);
    }

    INLINE StandardElement zero() const {
        auto zero_val = ring.zero();
        return FP2<decltype(zero_val), decltype(zero_val)>(zero_val, zero_val);
    }

    static INLINE StandardElement one(const Ring &ring) {
        auto one = ring.one();
        auto zero = ring.zero();
        return FP2<decltype(one), decltype(zero)>(one, zero);
    }

    INLINE StandardElement one() const {
        auto one_val = ring.one();
        auto zero_val = ring.zero();
        return FP2<decltype(one_val), decltype(zero_val)>(one_val, zero_val);
    }

    template<class FP2Type>
    INLINE auto mul_by_nonresidue(const FP2Type &fp2) const {
        // Multiply by (u + 1)
        auto result_x = fp2.x() + ring.negate(fp2.y());
        auto result_y = fp2.x() + fp2.y();
        return FP2<decltype(result_x), decltype(result_y)>(result_x, result_y);
    }

    template<class FP2Type>
    INLINE auto mul_3b(const FP2Type &fp2) const {
        // BLS12-381: b = 4(1+u), 3b = 12(1+u)
        return mul_by_nonresidue(fp2) * _12;
    }

    INLINE auto negative_range_offset() const {
        auto offset = ring.negative_range_offset();
        return FP2<decltype(offset), decltype(offset)>(offset, offset);
    }

    INLINE auto mul(const StandardElement &a, const StandardElement &b) const {
        return batch_reduce_expand(a.prep_left(ring) * b.prep(ring));
    }

    INLINE auto square(const StandardElement &a) const {
        return batch_reduce_expand(a.prep_left(ring).square());
    }

    template<class X, class Y>
    INLINE auto prep(const FP2<X, Y> &fp2) const {
        return fp2.prep(ring);
    }

    template<class X, class Y>
    INLINE auto prep_left(const FP2<X, Y> &fp2) const {
        return fp2.prep_left(ring);
    }

    template<class X, class Y>
    INLINE auto standard_negate(const FP2<X, Y> &fp2) const {
        return FP2<decltype(ring.standard_negate(fp2.x())), decltype(ring.standard_negate(fp2.y()))>(ring.standard_negate(fp2.x()), ring.standard_negate(fp2.y()));
    }

    template<class X, class Y>
    INLINE auto negate(const FP2<X, Y> &fp2) const {
        return FP2<decltype(ring.negate(fp2.x())), decltype(ring.negate(fp2.y()))>(ring.negate(fp2.x()), ring.negate(fp2.y()));
    }
};

}  // namespace blst::avx
