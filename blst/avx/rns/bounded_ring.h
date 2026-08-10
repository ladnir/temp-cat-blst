#pragma once
#include "ring_element.h"
#include "elementwise_reduce.h"
#include "change_base.h"
#include "../cpu/vector/multiplication.hpp"
#include "../runtime_constants.h"
#include "../../common/bytes.h"
#include "../../fields/Fp.h"
#include <array>
#include <tuple>
#include <utility>
#include <cstddef>
#include <type_traits>
#include <string>
#include <iostream>
#include "../common/preprocess.h"


// The maximum negative and positive values that are needed.  These occur in the Y coordinate of doubling_step.
namespace blst::avx {

template<int modulus_bits=381, int limbs=8, int element_bits=52, int RNS_MAX_NEG=-1932, int RNS_MAX_POS=2377, int LOG_MULTIPLES=4>
class BoundedRing {
    
    // -1 to account for each moduli being slighly less than element_bits
    static constexpr int RNS_BITS = element_bits * limbs - 1;

    // We store values only in the positive range, so we need a range of [0, MAX_POS - MAX_NEG]
    static constexpr int RNS_MAX_MULTIPLE =  (RNS_MAX_POS - RNS_MAX_NEG);
    static_assert(RNS_MAX_NEG <= 0, "RNS_MAX_NEG must be negative");
    static_assert(RNS_MAX_POS > 0, "RNS_MAX_POS must be positive");
    
    // We use a redundnace of 3, so we multiply by 9 (for products) to get the total range.
    static constexpr int RNS_REDUNDANT_MULTIPLE = 9 * RNS_MAX_MULTIPLE;
    static constexpr int required_redundance_bits = ceil_log2(RNS_REDUNDANT_MULTIPLE);
    static_assert(modulus_bits + required_redundance_bits <= RNS_BITS, "Modulus bits + required redundance must be less than or equal to RNS bits (conservative bound)");

    public:
    
    // In theory we could pass these reducers to the IntRNS2System constructor (also would fix the hardcoded 52 bit operations for the IntRNS2System)
    const IntRNS2System<modulus_bits, limbs> sys;
    const MontgomeryReduce2<limbs, element_bits, LOG_MULTIPLES> montgomery_reduce_m1;
    const MontgomeryReduce2<limbs, element_bits, LOG_MULTIPLES> montgomery_reduce_m2;

    // A 52-bit or less element, which can be used directly in 52-bit multiplication operations.
    constexpr static uint64_t MULT_OK_BOUND = 1ull << (52 - element_bits);
    // The output of the change base operation; if possible, we skip the reduction and output with a redundancy of 2
    constexpr static uint64_t STANDARD_BOUND = (MULT_OK_BOUND < 2ull) ? MULT_OK_BOUND : 2ull;
    // To use array types, we recast everything to some maximum number of additions.  You can probably make this bigger if needed.
    static constexpr uint64_t MAX_ADD = 800;
    using FullyReducedElement = BoundedElement<Bounds<0, 1>, limbs, element_bits>;
    using StandardBoundedElement = BoundedElement<Bounds<0, STANDARD_BOUND>, limbs, element_bits>;
    using MultOkBoundedElement = BoundedElement<Bounds<0, MULT_OK_BOUND>, limbs, element_bits>;
    // A fully reduced element, the output of RNS reduction.
    using StandardElement = ExpandedElement<StandardBoundedElement, Bounds<0, 1>>;
    template<int64_t A>
    using ReadyElement = ReadyToReduce<WideElement<Bounds<0, A>, limbs, element_bits>, MultOkBoundedElement>;

    template<int64_t rns_bound>
    using OffsetElement = ExpandedElement<FullyReducedElement, Bounds<rns_bound, rns_bound>>;

    template<int64_t ele_bound, int64_t rns_bound>
    using WideOffsetElement = ExpandedElement<WideElement<Bounds<ele_bound, ele_bound+1>, limbs, element_bits>, Bounds<rns_bound, rns_bound>>;

    template<int64_t upper_bound>
    using ReducedElement = ExpandedElement<StandardBoundedElement, Bounds<0, upper_bound>>;

    template<int64_t upper_bound>
    using SmallReducedElement = SmallElement<StandardBoundedElement, Bounds<0, upper_bound>>;

    using SmallStandardElement = SmallReducedElement<1>;

    using RawReduceOutput = SmallElement<Bounds<0, MAX_ADD + 2*limbs*MULT_OK_BOUND>, Bounds<0, 1>>;

    private:

    static INLINE StandardElement make_standard_element(
        const std::array<uint64_t, limbs>& m1,
        const std::array<uint64_t, limbs>& m2)
    {
        return StandardElement(StandardBoundedElement(AVXVector<limbs>(m1)), StandardBoundedElement(AVXVector<limbs>(m2)));
    }

    template<int64_t rns_bound>
    static INLINE OffsetElement<rns_bound> make_offset_element(
        const std::array<uint64_t, limbs>& m1,
        const std::array<uint64_t, limbs>& m2)
    {
        return OffsetElement<rns_bound>(FullyReducedElement(AVXVector<limbs>(m1)), FullyReducedElement(AVXVector<limbs>(m2)));
    }

    static INLINE AVXVector<limbs> scalar_fp_to_digits(const blst::Fp& value)
    {
        static_assert(limbs == 8, "The baked AVX field bridge currently supports only the 8-limb AVX layouts.");

        blst::Fp canonical = {};
        blst::Fp::from_mont(canonical, value);

        const uint64_t* limbs_64 = canonical.l;
        std::array<uint64_t, 8> digits = {};
        digits[0] = limbs_64[0] & ((1ULL << 52) - 1);
        digits[1] = ((limbs_64[0] >> 52) | (limbs_64[1] << 12)) & ((1ULL << 52) - 1);
        digits[2] = ((limbs_64[1] >> 40) | (limbs_64[2] << 24)) & ((1ULL << 52) - 1);
        digits[3] = ((limbs_64[2] >> 28) | (limbs_64[3] << 36)) & ((1ULL << 52) - 1);
        digits[4] = ((limbs_64[3] >> 16) | (limbs_64[4] << 48)) & ((1ULL << 52) - 1);
        digits[5] = (limbs_64[4] >> 4) & ((1ULL << 52) - 1);
        digits[6] = ((limbs_64[4] >> 56) | (limbs_64[5] << 8)) & ((1ULL << 52) - 1);
        digits[7] = (limbs_64[5] >> 44) & ((1ULL << 52) - 1);
        return AVXVector<limbs>(digits);
    }


    // One in RNS
    StandardElement one_;
    // Used to make the largest possible negative value positive.
    // WideOffsetElement<0, -RNS_MAX_NEG> max_neg_offset;
    OffsetElement<-RNS_MAX_NEG> max_neg_offset;
    // the modulus in RNS
    OffsetElement<1> modulus_in_rns;

    OffsetElement<6> six_modulus_in_rns;

    public:

    BoundedRing()
        :
        // 2 * RNS_MAX to allow for negative range
        // 9p * max redundance for elementwise reduction
          sys(make_IntRNS2SystemBaked<modulus_bits, limbs>(runtime_constants<element_bits>()))
         // There reducers differ in that they return without conditional correction (this should help when element_bits = 50, but has no advantage at 52)
        , montgomery_reduce_m1(sys.moduli1)
        , montgomery_reduce_m2(sys.moduli2)
        // To and from operate on 52 bit chunks independently of element bits, we let the system handle the digit calculation.
    {
        const auto& constants = runtime_constants<element_bits>();
        one_ = make_standard_element(constants.one_m1, constants.one_m2);
        max_neg_offset = make_offset_element<-RNS_MAX_NEG>(constants.max_neg_offset_m1, constants.max_neg_offset_m2);
        modulus_in_rns = make_offset_element<1>(constants.modulus_in_rns_m1, constants.modulus_in_rns_m2);
        six_modulus_in_rns = make_offset_element<6>(constants.six_modulus_in_rns_m1, constants.six_modulus_in_rns_m2);
    }

    static INLINE StandardElement zero() {
        return StandardElement{};
    }

    INLINE StandardElement one() const {
        return one_;
    }

    INLINE const auto& max_neg_offset_value() const {
        return max_neg_offset;
    }

    INLINE const auto& modulus_in_rns_value() const {
        return modulus_in_rns;
    }

    INLINE const auto& six_modulus_in_rns_value() const {
        return six_modulus_in_rns;
    }
    
    static INLINE auto wrap(const std::pair<AVXVector<limbs>, AVXVector<limbs>> &elements) {
        return StandardElement(
            StandardBoundedElement(elements.first), 
            StandardBoundedElement(elements.second)
        );
    };

    static INLINE auto wrap(const std::pair<std::array<uint64_t, limbs>, std::array<uint64_t, limbs>> &elements) {
        return wrap(std::make_pair(AVXVector<limbs>(elements.first), AVXVector<limbs>(elements.second)));
    }

    template<int bytes>
    INLINE auto from_hex(const char(&lit)[bytes]) const {
        blst::Fp canonical = {};
        const char* hex = lit;
        if (bytes >= 3 && lit[0] == '0' && (lit[1] == 'x' || lit[1] == 'X')) {
            hex = lit + 2;
        }
        const bool parsed =
            blst::limbs_from_hexascii(canonical.l, sizeof(canonical), hex);
        BLST_ASSERT(parsed);
        (void)parsed;
        blst::Fp mont = {};
        blst::Fp::to_mont(mont, canonical);
        return convert_to_rns_batch<1>(std::array<AVXVector<limbs>, 1>{scalar_fp_to_digits(mont)})[0];
    }

    template<int lower_bound, int upper_bound, class RNSBounds>
    INLINE auto recast(const ExpandedElement<BoundedElement<Bounds<lower_bound, upper_bound>, limbs, element_bits>, RNSBounds> &element) const {
        auto m1_recast = element.m1.template recast<lower_bound, upper_bound>();
        auto m2_recast = element.m2.template recast<lower_bound, upper_bound>();
        return ExpandedElement<decltype(m1_recast), RNSBounds>(m1_recast, m2_recast);
    }

    template<int multiple=1>
    INLINE auto elementwise_moduli() const {
        static_assert(multiple >= 1, "Multiple must be at least 1");
        // We only store log2(multiple) for efficiency.
        constexpr int required_multiple = ceil_log2(multiple);
        auto m1_moduli = montgomery_reduce_m1.template moduli<required_multiple>();
        auto m2_moduli = montgomery_reduce_m2.template moduli<required_multiple>();
        // This is 0 in the RNS system, so we can use the Bounds<0, 0> type.
        return ExpandedElement<decltype(m1_moduli), Bounds<0, 0>>(m1_moduli, m2_moduli);
    }

    template<int64_t lower, int64_t upper, class RNSBoundsType>
    INLINE auto wideify(const ExpandedElement<BoundedElement<Bounds<lower, upper>, limbs, element_bits>, RNSBoundsType> &element) const {
        auto recast_zero = BoundedElement<Bounds<0, 0>, limbs, element_bits>::zero().template recast<lower, upper>();
        using WideElementType = WideElement<Bounds<lower, upper>, limbs, element_bits>;
        auto m1 = WideElementType(recast_zero, element.m1);
        auto m2 = WideElementType(recast_zero, element.m2);
        return ExpandedElement<WideElementType, RNSBoundsType>(m1, m2);
    }
    
    INLINE auto negative_range_offset() const {
        return wideify(max_neg_offset);
    }

    // Prepare an element for a multiplication - elementwise reduce to make 52 bits.
    template<int64_t lower, int64_t upper, class RNSBoundsType, int to=MULT_OK_BOUND>
    INLINE auto prep(const ExpandedElement<BoundedElement<Bounds<lower, upper>, limbs, element_bits>, RNSBoundsType> &element) const {
        auto m1_reduced = montgomery_reduce_m1.template reduce_auto<to, upper, lower>(element.m1);
        auto m2_reduced = montgomery_reduce_m2.template reduce_auto<to, upper, lower>(element.m2);
        return ExpandedElement<decltype(m1_reduced), RNSBoundsType>(m1_reduced, m2_reduced);
    }

    template<class ElementType>
    INLINE auto prep_left(const ElementType &element) const {
        return prep(element);
    }

    template<int64_t lower, int64_t upper, class RNSBoundsType, int to=MULT_OK_BOUND>
    INLINE auto prep(const SmallElement<BoundedElement<Bounds<lower, upper>, limbs, element_bits>, RNSBoundsType> &element) const {
        auto m2_reduced = montgomery_reduce_m2.template reduce_auto<to, upper, lower>(element.element);
        return SmallElement<decltype(m2_reduced), RNSBoundsType>(m2_reduced);
    }

    template<int64_t lower, int64_t upper, class RNSBoundsType>
    INLINE auto prep_standard(const ExpandedElement<BoundedElement<Bounds<lower, upper>, limbs, element_bits>, RNSBoundsType> &element) const {
        return prep<lower, upper, RNSBoundsType, STANDARD_BOUND>(element);
    }

    template<int64_t lower, int64_t upper, class RNSBoundsType>
    INLINE auto prep_standard(const SmallElement<BoundedElement<Bounds<lower, upper>, limbs, element_bits>, RNSBoundsType> &element) const {
        return prep<lower, upper, RNSBoundsType, STANDARD_BOUND>(element);
    }
    
    template<int64_t upper_bound, int64_t lower_bound, class RNSBoundsType>
    INLINE auto negate(const ExpandedElement<BoundedElement<Bounds<lower_bound, upper_bound>, limbs, element_bits>, RNSBoundsType> &a) const {
        auto negation_offset = elementwise_moduli<upper_bound>();
        return negation_offset - a;
    }

    template<int64_t lower, int64_t upper, class RNSBoundsType>
    INLINE auto negate(const SmallElement<BoundedElement<Bounds<lower, upper>, limbs, element_bits>, RNSBoundsType> &a) const {
        auto negation_offset = elementwise_moduli<upper>();
        auto m2_negation_offset = negation_offset.m2;
        auto ns = SmallElement<decltype(m2_negation_offset), decltype(negation_offset.rns_bounds)>(m2_negation_offset);
        return ns - a;
    }

    // not as efficient as delaying negation but useful for making fp12 have a consistent input and output type
    template<int64_t lower, int64_t upper>
    INLINE StandardElement standard_negate(const ExpandedElement<BoundedElement<Bounds<lower, upper>, limbs, element_bits>, Bounds<0, 1>> &a) const {
        auto negated = negate(a) + modulus_in_rns;
        return prep<negated.m1.bounds.lower, negated.m1.bounds.upper, Bounds<0, 1>, STANDARD_BOUND>(negated);
    }

    template<int64_t lower, int64_t upper>
    INLINE auto standard_negate(const ExpandedElement<BoundedElement<Bounds<lower, upper>, limbs, element_bits>, Bounds<0, 6>> &a) const {
        auto negated = negate(a) + six_modulus_in_rns;
        return prep<negated.m1.bounds.lower, negated.m1.bounds.upper, Bounds<0, 6>, STANDARD_BOUND>(negated);
    }

    template<int64_t lower, int64_t upper>
    INLINE SmallElement<StandardBoundedElement, Bounds<0, 1>> standard_negate(const SmallElement<BoundedElement<Bounds<lower, upper>, limbs, element_bits>, Bounds<0, 1>> &a) const {
        auto modulus_in_rns_m2 = modulus_in_rns.m2;
        auto modulus_in_rns_s = SmallElement<decltype(modulus_in_rns_m2), decltype(modulus_in_rns.rns_bounds)>(modulus_in_rns_m2);
        auto negated = negate(a) + modulus_in_rns_s;
        return prep<negated.m2.bounds.lower, negated.m2.bounds.upper, Bounds<0, 1>, STANDARD_BOUND>(negated);
    }

    template<class ElementType, int64_t rns_lower, int64_t rns_upper>
    INLINE auto complete(const ExpandedElement<ElementType, Bounds<rns_lower, rns_upper>> &delayed_product) const {
        // automatically add negation offset if needed
        if constexpr(rns_lower < 0) {
            return negative_range_offset() + delayed_product;
        } else {
            // automatically complete delayed products if needed
            return delayed_product.complete();
        }
    }
    
    template<int64_t A, class ElementType, class RNSBoundsType>
    INLINE ReadyElement<A> ready(const ExpandedElement<ElementType, RNSBoundsType> &e) const {
        auto element = complete(e);
        using CompleteElementType = decltype(element);
        using LowType = decltype(element.m2.low);
        using HighType = decltype(element.m2.high);
        constexpr auto rns_lower_bound = CompleteElementType::rns_bounds.lower;
        constexpr auto rns_upper_bound = CompleteElementType::rns_bounds.upper;
        constexpr auto low_lower_bound = LowType::bounds.lower;
        constexpr auto high_lower_bound = HighType::bounds.lower;
        static_assert(rns_lower_bound >= 0, "RNS lower bound must be positive");
        static_assert(rns_upper_bound <= RNS_MAX_MULTIPLE, "RNS upper bound must be less than or equal to RNS_MAX_MULTIPLE");
        auto m2_recast_low = element.m2.low.template recast<0, A>();
        auto m2_recast_high = element.m2.high.template recast<0, A>();
        auto m2_recast = WideElement<Bounds<0, A>, limbs, element_bits>(m2_recast_high, m2_recast_low);
        static_assert(low_lower_bound >= 0);
        static_assert(high_lower_bound >= 0);
        auto m1_mont_ele_reduced = montgomery_reduce_m1.template reduce_auto<MULT_OK_BOUND>(element.m1);
        auto m1_recast = m1_mont_ele_reduced.template recast<0, MULT_OK_BOUND>();
        return ReadyToReduce<WideElement<Bounds<0, A>, limbs, element_bits>, MultOkBoundedElement>(m2_recast, m1_recast);
    }

    // Expand, maintains RNS bound and allows for multiplications
    template<int batch, int64_t RNS_Bound>
    INLINE std::array<ExpandedElement<StandardBoundedElement, Bounds<0, RNS_Bound>>, batch>
    batch_expand(const std::array<SmallElement<StandardBoundedElement, Bounds<0, RNS_Bound>>, batch> &c_m2) const {
        auto wide_zero = WideZero<limbs, element_bits>();
        using ZeroAccElement = ReadyToReduce<decltype(wide_zero), StandardBoundedElement>;
        std::array<ZeroAccElement, batch> zero_acc_elements;
        for (int i = 0; i < batch; i++) {
            zero_acc_elements[i] = ZeroAccElement(wide_zero, c_m2[i].element);
        }
        auto c_m1_acc = batch_change_base<batch>(zero_acc_elements, sys.intrns2.r2);
        std::array<ExpandedElement<StandardBoundedElement, Bounds<0, RNS_Bound>>, batch> result;
        for (int i = 0; i < batch; i++) {
            result[i].m1 = montgomery_reduce_m1.template reduce_auto<STANDARD_BOUND>(c_m1_acc[i]);
            result[i].m2 = c_m2[i].element;
        }
        return result;
    }

    template<int batch>
    INLINE std::array<StandardBoundedElement, batch>
    batch_expand(const std::array<StandardBoundedElement, batch> &c_m2) const {
        auto wide_zero = WideZero<limbs, element_bits>();
        using ZeroAccElement = ReadyToReduce<decltype(wide_zero), StandardBoundedElement>;
        std::array<ZeroAccElement, batch> zero_acc_elements;
        for (int i = 0; i < batch; i++) {
            zero_acc_elements[i] = ZeroAccElement(wide_zero, c_m2[i]);
        }
        auto c_m1_acc = batch_change_base<batch>(zero_acc_elements, sys.intrns2.r2);
        std::array<StandardBoundedElement, batch> result;
        for (int i = 0; i < batch; i++) {
            result[i] = montgomery_reduce_m1.template reduce_auto<STANDARD_BOUND>(c_m1_acc[i]);
        }
        return result;
    }

    template<class ElementType, class RNSBoundsType>
    INLINE auto small_m2(const ExpandedElement<ElementType, RNSBoundsType>& element) const {
        return SmallElement<ElementType, RNSBoundsType>(element.m2);
    }


    template<typename... Args, std::size_t... I>
    INLINE auto make_expand_array_impl(const std::tuple<Args...>& args, std::index_sequence<I...>) const {
        return std::array<StandardBoundedElement, sizeof...(Args)>{std::get<I>(args).element...};
    }

    template<class E, class Expanded>
    INLINE auto reattach_expand_bounds(const E& e, const Expanded& expanded) const {
        return ExpandedElement<StandardBoundedElement, std::remove_cv_t<decltype(e.rns_bounds)>>(expanded, e);
    }

    template<typename... Args, std::size_t... I>
    INLINE auto reattach_expand_tuple_impl(
        const std::tuple<Args...>& args,
        const std::array<StandardBoundedElement, sizeof...(Args)>& expanded,
        std::index_sequence<I...>) const {
        return std::tuple<decltype(reattach_expand_bounds(std::get<I>(args), expanded[I]))...>{
            reattach_expand_bounds(std::get<I>(args), expanded[I])...};
    }

    template<typename... Args>
    INLINE auto batch_expand(const Args&... args) const {
        constexpr int batch = static_cast<int>(sizeof...(Args));
        static_assert(batch > 0, "batch_expand requires at least one argument");
        const auto tuple = std::tie(args...);
        const auto data = make_expand_array_impl(tuple, std::make_index_sequence<sizeof...(Args)>{});
        const auto expanded = batch_expand<batch>(data);
        return reattach_expand_tuple_impl(tuple, expanded, std::make_index_sequence<sizeof...(Args)>{});
    }

    // Reduce, sets RNS bound to 1
    template<int batch, int64_t A=MAX_ADD>
    INLINE auto batch_reduce(const std::array<ReadyElement<A>, batch> &c_in) const {
        auto c_m2_acc = batch_change_base<batch>(c_in, sys.intrns2.r1);
        // We split reduce wide and reduce auto to enable scalar operations on unexpanded form.
        using MontHalfType = decltype(montgomery_reduce_m2.reduce_wide(c_m2_acc[0]));
        // Apply RNS Bounds manually since now RNS reduced
        using OutType = SmallElement<MontHalfType, Bounds<0, 1>>;
        std::array<OutType, batch> result;
        for (int i = 0; i < batch; i++) {
            result[i] = OutType(montgomery_reduce_m2.reduce_wide(c_m2_acc[i]));
        }
        return result;
    }

    template<class ElementType, class RNSBoundsType>
    INLINE auto prep_expand(const SmallElement<ElementType, RNSBoundsType> &e) const {
        return SmallElement<StandardBoundedElement, RNSBoundsType>(montgomery_reduce_m2.template reduce_auto<STANDARD_BOUND>(e.element));
    }

    // Helper to do reduction followed by expansion.
    template<int batch>
    INLINE std::array<StandardElement, batch> batch_reduce_expand(const std::array<ReadyElement<MAX_ADD>, batch> &elements) const {
        auto cm2_acc = batch_reduce<batch>(elements);
        std::array<SmallStandardElement, batch> c_m2;
        for (int i = 0; i < batch; i++) {
            c_m2[i] = SmallElement<StandardBoundedElement, Bounds<0, 1>>(montgomery_reduce_m2.template reduce_auto<STANDARD_BOUND>(cm2_acc[i].element));
        }
        return batch_expand<batch>(c_m2);
    }



    template<int batch>
    INLINE std::array<StandardElement, batch> convert_to_rns_batch(const std::array<AVXVector<limbs>, batch> &c_m2) const {
        std::array<AVXVector<limbs>, batch> c_m2_hi{};
        std::array<AVXVector<limbs>, batch> c_m2_lo{};
        sys.convert_to.conversion_matrix.template rns_reduce_raw_batch<batch>(&c_m2, c_m2_hi, c_m2_lo);
        std::array<SmallStandardElement, batch> result;
        using WideBoundedElement = BoundedElement<Bounds<0, MAX_ADD>, limbs, element_bits>;
        for (int i = 0; i < batch; i++) {
            result[i] = SmallElement<StandardBoundedElement, Bounds<0, 1>>(montgomery_reduce_m2.template reduce_auto<STANDARD_BOUND>(WideElement<Bounds<0, MAX_ADD>, limbs, element_bits>(WideBoundedElement(c_m2_hi[i]), WideBoundedElement(c_m2_lo[i]))));
        }
        return batch_expand<batch>(result);
    }
    

    // Helper to convert variadic arguments to array using index_sequence
    template<int64_t A, typename... Args, std::size_t... I>
    INLINE auto make_ready_array_impl(const std::tuple<Args...> &args, std::index_sequence<I...>) const {
        return std::array<ReadyElement<A>, sizeof...(Args)>{ready<A>(std::get<I>(args))...};
    }

    template<int64_t A, typename... Args, std::size_t... I>
    INLINE auto make_prep_array_impl(const std::tuple<Args...> &args, std::index_sequence<I...>) const {
        return std::array<SmallStandardElement, sizeof...(Args)>{prep(std::get<I>(args))...};
    }

     // Variadic template version that captures all the different arities
     template<int64_t A, typename... Args>
     INLINE auto batch_ready(const Args&... args) const {
         constexpr int batch = static_cast<int>(sizeof...(Args));
         static_assert(batch > 0, "batch_reduce_expand requires at least one argument");
         auto tuple = std::tie(args...);
         auto ready_array = make_ready_array_impl<A>(tuple, std::make_index_sequence<sizeof...(Args)>{});
         return ready_array;
    }

    // Variadic template version that captures all the different arities
    template<typename... Args>
    INLINE auto batch_reduce_expand(const Args&... args) const {
        constexpr int batch = static_cast<int>(sizeof...(Args));
        auto reduced = batch_ready<MAX_ADD>(args...);
        return batch_reduce_expand<batch>(reduced);
    }

    template<typename... Args>
    INLINE auto batch_reduce(const Args&... args) const {
        constexpr int batch = static_cast<int>(sizeof...(Args));
        auto reduced = batch_ready<MULT_OK_BOUND*MULT_OK_BOUND>(args...);
        return batch_reduce<batch>(reduced);
    }

    INLINE StandardElement modmul(const StandardElement &a, const StandardElement &b) const {
        auto product = (a * b);
        auto ready_product = ready<MAX_ADD>(product);
        std::array<ReadyElement<MAX_ADD>, 1> ready_elements = {ready_product};
        return batch_reduce_expand<1>(ready_elements)[0];
    }

    template<int batch>
    INLINE std::array<StandardElement, batch> batch_modmul(const std::array<StandardElement, batch> &a, const std::array<StandardElement, batch> &b) const {
        std::array<ReadyElement<MAX_ADD>, batch> ready_elements;
        for (int i = 0; i < batch; i++) {
            ready_elements[i] = ready<MAX_ADD>(a[i] * b[i]);
        }
        return batch_reduce_expand<batch>(ready_elements);
    }

    // This probably should be on the curve class
    template<class Element>
    INLINE auto mul_3b(const Element &x) const {
        return x * std::integral_constant<int, 12>();
    }


    // Find a number matching the bounds.

};

}  // namespace blst::avx
