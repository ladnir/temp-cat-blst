#include "runtime_constants_validation_support.h"

#include "blst/avx/cpu/vector/multiplication.hpp"
#include "blst/avx/fields/fp12.h"
#include "blst/avx/rns/bounded_ring.h"
#include "blst/avx/runtime_constants.h"

#include <array>
#include <stdexcept>
#include <string>

namespace {

using blst::avx::BigInt;
using blst::avx::FP12Ring;
using blst::avx::RuntimeConstantsView;
using blst::avx::make_IntRNS2SystemMontgomery;
using blst::avx::runtime_constants;
using blst::avx::uint128_t;

constexpr char kBls12_381BaseModulus[] =
    "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";

template<class T>
void expect_equal(const T& actual, const T& expected, const char* label)
{
    if (!(actual == expected)) {
        throw std::runtime_error(std::string("runtime constant mismatch: ") + label);
    }
}

template<class Element>
std::array<uint64_t, 8> store_bounded(const Element& element)
{
    std::array<uint64_t, 8> out = {};
    element.data.store(out.data());
    return out;
}

template<class Element>
std::pair<std::array<uint64_t, 8>, std::array<uint64_t, 8>> store_expanded(const Element& element)
{
    return {store_bounded(element.m1), store_bounded(element.m2)};
}

template<class Ring>
struct StoredWideOffset {
    std::array<uint64_t, 8> m1_high;
    std::array<uint64_t, 8> m1_low;
    std::array<uint64_t, 8> m2_high;
    std::array<uint64_t, 8> m2_low;
};

template<class Ring, class Element>
StoredWideOffset<Ring> store_wide_offset(const Element& element)
{
    return {
        store_bounded(element.m1.high),
        store_bounded(element.m1.low),
        store_bounded(element.m2.high),
        store_bounded(element.m2.low),
    };
}

template<class Ring>
StoredWideOffset<Ring> expected_fp12_offset()
{
    const BigInt target(kBls12_381BaseModulus, 10);
    const auto sys = make_IntRNS2SystemMontgomery<381, 8>(target, Ring::StandardBoundedElement::ELEMENT_BITS);

    constexpr int offset_rns_bound =
        10 * Ring::STANDARD_BOUND * Ring::STANDARD_BOUND * 36;
    constexpr int offset_ele_bound =
        10 * Ring::STANDARD_BOUND * ((Ring::STANDARD_BOUND > Ring::MULT_OK_BOUND) ? Ring::STANDARD_BOUND : Ring::MULT_OK_BOUND);

    const BigInt desired_value = BigInt(9) * BigInt(offset_rns_bound) * target * target;
    const auto desired_value_rns = sys.to_mont_exact(desired_value, true);

    StoredWideOffset<Ring> out = {};
    const uint64_t low_part_min = (1ull << 52) * offset_ele_bound;
    const uint64_t low_mask = (1ull << 52) - 1;

    for (int j = 0; j < 8; ++j) {
        const uint128_t val_m1 = uint128_t(offset_ele_bound) * uint128_t(sys.moduli1[j]) * uint128_t(sys.moduli1[j]);
        const uint128_t val_m2 = uint128_t(offset_ele_bound) * uint128_t(sys.moduli2[j]) * uint128_t(sys.moduli2[j]);

        out.m1_low[j] = static_cast<uint64_t>(val_m1 & low_mask);
        out.m2_low[j] = static_cast<uint64_t>(val_m2 & low_mask);
        out.m1_high[j] = static_cast<uint64_t>(val_m1 >> 52);
        out.m2_high[j] = static_cast<uint64_t>(val_m2 >> 52);

        out.m1_low[j] += desired_value_rns.first[j];
        out.m2_low[j] += desired_value_rns.second[j];
        out.m1_low[j] += sys.moduli1[j] * offset_ele_bound;
        out.m2_low[j] += sys.moduli2[j] * offset_ele_bound;

        while (out.m1_low[j] < low_part_min) {
            out.m1_low[j] += sys.moduli1[j];
        }
        while (out.m2_low[j] < low_part_min) {
            out.m2_low[j] += sys.moduli2[j];
        }
    }

    return out;
}

template<int element_bits>
void validate_runtime_constants_impl()
{
    using Ring = blst::avx::BoundedRing<381, 8, element_bits>;

    const BigInt target(kBls12_381BaseModulus, 10);
    const auto sys = make_IntRNS2SystemMontgomery<381, 8>(target, element_bits);
    const RuntimeConstantsView& constants = runtime_constants<element_bits>();

    expect_equal(constants.moduli1, sys.moduli1, "moduli1");
    expect_equal(constants.moduli2, sys.moduli2, "moduli2");

    std::array<std::array<uint64_t, 8>, 8> matrix = {};
    std::array<uint64_t, 8> correction = {};
    std::array<uint64_t, 8> sqe = {};
    std::array<uint64_t, 8> reducer_values = {};

    sys.intrns2.m1.get_moduli(reducer_values);
    expect_equal(constants.reducer1_moduli, reducer_values, "reducer1_moduli");
    sys.intrns2.m1.get_mont_factor(reducer_values);
    expect_equal(constants.reducer1_mont_factor, reducer_values, "reducer1_mont_factor");
    sys.intrns2.m1.get_t_squared(reducer_values);
    expect_equal(constants.reducer1_t_squared, reducer_values, "reducer1_t_squared");

    sys.intrns2.m2.get_moduli(reducer_values);
    expect_equal(constants.reducer2_moduli, reducer_values, "reducer2_moduli");
    sys.intrns2.m2.get_mont_factor(reducer_values);
    expect_equal(constants.reducer2_mont_factor, reducer_values, "reducer2_mont_factor");
    sys.intrns2.m2.get_t_squared(reducer_values);
    expect_equal(constants.reducer2_t_squared, reducer_values, "reducer2_t_squared");

    sys.intrns2.r1.get_rns_values(matrix);
    expect_equal(constants.r1_rns_mat, matrix, "r1_rns_mat");
    sys.intrns2.r1.get_correction(correction);
    expect_equal(constants.r1_correction, correction, "r1_correction");
    sys.intrns2.r1.get_shifted_quotient_estimations(sqe);
    expect_equal(constants.r1_shifted_quotient_estimations, sqe, "r1_shifted_quotient_estimations");

    sys.intrns2.r2.get_rns_values(matrix);
    expect_equal(constants.r2_rns_mat, matrix, "r2_rns_mat");
    sys.intrns2.r2.get_correction(correction);
    expect_equal(constants.r2_correction, correction, "r2_correction");
    sys.intrns2.r2.get_shifted_quotient_estimations(sqe);
    expect_equal(constants.r2_shifted_quotient_estimations, sqe, "r2_shifted_quotient_estimations");

    sys.convert_to.get_conversion_matrix().get_rns_values(matrix);
    expect_equal(constants.convert_to_rns_mat, matrix, "convert_to_rns_mat");
    sys.convert_to.get_conversion_matrix().get_correction(correction);
    expect_equal(constants.convert_to_correction, correction, "convert_to_correction");
    sys.convert_to.get_conversion_matrix().get_shifted_quotient_estimations(sqe);
    expect_equal(constants.convert_to_shifted_quotient_estimations, sqe, "convert_to_shifted_quotient_estimations");

    sys.convert_from.get_conversion_matrix().get_rns_values(matrix);
    expect_equal(constants.convert_from_rns_mat, matrix, "convert_from_rns_mat");
    sys.convert_from.get_conversion_matrix().get_correction(correction);
    expect_equal(constants.convert_from_correction, correction, "convert_from_correction");
    sys.convert_from.get_conversion_matrix().get_shifted_quotient_estimations(sqe);
    expect_equal(constants.convert_from_shifted_quotient_estimations, sqe, "convert_from_shifted_quotient_estimations");

    const auto one_expected = sys.to_mont_exact(BigInt(1));
    expect_equal(constants.one_m1, one_expected.first, "one_m1");
    expect_equal(constants.one_m2, one_expected.second, "one_m2");

    const auto modulus_expected = sys.to_mont_exact(BigInt(3) * target);
    expect_equal(constants.modulus_in_rns_m1, modulus_expected.first, "modulus_in_rns_m1");
    expect_equal(constants.modulus_in_rns_m2, modulus_expected.second, "modulus_in_rns_m2");

    const auto six_modulus_expected = sys.to_mont_exact(BigInt(18) * target);
    expect_equal(constants.six_modulus_in_rns_m1, six_modulus_expected.first, "six_modulus_in_rns_m1");
    expect_equal(constants.six_modulus_in_rns_m2, six_modulus_expected.second, "six_modulus_in_rns_m2");

    const auto max_neg_offset_expected = sys.to_mont_exact(BigInt(1932 * 9) * target * target, true);
    expect_equal(constants.max_neg_offset_m1, max_neg_offset_expected.first, "max_neg_offset_m1");
    expect_equal(constants.max_neg_offset_m2, max_neg_offset_expected.second, "max_neg_offset_m2");

    const auto fp12_offset_expected = expected_fp12_offset<Ring>();
    expect_equal(constants.fp12_offset_m1_high, fp12_offset_expected.m1_high, "fp12_offset_m1_high");
    expect_equal(constants.fp12_offset_m1_low, fp12_offset_expected.m1_low, "fp12_offset_m1_low");
    expect_equal(constants.fp12_offset_m2_high, fp12_offset_expected.m2_high, "fp12_offset_m2_high");
    expect_equal(constants.fp12_offset_m2_low, fp12_offset_expected.m2_low, "fp12_offset_m2_low");

    const Ring ring;
    const auto one_actual = store_expanded(ring.one());
    expect_equal(constants.one_m1, one_actual.first, "ring.one().m1");
    expect_equal(constants.one_m2, one_actual.second, "ring.one().m2");

    const auto max_neg_actual = store_expanded(ring.max_neg_offset_value());
    expect_equal(constants.max_neg_offset_m1, max_neg_actual.first, "ring.max_neg_offset_value().m1");
    expect_equal(constants.max_neg_offset_m2, max_neg_actual.second, "ring.max_neg_offset_value().m2");

    const auto modulus_actual = store_expanded(ring.modulus_in_rns_value());
    expect_equal(constants.modulus_in_rns_m1, modulus_actual.first, "ring.modulus_in_rns_value().m1");
    expect_equal(constants.modulus_in_rns_m2, modulus_actual.second, "ring.modulus_in_rns_value().m2");

    const auto six_modulus_actual = store_expanded(ring.six_modulus_in_rns_value());
    expect_equal(constants.six_modulus_in_rns_m1, six_modulus_actual.first, "ring.six_modulus_in_rns_value().m1");
    expect_equal(constants.six_modulus_in_rns_m2, six_modulus_actual.second, "ring.six_modulus_in_rns_value().m2");

    const FP12Ring<Ring> fp12_ring(ring);
    const auto fp12_actual = store_wide_offset<Ring>(fp12_ring.offset_value());
    expect_equal(constants.fp12_offset_m1_high, fp12_actual.m1_high, "fp12_ring.offset_value().m1_high");
    expect_equal(constants.fp12_offset_m1_low, fp12_actual.m1_low, "fp12_ring.offset_value().m1_low");
    expect_equal(constants.fp12_offset_m2_high, fp12_actual.m2_high, "fp12_ring.offset_value().m2_high");
    expect_equal(constants.fp12_offset_m2_low, fp12_actual.m2_low, "fp12_ring.offset_value().m2_low");
}

}  // namespace

void validate_runtime_constants_52()
{
    validate_runtime_constants_impl<52>();
}

void validate_runtime_constants_50()
{
    validate_runtime_constants_impl<50>();
}
