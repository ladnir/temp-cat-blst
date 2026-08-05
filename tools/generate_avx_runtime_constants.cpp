#include "blst/avx/fields/fp12.h"
#include "blst/avx/fields/conversion_inversion.h"
#include "blst/avx/rns/bounded_ring.h"
#include "blst/scalar/curve/GT.h"

#include <array>
#include <iostream>
#include <string_view>

namespace {

template <typename T, std::size_t N>
void print_array(const std::array<T, N>& values, std::string_view name)
{
    std::cout << "inline constexpr std::array<uint64_t, " << N << "> " << name << " = {";
    for (std::size_t i = 0; i < N; ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << static_cast<unsigned long long>(values[i]) << "ULL";
    }
    std::cout << "}};\n";
}

template <typename T, std::size_t Rows, std::size_t Cols>
void print_matrix(const std::array<std::array<T, Cols>, Rows>& values, std::string_view name)
{
    std::cout << "inline constexpr std::array<std::array<uint64_t, " << Cols << ">, " << Rows << "> " << name << " = {{";
    for (std::size_t row = 0; row < Rows; ++row) {
        if (row != 0) {
            std::cout << ", ";
        }
        std::cout << "{";
        for (std::size_t col = 0; col < Cols; ++col) {
            if (col != 0) {
                std::cout << ", ";
            }
            std::cout << static_cast<unsigned long long>(values[row][col]) << "ULL";
        }
        std::cout << "}";
    }
    std::cout << "};\n";
}

template <class Element>
void print_standard_element(const Element& value, std::string_view name)
{
    print_array(value.m1.to_unsigned_array(), std::string(name) + "_m1");
    print_array(value.m2.to_unsigned_array(), std::string(name) + "_m2");
}

template <class Element>
void print_fp12_constant(const Element& value, std::string_view name)
{
    std::cout << "inline constexpr Fp12Constant " << name << " = {{\n";
    for (const auto& coefficient : value) {
        const auto m1 = coefficient.m1.to_unsigned_array();
        const auto m2 = coefficient.m2.to_unsigned_array();
        std::cout << "    {{{";
        for (std::size_t i = 0; i < m1.size(); ++i) {
            std::cout << (i == 0 ? "" : ", ")
                      << static_cast<unsigned long long>(m1[i]) << "ULL";
        }
        std::cout << "}, {";
        for (std::size_t i = 0; i < m2.size(); ++i) {
            std::cout << (i == 0 ? "" : ", ")
                      << static_cast<unsigned long long>(m2[i]) << "ULL";
        }
        std::cout << "}}},\n";
    }
    std::cout << "}};\n";
}

template <class Ring, class WideElement>
void print_wide_element(const WideElement& value, std::string_view name)
{
    print_array(value.m1.high.to_unsigned_array(), std::string(name) + "_m1_high");
    print_array(value.m1.low.to_unsigned_array(), std::string(name) + "_m1_low");
    print_array(value.m2.high.to_unsigned_array(), std::string(name) + "_m2_high");
    print_array(value.m2.low.to_unsigned_array(), std::string(name) + "_m2_low");
}

template <class Ring>
void print_matrix_data(
    const auto& matrix,
    std::string_view prefix)
{
    constexpr std::size_t kLimbs = Ring::StandardBoundedElement::LIMBS;
    std::array<std::array<uint64_t, kLimbs>, kLimbs> rns_mat = {};
    std::array<uint64_t, kLimbs> correction = {};
    std::array<uint64_t, kLimbs> sqe = {};
    matrix.get_rns_mat(rns_mat);
    matrix.get_correction(correction);
    matrix.get_shifted_quotient_estimations(sqe);
    print_matrix(rns_mat, std::string(prefix) + "_rns_mat");
    print_array(correction, std::string(prefix) + "_correction");
    print_array(sqe, std::string(prefix) + "_shifted_quotient_estimations");
}

template <class Ring>
void print_config(std::string_view ns_name)
{
    constexpr std::size_t limbs = Ring::StandardBoundedElement::LIMBS;
    const Ring ring(blst::avx::BigInt(
        "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787",
        10));
    blst::avx::FP12Ring<Ring> fp12(ring);

    std::cout << "namespace " << ns_name << " {\n";

    print_array(ring.sys.moduli1, "moduli1");
    print_array(ring.sys.moduli2, "moduli2");

    std::array<uint64_t, limbs> reducer1_moduli = {};
    std::array<uint64_t, limbs> reducer1_mont = {};
    std::array<uint64_t, limbs> reducer1_t_squared = {};
    ring.montgomery_reduce_m1.get_moduli(reducer1_moduli);
    ring.montgomery_reduce_m1.get_mont_factor(reducer1_mont);
    ring.montgomery_reduce_m1.get_t_squared(reducer1_t_squared);
    print_array(reducer1_moduli, "reducer1_moduli");
    print_array(reducer1_mont, "reducer1_mont_factor");
    print_array(reducer1_t_squared, "reducer1_t_squared");

    std::array<uint64_t, limbs> reducer2_moduli = {};
    std::array<uint64_t, limbs> reducer2_mont = {};
    std::array<uint64_t, limbs> reducer2_t_squared = {};
    ring.montgomery_reduce_m2.get_moduli(reducer2_moduli);
    ring.montgomery_reduce_m2.get_mont_factor(reducer2_mont);
    ring.montgomery_reduce_m2.get_t_squared(reducer2_t_squared);
    print_array(reducer2_moduli, "reducer2_moduli");
    print_array(reducer2_mont, "reducer2_mont_factor");
    print_array(reducer2_t_squared, "reducer2_t_squared");

    print_matrix_data<Ring>(ring.sys.intrns2.r1, "r1");
    print_matrix_data<Ring>(ring.sys.intrns2.r2, "r2");
    print_matrix_data<Ring>(ring.sys.convert_to.conversion_matrix, "convert_to");
    print_matrix_data<Ring>(ring.sys.convert_from.conversion_matrix, "convert_from");

    print_standard_element(ring.one(), "one");
    print_standard_element(ring.modulus_in_rns_value(), "modulus_in_rns");
    print_standard_element(ring.six_modulus_in_rns_value(), "six_modulus_in_rns");
    print_standard_element(ring.max_neg_offset_value(), "max_neg_offset");

    print_wide_element<Ring>(fp12.offset_value(), "fp12_offset");

    if constexpr (Ring::StandardBoundedElement::ELEMENT_BITS == 52) {
        print_fp12_constant(
            blst::avx::avx_blst_direct_bridge::scalar_fp12_to_ring(
                ring, blst::scalar::GT::generator().data),
            "gt_generator");
    }

    std::cout << "}  // namespace " << ns_name << "\n\n";
}

}  // namespace

int main()
{
    std::cout << "#pragma once\n\n";
    std::cout << "#include <array>\n";
    std::cout << "#include <cstdint>\n\n";
    std::cout << "namespace blst::avx::generated {\n\n";
    std::cout << "using Fp12Constant = std::array<std::array<std::array<uint64_t, 8>, 2>, 12>;\n\n";
    print_config<blst::avx::BoundedRing<381, 8, 52>>("fp52");
    print_config<blst::avx::BoundedRing<381, 8, 50>>("fp50");
    std::cout << "}  // namespace blst::avx::generated\n";
    return 0;
}
