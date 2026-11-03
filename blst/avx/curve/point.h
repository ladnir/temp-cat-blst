/*
 * Copyright Supranational LLC and CatBLST contributors
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 *
 * Imported with the VROOM AVX backend and adapted to CatBLST's homogeneous
 * point representation, selected-backend API, and point-domain contracts.
 */
#pragma once

#include "ec.h"
#include "matrixnok_g1.h"
#include "matrixnok_g2.h"
#include "scalar_mult.h"
#include "../generated_matrixnok_curve_constants.hpp"
#include "blst/avx/fields/conversion_inversion.h"
#include "blst/scalar/curve/point.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>

namespace blst::avx {

namespace avx_curve_detail {

using CurveRing = matrixnok::FpRing;
using CurveFp2Ring = matrixnok::G2Fp2Ring;
using CurveG1 = matrixnok::G1Curve;
using CurveG2 = matrixnok::G2Curve;
using CurveG1Projective = matrixnok::G1Projective;
using CurveG1Affine = matrixnok::G1Affine;
using CurveG2Projective = matrixnok::G2Projective;
using CurveG2Affine = matrixnok::G2Affine;
using CurveG1GLV = matrixnok::G1GLV;
using CurveG2GLS = matrixnok::G2GLS;

constexpr AVXVector<CurveRing::limbs> literal_vector(
    const std::array<uint64_t, CurveRing::limbs>& lanes)
{
    return AVXVector<CurveRing::limbs>(_vv512_const_u64(
        lanes[0], lanes[1], lanes[2], lanes[3],
        lanes[4], lanes[5], lanes[6], lanes[7]));
}

constexpr CurveRing::StandardElement literal_standard_element(
    const std::array<uint64_t, CurveRing::limbs>& m1,
    const std::array<uint64_t, CurveRing::limbs>& m2)
{
    return CurveRing::StandardElement(
        CurveRing::StandardBoundedElement(literal_vector(m1)),
        CurveRing::StandardBoundedElement(literal_vector(m2)));
}

constexpr CurveFp2Ring::StandardElement literal_fp2_element(
    const std::array<uint64_t, CurveRing::limbs>& x_m1,
    const std::array<uint64_t, CurveRing::limbs>& x_m2,
    const std::array<uint64_t, CurveRing::limbs>& y_m1,
    const std::array<uint64_t, CurveRing::limbs>& y_m2)
{
    return CurveFp2Ring::StandardElement(
        literal_standard_element(x_m1, x_m2),
        literal_standard_element(y_m1, y_m2));
}

inline constexpr CurveRing kCurveRing{};
inline constexpr CurveFp2Ring kCurveFp2Ring{kCurveRing};
inline constexpr CurveG1 kCurveG1{};
inline constexpr CurveG2 kCurveG2{};
inline constexpr CurveG1GLV kCurveG1Glv(
    kCurveG1,
    literal_standard_element(
        matrixnok::generated::fp12_frob1_4_m1,
        matrixnok::generated::fp12_frob1_4_m2));
inline constexpr CurveG2GLS kCurveG2Gls(
    literal_standard_element(
        matrixnok::generated::fp12_frob1_3_m1,
        matrixnok::generated::fp12_frob1_3_m2),
    literal_standard_element(
        matrixnok::generated::fp12_frob1_8_m1,
        matrixnok::generated::fp12_frob1_8_m2),
    literal_standard_element(
        matrixnok::generated::fp12_frob1_7_m1,
        matrixnok::generated::fp12_frob1_7_m2));
inline constexpr CurveFp2Ring::StandardElement kPsiFrobeniusX(
    literal_standard_element(
        matrixnok::generated::aprime_e2_x_m1,
        matrixnok::generated::aprime_e2_x_m2),
    kCurveG2Gls.psi_x_i);
inline constexpr CurveFp2Ring::StandardElement kPsiFrobeniusY(
    kCurveG2Gls.psi_y_r,
    kCurveG2Gls.psi_y_i);
inline constexpr CurveRing::StandardElement kAprimeE1 =
    literal_standard_element(
        matrixnok::generated::aprime_e1_m1,
        matrixnok::generated::aprime_e1_m2);
inline constexpr CurveFp2Ring::StandardElement kAprimeE2 =
    literal_fp2_element(
        matrixnok::generated::aprime_e2_x_m1,
        matrixnok::generated::aprime_e2_x_m2,
        matrixnok::generated::aprime_e2_y_m1,
        matrixnok::generated::aprime_e2_y_m2);

inline const CurveRing& curve_ring()
{
    return kCurveRing;
}

inline const CurveFp2Ring& curve_fp2_ring()
{
    return kCurveFp2Ring;
}

inline const CurveG1& curve_g1()
{
    return kCurveG1;
}

inline const CurveG2& curve_g2()
{
    return kCurveG2;
}

inline const CurveG1GLV& curve_g1_glv()
{
    return kCurveG1Glv;
}

inline const CurveG2GLS& curve_g2_gls()
{
    return kCurveG2Gls;
}

inline bool is_zero(const CurveRing::StandardElement& value)
{
    return value.is_zero() != 0;
}

inline bool is_zero(const CurveFp2Ring::StandardElement& value)
{
    return value.is_zero() != 0;
}

inline blst::Fp ring_fp_to_scalar(const CurveRing::StandardElement& value)
{
    return avx_blst_direct_bridge::ring_to_scalar_fp(curve_ring(), value);
}

inline CurveRing::StandardElement scalar_fp_to_ring(const blst::Fp& value)
{
    return avx_blst_direct_bridge::scalar_fp_to_ring(curve_ring(), value);
}

inline blst::Fp2 ring_fp2_to_scalar(const CurveFp2Ring::StandardElement& value)
{
    blst::Fp2 out = {};
    out[0] = ring_fp_to_scalar(value.x());
    out[1] = ring_fp_to_scalar(value.y());
    return out;
}

inline CurveFp2Ring::StandardElement scalar_fp2_to_ring(const blst::Fp2& value)
{
    return CurveFp2Ring::StandardElement(
        scalar_fp_to_ring(value[0]),
        scalar_fp_to_ring(value[1]));
}

inline blst::scalar::G1AffinePoint curve_g1_affine_to_scalar(const CurveG1Affine& point)
{
    blst::scalar::G1AffinePoint out = {};
    out.x = ring_fp_to_scalar(point.x);
    out.y = ring_fp_to_scalar(point.y);
    return out;
}

inline blst::scalar::G1Point curve_g1_to_scalar(const CurveG1Projective& point)
{
    // The AVX curve uses homogeneous coordinates (X/Z, Y/Z), while the
    // scalar curve uses Jacobian coordinates (X/Z^2, Y/Z^3). Convert without
    // an inversion: (Xh*Zh, Yh*Zh^2, Zh) represents the same affine point.
    const auto z_squared = curve_ring().modmul(point.z, point.z);
    const auto jacobian_x = curve_ring().modmul(point.x, point.z);
    const auto jacobian_y = curve_ring().modmul(point.y, z_squared);
    blst::scalar::G1Point out = {};
    out.x = ring_fp_to_scalar(jacobian_x);
    out.y = ring_fp_to_scalar(jacobian_y);
    out.z = ring_fp_to_scalar(point.z);
    return out;
}

// Public observations must compare field values, not raw redundant RNS lanes.
// Keep is_zero() above for internal owners whose representation invariant is
// exact zero; this boundary uses the proved canonical export.
inline bool_t is_zero_mod_field(const CurveRing::StandardElement& value)
{
    const blst::Fp scalar = ring_fp_to_scalar(value);
    return blst::Fp::equal(scalar, blst::Fp::zero());
}

inline bool_t is_zero_mod_field(const CurveFp2Ring::StandardElement& value)
{
    const blst::Fp2 scalar = ring_fp2_to_scalar(value);
    return blst::Fp::equal(scalar[0], blst::Fp::zero()) &
        blst::Fp::equal(scalar[1], blst::Fp::zero());
}

inline bool_t equal_mod_field(
    const CurveRing::StandardElement& left,
    const CurveRing::StandardElement& right)
{
    return blst::Fp::equal(
        ring_fp_to_scalar(left), ring_fp_to_scalar(right));
}

inline bool_t equal_mod_field(
    const CurveFp2Ring::StandardElement& left,
    const CurveFp2Ring::StandardElement& right)
{
    return equal_mod_field(left.x(), right.x()) &
        equal_mod_field(left.y(), right.y());
}

inline blst::Fp canonical_fp(const CurveRing::StandardElement& value)
{
    const blst::Fp montgomery = ring_fp_to_scalar(value);
    blst::Fp canonical = {};
    blst::Fp::from_mont(canonical, montgomery);
    return canonical;
}

inline blst::Fp2 canonical_fp2(const CurveFp2Ring::StandardElement& value)
{
    const blst::Fp2 montgomery = ring_fp2_to_scalar(value);
    blst::Fp2 canonical = {};
    blst::Fp::from_mont(canonical[0], montgomery[0]);
    blst::Fp::from_mont(canonical[1], montgomery[1]);
    return canonical;
}

inline void canonical_fp_to_be(
    unsigned char out[48], const blst::Fp& canonical)
{
    be_bytes_from_limbs(out, canonical.l, sizeof(canonical));
}

inline bool decode_fp_be(
    CurveRing::StandardElement& out,
    const unsigned char in[48],
    bool clear_top_three_bits)
{
    blst::Fp canonical = {};
    limbs_from_be_bytes(canonical.l, in, sizeof(canonical));
    if (clear_top_three_bits) {
        canonical[blst::Fp::num_limb() - 1] &=
            (static_cast<limb_t>(0) - 1) >> 3;
    }

    blst::Fp reduced = {};
    blst::Fp::add(reduced, canonical, blst::Fp::zero());
    if (!vec_is_equal(reduced.l, canonical.l, sizeof(canonical)))
        return false;

    blst::Fp montgomery = {};
    blst::Fp::mul(montgomery, canonical, blst::get_BLST_RR());
    out = scalar_fp_to_ring(montgomery);
    return true;
}

inline void serialize_fp2_be(
    unsigned char out[96],
    const CurveFp2Ring::StandardElement& value)
{
    const blst::Fp2 canonical = canonical_fp2(value);
    canonical_fp_to_be(out, canonical[1]);
    canonical_fp_to_be(out + 48, canonical[0]);
}

inline bool decode_fp2_be(
    CurveFp2Ring::StandardElement& out,
    const unsigned char in[96],
    bool clear_top_three_bits)
{
    CurveRing::StandardElement imaginary = {};
    CurveRing::StandardElement real = {};
    if (!decode_fp_be(imaginary, in, clear_top_three_bits) ||
        !decode_fp_be(real, in + 48, false)) {
        return false;
    }
    out = CurveFp2Ring::StandardElement(real, imaginary);
    return true;
}

inline limb_t fp_sign(const CurveRing::StandardElement& value)
{
    return blst::Fp::sgn0_pty_mont(ring_fp_to_scalar(value));
}

inline limb_t fp2_sign(const CurveFp2Ring::StandardElement& value)
{
    return blst::Fp2::sgn0_pty_mont(ring_fp2_to_scalar(value));
}

inline bool recover_g1_y(
    CurveRing::StandardElement& out,
    const CurveRing::StandardElement& x)
{
    const blst::Fp scalar_x = ring_fp_to_scalar(x);
    blst::Fp scalar_y = {};
    blst::Fp::sqr(scalar_y, scalar_x);
    blst::Fp::mul(scalar_y, scalar_y, scalar_x);
    blst::Fp::add(scalar_y, scalar_y, blst::scalar::get_BLST_B_E1());
    if (!blst::Fp::sqrt(scalar_y, scalar_y))
        return false;
    out = scalar_fp_to_ring(scalar_y);
    return true;
}

inline bool recover_g2_y(
    CurveFp2Ring::StandardElement& out,
    const CurveFp2Ring::StandardElement& x)
{
    const blst::Fp2 scalar_x = ring_fp2_to_scalar(x);
    blst::Fp2 scalar_y = {};
    blst::Fp2::sqr(scalar_y, scalar_x);
    blst::Fp2::mul(scalar_y, scalar_y, scalar_x);
    blst::Fp2::add(scalar_y, scalar_y, blst::scalar::get_BLST_B_E2());
    if (!blst::Fp2::sqrt(scalar_y, scalar_y))
        return false;
    out = scalar_fp2_to_ring(scalar_y);
    return true;
}

inline blst::scalar::G2AffinePoint curve_g2_affine_to_scalar(const CurveG2Affine& point)
{
    blst::scalar::G2AffinePoint out = {};
    out.x = ring_fp2_to_scalar(point.x);
    out.y = ring_fp2_to_scalar(point.y);
    return out;
}

inline blst::scalar::G2Point curve_g2_to_scalar(const CurveG2Projective& point)
{
    const auto z_squared = curve_fp2_ring().mul(point.z, point.z);
    const auto jacobian_x = curve_fp2_ring().mul(point.x, point.z);
    const auto jacobian_y = curve_fp2_ring().mul(point.y, z_squared);
    blst::scalar::G2Point out = {};
    out.x = ring_fp2_to_scalar(jacobian_x);
    out.y = ring_fp2_to_scalar(jacobian_y);
    out.z = ring_fp2_to_scalar(point.z);
    return out;
}

inline std::array<uint8_t, 32> scalar_to_fixed_le_bytes(const byte* scalar, size_t nbits)
{
    std::array<uint8_t, 32> out = {};
    pow256 canonical = {};
    blst::scalar::canonical_scalar_256(canonical, scalar, nbits);
    std::memcpy(out.data(), canonical, sizeof(canonical));
    vec_zero(canonical, sizeof(canonical));
    return out;
}

inline std::array<uint8_t, 32> scalar_to_fixed_le_bytes_vartime(
    const byte* scalar, size_t nbits)
{
    std::array<uint8_t, 32> out = {};
    blst::scalar::canonical_scalar_256_vartime(out.data(), scalar, nbits);
    return out;
}

inline CurveG1Affine g1_affine_infinity()
{
    return {curve_ring().zero(), curve_ring().zero()};
}

inline CurveG2Affine g2_affine_infinity()
{
    return {curve_fp2_ring().zero(), curve_fp2_ring().zero()};
}

inline CurveRing::StandardElement curve_inverse_vartime(
    const CurveRing::StandardElement& value)
{
    return invert_fp_vartime_via_scalar(curve_ring(), value);
}

inline CurveFp2Ring::StandardElement curve_fp2_inverse_vartime(
    const CurveFp2Ring::StandardElement& value)
{
    blst::Fp2 scalar = ring_fp2_to_scalar(value);
    blst::Fp2 inverse = {};
    blst::Fp2::reciprocal(inverse, scalar);
    return scalar_fp2_to_ring(inverse);
}

inline CurveRing::StandardElement curve_secret_inverse(
    const CurveRing::StandardElement& value, bool_t& is_zero_value)
{
    blst::Fp scalar = ring_fp_to_scalar(value);
    is_zero_value = blst::Fp::equal(scalar, blst::Fp::zero());
    blst::Fp inverse = {};
    blst::Fp::flt_reciprocal(inverse, scalar);
    const auto out = scalar_fp_to_ring(inverse);
    vec_zero(&scalar, sizeof(scalar));
    vec_zero(&inverse, sizeof(inverse));
    return out;
}

inline CurveFp2Ring::StandardElement curve_secret_inverse(
    const CurveFp2Ring::StandardElement& value, bool_t& is_zero_value)
{
    blst::Fp2 scalar = ring_fp2_to_scalar(value);
    is_zero_value = blst::Fp::equal(scalar[0], blst::Fp::zero()) &
        blst::Fp::equal(scalar[1], blst::Fp::zero());
    blst::Fp2 inverse = {};
    blst::Fp2::flt_reciprocal(inverse, scalar);
    const auto out = scalar_fp2_to_ring(inverse);
    vec_zero(&scalar, sizeof(scalar));
    vec_zero(&inverse, sizeof(inverse));
    return out;
}

inline CurveG1Projective normalize_secret_projective(const CurveG1Projective& point)
{
    bool_t is_infinity = 0;
    const auto z_inverse = curve_secret_inverse(point.z, is_infinity);
    auto [x, y] = curve_ring().batch_reduce_expand(
        curve_ring().prep_left(point.x) * z_inverse,
        curve_ring().prep_left(point.y) * z_inverse);
    const auto z = curve_ring().one().ctime_select(
        !is_infinity, curve_ring().zero());
    return {x, y, z};
}

inline CurveG2Projective normalize_secret_projective(const CurveG2Projective& point)
{
    bool_t is_infinity = 0;
    const auto z_inverse = curve_secret_inverse(point.z, is_infinity);
    const auto x = curve_fp2_ring().mul(point.x, z_inverse);
    const auto y = curve_fp2_ring().mul(point.y, z_inverse);
    const auto z = curve_fp2_ring().one().ctime_select(
        !is_infinity, curve_fp2_ring().zero());
    return {x, y, z};
}

inline CurveG1Affine projective_to_affine(const CurveG1Projective& point)
{
    if (is_zero_mod_field(point.z)) {
        return g1_affine_infinity();
    }

    const auto z_inv = curve_inverse_vartime(point.z);
    auto [x, y] = curve_ring().batch_reduce_expand(
        curve_ring().prep_left(point.x) * z_inv,
        curve_ring().prep_left(point.y) * z_inv);
    return {x, y};
}

inline CurveG2Affine projective_to_affine(const CurveG2Projective& point)
{
    if (is_zero_mod_field(point.z)) {
        return g2_affine_infinity();
    }

    const auto z_inv = curve_fp2_inverse_vartime(point.z);
    return {
        curve_fp2_ring().mul(point.x, z_inv),
        curve_fp2_ring().mul(point.y, z_inv),
    };
}

inline CurveG1Projective affine_to_projective(const CurveG1Affine& point)
{
    if (is_zero_mod_field(point.x) & is_zero_mod_field(point.y)) {
        return curve_g1().zero(curve_ring());
    }
    return {point.x, point.y, curve_ring().one()};
}

inline CurveG2Projective affine_to_projective(const CurveG2Affine& point)
{
    if (is_zero_mod_field(point.x) & is_zero_mod_field(point.y)) {
        return curve_g2().zero(curve_fp2_ring());
    }
    return {point.x, point.y, curve_fp2_ring().one()};
}

inline bool_t on_curve(const CurveG1Affine& point)
{
    if (is_zero_mod_field(point.x) & is_zero_mod_field(point.y))
        return 1;

    const auto xx = curve_ring().modmul(point.x, point.x);
    const auto xxx = curve_ring().modmul(xx, point.x);
    const auto yy = curve_ring().modmul(point.y, point.y);
    const auto two = curve_ring().one() * std::integral_constant<int, 2>{};
    const auto four = two * std::integral_constant<int, 2>{};
    const auto rhs = NormalizePointCoord(
        xxx + four, curve_ring());
    return equal_mod_field(yy, rhs);
}

inline bool_t on_curve(const CurveG2Affine& point)
{
    if (is_zero_mod_field(point.x) & is_zero_mod_field(point.y))
        return 1;

    const auto xx = curve_fp2_ring().square(point.x);
    const auto xxx = curve_fp2_ring().mul(xx, point.x);
    const auto yy = curve_fp2_ring().square(point.y);
    const auto b = curve_fp2_ring().mul_by_nonresidue(
        curve_fp2_ring().one());
    const auto two_b = b * std::integral_constant<int, 2>{};
    const auto four_b = two_b * std::integral_constant<int, 2>{};
    const auto rhs = NormalizePointCoord(xxx + four_b, curve_fp2_ring());
    return equal_mod_field(yy, rhs);
}

inline bool_t on_curve(
    const CurveG1Projective& point, bool_t& is_infinity)
{
    is_infinity = is_zero_mod_field(point.z);
    if (is_infinity)
        return 1;

    const auto xx = curve_ring().modmul(point.x, point.x);
    const auto yy = curve_ring().modmul(point.y, point.y);
    const auto zz = curve_ring().modmul(point.z, point.z);
    const auto xxx = curve_ring().modmul(xx, point.x);
    const auto yyz = curve_ring().modmul(yy, point.z);
    const auto zzz = curve_ring().modmul(zz, point.z);
    const auto two_zzz = zzz * std::integral_constant<int, 2>{};
    const auto four_zzz = two_zzz * std::integral_constant<int, 2>{};
    const auto rhs = NormalizePointCoord(
        xxx + four_zzz, curve_ring());
    return equal_mod_field(yyz, rhs);
}

inline bool_t on_curve(const CurveG1Projective& point)
{
    bool_t is_infinity = 0;
    return on_curve(point, is_infinity);
}

inline bool_t on_curve(
    const CurveG2Projective& point, bool_t& is_infinity)
{
    is_infinity = is_zero_mod_field(point.z);
    if (is_infinity)
        return 1;

    const auto xx = curve_fp2_ring().square(point.x);
    const auto yy = curve_fp2_ring().square(point.y);
    const auto zz = curve_fp2_ring().square(point.z);
    const auto xxx = curve_fp2_ring().mul(xx, point.x);
    const auto yyz = curve_fp2_ring().mul(yy, point.z);
    const auto zzz = curve_fp2_ring().mul(zz, point.z);
    const auto bzzz = curve_fp2_ring().mul_by_nonresidue(zzz);
    const auto two_bzzz = bzzz * std::integral_constant<int, 2>{};
    const auto four_bzzz = two_bzzz * std::integral_constant<int, 2>{};
    const auto rhs = NormalizePointCoord(xxx + four_bzzz, curve_fp2_ring());
    return equal_mod_field(yyz, rhs);
}

inline bool_t on_curve(const CurveG2Projective& point)
{
    bool_t is_infinity = 0;
    return on_curve(point, is_infinity);
}

inline const CurveRing::StandardElement& aprime_e1()
{
    return kAprimeE1;
}

inline const CurveFp2Ring::StandardElement& aprime_e2()
{
    return kAprimeE2;
}

inline CurveFp2Ring::StandardElement fp2_conjugate(const CurveFp2Ring::StandardElement& value)
{
    return {value.x(), curve_ring().standard_negate(value.y())};
}

inline const CurveFp2Ring::StandardElement& psi_frobenius_x()
{
    return kPsiFrobeniusX;
}

inline const CurveFp2Ring::StandardElement& psi_frobenius_y()
{
    return kPsiFrobeniusY;
}

inline const std::array<byte, 32>& minus_z_bytes()
{
    static const std::array<byte, 32> value = {0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0xd2};
    return value;
}

}  // namespace avx_curve_detail

/** AVX G1 affine storage. The type does not encode subgroup-validation state.
 * Basic single-point operations require on-curve inputs; optimized
 * multiplication and advanced algorithms require prior G1 validation. */
struct G1AffinePoint
{
    using Field_type = avx_curve_detail::CurveRing::StandardElement;

    Field_type X;
    Field_type Y;

    constexpr G1AffinePoint() = default;
    constexpr G1AffinePoint(const Field_type& x, const Field_type& y) : X(x), Y(y) {}

    /// Prime-order G1 subgroup generator.
    static const G1AffinePoint& generator();
    /// Constant-time comparison of on-curve values; this is not a validator.
    static int equal(const G1AffinePoint& a, const G1AffinePoint& b);
    /// Test an on-curve value for infinity; this is not a validator.
    static bool_t is_inf(const G1AffinePoint& p);
    /// Validate the curve equation, accepting infinity.
    static bool_t on_curve(const G1AffinePoint& p);
    /// Validate both the curve equation and prime-order G1 membership.
    static bool_t in_subgroup(const G1AffinePoint& p);
    /// Constant-time conditional negation of an on-curve value.
    /// c must be zero or one; no subgroup validation occurs.
    static void cneg(G1AffinePoint& inOut, bool_t c);
    /// Write 96-byte raw big-endian X || Y for an on-curve public value.
    static limb_t serialize_BE(unsigned char* out, const G1AffinePoint& in);
    static void serialize(unsigned char* out, const G1AffinePoint& in);
    /// Parse 96-byte raw big-endian X || Y and validate range and curve.
    /// Success establishes on-curve status only; on failure out is unchanged.
    static BLST_ERROR deserialize_BE(G1AffinePoint& out, const unsigned char* in);
    /// Parse 96-byte IETF input. Success establishes canonical encoding and
    /// on-curve status only; on failure out is unchanged.
    static BLST_ERROR deserialize(G1AffinePoint& out, const unsigned char* in);
    /// Write 48-byte raw big-endian X and return y parity.
    static limb_t compress_BE(unsigned char* out, const G1AffinePoint& in);
    static void compress(unsigned char* out, const G1AffinePoint& in);
    /// Parse 48-byte raw big-endian X and reconstruct an on-curve value.
    /// Returns parity 0/1 on success or a negated BLST_ERROR in limb_t;
    /// success does not establish G1 membership and failure leaves out unchanged.
    static limb_t uncompress_BE(G1AffinePoint& out, const unsigned char* in);
    /// Parse 48-byte IETF compressed input. Success establishes canonical
    /// encoding and on-curve status only; on failure out is unchanged.
    static BLST_ERROR uncompress(G1AffinePoint& out, const unsigned char* in);
};

/** AVX G1 homogeneous projective storage `(X/Z,Y/Z)`. The type does not encode subgroup-validation state.
 * Unless stated otherwise, basic single-point methods accept on-curve points
 * and preserve, but do not establish, subgroup membership. */
struct G1Point
{
    using Field_type = avx_curve_detail::CurveRing::StandardElement;

    Field_type X;
    Field_type Y;
    Field_type Z;

    constexpr G1Point() = default;
    constexpr G1Point(const Field_type& x, const Field_type& y, const Field_type& z) : X(x), Y(y), Z(z) {}

    /// Complete add-or-double of on-curve points; subgroup membership is not required.
    static void dadd(G1Point& out, const G1Point& a, const G1Point& b);
    static void dadd(G1Point& out, const G1Point& a, const G1Point& b, const Field_type* a4);
    static void dadd(G1Point& out, const G1Point& a, const G1AffinePoint& b);
    /// Fast incomplete addition of on-curve distinct, non-opposite finite inputs.
    static void add(G1Point& out, const G1Point& a, const G1Point& b);
    static void add(G1Point& out, const G1Point& a, const G1AffinePoint& b);
    /// Double an on-curve point; subgroup membership is not required.
    static void doubling(G1Point& out, const G1Point& a);
    /// Constant-time scalar multiplication of a previously validated subgroup
    /// point. Release does not validate membership. For canonical scalars,
    /// timing is independent of contents at fixed public nbits; noncanonical
    /// inputs are reduced modulo r and may reveal that reduction was required.
    /// nbits must not exceed 256 and scalar may be null only for zero bits.
    /// Exact out/a aliasing is supported; partial overlap is not.
    static void mult(G1Point& out, const G1Point& a, const byte* scalar, size_t nbits);
    /// Optimized variable-time multiplication of a previously validated
    /// subgroup point using public scalar data. Release does not validate
    /// membership. nbits must not exceed 256 and scalar may be null only for
    /// zero bits. Exact out/a aliasing is supported; partial overlap is not.
    static void mult_vartime(G1Point& out, const G1Point& a, const byte* scalar, size_t nbits);
    /// Constant-time exact unsigned-scalar multiplication of an on-curve point.
    /// Subgroup membership is not required and no GLV reduction is performed.
    /// nbits must not exceed 256 and scalar may be null only for zero bits.
    /// Exact out/a aliasing is supported; partial overlap is not.
    static void unchecked_mult(G1Point& out, const G1Point& a, const byte* scalar, size_t nbits);
    /// Sign a previously validated subgroup hash point with a canonical,
    /// nonzero secret key accepted by keygen/sk_check. Timing is independent
    /// of SK; exact out/in aliasing is supported.
    static void sign(G1Point& out, const G1Point& in, const pow256 SK);
    /// Sign and optionally emit a 96-byte pre-decorated X || Y buffer and/or
    /// affine subgroup point. Finalize the byte buffer as described for to_pk2.
    /// Either output pointer may be null. Timing is independent of SK.
    static void sign2(unsigned char* out, G1AffinePoint* sig,
        const G1Point& hash, const pow256 SK);
    /// Derive the subgroup public key for a canonical, nonzero secret key.
    /// Timing is independent of SK.
    static void to_pk(G1Point& out, const pow256 SK);
    /// Optionally emit a 96-byte pre-decorated X || Y public-key buffer and/or
    /// affine subgroup point. The byte buffer is not yet a conforming encoding:
    /// set bit 7 and retain its first 48 bytes for compressed form, or clear bit
    /// 5 for uncompressed form. Either output pointer may be null. Timing is
    /// independent of SK.
    static void to_pk2(unsigned char* out, G1AffinePoint* PK, const pow256 SK);
    /// Constant-time conditional negation of an on-curve point.
    /// cbit must be zero or one; no subgroup validation occurs.
    static void cneg(G1Point& p, limb_t cbit);
    /// Negate an on-curve point; exact out/p aliasing is supported.
    static void neg(G1Point& out, const G1Point& p);
    /// Convert an on-curve public point to affine. This is variable-time;
    /// infinity maps to the all-zero affine value and subgroup state is preserved.
    static void to_affine(G1AffinePoint& out, const G1Point& in);
    /// Batch-normalize public points with one field inversion. Output and input
    /// lengths must match; prefix_products supplies at least one field element
    /// per point and may be retained by a caller-owned hot-path workspace.
    static void to_affine_batch(
        std::span<G1AffinePoint> out,
        std::span<const G1Point> in,
        std::span<Field_type> prefix_products);
    /// Convert an on-curve affine value to homogeneous projective storage;
    /// subgroup state is preserved and exact out/a overlap is not permitted.
    static void from_affine(G1Point& out, const G1AffinePoint& a);
    /// Legacy name for variable-time normalization of homogeneous projective
    /// storage. Subgroup state is preserved; exact out/in aliasing is supported.
    static void from_jacobian(G1Point& out, const G1Point& in);
    /// Constant-time conditional swap. Objects must be distinct and cbit 0/1.
    static void cswap(G1Point& a, G1Point& b, bool_t cbit);
    /// Constant-time conditional copy. Objects must be distinct and cbit 0/1.
    static void ccopy(G1Point& a, const G1Point& b, bool_t cbit);
    /// Constant-time comparison of on-curve values; this is not a validator.
    static limb_t equal(const G1Point& a, const G1Point& b);
    /// Prime-order G1 subgroup generator.
    static const G1Point& generator();
    /// Return the homogeneous projective infinity representation.
    static G1Point infinity();
    static G1Point identity();
    static bool_t is_identity(const G1Point& p);
    /// Test an on-curve value for infinity; this is not a validator.
    static bool_t is_inf(const G1Point& p);
    /// Validate the curve equation, accepting infinity.
    static bool_t on_curve(const G1Point& p);
    /// Validate both the curve equation and prime-order G1 membership.
    static bool_t in_subgroup(const G1Point& p);
    /// Write 96-byte raw big-endian affine coordinates for an on-curve public point.
    static limb_t serialize_BE(unsigned char* out, const G1Point& in);
    static void serialize(unsigned char* out, const G1Point& in);
    /// Write 48-byte raw big-endian X and return y parity.
    static limb_t compress_BE(unsigned char* out, const G1Point& in);
    static void compress(unsigned char* out, const G1Point& in);
};

/** AVX G2 affine storage. The type does not encode subgroup-validation state.
 * Basic single-point operations require on-curve inputs; optimized
 * multiplication and advanced algorithms require prior G2 validation. */
struct G2AffinePoint
{
    using Field_type = avx_curve_detail::CurveFp2Ring::StandardElement;

    Field_type X;
    Field_type Y;

    G2AffinePoint() = default;
    constexpr G2AffinePoint(const Field_type& x, const Field_type& y) : X(x), Y(y) {}

    /// Prime-order G2 subgroup generator.
    static const G2AffinePoint& generator();
    /// Constant-time comparison of on-curve values; this is not a validator.
    static int equal(const G2AffinePoint& a, const G2AffinePoint& b);
    /// Test an on-curve value for infinity; this is not a validator.
    static bool_t is_inf(const G2AffinePoint& p);
    /// Validate the curve equation, accepting infinity.
    static bool_t on_curve(const G2AffinePoint& p);
    /// Validate both the curve equation and prime-order G2 membership.
    static bool_t in_subgroup(const G2AffinePoint& p);
    /// Constant-time conditional negation of an on-curve value.
    /// c must be zero or one; no subgroup validation occurs.
    static void cneg(G2AffinePoint& inOut, bool_t c);
    /// Write 192-byte raw big-endian coordinates for an on-curve public value.
    static limb_t serialize_BE(unsigned char* out, const G2AffinePoint& in);
    static void serialize(unsigned char* out, const G2AffinePoint& in);
    /// Parse 192-byte raw big-endian coordinates and validate range and curve.
    /// Success establishes on-curve status only; on failure out is unchanged.
    static BLST_ERROR deserialize_BE(G2AffinePoint& out, const unsigned char* in);
    /// Parse 192-byte IETF input. Success establishes canonical encoding and
    /// on-curve status only; on failure out is unchanged.
    static BLST_ERROR deserialize(G2AffinePoint& out, const unsigned char* in);
    /// Write 96-byte raw big-endian X and return y parity.
    static limb_t compress_BE(unsigned char* out, const G2AffinePoint& in);
    static void compress(unsigned char* out, const G2AffinePoint& in);
    /// Parse 96-byte raw big-endian X and reconstruct an on-curve value.
    /// Returns parity 0/1 on success or a negated BLST_ERROR in limb_t;
    /// success does not establish G2 membership and failure leaves out unchanged.
    static limb_t uncompress_BE(G2AffinePoint& out, const unsigned char* in);
    /// Parse 96-byte IETF compressed input. Success establishes canonical
    /// encoding and on-curve status only; on failure out is unchanged.
    static BLST_ERROR uncompress(G2AffinePoint& out, const unsigned char* in);
};

/** AVX G2 homogeneous projective storage `(X/Z,Y/Z)`. The type does not encode subgroup-validation state.
 * Unless stated otherwise, basic single-point methods accept on-curve points
 * and preserve, but do not establish, subgroup membership. */
struct G2Point
{
    using Field_type = avx_curve_detail::CurveFp2Ring::StandardElement;

    Field_type X;
    Field_type Y;
    Field_type Z;

    G2Point() = default;
    constexpr G2Point(
        const Field_type& x, const Field_type& y, const Field_type& z)
        : X(x), Y(y), Z(z)
    {}

    /// Complete add-or-double of on-curve points; subgroup membership is not required.
    static void dadd(G2Point& out, const G2Point& a, const G2Point& b);
    static void dadd(G2Point& out, const G2Point& a, const G2Point& b, const Field_type* a4);
    static void dadd(G2Point& out, const G2Point& a, const G2AffinePoint& b);
    /// Fast incomplete addition of on-curve distinct, non-opposite finite inputs.
    static void add(G2Point& out, const G2Point& a, const G2Point& b);
    static void add(G2Point& out, const G2Point& a, const G2AffinePoint& b);
    /// Double an on-curve point; subgroup membership is not required.
    static void doubling(G2Point& out, const G2Point& a);
    /// Constant-time scalar multiplication of a previously validated subgroup
    /// point. Release does not validate membership. For canonical scalars,
    /// timing is independent of contents at fixed public nbits; noncanonical
    /// inputs are reduced modulo r and may reveal that reduction was required.
    /// nbits must not exceed 256 and scalar may be null only for zero bits.
    /// Exact out/a aliasing is supported; partial overlap is not.
    static void mult(G2Point& out, const G2Point& a, const byte* scalar, size_t nbits);
    /// Optimized variable-time multiplication of a previously validated
    /// subgroup point using public scalar data. Release does not validate
    /// membership. nbits must not exceed 256 and scalar may be null only for
    /// zero bits. Exact out/a aliasing is supported; partial overlap is not.
    static void mult_vartime(G2Point& out, const G2Point& a, const byte* scalar, size_t nbits);
    /// Constant-time exact unsigned-scalar multiplication of an on-curve point.
    /// Subgroup membership is not required and no GLS reduction is performed.
    /// nbits must not exceed 256 and scalar may be null only for zero bits.
    /// Exact out/a aliasing is supported; partial overlap is not.
    static void unchecked_mult(G2Point& out, const G2Point& a, const byte* scalar, size_t nbits);
    /// Sign a previously validated subgroup hash point with a canonical,
    /// nonzero secret key accepted by keygen/sk_check. Timing is independent
    /// of SK; exact out/in aliasing is supported.
    static void sign(G2Point& out, const G2Point& in, const pow256 SK);
    /// Sign and optionally emit a 192-byte pre-decorated X || Y buffer and/or
    /// affine subgroup point. Finalize the byte buffer as described for to_pk2.
    /// Either output pointer may be null. Timing is independent of SK.
    static void sign2(unsigned char* out, G2AffinePoint* sig,
        const G2Point& hash, const pow256 SK);
    /// Derive the subgroup public key for a canonical, nonzero secret key.
    /// Timing is independent of SK.
    static void to_pk(G2Point& out, const pow256 SK);
    /// Optionally emit a 192-byte pre-decorated X || Y public-key buffer and/or
    /// affine subgroup point. The byte buffer is not yet a conforming encoding:
    /// set bit 7 and retain its first 96 bytes for compressed form, or clear bit
    /// 5 for uncompressed form. Either output pointer may be null. Timing is
    /// independent of SK.
    static void to_pk2(unsigned char* out, G2AffinePoint* PK, const pow256 SK);
    /// Constant-time conditional negation of an on-curve point.
    /// cbit must be zero or one; no subgroup validation occurs.
    static void cneg(G2Point& p, limb_t cbit);
    /// Negate an on-curve point; exact out/p aliasing is supported.
    static void neg(G2Point& out, const G2Point& p);
    /// Convert an on-curve public point to affine. This is variable-time;
    /// infinity maps to the all-zero affine value and subgroup state is preserved.
    static void to_affine(G2AffinePoint& out, const G2Point& in);
    /// Batch-normalize public points with one Fp2 inversion. Output and input
    /// lengths must match; prefix_products supplies at least one field element
    /// per point and may be retained by a caller-owned hot-path workspace.
    static void to_affine_batch(
        std::span<G2AffinePoint> out,
        std::span<const G2Point> in,
        std::span<Field_type> prefix_products);
    /// Convert an on-curve affine value to homogeneous projective storage;
    /// subgroup state is preserved and exact out/a overlap is not permitted.
    static void from_affine(G2Point& out, const G2AffinePoint& a);
    /// Legacy name for variable-time normalization of homogeneous projective
    /// storage. Subgroup state is preserved; exact out/in aliasing is supported.
    static void from_jacobian(G2Point& out, const G2Point& in);
    /// Constant-time conditional swap. Objects must be distinct and cbit 0/1.
    static void cswap(G2Point& a, G2Point& b, bool_t cbit);
    /// Constant-time conditional copy. Objects must be distinct and cbit 0/1.
    static void ccopy(G2Point& a, const G2Point& b, bool_t cbit);
    /// Constant-time comparison of on-curve values; this is not a validator.
    static limb_t equal(const G2Point& a, const G2Point& b);
    /// Prime-order G2 subgroup generator.
    static const G2Point& generator();
    /// Return the homogeneous projective infinity representation.
    static G2Point infinity();
    static G2Point identity();
    static bool_t is_identity(const G2Point& p);
    /// Test an on-curve value for infinity; this is not a validator.
    static bool_t is_inf(const G2Point& p);
    /// Validate the curve equation, accepting infinity.
    static bool_t on_curve(const G2Point& p);
    /// Validate both the curve equation and prime-order G2 membership.
    static bool_t in_subgroup(const G2Point& p);
    /// Write 192-byte raw big-endian affine coordinates for an on-curve public point.
    static limb_t serialize_BE(unsigned char* out, const G2Point& in);
    static void serialize(unsigned char* out, const G2Point& in);
    /// Write 96-byte raw big-endian X and return y parity.
    static limb_t compress_BE(unsigned char* out, const G2Point& in);
    static void compress(unsigned char* out, const G2Point& in);
};

inline G1AffinePoint ctime_select_g1_affine(bool_t selection, const G1AffinePoint& when_true, const G1AffinePoint& when_false)
{
    return {
        when_true.X.ctime_select(selection, when_false.X),
        when_true.Y.ctime_select(selection, when_false.Y),
    };
}

inline G1Point ctime_select_g1(bool_t selection, const G1Point& when_true, const G1Point& when_false)
{
    return {
        when_true.X.ctime_select(selection, when_false.X),
        when_true.Y.ctime_select(selection, when_false.Y),
        when_true.Z.ctime_select(selection, when_false.Z),
    };
}

inline G2AffinePoint ctime_select_g2_affine(bool_t selection, const G2AffinePoint& when_true, const G2AffinePoint& when_false)
{
    return {
        when_true.X.ctime_select(selection, when_false.X),
        when_true.Y.ctime_select(selection, when_false.Y),
    };
}

inline G2Point ctime_select_g2(bool_t selection, const G2Point& when_true, const G2Point& when_false)
{
    return {
        when_true.X.ctime_select(selection, when_false.X),
        when_true.Y.ctime_select(selection, when_false.Y),
        when_true.Z.ctime_select(selection, when_false.Z),
    };
}

namespace avx_curve_detail
{

template<typename Group, size_t Window>
BLST_INLINE void unchecked_mult_window(
    Group& out, const Group& point, const byte* scalar, size_t bits)
{
    static_assert(Window == 4 || Window == 5);
    constexpr size_t table_size = size_t{1} << (Window - 1);
    Group table[table_size];
    table[0] = point;
    Group::doubling(table[1], point);
    for (size_t i = 2, j = 1; i < table_size; i += 2, ++j) {
        // Complete addition is required here: unchecked inputs can have low
        // order, so table entries may coincide or be opposites.
        Group::dadd(table[i], table[j], table[j - 1]);
        Group::doubling(table[i + 1], table[j]);
    }

    const auto gather = [&table](limb_t encoding) {
        const bool_t sign = (encoding >> Window) & 1;
        const limb_t index = encoding & ((limb_t{1} << Window) - 1);
        Group selected = Group::infinity();
        for (size_t i = 1; i <= table_size; ++i) {
            Group::ccopy(
                selected, table[i - 1],
                byte_is_zero(static_cast<byte>(i ^ index)));
        }
        Group::cneg(selected, sign);
        return selected;
    };

    size_t window = bits % Window;
    limb_t mask = (limb_t{1} << (window + 1)) - 1;
    size_t bit_index = bits - window;
    limb_t value = bit_index
        ? blst::scalar::get_wval(scalar, bit_index - 1, window + 1)
        : static_cast<limb_t>(scalar[0]) << 1;
    out = gather(blst::scalar::booth_encode(value & mask, Window));

    while (bit_index != 0) {
        for (size_t i = 0; i < Window; ++i)
            Group::doubling(out, out);

        bit_index -= Window;
        value = bit_index
            ? blst::scalar::get_wval(scalar, bit_index - 1, Window + 1)
            : static_cast<limb_t>(scalar[0]) << 1;
        Group row = gather(blst::scalar::booth_encode(
            value & ((limb_t{1} << (Window + 1)) - 1), Window));
        Group::dadd(out, out, row);
    }
}

template<typename Group>
BLST_INLINE void unchecked_mult(
    Group& out, const Group& point, const byte* scalar, size_t bits)
{
    BLST_ASSERT(bits <= 256);
    BLST_ASSERT(bits == 0 || scalar != nullptr);
    BLST_ASSERT(Group::on_curve(point));
    if (bits == 0) {
        out = Group::infinity();
    } else if (bits < 176) {
        unchecked_mult_window<Group, 4>(out, point, scalar, bits);
    } else {
        unchecked_mult_window<Group, 5>(out, point, scalar, bits);
    }
}

}  // namespace avx_curve_detail

inline bool_t g1_in_subgroup_relation(const G1Point& point);
inline bool_t g2_in_subgroup_relation(const G2Point& point);

inline const G1AffinePoint& G1AffinePoint::generator()
{
    static constexpr G1AffinePoint g(
        avx_curve_detail::literal_standard_element(
            matrixnok::generated::g1_generator_x_m1,
            matrixnok::generated::g1_generator_x_m2),
        avx_curve_detail::literal_standard_element(
            matrixnok::generated::g1_generator_y_m1,
            matrixnok::generated::g1_generator_y_m2));
    return g;
}

inline int G1AffinePoint::equal(const G1AffinePoint& a, const G1AffinePoint& b)
{
    return avx_curve_detail::equal_mod_field(a.X, b.X) &
        avx_curve_detail::equal_mod_field(a.Y, b.Y);
}

inline bool_t G1AffinePoint::is_inf(const G1AffinePoint& p)
{
    return avx_curve_detail::is_zero_mod_field(p.X) &
        avx_curve_detail::is_zero_mod_field(p.Y);
}

inline bool_t G1AffinePoint::on_curve(const G1AffinePoint& p)
{
    return avx_curve_detail::on_curve(
        avx_curve_detail::CurveG1Affine{p.X, p.Y});
}

inline bool_t G1AffinePoint::in_subgroup(const G1AffinePoint& p)
{
    if (!on_curve(p))
        return 0;
    G1Point projective = {};
    G1Point::from_affine(projective, p);
    return g1_in_subgroup_relation(projective);
}

inline void G1AffinePoint::cneg(G1AffinePoint& inOut, bool_t c)
{
    BLST_ASSERT(c <= 1);
    G1AffinePoint neg = inOut;
    neg.Y = avx_curve_detail::curve_ring().standard_negate(neg.Y);
    neg = ctime_select_g1_affine(is_inf(inOut), inOut, neg);
    inOut = ctime_select_g1_affine(c, neg, inOut);
}

inline limb_t G1AffinePoint::serialize_BE(unsigned char* out, const G1AffinePoint& in)
{
    const blst::Fp x = avx_curve_detail::canonical_fp(in.X);
    const blst::Fp y = avx_curve_detail::canonical_fp(in.Y);
    avx_curve_detail::canonical_fp_to_be(out, x);
    avx_curve_detail::canonical_fp_to_be(out + 48, y);
    return blst::Fp::sgn0_pty(y);
}

inline void G1AffinePoint::serialize(unsigned char* out, const G1AffinePoint& in)
{
    if (is_inf(in)) {
        bytes_zero(out, 96);
        out[0] = 0x40;
    }
    else {
        (void)serialize_BE(out, in);
    }
}

inline BLST_ERROR G1AffinePoint::deserialize_BE(G1AffinePoint& out, const unsigned char* in)
{
    G1AffinePoint decoded = {};
    if (!avx_curve_detail::decode_fp_be(decoded.X, in, true) ||
        !avx_curve_detail::decode_fp_be(decoded.Y, in + 48, false)) {
        return BLST_BAD_ENCODING;
    }
    if (!on_curve(decoded))
        return BLST_POINT_NOT_ON_CURVE;
    if (avx_curve_detail::is_zero_mod_field(decoded.X))
        return BLST_POINT_NOT_IN_GROUP;
    out = decoded;
    return BLST_SUCCESS;
}

inline BLST_ERROR G1AffinePoint::deserialize(G1AffinePoint& out, const unsigned char* in)
{
    const unsigned char in0 = in[0];
    if ((in0 & 0xe0) == 0)
        return deserialize_BE(out, in);
    if (in0 & 0x80)
        return uncompress(out, in);
    if (in0 & 0x40) {
        if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in + 1, 95)) {
            out = {};
            return BLST_SUCCESS;
        }
    }
    return BLST_BAD_ENCODING;
}

inline limb_t G1AffinePoint::compress_BE(unsigned char* out, const G1AffinePoint& in)
{
    avx_curve_detail::canonical_fp_to_be(
        out, avx_curve_detail::canonical_fp(in.X));
    return avx_curve_detail::fp_sign(in.Y);
}

inline void G1AffinePoint::compress(unsigned char* out, const G1AffinePoint& in)
{
    if (is_inf(in)) {
        bytes_zero(out, 48);
        out[0] = 0xc0;
    }
    else {
        const limb_t sign = compress_BE(out, in);
        out[0] |= static_cast<unsigned char>(0x80 | ((sign & 2) << 4));
    }
}

inline limb_t G1AffinePoint::uncompress_BE(G1AffinePoint& out, const unsigned char* in)
{
    G1AffinePoint decoded = {};
    if (!avx_curve_detail::decode_fp_be(decoded.X, in, true))
        return static_cast<limb_t>(0) - BLST_BAD_ENCODING;
    if (!avx_curve_detail::recover_g1_y(decoded.Y, decoded.X))
        return static_cast<limb_t>(0) - BLST_POINT_NOT_ON_CURVE;
    out = decoded;
    return avx_curve_detail::fp_sign(out.Y);
}

inline BLST_ERROR G1AffinePoint::uncompress(G1AffinePoint& out, const unsigned char* in)
{
    const unsigned char in0 = in[0];
    if ((in0 & 0x80) == 0)
        return BLST_BAD_ENCODING;
    if (in0 & 0x40) {
        if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in + 1, 47)) {
            out = {};
            return BLST_SUCCESS;
        }
        return BLST_BAD_ENCODING;
    }

    G1AffinePoint decoded = {};
    limb_t sign = uncompress_BE(decoded, in);
    if (sign > 3)
        return static_cast<BLST_ERROR>(0 - sign);
    sign = (sign >> 1) ^ ((in0 & 0x20) >> 5);
    cneg(decoded, sign);
    if (avx_curve_detail::is_zero_mod_field(decoded.X))
        return BLST_POINT_NOT_IN_GROUP;
    out = decoded;
    return BLST_SUCCESS;
}

inline void G1Point::dadd(G1Point& out, const G1Point& a, const G1Point& b) { dadd(out, a, b, nullptr); }
inline void G1Point::dadd(G1Point& out, const G1Point& a, const G1Point& b, const Field_type* a4)
{
    const auto raw = avx_curve_detail::curve_g1().dadd({a.X, a.Y, a.Z}, {b.X, b.Y, b.Z}, avx_curve_detail::curve_ring(), a4);
    out = {raw.x, raw.y, raw.z};
}
inline void G1Point::dadd(G1Point& out, const G1Point& a, const G1AffinePoint& b) { add(out, a, b); }

inline void G1Point::add(G1Point& out, const G1Point& a, const G1Point& b)
{
    const auto raw = avx_curve_detail::curve_g1().add_point({a.X, a.Y, a.Z}, {b.X, b.Y, b.Z}, avx_curve_detail::curve_ring());
    out = {raw.x, raw.y, raw.z};
}

inline void G1Point::add(G1Point& out, const G1Point& a, const G1AffinePoint& b)
{
    const auto raw = avx_curve_detail::curve_g1().add_mixed_point({a.X, a.Y, a.Z}, {b.X, b.Y}, avx_curve_detail::curve_ring());
    out = {raw.x, raw.y, raw.z};
}

inline void G1Point::doubling(G1Point& out, const G1Point& a)
{
    const auto raw = avx_curve_detail::curve_g1().double_point({a.X, a.Y, a.Z}, avx_curve_detail::curve_ring());
    out = {raw.x, raw.y, raw.z};
}

inline void G1Point::mult(G1Point& out, const G1Point& a, const byte* scalar, size_t nbits)
{
    BLST_ASSERT(nbits <= 256);
    BLST_ASSERT(nbits == 0 || scalar != nullptr);
    BLST_ASSERT(on_curve(a));
    BLST_ASSERT(in_subgroup(a));
    auto scalar_bytes = avx_curve_detail::scalar_to_fixed_le_bytes(scalar, nbits);
    const auto raw = avx_curve_detail::curve_g1_glv().multiply({a.X, a.Y, a.Z}, scalar_bytes, avx_curve_detail::curve_ring());
    out = {raw.x, raw.y, raw.z};
    vec_zero(scalar_bytes.data(), scalar_bytes.size());
}

inline void G1Point::mult_vartime(G1Point& out, const G1Point& a, const byte* scalar, size_t nbits)
{
    BLST_ASSERT(nbits <= 256);
    BLST_ASSERT(nbits == 0 || scalar != nullptr);
    BLST_ASSERT(on_curve(a));
    BLST_ASSERT(in_subgroup(a));
    auto scalar_bytes = avx_curve_detail::scalar_to_fixed_le_bytes_vartime(scalar, nbits);
    const auto raw = avx_curve_detail::curve_g1_glv().multiply_vartime(
        {a.X, a.Y, a.Z}, scalar_bytes, avx_curve_detail::curve_ring());
    out = {raw.x, raw.y, raw.z};
}

inline void G1Point::unchecked_mult(G1Point& out, const G1Point& a, const byte* scalar, size_t nbits)
{
    avx_curve_detail::unchecked_mult(out, a, scalar, nbits);
}

inline void G1Point::sign(G1Point& out, const G1Point& in, const pow256 SK)
{
    BLST_ASSERT(in_subgroup(in));
    BLST_ASSERT(blst::scalar::sk_check(SK));
    std::array<uint8_t, 32> scalar_bytes = {};
    std::memcpy(scalar_bytes.data(), SK, scalar_bytes.size());
    const auto product = avx_curve_detail::curve_g1_glv().multiply(
        {in.X, in.Y, in.Z}, scalar_bytes, avx_curve_detail::curve_ring());
    vec_zero(scalar_bytes.data(), scalar_bytes.size());
    const auto normalized = avx_curve_detail::normalize_secret_projective(product);
    out = {normalized.x, normalized.y, normalized.z};
}

inline void G1Point::sign2(
    unsigned char* out, G1AffinePoint* sig,
    const G1Point& hash, const pow256 SK)
{
    G1Point point = {};
    sign(point, hash, SK);
    const G1AffinePoint affine = {point.X, point.Y};
    if (sig != nullptr)
        *sig = affine;
    if (out != nullptr) {
        const limb_t parity = G1AffinePoint::serialize_BE(out, affine);
        out[0] |= (parity & 2) << 4;
        out[0] |= is_identity(point) << 6;
    }
}

inline void G1Point::to_pk(G1Point& out, const pow256 SK)
{
    sign(out, generator(), SK);
}

inline void G1Point::to_pk2(unsigned char* out, G1AffinePoint* PK, const pow256 SK)
{
    sign2(out, PK, generator(), SK);
}

inline void G1Point::cneg(G1Point& p, limb_t cbit)
{
    BLST_ASSERT(cbit <= 1);
    G1Point neg = p;
    neg.Y = avx_curve_detail::curve_ring().standard_negate(neg.Y);
    p = ctime_select_g1(cbit, neg, p);
}

inline void G1Point::neg(G1Point& out, const G1Point& p)
{
    out = p;
    out.Y = avx_curve_detail::curve_ring().standard_negate(out.Y);
}

inline void G1Point::to_affine(G1AffinePoint& out, const G1Point& in)
{
    const auto affine = avx_curve_detail::projective_to_affine({in.X, in.Y, in.Z});
    out = {affine.x, affine.y};
}

inline void G1Point::to_affine_batch(
    std::span<G1AffinePoint> out,
    std::span<const G1Point> in,
    std::span<Field_type> prefix_products)
{
    BLST_ASSERT(out.size() == in.size());
    BLST_ASSERT(prefix_products.size() >= in.size());
    if (in.empty()) {
        return;
    }

    const auto& ring = avx_curve_detail::curve_ring();
    Field_type product = ring.one();
    std::size_t nonzero_count = 0;
    for (std::size_t i = 0; i < in.size(); ++i) {
        prefix_products[i] = product;
        if (!avx_curve_detail::is_zero_mod_field(in[i].Z)) {
            product = ring.modmul(product, in[i].Z);
            ++nonzero_count;
        }
    }

    if (nonzero_count == 0) {
        const auto infinity = avx_curve_detail::g1_affine_infinity();
        std::fill(out.begin(), out.end(), G1AffinePoint{infinity.x, infinity.y});
        return;
    }

    Field_type running_inverse =
        avx_curve_detail::curve_inverse_vartime(product);
    for (std::size_t i = in.size(); i-- > 0;) {
        if (avx_curve_detail::is_zero_mod_field(in[i].Z)) {
            const auto infinity = avx_curve_detail::g1_affine_infinity();
            out[i] = {infinity.x, infinity.y};
            continue;
        }

        const Field_type z_inverse =
            ring.modmul(prefix_products[i], running_inverse);
        running_inverse = ring.modmul(running_inverse, in[i].Z);
        auto [x, y] = ring.batch_reduce_expand(
            ring.prep_left(in[i].X) * z_inverse,
            ring.prep_left(in[i].Y) * z_inverse);
        out[i] = {x, y};
    }
}

inline void G1Point::from_affine(G1Point& out, const G1AffinePoint& a)
{
    const auto point = avx_curve_detail::affine_to_projective({a.X, a.Y});
    out = {point.x, point.y, point.z};
}

inline void G1Point::from_jacobian(G1Point& out, const G1Point& in)
{
    G1AffinePoint affine = {};
    to_affine(affine, in);
    from_affine(out, affine);
}

inline void G1Point::cswap(G1Point& a, G1Point& b, bool_t cbit)
{
    BLST_ASSERT(&a != &b);
    BLST_ASSERT(cbit <= 1);
    const G1Point aa = a;
    const G1Point bb = b;
    a = ctime_select_g1(cbit, bb, aa);
    b = ctime_select_g1(cbit, aa, bb);
}

inline void G1Point::ccopy(G1Point& a, const G1Point& b, bool_t cbit)
{
    BLST_ASSERT(&a != &b);
    BLST_ASSERT(cbit <= 1);
    a = ctime_select_g1(cbit, b, a);
}

inline limb_t G1Point::equal(const G1Point& a, const G1Point& b)
{
    return ProjectiveEqual(
        avx_curve_detail::CurveG1Projective{a.X, a.Y, a.Z},
        avx_curve_detail::CurveG1Projective{b.X, b.Y, b.Z},
        avx_curve_detail::curve_ring(),
        [](const Field_type& value) {
            return avx_curve_detail::is_zero_mod_field(value);
        });
}

inline const G1Point& G1Point::generator()
{
    static constexpr G1Point g{
        avx_curve_detail::literal_standard_element(
            matrixnok::generated::g1_generator_x_m1,
            matrixnok::generated::g1_generator_x_m2),
        avx_curve_detail::literal_standard_element(
            matrixnok::generated::g1_generator_y_m1,
            matrixnok::generated::g1_generator_y_m2),
        avx_curve_detail::literal_standard_element(
            matrixnok::generated::one_m1,
            matrixnok::generated::one_m2)};
    return g;
}

inline G1Point G1Point::infinity()
{
    const auto raw = avx_curve_detail::curve_g1().zero(avx_curve_detail::curve_ring());
    return {raw.x, raw.y, raw.z};
}

inline G1Point G1Point::identity() { return infinity(); }
inline bool_t G1Point::is_identity(const G1Point& p)
{
    return avx_curve_detail::is_zero_mod_field(p.Z);
}
inline bool_t G1Point::is_inf(const G1Point& p) { return is_identity(p); }

inline bool_t G1Point::on_curve(const G1Point& p)
{
    return avx_curve_detail::on_curve(
        avx_curve_detail::CurveG1Projective{p.X, p.Y, p.Z});
}

inline bool_t G1Point::in_subgroup(const G1Point& p)
{
    return on_curve(p) && g1_in_subgroup_relation(p);
}

inline limb_t G1Point::serialize_BE(unsigned char* out, const G1Point& in)
{
    G1AffinePoint affine = {};
    to_affine(affine, in);
    return G1AffinePoint::serialize_BE(out, affine);
}

inline void G1Point::serialize(unsigned char* out, const G1Point& in)
{
    G1AffinePoint affine = {};
    to_affine(affine, in);
    G1AffinePoint::serialize(out, affine);
}

inline limb_t G1Point::compress_BE(unsigned char* out, const G1Point& in)
{
    G1AffinePoint affine = {};
    to_affine(affine, in);
    return G1AffinePoint::compress_BE(out, affine);
}

inline void G1Point::compress(unsigned char* out, const G1Point& in)
{
    G1AffinePoint affine = {};
    to_affine(affine, in);
    G1AffinePoint::compress(out, affine);
}

inline const G2AffinePoint& G2AffinePoint::generator()
{
    static constexpr G2AffinePoint g(
        avx_curve_detail::literal_fp2_element(
            matrixnok::generated::g2_generator_x0_m1,
            matrixnok::generated::g2_generator_x0_m2,
            matrixnok::generated::g2_generator_x1_m1,
            matrixnok::generated::g2_generator_x1_m2),
        avx_curve_detail::literal_fp2_element(
            matrixnok::generated::g2_generator_y0_m1,
            matrixnok::generated::g2_generator_y0_m2,
            matrixnok::generated::g2_generator_y1_m1,
            matrixnok::generated::g2_generator_y1_m2));
    return g;
}

inline int G2AffinePoint::equal(const G2AffinePoint& a, const G2AffinePoint& b)
{
    return avx_curve_detail::equal_mod_field(a.X, b.X) &
        avx_curve_detail::equal_mod_field(a.Y, b.Y);
}

inline bool_t G2AffinePoint::is_inf(const G2AffinePoint& p)
{
    return avx_curve_detail::is_zero_mod_field(p.X) &
        avx_curve_detail::is_zero_mod_field(p.Y);
}

inline bool_t G2AffinePoint::on_curve(const G2AffinePoint& p)
{
    return avx_curve_detail::on_curve(
        avx_curve_detail::CurveG2Affine{p.X, p.Y});
}

inline bool_t G2AffinePoint::in_subgroup(const G2AffinePoint& p)
{
    if (!on_curve(p))
        return 0;
    G2Point projective = {};
    G2Point::from_affine(projective, p);
    return g2_in_subgroup_relation(projective);
}

inline void G2AffinePoint::cneg(G2AffinePoint& inOut, bool_t c)
{
    BLST_ASSERT(c <= 1);
    G2AffinePoint neg = inOut;
    neg.Y = avx_curve_detail::curve_fp2_ring().standard_negate(neg.Y);
    neg = ctime_select_g2_affine(is_inf(inOut), inOut, neg);
    inOut = ctime_select_g2_affine(c, neg, inOut);
}

inline limb_t G2AffinePoint::serialize_BE(unsigned char* out, const G2AffinePoint& in)
{
    const blst::Fp2 x = avx_curve_detail::canonical_fp2(in.X);
    const blst::Fp2 y = avx_curve_detail::canonical_fp2(in.Y);
    avx_curve_detail::canonical_fp_to_be(out, x[1]);
    avx_curve_detail::canonical_fp_to_be(out + 48, x[0]);
    avx_curve_detail::canonical_fp_to_be(out + 96, y[1]);
    avx_curve_detail::canonical_fp_to_be(out + 144, y[0]);
    return blst::Fp2::sgn0_pty(y);
}

inline void G2AffinePoint::serialize(unsigned char* out, const G2AffinePoint& in)
{
    if (is_inf(in)) {
        bytes_zero(out, 192);
        out[0] = 0x40;
    }
    else {
        (void)serialize_BE(out, in);
    }
}

inline BLST_ERROR G2AffinePoint::deserialize_BE(G2AffinePoint& out, const unsigned char* in)
{
    G2AffinePoint decoded = {};
    if (!avx_curve_detail::decode_fp2_be(decoded.X, in, true) ||
        !avx_curve_detail::decode_fp2_be(decoded.Y, in + 96, false)) {
        return BLST_BAD_ENCODING;
    }
    if (!on_curve(decoded))
        return BLST_POINT_NOT_ON_CURVE;
    out = decoded;
    return BLST_SUCCESS;
}

inline BLST_ERROR G2AffinePoint::deserialize(G2AffinePoint& out, const unsigned char* in)
{
    const unsigned char in0 = in[0];
    if ((in0 & 0xe0) == 0)
        return deserialize_BE(out, in);
    if (in0 & 0x80)
        return uncompress(out, in);
    if (in0 & 0x40) {
        if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in + 1, 191)) {
            out = {};
            return BLST_SUCCESS;
        }
    }
    return BLST_BAD_ENCODING;
}

inline limb_t G2AffinePoint::compress_BE(unsigned char* out, const G2AffinePoint& in)
{
    avx_curve_detail::serialize_fp2_be(out, in.X);
    return avx_curve_detail::fp2_sign(in.Y);
}

inline void G2AffinePoint::compress(unsigned char* out, const G2AffinePoint& in)
{
    if (is_inf(in)) {
        bytes_zero(out, 96);
        out[0] = 0xc0;
    }
    else {
        const limb_t sign = compress_BE(out, in);
        out[0] |= static_cast<unsigned char>(0x80 | ((sign & 2) << 4));
    }
}

inline limb_t G2AffinePoint::uncompress_BE(G2AffinePoint& out, const unsigned char* in)
{
    G2AffinePoint decoded = {};
    if (!avx_curve_detail::decode_fp2_be(decoded.X, in, true))
        return static_cast<limb_t>(0) - BLST_BAD_ENCODING;
    if (!avx_curve_detail::recover_g2_y(decoded.Y, decoded.X))
        return static_cast<limb_t>(0) - BLST_POINT_NOT_ON_CURVE;
    out = decoded;
    return avx_curve_detail::fp2_sign(out.Y);
}

inline BLST_ERROR G2AffinePoint::uncompress(G2AffinePoint& out, const unsigned char* in)
{
    const unsigned char in0 = in[0];
    if ((in0 & 0x80) == 0)
        return BLST_BAD_ENCODING;
    if (in0 & 0x40) {
        if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in + 1, 95)) {
            out = {};
            return BLST_SUCCESS;
        }
        return BLST_BAD_ENCODING;
    }

    G2AffinePoint decoded = {};
    limb_t sign = uncompress_BE(decoded, in);
    if (sign > 3)
        return static_cast<BLST_ERROR>(0 - sign);
    sign = (sign >> 1) ^ ((in0 & 0x20) >> 5);
    cneg(decoded, sign);
    if (avx_curve_detail::is_zero_mod_field(decoded.X))
        return BLST_POINT_NOT_IN_GROUP;
    out = decoded;
    return BLST_SUCCESS;
}

inline void G2Point::dadd(G2Point& out, const G2Point& a, const G2Point& b) { dadd(out, a, b, nullptr); }
inline void G2Point::dadd(G2Point& out, const G2Point& a, const G2Point& b, const Field_type* a4)
{
    const auto raw = avx_curve_detail::curve_g2().dadd({a.X, a.Y, a.Z}, {b.X, b.Y, b.Z}, avx_curve_detail::curve_fp2_ring(), a4);
    out = {raw.x, raw.y, raw.z};
}
inline void G2Point::dadd(G2Point& out, const G2Point& a, const G2AffinePoint& b) { add(out, a, b); }

inline void G2Point::add(G2Point& out, const G2Point& a, const G2Point& b)
{
    const auto raw = avx_curve_detail::curve_g2().add_point({a.X, a.Y, a.Z}, {b.X, b.Y, b.Z}, avx_curve_detail::curve_fp2_ring());
    out = {raw.x, raw.y, raw.z};
}

inline void G2Point::add(G2Point& out, const G2Point& a, const G2AffinePoint& b)
{
    const auto raw = avx_curve_detail::curve_g2().add_mixed_point({a.X, a.Y, a.Z}, {b.X, b.Y}, avx_curve_detail::curve_fp2_ring());
    out = {raw.x, raw.y, raw.z};
}

inline void G2Point::doubling(G2Point& out, const G2Point& a)
{
    const auto raw = avx_curve_detail::curve_g2().double_point({a.X, a.Y, a.Z}, avx_curve_detail::curve_fp2_ring());
    out = {raw.x, raw.y, raw.z};
}

inline void G2Point::mult(G2Point& out, const G2Point& a, const byte* scalar, size_t nbits)
{
    BLST_ASSERT(nbits <= 256);
    BLST_ASSERT(nbits == 0 || scalar != nullptr);
    BLST_ASSERT(on_curve(a));
    BLST_ASSERT(in_subgroup(a));
    auto scalar_bytes = avx_curve_detail::scalar_to_fixed_le_bytes(scalar, nbits);
    const auto raw = avx_curve_detail::curve_g2_gls().multiply(
        {a.X, a.Y, a.Z},
        scalar_bytes,
        avx_curve_detail::curve_fp2_ring(),
        avx_curve_detail::curve_g2(),
        avx_curve_detail::curve_ring());
    out = {raw.x, raw.y, raw.z};
    vec_zero(scalar_bytes.data(), scalar_bytes.size());
}

inline void G2Point::mult_vartime(G2Point& out, const G2Point& a, const byte* scalar, size_t nbits)
{
    BLST_ASSERT(nbits <= 256);
    BLST_ASSERT(nbits == 0 || scalar != nullptr);
    BLST_ASSERT(on_curve(a));
    BLST_ASSERT(in_subgroup(a));
    auto scalar_bytes = avx_curve_detail::scalar_to_fixed_le_bytes_vartime(scalar, nbits);
    const auto raw = avx_curve_detail::curve_g2_gls().multiply_vartime(
        {a.X, a.Y, a.Z},
        scalar_bytes,
        avx_curve_detail::curve_fp2_ring(),
        avx_curve_detail::curve_g2(),
        avx_curve_detail::curve_ring());
    out = {raw.x, raw.y, raw.z};
}

inline void G2Point::unchecked_mult(G2Point& out, const G2Point& a, const byte* scalar, size_t nbits)
{
    avx_curve_detail::unchecked_mult(out, a, scalar, nbits);
}

inline void G2Point::sign(G2Point& out, const G2Point& in, const pow256 SK)
{
    BLST_ASSERT(in_subgroup(in));
    BLST_ASSERT(blst::scalar::sk_check(SK));
    std::array<uint8_t, 32> scalar_bytes = {};
    std::memcpy(scalar_bytes.data(), SK, scalar_bytes.size());
    const auto product = avx_curve_detail::curve_g2_gls().multiply(
        {in.X, in.Y, in.Z}, scalar_bytes,
        avx_curve_detail::curve_fp2_ring(), avx_curve_detail::curve_g2(),
        avx_curve_detail::curve_ring());
    vec_zero(scalar_bytes.data(), scalar_bytes.size());
    const auto normalized = avx_curve_detail::normalize_secret_projective(product);
    out = {normalized.x, normalized.y, normalized.z};
}

inline void G2Point::sign2(
    unsigned char* out, G2AffinePoint* sig,
    const G2Point& hash, const pow256 SK)
{
    G2Point point = {};
    sign(point, hash, SK);
    const G2AffinePoint affine = {point.X, point.Y};
    if (sig != nullptr)
        *sig = affine;
    if (out != nullptr) {
        const limb_t parity = G2AffinePoint::serialize_BE(out, affine);
        out[0] |= (parity & 2) << 4;
        out[0] |= is_identity(point) << 6;
    }
}

inline void G2Point::to_pk(G2Point& out, const pow256 SK)
{
    sign(out, generator(), SK);
}

inline void G2Point::to_pk2(unsigned char* out, G2AffinePoint* PK, const pow256 SK)
{
    sign2(out, PK, generator(), SK);
}

inline void G2Point::cneg(G2Point& p, limb_t cbit)
{
    BLST_ASSERT(cbit <= 1);
    G2Point neg = p;
    neg.Y = avx_curve_detail::curve_fp2_ring().standard_negate(neg.Y);
    p = ctime_select_g2(cbit, neg, p);
}

inline void G2Point::neg(G2Point& out, const G2Point& p)
{
    out = p;
    out.Y = avx_curve_detail::curve_fp2_ring().standard_negate(out.Y);
}

inline void G2Point::to_affine(G2AffinePoint& out, const G2Point& in)
{
    const auto affine = avx_curve_detail::projective_to_affine({in.X, in.Y, in.Z});
    out = {affine.x, affine.y};
}

inline void G2Point::to_affine_batch(
    std::span<G2AffinePoint> out,
    std::span<const G2Point> in,
    std::span<Field_type> prefix_products)
{
    BLST_ASSERT(out.size() == in.size());
    BLST_ASSERT(prefix_products.size() >= in.size());
    if (in.empty()) {
        return;
    }

    const auto& ring = avx_curve_detail::curve_fp2_ring();
    Field_type product = ring.one();
    std::size_t nonzero_count = 0;
    for (std::size_t i = 0; i < in.size(); ++i) {
        prefix_products[i] = product;
        if (!avx_curve_detail::is_zero_mod_field(in[i].Z)) {
            product = ring.mul(product, in[i].Z);
            ++nonzero_count;
        }
    }

    if (nonzero_count == 0) {
        const auto infinity = avx_curve_detail::g2_affine_infinity();
        std::fill(out.begin(), out.end(), G2AffinePoint{infinity.x, infinity.y});
        return;
    }

    Field_type running_inverse =
        avx_curve_detail::curve_fp2_inverse_vartime(product);
    for (std::size_t i = in.size(); i-- > 0;) {
        if (avx_curve_detail::is_zero_mod_field(in[i].Z)) {
            const auto infinity = avx_curve_detail::g2_affine_infinity();
            out[i] = {infinity.x, infinity.y};
            continue;
        }

        const Field_type z_inverse = ring.mul(prefix_products[i], running_inverse);
        running_inverse = ring.mul(running_inverse, in[i].Z);
        out[i] = {
            ring.mul(in[i].X, z_inverse),
            ring.mul(in[i].Y, z_inverse),
        };
    }
}

inline void G2Point::from_affine(G2Point& out, const G2AffinePoint& a)
{
    const auto point = avx_curve_detail::affine_to_projective({a.X, a.Y});
    out = {point.x, point.y, point.z};
}

inline void G2Point::from_jacobian(G2Point& out, const G2Point& in)
{
    G2AffinePoint affine = {};
    to_affine(affine, in);
    from_affine(out, affine);
}

inline void G2Point::cswap(G2Point& a, G2Point& b, bool_t cbit)
{
    BLST_ASSERT(&a != &b);
    BLST_ASSERT(cbit <= 1);
    const G2Point aa = a;
    const G2Point bb = b;
    a = ctime_select_g2(cbit, bb, aa);
    b = ctime_select_g2(cbit, aa, bb);
}

inline void G2Point::ccopy(G2Point& a, const G2Point& b, bool_t cbit)
{
    BLST_ASSERT(&a != &b);
    BLST_ASSERT(cbit <= 1);
    a = ctime_select_g2(cbit, b, a);
}

inline limb_t G2Point::equal(const G2Point& a, const G2Point& b)
{
    return ProjectiveEqual(
        avx_curve_detail::CurveG2Projective{a.X, a.Y, a.Z},
        avx_curve_detail::CurveG2Projective{b.X, b.Y, b.Z},
        avx_curve_detail::curve_fp2_ring(),
        [](const Field_type& value) {
            return avx_curve_detail::is_zero_mod_field(value);
        });
}

inline const G2Point& G2Point::generator()
{
    static constexpr G2Point g{
        avx_curve_detail::literal_fp2_element(
            matrixnok::generated::g2_generator_x0_m1,
            matrixnok::generated::g2_generator_x0_m2,
            matrixnok::generated::g2_generator_x1_m1,
            matrixnok::generated::g2_generator_x1_m2),
        avx_curve_detail::literal_fp2_element(
            matrixnok::generated::g2_generator_y0_m1,
            matrixnok::generated::g2_generator_y0_m2,
            matrixnok::generated::g2_generator_y1_m1,
            matrixnok::generated::g2_generator_y1_m2),
        avx_curve_detail::literal_fp2_element(
            matrixnok::generated::one_m1,
            matrixnok::generated::one_m2,
            matrixnok::generated::aprime_e2_x_m1,
            matrixnok::generated::aprime_e2_x_m2)};
    return g;
}

inline G2Point G2Point::infinity()
{
    const auto raw = avx_curve_detail::curve_g2().zero(avx_curve_detail::curve_fp2_ring());
    return {raw.x, raw.y, raw.z};
}

inline G2Point G2Point::identity() { return infinity(); }
inline bool_t G2Point::is_identity(const G2Point& p)
{
    return avx_curve_detail::is_zero_mod_field(p.Z);
}
inline bool_t G2Point::is_inf(const G2Point& p) { return is_identity(p); }

inline bool_t G2Point::on_curve(const G2Point& p)
{
    return avx_curve_detail::on_curve(
        avx_curve_detail::CurveG2Projective{p.X, p.Y, p.Z});
}

inline bool_t G2Point::in_subgroup(const G2Point& p)
{
    return on_curve(p) && g2_in_subgroup_relation(p);
}

inline limb_t G2Point::serialize_BE(unsigned char* out, const G2Point& in)
{
    G2AffinePoint affine = {};
    to_affine(affine, in);
    return G2AffinePoint::serialize_BE(out, affine);
}

inline void G2Point::serialize(unsigned char* out, const G2Point& in)
{
    G2AffinePoint affine = {};
    to_affine(affine, in);
    G2AffinePoint::serialize(out, affine);
}

inline limb_t G2Point::compress_BE(unsigned char* out, const G2Point& in)
{
    G2AffinePoint affine = {};
    to_affine(affine, in);
    return G2AffinePoint::compress_BE(out, affine);
}

inline void G2Point::compress(unsigned char* out, const G2Point& in)
{
    G2AffinePoint affine = {};
    to_affine(affine, in);
    G2AffinePoint::compress(out, affine);
}

inline const G1Point::Field_type& get_BLST_Aprime_E1()
{
    return avx_curve_detail::aprime_e1();
}

inline const G2Point::Field_type& get_BLST_Aprime_E2()
{
    return avx_curve_detail::aprime_e2();
}

inline const G1Point& get_BLST_G1()
{
    return G1Point::generator();
}

inline const G2Point& get_BLST_G2()
{
    return G2Point::generator();
}

inline const G1Point& get_BLST_NEG_G1()
{
    static constexpr G1Point neg{
        avx_curve_detail::literal_standard_element(
            matrixnok::generated::g1_generator_x_m1,
            matrixnok::generated::g1_generator_x_m2),
        avx_curve_detail::literal_standard_element(
            matrixnok::generated::g1_generator_neg_y_m1,
            matrixnok::generated::g1_generator_neg_y_m2),
        avx_curve_detail::literal_standard_element(
            matrixnok::generated::one_m1,
            matrixnok::generated::one_m2)};
    return neg;
}

inline const G2Point& get_BLST_NEG_G2()
{
    static constexpr G2Point neg{
        avx_curve_detail::literal_fp2_element(
            matrixnok::generated::g2_generator_x0_m1,
            matrixnok::generated::g2_generator_x0_m2,
            matrixnok::generated::g2_generator_x1_m1,
            matrixnok::generated::g2_generator_x1_m2),
        avx_curve_detail::literal_fp2_element(
            matrixnok::generated::g2_generator_neg_y0_m1,
            matrixnok::generated::g2_generator_neg_y0_m2,
            matrixnok::generated::g2_generator_neg_y1_m1,
            matrixnok::generated::g2_generator_neg_y1_m2),
        avx_curve_detail::literal_fp2_element(
            matrixnok::generated::one_m1,
            matrixnok::generated::one_m2,
            matrixnok::generated::aprime_e2_x_m1,
            matrixnok::generated::aprime_e2_x_m2)};
    return neg;
}

inline void psi(G2Point& out, const G2Point& in)
{
    out.X = avx_curve_detail::curve_fp2_ring().mul(
        avx_curve_detail::fp2_conjugate(in.X),
        avx_curve_detail::psi_frobenius_x());
    out.Y = avx_curve_detail::curve_fp2_ring().mul(
        avx_curve_detail::fp2_conjugate(in.Y),
        avx_curve_detail::psi_frobenius_y());
    out.Z = avx_curve_detail::fp2_conjugate(in.Z);
}

inline void print(const char* name, const G1Point& p)
{
    blst::scalar::print(
        name,
        avx_curve_detail::curve_g1_to_scalar({p.X, p.Y, p.Z}));
}

inline void print(const char* name, const G2Point& p)
{
    blst::scalar::print(
        name,
        avx_curve_detail::curve_g2_to_scalar({p.X, p.Y, p.Z}));
}

inline void print(const char* name, const G1AffinePoint& p)
{
    blst::scalar::print(
        name,
        avx_curve_detail::curve_g1_affine_to_scalar({p.X, p.Y}));
}

inline void print(const char* name, const G2AffinePoint& p)
{
    blst::scalar::print(
        name,
        avx_curve_detail::curve_g2_affine_to_scalar({p.X, p.Y}));
}

template<typename Group>
/// Multiply an on-curve point by `-z` using the fixed BLS12-381 chain.
/// `out` and `in` must be distinct; subgroup membership is not required.
inline void times_minus_z(Group& out, const Group& in);

template<typename Group>
inline void add_n_dbl(Group& out, const Group& p, size_t n)
{
    Group::dadd(out, out, p);
    while (n-- != 0) {
        Group::doubling(out, out);
    }
}

template<>
inline void times_minus_z<G1Point>(G1Point& out, const G1Point& in)
{
    BLST_ASSERT(&out != &in);
    G1Point::doubling(out, in);
    add_n_dbl(out, in, 2);
    add_n_dbl(out, in, 3);
    add_n_dbl(out, in, 9);
    add_n_dbl(out, in, 32);
    add_n_dbl(out, in, 16);
}

template<>
inline void times_minus_z<G2Point>(G2Point& out, const G2Point& in)
{
    BLST_ASSERT(&out != &in);
    G2Point::doubling(out, in);
    add_n_dbl(out, in, 2);
    add_n_dbl(out, in, 3);
    add_n_dbl(out, in, 9);
    add_n_dbl(out, in, 32);
    add_n_dbl(out, in, 16);
}

inline void sigma_squared(G1Point& out, const G1Point& in)
{
    out = in;
    out.X = avx_curve_detail::curve_ring().modmul(
        in.X, avx_curve_detail::curve_g1_glv().beta);
}

inline bool_t g1_in_subgroup_relation(const G1Point& point)
{
    G1Point t0 = {};
    G1Point t1 = {};

    times_minus_z(t0, point);
    times_minus_z(t1, t0);
    G1Point::cneg(t1, 1);

    sigma_squared(t0, point);

    return G1Point::equal(t0, t1);
}

inline bool_t g2_in_subgroup_relation(const G2Point& point)
{
    G2Point t0 = {};
    G2Point t1 = {};

    psi(t0, point);
    times_minus_z(t1, point);
    G2Point::cneg(t1, 1);

    return G2Point::equal(t0, t1);
}

}  // namespace blst::avx
