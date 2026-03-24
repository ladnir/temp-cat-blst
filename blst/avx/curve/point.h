#pragma once

#include "ec.h"
#include "scalar_mult.h"
#include "blst/avx/fields/conversion_inversion.h"
#include "blst/avx/fields/fp2.h"
#include "blst/avx/rns/bounded_ring.h"
#include "blst/scalar/curve/point.h"

#include <array>
#include <cstddef>
#include <cstring>

namespace blst::avx {

namespace avx_curve_detail {

using CurveRing = BoundedRing<381, 8, 52, -1932, 2377, 12>;
using CurveFp2Ring = FP2Ring<CurveRing>;
using CurveG1 = G1<CurveRing>;
using CurveG2 = G2<CurveRing>;
using CurveG1Projective = ProjectivePoint<typename CurveRing::StandardElement>;
using CurveG1Affine = AffinePoint<typename CurveRing::StandardElement>;
using CurveG2Projective = ProjectivePoint<typename CurveFp2Ring::StandardElement>;
using CurveG2Affine = AffinePoint<typename CurveFp2Ring::StandardElement>;

inline const BigInt& base_modulus()
{
    static const BigInt modulus(
        "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787",
        10);
    return modulus;
}

inline const CurveRing& curve_ring()
{
    static const CurveRing ring(base_modulus());
    return ring;
}

inline const CurveFp2Ring& curve_fp2_ring()
{
    static const CurveFp2Ring fp2ring(curve_ring());
    return fp2ring;
}

inline const CurveG1& curve_g1()
{
    static const CurveG1 curve;
    return curve;
}

inline const CurveG2& curve_g2()
{
    static const CurveG2 curve;
    return curve;
}

inline const GLVMult<CurveRing>& curve_g1_glv()
{
    static const GLVMult<CurveRing> mult(curve_g1(), curve_ring());
    return mult;
}

inline const GLSMult<CurveRing>& curve_g2_gls()
{
    static const GLSMult<CurveRing> mult(curve_ring());
    return mult;
}

inline bool is_zero(const CurveRing::StandardElement& value)
{
    return curve_ring().to_bigint(value) == BigInt(0);
}

inline bool is_zero(const CurveFp2Ring::StandardElement& value)
{
    return is_zero(value.x()) && is_zero(value.y());
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
    out.X = ring_fp_to_scalar(point.x);
    out.Y = ring_fp_to_scalar(point.y);
    return out;
}

inline blst::scalar::G1Point curve_g1_to_scalar(const CurveG1Projective& point)
{
    blst::scalar::G1Point out = {};
    out.X = ring_fp_to_scalar(point.x);
    out.Y = ring_fp_to_scalar(point.y);
    out.Z = ring_fp_to_scalar(point.z);
    return out;
}

inline blst::scalar::G2AffinePoint curve_g2_affine_to_scalar(const CurveG2Affine& point)
{
    blst::scalar::G2AffinePoint out = {};
    out.X = ring_fp2_to_scalar(point.x);
    out.Y = ring_fp2_to_scalar(point.y);
    return out;
}

inline blst::scalar::G2Point curve_g2_to_scalar(const CurveG2Projective& point)
{
    blst::scalar::G2Point out = {};
    out.X = ring_fp2_to_scalar(point.x);
    out.Y = ring_fp2_to_scalar(point.y);
    out.Z = ring_fp2_to_scalar(point.z);
    return out;
}

inline std::array<uint8_t, 32> scalar_to_fixed_le_bytes(const byte* scalar, size_t nbits)
{
    BLST_ASSERT(nbits <= 256);
    std::array<uint8_t, 32> out = {};
    const size_t nbytes = (nbits + 7) / 8;
    if (nbytes != 0) {
        std::memcpy(out.data(), scalar, nbytes);
        if ((nbits & 7u) != 0u) {
            out[nbytes - 1] &= static_cast<uint8_t>((1u << (nbits & 7u)) - 1u);
        }
    }
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

inline CurveRing::StandardElement curve_inverse(const CurveRing::StandardElement& value)
{
    auto mutable_value = value;
    return invert_via_blst(curve_ring(), mutable_value);
}

inline CurveFp2Ring::StandardElement curve_fp2_inverse(const CurveFp2Ring::StandardElement& value)
{
    blst::Fp2 scalar = ring_fp2_to_scalar(value);
    blst::Fp2 inverse = {};
    blst::Fp2::reciprocal(inverse, scalar);
    return scalar_fp2_to_ring(inverse);
}

inline CurveG1Affine projective_to_affine(const CurveG1Projective& point)
{
    if (is_zero(point.z)) {
        return g1_affine_infinity();
    }

    const auto z_inv = curve_inverse(point.z);
    auto [x, y] = curve_ring().batch_reduce_expand(
        curve_ring().prep_left(point.x) * z_inv,
        curve_ring().prep_left(point.y) * z_inv);
    return {x, y};
}

inline CurveG2Affine projective_to_affine(const CurveG2Projective& point)
{
    if (is_zero(point.z)) {
        return g2_affine_infinity();
    }

    const auto z_inv = curve_fp2_inverse(point.z);
    return {
        curve_fp2_ring().mul(point.x, z_inv)[0],
        curve_fp2_ring().mul(point.y, z_inv)[0],
    };
}

inline CurveG1Projective affine_to_projective(const CurveG1Affine& point)
{
    if (is_zero(point.x) && is_zero(point.y)) {
        return curve_g1().zero(curve_ring());
    }
    return {point.x, point.y, curve_ring().one()};
}

inline CurveG2Projective affine_to_projective(const CurveG2Affine& point)
{
    if (is_zero(point.x) && is_zero(point.y)) {
        return curve_g2().zero(curve_fp2_ring());
    }
    return {point.x, point.y, curve_fp2_ring().one()};
}

inline const CurveRing::StandardElement& aprime_e1()
{
    static const CurveRing::StandardElement value =
        scalar_fp_to_ring(blst::scalar::get_BLST_Aprime_E1());
    return value;
}

inline const CurveFp2Ring::StandardElement& aprime_e2()
{
    static const CurveFp2Ring::StandardElement value =
        scalar_fp2_to_ring(blst::scalar::get_BLST_Aprime_E2());
    return value;
}

inline CurveFp2Ring::StandardElement fp2_conjugate(const CurveFp2Ring::StandardElement& value)
{
    return {value.x(), curve_ring().standard_negate(value.y())};
}

inline const CurveFp2Ring::StandardElement& psi_frobenius_x()
{
    static const CurveFp2Ring::StandardElement value(
        curve_ring().zero(),
        curve_g2_gls().psi_x_i);
    return value;
}

inline const CurveFp2Ring::StandardElement& psi_frobenius_y()
{
    static const CurveFp2Ring::StandardElement value(
        curve_g2_gls().psi_y_r,
        curve_g2_gls().psi_y_i);
    return value;
}

inline const std::array<byte, 32>& minus_z_bytes()
{
    static const std::array<byte, 32> value = {0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0xd2};
    return value;
}

}  // namespace avx_curve_detail

struct G1AffinePoint
{
    using Field_type = avx_curve_detail::CurveRing::StandardElement;

    Field_type X;
    Field_type Y;

    constexpr G1AffinePoint() = default;
    constexpr G1AffinePoint(const Field_type& x, const Field_type& y) : X(x), Y(y) {}

    static const G1AffinePoint& generator();
    static int equal(const G1AffinePoint& a, const G1AffinePoint& b);
    static bool_t is_inf(const G1AffinePoint& p);
    static bool_t on_curve(const G1AffinePoint& p);
    static bool_t in_subgroup(const G1AffinePoint& p);
    static void cneg(G1AffinePoint& inOut, bool_t c);
    static limb_t serialize_BE(unsigned char* out, const G1AffinePoint& in);
    static void serialize(unsigned char* out, const G1AffinePoint& in);
    static BLST_ERROR deserialize_BE(G1AffinePoint& out, const unsigned char* in);
    static BLST_ERROR deserialize(G1AffinePoint& out, const unsigned char* in);
    static limb_t compress_BE(unsigned char* out, const G1AffinePoint& in);
    static void compress(unsigned char* out, const G1AffinePoint& in);
    static limb_t uncompress_BE(G1AffinePoint& out, const unsigned char* in);
    static BLST_ERROR uncompress(G1AffinePoint& out, const unsigned char* in);
};

struct G1Point
{
    using Field_type = avx_curve_detail::CurveRing::StandardElement;

    Field_type X;
    Field_type Y;
    Field_type Z;

    constexpr G1Point() = default;
    constexpr G1Point(const Field_type& x, const Field_type& y, const Field_type& z) : X(x), Y(y), Z(z) {}

    static void dadd(G1Point& out, const G1Point& a, const G1Point& b);
    static void dadd(G1Point& out, const G1Point& a, const G1Point& b, const Field_type* a4);
    static void dadd(G1Point& out, const G1Point& a, const G1AffinePoint& b);
    static void add(G1Point& out, const G1Point& a, const G1Point& b);
    static void add(G1Point& out, const G1Point& a, const G1AffinePoint& b);
    static void doubling(G1Point& out, const G1Point& a);
    static void mult(G1Point& out, const G1Point& a, const byte* scalar, size_t nbits);
    static void sign(G1Point& out, const G1Point& in, const pow256 SK);
    static void to_pk2(unsigned char* out, G1AffinePoint* PK, const pow256 SK);
    static void cneg(G1Point& p, limb_t cbit);
    static void neg(G1Point& out, const G1Point& p);
    static void to_affine(G1AffinePoint& out, const G1Point& in);
    static void from_affine(G1Point& out, const G1AffinePoint& a);
    static void from_jacobian(G1Point& out, const G1Point& in);
    static void cswap(G1Point& a, G1Point& b, bool_t cbit);
    static void ccopy(G1Point& a, const G1Point& b, bool_t cbit);
    static limb_t equal(const G1Point& a, const G1Point& b);
    static const G1Point& generator();
    static G1Point infinity();
    static G1Point identity();
    static bool_t is_identity(const G1Point& p);
    static bool_t is_inf(const G1Point& p);
    static bool_t on_curve(const G1Point& p);
    static bool_t in_subgroup(const G1Point& p);
    static limb_t serialize_BE(unsigned char* out, const G1Point& in);
    static void serialize(unsigned char* out, const G1Point& in);
    static limb_t compress_BE(unsigned char* out, const G1Point& in);
    static void compress(unsigned char* out, const G1Point& in);
};

struct G2AffinePoint
{
    using Field_type = avx_curve_detail::CurveFp2Ring::StandardElement;

    Field_type X;
    Field_type Y;

    constexpr G2AffinePoint() = default;
    constexpr G2AffinePoint(const Field_type& x, const Field_type& y) : X(x), Y(y) {}

    static const G2AffinePoint& generator();
    static int equal(const G2AffinePoint& a, const G2AffinePoint& b);
    static bool_t is_inf(const G2AffinePoint& p);
    static bool_t on_curve(const G2AffinePoint& p);
    static bool_t in_subgroup(const G2AffinePoint& p);
    static void cneg(G2AffinePoint& inOut, bool_t c);
    static limb_t serialize_BE(unsigned char* out, const G2AffinePoint& in);
    static void serialize(unsigned char* out, const G2AffinePoint& in);
    static BLST_ERROR deserialize_BE(G2AffinePoint& out, const unsigned char* in);
    static BLST_ERROR deserialize(G2AffinePoint& out, const unsigned char* in);
    static limb_t compress_BE(unsigned char* out, const G2AffinePoint& in);
    static void compress(unsigned char* out, const G2AffinePoint& in);
    static limb_t uncompress_BE(G2AffinePoint& out, const unsigned char* in);
    static BLST_ERROR uncompress(G2AffinePoint& out, const unsigned char* in);
};

struct G2Point
{
    using Field_type = avx_curve_detail::CurveFp2Ring::StandardElement;

    Field_type X;
    Field_type Y;
    Field_type Z;

    constexpr G2Point() = default;
    constexpr G2Point(const Field_type& x, const Field_type& y, const Field_type& z) : X(x), Y(y), Z(z) {}

    static void dadd(G2Point& out, const G2Point& a, const G2Point& b);
    static void dadd(G2Point& out, const G2Point& a, const G2Point& b, const Field_type* a4);
    static void dadd(G2Point& out, const G2Point& a, const G2AffinePoint& b);
    static void add(G2Point& out, const G2Point& a, const G2Point& b);
    static void add(G2Point& out, const G2Point& a, const G2AffinePoint& b);
    static void doubling(G2Point& out, const G2Point& a);
    static void mult(G2Point& out, const G2Point& a, const byte* scalar, size_t nbits);
    static void sign(G2Point& out, const G2Point& in, const pow256 SK);
    static void to_pk2(unsigned char* out, G2AffinePoint* PK, const pow256 SK);
    static void cneg(G2Point& p, limb_t cbit);
    static void neg(G2Point& out, const G2Point& p);
    static void to_affine(G2AffinePoint& out, const G2Point& in);
    static void from_affine(G2Point& out, const G2AffinePoint& a);
    static void from_jacobian(G2Point& out, const G2Point& in);
    static void cswap(G2Point& a, G2Point& b, bool_t cbit);
    static void ccopy(G2Point& a, const G2Point& b, bool_t cbit);
    static limb_t equal(const G2Point& a, const G2Point& b);
    static const G2Point& generator();
    static G2Point infinity();
    static G2Point identity();
    static bool_t is_identity(const G2Point& p);
    static bool_t is_inf(const G2Point& p);
    static bool_t on_curve(const G2Point& p);
    static bool_t in_subgroup(const G2Point& p);
    static limb_t serialize_BE(unsigned char* out, const G2Point& in);
    static void serialize(unsigned char* out, const G2Point& in);
    static limb_t compress_BE(unsigned char* out, const G2Point& in);
    static void compress(unsigned char* out, const G2Point& in);
};

inline const G1AffinePoint& G1AffinePoint::generator()
{
    static const G1AffinePoint g = [] {
        const blst::scalar::G1AffinePoint& scalar_g = blst::scalar::G1AffinePoint::generator();
        return G1AffinePoint(
            avx_curve_detail::scalar_fp_to_ring(scalar_g.X),
            avx_curve_detail::scalar_fp_to_ring(scalar_g.Y));
    }();
    return g;
}

inline int G1AffinePoint::equal(const G1AffinePoint& a, const G1AffinePoint& b)
{
    const blst::scalar::G1AffinePoint sa = avx_curve_detail::curve_g1_affine_to_scalar({a.X, a.Y});
    const blst::scalar::G1AffinePoint sb = avx_curve_detail::curve_g1_affine_to_scalar({b.X, b.Y});
    return blst::scalar::G1AffinePoint::equal(sa, sb);
}

inline bool_t G1AffinePoint::is_inf(const G1AffinePoint& p)
{
    return avx_curve_detail::is_zero(p.X) && avx_curve_detail::is_zero(p.Y);
}

inline bool_t G1AffinePoint::on_curve(const G1AffinePoint& p)
{
    return blst::scalar::G1AffinePoint::on_curve(avx_curve_detail::curve_g1_affine_to_scalar({p.X, p.Y}));
}

inline bool_t G1AffinePoint::in_subgroup(const G1AffinePoint& p)
{
    return blst::scalar::G1AffinePoint::in_subgroup(avx_curve_detail::curve_g1_affine_to_scalar({p.X, p.Y}));
}

inline void G1AffinePoint::cneg(G1AffinePoint& inOut, bool_t c)
{
    if (c) {
        inOut.Y = avx_curve_detail::curve_ring().standard_negate(inOut.Y);
    }
}

inline limb_t G1AffinePoint::serialize_BE(unsigned char* out, const G1AffinePoint& in)
{
    return blst::scalar::G1AffinePoint::serialize_BE(out, avx_curve_detail::curve_g1_affine_to_scalar({in.X, in.Y}));
}

inline void G1AffinePoint::serialize(unsigned char* out, const G1AffinePoint& in)
{
    blst::scalar::G1AffinePoint::serialize(out, avx_curve_detail::curve_g1_affine_to_scalar({in.X, in.Y}));
}

inline BLST_ERROR G1AffinePoint::deserialize_BE(G1AffinePoint& out, const unsigned char* in)
{
    blst::scalar::G1AffinePoint scalar = {};
    const BLST_ERROR err = blst::scalar::G1AffinePoint::deserialize_BE(scalar, in);
    if (err == BLST_SUCCESS) {
        out = G1AffinePoint(
            avx_curve_detail::scalar_fp_to_ring(scalar.X),
            avx_curve_detail::scalar_fp_to_ring(scalar.Y));
    }
    return err;
}

inline BLST_ERROR G1AffinePoint::deserialize(G1AffinePoint& out, const unsigned char* in)
{
    blst::scalar::G1AffinePoint scalar = {};
    const BLST_ERROR err = blst::scalar::G1AffinePoint::deserialize(scalar, in);
    if (err == BLST_SUCCESS) {
        out = G1AffinePoint(
            avx_curve_detail::scalar_fp_to_ring(scalar.X),
            avx_curve_detail::scalar_fp_to_ring(scalar.Y));
    }
    return err;
}

inline limb_t G1AffinePoint::compress_BE(unsigned char* out, const G1AffinePoint& in)
{
    return blst::scalar::G1AffinePoint::compress_BE(out, avx_curve_detail::curve_g1_affine_to_scalar({in.X, in.Y}));
}

inline void G1AffinePoint::compress(unsigned char* out, const G1AffinePoint& in)
{
    blst::scalar::G1AffinePoint::compress(out, avx_curve_detail::curve_g1_affine_to_scalar({in.X, in.Y}));
}

inline limb_t G1AffinePoint::uncompress_BE(G1AffinePoint& out, const unsigned char* in)
{
    blst::scalar::G1AffinePoint scalar = {};
    const limb_t parity = blst::scalar::G1AffinePoint::uncompress_BE(scalar, in);
    out = G1AffinePoint(
        avx_curve_detail::scalar_fp_to_ring(scalar.X),
        avx_curve_detail::scalar_fp_to_ring(scalar.Y));
    return parity;
}

inline BLST_ERROR G1AffinePoint::uncompress(G1AffinePoint& out, const unsigned char* in)
{
    blst::scalar::G1AffinePoint scalar = {};
    const BLST_ERROR err = blst::scalar::G1AffinePoint::uncompress(scalar, in);
    if (err == BLST_SUCCESS) {
        out = G1AffinePoint(
            avx_curve_detail::scalar_fp_to_ring(scalar.X),
            avx_curve_detail::scalar_fp_to_ring(scalar.Y));
    }
    return err;
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
    const auto scalar_bytes = avx_curve_detail::scalar_to_fixed_le_bytes(scalar, nbits);
    const auto raw = avx_curve_detail::curve_g1_glv().multiply({a.X, a.Y, a.Z}, scalar_bytes, avx_curve_detail::curve_ring());
    out = {raw.x, raw.y, raw.z};
}

inline void G1Point::sign(G1Point& out, const G1Point& in, const pow256 SK)
{
    G1AffinePoint affine = {};
    to_affine(affine, in);

    blst::scalar::G1AffinePoint scalar_affine = {};
    scalar_affine.X = avx_curve_detail::ring_fp_to_scalar(affine.X);
    scalar_affine.Y = avx_curve_detail::ring_fp_to_scalar(affine.Y);

    blst::scalar::G1Point scalar_in = {};
    blst::scalar::G1Point::from_affine(scalar_in, scalar_affine);
    blst::scalar::G1Point scalar_out = {};
    blst::scalar::G1Point::sign(scalar_out, scalar_in, SK);
    blst::scalar::G1AffinePoint scalar_out_affine = {};
    blst::scalar::G1Point::to_affine(scalar_out_affine, scalar_out);
    G1AffinePoint out_affine = {
        avx_curve_detail::scalar_fp_to_ring(scalar_out_affine.X),
        avx_curve_detail::scalar_fp_to_ring(scalar_out_affine.Y),
    };
    from_affine(out, out_affine);
}

inline void G1Point::to_pk2(unsigned char* out, G1AffinePoint* PK, const pow256 SK)
{
    blst::scalar::G1AffinePoint scalar_pk = {};
    blst::scalar::G1AffinePoint* scalar_pk_ptr = PK != nullptr ? &scalar_pk : nullptr;
    blst::scalar::G1Point::to_pk2(out, scalar_pk_ptr, SK);
    if (PK != nullptr) {
        *PK = {
            avx_curve_detail::scalar_fp_to_ring(scalar_pk.X),
            avx_curve_detail::scalar_fp_to_ring(scalar_pk.Y),
        };
    }
}

inline void G1Point::cneg(G1Point& p, limb_t cbit)
{
    if (cbit) {
        p.Y = avx_curve_detail::curve_ring().standard_negate(p.Y);
    }
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

inline void G1Point::cswap(G1Point& a, G1Point& b, bool_t cbit) { vec_cswap(&a, &b, sizeof(G1Point), cbit); }
inline void G1Point::ccopy(G1Point& a, const G1Point& b, bool_t cbit) { vec_select(&a, &b, &a, sizeof(G1Point), cbit); }

inline limb_t G1Point::equal(const G1Point& a, const G1Point& b)
{
    G1AffinePoint aa = {};
    G1AffinePoint bb = {};
    to_affine(aa, a);
    to_affine(bb, b);
    return G1AffinePoint::equal(aa, bb);
}

inline const G1Point& G1Point::generator()
{
    static const G1Point g = [] {
        G1Point out = {};
        from_affine(out, G1AffinePoint::generator());
        return out;
    }();
    return g;
}

inline G1Point G1Point::infinity()
{
    const auto raw = avx_curve_detail::curve_g1().zero(avx_curve_detail::curve_ring());
    return {raw.x, raw.y, raw.z};
}

inline G1Point G1Point::identity() { return infinity(); }
inline bool_t G1Point::is_identity(const G1Point& p) { return avx_curve_detail::is_zero(p.Z); }
inline bool_t G1Point::is_inf(const G1Point& p) { return is_identity(p); }

inline bool_t G1Point::on_curve(const G1Point& p)
{
    G1AffinePoint affine = {};
    to_affine(affine, p);
    return G1AffinePoint::on_curve(affine);
}

inline bool_t G1Point::in_subgroup(const G1Point& p)
{
    G1AffinePoint affine = {};
    to_affine(affine, p);
    return G1AffinePoint::in_subgroup(affine);
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
    static const G2AffinePoint g = [] {
        const blst::scalar::G2AffinePoint& scalar_g = blst::scalar::G2AffinePoint::generator();
        return G2AffinePoint(
            avx_curve_detail::scalar_fp2_to_ring(scalar_g.X),
            avx_curve_detail::scalar_fp2_to_ring(scalar_g.Y));
    }();
    return g;
}

inline int G2AffinePoint::equal(const G2AffinePoint& a, const G2AffinePoint& b)
{
    const blst::scalar::G2AffinePoint sa = avx_curve_detail::curve_g2_affine_to_scalar({a.X, a.Y});
    const blst::scalar::G2AffinePoint sb = avx_curve_detail::curve_g2_affine_to_scalar({b.X, b.Y});
    return blst::scalar::G2AffinePoint::equal(sa, sb);
}

inline bool_t G2AffinePoint::is_inf(const G2AffinePoint& p)
{
    return avx_curve_detail::is_zero(p.X) && avx_curve_detail::is_zero(p.Y);
}

inline bool_t G2AffinePoint::on_curve(const G2AffinePoint& p)
{
    return blst::scalar::G2AffinePoint::on_curve(avx_curve_detail::curve_g2_affine_to_scalar({p.X, p.Y}));
}

inline bool_t G2AffinePoint::in_subgroup(const G2AffinePoint& p)
{
    return blst::scalar::G2AffinePoint::in_subgroup(avx_curve_detail::curve_g2_affine_to_scalar({p.X, p.Y}));
}

inline void G2AffinePoint::cneg(G2AffinePoint& inOut, bool_t c)
{
    if (c) {
        inOut.Y = avx_curve_detail::curve_fp2_ring().standard_negate(inOut.Y);
    }
}

inline limb_t G2AffinePoint::serialize_BE(unsigned char* out, const G2AffinePoint& in)
{
    return blst::scalar::G2AffinePoint::serialize_BE(out, avx_curve_detail::curve_g2_affine_to_scalar({in.X, in.Y}));
}

inline void G2AffinePoint::serialize(unsigned char* out, const G2AffinePoint& in)
{
    blst::scalar::G2AffinePoint::serialize(out, avx_curve_detail::curve_g2_affine_to_scalar({in.X, in.Y}));
}

inline BLST_ERROR G2AffinePoint::deserialize_BE(G2AffinePoint& out, const unsigned char* in)
{
    blst::scalar::G2AffinePoint scalar = {};
    const BLST_ERROR err = blst::scalar::G2AffinePoint::deserialize_BE(scalar, in);
    if (err == BLST_SUCCESS) {
        out = G2AffinePoint(
            avx_curve_detail::scalar_fp2_to_ring(scalar.X),
            avx_curve_detail::scalar_fp2_to_ring(scalar.Y));
    }
    return err;
}

inline BLST_ERROR G2AffinePoint::deserialize(G2AffinePoint& out, const unsigned char* in)
{
    blst::scalar::G2AffinePoint scalar = {};
    const BLST_ERROR err = blst::scalar::G2AffinePoint::deserialize(scalar, in);
    if (err == BLST_SUCCESS) {
        out = G2AffinePoint(
            avx_curve_detail::scalar_fp2_to_ring(scalar.X),
            avx_curve_detail::scalar_fp2_to_ring(scalar.Y));
    }
    return err;
}

inline limb_t G2AffinePoint::compress_BE(unsigned char* out, const G2AffinePoint& in)
{
    return blst::scalar::G2AffinePoint::compress_BE(out, avx_curve_detail::curve_g2_affine_to_scalar({in.X, in.Y}));
}

inline void G2AffinePoint::compress(unsigned char* out, const G2AffinePoint& in)
{
    blst::scalar::G2AffinePoint::compress(out, avx_curve_detail::curve_g2_affine_to_scalar({in.X, in.Y}));
}

inline limb_t G2AffinePoint::uncompress_BE(G2AffinePoint& out, const unsigned char* in)
{
    blst::scalar::G2AffinePoint scalar = {};
    const limb_t parity = blst::scalar::G2AffinePoint::uncompress_BE(scalar, in);
    out = G2AffinePoint(
        avx_curve_detail::scalar_fp2_to_ring(scalar.X),
        avx_curve_detail::scalar_fp2_to_ring(scalar.Y));
    return parity;
}

inline BLST_ERROR G2AffinePoint::uncompress(G2AffinePoint& out, const unsigned char* in)
{
    blst::scalar::G2AffinePoint scalar = {};
    const BLST_ERROR err = blst::scalar::G2AffinePoint::uncompress(scalar, in);
    if (err == BLST_SUCCESS) {
        out = G2AffinePoint(
            avx_curve_detail::scalar_fp2_to_ring(scalar.X),
            avx_curve_detail::scalar_fp2_to_ring(scalar.Y));
    }
    return err;
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
    const auto scalar_bytes = avx_curve_detail::scalar_to_fixed_le_bytes(scalar, nbits);
    const auto raw = avx_curve_detail::curve_g2_gls().multiply(
        {a.X, a.Y, a.Z},
        scalar_bytes,
        avx_curve_detail::curve_fp2_ring(),
        avx_curve_detail::curve_g2(),
        avx_curve_detail::curve_ring());
    out = {raw.x, raw.y, raw.z};
}

inline void G2Point::sign(G2Point& out, const G2Point& in, const pow256 SK)
{
    G2AffinePoint affine = {};
    to_affine(affine, in);

    blst::scalar::G2AffinePoint scalar_affine = {};
    scalar_affine.X = avx_curve_detail::ring_fp2_to_scalar(affine.X);
    scalar_affine.Y = avx_curve_detail::ring_fp2_to_scalar(affine.Y);

    blst::scalar::G2Point scalar_in = {};
    blst::scalar::G2Point::from_affine(scalar_in, scalar_affine);
    blst::scalar::G2Point scalar_out = {};
    blst::scalar::G2Point::sign(scalar_out, scalar_in, SK);
    blst::scalar::G2AffinePoint scalar_out_affine = {};
    blst::scalar::G2Point::to_affine(scalar_out_affine, scalar_out);
    G2AffinePoint out_affine = {
        avx_curve_detail::scalar_fp2_to_ring(scalar_out_affine.X),
        avx_curve_detail::scalar_fp2_to_ring(scalar_out_affine.Y),
    };
    from_affine(out, out_affine);
}

inline void G2Point::to_pk2(unsigned char* out, G2AffinePoint* PK, const pow256 SK)
{
    blst::scalar::G2AffinePoint scalar_pk = {};
    blst::scalar::G2AffinePoint* scalar_pk_ptr = PK != nullptr ? &scalar_pk : nullptr;
    blst::scalar::G2Point::to_pk2(out, scalar_pk_ptr, SK);
    if (PK != nullptr) {
        *PK = {
            avx_curve_detail::scalar_fp2_to_ring(scalar_pk.X),
            avx_curve_detail::scalar_fp2_to_ring(scalar_pk.Y),
        };
    }
}

inline void G2Point::cneg(G2Point& p, limb_t cbit)
{
    if (cbit) {
        p.Y = avx_curve_detail::curve_fp2_ring().standard_negate(p.Y);
    }
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

inline void G2Point::cswap(G2Point& a, G2Point& b, bool_t cbit) { vec_cswap(&a, &b, sizeof(G2Point), cbit); }
inline void G2Point::ccopy(G2Point& a, const G2Point& b, bool_t cbit) { vec_select(&a, &b, &a, sizeof(G2Point), cbit); }

inline limb_t G2Point::equal(const G2Point& a, const G2Point& b)
{
    G2AffinePoint aa = {};
    G2AffinePoint bb = {};
    to_affine(aa, a);
    to_affine(bb, b);
    return G2AffinePoint::equal(aa, bb);
}

inline const G2Point& G2Point::generator()
{
    static const G2Point g = [] {
        G2Point out = {};
        from_affine(out, G2AffinePoint::generator());
        return out;
    }();
    return g;
}

inline G2Point G2Point::infinity()
{
    const auto raw = avx_curve_detail::curve_g2().zero(avx_curve_detail::curve_fp2_ring());
    return {raw.x, raw.y, raw.z};
}

inline G2Point G2Point::identity() { return infinity(); }
inline bool_t G2Point::is_identity(const G2Point& p) { return avx_curve_detail::is_zero(p.Z); }
inline bool_t G2Point::is_inf(const G2Point& p) { return is_identity(p); }

inline bool_t G2Point::on_curve(const G2Point& p)
{
    G2AffinePoint affine = {};
    to_affine(affine, p);
    return G2AffinePoint::on_curve(affine);
}

inline bool_t G2Point::in_subgroup(const G2Point& p)
{
    G2AffinePoint affine = {};
    to_affine(affine, p);
    return G2AffinePoint::in_subgroup(affine);
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
    static const G1Point neg = [] {
        G1Point out = G1Point::generator();
        G1Point::cneg(out, 1);
        return out;
    }();
    return neg;
}

inline const G2Point& get_BLST_NEG_G2()
{
    static const G2Point neg = [] {
        G2Point out = G2Point::generator();
        G2Point::cneg(out, 1);
        return out;
    }();
    return neg;
}

inline void psi(G2Point& out, const G2Point& in)
{
    out.X = avx_curve_detail::curve_fp2_ring().mul(
        avx_curve_detail::fp2_conjugate(in.X),
        avx_curve_detail::psi_frobenius_x())[0];
    out.Y = avx_curve_detail::curve_fp2_ring().mul(
        avx_curve_detail::fp2_conjugate(in.Y),
        avx_curve_detail::psi_frobenius_y())[0];
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
inline void times_minus_z(Group& out, const Group& in);

template<>
inline void times_minus_z<G1Point>(G1Point& out, const G1Point& in)
{
    G1Point::mult(out, in, avx_curve_detail::minus_z_bytes().data(), 64);
}

template<>
inline void times_minus_z<G2Point>(G2Point& out, const G2Point& in)
{
    G2Point::mult(out, in, avx_curve_detail::minus_z_bytes().data(), 64);
}

}  // namespace blst::avx
