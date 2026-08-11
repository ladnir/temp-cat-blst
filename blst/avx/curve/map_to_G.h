#pragma once

#include "blst/avx/curve/point.h"
#include "blst/scalar/curve/map_to_G.h"

namespace blst::avx {

/*
 * TODO(big, backend-defaulting): These functions are explicit temporary
 * compatibility wrappers. They convert AVX inputs to scalar, call the scalar
 * map-to-group implementation, and convert the result back to AVX with
 * inversion-free homogeneous/Jacobian transformations.
 *
 * They are intentionally not native AVX implementations.
 *
 * Real AVX implementations are still required for:
 * - map_to_isogenous<G1Point>
 * - map_to_isogenous<G2Point>
 * - isogeny_map_to<G1Point>
 * - isogeny_map_to<G2Point>
 * - clear_cofactor<G1Point>
 * - clear_cofactor<G2Point>
 * - map_to_curve<G1Point>
 * - map_to_curve<G2Point>
 * - Encode_to_curve<G1Point>
 * - Encode_to_curve<G2Point>
 * - Hash_to_curve<G1Point>
 * - Hash_to_curve<G2Point>
 *
 * The real native path also needs missing AVX field primitives used by SSWU,
 * especially sgn0- and recip_sqrt-compatible semantics for Fp/Fp2.
 */

namespace avx_map_to_G_detail {

inline scalar::G1Point to_scalar(const G1Point& point)
{
    return avx_curve_detail::curve_g1_to_scalar({point.X, point.Y, point.Z});
}

inline scalar::G2Point to_scalar(const G2Point& point)
{
    return avx_curve_detail::curve_g2_to_scalar({point.X, point.Y, point.Z});
}

inline G1Point from_scalar(const scalar::G1Point& point)
{
    // Scalar infinity is the all-zero Jacobian value. MatrixNoK accepts any
    // Z=0 representative in addition, but its multiplication kernels require
    // the canonical (0, 1, 0) identity produced by G1Point::infinity().
    if (scalar::G1Point::is_inf(point) != 0) {
        return G1Point::infinity();
    }

    const auto x = avx_curve_detail::scalar_fp_to_ring(point.x);
    const auto y = avx_curve_detail::scalar_fp_to_ring(point.y);
    const auto z = avx_curve_detail::scalar_fp_to_ring(point.z);
    const auto z_squared = avx_curve_detail::curve_ring().modmul(z, z);

    // Jacobian (X/Z^2, Y/Z^3) to homogeneous (X/Z, Y/Z), without inversion.
    return {
        avx_curve_detail::curve_ring().modmul(x, z),
        y,
        avx_curve_detail::curve_ring().modmul(z_squared, z),
    };
}

inline G2Point from_scalar(const scalar::G2Point& point)
{
    if (scalar::G2Point::is_inf(point) != 0) {
        return G2Point::infinity();
    }

    const auto x = avx_curve_detail::scalar_fp2_to_ring(point.x);
    const auto y = avx_curve_detail::scalar_fp2_to_ring(point.y);
    const auto z = avx_curve_detail::scalar_fp2_to_ring(point.z);
    const auto z_squared = avx_curve_detail::curve_fp2_ring().mul(z, z);

    return {
        avx_curve_detail::curve_fp2_ring().mul(x, z),
        y,
        avx_curve_detail::curve_fp2_ring().mul(z_squared, z),
    };
}

inline scalar::Fp to_scalar(const G1Point::Field_type& value)
{
    return avx_curve_detail::ring_fp_to_scalar(value);
}

inline scalar::Fp2 to_scalar(const G2Point::Field_type& value)
{
    return avx_curve_detail::ring_fp2_to_scalar(value);
}

inline G1Point::Field_type from_scalar(const scalar::Fp& value)
{
    return avx_curve_detail::scalar_fp_to_ring(value);
}

inline G2Point::Field_type from_scalar(const scalar::Fp2& value)
{
    return avx_curve_detail::scalar_fp2_to_ring(value);
}

}  // namespace avx_map_to_G_detail

template<typename Group>
BLST_INLINE void isogeny_map_to(Group& out, const Group& p);

template<>
BLST_INLINE void isogeny_map_to<G1Point>(G1Point& out, const G1Point& p)
{
    scalar::G1Point scalar_out = {};
    scalar::isogeny_map_to<scalar::G1Point>(scalar_out, avx_map_to_G_detail::to_scalar(p));
    out = avx_map_to_G_detail::from_scalar(scalar_out);
}

template<>
BLST_INLINE void isogeny_map_to<G2Point>(G2Point& out, const G2Point& p)
{
    scalar::G2Point scalar_out = {};
    scalar::isogeny_map_to<scalar::G2Point>(scalar_out, avx_map_to_G_detail::to_scalar(p));
    out = avx_map_to_G_detail::from_scalar(scalar_out);
}

template<typename Group>
BLST_INLINE void map_to_isogenous(Group& p, const typename Group::Field_type& u);

template<>
BLST_INLINE void map_to_isogenous<G1Point>(G1Point& p, const G1Point::Field_type& u)
{
    scalar::G1Point scalar_p = {};
    scalar::map_to_isogenous<scalar::G1Point>(scalar_p, avx_map_to_G_detail::to_scalar(u));
    p = avx_map_to_G_detail::from_scalar(scalar_p);
}

template<>
BLST_INLINE void map_to_isogenous<G2Point>(G2Point& p, const G2Point::Field_type& u)
{
    scalar::G2Point scalar_p = {};
    scalar::map_to_isogenous<scalar::G2Point>(scalar_p, avx_map_to_G_detail::to_scalar(u));
    p = avx_map_to_G_detail::from_scalar(scalar_p);
}

template<typename Group>
BLST_INLINE void clear_cofactor(Group& out, const Group& p);

template<>
BLST_INLINE void clear_cofactor<G1Point>(G1Point& out, const G1Point& p)
{
    scalar::G1Point scalar_out = {};
    scalar::clear_cofactor<scalar::G1Point>(scalar_out, avx_map_to_G_detail::to_scalar(p));
    out = avx_map_to_G_detail::from_scalar(scalar_out);
}

template<>
BLST_INLINE void clear_cofactor<G2Point>(G2Point& out, const G2Point& p)
{
    scalar::G2Point scalar_out = {};
    scalar::clear_cofactor<scalar::G2Point>(scalar_out, avx_map_to_G_detail::to_scalar(p));
    out = avx_map_to_G_detail::from_scalar(scalar_out);
}

template<typename Group>
BLST_INLINE void map_to_curve(
    Group& p,
    const typename Group::Field_type& u,
    const typename Group::Field_type* v);

template<>
BLST_INLINE void map_to_curve<G1Point>(G1Point& p, const G1Point::Field_type& u, const G1Point::Field_type* v)
{
    scalar::G1Point scalar_p = {};
    scalar::Fp scalar_v = {};
    const scalar::Fp* scalar_v_ptr = nullptr;
    if (v != nullptr) {
        scalar_v = avx_map_to_G_detail::to_scalar(*v);
        scalar_v_ptr = &scalar_v;
    }
    scalar::map_to_curve<scalar::G1Point>(scalar_p, avx_map_to_G_detail::to_scalar(u), scalar_v_ptr);
    p = avx_map_to_G_detail::from_scalar(scalar_p);
}

template<>
BLST_INLINE void map_to_curve<G2Point>(G2Point& p, const G2Point::Field_type& u, const G2Point::Field_type* v)
{
    scalar::G2Point scalar_p = {};
    scalar::Fp2 scalar_v = {};
    const scalar::Fp2* scalar_v_ptr = nullptr;
    if (v != nullptr) {
        scalar_v = avx_map_to_G_detail::to_scalar(*v);
        scalar_v_ptr = &scalar_v;
    }
    scalar::map_to_curve<scalar::G2Point>(scalar_p, avx_map_to_G_detail::to_scalar(u), scalar_v_ptr);
    p = avx_map_to_G_detail::from_scalar(scalar_p);
}

/**
 * Encode a message to the selected G1/G2 subgroup with the RFC 9380
 * XMD:SHA-256 SSWU nonuniform map. Byte pointers must cover their stated
 * lengths and may be null only for zero length. Allocation-free and
 * nonthrowing for valid inputs; no Release input validation is performed.
 */
template<typename Group>
BLST_INLINE void Encode_to_curve(
    Group& p,
    const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len);

template<>
BLST_INLINE void Encode_to_curve<G1Point>(
    G1Point& p,
    const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len)
{
    scalar::G1Point scalar_p = {};
    scalar::Encode_to_curve<scalar::G1Point>(scalar_p, msg, msg_len, DST, DST_len, aug, aug_len);
    p = avx_map_to_G_detail::from_scalar(scalar_p);
}

template<>
BLST_INLINE void Encode_to_curve<G2Point>(
    G2Point& p,
    const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len)
{
    scalar::G2Point scalar_p = {};
    scalar::Encode_to_curve<scalar::G2Point>(scalar_p, msg, msg_len, DST, DST_len, aug, aug_len);
    p = avx_map_to_G_detail::from_scalar(scalar_p);
}

/**
 * Hash a message to the selected G1/G2 subgroup with the RFC 9380
 * XMD:SHA-256 SSWU random-oracle map. Pointer and validation contracts match
 * Encode_to_curve.
 */
template<typename Group>
BLST_INLINE void Hash_to_curve(
    Group& p,
    const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len);

template<>
BLST_INLINE void Hash_to_curve<G1Point>(
    G1Point& p,
    const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len)
{
    scalar::G1Point scalar_p = {};
    scalar::Hash_to_curve<scalar::G1Point>(scalar_p, msg, msg_len, DST, DST_len, aug, aug_len);
    p = avx_map_to_G_detail::from_scalar(scalar_p);
}

template<>
BLST_INLINE void Hash_to_curve<G2Point>(
    G2Point& p,
    const unsigned char* msg, size_t msg_len,
    const unsigned char* DST, size_t DST_len,
    const unsigned char* aug, size_t aug_len)
{
    scalar::G2Point scalar_p = {};
    scalar::Hash_to_curve<scalar::G2Point>(scalar_p, msg, msg_len, DST, DST_len, aug, aug_len);
    p = avx_map_to_G_detail::from_scalar(scalar_p);
}

}  // namespace blst::avx
