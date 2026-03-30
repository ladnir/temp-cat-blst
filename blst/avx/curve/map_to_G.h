#pragma once

#include "blst/avx/curve/point.h"
#include "blst/scalar/curve/map_to_G.h"

namespace blst::avx {

/*
 * TODO(big, backend-defaulting): These functions are explicit temporary
 * compatibility wrappers. They convert AVX inputs to scalar, call the scalar
 * map-to-group implementation, and convert the result back to AVX.
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
    G1AffinePoint affine = {};
    G1Point::to_affine(affine, point);

    scalar::G1AffinePoint scalar_affine = {};
    scalar_affine.X = avx_curve_detail::ring_fp_to_scalar(affine.X);
    scalar_affine.Y = avx_curve_detail::ring_fp_to_scalar(affine.Y);

    scalar::G1Point out = {};
    scalar::G1Point::from_affine(out, scalar_affine);
    return out;
}

inline scalar::G2Point to_scalar(const G2Point& point)
{
    G2AffinePoint affine = {};
    G2Point::to_affine(affine, point);

    scalar::G2AffinePoint scalar_affine = {};
    scalar_affine.X = avx_curve_detail::ring_fp2_to_scalar(affine.X);
    scalar_affine.Y = avx_curve_detail::ring_fp2_to_scalar(affine.Y);

    scalar::G2Point out = {};
    scalar::G2Point::from_affine(out, scalar_affine);
    return out;
}

inline G1Point from_scalar(const scalar::G1Point& point)
{
    scalar::G1AffinePoint scalar_affine = {};
    scalar::G1Point::to_affine(scalar_affine, point);

    G1AffinePoint affine = {
        avx_curve_detail::scalar_fp_to_ring(scalar_affine.X),
        avx_curve_detail::scalar_fp_to_ring(scalar_affine.Y),
    };

    G1Point out = {};
    G1Point::from_affine(out, affine);
    return out;
}

inline G2Point from_scalar(const scalar::G2Point& point)
{
    scalar::G2AffinePoint scalar_affine = {};
    scalar::G2Point::to_affine(scalar_affine, point);

    G2AffinePoint affine = {
        avx_curve_detail::scalar_fp2_to_ring(scalar_affine.X),
        avx_curve_detail::scalar_fp2_to_ring(scalar_affine.Y),
    };

    G2Point out = {};
    G2Point::from_affine(out, affine);
    return out;
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
