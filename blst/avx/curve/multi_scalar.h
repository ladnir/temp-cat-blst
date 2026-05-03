#pragma once

#include "point.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

namespace blst::avx {

namespace avx_msm_detail {

inline size_t pippenger_window_size_fallback(size_t npoints)
{
    size_t wbits = 0;
    while (npoints >>= 1) {
        ++wbits;
    }

    if (wbits > 12) {
        return wbits - 3;
    }
    if (wbits > 8) {
        return wbits - 2;
    }
    if (wbits > 4) {
        return wbits - 1;
    }
    return wbits != 0 ? 2 : 1;
}

inline size_t g1_pippenger_window_size(size_t npoints)
{
    if (npoints <= 1) {
        return 1;
    }
    if (npoints <= 128) {
        return pippenger_window_size_fallback(npoints);
    }
    if (npoints <= 256) {
        return 6;
    }
    if (npoints <= 512) {
        return 7;
    }
    if (npoints <= 1024) {
        return 8;
    }
    if (npoints <= 4096) {
        return 9;
    }
    if (npoints <= 16384) {
        return 11;
    }
    if (npoints <= 32768) {
        return 12;
    }
    return pippenger_window_size_fallback(npoints);
}

inline size_t g2_pippenger_window_size(size_t npoints)
{
    if (npoints <= 1) {
        return 1;
    }
    if (npoints <= 128) {
        return pippenger_window_size_fallback(npoints);
    }
    if (npoints <= 256) {
        return 6;
    }
    if (npoints <= 512) {
        return 7;
    }
    if (npoints <= 2048) {
        return 8;
    }
    if (npoints <= 4096) {
        return 10;
    }
    if (npoints <= 8192) {
        return 11;
    }
    if (npoints <= 16384) {
        return 11;
    }
    if (npoints <= 32768) {
        return 12;
    }
    return pippenger_window_size_fallback(npoints);
}

inline size_t g1_glv_pippenger_window_size(size_t npoints)
{
    return pippenger_window_size_fallback(2 * npoints);
}

inline size_t g2_gls_pippenger_window_size(size_t npoints)
{
    if (npoints <= 1) {
        return 1;
    }
    if (npoints <= 128) {
        return pippenger_window_size_fallback(4 * npoints);
    }
    if (npoints <= 1024) {
        return 10;
    }
    if (npoints <= 4096) {
        return 11;
    }
    if (npoints <= 32768) {
        return 13;
    }
    return pippenger_window_size_fallback(4 * npoints);
}

inline size_t pippenger_window_size(size_t npoints)
{
    return std::max(g1_pippenger_window_size(npoints), g2_pippenger_window_size(npoints));
}

inline limb_t get_wval_limb(const byte* scalar, size_t bit0, size_t nbits)
{
    const size_t byte_index = bit0 >> 3;
    const size_t bit_shift = bit0 & 7u;
    const size_t bytes = (bit_shift + nbits + 7u) >> 3;
    limb_t out = scalar[byte_index];
    if (bytes > 1) {
        out |= static_cast<limb_t>(scalar[byte_index + 1]) << 8;
    }
    if (bytes > 2) {
        out |= static_cast<limb_t>(scalar[byte_index + 2]) << 16;
    }
    return out >> bit_shift;
}

inline limb_t booth_encode(limb_t wval, size_t wbits)
{
    const limb_t sign = wval >> wbits;
    limb_t mask = 0 - sign;
    BLST_launder(mask);

    const limb_t val = (wval + 1) >> 1;
    return (val ^ mask) - mask;
}

struct WindowExtract
{
    size_t missing_overlap;
    size_t extract_begin;
    size_t extract_size;
    limb_t extracted_mask;
    size_t encode_size;

    WindowExtract(size_t window_idx, size_t window_size, size_t nbits)
    {
        const size_t window_begin = window_idx * window_size;
        missing_overlap = window_begin == 0 ? 1 : 0;
        extract_begin = window_begin > 0 ? window_begin - 1 : 0;
        const size_t window_end = std::min(window_begin + window_size, nbits);
        extract_size = window_end - extract_begin;
        extracted_mask = (static_cast<limb_t>(1) << extract_size) - 1;
        encode_size = std::min(nbits - window_begin + 1, window_size);
    }

    limb_t extract(const byte* scalar) const
    {
        const limb_t wval = (get_wval_limb(scalar, extract_begin, extract_size) & extracted_mask)
            << missing_overlap;
        return booth_encode(wval, encode_size);
    }
};

inline G1AffinePoint maybe_negate(const G1AffinePoint& point, bool negate)
{
    if (!negate) {
        return point;
    }
    return {point.X, avx_curve_detail::curve_ring().standard_negate(point.Y)};
}

inline G2AffinePoint maybe_negate(const G2AffinePoint& point, bool negate)
{
    if (!negate) {
        return point;
    }
    return {point.X, avx_curve_detail::curve_fp2_ring().standard_negate(point.Y)};
}

using G1PippengerBucket = G1Point;
using G2PippengerBucket = G2Point;

template<class BucketPoint>
inline BucketPoint bucket_zero()
{
    return BucketPoint::zero();
}

inline G2Point bucket_zero_g2()
{
    return G2Point::infinity();
}

template<>
inline G1Point bucket_zero<G1Point>()
{
    return G1Point::infinity();
}

template<>
inline G2Point bucket_zero<G2Point>()
{
    return G2Point::infinity();
}

template<class BucketPoint, class AffinePoint>
inline void bucket_add_affine(BucketPoint& out, const BucketPoint& a, const AffinePoint& point, bool subtract)
{
    BucketPoint::dadd(out, a, point, subtract);
}

inline void bucket_add_affine(G2Point& out, const G2Point& a, const G2AffinePoint& point, bool subtract)
{
    if (subtract) {
        const G2AffinePoint signed_point = maybe_negate(point, true);
        G2Point::add(out, a, signed_point);
    } else {
        G2Point::add(out, a, point);
    }
}

inline void bucket_add_affine(G1Point& out, const G1Point& a, const G1AffinePoint& point, bool subtract)
{
    if (subtract) {
        const G1AffinePoint signed_point = maybe_negate(point, true);
        G1Point::add(out, a, signed_point);
    } else {
        G1Point::add(out, a, point);
    }
}

inline void bucket_set_affine(G1Point& out, const G1AffinePoint& point, bool subtract)
{
    if (subtract) {
        const G1AffinePoint signed_point = maybe_negate(point, true);
        G1Point::from_affine(out, signed_point);
    } else {
        G1Point::from_affine(out, point);
    }
}

inline void bucket_set_affine(G2Point& out, const G2AffinePoint& point, bool subtract)
{
    if (subtract) {
        const G2AffinePoint signed_point = maybe_negate(point, true);
        G2Point::from_affine(out, signed_point);
    } else {
        G2Point::from_affine(out, point);
    }
}

template<class BucketPoint>
inline void bucket_add_bucket(BucketPoint& out, const BucketPoint& a, const BucketPoint& b)
{
    BucketPoint::dadd(out, a, b);
}

inline void bucket_add_bucket(G2Point& out, const G2Point& a, const G2Point& b)
{
    G2Point::add(out, a, b);
}

inline void bucket_add_bucket(G1Point& out, const G1Point& a, const G1Point& b)
{
    G1Point::add(out, a, b);
}

inline void bucket_to_jacobian(G2Point& out, const G2Point& in)
{
    out = in;
}

inline void bucket_to_jacobian(G1Point& out, const G1Point& in)
{
    out = in;
}

template<class BucketPoint, class AffinePoint>
inline void point_bucket(
    BucketPoint* buckets,
    limb_t booth_idx,
    size_t wbits,
    const AffinePoint& point)
{
    const bool booth_sign = ((booth_idx >> wbits) & 1u) != 0;
    booth_idx &= (static_cast<limb_t>(1) << wbits) - 1;
    if (booth_idx == 0) {
        return;
    }

    bucket_add_affine(buckets[booth_idx - 1], buckets[booth_idx - 1], point, booth_sign);
}

template<class BucketPoint, class AffinePoint>
inline void point_bucket(
    BucketPoint* buckets,
    uint8_t* bucket_used,
    uint8_t bucket_tag,
    limb_t booth_idx,
    size_t wbits,
    const AffinePoint& point)
{
    const bool booth_sign = ((booth_idx >> wbits) & 1u) != 0;
    booth_idx &= (static_cast<limb_t>(1) << wbits) - 1;
    if (booth_idx == 0) {
        return;
    }

    const size_t bucket_idx = static_cast<size_t>(booth_idx - 1);
    if (bucket_used[bucket_idx] != bucket_tag) {
        bucket_set_affine(buckets[bucket_idx], point, booth_sign);
        bucket_used[bucket_idx] = bucket_tag;
    } else {
        bucket_add_affine(buckets[bucket_idx], buckets[bucket_idx], point, booth_sign);
    }
}

template<class BucketPoint>
inline void point_prefetch(const BucketPoint* buckets, limb_t booth_idx, size_t wbits)
{
    booth_idx &= (static_cast<limb_t>(1) << wbits) - 1;
    if (booth_idx == 0) {
        return;
    }
#if defined(__GNUC__) || defined(__clang__)
    __builtin_prefetch(&buckets[booth_idx - 1], 1, 1);
#elif defined(_MSC_VER)
    _mm_prefetch(reinterpret_cast<const char*>(&buckets[booth_idx - 1]), _MM_HINT_T0);
#endif
}

template<class ProjectivePoint, class BucketPoint>
inline void point_integrate_buckets(ProjectivePoint& out, BucketPoint* buckets, size_t wbits)
{
    size_t n = static_cast<size_t>(1) << wbits;
    BucketPoint acc = buckets[--n];
    BucketPoint ret = buckets[n];
    buckets[n] = bucket_zero<BucketPoint>();

    while (n--) {
        bucket_add_bucket(acc, acc, buckets[n]);
        bucket_add_bucket(ret, ret, acc);
        buckets[n] = bucket_zero<BucketPoint>();
    }
    bucket_to_jacobian(out, ret);
}

template<class ProjectivePoint, class BucketPoint>
inline void point_integrate_buckets(ProjectivePoint& out, BucketPoint* buckets, const uint8_t* bucket_used, uint8_t bucket_tag, size_t wbits)
{
    size_t n = static_cast<size_t>(1) << wbits;
    BucketPoint acc = bucket_zero<BucketPoint>();
    BucketPoint ret = bucket_zero<BucketPoint>();
    bool acc_nonzero = false;
    bool ret_nonzero = false;

    while (n--) {
        if (bucket_used[n] == bucket_tag) {
            if (acc_nonzero) {
                bucket_add_bucket(acc, acc, buckets[n]);
            } else {
                acc = buckets[n];
                acc_nonzero = true;
            }
        }
        if (acc_nonzero) {
            if (ret_nonzero) {
                bucket_add_bucket(ret, ret, acc);
            } else {
                ret = acc;
                ret_nonzero = true;
            }
        }
    }
    bucket_to_jacobian(out, ret);
}

template<class ProjectivePoint, class AffinePoint, class BucketPoint>
void points_tile_pippenger(
    ProjectivePoint& ret,
    const AffinePoint* points,
    const byte* scalars,
    size_t npoints,
    size_t nbits,
    BucketPoint* buckets,
    uint8_t* bucket_used,
    uint8_t bucket_tag,
    size_t window_idx,
    size_t window_size)
{
    const WindowExtract window(window_idx, window_size, nbits);
    const size_t bucket_count = static_cast<size_t>(1) << (window.encode_size - 1);
    if (bucket_used == nullptr) {
        std::fill_n(buckets, bucket_count, bucket_zero<BucketPoint>());
    }

    const size_t nbytes = (nbits + 7) / 8;
    limb_t encoding = window.extract(scalars);
    for (size_t i = 0; i + 1 < npoints; ++i) {
        const limb_t next_encoding = window.extract(scalars + (i + 1) * nbytes);
        point_prefetch(buckets, next_encoding, window.encode_size);
        if (bucket_used != nullptr) {
            point_bucket(buckets, bucket_used, bucket_tag, encoding, window.encode_size, points[i]);
        } else {
            point_bucket(buckets, encoding, window.encode_size, points[i]);
        }
        encoding = next_encoding;
    }
    if (bucket_used != nullptr) {
        point_bucket(buckets, bucket_used, bucket_tag, encoding, window.encode_size, points[npoints - 1]);
    } else {
        point_bucket(buckets, encoding, window.encode_size, points[npoints - 1]);
    }

    if (bucket_used != nullptr) {
        point_integrate_buckets(ret, buckets, bucket_used, bucket_tag, window.encode_size - 1);
    } else {
        point_integrate_buckets(ret, buckets, window.encode_size - 1);
    }
}

template<class ProjectivePoint, class AffinePoint, class BucketPoint>
void mult_pippenger(
    ProjectivePoint& ret,
    const AffinePoint* points,
    const byte* scalars,
    size_t npoints,
    size_t nbits,
    BucketPoint* buckets,
    uint8_t* bucket_used,
    size_t window_size)
{
    ret = ProjectivePoint::infinity();
    if (npoints == 0 || points == nullptr || scalars == nullptr || nbits == 0) {
        return;
    }

    window_size = window_size != 0 ? window_size : pippenger_window_size(npoints);
    const size_t num_windows = 1 + nbits / window_size;
    ProjectivePoint tile = ProjectivePoint::infinity();
    uint8_t bucket_tag = 1;
    if (bucket_used != nullptr) {
        std::fill_n(bucket_used, static_cast<size_t>(1) << (window_size - 1), 0);
    }

    for (size_t window_index = num_windows; window_index-- > 0;) {
        points_tile_pippenger(
            tile, points, scalars, npoints, nbits, buckets, bucket_used, bucket_tag, window_index, window_size);
        ProjectivePoint::add(ret, ret, tile);
        if (window_index != 0) {
            for (size_t i = 0; i < window_size; ++i) {
                ProjectivePoint::doubling(ret, ret);
            }
        }
        if (bucket_used != nullptr && ++bucket_tag == 0) {
            std::fill_n(bucket_used, static_cast<size_t>(1) << (window_size - 1), 0);
            bucket_tag = 1;
        }
    }
}

template<class ProjectivePoint, class AffinePoint, class BucketPoint>
void mult_pippenger_alloc(
    ProjectivePoint& ret,
    const AffinePoint* points,
    const byte* scalars,
    size_t npoints,
    size_t nbits,
    size_t window_size)
{
    window_size = window_size != 0 ? window_size : pippenger_window_size(npoints);
    std::vector<BucketPoint> buckets(static_cast<size_t>(1) << (window_size - 1));
    std::vector<uint8_t> bucket_used(static_cast<size_t>(1) << (window_size - 1));
    mult_pippenger(ret, points, scalars, npoints, nbits, buckets.data(), bucket_used.data(), window_size);
}

}  // namespace avx_msm_detail

inline size_t avx_pippenger_window_size(size_t npoints)
{
    return avx_msm_detail::pippenger_window_size(npoints);
}

inline size_t avx_g1_pippenger_window_size(size_t npoints)
{
    return avx_msm_detail::g1_pippenger_window_size(npoints);
}

inline size_t avx_g2_pippenger_window_size(size_t npoints)
{
    return avx_msm_detail::g2_pippenger_window_size(npoints);
}

inline size_t avx_g1_glv_pippenger_window_size(size_t npoints)
{
    return avx_msm_detail::g1_glv_pippenger_window_size(npoints);
}

inline size_t avx_g2_gls_pippenger_window_size(size_t npoints)
{
    return avx_msm_detail::g2_gls_pippenger_window_size(npoints);
}

inline size_t avx_g1_mult_pippenger_scratch_sizeof(size_t npoints)
{
    return sizeof(avx_msm_detail::G1PippengerBucket) << (avx_g1_pippenger_window_size(npoints) - 1);
}

inline size_t avx_g1_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints)
{
    return static_cast<size_t>(1) << (avx_g1_pippenger_window_size(npoints) - 1);
}

inline size_t avx_g1_mult_pippenger_glv_scratch_sizeof(size_t npoints)
{
    return sizeof(avx_msm_detail::G1PippengerBucket) << (avx_g1_glv_pippenger_window_size(npoints) - 1);
}

inline size_t avx_g1_mult_pippenger_glv_bucket_used_scratch_sizeof(size_t npoints)
{
    return static_cast<size_t>(1) << (avx_g1_glv_pippenger_window_size(npoints) - 1);
}

inline size_t avx_g1_mult_pippenger_glv_points_scratch_sizeof(size_t npoints)
{
    return sizeof(G1AffinePoint) * 2 * npoints;
}

inline size_t avx_g1_mult_pippenger_glv_scalars_scratch_sizeof(size_t npoints)
{
    return 32 * npoints;
}

inline size_t avx_g2_mult_pippenger_scratch_sizeof(size_t npoints)
{
    return sizeof(avx_msm_detail::G2PippengerBucket) << (avx_g2_pippenger_window_size(npoints) - 1);
}

inline size_t avx_g2_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints)
{
    return static_cast<size_t>(1) << (avx_g2_pippenger_window_size(npoints) - 1);
}

inline size_t avx_g2_mult_pippenger_gls_scratch_sizeof(size_t npoints)
{
    return sizeof(avx_msm_detail::G2PippengerBucket) << (avx_g2_gls_pippenger_window_size(npoints) - 1);
}

inline size_t avx_g2_mult_pippenger_gls_bucket_used_scratch_sizeof(size_t npoints)
{
    return static_cast<size_t>(1) << (avx_g2_gls_pippenger_window_size(npoints) - 1);
}

inline size_t avx_g2_mult_pippenger_gls_points_scratch_sizeof(size_t npoints)
{
    return sizeof(G2AffinePoint) * 4 * npoints;
}

inline size_t avx_g2_mult_pippenger_gls_scalars_scratch_sizeof(size_t npoints)
{
    return 32 * npoints;
}

using AvxG1PippengerBucket = avx_msm_detail::G1PippengerBucket;
using AvxG2PippengerBucket = avx_msm_detail::G2PippengerBucket;

inline void avx_g1_mult_pippenger(
    G1Point& ret,
    const G1AffinePoint* points,
    const byte* scalars,
    size_t npoints,
    size_t nbits,
    AvxG1PippengerBucket* buckets = nullptr,
    size_t window_size = 0,
    uint8_t* bucket_used = nullptr)
{
    window_size = window_size != 0 ? window_size : avx_g1_pippenger_window_size(npoints);
    if (buckets != nullptr) {
        if (bucket_used != nullptr) {
            avx_msm_detail::mult_pippenger(
                ret, points, scalars, npoints, nbits, buckets, bucket_used, window_size);
        } else {
            std::vector<uint8_t> local_bucket_used(static_cast<size_t>(1) << (window_size - 1));
            avx_msm_detail::mult_pippenger(
                ret, points, scalars, npoints, nbits, buckets, local_bucket_used.data(), window_size);
        }
    } else {
        avx_msm_detail::mult_pippenger_alloc<G1Point, G1AffinePoint, AvxG1PippengerBucket>(
            ret, points, scalars, npoints, nbits, window_size);
    }
}

inline void avx_g1_mult_pippenger_glv_prepare(
    G1AffinePoint* split_points,
    byte* split_scalars,
    const G1AffinePoint* points,
    const byte* scalars,
    size_t npoints)
{
    if (npoints == 0 || split_points == nullptr || split_scalars == nullptr || points == nullptr || scalars == nullptr) {
        return;
    }

    const auto& ring = avx_curve_detail::curve_ring();
    const auto& glv = avx_curve_detail::curve_g1_glv();
    constexpr size_t batch = 16;
    std::array<avx_curve_detail::CurveRing::StandardElement, batch> xs;
    std::array<avx_curve_detail::CurveRing::StandardElement, batch> betas;
    betas.fill(glv.beta);

    size_t i = 0;
    for (; i + batch <= npoints; i += batch) {
        for (size_t j = 0; j < batch; ++j) {
            const G1AffinePoint& point = points[i + j];
            split_points[2 * (i + j)] = point;
            xs[j] = point.X;
            split_points[2 * (i + j) + 1].Y = ring.standard_negate(point.Y);

            std::array<uint8_t, 32> scalar = {};
            std::memcpy(scalar.data(), scalars + 32 * (i + j), 32);
            auto parts = div_by_zz(scalar);
            std::memcpy(split_scalars + (2 * (i + j)) * 16, parts.first.data(), 16);
            std::memcpy(split_scalars + (2 * (i + j) + 1) * 16, parts.second.data(), 16);
        }
        const auto phis = ring.template batch_modmul<batch>(xs, betas);
        for (size_t j = 0; j < batch; ++j) {
            split_points[2 * (i + j) + 1].X = phis[j];
        }
    }
    for (; i < npoints; ++i) {
        split_points[2 * i] = points[i];
        split_points[2 * i + 1].X = ring.modmul(points[i].X, glv.beta);
        split_points[2 * i + 1].Y = ring.standard_negate(points[i].Y);

        std::array<uint8_t, 32> scalar = {};
        std::memcpy(scalar.data(), scalars + 32 * i, 32);
        auto parts = div_by_zz(scalar);
        std::memcpy(split_scalars + (2 * i) * 16, parts.first.data(), 16);
        std::memcpy(split_scalars + (2 * i + 1) * 16, parts.second.data(), 16);
    }
}

inline void avx_g1_mult_pippenger_glv_prepared(
    G1Point& ret,
    const G1AffinePoint* split_points,
    const byte* split_scalars,
    size_t npoints,
    AvxG1PippengerBucket* buckets = nullptr,
    size_t window_size = 0,
    uint8_t* bucket_used = nullptr)
{
    window_size = window_size != 0 ? window_size : avx_g1_glv_pippenger_window_size(npoints);
    if (buckets != nullptr) {
        if (bucket_used != nullptr) {
            avx_msm_detail::mult_pippenger(
                ret, split_points, split_scalars, 2 * npoints, 128, buckets, bucket_used, window_size);
        } else {
            std::vector<uint8_t> local_bucket_used(static_cast<size_t>(1) << (window_size - 1));
            avx_msm_detail::mult_pippenger(
                ret, split_points, split_scalars, 2 * npoints, 128, buckets, local_bucket_used.data(), window_size);
        }
    } else {
        avx_msm_detail::mult_pippenger_alloc<G1Point, G1AffinePoint, AvxG1PippengerBucket>(
            ret, split_points, split_scalars, 2 * npoints, 128, window_size);
    }
}

inline void avx_g1_mult_pippenger_glv(
    G1Point& ret,
    const G1AffinePoint* points,
    const byte* scalars,
    size_t npoints,
    AvxG1PippengerBucket* buckets = nullptr,
    size_t window_size = 0,
    uint8_t* bucket_used = nullptr)
{
    ret = G1Point::infinity();
    if (npoints == 0 || points == nullptr || scalars == nullptr) {
        return;
    }

    std::vector<G1AffinePoint> split_points(2 * npoints);
    std::vector<byte> split_scalars(2 * npoints * 16);
    avx_g1_mult_pippenger_glv_prepare(split_points.data(), split_scalars.data(), points, scalars, npoints);

    avx_g1_mult_pippenger_glv_prepared(
        ret, split_points.data(), split_scalars.data(), npoints, buckets, window_size, bucket_used);
}

inline G2AffinePoint avx_g2_minus_psi_affine(const G2AffinePoint& point)
{
    const auto& ring = avx_curve_detail::curve_ring();
    const auto& gls = avx_curve_detail::curve_g2_gls();
    auto x_r = point.X.y() * gls.psi_x_i;
    auto x_i = point.X.x() * gls.psi_x_i;
    auto common = point.Y.y() * gls.psi_y_r;
    auto y_r = common + point.Y.x() * gls.psi_y_i;
    auto y_i = common + point.Y.x() * gls.psi_y_r;
    auto [x_r_reduced, x_i_reduced, y_r_reduced, y_i_reduced] =
        ring.batch_reduce_expand(x_r, x_i, y_r, y_i);
    return {
        avx_curve_detail::CurveFp2Ring::StandardElement(x_r_reduced, x_i_reduced),
        avx_curve_detail::CurveFp2Ring::StandardElement(y_r_reduced, y_i_reduced)};
}

inline void avx_g2_minus_psi_affine_batch4(G2AffinePoint* out, const G2AffinePoint* point)
{
    const auto& ring = avx_curve_detail::curve_ring();
    const auto& gls = avx_curve_detail::curve_g2_gls();

    auto x0_r = point[0].X.y() * gls.psi_x_i;
    auto x0_i = point[0].X.x() * gls.psi_x_i;
    auto c0 = point[0].Y.y() * gls.psi_y_r;
    auto y0_r = c0 + point[0].Y.x() * gls.psi_y_i;
    auto y0_i = c0 + point[0].Y.x() * gls.psi_y_r;

    auto x1_r = point[1].X.y() * gls.psi_x_i;
    auto x1_i = point[1].X.x() * gls.psi_x_i;
    auto c1 = point[1].Y.y() * gls.psi_y_r;
    auto y1_r = c1 + point[1].Y.x() * gls.psi_y_i;
    auto y1_i = c1 + point[1].Y.x() * gls.psi_y_r;

    auto x2_r = point[2].X.y() * gls.psi_x_i;
    auto x2_i = point[2].X.x() * gls.psi_x_i;
    auto c2 = point[2].Y.y() * gls.psi_y_r;
    auto y2_r = c2 + point[2].Y.x() * gls.psi_y_i;
    auto y2_i = c2 + point[2].Y.x() * gls.psi_y_r;

    auto x3_r = point[3].X.y() * gls.psi_x_i;
    auto x3_i = point[3].X.x() * gls.psi_x_i;
    auto c3 = point[3].Y.y() * gls.psi_y_r;
    auto y3_r = c3 + point[3].Y.x() * gls.psi_y_i;
    auto y3_i = c3 + point[3].Y.x() * gls.psi_y_r;

    auto [
        x0_r_reduced, x0_i_reduced, y0_r_reduced, y0_i_reduced,
        x1_r_reduced, x1_i_reduced, y1_r_reduced, y1_i_reduced,
        x2_r_reduced, x2_i_reduced, y2_r_reduced, y2_i_reduced,
        x3_r_reduced, x3_i_reduced, y3_r_reduced, y3_i_reduced] =
        ring.batch_reduce_expand(
            x0_r, x0_i, y0_r, y0_i,
            x1_r, x1_i, y1_r, y1_i,
            x2_r, x2_i, y2_r, y2_i,
            x3_r, x3_i, y3_r, y3_i);

    out[0] = {
        avx_curve_detail::CurveFp2Ring::StandardElement(x0_r_reduced, x0_i_reduced),
        avx_curve_detail::CurveFp2Ring::StandardElement(y0_r_reduced, y0_i_reduced)};
    out[1] = {
        avx_curve_detail::CurveFp2Ring::StandardElement(x1_r_reduced, x1_i_reduced),
        avx_curve_detail::CurveFp2Ring::StandardElement(y1_r_reduced, y1_i_reduced)};
    out[2] = {
        avx_curve_detail::CurveFp2Ring::StandardElement(x2_r_reduced, x2_i_reduced),
        avx_curve_detail::CurveFp2Ring::StandardElement(y2_r_reduced, y2_i_reduced)};
    out[3] = {
        avx_curve_detail::CurveFp2Ring::StandardElement(x3_r_reduced, x3_i_reduced),
        avx_curve_detail::CurveFp2Ring::StandardElement(y3_r_reduced, y3_i_reduced)};
}

inline void avx_g2_mult_pippenger_gls_prepare(
    G2AffinePoint* split_points,
    byte* split_scalars,
    const G2AffinePoint* points,
    const byte* scalars,
    size_t npoints)
{
    if (npoints == 0 || split_points == nullptr || split_scalars == nullptr || points == nullptr || scalars == nullptr) {
        return;
    }

    constexpr size_t batch = 4;
    size_t i = 0;
    for (; i + batch <= npoints; i += batch) {
        std::array<G2AffinePoint, batch> psi1;
        std::array<G2AffinePoint, batch> psi2;
        std::array<G2AffinePoint, batch> psi3;
        avx_g2_minus_psi_affine_batch4(psi1.data(), points + i);
        avx_g2_minus_psi_affine_batch4(psi2.data(), psi1.data());
        avx_g2_minus_psi_affine_batch4(psi3.data(), psi2.data());

        for (size_t j = 0; j < batch; ++j) {
            const size_t idx = i + j;
            split_points[4 * idx] = points[idx];
            split_points[4 * idx + 1] = psi1[j];
            split_points[4 * idx + 2] = psi2[j];
            split_points[4 * idx + 3] = psi3[j];

            std::array<uint8_t, 32> scalar = {};
            std::memcpy(scalar.data(), scalars + 32 * idx, 32);
            auto parts_zz = div_by_zz(scalar);
            auto parts_low = div_by_z(parts_zz.first);
            auto parts_high = div_by_z(parts_zz.second);
            std::memcpy(split_scalars + (4 * idx) * 8, parts_low.first.data(), 8);
            std::memcpy(split_scalars + (4 * idx + 1) * 8, parts_low.second.data(), 8);
            std::memcpy(split_scalars + (4 * idx + 2) * 8, parts_high.first.data(), 8);
            std::memcpy(split_scalars + (4 * idx + 3) * 8, parts_high.second.data(), 8);
        }
    }

    for (; i < npoints; ++i) {
        split_points[4 * i] = points[i];
        split_points[4 * i + 1] = avx_g2_minus_psi_affine(split_points[4 * i]);
        split_points[4 * i + 2] = avx_g2_minus_psi_affine(split_points[4 * i + 1]);
        split_points[4 * i + 3] = avx_g2_minus_psi_affine(split_points[4 * i + 2]);

        std::array<uint8_t, 32> scalar = {};
        std::memcpy(scalar.data(), scalars + 32 * i, 32);
        auto parts_zz = div_by_zz(scalar);
        auto parts_low = div_by_z(parts_zz.first);
        auto parts_high = div_by_z(parts_zz.second);
        std::memcpy(split_scalars + (4 * i) * 8, parts_low.first.data(), 8);
        std::memcpy(split_scalars + (4 * i + 1) * 8, parts_low.second.data(), 8);
        std::memcpy(split_scalars + (4 * i + 2) * 8, parts_high.first.data(), 8);
        std::memcpy(split_scalars + (4 * i + 3) * 8, parts_high.second.data(), 8);
    }
}

inline void avx_g2_mult_pippenger_gls_prepared(
    G2Point& ret,
    const G2AffinePoint* split_points,
    const byte* split_scalars,
    size_t npoints,
    AvxG2PippengerBucket* buckets = nullptr,
    size_t window_size = 0,
    uint8_t* bucket_used = nullptr)
{
    window_size = window_size != 0 ? window_size : avx_g2_gls_pippenger_window_size(npoints);
    if (buckets != nullptr) {
        if (bucket_used != nullptr) {
            avx_msm_detail::mult_pippenger(
                ret, split_points, split_scalars, 4 * npoints, 64, buckets, bucket_used, window_size);
        } else {
            std::vector<uint8_t> local_bucket_used(static_cast<size_t>(1) << (window_size - 1));
            avx_msm_detail::mult_pippenger(
                ret, split_points, split_scalars, 4 * npoints, 64, buckets, local_bucket_used.data(), window_size);
        }
    } else {
        avx_msm_detail::mult_pippenger_alloc<G2Point, G2AffinePoint, AvxG2PippengerBucket>(
            ret, split_points, split_scalars, 4 * npoints, 64, window_size);
    }
}

inline void avx_g2_mult_pippenger_gls(
    G2Point& ret,
    const G2AffinePoint* points,
    const byte* scalars,
    size_t npoints,
    AvxG2PippengerBucket* buckets = nullptr,
    size_t window_size = 0,
    uint8_t* bucket_used = nullptr)
{
    ret = G2Point::infinity();
    if (npoints == 0 || points == nullptr || scalars == nullptr) {
        return;
    }

    std::vector<G2AffinePoint> split_points(4 * npoints);
    std::vector<byte> split_scalars(4 * npoints * 8);
    avx_g2_mult_pippenger_gls_prepare(split_points.data(), split_scalars.data(), points, scalars, npoints);
    avx_g2_mult_pippenger_gls_prepared(
        ret, split_points.data(), split_scalars.data(), npoints, buckets, window_size, bucket_used);
}

inline void avx_g2_mult_pippenger(
    G2Point& ret,
    const G2AffinePoint* points,
    const byte* scalars,
    size_t npoints,
    size_t nbits,
    AvxG2PippengerBucket* buckets = nullptr,
    size_t window_size = 0,
    uint8_t* bucket_used = nullptr)
{
    window_size = window_size != 0 ? window_size : avx_g2_pippenger_window_size(npoints);
    if (buckets != nullptr) {
        if (bucket_used != nullptr) {
            avx_msm_detail::mult_pippenger(
                ret, points, scalars, npoints, nbits, buckets, bucket_used, window_size);
        } else {
            std::vector<uint8_t> local_bucket_used(static_cast<size_t>(1) << (window_size - 1));
            avx_msm_detail::mult_pippenger(
                ret, points, scalars, npoints, nbits, buckets, local_bucket_used.data(), window_size);
        }
    } else {
        avx_msm_detail::mult_pippenger_alloc<G2Point, G2AffinePoint, AvxG2PippengerBucket>(
            ret, points, scalars, npoints, nbits, window_size);
    }
}

}  // namespace blst::avx
