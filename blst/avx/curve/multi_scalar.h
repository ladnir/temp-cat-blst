#pragma once

#include "GT.h"
#include "point.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace blst::avx {

namespace avx_msm_detail {

constexpr size_t kMaxScalarBits = 256;
constexpr size_t kMaxWindowBits = 16;
constexpr size_t kG2GlsDispatchMaxNpoints = 8192;

inline bool checked_mul_size(size_t a, size_t b, size_t& out)
{
    if (a != 0 && b > std::numeric_limits<size_t>::max() / a) {
        out = 0;
        return false;
    }
    out = a * b;
    return true;
}

inline bool checked_lshift_size(size_t value, size_t shift, size_t& out)
{
    if (shift >= sizeof(size_t) * 8 || value > (std::numeric_limits<size_t>::max() >> shift)) {
        out = 0;
        return false;
    }
    out = value << shift;
    return true;
}

inline bool window_size_is_valid(size_t window_size)
{
    return window_size != 0 && window_size <= kMaxWindowBits;
}

inline size_t clamp_window_size(size_t window_size)
{
    return std::min(window_size, kMaxWindowBits);
}

inline size_t pippenger_window_size_fallback(size_t npoints)
{
    size_t wbits = 0;
    while (npoints >>= 1) {
        ++wbits;
    }

    if (wbits > 12) {
        return clamp_window_size(wbits - 3);
    }
    if (wbits > 8) {
        return clamp_window_size(wbits - 2);
    }
    if (wbits > 4) {
        return clamp_window_size(wbits - 1);
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
    size_t split_npoints = 0;
    if (!checked_mul_size(npoints, 2, split_npoints)) {
        return kMaxWindowBits;
    }
    return pippenger_window_size_fallback(split_npoints);
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
    size_t split_npoints = 0;
    if (!checked_mul_size(npoints, 4, split_npoints)) {
        return kMaxWindowBits;
    }
    return pippenger_window_size_fallback(split_npoints);
}

inline size_t pippenger_window_size(size_t npoints)
{
    return std::max(g1_pippenger_window_size(npoints), g2_pippenger_window_size(npoints));
}

inline bool bucket_count_for_window(size_t window_size, size_t& out)
{
    if (!window_size_is_valid(window_size)) {
        out = 0;
        return false;
    }
    return checked_lshift_size(1, window_size - 1, out);
}

inline size_t bucket_count_or_zero(size_t window_size)
{
    size_t count = 0;
    return bucket_count_for_window(window_size, count) ? count : 0;
}

inline size_t bucket_scratch_sizeof(size_t npoints, size_t window_size, size_t bucket_size)
{
    if (npoints == 0) {
        return 0;
    }
    size_t count = 0;
    size_t bytes = 0;
    if (!bucket_count_for_window(window_size, count) || !checked_mul_size(count, bucket_size, bytes)) {
        return 0;
    }
    return bytes;
}

inline size_t bucket_used_scratch_sizeof(size_t npoints, size_t window_size)
{
    if (npoints == 0) {
        return 0;
    }
    return bucket_count_or_zero(window_size);
}

template<class AffinePoint>
inline bool point_count_is_addressable(size_t npoints)
{
    size_t bytes = 0;
    return checked_mul_size(npoints, sizeof(AffinePoint), bytes);
}

template<class AffinePoint>
inline bool scalar_count_is_addressable(size_t npoints, size_t nbits)
{
    if (nbits > kMaxScalarBits) {
        return false;
    }
    const size_t nbytes = (nbits + 7) / 8;
    size_t bytes = 0;
    return checked_mul_size(npoints, nbytes, bytes);
}

template<class AffinePoint>
inline bool msm_inputs_are_valid(size_t npoints, size_t nbits, size_t window_size)
{
    return window_size_is_valid(window_size) &&
        scalar_count_is_addressable<AffinePoint>(npoints, nbits) &&
        point_count_is_addressable<AffinePoint>(npoints);
}

inline bool split_count(size_t npoints, size_t factor, size_t& out)
{
    return checked_mul_size(npoints, factor, out);
}

inline size_t array_scratch_sizeof(size_t npoints, size_t factor, size_t element_size)
{
    size_t count = 0;
    size_t bytes = 0;
    if (npoints == 0 || !split_count(npoints, factor, count) || !checked_mul_size(count, element_size, bytes)) {
        return 0;
    }
    return bytes;
}

inline limb_t get_wval_limb(const byte* scalar, size_t bit0, size_t nbits)
{
    BLST_ASSERT(nbits <= 17);
    const size_t byte_index = bit0 >> 3;
    const size_t bit_shift = bit0 & 7u;
    const size_t bytes = (bit_shift + nbits + 7u) >> 3;
    BLST_ASSERT(bytes <= 3);
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
        BLST_ASSERT(window_size_is_valid(window_size));
        BLST_ASSERT(nbits <= kMaxScalarBits);
        const size_t window_begin = window_idx * window_size;
        missing_overlap = window_begin == 0 ? 1 : 0;
        extract_begin = window_begin > 0 ? window_begin - 1 : 0;
        const size_t window_end = std::min(window_begin + window_size, nbits);
        extract_size = window_end - extract_begin;
        BLST_ASSERT(extract_size != 0 && extract_size <= 17);
        extracted_mask = (static_cast<limb_t>(1) << extract_size) - 1;
        encode_size = std::min(nbits - window_begin + 1, window_size);
        BLST_ASSERT(encode_size != 0);
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
    G1AffinePoint out = point;
    G1AffinePoint::cneg(out, negate);
    return out;
}

inline G2AffinePoint maybe_negate(const G2AffinePoint& point, bool negate)
{
    G2AffinePoint out = point;
    G2AffinePoint::cneg(out, negate);
    return out;
}

using G1PippengerBucket = G1Point;
using G2PippengerBucket = G2Point;
using GtPippengerBucket = GT;

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

template<>
inline GT bucket_zero<GT>()
{
    return GT::identity();
}

template<class ProjectivePoint>
inline ProjectivePoint msm_identity()
{
    return ProjectivePoint::infinity();
}

template<>
inline GT msm_identity<GT>()
{
    return GT::identity();
}

template<class BucketPoint, class AffinePoint>
inline void bucket_add_affine(BucketPoint& out, const BucketPoint& a, const AffinePoint& point, bool subtract)
{
    BucketPoint::dadd(out, a, point, subtract);
}

inline void bucket_add_affine(G2Point& out, const G2Point& a, const G2AffinePoint& point, bool subtract)
{
    const G2AffinePoint signed_point = maybe_negate(point, subtract);
    G2Point::add(out, a, signed_point);
}

inline void bucket_add_affine(G1Point& out, const G1Point& a, const G1AffinePoint& point, bool subtract)
{
    const G1AffinePoint signed_point = maybe_negate(point, subtract);
    G1Point::add(out, a, signed_point);
}

inline void bucket_add_affine(GT& out, const GT& a, const GT& point, bool subtract)
{
    GT signed_point = point;
    GT::cneg(signed_point, subtract);
    GT::add(out, a, signed_point);
}

inline void bucket_set_affine(G1Point& out, const G1AffinePoint& point, bool subtract)
{
    const G1AffinePoint signed_point = maybe_negate(point, subtract);
    G1Point::from_affine(out, signed_point);
}

inline void bucket_set_affine(G2Point& out, const G2AffinePoint& point, bool subtract)
{
    const G2AffinePoint signed_point = maybe_negate(point, subtract);
    G2Point::from_affine(out, signed_point);
}

inline void bucket_set_affine(GT& out, const GT& point, bool subtract)
{
    out = point;
    GT::cneg(out, subtract);
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

inline void bucket_to_jacobian(GT& out, const GT& in)
{
    out = in;
}

template<class ProjectivePoint>
inline void msm_add(ProjectivePoint& out, const ProjectivePoint& a, const ProjectivePoint& b)
{
    ProjectivePoint::add(out, a, b);
}

inline void msm_double(G1Point& out, const G1Point& a)
{
    G1Point::doubling(out, a);
}

inline void msm_double(G2Point& out, const G2Point& a)
{
    G2Point::doubling(out, a);
}

inline void msm_double(GT& out, const GT& a)
{
    out.data = avx_gt_detail::narrow_to_standard(
        avx_gt_detail::fp12ring().fp12_flat_cyclotomic_sqr(a.data, avx_gt_detail::ring()));
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
    size_t n = 0;
    if (!checked_lshift_size(1, wbits, n) || n == 0) {
        out = msm_identity<ProjectivePoint>();
        BLST_ASSERT(false && "invalid AVX MSM bucket count");
        return;
    }
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
    size_t n = 0;
    if (!checked_lshift_size(1, wbits, n)) {
        out = msm_identity<ProjectivePoint>();
        BLST_ASSERT(false && "invalid AVX MSM bucket count");
        return;
    }
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
    size_t bucket_count = 0;
    if (!bucket_count_for_window(window.encode_size, bucket_count)) {
        ret = msm_identity<ProjectivePoint>();
        BLST_ASSERT(false && "invalid AVX MSM bucket count");
        return;
    }
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
    ret = msm_identity<ProjectivePoint>();
    if (npoints == 0 || points == nullptr || scalars == nullptr || nbits == 0) {
        return;
    }

    window_size = window_size != 0 ? window_size : pippenger_window_size(npoints);
    if (!msm_inputs_are_valid<AffinePoint>(npoints, nbits, window_size)) {
        BLST_ASSERT(false && "invalid AVX MSM shape");
        return;
    }
    const size_t num_windows = 1 + nbits / window_size;
    ProjectivePoint tile = msm_identity<ProjectivePoint>();
    uint8_t bucket_tag = 1;
    const size_t bucket_count = bucket_count_or_zero(window_size);
    if (bucket_count == 0) {
        BLST_ASSERT(false && "invalid AVX MSM bucket count");
        return;
    }
    if (bucket_used != nullptr) {
        std::fill_n(bucket_used, bucket_count, 0);
    }

    for (size_t window_index = num_windows; window_index-- > 0;) {
        points_tile_pippenger(
            tile, points, scalars, npoints, nbits, buckets, bucket_used, bucket_tag, window_index, window_size);
        msm_add(ret, ret, tile);
        if (window_index != 0) {
            for (size_t i = 0; i < window_size; ++i) {
                msm_double(ret, ret);
            }
        }
        if (bucket_used != nullptr && ++bucket_tag == 0) {
            std::fill_n(bucket_used, bucket_count, 0);
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
    const size_t bucket_count = bucket_count_or_zero(window_size);
    if (bucket_count == 0) {
        ret = msm_identity<ProjectivePoint>();
        BLST_ASSERT(false && "invalid AVX MSM bucket count");
        return;
    }
    std::vector<BucketPoint> buckets(bucket_count);
    std::vector<uint8_t> bucket_used(bucket_count);
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

inline size_t avx_gt_glv4_pippenger_window_size(size_t npoints)
{
    size_t split_npoints = 0;
    if (!avx_msm_detail::split_count(npoints, 4, split_npoints)) {
        return avx_msm_detail::kMaxWindowBits;
    }
    return avx_msm_detail::pippenger_window_size_fallback(split_npoints);
}

inline size_t avx_g1_mult_pippenger_plain_scratch_sizeof(size_t npoints)
{
    if (!avx_msm_detail::point_count_is_addressable<G1AffinePoint>(npoints)) {
        return 0;
    }
    return avx_msm_detail::bucket_scratch_sizeof(
        npoints, avx_g1_pippenger_window_size(npoints), sizeof(avx_msm_detail::G1PippengerBucket));
}

inline size_t avx_g1_mult_pippenger_plain_bucket_used_scratch_sizeof(size_t npoints)
{
    if (!avx_msm_detail::point_count_is_addressable<G1AffinePoint>(npoints)) {
        return 0;
    }
    return avx_msm_detail::bucket_used_scratch_sizeof(npoints, avx_g1_pippenger_window_size(npoints));
}

inline size_t avx_g1_mult_pippenger_glv_scratch_sizeof(size_t npoints)
{
    size_t split_npoints = 0;
    if (!avx_msm_detail::split_count(npoints, 2, split_npoints) ||
        !avx_msm_detail::point_count_is_addressable<G1AffinePoint>(split_npoints)) {
        return 0;
    }
    return avx_msm_detail::bucket_scratch_sizeof(
        split_npoints, avx_g1_glv_pippenger_window_size(npoints), sizeof(avx_msm_detail::G1PippengerBucket));
}

inline size_t avx_g1_mult_pippenger_glv_bucket_used_scratch_sizeof(size_t npoints)
{
    size_t split_npoints = 0;
    if (!avx_msm_detail::split_count(npoints, 2, split_npoints) ||
        !avx_msm_detail::point_count_is_addressable<G1AffinePoint>(split_npoints)) {
        return 0;
    }
    return avx_msm_detail::bucket_used_scratch_sizeof(split_npoints, avx_g1_glv_pippenger_window_size(npoints));
}

inline size_t avx_g1_mult_pippenger_glv_points_scratch_sizeof(size_t npoints)
{
    return avx_msm_detail::array_scratch_sizeof(npoints, 2, sizeof(G1AffinePoint));
}

inline size_t avx_g1_mult_pippenger_glv_scalars_scratch_sizeof(size_t npoints)
{
    return avx_msm_detail::array_scratch_sizeof(npoints, 1, 32);
}

inline size_t avx_g1_mult_pippenger_scratch_sizeof(size_t npoints)
{
    return std::max(
        avx_g1_mult_pippenger_plain_scratch_sizeof(npoints),
        avx_g1_mult_pippenger_glv_scratch_sizeof(npoints));
}

inline size_t avx_g1_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints)
{
    return std::max(
        avx_g1_mult_pippenger_plain_bucket_used_scratch_sizeof(npoints),
        avx_g1_mult_pippenger_glv_bucket_used_scratch_sizeof(npoints));
}

inline size_t avx_g2_mult_pippenger_plain_scratch_sizeof(size_t npoints)
{
    if (!avx_msm_detail::point_count_is_addressable<G2AffinePoint>(npoints)) {
        return 0;
    }
    return avx_msm_detail::bucket_scratch_sizeof(
        npoints, avx_g2_pippenger_window_size(npoints), sizeof(avx_msm_detail::G2PippengerBucket));
}

inline size_t avx_g2_mult_pippenger_plain_bucket_used_scratch_sizeof(size_t npoints)
{
    if (!avx_msm_detail::point_count_is_addressable<G2AffinePoint>(npoints)) {
        return 0;
    }
    return avx_msm_detail::bucket_used_scratch_sizeof(npoints, avx_g2_pippenger_window_size(npoints));
}

inline size_t avx_g2_mult_pippenger_gls_scratch_sizeof(size_t npoints)
{
    size_t split_npoints = 0;
    if (!avx_msm_detail::split_count(npoints, 4, split_npoints) ||
        !avx_msm_detail::point_count_is_addressable<G2AffinePoint>(split_npoints)) {
        return 0;
    }
    return avx_msm_detail::bucket_scratch_sizeof(
        split_npoints, avx_g2_gls_pippenger_window_size(npoints), sizeof(avx_msm_detail::G2PippengerBucket));
}

inline size_t avx_g2_mult_pippenger_gls_bucket_used_scratch_sizeof(size_t npoints)
{
    size_t split_npoints = 0;
    if (!avx_msm_detail::split_count(npoints, 4, split_npoints) ||
        !avx_msm_detail::point_count_is_addressable<G2AffinePoint>(split_npoints)) {
        return 0;
    }
    return avx_msm_detail::bucket_used_scratch_sizeof(split_npoints, avx_g2_gls_pippenger_window_size(npoints));
}

inline size_t avx_g2_mult_pippenger_gls_points_scratch_sizeof(size_t npoints)
{
    return avx_msm_detail::array_scratch_sizeof(npoints, 4, sizeof(G2AffinePoint));
}

inline size_t avx_g2_mult_pippenger_gls_scalars_scratch_sizeof(size_t npoints)
{
    return avx_msm_detail::array_scratch_sizeof(npoints, 1, 32);
}

inline size_t avx_g2_mult_pippenger_scratch_sizeof(size_t npoints)
{
    const size_t plain = avx_g2_mult_pippenger_plain_scratch_sizeof(npoints);
    if (npoints > avx_msm_detail::kG2GlsDispatchMaxNpoints) {
        return plain;
    }
    return std::max(plain, avx_g2_mult_pippenger_gls_scratch_sizeof(npoints));
}

inline size_t avx_g2_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints)
{
    const size_t plain = avx_g2_mult_pippenger_plain_bucket_used_scratch_sizeof(npoints);
    if (npoints > avx_msm_detail::kG2GlsDispatchMaxNpoints) {
        return plain;
    }
    return std::max(plain, avx_g2_mult_pippenger_gls_bucket_used_scratch_sizeof(npoints));
}

inline size_t avx_gt_mult_pippenger_glv4_scratch_sizeof(size_t npoints)
{
    size_t split_npoints = 0;
    if (!avx_msm_detail::split_count(npoints, 4, split_npoints) ||
        !avx_msm_detail::point_count_is_addressable<GT>(split_npoints)) {
        return 0;
    }
    return avx_msm_detail::bucket_scratch_sizeof(
        split_npoints, avx_gt_glv4_pippenger_window_size(npoints), sizeof(avx_msm_detail::GtPippengerBucket));
}

inline size_t avx_gt_mult_pippenger_glv4_bucket_used_scratch_sizeof(size_t npoints)
{
    size_t split_npoints = 0;
    if (!avx_msm_detail::split_count(npoints, 4, split_npoints) ||
        !avx_msm_detail::point_count_is_addressable<GT>(split_npoints)) {
        return 0;
    }
    return avx_msm_detail::bucket_used_scratch_sizeof(split_npoints, avx_gt_glv4_pippenger_window_size(npoints));
}

inline size_t avx_gt_mult_pippenger_glv4_points_scratch_sizeof(size_t npoints)
{
    return avx_msm_detail::array_scratch_sizeof(npoints, 4, sizeof(GT));
}

inline size_t avx_gt_mult_pippenger_glv4_scalars_scratch_sizeof(size_t npoints)
{
    return avx_msm_detail::array_scratch_sizeof(npoints, 4, 8);
}

inline size_t avx_gt_mult_pippenger_scratch_sizeof(size_t npoints)
{
    return avx_gt_mult_pippenger_glv4_scratch_sizeof(npoints);
}

inline size_t avx_gt_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints)
{
    return avx_gt_mult_pippenger_glv4_bucket_used_scratch_sizeof(npoints);
}

using AvxG1PippengerBucket = avx_msm_detail::G1PippengerBucket;
using AvxG2PippengerBucket = avx_msm_detail::G2PippengerBucket;
using AvxGtPippengerBucket = avx_msm_detail::GtPippengerBucket;

inline void avx_g1_mult_pippenger_glv(
    G1Point& ret,
    const G1AffinePoint* points,
    const byte* scalars,
    size_t npoints,
    AvxG1PippengerBucket* buckets,
    size_t window_size,
    uint8_t* bucket_used);

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
    if (nbits == 256 && window_size == 0) {
        avx_g1_mult_pippenger_glv(ret, points, scalars, npoints, buckets, 0, bucket_used);
        return;
    }

    window_size = window_size != 0 ? window_size : avx_g1_pippenger_window_size(npoints);
    const size_t bucket_count = avx_msm_detail::bucket_count_or_zero(window_size);
    if (bucket_count == 0) {
        ret = G1Point::infinity();
        BLST_ASSERT(false && "invalid AVX G1 MSM window size");
        return;
    }
    if (buckets != nullptr) {
        if (bucket_used != nullptr) {
            avx_msm_detail::mult_pippenger(
                ret, points, scalars, npoints, nbits, buckets, bucket_used, window_size);
        } else {
            std::vector<uint8_t> local_bucket_used(bucket_count);
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
    size_t split_npoints = 0;
    if (!avx_msm_detail::split_count(npoints, 2, split_npoints) ||
        !avx_msm_detail::point_count_is_addressable<G1AffinePoint>(npoints) ||
        !avx_msm_detail::scalar_count_is_addressable<G1AffinePoint>(npoints, 256)) {
        BLST_ASSERT(false && "invalid AVX G1 GLV prepare shape");
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

            std::array<uint8_t, 32> scalar = {};
            std::memcpy(scalar.data(), scalars + 32 * (i + j), 32);
            auto parts = div_by_zz(scalar);
            std::memcpy(split_scalars + (2 * (i + j)) * 16, parts.first.data(), 16);
            std::memcpy(split_scalars + (2 * (i + j) + 1) * 16, parts.second.data(), 16);
        }
        const auto phis = ring.template batch_modmul<batch>(xs, betas);
        for (size_t j = 0; j < batch; ++j) {
            G1AffinePoint phi = points[i + j];
            phi.X = phis[j];
            G1AffinePoint::cneg(phi, true);
            split_points[2 * (i + j) + 1] = phi;
        }
    }
    for (; i < npoints; ++i) {
        split_points[2 * i] = points[i];
        G1AffinePoint phi = points[i];
        phi.X = ring.modmul(points[i].X, glv.beta);
        G1AffinePoint::cneg(phi, true);
        split_points[2 * i + 1] = phi;

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
    size_t split_npoints = 0;
    const size_t bucket_count = avx_msm_detail::bucket_count_or_zero(window_size);
    if (!avx_msm_detail::split_count(npoints, 2, split_npoints) || bucket_count == 0) {
        ret = G1Point::infinity();
        BLST_ASSERT(false && "invalid AVX G1 GLV MSM shape");
        return;
    }
    if (buckets != nullptr) {
        if (bucket_used != nullptr) {
            avx_msm_detail::mult_pippenger(
                ret, split_points, split_scalars, split_npoints, 128, buckets, bucket_used, window_size);
        } else {
            std::vector<uint8_t> local_bucket_used(bucket_count);
            avx_msm_detail::mult_pippenger(
                ret, split_points, split_scalars, split_npoints, 128, buckets, local_bucket_used.data(), window_size);
        }
    } else {
        avx_msm_detail::mult_pippenger_alloc<G1Point, G1AffinePoint, AvxG1PippengerBucket>(
            ret, split_points, split_scalars, split_npoints, 128, window_size);
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

    size_t split_npoints = 0;
    size_t split_point_bytes = 0;
    size_t split_scalar_bytes = 0;
    if (!avx_msm_detail::split_count(npoints, 2, split_npoints) ||
        !avx_msm_detail::checked_mul_size(split_npoints, sizeof(G1AffinePoint), split_point_bytes) ||
        !avx_msm_detail::checked_mul_size(split_npoints, 16, split_scalar_bytes)) {
        BLST_ASSERT(false && "invalid AVX G1 GLV split count");
        return;
    }

    std::vector<G1AffinePoint> split_points(split_npoints);
    std::vector<byte> split_scalars(split_scalar_bytes);
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
    size_t split_npoints = 0;
    if (!avx_msm_detail::split_count(npoints, 4, split_npoints) ||
        !avx_msm_detail::point_count_is_addressable<G2AffinePoint>(npoints) ||
        !avx_msm_detail::scalar_count_is_addressable<G2AffinePoint>(npoints, 256)) {
        BLST_ASSERT(false && "invalid AVX G2 GLS prepare shape");
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
    size_t split_npoints = 0;
    const size_t bucket_count = avx_msm_detail::bucket_count_or_zero(window_size);
    if (!avx_msm_detail::split_count(npoints, 4, split_npoints) || bucket_count == 0) {
        ret = G2Point::infinity();
        BLST_ASSERT(false && "invalid AVX G2 GLS MSM shape");
        return;
    }
    if (buckets != nullptr) {
        if (bucket_used != nullptr) {
            avx_msm_detail::mult_pippenger(
                ret, split_points, split_scalars, split_npoints, 64, buckets, bucket_used, window_size);
        } else {
            std::vector<uint8_t> local_bucket_used(bucket_count);
            avx_msm_detail::mult_pippenger(
                ret, split_points, split_scalars, split_npoints, 64, buckets, local_bucket_used.data(), window_size);
        }
    } else {
        avx_msm_detail::mult_pippenger_alloc<G2Point, G2AffinePoint, AvxG2PippengerBucket>(
            ret, split_points, split_scalars, split_npoints, 64, window_size);
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

    size_t split_npoints = 0;
    size_t split_point_bytes = 0;
    size_t split_scalar_bytes = 0;
    if (!avx_msm_detail::split_count(npoints, 4, split_npoints) ||
        !avx_msm_detail::checked_mul_size(split_npoints, sizeof(G2AffinePoint), split_point_bytes) ||
        !avx_msm_detail::checked_mul_size(split_npoints, 8, split_scalar_bytes)) {
        BLST_ASSERT(false && "invalid AVX G2 GLS split count");
        return;
    }

    std::vector<G2AffinePoint> split_points(split_npoints);
    std::vector<byte> split_scalars(split_scalar_bytes);
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
    if (nbits == 256 && window_size == 0 && npoints <= avx_msm_detail::kG2GlsDispatchMaxNpoints) {
        avx_g2_mult_pippenger_gls(ret, points, scalars, npoints, buckets, 0, bucket_used);
        return;
    }

    window_size = window_size != 0 ? window_size : avx_g2_pippenger_window_size(npoints);
    const size_t bucket_count = avx_msm_detail::bucket_count_or_zero(window_size);
    if (bucket_count == 0) {
        ret = G2Point::infinity();
        BLST_ASSERT(false && "invalid AVX G2 MSM window size");
        return;
    }
    if (buckets != nullptr) {
        if (bucket_used != nullptr) {
            avx_msm_detail::mult_pippenger(
                ret, points, scalars, npoints, nbits, buckets, bucket_used, window_size);
        } else {
            std::vector<uint8_t> local_bucket_used(bucket_count);
            avx_msm_detail::mult_pippenger(
                ret, points, scalars, npoints, nbits, buckets, local_bucket_used.data(), window_size);
        }
    } else {
        avx_msm_detail::mult_pippenger_alloc<G2Point, G2AffinePoint, AvxG2PippengerBucket>(
            ret, points, scalars, npoints, nbits, window_size);
    }
}

inline void avx_gt_mult_pippenger_glv4_prepare(
    GT* split_points,
    byte* split_scalars,
    const GT* points,
    const byte* scalars,
    size_t npoints)
{
    if (npoints == 0 || split_points == nullptr || split_scalars == nullptr || points == nullptr || scalars == nullptr) {
        return;
    }
    size_t split_npoints = 0;
    if (!avx_msm_detail::split_count(npoints, 4, split_npoints) ||
        !avx_msm_detail::point_count_is_addressable<GT>(split_npoints) ||
        !avx_msm_detail::scalar_count_is_addressable<GT>(npoints, 256)) {
        BLST_ASSERT(false && "invalid AVX GT GLV4 prepare shape");
        return;
    }

    for (size_t i = 0; i < npoints; ++i) {
        const GT& point = points[i];
        std::array<uint64_t, 4> scalar = {};
        std::memcpy(scalar.data(), scalars + 32 * i, 32);

        int64_t parts[4] = {};
        avx_gt_glv4_detail::gt_split_scalar_4d(parts, scalar.data());

        split_points[4 * i] = point;
        GT::frobenius_map(split_points[4 * i + 1], point, 1);
        GT::frobenius_map(split_points[4 * i + 2], point, 2);
        GT::frobenius_map(split_points[4 * i + 3], point, 3);

        for (size_t j = 0; j < 4; ++j) {
            uint64_t magnitude = 0;
            if (parts[j] < 0) {
                GT::conjugate(split_points[4 * i + j]);
                magnitude = 0 - static_cast<uint64_t>(parts[j]);
            } else {
                magnitude = static_cast<uint64_t>(parts[j]);
            }
            std::memcpy(split_scalars + (4 * i + j) * 8, &magnitude, 8);
        }
    }
}

inline void avx_gt_mult_pippenger_glv4_prepared(
    GT& ret,
    const GT* split_points,
    const byte* split_scalars,
    size_t npoints,
    AvxGtPippengerBucket* buckets = nullptr,
    size_t window_size = 0,
    uint8_t* bucket_used = nullptr)
{
    window_size = window_size != 0 ? window_size : avx_gt_glv4_pippenger_window_size(npoints);
    size_t split_npoints = 0;
    const size_t bucket_count = avx_msm_detail::bucket_count_or_zero(window_size);
    if (!avx_msm_detail::split_count(npoints, 4, split_npoints) || bucket_count == 0) {
        ret = GT::identity();
        BLST_ASSERT(false && "invalid AVX GT GLV4 MSM shape");
        return;
    }
    if (buckets != nullptr) {
        if (bucket_used != nullptr) {
            avx_msm_detail::mult_pippenger(
                ret, split_points, split_scalars, split_npoints, 64, buckets, bucket_used, window_size);
        } else {
            std::vector<uint8_t> local_bucket_used(bucket_count);
            avx_msm_detail::mult_pippenger(
                ret, split_points, split_scalars, split_npoints, 64, buckets, local_bucket_used.data(), window_size);
        }
    } else {
        avx_msm_detail::mult_pippenger_alloc<GT, GT, AvxGtPippengerBucket>(
            ret, split_points, split_scalars, split_npoints, 64, window_size);
    }
}

inline void avx_gt_mult_pippenger_glv4(
    GT& ret,
    const GT* points,
    const byte* scalars,
    size_t npoints,
    AvxGtPippengerBucket* buckets = nullptr,
    size_t window_size = 0,
    uint8_t* bucket_used = nullptr,
    GT* split_points_scratch = nullptr,
    byte* split_scalars_scratch = nullptr)
{
    ret = GT::identity();
    if (npoints == 0 || points == nullptr || scalars == nullptr) {
        return;
    }

    size_t split_npoints = 0;
    size_t split_point_bytes = 0;
    size_t split_scalar_bytes = 0;
    if (!avx_msm_detail::split_count(npoints, 4, split_npoints) ||
        !avx_msm_detail::checked_mul_size(split_npoints, sizeof(GT), split_point_bytes) ||
        !avx_msm_detail::checked_mul_size(split_npoints, 8, split_scalar_bytes)) {
        BLST_ASSERT(false && "invalid AVX GT GLV4 split count");
        return;
    }

    if (split_points_scratch != nullptr && split_scalars_scratch != nullptr) {
        avx_gt_mult_pippenger_glv4_prepare(split_points_scratch, split_scalars_scratch, points, scalars, npoints);
        avx_gt_mult_pippenger_glv4_prepared(
            ret, split_points_scratch, split_scalars_scratch, npoints, buckets, window_size, bucket_used);
        return;
    }
    if (split_points_scratch != nullptr || split_scalars_scratch != nullptr) {
        BLST_ASSERT(false && "AVX GT GLV4 split scratch must provide both point and scalar buffers");
        return;
    }

    std::vector<GT> split_points(split_npoints);
    std::vector<byte> split_scalars(split_scalar_bytes);
    avx_gt_mult_pippenger_glv4_prepare(split_points.data(), split_scalars.data(), points, scalars, npoints);
    avx_gt_mult_pippenger_glv4_prepared(
        ret, split_points.data(), split_scalars.data(), npoints, buckets, window_size, bucket_used);
}

inline void avx_gt_mult_pippenger(
    GT& ret,
    const GT* points,
    const byte* scalars,
    size_t npoints,
    size_t nbits,
    AvxGtPippengerBucket* buckets = nullptr,
    size_t window_size = 0,
    uint8_t* bucket_used = nullptr)
{
    if (nbits != 256) {
        ret = GT::identity();
        BLST_ASSERT(false && "AVX GT MSM uses the GLV4 path and requires 256-bit scalars");
        return;
    }
    avx_gt_mult_pippenger_glv4(ret, points, scalars, npoints, buckets, window_size, bucket_used);
}

}  // namespace blst::avx
