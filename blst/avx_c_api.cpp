#include "bindings/blst_avx.h"

#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/multi_scalar.h"
#include "blst/scalar/curve/point.h"
#include "blst/scalar/fields/Fp12.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/map_to_G.h"
#include "blst/avx/curve/multi_scalar.h"
#include "blst/avx/curve/point.h"
#endif

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace {

using ScalarG1 = blst::scalar::G1Point;
using ScalarG1Affine = blst::scalar::G1AffinePoint;
using ScalarG2 = blst::scalar::G2Point;
using ScalarG2Affine = blst::scalar::G2AffinePoint;
using ScalarFp12 = blst::scalar::Fp12;

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
using AvxG1 = blst::avx::G1Point;
using AvxG1Affine = blst::avx::G1AffinePoint;
using AvxG2 = blst::avx::G2Point;
using AvxG2Affine = blst::avx::G2AffinePoint;
using AvxGt = blst::avx::GT;

static_assert(sizeof(blst_avx_p1) == sizeof(AvxG1));
static_assert(alignof(blst_avx_p1) >= alignof(AvxG1));
static_assert(sizeof(blst_avx_p1_affine) == sizeof(AvxG1Affine));
static_assert(alignof(blst_avx_p1_affine) >= alignof(AvxG1Affine));
static_assert(sizeof(blst_avx_p2) == sizeof(AvxG2));
static_assert(alignof(blst_avx_p2) >= alignof(AvxG2));
static_assert(sizeof(blst_avx_p2_affine) == sizeof(AvxG2Affine));
static_assert(alignof(blst_avx_p2_affine) >= alignof(AvxG2Affine));
static_assert(sizeof(blst_avx_fp12) == sizeof(AvxGt));
static_assert(alignof(blst_avx_fp12) >= alignof(AvxGt));
static_assert(std::is_standard_layout_v<AvxG1>);
static_assert(std::is_standard_layout_v<AvxG1Affine>);
static_assert(std::is_standard_layout_v<AvxG2>);
static_assert(std::is_standard_layout_v<AvxG2Affine>);
static_assert(std::is_standard_layout_v<AvxGt>);
static_assert(std::is_trivially_copyable_v<AvxG1>);
static_assert(std::is_trivially_copyable_v<AvxG1Affine>);
static_assert(std::is_trivially_copyable_v<AvxG2>);
static_assert(std::is_trivially_copyable_v<AvxG2Affine>);
static_assert(std::is_trivially_copyable_v<AvxGt>);

/*
 * The expert AVX C ABI intentionally accepts caller-allocated opaque shells
 * from bindings/blst_avx.h and reinterprets them as the internal AVX types
 * after checking size, alignment, and trivial layout properties. This is an
 * explicit ABI boundary choice for supported toolchains, not an accidental
 * internal aliasing pattern callers may extend beyond these entry points.
 */

inline AvxG1* as_avx(blst_avx_p1* point)
{
    return reinterpret_cast<AvxG1*>(point);
}

inline const AvxG1* as_avx(const blst_avx_p1* point)
{
    return reinterpret_cast<const AvxG1*>(point);
}

inline AvxG1Affine* as_avx(blst_avx_p1_affine* point)
{
    return reinterpret_cast<AvxG1Affine*>(point);
}

inline const AvxG1Affine* as_avx(const blst_avx_p1_affine* point)
{
    return reinterpret_cast<const AvxG1Affine*>(point);
}

inline AvxG2* as_avx(blst_avx_p2* point)
{
    return reinterpret_cast<AvxG2*>(point);
}

inline const AvxG2* as_avx(const blst_avx_p2* point)
{
    return reinterpret_cast<const AvxG2*>(point);
}

inline AvxG2Affine* as_avx(blst_avx_p2_affine* point)
{
    return reinterpret_cast<AvxG2Affine*>(point);
}

inline const AvxG2Affine* as_avx(const blst_avx_p2_affine* point)
{
    return reinterpret_cast<const AvxG2Affine*>(point);
}

inline AvxGt* as_avx(blst_avx_fp12* value)
{
    return reinterpret_cast<AvxGt*>(value);
}

inline const AvxGt* as_avx(const blst_avx_fp12* value)
{
    return reinterpret_cast<const AvxGt*>(value);
}

template<typename Point>
inline Point* aligned_scratch(blst::limb_t* scratch)
{
    auto address = reinterpret_cast<std::uintptr_t>(scratch);
    address = (address + alignof(Point) - 1) & ~(static_cast<std::uintptr_t>(alignof(Point)) - 1);
    return reinterpret_cast<Point*>(address);
}

template<typename Point>
inline size_t avx_pippenger_scratch_sizeof(const size_t npoints)
{
    if (npoints <= 1) {
        return 0;
    }
    const size_t window = blst::scalar::pippenger_window_size(npoints);
    return (sizeof(Point) << (window - 1)) + alignof(Point) - 1;
}

inline AvxG1 g1_add_complete(const AvxG1& a, const AvxG1& b)
{
    AvxG1 out = {};
    AvxG1::dadd(out, a, b);
    return out;
}

inline AvxG2 g2_add_complete(const AvxG2& a, const AvxG2& b)
{
    AvxG2 out = {};
    AvxG2::dadd(out, a, b);
    return out;
}

inline void g1_add_into_bucket(AvxG1& bucket, blst::limb_t booth_idx, const size_t wbits, const AvxG1& point)
{
    const bool subtract = ((booth_idx >> wbits) & 1) != 0;
    booth_idx &= (static_cast<blst::limb_t>(1) << wbits) - 1;
    if (booth_idx == 0) {
        return;
    }
    AvxG1 signed_point = point;
    if (subtract) {
        AvxG1::neg(signed_point, signed_point);
    }
    bucket = g1_add_complete(bucket, signed_point);
}

inline void g2_add_into_bucket(AvxG2& bucket, blst::limb_t booth_idx, const size_t wbits, const AvxG2& point)
{
    const bool subtract = ((booth_idx >> wbits) & 1) != 0;
    booth_idx &= (static_cast<blst::limb_t>(1) << wbits) - 1;
    if (booth_idx == 0) {
        return;
    }
    AvxG2 signed_point = point;
    if (subtract) {
        AvxG2::neg(signed_point, signed_point);
    }
    bucket = g2_add_complete(bucket, signed_point);
}

template<typename Point, typename AddComplete>
inline Point integrate_buckets(Point* buckets, const size_t wbits, AddComplete add_complete)
{
    size_t n = static_cast<size_t>(1) << wbits;
    Point acc = buckets[--n];
    Point out = buckets[n];
    buckets[n] = Point::identity();
    while (n--) {
        acc = add_complete(acc, buckets[n]);
        out = add_complete(out, acc);
        buckets[n] = Point::identity();
    }
    return out;
}

inline bool g1_pippenger(
    AvxG1& out,
    const AvxG1* points,
    const size_t npoints,
    const byte* scalars,
    const size_t nbits,
    blst::limb_t* scratch)
{
    if (npoints == 0) {
        out = AvxG1::identity();
        return true;
    }
    if (npoints == 1) {
        AvxG1::mult(out, points[0], scalars, nbits);
        return true;
    }
    if (scratch == nullptr) {
        return false;
    }

    const size_t scalar_bytes = (nbits + 7) / 8;
    const size_t window = blst::scalar::pippenger_window_size(npoints);
    AvxG1* buckets = aligned_scratch<AvxG1>(scratch);
    out = AvxG1::identity();

    const size_t num_windows = 1 + nbits / window;
    for (size_t window_index = num_windows; window_index-- > 0;) {
        blst::scalar::WindowExtract extract(window_index, window, nbits);
        const size_t bucket_count = static_cast<size_t>(1) << (extract.encode_size - 1);
        for (size_t i = 0; i < bucket_count; ++i) {
            buckets[i] = AvxG1::identity();
        }
        for (size_t i = 0; i < npoints; ++i) {
            const blst::limb_t digit = extract.extract(scalars + i * scalar_bytes);
            const blst::limb_t bucket_index = digit & ((static_cast<blst::limb_t>(1) << extract.encode_size) - 1);
            if (bucket_index != 0) {
                g1_add_into_bucket(buckets[bucket_index - 1], digit, extract.encode_size, points[i]);
            }
        }

        const AvxG1 tile = integrate_buckets<AvxG1>(buckets, extract.encode_size - 1, g1_add_complete);
        out = g1_add_complete(out, tile);
        for (size_t i = 0; i < (window_index != 0 ? window : 0); ++i) {
            AvxG1::doubling(out, out);
        }
    }
    return true;
}

inline bool g2_pippenger(
    AvxG2& out,
    const AvxG2* points,
    const size_t npoints,
    const byte* scalars,
    const size_t nbits,
    blst::limb_t* scratch)
{
    if (npoints == 0) {
        out = AvxG2::identity();
        return true;
    }
    if (npoints == 1) {
        AvxG2::mult(out, points[0], scalars, nbits);
        return true;
    }
    if (scratch == nullptr) {
        return false;
    }

    const size_t scalar_bytes = (nbits + 7) / 8;
    const size_t window = blst::scalar::pippenger_window_size(npoints);
    AvxG2* buckets = aligned_scratch<AvxG2>(scratch);
    out = AvxG2::identity();

    const size_t num_windows = 1 + nbits / window;
    for (size_t window_index = num_windows; window_index-- > 0;) {
        blst::scalar::WindowExtract extract(window_index, window, nbits);
        const size_t bucket_count = static_cast<size_t>(1) << (extract.encode_size - 1);
        for (size_t i = 0; i < bucket_count; ++i) {
            buckets[i] = AvxG2::identity();
        }
        for (size_t i = 0; i < npoints; ++i) {
            const blst::limb_t digit = extract.extract(scalars + i * scalar_bytes);
            const blst::limb_t bucket_index = digit & ((static_cast<blst::limb_t>(1) << extract.encode_size) - 1);
            if (bucket_index != 0) {
                g2_add_into_bucket(buckets[bucket_index - 1], digit, extract.encode_size, points[i]);
            }
        }

        const AvxG2 tile = integrate_buckets<AvxG2>(buckets, extract.encode_size - 1, g2_add_complete);
        out = g2_add_complete(out, tile);
        for (size_t i = 0; i < (window_index != 0 ? window : 0); ++i) {
            AvxG2::doubling(out, out);
        }
    }
    return true;
}

inline const ScalarG1* as_scalar(const blst_p1* point)
{
    return reinterpret_cast<const ScalarG1*>(point);
}

inline ScalarG1* as_scalar(blst_p1* point)
{
    return reinterpret_cast<ScalarG1*>(point);
}

inline const ScalarG1Affine* as_scalar(const blst_p1_affine* point)
{
    return reinterpret_cast<const ScalarG1Affine*>(point);
}

inline ScalarG1Affine* as_scalar(blst_p1_affine* point)
{
    return reinterpret_cast<ScalarG1Affine*>(point);
}

inline const ScalarG2* as_scalar(const blst_p2* point)
{
    return reinterpret_cast<const ScalarG2*>(point);
}

inline ScalarG2* as_scalar(blst_p2* point)
{
    return reinterpret_cast<ScalarG2*>(point);
}

inline const ScalarG2Affine* as_scalar(const blst_p2_affine* point)
{
    return reinterpret_cast<const ScalarG2Affine*>(point);
}

inline ScalarG2Affine* as_scalar(blst_p2_affine* point)
{
    return reinterpret_cast<ScalarG2Affine*>(point);
}

inline const ScalarFp12* as_scalar(const blst_fp12* value)
{
    return reinterpret_cast<const ScalarFp12*>(value);
}

inline ScalarFp12* as_scalar(blst_fp12* value)
{
    return reinterpret_cast<ScalarFp12*>(value);
}

inline bool is_aligned(const void* ptr, size_t alignment)
{
    return (reinterpret_cast<uintptr_t>(ptr) & (alignment - 1)) == 0;
}

inline AvxGt scalar_fp12_to_avx(const ScalarFp12& value)
{
    unsigned char serialized[12 * 48] = {};
    blst::scalar::GT::serialize(serialized, blst::scalar::GT(value));
    AvxGt out = {};
    AvxGt::deserialize(out, serialized);
    return out;
}

inline ScalarFp12 avx_fp12_to_scalar(const AvxGt& value)
{
    unsigned char serialized[12 * 48] = {};
    AvxGt::serialize(serialized, value);
    blst::scalar::GT scalar_gt = {};
    blst::scalar::GT::deserialize(scalar_gt, serialized);
    return scalar_gt.data;
}
#endif

} // namespace

extern "C" {

bool blst_avx_enabled(void)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return true;
#else
    return false;
#endif
}

size_t blst_avx_p1_sizeof(void) { return sizeof(blst_avx_p1); }
size_t blst_avx_p1_alignof(void) { return alignof(blst_avx_p1); }
size_t blst_avx_p1_affine_sizeof(void) { return sizeof(blst_avx_p1_affine); }
size_t blst_avx_p1_affine_alignof(void) { return alignof(blst_avx_p1_affine); }
size_t blst_avx_p2_sizeof(void) { return sizeof(blst_avx_p2); }
size_t blst_avx_p2_alignof(void) { return alignof(blst_avx_p2); }
size_t blst_avx_p2_affine_sizeof(void) { return sizeof(blst_avx_p2_affine); }
size_t blst_avx_p2_affine_alignof(void) { return alignof(blst_avx_p2_affine); }
size_t blst_avx_fp12_sizeof(void) { return sizeof(blst_avx_fp12); }
size_t blst_avx_fp12_alignof(void) { return alignof(blst_avx_fp12); }

bool blst_avx_p1_from_scalar(blst_avx_p1* out, const blst_p1* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    *as_avx(out) = blst::avx::avx_map_to_G_detail::from_scalar(*as_scalar(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_p1_from_avx(blst_p1* out, const blst_avx_p1* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    *as_scalar(out) = blst::avx::avx_map_to_G_detail::to_scalar(*as_avx(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_avx_p1_affine_from_scalar(blst_avx_p1_affine* out, const blst_p1_affine* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    *as_avx(out) = AvxG1Affine(
        blst::avx::avx_curve_detail::scalar_fp_to_ring(as_scalar(in)->X),
        blst::avx::avx_curve_detail::scalar_fp_to_ring(as_scalar(in)->Y));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_p1_affine_from_avx(blst_p1_affine* out, const blst_avx_p1_affine* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    as_scalar(out)->X = blst::avx::avx_curve_detail::ring_fp_to_scalar(as_avx(in)->X);
    as_scalar(out)->Y = blst::avx::avx_curve_detail::ring_fp_to_scalar(as_avx(in)->Y);
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_avx_p1_add(blst_avx_p1* out, const blst_avx_p1* a, const blst_avx_p1* b)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || a == nullptr || b == nullptr) {
        return false;
    }
    AvxG1::dadd(*as_avx(out), *as_avx(a), *as_avx(b));
    return true;
#else
    (void)out;
    (void)a;
    (void)b;
    return false;
#endif
}

bool blst_avx_p1_mult(blst_avx_p1* out, const blst_avx_p1* p, const byte* scalar, size_t nbits)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || p == nullptr || (scalar == nullptr && nbits != 0)) {
        return false;
    }
    AvxG1::mult(*as_avx(out), *as_avx(p), scalar, nbits);
    return true;
#else
    (void)out;
    (void)p;
    (void)scalar;
    (void)nbits;
    return false;
#endif
}

size_t blst_avx_p1_mult_pippenger_scratch_sizeof(size_t npoints)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return blst::avx::avx_g1_mult_pippenger_scratch_sizeof(npoints);
#else
    (void)npoints;
    return 0;
#endif
}

size_t blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return blst::avx::avx_g1_mult_pippenger_bucket_used_scratch_sizeof(npoints);
#else
    (void)npoints;
    return 0;
#endif
}

size_t blst_avx_p1_mult_pippenger_glv_points_scratch_sizeof(size_t npoints)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return blst::avx::avx_g1_mult_pippenger_glv_points_scratch_sizeof(npoints);
#else
    (void)npoints;
    return 0;
#endif
}

size_t blst_avx_p1_mult_pippenger_glv_scalars_scratch_sizeof(size_t npoints)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return blst::avx::avx_g1_mult_pippenger_glv_scalars_scratch_sizeof(npoints);
#else
    (void)npoints;
    return 0;
#endif
}

bool blst_avx_p1_mult_pippenger(
    blst_avx_p1* out,
    const blst_avx_p1_affine* points,
    const byte* scalars,
    size_t npoints,
    size_t nbits,
    void* scratch,
    size_t scratch_bytes,
    byte* bucket_used,
    size_t bucket_used_bytes)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || nbits > 256) {
        return false;
    }
    if (npoints == 0 || nbits == 0) {
        *as_avx(out) = AvxG1::infinity();
        return true;
    }
    const size_t required_scratch = nbits == 256
        ? blst::avx::avx_g1_mult_pippenger_scratch_sizeof(npoints)
        : blst::avx::avx_g1_mult_pippenger_plain_scratch_sizeof(npoints);
    const size_t required_bucket_used = nbits == 256
        ? blst::avx::avx_g1_mult_pippenger_bucket_used_scratch_sizeof(npoints)
        : blst::avx::avx_g1_mult_pippenger_plain_bucket_used_scratch_sizeof(npoints);
    if (required_scratch == 0 || required_bucket_used == 0 ||
        points == nullptr || scalars == nullptr || scratch == nullptr || bucket_used == nullptr ||
        scratch_bytes < required_scratch || bucket_used_bytes < required_bucket_used ||
        !is_aligned(scratch, alignof(blst::avx::AvxG1PippengerBucket))) {
        return false;
    }
    blst::avx::avx_g1_mult_pippenger(
        *as_avx(out),
        as_avx(points),
        scalars,
        npoints,
        nbits,
        reinterpret_cast<blst::avx::AvxG1PippengerBucket*>(scratch),
        0,
        bucket_used);
    return true;
#else
    (void)out;
    (void)points;
    (void)scalars;
    (void)npoints;
    (void)nbits;
    (void)scratch;
    (void)scratch_bytes;
    (void)bucket_used;
    (void)bucket_used_bytes;
    return false;
#endif
}

bool blst_avx_p1_mult_pippenger_glv_with_scratch(
    blst_avx_p1* out,
    const blst_avx_p1_affine* points,
    const byte* scalars,
    size_t npoints,
    void* scratch,
    size_t scratch_bytes,
    byte* bucket_used,
    size_t bucket_used_bytes,
    blst_avx_p1_affine* split_points,
    size_t split_points_bytes,
    byte* split_scalars,
    size_t split_scalars_bytes)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr) {
        return false;
    }
    if (npoints == 0) {
        *as_avx(out) = AvxG1::infinity();
        return true;
    }
    const size_t required_scratch = blst::avx::avx_g1_mult_pippenger_glv_scratch_sizeof(npoints);
    const size_t required_bucket_used =
        blst::avx::avx_g1_mult_pippenger_glv_bucket_used_scratch_sizeof(npoints);
    const size_t required_split_points =
        blst::avx::avx_g1_mult_pippenger_glv_points_scratch_sizeof(npoints);
    const size_t required_split_scalars =
        blst::avx::avx_g1_mult_pippenger_glv_scalars_scratch_sizeof(npoints);
    if (required_scratch == 0 || required_bucket_used == 0 ||
        required_split_points == 0 || required_split_scalars == 0 ||
        points == nullptr || scalars == nullptr || scratch == nullptr || bucket_used == nullptr ||
        split_points == nullptr || split_scalars == nullptr ||
        scratch_bytes < required_scratch || bucket_used_bytes < required_bucket_used ||
        split_points_bytes < required_split_points || split_scalars_bytes < required_split_scalars ||
        !is_aligned(scratch, alignof(blst::avx::AvxG1PippengerBucket)) ||
        !is_aligned(split_points, alignof(blst::avx::G1AffinePoint))) {
        return false;
    }
    blst::avx::avx_g1_mult_pippenger_glv_prepare(
        as_avx(split_points),
        split_scalars,
        as_avx(points),
        scalars,
        npoints);
    blst::avx::avx_g1_mult_pippenger_glv_prepared(
        *as_avx(out),
        as_avx(split_points),
        split_scalars,
        npoints,
        reinterpret_cast<blst::avx::AvxG1PippengerBucket*>(scratch),
        0,
        bucket_used);
    return true;
#else
    (void)out;
    (void)points;
    (void)scalars;
    (void)npoints;
    (void)scratch;
    (void)scratch_bytes;
    (void)bucket_used;
    (void)bucket_used_bytes;
    (void)split_points;
    (void)split_points_bytes;
    (void)split_scalars;
    (void)split_scalars_bytes;
    return false;
#endif
}

size_t blst_avx_p1_mult_pippenger_glv_window_scratch_sizeof(size_t npoints, size_t window)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    size_t split_npoints = 0;
    if (!blst::avx::avx_msm_detail::window_size_is_valid(window) ||
        !blst::avx::avx_msm_detail::split_count(npoints, 2, split_npoints) ||
        !blst::avx::avx_msm_detail::point_count_is_addressable<AvxG1Affine>(split_npoints)) {
        return 0;
    }
    return blst::avx::avx_msm_detail::bucket_scratch_sizeof(
        split_npoints, window, sizeof(blst::avx::AvxG1PippengerBucket));
#else
    (void)npoints;
    (void)window;
    return 0;
#endif
}

size_t blst_avx_p1_mult_pippenger_glv_window_bucket_used_scratch_sizeof(size_t npoints, size_t window)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    size_t split_npoints = 0;
    if (!blst::avx::avx_msm_detail::window_size_is_valid(window) ||
        !blst::avx::avx_msm_detail::split_count(npoints, 2, split_npoints) ||
        !blst::avx::avx_msm_detail::point_count_is_addressable<AvxG1Affine>(split_npoints)) {
        return 0;
    }
    return blst::avx::avx_msm_detail::bucket_used_scratch_sizeof(split_npoints, window);
#else
    (void)npoints;
    (void)window;
    return 0;
#endif
}

bool blst_avx_p1_mult_pippenger_glv_window_with_scratch(
    blst_avx_p1* out,
    const blst_avx_p1_affine* points,
    const byte* scalars,
    size_t npoints,
    size_t window,
    void* scratch,
    size_t scratch_bytes,
    byte* bucket_used,
    size_t bucket_used_bytes,
    blst_avx_p1_affine* split_points,
    size_t split_points_bytes,
    byte* split_scalars,
    size_t split_scalars_bytes)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr) {
        return false;
    }
    if (npoints == 0) {
        *as_avx(out) = AvxG1::infinity();
        return true;
    }
    const size_t required_scratch =
        blst_avx_p1_mult_pippenger_glv_window_scratch_sizeof(npoints, window);
    const size_t required_bucket_used =
        blst_avx_p1_mult_pippenger_glv_window_bucket_used_scratch_sizeof(npoints, window);
    const size_t required_split_points =
        blst::avx::avx_g1_mult_pippenger_glv_points_scratch_sizeof(npoints);
    const size_t required_split_scalars =
        blst::avx::avx_g1_mult_pippenger_glv_scalars_scratch_sizeof(npoints);
    if (required_scratch == 0 || required_bucket_used == 0 ||
        required_split_points == 0 || required_split_scalars == 0 ||
        points == nullptr || scalars == nullptr || scratch == nullptr || bucket_used == nullptr ||
        split_points == nullptr || split_scalars == nullptr ||
        scratch_bytes < required_scratch || bucket_used_bytes < required_bucket_used ||
        split_points_bytes < required_split_points || split_scalars_bytes < required_split_scalars ||
        !is_aligned(scratch, alignof(blst::avx::AvxG1PippengerBucket)) ||
        !is_aligned(split_points, alignof(blst::avx::G1AffinePoint))) {
        return false;
    }
    blst::avx::avx_g1_mult_pippenger_glv_prepare(
        as_avx(split_points),
        split_scalars,
        as_avx(points),
        scalars,
        npoints);
    blst::avx::avx_g1_mult_pippenger_glv_prepared(
        *as_avx(out),
        as_avx(split_points),
        split_scalars,
        npoints,
        reinterpret_cast<blst::avx::AvxG1PippengerBucket*>(scratch),
        window,
        bucket_used);
    return true;
#else
    (void)out;
    (void)points;
    (void)scalars;
    (void)npoints;
    (void)window;
    (void)scratch;
    (void)scratch_bytes;
    (void)bucket_used;
    (void)bucket_used_bytes;
    (void)split_points;
    (void)split_points_bytes;
    (void)split_scalars;
    (void)split_scalars_bytes;
    return false;
#endif
}

bool blst_avx_p1_to_affine(blst_avx_p1_affine* out, const blst_avx_p1* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    AvxG1::to_affine(*as_avx(out), *as_avx(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_avx_p1_from_affine(blst_avx_p1* out, const blst_avx_p1_affine* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    AvxG1::from_affine(*as_avx(out), *as_avx(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_avx_p1_in_g1(const blst_avx_p1* p)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return p != nullptr && AvxG1::in_subgroup(*as_avx(p));
#else
    (void)p;
    return false;
#endif
}

bool blst_avx_p1_is_inf(const blst_avx_p1* p)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return p != nullptr && AvxG1::is_inf(*as_avx(p));
#else
    (void)p;
    return false;
#endif
}

bool blst_avx_p1_is_equal(const blst_avx_p1* a, const blst_avx_p1* b)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return a != nullptr && b != nullptr && AvxG1::equal(*as_avx(a), *as_avx(b));
#else
    (void)a;
    (void)b;
    return false;
#endif
}

bool blst_avx_p1_compress(byte out[48], const blst_avx_p1* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    AvxG1::compress(out, *as_avx(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_avx_p1_serialize(byte out[96], const blst_avx_p1* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    AvxG1::serialize(out, *as_avx(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

BLST_ERROR blst_avx_p1_deserialize(blst_avx_p1* out, const byte in[96])
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return BLST_BAD_ENCODING;
    }
    AvxG1Affine affine = {};
    const blst::BLST_ERROR err = AvxG1Affine::deserialize(affine, in);
    if (err == blst::BLST_SUCCESS) {
        AvxG1::from_affine(*as_avx(out), affine);
    }
    return static_cast<BLST_ERROR>(err);
#else
    (void)out;
    (void)in;
    return BLST_BAD_ENCODING;
#endif
}

BLST_ERROR blst_avx_p1_uncompress(blst_avx_p1* out, const byte in[48])
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return BLST_BAD_ENCODING;
    }
    AvxG1Affine affine = {};
    const blst::BLST_ERROR err = AvxG1Affine::uncompress(affine, in);
    if (err == blst::BLST_SUCCESS) {
        AvxG1::from_affine(*as_avx(out), affine);
    }
    return static_cast<BLST_ERROR>(err);
#else
    (void)out;
    (void)in;
    return BLST_BAD_ENCODING;
#endif
}

size_t blst_avx_p1s_mult_pippenger_scratch_sizeof(size_t npoints)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return avx_pippenger_scratch_sizeof<AvxG1>(npoints);
#else
    (void)npoints;
    return 0;
#endif
}

bool blst_avx_p1s_mult_pippenger(
    blst_avx_p1* out,
    const blst_avx_p1* points,
    size_t npoints,
    const byte* scalars,
    size_t nbits,
    limb_t* scratch)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || (npoints != 0 && (points == nullptr || scalars == nullptr))) {
        return false;
    }
    return g1_pippenger(*as_avx(out), as_avx(points), npoints, scalars, nbits, scratch);
#else
    (void)out;
    (void)points;
    (void)npoints;
    (void)scalars;
    (void)nbits;
    (void)scratch;
    return false;
#endif
}

const blst_avx_p1* blst_avx_p1_generator(void)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    static const AvxG1 generator = AvxG1::generator();
    return reinterpret_cast<const blst_avx_p1*>(&generator);
#else
    return nullptr;
#endif
}

bool blst_avx_p2_from_scalar(blst_avx_p2* out, const blst_p2* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    *as_avx(out) = blst::avx::avx_map_to_G_detail::from_scalar(*as_scalar(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_p2_from_avx(blst_p2* out, const blst_avx_p2* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    *as_scalar(out) = blst::avx::avx_map_to_G_detail::to_scalar(*as_avx(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_avx_p2_affine_from_scalar(blst_avx_p2_affine* out, const blst_p2_affine* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    *as_avx(out) = AvxG2Affine(
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(as_scalar(in)->X),
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(as_scalar(in)->Y));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_p2_affine_from_avx(blst_p2_affine* out, const blst_avx_p2_affine* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    as_scalar(out)->X = blst::avx::avx_curve_detail::ring_fp2_to_scalar(as_avx(in)->X);
    as_scalar(out)->Y = blst::avx::avx_curve_detail::ring_fp2_to_scalar(as_avx(in)->Y);
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_avx_p2_add(blst_avx_p2* out, const blst_avx_p2* a, const blst_avx_p2* b)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || a == nullptr || b == nullptr) {
        return false;
    }
    AvxG2::dadd(*as_avx(out), *as_avx(a), *as_avx(b));
    return true;
#else
    (void)out;
    (void)a;
    (void)b;
    return false;
#endif
}

bool blst_avx_p2_cneg(blst_avx_p2* p, bool cbit)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (p == nullptr) {
        return false;
    }
    AvxG2::cneg(*as_avx(p), cbit ? 1 : 0);
    return true;
#else
    (void)p;
    (void)cbit;
    return false;
#endif
}

bool blst_avx_p2_mult(blst_avx_p2* out, const blst_avx_p2* p, const byte* scalar, size_t nbits)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || p == nullptr || (scalar == nullptr && nbits != 0)) {
        return false;
    }
    AvxG2::mult(*as_avx(out), *as_avx(p), scalar, nbits);
    return true;
#else
    (void)out;
    (void)p;
    (void)scalar;
    (void)nbits;
    return false;
#endif
}

size_t blst_avx_p2_mult_pippenger_scratch_sizeof(size_t npoints)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return blst::avx::avx_g2_mult_pippenger_scratch_sizeof(npoints);
#else
    (void)npoints;
    return 0;
#endif
}

size_t blst_avx_p2_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return blst::avx::avx_g2_mult_pippenger_bucket_used_scratch_sizeof(npoints);
#else
    (void)npoints;
    return 0;
#endif
}

bool blst_avx_p2_mult_pippenger(
    blst_avx_p2* out,
    const blst_avx_p2_affine* points,
    const byte* scalars,
    size_t npoints,
    size_t nbits,
    void* scratch,
    size_t scratch_bytes,
    byte* bucket_used,
    size_t bucket_used_bytes)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || nbits > 256) {
        return false;
    }
    if (npoints == 0 || nbits == 0) {
        *as_avx(out) = AvxG2::infinity();
        return true;
    }
    const size_t required_scratch = nbits == 256
        ? blst::avx::avx_g2_mult_pippenger_scratch_sizeof(npoints)
        : blst::avx::avx_g2_mult_pippenger_plain_scratch_sizeof(npoints);
    const size_t required_bucket_used = nbits == 256
        ? blst::avx::avx_g2_mult_pippenger_bucket_used_scratch_sizeof(npoints)
        : blst::avx::avx_g2_mult_pippenger_plain_bucket_used_scratch_sizeof(npoints);
    if (required_scratch == 0 || required_bucket_used == 0 ||
        points == nullptr || scalars == nullptr || scratch == nullptr || bucket_used == nullptr ||
        scratch_bytes < required_scratch || bucket_used_bytes < required_bucket_used ||
        !is_aligned(scratch, alignof(blst::avx::AvxG2PippengerBucket))) {
        return false;
    }
    blst::avx::avx_g2_mult_pippenger(
        *as_avx(out),
        as_avx(points),
        scalars,
        npoints,
        nbits,
        reinterpret_cast<blst::avx::AvxG2PippengerBucket*>(scratch),
        0,
        bucket_used);
    return true;
#else
    (void)out;
    (void)points;
    (void)scalars;
    (void)npoints;
    (void)nbits;
    (void)scratch;
    (void)scratch_bytes;
    (void)bucket_used;
    (void)bucket_used_bytes;
    return false;
#endif
}

bool blst_avx_p2_to_affine(blst_avx_p2_affine* out, const blst_avx_p2* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    AvxG2::to_affine(*as_avx(out), *as_avx(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_avx_p2_from_affine(blst_avx_p2* out, const blst_avx_p2_affine* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    AvxG2::from_affine(*as_avx(out), *as_avx(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_avx_p2_in_g2(const blst_avx_p2* p)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return p != nullptr && AvxG2::in_subgroup(*as_avx(p));
#else
    (void)p;
    return false;
#endif
}

bool blst_avx_p2_is_inf(const blst_avx_p2* p)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return p != nullptr && AvxG2::is_inf(*as_avx(p));
#else
    (void)p;
    return false;
#endif
}

bool blst_avx_p2_is_equal(const blst_avx_p2* a, const blst_avx_p2* b)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return a != nullptr && b != nullptr && AvxG2::equal(*as_avx(a), *as_avx(b));
#else
    (void)a;
    (void)b;
    return false;
#endif
}

bool blst_avx_p2_compress(byte out[96], const blst_avx_p2* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    AvxG2::compress(out, *as_avx(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

BLST_ERROR blst_avx_p2_uncompress(blst_avx_p2* out, const byte in[96])
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return BLST_BAD_ENCODING;
    }
    AvxG2Affine affine = {};
    const blst::BLST_ERROR err = AvxG2Affine::uncompress(affine, in);
    if (err == blst::BLST_SUCCESS) {
        AvxG2::from_affine(*as_avx(out), affine);
    }
    return static_cast<BLST_ERROR>(err);
#else
    (void)out;
    (void)in;
    return BLST_BAD_ENCODING;
#endif
}

size_t blst_avx_p2s_mult_pippenger_scratch_sizeof(size_t npoints)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return avx_pippenger_scratch_sizeof<AvxG2>(npoints);
#else
    (void)npoints;
    return 0;
#endif
}

bool blst_avx_p2s_mult_pippenger(
    blst_avx_p2* out,
    const blst_avx_p2* points,
    size_t npoints,
    const byte* scalars,
    size_t nbits,
    limb_t* scratch)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || (npoints != 0 && (points == nullptr || scalars == nullptr))) {
        return false;
    }
    return g2_pippenger(*as_avx(out), as_avx(points), npoints, scalars, nbits, scratch);
#else
    (void)out;
    (void)points;
    (void)npoints;
    (void)scalars;
    (void)nbits;
    (void)scratch;
    return false;
#endif
}

const blst_avx_p2* blst_avx_p2_generator(void)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    static const AvxG2 generator = AvxG2::generator();
    return reinterpret_cast<const blst_avx_p2*>(&generator);
#else
    return nullptr;
#endif
}

bool blst_avx_fp12_from_scalar(blst_avx_fp12* out, const blst_fp12* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    *as_avx(out) = scalar_fp12_to_avx(*as_scalar(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_fp12_from_avx(blst_fp12* out, const blst_avx_fp12* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    *as_scalar(out) = avx_fp12_to_scalar(*as_avx(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_avx_fp12_mul(blst_avx_fp12* ret, const blst_avx_fp12* a, const blst_avx_fp12* b)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (ret == nullptr || a == nullptr || b == nullptr) {
        return false;
    }
    AvxGt::add(*as_avx(ret), *as_avx(a), *as_avx(b));
    return true;
#else
    (void)ret;
    (void)a;
    (void)b;
    return false;
#endif
}

bool blst_avx_fp12_sqr(blst_avx_fp12* ret, const blst_avx_fp12* a)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (ret == nullptr || a == nullptr) {
        return false;
    }
    AvxGt::doubling(*as_avx(ret), *as_avx(a));
    return true;
#else
    (void)ret;
    (void)a;
    return false;
#endif
}

bool blst_avx_fp12_final_exp(blst_avx_fp12* ret, const blst_avx_fp12* a)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (ret == nullptr || a == nullptr) {
        return false;
    }
    AvxGt::final_mult(*as_avx(ret), *as_avx(a));
    return true;
#else
    (void)ret;
    (void)a;
    return false;
#endif
}

bool blst_avx_fp12_is_equal(const blst_avx_fp12* a, const blst_avx_fp12* b)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return a != nullptr && b != nullptr && AvxGt::equal(*as_avx(a), *as_avx(b));
#else
    (void)a;
    (void)b;
    return false;
#endif
}

bool blst_avx_fp12_in_group(const blst_avx_fp12* a)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return a != nullptr && AvxGt::in_group(*as_avx(a));
#else
    (void)a;
    return false;
#endif
}

const blst_avx_fp12* blst_avx_fp12_one(void)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    static const AvxGt one = AvxGt::identity();
    return reinterpret_cast<const blst_avx_fp12*>(&one);
#else
    return nullptr;
#endif
}

bool blst_avx_pairing(blst_avx_fp12* out, const blst_avx_p1* g1, const blst_avx_p2* g2)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || g1 == nullptr || g2 == nullptr) {
        return false;
    }
    AvxG1Affine g1_affine = {};
    AvxG2Affine g2_affine = {};
    AvxG1::to_affine(g1_affine, *as_avx(g1));
    AvxG2::to_affine(g2_affine, *as_avx(g2));

    AvxGt miller = {};
    AvxGt::miller_loop(miller, &g2_affine, &g1_affine);
    AvxGt::final_mult(*as_avx(out), miller);
    return true;
#else
    (void)out;
    (void)g1;
    (void)g2;
    return false;
#endif
}

} // extern "C"
