#include "blst/config.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#define BLST_AVX_USE_CPP_TYPES 1
#endif

#include "bindings/blst_avx.h"

#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/multi_scalar.h"
#include "blst/scalar/curve/point.h"
#include "blst/scalar/fields/Fp12.h"
#include "blst/sign/sha256.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/map_to_G.h"
#include "blst/avx/curve/multi_scalar.h"
#include "blst/avx/curve/point.h"
#endif

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <type_traits>

using blst::byte;
using blst::limb_t;

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
using AvxMillerResult = blst::avx::MillerResult;
using AvxPreparedG2Lines = blst::avx::PreparedLineTable;

static_assert(BLST_AVX_P1_SIZE == sizeof(AvxG1));
static_assert(BLST_AVX_ALIGNMENT >= alignof(AvxG1));
static_assert(BLST_AVX_P1_AFFINE_SIZE == sizeof(AvxG1Affine));
static_assert(BLST_AVX_ALIGNMENT >= alignof(AvxG1Affine));
static_assert(BLST_AVX_P2_SIZE == sizeof(AvxG2));
static_assert(BLST_AVX_ALIGNMENT >= alignof(AvxG2));
static_assert(BLST_AVX_P2_AFFINE_SIZE == sizeof(AvxG2Affine));
static_assert(BLST_AVX_ALIGNMENT >= alignof(AvxG2Affine));
static_assert(BLST_AVX_FP12_SIZE == sizeof(AvxGt));
static_assert(BLST_AVX_ALIGNMENT >= alignof(AvxGt));
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
static_assert(std::is_standard_layout_v<AvxPreparedG2Lines>);
static_assert(std::is_trivially_copyable_v<AvxPreparedG2Lines>);
static_assert(std::is_trivially_destructible_v<AvxPreparedG2Lines>);
static_assert(std::is_same_v<blst_avx_p1, AvxG1>);
static_assert(std::is_same_v<blst_avx_p1_affine, AvxG1Affine>);
static_assert(std::is_same_v<blst_avx_p2, AvxG2>);
static_assert(std::is_same_v<blst_avx_p2_affine, AvxG2Affine>);
static_assert(std::is_same_v<blst_avx_fp12, AvxGt>);

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

inline bool prepared_g2_storage_valid(const void* storage, const size_t storage_size)
{
    return storage != nullptr &&
        storage_size >= sizeof(AvxPreparedG2Lines) &&
        reinterpret_cast<std::uintptr_t>(storage) % alignof(AvxPreparedG2Lines) == 0;
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

constexpr size_t kBatchSubgroupSeedBytes = 32;
// The selected G1/G2 curves have no nontrivial 2-torsion. After conditioning
// on every other coefficient, at most one of {-1,0,+1} can cancel a fixed
// nonzero quotient component. Under the SHA-256 PRF model, the fresh seed below
// produces independent uniform rows, so a bad batch passes with probability at
// most (1/3)^81 < 2^-128.
constexpr size_t kBatchSubgroupRounds = 81;

inline void store_be64(byte out[8], std::uint64_t value)
{
    for (size_t i = 0; i < 8; ++i) {
        out[7 - i] = static_cast<byte>(value);
        value >>= 8;
    }
}

inline bool fill_random_seed(
    blst_avx_random_bytes_fn fill_random,
    void* random_context,
    byte seed[kBatchSubgroupSeedBytes]) noexcept
{
    if (fill_random == nullptr) {
        return false;
    }
#if defined(__cpp_exceptions) || defined(_CPPUNWIND)
    try {
        return fill_random(random_context, seed, kBatchSubgroupSeedBytes);
    }
    catch (...) {
        return false;
    }
#else
    return fill_random(random_context, seed, kBatchSubgroupSeedBytes);
#endif
}

class TernaryChallengeStream {
    static constexpr char domain[] =
        "CatBLST batch subgroup challenge v1";

    const byte* seed_;
    std::uint64_t count_;
    std::uint64_t block_ = 0;
    byte group_tag_;
    byte digest_[32] = {};
    size_t offset_ = sizeof(digest_);

    void refill()
    {
        byte metadata[17];
        metadata[0] = group_tag_;
        store_be64(metadata + 1, count_);
        store_be64(metadata + 9, block_++);

        blst::Sha256Context ctx;
        blst::sha256_init(&ctx);
        blst::sha256_update(&ctx, domain, sizeof(domain) - 1);
        blst::sha256_update(&ctx, seed_, kBatchSubgroupSeedBytes);
        blst::sha256_update(&ctx, metadata, sizeof(metadata));
        blst::sha256_final(digest_, &ctx);
        blst::sha256_clear(&ctx);
        blst::vec_zero(metadata, sizeof(metadata));
        offset_ = 0;
    }

public:
    TernaryChallengeStream(
        const byte seed[kBatchSubgroupSeedBytes],
        size_t count,
        byte group_tag)
        : seed_(seed), count_(static_cast<std::uint64_t>(count)),
          group_tag_(group_tag)
    {
        static_assert(sizeof(size_t) == sizeof(std::uint64_t));
    }

    byte next()
    {
        for (;;) {
            if (offset_ == sizeof(digest_)) {
                BLST_ASSERT(block_ != (std::numeric_limits<std::uint64_t>::max)());
                refill();
            }
            const byte candidate = digest_[offset_++];
            if (candidate < 255) {
                return static_cast<byte>(candidate % 3);
            }
        }
    }

    void clear() noexcept
    {
        blst::vec_zero(digest_, sizeof(digest_));
        seed_ = nullptr;
        count_ = 0;
        block_ = 0;
        group_tag_ = 0;
        offset_ = sizeof(digest_);
    }
};

template<typename Point, typename AddComplete, typename OnCurve, typename InSubgroup>
inline bool friendly_batch_subgroup_check_projective(
    const Point* const* points,
    const size_t count,
    blst_avx_random_bytes_fn fill_random,
    void* random_context,
    const Point& shift,
    const byte group_tag,
    AddComplete add_complete,
    OnCurve on_curve,
    InSubgroup in_subgroup)
{
    if (count == 0) {
        return true;
    }
    if (points == nullptr || fill_random == nullptr) {
        return false;
    }
    for (size_t i = 0; i < count; ++i) {
        if (points[i] == nullptr || !on_curve(*points[i])) {
            return false;
        }
    }

    byte seed[kBatchSubgroupSeedBytes] = {};
    if (!fill_random_seed(fill_random, random_context, seed)) {
        blst::vec_zero(seed, sizeof(seed));
        return false;
    }

    TernaryChallengeStream challenges(seed, count, group_tag);
    bool accepted = true;
    for (size_t round = 0; round < kBatchSubgroupRounds; ++round) {
        Point sum = shift;
        Point negative_sum = Point::identity();
        bool has_negative = false;
        for (size_t i = 0; i < count; ++i) {
            switch (challenges.next()) {
            case 0:
                break;
            case 1:
                sum = add_complete(sum, *points[i]);
                break;
            case 2:
                negative_sum = add_complete(negative_sum, *points[i]);
                has_negative = true;
                break;
            default:
                BLST_ASSERT(false);
                break;
            }
        }
        if (has_negative) {
            Point::neg(negative_sum, negative_sum);
            sum = add_complete(sum, negative_sum);
        }
        if (!in_subgroup(sum)) {
            accepted = false;
            break;
        }
    }

    challenges.clear();
    blst::vec_zero(seed, sizeof(seed));
    return accepted;
}

template<typename Point, typename AddComplete, typename OnCurve, typename InSubgroup>
inline bool batch_subgroup_check_projective(
    const Point* const* points,
    const size_t count,
    const byte* coefficients,
    const size_t width,
    const Point& shift,
    AddComplete add_complete,
    OnCurve on_curve,
    InSubgroup in_subgroup)
{
    if (!on_curve(shift)) {
        return false;
    }
    if (width == 0 || count == 0) {
        return true;
    }
    if (points == nullptr || coefficients == nullptr) {
        return false;
    }
    for (size_t i = 0; i < count; ++i) {
        if (points[i] == nullptr || !on_curve(*points[i])) {
            return false;
        }
    }
    for (size_t round = 0; round < width; ++round) {
        Point sum = shift;
        Point negative_sum = Point::identity();
        bool has_negative = false;
        const byte* row = coefficients + round * count;
        for (size_t i = 0; i < count; ++i) {
            switch (row[i]) {
            case 0:
                break;
            case 1:
                sum = add_complete(sum, *points[i]);
                break;
            case 2:
                negative_sum = add_complete(negative_sum, *points[i]);
                has_negative = true;
                break;
            default:
                return false;
            }
        }
        if (has_negative) {
            Point::neg(negative_sum, negative_sum);
            sum = add_complete(sum, negative_sum);
        }
        if (!in_subgroup(sum)) {
            return false;
        }
    }
    return true;
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
        AvxG1::mult_vartime(out, points[0], scalars, nbits);
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
        AvxG2::mult_vartime(out, points[0], scalars, nbits);
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
    return point;
}

inline ScalarG1* as_scalar(blst_p1* point)
{
    return point;
}

inline const ScalarG1Affine* as_scalar(const blst_p1_affine* point)
{
    return point;
}

inline ScalarG1Affine* as_scalar(blst_p1_affine* point)
{
    return point;
}

inline const ScalarG2* as_scalar(const blst_p2* point)
{
    return point;
}

inline ScalarG2* as_scalar(blst_p2* point)
{
    return point;
}

inline const ScalarG2Affine* as_scalar(const blst_p2_affine* point)
{
    return point;
}

inline ScalarG2Affine* as_scalar(blst_p2_affine* point)
{
    return point;
}

inline const ScalarFp12* as_scalar(const blst_fp12* value)
{
    return value;
}

inline ScalarFp12* as_scalar(blst_fp12* value)
{
    return value;
}

inline bool is_aligned(const void* ptr, size_t alignment)
{
    return (reinterpret_cast<uintptr_t>(ptr) & (alignment - 1)) == 0;
}

inline bool checked_add_size(size_t a, size_t b, size_t& out)
{
    if (b > static_cast<size_t>(-1) - a) {
        out = 0;
        return false;
    }
    out = a + b;
    return true;
}

struct ByteRange
{
    std::uintptr_t begin;
    std::uintptr_t end;
};

inline bool make_byte_range(const void* ptr, size_t bytes, ByteRange& out)
{
    if (bytes == 0) {
        out = {};
        return true;
    }
    if (ptr == nullptr) {
        return false;
    }
    const auto begin = reinterpret_cast<std::uintptr_t>(ptr);
    if (bytes > std::numeric_limits<std::uintptr_t>::max() - begin) {
        return false;
    }
    out = {begin, begin + bytes};
    return true;
}

template<size_t N>
inline bool byte_ranges_are_disjoint(
    const std::array<const void*, N>& pointers,
    const std::array<size_t, N>& sizes)
{
    std::array<ByteRange, N> ranges = {};
    for (size_t i = 0; i < N; ++i) {
        if (!make_byte_range(pointers[i], sizes[i], ranges[i])) {
            return false;
        }
    }
    for (size_t i = 0; i < N; ++i) {
        if (sizes[i] == 0) {
            continue;
        }
        for (size_t j = i + 1; j < N; ++j) {
            if (sizes[j] != 0 && ranges[i].begin < ranges[j].end && ranges[j].begin < ranges[i].end) {
                return false;
            }
        }
    }
    return true;
}

inline bool align_up_size(size_t value, size_t alignment, size_t& out)
{
    const size_t mask = alignment - 1;
    if ((alignment & mask) != 0 || value > static_cast<size_t>(-1) - mask) {
        out = 0;
        return false;
    }
    out = (value + mask) & ~mask;
    return true;
}

inline bool gt_msm_scratch_layout(
    size_t npoints,
    size_t& bucket_bytes,
    size_t& split_points_offset,
    size_t& split_scalars_offset,
    size_t& total_bytes)
{
    bucket_bytes = blst::avx::avx_gt_mult_pippenger_glv4_scratch_sizeof(npoints);
    const size_t split_point_bytes = blst::avx::avx_gt_mult_pippenger_glv4_points_scratch_sizeof(npoints);
    const size_t split_scalar_bytes = blst::avx::avx_gt_mult_pippenger_glv4_scalars_scratch_sizeof(npoints);
    if (bucket_bytes == 0 || split_point_bytes == 0 || split_scalar_bytes == 0) {
        split_points_offset = 0;
        split_scalars_offset = 0;
        total_bytes = 0;
        return false;
    }
    if (!align_up_size(bucket_bytes, alignof(AvxGt), split_points_offset)) {
        return false;
    }
    if (!checked_add_size(split_points_offset, split_point_bytes, split_scalars_offset)) {
        return false;
    }
    if (!checked_add_size(split_scalars_offset, split_scalar_bytes, total_bytes)) {
        return false;
    }
    return true;
}

inline bool g1_subgroup_msm_scratch_layout(
    size_t npoints,
    size_t& bucket_bytes,
    size_t& split_points_offset,
    size_t& split_scalars_offset,
    size_t& total_bytes)
{
    bucket_bytes = blst::avx::avx_g1_mult_pippenger_scratch_sizeof(npoints);
    const size_t split_point_bytes =
        blst::avx::avx_g1_mult_pippenger_glv_points_scratch_sizeof(npoints);
    const size_t split_scalar_bytes =
        blst::avx::avx_g1_mult_pippenger_glv_scalars_scratch_sizeof(npoints);
    if (bucket_bytes == 0 || split_point_bytes == 0 || split_scalar_bytes == 0) {
        split_points_offset = 0;
        split_scalars_offset = 0;
        total_bytes = 0;
        return false;
    }
    if (!align_up_size(bucket_bytes, alignof(AvxG1Affine), split_points_offset) ||
        !checked_add_size(split_points_offset, split_point_bytes, split_scalars_offset) ||
        !checked_add_size(split_scalars_offset, split_scalar_bytes, total_bytes)) {
        return false;
    }
    return true;
}

inline bool g2_subgroup_msm_scratch_layout(
    size_t npoints,
    size_t& bucket_bytes,
    size_t& split_points_offset,
    size_t& split_scalars_offset,
    size_t& total_bytes)
{
    bucket_bytes = blst::avx::avx_g2_mult_pippenger_scratch_sizeof(npoints);
    if (npoints > blst::avx::avx_msm_detail::kG2GlsDispatchMaxNpoints) {
        split_points_offset = 0;
        split_scalars_offset = 0;
        total_bytes = bucket_bytes;
        return bucket_bytes != 0;
    }
    const size_t split_point_bytes =
        blst::avx::avx_g2_mult_pippenger_gls_points_scratch_sizeof(npoints);
    const size_t split_scalar_bytes =
        blst::avx::avx_g2_mult_pippenger_gls_scalars_scratch_sizeof(npoints);
    if (bucket_bytes == 0 || split_point_bytes == 0 || split_scalar_bytes == 0) {
        split_points_offset = 0;
        split_scalars_offset = 0;
        total_bytes = 0;
        return false;
    }
    if (!align_up_size(bucket_bytes, alignof(AvxG2Affine), split_points_offset) ||
        !checked_add_size(split_points_offset, split_point_bytes, split_scalars_offset) ||
        !checked_add_size(split_scalars_offset, split_scalar_bytes, total_bytes)) {
        return false;
    }
    return true;
}

inline AvxGt scalar_fp12_to_avx(const ScalarFp12& value)
{
    return AvxGt(blst::avx::avx_gt_detail::from_scalar(value));
}

inline ScalarFp12 avx_fp12_to_scalar(const AvxGt& value)
{
    return blst::avx::avx_gt_detail::to_scalar(value.data);
}

inline void avx_fp12_mul(
    AvxGt& out, const AvxGt& a, const AvxGt& b)
{
    out.data = blst::avx::avx_gt_detail::fp12ring().fp12_flat_mul(
        a.data, b.data, blst::avx::avx_gt_detail::ring());
}

inline void avx_fp12_sqr(AvxGt& out, const AvxGt& a)
{
    out.data = blst::avx::avx_gt_detail::fp12ring().fp12_flat_sqr(
        a.data, blst::avx::avx_gt_detail::ring());
}

inline void avx_fp12_final_exp(AvxGt& out, const AvxGt& a)
{
    blst::avx::matrixnok::PairingInverter inverter;
    out.data = blst::avx::final_exp(
        a.data, inverter, blst::avx::avx_gt_detail::fp12ring(),
        blst::avx::avx_gt_detail::ring());
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
size_t blst_avx_p1_alignof(void) { return BLST_AVX_ALIGNMENT; }
size_t blst_avx_p1_affine_sizeof(void) { return sizeof(blst_avx_p1_affine); }
size_t blst_avx_p1_affine_alignof(void) { return BLST_AVX_ALIGNMENT; }
size_t blst_avx_p2_sizeof(void) { return sizeof(blst_avx_p2); }
size_t blst_avx_p2_alignof(void) { return BLST_AVX_ALIGNMENT; }
size_t blst_avx_p2_affine_sizeof(void) { return sizeof(blst_avx_p2_affine); }
size_t blst_avx_p2_affine_alignof(void) { return BLST_AVX_ALIGNMENT; }
size_t blst_avx_fp12_sizeof(void) { return sizeof(blst_avx_fp12); }
size_t blst_avx_fp12_alignof(void) { return BLST_AVX_ALIGNMENT; }

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
        blst::avx::avx_curve_detail::scalar_fp_to_ring(as_scalar(in)->x),
        blst::avx::avx_curve_detail::scalar_fp_to_ring(as_scalar(in)->y));
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
    as_scalar(out)->x = blst::avx::avx_curve_detail::ring_fp_to_scalar(as_avx(in)->X);
    as_scalar(out)->y = blst::avx::avx_curve_detail::ring_fp_to_scalar(as_avx(in)->Y);
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

bool blst_avx_p1_cneg(blst_avx_p1* p, bool cbit)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (p == nullptr) {
        return false;
    }
    AvxG1::cneg(*as_avx(p), cbit ? 1 : 0);
    return true;
#else
    (void)p;
    (void)cbit;
    return false;
#endif
}

bool blst_avx_p1_batch_subgroup_check(
    const blst_avx_p1* const* points,
    size_t count,
    const byte* coefficients,
    size_t width,
    const blst_avx_p1* shift)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (shift == nullptr) {
        return false;
    }
    return batch_subgroup_check_projective<AvxG1>(
        reinterpret_cast<const AvxG1* const*>(points),
        count,
        coefficients,
        width,
        *as_avx(shift),
        g1_add_complete,
        [](const AvxG1& point) { return AvxG1::on_curve(point); },
        [](const AvxG1& point) { return AvxG1::in_subgroup(point); });
#else
    (void)points;
    (void)count;
    (void)coefficients;
    (void)width;
    (void)shift;
    return false;
#endif
}

bool blst_avx_p1_mult(blst_avx_p1* out, const blst_avx_p1* p, const byte* scalar, size_t nbits)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || p == nullptr || nbits > 256 ||
        (scalar == nullptr && nbits != 0)) {
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
    size_t bucket_bytes = 0;
    size_t split_points_offset = 0;
    size_t split_scalars_offset = 0;
    size_t total_bytes = 0;
    if (!g1_subgroup_msm_scratch_layout(
            npoints, bucket_bytes, split_points_offset,
            split_scalars_offset, total_bytes)) {
        return 0;
    }
    return total_bytes;
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
    size_t bucket_bytes = 0;
    size_t split_points_offset = 0;
    size_t split_scalars_offset = 0;
    size_t required_scratch = 0;
    if (nbits == 256) {
        if (!g1_subgroup_msm_scratch_layout(
                npoints, bucket_bytes, split_points_offset,
                split_scalars_offset, required_scratch)) {
            return false;
        }
    } else {
        bucket_bytes = blst::avx::avx_g1_mult_pippenger_plain_scratch_sizeof(npoints);
        required_scratch = bucket_bytes;
    }
    const size_t required_bucket_used = nbits == 256
        ? blst::avx::avx_g1_mult_pippenger_bucket_used_scratch_sizeof(npoints)
        : blst::avx::avx_g1_mult_pippenger_plain_bucket_used_scratch_sizeof(npoints);
    size_t point_bytes = 0;
    size_t scalar_bytes = 0;
    const size_t scalar_stride = (nbits + 7) / 8;
    if (required_scratch == 0 || required_bucket_used == 0 ||
        points == nullptr || scalars == nullptr || scratch == nullptr || bucket_used == nullptr ||
        scratch_bytes < required_scratch || bucket_used_bytes < required_bucket_used ||
        !is_aligned(scratch, alignof(blst::avx::AvxG1PippengerBucket)) ||
        !blst::avx::avx_msm_detail::checked_mul_size(npoints, sizeof(AvxG1Affine), point_bytes) ||
        !blst::avx::avx_msm_detail::checked_mul_size(npoints, scalar_stride, scalar_bytes) ||
        !byte_ranges_are_disjoint<5>(
            {out, points, scalars, scratch, bucket_used},
            {sizeof(AvxG1), point_bytes, scalar_bytes, required_scratch, required_bucket_used})) {
        return false;
    }
#ifndef NDEBUG
    for (size_t i = 0; i < npoints; ++i) {
        BLST_ASSERT(AvxG1Affine::on_curve(as_avx(points)[i]));
        BLST_ASSERT(AvxG1Affine::in_subgroup(as_avx(points)[i]));
    }
#endif
    if (nbits == 256) {
        auto* scratch_bytes_base = reinterpret_cast<byte*>(scratch);
        auto* split_points = reinterpret_cast<AvxG1Affine*>(
            scratch_bytes_base + split_points_offset);
        byte* split_scalars = scratch_bytes_base + split_scalars_offset;
        blst::avx::avx_g1_mult_pippenger_glv_prepare(
            split_points, split_scalars, as_avx(points), scalars, npoints);
        blst::avx::avx_g1_mult_pippenger_glv_prepared(
            *as_avx(out), split_points, split_scalars, npoints,
            reinterpret_cast<blst::avx::AvxG1PippengerBucket*>(scratch),
            0, bucket_used);
    } else {
        blst::avx::avx_g1_mult_pippenger(
            *as_avx(out), as_avx(points), scalars, npoints, nbits,
            reinterpret_cast<blst::avx::AvxG1PippengerBucket*>(scratch),
            0, bucket_used);
    }
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

bool blst_avx_p1s_in_g1(
    const blst_avx_p1* const* points,
    size_t count,
    blst_avx_random_bytes_fn fill_random,
    void* random_context)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return friendly_batch_subgroup_check_projective<AvxG1>(
        reinterpret_cast<const AvxG1* const*>(points),
        count,
        fill_random,
        random_context,
        AvxG1::generator(),
        1,
        g1_add_complete,
        [](const AvxG1& point) { return AvxG1::on_curve(point); },
        [](const AvxG1& point) { return AvxG1::in_subgroup(point); });
#else
    (void)points;
    (void)count;
    (void)fill_random;
    (void)random_context;
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
    size_t point_bytes = 0;
    size_t scalar_bytes = 0;
    if (required_scratch == 0 || required_bucket_used == 0 ||
        required_split_points == 0 || required_split_scalars == 0 ||
        points == nullptr || scalars == nullptr || scratch == nullptr || bucket_used == nullptr ||
        split_points == nullptr || split_scalars == nullptr ||
        scratch_bytes < required_scratch || bucket_used_bytes < required_bucket_used ||
        split_points_bytes < required_split_points || split_scalars_bytes < required_split_scalars ||
        !is_aligned(scratch, alignof(blst::avx::AvxG1PippengerBucket)) ||
        !is_aligned(split_points, alignof(blst::avx::G1AffinePoint)) ||
        !blst::avx::avx_msm_detail::checked_mul_size(npoints, sizeof(AvxG1Affine), point_bytes) ||
        !blst::avx::avx_msm_detail::checked_mul_size(npoints, 32, scalar_bytes) ||
        !byte_ranges_are_disjoint<7>(
            {out, points, scalars, scratch, bucket_used, split_points, split_scalars},
            {sizeof(AvxG1), point_bytes, scalar_bytes, required_scratch,
                required_bucket_used, required_split_points, required_split_scalars})) {
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
    size_t point_bytes = 0;
    size_t scalar_bytes = 0;
    if (required_scratch == 0 || required_bucket_used == 0 ||
        required_split_points == 0 || required_split_scalars == 0 ||
        points == nullptr || scalars == nullptr || scratch == nullptr || bucket_used == nullptr ||
        split_points == nullptr || split_scalars == nullptr ||
        scratch_bytes < required_scratch || bucket_used_bytes < required_bucket_used ||
        split_points_bytes < required_split_points || split_scalars_bytes < required_split_scalars ||
        !is_aligned(scratch, alignof(blst::avx::AvxG1PippengerBucket)) ||
        !is_aligned(split_points, alignof(blst::avx::G1AffinePoint)) ||
        !blst::avx::avx_msm_detail::checked_mul_size(npoints, sizeof(AvxG1Affine), point_bytes) ||
        !blst::avx::avx_msm_detail::checked_mul_size(npoints, 32, scalar_bytes) ||
        !byte_ranges_are_disjoint<7>(
            {out, points, scalars, scratch, bucket_used, split_points, split_scalars},
            {sizeof(AvxG1), point_bytes, scalar_bytes, required_scratch,
                required_bucket_used, required_split_points, required_split_scalars})) {
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
    if (out == nullptr || nbits > 256 ||
        (npoints != 0 && nbits != 0 && (points == nullptr || scalars == nullptr))) {
        return false;
    }
    if (npoints == 0 || nbits == 0) {
        *as_avx(out) = AvxG1::identity();
        return true;
    }
    if (npoints > 1) {
        size_t point_bytes = 0;
        size_t scalar_bytes = 0;
        const size_t scratch_bytes = blst_avx_p1s_mult_pippenger_scratch_sizeof(npoints);
        const size_t scalar_stride = (nbits + 7) / 8;
        if (scratch == nullptr || scratch_bytes == 0 ||
            !blst::avx::avx_msm_detail::checked_mul_size(npoints, sizeof(AvxG1), point_bytes) ||
            !blst::avx::avx_msm_detail::checked_mul_size(npoints, scalar_stride, scalar_bytes) ||
            !byte_ranges_are_disjoint<4>(
                {out, points, scalars, scratch},
                {sizeof(AvxG1), point_bytes, scalar_bytes, scratch_bytes})) {
            return false;
        }
    }
#ifndef NDEBUG
    for (size_t i = 0; i < npoints; ++i) {
        BLST_ASSERT(AvxG1::on_curve(as_avx(points)[i]));
        BLST_ASSERT(AvxG1::in_subgroup(as_avx(points)[i]));
    }
#endif
    return g1_pippenger(
        *as_avx(out), as_avx(points), npoints, scalars, nbits, scratch);
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
    return reinterpret_cast<const blst_avx_p1*>(&AvxG1::generator());
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
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(as_scalar(in)->x),
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(as_scalar(in)->y));
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
    as_scalar(out)->x = blst::avx::avx_curve_detail::ring_fp2_to_scalar(as_avx(in)->X);
    as_scalar(out)->y = blst::avx::avx_curve_detail::ring_fp2_to_scalar(as_avx(in)->Y);
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

bool blst_avx_p2_batch_subgroup_check(
    const blst_avx_p2* const* points,
    size_t count,
    const byte* coefficients,
    size_t width,
    const blst_avx_p2* shift)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (shift == nullptr) {
        return false;
    }
    return batch_subgroup_check_projective<AvxG2>(
        reinterpret_cast<const AvxG2* const*>(points),
        count,
        coefficients,
        width,
        *as_avx(shift),
        g2_add_complete,
        [](const AvxG2& point) { return AvxG2::on_curve(point); },
        [](const AvxG2& point) { return AvxG2::in_subgroup(point); });
#else
    (void)points;
    (void)count;
    (void)coefficients;
    (void)width;
    (void)shift;
    return false;
#endif
}

bool blst_avx_p2_mult(blst_avx_p2* out, const blst_avx_p2* p, const byte* scalar, size_t nbits)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || p == nullptr || nbits > 256 ||
        (scalar == nullptr && nbits != 0)) {
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
    size_t bucket_bytes = 0;
    size_t split_points_offset = 0;
    size_t split_scalars_offset = 0;
    size_t total_bytes = 0;
    if (!g2_subgroup_msm_scratch_layout(
            npoints, bucket_bytes, split_points_offset,
            split_scalars_offset, total_bytes)) {
        return 0;
    }
    return total_bytes;
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
    const bool use_gls = nbits == 256 &&
        npoints <= blst::avx::avx_msm_detail::kG2GlsDispatchMaxNpoints;
    size_t bucket_bytes = 0;
    size_t split_points_offset = 0;
    size_t split_scalars_offset = 0;
    size_t required_scratch = 0;
    if (use_gls) {
        if (!g2_subgroup_msm_scratch_layout(
                npoints, bucket_bytes, split_points_offset,
                split_scalars_offset, required_scratch)) {
            return false;
        }
    } else {
        bucket_bytes = blst::avx::avx_g2_mult_pippenger_plain_scratch_sizeof(npoints);
        required_scratch = bucket_bytes;
    }
    const size_t required_bucket_used = use_gls
        ? blst::avx::avx_g2_mult_pippenger_bucket_used_scratch_sizeof(npoints)
        : blst::avx::avx_g2_mult_pippenger_plain_bucket_used_scratch_sizeof(npoints);
    size_t point_bytes = 0;
    size_t scalar_bytes = 0;
    const size_t scalar_stride = (nbits + 7) / 8;
    if (required_scratch == 0 || required_bucket_used == 0 ||
        points == nullptr || scalars == nullptr || scratch == nullptr || bucket_used == nullptr ||
        scratch_bytes < required_scratch || bucket_used_bytes < required_bucket_used ||
        !is_aligned(scratch, alignof(blst::avx::AvxG2PippengerBucket)) ||
        !blst::avx::avx_msm_detail::checked_mul_size(npoints, sizeof(AvxG2Affine), point_bytes) ||
        !blst::avx::avx_msm_detail::checked_mul_size(npoints, scalar_stride, scalar_bytes) ||
        !byte_ranges_are_disjoint<5>(
            {out, points, scalars, scratch, bucket_used},
            {sizeof(AvxG2), point_bytes, scalar_bytes, required_scratch, required_bucket_used})) {
        return false;
    }
#ifndef NDEBUG
    for (size_t i = 0; i < npoints; ++i) {
        BLST_ASSERT(AvxG2Affine::on_curve(as_avx(points)[i]));
        BLST_ASSERT(AvxG2Affine::in_subgroup(as_avx(points)[i]));
    }
#endif
    if (use_gls) {
        auto* scratch_bytes_base = reinterpret_cast<byte*>(scratch);
        auto* split_points = reinterpret_cast<AvxG2Affine*>(
            scratch_bytes_base + split_points_offset);
        byte* split_scalars = scratch_bytes_base + split_scalars_offset;
        blst::avx::avx_g2_mult_pippenger_gls_prepare(
            split_points, split_scalars, as_avx(points), scalars, npoints);
        blst::avx::avx_g2_mult_pippenger_gls_prepared(
            *as_avx(out), split_points, split_scalars, npoints,
            reinterpret_cast<blst::avx::AvxG2PippengerBucket*>(scratch),
            0, bucket_used);
    } else {
        blst::avx::avx_g2_mult_pippenger(
            *as_avx(out), as_avx(points), scalars, npoints, nbits,
            reinterpret_cast<blst::avx::AvxG2PippengerBucket*>(scratch),
            0, bucket_used);
    }
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

bool blst_avx_p2s_in_g2(
    const blst_avx_p2* const* points,
    size_t count,
    blst_avx_random_bytes_fn fill_random,
    void* random_context)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return friendly_batch_subgroup_check_projective<AvxG2>(
        reinterpret_cast<const AvxG2* const*>(points),
        count,
        fill_random,
        random_context,
        AvxG2::generator(),
        2,
        g2_add_complete,
        [](const AvxG2& point) { return AvxG2::on_curve(point); },
        [](const AvxG2& point) { return AvxG2::in_subgroup(point); });
#else
    (void)points;
    (void)count;
    (void)fill_random;
    (void)random_context;
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

bool blst_avx_p2_serialize(byte out[192], const blst_avx_p2* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    AvxG2::serialize(out, *as_avx(in));
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

BLST_ERROR blst_avx_p2_deserialize(blst_avx_p2* out, const byte in[192])
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return BLST_BAD_ENCODING;
    }
    AvxG2Affine affine = {};
    const blst::BLST_ERROR err = AvxG2Affine::deserialize(affine, in);
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
    if (out == nullptr || nbits > 256 ||
        (npoints != 0 && nbits != 0 && (points == nullptr || scalars == nullptr))) {
        return false;
    }
    if (npoints == 0 || nbits == 0) {
        *as_avx(out) = AvxG2::identity();
        return true;
    }
    if (npoints > 1) {
        size_t point_bytes = 0;
        size_t scalar_bytes = 0;
        const size_t scratch_bytes = blst_avx_p2s_mult_pippenger_scratch_sizeof(npoints);
        const size_t scalar_stride = (nbits + 7) / 8;
        if (scratch == nullptr || scratch_bytes == 0 ||
            !blst::avx::avx_msm_detail::checked_mul_size(npoints, sizeof(AvxG2), point_bytes) ||
            !blst::avx::avx_msm_detail::checked_mul_size(npoints, scalar_stride, scalar_bytes) ||
            !byte_ranges_are_disjoint<4>(
                {out, points, scalars, scratch},
                {sizeof(AvxG2), point_bytes, scalar_bytes, scratch_bytes})) {
            return false;
        }
    }
#ifndef NDEBUG
    for (size_t i = 0; i < npoints; ++i) {
        BLST_ASSERT(AvxG2::on_curve(as_avx(points)[i]));
        BLST_ASSERT(AvxG2::in_subgroup(as_avx(points)[i]));
    }
#endif
    return g2_pippenger(
        *as_avx(out), as_avx(points), npoints, scalars, nbits, scratch);
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
    return reinterpret_cast<const blst_avx_p2*>(&AvxG2::generator());
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

bool blst_avx_gt_serialize(byte out[12 * 48], const blst_avx_gt* in)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    AvxGt::serialize(out, *as_avx(in));
    return true;
#else
    (void)out;
    (void)in;
    return false;
#endif
}

bool blst_avx_gt_deserialize(blst_avx_gt* out, const byte in[12 * 48])
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || in == nullptr) {
        return false;
    }
    AvxGt decoded = {};
    if (!AvxGt::deserialize(decoded, in)) {
        return false;
    }
    *as_avx(out) = decoded;
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
    avx_fp12_mul(*as_avx(ret), *as_avx(a), *as_avx(b));
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
    avx_fp12_sqr(*as_avx(ret), *as_avx(a));
    return true;
#else
    (void)ret;
    (void)a;
    return false;
#endif
}

bool blst_avx_gt_inverse(blst_avx_gt* ret, const blst_avx_gt* a)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (ret == nullptr || a == nullptr) {
        return false;
    }
    AvxGt::neg(*as_avx(ret), *as_avx(a));
    return true;
#else
    (void)ret;
    (void)a;
    return false;
#endif
}

size_t blst_avx_gt_mult_pippenger_scratch_sizeof(size_t npoints)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    size_t bucket_bytes = 0;
    size_t split_points_offset = 0;
    size_t split_scalars_offset = 0;
    size_t total_bytes = 0;
    if (!gt_msm_scratch_layout(npoints, bucket_bytes, split_points_offset, split_scalars_offset, total_bytes)) {
        return 0;
    }
    return total_bytes;
#else
    (void)npoints;
    return 0;
#endif
}

size_t blst_avx_gt_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return blst::avx::avx_gt_mult_pippenger_glv4_bucket_used_scratch_sizeof(npoints);
#else
    (void)npoints;
    return 0;
#endif
}

bool blst_avx_gt_mult_pippenger(
    blst_avx_gt* out,
    const blst_avx_gt* points,
    const byte* scalars,
    size_t npoints,
    void* scratch,
    size_t scratch_bytes,
    byte* bucket_used,
    size_t bucket_used_bytes)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr) {
        return false;
    }
    if (npoints == 0) {
        *as_avx(out) = AvxGt::identity();
        return true;
    }

    size_t bucket_bytes = 0;
    size_t split_points_offset = 0;
    size_t split_scalars_offset = 0;
    size_t required_scratch = 0;
    const size_t required_bucket_used =
        blst::avx::avx_gt_mult_pippenger_glv4_bucket_used_scratch_sizeof(npoints);
    size_t point_bytes = 0;
    size_t scalar_bytes = 0;
    if (!gt_msm_scratch_layout(npoints, bucket_bytes, split_points_offset, split_scalars_offset, required_scratch) ||
        required_bucket_used == 0 ||
        points == nullptr || scalars == nullptr || scratch == nullptr || bucket_used == nullptr ||
        scratch_bytes < required_scratch || bucket_used_bytes < required_bucket_used ||
        !is_aligned(scratch, alignof(AvxGt)) ||
        !blst::avx::avx_msm_detail::checked_mul_size(npoints, sizeof(AvxGt), point_bytes) ||
        !blst::avx::avx_msm_detail::checked_mul_size(npoints, 32, scalar_bytes) ||
        !byte_ranges_are_disjoint<5>(
            {out, points, scalars, scratch, bucket_used},
            {sizeof(AvxGt), point_bytes, scalar_bytes, required_scratch, required_bucket_used})) {
        return false;
    }
#ifndef NDEBUG
    for (size_t i = 0; i < npoints; ++i)
        BLST_ASSERT(AvxGt::in_group(as_avx(points)[i]));
#endif

    auto* scratch_bytes_base = reinterpret_cast<byte*>(scratch);
    blst::avx::avx_gt_mult_pippenger_glv4(
        *as_avx(out),
        as_avx(points),
        scalars,
        npoints,
        reinterpret_cast<blst::avx::AvxGtPippengerBucket*>(scratch_bytes_base),
        0,
        bucket_used,
        reinterpret_cast<AvxGt*>(scratch_bytes_base + split_points_offset),
        scratch_bytes_base + split_scalars_offset);
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
    return false;
#endif
}

bool blst_avx_fp12_final_exp(blst_avx_fp12* ret, const blst_avx_fp12* a)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (ret == nullptr || a == nullptr) {
        return false;
    }
    avx_fp12_final_exp(*as_avx(ret), *as_avx(a));
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
    return reinterpret_cast<const blst_avx_fp12*>(&blst::avx::kGtIdentity);
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

    AvxMillerResult miller = {};
    AvxGt::miller_loop(miller, &g2_affine, &g1_affine);
    AvxGt::final_exp(*as_avx(out), miller);
    return true;
#else
    (void)out;
    (void)g1;
    (void)g2;
    return false;
#endif
}

bool blst_avx_pairing_n(
    blst_avx_fp12* out,
    const blst_avx_p1_affine* g1,
    const blst_avx_p2_affine* g2,
    size_t n)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || (n != 0 && (g1 == nullptr || g2 == nullptr))) {
        return false;
    }
    AvxMillerResult miller = {};
    AvxGt::miller_loop_n(
        miller,
        reinterpret_cast<const AvxG2Affine*>(g2),
        reinterpret_cast<const AvxG1Affine*>(g1),
        n);
    AvxGt::final_exp(*as_avx(out), miller);
    return true;
#else
    (void)out;
    (void)g1;
    (void)g2;
    (void)n;
    return false;
#endif
}

size_t blst_avx_prepared_g2_lines_sizeof(void)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return sizeof(AvxPreparedG2Lines);
#else
    return 0;
#endif
}

size_t blst_avx_prepared_g2_lines_alignof(void)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    return alignof(AvxPreparedG2Lines);
#else
    return 1;
#endif
}

bool blst_avx_prepared_g2_lines_init(
    void* storage,
    size_t storage_size,
    const blst_avx_p2* g2)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (!prepared_g2_storage_valid(storage, storage_size) || g2 == nullptr) {
        return false;
    }
    auto* prepared = ::new (storage) AvxPreparedG2Lines{};
    AvxG2Affine g2_affine = {};
    AvxG2::to_affine(g2_affine, *as_avx(g2));
    AvxGt::precompute_lines(prepared->lines, g2_affine);
    return true;
#else
    (void)storage;
    (void)storage_size;
    (void)g2;
    return false;
#endif
}

bool blst_avx_pairing_prepared_g2(
    blst_avx_fp12* out,
    const void* storage,
    size_t storage_size,
    const blst_avx_p1* g1)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    if (out == nullptr || g1 == nullptr ||
        !prepared_g2_storage_valid(storage, storage_size)) {
        return false;
    }
    const auto* prepared = reinterpret_cast<const AvxPreparedG2Lines*>(storage);
    AvxG1Affine g1_affine = {};
    AvxG1::to_affine(g1_affine, *as_avx(g1));
    AvxMillerResult miller = {};
    AvxGt::miller_loop_lines(miller, prepared->lines, g1_affine);
    AvxGt::final_exp(*as_avx(out), miller);
    return true;
#else
    (void)out;
    (void)storage;
    (void)storage_size;
    (void)g1;
    return false;
#endif
}

} // extern "C"
