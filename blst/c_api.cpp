#include "blst/backend/backend.h"
#include "blst/common/vect.h"
#include "blst/common/rb_tree.h"
#include "blst/curve/GT.h"
#include "blst/curve/map_to_G.h"
#include "blst/curve/point.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/hash_to_field.h"
#include "blst/sign/aggregate.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/multi_scalar.h"
#include "blst/scalar/curve/point.h"
#include "bindings/blst_scalar.h"
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/map_to_G.h"
#endif

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

using blst::byte;
using blst::limb_t;

namespace blst {
size_t blst_pairing_sizeof(void);
void blst_pairing_init(PAIRING* ctx, int hash_or_encode, const void* DST, size_t DST_len);
const void* blst_pairing_get_dst(const PAIRING* ctx);
void blst_pairing_commit(PAIRING* ctx);
BLST_ERROR blst_pairing_aggregate_pk_in_g2(PAIRING* ctx, const G2AffinePoint* PK, const G1AffinePoint* signature, const void* msg, size_t msg_len, const void* aug, size_t aug_len);
BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g2(PAIRING* ctx, const G2AffinePoint* PK, size_t pk_grpchk, const G1AffinePoint* signature, size_t sig_grpchk, const void* msg, size_t msg_len, const void* aug, size_t aug_len);
BLST_ERROR blst_pairing_mul_n_aggregate_pk_in_g2(PAIRING* ctx, const G2AffinePoint* PK, const G1AffinePoint* sig, const byte* scalar, size_t nbits, const void* msg, size_t msg_len, const void* aug, size_t aug_len);
BLST_ERROR blst_pairing_chk_n_mul_n_aggr_pk_in_g2(PAIRING* ctx, const G2AffinePoint* PK, size_t pk_grpchk, const G1AffinePoint* sig, size_t sig_grpchk, const byte* scalar, size_t nbits, const void* msg, size_t msg_len, const void* aug, size_t aug_len);
BLST_ERROR blst_pairing_aggregate_pk_in_g1(PAIRING* ctx, const G1AffinePoint* PK, const G2AffinePoint* signature, const void* msg, size_t msg_len, const void* aug, size_t aug_len);
BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g1(PAIRING* ctx, const G1AffinePoint* PK, size_t pk_grpchk, const G2AffinePoint* signature, size_t sig_grpchk, const void* msg, size_t msg_len, const void* aug, size_t aug_len);
BLST_ERROR blst_pairing_mul_n_aggregate_pk_in_g1(PAIRING* ctx, const G1AffinePoint* PK, const G2AffinePoint* sig, const byte* scalar, size_t nbits, const void* msg, size_t msg_len, const void* aug, size_t aug_len);
BLST_ERROR blst_pairing_chk_n_mul_n_aggr_pk_in_g1(PAIRING* ctx, const G1AffinePoint* PK, size_t pk_grpchk, const G2AffinePoint* sig, size_t sig_grpchk, const byte* scalar, size_t nbits, const void* msg, size_t msg_len, const void* aug, size_t aug_len);
BLST_ERROR blst_pairing_merge(PAIRING* ctx, const PAIRING* ctx1);
int blst_pairing_finalverify(const PAIRING* ctx, const Fp12* GTsig);
BLST_ERROR blst_aggregate_in_g1(G1Point& out, const G1Point* in, const unsigned char* zwire);
BLST_ERROR blst_aggregate_in_g2(G2Point& out, const G2Point* in, const unsigned char* zwire);
void blst_aggregated_in_g1(MillerResult& ret, const G1AffinePoint& sig);
void blst_aggregated_in_g2(MillerResult& ret, const G2AffinePoint& sig);
BLST_ERROR blst_core_verify_pk_in_g1(const G1AffinePoint* pk, const G2AffinePoint* signature, int hash_or_encode, const void* msg, size_t msg_len, const void* DST, size_t DST_len, const void* aug, size_t aug_len);
BLST_ERROR blst_core_verify_pk_in_g2(const G2AffinePoint* pk, const G1AffinePoint* signature, int hash_or_encode, const void* msg, size_t msg_len, const void* DST, size_t DST_len, const void* aug, size_t aug_len);
void blst_pairing_raw_aggregate(PAIRING* ctx, const G2AffinePoint* q, const G1AffinePoint* p);

void blst_keygen(pow256 SK, const void* IKM, size_t IKM_len, const void* info, size_t info_len) noexcept;
void blst_keygen_v3(pow256 SK, const void* IKM, size_t IKM_len, const void* info, size_t info_len) noexcept;
void blst_keygen_v4_5(pow256 SK, const void* IKM, size_t IKM_len, const void* salt, size_t salt_len, const void* info, size_t info_len) noexcept;
void blst_keygen_v5(pow256 SK, const void* IKM, size_t IKM_len, const void* salt, size_t salt_len, const void* info, size_t info_len) noexcept;
void blst_derive_master_eip2333(pow256 SK, const void* seed, size_t seed_len) noexcept;
void blst_derive_child_eip2333(pow256 SK, const pow256 parent_SK, unsigned int child_index) noexcept;
}

namespace blst::scalar {
template<typename Field>
void points_to_affine(AffinePoint<Field> dst[], const Point<Field>* const points[], size_t npoints);
}

namespace {

using BoundaryG1 = blst_p1;
using BoundaryG1Affine = blst_p1_affine;
using BoundaryG2 = blst_p2;
using BoundaryG2Affine = blst_p2_affine;
using BoundaryFp12 = blst_fp12;
using BoundaryFp6 = blst_fp6;

using PublicG1 = blst::G1Point;
using PublicG1Affine = blst::G1AffinePoint;
using PublicG2 = blst::G2Point;
using PublicG2Affine = blst::G2AffinePoint;
using PublicGt = blst::GT;
using PublicMillerResult = blst::MillerResult;
using BoundaryPairing = blst::PAIRING;

static_assert(std::is_same_v<BoundaryG1, blst::scalar::G1Point>);
static_assert(std::is_same_v<BoundaryG1Affine, blst::scalar::G1AffinePoint>);
static_assert(std::is_same_v<BoundaryG2, blst::scalar::G2Point>);
static_assert(std::is_same_v<BoundaryG2Affine, blst::scalar::G2AffinePoint>);
static_assert(std::is_same_v<BoundaryFp12, blst::scalar::Fp12>);
static_assert(std::is_standard_layout_v<BoundaryG1>);
static_assert(std::is_standard_layout_v<BoundaryG1Affine>);
static_assert(std::is_standard_layout_v<BoundaryG2>);
static_assert(std::is_standard_layout_v<BoundaryG2Affine>);
static_assert(std::is_standard_layout_v<BoundaryFp12>);
static_assert(std::is_trivially_copyable_v<BoundaryG1>);
static_assert(std::is_trivially_copyable_v<BoundaryG1Affine>);
static_assert(std::is_trivially_copyable_v<BoundaryG2>);
static_assert(std::is_trivially_copyable_v<BoundaryG2Affine>);
static_assert(std::is_trivially_copyable_v<BoundaryFp12>);
inline PublicGt boundary_to_public_gt(const BoundaryFp12& value)
{
#if BLST_DEFAULT_AVX
    return PublicGt(blst::avx::avx_gt_detail::from_scalar(value));
#else
    return PublicGt(value);
#endif
}

inline BoundaryFp12 public_to_boundary_gt(const PublicGt& value)
{
#if BLST_DEFAULT_AVX
    return blst::avx::avx_gt_detail::to_scalar(value.data);
#else
    return value.data;
#endif
}

inline PublicG1 boundary_to_public_g1(const BoundaryG1& point)
{
#if BLST_DEFAULT_AVX
    return blst::avx::avx_map_to_G_detail::from_scalar(point);
#else
    return point;
#endif
}

inline BoundaryG1 public_to_boundary_g1(const PublicG1& point)
{
#if BLST_DEFAULT_AVX
    return blst::avx::avx_map_to_G_detail::to_scalar(point);
#else
    return point;
#endif
}

inline PublicG2 boundary_to_public_g2(const BoundaryG2& point)
{
#if BLST_DEFAULT_AVX
    return blst::avx::avx_map_to_G_detail::from_scalar(point);
#else
    return point;
#endif
}

inline BoundaryG2 public_to_boundary_g2(const PublicG2& point)
{
#if BLST_DEFAULT_AVX
    return blst::avx::avx_map_to_G_detail::to_scalar(point);
#else
    return point;
#endif
}

inline PublicG1Affine boundary_to_public_g1_affine(const BoundaryG1Affine& point)
{
#if BLST_DEFAULT_AVX
    return PublicG1Affine(
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.x),
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.y));
#else
    return point;
#endif
}

inline BoundaryG1Affine public_to_boundary_g1_affine(const PublicG1Affine& point)
{
#if BLST_DEFAULT_AVX
    return BoundaryG1Affine{
        blst::avx::avx_curve_detail::ring_fp_to_scalar(point.X),
        blst::avx::avx_curve_detail::ring_fp_to_scalar(point.Y),
    };
#else
    return point;
#endif
}

inline PublicG2Affine boundary_to_public_g2_affine(const BoundaryG2Affine& point)
{
#if BLST_DEFAULT_AVX
    return PublicG2Affine(
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.x),
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.y));
#else
    return point;
#endif
}

inline BoundaryG2Affine public_to_boundary_g2_affine(const PublicG2Affine& point)
{
#if BLST_DEFAULT_AVX
    return BoundaryG2Affine{
        blst::avx::avx_curve_detail::ring_fp2_to_scalar(point.X),
        blst::avx::avx_curve_detail::ring_fp2_to_scalar(point.Y),
    };
#else
    return point;
#endif
}

inline BoundaryPairing* as_pairing(blst_pairing* ctx)
{
    return reinterpret_cast<BoundaryPairing*>(ctx);
}

inline const BoundaryPairing* as_pairing(const blst_pairing* ctx)
{
    return reinterpret_cast<const BoundaryPairing*>(ctx);
}

template<typename Group, typename GroupZZ>
inline void c_api_point_integrate_buckets(Group& out, GroupZZ buckets[], size_t wbits)
{
    GroupZZ ret = {};
    GroupZZ acc = {};
    size_t n = static_cast<size_t>(1) << wbits;

    blst::vec_copy(&acc, &buckets[--n], sizeof(acc));
    blst::vec_copy(&ret, &buckets[n], sizeof(ret));
    blst::vec_zero(&buckets[n], sizeof(buckets[n]));

    while (n--) {
        GroupZZ::dadd(acc, acc, buckets[n]);
        GroupZZ::dadd(ret, ret, acc);
        blst::vec_zero(&buckets[n], sizeof(buckets[n]));
    }
    blst::scalar::msm_to_jacobian(out, ret);
}

inline bool c_api_valid_wbits(size_t wbits)
{
    return blst::scalar::detail::valid_wbits(wbits);
}

inline bool c_api_valid_pippenger_window(size_t window)
{
    return blst::scalar::detail::valid_pippenger_window(window);
}

inline PublicMillerResult boundary_to_public_miller(const BoundaryFp12& value)
{
#if BLST_DEFAULT_AVX
    return PublicMillerResult(blst::avx::avx_gt_detail::from_scalar(value));
#else
    return PublicMillerResult(value);
#endif
}

inline BoundaryFp12 public_miller_to_boundary(const PublicMillerResult& value)
{
#if BLST_DEFAULT_AVX
    return blst::avx::avx_gt_detail::to_scalar(value.data);
#else
    return value.data;
#endif
}

template<typename T>
inline bool c_api_has_valid_array_head(const T* const items[], size_t n)
{
    return n == 0 || (items != nullptr && items[0] != nullptr);
}

template<typename T, typename U>
inline bool c_api_has_valid_pair_heads(const T* const xs[], const U* const ys[], size_t n)
{
    return n == 0 || (xs != nullptr && ys != nullptr && xs[0] != nullptr && ys[0] != nullptr);
}

template<typename GroupZZ>
inline size_t c_api_pippenger_fr_scratch_sizeof(size_t npoints)
{
    static_assert(alignof(GroupZZ) <= alignof(blst::limb_t));
    static_assert(sizeof(blst::pow256) % alignof(GroupZZ) == 0);

    if (npoints == 0) {
        return 0;
    }

    constexpr size_t scalar_size = sizeof(blst::pow256);
    constexpr size_t max_size = std::numeric_limits<size_t>::max();
    if (npoints > max_size / scalar_size) {
        return 0;
    }
    const size_t scalar_bytes = npoints * scalar_size;
    if (npoints == 1) {
        return scalar_bytes;
    }

    const size_t window = blst::scalar::pippenger_window_size(npoints);
    if (window == 0 || window - 1 >= std::numeric_limits<size_t>::digits) {
        return 0;
    }
    const size_t bucket_count = static_cast<size_t>(1) << (window - 1);
    if (bucket_count > (max_size - scalar_bytes) / sizeof(GroupZZ)) {
        return 0;
    }
    return scalar_bytes + bucket_count * sizeof(GroupZZ);
}

inline bool c_api_scalars_from_fr(
    blst::byte* out,
    const blst::Fr scalars[],
    size_t npoints)
{
    if (npoints == 0) {
        return true;
    }
    if (out == nullptr || scalars == nullptr) {
        return false;
    }
    if (npoints > std::numeric_limits<size_t>::max() / sizeof(blst::pow256)) {
        return false;
    }
    for (size_t i = 0; i < npoints; ++i) {
        blst::Fr::to_scaler(out + i * sizeof(blst::pow256), scalars[i]);
    }
    return true;
}

template<typename Group, typename AffineGroup, typename GroupZZ>
inline bool c_api_pippenger_fr(
    Group* ret,
    const AffineGroup points[],
    const blst::Fr scalars[],
    size_t npoints,
    blst::limb_t* scratch,
    size_t scratch_size)
{
    if (ret == nullptr) {
        return false;
    }
    if (npoints == 0) {
        *ret = {};
        return true;
    }

    const size_t required = c_api_pippenger_fr_scratch_sizeof<GroupZZ>(npoints);
    if (required == 0 || points == nullptr || scalars == nullptr ||
        scratch == nullptr || scratch_size < required) {
        return false;
    }

    auto* canonical = reinterpret_cast<blst::byte*>(scratch);
    const bool converted = c_api_scalars_from_fr(canonical, scalars, npoints);
    if (!converted) {
        return false;
    }

    if (npoints == 1) {
        Group point = {};
        Group::from_affine(point, points[0]);
        Group::mult(*ret, point, canonical, 255);
        return true;
    }

    const AffineGroup* point_descriptors[] = {points, nullptr};
    const blst::byte* scalar_descriptors[] = {canonical, nullptr};
    auto* buckets = reinterpret_cast<GroupZZ*>(
        canonical + npoints * sizeof(blst::pow256));
    blst::scalar::msm_mult_pippenger<Group, AffineGroup, GroupZZ>(
        *ret,
        point_descriptors,
        npoints,
        scalar_descriptors,
        255,
        buckets,
        0);
    return true;
}

template<typename Field, size_t BatchSize>
inline size_t c_api_batch_affine_algorithm_scratch_sizeof(size_t window)
{
    return blst::scalar::msm_mult_batch_affine_scratch_sizeof<
        Field, BatchSize>(window);
}

template<typename Field>
inline size_t c_api_batch_affine_fr_scratch_sizeof(
    size_t npoints,
    size_t window,
    size_t batch_size)
{
    if (npoints == 0) {
        return 0;
    }
    if (npoints > std::numeric_limits<std::uint32_t>::max()) {
        return 0;
    }

    size_t algorithm_bytes = 0;
    switch (batch_size) {
    case 128:
        algorithm_bytes =
            c_api_batch_affine_algorithm_scratch_sizeof<Field, 128>(
                window);
        break;
    case 384:
        algorithm_bytes =
            c_api_batch_affine_algorithm_scratch_sizeof<Field, 384>(
                window);
        break;
    default:
        return 0;
    }
    if (algorithm_bytes == 0) {
        return 0;
    }

    constexpr size_t scalar_size = sizeof(blst::pow256);
    size_t scalar_bytes;
    if (!blst::scalar::detail::checked_mul_size(
            npoints, scalar_size, scalar_bytes)) {
        return 0;
    }
    if (npoints == 1) {
        return scalar_bytes;
    }
    if (algorithm_bytes >
        std::numeric_limits<size_t>::max() - scalar_bytes) {
        return 0;
    }
    return scalar_bytes + algorithm_bytes;
}

template<typename Field, typename Group, typename AffineGroup>
inline bool c_api_batch_affine_fr(
    Group* ret,
    const AffineGroup points[],
    const blst::Fr scalars[],
    size_t npoints,
    blst::limb_t* scratch,
    size_t scratch_size,
    size_t window,
    size_t batch_size)
{
    if (ret == nullptr) {
        return false;
    }
    if (npoints == 0) {
        *ret = {};
        return true;
    }

    const size_t required =
        c_api_batch_affine_fr_scratch_sizeof<Field>(
            npoints, window, batch_size);
    if (required == 0 || points == nullptr || scalars == nullptr ||
        scratch == nullptr || scratch_size < required) {
        return false;
    }

    auto* canonical = reinterpret_cast<blst::byte*>(scratch);
    if (!c_api_scalars_from_fr(canonical, scalars, npoints)) {
        return false;
    }

    if (npoints == 1) {
        Group point = {};
        Group::from_affine(point, points[0]);
        Group::mult(*ret, point, canonical, 255);
        return true;
    }

    const size_t scalar_bytes = npoints * sizeof(blst::pow256);
    auto* algorithm_scratch =
        reinterpret_cast<blst::limb_t*>(canonical + scalar_bytes);
    const size_t algorithm_scratch_size = scratch_size - scalar_bytes;
    switch (batch_size) {
    case 128:
        return blst::scalar::msm_mult_batch_affine<Field, 128>(
            *ret, points, canonical, npoints, 255,
            algorithm_scratch, algorithm_scratch_size, window);
    case 384:
        return blst::scalar::msm_mult_batch_affine<Field, 384>(
            *ret, points, canonical, npoints, 255,
            algorithm_scratch, algorithm_scratch_size, window);
    default:
        return false;
    }
}

enum class CApiMsmFrAlgorithm
{
    Pippenger,
    BatchAffine128,
    BatchAffine384,
};

struct CApiMsmFrPlan
{
    CApiMsmFrAlgorithm algorithm;
    size_t window;
};

template<typename Field>
inline constexpr size_t c_api_batch_affine_measured_max_points()
{
    if constexpr (std::is_same_v<Field, blst::Fp>) {
        return size_t{1} << 18;
    }
    else {
        static_assert(std::is_same_v<Field, blst::Fp2>);
        return size_t{1} << 14;
    }
}

template<typename Field>
inline CApiMsmFrPlan c_api_msm_fr_plan(size_t npoints)
{
    constexpr size_t batch_threshold = size_t{1} << 12;
    if (npoints < batch_threshold ||
        npoints > c_api_batch_affine_measured_max_points<Field>()) {
        return {CApiMsmFrAlgorithm::Pippenger, 0};
    }

    size_t window;
    if (npoints <= (size_t{1} << 12)) {
        window = 10;
    }
    else if (npoints <= (size_t{1} << 14)) {
        window = 11;
    }
    else if (npoints <= (size_t{1} << 15)) {
        window = 12;
    }
    else if (npoints <= (size_t{1} << 16)) {
        window = 13;
    }
    else if (npoints <= (size_t{1} << 17)) {
        window = 14;
    }
    else {
        window = 15;
    }
    const CApiMsmFrAlgorithm algorithm =
        npoints >= (size_t{1} << 18)
            ? CApiMsmFrAlgorithm::BatchAffine384
            : CApiMsmFrAlgorithm::BatchAffine128;
    return {algorithm, window};
}

template<typename Field, typename GroupZZ>
inline size_t c_api_msm_fr_scratch_sizeof(size_t max_points)
{
    if (max_points == 0) {
        return 0;
    }
    const size_t pippenger_bytes =
        c_api_pippenger_fr_scratch_sizeof<GroupZZ>(max_points);
    if (pippenger_bytes == 0) {
        return 0;
    }

    const CApiMsmFrPlan plan = c_api_msm_fr_plan<Field>(max_points);
    size_t required = pippenger_bytes;
    if (plan.algorithm != CApiMsmFrAlgorithm::Pippenger) {
        const size_t batch_size =
            plan.algorithm == CApiMsmFrAlgorithm::BatchAffine384
                ? 384
                : 128;
        const size_t batch_bytes =
            c_api_batch_affine_fr_scratch_sizeof<Field>(
                max_points, plan.window, batch_size);
        if (batch_bytes == 0) {
            return 0;
        }
        required = std::max(required, batch_bytes);
    }

    // The query is a capacity query: storage for max_points must also run
    // every smaller input. Preserve the largest measured batch plan when
    // exact-size execution has conservatively returned to Pippenger.
    const size_t measured_max =
        c_api_batch_affine_measured_max_points<Field>();
    if (max_points > measured_max) {
        const CApiMsmFrPlan measured_plan =
            c_api_msm_fr_plan<Field>(measured_max);
        const size_t batch_size =
            measured_plan.algorithm == CApiMsmFrAlgorithm::BatchAffine384
                ? 384
                : 128;
        const size_t measured_batch_bytes =
            c_api_batch_affine_fr_scratch_sizeof<Field>(
                measured_max, measured_plan.window, batch_size);
        if (measured_batch_bytes == 0) {
            return 0;
        }
        required = std::max(required, measured_batch_bytes);
    }
    return required;
}

template<typename Field, typename Group, typename AffineGroup, typename GroupZZ>
inline bool c_api_msm_fr(
    Group* ret,
    const AffineGroup points[],
    const blst::Fr scalars[],
    size_t npoints,
    blst::limb_t* scratch,
    size_t scratch_size)
{
    if (ret == nullptr) {
        return false;
    }
    if (npoints == 0) {
        *ret = {};
        return true;
    }
    const size_t required =
        c_api_msm_fr_scratch_sizeof<Field, GroupZZ>(npoints);
    if (required == 0 || points == nullptr || scalars == nullptr ||
        scratch == nullptr || scratch_size < required) {
        return false;
    }

    const CApiMsmFrPlan plan = c_api_msm_fr_plan<Field>(npoints);
    if (plan.algorithm == CApiMsmFrAlgorithm::Pippenger) {
        return c_api_pippenger_fr<Group, AffineGroup, GroupZZ>(
            ret, points, scalars, npoints, scratch, scratch_size);
    }

    // The explicit batch kernel assumes non-identity inputs. Check that
    // boundary contract once, then keep the per-window kernel unchecked.
    for (size_t i = 0; i < npoints; ++i) {
        if (AffineGroup::is_inf(points[i])) {
            return c_api_pippenger_fr<Group, AffineGroup, GroupZZ>(
                ret, points, scalars, npoints, scratch, scratch_size);
        }
    }

    const size_t batch_size =
        plan.algorithm == CApiMsmFrAlgorithm::BatchAffine384 ? 384 : 128;
    return c_api_batch_affine_fr<Field>(
        ret, points, scalars, npoints, scratch, scratch_size,
        plan.window, batch_size);
}

template<typename AffineGroup, typename GroupZZ>
inline void c_api_point_bucket(GroupZZ buckets[], blst::limb_t booth_idx, size_t wbits, const AffineGroup& point)
{
    const bool booth_sign = ((booth_idx >> wbits) & 1) != 0;

    booth_idx &= (static_cast<blst::limb_t>(1) << wbits) - 1;
    if (booth_idx--) {
        blst::scalar::msm_add_or_sub(buckets[booth_idx], buckets[booth_idx], point, booth_sign);
    }
}

template<typename GroupZZ>
inline void c_api_point_prefetch(const GroupZZ buckets[], blst::limb_t booth_idx, size_t wbits)
{
    booth_idx &= (static_cast<blst::limb_t>(1) << wbits) - 1;
    if (booth_idx--) {
#ifndef __CUDA_ARCH__
        blst::vec_prefetch(&buckets[booth_idx], sizeof(buckets[booth_idx]));
#endif
    }
}

template<typename Group, typename AffineGroup, typename GroupZZ>
inline void c_api_points_tile_pippenger(
    Group& ret,
    const AffineGroup* const points[],
    size_t npoints,
    const blst::byte* const scalars[],
    size_t nbits,
    GroupZZ buckets[],
    size_t bit0,
    size_t window)
{
    ret = {};
    if (npoints == 0 || !c_api_has_valid_pair_heads(points, scalars, npoints) ||
        buckets == nullptr || !c_api_valid_pippenger_window(window) ||
        nbits == 0 || bit0 >= nbits) {
        return;
    }

    const blst::scalar::WindowExtract we =
        blst::scalar::WindowExtract::from_bit_offset(bit0, window, nbits);
    blst::vec_zero(buckets, sizeof(buckets[0]) << (we.encode_size - 1));
    size_t i = 0;
    const blst::byte* scalar = *scalars++;
    const AffineGroup* point = *points++;
    const size_t nbytes = (nbits + 7) / 8;

    auto extract_window = [&]() -> std::pair<const AffineGroup*, blst::limb_t> {
        const blst::limb_t encoding = we.extract(scalar);
        return {point, encoding};
    };
    auto advance_input = [&]() {
        scalar = *scalars ? *scalars++ : scalar + nbytes;
        point = *points ? *points++ : point + 1;
    };

    auto next = extract_window();
    ++i;

    for (; i < npoints; ++i) {
        auto current = next;
        advance_input();
        next = extract_window();
        c_api_point_prefetch(buckets, next.second, we.encode_size);
        c_api_point_bucket(buckets, current.second, we.encode_size, *current.first);
    }

    c_api_point_bucket(buckets, next.second, we.encode_size, *next.first);
    c_api_point_integrate_buckets(ret, buckets, we.encode_size - 1);
}

} // namespace

BLST_EXTERN_C_BEGIN_

bool blst_scalars_from_fr(
    blst_scalar out[],
    const blst_fr scalars[],
    size_t npoints)
{
    return c_api_scalars_from_fr(
        reinterpret_cast<blst::byte*>(out), scalars, npoints);
}

void blst_fp12_mul(BoundaryFp12* ret, const BoundaryFp12* a, const BoundaryFp12* b)
{
    BoundaryFp12::mul(*ret, *a, *b);
}

void blst_fp12_sqr(BoundaryFp12* ret, const BoundaryFp12* a)
{
    BoundaryFp12::sqr(*ret, *a);
}

void blst_fp12_cyclotomic_sqr(BoundaryFp12* ret, const BoundaryFp12* a)
{
    BoundaryFp12::cyclotomic_sqr(*ret, *a);
}

void blst_fp12_mul_by_xy00z0(BoundaryFp12* ret, const BoundaryFp12* a, const BoundaryFp6* xy00z0)
{
    BoundaryFp12::mul_by_xy00z0(*ret, *a, *xy00z0);
}

void blst_fp12_conjugate(BoundaryFp12* a)
{
    BoundaryFp12::conjugate(*a);
}

void blst_fp12_inverse(BoundaryFp12* ret, const BoundaryFp12* a)
{
    BoundaryFp12::inverse(*ret, *a);
}

bool blst_fp12_frobenius_map(BoundaryFp12* ret, const BoundaryFp12* a, size_t n)
{
    if (n < 1 || n > 3)
        return false;
    BoundaryFp12::frobenius_map(*ret, *a, n);
    return true;
}

bool blst_fp12_is_equal(const BoundaryFp12* a, const BoundaryFp12* b)
{
    return BoundaryFp12::equal(*a, *b) != 0;
}

bool blst_fp12_is_one(const BoundaryFp12* a)
{
    return BoundaryFp12::is_one(*a) != 0;
}

bool blst_fp12_in_group(const BoundaryFp12* a)
{
    return BoundaryFp12::in_group(*a) != 0;
}

const BoundaryFp12* blst_fp12_one(void)
{
    return &BoundaryFp12::one();
}

void blst_final_exp(BoundaryFp12* ret, const BoundaryFp12* f)
{
    PublicGt out = {};
    PublicGt::final_exp(out, boundary_to_public_miller(*f));
    *ret = public_to_boundary_gt(out);
}

const BoundaryFp12* blst_gt_identity(void)
{
    return blst_fp12_one();
}

const BoundaryFp12* blst_gt_generator(void)
{
    return &blst::scalar::GT::generator().data;
}

bool blst_gt_is_equal(const BoundaryFp12* a, const BoundaryFp12* b)
{
    return PublicGt::equal(boundary_to_public_gt(*a), boundary_to_public_gt(*b)) != 0;
}

bool blst_gt_is_identity(const BoundaryFp12* p)
{
    return PublicGt::is_identity(boundary_to_public_gt(*p)) != 0;
}

void blst_gt_add(BoundaryFp12* out, const BoundaryFp12* a, const BoundaryFp12* b)
{
    PublicGt result = {};
    PublicGt::add(result, boundary_to_public_gt(*a), boundary_to_public_gt(*b));
    *out = public_to_boundary_gt(result);
}

void blst_gt_sub(BoundaryFp12* out, const BoundaryFp12* a, const BoundaryFp12* b)
{
    PublicGt result = {};
    PublicGt::sub(result, boundary_to_public_gt(*a), boundary_to_public_gt(*b));
    *out = public_to_boundary_gt(result);
}

void blst_gt_neg(BoundaryFp12* out, const BoundaryFp12* p)
{
    PublicGt result = {};
    PublicGt::neg(result, boundary_to_public_gt(*p));
    *out = public_to_boundary_gt(result);
}

void blst_gt_double(BoundaryFp12* out, const BoundaryFp12* p)
{
    PublicGt result = {};
    PublicGt::doubling(result, boundary_to_public_gt(*p));
    *out = public_to_boundary_gt(result);
}

void blst_gt_mult(BoundaryFp12* out, const BoundaryFp12* p, const blst::byte* scalar, size_t nbits)
{
    PublicGt result = {};
    PublicGt::mult(result, boundary_to_public_gt(*p), scalar, nbits);
    *out = public_to_boundary_gt(result);
}

void blst_gt_mult_vartime(BoundaryFp12* out, const BoundaryFp12* p, const blst::byte* scalar, size_t nbits)
{
    PublicGt result = {};
    PublicGt::mult_vartime(result, boundary_to_public_gt(*p), scalar, nbits);
    *out = public_to_boundary_gt(result);
}

void blst_bendian_from_fp12(unsigned char out[48 * 12], const BoundaryFp12* a)
{
    BoundaryFp12::to_bendian(out, *a);
}

bool blst_fp12_finalverify(const BoundaryFp12* gt1, const BoundaryFp12* gt2)
{
    PublicMillerResult lhs = boundary_to_public_miller(*gt1);
    const PublicMillerResult rhs = boundary_to_public_miller(*gt2);
    PublicMillerResult::conjugate(lhs);
    PublicMillerResult::mul(lhs, lhs, rhs);
    PublicGt result = {};
    PublicGt::final_exp(result, lhs);
    return PublicGt::is_identity(result) != 0;
}

void blst_p1_add(BoundaryG1* out, const BoundaryG1* a, const BoundaryG1* b)
{
    PublicG1 result = {};
    PublicG1::add(result, boundary_to_public_g1(*a), boundary_to_public_g1(*b));
    *out = public_to_boundary_g1(result);
}

void blst_p1_sub(BoundaryG1* out, const BoundaryG1* a, const BoundaryG1* b)
{
    PublicG1 rhs = boundary_to_public_g1(*b);
    PublicG1::cneg(rhs, 1);
    PublicG1 result = {};
    PublicG1::dadd(result, boundary_to_public_g1(*a), rhs);
    *out = public_to_boundary_g1(result);
}

void blst_p1_add_or_double(BoundaryG1* out, const BoundaryG1* a, const BoundaryG1* b)
{
    PublicG1 result = {};
    PublicG1::dadd(result, boundary_to_public_g1(*a), boundary_to_public_g1(*b));
    *out = public_to_boundary_g1(result);
}

void blst_p1_add_affine(BoundaryG1* out, const BoundaryG1* a, const BoundaryG1Affine* b)
{
    PublicG1 result = {};
    PublicG1::add(result, boundary_to_public_g1(*a), boundary_to_public_g1_affine(*b));
    *out = public_to_boundary_g1(result);
}

void blst_p1_sub_affine(BoundaryG1* out, const BoundaryG1* a, const BoundaryG1Affine* b)
{
    PublicG1Affine rhs = boundary_to_public_g1_affine(*b);
    PublicG1Affine::cneg(rhs, 1);
    PublicG1 result = {};
    PublicG1::dadd(result, boundary_to_public_g1(*a), rhs);
    *out = public_to_boundary_g1(result);
}

void blst_p1_add_or_double_affine(BoundaryG1* out, const BoundaryG1* a, const BoundaryG1Affine* b)
{
    PublicG1 result = {};
    PublicG1::dadd(result, boundary_to_public_g1(*a), boundary_to_public_g1_affine(*b));
    *out = public_to_boundary_g1(result);
}

void blst_p1_double(BoundaryG1* out, const BoundaryG1* a)
{
    PublicG1 result = {};
    PublicG1::doubling(result, boundary_to_public_g1(*a));
    *out = public_to_boundary_g1(result);
}

void blst_p1_mult(BoundaryG1* out, const BoundaryG1* p, const blst::byte* scalar, size_t nbits)
{
    PublicG1 result = {};
    PublicG1::mult(result, boundary_to_public_g1(*p), scalar, nbits);
    *out = public_to_boundary_g1(result);
}

void blst_p1_mult_vartime(BoundaryG1* out, const BoundaryG1* p, const blst::byte* scalar, size_t nbits)
{
    PublicG1 result = {};
    PublicG1::mult_vartime(result, boundary_to_public_g1(*p), scalar, nbits);
    *out = public_to_boundary_g1(result);
}

void blst_p1_cneg(BoundaryG1* p, bool cbit)
{
    BoundaryG1::cneg(*p, cbit ? 1 : 0);
}

void blst_p1_affine_cneg(BoundaryG1Affine* p, bool cbit)
{
    BoundaryG1Affine::cneg(*p, cbit ? 1 : 0);
}

void blst_p1_to_affine(BoundaryG1Affine* out, const BoundaryG1* in)
{
    PublicG1Affine affine = {};
    PublicG1::to_affine(affine, boundary_to_public_g1(*in));
    *out = public_to_boundary_g1_affine(affine);
}

void blst_p1_from_affine(BoundaryG1* out, const BoundaryG1Affine* in)
{
    PublicG1 point = {};
    PublicG1::from_affine(point, boundary_to_public_g1_affine(*in));
    *out = public_to_boundary_g1(point);
}

void blst_p1_from_jacobian(BoundaryG1* out, const BoundaryG1* in)
{
    BoundaryG1::from_jacobian(*out, *in);
}

bool blst_p1_on_curve(const BoundaryG1* p)
{
    return BoundaryG1::on_curve(*p) != 0;
}

bool blst_p1_in_g1(const BoundaryG1* p)
{
    return BoundaryG1::in_subgroup(*p) != 0;
}

bool blst_p1_is_equal(const BoundaryG1* a, const BoundaryG1* b)
{
    return BoundaryG1::equal(*a, *b) != 0;
}

bool blst_p1_is_inf(const BoundaryG1* a)
{
    return BoundaryG1::is_inf(*a) != 0;
}

const BoundaryG1* blst_p1_generator(void)
{
    return &BoundaryG1::generator();
}

void blst_p1_unchecked_mult(BoundaryG1* out, const BoundaryG1* p, const blst::byte* scalar, size_t nbits)
{
	BoundaryG1 result = {};
	BoundaryG1::unchecked_mult(result, *p, scalar, nbits);
	*out = result;
}

bool blst_p1_affine_on_curve(const BoundaryG1Affine* p)
{
    return BoundaryG1Affine::on_curve(*p) != 0;
}

bool blst_p1_affine_in_g1(const BoundaryG1Affine* p)
{
    return BoundaryG1Affine::in_subgroup(*p) != 0;
}

bool blst_p1_affine_is_equal(const BoundaryG1Affine* a, const BoundaryG1Affine* b)
{
    return BoundaryG1Affine::equal(*a, *b) != 0;
}

bool blst_p1_affine_is_inf(const BoundaryG1Affine* a)
{
    return BoundaryG1Affine::is_inf(*a) != 0;
}

const BoundaryG1Affine* blst_p1_affine_generator(void)
{
    return &BoundaryG1Affine::generator();
}

void blst_p1s_to_affine(BoundaryG1Affine* dst, const BoundaryG1* const points[], size_t npoints)
{
    if (npoints == 0) {
        return;
    }
    if (!c_api_has_valid_array_head(points, npoints)) {
        for (size_t i = 0; i < npoints; ++i) {
            dst[i] = {};
        }
        return;
    }
    blst::scalar::points_to_affine<blst::Fp>(dst, points, npoints);
}

void blst_p1s_add(BoundaryG1* ret, const BoundaryG1Affine* const points[], size_t npoints)
{
    if (npoints == 0 || !c_api_has_valid_array_head(points, npoints)) {
        *ret = {};
        return;
    }
    blst::scalar::points_add<blst::Fp>(*ret, points, npoints);
}

size_t blst_p1s_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints)
{
    if (!c_api_valid_wbits(wbits)) {
        return 0;
    }
    return blst::scalar::points_mult_wbits_precompute_sizeof<blst::Fp>(wbits, npoints);
}

void blst_p1s_mult_wbits_precompute(
    BoundaryG1Affine table[],
    size_t wbits,
    const BoundaryG1Affine* const points[],
    size_t npoints)
{
    if (!c_api_valid_wbits(wbits) ||
        (npoints != 0 && table == nullptr) ||
        !c_api_has_valid_array_head(points, npoints)) {
        return;
    }
    blst::scalar::msm_precompute_wbits<blst::Fp>(table, wbits, points, npoints);
}

size_t blst_p1s_mult_wbits_scratch_sizeof(size_t npoints)
{
    return blst::scalar::points_mult_wbits_scratch_sizeof<blst::Fp>(npoints);
}

void blst_p1s_mult_wbits(
    BoundaryG1* ret,
    const BoundaryG1Affine table[],
    size_t wbits,
    size_t npoints,
    const blst::byte* const scalars[],
    size_t nbits,
    blst::limb_t* scratch)
{
    if (npoints == 0 || nbits == 0) {
        *ret = {};
        return;
    }
    if (nbits > 256 || !c_api_valid_wbits(wbits) || table == nullptr ||
        scalars == nullptr || scalars[0] == nullptr || scratch == nullptr) {
        *ret = {};
        return;
    }

    size_t table_len = 0;
    if (!blst::scalar::detail::checked_shl_size(npoints, wbits - 1, table_len)) {
        *ret = {};
        return;
    }

    using ScratchNode = blst::scalar::detail::BatchAffineNode<blst::Fp>;
    const size_t scratch_len = blst::scalar::points_mult_wbits_scratch_sizeof<blst::Fp>(npoints) / sizeof(ScratchNode);
    const blst::scalar::span<const BoundaryG1Affine> table_span(table, table_len);
    blst::scalar::detail::BatchAffineWorkspace<blst::Fp> scratch_workspace(
        scratch, scratch == nullptr ? 0 : scratch_len);
    blst::scalar::msm_mult_table<blst::Fp>(
        *ret, table_span, wbits, npoints, scalars, nbits, scratch_workspace);
}

size_t blst_p1s_mult_pippenger_scratch_sizeof(size_t npoints)
{
    return blst::scalar::msm_mult_pippenger_scratch_sizeof_<blst::scalar::G1PointZZ>(npoints);
}

size_t blst_p1s_mult_pippenger_window_scratch_sizeof(size_t window)
{
    return blst::scalar::msm_mult_pippenger_window_scratch_sizeof_<
        blst::scalar::G1PointZZ>(window);
}

void blst_p1s_mult_pippenger(
    BoundaryG1* ret,
    const BoundaryG1Affine* const points[],
    size_t npoints,
    const blst::byte* const scalars[],
    size_t nbits,
    blst::limb_t* scratch)
{
    if (npoints == 0 || nbits == 0) {
        *ret = {};
        return;
    }
    if (nbits > 256 ||
        !c_api_has_valid_pair_heads(points, scalars, npoints)) {
        *ret = {};
        return;
    }
    if (npoints == 1) {
        BoundaryG1 point = {};
        BoundaryG1::from_affine(point, *points[0]);
        BoundaryG1::mult_vartime(*ret, point, scalars[0], nbits);
        return;
    }
    if (scratch == nullptr) {
        *ret = {};
        return;
    }
    blst::scalar::msm_mult_pippenger<BoundaryG1, BoundaryG1Affine, blst::scalar::Pointxyzz<blst::Fp>>(
        *ret,
        points,
        npoints,
        scalars,
        nbits,
        reinterpret_cast<blst::scalar::Pointxyzz<blst::Fp>*>(scratch),
        0);
}

size_t blst_p1s_mult_pippenger_fr_scratch_sizeof(size_t npoints)
{
    return c_api_pippenger_fr_scratch_sizeof<blst::scalar::G1PointZZ>(npoints);
}

bool blst_p1s_mult_pippenger_fr(
    BoundaryG1* ret,
    const BoundaryG1Affine points[],
    const blst::Fr scalars[],
    size_t npoints,
    blst::limb_t* scratch,
    size_t scratch_size)
{
    return c_api_pippenger_fr<BoundaryG1, BoundaryG1Affine, blst::scalar::G1PointZZ>(
        ret, points, scalars, npoints, scratch, scratch_size);
}

size_t blst_p1s_mult_fr_scratch_sizeof(size_t max_points)
{
    return c_api_msm_fr_scratch_sizeof<
        blst::Fp, blst::scalar::G1PointZZ>(max_points);
}

bool blst_p1s_mult_fr(
    BoundaryG1* ret,
    const BoundaryG1Affine points[],
    const blst::Fr scalars[],
    size_t npoints,
    blst::limb_t* scratch,
    size_t scratch_size)
{
    return c_api_msm_fr<
        blst::Fp, BoundaryG1, BoundaryG1Affine,
        blst::scalar::G1PointZZ>(
        ret, points, scalars, npoints, scratch, scratch_size);
}

size_t blst_p1s_mult_batch_affine_fr_scratch_sizeof(
    size_t npoints,
    size_t window,
    size_t batch_size)
{
    return c_api_batch_affine_fr_scratch_sizeof<blst::Fp>(
        npoints, window, batch_size);
}

bool blst_p1s_mult_batch_affine_fr(
    BoundaryG1* ret,
    const BoundaryG1Affine points[],
    const blst::Fr scalars[],
    size_t npoints,
    blst::limb_t* scratch,
    size_t scratch_size,
    size_t window,
    size_t batch_size)
{
    return c_api_batch_affine_fr<blst::Fp>(
        ret, points, scalars, npoints, scratch, scratch_size,
        window, batch_size);
}

void blst_p1s_mult_pippenger_window(
    BoundaryG1* ret,
    const BoundaryG1Affine* const points[],
    size_t npoints,
    const blst::byte* const scalars[],
    size_t nbits,
    blst::limb_t* scratch,
    size_t window)
{
    if (npoints == 0 || nbits == 0) {
        *ret = {};
        return;
    }
    if (nbits > 256 || !c_api_valid_pippenger_window(window) ||
        !c_api_has_valid_pair_heads(points, scalars, npoints)) {
        *ret = {};
        return;
    }
    if (npoints == 1) {
        BoundaryG1 point = {};
        BoundaryG1::from_affine(point, *points[0]);
        BoundaryG1::mult_vartime(*ret, point, scalars[0], nbits);
        return;
    }
    if (scratch == nullptr) {
        *ret = {};
        return;
    }
    blst::scalar::msm_mult_pippenger<BoundaryG1, BoundaryG1Affine, blst::scalar::Pointxyzz<blst::Fp>>(
        *ret,
        points,
        npoints,
        scalars,
        nbits,
        reinterpret_cast<blst::scalar::Pointxyzz<blst::Fp>*>(scratch),
        window);
}

void blst_p1s_tile_pippenger(
    BoundaryG1* ret,
    const BoundaryG1Affine* const points[],
    size_t npoints,
    const blst::byte* const scalars[],
    size_t nbits,
    blst::limb_t* scratch,
    size_t bit0,
    size_t window)
{
    if (npoints == 0 || !c_api_has_valid_pair_heads(points, scalars, npoints) ||
        scratch == nullptr || !c_api_valid_pippenger_window(window) ||
        nbits == 0 || nbits > 256 || bit0 >= nbits) {
        *ret = {};
        return;
    }
    c_api_points_tile_pippenger(
        *ret,
        points,
        npoints,
        scalars,
        nbits,
        reinterpret_cast<blst::scalar::Pointxyzz<blst::Fp>*>(scratch),
        bit0,
        window);
}

void blst_p1_serialize(unsigned char out[96], const BoundaryG1* in)
{
    PublicG1::serialize(out, boundary_to_public_g1(*in));
}

void blst_p1_compress(unsigned char out[48], const BoundaryG1* in)
{
    PublicG1::compress(out, boundary_to_public_g1(*in));
}

void blst_p1_affine_serialize(unsigned char out[96], const BoundaryG1Affine* in)
{
    PublicG1Affine::serialize(out, boundary_to_public_g1_affine(*in));
}

void blst_p1_affine_compress(unsigned char out[48], const BoundaryG1Affine* in)
{
    PublicG1Affine::compress(out, boundary_to_public_g1_affine(*in));
}

blst::BLST_ERROR blst_p1_uncompress(BoundaryG1Affine* out, const unsigned char in[48])
{
    PublicG1Affine affine = {};
    const blst::BLST_ERROR err = PublicG1Affine::uncompress(affine, in);
    if (err == blst::BLST_SUCCESS) {
        *out = public_to_boundary_g1_affine(affine);
    }
    return err;
}

blst::BLST_ERROR blst_p1_deserialize(BoundaryG1Affine* out, const unsigned char in[96])
{
    PublicG1Affine affine = {};
    const blst::BLST_ERROR err = PublicG1Affine::deserialize(affine, in);
    if (err == blst::BLST_SUCCESS) {
        *out = public_to_boundary_g1_affine(affine);
    }
    return err;
}

void blst_p2_add(BoundaryG2* out, const BoundaryG2* a, const BoundaryG2* b)
{
    PublicG2 result = {};
    PublicG2::add(result, boundary_to_public_g2(*a), boundary_to_public_g2(*b));
    *out = public_to_boundary_g2(result);
}

void blst_p2_sub(BoundaryG2* out, const BoundaryG2* a, const BoundaryG2* b)
{
    PublicG2 rhs = boundary_to_public_g2(*b);
    PublicG2::cneg(rhs, 1);
    PublicG2 result = {};
    PublicG2::dadd(result, boundary_to_public_g2(*a), rhs);
    *out = public_to_boundary_g2(result);
}

void blst_p2_add_or_double(BoundaryG2* out, const BoundaryG2* a, const BoundaryG2* b)
{
    PublicG2 result = {};
    PublicG2::dadd(result, boundary_to_public_g2(*a), boundary_to_public_g2(*b));
    *out = public_to_boundary_g2(result);
}

void blst_p2_add_affine(BoundaryG2* out, const BoundaryG2* a, const BoundaryG2Affine* b)
{
    PublicG2 result = {};
    PublicG2::add(result, boundary_to_public_g2(*a), boundary_to_public_g2_affine(*b));
    *out = public_to_boundary_g2(result);
}

void blst_p2_sub_affine(BoundaryG2* out, const BoundaryG2* a, const BoundaryG2Affine* b)
{
    PublicG2Affine rhs = boundary_to_public_g2_affine(*b);
    PublicG2Affine::cneg(rhs, 1);
    PublicG2 result = {};
    PublicG2::dadd(result, boundary_to_public_g2(*a), rhs);
    *out = public_to_boundary_g2(result);
}

void blst_p2_add_or_double_affine(BoundaryG2* out, const BoundaryG2* a, const BoundaryG2Affine* b)
{
    PublicG2 result = {};
    PublicG2::dadd(result, boundary_to_public_g2(*a), boundary_to_public_g2_affine(*b));
    *out = public_to_boundary_g2(result);
}

void blst_p2_double(BoundaryG2* out, const BoundaryG2* a)
{
    PublicG2 result = {};
    PublicG2::doubling(result, boundary_to_public_g2(*a));
    *out = public_to_boundary_g2(result);
}

void blst_p2_mult(BoundaryG2* out, const BoundaryG2* p, const blst::byte* scalar, size_t nbits)
{
    PublicG2 result = {};
    PublicG2::mult(result, boundary_to_public_g2(*p), scalar, nbits);
    *out = public_to_boundary_g2(result);
}

void blst_p2_mult_vartime(BoundaryG2* out, const BoundaryG2* p, const blst::byte* scalar, size_t nbits)
{
    PublicG2 result = {};
    PublicG2::mult_vartime(result, boundary_to_public_g2(*p), scalar, nbits);
    *out = public_to_boundary_g2(result);
}

void blst_p2_cneg(BoundaryG2* p, bool cbit)
{
    BoundaryG2::cneg(*p, cbit ? 1 : 0);
}

void blst_p2_affine_cneg(BoundaryG2Affine* p, bool cbit)
{
    BoundaryG2Affine::cneg(*p, cbit ? 1 : 0);
}

void blst_p2_to_affine(BoundaryG2Affine* out, const BoundaryG2* in)
{
    PublicG2Affine affine = {};
    PublicG2::to_affine(affine, boundary_to_public_g2(*in));
    *out = public_to_boundary_g2_affine(affine);
}

void blst_p2_from_affine(BoundaryG2* out, const BoundaryG2Affine* in)
{
    PublicG2 point = {};
    PublicG2::from_affine(point, boundary_to_public_g2_affine(*in));
    *out = public_to_boundary_g2(point);
}

void blst_p2_from_jacobian(BoundaryG2* out, const BoundaryG2* in)
{
    BoundaryG2::from_jacobian(*out, *in);
}

bool blst_p2_on_curve(const BoundaryG2* p)
{
    return BoundaryG2::on_curve(*p) != 0;
}

bool blst_p2_in_g2(const BoundaryG2* p)
{
    return BoundaryG2::in_subgroup(*p) != 0;
}

bool blst_p2_is_equal(const BoundaryG2* a, const BoundaryG2* b)
{
    return BoundaryG2::equal(*a, *b) != 0;
}

bool blst_p2_is_inf(const BoundaryG2* a)
{
    return BoundaryG2::is_inf(*a) != 0;
}

const BoundaryG2* blst_p2_generator(void)
{
    return &BoundaryG2::generator();
}

void blst_p2_unchecked_mult(BoundaryG2* out, const BoundaryG2* p, const blst::byte* scalar, size_t nbits)
{
	BoundaryG2 result = {};
	BoundaryG2::unchecked_mult(result, *p, scalar, nbits);
	*out = result;
}

bool blst_p2_affine_on_curve(const BoundaryG2Affine* p)
{
    return BoundaryG2Affine::on_curve(*p) != 0;
}

bool blst_p2_affine_in_g2(const BoundaryG2Affine* p)
{
    return BoundaryG2Affine::in_subgroup(*p) != 0;
}

bool blst_p2_affine_is_equal(const BoundaryG2Affine* a, const BoundaryG2Affine* b)
{
    return BoundaryG2Affine::equal(*a, *b) != 0;
}

bool blst_p2_affine_is_inf(const BoundaryG2Affine* a)
{
    return BoundaryG2Affine::is_inf(*a) != 0;
}

const BoundaryG2Affine* blst_p2_affine_generator(void)
{
    return &BoundaryG2Affine::generator();
}

void blst_p2s_to_affine(BoundaryG2Affine* dst, const BoundaryG2* const points[], size_t npoints)
{
    if (npoints == 0) {
        return;
    }
    if (!c_api_has_valid_array_head(points, npoints)) {
        for (size_t i = 0; i < npoints; ++i) {
            dst[i] = {};
        }
        return;
    }
    blst::scalar::points_to_affine<blst::Fp2>(dst, points, npoints);
}

void blst_p2s_add(BoundaryG2* ret, const BoundaryG2Affine* const points[], size_t npoints)
{
    if (npoints == 0 || !c_api_has_valid_array_head(points, npoints)) {
        *ret = {};
        return;
    }
    blst::scalar::points_add<blst::Fp2>(*ret, points, npoints);
}

size_t blst_p2s_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints)
{
    if (!c_api_valid_wbits(wbits)) {
        return 0;
    }
    return blst::scalar::points_mult_wbits_precompute_sizeof<blst::Fp2>(wbits, npoints);
}

void blst_p2s_mult_wbits_precompute(
    BoundaryG2Affine table[],
    size_t wbits,
    const BoundaryG2Affine* const points[],
    size_t npoints)
{
    if (!c_api_valid_wbits(wbits) ||
        (npoints != 0 && table == nullptr) ||
        !c_api_has_valid_array_head(points, npoints)) {
        return;
    }
    blst::scalar::msm_precompute_wbits<blst::Fp2>(table, wbits, points, npoints);
}

size_t blst_p2s_mult_wbits_scratch_sizeof(size_t npoints)
{
    return blst::scalar::points_mult_wbits_scratch_sizeof<blst::Fp2>(npoints);
}

void blst_p2s_mult_wbits(
    BoundaryG2* ret,
    const BoundaryG2Affine table[],
    size_t wbits,
    size_t npoints,
    const blst::byte* const scalars[],
    size_t nbits,
    blst::limb_t* scratch)
{
    if (npoints == 0 || nbits == 0) {
        *ret = {};
        return;
    }
    if (nbits > 256 || !c_api_valid_wbits(wbits) || table == nullptr ||
        scalars == nullptr || scalars[0] == nullptr || scratch == nullptr) {
        *ret = {};
        return;
    }

    size_t table_len = 0;
    if (!blst::scalar::detail::checked_shl_size(npoints, wbits - 1, table_len)) {
        *ret = {};
        return;
    }

    using ScratchNode = blst::scalar::detail::BatchAffineNode<blst::Fp2>;
    const size_t scratch_len = blst::scalar::points_mult_wbits_scratch_sizeof<blst::Fp2>(npoints) / sizeof(ScratchNode);
    const blst::scalar::span<const BoundaryG2Affine> table_span(table, table_len);
    blst::scalar::detail::BatchAffineWorkspace<blst::Fp2> scratch_workspace(
        scratch, scratch == nullptr ? 0 : scratch_len);
    blst::scalar::msm_mult_table<blst::Fp2>(
        *ret, table_span, wbits, npoints, scalars, nbits, scratch_workspace);
}

size_t blst_p2s_mult_pippenger_scratch_sizeof(size_t npoints)
{
    return blst::scalar::msm_mult_pippenger_scratch_sizeof_<blst::scalar::G2PointZZ>(npoints);
}

size_t blst_p2s_mult_pippenger_window_scratch_sizeof(size_t window)
{
    return blst::scalar::msm_mult_pippenger_window_scratch_sizeof_<
        blst::scalar::G2PointZZ>(window);
}

void blst_p2s_mult_pippenger(
    BoundaryG2* ret,
    const BoundaryG2Affine* const points[],
    size_t npoints,
    const blst::byte* const scalars[],
    size_t nbits,
    blst::limb_t* scratch)
{
    if (npoints == 0 || nbits == 0) {
        *ret = {};
        return;
    }
    if (nbits > 256 ||
        !c_api_has_valid_pair_heads(points, scalars, npoints)) {
        *ret = {};
        return;
    }
    if (npoints == 1) {
        BoundaryG2 point = {};
        BoundaryG2::from_affine(point, *points[0]);
        BoundaryG2::mult_vartime(*ret, point, scalars[0], nbits);
        return;
    }
    if (scratch == nullptr) {
        *ret = {};
        return;
    }
    blst::scalar::msm_mult_pippenger<BoundaryG2, BoundaryG2Affine, blst::scalar::Pointxyzz<blst::Fp2>>(
        *ret,
        points,
        npoints,
        scalars,
        nbits,
        reinterpret_cast<blst::scalar::Pointxyzz<blst::Fp2>*>(scratch),
        0);
}

void blst_p2s_mult_pippenger_window(
    BoundaryG2* ret,
    const BoundaryG2Affine* const points[],
    size_t npoints,
    const blst::byte* const scalars[],
    size_t nbits,
    blst::limb_t* scratch,
    size_t window)
{
    if (npoints == 0 || nbits == 0) {
        *ret = {};
        return;
    }
    if (nbits > 256 || !c_api_valid_pippenger_window(window) ||
        !c_api_has_valid_pair_heads(points, scalars, npoints)) {
        *ret = {};
        return;
    }
    if (npoints == 1) {
        BoundaryG2 point = {};
        BoundaryG2::from_affine(point, *points[0]);
        BoundaryG2::mult_vartime(*ret, point, scalars[0], nbits);
        return;
    }
    if (scratch == nullptr) {
        *ret = {};
        return;
    }
    blst::scalar::msm_mult_pippenger<
        BoundaryG2, BoundaryG2Affine,
        blst::scalar::Pointxyzz<blst::Fp2>>(
        *ret, points, npoints, scalars, nbits,
        reinterpret_cast<blst::scalar::Pointxyzz<blst::Fp2>*>(scratch),
        window);
}

size_t blst_p2s_mult_pippenger_fr_scratch_sizeof(size_t npoints)
{
    return c_api_pippenger_fr_scratch_sizeof<blst::scalar::G2PointZZ>(npoints);
}

bool blst_p2s_mult_pippenger_fr(
    BoundaryG2* ret,
    const BoundaryG2Affine points[],
    const blst::Fr scalars[],
    size_t npoints,
    blst::limb_t* scratch,
    size_t scratch_size)
{
    return c_api_pippenger_fr<BoundaryG2, BoundaryG2Affine, blst::scalar::G2PointZZ>(
        ret, points, scalars, npoints, scratch, scratch_size);
}

size_t blst_p2s_mult_fr_scratch_sizeof(size_t max_points)
{
    return c_api_msm_fr_scratch_sizeof<
        blst::Fp2, blst::scalar::G2PointZZ>(max_points);
}

bool blst_p2s_mult_fr(
    BoundaryG2* ret,
    const BoundaryG2Affine points[],
    const blst::Fr scalars[],
    size_t npoints,
    blst::limb_t* scratch,
    size_t scratch_size)
{
    return c_api_msm_fr<
        blst::Fp2, BoundaryG2, BoundaryG2Affine,
        blst::scalar::G2PointZZ>(
        ret, points, scalars, npoints, scratch, scratch_size);
}

size_t blst_p2s_mult_batch_affine_fr_scratch_sizeof(
    size_t npoints,
    size_t window,
    size_t batch_size)
{
    return c_api_batch_affine_fr_scratch_sizeof<blst::Fp2>(
        npoints, window, batch_size);
}

bool blst_p2s_mult_batch_affine_fr(
    BoundaryG2* ret,
    const BoundaryG2Affine points[],
    const blst::Fr scalars[],
    size_t npoints,
    blst::limb_t* scratch,
    size_t scratch_size,
    size_t window,
    size_t batch_size)
{
    return c_api_batch_affine_fr<blst::Fp2>(
        ret, points, scalars, npoints, scratch, scratch_size,
        window, batch_size);
}

void blst_p2s_tile_pippenger(
    BoundaryG2* ret,
    const BoundaryG2Affine* const points[],
    size_t npoints,
    const blst::byte* const scalars[],
    size_t nbits,
    blst::limb_t* scratch,
    size_t bit0,
    size_t window)
{
    if (npoints == 0 || !c_api_has_valid_pair_heads(points, scalars, npoints) ||
        scratch == nullptr || !c_api_valid_pippenger_window(window) ||
        nbits == 0 || nbits > 256 || bit0 >= nbits) {
        *ret = {};
        return;
    }
    c_api_points_tile_pippenger(
        *ret,
        points,
        npoints,
        scalars,
        nbits,
        reinterpret_cast<blst::scalar::Pointxyzz<blst::Fp2>*>(scratch),
        bit0,
        window);
}

void blst_p2_serialize(unsigned char out[192], const BoundaryG2* in)
{
    PublicG2::serialize(out, boundary_to_public_g2(*in));
}

void blst_p2_compress(unsigned char out[96], const BoundaryG2* in)
{
    PublicG2::compress(out, boundary_to_public_g2(*in));
}

void blst_p2_affine_serialize(unsigned char out[192], const BoundaryG2Affine* in)
{
    PublicG2Affine::serialize(out, boundary_to_public_g2_affine(*in));
}

void blst_p2_affine_compress(unsigned char out[96], const BoundaryG2Affine* in)
{
    PublicG2Affine::compress(out, boundary_to_public_g2_affine(*in));
}

blst::BLST_ERROR blst_p2_uncompress(BoundaryG2Affine* out, const unsigned char in[96])
{
    PublicG2Affine affine = {};
    const blst::BLST_ERROR err = PublicG2Affine::uncompress(affine, in);
    if (err == blst::BLST_SUCCESS) {
        *out = public_to_boundary_g2_affine(affine);
    }
    return err;
}

blst::BLST_ERROR blst_p2_deserialize(BoundaryG2Affine* out, const unsigned char in[192])
{
    PublicG2Affine affine = {};
    const blst::BLST_ERROR err = PublicG2Affine::deserialize(affine, in);
    if (err == blst::BLST_SUCCESS) {
        *out = public_to_boundary_g2_affine(affine);
    }
    return err;
}

void blst_miller_loop(BoundaryFp12* ret, const BoundaryG2Affine* Q, const BoundaryG1Affine* P)
{
    PublicMillerResult out = {};
    const PublicG2Affine q = Q == nullptr ? PublicG2Affine::generator() : boundary_to_public_g2_affine(*Q);
    const PublicG1Affine p = P == nullptr ? PublicG1Affine::generator() : boundary_to_public_g1_affine(*P);
    PublicGt::miller_loop(out, &q, &p);
    *ret = public_miller_to_boundary(out);
}

void blst_miller_loop_n(
    BoundaryFp12* ret,
    const BoundaryG2Affine* const Qs[],
    const BoundaryG1Affine* const Ps[],
    size_t n)
{
    if (n == 0) {
        *ret = public_miller_to_boundary(PublicMillerResult::identity());
        return;
    }
    if (!c_api_has_valid_pair_heads(Qs, Ps, n)) {
        *ret = {};
        return;
    }

#if BLST_DEFAULT_AVX
    constexpr size_t batch_capacity = blst::avx::avx_gt_detail::kMillerLoopNMax;
    PublicG2Affine q_values[batch_capacity] = {};
    PublicG1Affine p_values[batch_capacity] = {};
    PublicMillerResult out = PublicMillerResult::identity();
    const BoundaryG2Affine* current_q = nullptr;
    const BoundaryG1Affine* current_p = nullptr;
    const BoundaryG2Affine* const* q_iter = Qs;
    const BoundaryG1Affine* const* p_iter = Ps;
    size_t batch_size = 0;

    for (size_t i = 0; i < n; ++i) {
        current_q = *q_iter != nullptr ? *q_iter++ : current_q + 1;
        current_p = *p_iter != nullptr ? *p_iter++ : current_p + 1;
        q_values[batch_size] = boundary_to_public_g2_affine(*current_q);
        p_values[batch_size] = boundary_to_public_g1_affine(*current_p);
        ++batch_size;

        if (batch_size == batch_capacity || i + 1 == n) {
            PublicMillerResult batch = {};
            PublicGt::miller_loop_n(batch, q_values, p_values, batch_size);
            PublicMillerResult::mul(out, out, batch);
            batch_size = 0;
        }
    }
#else
    static_assert(std::is_same_v<PublicG2Affine, BoundaryG2Affine>);
    static_assert(std::is_same_v<PublicG1Affine, BoundaryG1Affine>);
    PublicMillerResult out = {};
    PublicGt::miller_loop_n(out, Qs, Ps, n);
#endif
    *ret = public_miller_to_boundary(out);
}

void blst_precompute_lines(BoundaryFp6 Qlines[68], const BoundaryG2Affine* Q)
{
    const PublicG2Affine q = Q == nullptr ? PublicG2Affine::generator() : boundary_to_public_g2_affine(*Q);
    PublicGt::precompute_lines(Qlines, q);
}

void blst_miller_loop_lines(BoundaryFp12* ret, const BoundaryFp6 Qlines[68], const BoundaryG1Affine* P)
{
    PublicMillerResult out = {};
    const PublicG1Affine p = P == nullptr ? PublicG1Affine::generator() : boundary_to_public_g1_affine(*P);
    PublicGt::miller_loop_lines(out, Qlines, p);
    *ret = public_miller_to_boundary(out);
}

size_t blst_pairing_sizeof(void)
{
    return blst::blst_pairing_sizeof();
}

size_t blst_pairing_alignof(void)
{
    return alignof(BoundaryPairing);
}

void blst_pairing_init(blst_pairing* new_ctx, bool hash_or_encode, const byte* DST, size_t DST_len)
{
    blst::blst_pairing_init(as_pairing(new_ctx), hash_or_encode ? 1 : 0, DST, DST_len);
}

const byte* blst_pairing_get_dst(const blst_pairing* ctx)
{
    return reinterpret_cast<const unsigned char*>(blst::blst_pairing_get_dst(as_pairing(ctx)));
}

void blst_pairing_commit(blst_pairing* ctx)
{
    blst::blst_pairing_commit(as_pairing(ctx));
}

blst::BLST_ERROR blst_pairing_aggregate_pk_in_g2(
    blst_pairing* ctx,
    const BoundaryG2Affine* PK,
    const BoundaryG1Affine* signature,
    const unsigned char* msg,
    size_t msg_len,
    const unsigned char* aug,
    size_t aug_len)
{
    PublicG2Affine pk = {};
    PublicG1Affine sig = {};
    const PublicG2Affine* pk_ptr = PK == nullptr ? nullptr : &(pk = boundary_to_public_g2_affine(*PK));
    const PublicG1Affine* sig_ptr = signature == nullptr ? nullptr : &(sig = boundary_to_public_g1_affine(*signature));
    return blst::blst_pairing_aggregate_pk_in_g2(as_pairing(ctx), pk_ptr, sig_ptr, msg, msg_len, aug, aug_len);
}

blst::BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g2(
    blst_pairing* ctx,
    const BoundaryG2Affine* PK,
    bool pk_grpchk,
    const BoundaryG1Affine* signature,
    bool sig_grpchk,
    const unsigned char* msg,
    size_t msg_len,
    const unsigned char* aug,
    size_t aug_len)
{
    PublicG2Affine pk = {};
    PublicG1Affine sig = {};
    const PublicG2Affine* pk_ptr = PK == nullptr ? nullptr : &(pk = boundary_to_public_g2_affine(*PK));
    const PublicG1Affine* sig_ptr = signature == nullptr ? nullptr : &(sig = boundary_to_public_g1_affine(*signature));
    return blst::blst_pairing_chk_n_aggr_pk_in_g2(
        as_pairing(ctx),
        pk_ptr,
        pk_grpchk ? 1 : 0,
        sig_ptr,
        sig_grpchk ? 1 : 0,
        msg,
        msg_len,
        aug,
        aug_len);
}

blst::BLST_ERROR blst_pairing_mul_n_aggregate_pk_in_g2(
    blst_pairing* ctx,
    const BoundaryG2Affine* PK,
    const BoundaryG1Affine* sig,
    const unsigned char* scalar,
    size_t nbits,
    const unsigned char* msg,
    size_t msg_len,
    const unsigned char* aug,
    size_t aug_len)
{
    PublicG2Affine pk = {};
    PublicG1Affine signature = {};
    const PublicG2Affine* pk_ptr = PK == nullptr ? nullptr : &(pk = boundary_to_public_g2_affine(*PK));
    const PublicG1Affine* sig_ptr = sig == nullptr ? nullptr : &(signature = boundary_to_public_g1_affine(*sig));
    return blst::blst_pairing_mul_n_aggregate_pk_in_g2(
        as_pairing(ctx),
        pk_ptr,
        sig_ptr,
        scalar,
        nbits,
        msg,
        msg_len,
        aug,
        aug_len);
}

blst::BLST_ERROR blst_pairing_chk_n_mul_n_aggr_pk_in_g2(
    blst_pairing* ctx,
    const BoundaryG2Affine* PK,
    bool pk_grpchk,
    const BoundaryG1Affine* sig,
    bool sig_grpchk,
    const unsigned char* scalar,
    size_t nbits,
    const unsigned char* msg,
    size_t msg_len,
    const unsigned char* aug,
    size_t aug_len)
{
    PublicG2Affine pk = {};
    PublicG1Affine signature = {};
    const PublicG2Affine* pk_ptr = PK == nullptr ? nullptr : &(pk = boundary_to_public_g2_affine(*PK));
    const PublicG1Affine* sig_ptr = sig == nullptr ? nullptr : &(signature = boundary_to_public_g1_affine(*sig));
    return blst::blst_pairing_chk_n_mul_n_aggr_pk_in_g2(
        as_pairing(ctx),
        pk_ptr,
        pk_grpchk ? 1 : 0,
        sig_ptr,
        sig_grpchk ? 1 : 0,
        scalar,
        nbits,
        msg,
        msg_len,
        aug,
        aug_len);
}

blst::BLST_ERROR blst_pairing_aggregate_pk_in_g1(
    blst_pairing* ctx,
    const BoundaryG1Affine* PK,
    const BoundaryG2Affine* signature,
    const unsigned char* msg,
    size_t msg_len,
    const unsigned char* aug,
    size_t aug_len)
{
    PublicG1Affine pk = {};
    PublicG2Affine sig = {};
    const PublicG1Affine* pk_ptr = PK == nullptr ? nullptr : &(pk = boundary_to_public_g1_affine(*PK));
    const PublicG2Affine* sig_ptr = signature == nullptr ? nullptr : &(sig = boundary_to_public_g2_affine(*signature));
    return blst::blst_pairing_aggregate_pk_in_g1(as_pairing(ctx), pk_ptr, sig_ptr, msg, msg_len, aug, aug_len);
}

blst::BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g1(
    blst_pairing* ctx,
    const BoundaryG1Affine* PK,
    bool pk_grpchk,
    const BoundaryG2Affine* signature,
    bool sig_grpchk,
    const unsigned char* msg,
    size_t msg_len,
    const unsigned char* aug,
    size_t aug_len)
{
    PublicG1Affine pk = {};
    PublicG2Affine sig = {};
    const PublicG1Affine* pk_ptr = PK == nullptr ? nullptr : &(pk = boundary_to_public_g1_affine(*PK));
    const PublicG2Affine* sig_ptr = signature == nullptr ? nullptr : &(sig = boundary_to_public_g2_affine(*signature));
    return blst::blst_pairing_chk_n_aggr_pk_in_g1(
        as_pairing(ctx),
        pk_ptr,
        pk_grpchk ? 1 : 0,
        sig_ptr,
        sig_grpchk ? 1 : 0,
        msg,
        msg_len,
        aug,
        aug_len);
}

blst::BLST_ERROR blst_pairing_mul_n_aggregate_pk_in_g1(
    blst_pairing* ctx,
    const BoundaryG1Affine* PK,
    const BoundaryG2Affine* sig,
    const unsigned char* scalar,
    size_t nbits,
    const unsigned char* msg,
    size_t msg_len,
    const unsigned char* aug,
    size_t aug_len)
{
    PublicG1Affine pk = {};
    PublicG2Affine signature = {};
    const PublicG1Affine* pk_ptr = PK == nullptr ? nullptr : &(pk = boundary_to_public_g1_affine(*PK));
    const PublicG2Affine* sig_ptr = sig == nullptr ? nullptr : &(signature = boundary_to_public_g2_affine(*sig));
    return blst::blst_pairing_mul_n_aggregate_pk_in_g1(
        as_pairing(ctx),
        pk_ptr,
        sig_ptr,
        scalar,
        nbits,
        msg,
        msg_len,
        aug,
        aug_len);
}

blst::BLST_ERROR blst_pairing_chk_n_mul_n_aggr_pk_in_g1(
    blst_pairing* ctx,
    const BoundaryG1Affine* PK,
    bool pk_grpchk,
    const BoundaryG2Affine* sig,
    bool sig_grpchk,
    const unsigned char* scalar,
    size_t nbits,
    const unsigned char* msg,
    size_t msg_len,
    const unsigned char* aug,
    size_t aug_len)
{
    PublicG1Affine pk = {};
    PublicG2Affine signature = {};
    const PublicG1Affine* pk_ptr = PK == nullptr ? nullptr : &(pk = boundary_to_public_g1_affine(*PK));
    const PublicG2Affine* sig_ptr = sig == nullptr ? nullptr : &(signature = boundary_to_public_g2_affine(*sig));
    return blst::blst_pairing_chk_n_mul_n_aggr_pk_in_g1(
        as_pairing(ctx),
        pk_ptr,
        pk_grpchk ? 1 : 0,
        sig_ptr,
        sig_grpchk ? 1 : 0,
        scalar,
        nbits,
        msg,
        msg_len,
        aug,
        aug_len);
}

blst::BLST_ERROR blst_pairing_merge(blst_pairing* ctx, const blst_pairing* ctx1)
{
    return blst::blst_pairing_merge(as_pairing(ctx), as_pairing(ctx1));
}

bool blst_pairing_finalverify(const blst_pairing* ctx, const BoundaryFp12* gtsig)
{
    return blst::blst_pairing_finalverify(as_pairing(ctx), gtsig) != 0;
}

blst::BLST_ERROR blst_aggregate_in_g1(BoundaryG1* out, const BoundaryG1* in, const unsigned char* zwire)
{
    PublicG1 public_out = {};
    const PublicG1 public_in = in == nullptr ? PublicG1::identity() : boundary_to_public_g1(*in);
    const blst::BLST_ERROR err = blst::blst_aggregate_in_g1(public_out, in == nullptr ? nullptr : &public_in, zwire);
    if (err == blst::BLST_SUCCESS) {
        *out = public_to_boundary_g1(public_out);
    }
    return err;
}

blst::BLST_ERROR blst_aggregate_in_g2(BoundaryG2* out, const BoundaryG2* in, const unsigned char* zwire)
{
    PublicG2 public_out = {};
    const PublicG2 public_in = in == nullptr ? PublicG2::identity() : boundary_to_public_g2(*in);
    const blst::BLST_ERROR err = blst::blst_aggregate_in_g2(public_out, in == nullptr ? nullptr : &public_in, zwire);
    if (err == blst::BLST_SUCCESS) {
        *out = public_to_boundary_g2(public_out);
    }
    return err;
}

void blst_aggregated_in_g1(BoundaryFp12* out, const BoundaryG1Affine* signature)
{
    PublicMillerResult gt = {};
    const PublicG1Affine sig = boundary_to_public_g1_affine(*signature);
    blst::blst_aggregated_in_g1(gt, sig);
    *out = public_miller_to_boundary(gt);
}

void blst_aggregated_in_g2(BoundaryFp12* out, const BoundaryG2Affine* signature)
{
    PublicMillerResult gt = {};
    const PublicG2Affine sig = boundary_to_public_g2_affine(*signature);
    blst::blst_aggregated_in_g2(gt, sig);
    *out = public_miller_to_boundary(gt);
}

blst::BLST_ERROR blst_core_verify_pk_in_g1(
    const BoundaryG1Affine* pk,
    const BoundaryG2Affine* signature,
    bool hash_or_encode,
    const unsigned char* msg,
    size_t msg_len,
    const unsigned char* DST,
    size_t DST_len,
    const unsigned char* aug,
    size_t aug_len)
{
    const PublicG1Affine public_pk = boundary_to_public_g1_affine(*pk);
    const PublicG2Affine public_sig = boundary_to_public_g2_affine(*signature);
    return blst::blst_core_verify_pk_in_g1(
        &public_pk,
        &public_sig,
        hash_or_encode ? 1 : 0,
        msg,
        msg_len,
        DST,
        DST_len,
        aug,
        aug_len);
}

blst::BLST_ERROR blst_core_verify_pk_in_g2(
    const BoundaryG2Affine* pk,
    const BoundaryG1Affine* signature,
    bool hash_or_encode,
    const unsigned char* msg,
    size_t msg_len,
    const unsigned char* DST,
    size_t DST_len,
    const unsigned char* aug,
    size_t aug_len)
{
    const PublicG2Affine public_pk = boundary_to_public_g2_affine(*pk);
    const PublicG1Affine public_sig = boundary_to_public_g1_affine(*signature);
    return blst::blst_core_verify_pk_in_g2(
        &public_pk,
        &public_sig,
        hash_or_encode ? 1 : 0,
        msg,
        msg_len,
        DST,
        DST_len,
        aug,
        aug_len);
}

void blst_pairing_raw_aggregate(blst_pairing* ctx, const BoundaryG2Affine* q, const BoundaryG1Affine* p)
{
    const PublicG2Affine public_q = boundary_to_public_g2_affine(*q);
    const PublicG1Affine public_p = boundary_to_public_g1_affine(*p);
    blst::blst_pairing_raw_aggregate(as_pairing(ctx), &public_q, &public_p);
}

void blst_pairing_get_miller_result(BoundaryFp12* out, blst_pairing* ctx)
{
    blst::blst_pairing_get_miller_result(out, as_pairing(ctx));
}

void blst_keygen(blst_scalar* out_SK, const byte* IKM, size_t IKM_len, const byte* info, size_t info_len)
{
    BLST_ASSERT(out_SK != nullptr);
    BLST_ASSERT(IKM != nullptr || IKM_len < 32);
    blst::blst_keygen(out_SK->b, IKM, IKM_len, info, info_len);
}

void blst_keygen_v3(blst_scalar* out_SK, const byte* IKM, size_t IKM_len, const byte* info, size_t info_len)
{
    BLST_ASSERT(out_SK != nullptr);
    BLST_ASSERT(IKM != nullptr || IKM_len < 32);
    blst::blst_keygen_v3(out_SK->b, IKM, IKM_len, info, info_len);
}

void blst_keygen_v4_5(
    blst_scalar* out_SK,
    const byte* IKM,
    size_t IKM_len,
    const byte* salt,
    size_t salt_len,
    const byte* info,
    size_t info_len)
{
    BLST_ASSERT(out_SK != nullptr);
    BLST_ASSERT(IKM != nullptr || IKM_len < 32);
    blst::blst_keygen_v4_5(out_SK->b, IKM, IKM_len, salt, salt_len, info, info_len);
}

void blst_keygen_v5(
    blst_scalar* out_SK,
    const byte* IKM,
    size_t IKM_len,
    const byte* salt,
    size_t salt_len,
    const byte* info,
    size_t info_len)
{
    BLST_ASSERT(out_SK != nullptr);
    BLST_ASSERT(IKM != nullptr || IKM_len < 32);
    blst::blst_keygen_v5(out_SK->b, IKM, IKM_len, salt, salt_len, info, info_len);
}

void blst_derive_master_eip2333(blst_scalar* out_SK, const byte* IKM, size_t IKM_len)
{
    BLST_ASSERT(out_SK != nullptr);
    BLST_ASSERT(IKM != nullptr || IKM_len < 32);
    blst::blst_derive_master_eip2333(out_SK->b, IKM, IKM_len);
}

void blst_derive_child_eip2333(blst_scalar* out_SK, const blst_scalar* SK, uint32_t child_index)
{
    BLST_ASSERT(out_SK != nullptr);
    BLST_ASSERT(SK != nullptr);
    blst::blst_derive_child_eip2333(out_SK->b, SK->b, child_index);
}

size_t blst_uniq_sizeof(size_t capacity)
{
    return blst::blst_uniq_sizeof(capacity);
}

size_t blst_uniq_alignof(void)
{
    return blst::blst_uniq_alignof();
}

bool blst_uniq_init(blst_uniq* tree, size_t capacity)
{
    return blst::blst_uniq_init(tree, capacity);
}

BLST_UNIQ_RESULT blst_uniq_insert(blst_uniq* tree, const byte* msg, size_t len)
{
    return static_cast<BLST_UNIQ_RESULT>(blst::blst_uniq_insert(tree, msg, len));
}

void blst_expand_message_xmd(blst::byte* out, size_t out_len, const blst::byte* msg, size_t msg_len, const blst::byte* DST, size_t DST_len)
{
    blst::scalar::blst_expand_message_xmd(out, out_len, msg, msg_len, DST, DST_len);
}

size_t blst_p1_sizeof(void)
{
    return sizeof(BoundaryG1);
}

size_t blst_p1_affine_sizeof(void)
{
    return sizeof(BoundaryG1Affine);
}

size_t blst_p2_sizeof(void)
{
    return sizeof(BoundaryG2);
}

size_t blst_p2_affine_sizeof(void)
{
    return sizeof(BoundaryG2Affine);
}

size_t blst_fp12_sizeof(void)
{
    return sizeof(BoundaryFp12);
}

void blst_encode_to_g1(
    BoundaryG1* out,
    const blst::byte* msg,
    size_t msg_len,
    const blst::byte* DST,
    size_t DST_len,
    const blst::byte* aug,
    size_t aug_len)
{
    BoundaryG1 scalar_out = {};
    blst::scalar::Encode_to_curve(scalar_out, msg, msg_len, DST, DST_len, aug, aug_len);
    *out = scalar_out;
}

void blst_hash_to_g1(
    BoundaryG1* out,
    const blst::byte* msg,
    size_t msg_len,
    const blst::byte* DST,
    size_t DST_len,
    const blst::byte* aug,
    size_t aug_len)
{
    BoundaryG1 scalar_out = {};
    blst::scalar::Hash_to_curve(scalar_out, msg, msg_len, DST, DST_len, aug, aug_len);
    *out = scalar_out;
}

void blst_encode_to_g2(
    BoundaryG2* out,
    const blst::byte* msg,
    size_t msg_len,
    const blst::byte* DST,
    size_t DST_len,
    const blst::byte* aug,
    size_t aug_len)
{
    BoundaryG2 scalar_out = {};
    blst::scalar::Encode_to_curve(scalar_out, msg, msg_len, DST, DST_len, aug, aug_len);
    *out = scalar_out;
}

void blst_hash_to_g2(
    BoundaryG2* out,
    const blst::byte* msg,
    size_t msg_len,
    const blst::byte* DST,
    size_t DST_len,
    const blst::byte* aug,
    size_t aug_len)
{
    BoundaryG2 scalar_out = {};
    blst::scalar::Hash_to_curve(scalar_out, msg, msg_len, DST, DST_len, aug, aug_len);
    *out = scalar_out;
}

void blst_map_to_g1(BoundaryG1* out, const blst::Fp* u, const blst::Fp* v)
{
    blst::scalar::map_to_curve(*out, *u, v);
}

void blst_map_to_g2(BoundaryG2* out, const blst::Fp2* u, const blst::Fp2* v)
{
    blst::scalar::map_to_curve(*out, *u, v);
}

void blst_sk_to_pk_in_g1(BoundaryG1* out_pk, const blst_scalar* SK)
{
    PublicG1 public_out = {};
    PublicG1::to_pk(public_out, SK->b);
    *out_pk = public_to_boundary_g1(public_out);
}

void blst_sign_pk_in_g1(BoundaryG2* out_sig, const BoundaryG2* hash, const blst_scalar* SK)
{
    PublicG2 public_out = {};
    PublicG2::sign(public_out, boundary_to_public_g2(*hash), SK->b);
    *out_sig = public_to_boundary_g2(public_out);
}

void blst_sk_to_pk_in_g2(BoundaryG2* out_pk, const blst_scalar* SK)
{
    PublicG2 public_out = {};
    PublicG2::to_pk(public_out, SK->b);
    *out_pk = public_to_boundary_g2(public_out);
}

void blst_sign_pk_in_g2(BoundaryG1* out_sig, const BoundaryG1* hash, const blst_scalar* SK)
{
    PublicG1 public_out = {};
    PublicG1::sign(public_out, boundary_to_public_g1(*hash), SK->b);
    *out_sig = public_to_boundary_g1(public_out);
}

void blst_sk_to_pk2_in_g1(byte out[96], BoundaryG1Affine* out_pk, const blst_scalar* SK)
{
    PublicG1Affine public_pk = {};
    PublicG1Affine* public_pk_ptr = out_pk == nullptr ? nullptr : &public_pk;
    PublicG1::to_pk2(out, public_pk_ptr, SK->b);
    if (out_pk != nullptr) {
        *out_pk = public_to_boundary_g1_affine(public_pk);
    }
}

void blst_sign_pk2_in_g1(byte out[192], BoundaryG2Affine* out_sig, const BoundaryG2* hash, const blst_scalar* SK)
{
    PublicG2 public_point = {};
    PublicG2::sign(public_point, boundary_to_public_g2(*hash), SK->b);
    if (out_sig != nullptr) {
        PublicG2Affine public_sig = {};
        PublicG2::to_affine(public_sig, public_point);
        *out_sig = public_to_boundary_g2_affine(public_sig);
    }
    if (out != nullptr) {
        PublicG2::serialize(out, public_point);
        PublicG2Affine affine = {};
        PublicG2::to_affine(affine, public_point);
        out[0] |= (BoundaryG2Affine::compress_BE(out, public_to_boundary_g2_affine(affine)) & 2) << 4;
        out[0] |= PublicG2::is_inf(public_point) << 6;
    }
}

void blst_sk_to_pk2_in_g2(byte out[192], BoundaryG2Affine* out_pk, const blst_scalar* SK)
{
    PublicG2Affine public_pk = {};
    PublicG2Affine* public_pk_ptr = out_pk == nullptr ? nullptr : &public_pk;
    PublicG2::to_pk2(out, public_pk_ptr, SK->b);
    if (out_pk != nullptr) {
        *out_pk = public_to_boundary_g2_affine(public_pk);
    }
}

void blst_sign_pk2_in_g2(byte out[96], BoundaryG1Affine* out_sig, const BoundaryG1* hash, const blst_scalar* SK)
{
    PublicG1 public_point = {};
    PublicG1::sign(public_point, boundary_to_public_g1(*hash), SK->b);
    if (out_sig != nullptr) {
        PublicG1Affine public_sig = {};
        PublicG1::to_affine(public_sig, public_point);
        *out_sig = public_to_boundary_g1_affine(public_sig);
    }
    if (out != nullptr) {
        PublicG1::serialize(out, public_point);
        PublicG1Affine affine = {};
        PublicG1::to_affine(affine, public_point);
        out[0] |= (BoundaryG1Affine::compress_BE(out, public_to_boundary_g1_affine(affine)) & 2) << 4;
        out[0] |= PublicG1::is_inf(public_point) << 6;
    }
}

BLST_EXTERN_C_END_
