#include "blst/backend/backend.h"
#include "blst/common/vect.h"
#include "blst/curve/GT.h"
#include "blst/curve/map_to_G.h"
#include "blst/curve/point.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/hash_to_field.h"
#include "blst/sign/aggregate.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/multi_scalar.h"
#include "blst/scalar/curve/point.h"
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/map_to_G.h"
#endif

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace blst {
struct rb_tree;
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
void blst_aggregated_in_g1(GT& ret, const G1AffinePoint& sig);
void blst_aggregated_in_g2(GT& ret, const G2AffinePoint& sig);
BLST_ERROR blst_core_verify_pk_in_g1(const G1AffinePoint* pk, const G2AffinePoint* signature, int hash_or_encode, const void* msg, size_t msg_len, const void* DST, size_t DST_len, const void* aug, size_t aug_len);
BLST_ERROR blst_core_verify_pk_in_g2(const G2AffinePoint* pk, const G1AffinePoint* signature, int hash_or_encode, const void* msg, size_t msg_len, const void* DST, size_t DST_len, const void* aug, size_t aug_len);
void blst_pairing_raw_aggregate(PAIRING* ctx, const G2AffinePoint* q, const G1AffinePoint* p);

void blst_keygen(pow256 SK, const void* IKM, size_t IKM_len, const void* info, size_t info_len);
void blst_keygen_v3(pow256 SK, const void* IKM, size_t IKM_len, const void* info, size_t info_len);
void blst_keygen_v4_5(pow256 SK, const void* IKM, size_t IKM_len, const void* salt, size_t salt_len, const void* info, size_t info_len);
void blst_keygen_v5(pow256 SK, const void* IKM, size_t IKM_len, const void* salt, size_t salt_len, const void* info, size_t info_len);
void blst_derive_master_eip2333(pow256 SK, const void* seed, size_t seed_len);
void blst_derive_child_eip2333(pow256 SK, const pow256 parent_SK, unsigned int child_index);
size_t blst_uniq_sizeof(size_t n_nodes);
void blst_uniq_init(rb_tree* tree);
int blst_uniq_test(rb_tree* tree, const void* data, size_t len);
}

namespace blst::scalar {
template<typename Field>
void points_to_affine(AffinePoint<Field> dst[], const Point<Field>* const points[], size_t npoints);
}

namespace {

using BoundaryG1 = blst::scalar::G1Point;
using BoundaryG1Affine = blst::scalar::G1AffinePoint;
using BoundaryG2 = blst::scalar::G2Point;
using BoundaryG2Affine = blst::scalar::G2AffinePoint;
using BoundaryFp12 = blst::Fp12;
using BoundaryFp6 = blst::Fp6;

using PublicG1 = blst::G1Point;
using PublicG1Affine = blst::G1AffinePoint;
using PublicG2 = blst::G2Point;
using PublicG2Affine = blst::G2AffinePoint;
using PublicGt = blst::GT;
using BoundaryPairing = blst::PAIRING;

struct BlstP1Layout
{
    blst::Fp x;
    blst::Fp y;
    blst::Fp z;
};

struct BlstP1AffineLayout
{
    blst::Fp x;
    blst::Fp y;
};

struct BlstP2Layout
{
    blst::Fp2 x;
    blst::Fp2 y;
    blst::Fp2 z;
};

struct BlstP2AffineLayout
{
    blst::Fp2 x;
    blst::Fp2 y;
};

struct BlstFp12Layout
{
    blst::Fp6 fp6[2];
};

static_assert(sizeof(BlstP1Layout) == sizeof(BoundaryG1));
static_assert(alignof(BlstP1Layout) == alignof(BoundaryG1));
static_assert(sizeof(BlstP1AffineLayout) == sizeof(BoundaryG1Affine));
static_assert(alignof(BlstP1AffineLayout) == alignof(BoundaryG1Affine));
static_assert(sizeof(BlstP2Layout) == sizeof(BoundaryG2));
static_assert(alignof(BlstP2Layout) == alignof(BoundaryG2));
static_assert(sizeof(BlstP2AffineLayout) == sizeof(BoundaryG2Affine));
static_assert(alignof(BlstP2AffineLayout) == alignof(BoundaryG2Affine));
static_assert(sizeof(BlstFp12Layout) == sizeof(BoundaryFp12));
static_assert(alignof(BlstFp12Layout) == alignof(BoundaryFp12));
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
    unsigned char serialized[12 * 48] = {};
    blst::scalar::GT::serialize(serialized, blst::scalar::GT(value));
    PublicGt out = {};
    PublicGt::deserialize(out, serialized);
    return out;
}

inline BoundaryFp12 public_to_boundary_gt(const PublicGt& value)
{
    unsigned char serialized[12 * 48] = {};
    PublicGt::serialize(serialized, value);
    blst::scalar::GT scalar_gt = {};
    blst::scalar::GT::deserialize(scalar_gt, serialized);
    return scalar_gt.data;
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
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.X),
        blst::avx::avx_curve_detail::scalar_fp_to_ring(point.Y));
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
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.X),
        blst::avx::avx_curve_detail::scalar_fp2_to_ring(point.Y));
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

inline BoundaryPairing* as_pairing(void* ctx)
{
    return reinterpret_cast<BoundaryPairing*>(ctx);
}

inline const BoundaryPairing* as_pairing(const void* ctx)
{
    return reinterpret_cast<const BoundaryPairing*>(ctx);
}

inline blst::rb_tree* as_uniq(void* tree)
{
    return reinterpret_cast<blst::rb_tree*>(tree);
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
    if (npoints == 0 || !c_api_has_valid_pair_heads(points, scalars, npoints)) {
        return;
    }

    BLST_ASSERT(npoints > 1);

    const size_t window_idx = bit0 / window;
    const blst::scalar::WindowExtract we(window_idx, window, nbits);
    size_t i = 0;
    const blst::byte* scalar = *scalars++;
    const AffineGroup* point = *points++;
    const size_t nbytes = (nbits + 7) / 8;

    auto extract_window = [&]() -> std::pair<const AffineGroup*, blst::limb_t> {
        const blst::limb_t encoding = we.extract(scalar);
        const AffineGroup* current = point;
        scalar = *scalars ? *scalars++ : scalar + nbytes;
        point = *points ? *points++ : point + 1;
        return {current, encoding};
    };

    auto next = extract_window();
    ++i;

    for (; i < npoints; ++i) {
        auto current = next;
        next = extract_window();
        c_api_point_prefetch(buckets, next.second, we.encode_size);
        c_api_point_bucket(buckets, current.second, we.encode_size, *current.first);
    }

    c_api_point_bucket(buckets, next.second, we.encode_size, *next.first);
    c_api_point_integrate_buckets(ret, buckets, we.encode_size - 1);
}

template<typename Group, typename AffineGroup>
inline void c_api_points_add(Group& out, const AffineGroup* const points[], size_t npoints)
{
    out = {};
    if (npoints == 0 || !c_api_has_valid_array_head(points, npoints)) {
        return;
    }

    const AffineGroup* point = nullptr;
    for (size_t i = 0; i < npoints; ++i) {
        point = *points ? *points++ : point + 1;
        Group::dadd(out, out, *point);
    }
}

} // namespace

BLST_EXTERN_C_BEGIN_

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

void blst_fp12_frobenius_map(BoundaryFp12* ret, const BoundaryFp12* a, size_t n)
{
    BoundaryFp12::frobenius_map(*ret, *a, n);
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
    PublicGt::final_mult(out, boundary_to_public_gt(*f));
    *ret = public_to_boundary_gt(out);
}

void blst_bendian_from_fp12(unsigned char out[48 * 12], const BoundaryFp12* a)
{
    BoundaryFp12::to_bendian(out, *a);
}

bool blst_fp12_finalverify(const BoundaryFp12* gt1, const BoundaryFp12* gt2)
{
    PublicGt lhs = boundary_to_public_gt(*gt1);
    const PublicGt rhs = boundary_to_public_gt(*gt2);
    PublicGt::conjugate(lhs);
    PublicGt::add(lhs, lhs, rhs);
    PublicGt::final_mult(lhs, lhs);
    return PublicGt::is_identity(lhs) != 0;
}

void blst_p1_add(BoundaryG1* out, const BoundaryG1* a, const BoundaryG1* b)
{
    PublicG1 result = {};
    PublicG1::add(result, boundary_to_public_g1(*a), boundary_to_public_g1(*b));
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

void blst_p1_cneg(BoundaryG1* p, bool cbit)
{
    BoundaryG1::cneg(*p, cbit ? 1 : 0);
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
    blst_p1_mult(out, p, scalar, nbits);
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
    c_api_points_add(*ret, points, npoints);
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
    if (!c_api_valid_wbits(wbits) || !c_api_has_valid_array_head(points, npoints)) {
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
    if (npoints == 0 || !c_api_valid_wbits(wbits) || scalars == nullptr || scalars[0] == nullptr) {
        *ret = {};
        return;
    }

    using ScratchPoint = blst::scalar::Point<blst::Fp>;
    const size_t scratch_len = blst::scalar::points_mult_wbits_scratch_sizeof<blst::Fp>(npoints) / sizeof(ScratchPoint);
    const blst::scalar::span<const BoundaryG1Affine> table_span(table, npoints << (wbits - 1));
    blst::scalar::span<ScratchPoint> scratch_span(
        reinterpret_cast<ScratchPoint*>(scratch),
        scratch == nullptr ? 0 : scratch_len);
    blst::scalar::msm_mult_table<blst::Fp>(*ret, table_span, wbits, npoints, scalars, nbits, scratch_span);
}

size_t blst_p1s_mult_pippenger_scratch_sizeof(size_t npoints)
{
    return sizeof(blst::scalar::G1PointZZ) << (blst::scalar::pippenger_window_size(npoints) - 1);
}

void blst_p1s_mult_pippenger(
    BoundaryG1* ret,
    const BoundaryG1Affine* const points[],
    size_t npoints,
    const blst::byte* const scalars[],
    size_t nbits,
    blst::limb_t* scratch)
{
    if (npoints == 0 || !c_api_has_valid_pair_heads(points, scalars, npoints)) {
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
    if (npoints == 0 || !c_api_has_valid_pair_heads(points, scalars, npoints)) {
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

void blst_p2_cneg(BoundaryG2* p, bool cbit)
{
    BoundaryG2::cneg(*p, cbit ? 1 : 0);
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
    blst_p2_mult(out, p, scalar, nbits);
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
    c_api_points_add(*ret, points, npoints);
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
    if (!c_api_valid_wbits(wbits) || !c_api_has_valid_array_head(points, npoints)) {
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
    if (npoints == 0 || !c_api_valid_wbits(wbits) || scalars == nullptr || scalars[0] == nullptr) {
        *ret = {};
        return;
    }

    using ScratchPoint = blst::scalar::Point<blst::Fp2>;
    const size_t scratch_len = blst::scalar::points_mult_wbits_scratch_sizeof<blst::Fp2>(npoints) / sizeof(ScratchPoint);
    const blst::scalar::span<const BoundaryG2Affine> table_span(table, npoints << (wbits - 1));
    blst::scalar::span<ScratchPoint> scratch_span(
        reinterpret_cast<ScratchPoint*>(scratch),
        scratch == nullptr ? 0 : scratch_len);
    blst::scalar::msm_mult_table<blst::Fp2>(*ret, table_span, wbits, npoints, scalars, nbits, scratch_span);
}

size_t blst_p2s_mult_pippenger_scratch_sizeof(size_t npoints)
{
    return sizeof(blst::scalar::G2PointZZ) << (blst::scalar::pippenger_window_size(npoints) - 1);
}

void blst_p2s_mult_pippenger(
    BoundaryG2* ret,
    const BoundaryG2Affine* const points[],
    size_t npoints,
    const blst::byte* const scalars[],
    size_t nbits,
    blst::limb_t* scratch)
{
    if (npoints == 0 || !c_api_has_valid_pair_heads(points, scalars, npoints)) {
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
    if (npoints == 0 || !c_api_has_valid_pair_heads(points, scalars, npoints)) {
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
    PublicGt out = {};
    const PublicG2Affine q = Q == nullptr ? PublicG2Affine::generator() : boundary_to_public_g2_affine(*Q);
    const PublicG1Affine p = P == nullptr ? PublicG1Affine::generator() : boundary_to_public_g1_affine(*P);
    PublicGt::miller_loop(out, &q, &p);
    *ret = public_to_boundary_gt(out);
}

void blst_miller_loop_n(
    BoundaryFp12* ret,
    const BoundaryG2Affine* const Qs[],
    const BoundaryG1Affine* const Ps[],
    size_t n)
{
    if (n == 0 || !c_api_has_valid_pair_heads(Qs, Ps, n)) {
        *ret = public_to_boundary_gt(PublicGt::identity());
        return;
    }

    std::vector<PublicG2Affine> q_values(n);
    std::vector<PublicG1Affine> p_values(n);
    const BoundaryG2Affine* current_q = nullptr;
    const BoundaryG1Affine* current_p = nullptr;
    const BoundaryG2Affine* const* q_iter = Qs;
    const BoundaryG1Affine* const* p_iter = Ps;

    for (size_t i = 0; i < n; ++i) {
        current_q = *q_iter != nullptr ? *q_iter++ : current_q + 1;
        current_p = *p_iter != nullptr ? *p_iter++ : current_p + 1;
        q_values[i] = boundary_to_public_g2_affine(*current_q);
        p_values[i] = boundary_to_public_g1_affine(*current_p);
    }

    PublicGt out = {};
    PublicGt::miller_loop_n(out, q_values.data(), p_values.data(), n);
    *ret = public_to_boundary_gt(out);
}

void blst_precompute_lines(BoundaryFp6 Qlines[68], const BoundaryG2Affine* Q)
{
    const PublicG2Affine q = Q == nullptr ? PublicG2Affine::generator() : boundary_to_public_g2_affine(*Q);
    PublicGt::precompute_lines(Qlines, q);
}

void blst_miller_loop_lines(BoundaryFp12* ret, const BoundaryFp6 Qlines[68], const BoundaryG1Affine* P)
{
    PublicGt out = {};
    const PublicG1Affine p = P == nullptr ? PublicG1Affine::generator() : boundary_to_public_g1_affine(*P);
    PublicGt::miller_loop_lines(out, Qlines, p);
    *ret = public_to_boundary_gt(out);
}

size_t blst_pairing_sizeof(void)
{
    return blst::blst_pairing_sizeof();
}

size_t blst_pairing_alignof(void)
{
    return alignof(BoundaryPairing);
}

void blst_pairing_init(void* new_ctx, bool hash_or_encode, const unsigned char* DST, size_t DST_len)
{
    blst::blst_pairing_init(as_pairing(new_ctx), hash_or_encode ? 1 : 0, DST, DST_len);
}

const unsigned char* blst_pairing_get_dst(const void* ctx)
{
    return reinterpret_cast<const unsigned char*>(blst::blst_pairing_get_dst(as_pairing(ctx)));
}

void blst_pairing_commit(void* ctx)
{
    blst::blst_pairing_commit(as_pairing(ctx));
}

blst::BLST_ERROR blst_pairing_aggregate_pk_in_g2(
    void* ctx,
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
    void* ctx,
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
    void* ctx,
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
    void* ctx,
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
    void* ctx,
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
    void* ctx,
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
    void* ctx,
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
    void* ctx,
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

blst::BLST_ERROR blst_pairing_merge(void* ctx, const void* ctx1)
{
    return blst::blst_pairing_merge(as_pairing(ctx), as_pairing(ctx1));
}

bool blst_pairing_finalverify(const void* ctx, const BoundaryFp12* gtsig)
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
    PublicGt gt = {};
    const PublicG1Affine sig = boundary_to_public_g1_affine(*signature);
    blst::blst_aggregated_in_g1(gt, sig);
    *out = public_to_boundary_gt(gt);
}

void blst_aggregated_in_g2(BoundaryFp12* out, const BoundaryG2Affine* signature)
{
    PublicGt gt = {};
    const PublicG2Affine sig = boundary_to_public_g2_affine(*signature);
    blst::blst_aggregated_in_g2(gt, sig);
    *out = public_to_boundary_gt(gt);
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

void blst_pairing_raw_aggregate(void* ctx, const BoundaryG2Affine* q, const BoundaryG1Affine* p)
{
    const PublicG2Affine public_q = boundary_to_public_g2_affine(*q);
    const PublicG1Affine public_p = boundary_to_public_g1_affine(*p);
    blst::blst_pairing_raw_aggregate(as_pairing(ctx), &public_q, &public_p);
}

BoundaryFp12 blst_pairing_as_fp12(void* ctx)
{
    return blst::blst_pairing_as_fp12(as_pairing(ctx));
}

void blst_pairing_as_fp12_to(BoundaryFp12* out, void* ctx)
{
    blst::blst_pairing_as_fp12_to(out, as_pairing(ctx));
}

void blst_keygen(blst::byte* out_SK, const blst::byte* IKM, size_t IKM_len, const blst::byte* info, size_t info_len)
{
    blst::blst_keygen(out_SK, IKM, IKM_len, info, info_len);
}

void blst_keygen_v3(blst::byte* out_SK, const blst::byte* IKM, size_t IKM_len, const blst::byte* info, size_t info_len)
{
    blst::blst_keygen_v3(out_SK, IKM, IKM_len, info, info_len);
}

void blst_keygen_v4_5(
    blst::byte* out_SK,
    const blst::byte* IKM,
    size_t IKM_len,
    const blst::byte* salt,
    size_t salt_len,
    const blst::byte* info,
    size_t info_len)
{
    blst::blst_keygen_v4_5(out_SK, IKM, IKM_len, salt, salt_len, info, info_len);
}

void blst_keygen_v5(
    blst::byte* out_SK,
    const blst::byte* IKM,
    size_t IKM_len,
    const blst::byte* salt,
    size_t salt_len,
    const blst::byte* info,
    size_t info_len)
{
    blst::blst_keygen_v5(out_SK, IKM, IKM_len, salt, salt_len, info, info_len);
}

void blst_derive_master_eip2333(blst::byte* out_SK, const blst::byte* IKM, size_t IKM_len)
{
    blst::blst_derive_master_eip2333(out_SK, IKM, IKM_len);
}

void blst_derive_child_eip2333(blst::byte* out_SK, const blst::byte* SK, unsigned int child_index)
{
    blst::blst_derive_child_eip2333(out_SK, SK, child_index);
}

size_t blst_uniq_sizeof(size_t n_nodes)
{
    return blst::blst_uniq_sizeof(n_nodes);
}

void blst_uniq_init(void* tree)
{
    blst::blst_uniq_init(as_uniq(tree));
}

bool blst_uniq_test(void* tree, const blst::byte* msg, size_t len)
{
    return blst::blst_uniq_test(as_uniq(tree), msg, len) != 0;
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
    PublicG1 public_out = {};
    blst::Encode_to_curve(public_out, msg, msg_len, DST, DST_len, aug, aug_len);
    *out = public_to_boundary_g1(public_out);
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
    PublicG1 public_out = {};
    blst::Hash_to_curve(public_out, msg, msg_len, DST, DST_len, aug, aug_len);
    *out = public_to_boundary_g1(public_out);
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
    PublicG2 public_out = {};
    blst::Encode_to_curve(public_out, msg, msg_len, DST, DST_len, aug, aug_len);
    *out = public_to_boundary_g2(public_out);
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
    PublicG2 public_out = {};
    blst::Hash_to_curve(public_out, msg, msg_len, DST, DST_len, aug, aug_len);
    *out = public_to_boundary_g2(public_out);
}

void blst_map_to_g1(BoundaryG1* out, const blst::Fp* u, const blst::Fp* v)
{
    blst::scalar::map_to_curve(*out, *u, v);
}

void blst_map_to_g2(BoundaryG2* out, const blst::Fp2* u, const blst::Fp2* v)
{
    blst::scalar::map_to_curve(*out, *u, v);
}

void blst_sk_to_pk_in_g1(BoundaryG1* out_pk, const blst::byte* SK)
{
    PublicG1 public_out = {};
    PublicG1::to_pk(public_out, SK);
    *out_pk = public_to_boundary_g1(public_out);
}

void blst_sign_pk_in_g1(BoundaryG2* out_sig, const BoundaryG2* hash, const blst::byte* SK)
{
    PublicG2 public_out = {};
    PublicG2::sign(public_out, boundary_to_public_g2(*hash), SK);
    *out_sig = public_to_boundary_g2(public_out);
}

void blst_sk_to_pk_in_g2(BoundaryG2* out_pk, const blst::byte* SK)
{
    PublicG2 public_out = {};
    PublicG2::to_pk(public_out, SK);
    *out_pk = public_to_boundary_g2(public_out);
}

void blst_sign_pk_in_g2(BoundaryG1* out_sig, const BoundaryG1* hash, const blst::byte* SK)
{
    PublicG1 public_out = {};
    PublicG1::sign(public_out, boundary_to_public_g1(*hash), SK);
    *out_sig = public_to_boundary_g1(public_out);
}

void blst_sk_to_pk2_in_g1(unsigned char out[96], BoundaryG1Affine* out_pk, const blst::byte* SK)
{
    PublicG1Affine public_pk = {};
    PublicG1Affine* public_pk_ptr = out_pk == nullptr ? nullptr : &public_pk;
    PublicG1::to_pk2(out, public_pk_ptr, SK);
    if (out_pk != nullptr) {
        *out_pk = public_to_boundary_g1_affine(public_pk);
    }
}

void blst_sign_pk2_in_g1(unsigned char out[192], BoundaryG2Affine* out_sig, const BoundaryG2* hash, const blst::byte* SK)
{
    PublicG2 public_point = {};
    PublicG2::sign(public_point, boundary_to_public_g2(*hash), SK);
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

void blst_sk_to_pk2_in_g2(unsigned char out[192], BoundaryG2Affine* out_pk, const blst::byte* SK)
{
    PublicG2Affine public_pk = {};
    PublicG2Affine* public_pk_ptr = out_pk == nullptr ? nullptr : &public_pk;
    PublicG2::to_pk2(out, public_pk_ptr, SK);
    if (out_pk != nullptr) {
        *out_pk = public_to_boundary_g2_affine(public_pk);
    }
}

void blst_sign_pk2_in_g2(unsigned char out[96], BoundaryG1Affine* out_sig, const BoundaryG1* hash, const blst::byte* SK)
{
    PublicG1 public_point = {};
    PublicG1::sign(public_point, boundary_to_public_g1(*hash), SK);
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
