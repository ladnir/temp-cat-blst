#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fields/Fp12.h"
#include "blst/scalar/curve/GT.h"
#include "blst/scalar/curve/point.h"
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/curve/map_to_G.h"
#endif

#include <cstddef>
#include <type_traits>
#include <vector>

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

} // namespace

BLST_EXTERN_C_BEGIN_

void blst_fp12_mul(BoundaryFp12* ret, const BoundaryFp12* a, const BoundaryFp12* b)
{
    BoundaryFp12::mul(*ret, *a, *b);
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

void blst_p1_mult(BoundaryG1* out, const BoundaryG1* p, const byte* scalar, size_t nbits)
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

void blst_p2_mult(BoundaryG2* out, const BoundaryG2* p, const byte* scalar, size_t nbits)
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
    std::vector<PublicG2Affine> q_values(n);
    std::vector<PublicG1Affine> p_values(n);
    const BoundaryG2Affine* current_q = nullptr;
    const BoundaryG1Affine* current_p = nullptr;

    for (size_t i = 0; i < n; ++i) {
        current_q = Qs[i] != nullptr ? Qs[i] : current_q + 1;
        current_p = Ps[i] != nullptr ? Ps[i] : current_p + 1;
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

BLST_EXTERN_C_END_
