#include "blst/curve/pairing_avx.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include <array>

#include "blst/backend/avx512/avx_layout.h"
#include "blst/curve/G2_avx.h"
#include "blst/fields/Fp12_avx.h"

namespace blst::avx
{
namespace
{
constexpr std::size_t kMillerLoopNMax = 16;

inline void store_u64_lanes(std::array<u64, 8>& lanes, const cc512i v)
{
    _cc512_storeu_si512(lanes.data(), v);
}

#if defined(NDEBUG)
#define BLST_ASSERT_NO_OVERLAP_PTR(a, an, b, bn) ((void)0)
#else
inline bool avx_ptr_ranges_overlap(const void* a, std::size_t a_bytes, const void* b, std::size_t b_bytes)
{
    const auto* ap = static_cast<const u8*>(a);
    const auto* bp = static_cast<const u8*>(b);
    return (ap < (bp + b_bytes)) && (bp < (ap + a_bytes));
}
#define BLST_ASSERT_NO_OVERLAP_PTR(a, an, b, bn) \
    BLST_ASSERT(!avx_ptr_ranges_overlap((a), (an), (b), (bn)))
#endif

void line_by_px2(Fp6& line, const G1AffinePoint& px2)
{
    Fp::mul(line[1][0], line[1][0], px2.X);
    Fp::mul(line[1][1], line[1][1], px2.X);

    Fp::mul(line[2][0], line[2][0], px2.Y);
    Fp::mul(line[2][1], line[2][1], px2.Y);
}

void start_dbl_n(
    GT& ret,
    G2Point* BLST_restrict t,
    const G1AffinePoint* BLST_restrict px2,
    std::size_t n)
{
    Fp6 line{};

    line_dbl_vector(line, t[0], t[0]);
    line_by_px2(line, px2[0]);

    vec_zero(&ret, sizeof(GT));
    vec_copy(&ret.data[0][0], &line[0], 2 * sizeof(Fp2));
    vec_copy(&ret.data[1][1], &line[2], sizeof(Fp2));

    for (std::size_t i = 1; i < n; ++i) {
        line_dbl_vector(line, t[i], t[i]);
        line_by_px2(line, px2[i]);
        mul_by_xy00z0_fp12_vector(ret.data, ret.data, line);
    }
}

void add_n_dbl_n_generic(
    GT& ret,
    G2Point* BLST_restrict t,
    const G2AffinePoint* BLST_restrict q,
    const G1AffinePoint* BLST_restrict px2,
    std::size_t n,
    std::size_t k)
{
    Fp6 line{};

    for (std::size_t i = 0; i < n; ++i) {
        line_add_vector(line, t[i], t[i], q[i]);
        line_by_px2(line, px2[i]);
        mul_by_xy00z0_fp12_vector(ret.data, ret.data, line);
    }

    while (k--) {
        sqr_fp12_vector(ret.data, ret.data);
        for (std::size_t i = 0; i < n; ++i) {
            line_dbl_vector(line, t[i], t[i]);
            line_by_px2(line, px2[i]);
            mul_by_xy00z0_fp12_vector(ret.data, ret.data, line);
        }
    }
}

void add_n_dbl_n_vector(
    Fp12& ret,
    G2Point* BLST_restrict T,
    const G2AffinePoint* BLST_restrict Q,
    const G1AffinePoint* BLST_restrict Px2,
    std::size_t k)
{
    BLST_ASSERT(T != nullptr);
    BLST_ASSERT(Q != nullptr);
    BLST_ASSERT(Px2 != nullptr);
    BLST_ASSERT_NO_OVERLAP_PTR(T, sizeof(*T), Q, sizeof(*Q));
    BLST_ASSERT_NO_OVERLAP_PTR(T, sizeof(*T), Px2, sizeof(*Px2));
    BLST_ASSERT_NO_OVERLAP_PTR(Q, sizeof(*Q), Px2, sizeof(*Px2));

    std::array<std::array<cc512i, SWORDS>, 3> t;
    std::array<std::array<cc512i, SWORDS / 2>, 4> s;
    int i;

    fp2_2x2x2w X1Y1, Z1Y2, X2, l0Y3, l1, X3, Z3;
    const cc512i hh = VSET(7, 6, 5, 4, 7, 6, 5, 4);

    for (i = 0; i < SWORDS / 2; ++i) {
        s[0][i] = VSET(
            static_cast<i64>(T->X[1][i + (SWORDS / 2)]), static_cast<i64>(T->X[1][i]),
            static_cast<i64>(T->X[0][i + (SWORDS / 2)]), static_cast<i64>(T->X[0][i]),
            static_cast<i64>(T->Y[1][i + (SWORDS / 2)]), static_cast<i64>(T->Y[1][i]),
            static_cast<i64>(T->Y[0][i + (SWORDS / 2)]), static_cast<i64>(T->Y[0][i]));
        s[1][i] = VSET(
            static_cast<i64>(T->Z[1][i + (SWORDS / 2)]), static_cast<i64>(T->Z[1][i]),
            static_cast<i64>(T->Z[0][i + (SWORDS / 2)]), static_cast<i64>(T->Z[0][i]),
            static_cast<i64>(Q->Y[1][i + (SWORDS / 2)]), static_cast<i64>(Q->Y[1][i]),
            static_cast<i64>(Q->Y[0][i + (SWORDS / 2)]), static_cast<i64>(Q->Y[0][i]));
        s[2][i] = VSET(
            static_cast<i64>(Q->X[1][i + (SWORDS / 2)]), static_cast<i64>(Q->X[1][i]),
            static_cast<i64>(Q->X[0][i + (SWORDS / 2)]), static_cast<i64>(Q->X[0][i]),
            static_cast<i64>(Q->X[1][i + (SWORDS / 2)]), static_cast<i64>(Q->X[1][i]),
            static_cast<i64>(Q->X[0][i + (SWORDS / 2)]), static_cast<i64>(Q->X[0][i]));
    }

    conv_64to48_fp_4x2w(X1Y1, s[0].data());
    conv_64to48_fp_4x2w(Z1Y2, s[1].data());
    conv_64to48_fp_4x2w(X2, s[2].data());

    line_add_vec_v1(l0Y3, l1, X3, Z3, X1Y1, Z1Y2, X2);

    perm_var_hl(X3.data(), X3.data(), hh);
    perm_var_hl(Z3.data(), Z3.data(), hh);

    fp2_2x2x2w l12, _Px2;
    blend_0x0F_hl(l12.data(), Z3.data(), l1.data());
    for (i = 0; i < SWORDS / 2; ++i) {
        s[0][i] = VSET(
            static_cast<i64>(Px2->Y[i + (SWORDS / 2)]), static_cast<i64>(Px2->Y[i]),
            static_cast<i64>(Px2->Y[i + (SWORDS / 2)]), static_cast<i64>(Px2->Y[i]),
            static_cast<i64>(Px2->X[i + (SWORDS / 2)]), static_cast<i64>(Px2->X[i]),
            static_cast<i64>(Px2->X[i + (SWORDS / 2)]), static_cast<i64>(Px2->X[i]));
    }
    conv_64to48_fp_4x2w(_Px2, s[0].data());

    line_by_Px2_2x2x2w(l12, l12, _Px2);

    fp2_4x2x1w a01, a2, b01, b4, r0;
    fp2_2x2x2w r1;
    const cc512i m0 = VSET(6, 4, 6, 4, 6, 4, 6, 4);
    const cc512i m1 = VSET(7, 5, 7, 5, 7, 5, 7, 5);
    const cc512i m2 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m3 = VSET(6, 4, 2, 0, 6, 4, 2, 0);
    const cc512i m4 = VSET(7, 5, 3, 1, 7, 5, 3, 1);

    for (i = 0; i < SWORDS; ++i) {
        t[0][i] = VSET(
            static_cast<i64>(ret[1][1][1][i]), static_cast<i64>(ret[1][1][0][i]),
            static_cast<i64>(ret[1][0][1][i]), static_cast<i64>(ret[1][0][0][i]),
            static_cast<i64>(ret[0][1][1][i]), static_cast<i64>(ret[0][1][0][i]),
            static_cast<i64>(ret[0][0][1][i]), static_cast<i64>(ret[0][0][0][i]));
        t[1][i] = VSET(
            static_cast<i64>(ret[1][2][1][i]), static_cast<i64>(ret[1][2][0][i]),
            static_cast<i64>(ret[1][2][1][i]), static_cast<i64>(ret[1][2][0][i]),
            static_cast<i64>(ret[0][2][1][i]), static_cast<i64>(ret[0][2][0][i]),
            static_cast<i64>(ret[0][2][1][i]), static_cast<i64>(ret[0][2][0][i]));
    }
    conv_64to48_fp_8x1w(a01, t[0].data());
    conv_64to48_fp_8x1w(a2, t[1].data());

    carryp_fp_4x2w(l12);
    carryp_fp_4x2w(l0Y3);

    for (i = 0; i < VWORDS; ++i) {
        b4[i] = VPERMV(m0, l12[i]);
        b4[i + VWORDS] = VPERMV(m1, l12[i]);
    }
    perm_var_hl(l12.data(), l12.data(), m2);
    blend_0x0F_hl(l12.data(), l12.data(), l0Y3.data());
    for (i = 0; i < VWORDS; ++i) {
        b01[i] = VPERMV(m3, l12[i]);
        b01[i + VWORDS] = VPERMV(m4, l12[i]);
    }

    mul_by_xy00z0_fp12_vec_v1(r0, r1, a01, a2, b01, b4);

    fp2_4x2x1w a0, a1, _r0, _r1;
    fp2_4x2x1w X1Y1Z1, _l0, _l12, _X3, _Y3, _Z3, __Px2;
    fp2_2x2x2w t0;
    const cc512i m5 = VSET(0, 0, 1, 0, 3, 2, 0, 0);
    const cc512i m6 = VSET(4, 4, 0, 0, 0, 0, 2, 0);
    const cc512i m7 = VSET(5, 5, 0, 0, 0, 0, 3, 1);
    const cc512i m8 = VSET(6, 6, 6, 4, 0, 0, 0, 0);
    const cc512i m9 = VSET(7, 7, 7, 5, 0, 0, 0, 0);
    const cc512i m10 = VSET(5, 4, 5, 4, 5, 4, 5, 4);
    const cc512i m11 = VSET(6, 4, 6, 4, 2, 0, 2, 0);
    const cc512i m12 = VSET(7, 5, 7, 5, 3, 1, 3, 1);
    const cc512i m13 = VSET(3, 2, 3, 2, 3, 2, 3, 2);
    const cc512i m14 = VSET(5, 4, 5, 4, 5, 4, 5, 4);
    const cc512i m15 = VSET(0, 0, 7, 6, 1, 0, 0, 0);
    const cc512i m16 = VSET(7, 6, 7, 6, 7, 6, 7, 6);
    int first_iteration = 1;

    for (i = 0; i < SWORDS; ++i) {
        t[0][i] = VSET(
            static_cast<i64>(Px2->X[i]), static_cast<i64>(Px2->X[i]),
            static_cast<i64>(Px2->X[i]), static_cast<i64>(Px2->X[i]),
            static_cast<i64>(Px2->Y[i]), static_cast<i64>(Px2->Y[i]),
            static_cast<i64>(Px2->Y[i]), static_cast<i64>(Px2->Y[i]));
    }
    conv_64to48_fp_8x1w(__Px2, t[0].data());

    while (k--) {
        perm_var(a0.data(), r0.data(), m5);
        perm_var(a1.data(), r0.data(), m2);
        for (i = 0; i < VWORDS; ++i) {
            _r0[i] = VPERMV(m6, r1[i]);
            _r0[i + VWORDS] = VPERMV(m7, r1[i]);
            _r1[i] = VPERMV(m8, r1[i]);
            _r1[i + VWORDS] = VPERMV(m9, r1[i]);
        }
        blend_0xC3(a0.data(), a0.data(), _r0.data());
        blend_0x0F(a1.data(), _r1.data(), a1.data());

        sqr_fp12_vec_v1(_r0, _r1, a0, a1);

        if (first_iteration) {
            blend_0x0F_hl(t0.data(), X3.data(), Z3.data());
            for (i = 0; i < VWORDS; ++i) {
                X1Y1Z1[i] = VPERMV(m11, t0[i]);
                X1Y1Z1[i + VWORDS] = VPERMV(m12, t0[i]);
            }
            for (i = 0; i < VWORDS; ++i) {
                _Y3[i] = VPERMV(m0, l0Y3[i]);
                _Y3[i + VWORDS] = VPERMV(m1, l0Y3[i]);
            }
            blend_0xC3(X1Y1Z1.data(), X1Y1Z1.data(), _Y3.data());
            first_iteration = 0;
        } else {
            blend_0x0F(X1Y1Z1.data(), _X3.data(), _Z3.data());
            perm_var(X1Y1Z1.data(), X1Y1Z1.data(), m15);
            perm_var(_Y3.data(), _Y3.data(), m16);
            blend_0xC3(X1Y1Z1.data(), X1Y1Z1.data(), _Y3.data());
        }

        line_dbl_vec_v2(_l0, _l12, _X3, _Y3, _Z3, X1Y1Z1);

        line_by_Px2_4x2x1w(_l12, _l12, __Px2);

        perm_var(b01.data(), _l12.data(), m14);
        perm_var(b4.data(), _l0.data(), m13);
        blend_0x33(b01.data(), b01.data(), b4.data());
        perm_var(b4.data(), _l12.data(), m13);

        perm_var(a01.data(), _r1.data(), m2);
        blend_0x0F(a01.data(), a01.data(), _r0.data());
        perm_var(_r0.data(), _r0.data(), m10);
        perm_var(_r1.data(), _r1.data(), m10);
        blend_0x0F(a2.data(), _r1.data(), _r0.data());

        mul_by_xy00z0_fp12_vec_v1(r0, r1, a01, a2, b01, b4);
    }

    carryp_fp_4x2w(r1);
    conv_48to64_fp_8x1w(t[0].data(), r0);
    conv_48to64_fp_4x2w(s[0].data(), r1);

    for (i = 0; i < SWORDS; ++i) {
        std::array<u64, 8> lanes;
        store_u64_lanes(lanes, t[0][i]);
        ret[0][2][0][i] = lanes[0];
        ret[0][2][1][i] = lanes[1];
        ret[0][1][0][i] = lanes[2];
        ret[0][1][1][i] = lanes[3];
        ret[1][0][0][i] = lanes[4];
        ret[1][0][1][i] = lanes[5];
        ret[1][1][0][i] = lanes[6];
        ret[1][1][1][i] = lanes[7];
    }

    for (i = 0; i < SWORDS / 2; ++i) {
        std::array<u64, 8> lanes;
        store_u64_lanes(lanes, s[0][i]);
        ret[0][0][0][i] = lanes[0];
        ret[0][0][0][i + (SWORDS / 2)] = lanes[1];
        ret[0][0][1][i] = lanes[2];
        ret[0][0][1][i + (SWORDS / 2)] = lanes[3];
        ret[1][2][0][i] = lanes[4];
        ret[1][2][0][i + (SWORDS / 2)] = lanes[5];
        ret[1][2][1][i] = lanes[6];
        ret[1][2][1][i + (SWORDS / 2)] = lanes[7];
    }

    conv_48to64_fp_8x1w(t[0].data(), _X3);
    conv_48to64_fp_8x1w(t[1].data(), _Y3);
    conv_48to64_fp_8x1w(t[2].data(), _Z3);

    for (i = 0; i < SWORDS; ++i) {
        std::array<u64, 8> lanes0;
        std::array<u64, 8> lanes1;
        std::array<u64, 8> lanes2;
        store_u64_lanes(lanes0, t[0][i]);
        store_u64_lanes(lanes1, t[1][i]);
        store_u64_lanes(lanes2, t[2][i]);
        T->X[0][i] = lanes0[6];
        T->X[1][i] = lanes0[7];
        T->Y[0][i] = lanes1[6];
        T->Y[1][i] = lanes1[7];
        T->Z[0][i] = lanes2[0];
        T->Z[1][i] = lanes2[1];
    }
}

void add_n_dbl_n(
    GT& ret,
    G2Point* BLST_restrict t,
    const G2AffinePoint* BLST_restrict q,
    const G1AffinePoint* BLST_restrict px2,
    std::size_t n,
    std::size_t k)
{
    if (n == 1) {
        add_n_dbl_n_vector(ret.data, t, q, px2, k);
        return;
    }

    add_n_dbl_n_generic(ret, t, q, px2, n, k);
}

void mul_n_sqr_vector(Fp12& ret, const Fp12& a, std::size_t n)
{
    fp2_4x2x1w ab0, ab1, ab2;
    fp2_2x2x2w r001, r101, r2;
    fp4_1x2x2x2w _a;
    fp4_2x2x2x1w _bc;
    std::array<std::array<cc512i, SWORDS>, 3> t;
    std::array<cc512i, SWORDS / 2> s;
    int i;

    for (i = 0; i < SWORDS; ++i) {
        t[0][i] = VSET(
            static_cast<i64>(ret[1][0][1][i]), static_cast<i64>(ret[1][0][0][i]),
            static_cast<i64>(a[1][0][1][i]), static_cast<i64>(a[1][0][0][i]),
            static_cast<i64>(ret[0][0][1][i]), static_cast<i64>(ret[0][0][0][i]),
            static_cast<i64>(a[0][0][1][i]), static_cast<i64>(a[0][0][0][i]));
        t[1][i] = VSET(
            static_cast<i64>(ret[1][1][1][i]), static_cast<i64>(ret[1][1][0][i]),
            static_cast<i64>(a[1][1][1][i]), static_cast<i64>(a[1][1][0][i]),
            static_cast<i64>(ret[0][1][1][i]), static_cast<i64>(ret[0][1][0][i]),
            static_cast<i64>(a[0][1][1][i]), static_cast<i64>(a[0][1][0][i]));
        t[2][i] = VSET(
            static_cast<i64>(ret[1][2][1][i]), static_cast<i64>(ret[1][2][0][i]),
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(ret[0][2][1][i]), static_cast<i64>(ret[0][2][0][i]),
            static_cast<i64>(a[0][2][1][i]), static_cast<i64>(a[0][2][0][i]));
    }

    conv_64to48_fp_8x1w(ab0, t[0].data());
    conv_64to48_fp_8x1w(ab1, t[1].data());
    conv_64to48_fp_8x1w(ab2, t[2].data());

    mul_fp12_vec_v3(r001, r101, r2, ab0, ab1, ab2);

    carryp_fp_4x2w(r001);
    carryp_fp_4x2w(r101);
    carryp_fp_4x2w(r2);

    blend_0x0F_hl(_a.data(), r101.data(), r001.data());

    for (i = 0; i < VWORDS; ++i) {
        const u64* lanes_r2 = reinterpret_cast<const u64*>(&r2[i]);
        const u64* lanes_r001 = reinterpret_cast<const u64*>(&r001[i]);
        const u64* lanes_r101 = reinterpret_cast<const u64*>(&r101[i]);
        _bc[i] = VSET(
            static_cast<i64>(lanes_r2[6]),
            static_cast<i64>(lanes_r2[4]),
            static_cast<i64>(lanes_r001[6]),
            static_cast<i64>(lanes_r001[4]),
            static_cast<i64>(lanes_r2[2]),
            static_cast<i64>(lanes_r2[0]),
            static_cast<i64>(lanes_r101[2]),
            static_cast<i64>(lanes_r101[0]));
        _bc[i + VWORDS] = VSET(
            static_cast<i64>(lanes_r2[7]),
            static_cast<i64>(lanes_r2[5]),
            static_cast<i64>(lanes_r001[7]),
            static_cast<i64>(lanes_r001[5]),
            static_cast<i64>(lanes_r2[3]),
            static_cast<i64>(lanes_r2[1]),
            static_cast<i64>(lanes_r101[3]),
            static_cast<i64>(lanes_r101[1]));
    }

    while (n--) {
        cyclotomic_sqr_fp12_vec_v1(_a, _bc, _a, _bc);
    }

    carryp_fp_4x2w(_a);
    conv_48to64_fp_4x2w(s.data(), _a);
    for (i = 0; i < SWORDS / 2; ++i) {
        const u64* lanes = reinterpret_cast<const u64*>(&s[i]);
        ret[0][0][0][i] = lanes[0];
        ret[0][0][0][i + (SWORDS / 2)] = lanes[1];
        ret[0][0][1][i] = lanes[2];
        ret[0][0][1][i + (SWORDS / 2)] = lanes[3];
        ret[1][1][0][i] = lanes[4];
        ret[1][1][0][i + (SWORDS / 2)] = lanes[5];
        ret[1][1][1][i] = lanes[6];
        ret[1][1][1][i + (SWORDS / 2)] = lanes[7];
    }

    conv_48to64_fp_8x1w(t[0].data(), _bc);
    for (i = 0; i < SWORDS; ++i) {
        const u64* lanes = reinterpret_cast<const u64*>(&t[0][i]);
        ret[1][0][0][i] = lanes[0];
        ret[1][0][1][i] = lanes[1];
        ret[0][2][0][i] = lanes[2];
        ret[0][2][1][i] = lanes[3];
        ret[0][1][0][i] = lanes[4];
        ret[0][1][1][i] = lanes[5];
        ret[1][2][0][i] = lanes[6];
        ret[1][2][1][i] = lanes[7];
    }
}

void mul_n_sqr(Fp12& ret, const Fp12& a, std::size_t n)
{
    mul_n_sqr_vector(ret, a, n);
}

void raise_to_z_div_by_2(Fp12& ret, const Fp12& a)
{
    cyclotomic_sqr_fp12_vector(ret, a);
    mul_n_sqr(ret, a, 2);
    mul_n_sqr(ret, a, 3);
    mul_n_sqr(ret, a, 9);
    mul_n_sqr(ret, a, 32);
    mul_n_sqr(ret, a, 16 - 1);
    Fp12::conjugate(ret);
}

void raise_to_z(Fp12& ret, const Fp12& a)
{
    raise_to_z_div_by_2(ret, a);
    cyclotomic_sqr_fp12_vector(ret, ret);
}
} // namespace

void miller_loop_n(
    GT& ret,
    const G2AffinePoint* BLST_restrict q,
    const G1AffinePoint* BLST_restrict p,
    std::size_t n)
{
    if (n == 0) {
        vec_copy(&ret.data, get_BLST_Rx().p12, sizeof(Fp12));
        return;
    }

    BLST_ASSERT(q != nullptr);
    BLST_ASSERT(p != nullptr);
    if ((q == nullptr) || (p == nullptr)) {
        return;
    }
    BLST_ASSERT_NO_OVERLAP_PTR(&ret, sizeof(ret), q, sizeof(*q) * n);
    BLST_ASSERT_NO_OVERLAP_PTR(&ret, sizeof(ret), p, sizeof(*p) * n);
    BLST_ASSERT_NO_OVERLAP_PTR(q, sizeof(*q) * n, p, sizeof(*p) * n);

    if (n > kMillerLoopNMax) {
        return;
    }

    if ((n == 1) && (vec_is_zero(&q[0], sizeof(q[0])) |
                     vec_is_zero(&p[0], sizeof(p[0])))) {
        vec_copy(&ret.data, get_BLST_Rx().p12, sizeof(Fp12));
        return;
    }

    std::array<G2Point, kMillerLoopNMax> t;
    std::array<G1AffinePoint, kMillerLoopNMax> px2;

    for (std::size_t i = 0; i < n; ++i) {
        Fp::add(px2[i].X, p[i].X, p[i].X);
        Fp::neg(px2[i].X, px2[i].X);
        Fp::add(px2[i].Y, p[i].Y, p[i].Y);

        vec_copy(&t[i].X, &q[i].X, 2 * sizeof(t[i].X));
        vec_copy(&t[i].Z, &get_BLST_Rx().p2, sizeof(t[i].Z));
    }

    start_dbl_n(ret, t.data(), px2.data(), n);
    add_n_dbl_n(ret, t.data(), q, px2.data(), n, 2);
    add_n_dbl_n(ret, t.data(), q, px2.data(), n, 3);
    add_n_dbl_n(ret, t.data(), q, px2.data(), n, 9);
    add_n_dbl_n(ret, t.data(), q, px2.data(), n, 32);
    add_n_dbl_n(ret, t.data(), q, px2.data(), n, 16);
    Fp12::conjugate(ret.data);
}

void miller_loop_n(
    GT& ret,
    std::span<const G2AffinePoint> q,
    std::span<const G1AffinePoint> p)
{
    BLST_ASSERT(q.size() == p.size());
    if (q.size() != p.size()) {
        return;
    }

    const auto* BLST_restrict q_ptr = q.data();
    const auto* BLST_restrict p_ptr = p.data();
    miller_loop_n(ret, q_ptr, p_ptr, q.size());
}

void final_exp(GT& ret, const GT& f)
{
    Fp12 y0{};
    Fp12 y1{};
    Fp12 y2{};
    Fp12 y3{};

    vec_copy(&y1, &f.data, sizeof(y1));
    Fp12::conjugate(y1);
    Fp12::inverse(y2, f.data);
    mul_fp12_vector(ret.data, y1, y2);
    Fp12::frobenius_map(y2, ret.data, 2);
    mul_fp12_vector(ret.data, ret.data, y2);

    cyclotomic_sqr_fp12_vector(y0, ret.data);
    raise_to_z(y1, y0);
    raise_to_z_div_by_2(y2, y1);
    vec_copy(&y3, &ret.data, sizeof(y3));
    Fp12::conjugate(y3);
    mul_fp12_vector(y1, y1, y3);
    Fp12::conjugate(y1);
    mul_fp12_vector(y1, y1, y2);
    raise_to_z(y2, y1);
    raise_to_z(y3, y2);
    Fp12::conjugate(y1);
    mul_fp12_vector(y3, y3, y1);
    Fp12::conjugate(y1);
    Fp12::frobenius_map(y1, y1, 3);
    Fp12::frobenius_map(y2, y2, 2);
    mul_fp12_vector(y1, y1, y2);
    raise_to_z(y2, y3);
    mul_fp12_vector(y2, y2, y0);
    mul_fp12_vector(y2, y2, ret.data);
    mul_fp12_vector(y1, y1, y2);
    Fp12::frobenius_map(y2, y3, 1);
    mul_fp12_vector(ret.data, y1, y2);
}

void optimal_ate_pairing(
    GT& ret,
    const G2AffinePoint* BLST_restrict q,
    const G1AffinePoint* BLST_restrict p,
    std::size_t n)
{
    GT f{};
    miller_loop_n(f, q, p, n);
    final_exp(ret, f);
}

void optimal_ate_pairing(
    GT& ret,
    std::span<const G2AffinePoint> q,
    std::span<const G1AffinePoint> p)
{
    BLST_ASSERT(q.size() == p.size());
    if (q.size() != p.size()) {
        return;
    }

    const auto* BLST_restrict q_ptr = q.data();
    const auto* BLST_restrict p_ptr = p.data();
    optimal_ate_pairing(ret, q_ptr, p_ptr, q.size());
}

} // namespace blst::avx

#endif
