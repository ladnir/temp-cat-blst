#include <gtest/gtest.h>

#include <array>

#include "InsecureRng.h"
#include "blst/backend/avx512/avx_layout.h"
#include "blst/curve/point.h"
#include "blst/curve/G2_avx.h"
#include "blst/fields/Fp2.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

using namespace blst;

namespace
{
struct LineEval
{
    Fp2 c0;
    Fp2 c1;
    Fp2 c2;
};

struct JacobianPoint
{
    Fp2 X;
    Fp2 Y;
    Fp2 Z;
};

struct AffineInputPoint
{
    Fp2 X;
    Fp2 Y;
};

void rand_redc_fp(Fp& v, InsecureRng& rng)
{
    for (u64 i = 0; i < v.num_limb(); ++i) {
        v[i] = rng.next_u64();
    }
    while (Fp::is_redc(v) == 0) {
        Fp::rshift(v, v, 1);
    }
}

void rand_redc_fp2(Fp2& v, InsecureRng& rng)
{
    rand_redc_fp(v[0], rng);
    rand_redc_fp(v[1], rng);
}

void fill_edge_fp2(Fp2& v, int idx)
{
    Fp two;
    Fp p_minus_one;
    Fp::add(two, Fp::one(), Fp::one());
    Fp::sub(p_minus_one, Fp::zero(), Fp::one());
    const Fp* cycle[4] = {&Fp::zero(), &Fp::one(), &two, &p_minus_one};
    v[0] = *cycle[(idx + 0) & 3];
    v[1] = *cycle[(idx + 1) & 3];
}

void make_point(JacobianPoint& p, InsecureRng& rng, int edge_idx)
{
    if (edge_idx >= 0) {
        fill_edge_fp2(p.X, edge_idx + 0);
        fill_edge_fp2(p.Y, edge_idx + 1);
        fill_edge_fp2(p.Z, edge_idx + 2);
        return;
    }
    rand_redc_fp2(p.X, rng);
    rand_redc_fp2(p.Y, rng);
    rand_redc_fp2(p.Z, rng);
}

void make_affine(AffineInputPoint& p, InsecureRng& rng, int edge_idx)
{
    if (edge_idx >= 0) {
        fill_edge_fp2(p.X, edge_idx + 0);
        fill_edge_fp2(p.Y, edge_idx + 1);
        return;
    }
    rand_redc_fp2(p.X, rng);
    rand_redc_fp2(p.Y, rng);
}

void scalar_u64_le(byte* out, u64 v)
{
    for (int i = 0; i < 8; ++i) {
        out[i] = static_cast<byte>(v >> (8 * i));
    }
}

void make_valid_g2_jacobian_from_u64(JacobianPoint& out, u64 scalar_u64)
{
    std::array<byte, 8> s = {};
    scalar_u64_le(s.data(), scalar_u64 == 0 ? 1 : scalar_u64);

    G2Point p{};
    G2Point::mult(p, G2Point::generator(), s.data(), 64);
    if (G2Point::is_inf(p)) {
        p = G2Point::generator();
    }

    out.X = p.X;
    out.Y = p.Y;
    out.Z = p.Z;
}

void make_valid_g2_affine_from_u64(AffineInputPoint& out, u64 scalar_u64)
{
    JacobianPoint j{};
    make_valid_g2_jacobian_from_u64(j, scalar_u64);

    G2Point p{};
    p.X = j.X;
    p.Y = j.Y;
    p.Z = j.Z;

    G2AffinePoint a{};
    G2Point::to_affine(a, p);
    out.X = a.X;
    out.Y = a.Y;
}

bool is_line_add_nonexceptional(const JacobianPoint& r, const AffineInputPoint& q)
{
    if (vec_is_zero(&r.Z, sizeof(r.Z))) {
        return false;
    }
    if (vec_is_zero(&q.X, sizeof(q.X)) && vec_is_zero(&q.Y, sizeof(q.Y))) {
        return false;
    }

    Fp2 z1z1{}, u2{}, h{};
    Fp2::sqr(z1z1, r.Z);
    Fp2::mul(u2, q.X, z1z1);
    Fp2::sub(h, u2, r.X);
    return vec_is_zero(&h, sizeof(h)) == 0;
}

void make_valid_line_add_inputs(JacobianPoint& r, AffineInputPoint& q, InsecureRng& rng, int iter)
{
    for (int tries = 0; tries < 256; ++tries) {
        const u64 sr = (iter < 16)
            ? static_cast<u64>(iter + 1 + tries)
            : (rng.next_u64() | 1ULL);
        const u64 sq = (iter < 16)
            ? static_cast<u64>(iter + 33 + tries)
            : (rng.next_u64() | 1ULL);

        make_valid_g2_jacobian_from_u64(r, sr);
        make_valid_g2_affine_from_u64(q, sq);

        if (is_line_add_nonexceptional(r, q)) {
            return;
        }
    }

    make_valid_g2_jacobian_from_u64(r, 1);
    make_valid_g2_affine_from_u64(q, 2);
}

bool fp2_equal(const Fp2& a, const Fp2& b)
{
    return Fp::equal(a[0], b[0]) && Fp::equal(a[1], b[1]);
}

void expect_fp2_eq(const Fp2& got, const Fp2& expected)
{
    EXPECT_TRUE(fp2_equal(got, expected));
}

Fp fp_from_48_channel_8x1w(const fp_8x1w a, int ch, bool carry)
{
    u64 r48[NWORDS] = {};
    u64 r64[SWORDS] = {};
    Fp out{};
    get_channel_8x1w(r48, a, ch);
    if (carry) {
        carryp_mpi48(r48);
    }
    conv_48to64_mpi(r64, r48, static_cast<int>(SWORDS), static_cast<int>(NWORDS));
    for (u64 i = 0; i < SWORDS; ++i) {
        out[i] = r64[i];
    }
    return out;
}

Fp fp_from_48_channel_4x2w(const fp_4x2w a, int ch, bool carry)
{
    u64 r48[NWORDS] = {};
    u64 r64[SWORDS] = {};
    Fp out{};
    get_channel_4x2w(r48, a, ch);
    if (carry) {
        carryp_mpi48(r48);
    }
    conv_48to64_mpi(r64, r48, static_cast<int>(SWORDS), static_cast<int>(NWORDS));
    for (u64 i = 0; i < SWORDS; ++i) {
        out[i] = r64[i];
    }
    return out;
}

Fp2 fp2_from_8x1w_channels(const fp_8x1w a, int ch0, int ch1, bool carry)
{
    Fp2 out{};
    out[0] = fp_from_48_channel_8x1w(a, ch0, carry);
    out[1] = fp_from_48_channel_8x1w(a, ch1, carry);
    return out;
}

Fp2 fp2_from_4x2w_channels(const fp_4x2w a, int ch0, int ch1, bool carry)
{
    Fp2 out{};
    out[0] = fp_from_48_channel_4x2w(a, ch0, carry);
    out[1] = fp_from_48_channel_4x2w(a, ch1, carry);
    return out;
}

void line_dbl_scalar(LineEval& line, JacobianPoint& T, const JacobianPoint& Q)
{
    Fp2 ZZ{}, A{}, B{}, C{}, D{}, E{}, F{};

    Fp2::sqr(A, Q.X);
    Fp2::sqr(B, Q.Y);
    Fp2::sqr(ZZ, Q.Z);
    Fp2::sqr(C, B);

    Fp2::add(D, Q.X, B);
    Fp2::sqr(D, D);
    Fp2::sub(D, D, A);
    Fp2::sub(D, D, C);
    Fp2::add(D, D, D);

    Fp2::mul_by_3(E, A);
    Fp2::sqr(F, E);

    Fp2::add(line.c0, E, Q.X);

    Fp2::sub(T.X, F, D);
    Fp2::sub(T.X, T.X, D);

    Fp2::add(T.Z, Q.Y, Q.Z);
    Fp2::sqr(T.Z, T.Z);
    Fp2::sub(T.Z, T.Z, B);
    Fp2::sub(T.Z, T.Z, ZZ);

    Fp2::mul_by_8(C, C);
    Fp2::sub(T.Y, D, T.X);
    Fp2::mul(T.Y, T.Y, E);
    Fp2::sub(T.Y, T.Y, C);

    Fp2::sqr(line.c0, line.c0);
    Fp2::sub(line.c0, line.c0, A);
    Fp2::sub(line.c0, line.c0, F);
    Fp2::lshift(B, B, 2);
    Fp2::sub(line.c0, line.c0, B);

    Fp2::mul(line.c1, E, ZZ);
    Fp2::mul(line.c2, T.Z, ZZ);
}

void line_add_scalar(LineEval& line, JacobianPoint& T, const JacobianPoint& R, const AffineInputPoint& Q)
{
    Fp2 Z1Z1{}, U2{}, S2{}, H{}, HH{}, I{}, J{}, V{};
    Fp2 r{};

    Fp2::sqr(Z1Z1, R.Z);
    Fp2::mul(U2, Q.X, Z1Z1);

    Fp2::mul(S2, Q.Y, R.Z);
    Fp2::mul(S2, S2, Z1Z1);

    Fp2::sub(H, U2, R.X);

    Fp2::sqr(HH, H);
    Fp2::add(I, HH, HH);
    Fp2::add(I, I, I);

    Fp2::mul(J, H, I);

    Fp2::sub(r, S2, R.Y);
    Fp2::add(r, r, r);

    Fp2::mul(V, R.X, I);

    Fp2::sqr(T.X, r);
    Fp2::sub(T.X, T.X, J);
    Fp2::sub(T.X, T.X, V);
    Fp2::sub(T.X, T.X, V);

    Fp2::mul(J, J, R.Y);
    Fp2::sub(T.Y, V, T.X);
    Fp2::mul(T.Y, T.Y, r);
    Fp2::sub(T.Y, T.Y, J);
    Fp2::sub(T.Y, T.Y, J);

    Fp2::add(T.Z, R.Z, H);
    Fp2::sqr(T.Z, T.Z);
    Fp2::sub(T.Z, T.Z, Z1Z1);
    Fp2::sub(T.Z, T.Z, HH);

    Fp2::mul(I, r, Q.X);
    Fp2::mul(J, Q.Y, T.Z);
    Fp2::sub(I, I, J);
    Fp2::add(line.c0, I, I);
    line.c1 = r;
    line.c2 = T.Z;
}

void scale_line_by_px2_scalar(LineEval& line, const Fp& px, const Fp& py)
{
    Fp::mul(line.c1[0], line.c1[0], px);
    Fp::mul(line.c1[1], line.c1[1], px);
    Fp::mul(line.c2[0], line.c2[0], py);
    Fp::mul(line.c2[1], line.c2[1], py);
}

void pack_line_dbl_v2_input(fp2_4x2x1w& X1Y1Z1, const JacobianPoint& Q)
{
    cc512i t[SWORDS] = {};
    for (u64 i = 0; i < SWORDS; ++i) {
        t[i] = VSET(
            static_cast<i64>(Q.Y[1][i]), static_cast<i64>(Q.Y[0][i]),
            static_cast<i64>(Q.X[1][i]), static_cast<i64>(Q.X[0][i]),
            static_cast<i64>(Q.Z[1][i]), static_cast<i64>(Q.Z[0][i]),
            static_cast<i64>(Q.Y[1][i]), static_cast<i64>(Q.Y[0][i]));
    }
    conv_64to48_fp_8x1w(X1Y1Z1, t);
}

void pack_line_dbl_v1_inputs(fp2_2x2x2w& X1Y1, fp2_2x2x2w& Z1, const JacobianPoint& Q)
{
    cc512i s0[SWORDS / 2] = {};
    cc512i s1[SWORDS / 2] = {};
    for (u64 i = 0; i < SWORDS / 2; ++i) {
        s0[i] = VSET(
            static_cast<i64>(Q.Y[1][i + (SWORDS / 2)]), static_cast<i64>(Q.Y[1][i]),
            static_cast<i64>(Q.Y[0][i + (SWORDS / 2)]), static_cast<i64>(Q.Y[0][i]),
            static_cast<i64>(Q.X[1][i + (SWORDS / 2)]), static_cast<i64>(Q.X[1][i]),
            static_cast<i64>(Q.X[0][i + (SWORDS / 2)]), static_cast<i64>(Q.X[0][i]));
        s1[i] = VSET(
            static_cast<i64>(Q.Z[1][i + (SWORDS / 2)]), static_cast<i64>(Q.Z[1][i]),
            static_cast<i64>(Q.Z[0][i + (SWORDS / 2)]), static_cast<i64>(Q.Z[0][i]),
            static_cast<i64>(Q.Z[1][i + (SWORDS / 2)]), static_cast<i64>(Q.Z[1][i]),
            static_cast<i64>(Q.Z[0][i + (SWORDS / 2)]), static_cast<i64>(Q.Z[0][i]));
    }
    conv_64to48_fp_4x2w(X1Y1, s0);
    conv_64to48_fp_4x2w(Z1, s1);
}

void pack_line_add_v1_inputs(
    fp2_2x2x2w& X1Y1,
    fp2_2x2x2w& Z1Y2,
    fp2_2x2x2w& X2,
    const JacobianPoint& R,
    const AffineInputPoint& Q)
{
    cc512i s0[SWORDS / 2] = {};
    cc512i s1[SWORDS / 2] = {};
    cc512i s2[SWORDS / 2] = {};
    for (u64 i = 0; i < SWORDS / 2; ++i) {
        s0[i] = VSET(
            static_cast<i64>(R.X[1][i + (SWORDS / 2)]), static_cast<i64>(R.X[1][i]),
            static_cast<i64>(R.X[0][i + (SWORDS / 2)]), static_cast<i64>(R.X[0][i]),
            static_cast<i64>(R.Y[1][i + (SWORDS / 2)]), static_cast<i64>(R.Y[1][i]),
            static_cast<i64>(R.Y[0][i + (SWORDS / 2)]), static_cast<i64>(R.Y[0][i]));
        s1[i] = VSET(
            static_cast<i64>(R.Z[1][i + (SWORDS / 2)]), static_cast<i64>(R.Z[1][i]),
            static_cast<i64>(R.Z[0][i + (SWORDS / 2)]), static_cast<i64>(R.Z[0][i]),
            static_cast<i64>(Q.Y[1][i + (SWORDS / 2)]), static_cast<i64>(Q.Y[1][i]),
            static_cast<i64>(Q.Y[0][i + (SWORDS / 2)]), static_cast<i64>(Q.Y[0][i]));
        s2[i] = VSET(
            static_cast<i64>(Q.X[1][i + (SWORDS / 2)]), static_cast<i64>(Q.X[1][i]),
            static_cast<i64>(Q.X[0][i + (SWORDS / 2)]), static_cast<i64>(Q.X[0][i]),
            static_cast<i64>(Q.X[1][i + (SWORDS / 2)]), static_cast<i64>(Q.X[1][i]),
            static_cast<i64>(Q.X[0][i + (SWORDS / 2)]), static_cast<i64>(Q.X[0][i]));
    }
    conv_64to48_fp_4x2w(X1Y1, s0);
    conv_64to48_fp_4x2w(Z1Y2, s1);
    conv_64to48_fp_4x2w(X2, s2);
}

void pack_px2_4x2x1w(fp2_4x2x1w& Px2, const Fp& X, const Fp& Y)
{
    cc512i t[SWORDS] = {};
    for (u64 i = 0; i < SWORDS; ++i) {
        t[i] = VSET(
            static_cast<i64>(X[i]), static_cast<i64>(X[i]), static_cast<i64>(X[i]), static_cast<i64>(X[i]),
            static_cast<i64>(Y[i]), static_cast<i64>(Y[i]), static_cast<i64>(Y[i]), static_cast<i64>(Y[i]));
    }
    conv_64to48_fp_8x1w(Px2, t);
}

void pack_px2_2x2x2w(fp2_2x2x2w& Px2, const Fp& X, const Fp& Y)
{
    cc512i t[SWORDS / 2] = {};
    for (u64 i = 0; i < SWORDS / 2; ++i) {
        t[i] = VSET(
            static_cast<i64>(Y[i + (SWORDS / 2)]), static_cast<i64>(Y[i]),
            static_cast<i64>(Y[i + (SWORDS / 2)]), static_cast<i64>(Y[i]),
            static_cast<i64>(X[i + (SWORDS / 2)]), static_cast<i64>(X[i]),
            static_cast<i64>(X[i + (SWORDS / 2)]), static_cast<i64>(X[i]));
    }
    conv_64to48_fp_4x2w(Px2, t);
}

void decode_line_dbl_v2(
    LineEval& line,
    JacobianPoint& T,
    const fp2_4x2x1w l0,
    const fp2_4x2x1w l12,
    const fp2_4x2x1w X3,
    const fp2_4x2x1w Y3,
    const fp2_4x2x1w Z3)
{
    line.c0 = fp2_from_8x1w_channels(l0, 2, 3, false);
    line.c1 = fp2_from_8x1w_channels(l12, 4, 5, false);
    line.c2 = fp2_from_8x1w_channels(l12, 2, 3, false);
    T.X = fp2_from_8x1w_channels(X3, 6, 7, false);
    T.Y = fp2_from_8x1w_channels(Y3, 6, 7, false);
    T.Z = fp2_from_8x1w_channels(Z3, 0, 1, false);
}

void decode_line_dbl_v1(
    LineEval& line,
    JacobianPoint& T,
    const fp2_2x2x2w l0Y3,
    const fp2_2x2x2w l1,
    const fp2_2x2x2w l2,
    const fp2_2x2x2w X3,
    const fp2_2x2x2w Z3)
{
    line.c0 = fp2_from_4x2w_channels(l0Y3, 0, 2, true);
    line.c1 = fp2_from_4x2w_channels(l1, 0, 2, true);
    line.c2 = fp2_from_4x2w_channels(l2, 4, 6, true);
    T.X = fp2_from_4x2w_channels(X3, 4, 6, true);
    T.Y = fp2_from_4x2w_channels(l0Y3, 4, 6, true);
    T.Z = fp2_from_4x2w_channels(Z3, 4, 6, true);
}

void decode_line_add_v1(
    LineEval& line,
    JacobianPoint& T,
    const fp2_2x2x2w l0Y3,
    const fp2_2x2x2w l1,
    const fp2_2x2x2w X3,
    const fp2_2x2x2w Z3)
{
    line.c0 = fp2_from_4x2w_channels(l0Y3, 0, 2, true);
    line.c1 = fp2_from_4x2w_channels(l1, 0, 2, true);
    line.c2 = fp2_from_4x2w_channels(Z3, 4, 6, true);
    T.X = fp2_from_4x2w_channels(X3, 4, 6, true);
    T.Y = fp2_from_4x2w_channels(l0Y3, 4, 6, true);
    T.Z = fp2_from_4x2w_channels(Z3, 4, 6, true);
}

} // namespace

TEST(G2Avx, LineDblV2MatchesScalar)
{
    InsecureRng rng;

    for (int iter = 0; iter < 112; ++iter) {
        JacobianPoint q{};
        make_point(q, rng, (iter < 16) ? iter : -1);

        LineEval expected_line{};
        JacobianPoint expected_t{};
        line_dbl_scalar(expected_line, expected_t, q);

        fp2_4x2x1w x1y1z1{}, l0{}, l12{}, x3{}, y3{}, z3{};
        pack_line_dbl_v2_input(x1y1z1, q);
        line_dbl_vec_v2(l0, l12, x3, y3, z3, x1y1z1);

        LineEval got_line{};
        JacobianPoint got_t{};
        decode_line_dbl_v2(got_line, got_t, l0, l12, x3, y3, z3);

        expect_fp2_eq(got_line.c0, expected_line.c0);
        expect_fp2_eq(got_line.c1, expected_line.c1);
        expect_fp2_eq(got_line.c2, expected_line.c2);
        expect_fp2_eq(got_t.X, expected_t.X);
        expect_fp2_eq(got_t.Y, expected_t.Y);
        expect_fp2_eq(got_t.Z, expected_t.Z);
    }
}

TEST(G2Avx, LineDblV1MatchesScalar)
{
    InsecureRng rng;

    for (int iter = 0; iter < 112; ++iter) {
        JacobianPoint q{};
        make_point(q, rng, (iter < 16) ? iter : -1);

        LineEval expected_line{};
        JacobianPoint expected_t{};
        line_dbl_scalar(expected_line, expected_t, q);

        fp2_2x2x2w x1y1{}, z1{}, l0y3{}, l1{}, l2{}, x3{}, z3{};
        pack_line_dbl_v1_inputs(x1y1, z1, q);
        line_dbl_vec_v1(l0y3, l1, l2, x3, z3, x1y1, z1);

        LineEval got_line{};
        JacobianPoint got_t{};
        decode_line_dbl_v1(got_line, got_t, l0y3, l1, l2, x3, z3);

        expect_fp2_eq(got_line.c0, expected_line.c0);
        expect_fp2_eq(got_line.c1, expected_line.c1);
        expect_fp2_eq(got_line.c2, expected_line.c2);
        expect_fp2_eq(got_t.X, expected_t.X);
        expect_fp2_eq(got_t.Y, expected_t.Y);
        expect_fp2_eq(got_t.Z, expected_t.Z);
    }
}

TEST(G2Avx, LineDblVariantsAgree)
{
    InsecureRng rng;

    for (int iter = 0; iter < 112; ++iter) {
        JacobianPoint q{};
        make_point(q, rng, (iter < 16) ? iter : -1);

        fp2_4x2x1w x1y1z1{}, l0_v2{}, l12_v2{}, x3_v2{}, y3_v2{}, z3_v2{};
        pack_line_dbl_v2_input(x1y1z1, q);
        line_dbl_vec_v2(l0_v2, l12_v2, x3_v2, y3_v2, z3_v2, x1y1z1);

        LineEval line_v2{};
        JacobianPoint t_v2{};
        decode_line_dbl_v2(line_v2, t_v2, l0_v2, l12_v2, x3_v2, y3_v2, z3_v2);

        fp2_2x2x2w x1y1{}, z1{}, l0y3_v1{}, l1_v1{}, l2_v1{}, x3_v1{}, z3_v1{};
        pack_line_dbl_v1_inputs(x1y1, z1, q);
        line_dbl_vec_v1(l0y3_v1, l1_v1, l2_v1, x3_v1, z3_v1, x1y1, z1);

        LineEval line_v1{};
        JacobianPoint t_v1{};
        decode_line_dbl_v1(line_v1, t_v1, l0y3_v1, l1_v1, l2_v1, x3_v1, z3_v1);

        expect_fp2_eq(line_v2.c0, line_v1.c0);
        expect_fp2_eq(line_v2.c1, line_v1.c1);
        expect_fp2_eq(line_v2.c2, line_v1.c2);
        expect_fp2_eq(t_v2.X, t_v1.X);
        expect_fp2_eq(t_v2.Y, t_v1.Y);
        expect_fp2_eq(t_v2.Z, t_v1.Z);
    }
}

TEST(G2Avx, LineAddV1MatchesScalar)
{
    InsecureRng rng;

    for (int iter = 0; iter < 112; ++iter) {
        JacobianPoint r{};
        AffineInputPoint q{};
        make_valid_line_add_inputs(r, q, rng, iter);

        LineEval expected_line{};
        JacobianPoint expected_t{};
        line_add_scalar(expected_line, expected_t, r, q);

        fp2_2x2x2w x1y1{}, z1y2{}, x2{}, l0y3{}, l1{}, x3{}, z3{};
        pack_line_add_v1_inputs(x1y1, z1y2, x2, r, q);
        line_add_vec_v1(l0y3, l1, x3, z3, x1y1, z1y2, x2);

        LineEval got_line{};
        JacobianPoint got_t{};
        decode_line_add_v1(got_line, got_t, l0y3, l1, x3, z3);

        expect_fp2_eq(got_line.c0, expected_line.c0);
        expect_fp2_eq(got_line.c1, expected_line.c1);
        expect_fp2_eq(got_line.c2, expected_line.c2);
        expect_fp2_eq(got_t.X, expected_t.X);
        expect_fp2_eq(got_t.Y, expected_t.Y);
        expect_fp2_eq(got_t.Z, expected_t.Z);
    }
}

TEST(G2Avx, LineDblVectorMatchesScalar)
{
    InsecureRng rng;

    for (int iter = 0; iter < 112; ++iter) {
        JacobianPoint q{};
        make_point(q, rng, (iter < 16) ? iter : -1);

        LineEval expected_line{};
        JacobianPoint expected_t{};
        line_dbl_scalar(expected_line, expected_t, q);

        Fp6 got_line{};
        G2Point q_point{};
        G2Point got_t{};
        q_point.X = q.X;
        q_point.Y = q.Y;
        q_point.Z = q.Z;

        line_dbl_vector(got_line, got_t, q_point);

        expect_fp2_eq(got_line[0], expected_line.c0);
        expect_fp2_eq(got_line[1], expected_line.c1);
        expect_fp2_eq(got_line[2], expected_line.c2);
        expect_fp2_eq(got_t.X, expected_t.X);
        expect_fp2_eq(got_t.Y, expected_t.Y);
        expect_fp2_eq(got_t.Z, expected_t.Z);
    }
}

TEST(G2Avx, LineAddVectorMatchesScalar)
{
    InsecureRng rng;

    for (int iter = 0; iter < 112; ++iter) {
        JacobianPoint r{};
        AffineInputPoint q{};
        make_valid_line_add_inputs(r, q, rng, iter);

        LineEval expected_line{};
        JacobianPoint expected_t{};
        line_add_scalar(expected_line, expected_t, r, q);

        Fp6 got_line{};
        G2Point r_point{};
        G2AffinePoint q_point{};
        G2Point got_t{};
        r_point.X = r.X;
        r_point.Y = r.Y;
        r_point.Z = r.Z;
        q_point.X = q.X;
        q_point.Y = q.Y;

        line_add_vector(got_line, got_t, r_point, q_point);

        expect_fp2_eq(got_line[0], expected_line.c0);
        expect_fp2_eq(got_line[1], expected_line.c1);
        expect_fp2_eq(got_line[2], expected_line.c2);
        expect_fp2_eq(got_t.X, expected_t.X);
        expect_fp2_eq(got_t.Y, expected_t.Y);
        expect_fp2_eq(got_t.Z, expected_t.Z);
    }
}

TEST(G2Avx, LineByPx2_4x2x1w)
{
    InsecureRng rng;

    for (int iter = 0; iter < 96; ++iter) {
        JacobianPoint q{};
        make_point(q, rng, (iter < 16) ? iter : -1);

        Fp px{};
        Fp py{};
        if (iter < 16) {
            px = (iter & 1) ? Fp::one() : Fp::zero();
            py = (iter & 2) ? Fp::one() : Fp::zero();
        } else {
            rand_redc_fp(px, rng);
            rand_redc_fp(py, rng);
        }

        LineEval expected_line{};
        JacobianPoint expected_t{};
        line_dbl_scalar(expected_line, expected_t, q);
        scale_line_by_px2_scalar(expected_line, px, py);

        fp2_4x2x1w x1y1z1{}, l0{}, l12{}, x3{}, y3{}, z3{}, px2{};
        pack_line_dbl_v2_input(x1y1z1, q);
        line_dbl_vec_v2(l0, l12, x3, y3, z3, x1y1z1);

        pack_px2_4x2x1w(px2, px, py);
        line_by_Px2_4x2x1w(l12, l12, px2);

        Fp2 got_l1 = fp2_from_8x1w_channels(l12, 4, 5, false);
        Fp2 got_l2 = fp2_from_8x1w_channels(l12, 2, 3, false);
        expect_fp2_eq(got_l1, expected_line.c1);
        expect_fp2_eq(got_l2, expected_line.c2);
    }
}

TEST(G2Avx, LineByPx2_2x2x2w)
{
    InsecureRng rng;

    for (int iter = 0; iter < 96; ++iter) {
        JacobianPoint q{};
        make_point(q, rng, (iter < 16) ? iter : -1);

        Fp px{};
        Fp py{};
        if (iter < 16) {
            px = (iter & 1) ? Fp::one() : Fp::zero();
            py = (iter & 2) ? Fp::one() : Fp::zero();
        } else {
            rand_redc_fp(px, rng);
            rand_redc_fp(py, rng);
        }

        LineEval expected_line{};
        JacobianPoint expected_t{};
        line_dbl_scalar(expected_line, expected_t, q);
        scale_line_by_px2_scalar(expected_line, px, py);

        fp2_2x2x2w x1y1{}, z1{}, l0y3{}, l1{}, l2{}, x3{}, z3{}, l12{}, px2{};
        pack_line_dbl_v1_inputs(x1y1, z1, q);
        line_dbl_vec_v1(l0y3, l1, l2, x3, z3, x1y1, z1);
        blend_0x0F_hl(l12.data(), l2.data(), l1.data());

        pack_px2_2x2x2w(px2, px, py);
        line_by_Px2_2x2x2w(l12, l12, px2);

        Fp2 got_l1 = fp2_from_4x2w_channels(l12, 0, 2, true);
        Fp2 got_l2 = fp2_from_4x2w_channels(l12, 4, 6, true);
        expect_fp2_eq(got_l1, expected_line.c1);
        expect_fp2_eq(got_l2, expected_line.c2);
    }
}

#else

TEST(G2Avx, Disabled)
{
    SUCCEED();
}

#endif
