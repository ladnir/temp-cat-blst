#include "blst/curve/pairing_avx.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include "blst/curve/G2_avx.h"
#include "blst/fields/Fp12_avx.h"

namespace blst::avx
{
namespace
{
constexpr std::size_t kMillerLoopNMax = 16;

void line_by_px2(Fp6& line, const G1AffinePoint& px2)
{
    Fp::mul(line[1][0], line[1][0], px2.X);
    Fp::mul(line[1][1], line[1][1], px2.X);

    Fp::mul(line[2][0], line[2][0], px2.Y);
    Fp::mul(line[2][1], line[2][1], px2.Y);
}

void start_dbl_n(
    GT& ret,
    G2Point t[],
    const G1AffinePoint px2[],
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

void add_n_dbl_n(
    GT& ret,
    G2Point t[],
    const G2AffinePoint q[],
    const G1AffinePoint px2[],
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

void mul_n_sqr(Fp12& ret, const Fp12& a, std::size_t n)
{
    mul_fp12_vector(ret, ret, a);
    while (n--) {
        cyclotomic_sqr_fp12_vector(ret, ret);
    }
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
    const G2AffinePoint q[],
    const G1AffinePoint p[],
    std::size_t n)
{
    if (n > kMillerLoopNMax) {
        return;
    }

    if ((n == 1) && (vec_is_zero(&q[0], sizeof(q[0])) |
                     vec_is_zero(&p[0], sizeof(p[0])))) {
        vec_copy(&ret.data, get_BLST_Rx().p12, sizeof(Fp12));
        return;
    }

    G2Point t[kMillerLoopNMax];
    G1AffinePoint px2[kMillerLoopNMax];

    for (std::size_t i = 0; i < n; ++i) {
        Fp::add(px2[i].X, p[i].X, p[i].X);
        Fp::neg(px2[i].X, px2[i].X);
        Fp::add(px2[i].Y, p[i].Y, p[i].Y);

        vec_copy(&t[i].X, &q[i].X, 2 * sizeof(t[i].X));
        vec_copy(&t[i].Z, &get_BLST_Rx().p2, sizeof(t[i].Z));
    }

    start_dbl_n(ret, t, px2, n);
    add_n_dbl_n(ret, t, q, px2, n, 2);
    add_n_dbl_n(ret, t, q, px2, n, 3);
    add_n_dbl_n(ret, t, q, px2, n, 9);
    add_n_dbl_n(ret, t, q, px2, n, 32);
    add_n_dbl_n(ret, t, q, px2, n, 16);
    Fp12::conjugate(ret.data);
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
    const G2AffinePoint q[],
    const G1AffinePoint p[],
    std::size_t n)
{
    GT f{};
    miller_loop_n(f, q, p, n);
    final_exp(ret, f);
}

} // namespace blst::avx

#endif
