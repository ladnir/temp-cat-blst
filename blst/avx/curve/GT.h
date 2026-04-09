#pragma once

#include "blst/avx/fields/conversion_inversion.h"
#include "blst/avx/fields/final_exponentiation.h"
#include "blst/avx/fields/fp12.h"
#include "blst/avx/curve/point.h"
#include "blst/avx/gt/glv4.h"
#include "blst/avx/pairing/miller.h"
#include "blst/avx/rns/bounded_ring.h"
#include "blst/scalar/curve/GT.h"

#include <array>
#include <cstddef>
#include <utility>

namespace blst::avx {

namespace avx_gt_detail {

inline constexpr size_t kMillerLoopNMax = 16;

using Ring = BoundedRing<>;
using FP12 = FP12Ring<Ring>;
using PairingFp2 = FP2Ring<Ring>;
using PairingFp2Element = PairingFp2::StandardElement;
using PairingFp2Prepared = FP2N<typename Ring::StandardElement, typename Ring::StandardElement, typename Ring::StandardElement>;
using PairingFp2LeftPrepared = decltype(std::declval<PairingFp2Element>().prep_left(std::declval<const Ring&>()));

inline const Ring& ring()
{
    static const Ring value;
    return value;
}

inline const PairingFp2& fp2ring()
{
    static const PairingFp2 value(ring());
    return value;
}

inline const FP12& fp12ring()
{
    static const FP12 value(ring());
    return value;
}

inline std::array<uint64_t, 4> scalar_to_u64_le(const byte* scalar, size_t nbits)
{
    BLST_ASSERT(nbits <= 256);
    std::array<uint64_t, 4> out = {};
    const size_t nbytes = (nbits + 7) / 8;
    for (size_t i = 0; i < nbytes; ++i) {
        out[i / 8] |= static_cast<uint64_t>(scalar[i]) << ((i % 8) * 8);
    }
    return out;
}

inline blst::Fp12 to_scalar(const typename FP12::StandardElement& value)
{
    return avx_blst_direct_bridge::ring_fp12_to_scalar(ring(), value);
}

inline typename FP12::StandardElement from_scalar(const blst::Fp12& value)
{
    return avx_blst_direct_bridge::scalar_fp12_to_ring(ring(), value);
}

template<class Element>
inline typename FP12::StandardElement to_standard(const Element& value)
{
    return fp12ring().fp12_flat_mul(fp12ring().one(ring()), value, ring());
}

inline PairingFp2Prepared prep_left_helper(const PairingFp2Element& point)
{
    return PairingFp2Prepared(point.x(), point.y(), ring().standard_negate(point.y()));
}

template<class T>
inline PairingFp2Element normalize_fp2(const T& value)
{
    return fp2ring().batch_reduce_expand(value.prep_left(ring()).mul_fp(ring().one()))[0];
}

struct PreparedLine
{
    PairingFp2Element c = {};
    PairingFp2LeftPrepared x = PairingFp2Element().prep_left(ring());
    PairingFp2LeftPrepared y = PairingFp2Element().prep_left(ring());
};

struct alignas(64) PreparedLineTable
{
    PreparedLine lines[68];
};

inline typename FP12::StandardElement make_line_value(
    const PairingFp2Element& c,
    const PairingFp2Element& x,
    const PairingFp2Element& y)
{
    return {
        c.x(),
        c.y(),
        x.x(),
        x.y(),
        ring().zero(),
        ring().zero(),
        ring().zero(),
        ring().zero(),
        y.x(),
        y.y(),
        ring().zero(),
        ring().zero(),
    };
}

inline G1AffinePoint make_px2(const G1AffinePoint& point)
{
    G1AffinePoint out = {};
    out.X = ring().standard_negate(NormalizePointCoord(point.X + point.X, ring()));
    out.Y = NormalizePointCoord(point.Y + point.Y, ring());
    return out;
}

inline auto scale_line_xy(const PreparedLine& in, const G1AffinePoint& p)
{
    return fp2ring().batch_reduce_expand(
        in.x.mul_fp(p.X),
        in.y.mul_fp(p.Y));
}

inline void apply_prepared_line(
    typename FP12::StandardElement& data,
    const PreparedLine& line,
    const G1AffinePoint& p)
{
    const auto scaled = scale_line_xy(line, p);
    fp12ring().flat_mulxy00z0(data, data, line.c, scaled[0], scaled[1], ring());
}

inline void line_add(PreparedLine& line, G2Point& T, const G2Point& R, const G2AffinePoint& Q)
{
    static constexpr auto kThree = std::integral_constant<int, 3>{};

    const PairingFp2Prepared qx = prep_left_helper(Q.X);
    const PairingFp2Prepared qy = prep_left_helper(Q.Y);

    const auto x2x1 = qx * R.X;
    const auto x2y1 = qx * R.Y;
    const auto x2z1 = qx * R.Z;
    const auto y2x1 = qy * R.X;
    const auto y2y1 = qy * R.Y;
    const auto y2z1 = qy * R.Z;

    const auto [x2x1_reduced, x2y1_reduced, x2z1_reduced, y2x1_reduced, y2y1_reduced, y2z1_reduced] =
        fp2ring().batch_reduce_expand(x2x1, x2y1, x2z1, y2x1, y2y1, y2z1);

    line.c = normalize_fp2(x2y1_reduced - y2x1_reduced);
    line.x = normalize_fp2(y2z1_reduced - R.Y).prep_left(ring());
    line.y = normalize_fp2(R.X - x2z1_reduced).prep_left(ring());

    auto x2y1_plus_y2x1 = (x2y1_reduced + y2x1_reduced).prep_left(ring());
    const auto y2z1_plus_y1 = (y2z1_reduced + R.Y).prep(ring());
    const auto x2z1_plus_x1 = x2z1_reduced + R.X;
    const auto three_x2x1 = (x2x1_reduced * kThree).prep(ring());
    const auto b3z1 = fp2ring().mul_3b(R.Z).prep(ring());
    auto y2y1_plus_b3z1 = (y2y1_reduced + b3z1).prep_left(ring());
    const auto y2y1_minus_b3z1 = (y2y1_reduced - b3z1).prep(ring());
    auto b3_x2z1_plus_x1 = fp2ring().mul_3b(x2z1_plus_x1).prep_left(ring());

    const auto t2_19 = x2y1_plus_y2x1 * y2y1_minus_b3z1;
    const auto x3_20 = t2_19.mul_sub_accumulate(b3_x2z1_plus_x1, y2z1_plus_y1, ring());
    const auto y3_21 = b3_x2z1_plus_x1 * three_x2x1;
    const auto y3_23 = y3_21.mul_accumulate(y2y1_plus_b3z1, y2y1_minus_b3z1);
    const auto z3_25 = y2y1_plus_b3z1 * y2z1_plus_y1;
    const auto z3_26 = z3_25.mul_accumulate(x2y1_plus_y2x1, three_x2x1);

    const auto [x3, y3, z3] = fp2ring().batch_reduce_expand(x3_20, y3_23, z3_26);
    T.X = x3;
    T.Y = y3;
    T.Z = z3;
}

inline void line_dbl(PreparedLine& line, G2Point& T, const G2Point& Q)
{
    static constexpr auto kTwo = std::integral_constant<int, 2>{};
    static constexpr auto kThree = std::integral_constant<int, 3>{};
    static constexpr auto kEight = std::integral_constant<int, 8>{};

    const auto px = Q.X.prep_left(ring());
    const auto py = Q.Y.prep_left(ring());
    const auto pz = Q.Z.prep_left(ring());

    const auto yy = py.square();
    const auto zz = pz.square();
    const auto xx = px.square();
    const auto xy = px * py.to_fp2();
    const auto yz = py * pz.to_fp2();

    const auto [xx_reduced, yy_reduced, zz_reduced, xy_reduced, yz_reduced] =
        fp2ring().batch_reduce_expand(xx, yy, zz, xy, yz);

    const auto eight_yy = (yy_reduced * kEight).prep_left(ring());
    const auto b3_zz = fp2ring().mul_3b(zz_reduced).prep(ring());
    const auto yy_plus_b3_zz = (yy_reduced + b3_zz).prep(ring());
    const auto b9_zz = b3_zz * kThree;
    const auto yy_minus_b9_zz = (yy_reduced - b9_zz).prep_left(ring());
    const auto two_xy = (xy_reduced * kTwo).prep(ring());

    line.c = normalize_fp2(b3_zz - yy_reduced);
    line.x = normalize_fp2(xx_reduced * kThree).prep_left(ring());
    line.y = fp2ring().standard_negate(normalize_fp2(yz_reduced * kTwo)).prep_left(ring());

    const auto y3_14 = yy_minus_b9_zz * yy_plus_b3_zz;
    const auto x3_18 = yy_minus_b9_zz * two_xy;
    const auto z3_10 = eight_yy * yz_reduced;
    const auto y3_15 = y3_14.mul_accumulate(eight_yy, b3_zz);

    const auto [x3, y3, z3] = fp2ring().batch_reduce_expand(x3_18, y3_15, z3_10);
    T.X = x3;
    T.Y = y3;
    T.Z = z3;
}

inline void pre_add_n_dbl(PreparedLine lines[], G2Point& t, const G2AffinePoint& q, size_t n)
{
    line_add(lines[0], t, t, q);
    ++lines;
    while (n-- != 0) {
        line_dbl(lines[0], t, t);
        ++lines;
    }
}

inline void apply_add_n_dbl(
    typename FP12::StandardElement& data,
    const PreparedLine lines[],
    const G1AffinePoint& p,
    size_t n)
{
    apply_prepared_line(data, lines[0], p);
    ++lines;

    while (n-- != 0) {
        data = fp12ring().fp12_flat_sqr(data, ring());
        apply_prepared_line(data, lines[0], p);
        ++lines;
    }
}

inline void start_dbl_n(
    typename FP12::StandardElement& data,
    G2Point T[],
    const G1AffinePoint P[],
    size_t n)
{
    BLST_ASSERT(n != 0);

    PreparedLine line = {};
    line_dbl(line, T[0], T[0]);
    const auto first = scale_line_xy(line, P[0]);
    data = make_line_value(line.c, first[0], first[1]);

    for (size_t i = 1; i < n; ++i) {
        line_dbl(line, T[i], T[i]);
        apply_prepared_line(data, line, P[i]);
    }
}

inline void add_n_dbl_n(
    typename FP12::StandardElement& data,
    G2Point T[],
    const G2AffinePoint Q[],
    const G1AffinePoint P[],
    size_t n,
    size_t k)
{
    PreparedLine line = {};

    for (size_t i = 0; i < n; ++i) {
        line_add(line, T[i], T[i], Q[i]);
        apply_prepared_line(data, line, P[i]);
    }

    while (k-- != 0) {
        data = fp12ring().fp12_flat_sqr(data, ring());
        for (size_t i = 0; i < n; ++i) {
            line_dbl(line, T[i], T[i]);
            apply_prepared_line(data, line, P[i]);
        }
    }
}

}  // namespace avx_gt_detail

using PreparedLine = avx_gt_detail::PreparedLine;
using PreparedLineTable = avx_gt_detail::PreparedLineTable;

struct GT
{
    using Field_type = avx_gt_detail::FP12::StandardElement;

    Field_type data = {};

    constexpr GT() = default;
    constexpr GT(const Field_type& value) : data(value) {}

    auto operator[](size_t idx) -> decltype(data[idx]) { return data[idx]; }
    auto operator[](size_t idx) const -> decltype(data[idx]) { return data[idx]; }

    static void add(GT& out, const GT& a, const GT& b)
    {
        out.data = avx_gt_detail::fp12ring().fp12_flat_mul(a.data, b.data, avx_gt_detail::ring());
    }

    static void dadd(GT& out, const GT& a, const GT& b)
    {
        add(out, a, b);
    }

    static void sub(GT& out, const GT& a, const GT& b)
    {
        const auto b_inv = avx_gt_detail::fp12ring().fp12_flat_conjugate(b.data, avx_gt_detail::ring());
        out.data = avx_gt_detail::fp12ring().fp12_flat_mul(a.data, b_inv, avx_gt_detail::ring());
    }

    static void doubling(GT& out, const GT& a)
    {
        out.data = avx_gt_detail::to_standard(
            avx_gt_detail::fp12ring().fp12_flat_cyclotomic_sqr(a.data, avx_gt_detail::ring()));
    }

    static void add_n_dbl(GT& out, const GT& a, size_t n)
    {
        add(out, out, a);
        for (size_t i = 0; i < n; ++i) {
            doubling(out, out);
        }
    }

    static void mult(GT& out, const GT& a, const byte* scalar, size_t nbits)
    {
        const auto words = avx_gt_detail::scalar_to_u64_le(scalar, nbits);
        if (nbits == 256) {
            const auto props = avx_glv4_compile(words.data());
            out.data = avx_cyclotomic_exp_glv4(a.data, props, avx_gt_detail::fp12ring(), avx_gt_detail::ring());
            return;
        }
        out.data = avx_cyclotomic_exp_basic(a.data, words.data(), avx_gt_detail::fp12ring(), avx_gt_detail::ring());
    }

    static void cneg(GT& p, limb_t cbit)
    {
        if (cbit) {
            p.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(p.data, avx_gt_detail::ring());
        }
    }

    static void neg(GT& out, const GT& p)
    {
        out.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(p.data, avx_gt_detail::ring());
    }

    static void cswap(GT& a, GT& b, bool_t cbit)
    {
        vec_cswap(&a, &b, sizeof(GT), cbit);
    }

    static void ccopy(GT& a, const GT& b, bool_t cbit)
    {
        vec_select(&a, &b, &a, sizeof(GT), cbit);
    }

    static limb_t equal(const GT& a, const GT& b)
    {
        return blst::scalar::GT::equal(
            blst::scalar::GT(avx_gt_detail::to_scalar(a.data)),
            blst::scalar::GT(avx_gt_detail::to_scalar(b.data)));
    }

    static const GT& generator()
    {
        static const GT g(identity());
        return g;
    }

    static GT identity()
    {
        return GT(avx_gt_detail::fp12ring().one(avx_gt_detail::ring()));
    }

    static int is_identity(const GT& p)
    {
        return equal(p, identity());
    }

    static bool_t in_group(const GT& p)
    {
        return blst::scalar::GT::in_group(blst::scalar::GT(avx_gt_detail::to_scalar(p.data)));
    }

    static bool_t is_cyclotomic(const GT& p)
    {
        return blst::scalar::GT::is_cyclotomic(blst::scalar::GT(avx_gt_detail::to_scalar(p.data)));
    }

    static void serialize(unsigned char serialized[12 * 48], const GT& p)
    {
        blst::scalar::GT::serialize(serialized, blst::scalar::GT(avx_gt_detail::to_scalar(p.data)));
    }

    static void deserialize(GT& out, const unsigned char serialized[12 * 48])
    {
        blst::scalar::GT scalar = {};
        blst::scalar::GT::deserialize(scalar, serialized);
        out.data = avx_gt_detail::from_scalar(scalar.data);
    }

    static void final_mult(GT& out, const GT& f)
    {
        FP12RNS_BLSTInverter<avx_gt_detail::Ring> inverter;
        out.data = final_exp(f.data, inverter, avx_gt_detail::fp12ring(), avx_gt_detail::ring());
    }

    static void miller_loop(GT& ret, const G2AffinePoint* Q, const G1AffinePoint* P)
    {
        const G2AffinePoint& q = Q ? *Q : G2AffinePoint::generator();
        const G1AffinePoint& p = P ? *P : G1AffinePoint::generator();

        if (G2AffinePoint::is_inf(q) || G1AffinePoint::is_inf(p)) {
            ret = identity();
            return;
        }

        const FP2Ring<avx_gt_detail::Ring> fp2ring(avx_gt_detail::ring());
        const Miller<avx_gt_detail::Ring> miller(
            {p.X, p.Y},
            {q.X, q.Y},
            avx_gt_detail::ring());
        ret.data = miller.miller_loop(
            avx_gt_detail::ring(),
            fp2ring,
            avx_gt_detail::fp12ring());
    }

    static void miller_loop_n(
        GT& out,
        const G2AffinePoint Q[],
        const G1AffinePoint P[],
        size_t n)
    {
        out = identity();
        if (n == 0) {
            return;
        }

        size_t offset = 0;
        while (offset < n) {
            G2Point T[avx_gt_detail::kMillerLoopNMax] = {};
            G2AffinePoint q_batch[avx_gt_detail::kMillerLoopNMax] = {};
            G1AffinePoint p_batch[avx_gt_detail::kMillerLoopNMax] = {};
            size_t batch = 0;

            while (offset < n && batch < avx_gt_detail::kMillerLoopNMax) {
                if (!G2AffinePoint::is_inf(Q[offset]) && !G1AffinePoint::is_inf(P[offset])) {
                    q_batch[batch] = Q[offset];
                    p_batch[batch] = P[offset];
                    T[batch].X = Q[offset].X;
                    T[batch].Y = Q[offset].Y;
                    T[batch].Z = avx_gt_detail::fp2ring().one();
                    ++batch;
                }
                ++offset;
            }

            if (batch == 0) {
                continue;
            }

            typename GT::Field_type chunk = {};
            avx_gt_detail::start_dbl_n(chunk, T, p_batch, batch);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 2);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 3);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 9);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 32);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 16);
            chunk = avx_gt_detail::fp12ring().fp12_flat_conjugate(chunk, avx_gt_detail::ring());
            out.data = avx_gt_detail::fp12ring().fp12_flat_mul(
                out.data,
                chunk,
                avx_gt_detail::ring());
        }
    }

    static void miller_loop_n(
        GT& out,
        const G2AffinePoint* const Qs[],
        const G1AffinePoint* const Ps[],
        size_t n)
    {
        out = identity();
        if (n == 0) {
            return;
        }
        if (Qs == nullptr || Ps == nullptr || Qs[0] == nullptr || Ps[0] == nullptr) {
            return;
        }

        const G2AffinePoint* qptr = nullptr;
        const G1AffinePoint* pptr = nullptr;
        G2Point T[avx_gt_detail::kMillerLoopNMax] = {};
        G2AffinePoint q_batch[avx_gt_detail::kMillerLoopNMax] = {};
        G1AffinePoint p_batch[avx_gt_detail::kMillerLoopNMax] = {};
        size_t batch = 0;

        auto flush_batch = [&]() {
            if (batch == 0) {
                return;
            }
            typename GT::Field_type chunk = {};
            avx_gt_detail::start_dbl_n(chunk, T, p_batch, batch);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 2);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 3);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 9);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 32);
            avx_gt_detail::add_n_dbl_n(chunk, T, q_batch, p_batch, batch, 16);
            chunk = avx_gt_detail::fp12ring().fp12_flat_conjugate(chunk, avx_gt_detail::ring());
            out.data = avx_gt_detail::fp12ring().fp12_flat_mul(
                out.data,
                chunk,
                avx_gt_detail::ring());
            batch = 0;
        };

        for (size_t i = 0; i < n; ++i) {
            qptr = Qs[i] ? Qs[i] : qptr + 1;
            pptr = Ps[i] ? Ps[i] : pptr + 1;
            BLST_ASSERT(qptr != nullptr);
            BLST_ASSERT(pptr != nullptr);

            if (G2AffinePoint::is_inf(*qptr) || G1AffinePoint::is_inf(*pptr)) {
                continue;
            }

            q_batch[batch] = *qptr;
            p_batch[batch] = *pptr;
            T[batch].X = qptr->X;
            T[batch].Y = qptr->Y;
            T[batch].Z = avx_gt_detail::fp2ring().one();
            ++batch;

            if (batch == avx_gt_detail::kMillerLoopNMax) {
                flush_batch();
            }
        }

        flush_batch();
    }

    static void precompute_lines(PreparedLine Qlines[68], const G2AffinePoint& Q)
    {
        G2Point t = {};
        t.X = Q.X;
        t.Y = Q.Y;
        t.Z = avx_gt_detail::fp2ring().one();

        avx_gt_detail::line_dbl(Qlines[0], t, t);
        avx_gt_detail::pre_add_n_dbl(&Qlines[1], t, Q, 2);
        avx_gt_detail::pre_add_n_dbl(&Qlines[4], t, Q, 3);
        avx_gt_detail::pre_add_n_dbl(&Qlines[8], t, Q, 9);
        avx_gt_detail::pre_add_n_dbl(&Qlines[18], t, Q, 32);
        avx_gt_detail::pre_add_n_dbl(&Qlines[51], t, Q, 16);
    }

    static void precompute_lines(blst::scalar::Fp6 Qlines[68], const G2AffinePoint& Q)
    {
        blst::scalar::G2AffinePoint scalar_q = {};
        scalar_q.X = avx_curve_detail::ring_fp2_to_scalar(Q.X);
        scalar_q.Y = avx_curve_detail::ring_fp2_to_scalar(Q.Y);
        blst::scalar::GT::precompute_lines(Qlines, scalar_q);
    }

    static void miller_loop_lines(
        GT& ret,
        const PreparedLine Qlines[68],
        const G1AffinePoint& P)
    {
        const auto first_line = avx_gt_detail::scale_line_xy(Qlines[0], P);
        ret.data = avx_gt_detail::make_line_value(Qlines[0].c, first_line[0], first_line[1]);
        avx_gt_detail::apply_add_n_dbl(ret.data, &Qlines[1], P, 2);
        avx_gt_detail::apply_add_n_dbl(ret.data, &Qlines[4], P, 3);
        avx_gt_detail::apply_add_n_dbl(ret.data, &Qlines[8], P, 9);
        avx_gt_detail::apply_add_n_dbl(ret.data, &Qlines[18], P, 32);
        avx_gt_detail::apply_add_n_dbl(ret.data, &Qlines[51], P, 16);
        ret.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(ret.data, avx_gt_detail::ring());
    }

    static void miller_loop_lines(
        GT& ret,
        const blst::scalar::Fp6 Qlines[68],
        const G1AffinePoint& P)
    {
        blst::scalar::G1AffinePoint scalar_p = {};
        scalar_p.X = avx_curve_detail::ring_fp_to_scalar(P.X);
        scalar_p.Y = avx_curve_detail::ring_fp_to_scalar(P.Y);

        blst::scalar::GT scalar_ret = {};
        blst::scalar::GT::miller_loop_lines(scalar_ret, Qlines, scalar_p);
        ret.data = avx_gt_detail::from_scalar(scalar_ret.data);
    }

    static void mult_by_z(GT& out, const GT& a)
    {
        out.data = avx_gt_detail::to_standard(
            raise_to_z(a.data, avx_gt_detail::fp12ring(), avx_gt_detail::ring()));
    }

    static void mult_by_z_div_by_2(GT& out, const GT& a)
    {
        out.data = avx_gt_detail::to_standard(
            raise_to_z_div_by_2(a.data, avx_gt_detail::fp12ring(), avx_gt_detail::ring()));
    }

    static void frobenius_map(GT& out, const GT& a, size_t n)
    {
        switch (n) {
        case 1:
            out.data = avx_gt_detail::fp12ring().fp12_flat_frobenius_map1(a.data, avx_gt_detail::ring());
            break;
        case 2:
            out.data = avx_gt_detail::fp12ring().fp12_flat_frobenius_map2(a.data, avx_gt_detail::ring());
            break;
        case 3:
            out.data = avx_gt_detail::fp12ring().fp12_flat_frobenius_map3(a.data, avx_gt_detail::ring());
            break;
        default:
            BLST_ASSERT(false && "Unsupported Frobenius power");
        }
    }

    static void conjugate(GT& a)
    {
        a.data = avx_gt_detail::fp12ring().fp12_flat_conjugate(a.data, avx_gt_detail::ring());
    }
};

}  // namespace blst::avx
