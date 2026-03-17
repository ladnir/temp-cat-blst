#include "blst/fields/Fp6_avx.h"
#include "blst/backend/avx512/avx_layout.h"
#include "blst/fields/Fp12.h"

/*
 * WARNING: This file is intentionally excluded from the active build during the
 * VROOM backend reset. It is retained only as reference and must not be wired
 * back into active code until the post-reset AVX API is defined.
 */

#include <array>

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
namespace
{
    inline void shuf_01_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        shuf_01(r[0].data(), a[0].data());
        shuf_01(r[1].data(), a[1].data());
    }

    inline void blend_0x11_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
    {
        blend_0x11(r[0].data(), a[0].data(), b[0].data());
        blend_0x11(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x33_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
    {
        blend_0x33(r[0].data(), a[0].data(), b[0].data());
        blend_0x33(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x55_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
    {
        blend_0x55(r[0].data(), a[0].data(), b[0].data());
        blend_0x55(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x77_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
    {
        blend_0x77(r[0].data(), a[0].data(), b[0].data());
        blend_0x77(r[1].data(), a[1].data(), b[1].data());
    }

    inline void perm_0123_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        perm_0123(r[0].data(), a[0].data());
        perm_0123(r[1].data(), a[1].data());
    }

    inline void perm_0300_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        perm_0300(r[0].data(), a[0].data());
        perm_0300(r[1].data(), a[1].data());
    }

    inline void perm_1032_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        perm_1032(r[0].data(), a[0].data());
        perm_1032(r[1].data(), a[1].data());
    }

    inline void perm_1200_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        perm_1200(r[0].data(), a[0].data());
        perm_1200(r[1].data(), a[1].data());
    }

    inline void perm_2000_Fp2Pack8x1x1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a)
    {
        perm_2000(r[0].data(), a[0].data());
        perm_2000(r[1].data(), a[1].data());
    }

    inline void shuf_00_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        shuf_00_dl(r[0].data(), a[0].data());
        shuf_00_dl(r[1].data(), a[1].data());
    }

    inline void shuf_01_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        shuf_01_dl(r[0].data(), a[0].data());
        shuf_01_dl(r[1].data(), a[1].data());
    }

    inline void shuf_11_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        shuf_11_dl(r[0].data(), a[0].data());
        shuf_11_dl(r[1].data(), a[1].data());
    }

    inline void blend_0x33_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
    {
        blend_0x33_dl(r[0].data(), a[0].data(), b[0].data());
        blend_0x33_dl(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x55_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
    {
        blend_0x55_dl(r[0].data(), a[0].data(), b[0].data());
        blend_0x55_dl(r[1].data(), a[1].data(), b[1].data());
    }

    inline void perm_2323_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        perm_2323_dl(r[0].data(), a[0].data());
        perm_2323_dl(r[1].data(), a[1].data());
    }

    inline void perm_3000_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        perm_3000_dl(r[0].data(), a[0].data());
        perm_3000_dl(r[1].data(), a[1].data());
    }
} // namespace

void Fp6Pack8x1::from_aos(Fp6Pack8x1& r, const Fp6* in8)
{
    BLST_ASSERT(in8 != nullptr);

    std::array<Fp2, 8> c0{};
    std::array<Fp2, 8> c1{};
    std::array<Fp2, 8> c2{};

    for (std::size_t i = 0; i < 8; ++i) {
        c0[i] = in8[i][0];
        c1[i] = in8[i][1];
        c2[i] = in8[i][2];
    }

    Fp2Pack8x1x1::from_aos(r.c[0], c0.data());
    Fp2Pack8x1x1::from_aos(r.c[1], c1.data());
    Fp2Pack8x1x1::from_aos(r.c[2], c2.data());
}

void Fp6Pack8x1::to_aos(Fp6* out8, const Fp6Pack8x1& a)
{
    BLST_ASSERT(out8 != nullptr);

    std::array<Fp2, 8> c0{};
    std::array<Fp2, 8> c1{};
    std::array<Fp2, 8> c2{};

    Fp2Pack8x1x1::to_aos(c0.data(), a.c[0]);
    Fp2Pack8x1x1::to_aos(c1.data(), a.c[1]);
    Fp2Pack8x1x1::to_aos(c2.data(), a.c[2]);

    for (std::size_t i = 0; i < 8; ++i) {
        out8[i][0] = c0[i];
        out8[i][1] = c1[i];
        out8[i][2] = c2[i];
    }
}

void Fp6Pack8x1::from_aos_partial(
    Fp6Pack8x1& r,
    const Fp6* in,
    const std::size_t n_active,
    const Fp6& fill)
{
    BLST_ASSERT(n_active <= 8);
    BLST_ASSERT(in != nullptr || n_active == 0);

    std::array<Fp6, 8> tmp{};
    for (std::size_t i = 0; i < n_active; ++i) {
        tmp[i] = in[i];
    }
    for (std::size_t i = n_active; i < 8; ++i) {
        tmp[i] = fill;
    }

    from_aos(r, tmp.data());
}

void Fp6Pack8x1::to_aos_partial(Fp6* out, const std::size_t n_active, const Fp6Pack8x1& a)
{
    BLST_ASSERT(n_active <= 8);
    BLST_ASSERT(out != nullptr || n_active == 0);

    std::array<Fp6, 8> tmp{};
    to_aos(tmp.data(), a);
    for (std::size_t i = 0; i < n_active; ++i) {
        out[i] = tmp[i];
    }
}

void Fp6Pack8x1::add(Fp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    Fp2Pack8x1x1::add(r.c[0], a.c[0], b.c[0]);
    Fp2Pack8x1x1::add(r.c[1], a.c[1], b.c[1]);
    Fp2Pack8x1x1::add(r.c[2], a.c[2], b.c[2]);
}

void Fp6Pack8x1::sub(Fp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    Fp2Pack8x1x1::sub(r.c[0], a.c[0], b.c[0]);
    Fp2Pack8x1x1::sub(r.c[1], a.c[1], b.c[1]);
    Fp2Pack8x1x1::sub(r.c[2], a.c[2], b.c[2]);
}

void Fp6Pack8x1::mul(UFp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    Fp2Pack8x1x1 aa{}, bb{};
    UFp2Pack8x1x1 t0{}, t1{}, t2{}, t3{};

    UFp2Pack8x1x1::mul(t0, a.c[0], b.c[0]);
    UFp2Pack8x1x1::mul(t1, a.c[1], b.c[1]);
    UFp2Pack8x1x1::mul(t2, a.c[2], b.c[2]);

    // r.c[0] = ((a1+a2)*(b1+b2) - a1*b1 - a2*b2)*(u+1) + a0*b0
    Fp2Pack8x1x1::add(aa, a.c[1], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[1], b.c[2]);
    UFp2Pack8x1x1::mul(r.c[0], aa, bb);
    UFp2Pack8x1x1::sub(r.c[0], r.c[0], t1);
    UFp2Pack8x1x1::sub(r.c[0], r.c[0], t2);
    UFp2Pack8x1x1::mul_by_u_plus_1(t3, r.c[0]);
    UFp2Pack8x1x1::add(r.c[0], t3, t0);

    // r.c[1] = (a0+a1)*(b0+b1) - a0*b0 - a1*b1 + (u+1)*a2*b2
    Fp2Pack8x1x1::add(aa, a.c[0], a.c[1]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[1]);
    UFp2Pack8x1x1::mul(r.c[1], aa, bb);
    UFp2Pack8x1x1::sub(r.c[1], r.c[1], t0);
    UFp2Pack8x1x1::sub(r.c[1], r.c[1], t1);
    UFp2Pack8x1x1::mul_by_u_plus_1(t3, t2);
    UFp2Pack8x1x1::add(r.c[1], r.c[1], t3);

    // r.c[2] = (a0+a2)*(b0+b2) - a0*b0 - a2*b2 + a1*b1
    Fp2Pack8x1x1::add(aa, a.c[0], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[2]);
    UFp2Pack8x1x1::mul(r.c[2], aa, bb);
    UFp2Pack8x1x1::sub(r.c[2], r.c[2], t0);
    UFp2Pack8x1x1::sub(r.c[2], r.c[2], t2);
    UFp2Pack8x1x1::add(r.c[2], r.c[2], t1);
}

void Fp6Pack8x1::mul(Fp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    UFp6Pack8x1 wide{};
    Fp6Pack8x1::mul(wide, a, b);
    UFp6Pack8x1::redc(r, wide);
}

void Fp6Pack8x1::sqr(UFp6Pack8x1& r, const Fp6Pack8x1& a)
{
    Fp2Pack8x1x1 t0{};
    UFp2Pack8x1x1 s0{}, m01{}, m12{}, s2{}, rx{};

    UFp2Pack8x1x1::mul(s0, a.c[0], a.c[0]);

    UFp2Pack8x1x1::mul(m01, a.c[0], a.c[1]);
    UFp2Pack8x1x1::add(m01, m01, m01);

    UFp2Pack8x1x1::mul(m12, a.c[1], a.c[2]);
    UFp2Pack8x1x1::add(m12, m12, m12);

    UFp2Pack8x1x1::mul(s2, a.c[2], a.c[2]);

    // r.c[2] = (a0+a1+a2)^2 - a0^2 - a2^2 - 2*(a0*a1) - 2*(a1*a2)
    Fp2Pack8x1x1::add(t0, a.c[2], a.c[1]);
    Fp2Pack8x1x1::add(t0, t0, a.c[0]);
    UFp2Pack8x1x1::mul(rx, t0, t0);
    UFp2Pack8x1x1::sub(r.c[2], rx, s0);
    UFp2Pack8x1x1::sub(r.c[2], r.c[2], s2);
    UFp2Pack8x1x1::sub(r.c[2], r.c[2], m01);
    UFp2Pack8x1x1::sub(r.c[2], r.c[2], m12);

    // r.c[0] = a0^2 + 2*(a1*a2)*(u+1)
    UFp2Pack8x1x1::mul_by_u_plus_1(rx, m12);
    UFp2Pack8x1x1::add(r.c[0], rx, s0);

    // r.c[1] = a2^2*(u+1) + 2*(a0*a1)
    UFp2Pack8x1x1::mul_by_u_plus_1(rx, s2);
    UFp2Pack8x1x1::add(r.c[1], rx, m01);
}

void Fp6Pack8x1::sqr(Fp6Pack8x1& r, const Fp6Pack8x1& a)
{
    UFp6Pack8x1 wide{};
    Fp6Pack8x1::sqr(wide, a);
    UFp6Pack8x1::redc(r, wide);
}

void Fp6Pack8x1::mul_by_nonresidue(Fp6Pack8x1& r, const Fp6Pack8x1& a)
{
    Fp2Pack8x1x1 t0 = a.c[0];
    Fp2Pack8x1x1 t1 = a.c[1];
    Fp2Pack8x1x1 t2{};
    Fp2Pack8x1x1::mul_by_u_plus_1(t2, a.c[2]);
    r.c[0] = t2;
    r.c[1] = t0;
    r.c[2] = t1;
}

void Fp6Pack8x1::select(Fp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b, const u8 lane_mask)
{
    Fp2Pack8x1x1::select(r.c[0], a.c[0], b.c[0], lane_mask);
    Fp2Pack8x1x1::select(r.c[1], a.c[1], b.c[1], lane_mask);
    Fp2Pack8x1x1::select(r.c[2], a.c[2], b.c[2], lane_mask);
}

void UFp6Pack8x1::redc(Fp6Pack8x1& r, const UFp6Pack8x1& a)
{
    UFp2Pack8x1x1::redc(r.c[0], a.c[0]);
    UFp2Pack8x1x1::redc(r.c[1], a.c[1]);
    UFp2Pack8x1x1::redc(r.c[2], a.c[2]);
}

void UFp6Pack8x1::add(UFp6Pack8x1& r, const UFp6Pack8x1& a, const UFp6Pack8x1& b)
{
    UFp2Pack8x1x1::add(r.c[0], a.c[0], b.c[0]);
    UFp2Pack8x1x1::add(r.c[1], a.c[1], b.c[1]);
    UFp2Pack8x1x1::add(r.c[2], a.c[2], b.c[2]);
}

void UFp6Pack8x1::sub(UFp6Pack8x1& r, const UFp6Pack8x1& a, const UFp6Pack8x1& b)
{
    UFp2Pack8x1x1::sub(r.c[0], a.c[0], b.c[0]);
    UFp2Pack8x1x1::sub(r.c[1], a.c[1], b.c[1]);
    UFp2Pack8x1x1::sub(r.c[2], a.c[2], b.c[2]);
}

void UFp6Pack8x1::mul_by_nonresidue(UFp6Pack8x1& r, const UFp6Pack8x1& a)
{
    UFp2Pack8x1x1 t0 = a.c[0];
    UFp2Pack8x1x1 t1 = a.c[1];
    UFp2Pack8x1x1 t2{};
    UFp2Pack8x1x1::mul_by_u_plus_1(t2, a.c[2]);
    r.c[0] = t2;
    r.c[1] = t0;
    r.c[2] = t1;
}

void UFp6Pack8x1::select(UFp6Pack8x1& r, const UFp6Pack8x1& a, const UFp6Pack8x1& b, const u8 lane_mask)
{
    UFp2Pack8x1x1::select(r.c[0], a.c[0], b.c[0], lane_mask);
    UFp2Pack8x1x1::select(r.c[1], a.c[1], b.c[1], lane_mask);
    UFp2Pack8x1x1::select(r.c[2], a.c[2], b.c[2], lane_mask);
}

void Fp6Pack4x2x1x1::mul(
    UFp2Pack8x1x1& r01,
    UFp2Pack8x1x1& r2,
    const Fp2Pack8x1x1& ab0,
    const Fp2Pack8x1x1& ab1,
    const Fp2Pack8x1x1& ab2)
{
    Fp2Pack8x1x1 t0, t1, t2, t3, t4;
    UFp2Pack8x1x1 tt0, tt1, tt2, tt3, tt4, tt5;

    Fp2Pack8x1x1::add(t0, ab1, ab2);
    Fp2Pack8x1x1::add(t1, ab0, ab1);
    Fp2Pack8x1x1::add(t2, ab0, ab2);
    shuf_01_Fp2Pack8x1x1(t3, ab1);
    blend_0x55_Fp2Pack8x1x1(t3, t3, ab0);
    shuf_01_Fp2Pack8x1x1(t4, ab0);
    blend_0x55_Fp2Pack8x1x1(t4, ab1, t4);
    UFp2Pack8x1x1::mul(tt0, t3, t4);
    shuf_01_Fp2Pack8x1x1(t3, ab2);
    blend_0x55_Fp2Pack8x1x1(t3, t3, t0);
    shuf_01_Fp2Pack8x1x1(t4, t0);
    blend_0x55_Fp2Pack8x1x1(t4, ab2, t4);
    UFp2Pack8x1x1::mul(tt1, t3, t4);
    shuf_01_Fp2Pack8x1x1(t3, t1);
    blend_0x55_Fp2Pack8x1x1(t3, t3, t2);
    shuf_01_Fp2Pack8x1x1(t4, t2);
    blend_0x55_Fp2Pack8x1x1(t4, t1, t4);
    UFp2Pack8x1x1::mul(tt2, t3, t4);
    blend_0x55_UFp2Pack8x1x1(tt3, tt2, tt1);
    shuf_11_UFp2Pack8x1x1(tt4, tt0);
    UFp2Pack8x1x1::sub(tt3, tt3, tt4);
    shuf_01_UFp2Pack8x1x1(tt2, tt2);
    blend_0x55_UFp2Pack8x1x1(tt2, tt2, tt3);
    shuf_11_UFp2Pack8x1x1(tt1, tt1);
    UFp2Pack8x1x1::sub(tt2, tt2, tt1);
    shuf_11_UFp2Pack8x1x1(tt5, tt2);
    blend_0x55_UFp2Pack8x1x1(tt5, tt3, tt5);
    shuf_00_UFp2Pack8x1x1(tt0, tt0);
    UFp2Pack8x1x1::sub(tt5, tt5, tt0);
    blend_0x55_UFp2Pack8x1x1(tt2, tt1, tt2);
    UFp2Pack8x1x1::mul_by_u_plus_1(tt1, tt2);
    blend_0x55_UFp2Pack8x1x1(tt0, tt5, tt0);
    UFp2Pack8x1x1::add(r01, tt0, tt1);
    UFp2Pack8x1x1::add(r2, tt5, tt4);
}

void Fp6Pack2x2x2x1::mul(
    UFp2Pack4x2x1& r01,
    UFp2Pack4x2x1& r2,
    const Fp2Pack4x2x1& ab0,
    const Fp2Pack4x2x1& ab1,
    const Fp2Pack4x2x1& ab2)
{
    Fp2Pack4x2x1 t0, t1, t2, t3, t4;
    UFp2Pack4x2x1 tt0, tt1, tt2, tt3, tt4, tt5;

    Fp2Pack4x2x1::add(t0, ab1, ab2);
    Fp2Pack4x2x1::add(t1, ab0, ab1);
    Fp2Pack4x2x1::add(t2, ab0, ab2);
    perm_1032(t3.data(), ab1.data());
    blend_0x33(t3.data(), t3.data(), ab0.data());
    perm_1032(t4.data(), ab0.data());
    blend_0x33(t4.data(), ab1.data(), t4.data());
    UFp2Pack4x2x1::mul(tt0, t3, t4);
    perm_1032(t3.data(), ab2.data());
    blend_0x33(t3.data(), t3.data(), t0.data());
    perm_1032(t4.data(), t0.data());
    blend_0x33(t4.data(), ab2.data(), t4.data());
    UFp2Pack4x2x1::mul(tt1, t3, t4);
    perm_1032(t3.data(), t1.data());
    blend_0x33(t3.data(), t3.data(), t2.data());
    perm_1032(t4.data(), t2.data());
    blend_0x33(t4.data(), t1.data(), t4.data());
    UFp2Pack4x2x1::mul(tt2, t3, t4);
    blend_0x33_dl(tt3.data(), tt2.data(), tt1.data());
    perm_3232_dl(tt4.data(), tt0.data());
    UFp2Pack4x2x1::sub(tt3, tt3, tt4);
    perm_1032_dl(tt2.data(), tt2.data());
    blend_0x33_dl(tt2.data(), tt2.data(), tt3.data());
    perm_3232_dl(tt1.data(), tt1.data());
    UFp2Pack4x2x1::sub(tt2, tt2, tt1);
    perm_3232_dl(tt5.data(), tt2.data());
    blend_0x33_dl(tt5.data(), tt3.data(), tt5.data());
    perm_1010_dl(tt0.data(), tt0.data());
    UFp2Pack4x2x1::sub(tt5, tt5, tt0);
    blend_0x33_dl(tt2.data(), tt1.data(), tt2.data());
    UFp2Pack4x2x1::mul_by_u_plus_1(tt1, tt2);
    blend_0x33_dl(tt0.data(), tt5.data(), tt0.data());
    UFp2Pack4x2x1::add(r01, tt0, tt1);
    UFp2Pack4x2x1::add(r2, tt5, tt4);
}

void Fp6Pack1x2x2x2::mul(
    UFp2Pack2x2x2& r01,
    UFp2Pack2x2x2& r2,
    const Fp2Pack2x2x2& ab0,
    const Fp2Pack2x2x2& ab1,
    const Fp2Pack2x2x2& ab2)
{
    Fp2Pack2x2x2 t0, t1, t2, t3, t4;
    UFp2Pack2x2x2 tt0, tt1, tt2, tt3, tt4, tt5;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(7, 6, 5, 4, 7, 6, 5, 4);
    const cc512i m2 = VSET(3, 2, 1, 0, 3, 2, 1, 0);

    Fp2Pack2x2x2::add(t0, ab1, ab2);
    Fp2Pack2x2x2::add(t1, ab0, ab1);
    Fp2Pack2x2x2::add(t2, ab0, ab2);
    perm_var_hl(t3.data(), ab1.data(), m0);
    blend_0x0F_hl(t3.data(), t3.data(), ab0.data());
    perm_var_hl(t4.data(), ab0.data(), m0);
    blend_0x0F_hl(t4.data(), ab1.data(), t4.data());
    UFp2Pack2x2x2::mul(tt0, t3, t4);
    perm_var_hl(t3.data(), ab2.data(), m0);
    blend_0x0F_hl(t3.data(), t3.data(), t0.data());
    perm_var_hl(t4.data(), t0.data(), m0);
    blend_0x0F_hl(t4.data(), ab2.data(), t4.data());
    UFp2Pack2x2x2::mul(tt1, t3, t4);
    perm_var_hl(t3.data(), t1.data(), m0);
    blend_0x0F_hl(t3.data(), t3.data(), t2.data());
    perm_var_hl(t4.data(), t2.data(), m0);
    blend_0x0F_hl(t4.data(), t1.data(), t4.data());
    UFp2Pack2x2x2::mul(tt2, t3, t4);
    blend_0x0F_vl(tt3.data(), tt2.data(), tt1.data());
    perm_var_vl(tt4.data(), tt0.data(), m1);
    UFp2Pack2x2x2::sub(tt3, tt3, tt4);
    perm_var_vl(tt2.data(), tt2.data(), m0);
    blend_0x0F_vl(tt2.data(), tt2.data(), tt3.data());
    perm_var_vl(tt1.data(), tt1.data(), m1);
    UFp2Pack2x2x2::sub(tt2, tt2, tt1);
    perm_var_vl(tt5.data(), tt2.data(), m1);
    blend_0x0F_vl(tt5.data(), tt3.data(), tt5.data());
    perm_var_vl(tt0.data(), tt0.data(), m2);
    UFp2Pack2x2x2::sub(tt5, tt5, tt0);
    blend_0x0F_vl(tt2.data(), tt1.data(), tt2.data());
    UFp2Pack2x2x2::mul_by_u_plus_1(tt1, tt2);
    blend_0x0F_vl(tt0.data(), tt5.data(), tt0.data());
    UFp2Pack2x2x2::add(r01, tt0, tt1);
    UFp2Pack2x2x2::add(r2, tt5, tt4);
}

void Fp6Pack2x2x2x1::mul_by_xy00z0(
    UFp2Pack4x2x1& r01,
    UFp2Pack4x2x1& r2,
    UFp2Pack4x2x1& r3,
    UFp2Pack4x2x1& r45,
    const Fp2Pack4x2x1& a01,
    const Fp2Pack4x2x1& a2,
    const Fp2Pack4x2x1& b01,
    const Fp2Pack4x2x1& b4)
{
    Fp2Pack4x2x1 t0, t1;
    UFp2Pack4x2x1 tt0, tt1, tt2;

    perm_1032(t0.data(), b01.data());
    blend_0x33(t1.data(), b4.data(), t0.data());
    UFp2Pack4x2x1::mul(tt0, t1, a2);
    UFp2Pack4x2x1::mul_by_u_plus_1(r3, tt0);
    UFp2Pack4x2x1::mul(r45, a01, b4);
    UFp2Pack4x2x1::mul(tt0, a01, b01);
    blend_0x33(t0.data(), t0.data(), a01.data());
    perm_1032(t1.data(), a01.data());
    blend_0x33(t1.data(), b01.data(), t1.data());
    Fp2Pack4x2x1::add(t0, t0, t1);
    blend_0x33(t1.data(), t0.data(), b01.data());
    perm_1032(t0.data(), t0.data());
    blend_0x33(t0.data(), t0.data(), a2.data());
    UFp2Pack4x2x1::mul(tt1, t1, t0);
    perm_1032_dl(tt2.data(), tt0.data());
    UFp2Pack4x2x1::sax2(r2, tt1, tt2);
    blend_0x33_dl(tt1.data(), r2.data(), r3.data());
    UFp2Pack4x2x1::sax2(r01, tt1, tt0);
}

void Fp6Pack4x2x1x1::mul_by_xy00z0(
    UFp2Pack8x1x1& r01,
    UFp2Pack8x1x1& r23,
    UFp2Pack8x1x1& r45,
    const Fp2Pack8x1x1& a,
    const Fp2Pack8x1x1& b)
{
    Fp2Pack8x1x1 t0, t1;
    UFp2Pack8x1x1 tt0, tt1, tt2;

    perm_2000_Fp2Pack8x1x1(t0, b);
    UFp2Pack8x1x1::mul(r45, a, t0);
    blend_0x33_Fp2Pack8x1x1(t0, b, a);
    perm_0300_Fp2Pack8x1x1(t1, t0);
    perm_1200_Fp2Pack8x1x1(t0, t0);
    Fp2Pack8x1x1::add(t0, t0, t1);
    blend_0x77_Fp2Pack8x1x1(t1, t0, a);
    perm_0123_Fp2Pack8x1x1(t1, t1);
    perm_1032_Fp2Pack8x1x1(t0, t0);
    blend_0x11_Fp2Pack8x1x1(t0, b, t0);
    UFp2Pack8x1x1::mul(tt2, t0, t1);
    blend_0x33_UFp2Pack8x1x1(tt1, r45, tt2);
    perm_2323_UFp2Pack8x1x1(tt0, tt2);
    blend_0x33_UFp2Pack8x1x1(tt0, tt1, tt0);
    blend_0x33_dl(tt1[1].data(), tt1[0].data(), tt1[1].data());
    blend_0x33_dl(tt0[0].data(), tt0[1].data(), tt0[0].data());
    UFpPack8x1::ssasx2(r23[0], tt1[0], tt0[0]);
    UFpPack8x1::aaasx2(r23[1], tt1[1], tt0[1]);
    perm_3000_UFp2Pack8x1x1(tt0, r23);
    UFp2Pack8x1x1::asx4(r01, tt0, tt2);
}

} // namespace blst

#endif
