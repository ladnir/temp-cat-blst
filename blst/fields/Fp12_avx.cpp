#include "blst/fields/Fp12_avx.h"
#include "blst/backend/avx512/avx_layout.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include <array>

namespace blst
{
namespace
{
    BLST_FORCEINLINE void sqr_fp4_pack(
        Fp2Pack8x1x1& r0,
        Fp2Pack8x1x1& r1,
        const Fp2Pack8x1x1& a0,
        const Fp2Pack8x1x1& a1)
    {
        Fp2Pack8x1x1 t0{};
        Fp2Pack8x1x1 t1{};
        Fp2Pack8x1x1 t2{};

        Fp2Pack8x1x1::sqr(t0, a0);
        Fp2Pack8x1x1::sqr(t1, a1);

        Fp2Pack8x1x1::add(r1, a0, a1);
        Fp2Pack8x1x1::mul_by_u_plus_1(r0, t1);
        Fp2Pack8x1x1::add(r0, r0, t0);

        Fp2Pack8x1x1::sqr(t2, r1);
        Fp2Pack8x1x1::sub(t2, t2, t0);
        Fp2Pack8x1x1::sub(r1, t2, t1);
    }

    BLST_FORCEINLINE Fp2 make_fp2_real(const Fp& real)
    {
        Fp2 out{};
        out[0] = real;
        out[1] = Fp::zero();
        return out;
    }

    BLST_FORCEINLINE const std::array<Fp2Pack8x1x1, 4>& fp6_c1_coeff_packs()
    {
        static const std::array<Fp2Pack8x1x1, 4> kPacks = []() {
            const std::array<Fp2, 4> kCoeff = {
                Fp2::one(),
                Fp2{Fp{0},
                    Fp{TO_LIMB_T(0xcd03c9e48671f071), TO_LIMB_T(0x5dab22461fcda5d2),
                       TO_LIMB_T(0x587042afd3851b95), TO_LIMB_T(0x8eb60ebe01bacb9e),
                       TO_LIMB_T(0x03f97d6e83d050d2), TO_LIMB_T(0x18f0206554638741)}},
                Fp2{Fp{TO_LIMB_T(0x30f1361b798a64e8), TO_LIMB_T(0xf3b8ddab7ece5a2a),
                       TO_LIMB_T(0x16a8ca3ac61577f7), TO_LIMB_T(0xc26a2ff874fd029b),
                       TO_LIMB_T(0x3636b76660701c6e), TO_LIMB_T(0x051ba4ab241b6160)}},
                Fp2{Fp{0}, Fp{ONE_MONT_P_LIMBS}},
            };

            std::array<Fp2Pack8x1x1, 4> packs{};
            std::array<Fp2, 8> lanes{};
            for (std::size_t i = 0; i < packs.size(); ++i) {
                lanes.fill(kCoeff[i]);
                Fp2Pack8x1x1::from_aos(packs[i], lanes.data());
            }
            return packs;
        }();
        return kPacks;
    }

    BLST_FORCEINLINE const std::array<Fp2Pack8x1x1, 4>& fp6_c2_coeff_packs()
    {
        static const std::array<Fp2Pack8x1x1, 4> kPacks = []() {
            const std::array<Fp2, 4> kCoeff = {
                Fp2::one(),
                make_fp2_real(Fp{TO_LIMB_T(0x890dc9e4867545c3), TO_LIMB_T(0x2af322533285a5d5),
                                 TO_LIMB_T(0x50880866309b7e2c), TO_LIMB_T(0xa20d1b8c7e881024),
                                 TO_LIMB_T(0x14e4f04fe2db9068), TO_LIMB_T(0x14e56d3f1564853a)}),
                make_fp2_real(Fp{TO_LIMB_T(0xcd03c9e48671f071), TO_LIMB_T(0x5dab22461fcda5d2),
                                 TO_LIMB_T(0x587042afd3851b95), TO_LIMB_T(0x8eb60ebe01bacb9e),
                                 TO_LIMB_T(0x03f97d6e83d050d2), TO_LIMB_T(0x18f0206554638741)}),
                make_fp2_real(Fp{TO_LIMB_T(0x43f5fffffffcaaae), TO_LIMB_T(0x32b7fff2ed47fffd),
                                 TO_LIMB_T(0x07e83a49a2e99d69), TO_LIMB_T(0xeca8f3318332bb7a),
                                 TO_LIMB_T(0xef148d1ea0f4c069), TO_LIMB_T(0x040ab3263eff0206)}),
            };

            std::array<Fp2Pack8x1x1, 4> packs{};
            std::array<Fp2, 8> lanes{};
            for (std::size_t i = 0; i < packs.size(); ++i) {
                lanes.fill(kCoeff[i]);
                Fp2Pack8x1x1::from_aos(packs[i], lanes.data());
            }
            return packs;
        }();
        return kPacks;
    }

    BLST_FORCEINLINE const std::array<Fp2Pack8x1x1, 4>& fp12_c1_coeff_packs()
    {
        static const std::array<Fp2Pack8x1x1, 4> kPacks = []() {
            const std::array<Fp2, 4> kCoeff = {
                Fp2::one(),
                Fp2{
                    Fp{TO_LIMB_T(0x07089552b319d465), TO_LIMB_T(0xc6695f92b50a8313),
                       TO_LIMB_T(0x97e83cccd117228f), TO_LIMB_T(0xa35baecab2dc29ee),
                       TO_LIMB_T(0x1ce393ea5daace4d), TO_LIMB_T(0x08f2220fb0fb66eb)},
                    Fp{TO_LIMB_T(0xb2f66aad4ce5d646), TO_LIMB_T(0x5842a06bfc497cec),
                       TO_LIMB_T(0xcf4895d42599d394), TO_LIMB_T(0xc11b9cba40a8e8d0),
                       TO_LIMB_T(0x2e3813cbe5a0de89), TO_LIMB_T(0x110eefda88847faf)}},
                Fp2{Fp{TO_LIMB_T(0xecfb361b798dba3a), TO_LIMB_T(0xc100ddb891865a2c),
                       TO_LIMB_T(0x0ec08ff1232bda8e), TO_LIMB_T(0xd5c13cc6f1ca4721),
                       TO_LIMB_T(0x47222a47bf7b5c04), TO_LIMB_T(0x0110f184e51c5f59)}},
                Fp2{
                    Fp{TO_LIMB_T(0x3e2f585da55c9ad1), TO_LIMB_T(0x4294213d86c18183),
                       TO_LIMB_T(0x382844c88b623732), TO_LIMB_T(0x92ad2afd19103e18),
                       TO_LIMB_T(0x1d794e4fac7cf0b9), TO_LIMB_T(0x0bd592fc7d825ec8)},
                    Fp{TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),
                       TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),
                       TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2)}},
            };

            std::array<Fp2Pack8x1x1, 4> packs{};
            std::array<Fp2, 8> lanes{};
            for (std::size_t i = 0; i < packs.size(); ++i) {
                lanes.fill(kCoeff[i]);
                Fp2Pack8x1x1::from_aos(packs[i], lanes.data());
            }
            return packs;
        }();
        return kPacks;
    }

    BLST_FORCEINLINE void frobenius_map_fp2_pack(
        Fp2Pack8x1x1& r,
        const Fp2Pack8x1x1& a,
        const std::size_t n_mod4)
    {
        if ((n_mod4 & 1u) == 0) {
            r = a;
            return;
        }
        Fp2Pack8x1x1::conjugate(r, a);
    }

    BLST_FORCEINLINE void frobenius_map_fp6_pack(
        Fp6Pack8x1& r,
        const Fp6Pack8x1& a,
        const std::size_t n_mod4)
    {
        if (n_mod4 == 0) {
            r = a;
            return;
        }

        frobenius_map_fp2_pack(r.c[0], a.c[0], n_mod4);
        frobenius_map_fp2_pack(r.c[1], a.c[1], n_mod4);
        frobenius_map_fp2_pack(r.c[2], a.c[2], n_mod4);

        Fp2Pack8x1x1 t1{};
        Fp2Pack8x1x1 t2{};
        Fp2Pack8x1x1::mul(t1, r.c[1], fp6_c1_coeff_packs()[n_mod4]);
        Fp2Pack8x1x1::mul(t2, r.c[2], fp6_c2_coeff_packs()[n_mod4]);
        r.c[1] = t1;
        r.c[2] = t2;
    }

    inline void shuf_01_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        shuf_01_dl(r[0].data(), a[0].data());
        shuf_01_dl(r[1].data(), a[1].data());
    }

    inline void blend_0x01_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
    {
        blend_0x01_dl(r[0].data(), a[0].data(), b[0].data());
        blend_0x01_dl(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x03_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
    {
        blend_0x03_dl(r[0].data(), a[0].data(), b[0].data());
        blend_0x03_dl(r[1].data(), a[1].data(), b[1].data());
    }

    inline void blend_0x33_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
    {
        blend_0x33_dl(r[0].data(), a[0].data(), b[0].data());
        blend_0x33_dl(r[1].data(), a[1].data(), b[1].data());
    }

    inline void perm_1032_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        perm_1032_dl(r[0].data(), a[0].data());
        perm_1032_dl(r[1].data(), a[1].data());
    }

    inline void perm_var_UFp2Pack8x1x1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const cc512i idx)
    {
        perm_var_dl(r[0].data(), a[0].data(), idx);
        perm_var_dl(r[1].data(), a[1].data(), idx);
    }

    inline void conv_sltohl(cc512i* r, const cc512i* a)
    {
        blend_0x55_hl(r, &a[VWORDS], a);
    }

    template<class RPackT, class APackT>
    inline void conv_sltohl(RPackT& r, const APackT& a)
    {
        conv_sltohl(r.data(), a.data());
    }

    inline void conv_dltovl(cc512i* r, const cc512i* a)
    {
        const cc512i m0 = VSET(0, FMASK, 0, FMASK, 0, FMASK, 0, FMASK);
        for (u64 i = 0; i < 2 * VWORDS; ++i) {
            r[i] = VAND(a[i], m0);
        }
        blend_0x55_hl(&r[2 * VWORDS], &a[3 * VWORDS], &a[2 * VWORDS]);
    }

    template<class RPackT, class APackT>
    inline void conv_dltovl(RPackT& r, const APackT& a)
    {
        conv_dltovl(r.data(), a.data());
    }

} // namespace

void Fp12Pack8x1::from_aos(Fp12Pack8x1& r, const Fp12* in8)
{
    BLST_ASSERT(in8 != nullptr);

    std::array<Fp6, 8> c0{};
    std::array<Fp6, 8> c1{};

    for (std::size_t i = 0; i < 8; ++i) {
        c0[i] = in8[i][0];
        c1[i] = in8[i][1];
    }

    Fp6Pack8x1::from_aos(r.c[0], c0.data());
    Fp6Pack8x1::from_aos(r.c[1], c1.data());
}

void Fp12Pack8x1::to_aos(Fp12* out8, const Fp12Pack8x1& a)
{
    BLST_ASSERT(out8 != nullptr);

    std::array<Fp6, 8> c0{};
    std::array<Fp6, 8> c1{};

    Fp6Pack8x1::to_aos(c0.data(), a.c[0]);
    Fp6Pack8x1::to_aos(c1.data(), a.c[1]);

    for (std::size_t i = 0; i < 8; ++i) {
        out8[i][0] = c0[i];
        out8[i][1] = c1[i];
    }
}

void Fp12Pack8x1::from_aos_partial(
    Fp12Pack8x1& r,
    const Fp12* in,
    const std::size_t n_active,
    const Fp12& fill)
{
    BLST_ASSERT(n_active <= 8);
    BLST_ASSERT(in != nullptr || n_active == 0);

    std::array<Fp12, 8> tmp{};
    for (std::size_t i = 0; i < n_active; ++i) {
        tmp[i] = in[i];
    }
    for (std::size_t i = n_active; i < 8; ++i) {
        tmp[i] = fill;
    }

    from_aos(r, tmp.data());
}

void Fp12Pack8x1::to_aos_partial(Fp12* out, const std::size_t n_active, const Fp12Pack8x1& a)
{
    BLST_ASSERT(n_active <= 8);
    BLST_ASSERT(out != nullptr || n_active == 0);

    std::array<Fp12, 8> tmp{};
    to_aos(tmp.data(), a);
    for (std::size_t i = 0; i < n_active; ++i) {
        out[i] = tmp[i];
    }
}

void Fp12Pack8x1::add(Fp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b)
{
    Fp6Pack8x1::add(r.c[0], a.c[0], b.c[0]);
    Fp6Pack8x1::add(r.c[1], a.c[1], b.c[1]);
}

void Fp12Pack8x1::sub(Fp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b)
{
    Fp6Pack8x1::sub(r.c[0], a.c[0], b.c[0]);
    Fp6Pack8x1::sub(r.c[1], a.c[1], b.c[1]);
}

void Fp12Pack8x1::mul(UFp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b)
{
    UFp6Pack8x1 t0{};
    UFp6Pack8x1 t1{};
    UFp6Pack8x1 t2{};
    Fp6Pack8x1 t3{};
    Fp6Pack8x1 t4{};

    Fp6Pack8x1::mul(t0, a.c[0], b.c[0]);
    Fp6Pack8x1::mul(t1, a.c[1], b.c[1]);

    // r.c[1] = (a0+a1)*(b0+b1) - a0*b0 - a1*b1
    Fp6Pack8x1::add(t3, a.c[0], a.c[1]);
    Fp6Pack8x1::add(t4, b.c[0], b.c[1]);
    Fp6Pack8x1::mul(r.c[1], t3, t4);
    UFp6Pack8x1::sub(r.c[1], r.c[1], t0);
    UFp6Pack8x1::sub(r.c[1], r.c[1], t1);

    // r.c[0] = a0*b0 + v*(a1*b1)
    UFp6Pack8x1::mul_by_nonresidue(t2, t1);
    UFp6Pack8x1::add(r.c[0], t0, t2);
}

void Fp12Pack8x1::mul(Fp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b)
{
    UFp12Pack8x1 wide{};
    Fp12Pack8x1::mul(wide, a, b);
    UFp12Pack8x1::redc(r, wide);
}

void Fp12Pack8x1::sqr(UFp12Pack8x1& r, const Fp12Pack8x1& a)
{
    UFp6Pack8x1 s0{};
    UFp6Pack8x1 s1{};
    UFp6Pack8x1 m01{};
    UFp6Pack8x1 t0{};

    Fp6Pack8x1::sqr(s0, a.c[0]);
    Fp6Pack8x1::sqr(s1, a.c[1]);
    Fp6Pack8x1::mul(m01, a.c[0], a.c[1]);
    UFp6Pack8x1::add(r.c[1], m01, m01);
    UFp6Pack8x1::mul_by_nonresidue(t0, s1);
    UFp6Pack8x1::add(r.c[0], s0, t0);
}

void Fp12Pack8x1::sqr(Fp12Pack8x1& r, const Fp12Pack8x1& a)
{
    UFp12Pack8x1 wide{};
    Fp12Pack8x1::sqr(wide, a);
    UFp12Pack8x1::redc(r, wide);
}

void Fp12Pack8x1::cyclotomic_sqr(Fp12Pack8x1& r, const Fp12Pack8x1& a)
{
    Fp2Pack8x1x1 t0_0{};
    Fp2Pack8x1x1 t0_1{};
    Fp2Pack8x1x1 t1_0{};
    Fp2Pack8x1x1 t1_1{};
    Fp2Pack8x1x1 t2_0{};
    Fp2Pack8x1x1 t2_1{};
    Fp2Pack8x1x1 tmp{};

    // t0 = sqr_fp4(a00, a11)
    sqr_fp4_pack(t0_0, t0_1, a.c[0].c[0], a.c[1].c[1]);

    // ret00 = 2*(t0_0 - a00) + t0_0
    Fp2Pack8x1x1::sub(r.c[0].c[0], t0_0, a.c[0].c[0]);
    Fp2Pack8x1x1::add(r.c[0].c[0], r.c[0].c[0], r.c[0].c[0]);
    Fp2Pack8x1x1::add(r.c[0].c[0], r.c[0].c[0], t0_0);

    // ret11 = 2*(t0_1 + a11) + t0_1
    Fp2Pack8x1x1::add(r.c[1].c[1], t0_1, a.c[1].c[1]);
    Fp2Pack8x1x1::add(r.c[1].c[1], r.c[1].c[1], r.c[1].c[1]);
    Fp2Pack8x1x1::add(r.c[1].c[1], r.c[1].c[1], t0_1);

    // Must be computed before overwriting r.c[*] when r aliases a.
    // t2 = sqr_fp4(a01, a12)
    sqr_fp4_pack(t2_0, t2_1, a.c[0].c[1], a.c[1].c[2]);

    // t1 = sqr_fp4(a10, a02)
    sqr_fp4_pack(t1_0, t1_1, a.c[1].c[0], a.c[0].c[2]);

    // ret01 = 2*(t1_0 - a01) + t1_0
    Fp2Pack8x1x1::sub(r.c[0].c[1], t1_0, a.c[0].c[1]);
    Fp2Pack8x1x1::add(r.c[0].c[1], r.c[0].c[1], r.c[0].c[1]);
    Fp2Pack8x1x1::add(r.c[0].c[1], r.c[0].c[1], t1_0);

    // ret12 = 2*(t1_1 + a12) + t1_1
    Fp2Pack8x1x1::add(r.c[1].c[2], t1_1, a.c[1].c[2]);
    Fp2Pack8x1x1::add(r.c[1].c[2], r.c[1].c[2], r.c[1].c[2]);
    Fp2Pack8x1x1::add(r.c[1].c[2], r.c[1].c[2], t1_1);

    // ret02 = 2*(t2_0 - a02) + t2_0
    Fp2Pack8x1x1::sub(r.c[0].c[2], t2_0, a.c[0].c[2]);
    Fp2Pack8x1x1::add(r.c[0].c[2], r.c[0].c[2], r.c[0].c[2]);
    Fp2Pack8x1x1::add(r.c[0].c[2], r.c[0].c[2], t2_0);

    // ret10 = 2*((u+1)*t2_1 + a10) + (u+1)*t2_1
    Fp2Pack8x1x1::mul_by_u_plus_1(tmp, t2_1);
    Fp2Pack8x1x1::add(r.c[1].c[0], tmp, a.c[1].c[0]);
    Fp2Pack8x1x1::add(r.c[1].c[0], r.c[1].c[0], r.c[1].c[0]);
    Fp2Pack8x1x1::add(r.c[1].c[0], r.c[1].c[0], tmp);
}

void Fp12Pack8x1::frobenius_map(Fp12Pack8x1& r, const Fp12Pack8x1& a, const std::size_t n)
{
    const std::size_t n_mod4 = n & 3u;
    if (n_mod4 == 0) {
        r = a;
        return;
    }

    frobenius_map_fp6_pack(r.c[0], a.c[0], n_mod4);
    frobenius_map_fp6_pack(r.c[1], a.c[1], n_mod4);

    Fp2Pack8x1x1 c0{};
    Fp2Pack8x1x1 c1{};
    Fp2Pack8x1x1 c2{};
    const Fp2Pack8x1x1& coeff = fp12_c1_coeff_packs()[n_mod4];
    Fp2Pack8x1x1::mul(c0, r.c[1].c[0], coeff);
    Fp2Pack8x1x1::mul(c1, r.c[1].c[1], coeff);
    Fp2Pack8x1x1::mul(c2, r.c[1].c[2], coeff);
    r.c[1].c[0] = c0;
    r.c[1].c[1] = c1;
    r.c[1].c[2] = c2;
}

void Fp12Pack8x1::conjugate(Fp12Pack8x1& r, const Fp12Pack8x1& a)
{
    Fp6Pack8x1 zero{};
    r.c[0] = a.c[0];
    Fp6Pack8x1::sub(r.c[1], zero, a.c[1]);
}

void Fp12Pack8x1::select(Fp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b, const u8 lane_mask)
{
    Fp6Pack8x1::select(r.c[0], a.c[0], b.c[0], lane_mask);
    Fp6Pack8x1::select(r.c[1], a.c[1], b.c[1], lane_mask);
}

void UFp12Pack8x1::redc(Fp12Pack8x1& r, const UFp12Pack8x1& a)
{
    UFp6Pack8x1::redc(r.c[0], a.c[0]);
    UFp6Pack8x1::redc(r.c[1], a.c[1]);
}

void UFp12Pack8x1::select(UFp12Pack8x1& r, const UFp12Pack8x1& a, const UFp12Pack8x1& b, const u8 lane_mask)
{
    UFp6Pack8x1::select(r.c[0], a.c[0], b.c[0], lane_mask);
    UFp6Pack8x1::select(r.c[1], a.c[1], b.c[1], lane_mask);
}

void Fp12CyclotomicState::sqr_v1(Fp12CyclotomicState& r, const Fp12CyclotomicState& in)
{
    Fp4Pack1x2x2x2& ra = r.a;
    Fp12CyclotomicBc& rbc = r.bc;
    const Fp4Pack1x2x2x2& a = in.a;
    const Fp12CyclotomicBc& bc = in.bc;
    Fp4Pack1x2x2x2 ta;
    Fp4Pack2x2x2x1 tbc, t0;
    const cc512i m0 = VSET(3, 2, 1, 0, 5, 4, 7, 6);

    Fp4Pack1x2x2x2::sqr_v1(ta, a);
    Fp2Pack2x2x2::as(ra, ta, a);
    Fp2Pack2x2x2::add(ra, ra, ra);
    Fp2Pack2x2x2::add(ra, ra, ta);

    Fp4Pack2x2x2x1::sqr_v1(tbc, bc);
    Fp2Pack4x2x1::mul_by_u_plus_1(t0, tbc);
    blend_0xC0(tbc.data(), tbc.data(), t0.data());
    perm_var(tbc.data(), tbc.data(), m0);
    Fp2Pack4x2x1::assa(rbc, tbc, bc);
    Fp2Pack4x2x1::add(rbc, rbc, rbc);
    Fp2Pack4x2x1::add(rbc, rbc, tbc);
}

void Fp12CyclotomicState::sqr_v2(Fp12CyclotomicState& r, const Fp12CyclotomicState& in)
{
    Fp4Pack1x2x2x2& ra = r.a;
    Fp12CyclotomicBc& rbc = r.bc;
    const Fp4Pack1x2x2x2& a = in.a;
    const Fp12CyclotomicBc& bc = in.bc;
    Fp4Pack1x2x2x2 ta;
    Fp4Pack2x2x2x1 tbc, t0;
    const cc512i m0 = VSET(3, 2, 1, 0, 5, 4, 7, 6);

    Fp4Pack1x2x2x2::sqr_v2(ta, a);
    Fp2Pack2x2x2::as(ra, ta, a);
    Fp2Pack2x2x2::add(ra, ra, ra);
    Fp2Pack2x2x2::add(ra, ra, ta);

    Fp4Pack2x2x2x1::sqr_v2(tbc, bc);
    Fp2Pack4x2x1::mul_by_u_plus_1(t0, tbc);
    blend_0xC0(tbc.data(), tbc.data(), t0.data());
    perm_var(tbc.data(), tbc.data(), m0);
    Fp2Pack4x2x1::assa(rbc, tbc, bc);
    Fp2Pack4x2x1::add(rbc, rbc, rbc);
    Fp2Pack4x2x1::add(rbc, rbc, tbc);
}

void Fp12CyclotomicBc::compressed_sqr_v1(Fp12CyclotomicBc& rbc, const Fp12CyclotomicBc& bc)
{
    Fp4Pack2x2x2x1 tbc, t0;
    const cc512i m0 = VSET(3, 2, 1, 0, 5, 4, 7, 6);

    Fp4Pack2x2x2x1::sqr_v1(tbc, bc);
    Fp2Pack4x2x1::mul_by_u_plus_1(t0, tbc);
    blend_0xC0(tbc.data(), tbc.data(), t0.data());
    perm_var(tbc.data(), tbc.data(), m0);
    Fp2Pack4x2x1::assa(rbc, tbc, bc);
    Fp2Pack4x2x1::add(rbc, rbc, rbc);
    Fp2Pack4x2x1::add(rbc, rbc, tbc);
}

void Fp12MulInput8x1x1::mul_v1(
    Fp2Pack4x2x1& r01,
    Fp2Pack4x2x1& r2,
    const Fp12MulInput8x1x1& ab)
{
    const Fp2Pack8x1x1& ab0 = ab[0];
    const Fp2Pack8x1x1& ab1 = ab[1];
    const Fp2Pack8x1x1& ab2 = ab[2];
    UFp2Pack8x1x1 tt01, tt2, tt3, tt4;
    UFp2Pack4x2x1 ss0, ss1, ss2, ss3;
    Fp2Pack4x2x1 s0, s1;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 2, 1, 0, 5, 4, 7, 6);

    Fp6Pack4x2x1x1::mul(tt01, tt2, ab0, ab1, ab2);
    shuf_01_UFp2Pack8x1x1(tt3, tt01);
    shuf_01_UFp2Pack8x1x1(tt4, tt2);
    blend_0x01_UFp2Pack8x1x1(tt4, tt4, tt3);
    blend_0x01_UFp2Pack8x1x1(tt2, tt2, tt01);
    blend_0x03_UFp2Pack8x1x1(tt01, tt01, tt4);
    shuf_01_UFp2Pack8x1x1(tt3, tt01);
    blend_0x55_dl(ss0.data(), tt3[1].data(), tt01[0].data());
    blend_0x55_dl(ss1.data(), tt01[1].data(), tt3[0].data());
    perm_1032_dl(ss2.data(), ss0.data());
    blend_0x33_dl(ss2.data(), ss1.data(), ss2.data());
    perm_var_dl(ss2.data(), ss2.data(), m0);
    perm_1032_dl(ss3.data(), ss1.data());
    blend_0x33_dl(ss3.data(), ss3.data(), ss0.data());
    UFp2Pack4x2x1::add(ss2, ss2, ss3);
    UFp2Pack4x2x1::redc(r01, ss2);

    shuf_01_dl(tt3[0].data(), tt2[1].data());
    blend_0x55_dl(ss0.data(), tt3[0].data(), tt2[0].data());
    UFp2Pack4x2x1::redc(s0, ss0);
    Fp2Pack4x2x1::mul_by_u_plus_1(s1, s0);
    blend_0xC0(s1.data(), s0.data(), s1.data());
    perm_var(s0.data(), s1.data(), m1);
    Fp2Pack4x2x1::add(r2, s0, s1);
}

void Fp12MulInput8x1x1::mul_v2(
    Fp2Pack4x2x1& r01,
    Fp2Pack2x2x2& r2,
    const Fp12MulInput8x1x1& ab)
{
    const Fp2Pack8x1x1& ab0 = ab[0];
    const Fp2Pack8x1x1& ab1 = ab[1];
    const Fp2Pack8x1x1& ab2 = ab[2];
    UFp2Pack8x1x1 tt01, tt2, tt3, tt4;
    UFp2Pack4x2x1 ss0, ss1, ss2, ss3;
    Fp2Pack4x2x1 s0, s1;
    UFp2Pack2x2x2 hh0, hh1, hh2;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 3, 2, 2, 1, 1, 0, 0);
    const cc512i m2 = VSET(5, 5, 4, 4, 7, 7, 6, 6);

    Fp6Pack4x2x1x1::mul(tt01, tt2, ab0, ab1, ab2);
    shuf_01_UFp2Pack8x1x1(tt3, tt01);
    shuf_01_UFp2Pack8x1x1(tt4, tt2);
    blend_0x01_UFp2Pack8x1x1(tt4, tt4, tt3);
    blend_0x01_UFp2Pack8x1x1(tt2, tt2, tt01);
    blend_0x03_UFp2Pack8x1x1(tt01, tt01, tt4);
    shuf_01_UFp2Pack8x1x1(tt3, tt01);
    blend_0x55_dl(ss0.data(), tt3[1].data(), tt01[0].data());
    blend_0x55_dl(ss1.data(), tt01[1].data(), tt3[0].data());
    perm_1032_dl(ss2.data(), ss0.data());
    blend_0x33_dl(ss2.data(), ss1.data(), ss2.data());
    perm_var_dl(ss2.data(), ss2.data(), m0);
    perm_1032_dl(ss3.data(), ss1.data());
    blend_0x33_dl(ss3.data(), ss3.data(), ss0.data());
    UFp2Pack4x2x1::add(ss2, ss2, ss3);
    UFp2Pack4x2x1::redc(r01, ss2);

    shuf_01_dl(tt3[0].data(), tt2[1].data());
    blend_0x55_dl(ss0.data(), tt3[0].data(), tt2[0].data());
    perm_var_dl(ss1.data(), ss0.data(), m1);
    conv_dltovl(hh0, ss1);
    perm_var_dl(ss1.data(), ss0.data(), m2);
    conv_dltovl(hh1, ss1);
    UFp2Pack2x2x2::mul_by_u_plus_1(hh2, hh1);
    blend_0x0F_vl(hh1.data(), hh1.data(), hh2.data());
    UFp2Pack2x2x2::add(hh0, hh0, hh1);
    UFp2Pack2x2x2::redc(r2, hh0);
}

void Fp12MulInput4x2x1::mul_v3(
    Fp2Pack2x2x2& r001,
    Fp2Pack2x2x2& r101,
    Fp2Pack2x2x2& r2,
    const Fp12MulInput4x2x1& ab)
{
    const Fp2Pack4x2x1& ab0 = ab[0];
    const Fp2Pack4x2x1& ab1 = ab[1];
    const Fp2Pack4x2x1& ab2 = ab[2];
    UFp2Pack4x2x1 tt0, tt1, tt2, tt3, tt4, tt5;
    Fp2Pack4x2x1 t0;
    UFp2Pack2x2x2 ss0, ss1, ss2, ss3, ss4, ss5;
    Fp2Pack2x2x2 s0, s1, s2, s3;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 3, 2, 2, 1, 1, 0, 0);
    const cc512i m2 = VSET(7, 7, 6, 6, 7, 7, 6, 6);
    const cc512i m3 = VSET(7, 7, 6, 6, 5, 5, 4, 4);

    Fp6Pack2x2x2x1::mul(tt0, tt1, ab0, ab1, ab2);

    perm_1100_dl(tt2.data(), tt0.data());
    perm_1100_dl(tt1.data(), tt1.data());
    perm_var_dl(tt3.data(), tt1.data(), m0);
    blend_0x0F_dl(tt4.data(), tt2.data(), tt3.data());
    perm_var_dl(tt5.data(), tt0.data(), m1);
    conv_dltovl(ss0, tt4);
    UFp2Pack2x2x2::mul_by_u_plus_1(ss1, ss0);
    blend_0x0F_vl(ss1.data(), ss0.data(), ss1.data());
    conv_dltovl(ss4, tt5);
    UFp2Pack2x2x2::add(ss1, ss4, ss1);
    UFp2Pack2x2x2::redc(r001, ss1);
    conv_dltovl(ss2, tt1);
    perm_var_dl(tt5.data(), tt0.data(), m2);
    conv_dltovl(ss1, tt5);

    perm_var(t0.data(), ab0.data(), m1);
    conv_sltohl(s0, t0);
    perm_var(t0.data(), ab0.data(), m3);
    conv_sltohl(s1, t0);
    Fp2Pack2x2x2::add(s0, s0, s1);
    perm_var(t0.data(), ab1.data(), m1);
    conv_sltohl(s1, t0);
    perm_var(t0.data(), ab1.data(), m3);
    conv_sltohl(s2, t0);
    Fp2Pack2x2x2::add(s1, s1, s2);
    perm_var(t0.data(), ab2.data(), m1);
    conv_sltohl(s2, t0);
    perm_var(t0.data(), ab2.data(), m3);
    conv_sltohl(s3, t0);
    Fp2Pack2x2x2::add(s2, s2, s3);
    Fp6Pack1x2x2x2::mul(ss3, ss5, s0, s1, s2);
    UFp2Pack2x2x2::sub(ss3, ss3, ss4);
    perm_var_dl(tt5.data(), tt0.data(), m3);
    conv_dltovl(ss4, tt5);
    UFp2Pack2x2x2::sub(ss3, ss3, ss4);
    UFp2Pack2x2x2::redc(r101, ss3);
    UFp2Pack2x2x2::sub(ss5, ss5, ss2);
    UFp2Pack2x2x2::add(ss1, ss2, ss1);
    UFp2Pack2x2x2::sub(ss5, ss5, ss0);
    perm_var_vl(ss5.data(), ss5.data(), m0);
    blend_0x0F_vl(ss1.data(), ss5.data(), ss1.data());
    UFp2Pack2x2x2::redc(r2, ss1);
}

void Fp12MulInput4x2x1::mul_v4(
    Fp2Pack4x2x1& r0,
    Fp2Pack2x2x2& r101,
    Fp2Pack2x2x2& r12,
    const Fp12MulInput4x2x1& ab)
{
    const Fp2Pack4x2x1& ab0 = ab[0];
    const Fp2Pack4x2x1& ab1 = ab[1];
    const Fp2Pack4x2x1& ab2 = ab[2];
    UFp2Pack4x2x1 tt0, tt1, tt2, tt3, tt4, tt5;
    Fp2Pack4x2x1 t0;
    UFp2Pack2x2x2 ss0, ss1, ss2, ss3, ss4;
    Fp2Pack2x2x2 s0, s1, s2, s3;
    const cc512i m0 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 3, 2, 2, 1, 1, 0, 0);
    const cc512i m2 = VSET(1, 0, 7, 6, 5, 4, 3, 2);
    const cc512i m3 = VSET(7, 7, 6, 6, 5, 5, 4, 4);

    Fp6Pack2x2x2x1::mul(tt0, tt1, ab0, ab1, ab2);

    perm_var_dl(tt2.data(), tt1.data(), m0);
    blend_0x0F_dl(tt3.data(), tt2.data(), tt0.data());
    UFp2Pack4x2x1::mul_by_u_plus_1(tt4, tt2);
    perm_var_dl(tt2.data(), tt0.data(), m2);
    blend_0x03_dl(tt2.data(), tt2.data(), tt4.data());
    UFp2Pack4x2x1::add(tt3, tt3, tt2);
    UFp2Pack4x2x1::redc(r0, tt3);

    perm_1100_dl(tt2.data(), tt0.data());
    perm_1100_dl(tt1.data(), tt1.data());
    perm_var_dl(tt3.data(), tt1.data(), m0);
    blend_0x0F_dl(tt4.data(), tt2.data(), tt3.data());
    perm_var_dl(tt5.data(), tt0.data(), m1);
    conv_dltovl(ss0, tt4);
    conv_dltovl(ss4, tt5);
    conv_dltovl(ss2, tt1);

    perm_var(t0.data(), ab0.data(), m1);
    conv_sltohl(s0, t0);
    perm_var(t0.data(), ab0.data(), m3);
    conv_sltohl(s1, t0);
    Fp2Pack2x2x2::add(s0, s0, s1);
    perm_var(t0.data(), ab1.data(), m1);
    conv_sltohl(s1, t0);
    perm_var(t0.data(), ab1.data(), m3);
    conv_sltohl(s2, t0);
    Fp2Pack2x2x2::add(s1, s1, s2);
    perm_var(t0.data(), ab2.data(), m1);
    conv_sltohl(s2, t0);
    perm_var(t0.data(), ab2.data(), m3);
    conv_sltohl(s3, t0);
    Fp2Pack2x2x2::add(s2, s2, s3);
    Fp6Pack1x2x2x2::mul(ss3, ss1, s0, s1, s2);
    UFp2Pack2x2x2::sub(ss3, ss3, ss4);
    perm_var_dl(tt5.data(), tt0.data(), m3);
    conv_dltovl(ss4, tt5);
    UFp2Pack2x2x2::sub(ss3, ss3, ss4);
    UFp2Pack2x2x2::redc(r101, ss3);
    UFp2Pack2x2x2::sub(ss1, ss1, ss2);
    UFp2Pack2x2x2::sub(ss1, ss1, ss0);
    UFp2Pack2x2x2::redc(r12, ss1);
}

void Fp12SparseInput4x2x1::mul_v1(
    Fp2Pack4x2x1& r0,
    Fp2Pack2x2x2& r1,
    const Fp12SparseInput4x2x1& in)
{
    const Fp2Pack4x2x1& a01 = in.a01;
    const Fp2Pack4x2x1& a2 = in.a2;
    const Fp2Pack4x2x1& b01 = in.b01;
    const Fp2Pack4x2x1& b4 = in.b4;
    UFp2Pack4x2x1 tt0, tt1, tt2, tt3, tt4, tt5;
    UFp2Pack2x2x2 ss0, ss1, ss2;
    const cc512i m0 = VSET(1, 0, 3, 2, 7, 6, 5, 4);
    const cc512i m1 = VSET(3, 3, 2, 2, 3, 3, 2, 2);
    const cc512i m2 = VSET(7, 7, 6, 6, 7, 7, 6, 6);

    Fp6Pack2x2x2x1::mul_by_xy00z0(tt0, tt1, tt2, tt3, a01, a2, b01, b4);

    blend_0x03_dl(tt5.data(), tt0.data(), tt1.data());
    blend_0x33_dl(tt4.data(), tt2.data(), tt3.data());
    perm_var_dl(tt4.data(), tt4.data(), m0);
    UFp2Pack4x2x1::add(tt5, tt5, tt4);
    UFp2Pack4x2x1::redc(r0, tt5);

    perm_1100_dl(tt0.data(), tt0.data());
    perm_var_dl(tt4.data(), tt3.data(), m1);
    blend_0x0F_dl(tt4.data(), tt4.data(), tt0.data());
    conv_dltovl(ss0, tt4);
    perm_1100_dl(tt1.data(), tt1.data());
    perm_var_dl(tt3.data(), tt3.data(), m2);
    blend_0x0F_dl(tt1.data(), tt1.data(), tt3.data());
    conv_dltovl(ss1, tt1);
    UFp2Pack2x2x2::mul_by_u_plus_1(ss2, ss1);
    blend_0x0F_vl(ss1.data(), ss1.data(), ss2.data());
    UFp2Pack2x2x2::add(ss1, ss1, ss0);
    UFp2Pack2x2x2::redc(r1, ss1);
}

void Fp12SparseInput8x1x1::mul_v2(
    Fp2Pack4x2x1& r0,
    Fp2Pack2x2x2& r1,
    const Fp12SparseInput8x1x1& in)
{
    const Fp2Pack8x1x1& a = in.a;
    const Fp2Pack8x1x1& b = in.b;
    UFp2Pack8x1x1 tt0, tt1, tt2, tt3;
    UFp2Pack4x2x1 ss0, ss1;
    UFp2Pack2x2x2 hh0, hh1, hh2;
    const cc512i m0 = VSET(6, 6, 7, 7, 2, 2, 1, 1);
    const cc512i m1 = VSET(0, 0, 2, 2, 6, 6, 4, 4);
    const cc512i m2 = VSET(1, 1, 1, 1, 3, 3, 3, 3);
    const cc512i m3 = VSET(5, 5, 5, 5, 7, 7, 7, 7);

    Fp6Pack4x2x1x1::mul_by_xy00z0(tt0, tt1, tt2, a, b);

    blend_0x33_UFp2Pack8x1x1(tt3, tt0, tt1);
    perm_var_UFp2Pack8x1x1(tt3, tt3, m0);
    blend_0x55_dl(ss0.data(), tt3[1].data(), tt3[0].data());
    blend_0x33_UFp2Pack8x1x1(tt3, tt1, tt2);
    perm_var_UFp2Pack8x1x1(tt3, tt3, m1);
    blend_0x55_dl(ss1.data(), tt3[1].data(), tt3[0].data());
    UFp2Pack4x2x1::add(ss0, ss0, ss1);
    UFp2Pack4x2x1::redc(r0, ss0);

    blend_0x33_UFp2Pack8x1x1(tt3, tt0, tt2);
    perm_var_UFp2Pack8x1x1(tt3, tt3, m2);
    blend_0x33_dl(ss0.data(), tt3[1].data(), tt3[0].data());
    conv_dltovl(hh0, ss0);
    perm_1032_UFp2Pack8x1x1(tt3, tt2);
    blend_0x33_UFp2Pack8x1x1(tt3, tt3, tt1);
    perm_var_UFp2Pack8x1x1(tt3, tt3, m3);
    blend_0x33_dl(ss1.data(), tt3[1].data(), tt3[0].data());
    conv_dltovl(hh1, ss1);
    UFp2Pack2x2x2::mul_by_u_plus_1(hh2, hh1);
    blend_0x0F_vl(hh1.data(), hh1.data(), hh2.data());
    UFp2Pack2x2x2::add(hh1, hh1, hh0);
    UFp2Pack2x2x2::redc(r1, hh1);
}

void Fp12SquareInput4x2x1::sqr_v1(
    Fp2Pack4x2x1& r0,
    Fp2Pack4x2x1& r1,
    const Fp12SquareInput4x2x1& in)
{
    const Fp2Pack4x2x1& a0 = in.a0;
    const Fp2Pack4x2x1& a1 = in.a1;
    Fp2Pack4x2x1 t0, t1, t2, t3, t4;
    Fp2Pack4x2x1 s0, s1, s2;
    UFp2Pack4x2x1 tt0, tt1;
    const cc512i m0 = VSET(5, 4, 3, 2, 1, 0, 7, 6);
    const cc512i m1 = VSET(3, 2, 1, 0, 7, 6, 5, 4);
    const cc512i m2 = VSET(5, 4, 1, 0, 7, 6, 3, 2);
    const cc512i m3 = VSET(7, 6, 3, 2, 5, 4, 1, 0);
    const cc512i m4 = VSET(4, 4, 1, 0, 7, 6, 3, 2);
    const cc512i m5 = VSET(5, 5, 5, 4, 7, 6, 3, 2);

    FpPack8x1::asax6(t0, a0, a1);
    blend_0x0F(t1.data(), t0.data(), a1.data());
    perm_var(t1.data(), t1.data(), m0);
    Fp2Pack4x2x1::add(t1, a0, t1);
    perm_var(t2.data(), a1.data(), m1);
    blend_0x0F(t2.data(), t2.data(), a0.data());
    perm_var(s0.data(), t2.data(), m2);
    perm_1032(t3.data(), t0.data());
    blend_0x03(t3.data(), t3.data(), t1.data());
    blend_0x0F(s0.data(), s0.data(), t3.data());
    perm_var(s1.data(), t2.data(), m3);
    perm_1032(t3.data(), t1.data());
    blend_0x03(t3.data(), t0.data(), t3.data());
    blend_0x0F(s1.data(), s1.data(), t3.data());
    perm_1032(s2.data(), a1.data());
    blend_0xC0(s2.data(), a0.data(), s2.data());
    perm_var(t2.data(), t0.data(), m3);
    perm_var(t3.data(), t1.data(), m1);
    blend_0x03(t2.data(), t2.data(), t3.data());
    blend_0x0F(s2.data(), s2.data(), t2.data());

    Fp6Pack2x2x2x1::mul(tt0, tt1, s0, s1, s2);
    UFp2Pack4x2x1::redc(t0, tt0);
    UFp2Pack4x2x1::redc(t1, tt1);

    perm_var(t2.data(), t0.data(), m1);
    blend_0x0F(t2.data(), t1.data(), t2.data());
    Fp2Pack4x2x1::add(r1, t2, t2);

    perm_var(t4.data(), t1.data(), m4);
    blend_0x0F(t4.data(), t4.data(), t0.data());
    perm_var(t3.data(), t1.data(), m5);
    blend_0x0F(t3.data(), t3.data(), t2.data());

    FpPack8x1::assx6(t3, t4, t3);
    blend_0x0F(t4.data(), t3.data(), t2.data());
    perm_var(t4.data(), t4.data(), m0);
    Fp2Pack4x2x1::sub(r0, t3, t4);
}

void mul_fp12_vector(Fp12& ret, const Fp12& a, const Fp12& b)
{
    Fp12MulInput4x2x1 ab;
    Fp2Pack2x2x2 r001, r101, r2;
    std::array<cc512i, SWORDS> in0, in1, in2;
    std::array<cc512i, SWORDS / 2> out0, out1, out2;

    for (u64 i = 0; i < SWORDS; ++i) {
        in0[i] = VSET(
            static_cast<i64>(b[1][0][1][i]), static_cast<i64>(b[1][0][0][i]),
            static_cast<i64>(a[1][0][1][i]), static_cast<i64>(a[1][0][0][i]),
            static_cast<i64>(b[0][0][1][i]), static_cast<i64>(b[0][0][0][i]),
            static_cast<i64>(a[0][0][1][i]), static_cast<i64>(a[0][0][0][i]));
        in1[i] = VSET(
            static_cast<i64>(b[1][1][1][i]), static_cast<i64>(b[1][1][0][i]),
            static_cast<i64>(a[1][1][1][i]), static_cast<i64>(a[1][1][0][i]),
            static_cast<i64>(b[0][1][1][i]), static_cast<i64>(b[0][1][0][i]),
            static_cast<i64>(a[0][1][1][i]), static_cast<i64>(a[0][1][0][i]));
        in2[i] = VSET(
            static_cast<i64>(b[1][2][1][i]), static_cast<i64>(b[1][2][0][i]),
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(b[0][2][1][i]), static_cast<i64>(b[0][2][0][i]),
            static_cast<i64>(a[0][2][1][i]), static_cast<i64>(a[0][2][0][i]));
    }

    FpPack8x1::conv_64to48(ab[0], in0.data());
    FpPack8x1::conv_64to48(ab[1], in1.data());
    FpPack8x1::conv_64to48(ab[2], in2.data());

    Fp12MulInput4x2x1::mul_v3(r001, r101, r2, ab);

    FpPack4x2::carryp(r001);
    FpPack4x2::carryp(r101);
    FpPack4x2::carryp(r2);

    FpPack4x2::conv_48to64(out0.data(), r001);
    FpPack4x2::conv_48to64(out1.data(), r101);
    FpPack4x2::conv_48to64(out2.data(), r2);

    for (u64 i = 0; i < SWORDS / 2; ++i) {
        const u64* lanes0 = reinterpret_cast<const u64*>(&out0[i]);
        const u64* lanes1 = reinterpret_cast<const u64*>(&out1[i]);
        const u64* lanes2 = reinterpret_cast<const u64*>(&out2[i]);

        ret[0][0][0][i] = lanes0[0];
        ret[0][0][0][i + (SWORDS / 2)] = lanes0[1];
        ret[0][0][1][i] = lanes0[2];
        ret[0][0][1][i + (SWORDS / 2)] = lanes0[3];
        ret[0][1][0][i] = lanes0[4];
        ret[0][1][0][i + (SWORDS / 2)] = lanes0[5];
        ret[0][1][1][i] = lanes0[6];
        ret[0][1][1][i + (SWORDS / 2)] = lanes0[7];

        ret[1][0][0][i] = lanes1[0];
        ret[1][0][0][i + (SWORDS / 2)] = lanes1[1];
        ret[1][0][1][i] = lanes1[2];
        ret[1][0][1][i + (SWORDS / 2)] = lanes1[3];
        ret[1][1][0][i] = lanes1[4];
        ret[1][1][0][i + (SWORDS / 2)] = lanes1[5];
        ret[1][1][1][i] = lanes1[6];
        ret[1][1][1][i + (SWORDS / 2)] = lanes1[7];

        ret[0][2][0][i] = lanes2[0];
        ret[0][2][0][i + (SWORDS / 2)] = lanes2[1];
        ret[0][2][1][i] = lanes2[2];
        ret[0][2][1][i + (SWORDS / 2)] = lanes2[3];
        ret[1][2][0][i] = lanes2[4];
        ret[1][2][0][i + (SWORDS / 2)] = lanes2[5];
        ret[1][2][1][i] = lanes2[6];
        ret[1][2][1][i + (SWORDS / 2)] = lanes2[7];
    }
}

void mul_by_xy00z0_fp12_vector(Fp12& ret, const Fp12& a, const Fp6& xy00z0)
{
    Fp2Pack4x2x1 r0;
    Fp12SparseInput4x2x1 in;
    Fp2Pack2x2x2 r1;
    std::array<cc512i, SWORDS> in0, in1, in2, in3;
    std::array<cc512i, SWORDS> out0;
    std::array<cc512i, SWORDS / 2> out1;

    for (u64 i = 0; i < SWORDS; ++i) {
        in0[i] = VSET(
            static_cast<i64>(a[1][1][1][i]), static_cast<i64>(a[1][1][0][i]),
            static_cast<i64>(a[1][0][1][i]), static_cast<i64>(a[1][0][0][i]),
            static_cast<i64>(a[0][1][1][i]), static_cast<i64>(a[0][1][0][i]),
            static_cast<i64>(a[0][0][1][i]), static_cast<i64>(a[0][0][0][i]));
        in1[i] = VSET(
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(a[0][2][1][i]), static_cast<i64>(a[0][2][0][i]),
            static_cast<i64>(a[0][2][1][i]), static_cast<i64>(a[0][2][0][i]));
        in2[i] = VSET(
            static_cast<i64>(xy00z0[1][1][i]), static_cast<i64>(xy00z0[1][0][i]),
            static_cast<i64>(xy00z0[0][1][i]), static_cast<i64>(xy00z0[0][0][i]),
            static_cast<i64>(xy00z0[1][1][i]), static_cast<i64>(xy00z0[1][0][i]),
            static_cast<i64>(xy00z0[0][1][i]), static_cast<i64>(xy00z0[0][0][i]));
        in3[i] = VSET(
            static_cast<i64>(xy00z0[2][1][i]), static_cast<i64>(xy00z0[2][0][i]),
            static_cast<i64>(xy00z0[2][1][i]), static_cast<i64>(xy00z0[2][0][i]),
            static_cast<i64>(xy00z0[2][1][i]), static_cast<i64>(xy00z0[2][0][i]),
            static_cast<i64>(xy00z0[2][1][i]), static_cast<i64>(xy00z0[2][0][i]));
    }

    FpPack8x1::conv_64to48(in.a01, in0.data());
    FpPack8x1::conv_64to48(in.a2, in1.data());
    FpPack8x1::conv_64to48(in.b01, in2.data());
    FpPack8x1::conv_64to48(in.b4, in3.data());

    Fp12SparseInput4x2x1::mul_v1(r0, r1, in);

    FpPack4x2::carryp(r1);

    FpPack8x1::conv_48to64(out0.data(), r0);
    FpPack4x2::conv_48to64(out1.data(), r1);

    for (u64 i = 0; i < SWORDS; ++i) {
        const u64* lanes = reinterpret_cast<const u64*>(&out0[i]);
        ret[0][2][0][i] = lanes[0];
        ret[0][2][1][i] = lanes[1];
        ret[0][1][0][i] = lanes[2];
        ret[0][1][1][i] = lanes[3];
        ret[1][0][0][i] = lanes[4];
        ret[1][0][1][i] = lanes[5];
        ret[1][1][0][i] = lanes[6];
        ret[1][1][1][i] = lanes[7];
    }

    for (u64 i = 0; i < SWORDS / 2; ++i) {
        const u64* lanes = reinterpret_cast<const u64*>(&out1[i]);
        ret[0][0][0][i] = lanes[0];
        ret[0][0][0][i + (SWORDS / 2)] = lanes[1];
        ret[0][0][1][i] = lanes[2];
        ret[0][0][1][i + (SWORDS / 2)] = lanes[3];
        ret[1][2][0][i] = lanes[4];
        ret[1][2][0][i + (SWORDS / 2)] = lanes[5];
        ret[1][2][1][i] = lanes[6];
        ret[1][2][1][i + (SWORDS / 2)] = lanes[7];
    }
}

void sqr_fp12_vector(Fp12& ret, const Fp12& a)
{
    Fp2Pack4x2x1 r0, r1;
    Fp12SquareInput4x2x1 in;
    std::array<cc512i, SWORDS> in0, in1, out0, out1;

    for (u64 i = 0; i < SWORDS; ++i) {
        in0[i] = VSET(
            static_cast<i64>(a[1][2][0][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(a[0][2][1][i]), static_cast<i64>(a[0][2][0][i]),
            static_cast<i64>(a[0][1][1][i]), static_cast<i64>(a[0][1][0][i]),
            static_cast<i64>(a[0][0][1][i]), static_cast<i64>(a[0][0][0][i]));
        in1[i] = VSET(
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][1][i]),
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(a[1][1][1][i]), static_cast<i64>(a[1][1][0][i]),
            static_cast<i64>(a[1][0][1][i]), static_cast<i64>(a[1][0][0][i]));
    }

    FpPack8x1::conv_64to48(in.a0, in0.data());
    FpPack8x1::conv_64to48(in.a1, in1.data());

    Fp12SquareInput4x2x1::sqr_v1(r0, r1, in);

    FpPack8x1::conv_48to64(out0.data(), r0);
    FpPack8x1::conv_48to64(out1.data(), r1);

    for (u64 i = 0; i < SWORDS; ++i) {
        const u64* lanes0 = reinterpret_cast<const u64*>(&out0[i]);
        const u64* lanes1 = reinterpret_cast<const u64*>(&out1[i]);
        ret[0][0][0][i] = lanes0[0];
        ret[0][0][1][i] = lanes0[1];
        ret[0][1][0][i] = lanes0[2];
        ret[0][1][1][i] = lanes0[3];
        ret[0][2][0][i] = lanes0[4];
        ret[0][2][1][i] = lanes0[5];
        ret[1][0][0][i] = lanes1[0];
        ret[1][0][1][i] = lanes1[1];
        ret[1][1][0][i] = lanes1[2];
        ret[1][1][1][i] = lanes1[3];
        ret[1][2][0][i] = lanes1[4];
        ret[1][2][1][i] = lanes1[5];
    }
}

void cyclotomic_sqr_fp12_vector(Fp12& ret, const Fp12& a)
{
    Fp12CyclotomicState in, out;
    std::array<cc512i, SWORDS / 2> lanes_a;
    std::array<cc512i, SWORDS> lanes_bc;

    for (u64 i = 0; i < SWORDS / 2; ++i) {
        lanes_a[i] = VSET(
            static_cast<i64>(a[1][1][1][i + (SWORDS / 2)]), static_cast<i64>(a[1][1][1][i]),
            static_cast<i64>(a[1][1][0][i + (SWORDS / 2)]), static_cast<i64>(a[1][1][0][i]),
            static_cast<i64>(a[0][0][1][i + (SWORDS / 2)]), static_cast<i64>(a[0][0][1][i]),
            static_cast<i64>(a[0][0][0][i + (SWORDS / 2)]), static_cast<i64>(a[0][0][0][i]));
    }
    FpPack4x2::conv_64to48(in.a, lanes_a.data());

    for (u64 i = 0; i < SWORDS; ++i) {
        lanes_bc[i] = VSET(
            static_cast<i64>(a[1][2][1][i]), static_cast<i64>(a[1][2][0][i]),
            static_cast<i64>(a[0][1][1][i]), static_cast<i64>(a[0][1][0][i]),
            static_cast<i64>(a[0][2][1][i]), static_cast<i64>(a[0][2][0][i]),
            static_cast<i64>(a[1][0][1][i]), static_cast<i64>(a[1][0][0][i]));
    }
    FpPack8x1::conv_64to48(in.bc, lanes_bc.data());

    Fp12CyclotomicState::sqr_v1(out, in);

    FpPack4x2::carryp(out.a);
    FpPack4x2::conv_48to64(lanes_a.data(), out.a);
    for (u64 i = 0; i < SWORDS / 2; ++i) {
        const u64* lanes = reinterpret_cast<const u64*>(&lanes_a[i]);
        ret[0][0][0][i] = lanes[0];
        ret[0][0][0][i + (SWORDS / 2)] = lanes[1];
        ret[0][0][1][i] = lanes[2];
        ret[0][0][1][i + (SWORDS / 2)] = lanes[3];
        ret[1][1][0][i] = lanes[4];
        ret[1][1][0][i + (SWORDS / 2)] = lanes[5];
        ret[1][1][1][i] = lanes[6];
        ret[1][1][1][i + (SWORDS / 2)] = lanes[7];
    }

    FpPack8x1::conv_48to64(lanes_bc.data(), out.bc);
    for (u64 i = 0; i < SWORDS; ++i) {
        const u64* lanes = reinterpret_cast<const u64*>(&lanes_bc[i]);
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

} // namespace blst

#endif
