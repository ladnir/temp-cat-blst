#pragma once

#include <array>
#include <cstddef>

#include "blst/fields/Fp.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fp_avx.h"
#include "blst/fields/Fp6_avx.h"

namespace blst::bench_lufp_v2
{
enum class LUFpV2RenormClass
{
    ExactAdd,
    Fp6MulCommon,
};

struct LUFpPack8x1V2 : AvxPack<2 * NWORDS>
{
    static void add(LUFpPack8x1V2& r, const LUFpPack8x1V2& a, const LUFpPack8x1V2& b)
    {
        for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
            r[limb] = VADD(a[limb], b[limb]);
        }
    }

    static void sub(LUFpPack8x1V2& r, const LUFpPack8x1V2& a, const LUFpPack8x1V2& b)
    {
        for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
            r[limb] = VSUB(a[limb], b[limb]);
        }
    }

    template <LUFpV2RenormClass C>
    static void renorm(UFpPack8x1& r, const LUFpPack8x1V2& a);
};

inline auto make_scaled_mod_plus_offset(const u64 scale, const u64 offset) -> std::array<u64, NWORDS>
{
    std::array<u64, NWORDS> out{};
    u64 carry = offset;
    for (std::size_t i = 0; i < NWORDS; ++i) {
        const u64 accum = P48[i] * scale + carry;
        out[i] = accum & BMASK;
        carry = accum >> BRADIX;
    }
    return out;
}

template <u64 Scale>
inline auto scaled_p48() -> std::array<u64, NWORDS>
{
    return make_scaled_mod_plus_offset(Scale, 0);
}

inline auto broadcast_mod_vec(const std::array<u64, NWORDS>& mod) -> std::array<cc512i, NWORDS>
{
    std::array<cc512i, NWORDS> out{};
    for (std::size_t i = 0; i < NWORDS; ++i) {
        out[i] = VSET1(static_cast<i64>(mod[i]));
    }
    return out;
}

inline const std::array<u64, NWORDS> kV2P48X1 = make_scaled_mod_plus_offset(1, 0);
inline const std::array<u64, NWORDS> kV2Bias3PPlus3 = make_scaled_mod_plus_offset(3, 3);

inline const std::array<cc512i, NWORDS> kV2P48X1Vec = broadcast_mod_vec(kV2P48X1);
inline const std::array<cc512i, NWORDS> kV2Bias3PPlus3Vec = broadcast_mod_vec(kV2Bias3PPlus3);

BLST_FORCEINLINE auto v2_mask_gt_epi64(const cc512i a, const cc512i b) -> u8
{
#if BLST_AVX_NATIVE
    return static_cast<u8>(_mm512_cmpgt_epi64_mask(to_native_512(a), to_native_512(b)));
#else
    std::array<u64, 8> lhs{};
    std::array<u64, 8> rhs{};
    _cc512_storeu_si512(lhs.data(), a);
    _cc512_storeu_si512(rhs.data(), b);
    u8 mask = 0;
    for (u8 lane = 0; lane < 8; ++lane) {
        if (static_cast<i64>(lhs[lane]) > static_cast<i64>(rhs[lane])) {
            mask |= static_cast<u8>(1u << lane);
        }
    }
    return mask;
#endif
}

BLST_FORCEINLINE void v2_copy_limbs(cc512i* dst, const LUFpPack8x1V2& src)
{
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        dst[i] = src[i];
    }
}

BLST_FORCEINLINE void v2_store_wide(UFpPack8x1& dst, const cc512i* limbs)
{
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        dst[i] = limbs[i];
    }
}

BLST_FORCEINLINE void v2_add_bias_neg_p(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        r[i] = VADD(a[i], b[i]);
    }
    for (std::size_t i = 0; i < NWORDS; ++i) {
        r[NWORDS + i] = VSUB(VADD(a[NWORDS + i], b[NWORDS + i]), VSET1(static_cast<i64>(P48[i])));
    }
}

BLST_FORCEINLINE void v2_bias_neg_p(UFpPack8x1& r, const LUFpPack8x1V2& a)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        r[i] = a[i];
    }
    for (std::size_t i = 0; i < NWORDS; ++i) {
        r[NWORDS + i] = VSUB(a[NWORDS + i], VSET1(static_cast<i64>(P48[i])));
    }
}

BLST_FORCEINLINE auto v2_neg_mask(const cc512i* limbs) -> u8;

BLST_FORCEINLINE void v2_carry_range(cc512i* limbs, const std::size_t begin, const std::size_t end)
{
    const cc512i bmask = VSET1(BMASK);
    for (std::size_t i = begin; i < end; ++i) {
        limbs[i + 1] = VADD(limbs[i + 1], VSRA(limbs[i], BRADIX));
        limbs[i] = VAND(limbs[i], bmask);
    }
}

BLST_FORCEINLINE auto v2_carry_range_with_sign(
    cc512i* limbs, const std::size_t begin, const std::size_t end) -> u8
{
    v2_carry_range(limbs, begin, end);
    return v2_neg_mask(limbs);
}

BLST_FORCEINLINE auto v2_carry_from_biased(UFpPack8x1& r) -> cc512i
{
    const cc512i bmask = VSET1(BMASK);

    r[1] = VADD(r[1], VSRA(r[0], BRADIX)); r[0] = VAND(r[0], bmask);
    r[2] = VADD(r[2], VSRA(r[1], BRADIX)); r[1] = VAND(r[1], bmask);
    r[3] = VADD(r[3], VSRA(r[2], BRADIX)); r[2] = VAND(r[2], bmask);
    r[4] = VADD(r[4], VSRA(r[3], BRADIX)); r[3] = VAND(r[3], bmask);
    r[5] = VADD(r[5], VSRA(r[4], BRADIX)); r[4] = VAND(r[4], bmask);
    r[6] = VADD(r[6], VSRA(r[5], BRADIX)); r[5] = VAND(r[5], bmask);
    r[7] = VADD(r[7], VSRA(r[6], BRADIX)); r[6] = VAND(r[6], bmask);
    r[8] = VADD(r[8], VSRA(r[7], BRADIX)); r[7] = VAND(r[7], bmask);
    cc512i t0 = VADD(r[9], VSRA(r[8], BRADIX));
    t0 = VADD(r[10], VSRA(t0, BRADIX));
    t0 = VADD(r[11], VSRA(t0, BRADIX));
    t0 = VADD(r[12], VSRA(t0, BRADIX));
    t0 = VADD(r[13], VSRA(t0, BRADIX));
    t0 = VADD(r[14], VSRA(t0, BRADIX));
    t0 = VADD(r[15], VSRA(t0, BRADIX));
    return VSRA(t0, 63);
}

BLST_FORCEINLINE void v2_repair_one_step_from_sign(UFpPack8x1& r, const cc512i smask)
{
    const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
    const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
    const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
    const cc512i bmask = VSET1(BMASK);

    r[8] = VADD(r[8], VAND(p0, smask)); r[9] = VADD(r[9], VAND(p1, smask));
    r[10] = VADD(r[10], VAND(p2, smask)); r[11] = VADD(r[11], VAND(p3, smask));
    r[12] = VADD(r[12], VAND(p4, smask)); r[13] = VADD(r[13], VAND(p5, smask));
    r[14] = VADD(r[14], VAND(p6, smask)); r[15] = VADD(r[15], VAND(p7, smask));

    r[9] = VADD(r[9], VSRA(r[8], BRADIX)); r[8] = VAND(r[8], bmask);
    r[10] = VADD(r[10], VSRA(r[9], BRADIX)); r[9] = VAND(r[9], bmask);
    r[11] = VADD(r[11], VSRA(r[10], BRADIX)); r[10] = VAND(r[10], bmask);
    r[12] = VADD(r[12], VSRA(r[11], BRADIX)); r[11] = VAND(r[11], bmask);
    r[13] = VADD(r[13], VSRA(r[12], BRADIX)); r[12] = VAND(r[12], bmask);
    r[14] = VADD(r[14], VSRA(r[13], BRADIX)); r[13] = VAND(r[13], bmask);
    r[15] = VADD(r[15], VSRA(r[14], BRADIX)); r[14] = VAND(r[14], bmask);
    r[15] = VAND(r[15], bmask);
}

template <u64 Scale>
BLST_FORCEINLINE auto v2_repair_sign_step(UFpPack8x1& r, const cc512i smask) -> cc512i
{
    static const auto kP = scaled_p48<Scale>();
    const cc512i p0 = VSET1(kP[0]), p1 = VSET1(kP[1]), p2 = VSET1(kP[2]);
    const cc512i p3 = VSET1(kP[3]), p4 = VSET1(kP[4]), p5 = VSET1(kP[5]);
    const cc512i p6 = VSET1(kP[6]), p7 = VSET1(kP[7]);
    const cc512i bmask = VSET1(BMASK);

    r[8] = VADD(r[8], VAND(p0, smask)); r[9] = VADD(r[9], VAND(p1, smask));
    r[10] = VADD(r[10], VAND(p2, smask)); r[11] = VADD(r[11], VAND(p3, smask));
    r[12] = VADD(r[12], VAND(p4, smask)); r[13] = VADD(r[13], VAND(p5, smask));
    r[14] = VADD(r[14], VAND(p6, smask)); r[15] = VADD(r[15], VAND(p7, smask));

    r[9] = VADD(r[9], VSRA(r[8], BRADIX)); r[8] = VAND(r[8], bmask);
    r[10] = VADD(r[10], VSRA(r[9], BRADIX)); r[9] = VAND(r[9], bmask);
    r[11] = VADD(r[11], VSRA(r[10], BRADIX)); r[10] = VAND(r[10], bmask);
    r[12] = VADD(r[12], VSRA(r[11], BRADIX)); r[11] = VAND(r[11], bmask);
    r[13] = VADD(r[13], VSRA(r[12], BRADIX)); r[12] = VAND(r[12], bmask);
    r[14] = VADD(r[14], VSRA(r[13], BRADIX)); r[13] = VAND(r[13], bmask);
    r[15] = VADD(r[15], VSRA(r[14], BRADIX)); r[14] = VAND(r[14], bmask);
    r[15] = VAND(r[15], bmask);

    return VSRA(r[15], 63);
}

BLST_FORCEINLINE void v2_high_sub_const(cc512i* limbs, const std::array<cc512i, NWORDS>& c)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        limbs[NWORDS + i] = VSUB(limbs[NWORDS + i], c[i]);
    }
}

BLST_FORCEINLINE void v2_high_add_const_masked(
    cc512i* limbs, const std::array<cc512i, NWORDS>& c, const u8 lane_mask)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        limbs[NWORDS + i] = VMADD(limbs[NWORDS + i], lane_mask, limbs[NWORDS + i], c[i]);
    }
}

BLST_FORCEINLINE auto v2_neg_mask(const cc512i* limbs) -> u8
{
    return v2_mask_gt_epi64(VZERO(), limbs[(2 * NWORDS) - 1]);
}

template <u64 Scale>
BLST_FORCEINLINE auto v2_sign_step(cc512i* limbs, const u8 neg_mask) -> u8
{
    const auto& kPVec = []() -> const std::array<cc512i, NWORDS>& {
        if constexpr (Scale == 1) {
            return kV2P48X1Vec;
        } else {
            static const auto kScaled = broadcast_mod_vec(make_scaled_mod_plus_offset(Scale, 0));
            return kScaled;
        }
    }();
    v2_high_add_const_masked(limbs, kPVec, neg_mask);
    return v2_carry_range_with_sign(limbs, NWORDS, (2 * NWORDS) - 1);
}

template <u64... Scales>
BLST_FORCEINLINE auto v2_sign_step_sequence(cc512i* limbs, u8 sign_mask) -> u8
{
    if constexpr (sizeof...(Scales) > 0) {
        ((sign_mask = v2_sign_step<Scales>(limbs, sign_mask)), ...);
    }
    return sign_mask;
}

template <>
inline void LUFpPack8x1V2::renorm<LUFpV2RenormClass::ExactAdd>(UFpPack8x1& r, const LUFpPack8x1V2& a)
{
    v2_bias_neg_p(r, a);
    const cc512i smask = v2_carry_from_biased(r);
    v2_repair_one_step_from_sign(r, smask);
}

BLST_FORCEINLINE void bench_copy_as_v2(LUFpPack8x1V2& dst, const UFpPack8x1& src)
{
    for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
        dst[limb] = src[limb];
    }
}

struct LUFp2Pack8x1x1V2
{
    LUFpPack8x1V2 c[2];

    BLST_FORCEINLINE LUFpPack8x1V2& operator[](std::size_t i) { return c[i]; }
    BLST_FORCEINLINE const LUFpPack8x1V2& operator[](std::size_t i) const { return c[i]; }

    static void add(LUFp2Pack8x1x1V2& r, const LUFp2Pack8x1x1V2& a, const LUFp2Pack8x1x1V2& b)
    {
        LUFpPack8x1V2::add(r[0], a[0], b[0]);
        LUFpPack8x1V2::add(r[1], a[1], b[1]);
    }

    static void sub(LUFp2Pack8x1x1V2& r, const LUFp2Pack8x1x1V2& a, const LUFp2Pack8x1x1V2& b)
    {
        LUFpPack8x1V2::sub(r[0], a[0], b[0]);
        LUFpPack8x1V2::sub(r[1], a[1], b[1]);
    }

    static void mul_by_u_plus_1(LUFp2Pack8x1x1V2& r, const LUFp2Pack8x1x1V2& a)
    {
        LUFpPack8x1V2::sub(r[0], a[0], a[1]);
        LUFpPack8x1V2::add(r[1], a[0], a[1]);
    }

    template <LUFpV2RenormClass C>
    static void renorm(UFp2Pack8x1x1& r, const LUFp2Pack8x1x1V2& a)
    {
        LUFpPack8x1V2::renorm<C>(r[0], a[0]);
        LUFpPack8x1V2::renorm<C>(r[1], a[1]);
    }
};

template <>
inline void LUFpPack8x1V2::renorm<LUFpV2RenormClass::Fp6MulCommon>(UFpPack8x1& r, const LUFpPack8x1V2& a)
{
    cc512i limbs[2 * NWORDS]{};
    v2_copy_limbs(limbs, a);
    v2_high_sub_const(limbs, kV2Bias3PPlus3Vec);
    const u8 sign_mask = v2_carry_range_with_sign(limbs, 0, (2 * NWORDS) - 1);
    v2_sign_step_sequence<4, 2, 1, 1>(limbs, sign_mask);
    v2_store_wide(r, limbs);
}

BLST_FORCEINLINE void bench_copy_as_v2(LUFp2Pack8x1x1V2& dst, const UFp2Pack8x1x1& src)
{
    bench_copy_as_v2(dst[0], src[0]);
    bench_copy_as_v2(dst[1], src[1]);
}

BLST_FORCEINLINE void mul_fp6_pack_lufpv2_combine_only_raw(
    LUFp2Pack8x1x1V2& c0,
    LUFp2Pack8x1x1V2& c1,
    LUFp2Pack8x1x1V2& c2,
    const LUFp2Pack8x1x1V2& t0,
    const LUFp2Pack8x1x1V2& t1,
    const LUFp2Pack8x1x1V2& t2,
    const LUFp2Pack8x1x1V2& w0,
    const LUFp2Pack8x1x1V2& w1,
    const LUFp2Pack8x1x1V2& w2)
{
    LUFp2Pack8x1x1V2 t3{};

    c0 = w0;
    LUFp2Pack8x1x1V2::sub(c0, c0, t1);
    LUFp2Pack8x1x1V2::sub(c0, c0, t2);
    LUFp2Pack8x1x1V2::mul_by_u_plus_1(t3, c0);
    LUFp2Pack8x1x1V2::add(c0, t3, t0);

    c1 = w1;
    LUFp2Pack8x1x1V2::sub(c1, c1, t0);
    LUFp2Pack8x1x1V2::sub(c1, c1, t1);
    t3 = t2;
    LUFp2Pack8x1x1V2::mul_by_u_plus_1(t3, t3);
    LUFp2Pack8x1x1V2::add(c1, c1, t3);

    c2 = w2;
    LUFp2Pack8x1x1V2::sub(c2, c2, t0);
    LUFp2Pack8x1x1V2::sub(c2, c2, t2);
    LUFp2Pack8x1x1V2::add(c2, c2, t1);
}
}
