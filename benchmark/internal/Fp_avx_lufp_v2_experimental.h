#pragma once

#include "benchmark/internal/Fp_avx_lufp_v2_impl.h"

namespace blst::bench_lufp_v2
{
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
