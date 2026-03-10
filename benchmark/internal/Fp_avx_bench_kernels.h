#pragma once

#include <array>

#include <benchmark/benchmark.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fp2_avx.h"
#include "blst/fields/Fp6_avx.h"
#include "blst/fields/Fp_avx.h"

namespace blst::bench_internal
{
template <u64 Scale>
BLST_FORCEINLINE auto scaled_p48() -> std::array<u64, NWORDS>
{
    std::array<u64, NWORDS> out{};
    u64 carry = 0;
    for (std::size_t i = 0; i < NWORDS; ++i) {
        const u64 accum = P48[i] * Scale + carry;
        out[i] = accum & BMASK;
        carry = accum >> BRADIX;
    }
    return out;
}

inline void random_redc_fp(Fp& v, std::mt19937_64& gen)
{
    for (u64 i = 0; i < v.num_limb(); ++i) {
        v[i] = gen();
    }
    while (Fp::is_redc(v) == 0) {
        Fp::rshift(v, v, 1);
    }
}

inline void random_redc_fp2(Fp2& v, std::mt19937_64& gen)
{
    random_redc_fp(v[0], gen);
    random_redc_fp(v[1], gen);
}

inline void random_redc_fp6(Fp6& v, std::mt19937_64& gen)
{
    random_redc_fp2(v.data[0], gen);
    random_redc_fp2(v.data[1], gen);
    random_redc_fp2(v.data[2], gen);
}

struct Fp8x1WideBenchData
{
    explicit Fp8x1WideBenchData(const std::size_t n_packs)
        : n_packs_(n_packs)
        , wide_a(n_packs_)
        , wide_b(n_packs_)
        , wide_out(n_packs_)
        , lazy_out(n_packs_)
    {
        std::mt19937_64 gen(0x5546708F1A2D3E6CULL + static_cast<u64>(n_packs_));

        for (std::size_t i = 0; i < n_packs_; ++i) {
            std::array<Fp, 8> lane_a{};
            std::array<Fp, 8> lane_b{};
            std::array<Fp, 8> lane_c{};
            std::array<Fp, 8> lane_d{};
            FpPack8x1 pack_a{}, pack_b{}, pack_c{}, pack_d{};

            for (std::size_t lane = 0; lane < 8; ++lane) {
                random_redc_fp(lane_a[lane], gen);
                random_redc_fp(lane_b[lane], gen);
                random_redc_fp(lane_c[lane], gen);
                random_redc_fp(lane_d[lane], gen);
            }

            FpPack8x1::from_aos(pack_a, lane_a.data());
            FpPack8x1::from_aos(pack_b, lane_b.data());
            FpPack8x1::from_aos(pack_c, lane_c.data());
            FpPack8x1::from_aos(pack_d, lane_d.data());

            FpPack8x1::mul(wide_a[i], pack_a, pack_b);
            FpPack8x1::mul(wide_b[i], pack_c, pack_d);
        }
    }

    const std::size_t n_packs_;
    std::vector<UFpPack8x1> wide_a;
    std::vector<UFpPack8x1> wide_b;
    std::vector<UFpPack8x1> wide_out;
    std::vector<LUFpPack8x1> lazy_out;
};

struct Fp6BenchData
{
    explicit Fp6BenchData(const std::size_t packs)
        : n_packs(packs)
        , n_tuples(packs * 8)
        , a(n_tuples)
        , b(n_tuples)
        , out_wide(n_packs)
    {
        std::mt19937_64 gen(0xF6001234ULL + static_cast<u64>(packs));
        for (std::size_t i = 0; i < n_tuples; ++i) {
            random_redc_fp6(a[i], gen);
            random_redc_fp6(b[i], gen);
        }
        a_pack.resize(n_packs);
        b_pack.resize(n_packs);
        for (std::size_t p = 0; p < n_packs; ++p) {
            Fp6Pack8x1::from_aos(a_pack[p], &a[p * 8]);
            Fp6Pack8x1::from_aos(b_pack[p], &b[p * 8]);
        }
    }

    const std::size_t n_packs;
    const std::size_t n_tuples;
    std::vector<Fp6> a;
    std::vector<Fp6> b;
    std::vector<Fp6Pack8x1> a_pack;
    std::vector<Fp6Pack8x1> b_pack;
    std::vector<UFp6Pack8x1> out_wide;
};

inline void consume_wide(const std::vector<UFpPack8x1>& data)
{
    if (data.empty()) {
        return;
    }
    benchmark::DoNotOptimize(data.front()[0]);
    benchmark::DoNotOptimize(data.back()[(2 * NWORDS) - 1]);
    benchmark::ClobberMemory();
}

inline void consume_lazy(const std::vector<LUFpPack8x1>& data)
{
    if (data.empty()) {
        return;
    }
    benchmark::DoNotOptimize(data.front()[0]);
    benchmark::DoNotOptimize(data.back()[(2 * NWORDS) - 1]);
    benchmark::ClobberMemory();
}

inline void consume_wide_fp6(const std::vector<UFp6Pack8x1>& data)
{
    if (data.empty()) {
        return;
    }
    benchmark::DoNotOptimize(data.front().c[0][0][0]);
    benchmark::DoNotOptimize(data.back().c[2][1][NWORDS - 1]);
    benchmark::ClobberMemory();
}

struct BenchLUFp2Pack8x1x1
{
    std::array<LUFpPack8x1, 2> c;

    BLST_FORCEINLINE LUFpPack8x1& operator[](std::size_t i) { return c[i]; }
    BLST_FORCEINLINE const LUFpPack8x1& operator[](std::size_t i) const { return c[i]; }
};

BLST_FORCEINLINE void lufp8x1_bias_p(LUFpPack8x1& r, const LUFpPack8x1& a)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        r[i] = a[i];
    }
    for (std::size_t i = 0; i < NWORDS; ++i) {
        r[NWORDS + i] = VADD(a[NWORDS + i], VSET1(static_cast<i64>(P48[i])));
    }
}

BLST_FORCEINLINE void lufp8x1_bias_neg_p(UFpPack8x1& r, const LUFpPack8x1& a)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        r[i] = a[i];
    }
    for (std::size_t i = 0; i < NWORDS; ++i) {
        r[NWORDS + i] = VSUB(a[NWORDS + i], VSET1(static_cast<i64>(P48[i])));
    }
}

BLST_FORCEINLINE void lufp8x1_renorm_exact_add(UFpPack8x1& r, const LUFpPack8x1& a)
{
    const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
    const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
    const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
    const cc512i bmask = VSET1(BMASK);

    cc512i r0 = a[0], r1 = a[1], r2 = a[2], r3 = a[3];
    cc512i r4 = a[4], r5 = a[5], r6 = a[6], r7 = a[7];
    cc512i r8 = VSUB(a[8], p0), r9 = VSUB(a[9], p1);
    cc512i r10 = VSUB(a[10], p2), r11 = VSUB(a[11], p3);
    cc512i r12 = VSUB(a[12], p4), r13 = VSUB(a[13], p5);
    cc512i r14 = VSUB(a[14], p6), r15 = VSUB(a[15], p7);

    r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
    r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
    r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
    r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
    r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
    r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
    r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
    r8 = VADD(r8, VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
    cc512i t0 = VADD(r9, VSRA(r8, BRADIX));
    t0 = VADD(r10, VSRA(t0, BRADIX));
    t0 = VADD(r11, VSRA(t0, BRADIX));
    t0 = VADD(r12, VSRA(t0, BRADIX));
    t0 = VADD(r13, VSRA(t0, BRADIX));
    t0 = VADD(r14, VSRA(t0, BRADIX));
    t0 = VADD(r15, VSRA(t0, BRADIX));

    const cc512i smask = VSRA(t0, 63);
    r8 = VADD(r8, VAND(p0, smask)); r9 = VADD(r9, VAND(p1, smask));
    r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
    r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
    r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

    r9 = VADD(r9, VSRA(r8, BRADIX)); r8 = VAND(r8, bmask);
    r10 = VADD(r10, VSRA(r9, BRADIX)); r9 = VAND(r9, bmask);
    r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
    r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
    r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
    r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
    r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
    r15 = VAND(r15, bmask);

    r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
    r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
    r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
    r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

BLST_FORCEINLINE void lufp8x1_finish_exact_add_from_biased(UFpPack8x1& r)
{
    const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
    const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
    const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
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

    const cc512i smask = VSRA(t0, 63);
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

BLST_FORCEINLINE auto lufp8x1_carry_from_biased(UFpPack8x1& r) -> cc512i
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

BLST_FORCEINLINE void lufp8x1_repair_one_step_from_sign(UFpPack8x1& r, const cc512i smask)
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
BLST_FORCEINLINE auto lufp8x1_repair_sign_step(UFpPack8x1& r, const cc512i smask) -> cc512i
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

BLST_FORCEINLINE void lufp2_add(BenchLUFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    LUFpPack8x1::add(r[0], a[0], b[0]);
    LUFpPack8x1::add(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void lufp2_add(BenchLUFp2Pack8x1x1& r, const BenchLUFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    LUFpPack8x1::add(r[0], a[0], b[0]);
    LUFpPack8x1::add(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void lufp2_sub(BenchLUFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    LUFpPack8x1::sub(r[0], a[0], b[0]);
    LUFpPack8x1::sub(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void lufp2_sub(BenchLUFp2Pack8x1x1& r, const BenchLUFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    LUFpPack8x1::sub(r[0], a[0], b[0]);
    LUFpPack8x1::sub(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void lufp2_sub(BenchLUFp2Pack8x1x1& r, const BenchLUFp2Pack8x1x1& a, const BenchLUFp2Pack8x1x1& b)
{
    LUFpPack8x1::sub(r[0], a[0], b[0]);
    LUFpPack8x1::sub(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void lufp2_bias_p(BenchLUFp2Pack8x1x1& r, const BenchLUFp2Pack8x1x1& a)
{
    lufp8x1_bias_p(r[0], a[0]);
    lufp8x1_bias_p(r[1], a[1]);
}

BLST_FORCEINLINE void lufp2_normalize(UFp2Pack8x1x1& r, const BenchLUFp2Pack8x1x1& a)
{
    LUFpPack8x1::normalize<LUFpNormalizeBound::HiNeg1Pos2>(r[0], a[0]);
    LUFpPack8x1::normalize<LUFpNormalizeBound::HiNeg1Pos2>(r[1], a[1]);
}

BLST_FORCEINLINE void ufp8x1_add_no_repair(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        r[i] = VADD(a[i], b[i]);
    }
}

BLST_FORCEINLINE void ufp8x1_sub_no_repair(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        r[i] = VSUB(a[i], b[i]);
    }
}

BLST_FORCEINLINE void ufp2_add_no_repair(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    ufp8x1_add_no_repair(r[0], a[0], b[0]);
    ufp8x1_add_no_repair(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void ufp2_sub_no_repair(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    ufp8x1_sub_no_repair(r[0], a[0], b[0]);
    ufp8x1_sub_no_repair(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void ufp2_mul_by_u_plus_1_no_repair(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
{
    ufp8x1_sub_no_repair(r[0], a[0], a[1]);
    ufp8x1_add_no_repair(r[1], a[0], a[1]);
}

BLST_FORCEINLINE void lufp8x1_add_no_normalize(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    LUFpPack8x1 lazy{};
    LUFpPack8x1::add(lazy, a, b);
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        r[i] = lazy[i];
    }
}

BLST_FORCEINLINE void lufp8x1_sub_no_normalize(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    LUFpPack8x1 lazy{};
    LUFpPack8x1::sub(lazy, a, b);
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        r[i] = lazy[i];
    }
}

BLST_FORCEINLINE void lufp2_add_no_normalize(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    lufp8x1_add_no_normalize(r[0], a[0], b[0]);
    lufp8x1_add_no_normalize(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void lufp2_sub_no_normalize(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
{
    lufp8x1_sub_no_normalize(r[0], a[0], b[0]);
    lufp8x1_sub_no_normalize(r[1], a[1], b[1]);
}

BLST_FORCEINLINE void lufp2_mul_by_u_plus_1_no_normalize(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
{
    lufp8x1_sub_no_normalize(r[0], a[0], a[1]);
    lufp8x1_add_no_normalize(r[1], a[0], a[1]);
}

BLST_FORCEINLINE void ufp2_mul_no_repair(UFp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b)
{
    FpPack8x1 t0, t1;
    UFpPack8x1 tt0, tt1, tt2;

    FpPack8x1::mul(tt0, a[0], b[0]);
    FpPack8x1::mul(tt1, a[1], b[1]);
    FpPack8x1::add(t0, a[0], a[1]);
    FpPack8x1::add(t1, b[0], b[1]);
    FpPack8x1::mul(tt2, t0, t1);
    ufp8x1_sub_no_repair(r[0], tt0, tt1);
    ufp8x1_sub_no_repair(tt2, tt2, tt0);
    ufp8x1_sub_no_repair(r[1], tt2, tt1);
}

inline void fp6_combine_lufp_from_terms(
    UFp6Pack8x1& r,
    const UFp2Pack8x1x1& t0,
    const UFp2Pack8x1x1& t1,
    const UFp2Pack8x1x1& t2,
    const UFp2Pack8x1x1& w0,
    const UFp2Pack8x1x1& w1,
    const UFp2Pack8x1x1& w2)
{
    UFp2Pack8x1x1 strict0{}, strict1{};
    BenchLUFp2Pack8x1x1 lazy0{}, lazy1{}, lazy2{};

    lufp2_sub(lazy0, w0, t1);
    lufp2_sub(lazy1, lazy0, t2);
    lufp2_bias_p(lazy2, lazy1);
    lufp2_normalize(strict0, lazy2);
    UFp2Pack8x1x1::mul_by_u_plus_1(strict1, strict0);
    UFp2Pack8x1x1::add(r.c[0], strict1, t0);

    lufp2_sub(lazy0, w1, t0);
    lufp2_sub(lazy1, lazy0, t1);
    lufp2_bias_p(lazy2, lazy1);
    lufp2_normalize(strict0, lazy2);
    UFp2Pack8x1x1::mul_by_u_plus_1(strict1, t2);
    UFp2Pack8x1x1::add(r.c[1], strict0, strict1);

    lufp2_sub(lazy0, w2, t0);
    lufp2_sub(lazy1, lazy0, t2);
    lufp2_bias_p(lazy2, lazy1);
    lufp2_normalize(strict0, lazy2);
    UFp2Pack8x1x1::add(r.c[2], strict0, t1);
}

inline void fp6_mul_lufp_combine_only(UFp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    Fp2Pack8x1x1 aa{}, bb{};
    UFp2Pack8x1x1 t0{}, t1{}, t2{}, w0{}, w1{}, w2{};

    UFp2Pack8x1x1::mul(t0, a.c[0], b.c[0]);
    UFp2Pack8x1x1::mul(t1, a.c[1], b.c[1]);
    UFp2Pack8x1x1::mul(t2, a.c[2], b.c[2]);

    Fp2Pack8x1x1::add(aa, a.c[1], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[1], b.c[2]);
    UFp2Pack8x1x1::mul(w0, aa, bb);

    Fp2Pack8x1x1::add(aa, a.c[0], a.c[1]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[1]);
    UFp2Pack8x1x1::mul(w1, aa, bb);

    Fp2Pack8x1x1::add(aa, a.c[0], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[2]);
    UFp2Pack8x1x1::mul(w2, aa, bb);

    fp6_combine_lufp_from_terms(r, t0, t1, t2, w0, w1, w2);
}

inline void fp6_combine_no_repair_from_terms(
    UFp6Pack8x1& r,
    const UFp2Pack8x1x1& t0,
    const UFp2Pack8x1x1& t1,
    const UFp2Pack8x1x1& t2,
    const UFp2Pack8x1x1& w0,
    const UFp2Pack8x1x1& w1,
    const UFp2Pack8x1x1& w2)
{
    UFp2Pack8x1x1 t3{};

    r.c[0] = w0;
    ufp2_sub_no_repair(r.c[0], r.c[0], t1);
    ufp2_sub_no_repair(r.c[0], r.c[0], t2);
    ufp2_mul_by_u_plus_1_no_repair(t3, r.c[0]);
    ufp2_add_no_repair(r.c[0], t3, t0);

    r.c[1] = w1;
    ufp2_sub_no_repair(r.c[1], r.c[1], t0);
    ufp2_sub_no_repair(r.c[1], r.c[1], t1);
    ufp2_mul_by_u_plus_1_no_repair(t3, t2);
    ufp2_add_no_repair(r.c[1], r.c[1], t3);

    r.c[2] = w2;
    ufp2_sub_no_repair(r.c[2], r.c[2], t0);
    ufp2_sub_no_repair(r.c[2], r.c[2], t2);
    ufp2_add_no_repair(r.c[2], r.c[2], t1);
}

inline void fp6_mul_baseline_no_repair_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    Fp2Pack8x1x1 aa{}, bb{};
    UFp2Pack8x1x1 t0{}, t1{}, t2{}, w0{}, w1{}, w2{};

    ufp2_mul_no_repair(t0, a.c[0], b.c[0]);
    ufp2_mul_no_repair(t1, a.c[1], b.c[1]);
    ufp2_mul_no_repair(t2, a.c[2], b.c[2]);

    Fp2Pack8x1x1::add(aa, a.c[1], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[1], b.c[2]);
    ufp2_mul_no_repair(w0, aa, bb);

    Fp2Pack8x1x1::add(aa, a.c[0], a.c[1]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[1]);
    ufp2_mul_no_repair(w1, aa, bb);

    Fp2Pack8x1x1::add(aa, a.c[0], a.c[2]);
    Fp2Pack8x1x1::add(bb, b.c[0], b.c[2]);
    ufp2_mul_no_repair(w2, aa, bb);

    fp6_combine_no_repair_from_terms(r, t0, t1, t2, w0, w1, w2);
}
}
