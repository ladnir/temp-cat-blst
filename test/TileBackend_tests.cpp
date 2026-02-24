#include <gtest/gtest.h>

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/fields/Fp12.h"
#include "blst/curve/point.h"
#include "blst/backend/backend.h"
#include "blst/backend/avx512/tile_ops.h"

using namespace blst;

namespace blst
{
BLST_HD void line_add(Fp6& line, G2Point& T, const G2Point& R, const G2AffinePoint& Q);
BLST_HD void line_dbl(Fp6& line, G2Point& T, const G2Point& Q);
}

TEST(TileBackend, AvailabilityMacroMatchesConfig)
{
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
    EXPECT_EQ(BLST_TILE_BACKEND_AVAILABLE, 1);
#else
    EXPECT_EQ(BLST_TILE_BACKEND_AVAILABLE, 0);
#endif
}

#if BLST_TILE_BACKEND_AVAILABLE
namespace {
constexpr u64 kAvxWidth = kTileWidth;
constexpr u64 kAvxLimbs = NLIMBS(384);
constexpr u64 kAvxSeedX = 0x41;
constexpr u64 kAvxSeedY = 0x59;
constexpr u64 kAvxLaneMulX = 0x17;
constexpr u64 kAvxLaneMulY = 0x13;
constexpr u64 kAvxLimbMulX = 0x05;
constexpr u64 kAvxLimbMulY = 0x07;
constexpr u64 kFp2SeedX0 = 0x61;
constexpr u64 kFp2SeedX1 = 0x71;
constexpr u64 kFp2SeedY0 = 0x81;
constexpr u64 kFp2SeedY1 = 0x91;
constexpr u64 kFp2LaneMulX0 = 0x11;
constexpr u64 kFp2LaneMulX1 = 0x13;
constexpr u64 kFp2LaneMulY0 = 0x17;
constexpr u64 kFp2LaneMulY1 = 0x1B;
constexpr u64 kFp2LimbMulX0 = 0x03;
constexpr u64 kFp2LimbMulX1 = 0x05;
constexpr u64 kFp2LimbMulY0 = 0x07;
constexpr u64 kFp2LimbMulY1 = 0x09;
constexpr u64 kFp12SeedX = 0xA1;
constexpr u64 kFp12SeedY = 0xB1;
constexpr u64 kFp12LaneMulX = 0x0D;
constexpr u64 kFp12LaneMulY = 0x11;
constexpr u64 kFp12IMulX = 0x09;
constexpr u64 kFp12IMulY = 0x07;
constexpr u64 kFp12JMulX = 0x05;
constexpr u64 kFp12JMulY = 0x03;
constexpr u64 kFp12KMulX = 0x03;
constexpr u64 kFp12KMulY = 0x05;
constexpr u64 kFp12LimbMulX = 0x01;
constexpr u64 kFp12LimbMulY = 0x02;
constexpr u64 kTileMemStridePad = 3;
constexpr limb_t kTileMemMask = ((((limb_t)1) << kTileLimbBits) - 1);
constexpr u64 kIntrinSeedA = 0x0123456789ABCDEFULL;
constexpr u64 kIntrinSeedB = 0x0FEDCBA987654321ULL;
constexpr u64 kIntrinLaneStepA = 0x1111111111111111ULL;
constexpr u64 kIntrinLaneStepB = 0x0101010101010101ULL;
constexpr u64 kLineSeed = 0xD1;
constexpr u64 kLineLaneMul = 0x09;
constexpr u64 kLineIMul = 0x07;
constexpr u64 kLineKMul = 0x05;
constexpr u64 kLineLimbMul = 0x03;
#if defined(BLST_AVX512_PORTABLE) && !defined(BLST_AVX512)
constexpr u64 kIfmaPortableIterations = 64;
constexpr u64 kIfmaPortableSeed = 0x9E3779B97F4A7C15ULL;
constexpr u64 kIfmaPortableWordMask = 0xFFFFFFFFFFFFULL;
constexpr u64 kFpLimbs = NLIMBS(384);
constexpr u64 kLinePortableFuzzIterations = 32;
#endif

BLST_FORCEINLINE void set_fp_tile_lane_avx(FpTileRegAvx512& reg, u64 lane, const Fp& x)
{
    fp_lane_pack(reg, lane, x);
}

BLST_FORCEINLINE auto get_fp_tile_lane_avx(const FpTileRegAvx512& reg, u64 lane) -> Fp
{
    Fp out;
    fp_lane_unpack(out, reg, lane);
    return out;
}

BLST_FORCEINLINE void set_fp2_tile_lane_avx(Fp2TileRegAvx512& reg, u64 lane, const Fp2& x)
{
    fp_lane_pack(reg[0], lane, x[0]);
    fp_lane_pack(reg[1], lane, x[1]);
}

BLST_FORCEINLINE auto get_fp2_tile_lane_avx(const Fp2TileRegAvx512& reg, u64 lane) -> Fp2
{
    Fp2 out;
    fp_lane_unpack(out[0], reg[0], lane);
    fp_lane_unpack(out[1], reg[1], lane);
    return out;
}

BLST_FORCEINLINE void set_fp6_tile_lane_avx(Fp6TileRegAvx512& reg, u64 lane, const Fp6& x)
{
    for (u64 i = 0; i < 3; ++i) {
        fp_lane_pack(reg[i][0], lane, x[i][0]);
        fp_lane_pack(reg[i][1], lane, x[i][1]);
    }
}

BLST_FORCEINLINE auto get_fp6_tile_lane_avx(const Fp6TileRegAvx512& reg, u64 lane) -> Fp6
{
    Fp6 out;
    for (u64 i = 0; i < 3; ++i) {
        fp_lane_unpack(out[i][0], reg[i][0], lane);
        fp_lane_unpack(out[i][1], reg[i][1], lane);
    }
    return out;
}

BLST_FORCEINLINE void set_fp12_tile_lane_avx(Fp12TileRegAvx512& reg, u64 lane, const Fp12& x)
{
    for (u64 i = 0; i < 2; ++i) {
        for (u64 j = 0; j < 3; ++j) {
            for (u64 k = 0; k < 2; ++k) {
                fp_lane_pack(reg[i][j][k], lane, x[i][j][k]);
            }
        }
    }
}

BLST_FORCEINLINE auto get_fp12_tile_lane_avx(const Fp12TileRegAvx512& reg, u64 lane) -> Fp12
{
    Fp12 out;
    for (u64 i = 0; i < 2; ++i) {
        for (u64 j = 0; j < 3; ++j) {
            for (u64 k = 0; k < 2; ++k) {
                fp_lane_unpack(out[i][j][k], reg[i][j][k], lane);
            }
        }
    }
    return out;
}

BLST_FORCEINLINE void fp12_add_ref(Fp12& out, const Fp12& a, const Fp12& b)
{
    for (u64 i = 0; i < 2; ++i) {
        for (u64 j = 0; j < 3; ++j) {
            Fp2::add(out[i][j], a[i][j], b[i][j]);
        }
    }
}

BLST_FORCEINLINE void fp12_sub_ref(Fp12& out, const Fp12& a, const Fp12& b)
{
    for (u64 i = 0; i < 2; ++i) {
        for (u64 j = 0; j < 3; ++j) {
            Fp2::sub(out[i][j], a[i][j], b[i][j]);
        }
    }
}

BLST_FORCEINLINE void fp6_mul_ref(Fp6& out, const Fp6& a, const Fp6& b)
{
    Fp2 t0;
    Fp2 t1;
    Fp2 t2;
    Fp2 t3;
    Fp2 t4;
    Fp2 t5;

    Fp2::mul(t0, a[0], b[0]);
    Fp2::mul(t1, a[1], b[1]);
    Fp2::mul(t2, a[2], b[2]);

    Fp2::add(t4, a[1], a[2]);
    Fp2::add(t5, b[1], b[2]);
    Fp2::mul(t3, t4, t5);
    Fp2::sub(t3, t3, t1);
    Fp2::sub(t3, t3, t2);
    Fp2::mul_by_u_plus_1(t3, t3);
    Fp2::add(out[0], t3, t0);

    Fp2::add(t4, a[0], a[1]);
    Fp2::add(t5, b[0], b[1]);
    Fp2::mul(out[1], t4, t5);
    Fp2::sub(out[1], out[1], t0);
    Fp2::sub(out[1], out[1], t1);
    Fp2::mul_by_u_plus_1(t3, t2);
    Fp2::add(out[1], out[1], t3);

    Fp2::add(t4, a[0], a[2]);
    Fp2::add(t5, b[0], b[2]);
    Fp2::mul(out[2], t4, t5);
    Fp2::sub(out[2], out[2], t0);
    Fp2::sub(out[2], out[2], t2);
    Fp2::add(out[2], out[2], t1);
}

BLST_FORCEINLINE void fp6_mul_by_0y0_ref(Fp6& out, const Fp6& a, const Fp2& y)
{
    Fp2 t;
    Fp2::mul(t, a[2], y);
    Fp2::mul(out[2], a[1], y);
    Fp2::mul(out[1], a[0], y);
    Fp2::mul_by_u_plus_1(out[0], t);
}

BLST_FORCEINLINE void fp6_mul_by_xy0_ref(Fp6& out, const Fp6& a, const Fp6& b)
{
    Fp2 t0;
    Fp2 t1;
    Fp2 t3;
    Fp2 t4;
    Fp2 t5;

    Fp2::mul(t0, a[0], b[0]);
    Fp2::mul(t1, a[1], b[1]);

    Fp2::mul(t3, a[2], b[1]);
    Fp2::mul_by_u_plus_1(t3, t3);
    Fp2::add(out[0], t3, t0);

    Fp2::add(t4, a[0], a[1]);
    Fp2::add(t5, b[0], b[1]);
    Fp2::mul(out[1], t4, t5);
    Fp2::sub(out[1], out[1], t0);
    Fp2::sub(out[1], out[1], t1);

    Fp2::mul(out[2], a[2], b[0]);
    Fp2::add(out[2], out[2], t1);
}

BLST_FORCEINLINE auto fp6_equal_ref(const Fp6& a, const Fp6& b) -> bool
{
    for (u64 i = 0; i < 3; ++i) {
        for (u64 j = 0; j < 2; ++j) {
            for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                if (a[i][j][limb] != b[i][j][limb]) {
                    return false;
                }
            }
        }
    }
    return true;
}

BLST_FORCEINLINE auto make_cc512_from_lanes(const u64 lanes[8]) -> cc512i
{
    return _cc512_loadu_si512(static_cast<const void*>(lanes));
}

BLST_FORCEINLINE void store_cc512_to_lanes(u64 lanes[8], const cc512i& v)
{
    _cc512_storeu_si512(static_cast<void*>(lanes), v);
}

BLST_FORCEINLINE auto madd52_lo_ref(u64 acc, u64 a, u64 b) -> u64
{
    constexpr u64 kMask26 = ((u64)1 << 26) - 1;
    constexpr u64 kMask52 = ((u64)1 << 52) - 1;

    const u64 ax = a & kMask52;
    const u64 bx = b & kMask52;
    const u64 a0 = ax & kMask26;
    const u64 a1 = ax >> 26;
    const u64 b0 = bx & kMask26;
    const u64 b1 = bx >> 26;

    const u64 p0 = a0 * b0;
    const u64 p1 = a0 * b1 + a1 * b0;
    const u64 coeff0 = p0 & kMask26;
    const u64 coeff1 = ((p0 >> 26) + (p1 & kMask26)) & kMask26;
    const u64 low52 = coeff0 | (coeff1 << 26);
    return acc + low52;
}

BLST_FORCEINLINE auto madd52_hi_ref(u64 acc, u64 a, u64 b) -> u64
{
    constexpr u64 kMask26 = ((u64)1 << 26) - 1;
    constexpr u64 kMask52 = ((u64)1 << 52) - 1;

    const u64 ax = a & kMask52;
    const u64 bx = b & kMask52;
    const u64 a0 = ax & kMask26;
    const u64 a1 = ax >> 26;
    const u64 b0 = bx & kMask26;
    const u64 b1 = bx >> 26;

    const u64 p0 = a0 * b0;
    const u64 p1 = a0 * b1 + a1 * b0;
    const u64 p2 = a1 * b1;

    u64 coeff1 = (p0 >> 26) + (p1 & kMask26);
    const u64 carry1 = coeff1 >> 26;
    coeff1 &= kMask26;

    u64 coeff2 = (p1 >> 26) + (p2 & kMask26) + carry1;
    const u64 carry2 = coeff2 >> 26;
    coeff2 &= kMask26;

    const u64 coeff3 = (p2 >> 26) + carry2;
    const u64 high52 = coeff2 | ((coeff3 & kMask26) << 26);
    return acc + high52;
}

BLST_FORCEINLINE auto next_test_u64(u64& state) -> u64
{
    constexpr u64 kMul = 6364136223846793005ULL;
    constexpr u64 kAdd = 1442695040888963407ULL;
    state = state * kMul + kAdd;
    return state;
}

BLST_FORCEINLINE void fill_random_redc_fp_any(Fp& out, u64& state)
{
    for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
        out[limb] = static_cast<limb_t>(next_test_u64(state));
    }
    while (Fp::is_redc(out) == 0) {
        Fp::rshift(out, out, 1);
    }
}

BLST_FORCEINLINE void fill_random_redc_fp2_any(Fp2& out, u64& state)
{
    fill_random_redc_fp_any(out[0], state);
    fill_random_redc_fp_any(out[1], state);
}

BLST_FORCEINLINE void fill_random_redc_fp6_any(Fp6& out, u64& state)
{
    fill_random_redc_fp2_any(out[0], state);
    fill_random_redc_fp2_any(out[1], state);
    fill_random_redc_fp2_any(out[2], state);
}

BLST_FORCEINLINE void fill_random_redc_fp12_any(Fp12& out, u64& state)
{
    fill_random_redc_fp6_any(out[0], state);
    fill_random_redc_fp6_any(out[1], state);
}

#if defined(BLST_AVX512_PORTABLE) && !defined(BLST_AVX512)
BLST_FORCEINLINE auto next_lcg64(u64& state) -> u64
{
    constexpr u64 kMul = 6364136223846793005ULL;
    constexpr u64 kAdd = 1442695040888963407ULL;
    state = state * kMul + kAdd;
    return state;
}

BLST_FORCEINLINE void fill_random_redc_fp(Fp& out, u64& state)
{
    for (u64 limb = 0; limb < kFpLimbs; ++limb) {
        out[limb] = static_cast<limb_t>(next_lcg64(state));
    }
    while (Fp::is_redc(out) == 0) {
        Fp::rshift(out, out, 1);
    }
}

BLST_FORCEINLINE void fill_words48_from_fp(u64 words[8], const Fp& x)
{
    fp_to_words48(words, x);
    for (u64 limb = 0; limb < 8; ++limb) {
        words[limb] &= kIfmaPortableWordMask;
    }
}

BLST_FORCEINLINE void set_reg_from_words(cc512i reg[8], const u64 words_by_limb[8][8])
{
    for (u64 limb = 0; limb < 8; ++limb) {
        for (u64 lane = 0; lane < 8; ++lane) {
            reg[limb].lane[lane] = words_by_limb[limb][lane] & kIfmaPortableWordMask;
        }
    }
}

BLST_FORCEINLINE void get_words_from_reg(u64 words_by_limb[8][8], const cc512i reg[8])
{
    for (u64 limb = 0; limb < 8; ++limb) {
        for (u64 lane = 0; lane < 8; ++lane) {
            words_by_limb[limb][lane] = reg[limb].lane[lane] & kIfmaPortableWordMask;
        }
    }
}

BLST_FORCEINLINE void fill_random_redc_fp2(Fp2& out, u64& state)
{
    fill_random_redc_fp(out[0], state);
    fill_random_redc_fp(out[1], state);
}

BLST_FORCEINLINE void fill_random_point_like(G2Point& out, u64& state)
{
    fill_random_redc_fp2(out.X, state);
    fill_random_redc_fp2(out.Y, state);
    fill_random_redc_fp2(out.Z, state);
}

BLST_FORCEINLINE void fill_random_affine_like(G2AffinePoint& out, u64& state)
{
    fill_random_redc_fp2(out.X, state);
    fill_random_redc_fp2(out.Y, state);
}
#endif
} // namespace

TEST(TileBackend, AvxFpMulSqrMontMatchesScalar)
{
    FpTileRegAvx512 a_reg;
    FpTileRegAvx512 b_reg;
    FpTileRegAvx512 mul_reg;
    FpTileRegAvx512 sqr_reg;

    Fp a_vals[kAvxWidth];
    Fp b_vals[kAvxWidth];

    // Deterministic per-lane values keep tests reproducible and easy to diff.
    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        Fp x;
        Fp y;
        for (u64 i = 0; i < kAvxLimbs; ++i) {
            x[i] = static_cast<limb_t>(kAvxSeedX + lane * kAvxLaneMulX + i * kAvxLimbMulX);
            y[i] = static_cast<limb_t>(kAvxSeedY + lane * kAvxLaneMulY + i * kAvxLimbMulY);
        }
        while (Fp::is_redc(x) == 0) {
            Fp::rshift(x, x, 1);
        }
        while (Fp::is_redc(y) == 0) {
            Fp::rshift(y, y, 1);
        }
        a_vals[lane] = x;
        b_vals[lane] = y;
        set_fp_tile_lane_avx(a_reg, lane, a_vals[lane]);
        set_fp_tile_lane_avx(b_reg, lane, b_vals[lane]);
    }

    fp_tile_mul_mont(mul_reg, a_reg, b_reg, get_BLST_P().data, BLST_p0);
    fp_tile_sqr_mont(sqr_reg, a_reg, get_BLST_P().data, BLST_p0);

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        Fp expected_mul;
        Fp expected_sqr;
        Fp::mul(expected_mul, a_vals[lane], b_vals[lane]);
        Fp::sqr(expected_sqr, a_vals[lane]);
        Fp got_mul = get_fp_tile_lane_avx(mul_reg, lane);
        Fp got_sqr = get_fp_tile_lane_avx(sqr_reg, lane);
        for (u64 i = 0; i < kAvxLimbs; ++i) {
            EXPECT_EQ(got_mul[i], expected_mul[i]);
            EXPECT_EQ(got_sqr[i], expected_sqr[i]);
        }
    }
}

TEST(TileBackend, AvxFp2MulSqrMontMatchesScalar)
{
    Fp2TileRegAvx512 a_reg;
    Fp2TileRegAvx512 b_reg;
    Fp2TileRegAvx512 mul_reg;
    Fp2TileRegAvx512 sqr_reg;

    Fp2 a_vals[kAvxWidth];
    Fp2 b_vals[kAvxWidth];

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        Fp x0;
        Fp x1;
        Fp y0;
        Fp y1;
        for (u64 i = 0; i < kAvxLimbs; ++i) {
            x0[i] = static_cast<limb_t>(kFp2SeedX0 + lane * kFp2LaneMulX0 + i * kFp2LimbMulX0);
            x1[i] = static_cast<limb_t>(kFp2SeedX1 + lane * kFp2LaneMulX1 + i * kFp2LimbMulX1);
            y0[i] = static_cast<limb_t>(kFp2SeedY0 + lane * kFp2LaneMulY0 + i * kFp2LimbMulY0);
            y1[i] = static_cast<limb_t>(kFp2SeedY1 + lane * kFp2LaneMulY1 + i * kFp2LimbMulY1);
        }
        while (Fp::is_redc(x0) == 0) { Fp::rshift(x0, x0, 1); }
        while (Fp::is_redc(x1) == 0) { Fp::rshift(x1, x1, 1); }
        while (Fp::is_redc(y0) == 0) { Fp::rshift(y0, y0, 1); }
        while (Fp::is_redc(y1) == 0) { Fp::rshift(y1, y1, 1); }
        a_vals[lane][0] = x0; a_vals[lane][1] = x1;
        b_vals[lane][0] = y0; b_vals[lane][1] = y1;
        set_fp2_tile_lane_avx(a_reg, lane, a_vals[lane]);
        set_fp2_tile_lane_avx(b_reg, lane, b_vals[lane]);
    }

    fp2_tile_mul_mont(mul_reg, a_reg, b_reg, get_BLST_P().data, BLST_p0);
    fp2_tile_sqr_mont(sqr_reg, a_reg, get_BLST_P().data, BLST_p0);

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        Fp2 expected_mul;
        Fp2 expected_sqr;
        Fp2::mul(expected_mul, a_vals[lane], b_vals[lane]);
        Fp2::sqr(expected_sqr, a_vals[lane]);
        Fp2 got_mul = get_fp2_tile_lane_avx(mul_reg, lane);
        Fp2 got_sqr = get_fp2_tile_lane_avx(sqr_reg, lane);
        for (u64 c = 0; c < 2; ++c) {
            for (u64 i = 0; i < kAvxLimbs; ++i) {
                EXPECT_EQ(got_mul[c][i], expected_mul[c][i]);
                EXPECT_EQ(got_sqr[c][i], expected_sqr[c][i]);
            }
        }
    }
}

TEST(TileBackend, AvxFp12AddSubConjugateMatchesScalar)
{
    Fp12TileRegAvx512 a_reg;
    Fp12TileRegAvx512 b_reg;
    Fp12TileRegAvx512 add_reg;
    Fp12TileRegAvx512 sub_reg;
    Fp12TileRegAvx512 conj_reg;

    Fp12 a_vals[kAvxWidth];
    Fp12 b_vals[kAvxWidth];

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        for (u64 i = 0; i < 2; ++i) {
            for (u64 j = 0; j < 3; ++j) {
                for (u64 k = 0; k < 2; ++k) {
                    Fp x;
                    Fp y;
                    for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                        x[limb] = static_cast<limb_t>(kFp12SeedX + lane * kFp12LaneMulX + i * kFp12IMulX + j * kFp12JMulX + k * kFp12KMulX + limb * kFp12LimbMulX);
                        y[limb] = static_cast<limb_t>(kFp12SeedY + lane * kFp12LaneMulY + i * kFp12IMulY + j * kFp12JMulY + k * kFp12KMulY + limb * kFp12LimbMulY);
                    }
                    while (Fp::is_redc(x) == 0) { Fp::rshift(x, x, 1); }
                    while (Fp::is_redc(y) == 0) { Fp::rshift(y, y, 1); }
                    a_vals[lane][i][j][k] = x;
                    b_vals[lane][i][j][k] = y;
                }
            }
        }
        set_fp12_tile_lane_avx(a_reg, lane, a_vals[lane]);
        set_fp12_tile_lane_avx(b_reg, lane, b_vals[lane]);
    }

    fp12_tile_add_mod(add_reg, a_reg, b_reg, get_BLST_P().data);
    fp12_tile_sub_mod(sub_reg, a_reg, b_reg, get_BLST_P().data);
    fp12_tile_conjugate(conj_reg, a_reg, get_BLST_P().data);

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        Fp12 expected_add;
        Fp12 expected_sub;
        Fp12 expected_conj = a_vals[lane];
        fp12_add_ref(expected_add, a_vals[lane], b_vals[lane]);
        fp12_sub_ref(expected_sub, a_vals[lane], b_vals[lane]);
        Fp12::conjugate(expected_conj);

        Fp12 got_add = get_fp12_tile_lane_avx(add_reg, lane);
        Fp12 got_sub = get_fp12_tile_lane_avx(sub_reg, lane);
        Fp12 got_conj = get_fp12_tile_lane_avx(conj_reg, lane);

        for (u64 i = 0; i < 2; ++i) {
            for (u64 j = 0; j < 3; ++j) {
                for (u64 k = 0; k < 2; ++k) {
                    for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                        EXPECT_EQ(got_add[i][j][k][limb], expected_add[i][j][k][limb]);
                        EXPECT_EQ(got_sub[i][j][k][limb], expected_sub[i][j][k][limb]);
                        EXPECT_EQ(got_conj[i][j][k][limb], expected_conj[i][j][k][limb]);
                    }
                }
            }
        }
    }
}

TEST(TileBackend, AvxFp12CyclotomicAndSparseMulMatchScalar)
{
    Fp12TileRegAvx512 a_reg;
    Fp12TileRegAvx512 cyc_reg;
    Fp12TileRegAvx512 sparse_mul_reg;
    Fp6TileRegAvx512 sparse_reg;

    Fp12 a_vals[kAvxWidth];
    Fp6 sparse_vals[kAvxWidth];

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        for (u64 i = 0; i < 2; ++i) {
            for (u64 j = 0; j < 3; ++j) {
                for (u64 k = 0; k < 2; ++k) {
                    Fp x;
                    for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                        x[limb] = static_cast<limb_t>(kFp12SeedX + lane * kFp12LaneMulX + i * kFp12IMulX + j * kFp12JMulX + k * kFp12KMulX + limb * kFp12LimbMulX);
                    }
                    while (Fp::is_redc(x) == 0) { Fp::rshift(x, x, 1); }
                    a_vals[lane][i][j][k] = x;
                }
            }
        }

        for (u64 j = 0; j < 3; ++j) {
            for (u64 k = 0; k < 2; ++k) {
                Fp s;
                for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                    s[limb] = static_cast<limb_t>(0xC1 + lane * 0x0F + j * 0x07 + k * 0x05 + limb * 0x03);
                }
                while (Fp::is_redc(s) == 0) { Fp::rshift(s, s, 1); }
                sparse_vals[lane][j][k] = s;
            }
        }

        set_fp12_tile_lane_avx(a_reg, lane, a_vals[lane]);
        set_fp6_tile_lane_avx(sparse_reg, lane, sparse_vals[lane]);
    }

    fp12_tile_cyclotomic_sqr(cyc_reg, a_reg);
    fp12_tile_mul_by_xy00z0(sparse_mul_reg, a_reg, sparse_reg);

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        Fp12 expected_cyc;
        Fp12 expected_sparse;
        Fp12::cyclotomic_sqr(expected_cyc, a_vals[lane]);
        Fp12::mul_by_xy00z0(expected_sparse, a_vals[lane], sparse_vals[lane]);

        const Fp12 got_cyc = get_fp12_tile_lane_avx(cyc_reg, lane);
        const Fp12 got_sparse = get_fp12_tile_lane_avx(sparse_mul_reg, lane);

        for (u64 i = 0; i < 2; ++i) {
            for (u64 j = 0; j < 3; ++j) {
                for (u64 k = 0; k < 2; ++k) {
                    for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                        EXPECT_EQ(got_cyc[i][j][k][limb], expected_cyc[i][j][k][limb]);
                        EXPECT_EQ(got_sparse[i][j][k][limb], expected_sparse[i][j][k][limb]);
                    }
                }
            }
        }
    }
}

TEST(TileBackend, AvxFp6MulAndSparseHelpersMatchScalar)
{
    constexpr u64 kIterations = 32;
    u64 state = 0x6A09E667F3BCC909ULL;

    Fp6TileRegAvx512 a_reg;
    Fp6TileRegAvx512 b_reg;
    Fp6TileRegAvx512 mul_reg;
    Fp6TileRegAvx512 xy0_reg;
    Fp6TileRegAvx512 y0_reg;

    for (u64 iter = 0; iter < kIterations; ++iter) {
        Fp6 a_vals[kAvxWidth];
        Fp6 b_vals[kAvxWidth];
        Fp6 b_xy0_vals[kAvxWidth];

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            fill_random_redc_fp6_any(a_vals[lane], state);
            fill_random_redc_fp6_any(b_vals[lane], state);
            fill_random_redc_fp6_any(b_xy0_vals[lane], state);
            b_xy0_vals[lane][2] = Fp2::zero();

            set_fp6_tile_lane_avx(a_reg, lane, a_vals[lane]);
            set_fp6_tile_lane_avx(b_reg, lane, b_vals[lane]);
        }

        fp6_tile_mul_mont(mul_reg, a_reg, b_reg, get_BLST_P().data, BLST_p0);

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            set_fp6_tile_lane_avx(b_reg, lane, b_xy0_vals[lane]);
        }
        fp6_tile_mul_by_xy0(xy0_reg, a_reg, b_reg, get_BLST_P().data, BLST_p0);

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            set_fp6_tile_lane_avx(b_reg, lane, b_vals[lane]);
        }
        fp6_tile_mul_by_0y0(y0_reg, a_reg, b_reg[2], get_BLST_P().data, BLST_p0);

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            Fp6 expected_mul;
            Fp6 expected_xy0;
            Fp6 expected_0y0;

            fp6_mul_ref(expected_mul, a_vals[lane], b_vals[lane]);
            fp6_mul_by_xy0_ref(expected_xy0, a_vals[lane], b_xy0_vals[lane]);
            fp6_mul_by_0y0_ref(expected_0y0, a_vals[lane], b_vals[lane][2]);

            const Fp6 got_mul = get_fp6_tile_lane_avx(mul_reg, lane);
            const Fp6 got_xy0 = get_fp6_tile_lane_avx(xy0_reg, lane);
            const Fp6 got_0y0 = get_fp6_tile_lane_avx(y0_reg, lane);

            EXPECT_EQ(fp6_equal_ref(got_mul, expected_mul), true);
            EXPECT_EQ(fp6_equal_ref(got_xy0, expected_xy0), true);
            EXPECT_EQ(fp6_equal_ref(got_0y0, expected_0y0), true);
        }
    }
}

TEST(TileBackend, AvxFp12MulSqrConjMulMatchScalar)
{
    Fp12TileRegAvx512 a_reg;
    Fp12TileRegAvx512 b_reg;
    Fp12TileRegAvx512 mul_reg;
    Fp12TileRegAvx512 sqr_reg;
    Fp12TileRegAvx512 mul_conj_reg;

    Fp12 a_vals[kAvxWidth];
    Fp12 b_vals[kAvxWidth];

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        for (u64 i = 0; i < 2; ++i) {
            for (u64 j = 0; j < 3; ++j) {
                for (u64 k = 0; k < 2; ++k) {
                    Fp x;
                    Fp y;
                    for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                        x[limb] = static_cast<limb_t>(kFp12SeedX + lane * kFp12LaneMulX + i * kFp12IMulX + j * kFp12JMulX + k * kFp12KMulX + limb * kFp12LimbMulX);
                        y[limb] = static_cast<limb_t>(kFp12SeedY + lane * kFp12LaneMulY + i * kFp12IMulY + j * kFp12JMulY + k * kFp12KMulY + limb * kFp12LimbMulY);
                    }
                    while (Fp::is_redc(x) == 0) { Fp::rshift(x, x, 1); }
                    while (Fp::is_redc(y) == 0) { Fp::rshift(y, y, 1); }
                    a_vals[lane][i][j][k] = x;
                    b_vals[lane][i][j][k] = y;
                }
            }
        }
        set_fp12_tile_lane_avx(a_reg, lane, a_vals[lane]);
        set_fp12_tile_lane_avx(b_reg, lane, b_vals[lane]);
    }

    fp12_tile_mul_mont(mul_reg, a_reg, b_reg, get_BLST_P().data, BLST_p0);
    fp12_tile_sqr_mont(sqr_reg, a_reg, get_BLST_P().data, BLST_p0);
    fp12_tile_mul_mont(mul_conj_reg, a_reg, b_reg, true, get_BLST_P().data, BLST_p0);

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        Fp12 expected_mul;
        Fp12 expected_sqr;
        Fp12 expected_mul_conj;
        Fp12::mul(expected_mul, a_vals[lane], b_vals[lane]);
        Fp12::sqr(expected_sqr, a_vals[lane]);
        Fp12::mul(expected_mul_conj, a_vals[lane], b_vals[lane], true);

        const Fp12 got_mul = get_fp12_tile_lane_avx(mul_reg, lane);
        const Fp12 got_sqr = get_fp12_tile_lane_avx(sqr_reg, lane);
        const Fp12 got_mul_conj = get_fp12_tile_lane_avx(mul_conj_reg, lane);

        for (u64 i = 0; i < 2; ++i) {
            for (u64 j = 0; j < 3; ++j) {
                for (u64 k = 0; k < 2; ++k) {
                    for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                        EXPECT_EQ(got_mul[i][j][k][limb], expected_mul[i][j][k][limb]);
                        EXPECT_EQ(got_sqr[i][j][k][limb], expected_sqr[i][j][k][limb]);
                        EXPECT_EQ(got_mul_conj[i][j][k][limb], expected_mul_conj[i][j][k][limb]);
                    }
                }
            }
        }
    }
}

TEST(TileBackend, AvxFp12MulSqrConjMulRandomizedMatchesScalar)
{
    constexpr u64 kIterations = 32;
    u64 state = 0xA5A5F0F00D0D1234ULL;
    Fp12TileRegAvx512 a_reg;
    Fp12TileRegAvx512 b_reg;
    Fp12TileRegAvx512 mul_reg;
    Fp12TileRegAvx512 sqr_reg;
    Fp12TileRegAvx512 mul_conj_reg;

    for (u64 iter = 0; iter < kIterations; ++iter) {
        Fp12 a_vals[kAvxWidth];
        Fp12 b_vals[kAvxWidth];

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            fill_random_redc_fp12_any(a_vals[lane], state);
            fill_random_redc_fp12_any(b_vals[lane], state);
            set_fp12_tile_lane_avx(a_reg, lane, a_vals[lane]);
            set_fp12_tile_lane_avx(b_reg, lane, b_vals[lane]);
        }

        fp12_tile_mul_mont(mul_reg, a_reg, b_reg, get_BLST_P().data, BLST_p0);
        fp12_tile_sqr_mont(sqr_reg, a_reg, get_BLST_P().data, BLST_p0);
        fp12_tile_mul_mont(mul_conj_reg, a_reg, b_reg, true, get_BLST_P().data, BLST_p0);

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            Fp12 expected_mul;
            Fp12 expected_sqr;
            Fp12 expected_mul_conj;
            Fp12::mul(expected_mul, a_vals[lane], b_vals[lane]);
            Fp12::sqr(expected_sqr, a_vals[lane]);
            Fp12::mul(expected_mul_conj, a_vals[lane], b_vals[lane], true);

            const Fp12 got_mul = get_fp12_tile_lane_avx(mul_reg, lane);
            const Fp12 got_sqr = get_fp12_tile_lane_avx(sqr_reg, lane);
            const Fp12 got_mul_conj = get_fp12_tile_lane_avx(mul_conj_reg, lane);

            EXPECT_EQ(Fp12::equal(got_mul, expected_mul), 1);
            EXPECT_EQ(Fp12::equal(got_sqr, expected_sqr), 1);
            EXPECT_EQ(Fp12::equal(got_mul_conj, expected_mul_conj), 1);
        }
    }
}

TEST(TileBackend, AvxFp12SparseMulRandomizedAndInPlaceMatchesScalar)
{
    constexpr u64 kIterations = 32;
    u64 state = 0x13579BDF2468ACE0ULL;
    Fp12TileRegAvx512 a_reg;
    Fp12TileRegAvx512 out_reg;
    Fp12TileRegAvx512 in_place_reg;
    Fp6TileRegAvx512 sparse_reg;

    for (u64 iter = 0; iter < kIterations; ++iter) {
        Fp12 a_vals[kAvxWidth];
        Fp6 sparse_vals[kAvxWidth];

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            fill_random_redc_fp12_any(a_vals[lane], state);
            fill_random_redc_fp6_any(sparse_vals[lane], state);
            set_fp12_tile_lane_avx(a_reg, lane, a_vals[lane]);
            set_fp6_tile_lane_avx(sparse_reg, lane, sparse_vals[lane]);
        }

        fp12_tile_mul_by_xy00z0(out_reg, a_reg, sparse_reg);
        in_place_reg = a_reg;
        fp12_tile_mul_by_xy00z0(in_place_reg, in_place_reg, sparse_reg);

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            Fp12 expected;
            Fp12::mul_by_xy00z0(expected, a_vals[lane], sparse_vals[lane]);

            const Fp12 got = get_fp12_tile_lane_avx(out_reg, lane);
            const Fp12 got_in_place = get_fp12_tile_lane_avx(in_place_reg, lane);

            EXPECT_EQ(Fp12::equal(got, expected), 1);
            EXPECT_EQ(Fp12::equal(got_in_place, expected), 1);
        }
    }
}

TEST(TileBackend, AvxFp12MulAliasInPlaceMatchesScalar)
{
    constexpr u64 kIterations = 32;
    u64 state = 0xCAFEBABE12345678ULL;
    Fp12TileRegAvx512 a_reg;
    Fp12TileRegAvx512 b_reg;

    for (u64 iter = 0; iter < kIterations; ++iter) {
        Fp12 a_vals[kAvxWidth];
        Fp12 b_vals[kAvxWidth];
        Fp12TileRegAvx512 out_ref;
        Fp12TileRegAvx512 out_a_alias;
        Fp12TileRegAvx512 out_b_alias;
        Fp12TileRegAvx512 out_conj_ref;
        Fp12TileRegAvx512 out_conj_a_alias;

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            fill_random_redc_fp12_any(a_vals[lane], state);
            fill_random_redc_fp12_any(b_vals[lane], state);
            set_fp12_tile_lane_avx(a_reg, lane, a_vals[lane]);
            set_fp12_tile_lane_avx(b_reg, lane, b_vals[lane]);
        }

        // Baseline outputs
        fp12_tile_mul_mont(out_ref, a_reg, b_reg, get_BLST_P().data, BLST_p0);
        fp12_tile_mul_mont(out_conj_ref, a_reg, b_reg, true, get_BLST_P().data, BLST_p0);

        // out == a
        out_a_alias = a_reg;
        fp12_tile_mul_mont(out_a_alias, out_a_alias, b_reg, get_BLST_P().data, BLST_p0);

        // out == b
        out_b_alias = b_reg;
        fp12_tile_mul_mont(out_b_alias, a_reg, out_b_alias, get_BLST_P().data, BLST_p0);

        // out == a with conj path
        out_conj_a_alias = a_reg;
        fp12_tile_mul_mont(out_conj_a_alias, out_conj_a_alias, b_reg, true, get_BLST_P().data, BLST_p0);

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            const Fp12 expected = get_fp12_tile_lane_avx(out_ref, lane);
            const Fp12 got_a_alias = get_fp12_tile_lane_avx(out_a_alias, lane);
            const Fp12 got_b_alias = get_fp12_tile_lane_avx(out_b_alias, lane);
            const Fp12 expected_conj = get_fp12_tile_lane_avx(out_conj_ref, lane);
            const Fp12 got_conj_a_alias = get_fp12_tile_lane_avx(out_conj_a_alias, lane);

            EXPECT_EQ(Fp12::equal(got_a_alias, expected), 1);
            EXPECT_EQ(Fp12::equal(got_b_alias, expected), 1);
            EXPECT_EQ(Fp12::equal(got_conj_a_alias, expected_conj), 1);
        }
    }
}

TEST(TileBackend, AvxFp12SqrAliasInPlaceMatchesScalar)
{
    constexpr u64 kIterations = 32;
    u64 state = 0x0F0E0D0C0B0A0908ULL;
    Fp12TileRegAvx512 a_reg;

    for (u64 iter = 0; iter < kIterations; ++iter) {
        Fp12 a_vals[kAvxWidth];
        Fp12TileRegAvx512 out_ref;
        Fp12TileRegAvx512 out_alias;

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            fill_random_redc_fp12_any(a_vals[lane], state);
            set_fp12_tile_lane_avx(a_reg, lane, a_vals[lane]);
        }

        fp12_tile_sqr_mont(out_ref, a_reg, get_BLST_P().data, BLST_p0);
        out_alias = a_reg;
        fp12_tile_sqr_mont(out_alias, out_alias, get_BLST_P().data, BLST_p0);

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            const Fp12 expected = get_fp12_tile_lane_avx(out_ref, lane);
            const Fp12 got = get_fp12_tile_lane_avx(out_alias, lane);
            EXPECT_EQ(Fp12::equal(got, expected), 1);
        }
    }
}

TEST(TileBackend, AvxFpTileLoadStoreAndMemZero)
{
    constexpr u64 stride = kTileWidth + kTileMemStridePad;
    constexpr u64 limbs = FpTileMemViewAvx512::kLimbs;

    alignas(64) limb_t in_mem[limbs * stride];
    alignas(64) limb_t out_mem[limbs * stride];
    FpTileRegAvx512 reg;

    for (u64 i = 0; i < limbs * stride; ++i) {
        in_mem[i] = static_cast<limb_t>(0xABC000000000ull + i);
        out_mem[i] = static_cast<limb_t>(0xFFFFFFFFFFFFFFFFull);
    }

    FpTileMemViewAvx512 in_view;
    in_view.limbs = in_mem;
    in_view.limb_stride = stride;
    FpTileMemViewAvx512 out_view;
    out_view.limbs = out_mem;
    out_view.limb_stride = stride;

    fp_tile_load(reg, in_view);
    fp_tile_store(out_view, reg);

    for (u64 limb = 0; limb < limbs; ++limb) {
        for (u64 lane = 0; lane < kTileWidth; ++lane) {
            const u64 idx = limb * stride + lane;
            EXPECT_EQ(out_mem[idx], (in_mem[idx] & kTileMemMask));
        }
    }

    fp_tile_mem_zero(out_view);
    for (u64 limb = 0; limb < limbs; ++limb) {
        for (u64 lane = 0; lane < kTileWidth; ++lane) {
            const u64 idx = limb * stride + lane;
            EXPECT_EQ(out_mem[idx], static_cast<limb_t>(0));
        }
    }
}

TEST(TileBackend, AvxFp2TileLoadStoreAndMemZero)
{
    constexpr u64 stride = kTileWidth + kTileMemStridePad;
    constexpr u64 limbs = FpTileMemViewAvx512::kLimbs;

    alignas(64) limb_t in0[limbs * stride];
    alignas(64) limb_t in1[limbs * stride];
    alignas(64) limb_t out0[limbs * stride];
    alignas(64) limb_t out1[limbs * stride];
    Fp2TileRegAvx512 reg;

    for (u64 i = 0; i < limbs * stride; ++i) {
        in0[i] = static_cast<limb_t>(0x123450000000ull + i);
        in1[i] = static_cast<limb_t>(0x789AB0000000ull + i);
        out0[i] = static_cast<limb_t>(0xFFFFFFFFFFFFFFFFull);
        out1[i] = static_cast<limb_t>(0xFFFFFFFFFFFFFFFFull);
    }

    Fp2TileMemViewAvx512 in_view;
    in_view[0].limbs = in0;
    in_view[0].limb_stride = stride;
    in_view[1].limbs = in1;
    in_view[1].limb_stride = stride;

    Fp2TileMemViewAvx512 out_view;
    out_view[0].limbs = out0;
    out_view[0].limb_stride = stride;
    out_view[1].limbs = out1;
    out_view[1].limb_stride = stride;

    fp2_tile_load(reg, in_view);
    fp2_tile_store(out_view, reg);

    for (u64 limb = 0; limb < limbs; ++limb) {
        for (u64 lane = 0; lane < kTileWidth; ++lane) {
            const u64 idx = limb * stride + lane;
            EXPECT_EQ(out0[idx], (in0[idx] & kTileMemMask));
            EXPECT_EQ(out1[idx], (in1[idx] & kTileMemMask));
        }
    }

    fp2_tile_mem_zero(out_view);
    for (u64 limb = 0; limb < limbs; ++limb) {
        for (u64 lane = 0; lane < kTileWidth; ++lane) {
            const u64 idx = limb * stride + lane;
            EXPECT_EQ(out0[idx], static_cast<limb_t>(0));
            EXPECT_EQ(out1[idx], static_cast<limb_t>(0));
        }
    }
}

TEST(TileBackend, AvxFp12TileMemZero)
{
    constexpr u64 stride = kTileWidth + kTileMemStridePad;
    constexpr u64 limbs = FpTileMemViewAvx512::kLimbs;
    constexpr u64 rows = 2 * 3 * 2 * limbs;

    alignas(64) limb_t mem[rows * stride];
    for (u64 i = 0; i < rows * stride; ++i) {
        mem[i] = static_cast<limb_t>(0xEEEE00000000ull + i);
    }

    Fp12TileMemViewAvx512 view;
    u64 row_base = 0;
    for (u64 i = 0; i < 2; ++i) {
        for (u64 j = 0; j < 3; ++j) {
            for (u64 k = 0; k < 2; ++k) {
                view[i][j][k].limbs = mem + row_base * stride;
                view[i][j][k].limb_stride = stride;
                row_base += limbs;
            }
        }
    }

    fp12_tile_mem_zero(view);

    for (u64 row = 0; row < rows; ++row) {
        for (u64 lane = 0; lane < kTileWidth; ++lane) {
            EXPECT_EQ(mem[row * stride + lane], static_cast<limb_t>(0));
        }
    }
}

TEST(TileBackend, AvxIntrinDispatchPrimitivesMatchScalar)
{
    alignas(64) u64 a_lanes[8];
    alignas(64) u64 b_lanes[8];
    alignas(64) u64 out_lanes[8];
    alignas(64) u64 mem_in[8];
    alignas(64) u64 mem_out[8];
    u64 unaligned_buf[16];
    u64* const unaligned_ptr = unaligned_buf + 1;

    for (u64 lane = 0; lane < 8; ++lane) {
        a_lanes[lane] = kIntrinSeedA + lane * kIntrinLaneStepA;
        b_lanes[lane] = kIntrinSeedB + lane * kIntrinLaneStepB;
        mem_in[lane] = a_lanes[lane] ^ b_lanes[lane];
        mem_out[lane] = 0;
        unaligned_ptr[lane] = b_lanes[lane];
    }

    const cc512i a = make_cc512_from_lanes(a_lanes);
    const cc512i b = make_cc512_from_lanes(b_lanes);

    // setzero/set1
    const cc512i z = _cc512_setzero_si512();
    store_cc512_to_lanes(out_lanes, z);
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], static_cast<u64>(0));
    }

    const cc512i one = _cc512_set1_epi64(1);
    store_cc512_to_lanes(out_lanes, one);
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], static_cast<u64>(1));
    }

    // add/sub/and/or
    store_cc512_to_lanes(out_lanes, _cc512_add_epi64(a, b));
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], a_lanes[lane] + b_lanes[lane]);
    }

    store_cc512_to_lanes(out_lanes, _cc512_sub_epi64(a, b));
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], a_lanes[lane] - b_lanes[lane]);
    }

    store_cc512_to_lanes(out_lanes, _cc512_and_si512(a, b));
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], a_lanes[lane] & b_lanes[lane]);
    }

    store_cc512_to_lanes(out_lanes, _cc512_or_si512(a, b));
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], a_lanes[lane] | b_lanes[lane]);
    }

    store_cc512_to_lanes(out_lanes, _cc512_xor_si512(a, b));
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], a_lanes[lane] ^ b_lanes[lane]);
    }

    // shifts
    constexpr int kShift = 7;
    store_cc512_to_lanes(out_lanes, _cc512_slli_epi64(a, kShift));
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], a_lanes[lane] << kShift);
    }

    store_cc512_to_lanes(out_lanes, _cc512_srli_epi64(a, kShift));
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], a_lanes[lane] >> kShift);
    }

    store_cc512_to_lanes(out_lanes, _cc512_srai_epi64(a, kShift));
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], static_cast<u64>(static_cast<i64>(a_lanes[lane]) >> kShift));
    }

    alignas(64) u64 shift_lanes[8];
    for (u64 lane = 0; lane < 8; ++lane) {
        shift_lanes[lane] = lane + 1;
    }
    const cc512i sh = make_cc512_from_lanes(shift_lanes);
    store_cc512_to_lanes(out_lanes, _cc512_sllv_epi64(a, sh));
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], a_lanes[lane] << (shift_lanes[lane] & 63));
    }

    // load/store aligned + unaligned
    const cc512i la = _cc512_load_si512(static_cast<const void*>(mem_in));
    _cc512_store_si512(static_cast<void*>(mem_out), la);
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(mem_out[lane], mem_in[lane]);
    }

    const cc512i lu = _cc512_loadu_si512(static_cast<const void*>(unaligned_ptr));
    _cc512_storeu_si512(static_cast<void*>(mem_out), lu);
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(mem_out[lane], unaligned_ptr[lane]);
    }

    // set/permute/blend
    const cc512i setv = _cc512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
    store_cc512_to_lanes(out_lanes, setv);
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], lane);
    }

    alignas(64) u64 idx_lanes[8] = {7, 6, 5, 4, 3, 2, 1, 0};
    const cc512i idx = make_cc512_from_lanes(idx_lanes);
    store_cc512_to_lanes(out_lanes, _cc512_permutexvar_epi64(idx, a));
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(out_lanes[lane], a_lanes[idx_lanes[lane] & 7]);
    }

    const u64 blend_mask = 0b10101100;
    store_cc512_to_lanes(out_lanes, _cc512_mask_blend_epi64(blend_mask, a, b));
    for (u64 lane = 0; lane < 8; ++lane) {
        const u64 use_b = (blend_mask >> lane) & 1;
        EXPECT_EQ(out_lanes[lane], use_b ? b_lanes[lane] : a_lanes[lane]);
    }

    // mask add/sub
    const u64 arith_mask = 0b01101001;
    store_cc512_to_lanes(out_lanes, _cc512_mask_add_epi64(a, arith_mask, a, b));
    for (u64 lane = 0; lane < 8; ++lane) {
        const u64 use_add = (arith_mask >> lane) & 1;
        EXPECT_EQ(out_lanes[lane], use_add ? (a_lanes[lane] + b_lanes[lane]) : a_lanes[lane]);
    }

    store_cc512_to_lanes(out_lanes, _cc512_mask_sub_epi64(a, arith_mask, a, b));
    for (u64 lane = 0; lane < 8; ++lane) {
        const u64 use_sub = (arith_mask >> lane) & 1;
        EXPECT_EQ(out_lanes[lane], use_sub ? (a_lanes[lane] - b_lanes[lane]) : a_lanes[lane]);
    }

    // shuffle/permutex immediate and masked variants
    constexpr int kShufImm = 0x1B;
    const cc512i shuf = _cc512_shuffle_epi32<kShufImm>(a);
    alignas(64) u32 in32[16];
    alignas(64) u32 shuf32[16];
    for (u64 lane = 0; lane < 8; ++lane) {
        in32[lane * 2] = static_cast<u32>(a_lanes[lane] & 0xffffffffu);
        in32[lane * 2 + 1] = static_cast<u32>(a_lanes[lane] >> 32);
    }
    store_cc512_to_lanes(out_lanes, shuf);
    for (u64 lane = 0; lane < 8; ++lane) {
        shuf32[lane * 2] = static_cast<u32>(out_lanes[lane] & 0xffffffffu);
        shuf32[lane * 2 + 1] = static_cast<u32>(out_lanes[lane] >> 32);
    }
    for (u64 block = 0; block < 4; ++block) {
        const u64 base = block * 4;
        for (u64 j = 0; j < 4; ++j) {
            const u64 src = (static_cast<u64>(kShufImm) >> (2 * j)) & 3u;
            EXPECT_EQ(shuf32[base + j], in32[base + src]);
        }
    }

    const u64 shuf_mask = 0xA55Au;
    store_cc512_to_lanes(out_lanes, _cc512_maskz_shuffle_epi32<kShufImm>(shuf_mask, a));
    for (u64 lane = 0; lane < 8; ++lane) {
        shuf32[lane * 2] = static_cast<u32>(out_lanes[lane] & 0xffffffffu);
        shuf32[lane * 2 + 1] = static_cast<u32>(out_lanes[lane] >> 32);
    }
    for (u64 block = 0; block < 4; ++block) {
        const u64 base = block * 4;
        for (u64 j = 0; j < 4; ++j) {
            const u64 src = (static_cast<u64>(kShufImm) >> (2 * j)) & 3u;
            const u32 expected = ((shuf_mask >> (base + j)) & 1u) ? in32[base + src] : 0u;
            EXPECT_EQ(shuf32[base + j], expected);
        }
    }

    constexpr int kPermImm = 0x39;
    store_cc512_to_lanes(out_lanes, _cc512_permutex_epi64<kPermImm>(a));
    for (u64 block = 0; block < 2; ++block) {
        const u64 base = block * 4;
        for (u64 lane = 0; lane < 4; ++lane) {
            const u64 src = (static_cast<u64>(kPermImm) >> (2 * lane)) & 3u;
            EXPECT_EQ(out_lanes[base + lane], a_lanes[base + src]);
        }
    }

    const u64 perm_mask = 0b10110101;
    store_cc512_to_lanes(out_lanes, _cc512_maskz_permutex_epi64<kPermImm>(perm_mask, a));
    for (u64 block = 0; block < 2; ++block) {
        const u64 base = block * 4;
        for (u64 lane = 0; lane < 4; ++lane) {
            const u64 src = (static_cast<u64>(kPermImm) >> (2 * lane)) & 3u;
            const u64 expected = ((perm_mask >> (base + lane)) & 1u) ? a_lanes[base + src] : 0u;
            EXPECT_EQ(out_lanes[base + lane], expected);
        }
    }

    // extract helpers
    const cc256i low = _cc512_extracti64x4_epi64<0>(a);
    const cc256i high = _cc512_extracti64x4_epi64<1>(a);
    EXPECT_EQ(_cc256_extract_epi64<0>(low), a_lanes[0]);
    EXPECT_EQ(_cc256_extract_epi64<0>(high), a_lanes[4]);
    EXPECT_EQ(_cc256_extract_epi64<1>(low), a_lanes[1]);
    EXPECT_EQ(_cc256_extract_epi64<1>(high), a_lanes[5]);
    EXPECT_EQ(_cc256_extract_epi64<2>(low), a_lanes[2]);
    EXPECT_EQ(_cc256_extract_epi64<2>(high), a_lanes[6]);
    EXPECT_EQ(_cc256_extract_epi64<3>(low), a_lanes[3]);
    EXPECT_EQ(_cc256_extract_epi64<3>(high), a_lanes[7]);

    // madd52 low/high semantics
    const cc512i mlo = _cc512_madd52lo_epu64(a, a, b);
    const cc512i mhi = _cc512_madd52hi_epu64(a, a, b);
    alignas(64) u64 lo_lanes[8];
    alignas(64) u64 hi_lanes[8];
    store_cc512_to_lanes(lo_lanes, mlo);
    store_cc512_to_lanes(hi_lanes, mhi);
    for (u64 lane = 0; lane < 8; ++lane) {
        EXPECT_EQ(lo_lanes[lane], madd52_lo_ref(a_lanes[lane], a_lanes[lane], b_lanes[lane]));
        EXPECT_EQ(hi_lanes[lane], madd52_hi_ref(a_lanes[lane], a_lanes[lane], b_lanes[lane]));
    }
}

#if defined(BLST_AVX512_PORTABLE) && !defined(BLST_AVX512)
TEST(TileBackend, AvxPortableWords48RoundTripMatchesFp)
{
    u64 state = kIfmaPortableSeed;
    for (u64 iter = 0; iter < kIfmaPortableIterations; ++iter) {
        for (u64 lane = 0; lane < 8; ++lane) {
            Fp x;
            Fp y;
            u64 words[8];
            fill_random_redc_fp(x, state);
            fill_words48_from_fp(words, x);
            for (u64 limb = 0; limb < 8; ++limb) {
                EXPECT_LE(words[limb], kIfmaPortableWordMask);
            }
            words48_to_fp(y, words);
            EXPECT_EQ(Fp::equal(x, y), 1);
        }
    }
}

TEST(TileBackend, AvxPortableConv64And48RoundTrip)
{
    u64 state = kIfmaPortableSeed ^ 0xA5A5A5A5A5A5A5A5ULL;
    for (u64 iter = 0; iter < kIfmaPortableIterations; ++iter) {
        Fp lane_vals[8];
        u64 in_words[8][8] = {};
        cc512i in_reg[8] = {};
        cc512i mid_reg[6] = {};
        cc512i out_reg[8] = {};
        u64 out_words[8][8] = {};

        for (u64 lane = 0; lane < 8; ++lane) {
            u64 words[8];
            fill_random_redc_fp(lane_vals[lane], state);
            fill_words48_from_fp(words, lane_vals[lane]);
            for (u64 limb = 0; limb < 8; ++limb) {
                in_words[limb][lane] = words[limb];
            }
        }

        set_reg_from_words(in_reg, in_words);
        conv_48to64_fp_8x1w(mid_reg, in_reg);
        conv_64to48_fp_8x1w(out_reg, mid_reg);
        get_words_from_reg(out_words, out_reg);

        for (u64 limb = 0; limb < 8; ++limb) {
            for (u64 lane = 0; lane < 8; ++lane) {
                EXPECT_EQ(out_words[limb][lane], in_words[limb][lane]);
            }
        }

        for (u64 lane = 0; lane < 8; ++lane) {
            Fp decoded;
            u64 words[8];
            for (u64 limb = 0; limb < 8; ++limb) {
                words[limb] = out_words[limb][lane];
            }
            words48_to_fp(decoded, words);
            EXPECT_EQ(Fp::equal(decoded, lane_vals[lane]), 1);
        }
    }
}

TEST(TileBackend, AvxPortableIfmaFpOpsMatchScalarRandomized)
{
    u64 state = kIfmaPortableSeed ^ 0xC3C3C3C3C3C3C3C3ULL;
    for (u64 iter = 0; iter < kIfmaPortableIterations; ++iter) {
        Fp a_vals[8];
        Fp b_vals[8];
        Fp exp_add[8];
        Fp exp_sub[8];
        Fp exp_mul[8];

        u64 a_words[8][8] = {};
        u64 b_words[8][8] = {};
        cc512i a_reg[8] = {};
        cc512i b_reg[8] = {};
        cc512i add_reg[8] = {};
        cc512i sub_reg[8] = {};
        cc512i mul_reg[8] = {};
        u64 add_words[8][8] = {};
        u64 sub_words[8][8] = {};
        u64 mul_words[8][8] = {};

        for (u64 lane = 0; lane < 8; ++lane) {
            u64 lane_aw[8];
            u64 lane_bw[8];
            fill_random_redc_fp(a_vals[lane], state);
            fill_random_redc_fp(b_vals[lane], state);
            fill_words48_from_fp(lane_aw, a_vals[lane]);
            fill_words48_from_fp(lane_bw, b_vals[lane]);
            for (u64 limb = 0; limb < 8; ++limb) {
                a_words[limb][lane] = lane_aw[limb];
                b_words[limb][lane] = lane_bw[limb];
            }
            Fp::add(exp_add[lane], a_vals[lane], b_vals[lane]);
            Fp::sub(exp_sub[lane], a_vals[lane], b_vals[lane]);
            Fp::mul(exp_mul[lane], a_vals[lane], b_vals[lane]);
        }

        set_reg_from_words(a_reg, a_words);
        set_reg_from_words(b_reg, b_words);
        add_fp_8x1w(add_reg, a_reg, b_reg);
        sub_fp_8x1w(sub_reg, a_reg, b_reg);
        mul_redc_fp_8x1w(mul_reg, a_reg, b_reg);
        get_words_from_reg(add_words, add_reg);
        get_words_from_reg(sub_words, sub_reg);
        get_words_from_reg(mul_words, mul_reg);

        for (u64 lane = 0; lane < 8; ++lane) {
            Fp got_add;
            Fp got_sub;
            Fp got_mul;
            u64 lane_add[8];
            u64 lane_sub[8];
            u64 lane_mul[8];
            for (u64 limb = 0; limb < 8; ++limb) {
                lane_add[limb] = add_words[limb][lane];
                lane_sub[limb] = sub_words[limb][lane];
                lane_mul[limb] = mul_words[limb][lane];
            }
            words48_to_fp(got_add, lane_add);
            words48_to_fp(got_sub, lane_sub);
            words48_to_fp(got_mul, lane_mul);
            EXPECT_EQ(Fp::equal(got_add, exp_add[lane]), 1);
            EXPECT_EQ(Fp::equal(got_sub, exp_sub[lane]), 1);
            EXPECT_EQ(Fp::equal(got_mul, exp_mul[lane]), 1);
        }
    }
}

TEST(TileBackend, AvxPortableLineDblVecV1RandomizedMatchesScalar)
{
    u64 state = kIfmaPortableSeed ^ 0x1111222233334444ULL;
    for (u64 iter = 0; iter < kLinePortableFuzzIterations; ++iter) {
        Fp2TileRegAvx512 x1_reg;
        Fp2TileRegAvx512 y1_reg;
        Fp2TileRegAvx512 z1_reg;
        Fp6TileRegAvx512 line_reg;
        Fp2TileRegAvx512 x3_reg;
        Fp2TileRegAvx512 y3_reg;
        Fp2TileRegAvx512 z3_reg;
        G2Point q_vals[kAvxWidth];

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            fill_random_point_like(q_vals[lane], state);
            set_fp2_tile_lane_avx(x1_reg, lane, q_vals[lane].X);
            set_fp2_tile_lane_avx(y1_reg, lane, q_vals[lane].Y);
            set_fp2_tile_lane_avx(z1_reg, lane, q_vals[lane].Z);
        }

        line_dbl_vec_v1(line_reg, x3_reg, y3_reg, z3_reg, x1_reg, y1_reg, z1_reg);

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            G2Point expected_t;
            Fp6 expected_line;
            blst::line_dbl(expected_line, expected_t, q_vals[lane]);
            const Fp6 got_line = get_fp6_tile_lane_avx(line_reg, lane);
            const Fp2 got_x3 = get_fp2_tile_lane_avx(x3_reg, lane);
            const Fp2 got_y3 = get_fp2_tile_lane_avx(y3_reg, lane);
            const Fp2 got_z3 = get_fp2_tile_lane_avx(z3_reg, lane);
            for (u64 i = 0; i < 3; ++i) {
                for (u64 k = 0; k < 2; ++k) {
                    EXPECT_EQ(Fp::equal(got_line[i][k], expected_line[i][k]), 1);
                }
            }
            for (u64 k = 0; k < 2; ++k) {
                EXPECT_EQ(Fp::equal(got_x3[k], expected_t.X[k]), 1);
                EXPECT_EQ(Fp::equal(got_y3[k], expected_t.Y[k]), 1);
                EXPECT_EQ(Fp::equal(got_z3[k], expected_t.Z[k]), 1);
            }
        }
    }
}

TEST(TileBackend, AvxPortableLineDblVecV2RandomizedMatchesV1)
{
    u64 state = kIfmaPortableSeed ^ 0x5555666677778888ULL;
    for (u64 iter = 0; iter < kLinePortableFuzzIterations; ++iter) {
        Fp2TileRegAvx512 x1_reg;
        Fp2TileRegAvx512 y1_reg;
        Fp2TileRegAvx512 z1_reg;
        Fp6TileRegAvx512 line_v1;
        Fp6TileRegAvx512 line_v2;
        Fp2TileRegAvx512 x3_v1;
        Fp2TileRegAvx512 y3_v1;
        Fp2TileRegAvx512 z3_v1;
        Fp2TileRegAvx512 x3_v2;
        Fp2TileRegAvx512 y3_v2;
        Fp2TileRegAvx512 z3_v2;

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            G2Point q;
            fill_random_point_like(q, state);
            set_fp2_tile_lane_avx(x1_reg, lane, q.X);
            set_fp2_tile_lane_avx(y1_reg, lane, q.Y);
            set_fp2_tile_lane_avx(z1_reg, lane, q.Z);
        }

        line_dbl_vec_v1(line_v1, x3_v1, y3_v1, z3_v1, x1_reg, y1_reg, z1_reg);
        line_dbl_vec_v2(line_v2, x3_v2, y3_v2, z3_v2, x1_reg, y1_reg, z1_reg);

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            const Fp6 line_expected = get_fp6_tile_lane_avx(line_v1, lane);
            const Fp6 line_got = get_fp6_tile_lane_avx(line_v2, lane);
            const Fp2 x_expected = get_fp2_tile_lane_avx(x3_v1, lane);
            const Fp2 y_expected = get_fp2_tile_lane_avx(y3_v1, lane);
            const Fp2 z_expected = get_fp2_tile_lane_avx(z3_v1, lane);
            const Fp2 x_got = get_fp2_tile_lane_avx(x3_v2, lane);
            const Fp2 y_got = get_fp2_tile_lane_avx(y3_v2, lane);
            const Fp2 z_got = get_fp2_tile_lane_avx(z3_v2, lane);
            for (u64 i = 0; i < 3; ++i) {
                for (u64 k = 0; k < 2; ++k) {
                    EXPECT_EQ(Fp::equal(line_expected[i][k], line_got[i][k]), 1);
                }
            }
            for (u64 k = 0; k < 2; ++k) {
                EXPECT_EQ(Fp::equal(x_expected[k], x_got[k]), 1);
                EXPECT_EQ(Fp::equal(y_expected[k], y_got[k]), 1);
                EXPECT_EQ(Fp::equal(z_expected[k], z_got[k]), 1);
            }
        }
    }
}

TEST(TileBackend, AvxPortableLineAddVecV1RandomizedMatchesScalar)
{
    u64 state = kIfmaPortableSeed ^ 0x9999AAAABBBBCCCCULL;
    for (u64 iter = 0; iter < kLinePortableFuzzIterations; ++iter) {
        Fp2TileRegAvx512 x1_reg;
        Fp2TileRegAvx512 y1_reg;
        Fp2TileRegAvx512 z1_reg;
        Fp2TileRegAvx512 x2_reg;
        Fp2TileRegAvx512 y2_reg;
        Fp6TileRegAvx512 line_reg;
        Fp2TileRegAvx512 x3_reg;
        Fp2TileRegAvx512 y3_reg;
        Fp2TileRegAvx512 z3_reg;
        G2Point r_vals[kAvxWidth];
        G2AffinePoint q_vals[kAvxWidth];

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            fill_random_point_like(r_vals[lane], state);
            fill_random_affine_like(q_vals[lane], state);
            set_fp2_tile_lane_avx(x1_reg, lane, r_vals[lane].X);
            set_fp2_tile_lane_avx(y1_reg, lane, r_vals[lane].Y);
            set_fp2_tile_lane_avx(z1_reg, lane, r_vals[lane].Z);
            set_fp2_tile_lane_avx(x2_reg, lane, q_vals[lane].X);
            set_fp2_tile_lane_avx(y2_reg, lane, q_vals[lane].Y);
        }

        line_add_vec_v1(line_reg, x3_reg, y3_reg, z3_reg, x1_reg, y1_reg, z1_reg, x2_reg, y2_reg);

        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
            G2Point expected_t;
            Fp6 expected_line;
            blst::line_add(expected_line, expected_t, r_vals[lane], q_vals[lane]);
            const Fp6 got_line = get_fp6_tile_lane_avx(line_reg, lane);
            const Fp2 got_x3 = get_fp2_tile_lane_avx(x3_reg, lane);
            const Fp2 got_y3 = get_fp2_tile_lane_avx(y3_reg, lane);
            const Fp2 got_z3 = get_fp2_tile_lane_avx(z3_reg, lane);
            for (u64 i = 0; i < 3; ++i) {
                for (u64 k = 0; k < 2; ++k) {
                    EXPECT_EQ(Fp::equal(got_line[i][k], expected_line[i][k]), 1);
                }
            }
            for (u64 k = 0; k < 2; ++k) {
                EXPECT_EQ(Fp::equal(got_x3[k], expected_t.X[k]), 1);
                EXPECT_EQ(Fp::equal(got_y3[k], expected_t.Y[k]), 1);
                EXPECT_EQ(Fp::equal(got_z3[k], expected_t.Z[k]), 1);
            }
        }
    }
}
#endif

TEST(TileBackend, AvxLineByPx2KernelMatchesScalar)
{
    Fp6TileRegAvx512 line_reg;
    Fp6TileRegAvx512 out_reg_a;
    Fp6TileRegAvx512 out_reg_b;
    FpTileRegAvx512 px2_x_reg;
    FpTileRegAvx512 px2_y_reg;

    Fp6 line_vals[kAvxWidth];
    Fp px2_x_vals[kAvxWidth];
    Fp px2_y_vals[kAvxWidth];

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        for (u64 i = 0; i < 3; ++i) {
            for (u64 k = 0; k < 2; ++k) {
                Fp v;
                for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                    v[limb] = static_cast<limb_t>(
                        kLineSeed + lane * kLineLaneMul + i * kLineIMul + k * kLineKMul + limb * kLineLimbMul);
                }
                while (Fp::is_redc(v) == 0) { Fp::rshift(v, v, 1); }
                line_vals[lane][i][k] = v;
            }
        }

        Fp x;
        Fp y;
        for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
            x[limb] = static_cast<limb_t>(0x21 + lane * 0x0B + limb * 0x03);
            y[limb] = static_cast<limb_t>(0x31 + lane * 0x0D + limb * 0x05);
        }
        while (Fp::is_redc(x) == 0) { Fp::rshift(x, x, 1); }
        while (Fp::is_redc(y) == 0) { Fp::rshift(y, y, 1); }
        px2_x_vals[lane] = x;
        px2_y_vals[lane] = y;

        set_fp6_tile_lane_avx(line_reg, lane, line_vals[lane]);
        set_fp_tile_lane_avx(px2_x_reg, lane, px2_x_vals[lane]);
        set_fp_tile_lane_avx(px2_y_reg, lane, px2_y_vals[lane]);
    }

    line_by_Px2_4x2x1w(out_reg_a, line_reg, px2_x_reg, px2_y_reg);
    line_by_Px2_2x2x2w(out_reg_b, line_reg, px2_x_reg, px2_y_reg);

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        Fp6 expected;
        expected[0] = line_vals[lane][0];
        Fp::mul(expected[1][0], line_vals[lane][1][0], px2_x_vals[lane]);
        Fp::mul(expected[1][1], line_vals[lane][1][1], px2_x_vals[lane]);
        Fp::mul(expected[2][0], line_vals[lane][2][0], px2_y_vals[lane]);
        Fp::mul(expected[2][1], line_vals[lane][2][1], px2_y_vals[lane]);

        const Fp6 got_a = get_fp6_tile_lane_avx(out_reg_a, lane);
        const Fp6 got_b = get_fp6_tile_lane_avx(out_reg_b, lane);
        for (u64 i = 0; i < 3; ++i) {
            for (u64 k = 0; k < 2; ++k) {
                for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                    EXPECT_EQ(got_a[i][k][limb], expected[i][k][limb]);
                    EXPECT_EQ(got_b[i][k][limb], expected[i][k][limb]);
                }
            }
        }
    }
}

TEST(TileBackend, AvxLineByPx2VectorWrapperMatchesKernel)
{
    Fp6TileRegAvx512 line_reg;
    Fp6TileRegAvx512 out_ref;
    Fp6TileRegAvx512 out_vec;
    FpTileRegAvx512 px2_x_reg;
    FpTileRegAvx512 px2_y_reg;

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        Fp6 line_val;
        Fp px2_x;
        Fp px2_y;
        for (u64 i = 0; i < 3; ++i) {
            for (u64 k = 0; k < 2; ++k) {
                for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                    line_val[i][k][limb] = static_cast<limb_t>(0xC1 + lane * 0x09 + i * 0x07 + k * 0x05 + limb * 0x03);
                }
                while (Fp::is_redc(line_val[i][k]) == 0) { Fp::rshift(line_val[i][k], line_val[i][k], 1); }
            }
        }
        for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
            px2_x[limb] = static_cast<limb_t>(0xD1 + lane * 0x0B + limb * 0x03);
            px2_y[limb] = static_cast<limb_t>(0xE1 + lane * 0x0D + limb * 0x05);
        }
        while (Fp::is_redc(px2_x) == 0) { Fp::rshift(px2_x, px2_x, 1); }
        while (Fp::is_redc(px2_y) == 0) { Fp::rshift(px2_y, px2_y, 1); }

        set_fp6_tile_lane_avx(line_reg, lane, line_val);
        set_fp_tile_lane_avx(px2_x_reg, lane, px2_x);
        set_fp_tile_lane_avx(px2_y_reg, lane, px2_y);
    }

    line_by_Px2_4x2x1w(out_ref, line_reg, px2_x_reg, px2_y_reg);
    line_by_Px2_vector(out_vec, line_reg, px2_x_reg, px2_y_reg);

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        const Fp6 expected = get_fp6_tile_lane_avx(out_ref, lane);
        const Fp6 got = get_fp6_tile_lane_avx(out_vec, lane);
        for (u64 i = 0; i < 3; ++i) {
            for (u64 k = 0; k < 2; ++k) {
                EXPECT_EQ(Fp::equal(expected[i][k], got[i][k]), 1);
            }
        }
    }
}

TEST(TileBackend, AvxLineDblVecV1KernelMatchesScalar)
{
    Fp2TileRegAvx512 x1_reg;
    Fp2TileRegAvx512 y1_reg;
    Fp2TileRegAvx512 z1_reg;
    Fp6TileRegAvx512 line_reg;
    Fp2TileRegAvx512 x3_reg;
    Fp2TileRegAvx512 y3_reg;
    Fp2TileRegAvx512 z3_reg;
    G2Point q_vals[kAvxWidth];

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        for (u64 c = 0; c < 2; ++c) {
            for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                q_vals[lane].X[c][limb] = static_cast<limb_t>(0x41 + lane * 0x0B + c * 0x07 + limb * 0x03);
                q_vals[lane].Y[c][limb] = static_cast<limb_t>(0x51 + lane * 0x0D + c * 0x05 + limb * 0x03);
                q_vals[lane].Z[c][limb] = static_cast<limb_t>(0x61 + lane * 0x11 + c * 0x03 + limb * 0x05);
            }
            while (Fp::is_redc(q_vals[lane].X[c]) == 0) { Fp::rshift(q_vals[lane].X[c], q_vals[lane].X[c], 1); }
            while (Fp::is_redc(q_vals[lane].Y[c]) == 0) { Fp::rshift(q_vals[lane].Y[c], q_vals[lane].Y[c], 1); }
            while (Fp::is_redc(q_vals[lane].Z[c]) == 0) { Fp::rshift(q_vals[lane].Z[c], q_vals[lane].Z[c], 1); }
        }
        set_fp2_tile_lane_avx(x1_reg, lane, q_vals[lane].X);
        set_fp2_tile_lane_avx(y1_reg, lane, q_vals[lane].Y);
        set_fp2_tile_lane_avx(z1_reg, lane, q_vals[lane].Z);
    }

    line_dbl_vec_v1(line_reg, x3_reg, y3_reg, z3_reg, x1_reg, y1_reg, z1_reg);

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        G2Point expected_t;
        Fp6 expected_line;
        blst::line_dbl(expected_line, expected_t, q_vals[lane]);
        const Fp6 got_line = get_fp6_tile_lane_avx(line_reg, lane);
        const Fp2 got_x3 = get_fp2_tile_lane_avx(x3_reg, lane);
        const Fp2 got_y3 = get_fp2_tile_lane_avx(y3_reg, lane);
        const Fp2 got_z3 = get_fp2_tile_lane_avx(z3_reg, lane);
        for (u64 i = 0; i < 3; ++i) {
            for (u64 k = 0; k < 2; ++k) {
                EXPECT_EQ(Fp::equal(got_line[i][k], expected_line[i][k]), 1);
            }
        }
        for (u64 k = 0; k < 2; ++k) {
            EXPECT_EQ(Fp::equal(got_x3[k], expected_t.X[k]), 1);
            EXPECT_EQ(Fp::equal(got_y3[k], expected_t.Y[k]), 1);
            EXPECT_EQ(Fp::equal(got_z3[k], expected_t.Z[k]), 1);
        }
    }
}

TEST(TileBackend, AvxLineDblVecV2KernelMatchesV1)
{
    Fp2TileRegAvx512 x1_reg;
    Fp2TileRegAvx512 y1_reg;
    Fp2TileRegAvx512 z1_reg;
    Fp6TileRegAvx512 line_v1;
    Fp6TileRegAvx512 line_v2;
    Fp2TileRegAvx512 x3_v1;
    Fp2TileRegAvx512 y3_v1;
    Fp2TileRegAvx512 z3_v1;
    Fp2TileRegAvx512 x3_v2;
    Fp2TileRegAvx512 y3_v2;
    Fp2TileRegAvx512 z3_v2;

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        Fp2 x1;
        Fp2 y1;
        Fp2 z1;
        for (u64 c = 0; c < 2; ++c) {
            for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                x1[c][limb] = static_cast<limb_t>(0x121 + lane * 0x09 + c * 0x05 + limb * 0x03);
                y1[c][limb] = static_cast<limb_t>(0x151 + lane * 0x0B + c * 0x07 + limb * 0x03);
                z1[c][limb] = static_cast<limb_t>(0x181 + lane * 0x0D + c * 0x05 + limb * 0x05);
            }
            while (Fp::is_redc(x1[c]) == 0) { Fp::rshift(x1[c], x1[c], 1); }
            while (Fp::is_redc(y1[c]) == 0) { Fp::rshift(y1[c], y1[c], 1); }
            while (Fp::is_redc(z1[c]) == 0) { Fp::rshift(z1[c], z1[c], 1); }
        }
        set_fp2_tile_lane_avx(x1_reg, lane, x1);
        set_fp2_tile_lane_avx(y1_reg, lane, y1);
        set_fp2_tile_lane_avx(z1_reg, lane, z1);
    }

    line_dbl_vec_v1(line_v1, x3_v1, y3_v1, z3_v1, x1_reg, y1_reg, z1_reg);
    line_dbl_vec_v2(line_v2, x3_v2, y3_v2, z3_v2, x1_reg, y1_reg, z1_reg);

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        const Fp6 line_expected = get_fp6_tile_lane_avx(line_v1, lane);
        const Fp6 line_got = get_fp6_tile_lane_avx(line_v2, lane);
        const Fp2 x_expected = get_fp2_tile_lane_avx(x3_v1, lane);
        const Fp2 y_expected = get_fp2_tile_lane_avx(y3_v1, lane);
        const Fp2 z_expected = get_fp2_tile_lane_avx(z3_v1, lane);
        const Fp2 x_got = get_fp2_tile_lane_avx(x3_v2, lane);
        const Fp2 y_got = get_fp2_tile_lane_avx(y3_v2, lane);
        const Fp2 z_got = get_fp2_tile_lane_avx(z3_v2, lane);
        for (u64 i = 0; i < 3; ++i) {
            for (u64 k = 0; k < 2; ++k) {
                EXPECT_EQ(Fp::equal(line_expected[i][k], line_got[i][k]), 1);
            }
        }
        for (u64 k = 0; k < 2; ++k) {
            EXPECT_EQ(Fp::equal(x_expected[k], x_got[k]), 1);
            EXPECT_EQ(Fp::equal(y_expected[k], y_got[k]), 1);
            EXPECT_EQ(Fp::equal(z_expected[k], z_got[k]), 1);
        }
    }
}

TEST(TileBackend, AvxLineAddVecV1KernelMatchesScalar)
{
    Fp2TileRegAvx512 x1_reg;
    Fp2TileRegAvx512 y1_reg;
    Fp2TileRegAvx512 z1_reg;
    Fp2TileRegAvx512 x2_reg;
    Fp2TileRegAvx512 y2_reg;
    Fp6TileRegAvx512 line_reg;
    Fp2TileRegAvx512 x3_reg;
    Fp2TileRegAvx512 y3_reg;
    Fp2TileRegAvx512 z3_reg;
    G2Point r_vals[kAvxWidth];
    G2AffinePoint q_vals[kAvxWidth];

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        for (u64 c = 0; c < 2; ++c) {
            for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                r_vals[lane].X[c][limb] = static_cast<limb_t>(0x71 + lane * 0x0B + c * 0x07 + limb * 0x03);
                r_vals[lane].Y[c][limb] = static_cast<limb_t>(0x81 + lane * 0x0D + c * 0x05 + limb * 0x03);
                r_vals[lane].Z[c][limb] = static_cast<limb_t>(0x91 + lane * 0x11 + c * 0x03 + limb * 0x05);
                q_vals[lane].X[c][limb] = static_cast<limb_t>(0xA1 + lane * 0x13 + c * 0x05 + limb * 0x03);
                q_vals[lane].Y[c][limb] = static_cast<limb_t>(0xB1 + lane * 0x17 + c * 0x03 + limb * 0x05);
            }
            while (Fp::is_redc(r_vals[lane].X[c]) == 0) { Fp::rshift(r_vals[lane].X[c], r_vals[lane].X[c], 1); }
            while (Fp::is_redc(r_vals[lane].Y[c]) == 0) { Fp::rshift(r_vals[lane].Y[c], r_vals[lane].Y[c], 1); }
            while (Fp::is_redc(r_vals[lane].Z[c]) == 0) { Fp::rshift(r_vals[lane].Z[c], r_vals[lane].Z[c], 1); }
            while (Fp::is_redc(q_vals[lane].X[c]) == 0) { Fp::rshift(q_vals[lane].X[c], q_vals[lane].X[c], 1); }
            while (Fp::is_redc(q_vals[lane].Y[c]) == 0) { Fp::rshift(q_vals[lane].Y[c], q_vals[lane].Y[c], 1); }
        }
        set_fp2_tile_lane_avx(x1_reg, lane, r_vals[lane].X);
        set_fp2_tile_lane_avx(y1_reg, lane, r_vals[lane].Y);
        set_fp2_tile_lane_avx(z1_reg, lane, r_vals[lane].Z);
        set_fp2_tile_lane_avx(x2_reg, lane, q_vals[lane].X);
        set_fp2_tile_lane_avx(y2_reg, lane, q_vals[lane].Y);
    }

    line_add_vec_v1(line_reg, x3_reg, y3_reg, z3_reg, x1_reg, y1_reg, z1_reg, x2_reg, y2_reg);

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        G2Point expected_t;
        Fp6 expected_line;
        blst::line_add(expected_line, expected_t, r_vals[lane], q_vals[lane]);
        const Fp6 got_line = get_fp6_tile_lane_avx(line_reg, lane);
        const Fp2 got_x3 = get_fp2_tile_lane_avx(x3_reg, lane);
        const Fp2 got_y3 = get_fp2_tile_lane_avx(y3_reg, lane);
        const Fp2 got_z3 = get_fp2_tile_lane_avx(z3_reg, lane);
        for (u64 i = 0; i < 3; ++i) {
            for (u64 k = 0; k < 2; ++k) {
                EXPECT_EQ(Fp::equal(got_line[i][k], expected_line[i][k]), 1);
            }
        }
        for (u64 k = 0; k < 2; ++k) {
            EXPECT_EQ(Fp::equal(got_x3[k], expected_t.X[k]), 1);
            EXPECT_EQ(Fp::equal(got_y3[k], expected_t.Y[k]), 1);
            EXPECT_EQ(Fp::equal(got_z3[k], expected_t.Z[k]), 1);
        }
    }
}

TEST(TileBackend, AvxLineVectorWrappersMatchKernelEntrypoints)
{
    Fp2TileRegAvx512 x1_reg;
    Fp2TileRegAvx512 y1_reg;
    Fp2TileRegAvx512 z1_reg;
    Fp2TileRegAvx512 x2_reg;
    Fp2TileRegAvx512 y2_reg;
    Fp6TileRegAvx512 line_dbl_ref;
    Fp6TileRegAvx512 line_dbl_vec;
    Fp2TileRegAvx512 x3_dbl_ref;
    Fp2TileRegAvx512 y3_dbl_ref;
    Fp2TileRegAvx512 z3_dbl_ref;
    Fp2TileRegAvx512 x3_dbl_vec;
    Fp2TileRegAvx512 y3_dbl_vec;
    Fp2TileRegAvx512 z3_dbl_vec;
    Fp6TileRegAvx512 line_add_ref;
    Fp6TileRegAvx512 line_add_vec;
    Fp2TileRegAvx512 x3_add_ref;
    Fp2TileRegAvx512 y3_add_ref;
    Fp2TileRegAvx512 z3_add_ref;
    Fp2TileRegAvx512 x3_add_vec;
    Fp2TileRegAvx512 y3_add_vec;
    Fp2TileRegAvx512 z3_add_vec;

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        Fp2 x1;
        Fp2 y1;
        Fp2 z1;
        Fp2 x2;
        Fp2 y2;
        for (u64 c = 0; c < 2; ++c) {
            for (u64 limb = 0; limb < kAvxLimbs; ++limb) {
                x1[c][limb] = static_cast<limb_t>(0x201 + lane * 0x09 + c * 0x05 + limb * 0x03);
                y1[c][limb] = static_cast<limb_t>(0x221 + lane * 0x0B + c * 0x07 + limb * 0x05);
                z1[c][limb] = static_cast<limb_t>(0x241 + lane * 0x0D + c * 0x03 + limb * 0x05);
                x2[c][limb] = static_cast<limb_t>(0x261 + lane * 0x0F + c * 0x05 + limb * 0x03);
                y2[c][limb] = static_cast<limb_t>(0x281 + lane * 0x11 + c * 0x03 + limb * 0x05);
            }
            while (Fp::is_redc(x1[c]) == 0) { Fp::rshift(x1[c], x1[c], 1); }
            while (Fp::is_redc(y1[c]) == 0) { Fp::rshift(y1[c], y1[c], 1); }
            while (Fp::is_redc(z1[c]) == 0) { Fp::rshift(z1[c], z1[c], 1); }
            while (Fp::is_redc(x2[c]) == 0) { Fp::rshift(x2[c], x2[c], 1); }
            while (Fp::is_redc(y2[c]) == 0) { Fp::rshift(y2[c], y2[c], 1); }
        }
        set_fp2_tile_lane_avx(x1_reg, lane, x1);
        set_fp2_tile_lane_avx(y1_reg, lane, y1);
        set_fp2_tile_lane_avx(z1_reg, lane, z1);
        set_fp2_tile_lane_avx(x2_reg, lane, x2);
        set_fp2_tile_lane_avx(y2_reg, lane, y2);
    }

    line_dbl_vec_v2(line_dbl_ref, x3_dbl_ref, y3_dbl_ref, z3_dbl_ref, x1_reg, y1_reg, z1_reg);
    line_dbl_vector(line_dbl_vec, x3_dbl_vec, y3_dbl_vec, z3_dbl_vec, x1_reg, y1_reg, z1_reg);
    line_add_vec_v1(line_add_ref, x3_add_ref, y3_add_ref, z3_add_ref, x1_reg, y1_reg, z1_reg, x2_reg, y2_reg);
    line_add_vector(line_add_vec, x3_add_vec, y3_add_vec, z3_add_vec, x1_reg, y1_reg, z1_reg, x2_reg, y2_reg);

    for (u64 lane = 0; lane < kAvxWidth; ++lane) {
        const Fp6 dbl_expected = get_fp6_tile_lane_avx(line_dbl_ref, lane);
        const Fp6 dbl_got = get_fp6_tile_lane_avx(line_dbl_vec, lane);
        const Fp2 dbl_x_expected = get_fp2_tile_lane_avx(x3_dbl_ref, lane);
        const Fp2 dbl_y_expected = get_fp2_tile_lane_avx(y3_dbl_ref, lane);
        const Fp2 dbl_z_expected = get_fp2_tile_lane_avx(z3_dbl_ref, lane);
        const Fp2 dbl_x_got = get_fp2_tile_lane_avx(x3_dbl_vec, lane);
        const Fp2 dbl_y_got = get_fp2_tile_lane_avx(y3_dbl_vec, lane);
        const Fp2 dbl_z_got = get_fp2_tile_lane_avx(z3_dbl_vec, lane);
        const Fp6 add_expected = get_fp6_tile_lane_avx(line_add_ref, lane);
        const Fp6 add_got = get_fp6_tile_lane_avx(line_add_vec, lane);
        const Fp2 add_x_expected = get_fp2_tile_lane_avx(x3_add_ref, lane);
        const Fp2 add_y_expected = get_fp2_tile_lane_avx(y3_add_ref, lane);
        const Fp2 add_z_expected = get_fp2_tile_lane_avx(z3_add_ref, lane);
        const Fp2 add_x_got = get_fp2_tile_lane_avx(x3_add_vec, lane);
        const Fp2 add_y_got = get_fp2_tile_lane_avx(y3_add_vec, lane);
        const Fp2 add_z_got = get_fp2_tile_lane_avx(z3_add_vec, lane);

        for (u64 i = 0; i < 3; ++i) {
            for (u64 k = 0; k < 2; ++k) {
                EXPECT_EQ(Fp::equal(dbl_expected[i][k], dbl_got[i][k]), 1);
                EXPECT_EQ(Fp::equal(add_expected[i][k], add_got[i][k]), 1);
            }
        }

        for (u64 k = 0; k < 2; ++k) {
            EXPECT_EQ(Fp::equal(dbl_x_expected[k], dbl_x_got[k]), 1);
            EXPECT_EQ(Fp::equal(dbl_y_expected[k], dbl_y_got[k]), 1);
            EXPECT_EQ(Fp::equal(dbl_z_expected[k], dbl_z_got[k]), 1);
            EXPECT_EQ(Fp::equal(add_x_expected[k], add_x_got[k]), 1);
            EXPECT_EQ(Fp::equal(add_y_expected[k], add_y_got[k]), 1);
            EXPECT_EQ(Fp::equal(add_z_expected[k], add_z_got[k]), 1);
        }
    }
}

#endif
