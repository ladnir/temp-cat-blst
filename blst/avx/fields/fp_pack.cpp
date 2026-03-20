#include "blst/fields/Fp_avx.h"
#include "blst/fields/Fp.h"
#include <array>
#include <atomic>

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst::avx
{


inline constexpr const std::array< u64, NWORDS> P48X2 = {
    0xFFFFFFFF5556ULL, 0x62A7FFFF73FDULL, 0xEC483D57FFFDULL, 0xCE61A541ED61ULL,
    0x9709E70A257EULL, 0x869759AEC8EEULL, 0xCD3496374F6CULL, 0x340223D472FFULL,
};

inline constexpr const std::array< u64, NWORDS> P48X4 = {
    0xFFFFFFFEAAACULL, 0xC54FFFFEE7FBULL, 0xD8907AAFFFFAULL, 0x9CC34A83DAC3ULL,
    0x2E13CE144AFDULL, 0x0D2EB35D91DDULL, 0x9A692C6E9ED9ULL, 0x680447A8E5FFULL,
};

#ifndef BLST_HAVE_FP8X1_MUL_REDC_CORE_K2_ASM
#define BLST_HAVE_FP8X1_MUL_REDC_CORE_K2_ASM 0
#endif

#if BLST_HAVE_FP8X1_MUL_REDC_CORE_K2_ASM
extern "C" void blst_fp8x1_mul_redc_core_k2(
    LFpPack8x1* out0,
    LFpPack8x1* out1,
    const FpPack8x1* a0,
    const FpPack8x1* b0,
    const FpPack8x1* a1,
    const FpPack8x1* b1);
#endif

namespace
{
#if defined(NDEBUG)
void fp8x1_lfp_bounds_record(const LFpPack8x1&)
{
}
#else
std::atomic<u64> g_fp8x1_lfp_total_lanes{0};
std::atomic<u64> g_fp8x1_lfp_ge_p_lanes{0};
std::atomic<u64> g_fp8x1_lfp_max_limb{0};

BLST_FORCEINLINE bool fp8x1_lane_ge_p(const u64 limbs[NWORDS])
{
    for (int i = static_cast<int>(NWORDS) - 1; i >= 0; --i) {
        if (limbs[i] > P48[i]) {
            return true;
        }
        if (limbs[i] < P48[i]) {
            return false;
        }
    }
    return true; // equal to p
}

void fp8x1_lfp_bounds_record(const LFpPack8x1& v)
{
    std::array<std::array<u64, 8>, NWORDS> words{};
    for (std::size_t i = 0; i < NWORDS; ++i) {
        _cc512_storeu_si512(words[i].data(), v[i]);
    }

    u64 local_max_limb = 0;
    u64 ge_p_count = 0;
    for (std::size_t lane = 0; lane < 8; ++lane) {
        u64 lane_limbs[NWORDS];
        for (std::size_t i = 0; i < NWORDS; ++i) {
            const u64 limb = words[i][lane];
            BLST_ASSERT(limb <= BMASK);
            lane_limbs[i] = limb;
            if (limb > local_max_limb) {
                local_max_limb = limb;
            }
        }
        ge_p_count += fp8x1_lane_ge_p(lane_limbs) ? 1u : 0u;
    }

    g_fp8x1_lfp_total_lanes.fetch_add(8, std::memory_order_relaxed);
    g_fp8x1_lfp_ge_p_lanes.fetch_add(ge_p_count, std::memory_order_relaxed);

    u64 prev = g_fp8x1_lfp_max_limb.load(std::memory_order_relaxed);
    while (prev < local_max_limb &&
           !g_fp8x1_lfp_max_limb.compare_exchange_weak(
               prev, local_max_limb, std::memory_order_relaxed, std::memory_order_relaxed)) {
    }
}
#endif

BLST_FORCEINLINE i64 lufp_floor_div_radix(const i64 x)
{
    if (x >= 0) {
        return x >> BRADIX;
    }
    return -static_cast<i64>((static_cast<u64>(-x) + BMASK) >> BRADIX);
}

BLST_FORCEINLINE void lufp_scalar_signed_renorm_range(i64* limbs, const std::size_t begin, const std::size_t end)
{
    for (std::size_t i = begin; i < end; ++i) {
        const i64 carry = lufp_floor_div_radix(limbs[i]);
        limbs[i] -= (carry << BRADIX);
        limbs[i + 1] += carry;
    }
}

BLST_FORCEINLINE bool lufp_scalar_high_negative(const i64* limbs)
{
    return limbs[15] < 0;
}

BLST_FORCEINLINE bool lufp_scalar_high_ge_p(const i64* limbs)
{
    for (int i = static_cast<int>(NWORDS) - 1; i >= 0; --i) {
        const i64 limb = limbs[8 + i];
        const i64 mod = static_cast<i64>(P48[i]);
        if (limb > mod) {
            return true;
        }
        if (limb < mod) {
            return false;
        }
    }
    return true;
}

BLST_FORCEINLINE void lufp_scalar_high_add_p(i64* limbs)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        limbs[8 + i] += static_cast<i64>(P48[i]);
    }
    lufp_scalar_signed_renorm_range(limbs, 8, 15);
}

BLST_FORCEINLINE void lufp_scalar_high_sub_p(i64* limbs)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        limbs[8 + i] -= static_cast<i64>(P48[i]);
    }
    lufp_scalar_signed_renorm_range(limbs, 8, 15);
}

BLST_FORCEINLINE bool lufp_scalar_high_ge_mod(const i64* limbs, const u64* mod)
{
    for (int i = static_cast<int>(NWORDS) - 1; i >= 0; --i) {
        const i64 lhs = limbs[8 + i];
        const i64 rhs = static_cast<i64>(mod[i]);
        if (lhs > rhs) {
            return true;
        }
        if (lhs < rhs) {
            return false;
        }
    }
    return true;
}

BLST_FORCEINLINE void lufp_scalar_high_add_mod(i64* limbs, const u64* mod)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        limbs[8 + i] += static_cast<i64>(mod[i]);
    }
    lufp_scalar_signed_renorm_range(limbs, 8, 15);
}

BLST_FORCEINLINE void lufp_scalar_high_sub_mod(i64* limbs, const u64* mod)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        limbs[8 + i] -= static_cast<i64>(mod[i]);
    }
    lufp_scalar_signed_renorm_range(limbs, 8, 15);
}

BLST_FORCEINLINE void lufp_scalar_normalize_to_ufp(i64* limbs)
{
    lufp_scalar_signed_renorm_range(limbs, 0, 15);

    int guard = 0;
    while (lufp_scalar_high_negative(limbs)) {
        lufp_scalar_high_add_p(limbs);
        if (++guard > 64) {
            break;
        }
    }

    guard = 0;
    while (lufp_scalar_high_ge_p(limbs)) {
        lufp_scalar_high_sub_p(limbs);
        if (++guard > 64) {
            break;
        }
    }

    BLST_ASSERT(!lufp_scalar_high_negative(limbs));
    BLST_ASSERT(!lufp_scalar_high_ge_p(limbs));
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        BLST_ASSERT(limbs[i] >= 0);
        BLST_ASSERT(static_cast<u64>(limbs[i]) <= BMASK);
    }
}

template <class R, class A, class B, class Op>
BLST_FORCEINLINE void lufp8x1_binary_op(R& r, const A& a, const B& b, Op&& op)
{
    BLST_ASSERT_NO_OVERLAP_AVX(r, a);
    BLST_ASSERT_NO_OVERLAP_AVX(r, b);
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        r[i] = op(a[i], b[i]);
    }
}

BLST_FORCEINLINE void lufp8x1_signed_carry_range(cc512i* limbs, const std::size_t begin, const std::size_t end)
{
    for (std::size_t i = begin; i < end; ++i) {
        const cc512i carry = VSRA(limbs[i], BRADIX);
        limbs[i] = VSUB(limbs[i], VSHL(carry, BRADIX));
        limbs[i + 1] = VADD(limbs[i + 1], carry);
    }
}

BLST_FORCEINLINE auto lufp8x1_mask_eq_epi64(const cc512i a, const cc512i b) -> u8
{
#if BLST_AVX_NATIVE
    return static_cast<u8>(_mm512_cmpeq_epi64_mask(to_native_512(a), to_native_512(b)));
#else
    std::array<u64, 8> lhs{};
    std::array<u64, 8> rhs{};
    _cc512_storeu_si512(lhs.data(), a);
    _cc512_storeu_si512(rhs.data(), b);
    u8 mask = 0;
    for (u8 lane = 0; lane < 8; ++lane) {
        if (lhs[lane] == rhs[lane]) {
            mask |= static_cast<u8>(1u << lane);
        }
    }
    return mask;
#endif
}

BLST_FORCEINLINE auto lufp8x1_mask_gt_epi64(const cc512i a, const cc512i b) -> u8
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

BLST_FORCEINLINE auto lufp8x1_high_ge_p_mask(const cc512i* limbs) -> u8
{
    u8 eq_mask = 0xFF;
    u8 gt_mask = 0;

    for (int i = static_cast<int>(NWORDS) - 1; i >= 0; --i) {
        const cc512i mod = VSET1(static_cast<i64>(P48[i]));
        const u8 gt_i = lufp8x1_mask_gt_epi64(limbs[8 + i], mod);
        const u8 eq_i = lufp8x1_mask_eq_epi64(limbs[8 + i], mod);
        gt_mask |= static_cast<u8>(eq_mask & gt_i);
        eq_mask &= eq_i;
    }

    return static_cast<u8>(gt_mask | eq_mask);
}

BLST_FORCEINLINE auto lufp8x1_high_ge_mod_mask(const cc512i* limbs, const u64* mod) -> u8
{
    u8 eq_mask = 0xFF;
    u8 gt_mask = 0;

    for (int i = static_cast<int>(NWORDS) - 1; i >= 0; --i) {
        const cc512i mod_i = VSET1(static_cast<i64>(mod[i]));
        const u8 gt_i = lufp8x1_mask_gt_epi64(limbs[8 + i], mod_i);
        const u8 eq_i = lufp8x1_mask_eq_epi64(limbs[8 + i], mod_i);
        gt_mask |= static_cast<u8>(eq_mask & gt_i);
        eq_mask &= eq_i;
    }

    return static_cast<u8>(gt_mask | eq_mask);
}

BLST_FORCEINLINE void lufp8x1_high_add_p_masked(cc512i* limbs, const u8 lane_mask)
{
    limbs[8] = VMADD(limbs[8], lane_mask, limbs[8], VSET1(static_cast<i64>(P48[0])));
    limbs[9] = VMADD(limbs[9], lane_mask, limbs[9], VSET1(static_cast<i64>(P48[1])));
    limbs[10] = VMADD(limbs[10], lane_mask, limbs[10], VSET1(static_cast<i64>(P48[2])));
    limbs[11] = VMADD(limbs[11], lane_mask, limbs[11], VSET1(static_cast<i64>(P48[3])));
    limbs[12] = VMADD(limbs[12], lane_mask, limbs[12], VSET1(static_cast<i64>(P48[4])));
    limbs[13] = VMADD(limbs[13], lane_mask, limbs[13], VSET1(static_cast<i64>(P48[5])));
    limbs[14] = VMADD(limbs[14], lane_mask, limbs[14], VSET1(static_cast<i64>(P48[6])));
    limbs[15] = VMADD(limbs[15], lane_mask, limbs[15], VSET1(static_cast<i64>(P48[7])));
}

BLST_FORCEINLINE void lufp8x1_high_sub_p_masked(cc512i* limbs, const u8 lane_mask)
{
    limbs[8] = VMSUB(limbs[8], lane_mask, limbs[8], VSET1(static_cast<i64>(P48[0])));
    limbs[9] = VMSUB(limbs[9], lane_mask, limbs[9], VSET1(static_cast<i64>(P48[1])));
    limbs[10] = VMSUB(limbs[10], lane_mask, limbs[10], VSET1(static_cast<i64>(P48[2])));
    limbs[11] = VMSUB(limbs[11], lane_mask, limbs[11], VSET1(static_cast<i64>(P48[3])));
    limbs[12] = VMSUB(limbs[12], lane_mask, limbs[12], VSET1(static_cast<i64>(P48[4])));
    limbs[13] = VMSUB(limbs[13], lane_mask, limbs[13], VSET1(static_cast<i64>(P48[5])));
    limbs[14] = VMSUB(limbs[14], lane_mask, limbs[14], VSET1(static_cast<i64>(P48[6])));
    limbs[15] = VMSUB(limbs[15], lane_mask, limbs[15], VSET1(static_cast<i64>(P48[7])));
}

BLST_FORCEINLINE void lufp8x1_high_add_mod_masked(cc512i* limbs, const u8 lane_mask, const u64* mod)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        limbs[8 + i] = VMADD(
            limbs[8 + i], lane_mask, limbs[8 + i], VSET1(static_cast<i64>(mod[i])));
    }
}

BLST_FORCEINLINE void lufp8x1_high_sub_mod_masked(cc512i* limbs, const u8 lane_mask, const u64* mod)
{
    for (std::size_t i = 0; i < NWORDS; ++i) {
        limbs[8 + i] = VMSUB(
            limbs[8 + i], lane_mask, limbs[8 + i], VSET1(static_cast<i64>(mod[i])));
    }
}

#if !defined(NDEBUG)
void lufp8x1_debug_assert_hi_neg1_pos2(const cc512i* limbs)
{
    std::array<std::array<u64, 8>, 2 * NWORDS> words{};
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        _cc512_storeu_si512(words[i].data(), limbs[i]);
    }

    for (std::size_t lane = 0; lane < 8; ++lane) {
        i64 work[2 * NWORDS]{};
        for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
            work[i] = static_cast<i64>(words[i][lane]);
        }

        i64 lower[2 * NWORDS]{};
        i64 upper[2 * NWORDS]{};
        for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
            lower[i] = work[i];
            upper[i] = work[i];
        }

        lufp_scalar_high_add_p(lower);
        BLST_ASSERT(!lufp_scalar_high_negative(lower));

        if (lufp_scalar_high_ge_p(upper)) {
            lufp_scalar_high_sub_p(upper);
        }
        BLST_ASSERT(!lufp_scalar_high_ge_p(upper));
    }
}

void lufp8x1_debug_assert_hi_neg4_pos4(const cc512i* limbs)
{
    std::array<std::array<u64, 8>, 2 * NWORDS> words{};
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
        _cc512_storeu_si512(words[i].data(), limbs[i]);
    }

    for (std::size_t lane = 0; lane < 8; ++lane) {
        i64 work[2 * NWORDS]{};
        for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
            work[i] = static_cast<i64>(words[i][lane]);
        }

        i64 lower[2 * NWORDS]{};
        i64 upper[2 * NWORDS]{};
        for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
            lower[i] = work[i];
            upper[i] = work[i];
        }

        lufp_scalar_high_add_mod(lower, P48X4.data());
        if (lufp_scalar_high_negative(lower)) {
            lufp_scalar_high_add_p(lower);
        }
        BLST_ASSERT(!lufp_scalar_high_negative(lower));

        if (lufp_scalar_high_ge_mod(upper, P48X4.data())) {
            lufp_scalar_high_sub_mod(upper, P48X4.data());
        }
        if (lufp_scalar_high_ge_mod(upper, P48X2.data())) {
            lufp_scalar_high_sub_mod(upper, P48X2.data());
        }
        if (lufp_scalar_high_ge_p(upper)) {
            lufp_scalar_high_sub_p(upper);
        }
        BLST_ASSERT(!lufp_scalar_high_negative(lower));
        BLST_ASSERT(!lufp_scalar_high_ge_p(upper));
    }
}
#endif
} // namespace

void conv_64to48_fp(u64 *r, const u64 *a)
{
  r[0] =                  a[0]         & BMASK;
  r[1] = ((a[0] >> 48) | (a[1] << 16)) & BMASK;
  r[2] = ((a[1] >> 32) | (a[2] << 32)) & BMASK;
  r[3] =  (a[2] >> 16)                 & BMASK;
  r[4] =                  a[3]         & BMASK;
  r[5] = ((a[3] >> 48) | (a[4] << 16)) & BMASK;
  r[6] = ((a[4] >> 32) | (a[5] << 32)) & BMASK;
  r[7] =  (a[5] >> 16)                 & BMASK;
}

void FpPack8x1::conv_64to48(FpPack8x1& r, const cc512i* BLST_restrict a)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2];
  const cc512i a3 = a[3], a4 = a[4], a5 = a[5];
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7;

  r0 = VAND(                       a0,       bmask);
  r1 = VAND(VOR(VSHR(a0, 48), VSHL(a1, 16)), bmask);
  r2 = VAND(VOR(VSHR(a1, 32), VSHL(a2, 32)), bmask);
  r3 = VAND(    VSHR(a2, 16),                bmask);
  r4 = VAND(                       a3,       bmask);
  r5 = VAND(VOR(VSHR(a3, 48), VSHL(a4, 16)), bmask);
  r6 = VAND(VOR(VSHR(a4, 32), VSHL(a5, 32)), bmask);
  r7 = VAND(    VSHR(a5, 16),                bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
}

void FpPack4x2::conv_64to48(FpPack4x2& r, const cc512i* BLST_restrict a)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2];
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2, r3;

  r0 = VAND(                       a0,       bmask);
  r1 = VAND(VOR(VSHR(a0, 48), VSHL(a1, 16)), bmask);
  r2 = VAND(VOR(VSHR(a1, 32), VSHL(a2, 32)), bmask);
  r3 = VAND(    VSHR(a2, 16),                bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
}

void conv_64to48_mpi(u64 *r, const u64 *a, int rlen, int alen)
{
  int i, j, shr_pos, shl_pos;
  u64 word, temp;

  i = j = 0;
  shr_pos = 64; shl_pos = 0;
  temp = 0;
  while ((i < rlen) && (j < alen)) {
    word = ((temp >> shr_pos) | (a[j] << shl_pos));
    r[i] = (word & BMASK);
    shr_pos -= 16, shl_pos += 16;
    if ((shr_pos > 0) && (shl_pos < 64)) temp = a[j++];
    if (shr_pos <= 0) shr_pos += 64;
    if (shl_pos >= 64) shl_pos -= 64;
    // Any shift past 63 is undefined!
    if (shr_pos == 64) temp = 0;
    i++;
  }
  if (i < rlen) r[i++] = ((temp >> shr_pos) & BMASK);
  for (; i < rlen; i++) r[i] = 0;
}

void conv_48to64_fp(u64 *r, const u64 *a)
{
  r[0] = (a[1] << 48) | a[0];
  r[1] = (a[2] << 32) | a[1] >> 16;
  r[2] = (a[3] << 16) | a[2] >> 32;
  r[3] = (a[5] << 48) | a[4];
  r[4] = (a[6] << 32) | a[5] >> 16;
  r[5] = (a[7] << 16) | a[6] >> 32;  
}

void FpPack8x1::conv_48to64(cc512i* BLST_restrict r, const FpPack8x1& a)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  cc512i r0, r1, r2, r3, r4, r5;

  r0 = VOR(VSHL(a1, 48),      a0     );
  r1 = VOR(VSHL(a2, 32), VSHR(a1, 16));
  r2 = VOR(VSHL(a3, 16), VSHR(a2, 32));
  r3 = VOR(VSHL(a5, 48),      a4     );
  r4 = VOR(VSHL(a6, 32), VSHR(a5, 16));
  r5 = VOR(VSHL(a7, 16), VSHR(a6, 32));

  r[0] = r0; r[1] = r1; r[2] = r2; 
  r[3] = r3; r[4] = r4; r[5] = r5; 
}

void FpPack8x1::from_aos(FpPack8x1& r, const Fp* in8)
{
  BLST_ASSERT(in8 != nullptr);

  cc512i in64[SWORDS];
  for (u64 limb = 0; limb < SWORDS; ++limb) {
    std::array<u64, 8> lane_words{};
    for (u64 lane = 0; lane < 8; ++lane) {
      lane_words[lane] = static_cast<u64>(in8[lane][limb]);
    }
    in64[limb] = _cc512_loadu_si512(lane_words.data());
  }

  conv_64to48(r, in64);
}

void FpPack8x1::to_aos(Fp* out8, const FpPack8x1& a)
{
  BLST_ASSERT(out8 != nullptr);

  cc512i out64[SWORDS];
  conv_48to64(out64, a);

  for (u64 limb = 0; limb < SWORDS; ++limb) {
    std::array<u64, 8> lane_words{};
    _cc512_storeu_si512(lane_words.data(), out64[limb]);
    for (u64 lane = 0; lane < 8; ++lane) {
      out8[lane][limb] = lane_words[lane];
    }
  }
}

void FpPack8x1::from_aos_partial(FpPack8x1& r, const Fp* in, std::size_t n_active, const Fp& fill)
{
  BLST_ASSERT(n_active <= 8);
  BLST_ASSERT(in != nullptr || n_active == 0);

  std::array<Fp, 8> tmp{};
  for (std::size_t i = 0; i < n_active; ++i) {
    tmp[i] = in[i];
  }
  for (std::size_t i = n_active; i < 8; ++i) {
    tmp[i] = fill;
  }

  from_aos(r, tmp.data());
}

void FpPack8x1::to_aos_partial(Fp* out, std::size_t n_active, const FpPack8x1& a)
{
  BLST_ASSERT(n_active <= 8);
  BLST_ASSERT(out != nullptr || n_active == 0);

  std::array<Fp, 8> tmp{};
  to_aos(tmp.data(), a);
  for (std::size_t i = 0; i < n_active; ++i) {
    out[i] = tmp[i];
  }
}

void FpPack8x1::select(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b, const u8 lane_mask)
{
  const u64 mask = static_cast<u64>(lane_mask);
  for (u64 i = 0; i < NWORDS; ++i) {
    r[i] = VMBLEND(mask, a[i], b[i]);
  }
}

void FpPack4x2::conv_48to64(cc512i* BLST_restrict r, const FpPack4x2& a)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  cc512i r0, r1, r2;

  r0 = VOR(VSHL(a1, 48),      a0     );
  r1 = VOR(VSHL(a2, 32), VSHR(a1, 16));
  r2 = VOR(VSHL(a3, 16), VSHR(a2, 32));

  r[0] = r0; r[1] = r1; r[2] = r2; 
}

void UFpPack8x1::select(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b, const u8 lane_mask)
{
  const u64 mask = static_cast<u64>(lane_mask);
  for (u64 i = 0; i < 2 * NWORDS; ++i) {
    r[i] = VMBLEND(mask, a[i], b[i]);
  }
}

void LUFpPack8x1::select(LUFpPack8x1& r, const LUFpPack8x1& a, const LUFpPack8x1& b, const u8 lane_mask)
{
  const u64 mask = static_cast<u64>(lane_mask);
  for (u64 i = 0; i < 2 * NWORDS; ++i) {
    r[i] = VMBLEND(mask, a[i], b[i]);
  }
}

void LUFpPack8x1::normalize(UFpPack8x1& r, const LUFpPack8x1& a)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);

  std::array<std::array<u64, 8>, 2 * NWORDS> words{};
  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    _cc512_storeu_si512(words[i].data(), a[i]);
  }

  std::array<std::array<u64, 8>, 2 * NWORDS> out_words{};
  for (std::size_t lane = 0; lane < 8; ++lane) {
    i64 limbs[2 * NWORDS];
    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
      limbs[i] = static_cast<i64>(words[i][lane]);
    }

    lufp_scalar_normalize_to_ufp(limbs);

    for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
      out_words[i][lane] = static_cast<u64>(limbs[i]);
    }
  }

  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    r[i] = _cc512_loadu_si512(out_words[i].data());
  }
}

void LUFpPack8x1::normalize_hi_neg1_pos2(UFpPack8x1& r, const LUFpPack8x1& a)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);

  cc512i limbs[2 * NWORDS];
  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    limbs[i] = a[i];
  }

  lufp8x1_signed_carry_range(limbs, 0, (2 * NWORDS) - 1);

#if !defined(NDEBUG)
  lufp8x1_debug_assert_hi_neg1_pos2(limbs);
#endif

  const u8 neg_mask = lufp8x1_mask_gt_epi64(VZERO(), limbs[15]);
  lufp8x1_high_add_p_masked(limbs, neg_mask);
  lufp8x1_signed_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  const u8 ge_mask = lufp8x1_high_ge_p_mask(limbs);
  lufp8x1_high_sub_p_masked(limbs, ge_mask);
  lufp8x1_signed_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  const cc512i bmask = VSET1(BMASK);
  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    r[i] = VAND(limbs[i], bmask);
  }

#if !defined(NDEBUG)
  UFpPack8x1 ref{};
  LUFpPack8x1::normalize(ref, a);
  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    std::array<u64, 8> got{};
    std::array<u64, 8> want{};
    _cc512_storeu_si512(got.data(), r[i]);
    _cc512_storeu_si512(want.data(), ref[i]);
    BLST_ASSERT(got == want);
  }
#endif
}

void LUFpPack8x1::normalize_hi_neg4_pos4(UFpPack8x1& r, const LUFpPack8x1& a)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);

  cc512i limbs[2 * NWORDS];
  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    limbs[i] = a[i];
  }

  lufp8x1_signed_carry_range(limbs, 0, (2 * NWORDS) - 1);

#if !defined(NDEBUG)
  lufp8x1_debug_assert_hi_neg4_pos4(limbs);
#endif

  const u8 neg_mask = lufp8x1_mask_gt_epi64(VZERO(), limbs[15]);
  lufp8x1_high_add_mod_masked(limbs, neg_mask, P48X4.data());
  lufp8x1_signed_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  const u8 neg_fix_mask = lufp8x1_mask_gt_epi64(VZERO(), limbs[15]);
  lufp8x1_high_add_p_masked(limbs, neg_fix_mask);
  lufp8x1_signed_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  const u8 ge_4p_mask = lufp8x1_high_ge_mod_mask(limbs, P48X4.data());
  lufp8x1_high_sub_mod_masked(limbs, ge_4p_mask, P48X4.data());
  lufp8x1_signed_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  const u8 ge_2p_mask = lufp8x1_high_ge_mod_mask(limbs, P48X2.data());
  lufp8x1_high_sub_mod_masked(limbs, ge_2p_mask, P48X2.data());
  lufp8x1_signed_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  const u8 ge_p_mask = lufp8x1_high_ge_p_mask(limbs);
  lufp8x1_high_sub_p_masked(limbs, ge_p_mask);
  lufp8x1_signed_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  const cc512i bmask = VSET1(BMASK);
  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    r[i] = VAND(limbs[i], bmask);
  }

#if !defined(NDEBUG)
  UFpPack8x1 ref{};
  LUFpPack8x1::normalize(ref, a);
  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    std::array<u64, 8> got{};
    std::array<u64, 8> want{};
    _cc512_storeu_si512(got.data(), r[i]);
    _cc512_storeu_si512(want.data(), ref[i]);
    BLST_ASSERT(got == want);
  }
#endif
}

void conv_48to64_mpi(u64 *r, const u64 *a, int rlen, int alen)
{
  int i, j, bits_in_word, bits_to_shift;
  u64 word;

  i = j = 0;
  bits_in_word = bits_to_shift = 0;
  word = 0;
  while ((i < rlen) && (j < alen)) {
    word |= (a[j] << bits_in_word);
    bits_to_shift = (64 - bits_in_word);
    bits_in_word += 48;
    if (bits_in_word >= 64) {
      r[i++] = word;
      word = ((bits_to_shift > 0) ? (a[j] >> bits_to_shift) : 0);
      bits_in_word = ((bits_to_shift > 0) ? (48 - bits_to_shift) : 0);
    }
    j++;
  }
  if (i < rlen) r[i++] = word;
  for (; i < rlen; i++) r[i] = 0;
}

void carryp_mpi48(u64 *a)
{
  int i;

  for (i = 0; i < NWORDS-1; i++) {
    a[i+1] += a[i]>>BRADIX;
    a[i] &= BMASK;
  }
}

void fp8x1_lfp_bounds_reset()
{
#if !defined(NDEBUG)
    g_fp8x1_lfp_total_lanes.store(0, std::memory_order_relaxed);
    g_fp8x1_lfp_ge_p_lanes.store(0, std::memory_order_relaxed);
    g_fp8x1_lfp_max_limb.store(0, std::memory_order_relaxed);
#endif
}

LFpBoundsStats8x1 fp8x1_lfp_bounds_get()
{
#if defined(NDEBUG)
    return LFpBoundsStats8x1{0, 0, 0};
#else
    return LFpBoundsStats8x1{
        g_fp8x1_lfp_total_lanes.load(std::memory_order_relaxed),
        g_fp8x1_lfp_ge_p_lanes.load(std::memory_order_relaxed),
        g_fp8x1_lfp_max_limb.load(std::memory_order_relaxed),
    };
#endif
}

void FpPack4x2::carryp(FpPack4x2& a)
{
  cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i bmask = VSET1(BMASK);

  a1 = VADD(a1, VSRA(a0, BRADIX)); a0 = VAND(a0, bmask);
  a2 = VADD(a2, VSRA(a1, BRADIX)); a1 = VAND(a1, bmask);
  a3 = VADD(a3, VSRA(a2, BRADIX)); a2 = VAND(a2, bmask);

  a[0] = a0; a[1] = a1; a[2] = a2; a[3] = a3;
}

void FpPack8x1::get_channel(u64* BLST_restrict r, const FpPack8x1& a, const int ch) 
{
  int i;

  for(i = 0; i < NWORDS; i++) {
    r[i] = ((u64 *)&a[i])[ch];
  }
}

void FpPack4x2::get_channel(u64* BLST_restrict r, const FpPack4x2& a, const int ch) 
{
  int i;

  for(i = 0; i < VWORDS; i++) {
    r[i] = ((u64 *)&a[i])[ch];
    r[i+VWORDS] = ((u64 *)&a[i])[ch+1];
  }
}


void FpPack8x1::asx4(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask;
  cc512i t0, t1, t2, t3, t4, t5, t6, t7;

  // r = H+P | G | F+N | E | D+L | C | B+J | A
  r0  = VMADD(a0, 0xAA, a0, b0); r1  = VMADD(a1, 0xAA, a1, b1);
  r2  = VMADD(a2, 0xAA, a2, b2); r3  = VMADD(a3, 0xAA, a3, b3);
  r4  = VMADD(a4, 0xAA, a4, b4); r5  = VMADD(a5, 0xAA, a5, b5);
  r6  = VMADD(a6, 0xAA, a6, b6); r7  = VMADD(a7, 0xAA, a7, b7);

  // t = p | O | p | M | p | K | p | I
  t0 = VMBLEND(0xAA, b0, p0); t1 = VMBLEND(0xAA, b1, p1);
  t2 = VMBLEND(0xAA, b2, p2); t3 = VMBLEND(0xAA, b3, p3);
  t4 = VMBLEND(0xAA, b4, p4); t5 = VMBLEND(0xAA, b5, p5);
  t6 = VMBLEND(0xAA, b6, p6); t7 = VMBLEND(0xAA, b7, p7); 

  // r = H+P-p | G-O | F+N-p | E-M | D+L-p | C-K | B+J-p | A-I
  r0 = VSUB(r0, t0); r1 = VSUB(r1, t1);
  r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
  r4 = VSUB(r4, t4); r5 = VSUB(r5, t5);
  r6 = VSUB(r6, t6); r7 = VSUB(r7, t7);

  // get sign mask 
  t0 = VADD(r1, VSRA(r0, BRADIX));
  t0 = VADD(r2, VSRA(t0, BRADIX));
  t0 = VADD(r3, VSRA(t0, BRADIX));
  t0 = VADD(r4, VSRA(t0, BRADIX));
  t0 = VADD(r5, VSRA(t0, BRADIX));
  t0 = VADD(r6, VSRA(t0, BRADIX));
  t0 = VADD(r7, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)
  r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask)); 
  r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask)); 
  r4 = VADD(r4, VAND(p4, smask)); r5 = VADD(r5, VAND(p5, smask)); 
  r6 = VADD(r6, VAND(p6, smask)); r7 = VADD(r7, VAND(p7, smask)); 

  // carry propagation 
  r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
                                   r7 = VAND(r7, bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
}

void FpPack8x1::asax6(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask;
  cc512i t0, t1, t2, t3, t4, t5, t6, t7;

  // r = H+P | G | F+N | E+M | D+L | C+K | B+J | A+I
  r0  = VMADD(a0, 0xBF, a0, b0); r1  = VMADD(a1, 0xBF, a1, b1);
  r2  = VMADD(a2, 0xBF, a2, b2); r3  = VMADD(a3, 0xBF, a3, b3);
  r4  = VMADD(a4, 0xBF, a4, b4); r5  = VMADD(a5, 0xBF, a5, b5);
  r6  = VMADD(a6, 0xBF, a6, b6); r7  = VMADD(a7, 0xBF, a7, b7);

  // t = p | O | p | p | p | p | p | p
  t0 = VMBLEND(0xBF, b0, p0); t1 = VMBLEND(0xBF, b1, p1);
  t2 = VMBLEND(0xBF, b2, p2); t3 = VMBLEND(0xBF, b3, p3);
  t4 = VMBLEND(0xBF, b4, p4); t5 = VMBLEND(0xBF, b5, p5);
  t6 = VMBLEND(0xBF, b6, p6); t7 = VMBLEND(0xBF, b7, p7); 

  // r = H+P-p | G-O | F+N-p | E+M-p | D+L-p | C+K-p | B+J-p | A+I-p
  r0 = VSUB(r0, t0); r1 = VSUB(r1, t1);
  r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
  r4 = VSUB(r4, t4); r5 = VSUB(r5, t5);
  r6 = VSUB(r6, t6); r7 = VSUB(r7, t7);

  // get sign mask 
  t0 = VADD(r1, VSRA(r0, BRADIX));
  t0 = VADD(r2, VSRA(t0, BRADIX));
  t0 = VADD(r3, VSRA(t0, BRADIX));
  t0 = VADD(r4, VSRA(t0, BRADIX));
  t0 = VADD(r5, VSRA(t0, BRADIX));
  t0 = VADD(r6, VSRA(t0, BRADIX));
  t0 = VADD(r7, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)
  r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask)); 
  r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask)); 
  r4 = VADD(r4, VAND(p4, smask)); r5 = VADD(r5, VAND(p5, smask)); 
  r6 = VADD(r6, VAND(p6, smask)); r7 = VADD(r7, VAND(p7, smask)); 

  // carry propagation 
  r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
                                   r7 = VAND(r7, bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
}

void FpPack8x1::assx6(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask;
  cc512i t0, t1, t2, t3, t4, t5, t6, t7;

  // r = H+P | G | F | E | D | C | B | A
  r0  = VMADD(a0, 0x80, a0, b0); r1  = VMADD(a1, 0x80, a1, b1);
  r2  = VMADD(a2, 0x80, a2, b2); r3  = VMADD(a3, 0x80, a3, b3);
  r4  = VMADD(a4, 0x80, a4, b4); r5  = VMADD(a5, 0x80, a5, b5);
  r6  = VMADD(a6, 0x80, a6, b6); r7  = VMADD(a7, 0x80, a7, b7);

  // t = p | O | N | M | L | K | J | I
  t0 = VMBLEND(0x80, b0, p0); t1 = VMBLEND(0x80, b1, p1);
  t2 = VMBLEND(0x80, b2, p2); t3 = VMBLEND(0x80, b3, p3);
  t4 = VMBLEND(0x80, b4, p4); t5 = VMBLEND(0x80, b5, p5);
  t6 = VMBLEND(0x80, b6, p6); t7 = VMBLEND(0x80, b7, p7); 

  // r = H+P-p | G-O | F-N | E-M | D-L | C-K | B-J | A-I
  r0 = VSUB(r0, t0); r1 = VSUB(r1, t1);
  r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
  r4 = VSUB(r4, t4); r5 = VSUB(r5, t5);
  r6 = VSUB(r6, t6); r7 = VSUB(r7, t7);

  // get sign mask 
  t0 = VADD(r1, VSRA(r0, BRADIX));
  t0 = VADD(r2, VSRA(t0, BRADIX));
  t0 = VADD(r3, VSRA(t0, BRADIX));
  t0 = VADD(r4, VSRA(t0, BRADIX));
  t0 = VADD(r5, VSRA(t0, BRADIX));
  t0 = VADD(r6, VSRA(t0, BRADIX));
  t0 = VADD(r7, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)
  r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask)); 
  r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask)); 
  r4 = VADD(r4, VAND(p4, smask)); r5 = VADD(r5, VAND(p5, smask)); 
  r6 = VADD(r6, VAND(p6, smask)); r7 = VADD(r7, VAND(p7, smask)); 

  // carry propagation 
  r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
                                   r7 = VAND(r7, bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
}

void FpPack4x2::add(FpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask = VSET1(BMASK); 
  cc512i r0, r1, r2, r3, smask, t0;

  // r = a + b
  r0 = VADD(a0, b0); r1 = VADD(a1, b1); r2 = VADD(a2, b2); r3 = VADD(a3, b3);

  // r = r - p
  r0 = VSUB(r0, p0); r1 = VSUB(r1, p1); r2 = VSUB(r2, p2); r3 = VSUB(r3, p3);

  // get sign mask
  t0 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX));
  t0 = VMADD(r2, 0x55, r2, VSRA(t0, BRADIX));
  t0 = VMADD(r3, 0x55, r3, VSRA(t0, BRADIX));
  t0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0 = VMADD(r1, 0xAA, r1, VSRA(t0, BRADIX));
  t0 = VMADD(r2, 0xAA, r2, VSRA(t0, BRADIX));
  t0 = VMADD(r3, 0xAA, r3, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask)); 
  r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask)); 

  // carry propagation 
  r1 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX)); 
  r2 = VMADD(r2, 0x55, r2, VSRA(r1, BRADIX));
  r3 = VMADD(r3, 0x55, r3, VSRA(r2, BRADIX));
  r0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(r3, BRADIX)));
  r1 = VMADD(r1, 0xAA, r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VMADD(r2, 0xAA, r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VMADD(r3, 0xAA, r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r3 = VAND(r3, bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; 
}

void FpPack4x2::sub(FpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask = VSET1(BMASK); 
  cc512i r0, r1, r2, r3, smask, t0;

  // r = a - b
  r0 = VSUB(a0, b0); r1 = VSUB(a1, b1); r2 = VSUB(a2, b2); r3 = VSUB(a3, b3);

  // get sign mask
  t0 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX));
  t0 = VMADD(r2, 0x55, r2, VSRA(t0, BRADIX));
  t0 = VMADD(r3, 0x55, r3, VSRA(t0, BRADIX));
  t0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0 = VMADD(r1, 0xAA, r1, VSRA(t0, BRADIX));
  t0 = VMADD(r2, 0xAA, r2, VSRA(t0, BRADIX));
  t0 = VMADD(r3, 0xAA, r3, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask)); 
  r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask)); 

  // carry propagation 
  r1 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX)); 
  r2 = VMADD(r2, 0x55, r2, VSRA(r1, BRADIX));
  r3 = VMADD(r3, 0x55, r3, VSRA(r2, BRADIX));
  r0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(r3, BRADIX)));
  r1 = VMADD(r1, 0xAA, r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VMADD(r2, 0xAA, r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VMADD(r3, 0xAA, r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r3 = VAND(r3, bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; 
}

void FpPack4x2::asx2(FpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask = VSET1(BMASK); 
  cc512i r0, r1, r2, r3, smask, t0, t1, t2, t3; 

  // r = D'+H' | D+H | C' | C | B'+F' | B+F | A' | A 
  r0 = VMADD(a0, 0xCC, a0, b0); r1 = VMADD(a1, 0xCC, a1, b1);
  r2 = VMADD(a2, 0xCC, a2, b2); r3 = VMADD(a3, 0xCC, a3, b3);

  // t = p' | p | G' | G | p' | p | E' | E
  t0 = VMBLEND(0xCC, b0, p0); t1 = VMBLEND(0xCC, b1, p1);
  t2 = VMBLEND(0xCC, b2, p2); t3 = VMBLEND(0xCC, b3, p3);

  // r = D'+H'-p | D+H-p | C'-G' | C-G | B'+F'-p' | B+F-p | A'-E' | A-E
  r0 = VSUB(r0, t0); r1 = VSUB(r1, t1);
  r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);

  // get sign mask
  t0 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX));
  t0 = VMADD(r2, 0x55, r2, VSRA(t0, BRADIX));
  t0 = VMADD(r3, 0x55, r3, VSRA(t0, BRADIX));
  t0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0 = VMADD(r1, 0xAA, r1, VSRA(t0, BRADIX));
  t0 = VMADD(r2, 0xAA, r2, VSRA(t0, BRADIX));
  t0 = VMADD(r3, 0xAA, r3, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask)); 
  r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask)); 

  // carry propagation 
  r1 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX)); 
  r2 = VMADD(r2, 0x55, r2, VSRA(r1, BRADIX));
  r3 = VMADD(r3, 0x55, r3, VSRA(r2, BRADIX));
  r0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(r3, BRADIX)));
  r1 = VMADD(r1, 0xAA, r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VMADD(r2, 0xAA, r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VMADD(r3, 0xAA, r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r3 = VAND(r3, bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; 
}

void UFpPack8x1::aaasx2(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask;
  cc512i t0, t1, t2 , t3 , t4 , t5 , t6 , t7 ;

  // r = H+P | G+O | F+N | E | D+L | C+K | B+N | A
  r0  = VMADD(a0 , 0xEE, a0 , b0 ); r1  = VMADD(a1 , 0xEE, a1 , b1 );
  r2  = VMADD(a2 , 0xEE, a2 , b2 ); r3  = VMADD(a3 , 0xEE, a3 , b3 );
  r4  = VMADD(a4 , 0xEE, a4 , b4 ); r5  = VMADD(a5 , 0xEE, a5 , b5 );
  r6  = VMADD(a6 , 0xEE, a6 , b6 ); r7  = VMADD(a7 , 0xEE, a7 , b7 );
  r8  = VMADD(a8 , 0xEE, a8 , b8 ); r9  = VMADD(a9 , 0xEE, a9 , b9 );
  r10 = VMADD(a10, 0xEE, a10, b10); r11 = VMADD(a11, 0xEE, a11, b11);
  r12 = VMADD(a12, 0xEE, a12, b12); r13 = VMADD(a13, 0xEE, a13, b13);
  r14 = VMADD(a14, 0xEE, a14, b14); r15 = VMADD(a15, 0xEE, a15, b15);

  // t = p | p | p | M | p | p | p | I
  t0 = VMBLEND(0xEE, b8 , p0); t1 = VMBLEND(0xEE, b9 , p1);
  t2 = VMBLEND(0xEE, b10, p2); t3 = VMBLEND(0xEE, b11, p3);
  t4 = VMBLEND(0xEE, b12, p4); t5 = VMBLEND(0xEE, b13, p5);
  t6 = VMBLEND(0xEE, b14, p6); t7 = VMBLEND(0xEE, b15, p7); 

  // r = H+P | G+O | F+N | E-M | D+L | C+K | B+J | A-I
  r0 = VMSUB(r0, 0x11, r0, b0); r1 = VMSUB(r1, 0x11, r1, b1);
  r2 = VMSUB(r2, 0x11, r2, b2); r3 = VMSUB(r3, 0x11, r3, b3);
  r4 = VMSUB(r4, 0x11, r4, b4); r5 = VMSUB(r5, 0x11, r5, b5);
  r6 = VMSUB(r6, 0x11, r6, b6); r7 = VMSUB(r7, 0x11, r7, b7);
  // r = H+P-p | G+O-p | F+N-p | E-M | D+L-p | C+K-p | B+J-p | A-I
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);
  r12 = VSUB(r12, t4); r13 = VSUB(r13, t5);
  r14 = VSUB(r14, t6); r15 = VSUB(r15, t7);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void UFpPack8x1::aassx2(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask;
  cc512i t0, t1, t2 , t3 , t4 , t5 , t6 , t7 ;

  // r = H+P | G+O | F | E | D+L | C+K | B | A
  r0  = VMADD(a0 , 0xCC, a0 , b0 ); r1  = VMADD(a1 , 0xCC, a1 , b1 );
  r2  = VMADD(a2 , 0xCC, a2 , b2 ); r3  = VMADD(a3 , 0xCC, a3 , b3 );
  r4  = VMADD(a4 , 0xCC, a4 , b4 ); r5  = VMADD(a5 , 0xCC, a5 , b5 );
  r6  = VMADD(a6 , 0xCC, a6 , b6 ); r7  = VMADD(a7 , 0xCC, a7 , b7 );
  r8  = VMADD(a8 , 0xCC, a8 , b8 ); r9  = VMADD(a9 , 0xCC, a9 , b9 );
  r10 = VMADD(a10, 0xCC, a10, b10); r11 = VMADD(a11, 0xCC, a11, b11);
  r12 = VMADD(a12, 0xCC, a12, b12); r13 = VMADD(a13, 0xCC, a13, b13);
  r14 = VMADD(a14, 0xCC, a14, b14); r15 = VMADD(a15, 0xCC, a15, b15);

  // t = p | p | N | M | p | p | J | I
  t0 = VMBLEND(0xCC, b8 , p0); t1 = VMBLEND(0xCC, b9 , p1);
  t2 = VMBLEND(0xCC, b10, p2); t3 = VMBLEND(0xCC, b11, p3);
  t4 = VMBLEND(0xCC, b12, p4); t5 = VMBLEND(0xCC, b13, p5);
  t6 = VMBLEND(0xCC, b14, p6); t7 = VMBLEND(0xCC, b15, p7); 

  // r = H+P | G+O | F-N | E-M | D+L | C+K | B-J | A-I
  r0 = VMSUB(r0, 0x33, r0, b0); r1 = VMSUB(r1, 0x33, r1, b1);
  r2 = VMSUB(r2, 0x33, r2, b2); r3 = VMSUB(r3, 0x33, r3, b3);
  r4 = VMSUB(r4, 0x33, r4, b4); r5 = VMSUB(r5, 0x33, r5, b5);
  r6 = VMSUB(r6, 0x33, r6, b6); r7 = VMSUB(r7, 0x33, r7, b7);
  // r = H+P-p | G+O-p | F-N | E-M | D+L-p | C+K-p | B-J | A-I
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);
  r12 = VSUB(r12, t4); r13 = VSUB(r13, t5);
  r14 = VSUB(r14, t6); r15 = VSUB(r15, t7);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void UFpPack8x1::asx4(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask;
  cc512i t0, t1, t2 , t3 , t4 , t5 , t6 , t7 ;

  // r = H+P | G | F+N | E | D+L | C | B+J | A
  r0  = VMADD(a0 , 0xAA, a0 , b0 ); r1  = VMADD(a1 , 0xAA, a1 , b1 );
  r2  = VMADD(a2 , 0xAA, a2 , b2 ); r3  = VMADD(a3 , 0xAA, a3 , b3 );
  r4  = VMADD(a4 , 0xAA, a4 , b4 ); r5  = VMADD(a5 , 0xAA, a5 , b5 );
  r6  = VMADD(a6 , 0xAA, a6 , b6 ); r7  = VMADD(a7 , 0xAA, a7 , b7 );
  r8  = VMADD(a8 , 0xAA, a8 , b8 ); r9  = VMADD(a9 , 0xAA, a9 , b9 );
  r10 = VMADD(a10, 0xAA, a10, b10); r11 = VMADD(a11, 0xAA, a11, b11);
  r12 = VMADD(a12, 0xAA, a12, b12); r13 = VMADD(a13, 0xAA, a13, b13);
  r14 = VMADD(a14, 0xAA, a14, b14); r15 = VMADD(a15, 0xAA, a15, b15);

  // t = p | O | p | M | p | K | p | I
  t0 = VMBLEND(0xAA, b8 , p0); t1 = VMBLEND(0xAA, b9 , p1);
  t2 = VMBLEND(0xAA, b10, p2); t3 = VMBLEND(0xAA, b11, p3);
  t4 = VMBLEND(0xAA, b12, p4); t5 = VMBLEND(0xAA, b13, p5);
  t6 = VMBLEND(0xAA, b14, p6); t7 = VMBLEND(0xAA, b15, p7); 

  // r = H+P | G-O | F+N | E-M | D+L | C-K | B+J | A-I
  r0 = VMSUB(r0, 0x55, r0, b0); r1 = VMSUB(r1, 0x55, r1, b1);
  r2 = VMSUB(r2, 0x55, r2, b2); r3 = VMSUB(r3, 0x55, r3, b3);
  r4 = VMSUB(r4, 0x55, r4, b4); r5 = VMSUB(r5, 0x55, r5, b5);
  r6 = VMSUB(r6, 0x55, r6, b6); r7 = VMSUB(r7, 0x55, r7, b7);
  // r = H+P-p | G-O | F+N-p | E-M | D+L-p | C-K | B+J-p | A-I
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);
  r12 = VSUB(r12, t4); r13 = VSUB(r13, t5);
  r14 = VSUB(r14, t6); r15 = VSUB(r15, t7);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void UFpPack8x1::ssaax2(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask;
  cc512i t0, t1, t2 , t3 , t4 , t5 , t6 , t7 ;

  // r = H | G | F+N | E+M | D | C | B+J | A+I
  r0  = VMADD(a0 , 0x33, a0 , b0 ); r1  = VMADD(a1 , 0x33, a1 , b1 );
  r2  = VMADD(a2 , 0x33, a2 , b2 ); r3  = VMADD(a3 , 0x33, a3 , b3 );
  r4  = VMADD(a4 , 0x33, a4 , b4 ); r5  = VMADD(a5 , 0x33, a5 , b5 );
  r6  = VMADD(a6 , 0x33, a6 , b6 ); r7  = VMADD(a7 , 0x33, a7 , b7 );
  r8  = VMADD(a8 , 0x33, a8 , b8 ); r9  = VMADD(a9 , 0x33, a9 , b9 );
  r10 = VMADD(a10, 0x33, a10, b10); r11 = VMADD(a11, 0x33, a11, b11);
  r12 = VMADD(a12, 0x33, a12, b12); r13 = VMADD(a13, 0x33, a13, b13);
  r14 = VMADD(a14, 0x33, a14, b14); r15 = VMADD(a15, 0x33, a15, b15);

  // t = P | O | p | p | L | K | p | p
  t0 = VMBLEND(0x33, b8 , p0); t1 = VMBLEND(0x33, b9 , p1);
  t2 = VMBLEND(0x33, b10, p2); t3 = VMBLEND(0x33, b11, p3);
  t4 = VMBLEND(0x33, b12, p4); t5 = VMBLEND(0x33, b13, p5);
  t6 = VMBLEND(0x33, b14, p6); t7 = VMBLEND(0x33, b15, p7); 

  // r = H-P | G-O | F+N | E+M | D-L | C-K | B+J | A+I
  r0 = VMSUB(r0, 0xCC, r0, b0); r1 = VMSUB(r1, 0xCC, r1, b1);
  r2 = VMSUB(r2, 0xCC, r2, b2); r3 = VMSUB(r3, 0xCC, r3, b3);
  r4 = VMSUB(r4, 0xCC, r4, b4); r5 = VMSUB(r5, 0xCC, r5, b5);
  r6 = VMSUB(r6, 0xCC, r6, b6); r7 = VMSUB(r7, 0xCC, r7, b7);
  // r = H-P | G-O | F+N-p | E+M-p | D-L | C-K | B+J-p | A+I-p
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);
  r12 = VSUB(r12, t4); r13 = VSUB(r13, t5);
  r14 = VSUB(r14, t6); r15 = VSUB(r15, t7);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void UFpPack8x1::ssasx2(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask;
  cc512i t0, t1, t2 , t3 , t4 , t5 , t6 , t7 ;

  // r = H | G | F+N | E | D | C | B+N | A
  r0  = VMADD(a0 , 0x22, a0 , b0 ); r1  = VMADD(a1 , 0x22, a1 , b1 );
  r2  = VMADD(a2 , 0x22, a2 , b2 ); r3  = VMADD(a3 , 0x22, a3 , b3 );
  r4  = VMADD(a4 , 0x22, a4 , b4 ); r5  = VMADD(a5 , 0x22, a5 , b5 );
  r6  = VMADD(a6 , 0x22, a6 , b6 ); r7  = VMADD(a7 , 0x22, a7 , b7 );
  r8  = VMADD(a8 , 0x22, a8 , b8 ); r9  = VMADD(a9 , 0x22, a9 , b9 );
  r10 = VMADD(a10, 0x22, a10, b10); r11 = VMADD(a11, 0x22, a11, b11);
  r12 = VMADD(a12, 0x22, a12, b12); r13 = VMADD(a13, 0x22, a13, b13);
  r14 = VMADD(a14, 0x22, a14, b14); r15 = VMADD(a15, 0x22, a15, b15);

  // t = P | O | p | M | L | K | p | I
  t0 = VMBLEND(0x22, b8 , p0); t1 = VMBLEND(0x22, b9 , p1);
  t2 = VMBLEND(0x22, b10, p2); t3 = VMBLEND(0x22, b11, p3);
  t4 = VMBLEND(0x22, b12, p4); t5 = VMBLEND(0x22, b13, p5);
  t6 = VMBLEND(0x22, b14, p6); t7 = VMBLEND(0x22, b15, p7); 

  // r = H-P | G-O | F+N | E-M | D-L | C-K | B+J | A-I
  r0 = VMSUB(r0, 0xDD, r0, b0); r1 = VMSUB(r1, 0xDD, r1, b1);
  r2 = VMSUB(r2, 0xDD, r2, b2); r3 = VMSUB(r3, 0xDD, r3, b3);
  r4 = VMSUB(r4, 0xDD, r4, b4); r5 = VMSUB(r5, 0xDD, r5, b5);
  r6 = VMSUB(r6, 0xDD, r6, b6); r7 = VMSUB(r7, 0xDD, r7, b7);
  // r = H-P | G-O | F+N-p | E-M | D-L | C-K | B+J-p | A-I
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);
  r12 = VSUB(r12, t4); r13 = VSUB(r13, t5);
  r14 = VSUB(r14, t6); r15 = VSUB(r15, t7);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void UFpPack4x2::add(UFpPack4x2& r, const UFpPack4x2& a, const UFpPack4x2& b)
{
  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask0 = VSET(0, BMASK, 0, BMASK, 0, BMASK, 0, BMASK);
  const cc512i bmask1 = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11; 
  cc512i t0, t1, t2, t3, smask;

  // r = a + b
  r0  = VADD(a0 , b0 ); r1  = VADD(a1 , b1 );
  r2  = VADD(a2 , b2 ); r3  = VADD(a3 , b3 );
  r4  = VADD(a4 , b4 ); r5  = VADD(a5 , b5 );
  r6  = VADD(a6 , b6 ); r7  = VADD(a7 , b7 );
  r8  = VADD(a8 , b8 ); r9  = VADD(a9 , b9 );
  r10 = VADD(a10, b10); r11 = VADD(a11, b11);

  // r = r - p*2^384
  r8  = VSUB(r8 , p0); r9  = VSUB(r9 , p1);
  r10 = VSUB(r10, p2); r11 = VSUB(r11, p3);

  // get sign mask + carry propagation
  r1  = VMADD(r1 , 0x55, r1 , VSRA(r0, BRADIX)); 
  r4  = VMADD(r4 , 0x55, r4 , VSHUF<0x4E>(r0)); r0 = VAND(r0, bmask0);
  r2  = VMADD(r2 , 0x55, r2 , VSRA(r1, BRADIX)); 
  r5  = VMADD(r5 , 0x55, r5 , VSHUF<0x4E>(r1)); r1 = VAND(r1, bmask0);
  r3  = VMADD(r3 , 0x55, r3 , VSRA(r2, BRADIX)); 
  r6  = VMADD(r6 , 0x55, r6 , VSHUF<0x4E>(r2)); r2 = VAND(r2, bmask0);
  r4  = VMADD(r4 , 0x55, r4 , VSRA(r3, BRADIX)); 
  r7  = VMADD(r7 , 0x55, r7 , VSHUF<0x4E>(r3)); r3 = VAND(r3, bmask0);
  r5  = VMADD(r5 , 0x55, r5 , VSRA(r4, BRADIX)); 
  r8  = VMADD(r8 , 0x55, r8 , VSHUF<0x4E>(r4)); r4 = VAND(r4, bmask0);
  r6  = VMADD(r6 , 0x55, r6 , VSRA(r5, BRADIX)); 
  r9  = VMADD(r9 , 0x55, r9 , VSHUF<0x4E>(r5)); r5 = VAND(r5, bmask0);
  r7  = VMADD(r7 , 0x55, r7 , VSRA(r6, BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSHUF<0x4E>(r6)); r6 = VAND(r6, bmask0);
  r8  = VMADD(r8 , 0x55, r8 , VSRA(r7, BRADIX)); 
  r11 = VMADD(r11, 0x55, r11, VSHUF<0x4E>(r7)); r7 = VAND(r7, bmask0);
  t0  = VMADD(r9 , 0x55, r9 , VSRA(r8, BRADIX));
  t0  = VMADD(r10, 0x55, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0x55, r11, VSRA(t0, BRADIX));
  t0  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0  = VMADD(r9 , 0xAA, r9 , VSRA(t0, BRADIX));
  t0  = VMADD(r10, 0xAA, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0xAA, r11, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask)); 

  // carry propagation 
  r9  = VMADD(r9 , 0x55, r9 , VSRA(r8 , BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSRA(r9 , BRADIX));
  r11 = VMADD(r11, 0x55, r11, VSRA(r10, BRADIX));
  r8  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(r11, BRADIX)));
  r9  = VMADD(r9 , 0xAA, r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask1);
  r10 = VMADD(r10, 0xAA, r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask1);
  r11 = VMADD(r11, 0xAA, r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask1);
  r11 = VAND(r11, bmask1);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
}

void UFpPack4x2::sub(UFpPack4x2& r, const UFpPack4x2& a, const UFpPack4x2& b)
{
  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask0 = VSET(0, BMASK, 0, BMASK, 0, BMASK, 0, BMASK);
  const cc512i bmask1 = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11; 
  cc512i t0, t1, t2, t3, smask;

  // r = a - b
  r0  = VSUB(a0 , b0 ); r1  = VSUB(a1 , b1 );
  r2  = VSUB(a2 , b2 ); r3  = VSUB(a3 , b3 );
  r4  = VSUB(a4 , b4 ); r5  = VSUB(a5 , b5 );
  r6  = VSUB(a6 , b6 ); r7  = VSUB(a7 , b7 );
  r8  = VSUB(a8 , b8 ); r9  = VSUB(a9 , b9 );
  r10 = VSUB(a10, b10); r11 = VSUB(a11, b11);

  // get sign mask + carry propagation
  r1  = VMADD(r1 , 0x55, r1 , VSRA(r0, BRADIX)); 
  r4  = VMADD(r4 , 0x55, r4 , VSHUF<0x4E>(r0)); r0 = VAND(r0, bmask0);
  r2  = VMADD(r2 , 0x55, r2 , VSRA(r1, BRADIX)); 
  r5  = VMADD(r5 , 0x55, r5 , VSHUF<0x4E>(r1)); r1 = VAND(r1, bmask0);
  r3  = VMADD(r3 , 0x55, r3 , VSRA(r2, BRADIX)); 
  r6  = VMADD(r6 , 0x55, r6 , VSHUF<0x4E>(r2)); r2 = VAND(r2, bmask0);
  r4  = VMADD(r4 , 0x55, r4 , VSRA(r3, BRADIX)); 
  r7  = VMADD(r7 , 0x55, r7 , VSHUF<0x4E>(r3)); r3 = VAND(r3, bmask0);
  r5  = VMADD(r5 , 0x55, r5 , VSRA(r4, BRADIX)); 
  r8  = VMADD(r8 , 0x55, r8 , VSHUF<0x4E>(r4)); r4 = VAND(r4, bmask0);
  r6  = VMADD(r6 , 0x55, r6 , VSRA(r5, BRADIX)); 
  r9  = VMADD(r9 , 0x55, r9 , VSHUF<0x4E>(r5)); r5 = VAND(r5, bmask0);
  r7  = VMADD(r7 , 0x55, r7 , VSRA(r6, BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSHUF<0x4E>(r6)); r6 = VAND(r6, bmask0);
  r8  = VMADD(r8 , 0x55, r8 , VSRA(r7, BRADIX)); 
  r11 = VMADD(r11, 0x55, r11, VSHUF<0x4E>(r7)); r7 = VAND(r7, bmask0);
  t0  = VMADD(r9 , 0x55, r9 , VSRA(r8, BRADIX));
  t0  = VMADD(r10, 0x55, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0x55, r11, VSRA(t0, BRADIX));
  t0  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0  = VMADD(r9 , 0xAA, r9 , VSRA(t0, BRADIX));
  t0  = VMADD(r10, 0xAA, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0xAA, r11, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask)); 

  // carry propagation 
  r9  = VMADD(r9 , 0x55, r9 , VSRA(r8 , BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSRA(r9 , BRADIX));
  r11 = VMADD(r11, 0x55, r11, VSRA(r10, BRADIX));
  r8  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(r11, BRADIX)));
  r9  = VMADD(r9 , 0xAA, r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask1);
  r10 = VMADD(r10, 0xAA, r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask1);
  r11 = VMADD(r11, 0xAA, r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask1);
  r11 = VAND(r11, bmask1);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
}

void UFpPack4x2::asx2(UFpPack4x2& r, const UFpPack4x2& a, const UFpPack4x2& b)
{
  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask0 = VSET(0, BMASK, 0, BMASK, 0, BMASK, 0, BMASK);
  const cc512i bmask1 = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11; 
  cc512i t0, t1, t2, t3, smask;

  // r = D'+H' | D+H | C' | C | B'+F' | B+F | A' | A 
  r0  = VMADD(a0 , 0xCC, a0 , b0 ); r1  = VMADD(a1 , 0xCC, a1 , b1 );
  r2  = VMADD(a2 , 0xCC, a2 , b2 ); r3  = VMADD(a3 , 0xCC, a3 , b3 );
  r4  = VMADD(a4 , 0xCC, a4 , b4 ); r5  = VMADD(a5 , 0xCC, a5 , b5 );
  r6  = VMADD(a6 , 0xCC, a6 , b6 ); r7  = VMADD(a7 , 0xCC, a7 , b7 );
  r8  = VMADD(a8 , 0xCC, a8 , b8 ); r9  = VMADD(a9 , 0xCC, a9 , b9 );
  r10 = VMADD(a10, 0xCC, a10, b10); r11 = VMADD(a11, 0xCC, a11, b11);

  // t = p' | p | G' | G | p' | p | E' | E
  t0 = VMBLEND(0xCC, b8 , p0); t1 = VMBLEND(0xCC, b9 , p1);
  t2 = VMBLEND(0xCC, b10, p2); t3 = VMBLEND(0xCC, b11, p3);

  // r = D'+H' | D+H | C'-G' | C-G | B'+F' | B+F | A'-E' | A-E
  r0 = VMSUB(r0, 0x33, r0, b0); r1 = VMSUB(r1, 0x33, r1, b1);
  r2 = VMSUB(r2, 0x33, r2, b2); r3 = VMSUB(r3, 0x33, r3, b3);
  r4 = VMSUB(r4, 0x33, r4, b4); r5 = VMSUB(r5, 0x33, r5, b5);
  r6 = VMSUB(r6, 0x33, r6, b6); r7 = VMSUB(r7, 0x33, r7, b7);
  // r = D'+H'-p | D+H-p | C'-G' | C-G | B'+F'-p' | B+F-p | A'-E' | A-E
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);

  // get sign mask + carry propagation
  r1  = VMADD(r1 , 0x55, r1 , VSRA(r0, BRADIX)); 
  r4  = VMADD(r4 , 0x55, r4 , VSHUF<0x4E>(r0)); r0 = VAND(r0, bmask0);
  r2  = VMADD(r2 , 0x55, r2 , VSRA(r1, BRADIX)); 
  r5  = VMADD(r5 , 0x55, r5 , VSHUF<0x4E>(r1)); r1 = VAND(r1, bmask0);
  r3  = VMADD(r3 , 0x55, r3 , VSRA(r2, BRADIX)); 
  r6  = VMADD(r6 , 0x55, r6 , VSHUF<0x4E>(r2)); r2 = VAND(r2, bmask0);
  r4  = VMADD(r4 , 0x55, r4 , VSRA(r3, BRADIX)); 
  r7  = VMADD(r7 , 0x55, r7 , VSHUF<0x4E>(r3)); r3 = VAND(r3, bmask0);
  r5  = VMADD(r5 , 0x55, r5 , VSRA(r4, BRADIX)); 
  r8  = VMADD(r8 , 0x55, r8 , VSHUF<0x4E>(r4)); r4 = VAND(r4, bmask0);
  r6  = VMADD(r6 , 0x55, r6 , VSRA(r5, BRADIX)); 
  r9  = VMADD(r9 , 0x55, r9 , VSHUF<0x4E>(r5)); r5 = VAND(r5, bmask0);
  r7  = VMADD(r7 , 0x55, r7 , VSRA(r6, BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSHUF<0x4E>(r6)); r6 = VAND(r6, bmask0);
  r8  = VMADD(r8 , 0x55, r8 , VSRA(r7, BRADIX)); 
  r11 = VMADD(r11, 0x55, r11, VSHUF<0x4E>(r7)); r7 = VAND(r7, bmask0);
  t0  = VMADD(r9 , 0x55, r9 , VSRA(r8, BRADIX));
  t0  = VMADD(r10, 0x55, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0x55, r11, VSRA(t0, BRADIX));
  t0  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0  = VMADD(r9 , 0xAA, r9 , VSRA(t0, BRADIX));
  t0  = VMADD(r10, 0xAA, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0xAA, r11, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask)); 

  // carry propagation 
  r9  = VMADD(r9 , 0x55, r9 , VSRA(r8 , BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSRA(r9 , BRADIX));
  r11 = VMADD(r11, 0x55, r11, VSRA(r10, BRADIX));
  r8  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(r11, BRADIX)));
  r9  = VMADD(r9 , 0xAA, r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask1);
  r10 = VMADD(r10, 0xAA, r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask1);
  r11 = VMADD(r11, 0xAA, r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask1);
  r11 = VAND(r11, bmask1);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
}

void UFpPack4x2::aass(UFpPack4x2& r, const UFpPack4x2& a, const UFpPack4x2& b)
{
  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask0 = VSET(0, BMASK, 0, BMASK, 0, BMASK, 0, BMASK);
  const cc512i bmask1 = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11; 
  cc512i t0, t1, t2, t3, smask;

  // r = D'+H' | D+H | C'+G' | C+G' | B' | B | A' | A 
  r0  = VMADD(a0 , 0xF0, a0 , b0 ); r1  = VMADD(a1 , 0xF0, a1 , b1 );
  r2  = VMADD(a2 , 0xF0, a2 , b2 ); r3  = VMADD(a3 , 0xF0, a3 , b3 );
  r4  = VMADD(a4 , 0xF0, a4 , b4 ); r5  = VMADD(a5 , 0xF0, a5 , b5 );
  r6  = VMADD(a6 , 0xF0, a6 , b6 ); r7  = VMADD(a7 , 0xF0, a7 , b7 );
  r8  = VMADD(a8 , 0xF0, a8 , b8 ); r9  = VMADD(a9 , 0xF0, a9 , b9 );
  r10 = VMADD(a10, 0xF0, a10, b10); r11 = VMADD(a11, 0xF0, a11, b11);

  // t = p' | p | p' | p | F' | F | E' | E
  t0 = VMBLEND(0xF0, b8 , p0); t1 = VMBLEND(0xF0, b9 , p1);
  t2 = VMBLEND(0xF0, b10, p2); t3 = VMBLEND(0xF0, b11, p3);

  // r = D'+H' | D+H | C'+G' | C+G | B'-F' | B-F | A'-E' | A-E
  r0 = VMSUB(r0, 0x0F, r0, b0); r1 = VMSUB(r1, 0x0F, r1, b1);
  r2 = VMSUB(r2, 0x0F, r2, b2); r3 = VMSUB(r3, 0x0F, r3, b3);
  r4 = VMSUB(r4, 0x0F, r4, b4); r5 = VMSUB(r5, 0x0F, r5, b5);
  r6 = VMSUB(r6, 0x0F, r6, b6); r7 = VMSUB(r7, 0x0F, r7, b7);
  // r = D'+H'-p | D+H-p | C'+G'-p' | C+G-p | B'-F' | B-F | A'-E' | A-E
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);

  // get sign mask + carry propagation
  r1  = VMADD(r1 , 0x55, r1 , VSRA(r0, BRADIX)); 
  r4  = VMADD(r4 , 0x55, r4 , VSHUF<0x4E>(r0)); r0 = VAND(r0, bmask0);
  r2  = VMADD(r2 , 0x55, r2 , VSRA(r1, BRADIX)); 
  r5  = VMADD(r5 , 0x55, r5 , VSHUF<0x4E>(r1)); r1 = VAND(r1, bmask0);
  r3  = VMADD(r3 , 0x55, r3 , VSRA(r2, BRADIX)); 
  r6  = VMADD(r6 , 0x55, r6 , VSHUF<0x4E>(r2)); r2 = VAND(r2, bmask0);
  r4  = VMADD(r4 , 0x55, r4 , VSRA(r3, BRADIX)); 
  r7  = VMADD(r7 , 0x55, r7 , VSHUF<0x4E>(r3)); r3 = VAND(r3, bmask0);
  r5  = VMADD(r5 , 0x55, r5 , VSRA(r4, BRADIX)); 
  r8  = VMADD(r8 , 0x55, r8 , VSHUF<0x4E>(r4)); r4 = VAND(r4, bmask0);
  r6  = VMADD(r6 , 0x55, r6 , VSRA(r5, BRADIX)); 
  r9  = VMADD(r9 , 0x55, r9 , VSHUF<0x4E>(r5)); r5 = VAND(r5, bmask0);
  r7  = VMADD(r7 , 0x55, r7 , VSRA(r6, BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSHUF<0x4E>(r6)); r6 = VAND(r6, bmask0);
  r8  = VMADD(r8 , 0x55, r8 , VSRA(r7, BRADIX)); 
  r11 = VMADD(r11, 0x55, r11, VSHUF<0x4E>(r7)); r7 = VAND(r7, bmask0);
  t0  = VMADD(r9 , 0x55, r9 , VSRA(r8, BRADIX));
  t0  = VMADD(r10, 0x55, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0x55, r11, VSRA(t0, BRADIX));
  t0  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0  = VMADD(r9 , 0xAA, r9 , VSRA(t0, BRADIX));
  t0  = VMADD(r10, 0xAA, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0xAA, r11, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask)); 

  // carry propagation 
  r9  = VMADD(r9 , 0x55, r9 , VSRA(r8 , BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSRA(r9 , BRADIX));
  r11 = VMADD(r11, 0x55, r11, VSRA(r10, BRADIX));
  r8  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(r11, BRADIX)));
  r9  = VMADD(r9 , 0xAA, r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask1);
  r10 = VMADD(r10, 0xAA, r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask1);
  r11 = VMADD(r11, 0xAA, r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask1);
  r11 = VAND(r11, bmask1);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
}

void FpPack4x2::mul_v1(UFpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);
  BLST_ASSERT_NO_OVERLAP_AVX(r, b);

  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  cc512i z0 = VZERO(), z1 = VZERO(), z2  = VZERO(), z3  = VZERO();
  cc512i z4 = VZERO(), z5 = VZERO(), z6  = VZERO(), z7  = VZERO();
  cc512i z8 = VZERO(), z9 = VZERO(), z10 = VZERO(), z11 = VZERO();
  cc512i y0 = VZERO(), y1 = VZERO(), y2  = VZERO(), y3  = VZERO();
  cc512i y4 = VZERO(), y5 = VZERO(), y6  = VZERO(), y7  = VZERO();
  cc512i y8 = VZERO(), y9 = VZERO(), y10 = VZERO(), tb;

  tb = VSHUF<0x44>(b0);
  z0 = VMACLO(z0, tb, a0); z1 = VMACLO(z1, tb, a1);
  z2 = VMACLO(z2, tb, a2); z3 = VMACLO(z3, tb, a3);
  y0 = VMACHI(y0, tb, a0); y1 = VMACHI(y1, tb, a1);
  y2 = VMACHI(y2, tb, a2); y3 = VMACHI(y3, tb, a3);

  tb = VSHUF<0x44>(b1);
  z1 = VMACLO(z1, tb, a0); z2 = VMACLO(z2, tb, a1);
  z3 = VMACLO(z3, tb, a2); z4 = VMACLO(z4, tb, a3);
  y1 = VMACHI(y1, tb, a0); y2 = VMACHI(y2, tb, a1);
  y3 = VMACHI(y3, tb, a2); y4 = VMACHI(y4, tb, a3);

  tb = VSHUF<0x44>(b2);
  z2 = VMACLO(z2, tb, a0); z3 = VMACLO(z3, tb, a1);
  z4 = VMACLO(z4, tb, a2); z5 = VMACLO(z5, tb, a3);
  y2 = VMACHI(y2, tb, a0); y3 = VMACHI(y3, tb, a1);
  y4 = VMACHI(y4, tb, a2); y5 = VMACHI(y5, tb, a3);

  tb = VSHUF<0x44>(b3);
  z3 = VMACLO(z3, tb, a0); z4 = VMACLO(z4, tb, a1);
  z5 = VMACLO(z5, tb, a2); z6 = VMACLO(z6, tb, a3);
  y3 = VMACHI(y3, tb, a0); y4 = VMACHI(y4, tb, a1);
  y5 = VMACHI(y5, tb, a2); y6 = VMACHI(y6, tb, a3);

  tb = VSHUF<0xEE>(b0);
  z4 = VMACLO(z4, tb, a0); z5 = VMACLO(z5, tb, a1);
  z6 = VMACLO(z6, tb, a2); z7 = VMACLO(z7, tb, a3);
  y4 = VMACHI(y4, tb, a0); y5 = VMACHI(y5, tb, a1);
  y6 = VMACHI(y6, tb, a2); y7 = VMACHI(y7, tb, a3);

  tb = VSHUF<0xEE>(b1);
  z5 = VMACLO(z5, tb, a0); z6 = VMACLO(z6, tb, a1);
  z7 = VMACLO(z7, tb, a2); z8 = VMACLO(z8, tb, a3);
  y5 = VMACHI(y5, tb, a0); y6 = VMACHI(y6, tb, a1);
  y7 = VMACHI(y7, tb, a2); y8 = VMACHI(y8, tb, a3);

  tb = VSHUF<0xEE>(b2);
  z6 = VMACLO(z6, tb, a0); z7 = VMACLO(z7, tb, a1);
  z8 = VMACLO(z8, tb, a2); z9 = VMACLO(z9, tb, a3);
  y6 = VMACHI(y6, tb, a0); y7 = VMACHI(y7, tb, a1);
  y8 = VMACHI(y8, tb, a2); y9 = VMACHI(y9, tb, a3);

  tb = VSHUF<0xEE>(b3);
  z7 = VMACLO(z7, tb, a0); z8 = VMACLO(z8, tb, a1);
  z9 = VMACLO(z9, tb, a2); z10 = VMACLO(z10, tb, a3);
  y7 = VMACHI(y7, tb, a0); y8 = VMACHI(y8, tb, a1);
  y9 = VMACHI(y9, tb, a2); y10 = VMACHI(y10, tb, a3);

  z1  = VADD(z1 , VSHL(y0 , BALIGN));
  z2  = VADD(z2 , VSHL(y1 , BALIGN));
  z3  = VADD(z3 , VSHL(y2 , BALIGN));
  z4  = VADD(z4 , VSHL(y3 , BALIGN));
  z5  = VADD(z5 , VSHL(y4 , BALIGN));
  z6  = VADD(z6 , VSHL(y5 , BALIGN));
  z7  = VADD(z7 , VSHL(y6 , BALIGN));
  z8  = VADD(z8 , VSHL(y7 , BALIGN));
  z9  = VADD(z9 , VSHL(y8 , BALIGN));
  z10 = VADD(z10, VSHL(y9 , BALIGN));
  z11 = VADD(z11, VSHL(y10, BALIGN));

  r[0 ] = z0 ; r[1 ] = z1 ; r[2 ] = z2 ; r[3 ] = z3 ;
  r[4 ] = z4 ; r[5 ] = z5 ; r[6 ] = z6 ; r[7 ] = z7 ;
  r[8 ] = z8 ; r[9 ] = z9 ; r[10] = z10; r[11] = z11;
}

void FpPack4x2::mul(UFpPack4x2& r, const FpPack4x2& a, const FpPack4x2& b)
{
  FpPack4x2::mul_v1(r, a, b);
}

void UFpPack8x1::redc_core(LFpPack8x1& r, const UFpPack8x1& a)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);

  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  cc512i y0  = VZERO(), y1  = VZERO(), y2  = VZERO(), y3  = VZERO();
  cc512i y4  = VZERO(), y5  = VZERO(), y6  = VZERO(), y7  = VZERO();
  cc512i y8  = VZERO(), y9  = VZERO(), y10 = VZERO(), y11 = VZERO();
  cc512i y12 = VZERO(), y13 = VZERO(), y14 = VZERO();
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK), montw = VSET1(MONT_W_R48), zero = VZERO();
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, u;

  u  = VAND(VMACLO(zero, a0, montw), bmask);
  a0 = VMACLO(a0, u, p0); a1 = VMACLO(a1, u, p1); a2 = VMACLO(a2, u, p2);
  a3 = VMACLO(a3, u, p3); a4 = VMACLO(a4, u, p4); a5 = VMACLO(a5, u, p5);
  a6 = VMACLO(a6, u, p6); a7 = VMACLO(a7, u, p7);
  y0 = VMACHI(y0, u, p0); y1 = VMACHI(y1, u, p1); y2 = VMACHI(y2, u, p2);
  y3 = VMACHI(y3, u, p3); y4 = VMACHI(y4, u, p4); y5 = VMACHI(y5, u, p5);
  y6 = VMACHI(y6, u, p6); y7 = VMACHI(y7, u, p7);
  a1 = VADD(VADD(a1, VSRA(a0, BRADIX)), VSHL(y0, BALIGN));

  u  = VAND(VMACLO(zero, a1, montw), bmask);
  a1 = VMACLO(a1, u, p0); a2 = VMACLO(a2, u, p1); a3 = VMACLO(a3, u, p2);
  a4 = VMACLO(a4, u, p3); a5 = VMACLO(a5, u, p4); a6 = VMACLO(a6, u, p5);
  a7 = VMACLO(a7, u, p6); a8 = VMACLO(a8, u, p7);
  y1 = VMACHI(y1, u, p0); y2 = VMACHI(y2, u, p1); y3 = VMACHI(y3, u, p2);
  y4 = VMACHI(y4, u, p3); y5 = VMACHI(y5, u, p4); y6 = VMACHI(y6, u, p5);
  y7 = VMACHI(y7, u, p6); y8 = VMACHI(y8, u, p7);
  a2 = VADD(VADD(a2, VSRA(a1, BRADIX)), VSHL(y1, BALIGN));

  u  = VAND(VMACLO(zero, a2, montw), bmask);
  a2 = VMACLO(a2, u, p0); a3 = VMACLO(a3, u, p1); a4 = VMACLO(a4, u, p2);
  a5 = VMACLO(a5, u, p3); a6 = VMACLO(a6, u, p4); a7 = VMACLO(a7, u, p5);
  a8 = VMACLO(a8, u, p6); a9 = VMACLO(a9, u, p7);
  y2 = VMACHI(y2, u, p0); y3 = VMACHI(y3, u, p1); y4 = VMACHI(y4, u, p2);
  y5 = VMACHI(y5, u, p3); y6 = VMACHI(y6, u, p4); y7 = VMACHI(y7, u, p5);
  y8 = VMACHI(y8, u, p6); y9 = VMACHI(y9, u, p7);
  a3 = VADD(VADD(a3, VSRA(a2, BRADIX)), VSHL(y2, BALIGN));

  u  = VAND(VMACLO(zero, a3, montw), bmask);
  a3 = VMACLO(a3, u, p0); a4  = VMACLO(a4,  u, p1); a5 = VMACLO(a5, u, p2);
  a6 = VMACLO(a6, u, p3); a7  = VMACLO(a7,  u, p4); a8 = VMACLO(a8, u, p5);
  a9 = VMACLO(a9, u, p6); a10 = VMACLO(a10, u, p7);
  y3 = VMACHI(y3, u, p0); y4  = VMACHI(y4 , u, p1); y5 = VMACHI(y5, u, p2);
  y6 = VMACHI(y6, u, p3); y7  = VMACHI(y7 , u, p4); y8 = VMACHI(y8, u, p5);
  y9 = VMACHI(y9, u, p6); y10 = VMACHI(y10, u, p7);
  a4 = VADD(VADD(a4, VSRA(a3, BRADIX)), VSHL(y3, BALIGN));

  u   = VAND(VMACLO(zero, a4, montw), bmask);
  a4  = VMACLO(a4 , u, p0); a5  = VMACLO(a5 , u, p1); a6 = VMACLO(a6, u, p2);
  a7  = VMACLO(a7 , u, p3); a8  = VMACLO(a8 , u, p4); a9 = VMACLO(a9, u, p5);
  a10 = VMACLO(a10, u, p6); a11 = VMACLO(a11, u, p7);
  y4  = VMACHI(y4 , u, p0); y5  = VMACHI(y5 , u, p1); y6 = VMACHI(y6, u, p2);
  y7  = VMACHI(y7 , u, p3); y8  = VMACHI(y8 , u, p4); y9 = VMACHI(y9, u, p5);
  y10 = VMACHI(y10, u, p6); y11 = VMACHI(y11, u, p7);
  a5  = VADD(VADD(a5, VSRA(a4, BRADIX)), VSHL(y4, BALIGN));

  u   = VAND(VMACLO(zero, a5, montw), bmask);
  a5  = VMACLO(a5 , u, p0); a6  = VMACLO(a6 , u, p1); a7  = VMACLO(a7 , u, p2);
  a8  = VMACLO(a8 , u, p3); a9  = VMACLO(a9 , u, p4); a10 = VMACLO(a10, u, p5);
  a11 = VMACLO(a11, u, p6); a12 = VMACLO(a12, u, p7);
  y5  = VMACHI(y5 , u, p0); y6  = VMACHI(y6 , u, p1); y7  = VMACHI(y7 , u, p2);
  y8  = VMACHI(y8 , u, p3); y9  = VMACHI(y9 , u, p4); y10 = VMACHI(y10, u, p5);
  y11 = VMACHI(y11, u, p6); y12 = VMACHI(y12, u, p7);
  a6  = VADD(VADD(a6, VSRA(a5, BRADIX)), VSHL(y5, BALIGN));

  u   = VAND(VMACLO(zero, a6, montw), bmask);
  a6  = VMACLO(a6 , u, p0); a7  = VMACLO(a7 , u, p1); a8  = VMACLO(a8 , u, p2);
  a9  = VMACLO(a9 , u, p3); a10 = VMACLO(a10, u, p4); a11 = VMACLO(a11, u, p5);
  a12 = VMACLO(a12, u, p6); a13 = VMACLO(a13, u, p7);
  y6  = VMACHI(y6 , u, p0); y7  = VMACHI(y7 , u, p1); y8  = VMACHI(y8 , u, p2);
  y9  = VMACHI(y9 , u, p3); y10 = VMACHI(y10, u, p4); y11 = VMACHI(y11, u, p5);
  y12 = VMACHI(y12, u, p6); y13 = VMACHI(y13, u, p7);
  a7  = VADD(VADD(a7, VSRA(a6, BRADIX)), VSHL(y6, BALIGN));

  u   = VAND(VMACLO(zero, a7, montw), bmask);
  a7  = VMACLO(a7 , u, p0); a8  = VMACLO(a8 , u, p1); a9  = VMACLO(a9 , u, p2);
  a10 = VMACLO(a10, u, p3); a11 = VMACLO(a11, u, p4); a12 = VMACLO(a12, u, p5);
  a13 = VMACLO(a13, u, p6); a14 = VMACLO(a14, u, p7);
  y7  = VMACHI(y7 , u, p0); y8  = VMACHI(y8 , u, p1); y9  = VMACHI(y9 , u, p2);
  y10 = VMACHI(y10, u, p3); y11 = VMACHI(y11, u, p4); y12 = VMACHI(y12, u, p5);
  y13 = VMACHI(y13, u, p6); y14 = VMACHI(y14, u, p7);
  a8  = VADD(VADD(a8, VSRA(a7, BRADIX)), VSHL(y7, BALIGN));

  a9  = VADD(a9 , VSHL(y8 , BALIGN));
  a10 = VADD(a10, VSHL(y9 , BALIGN));
  a11 = VADD(a11, VSHL(y10, BALIGN));
  a12 = VADD(a12, VSHL(y11, BALIGN));
  a13 = VADD(a13, VSHL(y12, BALIGN));
  a14 = VADD(a14, VSHL(y13, BALIGN));
  a15 = VADD(a15, VSHL(y14, BALIGN));

  r0 = a8 ; r1 = a9 ; r2 = a10; r3 = a11;
  r4 = a12; r5 = a13; r6 = a14; r7 = a15;

  // carry propagation only; no final subtraction of p.
  r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
                                   r7 = VAND(r7, bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;

  fp8x1_lfp_bounds_record(r);
}

void UFpPack8x1::normalize(FpPack8x1& r, const LFpPack8x1& a)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);

  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0 = a[0], r1 = a[1], r2 = a[2], r3 = a[3];
  cc512i r4 = a[4], r5 = a[5], r6 = a[6], r7 = a[7];
  cc512i smask, t0;

  // r = r - p
  r0 = VSUB(r0, p0); r1 = VSUB(r1, p1); r2 = VSUB(r2, p2); r3 = VSUB(r3, p3);
  r4 = VSUB(r4, p4); r5 = VSUB(r5, p5); r6 = VSUB(r6, p6); r7 = VSUB(r7, p7);

  // get sign mask
  t0 = VADD(r1, VSRA(r0, BRADIX));
  t0 = VADD(r2, VSRA(t0, BRADIX));
  t0 = VADD(r3, VSRA(t0, BRADIX));
  t0 = VADD(r4, VSRA(t0, BRADIX));
  t0 = VADD(r5, VSRA(t0, BRADIX));
  t0 = VADD(r6, VSRA(t0, BRADIX));
  t0 = VADD(r7, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0
  // if r is negative,     smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)
  r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask));
  r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask));
  r4 = VADD(r4, VAND(p4, smask)); r5 = VADD(r5, VAND(p5, smask));
  r6 = VADD(r6, VAND(p6, smask)); r7 = VADD(r7, VAND(p7, smask));

  // carry propagation
  r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
                                   r7 = VAND(r7, bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
}

void UFpPack8x1::redc(FpPack8x1& r, const UFpPack8x1& a)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);

  LFpPack8x1 lazy{};
  UFpPack8x1::redc_core(lazy, a);
  UFpPack8x1::normalize(r, lazy);
}

void UFpPack4x2::redc(FpPack4x2& r, const UFpPack4x2& a)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);

  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i y0  = VZERO(), y1  = VZERO(), y2  = VZERO(), y3  = VZERO();
  cc512i y4  = VZERO(), y5  = VZERO(), y6  = VZERO(), y7  = VZERO();
  cc512i y8  = VZERO(), y9  = VZERO(), y10 = VZERO();
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask = VSET1(BMASK), montw = VSET1(MONT_W_R48), zero = VZERO();
  cc512i r0, r1, r2, r3, smask, u, t0;

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a0, montw), bmask));
  a0 = VMACLO(a0, u, p0); a1 = VMACLO(a1, u, p1);
  a2 = VMACLO(a2, u, p2); a3 = VMACLO(a3, u, p3);
  y0 = VMACHI(y0, u, p0); y1 = VMACHI(y1, u, p1);
  y2 = VMACHI(y2, u, p2); y3 = VMACHI(y3, u, p3);
  a4 = VMADD(a4, 0x55, a4, VSHUF<0x4E>(a0));
  a1 = VMADD(a1, 0x55, a1, VSRA (a0, BRADIX));
  a1 = VADD(a1, VSHL(y0, BALIGN));

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a1, montw), bmask));
  a1 = VMACLO(a1, u, p0); a2 = VMACLO(a2, u, p1);
  a3 = VMACLO(a3, u, p2); a4 = VMACLO(a4, u, p3);
  y1 = VMACHI(y1, u, p0); y2 = VMACHI(y2, u, p1);
  y3 = VMACHI(y3, u, p2); y4 = VMACHI(y4, u, p3);
  a5 = VMADD(a5, 0x55, a5, VSHUF<0x4E>(a1));
  a2 = VMADD(a2, 0x55, a2, VSRA (a1, BRADIX));
  a2 = VADD(a2, VSHL(y1, BALIGN));

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a2, montw), bmask));
  a2 = VMACLO(a2, u, p0); a3 = VMACLO(a3, u, p1);
  a4 = VMACLO(a4, u, p2); a5 = VMACLO(a5, u, p3);
  y2 = VMACHI(y2, u, p0); y3 = VMACHI(y3, u, p1);
  y4 = VMACHI(y4, u, p2); y5 = VMACHI(y5, u, p3);
  a6 = VMADD(a6, 0x55, a6, VSHUF<0x4E>(a2));
  a3 = VMADD(a3, 0x55, a3, VSRA (a2, BRADIX));
  a3 = VADD(a3, VSHL(y2, BALIGN));

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a3, montw), bmask));
  a3 = VMACLO(a3, u, p0); a4 = VMACLO(a4, u, p1);
  a5 = VMACLO(a5, u, p2); a6 = VMACLO(a6, u, p3);
  y3 = VMACHI(y3, u, p0); y4 = VMACHI(y4, u, p1);
  y5 = VMACHI(y5, u, p2); y6 = VMACHI(y6, u, p3);
  a7 = VMADD(a7, 0x55, a7, VSHUF<0x4E>(a3));
  a4 = VMADD(a4, 0x55, a4, VSRA (a3, BRADIX));
  a4 = VADD(a4, VSHL(y3, BALIGN));

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a4, montw), bmask));
  a4 = VMACLO(a4, u, p0); a5 = VMACLO(a5, u, p1);
  a6 = VMACLO(a6, u, p2); a7 = VMACLO(a7, u, p3);
  y4 = VMACHI(y4, u, p0); y5 = VMACHI(y5, u, p1);
  y6 = VMACHI(y6, u, p2); y7 = VMACHI(y7, u, p3);
  a8 = VMADD(a8, 0x55, a8, VSHUF<0x4E>(a4));
  a5 = VMADD(a5, 0x55, a5, VSRA (a4, BRADIX));
  a5 = VADD(a5, VSHL(y4, BALIGN));

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a5, montw), bmask));
  a5 = VMACLO(a5, u, p0); a6 = VMACLO(a6, u, p1);
  a7 = VMACLO(a7, u, p2); a8 = VMACLO(a8, u, p3);
  y5 = VMACHI(y5, u, p0); y6 = VMACHI(y6, u, p1);
  y7 = VMACHI(y7, u, p2); y8 = VMACHI(y8, u, p3);
  a9 = VMADD(a9, 0x55, a9, VSHUF<0x4E>(a5));
  a6 = VMADD(a6, 0x55, a6, VSRA (a5, BRADIX));
  a6 = VADD(a6, VSHL(y5, BALIGN));

  u   = VSHUF<0x44>(VAND(VMACLO(zero, a6, montw), bmask));
  a6  = VMACLO(a6, u, p0); a7  = VMACLO(a7, u, p1);
  a8  = VMACLO(a8, u, p2); a9  = VMACLO(a9, u, p3);
  y6  = VMACHI(y6, u, p0); y7  = VMACHI(y7, u, p1);
  y8  = VMACHI(y8, u, p2); y9  = VMACHI(y9, u, p3);
  a10 = VMADD(a10, 0x55, a10, VSHUF<0x4E>(a6));
  a7  = VMADD(a7 , 0x55, a7 , VSRA (a6, BRADIX));
  a7  = VADD(a7, VSHL(y6, BALIGN));

  u   = VSHUF<0x44>(VAND(VMACLO(zero, a7, montw), bmask));
  a7  = VMACLO(a7, u, p0); a8   = VMACLO(a8 , u, p1);
  a9  = VMACLO(a9, u, p2); a10  = VMACLO(a10, u, p3);
  y7  = VMACHI(y7, u, p0); y8   = VMACHI(y8 , u, p1);
  y9  = VMACHI(y9, u, p2); y10  = VMACHI(y10, u, p3);
  a11 = VMADD(a11, 0x55, a11, VSHUF<0x4E>(a7));
  a8  = VMADD(a8 , 0x55, a8 , VSRA (a7, BRADIX));
  a8  = VADD(a8, VSHL(y7, BALIGN));

  a9  = VADD(a9 , VSHL(y8 , BALIGN));
  a10 = VADD(a10, VSHL(y9 , BALIGN));
  a11 = VADD(a11, VSHL(y10, BALIGN));

  // final subtraction
  r0 = a8; r1 = a9; r2 = a10; r3 = a11;

  // r = r - p
  r0 = VSUB(r0, p0); r1 = VSUB(r1, p1); r2 = VSUB(r2, p2); r3 = VSUB(r3, p3);

  // get sign mask
  t0 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX));
  t0 = VMADD(r2, 0x55, r2, VSRA(t0, BRADIX));
  t0 = VMADD(r3, 0x55, r3, VSRA(t0, BRADIX));
  t0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0 = VMADD(r1, 0xAA, r1, VSRA(t0, BRADIX));
  t0 = VMADD(r2, 0xAA, r2, VSRA(t0, BRADIX));
  t0 = VMADD(r3, 0xAA, r3, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask)); 
  r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask)); 

  // carry propagation 
  r1 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX)); 
  r2 = VMADD(r2, 0x55, r2, VSRA(r1, BRADIX));
  r3 = VMADD(r3, 0x55, r3, VSRA(r2, BRADIX));
  r0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(r3, BRADIX)));
  r1 = VMADD(r1, 0xAA, r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VMADD(r2, 0xAA, r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VMADD(r3, 0xAA, r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r3 = VAND(r3, bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; 
}

namespace
{
BLST_FORCEINLINE void fp8x1_redc_core_pair_hand(
    LFpPack8x1& out0,
    LFpPack8x1& out1,
    const UFpPack8x1& in0,
    const UFpPack8x1& in1)
{
  BLST_ASSERT_NO_OVERLAP_AVX(out0, in0);
  BLST_ASSERT_NO_OVERLAP_AVX(out0, in1);
  BLST_ASSERT_NO_OVERLAP_AVX(out1, in0);
  BLST_ASSERT_NO_OVERLAP_AVX(out1, in1);

  cc512i a00 = in0[0], a01 = in0[1], a02 = in0[2], a03 = in0[3];
  cc512i a04 = in0[4], a05 = in0[5], a06 = in0[6], a07 = in0[7];
  cc512i a08 = in0[8], a09 = in0[9], a10 = in0[10], a11 = in0[11];
  cc512i a12 = in0[12], a13 = in0[13], a14 = in0[14], a15 = in0[15];

  cc512i b00 = in1[0], b01 = in1[1], b02 = in1[2], b03 = in1[3];
  cc512i b04 = in1[4], b05 = in1[5], b06 = in1[6], b07 = in1[7];
  cc512i b08 = in1[8], b09 = in1[9], b10 = in1[10], b11 = in1[11];
  cc512i b12 = in1[12], b13 = in1[13], b14 = in1[14], b15 = in1[15];

  cc512i y0 = VZERO(), y1 = VZERO(), y2 = VZERO(), y3 = VZERO();
  cc512i y4 = VZERO(), y5 = VZERO(), y6 = VZERO(), y7 = VZERO();
  cc512i z0 = VZERO(), z1 = VZERO(), z2 = VZERO(), z3 = VZERO();
  cc512i z4 = VZERO(), z5 = VZERO(), z6 = VZERO(), z7 = VZERO();

  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK), montw = VSET1(MONT_W_R48), zero = VZERO();
  cc512i u0, u1;

  u0 = VAND(VMACLO(zero, a00, montw), bmask);
  u1 = VAND(VMACLO(zero, b00, montw), bmask);
  a00 = VMACLO(a00, u0, p0); b00 = VMACLO(b00, u1, p0);
  a01 = VMACLO(a01, u0, p1); b01 = VMACLO(b01, u1, p1);
  a02 = VMACLO(a02, u0, p2); b02 = VMACLO(b02, u1, p2);
  a03 = VMACLO(a03, u0, p3); b03 = VMACLO(b03, u1, p3);
  a04 = VMACLO(a04, u0, p4); b04 = VMACLO(b04, u1, p4);
  a05 = VMACLO(a05, u0, p5); b05 = VMACLO(b05, u1, p5);
  a06 = VMACLO(a06, u0, p6); b06 = VMACLO(b06, u1, p6);
  a07 = VMACLO(a07, u0, p7); b07 = VMACLO(b07, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a01 = VADD(VADD(a01, VSRA(a00, BRADIX)), VSHL(y0, BALIGN));
  b01 = VADD(VADD(b01, VSRA(b00, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a01, montw), bmask);
  u1 = VAND(VMACLO(zero, b01, montw), bmask);
  a01 = VMACLO(a01, u0, p0); b01 = VMACLO(b01, u1, p0);
  a02 = VMACLO(a02, u0, p1); b02 = VMACLO(b02, u1, p1);
  a03 = VMACLO(a03, u0, p2); b03 = VMACLO(b03, u1, p2);
  a04 = VMACLO(a04, u0, p3); b04 = VMACLO(b04, u1, p3);
  a05 = VMACLO(a05, u0, p4); b05 = VMACLO(b05, u1, p4);
  a06 = VMACLO(a06, u0, p5); b06 = VMACLO(b06, u1, p5);
  a07 = VMACLO(a07, u0, p6); b07 = VMACLO(b07, u1, p6);
  a08 = VMACLO(a08, u0, p7); b08 = VMACLO(b08, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a02 = VADD(VADD(a02, VSRA(a01, BRADIX)), VSHL(y0, BALIGN));
  b02 = VADD(VADD(b02, VSRA(b01, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a02, montw), bmask);
  u1 = VAND(VMACLO(zero, b02, montw), bmask);
  a02 = VMACLO(a02, u0, p0); b02 = VMACLO(b02, u1, p0);
  a03 = VMACLO(a03, u0, p1); b03 = VMACLO(b03, u1, p1);
  a04 = VMACLO(a04, u0, p2); b04 = VMACLO(b04, u1, p2);
  a05 = VMACLO(a05, u0, p3); b05 = VMACLO(b05, u1, p3);
  a06 = VMACLO(a06, u0, p4); b06 = VMACLO(b06, u1, p4);
  a07 = VMACLO(a07, u0, p5); b07 = VMACLO(b07, u1, p5);
  a08 = VMACLO(a08, u0, p6); b08 = VMACLO(b08, u1, p6);
  a09 = VMACLO(a09, u0, p7); b09 = VMACLO(b09, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a03 = VADD(VADD(a03, VSRA(a02, BRADIX)), VSHL(y0, BALIGN));
  b03 = VADD(VADD(b03, VSRA(b02, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a03, montw), bmask);
  u1 = VAND(VMACLO(zero, b03, montw), bmask);
  a03 = VMACLO(a03, u0, p0); b03 = VMACLO(b03, u1, p0);
  a04 = VMACLO(a04, u0, p1); b04 = VMACLO(b04, u1, p1);
  a05 = VMACLO(a05, u0, p2); b05 = VMACLO(b05, u1, p2);
  a06 = VMACLO(a06, u0, p3); b06 = VMACLO(b06, u1, p3);
  a07 = VMACLO(a07, u0, p4); b07 = VMACLO(b07, u1, p4);
  a08 = VMACLO(a08, u0, p5); b08 = VMACLO(b08, u1, p5);
  a09 = VMACLO(a09, u0, p6); b09 = VMACLO(b09, u1, p6);
  a10 = VMACLO(a10, u0, p7); b10 = VMACLO(b10, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a04 = VADD(VADD(a04, VSRA(a03, BRADIX)), VSHL(y0, BALIGN));
  b04 = VADD(VADD(b04, VSRA(b03, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a04, montw), bmask);
  u1 = VAND(VMACLO(zero, b04, montw), bmask);
  a04 = VMACLO(a04, u0, p0); b04 = VMACLO(b04, u1, p0);
  a05 = VMACLO(a05, u0, p1); b05 = VMACLO(b05, u1, p1);
  a06 = VMACLO(a06, u0, p2); b06 = VMACLO(b06, u1, p2);
  a07 = VMACLO(a07, u0, p3); b07 = VMACLO(b07, u1, p3);
  a08 = VMACLO(a08, u0, p4); b08 = VMACLO(b08, u1, p4);
  a09 = VMACLO(a09, u0, p5); b09 = VMACLO(b09, u1, p5);
  a10 = VMACLO(a10, u0, p6); b10 = VMACLO(b10, u1, p6);
  a11 = VMACLO(a11, u0, p7); b11 = VMACLO(b11, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a05 = VADD(VADD(a05, VSRA(a04, BRADIX)), VSHL(y0, BALIGN));
  b05 = VADD(VADD(b05, VSRA(b04, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a05, montw), bmask);
  u1 = VAND(VMACLO(zero, b05, montw), bmask);
  a05 = VMACLO(a05, u0, p0); b05 = VMACLO(b05, u1, p0);
  a06 = VMACLO(a06, u0, p1); b06 = VMACLO(b06, u1, p1);
  a07 = VMACLO(a07, u0, p2); b07 = VMACLO(b07, u1, p2);
  a08 = VMACLO(a08, u0, p3); b08 = VMACLO(b08, u1, p3);
  a09 = VMACLO(a09, u0, p4); b09 = VMACLO(b09, u1, p4);
  a10 = VMACLO(a10, u0, p5); b10 = VMACLO(b10, u1, p5);
  a11 = VMACLO(a11, u0, p6); b11 = VMACLO(b11, u1, p6);
  a12 = VMACLO(a12, u0, p7); b12 = VMACLO(b12, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a06 = VADD(VADD(a06, VSRA(a05, BRADIX)), VSHL(y0, BALIGN));
  b06 = VADD(VADD(b06, VSRA(b05, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a06, montw), bmask);
  u1 = VAND(VMACLO(zero, b06, montw), bmask);
  a06 = VMACLO(a06, u0, p0); b06 = VMACLO(b06, u1, p0);
  a07 = VMACLO(a07, u0, p1); b07 = VMACLO(b07, u1, p1);
  a08 = VMACLO(a08, u0, p2); b08 = VMACLO(b08, u1, p2);
  a09 = VMACLO(a09, u0, p3); b09 = VMACLO(b09, u1, p3);
  a10 = VMACLO(a10, u0, p4); b10 = VMACLO(b10, u1, p4);
  a11 = VMACLO(a11, u0, p5); b11 = VMACLO(b11, u1, p5);
  a12 = VMACLO(a12, u0, p6); b12 = VMACLO(b12, u1, p6);
  a13 = VMACLO(a13, u0, p7); b13 = VMACLO(b13, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a07 = VADD(VADD(a07, VSRA(a06, BRADIX)), VSHL(y0, BALIGN));
  b07 = VADD(VADD(b07, VSRA(b06, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a07, montw), bmask);
  u1 = VAND(VMACLO(zero, b07, montw), bmask);
  a07 = VMACLO(a07, u0, p0); b07 = VMACLO(b07, u1, p0);
  a08 = VMACLO(a08, u0, p1); b08 = VMACLO(b08, u1, p1);
  a09 = VMACLO(a09, u0, p2); b09 = VMACLO(b09, u1, p2);
  a10 = VMACLO(a10, u0, p3); b10 = VMACLO(b10, u1, p3);
  a11 = VMACLO(a11, u0, p4); b11 = VMACLO(b11, u1, p4);
  a12 = VMACLO(a12, u0, p5); b12 = VMACLO(b12, u1, p5);
  a13 = VMACLO(a13, u0, p6); b13 = VMACLO(b13, u1, p6);
  a14 = VMACLO(a14, u0, p7); b14 = VMACLO(b14, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a08 = VADD(VADD(a08, VSRA(a07, BRADIX)), VSHL(y0, BALIGN));
  b08 = VADD(VADD(b08, VSRA(b07, BRADIX)), VSHL(z0, BALIGN));

  a09 = VADD(a09, VSHL(y1, BALIGN));
  a10 = VADD(a10, VSHL(y2, BALIGN));
  a11 = VADD(a11, VSHL(y3, BALIGN));
  a12 = VADD(a12, VSHL(y4, BALIGN));
  a13 = VADD(a13, VSHL(y5, BALIGN));
  a14 = VADD(a14, VSHL(y6, BALIGN));
  a15 = VADD(a15, VSHL(y7, BALIGN));

  b09 = VADD(b09, VSHL(z1, BALIGN));
  b10 = VADD(b10, VSHL(z2, BALIGN));
  b11 = VADD(b11, VSHL(z3, BALIGN));
  b12 = VADD(b12, VSHL(z4, BALIGN));
  b13 = VADD(b13, VSHL(z5, BALIGN));
  b14 = VADD(b14, VSHL(z6, BALIGN));
  b15 = VADD(b15, VSHL(z7, BALIGN));

  cc512i r00 = a08, r01 = a09, r02 = a10, r03 = a11;
  cc512i r04 = a12, r05 = a13, r06 = a14, r07 = a15;
  cc512i r10 = b08, r11 = b09, r12 = b10, r13 = b11;
  cc512i r14 = b12, r15 = b13, r16 = b14, r17 = b15;

  r01 = VADD(r01, VSRA(r00, BRADIX)); r00 = VAND(r00, bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r02 = VADD(r02, VSRA(r01, BRADIX)); r01 = VAND(r01, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r03 = VADD(r03, VSRA(r02, BRADIX)); r02 = VAND(r02, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r04 = VADD(r04, VSRA(r03, BRADIX)); r03 = VAND(r03, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r05 = VADD(r05, VSRA(r04, BRADIX)); r04 = VAND(r04, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
  r06 = VADD(r06, VSRA(r05, BRADIX)); r05 = VAND(r05, bmask);
  r16 = VADD(r16, VSRA(r15, BRADIX)); r15 = VAND(r15, bmask);
  r07 = VADD(r07, VSRA(r06, BRADIX)); r06 = VAND(r06, bmask); r07 = VAND(r07, bmask);
  r17 = VADD(r17, VSRA(r16, BRADIX)); r16 = VAND(r16, bmask); r17 = VAND(r17, bmask);

  out0[0] = r00; out0[1] = r01; out0[2] = r02; out0[3] = r03;
  out0[4] = r04; out0[5] = r05; out0[6] = r06; out0[7] = r07;
  out1[0] = r10; out1[1] = r11; out1[2] = r12; out1[3] = r13;
  out1[4] = r14; out1[5] = r15; out1[6] = r16; out1[7] = r17;

  fp8x1_lfp_bounds_record(out0);
  fp8x1_lfp_bounds_record(out1);
}

BLST_FORCEINLINE void fp8x1_redc_pair(
    FpPack8x1& out0,
    FpPack8x1& out1,
    const UFpPack8x1& in0,
    const UFpPack8x1& in1)
{
  LFpPack8x1 lazy0{};
  LFpPack8x1 lazy1{};
  fp8x1_redc_core_pair_hand(lazy0, lazy1, in0, in1);
  UFpPack8x1::normalize(out0, lazy0);
  UFpPack8x1::normalize(out1, lazy1);
}

BLST_FORCEINLINE void fp8x1_mul_redc_core_k2(
    LFpPack8x1& lazy0,
    LFpPack8x1& lazy1,
    const FpPack8x1& a0,
    const FpPack8x1& b0,
    const FpPack8x1& a1,
    const FpPack8x1& b1)
{
#if BLST_HAVE_FP8X1_MUL_REDC_CORE_K2_ASM
  blst_fp8x1_mul_redc_core_k2(&lazy0, &lazy1, &a0, &b0, &a1, &b1);
#else
  auto fp8x1_mul_redc_core_single_fused = [](
      LFpPack8x1& r,
      const FpPack8x1& a,
      const FpPack8x1& b) {
    BLST_ASSERT_NO_OVERLAP_AVX(r, a);
    BLST_ASSERT_NO_OVERLAP_AVX(r, b);

    const cc512i ia0 = a[0], ia1 = a[1], ia2 = a[2], ia3 = a[3];
    const cc512i ia4 = a[4], ia5 = a[5], ia6 = a[6], ia7 = a[7];
    const cc512i ib0 = b[0], ib1 = b[1], ib2 = b[2], ib3 = b[3];
    const cc512i ib4 = b[4], ib5 = b[5], ib6 = b[6], ib7 = b[7];
    cc512i ta0, ta1, ta2, ta3;
    cc512i tb0, tb1, tb2, tb3;
    cc512i z0  = VZERO(), z1  = VZERO(), z2  = VZERO(), z3  = VZERO();
    cc512i z4  = VZERO(), z5  = VZERO(), z6  = VZERO(), z7  = VZERO();
    cc512i z8  = VZERO(), z9  = VZERO(), z10 = VZERO(), z11 = VZERO();
    cc512i z12 = VZERO(), z13 = VZERO(), z14 = VZERO(), z15 = VZERO();
    cc512i y0  = VZERO(), y1  = VZERO(), y2  = VZERO(), y3  = VZERO();
    cc512i y4  = VZERO(), y5  = VZERO(), y6  = VZERO(), y7  = VZERO();
    cc512i y8  = VZERO(), y9  = VZERO(), y10 = VZERO(), y11 = VZERO();
    cc512i y12 = VZERO(), y13 = VZERO(), y14 = VZERO(), y15 = VZERO();
    cc512i m0  = VZERO(), m1  = VZERO(), m2  = VZERO(), m3  = VZERO();
    cc512i m4  = VZERO(), m5  = VZERO(), m6  = VZERO(), m7  = VZERO();

    z0 = VMACLO(z0, ia0, ib0);
    y0 = VMACHI(y0, ia0, ib0);
    y0 = VSHL(y0, BALIGN);

    z1 = VMACLO(y0, ia0, ib1); z1 = VMACLO(z1, ia1, ib0);
    y1 = VMACHI(y1, ia0, ib1); y1 = VMACHI(y1, ia1, ib0);
    y1 = VSHL(y1, BALIGN);

    z2 = VMACLO(y1, ia0, ib2); z2 = VMACLO(z2, ia1, ib1); z2 = VMACLO(z2, ia2, ib0);
    y2 = VMACHI(y2, ia0, ib2); y2 = VMACHI(y2, ia1, ib1); y2 = VMACHI(y2, ia2, ib0);
    y2 = VSHL(y2, BALIGN);

    z3 = VMACLO(y2, ia0, ib3); z3 = VMACLO(z3, ia1, ib2); z3 = VMACLO(z3, ia2, ib1);
    z3 = VMACLO(z3, ia3, ib0);
    y3 = VMACHI(y3, ia0, ib3); y3 = VMACHI(y3, ia1, ib2); y3 = VMACHI(y3, ia2, ib1);
    y3 = VMACHI(y3, ia3, ib0);
    y3 = VSHL(y3, BALIGN);

    z4 = VMACLO(y3, ia1, ib3); z4 = VMACLO(z4, ia2, ib2); z4 = VMACLO(z4, ia3, ib1);
    y4 = VMACHI(y4, ia1, ib3); y4 = VMACHI(y4, ia2, ib2); y4 = VMACHI(y4, ia3, ib1);
    y4 = VSHL(y4, BALIGN);

    z5 = VMACLO(y4, ia2, ib3); z5 = VMACLO(z5, ia3, ib2);
    y5 = VMACHI(y5, ia2, ib3); y5 = VMACHI(y5, ia3, ib2);
    y5 = VSHL(y5, BALIGN);

    z6 = VMACLO(y5, ia3, ib3);
    y6 = VMACHI(y6, ia3, ib3);
    y6 = VSHL(y6, BALIGN);
    z7 = y6;

    z8 = VMACLO(z8, ia4, ib4);
    y8 = VMACHI(y8, ia4, ib4);
    y8 = VSHL(y8, BALIGN);

    z9 = VMACLO(y8, ia4, ib5); z9 = VMACLO(z9, ia5, ib4);
    y9 = VMACHI(y9, ia4, ib5); y9 = VMACHI(y9, ia5, ib4);
    y9 = VSHL(y9, BALIGN);

    z10 = VMACLO(y9, ia4, ib6); z10 = VMACLO(z10, ia5, ib5); z10 = VMACLO(z10, ia6, ib4);
    y10 = VMACHI(y10, ia4, ib6); y10 = VMACHI(y10, ia5, ib5); y10 = VMACHI(y10, ia6, ib4);
    y10 = VSHL(y10, BALIGN);

    z11 = VMACLO(y10, ia4, ib7); z11 = VMACLO(z11, ia5, ib6);
    z11 = VMACLO(z11, ia6, ib5); z11 = VMACLO(z11, ia7, ib4);
    y11 = VMACHI(y11, ia4, ib7); y11 = VMACHI(y11, ia5, ib6);
    y11 = VMACHI(y11, ia6, ib5); y11 = VMACHI(y11, ia7, ib4);
    y11 = VSHL(y11, BALIGN);

    z12 = VMACLO(y11, ia5, ib7); z12 = VMACLO(z12, ia6, ib6); z12 = VMACLO(z12, ia7, ib5);
    y12 = VMACHI(y12, ia5, ib7); y12 = VMACHI(y12, ia6, ib6); y12 = VMACHI(y12, ia7, ib5);
    y12 = VSHL(y12, BALIGN);

    z13 = VMACLO(y12, ia6, ib7); z13 = VMACLO(z13, ia7, ib6);
    y13 = VMACHI(y13, ia6, ib7); y13 = VMACHI(y13, ia7, ib6);
    y13 = VSHL(y13, BALIGN);

    z14 = VMACLO(y13, ia7, ib7);
    y14 = VMACHI(y14, ia7, ib7);
    y14 = VSHL(y14, BALIGN);
    z15 = y14;

    ta0 = VADD(ia0, ia4); ta1 = VADD(ia1, ia5);
    ta2 = VADD(ia2, ia6); ta3 = VADD(ia3, ia7);
    tb0 = VADD(ib0, ib4); tb1 = VADD(ib1, ib5);
    tb2 = VADD(ib2, ib6); tb3 = VADD(ib3, ib7);

    y0 = y1 = y2 = y3 = y4 = y5 = y6 = y7 = VZERO();

    m0 = VMACLO(m0, ta0, tb0);
    y0 = VMACHI(y0, ta0, tb0);
    y0 = VSHL(y0, BALIGN);

    m1 = VMACLO(y0, ta0, tb1); m1 = VMACLO(m1, ta1, tb0);
    y1 = VMACHI(y1, ta0, tb1); y1 = VMACHI(y1, ta1, tb0);
    y1 = VSHL(y1, BALIGN);

    m2 = VMACLO(y1, ta0, tb2); m2 = VMACLO(m2, ta1, tb1); m2 = VMACLO(m2, ta2, tb0);
    y2 = VMACHI(y2, ta0, tb2); y2 = VMACHI(y2, ta1, tb1); y2 = VMACHI(y2, ta2, tb0);
    y2 = VSHL(y2, BALIGN);

    m3 = VMACLO(y2, ta0, tb3); m3 = VMACLO(m3, ta1, tb2);
    m3 = VMACLO(m3, ta2, tb1); m3 = VMACLO(m3, ta3, tb0);
    y3 = VMACHI(y3, ta0, tb3); y3 = VMACHI(y3, ta1, tb2);
    y3 = VMACHI(y3, ta2, tb1); y3 = VMACHI(y3, ta3, tb0);
    y3 = VSHL(y3, BALIGN);

    m4 = VMACLO(y3, ta1, tb3); m4 = VMACLO(m4, ta2, tb2); m4 = VMACLO(m4, ta3, tb1);
    y4 = VMACHI(y4, ta1, tb3); y4 = VMACHI(y4, ta2, tb2); y4 = VMACHI(y4, ta3, tb1);
    y4 = VSHL(y4, BALIGN);

    m5 = VMACLO(y4, ta2, tb3); m5 = VMACLO(m5, ta3, tb2);
    y5 = VMACHI(y5, ta2, tb3); y5 = VMACHI(y5, ta3, tb2);
    y5 = VSHL(y5, BALIGN);

    m6 = VMACLO(y5, ta3, tb3);
    y6 = VMACHI(y6, ta3, tb3);
    y6 = VSHL(y6, BALIGN);
    m7 = y6;

    m0 = VSUB(m0, VADD(z0, z8));   m1 = VSUB(m1, VADD(z1, z9));
    m2 = VSUB(m2, VADD(z2, z10));  m3 = VSUB(m3, VADD(z3, z11));
    m4 = VSUB(m4, VADD(z4, z12));  m5 = VSUB(m5, VADD(z5, z13));
    m6 = VSUB(m6, VADD(z6, z14));  m7 = VSUB(m7, VADD(z7, z15));

    z4  = VADD(z4, m0);   z5  = VADD(z5, m1);
    z6  = VADD(z6, m2);   z7  = VADD(z7, m3);
    z8  = VADD(z8, m4);   z9  = VADD(z9, m5);
    z10 = VADD(z10, m6);  z11 = VADD(z11, m7);

    cc512i w0 = z0, w1 = z1, w2 = z2, w3 = z3;
    cc512i w4 = z4, w5 = z5, w6 = z6, w7 = z7;
    cc512i w8 = z8, w9 = z9, w10 = z10, w11 = z11;
    cc512i w12 = z12, w13 = z13, w14 = z14, w15 = z15;

    y0 = y1 = y2 = y3 = y4 = y5 = y6 = y7 = VZERO();
    y8 = y9 = y10 = y11 = y12 = y13 = y14 = VZERO();

    const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
    const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
    const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
    const cc512i bmask = VSET1(BMASK), montw = VSET1(MONT_W_R48), zero = VZERO();
    cc512i u;

    u  = VAND(VMACLO(zero, w0, montw), bmask);
    w0 = VMACLO(w0, u, p0); w1 = VMACLO(w1, u, p1); w2 = VMACLO(w2, u, p2);
    w3 = VMACLO(w3, u, p3); w4 = VMACLO(w4, u, p4); w5 = VMACLO(w5, u, p5);
    w6 = VMACLO(w6, u, p6); w7 = VMACLO(w7, u, p7);
    y0 = VMACHI(y0, u, p0); y1 = VMACHI(y1, u, p1); y2 = VMACHI(y2, u, p2);
    y3 = VMACHI(y3, u, p3); y4 = VMACHI(y4, u, p4); y5 = VMACHI(y5, u, p5);
    y6 = VMACHI(y6, u, p6); y7 = VMACHI(y7, u, p7);
    w1 = VADD(VADD(w1, VSRA(w0, BRADIX)), VSHL(y0, BALIGN));

    u  = VAND(VMACLO(zero, w1, montw), bmask);
    w1 = VMACLO(w1, u, p0); w2 = VMACLO(w2, u, p1); w3 = VMACLO(w3, u, p2);
    w4 = VMACLO(w4, u, p3); w5 = VMACLO(w5, u, p4); w6 = VMACLO(w6, u, p5);
    w7 = VMACLO(w7, u, p6); w8 = VMACLO(w8, u, p7);
    y1 = VMACHI(y1, u, p0); y2 = VMACHI(y2, u, p1); y3 = VMACHI(y3, u, p2);
    y4 = VMACHI(y4, u, p3); y5 = VMACHI(y5, u, p4); y6 = VMACHI(y6, u, p5);
    y7 = VMACHI(y7, u, p6); y8 = VMACHI(y8, u, p7);
    w2 = VADD(VADD(w2, VSRA(w1, BRADIX)), VSHL(y1, BALIGN));

    u  = VAND(VMACLO(zero, w2, montw), bmask);
    w2 = VMACLO(w2, u, p0); w3 = VMACLO(w3, u, p1); w4 = VMACLO(w4, u, p2);
    w5 = VMACLO(w5, u, p3); w6 = VMACLO(w6, u, p4); w7 = VMACLO(w7, u, p5);
    w8 = VMACLO(w8, u, p6); w9 = VMACLO(w9, u, p7);
    y2 = VMACHI(y2, u, p0); y3 = VMACHI(y3, u, p1); y4 = VMACHI(y4, u, p2);
    y5 = VMACHI(y5, u, p3); y6 = VMACHI(y6, u, p4); y7 = VMACHI(y7, u, p5);
    y8 = VMACHI(y8, u, p6); y9 = VMACHI(y9, u, p7);
    w3 = VADD(VADD(w3, VSRA(w2, BRADIX)), VSHL(y2, BALIGN));

    u  = VAND(VMACLO(zero, w3, montw), bmask);
    w3 = VMACLO(w3, u, p0); w4 = VMACLO(w4, u, p1); w5 = VMACLO(w5, u, p2);
    w6 = VMACLO(w6, u, p3); w7 = VMACLO(w7, u, p4); w8 = VMACLO(w8, u, p5);
    w9 = VMACLO(w9, u, p6); w10 = VMACLO(w10, u, p7);
    y3 = VMACHI(y3, u, p0); y4 = VMACHI(y4, u, p1); y5 = VMACHI(y5, u, p2);
    y6 = VMACHI(y6, u, p3); y7 = VMACHI(y7, u, p4); y8 = VMACHI(y8, u, p5);
    y9 = VMACHI(y9, u, p6); y10 = VMACHI(y10, u, p7);
    w4 = VADD(VADD(w4, VSRA(w3, BRADIX)), VSHL(y3, BALIGN));

    u   = VAND(VMACLO(zero, w4, montw), bmask);
    w4  = VMACLO(w4, u, p0); w5  = VMACLO(w5, u, p1); w6 = VMACLO(w6, u, p2);
    w7  = VMACLO(w7, u, p3); w8  = VMACLO(w8, u, p4); w9 = VMACLO(w9, u, p5);
    w10 = VMACLO(w10, u, p6); w11 = VMACLO(w11, u, p7);
    y4  = VMACHI(y4, u, p0); y5  = VMACHI(y5, u, p1); y6 = VMACHI(y6, u, p2);
    y7  = VMACHI(y7, u, p3); y8  = VMACHI(y8, u, p4); y9 = VMACHI(y9, u, p5);
    y10 = VMACHI(y10, u, p6); y11 = VMACHI(y11, u, p7);
    w5 = VADD(VADD(w5, VSRA(w4, BRADIX)), VSHL(y4, BALIGN));

    u   = VAND(VMACLO(zero, w5, montw), bmask);
    w5  = VMACLO(w5, u, p0); w6  = VMACLO(w6, u, p1); w7  = VMACLO(w7, u, p2);
    w8  = VMACLO(w8, u, p3); w9  = VMACLO(w9, u, p4); w10 = VMACLO(w10, u, p5);
    w11 = VMACLO(w11, u, p6); w12 = VMACLO(w12, u, p7);
    y5  = VMACHI(y5, u, p0); y6  = VMACHI(y6, u, p1); y7  = VMACHI(y7, u, p2);
    y8  = VMACHI(y8, u, p3); y9  = VMACHI(y9, u, p4); y10 = VMACHI(y10, u, p5);
    y11 = VMACHI(y11, u, p6); y12 = VMACHI(y12, u, p7);
    w6 = VADD(VADD(w6, VSRA(w5, BRADIX)), VSHL(y5, BALIGN));

    u   = VAND(VMACLO(zero, w6, montw), bmask);
    w6  = VMACLO(w6, u, p0); w7  = VMACLO(w7, u, p1); w8  = VMACLO(w8, u, p2);
    w9  = VMACLO(w9, u, p3); w10 = VMACLO(w10, u, p4); w11 = VMACLO(w11, u, p5);
    w12 = VMACLO(w12, u, p6); w13 = VMACLO(w13, u, p7);
    y6  = VMACHI(y6, u, p0); y7  = VMACHI(y7, u, p1); y8  = VMACHI(y8, u, p2);
    y9  = VMACHI(y9, u, p3); y10 = VMACHI(y10, u, p4); y11 = VMACHI(y11, u, p5);
    y12 = VMACHI(y12, u, p6); y13 = VMACHI(y13, u, p7);
    w7 = VADD(VADD(w7, VSRA(w6, BRADIX)), VSHL(y6, BALIGN));

    u   = VAND(VMACLO(zero, w7, montw), bmask);
    w7  = VMACLO(w7, u, p0); w8  = VMACLO(w8, u, p1); w9  = VMACLO(w9, u, p2);
    w10 = VMACLO(w10, u, p3); w11 = VMACLO(w11, u, p4); w12 = VMACLO(w12, u, p5);
    w13 = VMACLO(w13, u, p6); w14 = VMACLO(w14, u, p7);
    y7  = VMACHI(y7, u, p0); y8  = VMACHI(y8, u, p1); y9  = VMACHI(y9, u, p2);
    y10 = VMACHI(y10, u, p3); y11 = VMACHI(y11, u, p4); y12 = VMACHI(y12, u, p5);
    y13 = VMACHI(y13, u, p6); y14 = VMACHI(y14, u, p7);
    w8 = VADD(VADD(w8, VSRA(w7, BRADIX)), VSHL(y7, BALIGN));

    w9  = VADD(w9, VSHL(y8, BALIGN));
    w10 = VADD(w10, VSHL(y9, BALIGN));
    w11 = VADD(w11, VSHL(y10, BALIGN));
    w12 = VADD(w12, VSHL(y11, BALIGN));
    w13 = VADD(w13, VSHL(y12, BALIGN));
    w14 = VADD(w14, VSHL(y13, BALIGN));
    w15 = VADD(w15, VSHL(y14, BALIGN));

    cc512i r0 = w8, r1 = w9, r2 = w10, r3 = w11;
    cc512i r4 = w12, r5 = w13, r6 = w14, r7 = w15;

    r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
    r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
    r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
    r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
    r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
    r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
    r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
                                   r7 = VAND(r7, bmask);

    r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
    r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
    fp8x1_lfp_bounds_record(r);
  };

  fp8x1_mul_redc_core_single_fused(lazy0, a0, b0);
  fp8x1_mul_redc_core_single_fused(lazy1, a1, b1);
#endif
}

BLST_FORCEINLINE void fp8x1_mul_redc_pair(
    FpPack8x1& out0,
    FpPack8x1& out1,
    const FpPack8x1& a0,
    const FpPack8x1& b0,
    const FpPack8x1& a1,
    const FpPack8x1& b1)
{
  LFpPack8x1 lazy0{};
  LFpPack8x1 lazy1{};
  fp8x1_mul_redc_core_k2(lazy0, lazy1, a0, b0, a1, b1);
  UFpPack8x1::normalize(out0, lazy0);
  UFpPack8x1::normalize(out1, lazy1);
}

template <int K>
BLST_FORCEINLINE void fp8x1_mul_block(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b)
{
  static_assert((K == 1) || (K == 2) || (K == 4), "fp8x1_mul_block supports K in {1,2,4}");

  if constexpr (K == 1) {
    FpPack8x1::mul(out[0], a[0], b[0]);
  } else if constexpr (K == 2) {
    FpPack8x1::mul(out[0], a[0], b[0]);
    FpPack8x1::mul(out[1], a[1], b[1]);
  } else {
    FpPack8x1::mul(out[0], a[0], b[0]);
    FpPack8x1::mul(out[1], a[1], b[1]);
    FpPack8x1::mul(out[2], a[2], b[2]);
    FpPack8x1::mul(out[3], a[3], b[3]);
  }
}

template <int K>
BLST_FORCEINLINE void fp8x1_redc_block(
    FpPack8x1* BLST_restrict out,
    const UFpPack8x1* BLST_restrict in)
{
  static_assert((K == 1) || (K == 2) || (K == 4), "fp8x1_redc_block supports K in {1,2,4}");

  if constexpr (K == 1) {
    UFpPack8x1::redc(out[0], in[0]);
  } else if constexpr (K == 2) {
    fp8x1_redc_pair(out[0], out[1], in[0], in[1]);
  } else {
    fp8x1_redc_pair(out[0], out[1], in[0], in[1]);
    fp8x1_redc_pair(out[2], out[3], in[2], in[3]);
  }
}

template <int K>
BLST_FORCEINLINE void fp8x1_mul_redc_block(
    FpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b)
{
  static_assert((K == 1) || (K == 2) || (K == 4), "fp8x1_mul_redc_block supports K in {1,2,4}");
  if constexpr (K == 1) {
    UFpPack8x1 wide{};
    FpPack8x1::mul(wide, a[0], b[0]);
    UFpPack8x1::redc(out[0], wide);
    return;
  }

  if constexpr (K == 2) {
    fp8x1_mul_redc_pair(out[0], out[1], a[0], b[0], a[1], b[1]);
  } else {
    fp8x1_mul_redc_pair(out[0], out[1], a[0], b[0], a[1], b[1]);
    fp8x1_mul_redc_pair(out[2], out[3], a[2], b[2], a[3], b[3]);
  }
}
} // namespace

template <int K>
BLST_AVX_NOINLINE void fp8x1_mul_k(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  static_assert((K == 1) || (K == 2) || (K == 4), "fp8x1_mul_k supports K in {1,2,4}");
  BLST_ASSERT(out != nullptr || n_blocks == 0);
  BLST_ASSERT(a != nullptr || n_blocks == 0);
  BLST_ASSERT(b != nullptr || n_blocks == 0);

  for (std::size_t i = 0; i < n_blocks; ++i) {
    fp8x1_mul_block<K>(out, a, b);
    out += K;
    a += K;
    b += K;
  }
}

template <int K>
BLST_AVX_NOINLINE void fp8x1_redc_k(
    FpPack8x1* BLST_restrict out,
    const UFpPack8x1* BLST_restrict in,
    const std::size_t n_blocks)
{
  static_assert((K == 1) || (K == 2) || (K == 4), "fp8x1_redc_k supports K in {1,2,4}");
  BLST_ASSERT(out != nullptr || n_blocks == 0);
  BLST_ASSERT(in != nullptr || n_blocks == 0);

  for (std::size_t i = 0; i < n_blocks; ++i) {
    fp8x1_redc_block<K>(out, in);
    out += K;
    in += K;
  }
}

template <int K>
BLST_AVX_NOINLINE void fp8x1_mul_redc_k(
    FpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  static_assert((K == 1) || (K == 2) || (K == 4), "fp8x1_mul_redc_k supports K in {1,2,4}");
  BLST_ASSERT(out != nullptr || n_blocks == 0);
  BLST_ASSERT(a != nullptr || n_blocks == 0);
  BLST_ASSERT(b != nullptr || n_blocks == 0);

  for (std::size_t i = 0; i < n_blocks; ++i) {
    fp8x1_mul_redc_block<K>(out, a, b);
    out += K;
    a += K;
    b += K;
  }
}

void fp8x1_mul_dispatch(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    std::size_t n_packs)
{
  BLST_ASSERT(out != nullptr || n_packs == 0);
  BLST_ASSERT(a != nullptr || n_packs == 0);
  BLST_ASSERT(b != nullptr || n_packs == 0);

  const std::size_t n4 = n_packs / 4;
  if (n4 != 0) {
    fp8x1_mul_k<4>(out, a, b, n4);
    out += n4 * 4;
    a += n4 * 4;
    b += n4 * 4;
    n_packs -= n4 * 4;
  }

  const std::size_t n2 = n_packs / 2;
  if (n2 != 0) {
    fp8x1_mul_k<2>(out, a, b, n2);
    out += n2 * 2;
    a += n2 * 2;
    b += n2 * 2;
    n_packs -= n2 * 2;
  }

  if (n_packs != 0) {
    fp8x1_mul_k<1>(out, a, b, n_packs);
  }
}

void fp8x1_redc_dispatch(
    FpPack8x1* BLST_restrict out,
    const UFpPack8x1* BLST_restrict in,
    std::size_t n_packs)
{
  BLST_ASSERT(out != nullptr || n_packs == 0);
  BLST_ASSERT(in != nullptr || n_packs == 0);

  const std::size_t n4 = n_packs / 4;
  if (n4 != 0) {
    fp8x1_redc_k<4>(out, in, n4);
    out += n4 * 4;
    in += n4 * 4;
    n_packs -= n4 * 4;
  }

  const std::size_t n2 = n_packs / 2;
  if (n2 != 0) {
    fp8x1_redc_k<2>(out, in, n2);
    out += n2 * 2;
    in += n2 * 2;
    n_packs -= n2 * 2;
  }

  if (n_packs != 0) {
    fp8x1_redc_k<1>(out, in, n_packs);
  }
}

void fp8x1_mul_redc_dispatch(
    FpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    std::size_t n_packs)
{
  BLST_ASSERT(out != nullptr || n_packs == 0);
  BLST_ASSERT(a != nullptr || n_packs == 0);
  BLST_ASSERT(b != nullptr || n_packs == 0);

  const std::size_t n4 = n_packs / 4;
  if (n4 != 0) {
    fp8x1_mul_redc_k<4>(out, a, b, n4);
    out += n4 * 4;
    a += n4 * 4;
    b += n4 * 4;
    n_packs -= n4 * 4;
  }

  const std::size_t n2 = n_packs / 2;
  if (n2 != 0) {
    fp8x1_mul_redc_k<2>(out, a, b, n2);
    out += n2 * 2;
    a += n2 * 2;
    b += n2 * 2;
    n_packs -= n2 * 2;
  }

  if (n_packs != 0) {
    fp8x1_mul_redc_k<1>(out, a, b, n_packs);
  }
}

template BLST_AVX_NOINLINE void fp8x1_mul_k<1>(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    std::size_t n_blocks);
template BLST_AVX_NOINLINE void fp8x1_mul_k<2>(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    std::size_t n_blocks);
template BLST_AVX_NOINLINE void fp8x1_mul_k<4>(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    std::size_t n_blocks);

template BLST_AVX_NOINLINE void fp8x1_redc_k<1>(
    FpPack8x1* BLST_restrict out,
    const UFpPack8x1* BLST_restrict in,
    std::size_t n_blocks);
template BLST_AVX_NOINLINE void fp8x1_redc_k<2>(
    FpPack8x1* BLST_restrict out,
    const UFpPack8x1* BLST_restrict in,
    std::size_t n_blocks);
template BLST_AVX_NOINLINE void fp8x1_redc_k<4>(
    FpPack8x1* BLST_restrict out,
    const UFpPack8x1* BLST_restrict in,
    std::size_t n_blocks);

template BLST_AVX_NOINLINE void fp8x1_mul_redc_k<1>(
    FpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    std::size_t n_blocks);
template BLST_AVX_NOINLINE void fp8x1_mul_redc_k<2>(
    FpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    std::size_t n_blocks);
template BLST_AVX_NOINLINE void fp8x1_mul_redc_k<4>(
    FpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    std::size_t n_blocks);

} // namespace blst

namespace blst
{
const std::array<u64, avx::NWORDS> P48 = avx::P48;
}

#endif
