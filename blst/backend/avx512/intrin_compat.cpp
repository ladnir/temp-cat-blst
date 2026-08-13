// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs

#include "blst/backend/avx512/intrin_compat.h"

#include <cstring>

#if defined(BLST_AVX512_PORTABLE)

namespace blst
{
#if defined(_MSC_VER)
# define BLST_CC_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define BLST_CC_NOINLINE __attribute__((noinline))
#else
# define BLST_CC_NOINLINE
#endif

BLST_CC_NOINLINE auto _cc512_set1_epi64(u64 x) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = x;
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_add_epi64(cc512i a, cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] + b.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_sub_epi64(cc512i a, cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] - b.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_and_si512(cc512i a, cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] & b.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_slli_epi64(cc512i a, int bits) -> cc512i
{
    cc512i out{};
    if (static_cast<unsigned>(bits) >= 64) {
        return out;
    }
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] << bits;
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_srli_epi64(cc512i a, int bits) -> cc512i
{
    cc512i out{};
    if (static_cast<unsigned>(bits) >= 64) {
        return out;
    }
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] >> bits;
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_loadu_si512(const void* p) -> cc512i
{
    cc512i out{};
    std::memcpy(out.lane, p, sizeof(out.lane));
    return out;
}

BLST_CC_NOINLINE void _cc512_storeu_si512(void* p, cc512i a)
{
    std::memcpy(p, a.lane, sizeof(a.lane));
}

BLST_CC_NOINLINE auto _cc512_permutexvar_epi64(cc512i idx, cc512i a) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[idx.lane[i] & 7];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_mask_blend_epi64(u64 mask, cc512i a, cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = ((mask >> i) & 1) != 0 ? b.lane[i] : a.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_madd52lo_epu64(cc512i acc, cc512i a, cc512i b) -> cc512i
{
    cc512i out{};
    constexpr u64 kMask26 = (u64{1} << 26) - 1;
    constexpr u64 kMask52 = (u64{1} << 52) - 1;
    for (u64 i = 0; i < 8; ++i) {
        const u64 ax = a.lane[i] & kMask52;
        const u64 bx = b.lane[i] & kMask52;
        const u64 a0 = ax & kMask26;
        const u64 a1 = ax >> 26;
        const u64 b0 = bx & kMask26;
        const u64 b1 = bx >> 26;
        const u64 p0 = a0 * b0;
        const u64 p1 = a0 * b1 + a1 * b0;
        const u64 coeff0 = p0 & kMask26;
        const u64 coeff1 = ((p0 >> 26) + (p1 & kMask26)) & kMask26;
        const u64 low52 = coeff0 | (coeff1 << 26);
        out.lane[i] = acc.lane[i] + low52;
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_madd52hi_epu64(cc512i acc, cc512i a, cc512i b) -> cc512i
{
    cc512i out{};
    constexpr u64 kMask26 = (u64{1} << 26) - 1;
    constexpr u64 kMask52 = (u64{1} << 52) - 1;
    for (u64 i = 0; i < 8; ++i) {
        const u64 ax = a.lane[i] & kMask52;
        const u64 bx = b.lane[i] & kMask52;
        const u64 a0 = ax & kMask26;
        const u64 a1 = ax >> 26;
        const u64 b0 = bx & kMask26;
        const u64 b1 = bx >> 26;
        const u64 p0 = a0 * b0;
        const u64 p1 = a0 * b1 + a1 * b0;
        const u64 p2 = a1 * b1;
        const u64 coeff1 = (p0 >> 26) + (p1 & kMask26);
        const u64 carry1 = coeff1 >> 26;
        u64 coeff2 = (p1 >> 26) + (p2 & kMask26) + carry1;
        const u64 carry2 = coeff2 >> 26;
        coeff2 &= kMask26;
        const u64 coeff3 = (p2 >> 26) + carry2;
        const u64 high52 = coeff2 | ((coeff3 & kMask26) << 26);
        out.lane[i] = acc.lane[i] + high52;
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_min_epu64(cc512i a, cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] < b.lane[i] ? a.lane[i] : b.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_cmpeq_epi64_mask(cc512i a, cc512i b) -> u64
{
    u64 mask = 0;
    for (u64 i = 0; i < 8; ++i) {
        if (a.lane[i] == b.lane[i]) {
            mask |= u64{1} << i;
        }
    }
    return mask;
}

BLST_CC_NOINLINE auto _cc512_extracti64x4_epi64_runtime(int imm8, cc512i a) -> cc256i
{
    cc256i out{};
    const u64 base = (static_cast<u64>(imm8) & 1) * 4;
    for (u64 i = 0; i < 4; ++i) {
        out.lane[i] = a.lane[base + i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc256_extract_epi64_runtime(int imm8, cc256i a) -> u64
{
    return a.lane[static_cast<u64>(imm8) & 3];
}

BLST_CC_NOINLINE auto _cc512_castsi512_si128(cc512i a) -> cc128i
{
    return cc128i{{a.lane[0], a.lane[1]}};
}

BLST_CC_NOINLINE auto _cc128_add_epi64(cc128i a, cc128i b) -> cc128i
{
    return cc128i{{a.lane[0] + b.lane[0], a.lane[1] + b.lane[1]}};
}

BLST_CC_NOINLINE auto _cc128_srli_epi64(cc128i a, int bits) -> cc128i
{
    if (static_cast<unsigned>(bits) >= 64) {
        return cc128i{};
    }
    return cc128i{{a.lane[0] >> bits, a.lane[1] >> bits}};
}

BLST_CC_NOINLINE auto _cc512_broadcastq_epi64(cc128i a) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[0];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc128_extract_epi64_runtime(int imm8, cc128i a) -> u64
{
    return a.lane[static_cast<u64>(imm8) & 1];
}

#undef BLST_CC_NOINLINE
} // namespace blst

#endif
