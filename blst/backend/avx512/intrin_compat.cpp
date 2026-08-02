#include "blst/backend/avx512/intrin_compat.h"

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

BLST_CC_NOINLINE auto _cc512_setzero_si512() -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = 0;
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_set1_epi64(long long x) -> cc512i
{
    cc512i out{};
    const u64 v = static_cast<u64>(x);
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = v;
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_add_epi64(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] + b.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_sub_epi64(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] - b.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_and_si512(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] & b.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_or_si512(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] | b.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_xor_si512(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] ^ b.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_slli_epi64(const cc512i a, int bits) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] << bits;
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_srli_epi64(const cc512i a, int bits) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] >> bits;
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_srai_epi64(const cc512i a, int bits) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const i64 v = static_cast<i64>(a.lane[i]);
        out.lane[i] = static_cast<u64>(v >> bits);
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_sllv_epi64(const cc512i a, const cc512i bits) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const u64 sh = bits.lane[i];
        out.lane[i] = sh < 64 ? a.lane[i] << sh : 0;
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_load_si512(const void* p) -> cc512i
{
    cc512i out{};
    const u64* src = static_cast<const u64*>(p);
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = src[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_loadu_si512(const void* p) -> cc512i
{
    cc512i out{};
    const u64* src = static_cast<const u64*>(p);
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = src[i];
    }
    return out;
}

BLST_CC_NOINLINE void _cc512_store_si512(void* p, const cc512i a)
{
    u64* dst = static_cast<u64*>(p);
    for (u64 i = 0; i < 8; ++i) {
        dst[i] = a.lane[i];
    }
}

BLST_CC_NOINLINE void _cc512_storeu_si512(void* p, const cc512i a)
{
    u64* dst = static_cast<u64*>(p);
    for (u64 i = 0; i < 8; ++i) {
        dst[i] = a.lane[i];
    }
}

BLST_CC_NOINLINE auto _cc512_set_epi64(
    long long x7, long long x6, long long x5, long long x4,
    long long x3, long long x2, long long x1, long long x0) -> cc512i
{
    cc512i out{};
    out.lane[0] = static_cast<u64>(x0);
    out.lane[1] = static_cast<u64>(x1);
    out.lane[2] = static_cast<u64>(x2);
    out.lane[3] = static_cast<u64>(x3);
    out.lane[4] = static_cast<u64>(x4);
    out.lane[5] = static_cast<u64>(x5);
    out.lane[6] = static_cast<u64>(x6);
    out.lane[7] = static_cast<u64>(x7);
    return out;
}

BLST_CC_NOINLINE auto _cc512_permutexvar_epi64(const cc512i idx, const cc512i a) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const u64 src = idx.lane[i] & 7;
        out.lane[i] = a.lane[src];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_mask_blend_epi64(u64 mask, const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const u64 take_b = (mask >> i) & 1;
        out.lane[i] = take_b ? b.lane[i] : a.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_madd52lo_epu64(const cc512i acc, const cc512i a, const cc512i b) -> cc512i
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

BLST_CC_NOINLINE auto _cc512_madd52hi_epu64(const cc512i acc, const cc512i a, const cc512i b) -> cc512i
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

BLST_CC_NOINLINE auto _cc512_mask_add_epi64(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const u64 use_add = (mask >> i) & 1;
        out.lane[i] = use_add ? (a.lane[i] + b.lane[i]) : src.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_mask_sub_epi64(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const u64 use_sub = (mask >> i) & 1;
        out.lane[i] = use_sub ? (a.lane[i] - b.lane[i]) : src.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_min_epu64(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = (a.lane[i] < b.lane[i]) ? a.lane[i] : b.lane[i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_cmpgt_epi64_mask(const cc512i a, const cc512i b) -> u64
{
    u64 mask = 0;
    for (u64 i = 0; i < 8; ++i) {
        if (static_cast<i64>(a.lane[i]) > static_cast<i64>(b.lane[i])) {
            mask |= (u64{1} << i);
        }
    }
    return mask;
}

BLST_CC_NOINLINE auto _cc512_cmpge_epi64_mask(const cc512i a, const cc512i b) -> u64
{
    u64 mask = 0;
    for (u64 i = 0; i < 8; ++i) {
        if (static_cast<i64>(a.lane[i]) >= static_cast<i64>(b.lane[i])) {
            mask |= (u64{1} << i);
        }
    }
    return mask;
}

BLST_CC_NOINLINE auto _cc512_shuffle_epi32_runtime(int imm8, const cc512i a) -> cc512i
{
    cc512i out{};
    u32 in32[16];
    u32 out32[16];
    for (u64 i = 0; i < 8; ++i) {
        in32[i * 2] = static_cast<u32>(a.lane[i] & 0xffffffffu);
        in32[i * 2 + 1] = static_cast<u32>(a.lane[i] >> 32);
    }
    for (u64 block = 0; block < 4; ++block) {
        const u64 base = block * 4;
        for (u64 j = 0; j < 4; ++j) {
            const u64 src = (static_cast<u64>(imm8) >> (2 * j)) & 3u;
            out32[base + j] = in32[base + src];
        }
    }
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = static_cast<u64>(out32[i * 2]) | (static_cast<u64>(out32[i * 2 + 1]) << 32);
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_maskz_shuffle_epi32_runtime(int imm8, u64 mask, const cc512i a) -> cc512i
{
    cc512i out{};
    const cc512i shuf = _cc512_shuffle_epi32_runtime(imm8, a);
    u32 shuf32[16];
    u32 out32[16];
    for (u64 i = 0; i < 8; ++i) {
        shuf32[i * 2] = static_cast<u32>(shuf.lane[i] & 0xffffffffu);
        shuf32[i * 2 + 1] = static_cast<u32>(shuf.lane[i] >> 32);
    }
    for (u64 i = 0; i < 16; ++i) {
        out32[i] = ((mask >> i) & 1u) ? shuf32[i] : 0u;
    }
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = static_cast<u64>(out32[i * 2]) | (static_cast<u64>(out32[i * 2 + 1]) << 32);
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_permutex_epi64_runtime(int imm8, const cc512i a) -> cc512i
{
    cc512i out{};
    for (u64 block = 0; block < 2; ++block) {
        const u64 base = block * 4;
        for (u64 lane = 0; lane < 4; ++lane) {
            const u64 src = (static_cast<u64>(imm8) >> (2 * lane)) & 3u;
            out.lane[base + lane] = a.lane[base + src];
        }
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_maskz_permutex_epi64_runtime(int imm8, u64 mask, const cc512i a) -> cc512i
{
    cc512i out{};
    const cc512i perm = _cc512_permutex_epi64_runtime(imm8, a);
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = ((mask >> i) & 1u) ? perm.lane[i] : 0u;
    }
    return out;
}

BLST_CC_NOINLINE auto _cc512_extracti64x4_epi64_runtime(int imm8, const cc512i a) -> cc256i
{
    cc256i out{};
    const u64 block = static_cast<u64>(imm8) & 1u;
    const u64 base = block * 4;
    for (u64 i = 0; i < 4; ++i) {
        out.lane[i] = a.lane[base + i];
    }
    return out;
}

BLST_CC_NOINLINE auto _cc256_extract_epi64_runtime(int imm8, const cc256i a) -> u64
{
    return a.lane[static_cast<u64>(imm8) & 3u];
}

#undef BLST_CC_NOINLINE
} // namespace blst

#endif
