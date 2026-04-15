#include "blst/backend/avx512/intrin_compat.h"

#if defined(BLST_AVX512_PORTABLE) && !defined(BLST_AVX512)

namespace blst
{
#if defined(_MSC_VER)
# define BLST_AVX_PORTABLE_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define BLST_AVX_PORTABLE_NOINLINE __attribute__((noinline))
#else
# define BLST_AVX_PORTABLE_NOINLINE
#endif

BLST_AVX_PORTABLE_NOINLINE auto setzero_si512_portable() -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = 0;
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto set1_epi64_portable(long long x) -> cc512i
{
    cc512i out{};
    const u64 v = static_cast<u64>(x);
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = v;
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto add_epi64_portable(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] + b.lane[i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto sub_epi64_portable(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] - b.lane[i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto and_si512_portable(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] & b.lane[i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto or_si512_portable(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] | b.lane[i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto xor_si512_portable(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] ^ b.lane[i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto slli_epi64_portable(const cc512i a, int bits) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] << bits;
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto srli_epi64_portable(const cc512i a, int bits) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = a.lane[i] >> bits;
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto srai_epi64_portable(const cc512i a, int bits) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const i64 v = static_cast<i64>(a.lane[i]);
        out.lane[i] = static_cast<u64>(v >> bits);
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto sllv_epi64_portable(const cc512i a, const cc512i bits) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const u64 sh = bits.lane[i] & 63;
        out.lane[i] = a.lane[i] << sh;
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto load_si512_portable(const void* p) -> cc512i
{
    cc512i out{};
    const u64* src = static_cast<const u64*>(p);
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = src[i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto loadu_si512_portable(const void* p) -> cc512i
{
    cc512i out{};
    const u64* src = static_cast<const u64*>(p);
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = src[i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE void store_si512_portable(void* p, const cc512i a)
{
    u64* dst = static_cast<u64*>(p);
    for (u64 i = 0; i < 8; ++i) {
        dst[i] = a.lane[i];
    }
}

BLST_AVX_PORTABLE_NOINLINE void storeu_si512_portable(void* p, const cc512i a)
{
    u64* dst = static_cast<u64*>(p);
    for (u64 i = 0; i < 8; ++i) {
        dst[i] = a.lane[i];
    }
}

BLST_AVX_PORTABLE_NOINLINE auto set_epi64_portable(
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

BLST_AVX_PORTABLE_NOINLINE auto permutexvar_epi64_portable(const cc512i idx, const cc512i a) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const u64 src = idx.lane[i] & 7;
        out.lane[i] = a.lane[src];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto mask_blend_epi64_portable(u64 mask, const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const u64 take_b = (mask >> i) & 1;
        out.lane[i] = take_b ? b.lane[i] : a.lane[i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto madd52lo_epu64_portable(const cc512i acc, const cc512i a, const cc512i b) -> cc512i
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

BLST_AVX_PORTABLE_NOINLINE auto madd52hi_epu64_portable(const cc512i acc, const cc512i a, const cc512i b) -> cc512i
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

        u64 coeff1 = (p0 >> 26) + (p1 & kMask26);
        const u64 carry1 = coeff1 >> 26;
        coeff1 &= kMask26;

        u64 coeff2 = (p1 >> 26) + (p2 & kMask26) + carry1;
        const u64 carry2 = coeff2 >> 26;
        coeff2 &= kMask26;

        const u64 coeff3 = (p2 >> 26) + carry2;
        const u64 high52 = coeff2 | ((coeff3 & kMask26) << 26);

        out.lane[i] = acc.lane[i] + high52;
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto mask_add_epi64_portable(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const u64 use_add = (mask >> i) & 1;
        out.lane[i] = use_add ? (a.lane[i] + b.lane[i]) : src.lane[i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto mask_sub_epi64_portable(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        const u64 use_sub = (mask >> i) & 1;
        out.lane[i] = use_sub ? (a.lane[i] - b.lane[i]) : src.lane[i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto min_epu64_portable(const cc512i a, const cc512i b) -> cc512i
{
    cc512i out{};
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = (a.lane[i] < b.lane[i]) ? a.lane[i] : b.lane[i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto cmpgt_epi64_mask_portable(const cc512i a, const cc512i b) -> u64
{
    u64 mask = 0;
    for (u64 i = 0; i < 8; ++i) {
        if (static_cast<i64>(a.lane[i]) > static_cast<i64>(b.lane[i])) {
            mask |= (u64{1} << i);
        }
    }
    return mask;
}

BLST_AVX_PORTABLE_NOINLINE auto cmpge_epi64_mask_portable(const cc512i a, const cc512i b) -> u64
{
    u64 mask = 0;
    for (u64 i = 0; i < 8; ++i) {
        if (static_cast<i64>(a.lane[i]) >= static_cast<i64>(b.lane[i])) {
            mask |= (u64{1} << i);
        }
    }
    return mask;
}

BLST_AVX_PORTABLE_NOINLINE auto shuffle_epi32_portable_runtime(int imm8, const cc512i a) -> cc512i
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

BLST_AVX_PORTABLE_NOINLINE auto maskz_shuffle_epi32_portable_runtime(int imm8, u64 mask, const cc512i a) -> cc512i
{
    cc512i out{};
    const cc512i shuf = shuffle_epi32_portable_runtime(imm8, a);
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

BLST_AVX_PORTABLE_NOINLINE auto permutex_epi64_portable_runtime(int imm8, const cc512i a) -> cc512i
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

BLST_AVX_PORTABLE_NOINLINE auto maskz_permutex_epi64_portable_runtime(int imm8, u64 mask, const cc512i a) -> cc512i
{
    cc512i out{};
    const cc512i perm = permutex_epi64_portable_runtime(imm8, a);
    for (u64 i = 0; i < 8; ++i) {
        out.lane[i] = ((mask >> i) & 1u) ? perm.lane[i] : 0u;
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto extracti64x4_epi64_portable_runtime(int imm8, const cc512i a) -> cc256i
{
    cc256i out{};
    const u64 block = static_cast<u64>(imm8) & 1u;
    const u64 base = block * 4;
    for (u64 i = 0; i < 4; ++i) {
        out.lane[i] = a.lane[base + i];
    }
    return out;
}

BLST_AVX_PORTABLE_NOINLINE auto extract_epi64_portable_runtime(int imm8, const cc256i a) -> u64
{
    return a.lane[static_cast<u64>(imm8) & 3u];
}

#undef BLST_AVX_PORTABLE_NOINLINE
} // namespace blst

#endif
