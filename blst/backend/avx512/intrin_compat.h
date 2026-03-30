#ifndef BLST_BACKEND_AVX512_INTRIN_COMPAT_H
#define BLST_BACKEND_AVX512_INTRIN_COMPAT_H

#include "blst/common/defines.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#if defined(BLST_AVX512)
#include <immintrin.h>
#endif

#if defined(BLST_AVX512)
#define BLST_AVX_NATIVE 1
#else
#define BLST_AVX_NATIVE 0
#endif

namespace blst
{
    // Logical 512-bit integer register used by tile code.
    struct cc512i {
#if BLST_AVX_NATIVE
        ::__m512i data;
#else
        u64 lane[8];
#endif
    };

    // Logical 256-bit integer register used only for extraction helpers.
    struct cc256i {
#if BLST_AVX_NATIVE
        ::__m256i data;
#else
        u64 lane[4];
#endif
    };

#if BLST_AVX_NATIVE
    BLST_FORCEINLINE auto from_native_512(const ::__m512i in) -> cc512i
    {
        cc512i out{};
        out.data = in;
        return out;
    }

    BLST_FORCEINLINE auto to_native_512(const cc512i in) -> ::__m512i
    {
        return in.data;
    }

    BLST_FORCEINLINE auto from_native_256(const ::__m256i in) -> cc256i
    {
        cc256i out{};
        out.data = in;
        return out;
    }

    BLST_FORCEINLINE auto to_native_256(const cc256i in) -> ::__m256i
    {
        return in.data;
    }
#endif

#if !BLST_AVX_NATIVE
    // Portable implementation helpers
    
    inline auto setzero_si512_portable() -> cc512i
    {
        cc512i out{};
        // Default init already zeros it, but to be explicit or if default changes:
        for (u64 i = 0; i < 8; ++i) out.lane[i] = 0;
        return out;
    }

    inline auto set1_epi64_portable(long long x) -> cc512i
    {
        cc512i out{};
        const u64 v = static_cast<u64>(x);
        for (u64 i = 0; i < 8; ++i) {
            out.lane[i] = v;
        }
        return out;
    }

    inline auto add_epi64_portable(const cc512i a, const cc512i b) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            out.lane[i] = a.lane[i] + b.lane[i];
        }
        return out;
    }

    inline auto sub_epi64_portable(const cc512i a, const cc512i b) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            out.lane[i] = a.lane[i] - b.lane[i];
        }
        return out;
    }

    inline auto and_si512_portable(const cc512i a, const cc512i b) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            out.lane[i] = a.lane[i] & b.lane[i];
        }
        return out;
    }

    inline auto or_si512_portable(const cc512i a, const cc512i b) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            out.lane[i] = a.lane[i] | b.lane[i];
        }
        return out;
    }

    inline auto xor_si512_portable(const cc512i a, const cc512i b) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            out.lane[i] = a.lane[i] ^ b.lane[i];
        }
        return out;
    }

    inline auto slli_epi64_portable(const cc512i a, int bits) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            out.lane[i] = a.lane[i] << bits;
        }
        return out;
    }

    inline auto srli_epi64_portable(const cc512i a, int bits) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            out.lane[i] = a.lane[i] >> bits;
        }
        return out;
    }

    inline auto srai_epi64_portable(const cc512i a, int bits) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            const i64 v = static_cast<i64>(a.lane[i]);
            out.lane[i] = static_cast<u64>(v >> bits);
        }
        return out;
    }

    inline auto sllv_epi64_portable(const cc512i a, const cc512i bits) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            const u64 sh = bits.lane[i] & 63;
            out.lane[i] = a.lane[i] << sh;
        }
        return out;
    }

    inline auto load_si512_portable(const void* p) -> cc512i
    {
        cc512i out{};
        const u64* src = static_cast<const u64*>(p);
        for (u64 i = 0; i < 8; ++i) {
            out.lane[i] = src[i];
        }
        return out;
    }

    inline auto loadu_si512_portable(const void* p) -> cc512i
    {
        cc512i out{};
        const u64* src = static_cast<const u64*>(p);
        for (u64 i = 0; i < 8; ++i) {
            out.lane[i] = src[i];
        }
        return out;
    }

    inline void store_si512_portable(void* p, const cc512i a)
    {
        u64* dst = static_cast<u64*>(p);
        for (u64 i = 0; i < 8; ++i) {
            dst[i] = a.lane[i];
        }
    }

    inline void storeu_si512_portable(void* p, const cc512i a)
    {
        u64* dst = static_cast<u64*>(p);
        for (u64 i = 0; i < 8; ++i) {
            dst[i] = a.lane[i];
        }
    }

    inline auto set_epi64_portable(
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

    inline auto permutexvar_epi64_portable(const cc512i idx, const cc512i a) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            const u64 src = idx.lane[i] & 7;
            out.lane[i] = a.lane[src];
        }
        return out;
    }

    inline auto mask_blend_epi64_portable(u64 mask, const cc512i a, const cc512i b) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            const u64 take_b = (mask >> i) & 1;
            out.lane[i] = take_b ? b.lane[i] : a.lane[i];
        }
        return out;
    }

    inline auto madd52lo_epu64_portable(const cc512i acc, const cc512i a, const cc512i b) -> cc512i
    {
        cc512i out{};
        constexpr u64 kMask26 = ((u64)1 << 26) - 1;
        constexpr u64 kMask52 = ((u64)1 << 52) - 1;
        // Exact emulation of madd52lo:
        // acc + low52( (a & ((1<<52)-1)) * (b & ((1<<52)-1)) )
        for (u64 i = 0; i < 8; ++i) {
            const u64 ax = a.lane[i] & kMask52;
            const u64 bx = b.lane[i] & kMask52;
            const u64 a0 = ax & kMask26;
            const u64 a1 = ax >> 26;
            const u64 b0 = bx & kMask26;
            const u64 b1 = bx >> 26;

            const u64 p0 = a0 * b0;
            const u64 p1 = a0 * b1 + a1 * b0;

            // Base-2^26 composition of the low 52 bits.
            const u64 coeff0 = p0 & kMask26;
            const u64 coeff1 = ((p0 >> 26) + (p1 & kMask26)) & kMask26;
            const u64 low52 = coeff0 | (coeff1 << 26);

            out.lane[i] = acc.lane[i] + low52;
        }
        return out;
    }

    inline auto madd52hi_epu64_portable(const cc512i acc, const cc512i a, const cc512i b) -> cc512i
    {
        cc512i out{};
        constexpr u64 kMask26 = ((u64)1 << 26) - 1;
        constexpr u64 kMask52 = ((u64)1 << 52) - 1;
        // Exact emulation of madd52hi:
        // acc + high52( (a & ((1<<52)-1)) * (b & ((1<<52)-1)) )
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

            // Product in base 2^26:
            // p = coeff0 + coeff1*2^26 + coeff2*2^52 + coeff3*2^78
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

    inline auto mask_add_epi64_portable(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            const u64 use_add = (mask >> i) & 1;
            out.lane[i] = use_add ? (a.lane[i] + b.lane[i]) : src.lane[i];
        }
        return out;
    }

    inline auto mask_sub_epi64_portable(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i
    {
        cc512i out{};
        for (u64 i = 0; i < 8; ++i) {
            const u64 use_sub = (mask >> i) & 1;
            out.lane[i] = use_sub ? (a.lane[i] - b.lane[i]) : src.lane[i];
        }
        return out;
    }

    template <int imm8>
    inline auto shuffle_epi32_portable(const cc512i a) -> cc512i
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

    template <int imm8>
    inline auto maskz_shuffle_epi32_portable(u64 mask, const cc512i a) -> cc512i
    {
        cc512i out{};
        const cc512i shuf = shuffle_epi32_portable<imm8>(a);
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

    template <int imm8>
    inline auto permutex_epi64_portable(const cc512i a) -> cc512i
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

    template <int imm8>
    inline auto maskz_permutex_epi64_portable(u64 mask, const cc512i a) -> cc512i
    {
        cc512i out{};
        const cc512i perm = permutex_epi64_portable<imm8>(a);
        for (u64 i = 0; i < 8; ++i) {
            out.lane[i] = ((mask >> i) & 1u) ? perm.lane[i] : 0u;
        }
        return out;
    }

    template <int imm8>
    inline auto extracti64x4_epi64_portable(const cc512i a) -> cc256i
    {
        cc256i out{};
        const u64 block = static_cast<u64>(imm8) & 1u;
        const u64 base = block * 4;
        for (u64 i = 0; i < 4; ++i) {
            out.lane[i] = a.lane[base + i];
        }
        return out;
    }

    template <int imm8>
    inline auto extract_epi64_portable(const cc256i a) -> u64
    {
        return a.lane[static_cast<u64>(imm8) & 3u];
    }
#endif // !BLST_AVX_NATIVE

    //
    // Primary API: Static Dispatch
    //
    
    // Set all lanes to zero
    BLST_FORCEINLINE auto setzero_si512() -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_setzero_si512());
#else
        return setzero_si512_portable();
#endif
    }

    // Broadcast a 64-bit integer to all lanes
    BLST_FORCEINLINE auto set1_epi64(long long x) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_set1_epi64(x));
#else
        return set1_epi64_portable(x);
#endif
    }

    // Add packed 64-bit integers
    BLST_FORCEINLINE auto add_epi64(const cc512i a, const cc512i b) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_add_epi64(to_native_512(a), to_native_512(b)));
#else
        return add_epi64_portable(a, b);
#endif
    }

    // Subtract packed 64-bit integers
    BLST_FORCEINLINE auto sub_epi64(const cc512i a, const cc512i b) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_sub_epi64(to_native_512(a), to_native_512(b)));
#else
        return sub_epi64_portable(a, b);
#endif
    }

    // Bitwise AND
    BLST_FORCEINLINE auto and_si512(const cc512i a, const cc512i b) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_and_si512(to_native_512(a), to_native_512(b)));
#else
        return and_si512_portable(a, b);
#endif
    }

    // Bitwise OR
    BLST_FORCEINLINE auto or_si512(const cc512i a, const cc512i b) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_or_si512(to_native_512(a), to_native_512(b)));
#else
        return or_si512_portable(a, b);
#endif
    }

    // Shift left logical (immediate)
    BLST_FORCEINLINE auto slli_epi64(const cc512i a, int bits) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_slli_epi64(to_native_512(a), bits));
#else
        return slli_epi64_portable(a, bits);
#endif
    }

    // Shift right logical (immediate)
    BLST_FORCEINLINE auto srli_epi64(const cc512i a, int bits) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_srli_epi64(to_native_512(a), bits));
#else
        return srli_epi64_portable(a, bits);
#endif
    }

    // Shift right arithmetic (immediate)
    BLST_FORCEINLINE auto srai_epi64(const cc512i a, int bits) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_srai_epi64(to_native_512(a), bits));
#else
        return srai_epi64_portable(a, bits);
#endif
    }

    // Load 512-bit integer (aligned)
    BLST_FORCEINLINE auto load_si512(const void* p) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_load_si512(p));
#else
        return load_si512_portable(p);
#endif
    }

    // Load 512-bit integer (unaligned)
    BLST_FORCEINLINE auto loadu_si512(const void* p) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_loadu_si512(p));
#else
        return loadu_si512_portable(p);
#endif
    }

    // Store 512-bit integer (aligned)
    BLST_FORCEINLINE void store_si512(void* p, const cc512i a)
    {
#if BLST_AVX_NATIVE
        _mm512_store_si512(p, to_native_512(a));
#else
        store_si512_portable(p, a);
#endif
    }

    // Store 512-bit integer (unaligned)
    BLST_FORCEINLINE void storeu_si512(void* p, const cc512i a)
    {
#if BLST_AVX_NATIVE
        _mm512_storeu_si512(p, to_native_512(a));
#else
        storeu_si512_portable(p, a);
#endif
    }

    // Multiply and Add Low 52-bits
    BLST_FORCEINLINE auto madd52lo_epu64(const cc512i acc, const cc512i a, const cc512i b) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_madd52lo_epu64(to_native_512(acc), to_native_512(a), to_native_512(b)));
#else
        return madd52lo_epu64_portable(acc, a, b);
#endif
    }

    // Multiply and Add High 52-bits
    BLST_FORCEINLINE auto madd52hi_epu64(const cc512i acc, const cc512i a, const cc512i b) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_madd52hi_epu64(to_native_512(acc), to_native_512(a), to_native_512(b)));
#else
        return madd52hi_epu64_portable(acc, a, b);
#endif
    }
    
    // --- Legacy _cc512_* wrappers mapping to standardized names or specific portable implementations ---

    BLST_FORCEINLINE auto _cc512_setzero_si512() -> cc512i { return setzero_si512(); }
    BLST_FORCEINLINE auto _cc512_set1_epi64(long long x) -> cc512i { return set1_epi64(x); }
    BLST_FORCEINLINE auto _cc512_add_epi64(const cc512i a, const cc512i b) -> cc512i { return add_epi64(a, b); }

    // Masked addition: src if mask bit 0, a+b if mask bit 1
    BLST_FORCEINLINE auto _cc512_mask_add_epi64(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_mask_add_epi64(to_native_512(src), static_cast<__mmask8>(mask), to_native_512(a), to_native_512(b)));
#else
        return mask_add_epi64_portable(src, mask, a, b);
#endif
    }

    BLST_FORCEINLINE auto _cc512_sub_epi64(const cc512i a, const cc512i b) -> cc512i { return sub_epi64(a, b); }

    // Masked subtraction: src if mask bit 0, a-b if mask bit 1
    BLST_FORCEINLINE auto _cc512_mask_sub_epi64(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_mask_sub_epi64(to_native_512(src), static_cast<__mmask8>(mask), to_native_512(a), to_native_512(b)));
#else
        return mask_sub_epi64_portable(src, mask, a, b);
#endif
    }

    BLST_FORCEINLINE auto _cc512_and_si512(const cc512i a, const cc512i b) -> cc512i { return and_si512(a, b); }
    BLST_FORCEINLINE auto _cc512_or_si512(const cc512i a, const cc512i b) -> cc512i { return or_si512(a, b); }

    // Bitwise XOR
    BLST_FORCEINLINE auto _cc512_xor_si512(const cc512i a, const cc512i b) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_xor_si512(to_native_512(a), to_native_512(b)));
#else
        return xor_si512_portable(a, b);
#endif
    }

    BLST_FORCEINLINE auto _cc512_slli_epi64(const cc512i a, int bits) -> cc512i { return slli_epi64(a, bits); }
    BLST_FORCEINLINE auto _cc512_srli_epi64(const cc512i a, int bits) -> cc512i { return srli_epi64(a, bits); }
    BLST_FORCEINLINE auto _cc512_srai_epi64(const cc512i a, int bits) -> cc512i { return srai_epi64(a, bits); }

    // Shift left logical (variable)
    BLST_FORCEINLINE auto _cc512_sllv_epi64(const cc512i a, const cc512i bits) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_sllv_epi64(to_native_512(a), to_native_512(bits)));
#else
        return sllv_epi64_portable(a, bits);
#endif
    }

    BLST_FORCEINLINE auto _cc512_load_si512(const void* p) -> cc512i { return load_si512(p); }
    BLST_FORCEINLINE auto _cc512_loadu_si512(const void* p) -> cc512i { return loadu_si512(p); }
    BLST_FORCEINLINE void _cc512_store_si512(void* p, const cc512i a) { store_si512(p, a); }
    BLST_FORCEINLINE void _cc512_storeu_si512(void* p, const cc512i a) { storeu_si512(p, a); }

    // Set 8 64-bit integers
    BLST_FORCEINLINE auto _cc512_set_epi64(
        long long x7, long long x6, long long x5, long long x4,
        long long x3, long long x2, long long x1, long long x0) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_set_epi64(x7, x6, x5, x4, x3, x2, x1, x0));
#else
        return set_epi64_portable(x7, x6, x5, x4, x3, x2, x1, x0);
#endif
    }

    // Permute using index vector
    BLST_FORCEINLINE auto _cc512_permutexvar_epi64(const cc512i idx, const cc512i a) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_permutexvar_epi64(to_native_512(idx), to_native_512(a)));
#else
        return permutexvar_epi64_portable(idx, a);
#endif
    }

    // Masked Blend
    BLST_FORCEINLINE auto _cc512_mask_blend_epi64(u64 mask, const cc512i a, const cc512i b) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_mask_blend_epi64(static_cast<__mmask8>(mask), to_native_512(a), to_native_512(b)));
#else
        return mask_blend_epi64_portable(mask, a, b);
#endif
    }

    BLST_FORCEINLINE auto _cc512_madd52lo_epu64(const cc512i acc, const cc512i a, const cc512i b) -> cc512i { return madd52lo_epu64(acc, a, b); }
    BLST_FORCEINLINE auto _cc512_madd52hi_epu64(const cc512i acc, const cc512i a, const cc512i b) -> cc512i { return madd52hi_epu64(acc, a, b); }

    // Shuffle 32-bit integers
    template <int imm8>
    BLST_FORCEINLINE auto _cc512_shuffle_epi32(const cc512i a) -> cc512i
    {
        static_assert((imm8 & ~0xFF) == 0, "imm8 must fit in 8 bits");
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_shuffle_epi32(to_native_512(a), static_cast<_MM_PERM_ENUM>(imm8)));
#else
        return shuffle_epi32_portable<imm8>(a);
#endif
    }

    // Masked shuffle 32-bit integers
    template <int imm8>
    BLST_FORCEINLINE auto _cc512_maskz_shuffle_epi32(u64 mask, const cc512i a) -> cc512i
    {
        static_assert((imm8 & ~0xFF) == 0, "imm8 must fit in 8 bits");
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_maskz_shuffle_epi32(static_cast<__mmask16>(mask), to_native_512(a), static_cast<_MM_PERM_ENUM>(imm8)));
#else
        return maskz_shuffle_epi32_portable<imm8>(mask, a);
#endif
    }

    // Permute 64-bit integers
    template <int imm8>
    BLST_FORCEINLINE auto _cc512_permutex_epi64(const cc512i a) -> cc512i
    {
        static_assert((imm8 & ~0xFF) == 0, "imm8 must fit in 8 bits");
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_permutex_epi64(to_native_512(a), imm8));
#else
        return permutex_epi64_portable<imm8>(a);
#endif
    }

    // Masked permute 64-bit integers
    template <int imm8>
    BLST_FORCEINLINE auto _cc512_maskz_permutex_epi64(u64 mask, const cc512i a) -> cc512i
    {
        static_assert((imm8 & ~0xFF) == 0, "imm8 must fit in 8 bits");
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_maskz_permutex_epi64(static_cast<__mmask8>(mask), to_native_512(a), imm8));
#else
        return maskz_permutex_epi64_portable<imm8>(mask, a);
#endif
    }

    // Extract 256-bit integer from 512-bit integer
    template <int imm8>
    BLST_FORCEINLINE auto _cc512_extracti64x4_epi64(const cc512i a) -> cc256i
    {
        static_assert(imm8 == 0 || imm8 == 1, "imm8 must be 0 or 1");
#if BLST_AVX_NATIVE
        return from_native_256(_mm512_extracti64x4_epi64(to_native_512(a), imm8));
#else
        return extracti64x4_epi64_portable<imm8>(a);
#endif
    }

    // Extract 64-bit integer from 256-bit integer
    template <int imm8>
    BLST_FORCEINLINE auto _cc256_extract_epi64(const cc256i a) -> u64
    {
        static_assert(imm8 >= 0 && imm8 < 4, "imm8 must be in [0,3]");
#if BLST_AVX_NATIVE
        return static_cast<u64>(_mm256_extract_epi64(to_native_256(a), imm8));
#else
        return extract_epi64_portable<imm8>(a);
#endif
    }
} // namespace blst::tile

#undef BLST_AVX_NATIVE

#endif

#endif // BLST_BACKEND_AVX512_INTRIN_COMPAT_H
