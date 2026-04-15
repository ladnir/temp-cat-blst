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
    auto setzero_si512_portable() -> cc512i;
    auto set1_epi64_portable(long long x) -> cc512i;
    auto add_epi64_portable(const cc512i a, const cc512i b) -> cc512i;
    auto sub_epi64_portable(const cc512i a, const cc512i b) -> cc512i;
    auto and_si512_portable(const cc512i a, const cc512i b) -> cc512i;
    auto or_si512_portable(const cc512i a, const cc512i b) -> cc512i;
    auto xor_si512_portable(const cc512i a, const cc512i b) -> cc512i;
    auto slli_epi64_portable(const cc512i a, int bits) -> cc512i;
    auto srli_epi64_portable(const cc512i a, int bits) -> cc512i;
    auto srai_epi64_portable(const cc512i a, int bits) -> cc512i;
    auto sllv_epi64_portable(const cc512i a, const cc512i bits) -> cc512i;
    auto load_si512_portable(const void* p) -> cc512i;
    auto loadu_si512_portable(const void* p) -> cc512i;
    void store_si512_portable(void* p, const cc512i a);
    void storeu_si512_portable(void* p, const cc512i a);
    auto set_epi64_portable(
        long long x7, long long x6, long long x5, long long x4,
        long long x3, long long x2, long long x1, long long x0) -> cc512i;
    auto permutexvar_epi64_portable(const cc512i idx, const cc512i a) -> cc512i;
    auto mask_blend_epi64_portable(u64 mask, const cc512i a, const cc512i b) -> cc512i;
    auto madd52lo_epu64_portable(const cc512i acc, const cc512i a, const cc512i b) -> cc512i;
    auto madd52hi_epu64_portable(const cc512i acc, const cc512i a, const cc512i b) -> cc512i;
    auto mask_add_epi64_portable(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i;
    auto mask_sub_epi64_portable(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i;
    auto min_epu64_portable(const cc512i a, const cc512i b) -> cc512i;
    auto cmpgt_epi64_mask_portable(const cc512i a, const cc512i b) -> u64;
    auto cmpge_epi64_mask_portable(const cc512i a, const cc512i b) -> u64;
    auto shuffle_epi32_portable_runtime(int imm8, const cc512i a) -> cc512i;
    auto maskz_shuffle_epi32_portable_runtime(int imm8, u64 mask, const cc512i a) -> cc512i;
    auto permutex_epi64_portable_runtime(int imm8, const cc512i a) -> cc512i;
    auto maskz_permutex_epi64_portable_runtime(int imm8, u64 mask, const cc512i a) -> cc512i;
    auto extracti64x4_epi64_portable_runtime(int imm8, const cc512i a) -> cc256i;
    auto extract_epi64_portable_runtime(int imm8, const cc256i a) -> u64;

    template <int imm8>
    BLST_FORCEINLINE auto shuffle_epi32_portable(const cc512i a) -> cc512i
    {
        return shuffle_epi32_portable_runtime(imm8, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto maskz_shuffle_epi32_portable(u64 mask, const cc512i a) -> cc512i
    {
        return maskz_shuffle_epi32_portable_runtime(imm8, mask, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto permutex_epi64_portable(const cc512i a) -> cc512i
    {
        return permutex_epi64_portable_runtime(imm8, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto maskz_permutex_epi64_portable(u64 mask, const cc512i a) -> cc512i
    {
        return maskz_permutex_epi64_portable_runtime(imm8, mask, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto extracti64x4_epi64_portable(const cc512i a) -> cc256i
    {
        return extracti64x4_epi64_portable_runtime(imm8, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto extract_epi64_portable(const cc256i a) -> u64
    {
        return extract_epi64_portable_runtime(imm8, a);
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

    BLST_FORCEINLINE auto _cc512_min_epu64(const cc512i a, const cc512i b) -> cc512i
    {
#if BLST_AVX_NATIVE
        return from_native_512(_mm512_min_epu64(to_native_512(a), to_native_512(b)));
#else
        return min_epu64_portable(a, b);
#endif
    }

    BLST_FORCEINLINE auto _cc512_cmpgt_epi64_mask(const cc512i a, const cc512i b) -> u64
    {
#if BLST_AVX_NATIVE
        return static_cast<u64>(_mm512_cmpgt_epi64_mask(to_native_512(a), to_native_512(b)));
#else
        return cmpgt_epi64_mask_portable(a, b);
#endif
    }

    BLST_FORCEINLINE auto _cc512_cmpge_epi64_mask(const cc512i a, const cc512i b) -> u64
    {
#if BLST_AVX_NATIVE
        return static_cast<u64>(_mm512_cmpge_epi64_mask(to_native_512(a), to_native_512(b)));
#else
        return cmpge_epi64_mask_portable(a, b);
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
