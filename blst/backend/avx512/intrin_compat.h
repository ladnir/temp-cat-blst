#ifndef BLST_BACKEND_AVX512_INTRIN_COMPAT_H
#define BLST_BACKEND_AVX512_INTRIN_COMPAT_H

#include "blst/common/defines.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#if defined(BLST_AVX512)
#include <immintrin.h>
#endif

namespace blst
{
#if defined(BLST_AVX512_PORTABLE)
    // Pure C/C++ representations used by the portable AVX-512 emulation.
    struct cc512i {
        u64 lane[8];
    };

    struct cc256i {
        u64 lane[4];
    };

    static_assert(sizeof(cc512i) == 64);
    static_assert(sizeof(cc256i) == 32);

    auto _cc512_setzero_si512() -> cc512i;
    auto _cc512_set1_epi64(long long x) -> cc512i;
    auto _cc512_add_epi64(const cc512i a, const cc512i b) -> cc512i;
    auto _cc512_sub_epi64(const cc512i a, const cc512i b) -> cc512i;
    auto _cc512_and_si512(const cc512i a, const cc512i b) -> cc512i;
    auto _cc512_or_si512(const cc512i a, const cc512i b) -> cc512i;
    auto _cc512_xor_si512(const cc512i a, const cc512i b) -> cc512i;
    auto _cc512_slli_epi64(const cc512i a, int bits) -> cc512i;
    auto _cc512_srli_epi64(const cc512i a, int bits) -> cc512i;
    auto _cc512_srai_epi64(const cc512i a, int bits) -> cc512i;
    auto _cc512_sllv_epi64(const cc512i a, const cc512i bits) -> cc512i;
    auto _cc512_load_si512(const void* p) -> cc512i;
    auto _cc512_loadu_si512(const void* p) -> cc512i;
    void _cc512_store_si512(void* p, const cc512i a);
    void _cc512_storeu_si512(void* p, const cc512i a);
    auto _cc512_set_epi64(
        long long x7, long long x6, long long x5, long long x4,
        long long x3, long long x2, long long x1, long long x0) -> cc512i;
    auto _cc512_permutexvar_epi64(const cc512i idx, const cc512i a) -> cc512i;
    auto _cc512_mask_blend_epi64(u64 mask, const cc512i a, const cc512i b) -> cc512i;
    auto _cc512_madd52lo_epu64(const cc512i acc, const cc512i a, const cc512i b) -> cc512i;
    auto _cc512_madd52hi_epu64(const cc512i acc, const cc512i a, const cc512i b) -> cc512i;
    auto _cc512_mask_add_epi64(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i;
    auto _cc512_mask_sub_epi64(const cc512i src, u64 mask, const cc512i a, const cc512i b) -> cc512i;
    auto _cc512_min_epu64(const cc512i a, const cc512i b) -> cc512i;
    auto _cc512_cmpgt_epi64_mask(const cc512i a, const cc512i b) -> u64;
    auto _cc512_cmpge_epi64_mask(const cc512i a, const cc512i b) -> u64;

    auto _cc512_shuffle_epi32_runtime(int imm8, const cc512i a) -> cc512i;
    auto _cc512_maskz_shuffle_epi32_runtime(int imm8, u64 mask, const cc512i a) -> cc512i;
    auto _cc512_permutex_epi64_runtime(int imm8, const cc512i a) -> cc512i;
    auto _cc512_maskz_permutex_epi64_runtime(int imm8, u64 mask, const cc512i a) -> cc512i;
    auto _cc512_extracti64x4_epi64_runtime(int imm8, const cc512i a) -> cc256i;
    auto _cc256_extract_epi64_runtime(int imm8, const cc256i a) -> u64;

    template <int imm8>
    BLST_FORCEINLINE auto _cc512_shuffle_epi32(const cc512i a) -> cc512i
    {
        static_assert((imm8 & ~0xff) == 0, "imm8 must fit in 8 bits");
        return _cc512_shuffle_epi32_runtime(imm8, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto _cc512_maskz_shuffle_epi32(u64 mask, const cc512i a) -> cc512i
    {
        static_assert((imm8 & ~0xff) == 0, "imm8 must fit in 8 bits");
        return _cc512_maskz_shuffle_epi32_runtime(imm8, mask, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto _cc512_permutex_epi64(const cc512i a) -> cc512i
    {
        static_assert((imm8 & ~0xff) == 0, "imm8 must fit in 8 bits");
        return _cc512_permutex_epi64_runtime(imm8, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto _cc512_maskz_permutex_epi64(u64 mask, const cc512i a) -> cc512i
    {
        static_assert((imm8 & ~0xff) == 0, "imm8 must fit in 8 bits");
        return _cc512_maskz_permutex_epi64_runtime(imm8, mask, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto _cc512_extracti64x4_epi64(const cc512i a) -> cc256i
    {
        static_assert(imm8 == 0 || imm8 == 1, "imm8 must be 0 or 1");
        return _cc512_extracti64x4_epi64_runtime(imm8, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto _cc256_extract_epi64(const cc256i a) -> u64
    {
        static_assert(imm8 >= 0 && imm8 < 4, "imm8 must be in [0,3]");
        return _cc256_extract_epi64_runtime(imm8, a);
    }
#endif

#if defined(BLST_AVX512)
#if defined(BLST_AVX512_PORTABLE)
    // Explicit bridges are for compatibility tests and intentional boundaries only.
    BLST_FORCEINLINE auto _cc512_from_mm512(const ::__m512i in) -> cc512i
    {
        cc512i out{};
        _mm512_storeu_si512(out.lane, in);
        return out;
    }

    BLST_FORCEINLINE auto _mm512_from_cc512(const cc512i in) -> ::__m512i
    {
        return _mm512_loadu_si512(in.lane);
    }

    BLST_FORCEINLINE auto _cc256_from_mm256(const ::__m256i in) -> cc256i
    {
        cc256i out{};
        _mm256_storeu_si256(reinterpret_cast<::__m256i*>(out.lane), in);
        return out;
    }

    BLST_FORCEINLINE auto _mm256_from_cc256(const cc256i in) -> ::__m256i
    {
        return _mm256_loadu_si256(reinterpret_cast<const ::__m256i*>(in.lane));
    }
#endif

    // Production-facing selected vector API. Native wins in combined builds.
    using vv512i = ::__m512i;
    using vv256i = ::__m256i;

    BLST_FORCEINLINE auto _vv512_setzero_si512() -> vv512i { return _mm512_setzero_si512(); }
    BLST_FORCEINLINE auto _vv512_set1_epi64(long long x) -> vv512i { return _mm512_set1_epi64(x); }
    BLST_FORCEINLINE auto _vv512_add_epi64(const vv512i a, const vv512i b) -> vv512i { return _mm512_add_epi64(a, b); }
    BLST_FORCEINLINE auto _vv512_sub_epi64(const vv512i a, const vv512i b) -> vv512i { return _mm512_sub_epi64(a, b); }
    BLST_FORCEINLINE auto _vv512_and_si512(const vv512i a, const vv512i b) -> vv512i { return _mm512_and_si512(a, b); }
    BLST_FORCEINLINE auto _vv512_or_si512(const vv512i a, const vv512i b) -> vv512i { return _mm512_or_si512(a, b); }
    BLST_FORCEINLINE auto _vv512_xor_si512(const vv512i a, const vv512i b) -> vv512i { return _mm512_xor_si512(a, b); }
    BLST_FORCEINLINE auto _vv512_slli_epi64(const vv512i a, int bits) -> vv512i { return _mm512_slli_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_srli_epi64(const vv512i a, int bits) -> vv512i { return _mm512_srli_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_srai_epi64(const vv512i a, int bits) -> vv512i { return _mm512_srai_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_sllv_epi64(const vv512i a, const vv512i bits) -> vv512i { return _mm512_sllv_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_load_si512(const void* p) -> vv512i { return _mm512_load_si512(p); }
    BLST_FORCEINLINE auto _vv512_loadu_si512(const void* p) -> vv512i { return _mm512_loadu_si512(p); }
    BLST_FORCEINLINE void _vv512_store_si512(void* p, const vv512i a) { _mm512_store_si512(p, a); }
    BLST_FORCEINLINE void _vv512_storeu_si512(void* p, const vv512i a) { _mm512_storeu_si512(p, a); }

    BLST_FORCEINLINE auto _vv512_set_epi64(
        long long x7, long long x6, long long x5, long long x4,
        long long x3, long long x2, long long x1, long long x0) -> vv512i
    {
        return _mm512_set_epi64(x7, x6, x5, x4, x3, x2, x1, x0);
    }

    BLST_FORCEINLINE auto _vv512_permutexvar_epi64(const vv512i idx, const vv512i a) -> vv512i
    {
        return _mm512_permutexvar_epi64(idx, a);
    }

    BLST_FORCEINLINE auto _vv512_mask_blend_epi64(u64 mask, const vv512i a, const vv512i b) -> vv512i
    {
        return _mm512_mask_blend_epi64(static_cast<__mmask8>(mask), a, b);
    }

    BLST_FORCEINLINE auto _vv512_madd52lo_epu64(const vv512i acc, const vv512i a, const vv512i b) -> vv512i
    {
        return _mm512_madd52lo_epu64(acc, a, b);
    }

    BLST_FORCEINLINE auto _vv512_madd52hi_epu64(const vv512i acc, const vv512i a, const vv512i b) -> vv512i
    {
        return _mm512_madd52hi_epu64(acc, a, b);
    }

    BLST_FORCEINLINE auto _vv512_mask_add_epi64(const vv512i src, u64 mask, const vv512i a, const vv512i b) -> vv512i
    {
        return _mm512_mask_add_epi64(src, static_cast<__mmask8>(mask), a, b);
    }

    BLST_FORCEINLINE auto _vv512_mask_sub_epi64(const vv512i src, u64 mask, const vv512i a, const vv512i b) -> vv512i
    {
        return _mm512_mask_sub_epi64(src, static_cast<__mmask8>(mask), a, b);
    }

    BLST_FORCEINLINE auto _vv512_min_epu64(const vv512i a, const vv512i b) -> vv512i { return _mm512_min_epu64(a, b); }

    BLST_FORCEINLINE auto _vv512_cmpgt_epi64_mask(const vv512i a, const vv512i b) -> u64
    {
        return static_cast<u64>(_mm512_cmpgt_epi64_mask(a, b));
    }

    BLST_FORCEINLINE auto _vv512_cmpge_epi64_mask(const vv512i a, const vv512i b) -> u64
    {
        return static_cast<u64>(_mm512_cmpge_epi64_mask(a, b));
    }

    template <int imm8>
    BLST_FORCEINLINE auto _vv512_shuffle_epi32(const vv512i a) -> vv512i
    {
        static_assert((imm8 & ~0xff) == 0, "imm8 must fit in 8 bits");
        return _mm512_shuffle_epi32(a, static_cast<_MM_PERM_ENUM>(imm8));
    }

    template <int imm8>
    BLST_FORCEINLINE auto _vv512_maskz_shuffle_epi32(u64 mask, const vv512i a) -> vv512i
    {
        static_assert((imm8 & ~0xff) == 0, "imm8 must fit in 8 bits");
        return _mm512_maskz_shuffle_epi32(static_cast<__mmask16>(mask), a, static_cast<_MM_PERM_ENUM>(imm8));
    }

    template <int imm8>
    BLST_FORCEINLINE auto _vv512_permutex_epi64(const vv512i a) -> vv512i
    {
        static_assert((imm8 & ~0xff) == 0, "imm8 must fit in 8 bits");
        return _mm512_permutex_epi64(a, imm8);
    }

    template <int imm8>
    BLST_FORCEINLINE auto _vv512_maskz_permutex_epi64(u64 mask, const vv512i a) -> vv512i
    {
        static_assert((imm8 & ~0xff) == 0, "imm8 must fit in 8 bits");
        return _mm512_maskz_permutex_epi64(static_cast<__mmask8>(mask), a, imm8);
    }

    template <int imm8>
    BLST_FORCEINLINE auto _vv512_extracti64x4_epi64(const vv512i a) -> vv256i
    {
        static_assert(imm8 == 0 || imm8 == 1, "imm8 must be 0 or 1");
        return _mm512_extracti64x4_epi64(a, imm8);
    }

    template <int imm8>
    BLST_FORCEINLINE auto _vv256_extract_epi64(const vv256i a) -> u64
    {
        static_assert(imm8 >= 0 && imm8 < 4, "imm8 must be in [0,3]");
        return static_cast<u64>(_mm256_extract_epi64(a, imm8));
    }
#elif defined(BLST_AVX512_PORTABLE)
    // Production-facing selected vector API when native AVX-512 is unavailable.
    using vv512i = cc512i;
    using vv256i = cc256i;

    BLST_FORCEINLINE auto _vv512_setzero_si512() -> vv512i { return _cc512_setzero_si512(); }
    BLST_FORCEINLINE auto _vv512_set1_epi64(long long x) -> vv512i { return _cc512_set1_epi64(x); }
    BLST_FORCEINLINE auto _vv512_add_epi64(const vv512i a, const vv512i b) -> vv512i { return _cc512_add_epi64(a, b); }
    BLST_FORCEINLINE auto _vv512_sub_epi64(const vv512i a, const vv512i b) -> vv512i { return _cc512_sub_epi64(a, b); }
    BLST_FORCEINLINE auto _vv512_and_si512(const vv512i a, const vv512i b) -> vv512i { return _cc512_and_si512(a, b); }
    BLST_FORCEINLINE auto _vv512_or_si512(const vv512i a, const vv512i b) -> vv512i { return _cc512_or_si512(a, b); }
    BLST_FORCEINLINE auto _vv512_xor_si512(const vv512i a, const vv512i b) -> vv512i { return _cc512_xor_si512(a, b); }
    BLST_FORCEINLINE auto _vv512_slli_epi64(const vv512i a, int bits) -> vv512i { return _cc512_slli_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_srli_epi64(const vv512i a, int bits) -> vv512i { return _cc512_srli_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_srai_epi64(const vv512i a, int bits) -> vv512i { return _cc512_srai_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_sllv_epi64(const vv512i a, const vv512i bits) -> vv512i { return _cc512_sllv_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_load_si512(const void* p) -> vv512i { return _cc512_load_si512(p); }
    BLST_FORCEINLINE auto _vv512_loadu_si512(const void* p) -> vv512i { return _cc512_loadu_si512(p); }
    BLST_FORCEINLINE void _vv512_store_si512(void* p, const vv512i a) { _cc512_store_si512(p, a); }
    BLST_FORCEINLINE void _vv512_storeu_si512(void* p, const vv512i a) { _cc512_storeu_si512(p, a); }
    BLST_FORCEINLINE auto _vv512_set_epi64(
        long long x7, long long x6, long long x5, long long x4,
        long long x3, long long x2, long long x1, long long x0) -> vv512i
    {
        return _cc512_set_epi64(x7, x6, x5, x4, x3, x2, x1, x0);
    }
    BLST_FORCEINLINE auto _vv512_permutexvar_epi64(const vv512i idx, const vv512i a) -> vv512i { return _cc512_permutexvar_epi64(idx, a); }
    BLST_FORCEINLINE auto _vv512_mask_blend_epi64(u64 mask, const vv512i a, const vv512i b) -> vv512i { return _cc512_mask_blend_epi64(mask, a, b); }
    BLST_FORCEINLINE auto _vv512_madd52lo_epu64(const vv512i acc, const vv512i a, const vv512i b) -> vv512i { return _cc512_madd52lo_epu64(acc, a, b); }
    BLST_FORCEINLINE auto _vv512_madd52hi_epu64(const vv512i acc, const vv512i a, const vv512i b) -> vv512i { return _cc512_madd52hi_epu64(acc, a, b); }
    BLST_FORCEINLINE auto _vv512_mask_add_epi64(const vv512i src, u64 mask, const vv512i a, const vv512i b) -> vv512i { return _cc512_mask_add_epi64(src, mask, a, b); }
    BLST_FORCEINLINE auto _vv512_mask_sub_epi64(const vv512i src, u64 mask, const vv512i a, const vv512i b) -> vv512i { return _cc512_mask_sub_epi64(src, mask, a, b); }
    BLST_FORCEINLINE auto _vv512_min_epu64(const vv512i a, const vv512i b) -> vv512i { return _cc512_min_epu64(a, b); }
    BLST_FORCEINLINE auto _vv512_cmpgt_epi64_mask(const vv512i a, const vv512i b) -> u64 { return _cc512_cmpgt_epi64_mask(a, b); }
    BLST_FORCEINLINE auto _vv512_cmpge_epi64_mask(const vv512i a, const vv512i b) -> u64 { return _cc512_cmpge_epi64_mask(a, b); }

    template <int imm8>
    BLST_FORCEINLINE auto _vv512_shuffle_epi32(const vv512i a) -> vv512i { return _cc512_shuffle_epi32<imm8>(a); }
    template <int imm8>
    BLST_FORCEINLINE auto _vv512_maskz_shuffle_epi32(u64 mask, const vv512i a) -> vv512i { return _cc512_maskz_shuffle_epi32<imm8>(mask, a); }
    template <int imm8>
    BLST_FORCEINLINE auto _vv512_permutex_epi64(const vv512i a) -> vv512i { return _cc512_permutex_epi64<imm8>(a); }
    template <int imm8>
    BLST_FORCEINLINE auto _vv512_maskz_permutex_epi64(u64 mask, const vv512i a) -> vv512i { return _cc512_maskz_permutex_epi64<imm8>(mask, a); }
    template <int imm8>
    BLST_FORCEINLINE auto _vv512_extracti64x4_epi64(const vv512i a) -> vv256i { return _cc512_extracti64x4_epi64<imm8>(a); }
    template <int imm8>
    BLST_FORCEINLINE auto _vv256_extract_epi64(const vv256i a) -> u64 { return _cc256_extract_epi64<imm8>(a); }
#endif
} // namespace blst

#endif

#endif // BLST_BACKEND_AVX512_INTRIN_COMPAT_H
