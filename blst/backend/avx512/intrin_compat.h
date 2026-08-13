// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// Selected-vector AVX-512 operations and their pure C/C++ emulation.

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
    /**
     * Pure C/C++ representations of the retained AVX-512 operation surface.
     *
     * `cc` means a C/C++ implementation. Scalar `cc` code is a supported
     * production backend. These vector `cc` operations specifically emulate
     * AVX-512 for developer testing and are not a deployment backend or a
     * constant-time code-generation claim.
     *
     * Lane zero is the least-significant lane used by Intel intrinsics. The
     * structures own ordinary lane arrays and never contain native vector
     * types. Every operation is allocation-free and nonthrowing by
     * construction.
     */
    struct cc512i {
        u64 lane[8];
    };

    struct cc256i {
        u64 lane[4];
    };

    struct cc128i {
        u64 lane[2];
    };

    static_assert(sizeof(cc512i) == 64);
    static_assert(sizeof(cc256i) == 32);
    static_assert(sizeof(cc128i) == 16);
    static_assert(alignof(cc512i) == alignof(u64));
    static_assert(alignof(cc256i) == alignof(u64));
    static_assert(alignof(cc128i) == alignof(u64));
    static_assert(std::is_standard_layout_v<cc512i>);
    static_assert(std::is_standard_layout_v<cc256i>);
    static_assert(std::is_standard_layout_v<cc128i>);
    static_assert(std::is_trivially_copyable_v<cc512i>);
    static_assert(std::is_trivially_copyable_v<cc256i>);
    static_assert(std::is_trivially_copyable_v<cc128i>);

    auto _cc512_set1_epi64(u64 x) -> cc512i;
    auto _cc512_add_epi64(cc512i a, cc512i b) -> cc512i;
    auto _cc512_sub_epi64(cc512i a, cc512i b) -> cc512i;
    auto _cc512_and_si512(cc512i a, cc512i b) -> cc512i;
    auto _cc512_slli_epi64(cc512i a, int bits) -> cc512i;
    auto _cc512_srli_epi64(cc512i a, int bits) -> cc512i;

    // `p` must name at least 64 readable or writable bytes. No alignment is
    // required, and the complete 64-byte object is transferred.
    auto _cc512_loadu_si512(const void* p) -> cc512i;
    void _cc512_storeu_si512(void* p, cc512i a);

    auto _cc512_permutexvar_epi64(cc512i idx, cc512i a) -> cc512i;
    auto _cc512_mask_blend_epi64(u64 mask, cc512i a, cc512i b) -> cc512i;
    auto _cc512_madd52lo_epu64(cc512i acc, cc512i a, cc512i b) -> cc512i;
    auto _cc512_madd52hi_epu64(cc512i acc, cc512i a, cc512i b) -> cc512i;
    auto _cc512_min_epu64(cc512i a, cc512i b) -> cc512i;
    auto _cc512_cmpeq_epi64_mask(cc512i a, cc512i b) -> u64;
    auto _cc512_extracti64x4_epi64_runtime(int imm8, cc512i a) -> cc256i;
    auto _cc256_extract_epi64_runtime(int imm8, cc256i a) -> u64;
    auto _cc512_castsi512_si128(cc512i a) -> cc128i;
    auto _cc128_add_epi64(cc128i a, cc128i b) -> cc128i;
    auto _cc128_srli_epi64(cc128i a, int bits) -> cc128i;
    auto _cc512_broadcastq_epi64(cc128i a) -> cc512i;
    auto _cc128_extract_epi64_runtime(int imm8, cc128i a) -> u64;

    template <int imm8>
    BLST_FORCEINLINE auto _cc512_extracti64x4_epi64(cc512i a) -> cc256i
    {
        static_assert(imm8 == 0 || imm8 == 1, "imm8 must be 0 or 1");
        return _cc512_extracti64x4_epi64_runtime(imm8, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto _cc256_extract_epi64(cc256i a) -> u64
    {
        static_assert(imm8 >= 0 && imm8 < 4, "imm8 must be in [0,3]");
        return _cc256_extract_epi64_runtime(imm8, a);
    }

    template <int imm8>
    BLST_FORCEINLINE auto _cc128_extract_epi64(cc128i a) -> u64
    {
        static_assert(imm8 == 0 || imm8 == 1, "imm8 must be 0 or 1");
        return _cc128_extract_epi64_runtime(imm8, a);
    }
#endif

#if defined(BLST_AVX512)
#if defined(BLST_AVX512_PORTABLE)
    // Explicit bridges are test-only comparisons between the independent
    // native and C/C++ domains. Production selected-vector code must not use
    // them.
    BLST_FORCEINLINE auto _cc512_from_mm512(::__m512i in) -> cc512i
    {
        cc512i out{};
        _mm512_storeu_si512(out.lane, in);
        return out;
    }

    BLST_FORCEINLINE auto _mm512_from_cc512(cc512i in) -> ::__m512i
    {
        return _mm512_loadu_si512(in.lane);
    }

    BLST_FORCEINLINE auto _cc256_from_mm256(::__m256i in) -> cc256i
    {
        cc256i out{};
        _mm256_storeu_si256(reinterpret_cast<::__m256i*>(out.lane), in);
        return out;
    }

    BLST_FORCEINLINE auto _mm256_from_cc256(cc256i in) -> ::__m256i
    {
        return _mm256_loadu_si256(reinterpret_cast<const ::__m256i*>(in.lane));
    }

    BLST_FORCEINLINE auto _cc128_from_mm128(::__m128i in) -> cc128i
    {
        cc128i out{};
        _mm_storeu_si128(reinterpret_cast<::__m128i*>(out.lane), in);
        return out;
    }

    BLST_FORCEINLINE auto _mm128_from_cc128(cc128i in) -> ::__m128i
    {
        return _mm_loadu_si128(reinterpret_cast<const ::__m128i*>(in.lane));
    }
#endif

    /**
     * Compile-time-selected vector operations used by production AVX code.
     * Native AVX-512 always wins when both configurations are enabled. The
     * caller must execute native code only on a CPU and OS context supporting
     * AVX-512F, AVX-512DQ, and AVX-512IFMA.
     *
     * Masks use their low eight bits. Immediate shifts follow Intel semantics:
     * counts outside [0,63] produce zero. Unaligned transfers require 64 valid
     * bytes. The wrappers add no dispatch, conversion, allocation, or exception
     * boundary.
     */
    static_assert(sizeof(long long) == sizeof(u64));
    using vv512i = ::__m512i;
    using vv256i = ::__m256i;
    using vv128i = ::__m128i;

    BLST_FORCEINLINE auto _vv512_set1_epi64(u64 x) -> vv512i
    {
        return _mm512_set1_epi64(std::bit_cast<long long>(x));
    }
    BLST_FORCEINLINE auto _vv512_add_epi64(vv512i a, vv512i b) -> vv512i { return _mm512_add_epi64(a, b); }
    BLST_FORCEINLINE auto _vv512_sub_epi64(vv512i a, vv512i b) -> vv512i { return _mm512_sub_epi64(a, b); }
    BLST_FORCEINLINE auto _vv512_and_si512(vv512i a, vv512i b) -> vv512i { return _mm512_and_si512(a, b); }
    BLST_FORCEINLINE auto _vv512_slli_epi64(vv512i a, int bits) -> vv512i { return _mm512_slli_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_srli_epi64(vv512i a, int bits) -> vv512i { return _mm512_srli_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_loadu_si512(const void* p) -> vv512i { return _mm512_loadu_si512(p); }
    BLST_FORCEINLINE void _vv512_storeu_si512(void* p, vv512i a) { _mm512_storeu_si512(p, a); }
    BLST_FORCEINLINE auto _vv512_permutexvar_epi64(vv512i idx, vv512i a) -> vv512i { return _mm512_permutexvar_epi64(idx, a); }
    BLST_FORCEINLINE auto _vv512_mask_blend_epi64(u64 mask, vv512i a, vv512i b) -> vv512i
    {
        return _mm512_mask_blend_epi64(static_cast<__mmask8>(mask), a, b);
    }
    BLST_FORCEINLINE auto _vv512_madd52lo_epu64(vv512i acc, vv512i a, vv512i b) -> vv512i
    {
        return _mm512_madd52lo_epu64(acc, a, b);
    }
    BLST_FORCEINLINE auto _vv512_madd52hi_epu64(vv512i acc, vv512i a, vv512i b) -> vv512i
    {
        return _mm512_madd52hi_epu64(acc, a, b);
    }
    BLST_FORCEINLINE auto _vv512_min_epu64(vv512i a, vv512i b) -> vv512i { return _mm512_min_epu64(a, b); }
    BLST_FORCEINLINE auto _vv512_cmpeq_epi64_mask(vv512i a, vv512i b) -> u64
    {
        return static_cast<u64>(_mm512_cmpeq_epi64_mask(a, b));
    }
    BLST_FORCEINLINE auto _vv512_castsi512_si128(vv512i a) -> vv128i { return _mm512_castsi512_si128(a); }
    BLST_FORCEINLINE auto _vv128_add_epi64(vv128i a, vv128i b) -> vv128i { return _mm_add_epi64(a, b); }
    BLST_FORCEINLINE auto _vv128_srli_epi64(vv128i a, int bits) -> vv128i { return _mm_srli_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_broadcastq_epi64(vv128i a) -> vv512i { return _mm512_broadcastq_epi64(a); }

    template <int imm8>
    BLST_FORCEINLINE auto _vv512_extracti64x4_epi64(vv512i a) -> vv256i
    {
        static_assert(imm8 == 0 || imm8 == 1, "imm8 must be 0 or 1");
        return _mm512_extracti64x4_epi64(a, imm8);
    }

    template <int imm8>
    BLST_FORCEINLINE auto _vv256_extract_epi64(vv256i a) -> u64
    {
        static_assert(imm8 >= 0 && imm8 < 4, "imm8 must be in [0,3]");
        return static_cast<u64>(_mm256_extract_epi64(a, imm8));
    }

    template <int imm8>
    BLST_FORCEINLINE auto _vv128_extract_epi64(vv128i a) -> u64
    {
        static_assert(imm8 == 0 || imm8 == 1, "imm8 must be 0 or 1");
        if constexpr (imm8 == 0) {
            return static_cast<u64>(_mm_cvtsi128_si64(a));
        } else {
            return static_cast<u64>(_mm_extract_epi64(a, 1));
        }
    }
#elif defined(BLST_AVX512_PORTABLE)
    // Selected-vector aliases for the developer-only AVX emulation build.
    using vv512i = cc512i;
    using vv256i = cc256i;
    using vv128i = cc128i;

    BLST_FORCEINLINE auto _vv512_set1_epi64(u64 x) -> vv512i { return _cc512_set1_epi64(x); }
    BLST_FORCEINLINE auto _vv512_add_epi64(vv512i a, vv512i b) -> vv512i { return _cc512_add_epi64(a, b); }
    BLST_FORCEINLINE auto _vv512_sub_epi64(vv512i a, vv512i b) -> vv512i { return _cc512_sub_epi64(a, b); }
    BLST_FORCEINLINE auto _vv512_and_si512(vv512i a, vv512i b) -> vv512i { return _cc512_and_si512(a, b); }
    BLST_FORCEINLINE auto _vv512_slli_epi64(vv512i a, int bits) -> vv512i { return _cc512_slli_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_srli_epi64(vv512i a, int bits) -> vv512i { return _cc512_srli_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_loadu_si512(const void* p) -> vv512i { return _cc512_loadu_si512(p); }
    BLST_FORCEINLINE void _vv512_storeu_si512(void* p, vv512i a) { _cc512_storeu_si512(p, a); }
    BLST_FORCEINLINE auto _vv512_permutexvar_epi64(vv512i idx, vv512i a) -> vv512i { return _cc512_permutexvar_epi64(idx, a); }
    BLST_FORCEINLINE auto _vv512_mask_blend_epi64(u64 mask, vv512i a, vv512i b) -> vv512i { return _cc512_mask_blend_epi64(mask, a, b); }
    BLST_FORCEINLINE auto _vv512_madd52lo_epu64(vv512i acc, vv512i a, vv512i b) -> vv512i { return _cc512_madd52lo_epu64(acc, a, b); }
    BLST_FORCEINLINE auto _vv512_madd52hi_epu64(vv512i acc, vv512i a, vv512i b) -> vv512i { return _cc512_madd52hi_epu64(acc, a, b); }
    BLST_FORCEINLINE auto _vv512_min_epu64(vv512i a, vv512i b) -> vv512i { return _cc512_min_epu64(a, b); }
    BLST_FORCEINLINE auto _vv512_cmpeq_epi64_mask(vv512i a, vv512i b) -> u64 { return _cc512_cmpeq_epi64_mask(a, b); }
    BLST_FORCEINLINE auto _vv512_castsi512_si128(vv512i a) -> vv128i { return _cc512_castsi512_si128(a); }
    BLST_FORCEINLINE auto _vv128_add_epi64(vv128i a, vv128i b) -> vv128i { return _cc128_add_epi64(a, b); }
    BLST_FORCEINLINE auto _vv128_srli_epi64(vv128i a, int bits) -> vv128i { return _cc128_srli_epi64(a, bits); }
    BLST_FORCEINLINE auto _vv512_broadcastq_epi64(vv128i a) -> vv512i { return _cc512_broadcastq_epi64(a); }

    template <int imm8>
    BLST_FORCEINLINE auto _vv512_extracti64x4_epi64(vv512i a) -> vv256i { return _cc512_extracti64x4_epi64<imm8>(a); }
    template <int imm8>
    BLST_FORCEINLINE auto _vv256_extract_epi64(vv256i a) -> u64 { return _cc256_extract_epi64<imm8>(a); }
    template <int imm8>
    BLST_FORCEINLINE auto _vv128_extract_epi64(vv128i a) -> u64 { return _cc128_extract_epi64<imm8>(a); }
#endif

    // Source-literal vector construction for generated immutable constants.
    // Keep this separate from the intrinsic-shaped runtime operation surface.
    constexpr auto _vv512_const_u64(
        u64 x0, u64 x1, u64 x2, u64 x3,
        u64 x4, u64 x5, u64 x6, u64 x7) -> vv512i
    {
#if defined(BLST_AVX512) && defined(_MSC_VER)
        vv512i out = {};
        out.m512i_u64[0] = x0;
        out.m512i_u64[1] = x1;
        out.m512i_u64[2] = x2;
        out.m512i_u64[3] = x3;
        out.m512i_u64[4] = x4;
        out.m512i_u64[5] = x5;
        out.m512i_u64[6] = x6;
        out.m512i_u64[7] = x7;
        return out;
#elif defined(BLST_AVX512)
        return vv512i{
            std::bit_cast<long long>(x0), std::bit_cast<long long>(x1),
            std::bit_cast<long long>(x2), std::bit_cast<long long>(x3),
            std::bit_cast<long long>(x4), std::bit_cast<long long>(x5),
            std::bit_cast<long long>(x6), std::bit_cast<long long>(x7)};
#else
        return vv512i{{x0, x1, x2, x3, x4, x5, x6, x7}};
#endif
    }
} // namespace blst

#endif

#endif // BLST_BACKEND_AVX512_INTRIN_COMPAT_H
