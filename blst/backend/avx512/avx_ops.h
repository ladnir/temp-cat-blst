#ifndef BLST_BACKEND_AVX512_AVX_OPS_H
#define BLST_BACKEND_AVX512_AVX_OPS_H

#include "blst/backend/avx512/intrin_compat.h"
#include "blst/common/defines.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    inline auto VSET(i64 x7, i64 x6, i64 x5, i64 x4, i64 x3, i64 x2, i64 x1, i64 x0) -> cc512i
    {
        return _cc512_set_epi64(x7, x6, x5, x4, x3, x2, x1, x0);
    }

    inline auto VZERO() -> cc512i
    {
        return _cc512_setzero_si512();
    }

    inline auto VSET1(i64 x) -> cc512i
    {
        return _cc512_set1_epi64(x);
    }

    inline auto VPERMV(const cc512i idx, const cc512i a) -> cc512i
    {
        return _cc512_permutexvar_epi64(idx, a);
    }

    template <int imm8>
    inline auto VPERM(const cc512i a) -> cc512i
    {
        return _cc512_permutex_epi64<imm8>(a);
    }

    template <int imm8>
    inline auto VZPERM(u64 mask, const cc512i a) -> cc512i
    {
        return _cc512_maskz_permutex_epi64<imm8>(mask, a);
    }

    template <int imm8>
    inline auto VSHUF(const cc512i a) -> cc512i
    {
        return _cc512_shuffle_epi32<imm8>(a);
    }


    template <int imm8>
    inline auto VZSHUF(u64 mask, const cc512i a) -> cc512i
    {
        return _cc512_maskz_shuffle_epi32<imm8>(mask, a);
    }

    inline auto VMBLEND(u64 mask, const cc512i a, const cc512i b) -> cc512i
    {
        return _cc512_mask_blend_epi64(mask, a, b);
    }

    inline auto VADD(const cc512i a, const cc512i b) -> cc512i
    {
        return _cc512_add_epi64(a, b);
    }

    inline auto VMADD(const cc512i w, u64 mask, const cc512i y, const cc512i z) -> cc512i
    {
        return _cc512_mask_add_epi64(w, mask, y, z);
    }

    inline auto VSUB(const cc512i a, const cc512i b) -> cc512i
    {
        return _cc512_sub_epi64(a, b);
    }

    inline auto VMSUB(const cc512i w, u64 mask, const cc512i y, const cc512i z) -> cc512i
    {
        return _cc512_mask_sub_epi64(w, mask, y, z);
    }

    inline auto VMACLO(const cc512i x, const cc512i y, const cc512i z) -> cc512i
    {
        return _cc512_madd52lo_epu64(x, y, z);
    }

    inline auto VMACHI(const cc512i x, const cc512i y, const cc512i z) -> cc512i
    {
        return _cc512_madd52hi_epu64(x, y, z);
    }

    inline auto VAND(const cc512i a, const cc512i b) -> cc512i
    {
        return _cc512_and_si512(a, b);
    }

    inline auto VOR(const cc512i a, const cc512i b) -> cc512i
    {
        return _cc512_or_si512(a, b);
    }

    inline auto VXOR(const cc512i a, const cc512i b) -> cc512i
    {
        return _cc512_xor_si512(a, b);
    }

    inline auto VSHR(const cc512i a, int bits) -> cc512i
    {
        return _cc512_srli_epi64(a, bits);
    }

    inline auto VSHL(const cc512i a, int bits) -> cc512i
    {
        return _cc512_slli_epi64(a, bits);
    }

    inline auto VSRA(const cc512i a, int bits) -> cc512i
    {
        return _cc512_srai_epi64(a, bits);
    }

    template <int imm8>
    inline auto VEXTR256(const cc512i a) -> cc256i
    {
        return _cc512_extracti64x4_epi64<imm8>(a);
    }

    inline auto VEXTR256(const cc512i a, int imm8) -> cc256i
    {
        switch (imm8) {
        case 0: return VEXTR256<0>(a);
        case 1: return VEXTR256<1>(a);
        default:
            BLST_ASSERT(false);
            return _cc512_extracti64x4_epi64<0>(a);
        }
    }

    template <int imm8>
    inline auto VEXTR64(const cc256i a) -> u64
    {
        return _cc256_extract_epi64<imm8>(a);
    }

    inline auto VEXTR64(const cc256i a, int imm8) -> u64
    {
        switch (imm8) {
        case 0: return VEXTR64<0>(a);
        case 1: return VEXTR64<1>(a);
        case 2: return VEXTR64<2>(a);
        case 3: return VEXTR64<3>(a);
        default:
            BLST_ASSERT(false);
            return 0;
        }
    }
}

#endif

#endif // BLST_BACKEND_AVX512_AVX_OPS_H
