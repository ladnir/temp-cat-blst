#ifndef BLST_FIELDS_FP_AVX_H
#define BLST_FIELDS_FP_AVX_H

#include "blst/common/defines.h"
#include "blst/backend/avx512/avx_ops.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    constexpr u64 BRADIX = 48;
    constexpr u64 NWORDS = 8;
    constexpr u64 VWORDS = 4;
    constexpr u64 SWORDS = 6;
    constexpr u64 BMASK = 0xFFFFFFFFFFFFULL;
    constexpr u64 FMASK = 0xFFFFFFFFFFFFFFFFULL;
    constexpr u64 BALIGN = 4;
    constexpr u64 MONT_W_R48 = 0xFFFCFFFCFFFDULL;

    extern const u64 P48[NWORDS];

    using fp_8x1w = cc512i[NWORDS];
    using fpx2_8x1w = cc512i[2 * NWORDS];
    using fp_4x2w = cc512i[VWORDS];
    using fpx2_4x2w = cc512i[3 * VWORDS];

    // Conversion/layout helpers
    void conv_64to48_fp(u64* r, const u64* a);
    void conv_64to48_fp_8x1w(fp_8x1w r, const cc512i* a);
    void conv_64to48_fp_4x2w(fp_4x2w r, const cc512i* a);
    void conv_64to48_mpi(u64* r, const u64* a, int rlen, int alen);

    void conv_48to64_fp(u64* r, const u64* a);
    void conv_48to64_fp_8x1w(cc512i* r, const fp_8x1w a);
    void conv_48to64_fp_4x2w(cc512i* r, const fp_4x2w a);
    void conv_48to64_mpi(u64* r, const u64* a, int rlen, int alen);

    void carryp_mpi48(u64* a);
    void carryp_fp_4x2w(fp_4x2w a);
    void get_channel_8x1w(u64* r, const fp_8x1w a, int ch);
    void get_channel_4x2w(u64* r, const fp_4x2w a, int ch);

    // Phase-1 Fp helpers
    void add_fp_8x1w(fp_8x1w r, const fp_8x1w a, const fp_8x1w b);
    void sub_fp_8x1w(fp_8x1w r, const fp_8x1w a, const fp_8x1w b);
    void asx4_fp_8x1w(fp_8x1w r, const fp_8x1w a, const fp_8x1w b);
    void asax6_fp_8x1w(fp_8x1w r, const fp_8x1w a, const fp_8x1w b);
    void assx6_fp_8x1w(fp_8x1w r, const fp_8x1w a, const fp_8x1w b);

    void add_fp_4x2w(fp_4x2w r, const fp_4x2w a, const fp_4x2w b);
    void sub_fp_4x2w(fp_4x2w r, const fp_4x2w a, const fp_4x2w b);
    void asx2_fp_4x2w(fp_4x2w r, const fp_4x2w a, const fp_4x2w b);

    void add_fpx2_8x1w(fpx2_8x1w r, const fpx2_8x1w a, const fpx2_8x1w b);
    void sub_fpx2_8x1w(fpx2_8x1w r, const fpx2_8x1w a, const fpx2_8x1w b);
    void aaasx2_fpx2_8x1w(fpx2_8x1w r, const fpx2_8x1w a, const fpx2_8x1w b);
    void aassx2_fpx2_8x1w(fpx2_8x1w r, const fpx2_8x1w a, const fpx2_8x1w b);
    void asx4_fpx2_8x1w(fpx2_8x1w r, const fpx2_8x1w a, const fpx2_8x1w b);
    void ssaax2_fpx2_8x1w(fpx2_8x1w r, const fpx2_8x1w a, const fpx2_8x1w b);
    void ssasx2_fpx2_8x1w(fpx2_8x1w r, const fpx2_8x1w a, const fpx2_8x1w b);

    void add_fpx2_4x2w(fpx2_4x2w r, const fpx2_4x2w a, const fpx2_4x2w b);
    void sub_fpx2_4x2w(fpx2_4x2w r, const fpx2_4x2w a, const fpx2_4x2w b);
    void asx2_fpx2_4x2w(fpx2_4x2w r, const fpx2_4x2w a, const fpx2_4x2w b);
    void aass_fpx2_4x2w(fpx2_4x2w r, const fpx2_4x2w a, const fpx2_4x2w b);

    // Phase-1 Fp kernels
    void mul_fpx2_8x1w_v1(fpx2_8x1w r, const fp_8x1w a, const fp_8x1w b);
    void mul_fpx2_8x1w_v2(fpx2_8x1w r, const fp_8x1w a, const fp_8x1w b);
    void mul_fpx2_8x1w_v3(fpx2_8x1w r, const fp_8x1w a, const fp_8x1w b);
    void mul_fpx2_8x1w_v4(fpx2_8x1w r, const fp_8x1w a, const fp_8x1w b);
    void mul_fpx2_4x2w_v1(fpx2_4x2w r, const fp_4x2w a, const fp_4x2w b);

    void redc_fpx2_8x1w(fp_8x1w r, const fpx2_8x1w a);
    void redc_fpx2_4x2w(fp_4x2w r, const fpx2_4x2w a);
}

#endif

#endif // BLST_FIELDS_FP_AVX_H
