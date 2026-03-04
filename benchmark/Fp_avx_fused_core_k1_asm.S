.intel_syntax noprefix
.text

# blst_fp8x1_mul_redc_core_k1_asm
#
# Purpose:
# - Full fused k=1 kernel for benchmark exploration.
# - Input: two reduced FpPack8x1 values (a, b)
# - Output: one LFpPack8x1 (lazy reduced, no final sub-p)
#
# SysV x86_64:
# - rdi: LFpPack8x1* out
# - rsi: const FpPack8x1* a
# - rdx: const FpPack8x1* b
#
# Notes:
# - Implements the same arithmetic structure as the C++ fused baseline:
#   karatsuba-style mul_v2 + redc_core (lazy) without normalize.
# - Uses only zmm registers (no stack spills by construction).

.macro FUSED_K1_CORE_BODY
    # zero helpers
    vpxorq zmm16, zmm16, zmm16

    # -------------------------------------------------------------------------
    # mul_v2 low half: zmm0..zmm7 from (a0..a3)*(b0..b3)
    # aL: zmm24..zmm27, bL: zmm28..zmm31, carry-high chain in zmm17
    # -------------------------------------------------------------------------
    vmovdqa64 zmm24, ZMMWORD PTR [rsi +   0]
    vmovdqa64 zmm25, ZMMWORD PTR [rsi +  64]
    vmovdqa64 zmm26, ZMMWORD PTR [rsi + 128]
    vmovdqa64 zmm27, ZMMWORD PTR [rsi + 192]
    vmovdqa64 zmm28, ZMMWORD PTR [rdx +   0]
    vmovdqa64 zmm29, ZMMWORD PTR [rdx +  64]
    vmovdqa64 zmm30, ZMMWORD PTR [rdx + 128]
    vmovdqa64 zmm31, ZMMWORD PTR [rdx + 192]

    vpxorq zmm0, zmm0, zmm0
    vpmadd52luq zmm0, zmm24, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm28
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm1, zmm17
    vpmadd52luq zmm1, zmm24, zmm29
    vpmadd52luq zmm1, zmm25, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm29
    vpmadd52huq zmm17, zmm25, zmm28
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm2, zmm17
    vpmadd52luq zmm2, zmm24, zmm30
    vpmadd52luq zmm2, zmm25, zmm29
    vpmadd52luq zmm2, zmm26, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm30
    vpmadd52huq zmm17, zmm25, zmm29
    vpmadd52huq zmm17, zmm26, zmm28
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm3, zmm17
    vpmadd52luq zmm3, zmm24, zmm31
    vpmadd52luq zmm3, zmm25, zmm30
    vpmadd52luq zmm3, zmm26, zmm29
    vpmadd52luq zmm3, zmm27, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm31
    vpmadd52huq zmm17, zmm25, zmm30
    vpmadd52huq zmm17, zmm26, zmm29
    vpmadd52huq zmm17, zmm27, zmm28
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm4, zmm17
    vpmadd52luq zmm4, zmm25, zmm31
    vpmadd52luq zmm4, zmm26, zmm30
    vpmadd52luq zmm4, zmm27, zmm29
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm25, zmm31
    vpmadd52huq zmm17, zmm26, zmm30
    vpmadd52huq zmm17, zmm27, zmm29
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm5, zmm17
    vpmadd52luq zmm5, zmm26, zmm31
    vpmadd52luq zmm5, zmm27, zmm30
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm26, zmm31
    vpmadd52huq zmm17, zmm27, zmm30
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm6, zmm17
    vpmadd52luq zmm6, zmm27, zmm31
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm27, zmm31
    vpsllq zmm17, zmm17, 4
    vmovdqa64 zmm7, zmm17

    # -------------------------------------------------------------------------
    # mul_v2 high half: zmm8..zmm15 from (a4..a7)*(b4..b7)
    # aH: zmm24..zmm27, bH: zmm28..zmm31, carry-high chain in zmm17
    # -------------------------------------------------------------------------
    vmovdqa64 zmm24, ZMMWORD PTR [rsi + 256]
    vmovdqa64 zmm25, ZMMWORD PTR [rsi + 320]
    vmovdqa64 zmm26, ZMMWORD PTR [rsi + 384]
    vmovdqa64 zmm27, ZMMWORD PTR [rsi + 448]
    vmovdqa64 zmm28, ZMMWORD PTR [rdx + 256]
    vmovdqa64 zmm29, ZMMWORD PTR [rdx + 320]
    vmovdqa64 zmm30, ZMMWORD PTR [rdx + 384]
    vmovdqa64 zmm31, ZMMWORD PTR [rdx + 448]

    vpxorq zmm8, zmm8, zmm8
    vpmadd52luq zmm8, zmm24, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm28
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm9, zmm17
    vpmadd52luq zmm9, zmm24, zmm29
    vpmadd52luq zmm9, zmm25, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm29
    vpmadd52huq zmm17, zmm25, zmm28
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm10, zmm17
    vpmadd52luq zmm10, zmm24, zmm30
    vpmadd52luq zmm10, zmm25, zmm29
    vpmadd52luq zmm10, zmm26, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm30
    vpmadd52huq zmm17, zmm25, zmm29
    vpmadd52huq zmm17, zmm26, zmm28
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm11, zmm17
    vpmadd52luq zmm11, zmm24, zmm31
    vpmadd52luq zmm11, zmm25, zmm30
    vpmadd52luq zmm11, zmm26, zmm29
    vpmadd52luq zmm11, zmm27, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm31
    vpmadd52huq zmm17, zmm25, zmm30
    vpmadd52huq zmm17, zmm26, zmm29
    vpmadd52huq zmm17, zmm27, zmm28
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm12, zmm17
    vpmadd52luq zmm12, zmm25, zmm31
    vpmadd52luq zmm12, zmm26, zmm30
    vpmadd52luq zmm12, zmm27, zmm29
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm25, zmm31
    vpmadd52huq zmm17, zmm26, zmm30
    vpmadd52huq zmm17, zmm27, zmm29
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm13, zmm17
    vpmadd52luq zmm13, zmm26, zmm31
    vpmadd52luq zmm13, zmm27, zmm30
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm26, zmm31
    vpmadd52huq zmm17, zmm27, zmm30
    vpsllq zmm17, zmm17, 4

    vmovdqa64 zmm14, zmm17
    vpmadd52luq zmm14, zmm27, zmm31
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm27, zmm31
    vpsllq zmm17, zmm17, 4
    vmovdqa64 zmm15, zmm17

    # -------------------------------------------------------------------------
    # middle product zM = (aL+aH)*(bL+bH) - zL - zH
    # ta: zmm24..zmm27, tb: zmm28..zmm31
    # m temp: zmm18, carry-high temp: zmm17, sum temp: zmm19
    # -------------------------------------------------------------------------
    vmovdqa64 zmm24, ZMMWORD PTR [rsi +   0]
    vpaddq zmm24, zmm24, ZMMWORD PTR [rsi + 256]
    vmovdqa64 zmm25, ZMMWORD PTR [rsi +  64]
    vpaddq zmm25, zmm25, ZMMWORD PTR [rsi + 320]
    vmovdqa64 zmm26, ZMMWORD PTR [rsi + 128]
    vpaddq zmm26, zmm26, ZMMWORD PTR [rsi + 384]
    vmovdqa64 zmm27, ZMMWORD PTR [rsi + 192]
    vpaddq zmm27, zmm27, ZMMWORD PTR [rsi + 448]

    vmovdqa64 zmm28, ZMMWORD PTR [rdx +   0]
    vpaddq zmm28, zmm28, ZMMWORD PTR [rdx + 256]
    vmovdqa64 zmm29, ZMMWORD PTR [rdx +  64]
    vpaddq zmm29, zmm29, ZMMWORD PTR [rdx + 320]
    vmovdqa64 zmm30, ZMMWORD PTR [rdx + 128]
    vpaddq zmm30, zmm30, ZMMWORD PTR [rdx + 384]
    vmovdqa64 zmm31, ZMMWORD PTR [rdx + 192]
    vpaddq zmm31, zmm31, ZMMWORD PTR [rdx + 448]

    # preserve original z4..z7 for m4..m7 correction terms
    vmovdqa64 zmm20, zmm4
    vmovdqa64 zmm21, zmm5
    vmovdqa64 zmm22, zmm6
    vmovdqa64 zmm23, zmm7

    # m0
    vpxorq zmm18, zmm18, zmm18
    vpmadd52luq zmm18, zmm24, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm28
    vpsllq zmm17, zmm17, 4
    vmovdqa64 zmm19, zmm0
    vpaddq zmm19, zmm19, zmm8
    vpsubq zmm18, zmm18, zmm19
    vpaddq zmm4, zmm4, zmm18

    # m1
    vmovdqa64 zmm18, zmm17
    vpmadd52luq zmm18, zmm24, zmm29
    vpmadd52luq zmm18, zmm25, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm29
    vpmadd52huq zmm17, zmm25, zmm28
    vpsllq zmm17, zmm17, 4
    vmovdqa64 zmm19, zmm1
    vpaddq zmm19, zmm19, zmm9
    vpsubq zmm18, zmm18, zmm19
    vpaddq zmm5, zmm5, zmm18

    # m2
    vmovdqa64 zmm18, zmm17
    vpmadd52luq zmm18, zmm24, zmm30
    vpmadd52luq zmm18, zmm25, zmm29
    vpmadd52luq zmm18, zmm26, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm30
    vpmadd52huq zmm17, zmm25, zmm29
    vpmadd52huq zmm17, zmm26, zmm28
    vpsllq zmm17, zmm17, 4
    vmovdqa64 zmm19, zmm2
    vpaddq zmm19, zmm19, zmm10
    vpsubq zmm18, zmm18, zmm19
    vpaddq zmm6, zmm6, zmm18

    # m3
    vmovdqa64 zmm18, zmm17
    vpmadd52luq zmm18, zmm24, zmm31
    vpmadd52luq zmm18, zmm25, zmm30
    vpmadd52luq zmm18, zmm26, zmm29
    vpmadd52luq zmm18, zmm27, zmm28
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm24, zmm31
    vpmadd52huq zmm17, zmm25, zmm30
    vpmadd52huq zmm17, zmm26, zmm29
    vpmadd52huq zmm17, zmm27, zmm28
    vpsllq zmm17, zmm17, 4
    vmovdqa64 zmm19, zmm3
    vpaddq zmm19, zmm19, zmm11
    vpsubq zmm18, zmm18, zmm19
    vpaddq zmm7, zmm7, zmm18

    # m4
    vmovdqa64 zmm18, zmm17
    vpmadd52luq zmm18, zmm25, zmm31
    vpmadd52luq zmm18, zmm26, zmm30
    vpmadd52luq zmm18, zmm27, zmm29
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm25, zmm31
    vpmadd52huq zmm17, zmm26, zmm30
    vpmadd52huq zmm17, zmm27, zmm29
    vpsllq zmm17, zmm17, 4
    vmovdqa64 zmm19, zmm20
    vpaddq zmm19, zmm19, zmm12
    vpsubq zmm18, zmm18, zmm19
    vpaddq zmm8, zmm8, zmm18

    # m5
    vmovdqa64 zmm18, zmm17
    vpmadd52luq zmm18, zmm26, zmm31
    vpmadd52luq zmm18, zmm27, zmm30
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm26, zmm31
    vpmadd52huq zmm17, zmm27, zmm30
    vpsllq zmm17, zmm17, 4
    vmovdqa64 zmm19, zmm21
    vpaddq zmm19, zmm19, zmm13
    vpsubq zmm18, zmm18, zmm19
    vpaddq zmm9, zmm9, zmm18

    # m6
    vmovdqa64 zmm18, zmm17
    vpmadd52luq zmm18, zmm27, zmm31
    vpxorq zmm17, zmm17, zmm17
    vpmadd52huq zmm17, zmm27, zmm31
    vpsllq zmm17, zmm17, 4
    vmovdqa64 zmm19, zmm22
    vpaddq zmm19, zmm19, zmm14
    vpsubq zmm18, zmm18, zmm19
    vpaddq zmm10, zmm10, zmm18

    # m7
    vmovdqa64 zmm18, zmm17
    vmovdqa64 zmm19, zmm23
    vpaddq zmm19, zmm19, zmm15
    vpsubq zmm18, zmm18, zmm19
    vpaddq zmm11, zmm11, zmm18

    # -------------------------------------------------------------------------
    # redc_core on zmm0..zmm15 -> lazy zmm8..zmm15
    # y0..y12 in zmm16..zmm28, montw/bmask in zmm29/zmm30, u/tmp in zmm31
    # -------------------------------------------------------------------------
    vpxorq zmm16, zmm16, zmm16
    vpxorq zmm17, zmm17, zmm17
    vpxorq zmm18, zmm18, zmm18
    vpxorq zmm19, zmm19, zmm19
    vpxorq zmm20, zmm20, zmm20
    vpxorq zmm21, zmm21, zmm21
    vpxorq zmm22, zmm22, zmm22
    vpxorq zmm23, zmm23, zmm23
    vpxorq zmm24, zmm24, zmm24
    vpxorq zmm25, zmm25, zmm25
    vpxorq zmm26, zmm26, zmm26
    vpxorq zmm27, zmm27, zmm27
    vpxorq zmm28, zmm28, zmm28
    vmovdqa64 zmm29, ZMMWORD PTR [rip + montw_vec]
    vmovdqa64 zmm30, ZMMWORD PTR [rip + bmask_vec]

    # round 0
    vpxorq zmm31, zmm31, zmm31
    vpmadd52luq zmm31, zmm0, zmm29
    vpandq zmm31, zmm31, zmm30
    vpmadd52luq zmm0, zmm31, ZMMWORD PTR [rip + p0_vec]
    vpmadd52luq zmm1, zmm31, ZMMWORD PTR [rip + p1_vec]
    vpmadd52luq zmm2, zmm31, ZMMWORD PTR [rip + p2_vec]
    vpmadd52luq zmm3, zmm31, ZMMWORD PTR [rip + p3_vec]
    vpmadd52luq zmm4, zmm31, ZMMWORD PTR [rip + p4_vec]
    vpmadd52luq zmm5, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52luq zmm6, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52luq zmm7, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpmadd52huq zmm16, zmm31, ZMMWORD PTR [rip + p0_vec]
    vpmadd52huq zmm17, zmm31, ZMMWORD PTR [rip + p1_vec]
    vpmadd52huq zmm18, zmm31, ZMMWORD PTR [rip + p2_vec]
    vpmadd52huq zmm19, zmm31, ZMMWORD PTR [rip + p3_vec]
    vpmadd52huq zmm20, zmm31, ZMMWORD PTR [rip + p4_vec]
    vpmadd52huq zmm21, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52huq zmm22, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52huq zmm23, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpsraq zmm31, zmm0, 48
    vpaddq zmm1, zmm1, zmm31
    vpsllq zmm31, zmm16, 4
    vpaddq zmm1, zmm1, zmm31
    vmovdqa64 zmm0, ZMMWORD PTR [rip + p0_vec]

    # round 1
    vpxorq zmm31, zmm31, zmm31
    vpmadd52luq zmm31, zmm1, zmm29
    vpandq zmm31, zmm31, zmm30
    vpmadd52luq zmm1, zmm31, zmm0
    vpmadd52luq zmm2, zmm31, ZMMWORD PTR [rip + p1_vec]
    vpmadd52luq zmm3, zmm31, ZMMWORD PTR [rip + p2_vec]
    vpmadd52luq zmm4, zmm31, ZMMWORD PTR [rip + p3_vec]
    vpmadd52luq zmm5, zmm31, ZMMWORD PTR [rip + p4_vec]
    vpmadd52luq zmm6, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52luq zmm7, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52luq zmm8, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpmadd52huq zmm17, zmm31, zmm0
    vpmadd52huq zmm18, zmm31, ZMMWORD PTR [rip + p1_vec]
    vpmadd52huq zmm19, zmm31, ZMMWORD PTR [rip + p2_vec]
    vpmadd52huq zmm20, zmm31, ZMMWORD PTR [rip + p3_vec]
    vpmadd52huq zmm21, zmm31, ZMMWORD PTR [rip + p4_vec]
    vpmadd52huq zmm22, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52huq zmm23, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52huq zmm24, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpsraq zmm31, zmm1, 48
    vpaddq zmm2, zmm2, zmm31
    vpsllq zmm31, zmm17, 4
    vpaddq zmm2, zmm2, zmm31
    vmovdqa64 zmm1, ZMMWORD PTR [rip + p1_vec]

    # round 2
    vpxorq zmm31, zmm31, zmm31
    vpmadd52luq zmm31, zmm2, zmm29
    vpandq zmm31, zmm31, zmm30
    vpmadd52luq zmm2, zmm31, zmm0
    vpmadd52luq zmm3, zmm31, zmm1
    vpmadd52luq zmm4, zmm31, ZMMWORD PTR [rip + p2_vec]
    vpmadd52luq zmm5, zmm31, ZMMWORD PTR [rip + p3_vec]
    vpmadd52luq zmm6, zmm31, ZMMWORD PTR [rip + p4_vec]
    vpmadd52luq zmm7, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52luq zmm8, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52luq zmm9, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpmadd52huq zmm18, zmm31, zmm0
    vpmadd52huq zmm19, zmm31, zmm1
    vpmadd52huq zmm20, zmm31, ZMMWORD PTR [rip + p2_vec]
    vpmadd52huq zmm21, zmm31, ZMMWORD PTR [rip + p3_vec]
    vpmadd52huq zmm22, zmm31, ZMMWORD PTR [rip + p4_vec]
    vpmadd52huq zmm23, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52huq zmm24, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52huq zmm25, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpsraq zmm31, zmm2, 48
    vpaddq zmm3, zmm3, zmm31
    vpsllq zmm31, zmm18, 4
    vpaddq zmm3, zmm3, zmm31
    vmovdqa64 zmm2, ZMMWORD PTR [rip + p2_vec]

    # round 3
    vpxorq zmm31, zmm31, zmm31
    vpmadd52luq zmm31, zmm3, zmm29
    vpandq zmm31, zmm31, zmm30
    vpmadd52luq zmm3, zmm31, zmm0
    vpmadd52luq zmm4, zmm31, zmm1
    vpmadd52luq zmm5, zmm31, zmm2
    vpmadd52luq zmm6, zmm31, ZMMWORD PTR [rip + p3_vec]
    vpmadd52luq zmm7, zmm31, ZMMWORD PTR [rip + p4_vec]
    vpmadd52luq zmm8, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52luq zmm9, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52luq zmm10, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpmadd52huq zmm19, zmm31, zmm0
    vpmadd52huq zmm20, zmm31, zmm1
    vpmadd52huq zmm21, zmm31, zmm2
    vpmadd52huq zmm22, zmm31, ZMMWORD PTR [rip + p3_vec]
    vpmadd52huq zmm23, zmm31, ZMMWORD PTR [rip + p4_vec]
    vpmadd52huq zmm24, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52huq zmm25, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52huq zmm26, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpsraq zmm31, zmm3, 48
    vpaddq zmm4, zmm4, zmm31
    vpsllq zmm31, zmm19, 4
    vpaddq zmm4, zmm4, zmm31
    vmovdqa64 zmm3, ZMMWORD PTR [rip + p3_vec]

    # round 4
    vpxorq zmm31, zmm31, zmm31
    vpmadd52luq zmm31, zmm4, zmm29
    vpandq zmm31, zmm31, zmm30
    vpmadd52luq zmm4, zmm31, zmm0
    vpmadd52luq zmm5, zmm31, zmm1
    vpmadd52luq zmm6, zmm31, zmm2
    vpmadd52luq zmm7, zmm31, zmm3
    vpmadd52luq zmm8, zmm31, ZMMWORD PTR [rip + p4_vec]
    vpmadd52luq zmm9, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52luq zmm10, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52luq zmm11, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpmadd52huq zmm20, zmm31, zmm0
    vpmadd52huq zmm21, zmm31, zmm1
    vpmadd52huq zmm22, zmm31, zmm2
    vpmadd52huq zmm23, zmm31, zmm3
    vpmadd52huq zmm24, zmm31, ZMMWORD PTR [rip + p4_vec]
    vpmadd52huq zmm25, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52huq zmm26, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52huq zmm27, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpsraq zmm31, zmm4, 48
    vpaddq zmm5, zmm5, zmm31
    vpsllq zmm31, zmm20, 4
    vpaddq zmm5, zmm5, zmm31
    vmovdqa64 zmm4, ZMMWORD PTR [rip + p4_vec]

    # round 5
    vpxorq zmm31, zmm31, zmm31
    vpmadd52luq zmm31, zmm5, zmm29
    vpandq zmm31, zmm31, zmm30
    vpmadd52luq zmm5, zmm31, zmm0
    vpmadd52luq zmm6, zmm31, zmm1
    vpmadd52luq zmm7, zmm31, zmm2
    vpmadd52luq zmm8, zmm31, zmm3
    vpmadd52luq zmm9, zmm31, zmm4
    vpmadd52luq zmm10, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52luq zmm11, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52luq zmm12, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpmadd52huq zmm21, zmm31, zmm0
    vpmadd52huq zmm22, zmm31, zmm1
    vpmadd52huq zmm23, zmm31, zmm2
    vpmadd52huq zmm24, zmm31, zmm3
    vpmadd52huq zmm25, zmm31, zmm4
    vpmadd52huq zmm26, zmm31, ZMMWORD PTR [rip + p5_vec]
    vpmadd52huq zmm27, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52huq zmm28, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpsraq zmm31, zmm5, 48
    vpaddq zmm6, zmm6, zmm31
    vpsllq zmm31, zmm21, 4
    vpaddq zmm6, zmm6, zmm31
    vmovdqa64 zmm5, ZMMWORD PTR [rip + p5_vec]

    # round 6
    vpxorq zmm31, zmm31, zmm31
    vpmadd52luq zmm31, zmm6, zmm29
    vpandq zmm31, zmm31, zmm30
    vpmadd52luq zmm6, zmm31, zmm0
    vpmadd52luq zmm7, zmm31, zmm1
    vpmadd52luq zmm8, zmm31, zmm2
    vpmadd52luq zmm9, zmm31, zmm3
    vpmadd52luq zmm10, zmm31, zmm4
    vpmadd52luq zmm11, zmm31, zmm5
    vpmadd52luq zmm12, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpmadd52luq zmm13, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpmadd52huq zmm22, zmm31, zmm0
    vpmadd52huq zmm23, zmm31, zmm1
    vpmadd52huq zmm24, zmm31, zmm2
    vpmadd52huq zmm25, zmm31, zmm3
    vpmadd52huq zmm26, zmm31, zmm4
    vpmadd52huq zmm27, zmm31, zmm5
    vpmadd52huq zmm28, zmm31, ZMMWORD PTR [rip + p6_vec]
    vpxorq zmm30, zmm30, zmm30
    vpmadd52huq zmm30, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpsraq zmm31, zmm6, 48
    vpaddq zmm7, zmm7, zmm31
    vpsllq zmm31, zmm22, 4
    vpaddq zmm7, zmm7, zmm31
    vpsllq zmm30, zmm30, 4
    vpaddq zmm14, zmm14, zmm30
    vmovdqa64 zmm6, ZMMWORD PTR [rip + p6_vec]
    vmovdqa64 zmm30, ZMMWORD PTR [rip + bmask_vec]

    # round 7
    vpxorq zmm31, zmm31, zmm31
    vpmadd52luq zmm31, zmm7, zmm29
    vpandq zmm31, zmm31, zmm30
    vpmadd52luq zmm7, zmm31, zmm0
    vpmadd52luq zmm8, zmm31, zmm1
    vpmadd52luq zmm9, zmm31, zmm2
    vpmadd52luq zmm10, zmm31, zmm3
    vpmadd52luq zmm11, zmm31, zmm4
    vpmadd52luq zmm12, zmm31, zmm5
    vpmadd52luq zmm13, zmm31, zmm6
    vpmadd52luq zmm14, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpmadd52huq zmm23, zmm31, zmm0
    vpmadd52huq zmm24, zmm31, zmm1
    vpmadd52huq zmm25, zmm31, zmm2
    vpmadd52huq zmm26, zmm31, zmm3
    vpmadd52huq zmm27, zmm31, zmm4
    vpmadd52huq zmm28, zmm31, zmm5
    vpxorq zmm30, zmm30, zmm30
    vpmadd52huq zmm30, zmm31, zmm6
    vpxorq zmm29, zmm29, zmm29
    vpmadd52huq zmm29, zmm31, ZMMWORD PTR [rip + p7_vec]
    vpsraq zmm31, zmm7, 48
    vpaddq zmm8, zmm8, zmm31
    vpsllq zmm31, zmm23, 4
    vpaddq zmm8, zmm8, zmm31
    vpsllq zmm30, zmm30, 4
    vpaddq zmm14, zmm14, zmm30
    vpsllq zmm29, zmm29, 4
    vpaddq zmm15, zmm15, zmm29
    vmovdqa64 zmm30, ZMMWORD PTR [rip + bmask_vec]

    # tail add of pending high terms
    vpsllq zmm31, zmm24, 4
    vpaddq zmm9, zmm9, zmm31
    vpsllq zmm31, zmm25, 4
    vpaddq zmm10, zmm10, zmm31
    vpsllq zmm31, zmm26, 4
    vpaddq zmm11, zmm11, zmm31
    vpsllq zmm31, zmm27, 4
    vpaddq zmm12, zmm12, zmm31
    vpsllq zmm31, zmm28, 4
    vpaddq zmm13, zmm13, zmm31

    # final carry propagation on zmm8..zmm15
    vpsraq zmm31, zmm8, 48
    vpaddq zmm9, zmm9, zmm31
    vpandq zmm8, zmm8, zmm30
    vpsraq zmm31, zmm9, 48
    vpaddq zmm10, zmm10, zmm31
    vpandq zmm9, zmm9, zmm30
    vpsraq zmm31, zmm10, 48
    vpaddq zmm11, zmm11, zmm31
    vpandq zmm10, zmm10, zmm30
    vpsraq zmm31, zmm11, 48
    vpaddq zmm12, zmm12, zmm31
    vpandq zmm11, zmm11, zmm30
    vpsraq zmm31, zmm12, 48
    vpaddq zmm13, zmm13, zmm31
    vpandq zmm12, zmm12, zmm30
    vpsraq zmm31, zmm13, 48
    vpaddq zmm14, zmm14, zmm31
    vpandq zmm13, zmm13, zmm30
    vpsraq zmm31, zmm14, 48
    vpaddq zmm15, zmm15, zmm31
    vpandq zmm14, zmm14, zmm30
    vpandq zmm15, zmm15, zmm30

    # store lazy output
    vmovdqa64 ZMMWORD PTR [rdi +   0], zmm8
    vmovdqa64 ZMMWORD PTR [rdi +  64], zmm9
    vmovdqa64 ZMMWORD PTR [rdi + 128], zmm10
    vmovdqa64 ZMMWORD PTR [rdi + 192], zmm11
    vmovdqa64 ZMMWORD PTR [rdi + 256], zmm12
    vmovdqa64 ZMMWORD PTR [rdi + 320], zmm13
    vmovdqa64 ZMMWORD PTR [rdi + 384], zmm14
    vmovdqa64 ZMMWORD PTR [rdi + 448], zmm15
.endm

.globl blst_fp8x1_mul_redc_core_k1_asm
.type blst_fp8x1_mul_redc_core_k1_asm, @function
blst_fp8x1_mul_redc_core_k1_asm:
    FUSED_K1_CORE_BODY
    ret

.size blst_fp8x1_mul_redc_core_k1_asm, .-blst_fp8x1_mul_redc_core_k1_asm

# blst_fp8x1_mul_redc_core_k1_asm_span
#
# SysV x86_64:
# - rdi: LFpPack8x1* out
# - rsi: const FpPack8x1* a
# - rdx: const FpPack8x1* b
# - rcx: size_t n_packs
.globl blst_fp8x1_mul_redc_core_k1_asm_span
.type blst_fp8x1_mul_redc_core_k1_asm_span, @function
blst_fp8x1_mul_redc_core_k1_asm_span:
    test rcx, rcx
    je .Lk1_span_done

.Lk1_span_loop:
    FUSED_K1_CORE_BODY
    add rdi, 512
    add rsi, 512
    add rdx, 512
    dec rcx
    jne .Lk1_span_loop

.Lk1_span_done:
    ret

.size blst_fp8x1_mul_redc_core_k1_asm_span, .-blst_fp8x1_mul_redc_core_k1_asm_span

.section .rodata
.p2align 6

bmask_vec:
    .quad 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF
    .quad 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF

montw_vec:
    .quad 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD
    .quad 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD

p0_vec:
    .quad 0xFFFFFFFFAAAB, 0xFFFFFFFFAAAB, 0xFFFFFFFFAAAB, 0xFFFFFFFFAAAB
    .quad 0xFFFFFFFFAAAB, 0xFFFFFFFFAAAB, 0xFFFFFFFFAAAB, 0xFFFFFFFFAAAB
p1_vec:
    .quad 0xB153FFFFB9FE, 0xB153FFFFB9FE, 0xB153FFFFB9FE, 0xB153FFFFB9FE
    .quad 0xB153FFFFB9FE, 0xB153FFFFB9FE, 0xB153FFFFB9FE, 0xB153FFFFB9FE
p2_vec:
    .quad 0xF6241EABFFFE, 0xF6241EABFFFE, 0xF6241EABFFFE, 0xF6241EABFFFE
    .quad 0xF6241EABFFFE, 0xF6241EABFFFE, 0xF6241EABFFFE, 0xF6241EABFFFE
p3_vec:
    .quad 0x6730D2A0F6B0, 0x6730D2A0F6B0, 0x6730D2A0F6B0, 0x6730D2A0F6B0
    .quad 0x6730D2A0F6B0, 0x6730D2A0F6B0, 0x6730D2A0F6B0, 0x6730D2A0F6B0
p4_vec:
    .quad 0x4B84F38512BF, 0x4B84F38512BF, 0x4B84F38512BF, 0x4B84F38512BF
    .quad 0x4B84F38512BF, 0x4B84F38512BF, 0x4B84F38512BF, 0x4B84F38512BF
p5_vec:
    .quad 0x434BACD76477, 0x434BACD76477, 0x434BACD76477, 0x434BACD76477
    .quad 0x434BACD76477, 0x434BACD76477, 0x434BACD76477, 0x434BACD76477
p6_vec:
    .quad 0xE69A4B1BA7B6, 0xE69A4B1BA7B6, 0xE69A4B1BA7B6, 0xE69A4B1BA7B6
    .quad 0xE69A4B1BA7B6, 0xE69A4B1BA7B6, 0xE69A4B1BA7B6, 0xE69A4B1BA7B6
p7_vec:
    .quad 0x1A0111EA397F, 0x1A0111EA397F, 0x1A0111EA397F, 0x1A0111EA397F
    .quad 0x1A0111EA397F, 0x1A0111EA397F, 0x1A0111EA397F, 0x1A0111EA397F

.att_syntax prefix
