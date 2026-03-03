#if defined(__x86_64__)

.intel_syntax noprefix
.text

.globl blst_redc_core_k2_asm
.type blst_redc_core_k2_asm, @function
.globl blst_redc_core_k2_asm_v2
.type blst_redc_core_k2_asm_v2, @function

.macro REDC_ROUND U, C, A0, A1, A2, A3, A4, A5, A6, A7, Y0, Y1, Y2, Y3, Y4, Y5, Y6, BASE, OFF
    vpxorq \U, \U, \U
    vpmadd52luq \U, \A0, ZMMWORD PTR [rip + montw_vec]
    vpandq \U, \U, ZMMWORD PTR [rip + bmask_vec]

    # new queue tail slot (logical y[r+6]) starts from zero every round
    vpxorq \Y6, \Y6, \Y6

    vmovdqa64 \C, ZMMWORD PTR [rip + p0_vec]
    vpmadd52luq \A0, \U, \C
    vpmadd52huq \Y0, \U, \C
    vmovdqa64 \C, ZMMWORD PTR [rip + p1_vec]
    vpmadd52luq \A1, \U, \C
    vpmadd52huq \Y1, \U, \C
    vmovdqa64 \C, ZMMWORD PTR [rip + p2_vec]
    vpmadd52luq \A2, \U, \C
    vpmadd52huq \Y2, \U, \C
    vmovdqa64 \C, ZMMWORD PTR [rip + p3_vec]
    vpmadd52luq \A3, \U, \C
    vpmadd52huq \Y3, \U, \C
    vmovdqa64 \C, ZMMWORD PTR [rip + p4_vec]
    vpmadd52luq \A4, \U, \C
    vpmadd52huq \Y4, \U, \C
    vmovdqa64 \C, ZMMWORD PTR [rip + p5_vec]
    vpmadd52luq \A5, \U, \C
    vpmadd52huq \Y5, \U, \C
    vmovdqa64 \C, ZMMWORD PTR [rip + p6_vec]
    vpmadd52luq \A6, \U, \C
    vpmadd52huq \Y6, \U, \C

    # a[r+1] += carry(a[r]) + (y[r] << 4)
    vpsraq \A0, \A0, 48
    vpaddq \A1, \A1, \A0
    vpsllq \Y0, \Y0, 4
    vpaddq \A1, \A1, \Y0

    # j=7: eager fold into next limb (this removes one persistent y slot)
    vmovdqa64 \A0, ZMMWORD PTR [\BASE + \OFF]
    vmovdqa64 \C, ZMMWORD PTR [rip + p7_vec]
    vpmadd52luq \A0, \U, \C
    vpxorq \Y0, \Y0, \Y0
    vpmadd52huq \Y0, \U, \C
    vpsllq \Y0, \Y0, 4
    vpaddq \A0, \A0, \Y0
.endm

.macro REDC_CARRY8 TMP, BM, R0, R1, R2, R3, R4, R5, R6, R7
    vpsraq \TMP, \R0, 48
    vpaddq \R1, \R1, \TMP
    vpandq \R0, \R0, \BM

    vpsraq \TMP, \R1, 48
    vpaddq \R2, \R2, \TMP
    vpandq \R1, \R1, \BM

    vpsraq \TMP, \R2, 48
    vpaddq \R3, \R3, \TMP
    vpandq \R2, \R2, \BM

    vpsraq \TMP, \R3, 48
    vpaddq \R4, \R4, \TMP
    vpandq \R3, \R3, \BM

    vpsraq \TMP, \R4, 48
    vpaddq \R5, \R5, \TMP
    vpandq \R4, \R4, \BM

    vpsraq \TMP, \R5, 48
    vpaddq \R6, \R6, \TMP
    vpandq \R5, \R5, \BM

    vpsraq \TMP, \R6, 48
    vpaddq \R7, \R7, \TMP
    vpandq \R6, \R6, \BM

    vpandq \R7, \R7, \BM
.endm

.macro REDC_CARRY8_X2 TMP0, TMP1, BM, A0, A1, A2, A3, A4, A5, A6, A7, B0, B1, B2, B3, B4, B5, B6, B7
    vpsraq \TMP0, \A0, 48
    vpaddq \A1, \A1, \TMP0
    vpandq \A0, \A0, \BM
    vpsraq \TMP1, \B0, 48
    vpaddq \B1, \B1, \TMP1
    vpandq \B0, \B0, \BM

    vpsraq \TMP0, \A1, 48
    vpaddq \A2, \A2, \TMP0
    vpandq \A1, \A1, \BM
    vpsraq \TMP1, \B1, 48
    vpaddq \B2, \B2, \TMP1
    vpandq \B1, \B1, \BM

    vpsraq \TMP0, \A2, 48
    vpaddq \A3, \A3, \TMP0
    vpandq \A2, \A2, \BM
    vpsraq \TMP1, \B2, 48
    vpaddq \B3, \B3, \TMP1
    vpandq \B2, \B2, \BM

    vpsraq \TMP0, \A3, 48
    vpaddq \A4, \A4, \TMP0
    vpandq \A3, \A3, \BM
    vpsraq \TMP1, \B3, 48
    vpaddq \B4, \B4, \TMP1
    vpandq \B3, \B3, \BM

    vpsraq \TMP0, \A4, 48
    vpaddq \A5, \A5, \TMP0
    vpandq \A4, \A4, \BM
    vpsraq \TMP1, \B4, 48
    vpaddq \B5, \B5, \TMP1
    vpandq \B4, \B4, \BM

    vpsraq \TMP0, \A5, 48
    vpaddq \A6, \A6, \TMP0
    vpandq \A5, \A5, \BM
    vpsraq \TMP1, \B5, 48
    vpaddq \B6, \B6, \TMP1
    vpandq \B5, \B5, \BM

    vpsraq \TMP0, \A6, 48
    vpaddq \A7, \A7, \TMP0
    vpandq \A6, \A6, \BM
    vpsraq \TMP1, \B6, 48
    vpaddq \B7, \B7, \TMP1
    vpandq \B6, \B6, \BM

    vpandq \A7, \A7, \BM
    vpandq \B7, \B7, \BM
.endm

.macro REDC_ROUND_BOTH6 U0, U1, C, T, A0, A1, A2, A3, A4, A5, A6, A7, B0, B1, B2, B3, B4, B5, B6, B7, YH, Y1, Y2, Y3, Y4, Y5, ZH, Z1, Z2, Z3, Z4, Z5, BASEA, OFFA, BASEB, OFFB
    # u0/u1 from pivot limbs
    vpxorq \U0, \U0, \U0
    vpbroadcastq \C, QWORD PTR [rip + montw_q]
    vpmadd52luq \U0, \A0, \C
    vpbroadcastq \C, QWORD PTR [rip + bmask_q]
    vpandq \U0, \U0, \C
    vpxorq \U1, \U1, \U1
    vpbroadcastq \C, QWORD PTR [rip + montw_q]
    vpmadd52luq \U1, \B0, \C
    vpbroadcastq \C, QWORD PTR [rip + bmask_q]
    vpandq \U1, \U1, \C

    # p0 low updates
    vpbroadcastq \C, QWORD PTR [rip + p0_q]
    vpmadd52luq \A0, \U0, \C
    vpmadd52luq \B0, \U1, \C

    # a[r+1] / b[r+1] += carry(pivot)
    vpsraq \T, \A0, 48
    vpaddq \A1, \A1, \T
    vpsraq \T, \B0, 48
    vpaddq \B1, \B1, \T

    # add old queue head
    vpsllq \T, \YH, 4
    vpaddq \A1, \A1, \T
    vpsllq \T, \ZH, 4
    vpaddq \B1, \B1, \T

    # add p0 high directly to next limb
    vpmadd52huq \T, \U0, \C
    vpsllq \T, \T, 4
    vpaddq \A1, \A1, \T
    vpmadd52huq \T, \U1, \C
    vpsllq \T, \T, 4
    vpaddq \B1, \B1, \T

    # queue head is consumed; tail slot is this physical head register
    vpxorq \YH, \YH, \YH
    vpxorq \ZH, \ZH, \ZH

    # p1..p6 low/high updates
    vpbroadcastq \C, QWORD PTR [rip + p1_q]
    vpmadd52luq \A1, \U0, \C
    vpmadd52huq \Y1, \U0, \C
    vpmadd52luq \B1, \U1, \C
    vpmadd52huq \Z1, \U1, \C

    vpbroadcastq \C, QWORD PTR [rip + p2_q]
    vpmadd52luq \A2, \U0, \C
    vpmadd52huq \Y2, \U0, \C
    vpmadd52luq \B2, \U1, \C
    vpmadd52huq \Z2, \U1, \C

    vpbroadcastq \C, QWORD PTR [rip + p3_q]
    vpmadd52luq \A3, \U0, \C
    vpmadd52huq \Y3, \U0, \C
    vpmadd52luq \B3, \U1, \C
    vpmadd52huq \Z3, \U1, \C

    vpbroadcastq \C, QWORD PTR [rip + p4_q]
    vpmadd52luq \A4, \U0, \C
    vpmadd52huq \Y4, \U0, \C
    vpmadd52luq \B4, \U1, \C
    vpmadd52huq \Z4, \U1, \C

    vpbroadcastq \C, QWORD PTR [rip + p5_q]
    vpmadd52luq \A5, \U0, \C
    vpmadd52huq \Y5, \U0, \C
    vpmadd52luq \B5, \U1, \C
    vpmadd52huq \Z5, \U1, \C

    vpbroadcastq \C, QWORD PTR [rip + p6_q]
    vpmadd52luq \A6, \U0, \C
    vpmadd52huq \YH, \U0, \C
    vpmadd52luq \B6, \U1, \C
    vpmadd52huq \ZH, \U1, \C

    # p7 eager fold into loaded next limb
    vmovdqa64 \A0, ZMMWORD PTR [\BASEA + \OFFA]
    vmovdqa64 \B0, ZMMWORD PTR [\BASEB + \OFFB]
    vpbroadcastq \C, QWORD PTR [rip + p7_q]
    vpmadd52luq \A0, \U0, \C
    vpmadd52luq \B0, \U1, \C
    vpmadd52huq \T, \U0, \C
    vpsllq \T, \T, 4
    vpaddq \A0, \A0, \T
    vpmadd52huq \T, \U1, \C
    vpsllq \T, \T, 4
    vpaddq \B0, \B0, \T
.endm

blst_redc_core_k2_asm:
    # stream 0 init: zmm0..zmm7 = a0..a7
    vmovdqa64 zmm0, [rdx + 0]
    vmovdqa64 zmm1, [rdx + 64]
    vmovdqa64 zmm2, [rdx + 128]
    vmovdqa64 zmm3, [rdx + 192]
    vmovdqa64 zmm4, [rdx + 256]
    vmovdqa64 zmm5, [rdx + 320]
    vmovdqa64 zmm6, [rdx + 384]
    vmovdqa64 zmm7, [rdx + 448]

    # stream 1 init: zmm8..zmm15 = b0..b7
    vmovdqa64 zmm8,  [rcx + 0]
    vmovdqa64 zmm9,  [rcx + 64]
    vmovdqa64 zmm10, [rcx + 128]
    vmovdqa64 zmm11, [rcx + 192]
    vmovdqa64 zmm12, [rcx + 256]
    vmovdqa64 zmm13, [rcx + 320]
    vmovdqa64 zmm14, [rcx + 384]
    vmovdqa64 zmm15, [rcx + 448]

    # lazy hi queues (7 slots each stream)
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
    vpxorq zmm29, zmm29, zmm29

    # Round 0
    REDC_ROUND zmm30, zmm31, zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7, zmm16, zmm17, zmm18, zmm19, zmm20, zmm21, zmm22, rdx, 512
    REDC_ROUND zmm31, zmm30, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15, zmm23, zmm24, zmm25, zmm26, zmm27, zmm28, zmm29, rcx, 512
    # Round 1
    REDC_ROUND zmm30, zmm31, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7, zmm0, zmm17, zmm18, zmm19, zmm20, zmm21, zmm22, zmm16, rdx, 576
    REDC_ROUND zmm31, zmm30, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15, zmm8, zmm24, zmm25, zmm26, zmm27, zmm28, zmm29, zmm23, rcx, 576
    # Round 2
    REDC_ROUND zmm30, zmm31, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7, zmm0, zmm1, zmm18, zmm19, zmm20, zmm21, zmm22, zmm16, zmm17, rdx, 640
    REDC_ROUND zmm31, zmm30, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15, zmm8, zmm9, zmm25, zmm26, zmm27, zmm28, zmm29, zmm23, zmm24, rcx, 640
    # Round 3
    REDC_ROUND zmm30, zmm31, zmm3, zmm4, zmm5, zmm6, zmm7, zmm0, zmm1, zmm2, zmm19, zmm20, zmm21, zmm22, zmm16, zmm17, zmm18, rdx, 704
    REDC_ROUND zmm31, zmm30, zmm11, zmm12, zmm13, zmm14, zmm15, zmm8, zmm9, zmm10, zmm26, zmm27, zmm28, zmm29, zmm23, zmm24, zmm25, rcx, 704
    # Round 4
    REDC_ROUND zmm30, zmm31, zmm4, zmm5, zmm6, zmm7, zmm0, zmm1, zmm2, zmm3, zmm20, zmm21, zmm22, zmm16, zmm17, zmm18, zmm19, rdx, 768
    REDC_ROUND zmm31, zmm30, zmm12, zmm13, zmm14, zmm15, zmm8, zmm9, zmm10, zmm11, zmm27, zmm28, zmm29, zmm23, zmm24, zmm25, zmm26, rcx, 768
    # Round 5
    REDC_ROUND zmm30, zmm31, zmm5, zmm6, zmm7, zmm0, zmm1, zmm2, zmm3, zmm4, zmm21, zmm22, zmm16, zmm17, zmm18, zmm19, zmm20, rdx, 832
    REDC_ROUND zmm31, zmm30, zmm13, zmm14, zmm15, zmm8, zmm9, zmm10, zmm11, zmm12, zmm28, zmm29, zmm23, zmm24, zmm25, zmm26, zmm27, rcx, 832
    # Round 6
    REDC_ROUND zmm30, zmm31, zmm6, zmm7, zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm22, zmm16, zmm17, zmm18, zmm19, zmm20, zmm21, rdx, 896
    REDC_ROUND zmm31, zmm30, zmm14, zmm15, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm29, zmm23, zmm24, zmm25, zmm26, zmm27, zmm28, rcx, 896
    # Round 7
    REDC_ROUND zmm30, zmm31, zmm7, zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm16, zmm17, zmm18, zmm19, zmm20, zmm21, zmm22, rdx, 960
    REDC_ROUND zmm31, zmm30, zmm15, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm23, zmm24, zmm25, zmm26, zmm27, zmm28, zmm29, rcx, 960

    # Tail adds (queue after round 7 corresponds to y8..y14)
    vpsllq zmm30, zmm17, 4
    vpaddq zmm1, zmm1, zmm30
    vpsllq zmm30, zmm18, 4
    vpaddq zmm2, zmm2, zmm30
    vpsllq zmm30, zmm19, 4
    vpaddq zmm3, zmm3, zmm30
    vpsllq zmm30, zmm20, 4
    vpaddq zmm4, zmm4, zmm30
    vpsllq zmm30, zmm21, 4
    vpaddq zmm5, zmm5, zmm30
    vpsllq zmm30, zmm22, 4
    vpaddq zmm6, zmm6, zmm30
    vpsllq zmm30, zmm16, 4
    vpaddq zmm7, zmm7, zmm30

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
    vpsllq zmm31, zmm29, 4
    vpaddq zmm14, zmm14, zmm31
    vpsllq zmm31, zmm23, 4
    vpaddq zmm15, zmm15, zmm31

    # Final carry/mask over 8 limbs each stream
    vmovdqa64 zmm31, ZMMWORD PTR [rip + bmask_vec]
    REDC_CARRY8_X2 zmm30, zmm29, zmm31, zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15

    # Store outputs
    vmovdqa64 [rdi + 0], zmm0
    vmovdqa64 [rdi + 64], zmm1
    vmovdqa64 [rdi + 128], zmm2
    vmovdqa64 [rdi + 192], zmm3
    vmovdqa64 [rdi + 256], zmm4
    vmovdqa64 [rdi + 320], zmm5
    vmovdqa64 [rdi + 384], zmm6
    vmovdqa64 [rdi + 448], zmm7

    vmovdqa64 [rsi + 0], zmm8
    vmovdqa64 [rsi + 64], zmm9
    vmovdqa64 [rsi + 128], zmm10
    vmovdqa64 [rsi + 192], zmm11
    vmovdqa64 [rsi + 256], zmm12
    vmovdqa64 [rsi + 320], zmm13
    vmovdqa64 [rsi + 384], zmm14
    vmovdqa64 [rsi + 448], zmm15

    ret

.size blst_redc_core_k2_asm, .-blst_redc_core_k2_asm

blst_redc_core_k2_asm_v2:
    # stream 0 init: zmm0..zmm7 = a0..a7
    vmovdqa64 zmm0, [rdx + 0]
    vmovdqa64 zmm1, [rdx + 64]
    vmovdqa64 zmm2, [rdx + 128]
    vmovdqa64 zmm3, [rdx + 192]
    vmovdqa64 zmm4, [rdx + 256]
    vmovdqa64 zmm5, [rdx + 320]
    vmovdqa64 zmm6, [rdx + 384]
    vmovdqa64 zmm7, [rdx + 448]

    # stream 1 init: zmm8..zmm15 = b0..b7
    vmovdqa64 zmm8,  [rcx + 0]
    vmovdqa64 zmm9,  [rcx + 64]
    vmovdqa64 zmm10, [rcx + 128]
    vmovdqa64 zmm11, [rcx + 192]
    vmovdqa64 zmm12, [rcx + 256]
    vmovdqa64 zmm13, [rcx + 320]
    vmovdqa64 zmm14, [rcx + 384]
    vmovdqa64 zmm15, [rcx + 448]

    # lazy hi queues: 6 slots per stream
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

    # Round 0
    REDC_ROUND_BOTH6 zmm28, zmm29, zmm30, zmm31, zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15, zmm16, zmm17, zmm18, zmm19, zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27, rdx, 512, rcx, 512
    # Round 1
    REDC_ROUND_BOTH6 zmm28, zmm29, zmm30, zmm31, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7, zmm0, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15, zmm8, zmm17, zmm18, zmm19, zmm20, zmm21, zmm16, zmm23, zmm24, zmm25, zmm26, zmm27, zmm22, rdx, 576, rcx, 576
    # Round 2
    REDC_ROUND_BOTH6 zmm28, zmm29, zmm30, zmm31, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7, zmm0, zmm1, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15, zmm8, zmm9, zmm18, zmm19, zmm20, zmm21, zmm16, zmm17, zmm24, zmm25, zmm26, zmm27, zmm22, zmm23, rdx, 640, rcx, 640
    # Round 3
    REDC_ROUND_BOTH6 zmm28, zmm29, zmm30, zmm31, zmm3, zmm4, zmm5, zmm6, zmm7, zmm0, zmm1, zmm2, zmm11, zmm12, zmm13, zmm14, zmm15, zmm8, zmm9, zmm10, zmm19, zmm20, zmm21, zmm16, zmm17, zmm18, zmm25, zmm26, zmm27, zmm22, zmm23, zmm24, rdx, 704, rcx, 704
    # Round 4
    REDC_ROUND_BOTH6 zmm28, zmm29, zmm30, zmm31, zmm4, zmm5, zmm6, zmm7, zmm0, zmm1, zmm2, zmm3, zmm12, zmm13, zmm14, zmm15, zmm8, zmm9, zmm10, zmm11, zmm20, zmm21, zmm16, zmm17, zmm18, zmm19, zmm26, zmm27, zmm22, zmm23, zmm24, zmm25, rdx, 768, rcx, 768
    # Round 5
    REDC_ROUND_BOTH6 zmm28, zmm29, zmm30, zmm31, zmm5, zmm6, zmm7, zmm0, zmm1, zmm2, zmm3, zmm4, zmm13, zmm14, zmm15, zmm8, zmm9, zmm10, zmm11, zmm12, zmm21, zmm16, zmm17, zmm18, zmm19, zmm20, zmm27, zmm22, zmm23, zmm24, zmm25, zmm26, rdx, 832, rcx, 832
    # Round 6
    REDC_ROUND_BOTH6 zmm28, zmm29, zmm30, zmm31, zmm6, zmm7, zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm14, zmm15, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm16, zmm17, zmm18, zmm19, zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27, rdx, 896, rcx, 896
    # Round 7
    REDC_ROUND_BOTH6 zmm28, zmm29, zmm30, zmm31, zmm7, zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm15, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm17, zmm18, zmm19, zmm20, zmm21, zmm16, zmm23, zmm24, zmm25, zmm26, zmm27, zmm22, rdx, 960, rcx, 960

    # Tail adds for 6-slot queue (post round-7 mapping: q0..q5 = z18,z19,z20,z21,z16,z17)
    vpsllq zmm31, zmm18, 4
    vpaddq zmm1, zmm1, zmm31
    vpsllq zmm31, zmm19, 4
    vpaddq zmm2, zmm2, zmm31
    vpsllq zmm31, zmm20, 4
    vpaddq zmm3, zmm3, zmm31
    vpsllq zmm31, zmm21, 4
    vpaddq zmm4, zmm4, zmm31
    vpsllq zmm31, zmm16, 4
    vpaddq zmm5, zmm5, zmm31
    vpsllq zmm31, zmm17, 4
    vpaddq zmm6, zmm6, zmm31

    # stream1 post round-7 mapping: q0..q5 = z24,z25,z26,z27,z22,z23
    vpsllq zmm31, zmm24, 4
    vpaddq zmm9, zmm9, zmm31
    vpsllq zmm31, zmm25, 4
    vpaddq zmm10, zmm10, zmm31
    vpsllq zmm31, zmm26, 4
    vpaddq zmm11, zmm11, zmm31
    vpsllq zmm31, zmm27, 4
    vpaddq zmm12, zmm12, zmm31
    vpsllq zmm31, zmm22, 4
    vpaddq zmm13, zmm13, zmm31
    vpsllq zmm31, zmm23, 4
    vpaddq zmm14, zmm14, zmm31

    # Final carry/mask over 8 limbs each stream
    vmovdqa64 zmm31, ZMMWORD PTR [rip + bmask_vec]
    REDC_CARRY8 zmm30, zmm31, zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7
    REDC_CARRY8 zmm30, zmm31, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15

    # Store outputs
    vmovdqa64 [rdi + 0], zmm0
    vmovdqa64 [rdi + 64], zmm1
    vmovdqa64 [rdi + 128], zmm2
    vmovdqa64 [rdi + 192], zmm3
    vmovdqa64 [rdi + 256], zmm4
    vmovdqa64 [rdi + 320], zmm5
    vmovdqa64 [rdi + 384], zmm6
    vmovdqa64 [rdi + 448], zmm7

    vmovdqa64 [rsi + 0], zmm8
    vmovdqa64 [rsi + 64], zmm9
    vmovdqa64 [rsi + 128], zmm10
    vmovdqa64 [rsi + 192], zmm11
    vmovdqa64 [rsi + 256], zmm12
    vmovdqa64 [rsi + 320], zmm13
    vmovdqa64 [rsi + 384], zmm14
    vmovdqa64 [rsi + 448], zmm15

    ret

.size blst_redc_core_k2_asm_v2, .-blst_redc_core_k2_asm_v2

.section .rodata
.p2align 6

montw_vec:
    .quad 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD
    .quad 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD, 0xFFFCFFFCFFFD
bmask_vec:
    .quad 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF
    .quad 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF
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

montw_q:
    .quad 0xFFFCFFFCFFFD
bmask_q:
    .quad 0xFFFFFFFFFFFF
p0_q:
    .quad 0xFFFFFFFFAAAB
p1_q:
    .quad 0xB153FFFFB9FE
p2_q:
    .quad 0xF6241EABFFFE
p3_q:
    .quad 0x6730D2A0F6B0
p4_q:
    .quad 0x4B84F38512BF
p5_q:
    .quad 0x434BACD76477
p6_q:
    .quad 0xE69A4B1BA7B6
p7_q:
    .quad 0x1A0111EA397F

.att_syntax prefix

#endif
