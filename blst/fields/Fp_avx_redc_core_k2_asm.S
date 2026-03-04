.intel_syntax noprefix
.text

# blst_fp8x1_redc_core_k2
#
# Purpose:
# - Montgomery reduction core for two independent UFpPack8x1 inputs in parallel.
# - Produces two LFpPack8x1 outputs (lazy reduced, carry-propagated, no final sub-p).
#
# Why this schedule:
# - This is the promoted "v3" schedule from benchmark tuning.
# - Keeps the arithmetic model identical to C++ redc_core (6-slot lazy-high queue).
# - Avoids stack spills while interleaving two streams.
# - Breaks a key temp dependency by using C as stream-1 high-fold temp in p0/p7 steps.
#
# Known bottlenecks (current state):
# - Backend bound from carry/shift/add dependency chains.
# - IFMA count is fixed by algorithm; gains come from dependency scheduling, not fewer IFMAs.
# - Further significant gains likely require deeper arithmetic/schedule changes or mul+redc fusion.
#
# Calling convention (SysV x86_64):
# - rdi: LFpPack8x1* out0
# - rsi: LFpPack8x1* out1
# - rdx: const UFpPack8x1* in0
# - rcx: const UFpPack8x1* in1

.globl blst_fp8x1_redc_core_k2
.type blst_fp8x1_redc_core_k2, @function

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

.macro REDC_ROUND_BOTH6_CTEMP U0, U1, C, T, A0, A1, A2, A3, A4, A5, A6, A7, B0, B1, B2, B3, B4, B5, B6, B7, YH, Y1, Y2, Y3, Y4, Y5, ZH, Z1, Z2, Z3, Z4, Z5, BASEA, OFFA, BASEB, OFFB
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

    # add p0 high directly to next limb; stream1 uses C to avoid temp reuse chain
    vpmadd52huq \T, \U0, \C
    vpsllq \T, \T, 4
    vpaddq \A1, \A1, \T
    vpmadd52huq \C, \U1, \C
    vpsllq \C, \C, 4
    vpaddq \B1, \B1, \C

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
    vpmadd52huq \C, \U1, \C
    vpsllq \C, \C, 4
    vpaddq \B0, \B0, \C
.endm

blst_fp8x1_redc_core_k2:
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
    REDC_ROUND_BOTH6_CTEMP zmm28, zmm29, zmm30, zmm31, zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15, zmm16, zmm17, zmm18, zmm19, zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27, rdx, 512, rcx, 512
    # Round 1
    REDC_ROUND_BOTH6_CTEMP zmm28, zmm29, zmm30, zmm31, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7, zmm0, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15, zmm8, zmm17, zmm18, zmm19, zmm20, zmm21, zmm16, zmm23, zmm24, zmm25, zmm26, zmm27, zmm22, rdx, 576, rcx, 576
    # Round 2
    REDC_ROUND_BOTH6_CTEMP zmm28, zmm29, zmm30, zmm31, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7, zmm0, zmm1, zmm10, zmm11, zmm12, zmm13, zmm14, zmm15, zmm8, zmm9, zmm18, zmm19, zmm20, zmm21, zmm16, zmm17, zmm24, zmm25, zmm26, zmm27, zmm22, zmm23, rdx, 640, rcx, 640
    # Round 3
    REDC_ROUND_BOTH6_CTEMP zmm28, zmm29, zmm30, zmm31, zmm3, zmm4, zmm5, zmm6, zmm7, zmm0, zmm1, zmm2, zmm11, zmm12, zmm13, zmm14, zmm15, zmm8, zmm9, zmm10, zmm19, zmm20, zmm21, zmm16, zmm17, zmm18, zmm25, zmm26, zmm27, zmm22, zmm23, zmm24, rdx, 704, rcx, 704
    # Round 4
    REDC_ROUND_BOTH6_CTEMP zmm28, zmm29, zmm30, zmm31, zmm4, zmm5, zmm6, zmm7, zmm0, zmm1, zmm2, zmm3, zmm12, zmm13, zmm14, zmm15, zmm8, zmm9, zmm10, zmm11, zmm20, zmm21, zmm16, zmm17, zmm18, zmm19, zmm26, zmm27, zmm22, zmm23, zmm24, zmm25, rdx, 768, rcx, 768
    # Round 5
    REDC_ROUND_BOTH6_CTEMP zmm28, zmm29, zmm30, zmm31, zmm5, zmm6, zmm7, zmm0, zmm1, zmm2, zmm3, zmm4, zmm13, zmm14, zmm15, zmm8, zmm9, zmm10, zmm11, zmm12, zmm21, zmm16, zmm17, zmm18, zmm19, zmm20, zmm27, zmm22, zmm23, zmm24, zmm25, zmm26, rdx, 832, rcx, 832
    # Round 6
    REDC_ROUND_BOTH6_CTEMP zmm28, zmm29, zmm30, zmm31, zmm6, zmm7, zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm14, zmm15, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm16, zmm17, zmm18, zmm19, zmm20, zmm21, zmm22, zmm23, zmm24, zmm25, zmm26, zmm27, rdx, 896, rcx, 896
    # Round 7
    REDC_ROUND_BOTH6_CTEMP zmm28, zmm29, zmm30, zmm31, zmm7, zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm15, zmm8, zmm9, zmm10, zmm11, zmm12, zmm13, zmm14, zmm17, zmm18, zmm19, zmm20, zmm21, zmm16, zmm23, zmm24, zmm25, zmm26, zmm27, zmm22, rdx, 960, rcx, 960

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

.size blst_fp8x1_redc_core_k2, .-blst_fp8x1_redc_core_k2

.section .rodata
.p2align 6

bmask_vec:
    .quad 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF
    .quad 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFF

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
