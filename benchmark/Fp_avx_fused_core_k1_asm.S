.intel_syntax noprefix
.text

# Exact instruction-stream clone of the C++ fused baseline core body.
# Source of truth: objdump of fp8x1_mul_redc_core_k1_fused_baseline in
# Fp_avx_fused_baseline_bench.cpp.o (Linux AVX-512 build).

.extern _ZN4blst3P48E

.globl blst_fp8x1_mul_redc_core_k1_asm
.type blst_fp8x1_mul_redc_core_k1_asm, @function
blst_fp8x1_mul_redc_core_k1_asm:
    # SysV x86_64:
    #   rdi = LFpPack8x1* out
    #   rsi = const FpPack8x1* in0
    #   rdx = const FpPack8x1* in1
    #
    # High-level flow (kept identical to C++ object clone):
    # 1) Fused multiply body (IFMA52 + Karatsuba-style cross-term correction).
    # 2) Montgomery reduction (8 fully unrolled REDC rounds).
    # 3) Final carry propagation + radix-48 masking, then store 8 limbs.
    push   rbp
    vpxor  xmm0,xmm0,xmm0
    vmovdqa64 zmm5,zmm0
    mov    rbp,rsp
    and    rsp,0xffffffffffffffc0
    sub    rsp,0x288
    # Stack scratch (compiler spill layout preserved intentionally):
    # +0x248 +0x208 +0x1c8 +0x188 +0x148 +0x108 +0x0c8 +0x088 +0x048 +0x008
    # -0x038 -0x078
    vmovdqa64 zmm31,zmm0
    vmovdqa64 zmm7,zmm0
    vmovdqa64 zmm9,zmm0
    vmovdqa64 zmm21,zmm0
    vmovdqa64 zmm23,zmm0
    # Phase A: fused mul accumulation.
    # The compiler schedule interleaves:
    # - low/high IFMA products for radix-48 limbs
    # - +4-bit shifts for hi halves (align to radix-48)
    # - Karatsuba correction terms between lower and upper 4-limb blocks
    # with frequent stack staging to keep IFMA chains fed.
    vmovdqa64 zmm1,ZMMWORD PTR [rdx+0xc0]
    vmovdqa64 zmm19,ZMMWORD PTR [rsi+0x40]
    vmovdqa64 zmm10,ZMMWORD PTR [rsi+0x80]
    vpmadd52huq zmm5,zmm19,zmm1
    vmovdqa64 zmm22,ZMMWORD PTR [rdx+0x80]
    vmovdqa64 zmm24,ZMMWORD PTR [rdx+0x40]
    vmovdqa64 zmm8,ZMMWORD PTR [rsi+0xc0]
    vmovdqa64 zmm30,ZMMWORD PTR [rsi+0x100]
    vpmadd52huq zmm5,zmm10,zmm22
    vmovdqa64 zmm3,ZMMWORD PTR [rdx+0x100]
    vmovdqa64 zmm4,ZMMWORD PTR [rdx+0x140]
    vpmadd52huq zmm31,zmm8,zmm1
    vmovdqa64 zmm29,ZMMWORD PTR [rsi+0x140]
    vpmadd52huq zmm5,zmm8,zmm24
    mov    rax,QWORD PTR _ZN4blst3P48E@GOTPCREL[rip]
    vmovdqa64 zmm28,ZMMWORD PTR [rsi+0x180]
    vpsllq zmm31,zmm31,0x4
    vmovdqa64 ZMMWORD PTR [rsp+0x1c8],zmm31
    vpsllq zmm5,zmm5,0x4
    vpmadd52luq zmm5,zmm10,zmm1
    vmovdqa64 zmm31,zmm0
    vpmadd52luq zmm31,zmm30,zmm3
    vpbroadcastq zmm18,QWORD PTR [rax]
    vpbroadcastq zmm17,QWORD PTR [rax+0x8]
    vpmadd52luq zmm5,zmm8,zmm22
    vpbroadcastq zmm16,QWORD PTR [rax+0x10]
    vpbroadcastq zmm15,QWORD PTR [rax+0x18]
    vpbroadcastq zmm14,QWORD PTR [rax+0x20]
    vpbroadcastq zmm13,QWORD PTR [rax+0x28]
    vmovdqa64 ZMMWORD PTR [rsp+0x248],zmm5
    vmovdqa64 zmm5,zmm0
    vpmadd52huq zmm5,zmm10,zmm1
    vpbroadcastq zmm12,QWORD PTR [rax+0x30]
    vpbroadcastq zmm11,QWORD PTR [rax+0x38]
    vmovdqa64 zmm27,ZMMWORD PTR [rsi+0x1c0]
    vmovdqa64 zmm26,ZMMWORD PTR [rdx]
    vpmadd52huq zmm5,zmm8,zmm22
    vmovdqa64 zmm2,ZMMWORD PTR [rdx+0x180]
    vmovdqa64 zmm6,ZMMWORD PTR [rdx+0x1c0]
    vmovdqa64 zmm20,ZMMWORD PTR [rsi]
    vmovdqa64 zmm25,zmm0
    vpsllq zmm5,zmm5,0x4
    vpmadd52luq zmm5,zmm8,zmm1
    vpmadd52huq zmm7,zmm20,zmm1
    vpmadd52huq zmm9,zmm20,zmm22
    vpmadd52huq zmm21,zmm20,zmm24
    vpmadd52huq zmm23,zmm20,zmm26
    vmovdqa64 ZMMWORD PTR [rsp+0x208],zmm5
    vmovdqa64 zmm5,zmm0
    vpmadd52huq zmm5,zmm30,zmm3
    vmovdqa64 ZMMWORD PTR [rsp+0x188],zmm31
    vpmadd52huq zmm7,zmm19,zmm22
    vpmadd52huq zmm9,zmm19,zmm24
    vpmadd52huq zmm21,zmm19,zmm26
    vpsllq zmm5,zmm5,0x4
    vpmadd52luq zmm5,zmm30,zmm4
    vpmadd52huq zmm7,zmm10,zmm24
    vpmadd52huq zmm9,zmm10,zmm26
    vpsllq zmm21,zmm21,0x4
    vpmadd52luq zmm21,zmm20,zmm22
    vpmadd52luq zmm5,zmm29,zmm3
    vpmadd52huq zmm7,zmm8,zmm26
    vpsllq zmm9,zmm9,0x4
    vpmadd52luq zmm9,zmm20,zmm1
    vpsllq zmm23,zmm23,0x4
    vmovdqa64 ZMMWORD PTR [rsp+0x148],zmm5
    vmovdqa64 zmm5,zmm0
    vpmadd52huq zmm5,zmm30,zmm4
    vpsllq zmm7,zmm7,0x4
    vpmadd52luq zmm7,zmm19,zmm1
    vpmadd52luq zmm9,zmm19,zmm22
    vpmadd52luq zmm25,zmm20,zmm26
    vpmadd52huq zmm5,zmm29,zmm3
    vpmadd52luq zmm23,zmm20,zmm24
    vpmadd52luq zmm7,zmm10,zmm22
    vpaddq zmm20,zmm20,zmm30
    vpaddq zmm22,zmm22,zmm2
    vpsllq zmm5,zmm5,0x4
    vpmadd52luq zmm5,zmm30,zmm2
    vpmadd52luq zmm7,zmm8,zmm24
    vpmadd52luq zmm21,zmm19,zmm24
    vpmadd52luq zmm9,zmm10,zmm24
    vpmadd52luq zmm23,zmm19,zmm26
    vpmadd52luq zmm5,zmm29,zmm4
    vpaddq zmm19,zmm19,zmm29
    vpmadd52luq zmm21,zmm10,zmm26
    vpmadd52luq zmm9,zmm8,zmm26
    vpaddq zmm10,zmm10,zmm28
    vmovdqa64 zmm31,zmm5
    vmovdqa64 zmm5,zmm0
    vpmadd52huq zmm5,zmm30,zmm2
    vpmadd52luq zmm31,zmm28,zmm3
    vpaddq zmm8,zmm8,zmm27
    vpmadd52huq zmm5,zmm29,zmm4
    vmovdqa64 ZMMWORD PTR [rsp+0x108],zmm31
    vmovdqa64 zmm31,zmm0
    vpmadd52huq zmm31,zmm29,zmm6
    vpmadd52huq zmm5,zmm28,zmm3
    vpsllq zmm5,zmm5,0x4
    vpmadd52luq zmm5,zmm30,zmm6
    vpmadd52luq zmm5,zmm29,zmm2
    vpmadd52luq zmm5,zmm28,zmm4
    vpmadd52luq zmm5,zmm27,zmm3
    vmovdqa64 ZMMWORD PTR [rsp+0xc8],zmm5
    vmovdqa64 zmm5,zmm0
    vpmadd52huq zmm5,zmm30,zmm6
    vmovdqa64 zmm30,zmm25
    vpmadd52huq zmm5,zmm29,zmm2
    vpmadd52huq zmm5,zmm28,zmm4
    vpmadd52huq zmm5,zmm27,zmm3
    vpaddq zmm3,zmm26,zmm3
    # Begin REDC setup constants:
    # zmm26 = montgomery constant (-p^{-1} mod 2^48).
    vpbroadcastq zmm26,QWORD PTR [rip+montw_q]
    vpsllq zmm5,zmm5,0x4
    vpmadd52luq zmm5,zmm29,zmm6
    vpmadd52luq zmm5,zmm28,zmm2
    vpmadd52luq zmm5,zmm27,zmm4
    vmovdqa64 ZMMWORD PTR [rsp+0x48],zmm5
    vmovdqa64 zmm5,zmm31
    vpmadd52huq zmm5,zmm28,zmm2
    vmovdqa64 zmm31,zmm0
    vpmadd52huq zmm31,zmm28,zmm6
    vpmadd52huq zmm5,zmm27,zmm4
    vpaddq zmm4,zmm24,zmm4
    vpaddq zmm24,zmm1,zmm6
    vmovdqa64 zmm1,zmm0
    vpmadd52huq zmm1,zmm20,zmm22
    vpsllq zmm5,zmm5,0x4
    vpmadd52luq zmm5,zmm28,zmm6
    vpmadd52huq zmm1,zmm19,zmm4
    vpmadd52luq zmm5,zmm27,zmm2
    vpmadd52huq zmm1,zmm10,zmm3
    vmovdqa64 ZMMWORD PTR [rsp+0x8],zmm5
    vmovdqa64 zmm5,zmm31
    vpmadd52huq zmm5,zmm27,zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52huq zmm2,zmm20,zmm3
    vpsllq zmm1,zmm1,0x4
    vpmadd52luq zmm1,zmm20,zmm24
    vpsllq zmm5,zmm5,0x4
    vmovdqa64 zmm31,zmm5
    vmovdqa64 zmm5,zmm0
    vpmadd52luq zmm31,zmm27,zmm6
    vpmadd52huq zmm5,zmm27,zmm6
    vmovdqa64 zmm27,zmm0
    vpmadd52huq zmm27,zmm20,zmm4
    vpmadd52luq zmm1,zmm19,zmm22
    vmovdqa64 zmm6,zmm0
    vpsllq zmm5,zmm5,0x4
    vmovdqa64 ZMMWORD PTR [rsp+0x88],zmm5
    vpmadd52huq zmm27,zmm19,zmm3
    vpsllq zmm5,zmm2,0x4
    vmovdqa64 zmm2,zmm0
    vpmadd52huq zmm2,zmm20,zmm24
    vpmadd52luq zmm5,zmm20,zmm4
    vpsllq zmm27,zmm27,0x4
    vpmadd52luq zmm27,zmm20,zmm22
    vpmadd52luq zmm1,zmm10,zmm4
    vpmadd52huq zmm2,zmm19,zmm22
    vpmadd52luq zmm6,zmm20,zmm3
    vpmadd52luq zmm5,zmm19,zmm3
    vpmadd52luq zmm27,zmm19,zmm4
    vpmadd52luq zmm1,zmm8,zmm3
    vpmadd52huq zmm2,zmm10,zmm4
    vmovdqa64 zmm20,zmm21
    vpaddq zmm21,zmm21,ZMMWORD PTR [rsp+0x108]
    vpmadd52luq zmm27,zmm10,zmm3
    vpaddq zmm1,zmm1,ZMMWORD PTR [rsp+0x1c8]
    vpmadd52huq zmm2,zmm8,zmm3
    vmovdqa64 zmm3,zmm0
    vpmadd52huq zmm3,zmm19,zmm24
    vpaddq zmm27,zmm27,ZMMWORD PTR [rsp+0x208]
    vpaddq zmm6,zmm6,zmm7
    vpsllq zmm2,zmm2,0x4
    vpmadd52luq zmm2,zmm19,zmm24
    vpmadd52huq zmm3,zmm10,zmm22
    vpaddq zmm5,zmm5,ZMMWORD PTR [rsp+0x248]
    vpmadd52luq zmm2,zmm10,zmm22
    vpmadd52huq zmm3,zmm8,zmm4
    vpmadd52luq zmm2,zmm8,zmm4
    vpsllq zmm4,zmm3,0x4
    vmovdqa64 zmm3,zmm0
    vpmadd52huq zmm3,zmm10,zmm24
    vpmadd52luq zmm4,zmm10,zmm24
    vmovdqa64 zmm10,zmm0
    vpmadd52huq zmm10,zmm8,zmm24
    vpaddq zmm2,zmm2,ZMMWORD PTR [rsp+0x188]
    vpmadd52huq zmm3,zmm8,zmm22
    vpmadd52luq zmm4,zmm8,zmm22
    vpsllq zmm10,zmm10,0x4
    vmovdqa64 ZMMWORD PTR [rsp-0x78],zmm10
    vpsllq zmm3,zmm3,0x4
    vpmadd52luq zmm3,zmm8,zmm24
    vmovdqa64 zmm8,zmm0
    # Phase B: 8 unrolled Montgomery REDC rounds.
    # Per round:
    #   u = (pivot * montw) & BMASK
    #   accum += u * p[0..7] using IFMA lo/hi
    #   fold carry from p0 low into next limb
    #
    # REDC round 0 (pivot in zmm25).
    vpmadd52luq zmm8,zmm25,zmm26
    # zmm10 = BMASK in all lanes.
    vmovdqa32 zmm10,ZMMWORD PTR [rip+bmask_vec]
    vpaddq zmm25,zmm25,ZMMWORD PTR [rsp+0x188]
    vmovdqa64 ZMMWORD PTR [rsp-0x38],zmm3
    vpsubq zmm3,zmm27,zmm21
    vpandd zmm8,zmm8,zmm10
    vpmadd52luq zmm20,zmm8,zmm16
    vpsubq zmm6,zmm6,zmm25
    vmovdqa64 zmm24,zmm23
    vpaddq zmm23,zmm23,ZMMWORD PTR [rsp+0x148]
    vpmadd52luq zmm30,zmm8,zmm18
    vmovdqa64 zmm22,zmm20
    vmovdqa64 zmm20,zmm9
    vpaddq zmm9,zmm9,ZMMWORD PTR [rsp+0xc8]
    vpsubq zmm5,zmm5,zmm23
    vpsubq zmm1,zmm1,zmm9
    vmovdqa64 zmm9,zmm0
    vpmadd52huq zmm9,zmm8,zmm18
    vpmadd52luq zmm24,zmm8,zmm17
    vpsraq zmm30,zmm30,0x30
    vpmadd52luq zmm20,zmm8,zmm15
    vpmadd52luq zmm6,zmm8,zmm14
    vmovdqa64 zmm29,zmm9
    vmovdqa64 zmm9,zmm0
    vpmadd52huq zmm9,zmm8,zmm17
    vpsllq zmm29,zmm29,0x4
    vpaddq zmm24,zmm24,zmm30
    vpmadd52luq zmm5,zmm8,zmm13
    vpmadd52luq zmm3,zmm8,zmm12
    vmovdqa64 zmm28,zmm9
    vmovdqa64 zmm9,zmm0
    vpmadd52huq zmm9,zmm8,zmm16
    vpmadd52luq zmm1,zmm8,zmm11
    vpaddq zmm24,zmm24,zmm29
    vmovdqa64 zmm29,ZMMWORD PTR [rsp+0x48]
    vmovdqa64 zmm30,ZMMWORD PTR [rsp+0x8]
    vmovdqa64 zmm27,zmm9
    vmovdqa64 zmm9,zmm0
    vpmadd52huq zmm9,zmm8,zmm15
    vpaddq zmm7,zmm7,zmm29
    vpsubq zmm2,zmm2,zmm7
    vmovdqa64 zmm7,zmm0
    vpaddq zmm4,zmm4,ZMMWORD PTR [rsp+0x148]
    vmovdqa64 zmm25,zmm9
    vmovdqa64 zmm9,zmm0
    vpmadd52huq zmm9,zmm8,zmm14
    vmovdqa64 zmm23,zmm9
    vmovdqa64 zmm9,zmm0
    vpmadd52huq zmm9,zmm8,zmm13
    vmovdqa64 zmm21,zmm9
    vmovdqa64 zmm9,zmm0
    vpmadd52huq zmm9,zmm8,zmm12
    vmovdqa64 zmm19,zmm9
    vmovdqa64 zmm9,zmm0
    vpmadd52huq zmm9,zmm8,zmm11
    vmovdqa64 zmm8,zmm0
    # REDC round 1 (pivot in zmm24).
    vpmadd52luq zmm8,zmm24,zmm26
    vpandd zmm8,zmm8,zmm10
    vpmadd52luq zmm24,zmm8,zmm18
    vpmadd52luq zmm22,zmm8,zmm17
    vpmadd52huq zmm28,zmm8,zmm18
    vpmadd52luq zmm20,zmm8,zmm16
    vpmadd52luq zmm6,zmm8,zmm15
    vpsraq zmm24,zmm24,0x30
    vpaddq zmm22,zmm22,zmm24
    vpsllq zmm28,zmm28,0x4
    vpmadd52luq zmm5,zmm8,zmm14
    vpmadd52luq zmm3,zmm8,zmm13
    vpmadd52luq zmm1,zmm8,zmm12
    vpmadd52luq zmm2,zmm8,zmm11
    vpmadd52huq zmm27,zmm8,zmm17
    vpmadd52huq zmm25,zmm8,zmm16
    vpmadd52huq zmm23,zmm8,zmm15
    vpmadd52huq zmm21,zmm8,zmm14
    vpmadd52huq zmm19,zmm8,zmm13
    vpmadd52huq zmm9,zmm8,zmm12
    vpmadd52huq zmm7,zmm8,zmm11
    vpaddq zmm22,zmm22,zmm28
    vmovdqa64 zmm8,zmm0
    # REDC round 2 (pivot in zmm22).
    vpmadd52luq zmm8,zmm22,zmm26
    vmovdqa64 zmm28,ZMMWORD PTR [rsp+0x248]
    vpaddq zmm24,zmm28,zmm30
    vpsubq zmm4,zmm4,zmm24
    vpandd zmm8,zmm8,zmm10
    vpmadd52luq zmm22,zmm8,zmm18
    vpmadd52luq zmm20,zmm8,zmm17
    vpmadd52huq zmm27,zmm8,zmm18
    vmovdqa64 zmm24,zmm0
    vpmadd52luq zmm6,zmm8,zmm16
    vpsraq zmm22,zmm22,0x30
    vpaddq zmm20,zmm20,zmm22
    vpsllq zmm27,zmm27,0x4
    vpmadd52luq zmm5,zmm8,zmm15
    vpmadd52luq zmm3,zmm8,zmm14
    vpmadd52luq zmm1,zmm8,zmm13
    vpmadd52luq zmm2,zmm8,zmm12
    vpmadd52luq zmm4,zmm8,zmm11
    vpmadd52huq zmm25,zmm8,zmm17
    vpmadd52huq zmm23,zmm8,zmm16
    vpmadd52huq zmm21,zmm8,zmm15
    vpmadd52huq zmm19,zmm8,zmm14
    vpmadd52huq zmm9,zmm8,zmm13
    vpmadd52huq zmm7,zmm8,zmm12
    vpmadd52huq zmm24,zmm8,zmm11
    vpaddq zmm20,zmm20,zmm27
    vmovdqa64 zmm8,zmm0
    # REDC round 3 (pivot in zmm20).
    vpmadd52luq zmm8,zmm20,zmm26
    vmovdqa64 zmm27,ZMMWORD PTR [rsp+0x108]
    vmovdqa64 zmm28,ZMMWORD PTR [rsp+0xc8]
    vpaddq zmm22,zmm27,ZMMWORD PTR [rsp-0x38]
    vpaddq zmm27,zmm31,ZMMWORD PTR [rsp+0x208]
    vpandd zmm8,zmm8,zmm10
    vpmadd52luq zmm20,zmm8,zmm18
    vpmadd52luq zmm6,zmm8,zmm17
    vpmadd52huq zmm25,zmm8,zmm18
    vpsubq zmm22,zmm22,zmm27
    vmovdqa64 zmm27,zmm0
    vpsraq zmm20,zmm20,0x30
    vpaddq zmm6,zmm6,zmm20
    vpsllq zmm25,zmm25,0x4
    vpmadd52luq zmm5,zmm8,zmm16
    vpmadd52luq zmm3,zmm8,zmm15
    vpmadd52luq zmm1,zmm8,zmm14
    vpmadd52luq zmm2,zmm8,zmm13
    vpmadd52luq zmm4,zmm8,zmm12
    vpmadd52luq zmm22,zmm8,zmm11
    vpmadd52huq zmm23,zmm8,zmm17
    vpmadd52huq zmm21,zmm8,zmm16
    vpmadd52huq zmm19,zmm8,zmm15
    vpmadd52huq zmm9,zmm8,zmm14
    vpmadd52huq zmm7,zmm8,zmm13
    vpmadd52huq zmm24,zmm8,zmm12
    vpmadd52huq zmm27,zmm8,zmm11
    vpaddq zmm6,zmm6,zmm25
    vmovdqa64 zmm8,zmm0
    # REDC round 4 (pivot in zmm6).
    vpmadd52luq zmm8,zmm6,zmm26
    vpaddq zmm20,zmm28,ZMMWORD PTR [rsp-0x78]
    vmovdqa64 zmm28,ZMMWORD PTR [rsp+0x1c8]
    vpaddq zmm25,zmm28,ZMMWORD PTR [rsp+0x88]
    vpandd zmm8,zmm8,zmm10
    vpmadd52luq zmm6,zmm8,zmm18
    vpmadd52luq zmm5,zmm8,zmm17
    vpmadd52huq zmm23,zmm8,zmm18
    vpsubq zmm20,zmm20,zmm25
    vmovdqa64 zmm28,zmm0
    vpsraq zmm6,zmm6,0x30
    vpaddq zmm5,zmm5,zmm6
    vpsllq zmm23,zmm23,0x4
    vpmadd52luq zmm3,zmm8,zmm16
    vpmadd52luq zmm1,zmm8,zmm15
    vpmadd52luq zmm2,zmm8,zmm14
    vpmadd52luq zmm4,zmm8,zmm13
    vpmadd52luq zmm22,zmm8,zmm12
    vpmadd52luq zmm20,zmm8,zmm11
    vpmadd52huq zmm21,zmm8,zmm17
    vpmadd52huq zmm19,zmm8,zmm16
    vpmadd52huq zmm9,zmm8,zmm15
    vpmadd52huq zmm7,zmm8,zmm14
    vpmadd52huq zmm24,zmm8,zmm13
    vpmadd52huq zmm27,zmm8,zmm12
    vpmadd52huq zmm28,zmm8,zmm11
    vpaddq zmm5,zmm5,zmm23
    vmovdqa64 zmm8,zmm0
    # REDC round 5 (pivot in zmm5).
    vpmadd52luq zmm8,zmm5,zmm26
    vmovdqa64 zmm6,zmm29
    vmovdqa64 zmm25,zmm28
    vmovdqa64 zmm23,zmm0
    vpandd zmm8,zmm8,zmm10
    vpmadd52luq zmm5,zmm8,zmm18
    vpmadd52luq zmm3,zmm8,zmm17
    vpmadd52huq zmm21,zmm8,zmm18
    vpmadd52luq zmm1,zmm8,zmm16
    vpmadd52huq zmm19,zmm8,zmm17
    vpsraq zmm5,zmm5,0x30
    vpaddq zmm3,zmm3,zmm5
    vpsllq zmm21,zmm21,0x4
    vpaddq zmm3,zmm3,zmm21
    vmovdqa64 zmm5,zmm0
    # REDC round 6 (pivot in zmm3).
    vpmadd52luq zmm5,zmm3,zmm26
    vpmadd52luq zmm2,zmm8,zmm15
    vpmadd52huq zmm9,zmm8,zmm16
    vpmadd52huq zmm7,zmm8,zmm15
    vpmadd52luq zmm4,zmm8,zmm14
    vpandd zmm5,zmm5,zmm10
    vpmadd52luq zmm3,zmm5,zmm18
    vpmadd52luq zmm1,zmm5,zmm17
    vpmadd52huq zmm19,zmm5,zmm18
    vpmadd52huq zmm24,zmm8,zmm14
    vpmadd52luq zmm22,zmm8,zmm13
    vpsraq zmm3,zmm3,0x30
    vpaddq zmm1,zmm1,zmm3
    vpsllq zmm19,zmm19,0x4
    vpaddq zmm1,zmm1,zmm19
    vmovdqa64 zmm3,zmm0
    # REDC round 7 (pivot in zmm1).
    vpmadd52luq zmm3,zmm1,zmm26
    vpmadd52huq zmm27,zmm8,zmm13
    vpmadd52luq zmm2,zmm5,zmm16
    vpmadd52huq zmm9,zmm5,zmm17
    vpmadd52luq zmm20,zmm8,zmm12
    vpmadd52huq zmm7,zmm5,zmm16
    vpmadd52huq zmm25,zmm8,zmm12
    vpmadd52luq zmm4,zmm5,zmm15
    vpandd zmm3,zmm3,zmm10
    vpmadd52luq zmm6,zmm8,zmm11
    vpmadd52huq zmm24,zmm5,zmm15
    vpmadd52huq zmm23,zmm8,zmm11
    vpmadd52luq zmm1,zmm3,zmm18
    vpmadd52luq zmm22,zmm5,zmm14
    vpmadd52luq zmm2,zmm3,zmm17
    vpmadd52huq zmm9,zmm3,zmm18
    vpmadd52huq zmm27,zmm5,zmm14
    vpmadd52huq zmm7,zmm3,zmm17
    vpmadd52luq zmm20,zmm5,zmm13
    vpmadd52luq zmm4,zmm3,zmm16
    vpmadd52huq zmm25,zmm5,zmm13
    vpmadd52huq zmm24,zmm3,zmm16
    vpmadd52luq zmm6,zmm5,zmm12
    vpmadd52luq zmm22,zmm3,zmm15
    vpmadd52huq zmm23,zmm5,zmm12
    vmovdqa64 zmm8,zmm30
    vpsraq zmm1,zmm1,0x30
    vpmadd52huq zmm27,zmm3,zmm15
    vpmadd52luq zmm8,zmm5,zmm11
    vmovdqa64 zmm21,zmm0
    vpsllq zmm9,zmm9,0x4
    vpaddq zmm2,zmm2,zmm1
    vpmadd52luq zmm20,zmm3,zmm14
    vpaddq zmm2,zmm2,zmm9
    vpmadd52huq zmm21,zmm5,zmm11
    vpsllq zmm7,zmm7,0x4
    vpmadd52huq zmm25,zmm3,zmm14
    vpsraq zmm1,zmm2,0x30
    vpaddq zmm4,zmm4,zmm7
    vpmadd52luq zmm6,zmm3,zmm13
    vpsllq zmm24,zmm24,0x4
    vpaddq zmm4,zmm4,zmm1
    vpmadd52huq zmm23,zmm3,zmm13
    vpsraq zmm1,zmm4,0x30
    vpaddq zmm22,zmm22,zmm24
    vpmadd52luq zmm8,zmm3,zmm12
    vpsllq zmm27,zmm27,0x4
    vpaddq zmm22,zmm22,zmm1
    vpmadd52huq zmm21,zmm3,zmm12
    vpsraq zmm1,zmm22,0x30
    vpaddq zmm20,zmm20,zmm27
    vpmadd52huq zmm0,zmm3,zmm11
    vpmadd52luq zmm31,zmm3,zmm11
    vpsllq zmm25,zmm25,0x4
    vpaddq zmm20,zmm20,zmm1
    # Phase C: terminal carry chain and radix-48 normalization.
    # Carries propagate through the remaining limbs, then each limb is masked
    # back to 48 bits to produce LFp output (lazy, not fully reduced mod p).
    vpsraq zmm1,zmm20,0x30
    vpaddq zmm6,zmm6,zmm25
    vpsllq zmm23,zmm23,0x4
    vpaddq zmm6,zmm6,zmm1
    vpsraq zmm1,zmm6,0x30
    vpaddq zmm8,zmm8,zmm23
    vpsllq zmm21,zmm21,0x4
    vpaddq zmm8,zmm8,zmm1
    vpsraq zmm1,zmm8,0x30
    vpaddq zmm3,zmm31,zmm21
    vpandd zmm2,zmm2,zmm10
    vpandd zmm4,zmm4,zmm10
    vpandd zmm22,zmm22,zmm10
    vpandd zmm20,zmm20,zmm10
    vpsllq zmm0,zmm0,0x4
    vpaddq zmm3,zmm3,zmm1
    # Store r[0..3].
    vmovdqa32 ZMMWORD PTR [rdi],zmm2
    vmovdqa32 ZMMWORD PTR [rdi+0x40],zmm4
    vmovdqa32 ZMMWORD PTR [rdi+0x80],zmm22
    vmovdqa32 ZMMWORD PTR [rdi+0xc0],zmm20
    vpaddq zmm0,zmm0,ZMMWORD PTR [rsp+0x88]
    vpsraq zmm1,zmm3,0x30
    vpaddq zmm0,zmm0,zmm1
    vpandd zmm6,zmm6,zmm10
    vpandd zmm8,zmm8,zmm10
    vpandd zmm3,zmm3,zmm10
    vpandd zmm0,zmm0,zmm10
    # Store r[4..7].
    vmovdqa32 ZMMWORD PTR [rdi+0x100],zmm6
    vmovdqa32 ZMMWORD PTR [rdi+0x140],zmm8
    vmovdqa32 ZMMWORD PTR [rdi+0x180],zmm3
    vmovdqa32 ZMMWORD PTR [rdi+0x1c0],zmm0
    # AVX/SSE domain cleanup.
    vzeroupper
    leave
    ret
.size blst_fp8x1_mul_redc_core_k1_asm, .-blst_fp8x1_mul_redc_core_k1_asm

// Compatibility alias used by fused baseline benchmark candidate wiring.
.globl blst_fp8x1_mul_redc_core_k1_candidate_asm
.set blst_fp8x1_mul_redc_core_k1_candidate_asm, blst_fp8x1_mul_redc_core_k1_asm

.globl blst_fp8x1_mul_redc_core_k1_asm_span
.type blst_fp8x1_mul_redc_core_k1_asm_span, @function
blst_fp8x1_mul_redc_core_k1_asm_span:
    # SysV x86_64:
    #   rdi = LFpPack8x1* out
    #   rsi = const FpPack8x1* a
    #   rdx = const FpPack8x1* b
    #   rcx = std::size_t n_packs
    test   rcx, rcx
    je     .Lspan_done

    push   rbx
    mov    r8, rdi
    mov    r9, rsi
    mov    r10, rdx
    mov    rbx, rcx

.Lspan_loop:
    mov    rdi, r8
    mov    rsi, r9
    mov    rdx, r10
    call   blst_fp8x1_mul_redc_core_k1_asm
    add    r8, 0x200
    add    r9, 0x200
    add    r10, 0x200
    dec    rbx
    jne    .Lspan_loop

    pop    rbx
.Lspan_done:
    ret
.size blst_fp8x1_mul_redc_core_k1_asm_span, .-blst_fp8x1_mul_redc_core_k1_asm_span

.section .rodata
.align 8
montw_q:
    # MONT_W_R48 = -p^{-1} mod 2^48.
    .quad 0x0000FFFCFFFCFFFD

.align 64
bmask_vec:
    # BMASK replicated across all 8 lanes.
    .quad 0x0000FFFFFFFFFFFF
    .quad 0x0000FFFFFFFFFFFF
    .quad 0x0000FFFFFFFFFFFF
    .quad 0x0000FFFFFFFFFFFF
    .quad 0x0000FFFFFFFFFFFF
    .quad 0x0000FFFFFFFFFFFF
    .quad 0x0000FFFFFFFFFFFF
    .quad 0x0000FFFFFFFFFFFF

.text
.att_syntax prefix
