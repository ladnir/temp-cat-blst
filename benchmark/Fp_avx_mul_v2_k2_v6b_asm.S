.text
.intel_syntax noprefix

.globl blst_fp8x1_mul_v2_k2_v6b_asm
.type blst_fp8x1_mul_v2_k2_v6b_asm, @function
.p2align 4,0x90
blst_fp8x1_mul_v2_k2_v6b_asm:
    test   rcx,rcx
    je .Lret
    xor    eax,eax
    vpxor  xmm0,xmm0,xmm0
    nop
.Lloop:
    vmovdqa64 zmm18,ZMMWORD PTR [rdx]
    vmovdqa64 zmm19,ZMMWORD PTR [rsi]
    vmovdqa64 zmm1,zmm0
    vpmadd52huq zmm1,zmm19,zmm18
    vmovdqa64 zmm15,ZMMWORD PTR [rdx+0x40]
    vmovdqa64 zmm16,ZMMWORD PTR [rsi+0x40]
    vmovdqa64 zmm2,zmm0
    vpmadd52luq zmm2,zmm19,zmm18
    vpsllq zmm1,zmm1,0x4
    vpmadd52luq zmm1,zmm19,zmm15
    vmovdqa64 zmm11,ZMMWORD PTR [rdx+0xc0]
    vmovdqa64 zmm13,ZMMWORD PTR [rdx+0x80]
    vmovdqa64 ZMMWORD PTR [rdi],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52luq zmm1,zmm16,zmm18
    vpmadd52huq zmm2,zmm19,zmm11
    vmovdqa64 zmm14,ZMMWORD PTR [rsi+0x80]
    vmovdqa64 zmm12,ZMMWORD PTR [rsi+0xc0]
    vmovdqa64 zmm9,ZMMWORD PTR [rdx+0x100]
    vmovdqa64 ZMMWORD PTR [rdi+0x40],zmm1
    vmovdqa64 zmm1,zmm0
    vpmadd52huq zmm1,zmm19,zmm13
    vpmadd52huq zmm2,zmm16,zmm13
    vmovdqa64 zmm10,ZMMWORD PTR [rsi+0x100]
    vmovdqa64 zmm5,ZMMWORD PTR [rdx+0x180]
    vmovdqa64 zmm7,ZMMWORD PTR [rdx+0x140]
    vpmadd52huq zmm1,zmm16,zmm15
    vpmadd52huq zmm2,zmm14,zmm15
    vmovdqa64 zmm8,ZMMWORD PTR [rsi+0x140]
    vmovdqa64 zmm6,ZMMWORD PTR [rsi+0x180]
    vmovdqa64 zmm4,ZMMWORD PTR [rdx+0x1c0]
    vpmadd52huq zmm1,zmm14,zmm18
    vpmadd52huq zmm2,zmm12,zmm18
    vmovdqa64 zmm17,zmm0
    vpmadd52huq zmm17,zmm19,zmm15
    vmovdqa64 zmm3,ZMMWORD PTR [rsi+0x1c0]
    vpsllq zmm1,zmm1,0x4
    vpmadd52luq zmm1,zmm19,zmm11
    vpsllq zmm2,zmm2,0x4
    vpmadd52luq zmm2,zmm19,zmm9
    vpmadd52huq zmm17,zmm16,zmm18
    inc    rax
    vpmadd52luq zmm1,zmm16,zmm13
    add    rdi,0x800
    vpmadd52luq zmm2,zmm16,zmm11
    vpsllq zmm17,zmm17,0x4
    vpmadd52luq zmm17,zmm19,zmm13
    vpmadd52luq zmm1,zmm14,zmm15
    add    rsi,0x400
    vpmadd52luq zmm2,zmm14,zmm13
    add    rdx,0x400
    vpmadd52luq zmm17,zmm16,zmm15
    vpmadd52luq zmm1,zmm12,zmm18
    vpmadd52luq zmm2,zmm12,zmm15
    vpmadd52luq zmm17,zmm14,zmm18
    vmovdqa64 ZMMWORD PTR [rdi-0x740],zmm1
    vmovdqa64 zmm1,zmm0
    vpmadd52huq zmm1,zmm19,zmm9
    vpmadd52luq zmm2,zmm10,zmm18
    vmovdqa64 ZMMWORD PTR [rdi-0x780],zmm17
    vmovdqa64 zmm17,zmm0
    vpmadd52huq zmm17,zmm19,zmm7
    vpmadd52huq zmm1,zmm16,zmm11
    vmovdqa64 ZMMWORD PTR [rdi-0x700],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52huq zmm2,zmm19,zmm5
    vpmadd52huq zmm17,zmm16,zmm9
    vpmadd52huq zmm1,zmm14,zmm13
    vpmadd52huq zmm2,zmm16,zmm7
    vpmadd52huq zmm17,zmm14,zmm11
    vpmadd52huq zmm1,zmm12,zmm15
    vpmadd52huq zmm2,zmm14,zmm9
    vpmadd52huq zmm17,zmm12,zmm13
    vpmadd52huq zmm1,zmm10,zmm18
    vpmadd52huq zmm2,zmm12,zmm11
    vpmadd52huq zmm17,zmm10,zmm15
    vpsllq zmm1,zmm1,0x4
    vpmadd52luq zmm1,zmm19,zmm7
    vpmadd52huq zmm2,zmm10,zmm13
    vpmadd52huq zmm17,zmm8,zmm18
    vpmadd52luq zmm1,zmm16,zmm9
    vpmadd52huq zmm2,zmm8,zmm15
    vpsllq zmm17,zmm17,0x4
    vpmadd52luq zmm17,zmm19,zmm5
    vpmadd52luq zmm1,zmm14,zmm11
    vpmadd52huq zmm2,zmm6,zmm18
    vpmadd52luq zmm17,zmm16,zmm7
    vpmadd52luq zmm1,zmm12,zmm13
    vpsllq zmm2,zmm2,0x4
    vpmadd52luq zmm2,zmm19,zmm4
    vpmadd52luq zmm17,zmm14,zmm9
    vpmadd52luq zmm1,zmm10,zmm15
    vpmadd52luq zmm2,zmm16,zmm5
    vpmadd52luq zmm17,zmm12,zmm11
    vpmadd52luq zmm1,zmm8,zmm18
    vpmadd52luq zmm2,zmm14,zmm7
    vpmadd52luq zmm17,zmm10,zmm13
    vmovdqa64 ZMMWORD PTR [rdi-0x6c0],zmm1
    vmovdqa64 zmm1,zmm0
    vpmadd52huq zmm1,zmm19,zmm4
    vpmadd52luq zmm2,zmm12,zmm9
    vpmadd52luq zmm17,zmm8,zmm15
    vpmadd52huq zmm1,zmm16,zmm5
    vpmadd52luq zmm2,zmm10,zmm11
    vpmadd52luq zmm17,zmm6,zmm18
    vpmadd52huq zmm1,zmm14,zmm7
    vpmadd52luq zmm2,zmm8,zmm13
    vmovdqa64 ZMMWORD PTR [rdi-0x680],zmm17
    vmovdqa64 zmm17,zmm0
    vpmadd52huq zmm1,zmm12,zmm9
    vpmadd52luq zmm2,zmm6,zmm15
    vpmadd52huq zmm1,zmm10,zmm11
    vpmadd52luq zmm2,zmm3,zmm18
    vpmadd52huq zmm1,zmm8,zmm13
    vmovdqa64 ZMMWORD PTR [rdi-0x640],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52huq zmm2,zmm16,zmm4
    vpmadd52huq zmm1,zmm6,zmm15
    vpmadd52huq zmm2,zmm14,zmm5
    vpmadd52huq zmm1,zmm3,zmm18
    vpmadd52huq zmm2,zmm12,zmm7
    vpsllq zmm1,zmm1,0x4
    vpmadd52luq zmm1,zmm16,zmm4
    vpmadd52huq zmm2,zmm10,zmm9
    vpmadd52luq zmm1,zmm14,zmm5
    vpmadd52huq zmm2,zmm8,zmm11
    vpmadd52luq zmm1,zmm12,zmm7
    vpmadd52huq zmm2,zmm6,zmm13
    vpmadd52luq zmm1,zmm10,zmm9
    vpmadd52huq zmm2,zmm3,zmm15
    vpmadd52luq zmm1,zmm8,zmm11
    vpsllq zmm2,zmm2,0x4
    vpmadd52luq zmm2,zmm14,zmm4
    vpmadd52luq zmm1,zmm6,zmm13
    vpmadd52luq zmm2,zmm12,zmm5
    vpmadd52luq zmm1,zmm3,zmm15
    vmovdqa64 zmm15,zmm0
    vpmadd52huq zmm15,zmm14,zmm4
    vpmadd52luq zmm2,zmm10,zmm7
    vmovdqa64 ZMMWORD PTR [rdi-0x600],zmm1
    vpmadd52huq zmm15,zmm12,zmm5
    vpmadd52luq zmm2,zmm8,zmm9
    vmovdqa64 zmm14,zmm15
    vpmadd52huq zmm14,zmm10,zmm7
    vpmadd52luq zmm2,zmm6,zmm11
    vpmadd52huq zmm14,zmm8,zmm9
    vpmadd52luq zmm2,zmm3,zmm13
    vpmadd52huq zmm14,zmm6,zmm11
    vmovdqa64 ZMMWORD PTR [rdi-0x5c0],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52huq zmm2,zmm12,zmm4
    vpmadd52huq zmm14,zmm3,zmm13
    vpmadd52huq zmm2,zmm10,zmm5
    vpsllq zmm1,zmm14,0x4
    vpmadd52luq zmm1,zmm12,zmm4
    vpmadd52huq zmm2,zmm8,zmm7
    vpmadd52luq zmm1,zmm10,zmm5
    vpmadd52huq zmm2,zmm6,zmm9
    vpmadd52luq zmm1,zmm8,zmm7
    vpmadd52huq zmm2,zmm3,zmm11
    vpmadd52luq zmm1,zmm6,zmm9
    vpsllq zmm2,zmm2,0x4
    vpmadd52luq zmm2,zmm10,zmm4
    vpmadd52luq zmm1,zmm3,zmm11
    vmovdqa64 zmm11,zmm0
    vpmadd52huq zmm11,zmm10,zmm4
    vpmadd52luq zmm2,zmm8,zmm5
    vmovdqa64 ZMMWORD PTR [rdi-0x580],zmm1
    vpmadd52huq zmm11,zmm8,zmm5
    vpmadd52luq zmm2,zmm6,zmm7
    vmovdqa64 zmm10,zmm11
    vpmadd52huq zmm10,zmm6,zmm7
    vpmadd52luq zmm2,zmm3,zmm9
    vpmadd52huq zmm10,zmm3,zmm9
    vmovdqa64 ZMMWORD PTR [rdi-0x540],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52huq zmm2,zmm8,zmm4
    vpsllq zmm1,zmm10,0x4
    vpmadd52luq zmm1,zmm8,zmm4
    vpmadd52huq zmm2,zmm6,zmm5
    vpmadd52luq zmm1,zmm6,zmm5
    vpmadd52huq zmm2,zmm3,zmm7
    vpmadd52luq zmm1,zmm3,zmm7
    vmovdqa64 zmm7,zmm0
    vpmadd52huq zmm7,zmm6,zmm4
    vpsllq zmm2,zmm2,0x4
    vpmadd52luq zmm2,zmm6,zmm4
    vmovdqa64 ZMMWORD PTR [rdi-0x500],zmm1
    vmovdqa64 zmm6,zmm7
    vpmadd52huq zmm6,zmm3,zmm5
    vpmadd52luq zmm2,zmm3,zmm5
    vpsllq zmm1,zmm6,0x4
    vmovdqa64 ZMMWORD PTR [rdi-0x4c0],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52luq zmm1,zmm3,zmm4
    vpmadd52huq zmm2,zmm3,zmm4
    vmovdqa64 ZMMWORD PTR [rdi-0x480],zmm1
    vpsllq zmm1,zmm2,0x4
    vmovdqa64 ZMMWORD PTR [rdi-0x440],zmm1
    vmovdqa64 zmm19,ZMMWORD PTR [rsi-0x200]
    vmovdqa64 zmm18,ZMMWORD PTR [rdx-0x200]
    vmovdqa64 zmm1,zmm0
    vpmadd52huq zmm1,zmm19,zmm18
    vmovdqa64 zmm15,ZMMWORD PTR [rdx-0x1c0]
    vmovdqa64 zmm16,ZMMWORD PTR [rsi-0x1c0]
    vmovdqa64 zmm2,zmm0
    vpmadd52luq zmm2,zmm19,zmm18
    vpsllq zmm1,zmm1,0x4
    vpmadd52luq zmm1,zmm19,zmm15
    vmovdqa64 zmm11,ZMMWORD PTR [rdx-0x140]
    vmovdqa64 zmm13,ZMMWORD PTR [rdx-0x180]
    vmovdqa64 ZMMWORD PTR [rdi-0x400],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52luq zmm1,zmm16,zmm18
    vpmadd52huq zmm2,zmm19,zmm11
    vmovdqa64 zmm14,ZMMWORD PTR [rsi-0x180]
    vmovdqa64 zmm12,ZMMWORD PTR [rsi-0x140]
    vmovdqa64 zmm9,ZMMWORD PTR [rdx-0x100]
    vmovdqa64 ZMMWORD PTR [rdi-0x3c0],zmm1
    vmovdqa64 zmm1,zmm0
    vpmadd52huq zmm1,zmm19,zmm13
    vpmadd52huq zmm2,zmm16,zmm13
    vmovdqa64 zmm10,ZMMWORD PTR [rsi-0x100]
    vmovdqa64 zmm5,ZMMWORD PTR [rdx-0x80]
    vmovdqa64 zmm7,ZMMWORD PTR [rdx-0xc0]
    vpmadd52huq zmm1,zmm16,zmm15
    vpmadd52huq zmm2,zmm14,zmm15
    vmovdqa64 zmm8,ZMMWORD PTR [rsi-0xc0]
    vmovdqa64 zmm6,ZMMWORD PTR [rsi-0x80]
    vmovdqa64 zmm4,ZMMWORD PTR [rdx-0x40]
    vpmadd52huq zmm1,zmm14,zmm18
    vpmadd52huq zmm2,zmm12,zmm18
    vpmadd52huq zmm17,zmm19,zmm15
    vmovdqa64 zmm3,ZMMWORD PTR [rsi-0x40]
    vpsllq zmm1,zmm1,0x4
    vpmadd52luq zmm1,zmm19,zmm11
    vpsllq zmm2,zmm2,0x4
    vpmadd52luq zmm2,zmm19,zmm9
    vpmadd52huq zmm17,zmm16,zmm18
    vpmadd52luq zmm1,zmm16,zmm13
    vpmadd52luq zmm2,zmm16,zmm11
    vpsllq zmm17,zmm17,0x4
    vpmadd52luq zmm17,zmm19,zmm13
    vpmadd52luq zmm1,zmm14,zmm15
    vpmadd52luq zmm2,zmm14,zmm13
    vpmadd52luq zmm17,zmm16,zmm15
    vpmadd52luq zmm1,zmm12,zmm18
    vpmadd52luq zmm2,zmm12,zmm15
    vpmadd52luq zmm17,zmm14,zmm18
    vmovdqa64 ZMMWORD PTR [rdi-0x340],zmm1
    vmovdqa64 zmm1,zmm0
    vpmadd52huq zmm1,zmm19,zmm9
    vpmadd52luq zmm2,zmm10,zmm18
    vmovdqa64 ZMMWORD PTR [rdi-0x380],zmm17
    vmovdqa64 zmm17,zmm0
    vpmadd52huq zmm17,zmm19,zmm7
    vpmadd52huq zmm1,zmm16,zmm11
    vmovdqa64 ZMMWORD PTR [rdi-0x300],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52huq zmm2,zmm19,zmm5
    vpmadd52huq zmm17,zmm16,zmm9
    vpmadd52huq zmm1,zmm14,zmm13
    vpmadd52huq zmm2,zmm16,zmm7
    vpmadd52huq zmm17,zmm14,zmm11
    vpmadd52huq zmm1,zmm12,zmm15
    vpmadd52huq zmm2,zmm14,zmm9
    vpmadd52huq zmm17,zmm12,zmm13
    vpmadd52huq zmm1,zmm10,zmm18
    vpmadd52huq zmm2,zmm12,zmm11
    vpmadd52huq zmm17,zmm10,zmm15
    vpsllq zmm1,zmm1,0x4
    vpmadd52luq zmm1,zmm19,zmm7
    vpmadd52huq zmm2,zmm10,zmm13
    vpmadd52huq zmm17,zmm8,zmm18
    vpmadd52luq zmm1,zmm16,zmm9
    vpmadd52huq zmm2,zmm8,zmm15
    vpsllq zmm17,zmm17,0x4
    vpmadd52luq zmm17,zmm19,zmm5
    vpmadd52luq zmm1,zmm14,zmm11
    vpmadd52huq zmm2,zmm6,zmm18
    vpmadd52luq zmm17,zmm16,zmm7
    vpmadd52luq zmm1,zmm12,zmm13
    vpsllq zmm2,zmm2,0x4
    vpmadd52luq zmm2,zmm19,zmm4
    vpmadd52luq zmm17,zmm14,zmm9
    vpmadd52luq zmm1,zmm10,zmm15
    vpmadd52luq zmm2,zmm16,zmm5
    vpmadd52luq zmm17,zmm12,zmm11
    vpmadd52luq zmm1,zmm8,zmm18
    vpmadd52luq zmm2,zmm14,zmm7
    vpmadd52luq zmm17,zmm10,zmm13
    vmovdqa64 ZMMWORD PTR [rdi-0x2c0],zmm1
    vmovdqa64 zmm1,zmm0
    vpmadd52huq zmm1,zmm19,zmm4
    vpmadd52luq zmm2,zmm12,zmm9
    vpmadd52luq zmm17,zmm8,zmm15
    vpmadd52huq zmm1,zmm16,zmm5
    vpmadd52luq zmm2,zmm10,zmm11
    vpmadd52luq zmm17,zmm6,zmm18
    vpmadd52huq zmm1,zmm14,zmm7
    vpmadd52luq zmm2,zmm8,zmm13
    vmovdqa64 ZMMWORD PTR [rdi-0x280],zmm17
    vpmadd52huq zmm1,zmm12,zmm9
    vpmadd52luq zmm2,zmm6,zmm15
    vpmadd52huq zmm1,zmm10,zmm11
    vpmadd52luq zmm2,zmm3,zmm18
    vpmadd52huq zmm1,zmm8,zmm13
    vmovdqa64 ZMMWORD PTR [rdi-0x240],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52huq zmm2,zmm16,zmm4
    vpmadd52huq zmm1,zmm6,zmm15
    vpmadd52huq zmm2,zmm14,zmm5
    vpmadd52huq zmm1,zmm3,zmm18
    vpmadd52huq zmm2,zmm12,zmm7
    vpsllq zmm1,zmm1,0x4
    vpmadd52luq zmm1,zmm16,zmm4
    vpmadd52huq zmm2,zmm10,zmm9
    vpmadd52luq zmm1,zmm14,zmm5
    vpmadd52huq zmm2,zmm8,zmm11
    vpmadd52luq zmm1,zmm12,zmm7
    vpmadd52huq zmm2,zmm6,zmm13
    vpmadd52luq zmm1,zmm10,zmm9
    vpmadd52huq zmm2,zmm3,zmm15
    vpmadd52luq zmm1,zmm8,zmm11
    vpsllq zmm2,zmm2,0x4
    vpmadd52luq zmm2,zmm14,zmm4
    vpmadd52luq zmm1,zmm6,zmm13
    vpmadd52luq zmm2,zmm12,zmm5
    vpmadd52luq zmm1,zmm3,zmm15
    vmovdqa64 zmm15,zmm0
    vpmadd52huq zmm15,zmm14,zmm4
    vpmadd52luq zmm2,zmm10,zmm7
    vmovdqa64 ZMMWORD PTR [rdi-0x200],zmm1
    vpmadd52huq zmm15,zmm12,zmm5
    vpmadd52luq zmm2,zmm8,zmm9
    vmovdqa64 zmm14,zmm15
    vpmadd52huq zmm14,zmm10,zmm7
    vpmadd52luq zmm2,zmm6,zmm11
    vpmadd52huq zmm14,zmm8,zmm9
    vpmadd52luq zmm2,zmm3,zmm13
    vpmadd52huq zmm14,zmm6,zmm11
    vmovdqa64 ZMMWORD PTR [rdi-0x1c0],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52huq zmm2,zmm12,zmm4
    vpmadd52huq zmm14,zmm3,zmm13
    vpmadd52huq zmm2,zmm10,zmm5
    vpsllq zmm1,zmm14,0x4
    vpmadd52luq zmm1,zmm12,zmm4
    vpmadd52huq zmm2,zmm8,zmm7
    vpmadd52luq zmm1,zmm10,zmm5
    vpmadd52huq zmm2,zmm6,zmm9
    vpmadd52luq zmm1,zmm8,zmm7
    vpmadd52huq zmm2,zmm3,zmm11
    vpmadd52luq zmm1,zmm6,zmm9
    vpsllq zmm2,zmm2,0x4
    vpmadd52luq zmm2,zmm10,zmm4
    vpmadd52luq zmm1,zmm3,zmm11
    vmovdqa64 zmm11,zmm0
    vpmadd52huq zmm11,zmm10,zmm4
    vpmadd52luq zmm2,zmm8,zmm5
    vmovdqa64 ZMMWORD PTR [rdi-0x180],zmm1
    vmovdqa64 zmm10,zmm11
    vpmadd52huq zmm10,zmm8,zmm5
    vpmadd52luq zmm2,zmm6,zmm7
    vpmadd52huq zmm10,zmm6,zmm7
    vpmadd52luq zmm2,zmm3,zmm9
    vpmadd52huq zmm10,zmm3,zmm9
    vmovdqa64 ZMMWORD PTR [rdi-0x140],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52huq zmm2,zmm8,zmm4
    vpsllq zmm1,zmm10,0x4
    vpmadd52luq zmm1,zmm8,zmm4
    vpmadd52huq zmm2,zmm6,zmm5
    vpmadd52luq zmm1,zmm6,zmm5
    vpmadd52huq zmm2,zmm3,zmm7
    vpmadd52luq zmm1,zmm3,zmm7
    vmovdqa64 zmm7,zmm0
    vpmadd52huq zmm7,zmm6,zmm4
    vpsllq zmm2,zmm2,0x4
    vpmadd52luq zmm2,zmm6,zmm4
    vmovdqa64 ZMMWORD PTR [rdi-0x100],zmm1
    vmovdqa64 zmm6,zmm7
    vpmadd52huq zmm6,zmm3,zmm5
    vpmadd52luq zmm2,zmm3,zmm5
    vpsllq zmm1,zmm6,0x4
    vmovdqa64 ZMMWORD PTR [rdi-0xc0],zmm2
    vmovdqa64 zmm2,zmm0
    vpmadd52luq zmm1,zmm3,zmm4
    vpmadd52huq zmm2,zmm3,zmm4
    vmovdqa64 ZMMWORD PTR [rdi-0x80],zmm1
    vpsllq zmm1,zmm2,0x4
    vmovdqa64 ZMMWORD PTR [rdi-0x40],zmm1
    cmp    rcx,rax
    jne .Lloop
    vzeroupper
.Lret:
    ret
.size blst_fp8x1_mul_v2_k2_v6b_asm, .-blst_fp8x1_mul_v2_k2_v6b_asm
