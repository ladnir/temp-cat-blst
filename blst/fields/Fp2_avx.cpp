#include "blst/fields/Fp2_avx.h"
#include "blst/backend/avx512/avx_layout.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{

    void add_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a, const fp2_8x1x1w& b)
    {
        FpPack8x1::add(r[0], a[0], b[0]);
        FpPack8x1::add(r[1], a[1], b[1]);
    }

    void add_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b)
    {
        FpPack8x1::add(r, a, b);
    }

    void sub_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b)
    {
        FpPack8x1::sub(r, a, b);
    }

    // a = < D1 | D0 | C1 | C0 | B1 | B0 | A1 | A0 >  
    // b = < H1 | H0 | G1 | G0 | F1 | F0 | E1 | E0 >
    // r = < D1+H1 | D0+H0 | C1+G1 | C0+G0 | B1-F1 | B0-F0 | A1-E1 | A0-E0 >
    void aass_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b)
    {
        const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
        const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
        const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
        const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
        const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
        const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
        const cc512i bmask = VSET1(BMASK);
        cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask;
        cc512i t0, t1, t2, t3, t4, t5, t6, t7;

        // r = D1+H1 | D0+H0 | C1+G1 | C0+G0 | B1 | B0 | A1 | A0
        r0 = VMADD(a0, 0xF0, a0, b0); r1 = VMADD(a1, 0xF0, a1, b1);
        r2 = VMADD(a2, 0xF0, a2, b2); r3 = VMADD(a3, 0xF0, a3, b3);
        r4 = VMADD(a4, 0xF0, a4, b4); r5 = VMADD(a5, 0xF0, a5, b5);
        r6 = VMADD(a6, 0xF0, a6, b6); r7 = VMADD(a7, 0xF0, a7, b7);

        // t = p | p | p | p | F1 | F0 | E1 | E0
        t0 = VMBLEND(0xF0, b0, p0); t1 = VMBLEND(0xF0, b1, p1);
        t2 = VMBLEND(0xF0, b2, p2); t3 = VMBLEND(0xF0, b3, p3);
        t4 = VMBLEND(0xF0, b4, p4); t5 = VMBLEND(0xF0, b5, p5);
        t6 = VMBLEND(0xF0, b6, p6); t7 = VMBLEND(0xF0, b7, p7);

        // r = D1+H1-p | D0+H0-p | C1+G1-p | C0+G0-p | B1-F1 | B0-F0 | A1-E1 | A0-E0
        r0 = VSUB(r0, t0); r1 = VSUB(r1, t1); r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
        r4 = VSUB(r4, t4); r5 = VSUB(r5, t5); r6 = VSUB(r6, t6); r7 = VSUB(r7, t7);

        // get sign mask 
        t0 = VADD(r1, VSRA(r0, BRADIX));
        t0 = VADD(r2, VSRA(t0, BRADIX));
        t0 = VADD(r3, VSRA(t0, BRADIX));
        t0 = VADD(r4, VSRA(t0, BRADIX));
        t0 = VADD(r5, VSRA(t0, BRADIX));
        t0 = VADD(r6, VSRA(t0, BRADIX));
        t0 = VADD(r7, VSRA(t0, BRADIX));

        // if r is non-negative, smask = all-0 
        // if r is     negative, smask = all-1
        smask = VSRA(t0, 63);
        // r = r + (p & smask)
        r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask));
        r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask));
        r4 = VADD(r4, VAND(p4, smask)); r5 = VADD(r5, VAND(p5, smask));
        r6 = VADD(r6, VAND(p6, smask)); r7 = VADD(r7, VAND(p7, smask));

        // carry propagation 
        r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
        r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
        r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
        r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
        r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
        r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
        r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
        r7 = VAND(r7, bmask);

        r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
        r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
    }

    // a = < D1 | D0 | C1 | C0 | B1 | B0 | A1 | A0 >  
    // b = < H1 | H0 | G1 | G0 | F1 | F0 | E1 | E0 >
    // r = < D1+H1 | D0+H0 | C1-G1 | C0-G0 | B1-F1 | B0-F0 | A1+E1 | A0+E0 >
    void assa_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b)
    {
        const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
        const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
        const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
        const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
        const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
        const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
        const cc512i bmask = VSET1(BMASK);
        cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask;
        cc512i t0, t1, t2, t3, t4, t5, t6, t7;

        // r = D1+H1 | D0+H0 | C1 | C0 | B1 | B0 | A1+E1 | A0+E0
        r0 = VMADD(a0, 0xC3, a0, b0); r1 = VMADD(a1, 0xC3, a1, b1);
        r2 = VMADD(a2, 0xC3, a2, b2); r3 = VMADD(a3, 0xC3, a3, b3);
        r4 = VMADD(a4, 0xC3, a4, b4); r5 = VMADD(a5, 0xC3, a5, b5);
        r6 = VMADD(a6, 0xC3, a6, b6); r7 = VMADD(a7, 0xC3, a7, b7);

        // t = p | p | G1 | G0 | F1 | F0 | p | p
        t0 = VMBLEND(0xC3, b0, p0); t1 = VMBLEND(0xC3, b1, p1);
        t2 = VMBLEND(0xC3, b2, p2); t3 = VMBLEND(0xC3, b3, p3);
        t4 = VMBLEND(0xC3, b4, p4); t5 = VMBLEND(0xC3, b5, p5);
        t6 = VMBLEND(0xC3, b6, p6); t7 = VMBLEND(0xC3, b7, p7);

        // r = D1+H1-p | D0+H0-p | C1-G1 | C0-G0 | B1-F1 | B0-F0 | A1+E1-p | A0+E0-p
        r0 = VSUB(r0, t0); r1 = VSUB(r1, t1); r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
        r4 = VSUB(r4, t4); r5 = VSUB(r5, t5); r6 = VSUB(r6, t6); r7 = VSUB(r7, t7);

        // get sign mask 
        t0 = VADD(r1, VSRA(r0, BRADIX));
        t0 = VADD(r2, VSRA(t0, BRADIX));
        t0 = VADD(r3, VSRA(t0, BRADIX));
        t0 = VADD(r4, VSRA(t0, BRADIX));
        t0 = VADD(r5, VSRA(t0, BRADIX));
        t0 = VADD(r6, VSRA(t0, BRADIX));
        t0 = VADD(r7, VSRA(t0, BRADIX));

        // if r is non-negative, smask = all-0 
        // if r is     negative, smask = all-1
        smask = VSRA(t0, 63);
        // r = r + (p & smask)
        r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask));
        r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask));
        r4 = VADD(r4, VAND(p4, smask)); r5 = VADD(r5, VAND(p5, smask));
        r6 = VADD(r6, VAND(p6, smask)); r7 = VADD(r7, VAND(p7, smask));

        // carry propagation 
        r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
        r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
        r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
        r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
        r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
        r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
        r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
        r7 = VAND(r7, bmask);

        r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
        r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
    }

    // a = < D1 | D0 | C1 | C0 | B1 | B0 | A1 | A0 >  
    // b = < H1 | H0 | G1 | G0 | F1 | F0 | E1 | E0 >
    // r = < D1-H1 | D0-H0 | C1+G1 | C0+G0 | B1-F1 | B0-F0 | A1-E1 | A0-E0 >
    void sass_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b)
    {
        const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
        const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
        const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
        const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
        const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
        const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
        const cc512i bmask = VSET1(BMASK);
        cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask;
        cc512i t0, t1, t2, t3, t4, t5, t6, t7;

        // r = D1 | D0 | C1+G1 | C0+G0 | B1 | B0 | A1 | A0
        r0 = VMADD(a0, 0x30, a0, b0); r1 = VMADD(a1, 0x30, a1, b1);
        r2 = VMADD(a2, 0x30, a2, b2); r3 = VMADD(a3, 0x30, a3, b3);
        r4 = VMADD(a4, 0x30, a4, b4); r5 = VMADD(a5, 0x30, a5, b5);
        r6 = VMADD(a6, 0x30, a6, b6); r7 = VMADD(a7, 0x30, a7, b7);

        // t = H1 | H0 | p | p | F1 | F0 | E1 | E0
        t0 = VMBLEND(0x30, b0, p0); t1 = VMBLEND(0x30, b1, p1);
        t2 = VMBLEND(0x30, b2, p2); t3 = VMBLEND(0x30, b3, p3);
        t4 = VMBLEND(0x30, b4, p4); t5 = VMBLEND(0x30, b5, p5);
        t6 = VMBLEND(0x30, b6, p6); t7 = VMBLEND(0x30, b7, p7);

        // r = D1-H1 | D0-H0 | C1+G1-p | C0+G0-p | B1-F1 | B0-F0 | A1-E1 | A0-E0
        r0 = VSUB(r0, t0); r1 = VSUB(r1, t1); r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
        r4 = VSUB(r4, t4); r5 = VSUB(r5, t5); r6 = VSUB(r6, t6); r7 = VSUB(r7, t7);

        // get sign mask 
        t0 = VADD(r1, VSRA(r0, BRADIX));
        t0 = VADD(r2, VSRA(t0, BRADIX));
        t0 = VADD(r3, VSRA(t0, BRADIX));
        t0 = VADD(r4, VSRA(t0, BRADIX));
        t0 = VADD(r5, VSRA(t0, BRADIX));
        t0 = VADD(r6, VSRA(t0, BRADIX));
        t0 = VADD(r7, VSRA(t0, BRADIX));

        // if r is non-negative, smask = all-0 
        // if r is     negative, smask = all-1
        smask = VSRA(t0, 63);
        // r = r + (p & smask)
        r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask));
        r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask));
        r4 = VADD(r4, VAND(p4, smask)); r5 = VADD(r5, VAND(p5, smask));
        r6 = VADD(r6, VAND(p6, smask)); r7 = VADD(r7, VAND(p7, smask));

        // carry propagation 
        r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
        r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
        r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
        r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
        r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
        r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
        r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
        r7 = VAND(r7, bmask);

        r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
        r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
    }

    // a = < D1 | D0 | C1 | C0 | B1 | B0 | A1 | A0 >  
    // b = < H1 | H0 | G1 | G0 | F1 | F0 | E1 | E0 >
    // r = < D1-H1 | D0-H0 | C1-G1 | C0-G0 | B1+F1 | B0+F0 | A1+E1 | A0+E0 >
    void ssaa_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b)
    {
        const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
        const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
        const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
        const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
        const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
        const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
        const cc512i bmask = VSET1(BMASK);
        cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask;
        cc512i t0, t1, t2, t3, t4, t5, t6, t7;

        // r = D1 | D0 | C1 | C0 | B1+F1 | B0+F0 | A1+E1 | A0+E0
        r0 = VMADD(a0, 0x0F, a0, b0); r1 = VMADD(a1, 0x0F, a1, b1);
        r2 = VMADD(a2, 0x0F, a2, b2); r3 = VMADD(a3, 0x0F, a3, b3);
        r4 = VMADD(a4, 0x0F, a4, b4); r5 = VMADD(a5, 0x0F, a5, b5);
        r6 = VMADD(a6, 0x0F, a6, b6); r7 = VMADD(a7, 0x0F, a7, b7);

        // t = H1 | H0 | G1 | G0 | p | p | p | p
        t0 = VMBLEND(0x0F, b0, p0); t1 = VMBLEND(0x0F, b1, p1);
        t2 = VMBLEND(0x0F, b2, p2); t3 = VMBLEND(0x0F, b3, p3);
        t4 = VMBLEND(0x0F, b4, p4); t5 = VMBLEND(0x0F, b5, p5);
        t6 = VMBLEND(0x0F, b6, p6); t7 = VMBLEND(0x0F, b7, p7);

        // r = D1-H1 | D0-H0 | C1-G1 | C0-G0 | B1+F1-p | B0+F0-p | A1+E1-p | A0+E0-p
        r0 = VSUB(r0, t0); r1 = VSUB(r1, t1); r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
        r4 = VSUB(r4, t4); r5 = VSUB(r5, t5); r6 = VSUB(r6, t6); r7 = VSUB(r7, t7);

        // get sign mask 
        t0 = VADD(r1, VSRA(r0, BRADIX));
        t0 = VADD(r2, VSRA(t0, BRADIX));
        t0 = VADD(r3, VSRA(t0, BRADIX));
        t0 = VADD(r4, VSRA(t0, BRADIX));
        t0 = VADD(r5, VSRA(t0, BRADIX));
        t0 = VADD(r6, VSRA(t0, BRADIX));
        t0 = VADD(r7, VSRA(t0, BRADIX));

        // if r is non-negative, smask = all-0 
        // if r is     negative, smask = all-1
        smask = VSRA(t0, 63);
        // r = r + (p & smask)
        r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask));
        r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask));
        r4 = VADD(r4, VAND(p4, smask)); r5 = VADD(r5, VAND(p5, smask));
        r6 = VADD(r6, VAND(p6, smask)); r7 = VADD(r7, VAND(p7, smask));

        // carry propagation 
        r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
        r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
        r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
        r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
        r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
        r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
        r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
        r7 = VAND(r7, bmask);

        r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
        r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
    }

    // r0 = a0 - a1
    // r1 = a0 + a1
    void mul_by_u_plus_1_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a)
    {
        fp2_4x2x1w t0;

        // a = A1 | A0 at Fp layer 
        shuf_01(t0.data(), a.data());         //     A0 |    A1
        FpPack8x1::asx4(r, a, t0);               //  A0+A1 | A0-A1
    }

    // r0 = a0*b0 - a1*b1
    // r1 = a0*b1 + a1*b0
    // schoolbook
    void mul_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b)
    {
        fp2_4x2x1w t0, t1;
        fp2x2_4x2x1w tt0, tt1;

        // a = A1 | A0 at Fp layer 
        // b = B1 | B0 at Fp layer
        shuf_00(t0.data(), a.data());         //          A0 |          A0 
        FpPack8x1::mul(tt0, t0, b);            //       A0*B1 |       A0*B0
        shuf_11(t0.data(), a.data());         //          A1 |          A1 
        shuf_01(t1.data(), b.data());         //          B0 |          B1 
        FpPack8x1::mul(tt1, t0, t1);           //       A1*B0 |       A1*B1
        UFpPack8x1::asx4(tt0, tt0, tt1);        // A0*B1+A1*B0 | A0*B0-A1*B1
        UFpPack8x1::redc(r, tt0);               // A0*B1+A1*B0 | A0*B0-A1*B1
    }

    // r0 = (a0 + a1)*(a0 - a1)
    // r1 = 2*a0*a1
    void sqr_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a)
    {
        fp2_4x2x1w t0, t1, t2;
        fp2x2_4x2x1w tt0;

        // a = A1 | A0 at Fp layer 
        shuf_00(t0.data(), a.data());         //      A0 |              A0
        shuf_01(t1.data(), a.data());         //      A0 |              A1 
        shuf_z1(t2.data(), a.data());         //       0 |              A1
        FpPack8x1::add(t0, t0, t1);              //    2*A0 |           A0+A1
        FpPack8x1::sub(t2, a, t2);               //      A1 |           A0-A1
        FpPack8x1::mul(tt0, t0, t2);           // 2*A0*A1 | (A0+A1)*(A0-A1)
        UFpPack8x1::redc(r, tt0);               // 2*A0*A1 | (A0+A1)*(A0-A1)
    }

    // r0 = a0*b0 - a1*b1
    // r1 = a0*b1 + a1*b0
    // schoolbook
    void mul_fp2_2x4x1w(fp2_2x4x1w& r, const fp2_2x4x1w& a, const fp2_2x4x1w& b)
    {
        fp2_2x4x1w t0, t1;
        fp2x2_2x4x1w tt0;

        // a = A1 | A0 | ... | ... at Fp layer
        // b = B1 | B0 | ... | ... at Fp layer
        perm_3322(t0.data(), a.data());       //        A1 |        A1 |   A0 |   A0
        perm_2332(t1.data(), b.data());       //        B0 |        B1 |   B1 |   B0
        FpPack8x1::mul(tt0, t0, t1);           //      A1B0 |      A1B1 | A0B1 | A0B0
        UFpPack8x1::redc(t0, tt0);              //      A1B0 |      A1B1 | A0B1 | A0B0
        perm_10zz(t1.data(), t0.data());      //      A0B1 |      A0B0 |    0 |    0 
        FpPack8x1::asx4(r, t1, t0);              // A0B1+A1B0 | A0B0-A1B1 |  ... |  ...
    }

    void add_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b)
    {
        FpPack4x2::add(r, a, b);
    }

    void sub_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b)
    {
        FpPack4x2::sub(r, a, b);
    }

    // a = < B1' | B1 | B0' | B0 | A1' | A1 | A0' | A0 >
    // b = < D1' | D1 | D0' | D0 | C1' | C1 | C0' | C0 >
    // r = < B1'+D1' | B1+D1 | B0'+D0' | B0+D0 | A1'-C1' | A1-C1 | A0'-C0' | A0-C0 >
    void as_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b)
    {
        const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
        const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
        const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
        const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
        const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
        const cc512i bmask = VSET1(BMASK);
        cc512i r0, r1, r2, r3, smask;
        cc512i t0, t1, t2, t3;

        // r =  B1'+D1' | B1+D1 | B0'+D0' | B0+D0 | A1' | A1 | A0' | A0 
        r0 = VMADD(a0, 0xF0, a0, b0); r1 = VMADD(a1, 0xF0, a1, b1);
        r2 = VMADD(a2, 0xF0, a2, b2); r3 = VMADD(a3, 0xF0, a3, b3);

        // t = p' | p | p' | p | C1' | C1 | C0' | C0
        t0 = VMBLEND(0xF0, b0, p0); t1 = VMBLEND(0xF0, b1, p1);
        t2 = VMBLEND(0xF0, b2, p2); t3 = VMBLEND(0xF0, b3, p3);

        // r = B1'+D1'-p' | B1+D1-p | B0'+D0'-p' | B0+D0-p | A1'-C1' | A1-C1 | A0'-C0' | A0-C0
        r0 = VSUB(r0, t0); r1 = VSUB(r1, t1); r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);

        // get sign mask
        t0 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX));
        t0 = VMADD(r2, 0x55, r2, VSRA(t0, BRADIX));
        t0 = VMADD(r3, 0x55, r3, VSRA(t0, BRADIX));
        t0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
        t0 = VMADD(r1, 0xAA, r1, VSRA(t0, BRADIX));
        t0 = VMADD(r2, 0xAA, r2, VSRA(t0, BRADIX));
        t0 = VMADD(r3, 0xAA, r3, VSRA(t0, BRADIX));

        // if r is non-negative, smask = all-0 
        // if r is     negative, smask = all-1
        smask = VSRA(t0, 63);
        smask = VSHUF<0xEE>(smask);

        // r = r + (p & smask)
        r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask));
        r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask));

        // carry propagation 
        r1 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX));
        r2 = VMADD(r2, 0x55, r2, VSRA(r1, BRADIX));
        r3 = VMADD(r3, 0x55, r3, VSRA(r2, BRADIX));
        r0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(r3, BRADIX)));
        r1 = VMADD(r1, 0xAA, r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
        r2 = VMADD(r2, 0xAA, r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
        r3 = VMADD(r3, 0xAA, r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
        r3 = VAND(r3, bmask);

        r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
    }

    // a = < B1' | B1 | B0' | B0 | A1' | A1 | A0' | A0 >
    // b = < D1' | D1 | D0' | D0 | C1' | C1 | C0' | C0 >
    // r = < B1'-D1' | B1-D1 | B0'-D0' | B0-D0 | A1'+C1' | A1+C1 | A0'+C0' | A0+C0 >
    void sa_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b)
    {
        const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
        const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
        const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
        const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
        const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
        const cc512i bmask = VSET1(BMASK);
        cc512i r0, r1, r2, r3, smask;
        cc512i t0, t1, t2, t3;

        // r =  B1' | B1 | B0' | B0 | A1'+C1' | A1+C1 | A0'+C0' | A0+C0 
        r0 = VMADD(a0, 0x0F, a0, b0); r1 = VMADD(a1, 0x0F, a1, b1);
        r2 = VMADD(a2, 0x0F, a2, b2); r3 = VMADD(a3, 0x0F, a3, b3);

        // t = D1' | D1 | D0' | D0 | p' | p | p' | p
        t0 = VMBLEND(0x0F, b0, p0); t1 = VMBLEND(0x0F, b1, p1);
        t2 = VMBLEND(0x0F, b2, p2); t3 = VMBLEND(0x0F, b3, p3);

        // r = B1'-D1' | B1-D1 | B0'-D0' | B0-D0 | A1'+C1'-p' | A1+C1-p' | A0'+C0'-p' | A0+C0-p'
        r0 = VSUB(r0, t0); r1 = VSUB(r1, t1); r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);

        // get sign mask
        t0 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX));
        t0 = VMADD(r2, 0x55, r2, VSRA(t0, BRADIX));
        t0 = VMADD(r3, 0x55, r3, VSRA(t0, BRADIX));
        t0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
        t0 = VMADD(r1, 0xAA, r1, VSRA(t0, BRADIX));
        t0 = VMADD(r2, 0xAA, r2, VSRA(t0, BRADIX));
        t0 = VMADD(r3, 0xAA, r3, VSRA(t0, BRADIX));

        // if r is non-negative, smask = all-0 
        // if r is     negative, smask = all-1
        smask = VSRA(t0, 63);
        smask = VSHUF<0xEE>(smask);

        // r = r + (p & smask)
        r0 = VADD(r0, VAND(p0, smask)); r1 = VADD(r1, VAND(p1, smask));
        r2 = VADD(r2, VAND(p2, smask)); r3 = VADD(r3, VAND(p3, smask));

        // carry propagation 
        r1 = VMADD(r1, 0x55, r1, VSRA(r0, BRADIX));
        r2 = VMADD(r2, 0x55, r2, VSRA(r1, BRADIX));
        r3 = VMADD(r3, 0x55, r3, VSRA(r2, BRADIX));
        r0 = VMADD(r0, 0xAA, r0, VZSHUF<0x4E>(0xCCCC, VSRA(r3, BRADIX)));
        r1 = VMADD(r1, 0xAA, r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
        r2 = VMADD(r2, 0xAA, r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
        r3 = VMADD(r3, 0xAA, r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
        r3 = VAND(r3, bmask);

        r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
    }

    // r0 = a0 - a1
    // r1 = a0 + a1
    void mul_by_u_plus_1_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a)
    {
        fp2_2x2x2w t0;

        // a = A1 | A0 at Fp layer 
        perm_1032_hl(t0.data(), a.data());    //     A0 |    A1
        FpPack4x2::asx2(r, a, t0);               //  A0+A1 | A0-A1
    }

    // r0 = a0*b0 - a1*b1
    // r1 = a0*b1 + a1*b0
    // schoolbook
    void mul_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b)
    {
        fp2_2x2x2w t0, t1;
        fp2x2_2x2x2w tt0, tt1;

        // a = A1 | A0 at Fp layer
        // b = B1 | B0 at Fp layer
        perm_1010_hl(t0.data(), a.data());    //          A0 |          A0
        FpPack4x2::mul(tt0, t0, b);            //       A0*B1 |       A0*B0
        perm_3232_hl(t0.data(), a.data());    //          A1 |          A1
        perm_1032_hl(t1.data(), b.data());    //          B0 |          B1
        FpPack4x2::mul(tt1, t0, t1);           //       A1*B0 |       A1*B1
        UFpPack4x2::asx2(tt0, tt0, tt1);        // A0*B1+A1*B0 | A0*B0-A1*B1
        UFpPack4x2::redc(r, tt0);               // A0*B1+A1*B0 | A0*B0-A1*B1
    }

    // r0 = (a0 + a1)*(a0 - a1)
    // r1 = 2*a0*a1
    void sqr_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a)
    {
        fp2_2x2x2w t0, t1, t2;
        fp2x2_2x2x2w tt0;

        // a = A1 | A0 at Fp layer
        perm_1010_hl(t0.data(), a.data());    //      A0 |              A0
        perm_1032_hl(t1.data(), a.data());    //      A0 |              A1 
        perm_zz32_hl(t2.data(), a.data());    //       0 |              A1
        FpPack4x2::add(t0, t0, t1);              //    2*A0 |           A0+A1
        FpPack4x2::sub(t2, a, t2);               //      A1 |           A0-A1
        FpPack4x2::mul(tt0, t0, t2);           // 2*A0*A1 | (A0+A1)*(A0-A1)
        UFpPack4x2::redc(r, tt0);               // 2*A0*A1 | (A0+A1)*(A0-A1)
    }

    // r0 = a0*b0 - a1*b1
    // r1 = a0*b1 + a1*b0
    // schoolbook
    void mul_fp2_1x4x2w(fp2_1x4x2w& r, const fp2_1x4x2w& a, const fp2_1x4x2w& b)
    {
        fp2_1x4x2w t0, t1;
        fp2x2_1x4x2w tt0;
        const cc512i m0 = VSET(7, 6, 7, 6, 5, 4, 5, 4);
        const cc512i m1 = VSET(5, 4, 7, 6, 7, 6, 5, 4);
        const cc512i m2 = VSET(3, 2, 1, 0, 7, 6, 5, 4);

        // a = A1 | A0 | ... | ... at Fp layer
        // b = B1 | B0 | ... | ... at Fp layer
        perm_var_hl(t0.data(), a.data(), m0); //          A1 |          A1 |    A0 |    A0
        perm_var_hl(t1.data(), b.data(), m1); //          B0 |          B1 |    B1 |    B0
        FpPack4x2::mul(tt0, t0, t1);           //       A1*B0 |       A1*B1 | A0*B1 | A0*B0
        UFpPack4x2::redc(t0, tt0);              //       A1*B0 |       A1*B1 | A0*B1 | A0*B0
        perm_var_hl(t1.data(), t0.data(), m2);//       A0*B1 |       A0*B0 |   ... |   ... 
        FpPack4x2::asx2(r, t1, t0);              // A0*B1+A1*B0 | A0*B0-A1*B1 |   ... |   ...
    }

    // ----------------------------------------------------------------------------
    // Fp2 double-length operations

    void add_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const fp2x2_8x1x1w& b)
    {
        UFpPack8x1::add(r[0], a[0], b[0]);
        UFpPack8x1::add(r[1], a[1], b[1]);
    }

    void sub_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const fp2x2_8x1x1w& b)
    {
        UFpPack8x1::sub(r[0], a[0], b[0]);
        UFpPack8x1::sub(r[1], a[1], b[1]);
    }

    void asx4_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const fp2x2_8x1x1w& b)
    {
        UFpPack8x1::asx4(r[0], a[0], b[0]);
        UFpPack8x1::asx4(r[1], a[1], b[1]);
    }

    // r0 = a0 - a1
    // r1 = a0 + a1
    void mul_by_u_plus_1_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a)
    {
        UFpPack8x1::sub(r[0], a[0], a[1]);      // a0-a1
        UFpPack8x1::add(r[1], a[0], a[1]);      // a0+a1
    }

    // r0 = a0*b0 - a1*b1
    // r1 = (a0 + a1)*(b0 + b1) - a0*b0 - a1*b1 = a0*b1 + a1*b0
    // Karatsuba
    void mul_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2_8x1x1w& a, const fp2_8x1x1w& b)
    {
        FpPack8x1 t0, t1;
        UFpPack8x1 tt0, tt1, tt2;

        FpPack8x1::mul(tt0, a[0], b[0]);       // a0*b0
        FpPack8x1::mul(tt1, a[1], b[1]);       // a1*b1
        FpPack8x1::add(t0, a[0], a[1]);          // a0+a1
        FpPack8x1::add(t1, b[0], b[1]);          // b0+b1
        FpPack8x1::mul(tt2, t0, t1);           // (a0+a1)*(b0+b1)
        UFpPack8x1::sub(r[0], tt0, tt1);        // a0*b0-a1*b1
        UFpPack8x1::sub(tt2, tt2, tt0);         // (a0+a1)*(b0+b1)-a0*b0
        UFpPack8x1::sub(r[1], tt2, tt1);        // (a0+a1)*(b0+b1)-a0*b0-a1*b1
    }

    void add_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2x2_4x2x1w& a, const fp2x2_4x2x1w& b)
    {
        UFpPack8x1::add(r, a, b);
    }

    void sub_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2x2_4x2x1w& a, const fp2x2_4x2x1w& b)
    {
        UFpPack8x1::sub(r, a, b);
    }

    void sax2_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2x2_4x2x1w& a, const fp2x2_4x2x1w& b)
    {
        UFpPack8x1::ssaax2(r, a, b);
    }

    // r0 = a0 - a1
    // r1 = a0 + a1
    void mul_by_u_plus_1_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2x2_4x2x1w& a)
    {
        fp2x2_4x2x1w tt0;

        // a = A1 | A0 at Fp layer 
        shuf_01_dl(tt0.data(), a.data());     //     A0 |    A1
        UFpPack8x1::asx4(r, a, tt0);            //  A0+A1 | A0-A1
    }

    // r0 = a0*b0 - a1*b1
    // r1 = a0*b1 + a1*b0
    // schoolbook
    void mul_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b)
    {
        fp2_4x2x1w t0, t1;
        fp2x2_4x2x1w tt0, tt1;

        // a = A1 | A0 at Fp layer 
        // b = B1 | B0 at Fp layer
        shuf_00(t0.data(), a.data());         //          A0 |          A0 
        FpPack8x1::mul(tt0, t0, b);            //       A0*B1 |       A0*B0
        shuf_11(t0.data(), a.data());         //          A1 |          A1 
        shuf_01(t1.data(), b.data());         //          B0 |          B1 
        FpPack8x1::mul(tt1, t0, t1);           //       A1*B0 |       A1*B1
        UFpPack8x1::asx4(r, tt0, tt1);          // A0*B1+A1*B0 | A0*B0-A1*B1
    }

    // r0 = (a0 + a1)*(a0 - a1)
    // r1 = 2*a0*a1
    void sqr_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2_4x2x1w& a)
    {
        fp2_4x2x1w t0, t1, t2;

        // a = A1 | A0 at Fp layer 
        shuf_00(t0.data(), a.data());         //      A0 |              A0
        shuf_01(t1.data(), a.data());         //      A0 |              A1 
        shuf_z1(t2.data(), a.data());         //       0 |              A1
        FpPack8x1::add(t0, t0, t1);              //    2*A0 |           A0+A1
        FpPack8x1::sub(t2, a, t2);               //      A1 |           A0-A1
        FpPack8x1::mul(r, t0, t2);             // 2*A0*A1 | (A0+A1)*(A0-A1) 
    }

    void redc_fp2x2_4x2x1w(fp2_4x2x1w& r, const fp2x2_4x2x1w& a)
    {
        UFpPack8x1::redc(r, a);
    }

    // r0 = a0*b0 - a1*b1
    // r1 = a0*b1 + a1*b0
    // schoolbook
    void mul_fp2x2_2x4x1w(fp2x2_2x4x1w& r, const fp2_2x4x1w& a, const fp2_2x4x1w& b)
    {
        fp2_2x4x1w t0, t1;
        fp2x2_2x4x1w tt0, tt1;

        // a = A1 | A0 | ... | ... at Fp layer
        // b = B1 | B0 | ... | ... at Fp layer
        perm_3322(t0.data(), a.data());       //          A1 |          A1 |   A0 |   A0
        perm_2332(t1.data(), b.data());       //          B0 |          B1 |   B1 |   B0
        FpPack8x1::mul(tt0, t0, t1);           //        A1B0 |        A1B1 | A0B1 | A0B0
        perm_10zz_dl(tt1.data(), tt0.data()); //       A0*B1 |       A0*B0 |    0 |    0 
        UFpPack8x1::asx4(r, tt1, tt0);          // A0*B1+A1*B0 | A0*B0-A1*B1 |  ... |  ...
    }

    void redc_fp2x2_2x2x2w(fp2_2x2x2w& r, const fp2x2_2x2x2w& a)
    {
        UFpPack4x2::redc(r, a);
    }

    void add_fp2x2_2x2x2w(fp2x2_2x2x2w& r, const fp2x2_2x2x2w& a, const fp2x2_2x2x2w& b)
    {
        UFpPack4x2::add(r, a, b);
    }

    void sub_fp2x2_2x2x2w(fp2x2_2x2x2w& r, const fp2x2_2x2x2w& a, const fp2x2_2x2x2w& b)
    {
        UFpPack4x2::sub(r, a, b);
    }

    // r0 = a0 - a1
    // r1 = a0 + a1
    void mul_by_u_plus_1_fp2x2_2x2x2w(fp2x2_2x2x2w& r, const fp2x2_2x2x2w& a)
    {
        fp2x2_2x2x2w t0;

        // a = A1 | A0 at Fp layer 
        perm_1032(t0.data(), a.data());           //     A0 |    A1
        perm_1032_hl(&t0[NWORDS], &a[NWORDS]);    //     A0 |    A1
        UFpPack4x2::asx2(r, a, t0);                 //  A0+A1 | A0-A1
    }

    // r0 = a0*b0 - a1*b1
    // r1 = a0*b1 + a1*b0
    // schoolbook
    void mul_fp2x2_2x2x2w(fp2x2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b)
    {
        fp2_2x2x2w t0, t1;
        fp2x2_2x2x2w tt0, tt1;

        // a = A1 | A0 at Fp layer 
        // b = B1 | B0 at Fp layer
        perm_1010_hl(t0.data(), a.data());    //          A0 |          A0
        FpPack4x2::mul(tt0, t0, b);            //       A0*B1 |       A0*B0
        perm_3232_hl(t0.data(), a.data());    //          A1 |          A1 
        perm_1032_hl(t1.data(), b.data());    //          B0 |          B1
        FpPack4x2::mul(tt1, t0, t1);           //       A1*B0 |       A1*B1
        UFpPack4x2::asx2(r, tt0, tt1);          // A0*B1+A1*B0 | A0*B0-A1*B1  
    }

    // r0 = (a0 + a1)*(a0 - a1)
    // r1 = 2*a0*a1
    void sqr_fp2x2_2x2x2w(fp2x2_2x2x2w& r, const fp2_2x2x2w& a)
    {
        fp2_2x2x2w t0, t1, t2;
        // a = A1 | A0 at Fp layer
        perm_1010_hl(t0.data(), a.data());    //      A0 |              A0
        perm_1032_hl(t1.data(), a.data());    //      A0 |              A1 
        perm_zz32_hl(t2.data(), a.data());    //       0 |              A1
        FpPack4x2::add(t0, t0, t1);              //    2*A0 |           A0+A1
        FpPack4x2::sub(t2, a, t2);               //      A1 |           A0-A1
        FpPack4x2::mul(r, t0, t2);             // 2*A0*A1 | (A0+A1)*(A0-A1)
    }

    // r0 = a0*b0 - a1*b1
    // r1 = a0*b1 + a1*b0
    // schoolbook
    void mul_fp2x2_1x4x2w(fp2x2_1x4x2w& r, const fp2_1x4x2w& a, const fp2_1x4x2w& b)
    {
        fp2_1x4x2w t0, t1;
        fp2x2_1x4x2w tt0, tt1;
        const cc512i m0 = VSET(7, 6, 7, 6, 5, 4, 5, 4);
        const cc512i m1 = VSET(5, 4, 7, 6, 7, 6, 5, 4);
        const cc512i m2 = VSET(3, 2, 1, 0, 7, 6, 5, 4);

        // a = A1 | A0 | ... | ... at Fp layer
        // b = B1 | B0 | ... | ... at Fp layer
        perm_var_hl(t0.data(), a.data(), m0);         //          A1 |          A1 |    A0 |    A0
        perm_var_hl(t1.data(), b.data(), m1);         //          B0 |          B1 |    B1 |    B0
        FpPack4x2::mul(tt0, t0, t1);                   //       A1*B0 |       A1*B1 | A0*B1 | A0*B0
        perm_var_vl(tt1.data(), tt0.data(), m2);      //       A0*B1 |       A0*B0 |   ... |   ... 
        UFpPack4x2::asx2(r, tt1, tt0);                  // A0*B1+A1*B0 | A0*B0-A1*B1 |   ... |   ...
    }

} // namespace blst

#endif
