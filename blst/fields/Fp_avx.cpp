#include "blst/fields/Fp_avx.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
const u64 P48[NWORDS] = {
    0xFFFFFFFFAAABULL, 0xB153FFFFB9FEULL, 0xF6241EABFFFEULL, 0x6730D2A0F6B0ULL,
    0x4B84F38512BFULL, 0x434BACD76477ULL, 0xE69A4B1BA7B6ULL, 0x1A0111EA397FULL,
};

void conv_64to48_fp(u64 *r, const u64 *a)
{
  r[0] =                  a[0]         & BMASK;
  r[1] = ((a[0] >> 48) | (a[1] << 16)) & BMASK;
  r[2] = ((a[1] >> 32) | (a[2] << 32)) & BMASK;
  r[3] =  (a[2] >> 16)                 & BMASK;
  r[4] =                  a[3]         & BMASK;
  r[5] = ((a[3] >> 48) | (a[4] << 16)) & BMASK;
  r[6] = ((a[4] >> 32) | (a[5] << 32)) & BMASK;
  r[7] =  (a[5] >> 16)                 & BMASK;
}

void conv_64to48_fp_8x1w(fp_8x1w& r, const cc512i* BLST_restrict a)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2];
  const cc512i a3 = a[3], a4 = a[4], a5 = a[5];
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7;

  r0 = VAND(                       a0,       bmask);
  r1 = VAND(VOR(VSHR(a0, 48), VSHL(a1, 16)), bmask);
  r2 = VAND(VOR(VSHR(a1, 32), VSHL(a2, 32)), bmask);
  r3 = VAND(    VSHR(a2, 16),                bmask);
  r4 = VAND(                       a3,       bmask);
  r5 = VAND(VOR(VSHR(a3, 48), VSHL(a4, 16)), bmask);
  r6 = VAND(VOR(VSHR(a4, 32), VSHL(a5, 32)), bmask);
  r7 = VAND(    VSHR(a5, 16),                bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
}

void conv_64to48_fp_4x2w(fp_4x2w& r, const cc512i* BLST_restrict a)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2];
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2, r3;

  r0 = VAND(                       a0,       bmask);
  r1 = VAND(VOR(VSHR(a0, 48), VSHL(a1, 16)), bmask);
  r2 = VAND(VOR(VSHR(a1, 32), VSHL(a2, 32)), bmask);
  r3 = VAND(    VSHR(a2, 16),                bmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
}

void conv_64to48_mpi(u64 *r, const u64 *a, int rlen, int alen)
{
  int i, j, shr_pos, shl_pos;
  u64 word, temp;

  i = j = 0;
  shr_pos = 64; shl_pos = 0;
  temp = 0;
  while ((i < rlen) && (j < alen)) {
    word = ((temp >> shr_pos) | (a[j] << shl_pos));
    r[i] = (word & BMASK);
    shr_pos -= 16, shl_pos += 16;
    if ((shr_pos > 0) && (shl_pos < 64)) temp = a[j++];
    if (shr_pos <= 0) shr_pos += 64;
    if (shl_pos >= 64) shl_pos -= 64;
    // Any shift past 63 is undefined!
    if (shr_pos == 64) temp = 0;
    i++;
  }
  if (i < rlen) r[i++] = ((temp >> shr_pos) & BMASK);
  for (; i < rlen; i++) r[i] = 0;
}

void conv_48to64_fp(u64 *r, const u64 *a)
{
  r[0] = (a[1] << 48) | a[0];
  r[1] = (a[2] << 32) | a[1] >> 16;
  r[2] = (a[3] << 16) | a[2] >> 32;
  r[3] = (a[5] << 48) | a[4];
  r[4] = (a[6] << 32) | a[5] >> 16;
  r[5] = (a[7] << 16) | a[6] >> 32;  
}

void conv_48to64_fp_8x1w(cc512i* BLST_restrict r, const fp_8x1w& a)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  cc512i r0, r1, r2, r3, r4, r5;

  r0 = VOR(VSHL(a1, 48),      a0     );
  r1 = VOR(VSHL(a2, 32), VSHR(a1, 16));
  r2 = VOR(VSHL(a3, 16), VSHR(a2, 32));
  r3 = VOR(VSHL(a5, 48),      a4     );
  r4 = VOR(VSHL(a6, 32), VSHR(a5, 16));
  r5 = VOR(VSHL(a7, 16), VSHR(a6, 32));

  r[0] = r0; r[1] = r1; r[2] = r2; 
  r[3] = r3; r[4] = r4; r[5] = r5; 
}

void conv_48to64_fp_4x2w(cc512i* BLST_restrict r, const fp_4x2w& a)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  cc512i r0, r1, r2;

  r0 = VOR(VSHL(a1, 48),      a0     );
  r1 = VOR(VSHL(a2, 32), VSHR(a1, 16));
  r2 = VOR(VSHL(a3, 16), VSHR(a2, 32));

  r[0] = r0; r[1] = r1; r[2] = r2; 
}

void conv_48to64_mpi(u64 *r, const u64 *a, int rlen, int alen)
{
  int i, j, bits_in_word, bits_to_shift;
  u64 word;

  i = j = 0;
  bits_in_word = bits_to_shift = 0;
  word = 0;
  while ((i < rlen) && (j < alen)) {
    word |= (a[j] << bits_in_word);
    bits_to_shift = (64 - bits_in_word);
    bits_in_word += 48;
    if (bits_in_word >= 64) {
      r[i++] = word;
      word = ((bits_to_shift > 0) ? (a[j] >> bits_to_shift) : 0);
      bits_in_word = ((bits_to_shift > 0) ? (48 - bits_to_shift) : 0);
    }
    j++;
  }
  if (i < rlen) r[i++] = word;
  for (; i < rlen; i++) r[i] = 0;
}

void carryp_mpi48(u64 *a)
{
  int i;

  for (i = 0; i < NWORDS-1; i++) {
    a[i+1] += a[i]>>BRADIX;
    a[i] &= BMASK;
  }
}

void carryp_fp_4x2w(fp_4x2w& a)
{
  cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i bmask = VSET1(BMASK);

  a1 = VADD(a1, VSRA(a0, BRADIX)); a0 = VAND(a0, bmask);
  a2 = VADD(a2, VSRA(a1, BRADIX)); a1 = VAND(a1, bmask);
  a3 = VADD(a3, VSRA(a2, BRADIX)); a2 = VAND(a2, bmask);

  a[0] = a0; a[1] = a1; a[2] = a2; a[3] = a3;
}

void get_channel_8x1w(u64* BLST_restrict r, const fp_8x1w& a, const int ch) 
{
  int i;

  for(i = 0; i < NWORDS; i++) {
    r[i] = ((u64 *)&a[i])[ch];
  }
}

void get_channel_4x2w(u64* BLST_restrict r, const fp_4x2w& a, const int ch) 
{
  int i;

  for(i = 0; i < VWORDS; i++) {
    r[i] = ((u64 *)&a[i])[ch];
    r[i+VWORDS] = ((u64 *)&a[i])[ch+1];
  }
}

void add_fp_8x1w(fp_8x1w& r, const fp_8x1w& a, const fp_8x1w& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask, t0;

  // r = a + b
  r0 = VADD(a0, b0); r1 = VADD(a1, b1); r2 = VADD(a2, b2); r3 = VADD(a3, b3);
  r4 = VADD(a4, b4); r5 = VADD(a5, b5); r6 = VADD(a6, b6); r7 = VADD(a7, b7);

  // r = r - p
  r0 = VSUB(r0, p0); r1 = VSUB(r1, p1); r2 = VSUB(r2, p2); r3 = VSUB(r3, p3);
  r4 = VSUB(r4, p4); r5 = VSUB(r5, p5); r6 = VSUB(r6, p6); r7 = VSUB(r7, p7);

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

void sub_fp_8x1w(fp_8x1w& r, const fp_8x1w& a, const fp_8x1w& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask, t0;

  // r = a - b
  r0 = VSUB(a0, b0); r1 = VSUB(a1, b1); r2 = VSUB(a2, b2); r3 = VSUB(a3, b3);
  r4 = VSUB(a4, b4); r5 = VSUB(a5, b5); r6 = VSUB(a6, b6); r7 = VSUB(a7, b7);

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

void asx4_fp_8x1w(fp_8x1w& r, const fp_8x1w& a, const fp_8x1w& b)
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

  // r = H+P | G | F+N | E | D+L | C | B+J | A
  r0  = VMADD(a0, 0xAA, a0, b0); r1  = VMADD(a1, 0xAA, a1, b1);
  r2  = VMADD(a2, 0xAA, a2, b2); r3  = VMADD(a3, 0xAA, a3, b3);
  r4  = VMADD(a4, 0xAA, a4, b4); r5  = VMADD(a5, 0xAA, a5, b5);
  r6  = VMADD(a6, 0xAA, a6, b6); r7  = VMADD(a7, 0xAA, a7, b7);

  // t = p | O | p | M | p | K | p | I
  t0 = VMBLEND(0xAA, b0, p0); t1 = VMBLEND(0xAA, b1, p1);
  t2 = VMBLEND(0xAA, b2, p2); t3 = VMBLEND(0xAA, b3, p3);
  t4 = VMBLEND(0xAA, b4, p4); t5 = VMBLEND(0xAA, b5, p5);
  t6 = VMBLEND(0xAA, b6, p6); t7 = VMBLEND(0xAA, b7, p7); 

  // r = H+P-p | G-O | F+N-p | E-M | D+L-p | C-K | B+J-p | A-I
  r0 = VSUB(r0, t0); r1 = VSUB(r1, t1);
  r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
  r4 = VSUB(r4, t4); r5 = VSUB(r5, t5);
  r6 = VSUB(r6, t6); r7 = VSUB(r7, t7);

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

void asax6_fp_8x1w(fp_8x1w& r, const fp_8x1w& a, const fp_8x1w& b)
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

  // r = H+P | G | F+N | E+M | D+L | C+K | B+J | A+I
  r0  = VMADD(a0, 0xBF, a0, b0); r1  = VMADD(a1, 0xBF, a1, b1);
  r2  = VMADD(a2, 0xBF, a2, b2); r3  = VMADD(a3, 0xBF, a3, b3);
  r4  = VMADD(a4, 0xBF, a4, b4); r5  = VMADD(a5, 0xBF, a5, b5);
  r6  = VMADD(a6, 0xBF, a6, b6); r7  = VMADD(a7, 0xBF, a7, b7);

  // t = p | O | p | p | p | p | p | p
  t0 = VMBLEND(0xBF, b0, p0); t1 = VMBLEND(0xBF, b1, p1);
  t2 = VMBLEND(0xBF, b2, p2); t3 = VMBLEND(0xBF, b3, p3);
  t4 = VMBLEND(0xBF, b4, p4); t5 = VMBLEND(0xBF, b5, p5);
  t6 = VMBLEND(0xBF, b6, p6); t7 = VMBLEND(0xBF, b7, p7); 

  // r = H+P-p | G-O | F+N-p | E+M-p | D+L-p | C+K-p | B+J-p | A+I-p
  r0 = VSUB(r0, t0); r1 = VSUB(r1, t1);
  r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
  r4 = VSUB(r4, t4); r5 = VSUB(r5, t5);
  r6 = VSUB(r6, t6); r7 = VSUB(r7, t7);

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

void assx6_fp_8x1w(fp_8x1w& r, const fp_8x1w& a, const fp_8x1w& b)
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

  // r = H+P | G | F | E | D | C | B | A
  r0  = VMADD(a0, 0x80, a0, b0); r1  = VMADD(a1, 0x80, a1, b1);
  r2  = VMADD(a2, 0x80, a2, b2); r3  = VMADD(a3, 0x80, a3, b3);
  r4  = VMADD(a4, 0x80, a4, b4); r5  = VMADD(a5, 0x80, a5, b5);
  r6  = VMADD(a6, 0x80, a6, b6); r7  = VMADD(a7, 0x80, a7, b7);

  // t = p | O | N | M | L | K | J | I
  t0 = VMBLEND(0x80, b0, p0); t1 = VMBLEND(0x80, b1, p1);
  t2 = VMBLEND(0x80, b2, p2); t3 = VMBLEND(0x80, b3, p3);
  t4 = VMBLEND(0x80, b4, p4); t5 = VMBLEND(0x80, b5, p5);
  t6 = VMBLEND(0x80, b6, p6); t7 = VMBLEND(0x80, b7, p7); 

  // r = H+P-p | G-O | F-N | E-M | D-L | C-K | B-J | A-I
  r0 = VSUB(r0, t0); r1 = VSUB(r1, t1);
  r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
  r4 = VSUB(r4, t4); r5 = VSUB(r5, t5);
  r6 = VSUB(r6, t6); r7 = VSUB(r7, t7);

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

void add_fp_4x2w(fp_4x2w& r, const fp_4x2w& a, const fp_4x2w& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask = VSET1(BMASK); 
  cc512i r0, r1, r2, r3, smask, t0;

  // r = a + b
  r0 = VADD(a0, b0); r1 = VADD(a1, b1); r2 = VADD(a2, b2); r3 = VADD(a3, b3);

  // r = r - p
  r0 = VSUB(r0, p0); r1 = VSUB(r1, p1); r2 = VSUB(r2, p2); r3 = VSUB(r3, p3);

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

void sub_fp_4x2w(fp_4x2w& r, const fp_4x2w& a, const fp_4x2w& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask = VSET1(BMASK); 
  cc512i r0, r1, r2, r3, smask, t0;

  // r = a - b
  r0 = VSUB(a0, b0); r1 = VSUB(a1, b1); r2 = VSUB(a2, b2); r3 = VSUB(a3, b3);

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

void asx2_fp_4x2w(fp_4x2w& r, const fp_4x2w& a, const fp_4x2w& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask = VSET1(BMASK); 
  cc512i r0, r1, r2, r3, smask, t0, t1, t2, t3; 

  // r = D'+H' | D+H | C' | C | B'+F' | B+F | A' | A 
  r0 = VMADD(a0, 0xCC, a0, b0); r1 = VMADD(a1, 0xCC, a1, b1);
  r2 = VMADD(a2, 0xCC, a2, b2); r3 = VMADD(a3, 0xCC, a3, b3);

  // t = p' | p | G' | G | p' | p | E' | E
  t0 = VMBLEND(0xCC, b0, p0); t1 = VMBLEND(0xCC, b1, p1);
  t2 = VMBLEND(0xCC, b2, p2); t3 = VMBLEND(0xCC, b3, p3);

  // r = D'+H'-p | D+H-p | C'-G' | C-G | B'+F'-p' | B+F-p | A'-E' | A-E
  r0 = VSUB(r0, t0); r1 = VSUB(r1, t1);
  r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);

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

void add_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask, t0;

  // r = a + b
  r0  = VADD(a0 , b0 ); r1  = VADD(a1 , b1 );
  r2  = VADD(a2 , b2 ); r3  = VADD(a3 , b3 );
  r4  = VADD(a4 , b4 ); r5  = VADD(a5 , b5 );
  r6  = VADD(a6 , b6 ); r7  = VADD(a7 , b7 );
  r8  = VADD(a8 , b8 ); r9  = VADD(a9 , b9 );
  r10 = VADD(a10, b10); r11 = VADD(a11, b11);
  r12 = VADD(a12, b12); r13 = VADD(a13, b13);
  r14 = VADD(a14, b14); r15 = VADD(a15, b15);

  // r = r - p*2^384
  r8  = VSUB(r8 , p0); r9  = VSUB(r9 , p1);
  r10 = VSUB(r10, p2); r11 = VSUB(r11, p3);
  r12 = VSUB(r12, p4); r13 = VSUB(r13, p5);
  r14 = VSUB(r14, p6); r15 = VSUB(r15, p7);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void sub_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask, t0;

  // r = a - b
  r0  = VSUB(a0 , b0 ); r1  = VSUB(a1 , b1 );
  r2  = VSUB(a2 , b2 ); r3  = VSUB(a3 , b3 );
  r4  = VSUB(a4 , b4 ); r5  = VSUB(a5 , b5 );
  r6  = VSUB(a6 , b6 ); r7  = VSUB(a7 , b7 );
  r8  = VSUB(a8 , b8 ); r9  = VSUB(a9 , b9 );
  r10 = VSUB(a10, b10); r11 = VSUB(a11, b11);
  r12 = VSUB(a12, b12); r13 = VSUB(a13, b13);
  r14 = VSUB(a14, b14); r15 = VSUB(a15, b15);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void aaasx2_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask;
  cc512i t0, t1, t2 , t3 , t4 , t5 , t6 , t7 ;

  // r = H+P | G+O | F+N | E | D+L | C+K | B+N | A
  r0  = VMADD(a0 , 0xEE, a0 , b0 ); r1  = VMADD(a1 , 0xEE, a1 , b1 );
  r2  = VMADD(a2 , 0xEE, a2 , b2 ); r3  = VMADD(a3 , 0xEE, a3 , b3 );
  r4  = VMADD(a4 , 0xEE, a4 , b4 ); r5  = VMADD(a5 , 0xEE, a5 , b5 );
  r6  = VMADD(a6 , 0xEE, a6 , b6 ); r7  = VMADD(a7 , 0xEE, a7 , b7 );
  r8  = VMADD(a8 , 0xEE, a8 , b8 ); r9  = VMADD(a9 , 0xEE, a9 , b9 );
  r10 = VMADD(a10, 0xEE, a10, b10); r11 = VMADD(a11, 0xEE, a11, b11);
  r12 = VMADD(a12, 0xEE, a12, b12); r13 = VMADD(a13, 0xEE, a13, b13);
  r14 = VMADD(a14, 0xEE, a14, b14); r15 = VMADD(a15, 0xEE, a15, b15);

  // t = p | p | p | M | p | p | p | I
  t0 = VMBLEND(0xEE, b8 , p0); t1 = VMBLEND(0xEE, b9 , p1);
  t2 = VMBLEND(0xEE, b10, p2); t3 = VMBLEND(0xEE, b11, p3);
  t4 = VMBLEND(0xEE, b12, p4); t5 = VMBLEND(0xEE, b13, p5);
  t6 = VMBLEND(0xEE, b14, p6); t7 = VMBLEND(0xEE, b15, p7); 

  // r = H+P | G+O | F+N | E-M | D+L | C+K | B+J | A-I
  r0 = VMSUB(r0, 0x11, r0, b0); r1 = VMSUB(r1, 0x11, r1, b1);
  r2 = VMSUB(r2, 0x11, r2, b2); r3 = VMSUB(r3, 0x11, r3, b3);
  r4 = VMSUB(r4, 0x11, r4, b4); r5 = VMSUB(r5, 0x11, r5, b5);
  r6 = VMSUB(r6, 0x11, r6, b6); r7 = VMSUB(r7, 0x11, r7, b7);
  // r = H+P-p | G+O-p | F+N-p | E-M | D+L-p | C+K-p | B+J-p | A-I
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);
  r12 = VSUB(r12, t4); r13 = VSUB(r13, t5);
  r14 = VSUB(r14, t6); r15 = VSUB(r15, t7);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void aassx2_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask;
  cc512i t0, t1, t2 , t3 , t4 , t5 , t6 , t7 ;

  // r = H+P | G+O | F | E | D+L | C+K | B | A
  r0  = VMADD(a0 , 0xCC, a0 , b0 ); r1  = VMADD(a1 , 0xCC, a1 , b1 );
  r2  = VMADD(a2 , 0xCC, a2 , b2 ); r3  = VMADD(a3 , 0xCC, a3 , b3 );
  r4  = VMADD(a4 , 0xCC, a4 , b4 ); r5  = VMADD(a5 , 0xCC, a5 , b5 );
  r6  = VMADD(a6 , 0xCC, a6 , b6 ); r7  = VMADD(a7 , 0xCC, a7 , b7 );
  r8  = VMADD(a8 , 0xCC, a8 , b8 ); r9  = VMADD(a9 , 0xCC, a9 , b9 );
  r10 = VMADD(a10, 0xCC, a10, b10); r11 = VMADD(a11, 0xCC, a11, b11);
  r12 = VMADD(a12, 0xCC, a12, b12); r13 = VMADD(a13, 0xCC, a13, b13);
  r14 = VMADD(a14, 0xCC, a14, b14); r15 = VMADD(a15, 0xCC, a15, b15);

  // t = p | p | N | M | p | p | J | I
  t0 = VMBLEND(0xCC, b8 , p0); t1 = VMBLEND(0xCC, b9 , p1);
  t2 = VMBLEND(0xCC, b10, p2); t3 = VMBLEND(0xCC, b11, p3);
  t4 = VMBLEND(0xCC, b12, p4); t5 = VMBLEND(0xCC, b13, p5);
  t6 = VMBLEND(0xCC, b14, p6); t7 = VMBLEND(0xCC, b15, p7); 

  // r = H+P | G+O | F-N | E-M | D+L | C+K | B-J | A-I
  r0 = VMSUB(r0, 0x33, r0, b0); r1 = VMSUB(r1, 0x33, r1, b1);
  r2 = VMSUB(r2, 0x33, r2, b2); r3 = VMSUB(r3, 0x33, r3, b3);
  r4 = VMSUB(r4, 0x33, r4, b4); r5 = VMSUB(r5, 0x33, r5, b5);
  r6 = VMSUB(r6, 0x33, r6, b6); r7 = VMSUB(r7, 0x33, r7, b7);
  // r = H+P-p | G+O-p | F-N | E-M | D+L-p | C+K-p | B-J | A-I
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);
  r12 = VSUB(r12, t4); r13 = VSUB(r13, t5);
  r14 = VSUB(r14, t6); r15 = VSUB(r15, t7);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void asx4_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask;
  cc512i t0, t1, t2 , t3 , t4 , t5 , t6 , t7 ;

  // r = H+P | G | F+N | E | D+L | C | B+J | A
  r0  = VMADD(a0 , 0xAA, a0 , b0 ); r1  = VMADD(a1 , 0xAA, a1 , b1 );
  r2  = VMADD(a2 , 0xAA, a2 , b2 ); r3  = VMADD(a3 , 0xAA, a3 , b3 );
  r4  = VMADD(a4 , 0xAA, a4 , b4 ); r5  = VMADD(a5 , 0xAA, a5 , b5 );
  r6  = VMADD(a6 , 0xAA, a6 , b6 ); r7  = VMADD(a7 , 0xAA, a7 , b7 );
  r8  = VMADD(a8 , 0xAA, a8 , b8 ); r9  = VMADD(a9 , 0xAA, a9 , b9 );
  r10 = VMADD(a10, 0xAA, a10, b10); r11 = VMADD(a11, 0xAA, a11, b11);
  r12 = VMADD(a12, 0xAA, a12, b12); r13 = VMADD(a13, 0xAA, a13, b13);
  r14 = VMADD(a14, 0xAA, a14, b14); r15 = VMADD(a15, 0xAA, a15, b15);

  // t = p | O | p | M | p | K | p | I
  t0 = VMBLEND(0xAA, b8 , p0); t1 = VMBLEND(0xAA, b9 , p1);
  t2 = VMBLEND(0xAA, b10, p2); t3 = VMBLEND(0xAA, b11, p3);
  t4 = VMBLEND(0xAA, b12, p4); t5 = VMBLEND(0xAA, b13, p5);
  t6 = VMBLEND(0xAA, b14, p6); t7 = VMBLEND(0xAA, b15, p7); 

  // r = H+P | G-O | F+N | E-M | D+L | C-K | B+J | A-I
  r0 = VMSUB(r0, 0x55, r0, b0); r1 = VMSUB(r1, 0x55, r1, b1);
  r2 = VMSUB(r2, 0x55, r2, b2); r3 = VMSUB(r3, 0x55, r3, b3);
  r4 = VMSUB(r4, 0x55, r4, b4); r5 = VMSUB(r5, 0x55, r5, b5);
  r6 = VMSUB(r6, 0x55, r6, b6); r7 = VMSUB(r7, 0x55, r7, b7);
  // r = H+P-p | G-O | F+N-p | E-M | D+L-p | C-K | B+J-p | A-I
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);
  r12 = VSUB(r12, t4); r13 = VSUB(r13, t5);
  r14 = VSUB(r14, t6); r15 = VSUB(r15, t7);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void ssaax2_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask;
  cc512i t0, t1, t2 , t3 , t4 , t5 , t6 , t7 ;

  // r = H | G | F+N | E+M | D | C | B+J | A+I
  r0  = VMADD(a0 , 0x33, a0 , b0 ); r1  = VMADD(a1 , 0x33, a1 , b1 );
  r2  = VMADD(a2 , 0x33, a2 , b2 ); r3  = VMADD(a3 , 0x33, a3 , b3 );
  r4  = VMADD(a4 , 0x33, a4 , b4 ); r5  = VMADD(a5 , 0x33, a5 , b5 );
  r6  = VMADD(a6 , 0x33, a6 , b6 ); r7  = VMADD(a7 , 0x33, a7 , b7 );
  r8  = VMADD(a8 , 0x33, a8 , b8 ); r9  = VMADD(a9 , 0x33, a9 , b9 );
  r10 = VMADD(a10, 0x33, a10, b10); r11 = VMADD(a11, 0x33, a11, b11);
  r12 = VMADD(a12, 0x33, a12, b12); r13 = VMADD(a13, 0x33, a13, b13);
  r14 = VMADD(a14, 0x33, a14, b14); r15 = VMADD(a15, 0x33, a15, b15);

  // t = P | O | p | p | L | K | p | p
  t0 = VMBLEND(0x33, b8 , p0); t1 = VMBLEND(0x33, b9 , p1);
  t2 = VMBLEND(0x33, b10, p2); t3 = VMBLEND(0x33, b11, p3);
  t4 = VMBLEND(0x33, b12, p4); t5 = VMBLEND(0x33, b13, p5);
  t6 = VMBLEND(0x33, b14, p6); t7 = VMBLEND(0x33, b15, p7); 

  // r = H-P | G-O | F+N | E+M | D-L | C-K | B+J | A+I
  r0 = VMSUB(r0, 0xCC, r0, b0); r1 = VMSUB(r1, 0xCC, r1, b1);
  r2 = VMSUB(r2, 0xCC, r2, b2); r3 = VMSUB(r3, 0xCC, r3, b3);
  r4 = VMSUB(r4, 0xCC, r4, b4); r5 = VMSUB(r5, 0xCC, r5, b5);
  r6 = VMSUB(r6, 0xCC, r6, b6); r7 = VMSUB(r7, 0xCC, r7, b7);
  // r = H-P | G-O | F+N-p | E+M-p | D-L | C-K | B+J-p | A+I-p
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);
  r12 = VSUB(r12, t4); r13 = VSUB(r13, t5);
  r14 = VSUB(r14, t6); r15 = VSUB(r15, t7);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void ssasx2_fpx2_8x1w(fpx2_8x1w& r, const fpx2_8x1w& a, const fpx2_8x1w& b)
{
  const cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  const cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  const cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  const cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  const cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  const cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  const cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK);
  cc512i r0, r1, r2 , r3 , r4 , r5 , r6 , r7 ;
  cc512i r8, r9, r10, r11, r12, r13, r14, r15, smask;
  cc512i t0, t1, t2 , t3 , t4 , t5 , t6 , t7 ;

  // r = H | G | F+N | E | D | C | B+N | A
  r0  = VMADD(a0 , 0x22, a0 , b0 ); r1  = VMADD(a1 , 0x22, a1 , b1 );
  r2  = VMADD(a2 , 0x22, a2 , b2 ); r3  = VMADD(a3 , 0x22, a3 , b3 );
  r4  = VMADD(a4 , 0x22, a4 , b4 ); r5  = VMADD(a5 , 0x22, a5 , b5 );
  r6  = VMADD(a6 , 0x22, a6 , b6 ); r7  = VMADD(a7 , 0x22, a7 , b7 );
  r8  = VMADD(a8 , 0x22, a8 , b8 ); r9  = VMADD(a9 , 0x22, a9 , b9 );
  r10 = VMADD(a10, 0x22, a10, b10); r11 = VMADD(a11, 0x22, a11, b11);
  r12 = VMADD(a12, 0x22, a12, b12); r13 = VMADD(a13, 0x22, a13, b13);
  r14 = VMADD(a14, 0x22, a14, b14); r15 = VMADD(a15, 0x22, a15, b15);

  // t = P | O | p | M | L | K | p | I
  t0 = VMBLEND(0x22, b8 , p0); t1 = VMBLEND(0x22, b9 , p1);
  t2 = VMBLEND(0x22, b10, p2); t3 = VMBLEND(0x22, b11, p3);
  t4 = VMBLEND(0x22, b12, p4); t5 = VMBLEND(0x22, b13, p5);
  t6 = VMBLEND(0x22, b14, p6); t7 = VMBLEND(0x22, b15, p7); 

  // r = H-P | G-O | F+N | E-M | D-L | C-K | B+J | A-I
  r0 = VMSUB(r0, 0xDD, r0, b0); r1 = VMSUB(r1, 0xDD, r1, b1);
  r2 = VMSUB(r2, 0xDD, r2, b2); r3 = VMSUB(r3, 0xDD, r3, b3);
  r4 = VMSUB(r4, 0xDD, r4, b4); r5 = VMSUB(r5, 0xDD, r5, b5);
  r6 = VMSUB(r6, 0xDD, r6, b6); r7 = VMSUB(r7, 0xDD, r7, b7);
  // r = H-P | G-O | F+N-p | E-M | D-L | C-K | B+J-p | A-I
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);
  r12 = VSUB(r12, t4); r13 = VSUB(r13, t5);
  r14 = VSUB(r14, t6); r15 = VSUB(r15, t7);

  // get sign mask
  r1 = VADD(r1 , VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
  r2 = VADD(r2 , VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
  r3 = VADD(r3 , VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
  r4 = VADD(r4 , VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
  r5 = VADD(r5 , VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
  r6 = VADD(r6 , VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
  r7 = VADD(r7 , VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
  r8 = VADD(r8 , VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
  t0 = VADD(r9 , VSRA(r8, BRADIX));
  t0 = VADD(r10, VSRA(t0, BRADIX));
  t0 = VADD(r11, VSRA(t0, BRADIX));
  t0 = VADD(r12, VSRA(t0, BRADIX));
  t0 = VADD(r13, VSRA(t0, BRADIX));
  t0 = VADD(r14, VSRA(t0, BRADIX));
  t0 = VADD(r15, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  // r = r + (p & smask)*2^384
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
  r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
  r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

  // carry propagation 
  r9  = VADD(r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask);
  r10 = VADD(r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
                                      r15 = VAND(r15, bmask);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
}

void add_fpx2_4x2w(fpx2_4x2w& r, const fpx2_4x2w& a, const fpx2_4x2w& b)
{
  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask0 = VSET(0, BMASK, 0, BMASK, 0, BMASK, 0, BMASK);
  const cc512i bmask1 = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11; 
  cc512i t0, t1, t2, t3, smask;

  // r = a + b
  r0  = VADD(a0 , b0 ); r1  = VADD(a1 , b1 );
  r2  = VADD(a2 , b2 ); r3  = VADD(a3 , b3 );
  r4  = VADD(a4 , b4 ); r5  = VADD(a5 , b5 );
  r6  = VADD(a6 , b6 ); r7  = VADD(a7 , b7 );
  r8  = VADD(a8 , b8 ); r9  = VADD(a9 , b9 );
  r10 = VADD(a10, b10); r11 = VADD(a11, b11);

  // r = r - p*2^384
  r8  = VSUB(r8 , p0); r9  = VSUB(r9 , p1);
  r10 = VSUB(r10, p2); r11 = VSUB(r11, p3);

  // get sign mask + carry propagation
  r1  = VMADD(r1 , 0x55, r1 , VSRA(r0, BRADIX)); 
  r4  = VMADD(r4 , 0x55, r4 , VSHUF<0x4E>(r0)); r0 = VAND(r0, bmask0);
  r2  = VMADD(r2 , 0x55, r2 , VSRA(r1, BRADIX)); 
  r5  = VMADD(r5 , 0x55, r5 , VSHUF<0x4E>(r1)); r1 = VAND(r1, bmask0);
  r3  = VMADD(r3 , 0x55, r3 , VSRA(r2, BRADIX)); 
  r6  = VMADD(r6 , 0x55, r6 , VSHUF<0x4E>(r2)); r2 = VAND(r2, bmask0);
  r4  = VMADD(r4 , 0x55, r4 , VSRA(r3, BRADIX)); 
  r7  = VMADD(r7 , 0x55, r7 , VSHUF<0x4E>(r3)); r3 = VAND(r3, bmask0);
  r5  = VMADD(r5 , 0x55, r5 , VSRA(r4, BRADIX)); 
  r8  = VMADD(r8 , 0x55, r8 , VSHUF<0x4E>(r4)); r4 = VAND(r4, bmask0);
  r6  = VMADD(r6 , 0x55, r6 , VSRA(r5, BRADIX)); 
  r9  = VMADD(r9 , 0x55, r9 , VSHUF<0x4E>(r5)); r5 = VAND(r5, bmask0);
  r7  = VMADD(r7 , 0x55, r7 , VSRA(r6, BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSHUF<0x4E>(r6)); r6 = VAND(r6, bmask0);
  r8  = VMADD(r8 , 0x55, r8 , VSRA(r7, BRADIX)); 
  r11 = VMADD(r11, 0x55, r11, VSHUF<0x4E>(r7)); r7 = VAND(r7, bmask0);
  t0  = VMADD(r9 , 0x55, r9 , VSRA(r8, BRADIX));
  t0  = VMADD(r10, 0x55, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0x55, r11, VSRA(t0, BRADIX));
  t0  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0  = VMADD(r9 , 0xAA, r9 , VSRA(t0, BRADIX));
  t0  = VMADD(r10, 0xAA, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0xAA, r11, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask)); 

  // carry propagation 
  r9  = VMADD(r9 , 0x55, r9 , VSRA(r8 , BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSRA(r9 , BRADIX));
  r11 = VMADD(r11, 0x55, r11, VSRA(r10, BRADIX));
  r8  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(r11, BRADIX)));
  r9  = VMADD(r9 , 0xAA, r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask1);
  r10 = VMADD(r10, 0xAA, r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask1);
  r11 = VMADD(r11, 0xAA, r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask1);
  r11 = VAND(r11, bmask1);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
}

void sub_fpx2_4x2w(fpx2_4x2w& r, const fpx2_4x2w& a, const fpx2_4x2w& b)
{
  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask0 = VSET(0, BMASK, 0, BMASK, 0, BMASK, 0, BMASK);
  const cc512i bmask1 = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11; 
  cc512i t0, t1, t2, t3, smask;

  // r = a - b
  r0  = VSUB(a0 , b0 ); r1  = VSUB(a1 , b1 );
  r2  = VSUB(a2 , b2 ); r3  = VSUB(a3 , b3 );
  r4  = VSUB(a4 , b4 ); r5  = VSUB(a5 , b5 );
  r6  = VSUB(a6 , b6 ); r7  = VSUB(a7 , b7 );
  r8  = VSUB(a8 , b8 ); r9  = VSUB(a9 , b9 );
  r10 = VSUB(a10, b10); r11 = VSUB(a11, b11);

  // get sign mask + carry propagation
  r1  = VMADD(r1 , 0x55, r1 , VSRA(r0, BRADIX)); 
  r4  = VMADD(r4 , 0x55, r4 , VSHUF<0x4E>(r0)); r0 = VAND(r0, bmask0);
  r2  = VMADD(r2 , 0x55, r2 , VSRA(r1, BRADIX)); 
  r5  = VMADD(r5 , 0x55, r5 , VSHUF<0x4E>(r1)); r1 = VAND(r1, bmask0);
  r3  = VMADD(r3 , 0x55, r3 , VSRA(r2, BRADIX)); 
  r6  = VMADD(r6 , 0x55, r6 , VSHUF<0x4E>(r2)); r2 = VAND(r2, bmask0);
  r4  = VMADD(r4 , 0x55, r4 , VSRA(r3, BRADIX)); 
  r7  = VMADD(r7 , 0x55, r7 , VSHUF<0x4E>(r3)); r3 = VAND(r3, bmask0);
  r5  = VMADD(r5 , 0x55, r5 , VSRA(r4, BRADIX)); 
  r8  = VMADD(r8 , 0x55, r8 , VSHUF<0x4E>(r4)); r4 = VAND(r4, bmask0);
  r6  = VMADD(r6 , 0x55, r6 , VSRA(r5, BRADIX)); 
  r9  = VMADD(r9 , 0x55, r9 , VSHUF<0x4E>(r5)); r5 = VAND(r5, bmask0);
  r7  = VMADD(r7 , 0x55, r7 , VSRA(r6, BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSHUF<0x4E>(r6)); r6 = VAND(r6, bmask0);
  r8  = VMADD(r8 , 0x55, r8 , VSRA(r7, BRADIX)); 
  r11 = VMADD(r11, 0x55, r11, VSHUF<0x4E>(r7)); r7 = VAND(r7, bmask0);
  t0  = VMADD(r9 , 0x55, r9 , VSRA(r8, BRADIX));
  t0  = VMADD(r10, 0x55, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0x55, r11, VSRA(t0, BRADIX));
  t0  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0  = VMADD(r9 , 0xAA, r9 , VSRA(t0, BRADIX));
  t0  = VMADD(r10, 0xAA, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0xAA, r11, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask)); 

  // carry propagation 
  r9  = VMADD(r9 , 0x55, r9 , VSRA(r8 , BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSRA(r9 , BRADIX));
  r11 = VMADD(r11, 0x55, r11, VSRA(r10, BRADIX));
  r8  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(r11, BRADIX)));
  r9  = VMADD(r9 , 0xAA, r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask1);
  r10 = VMADD(r10, 0xAA, r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask1);
  r11 = VMADD(r11, 0xAA, r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask1);
  r11 = VAND(r11, bmask1);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
}

void asx2_fpx2_4x2w(fpx2_4x2w& r, const fpx2_4x2w& a, const fpx2_4x2w& b)
{
  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask0 = VSET(0, BMASK, 0, BMASK, 0, BMASK, 0, BMASK);
  const cc512i bmask1 = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11; 
  cc512i t0, t1, t2, t3, smask;

  // r = D'+H' | D+H | C' | C | B'+F' | B+F | A' | A 
  r0  = VMADD(a0 , 0xCC, a0 , b0 ); r1  = VMADD(a1 , 0xCC, a1 , b1 );
  r2  = VMADD(a2 , 0xCC, a2 , b2 ); r3  = VMADD(a3 , 0xCC, a3 , b3 );
  r4  = VMADD(a4 , 0xCC, a4 , b4 ); r5  = VMADD(a5 , 0xCC, a5 , b5 );
  r6  = VMADD(a6 , 0xCC, a6 , b6 ); r7  = VMADD(a7 , 0xCC, a7 , b7 );
  r8  = VMADD(a8 , 0xCC, a8 , b8 ); r9  = VMADD(a9 , 0xCC, a9 , b9 );
  r10 = VMADD(a10, 0xCC, a10, b10); r11 = VMADD(a11, 0xCC, a11, b11);

  // t = p' | p | G' | G | p' | p | E' | E
  t0 = VMBLEND(0xCC, b8 , p0); t1 = VMBLEND(0xCC, b9 , p1);
  t2 = VMBLEND(0xCC, b10, p2); t3 = VMBLEND(0xCC, b11, p3);

  // r = D'+H' | D+H | C'-G' | C-G | B'+F' | B+F | A'-E' | A-E
  r0 = VMSUB(r0, 0x33, r0, b0); r1 = VMSUB(r1, 0x33, r1, b1);
  r2 = VMSUB(r2, 0x33, r2, b2); r3 = VMSUB(r3, 0x33, r3, b3);
  r4 = VMSUB(r4, 0x33, r4, b4); r5 = VMSUB(r5, 0x33, r5, b5);
  r6 = VMSUB(r6, 0x33, r6, b6); r7 = VMSUB(r7, 0x33, r7, b7);
  // r = D'+H'-p | D+H-p | C'-G' | C-G | B'+F'-p' | B+F-p | A'-E' | A-E
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);

  // get sign mask + carry propagation
  r1  = VMADD(r1 , 0x55, r1 , VSRA(r0, BRADIX)); 
  r4  = VMADD(r4 , 0x55, r4 , VSHUF<0x4E>(r0)); r0 = VAND(r0, bmask0);
  r2  = VMADD(r2 , 0x55, r2 , VSRA(r1, BRADIX)); 
  r5  = VMADD(r5 , 0x55, r5 , VSHUF<0x4E>(r1)); r1 = VAND(r1, bmask0);
  r3  = VMADD(r3 , 0x55, r3 , VSRA(r2, BRADIX)); 
  r6  = VMADD(r6 , 0x55, r6 , VSHUF<0x4E>(r2)); r2 = VAND(r2, bmask0);
  r4  = VMADD(r4 , 0x55, r4 , VSRA(r3, BRADIX)); 
  r7  = VMADD(r7 , 0x55, r7 , VSHUF<0x4E>(r3)); r3 = VAND(r3, bmask0);
  r5  = VMADD(r5 , 0x55, r5 , VSRA(r4, BRADIX)); 
  r8  = VMADD(r8 , 0x55, r8 , VSHUF<0x4E>(r4)); r4 = VAND(r4, bmask0);
  r6  = VMADD(r6 , 0x55, r6 , VSRA(r5, BRADIX)); 
  r9  = VMADD(r9 , 0x55, r9 , VSHUF<0x4E>(r5)); r5 = VAND(r5, bmask0);
  r7  = VMADD(r7 , 0x55, r7 , VSRA(r6, BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSHUF<0x4E>(r6)); r6 = VAND(r6, bmask0);
  r8  = VMADD(r8 , 0x55, r8 , VSRA(r7, BRADIX)); 
  r11 = VMADD(r11, 0x55, r11, VSHUF<0x4E>(r7)); r7 = VAND(r7, bmask0);
  t0  = VMADD(r9 , 0x55, r9 , VSRA(r8, BRADIX));
  t0  = VMADD(r10, 0x55, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0x55, r11, VSRA(t0, BRADIX));
  t0  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0  = VMADD(r9 , 0xAA, r9 , VSRA(t0, BRADIX));
  t0  = VMADD(r10, 0xAA, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0xAA, r11, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask)); 

  // carry propagation 
  r9  = VMADD(r9 , 0x55, r9 , VSRA(r8 , BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSRA(r9 , BRADIX));
  r11 = VMADD(r11, 0x55, r11, VSRA(r10, BRADIX));
  r8  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(r11, BRADIX)));
  r9  = VMADD(r9 , 0xAA, r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask1);
  r10 = VMADD(r10, 0xAA, r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask1);
  r11 = VMADD(r11, 0xAA, r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask1);
  r11 = VAND(r11, bmask1);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
}

void aass_fpx2_4x2w(fpx2_4x2w& r, const fpx2_4x2w& a, const fpx2_4x2w& b)
{
  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i b0  = b[0 ], b1  = b[1 ], b2  = b[2 ], b3  = b[3 ];
  cc512i b4  = b[4 ], b5  = b[5 ], b6  = b[6 ], b7  = b[7 ];
  cc512i b8  = b[8 ], b9  = b[9 ], b10 = b[10], b11 = b[11];
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask0 = VSET(0, BMASK, 0, BMASK, 0, BMASK, 0, BMASK);
  const cc512i bmask1 = VSET1(BMASK);
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11; 
  cc512i t0, t1, t2, t3, smask;

  // r = D'+H' | D+H | C'+G' | C+G' | B' | B | A' | A 
  r0  = VMADD(a0 , 0xF0, a0 , b0 ); r1  = VMADD(a1 , 0xF0, a1 , b1 );
  r2  = VMADD(a2 , 0xF0, a2 , b2 ); r3  = VMADD(a3 , 0xF0, a3 , b3 );
  r4  = VMADD(a4 , 0xF0, a4 , b4 ); r5  = VMADD(a5 , 0xF0, a5 , b5 );
  r6  = VMADD(a6 , 0xF0, a6 , b6 ); r7  = VMADD(a7 , 0xF0, a7 , b7 );
  r8  = VMADD(a8 , 0xF0, a8 , b8 ); r9  = VMADD(a9 , 0xF0, a9 , b9 );
  r10 = VMADD(a10, 0xF0, a10, b10); r11 = VMADD(a11, 0xF0, a11, b11);

  // t = p' | p | p' | p | F' | F | E' | E
  t0 = VMBLEND(0xF0, b8 , p0); t1 = VMBLEND(0xF0, b9 , p1);
  t2 = VMBLEND(0xF0, b10, p2); t3 = VMBLEND(0xF0, b11, p3);

  // r = D'+H' | D+H | C'+G' | C+G | B'-F' | B-F | A'-E' | A-E
  r0 = VMSUB(r0, 0x0F, r0, b0); r1 = VMSUB(r1, 0x0F, r1, b1);
  r2 = VMSUB(r2, 0x0F, r2, b2); r3 = VMSUB(r3, 0x0F, r3, b3);
  r4 = VMSUB(r4, 0x0F, r4, b4); r5 = VMSUB(r5, 0x0F, r5, b5);
  r6 = VMSUB(r6, 0x0F, r6, b6); r7 = VMSUB(r7, 0x0F, r7, b7);
  // r = D'+H'-p | D+H-p | C'+G'-p' | C+G-p | B'-F' | B-F | A'-E' | A-E
  r8  = VSUB(r8 , t0); r9  = VSUB(r9 , t1);
  r10 = VSUB(r10, t2); r11 = VSUB(r11, t3);

  // get sign mask + carry propagation
  r1  = VMADD(r1 , 0x55, r1 , VSRA(r0, BRADIX)); 
  r4  = VMADD(r4 , 0x55, r4 , VSHUF<0x4E>(r0)); r0 = VAND(r0, bmask0);
  r2  = VMADD(r2 , 0x55, r2 , VSRA(r1, BRADIX)); 
  r5  = VMADD(r5 , 0x55, r5 , VSHUF<0x4E>(r1)); r1 = VAND(r1, bmask0);
  r3  = VMADD(r3 , 0x55, r3 , VSRA(r2, BRADIX)); 
  r6  = VMADD(r6 , 0x55, r6 , VSHUF<0x4E>(r2)); r2 = VAND(r2, bmask0);
  r4  = VMADD(r4 , 0x55, r4 , VSRA(r3, BRADIX)); 
  r7  = VMADD(r7 , 0x55, r7 , VSHUF<0x4E>(r3)); r3 = VAND(r3, bmask0);
  r5  = VMADD(r5 , 0x55, r5 , VSRA(r4, BRADIX)); 
  r8  = VMADD(r8 , 0x55, r8 , VSHUF<0x4E>(r4)); r4 = VAND(r4, bmask0);
  r6  = VMADD(r6 , 0x55, r6 , VSRA(r5, BRADIX)); 
  r9  = VMADD(r9 , 0x55, r9 , VSHUF<0x4E>(r5)); r5 = VAND(r5, bmask0);
  r7  = VMADD(r7 , 0x55, r7 , VSRA(r6, BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSHUF<0x4E>(r6)); r6 = VAND(r6, bmask0);
  r8  = VMADD(r8 , 0x55, r8 , VSRA(r7, BRADIX)); 
  r11 = VMADD(r11, 0x55, r11, VSHUF<0x4E>(r7)); r7 = VAND(r7, bmask0);
  t0  = VMADD(r9 , 0x55, r9 , VSRA(r8, BRADIX));
  t0  = VMADD(r10, 0x55, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0x55, r11, VSRA(t0, BRADIX));
  t0  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(t0, BRADIX)));
  t0  = VMADD(r9 , 0xAA, r9 , VSRA(t0, BRADIX));
  t0  = VMADD(r10, 0xAA, r10, VSRA(t0, BRADIX));
  t0  = VMADD(r11, 0xAA, r11, VSRA(t0, BRADIX));

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(t0, 63);
  smask = VSHUF<0xEE>(smask);

  // r = r + (p & smask)
  r8  = VADD(r8 , VAND(p0, smask)); r9  = VADD(r9 , VAND(p1, smask));
  r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask)); 

  // carry propagation 
  r9  = VMADD(r9 , 0x55, r9 , VSRA(r8 , BRADIX)); 
  r10 = VMADD(r10, 0x55, r10, VSRA(r9 , BRADIX));
  r11 = VMADD(r11, 0x55, r11, VSRA(r10, BRADIX));
  r8  = VMADD(r8 , 0xAA, r8 , VZSHUF<0x4E>(0xCCCC, VSRA(r11, BRADIX)));
  r9  = VMADD(r9 , 0xAA, r9 , VSRA(r8 , BRADIX)); r8  = VAND(r8 , bmask1);
  r10 = VMADD(r10, 0xAA, r10, VSRA(r9 , BRADIX)); r9  = VAND(r9 , bmask1);
  r11 = VMADD(r11, 0xAA, r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask1);
  r11 = VAND(r11, bmask1);

  r[0 ] = r0 ; r[1 ] = r1 ; r[2 ] = r2 ; r[3 ] = r3 ;
  r[4 ] = r4 ; r[5 ] = r5 ; r[6 ] = r6 ; r[7 ] = r7 ;
  r[8 ] = r8 ; r[9 ] = r9 ; r[10] = r10; r[11] = r11;
}

void mul_fpx2_8x1w_v2(fpx2_8x1w& r, const fp_8x1w& a, const fp_8x1w& b)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);
  BLST_ASSERT_NO_OVERLAP_AVX(r, b);

  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
  const cc512i bmask = VSET1(BMASK);
  cc512i ta0, ta1, ta2, ta3;
  cc512i tb0, tb1, tb2, tb3;
  cc512i z0  = VZERO(), z1  = VZERO(), z2  = VZERO(), z3  = VZERO();
  cc512i z4  = VZERO(), z5  = VZERO(), z6  = VZERO(), z7  = VZERO();
  cc512i z8  = VZERO(), z9  = VZERO(), z10 = VZERO(), z11 = VZERO();
  cc512i z12 = VZERO(), z13 = VZERO(), z14 = VZERO(), z15 = VZERO();
  cc512i y0  = VZERO(), y1  = VZERO(), y2  = VZERO(), y3  = VZERO();
  cc512i y4  = VZERO(), y5  = VZERO(), y6  = VZERO(), y7  = VZERO();
  cc512i y8  = VZERO(), y9  = VZERO(), y10 = VZERO(), y11 = VZERO();
  cc512i y12 = VZERO(), y13 = VZERO(), y14 = VZERO(), y15 = VZERO();
  cc512i m0  = VZERO(), m1  = VZERO(), m2  = VZERO(), m3  = VZERO();
  cc512i m4  = VZERO(), m5  = VZERO(), m6  = VZERO(), m7  = VZERO();

  // compute zL(z0-z7) by aL(a0-a3) * bL(b0-b4)

  z0 = VMACLO(z0, a0, b0);
  y0 = VMACHI(y0, a0, b0);
  y0 = VSHL(y0, BALIGN);

  z1 = VMACLO(y0, a0, b1); z1 = VMACLO(z1, a1, b0);
  y1 = VMACHI(y1, a0, b1); y1 = VMACHI(y1, a1, b0);
  y1 = VSHL(y1, BALIGN);

  z2 = VMACLO(y1, a0, b2); z2 = VMACLO(z2, a1, b1); z2 = VMACLO(z2, a2, b0);
  y2 = VMACHI(y2, a0, b2); y2 = VMACHI(y2, a1, b1); y2 = VMACHI(y2, a2, b0);
  y2 = VSHL(y2, BALIGN);

  z3 = VMACLO(y2, a0, b3); z3 = VMACLO(z3, a1, b2); z3 = VMACLO(z3, a2, b1); 
  z3 = VMACLO(z3, a3, b0);
  y3 = VMACHI(y3, a0, b3); y3 = VMACHI(y3, a1, b2); y3 = VMACHI(y3, a2, b1); 
  y3 = VMACHI(y3, a3, b0);
  y3 = VSHL(y3, BALIGN);

  z4 = VMACLO(y3, a1, b3); z4 = VMACLO(z4, a2, b2); z4 = VMACLO(z4, a3, b1);
  y4 = VMACHI(y4, a1, b3); y4 = VMACHI(y4, a2, b2); y4 = VMACHI(y4, a3, b1);
  y4 = VSHL(y4, BALIGN);

  z5 = VMACLO(y4, a2, b3); z5 = VMACLO(z5, a3, b2);
  y5 = VMACHI(y5, a2, b3); y5 = VMACHI(y5, a3, b2);
  y5 = VSHL(y5, BALIGN);

  z6 = VMACLO(y5, a3, b3);
  y6 = VMACHI(y6, a3, b3);
  y6 = VSHL(y6, BALIGN);

  z7 = y6;

  // compute zH(z8-z15) by aH(a4-a7) * bH(b4-b7)

  z8 = VMACLO(z8, a4, b4);
  y8 = VMACHI(y8, a4, b4);
  y8 = VSHL(y8, BALIGN);

  z9 = VMACLO(y8, a4, b5); z9 = VMACLO(z9, a5, b4);
  y9 = VMACHI(y9, a4, b5); y9 = VMACHI(y9, a5, b4);
  y9 = VSHL(y9, BALIGN);

  z10 = VMACLO(y9,  a4, b6); z10 = VMACLO(z10, a5, b5); 
  z10 = VMACLO(z10, a6, b4);
  y10 = VMACHI(y10, a4, b6); y10 = VMACHI(y10, a5, b5); 
  y10 = VMACHI(y10, a6, b4);
  y10 = VSHL(y10, BALIGN);

  z11 = VMACLO(y10, a4, b7); z11 = VMACLO(z11, a5, b6);
  z11 = VMACLO(z11, a6, b5); z11 = VMACLO(z11, a7, b4);
  y11 = VMACHI(y11, a4, b7); y11 = VMACHI(y11, a5, b6);
  y11 = VMACHI(y11, a6, b5); y11 = VMACHI(y11, a7, b4);
  y11 = VSHL(y11, BALIGN);

  z12 = VMACLO(y11, a5, b7); z12 = VMACLO(z12, a6, b6);
  z12 = VMACLO(z12, a7, b5);
  y12 = VMACHI(y12, a5, b7); y12 = VMACHI(y12, a6, b6);
  y12 = VMACHI(y12, a7, b5);
  y12 = VSHL(y12, BALIGN);

  z13 = VMACLO(y12, a6, b7); z13 = VMACLO(z13, a7, b6);
  y13 = VMACHI(y13, a6, b7); y13 = VMACHI(y13, a7, b6);
  y13 = VSHL(y13, BALIGN);

  z14 = VMACLO(y13, a7, b7);
  y14 = VMACHI(y14, a7, b7);
  y14 = VSHL(y14, BALIGN);

  z15 = y14;

  // ta(ta0-ta3) = aL(a0-a3) + aH(a4-a7)
  ta0 = VADD(a0, a4); ta1 = VADD(a1, a5);
  ta2 = VADD(a2, a6); ta3 = VADD(a3, a7);

  // tb(tb0-tb3) = bL(b0-b3) + bH(b4-b7)
  tb0 = VADD(b0, b4); tb1 = VADD(b1, b5); 
  tb2 = VADD(b2, b6); tb3 = VADD(b3, b7);

  // zM = ta * tb - zL - zH 
  
  y0 = y1 = y2 = y3 = y4 = y5 = y6 = y7 = VZERO();

  m0 = VMACLO(m0, ta0, tb0);
  y0 = VMACHI(y0, ta0, tb0);
  y0 = VSHL(y0, BALIGN);

  m1 = VMACLO(y0, ta0, tb1); m1 = VMACLO(m1, ta1, tb0);
  y1 = VMACHI(y1, ta0, tb1); y1 = VMACHI(y1, ta1, tb0);
  y1 = VSHL(y1, BALIGN);

  m2 = VMACLO(y1, ta0, tb2); m2 = VMACLO(m2, ta1, tb1); 
  m2 = VMACLO(m2, ta2, tb0);
  y2 = VMACHI(y2, ta0, tb2); y2 = VMACHI(y2, ta1, tb1); 
  y2 = VMACHI(y2, ta2, tb0);
  y2 = VSHL(y2, BALIGN);

  m3 = VMACLO(y2, ta0, tb3); m3 = VMACLO(m3, ta1, tb2); 
  m3 = VMACLO(m3, ta2, tb1); m3 = VMACLO(m3, ta3, tb0);
  y3 = VMACHI(y3, ta0, tb3); y3 = VMACHI(y3, ta1, tb2); 
  y3 = VMACHI(y3, ta2, tb1); y3 = VMACHI(y3, ta3, tb0);
  y3 = VSHL(y3, BALIGN);

  m4 = VMACLO(y3, ta1, tb3); m4 = VMACLO(m4, ta2, tb2); 
  m4 = VMACLO(m4, ta3, tb1);
  y4 = VMACHI(y4, ta1, tb3); y4 = VMACHI(y4, ta2, tb2); 
  y4 = VMACHI(y4, ta3, tb1);
  y4 = VSHL(y4, BALIGN);

  m5 = VMACLO(y4, ta2, tb3); m5 = VMACLO(m5, ta3, tb2);
  y5 = VMACHI(y5, ta2, tb3); y5 = VMACHI(y5, ta3, tb2);
  y5 = VSHL(y5, BALIGN);

  m6 = VMACLO(y5, ta3, tb3);
  y6 = VMACHI(y6, ta3, tb3);
  y6 = VSHL(y6, BALIGN);

  m7 = y6;

  m0 = VSUB(m0, VADD(z0, z8 )); m1 = VSUB(m1, VADD(z1, z9 ));
  m2 = VSUB(m2, VADD(z2, z10)); m3 = VSUB(m3, VADD(z3, z11));
  m4 = VSUB(m4, VADD(z4, z12)); m5 = VSUB(m5, VADD(z5, z13));
  m6 = VSUB(m6, VADD(z6, z14)); m7 = VSUB(m7, VADD(z7, z15));

  // z = z + zM
  z4  = VADD(z4 , m0); z5  = VADD(z5 , m1);
  z6  = VADD(z6 , m2); z7  = VADD(z7 , m3);
  z8  = VADD(z8 , m4); z9  = VADD(z9 , m5);
  z10 = VADD(z10, m6); z11 = VADD(z11, m7);

  r[0 ] = z0 ; r[1 ] = z1 ; r[2 ] = z2 ; r[3 ] = z3 ;
  r[4 ] = z4 ; r[5 ] = z5 ; r[6 ] = z6 ; r[7 ] = z7 ;
  r[8 ] = z8 ; r[9 ] = z9 ; r[10] = z10; r[11] = z11;
  r[12] = z12; r[13] = z13; r[14] = z14; r[15] = z15; 
}

void mul_fpx2_8x1w(fpx2_8x1w& r, const fp_8x1w& a, const fp_8x1w& b)
{
  mul_fpx2_8x1w_v2(r, a, b);
}

void mul_fpx2_4x2w_v1(fpx2_4x2w& r, const fp_4x2w& a, const fp_4x2w& b)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);
  BLST_ASSERT_NO_OVERLAP_AVX(r, b);

  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  cc512i z0 = VZERO(), z1 = VZERO(), z2  = VZERO(), z3  = VZERO();
  cc512i z4 = VZERO(), z5 = VZERO(), z6  = VZERO(), z7  = VZERO();
  cc512i z8 = VZERO(), z9 = VZERO(), z10 = VZERO(), z11 = VZERO();
  cc512i y0 = VZERO(), y1 = VZERO(), y2  = VZERO(), y3  = VZERO();
  cc512i y4 = VZERO(), y5 = VZERO(), y6  = VZERO(), y7  = VZERO();
  cc512i y8 = VZERO(), y9 = VZERO(), y10 = VZERO(), tb;

  tb = VSHUF<0x44>(b0);
  z0 = VMACLO(z0, tb, a0); z1 = VMACLO(z1, tb, a1);
  z2 = VMACLO(z2, tb, a2); z3 = VMACLO(z3, tb, a3);
  y0 = VMACHI(y0, tb, a0); y1 = VMACHI(y1, tb, a1);
  y2 = VMACHI(y2, tb, a2); y3 = VMACHI(y3, tb, a3);

  tb = VSHUF<0x44>(b1);
  z1 = VMACLO(z1, tb, a0); z2 = VMACLO(z2, tb, a1);
  z3 = VMACLO(z3, tb, a2); z4 = VMACLO(z4, tb, a3);
  y1 = VMACHI(y1, tb, a0); y2 = VMACHI(y2, tb, a1);
  y3 = VMACHI(y3, tb, a2); y4 = VMACHI(y4, tb, a3);

  tb = VSHUF<0x44>(b2);
  z2 = VMACLO(z2, tb, a0); z3 = VMACLO(z3, tb, a1);
  z4 = VMACLO(z4, tb, a2); z5 = VMACLO(z5, tb, a3);
  y2 = VMACHI(y2, tb, a0); y3 = VMACHI(y3, tb, a1);
  y4 = VMACHI(y4, tb, a2); y5 = VMACHI(y5, tb, a3);

  tb = VSHUF<0x44>(b3);
  z3 = VMACLO(z3, tb, a0); z4 = VMACLO(z4, tb, a1);
  z5 = VMACLO(z5, tb, a2); z6 = VMACLO(z6, tb, a3);
  y3 = VMACHI(y3, tb, a0); y4 = VMACHI(y4, tb, a1);
  y5 = VMACHI(y5, tb, a2); y6 = VMACHI(y6, tb, a3);

  tb = VSHUF<0xEE>(b0);
  z4 = VMACLO(z4, tb, a0); z5 = VMACLO(z5, tb, a1);
  z6 = VMACLO(z6, tb, a2); z7 = VMACLO(z7, tb, a3);
  y4 = VMACHI(y4, tb, a0); y5 = VMACHI(y5, tb, a1);
  y6 = VMACHI(y6, tb, a2); y7 = VMACHI(y7, tb, a3);

  tb = VSHUF<0xEE>(b1);
  z5 = VMACLO(z5, tb, a0); z6 = VMACLO(z6, tb, a1);
  z7 = VMACLO(z7, tb, a2); z8 = VMACLO(z8, tb, a3);
  y5 = VMACHI(y5, tb, a0); y6 = VMACHI(y6, tb, a1);
  y7 = VMACHI(y7, tb, a2); y8 = VMACHI(y8, tb, a3);

  tb = VSHUF<0xEE>(b2);
  z6 = VMACLO(z6, tb, a0); z7 = VMACLO(z7, tb, a1);
  z8 = VMACLO(z8, tb, a2); z9 = VMACLO(z9, tb, a3);
  y6 = VMACHI(y6, tb, a0); y7 = VMACHI(y7, tb, a1);
  y8 = VMACHI(y8, tb, a2); y9 = VMACHI(y9, tb, a3);

  tb = VSHUF<0xEE>(b3);
  z7 = VMACLO(z7, tb, a0); z8 = VMACLO(z8, tb, a1);
  z9 = VMACLO(z9, tb, a2); z10 = VMACLO(z10, tb, a3);
  y7 = VMACHI(y7, tb, a0); y8 = VMACHI(y8, tb, a1);
  y9 = VMACHI(y9, tb, a2); y10 = VMACHI(y10, tb, a3);

  z1  = VADD(z1 , VSHL(y0 , BALIGN));
  z2  = VADD(z2 , VSHL(y1 , BALIGN));
  z3  = VADD(z3 , VSHL(y2 , BALIGN));
  z4  = VADD(z4 , VSHL(y3 , BALIGN));
  z5  = VADD(z5 , VSHL(y4 , BALIGN));
  z6  = VADD(z6 , VSHL(y5 , BALIGN));
  z7  = VADD(z7 , VSHL(y6 , BALIGN));
  z8  = VADD(z8 , VSHL(y7 , BALIGN));
  z9  = VADD(z9 , VSHL(y8 , BALIGN));
  z10 = VADD(z10, VSHL(y9 , BALIGN));
  z11 = VADD(z11, VSHL(y10, BALIGN));

  r[0 ] = z0 ; r[1 ] = z1 ; r[2 ] = z2 ; r[3 ] = z3 ;
  r[4 ] = z4 ; r[5 ] = z5 ; r[6 ] = z6 ; r[7 ] = z7 ;
  r[8 ] = z8 ; r[9 ] = z9 ; r[10] = z10; r[11] = z11;
}

void mul_fpx2_4x2w(fpx2_4x2w& r, const fp_4x2w& a, const fp_4x2w& b)
{
  mul_fpx2_4x2w_v1(r, a, b);
}

void redc_fpx2_8x1w(fp_8x1w& r, const fpx2_8x1w& a)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);

  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  cc512i y0  = VZERO(), y1  = VZERO(), y2  = VZERO(), y3  = VZERO();
  cc512i y4  = VZERO(), y5  = VZERO(), y6  = VZERO(), y7  = VZERO();
  cc512i y8  = VZERO(), y9  = VZERO(), y10 = VZERO(), y11 = VZERO();
  cc512i y12 = VZERO(), y13 = VZERO(), y14 = VZERO();
  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK), montw = VSET1(MONT_W_R48), zero = VZERO();
  cc512i r0, r1, r2, r3, r4, r5, r6, r7, smask, u, t0;

  u  = VAND(VMACLO(zero, a0, montw), bmask);
  a0 = VMACLO(a0, u, p0); a1 = VMACLO(a1, u, p1); a2 = VMACLO(a2, u, p2);
  a3 = VMACLO(a3, u, p3); a4 = VMACLO(a4, u, p4); a5 = VMACLO(a5, u, p5);
  a6 = VMACLO(a6, u, p6); a7 = VMACLO(a7, u, p7); 
  y0 = VMACHI(y0, u, p0); y1 = VMACHI(y1, u, p1); y2 = VMACHI(y2, u, p2);
  y3 = VMACHI(y3, u, p3); y4 = VMACHI(y4, u, p4); y5 = VMACHI(y5, u, p5);
  y6 = VMACHI(y6, u, p6); y7 = VMACHI(y7, u, p7); 
  a1 = VADD(VADD(a1, VSRA(a0, BRADIX)), VSHL(y0, BALIGN));

  u  = VAND(VMACLO(zero, a1, montw), bmask);
  a1 = VMACLO(a1, u, p0); a2 = VMACLO(a2, u, p1); a3 = VMACLO(a3, u, p2);
  a4 = VMACLO(a4, u, p3); a5 = VMACLO(a5, u, p4); a6 = VMACLO(a6, u, p5);
  a7 = VMACLO(a7, u, p6); a8 = VMACLO(a8, u, p7); 
  y1 = VMACHI(y1, u, p0); y2 = VMACHI(y2, u, p1); y3 = VMACHI(y3, u, p2);
  y4 = VMACHI(y4, u, p3); y5 = VMACHI(y5, u, p4); y6 = VMACHI(y6, u, p5);
  y7 = VMACHI(y7, u, p6); y8 = VMACHI(y8, u, p7); 
  a2 = VADD(VADD(a2, VSRA(a1, BRADIX)), VSHL(y1, BALIGN));

  u  = VAND(VMACLO(zero, a2, montw), bmask);
  a2 = VMACLO(a2, u, p0); a3 = VMACLO(a3, u, p1); a4 = VMACLO(a4, u, p2);
  a5 = VMACLO(a5, u, p3); a6 = VMACLO(a6, u, p4); a7 = VMACLO(a7, u, p5);
  a8 = VMACLO(a8, u, p6); a9 = VMACLO(a9, u, p7); 
  y2 = VMACHI(y2, u, p0); y3 = VMACHI(y3, u, p1); y4 = VMACHI(y4, u, p2);
  y5 = VMACHI(y5, u, p3); y6 = VMACHI(y6, u, p4); y7 = VMACHI(y7, u, p5);
  y8 = VMACHI(y8, u, p6); y9 = VMACHI(y9, u, p7); 
  a3 = VADD(VADD(a3, VSRA(a2, BRADIX)), VSHL(y2, BALIGN));

  u  = VAND(VMACLO(zero, a3, montw), bmask);
  a3 = VMACLO(a3, u, p0); a4  = VMACLO(a4,  u, p1); a5 = VMACLO(a5, u, p2);
  a6 = VMACLO(a6, u, p3); a7  = VMACLO(a7,  u, p4); a8 = VMACLO(a8, u, p5);
  a9 = VMACLO(a9, u, p6); a10 = VMACLO(a10, u, p7); 
  y3 = VMACHI(y3, u, p0); y4  = VMACHI(y4 , u, p1); y5 = VMACHI(y5, u, p2);
  y6 = VMACHI(y6, u, p3); y7  = VMACHI(y7 , u, p4); y8 = VMACHI(y8, u, p5);
  y9 = VMACHI(y9, u, p6); y10 = VMACHI(y10, u, p7); 
  a4 = VADD(VADD(a4, VSRA(a3, BRADIX)), VSHL(y3, BALIGN));

  u   = VAND(VMACLO(zero, a4, montw), bmask);
  a4  = VMACLO(a4 , u, p0); a5  = VMACLO(a5 , u, p1); a6 = VMACLO(a6, u, p2);
  a7  = VMACLO(a7 , u, p3); a8  = VMACLO(a8 , u, p4); a9 = VMACLO(a9, u, p5);
  a10 = VMACLO(a10, u, p6); a11 = VMACLO(a11, u, p7); 
  y4  = VMACHI(y4 , u, p0); y5  = VMACHI(y5 , u, p1); y6 = VMACHI(y6, u, p2);
  y7  = VMACHI(y7 , u, p3); y8  = VMACHI(y8 , u, p4); y9 = VMACHI(y9, u, p5);
  y10 = VMACHI(y10, u, p6); y11 = VMACHI(y11, u, p7); 
  a5  = VADD(VADD(a5, VSRA(a4, BRADIX)), VSHL(y4, BALIGN));

  u   = VAND(VMACLO(zero, a5, montw), bmask);
  a5  = VMACLO(a5 , u, p0); a6  = VMACLO(a6 , u, p1); a7  = VMACLO(a7 , u, p2);
  a8  = VMACLO(a8 , u, p3); a9  = VMACLO(a9 , u, p4); a10 = VMACLO(a10, u, p5);
  a11 = VMACLO(a11, u, p6); a12 = VMACLO(a12, u, p7); 
  y5  = VMACHI(y5 , u, p0); y6  = VMACHI(y6 , u, p1); y7  = VMACHI(y7 , u, p2);
  y8  = VMACHI(y8 , u, p3); y9  = VMACHI(y9 , u, p4); y10 = VMACHI(y10, u, p5);
  y11 = VMACHI(y11, u, p6); y12 = VMACHI(y12, u, p7); 
  a6  = VADD(VADD(a6, VSRA(a5, BRADIX)), VSHL(y5, BALIGN));  

  u   = VAND(VMACLO(zero, a6, montw), bmask);
  a6  = VMACLO(a6 , u, p0); a7  = VMACLO(a7 , u, p1); a8  = VMACLO(a8 , u, p2);
  a9  = VMACLO(a9 , u, p3); a10 = VMACLO(a10, u, p4); a11 = VMACLO(a11, u, p5);
  a12 = VMACLO(a12, u, p6); a13 = VMACLO(a13, u, p7); 
  y6  = VMACHI(y6 , u, p0); y7  = VMACHI(y7 , u, p1); y8  = VMACHI(y8 , u, p2);
  y9  = VMACHI(y9 , u, p3); y10 = VMACHI(y10, u, p4); y11 = VMACHI(y11, u, p5);
  y12 = VMACHI(y12, u, p6); y13 = VMACHI(y13, u, p7); 
  a7  = VADD(VADD(a7, VSRA(a6, BRADIX)), VSHL(y6, BALIGN));  

  u   = VAND(VMACLO(zero, a7, montw), bmask);
  a7  = VMACLO(a7 , u, p0); a8  = VMACLO(a8 , u, p1); a9  = VMACLO(a9 , u, p2);
  a10 = VMACLO(a10, u, p3); a11 = VMACLO(a11, u, p4); a12 = VMACLO(a12, u, p5);
  a13 = VMACLO(a13, u, p6); a14 = VMACLO(a14, u, p7); 
  y7  = VMACHI(y7 , u, p0); y8  = VMACHI(y8 , u, p1); y9  = VMACHI(y9 , u, p2);
  y10 = VMACHI(y10, u, p3); y11 = VMACHI(y11, u, p4); y12 = VMACHI(y12, u, p5);
  y13 = VMACHI(y13, u, p6); y14 = VMACHI(y14, u, p7); 
  a8  = VADD(VADD(a8, VSRA(a7, BRADIX)), VSHL(y7, BALIGN)); 

  a9  = VADD(a9 , VSHL(y8 , BALIGN));
  a10 = VADD(a10, VSHL(y9 , BALIGN));
  a11 = VADD(a11, VSHL(y10, BALIGN));
  a12 = VADD(a12, VSHL(y11, BALIGN));
  a13 = VADD(a13, VSHL(y12, BALIGN));
  a14 = VADD(a14, VSHL(y13, BALIGN));
  a15 = VADD(a15, VSHL(y14, BALIGN));

  // final subtraction 
  r0 = a8 ; r1 = a9 ; r2 = a10; r3 = a11;
  r4 = a12; r5 = a13; r6 = a14; r7 = a15;

  // r = r - p
  r0 = VSUB(r0, p0); r1 = VSUB(r1, p1); r2 = VSUB(r2, p2); r3 = VSUB(r3, p3);
  r4 = VSUB(r4, p4); r5 = VSUB(r5, p5); r6 = VSUB(r6, p6); r7 = VSUB(r7, p7);

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

void redc_fpx2_4x2w(fp_4x2w& r, const fpx2_4x2w& a)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);

  cc512i a0  = a[0 ], a1  = a[1 ], a2  = a[2 ], a3  = a[3 ];
  cc512i a4  = a[4 ], a5  = a[5 ], a6  = a[6 ], a7  = a[7 ];
  cc512i a8  = a[8 ], a9  = a[9 ], a10 = a[10], a11 = a[11];
  cc512i y0  = VZERO(), y1  = VZERO(), y2  = VZERO(), y3  = VZERO();
  cc512i y4  = VZERO(), y5  = VZERO(), y6  = VZERO(), y7  = VZERO();
  cc512i y8  = VZERO(), y9  = VZERO(), y10 = VZERO();
  const cc512i p0 = VSET(P48[4], P48[0], P48[4], P48[0], P48[4], P48[0], P48[4], P48[0]);
  const cc512i p1 = VSET(P48[5], P48[1], P48[5], P48[1], P48[5], P48[1], P48[5], P48[1]);
  const cc512i p2 = VSET(P48[6], P48[2], P48[6], P48[2], P48[6], P48[2], P48[6], P48[2]);
  const cc512i p3 = VSET(P48[7], P48[3], P48[7], P48[3], P48[7], P48[3], P48[7], P48[3]);
  const cc512i bmask = VSET1(BMASK), montw = VSET1(MONT_W_R48), zero = VZERO();
  cc512i r0, r1, r2, r3, smask, u, t0;

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a0, montw), bmask));
  a0 = VMACLO(a0, u, p0); a1 = VMACLO(a1, u, p1);
  a2 = VMACLO(a2, u, p2); a3 = VMACLO(a3, u, p3);
  y0 = VMACHI(y0, u, p0); y1 = VMACHI(y1, u, p1);
  y2 = VMACHI(y2, u, p2); y3 = VMACHI(y3, u, p3);
  a4 = VMADD(a4, 0x55, a4, VSHUF<0x4E>(a0));
  a1 = VMADD(a1, 0x55, a1, VSRA (a0, BRADIX));
  a1 = VADD(a1, VSHL(y0, BALIGN));

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a1, montw), bmask));
  a1 = VMACLO(a1, u, p0); a2 = VMACLO(a2, u, p1);
  a3 = VMACLO(a3, u, p2); a4 = VMACLO(a4, u, p3);
  y1 = VMACHI(y1, u, p0); y2 = VMACHI(y2, u, p1);
  y3 = VMACHI(y3, u, p2); y4 = VMACHI(y4, u, p3);
  a5 = VMADD(a5, 0x55, a5, VSHUF<0x4E>(a1));
  a2 = VMADD(a2, 0x55, a2, VSRA (a1, BRADIX));
  a2 = VADD(a2, VSHL(y1, BALIGN));

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a2, montw), bmask));
  a2 = VMACLO(a2, u, p0); a3 = VMACLO(a3, u, p1);
  a4 = VMACLO(a4, u, p2); a5 = VMACLO(a5, u, p3);
  y2 = VMACHI(y2, u, p0); y3 = VMACHI(y3, u, p1);
  y4 = VMACHI(y4, u, p2); y5 = VMACHI(y5, u, p3);
  a6 = VMADD(a6, 0x55, a6, VSHUF<0x4E>(a2));
  a3 = VMADD(a3, 0x55, a3, VSRA (a2, BRADIX));
  a3 = VADD(a3, VSHL(y2, BALIGN));

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a3, montw), bmask));
  a3 = VMACLO(a3, u, p0); a4 = VMACLO(a4, u, p1);
  a5 = VMACLO(a5, u, p2); a6 = VMACLO(a6, u, p3);
  y3 = VMACHI(y3, u, p0); y4 = VMACHI(y4, u, p1);
  y5 = VMACHI(y5, u, p2); y6 = VMACHI(y6, u, p3);
  a7 = VMADD(a7, 0x55, a7, VSHUF<0x4E>(a3));
  a4 = VMADD(a4, 0x55, a4, VSRA (a3, BRADIX));
  a4 = VADD(a4, VSHL(y3, BALIGN));

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a4, montw), bmask));
  a4 = VMACLO(a4, u, p0); a5 = VMACLO(a5, u, p1);
  a6 = VMACLO(a6, u, p2); a7 = VMACLO(a7, u, p3);
  y4 = VMACHI(y4, u, p0); y5 = VMACHI(y5, u, p1);
  y6 = VMACHI(y6, u, p2); y7 = VMACHI(y7, u, p3);
  a8 = VMADD(a8, 0x55, a8, VSHUF<0x4E>(a4));
  a5 = VMADD(a5, 0x55, a5, VSRA (a4, BRADIX));
  a5 = VADD(a5, VSHL(y4, BALIGN));

  u  = VSHUF<0x44>(VAND(VMACLO(zero, a5, montw), bmask));
  a5 = VMACLO(a5, u, p0); a6 = VMACLO(a6, u, p1);
  a7 = VMACLO(a7, u, p2); a8 = VMACLO(a8, u, p3);
  y5 = VMACHI(y5, u, p0); y6 = VMACHI(y6, u, p1);
  y7 = VMACHI(y7, u, p2); y8 = VMACHI(y8, u, p3);
  a9 = VMADD(a9, 0x55, a9, VSHUF<0x4E>(a5));
  a6 = VMADD(a6, 0x55, a6, VSRA (a5, BRADIX));
  a6 = VADD(a6, VSHL(y5, BALIGN));

  u   = VSHUF<0x44>(VAND(VMACLO(zero, a6, montw), bmask));
  a6  = VMACLO(a6, u, p0); a7  = VMACLO(a7, u, p1);
  a8  = VMACLO(a8, u, p2); a9  = VMACLO(a9, u, p3);
  y6  = VMACHI(y6, u, p0); y7  = VMACHI(y7, u, p1);
  y8  = VMACHI(y8, u, p2); y9  = VMACHI(y9, u, p3);
  a10 = VMADD(a10, 0x55, a10, VSHUF<0x4E>(a6));
  a7  = VMADD(a7 , 0x55, a7 , VSRA (a6, BRADIX));
  a7  = VADD(a7, VSHL(y6, BALIGN));

  u   = VSHUF<0x44>(VAND(VMACLO(zero, a7, montw), bmask));
  a7  = VMACLO(a7, u, p0); a8   = VMACLO(a8 , u, p1);
  a9  = VMACLO(a9, u, p2); a10  = VMACLO(a10, u, p3);
  y7  = VMACHI(y7, u, p0); y8   = VMACHI(y8 , u, p1);
  y9  = VMACHI(y9, u, p2); y10  = VMACHI(y10, u, p3);
  a11 = VMADD(a11, 0x55, a11, VSHUF<0x4E>(a7));
  a8  = VMADD(a8 , 0x55, a8 , VSRA (a7, BRADIX));
  a8  = VADD(a8, VSHL(y7, BALIGN));

  a9  = VADD(a9 , VSHL(y8 , BALIGN));
  a10 = VADD(a10, VSHL(y9 , BALIGN));
  a11 = VADD(a11, VSHL(y10, BALIGN));

  // final subtraction
  r0 = a8; r1 = a9; r2 = a10; r3 = a11;

  // r = r - p
  r0 = VSUB(r0, p0); r1 = VSUB(r1, p1); r2 = VSUB(r2, p2); r3 = VSUB(r3, p3);

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

} // namespace blst

#endif
