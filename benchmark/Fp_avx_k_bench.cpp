#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#if defined(BLST_AVX512)

#include <array>
#include <cstddef>
#include <random>
#include <vector>

#include "blst/fields/Fp.h"
#include "blst/fields/Fp_avx.h"

using namespace blst;

namespace
{
void random_redc_fp(Fp& v, std::mt19937_64& gen)
{
  for (u64 i = 0; i < v.num_limb(); ++i) {
    v[i] = gen();
  }
  while (Fp::is_redc(v) == 0) {
    Fp::rshift(v, v, 1);
  }
}

struct Fp8x1KBenchData
{
  explicit Fp8x1KBenchData(const std::size_t n_packs)
      : n_packs_(n_packs)
      , a(n_packs_)
      , b(n_packs_)
      , wide(n_packs_)
      , lazy(n_packs_)
      , out(n_packs_)
  {
    std::mt19937_64 gen(0xA6B1F81EULL + static_cast<u64>(n_packs_));

    for (std::size_t i = 0; i < n_packs_; ++i) {
      std::array<Fp, 8> lane_a{};
      std::array<Fp, 8> lane_b{};
      for (std::size_t lane = 0; lane < 8; ++lane) {
        random_redc_fp(lane_a[lane], gen);
        random_redc_fp(lane_b[lane], gen);
      }
      FpPack8x1::from_aos(a[i], lane_a.data());
      FpPack8x1::from_aos(b[i], lane_b.data());
      FpPack8x1::mul_v2(wide[i], a[i], b[i]);
      UFpPack8x1::redc_core(lazy[i], wide[i]);
    }
  }

  const std::size_t n_packs_;
  std::vector<FpPack8x1> a;
  std::vector<FpPack8x1> b;
  std::vector<UFpPack8x1> wide;
  std::vector<LFpPack8x1> lazy;
  std::vector<FpPack8x1> out;
};

template <int Chains>
struct Fp8x1MulRedcStreamData
{
  static_assert(Chains >= 1 && Chains <= 8, "Chains must be in [1,8].");

  Fp8x1MulRedcStreamData()
  {
    std::mt19937_64 gen(0xD4E5F617ULL + static_cast<u64>(Chains));
    for (int c = 0; c < Chains; ++c) {
      std::array<Fp, 8> lane_a{};
      std::array<Fp, 8> lane_b{};
      for (std::size_t lane = 0; lane < 8; ++lane) {
        random_redc_fp(lane_a[lane], gen);
        random_redc_fp(lane_b[lane], gen);
      }
      FpPack8x1::from_aos(cur[c], lane_a.data());
      FpPack8x1::from_aos(mul_by[c], lane_b.data());
    }
  }

  std::array<FpPack8x1, Chains> cur{};
  std::array<FpPack8x1, Chains> mul_by{};
  std::array<UFpPack8x1, Chains> wide{};
  std::array<FpPack8x1, Chains> out{};
};

template <int Chains>
struct Fp8x1RedcCoreStreamData
{
  static_assert(Chains >= 1 && Chains <= 8, "Chains must be in [1,8].");

  Fp8x1RedcCoreStreamData()
  {
    std::mt19937_64 gen(0xE71B93D4ULL + static_cast<u64>(Chains));
    for (int c = 0; c < Chains; ++c) {
      std::array<Fp, 8> lane_a{};
      std::array<Fp, 8> lane_b{};
      for (std::size_t lane = 0; lane < 8; ++lane) {
        random_redc_fp(lane_a[lane], gen);
        random_redc_fp(lane_b[lane], gen);
      }
      FpPack8x1 ap{}, bp{};
      FpPack8x1::from_aos(ap, lane_a.data());
      FpPack8x1::from_aos(bp, lane_b.data());
      FpPack8x1::mul_v2(wide[c], ap, bp);
    }
  }

  std::array<UFpPack8x1, Chains> wide{};
  std::array<LFpPack8x1, Chains> lazy{};
};

BLST_FORCEINLINE void redc_core_hand_interleaved_k2_once(
    LFpPack8x1& out0,
    LFpPack8x1& out1,
    const UFpPack8x1& in0,
    const UFpPack8x1& in1)
{
  cc512i a0 = in0[0], a1 = in0[1], a2 = in0[2], a3 = in0[3];
  cc512i a4 = in0[4], a5 = in0[5], a6 = in0[6], a7 = in0[7];
  cc512i a8 = in0[8], a9 = in0[9], a10 = in0[10], a11 = in0[11];
  cc512i a12 = in0[12], a13 = in0[13], a14 = in0[14], a15 = in0[15];

  cc512i b0 = in1[0], b1 = in1[1], b2 = in1[2], b3 = in1[3];
  cc512i b4 = in1[4], b5 = in1[5], b6 = in1[6], b7 = in1[7];
  cc512i b8 = in1[8], b9 = in1[9], b10 = in1[10], b11 = in1[11];
  cc512i b12 = in1[12], b13 = in1[13], b14 = in1[14], b15 = in1[15];

  cc512i y0 = VZERO(), y1 = VZERO(), y2 = VZERO(), y3 = VZERO();
  cc512i y4 = VZERO(), y5 = VZERO(), y6 = VZERO(), y7 = VZERO();
  cc512i y8 = VZERO(), y9 = VZERO(), y10 = VZERO(), y11 = VZERO();
  cc512i y12 = VZERO(), y13 = VZERO(), y14 = VZERO();

  cc512i z0 = VZERO(), z1 = VZERO(), z2 = VZERO(), z3 = VZERO();
  cc512i z4 = VZERO(), z5 = VZERO(), z6 = VZERO(), z7 = VZERO();
  cc512i z8 = VZERO(), z9 = VZERO(), z10 = VZERO(), z11 = VZERO();
  cc512i z12 = VZERO(), z13 = VZERO(), z14 = VZERO();

  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK), montw = VSET1(MONT_W_R48), zero = VZERO();

  cc512i u0, u1;

  u0 = VAND(VMACLO(zero, a0, montw), bmask);
  u1 = VAND(VMACLO(zero, b0, montw), bmask);
  a0 = VMACLO(a0, u0, p0); b0 = VMACLO(b0, u1, p0);
  a1 = VMACLO(a1, u0, p1); b1 = VMACLO(b1, u1, p1);
  a2 = VMACLO(a2, u0, p2); b2 = VMACLO(b2, u1, p2);
  a3 = VMACLO(a3, u0, p3); b3 = VMACLO(b3, u1, p3);
  a4 = VMACLO(a4, u0, p4); b4 = VMACLO(b4, u1, p4);
  a5 = VMACLO(a5, u0, p5); b5 = VMACLO(b5, u1, p5);
  a6 = VMACLO(a6, u0, p6); b6 = VMACLO(b6, u1, p6);
  a7 = VMACLO(a7, u0, p7); b7 = VMACLO(b7, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a1 = VADD(VADD(a1, VSRA(a0, BRADIX)), VSHL(y0, BALIGN));
  b1 = VADD(VADD(b1, VSRA(b0, BRADIX)), VSHL(z0, BALIGN));

  u0 = VAND(VMACLO(zero, a1, montw), bmask);
  u1 = VAND(VMACLO(zero, b1, montw), bmask);
  a1 = VMACLO(a1, u0, p0); b1 = VMACLO(b1, u1, p0);
  a2 = VMACLO(a2, u0, p1); b2 = VMACLO(b2, u1, p1);
  a3 = VMACLO(a3, u0, p2); b3 = VMACLO(b3, u1, p2);
  a4 = VMACLO(a4, u0, p3); b4 = VMACLO(b4, u1, p3);
  a5 = VMACLO(a5, u0, p4); b5 = VMACLO(b5, u1, p4);
  a6 = VMACLO(a6, u0, p5); b6 = VMACLO(b6, u1, p5);
  a7 = VMACLO(a7, u0, p6); b7 = VMACLO(b7, u1, p6);
  a8 = VMACLO(a8, u0, p7); b8 = VMACLO(b8, u1, p7);
  y1 = VMACHI(y1, u0, p0); z1 = VMACHI(z1, u1, p0);
  y2 = VMACHI(y2, u0, p1); z2 = VMACHI(z2, u1, p1);
  y3 = VMACHI(y3, u0, p2); z3 = VMACHI(z3, u1, p2);
  y4 = VMACHI(y4, u0, p3); z4 = VMACHI(z4, u1, p3);
  y5 = VMACHI(y5, u0, p4); z5 = VMACHI(z5, u1, p4);
  y6 = VMACHI(y6, u0, p5); z6 = VMACHI(z6, u1, p5);
  y7 = VMACHI(y7, u0, p6); z7 = VMACHI(z7, u1, p6);
  y8 = VMACHI(y8, u0, p7); z8 = VMACHI(z8, u1, p7);
  a2 = VADD(VADD(a2, VSRA(a1, BRADIX)), VSHL(y1, BALIGN));
  b2 = VADD(VADD(b2, VSRA(b1, BRADIX)), VSHL(z1, BALIGN));

  u0 = VAND(VMACLO(zero, a2, montw), bmask);
  u1 = VAND(VMACLO(zero, b2, montw), bmask);
  a2 = VMACLO(a2, u0, p0); b2 = VMACLO(b2, u1, p0);
  a3 = VMACLO(a3, u0, p1); b3 = VMACLO(b3, u1, p1);
  a4 = VMACLO(a4, u0, p2); b4 = VMACLO(b4, u1, p2);
  a5 = VMACLO(a5, u0, p3); b5 = VMACLO(b5, u1, p3);
  a6 = VMACLO(a6, u0, p4); b6 = VMACLO(b6, u1, p4);
  a7 = VMACLO(a7, u0, p5); b7 = VMACLO(b7, u1, p5);
  a8 = VMACLO(a8, u0, p6); b8 = VMACLO(b8, u1, p6);
  a9 = VMACLO(a9, u0, p7); b9 = VMACLO(b9, u1, p7);
  y2 = VMACHI(y2, u0, p0); z2 = VMACHI(z2, u1, p0);
  y3 = VMACHI(y3, u0, p1); z3 = VMACHI(z3, u1, p1);
  y4 = VMACHI(y4, u0, p2); z4 = VMACHI(z4, u1, p2);
  y5 = VMACHI(y5, u0, p3); z5 = VMACHI(z5, u1, p3);
  y6 = VMACHI(y6, u0, p4); z6 = VMACHI(z6, u1, p4);
  y7 = VMACHI(y7, u0, p5); z7 = VMACHI(z7, u1, p5);
  y8 = VMACHI(y8, u0, p6); z8 = VMACHI(z8, u1, p6);
  y9 = VMACHI(y9, u0, p7); z9 = VMACHI(z9, u1, p7);
  a3 = VADD(VADD(a3, VSRA(a2, BRADIX)), VSHL(y2, BALIGN));
  b3 = VADD(VADD(b3, VSRA(b2, BRADIX)), VSHL(z2, BALIGN));

  u0 = VAND(VMACLO(zero, a3, montw), bmask);
  u1 = VAND(VMACLO(zero, b3, montw), bmask);
  a3 = VMACLO(a3, u0, p0); b3 = VMACLO(b3, u1, p0);
  a4 = VMACLO(a4, u0, p1); b4 = VMACLO(b4, u1, p1);
  a5 = VMACLO(a5, u0, p2); b5 = VMACLO(b5, u1, p2);
  a6 = VMACLO(a6, u0, p3); b6 = VMACLO(b6, u1, p3);
  a7 = VMACLO(a7, u0, p4); b7 = VMACLO(b7, u1, p4);
  a8 = VMACLO(a8, u0, p5); b8 = VMACLO(b8, u1, p5);
  a9 = VMACLO(a9, u0, p6); b9 = VMACLO(b9, u1, p6);
  a10 = VMACLO(a10, u0, p7); b10 = VMACLO(b10, u1, p7);
  y3 = VMACHI(y3, u0, p0); z3 = VMACHI(z3, u1, p0);
  y4 = VMACHI(y4, u0, p1); z4 = VMACHI(z4, u1, p1);
  y5 = VMACHI(y5, u0, p2); z5 = VMACHI(z5, u1, p2);
  y6 = VMACHI(y6, u0, p3); z6 = VMACHI(z6, u1, p3);
  y7 = VMACHI(y7, u0, p4); z7 = VMACHI(z7, u1, p4);
  y8 = VMACHI(y8, u0, p5); z8 = VMACHI(z8, u1, p5);
  y9 = VMACHI(y9, u0, p6); z9 = VMACHI(z9, u1, p6);
  y10 = VMACHI(y10, u0, p7); z10 = VMACHI(z10, u1, p7);
  a4 = VADD(VADD(a4, VSRA(a3, BRADIX)), VSHL(y3, BALIGN));
  b4 = VADD(VADD(b4, VSRA(b3, BRADIX)), VSHL(z3, BALIGN));

  u0 = VAND(VMACLO(zero, a4, montw), bmask);
  u1 = VAND(VMACLO(zero, b4, montw), bmask);
  a4 = VMACLO(a4, u0, p0); b4 = VMACLO(b4, u1, p0);
  a5 = VMACLO(a5, u0, p1); b5 = VMACLO(b5, u1, p1);
  a6 = VMACLO(a6, u0, p2); b6 = VMACLO(b6, u1, p2);
  a7 = VMACLO(a7, u0, p3); b7 = VMACLO(b7, u1, p3);
  a8 = VMACLO(a8, u0, p4); b8 = VMACLO(b8, u1, p4);
  a9 = VMACLO(a9, u0, p5); b9 = VMACLO(b9, u1, p5);
  a10 = VMACLO(a10, u0, p6); b10 = VMACLO(b10, u1, p6);
  a11 = VMACLO(a11, u0, p7); b11 = VMACLO(b11, u1, p7);
  y4 = VMACHI(y4, u0, p0); z4 = VMACHI(z4, u1, p0);
  y5 = VMACHI(y5, u0, p1); z5 = VMACHI(z5, u1, p1);
  y6 = VMACHI(y6, u0, p2); z6 = VMACHI(z6, u1, p2);
  y7 = VMACHI(y7, u0, p3); z7 = VMACHI(z7, u1, p3);
  y8 = VMACHI(y8, u0, p4); z8 = VMACHI(z8, u1, p4);
  y9 = VMACHI(y9, u0, p5); z9 = VMACHI(z9, u1, p5);
  y10 = VMACHI(y10, u0, p6); z10 = VMACHI(z10, u1, p6);
  y11 = VMACHI(y11, u0, p7); z11 = VMACHI(z11, u1, p7);
  a5 = VADD(VADD(a5, VSRA(a4, BRADIX)), VSHL(y4, BALIGN));
  b5 = VADD(VADD(b5, VSRA(b4, BRADIX)), VSHL(z4, BALIGN));

  u0 = VAND(VMACLO(zero, a5, montw), bmask);
  u1 = VAND(VMACLO(zero, b5, montw), bmask);
  a5 = VMACLO(a5, u0, p0); b5 = VMACLO(b5, u1, p0);
  a6 = VMACLO(a6, u0, p1); b6 = VMACLO(b6, u1, p1);
  a7 = VMACLO(a7, u0, p2); b7 = VMACLO(b7, u1, p2);
  a8 = VMACLO(a8, u0, p3); b8 = VMACLO(b8, u1, p3);
  a9 = VMACLO(a9, u0, p4); b9 = VMACLO(b9, u1, p4);
  a10 = VMACLO(a10, u0, p5); b10 = VMACLO(b10, u1, p5);
  a11 = VMACLO(a11, u0, p6); b11 = VMACLO(b11, u1, p6);
  a12 = VMACLO(a12, u0, p7); b12 = VMACLO(b12, u1, p7);
  y5 = VMACHI(y5, u0, p0); z5 = VMACHI(z5, u1, p0);
  y6 = VMACHI(y6, u0, p1); z6 = VMACHI(z6, u1, p1);
  y7 = VMACHI(y7, u0, p2); z7 = VMACHI(z7, u1, p2);
  y8 = VMACHI(y8, u0, p3); z8 = VMACHI(z8, u1, p3);
  y9 = VMACHI(y9, u0, p4); z9 = VMACHI(z9, u1, p4);
  y10 = VMACHI(y10, u0, p5); z10 = VMACHI(z10, u1, p5);
  y11 = VMACHI(y11, u0, p6); z11 = VMACHI(z11, u1, p6);
  y12 = VMACHI(y12, u0, p7); z12 = VMACHI(z12, u1, p7);
  a6 = VADD(VADD(a6, VSRA(a5, BRADIX)), VSHL(y5, BALIGN));
  b6 = VADD(VADD(b6, VSRA(b5, BRADIX)), VSHL(z5, BALIGN));

  u0 = VAND(VMACLO(zero, a6, montw), bmask);
  u1 = VAND(VMACLO(zero, b6, montw), bmask);
  a6 = VMACLO(a6, u0, p0); b6 = VMACLO(b6, u1, p0);
  a7 = VMACLO(a7, u0, p1); b7 = VMACLO(b7, u1, p1);
  a8 = VMACLO(a8, u0, p2); b8 = VMACLO(b8, u1, p2);
  a9 = VMACLO(a9, u0, p3); b9 = VMACLO(b9, u1, p3);
  a10 = VMACLO(a10, u0, p4); b10 = VMACLO(b10, u1, p4);
  a11 = VMACLO(a11, u0, p5); b11 = VMACLO(b11, u1, p5);
  a12 = VMACLO(a12, u0, p6); b12 = VMACLO(b12, u1, p6);
  a13 = VMACLO(a13, u0, p7); b13 = VMACLO(b13, u1, p7);
  y6 = VMACHI(y6, u0, p0); z6 = VMACHI(z6, u1, p0);
  y7 = VMACHI(y7, u0, p1); z7 = VMACHI(z7, u1, p1);
  y8 = VMACHI(y8, u0, p2); z8 = VMACHI(z8, u1, p2);
  y9 = VMACHI(y9, u0, p3); z9 = VMACHI(z9, u1, p3);
  y10 = VMACHI(y10, u0, p4); z10 = VMACHI(z10, u1, p4);
  y11 = VMACHI(y11, u0, p5); z11 = VMACHI(z11, u1, p5);
  y12 = VMACHI(y12, u0, p6); z12 = VMACHI(z12, u1, p6);
  y13 = VMACHI(y13, u0, p7); z13 = VMACHI(z13, u1, p7);
  a7 = VADD(VADD(a7, VSRA(a6, BRADIX)), VSHL(y6, BALIGN));
  b7 = VADD(VADD(b7, VSRA(b6, BRADIX)), VSHL(z6, BALIGN));

  u0 = VAND(VMACLO(zero, a7, montw), bmask);
  u1 = VAND(VMACLO(zero, b7, montw), bmask);
  a7 = VMACLO(a7, u0, p0); b7 = VMACLO(b7, u1, p0);
  a8 = VMACLO(a8, u0, p1); b8 = VMACLO(b8, u1, p1);
  a9 = VMACLO(a9, u0, p2); b9 = VMACLO(b9, u1, p2);
  a10 = VMACLO(a10, u0, p3); b10 = VMACLO(b10, u1, p3);
  a11 = VMACLO(a11, u0, p4); b11 = VMACLO(b11, u1, p4);
  a12 = VMACLO(a12, u0, p5); b12 = VMACLO(b12, u1, p5);
  a13 = VMACLO(a13, u0, p6); b13 = VMACLO(b13, u1, p6);
  a14 = VMACLO(a14, u0, p7); b14 = VMACLO(b14, u1, p7);
  y7 = VMACHI(y7, u0, p0); z7 = VMACHI(z7, u1, p0);
  y8 = VMACHI(y8, u0, p1); z8 = VMACHI(z8, u1, p1);
  y9 = VMACHI(y9, u0, p2); z9 = VMACHI(z9, u1, p2);
  y10 = VMACHI(y10, u0, p3); z10 = VMACHI(z10, u1, p3);
  y11 = VMACHI(y11, u0, p4); z11 = VMACHI(z11, u1, p4);
  y12 = VMACHI(y12, u0, p5); z12 = VMACHI(z12, u1, p5);
  y13 = VMACHI(y13, u0, p6); z13 = VMACHI(z13, u1, p6);
  y14 = VMACHI(y14, u0, p7); z14 = VMACHI(z14, u1, p7);
  a8 = VADD(VADD(a8, VSRA(a7, BRADIX)), VSHL(y7, BALIGN));
  b8 = VADD(VADD(b8, VSRA(b7, BRADIX)), VSHL(z7, BALIGN));

  a9 = VADD(a9, VSHL(y8, BALIGN)); b9 = VADD(b9, VSHL(z8, BALIGN));
  a10 = VADD(a10, VSHL(y9, BALIGN)); b10 = VADD(b10, VSHL(z9, BALIGN));
  a11 = VADD(a11, VSHL(y10, BALIGN)); b11 = VADD(b11, VSHL(z10, BALIGN));
  a12 = VADD(a12, VSHL(y11, BALIGN)); b12 = VADD(b12, VSHL(z11, BALIGN));
  a13 = VADD(a13, VSHL(y12, BALIGN)); b13 = VADD(b13, VSHL(z12, BALIGN));
  a14 = VADD(a14, VSHL(y13, BALIGN)); b14 = VADD(b14, VSHL(z13, BALIGN));
  a15 = VADD(a15, VSHL(y14, BALIGN)); b15 = VADD(b15, VSHL(z14, BALIGN));

  cc512i ra0 = a8, ra1 = a9, ra2 = a10, ra3 = a11;
  cc512i ra4 = a12, ra5 = a13, ra6 = a14, ra7 = a15;
  cc512i rb0 = b8, rb1 = b9, rb2 = b10, rb3 = b11;
  cc512i rb4 = b12, rb5 = b13, rb6 = b14, rb7 = b15;

  ra1 = VADD(ra1, VSRA(ra0, BRADIX)); ra0 = VAND(ra0, bmask);
  rb1 = VADD(rb1, VSRA(rb0, BRADIX)); rb0 = VAND(rb0, bmask);
  ra2 = VADD(ra2, VSRA(ra1, BRADIX)); ra1 = VAND(ra1, bmask);
  rb2 = VADD(rb2, VSRA(rb1, BRADIX)); rb1 = VAND(rb1, bmask);
  ra3 = VADD(ra3, VSRA(ra2, BRADIX)); ra2 = VAND(ra2, bmask);
  rb3 = VADD(rb3, VSRA(rb2, BRADIX)); rb2 = VAND(rb2, bmask);
  ra4 = VADD(ra4, VSRA(ra3, BRADIX)); ra3 = VAND(ra3, bmask);
  rb4 = VADD(rb4, VSRA(rb3, BRADIX)); rb3 = VAND(rb3, bmask);
  ra5 = VADD(ra5, VSRA(ra4, BRADIX)); ra4 = VAND(ra4, bmask);
  rb5 = VADD(rb5, VSRA(rb4, BRADIX)); rb4 = VAND(rb4, bmask);
  ra6 = VADD(ra6, VSRA(ra5, BRADIX)); ra5 = VAND(ra5, bmask);
  rb6 = VADD(rb6, VSRA(rb5, BRADIX)); rb5 = VAND(rb5, bmask);
  ra7 = VADD(ra7, VSRA(ra6, BRADIX)); ra6 = VAND(ra6, bmask); ra7 = VAND(ra7, bmask);
  rb7 = VADD(rb7, VSRA(rb6, BRADIX)); rb6 = VAND(rb6, bmask); rb7 = VAND(rb7, bmask);

  out0[0] = ra0; out0[1] = ra1; out0[2] = ra2; out0[3] = ra3;
  out0[4] = ra4; out0[5] = ra5; out0[6] = ra6; out0[7] = ra7;
  out1[0] = rb0; out1[1] = rb1; out1[2] = rb2; out1[3] = rb3;
  out1[4] = rb4; out1[5] = rb5; out1[6] = rb6; out1[7] = rb7;
}

template <int K>
void BM_Fp8x1MulV2K(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * K;
  Fp8x1KBenchData data(n_packs);

  for (auto _ : state) {
    fp8x1_mul_v2_k<K>(data.wide.data(), data.a.data(), data.b.data(), n_blocks);
    benchmark::DoNotOptimize(data.wide.data());
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

template <int K>
void BM_Fp8x1RedcK(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * K;
  Fp8x1KBenchData data(n_packs);

  for (auto _ : state) {
    fp8x1_redc_k<K>(data.out.data(), data.wide.data(), n_blocks);
    benchmark::DoNotOptimize(data.out.data());
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

template <int K>
void BM_Fp8x1MulRedcK(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * K;
  Fp8x1KBenchData data(n_packs);

  for (auto _ : state) {
    fp8x1_mul_redc_k<K>(data.out.data(), data.a.data(), data.b.data(), n_blocks);
    benchmark::DoNotOptimize(data.out.data());
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1RedcCore(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1KBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      UFpPack8x1::redc_core(data.lazy[i], data.wide[i]);
    }
    benchmark::DoNotOptimize(data.lazy.data());
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1Normalize(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1KBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      UFpPack8x1::normalize(data.out[i], data.lazy[i]);
    }
    benchmark::DoNotOptimize(data.out.data());
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

template <int Chains>
void BM_Fp8x1MulRedc_SeqChains(benchmark::State& state)
{
  static_assert(Chains >= 1 && Chains <= 8, "SeqChains supports Chains in [1,8].");
  constexpr std::size_t kUnroll = 8;
  Fp8x1MulRedcStreamData<Chains> data{};

  for (auto _ : state) {
    if constexpr (Chains >= 1) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul_v2(data.wide[0], data.cur[0], data.mul_by[0]);
        UFpPack8x1::redc(data.out[0], data.wide[0]);
        data.cur[0] = data.out[0];
      }
    }
    if constexpr (Chains >= 2) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul_v2(data.wide[1], data.cur[1], data.mul_by[1]);
        UFpPack8x1::redc(data.out[1], data.wide[1]);
        data.cur[1] = data.out[1];
      }
    }
    if constexpr (Chains >= 3) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul_v2(data.wide[2], data.cur[2], data.mul_by[2]);
        UFpPack8x1::redc(data.out[2], data.wide[2]);
        data.cur[2] = data.out[2];
      }
    }
    if constexpr (Chains >= 4) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul_v2(data.wide[3], data.cur[3], data.mul_by[3]);
        UFpPack8x1::redc(data.out[3], data.wide[3]);
        data.cur[3] = data.out[3];
      }
    }
    if constexpr (Chains >= 5) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul_v2(data.wide[4], data.cur[4], data.mul_by[4]);
        UFpPack8x1::redc(data.out[4], data.wide[4]);
        data.cur[4] = data.out[4];
      }
    }
    if constexpr (Chains >= 6) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul_v2(data.wide[5], data.cur[5], data.mul_by[5]);
        UFpPack8x1::redc(data.out[5], data.wide[5]);
        data.cur[5] = data.out[5];
      }
    }
    if constexpr (Chains >= 7) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul_v2(data.wide[6], data.cur[6], data.mul_by[6]);
        UFpPack8x1::redc(data.out[6], data.wide[6]);
        data.cur[6] = data.out[6];
      }
    }
    if constexpr (Chains >= 8) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul_v2(data.wide[7], data.cur[7], data.mul_by[7]);
        UFpPack8x1::redc(data.out[7], data.wide[7]);
        data.cur[7] = data.out[7];
      }
    }
  }

  cc512i sink = VZERO();
  for (int c = 0; c < Chains; ++c) {
    sink = VADD(sink, data.cur[c][0]);
  }
  benchmark::DoNotOptimize(sink);

  const std::size_t ops_per_iter = static_cast<std::size_t>(Chains) * kUnroll;
  state.SetItemsProcessed(state.iterations() * ops_per_iter);
  state.counters["mul_redc_per_iter"] = static_cast<double>(ops_per_iter);
  state.counters["chains"] = static_cast<double>(Chains);
  state.counters["unroll"] = static_cast<double>(kUnroll);
}

template <int Chains>
void BM_Fp8x1MulRedc_HandInterleavedChains(benchmark::State& state)
{
  static_assert(Chains >= 1 && Chains <= 8, "HandInterleaved supports Chains in [1,8].");
  constexpr std::size_t kUnroll = 8;
  Fp8x1MulRedcStreamData<Chains> data{};

  for (auto _ : state) {
    for (std::size_t i = 0; i < kUnroll; ++i) {
      if constexpr (Chains >= 1) {
        FpPack8x1::mul_v2(data.wide[0], data.cur[0], data.mul_by[0]);
        UFpPack8x1::redc(data.out[0], data.wide[0]);
        data.cur[0] = data.out[0];
      }
      if constexpr (Chains >= 2) {
        FpPack8x1::mul_v2(data.wide[1], data.cur[1], data.mul_by[1]);
        UFpPack8x1::redc(data.out[1], data.wide[1]);
        data.cur[1] = data.out[1];
      }
      if constexpr (Chains >= 3) {
        FpPack8x1::mul_v2(data.wide[2], data.cur[2], data.mul_by[2]);
        UFpPack8x1::redc(data.out[2], data.wide[2]);
        data.cur[2] = data.out[2];
      }
      if constexpr (Chains >= 4) {
        FpPack8x1::mul_v2(data.wide[3], data.cur[3], data.mul_by[3]);
        UFpPack8x1::redc(data.out[3], data.wide[3]);
        data.cur[3] = data.out[3];
      }
      if constexpr (Chains >= 5) {
        FpPack8x1::mul_v2(data.wide[4], data.cur[4], data.mul_by[4]);
        UFpPack8x1::redc(data.out[4], data.wide[4]);
        data.cur[4] = data.out[4];
      }
      if constexpr (Chains >= 6) {
        FpPack8x1::mul_v2(data.wide[5], data.cur[5], data.mul_by[5]);
        UFpPack8x1::redc(data.out[5], data.wide[5]);
        data.cur[5] = data.out[5];
      }
      if constexpr (Chains >= 7) {
        FpPack8x1::mul_v2(data.wide[6], data.cur[6], data.mul_by[6]);
        UFpPack8x1::redc(data.out[6], data.wide[6]);
        data.cur[6] = data.out[6];
      }
      if constexpr (Chains >= 8) {
        FpPack8x1::mul_v2(data.wide[7], data.cur[7], data.mul_by[7]);
        UFpPack8x1::redc(data.out[7], data.wide[7]);
        data.cur[7] = data.out[7];
      }
    }
  }

  cc512i sink = VZERO();
  for (int c = 0; c < Chains; ++c) {
    sink = VADD(sink, data.cur[c][0]);
  }
  benchmark::DoNotOptimize(sink);

  const std::size_t ops_per_iter = static_cast<std::size_t>(Chains) * kUnroll;
  state.SetItemsProcessed(state.iterations() * ops_per_iter);
  state.counters["mul_redc_per_iter"] = static_cast<double>(ops_per_iter);
  state.counters["chains"] = static_cast<double>(Chains);
  state.counters["unroll"] = static_cast<double>(kUnroll);
}

template <int Chains>
void BM_Fp8x1RedcCore_SeqChains(benchmark::State& state)
{
  static_assert(Chains >= 1 && Chains <= 8, "SeqChains supports Chains in [1,8].");
  constexpr std::size_t kUnroll = 8;
  Fp8x1RedcCoreStreamData<Chains> data{};

  const cc512i one = VSET1(1);
  for (auto _ : state) {
    if constexpr (Chains >= 1) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        UFpPack8x1::redc_core(data.lazy[0], data.wide[0]);
        data.wide[0][0] = VADD(data.wide[0][0], VADD(data.lazy[0][0], one));
      }
    }
    if constexpr (Chains >= 2) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        UFpPack8x1::redc_core(data.lazy[1], data.wide[1]);
        data.wide[1][0] = VADD(data.wide[1][0], VADD(data.lazy[1][0], one));
      }
    }
    if constexpr (Chains >= 3) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        UFpPack8x1::redc_core(data.lazy[2], data.wide[2]);
        data.wide[2][0] = VADD(data.wide[2][0], VADD(data.lazy[2][0], one));
      }
    }
    if constexpr (Chains >= 4) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        UFpPack8x1::redc_core(data.lazy[3], data.wide[3]);
        data.wide[3][0] = VADD(data.wide[3][0], VADD(data.lazy[3][0], one));
      }
    }
    if constexpr (Chains >= 5) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        UFpPack8x1::redc_core(data.lazy[4], data.wide[4]);
        data.wide[4][0] = VADD(data.wide[4][0], VADD(data.lazy[4][0], one));
      }
    }
    if constexpr (Chains >= 6) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        UFpPack8x1::redc_core(data.lazy[5], data.wide[5]);
        data.wide[5][0] = VADD(data.wide[5][0], VADD(data.lazy[5][0], one));
      }
    }
    if constexpr (Chains >= 7) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        UFpPack8x1::redc_core(data.lazy[6], data.wide[6]);
        data.wide[6][0] = VADD(data.wide[6][0], VADD(data.lazy[6][0], one));
      }
    }
    if constexpr (Chains >= 8) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        UFpPack8x1::redc_core(data.lazy[7], data.wide[7]);
        data.wide[7][0] = VADD(data.wide[7][0], VADD(data.lazy[7][0], one));
      }
    }
  }

  cc512i sink = VZERO();
  for (int c = 0; c < Chains; ++c) {
    sink = VADD(sink, data.lazy[c][0]);
  }
  benchmark::DoNotOptimize(sink);

  const std::size_t ops_per_iter = static_cast<std::size_t>(Chains) * kUnroll;
  state.SetItemsProcessed(state.iterations() * ops_per_iter);
  state.counters["redc_core_per_iter"] = static_cast<double>(ops_per_iter);
  state.counters["chains"] = static_cast<double>(Chains);
  state.counters["unroll"] = static_cast<double>(kUnroll);
}

template <int Chains>
void BM_Fp8x1RedcCore_HandInterleavedChains(benchmark::State& state)
{
  static_assert((Chains == 2) || (Chains == 4) || (Chains == 8),
                "HandInterleaved benchmark is explicitly specialized for Chains in {2,4,8}.");
}

template <>
void BM_Fp8x1RedcCore_HandInterleavedChains<2>(benchmark::State& state)
{
  constexpr std::size_t kUnroll = 8;
  Fp8x1RedcCoreStreamData<2> data{};
  const cc512i one = VSET1(1);

  for (auto _ : state) {
    for (std::size_t i = 0; i < kUnroll; ++i) {
      redc_core_hand_interleaved_k2_once(data.lazy[0], data.lazy[1], data.wide[0], data.wide[1]);
      data.wide[0][0] = VADD(data.wide[0][0], VADD(data.lazy[0][0], one));
      data.wide[1][0] = VADD(data.wide[1][0], VADD(data.lazy[1][0], one));
    }
  }

  cc512i sink = VZERO();
  sink = VADD(sink, data.lazy[0][0]);
  sink = VADD(sink, data.lazy[1][0]);
  benchmark::DoNotOptimize(sink);

  const std::size_t ops_per_iter = 2 * kUnroll;
  state.SetItemsProcessed(state.iterations() * ops_per_iter);
  state.counters["redc_core_per_iter"] = static_cast<double>(ops_per_iter);
  state.counters["chains"] = 2.0;
  state.counters["unroll"] = static_cast<double>(kUnroll);
}

template <>
void BM_Fp8x1RedcCore_HandInterleavedChains<4>(benchmark::State& state)
{
  constexpr std::size_t kUnroll = 8;
  Fp8x1RedcCoreStreamData<4> data{};
  const cc512i one = VSET1(1);

  for (auto _ : state) {
    for (std::size_t i = 0; i < kUnroll; ++i) {
      redc_core_hand_interleaved_k2_once(data.lazy[0], data.lazy[1], data.wide[0], data.wide[1]);
      redc_core_hand_interleaved_k2_once(data.lazy[2], data.lazy[3], data.wide[2], data.wide[3]);
      data.wide[0][0] = VADD(data.wide[0][0], VADD(data.lazy[0][0], one));
      data.wide[1][0] = VADD(data.wide[1][0], VADD(data.lazy[1][0], one));
      data.wide[2][0] = VADD(data.wide[2][0], VADD(data.lazy[2][0], one));
      data.wide[3][0] = VADD(data.wide[3][0], VADD(data.lazy[3][0], one));
    }
  }

  cc512i sink = VZERO();
  sink = VADD(sink, data.lazy[0][0]);
  sink = VADD(sink, data.lazy[1][0]);
  sink = VADD(sink, data.lazy[2][0]);
  sink = VADD(sink, data.lazy[3][0]);
  benchmark::DoNotOptimize(sink);

  const std::size_t ops_per_iter = 4 * kUnroll;
  state.SetItemsProcessed(state.iterations() * ops_per_iter);
  state.counters["redc_core_per_iter"] = static_cast<double>(ops_per_iter);
  state.counters["chains"] = 4.0;
  state.counters["unroll"] = static_cast<double>(kUnroll);
}

template <>
void BM_Fp8x1RedcCore_HandInterleavedChains<8>(benchmark::State& state)
{
  constexpr std::size_t kUnroll = 8;
  Fp8x1RedcCoreStreamData<8> data{};
  const cc512i one = VSET1(1);

  for (auto _ : state) {
    for (std::size_t i = 0; i < kUnroll; ++i) {
      redc_core_hand_interleaved_k2_once(data.lazy[0], data.lazy[1], data.wide[0], data.wide[1]);
      redc_core_hand_interleaved_k2_once(data.lazy[2], data.lazy[3], data.wide[2], data.wide[3]);
      redc_core_hand_interleaved_k2_once(data.lazy[4], data.lazy[5], data.wide[4], data.wide[5]);
      redc_core_hand_interleaved_k2_once(data.lazy[6], data.lazy[7], data.wide[6], data.wide[7]);
      data.wide[0][0] = VADD(data.wide[0][0], VADD(data.lazy[0][0], one));
      data.wide[1][0] = VADD(data.wide[1][0], VADD(data.lazy[1][0], one));
      data.wide[2][0] = VADD(data.wide[2][0], VADD(data.lazy[2][0], one));
      data.wide[3][0] = VADD(data.wide[3][0], VADD(data.lazy[3][0], one));
      data.wide[4][0] = VADD(data.wide[4][0], VADD(data.lazy[4][0], one));
      data.wide[5][0] = VADD(data.wide[5][0], VADD(data.lazy[5][0], one));
      data.wide[6][0] = VADD(data.wide[6][0], VADD(data.lazy[6][0], one));
      data.wide[7][0] = VADD(data.wide[7][0], VADD(data.lazy[7][0], one));
    }
  }

  cc512i sink = VZERO();
  sink = VADD(sink, data.lazy[0][0]);
  sink = VADD(sink, data.lazy[1][0]);
  sink = VADD(sink, data.lazy[2][0]);
  sink = VADD(sink, data.lazy[3][0]);
  sink = VADD(sink, data.lazy[4][0]);
  sink = VADD(sink, data.lazy[5][0]);
  sink = VADD(sink, data.lazy[6][0]);
  sink = VADD(sink, data.lazy[7][0]);
  benchmark::DoNotOptimize(sink);

  const std::size_t ops_per_iter = 8 * kUnroll;
  state.SetItemsProcessed(state.iterations() * ops_per_iter);
  state.counters["redc_core_per_iter"] = static_cast<double>(ops_per_iter);
  state.counters["chains"] = 8.0;
  state.counters["unroll"] = static_cast<double>(kUnroll);
}

template <bool DoLo, bool DoHi>
void BM_IFMA52Latency_DepChain(benchmark::State& state)
{
  static_assert(DoLo || DoHi, "At least one IFMA op must be enabled.");
  constexpr std::size_t kLatencyUnroll = 8;

  const cc512i a = VSET1(0x0000FFFFFFFFFFFFULL);
  const cc512i b = VSET1(0x00000000FFFFFFFFULL);
  cc512i x = VSET1(1);

  for (auto _ : state) {
    for (std::size_t i = 0; i < kLatencyUnroll; ++i) {
      if constexpr (DoLo) {
        x = VMACLO(x, a, b);
      }
      if constexpr (DoHi) {
        x = VMACHI(x, a, b);
      }
    }
    benchmark::DoNotOptimize(x);
  }

  const std::size_t madd_per_iter =
      kLatencyUnroll * ((DoLo ? 1u : 0u) + (DoHi ? 1u : 0u));
  state.SetItemsProcessed(state.iterations() * madd_per_iter);
  state.counters["madd52_per_iter"] = static_cast<double>(madd_per_iter);
  state.counters["unroll"] = static_cast<double>(kLatencyUnroll);
}

template <int Streams, bool DoLo, bool DoHi>
void BM_IFMA52Throughput_Indep(benchmark::State& state)
{
  static_assert(Streams >= 1 && Streams <= 16, "Streams must be in [1,16].");
  static_assert(DoLo || DoHi, "At least one IFMA op must be enabled.");

  const cc512i a = VSET1(0x0000FFFFFFFFFFFFULL);
  const cc512i b = VSET1(0x00000000FFFFFFFFULL);

  cc512i x0 = VSET1(1), x1 = VSET1(2), x2 = VSET1(3), x3 = VSET1(4);
  cc512i x4 = VSET1(5), x5 = VSET1(6), x6 = VSET1(7), x7 = VSET1(8);
  cc512i x8 = VSET1(9), x9 = VSET1(10), x10 = VSET1(11), x11 = VSET1(12);
  cc512i x12 = VSET1(13), x13 = VSET1(14), x14 = VSET1(15), x15 = VSET1(16);

  for (auto _ : state) {
    if constexpr (DoLo) {
      if constexpr (Streams >= 1) { x0 = VMACLO(x0, a, b); }
      if constexpr (Streams >= 2) { x1 = VMACLO(x1, a, b); }
      if constexpr (Streams >= 3) { x2 = VMACLO(x2, a, b); }
      if constexpr (Streams >= 4) { x3 = VMACLO(x3, a, b); }
      if constexpr (Streams >= 5) { x4 = VMACLO(x4, a, b); }
      if constexpr (Streams >= 6) { x5 = VMACLO(x5, a, b); }
      if constexpr (Streams >= 7) { x6 = VMACLO(x6, a, b); }
      if constexpr (Streams >= 8) { x7 = VMACLO(x7, a, b); }
      if constexpr (Streams >= 9) { x8 = VMACLO(x8, a, b); }
      if constexpr (Streams >= 10) { x9 = VMACLO(x9, a, b); }
      if constexpr (Streams >= 11) { x10 = VMACLO(x10, a, b); }
      if constexpr (Streams >= 12) { x11 = VMACLO(x11, a, b); }
      if constexpr (Streams >= 13) { x12 = VMACLO(x12, a, b); }
      if constexpr (Streams >= 14) { x13 = VMACLO(x13, a, b); }
      if constexpr (Streams >= 15) { x14 = VMACLO(x14, a, b); }
      if constexpr (Streams >= 16) { x15 = VMACLO(x15, a, b); }
    }
    if constexpr (DoHi) {
      if constexpr (Streams >= 1) { x0 = VMACHI(x0, a, b); }
      if constexpr (Streams >= 2) { x1 = VMACHI(x1, a, b); }
      if constexpr (Streams >= 3) { x2 = VMACHI(x2, a, b); }
      if constexpr (Streams >= 4) { x3 = VMACHI(x3, a, b); }
      if constexpr (Streams >= 5) { x4 = VMACHI(x4, a, b); }
      if constexpr (Streams >= 6) { x5 = VMACHI(x5, a, b); }
      if constexpr (Streams >= 7) { x6 = VMACHI(x6, a, b); }
      if constexpr (Streams >= 8) { x7 = VMACHI(x7, a, b); }
      if constexpr (Streams >= 9) { x8 = VMACHI(x8, a, b); }
      if constexpr (Streams >= 10) { x9 = VMACHI(x9, a, b); }
      if constexpr (Streams >= 11) { x10 = VMACHI(x10, a, b); }
      if constexpr (Streams >= 12) { x11 = VMACHI(x11, a, b); }
      if constexpr (Streams >= 13) { x12 = VMACHI(x12, a, b); }
      if constexpr (Streams >= 14) { x13 = VMACHI(x13, a, b); }
      if constexpr (Streams >= 15) { x14 = VMACHI(x14, a, b); }
      if constexpr (Streams >= 16) { x15 = VMACHI(x15, a, b); }
    }
  }

  cc512i sink = VZERO();
  if constexpr (Streams >= 1) { sink = VADD(sink, x0); }
  if constexpr (Streams >= 2) { sink = VADD(sink, x1); }
  if constexpr (Streams >= 3) { sink = VADD(sink, x2); }
  if constexpr (Streams >= 4) { sink = VADD(sink, x3); }
  if constexpr (Streams >= 5) { sink = VADD(sink, x4); }
  if constexpr (Streams >= 6) { sink = VADD(sink, x5); }
  if constexpr (Streams >= 7) { sink = VADD(sink, x6); }
  if constexpr (Streams >= 8) { sink = VADD(sink, x7); }
  if constexpr (Streams >= 9) { sink = VADD(sink, x8); }
  if constexpr (Streams >= 10) { sink = VADD(sink, x9); }
  if constexpr (Streams >= 11) { sink = VADD(sink, x10); }
  if constexpr (Streams >= 12) { sink = VADD(sink, x11); }
  if constexpr (Streams >= 13) { sink = VADD(sink, x12); }
  if constexpr (Streams >= 14) { sink = VADD(sink, x13); }
  if constexpr (Streams >= 15) { sink = VADD(sink, x14); }
  if constexpr (Streams >= 16) { sink = VADD(sink, x15); }
  benchmark::DoNotOptimize(sink);

  const std::size_t madd_per_iter =
      static_cast<std::size_t>(Streams) * ((DoLo ? 1u : 0u) + (DoHi ? 1u : 0u));
  state.SetItemsProcessed(state.iterations() * madd_per_iter);
  state.counters["madd52_per_iter"] = static_cast<double>(madd_per_iter);
}

template <int Chains>
void BM_IFMA52LoHi_SeqChains(benchmark::State& state)
{
  static_assert(Chains >= 2 && Chains <= 8, "SeqChains PoC supports Chains in [2,8].");
  constexpr std::size_t kUnroll = 8;

  const cc512i a = VSET1(0x0000FFFFFFFFFFFFULL);
  const cc512i b = VSET1(0x00000000FFFFFFFFULL);
  cc512i x0 = VSET1(1), x1 = VSET1(2), x2 = VSET1(3), x3 = VSET1(4);
  cc512i x4 = VSET1(5), x5 = VSET1(6), x6 = VSET1(7), x7 = VSET1(8);

  for (auto _ : state) {
    for (std::size_t i = 0; i < kUnroll; ++i) {
      x0 = VMACLO(x0, a, b);
      x0 = VMACHI(x0, a, b);
    }
    for (std::size_t i = 0; i < kUnroll; ++i) {
      x1 = VMACLO(x1, a, b);
      x1 = VMACHI(x1, a, b);
    }
    if constexpr (Chains >= 3) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        x2 = VMACLO(x2, a, b);
        x2 = VMACHI(x2, a, b);
      }
    }
    if constexpr (Chains >= 4) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        x3 = VMACLO(x3, a, b);
        x3 = VMACHI(x3, a, b);
      }
    }
    if constexpr (Chains >= 5) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        x4 = VMACLO(x4, a, b);
        x4 = VMACHI(x4, a, b);
      }
    }
    if constexpr (Chains >= 6) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        x5 = VMACLO(x5, a, b);
        x5 = VMACHI(x5, a, b);
      }
    }
    if constexpr (Chains >= 7) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        x6 = VMACLO(x6, a, b);
        x6 = VMACHI(x6, a, b);
      }
    }
    if constexpr (Chains >= 8) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        x7 = VMACLO(x7, a, b);
        x7 = VMACHI(x7, a, b);
      }
    }
  }

  cc512i sink = VADD(x0, x1);
  if constexpr (Chains >= 3) {
    sink = VADD(sink, x2);
  }
  if constexpr (Chains >= 4) {
    sink = VADD(sink, x3);
  }
  if constexpr (Chains >= 5) {
    sink = VADD(sink, x4);
  }
  if constexpr (Chains >= 6) {
    sink = VADD(sink, x5);
  }
  if constexpr (Chains >= 7) {
    sink = VADD(sink, x6);
  }
  if constexpr (Chains >= 8) {
    sink = VADD(sink, x7);
  }
  benchmark::DoNotOptimize(sink);

  const std::size_t madd_per_iter = static_cast<std::size_t>(Chains) * kUnroll * 2;
  state.SetItemsProcessed(state.iterations() * madd_per_iter);
  state.counters["madd52_per_iter"] = static_cast<double>(madd_per_iter);
  state.counters["chains"] = static_cast<double>(Chains);
  state.counters["unroll"] = static_cast<double>(kUnroll);
}

template <int Chains>
void BM_IFMA52LoHi_HandInterleavedChains(benchmark::State& state)
{
  static_assert(Chains >= 2 && Chains <= 8, "HandInterleaved PoC supports Chains in [2,8].");
  constexpr std::size_t kUnroll = 8;

  const cc512i a = VSET1(0x0000FFFFFFFFFFFFULL);
  const cc512i b = VSET1(0x00000000FFFFFFFFULL);
  cc512i x0 = VSET1(1), x1 = VSET1(2), x2 = VSET1(3), x3 = VSET1(4);
  cc512i x4 = VSET1(5), x5 = VSET1(6), x6 = VSET1(7), x7 = VSET1(8);

  for (auto _ : state) {
    for (std::size_t i = 0; i < kUnroll; ++i) {
      // Hand-stamped cross-chain issue pattern.
      x0 = VMACLO(x0, a, b);
      x1 = VMACLO(x1, a, b);
      if constexpr (Chains >= 3) { x2 = VMACLO(x2, a, b); }
      if constexpr (Chains >= 4) { x3 = VMACLO(x3, a, b); }
      if constexpr (Chains >= 5) { x4 = VMACLO(x4, a, b); }
      if constexpr (Chains >= 6) { x5 = VMACLO(x5, a, b); }
      if constexpr (Chains >= 7) { x6 = VMACLO(x6, a, b); }
      if constexpr (Chains >= 8) { x7 = VMACLO(x7, a, b); }

      x0 = VMACHI(x0, a, b);
      x1 = VMACHI(x1, a, b);
      if constexpr (Chains >= 3) { x2 = VMACHI(x2, a, b); }
      if constexpr (Chains >= 4) { x3 = VMACHI(x3, a, b); }
      if constexpr (Chains >= 5) { x4 = VMACHI(x4, a, b); }
      if constexpr (Chains >= 6) { x5 = VMACHI(x5, a, b); }
      if constexpr (Chains >= 7) { x6 = VMACHI(x6, a, b); }
      if constexpr (Chains >= 8) { x7 = VMACHI(x7, a, b); }
    }
  }

  cc512i sink = VADD(x0, x1);
  if constexpr (Chains >= 3) {
    sink = VADD(sink, x2);
  }
  if constexpr (Chains >= 4) {
    sink = VADD(sink, x3);
  }
  if constexpr (Chains >= 5) {
    sink = VADD(sink, x4);
  }
  if constexpr (Chains >= 6) {
    sink = VADD(sink, x5);
  }
  if constexpr (Chains >= 7) {
    sink = VADD(sink, x6);
  }
  if constexpr (Chains >= 8) {
    sink = VADD(sink, x7);
  }
  benchmark::DoNotOptimize(sink);

  const std::size_t madd_per_iter = static_cast<std::size_t>(Chains) * kUnroll * 2;
  state.SetItemsProcessed(state.iterations() * madd_per_iter);
  state.counters["madd52_per_iter"] = static_cast<double>(madd_per_iter);
  state.counters["chains"] = static_cast<double>(Chains);
  state.counters["unroll"] = static_cast<double>(kUnroll);
}

void BM_IFMA52Latency_Lo_DepChain(benchmark::State& state)
{
  BM_IFMA52Latency_DepChain<true, false>(state);
}

void BM_IFMA52Latency_Hi_DepChain(benchmark::State& state)
{
  BM_IFMA52Latency_DepChain<false, true>(state);
}

void BM_IFMA52Latency_LoHi_DepChain(benchmark::State& state)
{
  BM_IFMA52Latency_DepChain<true, true>(state);
}

void BM_IFMA52Throughput_Lo_Indep4(benchmark::State& state) { BM_IFMA52Throughput_Indep<4, true, false>(state); }
void BM_IFMA52Throughput_Lo_Indep8(benchmark::State& state) { BM_IFMA52Throughput_Indep<8, true, false>(state); }
void BM_IFMA52Throughput_Lo_Indep12(benchmark::State& state) { BM_IFMA52Throughput_Indep<12, true, false>(state); }
void BM_IFMA52Throughput_Lo_Indep16(benchmark::State& state) { BM_IFMA52Throughput_Indep<16, true, false>(state); }

void BM_IFMA52Throughput_Hi_Indep4(benchmark::State& state) { BM_IFMA52Throughput_Indep<4, false, true>(state); }
void BM_IFMA52Throughput_Hi_Indep8(benchmark::State& state) { BM_IFMA52Throughput_Indep<8, false, true>(state); }
void BM_IFMA52Throughput_Hi_Indep12(benchmark::State& state) { BM_IFMA52Throughput_Indep<12, false, true>(state); }
void BM_IFMA52Throughput_Hi_Indep16(benchmark::State& state) { BM_IFMA52Throughput_Indep<16, false, true>(state); }

void BM_IFMA52Throughput_LoHi_Indep4(benchmark::State& state) { BM_IFMA52Throughput_Indep<4, true, true>(state); }
void BM_IFMA52Throughput_LoHi_Indep8(benchmark::State& state) { BM_IFMA52Throughput_Indep<8, true, true>(state); }
void BM_IFMA52Throughput_LoHi_Indep12(benchmark::State& state) { BM_IFMA52Throughput_Indep<12, true, true>(state); }
void BM_IFMA52Throughput_LoHi_Indep16(benchmark::State& state) { BM_IFMA52Throughput_Indep<16, true, true>(state); }

BENCHMARK_TEMPLATE(BM_Fp8x1MulV2K, 1)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1MulV2K, 2)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1MulV2K, 4)->Arg(256);

BENCHMARK_TEMPLATE(BM_Fp8x1RedcK, 1)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1RedcK, 2)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1RedcK, 4)->Arg(256);

BENCHMARK_TEMPLATE(BM_Fp8x1MulRedcK, 1)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1MulRedcK, 2)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1MulRedcK, 4)->Arg(256);

BENCHMARK(BM_Fp8x1RedcCore)->Arg(256);
BENCHMARK(BM_Fp8x1Normalize)->Arg(256);

BENCHMARK(BM_IFMA52Latency_Lo_DepChain);
BENCHMARK(BM_IFMA52Latency_Hi_DepChain);
BENCHMARK(BM_IFMA52Latency_LoHi_DepChain);

BENCHMARK(BM_IFMA52Throughput_Lo_Indep4);
BENCHMARK(BM_IFMA52Throughput_Lo_Indep8);
BENCHMARK(BM_IFMA52Throughput_Lo_Indep12);
BENCHMARK(BM_IFMA52Throughput_Lo_Indep16);

BENCHMARK(BM_IFMA52Throughput_Hi_Indep4);
BENCHMARK(BM_IFMA52Throughput_Hi_Indep8);
BENCHMARK(BM_IFMA52Throughput_Hi_Indep12);
BENCHMARK(BM_IFMA52Throughput_Hi_Indep16);

BENCHMARK(BM_IFMA52Throughput_LoHi_Indep4);
BENCHMARK(BM_IFMA52Throughput_LoHi_Indep8);
BENCHMARK(BM_IFMA52Throughput_LoHi_Indep12);
BENCHMARK(BM_IFMA52Throughput_LoHi_Indep16);

BENCHMARK_TEMPLATE(BM_IFMA52LoHi_SeqChains, 2);
BENCHMARK_TEMPLATE(BM_IFMA52LoHi_HandInterleavedChains, 2);
BENCHMARK_TEMPLATE(BM_IFMA52LoHi_SeqChains, 4);
BENCHMARK_TEMPLATE(BM_IFMA52LoHi_HandInterleavedChains, 4);
BENCHMARK_TEMPLATE(BM_IFMA52LoHi_SeqChains, 8);
BENCHMARK_TEMPLATE(BM_IFMA52LoHi_HandInterleavedChains, 8);

BENCHMARK_TEMPLATE(BM_Fp8x1MulRedc_SeqChains, 1);
BENCHMARK_TEMPLATE(BM_Fp8x1MulRedc_HandInterleavedChains, 1);
BENCHMARK_TEMPLATE(BM_Fp8x1MulRedc_SeqChains, 2);
BENCHMARK_TEMPLATE(BM_Fp8x1MulRedc_HandInterleavedChains, 2);
BENCHMARK_TEMPLATE(BM_Fp8x1MulRedc_SeqChains, 4);
BENCHMARK_TEMPLATE(BM_Fp8x1MulRedc_HandInterleavedChains, 4);
BENCHMARK_TEMPLATE(BM_Fp8x1MulRedc_SeqChains, 8);
BENCHMARK_TEMPLATE(BM_Fp8x1MulRedc_HandInterleavedChains, 8);

BENCHMARK_TEMPLATE(BM_Fp8x1RedcCore_SeqChains, 1);
BENCHMARK_TEMPLATE(BM_Fp8x1RedcCore_SeqChains, 2);
BENCHMARK_TEMPLATE(BM_Fp8x1RedcCore_HandInterleavedChains, 2);
BENCHMARK_TEMPLATE(BM_Fp8x1RedcCore_SeqChains, 4);
BENCHMARK_TEMPLATE(BM_Fp8x1RedcCore_HandInterleavedChains, 4);
BENCHMARK_TEMPLATE(BM_Fp8x1RedcCore_SeqChains, 8);
BENCHMARK_TEMPLATE(BM_Fp8x1RedcCore_HandInterleavedChains, 8);

} // namespace

#endif // BLST_AVX512
