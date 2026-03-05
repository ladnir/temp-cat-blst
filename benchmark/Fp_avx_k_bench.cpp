#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#if defined(BLST_AVX512)

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "blst/fields/Fp.h"
#include "blst/fields/Fp_avx.h"

using namespace blst;

#if defined(__linux__) && defined(__x86_64__) && defined(BLST_ASM) && defined(BLST_AVX512) && !defined(BLST_AVX512_PORTABLE)
extern "C" void blst_fp8x1_redc_core_k2(
    LFpPack8x1* out0,
    LFpPack8x1* out1,
    const UFpPack8x1* in0,
    const UFpPack8x1* in1);
#define BLST_HAVE_REDC_CORE_K2_ASM 1
#else
#define BLST_HAVE_REDC_CORE_K2_ASM 0
#endif

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
      for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
        delta[c][limb] = VSET1(gen() & BMASK);
      }
    }
  }

  std::array<UFpPack8x1, Chains> wide{};
  std::array<UFpPack8x1, Chains> delta{};
  std::array<LFpPack8x1, Chains> lazy{};
};

BLST_FORCEINLINE void redc_core_hand_interleaved_k2_once(
    LFpPack8x1& out0,
    LFpPack8x1& out1,
    const UFpPack8x1& in0,
    const UFpPack8x1& in1)
{
  cc512i a00 = in0[0], a01 = in0[1], a02 = in0[2], a03 = in0[3];
  cc512i a04 = in0[4], a05 = in0[5], a06 = in0[6], a07 = in0[7];
  cc512i a08 = in0[8], a09 = in0[9], a10 = in0[10], a11 = in0[11];
  cc512i a12 = in0[12], a13 = in0[13], a14 = in0[14], a15 = in0[15];

  cc512i b00 = in1[0], b01 = in1[1], b02 = in1[2], b03 = in1[3];
  cc512i b04 = in1[4], b05 = in1[5], b06 = in1[6], b07 = in1[7];
  cc512i b08 = in1[8], b09 = in1[9], b10 = in1[10], b11 = in1[11];
  cc512i b12 = in1[12], b13 = in1[13], b14 = in1[14], b15 = in1[15];

  cc512i y0 = VZERO(), y1 = VZERO(), y2 = VZERO(), y3 = VZERO();
  cc512i y4 = VZERO(), y5 = VZERO(), y6 = VZERO(), y7 = VZERO();
  cc512i z0 = VZERO(), z1 = VZERO(), z2 = VZERO(), z3 = VZERO();
  cc512i z4 = VZERO(), z5 = VZERO(), z6 = VZERO(), z7 = VZERO();

  const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
  const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
  const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
  const cc512i bmask = VSET1(BMASK), montw = VSET1(MONT_W_R48), zero = VZERO();
  cc512i u0, u1;

  u0 = VAND(VMACLO(zero, a00, montw), bmask);
  u1 = VAND(VMACLO(zero, b00, montw), bmask);
  a00 = VMACLO(a00, u0, p0); b00 = VMACLO(b00, u1, p0);
  a01 = VMACLO(a01, u0, p1); b01 = VMACLO(b01, u1, p1);
  a02 = VMACLO(a02, u0, p2); b02 = VMACLO(b02, u1, p2);
  a03 = VMACLO(a03, u0, p3); b03 = VMACLO(b03, u1, p3);
  a04 = VMACLO(a04, u0, p4); b04 = VMACLO(b04, u1, p4);
  a05 = VMACLO(a05, u0, p5); b05 = VMACLO(b05, u1, p5);
  a06 = VMACLO(a06, u0, p6); b06 = VMACLO(b06, u1, p6);
  a07 = VMACLO(a07, u0, p7); b07 = VMACLO(b07, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a01 = VADD(VADD(a01, VSRA(a00, BRADIX)), VSHL(y0, BALIGN));
  b01 = VADD(VADD(b01, VSRA(b00, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a01, montw), bmask);
  u1 = VAND(VMACLO(zero, b01, montw), bmask);
  a01 = VMACLO(a01, u0, p0); b01 = VMACLO(b01, u1, p0);
  a02 = VMACLO(a02, u0, p1); b02 = VMACLO(b02, u1, p1);
  a03 = VMACLO(a03, u0, p2); b03 = VMACLO(b03, u1, p2);
  a04 = VMACLO(a04, u0, p3); b04 = VMACLO(b04, u1, p3);
  a05 = VMACLO(a05, u0, p4); b05 = VMACLO(b05, u1, p4);
  a06 = VMACLO(a06, u0, p5); b06 = VMACLO(b06, u1, p5);
  a07 = VMACLO(a07, u0, p6); b07 = VMACLO(b07, u1, p6);
  a08 = VMACLO(a08, u0, p7); b08 = VMACLO(b08, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a02 = VADD(VADD(a02, VSRA(a01, BRADIX)), VSHL(y0, BALIGN));
  b02 = VADD(VADD(b02, VSRA(b01, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a02, montw), bmask);
  u1 = VAND(VMACLO(zero, b02, montw), bmask);
  a02 = VMACLO(a02, u0, p0); b02 = VMACLO(b02, u1, p0);
  a03 = VMACLO(a03, u0, p1); b03 = VMACLO(b03, u1, p1);
  a04 = VMACLO(a04, u0, p2); b04 = VMACLO(b04, u1, p2);
  a05 = VMACLO(a05, u0, p3); b05 = VMACLO(b05, u1, p3);
  a06 = VMACLO(a06, u0, p4); b06 = VMACLO(b06, u1, p4);
  a07 = VMACLO(a07, u0, p5); b07 = VMACLO(b07, u1, p5);
  a08 = VMACLO(a08, u0, p6); b08 = VMACLO(b08, u1, p6);
  a09 = VMACLO(a09, u0, p7); b09 = VMACLO(b09, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a03 = VADD(VADD(a03, VSRA(a02, BRADIX)), VSHL(y0, BALIGN));
  b03 = VADD(VADD(b03, VSRA(b02, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a03, montw), bmask);
  u1 = VAND(VMACLO(zero, b03, montw), bmask);
  a03 = VMACLO(a03, u0, p0); b03 = VMACLO(b03, u1, p0);
  a04 = VMACLO(a04, u0, p1); b04 = VMACLO(b04, u1, p1);
  a05 = VMACLO(a05, u0, p2); b05 = VMACLO(b05, u1, p2);
  a06 = VMACLO(a06, u0, p3); b06 = VMACLO(b06, u1, p3);
  a07 = VMACLO(a07, u0, p4); b07 = VMACLO(b07, u1, p4);
  a08 = VMACLO(a08, u0, p5); b08 = VMACLO(b08, u1, p5);
  a09 = VMACLO(a09, u0, p6); b09 = VMACLO(b09, u1, p6);
  a10 = VMACLO(a10, u0, p7); b10 = VMACLO(b10, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a04 = VADD(VADD(a04, VSRA(a03, BRADIX)), VSHL(y0, BALIGN));
  b04 = VADD(VADD(b04, VSRA(b03, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a04, montw), bmask);
  u1 = VAND(VMACLO(zero, b04, montw), bmask);
  a04 = VMACLO(a04, u0, p0); b04 = VMACLO(b04, u1, p0);
  a05 = VMACLO(a05, u0, p1); b05 = VMACLO(b05, u1, p1);
  a06 = VMACLO(a06, u0, p2); b06 = VMACLO(b06, u1, p2);
  a07 = VMACLO(a07, u0, p3); b07 = VMACLO(b07, u1, p3);
  a08 = VMACLO(a08, u0, p4); b08 = VMACLO(b08, u1, p4);
  a09 = VMACLO(a09, u0, p5); b09 = VMACLO(b09, u1, p5);
  a10 = VMACLO(a10, u0, p6); b10 = VMACLO(b10, u1, p6);
  a11 = VMACLO(a11, u0, p7); b11 = VMACLO(b11, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a05 = VADD(VADD(a05, VSRA(a04, BRADIX)), VSHL(y0, BALIGN));
  b05 = VADD(VADD(b05, VSRA(b04, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a05, montw), bmask);
  u1 = VAND(VMACLO(zero, b05, montw), bmask);
  a05 = VMACLO(a05, u0, p0); b05 = VMACLO(b05, u1, p0);
  a06 = VMACLO(a06, u0, p1); b06 = VMACLO(b06, u1, p1);
  a07 = VMACLO(a07, u0, p2); b07 = VMACLO(b07, u1, p2);
  a08 = VMACLO(a08, u0, p3); b08 = VMACLO(b08, u1, p3);
  a09 = VMACLO(a09, u0, p4); b09 = VMACLO(b09, u1, p4);
  a10 = VMACLO(a10, u0, p5); b10 = VMACLO(b10, u1, p5);
  a11 = VMACLO(a11, u0, p6); b11 = VMACLO(b11, u1, p6);
  a12 = VMACLO(a12, u0, p7); b12 = VMACLO(b12, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a06 = VADD(VADD(a06, VSRA(a05, BRADIX)), VSHL(y0, BALIGN));
  b06 = VADD(VADD(b06, VSRA(b05, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a06, montw), bmask);
  u1 = VAND(VMACLO(zero, b06, montw), bmask);
  a06 = VMACLO(a06, u0, p0); b06 = VMACLO(b06, u1, p0);
  a07 = VMACLO(a07, u0, p1); b07 = VMACLO(b07, u1, p1);
  a08 = VMACLO(a08, u0, p2); b08 = VMACLO(b08, u1, p2);
  a09 = VMACLO(a09, u0, p3); b09 = VMACLO(b09, u1, p3);
  a10 = VMACLO(a10, u0, p4); b10 = VMACLO(b10, u1, p4);
  a11 = VMACLO(a11, u0, p5); b11 = VMACLO(b11, u1, p5);
  a12 = VMACLO(a12, u0, p6); b12 = VMACLO(b12, u1, p6);
  a13 = VMACLO(a13, u0, p7); b13 = VMACLO(b13, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a07 = VADD(VADD(a07, VSRA(a06, BRADIX)), VSHL(y0, BALIGN));
  b07 = VADD(VADD(b07, VSRA(b06, BRADIX)), VSHL(z0, BALIGN));
  y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = y5; y5 = y6; y6 = y7; y7 = zero;
  z0 = z1; z1 = z2; z2 = z3; z3 = z4; z4 = z5; z5 = z6; z6 = z7; z7 = zero;

  u0 = VAND(VMACLO(zero, a07, montw), bmask);
  u1 = VAND(VMACLO(zero, b07, montw), bmask);
  a07 = VMACLO(a07, u0, p0); b07 = VMACLO(b07, u1, p0);
  a08 = VMACLO(a08, u0, p1); b08 = VMACLO(b08, u1, p1);
  a09 = VMACLO(a09, u0, p2); b09 = VMACLO(b09, u1, p2);
  a10 = VMACLO(a10, u0, p3); b10 = VMACLO(b10, u1, p3);
  a11 = VMACLO(a11, u0, p4); b11 = VMACLO(b11, u1, p4);
  a12 = VMACLO(a12, u0, p5); b12 = VMACLO(b12, u1, p5);
  a13 = VMACLO(a13, u0, p6); b13 = VMACLO(b13, u1, p6);
  a14 = VMACLO(a14, u0, p7); b14 = VMACLO(b14, u1, p7);
  y0 = VMACHI(y0, u0, p0); z0 = VMACHI(z0, u1, p0);
  y1 = VMACHI(y1, u0, p1); z1 = VMACHI(z1, u1, p1);
  y2 = VMACHI(y2, u0, p2); z2 = VMACHI(z2, u1, p2);
  y3 = VMACHI(y3, u0, p3); z3 = VMACHI(z3, u1, p3);
  y4 = VMACHI(y4, u0, p4); z4 = VMACHI(z4, u1, p4);
  y5 = VMACHI(y5, u0, p5); z5 = VMACHI(z5, u1, p5);
  y6 = VMACHI(y6, u0, p6); z6 = VMACHI(z6, u1, p6);
  y7 = VMACHI(y7, u0, p7); z7 = VMACHI(z7, u1, p7);
  a08 = VADD(VADD(a08, VSRA(a07, BRADIX)), VSHL(y0, BALIGN));
  b08 = VADD(VADD(b08, VSRA(b07, BRADIX)), VSHL(z0, BALIGN));

  a09 = VADD(a09, VSHL(y1, BALIGN));
  a10 = VADD(a10, VSHL(y2, BALIGN));
  a11 = VADD(a11, VSHL(y3, BALIGN));
  a12 = VADD(a12, VSHL(y4, BALIGN));
  a13 = VADD(a13, VSHL(y5, BALIGN));
  a14 = VADD(a14, VSHL(y6, BALIGN));
  a15 = VADD(a15, VSHL(y7, BALIGN));

  b09 = VADD(b09, VSHL(z1, BALIGN));
  b10 = VADD(b10, VSHL(z2, BALIGN));
  b11 = VADD(b11, VSHL(z3, BALIGN));
  b12 = VADD(b12, VSHL(z4, BALIGN));
  b13 = VADD(b13, VSHL(z5, BALIGN));
  b14 = VADD(b14, VSHL(z6, BALIGN));
  b15 = VADD(b15, VSHL(z7, BALIGN));

  cc512i r00 = a08, r01 = a09, r02 = a10, r03 = a11;
  cc512i r04 = a12, r05 = a13, r06 = a14, r07 = a15;
  cc512i r10 = b08, r11 = b09, r12 = b10, r13 = b11;
  cc512i r14 = b12, r15 = b13, r16 = b14, r17 = b15;

  r01 = VADD(r01, VSRA(r00, BRADIX)); r00 = VAND(r00, bmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
  r02 = VADD(r02, VSRA(r01, BRADIX)); r01 = VAND(r01, bmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
  r03 = VADD(r03, VSRA(r02, BRADIX)); r02 = VAND(r02, bmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
  r04 = VADD(r04, VSRA(r03, BRADIX)); r03 = VAND(r03, bmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
  r05 = VADD(r05, VSRA(r04, BRADIX)); r04 = VAND(r04, bmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
  r06 = VADD(r06, VSRA(r05, BRADIX)); r05 = VAND(r05, bmask);
  r16 = VADD(r16, VSRA(r15, BRADIX)); r15 = VAND(r15, bmask);
  r07 = VADD(r07, VSRA(r06, BRADIX)); r06 = VAND(r06, bmask); r07 = VAND(r07, bmask);
  r17 = VADD(r17, VSRA(r16, BRADIX)); r16 = VAND(r16, bmask); r17 = VAND(r17, bmask);

  out0[0] = r00; out0[1] = r01; out0[2] = r02; out0[3] = r03;
  out0[4] = r04; out0[5] = r05; out0[6] = r06; out0[7] = r07;
  out1[0] = r10; out1[1] = r11; out1[2] = r12; out1[3] = r13;
  out1[4] = r14; out1[5] = r15; out1[6] = r16; out1[7] = r17;
}

template <typename Fn>
void redc_core_k2_selfcheck_impl(Fn&& fn, const char* tag)
{
#ifndef NDEBUG
  std::mt19937_64 gen(0xC001D00DULL);
  constexpr std::size_t kCases = 64;
  for (std::size_t t = 0; t < kCases; ++t) {
    std::array<Fp, 8> lanes_a0{};
    std::array<Fp, 8> lanes_b0{};
    std::array<Fp, 8> lanes_a1{};
    std::array<Fp, 8> lanes_b1{};
    for (std::size_t i = 0; i < 8; ++i) {
      random_redc_fp(lanes_a0[i], gen);
      random_redc_fp(lanes_b0[i], gen);
      random_redc_fp(lanes_a1[i], gen);
      random_redc_fp(lanes_b1[i], gen);
    }

    FpPack8x1 a0{}, b0{}, a1{}, b1{};
    FpPack8x1::from_aos(a0, lanes_a0.data());
    FpPack8x1::from_aos(b0, lanes_b0.data());
    FpPack8x1::from_aos(a1, lanes_a1.data());
    FpPack8x1::from_aos(b1, lanes_b1.data());

    UFpPack8x1 wide0{};
    UFpPack8x1 wide1{};
    FpPack8x1::mul_v2(wide0, a0, b0);
    FpPack8x1::mul_v2(wide1, a1, b1);

    LFpPack8x1 got0{};
    LFpPack8x1 got1{};
    fn(got0, got1, wide0, wide1);

    LFpPack8x1 exp0{};
    LFpPack8x1 exp1{};
    UFpPack8x1::redc_core(exp0, wide0);
    UFpPack8x1::redc_core(exp1, wide1);

    FpPack8x1 got0_norm{}, got1_norm{}, exp0_norm{}, exp1_norm{};
    UFpPack8x1::normalize(got0_norm, got0);
    UFpPack8x1::normalize(got1_norm, got1);
    UFpPack8x1::normalize(exp0_norm, exp0);
    UFpPack8x1::normalize(exp1_norm, exp1);

    std::array<Fp, 8> got0_lanes{}, got1_lanes{}, exp0_lanes{}, exp1_lanes{};
    FpPack8x1::to_aos(got0_lanes.data(), got0_norm);
    FpPack8x1::to_aos(got1_lanes.data(), got1_norm);
    FpPack8x1::to_aos(exp0_lanes.data(), exp0_norm);
    FpPack8x1::to_aos(exp1_lanes.data(), exp1_norm);

    for (std::size_t lane = 0; lane < 8; ++lane) {
      for (std::size_t limb = 0; limb < Fp::num_limb(); ++limb) {
        if (got0_lanes[lane][limb] != exp0_lanes[lane][limb]) {
          std::fprintf(
              stderr,
              "%s selfcheck failed (case=%zu stream=0 lane=%zu limb=%zu got=%llx exp=%llx)\n",
              tag,
              t,
              lane,
              limb,
              static_cast<unsigned long long>(got0_lanes[lane][limb]),
              static_cast<unsigned long long>(exp0_lanes[lane][limb]));
          std::abort();
        }
        if (got1_lanes[lane][limb] != exp1_lanes[lane][limb]) {
          std::fprintf(
              stderr,
              "%s selfcheck failed (case=%zu stream=1 lane=%zu limb=%zu got=%llx exp=%llx)\n",
              tag,
              t,
              lane,
              limb,
              static_cast<unsigned long long>(got1_lanes[lane][limb]),
              static_cast<unsigned long long>(exp1_lanes[lane][limb]));
          std::abort();
        }
      }
    }
  }
#else
  (void)fn;
  (void)tag;
#endif
}

void redc_core_hand_k2_selfcheck()
{
  redc_core_k2_selfcheck_impl(redc_core_hand_interleaved_k2_once, "redc_core_hand_k2");
}

#if BLST_HAVE_REDC_CORE_K2_ASM
BLST_FORCEINLINE void redc_core_asm_k2_once(
    LFpPack8x1& out0,
    LFpPack8x1& out1,
    const UFpPack8x1& in0,
    const UFpPack8x1& in1)
{
  blst_fp8x1_redc_core_k2(&out0, &out1, &in0, &in1);
}

void redc_core_asm_k2_selfcheck()
{
  redc_core_k2_selfcheck_impl(redc_core_asm_k2_once, "redc_core_asm_k2");
}
#endif

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
  static bool checked = false;
  if (!checked) {
    redc_core_hand_k2_selfcheck();
    checked = true;
  }

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

  for (auto _ : state) {
    for (int c = 0; c < Chains; ++c) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        UFpPack8x1::redc_core(data.lazy[c], data.wide[c]);
        for (std::size_t limb = 0; limb < NWORDS; ++limb) {
          data.wide[c][limb] =
              VADD(data.wide[c][limb], VADD(data.lazy[c][limb], data.delta[c][limb]));
          data.wide[c][limb + NWORDS] = VADD(
              data.wide[c][limb + NWORDS],
              VADD(VSHL(data.lazy[c][limb], BALIGN), data.delta[c][limb + NWORDS]));
        }
      }
    }
  }

  cc512i sink = VZERO();
  for (int c = 0; c < Chains; ++c) {
    for (std::size_t limb = 0; limb < NWORDS; ++limb) {
      sink = VADD(sink, data.lazy[c][limb]);
    }
  }
  benchmark::DoNotOptimize(sink);
  benchmark::ClobberMemory();

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
  static bool checked = false;
  if (!checked) {
    redc_core_hand_k2_selfcheck();
    checked = true;
  }

  for (auto _ : state) {
    for (std::size_t i = 0; i < kUnroll; ++i) {
      redc_core_hand_interleaved_k2_once(data.lazy[0], data.lazy[1], data.wide[0], data.wide[1]);
      for (std::size_t limb = 0; limb < NWORDS; ++limb) {
        data.wide[0][limb] =
            VADD(data.wide[0][limb], VADD(data.lazy[0][limb], data.delta[0][limb]));
        data.wide[0][limb + NWORDS] =
            VADD(data.wide[0][limb + NWORDS],
                 VADD(VSHL(data.lazy[0][limb], BALIGN), data.delta[0][limb + NWORDS]));
        data.wide[1][limb] =
            VADD(data.wide[1][limb], VADD(data.lazy[1][limb], data.delta[1][limb]));
        data.wide[1][limb + NWORDS] =
            VADD(data.wide[1][limb + NWORDS],
                 VADD(VSHL(data.lazy[1][limb], BALIGN), data.delta[1][limb + NWORDS]));
      }
    }
  }

  cc512i sink = VZERO();
  for (std::size_t limb = 0; limb < NWORDS; ++limb) {
    sink = VADD(sink, data.lazy[0][limb]);
    sink = VADD(sink, data.lazy[1][limb]);
  }
  benchmark::DoNotOptimize(sink);
  benchmark::ClobberMemory();

  const std::size_t ops_per_iter = 2 * kUnroll;
  state.SetItemsProcessed(state.iterations() * ops_per_iter);
  state.counters["redc_core_per_iter"] = static_cast<double>(ops_per_iter);
  state.counters["chains"] = 2.0;
  state.counters["unroll"] = static_cast<double>(kUnroll);
}

#if BLST_HAVE_REDC_CORE_K2_ASM
void BM_Fp8x1RedcCore_AsmK2(benchmark::State& state)
{
  constexpr std::size_t kUnroll = 8;
  Fp8x1RedcCoreStreamData<2> data{};
  static bool checked = false;
  if (!checked) {
    redc_core_asm_k2_selfcheck();
    checked = true;
  }

  for (auto _ : state) {
    for (std::size_t i = 0; i < kUnroll; ++i) {
      redc_core_asm_k2_once(data.lazy[0], data.lazy[1], data.wide[0], data.wide[1]);
      for (std::size_t limb = 0; limb < NWORDS; ++limb) {
        data.wide[0][limb] =
            VADD(data.wide[0][limb], VADD(data.lazy[0][limb], data.delta[0][limb]));
        data.wide[0][limb + NWORDS] =
            VADD(data.wide[0][limb + NWORDS],
                 VADD(VSHL(data.lazy[0][limb], BALIGN), data.delta[0][limb + NWORDS]));
        data.wide[1][limb] =
            VADD(data.wide[1][limb], VADD(data.lazy[1][limb], data.delta[1][limb]));
        data.wide[1][limb + NWORDS] =
            VADD(data.wide[1][limb + NWORDS],
                 VADD(VSHL(data.lazy[1][limb], BALIGN), data.delta[1][limb + NWORDS]));
      }
    }
  }

  cc512i sink = VZERO();
  for (std::size_t limb = 0; limb < NWORDS; ++limb) {
    sink = VADD(sink, data.lazy[0][limb]);
    sink = VADD(sink, data.lazy[1][limb]);
  }
  benchmark::DoNotOptimize(sink);
  benchmark::ClobberMemory();

  const std::size_t ops_per_iter = 2 * kUnroll;
  state.SetItemsProcessed(state.iterations() * ops_per_iter);
  state.counters["redc_core_per_iter"] = static_cast<double>(ops_per_iter);
  state.counters["chains"] = 2.0;
  state.counters["unroll"] = static_cast<double>(kUnroll);
}

#endif

template <>
void BM_Fp8x1RedcCore_HandInterleavedChains<4>(benchmark::State& state)
{
  constexpr std::size_t kUnroll = 8;
  Fp8x1RedcCoreStreamData<4> data{};

  for (auto _ : state) {
    for (std::size_t i = 0; i < kUnroll; ++i) {
      redc_core_hand_interleaved_k2_once(data.lazy[0], data.lazy[1], data.wide[0], data.wide[1]);
      redc_core_hand_interleaved_k2_once(data.lazy[2], data.lazy[3], data.wide[2], data.wide[3]);
      for (std::size_t limb = 0; limb < NWORDS; ++limb) {
        data.wide[0][limb] =
            VADD(data.wide[0][limb], VADD(data.lazy[0][limb], data.delta[0][limb]));
        data.wide[0][limb + NWORDS] =
            VADD(data.wide[0][limb + NWORDS],
                 VADD(VSHL(data.lazy[0][limb], BALIGN), data.delta[0][limb + NWORDS]));
        data.wide[1][limb] =
            VADD(data.wide[1][limb], VADD(data.lazy[1][limb], data.delta[1][limb]));
        data.wide[1][limb + NWORDS] =
            VADD(data.wide[1][limb + NWORDS],
                 VADD(VSHL(data.lazy[1][limb], BALIGN), data.delta[1][limb + NWORDS]));
        data.wide[2][limb] =
            VADD(data.wide[2][limb], VADD(data.lazy[2][limb], data.delta[2][limb]));
        data.wide[2][limb + NWORDS] =
            VADD(data.wide[2][limb + NWORDS],
                 VADD(VSHL(data.lazy[2][limb], BALIGN), data.delta[2][limb + NWORDS]));
        data.wide[3][limb] =
            VADD(data.wide[3][limb], VADD(data.lazy[3][limb], data.delta[3][limb]));
        data.wide[3][limb + NWORDS] =
            VADD(data.wide[3][limb + NWORDS],
                 VADD(VSHL(data.lazy[3][limb], BALIGN), data.delta[3][limb + NWORDS]));
      }
    }
  }

  cc512i sink = VZERO();
  for (std::size_t limb = 0; limb < NWORDS; ++limb) {
    sink = VADD(sink, data.lazy[0][limb]);
    sink = VADD(sink, data.lazy[1][limb]);
    sink = VADD(sink, data.lazy[2][limb]);
    sink = VADD(sink, data.lazy[3][limb]);
  }
  benchmark::DoNotOptimize(sink);
  benchmark::ClobberMemory();

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

  for (auto _ : state) {
    for (std::size_t i = 0; i < kUnroll; ++i) {
      redc_core_hand_interleaved_k2_once(data.lazy[0], data.lazy[1], data.wide[0], data.wide[1]);
      redc_core_hand_interleaved_k2_once(data.lazy[2], data.lazy[3], data.wide[2], data.wide[3]);
      redc_core_hand_interleaved_k2_once(data.lazy[4], data.lazy[5], data.wide[4], data.wide[5]);
      redc_core_hand_interleaved_k2_once(data.lazy[6], data.lazy[7], data.wide[6], data.wide[7]);
      for (std::size_t limb = 0; limb < NWORDS; ++limb) {
        data.wide[0][limb] =
            VADD(data.wide[0][limb], VADD(data.lazy[0][limb], data.delta[0][limb]));
        data.wide[0][limb + NWORDS] =
            VADD(data.wide[0][limb + NWORDS],
                 VADD(VSHL(data.lazy[0][limb], BALIGN), data.delta[0][limb + NWORDS]));
        data.wide[1][limb] =
            VADD(data.wide[1][limb], VADD(data.lazy[1][limb], data.delta[1][limb]));
        data.wide[1][limb + NWORDS] =
            VADD(data.wide[1][limb + NWORDS],
                 VADD(VSHL(data.lazy[1][limb], BALIGN), data.delta[1][limb + NWORDS]));
        data.wide[2][limb] =
            VADD(data.wide[2][limb], VADD(data.lazy[2][limb], data.delta[2][limb]));
        data.wide[2][limb + NWORDS] =
            VADD(data.wide[2][limb + NWORDS],
                 VADD(VSHL(data.lazy[2][limb], BALIGN), data.delta[2][limb + NWORDS]));
        data.wide[3][limb] =
            VADD(data.wide[3][limb], VADD(data.lazy[3][limb], data.delta[3][limb]));
        data.wide[3][limb + NWORDS] =
            VADD(data.wide[3][limb + NWORDS],
                 VADD(VSHL(data.lazy[3][limb], BALIGN), data.delta[3][limb + NWORDS]));
        data.wide[4][limb] =
            VADD(data.wide[4][limb], VADD(data.lazy[4][limb], data.delta[4][limb]));
        data.wide[4][limb + NWORDS] =
            VADD(data.wide[4][limb + NWORDS],
                 VADD(VSHL(data.lazy[4][limb], BALIGN), data.delta[4][limb + NWORDS]));
        data.wide[5][limb] =
            VADD(data.wide[5][limb], VADD(data.lazy[5][limb], data.delta[5][limb]));
        data.wide[5][limb + NWORDS] =
            VADD(data.wide[5][limb + NWORDS],
                 VADD(VSHL(data.lazy[5][limb], BALIGN), data.delta[5][limb + NWORDS]));
        data.wide[6][limb] =
            VADD(data.wide[6][limb], VADD(data.lazy[6][limb], data.delta[6][limb]));
        data.wide[6][limb + NWORDS] =
            VADD(data.wide[6][limb + NWORDS],
                 VADD(VSHL(data.lazy[6][limb], BALIGN), data.delta[6][limb + NWORDS]));
        data.wide[7][limb] =
            VADD(data.wide[7][limb], VADD(data.lazy[7][limb], data.delta[7][limb]));
        data.wide[7][limb + NWORDS] =
            VADD(data.wide[7][limb + NWORDS],
                 VADD(VSHL(data.lazy[7][limb], BALIGN), data.delta[7][limb + NWORDS]));
      }
    }
  }

  cc512i sink = VZERO();
  for (std::size_t limb = 0; limb < NWORDS; ++limb) {
    sink = VADD(sink, data.lazy[0][limb]);
    sink = VADD(sink, data.lazy[1][limb]);
    sink = VADD(sink, data.lazy[2][limb]);
    sink = VADD(sink, data.lazy[3][limb]);
    sink = VADD(sink, data.lazy[4][limb]);
    sink = VADD(sink, data.lazy[5][limb]);
    sink = VADD(sink, data.lazy[6][limb]);
    sink = VADD(sink, data.lazy[7][limb]);
  }
  benchmark::DoNotOptimize(sink);
  benchmark::ClobberMemory();

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

#if BLST_HAVE_REDC_CORE_K2_ASM
BENCHMARK(BM_Fp8x1RedcCore_AsmK2);
#endif

} // namespace

#endif // BLST_AVX512
