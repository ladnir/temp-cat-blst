#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#if defined(BLST_AVX512)

#include <array>
#include <cstddef>
#include <mutex>
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

struct Fp8x1FusedBaselineData
{
  explicit Fp8x1FusedBaselineData(const std::size_t n_packs)
      : n_packs_(n_packs)
      , a(n_packs_)
      , b(n_packs_)
      , lazy(n_packs_)
      , out(n_packs_)
  {
    std::mt19937_64 gen(0xB17A51E0ULL + static_cast<u64>(n_packs_));
    for (std::size_t i = 0; i < n_packs_; ++i) {
      std::array<Fp, 8> lane_a{};
      std::array<Fp, 8> lane_b{};
      for (std::size_t lane = 0; lane < 8; ++lane) {
        random_redc_fp(lane_a[lane], gen);
        random_redc_fp(lane_b[lane], gen);
      }
      FpPack8x1::from_aos(a[i], lane_a.data());
      FpPack8x1::from_aos(b[i], lane_b.data());
    }
  }

  const std::size_t n_packs_;
  std::vector<FpPack8x1> a;
  std::vector<FpPack8x1> b;
  std::vector<LFpPack8x1> lazy;
  std::vector<FpPack8x1> out;
};

struct Fp8x1MontConsts
{
  cc512i p0;
  cc512i p1;
  cc512i p2;
  cc512i p3;
  cc512i p4;
  cc512i p5;
  cc512i p6;
  cc512i p7;
  cc512i bmask;
  cc512i montw;
  cc512i zero;
};

bool fp8x1_pack_equal(const FpPack8x1& lhs, const FpPack8x1& rhs)
{
  std::array<Fp, 8> l{};
  std::array<Fp, 8> r{};
  FpPack8x1::to_aos(l.data(), lhs);
  FpPack8x1::to_aos(r.data(), rhs);
  for (std::size_t lane = 0; lane < 8; ++lane) {
    if (!Fp::equal(l[lane], r[lane])) {
      return false;
    }
  }
  return true;
}

BLST_FORCEINLINE void fp8x1_mul_redc_core_k1_fused_baseline(
    LFpPack8x1& r,
    const FpPack8x1& a,
    const FpPack8x1& b,
    const Fp8x1MontConsts& c)
{
  BLST_ASSERT_NO_OVERLAP_AVX(r, a);
  BLST_ASSERT_NO_OVERLAP_AVX(r, b);

  const cc512i ia0 = a[0], ia1 = a[1], ia2 = a[2], ia3 = a[3];
  const cc512i ia4 = a[4], ia5 = a[5], ia6 = a[6], ia7 = a[7];
  const cc512i ib0 = b[0], ib1 = b[1], ib2 = b[2], ib3 = b[3];
  const cc512i ib4 = b[4], ib5 = b[5], ib6 = b[6], ib7 = b[7];
  cc512i ta0, ta1, ta2, ta3;
  cc512i tb0, tb1, tb2, tb3;
  cc512i z0  = VZERO(), z1  = VZERO(), z2  = VZERO(), z3  = VZERO();
  cc512i z4  = VZERO(), z5  = VZERO(), z6  = VZERO(), z7  = VZERO();
  cc512i z8  = VZERO(), z9  = VZERO(), z10 = VZERO(), z11 = VZERO();
  cc512i z12 = VZERO(), z13 = VZERO(), z14 = VZERO(), z15 = VZERO();
  cc512i y0  = VZERO(), y1  = VZERO(), y2  = VZERO(), y3  = VZERO();
  cc512i y4  = VZERO(), y5  = VZERO(), y6  = VZERO(), y7  = VZERO();
  cc512i y8  = VZERO(), y9  = VZERO(), y10 = VZERO(), y11 = VZERO();
  cc512i y12 = VZERO(), y13 = VZERO(), y14 = VZERO();
  cc512i m0  = VZERO(), m1  = VZERO(), m2  = VZERO(), m3  = VZERO();
  cc512i m4  = VZERO(), m5  = VZERO(), m6  = VZERO(), m7  = VZERO();

  z0 = VMACLO(z0, ia0, ib0);
  y0 = VMACHI(y0, ia0, ib0);
  y0 = VSHL(y0, BALIGN);

  z1 = VMACLO(y0, ia0, ib1); z1 = VMACLO(z1, ia1, ib0);
  y1 = VMACHI(y1, ia0, ib1); y1 = VMACHI(y1, ia1, ib0);
  y1 = VSHL(y1, BALIGN);

  z2 = VMACLO(y1, ia0, ib2); z2 = VMACLO(z2, ia1, ib1); z2 = VMACLO(z2, ia2, ib0);
  y2 = VMACHI(y2, ia0, ib2); y2 = VMACHI(y2, ia1, ib1); y2 = VMACHI(y2, ia2, ib0);
  y2 = VSHL(y2, BALIGN);

  z3 = VMACLO(y2, ia0, ib3); z3 = VMACLO(z3, ia1, ib2); z3 = VMACLO(z3, ia2, ib1);
  z3 = VMACLO(z3, ia3, ib0);
  y3 = VMACHI(y3, ia0, ib3); y3 = VMACHI(y3, ia1, ib2); y3 = VMACHI(y3, ia2, ib1);
  y3 = VMACHI(y3, ia3, ib0);
  y3 = VSHL(y3, BALIGN);

  z4 = VMACLO(y3, ia1, ib3); z4 = VMACLO(z4, ia2, ib2); z4 = VMACLO(z4, ia3, ib1);
  y4 = VMACHI(y4, ia1, ib3); y4 = VMACHI(y4, ia2, ib2); y4 = VMACHI(y4, ia3, ib1);
  y4 = VSHL(y4, BALIGN);

  z5 = VMACLO(y4, ia2, ib3); z5 = VMACLO(z5, ia3, ib2);
  y5 = VMACHI(y5, ia2, ib3); y5 = VMACHI(y5, ia3, ib2);
  y5 = VSHL(y5, BALIGN);

  z6 = VMACLO(y5, ia3, ib3);
  y6 = VMACHI(y6, ia3, ib3);
  y6 = VSHL(y6, BALIGN);
  z7 = y6;

  z8 = VMACLO(z8, ia4, ib4);
  y8 = VMACHI(y8, ia4, ib4);
  y8 = VSHL(y8, BALIGN);

  z9 = VMACLO(y8, ia4, ib5); z9 = VMACLO(z9, ia5, ib4);
  y9 = VMACHI(y9, ia4, ib5); y9 = VMACHI(y9, ia5, ib4);
  y9 = VSHL(y9, BALIGN);

  z10 = VMACLO(y9, ia4, ib6); z10 = VMACLO(z10, ia5, ib5); z10 = VMACLO(z10, ia6, ib4);
  y10 = VMACHI(y10, ia4, ib6); y10 = VMACHI(y10, ia5, ib5); y10 = VMACHI(y10, ia6, ib4);
  y10 = VSHL(y10, BALIGN);

  z11 = VMACLO(y10, ia4, ib7); z11 = VMACLO(z11, ia5, ib6);
  z11 = VMACLO(z11, ia6, ib5); z11 = VMACLO(z11, ia7, ib4);
  y11 = VMACHI(y11, ia4, ib7); y11 = VMACHI(y11, ia5, ib6);
  y11 = VMACHI(y11, ia6, ib5); y11 = VMACHI(y11, ia7, ib4);
  y11 = VSHL(y11, BALIGN);

  z12 = VMACLO(y11, ia5, ib7); z12 = VMACLO(z12, ia6, ib6); z12 = VMACLO(z12, ia7, ib5);
  y12 = VMACHI(y12, ia5, ib7); y12 = VMACHI(y12, ia6, ib6); y12 = VMACHI(y12, ia7, ib5);
  y12 = VSHL(y12, BALIGN);

  z13 = VMACLO(y12, ia6, ib7); z13 = VMACLO(z13, ia7, ib6);
  y13 = VMACHI(y13, ia6, ib7); y13 = VMACHI(y13, ia7, ib6);
  y13 = VSHL(y13, BALIGN);

  z14 = VMACLO(y13, ia7, ib7);
  y14 = VMACHI(y14, ia7, ib7);
  y14 = VSHL(y14, BALIGN);
  z15 = y14;

  ta0 = VADD(ia0, ia4); ta1 = VADD(ia1, ia5);
  ta2 = VADD(ia2, ia6); ta3 = VADD(ia3, ia7);
  tb0 = VADD(ib0, ib4); tb1 = VADD(ib1, ib5);
  tb2 = VADD(ib2, ib6); tb3 = VADD(ib3, ib7);

  y0 = y1 = y2 = y3 = y4 = y5 = y6 = y7 = VZERO();

  m0 = VMACLO(m0, ta0, tb0);
  y0 = VMACHI(y0, ta0, tb0);
  y0 = VSHL(y0, BALIGN);

  m1 = VMACLO(y0, ta0, tb1); m1 = VMACLO(m1, ta1, tb0);
  y1 = VMACHI(y1, ta0, tb1); y1 = VMACHI(y1, ta1, tb0);
  y1 = VSHL(y1, BALIGN);

  m2 = VMACLO(y1, ta0, tb2); m2 = VMACLO(m2, ta1, tb1); m2 = VMACLO(m2, ta2, tb0);
  y2 = VMACHI(y2, ta0, tb2); y2 = VMACHI(y2, ta1, tb1); y2 = VMACHI(y2, ta2, tb0);
  y2 = VSHL(y2, BALIGN);

  m3 = VMACLO(y2, ta0, tb3); m3 = VMACLO(m3, ta1, tb2);
  m3 = VMACLO(m3, ta2, tb1); m3 = VMACLO(m3, ta3, tb0);
  y3 = VMACHI(y3, ta0, tb3); y3 = VMACHI(y3, ta1, tb2);
  y3 = VMACHI(y3, ta2, tb1); y3 = VMACHI(y3, ta3, tb0);
  y3 = VSHL(y3, BALIGN);

  m4 = VMACLO(y3, ta1, tb3); m4 = VMACLO(m4, ta2, tb2); m4 = VMACLO(m4, ta3, tb1);
  y4 = VMACHI(y4, ta1, tb3); y4 = VMACHI(y4, ta2, tb2); y4 = VMACHI(y4, ta3, tb1);
  y4 = VSHL(y4, BALIGN);

  m5 = VMACLO(y4, ta2, tb3); m5 = VMACLO(m5, ta3, tb2);
  y5 = VMACHI(y5, ta2, tb3); y5 = VMACHI(y5, ta3, tb2);
  y5 = VSHL(y5, BALIGN);

  m6 = VMACLO(y5, ta3, tb3);
  y6 = VMACHI(y6, ta3, tb3);
  y6 = VSHL(y6, BALIGN);
  m7 = y6;

  m0 = VSUB(m0, VADD(z0, z8));   m1 = VSUB(m1, VADD(z1, z9));
  m2 = VSUB(m2, VADD(z2, z10));  m3 = VSUB(m3, VADD(z3, z11));
  m4 = VSUB(m4, VADD(z4, z12));  m5 = VSUB(m5, VADD(z5, z13));
  m6 = VSUB(m6, VADD(z6, z14));  m7 = VSUB(m7, VADD(z7, z15));

  z4  = VADD(z4, m0);   z5  = VADD(z5, m1);
  z6  = VADD(z6, m2);   z7  = VADD(z7, m3);
  z8  = VADD(z8, m4);   z9  = VADD(z9, m5);
  z10 = VADD(z10, m6);  z11 = VADD(z11, m7);

  y0 = y1 = y2 = y3 = y4 = y5 = y6 = y7 = c.zero;
  y8 = y9 = y10 = y11 = y12 = y13 = y14 = c.zero;
  cc512i u;

  u  = VAND(VMACLO(c.zero, z0, c.montw), c.bmask);
  z0 = VMACLO(z0, u, c.p0); z1 = VMACLO(z1, u, c.p1); z2 = VMACLO(z2, u, c.p2);
  z3 = VMACLO(z3, u, c.p3); z4 = VMACLO(z4, u, c.p4); z5 = VMACLO(z5, u, c.p5);
  z6 = VMACLO(z6, u, c.p6); z7 = VMACLO(z7, u, c.p7);
  y0 = VMACHI(y0, u, c.p0); y1 = VMACHI(y1, u, c.p1); y2 = VMACHI(y2, u, c.p2);
  y3 = VMACHI(y3, u, c.p3); y4 = VMACHI(y4, u, c.p4); y5 = VMACHI(y5, u, c.p5);
  y6 = VMACHI(y6, u, c.p6); y7 = VMACHI(y7, u, c.p7);
  z1 = VADD(VADD(z1, VSRA(z0, BRADIX)), VSHL(y0, BALIGN));

  u  = VAND(VMACLO(c.zero, z1, c.montw), c.bmask);
  z1 = VMACLO(z1, u, c.p0); z2 = VMACLO(z2, u, c.p1); z3 = VMACLO(z3, u, c.p2);
  z4 = VMACLO(z4, u, c.p3); z5 = VMACLO(z5, u, c.p4); z6 = VMACLO(z6, u, c.p5);
  z7 = VMACLO(z7, u, c.p6); z8 = VMACLO(z8, u, c.p7);
  y1 = VMACHI(y1, u, c.p0); y2 = VMACHI(y2, u, c.p1); y3 = VMACHI(y3, u, c.p2);
  y4 = VMACHI(y4, u, c.p3); y5 = VMACHI(y5, u, c.p4); y6 = VMACHI(y6, u, c.p5);
  y7 = VMACHI(y7, u, c.p6); y8 = VMACHI(y8, u, c.p7);
  z2 = VADD(VADD(z2, VSRA(z1, BRADIX)), VSHL(y1, BALIGN));

  u  = VAND(VMACLO(c.zero, z2, c.montw), c.bmask);
  z2 = VMACLO(z2, u, c.p0); z3 = VMACLO(z3, u, c.p1); z4 = VMACLO(z4, u, c.p2);
  z5 = VMACLO(z5, u, c.p3); z6 = VMACLO(z6, u, c.p4); z7 = VMACLO(z7, u, c.p5);
  z8 = VMACLO(z8, u, c.p6); z9 = VMACLO(z9, u, c.p7);
  y2 = VMACHI(y2, u, c.p0); y3 = VMACHI(y3, u, c.p1); y4 = VMACHI(y4, u, c.p2);
  y5 = VMACHI(y5, u, c.p3); y6 = VMACHI(y6, u, c.p4); y7 = VMACHI(y7, u, c.p5);
  y8 = VMACHI(y8, u, c.p6); y9 = VMACHI(y9, u, c.p7);
  z3 = VADD(VADD(z3, VSRA(z2, BRADIX)), VSHL(y2, BALIGN));

  u  = VAND(VMACLO(c.zero, z3, c.montw), c.bmask);
  z3 = VMACLO(z3, u, c.p0); z4 = VMACLO(z4, u, c.p1); z5 = VMACLO(z5, u, c.p2);
  z6 = VMACLO(z6, u, c.p3); z7 = VMACLO(z7, u, c.p4); z8 = VMACLO(z8, u, c.p5);
  z9 = VMACLO(z9, u, c.p6); z10 = VMACLO(z10, u, c.p7);
  y3 = VMACHI(y3, u, c.p0); y4 = VMACHI(y4, u, c.p1); y5 = VMACHI(y5, u, c.p2);
  y6 = VMACHI(y6, u, c.p3); y7 = VMACHI(y7, u, c.p4); y8 = VMACHI(y8, u, c.p5);
  y9 = VMACHI(y9, u, c.p6); y10 = VMACHI(y10, u, c.p7);
  z4 = VADD(VADD(z4, VSRA(z3, BRADIX)), VSHL(y3, BALIGN));

  u   = VAND(VMACLO(c.zero, z4, c.montw), c.bmask);
  z4  = VMACLO(z4, u, c.p0); z5  = VMACLO(z5, u, c.p1); z6 = VMACLO(z6, u, c.p2);
  z7  = VMACLO(z7, u, c.p3); z8  = VMACLO(z8, u, c.p4); z9 = VMACLO(z9, u, c.p5);
  z10 = VMACLO(z10, u, c.p6); z11 = VMACLO(z11, u, c.p7);
  y4  = VMACHI(y4, u, c.p0); y5  = VMACHI(y5, u, c.p1); y6 = VMACHI(y6, u, c.p2);
  y7  = VMACHI(y7, u, c.p3); y8  = VMACHI(y8, u, c.p4); y9 = VMACHI(y9, u, c.p5);
  y10 = VMACHI(y10, u, c.p6); y11 = VMACHI(y11, u, c.p7);
  z5 = VADD(VADD(z5, VSRA(z4, BRADIX)), VSHL(y4, BALIGN));

  u   = VAND(VMACLO(c.zero, z5, c.montw), c.bmask);
  z5  = VMACLO(z5, u, c.p0); z6  = VMACLO(z6, u, c.p1); z7  = VMACLO(z7, u, c.p2);
  z8  = VMACLO(z8, u, c.p3); z9  = VMACLO(z9, u, c.p4); z10 = VMACLO(z10, u, c.p5);
  z11 = VMACLO(z11, u, c.p6); z12 = VMACLO(z12, u, c.p7);
  y5  = VMACHI(y5, u, c.p0); y6  = VMACHI(y6, u, c.p1); y7  = VMACHI(y7, u, c.p2);
  y8  = VMACHI(y8, u, c.p3); y9  = VMACHI(y9, u, c.p4); y10 = VMACHI(y10, u, c.p5);
  y11 = VMACHI(y11, u, c.p6); y12 = VMACHI(y12, u, c.p7);
  z6 = VADD(VADD(z6, VSRA(z5, BRADIX)), VSHL(y5, BALIGN));

  u   = VAND(VMACLO(c.zero, z6, c.montw), c.bmask);
  z6  = VMACLO(z6, u, c.p0); z7  = VMACLO(z7, u, c.p1); z8  = VMACLO(z8, u, c.p2);
  z9  = VMACLO(z9, u, c.p3); z10 = VMACLO(z10, u, c.p4); z11 = VMACLO(z11, u, c.p5);
  z12 = VMACLO(z12, u, c.p6); z13 = VMACLO(z13, u, c.p7);
  y6  = VMACHI(y6, u, c.p0); y7  = VMACHI(y7, u, c.p1); y8  = VMACHI(y8, u, c.p2);
  y9  = VMACHI(y9, u, c.p3); y10 = VMACHI(y10, u, c.p4); y11 = VMACHI(y11, u, c.p5);
  y12 = VMACHI(y12, u, c.p6); y13 = VMACHI(y13, u, c.p7);
  z7 = VADD(VADD(z7, VSRA(z6, BRADIX)), VSHL(y6, BALIGN));

  u   = VAND(VMACLO(c.zero, z7, c.montw), c.bmask);
  z7  = VMACLO(z7, u, c.p0); z8  = VMACLO(z8, u, c.p1); z9  = VMACLO(z9, u, c.p2);
  z10 = VMACLO(z10, u, c.p3); z11 = VMACLO(z11, u, c.p4); z12 = VMACLO(z12, u, c.p5);
  z13 = VMACLO(z13, u, c.p6); z14 = VMACLO(z14, u, c.p7);
  y7  = VMACHI(y7, u, c.p0); y8  = VMACHI(y8, u, c.p1); y9  = VMACHI(y9, u, c.p2);
  y10 = VMACHI(y10, u, c.p3); y11 = VMACHI(y11, u, c.p4); y12 = VMACHI(y12, u, c.p5);
  y13 = VMACHI(y13, u, c.p6); y14 = VMACHI(y14, u, c.p7);
  z8 = VADD(VADD(z8, VSRA(z7, BRADIX)), VSHL(y7, BALIGN));

  z9  = VADD(z9, VSHL(y8, BALIGN));
  z10 = VADD(z10, VSHL(y9, BALIGN));
  z11 = VADD(z11, VSHL(y10, BALIGN));
  z12 = VADD(z12, VSHL(y11, BALIGN));
  z13 = VADD(z13, VSHL(y12, BALIGN));
  z14 = VADD(z14, VSHL(y13, BALIGN));
  z15 = VADD(z15, VSHL(y14, BALIGN));

  z9  = VADD(z9, VSRA(z8, BRADIX)); z8  = VAND(z8, c.bmask);
  z10 = VADD(z10, VSRA(z9, BRADIX)); z9  = VAND(z9, c.bmask);
  z11 = VADD(z11, VSRA(z10, BRADIX)); z10 = VAND(z10, c.bmask);
  z12 = VADD(z12, VSRA(z11, BRADIX)); z11 = VAND(z11, c.bmask);
  z13 = VADD(z13, VSRA(z12, BRADIX)); z12 = VAND(z12, c.bmask);
  z14 = VADD(z14, VSRA(z13, BRADIX)); z13 = VAND(z13, c.bmask);
  z15 = VADD(z15, VSRA(z14, BRADIX)); z14 = VAND(z14, c.bmask);
                                    z15 = VAND(z15, c.bmask);

  r[0] = z8; r[1] = z9; r[2] = z10; r[3] = z11;
  r[4] = z12; r[5] = z13; r[6] = z14; r[7] = z15;
}

bool fused_k1_selfcheck()
{
  constexpr std::size_t kPacks = 32;
  Fp8x1FusedBaselineData data(kPacks);
  std::vector<FpPack8x1> ref(kPacks);
  std::vector<FpPack8x1> got(kPacks);

  const Fp8x1MontConsts c{
      VSET1(P48[0]), VSET1(P48[1]), VSET1(P48[2]), VSET1(P48[3]),
      VSET1(P48[4]), VSET1(P48[5]), VSET1(P48[6]), VSET1(P48[7]),
      VSET1(BMASK), VSET1(MONT_W_R48), VZERO()};

  fp8x1_mul_redc_k<1>(ref.data(), data.a.data(), data.b.data(), kPacks);
  for (std::size_t i = 0; i < kPacks; ++i) {
    fp8x1_mul_redc_core_k1_fused_baseline(data.lazy[i], data.a[i], data.b[i], c);
    UFpPack8x1::normalize(got[i], data.lazy[i]);
  }

  for (std::size_t i = 0; i < kPacks; ++i) {
    if (!fp8x1_pack_equal(got[i], ref[i])) {
      return false;
    }
  }
  return true;
}

void fused_k1_require_selfcheck(benchmark::State& state)
{
  static std::once_flag once;
  static bool ok = true;
  std::call_once(once, []() { ok = fused_k1_selfcheck(); });
  if (!ok) {
    state.SkipWithError("Fused K1 self-check failed (mismatch vs fp8x1_mul_redc_k<1>).");
  }
}

void BM_Fp8x1MulRedcK1_Current(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1FusedBaselineData data(n_packs);

  for (auto _ : state) {
    fp8x1_mul_redc_k<1>(data.out.data(), data.a.data(), data.b.data(), n_packs);
    benchmark::DoNotOptimize(data.out.data());
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1MulRedcCoreK1_FusedBaseline(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1FusedBaselineData data(n_packs);
  const Fp8x1MontConsts c{
      VSET1(P48[0]), VSET1(P48[1]), VSET1(P48[2]), VSET1(P48[3]),
      VSET1(P48[4]), VSET1(P48[5]), VSET1(P48[6]), VSET1(P48[7]),
      VSET1(BMASK), VSET1(MONT_W_R48), VZERO()};
  fused_k1_require_selfcheck(state);
  if (state.skipped()) {
    return;
  }

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      fp8x1_mul_redc_core_k1_fused_baseline(data.lazy[i], data.a[i], data.b[i], c);
    }
    benchmark::DoNotOptimize(data.lazy.data());
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1MulRedcK1_FusedBaseline(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1FusedBaselineData data(n_packs);
  const Fp8x1MontConsts c{
      VSET1(P48[0]), VSET1(P48[1]), VSET1(P48[2]), VSET1(P48[3]),
      VSET1(P48[4]), VSET1(P48[5]), VSET1(P48[6]), VSET1(P48[7]),
      VSET1(BMASK), VSET1(MONT_W_R48), VZERO()};
  fused_k1_require_selfcheck(state);
  if (state.skipped()) {
    return;
  }

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      fp8x1_mul_redc_core_k1_fused_baseline(data.lazy[i], data.a[i], data.b[i], c);
      UFpPack8x1::normalize(data.out[i], data.lazy[i]);
    }
    benchmark::DoNotOptimize(data.out.data());
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}
} // namespace

BENCHMARK(BM_Fp8x1MulRedcK1_Current)->Arg(256);
BENCHMARK(BM_Fp8x1MulRedcCoreK1_FusedBaseline)->Arg(256);
BENCHMARK(BM_Fp8x1MulRedcK1_FusedBaseline)->Arg(256);

#endif
