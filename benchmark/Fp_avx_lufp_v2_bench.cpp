#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#include <array>
#include <cstddef>
#include <cstring>
#include <random>
#include <vector>

#include "blst/fields/Fp.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fp_avx.h"
#include "blst/fields/Fp6_avx.h"

using namespace blst;

namespace
{
enum class LUFpV2RenormClass
{
  ExactAdd,
  Fp6MulCommon,
};

void random_redc_fp(Fp& v, std::mt19937_64& gen)
{
  for (u64 i = 0; i < v.num_limb(); ++i) {
    v[i] = gen();
  }
  while (Fp::is_redc(v) == 0) {
    Fp::rshift(v, v, 1);
  }
}

struct LUFpPack8x1V2 : AvxPack<2 * NWORDS>
{
  static void add(LUFpPack8x1V2& r, const LUFpPack8x1V2& a, const LUFpPack8x1V2& b)
  {
    for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
      r[limb] = VADD(a[limb], b[limb]);
    }
  }

  static void sub(LUFpPack8x1V2& r, const LUFpPack8x1V2& a, const LUFpPack8x1V2& b)
  {
    for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
      r[limb] = VSUB(a[limb], b[limb]);
    }
  }

  template <LUFpV2RenormClass C>
  static void renorm(UFpPack8x1& r, const LUFpPack8x1V2& a);
};

struct LUFp2Pack8x1x1V2
{
  LUFpPack8x1V2 c[2];

  BLST_FORCEINLINE LUFpPack8x1V2& operator[](std::size_t i) { return c[i]; }
  BLST_FORCEINLINE const LUFpPack8x1V2& operator[](std::size_t i) const { return c[i]; }

  static void add(LUFp2Pack8x1x1V2& r, const LUFp2Pack8x1x1V2& a, const LUFp2Pack8x1x1V2& b)
  {
    LUFpPack8x1V2::add(r[0], a[0], b[0]);
    LUFpPack8x1V2::add(r[1], a[1], b[1]);
  }

  static void sub(LUFp2Pack8x1x1V2& r, const LUFp2Pack8x1x1V2& a, const LUFp2Pack8x1x1V2& b)
  {
    LUFpPack8x1V2::sub(r[0], a[0], b[0]);
    LUFpPack8x1V2::sub(r[1], a[1], b[1]);
  }

  static void mul_by_u_plus_1(LUFp2Pack8x1x1V2& r, const LUFp2Pack8x1x1V2& a)
  {
    LUFpPack8x1V2::sub(r[0], a[0], a[1]);
    LUFpPack8x1V2::add(r[1], a[0], a[1]);
  }

  template <LUFpV2RenormClass C>
  static void renorm(UFp2Pack8x1x1& r, const LUFp2Pack8x1x1V2& a)
  {
    LUFpPack8x1V2::renorm<C>(r[0], a[0]);
    LUFpPack8x1V2::renorm<C>(r[1], a[1]);
  }
};

std::array<u64, NWORDS> make_scaled_mod_plus_offset(const u64 scale, const u64 offset)
{
  std::array<u64, NWORDS> out{};
  u64 carry = offset;
  for (std::size_t i = 0; i < NWORDS; ++i) {
    const unsigned __int128 accum =
        static_cast<unsigned __int128>(P48[i]) * scale + carry;
    out[i] = static_cast<u64>(accum) & BMASK;
    carry = static_cast<u64>(accum >> BRADIX);
  }
  return out;
}

std::array<cc512i, NWORDS> broadcast_mod_vec(const std::array<u64, NWORDS>& mod)
{
  std::array<cc512i, NWORDS> out{};
  for (std::size_t i = 0; i < NWORDS; ++i) {
    out[i] = VSET1(static_cast<i64>(mod[i]));
  }
  return out;
}

const std::array<u64, NWORDS> kV2P48X1 = make_scaled_mod_plus_offset(1, 0);
const std::array<u64, NWORDS> kV2P48X2 = make_scaled_mod_plus_offset(2, 0);
const std::array<u64, NWORDS> kV2P48X3 = make_scaled_mod_plus_offset(3, 0);
const std::array<u64, NWORDS> kV2P48X4 = make_scaled_mod_plus_offset(4, 0);
const std::array<u64, NWORDS> kV2Bias3PPlus3 = make_scaled_mod_plus_offset(3, 3);

const std::array<cc512i, NWORDS> kV2P48X1Vec = broadcast_mod_vec(kV2P48X1);
const std::array<cc512i, NWORDS> kV2P48X2Vec = broadcast_mod_vec(kV2P48X2);
const std::array<cc512i, NWORDS> kV2P48X4Vec = broadcast_mod_vec(kV2P48X4);
const std::array<cc512i, NWORDS> kV2Bias3PPlus3Vec = broadcast_mod_vec(kV2Bias3PPlus3);

BLST_FORCEINLINE auto v2_mask_gt_epi64(const cc512i a, const cc512i b) -> u8
{
#if BLST_AVX_NATIVE
  return static_cast<u8>(_mm512_cmpgt_epi64_mask(to_native_512(a), to_native_512(b)));
#else
  std::array<u64, 8> lhs{};
  std::array<u64, 8> rhs{};
  _cc512_storeu_si512(lhs.data(), a);
  _cc512_storeu_si512(rhs.data(), b);
  u8 mask = 0;
  for (u8 lane = 0; lane < 8; ++lane) {
    if (static_cast<i64>(lhs[lane]) > static_cast<i64>(rhs[lane])) {
      mask |= static_cast<u8>(1u << lane);
    }
  }
  return mask;
}

BLST_FORCEINLINE void v2_copy_limbs(cc512i* dst, const LUFpPack8x1V2& src)
{
  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    dst[i] = src[i];
  }
}

BLST_FORCEINLINE void v2_store_wide(UFpPack8x1& dst, const cc512i* limbs)
{
  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    dst[i] = limbs[i];
  }
}

BLST_FORCEINLINE void v2_carry_range(cc512i* limbs, const std::size_t begin, const std::size_t end)
{
  const cc512i bmask = VSET1(BMASK);
  for (std::size_t i = begin; i < end; ++i) {
    limbs[i + 1] = VADD(limbs[i + 1], VSRA(limbs[i], BRADIX));
    limbs[i] = VAND(limbs[i], bmask);
  }
}

BLST_FORCEINLINE void v2_high_sub_const(cc512i* limbs, const std::array<cc512i, NWORDS>& c)
{
  for (std::size_t i = 0; i < NWORDS; ++i) {
    limbs[NWORDS + i] = VSUB(limbs[NWORDS + i], c[i]);
  }
}

BLST_FORCEINLINE void v2_high_add_const_masked(
    cc512i* limbs, const std::array<cc512i, NWORDS>& c, const u8 lane_mask)
{
  for (std::size_t i = 0; i < NWORDS; ++i) {
    limbs[NWORDS + i] = VMADD(limbs[NWORDS + i], lane_mask, limbs[NWORDS + i], c[i]);
  }
}

BLST_FORCEINLINE auto v2_neg_mask(const cc512i* limbs) -> u8
{
  return v2_mask_gt_epi64(VZERO(), limbs[(2 * NWORDS) - 1]);
}

BLST_FORCEINLINE void bench_make_v2_hi_neg4_pos4(
    LUFpPack8x1V2& out,
    const UFpPack8x1& w0,
    const UFpPack8x1& w1,
    const UFpPack8x1& w2,
    const UFpPack8x1& w3,
    const UFpPack8x1& w4,
    const UFpPack8x1& w5,
    const UFpPack8x1& w6,
    const bool plus_branch)
{
  LUFpPack8x1V2 pre0{};
  LUFpPack8x1V2 pre1{};
  LUFpPack8x1V2 tmp{};

  LUFpPack8x1V2::sub(pre0, reinterpret_cast<const LUFpPack8x1V2&>(w0), reinterpret_cast<const LUFpPack8x1V2&>(w1));
  LUFpPack8x1V2::sub(pre0, pre0, reinterpret_cast<const LUFpPack8x1V2&>(w2));

  LUFpPack8x1V2::sub(pre1, reinterpret_cast<const LUFpPack8x1V2&>(w3), reinterpret_cast<const LUFpPack8x1V2&>(w4));
  LUFpPack8x1V2::sub(pre1, pre1, reinterpret_cast<const LUFpPack8x1V2&>(w5));

  if (plus_branch) {
    LUFpPack8x1V2::add(tmp, pre0, pre1);
  } else {
    LUFpPack8x1V2::sub(tmp, pre0, pre1);
  }

  LUFpPack8x1V2::add(out, tmp, reinterpret_cast<const LUFpPack8x1V2&>(w6));
}

template <>
void LUFpPack8x1V2::renorm<LUFpV2RenormClass::ExactAdd>(UFpPack8x1& r, const LUFpPack8x1V2& a)
{
  cc512i limbs[2 * NWORDS]{};
  v2_copy_limbs(limbs, a);

  v2_high_sub_const(limbs, kV2P48X1Vec);
  v2_carry_range(limbs, 0, (2 * NWORDS) - 1);

  const u8 neg_mask = v2_neg_mask(limbs);
  v2_high_add_const_masked(limbs, kV2P48X1Vec, neg_mask);
  v2_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  v2_store_wide(r, limbs);
}

template <>
void LUFpPack8x1V2::renorm<LUFpV2RenormClass::Fp6MulCommon>(UFpPack8x1& r, const LUFpPack8x1V2& a)
{
  cc512i limbs[2 * NWORDS]{};
  v2_copy_limbs(limbs, a);

  // Bias by B = 3p + 3 so the carried class lands inside [-8p, p).
  v2_high_sub_const(limbs, kV2Bias3PPlus3Vec);
  v2_carry_range(limbs, 0, (2 * NWORDS) - 1);

  // Sign-driven ladder for [-8p, p): +4p, +2p, +p, then a final +p cleanup.
  u8 neg_mask = v2_neg_mask(limbs);
  v2_high_add_const_masked(limbs, kV2P48X4Vec, neg_mask);
  v2_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  neg_mask = v2_neg_mask(limbs);
  v2_high_add_const_masked(limbs, kV2P48X2Vec, neg_mask);
  v2_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  neg_mask = v2_neg_mask(limbs);
  v2_high_add_const_masked(limbs, kV2P48X1Vec, neg_mask);
  v2_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  neg_mask = v2_neg_mask(limbs);
  v2_high_add_const_masked(limbs, kV2P48X1Vec, neg_mask);
  v2_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);

  v2_store_wide(r, limbs);
}

BLST_FORCEINLINE void bench_copy_as_v2(LUFpPack8x1V2& dst, const UFpPack8x1& src)
{
  for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
    dst[limb] = src[limb];
  }
}

BLST_FORCEINLINE void bench_copy_as_v2(LUFp2Pack8x1x1V2& dst, const UFp2Pack8x1x1& src)
{
  bench_copy_as_v2(dst[0], src[0]);
  bench_copy_as_v2(dst[1], src[1]);
}

BLST_FORCEINLINE void bench_consume_fp6_pack(std::vector<Fp6Pack8x1>& data)
{
  if (data.empty()) {
    return;
  }
  benchmark::DoNotOptimize(data.front().c[0][0][0]);
  benchmark::DoNotOptimize(data.back().c[2][1][NWORDS - 1]);
  benchmark::ClobberMemory();
}

struct Fp8x1LUFpV2BenchData
{
  explicit Fp8x1LUFpV2BenchData(const std::size_t n_packs)
      : n_packs_(n_packs)
      , wide_a(n_packs_)
      , wide_b(n_packs_)
      , lazy_a(n_packs_)
      , lazy_b(n_packs_)
      , lazy_tmp(n_packs_)
      , lazy_fp6_common(n_packs_)
      , wide_out(n_packs_)
  {
    std::mt19937_64 gen(0xD2C2B3A41516ULL + static_cast<u64>(n_packs_));

    for (std::size_t i = 0; i < n_packs_; ++i) {
      std::array<Fp, 8> lane_a{};
      std::array<Fp, 8> lane_b{};
      std::array<Fp, 8> lane_c{};
      std::array<Fp, 8> lane_d{};
      FpPack8x1 pack_a{}, pack_b{}, pack_c{}, pack_d{};

      for (std::size_t lane = 0; lane < 8; ++lane) {
        random_redc_fp(lane_a[lane], gen);
        random_redc_fp(lane_b[lane], gen);
        random_redc_fp(lane_c[lane], gen);
        random_redc_fp(lane_d[lane], gen);
      }

      FpPack8x1::from_aos(pack_a, lane_a.data());
      FpPack8x1::from_aos(pack_b, lane_b.data());
      FpPack8x1::from_aos(pack_c, lane_c.data());
      FpPack8x1::from_aos(pack_d, lane_d.data());
      FpPack8x1::mul(wide_a[i], pack_a, pack_b);
      FpPack8x1::mul(wide_b[i], pack_c, pack_d);
      bench_copy_as_v2(lazy_a[i], wide_a[i]);
      bench_copy_as_v2(lazy_b[i], wide_b[i]);

      UFpPack8x1 wide_c{}, wide_d{}, wide_e{}, wide_f{}, wide_g{};
      FpPack8x1::mul(wide_c, pack_a, pack_c);
      FpPack8x1::mul(wide_d, pack_b, pack_d);
      FpPack8x1::mul(wide_e, pack_a, pack_d);
      FpPack8x1::mul(wide_f, pack_b, pack_c);
      FpPack8x1::mul(wide_g, pack_a, pack_a);
      bench_make_v2_hi_neg4_pos4(
          lazy_fp6_common[i],
          wide_a[i],
          wide_b[i],
          wide_c,
          wide_d,
          wide_e,
          wide_f,
          wide_g,
          (i & 1u) == 0u);
    }
  }

  const std::size_t n_packs_;
  std::vector<UFpPack8x1> wide_a;
  std::vector<UFpPack8x1> wide_b;
  std::vector<LUFpPack8x1V2> lazy_a;
  std::vector<LUFpPack8x1V2> lazy_b;
  std::vector<LUFpPack8x1V2> lazy_tmp;
  std::vector<LUFpPack8x1V2> lazy_fp6_common;
  std::vector<UFpPack8x1> wide_out;
};

struct Fp6Pack8x1V2BenchData
{
  explicit Fp6Pack8x1V2BenchData(const std::size_t n_packs)
      : n_packs_(n_packs)
      , a_scalar(n_packs_ * 8)
      , b_scalar(n_packs_ * 8)
      , a_pack(n_packs_)
      , b_pack(n_packs_)
      , out_pack(n_packs_)
  {
    std::mt19937_64 gen(0x8C4F7D2E19B3ULL + static_cast<u64>(n_packs_));

    for (std::size_t i = 0; i < n_packs_ * 8; ++i) {
      for (std::size_t c = 0; c < 3; ++c) {
        for (std::size_t d = 0; d < 2; ++d) {
          random_redc_fp(a_scalar[i][c][d], gen);
          random_redc_fp(b_scalar[i][c][d], gen);
        }
      }
    }

    for (std::size_t p = 0; p < n_packs_; ++p) {
      Fp6Pack8x1::from_aos(a_pack[p], a_scalar.data() + (8 * p));
      Fp6Pack8x1::from_aos(b_pack[p], b_scalar.data() + (8 * p));
    }
  }

  const std::size_t n_packs_;
  std::vector<Fp6> a_scalar;
  std::vector<Fp6> b_scalar;
  std::vector<Fp6Pack8x1> a_pack;
  std::vector<Fp6Pack8x1> b_pack;
  std::vector<Fp6Pack8x1> out_pack;
};

struct Fp6WideCombineV2BenchData
{
  explicit Fp6WideCombineV2BenchData(const std::size_t n_packs)
      : n_packs_(n_packs)
      , a_scalar(n_packs_ * 8)
      , b_scalar(n_packs_ * 8)
      , t0_pack(n_packs_)
      , t1_pack(n_packs_)
      , t2_pack(n_packs_)
      , w0_pack(n_packs_)
      , w1_pack(n_packs_)
      , w2_pack(n_packs_)
      , vt0_pack(n_packs_)
      , vt1_pack(n_packs_)
      , vt2_pack(n_packs_)
      , vw0_pack(n_packs_)
      , vw1_pack(n_packs_)
      , vw2_pack(n_packs_)
      , lazy_c0_pack(n_packs_)
      , lazy_c1_pack(n_packs_)
      , lazy_c2_pack(n_packs_)
      , out_wide_pack(n_packs_)
  {
    std::mt19937_64 gen(0x91F63A4DE8B2ULL + static_cast<u64>(n_packs_));

    for (std::size_t i = 0; i < n_packs_ * 8; ++i) {
      for (std::size_t c = 0; c < 3; ++c) {
        for (std::size_t d = 0; d < 2; ++d) {
          random_redc_fp(a_scalar[i][c][d], gen);
          random_redc_fp(b_scalar[i][c][d], gen);
        }
      }
    }

    for (std::size_t p = 0; p < n_packs_; ++p) {
      Fp6Pack8x1 a_pack{};
      Fp6Pack8x1 b_pack{};
      Fp2Pack8x1x1 aa{}, bb{};

      Fp6Pack8x1::from_aos(a_pack, a_scalar.data() + (8 * p));
      Fp6Pack8x1::from_aos(b_pack, b_scalar.data() + (8 * p));

      UFp2Pack8x1x1::mul(t0_pack[p], a_pack.c[0], b_pack.c[0]);
      UFp2Pack8x1x1::mul(t1_pack[p], a_pack.c[1], b_pack.c[1]);
      UFp2Pack8x1x1::mul(t2_pack[p], a_pack.c[2], b_pack.c[2]);

      Fp2Pack8x1x1::add(aa, a_pack.c[1], a_pack.c[2]);
      Fp2Pack8x1x1::add(bb, b_pack.c[1], b_pack.c[2]);
      UFp2Pack8x1x1::mul(w0_pack[p], aa, bb);

      Fp2Pack8x1x1::add(aa, a_pack.c[0], a_pack.c[1]);
      Fp2Pack8x1x1::add(bb, b_pack.c[0], b_pack.c[1]);
      UFp2Pack8x1x1::mul(w1_pack[p], aa, bb);

      Fp2Pack8x1x1::add(aa, a_pack.c[0], a_pack.c[2]);
      Fp2Pack8x1x1::add(bb, b_pack.c[0], b_pack.c[2]);
      UFp2Pack8x1x1::mul(w2_pack[p], aa, bb);

      bench_copy_as_v2(vt0_pack[p], t0_pack[p]);
      bench_copy_as_v2(vt1_pack[p], t1_pack[p]);
      bench_copy_as_v2(vt2_pack[p], t2_pack[p]);
      bench_copy_as_v2(vw0_pack[p], w0_pack[p]);
      bench_copy_as_v2(vw1_pack[p], w1_pack[p]);
      bench_copy_as_v2(vw2_pack[p], w2_pack[p]);

      LUFp2Pack8x1x1V2 t3{};

      lazy_c0_pack[p] = vw0_pack[p];
      LUFp2Pack8x1x1V2::sub(lazy_c0_pack[p], lazy_c0_pack[p], vt1_pack[p]);
      LUFp2Pack8x1x1V2::sub(lazy_c0_pack[p], lazy_c0_pack[p], vt2_pack[p]);
      LUFp2Pack8x1x1V2::mul_by_u_plus_1(t3, lazy_c0_pack[p]);
      LUFp2Pack8x1x1V2::add(lazy_c0_pack[p], t3, vt0_pack[p]);

      lazy_c1_pack[p] = vw1_pack[p];
      LUFp2Pack8x1x1V2::sub(lazy_c1_pack[p], lazy_c1_pack[p], vt0_pack[p]);
      LUFp2Pack8x1x1V2::sub(lazy_c1_pack[p], lazy_c1_pack[p], vt1_pack[p]);
      t3 = vt2_pack[p];
      LUFp2Pack8x1x1V2::mul_by_u_plus_1(t3, t3);
      LUFp2Pack8x1x1V2::add(lazy_c1_pack[p], lazy_c1_pack[p], t3);

      lazy_c2_pack[p] = vw2_pack[p];
      LUFp2Pack8x1x1V2::sub(lazy_c2_pack[p], lazy_c2_pack[p], vt0_pack[p]);
      LUFp2Pack8x1x1V2::sub(lazy_c2_pack[p], lazy_c2_pack[p], vt2_pack[p]);
      LUFp2Pack8x1x1V2::add(lazy_c2_pack[p], lazy_c2_pack[p], vt1_pack[p]);
    }
  }

  const std::size_t n_packs_;
  std::vector<Fp6> a_scalar;
  std::vector<Fp6> b_scalar;
  std::vector<UFp2Pack8x1x1> t0_pack;
  std::vector<UFp2Pack8x1x1> t1_pack;
  std::vector<UFp2Pack8x1x1> t2_pack;
  std::vector<UFp2Pack8x1x1> w0_pack;
  std::vector<UFp2Pack8x1x1> w1_pack;
  std::vector<UFp2Pack8x1x1> w2_pack;
  std::vector<LUFp2Pack8x1x1V2> vt0_pack;
  std::vector<LUFp2Pack8x1x1V2> vt1_pack;
  std::vector<LUFp2Pack8x1x1V2> vt2_pack;
  std::vector<LUFp2Pack8x1x1V2> vw0_pack;
  std::vector<LUFp2Pack8x1x1V2> vw1_pack;
  std::vector<LUFp2Pack8x1x1V2> vw2_pack;
  std::vector<LUFp2Pack8x1x1V2> lazy_c0_pack;
  std::vector<LUFp2Pack8x1x1V2> lazy_c1_pack;
  std::vector<LUFp2Pack8x1x1V2> lazy_c2_pack;
  std::vector<UFp6Pack8x1> out_wide_pack;
};

void mul_fp6_pack_lufpv2_combine_only_raw(
    LUFp2Pack8x1x1V2& c0,
    LUFp2Pack8x1x1V2& c1,
    LUFp2Pack8x1x1V2& c2,
    const LUFp2Pack8x1x1V2& t0,
    const LUFp2Pack8x1x1V2& t1,
    const LUFp2Pack8x1x1V2& t2,
    const LUFp2Pack8x1x1V2& w0,
    const LUFp2Pack8x1x1V2& w1,
    const LUFp2Pack8x1x1V2& w2)
{
  LUFp2Pack8x1x1V2 t3{};

  c0 = w0;
  LUFp2Pack8x1x1V2::sub(c0, c0, t1);
  LUFp2Pack8x1x1V2::sub(c0, c0, t2);
  LUFp2Pack8x1x1V2::mul_by_u_plus_1(t3, c0);
  LUFp2Pack8x1x1V2::add(c0, t3, t0);

  c1 = w1;
  LUFp2Pack8x1x1V2::sub(c1, c1, t0);
  LUFp2Pack8x1x1V2::sub(c1, c1, t1);
  t3 = t2;
  LUFp2Pack8x1x1V2::mul_by_u_plus_1(t3, t3);
  LUFp2Pack8x1x1V2::add(c1, c1, t3);

  c2 = w2;
  LUFp2Pack8x1x1V2::sub(c2, c2, t0);
  LUFp2Pack8x1x1V2::sub(c2, c2, t2);
  LUFp2Pack8x1x1V2::add(c2, c2, t1);
}

void mul_fp6_pack_lufpv2_combine_only_wide(
    UFp6Pack8x1& r,
    const UFp2Pack8x1x1& t0,
    const UFp2Pack8x1x1& t1,
    const UFp2Pack8x1x1& t2,
    const UFp2Pack8x1x1& w0,
    const UFp2Pack8x1x1& w1,
    const UFp2Pack8x1x1& w2)
{
  LUFp2Pack8x1x1V2 lazy{}, t3{}, vt0{}, vt1{}, vt2{};
  bench_copy_as_v2(vt0, t0);
  bench_copy_as_v2(vt1, t1);
  bench_copy_as_v2(vt2, t2);

  bench_copy_as_v2(lazy, w0);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt1);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt2);
  LUFp2Pack8x1x1V2::mul_by_u_plus_1(t3, lazy);
  LUFp2Pack8x1x1V2::add(lazy, t3, vt0);
  LUFp2Pack8x1x1V2::renorm<LUFpV2RenormClass::Fp6MulCommon>(r.c[0], lazy);

  bench_copy_as_v2(lazy, w1);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt0);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt1);
  t3 = vt2;
  LUFp2Pack8x1x1V2::mul_by_u_plus_1(t3, t3);
  LUFp2Pack8x1x1V2::add(lazy, lazy, t3);
  LUFp2Pack8x1x1V2::renorm<LUFpV2RenormClass::Fp6MulCommon>(r.c[1], lazy);

  bench_copy_as_v2(lazy, w2);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt0);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt2);
  LUFp2Pack8x1x1V2::add(lazy, lazy, vt1);
  LUFp2Pack8x1x1V2::renorm<LUFpV2RenormClass::Fp6MulCommon>(r.c[2], lazy);
}

void mul_fp6_pack_lufpv2_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
  Fp2Pack8x1x1 aa{}, bb{};
  UFp2Pack8x1x1 t0{}, t1{}, t2{}, wide{};
  LUFp2Pack8x1x1V2 lazy{}, t3{}, vt0{}, vt1{}, vt2{};

  UFp2Pack8x1x1::mul(t0, a.c[0], b.c[0]);
  UFp2Pack8x1x1::mul(t1, a.c[1], b.c[1]);
  UFp2Pack8x1x1::mul(t2, a.c[2], b.c[2]);
  bench_copy_as_v2(vt0, t0);
  bench_copy_as_v2(vt1, t1);
  bench_copy_as_v2(vt2, t2);

  Fp2Pack8x1x1::add(aa, a.c[1], a.c[2]);
  Fp2Pack8x1x1::add(bb, b.c[1], b.c[2]);
  UFp2Pack8x1x1::mul(wide, aa, bb);
  bench_copy_as_v2(lazy, wide);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt1);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt2);
  LUFp2Pack8x1x1V2::mul_by_u_plus_1(t3, lazy);
  LUFp2Pack8x1x1V2::add(lazy, t3, vt0);
  LUFp2Pack8x1x1V2::renorm<LUFpV2RenormClass::Fp6MulCommon>(r.c[0], lazy);

  Fp2Pack8x1x1::add(aa, a.c[0], a.c[1]);
  Fp2Pack8x1x1::add(bb, b.c[0], b.c[1]);
  UFp2Pack8x1x1::mul(wide, aa, bb);
  bench_copy_as_v2(lazy, wide);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt0);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt1);
  t3 = vt2;
  LUFp2Pack8x1x1V2::mul_by_u_plus_1(t3, t3);
  LUFp2Pack8x1x1V2::add(lazy, lazy, t3);
  LUFp2Pack8x1x1V2::renorm<LUFpV2RenormClass::Fp6MulCommon>(r.c[1], lazy);

  Fp2Pack8x1x1::add(aa, a.c[0], a.c[2]);
  Fp2Pack8x1x1::add(bb, b.c[0], b.c[2]);
  UFp2Pack8x1x1::mul(wide, aa, bb);
  bench_copy_as_v2(lazy, wide);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt0);
  LUFp2Pack8x1x1V2::sub(lazy, lazy, vt2);
  LUFp2Pack8x1x1V2::add(lazy, lazy, vt1);
  LUFp2Pack8x1x1V2::renorm<LUFpV2RenormClass::Fp6MulCommon>(r.c[2], lazy);
}

void mul_fp6_pack_lufpv2(Fp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
  UFp6Pack8x1 wide{};
  mul_fp6_pack_lufpv2_wide(wide, a, b);
  UFp6Pack8x1::redc(r, wide);
}

BLST_FORCEINLINE void bench_consume_wide_v2(const std::vector<UFpPack8x1>& data)
{
  if (data.empty()) {
    return;
  }
  benchmark::DoNotOptimize(data.front()[0]);
  benchmark::DoNotOptimize(data.back()[(2 * NWORDS) - 1]);
  benchmark::ClobberMemory();
}

BLST_FORCEINLINE void bench_consume_luwide_v2(const std::vector<LUFpPack8x1V2>& data)
{
  if (data.empty()) {
    return;
  }
  benchmark::DoNotOptimize(data.front()[0]);
  benchmark::DoNotOptimize(data.back()[(2 * NWORDS) - 1]);
  benchmark::ClobberMemory();
}

} // namespace

#endif
