#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include <array>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <random>
#include <vector>

#include "blst/fields/Fp.h"
#include "blst/fields/Fp_avx.h"

using namespace blst;

#if defined(_MSC_VER)
#define BLST_BENCH_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define BLST_BENCH_NOINLINE __attribute__((noinline))
#else
#define BLST_BENCH_NOINLINE
#endif

#define BLST_HAVE_MUL_V2_K2_V6B_ASM 0

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
      FpPack8x1::mul(wide[i], a[i], b[i]);
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

BLST_FORCEINLINE void bench_copy_as_lufp(LUFpPack8x1& dst, const UFpPack8x1& src)
{
  for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
    dst[limb] = src[limb];
  }
}

BLST_FORCEINLINE void bench_make_lufp_hi_neg4_pos4(
    LUFpPack8x1& out,
    const UFpPack8x1& w0,
    const UFpPack8x1& w1,
    const UFpPack8x1& w2,
    const UFpPack8x1& w3,
    const UFpPack8x1& w4,
    const UFpPack8x1& w5,
    const UFpPack8x1& w6,
    const bool plus_branch)
{
  LUFpPack8x1 pre0{};
  LUFpPack8x1 pre1{};
  LUFpPack8x1 tmp{};

  LUFpPack8x1::sub(pre0, w0, w1);
  LUFpPack8x1::sub(pre0, pre0, w2);

  LUFpPack8x1::sub(pre1, w3, w4);
  LUFpPack8x1::sub(pre1, pre1, w5);

  if (plus_branch) {
    LUFpPack8x1::add(tmp, pre0, pre1);
  } else {
    LUFpPack8x1::sub(tmp, pre0, pre1);
  }

  LUFpPack8x1::add(out, tmp, w6);
}

BLST_FORCEINLINE void bench_copy_wide(UFpPack8x1& dst, const UFpPack8x1& src)
{
  for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
    dst[limb] = src[limb];
  }
}

BLST_FORCEINLINE void bench_add_wide_raw(UFpPack8x1& dst, const UFpPack8x1& a, const UFpPack8x1& b)
{
  for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
    dst[limb] = VADD(a[limb], b[limb]);
  }
}

std::array<u64, NWORDS> bench_make_scaled_mod(const u64 scale)
{
  std::array<u64, NWORDS> out{};
  u64 carry = 0;
  for (std::size_t i = 0; i < NWORDS; ++i) {
    const unsigned __int128 accum =
        static_cast<unsigned __int128>(P48[i]) * scale + carry;
    out[i] = static_cast<u64>(accum) & BMASK;
    carry = static_cast<u64>(accum >> BRADIX);
  }
  return out;
}

const std::array<u64, NWORDS> kBenchP48X1 = bench_make_scaled_mod(1);
const std::array<u64, NWORDS> kBenchP48X2 = bench_make_scaled_mod(2);
const std::array<u64, NWORDS> kBenchP48X3 = bench_make_scaled_mod(3);
const std::array<u64, NWORDS> kBenchP48X4 = bench_make_scaled_mod(4);
const std::array<u64, NWORDS> kBenchP48X7 = bench_make_scaled_mod(7);

std::array<cc512i, NWORDS> bench_broadcast_mod(const std::array<u64, NWORDS>& mod)
{
  std::array<cc512i, NWORDS> out{};
  for (std::size_t i = 0; i < NWORDS; ++i) {
    out[i] = VSET1(static_cast<i64>(mod[i]));
  }
  return out;
}

std::array<cc512i, NWORDS> bench_make_lane_varying_vec(
    const std::array<u64, NWORDS>& base, const u64 lane_step)
{
  std::array<cc512i, NWORDS> out{};
  for (std::size_t i = 0; i < NWORDS; ++i) {
    const u64 b = base[i];
    out[i] = VSET(
        static_cast<i64>(b + 7u * lane_step),
        static_cast<i64>(b + 6u * lane_step),
        static_cast<i64>(b + 5u * lane_step),
        static_cast<i64>(b + 4u * lane_step),
        static_cast<i64>(b + 3u * lane_step),
        static_cast<i64>(b + 2u * lane_step),
        static_cast<i64>(b + 1u * lane_step),
        static_cast<i64>(b + 0u * lane_step));
  }
  return out;
}

const std::array<cc512i, NWORDS> kBenchP48X1Vec = bench_broadcast_mod(kBenchP48X1);
const std::array<cc512i, NWORDS> kBenchP48X2Vec = bench_broadcast_mod(kBenchP48X2);
const std::array<cc512i, NWORDS> kBenchP48X3Vec = bench_broadcast_mod(kBenchP48X3);
const std::array<cc512i, NWORDS> kBenchP48X4Vec = bench_broadcast_mod(kBenchP48X4);
const std::array<cc512i, NWORDS> kBenchP48X7Vec = bench_broadcast_mod(kBenchP48X7);
const std::array<u64, NWORDS> kBenchCompareSeed64 = {
    0x000000001111ULL, 0x000000022222ULL, 0x000000033333ULL, 0x000000044444ULL,
    0x000000055555ULL, 0x000000066666ULL, 0x000000077777ULL, 0x000000088888ULL,
};
const std::array<u64, NWORDS> kBenchCompareKick64 = {
    0x000000000123ULL, 0x000000004567ULL, 0x000000089ABCULL, 0x0000000DEF01ULL,
    0x000000123456ULL, 0x0000000789ABULL, 0x0000000CDEF0ULL, 0x00000013579BULL,
};
const std::array<cc512i, NWORDS> kBenchCompareSeedVec = bench_make_lane_varying_vec(kBenchCompareSeed64, 0x11ULL);
const std::array<cc512i, NWORDS> kBenchCompareKickVec = bench_make_lane_varying_vec(kBenchCompareKick64, 0x101ULL);

template <u64 Scale>
BLST_FORCEINLINE const std::array<cc512i, NWORDS>& bench_mod_vec();

template <>
BLST_FORCEINLINE const std::array<cc512i, NWORDS>& bench_mod_vec<1>() { return kBenchP48X1Vec; }
template <>
BLST_FORCEINLINE const std::array<cc512i, NWORDS>& bench_mod_vec<2>() { return kBenchP48X2Vec; }
template <>
BLST_FORCEINLINE const std::array<cc512i, NWORDS>& bench_mod_vec<3>() { return kBenchP48X3Vec; }
template <>
BLST_FORCEINLINE const std::array<cc512i, NWORDS>& bench_mod_vec<4>() { return kBenchP48X4Vec; }
template <>
BLST_FORCEINLINE const std::array<cc512i, NWORDS>& bench_mod_vec<7>() { return kBenchP48X7Vec; }

BLST_FORCEINLINE void bench_mutate_compare_high(cc512i* limbs)
{
  cc512i prev = limbs[(2 * NWORDS) - 1];
  for (std::size_t i = 0; i < NWORDS; ++i) {
    const cc512i cur = limbs[NWORDS + i];
    limbs[NWORDS + i] = VXOR(VADD(cur, kBenchCompareKickVec[i]), prev);
    prev = cur;
  }
}

BLST_FORCEINLINE void bench_seed_compare_high(cc512i* limbs)
{
  for (std::size_t i = 0; i < NWORDS; ++i) {
    limbs[NWORDS + i] = kBenchCompareSeedVec[i];
  }
}

BLST_FORCEINLINE void bench_signed_carry_range(cc512i* limbs, const std::size_t begin, const std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i) {
    const cc512i carry = VSRA(limbs[i], BRADIX);
    limbs[i] = VSUB(limbs[i], VSHL(carry, BRADIX));
    limbs[i + 1] = VADD(limbs[i + 1], carry);
  }
}

BLST_FORCEINLINE auto bench_mask_eq_epi64(const cc512i a, const cc512i b) -> u8
{
#if BLST_AVX_NATIVE
  return static_cast<u8>(_mm512_cmpeq_epi64_mask(to_native_512(a), to_native_512(b)));
#else
  std::array<u64, 8> lhs{};
  std::array<u64, 8> rhs{};
  _cc512_storeu_si512(lhs.data(), a);
  _cc512_storeu_si512(rhs.data(), b);
  u8 mask = 0;
  for (u8 lane = 0; lane < 8; ++lane) {
    if (lhs[lane] == rhs[lane]) {
      mask |= static_cast<u8>(1u << lane);
    }
  }
  return mask;
#endif
}

BLST_FORCEINLINE auto bench_mask_gt_epi64(const cc512i a, const cc512i b) -> u8
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
#endif
}

struct BenchGtEq
{
  u8 gt;
  u8 eq;
};

struct BenchThreeWayBuckets
{
  u8 ge_p;
  u8 ge_2p;
  u8 ge_3p;
};

BLST_FORCEINLINE auto bench_gt_eq_pair(const cc512i a, const cc512i b) -> BenchGtEq
{
#if BLST_AVX_NATIVE
  return BenchGtEq{
      static_cast<u8>(_mm512_cmpgt_epi64_mask(to_native_512(a), to_native_512(b))),
      static_cast<u8>(_mm512_cmpeq_epi64_mask(to_native_512(a), to_native_512(b))),
  };
#else
  return BenchGtEq{
      bench_mask_gt_epi64(a, b),
      bench_mask_eq_epi64(a, b),
  };
#endif
}

BLST_FORCEINLINE auto bench_gt_eq_combine_hi_lo(const BenchGtEq hi, const BenchGtEq lo) -> BenchGtEq
{
#if BLST_AVX_NATIVE
  const __mmask8 hi_gt = static_cast<__mmask8>(hi.gt);
  const __mmask8 hi_eq = static_cast<__mmask8>(hi.eq);
  const __mmask8 lo_gt = static_cast<__mmask8>(lo.gt);
  const __mmask8 lo_eq = static_cast<__mmask8>(lo.eq);
  return BenchGtEq{
      static_cast<u8>(_kor_mask8(hi_gt, _kand_mask8(hi_eq, lo_gt))),
      static_cast<u8>(_kand_mask8(hi_eq, lo_eq)),
  };
#else
  return BenchGtEq{
      static_cast<u8>(hi.gt | (hi.eq & lo.gt)),
      static_cast<u8>(hi.eq & lo.eq),
  };
#endif
}

template <u64 Scale>
BLST_FORCEINLINE auto bench_high_ge_mod_mask_prefix_lex(const cc512i* limbs) -> u8
{
  const auto& mod_vec = bench_mod_vec<Scale>();

  const BenchGtEq c7 = bench_gt_eq_pair(limbs[NWORDS + 7], mod_vec[7]);
  const BenchGtEq c6 = bench_gt_eq_pair(limbs[NWORDS + 6], mod_vec[6]);
  const BenchGtEq c5 = bench_gt_eq_pair(limbs[NWORDS + 5], mod_vec[5]);
  const BenchGtEq c4 = bench_gt_eq_pair(limbs[NWORDS + 4], mod_vec[4]);
  const BenchGtEq c3 = bench_gt_eq_pair(limbs[NWORDS + 3], mod_vec[3]);
  const BenchGtEq c2 = bench_gt_eq_pair(limbs[NWORDS + 2], mod_vec[2]);
  const BenchGtEq c1 = bench_gt_eq_pair(limbs[NWORDS + 1], mod_vec[1]);
  const BenchGtEq c0 = bench_gt_eq_pair(limbs[NWORDS + 0], mod_vec[0]);

  const BenchGtEq c76 = bench_gt_eq_combine_hi_lo(c7, c6);
  const BenchGtEq c54 = bench_gt_eq_combine_hi_lo(c5, c4);
  const BenchGtEq c32 = bench_gt_eq_combine_hi_lo(c3, c2);
  const BenchGtEq c10 = bench_gt_eq_combine_hi_lo(c1, c0);

  const BenchGtEq c7654 = bench_gt_eq_combine_hi_lo(c76, c54);
  const BenchGtEq c3210 = bench_gt_eq_combine_hi_lo(c32, c10);
  const BenchGtEq all = bench_gt_eq_combine_hi_lo(c7654, c3210);

#if BLST_AVX_NATIVE
  return static_cast<u8>(_kor_mask8(static_cast<__mmask8>(all.gt), static_cast<__mmask8>(all.eq)));
#else
  return static_cast<u8>(all.gt | all.eq);
#endif
}

BLST_FORCEINLINE auto bench_high_ge_p_2p_3p_masks_fused(const cc512i* limbs)
    -> BenchThreeWayBuckets
{
  const BenchGtEq p7 = bench_gt_eq_pair(limbs[NWORDS + 7], kBenchP48X1Vec[7]);
  const BenchGtEq p6 = bench_gt_eq_pair(limbs[NWORDS + 6], kBenchP48X1Vec[6]);
  const BenchGtEq p5 = bench_gt_eq_pair(limbs[NWORDS + 5], kBenchP48X1Vec[5]);
  const BenchGtEq p4 = bench_gt_eq_pair(limbs[NWORDS + 4], kBenchP48X1Vec[4]);
  const BenchGtEq p3 = bench_gt_eq_pair(limbs[NWORDS + 3], kBenchP48X1Vec[3]);
  const BenchGtEq p2 = bench_gt_eq_pair(limbs[NWORDS + 2], kBenchP48X1Vec[2]);
  const BenchGtEq p1 = bench_gt_eq_pair(limbs[NWORDS + 1], kBenchP48X1Vec[1]);
  const BenchGtEq p0 = bench_gt_eq_pair(limbs[NWORDS + 0], kBenchP48X1Vec[0]);

  const BenchGtEq p27 = bench_gt_eq_pair(limbs[NWORDS + 7], kBenchP48X2Vec[7]);
  const BenchGtEq p26 = bench_gt_eq_pair(limbs[NWORDS + 6], kBenchP48X2Vec[6]);
  const BenchGtEq p25 = bench_gt_eq_pair(limbs[NWORDS + 5], kBenchP48X2Vec[5]);
  const BenchGtEq p24 = bench_gt_eq_pair(limbs[NWORDS + 4], kBenchP48X2Vec[4]);
  const BenchGtEq p23 = bench_gt_eq_pair(limbs[NWORDS + 3], kBenchP48X2Vec[3]);
  const BenchGtEq p22 = bench_gt_eq_pair(limbs[NWORDS + 2], kBenchP48X2Vec[2]);
  const BenchGtEq p21 = bench_gt_eq_pair(limbs[NWORDS + 1], kBenchP48X2Vec[1]);
  const BenchGtEq p20 = bench_gt_eq_pair(limbs[NWORDS + 0], kBenchP48X2Vec[0]);

  const BenchGtEq p37 = bench_gt_eq_pair(limbs[NWORDS + 7], kBenchP48X3Vec[7]);
  const BenchGtEq p36 = bench_gt_eq_pair(limbs[NWORDS + 6], kBenchP48X3Vec[6]);
  const BenchGtEq p35 = bench_gt_eq_pair(limbs[NWORDS + 5], kBenchP48X3Vec[5]);
  const BenchGtEq p34 = bench_gt_eq_pair(limbs[NWORDS + 4], kBenchP48X3Vec[4]);
  const BenchGtEq p33 = bench_gt_eq_pair(limbs[NWORDS + 3], kBenchP48X3Vec[3]);
  const BenchGtEq p32 = bench_gt_eq_pair(limbs[NWORDS + 2], kBenchP48X3Vec[2]);
  const BenchGtEq p31 = bench_gt_eq_pair(limbs[NWORDS + 1], kBenchP48X3Vec[1]);
  const BenchGtEq p30 = bench_gt_eq_pair(limbs[NWORDS + 0], kBenchP48X3Vec[0]);

  const BenchGtEq p76 = bench_gt_eq_combine_hi_lo(p7, p6);
  const BenchGtEq p54 = bench_gt_eq_combine_hi_lo(p5, p4);
  const BenchGtEq p32g = bench_gt_eq_combine_hi_lo(p3, p2);
  const BenchGtEq p10 = bench_gt_eq_combine_hi_lo(p1, p0);
  const BenchGtEq p7654 = bench_gt_eq_combine_hi_lo(p76, p54);
  const BenchGtEq p3210 = bench_gt_eq_combine_hi_lo(p32g, p10);
  const BenchGtEq pall = bench_gt_eq_combine_hi_lo(p7654, p3210);

  const BenchGtEq p276 = bench_gt_eq_combine_hi_lo(p27, p26);
  const BenchGtEq p254 = bench_gt_eq_combine_hi_lo(p25, p24);
  const BenchGtEq p232 = bench_gt_eq_combine_hi_lo(p23, p22);
  const BenchGtEq p210 = bench_gt_eq_combine_hi_lo(p21, p20);
  const BenchGtEq p27654 = bench_gt_eq_combine_hi_lo(p276, p254);
  const BenchGtEq p23210 = bench_gt_eq_combine_hi_lo(p232, p210);
  const BenchGtEq p2all = bench_gt_eq_combine_hi_lo(p27654, p23210);

  const BenchGtEq p376 = bench_gt_eq_combine_hi_lo(p37, p36);
  const BenchGtEq p354 = bench_gt_eq_combine_hi_lo(p35, p34);
  const BenchGtEq p332 = bench_gt_eq_combine_hi_lo(p33, p32);
  const BenchGtEq p310 = bench_gt_eq_combine_hi_lo(p31, p30);
  const BenchGtEq p37654 = bench_gt_eq_combine_hi_lo(p376, p354);
  const BenchGtEq p33210 = bench_gt_eq_combine_hi_lo(p332, p310);
  const BenchGtEq p3all = bench_gt_eq_combine_hi_lo(p37654, p33210);

#if BLST_AVX_NATIVE
  return BenchThreeWayBuckets{
      static_cast<u8>(_kor_mask8(static_cast<__mmask8>(pall.gt), static_cast<__mmask8>(pall.eq))),
      static_cast<u8>(_kor_mask8(static_cast<__mmask8>(p2all.gt), static_cast<__mmask8>(p2all.eq))),
      static_cast<u8>(_kor_mask8(static_cast<__mmask8>(p3all.gt), static_cast<__mmask8>(p3all.eq))),
  };
#else
  return BenchThreeWayBuckets{
      static_cast<u8>(pall.gt | pall.eq),
      static_cast<u8>(p2all.gt | p2all.eq),
      static_cast<u8>(p3all.gt | p3all.eq),
  };
#endif
}

template <u64 Scale>
BLST_FORCEINLINE void bench_high_add_mod_masked(cc512i* limbs, const u8 lane_mask)
{
  const auto& mod_vec = bench_mod_vec<Scale>();
  for (std::size_t i = 0; i < NWORDS; ++i) {
    limbs[NWORDS + i] = VMADD(
        limbs[NWORDS + i], lane_mask, limbs[NWORDS + i], mod_vec[i]);
  }
}

template <u64 Scale>
BLST_FORCEINLINE void bench_high_sub_mod_masked(cc512i* limbs, const u8 lane_mask)
{
  const auto& mod_vec = bench_mod_vec<Scale>();
  for (std::size_t i = 0; i < NWORDS; ++i) {
    limbs[NWORDS + i] = VMSUB(
        limbs[NWORDS + i], lane_mask, limbs[NWORDS + i], mod_vec[i]);
  }
}

struct Fp8x1LUFpBenchData
{
  explicit Fp8x1LUFpBenchData(const std::size_t n_packs)
      : n_packs_(n_packs)
      , lazy_a(n_packs_)
      , lazy_b(n_packs_)
      , lazy_in(n_packs_)
      , lazy_in_hi_neg4_pos4(n_packs_)
      , lazy_out(n_packs_)
  {
    std::mt19937_64 gen(0x6C7546708F1A2D3EULL + static_cast<u64>(n_packs_));

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
      UFpPack8x1 wide_a_i{};
      UFpPack8x1 wide_b_i{};
      UFpPack8x1 wide_c_i{};
      UFpPack8x1 wide_d_i{};
      UFpPack8x1 wide_e_i{};
      UFpPack8x1 wide_f_i{};
      UFpPack8x1 wide_g_i{};
      FpPack8x1::mul(wide_a_i, pack_a, pack_b);
      FpPack8x1::mul(wide_b_i, pack_c, pack_d);
      FpPack8x1::mul(wide_c_i, pack_a, pack_c);
      FpPack8x1::mul(wide_d_i, pack_b, pack_d);
      FpPack8x1::mul(wide_e_i, pack_a, pack_d);
      FpPack8x1::mul(wide_f_i, pack_b, pack_c);
      FpPack8x1::mul(wide_g_i, pack_a, pack_a);

      bench_copy_as_lufp(lazy_a[i], wide_a_i);
      bench_copy_as_lufp(lazy_b[i], wide_b_i);

      LUFpPack8x1::sub(lazy_in[i], wide_a_i, wide_b_i);
      bench_make_lufp_hi_neg4_pos4(
          lazy_in_hi_neg4_pos4[i],
          wide_a_i,
          wide_b_i,
          wide_c_i,
          wide_d_i,
          wide_e_i,
          wide_f_i,
          wide_g_i,
          (i & 1u) == 0u);
    }
  }

  const std::size_t n_packs_;
  std::vector<LUFpPack8x1> lazy_a;
  std::vector<LUFpPack8x1> lazy_b;
  std::vector<LUFpPack8x1> lazy_in;
  std::vector<LUFpPack8x1> lazy_in_hi_neg4_pos4;
  std::vector<LUFpPack8x1> lazy_out;
};

struct Fp8x1UFpBenchData
{
  explicit Fp8x1UFpBenchData(const std::size_t n_packs)
      : n_packs_(n_packs)
      , wide_a(n_packs_)
      , wide_b(n_packs_)
      , wide_out(n_packs_)
  {
    std::mt19937_64 gen(0x5546708F1A2D3E6CULL + static_cast<u64>(n_packs_));

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
    }
  }

  const std::size_t n_packs_;
  std::vector<UFpPack8x1> wide_a;
  std::vector<UFpPack8x1> wide_b;
  std::vector<UFpPack8x1> wide_out;
};

struct Fp8x1WidePairBenchData
{
  explicit Fp8x1WidePairBenchData(const std::size_t n_packs)
      : n_packs_(n_packs)
      , wide_a(n_packs_)
      , wide_b(n_packs_)
      , lazy_a(n_packs_)
      , lazy_b(n_packs_)
      , wide_out(n_packs_)
      , lazy_out(n_packs_)
  {
    std::mt19937_64 gen(0x8D13D4375BA14E29ULL + static_cast<u64>(n_packs_));

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
      bench_copy_as_lufp(lazy_a[i], wide_a[i]);
      bench_copy_as_lufp(lazy_b[i], wide_b[i]);
    }
  }

  const std::size_t n_packs_;
  std::vector<UFpPack8x1> wide_a;
  std::vector<UFpPack8x1> wide_b;
  std::vector<LUFpPack8x1> lazy_a;
  std::vector<LUFpPack8x1> lazy_b;
  std::vector<UFpPack8x1> wide_out;
  std::vector<LUFpPack8x1> lazy_out;
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
      FpPack8x1::mul(wide[c], ap, bp);
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
    FpPack8x1::mul(wide0, a0, b0);
    FpPack8x1::mul(wide1, a1, b1);

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

using MulV2K2Fn = void (*)(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    std::size_t n_blocks);

BLST_BENCH_NOINLINE void mul_v2_k2_cpp_baseline_once(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  fp8x1_mul_k<2>(out, a, b, n_blocks);
}

BLST_BENCH_NOINLINE void mul_v2_k2_stream_low_high_once(
    UFpPack8x1& r,
    const FpPack8x1& a,
    const FpPack8x1& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];

  cc512i z0  = VZERO(), z1  = VZERO(), z2  = VZERO(), z3  = VZERO();
  cc512i z4  = VZERO(), z5  = VZERO(), z6  = VZERO(), z7  = VZERO();
  cc512i z8  = VZERO(), z9  = VZERO(), z10 = VZERO(), z11 = VZERO();
  cc512i z12 = VZERO(), z13 = VZERO(), z14 = VZERO(), z15 = VZERO();
  cc512i y0  = VZERO(), y1  = VZERO(), y2  = VZERO(), y3  = VZERO();
  cc512i y4  = VZERO(), y5  = VZERO(), y6  = VZERO();
  cc512i y8  = VZERO(), y9  = VZERO(), y10 = VZERO(), y11 = VZERO();
  cc512i y12 = VZERO(), y13 = VZERO(), y14 = VZERO();

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

  z8 = VMACLO(z8, a4, b4);
  y8 = VMACHI(y8, a4, b4);
  y8 = VSHL(y8, BALIGN);

  z9 = VMACLO(y8, a4, b5); z9 = VMACLO(z9, a5, b4);
  y9 = VMACHI(y9, a4, b5); y9 = VMACHI(y9, a5, b4);
  y9 = VSHL(y9, BALIGN);

  z10 = VMACLO(y9, a4, b6); z10 = VMACLO(z10, a5, b5); z10 = VMACLO(z10, a6, b4);
  y10 = VMACHI(y10, a4, b6); y10 = VMACHI(y10, a5, b5); y10 = VMACHI(y10, a6, b4);
  y10 = VSHL(y10, BALIGN);

  z11 = VMACLO(y10, a4, b7); z11 = VMACLO(z11, a5, b6);
  z11 = VMACLO(z11, a6, b5); z11 = VMACLO(z11, a7, b4);
  y11 = VMACHI(y11, a4, b7); y11 = VMACHI(y11, a5, b6);
  y11 = VMACHI(y11, a6, b5); y11 = VMACHI(y11, a7, b4);
  y11 = VSHL(y11, BALIGN);

  z12 = VMACLO(y11, a5, b7); z12 = VMACLO(z12, a6, b6); z12 = VMACLO(z12, a7, b5);
  y12 = VMACHI(y12, a5, b7); y12 = VMACHI(y12, a6, b6); y12 = VMACHI(y12, a7, b5);
  y12 = VSHL(y12, BALIGN);

  z13 = VMACLO(y12, a6, b7); z13 = VMACLO(z13, a7, b6);
  y13 = VMACHI(y13, a6, b7); y13 = VMACHI(y13, a7, b6);
  y13 = VSHL(y13, BALIGN);

  z14 = VMACLO(y13, a7, b7);
  y14 = VMACHI(y14, a7, b7);
  y14 = VSHL(y14, BALIGN);
  z15 = y14;

  r[0] = z0; r[1] = z1; r[2] = z2; r[3] = z3;
  r[4] = z4; r[5] = z5; r[6] = z6; r[7] = z7;
  r[8] = z8; r[9] = z9; r[10] = z10; r[11] = z11;
  r[12] = z12; r[13] = z13; r[14] = z14; r[15] = z15;
}

BLST_BENCH_NOINLINE void mul_v2_k2_stream_middle_once(
    UFpPack8x1& r,
    const FpPack8x1& a,
    const FpPack8x1& b)
{
  const cc512i ta0 = VADD(a[0], a[4]), ta1 = VADD(a[1], a[5]);
  const cc512i ta2 = VADD(a[2], a[6]), ta3 = VADD(a[3], a[7]);
  const cc512i tb0 = VADD(b[0], b[4]), tb1 = VADD(b[1], b[5]);
  const cc512i tb2 = VADD(b[2], b[6]), tb3 = VADD(b[3], b[7]);

  const cc512i z4 = r[4], z5 = r[5], z6 = r[6], z7 = r[7];

  cc512i y0 = VZERO(), y1 = VZERO(), y2 = VZERO(), y3 = VZERO();
  cc512i y4 = VZERO(), y5 = VZERO(), y6 = VZERO();
  cc512i m0 = VZERO(), m1 = VZERO(), m2 = VZERO(), m3 = VZERO();
  cc512i m4 = VZERO(), m5 = VZERO(), m6 = VZERO(), m7 = VZERO();

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

  m0 = VSUB(m0, VADD(r[0], r[8]));
  m1 = VSUB(m1, VADD(r[1], r[9]));
  m2 = VSUB(m2, VADD(r[2], r[10]));
  m3 = VSUB(m3, VADD(r[3], r[11]));
  m4 = VSUB(m4, VADD(z4, r[12]));
  m5 = VSUB(m5, VADD(z5, r[13]));
  m6 = VSUB(m6, VADD(z6, r[14]));
  m7 = VSUB(m7, VADD(z7, r[15]));

  r[4] = VADD(r[4], m0); r[5] = VADD(r[5], m1);
  r[6] = VADD(r[6], m2); r[7] = VADD(r[7], m3);
  r[8] = VADD(r[8], m4); r[9] = VADD(r[9], m5);
  r[10] = VADD(r[10], m6); r[11] = VADD(r[11], m7);
}

BLST_BENCH_NOINLINE void mul_v2_k2_cpp_candidate_v6b_once(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks);
extern "C" BLST_BENCH_NOINLINE void blst_fp8x1_mul_v2_k2_v6b_cpp_entry(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    std::size_t n_blocks);

BLST_BENCH_NOINLINE void mul_v2_k2_cpp_candidate_v2_once(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  for (std::size_t i = 0; i < n_blocks; ++i) {
    // Stage-ready schedule: each stream advances once prerequisites are complete.
    mul_v2_k2_stream_low_high_once(out[0], a[0], b[0]);
    mul_v2_k2_stream_low_high_once(out[1], a[1], b[1]);
    mul_v2_k2_stream_middle_once(out[0], a[0], b[0]);
    mul_v2_k2_stream_middle_once(out[1], a[1], b[1]);
    out += 2;
    a += 2;
    b += 2;
  }
}

BLST_BENCH_NOINLINE void mul_v2_k2_cpp_candidate_once(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  // Promote current best candidate: v6b pair-carry schoolbook.
  mul_v2_k2_cpp_candidate_v6b_once(out, a, b, n_blocks);
}

BLST_BENCH_NOINLINE void mul_v2_k2_stream_full_once(
    UFpPack8x1& r,
    const FpPack8x1& a,
    const FpPack8x1& b)
{
  {
    const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
    const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
    const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
    const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];

    cc512i z0  = VZERO(), z1  = VZERO(), z2  = VZERO(), z3  = VZERO();
    cc512i z4  = VZERO(), z5  = VZERO(), z6  = VZERO(), z7  = VZERO();
    cc512i z8  = VZERO(), z9  = VZERO(), z10 = VZERO(), z11 = VZERO();
    cc512i z12 = VZERO(), z13 = VZERO(), z14 = VZERO(), z15 = VZERO();
    cc512i y0  = VZERO(), y1  = VZERO(), y2  = VZERO(), y3  = VZERO();
    cc512i y4  = VZERO(), y5  = VZERO(), y6  = VZERO();
    cc512i y8  = VZERO(), y9  = VZERO(), y10 = VZERO(), y11 = VZERO();
    cc512i y12 = VZERO(), y13 = VZERO(), y14 = VZERO();

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

    z8 = VMACLO(z8, a4, b4);
    y8 = VMACHI(y8, a4, b4);
    y8 = VSHL(y8, BALIGN);

    z9 = VMACLO(y8, a4, b5); z9 = VMACLO(z9, a5, b4);
    y9 = VMACHI(y9, a4, b5); y9 = VMACHI(y9, a5, b4);
    y9 = VSHL(y9, BALIGN);

    z10 = VMACLO(y9, a4, b6); z10 = VMACLO(z10, a5, b5); z10 = VMACLO(z10, a6, b4);
    y10 = VMACHI(y10, a4, b6); y10 = VMACHI(y10, a5, b5); y10 = VMACHI(y10, a6, b4);
    y10 = VSHL(y10, BALIGN);

    z11 = VMACLO(y10, a4, b7); z11 = VMACLO(z11, a5, b6);
    z11 = VMACLO(z11, a6, b5); z11 = VMACLO(z11, a7, b4);
    y11 = VMACHI(y11, a4, b7); y11 = VMACHI(y11, a5, b6);
    y11 = VMACHI(y11, a6, b5); y11 = VMACHI(y11, a7, b4);
    y11 = VSHL(y11, BALIGN);

    z12 = VMACLO(y11, a5, b7); z12 = VMACLO(z12, a6, b6); z12 = VMACLO(z12, a7, b5);
    y12 = VMACHI(y12, a5, b7); y12 = VMACHI(y12, a6, b6); y12 = VMACHI(y12, a7, b5);
    y12 = VSHL(y12, BALIGN);

    z13 = VMACLO(y12, a6, b7); z13 = VMACLO(z13, a7, b6);
    y13 = VMACHI(y13, a6, b7); y13 = VMACHI(y13, a7, b6);
    y13 = VSHL(y13, BALIGN);

    z14 = VMACLO(y13, a7, b7);
    y14 = VMACHI(y14, a7, b7);
    y14 = VSHL(y14, BALIGN);
    z15 = y14;

    r[0] = z0; r[1] = z1; r[2] = z2; r[3] = z3;
    r[4] = z4; r[5] = z5; r[6] = z6; r[7] = z7;
    r[8] = z8; r[9] = z9; r[10] = z10; r[11] = z11;
    r[12] = z12; r[13] = z13; r[14] = z14; r[15] = z15;
  }

  {
    const cc512i ta0 = VADD(a[0], a[4]), ta1 = VADD(a[1], a[5]);
    const cc512i ta2 = VADD(a[2], a[6]), ta3 = VADD(a[3], a[7]);
    const cc512i tb0 = VADD(b[0], b[4]), tb1 = VADD(b[1], b[5]);
    const cc512i tb2 = VADD(b[2], b[6]), tb3 = VADD(b[3], b[7]);

    const cc512i z4 = r[4], z5 = r[5], z6 = r[6], z7 = r[7];

    cc512i y0 = VZERO(), y1 = VZERO(), y2 = VZERO(), y3 = VZERO();
    cc512i y4 = VZERO(), y5 = VZERO(), y6 = VZERO();
    cc512i m0 = VZERO(), m1 = VZERO(), m2 = VZERO(), m3 = VZERO();
    cc512i m4 = VZERO(), m5 = VZERO(), m6 = VZERO(), m7 = VZERO();

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

    m0 = VSUB(m0, VADD(r[0], r[8]));
    m1 = VSUB(m1, VADD(r[1], r[9]));
    m2 = VSUB(m2, VADD(r[2], r[10]));
    m3 = VSUB(m3, VADD(r[3], r[11]));
    m4 = VSUB(m4, VADD(z4, r[12]));
    m5 = VSUB(m5, VADD(z5, r[13]));
    m6 = VSUB(m6, VADD(z6, r[14]));
    m7 = VSUB(m7, VADD(z7, r[15]));

    r[4] = VADD(r[4], m0); r[5] = VADD(r[5], m1);
    r[6] = VADD(r[6], m2); r[7] = VADD(r[7], m3);
    r[8] = VADD(r[8], m4); r[9] = VADD(r[9], m5);
    r[10] = VADD(r[10], m6); r[11] = VADD(r[11], m7);
  }
}

BLST_BENCH_NOINLINE void mul_v2_k2_cpp_candidate_v3_once(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  for (std::size_t i = 0; i < n_blocks; ++i) {
    // v3: per-stream full-stage helper (low/high + middle fused in one call).
    mul_v2_k2_stream_full_once(out[0], a[0], b[0]);
    mul_v2_k2_stream_full_once(out[1], a[1], b[1]);
    out += 2;
    a += 2;
    b += 2;
  }
}

BLST_BENCH_NOINLINE void mul_v2_k2_cpp_candidate_v4_once(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  if (n_blocks == 0) {
    return;
  }

  UFpPack8x1* out_cur = out;
  const FpPack8x1* a_cur = a;
  const FpPack8x1* b_cur = b;

  // Prime the first block low/high products.
  mul_v2_k2_stream_low_high_once(out_cur[0], a_cur[0], b_cur[0]);
  mul_v2_k2_stream_low_high_once(out_cur[1], a_cur[1], b_cur[1]);

  // Wavefront schedule:
  // next-block low/high for each stream, then current-block middle for that stream.
  for (std::size_t block = 0; (block + 1) < n_blocks; ++block) {
    UFpPack8x1* out_next = out_cur + 2;
    const FpPack8x1* a_next = a_cur + 2;
    const FpPack8x1* b_next = b_cur + 2;

    mul_v2_k2_stream_low_high_once(out_next[0], a_next[0], b_next[0]);
    mul_v2_k2_stream_middle_once(out_cur[0], a_cur[0], b_cur[0]);

    mul_v2_k2_stream_low_high_once(out_next[1], a_next[1], b_next[1]);
    mul_v2_k2_stream_middle_once(out_cur[1], a_cur[1], b_cur[1]);

    out_cur = out_next;
    a_cur = a_next;
    b_cur = b_next;
  }

  // Drain: finalize middle stage for the last primed block.
  mul_v2_k2_stream_middle_once(out_cur[0], a_cur[0], b_cur[0]);
  mul_v2_k2_stream_middle_once(out_cur[1], a_cur[1], b_cur[1]);
}

BLST_BENCH_NOINLINE void mul_v2_k2_cpp_candidate_v5_once(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  for (std::size_t i = 0; i < n_blocks; ++i) {
    // Hybrid lockstep schedule: keep split Karatsuba stages, but run stream-local
    // low/high -> middle before switching to the second stream.
    mul_v2_k2_stream_low_high_once(out[0], a[0], b[0]);
    mul_v2_k2_stream_middle_once(out[0], a[0], b[0]);
    mul_v2_k2_stream_low_high_once(out[1], a[1], b[1]);
    mul_v2_k2_stream_middle_once(out[1], a[1], b[1]);
    out += 2;
    a += 2;
    b += 2;
  }
}

BLST_FORCEINLINE void mul_v2_k2_schoolbook_stream_once(
    UFpPack8x1& r,
    const FpPack8x1& a,
    const FpPack8x1& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];

  cc512i carry = VZERO();
  cc512i lo = VZERO();
  cc512i hi = VZERO();

#define BLST_SB_COL_BEGIN() \
  lo = carry;               \
  hi = VZERO()
#define BLST_SB_TERM(aa, bb)      \
  lo = VMACLO(lo, (aa), (bb));    \
  hi = VMACHI(hi, (aa), (bb))
#define BLST_SB_COL_END(idx) \
  r[(idx)] = lo;             \
  carry = VSHL(hi, BALIGN)

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a0, b0);
  BLST_SB_COL_END(0);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a0, b1); BLST_SB_TERM(a1, b0);
  BLST_SB_COL_END(1);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a0, b2); BLST_SB_TERM(a1, b1); BLST_SB_TERM(a2, b0);
  BLST_SB_COL_END(2);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a0, b3); BLST_SB_TERM(a1, b2); BLST_SB_TERM(a2, b1); BLST_SB_TERM(a3, b0);
  BLST_SB_COL_END(3);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a0, b4); BLST_SB_TERM(a1, b3); BLST_SB_TERM(a2, b2); BLST_SB_TERM(a3, b1);
  BLST_SB_TERM(a4, b0);
  BLST_SB_COL_END(4);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a0, b5); BLST_SB_TERM(a1, b4); BLST_SB_TERM(a2, b3); BLST_SB_TERM(a3, b2);
  BLST_SB_TERM(a4, b1); BLST_SB_TERM(a5, b0);
  BLST_SB_COL_END(5);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a0, b6); BLST_SB_TERM(a1, b5); BLST_SB_TERM(a2, b4); BLST_SB_TERM(a3, b3);
  BLST_SB_TERM(a4, b2); BLST_SB_TERM(a5, b1); BLST_SB_TERM(a6, b0);
  BLST_SB_COL_END(6);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a0, b7); BLST_SB_TERM(a1, b6); BLST_SB_TERM(a2, b5); BLST_SB_TERM(a3, b4);
  BLST_SB_TERM(a4, b3); BLST_SB_TERM(a5, b2); BLST_SB_TERM(a6, b1); BLST_SB_TERM(a7, b0);
  BLST_SB_COL_END(7);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a1, b7); BLST_SB_TERM(a2, b6); BLST_SB_TERM(a3, b5); BLST_SB_TERM(a4, b4);
  BLST_SB_TERM(a5, b3); BLST_SB_TERM(a6, b2); BLST_SB_TERM(a7, b1);
  BLST_SB_COL_END(8);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a2, b7); BLST_SB_TERM(a3, b6); BLST_SB_TERM(a4, b5); BLST_SB_TERM(a5, b4);
  BLST_SB_TERM(a6, b3); BLST_SB_TERM(a7, b2);
  BLST_SB_COL_END(9);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a3, b7); BLST_SB_TERM(a4, b6); BLST_SB_TERM(a5, b5); BLST_SB_TERM(a6, b4);
  BLST_SB_TERM(a7, b3);
  BLST_SB_COL_END(10);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a4, b7); BLST_SB_TERM(a5, b6); BLST_SB_TERM(a6, b5); BLST_SB_TERM(a7, b4);
  BLST_SB_COL_END(11);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a5, b7); BLST_SB_TERM(a6, b6); BLST_SB_TERM(a7, b5);
  BLST_SB_COL_END(12);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a6, b7); BLST_SB_TERM(a7, b6);
  BLST_SB_COL_END(13);

  BLST_SB_COL_BEGIN();
  BLST_SB_TERM(a7, b7);
  BLST_SB_COL_END(14);

  r[15] = carry;

#undef BLST_SB_COL_BEGIN
#undef BLST_SB_TERM
#undef BLST_SB_COL_END
}

BLST_BENCH_NOINLINE void mul_v2_k2_cpp_candidate_v6_once(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  for (std::size_t i = 0; i < n_blocks; ++i) {
    // v6: full schoolbook (more IFMA, more regular dependency graph) for A/B testing.
    mul_v2_k2_schoolbook_stream_once(out[0], a[0], b[0]);
    mul_v2_k2_schoolbook_stream_once(out[1], a[1], b[1]);
    out += 2;
    a += 2;
    b += 2;
  }
}

BLST_FORCEINLINE void mul_v2_k2_schoolbook_stream_paircarry_once(
    UFpPack8x1& r,
    const FpPack8x1& a,
    const FpPack8x1& b)
{
  const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  const cc512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  const cc512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];

  cc512i carry = VZERO();
  cc512i lo0 = VZERO(), hi0 = VZERO();
  cc512i lo1 = VZERO(), hi1 = VZERO();

#define BLST_SB_PC_PAIR_BEGIN() \
  lo0 = carry;                  \
  hi0 = VZERO()
#define BLST_SB_PC_TERM0(aa, bb)    \
  lo0 = VMACLO(lo0, (aa), (bb));    \
  hi0 = VMACHI(hi0, (aa), (bb))
#define BLST_SB_PC_PAIR_MID()  \
  lo1 = VSHL(hi0, BALIGN);     \
  hi1 = VZERO()
#define BLST_SB_PC_TERM1(aa, bb)    \
  lo1 = VMACLO(lo1, (aa), (bb));    \
  hi1 = VMACHI(hi1, (aa), (bb))
#define BLST_SB_PC_PAIR_END(idx0, idx1) \
  r[(idx0)] = lo0;                     \
  r[(idx1)] = lo1;                     \
  carry = VSHL(hi1, BALIGN)

  BLST_SB_PC_PAIR_BEGIN();
  BLST_SB_PC_TERM0(a0, b0);
  BLST_SB_PC_PAIR_MID();
  BLST_SB_PC_TERM1(a0, b1); BLST_SB_PC_TERM1(a1, b0);
  BLST_SB_PC_PAIR_END(0, 1);

  BLST_SB_PC_PAIR_BEGIN();
  BLST_SB_PC_TERM0(a0, b2); BLST_SB_PC_TERM0(a1, b1); BLST_SB_PC_TERM0(a2, b0);
  BLST_SB_PC_PAIR_MID();
  BLST_SB_PC_TERM1(a0, b3); BLST_SB_PC_TERM1(a1, b2); BLST_SB_PC_TERM1(a2, b1); BLST_SB_PC_TERM1(a3, b0);
  BLST_SB_PC_PAIR_END(2, 3);

  BLST_SB_PC_PAIR_BEGIN();
  BLST_SB_PC_TERM0(a0, b4); BLST_SB_PC_TERM0(a1, b3); BLST_SB_PC_TERM0(a2, b2); BLST_SB_PC_TERM0(a3, b1);
  BLST_SB_PC_TERM0(a4, b0);
  BLST_SB_PC_PAIR_MID();
  BLST_SB_PC_TERM1(a0, b5); BLST_SB_PC_TERM1(a1, b4); BLST_SB_PC_TERM1(a2, b3); BLST_SB_PC_TERM1(a3, b2);
  BLST_SB_PC_TERM1(a4, b1); BLST_SB_PC_TERM1(a5, b0);
  BLST_SB_PC_PAIR_END(4, 5);

  BLST_SB_PC_PAIR_BEGIN();
  BLST_SB_PC_TERM0(a0, b6); BLST_SB_PC_TERM0(a1, b5); BLST_SB_PC_TERM0(a2, b4); BLST_SB_PC_TERM0(a3, b3);
  BLST_SB_PC_TERM0(a4, b2); BLST_SB_PC_TERM0(a5, b1); BLST_SB_PC_TERM0(a6, b0);
  BLST_SB_PC_PAIR_MID();
  BLST_SB_PC_TERM1(a0, b7); BLST_SB_PC_TERM1(a1, b6); BLST_SB_PC_TERM1(a2, b5); BLST_SB_PC_TERM1(a3, b4);
  BLST_SB_PC_TERM1(a4, b3); BLST_SB_PC_TERM1(a5, b2); BLST_SB_PC_TERM1(a6, b1); BLST_SB_PC_TERM1(a7, b0);
  BLST_SB_PC_PAIR_END(6, 7);

  BLST_SB_PC_PAIR_BEGIN();
  BLST_SB_PC_TERM0(a1, b7); BLST_SB_PC_TERM0(a2, b6); BLST_SB_PC_TERM0(a3, b5); BLST_SB_PC_TERM0(a4, b4);
  BLST_SB_PC_TERM0(a5, b3); BLST_SB_PC_TERM0(a6, b2); BLST_SB_PC_TERM0(a7, b1);
  BLST_SB_PC_PAIR_MID();
  BLST_SB_PC_TERM1(a2, b7); BLST_SB_PC_TERM1(a3, b6); BLST_SB_PC_TERM1(a4, b5); BLST_SB_PC_TERM1(a5, b4);
  BLST_SB_PC_TERM1(a6, b3); BLST_SB_PC_TERM1(a7, b2);
  BLST_SB_PC_PAIR_END(8, 9);

  BLST_SB_PC_PAIR_BEGIN();
  BLST_SB_PC_TERM0(a3, b7); BLST_SB_PC_TERM0(a4, b6); BLST_SB_PC_TERM0(a5, b5); BLST_SB_PC_TERM0(a6, b4);
  BLST_SB_PC_TERM0(a7, b3);
  BLST_SB_PC_PAIR_MID();
  BLST_SB_PC_TERM1(a4, b7); BLST_SB_PC_TERM1(a5, b6); BLST_SB_PC_TERM1(a6, b5); BLST_SB_PC_TERM1(a7, b4);
  BLST_SB_PC_PAIR_END(10, 11);

  BLST_SB_PC_PAIR_BEGIN();
  BLST_SB_PC_TERM0(a5, b7); BLST_SB_PC_TERM0(a6, b6); BLST_SB_PC_TERM0(a7, b5);
  BLST_SB_PC_PAIR_MID();
  BLST_SB_PC_TERM1(a6, b7); BLST_SB_PC_TERM1(a7, b6);
  BLST_SB_PC_PAIR_END(12, 13);

  // Final pair: column 14 has one term; column 15 is carry only.
  BLST_SB_PC_PAIR_BEGIN();
  BLST_SB_PC_TERM0(a7, b7);
  r[14] = lo0;
  r[15] = VSHL(hi0, BALIGN);

#undef BLST_SB_PC_PAIR_BEGIN
#undef BLST_SB_PC_TERM0
#undef BLST_SB_PC_PAIR_MID
#undef BLST_SB_PC_TERM1
#undef BLST_SB_PC_PAIR_END
}

BLST_BENCH_NOINLINE void mul_v2_k2_cpp_candidate_v6b_once(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  for (std::size_t i = 0; i < n_blocks; ++i) {
    // v6b: full-schoolbook with 2-column carry cadence (pair-carry).
    mul_v2_k2_schoolbook_stream_paircarry_once(out[0], a[0], b[0]);
    mul_v2_k2_schoolbook_stream_paircarry_once(out[1], a[1], b[1]);
    out += 2;
    a += 2;
    b += 2;
  }
}

BLST_BENCH_NOINLINE void mul_v2_k1_cpp_candidate_v6b_once(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  for (std::size_t i = 0; i < n_blocks; ++i) {
    // k=1 benchmark-only candidate: same pair-carry schoolbook inner style as v6b.
    mul_v2_k2_schoolbook_stream_paircarry_once(out[0], a[0], b[0]);
    ++out;
    ++a;
    ++b;
  }
}

extern "C" BLST_BENCH_NOINLINE void blst_fp8x1_mul_v2_k2_v6b_cpp_entry(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  mul_v2_k2_cpp_candidate_v6b_once(out, a, b, n_blocks);
}

BLST_FORCEINLINE void mul_v2_k2_schoolbook_dual_interleaved_once(
    UFpPack8x1& r0,
    UFpPack8x1& r1,
    const FpPack8x1& a0,
    const FpPack8x1& b0,
    const FpPack8x1& a1,
    const FpPack8x1& b1)
{
  cc512i carry0 = VZERO(), carry1 = VZERO();
  cc512i lo0 = VZERO(), hi0 = VZERO();
  cc512i lo1 = VZERO(), hi1 = VZERO();

#define BLST_SB2_COL_BEGIN() \
  lo0 = carry0;              \
  hi0 = VZERO();             \
  lo1 = carry1;              \
  hi1 = VZERO()
#define BLST_SB2_TERM(i, j)             \
  lo0 = VMACLO(lo0, a0[(i)], b0[(j)]);  \
  hi0 = VMACHI(hi0, a0[(i)], b0[(j)]);  \
  lo1 = VMACLO(lo1, a1[(i)], b1[(j)]);  \
  hi1 = VMACHI(hi1, a1[(i)], b1[(j)])
#define BLST_SB2_COL_END(idx)    \
  r0[(idx)] = lo0;               \
  carry0 = VSHL(hi0, BALIGN);    \
  r1[(idx)] = lo1;               \
  carry1 = VSHL(hi1, BALIGN)

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(0, 0);
  BLST_SB2_COL_END(0);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(0, 1); BLST_SB2_TERM(1, 0);
  BLST_SB2_COL_END(1);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(0, 2); BLST_SB2_TERM(1, 1); BLST_SB2_TERM(2, 0);
  BLST_SB2_COL_END(2);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(0, 3); BLST_SB2_TERM(1, 2); BLST_SB2_TERM(2, 1); BLST_SB2_TERM(3, 0);
  BLST_SB2_COL_END(3);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(0, 4); BLST_SB2_TERM(1, 3); BLST_SB2_TERM(2, 2); BLST_SB2_TERM(3, 1);
  BLST_SB2_TERM(4, 0);
  BLST_SB2_COL_END(4);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(0, 5); BLST_SB2_TERM(1, 4); BLST_SB2_TERM(2, 3); BLST_SB2_TERM(3, 2);
  BLST_SB2_TERM(4, 1); BLST_SB2_TERM(5, 0);
  BLST_SB2_COL_END(5);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(0, 6); BLST_SB2_TERM(1, 5); BLST_SB2_TERM(2, 4); BLST_SB2_TERM(3, 3);
  BLST_SB2_TERM(4, 2); BLST_SB2_TERM(5, 1); BLST_SB2_TERM(6, 0);
  BLST_SB2_COL_END(6);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(0, 7); BLST_SB2_TERM(1, 6); BLST_SB2_TERM(2, 5); BLST_SB2_TERM(3, 4);
  BLST_SB2_TERM(4, 3); BLST_SB2_TERM(5, 2); BLST_SB2_TERM(6, 1); BLST_SB2_TERM(7, 0);
  BLST_SB2_COL_END(7);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(1, 7); BLST_SB2_TERM(2, 6); BLST_SB2_TERM(3, 5); BLST_SB2_TERM(4, 4);
  BLST_SB2_TERM(5, 3); BLST_SB2_TERM(6, 2); BLST_SB2_TERM(7, 1);
  BLST_SB2_COL_END(8);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(2, 7); BLST_SB2_TERM(3, 6); BLST_SB2_TERM(4, 5); BLST_SB2_TERM(5, 4);
  BLST_SB2_TERM(6, 3); BLST_SB2_TERM(7, 2);
  BLST_SB2_COL_END(9);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(3, 7); BLST_SB2_TERM(4, 6); BLST_SB2_TERM(5, 5); BLST_SB2_TERM(6, 4);
  BLST_SB2_TERM(7, 3);
  BLST_SB2_COL_END(10);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(4, 7); BLST_SB2_TERM(5, 6); BLST_SB2_TERM(6, 5); BLST_SB2_TERM(7, 4);
  BLST_SB2_COL_END(11);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(5, 7); BLST_SB2_TERM(6, 6); BLST_SB2_TERM(7, 5);
  BLST_SB2_COL_END(12);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(6, 7); BLST_SB2_TERM(7, 6);
  BLST_SB2_COL_END(13);

  BLST_SB2_COL_BEGIN();
  BLST_SB2_TERM(7, 7);
  BLST_SB2_COL_END(14);

  r0[15] = carry0;
  r1[15] = carry1;

#undef BLST_SB2_COL_BEGIN
#undef BLST_SB2_TERM
#undef BLST_SB2_COL_END
}

BLST_BENCH_NOINLINE void mul_v2_k2_cpp_candidate_v6a_once(
    UFpPack8x1* BLST_restrict out,
    const FpPack8x1* BLST_restrict a,
    const FpPack8x1* BLST_restrict b,
    const std::size_t n_blocks)
{
  for (std::size_t i = 0; i < n_blocks; ++i) {
    // v6a: full-schoolbook with dual-stream diagonal interleave.
    mul_v2_k2_schoolbook_dual_interleaved_once(out[0], out[1], a[0], b[0], a[1], b[1]);
    out += 2;
    a += 2;
    b += 2;
  }
}

void mul_v2_k2_selfcheck_impl(MulV2K2Fn fn, const char* tag)
{
  constexpr std::size_t kCases = 32;
  constexpr std::size_t kBlocks = 7;
  constexpr std::size_t kPacks = kBlocks * 2;

  std::mt19937_64 gen(0xB61AF1330C5D119FULL);
  std::array<FpPack8x1, kPacks> a{};
  std::array<FpPack8x1, kPacks> b{};
  std::array<UFpPack8x1, kPacks> got{};
  std::array<UFpPack8x1, kPacks> exp{};

  for (std::size_t t = 0; t < kCases; ++t) {
    for (std::size_t i = 0; i < kPacks; ++i) {
      std::array<Fp, 8> lane_a{};
      std::array<Fp, 8> lane_b{};
      for (std::size_t lane = 0; lane < 8; ++lane) {
        random_redc_fp(lane_a[lane], gen);
        random_redc_fp(lane_b[lane], gen);
      }
      FpPack8x1::from_aos(a[i], lane_a.data());
      FpPack8x1::from_aos(b[i], lane_b.data());
    }

    mul_v2_k2_cpp_baseline_once(exp.data(), a.data(), b.data(), kBlocks);
    fn(got.data(), a.data(), b.data(), kBlocks);

    for (std::size_t i = 0; i < kPacks; ++i) {
      FpPack8x1 got_redc{}, exp_redc{};
      UFpPack8x1::redc(got_redc, got[i]);
      UFpPack8x1::redc(exp_redc, exp[i]);

      std::array<Fp, 8> got_lanes{}, exp_lanes{};
      FpPack8x1::to_aos(got_lanes.data(), got_redc);
      FpPack8x1::to_aos(exp_lanes.data(), exp_redc);

      for (std::size_t lane = 0; lane < 8; ++lane) {
        for (std::size_t limb = 0; limb < Fp::num_limb(); ++limb) {
          if (got_lanes[lane][limb] != exp_lanes[lane][limb]) {
            std::fprintf(
                stderr,
                "%s selfcheck failed (case=%zu pack=%zu lane=%zu limb=%zu got=%llx exp=%llx)\n",
                tag,
                t,
                i,
                lane,
                limb,
                static_cast<unsigned long long>(got_lanes[lane][limb]),
                static_cast<unsigned long long>(exp_lanes[lane][limb]));
            std::abort();
          }
        }
      }
    }
  }
}

void mul_v2_k1_selfcheck_impl(MulV2K2Fn fn, const char* tag)
{
  constexpr std::size_t kCases = 32;
  constexpr std::size_t kBlocks = 11;
  constexpr std::size_t kPacks = kBlocks;

  std::mt19937_64 gen(0x91D7B0AF3C2E44D1ULL);
  std::array<FpPack8x1, kPacks> a{};
  std::array<FpPack8x1, kPacks> b{};
  std::array<UFpPack8x1, kPacks> got{};
  std::array<UFpPack8x1, kPacks> exp{};

  for (std::size_t t = 0; t < kCases; ++t) {
    for (std::size_t i = 0; i < kPacks; ++i) {
      std::array<Fp, 8> lane_a{};
      std::array<Fp, 8> lane_b{};
      for (std::size_t lane = 0; lane < 8; ++lane) {
        random_redc_fp(lane_a[lane], gen);
        random_redc_fp(lane_b[lane], gen);
      }
      FpPack8x1::from_aos(a[i], lane_a.data());
      FpPack8x1::from_aos(b[i], lane_b.data());
    }

    fp8x1_mul_k<1>(exp.data(), a.data(), b.data(), kBlocks);
    fn(got.data(), a.data(), b.data(), kBlocks);

    for (std::size_t i = 0; i < kPacks; ++i) {
      FpPack8x1 got_redc{}, exp_redc{};
      UFpPack8x1::redc(got_redc, got[i]);
      UFpPack8x1::redc(exp_redc, exp[i]);

      std::array<Fp, 8> got_lanes{}, exp_lanes{};
      FpPack8x1::to_aos(got_lanes.data(), got_redc);
      FpPack8x1::to_aos(exp_lanes.data(), exp_redc);

      for (std::size_t lane = 0; lane < 8; ++lane) {
        for (std::size_t limb = 0; limb < Fp::num_limb(); ++limb) {
          if (got_lanes[lane][limb] != exp_lanes[lane][limb]) {
            std::fprintf(
                stderr,
                "%s selfcheck failed (case=%zu pack=%zu lane=%zu limb=%zu got=%llx exp=%llx)\n",
                tag,
                t,
                i,
                lane,
                limb,
                static_cast<unsigned long long>(got_lanes[lane][limb]),
                static_cast<unsigned long long>(exp_lanes[lane][limb]));
            std::abort();
          }
        }
      }
    }
  }
}

void mul_v2_k1_v6b_selfcheck_once()
{
  static std::once_flag once;
  std::call_once(once, []() {
    mul_v2_k1_selfcheck_impl(mul_v2_k1_cpp_candidate_v6b_once, "mul_v2_k1_cpp_candidate_v6b");
  });
}

void mul_v2_k2_selfcheck_once()
{
  static std::once_flag once;
  std::call_once(once, []() {
    mul_v2_k2_selfcheck_impl(mul_v2_k2_cpp_candidate_once, "mul_v2_k2_cpp_candidate");
  });
}

void mul_v2_k2_v3_selfcheck_once()
{
  static std::once_flag once;
  std::call_once(once, []() {
    mul_v2_k2_selfcheck_impl(mul_v2_k2_cpp_candidate_v3_once, "mul_v2_k2_cpp_candidate_v3");
  });
}

void mul_v2_k2_v4_selfcheck_once()
{
  static std::once_flag once;
  std::call_once(once, []() {
    mul_v2_k2_selfcheck_impl(mul_v2_k2_cpp_candidate_v4_once, "mul_v2_k2_cpp_candidate_v4");
  });
}

void mul_v2_k2_v5_selfcheck_once()
{
  static std::once_flag once;
  std::call_once(once, []() {
    mul_v2_k2_selfcheck_impl(mul_v2_k2_cpp_candidate_v5_once, "mul_v2_k2_cpp_candidate_v5");
  });
}

void mul_v2_k2_v6_selfcheck_once()
{
  static std::once_flag once;
  std::call_once(once, []() {
    mul_v2_k2_selfcheck_impl(mul_v2_k2_cpp_candidate_v6_once, "mul_v2_k2_cpp_candidate_v6");
  });
}

void mul_v2_k2_v6a_selfcheck_once()
{
  static std::once_flag once;
  std::call_once(once, []() {
    mul_v2_k2_selfcheck_impl(mul_v2_k2_cpp_candidate_v6a_once, "mul_v2_k2_cpp_candidate_v6a");
  });
}

void mul_v2_k2_v6b_selfcheck_once()
{
  static std::once_flag once;
  std::call_once(once, []() {
    mul_v2_k2_selfcheck_impl(mul_v2_k2_cpp_candidate_v6b_once, "mul_v2_k2_cpp_candidate_v6b");
  });
}

#if BLST_HAVE_MUL_V2_K2_V6B_ASM
void mul_v2_k2_v6b_asm_selfcheck_once()
{
  static std::once_flag once;
  std::call_once(once, []() {
    mul_v2_k2_selfcheck_impl(blst_fp8x1_mul_v2_k2_v6b_asm, "blst_fp8x1_mul_v2_k2_v6b_asm");
  });
}
#endif

BLST_FORCEINLINE void bench_consume_wide(std::vector<UFpPack8x1>& wide)
{
  if (wide.empty()) {
    return;
  }
  benchmark::DoNotOptimize(wide.front()[0]);
  benchmark::DoNotOptimize(wide.back()[2 * NWORDS - 1]);
  benchmark::ClobberMemory();
}

BLST_FORCEINLINE void bench_consume_lazy(std::vector<LFpPack8x1>& lazy)
{
  if (lazy.empty()) {
    return;
  }
  benchmark::DoNotOptimize(lazy.front()[0]);
  benchmark::DoNotOptimize(lazy.back()[NWORDS - 1]);
  benchmark::ClobberMemory();
}

BLST_FORCEINLINE void bench_consume_luwide(std::vector<LUFpPack8x1>& lazy)
{
  if (lazy.empty()) {
    return;
  }
  benchmark::DoNotOptimize(lazy.front()[0]);
  benchmark::DoNotOptimize(lazy.back()[2 * NWORDS - 1]);
  benchmark::ClobberMemory();
}

BLST_FORCEINLINE void bench_consume_out(std::vector<FpPack8x1>& out)
{
  if (out.empty()) {
    return;
  }
  benchmark::DoNotOptimize(out.front()[0]);
  benchmark::DoNotOptimize(out.back()[NWORDS - 1]);
  benchmark::ClobberMemory();
}

void BM_Fp8x1MulV2K2_CppCandidate(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * 2;
  Fp8x1KBenchData data(n_packs);
  mul_v2_k2_selfcheck_once();

  for (auto _ : state) {
    mul_v2_k2_cpp_candidate_once(data.wide.data(), data.a.data(), data.b.data(), n_blocks);
    bench_consume_wide(data.wide);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1MulV2K2_CppCandidateV3(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * 2;
  Fp8x1KBenchData data(n_packs);
  mul_v2_k2_v3_selfcheck_once();

  for (auto _ : state) {
    mul_v2_k2_cpp_candidate_v3_once(data.wide.data(), data.a.data(), data.b.data(), n_blocks);
    bench_consume_wide(data.wide);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1MulV2K2_CppCandidateV4(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * 2;
  Fp8x1KBenchData data(n_packs);
  mul_v2_k2_v4_selfcheck_once();

  for (auto _ : state) {
    mul_v2_k2_cpp_candidate_v4_once(data.wide.data(), data.a.data(), data.b.data(), n_blocks);
    bench_consume_wide(data.wide);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1MulV2K2_CppCandidateV5(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * 2;
  Fp8x1KBenchData data(n_packs);
  mul_v2_k2_v5_selfcheck_once();

  for (auto _ : state) {
    mul_v2_k2_cpp_candidate_v5_once(data.wide.data(), data.a.data(), data.b.data(), n_blocks);
    bench_consume_wide(data.wide);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1MulV2K2_CppCandidateV6(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * 2;
  Fp8x1KBenchData data(n_packs);
  mul_v2_k2_v6_selfcheck_once();

  for (auto _ : state) {
    mul_v2_k2_cpp_candidate_v6_once(data.wide.data(), data.a.data(), data.b.data(), n_blocks);
    bench_consume_wide(data.wide);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1MulV2K2_CppCandidateV6A(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * 2;
  Fp8x1KBenchData data(n_packs);
  mul_v2_k2_v6a_selfcheck_once();

  for (auto _ : state) {
    mul_v2_k2_cpp_candidate_v6a_once(data.wide.data(), data.a.data(), data.b.data(), n_blocks);
    bench_consume_wide(data.wide);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1MulV2K2_CppCandidateV6B(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * 2;
  Fp8x1KBenchData data(n_packs);
  mul_v2_k2_v6b_selfcheck_once();

  for (auto _ : state) {
    mul_v2_k2_cpp_candidate_v6b_once(data.wide.data(), data.a.data(), data.b.data(), n_blocks);
    bench_consume_wide(data.wide);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

#if BLST_HAVE_MUL_V2_K2_V6B_ASM
void BM_Fp8x1MulV2K2_CppCandidateV6BAsm(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * 2;
  Fp8x1KBenchData data(n_packs);
  mul_v2_k2_v6b_asm_selfcheck_once();

  for (auto _ : state) {
    blst_fp8x1_mul_v2_k2_v6b_asm(data.wide.data(), data.a.data(), data.b.data(), n_blocks);
    bench_consume_wide(data.wide);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}
#endif

void BM_Fp8x1MulV2K1_CppCandidateV6B(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks;
  Fp8x1KBenchData data(n_packs);
  mul_v2_k1_v6b_selfcheck_once();

  for (auto _ : state) {
    mul_v2_k1_cpp_candidate_v6b_once(data.wide.data(), data.a.data(), data.b.data(), n_blocks);
    bench_consume_wide(data.wide);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

template <int K>
void BM_Fp8x1MulV2K(benchmark::State& state)
{
  const std::size_t n_blocks = static_cast<std::size_t>(state.range(0));
  const std::size_t n_packs = n_blocks * K;
  Fp8x1KBenchData data(n_packs);

  for (auto _ : state) {
    fp8x1_mul_k<K>(data.wide.data(), data.a.data(), data.b.data(), n_blocks);
    bench_consume_wide(data.wide);
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
    bench_consume_out(data.out);
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
    bench_consume_out(data.out);
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
    bench_consume_lazy(data.lazy);
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
    bench_consume_out(data.out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1UFpAdd(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1UFpBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      UFpPack8x1::add(data.wide_out[i], data.wide_a[i], data.wide_b[i]);
    }
    bench_consume_wide(data.wide_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1WideCopy(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1UFpBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      bench_copy_wide(data.wide_out[i], data.wide_a[i]);
    }
    bench_consume_wide(data.wide_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1WideAddRaw(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1UFpBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      bench_add_wide_raw(data.wide_out[i], data.wide_a[i], data.wide_b[i]);
    }
    bench_consume_wide(data.wide_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1UFpSub(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1UFpBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      UFpPack8x1::sub(data.wide_out[i], data.wide_a[i], data.wide_b[i]);
    }
    bench_consume_wide(data.wide_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1LUFpAdd(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1LUFpBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      LUFpPack8x1::add(data.lazy_out[i], data.lazy_a[i], data.lazy_b[i]);
    }
    bench_consume_luwide(data.lazy_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1LUFpSub(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1LUFpBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      LUFpPack8x1::sub(data.lazy_out[i], data.lazy_a[i], data.lazy_b[i]);
    }
    bench_consume_luwide(data.lazy_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

template <int Rounds>
void BM_Fp8x1UFpAddChain(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1WidePairBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      UFpPack8x1 x = data.wide_a[i];
      const UFpPack8x1 y = data.wide_b[i];
      for (int round = 0; round < Rounds; ++round) {
        UFpPack8x1::add(x, x, y);
      }
      data.wide_out[i] = x;
    }
    bench_consume_wide(data.wide_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs * static_cast<std::size_t>(Rounds));
}

template <int Rounds>
void BM_Fp8x1LUFpAddChain(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1WidePairBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      LUFpPack8x1 x = data.lazy_a[i];
      const LUFpPack8x1 y = data.lazy_b[i];
      for (int round = 0; round < Rounds; ++round) {
        LUFpPack8x1::add(x, x, y);
      }
      data.lazy_out[i] = x;
    }
    bench_consume_luwide(data.lazy_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs * static_cast<std::size_t>(Rounds));
}

template <int Rounds>
void BM_Fp8x1UFpAddSubChain(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1WidePairBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      UFpPack8x1 x = data.wide_a[i];
      UFpPack8x1 y = data.wide_b[i];
      for (int round = 0; round < Rounds; ++round) {
        UFpPack8x1::add(x, x, y);
        UFpPack8x1::sub(y, x, y);
      }
      data.wide_out[i] = x;
    }
    bench_consume_wide(data.wide_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs * static_cast<std::size_t>(2 * Rounds));
}

template <int Rounds>
void BM_Fp8x1LUFpAddSubChain(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1WidePairBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      LUFpPack8x1 x = data.lazy_a[i];
      LUFpPack8x1 y = data.lazy_b[i];
      for (int round = 0; round < Rounds; ++round) {
        LUFpPack8x1::add(x, x, y);
        LUFpPack8x1::sub(y, x, y);
      }
      data.lazy_out[i] = x;
    }
    bench_consume_luwide(data.lazy_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs * static_cast<std::size_t>(2 * Rounds));
}

template <LUFpRenormClass B>
void BM_Fp8x1LUFpRenorm(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1LUFpBenchData data(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      LUFpPack8x1::renorm<B>(data.lazy_out[i], data.lazy_in[i]);
    }
    bench_consume_luwide(data.lazy_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

template <int Rounds>
void BM_Fp8x1CarryKickOnlyChain(benchmark::State& state)
{
  cc512i limbs[2 * NWORDS]{};
  const cc512i kick = VSET1(static_cast<i64>(1ULL << BRADIX));

  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    limbs[i] = VZERO();
  }

  for (auto _ : state) {
    for (int round = 0; round < Rounds; ++round) {
      for (std::size_t i = 0; i < (2 * NWORDS) - 1; ++i) {
        limbs[i] = VADD(limbs[i], kick);
      }
    }
  }

  benchmark::DoNotOptimize(limbs[0]);
  benchmark::DoNotOptimize(limbs[(2 * NWORDS) - 1]);
  state.SetItemsProcessed(state.iterations() * static_cast<std::size_t>(Rounds));
}

template <int Rounds>
void BM_Fp8x1CarryChain(benchmark::State& state)
{
  cc512i limbs[2 * NWORDS]{};
  const cc512i kick = VSET1(static_cast<i64>(1ULL << BRADIX));

  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    limbs[i] = VZERO();
  }

  for (auto _ : state) {
    for (int round = 0; round < Rounds; ++round) {
      for (std::size_t i = 0; i < (2 * NWORDS) - 1; ++i) {
        limbs[i] = VADD(limbs[i], kick);
      }
      bench_signed_carry_range(limbs, 0, (2 * NWORDS) - 1);
    }
  }

  benchmark::DoNotOptimize(limbs[0]);
  benchmark::DoNotOptimize(limbs[(2 * NWORDS) - 1]);
  state.SetItemsProcessed(state.iterations() * static_cast<std::size_t>(Rounds));
}

template <int Steps>
BLST_FORCEINLINE void bench_repair_round(cc512i* limbs)
{
  if constexpr (Steps == 1) {
    bench_high_add_mod_masked<1>(limbs, 0xFF);
    const u8 ge_p = bench_high_ge_mod_mask_prefix_lex<1>(limbs);
    bench_high_sub_mod_masked<1>(limbs, ge_p);
  } else if constexpr (Steps == 2) {
    bench_high_add_mod_masked<3>(limbs, 0xFF);
    const u8 ge_2p = bench_high_ge_mod_mask_prefix_lex<2>(limbs);
    bench_high_sub_mod_masked<2>(limbs, ge_2p);
    const u8 ge_p = bench_high_ge_mod_mask_prefix_lex<1>(limbs);
    bench_high_sub_mod_masked<1>(limbs, ge_p);
  } else if constexpr (Steps == 4) {
    bench_high_add_mod_masked<7>(limbs, 0xFF);
    const u8 ge_4p = bench_high_ge_mod_mask_prefix_lex<4>(limbs);
    bench_high_sub_mod_masked<4>(limbs, ge_4p);
    const u8 ge_2p = bench_high_ge_mod_mask_prefix_lex<2>(limbs);
    bench_high_sub_mod_masked<2>(limbs, ge_2p);
    const u8 ge_p = bench_high_ge_mod_mask_prefix_lex<1>(limbs);
    bench_high_sub_mod_masked<1>(limbs, ge_p);
  } else {
    static_assert(Steps == 3, "Repair chain supports 1, 2, 3, or 4 steps.");
    bench_high_sub_mod_masked<1>(limbs, 0xFF);
    const u8 neg_mask = bench_mask_gt_epi64(VZERO(), limbs[(2 * NWORDS) - 1]);
    bench_high_add_mod_masked<4>(limbs, neg_mask);
    const u8 ge_2p = bench_high_ge_mod_mask_prefix_lex<2>(limbs);
    bench_high_sub_mod_masked<2>(limbs, ge_2p);
    const u8 ge_p = bench_high_ge_mod_mask_prefix_lex<1>(limbs);
    bench_high_sub_mod_masked<1>(limbs, ge_p);
  }
}

BLST_FORCEINLINE void bench_repair_round_fused_2step(cc512i* limbs)
{
  bench_high_add_mod_masked<3>(limbs, 0xFF);

  const BenchThreeWayBuckets ge = bench_high_ge_p_2p_3p_masks_fused(limbs);
  const u8 sub_3p = ge.ge_3p;
  const u8 sub_2p = static_cast<u8>(ge.ge_2p & static_cast<u8>(~ge.ge_3p));
  const u8 sub_1p = static_cast<u8>(ge.ge_p & static_cast<u8>(~ge.ge_2p));

  bench_high_sub_mod_masked<3>(limbs, sub_3p);
  bench_high_sub_mod_masked<2>(limbs, sub_2p);
  bench_high_sub_mod_masked<1>(limbs, sub_1p);
}

template <int Steps>
BLST_FORCEINLINE void bench_repair_apply_round(cc512i* limbs, const u8 lane_mask)
{
  if constexpr (Steps == 1) {
    bench_high_add_mod_masked<1>(limbs, lane_mask);
    bench_high_sub_mod_masked<1>(limbs, lane_mask);
  } else if constexpr (Steps == 2) {
    bench_high_add_mod_masked<3>(limbs, lane_mask);
    bench_high_sub_mod_masked<2>(limbs, lane_mask);
    bench_high_sub_mod_masked<1>(limbs, lane_mask);
  } else if constexpr (Steps == 4) {
    bench_high_add_mod_masked<7>(limbs, lane_mask);
    bench_high_sub_mod_masked<4>(limbs, lane_mask);
    bench_high_sub_mod_masked<2>(limbs, lane_mask);
    bench_high_sub_mod_masked<1>(limbs, lane_mask);
  } else {
    static_assert(Steps == 3, "Repair apply chain supports 1, 2, 3, or 4 steps.");
    bench_high_sub_mod_masked<1>(limbs, lane_mask);
    bench_high_add_mod_masked<4>(limbs, lane_mask);
    bench_high_sub_mod_masked<2>(limbs, lane_mask);
    bench_high_sub_mod_masked<1>(limbs, lane_mask);
  }
}

template <int Steps, int Rounds>
void BM_Fp8x1RepairChainImpl(benchmark::State& state)
{
  cc512i limbs[2 * NWORDS]{};

  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    limbs[i] = VZERO();
  }
  bench_seed_compare_high(limbs);

  for (auto _ : state) {
    for (int round = 0; round < Rounds; ++round) {
      bench_repair_round<Steps>(limbs);
    }
  }

  benchmark::DoNotOptimize(limbs[NWORDS]);
  benchmark::DoNotOptimize(limbs[(2 * NWORDS) - 1]);
  state.SetItemsProcessed(state.iterations() * static_cast<std::size_t>(Rounds));
}

template <int Steps, int Rounds>
void BM_Fp8x1RepairChainPrefixLex(benchmark::State& state)
{
  BM_Fp8x1RepairChainImpl<Steps, Rounds>(state);
}

template <int Rounds>
void BM_Fp8x1RepairChainFused3Way(benchmark::State& state)
{
  cc512i limbs[2 * NWORDS]{};

  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    limbs[i] = VZERO();
  }
  bench_seed_compare_high(limbs);

  for (auto _ : state) {
    for (int round = 0; round < Rounds; ++round) {
      bench_repair_round_fused_2step(limbs);
    }
  }

  benchmark::DoNotOptimize(limbs[NWORDS]);
  benchmark::DoNotOptimize(limbs[(2 * NWORDS) - 1]);
  state.SetItemsProcessed(state.iterations() * static_cast<std::size_t>(Rounds));
}

template <int Steps, int Rounds>
void BM_Fp8x1RepairApplyChain(benchmark::State& state)
{
  cc512i limbs[2 * NWORDS]{};
  u8 lane_mask = 0xA5;

  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    limbs[i] = VZERO();
  }

  for (auto _ : state) {
    for (int round = 0; round < Rounds; ++round) {
      lane_mask = static_cast<u8>((lane_mask << 1) | (lane_mask >> 7));
      bench_repair_apply_round<Steps>(limbs, lane_mask);
    }
  }

  benchmark::DoNotOptimize(lane_mask);
  benchmark::DoNotOptimize(limbs[NWORDS]);
  benchmark::DoNotOptimize(limbs[(2 * NWORDS) - 1]);
  state.SetItemsProcessed(state.iterations() * static_cast<std::size_t>(Rounds));
}

template <int Rounds>
void BM_Fp8x1RepairNegMaskChain(benchmark::State& state)
{
  cc512i top = VSET1(-1);
  u8 mask_acc = 0;

  for (auto _ : state) {
    for (int round = 0; round < Rounds; ++round) {
      mask_acc ^= bench_mask_gt_epi64(VZERO(), top);
    }
  }

  benchmark::DoNotOptimize(mask_acc);
  benchmark::DoNotOptimize(top);
  state.SetItemsProcessed(state.iterations() * static_cast<std::size_t>(Rounds));
}

template <int Rounds>
void BM_Fp8x1RepairHighKickChain(benchmark::State& state)
{
  cc512i limbs[2 * NWORDS]{};
  const cc512i kick = VSET1(1);

  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    limbs[i] = VZERO();
  }

  for (auto _ : state) {
    for (int round = 0; round < Rounds; ++round) {
      limbs[NWORDS] = VADD(limbs[NWORDS], kick);
    }
  }

  benchmark::DoNotOptimize(limbs[NWORDS]);
  state.SetItemsProcessed(state.iterations() * static_cast<std::size_t>(Rounds));
}

template <int Compares, int Rounds>
void BM_Fp8x1RepairCompareChainImpl(benchmark::State& state)
{
  cc512i limbs[2 * NWORDS]{};
  u8 mask_acc = 0;

  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    limbs[i] = VZERO();
  }
  bench_seed_compare_high(limbs);

  for (auto _ : state) {
    for (int round = 0; round < Rounds; ++round) {
      bench_mutate_compare_high(limbs);
      if constexpr (Compares >= 1) {
        mask_acc ^= bench_high_ge_mod_mask_prefix_lex<1>(limbs);
      }
      if constexpr (Compares >= 2) {
        mask_acc ^= bench_high_ge_mod_mask_prefix_lex<2>(limbs);
      }
    }
  }

  benchmark::DoNotOptimize(mask_acc);
  benchmark::DoNotOptimize(limbs[NWORDS]);
  state.SetItemsProcessed(state.iterations() * static_cast<std::size_t>(Rounds));
}

template <int Compares, int Rounds>
void BM_Fp8x1RepairCompareChainPrefixLex(benchmark::State& state)
{
  BM_Fp8x1RepairCompareChainImpl<Compares, Rounds>(state);
}

template <int Rounds>
void BM_Fp8x1RepairCompareChainFused3Way(benchmark::State& state)
{
  cc512i limbs[2 * NWORDS]{};
  u8 mask_acc = 0;

  for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
    limbs[i] = VZERO();
  }
  bench_seed_compare_high(limbs);

  for (auto _ : state) {
    for (int round = 0; round < Rounds; ++round) {
      bench_mutate_compare_high(limbs);
      const BenchThreeWayBuckets ge = bench_high_ge_p_2p_3p_masks_fused(limbs);
      mask_acc ^= ge.ge_p;
      mask_acc ^= ge.ge_2p;
      mask_acc ^= ge.ge_3p;
    }
  }

  benchmark::DoNotOptimize(mask_acc);
  benchmark::DoNotOptimize(limbs[NWORDS]);
  state.SetItemsProcessed(state.iterations() * static_cast<std::size_t>(Rounds));
}

void BM_Fp8x1LUFpNormalize(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1LUFpBenchData data(n_packs);
  std::vector<UFpPack8x1> wide_out(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      LUFpPack8x1::normalize(wide_out[i], data.lazy_in[i]);
    }
    bench_consume_wide(wide_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1LUFpNormalizeHiNeg4Pos4(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1LUFpBenchData data(n_packs);
  std::vector<UFpPack8x1> wide_out(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      LUFpPack8x1::normalize<LUFpNormalizeBound::HiNeg4Pos4>(
          wide_out[i], data.lazy_in_hi_neg4_pos4[i]);
    }
    bench_consume_wide(wide_out);
  }

  state.SetItemsProcessed(state.iterations() * n_packs);
}

void BM_Fp8x1LUFpNormalizeGenericHiNeg4Pos4(benchmark::State& state)
{
  const std::size_t n_packs = static_cast<std::size_t>(state.range(0));
  Fp8x1LUFpBenchData data(n_packs);
  std::vector<UFpPack8x1> wide_out(n_packs);

  for (auto _ : state) {
    for (std::size_t i = 0; i < n_packs; ++i) {
      LUFpPack8x1::normalize(wide_out[i], data.lazy_in_hi_neg4_pos4[i]);
    }
    bench_consume_wide(wide_out);
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
        FpPack8x1::mul(data.wide[0], data.cur[0], data.mul_by[0]);
        UFpPack8x1::redc(data.out[0], data.wide[0]);
        data.cur[0] = data.out[0];
      }
    }
    if constexpr (Chains >= 2) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul(data.wide[1], data.cur[1], data.mul_by[1]);
        UFpPack8x1::redc(data.out[1], data.wide[1]);
        data.cur[1] = data.out[1];
      }
    }
    if constexpr (Chains >= 3) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul(data.wide[2], data.cur[2], data.mul_by[2]);
        UFpPack8x1::redc(data.out[2], data.wide[2]);
        data.cur[2] = data.out[2];
      }
    }
    if constexpr (Chains >= 4) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul(data.wide[3], data.cur[3], data.mul_by[3]);
        UFpPack8x1::redc(data.out[3], data.wide[3]);
        data.cur[3] = data.out[3];
      }
    }
    if constexpr (Chains >= 5) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul(data.wide[4], data.cur[4], data.mul_by[4]);
        UFpPack8x1::redc(data.out[4], data.wide[4]);
        data.cur[4] = data.out[4];
      }
    }
    if constexpr (Chains >= 6) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul(data.wide[5], data.cur[5], data.mul_by[5]);
        UFpPack8x1::redc(data.out[5], data.wide[5]);
        data.cur[5] = data.out[5];
      }
    }
    if constexpr (Chains >= 7) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul(data.wide[6], data.cur[6], data.mul_by[6]);
        UFpPack8x1::redc(data.out[6], data.wide[6]);
        data.cur[6] = data.out[6];
      }
    }
    if constexpr (Chains >= 8) {
      for (std::size_t i = 0; i < kUnroll; ++i) {
        FpPack8x1::mul(data.wide[7], data.cur[7], data.mul_by[7]);
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
        FpPack8x1::mul(data.wide[0], data.cur[0], data.mul_by[0]);
        UFpPack8x1::redc(data.out[0], data.wide[0]);
        data.cur[0] = data.out[0];
      }
      if constexpr (Chains >= 2) {
        FpPack8x1::mul(data.wide[1], data.cur[1], data.mul_by[1]);
        UFpPack8x1::redc(data.out[1], data.wide[1]);
        data.cur[1] = data.out[1];
      }
      if constexpr (Chains >= 3) {
        FpPack8x1::mul(data.wide[2], data.cur[2], data.mul_by[2]);
        UFpPack8x1::redc(data.out[2], data.wide[2]);
        data.cur[2] = data.out[2];
      }
      if constexpr (Chains >= 4) {
        FpPack8x1::mul(data.wide[3], data.cur[3], data.mul_by[3]);
        UFpPack8x1::redc(data.out[3], data.wide[3]);
        data.cur[3] = data.out[3];
      }
      if constexpr (Chains >= 5) {
        FpPack8x1::mul(data.wide[4], data.cur[4], data.mul_by[4]);
        UFpPack8x1::redc(data.out[4], data.wide[4]);
        data.cur[4] = data.out[4];
      }
      if constexpr (Chains >= 6) {
        FpPack8x1::mul(data.wide[5], data.cur[5], data.mul_by[5]);
        UFpPack8x1::redc(data.out[5], data.wide[5]);
        data.cur[5] = data.out[5];
      }
      if constexpr (Chains >= 7) {
        FpPack8x1::mul(data.wide[6], data.cur[6], data.mul_by[6]);
        UFpPack8x1::redc(data.out[6], data.wide[6]);
        data.cur[6] = data.out[6];
      }
      if constexpr (Chains >= 8) {
        FpPack8x1::mul(data.wide[7], data.cur[7], data.mul_by[7]);
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
BENCHMARK(BM_Fp8x1MulV2K1_CppCandidateV6B)->Arg(256);
BENCHMARK(BM_Fp8x1MulV2K2_CppCandidate)->Arg(256);
BENCHMARK(BM_Fp8x1MulV2K2_CppCandidateV3)->Arg(256);
BENCHMARK(BM_Fp8x1MulV2K2_CppCandidateV4)->Arg(256);
BENCHMARK(BM_Fp8x1MulV2K2_CppCandidateV5)->Arg(256);
BENCHMARK(BM_Fp8x1MulV2K2_CppCandidateV6)->Arg(256);
BENCHMARK(BM_Fp8x1MulV2K2_CppCandidateV6A)->Arg(256);
BENCHMARK(BM_Fp8x1MulV2K2_CppCandidateV6B)->Arg(256);
#if BLST_HAVE_MUL_V2_K2_V6B_ASM
BENCHMARK(BM_Fp8x1MulV2K2_CppCandidateV6BAsm)->Arg(256);
#endif

BENCHMARK_TEMPLATE(BM_Fp8x1RedcK, 1)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1RedcK, 2)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1RedcK, 4)->Arg(256);

BENCHMARK_TEMPLATE(BM_Fp8x1MulRedcK, 1)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1MulRedcK, 2)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1MulRedcK, 4)->Arg(256);

BENCHMARK(BM_Fp8x1RedcCore)->Arg(256);
BENCHMARK(BM_Fp8x1Normalize)->Arg(256);
BENCHMARK(BM_Fp8x1WideCopy)->Arg(256);
BENCHMARK(BM_Fp8x1WideAddRaw)->Arg(256);
BENCHMARK(BM_Fp8x1UFpAdd)->Arg(256);
BENCHMARK(BM_Fp8x1UFpSub)->Arg(256);
BENCHMARK(BM_Fp8x1LUFpAdd)->Arg(256);
BENCHMARK(BM_Fp8x1LUFpSub)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1UFpAddChain, 16)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1LUFpAddChain, 16)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1UFpAddSubChain, 16)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1LUFpAddSubChain, 16)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1LUFpRenorm, LUFpRenormClass::Tight)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1LUFpRenorm, LUFpRenormClass::Medium)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1LUFpRenorm, LUFpRenormClass::Loose)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp8x1CarryKickOnlyChain, 128);
BENCHMARK_TEMPLATE(BM_Fp8x1CarryChain, 128);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairChainPrefixLex, 1, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairChainPrefixLex, 2, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairChainPrefixLex, 3, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairChainPrefixLex, 4, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairChainFused3Way, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairApplyChain, 1, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairApplyChain, 2, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairApplyChain, 3, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairApplyChain, 4, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairNegMaskChain, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairHighKickChain, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairCompareChainPrefixLex, 1, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairCompareChainPrefixLex, 2, 256);
BENCHMARK_TEMPLATE(BM_Fp8x1RepairCompareChainFused3Way, 256);
BENCHMARK(BM_Fp8x1LUFpNormalize)->Arg(256);
BENCHMARK(BM_Fp8x1LUFpNormalizeHiNeg4Pos4)->Arg(256);
BENCHMARK(BM_Fp8x1LUFpNormalizeGenericHiNeg4Pos4)->Arg(256);

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

BENCHMARK_TEMPLATE(BM_Fp8x1RedcCore_HandInterleavedChains, 2);

} // namespace

#endif // BLST_AVX512 || BLST_AVX512_PORTABLE
