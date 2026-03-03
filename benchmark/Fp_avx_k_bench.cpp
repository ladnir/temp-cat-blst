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

#define BLST_APPLY_C2(STMT)                                                                 \
  do {                                                                                       \
    { constexpr std::size_t c = 0; STMT }                                                   \
    { constexpr std::size_t c = 1; STMT }                                                   \
  } while (0)

#define BLST_APPLY_C4(STMT)                                                                 \
  do {                                                                                       \
    { constexpr std::size_t c = 0; STMT }                                                   \
    { constexpr std::size_t c = 1; STMT }                                                   \
    { constexpr std::size_t c = 2; STMT }                                                   \
    { constexpr std::size_t c = 3; STMT }                                                   \
  } while (0)

#define BLST_APPLY_C8(STMT)                                                                 \
  do {                                                                                       \
    { constexpr std::size_t c = 0; STMT }                                                   \
    { constexpr std::size_t c = 1; STMT }                                                   \
    { constexpr std::size_t c = 2; STMT }                                                   \
    { constexpr std::size_t c = 3; STMT }                                                   \
    { constexpr std::size_t c = 4; STMT }                                                   \
    { constexpr std::size_t c = 5; STMT }                                                   \
    { constexpr std::size_t c = 6; STMT }                                                   \
    { constexpr std::size_t c = 7; STMT }                                                   \
  } while (0)

#define BLST_REDC_CORE_STEP(APPLY, IDX)                                                     \
  do {                                                                                       \
    APPLY({ u[c] = VAND(VMACLO(zero, a[(IDX)][c], montw), bmask); });                       \
    APPLY({                                                                                  \
      a[(IDX) + 0][c] = VMACLO(a[(IDX) + 0][c], u[c], p0);                                  \
      a[(IDX) + 1][c] = VMACLO(a[(IDX) + 1][c], u[c], p1);                                  \
      a[(IDX) + 2][c] = VMACLO(a[(IDX) + 2][c], u[c], p2);                                  \
      a[(IDX) + 3][c] = VMACLO(a[(IDX) + 3][c], u[c], p3);                                  \
      a[(IDX) + 4][c] = VMACLO(a[(IDX) + 4][c], u[c], p4);                                  \
      a[(IDX) + 5][c] = VMACLO(a[(IDX) + 5][c], u[c], p5);                                  \
      a[(IDX) + 6][c] = VMACLO(a[(IDX) + 6][c], u[c], p6);                                  \
      a[(IDX) + 7][c] = VMACLO(a[(IDX) + 7][c], u[c], p7);                                  \
    });                                                                                      \
    APPLY({                                                                                  \
      y[(IDX) + 0][c] = VMACHI(y[(IDX) + 0][c], u[c], p0);                                  \
      y[(IDX) + 1][c] = VMACHI(y[(IDX) + 1][c], u[c], p1);                                  \
      y[(IDX) + 2][c] = VMACHI(y[(IDX) + 2][c], u[c], p2);                                  \
      y[(IDX) + 3][c] = VMACHI(y[(IDX) + 3][c], u[c], p3);                                  \
      y[(IDX) + 4][c] = VMACHI(y[(IDX) + 4][c], u[c], p4);                                  \
      y[(IDX) + 5][c] = VMACHI(y[(IDX) + 5][c], u[c], p5);                                  \
      y[(IDX) + 6][c] = VMACHI(y[(IDX) + 6][c], u[c], p6);                                  \
      y[(IDX) + 7][c] = VMACHI(y[(IDX) + 7][c], u[c], p7);                                  \
    });                                                                                      \
    APPLY({                                                                                  \
      a[(IDX) + 1][c] =                                                                      \
          VADD(VADD(a[(IDX) + 1][c], VSRA(a[(IDX) + 0][c], BRADIX)), VSHL(y[(IDX) + 0][c], BALIGN)); \
    });                                                                                      \
  } while (0)

#define BLST_REDC_CORE_FINALIZE(APPLY)                                                      \
  do {                                                                                       \
    APPLY({ a[9][c] = VADD(a[9][c], VSHL(y[8][c], BALIGN)); });                             \
    APPLY({ a[10][c] = VADD(a[10][c], VSHL(y[9][c], BALIGN)); });                           \
    APPLY({ a[11][c] = VADD(a[11][c], VSHL(y[10][c], BALIGN)); });                          \
    APPLY({ a[12][c] = VADD(a[12][c], VSHL(y[11][c], BALIGN)); });                          \
    APPLY({ a[13][c] = VADD(a[13][c], VSHL(y[12][c], BALIGN)); });                          \
    APPLY({ a[14][c] = VADD(a[14][c], VSHL(y[13][c], BALIGN)); });                          \
    APPLY({ a[15][c] = VADD(a[15][c], VSHL(y[14][c], BALIGN)); });                          \
                                                                                             \
    APPLY({                                                                                  \
      r[0][c] = a[8][c]; r[1][c] = a[9][c]; r[2][c] = a[10][c]; r[3][c] = a[11][c];       \
      r[4][c] = a[12][c]; r[5][c] = a[13][c]; r[6][c] = a[14][c]; r[7][c] = a[15][c];     \
    });                                                                                      \
                                                                                             \
    APPLY({ r[1][c] = VADD(r[1][c], VSRA(r[0][c], BRADIX)); r[0][c] = VAND(r[0][c], bmask); }); \
    APPLY({ r[2][c] = VADD(r[2][c], VSRA(r[1][c], BRADIX)); r[1][c] = VAND(r[1][c], bmask); }); \
    APPLY({ r[3][c] = VADD(r[3][c], VSRA(r[2][c], BRADIX)); r[2][c] = VAND(r[2][c], bmask); }); \
    APPLY({ r[4][c] = VADD(r[4][c], VSRA(r[3][c], BRADIX)); r[3][c] = VAND(r[3][c], bmask); }); \
    APPLY({ r[5][c] = VADD(r[5][c], VSRA(r[4][c], BRADIX)); r[4][c] = VAND(r[4][c], bmask); }); \
    APPLY({ r[6][c] = VADD(r[6][c], VSRA(r[5][c], BRADIX)); r[5][c] = VAND(r[5][c], bmask); }); \
    APPLY({                                                                                  \
      r[7][c] = VADD(r[7][c], VSRA(r[6][c], BRADIX));                                       \
      r[6][c] = VAND(r[6][c], bmask);                                                        \
      r[7][c] = VAND(r[7][c], bmask);                                                        \
    });                                                                                      \
  } while (0)

template <int Chains>
BLST_FORCEINLINE void redc_core_hand_interleaved_once(
    std::array<LFpPack8x1, Chains>& out,
    const std::array<UFpPack8x1, Chains>& in);

#define BLST_DEFINE_REDC_CORE_HAND_SPECIALIZATION(CHAINS, APPLY)                            \
  template <>                                                                                \
  BLST_FORCEINLINE void redc_core_hand_interleaved_once<CHAINS>(                            \
      std::array<LFpPack8x1, CHAINS>& out,                                                  \
      const std::array<UFpPack8x1, CHAINS>& in)                                             \
  {                                                                                          \
    std::array<std::array<cc512i, CHAINS>, 16> a{};                                         \
    std::array<std::array<cc512i, CHAINS>, 15> y{};                                         \
    std::array<std::array<cc512i, CHAINS>, 8> r{};                                          \
    std::array<cc512i, CHAINS> u{};                                                         \
                                                                                             \
    const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);                \
    const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);                \
    const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);                                    \
    const cc512i bmask = VSET1(BMASK), montw = VSET1(MONT_W_R48), zero = VZERO();          \
                                                                                             \
    for (std::size_t i = 0; i < 16; ++i) {                                                  \
      APPLY({ a[i][c] = in[c][i]; });                                                       \
    }                                                                                        \
    for (std::size_t i = 0; i < 15; ++i) {                                                  \
      APPLY({ y[i][c] = VZERO(); });                                                        \
    }                                                                                        \
                                                                                             \
    BLST_REDC_CORE_STEP(APPLY, 0);                                                           \
    BLST_REDC_CORE_STEP(APPLY, 1);                                                           \
    BLST_REDC_CORE_STEP(APPLY, 2);                                                           \
    BLST_REDC_CORE_STEP(APPLY, 3);                                                           \
    BLST_REDC_CORE_STEP(APPLY, 4);                                                           \
    BLST_REDC_CORE_STEP(APPLY, 5);                                                           \
    BLST_REDC_CORE_STEP(APPLY, 6);                                                           \
    BLST_REDC_CORE_STEP(APPLY, 7);                                                           \
    BLST_REDC_CORE_FINALIZE(APPLY);                                                          \
                                                                                             \
    for (std::size_t i = 0; i < 8; ++i) {                                                   \
      APPLY({ out[c][i] = r[i][c]; });                                                      \
    }                                                                                        \
  }

BLST_DEFINE_REDC_CORE_HAND_SPECIALIZATION(2, BLST_APPLY_C2);
BLST_DEFINE_REDC_CORE_HAND_SPECIALIZATION(4, BLST_APPLY_C4);
BLST_DEFINE_REDC_CORE_HAND_SPECIALIZATION(8, BLST_APPLY_C8);

#undef BLST_DEFINE_REDC_CORE_HAND_SPECIALIZATION
#undef BLST_REDC_CORE_FINALIZE
#undef BLST_REDC_CORE_STEP

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
      redc_core_hand_interleaved_once<2>(data.lazy, data.wide);
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
      redc_core_hand_interleaved_once<4>(data.lazy, data.wide);
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
      redc_core_hand_interleaved_once<8>(data.lazy, data.wide);
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

#undef BLST_APPLY_C2
#undef BLST_APPLY_C4
#undef BLST_APPLY_C8

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
