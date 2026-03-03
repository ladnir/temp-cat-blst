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

} // namespace

#endif // BLST_AVX512
