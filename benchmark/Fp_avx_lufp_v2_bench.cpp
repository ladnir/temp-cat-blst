#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#include <cstddef>
#include <random>
#include <vector>

#include "blst/fields/Fp.h"
#include "benchmark/internal/Fp_avx_lufp_v2_impl.h"

using namespace blst;
using namespace blst::bench_lufp_v2;

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

BLST_FORCEINLINE void bench_consume_wide_v2(const std::vector<UFpPack8x1>& data)
{
    if (data.empty()) {
        return;
    }
    benchmark::DoNotOptimize(data.front()[0]);
    benchmark::DoNotOptimize(data.back()[(2 * NWORDS) - 1]);
    benchmark::ClobberMemory();
}

BLST_FORCEINLINE void bench_consume_luwide_v2(const std::vector<UFpPack8x1>& data)
{
    if (data.empty()) {
        return;
    }
    benchmark::DoNotOptimize(data.front()[0]);
    benchmark::DoNotOptimize(data.back()[(2 * NWORDS) - 1]);
    benchmark::ClobberMemory();
}

struct Fp8x1LUFpV2BenchData
{
    explicit Fp8x1LUFpV2BenchData(const std::size_t n_packs)
        : n_packs_(n_packs)
        , wide_a(n_packs_)
        , wide_b(n_packs_)
        , lazy_tmp(n_packs_)
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
        }
    }

    const std::size_t n_packs_;
    std::vector<UFpPack8x1> wide_a;
    std::vector<UFpPack8x1> wide_b;
    std::vector<UFpPack8x1> lazy_tmp;
    std::vector<UFpPack8x1> wide_out;
};

enum class LUFpAddStage
{
    RawAdd,
    RawAddBias,
    RawAddBiasCarry,
    Full,
    FullSign4,
};

template <LUFpAddStage Stage>
BLST_FORCEINLINE void lufp_add_accounting_v2(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    if constexpr (Stage == LUFpAddStage::RawAdd) {
        r = a;
        lufp_v2_add(r, r, b);
    } else {
        UFpPack8x1 wide{};
        v2_add_bias_neg_p(wide, a, b);

        if constexpr (Stage >= LUFpAddStage::RawAddBiasCarry) {
            const cc512i smask = v2_carry_from_biased(wide);
            benchmark::DoNotOptimize(smask);

            if constexpr (Stage == LUFpAddStage::Full) {
                v2_repair_one_step_from_sign(wide, smask);
            } else if constexpr (Stage == LUFpAddStage::FullSign4) {
                cc512i s = v2_repair_sign_step<4>(wide, smask);
                s = v2_repair_sign_step<2>(wide, s);
                s = v2_repair_sign_step<1>(wide, s);
                benchmark::DoNotOptimize(v2_repair_sign_step<1>(wide, s));
            }
        }

        r = wide;
    }
}

template <LUFpAddStage Stage>
void BM_LUFpAddAccounting(benchmark::State& state)
{
    constexpr std::size_t kRounds = 16;
    Fp8x1LUFpV2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            UFpPack8x1 acc = data.wide_a[p];
            for (std::size_t round = 0; round < kRounds; ++round) {
                lufp_add_accounting_v2<Stage>(acc, acc, data.wide_b[p]);
            }
            data.wide_out[p] = acc;
        }
        bench_consume_wide_v2(data.wide_out);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_) *
                            static_cast<int64_t>(kRounds));
}

void BM_LUFpAddChain(benchmark::State& state)
{
    constexpr std::size_t kRounds = 16;
    Fp8x1LUFpV2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            UFpPack8x1 acc = data.wide_a[p];
            for (std::size_t round = 0; round < kRounds; ++round) {
                lufp_v2_add(acc, acc, data.wide_b[p]);
            }
            data.lazy_tmp[p] = acc;
        }
        bench_consume_luwide_v2(data.lazy_tmp);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_) *
                            static_cast<int64_t>(kRounds));
}

BENCHMARK(BM_LUFpAddChain)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpAddAccounting, LUFpAddStage::RawAdd)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpAddAccounting, LUFpAddStage::RawAddBias)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpAddAccounting, LUFpAddStage::RawAddBiasCarry)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpAddAccounting, LUFpAddStage::Full)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpAddAccounting, LUFpAddStage::FullSign4)->Arg(256);
}
