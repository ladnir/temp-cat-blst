#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include <array>
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

BLST_FORCEINLINE void consume_wide(const std::vector<UFpPack8x1>& data)
{
    if (data.empty()) {
        return;
    }
    benchmark::DoNotOptimize(data.front()[0]);
    benchmark::DoNotOptimize(data.back()[(2 * NWORDS) - 1]);
    benchmark::ClobberMemory();
}

struct LUFpV2ExactAddRenormBenchData
{
    explicit LUFpV2ExactAddRenormBenchData(const std::size_t n_packs)
        : n_packs_(n_packs)
        , lazy_in(n_packs_)
        , wide_out(n_packs_)
    {
        std::mt19937_64 gen(0xA6DDEB4A5512C4B1ULL + static_cast<u64>(n_packs_));

        for (std::size_t i = 0; i < n_packs_; ++i) {
            std::array<Fp, 8> lane_a{};
            std::array<Fp, 8> lane_b{};
            std::array<Fp, 8> lane_c{};
            std::array<Fp, 8> lane_d{};
            FpPack8x1 pack_a{}, pack_b{}, pack_c{}, pack_d{};
            UFpPack8x1 wide_a{}, wide_b{};

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

            FpPack8x1::mul(wide_a, pack_a, pack_b);
            FpPack8x1::mul(wide_b, pack_c, pack_d);

            bench_copy_as_v2(lazy_in[i], wide_a);
            LUFpPack8x1V2::add(lazy_in[i], lazy_in[i], reinterpret_cast<const LUFpPack8x1V2&>(wide_b));
        }
    }

    const std::size_t n_packs_;
    std::vector<LUFpPack8x1V2> lazy_in;
    std::vector<UFpPack8x1> wide_out;
};

enum class LUFpV2ExactAddRenormStage
{
    Bias,
    BiasCarry,
    Full,
    FullSign4,
};

template <LUFpV2ExactAddRenormStage Stage>
BLST_FORCEINLINE void lufp_v2_exact_add_renorm_accounting(UFpPack8x1& r, const LUFpPack8x1V2& a)
{
    v2_bias_neg_p(r, a);

    if constexpr (Stage >= LUFpV2ExactAddRenormStage::BiasCarry) {
        const cc512i smask = v2_carry_from_biased(r);
        benchmark::DoNotOptimize(smask);

        if constexpr (Stage == LUFpV2ExactAddRenormStage::Full) {
            v2_repair_one_step_from_sign(r, smask);
        } else if constexpr (Stage == LUFpV2ExactAddRenormStage::FullSign4) {
            cc512i s = v2_repair_sign_step<4>(r, smask);
            s = v2_repair_sign_step<2>(r, s);
            s = v2_repair_sign_step<1>(r, s);
            benchmark::DoNotOptimize(v2_repair_sign_step<1>(r, s));
        }
    }
}

template <LUFpV2ExactAddRenormStage Stage>
void BM_LUFpV2ExactAddRenormAccounting(benchmark::State& state)
{
    constexpr std::size_t kRounds = 16;
    LUFpV2ExactAddRenormBenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            UFpPack8x1 out{};
            for (std::size_t round = 0; round < kRounds; ++round) {
                lufp_v2_exact_add_renorm_accounting<Stage>(out, data.lazy_in[p]);
            }
            data.wide_out[p] = out;
        }
        consume_wide(data.wide_out);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_) *
                            static_cast<int64_t>(kRounds));
}

BENCHMARK_TEMPLATE(BM_LUFpV2ExactAddRenormAccounting, LUFpV2ExactAddRenormStage::Bias)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2ExactAddRenormAccounting, LUFpV2ExactAddRenormStage::BiasCarry)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2ExactAddRenormAccounting, LUFpV2ExactAddRenormStage::Full)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2ExactAddRenormAccounting, LUFpV2ExactAddRenormStage::FullSign4)->Arg(256);
}

#endif
