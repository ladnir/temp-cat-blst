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

BLST_FORCEINLINE void bench_consume_fp6_wide(const std::vector<UFp6Pack8x1>& data)
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
        , out_wide(n_packs_)
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
    std::vector<UFp6Pack8x1> out_wide;
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
        LUFpPack8x1V2 lazy{};
        bench_copy_as_v2(reinterpret_cast<LUFpPack8x1V2&>(lazy), a);
        LUFpPack8x1V2::add(lazy, lazy, reinterpret_cast<const LUFpPack8x1V2&>(b));
        for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
            r[i] = lazy[i];
        }
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
            LUFpPack8x1V2 acc{};
            bench_copy_as_v2(acc, data.wide_a[p]);
            for (std::size_t round = 0; round < kRounds; ++round) {
                LUFpPack8x1V2::add(acc, acc, data.lazy_b[p]);
            }
            data.lazy_tmp[p] = acc;
        }
        bench_consume_luwide_v2(data.lazy_tmp);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_) *
                            static_cast<int64_t>(kRounds));
}

enum class LUFpFp6MulStage
{
    PreAdd,
    PreAddMul,
    PreAddMulRawCombine,
    Full,
};

template <LUFpFp6MulStage Stage>
BLST_FORCEINLINE void lufp_fp6_mul_accounting_v2(Fp6Pack8x1V2BenchData& data)
{
    for (std::size_t p = 0; p < data.n_packs_; ++p) {
        Fp2Pack8x1x1 aa0{}, aa1{}, aa2{};
        Fp2Pack8x1x1 bb0{}, bb1{}, bb2{};
        UFp2Pack8x1x1 t0{}, t1{}, t2{};
        UFp2Pack8x1x1 w0{}, w1{}, w2{};

        Fp2Pack8x1x1::add(aa0, data.a_pack[p].c[1], data.a_pack[p].c[2]);
        Fp2Pack8x1x1::add(bb0, data.b_pack[p].c[1], data.b_pack[p].c[2]);
        Fp2Pack8x1x1::add(aa1, data.a_pack[p].c[0], data.a_pack[p].c[1]);
        Fp2Pack8x1x1::add(bb1, data.b_pack[p].c[0], data.b_pack[p].c[1]);
        Fp2Pack8x1x1::add(aa2, data.a_pack[p].c[0], data.a_pack[p].c[2]);
        Fp2Pack8x1x1::add(bb2, data.b_pack[p].c[0], data.b_pack[p].c[2]);

        if constexpr (Stage >= LUFpFp6MulStage::PreAddMul) {
            UFp2Pack8x1x1::mul(t0, data.a_pack[p].c[0], data.b_pack[p].c[0]);
            UFp2Pack8x1x1::mul(t1, data.a_pack[p].c[1], data.b_pack[p].c[1]);
            UFp2Pack8x1x1::mul(t2, data.a_pack[p].c[2], data.b_pack[p].c[2]);
            UFp2Pack8x1x1::mul(w0, aa0, bb0);
            UFp2Pack8x1x1::mul(w1, aa1, bb1);
            UFp2Pack8x1x1::mul(w2, aa2, bb2);

            if constexpr (Stage >= LUFpFp6MulStage::PreAddMulRawCombine) {
                LUFp2Pack8x1x1V2 c0{}, c1{}, c2{};
                LUFp2Pack8x1x1V2 vt0{}, vt1{}, vt2{}, vw0{}, vw1{}, vw2{};
                bench_copy_as_v2(vt0, t0);
                bench_copy_as_v2(vt1, t1);
                bench_copy_as_v2(vt2, t2);
                bench_copy_as_v2(vw0, w0);
                bench_copy_as_v2(vw1, w1);
                bench_copy_as_v2(vw2, w2);
                mul_fp6_pack_lufpv2_combine_only_raw(c0, c1, c2, vt0, vt1, vt2, vw0, vw1, vw2);

                if constexpr (Stage >= LUFpFp6MulStage::Full) {
                    LUFp2Pack8x1x1V2::renorm<LUFpV2RenormClass::Fp6MulCommon>(data.out_wide[p].c[0], c0);
                    LUFp2Pack8x1x1V2::renorm<LUFpV2RenormClass::Fp6MulCommon>(data.out_wide[p].c[1], c1);
                    LUFp2Pack8x1x1V2::renorm<LUFpV2RenormClass::Fp6MulCommon>(data.out_wide[p].c[2], c2);
                } else {
                    benchmark::DoNotOptimize(c0[0][0]);
                    benchmark::DoNotOptimize(c2[1][NWORDS - 1]);
                }
            } else {
                benchmark::DoNotOptimize(t0[0][0]);
                benchmark::DoNotOptimize(w2[1][NWORDS - 1]);
            }
        } else {
            benchmark::DoNotOptimize(aa0[0][0]);
            benchmark::DoNotOptimize(bb2[1][NWORDS - 1]);
        }
    }
}

template <LUFpFp6MulStage Stage>
void BM_Fp6MulLufpAccounting(benchmark::State& state)
{
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        lufp_fp6_mul_accounting_v2<Stage>(data);
        if constexpr (Stage >= LUFpFp6MulStage::Full) {
            bench_consume_fp6_wide(data.out_wide);
        } else {
            benchmark::ClobberMemory();
        }
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

void BM_Fp6MulLufpWide(benchmark::State& state)
{
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        lufp_fp6_mul_accounting_v2<LUFpFp6MulStage::Full>(data);
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

BENCHMARK(BM_LUFpAddChain)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpAddAccounting, LUFpAddStage::RawAdd)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpAddAccounting, LUFpAddStage::RawAddBias)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpAddAccounting, LUFpAddStage::RawAddBiasCarry)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpAddAccounting, LUFpAddStage::Full)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpAddAccounting, LUFpAddStage::FullSign4)->Arg(256);

BENCHMARK(BM_Fp6MulLufpWide)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAdd)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAddMul)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAddMulRawCombine)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::Full)->Arg(256);
}
