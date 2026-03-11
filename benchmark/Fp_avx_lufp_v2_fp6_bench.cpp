#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "blst/fields/Fp.h"
#include "benchmark/internal/Fp_avx_lufp_v2_impl.h"

using namespace blst;
using namespace blst::bench_lufp_v2;

namespace
{
bool fp2_equal(const Fp2& a, const Fp2& b)
{
    return Fp::equal(a[0], b[0]) && Fp::equal(a[1], b[1]);
}

bool fp6_equal(const Fp6& a, const Fp6& b)
{
    return fp2_equal(a[0], b[0]) && fp2_equal(a[1], b[1]) && fp2_equal(a[2], b[2]);
}

void random_redc_fp(Fp& v, std::mt19937_64& gen)
{
    for (u64 i = 0; i < v.num_limb(); ++i) {
        v[i] = gen();
    }
    while (Fp::is_redc(v) == 0) {
        Fp::rshift(v, v, 1);
    }
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

//struct LUFp6Pack8x1Raw
//{
//    LUFp2Pack8x1x1 c[3];
//};

BLST_FORCEINLINE void prepare_lufp_v2_fp6_raw(UFp6Pack8x1& out, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
{
    Fp2Pack8x1x1 aa0{}, aa1{}, aa2{};
    Fp2Pack8x1x1 bb0{}, bb1{}, bb2{};
    UFp2Pack8x1x1 t0{}, t1{}, t2{};
    UFp2Pack8x1x1 w0{}, w1{}, w2{};

    Fp2Pack8x1x1::add(aa0, a.c[1], a.c[2]);
    Fp2Pack8x1x1::add(bb0, b.c[1], b.c[2]);
    Fp2Pack8x1x1::add(aa1, a.c[0], a.c[1]);
    Fp2Pack8x1x1::add(bb1, b.c[0], b.c[1]);
    Fp2Pack8x1x1::add(aa2, a.c[0], a.c[2]);
    Fp2Pack8x1x1::add(bb2, b.c[0], b.c[2]);

    UFp2Pack8x1x1::mul(t0, a.c[0], b.c[0]);
    UFp2Pack8x1x1::mul(t1, a.c[1], b.c[1]);
    UFp2Pack8x1x1::mul(t2, a.c[2], b.c[2]);
    UFp2Pack8x1x1::mul(w0, aa0, bb0);
    UFp2Pack8x1x1::mul(w1, aa1, bb1);
    UFp2Pack8x1x1::mul(w2, aa2, bb2);

    mul_fp6_pack_lufpv2_combine_only_raw(
        out.c[0],
        out.c[1],
        out.c[2],
        t0,
        t1,
        t2,
        w0,
        w1,
        w2);
}

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

void validate_lufp_v2_fp6()
{
    Fp6Pack8x1V2BenchData data(1);
    std::array<Fp6, 8> got{};
    std::array<Fp6, 8> ref{};

    for (std::size_t p = 0; p < data.n_packs_; ++p) {
        Fp6Pack8x1 wide_redc{};
        Fp6Pack8x1 accounting_redc{};
        Fp6Pack8x1 ref_mul{};
        UFp6Pack8x1 wide{};
        UFp6Pack8x1 accounting{};

        mul_fp6_pack_lufpv2_wide(wide, data.a_pack[p], data.b_pack[p]);
        UFp6Pack8x1::redc(wide_redc, wide);

        prepare_lufp_v2_fp6_raw(accounting, data.a_pack[p], data.b_pack[p]);
        lufp2_v2_renorm_fp6_mul_common(accounting.c[0]);
        lufp2_v2_renorm_fp6_mul_common(accounting.c[1]);
        lufp2_v2_renorm_fp6_mul_common(accounting.c[2]);
        UFp6Pack8x1::redc(accounting_redc, accounting);

        Fp6Pack8x1::mul(ref_mul, data.a_pack[p], data.b_pack[p]);

        Fp6Pack8x1::to_aos(got.data(), wide_redc);
        Fp6Pack8x1::to_aos(ref.data(), ref_mul);
        for (std::size_t lane = 0; lane < got.size(); ++lane) {
            if (!fp6_equal(got[lane], ref[lane])) {
                std::fprintf(stderr, "Fp6 LUFp v2 wide validation failed at pack %zu lane %zu\n", p, lane);
                std::abort();
            }
        }

        Fp6Pack8x1::to_aos(got.data(), accounting_redc);
        for (std::size_t lane = 0; lane < got.size(); ++lane) {
            if (!fp6_equal(got[lane], ref[lane])) {
                std::fprintf(stderr, "Fp6 LUFp v2 accounting validation failed at pack %zu lane %zu\n", p, lane);
                std::abort();
            }
        }
    }
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
                if constexpr (Stage >= LUFpFp6MulStage::Full) {
                    mul_fp6_pack_lufpv2_combine_only_raw(
                        data.out_wide[p].c[0],
                        data.out_wide[p].c[1],
                        data.out_wide[p].c[2],
                        t0,
                        t1,
                        t2,
                        w0,
                        w1,
                        w2);
                    lufp2_v2_renorm_fp6_mul_common(data.out_wide[p].c[0]);
                    lufp2_v2_renorm_fp6_mul_common(data.out_wide[p].c[1]);
                    lufp2_v2_renorm_fp6_mul_common(data.out_wide[p].c[2]);
                } else {
                    UFp2Pack8x1x1 c0{}, c1{}, c2{};
                    mul_fp6_pack_lufpv2_combine_only_raw(
                        c0,
                        c1,
                        c2,
                        t0,
                        t1,
                        t2,
                        w0,
                        w1,
                        w2);
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
    validate_lufp_v2_fp6();
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

enum class LUFpV2Fp6MulCommonRenormStage
{
    Bias,
    BiasCarry,
    BiasCarrySign4,
    BiasCarrySign42,
    BiasCarrySign421,
    BiasCarrySign4211,
    Full,
};

template <LUFpV2Fp6MulCommonRenormStage Stage>
BLST_FORCEINLINE void lufp_v2_fp6_mul_common_renorm_accounting(UFpPack8x1& in_out)
{
    v2_high_sub_const(in_out.data(), kV2Bias3PPlus3Vec);

    if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::BiasCarry) {
        u8 sign_mask = v2_carry_range_with_sign(in_out.data(), 0, (2 * NWORDS) - 1);

        if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::BiasCarrySign4) {
            sign_mask = v2_sign_step<4>(in_out.data(), sign_mask);

            if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::BiasCarrySign42) {
                sign_mask = v2_sign_step<2>(in_out.data(), sign_mask);

                if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::BiasCarrySign421) {
                    sign_mask = v2_sign_step<1>(in_out.data(), sign_mask);

                    if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::BiasCarrySign4211) {
                        benchmark::DoNotOptimize(v2_sign_step<1>(in_out.data(), sign_mask));
                    }
                }
            }
        }
    }

    if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::Full) {
        benchmark::DoNotOptimize(in_out[0]);
        benchmark::DoNotOptimize(in_out[(2 * NWORDS) - 1]);
    } else {
        benchmark::DoNotOptimize(in_out[0]);
        benchmark::DoNotOptimize(in_out[(2 * NWORDS) - 1]);
    }
}

template <LUFpV2Fp6MulCommonRenormStage Stage>
void BM_LUFpV2Fp6MulCommonRenormAccounting(benchmark::State& state)
{
    validate_lufp_v2_fp6();
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        state.PauseTiming();
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            prepare_lufp_v2_fp6_raw(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        state.ResumeTiming();

        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[0][0]);
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[0][1]);
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[1][0]);
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[1][1]);
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[2][0]);
            lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[2][1]);
        }

        benchmark::ClobberMemory();
        if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::Full) {
            bench_consume_fp6_wide(data.out_wide);
        }
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_) *
                            static_cast<int64_t>(6));
}

void BM_Fp6MulLufpWide(benchmark::State& state)
{
    validate_lufp_v2_fp6();
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            mul_fp6_pack_lufpv2_wide(data.out_wide[p], data.a_pack[p], data.b_pack[p]);
        }
        bench_consume_fp6_wide(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_));
}

BENCHMARK(BM_Fp6MulLufpWide)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAdd)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAddMul)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAddMulRawCombine)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::Full)->Arg(256);

BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::Bias)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::BiasCarry)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::BiasCarrySign4)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::BiasCarrySign42)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::BiasCarrySign421)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::BiasCarrySign4211)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::Full)->Arg(256);
}
