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

BLST_FORCEINLINE void bench_consume_fp6_wide(const std::vector<UFp6Pack8x1>& data)
{
    if (data.empty()) {
        return;
    }
    benchmark::DoNotOptimize(data.front().c[0][0][0]);
    benchmark::DoNotOptimize(data.back().c[2][1][NWORDS - 1]);
    benchmark::ClobberMemory();
}

struct LUFp6Pack8x1V2Raw
{
    LUFp2Pack8x1x1V2 c[3];
};

static_assert(sizeof(UFp2Pack8x1x1) == sizeof(LUFp2Pack8x1x1V2));
static_assert(alignof(UFp2Pack8x1x1) == alignof(LUFp2Pack8x1x1V2));

BLST_FORCEINLINE const LUFp2Pack8x1x1V2& bench_alias_as_v2(const UFp2Pack8x1x1& src)
{
    return reinterpret_cast<const LUFp2Pack8x1x1V2&>(src);
}

struct Fp6Pack8x1V2BenchData
{
    explicit Fp6Pack8x1V2BenchData(const std::size_t n_packs)
        : n_packs_(n_packs)
        , a_scalar(n_packs_ * 8)
        , b_scalar(n_packs_ * 8)
        , a_pack(n_packs_)
        , b_pack(n_packs_)
        , raw_lazy_out(n_packs_)
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
            Fp2Pack8x1x1 aa0{}, aa1{}, aa2{};
            Fp2Pack8x1x1 bb0{}, bb1{}, bb2{};
            UFp2Pack8x1x1 t0{}, t1{}, t2{};
            UFp2Pack8x1x1 w0{}, w1{}, w2{};

            Fp6Pack8x1::from_aos(a_pack[p], a_scalar.data() + (8 * p));
            Fp6Pack8x1::from_aos(b_pack[p], b_scalar.data() + (8 * p));

            Fp2Pack8x1x1::add(aa0, a_pack[p].c[1], a_pack[p].c[2]);
            Fp2Pack8x1x1::add(bb0, b_pack[p].c[1], b_pack[p].c[2]);
            Fp2Pack8x1x1::add(aa1, a_pack[p].c[0], a_pack[p].c[1]);
            Fp2Pack8x1x1::add(bb1, b_pack[p].c[0], b_pack[p].c[1]);
            Fp2Pack8x1x1::add(aa2, a_pack[p].c[0], a_pack[p].c[2]);
            Fp2Pack8x1x1::add(bb2, b_pack[p].c[0], b_pack[p].c[2]);

            UFp2Pack8x1x1::mul(t0, a_pack[p].c[0], b_pack[p].c[0]);
            UFp2Pack8x1x1::mul(t1, a_pack[p].c[1], b_pack[p].c[1]);
            UFp2Pack8x1x1::mul(t2, a_pack[p].c[2], b_pack[p].c[2]);
            UFp2Pack8x1x1::mul(w0, aa0, bb0);
            UFp2Pack8x1x1::mul(w1, aa1, bb1);
            UFp2Pack8x1x1::mul(w2, aa2, bb2);

            mul_fp6_pack_lufpv2_combine_only_raw(
                raw_lazy_out[p].c[0],
                raw_lazy_out[p].c[1],
                raw_lazy_out[p].c[2],
                bench_alias_as_v2(t0),
                bench_alias_as_v2(t1),
                bench_alias_as_v2(t2),
                bench_alias_as_v2(w0),
                bench_alias_as_v2(w1),
                bench_alias_as_v2(w2));
        }
    }

    const std::size_t n_packs_;
    std::vector<Fp6> a_scalar;
    std::vector<Fp6> b_scalar;
    std::vector<Fp6Pack8x1> a_pack;
    std::vector<Fp6Pack8x1> b_pack;
    std::vector<LUFp6Pack8x1V2Raw> raw_lazy_out;
    std::vector<UFp6Pack8x1> out_wide;
};

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
                mul_fp6_pack_lufpv2_combine_only_raw(
                    c0,
                    c1,
                    c2,
                    bench_alias_as_v2(t0),
                    bench_alias_as_v2(t1),
                    bench_alias_as_v2(t2),
                    bench_alias_as_v2(w0),
                    bench_alias_as_v2(w1),
                    bench_alias_as_v2(w2));

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

enum class LUFpV2Fp6MulCommonRenormStage
{
    CopyBias,
    CopyBiasCarry,
    CopyBiasCarrySign4,
    CopyBiasCarrySign42,
    CopyBiasCarrySign421,
    CopyBiasCarrySign4211,
    Full,
};

template <LUFpV2Fp6MulCommonRenormStage Stage>
BLST_FORCEINLINE void lufp_v2_fp6_mul_common_renorm_accounting(UFpPack8x1& r, const LUFpPack8x1V2& a)
{
    cc512i limbs[2 * NWORDS]{};
    v2_copy_limbs(limbs, a);
    v2_high_sub_const(limbs, kV2Bias3PPlus3Vec);

    if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::CopyBiasCarry) {
        u8 sign_mask = v2_carry_range_with_sign(limbs, 0, (2 * NWORDS) - 1);

        if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign4) {
            sign_mask = v2_sign_step<4>(limbs, sign_mask);

            if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign42) {
                sign_mask = v2_sign_step<2>(limbs, sign_mask);

                if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign421) {
                    sign_mask = v2_sign_step<1>(limbs, sign_mask);

                    if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign4211) {
                        benchmark::DoNotOptimize(v2_sign_step<1>(limbs, sign_mask));
                    }
                }
            }
        }
    }

    if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::Full) {
        v2_store_wide(r, limbs);
    } else {
        benchmark::DoNotOptimize(limbs[0]);
        benchmark::DoNotOptimize(limbs[(2 * NWORDS) - 1]);
    }
}

template <LUFpV2Fp6MulCommonRenormStage Stage>
void BM_LUFpV2Fp6MulCommonRenormAccounting(benchmark::State& state)
{
    Fp6Pack8x1V2BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::Full) {
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[0][0], data.raw_lazy_out[p].c[0][0]);
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[0][1], data.raw_lazy_out[p].c[0][1]);
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[1][0], data.raw_lazy_out[p].c[1][0]);
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[1][1], data.raw_lazy_out[p].c[1][1]);
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[2][0], data.raw_lazy_out[p].c[2][0]);
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(data.out_wide[p].c[2][1], data.raw_lazy_out[p].c[2][1]);
            } else {
                UFpPack8x1 out{};
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(out, data.raw_lazy_out[p].c[0][0]);
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(out, data.raw_lazy_out[p].c[0][1]);
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(out, data.raw_lazy_out[p].c[1][0]);
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(out, data.raw_lazy_out[p].c[1][1]);
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(out, data.raw_lazy_out[p].c[2][0]);
                lufp_v2_fp6_mul_common_renorm_accounting<Stage>(out, data.raw_lazy_out[p].c[2][1]);
            }
        }

        if constexpr (Stage >= LUFpV2Fp6MulCommonRenormStage::Full) {
            bench_consume_fp6_wide(data.out_wide);
        } else {
            benchmark::ClobberMemory();
        }
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_) *
                            static_cast<int64_t>(6));
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

BENCHMARK(BM_Fp6MulLufpWide)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAdd)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAddMul)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::PreAddMulRawCombine)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulLufpAccounting, LUFpFp6MulStage::Full)->Arg(256);

BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::CopyBias)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::CopyBiasCarry)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign4)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign42)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign421)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign4211)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpV2Fp6MulCommonRenormAccounting, LUFpV2Fp6MulCommonRenormStage::Full)->Arg(256);
}
