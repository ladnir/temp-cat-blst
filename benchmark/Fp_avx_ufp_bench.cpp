#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include "benchmark/internal/Fp_avx_bench_kernels.h"

using namespace blst;
using namespace blst::bench_internal;

namespace
{
enum class UFpAddStage
{
    RawAdd,
    RawAddBias,
    RawAddBiasCarry,
    Full,
};

template <UFpAddStage Stage>
BLST_FORCEINLINE void ufp_add_accounting(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
{
    const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
    const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
    const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
    const cc512i bmask = VSET1(BMASK);

    cc512i r0 = VADD(a[0], b[0]), r1 = VADD(a[1], b[1]);
    cc512i r2 = VADD(a[2], b[2]), r3 = VADD(a[3], b[3]);
    cc512i r4 = VADD(a[4], b[4]), r5 = VADD(a[5], b[5]);
    cc512i r6 = VADD(a[6], b[6]), r7 = VADD(a[7], b[7]);
    cc512i r8 = VADD(a[8], b[8]), r9 = VADD(a[9], b[9]);
    cc512i r10 = VADD(a[10], b[10]), r11 = VADD(a[11], b[11]);
    cc512i r12 = VADD(a[12], b[12]), r13 = VADD(a[13], b[13]);
    cc512i r14 = VADD(a[14], b[14]), r15 = VADD(a[15], b[15]);

    if constexpr (Stage >= UFpAddStage::RawAddBias) {
        r8 = VSUB(r8, p0); r9 = VSUB(r9, p1);
        r10 = VSUB(r10, p2); r11 = VSUB(r11, p3);
        r12 = VSUB(r12, p4); r13 = VSUB(r13, p5);
        r14 = VSUB(r14, p6); r15 = VSUB(r15, p7);
    }

    cc512i smask = VZERO();
    if constexpr (Stage >= UFpAddStage::RawAddBiasCarry) {
        r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, bmask);
        r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, bmask);
        r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, bmask);
        r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, bmask);
        r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, bmask);
        r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, bmask);
        r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, bmask);
        r8 = VADD(r8, VSRA(r7, BRADIX)); r7 = VAND(r7, bmask);
        cc512i t0 = VADD(r9, VSRA(r8, BRADIX));
        t0 = VADD(r10, VSRA(t0, BRADIX));
        t0 = VADD(r11, VSRA(t0, BRADIX));
        t0 = VADD(r12, VSRA(t0, BRADIX));
        t0 = VADD(r13, VSRA(t0, BRADIX));
        t0 = VADD(r14, VSRA(t0, BRADIX));
        t0 = VADD(r15, VSRA(t0, BRADIX));
        smask = VSRA(t0, 63);

        if constexpr (Stage >= UFpAddStage::Full) {
            r8 = VADD(r8, VAND(p0, smask)); r9 = VADD(r9, VAND(p1, smask));
            r10 = VADD(r10, VAND(p2, smask)); r11 = VADD(r11, VAND(p3, smask));
            r12 = VADD(r12, VAND(p4, smask)); r13 = VADD(r13, VAND(p5, smask));
            r14 = VADD(r14, VAND(p6, smask)); r15 = VADD(r15, VAND(p7, smask));

            r9 = VADD(r9, VSRA(r8, BRADIX)); r8 = VAND(r8, bmask);
            r10 = VADD(r10, VSRA(r9, BRADIX)); r9 = VAND(r9, bmask);
            r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, bmask);
            r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, bmask);
            r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, bmask);
            r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, bmask);
            r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, bmask);
            r15 = VAND(r15, bmask);
        }
    }

    r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
    r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
    r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
    r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15;
    benchmark::DoNotOptimize(smask);
}

template <UFpAddStage Stage>
void BM_UFpAddAccounting(benchmark::State& state)
{
    constexpr std::size_t kRounds = 16;
    Fp8x1WideBenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            UFpPack8x1 acc = data.wide_a[p];
            for (std::size_t round = 0; round < kRounds; ++round) {
                ufp_add_accounting<Stage>(acc, acc, data.wide_b[p]);
            }
            data.wide_out[p] = acc;
        }
        benchmark::DoNotOptimize(data.wide_out.data());
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_) *
                            static_cast<int64_t>(kRounds));
}

void BM_UFpAddChain(benchmark::State& state)
{
    constexpr std::size_t kRounds = 16;
    Fp8x1WideBenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            UFpPack8x1 acc = data.wide_a[p];
            for (std::size_t round = 0; round < kRounds; ++round) {
                UFpPack8x1::add(acc, acc, data.wide_b[p]);
            }
            data.wide_out[p] = acc;
        }
        consume_wide(data.wide_out);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs_) *
                            static_cast<int64_t>(kRounds));
}

enum class UFpFp6MulStage
{
    PreAdd,
    PreAddMul,
    Full,
};

template <UFpFp6MulStage Stage>
BLST_FORCEINLINE void ufp_fp6_mul_accounting(Fp6BenchData& data)
{
    for (std::size_t p = 0; p < data.n_packs; ++p) {
        Fp2Pack8x1x1 aa0{}, aa1{}, aa2{};
        Fp2Pack8x1x1 bb0{}, bb1{}, bb2{};
        UFp2Pack8x1x1 t0{}, t1{}, t2{};
        UFp2Pack8x1x1 w0{}, w1{}, w2{};
        UFp2Pack8x1x1 t3{};
        UFp6Pack8x1 out{};

        Fp2Pack8x1x1::add(aa0, data.a_pack[p].c[1], data.a_pack[p].c[2]);
        Fp2Pack8x1x1::add(bb0, data.b_pack[p].c[1], data.b_pack[p].c[2]);
        Fp2Pack8x1x1::add(aa1, data.a_pack[p].c[0], data.a_pack[p].c[1]);
        Fp2Pack8x1x1::add(bb1, data.b_pack[p].c[0], data.b_pack[p].c[1]);
        Fp2Pack8x1x1::add(aa2, data.a_pack[p].c[0], data.a_pack[p].c[2]);
        Fp2Pack8x1x1::add(bb2, data.b_pack[p].c[0], data.b_pack[p].c[2]);

        if constexpr (Stage >= UFpFp6MulStage::PreAddMul) {
            UFp2Pack8x1x1::mul(t0, data.a_pack[p].c[0], data.b_pack[p].c[0]);
            UFp2Pack8x1x1::mul(t1, data.a_pack[p].c[1], data.b_pack[p].c[1]);
            UFp2Pack8x1x1::mul(t2, data.a_pack[p].c[2], data.b_pack[p].c[2]);
            UFp2Pack8x1x1::mul(w0, aa0, bb0);
            UFp2Pack8x1x1::mul(w1, aa1, bb1);
            UFp2Pack8x1x1::mul(w2, aa2, bb2);

            if constexpr (Stage >= UFpFp6MulStage::Full) {
                out.c[0] = w0;
                UFp2Pack8x1x1::sub(out.c[0], out.c[0], t1);
                UFp2Pack8x1x1::sub(out.c[0], out.c[0], t2);
                UFp2Pack8x1x1::mul_by_u_plus_1(t3, out.c[0]);
                UFp2Pack8x1x1::add(out.c[0], t3, t0);

                out.c[1] = w1;
                UFp2Pack8x1x1::sub(out.c[1], out.c[1], t0);
                UFp2Pack8x1x1::sub(out.c[1], out.c[1], t1);
                UFp2Pack8x1x1::mul_by_u_plus_1(t3, t2);
                UFp2Pack8x1x1::add(out.c[1], out.c[1], t3);

                out.c[2] = w2;
                UFp2Pack8x1x1::sub(out.c[2], out.c[2], t0);
                UFp2Pack8x1x1::sub(out.c[2], out.c[2], t2);
                UFp2Pack8x1x1::add(out.c[2], out.c[2], t1);
                data.out_wide[p] = out;
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

template <UFpFp6MulStage Stage>
void BM_Fp6MulUfpAccounting(benchmark::State& state)
{
    Fp6BenchData data(static_cast<std::size_t>(state.range(0)));

    for (auto _ : state) {
        ufp_fp6_mul_accounting<Stage>(data);
        if constexpr (Stage >= UFpFp6MulStage::Full) {
            consume_wide_fp6(data.out_wide);
        } else {
            benchmark::ClobberMemory();
        }
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs));
}

void BM_Fp6MulUfpWide(benchmark::State& state)
{
    Fp6BenchData data(static_cast<std::size_t>(state.range(0)));
    for (auto _ : state) {
        ufp_fp6_mul_accounting<UFpFp6MulStage::Full>(data);
        consume_wide_fp6(data.out_wide);
    }

    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) *
                            static_cast<int64_t>(data.n_packs));
}

BENCHMARK(BM_UFpAddChain)->Arg(256);
BENCHMARK_TEMPLATE(BM_UFpAddAccounting, UFpAddStage::RawAdd)->Arg(256);
BENCHMARK_TEMPLATE(BM_UFpAddAccounting, UFpAddStage::RawAddBias)->Arg(256);
BENCHMARK_TEMPLATE(BM_UFpAddAccounting, UFpAddStage::RawAddBiasCarry)->Arg(256);
BENCHMARK_TEMPLATE(BM_UFpAddAccounting, UFpAddStage::Full)->Arg(256);

BENCHMARK(BM_Fp6MulUfpWide)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulUfpAccounting, UFpFp6MulStage::PreAdd)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulUfpAccounting, UFpFp6MulStage::PreAddMul)->Arg(256);
BENCHMARK_TEMPLATE(BM_Fp6MulUfpAccounting, UFpFp6MulStage::Full)->Arg(256);
}

#endif
