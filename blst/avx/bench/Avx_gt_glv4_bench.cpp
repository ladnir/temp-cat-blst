#include <benchmark/benchmark.h>

#include "blst/avx/gt/glv4.h"
#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_avx_test::AvxFp12Ring;
using blst_avx_test::AvxRing;

static void BM_AvxGTCyclotomicExp_Basic_52(benchmark::State& state)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());
    const AvxFp12Ring fp12ring(ring);
    const auto base = blst_avx_test::make_avx_cyclotomic_fp12_from_seed(ring, 42);
    const auto scalar_limbs = blst_avx_test::seeded_scalar_limbs_le(123);

    benchmark::DoNotOptimize(base);
    benchmark::DoNotOptimize(scalar_limbs);

    for (auto _ : state) {
        const auto result = avx_cyclotomic_exp_basic(base, scalar_limbs.data(), fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("AvxGTCyclotomicExp_Basic_52");
}

static void BM_AvxGTCyclotomicExp_GLV4_52(benchmark::State& state)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());
    const AvxFp12Ring fp12ring(ring);
    const auto base = blst_avx_test::make_avx_cyclotomic_fp12_from_seed(ring, 42);
    const auto scalar_limbs = blst_avx_test::seeded_scalar_limbs_le(123);

    benchmark::DoNotOptimize(base);
    benchmark::DoNotOptimize(scalar_limbs);

    for (auto _ : state) {
        const auto result = avx_cyclotomic_exp_glv4(base, scalar_limbs.data(), fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("AvxGTCyclotomicExp_GLV4_52");
}

static void BM_AvxGTCyclotomicExp_GLV4_Constant_52(benchmark::State& state)
{
    const AvxRing ring(blst_avx_test::modulus_bigint());
    const AvxFp12Ring fp12ring(ring);
    const auto base = blst_avx_test::make_avx_cyclotomic_fp12_from_seed(ring, 42);
    const auto scalar_limbs = blst_avx_test::seeded_scalar_limbs_le(123);
    const AvxGLV4Exp props = avx_glv4_compile(scalar_limbs.data());

    benchmark::DoNotOptimize(base);
    benchmark::DoNotOptimize(props);

    for (auto _ : state) {
        const auto result = avx_cyclotomic_exp_glv4(base, props, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("AvxGTCyclotomicExp_GLV4_Constant_52");
}

BENCHMARK(BM_AvxGTCyclotomicExp_Basic_52);
BENCHMARK(BM_AvxGTCyclotomicExp_GLV4_52);
BENCHMARK(BM_AvxGTCyclotomicExp_GLV4_Constant_52);

}  // namespace
