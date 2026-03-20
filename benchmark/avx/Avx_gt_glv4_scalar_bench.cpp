#include <benchmark/benchmark.h>

#include "test/avx/avx_test_utils.h"

using namespace blst::avx;

namespace {

static void BM_GTScalarCyclotomicExp_Basic(benchmark::State& state)
{
    const blst::Fp12 base = blst_avx_test::make_scalar_cyclotomic_fp12_from_seed(42);
    const auto scalar_bytes = blst_avx_test::seeded_scalar_bytes_le(123);

    benchmark::DoNotOptimize(base);
    benchmark::DoNotOptimize(scalar_bytes);

    for (auto _ : state) {
        blst::Fp12 result = {};
        blst::Fp12::cyclotomic_exp_basic(result, base, scalar_bytes.data(), 256);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("GTScalarCyclotomicExp_Basic");
}

static void BM_GTScalarCyclotomicExp_GLV4(benchmark::State& state)
{
    const blst::Fp12 base = blst_avx_test::make_scalar_cyclotomic_fp12_from_seed(42);
    const auto scalar_limbs = blst_avx_test::seeded_scalar_limbs_le(123);

    benchmark::DoNotOptimize(base);
    benchmark::DoNotOptimize(scalar_limbs);

    for (auto _ : state) {
        blst::Fp12 result = {};
        blst::Fp12::cyclotomic_exp_glv4(result, base, scalar_limbs.data());
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("GTScalarCyclotomicExp_GLV4");
}

static void BM_GTScalarCyclotomicExp_GLV4_Constant(benchmark::State& state)
{
    const blst::Fp12 base = blst_avx_test::make_scalar_cyclotomic_fp12_from_seed(42);
    const auto scalar_limbs = blst_avx_test::seeded_scalar_limbs_le(123);
    const blst::GLV4_Exp props = blst::Fp12::glv4_compile(scalar_limbs.data());

    benchmark::DoNotOptimize(base);
    benchmark::DoNotOptimize(props);

    for (auto _ : state) {
        blst::Fp12 result = {};
        blst::Fp12::cyclotomic_exp_glv4(result, base, props);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("GTScalarCyclotomicExp_GLV4_Constant");
}

BENCHMARK(BM_GTScalarCyclotomicExp_Basic);
BENCHMARK(BM_GTScalarCyclotomicExp_GLV4);
BENCHMARK(BM_GTScalarCyclotomicExp_GLV4_Constant);

}  // namespace
