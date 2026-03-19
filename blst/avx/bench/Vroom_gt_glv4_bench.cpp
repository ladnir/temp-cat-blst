#include <benchmark/benchmark.h>

#include "blst/avx/gt/glv4.h"
#include "blst/avx/tests/avx_test_utils.h"

namespace {

using blst_vroom_test::VroomFp12Ring;
using blst_vroom_test::VroomRing;

static void BM_VroomGTCyclotomicExp_Basic_52(benchmark::State& state)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const auto base = blst_vroom_test::make_vroom_cyclotomic_fp12_from_seed(ring, 42);
    const auto scalar_limbs = blst_vroom_test::seeded_scalar_limbs_le(123);

    benchmark::DoNotOptimize(base);
    benchmark::DoNotOptimize(scalar_limbs);

    for (auto _ : state) {
        const auto result = vroom_cyclotomic_exp_basic(base, scalar_limbs.data(), fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("VroomGTCyclotomicExp_Basic_52");
}

static void BM_VroomGTCyclotomicExp_GLV4_52(benchmark::State& state)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const auto base = blst_vroom_test::make_vroom_cyclotomic_fp12_from_seed(ring, 42);
    const auto scalar_limbs = blst_vroom_test::seeded_scalar_limbs_le(123);

    benchmark::DoNotOptimize(base);
    benchmark::DoNotOptimize(scalar_limbs);

    for (auto _ : state) {
        const auto result = vroom_cyclotomic_exp_glv4(base, scalar_limbs.data(), fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("VroomGTCyclotomicExp_GLV4_52");
}

static void BM_VroomGTCyclotomicExp_GLV4_Constant_52(benchmark::State& state)
{
    const VroomRing ring(blst_vroom_test::modulus_bigint());
    const VroomFp12Ring fp12ring(ring);
    const auto base = blst_vroom_test::make_vroom_cyclotomic_fp12_from_seed(ring, 42);
    const auto scalar_limbs = blst_vroom_test::seeded_scalar_limbs_le(123);
    const VroomGLV4Exp props = vroom_glv4_compile(scalar_limbs.data());

    benchmark::DoNotOptimize(base);
    benchmark::DoNotOptimize(props);

    for (auto _ : state) {
        const auto result = vroom_cyclotomic_exp_glv4(base, props, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("VroomGTCyclotomicExp_GLV4_Constant_52");
}

BENCHMARK(BM_VroomGTCyclotomicExp_Basic_52);
BENCHMARK(BM_VroomGTCyclotomicExp_GLV4_52);
BENCHMARK(BM_VroomGTCyclotomicExp_GLV4_Constant_52);

}  // namespace
