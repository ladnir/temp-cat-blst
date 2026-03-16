#include <benchmark/benchmark.h>

#include <array>

#include "blst/crns_batch/VroomLiteral.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fp12_avx.h"
#include "benchmark/internal/Fp_avx_lufp_v2_impl.h"
#include "test/InsecureRng.h"

#if defined(BLST_CRNS_BATCH_VROOM_VENDOR) && (defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE))

using namespace blst;

namespace
{
void rand_redc_fp(Fp& v, InsecureRng& rng)
{
    for (u64 i = 0; i < v.num_limb(); ++i) {
        v[i] = rng.next_u64();
    }
    while (Fp::is_redc(v) == 0) {
        Fp::rshift(v, v, 1);
    }
}

void rand_redc_fp2(Fp2& v, InsecureRng& rng)
{
    rand_redc_fp(v[0], rng);
    rand_redc_fp(v[1], rng);
}

void rand_redc_fp12(Fp12& v, InsecureRng& rng)
{
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            rand_redc_fp2(v[i][j], rng);
        }
    }
}

bool fp12_equal(const Fp12& a, const Fp12& b)
{
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 2; ++k) {
                if (!Fp::equal(a[i][j][k], b[i][j][k])) {
                    return false;
                }
            }
        }
    }
    return true;
}

struct CompareFp12Data
{
    CompareFp12Data()
    {
        InsecureRng rng(0x78563412u);
        for (size_t i = 0; i < a_scalar.size(); ++i) {
            rand_redc_fp12(a_scalar[i], rng);
            rand_redc_fp12(b_scalar[i], rng);
            a_vroom[i] = blst::crns_batch::vroom_literal::fp12_from_blst(a_scalar[i]);
            b_vroom[i] = blst::crns_batch::vroom_literal::fp12_from_blst(b_scalar[i]);
        }
        Fp12Pack8x1::from_aos(a_pack, a_scalar.data());
        Fp12Pack8x1::from_aos(b_pack, b_scalar.data());
    }

    std::array<Fp12, 8> a_scalar{};
    std::array<Fp12, 8> b_scalar{};
    std::array<Fp12, 8> out_scalar{};
    std::array<Fp12, 8> out_pack_aos{};
    std::array<Fp12, 8> out_vroom{};
    std::array<blst::crns_batch::vroom_literal::VroomFp12Element, 8> a_vroom{};
    std::array<blst::crns_batch::vroom_literal::VroomFp12Element, 8> b_vroom{};
    std::array<blst::crns_batch::vroom_literal::VroomFp12Element, 8> out_vroom_native{};
    Fp12Pack8x1 a_pack{};
    Fp12Pack8x1 b_pack{};
    Fp12Pack8x1 out_pack{};
};

void validate_fp12_compare_paths()
{
    static const bool validated = []() {
        CompareFp12Data data{};
        UFp12Pack8x1 wide{};

        for (size_t i = 0; i < data.a_scalar.size(); ++i) {
            Fp12::mul_scalar(data.out_scalar[i], data.a_scalar[i], data.b_scalar[i]);
        }

        blst::bench_lufp_v2::mul_fp12_pack_lufpv2_wide(wide, data.a_pack, data.b_pack);
        UFp12Pack8x1::redc(data.out_pack, wide);
        Fp12Pack8x1::to_aos(data.out_pack_aos.data(), data.out_pack);

        for (size_t i = 0; i < data.a_vroom.size(); ++i) {
            data.out_vroom[i] = blst::crns_batch::vroom_literal::fp12_to_blst(
                blst::crns_batch::vroom_literal::fp12_mul(data.a_vroom[i], data.b_vroom[i]));
        }

        for (size_t i = 0; i < data.out_scalar.size(); ++i) {
            if (!fp12_equal(data.out_scalar[i], data.out_pack_aos[i])) {
                return false;
            }
            if (!fp12_equal(data.out_scalar[i], data.out_vroom[i])) {
                return false;
            }
        }
        return true;
    }();

    if (!validated) {
        std::abort();
    }
    benchmark::DoNotOptimize(validated);
}

static void BM_Fp12Compare_ScalarAsmX8(benchmark::State& state)
{
    validate_fp12_compare_paths();
    CompareFp12Data data{};

    for (auto _ : state) {
        for (size_t i = 0; i < data.a_scalar.size(); ++i) {
            Fp12::mul_scalar(data.out_scalar[i], data.a_scalar[i], data.b_scalar[i]);
        }
        benchmark::DoNotOptimize(data.out_scalar);
    }
    state.SetItemsProcessed(state.iterations() * 8);
}

static void BM_Fp12Compare_LUFpV2WideX8(benchmark::State& state)
{
    validate_fp12_compare_paths();
    CompareFp12Data data{};

    for (auto _ : state) {
        UFp12Pack8x1 wide{};
        blst::bench_lufp_v2::mul_fp12_pack_lufpv2_wide(wide, data.a_pack, data.b_pack);
        UFp12Pack8x1::redc(data.out_pack, wide);
        benchmark::DoNotOptimize(data.out_pack);
    }
    state.SetItemsProcessed(state.iterations() * 8);
}

static void BM_Fp12Compare_VroomLiteralX8(benchmark::State& state)
{
    validate_fp12_compare_paths();
    CompareFp12Data data{};

    for (auto _ : state) {
        for (size_t i = 0; i < data.a_vroom.size(); ++i) {
            data.out_vroom_native[i] =
                blst::crns_batch::vroom_literal::fp12_mul(data.a_vroom[i], data.b_vroom[i]);
        }
        benchmark::DoNotOptimize(data.out_vroom_native);
    }
    state.SetItemsProcessed(state.iterations() * 8);
}

BENCHMARK(BM_Fp12Compare_ScalarAsmX8);
BENCHMARK(BM_Fp12Compare_LUFpV2WideX8);
BENCHMARK(BM_Fp12Compare_VroomLiteralX8);
} // namespace

#endif
