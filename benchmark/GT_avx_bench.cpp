#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#if defined(BLST_AVX512)

#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

#include "blst/curve/GT.h"
#include "blst/curve/GT_avx.h"
#include "blst/fields/Fr.h"

using namespace blst;

namespace
{
constexpr std::size_t kStride = 32;
constexpr std::size_t kNBits = 256;

u64 load_le_u64(const byte* in)
{
    u64 v = 0;
    for (std::size_t i = 0; i < sizeof(v); ++i) {
        v |= static_cast<u64>(in[i]) << (8 * i);
    }
    return v;
}

void load_scalar_limbs_256(limb_t out[4], const byte* scalar_32)
{
    for (std::size_t i = 0; i < 4; ++i) {
        out[i] = static_cast<limb_t>(load_le_u64(scalar_32 + i * sizeof(u64)));
    }
}

bool scalar_ge_group_order(const byte* scalar_32)
{
    const Fr& order = get_BLST_r();
    for (int limb = 3; limb >= 0; --limb) {
        const u64 s = load_le_u64(scalar_32 + limb * 8);
        const u64 r = order[limb];
        if (s > r) {
            return true;
        }
        if (s < r) {
            return false;
        }
    }
    return true;
}

void scalar_rshift_1(byte* scalar_32)
{
    byte carry = 0;
    for (int i = 31; i >= 0; --i) {
        const byte next_carry = static_cast<byte>(scalar_32[i] & 1u);
        scalar_32[i] = static_cast<byte>((scalar_32[i] >> 1) | (carry << 7));
        carry = next_carry;
    }
}

void random_reduced_scalar_32(byte* out, std::mt19937_64& gen)
{
    for (std::size_t i = 0; i < kStride; i += sizeof(u64)) {
        const u64 x = gen();
        for (std::size_t j = 0; j < sizeof(u64); ++j) {
            out[i + j] = static_cast<byte>((x >> (8 * j)) & 0xFFu);
        }
    }
    while (scalar_ge_group_order(out)) {
        scalar_rshift_1(out);
    }
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

void random_redc_fp2(Fp2& v, std::mt19937_64& gen)
{
    random_redc_fp(v[0], gen);
    random_redc_fp(v[1], gen);
}

void random_gt(GT& out, std::mt19937_64& gen)
{
    Fp12 seed{};
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            random_redc_fp2(seed[i][j], gen);
        }
    }
    Fp12::final_exp(out.data, seed);
}

struct GTBenchData
{
    explicit GTBenchData(const std::size_t packs)
        : n_packs(packs)
        , n_tuples(packs * 8)
        , a(n_tuples)
        , b(n_tuples)
        , out_scalar(n_tuples)
        , scalars(n_tuples * kStride)
        , scalar_props(n_tuples)
        , a_pack(n_packs)
        , b_pack(n_packs)
        , out_pack(n_packs)
        , pack_glv4(n_packs)
    {
        std::mt19937_64 gen(0xB17A5EED1234ULL + static_cast<u64>(packs));

        for (std::size_t i = 0; i < n_tuples; ++i) {
            random_gt(a[i], gen);
            random_gt(b[i], gen);
            random_reduced_scalar_32(&scalars[i * kStride], gen);

            limb_t scalar_limbs[4]{};
            load_scalar_limbs_256(scalar_limbs, &scalars[i * kStride]);
            scalar_props[i] = Fp12::glv4_compile(scalar_limbs);
        }

        for (std::size_t p = 0; p < n_packs; ++p) {
            GTPack8x1::from_aos(a_pack[p], &a[p * 8]);
            GTPack8x1::from_aos(b_pack[p], &b[p * 8]);

            const bool ok =
                GTPack8x1::compile_glv4(pack_glv4[p], &scalars[p * 8 * kStride], kStride, kNBits);
            BLST_ASSERT(ok);
        }
    }

    const std::size_t n_packs;
    const std::size_t n_tuples;
    std::vector<GT> a;
    std::vector<GT> b;
    std::vector<GT> out_scalar;
    std::vector<byte> scalars;
    std::vector<GLV4_Exp> scalar_props;
    std::vector<GTPack8x1> a_pack;
    std::vector<GTPack8x1> b_pack;
    std::vector<GTPack8x1> out_pack;
    std::vector<GTGlv4Pack8x1> pack_glv4;
};

void BM_GTAdd_Scalar(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);

    for (auto _ : state) {
        for (std::size_t i = 0; i < data.n_tuples; ++i) {
            GT::add(data.out_scalar[i], data.a[i], data.b[i]);
        }
        benchmark::DoNotOptimize(data.out_scalar.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTAdd_Scalar)->Arg(1)->Arg(8)->Arg(32);

void BM_GTAdd_Pack(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);

    for (auto _ : state) {
        GTPack8x1::add(
            std::span<GTPack8x1>(data.out_pack.data(), data.n_packs),
            std::span<const GTPack8x1>(data.a_pack.data(), data.n_packs),
            std::span<const GTPack8x1>(data.b_pack.data(), data.n_packs));
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTAdd_Pack)->Arg(1)->Arg(8)->Arg(32);

void BM_GTDoubling_Scalar(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);

    for (auto _ : state) {
        for (std::size_t i = 0; i < data.n_tuples; ++i) {
            GT::doubling(data.out_scalar[i], data.a[i]);
        }
        benchmark::DoNotOptimize(data.out_scalar.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTDoubling_Scalar)->Arg(1)->Arg(8)->Arg(32);

void BM_GTDoubling_Pack(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);

    for (auto _ : state) {
        GTPack8x1::doubling(
            std::span<GTPack8x1>(data.out_pack.data(), data.n_packs),
            std::span<const GTPack8x1>(data.a_pack.data(), data.n_packs));
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTDoubling_Pack)->Arg(1)->Arg(8)->Arg(32);

void BM_GTScalarMult_Scalar(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);

    for (auto _ : state) {
        for (std::size_t i = 0; i < data.n_tuples; ++i) {
            GT::mult(data.out_scalar[i], data.a[i], &data.scalars[i * kStride], kNBits);
        }
        benchmark::DoNotOptimize(data.out_scalar.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTScalarMult_Scalar)->Arg(1)->Arg(4);

void BM_GTScalarMult_Pack(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);

    for (auto _ : state) {
        GTPack8x1::mult(
            std::span<GTPack8x1>(data.out_pack.data(), data.n_packs),
            std::span<const GTPack8x1>(data.a_pack.data(), data.n_packs),
            std::span<const byte>(data.scalars.data(), data.scalars.size()),
            kStride,
            kNBits);
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTScalarMult_Pack)->Arg(1)->Arg(4);

void BM_GTScalarMult_ScalarPrecompiled(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);

    for (auto _ : state) {
        for (std::size_t i = 0; i < data.n_tuples; ++i) {
            Fp12::cyclotomic_exp_glv4(data.out_scalar[i].data, data.a[i].data, data.scalar_props[i]);
        }
        benchmark::DoNotOptimize(data.out_scalar.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTScalarMult_ScalarPrecompiled)->Arg(1)->Arg(4);

void BM_GTScalarMult_PackPrecompiled(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);

    for (auto _ : state) {
        GTPack8x1::mult(
            std::span<GTPack8x1>(data.out_pack.data(), data.n_packs),
            std::span<const GTPack8x1>(data.a_pack.data(), data.n_packs),
            std::span<const GTGlv4Pack8x1>(data.pack_glv4.data(), data.n_packs));
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTScalarMult_PackPrecompiled)->Arg(1)->Arg(4);

} // namespace

#endif // BLST_AVX512
