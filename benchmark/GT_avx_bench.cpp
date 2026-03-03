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
        , a_pack(n_packs)
        , b_pack(n_packs)
        , out_pack(n_packs)
    {
        std::mt19937_64 gen(0xB17A5EED1234ULL + static_cast<u64>(packs));

        for (std::size_t i = 0; i < n_tuples; ++i) {
            random_gt(a[i], gen);
            random_gt(b[i], gen);
            random_reduced_scalar_32(&scalars[i * kStride], gen);
        }

        for (std::size_t p = 0; p < n_packs; ++p) {
            GTPack8x1::from_aos(a_pack[p], &a[p * 8]);
            GTPack8x1::from_aos(b_pack[p], &b[p * 8]);
        }
    }

    const std::size_t n_packs;
    const std::size_t n_tuples;
    std::vector<GT> a;
    std::vector<GT> b;
    std::vector<GT> out_scalar;
    std::vector<byte> scalars;
    std::vector<GTPack8x1> a_pack;
    std::vector<GTPack8x1> b_pack;
    std::vector<GTPack8x1> out_pack;
};

GTPack8x1 gt_identity_pack_bench()
{
    static const GTPack8x1 kIdentity = []() {
        std::array<GT, 8> id{};
        for (std::size_t i = 0; i < id.size(); ++i) {
            id[i] = GT::identity();
        }
        GTPack8x1 out{};
        GTPack8x1::from_aos(out, id.data());
        return out;
    }();
    return kIdentity;
}

struct GTGlv4KernelPack
{
    std::array<Fp12Pack8x1, 16> table{};
    Fp12Pack8x1 init{};
};

struct GTGlv4DecodedPack
{
    struct StepMulSparse
    {
        u8 count = 0;
        std::array<u8, 8> idx{};
        std::array<u8, 8> lane_mask{};
    };

    std::size_t max_steps = 0;
    std::array<u8, 65> sqr_mask{};
    std::array<StepMulSparse, 65> mul_terms{};
};

void gt_glv4_build_table_init(GTGlv4KernelPack& kernel, const GTPack8x1& base, const GTGlv4Pack8x1& plan)
{
    std::array<Fp12Pack8x1, 4> basis{};
    basis[0] = base.c;
    Fp12Pack8x1::frobenius_map(basis[1], base.c, 1);
    Fp12Pack8x1::frobenius_map(basis[2], base.c, 2);
    Fp12Pack8x1::frobenius_map(basis[3], base.c, 3);

    for (std::size_t basis_idx = 0; basis_idx < basis.size(); ++basis_idx) {
        u8 conj_mask = 0;
        for (std::size_t lane = 0; lane < 8; ++lane) {
            const u8 sign_bit = static_cast<u8>((plan.props[lane].sign_mask >> basis_idx) & 1u);
            conj_mask |= static_cast<u8>(sign_bit << lane);
        }
        if (conj_mask != 0) {
            Fp12Pack8x1 conj{};
            Fp12Pack8x1::conjugate(conj, basis[basis_idx]);
            Fp12Pack8x1::select(basis[basis_idx], basis[basis_idx], conj, conj_mask);
        }
    }

    kernel.table[0] = gt_identity_pack_bench().c;
    kernel.table[1] = basis[0];
    kernel.table[2] = basis[1];
    kernel.table[4] = basis[2];
    kernel.table[8] = basis[3];

    Fp12Pack8x1::mul(kernel.table[3], kernel.table[1], kernel.table[2]);
    Fp12Pack8x1::mul(kernel.table[5], kernel.table[1], kernel.table[4]);
    Fp12Pack8x1::mul(kernel.table[6], kernel.table[2], kernel.table[4]);
    Fp12Pack8x1::mul(kernel.table[7], kernel.table[3], kernel.table[4]);
    Fp12Pack8x1::mul(kernel.table[9], kernel.table[1], kernel.table[8]);
    Fp12Pack8x1::mul(kernel.table[10], kernel.table[2], kernel.table[8]);
    Fp12Pack8x1::mul(kernel.table[11], kernel.table[3], kernel.table[8]);
    Fp12Pack8x1::mul(kernel.table[12], kernel.table[4], kernel.table[8]);
    Fp12Pack8x1::mul(kernel.table[13], kernel.table[5], kernel.table[8]);
    Fp12Pack8x1::mul(kernel.table[14], kernel.table[6], kernel.table[8]);
    Fp12Pack8x1::mul(kernel.table[15], kernel.table[7], kernel.table[8]);

    kernel.init = gt_identity_pack_bench().c;
    for (u8 idx = 1; idx < 16; ++idx) {
        u8 lane_mask = 0;
        for (std::size_t lane = 0; lane < 8; ++lane) {
            lane_mask |= static_cast<u8>((plan.init_idx[lane] == idx ? 1u : 0u) << lane);
        }
        if (lane_mask != 0) {
            Fp12Pack8x1::select(kernel.init, kernel.init, kernel.table[idx], lane_mask);
        }
    }
}

void gt_glv4_apply_loop(Fp12Pack8x1& out, const GTGlv4KernelPack& kernel, const GTGlv4Pack8x1& plan)
{
    Fp12Pack8x1 sqr_candidate{};
    Fp12Pack8x1 mul_candidate{};
    Fp12Pack8x1 mul_operand{};
    for (std::size_t step = 1; step <= plan.max_steps; ++step) {
        u8 sqr_mask = 0;
        std::array<u8, 16> mul_masks{};
        for (std::size_t lane = 0; lane < 8; ++lane) {
            const bool lane_active = (plan.init_idx[lane] != 0) && (step <= static_cast<std::size_t>(plan.total_steps[lane]));
            sqr_mask |= static_cast<u8>((lane_active ? 1u : 0u) << lane);
            if (lane_active) {
                const u8 idx = plan.mul_idx_at_step[lane][step];
                mul_masks[idx] |= static_cast<u8>(1u << lane);
            }
        }
        if (sqr_mask != 0) {
            Fp12Pack8x1::cyclotomic_sqr(sqr_candidate, out);
            Fp12Pack8x1::select(out, out, sqr_candidate, sqr_mask);
        }

        u8 nonzero_mul_mask = 0;
        for (u8 idx = 1; idx < 16; ++idx) {
            nonzero_mul_mask |= mul_masks[idx];
        }
        if (nonzero_mul_mask != 0) {
            mul_operand = kernel.table[0];
            for (u8 idx = 1; idx < 16; ++idx) {
                if (mul_masks[idx] != 0) {
                    Fp12Pack8x1::select(mul_operand, mul_operand, kernel.table[idx], mul_masks[idx]);
                }
            }
            Fp12Pack8x1::mul(mul_candidate, out, mul_operand);
            out = mul_candidate;
        }
    }
}

void gt_glv4_predecode(GTGlv4DecodedPack& decoded, const GTGlv4Pack8x1& plan)
{
    decoded.max_steps = static_cast<std::size_t>(plan.max_steps);
    for (std::size_t step = 1; step <= decoded.max_steps; ++step) {
        u8 sqr_mask = 0;
        auto& mul_terms = decoded.mul_terms[step];
        mul_terms.count = 0;
        for (std::size_t lane = 0; lane < 8; ++lane) {
            const bool lane_active = (plan.init_idx[lane] != 0) && (step <= static_cast<std::size_t>(plan.total_steps[lane]));
            sqr_mask |= static_cast<u8>((lane_active ? 1u : 0u) << lane);
            if (lane_active) {
                const u8 idx = plan.mul_idx_at_step[lane][step];
                if (idx != 0) {
                    const u8 lane_bit = static_cast<u8>(1u << lane);
                    bool merged = false;
                    for (u8 i = 0; i < mul_terms.count; ++i) {
                        if (mul_terms.idx[i] == idx) {
                            mul_terms.lane_mask[i] |= lane_bit;
                            merged = true;
                            break;
                        }
                    }
                    if (!merged) {
                        BLST_ASSERT(mul_terms.count < 8);
                        mul_terms.idx[mul_terms.count] = idx;
                        mul_terms.lane_mask[mul_terms.count] = lane_bit;
                        ++mul_terms.count;
                    }
                }
            }
        }

        decoded.sqr_mask[step] = sqr_mask;
    }
}

void gt_glv4_apply_loop_predecoded(Fp12Pack8x1& out, const GTGlv4KernelPack& kernel, const GTGlv4DecodedPack& decoded)
{
    Fp12Pack8x1 sqr_candidate{};
    Fp12Pack8x1 mul_candidate{};

    for (std::size_t step = 1; step <= decoded.max_steps; ++step) {
        const u8 sqr_mask = decoded.sqr_mask[step];
        if (sqr_mask != 0) {
            Fp12Pack8x1::cyclotomic_sqr(sqr_candidate, out);
            Fp12Pack8x1::select(out, out, sqr_candidate, sqr_mask);
        }

        const auto& mul_terms = decoded.mul_terms[step];
        if (mul_terms.count != 0) {
            Fp12Pack8x1 mul_operand = kernel.table[0];
            for (u8 i = 0; i < mul_terms.count; ++i) {
                Fp12Pack8x1::select(
                    mul_operand, mul_operand, kernel.table[mul_terms.idx[i]], mul_terms.lane_mask[i]);
            }
            Fp12Pack8x1::mul(mul_candidate, out, mul_operand);
            out = mul_candidate;
        }
    }
}

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

void BM_GTScalarMult_PackGLV4_PrecompiledExec(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);
    std::vector<GTGlv4Pack8x1> plans(data.n_packs);
    const bool ok = GTPack8x1::compile_glv4(
        std::span<GTGlv4Pack8x1>(plans.data(), plans.size()),
        std::span<const byte>(data.scalars.data(), data.scalars.size()),
        kStride,
        kNBits);
    if (!ok) {
        state.SkipWithError("compile_glv4 failed");
        return;
    }

    for (auto _ : state) {
        GTPack8x1::mult(
            std::span<GTPack8x1>(data.out_pack.data(), data.n_packs),
            std::span<const GTPack8x1>(data.a_pack.data(), data.n_packs),
            std::span<const GTGlv4Pack8x1>(plans.data(), plans.size()));
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTScalarMult_PackGLV4_PrecompiledExec)->Arg(1)->Arg(4);

void BM_GTScalarMult_PackGLV4_TableInitOnly(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);
    std::vector<GTGlv4Pack8x1> plans(data.n_packs);
    std::vector<GTGlv4KernelPack> kernels(data.n_packs);

    const bool ok = GTPack8x1::compile_glv4(
        std::span<GTGlv4Pack8x1>(plans.data(), plans.size()),
        std::span<const byte>(data.scalars.data(), data.scalars.size()),
        kStride,
        kNBits);
    if (!ok) {
        state.SkipWithError("compile_glv4 failed");
        return;
    }

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            gt_glv4_build_table_init(kernels[p], data.a_pack[p], plans[p]);
        }
        benchmark::DoNotOptimize(kernels.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTScalarMult_PackGLV4_TableInitOnly)->Arg(1)->Arg(4);

void BM_GTScalarMult_PackGLV4_LoopOnly(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);
    std::vector<GTGlv4Pack8x1> plans(data.n_packs);
    std::vector<GTGlv4KernelPack> kernels(data.n_packs);

    const bool ok = GTPack8x1::compile_glv4(
        std::span<GTGlv4Pack8x1>(plans.data(), plans.size()),
        std::span<const byte>(data.scalars.data(), data.scalars.size()),
        kStride,
        kNBits);
    if (!ok) {
        state.SkipWithError("compile_glv4 failed");
        return;
    }
    for (std::size_t p = 0; p < data.n_packs; ++p) {
        gt_glv4_build_table_init(kernels[p], data.a_pack[p], plans[p]);
    }

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            data.out_pack[p].c = kernels[p].init;
            gt_glv4_apply_loop(data.out_pack[p].c, kernels[p], plans[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTScalarMult_PackGLV4_LoopOnly)->Arg(1)->Arg(4);

void BM_GTScalarMult_PackGLV4_LoopOnly_Predecoded(benchmark::State& state)
{
    const std::size_t packs = static_cast<std::size_t>(state.range(0));
    GTBenchData data(packs);
    std::vector<GTGlv4Pack8x1> plans(data.n_packs);
    std::vector<GTGlv4KernelPack> kernels(data.n_packs);
    std::vector<GTGlv4DecodedPack> decoded(data.n_packs);

    const bool ok = GTPack8x1::compile_glv4(
        std::span<GTGlv4Pack8x1>(plans.data(), plans.size()),
        std::span<const byte>(data.scalars.data(), data.scalars.size()),
        kStride,
        kNBits);
    if (!ok) {
        state.SkipWithError("compile_glv4 failed");
        return;
    }
    for (std::size_t p = 0; p < data.n_packs; ++p) {
        gt_glv4_build_table_init(kernels[p], data.a_pack[p], plans[p]);
        gt_glv4_predecode(decoded[p], plans[p]);
    }

    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs; ++p) {
            data.out_pack[p].c = kernels[p].init;
            gt_glv4_apply_loop_predecoded(data.out_pack[p].c, kernels[p], decoded[p]);
        }
        benchmark::DoNotOptimize(data.out_pack.data());
    }
    state.SetItemsProcessed(state.iterations() * data.n_tuples);
}
BENCHMARK(BM_GTScalarMult_PackGLV4_LoopOnly_Predecoded)->Arg(1)->Arg(4);

} // namespace

#endif // BLST_AVX512
