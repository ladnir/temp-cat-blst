#include "blst/common/defines.h"
#include <benchmark/benchmark.h>

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

#include <array>

#include "benchmark/internal/Fp_avx_bench_kernels.h"

using namespace blst;
using namespace blst::bench_internal;

namespace
{
std::array<u64, NWORDS> make_scaled_mod(const u64 scale)
{
    std::array<u64, NWORDS> out{};
    u64 carry = 0;
    for (std::size_t i = 0; i < NWORDS; ++i) {
        const u64 accum = P48[i] * scale + carry;
        out[i] = accum & BMASK;
        carry = accum >> BRADIX;
    }
    return out;
}

std::array<cc512i, NWORDS> broadcast_mod(const std::array<u64, NWORDS>& mod)
{
    std::array<cc512i, NWORDS> out{};
    for (std::size_t i = 0; i < NWORDS; ++i) {
        out[i] = VSET1(static_cast<i64>(mod[i]));
    }
    return out;
}

const std::array<cc512i, NWORDS> kP1 = broadcast_mod(make_scaled_mod(1));
const std::array<cc512i, NWORDS> kP2 = broadcast_mod(make_scaled_mod(2));
const std::array<cc512i, NWORDS> kP4 = broadcast_mod(make_scaled_mod(4));

struct BenchGtEq
{
    u8 gt;
    u8 eq;
};

BLST_FORCEINLINE auto mask_eq_epi64(const cc512i a, const cc512i b) -> u8
{
#if BLST_AVX_NATIVE
    return static_cast<u8>(_mm512_cmpeq_epi64_mask(to_native_512(a), to_native_512(b)));
#else
    std::array<u64, 8> lhs{};
    std::array<u64, 8> rhs{};
    _cc512_storeu_si512(lhs.data(), a);
    _cc512_storeu_si512(rhs.data(), b);
    u8 mask = 0;
    for (u8 lane = 0; lane < 8; ++lane) {
        if (lhs[lane] == rhs[lane]) {
            mask |= static_cast<u8>(1u << lane);
        }
    }
    return mask;
#endif
}

BLST_FORCEINLINE auto mask_gt_epi64(const cc512i a, const cc512i b) -> u8
{
#if BLST_AVX_NATIVE
    return static_cast<u8>(_mm512_cmpgt_epi64_mask(to_native_512(a), to_native_512(b)));
#else
    std::array<u64, 8> lhs{};
    std::array<u64, 8> rhs{};
    _cc512_storeu_si512(lhs.data(), a);
    _cc512_storeu_si512(rhs.data(), b);
    u8 mask = 0;
    for (u8 lane = 0; lane < 8; ++lane) {
        if (static_cast<i64>(lhs[lane]) > static_cast<i64>(rhs[lane])) {
            mask |= static_cast<u8>(1u << lane);
        }
    }
    return mask;
#endif
}

BLST_FORCEINLINE auto gt_eq_pair(const cc512i a, const cc512i b) -> BenchGtEq
{
    return BenchGtEq{mask_gt_epi64(a, b), mask_eq_epi64(a, b)};
}

BLST_FORCEINLINE auto gt_eq_combine_hi_lo(const BenchGtEq hi, const BenchGtEq lo) -> BenchGtEq
{
#if BLST_AVX_NATIVE
    const __mmask8 hi_gt = static_cast<__mmask8>(hi.gt);
    const __mmask8 hi_eq = static_cast<__mmask8>(hi.eq);
    const __mmask8 lo_gt = static_cast<__mmask8>(lo.gt);
    const __mmask8 lo_eq = static_cast<__mmask8>(lo.eq);
    return BenchGtEq{
        static_cast<u8>(_kor_mask8(hi_gt, _kand_mask8(hi_eq, lo_gt))),
        static_cast<u8>(_kand_mask8(hi_eq, lo_eq)),
    };
#else
    return BenchGtEq{
        static_cast<u8>(hi.gt | (hi.eq & lo.gt)),
        static_cast<u8>(hi.eq & lo.eq),
    };
#endif
}

template <u64 Scale>
BLST_FORCEINLINE auto high_ge_mod_mask_prefix_lex(const cc512i* limbs) -> u8
{
    const auto& mod_vec = Scale == 1 ? kP1 : (Scale == 2 ? kP2 : kP4);
    const BenchGtEq c7 = gt_eq_pair(limbs[NWORDS + 7], mod_vec[7]);
    const BenchGtEq c6 = gt_eq_pair(limbs[NWORDS + 6], mod_vec[6]);
    const BenchGtEq c5 = gt_eq_pair(limbs[NWORDS + 5], mod_vec[5]);
    const BenchGtEq c4 = gt_eq_pair(limbs[NWORDS + 4], mod_vec[4]);
    const BenchGtEq c3 = gt_eq_pair(limbs[NWORDS + 3], mod_vec[3]);
    const BenchGtEq c2 = gt_eq_pair(limbs[NWORDS + 2], mod_vec[2]);
    const BenchGtEq c1 = gt_eq_pair(limbs[NWORDS + 1], mod_vec[1]);
    const BenchGtEq c0 = gt_eq_pair(limbs[NWORDS + 0], mod_vec[0]);

    const BenchGtEq c76 = gt_eq_combine_hi_lo(c7, c6);
    const BenchGtEq c54 = gt_eq_combine_hi_lo(c5, c4);
    const BenchGtEq c32 = gt_eq_combine_hi_lo(c3, c2);
    const BenchGtEq c10 = gt_eq_combine_hi_lo(c1, c0);
    const BenchGtEq c7654 = gt_eq_combine_hi_lo(c76, c54);
    const BenchGtEq c3210 = gt_eq_combine_hi_lo(c32, c10);
    const BenchGtEq all = gt_eq_combine_hi_lo(c7654, c3210);

#if BLST_AVX_NATIVE
    return static_cast<u8>(_kor_mask8(static_cast<__mmask8>(all.gt), static_cast<__mmask8>(all.eq)));
#else
    return static_cast<u8>(all.gt | all.eq);
#endif
}

template <u64 Scale>
BLST_FORCEINLINE void high_sub_mod_masked(cc512i* limbs, const u8 lane_mask)
{
    const auto& mod_vec = Scale == 1 ? kP1 : (Scale == 2 ? kP2 : kP4);
    for (std::size_t i = 0; i < NWORDS; ++i) {
        limbs[NWORDS + i] = VMSUB(limbs[NWORDS + i], lane_mask, limbs[NWORDS + i], mod_vec[i]);
    }
}

BLST_FORCEINLINE void signed_carry_range(cc512i* limbs, const std::size_t begin, const std::size_t end)
{
    const cc512i bmask = VSET1(BMASK);
    for (std::size_t i = begin; i < end; ++i) {
        limbs[i + 1] = VADD(limbs[i + 1], VSRA(limbs[i], BRADIX));
        limbs[i] = VAND(limbs[i], bmask);
    }
}

struct RenormBenchData
{
    explicit RenormBenchData(const std::size_t n_packs)
        : n_packs_(n_packs)
        , in(n_packs_)
        , out(n_packs_)
    {
        std::mt19937_64 gen(0x6C7546708F1A2D3EULL + static_cast<u64>(n_packs_));
        for (std::size_t i = 0; i < n_packs_; ++i) {
            std::array<Fp, 8> lane_a{};
            std::array<Fp, 8> lane_b{};
            FpPack8x1 pack_a{}, pack_b{};
            UFpPack8x1 wide_a{}, wide_b{};
            for (std::size_t lane = 0; lane < 8; ++lane) {
                random_redc_fp(lane_a[lane], gen);
                random_redc_fp(lane_b[lane], gen);
            }
            FpPack8x1::from_aos(pack_a, lane_a.data());
            FpPack8x1::from_aos(pack_b, lane_b.data());
            FpPack8x1::mul(wide_a, pack_a, pack_a);
            FpPack8x1::mul(wide_b, pack_b, pack_b);
            LUFpPack8x1::sub(in[i], wide_a, wide_b);
        }
    }

    const std::size_t n_packs_;
    std::vector<LUFpPack8x1> in;
    std::vector<UFpPack8x1> out;
};

void BM_LUFpNormalizeHiNeg1Pos2(benchmark::State& state)
{
    RenormBenchData data(static_cast<std::size_t>(state.range(0)));
    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            LUFpPack8x1::normalize<LUFpNormalizeBound::HiNeg1Pos2>(data.out[p], data.in[p]);
        }
        consume_wide(data.out);
    }
}

void BM_LUFpCarryChain(benchmark::State& state)
{
    constexpr std::size_t kRounds = 128;
    RenormBenchData data(static_cast<std::size_t>(state.range(0)));
    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            cc512i limbs[2 * NWORDS]{};
            for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
                limbs[i] = data.in[p][i];
            }
            for (std::size_t round = 0; round < kRounds; ++round) {
                signed_carry_range(limbs, 0, (2 * NWORDS) - 1);
            }
            for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
                data.out[p][i] = limbs[i];
            }
        }
        consume_wide(data.out);
    }
}

template <std::size_t Steps>
void BM_LUFpRepairPrefixLex(benchmark::State& state)
{
    constexpr std::size_t kRounds = 256;
    RenormBenchData data(static_cast<std::size_t>(state.range(0)));
    for (auto _ : state) {
        for (std::size_t p = 0; p < data.n_packs_; ++p) {
            cc512i limbs[2 * NWORDS]{};
            for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
                limbs[i] = data.in[p][i];
            }
            for (std::size_t round = 0; round < kRounds; ++round) {
                signed_carry_range(limbs, 0, (2 * NWORDS) - 1);
                if constexpr (Steps >= 1) {
                    const u8 ge1 = high_ge_mod_mask_prefix_lex<1>(limbs);
                    high_sub_mod_masked<1>(limbs, ge1);
                    signed_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);
                }
                if constexpr (Steps >= 2) {
                    const u8 ge2 = high_ge_mod_mask_prefix_lex<2>(limbs);
                    high_sub_mod_masked<2>(limbs, ge2);
                    signed_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);
                }
                if constexpr (Steps >= 4) {
                    const u8 ge4 = high_ge_mod_mask_prefix_lex<4>(limbs);
                    high_sub_mod_masked<4>(limbs, ge4);
                    signed_carry_range(limbs, NWORDS, (2 * NWORDS) - 1);
                }
            }
            for (std::size_t i = 0; i < 2 * NWORDS; ++i) {
                data.out[p][i] = limbs[i];
            }
        }
        consume_wide(data.out);
    }
}

BENCHMARK(BM_LUFpNormalizeHiNeg1Pos2)->Arg(256);
BENCHMARK(BM_LUFpCarryChain)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpRepairPrefixLex, 1)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpRepairPrefixLex, 2)->Arg(256);
BENCHMARK_TEMPLATE(BM_LUFpRepairPrefixLex, 4)->Arg(256);
}

#endif
