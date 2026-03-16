#include <benchmark/benchmark.h>

#include <array>
#include <random>

#include "blst/crns_batch/VroomLiteral.h"
#include "blst/crns_batch/vroom_port/src/bounded_ring.hpp"
#include "blst/crns_batch/vroom_port/src/fp12.hpp"
#include "test/InsecureRng.h"

#if defined(BLST_CRNS_BATCH_VROOM_VENDOR)

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

template<class Ring>
auto generate_random_fp12(const Ring& ring, uint64_t seed) -> typename FP12Ring<Ring>::StandardElement
{
    typename FP12Ring<Ring>::StandardElement result{};
    const BigInt p(blst::crns_batch::vroom_literal::kBls12381ModulusDec, 10);
    volatile uint64_t s = seed;
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(s);

    for (size_t i = 0; i < result.size(); ++i) {
        BigInt val(rng.get_z_range(p.get_mpz()));
        result[i] = ring.from_bigint(val);
    }

    return result;
}

static void BM_VroomVendorFpMulPacked(benchmark::State& state)
{
    InsecureRng rng(0x42424242u);
    std::array<Fp, 8> a{};
    std::array<Fp, 8> b{};
    for (size_t lane = 0; lane < 8; ++lane) {
        rand_redc_fp(a[lane], rng);
        rand_redc_fp(b[lane], rng);
    }

    blst::crns_batch::vroom_literal::VroomFpPack8 pa{};
    blst::crns_batch::vroom_literal::VroomFpPack8 pb{};
    blst::crns_batch::vroom_literal::from_aos(pa, a.data());
    blst::crns_batch::vroom_literal::from_aos(pb, b.data());

    for (auto _ : state) {
        auto pr = blst::crns_batch::vroom_literal::mul(pa, pb);
        benchmark::DoNotOptimize(pr);
    }
    state.SetItemsProcessed(state.iterations() * 8);
}

static void BM_VroomVendorBatchModMul8Direct(benchmark::State& state)
{
    const auto& ring = blst::crns_batch::vroom_literal::ring();
    InsecureRng rng(0x51525354u);
    std::array<blst::crns_batch::vroom_literal::VroomFpElement, 8> a{};
    std::array<blst::crns_batch::vroom_literal::VroomFpElement, 8> b{};

    for (size_t lane = 0; lane < 8; ++lane) {
        Fp af{};
        Fp bf{};
        rand_redc_fp(af, rng);
        rand_redc_fp(bf, rng);
        a[lane] = ring.from_bigint(blst::crns_batch::vroom_literal::fp_to_bigint(af));
        b[lane] = ring.from_bigint(blst::crns_batch::vroom_literal::fp_to_bigint(bf));
    }

    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        auto result = ring.batch_modmul<8>(a, b);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations() * 8);
}

static void BM_VroomVendorFP12MulDirect(benchmark::State& state)
{
    const auto& ring = blst::crns_batch::vroom_literal::ring();
    const auto& fp12ring = blst::crns_batch::vroom_literal::fp12_ring();

    volatile uint64_t seed_a = 42;
    volatile uint64_t seed_b = 123;
    auto a = generate_random_fp12(ring, seed_a);
    auto b = generate_random_fp12(ring, seed_b);

    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_mul(a, b, ring);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_VroomVendorFP12MulWrapped(benchmark::State& state)
{
    InsecureRng rng(0x12344321u);
    Fp12 a{};
    Fp12 b{};
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            rand_redc_fp(a[i][j][0], rng);
            rand_redc_fp(a[i][j][1], rng);
            rand_redc_fp(b[i][j][0], rng);
            rand_redc_fp(b[i][j][1], rng);
        }
    }

    const auto va = blst::crns_batch::vroom_literal::fp12_from_blst(a);
    const auto vb = blst::crns_batch::vroom_literal::fp12_from_blst(b);

    for (auto _ : state) {
        auto result = blst::crns_batch::vroom_literal::fp12_mul(va, vb);
        benchmark::DoNotOptimize(result);
    }
}

BENCHMARK(BM_VroomVendorFpMulPacked);
BENCHMARK(BM_VroomVendorBatchModMul8Direct);
BENCHMARK(BM_VroomVendorFP12MulDirect);
BENCHMARK(BM_VroomVendorFP12MulWrapped);
} // namespace

#endif
