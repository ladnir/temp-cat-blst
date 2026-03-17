#include <benchmark/benchmark.h>
#include <random>
#include <array>
#include <cstdint>
#include "fp12.hpp"
#include "fp2.hpp"
#include "bounded_ring.hpp"
#include "final_exponentiation.hpp"
#include "miller.hpp"
#include "inversion.hpp"
#include "conversion_inversion.hpp"
#include "pairing.hpp"
#include "ec.hpp"
#include "scalar_mult.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include "../test_data/test_values.hpp"

// G1 and G2 Generator values for Miller loop benchmarks
const char* g1_gen_x = "3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507";
const char* g1_gen_y = "1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569";
const char* g2_gen_x_c0 = "352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160";
const char* g2_gen_x_c1 = "3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758";
const char* g2_gen_y_c0 = "1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905";
const char* g2_gen_y_c1 = "927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582";

// BLS12-381 scalar field modulus r
const char* bls12_381_scalar_modulus_hex = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";

// Helper function to convert BigInt to byte array (little-endian)
static void bigint_to_bytes_le(uint8_t *bytes, const BigInt &value, size_t len) {
    BigInt temp = value;
    for (size_t i = 0; i < len; i++) {
        bytes[i] = static_cast<uint8_t>(temp.to_ulong() & 0xff);
        temp = temp >> 8;
    }
}

// Helper function to generate a random FP12 element
// Uses volatile copy of seed to prevent compiler optimization of the random generation
template<class Ring>
typename FP12Ring<Ring>::StandardElement generate_random_fp12(const Ring &ring, uint64_t seed) {
    typename FP12Ring<Ring>::StandardElement result;
    BigInt p = BigInt(modulus, 10);
    volatile uint64_t s = seed;
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(s);
    
    for (int i = 0; i < 12; i++) {
        // Generate random value less than modulus
        BigInt val(rng.get_z_range(p.get_mpz()));
        result[i] = ring.from_bigint(val);
    }
    
    return result;
}

static void BM_FP12_Mul(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    
    // Generate random inputs once, but use volatile to prevent optimization
    volatile uint64_t seed_a = 42;
    volatile uint64_t seed_b = 123;
    auto a = generate_random_fp12(ring, seed_a);
    auto b = generate_random_fp12(ring, seed_b);
    
    // Make inputs appear as black box to compiler
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);
    
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_mul(a, b, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("FP12_Mul");
}

static void BM_FP12_Sqr(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12(ring, seed);
    benchmark::DoNotOptimize(a);
    
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_sqr(a, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("FP12_Sqr");
}

static void BM_FP12_CyclotomicSqr(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12(ring, seed);
    benchmark::DoNotOptimize(a);
    
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_cyclotomic_sqr(a, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("FP12_CyclotomicSqr");
}

static void BM_FP12_Conjugate(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12(ring, seed);
    benchmark::DoNotOptimize(a);
    
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_conjugate(a, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("FP12_Conjugate");
}

static void BM_FP12_FrobeniusMap1(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12(ring, seed);
    benchmark::DoNotOptimize(a);
    
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_frobenius_map1(a, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("FP12_FrobeniusMap1");
}

static void BM_FP12_FrobeniusMap2(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12(ring, seed);
    benchmark::DoNotOptimize(a);
    
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_frobenius_map2(a, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("FP12_FrobeniusMap2");
}

static void BM_FP12_FrobeniusMap3(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12(ring, seed);
    benchmark::DoNotOptimize(a);
    
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_frobenius_map3(a, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("FP12_FrobeniusMap3");
}

static void BM_FP12_FinalExpAfterInverse(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    
    volatile uint64_t seed_f = 42;
    volatile uint64_t seed_f_inv = 123;
    auto f = generate_random_fp12(ring, seed_f);
    auto f_inverse = generate_random_fp12(ring, seed_f_inv);
    benchmark::DoNotOptimize(f);
    benchmark::DoNotOptimize(f_inverse);
    
    for (auto _ : state) {
        auto result = final_exp_after_inverse(f, f_inverse, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("FP12_FinalExpAfterInverse");
}

static void BM_FP12_MulByXy0z00(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    
    volatile uint64_t seed_x = 42;
    auto x = generate_random_fp12(ring, seed_x);
    
    // Generate random sparse element y with 6 components
    // y is std::array<typename Ring::StandardElement, 6>
    using RingType = decltype(ring);
    std::array<typename RingType::StandardElement, 6> y;
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(123);
    for (int i = 0; i < 6; i++) {
        BigInt val(rng.get_z_range(p.get_mpz()));
        y[i] = ring.from_bigint(val);
    }
    auto c_term = FP2<typename RingType::StandardElement, typename RingType::StandardElement>(y[0], y[1]);
    auto x_term = FP2<typename RingType::StandardElement, typename RingType::StandardElement>(y[2], y[3]);
    auto y_term = FP2<typename RingType::StandardElement, typename RingType::StandardElement>(y[4], y[5]);
    
    benchmark::DoNotOptimize(x);
    benchmark::DoNotOptimize(y);
    
    for (auto _ : state) {
        typename FP12Ring<RingType>::StandardElement result;
        fp12ring.flat_mulxy00z0(result, x, c_term, x_term, y_term, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("FP12_MulByXy0z00");
}

static void BM_Inverse(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    BLS381AddChainInversion<BoundedRing<>> inverter(ring);
    
    // Generate random field element
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    BigInt random_val(rng.get_z_range(p.get_mpz()));
    if (random_val == BigInt(0)) {
        random_val = BigInt(1); // Avoid zero
    }
    
    auto element = ring.from_bigint(random_val);
    benchmark::DoNotOptimize(element);
    
    for (auto _ : state) {
        auto result = inverter.invert(element, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("Inverse");
}

static void BM_MillerLoop_DoubleStep(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP2Ring<BoundedRing<>> fp2ring(ring);
    FP12Ring<BoundedRing<>> fp12ring(ring);
    
    // Set up G1 and G2 points
    BigInt g1_x(g1_gen_x, 10);
    BigInt g1_y(g1_gen_y, 10);
    AffinePoint<typename BoundedRing<>::StandardElement> E1_point{
        ring.from_bigint(g1_x),
        ring.from_bigint(g1_y)
    };
    
    BigInt g2_x_c0(g2_gen_x_c0, 10);
    BigInt g2_x_c1(g2_gen_x_c1, 10);
    BigInt g2_y_c0(g2_gen_y_c0, 10);
    BigInt g2_y_c1(g2_gen_y_c1, 10);
    AffinePoint<FP2Ring<BoundedRing<>>::StandardElement> E2_point{
        FP2Ring<BoundedRing<>>::StandardElement(
            ring.from_bigint(g2_x_c0),
            ring.from_bigint(g2_x_c1)
        ),
        FP2Ring<BoundedRing<>>::StandardElement(
            ring.from_bigint(g2_y_c0),
            ring.from_bigint(g2_y_c1)
        )
    };
    
    // Create Miller object
    Miller<BoundedRing<>> miller(E1_point, E2_point, ring);
    
    // Set up projective point and f for double step
    ProjectivePoint<FP2Ring<BoundedRing<>>::StandardElement> P{
        E2_point.x,
        E2_point.y,
        FP2Ring<BoundedRing<>>::one(ring)
    };
    auto f = fp12ring.one(ring);
    
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(f);
    
    for (auto _ : state) {
        P = miller.double_step(P, f, ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(P);
        benchmark::DoNotOptimize(f);
    }
    
    state.SetLabel("MillerLoop_DoubleStep");
}

static void BM_MillerLoop_AddStep(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP2Ring<BoundedRing<>> fp2ring(ring);
    FP12Ring<BoundedRing<>> fp12ring(ring);
    
    // Set up G1 and G2 points
    BigInt g1_x(g1_gen_x, 10);
    BigInt g1_y(g1_gen_y, 10);
    AffinePoint<typename BoundedRing<>::StandardElement> E1_point{
        ring.from_bigint(g1_x),
        ring.from_bigint(g1_y)
    };
    
    BigInt g2_x_c0(g2_gen_x_c0, 10);
    BigInt g2_x_c1(g2_gen_x_c1, 10);
    BigInt g2_y_c0(g2_gen_y_c0, 10);
    BigInt g2_y_c1(g2_gen_y_c1, 10);
    AffinePoint<FP2Ring<BoundedRing<>>::StandardElement> E2_point{
        FP2Ring<BoundedRing<>>::StandardElement(
            ring.from_bigint(g2_x_c0),
            ring.from_bigint(g2_x_c1)
        ),
        FP2Ring<BoundedRing<>>::StandardElement(
            ring.from_bigint(g2_y_c0),
            ring.from_bigint(g2_y_c1)
        )
    };
    
    // Create Miller object
    Miller<BoundedRing<>> miller(E1_point, E2_point, ring);
    
    // Set up projective point (after one double step) and f for add step
    ProjectivePoint<FP2Ring<BoundedRing<>>::StandardElement> P{
        E2_point.x,
        E2_point.y,
        FP2Ring<BoundedRing<>>::one(ring)
    };
    auto f = fp12ring.one(ring);
    // Do one double step to get a valid state for add step
    P = miller.double_step(P, f, ring, fp2ring, fp12ring);
    
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(f);
    
    for (auto _ : state) {
        P = miller.add_step(P, f, ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(P);
        benchmark::DoNotOptimize(f);
    }
    
    state.SetLabel("MillerLoop_AddStep");
}

static void BM_MillerLoop(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP2Ring<BoundedRing<>> fp2ring(ring);
    FP12Ring<BoundedRing<>> fp12ring(ring);
    
    // Set up G1 and G2 points
    BigInt g1_x(g1_gen_x, 10);
    BigInt g1_y(g1_gen_y, 10);
    AffinePoint<typename BoundedRing<>::StandardElement> E1_point{
        ring.from_bigint(g1_x),
        ring.from_bigint(g1_y)
    };
    
    BigInt g2_x_c0(g2_gen_x_c0, 10);
    BigInt g2_x_c1(g2_gen_x_c1, 10);
    BigInt g2_y_c0(g2_gen_y_c0, 10);
    BigInt g2_y_c1(g2_gen_y_c1, 10);
    AffinePoint<FP2Ring<BoundedRing<>>::StandardElement> E2_point{
        FP2Ring<BoundedRing<>>::StandardElement(
            ring.from_bigint(g2_x_c0),
            ring.from_bigint(g2_x_c1)
        ),
        FP2Ring<BoundedRing<>>::StandardElement(
            ring.from_bigint(g2_y_c0),
            ring.from_bigint(g2_y_c1)
        )
    };
    
    // Create Miller object
    Miller<BoundedRing<>> miller(E1_point, E2_point, ring);
    
    for (auto _ : state) {
        auto result = miller.miller_loop(ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("MillerLoop");
}

// Benchmark FP12 inversion using BLST-backed inverter
static void BM_FP12_Inverse_BLST(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);

    volatile uint64_t seed = 42;
    auto a = generate_random_fp12(ring, seed);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        auto inv = invert_fp12_via_blst(ring, a);
        benchmark::DoNotOptimize(inv);
    }

    state.SetLabel("FP12_Inverse_BLST");
}

// Benchmark FP12 inversion using RNS flat inverse with BLST for single-Fp inverse
static void BM_FP12_Inverse_RNS_BLST(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    FP12RNS_BLSTInverter<BoundedRing<>> inverter;

    volatile uint64_t seed = 42;
    auto a = generate_random_fp12(ring, seed);
    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        auto inv = inverter.invert(a, fp12ring, ring);
        benchmark::DoNotOptimize(inv);
    }

    state.SetLabel("FP12_Inverse_RNS_BLST");
}

// Benchmark final exponentiation using BLST-backed FP12 inverter
static void BM_FP12_FinalExp_BLST_Inverter(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    FP12BLSTInverter<BoundedRing<>> inverter;

    volatile uint64_t seed_f = 42;
    auto f = generate_random_fp12(ring, seed_f);
    benchmark::DoNotOptimize(f);

    for (auto _ : state) {
        auto result = final_exp(f, inverter, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FinalExp_BLST_Inverter");
}

// Benchmark final exponentiation using RNS flat inverse with BLST for single-Fp inverse
static void BM_FP12_FinalExp_RNS_BLST_Inverter(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP12Ring fp12ring(ring);
    FP12RNS_BLSTInverter<BoundedRing<>> inverter;

    volatile uint64_t seed_f = 42;
    auto f = generate_random_fp12(ring, seed_f);
    benchmark::DoNotOptimize(f);

    for (auto _ : state) {
        auto result = final_exp(f, inverter, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP12_FinalExp_RNS_BLST_Inverter");
}

// Benchmark full pairing (Miller loop + final exponentiation) using BLST-backed FP12 inverter
static void BM_Pairing_BLST_Inverter(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP2Ring<BoundedRing<>> fp2ring(ring);
    FP12Ring<BoundedRing<>> fp12ring(ring);

    // Set up G1 and G2 points (generators)
    BigInt g1_x(g1_gen_x, 10);
    BigInt g1_y(g1_gen_y, 10);
    AffinePoint<typename BoundedRing<>::StandardElement> E1_point{
        ring.from_bigint(g1_x),
        ring.from_bigint(g1_y)
    };

    BigInt g2_x_c0(g2_gen_x_c0, 10);
    BigInt g2_x_c1(g2_gen_x_c1, 10);
    BigInt g2_y_c0(g2_gen_y_c0, 10);
    BigInt g2_y_c1(g2_gen_y_c1, 10);
    AffinePoint<FP2Ring<BoundedRing<>>::StandardElement> E2_point{
        FP2Ring<BoundedRing<>>::StandardElement(
            ring.from_bigint(g2_x_c0),
            ring.from_bigint(g2_x_c1)
        ),
        FP2Ring<BoundedRing<>>::StandardElement(
            ring.from_bigint(g2_y_c0),
            ring.from_bigint(g2_y_c1)
        )
    };

    FP12BLSTInverter<BoundedRing<>> inverter;

    for (auto _ : state) {
        auto result = pairing(E1_point, E2_point, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_BLST_Inverter");
}

// Benchmark full pairing using RNS flat inverse with BLST for single-Fp inverse
static void BM_Pairing_RNS_BLST_Inverter(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP2Ring<BoundedRing<>> fp2ring(ring);
    FP12Ring<BoundedRing<>> fp12ring(ring);

    BigInt g1_x(g1_gen_x, 10);
    BigInt g1_y(g1_gen_y, 10);
    AffinePoint<typename BoundedRing<>::StandardElement> E1_point{
        ring.from_bigint(g1_x),
        ring.from_bigint(g1_y)
    };

    BigInt g2_x_c0(g2_gen_x_c0, 10);
    BigInt g2_x_c1(g2_gen_x_c1, 10);
    BigInt g2_y_c0(g2_gen_y_c0, 10);
    BigInt g2_y_c1(g2_gen_y_c1, 10);
    AffinePoint<FP2Ring<BoundedRing<>>::StandardElement> E2_point{
        FP2Ring<BoundedRing<>>::StandardElement(
            ring.from_bigint(g2_x_c0),
            ring.from_bigint(g2_x_c1)
        ),
        FP2Ring<BoundedRing<>>::StandardElement(
            ring.from_bigint(g2_y_c0),
            ring.from_bigint(g2_y_c1)
        )
    };

    FP12RNS_BLSTInverter<BoundedRing<>> inverter;

    for (auto _ : state) {
        auto result = pairing(E1_point, E2_point, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("Pairing_RNS_BLST_Inverter");
}

static void BM_BatchReduceExpand_1(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    
    // Generate 1 random element
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    
    std::array<typename BoundedRing<>::StandardElement, 1> a_elements;
    std::array<typename BoundedRing<>::StandardElement, 1> b_elements;
    for (int i = 0; i < 1; i++) {
        BigInt val_a(rng.get_z_range(p.get_mpz()));
        BigInt val_b(rng.get_z_range(p.get_mpz()));
        a_elements[i] = ring.from_bigint(val_a);
        b_elements[i] = ring.from_bigint(val_b);
    }
    
    // Create products and convert to ReadyElement array
    std::array<typename BoundedRing<>::ReadyElement<BoundedRing<>::MAX_ADD>, 1> ready_elements;
    for (int i = 0; i < 1; i++) {
        auto product = a_elements[i] * b_elements[i];
        ready_elements[i] = ring.template ready<BoundedRing<>::MAX_ADD>(product);
    }
    
    benchmark::DoNotOptimize(ready_elements);
    
    for (auto _ : state) {
        auto result = ring.batch_reduce_expand<1>(ready_elements);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("BatchReduceExpand_1");
}

static void BM_BatchReduceExpand_2(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    
    // Generate 2 random elements
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    
    std::array<typename BoundedRing<>::StandardElement, 2> a_elements;
    std::array<typename BoundedRing<>::StandardElement, 2> b_elements;
    for (int i = 0; i < 2; i++) {
        BigInt val_a(rng.get_z_range(p.get_mpz()));
        BigInt val_b(rng.get_z_range(p.get_mpz()));
        a_elements[i] = ring.from_bigint(val_a);
        b_elements[i] = ring.from_bigint(val_b);
    }
    
    // Create products and convert to ReadyElement array
    std::array<typename BoundedRing<>::ReadyElement<BoundedRing<>::MAX_ADD>, 2> ready_elements;
    for (int i = 0; i < 2; i++) {
        auto product = a_elements[i] * b_elements[i];
        ready_elements[i] = ring.template ready<BoundedRing<>::MAX_ADD>(product);
    }
    
    benchmark::DoNotOptimize(ready_elements);
    
    for (auto _ : state) {
        auto result = ring.batch_reduce_expand<2>(ready_elements);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("BatchReduceExpand_2");
}

static void BM_BatchReduceExpand_3(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    
    // Generate 3 random elements
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    
    std::array<typename BoundedRing<>::StandardElement, 3> a_elements;
    std::array<typename BoundedRing<>::StandardElement, 3> b_elements;
    for (int i = 0; i < 3; i++) {
        BigInt val_a(rng.get_z_range(p.get_mpz()));
        BigInt val_b(rng.get_z_range(p.get_mpz()));
        a_elements[i] = ring.from_bigint(val_a);
        b_elements[i] = ring.from_bigint(val_b);
    }
    
    // Create products and convert to ReadyElement array
    std::array<typename BoundedRing<>::ReadyElement<BoundedRing<>::MAX_ADD>, 3> ready_elements;
    for (int i = 0; i < 3; i++) {
        auto product = a_elements[i] * b_elements[i];
        ready_elements[i] = ring.template ready<BoundedRing<>::MAX_ADD>(product);
    }
    
    benchmark::DoNotOptimize(ready_elements);
    
    for (auto _ : state) {
        auto result = ring.batch_reduce_expand<3>(ready_elements);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("BatchReduceExpand_3");
}

static void BM_BatchReduceExpand_4(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    
    // Generate 4 random elements
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    
    std::array<typename BoundedRing<>::StandardElement, 4> a_elements;
    std::array<typename BoundedRing<>::StandardElement, 4> b_elements;
    for (int i = 0; i < 4; i++) {
        BigInt val_a(rng.get_z_range(p.get_mpz()));
        BigInt val_b(rng.get_z_range(p.get_mpz()));
        a_elements[i] = ring.from_bigint(val_a);
        b_elements[i] = ring.from_bigint(val_b);
    }
    
    // Create products and convert to ReadyElement array
    std::array<typename BoundedRing<>::ReadyElement<BoundedRing<>::MAX_ADD>, 4> ready_elements;
    for (int i = 0; i < 4; i++) {
        auto product = a_elements[i] * b_elements[i];
        ready_elements[i] = ring.template ready<BoundedRing<>::MAX_ADD>(product);
    }
    
    benchmark::DoNotOptimize(ready_elements);
    
    for (auto _ : state) {
        auto result = ring.batch_reduce_expand<4>(ready_elements);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("BatchReduceExpand_4");
}

static void BM_BatchReduceExpand_5(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    
    // Generate 5 random elements
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    
    std::array<typename BoundedRing<>::StandardElement, 5> a_elements;
    std::array<typename BoundedRing<>::StandardElement, 5> b_elements;
    for (int i = 0; i < 5; i++) {
        BigInt val_a(rng.get_z_range(p.get_mpz()));
        BigInt val_b(rng.get_z_range(p.get_mpz()));
        a_elements[i] = ring.from_bigint(val_a);
        b_elements[i] = ring.from_bigint(val_b);
    }
    
    // Create products and convert to ReadyElement array
    std::array<typename BoundedRing<>::ReadyElement<BoundedRing<>::MAX_ADD>, 5> ready_elements;
    for (int i = 0; i < 5; i++) {
        auto product = a_elements[i] * b_elements[i];
        ready_elements[i] = ring.template ready<BoundedRing<>::MAX_ADD>(product);
    }
    
    benchmark::DoNotOptimize(ready_elements);
    
    for (auto _ : state) {
        auto result = ring.batch_reduce_expand<5>(ready_elements);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("BatchReduceExpand_5");
}

static void BM_BatchReduceExpand_6(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    
    // Generate 6 random elements
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    
    std::array<typename BoundedRing<>::StandardElement, 6> a_elements;
    std::array<typename BoundedRing<>::StandardElement, 6> b_elements;
    for (int i = 0; i < 6; i++) {
        BigInt val_a(rng.get_z_range(p.get_mpz()));
        BigInt val_b(rng.get_z_range(p.get_mpz()));
        a_elements[i] = ring.from_bigint(val_a);
        b_elements[i] = ring.from_bigint(val_b);
    }
    
    // Create products and convert to ReadyElement array
    std::array<typename BoundedRing<>::ReadyElement<BoundedRing<>::MAX_ADD>, 6> ready_elements;
    for (int i = 0; i < 6; i++) {
        auto product = a_elements[i] * b_elements[i];
        ready_elements[i] = ring.template ready<BoundedRing<>::MAX_ADD>(product);
    }
    
    benchmark::DoNotOptimize(ready_elements);
    
    for (auto _ : state) {
        auto result = ring.batch_reduce_expand<6>(ready_elements);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("BatchReduceExpand_6");
}

static void BM_BatchReduceExpand_12(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    
    // Generate 12 random elements
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    
    std::array<typename BoundedRing<>::StandardElement, 12> a_elements;
    std::array<typename BoundedRing<>::StandardElement, 12> b_elements;
    for (int i = 0; i < 12; i++) {
        BigInt val_a(rng.get_z_range(p.get_mpz()));
        BigInt val_b(rng.get_z_range(p.get_mpz()));
        a_elements[i] = ring.from_bigint(val_a);
        b_elements[i] = ring.from_bigint(val_b);
    }
    
    // Create products and convert to ReadyElement array
    std::array<typename BoundedRing<>::ReadyElement<BoundedRing<>::MAX_ADD>, 12> ready_elements;
    for (int i = 0; i < 12; i++) {
        auto product = a_elements[i] * b_elements[i];
        ready_elements[i] = ring.template ready<BoundedRing<>::MAX_ADD>(product);
    }
    
    benchmark::DoNotOptimize(ready_elements);
    
    for (auto _ : state) {
        auto result = ring.batch_reduce_expand<12>(ready_elements);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("BatchReduceExpand_12");
}

// Benchmark FP2 multiplication using FP2Ring
static void BM_FP2_Mul(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP2Ring<BoundedRing<>> fp2ring(ring);

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);

    BigInt ax(rng.get_z_range(p.get_mpz()));
    BigInt ay(rng.get_z_range(p.get_mpz()));
    BigInt bx(rng.get_z_range(p.get_mpz()));
    BigInt by(rng.get_z_range(p.get_mpz()));

    auto a = FP2Ring<BoundedRing<>>::StandardElement(
        ring.from_bigint(ax),
        ring.from_bigint(ay)
    );
    auto b = FP2Ring<BoundedRing<>>::StandardElement(
        ring.from_bigint(bx),
        ring.from_bigint(by)
    );

    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        auto result = fp2ring.mul(a, b);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP2_Mul");
}

// Benchmark FP2 squaring using FP2Ring
static void BM_FP2_Sqr(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    FP2Ring<BoundedRing<>> fp2ring(ring);

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);

    BigInt ax(rng.get_z_range(p.get_mpz()));
    BigInt ay(rng.get_z_range(p.get_mpz()));

    auto a = FP2Ring<BoundedRing<>>::StandardElement(
        ring.from_bigint(ax),
        ring.from_bigint(ay)
    );

    benchmark::DoNotOptimize(a);

    for (auto _ : state) {
        auto result = fp2ring.square(a);
        benchmark::DoNotOptimize(result);
    }

    state.SetLabel("FP2_Sqr");
}
static void BM_ModMul(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing ring(p);
    
    // Generate random field elements
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    
    BigInt val_a(rng.get_z_range(p.get_mpz()));
    BigInt val_b(rng.get_z_range(p.get_mpz()));
    
    auto a = ring.from_bigint(val_a);
    auto b = ring.from_bigint(val_b);
    
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);
    
    for (auto _ : state) {
        auto result = ring.modmul(a, b);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("ModMul");
}

// Helper to generate random G1 point (just random field elements)
template<class Ring>
typename G1<Ring>::ProjPoint generate_random_g1_projective(const Ring &ring, uint64_t seed) {
    BigInt p = BigInt(modulus, 10);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    BigInt x(rng.get_z_range(p.get_mpz()));
    BigInt y(rng.get_z_range(p.get_mpz()));
    return typename G1<Ring>::ProjPoint(ring.from_bigint(x), ring.from_bigint(y), ring.one());
}

// Helper to generate random G1 affine point
template<class Ring>
AffinePoint<typename Ring::StandardElement> generate_random_g1_affine(const Ring &ring, uint64_t seed) {
    BigInt p = BigInt(modulus, 10);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    BigInt x(rng.get_z_range(p.get_mpz()));
    BigInt y(rng.get_z_range(p.get_mpz()));
    return AffinePoint<typename Ring::StandardElement>{ring.from_bigint(x), ring.from_bigint(y)};
}

// Helper to generate random G2 point (just random FP2 field elements)
template<class Ring>
typename G2<Ring>::ProjPoint generate_random_g2_projective(const FP2Ring<Ring> &fp2ring, const Ring &ring, uint64_t seed) {
    BigInt p = BigInt(modulus, 10);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    BigInt x0(rng.get_z_range(p.get_mpz()));
    BigInt x1(rng.get_z_range(p.get_mpz()));
    BigInt y0(rng.get_z_range(p.get_mpz()));
    BigInt y1(rng.get_z_range(p.get_mpz()));
    return typename G2<Ring>::ProjPoint(
        typename FP2Ring<Ring>::StandardElement(ring.from_bigint(x0), ring.from_bigint(x1)),
        typename FP2Ring<Ring>::StandardElement(ring.from_bigint(y0), ring.from_bigint(y1)),
        fp2ring.one()
    );
}

// Helper to generate random G2 affine point
template<class Ring>
AffinePoint<typename FP2Ring<Ring>::StandardElement> generate_random_g2_affine(const Ring &ring, uint64_t seed) {
    BigInt p = BigInt(modulus, 10);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    BigInt x0(rng.get_z_range(p.get_mpz()));
    BigInt x1(rng.get_z_range(p.get_mpz()));
    BigInt y0(rng.get_z_range(p.get_mpz()));
    BigInt y1(rng.get_z_range(p.get_mpz()));
    FP2Ring<Ring> fp2ring(ring);
    return AffinePoint<typename FP2Ring<Ring>::StandardElement>{
        typename FP2Ring<Ring>::StandardElement(ring.from_bigint(x0), ring.from_bigint(x1)),
        typename FP2Ring<Ring>::StandardElement(ring.from_bigint(y0), ring.from_bigint(y1))
    };
}

static void BM_G1_Add(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing<381, 8, 52, -1932, 2377, 12> ring(p);
    G1<decltype(ring)> g1;
    
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    auto P = generate_random_g1_projective(ring, seed1);
    auto Q = generate_random_g1_projective(ring, seed2);
    
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(Q);
    
    for (auto _ : state) {
        auto result = g1.add_point(P, Q, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("G1_Add");
}

static void BM_G1_MixedAdd(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing<381, 8, 52, -1932, 2377, 12> ring(p);
    G1<decltype(ring)> g1;
    
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    auto P = generate_random_g1_projective(ring, seed1);
    auto Q_affine = generate_random_g1_affine(ring, seed2);
    
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(Q_affine);
    
    for (auto _ : state) {
        auto result = g1.add_mixed_point(P, Q_affine, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("G1_MixedAdd");
}

static void BM_G1_Double(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing<381, 8, 52, -1932, 2377, 12> ring(p);
    G1<decltype(ring)> g1;
    
    volatile uint64_t seed = 42;
    auto P = generate_random_g1_projective(ring, seed);
    
    benchmark::DoNotOptimize(P);
    
    for (auto _ : state) {
        auto result = g1.double_point(P, ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("G1_Double");
}

static void BM_G2_Add(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing<381, 8, 52, -1932, 2377, 12> ring(p);
    FP2Ring<decltype(ring)> fp2ring(ring);
    G2<decltype(ring)> g2;
    
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    auto P = generate_random_g2_projective(fp2ring, ring, seed1);
    auto Q = generate_random_g2_projective(fp2ring, ring, seed2);
    
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(Q);
    
    for (auto _ : state) {
        auto result = g2.add_point(P, Q, fp2ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("G2_Add");
}

static void BM_G2_MixedAdd(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing<381, 8, 52, -1932, 2377, 12> ring(p);
    FP2Ring<decltype(ring)> fp2ring(ring);
    G2<decltype(ring)> g2;
    
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    auto P = generate_random_g2_projective(fp2ring, ring, seed1);
    auto Q_affine = generate_random_g2_affine(ring, seed2);
    
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(Q_affine);
    
    for (auto _ : state) {
        auto result = g2.add_mixed_point(P, Q_affine, fp2ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("G2_MixedAdd");
}

static void BM_G2_Double(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing<381, 8, 52, -1932, 2377, 12> ring(p);
    FP2Ring<decltype(ring)> fp2ring(ring);
    G2<decltype(ring)> g2;
    
    volatile uint64_t seed = 42;
    auto P = generate_random_g2_projective(fp2ring, ring, seed);
    
    benchmark::DoNotOptimize(P);
    
    for (auto _ : state) {
        auto result = g2.double_point(P, fp2ring);
        benchmark::DoNotOptimize(result);
    }
    
    state.SetLabel("G2_Double");
}

// Benchmark GLV scalar multiplication on G1 (byte array version, 52-bit configuration)
static void BM_G1_ScalarMult_GLV(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing<381, 8, 52, -1932, 2377, 12> ring(p);
    using RingType = decltype(ring);
    G1<RingType> g1_curve;
    GLVMult<RingType> glv_mult(g1_curve, ring);

    // Generate random projective point
    volatile uint64_t seed_point = 42;
    auto P = generate_random_g1_projective(ring, seed_point);

    // Generate random scalar in [0, r)
    BigInt r_scalar(bls12_381_scalar_modulus_hex, 16);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(123);
    BigInt scalar_bigint(rng.get_z_range(r_scalar.get_mpz()));

    // Convert scalar to 32-byte little-endian array
    std::array<uint8_t, 32> scalar_bytes = {0};
    bigint_to_bytes_le(scalar_bytes.data(), scalar_bigint, 32);

    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(scalar_bytes);

    for (auto _ : state) {
        auto R = glv_mult.multiply(P, scalar_bytes, ring);
        benchmark::DoNotOptimize(R);
    }

    state.SetLabel("G1_ScalarMult_GLV");
}

// Benchmark GLS scalar multiplication on G2 (byte array version, 52-bit configuration)
static void BM_G2_ScalarMult_GLS(benchmark::State& state) {
    BigInt p = BigInt(modulus, 10);
    BoundedRing<381, 8, 52, -1932, 2377, 12> ring(p);
    using RingType = decltype(ring);
    using FP2RingType = FP2Ring<RingType>;

    FP2RingType fp2ring(ring);
    GLSMult<RingType> gls_mult(ring);
    typename GLSMult<RingType>::G2Curve g2_curve;

    // Generate random projective point
    volatile uint64_t seed_point = 42;
    auto P = generate_random_g2_projective(fp2ring, ring, seed_point);

    // Generate random scalar in [0, r)
    BigInt r_scalar(bls12_381_scalar_modulus_hex, 16);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(321);
    BigInt scalar_bigint(rng.get_z_range(r_scalar.get_mpz()));

    // Convert scalar to 32-byte little-endian array
    std::array<uint8_t, 32> scalar_bytes = {0};
    bigint_to_bytes_le(scalar_bytes.data(), scalar_bigint, 32);

    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(scalar_bytes);

    for (auto _ : state) {
        auto R = gls_mult.multiply(P, scalar_bytes, fp2ring, g2_curve, ring);
        benchmark::DoNotOptimize(R);
    }

    state.SetLabel("G2_ScalarMult_GLS");
}

BENCHMARK(BM_ModMul);
BENCHMARK(BM_BatchReduceExpand_1);
BENCHMARK(BM_BatchReduceExpand_2);
BENCHMARK(BM_BatchReduceExpand_3);
BENCHMARK(BM_BatchReduceExpand_4);
BENCHMARK(BM_BatchReduceExpand_5);
BENCHMARK(BM_BatchReduceExpand_6);
BENCHMARK(BM_BatchReduceExpand_12);
BENCHMARK(BM_FP2_Mul);
BENCHMARK(BM_FP2_Sqr);
BENCHMARK(BM_G1_Add);
BENCHMARK(BM_G1_MixedAdd);
BENCHMARK(BM_G1_Double);
BENCHMARK(BM_G2_Add);
BENCHMARK(BM_G2_MixedAdd);
BENCHMARK(BM_G2_Double);
BENCHMARK(BM_G1_ScalarMult_GLV);
BENCHMARK(BM_G2_ScalarMult_GLS);
BENCHMARK(BM_FP12_Mul);
BENCHMARK(BM_FP12_Sqr);
BENCHMARK(BM_FP12_CyclotomicSqr);
BENCHMARK(BM_FP12_Conjugate);
BENCHMARK(BM_FP12_FrobeniusMap1);
BENCHMARK(BM_FP12_FrobeniusMap2);
BENCHMARK(BM_FP12_FrobeniusMap3);
BENCHMARK(BM_FP12_MulByXy0z00);
BENCHMARK(BM_FP12_FinalExpAfterInverse);
BENCHMARK(BM_Inverse);
BENCHMARK(BM_MillerLoop_DoubleStep);
BENCHMARK(BM_MillerLoop_AddStep);
BENCHMARK(BM_MillerLoop);
BENCHMARK(BM_FP12_Inverse_BLST);
BENCHMARK(BM_FP12_Inverse_RNS_BLST);
BENCHMARK(BM_FP12_FinalExp_BLST_Inverter);
BENCHMARK(BM_FP12_FinalExp_RNS_BLST_Inverter);
BENCHMARK(BM_Pairing_BLST_Inverter);
BENCHMARK(BM_Pairing_RNS_BLST_Inverter);

BENCHMARK_MAIN();

