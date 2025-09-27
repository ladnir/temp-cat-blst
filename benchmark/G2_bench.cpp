#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>

// Include the actual BLST curve and field headers
#include "blst/curve/e2.h"
#include "blst/curve/map_to_g2.h"

using namespace blst;

// Helper to generate a random 256-bit scalar
static void generate_random_scalar(byte* scalar) {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint64_t> dis;

    uint64_t* as_limbs = reinterpret_cast<uint64_t*>(scalar);
    for (size_t i = 0; i < 256 / 64; ++i) {
        as_limbs[i] = dis(gen);
    }
}

// Helper to generate a random G2 point by multiplying the generator
static void generate_random_g2_point(POINTonE2& p) {
    byte scalar[32];
    generate_random_scalar(scalar);
    blst_p2_mult(p, blst_p2_generator(), scalar, 256);
}

// ========== POINT ADDITION BENCHMARKS ==========

// Benchmark for G2 point addition (Jacobian)
static void BM_G2_Add(benchmark::State& state) {
    POINTonE2 a, b, c;
    generate_random_g2_point(a);
    generate_random_g2_point(b);

    for (auto _ : state) {
        blst_p2_add(c, a, b);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Add);

// Benchmark for G2 point addition with an affine point
static void BM_G2_AddAffine(benchmark::State& state) {
    POINTonE2 a, c;
    POINTonE2_affine b_affine;
    generate_random_g2_point(a);

    POINTonE2 b;
    generate_random_g2_point(b);
    blst_p2_to_affine(b_affine, b);

    for (auto _ : state) {
        blst_p2_add_affine(c, a, b_affine);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_AddAffine);

// ========== POINT DOUBLING BENCHMARK ==========

// Benchmark for G2 point doubling
static void BM_G2_Double(benchmark::State& state) {
    POINTonE2 a, c;
    generate_random_g2_point(a);

    for (auto _ : state) {
        blst_p2_double(c, a);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Double);

// ========== SCALAR MULTIPLICATION BENCHMARK ==========

// Benchmark for G2 scalar multiplication
static void BM_G2_ScalarMult(benchmark::State& state) {
    POINTonE2 p, result;
    byte scalar[32];

    p = blst_p2_generator();
    generate_random_scalar(scalar);

    for (auto _ : state) {
        blst_p2_mult(result, p, scalar, 256);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_ScalarMult);

// ========== COORDINATE CONVERSION BENCHMARKS ==========

// Benchmark for converting from Jacobian to Affine coordinates
static void BM_G2_ToAffine(benchmark::State& state) {
    POINTonE2 p;
    POINTonE2_affine p_affine;
    generate_random_g2_point(p);

    for (auto _ : state) {
        blst_p2_to_affine(p_affine, p);
        benchmark::DoNotOptimize(p_affine);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_ToAffine);

// Benchmark for converting from Affine to Jacobian coordinates
static void BM_G2_FromAffine(benchmark::State& state) {
    POINTonE2 p;
    POINTonE2_affine p_affine = blst_p2_affine_generator();

    for (auto _ : state) {
        blst_p2_from_affine(p, p_affine);
        benchmark::DoNotOptimize(p);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_FromAffine);

// ========== (DE)SERIALIZATION & (DE)COMPRESSION BENCHMARKS ==========

// Benchmark for G2 point serialization
static void BM_G2_Serialize(benchmark::State& state) {
    POINTonE2 p;
    generate_random_g2_point(p);
    unsigned char out[192];

    for (auto _ : state) {
        blst_p2_serialize(out, p);
        benchmark::DoNotOptimize(out);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Serialize);

// Benchmark for G2 point deserialization
static void BM_G2_Deserialize(benchmark::State& state) {
    POINTonE2 p = blst_p2_generator();
    unsigned char serialized[192];
    blst_p2_serialize(serialized, p);
    POINTonE2_affine deserialized;

    for (auto _ : state) {
        blst_p2_deserialize(deserialized, serialized);
        benchmark::DoNotOptimize(deserialized);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Deserialize);

// Benchmark for G2 point compression
static void BM_G2_Compress(benchmark::State& state) {
    POINTonE2 p;
    generate_random_g2_point(p);
    unsigned char out[96];

    for (auto _ : state) {
        blst_p2_compress(out, p);
        benchmark::DoNotOptimize(out);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Compress);

// Benchmark for G2 point uncompression
static void BM_G2_Uncompress(benchmark::State& state) {
    POINTonE2 p = blst_p2_generator();
    unsigned char compressed[96];
    blst_p2_compress(compressed, p);
    POINTonE2_affine uncompressed;

    for (auto _ : state) {
        blst_p2_uncompress(uncompressed, compressed);
        benchmark::DoNotOptimize(uncompressed);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Uncompress);

// ========== HASH TO CURVE BENCHMARK ==========

// Benchmark for hashing to G2
static void BM_G2_HashToCurve(benchmark::State& state) {
    POINTonE2 p;
    unsigned char msg[] = "sample message for hashing to curve";
    unsigned char dst[] = "BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_";

    for (auto _ : state) {
        blst_hash_to_g2(p, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);
        benchmark::DoNotOptimize(p);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_HashToCurve);

// ========== KEY GENERATION BENCHMARK ==========

// Benchmark for secret key to public key generation
static void BM_G2_SkToPk(benchmark::State& state) {
    POINTonE2 pk;
    byte sk[32];
    generate_random_scalar(sk);

    for (auto _ : state) {
        blst_sk_to_pk_in_g2(pk, sk);
        benchmark::DoNotOptimize(pk);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_SkToPk);