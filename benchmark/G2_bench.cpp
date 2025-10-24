#include <benchmark/benchmark.h>
#include <random>
#include <cstring>
#include <stdint.h>

// Include the actual BLST curve and field headers
#include "blst/curve/point.h"

#include "blst/curve/map_to_G.h"

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
static void generate_random_g2_point(G2Point& p) {
    byte scalar[32];
    generate_random_scalar(scalar);
    G2Point::mult(p, G2Point::generator(), scalar, 256);
}

// ========== POINT ADDITION BENCHMARKS ==========

// Benchmark for G2 point addition (Jacobian)
static void BM_G2_Add(benchmark::State& state) {
    G2Point a, b, c;
    generate_random_g2_point(a);
    generate_random_g2_point(b);

    for (auto _ : state) {
        G2Point::add(c, a, b);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Add);

// Benchmark for G2 point addition with an affine point
static void BM_G2_AddAffine(benchmark::State& state) {
    G2Point a, c;
    G2AffinePoint b_affine;
    generate_random_g2_point(a);

    G2Point b;
    generate_random_g2_point(b);
    G2Point::to_affine(b_affine, b);

    for (auto _ : state) {
        G2Point::add(c, a, b_affine);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_AddAffine);

// ========== POINT DOUBLING BENCHMARK ==========

// Benchmark for G2 point doubling
static void BM_G2_Double(benchmark::State& state) {
    G2Point a, c;
    generate_random_g2_point(a);

    for (auto _ : state) {
        G2Point::doubling(c, a);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Double);

// ========== SCALAR MULTIPLICATION BENCHMARK ==========

// Benchmark for G2 scalar multiplication
static void BM_G2_ScalarMult(benchmark::State& state) {
    G2Point p, result;
    byte scalar[32];

    p = G2Point::generator();
    generate_random_scalar(scalar);

    for (auto _ : state) {
        G2Point::mult(result, p, scalar, 256);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_ScalarMult);

// ========== COORDINATE CONVERSION BENCHMARKS ==========

// Benchmark for converting from Jacobian to Affine coordinates
static void BM_G2_ToAffine(benchmark::State& state) {
    G2Point p;
    G2AffinePoint p_affine;
    generate_random_g2_point(p);

    for (auto _ : state) {
        G2Point::to_affine(p_affine, p);
        benchmark::DoNotOptimize(p_affine);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_ToAffine);

// Benchmark for converting from Affine to Jacobian coordinates
static void BM_G2_FromAffine(benchmark::State& state) {
    G2Point p;
    G2AffinePoint p_affine = G2AffinePoint::generator();

    for (auto _ : state) {
        G2Point::from_affine(p, p_affine);
        benchmark::DoNotOptimize(p);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_FromAffine);

// ========== (DE)SERIALIZATION & (DE)COMPRESSION BENCHMARKS ==========

// Benchmark for G2 point serialization
static void BM_G2_Serialize(benchmark::State& state) {
    G2Point p;
    generate_random_g2_point(p);
    unsigned char out[192];

    for (auto _ : state) {
        G2Point::serialize(out, p);
        benchmark::DoNotOptimize(out);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Serialize);

// Benchmark for G2 point deserialization
static void BM_G2_Deserialize(benchmark::State& state) {
    G2Point p = G2Point::generator();
    unsigned char serialized[192];
    G2Point::serialize(serialized, p);
    G2AffinePoint deserialized;

    for (auto _ : state) {
        G2AffinePoint::deserialize(deserialized, serialized);
        benchmark::DoNotOptimize(deserialized);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Deserialize);

// Benchmark for G2 point compression
static void BM_G2_Compress(benchmark::State& state) {
    G2Point p;
    generate_random_g2_point(p);
    unsigned char out[96];

    for (auto _ : state) {
        G2Point::compress(out, p);
        benchmark::DoNotOptimize(out);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Compress);

// Benchmark for G2 point uncompression
static void BM_G2_Uncompress(benchmark::State& state) {
    G2Point p = G2Point::generator();
    unsigned char compressed[96];
    G2Point::compress(compressed, p);
    G2AffinePoint uncompressed;

    for (auto _ : state) {
        G2AffinePoint::uncompress(uncompressed, compressed);
        benchmark::DoNotOptimize(uncompressed);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_Uncompress);

// ========== HASH TO CURVE BENCHMARK ==========

// Benchmark for hashing to G2
static void BM_G2_HashToCurve(benchmark::State& state) {
    G2Point p;
    unsigned char msg[] = "sample message for hashing to curve";
    unsigned char dst[] = "BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_";

    for (auto _ : state) {
        Hash_to_curve<G2Point>(p, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);
        benchmark::DoNotOptimize(p);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_HashToCurve);

// ========== KEY GENERATION BENCHMARK ==========

// Benchmark for secret key to public key generation
static void BM_G2_SkToPk(benchmark::State& state) {
    G2Point pk;
    byte sk[32];
    generate_random_scalar(sk);

    for (auto _ : state) {
        G2Point::to_pk(pk, sk);
        benchmark::DoNotOptimize(pk);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G2_SkToPk);