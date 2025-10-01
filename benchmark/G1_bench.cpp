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

// Helper to generate a random G1 point by multiplying the generator
static void generate_random_g1_point(G1Point& p) {
    byte scalar[32];
    generate_random_scalar(scalar);
    G1Point::mult(p, G1Point::generator(), scalar, 256);
}

// ========== POINT ADDITION BENCHMARKS ==========

// Benchmark for G1 point addition (Jacobian)
static void BM_G1_Add(benchmark::State& state) {
    G1Point a, b, c;
    generate_random_g1_point(a);
    generate_random_g1_point(b);

    for (auto _ : state) {
        G1Point::add(c, a, b);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_Add);

// Benchmark for G1 point addition with an affine point
static void BM_G1_AddAffine(benchmark::State& state) {
    G1Point a, c;
    G1AffinePoint b_affine;
    generate_random_g1_point(a);

    G1Point b;
    generate_random_g1_point(b);
    G1Point::to_affine(b_affine, b);

    for (auto _ : state) {
        G1Point::add_affine(c, a, b_affine);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_AddAffine);

// ========== POINT DOUBLING BENCHMARK ==========

// Benchmark for G1 point doubling
static void BM_G1_Double(benchmark::State& state) {
    G1Point a, c;
    generate_random_g1_point(a);

    for (auto _ : state) {
        G1Point::doubling(c, a);
        benchmark::DoNotOptimize(c);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_Double);

// ========== SCALAR MULTIPLICATION BENCHMARK ==========

// Benchmark for G1 scalar multiplication
static void BM_G1_ScalarMult(benchmark::State& state) {
    G1Point p, result;
    byte scalar[32];

    p = G1Point::generator();
    generate_random_scalar(scalar);

    for (auto _ : state) {
        G1Point::mult(result, p, scalar, 256);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_ScalarMult);

// ========== COORDINATE CONVERSION BENCHMARKS ==========

// Benchmark for converting from Jacobian to Affine coordinates
static void BM_G1_ToAffine(benchmark::State& state) {
    G1Point p;
    G1AffinePoint p_affine;
    generate_random_g1_point(p);

    for (auto _ : state) {
        G1Point::to_affine(p_affine, p);
        benchmark::DoNotOptimize(p_affine);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_ToAffine);

// Benchmark for converting from Affine to Jacobian coordinates
static void BM_G1_FromAffine(benchmark::State& state) {
    G1Point p;
    G1AffinePoint p_affine = G1AffinePoint::generator();

    for (auto _ : state) {
        blst_p1_from_affine(p, p_affine);
        benchmark::DoNotOptimize(p);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_FromAffine);

// ========== (DE)SERIALIZATION & (DE)COMPRESSION BENCHMARKS ==========

// Benchmark for G1 point serialization
static void BM_G1_Serialize(benchmark::State& state) {
    G1Point p;
    generate_random_g1_point(p);
    unsigned char out[96];

    for (auto _ : state) {
        G1Point::serialize(out, p);
        benchmark::DoNotOptimize(out);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_Serialize);

// Benchmark for G1 point deserialization
static void BM_G1_Deserialize(benchmark::State& state) {
    G1Point p = G1Point::generator();
    unsigned char serialized[96];
    G1Point::serialize(serialized, p);
    G1AffinePoint deserialized;

    for (auto _ : state) {
        G1AffinePoint::deserialize(deserialized, serialized);
        benchmark::DoNotOptimize(deserialized);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_Deserialize);

// Benchmark for G1 point compression
static void BM_G1_Compress(benchmark::State& state) {
    G1Point p;
    generate_random_g1_point(p);
    unsigned char out[48];

    for (auto _ : state) {
        G1Point::compress(out, p);
        benchmark::DoNotOptimize(out);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_Compress);

// Benchmark for G1 point uncompression
static void BM_G1_Uncompress(benchmark::State& state) {
    G1Point p = G1Point::generator();
    unsigned char compressed[48];
    G1Point::compress(compressed, p);
    G1AffinePoint uncompressed;

    for (auto _ : state) {
        G1AffinePoint::uncompress(uncompressed, compressed);
        benchmark::DoNotOptimize(uncompressed);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_Uncompress);

// ========== HASH TO CURVE BENCHMARK ==========

// Benchmark for hashing to G1
static void BM_G1_HashToCurve(benchmark::State& state) {
    G1Point p;
    unsigned char msg[] = "sample message for hashing to curve";
    unsigned char dst[] = "BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_NUL_";

    for (auto _ : state) {
        Hash_to_curve<G1Point>(p, msg, sizeof(msg), dst, sizeof(dst), NULL, 0);
        benchmark::DoNotOptimize(p);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_HashToCurve);

// ========== KEY GENERATION BENCHMARK ==========

// Benchmark for secret key to public key generation
static void BM_G1_SkToPk(benchmark::State& state) {
    G1Point pk;
    byte sk[32];
    generate_random_scalar(sk);

    for (auto _ : state) {
        G1Point::to_pk(pk, sk);
        benchmark::DoNotOptimize(pk);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_G1_SkToPk);

