#include <iostream>
#include <vector>
#include <gtest/gtest.h>
#include "blst/fields/Fr.h"
#include "blst/curve/point.h"
#include "blst/curve/GT.h"
#include "blst/fft/primitveRoot.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/middle_product.h"
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/fft/gt_context.h"
#endif
using namespace blst;
/**
 * Naive polynomial multiplication (discrete convolution)
 * @tparam G The Group type for coefficients of A (e.g., G1)
 * @tparam G2 The Group type for coefficients of B (e.g., G2)
 * @tparam GT The Target Group type for result (e.g., GT)
 */
template <typename G, typename G2, typename GT, typename CoeffCtx = DefaultCoeffCtx>
std::vector<GT> naive_multiply(const std::vector<G>& A, const std::vector<G2>& B) {
    size_t m = A.size();
    size_t n = B.size();

    // The result polynomial has size m + n - 1
    // Initialize with group identity [0]
    std::vector<GT> expected(m + n - 1);

    // O(m * n) nested loops to compute full convolution
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            // Multiplication and accumulation: 
            // result[i+j] = sum(A[i] * B[j])
            auto p = CoeffCtx::mult(A[i], B[j]);
            expected[i + j] = CoeffCtx::plus(expected[i + j], p);
            //expected[i + j] = expected[i + j] + (A[i] * B[j]);
        }
    }

    return expected;
}

//
//TEST(Fft, basic) {
//
//    auto n = 4;
//    size_t N = n*2; // FFT size (2*B) to avoid aliasing [cite: 150]
//
//    std::vector<F65537> A = { 1, 2, 3, 4 }; // Padded keys [cite: 149]
//    std::vector<F65537> B = { 5, 6, 7, 8 }; // Padded CRS [cite: 149]
//
//    // 1. Naive Convolution (The "B^2" bottleneck) [cite: 36, 122]
//    std::vector<F65537> expected = naive_multiply<F65537, F65537, F65537>(A, B);
//	expected.resize(N); // Pad to N
//
//    auto prim = primRootOfUnity<F65537>(N);
//    std::vector<F65537> roots(N);
//    precomputeRootsOfUnity< F65537>(prim, roots);
//    
//    std::vector<F65537> inv_roots(N);
//    for (size_t i = 0; i < N; ++i) {
//        inv_roots[i] = roots[i].inverse();
//    }
//    // 2. FFT Convolution 
//    // (Assume roots are precomputed for simplicity in this test)
//
//    A.resize(N);
//    B.resize(N);
//
//    fft(A, roots, false);
//    fft(B, roots, false);
//
//    std::vector<F65537> C(N);
//    for (int i = 0; i < N; ++i) C[i] = A[i] * B[i]; // Pointwise pairing 
//
//    fft(C, inv_roots, true); // Inverse FFT 
//
//    //std::cout << "FFT Convolution Result (first 4): ";
//    //for (int i = 0; i < 4; ++i) std::cout << C[i] << " ";
//    //std::cout << "\nNaive Result (first 4): ";
//    //for (int i = 0; i < 4; ++i) std::cout << expected[i] << " ";
//
//    EXPECT_EQ(C, expected);
//}


TEST(Fft, BlstFrRoots) {
    auto one = blst::Fr::one();
    blst::Fr negOne;
    blst::Fr::neg(negOne, one);

    EXPECT_TRUE(blst::Fr::equal(blst::primRootOfUnity<blst::Fr>(1), blst::Fr::one()));

    // Verify roots for powers of 2 up to 2^32
    for (int k = 0; k <= 32; ++k) {
        u64 n = 1ULL << k;
        auto w = primRootOfUnity<blst::Fr>(n);

        // Ensure the value is a valid field element (reduced and non-zero)
        EXPECT_TRUE(blst::Fr::check_mod(w)) << "Root is not valid/reduced for k=" << k;

        if (k == 0) {
            EXPECT_TRUE(blst::Fr::equal(w, one)) 
                << " did not get 1 at k=0";
            continue;
        }

        // To verify w is a primitive n-th root of unity (where n = 2^k), 
        // we check that w^(n/2) == -1.
        // We compute w^(n/2) by squaring "w" (k-1) times.
        blst::Fr accum = w;
        for (int i = 0; i < k - 1; ++i) {
            blst::Fr::sqr(accum, accum);
        }

        EXPECT_TRUE(blst::Fr::equal(accum, negOne))
            << "For n=2^" << k << ", w^(n/2) != -1 (not primitive)";

        // Finally check w^n == 1 by squaring one more time
        blst::Fr::sqr(accum, accum);
        EXPECT_TRUE(blst::Fr::equal(accum, one))
            << "For n=2^" << k << ", w^n != 1 (not a root of unity)";

        EXPECT_TRUE(isPrimRootOfUnity<blst::Fr>(n, w))
			<< "isPrimRootOfUnity failed for n=2^" << k;
    }
}


Fr mf(u64 i)
{
    Fr r;
    Fr::from_uint64(r, i);
    return r;
}
//
//
//struct FrCoeffCtx
//{
//    static auto plus(const Fr& a, const Fr& b) {
//        Fr r;
//        Fr::add(r, a, b);
//        return r;
//    }
//    static auto minus(const Fr& a, const Fr& b) {
//        Fr r;
//        Fr::sub(r, a, b);
//        return r;
//    }
//
//    static auto mult(const Fr& a, const Fr& b) {
//        Fr r;
//        Fr::mul(r, a, b);
//        return r;
//    }
//    static auto inverse(const Fr& a) {
//        Fr r;
//        Fr::reciprocal(r, a);
//        return r;
//    }
//
//    template<typename F>
//    static auto from_u64(u64 a) {
//        return mf(a);
//    }
//
//};

std::ostream & operator<<(std::ostream& os, const Fr& f)
{
    char buf[65];
    Fr::to_hexascii(buf, f);
    os << buf;
    return os;
}

TEST(Fft, FrFft) {
    auto n = 4;
    size_t N = n * 2; // FFT size

    // 1. Setup data using the wrapper
    std::vector<Fr> A = { mf(1), mf(2), mf(3), mf(4) };
    std::vector<Fr> B = { mf(5), mf(6), mf(7), mf(8) };

    // 2. Naive Convolution for expected result
    std::vector<Fr> expected = naive_multiply<Fr, Fr, Fr, FrCoeffCtx>(A, B);
    Fr zero;
    memset(&zero, 0, sizeof(Fr));
    expected.resize(N, zero); // Pad to N

    // 3. FFT Setup
    // Get primitive root using the existing specialization for blst::Fr
    Fr prim = primRootOfUnity<Fr>(N);

    //EXPECT_TRUE(isPrimRootOfUnity<Fr>(N, prim));
        
    // Manually precompute roots to avoid depending on BasicFp generic implementation details
    std::vector<Fr> roots(N);
    Fr current_root = mf(1);
    for (size_t i = 0; i < N; ++i) {
        roots[i] = current_root;
        Fr::mul(current_root, current_root, prim);
    }

    std::vector<Fr> inv_roots(N);
    for (size_t i = 0; i < N; ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }

    // 4. Perform FFT
    A.resize(N, zero);
    B.resize(N, zero);

    fft<Fr,Fr, FrCoeffCtx>(A, roots, false);
    fft<Fr, Fr, FrCoeffCtx>(B, roots, false);

    std::vector<Fr> C(N);
    for (size_t i = 0; i < N; ++i) {
        Fr::mul(C[i], A[i], B[i]); // Pointwise multiplication
    }

    fft<Fr, Fr, FrCoeffCtx>(C, inv_roots, true); // Inverse FFT

    if (C != expected)
    {
        for (u64 i = 0; i < N; ++i)
        {
            std::cout << "exp[" << i << "] = " << expected[i] << std::endl;
        }
        for (u64 i = 0; i < N; ++i)
        {
            std::cout << "act[" << i << "] = " << C[i] << std::endl;
        }
        EXPECT_EQ(C, expected);
    }
}

// Helper to create G1 point = generator * scalar
G1Point g1(u64 scalar) {
    G1Point r;
    Fr s;
    Fr::from_uint64(s, scalar);
    Fr s_non_mont;
    Fr::from_mont(s_non_mont, s);
    G1Point::mult(r, G1Point::generator(), reinterpret_cast<const byte*>(&s_non_mont), 256);
    return r;
}

// Context for naive multiply: Fr * G1 -> G1
struct FrG1CoeffCtx
{
    static auto plus(const G1Point& a, const G1Point& b) {
        G1Point r;
        G1Point::dadd(r, a, b);
        return r;
    }

    static auto mult(const Fr& a, const G1Point& b) {
        G1Point r;
        Fr a_non_mont;
        Fr::from_mont(a_non_mont, a);
        G1Point::mult(r, b, reinterpret_cast<const byte*>(&a_non_mont), 256);
        return r;
    }
};

TEST(Fft, G1Fft) {
    auto n = 4;
    size_t N = n * 2; // FFT size

    // 1. Setup data: G1 points (A) and Fr scalars (B)
    // A[i] = g1(i+1) where g1(x) = generator * x
    std::vector<G1Point> A = { g1(1), g1(2), g1(3), g1(4) };
    std::vector<Fr> B = { mf(5), mf(6), mf(7), mf(8) };

    // 2. Naive Convolution for expected result
    // Computes sum of B[i] * A[j] for i+j = k
    // Result is in G1 since Fr * G1 -> G1
    std::vector<G1Point> expected = naive_multiply<Fr, G1Point, G1Point, FrG1CoeffCtx>(B, A);

    // Pad expected to N with identity (point at infinity)
    G1Point identity = G1Point::infinity();
    expected.resize(N, identity);

    // 3. FFT Setup
    Fr prim = primRootOfUnity<Fr>(N);

    // Precompute roots
    std::vector<Fr> roots(N);
    Fr current_root = Fr::one();
    for (size_t i = 0; i < N; ++i) {
        roots[i] = current_root;
        Fr::mul(current_root, current_root, prim);
    }

    std::vector<Fr> inv_roots(N);
    for (size_t i = 0; i < N; ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }

    // 4. Pad input arrays
    A.resize(N, identity);
    Fr zero;
    memset(&zero, 0, sizeof(Fr));
    B.resize(N, zero);

    // 5. Perform FFT
    // FFT on G1 points (in the exponent)
    fft<G1Point, Fr, G1CoeffCtx>(A, roots, false);
    // FFT on Fr scalars
    fft<Fr, Fr, FrCoeffCtx>(B, roots, false);

    // 6. Pointwise multiplication: B[i] * A[i] -> G1
    std::vector<G1Point> C(N);
    for (size_t i = 0; i < N; ++i) {
        C[i] = G1CoeffCtx::mult(A[i], B[i]);
    }

    // 7. Inverse FFT on G1 points
    fft<G1Point, Fr, G1CoeffCtx>(C, inv_roots, true);

    // 8. Compare results
    bool match = true;
    for (size_t i = 0; i < N; ++i) {
        if (!G1Point::equal(C[i], expected[i])) {
            match = false;
            break;
        }
    }

    if (!match) {
        std::cout << "G1 FFT Convolution FAILED\n";
        for (size_t i = 0; i < N; ++i) {
            byte exp_bytes[96], act_bytes[96];
            G1Point::serialize(exp_bytes, expected[i]);
            G1Point::serialize(act_bytes, C[i]);
            std::cout << "exp[" << i << "]: ";
            for (int j = 0; j < 12; ++j) printf("%02x", exp_bytes[j]);
            std::cout << "...\nact[" << i << "]: ";
            for (int j = 0; j < 12; ++j) printf("%02x", act_bytes[j]);
            std::cout << "...\n";
        }
    }

    EXPECT_TRUE(match);
}

// Helper to create G2 point = generator * scalar
G2Point g2(u64 scalar) {
    G2Point r;
    Fr s;
    Fr::from_uint64(s, scalar);
    Fr s_non_mont;
    Fr::from_mont(s_non_mont, s);
    G2Point::mult(r, G2Point::generator(), reinterpret_cast<const byte*>(&s_non_mont), 256);
    return r;
}

// Context for naive multiply: Fr * G2 -> G2
struct FrG2CoeffCtx
{
    static auto plus(const G2Point& a, const G2Point& b) {
        G2Point r;
        G2Point::dadd(r, a, b);
        return r;
    }

    static auto mult(const Fr& a, const G2Point& b) {
        G2Point r;
        Fr a_non_mont;
        Fr::from_mont(a_non_mont, a);
        G2Point::mult(r, b, reinterpret_cast<const byte*>(&a_non_mont), 256);
        return r;
    }
};

TEST(Fft, G2Fft) {
    auto n = 4;
    size_t N = n * 2; // FFT size

    // 1. Setup data: G2 points (A) and Fr scalars (B)
    // A[i] = g2(i+1) where g2(x) = generator * x
    std::vector<G2Point> A = { g2(1), g2(2), g2(3), g2(4) };
    std::vector<Fr> B = { mf(5), mf(6), mf(7), mf(8) };

    // 2. Naive Convolution for expected result
    // Computes sum of B[i] * A[j] for i+j = k
    // Result is in G2 since Fr * G2 -> G2
    std::vector<G2Point> expected = naive_multiply<Fr, G2Point, G2Point, FrG2CoeffCtx>(B, A);

    // Pad expected to N with identity (point at infinity)
    G2Point identity = G2Point::infinity();
    expected.resize(N, identity);

    // 3. FFT Setup
    Fr prim = primRootOfUnity<Fr>(N);

    // Precompute roots
    std::vector<Fr> roots(N);
    Fr current_root = Fr::one();
    for (size_t i = 0; i < N; ++i) {
        roots[i] = current_root;
        Fr::mul(current_root, current_root, prim);
    }

    std::vector<Fr> inv_roots(N);
    for (size_t i = 0; i < N; ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }

    // 4. Pad input arrays
    A.resize(N, identity);
    Fr zero;
    memset(&zero, 0, sizeof(Fr));
    B.resize(N, zero);

    // 5. Perform FFT
    // FFT on G2 points (in the exponent)
    fft<G2Point, Fr, G2CoeffCtx>(A, roots, false);
    // FFT on Fr scalars
    fft<Fr, Fr, FrCoeffCtx>(B, roots, false);

    // 6. Pointwise multiplication: B[i] * A[i] -> G2
    std::vector<G2Point> C(N);
    for (size_t i = 0; i < N; ++i) {
        C[i] = G2CoeffCtx::mult(A[i], B[i]);
    }

    // 7. Inverse FFT on G2 points
    fft<G2Point, Fr, G2CoeffCtx>(C, inv_roots, true);

    // 8. Compare results
    bool match = true;
    for (size_t i = 0; i < N; ++i) {
        if (!G2Point::equal(C[i], expected[i])) {
            match = false;
            break;
        }
    }

    if (!match) {
        std::cout << "G2 FFT Convolution FAILED\n";
        for (size_t i = 0; i < N; ++i) {
            byte exp_bytes[192], act_bytes[192];
            G2Point::serialize(exp_bytes, expected[i]);
            G2Point::serialize(act_bytes, C[i]);
            std::cout << "exp[" << i << "]: ";
            for (int j = 0; j < 12; ++j) printf("%02x", exp_bytes[j]);
            std::cout << "...\nact[" << i << "]: ";
            for (int j = 0; j < 12; ++j) printf("%02x", act_bytes[j]);
            std::cout << "...\n";
        }
    }

    EXPECT_TRUE(match);
}




// Pairing function: G1 x G2 -> GT
GT pair(const G1Point& p, const G2Point& q) {
    G1AffinePoint p_aff;
    G1Point::to_affine(p_aff, p);
    G2AffinePoint q_aff;
    G2Point::to_affine(q_aff, q);

    GT result;
    GT::miller_loop(result, &q_aff, &p_aff);
    GT::final_mult(result, result);
    return result;
}

// Context for naive multiply: pairing(G1, G2) -> GT
struct G1G2PairingCoeffCtx
{
    static auto plus(const GT& a, const GT& b) {
        GT r;
        GT::add(r, a, b);
        return r;
    }

    static auto mult(const G1Point& a, const G2Point& b) {
        return pair(a, b);
    }

    // Provide identity for proper initialization
    static GT identity() {
        return GT::identity();
    }
};
// Test that basic pairing bilinearity works
TEST(Fft, PairingBilinearity) {
    // Test: e(aP, bQ) = e(P, Q)^(ab)
    // Using the pattern from Pairing_tests.cpp

    G1Point P = G1Point::generator();
    G2Point Q = G2Point::generator();

    // Use a=2, b=3 so ab=6
    // Compute 2*P and 3*Q
    G1Point P2;
    G1Point::doubling(P2, P);  // 2*P

    G2Point Q3;
    G2Point::doubling(Q3, Q);   // 2*Q
    G2Point::dadd(Q3, Q3, Q);   // 3*Q

    // Convert to affine
    G1AffinePoint P2_aff;
    G2AffinePoint Q3_aff;
    G1Point::to_affine(P2_aff, P2);
    G2Point::to_affine(Q3_aff, Q3);

    G1AffinePoint P_aff;
    G2AffinePoint Q_aff;
    G1Point::to_affine(P_aff, P);
    G2Point::to_affine(Q_aff, Q);

    // Compute e(2P, 3Q)
    GT left;
    GT::miller_loop(left, &Q3_aff, &P2_aff);
    GT::final_mult(left, left);

    // Compute e(P, Q)^6 by repeated doubling/multiplication
    // e(P,Q)^6 = e(P,Q)^4 * e(P,Q)^2 = ((e(P,Q)^2)^2) * (e(P,Q)^2)
    GT base;
    GT::miller_loop(base, &Q_aff, &P_aff);
    GT::final_mult(base, base);

    // e(P,Q)^2
    GT base2;
    GT::doubling(base2, base);

    // e(P,Q)^4
    GT base4;
    GT::doubling(base4, base2);

    // e(P,Q)^6 = e(P,Q)^4 * e(P,Q)^2
    GT right;
    GT::add(right, base4, base2);

    EXPECT_TRUE(GT::equal(left, right)) << "Pairing bilinearity failed: e(2P,3Q) != e(P,Q)^6";
}

// Test GT FFT round-trip (FFT then iFFT should give back original)
TEST(Fft, GTFftRoundTrip) {
    size_t N = 16;

    // Create some GT elements via pairings - use the correct pairing pattern
    G1AffinePoint P_aff;
    G2AffinePoint Q_aff;
    G1Point::to_affine(P_aff, G1Point::generator());
    G2Point::to_affine(Q_aff, G2Point::generator());

    // Base pairing
    GT base;
    GT::miller_loop(base, &Q_aff, &P_aff);
    GT::final_mult(base, base);

    std::vector<GT> original(N);
    original[0] = GT::identity();

    // Create distinct GT elements by repeated "addition" (Fp12 multiplication)
    GT accum = base;
    for (size_t i = 1; i < N; ++i) {
        original[i] = accum;
        GT::add(accum, accum, base);  // accum = accum * base in Fp12
    }

    // Copy for FFT
    std::vector<GT> data = original;

    // FFT Setup
    Fr prim = primRootOfUnity<Fr>(N);
    std::vector<Fr> roots(N);
    Fr current_root = Fr::one();
    for (size_t i = 0; i < N; ++i) {
        roots[i] = current_root;
        Fr::mul(current_root, current_root, prim);
    }

    std::vector<Fr> inv_roots(N);
    for (size_t i = 0; i < N; ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }

    // FFT forward
    fft<GT, Fr, GTCoeffCtx>(data, roots, false);

    // FFT inverse
    fft<GT, Fr, GTCoeffCtx>(data, inv_roots, true);

    // Compare
    bool match = true;
    for (size_t i = 0; i < N; ++i) {
        if (!GT::equal(data[i], original[i])) {
            match = false;
            std::cout << "Mismatch at index " << i << std::endl;
            break;
        }
    }

    EXPECT_TRUE(match) << "GT FFT round-trip failed";
}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
TEST(Fft, GTAvx52InverseFftMatchesScalar) {
    const size_t N = 16;

    G1AffinePoint p_aff = {};
    G2AffinePoint q_aff = {};
    G1Point::to_affine(p_aff, G1Point::generator());
    G2Point::to_affine(q_aff, G2Point::generator());

    GT base = {};
    GT::miller_loop(base, &q_aff, &p_aff);
    GT::final_mult(base, base);

    std::vector<GT> time_domain(N);
    time_domain[0] = GT::identity();
    GT accum = base;
    for (size_t i = 1; i < N; ++i) {
        time_domain[i] = accum;
        GT::add(accum, accum, base);
    }

    Fr prim = primRootOfUnity<Fr>(N);
    std::vector<Fr> roots(N);
    Fr current_root = Fr::one();
    for (size_t i = 0; i < N; ++i) {
        roots[i] = current_root;
        Fr::mul(current_root, current_root, prim);
    }

    std::vector<Fr> inv_roots(N);
    for (size_t i = 0; i < N; ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }

    std::vector<GT> scalar_freq_domain = time_domain;
    fft<GT, Fr, GTCoeffCtx>(scalar_freq_domain, roots, false);

    std::vector<blst::avx::AvxGt52> avx_freq_domain =
        blst::avx::scalar_gt_vector_to_avx52(std::span<const GT>(scalar_freq_domain));
    blst::avx::AvxGt52PreparedRoots prepared_inv_roots{std::span<const Fr>(inv_roots)};
    prepared_inv_roots.set_scale_from_size(N);
    blst::avx::AvxGt52CoeffCtx::ScopedRoots scoped_roots(prepared_inv_roots);

    fft<blst::avx::AvxGt52, Fr, blst::avx::AvxGt52CoeffCtx>(avx_freq_domain, inv_roots, true);

    const std::vector<GT> avx_round_trip =
        blst::avx::avx52_vector_to_scalar_gt(std::span<const blst::avx::AvxGt52>(avx_freq_domain));

    bool match = true;
    for (size_t i = 0; i < N; ++i) {
        if (!GT::equal(avx_round_trip[i], time_domain[i])) {
            match = false;
            std::cout << "AVX GT inverse FFT mismatch at index " << i << std::endl;
            break;
        }
    }

    EXPECT_TRUE(match) << "AVX GT inverse FFT round-trip failed";
}
#endif

// Test that GT scalar multiplication (exponentiation) distributes correctly
TEST(Fft, GTScalarMult) {
    GT a = pair(g1(2), g2(3));
    GT b = pair(g1(4), g2(5));

    Fr scalar;
    Fr::from_uint64(scalar, 7);

    // (a + b) * scalar should equal a*scalar + b*scalar
    // In GT terms: (a * b)^scalar = a^scalar * b^scalar
    GT sum;
    GT::add(sum, a, b);
    GT left = GTCoeffCtx::mult(sum, scalar);

    GT a_scaled = GTCoeffCtx::mult(a, scalar);
    GT b_scaled = GTCoeffCtx::mult(b, scalar);
    GT right;
    GT::add(right, a_scaled, b_scaled);

    EXPECT_TRUE(GT::equal(left, right)) << "GT scalar multiplication distributivity failed";
}

TEST(Fft, G1G2PairingFft) {
    size_t n = 8;
    size_t N = n * 2; // FFT size

    // 1. Setup data: G1 points (A) and G2 points (B)
    //std::vector<G1Point> A = { g1(1), g1(2), g1(3), g1(4) };
    //std::vector<G2Point> B = { g2(5), g2(6), g2(7), g2(8) };

    std::vector<G1Point> A;
    std::vector<G2Point> B;
    while (A.size() != n)
    {
        A.push_back(g1(rand()));
        B.push_back(g2(rand()));
    }

    // 2. Naive Convolution for expected result
    // Computes sum of pairing(A[i], B[j]) for i+j = k
    // Result is in GT since pairing(G1, G2) -> GT
    // Must initialize with GT identity (Fp12::one()), not default zeros
    size_t m = A.size();
    size_t nn = B.size();
    std::vector<GT> expected(m + nn - 1, GT::identity());
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < nn; ++j) {
            GT p = pair(A[i], B[j]);
            GT::add(expected[i + j], expected[i + j], p);
        }
    }

    // Pad expected to N with identity
    GT gt_identity = GT::identity();
    expected.resize(N, gt_identity);

    // 3. FFT Setup
    Fr prim = primRootOfUnity<Fr>(N);

    // Precompute roots
    std::vector<Fr> roots(N);
    Fr current_root = Fr::one();
    for (size_t i = 0; i < N; ++i) {
        roots[i] = current_root;
        Fr::mul(current_root, current_root, prim);
    }

    std::vector<Fr> inv_roots(N);
    for (size_t i = 0; i < N; ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }

    // 4. Pad input arrays
    G1Point g1_identity = G1Point::infinity();
    G2Point g2_identity = G2Point::infinity();
    A.resize(N, g1_identity);
    B.resize(N, g2_identity);

    // 5. Perform FFT on G1 and G2 points (in the exponent)
    fft<G1Point, Fr, G1CoeffCtx>(A, roots, false);
    fft<G2Point, Fr, G2CoeffCtx>(B, roots, false);

    // 6. Pointwise pairing: e(A[i], B[i]) -> GT
    std::vector<GT> C(N);
    for (size_t i = 0; i < N; ++i) {
        C[i] = pair(A[i], B[i]);
    }

    // 7. Inverse FFT on GT elements
    fft<GT, Fr, GTCoeffCtx>(C, inv_roots, true);

    // 8. Compare results
    bool match = true;
    for (size_t i = 0; i < N; ++i) {
        if (!GT::equal(C[i], expected[i])) {
            match = false;
            break;
        }
    }

    if (!match) {
        std::cout << "G1G2 Pairing FFT Convolution FAILED\n";
        for (size_t i = 0; i < N; ++i) {
            byte exp_bytes[576], act_bytes[576];
            GT::serialize(exp_bytes, expected[i]);
            GT::serialize(act_bytes, C[i]);
            std::cout << "exp[" << i << "]: ";
            for (int j = 0; j < 12; ++j) printf("%02x", exp_bytes[j]);
            std::cout << "...\nact[" << i << "]: ";
            for (int j = 0; j < 12; ++j) printf("%02x", act_bytes[j]);
            std::cout << "...\n";
        }
    }

    EXPECT_TRUE(match);
}

// Context for poly_mult_fft with G1, G2 -> GT via pairing
struct PairingCoeffCtx
{
    // G1 operations
    template<typename G, typename F>
    static auto plus(const G& a, const F& b) {
        G r;
        G::dadd(r, a, b);
        return r;
    }

    template<typename G, typename F>
    static auto minus(const G& a, const F& b) {
        G r;
        F neg_b;
        F::neg(neg_b, b);
        G::dadd(r, a, neg_b);
        return r;
    }

    // Scalar multiplication for G1/G2/GT
    template<typename G, typename F>
    static auto mult(const G& a, const F& b) {

        if constexpr (std::is_same<G, F>::value) {
            // Scalar * Scalar case (Fr * Fr)
            G r;
            G::mul(r, a, b);
            return r;
        }
        else {
            G r;
            F b_non_mont;
            F::from_mont(b_non_mont, b);
            G::mult(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
            return r;
        }
    }

    // Pairing multiplication: G1 x G2 -> GT
    static auto mult(const G1Point& a, const G2Point& b) {
        return pair(a, b);
    }

    static auto inverse(const Fr& a) {
        Fr r;
        Fr::reciprocal(r, a);
        return r;
    }

    template<typename F>
    static auto from_u64(u64 a) {
        return mf(a);
    }
};
TEST(Fft, PolyMultFft) {
    size_t n = 8;
    size_t N = n * 2; // FFT size

    // 1. Setup data: G1 points (A) and G2 points (B)
    std::vector<G1Point> A;
    std::vector<G2Point> B;
    for (size_t i = 0; i < n; ++i) {
        A.push_back(g1(rand()));
        B.push_back(g2(rand()));
    }

    // 2. Naive Convolution for expected result
    size_t m = A.size();
    size_t nn = B.size();
    std::vector<GT> expected(m + nn - 1, GT::identity());
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < nn; ++j) {
            GT p = pair(A[i], B[j]);
            GT::add(expected[i + j], expected[i + j], p);
        }
    }
    expected.resize(N, GT::identity());

    // 3. FFT Setup
    Fr prim = primRootOfUnity<Fr>(N);

    std::vector<Fr> roots(N);
    Fr current_root = Fr::one();
    for (size_t i = 0; i < N; ++i) {
        roots[i] = current_root;
        Fr::mul(current_root, current_root, prim);
    }

    std::vector<Fr> inv_roots(N);
    for (size_t i = 0; i < N; ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }

    // 4. Pad input arrays to FFT size
    A.resize(N, G1Point::infinity());
    B.resize(N, G2Point::infinity());
    std::vector<GT> C(N, GT::identity());

    std::vector<G2Point> fftB = B;
	// FFT on G2 points
	fft<G2Point, Fr, G2CoeffCtx>(fftB, roots, false);
    
    // 5. Use poly_mult_fft
    poly_mult_fft<G1Point, G2Point, GT, Fr, PairingCoeffCtx>(A, fftB, true, C, roots, inv_roots);

    // 6. Compare results
    bool match = true;
    for (size_t i = 0; i < N; ++i) {
        if (!GT::equal(C[i], expected[i])) {
            match = false;
            break;
        }
    }

    if (!match) {
        std::cout << "PolyMultFft FAILED\n";
        for (size_t i = 0; i < N; ++i) {
            byte exp_bytes[576], act_bytes[576];
            GT::serialize(exp_bytes, expected[i]);
            GT::serialize(act_bytes, C[i]);
            std::cout << "exp[" << i << "]: ";
            for (int j = 0; j < 12; ++j) printf("%02x", exp_bytes[j]);
            std::cout << "...\nact[" << i << "]: ";
            for (int j = 0; j < 12; ++j) printf("%02x", act_bytes[j]);
            std::cout << "...\n";
        }
    }

    EXPECT_TRUE(match);
}


#if 0  // Legacy middle_product_fft/ifft helper coverage (API removed)
// Test middle_product_fft against naive convolution
TEST(Fft, MiddleProductFft) {
    std::size_t n = 8;
    std::size_t m = 16;  // FFT size, must be >= 2n

    // 1. Setup data: G1 points (A) with n elements, G2 points (B) with 2n elements
    std::vector<G1Point> A_orig;
    std::vector<G2Point> B_orig;

    srand(42);
    for (std::size_t i = 0; i < n; ++i) {
        A_orig.push_back(g1(rand() % 1000 + 1));
    }
    for (std::size_t i = 0; i < 2 * n; ++i) {
        B_orig.push_back(g2(rand() % 1000 + 1));
    }

    // 2. Naive Convolution for expected result
    std::size_t full_size = n + 2 * n - 1;
    std::vector<GT> full_product(full_size, GT::identity());
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < 2 * n; ++j) {
            GT p = pair(A_orig[i], B_orig[j]);
            GT::add(full_product[i + j], full_product[i + j], p);
        }
    }

    // Expected middle n coefficients: c[n], c[n+1], ..., c[2n-1]
    std::vector<GT> expected(n);
    for (std::size_t k = 0; k < n; ++k) {
        expected[k] = full_product[n + k];
    }

    // 3. FFT Setup
    Fr prim = primRootOfUnity<Fr>(m);
    std::vector<Fr> roots(m);
    Fr current_root = Fr::one();
    for (std::size_t i = 0; i < m; ++i) {
        roots[i] = current_root;
        Fr::mul(current_root, current_root, prim);
    }
    std::vector<Fr> inv_roots(m);
    for (std::size_t i = 0; i < m; ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }

    // 4. Prepare buffers (caller allocates everything)
    std::vector<G1Point> A(m, G1Point::infinity());
    std::vector<G2Point> B0(m, G2Point::infinity());
    std::vector<G2Point> B1(m, G2Point::infinity());
    std::vector<GT> out(n, GT::identity());
    std::vector<GT> C_scratch(m, GT::identity());

    // Copy input data
    for (std::size_t i = 0; i < n; ++i) {
        A[i] = A_orig[i];
        B0[i] = B_orig[i];        // B0 = B[0..n-1]
        B1[i] = B_orig[n + i];    // B1 = B[n..2n-1]
    }

    // 5. Call middle_product_fft
    middle_product_fft<G1Point, G2Point, GT, Fr, PairingCoeffCtx>(
        std::span<G1Point>(A), std::span<G2Point>(B0), std::span<G2Point>(B1),
        false,
        n, m,
        std::span<GT>(out), std::span<GT>(C_scratch),
        std::span<const Fr>(roots.data(), m), std::span<const Fr>(inv_roots.data(), m),
        prim,
        GT::identity(), G1Point::infinity(), G2Point::infinity());

    // 6. Compare results
    bool match = true;
    for (std::size_t i = 0; i < n; ++i) {
        if (!GT::equal(out[i], expected[i])) {
            match = false;
            break;
        }
    }

    if (!match) {
        std::cout << "MiddleProductFft FAILED\n";
        for (std::size_t i = 0; i < n; ++i) {
            byte exp_bytes[576], act_bytes[576];
            GT::serialize(exp_bytes, expected[i]);
            GT::serialize(act_bytes, out[i]);
            std::cout << "exp[" << i << "] (c[" << (n + i) << "]): ";
            for (int j = 0; j < 12; ++j) printf("%02x", exp_bytes[j]);
            std::cout << "...\nact[" << i << "]: ";
            for (int j = 0; j < 12; ++j) printf("%02x", act_bytes[j]);
            std::cout << "...\n";
        }
    }

    EXPECT_TRUE(match);
}

// Test middle_product_fft with different sizes
TEST(Fft, MiddleProductFftSizes) {
    std::vector<std::size_t> test_sizes = { 2, 4, 8, 16 };

    for (std::size_t n : test_sizes) {
        std::size_t m = 1;
        while (m < 2 * n) m <<= 1;

        // Setup data
        std::vector<G1Point> A_orig;
        std::vector<G2Point> B_orig;

        srand(static_cast<unsigned int>(n * 123));
        for (std::size_t i = 0; i < n; ++i) {
            A_orig.push_back(g1(rand() % 100 + 1));
        }
        for (std::size_t i = 0; i < 2 * n; ++i) {
            B_orig.push_back(g2(rand() % 100 + 1));
        }

        // Naive computation
        std::vector<GT> expected(n, GT::identity());
        for (std::size_t k = 0; k < n; ++k) {
            for (std::size_t i = 0; i < n; ++i) {
                std::size_t j = n + k - i;
                if (j < 2 * n) {
                    GT p = pair(A_orig[i], B_orig[j]);
                    GT::add(expected[k], expected[k], p);
                }
            }
        }

        // FFT Setup
        Fr prim = primRootOfUnity<Fr>(m);
        std::vector<Fr> roots(m), inv_roots(m);
        Fr current = Fr::one();
        for (std::size_t i = 0; i < m; ++i) {
            roots[i] = current;
            Fr::mul(current, current, prim);
        }
        for (std::size_t i = 0; i < m; ++i) {
            Fr::reciprocal(inv_roots[i], roots[i]);
        }

        // Prepare buffers
        std::vector<G1Point> A(m, G1Point::infinity());
        std::vector<G2Point> B0(m, G2Point::infinity());
        std::vector<G2Point> B1(m, G2Point::infinity());
        std::vector<GT> out(n, GT::identity());
        std::vector<GT> C_scratch(m, GT::identity());

        for (std::size_t i = 0; i < n; ++i) {
            A[i] = A_orig[i];
            B0[i] = B_orig[i];
            B1[i] = B_orig[n + i];
        }

        // Call
        middle_product_fft<G1Point, G2Point, GT, Fr, PairingCoeffCtx>(
            std::span<G1Point>(A), std::span<G2Point>(B0), std::span<G2Point>(B1),
            false,
            n, m,
            std::span<GT>(out), std::span<GT>(C_scratch),
            std::span<const Fr>(roots.data(), m), std::span<const Fr>(inv_roots.data(), m),
            prim,
            GT::identity(), G1Point::infinity(), G2Point::infinity());

        // Compare
        bool match = true;
        for (std::size_t i = 0; i < n; ++i) {
            if (!GT::equal(out[i], expected[i])) {
                match = false;
                std::cout << "Mismatch at n=" << n << ", index=" << i << std::endl;
                break;
            }
        }

        EXPECT_TRUE(match) << "MiddleProductFft failed for n=" << n;
    }
}

// Test that middle_product_fft matches extracting middle from full product
TEST(Fft, MiddleProductVsFullProduct) {
    std::size_t n = 4;
    std::size_t m = 8;

    // Setup data
    std::vector<G1Point> A_orig;
    std::vector<G2Point> B_orig;

    srand(999);
    for (std::size_t i = 0; i < n; ++i) {
        A_orig.push_back(g1(rand() % 50 + 1));
    }
    for (std::size_t i = 0; i < 2 * n; ++i) {
        B_orig.push_back(g2(rand() % 50 + 1));
    }

    // Compute full product via FFT
    std::size_t full_fft_size = 1;
    while (full_fft_size < 3 * n - 1) full_fft_size <<= 1;

    Fr prim_full = primRootOfUnity<Fr>(full_fft_size);
    std::vector<Fr> roots_full(full_fft_size), inv_roots_full(full_fft_size);
    Fr current = Fr::one();
    for (std::size_t i = 0; i < full_fft_size; ++i) {
        roots_full[i] = current;
        Fr::mul(current, current, prim_full);
    }
    for (std::size_t i = 0; i < full_fft_size; ++i) {
        Fr::reciprocal(inv_roots_full[i], roots_full[i]);
    }

    std::vector<G1Point> A_full(full_fft_size, G1Point::infinity());
    std::vector<G2Point> B_full(full_fft_size, G2Point::infinity());
    for (std::size_t i = 0; i < n; ++i) A_full[i] = A_orig[i];
    for (std::size_t i = 0; i < 2 * n; ++i) B_full[i] = B_orig[i];

    fft<G1Point, Fr, G1CoeffCtx>(A_full, roots_full, false);
    fft<G2Point, Fr, G2CoeffCtx>(B_full, roots_full, false);

    std::vector<GT> C_full(full_fft_size);
    for (std::size_t i = 0; i < full_fft_size; ++i) {
        C_full[i] = pair(A_full[i], B_full[i]);
    }
    fft<GT, Fr, GTCoeffCtx>(C_full, inv_roots_full, true);

    // Expected from full product
    std::vector<GT> expected(n);
    for (std::size_t k = 0; k < n; ++k) {
        expected[k] = C_full[n + k];
    }

    // Compute via middle_product_fft
    Fr prim = primRootOfUnity<Fr>(m);
    std::vector<Fr> roots(m), inv_roots(m);
    current = Fr::one();
    for (std::size_t i = 0; i < m; ++i) {
        roots[i] = current;
        Fr::mul(current, current, prim);
    }
    for (std::size_t i = 0; i < m; ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }

    std::vector<G1Point> A(m, G1Point::infinity());
    std::vector<G2Point> B0(m, G2Point::infinity());
    std::vector<G2Point> B1(m, G2Point::infinity());
    std::vector<GT> out(n, GT::identity());
    std::vector<GT> C_scratch(m, GT::identity());

    for (std::size_t i = 0; i < n; ++i) {
        A[i] = A_orig[i];
        B0[i] = B_orig[i];
        B1[i] = B_orig[n + i];
    }

    middle_product_fft<G1Point, G2Point, GT, Fr, PairingCoeffCtx>(
        std::span<G1Point>(A), std::span<G2Point>(B0), std::span<G2Point>(B1),
        false,
        n, m,
        std::span<GT>(out), std::span<GT>(C_scratch),
        std::span<const Fr>(roots.data(), m), std::span<const Fr>(inv_roots.data(), m),
        prim,
        GT::identity(), G1Point::infinity(), G2Point::infinity());

    // Compare
    bool match = true;
    for (std::size_t i = 0; i < n; ++i) {
        if (!GT::equal(out[i], expected[i])) {
            match = false;
            std::cout << "Mismatch at index " << i << std::endl;
            break;
        }
    }

    EXPECT_TRUE(match) << "middle_product_fft does not match full product extraction";
}









// Test truncated inverse FFT correctness against full inverse FFT
TEST(Fft, TruncatedInverseFftCorrectness) {
    std::size_t m = 16;

    std::vector<std::size_t> test_lengths = { 1, 2, 4, 8, 12, 15, 16 };

    for (std::size_t L : test_lengths) {
        std::vector<GT> freq_full(m);
        std::vector<GT> freq_trunc(m);

        srand(42 + static_cast<unsigned int>(L));
        for (std::size_t i = 0; i < m; ++i) {
            freq_full[i] = pair(g1(rand() % 100 + 1), g2(rand() % 100 + 1));
            freq_trunc[i] = freq_full[i];
        }

        Fr prim = primRootOfUnity<Fr>(m);
        Fr m_inv = Fr::from_uint64(m);
        Fr::reciprocal(m_inv, m_inv);

        // Setup inv_roots (reciprocals of forward roots)
        std::vector<Fr> inv_roots(m);
        Fr current = Fr::one();
        for (std::size_t i = 0; i < m; ++i) {
            Fr::reciprocal(inv_roots[i], current);
            Fr::mul(current, current, prim);
        }

        // Full inverse FFT using standard fft with invert=true
        fft<GT, Fr, GTCoeffCtx>(freq_full, inv_roots, true);

        // Truncated inverse FFT
        fft_inverse_truncated<GT, Fr, GTCoeffCtx>(
            std::span<GT>(freq_trunc), m, L,
            std::span<const Fr>(inv_roots), m_inv);

        bool match = true;
        for (std::size_t i = 0; i < L; ++i) {
            if (!GT::equal(freq_full[i], freq_trunc[i])) {
                match = false;
                std::cout << "Mismatch at L=" << L << ", index=" << i << std::endl;
                break;
            }
        }

        EXPECT_TRUE(match) << "Truncated inverse FFT failed for L=" << L;
    }
}

// Test full inverse matches standard
TEST(Fft, FullInverseFftCorrectness) {
    std::size_t m = 16;

    std::vector<GT> freq_std(m);
    std::vector<GT> freq_new(m);

    srand(123);
    for (std::size_t i = 0; i < m; ++i) {
        freq_std[i] = pair(g1(rand() % 100 + 1), g2(rand() % 100 + 1));
        freq_new[i] = freq_std[i];
    }

    Fr prim = primRootOfUnity<Fr>(m);
    Fr m_inv = Fr::from_uint64(m);
    Fr::reciprocal(m_inv, m_inv);

    std::vector<Fr> inv_roots(m);
    Fr current = Fr::one();
    for (std::size_t i = 0; i < m; ++i) {
        Fr::reciprocal(inv_roots[i], current);
        Fr::mul(current, current, prim);
    }

    // Standard inverse
    fft<GT, Fr, GTCoeffCtx>(freq_std, inv_roots, true);

    // New inverse function
    fft_inverse_full<GT, Fr, GTCoeffCtx>(
        std::span<GT>(freq_new), m,
        std::span<const Fr>(inv_roots), m_inv);

    bool match = true;
    for (std::size_t i = 0; i < m; ++i) {
        if (!GT::equal(freq_std[i], freq_new[i])) {
            match = false;
            std::cout << "Mismatch at index " << i << std::endl;
            break;
        }
    }

    EXPECT_TRUE(match) << "Full inverse FFT does not match standard inverse";
}
#endif

TEST(Fft, TruncatedForwardPrefixMatchesFullFr)
{
    const std::size_t m = 32;
    Fr prim = primRootOfUnity<Fr>(m);
    std::vector<Fr> roots(m);
    Fr cur = Fr::one();
    for (std::size_t i = 0; i < m; ++i) {
        roots[i] = cur;
        Fr::mul(cur, cur, prim);
    }

    const std::vector<std::size_t> xns = { 1, 5, 13, 24, 32 };
    for (std::size_t xn : xns) {
        std::vector<Fr> full(m);
        std::vector<Fr> trunc(m);
        Fr zero;
        memset(&zero, 0, sizeof(Fr));

        for (std::size_t i = 0; i < m; ++i) {
            if (i < xn) {
                full[i] = mf(static_cast<u64>(i + 3));
            }
            else {
                full[i] = zero;
            }
        }
        trunc = full;

        fft<Fr, Fr, FrCoeffCtx>(full, roots, false);
        fft_forward_truncated_prefix<Fr, Fr, FrCoeffCtx>(trunc, roots, xn);

        EXPECT_EQ(full, trunc) << "Fr truncated forward mismatch at xn=" << xn;
    }
}

TEST(Fft, TruncatedForwardPrefixMatchesFullG1)
{
    const std::size_t m = 32;
    Fr prim = primRootOfUnity<Fr>(m);
    std::vector<Fr> roots(m);
    Fr cur = Fr::one();
    for (std::size_t i = 0; i < m; ++i) {
        roots[i] = cur;
        Fr::mul(cur, cur, prim);
    }

    const std::vector<std::size_t> xns = { 1, 5, 13, 24, 32 };
    for (std::size_t xn : xns) {
        std::vector<G1Point> full(m, G1Point::infinity());
        std::vector<G1Point> trunc(m, G1Point::infinity());

        for (std::size_t i = 0; i < xn; ++i) {
            full[i] = g1(static_cast<u64>(i + 2));
        }
        trunc = full;

        fft<G1Point, Fr, G1CoeffCtx>(full, roots, false);
        fft_forward_truncated_prefix<G1Point, Fr, G1CoeffCtx>(trunc, roots, xn);

        bool match = true;
        for (std::size_t i = 0; i < m; ++i) {
            if (!G1Point::equal(full[i], trunc[i])) {
                match = false;
                break;
            }
        }
        EXPECT_TRUE(match) << "G1 truncated forward mismatch at xn=" << xn;
    }
}

TEST(Fft, InverseTruncatedRangeMatchesFullMiddleWindow)
{
    const std::size_t m = 64;
    Fr prim = primRootOfUnity<Fr>(m);
    std::vector<Fr> roots(m);
    Fr cur = Fr::one();
    for (std::size_t i = 0; i < m; ++i) {
        roots[i] = cur;
        Fr::mul(cur, cur, prim);
    }
    std::vector<Fr> inv_roots(m);
    for (std::size_t i = 0; i < m; ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }

    std::vector<GT> freq_full(m, GT::identity());
    srand(12345);
    for (std::size_t i = 0; i < m; ++i) {
        freq_full[i] = pair(g1(rand() % 100 + 1), g2(rand() % 100 + 1));
    }
    std::vector<GT> freq_trunc = freq_full;

    const std::size_t start = 17;
    const std::size_t len = 21;

    fft<GT, Fr, GTCoeffCtx>(freq_full, inv_roots, true);
    fft_inverse_truncated_range<GT, Fr, GTCoeffCtx>(
        std::span<GT>(freq_trunc),
        std::span<const Fr>(inv_roots),
        start,
        len);

    bool match = true;
    for (std::size_t i = 0; i < len; ++i) {
        const std::size_t idx = start + i;
        if (!GT::equal(freq_full[idx], freq_trunc[idx])) {
            match = false;
            break;
        }
    }
    EXPECT_TRUE(match) << "Inverse range truncation mismatch";
}

#if 0  // Legacy middle_product_fft_truncated coverage (API removed)
// Test middle_product_fft_truncated matches middle_product_fft
TEST(Fft, MiddleProductFftTruncatedCorrectness) {
    std::size_t n = 8;
    std::size_t m = 16;

    std::vector<G1Point> A_orig;
    std::vector<G2Point> B_orig;

    srand(42);
    for (std::size_t i = 0; i < n; ++i) {
        A_orig.push_back(g1(rand() % 1000 + 1));
    }
    for (std::size_t i = 0; i < 2 * n; ++i) {
        B_orig.push_back(g2(rand() % 1000 + 1));
    }

    // FFT Setup
    Fr prim = primRootOfUnity<Fr>(m);
    Fr m_inv = Fr::from_uint64(m);
    Fr::reciprocal(m_inv, m_inv);

    std::vector<Fr> roots(m);
    Fr current = Fr::one();
    for (std::size_t i = 0; i < m; ++i) {
        roots[i] = current;
        Fr::mul(current, current, prim);
    }
    std::vector<Fr> inv_roots(m);
    for (std::size_t i = 0; i < m; ++i) {
        Fr::reciprocal(inv_roots[i], roots[i]);
    }

    // Standard version
    std::vector<G1Point> A_std(m, G1Point::infinity());
    std::vector<G2Point> B0_std(m, G2Point::infinity());
    std::vector<G2Point> B1_std(m, G2Point::infinity());
    std::vector<GT> out_std(n, GT::identity());
    std::vector<GT> C_std(m, GT::identity());

    for (std::size_t i = 0; i < n; ++i) {
        A_std[i] = A_orig[i];
        B0_std[i] = B_orig[i];
        B1_std[i] = B_orig[n + i];
    }

    middle_product_fft<G1Point, G2Point, GT, Fr, PairingCoeffCtx>(
        std::span<G1Point>(A_std), std::span<G2Point>(B0_std), std::span<G2Point>(B1_std),
        false,
        n, m,
        std::span<GT>(out_std), std::span<GT>(C_std),
        std::span<const Fr>(roots), std::span<const Fr>(inv_roots),
        prim,
        GT::identity(), G1Point::infinity(), G2Point::infinity());

    // Truncated version
    std::vector<G1Point> A_trunc(m, G1Point::infinity());
    std::vector<G2Point> B0_trunc(m, G2Point::infinity());
    std::vector<G2Point> B1_trunc(m, G2Point::infinity());
    std::vector<GT> out_trunc(n, GT::identity());
    std::vector<GT> C_trunc(m, GT::identity());

    for (std::size_t i = 0; i < n; ++i) {
        A_trunc[i] = A_orig[i];
        B0_trunc[i] = B_orig[i];
        B1_trunc[i] = B_orig[n + i];
    }

    middle_product_fft_truncated<G1Point, G2Point, GT, Fr, PairingCoeffCtx>(
        std::span<G1Point>(A_trunc), std::span<G2Point>(B0_trunc), std::span<G2Point>(B1_trunc),
        false,
        n, m,
        std::span<GT>(out_trunc), std::span<GT>(C_trunc),
        std::span<const Fr>(roots), std::span<const Fr>(inv_roots),
        prim, m_inv,
        GT::identity(), G1Point::infinity(), G2Point::infinity());

    // Compare
    bool match = true;
    for (std::size_t i = 0; i < n; ++i) {
        if (!GT::equal(out_std[i], out_trunc[i])) {
            match = false;
            std::cout << "Mismatch at index " << i << std::endl;
            break;
        }
    }

    EXPECT_TRUE(match) << "Truncated middle product does not match standard version";
}
#endif

static Fr fr_zero()
{
    Fr z;
    memset(&z, 0, sizeof(Fr));
    return z;
}

static std::vector<Fr> middle_product_oracle_fr(
    const std::vector<Fr>& a,
    const std::vector<Fr>& b)
{
    const std::size_t n = a.size();
    EXPECT_EQ(b.size(), 2 * n);

    std::vector<Fr> out(n, fr_zero());
    for (std::size_t k = 0; k < n; ++k) {
        Fr acc = fr_zero();
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t j = n + k - i;
            if (j < (2 * n)) {
                Fr term;
                Fr::mul(term, a[i], b[j]);
                Fr::add(acc, acc, term);
            }
        }
        out[k] = acc;
    }
    return out;
}

static std::vector<G1Point> points_from_scalars_g1(const std::vector<u64>& scalars)
{
    std::vector<G1Point> points;
    points.reserve(scalars.size());
    for (u64 s : scalars) {
        points.push_back(g1(s));
    }
    return points;
}

static std::vector<G2Point> points_from_scalars_g2(const std::vector<u64>& scalars)
{
    std::vector<G2Point> points;
    points.reserve(scalars.size());
    for (u64 s : scalars) {
        points.push_back(g2(s));
    }
    return points;
}

static std::vector<Fr> fr_from_scalars(const std::vector<u64>& scalars)
{
    std::vector<Fr> out;
    out.reserve(scalars.size());
    for (u64 s : scalars) {
        out.push_back(mf(s));
    }
    return out;
}

TEST(Fft, MiddleProductBasicMatchesFieldOracle)
{
    const std::vector<std::size_t> test_sizes = { 3, 5, 9 };

    G1Point g1_gen = g1(1);
    G2Point g2_gen = g2(1);
    GT gt_base = pair(g1_gen, g2_gen);

    for (std::size_t n : test_sizes) {
        std::vector<u64> a_scalars;
        std::vector<u64> b_scalars;
        a_scalars.reserve(n);
        b_scalars.reserve(2 * n);

        for (std::size_t i = 0; i < n; ++i) {
            a_scalars.push_back(static_cast<u64>((7 * i + 3) % 97 + 1));
        }
        for (std::size_t i = 0; i < 2 * n; ++i) {
            b_scalars.push_back(static_cast<u64>((11 * i + 5) % 101 + 1));
        }

        const std::vector<G1Point> a_points = points_from_scalars_g1(a_scalars);
        const std::vector<G2Point> b_points = points_from_scalars_g2(b_scalars);
        const std::vector<Fr> a_fr = fr_from_scalars(a_scalars);
        const std::vector<Fr> b_fr = fr_from_scalars(b_scalars);
        const std::vector<Fr> expected_fr = middle_product_oracle_fr(a_fr, b_fr);

        std::vector<GT> out(n, GT::identity());
        middle_product_basic(
            std::span<const G1Point>(a_points),
            std::span<const G2Point>(b_points),
            std::span<GT>(out));

        bool match = true;
        for (std::size_t i = 0; i < n; ++i) {
            GT expected = GTCoeffCtx::mult(gt_base, expected_fr[i]);
            if (!GT::equal(out[i], expected)) {
                match = false;
                break;
            }
        }

        EXPECT_TRUE(match) << "Basic middle product mismatch at n=" << n;
    }
}

TEST(Fft, MiddleProductAdmissibleRoundingBoundaryCases)
{
    const std::size_t fft_small = 512;
    const std::size_t fft_large = 4096;

    const std::vector<std::size_t> inputs = { 1, 2, 15, 16, 17, 63, 64, 65, 255, 256, 257 };
    for (std::size_t x : inputs) {
        const std::size_t rounded = fft_admissible_size(x, fft_small);
        EXPECT_GE(rounded, static_cast<std::size_t>(1));
        EXPECT_LE(rounded, fft_small);
        EXPECT_EQ(rounded % 16, static_cast<std::size_t>(0));
        if (rounded != fft_small) {
            EXPECT_GE(rounded, x);
        }
    }

    // Near-full thresholds: small uses n - n/8, large uses n - n/16.
    EXPECT_EQ(fft_admissible_size(433, fft_small), static_cast<std::size_t>(448));
    EXPECT_EQ(fft_admissible_size(449, fft_small), fft_small);
    EXPECT_EQ(fft_admissible_size(3839, fft_large), static_cast<std::size_t>(3840));
    EXPECT_EQ(fft_admissible_size(3841, fft_large), fft_large);
}

TEST(Fft, MiddleProductMatchesFieldOracle)
{
    const std::vector<std::size_t> test_sizes = { 3, 5, 9, 24 };

    G1Point g1_gen = g1(1);
    G2Point g2_gen = g2(1);
    GT gt_base = pair(g1_gen, g2_gen);

    for (std::size_t n : test_sizes) {
        std::vector<u64> a_scalars;
        std::vector<u64> b_scalars;
        a_scalars.reserve(n);
        b_scalars.reserve(2 * n);

        for (std::size_t i = 0; i < n; ++i) {
            a_scalars.push_back(static_cast<u64>((59 * i + 17) % 233 + 1));
        }
        for (std::size_t i = 0; i < 2 * n; ++i) {
            b_scalars.push_back(static_cast<u64>((61 * i + 19) % 239 + 1));
        }

        const std::vector<G1Point> a_points = points_from_scalars_g1(a_scalars);
        const std::vector<G2Point> b_points = points_from_scalars_g2(b_scalars);
        const std::vector<Fr> a_fr = fr_from_scalars(a_scalars);
        const std::vector<Fr> b_fr = fr_from_scalars(b_scalars);
        const std::vector<Fr> expected_fr = middle_product_oracle_fr(a_fr, b_fr);

        std::vector<GT> out(n, GT::identity());
        middle_product(
            std::span<const G1Point>(a_points),
            std::span<const G2Point>(b_points),
            std::span<GT>(out));

        bool match = true;
        for (std::size_t i = 0; i < n; ++i) {
            GT expected = GTCoeffCtx::mult(gt_base, expected_fr[i]);
            if (!GT::equal(out[i], expected)) {
                match = false;
                break;
            }
        }

        EXPECT_TRUE(match) << "Reduced middle product mismatch at n=" << n;
    }
}

TEST(Fft, MiddleProductMatchesBasic)
{
    const std::vector<std::size_t> test_sizes = { 4, 6, 10, 24 };

    for (std::size_t n : test_sizes) {
        std::vector<u64> a_scalars;
        std::vector<u64> b_scalars;
        a_scalars.reserve(n);
        b_scalars.reserve(2 * n);

        for (std::size_t i = 0; i < n; ++i) {
            a_scalars.push_back(static_cast<u64>((67 * i + 23) % 241 + 1));
        }
        for (std::size_t i = 0; i < 2 * n; ++i) {
            b_scalars.push_back(static_cast<u64>((71 * i + 29) % 251 + 1));
        }

        const std::vector<G1Point> a_points = points_from_scalars_g1(a_scalars);
        const std::vector<G2Point> b_points = points_from_scalars_g2(b_scalars);

        std::vector<GT> expected(n, GT::identity());
        std::vector<GT> reduced_raw(n, GT::identity());
        std::vector<GT> reduced_pre(n, GT::identity());

        middle_product_basic(
            std::span<const G1Point>(a_points),
            std::span<const G2Point>(b_points),
            std::span<GT>(expected));

        middle_product(
            std::span<const G1Point>(a_points),
            std::span<const G2Point>(b_points),
            std::span<GT>(reduced_raw));

        const auto precomputed_b = precompute_middle_product_b(
            std::span<const G2Point>(b_points), n);
        middle_product(
            std::span<const G1Point>(a_points),
            precomputed_b,
            std::span<GT>(reduced_pre));

        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_TRUE(GT::equal(expected[i], reduced_raw[i]))
                << "Reduced raw mismatch vs basic at n=" << n << ", i=" << i;
            EXPECT_TRUE(GT::equal(reduced_raw[i], reduced_pre[i]))
                << "Reduced precompute mismatch at n=" << n << ", i=" << i;
        }
    }
}

TEST(Fft, MiddleProductPrecomputedMatchesRaw)
{
    const std::vector<std::size_t> test_sizes = { 4, 6, 10 };

    for (std::size_t n : test_sizes) {
        std::vector<u64> a_scalars;
        std::vector<u64> b_scalars;
        a_scalars.reserve(n);
        b_scalars.reserve(2 * n);

        for (std::size_t i = 0; i < n; ++i) {
            a_scalars.push_back(static_cast<u64>((19 * i + 1) % 137 + 1));
        }
        for (std::size_t i = 0; i < 2 * n; ++i) {
            b_scalars.push_back(static_cast<u64>((23 * i + 2) % 149 + 1));
        }

        const std::vector<G1Point> a_points = points_from_scalars_g1(a_scalars);
        const std::vector<G2Point> b_points = points_from_scalars_g2(b_scalars);

        std::vector<GT> raw(n, GT::identity());
        std::vector<GT> precomputed(n, GT::identity());

        middle_product(
            std::span<const G1Point>(a_points),
            std::span<const G2Point>(b_points),
            std::span<GT>(raw));

        const auto precomputed_b = precompute_middle_product_b(
            std::span<const G2Point>(b_points), n);
        middle_product(
            std::span<const G1Point>(a_points),
            precomputed_b,
            std::span<GT>(precomputed));

        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_TRUE(GT::equal(raw[i], precomputed[i]))
                << "Reduced precompute mismatch at n=" << n << ", i=" << i;
        }
    }
}

TEST(Fft, MiddleProductWorkspaceMatchesDefault)
{
    const std::vector<std::size_t> test_sizes = { 6, 10, 24 };

    for (std::size_t n : test_sizes) {
        std::vector<u64> a_scalars;
        std::vector<u64> b_scalars;
        a_scalars.reserve(n);
        b_scalars.reserve(2 * n);

        for (std::size_t i = 0; i < n; ++i) {
            a_scalars.push_back(static_cast<u64>((73 * i + 31) % 257 + 1));
        }
        for (std::size_t i = 0; i < 2 * n; ++i) {
            b_scalars.push_back(static_cast<u64>((79 * i + 37) % 263 + 1));
        }

        const std::vector<G1Point> a_points = points_from_scalars_g1(a_scalars);
        const std::vector<G2Point> b_points = points_from_scalars_g2(b_scalars);
        const auto precomputed_b = precompute_middle_product_b(
            std::span<const G2Point>(b_points), n);

        std::vector<GT> out_default(n, GT::identity());
        std::vector<GT> out_ws(n, GT::identity());
        MiddleProductWorkspace workspace;

        middle_product(
            std::span<const G1Point>(a_points),
            precomputed_b,
            std::span<GT>(out_default));

        middle_product(
            std::span<const G1Point>(a_points),
            precomputed_b,
            std::span<GT>(out_ws),
            workspace);

        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_TRUE(GT::equal(out_default[i], out_ws[i]))
                << "Reduced workspace mismatch at n=" << n << ", i=" << i;
        }
    }
}

TEST(Fft, MiddleProductMultipairMatchesFft)
{
    const std::vector<std::size_t> test_sizes = { 6, 10, 24 };

    for (std::size_t n : test_sizes) {
        std::vector<u64> a_scalars;
        std::vector<u64> b_scalars;
        a_scalars.reserve(n);
        b_scalars.reserve(2 * n);

        for (std::size_t i = 0; i < n; ++i) {
            a_scalars.push_back(static_cast<u64>((83 * i + 41) % 269 + 1));
        }
        for (std::size_t i = 0; i < 2 * n; ++i) {
            b_scalars.push_back(static_cast<u64>((89 * i + 43) % 271 + 1));
        }

        const std::vector<G1Point> a_points = points_from_scalars_g1(a_scalars);
        const std::vector<G2Point> b_points = points_from_scalars_g2(b_scalars);
        std::vector<GT> out_fft(n, GT::identity());
        std::vector<GT> out_multipair(n, GT::identity());

        middle_product(
            std::span<const G1Point>(a_points),
            std::span<const G2Point>(b_points),
            std::span<GT>(out_fft),
            MiddleProductAlgorithm::Fft);
        middle_product(
            std::span<const G1Point>(a_points),
            std::span<const G2Point>(b_points),
            std::span<GT>(out_multipair),
            MiddleProductAlgorithm::Multipair);

        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_TRUE(GT::equal(out_fft[i], out_multipair[i]))
                << "Reduced multipair mismatch at n=" << n << ", i=" << i;
        }
    }
}
