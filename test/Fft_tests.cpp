#include <iostream>
#include <vector>
#include "blst/fft/ctFft.h"
#include <gtest/gtest.h>
#include "BasicFp.h"
#include "blst/fields/Fr.h"
#include "blst/curve/point.h"
#include "blst/curve/GT.h"

using F65537 = BasicFp<65537, u32, u64>;
//struct F65537 {
//    long long value;
//    static const long long MOD = 65537; // For testing only
//
//    F65537(long long v = 0) : value((v% MOD + MOD) % MOD) {}
//
//    F65537 operator+(const F65537& other) const { return F65537(value + other.value); }
//    F65537 operator-(const F65537& other) const { return F65537(value - other.value); }
//    F65537 operator*(const F65537& other) const { return F65537(value * other.value); }
//
//    F65537& operator=(const F65537& other) { value = other.value; return *this; }
//
//    // Correct modular inverse for the iFFT scaling factor [cite: 154]
//    F65537 inv() const {
//        return pow(MOD - 2);
//    }
//
//    F65537 pow(long long exp) const {
//        F65537 res(1);
//        auto base = *this;
//
//        while (exp > 0) {
//            if (exp % 2 == 1) res = res * base;
//            base = base * base;
//            exp /= 2;
//        }
//        return res;
//    }
//};

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


TEST(Fft, basic) {

    auto n = 4;
    size_t N = n*2; // FFT size (2*B) to avoid aliasing [cite: 150]

    std::vector<F65537> A = { 1, 2, 3, 4 }; // Padded keys [cite: 149]
    std::vector<F65537> B = { 5, 6, 7, 8 }; // Padded CRS [cite: 149]

    // 1. Naive Convolution (The "B^2" bottleneck) [cite: 36, 122]
    std::vector<F65537> expected = naive_multiply<F65537, F65537, F65537>(A, B);
	expected.resize(N); // Pad to N

    auto prim = primRootOfUnity<F65537>(N);
    std::vector<F65537> roots(N);
    precomputeRootsOfUnity< F65537>(prim, roots);
    
    std::vector<F65537> inv_roots(N);
    for (size_t i = 0; i < N; ++i) {
        inv_roots[i] = roots[i].inverse();
    }
    // 2. FFT Convolution 
    // (Assume roots are precomputed for simplicity in this test)

    A.resize(N);
    B.resize(N);

    fft(A, roots, false);
    fft(B, roots, false);

    std::vector<F65537> C(N);
    for (int i = 0; i < N; ++i) C[i] = A[i] * B[i]; // Pointwise pairing 

    fft(C, inv_roots, true); // Inverse FFT 

    std::cout << "FFT Convolution Result (first 4): ";
    for (int i = 0; i < 4; ++i) std::cout << C[i] << " ";
    std::cout << "\nNaive Result (first 4): ";
    for (int i = 0; i < 4; ++i) std::cout << expected[i] << " ";

    EXPECT_EQ(C, expected);
}


// Explicit specialization of primRootOfUnity for blst::Fr.
// This bypasses the random search in findGenerator by using the known
// structure of the BLS12-381 scalar field (2-adicity = 32, generator = 7).
template<>
blst::Fr primRootOfUnity<blst::Fr>(u64 n)
{
    // The explicit table cache is initialized exactly once (thread-safe).
    static const std::array<blst::Fr, 33> roots = []() {
        std::array<blst::Fr, 33> table; // Stores primitive 2^k-th roots for k=0..32

		table[32][0] = 0xb9b58d8c5f0e466a;
		table[32][1] = 0x5b1b4c801819d7ec;
		table[32][2] = 0x0af53ae352a31e64;
		table[32][3] = 0x5bf3adda19e9b27b;

        // Compute 2^k roots by repeated squaring
        for (int i = 31; i >= 0; --i) {
            blst::Fr::sqr(table[i], table[i + 1]);
        }
        return table;

        }();

    // Verify n is a power of 2
    u64 k = 0;
    if (n > 0 && (n & (n - 1)) == 0) {
        u64 v = n;
        while (v >>= 1) k++;
    }
    else {
        throw std::runtime_error("primRootOfUnity: n must be a power of 2");
    }

    if (k >= roots.size()) {
        throw std::runtime_error("primRootOfUnity: n exceeds 2-adicity (2^32)");
    }

    // roots[k] is the primitive 2^k-th root of unity
    return roots[k];
}

template<>
inline bool isPrimRootOfUnity<Fr>(std::span<Factor> factors, const Fr& u)
{
    auto p = Fr::modulus();
    u64 n = 1;
    for (auto fe : factors)
    {
        for (u64 i = 0; i < fe.mExp; ++i)
            n *= fe.mFactor;
    }

    //Fr pMinusOne, N, pMinusOneModN;
    //Fr::sub(pMinusOne, p, Fr::one());
    //Fr::from_uint64(N, n);
    //Fr::mod(pMinusOneModN, pMinusOneMod, N);

    //if (pMinusOne % N)
    //    throw RTE_LOC;

    // make suer u is in Fp*
    if (u == Fr::zero())
        return false;

    Fr uPowN;
    Fr::exp(uPowN, u, (byte*)&n, std::bit_width(n));
    // make sure u is a root of unity.
    if (uPowN != Fr::one())
        return false;

    // check that u is a primitive root of unity.
    for (u64 i = 1; i < factors.size() - 1; ++i)
    {
        auto ni = n / factors[i].mFactor;
        Fr uPowNi;
        Fr::exp(uPowNi, u, (byte*)&ni, std::bit_width(ni));
        if (uPowNi == Fr::one())
        {
            return false;
        }
    }
    return true;
}

// returns true if f is a primitive root of unity.
template<>
inline bool isPrimRootOfUnity<Fr>(u64 n, const Fr& u)
{
    auto factors = uniqueFactor(n);
    return isPrimRootOfUnity<Fr>(factors, u);
}

TEST(Fft, BlstFrRoots) {
    auto one = blst::Fr::one();
    blst::Fr negOne;
    blst::Fr::neg(negOne, one);

    EXPECT_TRUE(blst::Fr::equal(primRootOfUnity<blst::Fr>(1), blst::Fr::one()));

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


struct FrCoeffCtx
{
    static auto plus(const Fr& a, const Fr& b) {
        Fr r;
        Fr::add(r, a, b);
        return r;
    }
    static auto minus(const Fr& a, const Fr& b) {
        Fr r;
        Fr::sub(r, a, b);
        return r;
    }

    static auto mult(const Fr& a, const Fr& b) {
        Fr r;
        Fr::mul(r, a, b);
        return r;
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

// Context for FFT operations on G1 points with Fr scalars (FFT in the exponent)
struct G1CoeffCtx
{
    // G1 + G1 -> G1 (point addition)
    static auto plus(const G1Point& a, const G1Point& b) {
        G1Point r;
        G1Point::dadd(r, a, b);
        return r;
    }

    // G1 - G1 -> G1 (point subtraction)
    static auto minus(const G1Point& a, const G1Point& b) {
        G1Point neg_b;
        G1Point::neg(neg_b, b);
        G1Point r;
        G1Point::dadd(r, a, neg_b);
        return r;
    }

    // G1 * Fr -> G1 (scalar multiplication)
    static auto mult(const G1Point& a, const Fr& b) {
        G1Point r;
        // Convert Fr from Montgomery form to scalar bytes for mult
        Fr b_non_mont;
        Fr::from_mont(b_non_mont, b);
        G1Point::mult(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
        return r;
    }

    // Inverse of Fr scalar (for iFFT scaling)
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

// Context for FFT operations on G2 points with Fr scalars (FFT in the exponent)
struct G2CoeffCtx
{
    // G2 + G2 -> G2 (point addition)
    static auto plus(const G2Point& a, const G2Point& b) {
        G2Point r;
        G2Point::dadd(r, a, b);
        return r;
    }

    // G2 - G2 -> G2 (point subtraction)
    static auto minus(const G2Point& a, const G2Point& b) {
        G2Point neg_b;
        G2Point::neg(neg_b, b);
        G2Point r;
        G2Point::dadd(r, a, neg_b);
        return r;
    }

    // G2 * Fr -> G2 (scalar multiplication)
    static auto mult(const G2Point& a, const Fr& b) {
        G2Point r;
        // Convert Fr from Montgomery form to scalar bytes for mult
        Fr b_non_mont;
        Fr::from_mont(b_non_mont, b);
        G2Point::mult(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
        return r;
    }

    // Inverse of Fr scalar (for iFFT scaling)
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




// Context for FFT operations on GT elements with Fr scalars
struct GTCoeffCtx
{
    // GT + GT -> GT (group addition = Fp12 multiplication)
    static auto plus(const GT& a, const GT& b) {
        GT r;
        GT::add(r, a, b);
        return r;
    }

    // GT - GT -> GT (group subtraction = Fp12 division)
    static auto minus(const GT& a, const GT& b) {
        GT r;
        GT::sub(r, a, b);
        return r;
    }

    // GT * Fr -> GT (scalar multiplication = Fp12 exponentiation)
    static auto mult(const GT& a, const Fr& b) {
        GT r;
        Fr b_non_mont;
        Fr::from_mont(b_non_mont, b);
        GT::mult(r, a, reinterpret_cast<const byte*>(&b_non_mont), 256);
        return r;
    }

    // Inverse of Fr scalar (for iFFT scaling)
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
    size_t N = 8;

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
    auto n = 4;
    size_t N = n * 2; // FFT size

    // 1. Setup data: G1 points (A) and G2 points (B)
    std::vector<G1Point> A = { g1(1), g1(2), g1(3), g1(4) };
    std::vector<G2Point> B = { g2(5), g2(6), g2(7), g2(8) };

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