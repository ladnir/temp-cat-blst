#include <iostream>
#include <vector>
#include <gtest/gtest.h>

#include "blst/fft/ctFft.h"
#include "blst/fft/ctFft_cuda.h"
#include "blst/fields/Fr.h"
#include "blst/fft/primitveRoot.h"
#include "blst/fft/context.h"
#include "blst/curve/point.h"
#include "blst/curve/GT.h"

#include "blst/common/defines.h"
#include <array>
#include <vector>
#include <random>
#include <span>
#include <exception>
using namespace blst;

namespace
{

	// Helper to create Fr from u64
	Fr make_fr(u64 i)
	{
		Fr r;
		Fr::from_uint64(r, i);
		return r;
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

	// Pairing function for GT tests
	GT pair_points(const G1Point& p, const G2Point& q) {
		G1AffinePoint p_aff;
		G1Point::to_affine(p_aff, p);
		G2AffinePoint q_aff;
		G2Point::to_affine(q_aff, q);

		GT result;
		GT::miller_loop(result, &q_aff, &p_aff);
		GT::final_mult(result, result);
		return result;
	}
}

#ifdef BLST_CUDA

TEST(FftCuda, BasicRoundTrip) {
	size_t N = 8;

	// Create test data
	std::vector<Fr> original(N);
	for (size_t i = 0; i < N; ++i) {
		original[i] = Fr::from_uint64(i + 1);
	}

	// Copy for CUDA test
	std::vector<Fr> data = original;

	// Get primitive root and precompute roots
	Fr prim = primRootOfUnity<Fr>(N);
	std::vector<Fr> roots(N);
	Fr current = Fr::one();
	for (size_t i = 0; i < N; ++i) {
		roots[i] = current;
		Fr::mul(current, current, prim);
	}

	std::vector<Fr> inv_roots(N);
	for (size_t i = 0; i < N; ++i) {
		Fr::reciprocal(inv_roots[i], roots[i]);
	}

	// Forward FFT on GPU
	fft_cuda(data, roots, false);

	// Inverse FFT on GPU  
	fft_cuda(data, inv_roots, true);

	// Compare with original
	bool match = true;
	for (size_t i = 0; i < N; ++i) {
		if (!Fr::equal(data[i], original[i])) {
			match = false;
			std::cout << "Mismatch at index " << i << std::endl;
			break;
		}
	}

	EXPECT_TRUE(match) << "CUDA FFT round-trip failed";
}

TEST(FftCuda, CompareWithCpu) {
	size_t N = 16;

	// Create test data
	std::vector<Fr> cpu_data(N);
	for (size_t i = 0; i < N; ++i) {
		cpu_data[i] = Fr::from_uint64((i + 1) * 7);
	}
	std::vector<Fr> gpu_data = cpu_data;

	// Get roots
	Fr prim = primRootOfUnity<Fr>(N);
	std::vector<Fr> roots(N);
	Fr current = Fr::one();
	for (size_t i = 0; i < N; ++i) {
		roots[i] = current;
		Fr::mul(current, current, prim);
	}

	fft<Fr, Fr, FrCoeffCtx>(cpu_data, roots, false);

	// GPU FFT
	fft_cuda(gpu_data, roots, false);

	// Compare results
	bool match = true;
	for (size_t i = 0; i < N; ++i) {
		if (!Fr::equal(cpu_data[i], gpu_data[i])) {
			match = false;
			char cpu_buf[65], gpu_buf[65];
			Fr::to_hexascii(cpu_buf, cpu_data[i]);
			Fr::to_hexascii(gpu_buf, gpu_data[i]);
			std::cout << "Mismatch at index " << i << std::endl;
			std::cout << "  CPU: " << cpu_buf << std::endl;
			std::cout << "  GPU: " << gpu_buf << std::endl;
			break;
		}
	}

	EXPECT_TRUE(match) << "CUDA FFT does not match CPU FFT";
}

TEST(FftCuda, Convolution) {
	size_t n = 4;
	size_t N = n * 2;

	// Input polynomials
	std::vector<Fr> A = { Fr::from_uint64(1), Fr::from_uint64(2), Fr::from_uint64(3), Fr::from_uint64(4) };
	std::vector<Fr> B = { Fr::from_uint64(5), Fr::from_uint64(6), Fr::from_uint64(7), Fr::from_uint64(8) };

	// Compute expected result via naive convolution
	std::vector<Fr> expected(N);
	Fr zero;
	memset(&zero, 0, sizeof(Fr));
	for (size_t i = 0; i < N; ++i) expected[i] = zero;

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			Fr prod;
			Fr::mul(prod, A[i], B[j]);
			Fr::add(expected[i + j], expected[i + j], prod);
		}
	}

	// Pad inputs
	A.resize(N, zero);
	B.resize(N, zero);

	// Get roots
	Fr prim = primRootOfUnity<Fr>(N);
	std::vector<Fr> roots(N);
	Fr current = Fr::one();
	for (size_t i = 0; i < N; ++i) {
		roots[i] = current;
		Fr::mul(current, current, prim);
	}

	std::vector<Fr> inv_roots(N);
	for (size_t i = 0; i < N; ++i) {
		Fr::reciprocal(inv_roots[i], roots[i]);
	}

	// FFT convolution on GPU
	fft_cuda(A, roots, false);
	fft_cuda(B, roots, false);

	// Pointwise multiplication
	std::vector<Fr> C(N);
	for (size_t i = 0; i < N; ++i) {
		Fr::mul(C[i], A[i], B[i]);
	}

	// Inverse FFT
	fft_cuda(C, inv_roots, true);

	// Compare
	bool match = true;
	for (size_t i = 0; i < N; ++i) {
		if (!Fr::equal(C[i], expected[i])) {
			match = false;
			break;
		}
	}

	EXPECT_TRUE(match) << "CUDA FFT convolution failed";
}

// ==================== G1 CUDA FFT Tests ====================

TEST(FftCuda, G1RoundTrip) {
	size_t N = 8;

	// Create test G1 points
	std::vector<G1Point> original(N);
	for (size_t i = 0; i < N; ++i) {
		original[i] = g1(i + 1);
	}
	std::vector<G1Point> data = original;

	// Get roots
	Fr prim = primRootOfUnity<Fr>(N);
	std::vector<Fr> roots(N);
	Fr current = Fr::one();
	for (size_t i = 0; i < N; ++i) {
		roots[i] = current;
		Fr::mul(current, current, prim);
	}

	std::vector<Fr> inv_roots(N);
	for (size_t i = 0; i < N; ++i) {
		Fr::reciprocal(inv_roots[i], roots[i]);
	}

	// Forward FFT on GPU
	fft_cuda(data, roots, false);

	// Inverse FFT on GPU
	fft_cuda(data, inv_roots, true);

	// Compare with original
	bool match = true;
	for (size_t i = 0; i < N; ++i) {
		if (!G1Point::equal(data[i], original[i])) {
			match = false;
			std::cout << "G1 Mismatch at index " << i << std::endl;
			break;
		}
	}

	EXPECT_TRUE(match) << "CUDA G1 FFT round-trip failed";
}

TEST(FftCuda, G1CompareWithCpu) {
	size_t N = 8;

	// Create test data
	std::vector<G1Point> cpu_data(N);
	for (size_t i = 0; i < N; ++i) {
		cpu_data[i] = g1((i + 1) * 3);
	}
	std::vector<G1Point> gpu_data = cpu_data;

	// Get roots
	Fr prim = primRootOfUnity<Fr>(N);
	std::vector<Fr> roots(N);
	Fr current = Fr::one();
	for (size_t i = 0; i < N; ++i) {
		roots[i] = current;
		Fr::mul(current, current, prim);
	}

	// CPU FFT
	fft<G1Point, Fr, G1CoeffCtx>(cpu_data, roots, false);

	// GPU FFT
	fft_cuda(gpu_data, roots, false);

	// Compare results
	bool match = true;
	for (size_t i = 0; i < N; ++i) {
		if (!G1Point::equal(cpu_data[i], gpu_data[i])) {
			match = false;
			std::cout << "G1 Mismatch at index " << i << std::endl;
			break;
		}
	}

	EXPECT_TRUE(match) << "CUDA G1 FFT does not match CPU FFT";
}

// ==================== G2 CUDA FFT Tests ====================

TEST(FftCuda, G2RoundTrip) {
	size_t N = 8;

	// Create test G2 points
	std::vector<G2Point> original(N);
	for (size_t i = 0; i < N; ++i) {
		original[i] = g2(i + 1);
	}
	std::vector<G2Point> data = original;

	// Get roots
	Fr prim = primRootOfUnity<Fr>(N);
	std::vector<Fr> roots(N);
	Fr current = Fr::one();
	for (size_t i = 0; i < N; ++i) {
		roots[i] = current;
		Fr::mul(current, current, prim);
	}

	std::vector<Fr> inv_roots(N);
	for (size_t i = 0; i < N; ++i) {
		Fr::reciprocal(inv_roots[i], roots[i]);
	}

	// Forward FFT on GPU
	fft_cuda(data, roots, false);

	// Inverse FFT on GPU
	fft_cuda(data, inv_roots, true);

	// Compare with original
	bool match = true;
	for (size_t i = 0; i < N; ++i) {
		if (!G2Point::equal(data[i], original[i])) {
			match = false;
			std::cout << "G2 Mismatch at index " << i << std::endl;
			break;
		}
	}

	EXPECT_TRUE(match) << "CUDA G2 FFT round-trip failed";
}

TEST(FftCuda, G2CompareWithCpu) {
	size_t N = 8;

	// Create test data
	std::vector<G2Point> cpu_data(N);
	for (size_t i = 0; i < N; ++i) {
		cpu_data[i] = g2((i + 1) * 5);
	}
	std::vector<G2Point> gpu_data = cpu_data;

	// Get roots
	Fr prim = primRootOfUnity<Fr>(N);
	std::vector<Fr> roots(N);
	Fr current = Fr::one();
	for (size_t i = 0; i < N; ++i) {
		roots[i] = current;
		Fr::mul(current, current, prim);
	}

	// CPU FFT
	fft<G2Point, Fr, G2CoeffCtx>(cpu_data, roots, false);

	// GPU FFT
	fft_cuda(gpu_data, roots, false);

	// Compare results
	bool match = true;
	for (size_t i = 0; i < N; ++i) {
		if (!G2Point::equal(cpu_data[i], gpu_data[i])) {
			match = false;
			std::cout << "G2 Mismatch at index " << i << std::endl;
			break;
		}
	}

	EXPECT_TRUE(match) << "CUDA G2 FFT does not match CPU FFT";
}

// ==================== GT CUDA FFT Tests ====================

TEST(FftCuda, GTRoundTrip) {
	size_t N = 8;

	// Create test GT elements via pairings
	G1AffinePoint P_aff;
	G2AffinePoint Q_aff;
	G1Point::to_affine(P_aff, G1Point::generator());
	G2Point::to_affine(Q_aff, G2Point::generator());

	GT base;
	GT::miller_loop(base, &Q_aff, &P_aff);
	GT::final_mult(base, base);

	std::vector<GT> original(N);
	original[0] = GT::identity();
	GT accum = base;
	for (size_t i = 1; i < N; ++i) {
		original[i] = accum;
		GT::add(accum, accum, base);
	}

	std::vector<GT> data = original;

	// Get roots
	Fr prim = primRootOfUnity<Fr>(N);
	std::vector<Fr> roots(N);
	Fr current = Fr::one();
	for (size_t i = 0; i < N; ++i) {
		roots[i] = current;
		Fr::mul(current, current, prim);
	}

	std::vector<Fr> inv_roots(N);
	for (size_t i = 0; i < N; ++i) {
		Fr::reciprocal(inv_roots[i], roots[i]);
	}

	// Forward FFT on GPU
	fft_cuda(data, roots, false);

	// Inverse FFT on GPU
	fft_cuda(data, inv_roots, true);

	// Compare with original
	bool match = true;
	for (size_t i = 0; i < N; ++i) {
		if (!GT::equal(data[i], original[i])) {
			match = false;
			std::cout << "GT Mismatch at index " << i << std::endl;
			break;
		}
	}

	EXPECT_TRUE(match) << "CUDA GT FFT round-trip failed";
}

TEST(FftCuda, GTCompareWithCpu) {
	size_t N = 8;

	// Create test GT elements
	G1AffinePoint P_aff;
	G2AffinePoint Q_aff;
	G1Point::to_affine(P_aff, G1Point::generator());
	G2Point::to_affine(Q_aff, G2Point::generator());

	GT base;
	GT::miller_loop(base, &Q_aff, &P_aff);
	GT::final_mult(base, base);

	std::vector<GT> cpu_data(N);
	cpu_data[0] = GT::identity();
	GT accum = base;
	for (size_t i = 1; i < N; ++i) {
		cpu_data[i] = accum;
		GT::add(accum, accum, base);
	}

	std::vector<GT> gpu_data = cpu_data;

	// Get roots
	Fr prim = primRootOfUnity<Fr>(N);
	std::vector<Fr> roots(N);
	Fr current = Fr::one();
	for (size_t i = 0; i < N; ++i) {
		roots[i] = current;
		Fr::mul(current, current, prim);
	}

	// CPU FFT
	fft<GT, Fr, GTCoeffCtx>(cpu_data, roots, false);

	// GPU FFT
	fft_cuda(gpu_data, roots, false);

	// Compare results
	bool match = true;
	for (size_t i = 0; i < N; ++i) {
		if (!GT::equal(cpu_data[i], gpu_data[i])) {
			match = false;
			std::cout << "GT Mismatch at index " << i << std::endl;
			break;
		}
	}

	EXPECT_TRUE(match) << "CUDA GT FFT does not match CPU FFT";
}

// ==================== Mixed Convolution Test (G1 * Fr) ====================

TEST(FftCuda, G1FrConvolution) {
	size_t n = 4;
	size_t N = n * 2;

	// G1 points and Fr scalars
	std::vector<G1Point> A = { g1(1), g1(2), g1(3), g1(4) };
	std::vector<Fr> B = { make_fr(5), make_fr(6), make_fr(7), make_fr(8) };

	// Compute expected via naive convolution
	G1Point identity = G1Point::infinity();
	std::vector<G1Point> expected(N, identity);

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			G1Point prod = G1CoeffCtx::mult(A[i], B[j]);
			G1Point::dadd(expected[i + j], expected[i + j], prod);
		}
	}

	// Pad inputs
	Fr zero;
	memset(&zero, 0, sizeof(Fr));
	A.resize(N, identity);
	B.resize(N, zero);

	// Get roots
	Fr prim = primRootOfUnity<Fr>(N);
	std::vector<Fr> roots(N);
	Fr current = Fr::one();
	for (size_t i = 0; i < N; ++i) {
		roots[i] = current;
		Fr::mul(current, current, prim);
	}

	std::vector<Fr> inv_roots(N);
	for (size_t i = 0; i < N; ++i) {
		Fr::reciprocal(inv_roots[i], roots[i]);
	}

	// FFT on G1 (GPU) and Fr (GPU)
	fft_cuda(A, roots, false);
	fft_cuda(B, roots, false);

	// Pointwise multiplication (on host for simplicity)
	std::vector<G1Point> C(N);
	for (size_t i = 0; i < N; ++i) {
		C[i] = G1CoeffCtx::mult(A[i], B[i]);
	}

	// Inverse FFT on G1
	fft_cuda(C, inv_roots, true);

	// Compare
	bool match = true;
	for (size_t i = 0; i < N; ++i) {
		if (!G1Point::equal(C[i], expected[i])) {
			match = false;
			std::cout << "G1*Fr Convolution mismatch at index " << i << std::endl;
			break;
		}
	}

	EXPECT_TRUE(match) << "CUDA G1*Fr FFT convolution failed";
}

#endif // BLST_CUDA