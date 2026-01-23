
#pragma once
#include <vector>
#include <span>
#include "blst/common/defines.h"

// Standard Bit-Reversal for in-place FFT
template<typename F>
void bit_reverse(std::span<F> a) {
	int n = a.size();
	for (int i = 1, j = 0; i < n; i++) {
		int bit = n >> 1;
		for (; j & bit; bit >>= 1) j ^= bit;
		j ^= bit;
		if (i < j) std::swap(a[i], a[j]);
	}
}

struct DefaultCoeffCtx
{
	template<typename G, typename F>
	static auto plus(const G& a, const F& b) {
		return a + b;
	}
	template<typename G, typename F>
	static auto minus(const G& a, const F& b) {
		return a - b;
	}
	template<typename G, typename F>
	static auto mult(const G& a, const F& b) {
		return a * b;
	}
	template<typename F>
	static auto inverse(const F& a) {
		return a.inverse();
	}

	template<typename F>
	static auto from_u64(std::uint64_t a) {
		return F(a);
	}

};

template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft(std::span<G> a, std::span<const F> roots, bool invert) {
	int n = a.size();
	bit_reverse(a);

	for (int len = 2; len <= n; len <<= 1) {
		int half = len / 2;
		for (int i = 0; i < n; i += len) {
			for (int j = 0; j < half; j++) {
				// In a real implementation, you'd precompute or select roots 
				// based on 'len' and 'invert'
				F w = roots[j * (n / len)];


				G u = a[i + j];
				G t = CoeffCtx::mult(a[i + j + half], w);
				a[i + j] = CoeffCtx::plus(u, t);
				a[i + j + half] = CoeffCtx::minus(u, t);
			}
		}
	}

	// Integrated scaling for inverse FFT 
	if (invert) {
		F inv_n = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(n));
		for (G& x : a) {
			x = CoeffCtx::mult(x, inv_n); // Scale every group element by 1/n
		}
	}
}

// Computes C = A * B via FFT-based polynomial multiplication
// All vectors must be pre-sized to n (the FFT size, must be power of 2)
// Caller is responsible for padding inputs and extracting desired output coefficients
//
// Template parameters:
//   G1 - coefficient type for polynomial A
//   G2 - coefficient type for polynomial B  
//   GT - coefficient type for product polynomial C
//   F  - type for FFT roots of unity
//   CoeffCtx - context providing arithmetic operations
template <typename G1, typename G2, typename GT, typename F, typename CoeffCtx = DefaultCoeffCtx>
void poly_mult_fft(
	std::span<G1> A,
	std::span<G2> B,
	bool skipBFft,
	std::span<GT> C,
	std::span<const F> roots,
	std::span<const F> inv_roots)
{
	std::size_t n = A.size();

	// Forward FFT on both polynomials (in-place)
	fft<G1, F, CoeffCtx>(A, roots, false);

	if (skipBFft == false)
		fft<G2, F, CoeffCtx>(B, roots, false);

	// Pointwise multiplication in evaluation domain
	for (std::size_t i = 0; i < n; ++i) {
		C[i] = CoeffCtx::mult(A[i], B[i]);
	}

	// Inverse FFT to get product coefficients
	fft<GT, F, CoeffCtx>(C, inv_roots, true);
}

//// Computes the middle n coefficients of C = A * B
//// A has n coefficients (type G1), B has 2n coefficients (type G2)
//// Returns coefficients c[n-1], c[n], ..., c[2n-2] (0-indexed)
//// which are the 1-indexed coefficients n, n+1, ..., 2n-1
////
//// Template parameters:
////   G1 - coefficient type for polynomial A
////   G2 - coefficient type for polynomial B  
////   GT - coefficient type for product polynomial C
////   F  - type for FFT roots of unity
////   CoeffCtx - context providing arithmetic operations
//template <typename G1, typename G2, typename GT, typename F, typename CoeffCtx = DefaultCoeffCtx>
//std::vector<GT> poly_mult_middle(
//	const std::vector<G1>& A,
//	const std::vector<G2>& B,
//	bool skipFftB,
//	const std::vector<F>& roots,
//	const std::vector<F>& inv_roots)
//{
//	std::size_t n = A.size();
//	// B should have size 2n, product has size 3n-1
//	// We need FFT size >= 3n-1, use next power of 2
//
//	std::size_t fft_size = 1;
//	while (fft_size < 3 * n - 1) {
//		fft_size <<= 1;
//	}
//
//	// Pad A to fft_size (A has n elements)
//	std::vector<G1> a_padded(fft_size);
//	for (std::size_t i = 0; i < n; ++i) {
//		a_padded[i] = CoeffCtx::template from_g1<GT>(A[i]);
//	}
//	// Remaining elements are zero-initialized
//
//	// Pad B to fft_size (B has 2n elements)
//	std::vector<GT> b_padded(fft_size);
//	for (std::size_t i = 0; i < B.size(); ++i) {
//		b_padded[i] = CoeffCtx::template from_g2<GT>(B[i]);
//	}
//
//	// Prepare output vector for product coefficients
//	std::vector<GT> c_eval(fft_size);
//	// Perform polynomial multiplication via FFT
//	poly_mult_fft<G1, G2, GT, F, CoeffCtx>(a_padded, b_padded, c_eval, roots, inv_roots);
//
//	// Extract middle n coefficients: indices n-1 to 2n-2 (0-indexed)
//	// These correspond to 1-indexed coefficients n, n+1, ..., 2n-1
//	std::vector<GT> result(n);
//	for (std::size_t i = 0; i < n; ++i) {
//		result[i] = c_eval[n - 1 + i];
//	}
//
//	return result;
//}





// Derive NTT of reversed-and-padded vector from NTT of original padded vector
// Given A_hat = NTT(A_pad) where A_pad = [a[0], a[1], ..., a[n-1], 0, 0, ..., 0] (length m)
// Computes A_rev_hat = NTT(A_rev_pad) where A_rev_pad = [a[n-1], a[n-2], ..., a[0], 0, ..., 0]
//
// Identity: A_rev_hat[k] = ω^((n-1)*k) * A_hat[(-k) mod m]
//
// Parameters:
//   A_hat     - input: NTT of padded A, size m (in natural frequency order)
//   A_rev_hat - output: NTT of reversed padded A, size m
//   n         - original polynomial size (before padding)
//   m         - FFT size (power of 2)
//   omega     - primitive m-th root of unity
//
// Cost: m scalar multiplications in G, O(m) memory reads/writes
// Saves: one full forward NTT over G of length m
template <typename G, typename F, typename CoeffCtx>
void derive_ntt_of_reversed(
	std::span<const G> A_hat,
	std::span<G> A_rev_hat,
	std::size_t n,
	std::size_t m,
	const F& omega,
	const G& identity_g)
{
	BLST_ASSERT(A_hat.size() >= m);
	BLST_ASSERT(A_rev_hat.size() >= m);
	BLST_ASSERT(n <= m);
	BLST_ASSERT((m & (m - 1)) == 0);

	// Compute γ = ω^(n-1)
	// We need to compute this via repeated squaring or use precomputed roots
	F gamma = omega;
	{
		// Compute ω^(n-1) by repeated squaring
		std::size_t exp = n - 1;
		F base = omega;
		gamma = CoeffCtx::template from_u64<F>(1);
		while (exp > 0) {
			if (exp & 1) {
				gamma = CoeffCtx::mult(gamma, base);
			}
			base = CoeffCtx::mult(base, base);
			exp >>= 1;
		}
	}

	// phase accumulator, starts at ω^0 = 1
	F phase = CoeffCtx::template from_u64<F>(1);

	for (std::size_t k = 0; k < m; ++k) {
		// src = (-k) mod m = (m - k) mod m
		std::size_t src = (k == 0) ? 0 : (m - k);

		// A_rev_hat[k] = phase * A_hat[src] where phase = ω^((n-1)*k)
		A_rev_hat[k] = CoeffCtx::mult(A_hat[src], phase);

		// Update phase: phase *= γ
		phase = CoeffCtx::mult(phase, gamma);
	}
}


// Overload that computes in-place (uses A_hat as both input and output)
// Note: This requires a temporary buffer or careful ordering since we read from
// indices that may be overwritten. For safety, allocates a temp copy.
template <typename G, typename F, typename CoeffCtx>
void derive_ntt_of_reversed_inplace(
	std::span<G> A_hat,
	std::size_t n,
	std::size_t m,
	const F& omega)
{
	// Need temp storage since src and dst indices overlap
	std::vector<G> temp(m);
	derive_ntt_of_reversed<G, F, CoeffCtx>(
		const_span<G>(A_hat.data(), m),
		std::span<G>(temp),
		n, m, omega, A_hat[0] /* just need type, not used */);

	for (std::size_t i = 0; i < m; ++i) {
		A_hat[i] = temp[i];
	}
}


// Computes the middle n coefficients of C = A * B using the reversal trick
// A has n coefficients (type G1), B has 2n coefficients split as B0 and B1
// Writes coefficients c[n], c[n+1], ..., c[2n-1] to out (0-indexed)
//
// OPTIMIZATION: Uses frequency reversal + phase trick to eliminate one G1 NTT
// Instead of computing NTT(rev(A)) directly, we derive it from NTT(A):
//   NTT(rev(A))[k] = ω^((n-1)*k) * NTT(A)[(-k) mod m]
//
// This saves one full G1 forward NTT (m log m scalar mults in G1)
// and replaces it with m scalar mults in G1 (linear cost).
//
// IMPORTANT: This function modifies A, B0, B1 in-place. Caller must copy if needed.
//
// Parameters:
//   A         - input polynomial, size m, input in [0,n), rest used as scratch, MODIFIED
//   B0        - input polynomial (low half of B), size m, input in [0,n), MODIFIED
//   B1        - input polynomial (high half of B), size m, input in [0,n), MODIFIED
//   n         - logical size of A, B0, B1 each have n valid coefficients
//   m         - FFT size, must be power of 2 and >= 2n
//   out       - output buffer, size >= n
//   C_scratch - scratch buffer for FFT results, size m
//   roots     - precomputed roots of unity, size m
//   inv_roots - precomputed inverse roots, size m
//   omega     - primitive m-th root of unity (roots[1])
//   identity_gt - additive identity for GT
//   identity_g1 - additive identity for G1
//   identity_g2 - additive identity for G2
template <typename G1, typename G2, typename GT, typename F, typename CoeffCtx>
void middle_product_fft(
	std::span<G1> A,
	std::span<G2> B0,
	std::span<G2> B1,
	bool skipFftB,
	std::size_t n,
	std::size_t m,
	std::span<GT> out,
	std::span<GT> C_scratch,
	std::span<const F> roots,
	std::span<const F> inv_roots,
	const F& omega,
	const GT& identity_gt,
	const G1& identity_g1,
	const G2& identity_g2)
{
	BLST_ASSERT(A.size() >= m);
	BLST_ASSERT(B0.size() >= m);
	BLST_ASSERT(B1.size() >= m);
	BLST_ASSERT(out.size() >= n);
	BLST_ASSERT(C_scratch.size() >= m);
	BLST_ASSERT(roots.size() >= m);
	BLST_ASSERT(inv_roots.size() >= m);
	BLST_ASSERT(m >= 2 * n);
	BLST_ASSERT((m & (m - 1)) == 0);

	// === Phase 1: Compute p = low_n(A * B1) ===

	// Pad A and B1 with identity
	for (std::size_t i = n; i < m; ++i) {
		A[i] = identity_g1;
		B1[i] = identity_g2;
	}

	// Forward FFT on A
	fft<G1, F, CoeffCtx>(A, roots, false);

	if (skipFftB == false)
	{
		// Forward FFT on B1
		fft<G2, F, CoeffCtx>(B1, roots, false);
	}

	// Pointwise multiply -> C_scratch
	for (std::size_t t = 0; t < m; ++t) {
		C_scratch[t] = CoeffCtx::mult(A[t], B1[t]);
	}

	// Inverse FFT to get (A * B1) coefficients
	fft<GT, F, CoeffCtx>(C_scratch, inv_roots, true);

	// Save p[k] = C_scratch[k] for k = 0..n-1 into out
	for (std::size_t k = 0; k < n; ++k) {
		out[k] = C_scratch[k];
	}

	// === Phase 2: Compute q = high_n(A * B0) via reversal ===

	std::vector<G1> ARev(m);

	// Derive NTT(rev(A)) from NTT(A) using the phase trick
	// A_rev_hat[k] = ω^((n-1)*k) * A_hat[(-k) mod m]
	derive_ntt_of_reversed<G1, F, CoeffCtx>(
		std::span<const G1>(A.data(), m),
		ARev,
		n, m, omega, identity_g1);

	if (skipFftB == false)
	{
		// Reverse B0 in-place and pad
		for (std::size_t i = 0; i < n / 2; ++i) {
			std::swap(B0[i], B0[n - 1 - i]);
		}
		for (std::size_t i = n; i < m; ++i) {
			B0[i] = identity_g2;
		}

		// Forward FFT on rev(B0)
		fft<G2, F, CoeffCtx>(B0, roots, false);
	}

	// Pointwise multiply
	for (std::size_t t = 0; t < m; ++t) {
		C_scratch[t] = CoeffCtx::mult(ARev[t], B0[t]);
	}

	// Inverse FFT
	fft<GT, F, CoeffCtx>(C_scratch, inv_roots, true);
	// (A*B0)[n+k] = D[n-2-k] for k = 0..n-2

	// Combine: out[k] = p[k] + q[k]
	for (std::size_t k = 0; k < n - 1; ++k) {
		out[k] = CoeffCtx::plus(out[k], C_scratch[n - 2 - k]);
	}
	// out[n-1] already has p[n-1], q[n-1] = 0, so no change needed
}

// Truncated Inverse NTT (DIF form) - computes only the first L outputs
// 
// IMPORTANT: DIF produces output in bit-reversed order. To get natural order output,
// we need to account for this. However, since we're truncating, we use a different
// approach: we implement the transpose of the forward DIT FFT.
//
// For inverse FFT with natural-order output, we use the same structure as forward
// but with inverse roots and post-normalization.
template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft_inverse_truncated(
	std::span<G> X,
	std::size_t m,
	std::size_t L,
	std::span<const F> inv_roots,
	const F& m_inv)
{
	BLST_ASSERT(X.size() >= m);
	BLST_ASSERT(L <= m);
	BLST_ASSERT(inv_roots.size() >= m);
	BLST_ASSERT((m & (m - 1)) == 0);  // power of 2

	if (L == 0) return;

	// For truncated inverse, we use the standard DIT structure (same as forward fft)
	// but with inv_roots. The key insight is that the standard fft with inv_roots
	// and post-scaling gives the inverse.
	//
	// We can truncate by not computing butterflies that only affect indices >= L
	// in the FINAL output. However, due to bit-reversal at the start, the mapping
	// between input indices and output indices is complex.
	//
	// For correctness, we do the full computation but this can be optimized later.
	// The bit-reversal permutation means early truncation is non-trivial.

	int n = static_cast<int>(m);

	// Bit-reversal permutation (same as forward)
	bit_reverse(X.subspan(0, m));

	// DIT butterflies with inverse roots
	for (int len = 2; len <= n; len <<= 1) {
		int half = len / 2;
		for (int i = 0; i < n; i += len) {
			for (int j = 0; j < half; j++) {
				F w = inv_roots[j * (n / len)];

				G u = X[i + j];
				G t = CoeffCtx::mult(X[i + j + half], w);
				X[i + j] = CoeffCtx::plus(u, t);
				X[i + j + half] = CoeffCtx::minus(u, t);
			}
		}
	}

	// Final normalization - only for prefix we care about
	for (std::size_t k = 0; k < L; ++k) {
		X[k] = CoeffCtx::mult(X[k], m_inv);
	}
}

// Full inverse using same structure as forward (for benchmarking comparison)
template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft_inverse_full(
	std::span<G> X,
	std::size_t m,
	std::span<const F> inv_roots,
	const F& m_inv)
{
	BLST_ASSERT(X.size() >= m);
	BLST_ASSERT(inv_roots.size() >= m);
	BLST_ASSERT((m & (m - 1)) == 0);

	int n = static_cast<int>(m);

	bit_reverse(X.subspan(0, m));

	for (int len = 2; len <= n; len <<= 1) {
		int half = len / 2;
		for (int i = 0; i < n; i += len) {
			for (int j = 0; j < half; j++) {
				F w = inv_roots[j * (n / len)];

				G u = X[i + j];
				G t = CoeffCtx::mult(X[i + j + half], w);
				X[i + j] = CoeffCtx::plus(u, t);
				X[i + j + half] = CoeffCtx::minus(u, t);
			}
		}
	}

	for (std::size_t k = 0; k < m; ++k) {
		X[k] = CoeffCtx::mult(X[k], m_inv);
	}
}

// Middle product using truncated inverse FFT
// Note: This version currently does full computation but only normalizes the prefix.
// True truncation with DIT structure is complex due to bit-reversal permutation.
template <typename G1, typename G2, typename GT, typename F, typename CoeffCtx>
void middle_product_fft_truncated(
	std::span<G1> A,
	std::span<G2> B0,
	std::span<G2> B1,
	bool skipFftB,
	std::size_t n,
	std::size_t m,
	std::span<GT> out,
	std::span<GT> C_scratch,
	std::span<const F> roots,
	std::span<const F> inv_roots,
	const F& omega,
	const F& m_inv,
	const GT& identity_gt,
	const G1& identity_g1,
	const G2& identity_g2)
{
	BLST_ASSERT(A.size() >= m);
	BLST_ASSERT(B0.size() >= m);
	BLST_ASSERT(B1.size() >= m);
	BLST_ASSERT(out.size() >= n);
	BLST_ASSERT(C_scratch.size() >= m);
	BLST_ASSERT(roots.size() >= m);
	BLST_ASSERT(inv_roots.size() >= m);
	BLST_ASSERT(m >= 2 * n);
	BLST_ASSERT((m & (m - 1)) == 0);

	// === Phase 1: Compute p = low_n(A * B1) ===

	for (std::size_t i = n; i < m; ++i) {
		A[i] = identity_g1;
		B1[i] = identity_g2;
	}

	fft<G1, F, CoeffCtx>(A, roots, false);
	if (skipFftB == false)
		fft<G2, F, CoeffCtx>(B1, roots, false);

	for (std::size_t t = 0; t < m; ++t) {
		C_scratch[t] = CoeffCtx::mult(A[t], B1[t]);
	}

	// Truncated inverse FFT: only need first n outputs
	fft_inverse_truncated<GT, F, CoeffCtx>(C_scratch, m, n, inv_roots, m_inv);

	for (std::size_t k = 0; k < n; ++k) {
		out[k] = C_scratch[k];
	}

	// === Phase 2: Compute q = high_n(A * B0) via reversal ===

	std::vector<G1> ARev(m);
	derive_ntt_of_reversed<G1, F, CoeffCtx>(
		std::span<const G1>(A.data(), m),
		ARev,
		n, m, omega, identity_g1);

	if (skipFftB == false)
	{
		for (std::size_t i = 0; i < n / 2; ++i) {
			std::swap(B0[i], B0[n - 1 - i]);
		}
		for (std::size_t i = n; i < m; ++i) {
			B0[i] = identity_g2;
		}

		fft<G2, F, CoeffCtx>(B0, roots, false);
	}

	for (std::size_t t = 0; t < m; ++t) {
		C_scratch[t] = CoeffCtx::mult(ARev[t], B0[t]);
	}

	// Truncated inverse FFT: only need first n-1 outputs
	fft_inverse_truncated<GT, F, CoeffCtx>(C_scratch, m, n - 1, inv_roots, m_inv);

	for (std::size_t k = 0; k < n - 1; ++k) {
		out[k] = CoeffCtx::plus(out[k], C_scratch[n - 2 - k]);
	}
}