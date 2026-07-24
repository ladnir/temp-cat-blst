
#pragma once
#include <algorithm>
#include <vector>
#include <span>
#include "blst/common/defines.h"
#include "blst/common/timer.h"

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

template <bool InverseFromForwardRoots, typename G, typename F, typename CoeffCtx>
void fft_butterflies(std::span<G> a, std::span<const F> roots) {
	int n = a.size();
	bit_reverse(a);

	for (int len = 2; len <= n; len <<= 1) {
		int half = len / 2;
		for (int i = 0; i < n; i += len) {

			// j == 0 case, this is 20% of the work
			// and we get to skip the multiplication
			// because w == 1.
			{ 
				G u = a[i];
				G t = a[i + half];
				a[i] = CoeffCtx::plus(u, t);
				a[i + half] = CoeffCtx::minus(u, t);
			}

			for (int j = 1; j < half; j++) {
				int root_index = j * (n / len);
				if constexpr (InverseFromForwardRoots) {
					root_index = n - root_index;
				}
				const F& w = roots[root_index];

				G u = a[i + j];
				G t = CoeffCtx::mult(a[i + j + half], w);
				a[i + j] = CoeffCtx::plus(u, t);
				a[i + j + half] = CoeffCtx::minus(u, t);
			}
		}
	}
}

template <typename G, typename F, typename CoeffCtx>
void fft_scale(std::span<G> a, const F& inverse_size) {
	for (G& x : a) {
		x = CoeffCtx::mult(x, inverse_size);
	}
}

template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft(std::span<G> a, std::span<const F> roots, bool invert) {
	fft_butterflies<false, G, F, CoeffCtx>(a, roots);

	// Integrated scaling for inverse FFT.
	if (invert) {
		F inv_n = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(a.size()));
		fft_scale<G, F, CoeffCtx>(a, inv_n);
	}
}

// Runs a transform from one forward root table. The inverse path selects
// reciprocal powers by index at compile time inside the butterfly loop and
// uses the caller's prepared 1/n scale.
template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft_with_forward_roots(
	std::span<G> a,
	std::span<const F> forward_roots,
	bool invert,
	const F& inverse_size) {
	if (invert) {
		fft_butterflies<true, G, F, CoeffCtx>(a, forward_roots);
		fft_scale<G, F, CoeffCtx>(a, inverse_size);
	} else {
		fft_butterflies<false, G, F, CoeffCtx>(a, forward_roots);
	}
}

// template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
// std::span<G> fft_stockham(
// 	std::span<const G> in,
// 	std::span<G> out0,
// 	std::span<G> out1,
// 	std::span<const F> roots,
// 	bool invert
// ) {
// 	int n = static_cast<int>(in.size());
// 	BLST_ASSERT(out0.size() == static_cast<std::size_t>(n));
// 	BLST_ASSERT(out1.size() == static_cast<std::size_t>(n));
// 	BLST_ASSERT((n & (n - 1)) == 0);

// 	// Ping-pong buffers.
// 	std::span<const G> src = in;
// 	std::span<G> dst = out0;

// 	// Stages: m = 2, 4, 8, ..., n
// 	for (int m = 2; m <= n; m <<= 1) {
// 		int mh = m >> 1;
// 		int r = n / m;           // "radix" in Stockham indexing
// 		int stride = r;          // twiddle index step = n/m

// 		// Core Stockham mapping:
// 		// Input index:  k*m + j
// 		// Output index: j*r + k          (and + n/2 for the "upper" output)
// 		for (int k = 0; k < r; k++) {
// 			int base_in = k * m;

// 			{ // j == 0 case: w == 1
// 				G u = src[base_in + 0];
// 				G v = src[base_in + mh];
// 				int base_out = 0 * r + k;
// 				dst[base_out] = CoeffCtx::plus(u, v);
// 				dst[base_out + (n >> 1)] = CoeffCtx::minus(u, v);
// 			}

// 			for (int j = 1; j < mh; j++) {
// 				G u = src[base_in + j];
// 				G v = src[base_in + j + mh];

// 				// Twiddle w = ω^(j * n/m) = roots[j * r]
// 				// If you want inverse, pass inv_roots as roots.
// 				F w = roots[j * stride];

// 				// Short-circuit w == 1 by index: j == 0 => w == roots[0] == 1.
// 				G t = CoeffCtx::mult(v, w);

// 				int base_out = j * r + k;

// 				dst[base_out] = CoeffCtx::plus(u, t);
// 				dst[base_out + (n >> 1)] = CoeffCtx::minus(u, t);
// 			}
// 		}

// 		// Swap src/dst for next stage.
// 		if (dst.data() == out0.data()) {
// 			src = out0;
// 			dst = out1;
// 		}
// 		else {
// 			src = out1;
// 			dst = out0;
// 		}
// 	}

// 	if (src.data() == out0.data())
// 	{
// 		dst = out0;
// 	}
// 	else
// 	{
// 		dst = out1;
// 	}

// 	//// After the final stage, the result lives in `src`.
// 	//// Ensure it ends up in `out`.
// 	//if (src.data() != out.data()) {
// 	//	for (int i = 0; i < n; i++) {
// 	//		out[i] = src[i];
// 	//	}
// 	//}

// 	// Scaling for inverse FFT: multiply every output by 1/n.
// 	// (As before: this assumes you provided inverse roots when invert==true.)
// 	if (invert) {
// 		F inv_n = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(static_cast<std::uint64_t>(n)));
// 		for (int i = 0; i < n; i++) {
// 			dst[i] = CoeffCtx::mult(src[i], inv_n);
// 		}
// 	}

// 	return dst;
// }


//template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
//void fft_dif(std::span<G> a, std::span<const F> roots, bool invert) {
//	int n = static_cast<int>(a.size());
//	BLST_ASSERT(n > 0);
//	BLST_ASSERT((n & (n - 1)) == 0);
//
//	int total = 0;
//	int zeroCount = 0;
//	int oneCount = 0;
//	int negOneCount = 0;
//
//	F one = CoeffCtx::template from_u64<F>(1);
//	F zero = CoeffCtx::template from_u64<F>(0);
//	F negOne;
//	F::neg(negOne, one);
//
//	// DIF: stages go from large blocks down to small blocks.
//	// No bit-reverse at the start; do it at the end.
//	for (int len = n; len >= 2; len >>= 1) {
//		int half = len / 2;
//		int stride = n / len;
//
//		for (int i = 0; i < n; i += len) {
//			for (int j = 0; j < half; j++) {
//				// Inputs
//				G u = a[i + j];
//				G v = a[i + j + half];
//
//				// Butterfly core
//				// a0 = u + v
//				// a1 = (u - v) * w
//				G sum = CoeffCtx::plus(u, v);
//				G diff = CoeffCtx::minus(u, v);
//
//				F w = roots[j * stride];
//
//				++total;
//
//				a[i + j] = sum;
//
//				// Hot-path short-circuits: these matter a lot for you.
//				if (F::equal(w, one)) {
//					oneCount++;
//					a[i + j + half] = diff;
//				}
//				else if (F::equal(w, negOne)) {
//					// Multiply by -1 in the scalar field == group negation.
//					negOneCount++;
//					a[i + j + half] = CoeffCtx::neg(diff);
//				}
//				//else if (F::equal(w, zero)) {
//				//	// Should not occur for roots of unity in a correct NTT, but kept for completeness.
//				//	a[i + j + half] = G::infinity();
//				//}
//				else {
//					a[i + j + half] = CoeffCtx::mult(diff, w);
//				}
//			}
//		}
//	}
//
//	//printf(
//	//	"FFT(DIF) total butterflies: %d, special cases: %d (zero: %d, one: %d, negOne: %d)\n",
//	//	total,
//	//	zeroCount + oneCount + negOneCount,
//	//	zeroCount,
//	//	oneCount,
//	//	negOneCount
//	//);
//
//	// Bit-reverse at the end for DIF to match natural-order output
//	bit_reverse(a);
//
//	// Scaling for inverse FFT (same as your code)
//	if (invert) {
//		F inv_n = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(n));
//		for (G& x : a) {
//			x = CoeffCtx::mult(x, inv_n);
//		}
//	}
//}


// ============================================================================
// Radix-4 FFT (DIT, in-place)
// 
// Reduces number of stages by half compared to radix-2.
// Requires N to be a power of 4 for pure radix-4, or handles mixed radix.
// ============================================================================
template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft_radix4(std::span<G> a, std::span<const F> roots, bool invert) {
	std::size_t n = a.size();
	BLST_ASSERT((n & (n - 1)) == 0); // n must be power of 2

	bit_reverse(a);

	// Precompute i4 = ω^(N/4), which satisfies i4^2 = -1 in the field
	// For NTT, this is the 4th root of unity
	F i4 = roots[n / 4];

	// Handle case where n is not a power of 4: do one radix-2 stage first
	std::size_t start_len = 4;
	if ((n & 3) != 0) {
		// n is 2 mod 4, do a radix-2 stage
		for (std::size_t i = 0; i < n; i += 2) {
			G u = a[i];
			G v = a[i + 1];
			a[i] = CoeffCtx::plus(u, v);
			a[i + 1] = CoeffCtx::minus(u, v);
		}
		start_len = 4;
	}
	else if (n == 2) {
		// Special case: just radix-2
		G u = a[0];
		G v = a[1];
		a[0] = CoeffCtx::plus(u, v);
		a[1] = CoeffCtx::minus(u, v);
		if (invert) {
			F inv_n = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(n));
			a[0] = CoeffCtx::mult(a[0], inv_n);
			a[1] = CoeffCtx::mult(a[1], inv_n);
		}
		return;
	}

	// Main radix-4 stages
	for (std::size_t m = start_len; m <= n; m <<= 2) {
		std::size_t m4 = m / 4;
		std::size_t stride = n / m;

		for (std::size_t k = 0; k < n; k += m) {
			for (std::size_t j = 0; j < m4; j++) {
				// Twiddle factors: wj = ω^(j * stride)
				F wj = roots[j * stride];
				F wj2 = CoeffCtx::mult(wj, wj);     // ω^(2j * stride)
				F wj3 = CoeffCtx::mult(wj2, wj);    // ω^(3j * stride)

				// Gather butterfly inputs
				G a0 = a[k + j];
				G a1 = a[k + j + m4];
				G a2 = a[k + j + 2 * m4];
				G a3 = a[k + j + 3 * m4];

				// Apply twiddles to legs 1, 2, 3
				G b0 = a0;
				G b1 = CoeffCtx::mult(a1, wj);
				G b2 = CoeffCtx::mult(a2, wj2);
				G b3 = CoeffCtx::mult(a3, wj3);

				// 4-point DFT kernel
				// s0 = b0 + b2, d0 = b0 - b2
				// s1 = b1 + b3, d1 = b1 - b3
				// t = d1 * i4 (where i4^2 = -1)
				G s0 = CoeffCtx::plus(b0, b2);
				G d0 = CoeffCtx::minus(b0, b2);
				G s1 = CoeffCtx::plus(b1, b3);
				G d1 = CoeffCtx::minus(b1, b3);
				G t = CoeffCtx::mult(d1, i4);

				// Output
				a[k + j] = CoeffCtx::plus(s0, s1);           // y0
				a[k + j + m4] = CoeffCtx::minus(d0, t);      // y1
				a[k + j + 2 * m4] = CoeffCtx::minus(s0, s1); // y2
				a[k + j + 3 * m4] = CoeffCtx::plus(d0, t);   // y3
			}
		}
	}

	// Scaling for inverse FFT
	if (invert) {
		F inv_n = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(n));
		for (G& x : a) {
			x = CoeffCtx::mult(x, inv_n);
		}
	}
}

// ============================================================================
// Radix-8 FFT (DIT, in-place)
// 
// Further reduces stages compared to radix-4.
// Requires N to be a power of 8 for pure radix-8, or handles mixed radix.
// ============================================================================
template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft_radix8(std::span<G> a, std::span<const F> roots, bool invert) {
	std::size_t n = a.size();
	BLST_ASSERT((n & (n - 1)) == 0); // n must be power of 2

	bit_reverse(a);

	// Precompute 8th roots needed for the radix-8 kernel
	// z1 = ω^(N/8), z2 = ω^(N/4), z3 = ω^(3N/8)
	F z1 = roots[n / 8];                      // ω^(N/8) - primitive 8th root
	F z2 = roots[n / 4];                      // ω^(N/4) - i4, the 4th root
	F z3 = CoeffCtx::mult(z1, z2);            // ω^(3N/8)

	// Handle sizes that aren't powers of 8 by doing smaller radix stages first
	std::size_t start_len = 8;
	std::size_t log2n = 0;
	for (std::size_t tmp = n; tmp > 1; tmp >>= 1) log2n++;

	// Do radix-2 or radix-4 stages as needed to reach a power of 8
	std::size_t remainder = log2n % 3;

	if (remainder == 1) {
		// One radix-2 stage
		for (std::size_t i = 0; i < n; i += 2) {
			G u = a[i];
			G v = a[i + 1];
			a[i] = CoeffCtx::plus(u, v);
			a[i + 1] = CoeffCtx::minus(u, v);
		}
		start_len = 4;

		// Then one radix-4 stage
		if (n >= 4) {
			std::size_t m = 4;
			std::size_t m4 = 1;
			std::size_t stride = n / m;
			F i4 = z2;

			for (std::size_t k = 0; k < n; k += m) {
				for (std::size_t j = 0; j < m4; j++) {
					F wj = roots[j * stride];
					F wj2 = CoeffCtx::mult(wj, wj);
					F wj3 = CoeffCtx::mult(wj2, wj);

					G a0 = a[k + j];
					G a1 = a[k + j + m4];
					G a2 = a[k + j + 2 * m4];
					G a3 = a[k + j + 3 * m4];

					G b0 = a0;
					G b1 = CoeffCtx::mult(a1, wj);
					G b2 = CoeffCtx::mult(a2, wj2);
					G b3 = CoeffCtx::mult(a3, wj3);

					G s0 = CoeffCtx::plus(b0, b2);
					G d0 = CoeffCtx::minus(b0, b2);
					G s1 = CoeffCtx::plus(b1, b3);
					G d1 = CoeffCtx::minus(b1, b3);
					G t = CoeffCtx::mult(d1, i4);

					a[k + j] = CoeffCtx::plus(s0, s1);
					a[k + j + m4] = CoeffCtx::minus(d0, t);
					a[k + j + 2 * m4] = CoeffCtx::minus(s0, s1);
					a[k + j + 3 * m4] = CoeffCtx::plus(d0, t);
				}
			}
		}
		start_len = 8;
	}
	else if (remainder == 2) {
		// One radix-4 stage (4 = 2^2)
		std::size_t m = 4;
		std::size_t m4 = 1;
		std::size_t stride = n / m;
		F i4 = z2;

		for (std::size_t k = 0; k < n; k += m) {
			for (std::size_t j = 0; j < m4; j++) {
				F wj = roots[j * stride];
				F wj2 = CoeffCtx::mult(wj, wj);
				F wj3 = CoeffCtx::mult(wj2, wj);

				G a0 = a[k + j];
				G a1 = a[k + j + m4];
				G a2 = a[k + j + 2 * m4];
				G a3 = a[k + j + 3 * m4];

				G b0 = a0;
				G b1 = CoeffCtx::mult(a1, wj);
				G b2 = CoeffCtx::mult(a2, wj2);
				G b3 = CoeffCtx::mult(a3, wj3);

				G s0 = CoeffCtx::plus(b0, b2);
				G d0 = CoeffCtx::minus(b0, b2);
				G s1 = CoeffCtx::plus(b1, b3);
				G d1 = CoeffCtx::minus(b1, b3);
				G t = CoeffCtx::mult(d1, i4);

				a[k + j] = CoeffCtx::plus(s0, s1);
				a[k + j + m4] = CoeffCtx::minus(d0, t);
				a[k + j + 2 * m4] = CoeffCtx::minus(s0, s1);
				a[k + j + 3 * m4] = CoeffCtx::plus(d0, t);
			}
		}
		start_len = 8;
	}

	// Main radix-8 stages
	for (std::size_t m = start_len; m <= n; m <<= 3) {
		std::size_t m8 = m / 8;
		std::size_t stride = n / m;

		for (std::size_t k = 0; k < n; k += m) {
			for (std::size_t j = 0; j < m8; j++) {
				// Twiddle factor wj = ω^(j * stride)
				F wj = roots[j * stride];

				// Compute powers of wj for legs 1-7
				F wj2 = CoeffCtx::mult(wj, wj);
				F wj3 = CoeffCtx::mult(wj2, wj);
				F wj4 = CoeffCtx::mult(wj2, wj2);
				F wj5 = CoeffCtx::mult(wj4, wj);
				F wj6 = CoeffCtx::mult(wj3, wj3);
				F wj7 = CoeffCtx::mult(wj6, wj);

				// Gather butterfly inputs
				G a0 = a[k + j];
				G a1 = a[k + j + m8];
				G a2 = a[k + j + 2 * m8];
				G a3 = a[k + j + 3 * m8];
				G a4 = a[k + j + 4 * m8];
				G a5 = a[k + j + 5 * m8];
				G a6 = a[k + j + 6 * m8];
				G a7 = a[k + j + 7 * m8];

				// Apply twiddles
				G b0 = a0;
				G b1 = CoeffCtx::mult(a1, wj);
				G b2 = CoeffCtx::mult(a2, wj2);
				G b3 = CoeffCtx::mult(a3, wj3);
				G b4 = CoeffCtx::mult(a4, wj4);
				G b5 = CoeffCtx::mult(a5, wj5);
				G b6 = CoeffCtx::mult(a6, wj6);
				G b7 = CoeffCtx::mult(a7, wj7);

				// 8-point DFT via two 4-point DFTs
				// Evens: e0=b0, e1=b2, e2=b4, e3=b6
				// Odds:  o0=b1, o1=b3, o2=b5, o3=b7

				// DFT4 of evens
				G es0 = CoeffCtx::plus(b0, b4);
				G ed0 = CoeffCtx::minus(b0, b4);
				G es1 = CoeffCtx::plus(b2, b6);
				G ed1 = CoeffCtx::minus(b2, b6);
				G et = CoeffCtx::mult(ed1, z2); // ed1 * i4

				G E0 = CoeffCtx::plus(es0, es1);
				G E1 = CoeffCtx::minus(ed0, et);
				G E2 = CoeffCtx::minus(es0, es1);
				G E3 = CoeffCtx::plus(ed0, et);

				// DFT4 of odds
				G os0 = CoeffCtx::plus(b1, b5);
				G od0 = CoeffCtx::minus(b1, b5);
				G os1 = CoeffCtx::plus(b3, b7);
				G od1 = CoeffCtx::minus(b3, b7);
				G ot = CoeffCtx::mult(od1, z2); // od1 * i4

				G O0 = CoeffCtx::plus(os0, os1);
				G O1 = CoeffCtx::minus(od0, ot);
				G O2 = CoeffCtx::minus(os0, os1);
				G O3 = CoeffCtx::plus(od0, ot);

				// Combine with 8th roots: y[k] = E[k] + O[k] * ζ^k
				// y[k+4] = E[k] - O[k] * ζ^k
				G O0z0 = O0;                          // O0 * ζ^0 = O0
				G O1z1 = CoeffCtx::mult(O1, z1);      // O1 * ζ^1
				G O2z2 = CoeffCtx::mult(O2, z2);      // O2 * ζ^2 = O2 * i4
				G O3z3 = CoeffCtx::mult(O3, z3);      // O3 * ζ^3

				a[k + j] = CoeffCtx::plus(E0, O0z0);
				a[k + j + m8] = CoeffCtx::plus(E1, O1z1);
				a[k + j + 2 * m8] = CoeffCtx::plus(E2, O2z2);
				a[k + j + 3 * m8] = CoeffCtx::plus(E3, O3z3);
				a[k + j + 4 * m8] = CoeffCtx::minus(E0, O0z0);
				a[k + j + 5 * m8] = CoeffCtx::minus(E1, O1z1);
				a[k + j + 6 * m8] = CoeffCtx::minus(E2, O2z2);
				a[k + j + 7 * m8] = CoeffCtx::minus(E3, O3z3);
			}
		}
	}

	// Scaling for inverse FFT
	if (invert) {
		F inv_n = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(n));
		for (G& x : a) {
			x = CoeffCtx::mult(x, inv_n);
		}
	}
}


// // ============================================================================
// // Adaptive FFT - selects best radix based on input size
// // ============================================================================
// template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
// void fft_adaptive(std::span<G> a, std::span<const F> roots, bool invert) {
// 	std::size_t n = a.size();

// 	if (n < 8) {
// 		// Fall back to radix-2 for very small sizes
// 		fft<G, F, CoeffCtx>(a, roots, invert);
// 	}
// 	else if (n >= 64) {
// 		// Use radix-8 for larger sizes (fewer stages)
// 		fft_radix8<G, F, CoeffCtx>(a, roots, invert);
// 	}
// 	else {
// 		// Use radix-4 for medium sizes
// 		fft_radix4<G, F, CoeffCtx>(a, roots, invert);
// 	}
// }

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
template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft_forward_truncated_prefix(
        std::span<G> a,
        std::span<const F> roots,
        std::size_t xn,
        const struct FftTruncatedForwardPlan& plan,
        blst::TimingStats* timings = nullptr,
        std::string_view metric_prefix = "fft_trunc");

template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft_forward_truncated_prefix(
        std::span<G> a,
        std::span<const F> roots,
        std::size_t xn,
        blst::TimingStats* timings = nullptr,
        std::string_view metric_prefix = "fft_trunc");

struct FftTruncatedForwardPlan
{
        struct StageWork
        {
                std::vector<std::size_t> left_j0;
                std::vector<std::size_t> left_nonzero;
                std::vector<std::size_t> twiddle_nonzero;
        };

        std::size_t n = 0;
        std::size_t xn = 0;
        std::size_t log2_n = 0;
        std::size_t butterflies_total = 0;
        std::size_t butterflies_executed = 0;
        std::size_t butterflies_skipped = 0;
        std::vector<StageWork> stage_work;
};

struct FftInverseTruncatedRangePlan
{
        std::size_t n = 0;
        std::size_t start = 0;
        std::size_t length = 0;
        std::size_t stages = 0;
        std::size_t butterflies_total = 0;
        std::size_t butterflies_executed = 0;
        std::size_t butterflies_skipped = 0;
        std::vector<std::vector<std::size_t>> stage_left;
};

inline std::size_t fft_log2_size(std::size_t n)
{
        std::size_t log2_n = 0;
        for (std::size_t t = n; t > 1; t >>= 1) {
                ++log2_n;
        }
        return log2_n;
}

inline FftTruncatedForwardPlan make_fft_truncated_forward_plan(std::size_t n, std::size_t xn)
{
        BLST_ASSERT(n > 0);
        BLST_ASSERT((n & (n - 1)) == 0);
        BLST_ASSERT(xn > 0);
        BLST_ASSERT(xn <= n);

        FftTruncatedForwardPlan plan;
        plan.n = n;
        plan.xn = xn;
        plan.log2_n = fft_log2_size(n);
        plan.butterflies_total = (n / 2) * plan.log2_n;

        if (xn == n) {
                plan.butterflies_executed = plan.butterflies_total;
                return plan;
        }

        plan.stage_work.resize(plan.log2_n);
        std::vector<unsigned char> active(n, 0);
        std::vector<unsigned char> next_active(n, 0);
        for (std::size_t i = 0; i < xn; ++i) {
                active[i] = 1;
        }
        bit_reverse<unsigned char>(std::span<unsigned char>(active));

        std::size_t stage_idx = 0;
        for (std::size_t len = 2; len <= n; len <<= 1, ++stage_idx) {
                const std::size_t half = len >> 1;
                const std::size_t stride = n / len;
                auto& stage = plan.stage_work[stage_idx];
                stage.left_j0.reserve(n / len);
                stage.left_nonzero.reserve((n / 2) - (n / len));
                stage.twiddle_nonzero.reserve((n / 2) - (n / len));
                std::fill(next_active.begin(), next_active.end(), static_cast<unsigned char>(0));

                for (std::size_t i = 0; i < n; i += len) {
                        for (std::size_t j = 0; j < half; ++j) {
                                const std::size_t left = i + j;
                                const std::size_t right = left + half;
                                if (!(active[left] || active[right])) {
                                        continue;
                                }
                                if (j == 0) {
                                        stage.left_j0.push_back(left);
                                }
                                else {
                                        stage.left_nonzero.push_back(left);
                                        stage.twiddle_nonzero.push_back(j * stride);
                                }
                                next_active[left] = 1;
                                next_active[right] = 1;
                        }
                }
                plan.butterflies_executed += stage.left_j0.size() + stage.left_nonzero.size();
                active.swap(next_active);
        }

        plan.butterflies_skipped = plan.butterflies_total - plan.butterflies_executed;
        return plan;
}

inline FftInverseTruncatedRangePlan make_fft_inverse_truncated_range_plan(
        std::size_t n,
        std::size_t start,
        std::size_t length)
{
        BLST_ASSERT(n > 0);
        BLST_ASSERT((n & (n - 1)) == 0);
        BLST_ASSERT(length > 0);
        BLST_ASSERT(start < n);
        BLST_ASSERT(start + length <= n);

        FftInverseTruncatedRangePlan plan;
        plan.n = n;
        plan.start = start;
        plan.length = length;
        plan.stages = fft_log2_size(n);
        plan.butterflies_total = (n / 2) * plan.stages;
        plan.stage_left.resize(plan.stages);

        std::vector<unsigned char> out_need(n, 0);
        std::vector<unsigned char> in_need(n, 0);
        for (std::size_t i = 0; i < length; ++i) {
                out_need[start + i] = 1;
        }

        for (std::size_t k = plan.stages; k > 0; --k) {
                const std::size_t len = std::size_t(1) << k;
                const std::size_t half = len >> 1;
                auto& lefts = plan.stage_left[k - 1];
                lefts.reserve(n / 2);
                std::fill(in_need.begin(), in_need.end(), static_cast<unsigned char>(0));

                for (std::size_t i = 0; i < n; i += len) {
                        for (std::size_t j = 0; j < half; ++j) {
                                const std::size_t left = i + j;
                                const std::size_t right = left + half;
                                if (!(out_need[left] || out_need[right])) {
                                        continue;
                                }
                                lefts.push_back(left);
                                in_need[left] = 1;
                                in_need[right] = 1;
                        }
                }

                plan.butterflies_executed += lefts.size();
                out_need.swap(in_need);
        }

        plan.butterflies_skipped = plan.butterflies_total - plan.butterflies_executed;
        return plan;
}

template <typename G1, typename G2, typename GT, typename F, typename CoeffCtx = DefaultCoeffCtx>
void poly_mult_fft(
	std::span<G1> A,
	std::span<G2> B,
	bool skipBFft,
	std::span<GT> C,
	std::span<const F> roots,
	std::span<const F> inv_roots,
	blst::TimingStats* timings = nullptr,
	std::size_t xnA = 0,
	std::size_t xnB = 0)
{
	std::size_t n = A.size();
	if (xnA == 0) {
		xnA = n;
	}
	if (xnB == 0) {
		xnB = n;
	}
	BLST_ASSERT(xnA <= n);
	BLST_ASSERT(xnB <= n);

	// Forward FFT on both polynomials (in-place)
	{
		blst::ScopedTimer timer(timings, "poly_mult.fft_a_ns");
		fft_forward_truncated_prefix<G1, F, CoeffCtx>(A, roots, xnA, timings, "poly_mult.fft_a");
	}

	if (skipBFft == false) {
		blst::ScopedTimer timer(timings, "poly_mult.fft_b_ns");
		fft_forward_truncated_prefix<G2, F, CoeffCtx>(B, roots, xnB, timings, "poly_mult.fft_b");
	}

	// Pointwise multiplication in evaluation domain
	{
		blst::ScopedTimer timer(timings, "poly_mult.pointwise_ns");
		for (std::size_t i = 0; i < n; ++i) {
			C[i] = CoeffCtx::mult(A[i], B[i]);
		}
	}

	// Inverse FFT to get product coefficients
	{
		blst::ScopedTimer timer(timings, "poly_mult.ifft_ns");
		fft<GT, F, CoeffCtx>(C, inv_roots, true);
	}
}

template <typename G, typename F, typename CoeffCtx>
void fft_forward_truncated_prefix(
        std::span<G> a,
        std::span<const F> roots,
        std::size_t xn,
        const FftTruncatedForwardPlan& plan,
        blst::TimingStats* timings,
        std::string_view metric_prefix) {
        const std::size_t n = a.size();
        BLST_ASSERT(n > 0);
        BLST_ASSERT((n & (n - 1)) == 0);
        BLST_ASSERT(roots.size() >= n);
        BLST_ASSERT(xn > 0);
        BLST_ASSERT(xn <= n);
        BLST_ASSERT(plan.n == n);
        BLST_ASSERT(plan.xn == xn);

        if (xn == n) {
                fft<G, F, CoeffCtx>(a, roots, false);
                if (timings != nullptr) {
                        timings->add(std::string(metric_prefix) + ".xn", xn);
                        timings->add(std::string(metric_prefix) + ".butterflies_total", plan.butterflies_total);
                        timings->add(std::string(metric_prefix) + ".butterflies_executed", plan.butterflies_executed);
                        timings->add(std::string(metric_prefix) + ".butterflies_skipped", plan.butterflies_skipped);
                }
                return;
        }

        bit_reverse(a);

        for (std::size_t stage = 0, len = 2; len <= n; len <<= 1, ++stage) {
                const std::size_t half = len >> 1;
                const auto& stage_work = plan.stage_work[stage];

                for (std::size_t idx = 0; idx < stage_work.left_j0.size(); ++idx) {
                        const std::size_t left = stage_work.left_j0[idx];
                        const std::size_t right = left + half;
                        G u = a[left];
                        G t = a[right];
                        a[left] = CoeffCtx::plus(u, t);
                        a[right] = CoeffCtx::minus(u, t);
                }

                for (std::size_t idx = 0; idx < stage_work.left_nonzero.size(); ++idx) {
                        const std::size_t left = stage_work.left_nonzero[idx];
                        const std::size_t right = left + half;
                        const F& w = roots[stage_work.twiddle_nonzero[idx]];
                        G u = a[left];
                        G t = CoeffCtx::mult(a[right], w);
                        a[left] = CoeffCtx::plus(u, t);
                        a[right] = CoeffCtx::minus(u, t);
                }
        }

        if (timings != nullptr) {
                timings->add(std::string(metric_prefix) + ".xn", xn);
                timings->add(std::string(metric_prefix) + ".butterflies_total", plan.butterflies_total);
                timings->add(std::string(metric_prefix) + ".butterflies_executed", plan.butterflies_executed);
                timings->add(std::string(metric_prefix) + ".butterflies_skipped", plan.butterflies_skipped);
        }
}

template <typename G, typename F, typename CoeffCtx>
void fft_forward_truncated_prefix(
        std::span<G> a,
        std::span<const F> roots,
        std::size_t xn,
        blst::TimingStats* timings,
        std::string_view metric_prefix)
{
        const FftTruncatedForwardPlan plan = make_fft_truncated_forward_plan(a.size(), xn);
        fft_forward_truncated_prefix<G, F, CoeffCtx>(a, roots, xn, plan, timings, metric_prefix);
}

template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft_inverse_truncated_range(
        std::span<G> x,
        std::span<const F> inv_roots,
        std::size_t start,
        std::size_t length,
        const FftInverseTruncatedRangePlan& plan,
        blst::TimingStats* timings = nullptr,
        std::string_view metric_prefix = "ifft_trunc_range")
{
        const std::size_t n = x.size();
        BLST_ASSERT(n > 0);
        BLST_ASSERT((n & (n - 1)) == 0);
        BLST_ASSERT(inv_roots.size() >= n);
        BLST_ASSERT(length > 0);
        BLST_ASSERT(start < n);
        BLST_ASSERT(start + length <= n);
        BLST_ASSERT(plan.n == n);
        BLST_ASSERT(plan.start == start);
        BLST_ASSERT(plan.length == length);

        bit_reverse(x);

        for (std::size_t stage = 0, len = 2; len <= n; len <<= 1, ++stage) {
                const std::size_t half = len >> 1;
                const std::size_t stride = n / len;
                const auto& lefts = plan.stage_left[stage];

                for (std::size_t idx = 0; idx < lefts.size(); ++idx) {
                        const std::size_t left = lefts[idx];
                        const std::size_t right = left + half;
                        const std::size_t j = left & (half - 1);

                        if (j == 0) {
                                G u = x[left];
                                G t = x[right];
                                x[left] = CoeffCtx::plus(u, t);
                                x[right] = CoeffCtx::minus(u, t);
                        }
                        else {
                                const F& w = inv_roots[j * stride];
                                G u = x[left];
                                G t = CoeffCtx::mult(x[right], w);
                                x[left] = CoeffCtx::plus(u, t);
                                x[right] = CoeffCtx::minus(u, t);
                        }
                }
        }

        F inv_n = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(n));
        for (std::size_t i = 0; i < length; ++i) {
                const std::size_t idx = start + i;
                x[idx] = CoeffCtx::mult(x[idx], inv_n);
        }

        if (timings != nullptr) {
                timings->add(std::string(metric_prefix) + ".start", start);
                timings->add(std::string(metric_prefix) + ".length", length);
                timings->add(std::string(metric_prefix) + ".butterflies_total", plan.butterflies_total);
                timings->add(std::string(metric_prefix) + ".butterflies_executed", plan.butterflies_executed);
                timings->add(std::string(metric_prefix) + ".butterflies_skipped", plan.butterflies_skipped);
        }
}

template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft_inverse_truncated_range(
        std::span<G> x,
        std::span<const F> inv_roots,
	std::size_t start,
	std::size_t length,
	blst::TimingStats* timings = nullptr,
	std::string_view metric_prefix = "ifft_trunc_range")
{
	const std::size_t n = x.size();
	BLST_ASSERT(n > 0);
	BLST_ASSERT((n & (n - 1)) == 0);
	BLST_ASSERT(inv_roots.size() >= n);
	BLST_ASSERT(length > 0);
	BLST_ASSERT(start < n);
	BLST_ASSERT(start + length <= n);

	const FftInverseTruncatedRangePlan plan =
			make_fft_inverse_truncated_range_plan(n, start, length);
	fft_inverse_truncated_range<G, F, CoeffCtx>(
			x, inv_roots, start, length, plan, timings, metric_prefix);
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


// // Computes the middle n coefficients of C = A * B using the reversal trick
// // A has n coefficients (type G1), B has 2n coefficients split as B0 and B1
// // Writes coefficients c[n], c[n+1], ..., c[2n-1] to out (0-indexed)
// //
// // OPTIMIZATION: Uses frequency reversal + phase trick to eliminate one G1 NTT
// // Instead of computing NTT(rev(A)) directly, we derive it from NTT(A):
// //   NTT(rev(A))[k] = ω^((n-1)*k) * NTT(A)[(-k) mod m]
// //
// // This saves one full G1 forward NTT (m log m scalar mults in G1)
// // and replaces it with m scalar mults in G1 (linear cost).
// //
// // IMPORTANT: This function modifies A, B0, B1 in-place. Caller must copy if needed.
// //
// // Parameters:
// //   A         - input polynomial, size m, input in [0,n), rest used as scratch, MODIFIED
// //   B0        - input polynomial (low half of B), size m, input in [0,n), MODIFIED
// //   B1        - input polynomial (high half of B), size m, input in [0,n), MODIFIED
// //   n         - logical size of A, B0, B1 each have n valid coefficients
// //   m         - FFT size, must be power of 2 and >= 2n
// //   out       - output buffer, size >= n
// //   C_scratch - scratch buffer for FFT results, size m
// //   roots     - precomputed roots of unity, size m
// //   inv_roots - precomputed inverse roots, size m
// //   omega     - primitive m-th root of unity (roots[1])
// //   identity_gt - additive identity for GT
// //   identity_g1 - additive identity for G1
// //   identity_g2 - additive identity for G2
// template <typename G1, typename G2, typename GT, typename F, typename CoeffCtx>
// void middle_product_fft(
// 	std::span<G1> A,
// 	std::span<G2> B0,
// 	std::span<G2> B1,
// 	bool skipFftB,
// 	std::size_t n,
// 	std::size_t m,
// 	std::span<GT> out,
// 	std::span<GT> C_scratch,
// 	std::span<const F> roots,
// 	std::span<const F> inv_roots,
// 	const F& omega,
// 	const GT& identity_gt,
// 	const G1& identity_g1,
// 	const G2& identity_g2,
// 	blst::TimingStats* timings = nullptr)
// {
// 	BLST_ASSERT(A.size() >= m);
// 	BLST_ASSERT(B0.size() >= m);
// 	BLST_ASSERT(B1.size() >= m);
// 	BLST_ASSERT(out.size() >= n);
// 	BLST_ASSERT(C_scratch.size() >= m);
// 	BLST_ASSERT(roots.size() >= m);
// 	BLST_ASSERT(inv_roots.size() >= m);
// 	BLST_ASSERT(m >= 2 * n);
// 	BLST_ASSERT((m & (m - 1)) == 0);

// 	// === Phase 1: Compute p = low_n(A * B1) ===

// 	// Pad A and B1 with identity
// 	{
// 		blst::ScopedTimer timer(timings, "middle.phase1_prepare_ns");
// 		for (std::size_t i = n; i < m; ++i) {
// 			A[i] = identity_g1;
// 			if (skipFftB == false) {
// 				B1[i] = identity_g2;
// 			}
// 		}
// 	}

// 	// Forward FFT on A
// 	{
// 		blst::ScopedTimer timer(timings, "middle.phase1_fft_a_ns");
// 		fft<G1, F, CoeffCtx>(A, roots, false);
// 	}

// 	if (skipFftB == false)
// 	{
// 		// Forward FFT on B1
// 		blst::ScopedTimer timer(timings, "middle.phase1_fft_b1_ns");
// 		fft<G2, F, CoeffCtx>(B1, roots, false);
// 	}

// 	// Pointwise multiply -> C_scratch
// 	{
// 		blst::ScopedTimer timer(timings, "middle.phase1_pointwise_ns");
// 		for (std::size_t t = 0; t < m; ++t) {
// 			C_scratch[t] = CoeffCtx::mult(A[t], B1[t]);
// 		}
// 	}

// 	// Inverse FFT to get (A * B1) coefficients
// 	{
// 		blst::ScopedTimer timer(timings, "middle.phase1_ifft_ns");
// 		fft<GT, F, CoeffCtx>(C_scratch, inv_roots, true);
// 	}

// 	// Save p[k] = C_scratch[k] for k = 0..n-1 into out
// 	{
// 		blst::ScopedTimer timer(timings, "middle.phase1_copy_ns");
// 		for (std::size_t k = 0; k < n; ++k) {
// 			out[k] = C_scratch[k];
// 		}
// 	}

// 	// === Phase 2: Compute q = high_n(A * B0) via reversal ===

// 	std::vector<G1> ARev(m);

// 	// Derive NTT(rev(A)) from NTT(A) using the phase trick
// 	// A_rev_hat[k] = ω^((n-1)*k) * A_hat[(-k) mod m]
// 	{
// 		blst::ScopedTimer timer(timings, "middle.phase2_derive_ntt_ns");
// 		derive_ntt_of_reversed<G1, F, CoeffCtx>(
// 			std::span<const G1>(A.data(), m),
// 			ARev,
// 			n, m, omega, identity_g1);
// 	}

// 	if (skipFftB == false)
// 	{
// 		// Reverse B0 in-place and pad
// 		{
// 			blst::ScopedTimer timer(timings, "middle.phase2_prepare_b0_ns");
// 			for (std::size_t i = 0; i < n / 2; ++i) {
// 				std::swap(B0[i], B0[n - 1 - i]);
// 			}
// 			for (std::size_t i = n; i < m; ++i) {
// 				B0[i] = identity_g2;
// 			}
// 		}
// 		// Forward FFT on rev(B0)
// 		{
// 			blst::ScopedTimer timer(timings, "middle.phase2_fft_b0_ns");
// 			fft<G2, F, CoeffCtx>(B0, roots, false);
// 		}
// 	}

// 	// Pointwise multiply
// 	{
// 		blst::ScopedTimer timer(timings, "middle.phase2_pointwise_ns");
// 		for (std::size_t t = 0; t < m; ++t) {
// 			C_scratch[t] = CoeffCtx::mult(ARev[t], B0[t]);
// 		}
// 	}

// 	// Inverse FFT
// 	{
// 		blst::ScopedTimer timer(timings, "middle.phase2_ifft_ns");
// 		fft<GT, F, CoeffCtx>(C_scratch, inv_roots, true);
// 	}
// 	// (A*B0)[n+k] = D[n-2-k] for k = 0..n-2

// 	// Combine: out[k] = p[k] + q[k]
// 	{
// 		blst::ScopedTimer timer(timings, "middle.phase2_combine_ns");
// 		for (std::size_t k = 0; k < n - 1; ++k) {
// 			out[k] = CoeffCtx::plus(out[k], C_scratch[n - 2 - k]);
// 		}
// 	}
// 	// out[n-1] already has p[n-1], q[n-1] = 0, so no change needed
// }

// // Truncated Inverse NTT (DIF form) - computes only the first L outputs
// // 
// // IMPORTANT: DIF produces output in bit-reversed order. To get natural order output,
// // we need to account for this. However, since we're truncating, we use a different
// // approach: we implement the transpose of the forward DIT FFT.
// //
// // For inverse FFT with natural-order output, we use the same structure as forward
// // but with inverse roots and post-normalization.
// template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
// void fft_inverse_truncated(
// 	std::span<G> X,
// 	std::size_t m,
// 	std::size_t L,
// 	std::span<const F> inv_roots,
// 	const F& m_inv)
// {
// 	BLST_ASSERT(X.size() >= m);
// 	BLST_ASSERT(L <= m);
// 	BLST_ASSERT(inv_roots.size() >= m);
// 	BLST_ASSERT((m & (m - 1)) == 0);  // power of 2

// 	if (L == 0) return;

// 	// For truncated inverse, we use the standard DIT structure (same as forward fft)
// 	// but with inv_roots. The key insight is that the standard fft with inv_roots
// 	// and post-scaling gives the inverse.
// 	//
// 	// We can truncate by not computing butterflies that only affect indices >= L
// 	// in the FINAL output. However, due to bit-reversal at the start, the mapping
// 	// between input indices and output indices is complex.
// 	//
// 	// For correctness, we do the full computation but this can be optimized later.
// 	// The bit-reversal permutation means early truncation is non-trivial.

// 	int n = static_cast<int>(m);

// 	// Bit-reversal permutation (same as forward)
// 	bit_reverse(X.subspan(0, m));

// 	// DIT butterflies with inverse roots
// 	for (int len = 2; len <= n; len <<= 1) {
// 		int half = len / 2;
// 		for (int i = 0; i < n; i += len) {
// 			for (int j = 0; j < half; j++) {
// 				F w = inv_roots[j * (n / len)];

// 				G u = X[i + j];
// 				G t = CoeffCtx::mult(X[i + j + half], w);
// 				X[i + j] = CoeffCtx::plus(u, t);
// 				X[i + j + half] = CoeffCtx::minus(u, t);
// 			}
// 		}
// 	}

// 	// Final normalization - only for prefix we care about
// 	for (std::size_t k = 0; k < L; ++k) {
// 		X[k] = CoeffCtx::mult(X[k], m_inv);
// 	}
// }

// // Full inverse using same structure as forward (for benchmarking comparison)
// template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
// void fft_inverse_full(
// 	std::span<G> X,
// 	std::size_t m,
// 	std::span<const F> inv_roots,
// 	const F& m_inv)
// {
// 	BLST_ASSERT(X.size() >= m);
// 	BLST_ASSERT(inv_roots.size() >= m);
// 	BLST_ASSERT((m & (m - 1)) == 0);

// 	int n = static_cast<int>(m);

// 	bit_reverse(X.subspan(0, m));

// 	for (int len = 2; len <= n; len <<= 1) {
// 		int half = len / 2;
// 		for (int i = 0; i < n; i += len) {
// 			for (int j = 0; j < half; j++) {
// 				F w = inv_roots[j * (n / len)];

// 				G u = X[i + j];
// 				G t = CoeffCtx::mult(X[i + j + half], w);
// 				X[i + j] = CoeffCtx::plus(u, t);
// 				X[i + j + half] = CoeffCtx::minus(u, t);
// 			}
// 		}
// 	}

// 	for (std::size_t k = 0; k < m; ++k) {
// 		X[k] = CoeffCtx::mult(X[k], m_inv);
// 	}
// }

// // Middle product using truncated inverse FFT
// // Note: This version currently does full computation but only normalizes the prefix.
// // True truncation with DIT structure is complex due to bit-reversal permutation.
// template <typename G1, typename G2, typename GT, typename F, typename CoeffCtx>
// void middle_product_fft_truncated(
// 	std::span<G1> A,
// 	std::span<G2> B0,
// 	std::span<G2> B1,
// 	bool skipFftB,
// 	std::size_t n,
// 	std::size_t m,
// 	std::span<GT> out,
// 	std::span<GT> C_scratch,
// 	std::span<const F> roots,
// 	std::span<const F> inv_roots,
// 	const F& omega,
// 	const F& m_inv,
// 	const GT& identity_gt,
// 	const G1& identity_g1,
// 	const G2& identity_g2,
// 	blst::TimingStats* timings = nullptr)
// {
// 	BLST_ASSERT(A.size() >= m);
// 	BLST_ASSERT(B0.size() >= m);
// 	BLST_ASSERT(B1.size() >= m);
// 	BLST_ASSERT(out.size() >= n);
// 	BLST_ASSERT(C_scratch.size() >= m);
// 	BLST_ASSERT(roots.size() >= m);
// 	BLST_ASSERT(inv_roots.size() >= m);
// 	BLST_ASSERT(m >= 2 * n);
// 	BLST_ASSERT((m & (m - 1)) == 0);

// 	// === Phase 1: Compute p = low_n(A * B1) ===

// 	{
// 		blst::ScopedTimer timer(timings, "middle_trunc.phase1_prepare_ns");
// 		for (std::size_t i = n; i < m; ++i) {
// 			A[i] = identity_g1;
// 			if (skipFftB == false) {
// 				B1[i] = identity_g2;
// 			}
// 		}
// 	}

// 	{
// 		blst::ScopedTimer timer(timings, "middle_trunc.phase1_fft_a_ns");
// 		fft<G1, F, CoeffCtx>(A, roots, false);
// 	}
// 	if (skipFftB == false) {
// 		blst::ScopedTimer timer(timings, "middle_trunc.phase1_fft_b1_ns");
// 		fft<G2, F, CoeffCtx>(B1, roots, false);
// 	}

// 	{
// 		blst::ScopedTimer timer(timings, "middle_trunc.phase1_pointwise_ns");
// 		for (std::size_t t = 0; t < m; ++t) {
// 			C_scratch[t] = CoeffCtx::mult(A[t], B1[t]);
// 		}
// 	}

// 	// Truncated inverse FFT: only need first n outputs
// 	{
// 		blst::ScopedTimer timer(timings, "middle_trunc.phase1_ifft_ns");
// 		fft_inverse_truncated<GT, F, CoeffCtx>(C_scratch, m, n, inv_roots, m_inv);
// 	}

// 	{
// 		blst::ScopedTimer timer(timings, "middle_trunc.phase1_copy_ns");
// 		for (std::size_t k = 0; k < n; ++k) {
// 			out[k] = C_scratch[k];
// 		}
// 	}

// 	// === Phase 2: Compute q = high_n(A * B0) via reversal ===

// 	std::vector<G1> ARev(m);
// 	{
// 		blst::ScopedTimer timer(timings, "middle_trunc.phase2_derive_ntt_ns");
// 		derive_ntt_of_reversed<G1, F, CoeffCtx>(
// 			std::span<const G1>(A.data(), m),
// 			ARev,
// 			n, m, omega, identity_g1);
// 	}

// 	if (skipFftB == false)
// 	{
// 		{
// 			blst::ScopedTimer timer(timings, "middle_trunc.phase2_prepare_b0_ns");
// 			for (std::size_t i = 0; i < n / 2; ++i) {
// 				std::swap(B0[i], B0[n - 1 - i]);
// 			}
// 			for (std::size_t i = n; i < m; ++i) {
// 				B0[i] = identity_g2;
// 			}
// 		}
// 		{
// 			blst::ScopedTimer timer(timings, "middle_trunc.phase2_fft_b0_ns");
// 			fft<G2, F, CoeffCtx>(B0, roots, false);
// 		}
// 	}

// 	{
// 		blst::ScopedTimer timer(timings, "middle_trunc.phase2_pointwise_ns");
// 		for (std::size_t t = 0; t < m; ++t) {
// 			C_scratch[t] = CoeffCtx::mult(ARev[t], B0[t]);
// 		}
// 	}

// 	// Truncated inverse FFT: only need first n-1 outputs
// 	{
// 		blst::ScopedTimer timer(timings, "middle_trunc.phase2_ifft_ns");
// 		fft_inverse_truncated<GT, F, CoeffCtx>(C_scratch, m, n - 1, inv_roots, m_inv);
// 	}

// 	{
// 		blst::ScopedTimer timer(timings, "middle_trunc.phase2_combine_ns");
// 		for (std::size_t k = 0; k < n - 1; ++k) {
// 			out[k] = CoeffCtx::plus(out[k], C_scratch[n - 2 - k]);
// 		}
// 	}
// }
