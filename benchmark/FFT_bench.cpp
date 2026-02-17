#include <benchmark/benchmark.h>

#include <vector>
#include <random>
#include <chrono>

#include "blst/fields/Fr.h"
#include "blst/curve/point.h"
#include "blst/curve/GT.h"
#include "blst/fft/primitveRoot.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"

#include "blst/common/defines.h"

#if defined(BLST_CUDA)
#include "blst/fft/ctFft_cuda.h"
#endif

using namespace blst;

namespace
{
	// Pairing function for benchmarks
	static GT pair_points(const G1Point& p, const G2Point& q) {
		G1AffinePoint p_aff;
		G1Point::to_affine(p_aff, p);
		G2AffinePoint q_aff;
		G2Point::to_affine(q_aff, q);

		GT result;
		GT::miller_loop(result, &q_aff, &p_aff);
		GT::final_mult(result, result);
		return result;
	}

	using BenchPairingCoeffCtx = CoeffCtx;

	static std::vector<Fr> make_roots(size_t N)
	{
		Fr prim = primRootOfUnity<Fr>(N);

		std::vector<Fr> roots(N);
		Fr cur = Fr::one();
		for (size_t i = 0; i < N; ++i) {
			roots[i] = cur;
			Fr::mul(cur, cur, prim);
		}
		return roots;
	}

	static std::vector<Fr> make_inv_roots(const std::vector<Fr>& roots)
	{
		std::vector<Fr> inv_roots(roots.size());
		for (size_t i = 0; i < roots.size(); ++i) {
			Fr::reciprocal(inv_roots[i], roots[i]);
		}
		return inv_roots;
	}

	static std::vector<G1Point> make_g1_data(size_t N, uint64_t seed)
	{
		std::mt19937_64 rng(seed);
		std::uniform_int_distribution<uint64_t> dist(1, (uint64_t)-1);

		std::vector<G1Point> data(N);
		for (size_t i = 0; i < N; ++i) {
			Fr s;
			Fr::from_uint64(s, dist(rng));
			Fr s_non_mont;
			Fr::from_mont(s_non_mont, s);
			G1Point::mult(data[i], G1Point::generator(), reinterpret_cast<const byte*>(&s_non_mont), 256);
		}
		return data;
	}

	static std::vector<G2Point> make_g2_data(size_t N, uint64_t seed)
	{
		std::mt19937_64 rng(seed);
		std::uniform_int_distribution<uint64_t> dist(1, (uint64_t)-1);

		std::vector<G2Point> data(N);
		for (size_t i = 0; i < N; ++i) {
			Fr s;
			Fr::from_uint64(s, dist(rng));
			Fr s_non_mont;
			Fr::from_mont(s_non_mont, s);
			G2Point::mult(data[i], G2Point::generator(), reinterpret_cast<const byte*>(&s_non_mont), 256);
		}
		return data;
	}

	static GT pair_gen()
	{
		G1AffinePoint p_aff;
		G2AffinePoint q_aff;
		G1Point::to_affine(p_aff, G1Point::generator());
		G2Point::to_affine(q_aff, G2Point::generator());

		GT r;
		GT::miller_loop(r, &q_aff, &p_aff);
		GT::final_mult(r, r);
		return r;
	}

	static std::vector<GT> make_gt_data(size_t N)
	{
		GT base = pair_gen();

		std::vector<GT> data(N);
		data[0] = GT::identity();

		GT acc = base;
		for (size_t i = 1; i < N; ++i) {
			data[i] = acc;
			GT::add(acc, acc, base);
		}
		return data;
	}

	// ============================================================================
	// Radix-2 FFT benchmarks (existing)
	// ============================================================================

	static void bench_g1_fft_cpu(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<G1Point> original = make_g1_data(N, /*seed*/ 0xC0FFEEULL);
		std::vector<G1Point> data = original;

		for (auto _ : state) {
			data = original;
			benchmark::DoNotOptimize(data.data());
			fft<G1Point, Fr, G1CoeffCtx>(data, roots, /*inverse*/ false);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["FFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	//static void bench_g1_fft_dif_cpu(benchmark::State& state)
	//{
	//	const size_t N = static_cast<size_t>(state.range(0));

	//	const std::vector<Fr> roots = make_roots(N);
	//	const std::vector<G1Point> original = make_g1_data(N, /*seed*/ 0xC0FFEEULL);
	//	std::vector<G1Point> data = original;

	//	for (auto _ : state) {
	//		data = original;
	//		benchmark::DoNotOptimize(data.data());
	//		fft_dif<G1Point, Fr, G1CoeffCtx>(data, roots, /*inverse*/ false);
	//		benchmark::ClobberMemory();
	//	}

	//	state.counters["N"] = static_cast<double>(N);
	//	state.counters["FFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	//}

	// ============================================================================
	// Radix-4 FFT benchmarks
	// ============================================================================

	static void bench_g1_fft_radix4_cpu(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<G1Point> original = make_g1_data(N, /*seed*/ 0xC0FFEEULL);
		std::vector<G1Point> data = original;

		for (auto _ : state) {
			data = original;
			benchmark::DoNotOptimize(data.data());
			fft_radix4<G1Point, Fr, G1CoeffCtx>(data, roots, /*inverse*/ false);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["FFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	static void bench_gt_ifft_radix4_cpu(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);

		const std::vector<GT> time_domain = make_gt_data(N);
		std::vector<GT> freq_domain = time_domain;
		fft<GT, Fr, GTCoeffCtx>(freq_domain, roots, /*inverse*/ false);

		std::vector<GT> data = freq_domain;

		for (auto _ : state) {
			data = freq_domain;
			benchmark::DoNotOptimize(data.data());
			fft_radix4<GT, Fr, GTCoeffCtx>(data, inv_roots, /*inverse*/ true);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["iFFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	// ============================================================================
	// Radix-8 FFT benchmarks
	// ============================================================================

	static void bench_g1_fft_radix8_cpu(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<G1Point> original = make_g1_data(N, /*seed*/ 0xC0FFEEULL);
		std::vector<G1Point> data = original;

		for (auto _ : state) {
			data = original;
			benchmark::DoNotOptimize(data.data());
			fft_radix8<G1Point, Fr, G1CoeffCtx>(data, roots, /*inverse*/ false);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["FFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	static void bench_gt_ifft_radix8_cpu(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);

		const std::vector<GT> time_domain = make_gt_data(N);
		std::vector<GT> freq_domain = time_domain;
		fft<GT, Fr, GTCoeffCtx>(freq_domain, roots, /*inverse*/ false);

		std::vector<GT> data = freq_domain;

		for (auto _ : state) {
			data = freq_domain;
			benchmark::DoNotOptimize(data.data());
			fft_radix8<GT, Fr, GTCoeffCtx>(data, inv_roots, /*inverse*/ true);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["iFFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}


#if defined(BLST_CUDA)
	static void bench_g1_fft_gpu(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<G1Point> original = make_g1_data(N, /*seed*/ 0xC0FFEEULL);
		std::vector<G1Point> data = original;

		for (auto _ : state) {
			data = original;
			benchmark::DoNotOptimize(data.data());
			fft_cuda(data, roots, /*inverse*/ false);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["FFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}
#endif

	static void bench_gt_ifft_cpu(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);

		// Put input in "frequency domain" first so we can time just the inverse FFT.
		const std::vector<GT> time_domain = make_gt_data(N);
		std::vector<GT> freq_domain = time_domain;
		fft<GT, Fr, GTCoeffCtx>(freq_domain, roots, /*inverse*/ false);

		std::vector<GT> data = freq_domain;

		for (auto _ : state) {
			data = freq_domain;
			benchmark::DoNotOptimize(data.data());
			fft<GT, Fr, GTCoeffCtx>(data, inv_roots, /*inverse*/ false);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["iFFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

#if defined(BLST_CUDA)
	static void bench_gt_ifft_gpu(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);

		const std::vector<GT> time_domain = make_gt_data(N);
		std::vector<GT> freq_domain = time_domain;
		fft_cuda(freq_domain, roots, /*inverse*/ false);

		std::vector<GT> data = freq_domain;

		for (auto _ : state) {
			data = freq_domain;
			benchmark::DoNotOptimize(data.data());
			fft_cuda(data, inv_roots, /*inverse*/ true);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["iFFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}
#endif



	// ============================================================================
	// Stockham FFT benchmarks
	// ============================================================================

	static void bench_g1_fft_stockham_cpu(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<G1Point> original = make_g1_data(N, /*seed*/ 0xC0FFEEULL);
		std::vector<G1Point> out0(N);
		std::vector<G1Point> out1(N);

		for (auto _ : state) {
			benchmark::DoNotOptimize(original.data());
			benchmark::DoNotOptimize(out0.data());
			benchmark::DoNotOptimize(out1.data());
			auto result = fft_stockham<G1Point, Fr, G1CoeffCtx>(
				std::span<const G1Point>(original),
				std::span<G1Point>(out0),
				std::span<G1Point>(out1),
				std::span<const Fr>(roots),
				/*invert*/ false);
			benchmark::DoNotOptimize(result.data());
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["FFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	static void bench_gt_ifft_stockham_cpu(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);

		const std::vector<GT> time_domain = make_gt_data(N);
		std::vector<GT> freq_domain = time_domain;
		fft<GT, Fr, GTCoeffCtx>(freq_domain, roots, /*inverse*/ false);

		std::vector<GT> out0(N);
		std::vector<GT> out1(N);

		for (auto _ : state) {
			benchmark::DoNotOptimize(freq_domain.data());
			benchmark::DoNotOptimize(out0.data());
			benchmark::DoNotOptimize(out1.data());
			auto result = fft_stockham<GT, Fr, GTCoeffCtx>(
				std::span<const GT>(freq_domain),
				std::span<GT>(out0),
				std::span<GT>(out1),
				std::span<const Fr>(inv_roots),
				/*invert*/ false);
			benchmark::DoNotOptimize(result.data());
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["iFFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}
































	// Benchmark middle_product_basic: computes full convolution of G1 * G2 -> GT
	static void bench_middle_product_basic(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const size_t N = 1 << log2ceil(n * 3); // FFT size for convolution of three size-n polynomials

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);

		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(n, /*seed*/ 0xCAFEBABEULL);

		// Pre-allocate buffers
		std::vector<G1Point> A(N, G1Point::infinity());
		std::vector<G2Point> B(N, G2Point::infinity());
		std::vector<GT> C(N, GT::identity());

		for (auto _ : state) {
			// Reset inputs
			for (size_t i = 0; i < n; ++i) {
				A[i] = A_orig[i];
				B[i] = B_orig[i];
			}
			for (size_t i = n; i < N; ++i) {
				A[i] = G1Point::infinity();
				B[i] = G2Point::infinity();
			}

			benchmark::DoNotOptimize(A.data());
			benchmark::DoNotOptimize(B.data());
			benchmark::DoNotOptimize(C.data());

			poly_mult_fft<G1Point, G2Point, GT, Fr, BenchPairingCoeffCtx>(
				A, B, /*skipBFft*/ true, C, roots, inv_roots);

			benchmark::ClobberMemory();
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["N"] = static_cast<double>(N);
		state.counters["ops/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	// Detailed benchmark breaking down bench_middle_product_basic into phases
	static void bench_middle_product_basic_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const size_t N = 1 << log2ceil(n * 3); // FFT size for convolution of three size-n polynomials

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);

		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(n, /*seed*/ 0xCAFEBABEULL);

		// Pre-allocate buffers
		std::vector<G1Point> A(N, G1Point::infinity());
		std::vector<G2Point> B(N, G2Point::infinity());
		std::vector<GT> C(N, GT::identity());

		// Accumulators for phase timings (in nanoseconds)
		double time_fft_g1 = 0;
		double time_fft_g2 = 0;
		double time_pairing = 0;
		double time_ifft_gt = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			state.PauseTiming();
			// Reset inputs
			for (size_t i = 0; i < n; ++i) {
				A[i] = A_orig[i];
				B[i] = B_orig[i];
			}
			for (size_t i = n; i < N; ++i) {
				A[i] = G1Point::infinity();
				B[i] = G2Point::infinity();
			}
			state.ResumeTiming();

			auto t0 = std::chrono::high_resolution_clock::now();

			// Forward FFT on G1
			fft<G1Point, Fr, G1CoeffCtx>(A, roots, false);
			auto t1 = std::chrono::high_resolution_clock::now();

			// Forward FFT on G2 (skipBFft=true in original, so this simulates skipBFft=false)
			//fft<G2Point, Fr, G2CoeffCtx>(B, roots, false);
			auto t2 = std::chrono::high_resolution_clock::now();

			// Pointwise pairing: e(A[i], B[i]) -> GT
			for (size_t i = 0; i < N; ++i) {
				C[i] = pair_points(A[i], B[i]);
			}
			auto t3 = std::chrono::high_resolution_clock::now();

			// Inverse FFT on GT
			fft<GT, Fr, GTCoeffCtx>(C, inv_roots, true);
			auto t4 = std::chrono::high_resolution_clock::now();

			benchmark::ClobberMemory();

			// Accumulate phase times
			time_fft_g1 += std::chrono::duration<double, std::nano>(t1 - t0).count();
			time_fft_g2 += std::chrono::duration<double, std::nano>(t2 - t1).count();
			time_pairing += std::chrono::duration<double, std::nano>(t3 - t2).count();
			time_ifft_gt += std::chrono::duration<double, std::nano>(t4 - t3).count();
			iterations++;
		}

		// Report phase timings
		double total = time_fft_g1 + time_fft_g2 + time_pairing + time_ifft_gt;

		state.counters["n"] = static_cast<double>(n);
		state.counters["N"] = static_cast<double>(N);

		// Absolute times in milliseconds
		state.counters["fft_g1_ms"] = benchmark::Counter(time_fft_g1 / 1e6 / iterations);
		state.counters["fft_g2_ms"] = benchmark::Counter(time_fft_g2 / 1e6 / iterations);
		state.counters["pairing_ms"] = benchmark::Counter(time_pairing / 1e6 / iterations);
		state.counters["ifft_gt_ms"] = benchmark::Counter(time_ifft_gt / 1e6 / iterations);

		//// Percentages
		//state.counters["fft_g1_pct"] = benchmark::Counter(100.0 * time_fft_g1 / total);
		//state.counters["fft_g2_pct"] = benchmark::Counter(100.0 * time_fft_g2 / total);
		//state.counters["pairing_pct"] = benchmark::Counter(100.0 * time_pairing / total);
		//state.counters["ifft_gt_pct"] = benchmark::Counter(100.0 * time_ifft_gt / total);
	}

	// Benchmark middle_product_fft: computes middle n coefficients of (size n) * (size 2n) convolution
	static void bench_middle_product_fft(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		size_t m = 1;
		while (m < 2 * n) m <<= 1; // FFT size

		const std::vector<Fr> roots = make_roots(m);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);
		const Fr prim = primRootOfUnity<Fr>(m);

		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);

		// Pre-allocate buffers
		std::vector<G1Point> A(m, G1Point::infinity());
		std::vector<G2Point> B0(m, G2Point::infinity());
		std::vector<G2Point> B1(m, G2Point::infinity());
		std::vector<GT> out(n, GT::identity());
		std::vector<GT> C_scratch(m, GT::identity());

		for (auto _ : state) {
			// Reset inputs
			for (size_t i = 0; i < n; ++i) {
				A[i] = A_orig[i];
				B0[i] = B_orig[i];
				B1[i] = B_orig[n + i];
			}
			for (size_t i = n; i < m; ++i) {
				A[i] = G1Point::infinity();
				B0[i] = G2Point::infinity();
				B1[i] = G2Point::infinity();
			}

			benchmark::DoNotOptimize(A.data());
			benchmark::DoNotOptimize(out.data());

			middle_product_fft<G1Point, G2Point, GT, Fr, BenchPairingCoeffCtx>(
				std::span<G1Point>(A), std::span<G2Point>(B0), std::span<G2Point>(B1),
				true,
				n, m,
				std::span<GT>(out), std::span<GT>(C_scratch),
				std::span<const Fr>(roots), std::span<const Fr>(inv_roots),
				prim,
				GT::identity(), G1Point::infinity(), G2Point::infinity());

			benchmark::ClobberMemory();
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(m);
		state.counters["ops/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	// Benchmark middle_product_fft_truncated: uses truncated inverse FFT
	static void bench_middle_product_fft_truncated(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		size_t m = 1;
		while (m < 2 * n) m <<= 1; // FFT size

		const std::vector<Fr> roots = make_roots(m);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);
		const Fr prim = primRootOfUnity<Fr>(m);
		Fr m_inv = Fr::from_uint64(m);
		Fr::reciprocal(m_inv, m_inv);

		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);

		// Pre-allocate buffers
		std::vector<G1Point> A(m, G1Point::infinity());
		std::vector<G2Point> B0(m, G2Point::infinity());
		std::vector<G2Point> B1(m, G2Point::infinity());
		std::vector<GT> out(n, GT::identity());
		std::vector<GT> C_scratch(m, GT::identity());

		for (auto _ : state) {
			// Reset inputs
			for (size_t i = 0; i < n; ++i) {
				A[i] = A_orig[i];
				B0[i] = B_orig[i];
				B1[i] = B_orig[n + i];
			}
			for (size_t i = n; i < m; ++i) {
				A[i] = G1Point::infinity();
				B0[i] = G2Point::infinity();
				B1[i] = G2Point::infinity();
			}

			benchmark::DoNotOptimize(A.data());
			benchmark::DoNotOptimize(out.data());

			middle_product_fft_truncated<G1Point, G2Point, GT, Fr, BenchPairingCoeffCtx>(
				std::span<G1Point>(A), std::span<G2Point>(B0), std::span<G2Point>(B1),
				true,
				n, m,
				std::span<GT>(out), std::span<GT>(C_scratch),
				std::span<const Fr>(roots), std::span<const Fr>(inv_roots),
				prim, m_inv,
				GT::identity(), G1Point::infinity(), G2Point::infinity());

			benchmark::ClobberMemory();
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(m);
		state.counters["ops/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	// Benchmark truncated inverse FFT in isolation
	static void bench_gt_ifft_truncated(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));
		const size_t L = N / 2; // Truncate to half

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);
		Fr m_inv = Fr::from_uint64(N);
		Fr::reciprocal(m_inv, m_inv);

		const std::vector<GT> time_domain = make_gt_data(N);
		std::vector<GT> freq_domain = time_domain;
		fft<GT, Fr, GTCoeffCtx>(freq_domain, roots, /*inverse*/ false);

		std::vector<GT> data = freq_domain;

		for (auto _ : state) {
			data = freq_domain;
			benchmark::DoNotOptimize(data.data());
			fft_inverse_truncated<GT, Fr, GTCoeffCtx>(
				std::span<GT>(data), N, L,
				std::span<const Fr>(inv_roots), m_inv);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["L"] = static_cast<double>(L);
		state.counters["iFFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}
}


// Detailed benchmark breaking down middle_product_fft into phases
static void bench_middle_product_fft_detailed(benchmark::State& state)
{
	const size_t n = static_cast<size_t>(state.range(0));
	size_t m = 1;
	while (m < 2 * n) m <<= 1;

	const std::vector<Fr> roots = make_roots(m);
	const std::vector<Fr> inv_roots = make_inv_roots(roots);
	const Fr prim = primRootOfUnity<Fr>(m);

	const std::vector<G1Point> A_orig = make_g1_data(n, 0xDEADBEEFULL);
	const std::vector<G2Point> B_orig = make_g2_data(2 * n, 0xCAFEBABEULL);

	// Pre-allocate all buffers
	std::vector<G1Point> AA(m * 2, G1Point::infinity());
	std::vector<G2Point> BB(m * 2, G2Point::infinity());
	//std::vector<G2Point> B0(m, G2Point::infinity());
	//std::vector<G2Point> B1(m, G2Point::infinity());
	std::vector<GT> out(n, GT::identity());
	std::vector<GT> C_scratch(m * 2, GT::identity());
	//std::vector<GT> C_scratch0(m, GT::identity());
	//std::vector<GT> C_scratch1(m, GT::identity());
	//std::vector<G1Point> ARev(m, G1Point::infinity());

	std::span<G1Point> A(AA.data(), m);
	std::span<G1Point> ARev(AA.data() + m, m);
	std::span<G2Point> B1(BB.data(), m);
	std::span<G2Point> B0(BB.data() + m, m);
	std::span<GT> C_scratch0(C_scratch.data(), m);
	std::span<GT> C_scratch1(C_scratch.data() + m, m);


	// Accumulators for phase timings (in nanoseconds)
	double time_fft_g1 = 0;
	double time_fft_g2 = 0;
	double time_pairing = 0;
	double time_ifft_gt = 0;
	double time_derive_rev = 0;
	double time_combine = 0;
	int64_t iterations = 0;

	for (auto _ : state) {
		state.PauseTiming();
		// Reset inputs
		for (size_t i = 0; i < n; ++i) {
			A[i] = A_orig[i];
			B0[i] = B_orig[i];
			B1[i] = B_orig[n + i];
		}
		for (size_t i = n; i < m; ++i) {
			A[i] = G1Point::infinity();
			B0[i] = G2Point::infinity();
			B1[i] = G2Point::infinity();
		}
		state.ResumeTiming();

		auto t0 = std::chrono::high_resolution_clock::now();

		// === Phase 1: FFT on G1 (A) ===
		fft<G1Point, Fr, G1CoeffCtx>(A, roots, false);
		auto t1 = std::chrono::high_resolution_clock::now();

		// === Phase 1: FFT on G2 (B1) ===
		//fft<G2Point, Fr, G2CoeffCtx>(B1, roots, false);
		//auto t2 = std::chrono::high_resolution_clock::now();

		// === Phase 1: Pointwise pairing ===
		//for (size_t t = 0; t < m; ++t) {
		//}
		auto t3 = std::chrono::high_resolution_clock::now();

		// === Phase 1: Inverse FFT on GT ===
		auto t4 = std::chrono::high_resolution_clock::now();

		// Save p[k] into out
		//for (size_t k = 0; k < n; ++k) {
		//	out[k] = C_scratch[k];
		//}

		// === Phase 2: Derive NTT of reversed A ===
		derive_ntt_of_reversed<G1Point, Fr, BenchPairingCoeffCtx>(
			A, ARev,
			n, m, prim, G1Point::infinity());
		auto t5 = std::chrono::high_resolution_clock::now();

		//// Reverse B0 in-place
		//for (size_t i = 0; i < n / 2; ++i) {
		//	std::swap(B0[i], B0[n - 1 - i]);
		//}

		//// === Phase 2: FFT on G2 (rev B0) ===
		//fft<G2Point, Fr, G2CoeffCtx>(B0, roots, false);
		//auto t6 = std::chrono::high_resolution_clock::now();

		// === Phase 2: Pointwise pairing ===
		//for (size_t t = 0; t < m; ++t) 
		//	C_scratch1[t] = pair_points(ARev[t], B0[t]);
		//for (size_t t = 0; t < m; ++t)
		//	C_scratch0[t] = pair_points(A[t], B1[t]);
		auto m2 = m * 2;
		for (size_t t = 0; t < m2; ++t)
			C_scratch[t] = pair_points(AA[t], BB[t]);

		auto t7 = std::chrono::high_resolution_clock::now();

		// === Phase 2: Inverse FFT on GT ===
		fft<GT, Fr, GTCoeffCtx>(C_scratch0, inv_roots, true);
		fft<GT, Fr, GTCoeffCtx>(C_scratch1, inv_roots, true);
		auto t8 = std::chrono::high_resolution_clock::now();

		// === Combine results ===
		for (size_t k = 0; k < n - 1; ++k) {
			out[k] = GTCoeffCtx::plus(C_scratch0[k], C_scratch1[n - 2 - k]);
		}
		auto t9 = std::chrono::high_resolution_clock::now();

		benchmark::ClobberMemory();

		// Accumulate phase times
		time_fft_g1 += std::chrono::duration<double, std::nano>(t1 - t0).count();
		//time_fft_g2 += std::chrono::duration<double, std::nano>(t2 - t1).count();
		//time_fft_g2 += std::chrono::duration<double, std::nano>(t6 - t5).count(); // second G2 FFT
		time_pairing += std::chrono::duration<double, std::nano>(t3 - t1).count();
		time_pairing += std::chrono::duration<double, std::nano>(t7 - t5).count(); // second pairing
		time_ifft_gt += std::chrono::duration<double, std::nano>(t4 - t3).count();
		time_ifft_gt += std::chrono::duration<double, std::nano>(t8 - t7).count(); // second iFFT
		time_derive_rev += std::chrono::duration<double, std::nano>(t5 - t4).count();
		time_combine += std::chrono::duration<double, std::nano>(t9 - t8).count();
		iterations++;
	}

	// Report phase timings as percentages and absolute times
	double total = time_fft_g1 + time_fft_g2 + time_pairing + time_ifft_gt + time_derive_rev + time_combine;

	state.counters["n"] = static_cast<double>(n);
	state.counters["m"] = static_cast<double>(m);
	state.counters["fft_g1_ms"] = benchmark::Counter(time_fft_g1 / 1e6 / iterations);
	state.counters["fft_g2_ms"] = benchmark::Counter(time_fft_g2 / 1e6 / iterations);
	state.counters["pairing_ms"] = benchmark::Counter(time_pairing / 1e6 / iterations);
	state.counters["ifft_gt_ms"] = benchmark::Counter(time_ifft_gt / 1e6 / iterations);
	state.counters["derive_rev_ms"] = benchmark::Counter(time_derive_rev / 1e6 / iterations);
	state.counters["combine_ms"] = benchmark::Counter(time_combine / 1e6 / iterations);
	//state.counters["fft_g1_pct"] = benchmark::Counter(100.0 * time_fft_g1 / total);
	//state.counters["fft_g2_pct"] = benchmark::Counter(100.0 * time_fft_g2 / total);
	//state.counters["pairing_pct"] = benchmark::Counter(100.0 * time_pairing / total);
	//state.counters["ifft_gt_pct"] = benchmark::Counter(100.0 * time_ifft_gt / total);
	//state.counters["derive_rev_pct"] = benchmark::Counter(100.0 * time_derive_rev / total);
	//state.counters["combine_pct"] = benchmark::Counter(100.0 * time_combine / total);
}

// Keep sizes as powers of two (roots exist up to 2^32 for Fr).
#define FFT_SIZES \
->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)->Arg(8192)->Arg(16384)

// Smaller sizes for expensive operations (middle product involves pairings)
#define MIDDLE_PROD_SIZES \
->Arg(96)->Arg(128)->Arg(160)->Arg(192)

// ============================================================================
// Radix-2 FFT benchmarks (baseline)
// ============================================================================

BENCHMARK(bench_g1_fft_cpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();


//BENCHMARK(bench_g1_fft_dif_cpu)
//FFT_SIZES
//->Unit(benchmark::kMillisecond)
//->UseRealTime();

BENCHMARK(bench_gt_ifft_cpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();


BENCHMARK(bench_g1_fft_stockham_cpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_gt_ifft_stockham_cpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

// ============================================================================
// Radix-4 FFT benchmarks
// ============================================================================

BENCHMARK(bench_g1_fft_radix4_cpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_gt_ifft_radix4_cpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

// ============================================================================
// Radix-8 FFT benchmarks
// ============================================================================

BENCHMARK(bench_g1_fft_radix8_cpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_gt_ifft_radix8_cpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

// ============================================================================
// Adaptive FFT benchmarks
// ============================================================================

//BENCHMARK(bench_g1_fft_adaptive_cpu)
//FFT_SIZES
//->Unit(benchmark::kMillisecond)
//->UseRealTime();
//
//BENCHMARK(bench_gt_ifft_adaptive_cpu)
//FFT_SIZES
//->Unit(benchmark::kMillisecond)
//->UseRealTime();

// ============================================================================
// Comparison benchmarks (all radix variants with same N)
// ============================================================================

//BENCHMARK(bench_g1_fft_comparison)
//->Args({ 1024, 2 })->Args({ 1024, 4 })->Args({ 1024, 8 })
//->Args({ 4096, 2 })->Args({ 4096, 4 })->Args({ 4096, 8 })
//->Args({ 16384, 2 })->Args({ 16384, 4 })->Args({ 16384, 8 })
//->Unit(benchmark::kMillisecond)
//->UseRealTime();
//
//BENCHMARK(bench_gt_ifft_comparison)
//->Args({ 1024, 2 })->Args({ 1024, 4 })->Args({ 1024, 8 })
//->Args({ 4096, 2 })->Args({ 4096, 4 })->Args({ 4096, 8 })
//->Args({ 16384, 2 })->Args({ 16384, 4 })->Args({ 16384, 8 })
//->Unit(benchmark::kMillisecond)
//->UseRealTime();

// ============================================================================
// Truncated and other FFT benchmarks
// ============================================================================

BENCHMARK(bench_gt_ifft_truncated)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

#if defined(BLST_CUDA)
BENCHMARK(bench_g1_fft_gpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_gt_ifft_gpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
#endif

// ============================================================================
// Middle product benchmarks
// ============================================================================

BENCHMARK(bench_middle_product_basic)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_fft)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_fft_truncated)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_fft_detailed)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_basic_detailed)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
//BENCHMARK_MAIN();
