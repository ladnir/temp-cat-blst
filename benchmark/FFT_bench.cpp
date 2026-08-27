#include <benchmark/benchmark.h>

#include <vector>
#include <random>
#include <chrono>

#include "blst/common/defines.h"
#include "blst/fields/Fr.h"
#include "blst/curve/point.h"
#include "blst/curve/GT.h"
#include "blst/fft/primitive_root.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/middle_product.h"
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "blst/avx/fft/gt_context.h"
#endif
#include "blst/common/timer.h"

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

		MillerResult miller;
		GT result;
		GT::miller_loop(miller, &q_aff, &p_aff);
		GT::final_exp(result, miller);
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

		MillerResult miller;
		GT r;
		GT::miller_loop(miller, &q_aff, &p_aff);
		GT::final_exp(r, miller);
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

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
	static std::vector<blst::avx::AvxGt52> make_avx_gt_data(const std::vector<GT>& scalar_data)
	{
		return blst::avx::scalar_gt_vector_to_avx52(std::span<const GT>(scalar_data));
	}
#endif

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
			fft<GT, Fr, GTCoeffCtx>(data, inv_roots, /*inverse*/ true);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["iFFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
	static void bench_gt_ifft_avx52_precomputed_cpu(benchmark::State& state)
	{
		const size_t N = static_cast<size_t>(state.range(0));

		const std::vector<Fr> roots = make_roots(N);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);

		const std::vector<GT> time_domain = make_gt_data(N);
		std::vector<GT> freq_domain_scalar = time_domain;
		fft<GT, Fr, GTCoeffCtx>(freq_domain_scalar, roots, /*inverse*/ false);

		const std::vector<blst::avx::AvxGt52> freq_domain = make_avx_gt_data(freq_domain_scalar);
		std::vector<blst::avx::AvxGt52> data = freq_domain;

		const auto precompute_begin = std::chrono::steady_clock::now();
		blst::avx::AvxGt52PreparedRoots prepared_inv_roots{std::span<const Fr>(inv_roots)};
		prepared_inv_roots.set_scale_from_size(N);
		const auto precompute_end = std::chrono::steady_clock::now();
		const double precompute_ms = std::chrono::duration<double, std::milli>(precompute_end - precompute_begin).count();

		blst::avx::AvxGt52CoeffCtx::ScopedRoots scoped_roots(prepared_inv_roots);

		for (auto _ : state) {
			data = freq_domain;
			benchmark::DoNotOptimize(data.data());
			fft<blst::avx::AvxGt52, Fr, blst::avx::AvxGt52CoeffCtx>(data, inv_roots, /*inverse*/ true);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(N);
		state.counters["root_glv4_precompute_ms"] = benchmark::Counter(precompute_ms);
		state.counters["iFFTs/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}
#endif

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
































	// Benchmark middle_product_basic: computes full convolution of G1 * G2 -> GT
	static void bench_middle_product_basic(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);
		std::vector<GT> out(n, GT::identity());
		const size_t m = middle_product_basic_fft_size(n);

		for (auto _ : state) {
			benchmark::DoNotOptimize(out.data());
			middle_product_basic(
				std::span<const G1Point>(A_orig),
				std::span<const G2Point>(B_orig),
				std::span<GT>(out));

			benchmark::ClobberMemory();
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(m);
		state.counters["ops/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	// Detailed benchmark breaking down bench_middle_product_basic into phases
	static void bench_middle_product_basic_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);
		std::vector<GT> out(n, GT::identity());
		const size_t m = middle_product_basic_fft_size(n);

		// Accumulators for phase timings (in nanoseconds)
		double time_fft_g1 = 0;
		double time_fft_g2 = 0;
		double time_pairing = 0;
		double time_ifft_gt = 0;
		double time_roots = 0;
		double time_extract = 0;
		double fft_a_total = 0;
		double fft_a_exec = 0;
		double fft_a_skip = 0;
		double fft_b_total = 0;
		double fft_b_exec = 0;
		double fft_b_skip = 0;
		double online_total = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			TimingStats iter_stats;
			middle_product_basic(
				std::span<const G1Point>(A_orig),
				std::span<const G2Point>(B_orig),
				std::span<GT>(out),
				&iter_stats);
			benchmark::ClobberMemory();

			time_fft_g1 += static_cast<double>(iter_stats.get_ns("poly_mult.fft_a_ns"));
			time_fft_g2 += static_cast<double>(iter_stats.get_ns("poly_mult.fft_b_ns"));
			time_pairing += static_cast<double>(iter_stats.get_ns("poly_mult.pointwise_ns"));
			time_ifft_gt += static_cast<double>(iter_stats.get_ns("poly_mult.ifft_ns"));
			time_roots += static_cast<double>(iter_stats.get_ns("middle_product_basic.roots_ns"));
			time_roots += static_cast<double>(iter_stats.get_ns("middle_product_basic.inv_roots_ns"));
			time_extract += static_cast<double>(iter_stats.get_ns("middle_product_basic.extract_ns"));
			fft_a_total += static_cast<double>(iter_stats.get_ns("poly_mult.fft_a.butterflies_total"));
			fft_a_exec += static_cast<double>(iter_stats.get_ns("poly_mult.fft_a.butterflies_executed"));
			fft_a_skip += static_cast<double>(iter_stats.get_ns("poly_mult.fft_a.butterflies_skipped"));
			fft_b_total += static_cast<double>(iter_stats.get_ns("poly_mult.fft_b.butterflies_total"));
			fft_b_exec += static_cast<double>(iter_stats.get_ns("poly_mult.fft_b.butterflies_executed"));
			fft_b_skip += static_cast<double>(iter_stats.get_ns("poly_mult.fft_b.butterflies_skipped"));
			online_total += static_cast<double>(iter_stats.get_ns("middle_product_basic.total_ns"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(m);

		// Absolute times in milliseconds
		state.counters["fft_g1_ms"] = benchmark::Counter(time_fft_g1 / 1e6 / iterations);
		state.counters["fft_g2_ms"] = benchmark::Counter(time_fft_g2 / 1e6 / iterations);
		state.counters["pairing_ms"] = benchmark::Counter(time_pairing / 1e6 / iterations);
		state.counters["ifft_gt_ms"] = benchmark::Counter(time_ifft_gt / 1e6 / iterations);
		state.counters["roots_ms"] = benchmark::Counter(time_roots / 1e6 / iterations);
		state.counters["extract_ms"] = benchmark::Counter(time_extract / 1e6 / iterations);
		state.counters["online_total_ms"] = benchmark::Counter(online_total / 1e6 / iterations);
		state.counters["online_compute_ms"] = benchmark::Counter((online_total - time_roots) / 1e6 / iterations);
		state.counters["fft_a_bfly_total"] = benchmark::Counter(fft_a_total / iterations);
		state.counters["fft_a_bfly_exec"] = benchmark::Counter(fft_a_exec / iterations);
		state.counters["fft_a_bfly_skip"] = benchmark::Counter(fft_a_skip / iterations);
		state.counters["fft_a_skip_ratio"] = benchmark::Counter(
			(fft_a_total > 0.0) ? (fft_a_skip / fft_a_total) : 0.0);
		state.counters["fft_b_bfly_total"] = benchmark::Counter(fft_b_total / iterations);
		state.counters["fft_b_bfly_exec"] = benchmark::Counter(fft_b_exec / iterations);
		state.counters["fft_b_bfly_skip"] = benchmark::Counter(fft_b_skip / iterations);
		state.counters["fft_b_skip_ratio"] = benchmark::Counter(
			(fft_b_total > 0.0) ? (fft_b_skip / fft_b_total) : 0.0);
	}

	// Detailed benchmark for Basic FFT path with NO truncation/pruning.
	// Uses full-size FFT for A and B (xn = fft_size) to provide a baseline.
	static void bench_middle_product_basic_no_trunc_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
			const size_t m = middle_product_basic_fft_size(n);
		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);
		const std::vector<Fr> roots = make_roots(m);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);

		// Offline precompute B with full FFT (no truncation).
		TimingStats pre_stats;
		std::vector<G2Point> B_fft(m, G2Point::infinity());
		{
			ScopedTimer timer(&pre_stats, "precompute_no_trunc.roots_ns");
		}
		for (size_t i = 0; i < 2 * n; ++i) {
			B_fft[i] = B_orig[i];
		}
		{
			ScopedTimer timer(&pre_stats, "precompute_no_trunc.fft_b_ns");
			fft<G2Point, Fr, G2CoeffCtx>(std::span<G2Point>(B_fft), std::span<const Fr>(roots), false);
		}

		std::vector<GT> out(n, GT::identity());

		double time_fft_g1 = 0;
		double time_pairing = 0;
		double time_ifft_gt = 0;
		double time_roots = 0;
		double time_extract = 0;
		double fft_a_total = 0;
		double fft_a_exec = 0;
		double fft_a_skip = 0;
		double online_total = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			std::vector<G1Point> A_fft(m, G1Point::infinity());
			for (size_t i = 0; i < n; ++i) {
				A_fft[i] = A_orig[i];
			}
			std::vector<G2Point> B_fft_copy = B_fft;
			std::vector<GT> C_fft(m, GT::identity());

			TimingStats iter_stats;
			{
				ScopedTimer timer(&iter_stats, "basic_no_trunc.roots_ns");
			}
			poly_mult_fft<G1Point, G2Point, GT, Fr, CoeffCtx>(
				std::span<G1Point>(A_fft),
				std::span<G2Point>(B_fft_copy),
				true,
				std::span<GT>(C_fft),
				std::span<const Fr>(roots),
				std::span<const Fr>(inv_roots),
				&iter_stats);

			{
				ScopedTimer timer(&iter_stats, "basic_no_trunc.extract_ns");
				for (size_t i = 0; i < n; ++i) {
					out[i] = C_fft[n + i];
				}
			}
			benchmark::ClobberMemory();

			time_fft_g1 += static_cast<double>(iter_stats.get_ns("poly_mult.fft_a_ns"));
			time_pairing += static_cast<double>(iter_stats.get_ns("poly_mult.pointwise_ns"));
			time_ifft_gt += static_cast<double>(iter_stats.get_ns("poly_mult.ifft_ns"));
			time_roots += static_cast<double>(iter_stats.get_ns("basic_no_trunc.roots_ns"));
			time_extract += static_cast<double>(iter_stats.get_ns("basic_no_trunc.extract_ns"));
			fft_a_total += static_cast<double>(iter_stats.get_ns("poly_mult.fft_a.butterflies_total"));
			fft_a_exec += static_cast<double>(iter_stats.get_ns("poly_mult.fft_a.butterflies_executed"));
			fft_a_skip += static_cast<double>(iter_stats.get_ns("poly_mult.fft_a.butterflies_skipped"));
			online_total += static_cast<double>(iter_stats.get_ns("poly_mult.fft_a_ns"));
			online_total += static_cast<double>(iter_stats.get_ns("poly_mult.pointwise_ns"));
			online_total += static_cast<double>(iter_stats.get_ns("poly_mult.ifft_ns"));
			online_total += static_cast<double>(iter_stats.get_ns("basic_no_trunc.extract_ns"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(m);
		state.counters["fft_g1_ms"] = benchmark::Counter(time_fft_g1 / 1e6 / iterations);
		state.counters["pairing_ms"] = benchmark::Counter(time_pairing / 1e6 / iterations);
		state.counters["ifft_gt_ms"] = benchmark::Counter(time_ifft_gt / 1e6 / iterations);
		state.counters["roots_ms"] = benchmark::Counter(time_roots / 1e6 / iterations);
		state.counters["extract_ms"] = benchmark::Counter(time_extract / 1e6 / iterations);
		state.counters["online_total_ms"] = benchmark::Counter(online_total / 1e6 / iterations);
		state.counters["online_compute_ms"] = benchmark::Counter((online_total - time_roots) / 1e6 / iterations);
		state.counters["fft_a_bfly_total"] = benchmark::Counter(fft_a_total / iterations);
		state.counters["fft_a_bfly_exec"] = benchmark::Counter(fft_a_exec / iterations);
		state.counters["fft_a_bfly_skip"] = benchmark::Counter(fft_a_skip / iterations);
		state.counters["fft_a_skip_ratio"] = benchmark::Counter(
			(fft_a_total > 0.0) ? (fft_a_skip / fft_a_total) : 0.0);
		state.counters["offline_precompute_once_ms"] = benchmark::Counter(
			(static_cast<double>(pre_stats.get_ns("precompute_no_trunc.roots_ns")) +
				static_cast<double>(pre_stats.get_ns("precompute_no_trunc.fft_b_ns"))) /
			1e6);
	}

        #if 0
        static void bench_middle_product_basic_ifft_trunc_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);
		TimingStats pre_stats;
		const auto precomputed_b = precompute_middle_product_b(
			std::span<const G2Point>(B_orig), n, MiddleProductAlgorithm::BasicFFT, &pre_stats);
		std::vector<GT> out(n, GT::identity());

		double time_fft_g1 = 0;
		double time_pairing = 0;
		double time_ifft_gt = 0;
		double time_roots = 0;
		double time_extract = 0;
		double fft_a_total = 0;
		double fft_a_exec = 0;
		double fft_a_skip = 0;
		double ifft_total = 0;
		double ifft_exec = 0;
		double ifft_skip = 0;
		double online_total = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			TimingStats iter_stats;
			middle_product_basic_ifft_truncated(
				std::span<const G1Point>(A_orig),
				precomputed_b,
				std::span<GT>(out),
				&iter_stats);
			benchmark::ClobberMemory();

			time_fft_g1 += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.fft_a_ns"));
			time_pairing += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.pointwise_ns"));
			time_ifft_gt += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.ifft_ns"));
			time_roots += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.roots_ns"));
			time_roots += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.inv_roots_ns"));
			time_extract += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.extract_ns"));
			fft_a_total += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.fft_a.butterflies_total"));
			fft_a_exec += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.fft_a.butterflies_executed"));
			fft_a_skip += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.fft_a.butterflies_skipped"));
			ifft_total += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.ifft.butterflies_total"));
			ifft_exec += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.ifft.butterflies_executed"));
			ifft_skip += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.ifft.butterflies_skipped"));
			online_total += static_cast<double>(iter_stats.get_ns("basic_ifft_trunc.total_ns"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(precomputed_b.fft_size);
		state.counters["fft_g1_ms"] = benchmark::Counter(time_fft_g1 / 1e6 / iterations);
		state.counters["pairing_ms"] = benchmark::Counter(time_pairing / 1e6 / iterations);
		state.counters["ifft_gt_ms"] = benchmark::Counter(time_ifft_gt / 1e6 / iterations);
		state.counters["roots_ms"] = benchmark::Counter(time_roots / 1e6 / iterations);
		state.counters["extract_ms"] = benchmark::Counter(time_extract / 1e6 / iterations);
		state.counters["online_total_ms"] = benchmark::Counter(online_total / 1e6 / iterations);
		state.counters["online_compute_ms"] = benchmark::Counter((online_total - time_roots) / 1e6 / iterations);
		state.counters["fft_a_bfly_total"] = benchmark::Counter(fft_a_total / iterations);
		state.counters["fft_a_bfly_exec"] = benchmark::Counter(fft_a_exec / iterations);
		state.counters["fft_a_bfly_skip"] = benchmark::Counter(fft_a_skip / iterations);
		state.counters["fft_a_skip_ratio"] = benchmark::Counter(
			(fft_a_total > 0.0) ? (fft_a_skip / fft_a_total) : 0.0);
		state.counters["ifft_bfly_total"] = benchmark::Counter(ifft_total / iterations);
		state.counters["ifft_bfly_exec"] = benchmark::Counter(ifft_exec / iterations);
		state.counters["ifft_bfly_skip"] = benchmark::Counter(ifft_skip / iterations);
		state.counters["ifft_skip_ratio"] = benchmark::Counter(
			(ifft_total > 0.0) ? (ifft_skip / ifft_total) : 0.0);
		state.counters["offline_precompute_once_ms"] = benchmark::Counter(
			(static_cast<double>(pre_stats.get_ns("precompute.roots_ns")) +
				static_cast<double>(pre_stats.get_ns("precompute.basic_fft_b_ns"))) /
			1e6);
	}

	static void bench_middle_product_basic_pruned_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);
		TimingStats pre_stats;
		const auto precomputed_b = precompute_middle_product_b(
			std::span<const G2Point>(B_orig), n, MiddleProductAlgorithm::BasicFFT, &pre_stats);
		std::vector<GT> out(n, GT::identity());
		MiddleProductBasicPrunedWorkspace workspace;

		double time_fft_g1 = 0;
		double time_pairing = 0;
		double time_ifft_gt = 0;
		double time_roots = 0;
		double time_extract = 0;
		double fft_a_total = 0;
		double fft_a_exec = 0;
		double fft_a_skip = 0;
		double ifft_total = 0;
		double ifft_exec = 0;
		double ifft_skip = 0;
		double online_total = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			TimingStats iter_stats;
			middle_product_basic_pruned(
				std::span<const G1Point>(A_orig),
				precomputed_b,
				std::span<GT>(out),
				workspace,
				&iter_stats);
			benchmark::ClobberMemory();

			time_fft_g1 += static_cast<double>(iter_stats.get_ns("basic_pruned.fft_a_ns"));
			time_pairing += static_cast<double>(iter_stats.get_ns("basic_pruned.pointwise_ns"));
			time_ifft_gt += static_cast<double>(iter_stats.get_ns("basic_pruned.ifft_ns"));
			time_roots += static_cast<double>(iter_stats.get_ns("basic_pruned.roots_ns"));
			time_roots += static_cast<double>(iter_stats.get_ns("basic_pruned.inv_roots_ns"));
			time_extract += static_cast<double>(iter_stats.get_ns("basic_pruned.extract_ns"));
			fft_a_total += static_cast<double>(iter_stats.get_ns("basic_pruned.fft_a.butterflies_total"));
			fft_a_exec += static_cast<double>(iter_stats.get_ns("basic_pruned.fft_a.butterflies_executed"));
			fft_a_skip += static_cast<double>(iter_stats.get_ns("basic_pruned.fft_a.butterflies_skipped"));
			ifft_total += static_cast<double>(iter_stats.get_ns("basic_pruned.ifft.butterflies_total"));
			ifft_exec += static_cast<double>(iter_stats.get_ns("basic_pruned.ifft.butterflies_executed"));
			ifft_skip += static_cast<double>(iter_stats.get_ns("basic_pruned.ifft.butterflies_skipped"));
			online_total += static_cast<double>(iter_stats.get_ns("basic_pruned.total_ns"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(precomputed_b.fft_size);
		state.counters["fft_g1_ms"] = benchmark::Counter(time_fft_g1 / 1e6 / iterations);
		state.counters["pairing_ms"] = benchmark::Counter(time_pairing / 1e6 / iterations);
		state.counters["ifft_gt_ms"] = benchmark::Counter(time_ifft_gt / 1e6 / iterations);
		state.counters["roots_ms"] = benchmark::Counter(time_roots / 1e6 / iterations);
		state.counters["extract_ms"] = benchmark::Counter(time_extract / 1e6 / iterations);
		state.counters["online_total_ms"] = benchmark::Counter(online_total / 1e6 / iterations);
		state.counters["online_compute_ms"] = benchmark::Counter((online_total - time_roots) / 1e6 / iterations);
		state.counters["fft_a_bfly_total"] = benchmark::Counter(fft_a_total / iterations);
		state.counters["fft_a_bfly_exec"] = benchmark::Counter(fft_a_exec / iterations);
		state.counters["fft_a_bfly_skip"] = benchmark::Counter(fft_a_skip / iterations);
		state.counters["fft_a_skip_ratio"] = benchmark::Counter(
			(fft_a_total > 0.0) ? (fft_a_skip / fft_a_total) : 0.0);
		state.counters["ifft_bfly_total"] = benchmark::Counter(ifft_total / iterations);
		state.counters["ifft_bfly_exec"] = benchmark::Counter(ifft_exec / iterations);
		state.counters["ifft_bfly_skip"] = benchmark::Counter(ifft_skip / iterations);
		state.counters["ifft_skip_ratio"] = benchmark::Counter(
			(ifft_total > 0.0) ? (ifft_skip / ifft_total) : 0.0);
		state.counters["offline_precompute_once_ms"] = benchmark::Counter(
			(static_cast<double>(pre_stats.get_ns("precompute.roots_ns")) +
				static_cast<double>(pre_stats.get_ns("precompute.basic_fft_b_ns"))) /
			1e6);
	}

	static void bench_middle_product_specialized_trunc_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);
		TimingStats pre_stats;
		const auto precomputed_b = precompute_middle_product_b(
			std::span<const G2Point>(B_orig), n, MiddleProductAlgorithm::SpecializedFFT, &pre_stats);
		std::vector<GT> out(n, GT::identity());

		double time_fft_g1 = 0;
		double time_pairing = 0;
		double time_ifft_gt = 0;
		double time_derive = 0;
		double time_roots = 0;
		double time_extract = 0;
		double fft_a_total = 0;
		double fft_a_exec = 0;
		double fft_a_skip = 0;
		double ifft1_total = 0;
		double ifft1_exec = 0;
		double ifft1_skip = 0;
		double ifft2_total = 0;
		double ifft2_exec = 0;
		double ifft2_skip = 0;
		double online_total = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			TimingStats iter_stats;
			middle_product_specialized_truncated(
				std::span<const G1Point>(A_orig),
				precomputed_b,
				std::span<GT>(out),
				&iter_stats);
			benchmark::ClobberMemory();

			time_fft_g1 += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase1_fft_a_ns"));
			time_pairing += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase1_pointwise_ns"));
			time_pairing += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase2_pointwise_ns"));
			time_ifft_gt += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase1_ifft_ns"));
			time_ifft_gt += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase2_ifft_ns"));
			time_derive += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase2_derive_ntt_ns"));
			time_roots += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.roots_ns"));
			time_roots += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.inv_roots_ns"));
			time_extract += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase1_copy_ns"));
			time_extract += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase2_combine_ns"));
			fft_a_total += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase1_fft_a.butterflies_total"));
			fft_a_exec += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase1_fft_a.butterflies_executed"));
			fft_a_skip += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase1_fft_a.butterflies_skipped"));
			ifft1_total += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase1_ifft.butterflies_total"));
			ifft1_exec += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase1_ifft.butterflies_executed"));
			ifft1_skip += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase1_ifft.butterflies_skipped"));
			ifft2_total += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase2_ifft.butterflies_total"));
			ifft2_exec += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase2_ifft.butterflies_executed"));
			ifft2_skip += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.phase2_ifft.butterflies_skipped"));
			online_total += static_cast<double>(iter_stats.get_ns("middle_spec_trunc.total_ns"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(precomputed_b.fft_size);
		state.counters["fft_g1_ms"] = benchmark::Counter(time_fft_g1 / 1e6 / iterations);
		state.counters["pairing_ms"] = benchmark::Counter(time_pairing / 1e6 / iterations);
		state.counters["ifft_gt_ms"] = benchmark::Counter(time_ifft_gt / 1e6 / iterations);
		state.counters["derive_rev_ms"] = benchmark::Counter(time_derive / 1e6 / iterations);
		state.counters["roots_ms"] = benchmark::Counter(time_roots / 1e6 / iterations);
		state.counters["extract_ms"] = benchmark::Counter(time_extract / 1e6 / iterations);
		state.counters["online_total_ms"] = benchmark::Counter(online_total / 1e6 / iterations);
		state.counters["online_compute_ms"] = benchmark::Counter((online_total - time_roots) / 1e6 / iterations);
		state.counters["fft_a_bfly_total"] = benchmark::Counter(fft_a_total / iterations);
		state.counters["fft_a_bfly_exec"] = benchmark::Counter(fft_a_exec / iterations);
		state.counters["fft_a_bfly_skip"] = benchmark::Counter(fft_a_skip / iterations);
		state.counters["fft_a_skip_ratio"] = benchmark::Counter(
			(fft_a_total > 0.0) ? (fft_a_skip / fft_a_total) : 0.0);
		state.counters["ifft1_bfly_total"] = benchmark::Counter(ifft1_total / iterations);
		state.counters["ifft1_bfly_exec"] = benchmark::Counter(ifft1_exec / iterations);
		state.counters["ifft1_bfly_skip"] = benchmark::Counter(ifft1_skip / iterations);
		state.counters["ifft1_skip_ratio"] = benchmark::Counter(
			(ifft1_total > 0.0) ? (ifft1_skip / ifft1_total) : 0.0);
		state.counters["ifft2_bfly_total"] = benchmark::Counter(ifft2_total / iterations);
		state.counters["ifft2_bfly_exec"] = benchmark::Counter(ifft2_exec / iterations);
		state.counters["ifft2_bfly_skip"] = benchmark::Counter(ifft2_skip / iterations);
		state.counters["ifft2_skip_ratio"] = benchmark::Counter(
			(ifft2_total > 0.0) ? (ifft2_skip / ifft2_total) : 0.0);
		state.counters["offline_precompute_once_ms"] = benchmark::Counter(
			(static_cast<double>(pre_stats.get_ns("precompute.roots_ns")) +
				static_cast<double>(pre_stats.get_ns("precompute.spec_fft_b1_ns")) +
				static_cast<double>(pre_stats.get_ns("precompute.spec_fft_b0_ns"))) /
			1e6);
	}

        #endif

        static void bench_middle_product_one_ifft_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);
		TimingStats pre_stats;
		const auto precomputed_b = precompute_middle_product_b(
			std::span<const G2Point>(B_orig), n, &pre_stats);
		std::vector<GT> out(n, GT::identity());
		MiddleProductWorkspace workspace;

		double time_fft_g1 = 0;
                double time_pairing = 0;
                double time_ifft_gt = 0;
                double time_fill = 0;
                double time_roots = 0;
                double time_extract = 0;
                double online_total = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			TimingStats iter_stats;
				middle_product(
					std::span<const G1Point>(A_orig),
					precomputed_b,
					std::span<GT>(out),
					workspace,
					&iter_stats);
			benchmark::ClobberMemory();

				time_fft_g1 += static_cast<double>(iter_stats.get_ns("middle_product_pre.fft_a_ns"));
	                        time_pairing += static_cast<double>(iter_stats.get_ns("middle_product_pre.pointwise_ns"));
	                        time_ifft_gt += static_cast<double>(iter_stats.get_ns("middle_product_pre.ifft_ns"));
	                        time_fill += static_cast<double>(iter_stats.get_ns("middle_product_pre.fill_a_fft_ns"));
	                        time_fill += static_cast<double>(iter_stats.get_ns("middle_product_pre.fill_c_fft_ns"));
	                        time_extract += static_cast<double>(iter_stats.get_ns("middle_product_pre.extract_ns"));
				online_total += static_cast<double>(iter_stats.get_ns("middle_product_pre.total_ns"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(precomputed_b.fft_size);
                state.counters["fft_g1_ms"] = benchmark::Counter(time_fft_g1 / 1e6 / iterations);
                state.counters["pairing_ms"] = benchmark::Counter(time_pairing / 1e6 / iterations);
                state.counters["ifft_gt_ms"] = benchmark::Counter(time_ifft_gt / 1e6 / iterations);
                state.counters["fill_ms"] = benchmark::Counter(time_fill / 1e6 / iterations);
                state.counters["roots_ms"] = benchmark::Counter(time_roots / 1e6 / iterations);
                state.counters["extract_ms"] = benchmark::Counter(time_extract / 1e6 / iterations);
		state.counters["online_total_ms"] = benchmark::Counter(online_total / 1e6 / iterations);
		state.counters["online_compute_ms"] = benchmark::Counter((online_total - time_roots) / 1e6 / iterations);
			state.counters["offline_precompute_once_ms"] = benchmark::Counter(
				(static_cast<double>(pre_stats.get_ns("middle_product_precompute.roots_ns")) +
					static_cast<double>(pre_stats.get_ns("middle_product_precompute.fft_b0_ns")) +
					static_cast<double>(pre_stats.get_ns("middle_product_precompute.fft_b1_ns")) +
					static_cast<double>(pre_stats.get_ns("middle_product_precompute.fold_ns"))) /
				1e6);
		}

        static void bench_middle_product_scalar_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);
		TimingStats pre_stats;
		const auto precomputed_b = precompute_middle_product_b_scalar(
			std::span<const G2Point>(B_orig), n, &pre_stats);
		std::vector<GT> out(n, GT::identity());
		MiddleProductWorkspaceScalar workspace;

		double time_fft_g1 = 0;
                double time_pairing = 0;
                double time_ifft_gt = 0;
                double time_fill = 0;
                double time_roots = 0;
                double time_extract = 0;
                double online_total = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			TimingStats iter_stats;
				middle_product_scalar(
					std::span<const G1Point>(A_orig),
					precomputed_b,
					std::span<GT>(out),
					workspace,
					&iter_stats);
			benchmark::ClobberMemory();

				time_fft_g1 += static_cast<double>(iter_stats.get_ns("middle_product_pre.fft_a_ns"));
	                        time_pairing += static_cast<double>(iter_stats.get_ns("middle_product_pre.pointwise_ns"));
	                        time_ifft_gt += static_cast<double>(iter_stats.get_ns("middle_product_pre.ifft_ns"));
	                        time_fill += static_cast<double>(iter_stats.get_ns("middle_product_pre.fill_a_fft_ns"));
	                        time_fill += static_cast<double>(iter_stats.get_ns("middle_product_pre.fill_c_fft_ns"));
	                        time_extract += static_cast<double>(iter_stats.get_ns("middle_product_pre.extract_ns"));
				online_total += static_cast<double>(iter_stats.get_ns("middle_product_pre.total_ns"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(precomputed_b.fft_size);
                state.counters["fft_g1_ms"] = benchmark::Counter(time_fft_g1 / 1e6 / iterations);
                state.counters["pairing_ms"] = benchmark::Counter(time_pairing / 1e6 / iterations);
                state.counters["ifft_gt_ms"] = benchmark::Counter(time_ifft_gt / 1e6 / iterations);
                state.counters["fill_ms"] = benchmark::Counter(time_fill / 1e6 / iterations);
                state.counters["roots_ms"] = benchmark::Counter(time_roots / 1e6 / iterations);
                state.counters["extract_ms"] = benchmark::Counter(time_extract / 1e6 / iterations);
		state.counters["online_total_ms"] = benchmark::Counter(online_total / 1e6 / iterations);
		state.counters["online_compute_ms"] = benchmark::Counter((online_total - time_roots) / 1e6 / iterations);
			state.counters["offline_precompute_once_ms"] = benchmark::Counter(
				(static_cast<double>(pre_stats.get_ns("middle_product_precompute.roots_ns")) +
					static_cast<double>(pre_stats.get_ns("middle_product_precompute.fft_b0_ns")) +
					static_cast<double>(pre_stats.get_ns("middle_product_precompute.fft_b1_ns")) +
					static_cast<double>(pre_stats.get_ns("middle_product_precompute.fold_ns"))) /
				1e6);
		}

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
	static void bench_middle_product_avx(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);
		const auto precomputed_b = precompute_middle_product_b_avx(
			std::span<const G2Point>(B_orig), n);
		std::vector<GT> out(n, GT::identity());
		MiddleProductWorkspaceAvx workspace;

		for (auto _ : state) {
			benchmark::DoNotOptimize(out.data());
			middle_product_avx(
				std::span<const G1Point>(A_orig),
				precomputed_b,
				std::span<GT>(out),
				workspace);
			benchmark::ClobberMemory();
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(precomputed_b.fft_size);
	}

	static void bench_middle_product_avx_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);
		TimingStats pre_stats;
		const auto precomputed_b = precompute_middle_product_b_avx(
			std::span<const G2Point>(B_orig), n, &pre_stats);
		std::vector<GT> out(n, GT::identity());
		MiddleProductWorkspaceAvx workspace;

		double time_fft_g1 = 0;
                double time_pairing = 0;
                double time_ifft_gt = 0;
                double time_fill = 0;
                double time_roots = 0;
                double time_extract = 0;
                double online_total = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			TimingStats iter_stats;
				middle_product_avx(
					std::span<const G1Point>(A_orig),
					precomputed_b,
					std::span<GT>(out),
					workspace,
					&iter_stats);
			benchmark::ClobberMemory();

				time_fft_g1 += static_cast<double>(iter_stats.get_ns("middle_product_pre.fft_a_ns"));
	                        time_pairing += static_cast<double>(iter_stats.get_ns("middle_product_pre.pointwise_ns"));
	                        time_ifft_gt += static_cast<double>(iter_stats.get_ns("middle_product_pre.ifft_ns"));
	                        time_fill += static_cast<double>(iter_stats.get_ns("middle_product_pre.fill_a_fft_ns"));
	                        time_fill += static_cast<double>(iter_stats.get_ns("middle_product_pre.fill_c_fft_ns"));
	                        time_extract += static_cast<double>(iter_stats.get_ns("middle_product_pre.extract_ns"));
				online_total += static_cast<double>(iter_stats.get_ns("middle_product_pre.total_ns"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(precomputed_b.fft_size);
                state.counters["fft_g1_ms"] = benchmark::Counter(time_fft_g1 / 1e6 / iterations);
                state.counters["pairing_ms"] = benchmark::Counter(time_pairing / 1e6 / iterations);
                state.counters["ifft_gt_ms"] = benchmark::Counter(time_ifft_gt / 1e6 / iterations);
                state.counters["fill_ms"] = benchmark::Counter(time_fill / 1e6 / iterations);
                state.counters["roots_ms"] = benchmark::Counter(time_roots / 1e6 / iterations);
                state.counters["extract_ms"] = benchmark::Counter(time_extract / 1e6 / iterations);
		state.counters["online_total_ms"] = benchmark::Counter(online_total / 1e6 / iterations);
		state.counters["online_compute_ms"] = benchmark::Counter((online_total - time_roots) / 1e6 / iterations);
			state.counters["offline_precompute_once_ms"] = benchmark::Counter(
				(static_cast<double>(pre_stats.get_ns("middle_product_precompute.roots_ns")) +
					static_cast<double>(pre_stats.get_ns("middle_product_precompute.fft_b0_ns")) +
					static_cast<double>(pre_stats.get_ns("middle_product_precompute.fft_b1_ns")) +
					static_cast<double>(pre_stats.get_ns("middle_product_precompute.fold_ns"))) /
				1e6);
		}
#endif

        // Benchmark middle_product_fft: computes middle n coefficients of (size n) * (size 2n) convolution
#if 0
        static void bench_middle_product_fft(benchmark::State& state)
        {
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<G1Point> A_orig = make_g1_data(n, /*seed*/ 0xDEADBEEFULL);
		const std::vector<G2Point> B_orig = make_g2_data(2 * n, /*seed*/ 0xCAFEBABEULL);
		const auto precomputed_b = precompute_middle_product_b(
			std::span<const G2Point>(B_orig), n, MiddleProductAlgorithm::SpecializedFFT);
		std::vector<GT> out(n, GT::identity());

		for (auto _ : state) {
			benchmark::DoNotOptimize(out.data());
			middle_product(
				std::span<const G1Point>(A_orig),
				precomputed_b,
				std::span<GT>(out));

			benchmark::ClobberMemory();
		}

		state.counters["n"] = static_cast<double>(n);
                state.counters["m"] = static_cast<double>(precomputed_b.fft_size);
                state.counters["ops/sec"] = benchmark::Counter(1, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
        }
#endif

	// Benchmark middle_product_fft_truncated: uses truncated inverse FFT
#if 0
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
#endif

	// Standalone G1 forward FFT in middle-product shape:
	// input has n non-zero coefficients embedded in m.
	static void bench_g1_fft_full_middle_shape_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
			const size_t m = middle_product_basic_fft_size(n);
		const size_t xn = fft_admissible_size(n, m);
		const std::vector<Fr> roots = make_roots(m);
		const std::vector<G1Point> a_orig = make_g1_data(n, /*seed*/ 0xA11CEULL);
		std::vector<G1Point> a(m, G1Point::infinity());

		double fft_ns = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			for (size_t i = 0; i < n; ++i) {
				a[i] = a_orig[i];
			}
			for (size_t i = n; i < m; ++i) {
				a[i] = G1Point::infinity();
			}

			benchmark::DoNotOptimize(a.data());
			TimingStats stats;
			{
				ScopedTimer timer(&stats, "g1_fft_full.total_ns");
				fft<G1Point, Fr, G1CoeffCtx>(std::span<G1Point>(a), std::span<const Fr>(roots), false);
			}
			benchmark::ClobberMemory();

			fft_ns += static_cast<double>(stats.get_ns("g1_fft_full.total_ns"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(m);
		state.counters["xn"] = static_cast<double>(xn);
		state.counters["fft_ms"] = benchmark::Counter(fft_ns / 1e6 / iterations);
	}

	// Standalone truncated G1 forward FFT in middle-product shape:
	// only the active prefix xn contributes work.
	static void bench_g1_fft_trunc_middle_shape_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
			const size_t m = middle_product_basic_fft_size(n);
		const size_t xn = fft_admissible_size(n, m);
		const std::vector<Fr> roots = make_roots(m);
		const std::vector<G1Point> a_orig = make_g1_data(n, /*seed*/ 0xA11CEULL);
		const FftTruncatedForwardPlan plan = make_fft_truncated_forward_plan(m, xn);
		std::vector<G1Point> a(m, G1Point::infinity());

		double fft_ns = 0;
		double bfly_total = 0;
		double bfly_exec = 0;
		double bfly_skip = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			for (size_t i = 0; i < n; ++i) {
				a[i] = a_orig[i];
			}
			for (size_t i = n; i < m; ++i) {
				a[i] = G1Point::infinity();
			}

			benchmark::DoNotOptimize(a.data());
			TimingStats stats;
			fft_forward_truncated_prefix<G1Point, Fr, G1CoeffCtx>(
				std::span<G1Point>(a),
				std::span<const Fr>(roots),
				xn,
				plan,
				&stats,
				"g1_fft_trunc");
			benchmark::ClobberMemory();

			fft_ns += static_cast<double>(stats.get_ns("g1_fft_trunc.total_ns"));
			bfly_total += static_cast<double>(stats.get_ns("g1_fft_trunc.butterflies_total"));
			bfly_exec += static_cast<double>(stats.get_ns("g1_fft_trunc.butterflies_executed"));
			bfly_skip += static_cast<double>(stats.get_ns("g1_fft_trunc.butterflies_skipped"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(m);
		state.counters["xn"] = static_cast<double>(xn);
		state.counters["fft_ms"] = benchmark::Counter(fft_ns / 1e6 / iterations);
		state.counters["bfly_total"] = benchmark::Counter(bfly_total / iterations);
		state.counters["bfly_exec"] = benchmark::Counter(bfly_exec / iterations);
		state.counters["bfly_skip"] = benchmark::Counter(bfly_skip / iterations);
		state.counters["bfly_skip_ratio"] = benchmark::Counter(
			(bfly_total > 0.0) ? (bfly_skip / bfly_total) : 0.0);
	}

	// Standalone GT inverse FFT in one-iFFT shape:
	// compare full iFFT(m) vs truncated-range iFFT for the used window [n, 2n).
	static void bench_gt_ifft_full_middle_shape_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
			const size_t m = middle_product_fft_size(n);
		const std::vector<Fr> roots = make_roots(m);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);
		const std::vector<GT> time_domain = make_gt_data(m);
		std::vector<GT> freq_domain = time_domain;
		fft<GT, Fr, GTCoeffCtx>(freq_domain, roots, /*inverse*/ false);
		std::vector<GT> data = freq_domain;

		double ifft_ns = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			data = freq_domain;
			benchmark::DoNotOptimize(data.data());
			TimingStats stats;
			{
				ScopedTimer timer(&stats, "gt_ifft_full.total_ns");
				fft<GT, Fr, GTCoeffCtx>(std::span<GT>(data), std::span<const Fr>(inv_roots), true);
			}
			benchmark::ClobberMemory();

			ifft_ns += static_cast<double>(stats.get_ns("gt_ifft_full.total_ns"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(m);
		state.counters["ifft_ms"] = benchmark::Counter(ifft_ns / 1e6 / iterations);
	}

	static void bench_gt_ifft_trunc_middle_shape_detailed(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
			const size_t m = middle_product_fft_size(n);
		const std::vector<Fr> roots = make_roots(m);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);
		const std::vector<GT> time_domain = make_gt_data(m);
		std::vector<GT> freq_domain = time_domain;
		fft<GT, Fr, GTCoeffCtx>(freq_domain, roots, /*inverse*/ false);
		std::vector<GT> data = freq_domain;

		double ifft_ns = 0;
		double bfly_total = 0;
		double bfly_exec = 0;
		double bfly_skip = 0;
		int64_t iterations = 0;

		for (auto _ : state) {
			data = freq_domain;
			benchmark::DoNotOptimize(data.data());
			TimingStats stats;
			{
				ScopedTimer timer(&stats, "gt_ifft_trunc.total_ns");
				fft_inverse_truncated_range<GT, Fr, GTCoeffCtx>(
					std::span<GT>(data),
					std::span<const Fr>(inv_roots),
					n,
					n,
					&stats,
					"gt_ifft_trunc");
			}
			benchmark::ClobberMemory();

			ifft_ns += static_cast<double>(stats.get_ns("gt_ifft_trunc.total_ns"));
			bfly_total += static_cast<double>(stats.get_ns("gt_ifft_trunc.butterflies_total"));
			bfly_exec += static_cast<double>(stats.get_ns("gt_ifft_trunc.butterflies_executed"));
			bfly_skip += static_cast<double>(stats.get_ns("gt_ifft_trunc.butterflies_skipped"));
			iterations++;
		}

		state.counters["n"] = static_cast<double>(n);
		state.counters["m"] = static_cast<double>(m);
		state.counters["ifft_ms"] = benchmark::Counter(ifft_ns / 1e6 / iterations);
		state.counters["bfly_total"] = benchmark::Counter(bfly_total / iterations);
		state.counters["bfly_exec"] = benchmark::Counter(bfly_exec / iterations);
		state.counters["bfly_skip"] = benchmark::Counter(bfly_skip / iterations);
		state.counters["bfly_skip_ratio"] = benchmark::Counter(
			(bfly_total > 0.0) ? (bfly_skip / bfly_total) : 0.0);
	}
}


// Detailed benchmark breaking down middle_product_fft into phases
#if 0
static void bench_middle_product_fft_detailed(benchmark::State& state)
{
	const size_t n = static_cast<size_t>(state.range(0));
	const std::vector<G1Point> A_orig = make_g1_data(n, 0xDEADBEEFULL);
	const std::vector<G2Point> B_orig = make_g2_data(2 * n, 0xCAFEBABEULL);
	TimingStats pre_stats;
	const auto precomputed_b = precompute_middle_product_b(
		std::span<const G2Point>(B_orig), n, MiddleProductAlgorithm::SpecializedFFT, &pre_stats);
	std::vector<GT> out(n, GT::identity());

	// Accumulators for phase timings (in nanoseconds)
	double time_fft_g1 = 0;
	double time_fft_g2 = 0; // online path uses precomputed B
	double time_pairing = 0;
	double time_ifft_gt = 0;
	double time_derive_rev = 0;
	double time_combine = 0;
	double time_roots = 0;
	int64_t iterations = 0;

	for (auto _ : state) {
		TimingStats iter_stats;
		middle_product(
			std::span<const G1Point>(A_orig),
			precomputed_b,
			std::span<GT>(out),
			&iter_stats);
		benchmark::ClobberMemory();

		time_fft_g1 += static_cast<double>(iter_stats.get_ns("middle.phase1_fft_a_ns"));
		time_fft_g2 += static_cast<double>(iter_stats.get_ns("middle.phase1_fft_b1_ns"));
		time_fft_g2 += static_cast<double>(iter_stats.get_ns("middle.phase2_fft_b0_ns"));
		time_pairing += static_cast<double>(iter_stats.get_ns("middle.phase1_pointwise_ns"));
		time_pairing += static_cast<double>(iter_stats.get_ns("middle.phase2_pointwise_ns"));
		time_ifft_gt += static_cast<double>(iter_stats.get_ns("middle.phase1_ifft_ns"));
		time_ifft_gt += static_cast<double>(iter_stats.get_ns("middle.phase2_ifft_ns"));
		time_derive_rev += static_cast<double>(iter_stats.get_ns("middle.phase2_derive_ntt_ns"));
		time_combine += static_cast<double>(iter_stats.get_ns("middle.phase2_combine_ns"));
		time_roots += static_cast<double>(iter_stats.get_ns("middle_pre.roots_ns"));
		time_roots += static_cast<double>(iter_stats.get_ns("middle_pre.inv_roots_ns"));
		iterations++;
	}

	state.counters["n"] = static_cast<double>(n);
	state.counters["m"] = static_cast<double>(precomputed_b.fft_size);
	state.counters["fft_g1_ms"] = benchmark::Counter(time_fft_g1 / 1e6 / iterations);
	state.counters["fft_g2_ms"] = benchmark::Counter(time_fft_g2 / 1e6 / iterations);
	state.counters["pairing_ms"] = benchmark::Counter(time_pairing / 1e6 / iterations);
	state.counters["ifft_gt_ms"] = benchmark::Counter(time_ifft_gt / 1e6 / iterations);
	state.counters["derive_rev_ms"] = benchmark::Counter(time_derive_rev / 1e6 / iterations);
	state.counters["combine_ms"] = benchmark::Counter(time_combine / 1e6 / iterations);
	state.counters["roots_ms"] = benchmark::Counter(time_roots / 1e6 / iterations);
        state.counters["offline_precompute_ms"] = benchmark::Counter(
                (static_cast<double>(pre_stats.get_ns("precompute.roots_ns")) +
                        static_cast<double>(pre_stats.get_ns("precompute.spec_fft_b1_ns")) +
                        static_cast<double>(pre_stats.get_ns("precompute.spec_fft_b0_ns"))) /
                1e6);
}
#endif

// Keep sizes as powers of two (roots exist up to 2^32 for Fr).
#define FFT_SIZES \
->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)->Arg(8192)->Arg(16384)

// Smaller sizes for expensive operations (middle product involves pairings)
#define MIDDLE_PROD_SIZES \
->Arg(2)->Arg(4)->Arg(8)->Arg(12)->Arg(16)->Arg(24)->Arg(32)->Arg(48)->Arg(64) \
->DenseRange(96, 512, 32)

// Medium sweep for detailed reports (non power-of-two n included).
#define MIDDLE_PROD_SWEEP_SIZES \
->Arg(2)->Arg(4)->Arg(8)->Arg(12)->Arg(16)->Arg(24)->Arg(32)->Arg(48)->Arg(64) \
->DenseRange(96, 512, 32)

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

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
BENCHMARK(bench_gt_ifft_avx52_precomputed_cpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
#endif


// ============================================================================
// Truncated and other FFT benchmarks
// ============================================================================

#if 0
BENCHMARK(bench_gt_ifft_truncated)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
#endif

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

#if 0
BENCHMARK(bench_middle_product_fft_truncated)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
#endif

BENCHMARK(bench_middle_product_basic_detailed)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_basic_no_trunc_detailed)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_scalar_detailed)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
BENCHMARK(bench_middle_product_avx_detailed)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_middle_product_avx)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
#endif

BENCHMARK(bench_middle_product_one_ifft_detailed)
MIDDLE_PROD_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_g1_fft_full_middle_shape_detailed)
MIDDLE_PROD_SWEEP_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_g1_fft_trunc_middle_shape_detailed)
MIDDLE_PROD_SWEEP_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_gt_ifft_full_middle_shape_detailed)
MIDDLE_PROD_SWEEP_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_gt_ifft_trunc_middle_shape_detailed)
MIDDLE_PROD_SWEEP_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
//BENCHMARK_MAIN();
