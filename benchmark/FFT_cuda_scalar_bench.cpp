#include <benchmark/benchmark.h>

#include <chrono>
#include <random>
#include <vector>

#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/primitveRoot.h"
#include "blst/fields/Fr.h"

#if defined(BLST_CUDA)
#include "blst/fft/ctFft_cuda.h"
#endif

using namespace blst;

namespace
{
	static std::vector<Fr> make_roots(size_t n)
	{
		Fr prim = primRootOfUnity<Fr>(n);
		std::vector<Fr> roots(n);
		Fr cur = Fr::one();
		for (size_t i = 0; i < n; ++i) {
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

	static std::vector<G1Point> make_g1_data(size_t n, uint64_t seed)
	{
		std::mt19937_64 rng(seed);
		std::uniform_int_distribution<uint64_t> dist(1, (uint64_t)-1);

		std::vector<G1Point> data(n);
		for (size_t i = 0; i < n; ++i) {
			Fr s;
			Fr::from_uint64(s, dist(rng));
			Fr s_non_mont;
			Fr::from_mont(s_non_mont, s);
			G1Point::mult(data[i], G1Point::generator(), reinterpret_cast<const byte*>(&s_non_mont), 256);
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

	static std::vector<GT> make_gt_data(size_t n)
	{
		GT base = pair_gen();
		std::vector<GT> data(n);
		data[0] = GT::identity();

		GT acc = base;
		for (size_t i = 1; i < n; ++i) {
			data[i] = acc;
			GT::add(acc, acc, base);
		}
		return data;
	}

	static void bench_g1_fft_cpu_scalar(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<Fr> roots = make_roots(n);
		const std::vector<G1Point> original = make_g1_data(n, 0xC0FFEEULL);
		std::vector<G1Point> data = original;

		for (auto _ : state) {
			data = original;
			benchmark::DoNotOptimize(data.data());
			fft<G1Point, Fr, G1CoeffCtx>(data, roots, false);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(n);
		state.counters["FFTs/sec"] = benchmark::Counter(
			1,
			benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	static void bench_gt_ifft_cpu_scalar(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<Fr> roots = make_roots(n);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);
		const std::vector<GT> time_domain = make_gt_data(n);
		std::vector<GT> freq_domain = time_domain;
		fft<GT, Fr, GTCoeffCtx>(freq_domain, roots, false);

		std::vector<GT> data = freq_domain;
		for (auto _ : state) {
			data = freq_domain;
			benchmark::DoNotOptimize(data.data());
			fft<GT, Fr, GTCoeffCtx>(data, inv_roots, true);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(n);
		state.counters["iFFTs/sec"] = benchmark::Counter(
			1,
			benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

#if defined(BLST_CUDA)
	static void bench_g1_fft_gpu_convenience(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<Fr> roots = make_roots(n);
		const std::vector<G1Point> original = make_g1_data(n, 0xC0FFEEULL);
		std::vector<G1Point> data = original;
		std::vector<G1Point> warmup = original;

		fft_cuda(warmup, roots, false);

		for (auto _ : state) {
			data = original;
			benchmark::DoNotOptimize(data.data());
			fft_cuda(data, roots, false);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(n);
		state.counters["FFTs/sec"] = benchmark::Counter(
			1,
			benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	static void bench_g1_fft_gpu_prepared(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<Fr> roots = make_roots(n);
		const std::vector<G1Point> original = make_g1_data(n, 0xC0FFEEULL);
		std::vector<G1Point> data = original;

		const auto prepare_begin = std::chrono::steady_clock::now();
		CudaFftPreparedRoots<G1Point> prepared_roots(roots);
		CudaFftWorkspace<G1Point> workspace;
		workspace.ensure(n);
		const auto prepare_end = std::chrono::steady_clock::now();
		const double prepared_once_ms =
			std::chrono::duration<double, std::milli>(prepare_end - prepare_begin).count();
		std::vector<G1Point> warmup = original;

		fft_cuda(warmup, prepared_roots, workspace, false);

		for (auto _ : state) {
			data = original;
			benchmark::DoNotOptimize(data.data());
			fft_cuda(data, prepared_roots, workspace, false);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(n);
		state.counters["prepared_once_ms"] = benchmark::Counter(prepared_once_ms);
		state.counters["FFTs/sec"] = benchmark::Counter(
			1,
			benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	static void bench_gt_ifft_gpu_convenience(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<Fr> roots = make_roots(n);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);
		const std::vector<GT> time_domain = make_gt_data(n);
		std::vector<GT> freq_domain = time_domain;
		fft_cuda(freq_domain, roots, false);
		std::vector<GT> warmup = freq_domain;

		fft_cuda(warmup, inv_roots, true);

		std::vector<GT> data = freq_domain;
		for (auto _ : state) {
			data = freq_domain;
			benchmark::DoNotOptimize(data.data());
			fft_cuda(data, inv_roots, true);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(n);
		state.counters["iFFTs/sec"] = benchmark::Counter(
			1,
			benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}

	static void bench_gt_ifft_gpu_prepared(benchmark::State& state)
	{
		const size_t n = static_cast<size_t>(state.range(0));
		const std::vector<Fr> roots = make_roots(n);
		const std::vector<Fr> inv_roots = make_inv_roots(roots);

		const auto prepare_begin = std::chrono::steady_clock::now();
		CudaFftPreparedRoots<GT> prepared_roots(roots);
		CudaFftPreparedRoots<GT> prepared_inv_roots(inv_roots);
		CudaFftWorkspace<GT> workspace;
		workspace.ensure(n);
		const auto prepare_end = std::chrono::steady_clock::now();
		const double prepared_once_ms =
			std::chrono::duration<double, std::milli>(prepare_end - prepare_begin).count();

		const std::vector<GT> time_domain = make_gt_data(n);
		std::vector<GT> freq_domain = time_domain;
		fft_cuda(freq_domain, prepared_roots, workspace, false);
		std::vector<GT> warmup = freq_domain;

		fft_cuda(warmup, prepared_inv_roots, workspace, true);

		std::vector<GT> data = freq_domain;
		for (auto _ : state) {
			data = freq_domain;
			benchmark::DoNotOptimize(data.data());
			fft_cuda(data, prepared_inv_roots, workspace, true);
			benchmark::ClobberMemory();
		}

		state.counters["N"] = static_cast<double>(n);
		state.counters["prepared_once_ms"] = benchmark::Counter(prepared_once_ms);
		state.counters["iFFTs/sec"] = benchmark::Counter(
			1,
			benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
	}
#endif
}

#define FFT_COMPARE_SIZES \
->Arg(1024)->Arg(2048)->Arg(4096)

BENCHMARK(bench_g1_fft_cpu_scalar)
FFT_COMPARE_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_gt_ifft_cpu_scalar)
FFT_COMPARE_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

#if defined(BLST_CUDA)
BENCHMARK(bench_g1_fft_gpu_convenience)
FFT_COMPARE_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_g1_fft_gpu_prepared)
FFT_COMPARE_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_gt_ifft_gpu_convenience)
FFT_COMPARE_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

BENCHMARK(bench_gt_ifft_gpu_prepared)
FFT_COMPARE_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
#endif
