#include <benchmark/benchmark.h>

#include <vector>
#include <random>

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
}

// Keep sizes as powers of two (roots exist up to 2^32 for Fr).
#define FFT_SIZES \
->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096)->Arg(8192)->Arg(16384)

BENCHMARK(bench_g1_fft_cpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

#if defined(BLST_CUDA)
BENCHMARK(bench_g1_fft_gpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
#endif

BENCHMARK(bench_gt_ifft_cpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();

#if defined(BLST_CUDA)
BENCHMARK(bench_gt_ifft_gpu)
FFT_SIZES
->Unit(benchmark::kMillisecond)
->UseRealTime();
#endif

//BENCHMARK_MAIN();//