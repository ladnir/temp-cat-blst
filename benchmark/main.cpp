
#include <benchmark/benchmark.h>
#include <cstring>
#include <iostream>

#ifdef BLST_CUDA
#include "blst/fields/Fp12.h"
//#include "blst/backend/cuda/mul_mont_cuda.h"
#include <thrust/device_vector.h>
#include "blst/curve/point.h"
#include <span>
#include "parallel_prefix_bench.h"
#include <chrono>
#include "sqr_impl.h"
#include "blst/fields/FpTile.h"
//BENCHMARK_MAIN();
//void FpMain();
//void Fp12Main();
using namespace blst;



BLST_GLOBAL void main_fp12_kernel
(
	Fp12* d_data,
	u64   n
)
{
	u64 idx;
	//u64 stride;

	idx = blockIdx.x * blockDim.x + threadIdx.x;

	//stride = gridDim.x * blockDim.x;
	//for (u64 i = idx; i < n; i += stride)
	//{
	if (idx < n) {
		Fp12 temp;
		cyclotomic_sqr(temp, d_data[idx], NoopIdx{});
		//Fp12::cyclotomic_sqr(temp, d_data[idx]);
		d_data[idx] = temp;
	}
}

BLST_GLOBAL void main_fp12_kernel(
	SoAFp12 d_data,
	u64   n)
{
	u64 idx;
	//u64 stride;

	idx = blockIdx.x * blockDim.x + threadIdx.x;

	//stride = gridDim.x * blockDim.x;
	//for (u64 i = idx; i < n; i += stride)
	//{
	if (idx < n) {
		//Fp12 temp;
		cyclotomic_sqr(d_data, d_data, BasicIdx{ idx });
		//Fp12::cyclotomic_sqr(temp, d_data[idx]);
		//d_data[idx] = temp;
	}
}

__launch_bounds__(64, 8)
BLST_GLOBAL void main_fp12Tile_kernel(
	bench::Fp12Tile* d_data,
	u64   n)
{
	u64 idx;
	//u64 stride;

	idx = blockIdx.x * blockDim.x + threadIdx.x;

	//stride = gridDim.x * blockDim.x;
	//for (u64 i = idx; i < n; i += stride)
	//{
	if (idx < n) {

#ifdef CHECK_RES
		Fp12 in;
		for (u64 i = 0; i < 2; ++i)
			for (u64 j = 0; j < 3; ++j)
				for (u64 k = 0; k < 2; ++k)
					for (u64 l = 0; l < 6; ++l)
					{
						auto& ss = d_data[idx / 32].data[i][j].data[k][l / 2];
						in[i][j][k][l] = l & 1 ? ss.y : ss.x;
					}
		
		Fp12::cyclotomic_sqr(in, in);
#endif

		cyclotomic_sqr(d_data[idx / 32], d_data[idx / 32], NoopIdx{});

#ifdef CHECK_RES
		Fp12 out;
		for (u64 i = 0; i < 2; ++i)
		{
			for (u64 j = 0; j < 3; ++j)
				for (u64 k = 0; k < 2; ++k)
					for (u64 l = 0; l < 6; ++l)
					{
						auto& ss = d_data[idx / 32].data[i][j].data[k][l / 2];
						out[i][j][k][l] = l & 1 ? ss.y : ss.x;
					}
		}

		if (Fp12::equal(in, out) == false)
		{
			printf("failed");
		}
#endif

		//Fp12::cyclotomic_sqr(temp, d_data[idx]);
		//d_data[idx] = temp;
	}
}


__launch_bounds__(128, 4)
BLST_GLOBAL void blst_fp12Tile_kernel(
	blst::Fp12Tile* d_data,
	u64   n)
{
	u64 idx;
	//u64 stride;

	idx = blockIdx.x * blockDim.x + threadIdx.x;

	//stride = gridDim.x * blockDim.x;
	//for (u64 i = idx; i < n; i += stride)
	//{
	if (idx < n) {
		//blst::Fp12Tile::cyclotomic_sqr(d_data[idx / 32], d_data[idx / 32]);
		//Fp12 temp;
		//cyclotomic_sqr(d_data[idx / 32], d_data[idx / 32], NoopIdx{});
		
		//auto tile = d_data[idx / 32].load();
		//cyclotomic_sqr(tile, tile, NoopIdx{});
		//d_data[idx / 32].store(tile);

		cyclotomic_sqr(d_data[idx / 32], d_data[idx / 32], NoopIdx{});

		//Fp12::cyclotomic_sqr(temp, d_data[idx]);
		//d_data[idx] = temp;
	}
}

void Fp12Main()
{
	u64 n = 1 << 18;
	bool AoS = false;
	if (AoS)
	{



		//std::vector<Fp12> h_data(n);
		thrust::device_vector<Fp12> d_data(n);

		// DIM
		dim3 block(64);
		dim3 grid((n + block.x - 1) / block.x);
		//cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, 0);
		//grid = dim3((unsigned)std::max<u64>(grid.x, (u64)sms * 8));

		auto begin = std::chrono::high_resolution_clock::now();

		for (u64 i = 0; i < 10; ++i)
		{

			main_fp12_kernel << <grid, block >> > (thrust::raw_pointer_cast(d_data.data()), n);
			// wait for GPU to finish
			cudaDeviceSynchronize();
		}

		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end - begin;
		printf("Time taken (AoS): %f s\n", diff.count());
		printf("Throughput: %f ops/sec\n", n / (diff.count() / 1000.0));

	}
	else
	{
		bool cuda = true;
		SoAFp12 d_data;
		for (auto i = 0; i < 2; ++i)
		{
			for (auto j = 0; j < 3; ++j)
			{
				for (auto k = 0; k < 2; ++k)
				{
#ifdef USE_LONG2
					for (auto l = 0; l < 3; ++l)
					{
						ulong2* ptr;
						if (cuda)
						{
							cudaMalloc(&ptr, n * sizeof(ulong2));
						}
						else
							ptr = (ulong2*)malloc(n * sizeof(ulong2));
						d_data[i][j][k][l] = ptr;
					}

#else
					for (auto l = 0; l < 6; ++l)
					{
						//d_data.data[i][j][k][l];
						limb_t* ptr;
						if (cuda)
						{
							cudaMalloc(&ptr, n * sizeof(limb_t));
						}
						else
							ptr = (limb_t*)malloc(n * sizeof(limb_t));
						//d_data[i][j][k][l] = std::span<limb_t>(ptr, n);
						d_data[i][j][k][l] = ptr;
					}
#endif
				}
				//d_data.data[i][j].resize(n);
			}
		}
		// DIM
		dim3 block(64);
		dim3 grid((n + block.x - 1) / block.x);

		for (u64 i = 0; i < 10; ++i)
		{
			if (cuda)
			{
				main_fp12_kernel << <grid, block >> > (d_data, n);
				auto error = cudaGetLastError();
				if (error != cudaSuccess)
				{
					printf("CUDA error: %s\n", cudaGetErrorString(error));
				}
				error = cudaDeviceSynchronize();
				if (error != cudaSuccess)
				{
					printf("CUDA sync error: %s\n", cudaGetErrorString(error));
				}
			}
			else
			{
				for (u64 j = 0; j < n; ++j)
				{
					cyclotomic_sqr(d_data, d_data, BasicIdx{ j });
				}
			}
		}

		for (auto i = 0; i < 2; ++i)
		{
			for (auto j = 0; j < 3; ++j)
			{
				for (auto k = 0; k < 2; ++k)
				{
					for (auto l = 0; l < 6; ++l)
					{
						//auto ptr = d_data[i][j][k][l].data();
						auto ptr = d_data[i][j][k][l];
						//d_data.data[i][j][k][l];
						if (cuda)
							cudaFree(ptr);
						else
							free(ptr);

					}
				}
				//d_data.data[i][j].resize(n);
			}
		}
	}
	printf("success\n");
}


void benchFp12TileMain()
{
	u64 n = 1 << 18;
	bool cuda = true;
	using Tile = bench::Fp12Tile;
	Tile* d_data;
	auto num_tiles = (n + 31) / 32;
	if (cuda)
		cudaMalloc(&d_data, num_tiles * sizeof(Tile));
	else
		d_data = (Tile*)malloc(num_tiles * sizeof(Tile));

	// DIM
	dim3 block(64);
	dim3 grid((n + block.x - 1) / block.x);

	auto begin = std::chrono::high_resolution_clock::now();

	for (u64 i = 0; i < 10; ++i)
	{
		if (cuda)
		{
			main_fp12Tile_kernel << <grid, block >> > (d_data, n);
			auto error = cudaGetLastError();
			if (error != cudaSuccess)
			{
				printf("CUDA error: %s\n", cudaGetErrorString(error));
			}
			error = cudaDeviceSynchronize();
			if (error != cudaSuccess)
			{
				printf("CUDA sync error: %s\n", cudaGetErrorString(error));
			}
		}
		else
		{
			//for (u64 j = 0; j < n; ++j)
			//{
			//	cyclotomic_sqr(d_data, d_data, BasicIdx{ j });
			//}
		}
	}
	auto end = std::chrono::high_resolution_clock::now();

	if (cuda)
		cudaFree(d_data);
	else
		free(d_data);

	printf("success\n");
	auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	printf("time %d ms", (int)dur);
}


void blstFp12TileMain()
{
	u64 n = 1 << 18;
	bool cuda = true;
	blst::Fp12Tile* d_data;
	auto num_tiles = (n + 31) / 32;
	if (cuda)
		cudaMalloc(&d_data, num_tiles * sizeof(blst::Fp12Tile));
	else
		d_data = (blst::Fp12Tile*)malloc(num_tiles * sizeof(blst::Fp12Tile));

	// DIM
	dim3 block(128);
	dim3 grid((n + block.x - 1) / block.x);
	auto begin = std::chrono::high_resolution_clock::now();

	for (u64 i = 0; i < 10; ++i)
	{
		if (cuda)
		{
			blst_fp12Tile_kernel << <grid, block >> > (d_data, n);
			auto error = cudaGetLastError();
			if (error != cudaSuccess)
			{
				printf("CUDA error: %s\n", cudaGetErrorString(error));
			}
			error = cudaDeviceSynchronize();
			if (error != cudaSuccess)
			{
				printf("CUDA sync error: %s\n", cudaGetErrorString(error));
			}
		}
		else
		{
			//for (u64 j = 0; j < n; ++j)
			//{
			//	cyclotomic_sqr(d_data, d_data, BasicIdx{ j });
			//}
		}
	}
	auto end = std::chrono::high_resolution_clock::now();

	if (cuda)
		cudaFree(d_data);
	else
		free(d_data);

	printf("success\n");
	auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	printf("time %d ms", (int)dur);
}

BLST_GLOBAL void g1_add_kernel(G1Point* ret, blst::G1Point* a, blst::G1Point* b, size_t n)
{
	u64 idx;
	//u64 stride;

	idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n) {
		G1Point::add(ret[idx], a[idx], b[idx]);
	}
}


void G1Main()
{
	u64 n = 1 << 18;


	//std::vector<Fp12> h_data(n);
	thrust::device_vector<G1Point> a(n);
	thrust::device_vector<G1Point> b(n);
	thrust::device_vector<G1Point> c(n);

	// DIM
	dim3 block(32);
	dim3 grid((n + block.x - 1) / block.x);

	for (u64 i = 0; i < 10; ++i)
	{

		g1_add_kernel << <grid, block >> > (
			thrust::raw_pointer_cast(c.data()),
			thrust::raw_pointer_cast(a.data()),
			thrust::raw_pointer_cast(b.data()),
			n);
		// wait for GPU to finish
		cudaDeviceSynchronize();
	}
}



// Optimization: Single Base Pointer + Launch Bounds
// Forces compiler to use fewer registers, increasing occupancy.
// 256 threads, min 4 blocks => max 64 registers/thread
template<typename Field>
__global__ void bench_fp_kernel(
	SoAFp a_ptr,
	SoAFp b_ptr,
	uint64_t n_elems
) {
	constexpr int NLIMBS = Field::num_limb();
	auto elem_idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (elem_idx > n_elems)
		return;
#ifdef USE_LONG2
	Fp aa;
	for (u64 i = 0; i < NLIMBS / 2; ++i)
	{
		ulong2 val = a_ptr[i][elem_idx];
		aa[2 * i + 0] = val.x;
		aa[2 * i + 1] = val.y;
	}
	Fp bb;
	for (u64 i = 0; i < NLIMBS / 2; ++i)
	{
		ulong2 val = b_ptr[i][elem_idx];
		bb[2 * i + 0] = val.x;
		bb[2 * i + 1] = val.y;
	}
#else
	Fp aa;
	for (u64 i = 0; i < NLIMBS; ++i)
		aa[i] = a_ptr[i][elem_idx];

	Fp bb;
	for (u64 i = 0; i < NLIMBS; ++i)
		bb[i] = b_ptr[i][elem_idx];
#endif

	Fp rr;

	for (u64 i = 0; i < 10; ++i)
	{
		ptx_mul_mont_u64<Field>(rr, aa, bb, NoopIdx{});
		//mul_mont_coop_tile<Field>(tile, rr, aa, bb);
		aa = rr;
	}

#ifdef USE_LONG2
	for (u64 i = 0; i < NLIMBS / 2; ++i)
	{
		ulong2 val;
		val.x = rr[2 * i + 0];
		val.y = rr[2 * i + 1];
		a_ptr[i][elem_idx] = val;
	}
#else
	//a_ptr[my_offset] = rr;
	for (u64 i = 0; i < NLIMBS; ++i)
		a_ptr[i][n_elems] = rr[i];

#endif
}

//// ==== __global__ kernels (1 element per thread, strided) ====
//BLST_GLOBAL //__launch_bounds__(128, 2)
//	void kernel_mul_mont_384_u64
//	(
//		limb_t* BLST_restrict out,
//		const limb_t* BLST_restrict a,
//		const limb_t* BLST_restrict b,
//		size_t n_elems
//	)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	size_t stride = gridDim.x * blockDim.x;
//	for (size_t i = idx; i < n_elems; i += stride)
//	{
//		const limb_t* pa = a + 6 * i;
//		const limb_t* pb = b + 6 * i;
//		limb_t* pr = out + 6 * i;
//		Fp tt;
//		ptx_mul_mont_u64<Fp>(pr, pa, pb/*, (u64*)Fp::modulus().data, BLST_p0*/);
//		//ptx_mul_mont_u64<Fp>(tt.data, pa, pb);
//		//ptx_mul_mont_u64<Fp>(tt.data, tt.data, pb);
//		//ptx_mul_mont_u64<Fp>(tt.data, tt.data, pb);
//		//ptx_mul_mont_u64<Fp>(tt.data, tt.data, pb);
//		//ptx_mul_mont_u64<Fp>(tt.data, tt.data, pb);
//		//ptx_mul_mont_u64<Fp>(tt.data, tt.data, pb);
//		//ptx_mul_mont_u64<Fp>(pr, tt.data, pb);
//	}
//}


//struct GpuBuf64 {
//	limb_t* d_a = nullptr;
//	limb_t* d_b = nullptr;
//	limb_t* d_r = nullptr;
//	size_t n = 0;
//};
//void FpMain()
//{
//	const size_t n = 1 << 18;
//	GpuBuf64 buf;

//	// host data
//	/*std::vector<limb_t> h_a(6 * n);
//	std::vector<limb_t> h_b(6 * n);*/

//	//fill_host_u64(h_a);
//	//fill_host_u64(h_b);

//	// device alloc + copy
//	buf.n = n;
//	cudaMalloc(&buf.d_a, 6 * n * sizeof(limb_t));
//	cudaMalloc(&buf.d_b, 6 * n * sizeof(limb_t));
//	cudaMalloc(&buf.d_r, 6 * n * sizeof(limb_t));

//	/*CUDA_OK(cudaMemcpy(buf.d_a, h_a.data(), 6 * n * sizeof(limb_t), cudaMemcpyHostToDevice));
//	CUDA_OK(cudaMemcpy(buf.d_b, h_b.data(), 6 * n * sizeof(limb_t), cudaMemcpyHostToDevice));*/

//	dim3 grid, block;
//	block.x = 32;
//	grid.x = (n + block.x - 1) / block.x;
//	//choose_launch(n, grid, block);

//	for (u64 i = 0; i < 100; ++i)
//	{
//		kernel_mul_mont_384_u64 << <grid, block >> > (buf.d_r, buf.d_a, buf.d_b, n);
//	}

//}




void pp_Fp() {
	using Field = blst::Fp;
	constexpr int NLIMBS = 6;
	u64 n = 512;

	// Allocate ONE contiguous block for 12 * N * 6 limbs
	size_t total_limbs = n * NLIMBS;
	size_t bytes = total_limbs * sizeof(limb_t);

	limb_t* d_a;
	limb_t* d_b;
	cudaMalloc(&d_a, bytes);
	cudaMemset(d_a, 0, bytes); // Initialize to avoid NaN issues if any
	cudaMalloc(&d_b, bytes);
	cudaMemset(d_b, 0, bytes); // Initialize to avoid NaN issues if any


	int ops_per_warp = 32 / NLIMBS;
	int num_warps = ceil_div(n, ops_per_warp);
	int warps_per_block = 16;

	int num_blocks = ceil_div(num_warps, warps_per_block);
	int threads_per_block = warps_per_block * 32;

	auto begin = std::chrono::high_resolution_clock::now();
	for (u64 i = 0; i < 10; ++i)
	{
		bench_fp_optimized_kernel<Field> << <num_blocks, threads_per_block >> > (d_a, d_b, n);
		auto ec = cudaDeviceSynchronize();
		if (ec != cudaSuccess)
		{
			std::cout << "cuda " << cudaGetErrorString(ec) << std::endl;
			return;
		}
		ec = cudaGetLastError();
		if (ec != cudaSuccess)
		{
			std::cout << "cuda <<<" << num_blocks << "," << threads_per_block << ">>> "
				<< cudaGetErrorString(ec) << std::endl;
			return;
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - begin;
	printf("bench_fp_optimized_kernel time: %f ms\n", elapsed.count());
	printf("Throughput: %f ops/sec\n", n / (elapsed.count() / 1000.0));

	cudaFree(d_a);
	cudaFree(d_b);
}


void FpMain() {
	using Field = blst::Fp;
	constexpr int NLIMBS = 6;
	u64 n = 512;

	// Allocate ONE contiguous block for 12 * N * 6 limbs
	size_t total_limbs = n;
	size_t bytes = total_limbs * sizeof(limb_t);

	SoAFp d_a;
	SoAFp d_b;
	for (u64 i = 0; i < NLIMBS; ++i)
		cudaMalloc(&d_a[i], bytes);
	for (u64 i = 0; i < NLIMBS; ++i)
		cudaMalloc(&d_b[i], bytes);


	int ops_per_warp = 32;
	int num_warps = ceil_div(n, ops_per_warp);
	int warps_per_block = 16;

	int num_blocks = ceil_div(num_warps, warps_per_block);
	int threads_per_block = warps_per_block * 32;

	auto begin = std::chrono::high_resolution_clock::now();
	for (u64 i = 0; i < 10; ++i)
	{
		bench_fp_kernel<Field> << <num_blocks, threads_per_block >> > (d_a, d_b, n);
		cudaDeviceSynchronize();
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - begin;
	printf("bench_fp_kernel time: %f ms\n", elapsed.count());
	printf("Throughput: %f ops/sec\n", n / (elapsed.count() / 1000.0));

	for (u64 i = 0; i < NLIMBS; ++i)
		cudaFree(d_a[i]);
	for (u64 i = 0; i < NLIMBS; ++i)
		cudaFree(d_b[i]);
}

void pp_Fp12() {
	using Field = blst::Fp;
	constexpr int NLIMBS = 6;
	u64 n = 1 << 18;

	// Allocate ONE contiguous block for 12 * N * 6 limbs
	size_t total_limbs = 12 * n * NLIMBS;
	size_t bytes = total_limbs * sizeof(limb_t);

	limb_t* d_data;
	cudaMalloc(&d_data, bytes);
	cudaMemset(d_data, 0, bytes); // Initialize to avoid NaN issues if any

	int threads = 256;
	int tiles_per_warp = 32 / NLIMBS;
	int ops_per_block = (threads / 32) * tiles_per_warp;
	int grids = (n + ops_per_block - 1) / ops_per_block;
	auto begin = std::chrono::high_resolution_clock::now();
	//for (auto _ : state) 
	for (u64 i = 0; i < 10; ++i)
	{
		bench_fp12_optimized_kernel<Field> << <grids, threads >> > (d_data, n);
		auto ec = cudaDeviceSynchronize();
		if (ec != cudaSuccess)
		{
			std::cout << "cuda" << cudaGetErrorString(ec) << std::endl;
			return;
		}

	}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - begin;
	printf("bench_fp12_optimized_kernel time: %f ms\n", elapsed.count());
	printf("Throughput: %f ops/sec\n", n / (elapsed.count() / 1000.0));

	cudaFree(d_data);
}

void mainFp12_Triad_NoShuffle();


void user(int argc, char** argv)
{
	//FpMain();
	if (argc > 2 && argv[2] && strcmp(argv[2], "--avxbls") == 0)
	{
		std::cout << "--avxbls is intentionally disabled during the AVX backend refactor." << std::endl;
		return;
	}

	if (argc > 2 && argv[2] && strcmp(argv[2], "--soa") == 0)
	{
		Fp12Main();
		//FpMain();
		return;
	}


	if (argc > 2 && argv[2] && strcmp(argv[2], "--pp") == 0)
	{
		//pp_Fp12();
		//pp_Fp();
		mainFp12_Triad_NoShuffle();
		return;
	}

	if (argc > 2 && argv[2] && strcmp(argv[2], "--tile") == 0)
	{
		blstFp12TileMain();
		return;
	}
	if (argc > 2 && argv[2] && strcmp(argv[2], "--benchtile") == 0)
	{
		benchFp12TileMain();
		return;
	}


	std::cout << "use --pp, --soa or --tile" << std::endl;
	//G1Main();
}

#else

void user(int argc, char** argv)
{
	if (argc > 2 && argv[2] && strcmp(argv[2], "--avxbls") == 0) {
		std::cout << "--avxbls is intentionally disabled during the AVX backend refactor." << std::endl;
		return;
	}

	(void)argc;
	(void)argv;
	std::cout << "BLST_CUDA is disabled." << std::endl;
}

#endif


int main(int argc, char** argv) {

	if (argc > 1 && argv[1] && strcmp(argv[1], "--user") == 0) {
		user(argc, argv);
		return 0;
	}

	char arg0_default[] = "benchmark";
	char* args_default = arg0_default;
	if (!argv) {

		argc = 1;
		argv = &args_default;
	}
	::benchmark::Initialize(&argc, argv);
	if (::benchmark::ReportUnrecognizedArguments(argc, argv))
		return 1;
	::benchmark::RunSpecifiedBenchmarks();
	::benchmark::Shutdown();
	return 0;
}
