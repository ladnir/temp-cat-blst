#pragma once

#include <benchmark/benchmark.h>
#include <thrust/device_vector.h>
#include "blst/fields/Fp12.h"
#include "blst/backend/cuda/mul_mont_cuda.h"
#include "blst/cuda/fields/parallel_prefix.h"
using namespace blst;

namespace {
    // =========================================================
    // 1. Data Structures
    // =========================================================

    // Register Resident Types (Same as before)
    struct CoopFp2 { limb_t c0; limb_t c1; };
    struct CoopFp12 { CoopFp2 c[6]; };

    // =========================================================
    // 2. Arithmetic Wrappers (Same as before)
    // =========================================================
    // ... [Copy fp2_add, fp2_sub, fp2_mul_u1, fp2_sqr, sqr_fp4 from previous code] ...
    // ... [They are unchanged] ...

    template<typename Field, int N>
    __device__ __forceinline__ void fp2_add(Tile<N> t, CoopFp2& r, const CoopFp2& a, const CoopFp2& b) {
        pp_add_coop_tile_aos<Field, N>(t, a.c0, b.c0, r.c0);
        pp_add_coop_tile_aos<Field, N>(t, a.c1, b.c1, r.c1);
    }

    template<typename Field, int N>
    __device__ __forceinline__ void fp2_sub(Tile<N> t, CoopFp2& r, const CoopFp2& a, const CoopFp2& b) {
        pp_sub_coop_tile_aos<Field, N>(t, a.c0, b.c0, r.c0);
        pp_sub_coop_tile_aos<Field, N>(t, a.c1, b.c1, r.c1);
    }

    template<typename Field, int N>
    __device__ __forceinline__ void fp2_mul_u1(Tile<N> t, CoopFp2& r, const CoopFp2& a) {
        limb_t t_val;
        pp_add_coop_tile_aos<Field, N>(t, a.c0, a.c1, t_val);
        pp_sub_coop_tile_aos<Field, N>(t, a.c0, a.c1, r.c0);
        r.c1 = t_val;
    }

    template<typename Field, int N>
    __device__ __forceinline__ void fp2_sqr(Tile<N> t, CoopFp2& r, const CoopFp2& a) {
        limb_t t0, t1, sum, r1;
        mul_mont_coop_tile<Field>(t, t0, a.c0, a.c0);
        mul_mont_coop_tile<Field>(t, t1, a.c1, a.c1);
        pp_add_coop_tile_aos<Field, N>(t, a.c0, a.c1, sum);
        mul_mont_coop_tile<Field>(t, r1, sum, sum);
        pp_sub_coop_tile_aos<Field, N>(t, r1, t0, r1);
        pp_sub_coop_tile_aos<Field, N>(t, r1, t1, r1);
        r.c1 = r1;
        pp_sub_coop_tile_aos<Field, N>(t, t0, t1, r.c0);
    }

    template<typename Field, int N>
    __device__ __forceinline__ void sqr_fp4(Tile<N> t, CoopFp2& r0, CoopFp2& r1, const CoopFp2& a0, const CoopFp2& a1) {
        CoopFp2 t0, t1, sum, t1_mul;
        fp2_sqr<Field,N>(t, t0, a0);
        fp2_sqr<Field,N>(t, t1, a1);
        fp2_add<Field,N>(t, sum, a0, a1);
        fp2_sqr<Field,N>(t, r1, sum);
        fp2_sub<Field,N>(t, r1, r1, t0);
        fp2_sub<Field,N>(t, r1, r1, t1);
        fp2_mul_u1<Field,N>(t, t1_mul, t1);
        fp2_add<Field,N>(t, r0, t0, t1_mul);
    }

    template<typename Field, int N>
    __device__ __forceinline__ void cyclotomic_sqr_slim(Tile<N> t, CoopFp12& val) {
        CoopFp2 t0_0, t0_1; 
        CoopFp2 t1_0, t1_1; 
        CoopFp2 t2_0, t2_1; 

        sqr_fp4<Field,N>(t, t0_0, t0_1, val.c[0], val.c[4]);
        sqr_fp4<Field,N>(t, t1_0, t1_1, val.c[3], val.c[2]);
        sqr_fp4<Field,N>(t, t2_0, t2_1, val.c[1], val.c[5]);

        CoopFp2 tmp;

        // r00
        fp2_sub<Field,N>(t, tmp, t0_0, val.c[0]);
        fp2_add<Field,N>(t, tmp, tmp, t0_0);
        fp2_add<Field,N>(t, val.c[0], tmp, tmp);

        // r11
        fp2_add<Field,N>(t, tmp, t0_1, val.c[4]);
        fp2_add<Field,N>(t, tmp, tmp, t0_1);
        fp2_add<Field,N>(t, val.c[4], tmp, tmp);

        // r01
        fp2_sub<Field,N>(t, tmp, t1_0, val.c[1]);
        fp2_add<Field,N>(t, tmp, tmp, t1_0);
        fp2_add<Field,N>(t, val.c[1], tmp, tmp);

        // r12
        fp2_add<Field,N>(t, tmp, t1_1, val.c[5]);
        fp2_add<Field,N>(t, tmp, tmp, t1_1);
        fp2_add<Field,N>(t, val.c[5], tmp, tmp);

        // r02
        fp2_sub<Field,N>(t, tmp, t2_0, val.c[2]);
        fp2_add<Field,N>(t, tmp, tmp, t2_0);
        fp2_add<Field,N>(t, val.c[2], tmp, tmp);

        // r10
        CoopFp2 t2_1_u;
        fp2_mul_u1<Field,N>(t, t2_1_u, t2_1);
        fp2_add<Field,N>(t, tmp, t2_1_u, val.c[3]);
        fp2_add<Field,N>(t, tmp, tmp, t2_1_u);
        fp2_add<Field,N>(t, val.c[3], tmp, tmp);
    }

    // =========================================================
    // 4. Kernel (Optimization Applied)
    // =========================================================

    // Optimization: Single Base Pointer + Launch Bounds
    // Forces compiler to use fewer registers, increasing occupancy.
    // 256 threads, min 4 blocks => max 64 registers/thread
    template<typename Field>
    __launch_bounds__(256, 3)
    __global__ void bench_fp12_optimized_kernel(
        limb_t* __restrict__ base_ptr, 
        uint64_t n_elems
    ) {
        constexpr int NLIMBS = Field::num_limb();
        Tile<NLIMBS> tile;
        if (!tile.is_full()) return;

        int warp_id = blockIdx.x * (blockDim.x / 32) + (threadIdx.x / 32);
        int elem_idx = warp_id * tile.tiles_per_warp() + tile.mGid;

        if (elem_idx >= n_elems) return;

        CoopFp12 reg;
        int rank = tile.thread_rank();
        
        // Stride: Each coordinate (c00, c01...) is separated by (N * 6) limbs
        // Or simpler: Flat SoA.
        // Array 0: c00 for all N elems.
        // Array 1: c01 for all N elems.
        // Offset between arrays = n_elems * NLIMBS
        size_t array_stride = n_elems * NLIMBS;
        size_t my_offset = elem_idx * NLIMBS + rank;

        // 1. Coalesced Load (Unrolled)
        #pragma unroll
        for(int k=0; k<6; ++k) {
            reg.c[k].c0 = base_ptr[ (2*k)   * array_stride + my_offset ];
            reg.c[k].c1 = base_ptr[ (2*k+1) * array_stride + my_offset ];
        }

        // 2. Compute
        cyclotomic_sqr_slim<Field, NLIMBS>(tile, reg);

        // 3. Coalesced Store
        #pragma unroll
        for(int k=0; k<6; ++k) {
            base_ptr[ (2*k)   * array_stride + my_offset ] = reg.c[k].c0;
            base_ptr[ (2*k+1) * array_stride + my_offset ] = reg.c[k].c1;
        }
    }


    // Optimization: Single Base Pointer + Launch Bounds
    // Forces compiler to use fewer registers, increasing occupancy.
    // 256 threads, min 4 blocks => max 64 registers/thread
    template<typename Field>
        __global__ void bench_fp_optimized_kernel(
            limb_t* __restrict__ a_ptr,
            limb_t* __restrict__ b_ptr,
            uint64_t n_elems
        ) {
        constexpr int NLIMBS = Field::num_limb();
        Tile<NLIMBS> tile;
        if (!tile.is_full()) return;

        int warp_id = blockIdx.x * (blockDim.x / 32) + (threadIdx.x / 32);
        int elem_idx = warp_id * tile.tiles_per_warp() + tile.mGid;

        if (elem_idx >= n_elems) return;

        int rank = tile.thread_rank();

        // Stride: Each coordinate (c00, c01...) is separated by (N * 6) limbs
        // Or simpler: Flat SoA.
        // Array 0: c00 for all N elems.
        // Array 1: c01 for all N elems.
        // Offset between arrays = n_elems * NLIMBS
        size_t my_offset = elem_idx * NLIMBS + rank;

		limb_t aa = a_ptr[my_offset];
		limb_t bb = b_ptr[my_offset];
		limb_t rr;

        for(u64 i = 0; i < 10; ++i)
        {
		    mul_mont_coop_tile<Field>(tile, rr, aa, bb);
			aa = rr;
        }

		a_ptr[my_offset] = rr;
    }



    // =========================================================
    // 5. Benchmark Harness
    // =========================================================

    void BenchFp12_Optimized(benchmark::State& state) {
        using Field = blst::Fp;
        constexpr int NLIMBS = 6;
        size_t n = 1 << 16; 
        
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

        for (auto _ : state) {
            bench_fp12_optimized_kernel<Field><<<grids, threads>>>(d_data, n);
            cudaDeviceSynchronize();
        }

        cudaFree(d_data);
    }
}

BENCHMARK(BenchFp12_Optimized)->Unit(benchmark::kMillisecond);;