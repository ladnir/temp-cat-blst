#include <benchmark/benchmark.h>
#include <thrust/device_vector.h>
#include "blst/fields/Fp12.h"
#include "sqr_impl.h" 

using namespace blst;

namespace {

    // =========================================================
    // 1. Scalar Math (Unchanged)
    // =========================================================
    using RegFp = Fp;
    struct RegFp2 { RegFp c0, c1; };
    struct RegFp4 { RegFp2 c0, c1; };

    __device__ __forceinline__ void scalar_fp_add(RegFp& r, const RegFp& a, const RegFp& b, const limb_t* P) { Fp__add<Fp>(r, a, b, BasicIdx{ 0 }); }
    __device__ __forceinline__ void scalar_fp_sub(RegFp& r, const RegFp& a, const RegFp& b, const limb_t* P) { Fp__sub<Fp>(r, a, b, BasicIdx{ 0 }); }
    __device__ __forceinline__ void scalar_fp_mul(RegFp& r, const RegFp& a, const RegFp& b, const limb_t* P, uint64_t n0) { ptx_mul_mont_u64<Fp>(r, a, b, BasicIdx{ 0 }); }

    __device__ __forceinline__ void scalar_fp2_add(RegFp2& r, const RegFp2& a, const RegFp2& b, const limb_t* P) { scalar_fp_add(r.c0, a.c0, b.c0, P); scalar_fp_add(r.c1, a.c1, b.c1, P); }
    __device__ __forceinline__ void scalar_fp2_sub(RegFp2& r, const RegFp2& a, const RegFp2& b, const limb_t* P) { scalar_fp_sub(r.c0, a.c0, b.c0, P); scalar_fp_sub(r.c1, a.c1, b.c1, P); }
    __device__ __forceinline__ void scalar_fp2_mul_u1(RegFp2& r, const RegFp2& a, const limb_t* P) { RegFp t; scalar_fp_add(t, a.c0, a.c1, P); scalar_fp_sub(r.c0, a.c0, a.c1, P); r.c1 = t; }

    __device__ __forceinline__ void scalar_fp2_sqr(RegFp2& r, const RegFp2& a, const limb_t* P, uint64_t n0) {
        RegFp t0, t1, sum, r1;
        scalar_fp_mul(t0, a.c0, a.c0, P, n0);
        scalar_fp_mul(t1, a.c1, a.c1, P, n0);
        scalar_fp_add(sum, a.c0, a.c1, P);
        scalar_fp_mul(r1, sum, sum, P, n0);
        scalar_fp_sub(r1, r1, t0, P);
        scalar_fp_sub(r1, r1, t1, P);
        r.c1 = r1;
        scalar_fp_sub(r.c0, t0, t1, P);
    }

    __device__ __forceinline__ void scalar_sqr_fp4(RegFp4& ret, const RegFp2& a0, const RegFp2& a1, const limb_t* P, uint64_t n0) {
        RegFp2 t0, t1, sum, t1_mul;
        scalar_fp2_sqr(t0, a0, P, n0);
        scalar_fp2_sqr(t1, a1, P, n0);
        scalar_fp2_add(sum, a0, a1, P);
        scalar_fp2_sqr(ret.c1, sum, P, n0);
        scalar_fp2_sub(ret.c1, ret.c1, t0, P);
        scalar_fp2_sub(ret.c1, ret.c1, t1, P);
        scalar_fp2_mul_u1(t1_mul, t1, P);
        scalar_fp2_add(ret.c0, t0, t1_mul, P);
    }

    __device__ __forceinline__ void shfl_fp2(RegFp2& dst, const RegFp2& src, int src_lane) {
#pragma unroll
        for (int i = 0; i < 6; ++i) dst.c0.data[i] = __shfl_sync(0xFFFFFFFF, src.c0.data[i], src_lane);
#pragma unroll
        for (int i = 0; i < 6; ++i) dst.c1.data[i] = __shfl_sync(0xFFFFFFFF, src.c1.data[i], src_lane);
    }

    // =========================================================
    // 3. Coalesced Kernel
    // =========================================================

    template<typename Field>
    __launch_bounds__(256, 3)
        __global__ void bench_fp12_triad_coalesced_kernel(
            limb_t* base_ptr,
            const limb_t* P, uint64_t n0
        ) {
        // --- 1. Linear ID ---
        // Every thread corresponds to one Task (or is padding).
        int global_tid = blockIdx.x * blockDim.x + threadIdx.x;
        int lane = threadIdx.x & 31;

        // --- 2. Triad Logic ---
        // Even though memory is linear, logic depends on Rank (0,1,2)
        int triad_rank = lane % 3;
        int triad_idx = lane / 3;

        // Skip padding threads (30, 31)
        if (lane >= 30) return;

        // SoA Layout: 24 Arrays (12 input limbs for A, 12 input limbs for B).
        // Stride is TOTAL_THREADS (including padding).
        // This ensures thread T reads address X, thread T+1 reads address X+1.
        size_t grid_stride = gridDim.x * blockDim.x;

        // --- 3. Coalesced Load ---
        RegFp2 in_A, in_B;

        // Loop unrolled over the 6 limbs of the Fp
        // input_A is limbs 0..5. input_B is limbs 6..11.
        // Each limb is in its own array separated by 'grid_stride'.

        limb_t* my_ptr = base_ptr + global_tid;

#pragma unroll
        for (int i = 0; i < 6; ++i) in_A.c0.data[i] = my_ptr[(0 * 6 + i) * grid_stride];
#pragma unroll
        for (int i = 0; i < 6; ++i) in_A.c1.data[i] = my_ptr[(1 * 6 + i) * grid_stride];

#pragma unroll
        for (int i = 0; i < 6; ++i) in_B.c0.data[i] = my_ptr[(2 * 6 + i) * grid_stride];
#pragma unroll
        for (int i = 0; i < 6; ++i) in_B.c1.data[i] = my_ptr[(3 * 6 + i) * grid_stride];

        // --- 4. Compute ---
        RegFp4 T_local;
        scalar_sqr_fp4(T_local, in_A, in_B, P, n0);

        // --- 5. Reconstruction Swap (Warp Shuffle) ---
        // Rank 1 needs data from Rank 2. Rank 2 needs data from Rank 1.
        // Base lane for this triad
        int base_lane = triad_idx * 3;

        RegFp2 rec_A = in_A;
        RegFp2 rec_B = in_B;

        if (triad_rank == 1) {
            shfl_fp2(rec_A, in_A, base_lane + 2); // Get a01 from R2
            shfl_fp2(rec_B, in_B, base_lane + 2); // Get a12 from R2
        }
        else if (triad_rank == 2) {
            shfl_fp2(rec_A, in_B, base_lane + 1); // Get a02 from R1
            shfl_fp2(rec_B, in_A, base_lane + 1); // Get a10 from R1
        }

        // --- 6. Reconstruct ---
        RegFp2 res_A, res_B, tmp;

        // Result A = 3*T.c0 - 2*rec_A
        scalar_fp2_sub(tmp, T_local.c0, rec_A, P);
        scalar_fp2_add(tmp, tmp, T_local.c0, P);
        scalar_fp2_add(res_A, tmp, tmp, P);

        // Result B = 3*T.c1 + 2*rec_B (with mul_u1 for rank 2)
        RegFp2 T_partB = T_local.c1;
        if (triad_rank == 2) scalar_fp2_mul_u1(T_partB, T_partB, P);

        scalar_fp2_add(tmp, T_partB, rec_B, P);
        scalar_fp2_add(tmp, tmp, T_partB, P);
        scalar_fp2_add(res_B, tmp, tmp, P);

        // --- 7. Coalesced Store ---
#pragma unroll
        for (int i = 0; i < 6; ++i) my_ptr[(0 * 6 + i) * grid_stride] = res_A.c0.data[i];
#pragma unroll
        for (int i = 0; i < 6; ++i) my_ptr[(1 * 6 + i) * grid_stride] = res_A.c1.data[i];

#pragma unroll
        for (int i = 0; i < 6; ++i) my_ptr[(2 * 6 + i) * grid_stride] = res_B.c0.data[i];
#pragma unroll
        for (int i = 0; i < 6; ++i) my_ptr[(3 * 6 + i) * grid_stride] = res_B.c1.data[i];
    }

    // =========================================================
    // 4. Host Repacking Helpers
    // =========================================================

    // Fp12 indices for tasks
    // T0 (Rank0): Inputs a00(0), a11(4).
    // T1 (Rank1): Inputs a10(3), a02(2).
    // T2 (Rank2): Inputs a01(1), a12(5).

    //// We iterate over input Fp12s and scatter them into the padded buffer.
    //void repack_input(
    //    const std::vector<Fp12>& h_in,
    //    std::vector<limb_t>& h_out,
    //    size_t n_elems, size_t total_threads
    //) {
    //    // Map Task->Input Index
    //    int map_A[] = { 0, 3, 1 }; // Rank 0->a00, Rank 1->a10, Rank 2->a01
    //    int map_B[] = { 4, 2, 5 }; // Rank 0->a11, Rank 1->a02, Rank 2->a12

    //    // We fill h_out as if it were on GPU: SoA with stride = total_threads.
    //    // h_out has size: 24 * total_threads
    //    // For each Fp12 'k', we write to 3 threads.

    //    for (size_t k = 0; k < n_elems; ++k) {
    //        // Which warp/triad is this?
    //        size_t triad_global = k;
    //        size_t warp_idx = triad_global / 10;
    //        size_t triad_in_warp = triad_global % 10;

    //        // Base thread ID for this triad
    //        size_t base_tid = warp_idx * 32 + triad_in_warp * 3;

    //        for (int r = 0; r < 3; ++r) {
    //            size_t tid = base_tid + r;

    //            // Get the Fp2 input pair for this rank
    //            const Fp2& fpA = (map_A[r] < 3) ? h_in[k].data[0].data[map_A[r]] : h_in[k].data[1].data[map_A[r] - 3];
    //            const Fp2& fpB = (map_B[r] < 3) ? h_in[k].data[0].data[map_B[r]] : h_in[k].data[1].data[map_B[r] - 3];

    //            // Write to h_out [Limb][TID]
    //            for (int i = 0; i < 6; ++i) {
    //                h_out[(0 * 6 + i) * total_threads + tid] = fpA.data[0][i]; // in_A.c0
    //                h_out[(1 * 6 + i) * total_threads + tid] = fpA.data[1][i]; // in_A.c1
    //                h_out[(2 * 6 + i) * total_threads + tid] = fpB.data[0][i]; // in_B.c0
    //                h_out[(3 * 6 + i) * total_threads + tid] = fpB.data[1][i]; // in_B.c1
    //            }
    //        }
    //    }
    //}

    void BenchFp12_Triad_Coalesced(benchmark::State& state) {
        using Field = blst::Fp;
        constexpr int NLIMBS = 6;
        size_t n = 1 << 16;

        // 1. Calculate Padded Size
        // 10 elems per warp (32 threads).
        size_t warps_needed = (n + 9) / 10;
        size_t total_threads = warps_needed * 32;
        size_t buffer_limbs = total_threads * 24; // 24 limbs per thread (2 input Fp2s)

        // 2. Alloc
        limb_t* d_data;
        cudaMalloc(&d_data, buffer_limbs * sizeof(limb_t));
        // Note: For real bench, we'd repack random data here. 
        // For pure throughput bench, uninitialized/zero is fine to measure cycles.
        cudaMemset(d_data, 0, buffer_limbs * sizeof(limb_t));

        limb_t* d_P;
        cudaMalloc(&d_P, NLIMBS * sizeof(limb_t));
        cudaMemcpy(d_P, Field::modulus().data, NLIMBS * sizeof(limb_t), cudaMemcpyHostToDevice);
        uint64_t n0 = Field::montgomery_constant();

        int threads = 256; // 8 Warps
        int blocks = (total_threads + threads - 1) / threads;

        for (auto _ : state) {
            bench_fp12_triad_coalesced_kernel<Field> << <blocks, threads >> > (d_data, d_P, n0);
            cudaDeviceSynchronize();
        }

        cudaFree(d_data);
        cudaFree(d_P);
    }
}

BENCHMARK(BenchFp12_Triad_Coalesced)->Unit(benchmark::kMillisecond);

void mainFp12_Triad_NoShuffle() {
    //using Field = blst::Fp;
    //constexpr int NLIMBS = 6;
    //size_t n = 1 << 18;

    //size_t total_limbs = 12 * n * NLIMBS;
    //size_t bytes = total_limbs * sizeof(limb_t);

    //limb_t* d_data;
    //cudaMalloc(&d_data, bytes);
    //cudaMemset(d_data, 0, bytes);

    //limb_t* d_P;
    //cudaMalloc(&d_P, NLIMBS * sizeof(limb_t));
    //cudaMemcpy(d_P, Field::modulus().data, NLIMBS * sizeof(limb_t), cudaMemcpyHostToDevice);
    //uint64_t n0 = Field::montgomery_constant();

    //int threads = 256;
    //int tiles_per_warp = 32 / 3; // 10 Triads per warp
    //int ops_per_block = (threads / 32) * tiles_per_warp;
    //int grids = (n + ops_per_block - 1) / ops_per_block;

    //for (int ii = 0; ii < 10; ++ii) 
    //{
    //    bench_fp12_triad_noshuffle_kernel<Field> << <grids, threads >> > (d_data, n, d_P, n0);
    //    cudaDeviceSynchronize();
    //}

    //cudaFree(d_data);
    //cudaFree(d_P);
    using Field = blst::Fp;
    constexpr int NLIMBS = 6;
    size_t n = 1 << 18;

    // 1. Calculate Padded Size
    // 10 elems per warp (32 threads).
    size_t warps_needed = (n + 9) / 10;
    size_t total_threads = warps_needed * 32;
    size_t buffer_limbs = total_threads * 24; // 24 limbs per thread (2 input Fp2s)

    // 2. Alloc
    limb_t* d_data;
    cudaMalloc(&d_data, buffer_limbs * sizeof(limb_t));
    // Note: For real bench, we'd repack random data here. 
    // For pure throughput bench, uninitialized/zero is fine to measure cycles.
    cudaMemset(d_data, 0, buffer_limbs * sizeof(limb_t));

    limb_t* d_P;
    cudaMalloc(&d_P, NLIMBS * sizeof(limb_t));
    cudaMemcpy(d_P, Field::modulus().data, NLIMBS * sizeof(limb_t), cudaMemcpyHostToDevice);
    uint64_t n0 = Field::montgomery_constant();

    int threads = 256; // 8 Warps
    int blocks = (total_threads + threads - 1) / threads;

    for (int ii = 0; ii < 10; ++ii) 
    {
        bench_fp12_triad_coalesced_kernel<Field> << <blocks, threads >> > (d_data, d_P, n0);
        cudaDeviceSynchronize();
    }

    cudaFree(d_data);
    cudaFree(d_P);
}