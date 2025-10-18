#pragma once

#include "multi_scalar.h"
#include <cuda_runtime.h>

#ifndef __INTELLISENSE__
#include <thrust/iterator/iterator_facade.h>

#include <cub/cub.cuh>
#include <cub/device/device_reduce.cuh >
#include <thrust/iterator/transform_iterator.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/logical.h>
#include <thrust/functional.h>
#endif
#include <cstdint>
#include <vector>
#include "blst/common/cuda_launch.h"

namespace blst
{

    template<typename T>
    using d2Span = cuda::std::mdspan<T, cuda::std::dextents<size_t, 2>>;

    template<typename T>
    using d3Span = cuda::std::mdspan<T, cuda::std::dextents<size_t, 3>>;


    // ---- Packed val: (idx<<1)|sign --------------------------------------
    struct IdxSign {
        uint32_t v;
        // the index of the point in the input array that this value corresponds to
        BLST_HD BLST_FORCEINLINE uint32_t idx()  const { return v >> 1; }

        // the sign bit: 0 = positive, 1 = negative
        BLST_HD BLST_FORCEINLINE uint32_t sign() const { return v & 1; }
    };

    // ========== Kernels ==========

    // taka flat array of scalars, extract the w-bit Booth encoded digits
    // the output order is column-major by digit. First row has the N LSB digits, ...
    BLST_GLOBAL
        void extract_digits_flat_kernel(const uint8_t* BLST_restrict scalars,
            int nbytes,           // bytes per scalar
            int N,                // number of scalars
            int nbits,            // significant bits per scalar
            int wbits,            // window size (e.g., 5). W=2^wbits
            int D,                // number of windows = ceil(nbits / wbits)
            u8* BLST_restrict digits,
            IdxSign* BLST_restrict index_signs)   // [0,N]|{0,1}^N*D
    {
        int tid = blockDim.x * blockIdx.x + threadIdx.x;
        int total = N * D;
        if (tid >= total)
            return;

        int i = tid / D;     // scalar index
        int d = tid % D;     // window index (0 = LSB window)

        auto outIdx = d * N + i;

        const uint8_t* s = scalars + i * nbytes;
        uint32_t digit = get_wval_limb(s, d * wbits, wbits);
        auto encode = booth_encode(digit, wbits);
        digits[outIdx] = encode & ((1u << wbits) - 1u);
        index_signs[outIdx] = IdxSign{ uint32_t((encode >> wbits) & 1u | (i << 1)) };
    }

    // taka flat array of scalars, extract the w-bit Booth encoded digits
    // the output order is column-major by digit. First row has the N LSB digits, ...
    cudaError extract_digits_flat(
        const uint8_t* BLST_restrict scalars,
        int nbytes,           // bytes per scalar
        int N,                // number of scalars
        int nbits,            // significant bits per scalar
        int wbits,            // window size (e.g., 5). W=2^wbits
        int D,                // number of windows = ceil(nbits / wbits)
        u8* BLST_restrict digit,	// {0,1}^N*D
        IdxSign* BLST_restrict idx_sign)   // [0,N]|{0,1}^N*D
    {
        if (N <= 0 || D <= 0)
            return cudaSuccess;

        int total = N * D;
        int blockSize = 256;
        int numBlocks = (total + blockSize - 1) / blockSize;
        extract_digits_flat_kernel << < numBlocks, blockSize >> > (
            scalars,
            nbytes,
            N,
            nbits,
            wbits,
            D,
            digit,
            idx_sign);

        return cudaGetLastError();
    }

    static inline size_t align_up(size_t n, size_t a)
    {
        return (n + (a - 1)) & ~(a - 1);
    }

    // Sort each of D digit slices of N elements with a single segmented sort.
	// given D digit arrays of N elements each, sort each digit array
	// we will sort key value pairs, where the key is the digit value (u8)
	// The value is the index of the point in the input array that this digit
	// corresponds to. The LSB of the value is reserved for the sign bit.
    //
    cudaError_t sort_digits_segmented(
        void* d_workspace,
        size_t& workspace_bytes,
        const u8* d_keys_input,          // [D*N], column-major by digit
        u8* d_keys_output,          // [D*N], column-major by digit
        const IdxSign* d_vals_input,          // [D*N]
        IdxSign* d_vals_output,          // [D*N]
        int          N,
        int          D,
        int          wbits,
        cudaStream_t stream = 0)
    {
        if (N <= 0 || D <= 0) return cudaSuccess;
        if (wbits < 0) return cudaErrorInvalidValue;
        if (wbits > 8) return cudaErrorInvalidValue;

        size_t count = static_cast<size_t>(N) * static_cast<size_t>(D);
        auto d_begins = thrust::make_transform_iterator(
            thrust::make_counting_iterator<u32>(0),
            [=] __host__ __device__(u32 i) { return i * N; }
        );
        auto d_ends = d_begins + 1;


        // 1) Query CUB temp storage
        size_t required = 0;
        cub::DeviceSegmentedRadixSort::SortPairs(
            nullptr,					       // void *d_temp_storage,
            required,						   // size_t &temp_storage_bytes,
            d_keys_input,					   // const KeyT *d_keys_in,
            d_keys_output,					   // KeyT *d_keys_out,
            d_vals_input,					   // const ValueT *d_values_in,
            d_vals_output,					   // ValueT *d_values_out,
            count,							   // ::cuda::std::int64_t num_items,
            D,								   // ::cuda::std::int64_t num_segments,
            d_begins ,			               // BeginOffsetIteratorT d_begin_offsets,
            d_ends, 			               // EndOffsetIteratorT d_end_offsets,
            0,								   // int begin_bit = 0,
            wbits,						       // int end_bit = sizeof(KeyT) * 8,
            stream							   // cudaStream_t stream = 0,
        );

        if (d_workspace == nullptr) {
            workspace_bytes = required;
            return cudaSuccess;
        }
        if (workspace_bytes < required) 
            return cudaErrorInvalidValue;

        // 2) Run segmented sort
        auto ec = cub::DeviceSegmentedRadixSort::SortPairs(
            d_workspace,
            workspace_bytes,
            d_keys_input,
            d_keys_output,
            d_vals_input,
            d_vals_output,
            count,
            D,
            d_begins,
            d_ends,
            0,
            wbits,
            stream);
        if (ec != cudaSuccess) 
            return ec;

        return cudaSuccess;
    }


    // cub will call operations on default-constructed Pointxyzz.
    // so we need to ensure that the default constructor is safe.
    // we also require alignment for Pointxyzz.
    template<typename Field>
    struct DefaultedPointxyzz : public Pointxyzz<Field> {

        BLST_DEVICE BLST_FORCEINLINE DefaultedPointxyzz()
        {
            memset(this, 0, sizeof(*this));
        }
        DefaultedPointxyzz(const DefaultedPointxyzz& other) = default;
        DefaultedPointxyzz& operator=(const DefaultedPointxyzz& other) = default;

    };


    // ---- Transform functor: packed -> Pointxyzz by gathering -------------
    template<class Field>
    struct PackedToPointXYZZ {
        const AffinePoint<Field>* BLST_restrict points; // device
        int num_points;
        BLST_DEVICE BLST_FORCEINLINE
            DefaultedPointxyzz<Field> operator()(IdxSign index_sign) const {

            assert(index_sign.idx() < num_points);
            AffinePoint<Field> p = points[index_sign.idx()];
            AffinePoint<Field>::cneg(p, index_sign.sign());
            DefaultedPointxyzz<Field> q;
            DefaultedPointxyzz<Field>::from_affine(q, p);
            return q;
        }
    };

    // ---- Point monoid for reduces/scans ----------------------------------
    // add: DefaultedPointxyzz x DefaultedPointxyzz -> DefaultedPointxyzz
    template<class Field>
    struct XYZZAdd {
        BLST_DEVICE BLST_FORCEINLINE
            DefaultedPointxyzz<Field> operator()(
                const DefaultedPointxyzz<Field>& a, 
                const DefaultedPointxyzz<Field>& b) const 
        {
            DefaultedPointxyzz<Field> r;
            DefaultedPointxyzz<Field>::dadd(r, a, b);
            return r;
        }
    };

    // Helpers: align up sizes / pointers
    BLST_HD BLST_INLINE static size_t align_up_size(
        size_t x,
        size_t a)
    {
        size_t mask = a - 1;
        return (x + mask) & ~mask;
    }

    BLST_HD BLST_INLINE static uint8_t* align_up_ptr(
        uint8_t* p,
        size_t a)
    {
        uintptr_t u = reinterpret_cast<uintptr_t>(p);
        uintptr_t aligned = (u + (a - 1)) & ~(uintptr_t)(a - 1);
        return reinterpret_cast<uint8_t*>(aligned);
    }

    // for each input key,val pair, scatter the val to d_groups[key-1].
    // we will assume that d_groups has been zeroed out before this kernel is launched.
    // any missing groups will remain at infinity (0).
    // we assume that keys are in [0..max_key], and that group 0 is skipped.
    template<class Field>
    BLST_GLOBAL void scatter_group_sums(
        const u16* d_group_keys, // we have num_keys entries, each in [0..max_key]
        const DefaultedPointxyzz<Field>* d_group_vals, // we have num_keys entries
        u64 num_keys,
        u64 max_key,
        Pointxyzz<Field>* d_groups)
    {
        u64 tid = blockDim.x * blockIdx.x + threadIdx.x;
        if (tid >= num_keys)
            return;
        u16 k = d_group_keys[tid];
        if (k == 0)
            return;
        if (k > max_key)
        {
            printf("%d > %d", (int)k, (int)max_key);
            assert(0 && "out of bounds key"); // should not happen
            return;
        }
        DefaultedPointxyzz<Field> acc = d_group_vals[tid];
        // scatter to group k-1
        d_groups[k - 1] = acc;

    }


    template<class Field, typename SrcIter>
    cudaError_t reduce_groups(
        void* d_workspace,
        size_t& workspace_bytes,
        const u16* d_keys,
        SrcIter&& src_iter,
        Pointxyzz<Field>* d_groups,
        int N,
        int W,
        cudaStream_t stream = 0
    ) {
        using XY = DefaultedPointxyzz<Field>;

        size_t group_count = W + 1;

#ifndef NDEBUG
        // make sure all keys are in range. We will launch a kernel to check this
        {
            bool keys_ok =
                thrust::all_of(
                    thrust::cuda::par.on(stream),
                    d_keys,
                    d_keys + N,
                    [=] BLST_DEVICE(u16 k)
            {
                return k <= static_cast<u16>(W);
            }
                );

            if (!keys_ok)
            {
                return cudaErrorInvalidValue;
            }
        }
#endif

        size_t temp_storage_bytes = 0;
        {
            //auto gather_iter = thrust::make_transform_iterator(
            //	d_point_index,
            //	PackedToPointXYZZ<Field>{ d_points, d_signs, 0 }
            //);

            u16* dummy_keys = nullptr;
            XY* dummy_vals = nullptr;
            size_t* dummy_size_ptr = nullptr;

            cudaError_t err = cub::DeviceReduce::ReduceByKey(
                nullptr,
                temp_storage_bytes,
                d_keys,
                dummy_keys,
                src_iter,
                dummy_vals,
                dummy_size_ptr,
                XYZZAdd<Field>{},
                N,
                stream
            );
            if (err != cudaSuccess) {
                return err;
            }

            size_t required = 0;
            required += temp_storage_bytes;
            required = align_up_size(required, alignof(XY));
            required += sizeof(XY) * group_count;

            required = align_up_size(required, alignof(u16));
            required += sizeof(u16) * group_count;

            required = align_up_size(required, alignof(size_t));
            required += sizeof(size_t);

            if (d_workspace == nullptr) {
                workspace_bytes = required;
                return cudaSuccess;
            }

            if (workspace_bytes < required) {
                workspace_bytes = required;
                return cudaErrorInvalidValue;
            }
        }

        uint8_t* ws = static_cast<uint8_t*>(d_workspace);

        void* d_temp_storage = ws;
        ws += temp_storage_bytes;

        ws = align_up_ptr(ws, alignof(XY));
        XY* d_group_vals = reinterpret_cast<XY*>(ws);
        ws += sizeof(XY) * group_count;

        ws = align_up_ptr(ws, alignof(u16));
        u16* d_group_keys = reinterpret_cast<u16*>(ws);
        ws += sizeof(u16) * group_count;
        cudaMemsetAsync(d_group_keys, 0, sizeof(u16) * group_count, stream);

        ws = align_up_ptr(ws, alignof(u64));
        size_t* d_size = (size_t*)ws;
        ws += sizeof(size_t);



        cudaError_t err = cub::DeviceReduce::ReduceByKey(
            d_temp_storage,
            temp_storage_bytes,
            d_keys,
            d_group_keys,
            src_iter,
            d_group_vals,
            d_size,
            XYZZAdd<Field>{},
            N,
            stream
        );
        if (err != cudaSuccess) {
            return err;
        }

        int blk = 256;
        int grd = static_cast<int>((W + 1 + blk - 1) / blk);

        // copy the groups to output, skipping group 0.
        // memset so that missing groups remain at infinity.
        cudaMemsetAsync(d_groups, 0, sizeof(Point<Field>) * W, stream);
        scatter_group_sums<Field> << <grd, blk, 0, stream >> > (
            d_group_keys,
            d_group_vals,
            W + 1,
            static_cast<u64>(W),
            d_groups);

        return cudaGetLastError();
    }

    template<class Field, class ValT>
    __global__ void scatter_group_sums_multi(
        const uint32_t* __restrict d_group_keys,
        const ValT* __restrict d_group_vals,
        const size_t* __restrict d_num_runs_ptr,
        int W,
        int M,
        Pointxyzz<Field>* __restrict d_groups)
    {
        size_t num_runs = *d_num_runs_ptr;
        size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
        if (tid >= num_runs)
            return;

        uint32_t gk = d_group_keys[tid];
        uint32_t wspan = static_cast<uint32_t>(W + 1);
        uint32_t m = gk / wspan;
        uint32_t k = gk % wspan;

        if (k == 0)
            return;

        size_t out_idx = static_cast<size_t>(m) * static_cast<size_t>(W) + static_cast<size_t>(k - 1);
        d_groups[out_idx] = static_cast<Pointxyzz<Field>>(d_group_vals[tid]);
    }

    // ===== Per-group reduce: keys/vals sorted by key =======================
    // 
    // Takes an input M vectors of size N, each with a key in [0..W].
    // Groups are defined by keys. Key 0 is skipped.
    // For each M, output W groups, each being the sum of all points
    // within the respective group.
    // 
    // @d_workspace: temp buffer, or nullptr to query size c
    // @workspace_bytes: in: size of d_workspace; out: required size if d_workspace is nullptr or too small
    // @d_keys: size M * N matrix iter with values in [0..W], in grouped order
    // @d_vals: size M * N matrix iter of points to be reduced.
    // @d_groups: size M * W vec of Point<Field>, the output groups, skipping group 0.
    // @M: number of independent vectors (digit vectors)
    // @N: number of keys/vals/points
    // @W: number of groups = 2^(wbits-1). keys are in [0..W], W+1 possible values. 0 is skipped.
    // @stream: cuda stream
    // 
    // This algorithm works as follows:
    //   1) load the points: we have N vals which are packed (idx<<1)|sign,
    //      for each val, we load points[idx] and apply sign to get DefaultedPointxyzz.
    //   2) for each group (key), we reduce the points in that group to get 
    //      the group sum. This will be done via  thrust::reduce_by_key.
    //   3) We scatter these to the output groups array. we will ignore group 0 
    //      and fill any missing groups with the point at infinity.
    //
    template<class Field, typename SrcIter>
    cudaError_t reduce_groups(
        void* d_workspace,
        size_t& workspace_bytes,
        const u16* d_keys,
        SrcIter&& src_iter,
        Pointxyzz<Field>* d_groups,
        int M,
        int N,
        int W,
        cudaStream_t stream)
    {
        using SrcItDecayed = typename std::decay<SrcIter>::type;
        using ValT = typename thrust::iterator_value<SrcItDecayed>::type;

        size_t total = static_cast<size_t>(M) * static_cast<size_t>(N);
        size_t group_count = static_cast<size_t>(M) * static_cast<size_t>(W + 1);

#ifndef NDEBUG
        {
            bool keys_ok =
                thrust::all_of(
                    thrust::cuda::par.on(stream),
                    d_keys,
                    d_keys + total,
                    [=] BLST_DEVICE(u16 k)
            {
                return k <= static_cast<u16>(W);
            });
            if (!keys_ok)
                return cudaErrorInvalidValue;
        }
#endif

        auto idx_first = thrust::make_counting_iterator<size_t>(0);
        auto global_key_iter =
            thrust::make_transform_iterator(
                idx_first,
                [=] BLST_DEVICE(size_t idx) -> uint32_t
        {
            u16 k = d_keys[idx];
            size_t m = idx / static_cast<size_t>(N);
            size_t wspan = static_cast<size_t>(W + 1);
            return static_cast<uint32_t>(m * wspan + static_cast<size_t>(k));
        });

        size_t temp_storage_bytes = 0;
        {
            uint32_t* dummy_keys = nullptr;
            ValT* dummy_vals = nullptr;
            size_t* dummy_num_runs = nullptr;

            cudaError_t err =
                cub::DeviceReduce::ReduceByKey(
                    nullptr,
                    temp_storage_bytes,
                    global_key_iter,
                    dummy_keys,
                    src_iter,
                    dummy_vals,
                    dummy_num_runs,
                    XYZZAdd<Field>{},
                    static_cast<int>(total),
                    stream);
            if (err != cudaSuccess)
                return err;

            size_t required = 0;
            required += temp_storage_bytes;
            required = align_up_size(required, alignof(ValT));
            required += sizeof(ValT) * group_count;
            required = align_up_size(required, alignof(uint32_t));
            required += sizeof(uint32_t) * group_count;
            required = align_up_size(required, alignof(size_t));
            required += sizeof(size_t);

            if (d_workspace == nullptr)
            {
                workspace_bytes = required;
                return cudaSuccess;
            }

            if (workspace_bytes < required)
            {
                workspace_bytes = required;
                return cudaErrorInvalidValue;
            }
        }

        uint8_t* ws = static_cast<uint8_t*>(d_workspace);

        void* d_temp_storage = ws;
        ws += temp_storage_bytes;

        ws = align_up_ptr(ws, alignof(ValT));
        ValT* d_group_vals = reinterpret_cast<ValT*>(ws);
        ws += sizeof(ValT) * group_count;

        ws = align_up_ptr(ws, alignof(uint32_t));
        uint32_t* d_group_keys = reinterpret_cast<uint32_t*>(ws);
        ws += sizeof(uint32_t) * group_count;
        cudaMemsetAsync(d_group_keys, 0, sizeof(uint32_t) * group_count, stream);

        ws = align_up_ptr(ws, alignof(size_t));
        size_t* d_num_runs = reinterpret_cast<size_t*>(ws);
        ws += sizeof(size_t);

        cudaError_t err =
            cub::DeviceReduce::ReduceByKey(
                d_temp_storage,
                temp_storage_bytes,
                global_key_iter,
                d_group_keys,
                src_iter,
                d_group_vals,
                d_num_runs,
                XYZZAdd<Field>{},
                static_cast<int>(total),
                stream);
        if (err != cudaSuccess)
            return err;

        size_t out_elems = static_cast<size_t>(M) * static_cast<size_t>(W);
        cudaMemsetAsync(d_groups, 0, sizeof(Pointxyzz<Field>) * out_elems, stream);

        int blk = 256;
        size_t scatter_threads = group_count;
        int grd = static_cast<int>((scatter_threads + blk - 1) / blk);

        scatter_group_sums_multi<Field, ValT>
            << <grd, blk, 0, stream >> > (
                d_group_keys,
                d_group_vals,
                d_num_runs,
                W,
                M,
                d_groups);

        return cudaGetLastError();
    }


	// we have D groups, each with W values. for each group, we want to compute
	// the integrated sum of the values in reverse order.
	// 
	// High->low running sum over buckets(xyzz form) : 
	// computes ?? d??S_d from pre-summed buckets and writes 
	// the window???s contribution as Jacobian.
    template<class Field>
    BLST_GLOBAL void integrate_buckets_xyzz_kernel(
        const Pointxyzz<Field>* BLST_restrict d_groups_xyzz,  // [D*W]
        int D,
        int W,
		int wbits,
        Point<Field>* BLST_restrict d_window_sums       // [M*D] Jacobian out
    )
    {
        int tid = blockDim.x * blockIdx.x + threadIdx.x;
        if (tid >= D)
            return;

        const Pointxyzz<Field>* buckets = d_groups_xyzz + (size_t)tid * (size_t)W;

        Point<Field> out;
		point_integrate_buckets(out, buckets, wbits);
        d_window_sums[tid] = out;
    }

    template<class Field>
    cudaError_t integrate_buckets_xyzz(
        Pointxyzz<Field>* d_groups_xyzz,  // [M*D*W]
        int D,
        int W,
		int wbits,
        Point<Field>* d_window_sums,      // [M*D]
        cudaStream_t stream = 0)
    {
        int blk = 256;
        int grd = (D + blk - 1) / blk;
        integrate_buckets_xyzz_kernel<Field> << <grd, blk, 0, stream >> > (
            d_groups_xyzz,
            D,
            W,
			wbits,
            d_window_sums);
        return cudaGetLastError();
    }

    template<class Field>
    BLST_GLOBAL void fold_windows_kernel(
        const Point<Field>* BLST_restrict d_window_sums, // [M*D], index = m*D + d
        int M,
        int D,
        int wbits,
        Point<Field>* BLST_restrict d_outputs            // [M]
    )
    {
        int m = blockDim.x * blockIdx.x + threadIdx.x;
        if (m >= M)
            return;

        Point<Field> acc;
        vec_zero(&acc, sizeof(acc));

        for (int d = D - 1; d >= 0; --d)
        {
            for (int j = 0; j < wbits; ++j)
                Point<Field>::doubling(acc, acc);

            const Point<Field>& win = d_window_sums[m * D + d];
            Point<Field> tmp;
            Point<Field>::dadd(tmp, acc, win, NULL);
            acc = tmp;
        }

        d_outputs[m] = acc;
    }

    template<class Field>
    cudaError_t fold_windows(
        const Point<Field>* d_window_sums, // [M*D]
        int M,
        int D,
        int wbits,
        Point<Field>* d_outputs,           // [M]
        cudaStream_t stream = 0)
    {
        int blk = 256;
        int grd = (M + blk - 1) / blk;
        fold_windows_kernel<Field> << <grd, blk, 0, stream >> > (
            d_window_sums,
            M,
            D,
            wbits,
            d_outputs);
        return cudaGetLastError();
    }

    template<class Field>
    cudaError_t msm_booth_batched_pipeline(
        void* d_workspace,                          // nullptr => size query
        size_t& workspace_bytes,                    // in/out
        const AffinePoint<Field>* d_points,         // [N]
        const uint8_t* d_scalars,                   // [M * N * nbytes], flat (instance-major or scalar-major is fine as long as your digit-extract matches)
        int nbytes,
        int nbits,
        int N,                                      // points per MSM
        int M,                                      // number of MSMs
        int wbits,                                  // window bits (Booth)
        /* final outputs: */
        Point<Field>* d_outputs,                    // [M]
        cudaStream_t stream = 0)
    {
        if (N <= 0 || M <= 0 || wbits <= 0)
            return cudaErrorInvalidValue;

        int D = (nbits + wbits - 1) / wbits;        // number of windows
        int W = 1 << (wbits - 1);                   // Booth bucket count (no zero)

        // ---- Workspace sizing for downstream steps ----
        size_t needed = 0;
        size_t sort_bytes = 0;

        // Per-digit segmented sort workspace:
        {
            size_t sort_ws = 0;
            // Assume you invoke sort per MSM or per packed layout; here we size for one (D,N) pack.
            cudaError_t ec = sort_digits_segmented(
                nullptr,
                sort_ws,
                (u8*)nullptr,
                (u8*)nullptr,
                (IdxSign*) nullptr,
                (IdxSign*) nullptr,
                N,
                D,
                wbits,
                stream);

            if (ec != cudaSuccess)
                return ec;
            needed += sort_ws;                       // upper bound; reuse per digit / per MSM in-place
        }

        // Per-digit group reduce workspace:
        size_t reduce_ws = 0;
        {
            // Size for one call of reduce_groups
            Point<Field>* dummy_groups = nullptr;
            cudaError_t ec = reduce_groups<Field>(
                nullptr,
                reduce_ws,
                (u16*)nullptr,
                (Pointxyzz<Field>*)nullptr,
                (Point<Field>*)nullptr,
                D * M,
                N,
                W,
                stream);
            if (ec != cudaSuccess && ec != cudaErrorInvalidValue)
                return ec;
            needed += reduce_ws;                     // also reused per (m,d)
        }

        // Storage for all (m,d) window sums:
        size_t window_sums_bytes = sizeof(Point<Field>) * (size_t)M * (size_t)D;
        needed += window_sums_bytes;

        if (d_workspace == nullptr)
        {
            workspace_bytes = needed;
            return cudaSuccess;
        }
        if (workspace_bytes < needed)
        {
            workspace_bytes = needed;
            return cudaErrorInvalidValue;
        }

        // ---- Carve workspace ----
        uint8_t* base = static_cast<uint8_t*>(d_workspace);
        uint8_t* ws = base;

        // temp sort buffer (reuse)
        void* sort_temp = ws;
        // (actual size is returned by the query above, but we re-query inside call and CUB will honor larger buffers)
        ws += workspace_bytes;                       // simple carve-forward; you can make this precise if you prefer

        // temp reduce buffer (reuse)
        void* reduce_temp = sort_temp;               // reuse same large region

        // all window contributions (Jacobian)
        Point<Field>* d_window_sums = reinterpret_cast<Point<Field>*>(base + (needed - window_sums_bytes));

        // ---- Pipeline ----
        // 

        d3Span<u8> d_digit(reinterpret_cast<u8*>(ws), D, M, N);    // [D][M][N]
        ws += sizeof(u8) * d_digit.size();
        d3Span<IdxSign> d_index_signs(reinterpret_cast<IdxSign*>(ws), D, M, N);    // [D][M][N]
        ws += sizeof(IdxSign) * d_index_signs.size();
        
        d3Span<u8> d_sorted_digit(reinterpret_cast<u8*>(ws), D, M, N);    // [D][M][N]
        ws += sizeof(u8) * d_sorted_digit.size();
        d3Span<IdxSign> d_sorted_index_signs(reinterpret_cast<IdxSign*>(ws), D, M, N);    // [D][M][N]
        ws += sizeof(IdxSign) * d_sorted_index_signs.size();

        // we have M*N scalars, each nbytes long.
        // for each we need to extract D digits of wbits each.
        // The outer most dimention is the digit index (0..D-1).
        // the middle dimension is the MSM index (0..M-1).
        // the inner most dimension is the scalar index (0..N-1).
        auto ec = extract_digits_flat(
            d_scalars,       // [M*N*nbytes]
            (M * N),         // total scalars
            nbytes,
            nbits,
            wbits,
            D,
            d_digit.data_handle(),    // [D*M*N]  digit
            d_index_signs.data_handle()     // [D*M*N]  index|sign
        );
        if (ec != cudaSuccess)
            return ec;

        // collapse the D*M dimension so that we can sort all digits at once.
        auto as_2d = [](auto& span3d) {
            return d2Span<typename std::remove_reference<decltype(span3d)>::type::value_type>(
                span3d.data_handle(),
                span3d.extent(0) * span3d.extent(1),// D * M
				span3d.extent(2)); // N
			};


        // for each of the D*M digit slices of N elements,
        // sort the digits and store the result in d_index_2d.
        // This will contain the gather
        ec = sort_digits_segmented(
            sort_temp,
            sort_bytes,		 // in/out; CUB will honor larger buffers
            d_digit.data_handle(),	        // [D*N] per MSM or [D*M*N] packed
            d_sorted_digit.data_handle(),	// [D*N] per MSM or [D*M*N] packed
            d_index_signs.data_handle(),   // D*M vectors, each of length N. contains src idx
            d_sorted_index_signs.data_handle(),   // D*M vectors, each of length N. contains src idx
            N, 
            D, 
            wbits,
            stream);
        if (ec != cudaSuccess)
            return ec;

        // gather the points into Pointxyzz, applying sign.
        auto gather_iter = thrust::make_transform_iterator(
            d_sorted_index_signs.data_handle(),
            PackedToPointXYZZ<Field>{ d_points, M * N }
        );

        d3Span<Pointxyzz<Field>> d_groups(reinterpret_cast<Pointxyzz<Field>*>(ws), D, M, W); // [D][M][W]
		ws += sizeof(Pointxyzz<Field>) * d_groups.size();


        // now reduce each of the D*M digit slices of N points to W buckets.
        ec = reduce_groups<Field>(
            reduce_temp,
            reduce_ws,       // in/out; CUB will honor larger buffers
            d_sorted_digit.data_handle(),      // [D*N] per MSM or [D*M*N] packed
            gather_iter,
            d_groups.data_handle(), // [M*D*W] (carve from workspace above)
            D*M,
            N,
            W,
            stream);
        if (ec != cudaSuccess)
            return ec;

        ec = integrate_buckets_xyzz<Field>(
            d_groups,
            D*M,
            W,
			wbits,
            d_window_sums,
            stream);
        if (ec != cudaSuccess)
            return ec;

        // 4) fold D window sums per MSM with wbits doublings in between
        ec = fold_windows<Field>(
            d_window_sums,
            M,
            D,
            wbits,
            d_outputs,
            stream);
        if (ec != cudaSuccess)
            return ec;

        //// 
        //// For each MSM instance m and each digit d:
        ////   1) run segmented sort for that digit's N pairs
        ////   2) reduce groups -> W buckets
        //// After you???ve produced all (m,d) strips of W buckets, call the integration kernel and the fold.

        //// NOTE: The following is pseudo-looping ??? adapt to your exact memory layout of d_digit_keys/vals.
        //for (int m = 0; m < M; ++m)
        //{
        //	for (int d = 0; d < D; ++d)
        //	{
        //		// Pointers to the (N) keys/vals for this (m,d).
        //		u16* keys_md = d_digit_keys + (size_t)((m * D + d) * N);
        //		uint32_t* vals_md = d_digit_vals + (size_t)((m * D + d) * N);

        //		// 1) segmented sort over this single segment is effectively a normal sort;
        //		// if you already produced grouped-by-key order, you can skip this.
        //		{
        //			size_t ws_bytes_dummy = workspace_bytes;
        //			auto ec = sort_digits_segmented(
        //				sort_temp,
        //				ws_bytes_dummy,
        //				keys_md,
        //				vals_md,
        //				N,
        //				1,         // one segment if operating per (m,d)
        //				wbits,
        //				stream);
        //			if (ec != cudaSuccess)
        //				return ec;
        //		}

        //		// 2) reduce to W buckets
        //		Pointxyzz<Field>* buckets_md = /* carve/point to */ nullptr;
        //		// You likely already allocated a [M*D*W] output buffer for buckets.
        //		// If not, you can reuse a staging buffer and immediately integrate for this (m,d).

        //		size_t ws_reduce = reduce_ws;
        //		auto ec2 = reduce_groups<Field>(
        //			reduce_temp,
        //			ws_reduce,
        //			keys_md,
        //			vals_md,
        //			d_points,
        //			buckets_md,
        //			N,
        //			W,
        //			stream);
        //		if (ec2 != cudaSuccess)
        //			return ec2;
        //	}
        //}

        //// 3) integrate all (m,d) bucket strips ??? window_sums[m*D + d]
        ////    Use either the Jacobian or XYZZ variant depending on what you wrote out.
        //{
        //	auto ec = integrate_buckets_xyzz<Field>(
        //		d_groups,
        //		M,
        //		D,
        //		W,
        //		d_window_sums,
        //		stream);
        //	if (ec != cudaSuccess)
        //		return ec;
        //}

        //// 4) fold D window sums per MSM with wbits doublings in between
        //{
        //	auto ec = fold_windows<Field>(
        //		d_window_sums,
        //		M,
        //		D,
        //		wbits,
        //		d_outputs,
        //		stream);
        //	if (ec != cudaSuccess)
        //		return ec;
        //}

        return cudaSuccess;
    }


}
