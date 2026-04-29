/*
 * Copyright BTX authors
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/Fp.h"
#include "blst/fields/Fp_cuda.h"
#include "blst/exports.h"
#include "blst/common/cuda_launch.h"
#include "blst/backend/cuda/mul_mont_cuda.h"
#include "blst/backend/cuda/add_mod_cuda.h"

#ifdef BLST_CUDA
#define BLST_FP_KERNELS
#endif

	namespace blst
{

	namespace
	{
		void choose_launch(size_t n, dim3& grid, dim3& block)
		{
			int sm_count = 0;
			cudaDeviceProp prop{};
			auto ec = cudaGetDeviceProperties(&prop, 0);
			if (ec != cudaSuccess)
			{
				// default values
				const int blk = 256;
				size_t blocks = (n + blk - 1) / blk;
				blocks = std::max(blocks, (size_t)16);
				grid = dim3((unsigned)blocks);
				block = dim3(blk);
				return;
			}
			else
			{
				sm_count = prop.multiProcessorCount;
				const int blk = 256;
				size_t blocks = (n + blk - 1) / blk;
				blocks = std::max(blocks, (size_t)sm_count * 4);
				grid = dim3((unsigned)blocks);
				block = dim3(blk);
			}
		}
	}

	// Vectorized addition of Fp elements (CPU)
	batched_error_code add_fp_vec_cc(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			Fp::add(*out, *a, *b);
			out++; a++; b++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_FP_KERNELS
	// kernel for addition of Fp elements (GPU)
	BLST_GLOBAL void add_fp_vec_ker(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		size_t stride = gridDim.x * blockDim.x;
		for (size_t i = idx; i < num; i += stride)
		{
			auto pa = a + i;
			auto pb = b + i;
			auto pr = out + i;
			Fp::add(*pr, *pa, *pb);
		}
	}

	// Vectorized addition of Fp elements (GPU)
	// a1 a2 ... a_num
	// b1 b2 ... b_num
	//+_______________
	// out1 ...  out_num
	batched_error_code add_fp_vec_cu(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		dim3 grid, block;
		choose_launch(num, grid, block);
		cudaStream_t stream = 0; // default stream
		add_fp_vec_ker BLST_LAUNCH(grid, block, 0, stream)(out, a, b, num);

		// sync
		auto ec = cudaStreamSynchronize(stream);
		if (ec != cudaSuccess)
			return cuda_to_batched(ec);

		ec = cudaGetLastError();
		return cuda_to_batched(ec);
	}
#endif

	// Vectorized addition of Fp elements
	batched_error_code add_fp_vec(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
#ifdef BLST_FP_KERNELS
		return add_fp_vec_cu(out, a, b, num);
#else
		return add_fp_vec_cc(out, a, b, num);
#endif
	}


	batched_error_code sub_fp_vec_cc(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;


		while (num--) {
			Fp::sub(*out, *a, *b);
			out++; a++; b++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_FP_KERNELS
	// kernel for subtraction of Fp elements (GPU)
	BLST_GLOBAL void sub_fp_vec_ker(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		size_t stride = gridDim.x * blockDim.x;
		for (size_t i = idx; i < num; i += stride)
		{
			auto pa = a + i;
			auto pb = b + i;
			auto pr = out + i;
			Fp::sub(*pr, *pa, *pb);
		}
	}

	// Vectorized subtraction of Fp elements (GPU)
	batched_error_code sub_fp_vec_cu(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		dim3 grid, block;
		choose_launch(num, grid, block);
		cudaStream_t stream = 0; // default stream
		sub_fp_vec_ker BLST_LAUNCH(grid, block, 0, stream)(out, a, b, num);

		// sync
		auto ec = cudaStreamSynchronize(stream);
		if (ec != cudaSuccess)
			return cuda_to_batched(ec);

		ec = cudaGetLastError();
		return cuda_to_batched(ec);
	}
#endif

	// Vectorized subtraction of Fp elements
	batched_error_code sub_fp_vec(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
#ifdef BLST_FP_KERNELS
		return sub_fp_vec_cu(out, a, b, num);
#else
		return sub_fp_vec_cc(out, a, b, num);
#endif
	}

	// Vectorized multiplication of Fp elements (CPU)
	batched_error_code mul_fp_vec_cc(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;


		while (num--) {
			Fp::mul(*out, *a, *b);
			out++; a++; b++;
		}
		return BATCHED_SUCCESS;
	}
#ifdef BLST_FP_KERNELS
	// kernel for multiplication of Fp elements (GPU)
	BLST_GLOBAL void mul_fp_vec_ker(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		size_t stride = gridDim.x * blockDim.x;
		for (size_t i = idx; i < num; i += stride)
		{
			auto pa = a + i;
			auto pb = b + i;
			auto pr = out + i;
			Fp::mul(*pr, *pa, *pb);
		}
	}
	// Vectorized multiplication of Fp elements (GPU)
	batched_error_code mul_fp_vec_cu(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		dim3 grid, block;
		choose_launch(num, grid, block);
		cudaStream_t stream = 0; // default stream
		mul_fp_vec_ker BLST_LAUNCH(grid, block, 0, stream)(out, a, b, num);

		// sync
		auto ec = cudaStreamSynchronize(stream);
		if (ec != cudaSuccess)
			return cuda_to_batched(ec);

		ec = cudaGetLastError();
		return cuda_to_batched(ec);


		//auto fut = launch_auto(mul_fp_vec_ker, num, out, a, b);
		//return cuda_to_batched(fut.get());
	}
#endif
	// Vectorized multiplication of Fp elements
	batched_error_code mul_fp_vec(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
#ifdef BLST_FP_KERNELS
		return mul_fp_vec_cu(out, a, b, num);
#else
		return mul_fp_vec_cc(out, a, b, num);
#endif
	}

	// Vectorized squaring of Fp elements (CPU)
	batched_error_code sqr_fp_vec_cc(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;


		while (num--) {
			Fp::sqr(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}
#ifdef BLST_FP_KERNELS
	// kernel for squaring of Fp elements (GPU)
	BLST_GLOBAL void sqr_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		size_t stride = gridDim.x * blockDim.x;
		for (size_t i = idx; i < num; i += stride)
		{
			auto pa = a + i;
			auto pr = out + i;
			Fp::sqr(*pr, *pa);
		}
	}
	// Vectorized squaring of Fp elements (GPU)
	batched_error_code sqr_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		dim3 grid, block;
		choose_launch(num, grid, block);
		cudaStream_t stream = 0; // default stream
		sqr_fp_vec_ker BLST_LAUNCH(grid, block, 0, stream)(out, a, num);

		// sync
		auto ec = cudaStreamSynchronize(stream);
		if (ec != cudaSuccess)
			return cuda_to_batched(ec);

		ec = cudaGetLastError();
		return cuda_to_batched(ec);
	}
#endif
	// Vectorized squaring of Fp elements
	batched_error_code sqr_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_FP_KERNELS
		return sqr_fp_vec_cu(out, a, num);
#else
		return sqr_fp_vec_cc(out, a, num);
#endif
	}

	// Vectorized negation of Fp elements (CPU)
	batched_error_code cneg_fp_vec_cc(Fp* out, const Fp* a, int flag, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;


		while (num--) {
			Fp::cneg(*out, *a, flag);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}
#ifdef BLST_FP_KERNELS
	// kernel for negation of Fp elements (GPU)
	BLST_GLOBAL void cneg_fp_vec_ker(Fp* out, const Fp* a, int flag, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		size_t stride = gridDim.x * blockDim.x;
		for (size_t i = idx; i < num; i += stride)
		{
			auto pa = a + i;
			auto pr = out + i;
			Fp::cneg(*pr, *pa, flag);
		}
	}
	// Vectorized negation of Fp elements (GPU)
	batched_error_code cneg_fp_vec_cu(Fp* out, const Fp* a, int flag, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		dim3 grid, block;
		choose_launch(num, grid, block);
		cudaStream_t stream = 0; // default stream
		cneg_fp_vec_ker BLST_LAUNCH(grid, block, 0, stream)(out, a, flag, num);

		// sync
		auto ec = cudaStreamSynchronize(stream);
		if (ec != cudaSuccess)
			return cuda_to_batched(ec);

		ec = cudaGetLastError();
		return cuda_to_batched(ec);
	}
#endif

	// Vectorized negation of Fp elements
	batched_error_code cneg_fp_vec(Fp* out, const Fp* a, int flag, size_t num)
	{
#ifdef BLST_FP_KERNELS
		return cneg_fp_vec_cu(out, a, flag, num);
#else
		return cneg_fp_vec_cc(out, a, flag, num);
#endif
	}



	// mult by 3

	// Vectorized multiplication of Fp elements by 3 (CPU)
	batched_error_code mul3_fp_vec_cc(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;


		while (num--) {
			Fp::mul_by_3(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}
#ifdef BLST_FP_KERNELS
	// kernel for multiplication of Fp elements by 3 (GPU)
	BLST_GLOBAL void mul3_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		size_t stride = gridDim.x * blockDim.x;
		for (size_t i = idx; i < num; i += stride)
		{
			auto pa = a + i;
			auto pr = out + i;
			Fp::mul_by_3(*pr, *pa);
		}
	}
	// Vectorized multiplication of Fp elements by 3 (GPU)
	batched_error_code mul3_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		dim3 grid, block;
		choose_launch(num, grid, block);
		cudaStream_t stream = 0; // default stream
		mul3_fp_vec_ker BLST_LAUNCH(grid, block, 0, stream)(out, a, num);

		// sync
		auto ec = cudaStreamSynchronize(stream);
		if (ec != cudaSuccess)
			return cuda_to_batched(ec);

		ec = cudaGetLastError();
		return cuda_to_batched(ec);
	}
#endif
	// Vectorized multiplication of Fp elements by 3
	batched_error_code mul3_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_FP_KERNELS
		return mul3_fp_vec_cu(out, a, num);
#else
		return mul3_fp_vec_cc(out, a, num);
#endif
	}

	// mul by 8 

	// Vectorized multiplication of Fp elements by 8 (CPU)
	batched_error_code mul8_fp_vec_cc(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;


		while (num--) {
			Fp::mul_by_8(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}
#ifdef BLST_FP_KERNELS
	// kernel for multiplication of Fp elements by 8 (GPU)
	BLST_GLOBAL void mul8_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		size_t stride = gridDim.x * blockDim.x;
		for (size_t i = idx; i < num; i += stride)
		{
			auto pa = a + i;
			auto pr = out + i;
			Fp::mul_by_8(*pr, *pa);
		}
	}
	// Vectorized multiplication of Fp elements by 8 (GPU)
	batched_error_code mul8_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		dim3 grid, block;
		choose_launch(num, grid, block);
		cudaStream_t stream = 0; // default stream
		mul8_fp_vec_ker BLST_LAUNCH(grid, block, 0, stream)(out, a, num);

		// sync
		auto ec = cudaStreamSynchronize(stream);
		if (ec != cudaSuccess)
			return cuda_to_batched(ec);

		ec = cudaGetLastError();
		return cuda_to_batched(ec);
	}
#endif
	// Vectorized multiplication of Fp elements by 8
	batched_error_code mul8_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_FP_KERNELS
		return mul8_fp_vec_cu(out, a, num);
#else
		return mul8_fp_vec_cc(out, a, num);
#endif
	}


	// div by 2

	// Vectorized division of Fp elements by 2 (CPU)
	batched_error_code div2_fp_vec_cc(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;


		while (num--) {
			Fp::div_by_2(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}
#ifdef BLST_FP_KERNELS
	// kernel for division of Fp elements by 2 (GPU)
	BLST_GLOBAL void div2_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		size_t stride = gridDim.x * blockDim.x;
		for (size_t i = idx; i < num; i += stride)
		{
			auto pa = a + i;
			auto pr = out + i;
			Fp::div_by_2(*pr, *pa);
		}
	}
	// Vectorized division of Fp elements by 2 (GPU)
	batched_error_code div2_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		dim3 grid, block;
		choose_launch(num, grid, block);
		cudaStream_t stream = 0; // default stream
		div2_fp_vec_ker BLST_LAUNCH(grid, block, 0, stream)(out, a, num);

		// sync
		auto ec = cudaStreamSynchronize(stream);
		if (ec != cudaSuccess)
			return cuda_to_batched(ec);

		ec = cudaGetLastError();
		return cuda_to_batched(ec);
	}
#endif
	// Vectorized division of Fp elements by 2
	batched_error_code div2_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_FP_KERNELS
		return div2_fp_vec_cu(out, a, num);
#else
		return div2_fp_vec_cc(out, a, num);
#endif
	}

	// lshift 

	// Vectorized left shift of Fp elements (CPU)
	batched_error_code lshift_fp_vec_cc(Fp* out, const Fp* a, size_t count, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;


		while (num--) {
			Fp::lshift(*out, *a, count);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_FP_KERNELS
	// kernel for left shift of Fp elements (GPU)
	BLST_GLOBAL void lshift_fp_vec_ker(Fp* out, const Fp* a, size_t count, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		size_t stride = gridDim.x * blockDim.x;
		for (size_t i = idx; i < num; i += stride)
		{
			auto pa = a + i;
			auto pr = out + i;
			Fp::lshift(*pr, *pa, count);
		}
	}
	// Vectorized left shift of Fp elements (GPU)
	batched_error_code lshift_fp_vec_cu(Fp* out, const Fp* a, size_t count, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		dim3 grid, block;
		choose_launch(num, grid, block);
		cudaStream_t stream = 0; // default stream
		lshift_fp_vec_ker BLST_LAUNCH(grid, block, 0, stream)(out, a, count, num);

		// sync
		auto ec = cudaStreamSynchronize(stream);
		if (ec != cudaSuccess)
			return cuda_to_batched(ec);

		ec = cudaGetLastError();
		return cuda_to_batched(ec);
	}
#endif
	// Vectorized left shift of Fp elements
	batched_error_code lshift_fp_vec(Fp* out, const Fp* a, size_t count, size_t num)
	{
#ifdef BLST_FP_KERNELS
		return lshift_fp_vec_cu(out, a, count, num);
#else
		return lshift_fp_vec_cc(out, a, count, num);
#endif
	}


	// to fp

	// Vectorized conversion from a regular representation to a montgomery Fp (CPU)
	batched_error_code to_fp_vec_cc(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;


		while (num--) {
			Fp::to_mont(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}
#ifdef BLST_FP_KERNELS
	// kernel for conversion from a regular representation to a montgomery Fp (GPU)
	BLST_GLOBAL void to_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		size_t stride = gridDim.x * blockDim.x;
		for (size_t i = idx; i < num; i += stride)
		{
			auto pa = a + i;
			auto pr = out + i;
			Fp::to_mont(*pr, *pa);
		}
	}
	// Vectorized conversion from a regular representation to a montgomery Fp (GPU)
	batched_error_code to_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		dim3 grid, block;
		choose_launch(num, grid, block);
		cudaStream_t stream = 0; // default stream
		to_fp_vec_ker BLST_LAUNCH(grid, block, 0, stream)(out, a, num);

		// sync
		auto ec = cudaStreamSynchronize(stream);
		if (ec != cudaSuccess)
			return cuda_to_batched(ec);

		ec = cudaGetLastError();
		return cuda_to_batched(ec);
	}
#endif
	// Vectorized conversion from a regular representation to a montgomery Fp
	batched_error_code to_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_FP_KERNELS
		return to_fp_vec_cu(out, a, num);
#else
		return to_fp_vec_cc(out, a, num);
#endif
	}

	// from fp

	// Vectorized conversion from a montgomery Fp to a regular representation (CPU)
	batched_error_code from_fp_vec_cc(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;


		while (num--) {
			Fp::from_mont(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}
#ifdef BLST_FP_KERNELS
	// kernel for conversion from a montgomery Fp to a regular representation (GPU)
	BLST_GLOBAL void from_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		size_t stride = gridDim.x * blockDim.x;
		for (size_t i = idx; i < num; i += stride)
		{
			auto pa = a + i;
			auto pr = out + i;
			Fp::from_mont(*pr, *pa);
		}
	}
	// Vectorized conversion from a montgomery Fp to a regular representation (GPU)
	batched_error_code from_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		dim3 grid, block;
		choose_launch(num, grid, block);
		cudaStream_t stream = 0; // default stream
		from_fp_vec_ker BLST_LAUNCH(grid, block, 0, stream)(out, a, num);

		// sync
		auto ec = cudaStreamSynchronize(stream);
		if (ec != cudaSuccess)
			return cuda_to_batched(ec);

		ec = cudaGetLastError();
		return cuda_to_batched(ec);
	}
#endif
	// Vectorized conversion from a montgomery Fp to a regular representation
	batched_error_code from_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_FP_KERNELS
		return from_fp_vec_cu(out, a, num);
#else
		return from_fp_vec_cc(out, a, num);
#endif
	}



}
