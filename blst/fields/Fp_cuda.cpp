/*
 * Copyright Category
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/Fp.h"
#include "blst/fields/Fp_cuda.h"
#include "blst/exports.h"
#include "blst/common/cuda_launch.h"

namespace blst
{

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

#ifdef BLST_CUDA
	// kernel for addition of Fp elements (GPU)
	GLOBAL_FUNC void add_fp_vec_ker(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp::add(out[idx], a[idx], b[idx]);
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

		auto fut = launch_auto(add_fp_vec_ker, num, out, a, b);
		return cuda_to_batched(fut.get());
	}
#endif

	// Vectorized addition of Fp elements
	batched_error_code add_fp_vec(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
#ifdef BLST_CUDA
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

#ifdef BLST_CUDA
	// kernel for subtraction of Fp elements (GPU)
	GLOBAL_FUNC void sub_fp_vec_ker(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp::sub(out[idx], a[idx], b[idx]);
		}
	}

	// Vectorized subtraction of Fp elements (GPU)
	batched_error_code sub_fp_vec_cu(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		auto fut = launch_auto(sub_fp_vec_ker, num, out, a, b);
		return cuda_to_batched(fut.get());
	}
#endif

	// Vectorized subtraction of Fp elements
	batched_error_code sub_fp_vec(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
#ifdef BLST_CUDA
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
#ifdef BLST_CUDA
	// kernel for multiplication of Fp elements (GPU)
	GLOBAL_FUNC void mul_fp_vec_ker(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp::mul(out[idx], a[idx], b[idx]);
		}
	}
	// Vectorized multiplication of Fp elements (GPU)
	batched_error_code mul_fp_vec_cu(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		auto fut = launch_auto(mul_fp_vec_ker, num, out, a, b);
		return cuda_to_batched(fut.get());
	}
#endif
	// Vectorized multiplication of Fp elements
	batched_error_code mul_fp_vec(Fp* out, const Fp* a, const Fp* b, size_t num)
	{
#ifdef BLST_CUDA
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
#ifdef BLST_CUDA
	// kernel for squaring of Fp elements (GPU)
	GLOBAL_FUNC void sqr_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp::sqr(out[idx], a[idx]);
		}
	}
	// Vectorized squaring of Fp elements (GPU)
	batched_error_code sqr_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto fut = launch_auto(sqr_fp_vec_ker, num, out, a);
		return cuda_to_batched(fut.get());
	}
#endif
	// Vectorized squaring of Fp elements
	batched_error_code sqr_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_CUDA
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
#ifdef BLST_CUDA
	// kernel for negation of Fp elements (GPU)
	GLOBAL_FUNC void cneg_fp_vec_ker(Fp* out, const Fp* a, int flag, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp::cneg(out[idx], a[idx], flag);
		}
	}
	// Vectorized negation of Fp elements (GPU)
	batched_error_code cneg_fp_vec_cu(Fp* out, const Fp* a, int flag, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto fut = launch_auto(cneg_fp_vec_ker, num, out, a, flag);
		return cuda_to_batched(fut.get());
	}
#endif

	// Vectorized negation of Fp elements
	batched_error_code cneg_fp_vec(Fp* out, const Fp* a, int flag, size_t num)
	{
#ifdef BLST_CUDA
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
#ifdef BLST_CUDA
	// kernel for multiplication of Fp elements by 3 (GPU)
	GLOBAL_FUNC void mul3_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp::mul_by_3(out[idx], a[idx]);
		}
	}
	// Vectorized multiplication of Fp elements by 3 (GPU)
	batched_error_code mul3_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto fut = launch_auto(mul3_fp_vec_ker, num, out, a);
		return cuda_to_batched(fut.get());
	}
#endif
	// Vectorized multiplication of Fp elements by 3
	batched_error_code mul3_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_CUDA
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
#ifdef BLST_CUDA
	// kernel for multiplication of Fp elements by 8 (GPU)
	GLOBAL_FUNC void mul8_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp::mul_by_8(out[idx], a[idx]);
		}
	}
	// Vectorized multiplication of Fp elements by 8 (GPU)
	batched_error_code mul8_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto fut = launch_auto(mul8_fp_vec_ker, num, out, a);
		return cuda_to_batched(fut.get());
	}
#endif
	// Vectorized multiplication of Fp elements by 8
	batched_error_code mul8_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_CUDA
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
#ifdef BLST_CUDA
	// kernel for division of Fp elements by 2 (GPU)
	GLOBAL_FUNC void div2_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp::div_by_2(out[idx], a[idx]);
		}
	}
	// Vectorized division of Fp elements by 2 (GPU)
	batched_error_code div2_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto fut = launch_auto(div2_fp_vec_ker, num, out, a);
		return cuda_to_batched(fut.get());
	}
#endif
	// Vectorized division of Fp elements by 2
	batched_error_code div2_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_CUDA
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

#ifdef BLST_CUDA
	// kernel for left shift of Fp elements (GPU)
	GLOBAL_FUNC void lshift_fp_vec_ker(Fp* out, const Fp* a, size_t count, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp::lshift(out[idx], a[idx], count);
		}
	}
	// Vectorized left shift of Fp elements (GPU)
	batched_error_code lshift_fp_vec_cu(Fp* out, const Fp* a, size_t count, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto fut = launch_auto(lshift_fp_vec_ker, num, out, a, count);
		return cuda_to_batched(fut.get());
	}
#endif
	// Vectorized left shift of Fp elements
	batched_error_code lshift_fp_vec(Fp* out, const Fp* a, size_t count, size_t num)
	{
#ifdef BLST_CUDA
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
#ifdef BLST_CUDA
	// kernel for conversion from a regular representation to a montgomery Fp (GPU)
	GLOBAL_FUNC void to_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp::to_mont(out[idx], a[idx]);
		}
	}
	// Vectorized conversion from a regular representation to a montgomery Fp (GPU)
	batched_error_code to_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto fut = launch_auto(to_fp_vec_ker, num, out, a);
		return cuda_to_batched(fut.get());
	}
#endif
	// Vectorized conversion from a regular representation to a montgomery Fp
	batched_error_code to_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_CUDA
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
#ifdef BLST_CUDA
	// kernel for conversion from a montgomery Fp to a regular representation (GPU)
	GLOBAL_FUNC void from_fp_vec_ker(Fp* out, const Fp* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp::from_mont(out[idx], a[idx]);
		}
	}
	// Vectorized conversion from a montgomery Fp to a regular representation (GPU)
	batched_error_code from_fp_vec_cu(Fp* out, const Fp* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto fut = launch_auto(from_fp_vec_ker, num, out, a);
		return cuda_to_batched(fut.get());
	}
#endif
	// Vectorized conversion from a montgomery Fp to a regular representation
	batched_error_code from_fp_vec(Fp* out, const Fp* a, size_t num)
	{
#ifdef BLST_CUDA
		return from_fp_vec_cu(out, a, num);
#else
		return from_fp_vec_cc(out, a, num);
#endif
	}



}

