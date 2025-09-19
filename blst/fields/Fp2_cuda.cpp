/*
 * Copyright Category
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/Fp2.h"
#include "blst/fields/Fp2_cuda.h"
#include "blst/exports.h"
#include "blst/common/errors.h"
#include "blst/common/cuda_launch.h"
 // printf for debugging
#include <stdio.h>

// Vectorized addition of Fp2 elements (CPU)
batched_error_code add_fp2_vec_cc(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

	while (num--) {
		add_fp2(*out, *a, *b);
		out++; a++; b++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for addition of Fp2 elements (GPU)
GLOBAL_FUNC void add_fp2_vec_ker(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		add_fp2(out[idx], a[idx], b[idx]);
	}
}

// Vectorized addition of Fp2 elements (GPU)
batched_error_code add_fp2_vec_cu(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

	auto fut = blst::launch_auto(add_fp2_vec_ker, num, out, a, b);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized addition of Fp2 elements
batched_error_code add_fp2_vec(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
#ifdef BLST_CUDA
	return add_fp2_vec_cu(out, a, b, num);
#else
	return add_fp2_vec_cc(out, a, b, num);
#endif
}

// Vectorized subtraction of Fp2 elements (CPU)
batched_error_code sub_fp2_vec_cc(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

	while (num--) {
		sub_fp2(*out, *a, *b);
		out++; a++; b++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for subtraction of Fp2 elements (GPU)
GLOBAL_FUNC void sub_fp2_vec_ker(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		sub_fp2(out[idx], a[idx], b[idx]);
	}
}

// Vectorized subtraction of Fp2 elements (GPU)
batched_error_code sub_fp2_vec_cu(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

	auto fut = blst::launch_auto(sub_fp2_vec_ker, num, out, a, b);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized subtraction of Fp2 elements
batched_error_code sub_fp2_vec(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
#ifdef BLST_CUDA
	return sub_fp2_vec_cu(out, a, b, num);
#else
	return sub_fp2_vec_cc(out, a, b, num);
#endif
}

// Vectorized multiplication of Fp2 elements (CPU)
batched_error_code mul_fp2_vec_cc(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

	while (num--) {
		mul_fp2(*out, *a, *b);
		out++; a++; b++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for multiplication of Fp2 elements (GPU)
GLOBAL_FUNC void mul_fp2_vec_ker(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		mul_fp2(out[idx], a[idx], b[idx]);
	}
}

// Vectorized multiplication of Fp2 elements (GPU)
batched_error_code mul_fp2_vec_cu(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

	auto fut = blst::launch_auto(mul_fp2_vec_ker, num, out, a, b);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized multiplication of Fp2 elements
batched_error_code mul_fp2_vec(vec384x* out, const vec384x* a, const vec384x* b, size_t num)
{
#ifdef BLST_CUDA
	return mul_fp2_vec_cu(out, a, b, num);
#else
	return mul_fp2_vec_cc(out, a, b, num);
#endif
}

// Vectorized squaring of Fp2 elements (CPU)
batched_error_code sqr_fp2_vec_cc(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	while (num--) {
		sqr_fp2(*out, *a);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for squaring of Fp2 elements (GPU)
GLOBAL_FUNC void sqr_fp2_vec_ker(vec384x* out, const vec384x* a, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		sqr_fp2(out[idx], a[idx]);
	}
}

// Vectorized squaring of Fp2 elements (GPU)
batched_error_code sqr_fp2_vec_cu(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	auto fut = blst::launch_auto(sqr_fp2_vec_ker, num, out, a);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized squaring of Fp2 elements
batched_error_code sqr_fp2_vec(vec384x* out, const vec384x* a, size_t num)
{
#ifdef BLST_CUDA
	return sqr_fp2_vec_cu(out, a, num);
#else
	return sqr_fp2_vec_cc(out, a, num);
#endif
}

// Vectorized conditional negation of Fp2 elements (CPU)
batched_error_code cneg_fp2_vec_cc(vec384x* out, const vec384x* a, int flag, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	while (num--) {
		cneg_fp2(*out, *a, flag);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for conditional negation of Fp2 elements (GPU)
GLOBAL_FUNC void cneg_fp2_vec_ker(vec384x* out, const vec384x* a, int flag, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		cneg_fp2(out[idx], a[idx], flag);
	}
}

// Vectorized conditional negation of Fp2 elements (GPU)
batched_error_code cneg_fp2_vec_cu(vec384x* out, const vec384x* a, int flag, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	auto fut = blst::launch_auto(cneg_fp2_vec_ker, num, out, a, flag);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized conditional negation of Fp2 elements
batched_error_code cneg_fp2_vec(vec384x* out, const vec384x* a, int flag, size_t num)
{
#ifdef BLST_CUDA
	return cneg_fp2_vec_cu(out, a, flag, num);
#else
	return cneg_fp2_vec_cc(out, a, flag, num);
#endif
}

// Vectorized multiplication of Fp2 elements by 3 (CPU)
batched_error_code mul3_fp2_vec_cc(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	while (num--) {
		mul_by_3_fp2(*out, *a);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for multiplication of Fp2 elements by 3 (GPU)
GLOBAL_FUNC void mul3_fp2_vec_ker(vec384x* out, const vec384x* a, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		mul_by_3_fp2(out[idx], a[idx]);
	}
}

// Vectorized multiplication of Fp2 elements by 3 (GPU)
batched_error_code mul3_fp2_vec_cu(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	auto fut = blst::launch_auto(mul3_fp2_vec_ker, num, out, a);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized multiplication of Fp2 elements by 3
batched_error_code mul3_fp2_vec(vec384x* out, const vec384x* a, size_t num)
{
#ifdef BLST_CUDA
	return mul3_fp2_vec_cu(out, a, num);
#else
	return mul3_fp2_vec_cc(out, a, num);
#endif
}

// Vectorized multiplication of Fp2 elements by 8 (CPU)
batched_error_code mul8_fp2_vec_cc(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	while (num--) {
		mul_by_8_fp2(*out, *a);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for multiplication of Fp2 elements by 8 (GPU)
GLOBAL_FUNC void mul8_fp2_vec_ker(vec384x* out, const vec384x* a, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		mul_by_8_fp2(out[idx], a[idx]);
	}
}

// Vectorized multiplication of Fp2 elements by 8 (GPU)
batched_error_code mul8_fp2_vec_cu(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	auto fut = blst::launch_auto(mul8_fp2_vec_ker, num, out, a);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized multiplication of Fp2 elements by 8
batched_error_code mul8_fp2_vec(vec384x* out, const vec384x* a, size_t num)
{
#ifdef BLST_CUDA
	return mul8_fp2_vec_cu(out, a, num);
#else
	return mul8_fp2_vec_cc(out, a, num);
#endif
}

// Vectorized left shift of Fp2 elements (CPU)
batched_error_code lshift_fp2_vec_cc(vec384x* out, const vec384x* a, size_t count, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	while (num--) {
		lshift_fp2(*out, *a, count);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for left shift of Fp2 elements (GPU)
GLOBAL_FUNC void lshift_fp2_vec_ker(vec384x* out, const vec384x* a, size_t count, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		lshift_fp2(out[idx], a[idx], count);
	}
}

// Vectorized left shift of Fp2 elements (GPU)
batched_error_code lshift_fp2_vec_cu(vec384x* out, const vec384x* a, size_t count, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	auto fut = blst::launch_auto(lshift_fp2_vec_ker, num, out, a, count);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized left shift of Fp2 elements
batched_error_code lshift_fp2_vec(vec384x* out, const vec384x* a, size_t count, size_t num)
{
#ifdef BLST_CUDA
	return lshift_fp2_vec_cu(out, a, count, num);
#else
	return lshift_fp2_vec_cc(out, a, count, num);
#endif
}

// Vectorized conjugation of Fp2 elements (CPU)
batched_error_code conjugate_fp2_vec_cc(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	while (num--) {
		conjugate_fp2(*out, *a);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for conjugation of Fp2 elements (GPU)
GLOBAL_FUNC void conjugate_fp2_vec_ker(vec384x* out, const vec384x* a, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		conjugate_fp2(out[idx], a[idx]);
	}
}

// Vectorized conjugation of Fp2 elements (GPU)
batched_error_code conjugate_fp2_vec_cu(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	auto fut = blst::launch_auto(conjugate_fp2_vec_ker, num, out, a);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized conjugation of Fp2 elements
batched_error_code conjugate_fp2_vec(vec384x* out, const vec384x* a, size_t num)
{
#ifdef BLST_CUDA
	return conjugate_fp2_vec_cu(out, a, num);
#else
	return conjugate_fp2_vec_cc(out, a, num);
#endif
}

// Vectorized inversion of Fp2 elements (CPU)
batched_error_code inverse_fp2_vec_cc(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	while (num--) {
		reciprocal_fp2(*out, *a);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for inversion of Fp2 elements (GPU)
GLOBAL_FUNC void inverse_fp2_vec_ker(vec384x* out, const vec384x* a, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		reciprocal_fp2(out[idx], a[idx]);
	}
}

// Vectorized inversion of Fp2 elements (GPU)
batched_error_code inverse_fp2_vec_cu(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	return cuda_to_batched(blst::launch_auto(inverse_fp2_vec_ker, num, out, a).get());
}
#endif

// Vectorized inversion of Fp2 elements
batched_error_code inverse_fp2_vec(vec384x* out, const vec384x* a, size_t num)
{
#ifdef BLST_CUDA
	return inverse_fp2_vec_cu(out, a, num);
#else
	return inverse_fp2_vec_cc(out, a, num);
#endif
}

// Vectorized square root of Fp2 elements (CPU)
batched_error_code sqrt_fp2_vec_cc(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

	while (num--) {
		sqrt_fp2(*out, *a);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for square root of Fp2 elements (GPU)
GLOBAL_FUNC void sqrt_fp2_vec_ker(vec384x* out, const vec384x* a, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		sqrt_fp2(out[idx], a[idx]);
	}
}

// Vectorized square root of Fp2 elements (GPU)
batched_error_code sqrt_fp2_vec_cu(vec384x* out, const vec384x* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	return cuda_to_batched(blst::launch_auto(sqrt_fp2_vec_ker, num, out, a).get());
}
#endif

// Vectorized square root of Fp2 elements
batched_error_code sqrt_fp2_vec(vec384x* out, const vec384x* a, size_t num)
{
#ifdef BLST_CUDA
	return sqrt_fp2_vec_cu(out, a, num);
#else
	return sqrt_fp2_vec_cc(out, a, num);
#endif
}

