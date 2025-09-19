/*
 * Copyright Category
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/Fr.h"
#include "blst/exports.h"
#include "blst/fields/Fr_cuda.h"
#include "blst/common/errors.h"
#include "blst/common/cuda_launch.h"

// ====================
// ====================
// ====================
// ====================
// ====================
// 
// Fr vectorized operations
// 
// ====================
// ====================
// ====================
// ====================
// ====================

// Vectorized addition of Fr elements (CPU)
batched_error_code add_fr_vec_cc(vec256* out, const vec256* a, const vec256* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_add(*out, *a, *b);
		out++; a++; b++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for addition of Fr elements (GPU)
GLOBAL_FUNC void add_fr_vec_ker(vec256* out, const vec256* a, const vec256* b, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_add(out[idx], a[idx], b[idx]);
	}
}

// Vectorized addition of Fr elements (GPU)
batched_error_code add_fr_vec_cu(vec256* out, const vec256* a, const vec256* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(add_fr_vec_ker, num, out, a, b);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized addition of Fr elements
batched_error_code add_fr_vec(vec256* out, const vec256* a, const vec256* b, size_t num)
{
#ifdef BLST_CUDA
	return add_fr_vec_cu(out, a, b, num);
#else
	return add_fr_vec_cc(out, a, b, num);
#endif
}

// Vectorized subtraction of Fr elements (CPU)
batched_error_code sub_fr_vec_cc(vec256* out, const vec256* a, const vec256* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_sub(*out, *a, *b);
		out++; a++; b++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for subtraction of Fr elements (GPU)
GLOBAL_FUNC void sub_fr_vec_ker(vec256* out, const vec256* a, const vec256* b, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_sub(out[idx], a[idx], b[idx]);
	}
}

// Vectorized subtraction of Fr elements (GPU)
batched_error_code sub_fr_vec_cu(vec256* out, const vec256* a, const vec256* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(sub_fr_vec_ker, num, out, a, b);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized subtraction of Fr elements
batched_error_code sub_fr_vec(vec256* out, const vec256* a, const vec256* b, size_t num)
{
#ifdef BLST_CUDA
	return sub_fr_vec_cu(out, a, b, num);
#else
	return sub_fr_vec_cc(out, a, b, num);
#endif
}

// Vectorized multiplication of Fr elements (CPU)
batched_error_code mul_fr_vec_cc(vec256* out, const vec256* a, const vec256* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_mul(*out, *a, *b);
		out++; a++; b++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for multiplication of Fr elements (GPU)
GLOBAL_FUNC void mul_fr_vec_ker(vec256* out, const vec256* a, const vec256* b, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_mul(out[idx], a[idx], b[idx]);
	}
}

// Vectorized multiplication of Fr elements (GPU)
batched_error_code mul_fr_vec_cu(vec256* out, const vec256* a, const vec256* b, size_t num)
{
	if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(mul_fr_vec_ker, num, out, a, b);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized multiplication of Fr elements
batched_error_code mul_fr_vec(vec256* out, const vec256* a, const vec256* b, size_t num)
{
#ifdef BLST_CUDA
	return mul_fr_vec_cu(out, a, b, num);
#else
	return mul_fr_vec_cc(out, a, b, num);
#endif
}

// Vectorized squaring of Fr elements (CPU)
batched_error_code sqr_fr_vec_cc(vec256* out, const vec256* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_sqr(*out, *a);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for squaring of Fr elements (GPU)
GLOBAL_FUNC void sqr_fr_vec_ker(vec256* out, const vec256* a, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_sqr(out[idx], a[idx]);
	}
}

// Vectorized squaring of Fr elements (GPU)
batched_error_code sqr_fr_vec_cu(vec256* out, const vec256* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(sqr_fr_vec_ker, num, out, a);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized squaring of Fr elements
batched_error_code sqr_fr_vec(vec256* out, const vec256* a, size_t num)
{
#ifdef BLST_CUDA
	return sqr_fr_vec_cu(out, a, num);
#else
	return sqr_fr_vec_cc(out, a, num);
#endif
}

// Vectorized conditional negation of Fr elements (CPU)
batched_error_code cneg_fr_vec_cc(vec256* out, const vec256* a, int flag, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_cneg(*out, *a, flag);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for conditional negation of Fr elements (GPU)
GLOBAL_FUNC void cneg_fr_vec_ker(vec256* out, const vec256* a, int flag, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_cneg(out[idx], a[idx], flag);
	}
}

// Vectorized conditional negation of Fr elements (GPU)
batched_error_code cneg_fr_vec_cu(vec256* out, const vec256* a, int flag, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(cneg_fr_vec_ker, num, out, a, flag);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized conditional negation of Fr elements
batched_error_code cneg_fr_vec(vec256* out, const vec256* a, int flag, size_t num)
{
#ifdef BLST_CUDA
	return cneg_fr_vec_cu(out, a, flag, num);
#else
	return cneg_fr_vec_cc(out, a, flag, num);
#endif
}

// Vectorized multiplication of Fr elements by 3 (CPU)
batched_error_code mul3_fr_vec_cc(vec256* out, const vec256* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_mul_by_3(*out, *a);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for multiplication of Fr elements by 3 (GPU)
GLOBAL_FUNC void mul3_fr_vec_ker(vec256* out, const vec256* a, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_mul_by_3(out[idx], a[idx]);
	}
}

// Vectorized multiplication of Fr elements by 3 (GPU)
batched_error_code mul3_fr_vec_cu(vec256* out, const vec256* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(mul3_fr_vec_ker, num, out, a);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized multiplication of Fr elements by 3
batched_error_code mul3_fr_vec(vec256* out, const vec256* a, size_t num)
{
#ifdef BLST_CUDA
	return mul3_fr_vec_cu(out, a, num);
#else
	return mul3_fr_vec_cc(out, a, num);
#endif
}

// Vectorized left shift of Fr elements (CPU)
batched_error_code lshift_fr_vec_cc(vec256* out, const vec256* a, size_t count, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_lshift(*out, *a, count);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for left shift of Fr elements (GPU)
GLOBAL_FUNC void lshift_fr_vec_ker(vec256* out, const vec256* a, size_t count, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_lshift(out[idx], a[idx], count);
	}
}

// Vectorized left shift of Fr elements (GPU)
batched_error_code lshift_fr_vec_cu(vec256* out, const vec256* a, size_t count, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(lshift_fr_vec_ker, num, out, a, count);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized left shift of Fr elements
batched_error_code lshift_fr_vec(vec256* out, const vec256* a, size_t count, size_t num)
{
#ifdef BLST_CUDA
	return lshift_fr_vec_cu(out, a, count, num);
#else
	return lshift_fr_vec_cc(out, a, count, num);
#endif
}

// Vectorized right shift of Fr elements (CPU)
batched_error_code rshift_fr_vec_cc(vec256* out, const vec256* a, size_t count, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_rshift(*out, *a, count);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for right shift of Fr elements (GPU)
GLOBAL_FUNC void rshift_fr_vec_ker(vec256* out, const vec256* a, size_t count, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_rshift(out[idx], a[idx], count);
	}
}

// Vectorized right shift of Fr elements (GPU)
batched_error_code rshift_fr_vec_cu(vec256* out, const vec256* a, size_t count, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(rshift_fr_vec_ker, num, out, a, count);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized right shift of Fr elements
batched_error_code rshift_fr_vec(vec256* out, const vec256* a, size_t count, size_t num)
{
#ifdef BLST_CUDA
	return rshift_fr_vec_cu(out, a, count, num);
#else
	return rshift_fr_vec_cc(out, a, count, num);
#endif
}

// Vectorized conversion from regular to Montgomery form for Fr elements (CPU)
batched_error_code to_fr_vec_cc(vec256* out, const vec256* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_to(*out, *a);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for conversion from regular to Montgomery form for Fr elements (GPU)
GLOBAL_FUNC void to_fr_vec_ker(vec256* out, const vec256* a, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_to(out[idx], a[idx]);
	}
}

// Vectorized conversion from regular to Montgomery form for Fr elements (GPU)
batched_error_code to_fr_vec_cu(vec256* out, const vec256* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(to_fr_vec_ker, num, out, a);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized conversion from regular to Montgomery form for Fr elements
batched_error_code to_fr_vec(vec256* out, const vec256* a, size_t num)
{
#ifdef BLST_CUDA
	return to_fr_vec_cu(out, a, num);
#else
	return to_fr_vec_cc(out, a, num);
#endif
}

// Vectorized conversion from Montgomery to regular form for Fr elements (CPU)
batched_error_code from_fr_vec_cc(vec256* out, const vec256* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_from(*out, *a);
		out++; a++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for conversion from Montgomery to regular form for Fr elements (GPU)
GLOBAL_FUNC void from_fr_vec_ker(vec256* out, const vec256* a, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_from(out[idx], a[idx]);
	}
}

// Vectorized conversion from Montgomery to regular form for Fr elements (GPU)
batched_error_code from_fr_vec_cu(vec256* out, const vec256* a, size_t num)
{
	if (!out || !a) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(from_fr_vec_ker, num, out, a);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized conversion from Montgomery to regular form for Fr elements
batched_error_code from_fr_vec(vec256* out, const vec256* a, size_t num)
{
#ifdef BLST_CUDA
	return from_fr_vec_cu(out, a, num);
#else
	return from_fr_vec_cc(out, a, num);
#endif
}

// Vectorized CT butterfly operations for Fr elements (CPU)
batched_error_code ct_bfly_fr_vec_cc(vec256* x0, vec256* x1, const vec256* twiddle, size_t num)
{
	if (!x0 || !x1 || !twiddle) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_ct_bfly(*x0, *x1, *twiddle);
		x0++; x1++; twiddle++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for CT butterfly operations for Fr elements (GPU)
GLOBAL_FUNC void ct_bfly_fr_vec_ker(vec256* x0, vec256* x1, const vec256* twiddle, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_ct_bfly(x0[idx], x1[idx], twiddle[idx]);
	}
}

// Vectorized CT butterfly operations for Fr elements (GPU)
batched_error_code ct_bfly_fr_vec_cu(vec256* x0, vec256* x1, const vec256* twiddle, size_t num)
{
	if (!x0 || !x1 || !twiddle) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(ct_bfly_fr_vec_ker, num, x0, x1, twiddle);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized CT butterfly operations for Fr elements
batched_error_code ct_bfly_fr_vec(vec256* x0, vec256* x1, const vec256* twiddle, size_t num)
{
#ifdef BLST_CUDA
	return ct_bfly_fr_vec_cu(x0, x1, twiddle, num);
#else
	return ct_bfly_fr_vec_cc(x0, x1, twiddle, num);
#endif
}

// Vectorized GS butterfly operations for Fr elements (CPU)
batched_error_code gs_bfly_fr_vec_cc(vec256* x0, vec256* x1, const vec256* twiddle, size_t num)
{
	if (!x0 || !x1 || !twiddle) return BATCHED_ERROR_NULL_POINTER;
	

	while (num--) {
		fr_gs_bfly(*x0, *x1, *twiddle);
		x0++; x1++; twiddle++;
	}
	return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// kernel for GS butterfly operations for Fr elements (GPU)
GLOBAL_FUNC void gs_bfly_fr_vec_ker(vec256* x0, vec256* x1, const vec256* twiddle, size_t num)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num) {
		fr_gs_bfly(x0[idx], x1[idx], twiddle[idx]);
	}
}

// Vectorized GS butterfly operations for Fr elements (GPU)
batched_error_code gs_bfly_fr_vec_cu(vec256* x0, vec256* x1, const vec256* twiddle, size_t num)
{
	if (!x0 || !x1 || !twiddle) return BATCHED_ERROR_NULL_POINTER;
	
	auto fut = blst::launch_auto(gs_bfly_fr_vec_ker, num, x0, x1, twiddle);
	return cuda_to_batched(fut.get());
}
#endif

// Vectorized GS butterfly operations for Fr elements
batched_error_code gs_bfly_fr_vec(vec256* x0, vec256* x1, 
	const vec256* twiddle, size_t num)
{
	#ifdef BLST_CUDA
	return gs_bfly_fr_vec_cu(x0, x1, twiddle, num);
#else
	return gs_bfly_fr_vec_cc(x0, x1, twiddle, num);
#endif
}


