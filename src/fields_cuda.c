///*
// * Copyright Supranational LLC
// * Licensed under the Apache License, Version 2.0, see LICENSE for details.
// * SPDX-License-Identifier: Apache-2.0
// */
//
//#include "fields.h"
//BLST_EXTERN_C_BEGIN
//
//
//// Vectorized addition of Fp elements (CPU)
//BLST_INLINE void add_fp_vec_cc(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//	while (num--) {
//		add_fp(*out, *a, *b);
//		out++; a++; b++;
//	}
//}
//
//#ifdef BLST_CUDA
//// kernel for addition of Fp elements (GPU)
//GLOBAL_FUNC BLST_INLINE void add_fp_vec_ker(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	if (idx < num) {
//		add_fp(out[idx], a[idx], b[idx]);
//	}
//}
//
//// Vectorized addition of Fp elements (GPU)
//BLST_INLINE void add_fp_vec_cu(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//	size_t threads_per_block = 256;
//	size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
//	add_fp_vec_ker << <num_blocks, threads_per_block >> > (out, a, b, num);
//	cudaDeviceSynchronize();
//}
//#endif
//
//// Vectorized addition of Fp elements
//void add_fp_vec(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//#ifdef BLST_CUDA
//	add_fp_vec_cu(out, a, b, num);
//#else
//	add_fp_vec_cc(out, a, b, num);
//#endif
//}
//
//
//BLST_INLINE void sub_fp_vec_cc(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//	while (num--) {
//		sub_fp(*out, *a, *b);
//		out++; a++; b++;
//	}
//}
//
//#ifdef BLST_CUDA
//// kernel for subtraction of Fp elements (GPU)
//GLOBAL_FUNC BLST_INLINE void sub_fp_vec_ker(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	if (idx < num) {
//		sub_fp(out[idx], a[idx], b[idx]);
//	}
//}
//
//// Vectorized subtraction of Fp elements (GPU)
//BLST_INLINE void sub_fp_vec_cu(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//	size_t threads_per_block = 256;
//	size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
//	sub_fp_vec_ker << <num_blocks, threads_per_block >> > (out, a, b, num);
//	cudaDeviceSynchronize();
//}
//#endif
//
//// Vectorized subtraction of Fp elements
//BLST_INLINE void sub_fp_vec(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//#ifdef BLST_CUDA
//	sub_fp_vec_cu(out, a, b, num);
//#else
//	sub_fp_vec_cc(out, a, b, num);
//#endif
//}
//
//// Vectorized multiplication of Fp elements (CPU)
//BLST_INLINE void mul_fp_vec_cc(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//	while (num--) {
//		mul_fp(*out, *a, *b);
//		out++; a++; b++;
//	}
//}
//#ifdef BLST_CUDA
//// kernel for multiplication of Fp elements (GPU)
//GLOBAL_FUNC BLST_INLINE void mul_fp_vec_ker(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	if (idx < num) {
//		mul_fp(out[idx], a[idx], b[idx]);
//	}
//}
//// Vectorized multiplication of Fp elements (GPU)
//BLST_INLINE void mul_fp_vec_cu(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//	size_t threads_per_block = 256;
//	size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
//	mul_fp_vec_ker << <num_blocks, threads_per_block >> > (out, a, b, num);
//	cudaDeviceSynchronize();
//}
//#endif
//// Vectorized multiplication of Fp elements
//BLST_INLINE void mul_fp_vec(vec384* out, const vec384* a, const vec384* b, size_t num)
//{
//#ifdef BLST_CUDA
//	mul_fp_vec_cu(out, a, b, num);
//#else
//	mul_fp_vec_cc(out, a, b, num);
//#endif
//}
//
//// Vectorized squaring of Fp elements (CPU)
//BLST_INLINE void sqr_fp_vec_cc(vec384* out, const vec384* a, size_t num)
//{
//	while (num--) {
//		sqr_fp(*out, *a);
//		out++; a++;
//	}
//}
//#ifdef BLST_CUDA
//// kernel for squaring of Fp elements (GPU)
//GLOBAL_FUNC BLST_INLINE void sqr_fp_vec_ker(vec384* out, const vec384* a, size_t num)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	if (idx < num) {
//		sqr_fp(out[idx], a[idx]);
//	}
//}
//// Vectorized squaring of Fp elements (GPU)
//BLST_INLINE void sqr_fp_vec_cu(vec384* out, const vec384* a, size_t num)
//{
//	size_t threads_per_block = 256;
//	size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
//	sqr_fp_vec_ker << <num_blocks, threads_per_block >> > (out, a, num);
//	cudaDeviceSynchronize();
//}
//#endif
//// Vectorized squaring of Fp elements
//BLST_INLINE void sqr_fp_vec(vec384* out, const vec384* a, size_t num)
//{
//#ifdef BLST_CUDA
//	sqr_fp_vec_cu(out, a, num);
//#else
//	sqr_fp_vec_cc(out, a, num);
//#endif
//}
//
//// Vectorized negation of Fp elements (CPU)
//BLST_INLINE void cneg_fp_vec_cc(vec384* out, const vec384* a, int flag, size_t num)
//{
//	while (num--) {
//		cneg_fp(*out, *a, flag);
//		out++; a++;
//	}
//}
//#ifdef BLST_CUDA
//// kernel for negation of Fp elements (GPU)
//GLOBAL_FUNC BLST_INLINE void cneg_fp_vec_ker(vec384* out, const vec384* a, int flag, size_t num)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	if (idx < num) {
//		cneg_fp(out[idx], a[idx], flag);
//	}
//}
//// Vectorized negation of Fp elements (GPU)
//BLST_INLINE void cneg_fp_vec_cu(vec384* out, const vec384* a, int flag, size_t num)
//{
//	size_t threads_per_block = 256;
//	size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
//	cneg_fp_vec_ker << <num_blocks, threads_per_block >> > (out, a, flag, num);
//	cudaDeviceSynchronize();
//}
//#endif
//
//// Vectorized negation of Fp elements
//BLST_INLINE void cneg_fp_vec(vec384* out, const vec384* a, int flag, size_t num)
//{
//#ifdef BLST_CUDA
//	cneg_fp_vec_cu(out, a, flag, num);
//#else
//	cneg_fp_vec_cc(out, a, flag, num);
//#endif
//}
//
//
//
//// mult by 3
//
//// Vectorized multiplication of Fp elements by 3 (CPU)
//BLST_INLINE void mul3_fp_vec_cc(vec384* out, const vec384* a, size_t num)
//{
//	while (num--) {
//		mul3_fp(*out, *a);
//		out++; a++;
//	}
//}
//#ifdef BLST_CUDA
//// kernel for multiplication of Fp elements by 3 (GPU)
//GLOBAL_FUNC BLST_INLINE void mul3_fp_vec_ker(vec384* out, const vec384* a, size_t num)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	if (idx < num) {
//		mul3_fp(out[idx], a[idx]);
//	}
//}
//// Vectorized multiplication of Fp elements by 3 (GPU)
//BLST_INLINE void mul3_fp_vec_cu(vec384* out, const vec384* a, size_t num)
//{
//	size_t threads_per_block = 256;
//	size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
//	mul3_fp_vec_ker << <num_blocks, threads_per_block >> > (out, a, num);
//	cudaDeviceSynchronize();
//}
//#endif
//// Vectorized multiplication of Fp elements by 3
//BLST_INLINE void mul3_fp_vec(vec384* out, const vec384* a, size_t num)
//{
//#ifdef BLST_CUDA
//	mul3_fp_vec_cu(out, a, num);
//#else
//	mul3_fp_vec_cc(out, a, num);
//#endif
//}
//
//// mul by 8 
//
//// Vectorized multiplication of Fp elements by 8 (CPU)
//BLST_INLINE void mul8_fp_vec_cc(vec384* out, const vec384* a, size_t num)
//{
//	while (num--) {
//		mul8_fp(*out, *a);
//		out++; a++;
//	}
//}
//#ifdef BLST_CUDA
//// kernel for multiplication of Fp elements by 8 (GPU)
//GLOBAL_FUNC BLST_INLINE void mul8_fp_vec_ker(vec384* out, const vec384* a, size_t num)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	if (idx < num) {
//		mul8_fp(out[idx], a[idx]);
//	}
//}
//// Vectorized multiplication of Fp elements by 8 (GPU)
//BLST_INLINE void mul8_fp_vec_cu(vec384* out, const vec384* a, size_t num)
//{
//	size_t threads_per_block = 256;
//	size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
//	mul8_fp_vec_ker << <num_blocks, threads_per_block >> > (out, a, num);
//	cudaDeviceSynchronize();
//}
//#endif
//// Vectorized multiplication of Fp elements by 8
//BLST_INLINE void mul8_fp_vec(vec384* out, const vec384* a, size_t num)
//{
//#ifdef BLST_CUDA
//	mul8_fp_vec_cu(out, a, num);
//#else
//	mul8_fp_vec_cc(out, a, num);
//#endif
//}
//
//
//// div by 2
//
//// Vectorized division of Fp elements by 2 (CPU)
//BLST_INLINE void div2_fp_vec_cc(vec384* out, const vec384* a
//	, size_t num)
//{
//	while (num--) {
//		div2_fp(*out, *a);
//		out++; a++;
//	}
//}
//#ifdef BLST_CUDA
//// kernel for division of Fp elements by 2 (GPU)
//GLOBAL_FUNC BLST_INLINE void div2_fp_vec_ker(vec384* out, const vec
//	384 * a, size_t num)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	if (idx < num) {
//		div2_fp(out[idx], a[idx]);
//	}
//}
//// Vectorized division of Fp elements by 2 (GPU)
//BLST_INLINE void div2_fp_vec_cu(vec384* out, const vec384*
//	a, size_t num)
//{
//	size_t threads_per_block = 256;
//	size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
//	div2_fp_vec_ker << <num_blocks, threads_per_block >> > (out, a, num);
//	cudaDeviceSynchronize();
//}
//#endif
//// Vectorized division of Fp elements by 2
//BLST_INLINE void div2_fp_vec(vec384* out, const vec384* a
//	, size_t num)
//{
//#ifdef BLST_CUDA
//	div2_fp_vec_cu(out, a, num);
//#else
//	div2_fp_vec_cc(out, a, num);
//#endif
//}
//
//// lshift 
//
//// Vectorized left shift of Fp elements (CPU)
//BLST_INLINE void lshift_fp_vec_cc(vec384* out, const vec384* a
//	, size_t count, size_t num)
//{
//	while (num--) {
//		lshift_fp(*out, *a, count);
//		out++; a++;
//	}
//}
//
//#ifdef BLST_CUDA
//// kernel for left shift of Fp elements (GPU)
//GLOBAL_FUNC BLST_INLINE void lshift_fp_vec_ker(vec384* out, const vec
//	384 * a, size_t count, size_t num)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	if (idx < num) {
//		lshift_fp(out[idx], a[idx], count);
//	}
//}
//// Vectorized left shift of Fp elements (GPU)
//BLST_INLINE void lshift_fp_vec_cu(vec384* out, const vec384*
//	a, size_t count, size_t num)
//{
//	size_t threads_per_block = 256;
//	size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
//	lshift_fp_vec_ker << <num_blocks, threads_per_block >> > (out, a, count, num);
//	cudaDeviceSynchronize();
//}
//#endif
//// Vectorized left shift of Fp elements
//BLST_INLINE void lshift_fp_vec(vec384* out, const vec384* a
//	, size_t count, size_t num)
//{
//#ifdef BLST_CUDA
//	lshift_fp_vec_cu(out, a, count, num);
//#else
//	lshift_fp_vec_cc(out, a, count, num);
//#endif
//}
//
//
//// to fp
//
//// Vectorized conversion from a regular representation to a montgomery Fp (CPU)
//BLST_INLINE void to_fp_vec_cc(vec384* out, const vec384* a, size_t num)
//{
//	while (num--) {
//		to_fp(*out, *a);
//		out++; a++;
//	}
//}
//#ifdef BLST_CUDA
//// kernel for conversion from a regular representation to a montgomery Fp (GPU)
//GLOBAL_FUNC BLST_INLINE void to_fp_vec_ker(vec384* out, const vec384
//	* a, size_t num)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	if (idx < num) {
//		to_fp(out[idx], a[idx]);
//	}
//}
//// Vectorized conversion from a regular representation to a montgomery Fp (GPU)
//BLST_INLINE void to_fp_vec_cu(vec384* out, const vec384* a
//	, size_t num)
//{
//	size_t threads_per_block = 256;
//	size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
//	to_fp_vec_ker << <num_blocks, threads_per_block >> > (out, a, num);
//	cudaDeviceSynchronize();
//}
//#endif
//// Vectorized conversion from a regular representation to a montgomery Fp
//BLST_INLINE void to_fp_vec(vec384* out, const vec384* a
//	, size_t num)
//{
//#ifdef BLST_CUDA
//	to_fp_vec_cu(out, a, num);
//#else
//	to_fp_vec_cc(out, a, num);
//#endif
//}
//
//// from fp
//
//// Vectorized conversion from a montgomery Fp to a regular representation (CPU)
//BLST_INLINE void from_fp_vec_cc(vec384* out, const vec384* a,
//	size_t num)
//{
//	while (num--) {
//		from_fp(*out, *a);
//		out++; a++;
//	}
//}
//#ifdef BLST_CUDA
//// kernel for conversion from a montgomery Fp to a regular representation (GPU)
//GLOBAL_FUNC BLST_INLINE void from_fp_vec_ker(vec384* out, const vec384
//	* a, size_t num)
//{
//	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//	if (idx < num) {
//		from_fp(out[idx], a[idx]);
//	}
//}
//// Vectorized conversion from a montgomery Fp to a regular representation (GPU)
//BLST_INLINE void from_fp_vec_cu(vec384* out, const vec384*
//	a, size_t num)
//{
//	size_t threads_per_block = 256;
//	size_t num_blocks = (num + threads_per_block - 1) / threads_per_block
//		;
//	from_fp_vec_ker << <num_blocks, threads_per_block >> > (out, a
//		, num);
//	cudaDeviceSynchronize();
//}
//#endif
//// Vectorized conversion from a montgomery Fp to a regular representation
//BLST_INLINE void from_fp_vec(vec384* out, const vec384* a,
//	size_t num)
//{
//#ifdef BLST_CUDA
//	from_fp_vec_cu(out, a, num);
//#else
//	from_fp_vec_cc(out, a, num);
//#endif
//}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//BLST_EXTERN_C_END
//
