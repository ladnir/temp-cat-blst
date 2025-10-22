/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/curve/point.h"
#include "blst/fields/Fp12.h"
#include "pairing.h"
#include "blst/curve/GT_cuda.h"
#include "blst/exports.h"
#include "blst/common/errors.h"
#include "blst/sign/aggregate.h"

#ifdef BLST_GT_KERNELS
#include "blst/common/cuda_launch.h"
#endif

namespace blst
{

	///////////////////////////////////////////
	// Core GT Operations

	// Vectorized multiplication of GT elements (CPU)
	batched_error_code mult_gt_vec_cc(Fp12* out, const Fp12* a, const Fp12* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			Fp12::mul(*out, *a, *b);
			out++; a++; b++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for multiplication of GT elements (GPU)
	__global__ void mult_gt_vec_ker(Fp12* out, const Fp12* a, const Fp12* b, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp12::mul(out[idx], a[idx], b[idx]);
		}
	}

	// Vectorized multiplication of GT elements (GPU)
	batched_error_code mult_gt_vec_cu(Fp12* out, const Fp12* a, const Fp12* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(mult_gt_vec_ker, num, out, a, b);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized multiplication of GT elements
	batched_error_code mult_gt_vec(Fp12* out, const Fp12* a, const Fp12* b, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return mult_gt_vec_cu(out, a, b, num);
#else
		return mult_gt_vec_cc(out, a, b, num);
#endif
	}

	// Vectorized squaring of GT elements (CPU)
	batched_error_code square_gt_vec_cc(Fp12* out, const Fp12* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			Fp12::sqr(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for squaring of GT elements (GPU)
	__global__ void square_gt_vec_ker(Fp12* out, const Fp12* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp12::sqr(out[idx], a[idx]);
		}
	}

	// Vectorized squaring of GT elements (GPU)
	batched_error_code square_gt_vec_cu(Fp12* out, const Fp12* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(square_gt_vec_ker, num, out, a);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized squaring of GT elements
	batched_error_code square_gt_vec(Fp12* out, const Fp12* a, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return square_gt_vec_cu(out, a, num);
#else
		return square_gt_vec_cc(out, a, num);
#endif
	}

	// Vectorized cyclotomic squaring of GT elements (CPU)
	batched_error_code cyclotomic_square_gt_vec_cc(Fp12* out, const Fp12* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			Fp12::cyclotomic_sqr(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for cyclotomic squaring of GT elements (GPU)
	__global__ void cyclotomic_square_gt_vec_ker(Fp12* out, const Fp12* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp12::cyclotomic_sqr(out[idx], a[idx]);
		}
	}

	// Vectorized cyclotomic squaring of GT elements (GPU)
	batched_error_code cyclotomic_square_gt_vec_cu(Fp12* out, const Fp12* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(cyclotomic_square_gt_vec_ker, num, out, a);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized cyclotomic squaring of GT elements
	batched_error_code cyclotomic_square_gt_vec(Fp12* out, const Fp12* a, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return cyclotomic_square_gt_vec_cu(out, a, num);
#else
		return cyclotomic_square_gt_vec_cc(out, a, num);
#endif
	}

	// Vectorized inversion of GT elements (CPU)
	batched_error_code inverse_gt_vec_cc(Fp12* out, const Fp12* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			Fp12::inverse(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for inversion of GT elements (GPU)
	__global__ void inverse_gt_vec_ker(Fp12* out, const Fp12* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp12::inverse(out[idx], a[idx]);
		}
	}

	// Vectorized inversion of GT elements (GPU)
	batched_error_code inverse_gt_vec_cu(Fp12* out, const Fp12* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(inverse_gt_vec_ker, num, out, a);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized inversion of GT elements
	batched_error_code inverse_gt_vec(Fp12* out, const Fp12* a, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return inverse_gt_vec_cu(out, a, num);
#else
		return inverse_gt_vec_cc(out, a, num);
#endif
	}

	// Vectorized conjugation of GT elements (CPU)
	batched_error_code conjugate_gt_vec_cc(Fp12* out, const Fp12* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			vec_copy(out, a, sizeof(Fp12));
			Fp12::conjugate(*out);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for conjugation of GT elements (GPU)
	__global__ void conjugate_gt_vec_ker(Fp12* out, const Fp12* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			vec_copy(&out[idx], &a[idx], sizeof(Fp12));
			Fp12::conjugate(out[idx]);
		}
	}

	// Vectorized conjugation of GT elements (GPU)
	batched_error_code conjugate_gt_vec_cu(Fp12* out, const Fp12* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(conjugate_gt_vec_ker, num, out, a);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized conjugation of GT elements
	batched_error_code conjugate_gt_vec(Fp12* out, const Fp12* a, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return conjugate_gt_vec_cu(out, a, num);
#else
		return conjugate_gt_vec_cc(out, a, num);
#endif
	}

	// Vectorized exponentiation of GT elements (CPU)
	batched_error_code pow_gt_vec_cc(Fp12* out, const Fp12* bases, const unsigned char* exponents,
		size_t exp_bits, size_t num)
	{
		if (!out || !bases || !exponents) return BATCHED_ERROR_NULL_POINTER;

		const size_t exp_bytes = (exp_bits + 7) / 8;

		while (num--) {
			Fp12::exp(*out, *bases, (limb_t*)exponents, exp_bits);
			out++; bases++; exponents += exp_bytes;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for exponentiation of GT elements (GPU)
	__global__ void pow_gt_vec_ker(Fp12* out, const Fp12* bases, const unsigned char* exponents,
		size_t exp_bits, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			const size_t exp_bytes = (exp_bits + 7) / 8;
			Fp12::exp(out[idx], bases[idx], (limb_t*)&exponents[idx * exp_bytes], exp_bits);
		}
	}

	// Vectorized exponentiation of GT elements (GPU)
	batched_error_code pow_gt_vec_cu(Fp12* out, const Fp12* bases, const unsigned char* exponents,
		size_t exp_bits, size_t num)
	{
		if (!out || !bases || !exponents) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(pow_gt_vec_ker, num, out, bases, exponents, exp_bits);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized exponentiation of GT elements
	batched_error_code pow_gt_vec(Fp12* out, const Fp12* bases, const unsigned char* exponents,
		size_t exp_bits, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return pow_gt_vec_cu(out, bases, exponents, exp_bits, num);
#else
		return pow_gt_vec_cc(out, bases, exponents, exp_bits, num);
#endif
	}

	///////////////////////////////////////////
	// Validation Operations

	// Vectorized check if GT elements equal identity (CPU)
	batched_error_code is_one_gt_vec_cc(bool* results, const Fp12* elements, size_t num)
	{
		if (!results || !elements) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			*results = Fp12::is_one(*elements);
			results++; elements++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for identity check of GT elements (GPU)
	__global__ void is_one_gt_vec_ker(bool* results, const Fp12* elements, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			results[idx] = Fp12::is_one(elements[idx]);
		}
	}

	// Vectorized check if GT elements equal identity (GPU)
	batched_error_code is_one_gt_vec_cu(bool* results, const Fp12* elements, size_t num)
	{
		if (!results || !elements) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(is_one_gt_vec_ker, num, results, elements);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized check if GT elements equal identity
	batched_error_code is_one_gt_vec(bool* results, const Fp12* elements, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return is_one_gt_vec_cu(results, elements, num);
#else
		return is_one_gt_vec_cc(results, elements, num);
#endif
	}

	// Vectorized GT element equality check (CPU)
	batched_error_code is_equal_gt_vec_cc(bool* results, const Fp12* a, const Fp12* b, size_t num)
	{
		if (!results || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			*results = Fp12::equal(*a, *b);
			results++; a++; b++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for GT element equality check (GPU)
	__global__ void is_equal_gt_vec_ker(bool* results, const Fp12* a, const Fp12* b, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			results[idx] = Fp12::equal(a[idx], b[idx]);
		}
	}

	// Vectorized GT element equality check (GPU)
	batched_error_code is_equal_gt_vec_cu(bool* results, const Fp12* a, const Fp12* b, size_t num)
	{
		if (!results || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(is_equal_gt_vec_ker, num, results, a, b);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized GT element equality check
	batched_error_code is_equal_gt_vec(bool* results, const Fp12* a, const Fp12* b, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return is_equal_gt_vec_cu(results, a, b, num);
#else
		return is_equal_gt_vec_cc(results, a, b, num);
#endif
	}

	// Vectorized GT group membership check (CPU)
	batched_error_code in_group_gt_vec_cc(bool* results, const Fp12* elements, size_t num)
	{
		if (!results || !elements) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			*results = Fp12::in_group(*elements);
			results++; elements++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for GT group membership check (GPU)
	__global__ void in_group_gt_vec_ker(bool* results, const Fp12* elements, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			results[idx] = Fp12::in_group(elements[idx]);
		}
	}

	// Vectorized GT group membership check (GPU)
	batched_error_code in_group_gt_vec_cu(bool* results, const Fp12* elements, size_t num)
	{
		if (!results || !elements) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(in_group_gt_vec_ker, num, results, elements);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized GT group membership check
	batched_error_code in_group_gt_vec(bool* results, const Fp12* elements, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return in_group_gt_vec_cu(results, elements, num);
#else
		return in_group_gt_vec_cc(results, elements, num);
#endif
	}

	///////////////////////////////////////////
	// Pairing-related Operations

	// Vectorized Miller loop computation (CPU)
	batched_error_code miller_loop_gt_vec_cc(Fp12* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		if (!out || !g2_points || !g1_points) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			blst_miller_loop(*out, g2_points, g1_points);
			out++; g2_points++; g1_points++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for Miller loop computation (GPU)
	__global__ void miller_loop_gt_vec_ker(Fp12* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			blst_miller_loop(out[idx], &g2_points[idx], &g1_points[idx]);
		}
	}

	// Vectorized Miller loop computation (GPU)
	batched_error_code miller_loop_gt_vec_cu(Fp12* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		if (!out || !g2_points || !g1_points) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(miller_loop_gt_vec_ker, num, out, g2_points, g1_points);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized Miller loop computation
	batched_error_code miller_loop_gt_vec(Fp12* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return miller_loop_gt_vec_cu(out, g2_points, g1_points, num);
#else
		return miller_loop_gt_vec_cc(out, g2_points, g1_points, num);
#endif
	}

	// Vectorized final exponentiation (CPU)
	batched_error_code final_exp_gt_vec_cc(Fp12* out, const Fp12* miller_results, size_t num)
	{
		if (!out || !miller_results) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			Fp12::final_exp(*out, *miller_results);
			out++; miller_results++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for final exponentiation (GPU)
	__global__ void final_exp_gt_vec_ker(Fp12* out, const Fp12* miller_results, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp12::final_exp(out[idx], miller_results[idx]);
		}
	}

	// Vectorized final exponentiation (GPU)
	batched_error_code final_exp_gt_vec_cu(Fp12* out, const Fp12* miller_results, size_t num)
	{
		if (!out || !miller_results) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(final_exp_gt_vec_ker, num, out, miller_results);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized final exponentiation
	batched_error_code final_exp_gt_vec(Fp12* out, const Fp12* miller_results, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return final_exp_gt_vec_cu(out, miller_results, num);
#else
		return final_exp_gt_vec_cc(out, miller_results, num);
#endif
	}

	// Vectorized complete pairing computation (CPU)
	batched_error_code pairing_gt_vec_cc(Fp12* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		if (!out || !g2_points || !g1_points) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			blst_miller_loop(*out, g2_points, g1_points);
			Fp12::final_exp(*out, *out);
			out++; g2_points++; g1_points++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for complete pairing computation (GPU)
	__global__ void pairing_gt_vec_ker(Fp12* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			blst_miller_loop(out[idx], &g2_points[idx], &g1_points[idx]);
			Fp12::final_exp(out[idx], out[idx]);
		}
	}

	// Vectorized complete pairing computation (GPU)
	batched_error_code pairing_gt_vec_cu(Fp12* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		if (!out || !g2_points || !g1_points) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(pairing_gt_vec_ker, num, out, g2_points, g1_points);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized complete pairing computation
	batched_error_code pairing_gt_vec(Fp12* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return pairing_gt_vec_cu(out, g2_points, g1_points, num);
#else
		return pairing_gt_vec_cc(out, g2_points, g1_points, num);
#endif
	}

	///////////////////////////////////////////
	// Serialization Operations

	// Vectorized GT element serialization (CPU)
	batched_error_code serialize_gt_vec_cc(unsigned char* out, const Fp12* elements, size_t num)
	{
		if (!out || !elements) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			Fp12::to_bendian(out, *elements);
			out += 576; elements++; // 12 * 48 bytes per Fp12 element
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for GT element serialization (GPU)
	__global__ void serialize_gt_vec_ker(unsigned char* out, const Fp12* elements, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp12::to_bendian(&out[idx * 576], elements[idx]);
		}
	}

	// Vectorized GT element serialization (GPU)
	batched_error_code serialize_gt_vec_cu(unsigned char* out, const Fp12* elements, size_t num)
	{
		if (!out || !elements) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(serialize_gt_vec_ker, num, out, elements);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized GT element serialization
	batched_error_code serialize_gt_vec(unsigned char* out, const Fp12* elements, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return serialize_gt_vec_cu(out, elements, num);
#else
		return serialize_gt_vec_cc(out, elements, num);
#endif
	}

	// Vectorized GT element deserialization (CPU)
	batched_error_code deserialize_gt_vec_cc(Fp12* out, const unsigned char* serialized, size_t num)
	{
		if (!out || !serialized) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			Fp12::from_bendian(*out, serialized);
			out++; serialized += 576;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for GT element deserialization (GPU)
	__global__ void deserialize_gt_vec_ker(Fp12* out, const unsigned char* serialized, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp12::from_bendian(out[idx], &serialized[idx * 576]);
		}
	}

	// Vectorized GT element deserialization (GPU)
	batched_error_code deserialize_gt_vec_cu(Fp12* out, const unsigned char* serialized, size_t num)
	{
		if (!out || !serialized) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(deserialize_gt_vec_ker, num, out, serialized);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized GT element deserialization
	batched_error_code deserialize_gt_vec(Fp12* out, const unsigned char* serialized, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return deserialize_gt_vec_cu(out, serialized, num);
#else
		return deserialize_gt_vec_cc(out, serialized, num);
#endif
	}

	///////////////////////////////////////////
	// Advanced Operations

//	// Multi-exponentiation for GT elements (CPU) - simplified version
//	// result = prod_i (bases[i] ^ exponents[i])
//	batched_error_code multi_exp_gt_cc(Fp12* result, const Fp12* bases,
//		const unsigned char* exponents, size_t exp_bits, size_t num)
//	{
//		if (!result || !bases || !exponents) return BATCHED_ERROR_NULL_POINTER;
//
//		// result = 1
//		Fp12::one_copy(*result);
//
//		const size_t exp_bytes = (exp_bits + 7) / 8;
//
//		// Simple multi-exponentiation: accumulate exponentiations
//		for (size_t i = 0; i < num; i++) {
//			Fp12 temp;
//
//			// temp = bases[i] ^ exponents[i]
//			Fp12::pow(temp, bases[i], (limb_t*)&exponents[i * exp_bytes], exp_bits);
//
//			// result *= temp
//			if (i == 0) {
//				vec_copy(*result, temp, sizeof(temp));
//			}
//			else {
//				Fp12::mul(*result, *result, temp);
//			}
//		}
//
//
//		return BATCHED_SUCCESS;
//	}
//#ifdef __INTELLISENSE__
//	// Disable IntelliSense warning about __shfl_down_sync unknown
//#define __shfl_down_sync(mask, var, offset) (var)
//#define __syncthreads()
//#endif
//
//#ifdef BLST_GT_KERNELS
//
//	// --- helpers ---------------------------------------------------------------
//
//#ifndef CUDA_OK
//#define CUDA_OK(stmt) do { \
//		cudaError_t _e = (stmt); \
//		if (_e != cudaSuccess) {\
//			printf("CUDA error %d at %s:%d\n", (int)_e, __FILE__, __LINE__);\
//			printf("  %s\n", cudaGetErrorString(_e));\
//			return cuda_to_batched(_e);\
//		} } while (0)
//#endif
//
//
//// ---------- helpers ----------
//	// --- add near the top of your TU ---
//	static inline void dbg_ptr(const void* p, const char* name) {
//#if !defined(NDEBUG)
//		cudaPointerAttributes a;
//		auto s = cudaPointerGetAttributes(&a, p);
//		if (s == cudaSuccess) {
//			fprintf(stderr, "[%s] type=%d (1=host,2=device,3=managed,4=unified)\n", name, (int)a.type);
//		}
//		else {
//			fprintf(stderr, "[%s] cudaPointerGetAttributes err=%d\n", name, (int)s);
//		}
//#endif
//	}
//	static inline void dbg_kernel_attrs(const void* k, const char* name) {
//#if !defined(NDEBUG)
//		cudaFuncAttributes fa{};
//		if (cudaFuncGetAttributes(&fa, k) == cudaSuccess) {
//			fprintf(stderr,
//				"[ker %s] regs=%d, staticSmem=%zu, local=%zu, maxTPB=%d, bin=%d, ptx=%d\n",
//				name, fa.numRegs, (size_t)fa.sharedSizeBytes, (size_t)fa.localSizeBytes,
//				fa.maxThreadsPerBlock, fa.binaryVersion, fa.ptxVersion);
//		}
//		else {
//			fprintf(stderr, "[ker %s] cudaFuncGetAttributes failed\n", name);
//		}
//#endif
//	}
//
//	struct fp12_words64 { static constexpr int kWords = (int)(sizeof(Fp12) / sizeof(uint64_t)); uint64_t w[kWords]; };
//
//	// shuffle directly from src to dst (only one 576B tmp owned by caller)
//	__device__ __forceinline__
//		void shfl_down_fp12(Fp12& dst, const Fp12& src, int offset, unsigned mask = 0xFFFFFFFFu)
//	{
//		// treat arrays as contiguous bytes; then as 64-bit lanes
//		const uint64_t* __restrict__ sw = reinterpret_cast<const uint64_t*>(src); // src decays
//		uint64_t* __restrict__       dw = reinterpret_cast<uint64_t*>(dst);       // dst decays
//#pragma unroll
//		for (int i = 0; i < fp12_words64::kWords; ++i) {
//			dw[i] = __shfl_down_sync(mask, sw[i], offset);
//		}
//	}
//
//	__device__ __forceinline__
//		void warp_reduce_mul(Fp12& acc)
//	{
//		unsigned mask = 0xFFFFFFFFu;
//		for (int off = 16; off > 0; off >>= 1) {
//			Fp12 tmp;                    // one 576B buffer reused per step
//			shfl_down_fp12(tmp, acc, off, mask);
//			Fp12::mul(acc, acc, tmp);
//		}
//	}
//
//	// Avoid bank conflicts with a tiny pad; keep alignment tame on MSVC.
//	//struct fp12_padded { Fp12 v; uint32_t _pad; };
//	// ===== conflict-free shared-memory record for fp12 =====
//
//// we want (sizeof(fp12_padded)/4) % 32 == 1
//	constexpr size_t kFP12 = sizeof(Fp12);       // you reported 576
//	constexpr size_t kStride4 = kFP12 / 4;             // 144
//	constexpr size_t kPad4 = (32 + 1 - (kStride4 % 32)) % 32; // (33 - 16) % 32 = 17
//	constexpr size_t kPadBytes = kPad4 * 4;            // 68
//
//	// force *small* member alignment so compiler won't add surprise tail padding
//#pragma pack(push, 4)
//	struct fp12_padded {
//		uint8_t bytes[kFP12];
//		uint8_t pad[kPadBytes];                        // 68 bytes -> total 644
//	};
//#pragma pack(pop)
//
//	static_assert(sizeof(fp12_padded) == kFP12 + kPadBytes, "unexpected fp12_padded size");
//	static_assert(((sizeof(fp12_padded) / 4) % 32) == 1, "bad shared-mem stride (bank conflicts)");
//
//	// Pack exponent bytes (LE) into 64b limb_ts (<= 1024 bits here; bump if needed)
//	template<int Maxlimb_ts = 16>
//	__device__ __forceinline__ int load_exp_le(uint64_t(&dst)[Maxlimb_ts], const uint8_t* src, size_t exp_bits) {
//		const size_t exp_bytes = (exp_bits + 7u) / 8u;
//		const int limb_t_cnt = (int)((exp_bits + 63u) / 64u);
//#pragma unroll
//		for (int i = 0; i < Maxlimb_ts; ++i) dst[i] = 0ull;
//		int i = 0;
//		for (; i + 8 <= (int)exp_bytes && (i / 8) < Maxlimb_ts; i += 8) { uint64_t v; memcpy(&v, src + i, 8); dst[i / 8] = v; }
//		if (i < (int)exp_bytes && (i / 8) < Maxlimb_ts) {
//			uint64_t v = 0; int tail = (int)exp_bytes - i; uint8_t* pv = (uint8_t*)&v;
//#pragma unroll
//			for (int b = 0; b < tail; ++b) pv[b] = src[i + b];
//			dst[i / 8] = v;
//		}
//		return limb_t_cnt;
//	}
//
//	// ---------- KERNELS (NO DYNAMIC SMEM) ----------
//
//	// One fp12 per warp goes into a fixed shared buffer (max 32 warps/block).
//	__global__ void multi_exp_gt_ker(Fp12* __restrict__ result,
//		const Fp12* __restrict__ bases,
//		const unsigned char* __restrict__ exponents,
//		size_t exp_bits, size_t num)
//	{
//		__shared__ fp12_padded warp_buf[32];  // static ~32 * (576+4) = ~18.9 KiB
//
//		const int tid = threadIdx.x;
//		const int lane = tid & 31;
//		const int wid = tid >> 5;
//		const int warps_per_block = (blockDim.x + 31) / 32;
//		const size_t idx = (size_t)blockIdx.x * blockDim.x + tid;
//		const size_t exp_bytes = (exp_bits + 7u) / 8u;
//
//		Fp12 acc;
//		if (idx < num) {
//			uint64_t limb_ts[16];
//			(void)load_exp_le(limb_ts, exponents + idx * exp_bytes, exp_bits);
//			Fp12::pow(acc, bases[idx], (const limb_t*)limb_ts, exp_bits);
//		}
//		else {
//			Fp12::one_copy(acc);
//		}
//
//		warp_reduce_mul(acc);
//
//		if (lane == 0)
//		{
//			//warp_buf[wid].v = acc;
//			memcpy(warp_buf[wid].bytes, acc, sizeof(acc));
//		}
//		__syncthreads();
//
//		if (wid == 0) {
//			Fp12 r; Fp12::one_copy(r);
//			if (lane < warps_per_block)
//			{
//				//r = warp_buf[lane].v;
//				memcpy(r, warp_buf[lane].bytes, sizeof(r));
//			}
//			warp_reduce_mul(r);
//			if (lane == 0) {
//				//result[blockIdx.x] = r;
//				vec_copy(result[blockIdx.x], r, sizeof(r));
//			}
//		}
//	}
//
//	// Reduce N partials down to one (grid-stride, no dyn smem)
//	__global__ void multi_exp_final_reduce_ker(Fp12* __restrict__ out,
//		const Fp12* __restrict__ in,
//		size_t n)
//	{
//		__shared__ fp12_padded warp_buf[32];
//
//		const int tid = threadIdx.x;
//		const int lane = tid & 31;
//		const int wid = tid >> 5;
//		const int warps_per_block = (blockDim.x + 31) / 32;
//
//		Fp12 acc; Fp12::one_copy(acc);
//		for (size_t i = tid; i < n; i += blockDim.x) {
//			Fp12 v /*= in[i]*/;
//			vec_copy(v, in[i], sizeof(v));
//			Fp12::mul(acc, acc, v);
//		}
//
//		warp_reduce_mul(acc);
//
//		if (lane == 0) {
//			//warp_buf[wid].v = acc;
//			memcpy(warp_buf[wid].bytes, acc, sizeof(acc));
//		}
//		__syncthreads();
//
//		if (wid == 0) {
//			Fp12 r; Fp12::one_copy(r);
//			if (lane < warps_per_block) {
//				// r = warp_buf[lane].v;
//				memcpy(r, warp_buf[lane].bytes, sizeof(r));
//			}
//			warp_reduce_mul(r);
//			if (lane == 0) {
//				//*out = r;
//				vec_copy(*out, r, sizeof(r));
//			}
//		}
//	}
//
//	// ---------- HOST API ----------
//	batched_error_code multi_exp_gt_cu(Fp12* result,
//		const Fp12* bases,
//		const unsigned char* exponents,
//		size_t exp_bits,
//		size_t num)
//	{
//		if (!result || !bases || !exponents) return BATCHED_ERROR_NULL_POINTER;
//
//		auto l = launch_auto(multi_exp_gt_ker, num,
//			(Fp12*) nullptr,
//			(Fp12*) nullptr,
//			(const unsigned char*)nullptr,
//			exp_bits);
//		auto e = l.get();
//		if( e != cudaSuccess ) {
//			return cuda_to_batched(e);
//		}
//		return BATCHED_SUCCESS;
//		//// Trivial cases
//		//if (num == 0) {
//		//	Fp12 one; Fp12::one_copy(one);
//		//	CUDA_OK(cudaMemcpy(result, &one, sizeof(one), cudaMemcpyHostToDevice));
//		//	return BATCHED_SUCCESS;
//		//}
//
//		//// Give threads some stack (the debug launcher also does 16 KiB, but we can be explicit)
//		////CUDA_OK(cudaDeviceSetLimit(cudaLimitStackSize, 32 * 1024));
//
//		//// Phase 1 config ----------------------------------------------------------
//		//// Suggest a block size from occupancy; clamp to at least 32
//		//int bs1 = suggest_block_size(multi_exp_gt_ker, /*dyn_smem=*/0);
//		//if (bs1 < 32) bs1 = 32;
//
//		//// If num is tiny, avoid oversizing the block
//		//if ((size_t)bs1 > num) bs1 = (int)num;
//		//bs1 = (bs1 + 31) / 32 * 32;
//
//		//int num_blocks = (int)((num + (size_t)bs1 - 1) / (size_t)bs1);
//		//if (num_blocks <= 0) num_blocks = 1;
//
//
//		//// Scratch for per-block partial reductions
//		//Fp12* d_partial = nullptr;
//		//CUDA_OK(cudaMalloc(&d_partial, (size_t)num_blocks * sizeof(Fp12)));
//
//		//// Launch phase 1 via the debug launcher
//		//{
//		//	cudaError_t e = launch_checked(
//		//		multi_exp_gt_ker,
//		//		dim3((u32)num_blocks, 1, 1),
//		//		dim3((u32)bs1, 1, 1),
//		//		/*dyn_smem=*/0,
//		//		/*stream=*/0,
//		//		"multi_exp_gt_ker",
//		//		d_partial, bases, exponents, exp_bits, num
//		//	);
//		//	if (e != cudaSuccess) {
//		//		cudaFree(d_partial);
//		//		return cuda_to_batched(e);
//		//	}
//		//}
//
//		//// Phase 2 (final reduction) ----------------------------------------------
//		//if (num_blocks == 1) {
//		//	// Nothing to reduce; just copy the single block result
//		//	CUDA_OK(cudaMemcpy(result, d_partial, sizeof(Fp12), cudaMemcpyDeviceToDevice));
//		//	CUDA_OK(cudaFree(d_partial));
//		//	return BATCHED_SUCCESS;
//		//}
//
//		//// Suggest block size for the final reduce kernel (grid=1)
//		//int bs2 = suggest_block_size(multi_exp_final_reduce_ker, /*dyn_smem=*/0);
//		//if (bs2 < 32) bs2 = 32;
//
//		//{
//		//	cudaError_t e = launch_checked(
//		//		multi_exp_final_reduce_ker,
//		//		dim3(1, 1, 1),
//		//		dim3((unsigned)bs2, 1, 1),
//		//		/*dyn_smem=*/0,
//		//		/*stream=*/0,
//		//		"multi_exp_final_reduce_ker",
//		//		result, d_partial, (size_t)num_blocks
//		//	);
//		//	if (e != cudaSuccess) {
//		//		cudaFree(d_partial);
//		//		return cuda_to_batched(e);
//		//	}
//		//}
//
//		//CUDA_OK(cudaFree(d_partial));
//		//return BATCHED_SUCCESS;
//	}
//
//#endif // BLST_CUDA
//
//
//
//		// Multi-exponentiation for GT elements
//	batched_error_code multi_exp_gt(Fp12* result, const Fp12* bases,
//		const unsigned char* exponents, size_t exp_bits, size_t num)
//	{
//#ifdef BLST_GT_KERNELS
//		return multi_exp_gt_cu(result, bases, exponents, exp_bits, num);
//#else
//		return multi_exp_gt_cc(result, bases, exponents, exp_bits, num);
//#endif
//	}

	// Vectorized finalverify operation (CPU)
	batched_error_code finalverify_gt_vec_cc(bool* results, const Fp12* gt1, const Fp12* gt2, size_t num)
	{
		if (!results || !gt1 || !gt2) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			*results = blst_fp12_finalverify(*gt1, *gt2);
			results++; gt1++; gt2++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for finalverify operation (GPU)
	__global__ void finalverify_gt_vec_ker(bool* results, const Fp12* gt1, const Fp12* gt2, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			results[idx] = blst_fp12_finalverify(gt1[idx], gt2[idx]);
		}
	}

	// Vectorized finalverify operation (GPU)
	batched_error_code finalverify_gt_vec_cu(bool* results, const Fp12* gt1, const Fp12* gt2, size_t num)
	{
		if (!results || !gt1 || !gt2) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(finalverify_gt_vec_ker, num, results, gt1, gt2);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized finalverify operation
	batched_error_code finalverify_gt_vec(bool* results, const Fp12* gt1, const Fp12* gt2, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return finalverify_gt_vec_cu(results, gt1, gt2, num);
#else
		return finalverify_gt_vec_cc(results, gt1, gt2, num);
#endif
	}

	// Batch pairing verification (CPU)
	batched_error_code pairing_verify_gt_vec_cc(bool* results,
		const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points,
		const Fp12* target_values,
		size_t num)
	{
		if (!results || !g2_points || !g1_points || !target_values) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			Fp12 computed;
			blst_miller_loop(computed, g2_points, g1_points);
			Fp12::final_exp(computed, computed);
			*results = Fp12::equal(computed, *target_values);
			results++; g2_points++; g1_points++; target_values++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for batch pairing verification (GPU)
	__global__ void pairing_verify_gt_vec_ker(bool* results,
		const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points,
		const Fp12* target_values,
		size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp12 computed;
			blst_miller_loop(computed, &g2_points[idx], &g1_points[idx]);
			Fp12::final_exp(computed, computed);
			results[idx] = Fp12::equal(computed, target_values[idx]);
		}
	}

	// Batch pairing verification (GPU)
	batched_error_code pairing_verify_gt_vec_cu(bool* results,
		const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points,
		const Fp12* target_values,
		size_t num)
	{
		if (!results || !g2_points || !g1_points || !target_values) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(pairing_verify_gt_vec_ker, num, results, g2_points, g1_points, target_values);
		return cuda_to_batched(future.get());
	}
#endif

	// Batch pairing verification
	batched_error_code pairing_verify_gt_vec(bool* results,
		const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points,
		const Fp12* target_values,
		size_t num)
	{
#ifdef BLST_GT_KERNELS
		return pairing_verify_gt_vec_cu(results, g2_points, g1_points, target_values, num);
#else
		return pairing_verify_gt_vec_cc(results, g2_points, g1_points, target_values, num);
#endif
	}

}