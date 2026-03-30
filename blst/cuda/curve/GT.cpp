/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/curve/point.h"
#include "blst/curve/GT.h"
#include "pairing.h"
#include "blst/curve/GT_cuda.h"
#include "blst/exports.h"
#include "blst/common/errors.h"
#include "blst/sign/aggregate.h"

#define BLST_GT_KERNELS

#ifdef BLST_GT_KERNELS
#include "blst/common/cuda_launch.h"
#endif

namespace blst
{

	///////////////////////////////////////////
	// Core GT Operations

	// Vectorized multiplication of GT elements (CPU)
	batched_error_code mult_gt_vec_cc(GT* out, const GT* a, const GT* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			GT::add(*out, *a, *b);
			out++; a++; b++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for multiplication of GT elements (GPU)
	__global__ void mult_gt_vec_ker(GT* out, const GT* a, const GT* b, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			GT::add(out[idx], a[idx], b[idx]);
		}
	}

	// Vectorized multiplication of GT elements (GPU)
	batched_error_code mult_gt_vec_cu(GT* out, const GT* a, const GT* b, size_t num)
	{
		if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(mult_gt_vec_ker, num, out, a, b);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized multiplication of GT elements
	batched_error_code mult_gt_vec(GT* out, const GT* a, const GT* b, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return mult_gt_vec_cu(out, a, b, num);
#else
		return mult_gt_vec_cc(out, a, b, num);
#endif
	}

	// Vectorized squaring of GT elements (CPU)
	batched_error_code square_gt_vec_cc(GT* out, const GT* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			GT::doubling(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for squaring of GT elements (GPU)
	__global__ void square_gt_vec_ker(GT* out, const GT* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			GT::doubling(out[idx], a[idx]);
		}
	}

	// Vectorized squaring of GT elements (GPU)
	batched_error_code square_gt_vec_cu(GT* out, const GT* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(square_gt_vec_ker, num, out, a);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized squaring of GT elements
	batched_error_code square_gt_vec(GT* out, const GT* a, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return square_gt_vec_cu(out, a, num);
#else
		return square_gt_vec_cc(out, a, num);
#endif
	}

	// Vectorized cyclotomic squaring of GT elements (CPU)
	batched_error_code cyclotomic_square_gt_vec_cc(GT* out, const GT* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			GT::doubling(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for cyclotomic squaring of GT elements (GPU)
	__global__ void cyclotomic_square_gt_vec_ker(GT* out, const GT* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			GT::doubling(out[idx], a[idx]);
		}
	}

	// Vectorized cyclotomic squaring of GT elements (GPU)
	batched_error_code cyclotomic_square_gt_vec_cu(GT* out, const GT* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(cyclotomic_square_gt_vec_ker, num, out, a);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized cyclotomic squaring of GT elements
	batched_error_code cyclotomic_square_gt_vec(GT* out, const GT* a, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return cyclotomic_square_gt_vec_cu(out, a, num);
#else
		return cyclotomic_square_gt_vec_cc(out, a, num);
#endif
	}

	// Vectorized inversion of GT elements (CPU)
	batched_error_code inverse_gt_vec_cc(GT* out, const GT* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			GT::neg(*out, *a);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for inversion of GT elements (GPU)
	__global__ void inverse_gt_vec_ker(GT* out, const GT* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			GT::neg(out[idx], a[idx]);
		}
	}

	// Vectorized inversion of GT elements (GPU)
	batched_error_code inverse_gt_vec_cu(GT* out, const GT* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(inverse_gt_vec_ker, num, out, a);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized inversion of GT elements
	batched_error_code inverse_gt_vec(GT* out, const GT* a, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return inverse_gt_vec_cu(out, a, num);
#else
		return inverse_gt_vec_cc(out, a, num);
#endif
	}

	// Vectorized conjugation of GT elements (CPU)
	batched_error_code conjugate_gt_vec_cc(GT* out, const GT* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			vec_copy(out, a, sizeof(GT));
			GT::conjugate(*out);
			out++; a++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for conjugation of GT elements (GPU)
	__global__ void conjugate_gt_vec_ker(GT* out, const GT* a, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			vec_copy(&out[idx], &a[idx], sizeof(GT));
			GT::conjugate(out[idx]);
		}
	}

	// Vectorized conjugation of GT elements (GPU)
	batched_error_code conjugate_gt_vec_cu(GT* out, const GT* a, size_t num)
	{
		if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(conjugate_gt_vec_ker, num, out, a);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized conjugation of GT elements
	batched_error_code conjugate_gt_vec(GT* out, const GT* a, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return conjugate_gt_vec_cu(out, a, num);
#else
		return conjugate_gt_vec_cc(out, a, num);
#endif
	}

	// Vectorized exponentiation of GT elements (CPU)
	batched_error_code pow_gt_vec_cc(GT* out, const GT* bases, const unsigned char* exponents,
		size_t exp_bits, size_t num)
	{
		if (!out || !bases || !exponents) return BATCHED_ERROR_NULL_POINTER;

		const size_t exp_bytes = (exp_bits + 7) / 8;

		while (num--) {
			GT::mult(*out, *bases, exponents, exp_bits);
			out++; bases++; exponents += exp_bytes;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for exponentiation of GT elements (GPU)
	__global__ void pow_gt_vec_ker(GT* out, const GT* bases, const unsigned char* exponents,
		size_t exp_bits, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			const size_t exp_bytes = (exp_bits + 7) / 8;
			GT::mult(out[idx], bases[idx], &exponents[idx * exp_bytes], exp_bits);
		}
	}

	// Vectorized exponentiation of GT elements (GPU)
	batched_error_code pow_gt_vec_cu(GT* out, const GT* bases, const unsigned char* exponents,
		size_t exp_bits, size_t num)
	{
		if (!out || !bases || !exponents) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(pow_gt_vec_ker, num, out, bases, exponents, exp_bits);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized exponentiation of GT elements
	batched_error_code pow_gt_vec(GT* out, const GT* bases, const unsigned char* exponents,
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
	batched_error_code is_one_gt_vec_cc(bool* results, const GT* elements, size_t num)
	{
		if (!results || !elements) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			*results = GT::is_identity(*elements);
			results++; elements++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for identity check of GT elements (GPU)
	__global__ void is_one_gt_vec_ker(bool* results, const GT* elements, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			results[idx] = GT::is_identity(elements[idx]);
		}
	}

	// Vectorized check if GT elements equal identity (GPU)
	batched_error_code is_one_gt_vec_cu(bool* results, const GT* elements, size_t num)
	{
		if (!results || !elements) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(is_one_gt_vec_ker, num, results, elements);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized check if GT elements equal identity
	batched_error_code is_one_gt_vec(bool* results, const GT* elements, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return is_one_gt_vec_cu(results, elements, num);
#else
		return is_one_gt_vec_cc(results, elements, num);
#endif
	}

	// Vectorized GT element equality check (CPU)
	batched_error_code is_equal_gt_vec_cc(bool* results, const GT* a, const GT* b, size_t num)
	{
		if (!results || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			*results = GT::equal(*a, *b);
			results++; a++; b++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for GT element equality check (GPU)
	__global__ void is_equal_gt_vec_ker(bool* results, const GT* a, const GT* b, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			results[idx] = GT::equal(a[idx], b[idx]);
		}
	}

	// Vectorized GT element equality check (GPU)
	batched_error_code is_equal_gt_vec_cu(bool* results, const GT* a, const GT* b, size_t num)
	{
		if (!results || !a || !b) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(is_equal_gt_vec_ker, num, results, a, b);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized GT element equality check
	batched_error_code is_equal_gt_vec(bool* results, const GT* a, const GT* b, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return is_equal_gt_vec_cu(results, a, b, num);
#else
		return is_equal_gt_vec_cc(results, a, b, num);
#endif
	}

	// Vectorized GT group membership check (CPU)
	batched_error_code in_group_gt_vec_cc(bool* results, const GT* elements, size_t num)
	{
		if (!results || !elements) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			*results = GT::in_group(*elements);
			results++; elements++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for GT group membership check (GPU)
	__global__ void in_group_gt_vec_ker(bool* results, const GT* elements, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			results[idx] = GT::in_group(elements[idx]);
		}
	}

	// Vectorized GT group membership check (GPU)
	batched_error_code in_group_gt_vec_cu(bool* results, const GT* elements, size_t num)
	{
		if (!results || !elements) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(in_group_gt_vec_ker, num, results, elements);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized GT group membership check
	batched_error_code in_group_gt_vec(bool* results, const GT* elements, size_t num)
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
	batched_error_code miller_loop_gt_vec_cc(GT* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		if (!out || !g2_points || !g1_points) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			GT::miller_loop(*out, g2_points, g1_points);
			out++; g2_points++; g1_points++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for Miller loop computation (GPU)
	__global__ void miller_loop_gt_vec_ker(GT* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			GT::miller_loop(out[idx], &g2_points[idx], &g1_points[idx]);
		}
	}

	// Vectorized Miller loop computation (GPU)
	batched_error_code miller_loop_gt_vec_cu(GT* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		if (!out || !g2_points || !g1_points) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(miller_loop_gt_vec_ker, num, out, g2_points, g1_points);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized Miller loop computation
	batched_error_code miller_loop_gt_vec(GT* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return miller_loop_gt_vec_cu(out, g2_points, g1_points, num);
#else
		return miller_loop_gt_vec_cc(out, g2_points, g1_points, num);
#endif
	}

	// Vectorized final exponentiation (CPU)
	batched_error_code final_exp_gt_vec_cc(GT* out, const GT* miller_results, size_t num)
	{
		if (!out || !miller_results) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			GT::final_mult(*out, *miller_results);
			out++; miller_results++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for final exponentiation (GPU)
	__global__ void final_exp_gt_vec_ker(GT* out, const GT* miller_results, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			GT::final_mult(out[idx], miller_results[idx]);
		}
	}

	// Vectorized final exponentiation (GPU)
	batched_error_code final_exp_gt_vec_cu(GT* out, const GT* miller_results, size_t num)
	{
		if (!out || !miller_results) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(final_exp_gt_vec_ker, num, out, miller_results);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized final exponentiation
	batched_error_code final_exp_gt_vec(GT* out, const GT* miller_results, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return final_exp_gt_vec_cu(out, miller_results, num);
#else
		return final_exp_gt_vec_cc(out, miller_results, num);
#endif
	}

	// Vectorized complete pairing computation (CPU)
	batched_error_code pairing_gt_vec_cc(GT* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		if (!out || !g2_points || !g1_points) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			GT::miller_loop(*out, g2_points, g1_points);
			GT::final_mult(*out, *out);
			out++; g2_points++; g1_points++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for complete pairing computation (GPU)
	__global__ void pairing_gt_vec_ker(GT* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			GT::miller_loop(out[idx], &g2_points[idx], &g1_points[idx]);
			GT::final_mult(out[idx], out[idx]);
		}
	}

	// Vectorized complete pairing computation (GPU)
	batched_error_code pairing_gt_vec_cu(GT* out, const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points, size_t num)
	{
		if (!out || !g2_points || !g1_points) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(pairing_gt_vec_ker, num, out, g2_points, g1_points);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized complete pairing computation
	batched_error_code pairing_gt_vec(GT* out, const G2AffinePoint* g2_points,
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
	batched_error_code serialize_gt_vec_cc(unsigned char* out, const GT* elements, size_t num)
	{
		if (!out || !elements) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			Fp12::to_bendian(out, elements->data);
			out += 576; elements++; // 12 * 48 bytes per GT element
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for GT element serialization (GPU)
	__global__ void serialize_gt_vec_ker(unsigned char* out, const GT* elements, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp12::to_bendian(&out[idx * 576], elements[idx].data);
		}
	}

	// Vectorized GT element serialization (GPU)
	batched_error_code serialize_gt_vec_cu(unsigned char* out, const GT* elements, size_t num)
	{
		if (!out || !elements) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(serialize_gt_vec_ker, num, out, elements);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized GT element serialization
	batched_error_code serialize_gt_vec(unsigned char* out, const GT* elements, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return serialize_gt_vec_cu(out, elements, num);
#else
		return serialize_gt_vec_cc(out, elements, num);
#endif
	}

	// Vectorized GT element deserialization (CPU)
	batched_error_code deserialize_gt_vec_cc(GT* out, const unsigned char* serialized, size_t num)
	{
		if (!out || !serialized) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			Fp12::from_bendian(out->data, serialized);
			out++; serialized += 576;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for GT element deserialization (GPU)
	__global__ void deserialize_gt_vec_ker(GT* out, const unsigned char* serialized, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			Fp12::from_bendian(out[idx].data, &serialized[idx * 576]);
		}
	}

	// Vectorized GT element deserialization (GPU)
	batched_error_code deserialize_gt_vec_cu(GT* out, const unsigned char* serialized, size_t num)
	{
		if (!out || !serialized) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(deserialize_gt_vec_ker, num, out, serialized);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized GT element deserialization
	batched_error_code deserialize_gt_vec(GT* out, const unsigned char* serialized, size_t num)
	{
#ifdef BLST_GT_KERNELS
		return deserialize_gt_vec_cu(out, serialized, num);
#else
		return deserialize_gt_vec_cc(out, serialized, num);
#endif
	}


	// Vectorized finalverify operation (CPU)
	batched_error_code finalverify_gt_vec_cc(bool* results, const GT* gt1, const GT* gt2, size_t num)
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
	__global__ void finalverify_gt_vec_ker(bool* results, const GT* gt1, const GT* gt2, size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			results[idx] = blst_fp12_finalverify(gt1[idx], gt2[idx]);
		}
	}

	// Vectorized finalverify operation (GPU)
	batched_error_code finalverify_gt_vec_cu(bool* results, const GT* gt1, const GT* gt2, size_t num)
	{
		if (!results || !gt1 || !gt2) return BATCHED_ERROR_NULL_POINTER;

		auto future = launch_auto(finalverify_gt_vec_ker, num, results, gt1, gt2);
		return cuda_to_batched(future.get());
	}
#endif

	// Vectorized finalverify operation
	batched_error_code finalverify_gt_vec(bool* results, const GT* gt1, const GT* gt2, size_t num)
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
		const GT* target_values,
		size_t num)
	{
		if (!results || !g2_points || !g1_points || !target_values) return BATCHED_ERROR_NULL_POINTER;

		while (num--) {
			GT computed;
			GT::miller_loop(computed, g2_points, g1_points);
			GT::final_mult(computed, computed);
			*results = GT::equal(computed, *target_values);
			results++; g2_points++; g1_points++; target_values++;
		}
		return BATCHED_SUCCESS;
	}

#ifdef BLST_GT_KERNELS
	// Kernel for batch pairing verification (GPU)
	__global__ void pairing_verify_gt_vec_ker(bool* results,
		const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points,
		const GT* target_values,
		size_t num)
	{
		size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < num) {
			GT computed;
			GT::miller_loop(computed, &g2_points[idx], &g1_points[idx]);
			GT::final_mult(computed, computed);
			results[idx] = GT::equal(computed, target_values[idx]);
		}
	}

	// Batch pairing verification (GPU)
	batched_error_code pairing_verify_gt_vec_cu(bool* results,
		const G2AffinePoint* g2_points,
		const G1AffinePoint* g1_points,
		const GT* target_values,
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
		const GT* target_values,
		size_t num)
	{
#ifdef BLST_GT_KERNELS
		return pairing_verify_gt_vec_cu(results, g2_points, g1_points, target_values, num);
#else
		return pairing_verify_gt_vec_cc(results, g2_points, g1_points, target_values, num);
#endif
	}

}