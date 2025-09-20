/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_FP2_CUDA_H
#define BLST_FP2_CUDA_H

#include "blst/fields/Fp2.h"
#include "blst/common/errors.h"


namespace blst
{

	///////////////////////////////////////////
	// Vectorized Fp2 operations

	// Vectorized addition of Fp2 elements. All inputs and outputs are in the Montgomery domain.
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code add_fp2_vec(vec384x* out, const vec384x* a, const vec384x* b, size_t num);

	// Vectorized subtraction of Fp2 elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code sub_fp2_vec(vec384x* out, const vec384x* a, const vec384x* b, size_t num);

	// Vectorized multiplication of Fp2 elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code mul_fp2_vec(vec384x* out, const vec384x* a, const vec384x* b, size_t num);

	// Vectorized squaring of Fp2 elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code sqr_fp2_vec(vec384x* out, const vec384x* a, size_t num);

	// Vectorized conditional negation of Fp2 elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code cneg_fp2_vec(vec384x* out, const vec384x* a, int flag, size_t num);

	// Vectorized multiplication of Fp2 elements by 3
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code mul3_fp2_vec(vec384x* out, const vec384x* a, size_t num);

	// Vectorized multiplication of Fp2 elements by 8
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code mul8_fp2_vec(vec384x* out, const vec384x* a, size_t num);

	// Vectorized left shift of Fp2 elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code lshift_fp2_vec(vec384x* out, const vec384x* a, size_t count, size_t num);

	// Vectorized conjugation of Fp2 elements (negates imaginary part)
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code conjugate_fp2_vec(vec384x* out, const vec384x* a, size_t num);

	// Vectorized inversion of Fp2 elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code inverse_fp2_vec(vec384x* out, const vec384x* a, size_t num);

	// Vectorized square root of Fp2 elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	// Returns BATCHED_SUCCESS if all elements have square roots, error code otherwise
	batched_error_code sqrt_fp2_vec(vec384x* out, const vec384x* a, size_t num);

}

#endif /* BLST_FIELDS_CUDA_H */