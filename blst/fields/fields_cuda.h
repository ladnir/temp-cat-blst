/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_FIELDS_CUDA_H
#define BLST_FIELDS_CUDA_H

#include "blst/fields/fields.h"
BLST_EXTERN_C_BEGIN

// Error codes for CUDA field operations - focused on most common errors
typedef enum {
	BATCHED_SUCCESS = 0,
	BATCHED_ERROR_NULL_POINTER,
	BATCHED_ERROR_INVALID_SIZE,
	BATCHED_ERROR_ILLEGAL_ADDRESS,
	BATCHED_ERROR_CUDA_MALLOC_FAILED,
	BATCHED_ERROR_CUDA_MEMCPY_FAILED,
	BATCHED_ERROR_CUDA_KERNEL_LAUNCH_FAILED,
	BATCHED_ERROR_CUDA_SYNC_FAILED,
	BATCHED_ERROR_CUDA_DEVICE_ERROR,
	BATCHED_ERROR_HOST_MALLOC_FAILED,
	BATCHED_ERROR_UNKNOWN
} batched_error_code;

// Convert error code to string
BLST_INLINE const char* batched_error_code_string(batched_error_code error)
{
	switch (error) {
	case BATCHED_SUCCESS:
		return "Success";
	case BATCHED_ERROR_NULL_POINTER:
		return "Invalid pointer arguments";
	case BATCHED_ERROR_INVALID_SIZE:
		return "Invalid number of elements";
	case BATCHED_ERROR_ILLEGAL_ADDRESS:
		return "Illegal memory access";
	case BATCHED_ERROR_CUDA_MALLOC_FAILED:
		return "CUDA memory allocation failed";
	case BATCHED_ERROR_CUDA_MEMCPY_FAILED:
		return "CUDA memory copy failed";
	case BATCHED_ERROR_CUDA_KERNEL_LAUNCH_FAILED:
		return "CUDA kernel launch failed";
	case BATCHED_ERROR_CUDA_SYNC_FAILED:
		return "CUDA device synchronization failed";
	case BATCHED_ERROR_CUDA_DEVICE_ERROR:
		return "CUDA device error";
	case BATCHED_ERROR_HOST_MALLOC_FAILED:
		return "Host memory allocation failed";
	case BATCHED_ERROR_UNKNOWN:
	default:
		return "Unknown error";
	}
}

///////////////////////////////////////////
// Vectorized Fp operations

// Vectorized addition of Fp elements. All inputs and outputs are in the Montgomery domain.
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code add_fp_vec(vec384* out, const vec384* a, const vec384* b, size_t num);

// Vectorized subtraction of Fp elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code sub_fp_vec(vec384* out, const vec384* a, const vec384* b, size_t num);

// Vectorized multiplication of Fp elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code mul_fp_vec(vec384* out, const vec384* a, const vec384* b, size_t num);

// Vectorized squaring of Fp elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code sqr_fp_vec(vec384* out, const vec384* a, size_t num);

// Vectorized negation of Fp elements
batched_error_code cneg_fp_vec(vec384* out, const vec384* a, int flag, size_t num);

// Vectorized multiplication of Fp elements by 3
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code mul3_fp_vec(vec384* out, const vec384* a, size_t num);

// Vectorized multiplication of Fp elements by 8
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code mul8_fp_vec(vec384* out, const vec384* a, size_t num);

// Vectorized division of Fp elements by 2
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code div2_fp_vec(vec384* out, const vec384* a
	, size_t num);


// Vectorized left shift of Fp elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code lshift_fp_vec(vec384* out, const vec384* a
	, size_t count, size_t num);


// Vectorized conversion from a regular representation to a montgomery Fp
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code to_fp_vec(vec384* out, const vec384* a
	, size_t num);


// Vectorized conversion from a montgomery Fp to a regular representation
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code from_fp_vec(vec384* out, const vec384* a,
	size_t num);


///////////////////////////////////////////
// Vectorized Fp operations


// Vectorized addition of Fr elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code add_fr_vec(vec256* out, const vec256* a, const vec256* b, size_t num);

// Vectorized subtraction of Fr elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code sub_fr_vec(vec256* out, const vec256* a, const vec256* b, size_t num);

// Vectorized multiplication of Fr elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code mul_fr_vec(vec256* out, const vec256* a, const vec256* b, size_t num);

// Vectorized squaring of Fr elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code sqr_fr_vec(vec256* out, const vec256* a, size_t num);

// Vectorized conditional negation of Fr elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code cneg_fr_vec(vec256* out, const vec256* a, int flag, size_t num);

// Vectorized multiplication of Fr elements by 3
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code mul3_fr_vec(vec256* out, const vec256* a, size_t num);

// Vectorized left shift of Fr elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code lshift_fr_vec(vec256* out, const vec256* a, size_t count, size_t num);

// Vectorized right shift of Fr elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code rshift_fr_vec(vec256* out, const vec256* a, size_t count, size_t num);

// Vectorized conversion from regular to Montgomery form for Fr elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code to_fr_vec(vec256* out, const vec256* a, size_t num);

// Vectorized conversion from Montgomery to regular form for Fr elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code from_fr_vec(vec256* out, const vec256* a, size_t num);

// Vectorized CT butterfly operations for Fr elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code ct_bfly_fr_vec(vec256* x0, vec256* x1, const vec256* twiddle, size_t num);

// Vectorized GS butterfly operations for Fr elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code gs_bfly_fr_vec(vec256* x0, vec256* x1, const vec256* twiddle, size_t num);



///////////////////////////////////////////
// Vectorized Fp2 operations


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


BLST_EXTERN_C_END
#endif /* BLST_FIELDS_CUDA_H */
