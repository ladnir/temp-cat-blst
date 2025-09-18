/*
 * Copyright Category
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_FP_CUDA_H
#define BLST_FP_CUDA_H

#include "blst/fields/Fp.h"
#include "blst/common/errors.h"

//BLST_EXTERN_C_BEGIN

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


//BLST_EXTERN_C_END
#endif /* BLST_FIELDS_CUDA_H */
