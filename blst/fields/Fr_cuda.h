/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_FR_CUDA_H
#define BLST_FR_CUDA_H

#include "blst/fields/Fr.h"
#include "blst/common/errors.h"
//BLST_EXTERN_C_BEGIN

namespace blst
{

	///////////////////////////////////////////
	// Vectorized Fp operations


	// Vectorized addition of Fr elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code add_fr_vec(Fr* out, const Fr* a, const Fr* b, size_t num);

	// Vectorized subtraction of Fr elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code sub_fr_vec(Fr* out, const Fr* a, const Fr* b, size_t num);

	// Vectorized multiplication of Fr elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code mul_fr_vec(Fr* out, const Fr* a, const Fr* b, size_t num);

	// Vectorized squaring of Fr elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code sqr_fr_vec(Fr* out, const Fr* a, size_t num);

	// Vectorized conditional negation of Fr elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code cneg_fr_vec(Fr* out, const Fr* a, int flag, size_t num);

	// Vectorized multiplication of Fr elements by 3
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code mul3_fr_vec(Fr* out, const Fr* a, size_t num);

	// Vectorized left shift of Fr elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code lshift_fr_vec(Fr* out, const Fr* a, size_t count, size_t num);

	// Vectorized right shift of Fr elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code rshift_fr_vec(Fr* out, const Fr* a, size_t count, size_t num);

	// Vectorized conversion from regular to Montgomery form for Fr elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code to_fr_vec(Fr* out, const Fr* a, size_t num);

	// Vectorized conversion from Montgomery to regular form for Fr elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code from_fr_vec(Fr* out, const Fr* a, size_t num);

	// Vectorized CT butterfly operations for Fr elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code ct_bfly_fr_vec(Fr* x0, Fr* x1, const Fr* twiddle, size_t num);

	// Vectorized GS butterfly operations for Fr elements
	// inputs must be allocated on the device if compiled with BLST_CUDA
	batched_error_code gs_bfly_fr_vec(Fr* x0, Fr* x1, const Fr* twiddle, size_t num);


}
//BLST_EXTERN_C_END
#endif /* BLST_FIELDS_CUDA_H */
