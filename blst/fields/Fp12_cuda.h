/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_FP12_CUDA_H
#define BLST_FP12_CUDA_H

#include "blst/fields/Fp12.h"
#include "blst/common/errors.h"

BLST_EXTERN_C_BEGIN

///////////////////////////////////////////
// Vectorized Fp12 operations

// Vectorized multiplication of Fp12 elements. All inputs and outputs are in the Montgomery domain.
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code mul_fp12_vec(vec384fp12* out, const vec384fp12* a, const vec384fp12* b, size_t num);

// Vectorized squaring of Fp12 elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code sqr_fp12_vec(vec384fp12* out, const vec384fp12* a, size_t num);

// Vectorized cyclotomic squaring of Fp12 elements (for elements in cyclotomic subgroup)
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code cyclotomic_sqr_fp12_vec(vec384fp12* out, const vec384fp12* a, size_t num);

// Vectorized conjugation of Fp12 elements (negates c1 component)
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code conjugate_fp12_vec(vec384fp12* out, const vec384fp12* a, size_t num);

// Vectorized inversion of Fp12 elements
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code inverse_fp12_vec(vec384fp12* out, const vec384fp12* a, size_t num);

// Vectorized Frobenius map of Fp12 elements
// inputs must be allocated on the device if compiled with BLST_CUDA
// power: the power of the Frobenius endomorphism (1, 2, or 3)
batched_error_code frobenius_map_fp12_vec(vec384fp12* out, const vec384fp12* a, size_t power, size_t num);

// Vectorized equality check for Fp12 elements
// inputs must be allocated on the device if compiled with BLST_CUDA
// out: array of bool results (true if a[i] == b[i])
batched_error_code is_equal_fp12_vec(bool* out, const vec384fp12* a, const vec384fp12* b, size_t num);

// Vectorized check if Fp12 elements are one
// inputs must be allocated on the device if compiled with BLST_CUDA
// out: array of bool results (true if a[i] == 1)
batched_error_code is_one_fp12_vec(bool* out, const vec384fp12* a, size_t num);

// Vectorized serialization of Fp12 elements to big-endian byte representation
// inputs must be allocated on the device if compiled with BLST_CUDA
// out: array of byte arrays, each 12*48 = 576 bytes
batched_error_code bendian_from_fp12_vec(unsigned char* out, const vec384fp12* a, size_t num);

// Vectorized deserialization of Fp12 elements from big-endian byte representation
// inputs must be allocated on the device if compiled with BLST_CUDA
// in: array of byte arrays, each 12*48 = 576 bytes
batched_error_code fp12_from_bendian_vec(vec384fp12* out, const unsigned char* in, size_t num);

// Vectorized multiplication by xy00z0 form (optimized for pairing operations)
// inputs must be allocated on the device if compiled with BLST_CUDA
// xy00z0: array of vec384fp6 elements representing the xy00z0 form
batched_error_code mul_by_xy00z0_fp12_vec(vec384fp12* out, const vec384fp12* a, const vec384fp6* xy00z0, size_t num);

BLST_EXTERN_C_END
#endif /* BLST_FP12_CUDA_H */