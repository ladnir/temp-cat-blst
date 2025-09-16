/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_G1_CUDA_H
#define BLST_G1_CUDA_H

#include "blst/curve/point.h"
#include "blst/curve/e1.h"
#include "blst/common/errors.h"

BLST_EXTERN_C_BEGIN

///////////////////////////////////////////
// Vectorized G1 operations

// Vectorized addition of G1 points. All inputs and outputs are in Jacobian coordinates.
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code add_g1_vec(POINTonE1* out, const POINTonE1* a, const POINTonE1* b, size_t num);

// Vectorized doubling of G1 points
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code double_g1_vec(POINTonE1* out, const POINTonE1* a, size_t num);

// Vectorized conditional negation of G1 points
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code cneg_g1_vec(POINTonE1* out, const POINTonE1* a, int flag, size_t num);

// Vectorized scalar multiplication of G1 points
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code mult_g1_vec(POINTonE1* out, const POINTonE1* points, const unsigned char* scalars, size_t scalar_bits, size_t num);

///////////////////////////////////////////
// G1 coordinate conversion operations

// Vectorized Jacobian to affine conversion
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code to_affine_g1_vec(POINTonE1_affine* out, const POINTonE1* a, size_t num);

// Vectorized affine to Jacobian conversion
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code from_affine_g1_vec(POINTonE1* out, const POINTonE1_affine* a, size_t num);

// Vectorized Jacobian normalization
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code from_jacobian_g1_vec(POINTonE1* out, const POINTonE1* a, size_t num);

///////////////////////////////////////////
// G1 validation operations

// Vectorized curve membership check
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code on_curve_g1_vec(bool* results, const POINTonE1* points, size_t num);

// Vectorized infinity check
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code is_inf_g1_vec(bool* results, const POINTonE1* points, size_t num);

// Vectorized point equality check
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code is_equal_g1_vec(bool* results, const POINTonE1* a, const POINTonE1* b, size_t num);

///////////////////////////////////////////
// G1 serialization operations

// Vectorized point serialization (uncompressed)
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code serialize_g1_vec(unsigned char* out, const POINTonE1* points, size_t num);

// Vectorized point deserialization (uncompressed)
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code deserialize_g1_vec(POINTonE1_affine* out, const unsigned char* serialized, size_t num);

// Vectorized point compression
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code compress_g1_vec(unsigned char* out, const POINTonE1* points, size_t num);

// Vectorized point decompression
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code uncompress_g1_vec(POINTonE1_affine* out, const unsigned char* compressed, size_t num);

///////////////////////////////////////////
// Advanced G1 operations

// Vectorized hash-to-curve
// inputs must be allocated on the device if compiled with BLST_CUDA
// msgs is a concatenated array of messages, msg_lens contains the length of each message
batched_error_code hash_to_g1_vec(POINTonE1* out, const unsigned char* msgs, const size_t* msg_lens,
    const unsigned char* dst, size_t dst_len, size_t num);

// Multi-scalar multiplication (MSM)
// inputs must be allocated on the device if compiled with BLST_CUDA
batched_error_code multi_scalar_mult_g1(POINTonE1* result, const POINTonE1* points,
    const unsigned char* scalars, size_t scalar_bits, size_t num);

BLST_EXTERN_C_END
#endif /* BLST_G1_CUDA_H */