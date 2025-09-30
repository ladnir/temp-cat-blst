/*
 * Copyright Category
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_G2_CUDA_H
#define BLST_G2_CUDA_H

#include "blst/curve/point.h"
#include "blst/curve/e2.h"
#include "blst/common/errors.h"

namespace blst
{


    ///////////////////////////////////////////
    // Vectorized G2 operations

    // Vectorized addition of G2 points. All inputs and outputs are in Jacobian coordinates.
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code add_g2_vec(G2Point* out, const G2Point* a, const G2Point* b, size_t num);

    // Vectorized doubling of G2 points
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code double_g2_vec(G2Point* out, const G2Point* a, size_t num);

    // Vectorized conditional negation of G2 points
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code cneg_g2_vec(G2Point* out, const G2Point* a, int flag, size_t num);

    // Vectorized scalar multiplication of G2 points
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code mult_g2_vec(G2Point* out, const G2Point* points, const unsigned char* scalars, size_t scalar_bits, size_t num);

    ///////////////////////////////////////////
    // G2 coordinate conversion operations

    // Vectorized Jacobian to affine conversion
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code to_affine_g2_vec(G2AffinePoint* out, const G2Point* a, size_t num);

    // Vectorized affine to Jacobian conversion
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code from_affine_g2_vec(G2Point* out, const G2AffinePoint* a, size_t num);

    // Vectorized Jacobian normalization
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code from_jacobian_g2_vec(G2Point* out, const G2Point* a, size_t num);

    ///////////////////////////////////////////
    // G2 validation operations

    // Vectorized curve membership check
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code on_curve_g2_vec(bool* results, const G2Point* points, size_t num);

    // Vectorized infinity check
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code is_inf_g2_vec(bool* results, const G2Point* points, size_t num);

    // Vectorized point equality check
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code is_equal_g2_vec(bool* results, const G2Point* a, const G2Point* b, size_t num);

    // Vectorized subgroup membership check (G2-specific)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code in_group_g2_vec(bool* results, const G2AffinePoint* points, size_t num);

    ///////////////////////////////////////////
    // G2 serialization operations

    // Vectorized point serialization (uncompressed, 192 bytes per point)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code serialize_g2_vec(unsigned char* out, const G2Point* points, size_t num);

    // Vectorized point deserialization (uncompressed, 192 bytes per point)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code deserialize_g2_vec(G2AffinePoint* out, const unsigned char* serialized, BLST_ERROR* errors, size_t num);

    // Vectorized point compression (96 bytes per point)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code compress_g2_vec(unsigned char* out, const G2Point* points, size_t num);

    // Vectorized point decompression (96 bytes per point)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code uncompress_g2_vec(G2AffinePoint* out, const unsigned char* compressed, BLST_ERROR* errors, size_t num);

    ///////////////////////////////////////////
    // Advanced G2 operations

    // Vectorized hash-to-curve for G2
    // inputs must be allocated on the device if compiled with BLST_CUDA
    // msgs is an array of message pointers, msg_lens contains the length of each message
    batched_error_code hash_to_g2_vec(G2Point* out, const unsigned char* const* msgs, const size_t* msg_lens,
        const unsigned char* dst, size_t dst_len, size_t num);

    // Multi-scalar multiplication (MSM) for G2
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code multi_scalar_mult_g2(G2Point* result, const G2Point* points,
        const unsigned char* scalars, size_t scalar_bits, size_t num);

    // Vectorized mixed addition (Jacobian + affine -> Jacobian)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code add_mixed_g2_vec(G2Point* out, const G2Point* a, const G2AffinePoint* b, size_t num);

    // Vectorized psi endomorphism (G2-specific Frobenius map)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code psi_g2_vec(G2Point* out, const G2Point* a, size_t num);

}

#endif /* BLST_G2_CUDA_H */