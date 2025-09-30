/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_GT_CUDA_H
#define BLST_GT_CUDA_H

#include "blst/curve/point.h"
#include "blst/curve/e2.h"
#include "blst/fields/Fp12.h"
#include "blst/common/errors.h"

namespace blst
{

    ///////////////////////////////////////////
    // Vectorized GT operations (Fp12 group operations)

    // Vectorized multiplication of GT elements
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code mult_gt_vec(Fp12* out, const Fp12* a, const Fp12* b, size_t num);

    // Vectorized squaring of GT elements
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code square_gt_vec(Fp12* out, const Fp12* a, size_t num);

    // Vectorized cyclotomic squaring of GT elements (optimized for GT group)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code cyclotomic_square_gt_vec(Fp12* out, const Fp12* a, size_t num);

    // Vectorized inversion of GT elements
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code inverse_gt_vec(Fp12* out, const Fp12* a, size_t num);

    // Vectorized conjugation of GT elements
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code conjugate_gt_vec(Fp12* out, const Fp12* a, size_t num);

    // Vectorized exponentiation of GT elements
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code pow_gt_vec(Fp12* out, const Fp12* bases, const unsigned char* exponents,
        size_t exp_bits, size_t num);

    ///////////////////////////////////////////
    // GT validation operations

    // Vectorized check if GT elements equal identity (one)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code is_one_gt_vec(bool* results, const Fp12* elements, size_t num);

    // Vectorized GT element equality check
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code is_equal_gt_vec(bool* results, const Fp12* a, const Fp12* b, size_t num);

    // Vectorized GT group membership check
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code in_group_gt_vec(bool* results, const Fp12* elements, size_t num);

    ///////////////////////////////////////////
    // GT pairing-related operations

    // Vectorized Miller loop computation
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code miller_loop_gt_vec(Fp12* out, const G2AffinePoint* g2_points,
        const G1AffinePoint* g1_points, size_t num);

    // Vectorized final exponentiation
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code final_exp_gt_vec(Fp12* out, const Fp12* miller_results, size_t num);

    // Vectorized complete pairing computation (Miller loop + final exp)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code pairing_gt_vec(Fp12* out, const G2AffinePoint* g2_points,
        const G1AffinePoint* g1_points, size_t num);

    ///////////////////////////////////////////
    // GT serialization operations

    // Vectorized GT element serialization (576 bytes per element - 12 * 48)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code serialize_gt_vec(unsigned char* out, const Fp12* elements, size_t num);

    // Vectorized GT element deserialization (576 bytes per element)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code deserialize_gt_vec(Fp12* out, const unsigned char* serialized, size_t num);

    ///////////////////////////////////////////
    // Advanced GT operations

    // Multi-exponentiation for GT elements (MSM for GT group)
    // inputs must be allocated on the device if compiled with BLST_CUDA
    //batched_error_code multi_exp_gt(Fp12* result, const Fp12* bases,
    //    const unsigned char* exponents, size_t exp_bits, size_t num);

    // Vectorized finalverify operation for signature verification
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code finalverify_gt_vec(bool* results, const Fp12* gt1, const Fp12* gt2, size_t num);

    // Batch pairing verification with precomputed lines
    // inputs must be allocated on the device if compiled with BLST_CUDA
    batched_error_code pairing_verify_gt_vec(bool* results,
        const G2AffinePoint* g2_points,
        const G1AffinePoint* g1_points,
        const Fp12* target_values,
        size_t num);

}

#endif /* BLST_GT_CUDA_H */