/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/curve/point.h"
#include "blst/curve/e2.h"
#include "blst/curve/G2_cuda.h"
#include "blst/curve/map_to_g2.h"
#include "blst/exports.h"
#include "blst/common/errors.h"

#ifdef BLST_CUDA
#include "blst/common/cuda_launch.h"
#endif

namespace blst
{

    ///////////////////////////////////////////
    // Core G2 Operations

    // Vectorized addition of G2 points (CPU)
    batched_error_code add_g2_vec_cc(G2Point* out, const G2Point* a, const G2Point* b, size_t num)
    {
        if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            blst_p2_add_or_double(*out, *a, *b);
            out++; a++; b++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for addition of G2 points (GPU)
    __global__ void add_g2_vec_ker(G2Point* out, const G2Point* a, const G2Point* b, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            blst_p2_add_or_double(out[idx], a[idx], b[idx]);
        }
    }

    // Vectorized addition of G2 points (GPU)
    batched_error_code add_g2_vec_cu(G2Point* out, const G2Point* a, const G2Point* b, size_t num)
    {
        if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(add_g2_vec_ker, num, out, a, b);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized addition of G2 points
    batched_error_code add_g2_vec(G2Point* out, const G2Point* a, const G2Point* b, size_t num)
    {
#ifdef BLST_CUDA
        return add_g2_vec_cu(out, a, b, num);
#else
        return add_g2_vec_cc(out, a, b, num);
#endif
    }

    // Vectorized doubling of G2 points (CPU)
    batched_error_code double_g2_vec_cc(G2Point* out, const G2Point* a, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            blst_p2_double(*out, *a);
            out++; a++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for doubling of G2 points (GPU)
    __global__ void double_g2_vec_ker(G2Point* out, const G2Point* a, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            blst_p2_double(out[idx], a[idx]);
        }
    }

    // Vectorized doubling of G2 points (GPU)
    batched_error_code double_g2_vec_cu(G2Point* out, const G2Point* a, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(double_g2_vec_ker, num, out, a);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized doubling of G2 points
    batched_error_code double_g2_vec(G2Point* out, const G2Point* a, size_t num)
    {
#ifdef BLST_CUDA
        return double_g2_vec_cu(out, a, num);
#else
        return double_g2_vec_cc(out, a, num);
#endif
    }

    // Vectorized conditional negation of G2 points (CPU)
    batched_error_code cneg_g2_vec_cc(G2Point* out, const G2Point* a, int flag, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            *out = *a;
            blst_p2_cneg(*out, flag);
            out++; a++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for conditional negation of G2 points (GPU)
    __global__ void cneg_g2_vec_ker(G2Point* out, const G2Point* a, int flag, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            out[idx] = a[idx];
            blst_p2_cneg(out[idx], flag);
        }
    }

    // Vectorized conditional negation of G2 points (GPU)
    batched_error_code cneg_g2_vec_cu(G2Point* out, const G2Point* a, int flag, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(cneg_g2_vec_ker, num, out, a, flag);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized conditional negation of G2 points
    batched_error_code cneg_g2_vec(G2Point* out, const G2Point* a, int flag, size_t num)
    {
#ifdef BLST_CUDA
        return cneg_g2_vec_cu(out, a, flag, num);
#else
        return cneg_g2_vec_cc(out, a, flag, num);
#endif
    }

    // Vectorized scalar multiplication of G2 points (CPU)
    batched_error_code mult_g2_vec_cc(G2Point* out, const G2Point* points, const unsigned char* scalars, size_t scalar_bits, size_t num)
    {
        if (!out || !points || !scalars) return BATCHED_ERROR_NULL_POINTER;

        const size_t scalar_bytes = (scalar_bits + 7) / 8;

        while (num--) {
            blst_p2_mult(*out, *points, scalars, scalar_bits);
            out++; points++; scalars += scalar_bytes;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for scalar multiplication of G2 points (GPU)
    __global__ void mult_g2_vec_ker(G2Point* out, const G2Point* points, const unsigned char* scalars, size_t scalar_bits, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            const size_t scalar_bytes = (scalar_bits + 7) / 8;

            blst_p2_mult(
                out[idx],
                points[idx],
                &scalars[idx * scalar_bytes],
                scalar_bits);
        }
    }

    // Vectorized scalar multiplication of G2 points (GPU)
    batched_error_code mult_g2_vec_cu(G2Point* out, const G2Point* points, const unsigned char* scalars, size_t scalar_bits, size_t num)
    {
        if (!out || !points || !scalars) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(mult_g2_vec_ker, num, out, points, scalars, scalar_bits);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized scalar multiplication of G2 points
    batched_error_code mult_g2_vec(G2Point* out, const G2Point* points, const unsigned char* scalars, size_t scalar_bits, size_t num)
    {
#ifdef BLST_CUDA
        return mult_g2_vec_cu(out, points, scalars, scalar_bits, num);
#else
        return mult_g2_vec_cc(out, points, scalars, scalar_bits, num);
#endif
    }

    ///////////////////////////////////////////
    // Coordinate Conversion Operations

    // Vectorized Jacobian to affine conversion (CPU)
    batched_error_code to_affine_g2_vec_cc(G2AffinePoint* out, const G2Point* a, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            blst_p2_to_affine(*out, *a);
            out++; a++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for Jacobian to affine conversion (GPU)
    __global__ void to_affine_g2_vec_ker(G2AffinePoint* out, const G2Point* a, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            blst_p2_to_affine(out[idx], a[idx]);
        }
    }

    // Vectorized Jacobian to affine conversion (GPU)
    batched_error_code to_affine_g2_vec_cu(G2AffinePoint* out, const G2Point* a, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(to_affine_g2_vec_ker, num, out, a);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized Jacobian to affine conversion
    batched_error_code to_affine_g2_vec(G2AffinePoint* out, const G2Point* a, size_t num)
    {
#ifdef BLST_CUDA
        return to_affine_g2_vec_cu(out, a, num);
#else
        return to_affine_g2_vec_cc(out, a, num);
#endif
    }

    // Vectorized affine to Jacobian conversion (CPU)
    batched_error_code from_affine_g2_vec_cc(G2Point* out, const G2AffinePoint* a, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            blst_p2_from_affine(*out, *a);
            out++; a++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for affine to Jacobian conversion (GPU)
    __global__ void from_affine_g2_vec_ker(G2Point* out, const G2AffinePoint* a, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            blst_p2_from_affine(out[idx], a[idx]);
        }
    }

    // Vectorized affine to Jacobian conversion (GPU)
    batched_error_code from_affine_g2_vec_cu(G2Point* out, const G2AffinePoint* a, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(from_affine_g2_vec_ker, num, out, a);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized affine to Jacobian conversion
    batched_error_code from_affine_g2_vec(G2Point* out, const G2AffinePoint* a, size_t num)
    {
#ifdef BLST_CUDA
        return from_affine_g2_vec_cu(out, a, num);
#else
        return from_affine_g2_vec_cc(out, a, num);
#endif
    }

    // Vectorized Jacobian normalization (CPU)
    batched_error_code from_jacobian_g2_vec_cc(G2Point* out, const G2Point* a, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            blst_p2_from_jacobian(*out, *a);
            out++; a++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for Jacobian normalization (GPU)
    __global__ void from_jacobian_g2_vec_ker(G2Point* out, const G2Point* a, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            blst_p2_from_jacobian(out[idx], a[idx]);
        }
    }

    // Vectorized Jacobian normalization (GPU)
    batched_error_code from_jacobian_g2_vec_cu(G2Point* out, const G2Point* a, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(from_jacobian_g2_vec_ker, num, out, a);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized Jacobian normalization
    batched_error_code from_jacobian_g2_vec(G2Point* out, const G2Point* a, size_t num)
    {
#ifdef BLST_CUDA
        return from_jacobian_g2_vec_cu(out, a, num);
#else
        return from_jacobian_g2_vec_cc(out, a, num);
#endif
    }

    ///////////////////////////////////////////
    // Validation Operations

    // Vectorized curve membership check (CPU)
    batched_error_code on_curve_g2_vec_cc(bool* results, const G2Point* points, size_t num)
    {
        if (!results || !points) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            *results = blst_p2_on_curve(*points);
            results++; points++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for curve membership check (GPU)
    __global__ void on_curve_g2_vec_ker(bool* results, const G2Point* points, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            results[idx] = blst_p2_on_curve(points[idx]);
        }
    }

    // Vectorized curve membership check (GPU)
    batched_error_code on_curve_g2_vec_cu(bool* results, const G2Point* points, size_t num)
    {
        if (!results || !points) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(on_curve_g2_vec_ker, num, results, points);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized curve membership check
    batched_error_code on_curve_g2_vec(bool* results, const G2Point* points, size_t num)
    {
#ifdef BLST_CUDA
        return on_curve_g2_vec_cu(results, points, num);
#else
        return on_curve_g2_vec_cc(results, points, num);
#endif
    }

    // Vectorized infinity check (CPU)
    batched_error_code is_inf_g2_vec_cc(bool* results, const G2Point* points, size_t num)
    {
        if (!results || !points) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            *results = blst_p2_is_inf(*points);
            results++; points++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for infinity check (GPU)
    __global__ void is_inf_g2_vec_ker(bool* results, const G2Point* points, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            results[idx] = blst_p2_is_inf(points[idx]);
        }
    }

    // Vectorized infinity check (GPU)
    batched_error_code is_inf_g2_vec_cu(bool* results, const G2Point* points, size_t num)
    {
        if (!results || !points) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(is_inf_g2_vec_ker, num, results, points);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized infinity check
    batched_error_code is_inf_g2_vec(bool* results, const G2Point* points, size_t num)
    {
#ifdef BLST_CUDA
        return is_inf_g2_vec_cu(results, points, num);
#else
        return is_inf_g2_vec_cc(results, points, num);
#endif
    }

    // Vectorized point equality check (CPU)
    batched_error_code is_equal_g2_vec_cc(bool* results, const G2Point* a, const G2Point* b, size_t num)
    {
        if (!results || !a || !b) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            *results = blst_p2_is_equal(*a, *b);
            results++; a++; b++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for point equality check (GPU)
    __global__ void is_equal_g2_vec_ker(bool* results, const G2Point* a, const G2Point* b, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            results[idx] = blst_p2_is_equal(a[idx], b[idx]);
        }
    }

    // Vectorized point equality check (GPU)
    batched_error_code is_equal_g2_vec_cu(bool* results, const G2Point* a, const G2Point* b, size_t num)
    {
        if (!results || !a || !b) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(is_equal_g2_vec_ker, num, results, a, b);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized point equality check
    batched_error_code is_equal_g2_vec(bool* results, const G2Point* a, const G2Point* b, size_t num)
    {
#ifdef BLST_CUDA
        return is_equal_g2_vec_cu(results, a, b, num);
#else
        return is_equal_g2_vec_cc(results, a, b, num);
#endif
    }

    // Forward declaration for G2 subgroup check
    //extern "C" int blst_p2_affine_in_g2(const G2AffinePoint* p);

    // Vectorized subgroup membership check (CPU)
    batched_error_code in_group_g2_vec_cc(bool* results, const G2AffinePoint* points, size_t num)
    {
        if (!results || !points) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            *results = blst_p2_affine_in_g2(*points);
            results++; points++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for subgroup membership check (GPU)
    __global__ void in_group_g2_vec_ker(bool* results, const G2AffinePoint* points, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            results[idx] = blst_p2_affine_in_g2(points[idx]);
        }
    }

    // Vectorized subgroup membership check (GPU)
    batched_error_code in_group_g2_vec_cu(bool* results, const G2AffinePoint* points, size_t num)
    {
        if (!results || !points) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(in_group_g2_vec_ker, num, results, points);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized subgroup membership check
    batched_error_code in_group_g2_vec(bool* results, const G2AffinePoint* points, size_t num)
    {
#ifdef BLST_CUDA
        return in_group_g2_vec_cu(results, points, num);
#else
        return in_group_g2_vec_cc(results, points, num);
#endif
    }

    ///////////////////////////////////////////
    // Serialization Operations

    // Vectorized point serialization (CPU)
    batched_error_code serialize_g2_vec_cc(unsigned char* out, const G2Point* points, size_t num)
    {
        if (!out || !points) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            blst_p2_serialize(out, *points);
            out += 192; points++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for point serialization (GPU)
    __global__ void serialize_g2_vec_ker(unsigned char* out, const G2Point* points, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            blst_p2_serialize(&out[idx * 192], points[idx]);
        }
    }

    // Vectorized point serialization (GPU)
    batched_error_code serialize_g2_vec_cu(unsigned char* out, const G2Point* points, size_t num)
    {
        if (!out || !points) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(serialize_g2_vec_ker, num, out, points);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized point serialization
    batched_error_code serialize_g2_vec(unsigned char* out, const G2Point* points, size_t num)
    {
#ifdef BLST_CUDA
        return serialize_g2_vec_cu(out, points, num);
#else
        return serialize_g2_vec_cc(out, points, num);
#endif
    }

    // Vectorized point deserialization (CPU)
    batched_error_code deserialize_g2_vec_cc(G2AffinePoint* out, const unsigned char* serialized, BLST_ERROR* errors, size_t num)
    {
        if (!out || !serialized || !errors) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            *errors = blst_p2_deserialize(*out, serialized);
            out++; serialized += 192; errors++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for point deserialization (GPU)
    __global__ void deserialize_g2_vec_ker(G2AffinePoint* out, const unsigned char* serialized, BLST_ERROR* errors, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            errors[idx] = blst_p2_deserialize(out[idx], &serialized[idx * 192]);
        }
    }

    // Vectorized point deserialization (GPU)
    batched_error_code deserialize_g2_vec_cu(G2AffinePoint* out, const unsigned char* serialized, BLST_ERROR* errors, size_t num)
    {
        if (!out || !serialized || !errors) return BATCHED_ERROR_NULL_POINTER;

        BLST_ERROR* d_errors;
        cudaError_t cuda_err = cudaMalloc(&d_errors, num * sizeof(BLST_ERROR));
        if (cuda_err != cudaSuccess) {
            return cuda_to_batched(cuda_err);
        }

        auto future = launch_auto(deserialize_g2_vec_ker, num, out, serialized, d_errors);
        batched_error_code sync_result = cuda_to_batched(future.get());
        if (sync_result != BATCHED_SUCCESS) {
            cudaFree(d_errors);
            return sync_result;
        }

        cuda_err = cudaMemcpy(errors, d_errors, num * sizeof(BLST_ERROR), cudaMemcpyDeviceToHost);
        cudaFree(d_errors);

        if (cuda_err != cudaSuccess) {
            return cuda_to_batched(cuda_err);
        }

        return BATCHED_SUCCESS;
    }
#endif

    // Vectorized point deserialization
    batched_error_code deserialize_g2_vec(G2AffinePoint* out, const unsigned char* serialized, BLST_ERROR* errors, size_t num)
    {
#ifdef BLST_CUDA
        return deserialize_g2_vec_cu(out, serialized, errors, num);
#else
        return deserialize_g2_vec_cc(out, serialized, errors, num);
#endif
    }

    // Vectorized point compression (CPU)
    batched_error_code compress_g2_vec_cc(unsigned char* out, const G2Point* points, size_t num)
    {
        if (!out || !points) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            blst_p2_compress(out, *points);
            out += 96; points++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for point compression (GPU)
    __global__ void compress_g2_vec_ker(unsigned char* out, const G2Point* points, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            blst_p2_compress(&out[idx * 96], points[idx]);
        }
    }

    // Vectorized point compression (GPU)
    batched_error_code compress_g2_vec_cu(unsigned char* out, const G2Point* points, size_t num)
    {
        if (!out || !points) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(compress_g2_vec_ker, num, out, points);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized point compression
    batched_error_code compress_g2_vec(unsigned char* out, const G2Point* points, size_t num)
    {
#ifdef BLST_CUDA
        return compress_g2_vec_cu(out, points, num);
#else
        return compress_g2_vec_cc(out, points, num);
#endif
    }

    // Vectorized point decompression (CPU)
    batched_error_code uncompress_g2_vec_cc(G2AffinePoint* out, const unsigned char* compressed, BLST_ERROR* errors, size_t num)
    {
        if (!out || !compressed || !errors) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            *errors = blst_p2_uncompress(*out, compressed);
            out++; compressed += 96; errors++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for point decompression (GPU)
    __global__ void uncompress_g2_vec_ker(G2AffinePoint* out, const unsigned char* compressed, BLST_ERROR* errors, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            errors[idx] = blst_p2_uncompress(out[idx], &compressed[idx * 96]);
        }
    }

    // Vectorized point decompression (GPU)
    batched_error_code uncompress_g2_vec_cu(G2AffinePoint* out, const unsigned char* compressed, BLST_ERROR* errors, size_t num)
    {
        if (!out || !compressed || !errors) return BATCHED_ERROR_NULL_POINTER;

        BLST_ERROR* d_errors;
        cudaError_t cuda_err = cudaMalloc(&d_errors, num * sizeof(BLST_ERROR));
        if (cuda_err != cudaSuccess) {
            return cuda_to_batched(cuda_err);
        }

        auto future = launch_auto(uncompress_g2_vec_ker, num, out, compressed, d_errors);
        batched_error_code sync_result = cuda_to_batched(future.get());
        if (sync_result != BATCHED_SUCCESS) {
            cudaFree(d_errors);
            return sync_result;
        }

        cuda_err = cudaMemcpy(errors, d_errors, num * sizeof(BLST_ERROR), cudaMemcpyDeviceToHost);
        cudaFree(d_errors);

        if (cuda_err != cudaSuccess) {
            return cuda_to_batched(cuda_err);
        }

        return BATCHED_SUCCESS;
    }
#endif

    // Vectorized point decompression
    batched_error_code uncompress_g2_vec(G2AffinePoint* out, const unsigned char* compressed, BLST_ERROR* errors, size_t num)
    {
#ifdef BLST_CUDA
        return uncompress_g2_vec_cu(out, compressed, errors, num);
#else
        return uncompress_g2_vec_cc(out, compressed, errors, num);
#endif
    }

    ///////////////////////////////////////////
    // Advanced Operations

    // Vectorized hash-to-curve (CPU)
    batched_error_code hash_to_g2_vec_cc(G2Point* out, const unsigned char* const* msgs, const size_t* msg_lens,
        const unsigned char* dst, size_t dst_len, size_t num)
    {
        if (!out || !msgs || !msg_lens || !dst) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            if (!msgs[i]) return BATCHED_ERROR_NULL_POINTER;
            Hash_to_curve<G2Point>(out[i], msgs[i], msg_lens[i], dst, dst_len, NULL, 0);
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for hash-to-curve (GPU)
    __global__ void hash_to_g2_vec_ker(G2Point* out, const unsigned char* const* msgs, const size_t* msg_lens,
        const unsigned char* dst, size_t dst_len, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            Hash_to_curve<G2Point>(out[idx], msgs[idx], msg_lens[idx], dst, dst_len, NULL, 0);
        }
    }

    // Vectorized hash-to-curve (GPU)
    batched_error_code hash_to_g2_vec_cu(G2Point* out, const unsigned char* const* msgs, const size_t* msg_lens,
        const unsigned char* dst, size_t dst_len, size_t num)
    {
        if (!out || !msgs || !msg_lens || !dst) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(hash_to_g2_vec_ker, num, out, msgs, msg_lens, dst, dst_len);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized hash-to-curve
    batched_error_code hash_to_g2_vec(G2Point* out, const unsigned char* const* msgs, const size_t* msg_lens,
        const unsigned char* dst, size_t dst_len, size_t num)
    {
#ifdef BLST_CUDA
        return hash_to_g2_vec_cu(out, msgs, msg_lens, dst, dst_len, num);
#else
        return hash_to_g2_vec_cc(out, msgs, msg_lens, dst, dst_len, num);
#endif
    }

    // Multi-scalar multiplication (CPU) - simplified version
    batched_error_code multi_scalar_mult_g2_cc(G2Point* result, const G2Point* points,
        const unsigned char* scalars, size_t scalar_bits, size_t num)
    {
        if (!result || !points || !scalars) return BATCHED_ERROR_NULL_POINTER;

        // Zero the result
        vec_zero(result, sizeof(*result));

        const size_t scalar_bytes = (scalar_bits + 7) / 8;

        // Simple MSM: accumulate scalar multiplications
        for (size_t i = 0; i < num; i++) {
            G2Point temp;
            blst_p2_mult(temp, points[i], &scalars[i * scalar_bytes], scalar_bits);
            if (i == 0) {
                *result = temp;
            }
            else {
                blst_p2_add_or_double(*result, *result, temp);
            }
        }

        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Multi-scalar multiplication (GPU) - simplified version
    batched_error_code multi_scalar_mult_g2_cu(G2Point* result, const G2Point* points,
        const unsigned char* scalars, size_t scalar_bits, size_t num)
    {
        // For now, just delegate to CPU implementation
        // A full GPU MSM would require more sophisticated algorithms like Pippenger's
        return multi_scalar_mult_g2_cc(result, points, scalars, scalar_bits, num);
    }
#endif

    // Multi-scalar multiplication
    batched_error_code multi_scalar_mult_g2(G2Point* result, const G2Point* points,
        const unsigned char* scalars, size_t scalar_bits, size_t num)
    {
#ifdef BLST_CUDA
        return multi_scalar_mult_g2_cu(result, points, scalars, scalar_bits, num);
#else
        return multi_scalar_mult_g2_cc(result, points, scalars, scalar_bits, num);
#endif
    }

    // Vectorized mixed addition (CPU)
    batched_error_code add_mixed_g2_vec_cc(G2Point* out, const G2Point* a, const G2AffinePoint* b, size_t num)
    {
        if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            blst_p2_add_affine(*out, *a, *b);
            out++; a++; b++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for mixed addition (GPU)
    __global__ void add_mixed_g2_vec_ker(G2Point* out, const G2Point* a, const G2AffinePoint* b, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            blst_p2_add_affine(out[idx], a[idx], b[idx]);
        }
    }

    // Vectorized mixed addition (GPU)
    batched_error_code add_mixed_g2_vec_cu(G2Point* out, const G2Point* a, const G2AffinePoint* b, size_t num)
    {
        if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(add_mixed_g2_vec_ker, num, out, a, b);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized mixed addition
    batched_error_code add_mixed_g2_vec(G2Point* out, const G2Point* a, const G2AffinePoint* b, size_t num)
    {
#ifdef BLST_CUDA
        return add_mixed_g2_vec_cu(out, a, b, num);
#else
        return add_mixed_g2_vec_cc(out, a, b, num);
#endif
    }


    // Vectorized psi endomorphism (CPU)
    batched_error_code psi_g2_vec_cc(G2Point* out, const G2Point* a, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        while (num--) {
            psi(*out, *a);
            out++; a++;
        }
        return BATCHED_SUCCESS;
    }

#ifdef BLST_CUDA
    // Kernel for psi endomorphism (GPU)
    __global__ void psi_g2_vec_ker(G2Point* out, const G2Point* a, size_t num)
    {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            psi(out[idx], a[idx]);
        }
    }

    // Vectorized psi endomorphism (GPU)
    batched_error_code psi_g2_vec_cu(G2Point* out, const G2Point* a, size_t num)
    {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto future = launch_auto(psi_g2_vec_ker, num, out, a);
        return cuda_to_batched(future.get());
    }
#endif

    // Vectorized psi endomorphism
    batched_error_code psi_g2_vec(G2Point* out, const G2Point* a, size_t num)
    {
#ifdef BLST_CUDA
        return psi_g2_vec_cu(out, a, num);
#else
        return psi_g2_vec_cc(out, a, num);
#endif
    }
}