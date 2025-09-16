/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/curve/point.h"
#include "blst/curve/e1.h"
#include "blst/curve/G1_cuda.h"
#include "blst/curve/map_to_g1.h"
#include "blst/exports.h"

BLST_EXTERN_C_BEGIN

///////////////////////////////////////////
// Core G1 Operations

// Vectorized addition of G1 points (CPU)
batched_error_code add_g1_vec_cc(POINTonE1* out, const POINTonE1* a, const POINTonE1* b, size_t num)
{
    if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        blst_p1_add_or_double(out, a, b);
        out++; a++; b++;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for addition of G1 points (GPU)
GLOBAL_FUNC void add_g1_vec_ker(POINTonE1* out, const POINTonE1* a, const POINTonE1* b, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        blst_p1_add_or_double(&out[idx], &a[idx], &b[idx]);
    }
}

// Vectorized addition of G1 points (GPU)
batched_error_code add_g1_vec_cu(POINTonE1* out, const POINTonE1* a, const POINTonE1* b, size_t num)
{
    if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    add_g1_vec_ker << <num_blocks, threads_per_block >> > (out, a, b, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized addition of G1 points
batched_error_code add_g1_vec(POINTonE1* out, const POINTonE1* a, const POINTonE1* b, size_t num)
{
#ifdef BLST_CUDA
    return add_g1_vec_cu(out, a, b, num);
#else
    return add_g1_vec_cc(out, a, b, num);
#endif
}

// Vectorized doubling of G1 points (CPU)
batched_error_code double_g1_vec_cc(POINTonE1* out, const POINTonE1* a, size_t num)
{
    if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        blst_p1_double(out, a);
        out++; a++;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for doubling of G1 points (GPU)
GLOBAL_FUNC void double_g1_vec_ker(POINTonE1* out, const POINTonE1* a, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        blst_p1_double(&out[idx], &a[idx]);
    }
}

// Vectorized doubling of G1 points (GPU)
batched_error_code double_g1_vec_cu(POINTonE1* out, const POINTonE1* a, size_t num)
{
    if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    double_g1_vec_ker << <num_blocks, threads_per_block >> > (out, a, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized doubling of G1 points
batched_error_code double_g1_vec(POINTonE1* out, const POINTonE1* a, size_t num)
{
#ifdef BLST_CUDA
    return double_g1_vec_cu(out, a, num);
#else
    return double_g1_vec_cc(out, a, num);
#endif
}

// Vectorized conditional negation of G1 points (CPU)
batched_error_code cneg_g1_vec_cc(POINTonE1* out, const POINTonE1* a, int flag, size_t num)
{
    if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        *out = *a;
        blst_p1_cneg(out, flag);
        out++; a++;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for conditional negation of G1 points (GPU)
GLOBAL_FUNC void cneg_g1_vec_ker(POINTonE1* out, const POINTonE1* a, int flag, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        out[idx] = a[idx];
        blst_p1_cneg(&out[idx], flag);
    }
}

// Vectorized conditional negation of G1 points (GPU)
batched_error_code cneg_g1_vec_cu(POINTonE1* out, const POINTonE1* a, int flag, size_t num)
{
    if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    cneg_g1_vec_ker << <num_blocks, threads_per_block >> > (out, a, flag, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized conditional negation of G1 points
batched_error_code cneg_g1_vec(POINTonE1* out, const POINTonE1* a, int flag, size_t num)
{
#ifdef BLST_CUDA
    return cneg_g1_vec_cu(out, a, flag, num);
#else
    return cneg_g1_vec_cc(out, a, flag, num);
#endif
}

// Vectorized scalar multiplication of G1 points (CPU)
batched_error_code mult_g1_vec_cc(POINTonE1* out, const POINTonE1* points, const unsigned char* scalars, size_t scalar_bits, size_t num)
{
    if (!out || !points || !scalars) return BATCHED_ERROR_NULL_POINTER;

    const size_t scalar_bytes = (scalar_bits + 7) / 8;

    while (num--) {
        blst_p1_mult(out, points, scalars, scalar_bits);
        out++; points++; scalars += scalar_bytes;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for scalar multiplication of G1 points (GPU)
GLOBAL_FUNC void mult_g1_vec_ker(POINTonE1* out, const POINTonE1* points, const unsigned char* scalars, size_t scalar_bits, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        const size_t scalar_bytes = (scalar_bits + 7) / 8;
        blst_p1_mult(&out[idx], &points[idx], &scalars[idx * scalar_bytes], scalar_bits);
    }
}

// Vectorized scalar multiplication of G1 points (GPU)
batched_error_code mult_g1_vec_cu(POINTonE1* out, const POINTonE1* points, const unsigned char* scalars, size_t scalar_bits, size_t num)
{
    if (!out || !points || !scalars) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 128; // Smaller blocks for scalar multiplication
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    mult_g1_vec_ker << <num_blocks, threads_per_block >> > (out, points, scalars, scalar_bits, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized scalar multiplication of G1 points
batched_error_code mult_g1_vec(POINTonE1* out, const POINTonE1* points, const unsigned char* scalars, size_t scalar_bits, size_t num)
{
#ifdef BLST_CUDA
    return mult_g1_vec_cu(out, points, scalars, scalar_bits, num);
#else
    return mult_g1_vec_cc(out, points, scalars, scalar_bits, num);
#endif
}

///////////////////////////////////////////
// Coordinate Conversion Operations

// Vectorized Jacobian to affine conversion (CPU)
batched_error_code to_affine_g1_vec_cc(POINTonE1_affine* out, const POINTonE1* a, size_t num)
{
    if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        blst_p1_to_affine(out, a);
        out++; a++;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for Jacobian to affine conversion (GPU)
GLOBAL_FUNC void to_affine_g1_vec_ker(POINTonE1_affine* out, const POINTonE1* a, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        blst_p1_to_affine(&out[idx], &a[idx]);
    }
}

// Vectorized Jacobian to affine conversion (GPU)
batched_error_code to_affine_g1_vec_cu(POINTonE1_affine* out, const POINTonE1* a, size_t num)
{
    if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    to_affine_g1_vec_ker << <num_blocks, threads_per_block >> > (out, a, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized Jacobian to affine conversion
batched_error_code to_affine_g1_vec(POINTonE1_affine* out, const POINTonE1* a, size_t num)
{
#ifdef BLST_CUDA
    return to_affine_g1_vec_cu(out, a, num);
#else
    return to_affine_g1_vec_cc(out, a, num);
#endif
}

// Vectorized affine to Jacobian conversion (CPU)
batched_error_code from_affine_g1_vec_cc(POINTonE1* out, const POINTonE1_affine* a, size_t num)
{
    if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        blst_p1_from_affine(out, a);
        out++; a++;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for affine to Jacobian conversion (GPU)
GLOBAL_FUNC void from_affine_g1_vec_ker(POINTonE1* out, const POINTonE1_affine* a, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        blst_p1_from_affine(&out[idx], &a[idx]);
    }
}

// Vectorized affine to Jacobian conversion (GPU)
batched_error_code from_affine_g1_vec_cu(POINTonE1* out, const POINTonE1_affine* a, size_t num)
{
    if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    from_affine_g1_vec_ker << <num_blocks, threads_per_block >> > (out, a, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized affine to Jacobian conversion
batched_error_code from_affine_g1_vec(POINTonE1* out, const POINTonE1_affine* a, size_t num)
{
#ifdef BLST_CUDA
    return from_affine_g1_vec_cu(out, a, num);
#else
    return from_affine_g1_vec_cc(out, a, num);
#endif
}

// Vectorized Jacobian normalization (CPU)
batched_error_code from_jacobian_g1_vec_cc(POINTonE1* out, const POINTonE1* a, size_t num)
{
    if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        blst_p1_from_jacobian(out, a);
        out++; a++;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for Jacobian normalization (GPU)
GLOBAL_FUNC void from_jacobian_g1_vec_ker(POINTonE1* out, const POINTonE1* a, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        blst_p1_from_jacobian(&out[idx], &a[idx]);
    }
}

// Vectorized Jacobian normalization (GPU)
batched_error_code from_jacobian_g1_vec_cu(POINTonE1* out, const POINTonE1* a, size_t num)
{
    if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    from_jacobian_g1_vec_ker << <num_blocks, threads_per_block >> > (out, a, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized Jacobian normalization
batched_error_code from_jacobian_g1_vec(POINTonE1* out, const POINTonE1* a, size_t num)
{
#ifdef BLST_CUDA
    return from_jacobian_g1_vec_cu(out, a, num);
#else
    return from_jacobian_g1_vec_cc(out, a, num);
#endif
}

///////////////////////////////////////////
// Validation Operations

// Vectorized curve membership check (CPU)
batched_error_code on_curve_g1_vec_cc(bool* results, const POINTonE1* points, size_t num)
{
    if (!results || !points) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        *results = blst_p1_on_curve(points);
        results++; points++;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for curve membership check (GPU)
GLOBAL_FUNC void on_curve_g1_vec_ker(bool* results, const POINTonE1* points, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        results[idx] = blst_p1_on_curve(&points[idx]);
    }
}

// Vectorized curve membership check (GPU)
batched_error_code on_curve_g1_vec_cu(bool* results, const POINTonE1* points, size_t num)
{
    if (!results || !points) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    on_curve_g1_vec_ker << <num_blocks, threads_per_block >> > (results, points, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized curve membership check
batched_error_code on_curve_g1_vec(bool* results, const POINTonE1* points, size_t num)
{
#ifdef BLST_CUDA
    return on_curve_g1_vec_cu(results, points, num);
#else
    return on_curve_g1_vec_cc(results, points, num);
#endif
}

// Vectorized infinity check (CPU)
batched_error_code is_inf_g1_vec_cc(bool* results, const POINTonE1* points, size_t num)
{
    if (!results || !points) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        *results = blst_p1_is_inf(points);
        results++; points++;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for infinity check (GPU)
GLOBAL_FUNC void is_inf_g1_vec_ker(bool* results, const POINTonE1* points, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        results[idx] = blst_p1_is_inf(&points[idx]);
    }
}

// Vectorized infinity check (GPU)
batched_error_code is_inf_g1_vec_cu(bool* results, const POINTonE1* points, size_t num)
{
    if (!results || !points) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    is_inf_g1_vec_ker << <num_blocks, threads_per_block >> > (results, points, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized infinity check
batched_error_code is_inf_g1_vec(bool* results, const POINTonE1* points, size_t num)
{
#ifdef BLST_CUDA
    return is_inf_g1_vec_cu(results, points, num);
#else
    return is_inf_g1_vec_cc(results, points, num);
#endif
}

// Vectorized point equality check (CPU)
batched_error_code is_equal_g1_vec_cc(bool* results, const POINTonE1* a, const POINTonE1* b, size_t num)
{
    if (!results || !a || !b) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        *results = blst_p1_is_equal(a, b);
        results++; a++; b++;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for point equality check (GPU)
GLOBAL_FUNC void is_equal_g1_vec_ker(bool* results, const POINTonE1* a, const POINTonE1* b, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        results[idx] = blst_p1_is_equal(&a[idx], &b[idx]);
    }
}

// Vectorized point equality check (GPU)
batched_error_code is_equal_g1_vec_cu(bool* results, const POINTonE1* a, const POINTonE1* b, size_t num)
{
    if (!results || !a || !b) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    is_equal_g1_vec_ker << <num_blocks, threads_per_block >> > (results, a, b, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized point equality check
batched_error_code is_equal_g1_vec(bool* results, const POINTonE1* a, const POINTonE1* b, size_t num)
{
#ifdef BLST_CUDA
    return is_equal_g1_vec_cu(results, a, b, num);
#else
    return is_equal_g1_vec_cc(results, a, b, num);
#endif
}

///////////////////////////////////////////
// Serialization Operations

// Vectorized point serialization (CPU)
batched_error_code serialize_g1_vec_cc(unsigned char* out, const POINTonE1* points, size_t num)
{
    if (!out || !points) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        blst_p1_serialize(out, points);
        out += 96; points++;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for point serialization (GPU)
GLOBAL_FUNC void serialize_g1_vec_ker(unsigned char* out, const POINTonE1* points, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        blst_p1_serialize(&out[idx * 96], &points[idx]);
    }
}

// Vectorized point serialization (GPU)
batched_error_code serialize_g1_vec_cu(unsigned char* out, const POINTonE1* points, size_t num)
{
    if (!out || !points) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    serialize_g1_vec_ker << <num_blocks, threads_per_block >> > (out, points, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized point serialization
batched_error_code serialize_g1_vec(unsigned char* out, const POINTonE1* points, size_t num)
{
#ifdef BLST_CUDA
    return serialize_g1_vec_cu(out, points, num);
#else
    return serialize_g1_vec_cc(out, points, num);
#endif
}

// Vectorized point deserialization (CPU)
batched_error_code deserialize_g1_vec_cc(POINTonE1_affine* out, const unsigned char* serialized, size_t num)
{
    if (!out || !serialized) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        BLST_ERROR err = blst_p1_deserialize(out, serialized);
        if (err != BLST_SUCCESS) return BATCHED_ERROR_POINT_NOT_ON_CURVE;
        out++; serialized += 96;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for point deserialization (GPU)
GLOBAL_FUNC void deserialize_g1_vec_ker(POINTonE1_affine* out, const unsigned char* serialized, size_t num, bool* errors)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        BLST_ERROR err = blst_p1_deserialize(&out[idx], &serialized[idx * 96]);
        errors[idx] = (err != BLST_SUCCESS);
    }
}

// Vectorized point deserialization (GPU)
batched_error_code deserialize_g1_vec_cu(POINTonE1_affine* out, const unsigned char* serialized, size_t num)
{
    if (!out || !serialized) return BATCHED_ERROR_NULL_POINTER;

    bool* d_errors;
    bool* h_errors = (bool*)malloc(num * sizeof(bool));
    if (!h_errors) return BATCHED_ERROR_CUDA_MALLOC_FAILED;

    cudaError_t cuda_err = cudaMalloc(&d_errors, num * sizeof(bool));
    if (cuda_err != cudaSuccess) {
        free(h_errors);
        return BATCHED_ERROR_CUDA_MALLOC_FAILED;
    }

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    deserialize_g1_vec_ker << <num_blocks, threads_per_block >> > (out, serialized, num, d_errors);

    batched_error_code sync_result = cuda_to_batched(cudaDeviceSynchronize());
    if (sync_result != BATCHED_SUCCESS) {
        cudaFree(d_errors);
        free(h_errors);
        return sync_result;
    }

    cuda_err = cudaMemcpy(h_errors, d_errors, num * sizeof(bool), cudaMemcpyDeviceToHost);
    cudaFree(d_errors);

    if (cuda_err != cudaSuccess) {
        free(h_errors);
        return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
    }

    // Check for any errors
    for (size_t i = 0; i < num; i++) {
        if (h_errors[i]) {
            free(h_errors);
            return BATCHED_ERROR_POINT_NOT_ON_CURVE;
        }
    }

    free(h_errors);
    return BATCHED_SUCCESS;
}
#endif

// Vectorized point deserialization
batched_error_code deserialize_g1_vec(POINTonE1_affine* out, const unsigned char* serialized, size_t num)
{
#ifdef BLST_CUDA
    return deserialize_g1_vec_cu(out, serialized, num);
#else
    return deserialize_g1_vec_cc(out, serialized, num);
#endif
}

// Vectorized point compression (CPU)
batched_error_code compress_g1_vec_cc(unsigned char* out, const POINTonE1* points, size_t num)
{
    if (!out || !points) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        blst_p1_compress(out, points);
        out += 48; points++;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for point compression (GPU)
GLOBAL_FUNC void compress_g1_vec_ker(unsigned char* out, const POINTonE1* points, size_t num)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        blst_p1_compress(&out[idx * 48], &points[idx]);
    }
}

// Vectorized point compression (GPU)
batched_error_code compress_g1_vec_cu(unsigned char* out, const POINTonE1* points, size_t num)
{
    if (!out || !points) return BATCHED_ERROR_NULL_POINTER;

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    compress_g1_vec_ker << <num_blocks, threads_per_block >> > (out, points, num);

    return cuda_to_batched(cudaDeviceSynchronize());
}
#endif

// Vectorized point compression
batched_error_code compress_g1_vec(unsigned char* out, const POINTonE1* points, size_t num)
{
#ifdef BLST_CUDA
    return compress_g1_vec_cu(out, points, num);
#else
    return compress_g1_vec_cc(out, points, num);
#endif
}

// Vectorized point decompression (CPU)
batched_error_code uncompress_g1_vec_cc(POINTonE1_affine* out, const unsigned char* compressed, size_t num)
{
    if (!out || !compressed) return BATCHED_ERROR_NULL_POINTER;

    while (num--) {
        BLST_ERROR err = blst_p1_uncompress(out, compressed);
        if (err != BLST_SUCCESS) return BATCHED_ERROR_POINT_NOT_ON_CURVE;
        out++; compressed += 48;
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for point decompression (GPU)
GLOBAL_FUNC void uncompress_g1_vec_ker(POINTonE1_affine* out, const unsigned char* compressed, size_t num, bool* errors)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        BLST_ERROR err = blst_p1_uncompress(&out[idx], &compressed[idx * 48]);
        errors[idx] = (err != BLST_SUCCESS);
    }
}

// Vectorized point decompression (GPU)
batched_error_code uncompress_g1_vec_cu(POINTonE1_affine* out, const unsigned char* compressed, size_t num)
{
    if (!out || !compressed) return BATCHED_ERROR_NULL_POINTER;

    bool* d_errors;
    bool* h_errors = (bool*)malloc(num * sizeof(bool));
    if (!h_errors) return BATCHED_ERROR_CUDA_MALLOC_FAILED;

    cudaError_t cuda_err = cudaMalloc(&d_errors, num * sizeof(bool));
    if (cuda_err != cudaSuccess) {
        free(h_errors);
        return BATCHED_ERROR_CUDA_MALLOC_FAILED;
    }

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    uncompress_g1_vec_ker << <num_blocks, threads_per_block >> > (out, compressed, num, d_errors);

    batched_error_code sync_result = cuda_to_batched(cudaDeviceSynchronize());
    if (sync_result != BATCHED_SUCCESS) {
        cudaFree(d_errors);
        free(h_errors);
        return sync_result;
    }

    cuda_err = cudaMemcpy(h_errors, d_errors, num * sizeof(bool), cudaMemcpyDeviceToHost);
    cudaFree(d_errors);

    if (cuda_err != cudaSuccess) {
        free(h_errors);
        return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
    }

    // Check for any errors
    for (size_t i = 0; i < num; i++) {
        if (h_errors[i]) {
            free(h_errors);
            return BATCHED_ERROR_POINT_NOT_ON_CURVE;
        }
    }

    free(h_errors);
    return BATCHED_SUCCESS;
}
#endif

// Vectorized point decompression
batched_error_code uncompress_g1_vec(POINTonE1_affine* out, const unsigned char* compressed, size_t num)
{
#ifdef BLST_CUDA
    return uncompress_g1_vec_cu(out, compressed, num);
#else
    return uncompress_g1_vec_cc(out, compressed, num);
#endif
}

///////////////////////////////////////////
// Advanced Operations

// Vectorized hash-to-curve (CPU)
batched_error_code hash_to_g1_vec_cc(POINTonE1* out, const unsigned char* msgs, const size_t* msg_lens,
    const unsigned char* dst, size_t dst_len, size_t num)
{
    if (!out || !msgs || !msg_lens || !dst) return BATCHED_ERROR_NULL_POINTER;

    const unsigned char* current_msg = msgs;
    for (size_t i = 0; i < num; i++) {
        blst_hash_to_g1(&out[i], current_msg, msg_lens[i], dst, dst_len, NULL, 0);
        current_msg += msg_lens[i];
    }
    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Kernel for hash-to-curve (GPU)
GLOBAL_FUNC void hash_to_g1_vec_ker(POINTonE1* out, const unsigned char* msgs, const size_t* msg_lens,
    const unsigned char* dst, size_t dst_len, size_t num, const size_t* msg_offsets)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        blst_hash_to_g1(&out[idx], &msgs[msg_offsets[idx]], msg_lens[idx], dst, dst_len, NULL, 0);
    }
}

// Vectorized hash-to-curve (GPU)
batched_error_code hash_to_g1_vec_cu(POINTonE1* out, const unsigned char* msgs, const size_t* msg_lens,
    const unsigned char* dst, size_t dst_len, size_t num)
{
    if (!out || !msgs || !msg_lens || !dst) return BATCHED_ERROR_NULL_POINTER;

    // Calculate message offsets for GPU kernel
    size_t* h_offsets = (size_t*)malloc(num * sizeof(size_t));
    if (!h_offsets) return BATCHED_ERROR_CUDA_MALLOC_FAILED;

    h_offsets[0] = 0;
    for (size_t i = 1; i < num; i++) {
        h_offsets[i] = h_offsets[i - 1] + msg_lens[i - 1];
    }

    size_t* d_offsets;
    cudaError_t cuda_err = cudaMalloc(&d_offsets, num * sizeof(size_t));
    if (cuda_err != cudaSuccess) {
        free(h_offsets);
        return BATCHED_ERROR_CUDA_MALLOC_FAILED;
    }

    cuda_err = cudaMemcpy(d_offsets, h_offsets, num * sizeof(size_t), cudaMemcpyHostToDevice);
    if (cuda_err != cudaSuccess) {
        cudaFree(d_offsets);
        free(h_offsets);
        return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
    }

    size_t threads_per_block = 256;
    size_t num_blocks = (num + threads_per_block - 1) / threads_per_block;
    hash_to_g1_vec_ker << <num_blocks, threads_per_block >> > (out, msgs, msg_lens, dst, dst_len, num, d_offsets);

    batched_error_code result = cuda_to_batched(cudaDeviceSynchronize());

    cudaFree(d_offsets);
    free(h_offsets);
    return result;
}
#endif

// Vectorized hash-to-curve
batched_error_code hash_to_g1_vec(POINTonE1* out, const unsigned char* msgs, const size_t* msg_lens,
    const unsigned char* dst, size_t dst_len, size_t num)
{
#ifdef BLST_CUDA
    return hash_to_g1_vec_cu(out, msgs, msg_lens, dst, dst_len, num);
#else
    return hash_to_g1_vec_cc(out, msgs, msg_lens, dst, dst_len, num);
#endif
}

// Multi-scalar multiplication (CPU) - simplified version
batched_error_code multi_scalar_mult_g1_cc(POINTonE1* result, const POINTonE1* points,
    const unsigned char* scalars, size_t scalar_bits, size_t num)
{
    if (!result || !points || !scalars) return BATCHED_ERROR_NULL_POINTER;

    // Zero the result
    vec_zero(result, sizeof(*result));

    const size_t scalar_bytes = (scalar_bits + 7) / 8;

    // Simple MSM: accumulate scalar multiplications
    for (size_t i = 0; i < num; i++) {
        POINTonE1 temp;
        blst_p1_mult(&temp, &points[i], &scalars[i * scalar_bytes], scalar_bits);
        if (i == 0) {
            *result = temp;
        }
        else {
            blst_p1_add_or_double(result, result, &temp);
        }
    }

    return BATCHED_SUCCESS;
}

#ifdef BLST_CUDA
// Multi-scalar multiplication (GPU) - simplified version
batched_error_code multi_scalar_mult_g1_cu(POINTonE1* result, const POINTonE1* points,
    const unsigned char* scalars, size_t scalar_bits, size_t num)
{
    // For now, just delegate to CPU implementation
    // A full GPU MSM would require more sophisticated algorithms like Pippenger's
    return multi_scalar_mult_g1_cc(result, points, scalars, scalar_bits, num);
}
#endif

// Multi-scalar multiplication
batched_error_code multi_scalar_mult_g1(POINTonE1* result, const POINTonE1* points,
    const unsigned char* scalars, size_t scalar_bits, size_t num)
{
#ifdef BLST_CUDA
    return multi_scalar_mult_g1_cu(result, points, scalars, scalar_bits, num);
#else
    return multi_scalar_mult_g1_cc(result, points, scalars, scalar_bits, num);
#endif
}

BLST_EXTERN_C_END