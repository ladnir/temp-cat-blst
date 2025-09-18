/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/Fp12_cuda.h"

#ifdef __CUDACC__

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

 // CUDA kernel for vectorized Fp12 multiplication
__global__ void mul_fp12_kernel(vec384fp12* out, const vec384fp12* a, const vec384fp12* b, size_t num) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        mul_fp12(out[idx], a[idx], b[idx]);
    }
}

// CUDA kernel for vectorized Fp12 squaring
__global__ void sqr_fp12_kernel(vec384fp12* out, const vec384fp12* a, size_t num) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        sqr_fp12(out[idx], a[idx]);
    }
}

// CUDA kernel for vectorized Fp12 cyclotomic squaring
__global__ void cyclotomic_sqr_fp12_kernel(vec384fp12* out, const vec384fp12* a, size_t num) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        cyclotomic_sqr_fp12(out[idx], a[idx]);
    }
}

// CUDA kernel for vectorized Fp12 conjugation
__global__ void conjugate_fp12_kernel(vec384fp12* out, const vec384fp12* a, size_t num) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        // Copy input to output first, then conjugate in-place
        for (int i = 0; i < sizeof(vec384fp12) / sizeof(limb_t); i++) {
            ((limb_t*)out[idx])[i] = ((const limb_t*)a[idx])[i];
        }
        conjugate_fp12(out[idx]);
    }
}

// CUDA kernel for vectorized Fp12 inversion
__global__ void inverse_fp12_kernel(vec384fp12* out, const vec384fp12* a, size_t num) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        inverse_fp12(out[idx], a[idx]);
    }
}

// CUDA kernel for vectorized Fp12 Frobenius map
__global__ void frobenius_map_fp12_kernel(vec384fp12* out, const vec384fp12* a, size_t power, size_t num) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        frobenius_map_fp12(out[idx], a[idx], power);
    }
}

// CUDA kernel for vectorized Fp12 equality check
__global__ void is_equal_fp12_kernel(bool* out, const vec384fp12* a, const vec384fp12* b, size_t num) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        out[idx] = blst_fp12_is_equal(a[idx], b[idx]);
    }
}

// CUDA kernel for vectorized Fp12 is_one check
__global__ void is_one_fp12_kernel(bool* out, const vec384fp12* a, size_t num) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        out[idx] = blst_fp12_is_one(a[idx]);
    }
}

// CUDA kernel for vectorized Fp12 serialization to big-endian
__global__ void bendian_from_fp12_kernel(unsigned char* out, const vec384fp12* a, size_t num) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        blst_bendian_from_fp12(&out[idx * 12 * 48], a[idx]);
    }
}

// CUDA kernel for vectorized Fp12 deserialization from big-endian
__global__ void fp12_from_bendian_kernel(vec384fp12* out, const unsigned char* in, size_t num) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        blst_fp12_from_bendian(out[idx], &in[idx * 12 * 48]);
    }
}

// CUDA kernel for vectorized Fp12 multiplication by xy00z0
__global__ void mul_by_xy00z0_fp12_kernel(vec384fp12* out, const vec384fp12* a, const vec384fp6* xy00z0, size_t num) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num) {
        mul_by_xy00z0_fp12(out[idx], a[idx], xy00z0[idx]);
    }
}


typedef enum {
    PTR_UNKNOWN = -1,
    PTR_UNREGISTERED_HOST = 0,
    PTR_HOST_PINNED = 1,
    PTR_DEVICE = 2,
    PTR_MANAGED = 3
} ptr_kind_t;

static inline ptr_kind_t cuda_ptr_kind(const void* p) {
    if (!p) return PTR_UNKNOWN;

    cudaPointerAttributes attr;
    cudaError_t err = cudaPointerGetAttributes(&attr, p);

    if (err != cudaSuccess) {
        // In case caller cares about sticky errors, clear it:
        cudaGetLastError();
        return PTR_UNKNOWN;
    }

#if CUDART_VERSION >= 10000
    // CUDA 10+ uses 'type'
    switch (attr.type) {
    case cudaMemoryTypeHost:       return PTR_HOST_PINNED;      // cudaHostAlloc/cudaHostRegister
    case cudaMemoryTypeDevice:     return PTR_DEVICE;           // cudaMalloc / cudaMallocAsync
    case cudaMemoryTypeManaged:    return PTR_MANAGED;          // cudaMallocManaged
    case cudaMemoryTypeUnregistered: return PTR_UNREGISTERED_HOST; // ordinary malloc/new
    default: return PTR_UNKNOWN;
    }
#else
    // CUDA 9.x and earlier
    if (attr.memoryType == cudaMemoryTypeDevice) return PTR_DEVICE;
    if (attr.isManaged)                           return PTR_MANAGED;
    // Older runtimes don’t distinguish pinned vs pageable well; assume host
    return PTR_HOST_PINNED; // or PTR_UNREGISTERED_HOST, depending on your needs
#endif
}

// Helper function to calculate optimal grid and block dimensions
static batched_error_code calculate_grid_block(size_t num_elements, dim3* grid, dim3* block) {
    //const int max_threads_per_block = 256; // Conservative choice for Fp12 complexity
    //const int max_blocks = 65535;

    //block->x = (num_elements < max_threads_per_block) ? num_elements : max_threads_per_block;
    //block->y = 1;
    //block->z = 1;

    //grid->x = (num_elements + block->x - 1) / block->x;
    //if (grid->x > max_blocks) grid->x = max_blocks;
    //grid->y = 1;
    //grid->z = 1;


    if (grid == NULL || block == NULL) {
        //printf("NULL grid/block pointers in %s\n", operation_name);
        return BATCHED_ERROR_NULL_POINTER;
    }

    if (num_elements == 0) {
        // Set minimal valid configuration for zero elements
        grid->x = grid->y = grid->z = 1;
        block->x = block->y = block->z = 1;
        return BATCHED_SUCCESS;
    }

    // Query device properties for optimal configuration
    cudaDeviceProp prop;
    cudaError_t err = cudaGetDeviceProperties(&prop, 0);

    int max_threads_per_block;
    int max_blocks;

    if (err != cudaSuccess) {
        // Use conservative defaults if query fails
        max_threads_per_block = 64;  // Very conservative for Fp12
        max_blocks = 32768;
        printf("Warning: Failed to query device properties,"
            " using conservative defaults\n");
    }
    else {
        // Use device-specific limits but cap them for Fp12 complexity
        max_threads_per_block = prop.maxThreadsPerBlock / 8; // Very conservative for Fp12
        if (max_threads_per_block > 128) max_threads_per_block = 128; // Hard cap for Fp12
        max_blocks = prop.maxGridSize[0];
        if (max_blocks > 65535) max_blocks = 65535; // Reasonable limit
    }

    // Calculate block size
    block->x = (num_elements < max_threads_per_block) ? (unsigned int)num_elements : max_threads_per_block;
    block->y = 1;
    block->z = 1;

    // Calculate grid size
    unsigned int blocks_needed = (unsigned int)((num_elements + block->x - 1) / block->x);
    grid->x = (blocks_needed > max_blocks) ? max_blocks : blocks_needed;
    grid->y = 1;
    grid->z = 1;

    // Final validation of calculated configuration
    if (grid->x == 0 || block->x == 0) {
        printf("Invalid configuration calculated : grid(%u,%u,%u) block(%u,%u,%u)\n",
            grid->x, grid->y, grid->z, block->x, block->y, block->z);
        return BATCHED_ERROR_INVALID_SIZE;
    }

    // Check for reasonable limits
    if (block->x > 1024 || grid->x > 65535) {
        printf("Configuration exceeds limits : grid(%u,%u,%u) block(%u,%u,%u)\n",
            grid->x, grid->y, grid->z, block->x, block->y, block->z);
        return BATCHED_ERROR_INVALID_SIZE;
    }

    printf("Launch config: grid(%u,%u,%u) block(%u,%u,%u) for %zu elements\n",
        grid->x, grid->y, grid->z, block->x, block->y, block->z, num_elements);

    return BATCHED_SUCCESS;
}

// Helper function to validate inputs
static batched_error_code validate_inputs(const void* out, const void* a, size_t num) {
    if (out == NULL || a == NULL) {
        return BATCHED_ERROR_NULL_POINTER;
    }
    ptr_kind_t out_kind = cuda_ptr_kind(out);
	ptr_kind_t a_kind = cuda_ptr_kind(a);
	if (out_kind != PTR_DEVICE && out_kind != PTR_MANAGED) {
        return BATCHED_ERROR_INVALID_POINTER;
	}
    if (a_kind != PTR_DEVICE && a_kind != PTR_MANAGED) {
		return BATCHED_ERROR_INVALID_POINTER;
	}
    
    return BATCHED_SUCCESS;
}

// Helper function to validate inputs for binary operations
static batched_error_code validate_binary_inputs(const void* out, const void* a, const void* b, size_t num) {
    if (out == NULL || a == NULL || b == NULL) {
        return BATCHED_ERROR_NULL_POINTER;
    }
	ptr_kind_t out_kind = cuda_ptr_kind(out);
	ptr_kind_t a_kind = cuda_ptr_kind(a);
	ptr_kind_t b_kind = cuda_ptr_kind(b);
	if (out_kind != PTR_DEVICE && out_kind != PTR_MANAGED) {
		return BATCHED_ERROR_INVALID_POINTER;
	}

    return BATCHED_SUCCESS;
}

#endif // __CUDACC__

// Implementation functions

batched_error_code mul_fp12_vec(vec384fp12* out, const vec384fp12* a, const vec384fp12* b, size_t num) {
    if (num == 0) return BATCHED_SUCCESS;

#ifdef __CUDACC__
    batched_error_code validation = validate_binary_inputs(out, a, b, num);
    if (validation != BATCHED_SUCCESS) return validation;

    dim3 grid, block;
    calculate_grid_block(num, &grid, &block);

    mul_fp12_kernel << <grid, block >> > (out, a, b, num);

    cudaError_t cuda_error = cudaGetLastError();
    if (cuda_error != cudaSuccess) {
		printf("cuda error in mul_fp12_vec: %s\n", cudaGetErrorString(cuda_error));
        return cuda_to_batched(cuda_error);
    }

    cuda_error = cudaDeviceSynchronize();
    return cuda_to_batched(cuda_error);
#else
    // CPU fallback
    batched_error_code validation = validate_binary_inputs(out, a, b, num);
    if (validation != BATCHED_SUCCESS) return validation;

    for (size_t i = 0; i < num; i++) {
        blst_fp12_mul(out[i], a[i], b[i]);
    }
    return BATCHED_SUCCESS;
#endif
}

batched_error_code sqr_fp12_vec(vec384fp12* out, const vec384fp12* a, size_t num) {
    if (num == 0) return BATCHED_SUCCESS;

#ifdef __CUDACC__
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    dim3 grid, block;
    calculate_grid_block(num, &grid, &block);

    sqr_fp12_kernel << <grid, block >> > (out, a, num);

    cudaError_t cuda_error = cudaGetLastError();
    if (cuda_error != cudaSuccess) {
        return cuda_to_batched(cuda_error);
    }

    cuda_error = cudaDeviceSynchronize();
    return cuda_to_batched(cuda_error);
#else
    // CPU fallback
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    for (size_t i = 0; i < num; i++) {
        blst_fp12_sqr(out[i], a[i]);
    }
    return BATCHED_SUCCESS;
#endif
}

batched_error_code cyclotomic_sqr_fp12_vec(vec384fp12* out, const vec384fp12* a, size_t num) {
    if (num == 0) return BATCHED_SUCCESS;

#ifdef __CUDACC__
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    dim3 grid, block;
    calculate_grid_block(num, &grid, &block);

    cyclotomic_sqr_fp12_kernel << <grid, block >> > (out, a, num);

    cudaError_t cuda_error = cudaGetLastError();
    if (cuda_error != cudaSuccess) {
        return cuda_to_batched(cuda_error);
    }

    cuda_error = cudaDeviceSynchronize();
    return cuda_to_batched(cuda_error);
#else
    // CPU fallback
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    for (size_t i = 0; i < num; i++) {
        blst_fp12_cyclotomic_sqr(out[i], a[i]);
    }
    return BATCHED_SUCCESS;
#endif
}

batched_error_code conjugate_fp12_vec(vec384fp12* out, const vec384fp12* a, size_t num) {
    if (num == 0) return BATCHED_SUCCESS;

#ifdef __CUDACC__
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    dim3 grid, block;
    calculate_grid_block(num, &grid, &block);

    conjugate_fp12_kernel << <grid, block >> > (out, a, num);

    cudaError_t cuda_error = cudaGetLastError();
    if (cuda_error != cudaSuccess) {
        return cuda_to_batched(cuda_error);
    }

    cuda_error = cudaDeviceSynchronize();
    return cuda_to_batched(cuda_error);
#else
    // CPU fallback
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    for (size_t i = 0; i < num; i++) {
        // Copy then conjugate
        memcpy(out[i], a[i], sizeof(vec384fp12));
        blst_fp12_conjugate(out[i]);
    }
    return BATCHED_SUCCESS;
#endif
}

batched_error_code inverse_fp12_vec(vec384fp12* out, const vec384fp12* a, size_t num) {
    if (num == 0) return BATCHED_SUCCESS;

#ifdef __CUDACC__
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    dim3 grid, block;
    calculate_grid_block(num, &grid, &block);

    inverse_fp12_kernel << <grid, block >> > (out, a, num);

    cudaError_t cuda_error = cudaGetLastError();
    if (cuda_error != cudaSuccess) {
        return cuda_to_batched(cuda_error);
    }

    cuda_error = cudaDeviceSynchronize();
    return cuda_to_batched(cuda_error);
#else
    // CPU fallback
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    for (size_t i = 0; i < num; i++) {
        blst_fp12_inverse(out[i], a[i]);
    }
    return BATCHED_SUCCESS;
#endif
}

batched_error_code frobenius_map_fp12_vec(vec384fp12* out, const vec384fp12* a, size_t power, size_t num) {
    if (num == 0) return BATCHED_SUCCESS;

#ifdef __CUDACC__
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    if (power == 0 || power > 3) {
        return BATCHED_ERROR_INVALID_SIZE;
    }

    dim3 grid, block;
    calculate_grid_block(num, &grid, &block);

    frobenius_map_fp12_kernel << <grid, block >> > (out, a, power, num);

    cudaError_t cuda_error = cudaGetLastError();
    if (cuda_error != cudaSuccess) {
        return cuda_to_batched(cuda_error);
    }

    cuda_error = cudaDeviceSynchronize();
    return cuda_to_batched(cuda_error);
#else
    // CPU fallback
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    if (power == 0 || power > 3) {
        return BATCHED_ERROR_INVALID_SIZE;
    }

    for (size_t i = 0; i < num; i++) {
        blst_fp12_frobenius_map(out[i], a[i], power);
    }
    return BATCHED_SUCCESS;
#endif
}

batched_error_code is_equal_fp12_vec(bool* out, const vec384fp12* a, const vec384fp12* b, size_t num) {
    if (num == 0) return BATCHED_SUCCESS;

#ifdef __CUDACC__
    batched_error_code validation = validate_binary_inputs(out, a, b, num);
    if (validation != BATCHED_SUCCESS) return validation;

    dim3 grid, block;
    calculate_grid_block(num, &grid, &block);

    is_equal_fp12_kernel << <grid, block >> > (out, a, b, num);

    cudaError_t cuda_error = cudaGetLastError();
    if (cuda_error != cudaSuccess) {
        return cuda_to_batched(cuda_error);
    }

    cuda_error = cudaDeviceSynchronize();
    return cuda_to_batched(cuda_error);
#else
    // CPU fallback
    batched_error_code validation = validate_binary_inputs(out, a, b, num);
    if (validation != BATCHED_SUCCESS) return validation;

    for (size_t i = 0; i < num; i++) {
        out[i] = blst_fp12_is_equal(a[i], b[i]);
    }
    return BATCHED_SUCCESS;
#endif
}

batched_error_code is_one_fp12_vec(bool* out, const vec384fp12* a, size_t num) {
    if (num == 0) return BATCHED_SUCCESS;

#ifdef __CUDACC__
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    dim3 grid, block;
    calculate_grid_block(num, &grid, &block);

    is_one_fp12_kernel << <grid, block >> > (out, a, num);

    cudaError_t cuda_error = cudaGetLastError();
    if (cuda_error != cudaSuccess) {
        return cuda_to_batched(cuda_error);
    }

    cuda_error = cudaDeviceSynchronize();
    return cuda_to_batched(cuda_error);
#else
    // CPU fallback
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    for (size_t i = 0; i < num; i++) {
        out[i] = blst_fp12_is_one(a[i]);
    }
    return BATCHED_SUCCESS;
#endif
}

batched_error_code bendian_from_fp12_vec(unsigned char* out, const vec384fp12* a, size_t num) {
    if (num == 0) return BATCHED_SUCCESS;

#ifdef __CUDACC__
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    dim3 grid, block;
    calculate_grid_block(num, &grid, &block);

    bendian_from_fp12_kernel << <grid, block >> > (out, a, num);

    cudaError_t cuda_error = cudaGetLastError();
    if (cuda_error != cudaSuccess) {
        return cuda_to_batched(cuda_error);
    }

    cuda_error = cudaDeviceSynchronize();
    return cuda_to_batched(cuda_error);
#else
    // CPU fallback
    batched_error_code validation = validate_inputs(out, a, num);
    if (validation != BATCHED_SUCCESS) return validation;

    for (size_t i = 0; i < num; i++) {
        blst_bendian_from_fp12(&out[i * 12 * 48], a[i]);
    }
    return BATCHED_SUCCESS;
#endif
}

batched_error_code fp12_from_bendian_vec(vec384fp12* out, const unsigned char* in, size_t num) {
    if (num == 0) return BATCHED_SUCCESS;

#ifdef __CUDACC__
    batched_error_code validation = validate_inputs(out, in, num);
    if (validation != BATCHED_SUCCESS) return validation;

    dim3 grid, block;
    calculate_grid_block(num, &grid, &block);

    fp12_from_bendian_kernel << <grid, block >> > (out, in, num);

    cudaError_t cuda_error = cudaGetLastError();
    if (cuda_error != cudaSuccess) {
        return cuda_to_batched(cuda_error);
    }

    cuda_error = cudaDeviceSynchronize();
    return cuda_to_batched(cuda_error);
#else
    // CPU fallback
    batched_error_code validation = validate_inputs(out, in, num);
    if (validation != BATCHED_SUCCESS) return validation;

    for (size_t i = 0; i < num; i++) {
        blst_fp12_from_bendian(out[i], &in[i * 12 * 48]);
    }
    return BATCHED_SUCCESS;
#endif
}

batched_error_code mul_by_xy00z0_fp12_vec(vec384fp12* out, const vec384fp12* a, const vec384fp6* xy00z0, size_t num) {
    if (num == 0) return BATCHED_SUCCESS;

#ifdef __CUDACC__
    batched_error_code validation = validate_binary_inputs(out, a, xy00z0, num);
    if (validation != BATCHED_SUCCESS) return validation;

    dim3 grid, block;
    calculate_grid_block(num, &grid, &block);

    mul_by_xy00z0_fp12_kernel << <grid, block >> > (out, a, xy00z0, num);

    cudaError_t cuda_error = cudaGetLastError();
    if (cuda_error != cudaSuccess) {
        return cuda_to_batched(cuda_error);
    }

    cuda_error = cudaDeviceSynchronize();
    return cuda_to_batched(cuda_error);
#else
    // CPU fallback
    batched_error_code validation = validate_binary_inputs(out, a, xy00z0, num);
    if (validation != BATCHED_SUCCESS) return validation;

    for (size_t i = 0; i < num; i++) {
        blst_fp12_mul_by_xy00z0(out[i], a[i], xy00z0[i]);
    }
    return BATCHED_SUCCESS;
#endif
}