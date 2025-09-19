/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef __BLS12_381_ASM_ERRORS_H__
#define __BLS12_381_ASM_ERRORS_H__

#ifdef BLST_CUDA
#include <cuda_runtime.h>
#endif


#include "blst/config.h"
typedef enum {
    BLST_SUCCESS = 0,
    BLST_BAD_ENCODING,
    BLST_POINT_NOT_ON_CURVE,
    BLST_POINT_NOT_IN_GROUP,
    BLST_AGGR_TYPE_MISMATCH,
    BLST_VERIFY_FAIL,
    BLST_PK_IS_INFINITY,
} BLST_ERROR;

// Error codes for CUDA field operations - focused on most common errors
typedef enum {
	BATCHED_SUCCESS = 0,
	BATCHED_ERROR_NULL_POINTER,
	BATCHED_ERROR_INVALID_POINTER,
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
	case BATCHED_ERROR_INVALID_POINTER:
		return "Invalid pointer (e.g. wrong host vs device)";
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

#ifdef BLST_CUDA

// Convert CUDA error to batched error code - focusing on most common errors
BLST_INLINE batched_error_code cuda_to_batched(cudaError_t error)
{
	switch (error) {
	case cudaSuccess:
		return BATCHED_SUCCESS;
	case cudaErrorInvalidValue:
		return BATCHED_ERROR_NULL_POINTER;
	case cudaErrorMemoryAllocation:
		return BATCHED_ERROR_CUDA_MALLOC_FAILED;
	case cudaErrorInvalidDevicePointer:
		return BATCHED_ERROR_NULL_POINTER;
	case cudaErrorInvalidHostPointer:
		return BATCHED_ERROR_NULL_POINTER;
	case cudaErrorInvalidMemcpyDirection:
		return BATCHED_ERROR_CUDA_MEMCPY_FAILED;
	case cudaErrorInvalidConfiguration:
		return BATCHED_ERROR_CUDA_KERNEL_LAUNCH_FAILED;
	case cudaErrorLaunchFailure:
		return BATCHED_ERROR_CUDA_KERNEL_LAUNCH_FAILED;
	case cudaErrorLaunchTimeout:
		return BATCHED_ERROR_CUDA_KERNEL_LAUNCH_FAILED;
	case cudaErrorLaunchOutOfResources:
		return BATCHED_ERROR_CUDA_KERNEL_LAUNCH_FAILED;
	case cudaErrorInvalidDevice:
		return BATCHED_ERROR_CUDA_DEVICE_ERROR;
	case cudaErrorDeviceUninitialized:
		return BATCHED_ERROR_CUDA_DEVICE_ERROR;
	case cudaErrorNoDevice:
		return BATCHED_ERROR_CUDA_DEVICE_ERROR;
	case cudaErrorNotReady:
		return BATCHED_ERROR_CUDA_SYNC_FAILED;
	case cudaErrorIllegalAddress:
		return BATCHED_ERROR_ILLEGAL_ADDRESS;
	case cudaErrorInvalidKernelImage:
		return BATCHED_ERROR_CUDA_KERNEL_LAUNCH_FAILED;
	case cudaErrorHostMemoryAlreadyRegistered:
		return BATCHED_ERROR_HOST_MALLOC_FAILED;
	case cudaErrorHostMemoryNotRegistered:
		return BATCHED_ERROR_HOST_MALLOC_FAILED;
	default:
		return BATCHED_ERROR_UNKNOWN;
	}
}
#endif
#endif
