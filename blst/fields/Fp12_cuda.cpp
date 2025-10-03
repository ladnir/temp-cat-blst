/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/fields/Fp12_cuda.h"
#include "blst/common/cuda_launch.h"

namespace blst
{

#ifdef BLST_CUDA

    // CUDA kernel for vectorized Fp12 multiplication
    GLOBAL_FUNC void mul_fp12_kernel(Fp12* out, const Fp12* a, const Fp12* b, size_t num) {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            Fp12::mul(out[idx], a[idx], b[idx]);
        }
    }

    // CUDA kernel for vectorized Fp12 squaring
    GLOBAL_FUNC void sqr_fp12_kernel(Fp12* out, const Fp12* a, size_t num) {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            Fp12::sqr(out[idx], a[idx]);
        }
    }

    // CUDA kernel for vectorized Fp12 cyclotomic squaring
    GLOBAL_FUNC void cyclotomic_sqr_fp12_kernel(Fp12* out, const Fp12* a, size_t num) {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            Fp12::cyclotomic_sqr(out[idx], a[idx]);
        }
    }

    // CUDA kernel for vectorized Fp12 conjugation
    GLOBAL_FUNC void conjugate_fp12_kernel(Fp12* out, const Fp12* a, size_t num) {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            // Copy input to output first, then conjugate in-place
            for (int i = 0; i < sizeof(Fp12) / sizeof(limb_t); i++) {
                ((limb_t*)out[idx].data)[i] = ((const limb_t*)a[idx].data)[i];
            }
            Fp12::conjugate(out[idx]);
        }
    }

    // CUDA kernel for vectorized Fp12 inversion
    GLOBAL_FUNC void inverse_fp12_kernel(Fp12* out, const Fp12* a, size_t num) {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            Fp12::inverse(out[idx], a[idx]);
        }
    }

    // CUDA kernel for vectorized Fp12 Frobenius map
    GLOBAL_FUNC void frobenius_map_fp12_kernel(Fp12* out, const Fp12* a, size_t power, size_t num) {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            Fp12::frobenius_map(out[idx], a[idx], power);
        }
    }

    // CUDA kernel for vectorized Fp12 equality check
    GLOBAL_FUNC void is_equal_fp12_kernel(bool* out, const Fp12* a, const Fp12* b, size_t num) {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            out[idx] = Fp12::equal(a[idx], b[idx]);
        }
    }

    // CUDA kernel for vectorized Fp12 is_one check
    GLOBAL_FUNC void is_one_fp12_kernel(bool* out, const Fp12* a, size_t num) {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            out[idx] = Fp12::is_one(a[idx]);
        }
    }

    // CUDA kernel for vectorized Fp12 serialization to big-endian
    GLOBAL_FUNC void bendian_from_fp12_kernel(unsigned char* out, const Fp12* a, size_t num) {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            Fp12::to_bendian(&out[idx * 12 * 48], a[idx]);
        }
    }

    // CUDA kernel for vectorized Fp12 deserialization from big-endian
    GLOBAL_FUNC void fp12_from_bendian_kernel(Fp12* out, const unsigned char* in, size_t num) {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            Fp12::from_bendian(out[idx], &in[idx * 12 * 48]);
        }
    }

    // CUDA kernel for vectorized Fp12 multiplication by xy00z0
    GLOBAL_FUNC void mul_by_xy00z0_fp12_kernel(Fp12* out, const Fp12* a, const Fp6* xy00z0, size_t num) {
        size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < num) {
            Fp12::mul_by_xy00z0(out[idx], a[idx], xy00z0[idx]);
        }
    }

#endif // BLST_CUDA

    // CPU fallback implementations
    batched_error_code mul_fp12_vec_cc(Fp12* out, const Fp12* a, const Fp12* b, size_t num) {
        if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            Fp12::mul(out[i], a[i], b[i]);
        }
        return BATCHED_SUCCESS;
    }

    batched_error_code sqr_fp12_vec_cc(Fp12* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            Fp12::sqr(out[i], a[i]);
        }
        return BATCHED_SUCCESS;
    }

    batched_error_code cyclotomic_sqr_fp12_vec_cc(Fp12* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            Fp12::cyclotomic_sqr(out[i], a[i]);
        }
        return BATCHED_SUCCESS;
    }

    batched_error_code conjugate_fp12_vec_cc(Fp12* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            memcpy(&out[i], &a[i], sizeof(Fp12));
            Fp12::conjugate(out[i]);
        }
        return BATCHED_SUCCESS;
    }

    batched_error_code inverse_fp12_vec_cc(Fp12* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            Fp12::inverse(out[i], a[i]);
        }
        return BATCHED_SUCCESS;
    }

    batched_error_code frobenius_map_fp12_vec_cc(Fp12* out, const Fp12* a, size_t power, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            Fp12::frobenius_map(out[i], a[i], power);
        }
        return BATCHED_SUCCESS;
    }

    batched_error_code is_equal_fp12_vec_cc(bool* out, const Fp12* a, const Fp12* b, size_t num) {
        if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            out[i] = Fp12::equal(a[i], b[i]);
        }
        return BATCHED_SUCCESS;
    }

    batched_error_code is_one_fp12_vec_cc(bool* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            out[i] = Fp12::is_one(a[i]);
        }
        return BATCHED_SUCCESS;
    }

    batched_error_code bendian_from_fp12_vec_cc(unsigned char* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            Fp12::to_bendian(&out[i * 12 * 48], a[i]);
        }
        return BATCHED_SUCCESS;
    }

    batched_error_code fp12_from_bendian_vec_cc(Fp12* out, const unsigned char* in, size_t num) {
        if (!out || !in) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            Fp12::from_bendian(out[i], &in[i * 12 * 48]);
        }
        return BATCHED_SUCCESS;
    }

    batched_error_code mul_by_xy00z0_fp12_vec_cc(Fp12* out, const Fp12* a, const Fp6* xy00z0, size_t num) {
        if (!out || !a || !xy00z0) return BATCHED_ERROR_NULL_POINTER;

        for (size_t i = 0; i < num; i++) {
            Fp12::mul_by_xy00z0(out[i], a[i], xy00z0[i]);
        }
        return BATCHED_SUCCESS;
    }

    // CUDA implementations
#ifdef BLST_CUDA
    batched_error_code mul_fp12_vec_cu(Fp12* out, const Fp12* a, const Fp12* b, size_t num) {
        if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

        auto fut = launch_auto(mul_fp12_kernel, num, out, a, b);
        return cuda_to_batched(fut.get());
    }

    batched_error_code sqr_fp12_vec_cu(Fp12* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto fut = launch_auto(sqr_fp12_kernel, num, out, a);
        return cuda_to_batched(fut.get());
    }

    batched_error_code cyclotomic_sqr_fp12_vec_cu(Fp12* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto fut = launch_auto(cyclotomic_sqr_fp12_kernel, num, out, a);
        return cuda_to_batched(fut.get());
    }

    batched_error_code conjugate_fp12_vec_cu(Fp12* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto fut = launch_auto(conjugate_fp12_kernel, num, out, a);
        return cuda_to_batched(fut.get());
    }

    batched_error_code inverse_fp12_vec_cu(Fp12* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto fut = launch_auto(inverse_fp12_kernel, num, out, a);
        return cuda_to_batched(fut.get());
    }

    batched_error_code frobenius_map_fp12_vec_cu(Fp12* out, const Fp12* a, size_t power, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto fut = launch_auto(frobenius_map_fp12_kernel, num, out, a, power);
        return cuda_to_batched(fut.get());
    }

    batched_error_code is_equal_fp12_vec_cu(bool* out, const Fp12* a, const Fp12* b, size_t num) {
        if (!out || !a || !b) return BATCHED_ERROR_NULL_POINTER;

        auto fut = launch_auto(is_equal_fp12_kernel, num, out, a, b);
        return cuda_to_batched(fut.get());
    }

    batched_error_code is_one_fp12_vec_cu(bool* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto fut = launch_auto(is_one_fp12_kernel, num, out, a);
        return cuda_to_batched(fut.get());
    }

    batched_error_code bendian_from_fp12_vec_cu(unsigned char* out, const Fp12* a, size_t num) {
        if (!out || !a) return BATCHED_ERROR_NULL_POINTER;

        auto fut = launch_auto(bendian_from_fp12_kernel, num, out, a);
        return cuda_to_batched(fut.get());
    }

    batched_error_code fp12_from_bendian_vec_cu(Fp12* out, const unsigned char* in, size_t num) {
        if (!out || !in) return BATCHED_ERROR_NULL_POINTER;

        auto fut = launch_auto(fp12_from_bendian_kernel, num, out, in);
        return cuda_to_batched(fut.get());
    }

    batched_error_code mul_by_xy00z0_fp12_vec_cu(Fp12* out, const Fp12* a, const Fp6* xy00z0, size_t num) {
        if (!out || !a || !xy00z0) return BATCHED_ERROR_NULL_POINTER;

        auto fut = launch_auto(mul_by_xy00z0_fp12_kernel, num, out, a, xy00z0);
        return cuda_to_batched(fut.get());
    }
#endif // BLST_CUDA

    // Public API functions
    batched_error_code mul_fp12_vec(Fp12* out, const Fp12* a, const Fp12* b, size_t num) {
        if (num == 0) return BATCHED_SUCCESS;

#ifdef BLST_CUDA
        return mul_fp12_vec_cu(out, a, b, num);
#else
        return mul_fp12_vec_cc(out, a, b, num);
#endif
    }

    batched_error_code sqr_fp12_vec(Fp12* out, const Fp12* a, size_t num) {
        if (num == 0) return BATCHED_SUCCESS;

#ifdef BLST_CUDA
        return sqr_fp12_vec_cu(out, a, num);
#else
        return sqr_fp12_vec_cc(out, a, num);
#endif
    }

    batched_error_code cyclotomic_sqr_fp12_vec(Fp12* out, const Fp12* a, size_t num) {
        if (num == 0) return BATCHED_SUCCESS;

#ifdef BLST_CUDA
        return cyclotomic_sqr_fp12_vec_cu(out, a, num);
#else
        return cyclotomic_sqr_fp12_vec_cc(out, a, num);
#endif
    }

    batched_error_code conjugate_fp12_vec(Fp12* out, const Fp12* a, size_t num) {
        if (num == 0) return BATCHED_SUCCESS;

#ifdef BLST_CUDA
        return conjugate_fp12_vec_cu(out, a, num);
#else
        return conjugate_fp12_vec_cc(out, a, num);
#endif
    }

    batched_error_code inverse_fp12_vec(Fp12* out, const Fp12* a, size_t num) {
        if (num == 0) return BATCHED_SUCCESS;

#ifdef BLST_CUDA
        return inverse_fp12_vec_cu(out, a, num);
#else
        return inverse_fp12_vec_cc(out, a, num);
#endif
    }

    batched_error_code frobenius_map_fp12_vec(Fp12* out, const Fp12* a, size_t power, size_t num) {
        if (num == 0) return BATCHED_SUCCESS;

#ifdef BLST_CUDA
        return frobenius_map_fp12_vec_cu(out, a, power, num);
#else
        return frobenius_map_fp12_vec_cc(out, a, power, num);
#endif
    }

    batched_error_code is_equal_fp12_vec(bool* out, const Fp12* a, const Fp12* b, size_t num) {
        if (num == 0) return BATCHED_SUCCESS;

#ifdef BLST_CUDA
        return is_equal_fp12_vec_cu(out, a, b, num);
#else
        return is_equal_fp12_vec_cc(out, a, b, num);
#endif
    }

    batched_error_code is_one_fp12_vec(bool* out, const Fp12* a, size_t num) {
        if (num == 0) return BATCHED_SUCCESS;

#ifdef BLST_CUDA
        return is_one_fp12_vec_cu(out, a, num);
#else
        return is_one_fp12_vec_cc(out, a, num);
#endif
    }

    batched_error_code bendian_from_fp12_vec(unsigned char* out, const Fp12* a, size_t num) {
        if (num == 0) return BATCHED_SUCCESS;

#ifdef BLST_CUDA
        return bendian_from_fp12_vec_cu(out, a, num);
#else
        return bendian_from_fp12_vec_cc(out, a, num);
#endif
    }

    batched_error_code fp12_from_bendian_vec(Fp12* out, const unsigned char* in, size_t num) {
        if (num == 0) return BATCHED_SUCCESS;

#ifdef BLST_CUDA
        return fp12_from_bendian_vec_cu(out, in, num);
#else
        return fp12_from_bendian_vec_cc(out, in, num);
#endif
    }

    batched_error_code mul_by_xy00z0_fp12_vec(Fp12* out, const Fp12* a, const Fp6* xy00z0, size_t num) {
        if (num == 0) return BATCHED_SUCCESS;

#ifdef BLST_CUDA
        return mul_by_xy00z0_fp12_vec_cu(out, a, xy00z0, num);
#else
        return mul_by_xy00z0_fp12_vec_cc(out, a, xy00z0, num);
#endif
    }
}