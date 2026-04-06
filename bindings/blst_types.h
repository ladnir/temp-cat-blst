/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef __BLST_TYPES_H__
#define __BLST_TYPES_H__

#ifdef __SIZE_TYPE__
typedef __SIZE_TYPE__ size_t;
#else
#include <stddef.h>
#endif

#if defined(__UINT8_TYPE__) && defined(__UINT32_TYPE__) && defined(__UINT64_TYPE__)
typedef __UINT8_TYPE__ uint8_t;
typedef __UINT32_TYPE__ uint32_t;
typedef __UINT64_TYPE__ uint64_t;
#else
#include <stdint.h>
#endif

#if !defined(__cplusplus) && (!defined(__STDC_VERSION__) || __STDC_VERSION__ < 202311L)
# if !defined(bool)
#  if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#   define bool _Bool
#  else
#   define bool int
#  endif
# endif
#endif

/**
 * Stable error codes returned by the plain C API.
 *
 * Invariants:
 * - Numeric values are part of the ABI.
 * - The same codes are used by the scalar, AVX, and middle-product headers.
 */
typedef enum {
    BLST_SUCCESS = 0,
    BLST_BAD_ENCODING,
    BLST_POINT_NOT_ON_CURVE,
    BLST_POINT_NOT_IN_GROUP,
    BLST_AGGR_TYPE_MISMATCH,
    BLST_VERIFY_FAIL,
    BLST_PK_IS_INFINITY,
    BLST_BAD_SCALAR,
} BLST_ERROR;

typedef uint8_t byte;
typedef uint64_t limb_t;

#if defined(_MSC_VER)
# define BLST_AUX_ALIGNAS(n) __declspec(align(n))
#elif defined(__GNUC__) || defined(__clang__)
# define BLST_AUX_ALIGNAS(n) __attribute__((aligned(n)))
#else
# define BLST_AUX_ALIGNAS(n)
#endif

/**
 * Scalar value encoded as 32 opaque bytes.
 *
 * Invariants:
 * - The byte layout is ABI-stable.
 * - Callers may store and move this value as plain POD bytes.
 */
typedef struct { byte b[256 / 8]; } blst_scalar;

/**
 * Scalar field element in Montgomery form.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - All arithmetic functions interpret values in Montgomery form.
 */
typedef struct { limb_t l[256 / 8 / sizeof(limb_t)]; } blst_fr;

/**
 * Base field element in Montgomery form.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - All arithmetic functions interpret values in Montgomery form.
 */
typedef struct { limb_t l[384 / 8 / sizeof(limb_t)]; } blst_fp;

/**
 * Quadratic extension field element in Montgomery form.
 */
typedef struct { blst_fp fp[2]; } blst_fp2;

/**
 * Sextic extension field element in Montgomery form.
 */
typedef struct { blst_fp2 fp2[3]; } blst_fp6;

/**
 * Twelfth-degree extension field element in Montgomery form.
 */
typedef struct { blst_fp6 fp6[2]; } blst_fp12;

/**
 * G1 projective point over Montgomery Fp coordinates.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - The point may be the point at infinity.
 */
typedef struct { blst_fp x, y, z; } blst_p1;

/**
 * G1 affine point over Montgomery Fp coordinates.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - The point may be the point at infinity.
 */
typedef struct { blst_fp x, y; } blst_p1_affine;

/**
 * G2 projective point over Montgomery Fp2 coordinates.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - The point may be the point at infinity.
 */
typedef struct { blst_fp2 x, y, z; } blst_p2;

/**
 * G2 affine point over Montgomery Fp2 coordinates.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - The point may be the point at infinity.
 */
typedef struct { blst_fp2 x, y; } blst_p2_affine;

#endif
