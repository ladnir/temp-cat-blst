/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef __BLST_AVX_H__
#define __BLST_AVX_H__

#include "blst_types.h"

#ifdef __cplusplus
extern "C" {
#endif

enum {
    BLST_AVX_ALIGNMENT = 64,
    BLST_AVX_P1_SIZE = 384,
    BLST_AVX_P1_AFFINE_SIZE = 256,
    BLST_AVX_P2_SIZE = 768,
    BLST_AVX_P2_AFFINE_SIZE = 512,
    BLST_AVX_FP12_SIZE = 1536,
};

/**
 * AVX-native POD shells.
 *
 * Invariants:
 * - These layouts are stable for the current expert AVX ABI.
 * - The byte arrays are intentionally opaque to C and Rust callers.
 * - Callers may allocate storage of this size and alignment and pass it back
 *   to the library, but must not inspect or mutate the internal
 *   representation except through these entry points.
 * - The library intentionally reinterprets these shells as its internal AVX
 *   types after validating size, alignment, and layout assumptions.
 * - The storage is aligned for native AVX builds; portable AVX builds may use
 *   a weaker internal alignment while remaining representable by these shells.
 */
typedef struct BLST_AUX_ALIGNAS(BLST_AVX_ALIGNMENT) { byte b[BLST_AVX_P1_SIZE]; } blst_avx_p1;
typedef struct BLST_AUX_ALIGNAS(BLST_AVX_ALIGNMENT) { byte b[BLST_AVX_P1_AFFINE_SIZE]; } blst_avx_p1_affine;
typedef struct BLST_AUX_ALIGNAS(BLST_AVX_ALIGNMENT) { byte b[BLST_AVX_P2_SIZE]; } blst_avx_p2;
typedef struct BLST_AUX_ALIGNAS(BLST_AVX_ALIGNMENT) { byte b[BLST_AVX_P2_AFFINE_SIZE]; } blst_avx_p2_affine;
typedef struct BLST_AUX_ALIGNAS(BLST_AVX_ALIGNMENT) { byte b[BLST_AVX_FP12_SIZE]; } blst_avx_fp12;

/**
 * Reports whether the current build contains the expert AVX API.
 *
 * Return value:
 * - `true` when the library was compiled with native or portable AVX support.
 * - `false` when all AVX entry points are inert stubs.
 */
bool blst_avx_enabled(void);

/** Returns the ABI size in bytes of `blst_avx_p1`. */
size_t blst_avx_p1_sizeof(void);
/** Returns the required ABI alignment of `blst_avx_p1`. */
size_t blst_avx_p1_alignof(void);
/** Returns the ABI size in bytes of `blst_avx_p1_affine`. */
size_t blst_avx_p1_affine_sizeof(void);
/** Returns the required ABI alignment of `blst_avx_p1_affine`. */
size_t blst_avx_p1_affine_alignof(void);
/** Returns the ABI size in bytes of `blst_avx_p2`. */
size_t blst_avx_p2_sizeof(void);
/** Returns the required ABI alignment of `blst_avx_p2`. */
size_t blst_avx_p2_alignof(void);
/** Returns the ABI size in bytes of `blst_avx_p2_affine`. */
size_t blst_avx_p2_affine_sizeof(void);
/** Returns the required ABI alignment of `blst_avx_p2_affine`. */
size_t blst_avx_p2_affine_alignof(void);
/** Returns the ABI size in bytes of `blst_avx_fp12`. */
size_t blst_avx_fp12_sizeof(void);
/** Returns the required ABI alignment of `blst_avx_fp12`. */
size_t blst_avx_fp12_alignof(void);

/**
 * Converts a scalar/public G1 projective point into the AVX G1 layout.
 *
 * Parameters:
 * - `out` receives the AVX-form point.
 * - `in` supplies the scalar/public point.
 *
 * Return value:
 * - `true` on success.
 * - `false` when AVX is unavailable or either pointer is null.
 */
bool blst_avx_p1_from_scalar(blst_avx_p1* out, const blst_p1* in);

/**
 * Converts an AVX G1 projective point into the scalar/public layout.
 *
 * Parameters:
 * - `out` receives the scalar/public point.
 * - `in` supplies the AVX-form point.
 *
 * Return value:
 * - `true` on success.
 * - `false` when AVX is unavailable or either pointer is null.
 */
bool blst_p1_from_avx(blst_p1* out, const blst_avx_p1* in);

/** Converts a scalar/public G1 affine point into the AVX affine layout. */
bool blst_avx_p1_affine_from_scalar(blst_avx_p1_affine* out, const blst_p1_affine* in);
/** Converts an AVX G1 affine point into the scalar/public affine layout. */
bool blst_p1_affine_from_avx(blst_p1_affine* out, const blst_avx_p1_affine* in);
/** Adds two AVX G1 projective points and writes the sum to `out`. */
bool blst_avx_p1_add(blst_avx_p1* out, const blst_avx_p1* a, const blst_avx_p1* b);

/**
 * Multiplies an AVX G1 projective point by an unsigned little-endian scalar.
 *
 * Parameters:
 * - `out` receives the product.
 * - `p` is the input point.
 * - `scalar` points to the scalar bytes.
 * - `nbits` is the number of meaningful low bits in `scalar`.
 */
bool blst_avx_p1_mult(blst_avx_p1* out, const blst_avx_p1* p, const byte* scalar, size_t nbits);
/** Converts an AVX G1 projective point to AVX affine form. */
bool blst_avx_p1_to_affine(blst_avx_p1_affine* out, const blst_avx_p1* in);
/** Converts an AVX G1 affine point to AVX projective form. */
bool blst_avx_p1_from_affine(blst_avx_p1* out, const blst_avx_p1_affine* in);
/** Returns `true` when the AVX G1 point is in the correct subgroup. */
bool blst_avx_p1_in_g1(const blst_avx_p1* p);
/** Returns `true` when the two AVX G1 points represent the same point. */
bool blst_avx_p1_is_equal(const blst_avx_p1* a, const blst_avx_p1* b);

/**
 * Returns the AVX G1 generator.
 *
 * Return value:
 * - Pointer to immutable static storage.
 * - Null only when AVX is unavailable.
 */
const blst_avx_p1* blst_avx_p1_generator(void);

/** Converts a scalar/public G2 projective point into the AVX G2 layout. */
bool blst_avx_p2_from_scalar(blst_avx_p2* out, const blst_p2* in);
/** Converts an AVX G2 projective point into the scalar/public layout. */
bool blst_p2_from_avx(blst_p2* out, const blst_avx_p2* in);
/** Converts a scalar/public G2 affine point into the AVX affine layout. */
bool blst_avx_p2_affine_from_scalar(blst_avx_p2_affine* out, const blst_p2_affine* in);
/** Converts an AVX G2 affine point into the scalar/public affine layout. */
bool blst_p2_affine_from_avx(blst_p2_affine* out, const blst_avx_p2_affine* in);
/** Adds two AVX G2 projective points and writes the sum to `out`. */
bool blst_avx_p2_add(blst_avx_p2* out, const blst_avx_p2* a, const blst_avx_p2* b);

/**
 * Multiplies an AVX G2 projective point by an unsigned little-endian scalar.
 *
 * Parameters mirror `blst_avx_p1_mult`.
 */
bool blst_avx_p2_mult(blst_avx_p2* out, const blst_avx_p2* p, const byte* scalar, size_t nbits);
/** Converts an AVX G2 projective point to AVX affine form. */
bool blst_avx_p2_to_affine(blst_avx_p2_affine* out, const blst_avx_p2* in);
/** Converts an AVX G2 affine point to AVX projective form. */
bool blst_avx_p2_from_affine(blst_avx_p2* out, const blst_avx_p2_affine* in);
/** Returns `true` when the AVX G2 point is in the correct subgroup. */
bool blst_avx_p2_in_g2(const blst_avx_p2* p);
/** Returns `true` when the two AVX G2 points represent the same point. */
bool blst_avx_p2_is_equal(const blst_avx_p2* a, const blst_avx_p2* b);

/**
 * Returns the AVX G2 generator.
 *
 * Return value:
 * - Pointer to immutable static storage.
 * - Null only when AVX is unavailable.
 */
const blst_avx_p2* blst_avx_p2_generator(void);

/** Converts a scalar/public Fp12 value into the AVX GT/Fp12 layout. */
bool blst_avx_fp12_from_scalar(blst_avx_fp12* out, const blst_fp12* in);
/** Converts an AVX GT/Fp12 value into the scalar/public layout. */
bool blst_fp12_from_avx(blst_fp12* out, const blst_avx_fp12* in);
/** Multiplies two AVX GT/Fp12 values and writes the product to `ret`. */
bool blst_avx_fp12_mul(blst_avx_fp12* ret, const blst_avx_fp12* a, const blst_avx_fp12* b);
/** Squares an AVX GT/Fp12 value and writes the result to `ret`. */
bool blst_avx_fp12_sqr(blst_avx_fp12* ret, const blst_avx_fp12* a);

/**
 * Applies the final exponentiation to an AVX Miller-loop value.
 *
 * Return value:
 * - `true` on success.
 * - `false` when AVX is unavailable or either pointer is null.
 */
bool blst_avx_fp12_final_exp(blst_avx_fp12* ret, const blst_avx_fp12* a);
/** Returns `true` when the two AVX GT/Fp12 values are equal. */
bool blst_avx_fp12_is_equal(const blst_avx_fp12* a, const blst_avx_fp12* b);
/** Returns `true` when the AVX GT/Fp12 value is in the target subgroup. */
bool blst_avx_fp12_in_group(const blst_avx_fp12* a);

/**
 * Returns the multiplicative identity in AVX GT/Fp12 form.
 *
 * Return value:
 * - Pointer to immutable static storage.
 * - Null only when AVX is unavailable.
 */
const blst_avx_fp12* blst_avx_fp12_one(void);

#ifdef __cplusplus
}
#endif

#endif
