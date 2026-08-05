/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BLST_AVX_H
#define BLST_AVX_H

#include "blst_types.h"

#if defined(__cplusplus) && defined(BLST_AVX_USE_CPP_TYPES) && \
    (defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE))
#include "blst/avx/curve/GT.h"
#include "blst/avx/curve/point.h"
#endif

#define BLST_AVX_ALIGNMENT 64
#define BLST_AVX_P1_SIZE 384
#define BLST_AVX_P1_AFFINE_SIZE 256
#define BLST_AVX_P2_SIZE 768
#define BLST_AVX_P2_AFFINE_SIZE 512
#define BLST_AVX_FP12_SIZE 1536

/**
 * AVX-native POD shells.
 *
 * Invariants:
 * - These layouts are stable for the current expert AVX ABI.
 * - The byte arrays are intentionally opaque to C and Rust callers.
 * - Callers may allocate storage of this size and alignment and pass it back
 *   to the library, but must not inspect or mutate the internal
 *   representation except through these entry points.
 * - C callers see opaque aligned byte storage. C++ AVX builds see the real
 *   internal AVX types, which avoids violating strict-aliasing/effective-type
 *   rules inside the implementation while preserving the C ABI layout.
 * - The storage is aligned for native AVX builds; portable AVX builds may use
 *   a weaker internal alignment while remaining representable by these shells.
 */
#if defined(__cplusplus) && defined(BLST_AVX_USE_CPP_TYPES) && \
    (defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE))
typedef blst::avx::G1Point blst_avx_p1;
typedef blst::avx::G1AffinePoint blst_avx_p1_affine;
typedef blst::avx::G2Point blst_avx_p2;
typedef blst::avx::G2AffinePoint blst_avx_p2_affine;
typedef blst::avx::GT blst_avx_fp12;
#else
typedef struct BLST_AUX_ALIGNAS(BLST_AVX_ALIGNMENT) { blst_byte b[BLST_AVX_P1_SIZE]; } blst_avx_p1;
typedef struct BLST_AUX_ALIGNAS(BLST_AVX_ALIGNMENT) { blst_byte b[BLST_AVX_P1_AFFINE_SIZE]; } blst_avx_p1_affine;
typedef struct BLST_AUX_ALIGNAS(BLST_AVX_ALIGNMENT) { blst_byte b[BLST_AVX_P2_SIZE]; } blst_avx_p2;
typedef struct BLST_AUX_ALIGNAS(BLST_AVX_ALIGNMENT) { blst_byte b[BLST_AVX_P2_AFFINE_SIZE]; } blst_avx_p2_affine;
typedef struct BLST_AUX_ALIGNAS(BLST_AVX_ALIGNMENT) { blst_byte b[BLST_AVX_FP12_SIZE]; } blst_avx_fp12;
#endif

/** Semantic alias for AVX Fp12 storage known to contain a GT element. */
typedef blst_avx_fp12 blst_avx_gt;

#ifdef __cplusplus
extern "C" {
#endif

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
 * Multiplies an AVX G1 projective subgroup point by an unsigned
 * little-endian scalar.
 *
 * `p` must already be on-curve and in G1; Release does not validate this
 * precondition. The scalar is interpreted modulo r. For canonical scalars,
 * runtime and memory access are independent of scalar contents for fixed
 * public `nbits`; noncanonical timing may reveal that reduction was required.
 *
 * Parameters:
 * - `out` receives the product.
 * - `p` is the input point.
 * - `scalar` points to the scalar bytes.
 * - `nbits` is the number of meaningful low bits in `scalar` and must not
 *   exceed 256. `scalar` may be null only when `nbits` is zero.
 */
bool blst_avx_p1_mult(blst_avx_p1* out, const blst_avx_p1* p, const blst_byte* scalar, size_t nbits);

/** Returns aggregate scratch bytes required by subgroup `blst_avx_p1_mult_pippenger`. */
size_t blst_avx_p1_mult_pippenger_scratch_sizeof(size_t npoints);
/** Returns byte-flag scratch bytes required by `blst_avx_p1_mult_pippenger`. */
size_t blst_avx_p1_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints);
/** Returns split-point scratch bytes required by allocation-free G1 GLV MSM. */
size_t blst_avx_p1_mult_pippenger_glv_points_scratch_sizeof(size_t npoints);
/** Returns split-scalar scratch bytes required by allocation-free G1 GLV MSM. */
size_t blst_avx_p1_mult_pippenger_glv_scalars_scratch_sizeof(size_t npoints);

/**
 * Variable-time subgroup-optimized Pippenger MSM over AVX G1 affine points.
 *
 * Every point must already be on-curve and in G1; Release does not validate
 * this precondition. Scalars contain at most 256 meaningful little-endian bits
 * each and are interpreted modulo r. The selected algorithm may use GLV.
 *
 * `scalars` is a contiguous array of `npoints` little-endian scalars, each
 * `(nbits+7)/8` bytes. `scratch` and `bucket_used` must point to caller-owned
 * storage of at least the sizes returned by the two scratch helpers above.
 */
bool blst_avx_p1_mult_pippenger(
    blst_avx_p1* out,
    const blst_avx_p1_affine* points,
    const blst_byte* scalars,
    size_t npoints,
    size_t nbits,
    void* scratch,
    size_t scratch_bytes,
    blst_byte* bucket_used,
    size_t bucket_used_bytes);

/**
 * Variable-time GLV Pippenger MSM over AVX G1 affine points with all scratch
 * supplied by the caller.
 *
 * Every point must already be on-curve and in G1; Release does not validate
 * this precondition. Each complete 256-bit scalar is reduced modulo r before
 * decomposition.
 *
 * `scalars` is a contiguous array of `npoints` 32-byte little-endian scalars.
 * `scratch`, `bucket_used`, `split_points`, and `split_scalars` must point to
 * caller-owned storage of at least the sizes returned by the scratch helpers.
 */
bool blst_avx_p1_mult_pippenger_glv_with_scratch(
    blst_avx_p1* out,
    const blst_avx_p1_affine* points,
    const blst_byte* scalars,
    size_t npoints,
    void* scratch,
    size_t scratch_bytes,
    blst_byte* bucket_used,
    size_t bucket_used_bytes,
    blst_avx_p1_affine* split_points,
    size_t split_points_bytes,
    blst_byte* split_scalars,
    size_t split_scalars_bytes);

/** Returns bucket scratch bytes needed by explicit-window AVX G1 GLV MSM. */
size_t blst_avx_p1_mult_pippenger_glv_window_scratch_sizeof(size_t npoints, size_t window);
/** Returns byte-flag scratch bytes needed by explicit-window AVX G1 GLV MSM. */
size_t blst_avx_p1_mult_pippenger_glv_window_bucket_used_scratch_sizeof(size_t npoints, size_t window);
/**
 * Variable-time GLV Pippenger MSM over AVX G1 affine points using an explicit
 * window and caller-owned scratch. Intended for benchmarking/tuning small MSMs.
 * Every point must already be on-curve and in G1; each 256-bit scalar is
 * reduced modulo r before GLV decomposition. Release does not validate point
 * membership.
 *
 * Current BTX-DKG measurements did not justify using this from production; the
 * auto-window GLV API remains the production path.
 */
bool blst_avx_p1_mult_pippenger_glv_window_with_scratch(
    blst_avx_p1* out,
    const blst_avx_p1_affine* points,
    const blst_byte* scalars,
    size_t npoints,
    size_t window,
    void* scratch,
    size_t scratch_bytes,
    blst_byte* bucket_used,
    size_t bucket_used_bytes,
    blst_avx_p1_affine* split_points,
    size_t split_points_bytes,
    blst_byte* split_scalars,
    size_t split_scalars_bytes);

/** Converts an AVX G1 projective point to AVX affine form. */
bool blst_avx_p1_to_affine(blst_avx_p1_affine* out, const blst_avx_p1* in);
/** Converts an AVX G1 affine point to AVX projective form. */
bool blst_avx_p1_from_affine(blst_avx_p1* out, const blst_avx_p1_affine* in);
/** Returns `true` when the AVX G1 point is on-curve and in the correct subgroup. */
bool blst_avx_p1_in_g1(const blst_avx_p1* p);
/** Returns `true` when the AVX G1 point is the point at infinity. */
bool blst_avx_p1_is_inf(const blst_avx_p1* p);
/** Returns `true` when the two AVX G1 points represent the same point. */
bool blst_avx_p1_is_equal(const blst_avx_p1* a, const blst_avx_p1* b);
/** Conditionally negates an AVX G1 projective point in place. */
bool blst_avx_p1_cneg(blst_avx_p1* p, bool cbit);
/**
 * Runs `width` rounds of the {-1,0,+1} random-linear-combination subgroup
 * check over AVX G1 projective points.
 *
 * `points` is an array of `count` point pointers. `coefficients` is row-major
 * with `width * count` bytes; coefficient values are 0, 1, or 2, where 2
 * means -1. `shift` is the subgroup random-shift element used to initialize
 * each round sum.
 */
bool blst_avx_p1_batch_subgroup_check(
    const blst_avx_p1* const* points,
    size_t count,
    const blst_byte* coefficients,
    size_t width,
    const blst_avx_p1* shift);
/** Compresses an AVX G1 projective point to the standard 48-byte encoding. */
bool blst_avx_p1_compress(blst_byte out[48], const blst_avx_p1* in);
/** Serializes an AVX G1 projective point to the standard 96-byte affine encoding. */
bool blst_avx_p1_serialize(blst_byte out[96], const blst_avx_p1* in);
/** Parses a standard 96-byte uncompressed G1 encoding into AVX projective form. */
BLST_ERROR blst_avx_p1_deserialize(blst_avx_p1* out, const blst_byte in[96]);
/** Parses a standard 48-byte compressed G1 encoding into AVX projective form. */
BLST_ERROR blst_avx_p1_uncompress(blst_avx_p1* out, const blst_byte in[48]);
/** Returns caller-owned scratch bytes needed by subgroup AVX G1 Pippenger MSM. */
size_t blst_avx_p1s_mult_pippenger_scratch_sizeof(size_t npoints);
/**
 * Computes variable-time subgroup G1 MSM.
 *
 * Inputs must already be on-curve G1 elements. Scalars contain at most 256
 * meaningful little-endian bits and are interpreted modulo r. Release checks
 * pointer/bit-count shape but does not validate point membership. The checked
 * one-point path requires no scratch.
 */
bool blst_avx_p1s_mult_pippenger(
    blst_avx_p1* out,
    const blst_avx_p1* points,
    size_t npoints,
    const blst_byte* scalars,
    size_t nbits,
    blst_limb_t* scratch);
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
/** Conditionally negates an AVX G2 projective point in place. */
bool blst_avx_p2_cneg(blst_avx_p2* p, bool cbit);
/**
 * Runs `width` rounds of the {-1,0,+1} random-linear-combination subgroup
 * check over AVX G2 projective points.
 *
 * Layout and coefficient conventions match `blst_avx_p1_batch_subgroup_check`.
 */
bool blst_avx_p2_batch_subgroup_check(
    const blst_avx_p2* const* points,
    size_t count,
    const blst_byte* coefficients,
    size_t width,
    const blst_avx_p2* shift);

/**
 * Multiplies an AVX G2 projective subgroup point by an unsigned
 * little-endian scalar.
 *
 * Preconditions, modulo-r scalar semantics, fixed-`nbits` timing, and
 * parameters mirror `blst_avx_p1_mult`.
 */
bool blst_avx_p2_mult(blst_avx_p2* out, const blst_avx_p2* p, const blst_byte* scalar, size_t nbits);

/** Returns aggregate scratch bytes required by subgroup `blst_avx_p2_mult_pippenger`. */
size_t blst_avx_p2_mult_pippenger_scratch_sizeof(size_t npoints);
/** Returns byte-flag scratch bytes required by `blst_avx_p2_mult_pippenger`. */
size_t blst_avx_p2_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints);

/**
 * Variable-time subgroup-optimized Pippenger MSM over AVX G2 affine points.
 *
 * Every point must already be on-curve and in G2; Release does not validate
 * this precondition. Scalars contain at most 256 meaningful little-endian bits
 * each and are interpreted modulo r. The selected algorithm may use GLS.
 *
 * `scalars` is a contiguous array of `npoints` little-endian scalars, each
 * `(nbits+7)/8` bytes. `scratch` and `bucket_used` must point to caller-owned
 * storage of at least the sizes returned by the two scratch helpers above.
 */
bool blst_avx_p2_mult_pippenger(
    blst_avx_p2* out,
    const blst_avx_p2_affine* points,
    const blst_byte* scalars,
    size_t npoints,
    size_t nbits,
    void* scratch,
    size_t scratch_bytes,
    blst_byte* bucket_used,
    size_t bucket_used_bytes);

/** Converts an AVX G2 projective point to AVX affine form. */
bool blst_avx_p2_to_affine(blst_avx_p2_affine* out, const blst_avx_p2* in);
/** Converts an AVX G2 affine point to AVX projective form. */
bool blst_avx_p2_from_affine(blst_avx_p2* out, const blst_avx_p2_affine* in);
/** Returns `true` when the AVX G2 point is on-curve and in the correct subgroup. */
bool blst_avx_p2_in_g2(const blst_avx_p2* p);
/** Returns `true` when the AVX G2 point is the point at infinity. */
bool blst_avx_p2_is_inf(const blst_avx_p2* p);
/** Returns `true` when the two AVX G2 points represent the same point. */
bool blst_avx_p2_is_equal(const blst_avx_p2* a, const blst_avx_p2* b);
/** Compresses an AVX G2 projective point to the standard 96-byte encoding. */
bool blst_avx_p2_compress(blst_byte out[96], const blst_avx_p2* in);
/** Serializes an AVX G2 projective point to the standard 192-byte affine encoding. */
bool blst_avx_p2_serialize(blst_byte out[192], const blst_avx_p2* in);
/** Parses a standard 96-byte compressed G2 encoding into AVX projective form. */
BLST_ERROR blst_avx_p2_uncompress(blst_avx_p2* out, const blst_byte in[96]);
/** Parses a standard 192-byte uncompressed G2 encoding into AVX projective form. */
BLST_ERROR blst_avx_p2_deserialize(blst_avx_p2* out, const blst_byte in[192]);
/** Returns caller-owned scratch bytes needed by subgroup AVX G2 Pippenger MSM. */
size_t blst_avx_p2s_mult_pippenger_scratch_sizeof(size_t npoints);
/**
 * Computes variable-time subgroup G2 MSM.
 *
 * Inputs must already be on-curve G2 elements. Scalars contain at most 256
 * meaningful little-endian bits and are interpreted modulo r. Release checks
 * pointer/bit-count shape but does not validate point membership. The checked
 * one-point path requires no scratch.
 */
bool blst_avx_p2s_mult_pippenger(
    blst_avx_p2* out,
    const blst_avx_p2* points,
    size_t npoints,
    const blst_byte* scalars,
    size_t nbits,
    blst_limb_t* scratch);
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

/** Returns aggregate scratch bytes required by `blst_avx_gt_mult_pippenger`. */
size_t blst_avx_gt_mult_pippenger_scratch_sizeof(size_t npoints);
/** Returns byte-flag scratch bytes required by `blst_avx_gt_mult_pippenger`. */
size_t blst_avx_gt_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints);

/**
 * Variable-time GLV4 Pippenger MSM over previously validated GT elements.
 *
 * Each scalar is a 256-bit little-endian integer interpreted modulo r.
 * Release does not validate GT membership. Scratch sizes come from the two
 * matching helpers above.
 */
bool blst_avx_gt_mult_pippenger(
    blst_avx_gt* out,
    const blst_avx_gt* points,
    const blst_byte* scalars,
    size_t npoints,
    void* scratch,
    size_t scratch_bytes,
    blst_byte* bucket_used,
    size_t bucket_used_bytes);

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

/**
 * Computes final_exponentiation(miller_loop(g2, g1)) in AVX-native form.
 */
bool blst_avx_pairing(blst_avx_fp12* out, const blst_avx_p1* g1, const blst_avx_p2* g2);

#ifdef __cplusplus
}
#endif

#endif
