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

/** Immutable prepared-G2 batch. Contents are owned by CatBLST. */
typedef struct blst_avx_prepared_g2_batch blst_avx_prepared_g2_batch;
/** Mutable scratch for evaluating a prepared-G2 batch. */
typedef struct blst_avx_prepared_g2_batch_workspace blst_avx_prepared_g2_batch_workspace;

/**
 * Fills exactly `out_len` bytes from a cryptographically secure random source.
 * Returns false on entropy failure. The callback must not throw across the C
 * boundary; CatBLST translates a C++ exception to failure when exceptions are
 * enabled.
 */
typedef bool (*blst_avx_random_bytes_fn)(
    void* context, blst_byte* out, size_t out_len);

/**
 * Expert AVX C-boundary rules.
 *
 * - No exception may cross this ABI. These entry points allocate no memory;
 *   exceptions raised by the random callback are translated to failure when
 *   the library is built with exception support.
 * - A `false` or non-success result does not publish a partial output. Unless
 *   a declaration says otherwise, output storage is unchanged on failure.
 * - Same-type arithmetic supports exact output/input aliasing. Partial overlap
 *   is never supported. Differently typed conversion buffers must not overlap.
 * - Scratch-size helpers return zero for empty work, an unavailable backend,
 *   an invalid tuning parameter, or an unrepresentable size. Callers must
 *   distinguish empty work before treating zero as failure.
 * - The opaque byte representation is not a serialization format. Only
 *   library-produced values and exact copies satisfy the representation
 *   invariants required by the operations below.
 */

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
 *
 * `in` must be on-curve. The conversion preserves subgroup membership but
 * neither requires nor establishes it.
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
 *
 * `in` must be on-curve. The conversion preserves subgroup membership but
 * neither requires nor establishes it.
 */
bool blst_p1_from_avx(blst_p1* out, const blst_avx_p1* in);

/** Converts an on-curve scalar/public G1 affine point into the AVX layout.
 * Subgroup membership is preserved but neither required nor established. */
bool blst_avx_p1_affine_from_scalar(blst_avx_p1_affine* out, const blst_p1_affine* in);
/** Converts an on-curve AVX G1 affine point into the scalar/public layout.
 * Subgroup membership is preserved but neither required nor established. */
bool blst_p1_affine_from_avx(blst_p1_affine* out, const blst_avx_p1_affine* in);
/** Completely adds two on-curve AVX G1 projective points. Subgroup membership
 * is not required or validated. */
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

/**
 * Constant-time MSM over AVX G1 projective subgroup points and secret scalars.
 * Computes `out = sum_i scalars[i] * points[i]`.
 *
 * Every point must already be on-curve and in G1; Release does not validate
 * this precondition. `scalars` contains `npoints` consecutive reduced native
 * Montgomery-form `blst_fr` values. Runtime and memory access are independent
 * of scalar contents for fixed public `npoints`.
 *
 * The selected kernel shares its doubling schedule in fixed 32-, 16-, and
 * 8-point chunks and uses a scalar-multiplication tail for shorter remainders.
 * It allocates nothing. For 8 or more points, `scratch` must be aligned to
 * `BLST_AVX_ALIGNMENT` and provide at least the number of bytes returned by
 * `blst_avx_p1s_mult_ct_scratch_sizeof`; the used prefix must be disjoint
 * from both inputs and `out`. Empty input returns the identity. Null `points`
 * or `scalars` is accepted only for empty input. Failure leaves `out`
 * unchanged.
 */
size_t blst_avx_p1s_mult_ct_scratch_sizeof(void);
bool blst_avx_p1s_mult_ct(
    blst_avx_p1* out,
    const blst_avx_p1* points,
    const blst_fr* scalars,
    size_t npoints,
    void* scratch,
    size_t scratch_bytes);

/**
 * Returns caller-owned scratch capacity for public-scalar native G1 MSM.
 *
 * Storage sized for `max_points` supports every input length through that
 * capacity. The selected implementation uses shared-doubling GLV batches from
 * two points through the small-input range and Pippenger for larger inputs.
 * A one-point capacity needs no scratch. Zero otherwise means that the
 * requested capacity is not representable.
 */
size_t blst_avx_p1s_mult_fr_vartime_scratch_sizeof(size_t max_points);

/**
 * Variable-time MSM over AVX G1 projective subgroup points and public scalars.
 * Computes `out = sum_i scalars[i] * points[i]`.
 *
 * Every point must already be on-curve and in G1; Release does not validate
 * this precondition. `scalars` contains `npoints` consecutive reduced native
 * Montgomery-form `blst_fr` values and must contain no secrets. The selected
 * implementation uses direct scalar-dependent table accesses and branches.
 *
 * The function allocates nothing. Nonzero scratch requirements must be
 * obtained from `blst_avx_p1s_mult_fr_vartime_scratch_sizeof`; `scratch` must
 * be aligned to `BLST_AVX_ALIGNMENT`, and its used prefix must be disjoint
 * from both inputs and `out`. Empty input returns the identity. Null `points`
 * or `scalars` is accepted only for empty input. Failure leaves `out`
 * unchanged.
 */
bool blst_avx_p1s_mult_fr_vartime(
    blst_avx_p1* out,
    const blst_avx_p1* points,
    const blst_fr* scalars,
    size_t npoints,
    void* scratch,
    size_t scratch_bytes);

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
 * For nonempty input, the output, point array, scalar array, required scratch
 * prefix, and required bucket-flag prefix must be pairwise disjoint. Invalid
 * shape, size, alignment, or overlap returns false without changing `out`.
 * The function allocates nothing and does not throw across the C boundary.
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
 * The output, both inputs, and all four required scratch prefixes must be
 * pairwise disjoint. Invalid shape, size, alignment, or overlap returns false
 * without changing `out`. The function allocates nothing and does not throw.
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
 * The disjoint live-range, failure-atomicity, allocation, and exception
 * contract is the same as `blst_avx_p1_mult_pippenger_glv_with_scratch`.
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

/** Converts an on-curve AVX G1 projective point to affine form. Subgroup
 * membership is preserved but neither required nor established. */
bool blst_avx_p1_to_affine(blst_avx_p1_affine* out, const blst_avx_p1* in);
/** Converts an on-curve AVX G1 affine point to projective form. Subgroup
 * membership is preserved but neither required nor established. */
bool blst_avx_p1_from_affine(blst_avx_p1* out, const blst_avx_p1_affine* in);
/** Returns `true` exactly when the candidate is on-curve and in G1. No prior
 * point-domain precondition applies. */
bool blst_avx_p1_in_g1(const blst_avx_p1* p);
/** Tests whether an on-curve AVX G1 point is infinity. This is not a validator. */
bool blst_avx_p1_is_inf(const blst_avx_p1* p);
/** Compares two on-curve AVX G1 points. This is not a validator. */
bool blst_avx_p1_is_equal(const blst_avx_p1* a, const blst_avx_p1* b);
/** Conditionally negates an on-curve AVX G1 point. Subgroup membership is not
 * required or validated. */
bool blst_avx_p1_cneg(blst_avx_p1* p, bool cbit);
/**
 * Runs `width` rounds of the {-1,0,+1} random-linear-combination subgroup
 * check over AVX G1 projective points.
 *
 * This is an expert primitive, not a misuse-resistant validation boundary.
 * Candidate points must be on-curve and `shift` must already be in G1. A
 * nonempty batch requires `width > 0`. Coefficient rows must be generated by a
 * sound unpredictable challenge process and must not be all zero; the function
 * cannot validate those conditions. Violating them can make invalid points pass.
 * General callers should use `blst_avx_p1s_in_g1`, which owns those choices.
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
/**
 * Tests whether every point belongs to G1 using an internally constructed
 * randomized batch check targeting a false-accept probability below 2^-128.
 *
 * `points` contains `count` on-curve or adversarial AVX projective candidates.
 * For a nonempty batch, `fill_random` must provide one fresh, unpredictable
 * 32-byte seed after the batch has been fixed, and the pointed-to values must
 * remain unchanged for the duration of the call. CatBLST derives independent
 * uniform ternary rows with domain-separated SHA-256 and owns the challenge
 * width and subgroup shift. The stated bound assumes SHA-256 behaves as a PRF
 * for that fresh seed and uses the fact that the selected curve has no
 * nontrivial 2-torsion. The routine allocates no memory. With AVX support,
 * empty input succeeds without invoking `fill_random`. Any unavailable AVX
 * backend, invalid shape, off-curve point, entropy failure, or failed subgroup
 * relation returns false.
 */
bool blst_avx_p1s_in_g1(
    const blst_avx_p1* const* points,
    size_t count,
    blst_avx_random_bytes_fn fill_random,
    void* random_context);
/** Compresses an on-curve AVX G1 projective point to 48 bytes. Subgroup
 * membership is not required or validated. */
bool blst_avx_p1_compress(blst_byte out[48], const blst_avx_p1* in);
/** Serializes an on-curve AVX G1 projective point to 96 bytes. Subgroup
 * membership is not required or validated. */
bool blst_avx_p1_serialize(blst_byte out[96], const blst_avx_p1* in);
/** Decodes either a 48-byte compressed or 96-byte uncompressed G1 encoding,
 * selected by its flags; the caller must provide the corresponding readable
 * size. Success establishes canonical encoding and on-curve status, not G1
 * membership. Infinity is accepted; on failure `out` is unchanged. */
BLST_ERROR blst_avx_p1_deserialize(blst_avx_p1* out, const blst_byte in[96]);
/** Decodes a 48-byte compressed G1 encoding. Success establishes canonical
 * encoding and on-curve status, not G1 membership. Infinity is accepted; on
 * failure `out` is unchanged. */
BLST_ERROR blst_avx_p1_uncompress(blst_avx_p1* out, const blst_byte in[48]);
/** Returns caller-owned scratch bytes needed by subgroup AVX G1 Pippenger MSM. */
size_t blst_avx_p1s_mult_pippenger_scratch_sizeof(size_t npoints);
/**
 * Computes variable-time subgroup G1 MSM.
 *
 * Inputs must already be prime-order G1 subgroup elements. Scalars contain at most 256
 * meaningful little-endian bits and are interpreted modulo r. Release checks
 * pointer/bit-count shape but does not validate point membership. The checked
 * one-point path requires no scratch and permits `out == points`; other input
 * and output ranges must be valid. For multiple points, output, points,
 * scalars, and the queried scratch range must be pairwise disjoint. Invalid
 * shape or overlap returns false without changing `out`. The function is
 * allocation-free and nonthrowing.
 */
bool blst_avx_p1s_mult_pippenger(
    blst_avx_p1* out,
    const blst_avx_p1* points,
    size_t npoints,
    const blst_byte* scalars,
    size_t nbits,
    blst_limb_t* scratch);
/**
 * Writes the canonical AVX G1 identity to caller-owned storage.
 *
 * `out` must be non-null and suitably aligned. The operation allocates
 * nothing and performs no scalar multiplication.
 */
bool blst_avx_p1_identity(blst_avx_p1* out);
/**
 * Returns the AVX G1 generator.
 *
 * Return value:
 * - Pointer to immutable static storage.
 * - Null only when AVX is unavailable.
 */
const blst_avx_p1* blst_avx_p1_generator(void);

/** Converts an on-curve scalar/public G2 projective point into the AVX layout.
 * Subgroup membership is preserved but neither required nor established. */
bool blst_avx_p2_from_scalar(blst_avx_p2* out, const blst_p2* in);
/** Converts an on-curve AVX G2 projective point into the scalar/public layout.
 * Subgroup membership is preserved but neither required nor established. */
bool blst_p2_from_avx(blst_p2* out, const blst_avx_p2* in);
/** Converts an on-curve scalar/public G2 affine point into the AVX layout.
 * Subgroup membership is preserved but neither required nor established. */
bool blst_avx_p2_affine_from_scalar(blst_avx_p2_affine* out, const blst_p2_affine* in);
/** Converts an on-curve AVX G2 affine point into the scalar/public layout.
 * Subgroup membership is preserved but neither required nor established. */
bool blst_p2_affine_from_avx(blst_p2_affine* out, const blst_avx_p2_affine* in);
/** Completely adds two on-curve AVX G2 projective points. Subgroup membership
 * is not required or validated. */
bool blst_avx_p2_add(blst_avx_p2* out, const blst_avx_p2* a, const blst_avx_p2* b);
/** Conditionally negates an on-curve AVX G2 point. Subgroup membership is not
 * required or validated. */
bool blst_avx_p2_cneg(blst_avx_p2* p, bool cbit);
/**
 * Runs `width` rounds of the {-1,0,+1} random-linear-combination subgroup
 * check over AVX G2 projective points.
 *
 * This has the same expert-only preconditions and challenge-soundness warning
 * as `blst_avx_p1_batch_subgroup_check`, with `shift` required to be in G2.
 * General callers should use `blst_avx_p2s_in_g2`.
 */
bool blst_avx_p2_batch_subgroup_check(
    const blst_avx_p2* const* points,
    size_t count,
    const blst_byte* coefficients,
    size_t width,
    const blst_avx_p2* shift);
/** Returns the caller-owned workspace size required by `blst_avx_p2s_in_g2`. */
size_t blst_avx_p2s_in_g2_workspace_sizeof(void);

/**
 * Tests whether every point belongs to the prime-order BLS12-381 G2 subgroup.
 *
 * For a nonempty batch, `fill_random` must provide one fresh, unpredictable
 * 32-byte seed after the batch has been fixed. The pointed-to values must
 * remain unchanged for the duration of the call. CatBLST derives 22 rows whose
 * coefficients are uniform in {-6,...,+6}. It checks one random linear
 * combination per row with a fixed G2 shift.
 *
 * Because 13 is the smallest prime divisor of the BLS12-381 G2 cofactor, a
 * fixed invalid batch passes one ideal row with probability at most 1/13.
 * Under the SHA-256 PRF assumption, one invocation therefore has false-accept
 * probability at most 13^-22 plus the PRF distinguishing advantage. The
 * concrete term 13^-22 is less than 2^-81.4.
 *
 * `workspace` must provide at least the queried byte count, must not overlap
 * any candidate point, and is modified. The function allocates no memory.
 * Empty input succeeds without entropy or workspace. An unavailable AVX
 * backend, invalid shape, off-curve point, entropy failure, insufficient
 * workspace, or failed subgroup relation returns false.
 */
bool blst_avx_p2s_in_g2(
    const blst_avx_p2* const* points,
    size_t count,
    blst_avx_random_bytes_fn fill_random,
    void* random_context,
    void* workspace,
    size_t workspace_bytes);

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
 * For nonempty input, the output, point array, scalar array, required scratch
 * prefix, and required bucket-flag prefix must be pairwise disjoint. Invalid
 * shape, size, alignment, or overlap returns false without changing `out`.
 * The function allocates nothing and does not throw across the C boundary.
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

/** Converts an on-curve AVX G2 projective point to affine form. Subgroup
 * membership is preserved but neither required nor established. */
bool blst_avx_p2_to_affine(blst_avx_p2_affine* out, const blst_avx_p2* in);
/** Converts an on-curve AVX G2 affine point to projective form. Subgroup
 * membership is preserved but neither required nor established. */
bool blst_avx_p2_from_affine(blst_avx_p2* out, const blst_avx_p2_affine* in);
/** Returns `true` exactly when the candidate is on-curve and in G2. No prior
 * point-domain precondition applies. */
bool blst_avx_p2_in_g2(const blst_avx_p2* p);
/** Tests whether an on-curve AVX G2 point is infinity. This is not a validator. */
bool blst_avx_p2_is_inf(const blst_avx_p2* p);
/** Compares two on-curve AVX G2 points. This is not a validator. */
bool blst_avx_p2_is_equal(const blst_avx_p2* a, const blst_avx_p2* b);
/** Compresses an on-curve AVX G2 projective point to 96 bytes. Subgroup
 * membership is not required or validated. */
bool blst_avx_p2_compress(blst_byte out[96], const blst_avx_p2* in);
/** Serializes an on-curve AVX G2 projective point to 192 bytes. Subgroup
 * membership is not required or validated. */
bool blst_avx_p2_serialize(blst_byte out[192], const blst_avx_p2* in);
/** Decodes a 96-byte compressed G2 encoding. Success establishes canonical
 * encoding and on-curve status, not G2 membership. Infinity is accepted; on
 * failure `out` is unchanged. */
BLST_ERROR blst_avx_p2_uncompress(blst_avx_p2* out, const blst_byte in[96]);
/** Decodes either a 96-byte compressed or 192-byte uncompressed G2 encoding,
 * selected by its flags; the caller must provide the corresponding readable
 * size. Success establishes canonical encoding and on-curve status, not G2
 * membership. Infinity is accepted; on failure `out` is unchanged. */
BLST_ERROR blst_avx_p2_deserialize(blst_avx_p2* out, const blst_byte in[192]);
/** Returns caller-owned scratch bytes needed by subgroup AVX G2 Pippenger MSM. */
size_t blst_avx_p2s_mult_pippenger_scratch_sizeof(size_t npoints);
/**
 * Computes variable-time subgroup G2 MSM.
 *
 * Inputs must already be prime-order G2 subgroup elements. Scalars contain at most 256
 * meaningful little-endian bits and are interpreted modulo r. Release checks
 * pointer/bit-count shape but does not validate point membership. The checked
 * one-point path requires no scratch and permits `out == points`; other input
 * and output ranges must be valid. For multiple points, output, points,
 * scalars, and the queried scratch range must be pairwise disjoint. Invalid
 * shape or overlap returns false without changing `out`. The function is
 * allocation-free and nonthrowing.
 */
bool blst_avx_p2s_mult_pippenger(
    blst_avx_p2* out,
    const blst_avx_p2* points,
    size_t npoints,
    const blst_byte* scalars,
    size_t nbits,
    blst_limb_t* scratch);
/**
 * Writes the canonical AVX G2 identity to caller-owned storage.
 *
 * `out` must be non-null and suitably aligned. The operation allocates
 * nothing and performs no scalar multiplication.
 */
bool blst_avx_p2_identity(blst_avx_p2* out);
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
/**
 * Serializes a validated AVX GT value as twelve canonical 48-byte big-endian
 * Fp coefficients. The coefficient order is Fp2 index, then Fp6 half, then Fp
 * component, exactly matching `blst_bendian_from_fp12`. `out` and `in` must be
 * non-null. No allocation is performed.
 */
bool blst_avx_gt_serialize(blst_byte out[576], const blst_avx_gt* in);
/**
 * Decodes the same twelve-coefficient order and requires every coefficient to
 * be canonical and the assembled Fp12 value to lie in GT. The identity is
 * accepted. Failure returns false and leaves `out` unchanged. Input bytes may
 * overlap `out`; both pointers must be non-null. Membership validation is
 * variable-time with public input. No allocation is performed.
 */
bool blst_avx_gt_deserialize(blst_avx_gt* out, const blst_byte in[576]);
/**
 * Constant-time scalar multiplication of a validated AVX GT element.
 *
 * `scalar` names `ceil(nbits/8)` readable little-endian bytes and may be null
 * only when `nbits` is zero. `nbits` must not exceed 256. Runtime and memory
 * access are independent of a canonical scalar's contents for fixed public
 * `nbits`; noncanonical scalars are accepted modulo r, but reduction timing
 * may reveal that reduction was required. `ret` may exactly alias `a`; partial
 * overlap is unsupported. The function allocates nothing.
 */
bool blst_avx_gt_mult(
    blst_avx_gt* ret,
    const blst_avx_gt* a,
    const blst_byte* scalar,
    size_t nbits);
/**
 * Multiplies two well-formed AVX Fp12 values and writes the generic field
 * product to `ret`. GT membership is neither required nor established.
 * `ret` may exactly alias either input; partial overlap is not supported.
 */
bool blst_avx_fp12_mul(blst_avx_fp12* ret, const blst_avx_fp12* a, const blst_avx_fp12* b);
/**
 * Squares a well-formed AVX Fp12 value with the dedicated generic field-square
 * kernel. GT membership is neither required nor established. `ret` may
 * exactly alias `a`; partial overlap is not supported.
 */
bool blst_avx_fp12_sqr(blst_avx_fp12* ret, const blst_avx_fp12* a);
/** Inverts a validated AVX target-group value. `ret` may alias `a`. */
bool blst_avx_gt_inverse(blst_avx_gt* ret, const blst_avx_gt* a);

/** Returns aggregate scratch bytes required by `blst_avx_gt_mult_pippenger`. */
size_t blst_avx_gt_mult_pippenger_scratch_sizeof(size_t npoints);
/** Returns byte-flag scratch bytes required by `blst_avx_gt_mult_pippenger`. */
size_t blst_avx_gt_mult_pippenger_bucket_used_scratch_sizeof(size_t npoints);

/**
 * Variable-time GLV4 Pippenger MSM over previously validated GT elements.
 *
 * Each scalar is a 256-bit little-endian integer interpreted modulo r.
 * Release does not validate GT membership. Scratch sizes come from the two
 * matching helpers above. For nonempty input, the output, point array, scalar
 * array, required aggregate-scratch prefix, and required bucket-flag prefix
 * must be pairwise disjoint. Invalid shape, size, alignment, or overlap returns
 * false without changing `out`. The function allocates nothing and does not
 * throw across the C boundary.
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
 * Preconditions:
 * - `a` represents a valid nonzero Fp12 value, normally produced by a Miller
 *   loop. No field, nonzero, or subgroup validation is performed.
 * - `ret` and `a` may be identical but must not partially overlap.
 *
 * Postcondition:
 * - On success, `ret` contains BLST's cubed final exponent in GT.
 * - Failure leaves `ret` unspecified.
 *
 * Return value:
 * - `true` on success.
 * - `false` only when AVX is unavailable or either pointer is null; this does
 *   not report semantic input validation.
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
 * Computes one unfinalized Miller-loop contribution in AVX-native form.
 *
 * Both inputs must be library-produced prime-order subgroup points; no curve
 * or subgroup validation is performed. Identity inputs produce the Miller
 * multiplicative identity. The operation allocates no memory. `out` may not
 * alias either input and is unspecified on failure.
 *
 * Return value:
 * - `true` on success.
 * - `false` only when AVX is unavailable or an argument is null; this does
 *   not report semantic input validation.
 */
bool blst_avx_miller_loop(
    blst_avx_fp12* out,
    const blst_avx_p1* g1,
    const blst_avx_p2* g2);

/**
 * Computes final_exponentiation(miller_loop(g2, g1)) in AVX-native form.
 * Both inputs must already be prime-order subgroup points; no point validation
 * is performed and behavior for non-subgroup inputs is not promised. The
 * boolean reports only null arguments or unavailable AVX support; it does not
 * report curve or subgroup validation.
 */
bool blst_avx_pairing(blst_avx_fp12* out, const blst_avx_p1* g1, const blst_avx_p2* g2);

/**
 * Computes one final exponentiation of the product of `n` Miller loops over
 * contiguous native affine arrays. Empty input returns the GT identity. Points
 * must already be subgroup validated. No allocation is performed.
 */
bool blst_avx_pairing_n(
    blst_avx_fp12* out,
    const blst_avx_p1_affine* g1,
    const blst_avx_p2_affine* g2,
    size_t n);

/** Returns caller-owned storage bytes required for one native G2 line table. */
size_t blst_avx_prepared_g2_lines_sizeof(void);
/** Returns the required alignment for one native G2 line table. */
size_t blst_avx_prepared_g2_lines_alignof(void);

/**
 * Initializes one native G2 Miller-line table in caller-owned storage.
 *
 * `storage` must be aligned to `blst_avx_prepared_g2_lines_alignof()` and
 * contain at least `storage_size` bytes. Construction is intended for the
 * cold public-parameter setup path; evaluation performs no allocation.
 * `g2` must be a library-produced, subgroup-validated native point.
 */
bool blst_avx_prepared_g2_lines_init(
    void* storage,
    size_t storage_size,
    const blst_avx_p2* g2);

/**
 * Computes final_exponentiation(miller_loop(prepared_g2, g1)).
 * `storage` must have completed `blst_avx_prepared_g2_lines_init`
 * successfully, and `g1` must be a subgroup-validated native point.
 */
bool blst_avx_pairing_prepared_g2(
    blst_avx_fp12* out,
    const void* storage,
    size_t storage_size,
    const blst_avx_p1* g1);

/**
 * Prepares one Miller-line table for every point in `g2[0..count)`.
 *
 * Construction is a cold-path allocating operation. Inputs must be
 * library-produced, subgroup-validated native points; infinity is accepted
 * and is recorded as the pairing identity. The returned object is immutable
 * and may be shared across threads. Returns null for unavailable AVX support,
 * invalid input, an empty batch, or allocation failure.
 */
blst_avx_prepared_g2_batch* blst_avx_prepared_g2_batch_new(
    const blst_avx_p2* g2,
    size_t count);

/** Returns the logical number of prepared G2 operands, or zero for null. */
size_t blst_avx_prepared_g2_batch_len(const blst_avx_prepared_g2_batch* batch);

/**
 * Returns an implementation-accounting estimate of owned bytes, including
 * line tables and identity flags but excluding allocator bookkeeping.
 */
size_t blst_avx_prepared_g2_batch_storage_bytes(
    const blst_avx_prepared_g2_batch* batch);

/** Releases one immutable prepared batch. Null is accepted. */
void blst_avx_prepared_g2_batch_free(blst_avx_prepared_g2_batch* batch);

/**
 * Allocates evaluation scratch sized for the complete prepared batch.
 * The workspace may subsequently evaluate any prefix of that batch without
 * allocating. A workspace is mutable, is not thread-safe, and must not be used
 * concurrently. Returns null for invalid input or allocation failure.
 */
blst_avx_prepared_g2_batch_workspace*
blst_avx_prepared_g2_batch_workspace_new(
    const blst_avx_prepared_g2_batch* batch);

/**
 * Allocates evaluation scratch for at most `capacity` pairings.
 *
 * `capacity` must be nonzero and no greater than the immutable batch length.
 * This form lets independent range workers retain only the scratch required by
 * their maximum range grain. It has the same ownership and thread-safety rules
 * as the full-batch constructor.
 */
blst_avx_prepared_g2_batch_workspace*
blst_avx_prepared_g2_batch_workspace_new_with_capacity(
    const blst_avx_prepared_g2_batch* batch,
    size_t capacity);

/** Releases one prepared-batch workspace. Null is accepted. */
void blst_avx_prepared_g2_batch_workspace_free(
    blst_avx_prepared_g2_batch_workspace* workspace);

/**
 * Computes `out[i] = e(g1, prepared_g2[i])` for `0 <= i < count`.
 *
 * `count` may be any prefix length through the prepared batch length. The
 * output has exactly `count` native GT elements and must not overlap `g1` or
 * either owned object. The workspace must have capacity for `count`. This
 * operation converts `g1` to affine once, batches final exponentiation across
 * the prefix, and allocates no memory. Inputs must already satisfy the native
 * point invariants. On failure, output contents are unspecified.
 */
bool blst_avx_prepared_g2_batch_pairing_prefix(
    blst_avx_gt* out,
    size_t count,
    const blst_avx_p1* g1,
    const blst_avx_prepared_g2_batch* batch,
    blst_avx_prepared_g2_batch_workspace* workspace);

/**
 * Computes `out[i] = e(g1, prepared_g2[indices[i]])` for `0 <= i < count`.
 *
 * Every index must be less than the prepared batch length. Repetition is
 * allowed. The output has exactly `count` native GT elements. This is the
 * sparse counterpart to the prefix evaluator and has the same ownership,
 * validation, and allocation-free workspace requirements.
 */
bool blst_avx_prepared_g2_batch_pairing_indices(
    blst_avx_gt* out,
    const size_t* indices,
    size_t count,
    const blst_avx_p1* g1,
    const blst_avx_prepared_g2_batch* batch,
    blst_avx_prepared_g2_batch_workspace* workspace);

/**
 * Computes `out[i] = e(g1, prepared_g2[begin + i])` for `0 <= i < count`.
 *
 * The contiguous range must fit inside the immutable batch and `count` must not
 * exceed the mutable workspace capacity. The operation converts `g1` to affine
 * once, batches final exponentiation across the range, allocates no memory, and
 * retains no state. Concurrent calls may share the batch but require distinct
 * workspaces and disjoint output ranges. Nonempty output must not overlap the
 * G1 input, prepared batch, or workspace. Inputs must already satisfy the
 * native point invariants. On failure, output contents are unspecified.
 */
bool blst_avx_prepared_g2_batch_pairing_range(
    blst_avx_gt* out,
    size_t begin,
    size_t count,
    const blst_avx_p1* g1,
    const blst_avx_prepared_g2_batch* batch,
    blst_avx_prepared_g2_batch_workspace* workspace);

/**
 * Multiplies Miller loops for per-pair G1 points and prepared G2 operands.
 *
 * The operation computes
 * `out = product_{i=0}^{count-1} miller(g1[i], prepared_g2[indices[i]])`.
 * `out` is an unfinalized Miller-loop value, not a GT element. The caller may
 * multiply products from disjoint ranges and apply one final exponentiation to
 * the aggregate.
 *
 * Preconditions:
 * - Every G1 point is a library-produced, subgroup-validated native point.
 * - Every index is less than the prepared batch length; repetition is allowed.
 * - For nonzero `count`, `g1` and `indices` each address `count` elements.
 * - `out` does not overlap either input array or the prepared batch.
 *
 * The prepared batch is immutable, so concurrent calls are permitted. The
 * function retains no cursor or task state, allocates no memory, and is
 * variable-time with respect to the public count, indices, and point
 * identities. Empty input returns the multiplicative identity. Rejected null
 * pointers, out-of-range indices, backend unavailability, or a caught
 * exception return false without changing a nonnull `out`.
 */
bool blst_avx_prepared_g2_batch_miller_loop_product_indices(
    blst_avx_fp12* out,
    const blst_avx_p1* g1,
    const size_t* indices,
    size_t count,
    const blst_avx_prepared_g2_batch* batch);

#ifdef __cplusplus
}
#endif

#endif
