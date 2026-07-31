/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef __BLST_SCALAR_H__
#define __BLST_SCALAR_H__

#include "blst_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct blst_opaque blst_pairing;
typedef struct blst_opaque blst_uniq;

/** Result of inserting a message into a `blst_uniq` tree. */
typedef enum BLST_UNIQ_RESULT {
    BLST_UNIQ_INSERTED = 0,
    BLST_UNIQ_DUPLICATE = 1,
    BLST_UNIQ_CAPACITY_EXCEEDED = 2,
    BLST_UNIQ_DEPTH_EXCEEDED = 3,
    BLST_UNIQ_INVALID_INPUT = 4,
} BLST_UNIQ_RESULT;

/**
 * Plain scalar/public C API.
 *
 * Parameter conventions:
 * - `ret` or `out` receives the result and must point to writable storage.
 * - Input pointers must be non-null unless the function explicitly documents
 *   nullable behavior.
 * - `a`, `b`, `p`, `q`, `PK`, and `signature` name read-only inputs.
 * - Field elements are interpreted in Montgomery form.
 * - Projective and affine points use the concrete POD layouts from
 *   `blst_types.h`.
 *
 * Return conventions:
 * - `void` functions report errors only through documented invariants.
 * - `bool` results report success or the queried predicate directly.
 * - `BLST_ERROR` reports decoding, subgroup, or verification failures.
 */

/**
 * Scalar byte helpers.
 *
 * Invariants:
 * - `blst_scalar_from_{le,be}_bytes` accept variable-length inputs.
 * - Zero-length input decodes to the zero scalar, clears `out`, and returns
 *   `false`.
 * - The checked arithmetic helpers return `false` when the resulting scalar is
 *   zero or otherwise rejected by the scalar validity contract.
 */
bool blst_scalar_fr_check(const blst_scalar* a);
bool blst_sk_check(const blst_scalar* a);
bool blst_sk_add_n_check(blst_scalar* out, const blst_scalar* a, const blst_scalar* b);
bool blst_sk_sub_n_check(blst_scalar* out, const blst_scalar* a, const blst_scalar* b);
bool blst_sk_mul_n_check(blst_scalar* out, const blst_scalar* a, const blst_scalar* b);
void blst_sk_inverse(blst_scalar* out, const blst_scalar* a);
bool blst_scalar_from_le_bytes(blst_scalar* out, const byte* in, size_t len);
bool blst_scalar_from_be_bytes(blst_scalar* out, const byte* in, size_t len);

/** Adds two Montgomery-form Fr values and writes the sum to `ret`. */
void blst_fr_add(blst_fr* ret, const blst_fr* a, const blst_fr* b);
/** Subtracts `b` from `a` in Fr and writes the difference to `ret`. */
void blst_fr_sub(blst_fr* ret, const blst_fr* a, const blst_fr* b);
/** Multiplies an Fr value by the integer constant 3. */
void blst_fr_mul_by_3(blst_fr* ret, const blst_fr* a);
/** Shifts an Fr value left by `count` bits inside its limb representation. */
void blst_fr_lshift(blst_fr* ret, const blst_fr* a, size_t count);
/** Shifts an Fr value right by `count` bits inside its limb representation. */
void blst_fr_rshift(blst_fr* ret, const blst_fr* a, size_t count);
/** Multiplies two Montgomery-form Fr values and writes the product to `ret`. */
void blst_fr_mul(blst_fr* ret, const blst_fr* a, const blst_fr* b);
/**
 * Raises `a` to a public little-endian exponent and writes the result to `ret`.
 *
 * A zero-bit exponent produces one. Otherwise `nbits` identifies the highest
 * set bit and `exponent` must contain at least `(nbits + 7) / 8` readable
 * bytes. Inputs and output may alias.
 */
void blst_fr_pow(blst_fr* ret, const blst_fr* a, const byte* exponent, size_t nbits);
/** Squares a Montgomery-form Fr value and writes the result to `ret`. */
void blst_fr_sqr(blst_fr* ret, const blst_fr* a);
/** Negates `a` in Fr when `flag` is true; otherwise copies `a` into `ret`. */
void blst_fr_cneg(blst_fr* ret, const blst_fr* a, bool flag);
/** Computes the Euclidean inverse of `a` in Fr and writes it to `ret`. */
void blst_fr_eucl_inverse(blst_fr* ret, const blst_fr* a);
/** Computes the multiplicative inverse of `a` in Fr and writes it to `ret`. */
void blst_fr_inverse(blst_fr* ret, const blst_fr* a);
/** Converts a four-limb unsigned integer into Montgomery-form Fr. */
void blst_fr_from_uint64(blst_fr* ret, const uint64_t a[4]);
/** Converts a Montgomery-form Fr value into four unsigned 64-bit limbs. */
void blst_uint64_from_fr(uint64_t ret[4], const blst_fr* a);
/** Converts a scalar byte representation into Montgomery-form Fr. */
void blst_fr_from_scalar(blst_fr* ret, const blst_scalar* a);
/** Converts a Montgomery-form Fr value into the scalar byte representation. */
void blst_scalar_from_fr(blst_scalar* ret, const blst_fr* a);
/** Bulk-converts contiguous Montgomery-form field elements into canonical little-endian scalars. Empty input succeeds without storage; nonempty input requires both arrays. */
bool blst_scalars_from_fr(blst_scalar out[], const blst_fr scalars[], size_t npoints);
/** Applies the constant-time butterfly step used by the Fr FFT routines. */
void blst_fr_ct_bfly(blst_fr* x0, blst_fr* x1, const blst_fr* twiddle);
/** Applies the Gentleman-Sande butterfly step used by the Fr FFT routines. */
void blst_fr_gs_bfly(blst_fr* x0, blst_fr* x1, const blst_fr* twiddle);
/** Converts an Fr value into Montgomery form and writes it to `ret`. */
void blst_fr_to(blst_fr* ret, const blst_fr* a);
/** Converts an Fr value out of Montgomery form and writes it to `ret`. */
void blst_fr_from(blst_fr* ret, const blst_fr* a);
/** Parses a hexadecimal ASCII string into Montgomery-form Fr. */
void blst_fr_from_hexascii(blst_fr* ret, const byte* hex);

/**
 * Arithmetic over Fp. Integer conversions use the exact array widths shown in
 * the signatures. Big-endian and little-endian byte helpers always operate on
 * 48-byte field encodings.
 */
/** Adds two Montgomery-form Fp values and writes the sum to `ret`. */
void blst_fp_add(blst_fp* ret, const blst_fp* a, const blst_fp* b);
/** Subtracts `b` from `a` in Fp and writes the difference to `ret`. */
void blst_fp_sub(blst_fp* ret, const blst_fp* a, const blst_fp* b);
/** Multiplies an Fp value by the integer constant 3. */
void blst_fp_mul_by_3(blst_fp* ret, const blst_fp* a);
/** Multiplies an Fp value by the integer constant 8. */
void blst_fp_mul_by_8(blst_fp* ret, const blst_fp* a);
/** Shifts an Fp value left by `count` bits inside its limb representation. */
void blst_fp_lshift(blst_fp* ret, const blst_fp* a, size_t count);
/** Multiplies two Montgomery-form Fp values and writes the product to `ret`. */
void blst_fp_mul(blst_fp* ret, const blst_fp* a, const blst_fp* b);
/** Squares a Montgomery-form Fp value and writes the result to `ret`. */
void blst_fp_sqr(blst_fp* ret, const blst_fp* a);
/** Negates `a` in Fp when `flag` is true; otherwise copies `a` into `ret`. */
void blst_fp_cneg(blst_fp* ret, const blst_fp* a, bool flag);
/** Computes the Euclidean inverse of `a` in Fp and writes it to `ret`. */
void blst_fp_eucl_inverse(blst_fp* ret, const blst_fp* a);
/** Computes the multiplicative inverse of `a` in Fp and writes it to `ret`. */
void blst_fp_inverse(blst_fp* ret, const blst_fp* a);
/** Attempts to compute an Fp square root; returns `true` only for quadratic residues. */
bool blst_fp_sqrt(blst_fp* ret, const blst_fp* a);
/** Converts a 12-word unsigned integer into Montgomery-form Fp. */
void blst_fp_from_uint32(blst_fp* ret, const uint32_t a[12]);
/** Converts a Montgomery-form Fp value into 12 unsigned 32-bit words. */
void blst_uint32_from_fp(uint32_t ret[12], const blst_fp* a);
/** Converts a six-limb unsigned integer into Montgomery-form Fp. */
void blst_fp_from_uint64(blst_fp* ret, const uint64_t a[6]);
/** Converts a Montgomery-form Fp value into six unsigned 64-bit limbs. */
void blst_uint64_from_fp(uint64_t ret[6], const blst_fp* a);
/** Decodes a 48-byte big-endian field element into Montgomery-form Fp. */
void blst_fp_from_bendian(blst_fp* ret, const byte a[48]);
/** Encodes a Montgomery-form Fp value into a 48-byte big-endian string. */
void blst_bendian_from_fp(byte ret[48], const blst_fp* a);
/** Decodes a 48-byte little-endian field element into Montgomery-form Fp. */
void blst_fp_from_lendian(blst_fp* ret, const byte a[48]);
/**
 * Reduces an arbitrary-length little-endian unsigned integer modulo Fp.
 * `in` must reference `len` readable bytes when `len` is nonzero.
 * When `len` is zero, `in` may be null and `ret` is set to field zero.
 */
void blst_fp_from_le_bytes(blst_fp* ret, const byte* in, size_t len);
/**
 * Reduces an arbitrary-length big-endian unsigned integer modulo Fp.
 * `in` must reference `len` readable bytes when `len` is nonzero.
 * When `len` is zero, `in` may be null and `ret` is set to field zero.
 */
void blst_fp_from_be_bytes(blst_fp* ret, const byte* in, size_t len);
/** Encodes a Montgomery-form Fp value into a 48-byte little-endian string. */
void blst_lendian_from_fp(byte ret[48], const blst_fp* a);
/** Converts an Fp value into Montgomery form and writes it to `ret`. */
void blst_fp_to(blst_fp* ret, const blst_fp* a);
/** Converts an Fp value out of Montgomery form and writes it to `ret`. */
void blst_fp_from(blst_fp* ret, const blst_fp* a);
/** Returns `true` when `a` is a quadratic residue in Fp. */
bool blst_fp_is_square(const blst_fp* a);
/** Parses a hexadecimal ASCII string into Montgomery-form Fp. */
void blst_fp_from_hexascii(blst_fp* ret, const byte* hex);

/** Adds two Montgomery-form Fp2 values and writes the sum to `ret`. */
void blst_fp2_add(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b);
/** Subtracts `b` from `a` in Fp2 and writes the difference to `ret`. */
void blst_fp2_sub(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b);
/** Multiplies an Fp2 value by the integer constant 3. */
void blst_fp2_mul_by_3(blst_fp2* ret, const blst_fp2* a);
/** Multiplies an Fp2 value by the integer constant 8. */
void blst_fp2_mul_by_8(blst_fp2* ret, const blst_fp2* a);
/** Shifts an Fp2 value left by `count` bits inside its limb representation. */
void blst_fp2_lshift(blst_fp2* ret, const blst_fp2* a, size_t count);
/** Multiplies two Montgomery-form Fp2 values and writes the product to `ret`. */
void blst_fp2_mul(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b);
/** Squares a Montgomery-form Fp2 value and writes the result to `ret`. */
void blst_fp2_sqr(blst_fp2* ret, const blst_fp2* a);
/** Negates `a` in Fp2 when `flag` is true; otherwise copies `a` into `ret`. */
void blst_fp2_cneg(blst_fp2* ret, const blst_fp2* a, bool flag);
/** Computes the Euclidean inverse of `a` in Fp2 and writes it to `ret`. */
void blst_fp2_eucl_inverse(blst_fp2* ret, const blst_fp2* a);
/** Computes the multiplicative inverse of `a` in Fp2 and writes it to `ret`. */
void blst_fp2_inverse(blst_fp2* ret, const blst_fp2* a);
/** Attempts to compute an Fp2 square root; returns `true` only for quadratic residues. */
bool blst_fp2_sqrt(blst_fp2* ret, const blst_fp2* a);
/** Returns `true` when `a` is a quadratic residue in Fp2. */
bool blst_fp2_is_square(const blst_fp2* a);

/**
 * Arithmetic over Fp12 and GT-style helpers.
 *
 * Invariants:
 * - `blst_fp12_one()` returns a pointer to immutable static storage.
 * - `blst_final_exp` and `blst_fp12_finalverify` interpret inputs as Miller
 *   loop outputs or GT elements, respectively.
 */
/** Multiplies two Fp12 values and writes the product to `ret`. */
void blst_fp12_mul(blst_fp12* ret, const blst_fp12* a, const blst_fp12* b);
/** Squares an Fp12 value and writes the result to `ret`. */
void blst_fp12_sqr(blst_fp12* ret, const blst_fp12* a);
/** Applies the cyclotomic squaring shortcut to `a` and writes the result to `ret`. */
void blst_fp12_cyclotomic_sqr(blst_fp12* ret, const blst_fp12* a);
/** Multiplies `a` by a sparse `xy00z0` Fp6 value and writes the result to `ret`. */
void blst_fp12_mul_by_xy00z0(blst_fp12* ret, const blst_fp12* a, const blst_fp6* xy00z0);
/** Replaces `a` with its Fp12 conjugate in place. */
void blst_fp12_conjugate(blst_fp12* a);
/** Computes the multiplicative inverse of `a` in Fp12 and writes it to `ret`. */
void blst_fp12_inverse(blst_fp12* ret, const blst_fp12* a);
/** Applies the `n`th Frobenius map to `a` and writes the result to `ret`. */
void blst_fp12_frobenius_map(blst_fp12* ret, const blst_fp12* a, size_t n);
/** Returns `true` when the two Fp12 values are equal. */
bool blst_fp12_is_equal(const blst_fp12* a, const blst_fp12* b);
/** Returns `true` when `a` is the multiplicative identity. */
bool blst_fp12_is_one(const blst_fp12* a);
/** Returns `true` when `a` lies in the target subgroup used by pairings. */
bool blst_fp12_in_group(const blst_fp12* a);
/** Returns a pointer to immutable static storage holding the Fp12 identity. */
const blst_fp12* blst_fp12_one(void);
/** Applies the pairing final exponentiation to `f` and writes the result to `ret`. */
void blst_final_exp(blst_fp12* ret, const blst_fp12* f);
/** Returns immutable static storage holding the additive GT identity. */
const blst_fp12* blst_gt_identity(void);
/** Returns immutable static storage holding pairing(G1 generator, G2 generator). */
const blst_fp12* blst_gt_generator(void);
/** Returns `true` when target-group elements `a` and `b` are equal. */
bool blst_gt_is_equal(const blst_fp12* a, const blst_fp12* b);
/** Returns `true` when target-group element `p` is the additive identity. */
bool blst_gt_is_identity(const blst_fp12* p);
/** Adds target-group elements `a` and `b`. Inputs and output may alias. */
void blst_gt_add(blst_fp12* out, const blst_fp12* a, const blst_fp12* b);
/** Subtracts target-group element `b` from `a`. Inputs and output may alias. */
void blst_gt_sub(blst_fp12* out, const blst_fp12* a, const blst_fp12* b);
/** Negates target-group element `p`. Input and output may alias. */
void blst_gt_neg(blst_fp12* out, const blst_fp12* p);
/** Doubles target-group element `p`. Input and output may alias. */
void blst_gt_double(blst_fp12* out, const blst_fp12* p);
/**
 * Multiplies target-group element `p` by a little-endian scalar.
 *
 * `p` must already lie in GT and `nbits` must not exceed 256. The operation
 * is cyclotomic exponentiation; it does not validate subgroup membership.
 */
void blst_gt_mult(blst_fp12* out, const blst_fp12* p, const byte* scalar, size_t nbits);
/** Returns `true` when `gt1` and `gt2` pass the final GT verification check. */
bool blst_fp12_finalverify(const blst_fp12* gt1, const blst_fp12* gt2);
/** Encodes an Fp12 value into a 576-byte big-endian string. */
void blst_bendian_from_fp12(byte out[48 * 12], const blst_fp12* a);

/**
 * G1 projective and affine helpers.
 *
 * Parameter interpretation:
 * - Scalar multiplication reads `scalar` as an unsigned little-endian bitstring
 *   with exactly `nbits` meaningful bits.
 * - Serialization writes 96-byte uncompressed or 48-byte compressed encodings.
 * - Generator accessors return immutable static storage.
 */
/** Adds two G1 projective points and writes the sum to `out`. */
void blst_p1_add(blst_p1* out, const blst_p1* a, const blst_p1* b);
/** Subtracts G1 projective point `b` from `a` and writes the result to `out`. */
void blst_p1_sub(blst_p1* out, const blst_p1* a, const blst_p1* b);
/** Adds `a` and `b`, doubling when they represent the same G1 point. */
void blst_p1_add_or_double(blst_p1* out, const blst_p1* a, const blst_p1* b);
/** Adds a projective G1 point and an affine G1 point and writes the sum to `out`. */
void blst_p1_add_affine(blst_p1* out, const blst_p1* a, const blst_p1_affine* b);
/** Subtracts affine G1 point `b` from projective `a` and writes the result to `out`. */
void blst_p1_sub_affine(blst_p1* out, const blst_p1* a, const blst_p1_affine* b);
/** Adds or doubles using a projective `a` and affine `b` G1 input. */
void blst_p1_add_or_double_affine(blst_p1* out, const blst_p1* a, const blst_p1_affine* b);
/** Doubles the G1 projective point `a` and writes the result to `out`. */
void blst_p1_double(blst_p1* out, const blst_p1* a);
/** Multiplies G1 point `p` by a little-endian scalar with `nbits` meaningful bits. */
void blst_p1_mult(blst_p1* out, const blst_p1* p, const byte* scalar, size_t nbits);
/** Negates `p` in place when `cbit` is true; otherwise leaves it unchanged. */
void blst_p1_cneg(blst_p1* p, bool cbit);
/** Conditionally negates affine G1 point `p` in place. */
void blst_p1_affine_cneg(blst_p1_affine* p, bool cbit);
/** Converts a G1 projective point into affine form. */
void blst_p1_to_affine(blst_p1_affine* out, const blst_p1* in);
/** Converts a G1 affine point into projective form. */
void blst_p1_from_affine(blst_p1* out, const blst_p1_affine* in);
/** Normalizes a G1 Jacobian point into the canonical projective form. */
void blst_p1_from_jacobian(blst_p1* out, const blst_p1* in);
/** Returns `true` when `p` lies on the G1 curve equation. */
bool blst_p1_on_curve(const blst_p1* p);
/** Returns `true` when `p` lies in the prime-order G1 subgroup. */
bool blst_p1_in_g1(const blst_p1* p);
/** Returns `true` when the two projective G1 points represent the same point. */
bool blst_p1_is_equal(const blst_p1* a, const blst_p1* b);
/** Returns `true` when `a` is the G1 point at infinity. */
bool blst_p1_is_inf(const blst_p1* a);
/** Returns a pointer to immutable static storage holding the G1 generator. */
const blst_p1* blst_p1_generator(void);
/** Multiplies `p` by `scalar` without subgroup safety checks. */
void blst_p1_unchecked_mult(blst_p1* out, const blst_p1* p, const byte* scalar, size_t nbits);
/** Returns `true` when affine G1 point `p` lies on the curve. */
bool blst_p1_affine_on_curve(const blst_p1_affine* p);
/** Returns `true` when affine G1 point `p` lies in the prime-order subgroup. */
bool blst_p1_affine_in_g1(const blst_p1_affine* p);
/** Returns `true` when the two affine G1 points represent the same point. */
bool blst_p1_affine_is_equal(const blst_p1_affine* a, const blst_p1_affine* b);
/** Returns `true` when affine G1 point `a` is the point at infinity. */
bool blst_p1_affine_is_inf(const blst_p1_affine* a);
/** Returns a pointer to immutable static storage holding the affine G1 generator. */
const blst_p1_affine* blst_p1_affine_generator(void);
/** Batch-converts `npoints` projective G1 inputs into affine outputs; when `npoints > 0`, `points[0]` must be non-null and later null entries mean contiguous continuation from the previous non-null element. */
void blst_p1s_to_affine(blst_p1_affine dst[], const blst_p1* const points[], size_t npoints);
/** Adds `npoints` affine G1 inputs and writes the total into `ret`; when `npoints > 0`, `points[0]` must be non-null and later null entries mean contiguous continuation from the previous non-null element. */
void blst_p1s_add(blst_p1* ret, const blst_p1_affine* const points[], size_t npoints);
/** Returns the table size needed for G1 wbits precomputation, or 0 when `wbits` is outside the supported range `[2, 14]` or the size computation would overflow. */
size_t blst_p1s_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints);
/** Precomputes a G1 wbits multiplication table for `npoints` affine inputs; `wbits` must lie in `[2, 14]` and `points[0]` must be non-null when `npoints > 0`. */
void blst_p1s_mult_wbits_precompute(blst_p1_affine table[], size_t wbits, const blst_p1_affine* const points[], size_t npoints);
/** Returns the scratch size needed by the G1 wbits MSM routine. */
size_t blst_p1s_mult_wbits_scratch_sizeof(size_t npoints);
/** Evaluates a G1 MSM from a wbits precompute table and caller-owned scratch; `wbits` must lie in `[2, 14]` and `scalars[0]` must be non-null when `npoints > 0`. */
void blst_p1s_mult_wbits(blst_p1* ret, const blst_p1_affine table[], size_t wbits, size_t npoints, const byte* const scalars[], size_t nbits, limb_t* scratch);
/** Returns the scratch size needed by the G1 Pippenger MSM routine. */
size_t blst_p1s_mult_pippenger_scratch_sizeof(size_t npoints);
/** Returns the scratch size needed by the G1 Pippenger MSM routine for an explicit window, or 0 when `window` is invalid. */
size_t blst_p1s_mult_pippenger_window_scratch_sizeof(size_t window);
/** Evaluates a G1 MSM with Pippenger using caller-owned scratch memory; when `npoints == 1`, the wrapper falls back to single-point multiplication, and when `npoints > 1`, `scratch`, `points[0]`, and `scalars[0]` must be non-null. Later null entries mean contiguous continuation from the previous non-null element. */
void blst_p1s_mult_pippenger(blst_p1* ret, const blst_p1_affine* const points[], size_t npoints, const byte* const scalars[], size_t nbits, limb_t* scratch);
/** Returns total scratch bytes for a contiguous G1 MSM from native Montgomery-form `blst_fr` scalars, including canonical scalar preparation and Pippenger buckets. Zero is valid for an empty MSM and otherwise reports an unrepresentable size. */
size_t blst_p1s_mult_pippenger_fr_scratch_sizeof(size_t npoints);
/** Evaluates a contiguous G1 MSM from native Montgomery-form `blst_fr` scalars. Scratch must be limb-aligned and contain at least the queried byte count. Empty input writes the identity and requires no input or scratch storage. */
bool blst_p1s_mult_pippenger_fr(blst_p1* ret, const blst_p1_affine points[], const blst_fr scalars[], size_t npoints, limb_t* scratch, size_t scratch_size);
/** Returns caller-owned scratch bytes for the selected contiguous G1 MSM implementation. This is a capacity query: storage for `max_points` supports every input length through that capacity. */
size_t blst_p1s_mult_fr_scratch_sizeof(size_t max_points);
/** Evaluates a contiguous G1 MSM from native Montgomery-form `blst_fr` scalars using the measured scalar-backend policy. Inputs must be prevalidated subgroup points; identity points are accepted and use the Pippenger fallback. */
bool blst_p1s_mult_fr(blst_p1* ret, const blst_p1_affine points[], const blst_fr scalars[], size_t npoints, limb_t* scratch, size_t scratch_size);
/** Returns total scratch bytes for the explicit G1 batch-affine MSM, including canonical scalar preparation. Supported batch sizes are 128 and 384; `window` must be non-zero and no larger than 32. Zero denotes either empty input or invalid/unrepresentable parameters. */
size_t blst_p1s_mult_batch_affine_fr_scratch_sizeof(size_t npoints, size_t window, size_t batch_size);
/** Evaluates a contiguous G1 MSM with the explicit batch-affine algorithm. Inputs are public, prevalidated, non-identity subgroup points and native Montgomery-form scalars. Scratch must be limb-aligned and contain at least the queried byte count. This routine performs no allocation and no input validation. */
bool blst_p1s_mult_batch_affine_fr(blst_p1* ret, const blst_p1_affine points[], const blst_fr scalars[], size_t npoints, limb_t* scratch, size_t scratch_size, size_t window, size_t batch_size);
/** Evaluates a G1 MSM with Pippenger using an explicit window size. */
void blst_p1s_mult_pippenger_window(blst_p1* ret, const blst_p1_affine* const points[], size_t npoints, const byte* const scalars[], size_t nbits, limb_t* scratch, size_t window);
/** Evaluates one tiled window of a G1 Pippenger MSM; when `npoints > 0`, `scratch`, `points[0]`, and `scalars[0]` must be non-null, `window` must be non-zero, and `bit0` must be less than `nbits`. Later null entries mean contiguous continuation from the previous non-null element. */
void blst_p1s_tile_pippenger(blst_p1* ret, const blst_p1_affine* const points[], size_t npoints, const byte* const scalars[], size_t nbits, limb_t* scratch, size_t bit0, size_t window);
/** Serializes a G1 projective point into a 96-byte uncompressed encoding. */
void blst_p1_serialize(byte out[96], const blst_p1* in);
/** Compresses a G1 projective point into a 48-byte encoding. */
void blst_p1_compress(byte out[48], const blst_p1* in);
/** Serializes an affine G1 point into a 96-byte uncompressed encoding. */
void blst_p1_affine_serialize(byte out[96], const blst_p1_affine* in);
/** Compresses an affine G1 point into a 48-byte encoding. */
void blst_p1_affine_compress(byte out[48], const blst_p1_affine* in);
/** Decodes a 48-byte compressed G1 encoding into `out`. */
BLST_ERROR blst_p1_uncompress(blst_p1_affine* out, const byte in[48]);
/** Decodes a 96-byte uncompressed G1 encoding into `out`. */
BLST_ERROR blst_p1_deserialize(blst_p1_affine* out, const byte in[96]);

/**
 * G2 projective and affine helpers.
 *
 * Parameter interpretation matches the G1 section, with 192-byte uncompressed
 * and 96-byte compressed encodings.
 */
/** Adds two G2 projective points and writes the sum to `out`. */
void blst_p2_add(blst_p2* out, const blst_p2* a, const blst_p2* b);
/** Subtracts G2 projective point `b` from `a` and writes the result to `out`. */
void blst_p2_sub(blst_p2* out, const blst_p2* a, const blst_p2* b);
/** Adds `a` and `b`, doubling when they represent the same G2 point. */
void blst_p2_add_or_double(blst_p2* out, const blst_p2* a, const blst_p2* b);
/** Adds a projective G2 point and an affine G2 point and writes the sum to `out`. */
void blst_p2_add_affine(blst_p2* out, const blst_p2* a, const blst_p2_affine* b);
/** Subtracts affine G2 point `b` from projective `a` and writes the result to `out`. */
void blst_p2_sub_affine(blst_p2* out, const blst_p2* a, const blst_p2_affine* b);
/** Adds or doubles using a projective `a` and affine `b` G2 input. */
void blst_p2_add_or_double_affine(blst_p2* out, const blst_p2* a, const blst_p2_affine* b);
/** Doubles the G2 projective point `a` and writes the result to `out`. */
void blst_p2_double(blst_p2* out, const blst_p2* a);
/** Multiplies G2 point `p` by a little-endian scalar with `nbits` meaningful bits. */
void blst_p2_mult(blst_p2* out, const blst_p2* p, const byte* scalar, size_t nbits);
/** Negates `p` in place when `cbit` is true; otherwise leaves it unchanged. */
void blst_p2_cneg(blst_p2* p, bool cbit);
/** Conditionally negates affine G2 point `p` in place. */
void blst_p2_affine_cneg(blst_p2_affine* p, bool cbit);
/** Converts a G2 projective point into affine form. */
void blst_p2_to_affine(blst_p2_affine* out, const blst_p2* in);
/** Converts a G2 affine point into projective form. */
void blst_p2_from_affine(blst_p2* out, const blst_p2_affine* in);
/** Normalizes a G2 Jacobian point into the canonical projective form. */
void blst_p2_from_jacobian(blst_p2* out, const blst_p2* in);
/** Returns `true` when `p` lies on the G2 curve equation. */
bool blst_p2_on_curve(const blst_p2* p);
/** Returns `true` when `p` lies in the prime-order G2 subgroup. */
bool blst_p2_in_g2(const blst_p2* p);
/** Returns `true` when the two projective G2 points represent the same point. */
bool blst_p2_is_equal(const blst_p2* a, const blst_p2* b);
/** Returns `true` when `a` is the G2 point at infinity. */
bool blst_p2_is_inf(const blst_p2* a);
/** Returns a pointer to immutable static storage holding the G2 generator. */
const blst_p2* blst_p2_generator(void);
/** Multiplies `p` by `scalar` without subgroup safety checks. */
void blst_p2_unchecked_mult(blst_p2* out, const blst_p2* p, const byte* scalar, size_t nbits);
/** Returns `true` when affine G2 point `p` lies on the curve. */
bool blst_p2_affine_on_curve(const blst_p2_affine* p);
/** Returns `true` when affine G2 point `p` lies in the prime-order subgroup. */
bool blst_p2_affine_in_g2(const blst_p2_affine* p);
/** Returns `true` when the two affine G2 points represent the same point. */
bool blst_p2_affine_is_equal(const blst_p2_affine* a, const blst_p2_affine* b);
/** Returns `true` when affine G2 point `a` is the point at infinity. */
bool blst_p2_affine_is_inf(const blst_p2_affine* a);
/** Returns a pointer to immutable static storage holding the affine G2 generator. */
const blst_p2_affine* blst_p2_affine_generator(void);
/** Batch-converts `npoints` projective G2 inputs into affine outputs; when `npoints > 0`, `points[0]` must be non-null and later null entries mean contiguous continuation from the previous non-null element. */
void blst_p2s_to_affine(blst_p2_affine dst[], const blst_p2* const points[], size_t npoints);
/** Adds `npoints` affine G2 inputs and writes the total into `ret`; when `npoints > 0`, `points[0]` must be non-null and later null entries mean contiguous continuation from the previous non-null element. */
void blst_p2s_add(blst_p2* ret, const blst_p2_affine* const points[], size_t npoints);
/** Returns the table size needed for G2 wbits precomputation, or 0 when `wbits` is outside the supported range `[2, 14]` or the size computation would overflow. */
size_t blst_p2s_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints);
/** Precomputes a G2 wbits multiplication table for `npoints` affine inputs; `wbits` must lie in `[2, 14]` and `points[0]` must be non-null when `npoints > 0`. */
void blst_p2s_mult_wbits_precompute(blst_p2_affine table[], size_t wbits, const blst_p2_affine* const points[], size_t npoints);
/** Returns the scratch size needed by the G2 wbits MSM routine. */
size_t blst_p2s_mult_wbits_scratch_sizeof(size_t npoints);
/** Evaluates a G2 MSM from a wbits precompute table and caller-owned scratch; `wbits` must lie in `[2, 14]` and `scalars[0]` must be non-null when `npoints > 0`. */
void blst_p2s_mult_wbits(blst_p2* ret, const blst_p2_affine table[], size_t wbits, size_t npoints, const byte* const scalars[], size_t nbits, limb_t* scratch);
/** Returns the scratch size needed by the G2 Pippenger MSM routine. */
size_t blst_p2s_mult_pippenger_scratch_sizeof(size_t npoints);
/** Returns the scratch size needed by the G2 Pippenger MSM routine for an explicit window, or 0 when `window` is invalid. */
size_t blst_p2s_mult_pippenger_window_scratch_sizeof(size_t window);
/** Evaluates a G2 MSM with Pippenger using caller-owned scratch memory; when `npoints == 1`, the wrapper falls back to single-point multiplication, and when `npoints > 1`, `scratch`, `points[0]`, and `scalars[0]` must be non-null. Later null entries mean contiguous continuation from the previous non-null element. */
void blst_p2s_mult_pippenger(blst_p2* ret, const blst_p2_affine* const points[], size_t npoints, const byte* const scalars[], size_t nbits, limb_t* scratch);
/** Evaluates a G2 MSM with Pippenger using an explicit window size. */
void blst_p2s_mult_pippenger_window(blst_p2* ret, const blst_p2_affine* const points[], size_t npoints, const byte* const scalars[], size_t nbits, limb_t* scratch, size_t window);
/** Returns total scratch bytes for a contiguous G2 MSM from native Montgomery-form `blst_fr` scalars, including canonical scalar preparation and Pippenger buckets. Zero is valid for an empty MSM and otherwise reports an unrepresentable size. */
size_t blst_p2s_mult_pippenger_fr_scratch_sizeof(size_t npoints);
/** Evaluates a contiguous G2 MSM from native Montgomery-form `blst_fr` scalars. Scratch must be limb-aligned and contain at least the queried byte count. Empty input writes the identity and requires no input or scratch storage. */
bool blst_p2s_mult_pippenger_fr(blst_p2* ret, const blst_p2_affine points[], const blst_fr scalars[], size_t npoints, limb_t* scratch, size_t scratch_size);
/** Returns caller-owned scratch bytes for the selected contiguous G2 MSM implementation. This is a capacity query: storage for `max_points` supports every input length through that capacity. */
size_t blst_p2s_mult_fr_scratch_sizeof(size_t max_points);
/** Evaluates a contiguous G2 MSM from native Montgomery-form `blst_fr` scalars using the measured scalar-backend policy. Inputs must be prevalidated subgroup points; identity points are accepted and use the Pippenger fallback. */
bool blst_p2s_mult_fr(blst_p2* ret, const blst_p2_affine points[], const blst_fr scalars[], size_t npoints, limb_t* scratch, size_t scratch_size);
/** Returns total scratch bytes for the explicit G2 batch-affine MSM, including canonical scalar preparation. Supported batch sizes are 128 and 384; `window` must be non-zero and no larger than 32. Zero denotes either empty input or invalid/unrepresentable parameters. */
size_t blst_p2s_mult_batch_affine_fr_scratch_sizeof(size_t npoints, size_t window, size_t batch_size);
/** Evaluates a contiguous G2 MSM with the explicit batch-affine algorithm. Inputs are public, prevalidated, non-identity subgroup points and native Montgomery-form scalars. Scratch must be limb-aligned and contain at least the queried byte count. This routine performs no allocation and no input validation. */
bool blst_p2s_mult_batch_affine_fr(blst_p2* ret, const blst_p2_affine points[], const blst_fr scalars[], size_t npoints, limb_t* scratch, size_t scratch_size, size_t window, size_t batch_size);
/** Evaluates one tiled window of a G2 Pippenger MSM; when `npoints > 0`, `scratch`, `points[0]`, and `scalars[0]` must be non-null, `window` must be non-zero, and `bit0` must be less than `nbits`. Later null entries mean contiguous continuation from the previous non-null element. */
void blst_p2s_tile_pippenger(blst_p2* ret, const blst_p2_affine* const points[], size_t npoints, const byte* const scalars[], size_t nbits, limb_t* scratch, size_t bit0, size_t window);
/** Serializes a G2 projective point into a 192-byte uncompressed encoding. */
void blst_p2_serialize(byte out[192], const blst_p2* in);
/** Compresses a G2 projective point into a 96-byte encoding. */
void blst_p2_compress(byte out[96], const blst_p2* in);
/** Serializes an affine G2 point into a 192-byte uncompressed encoding. */
void blst_p2_affine_serialize(byte out[192], const blst_p2_affine* in);
/** Compresses an affine G2 point into a 96-byte encoding. */
void blst_p2_affine_compress(byte out[96], const blst_p2_affine* in);
/** Decodes a 96-byte compressed G2 encoding into `out`. */
BLST_ERROR blst_p2_uncompress(blst_p2_affine* out, const byte in[96]);
/** Decodes a 192-byte uncompressed G2 encoding into `out`. */
BLST_ERROR blst_p2_deserialize(blst_p2_affine* out, const byte in[192]);

/**
 * Pairing and aggregate verification helpers.
 *
 * Invariants:
 * - `blst_pairing_sizeof()` and `blst_pairing_alignof()` describe the exact
 *   storage requirements for `blst_pairing`.
 * - `blst_pairing_init()` normally borrows `DST`; the caller must keep those
 *   bytes alive and unchanged while the context can still aggregate messages.
 * - When `DST` points exactly `blst_pairing_sizeof()` bytes past `ctx`, the
 *   implementation stores a relocatable self-relative reference. In that
 *   layout, the caller owns one allocation containing `[ctx | DST]`.
 * - `blst_pairing_get_dst()` returns the active borrowed DST pointer; its
 *   lifetime is governed by the preceding ownership rule.
 * - `msg` and `aug` are byte strings of the supplied lengths.
 */
/** Computes a Miller loop from affine `(Q, P)` inputs and writes it to `ret`. */
void blst_miller_loop(blst_fp12* ret, const blst_p2_affine* Q, const blst_p1_affine* P);
/**
 * Multiplies `n` Miller loops together from paired affine input arrays.
 *
 * Pointer-array contract:
 * - when `n == 0`, `ret` receives the Fp12 multiplicative identity
 * - when `n > 0`, `Qs[0]` and `Ps[0]` must be non-null
 * - a later null entry means "continue from the next contiguous element after
 *   the previous non-null entry"
 */
void blst_miller_loop_n(blst_fp12* ret, const blst_p2_affine* const Qs[], const blst_p1_affine* const Ps[], size_t n);
/** Precomputes the 68 line-function coefficients needed for repeated Miller loops with `Q`; when `Q` is null, the G2 generator is used. */
void blst_precompute_lines(blst_fp6 Qlines[68], const blst_p2_affine* Q);
/** Runs a Miller loop using a precomputed `Qlines` table and affine `P`; when `P` is null, the G1 generator is used. */
void blst_miller_loop_lines(blst_fp12* ret, const blst_fp6 Qlines[68], const blst_p1_affine* P);
/** Returns the exact storage size, in bytes, required for one `blst_pairing` object. */
size_t blst_pairing_sizeof(void);
/** Returns the required alignment for `blst_pairing` storage. */
size_t blst_pairing_alignof(void);
/** Initializes caller-owned `blst_pairing` storage with the given DST and hash/encode mode. */
void blst_pairing_init(blst_pairing* new_ctx, bool hash_or_encode, const byte* DST, size_t DST_len);
/** Returns the borrowed DST pointer currently stored in `ctx`. */
const byte* blst_pairing_get_dst(const blst_pairing* ctx);
/** Finalizes pending work in `ctx` before merge or verification. */
void blst_pairing_commit(blst_pairing* ctx);
/** Returns the current pairing accumulator as an Fp12 value by value. */
blst_fp12 blst_pairing_as_fp12(blst_pairing* ctx);
/** Writes the current pairing accumulator into caller-owned `out`. */
void blst_pairing_as_fp12_to(blst_fp12* out, blst_pairing* ctx);
/** Aggregates one raw affine `(q, p)` pair into the pairing context. */
void blst_pairing_raw_aggregate(blst_pairing* ctx, const blst_p2_affine* q, const blst_p1_affine* p);
/** Aggregates one PK-in-G2 signature verification item into `ctx`. */
BLST_ERROR blst_pairing_aggregate_pk_in_g2(blst_pairing* ctx, const blst_p2_affine* PK, const blst_p1_affine* signature, const byte* msg, size_t msg_len, const byte* aug, size_t aug_len);
/** Aggregates one PK-in-G2 verification item with optional subgroup checks. */
BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g2(blst_pairing* ctx, const blst_p2_affine* PK, bool pk_grpchk, const blst_p1_affine* signature, bool sig_grpchk, const byte* msg, size_t msg_len, const byte* aug, size_t aug_len);
/** Aggregates one PK-in-G2 verification item weighted by `scalar`. */
BLST_ERROR blst_pairing_mul_n_aggregate_pk_in_g2(blst_pairing* ctx, const blst_p2_affine* PK, const blst_p1_affine* sig, const byte* scalar, size_t nbits, const byte* msg, size_t msg_len, const byte* aug, size_t aug_len);
/** Aggregates one weighted PK-in-G2 verification item with optional subgroup checks. */
BLST_ERROR blst_pairing_chk_n_mul_n_aggr_pk_in_g2(blst_pairing* ctx, const blst_p2_affine* PK, bool pk_grpchk, const blst_p1_affine* sig, bool sig_grpchk, const byte* scalar, size_t nbits, const byte* msg, size_t msg_len, const byte* aug, size_t aug_len);
/** Aggregates one PK-in-G1 signature verification item into `ctx`. */
BLST_ERROR blst_pairing_aggregate_pk_in_g1(blst_pairing* ctx, const blst_p1_affine* PK, const blst_p2_affine* signature, const byte* msg, size_t msg_len, const byte* aug, size_t aug_len);
/** Aggregates one PK-in-G1 verification item with optional subgroup checks. */
BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g1(blst_pairing* ctx, const blst_p1_affine* PK, bool pk_grpchk, const blst_p2_affine* signature, bool sig_grpchk, const byte* msg, size_t msg_len, const byte* aug, size_t aug_len);
/** Aggregates one PK-in-G1 verification item weighted by `scalar`. */
BLST_ERROR blst_pairing_mul_n_aggregate_pk_in_g1(blst_pairing* ctx, const blst_p1_affine* PK, const blst_p2_affine* sig, const byte* scalar, size_t nbits, const byte* msg, size_t msg_len, const byte* aug, size_t aug_len);
/** Aggregates one weighted PK-in-G1 verification item with optional subgroup checks. */
BLST_ERROR blst_pairing_chk_n_mul_n_aggr_pk_in_g1(blst_pairing* ctx, const blst_p1_affine* PK, bool pk_grpchk, const blst_p2_affine* sig, bool sig_grpchk, const byte* scalar, size_t nbits, const byte* msg, size_t msg_len, const byte* aug, size_t aug_len);
/** Merges `ctx1` into `ctx`; both contexts must be compatible. */
BLST_ERROR blst_pairing_merge(blst_pairing* ctx, const blst_pairing* ctx1);
/** Returns `true` when `ctx` verifies against `gtsig`, or against the identity when `gtsig` is null. */
bool blst_pairing_finalverify(const blst_pairing* ctx, const blst_fp12* gtsig);
/** Aggregates one encoded or explicit G1 signature contribution into `out`. */
BLST_ERROR blst_aggregate_in_g1(blst_p1* out, const blst_p1* in, const byte* zwire);
/** Aggregates one encoded or explicit G2 signature contribution into `out`. */
BLST_ERROR blst_aggregate_in_g2(blst_p2* out, const blst_p2* in, const byte* zwire);
/** Computes the pairing-side aggregate value for a G1 affine signature. */
void blst_aggregated_in_g1(blst_fp12* out, const blst_p1_affine* signature);
/** Computes the pairing-side aggregate value for a G2 affine signature. */
void blst_aggregated_in_g2(blst_fp12* out, const blst_p2_affine* signature);
/** Verifies one PK-in-G1 signature against `msg`, `DST`, and optional `aug`. */
BLST_ERROR blst_core_verify_pk_in_g1(const blst_p1_affine* pk, const blst_p2_affine* signature, bool hash_or_encode, const byte* msg, size_t msg_len, const byte* DST, size_t DST_len, const byte* aug, size_t aug_len);
/** Verifies one PK-in-G2 signature against `msg`, `DST`, and optional `aug`. */
BLST_ERROR blst_core_verify_pk_in_g2(const blst_p2_affine* pk, const blst_p1_affine* signature, bool hash_or_encode, const byte* msg, size_t msg_len, const byte* DST, size_t DST_len, const byte* aug, size_t aug_len);

/**
 * Key derivation, hashing, mapping, and utility helpers.
 *
 * Invariants:
 * - `IKM`, `info`, `salt`, `msg`, `DST`, and `aug` are raw byte strings of the
 *   supplied lengths.
 * - `blst_uniq_sizeof()` returns the exact caller allocation size for the
 *   opaque `blst_uniq` object, or zero if the capacity is invalid.
 * - `blst_expand_message_xmd` implements the XMD expand step used by the hash
 *   to curve suite.
 */
/** Derives a secret key from `IKM` and optional `info` using the default keygen scheme; `IKM_len` must be at least 32 bytes or `out_SK` is cleared to zero. */
void blst_keygen(blst_scalar* out_SK, const byte* IKM, size_t IKM_len, const byte* info, size_t info_len);
/** Derives a secret key using the v3 keygen scheme; `IKM_len` must be at least 32 bytes or `out_SK` is cleared to zero. */
void blst_keygen_v3(blst_scalar* out_SK, const byte* IKM, size_t IKM_len, const byte* info, size_t info_len);
/** Derives a secret key using the v4.5 keygen scheme with explicit salt; `IKM_len` must be at least 32 bytes or `out_SK` is cleared to zero. */
void blst_keygen_v4_5(blst_scalar* out_SK, const byte* IKM, size_t IKM_len, const byte* salt, size_t salt_len, const byte* info, size_t info_len);
/** Derives a secret key using the v5 keygen scheme with explicit salt; `IKM_len` must be at least 32 bytes, `salt` must be non-null, and invalid input clears `out_SK` to zero. */
void blst_keygen_v5(blst_scalar* out_SK, const byte* IKM, size_t IKM_len, const byte* salt, size_t salt_len, const byte* info, size_t info_len);
/** Derives the EIP-2333 master secret key from `IKM`; `IKM_len` must be at least 32 bytes or `out_SK` is cleared to zero. */
void blst_derive_master_eip2333(blst_scalar* out_SK, const byte* IKM, size_t IKM_len);
/** Derives the EIP-2333 child secret key for `child_index` from `SK`. */
void blst_derive_child_eip2333(blst_scalar* out_SK, const blst_scalar* SK, uint32_t child_index);
/** Encodes `msg` to a deterministic G1 point using `DST` and optional `aug`. */
void blst_encode_to_g1(blst_p1* out, const byte* msg, size_t msg_len, const byte* DST, size_t DST_len, const byte* aug, size_t aug_len);
/** Hashes `msg` to a random-oracle G1 point using `DST` and optional `aug`. */
void blst_hash_to_g1(blst_p1* out, const byte* msg, size_t msg_len, const byte* DST, size_t DST_len, const byte* aug, size_t aug_len);
/** Encodes `msg` to a deterministic G2 point using `DST` and optional `aug`. */
void blst_encode_to_g2(blst_p2* out, const byte* msg, size_t msg_len, const byte* DST, size_t DST_len, const byte* aug, size_t aug_len);
/** Hashes `msg` to a random-oracle G2 point using `DST` and optional `aug`. */
void blst_hash_to_g2(blst_p2* out, const byte* msg, size_t msg_len, const byte* DST, size_t DST_len, const byte* aug, size_t aug_len);
/** Maps field elements `(u, v)` into a G1 point without hashing. */
void blst_map_to_g1(blst_p1* out, const blst_fp* u, const blst_fp* v);
/** Maps field elements `(u, v)` into a G2 point without hashing. */
void blst_map_to_g2(blst_p2* out, const blst_fp2* u, const blst_fp2* v);
/** Derives the projective G1 public key corresponding to `SK`. */
void blst_sk_to_pk_in_g1(blst_p1* out_pk, const blst_scalar* SK);
/** Signs a G2 hash point with `SK` and writes the projective signature to `out_sig`. */
void blst_sign_pk_in_g1(blst_p2* out_sig, const blst_p2* hash, const blst_scalar* SK);
/** Derives the projective G2 public key corresponding to `SK`. */
void blst_sk_to_pk_in_g2(blst_p2* out_pk, const blst_scalar* SK);
/** Signs a G1 hash point with `SK` and writes the projective signature to `out_sig`. */
void blst_sign_pk_in_g2(blst_p1* out_sig, const blst_p1* hash, const blst_scalar* SK);
/** Derives the affine G1 public key from `SK` and optionally serializes it into `out`. */
void blst_sk_to_pk2_in_g1(byte out[96], blst_p1_affine* out_pk, const blst_scalar* SK);
/** Signs a G2 hash point, writes the affine signature, and optionally serializes it into `out`. */
void blst_sign_pk2_in_g1(byte out[192], blst_p2_affine* out_sig, const blst_p2* hash, const blst_scalar* SK);
/** Derives the affine G2 public key from `SK` and optionally serializes it into `out`. */
void blst_sk_to_pk2_in_g2(byte out[192], blst_p2_affine* out_pk, const blst_scalar* SK);
/** Signs a G1 hash point, writes the affine signature, and optionally serializes it into `out`. */
void blst_sign_pk2_in_g2(byte out[96], blst_p1_affine* out_sig, const blst_p1* hash, const blst_scalar* SK);
/** Returns the exact allocation size required for a `blst_uniq` tree with `capacity` nodes; returns 0 for zero capacity or size overflow. */
size_t blst_uniq_sizeof(size_t capacity);
/** Returns the required alignment for caller-owned `blst_uniq` storage. */
size_t blst_uniq_alignof(void);
/** Initializes suitably aligned caller-owned storage for exactly `capacity` nodes; returns false for null or misaligned storage, or an invalid capacity. */
bool blst_uniq_init(blst_uniq* tree, size_t capacity);
/**
 * Inserts `msg[0..len)` and reports the exact outcome.
 *
 * `tree` must have been successfully initialized, and concurrent insertions
 * into the same tree require external synchronization.
 * `msg` may be null only when `len == 0`. The message bytes are borrowed, not
 * copied, and must remain readable and unchanged until the tree is no longer
 * used. A duplicate is still reported as such when the tree is at capacity.
 */
BLST_UNIQ_RESULT blst_uniq_insert(blst_uniq* tree, const byte* msg, size_t len);
/** Expands `msg` into `out_len` bytes using the XMD construction and `DST`; unsupported lengths or internal allocation failure clear `out` to zero. */
void blst_expand_message_xmd(byte* out, size_t out_len, const byte* msg, size_t msg_len, const byte* DST, size_t DST_len);
/** Returns `sizeof(blst_p1)` for ABI consumers that need runtime verification. */
size_t blst_p1_sizeof(void);
/** Returns `sizeof(blst_p1_affine)` for ABI consumers that need runtime verification. */
size_t blst_p1_affine_sizeof(void);
/** Returns `sizeof(blst_p2)` for ABI consumers that need runtime verification. */
size_t blst_p2_sizeof(void);
/** Returns `sizeof(blst_p2_affine)` for ABI consumers that need runtime verification. */
size_t blst_p2_affine_sizeof(void);
/** Returns `sizeof(blst_fp12)` for ABI consumers that need runtime verification. */
size_t blst_fp12_sizeof(void);

#ifdef __cplusplus
}
#endif

#endif
