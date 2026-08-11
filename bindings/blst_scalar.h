/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BLST_SCALAR_H
#define BLST_SCALAR_H

#include "blst_types.h"

#ifdef __cplusplus
extern "C" {
# define BLST_CXX_DEFAULT = 0
#else
# define BLST_CXX_DEFAULT
#endif

typedef struct blst_pairing_context blst_pairing;
typedef struct blst_uniq_context blst_uniq;

/** Result of inserting a message into a `blst_uniq` tree. */
typedef enum BLST_UNIQ_RESULT {
    BLST_UNIQ_INSERTED = 0,
    BLST_UNIQ_DUPLICATE = 1,
    BLST_UNIQ_CAPACITY_EXCEEDED = 2,
    BLST_UNIQ_DEPTH_EXCEEDED = 3,
    BLST_UNIQ_INVALID_INPUT = 4,
} BLST_UNIQ_RESULT;
typedef char blst_uniq_result_must_be_4_bytes[
    sizeof(BLST_UNIQ_RESULT) == 4 ? 1 : -1];

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
 * Computes SHA-256 over `msg[0..msg_len)` and writes exactly 32 digest bytes.
 *
 * `out` must point to 32 writable bytes. `msg` must cover `msg_len` readable
 * bytes and may be null only when `msg_len == 0`. The message must contain
 * fewer than 2^64 bits. Output may overlap the input because all message bytes
 * are consumed before the digest is published. This operation is
 * allocation-free and nonthrowing, and erases its locally owned context before
 * returning. Its control flow depends on public message length, not contents.
 */
void blst_sha256(blst_byte out[32], const void* msg, size_t msg_len);

/**
 * Scalar byte helpers.
 *
 * Invariants:
 * - `blst_scalar` is 32-byte little-endian storage and may contain any
 *   256-bit integer. `blst_scalar_fr_check` establishes `[0,r)` and
 *   `blst_sk_check` establishes the nonzero secret-key range `[1,r)`.
 *   `blst_fr` is the corresponding Montgomery-form arithmetic type.
 * - The fixed-width word and scalar/Fr conversions below permit disjoint or
 *   exactly equal input/output storage. Partial overlap is not supported.
 * - `blst_scalar_from_{le,be}_bytes` accept variable-length inputs.
 * - Zero-length input decodes to the zero scalar, clears `out`, and returns
 *   `false`.
 * - The checked arithmetic helpers return `false` when the resulting scalar is
 *   zero or otherwise rejected by the scalar validity contract.
 * - All pointers named by a call must be valid for the complete object or byte
 *   range described by its signature. These entry points allocate nothing and
 *   do not allow C++ exceptions to cross the C ABI.
 */
/** Converts eight conventional 32-bit words into scalar storage without
 * reduction.
 * Exact aliasing is allowed. */
void blst_scalar_from_uint32(blst_scalar* out, const uint32_t a[8]);
/** Converts scalar storage into eight conventional 32-bit words.
 * Exact aliasing is allowed. */
void blst_uint32_from_scalar(uint32_t out[8], const blst_scalar* a);
/** Converts four conventional 64-bit words into scalar storage without
 * reduction.
 * Exact aliasing is allowed. */
void blst_scalar_from_uint64(blst_scalar* out, const uint64_t a[4]);
/** Converts scalar storage into four conventional 64-bit words.
 * Exact aliasing is allowed. */
void blst_uint64_from_scalar(uint64_t out[4], const blst_scalar* a);
/** Imports an exact 32-byte big-endian scalar representation. */
void blst_scalar_from_bendian(blst_scalar* out, const blst_byte a[32]);
/** Exports an exact 32-byte big-endian scalar representation. */
void blst_bendian_from_scalar(blst_byte out[32], const blst_scalar* a);
/** Imports an exact 32-byte little-endian scalar representation. */
void blst_scalar_from_lendian(blst_scalar* out, const blst_byte a[32]);
/** Exports an exact 32-byte little-endian scalar representation. */
void blst_lendian_from_scalar(blst_byte out[32], const blst_scalar* a);
/** Returns `true` iff `a` encodes an integer in `[0,r)`. Zero is accepted. */
bool blst_scalar_fr_check(const blst_scalar* a);
/** Returns `true` iff `a` is a canonical secret key in `[1,r)`. */
bool blst_sk_check(const blst_scalar* a);
/** Computes `(a+b) mod r` for canonical secret-key inputs; returns `true` iff
 * the result is nonzero. `out` receives the modular result even when it is
 * zero. Inputs outside `[1,r)` are outside the supported contract. */
bool blst_sk_add_n_check(blst_scalar* out, const blst_scalar* a, const blst_scalar* b);
/** Computes `(a-b) mod r` for canonical secret-key inputs; returns `true` iff
 * the result is nonzero. `out` receives the modular result even when it is
 * zero. Inputs outside `[1,r)` are outside the supported contract. */
bool blst_sk_sub_n_check(blst_scalar* out, const blst_scalar* a, const blst_scalar* b);
/** Computes `(a*b) mod r` for canonical secret-key inputs; returns `true` iff
 * the result is nonzero. `out` receives the modular result even when it is
 * zero. Inputs outside `[1,r)` are outside the supported contract. */
bool blst_sk_mul_n_check(blst_scalar* out, const blst_scalar* a, const blst_scalar* b);
/** Computes the modular inverse of a canonical nonzero secret key.
 * `a` must lie in `[1,r)`; the function does not report invalid input. */
void blst_sk_inverse(blst_scalar* out, const blst_scalar* a);
/** Reduces a variable-length little-endian integer modulo r. Returns `true`
 * iff the result is a valid nonzero secret key. Empty input writes zero and
 * returns `false`; otherwise `in` must reference `len` readable bytes. */
bool blst_scalar_from_le_bytes(blst_scalar* out, const blst_byte* in, size_t len);
/** Reduces a variable-length big-endian integer modulo r. Returns `true` iff
 * the result is a valid nonzero secret key. Empty input writes zero and returns
 * `false`; otherwise `in` must reference `len` readable bytes. */
bool blst_scalar_from_be_bytes(blst_scalar* out, const blst_byte* in, size_t len);
/** Strictly parses a complete hexadecimal string into a scalar. */
BLST_NODISCARD bool blst_scalar_from_hexascii(blst_scalar* out, const blst_byte* hex);

/**
 * Arithmetic over Fr, the scalar field modulo r.
 *
 * Inputs are reduced Montgomery-form values unless a conversion explicitly
 * says otherwise. Ordinary arithmetic produces a reduced Montgomery-form
 * value, permits exact output aliasing with an input, and forbids partial
 * overlap. The functions allocate nothing and are nonthrowing. Euclidean
 * naming is retained for compatibility, but both Fr inverse entry points use
 * the selected constant-time inverse. Other operations in this family are
 * intended to be constant-time with respect to field values unless documented
 * otherwise.
 */
/** Adds two Montgomery-form Fr values and writes the sum to `ret`. */
void blst_fr_add(blst_fr* ret, const blst_fr* a, const blst_fr* b);
/** Subtracts `b` from `a` in Fr and writes the difference to `ret`. */
void blst_fr_sub(blst_fr* ret, const blst_fr* a, const blst_fr* b);
/** Multiplies an Fr value by the integer constant 3. */
void blst_fr_mul_by_3(blst_fr* ret, const blst_fr* a);
/** Multiplies an Fr value by `2^count` modulo r; `count` is public and zero copies `a`. */
void blst_fr_lshift(blst_fr* ret, const blst_fr* a, size_t count);
/** Divides an Fr value by `2^count` modulo r; `count` is public and zero copies `a`. */
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
void blst_fr_pow(blst_fr* ret, const blst_fr* a, const blst_byte* exponent, size_t nbits);
/** Squares a Montgomery-form Fr value and writes the result to `ret`. */
void blst_fr_sqr(blst_fr* ret, const blst_fr* a);
/** Constant-time conditional negation: negates `a` when `flag` is true and
 * otherwise copies it into `ret`. */
void blst_fr_cneg(blst_fr* ret, const blst_fr* a, bool flag);
/** Compatibility spelling for the selected constant-time inverse of nonzero
 * `a` in Fr. Zero is outside the supported domain and has no error result. */
void blst_fr_eucl_inverse(blst_fr* ret, const blst_fr* a);
/** Computes the constant-time multiplicative inverse of nonzero `a` in Fr.
 * Zero is outside the supported domain and has no error result. */
void blst_fr_inverse(blst_fr* ret, const blst_fr* a);
/** Converts a four-limb unsigned integer into Montgomery-form Fr. */
void blst_fr_from_uint64(blst_fr* ret, const uint64_t a[4]);
/** Converts a Montgomery-form Fr value into four unsigned 64-bit limbs. */
void blst_uint64_from_fr(uint64_t ret[4], const blst_fr* a);
/** Converts scalar storage in `[0,r)` into Montgomery-form Fr. Callers must
 * establish the range with `blst_scalar_fr_check`; no Release validation or
 * reduction is performed. Exact aliasing is allowed. */
void blst_fr_from_scalar(blst_fr* ret, const blst_scalar* a);
/** Converts reduced Montgomery-form Fr into canonical scalar bytes in `[0,r)`.
 * Exact aliasing is allowed. */
void blst_scalar_from_fr(blst_scalar* ret, const blst_fr* a);
/** Bulk-converts contiguous Montgomery-form field elements into canonical little-endian scalars. Empty input succeeds without storage; nonempty input requires both arrays. */
bool blst_scalars_from_fr(blst_scalar out[], const blst_fr scalars[], size_t npoints);
/** Applies the constant-time butterfly step used by the Fr FFT routines.
 * `x0` and `x1` must be distinct objects; `twiddle` may alias either. */
void blst_fr_ct_bfly(blst_fr* x0, blst_fr* x1, const blst_fr* twiddle);
/** Applies the constant-time Gentleman-Sande butterfly step. `x0` and `x1`
 * must be distinct objects; `twiddle` may alias either. */
void blst_fr_gs_bfly(blst_fr* x0, blst_fr* x1, const blst_fr* twiddle);
/** Converts an Fr value into Montgomery form and writes it to `ret`. */
void blst_fr_to(blst_fr* ret, const blst_fr* a);
/** Converts an Fr value out of Montgomery form and writes it to `ret`. */
void blst_fr_from(blst_fr* ret, const blst_fr* a);
/** Strictly parses a complete hexadecimal string into Montgomery-form Fr. */
BLST_NODISCARD bool blst_fr_from_hexascii(blst_fr* ret, const blst_byte* hex);

/**
 * Arithmetic over Fp. Integer conversions use the exact array widths shown in
 * the signatures. Big-endian and little-endian byte helpers always operate on
 * 48-byte field encodings. Inputs are reduced Montgomery-form values unless a
 * conversion explicitly says otherwise. Ordinary arithmetic produces a
 * reduced Montgomery-form value, permits exact output aliasing, forbids
 * partial overlap, allocates nothing, and is nonthrowing. Euclidean inversion
 * is variable-time; fixed-chain inversion and ordinary arithmetic are intended
 * to be constant-time with respect to field values.
 */
/** Adds two Montgomery-form Fp values and writes the sum to `ret`. */
void blst_fp_add(blst_fp* ret, const blst_fp* a, const blst_fp* b);
/** Subtracts `b` from `a` in Fp and writes the difference to `ret`. */
void blst_fp_sub(blst_fp* ret, const blst_fp* a, const blst_fp* b);
/** Multiplies an Fp value by the integer constant 3. */
void blst_fp_mul_by_3(blst_fp* ret, const blst_fp* a);
/** Multiplies an Fp value by the integer constant 8. */
void blst_fp_mul_by_8(blst_fp* ret, const blst_fp* a);
/** Multiplies an Fp value by `2^count` modulo p; `count` is public and zero copies `a`. */
void blst_fp_lshift(blst_fp* ret, const blst_fp* a, size_t count);
/** Multiplies two Montgomery-form Fp values and writes the product to `ret`. */
void blst_fp_mul(blst_fp* ret, const blst_fp* a, const blst_fp* b);
/** Squares a Montgomery-form Fp value and writes the result to `ret`. */
void blst_fp_sqr(blst_fp* ret, const blst_fp* a);
/** Constant-time conditional negation: negates `a` when `flag` is true and
 * otherwise copies it into `ret`. */
void blst_fp_cneg(blst_fp* ret, const blst_fp* a, bool flag);
/** Computes the variable-time Euclidean inverse of nonzero `a` in Fp.
 * Zero is outside the supported domain and has no error result. */
void blst_fp_eucl_inverse(blst_fp* ret, const blst_fp* a);
/** Computes the fixed-chain constant-time inverse of nonzero `a` in Fp.
 * Zero is outside the supported domain and has no error result. */
void blst_fp_inverse(blst_fp* ret, const blst_fp* a);
/** Attempts to compute an Fp square root; returns `true` and publishes a root
 * for quadratic residues, otherwise returns `false` and clears `ret`. */
bool blst_fp_sqrt(blst_fp* ret, const blst_fp* a);
/** Converts a 12-word unsigned integer into Montgomery-form Fp. */
void blst_fp_from_uint32(blst_fp* ret, const uint32_t a[12]);
/** Converts a Montgomery-form Fp value into 12 unsigned 32-bit words. */
void blst_uint32_from_fp(uint32_t ret[12], const blst_fp* a);
/** Converts a six-limb unsigned integer into Montgomery-form Fp. */
void blst_fp_from_uint64(blst_fp* ret, const uint64_t a[6]);
/** Converts a Montgomery-form Fp value into six unsigned 64-bit limbs. */
void blst_uint64_from_fp(uint64_t ret[6], const blst_fp* a);
/** Imports a 48-byte big-endian integer, reduces it modulo p, and writes
 * Montgomery-form Fp. This is not a canonical-decoding validator. */
void blst_fp_from_bendian(blst_fp* ret, const blst_byte a[48]);
/** Encodes a Montgomery-form Fp value into a 48-byte big-endian string. */
void blst_bendian_from_fp(blst_byte ret[48], const blst_fp* a);
/** Imports a 48-byte little-endian integer, reduces it modulo p, and writes
 * Montgomery-form Fp. This is not a canonical-decoding validator. */
void blst_fp_from_lendian(blst_fp* ret, const blst_byte a[48]);
/**
 * Reduces an arbitrary-length little-endian unsigned integer modulo Fp.
 * `in` must reference `len` readable bytes when `len` is nonzero.
 * When `len` is zero, `in` may be null and `ret` is set to field zero.
 */
void blst_fp_from_le_bytes(blst_fp* ret, const blst_byte* in, size_t len);
/**
 * Reduces an arbitrary-length big-endian unsigned integer modulo Fp.
 * `in` must reference `len` readable bytes when `len` is nonzero.
 * When `len` is zero, `in` may be null and `ret` is set to field zero.
 */
void blst_fp_from_be_bytes(blst_fp* ret, const blst_byte* in, size_t len);
/** Encodes a Montgomery-form Fp value into a 48-byte little-endian string. */
void blst_lendian_from_fp(blst_byte ret[48], const blst_fp* a);
/** Converts an Fp value into Montgomery form and writes it to `ret`. */
void blst_fp_to(blst_fp* ret, const blst_fp* a);
/** Converts an Fp value out of Montgomery form and writes it to `ret`. */
void blst_fp_from(blst_fp* ret, const blst_fp* a);
/** Returns `true` when `a` is a quadratic residue in Fp. */
bool blst_fp_is_square(const blst_fp* a);
/** Strictly parses a complete hexadecimal string into Montgomery-form Fp. */
BLST_NODISCARD bool blst_fp_from_hexascii(blst_fp* ret, const blst_byte* hex);

/**
 * Arithmetic over Fp2. Components are reduced Montgomery-form Fp values.
 * Ordinary arithmetic produces reduced Montgomery-form components, permits
 * exact output aliasing, forbids partial overlap, allocates nothing, and is
 * nonthrowing. Euclidean inversion is variable-time; fixed-chain inversion and
 * ordinary arithmetic are intended to be constant-time with respect to field
 * values.
 */
/** Adds two Montgomery-form Fp2 values and writes the sum to `ret`. */
void blst_fp2_add(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b);
/** Subtracts `b` from `a` in Fp2 and writes the difference to `ret`. */
void blst_fp2_sub(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b);
/** Multiplies an Fp2 value by the integer constant 3. */
void blst_fp2_mul_by_3(blst_fp2* ret, const blst_fp2* a);
/** Multiplies an Fp2 value by the integer constant 8. */
void blst_fp2_mul_by_8(blst_fp2* ret, const blst_fp2* a);
/** Multiplies an Fp2 value by `2^count` modulo p; `count` is public and zero copies `a`. */
void blst_fp2_lshift(blst_fp2* ret, const blst_fp2* a, size_t count);
/** Multiplies two Montgomery-form Fp2 values and writes the product to `ret`. */
void blst_fp2_mul(blst_fp2* ret, const blst_fp2* a, const blst_fp2* b);
/** Squares a Montgomery-form Fp2 value and writes the result to `ret`. */
void blst_fp2_sqr(blst_fp2* ret, const blst_fp2* a);
/** Constant-time conditional negation: negates `a` when `flag` is true and
 * otherwise copies it into `ret`. */
void blst_fp2_cneg(blst_fp2* ret, const blst_fp2* a, bool flag);
/** Computes the variable-time Euclidean inverse of nonzero `a` in Fp2.
 * Zero is outside the supported domain and has no error result. */
void blst_fp2_eucl_inverse(blst_fp2* ret, const blst_fp2* a);
/** Computes the fixed-chain constant-time inverse of nonzero `a` in Fp2.
 * Zero is outside the supported domain and has no error result. */
void blst_fp2_inverse(blst_fp2* ret, const blst_fp2* a);
/**
 * Attempts to compute an Fp2 square root and returns `true` only for quadratic
 * residues. On failure, `ret` contains an unspecified intermediate and must be
 * ignored. If `ret == a`, failure may overwrite the input.
 */
bool blst_fp2_sqrt(blst_fp2* ret, const blst_fp2* a);
/** Returns `true` when `a` is a quadratic residue in Fp2. */
bool blst_fp2_is_square(const blst_fp2* a);

/**
 * Arithmetic over Fp12 and GT-style helpers.
 *
 * Invariants:
 * - `blst_fp12_one()` returns a pointer to immutable static storage.
 * - `blst_final_exp` consumes one nonzero Miller-loop/general-Fp12 value and
 *   produces a GT element.
 * - `blst_fp12_finalverify` consumes two trusted Miller-loop contributions and
 *   tests whether their final-exponentiation images are equal.
 * - General Fp12 arithmetic accepts reduced Montgomery-form coefficients and
 *   permits exact output aliasing. Partial overlap is forbidden.
 * - Functions in this family allocate nothing and are nonthrowing.
 */
/** Multiplies two Fp12 values and writes the product to `ret`. */
void blst_fp12_mul(blst_fp12* ret, const blst_fp12* a, const blst_fp12* b);
/** Squares an Fp12 value and writes the result to `ret`. */
void blst_fp12_sqr(blst_fp12* ret, const blst_fp12* a);
/** Applies the cyclotomic squaring shortcut. `a` must lie in the cyclotomic
 * subgroup; Release does not validate this precondition. Use `blst_fp12_sqr`
 * for a general Fp12 value. */
void blst_fp12_cyclotomic_sqr(blst_fp12* ret, const blst_fp12* a);
/**
 * Multiply by packed `{x,y,z}` representing the sparse Fp12 element
 * `((x,y,0),(0,z,0))`. `xy00z0` is not an ordinary Fp6 multiplier.
 * Exact `ret == a` aliasing is supported; partial overlap is not.
 */
void blst_fp12_mul_by_xy00z0(blst_fp12* ret, const blst_fp12* a, const blst_fp6* xy00z0);
/** Replaces `a` with its Fp12 conjugate in place. */
void blst_fp12_conjugate(blst_fp12* a);
/** Computes the multiplicative inverse of nonzero `a` in Fp12. Zero is outside
 * the supported domain and has no error result. */
void blst_fp12_inverse(blst_fp12* ret, const blst_fp12* a);
/**
 * Applies the `n`th Frobenius map to `a` and writes the result to `ret`.
 * Only the powers 1, 2, and 3 used by CatBLST are supported. Returns `false`
 * and leaves `ret` unchanged for every other value of `n`.
 */
BLST_NODISCARD bool blst_fp12_frobenius_map(blst_fp12* ret, const blst_fp12* a, size_t n);
/** Returns `true` when the two Fp12 values are equal. */
bool blst_fp12_is_equal(const blst_fp12* a, const blst_fp12* b);
/** Returns `true` when `a` is the multiplicative identity. */
bool blst_fp12_is_one(const blst_fp12* a);
/** Returns `true` when `a` lies in the target subgroup used by pairings. */
bool blst_fp12_in_group(const blst_fp12* a);
/** Returns a pointer to immutable static storage holding the Fp12 identity. */
const blst_fp12* blst_fp12_one(void);
/** Applies pairing final exponentiation to a nonzero general Fp12/Miller-loop
 * value and produces a GT element. `ret` and `f` must be non-null; exact
 * aliasing is supported. Zero is outside the supported domain and is not
 * reported through a status result. */
void blst_final_exp(blst_fp12* ret, const blst_fp12* f);
/** Returns immutable static storage holding the additive GT identity. */
const blst_fp12* blst_gt_identity(void);
/** Returns immutable static storage holding pairing(G1 generator, G2 generator). */
const blst_fp12* blst_gt_generator(void);
/** Returns `true` when target-group elements `a` and `b` are equal. */
bool blst_gt_is_equal(const blst_fp12* a, const blst_fp12* b);
/** Returns `true` when target-group element `p` is the additive identity. */
bool blst_gt_is_identity(const blst_fp12* p);
/** Adds target-group elements `a` and `b`. Both inputs must already lie in GT;
 * Release does not validate them. Inputs and output may exactly alias. */
void blst_gt_add(blst_fp12* out, const blst_fp12* a, const blst_fp12* b);
/** Subtracts target-group element `b` from `a`. Both inputs must already lie in
 * GT; Release does not validate them. Inputs and output may exactly alias. */
void blst_gt_sub(blst_fp12* out, const blst_fp12* a, const blst_fp12* b);
/** Negates target-group element `p`. The input must already lie in GT; Release
 * does not validate it. Input and output may exactly alias. */
void blst_gt_neg(blst_fp12* out, const blst_fp12* p);
/** Doubles target-group element `p`. The input must already lie in GT; Release
 * does not validate it. Input and output may exactly alias. */
void blst_gt_double(blst_fp12* out, const blst_fp12* p);
/**
 * Multiplies target-group element `p` by a little-endian scalar.
 *
 * `p` must already lie in GT and `nbits` must not exceed 256. The scalar need
 * not be reduced modulo the group order. The operation is cyclotomic
 * exponentiation; it does not validate subgroup membership.
 *
 * For canonical scalars in `[0,r)`, runtime and memory access are independent
 * of scalar contents for a fixed public `nbits`. A noncanonical scalar is
 * accepted and reduced modulo r, but timing may reveal that reduction was
 * required. Exact `out == p` aliasing is supported; partial overlap is not.
 */
void blst_gt_mult(blst_fp12* out, const blst_fp12* p, const blst_byte* scalar, size_t nbits);
/** Optimized variable-time GT multiplication for public scalar data only.
 * Domain, scalar reduction, bounds, and aliasing match `blst_gt_mult`; Release
 * does not validate GT membership. */
void blst_gt_mult_vartime(blst_fp12* out, const blst_fp12* p, const blst_byte* scalar, size_t nbits);
/**
 * Returns `true` when trusted Miller-loop contributions `gt1` and `gt2` have
 * equal images after final exponentiation. The historical parameter names and
 * `blst_fp12` storage do not make the inputs validated GT elements; no input
 * validation is performed. Both pointers must be non-null.
 */
bool blst_fp12_finalverify(const blst_fp12* gt1, const blst_fp12* gt2);
/** Encodes an Fp12 value into a 576-byte big-endian string. */
void blst_bendian_from_fp12(blst_byte out[48 * 12], const blst_fp12* a);

/**
 * G1 projective and affine helpers.
 *
 * Parameter interpretation:
 * - Scalar multiplication reads `scalar` as an unsigned little-endian bitstring
 *   with exactly `nbits` meaningful bits. `nbits` must not exceed 256, and
 *   `scalar` may be null only when `nbits` is zero.
 * - The unsuffixed single-point multiply is constant-time for fixed public
 *   `nbits`; `_vartime` and all multi-scalar/Pippenger operations require
 *   public scalars.
 * - Serialization writes 96-byte uncompressed or 48-byte compressed encodings.
 * - Generator accessors return immutable static storage.
 * - Basic single-point arithmetic, conversion, observation, and encoding accept
 *   on-curve points; subgroup membership is not required unless stated.
 * - Batch helpers, optimized multiplication, MSM, pairing, and protocol
 *   operations require previously validated subgroup points.
 */
/** Adds on-curve G1 projective points using a fast incomplete formula, with infinity
 * supported. For finite inputs, "opposite" means the same affine x coordinate
 * and negated affine y coordinate, so their mathematical sum is infinity.
 * Equal or opposite finite inputs are outside this function's
 * contract; use `blst_p1_add_or_double` when those cases are possible.
 * Subgroup membership is not required or validated. */
void blst_p1_add(blst_p1* out, const blst_p1* a, const blst_p1* b);
/** Subtracts on-curve G1 projective point `b` from `a`. Subgroup membership is
 * not required or validated; `out` remains on-curve. */
void blst_p1_sub(blst_p1* out, const blst_p1* a, const blst_p1* b);
/** Complete addition of on-curve G1 points: doubles equal inputs and returns
 * infinity for opposites. Subgroup membership is not required or validated. */
void blst_p1_add_or_double(blst_p1* out, const blst_p1* a, const blst_p1* b);
/** Adds on-curve projective and affine G1 points using a fast incomplete formula, with
 * infinity supported. Finite inputs must be distinct and not opposites; use
 * `blst_p1_add_or_double_affine` when either case is possible. Subgroup
 * membership is not required or validated. */
void blst_p1_add_affine(blst_p1* out, const blst_p1* a, const blst_p1_affine* b);
/** Subtracts on-curve affine G1 point `b` from on-curve projective `a`.
 * Subgroup membership is not required or validated. */
void blst_p1_sub_affine(blst_p1* out, const blst_p1* a, const blst_p1_affine* b);
/** Complete mixed addition of on-curve G1 points, including equal and opposite
 * finite inputs. Subgroup membership is not required or validated. */
void blst_p1_add_or_double_affine(blst_p1* out, const blst_p1* a, const blst_p1_affine* b);
/** Doubles on-curve G1 projective point `a`. Subgroup membership is not
 * required or validated; `out` remains on-curve. */
void blst_p1_double(blst_p1* out, const blst_p1* a);
/** Multiplies subgroup point `p` by `scalar mod r`. For canonical scalars and fixed public `nbits`, the operation is constant-time. Noncanonical scalars are accepted, but timing may reveal that reduction was required. Release does not validate `p`. Exact `out == p` aliasing is supported; partial overlap is not. */
void blst_p1_mult(blst_p1* out, const blst_p1* p, const blst_byte* scalar, size_t nbits);
/** Optimized variable-time G1 multiplication for public scalar data only. `p` must already be an on-curve prime-order subgroup point; scalars are interpreted modulo r and Release does not validate `p`. */
void blst_p1_mult_vartime(blst_p1* out, const blst_p1* p, const blst_byte* scalar, size_t nbits);
/** Conditionally negates an on-curve point in constant time. Subgroup
 * membership is not required or validated. */
void blst_p1_cneg(blst_p1* p, bool cbit);
/** Conditionally negates an on-curve affine point in constant time. Subgroup
 * membership is not required or validated. */
void blst_p1_affine_cneg(blst_p1_affine* p, bool cbit);
/** Converts an on-curve projective point into affine form. Subgroup membership
 * is preserved but neither required nor established. */
void blst_p1_to_affine(blst_p1_affine* out, const blst_p1* in);
/** Converts an on-curve affine point into projective form. Subgroup membership
 * is preserved but neither required nor established. */
void blst_p1_from_affine(blst_p1* out, const blst_p1_affine* in);
/** Normalizes an on-curve Jacobian point. Subgroup membership is preserved but
 * neither required nor established. */
void blst_p1_from_jacobian(blst_p1* out, const blst_p1* in);
/** Returns `true` when `p` lies on the G1 curve equation. */
bool blst_p1_on_curve(const blst_p1* p);
/** Returns `true` when `p` is on-curve and lies in the prime-order G1 subgroup. */
bool blst_p1_in_g1(const blst_p1* p);
/** Compares two on-curve projective points. Subgroup membership is not required;
 * this is not a validator. */
bool blst_p1_is_equal(const blst_p1* a, const blst_p1* b);
/** Tests whether an on-curve point is infinity. Subgroup membership is not
 * required; this is not a validator. */
bool blst_p1_is_inf(const blst_p1* a);
/** Returns a pointer to immutable static storage holding the G1 generator. */
const blst_p1* blst_p1_generator(void);
/** Legacy upstream `unchecked_mult`: multiplies an on-curve `p` by the exact scalar in constant time without subgroup endomorphisms. `p` need not be in G1; Release does not validate curve membership. */
void blst_p1_unchecked_mult(blst_p1* out, const blst_p1* p, const blst_byte* scalar, size_t nbits);
/** Returns `true` when affine G1 point `p` lies on the curve. */
bool blst_p1_affine_on_curve(const blst_p1_affine* p);
/** Returns `true` when affine `p` is on-curve and lies in the prime-order G1 subgroup. */
bool blst_p1_affine_in_g1(const blst_p1_affine* p);
/** Compares two on-curve affine points. Subgroup membership is not required;
 * this is not a validator. */
bool blst_p1_affine_is_equal(const blst_p1_affine* a, const blst_p1_affine* b);
/** Tests whether an on-curve affine point is infinity. Subgroup membership is
 * not required; this is not a validator. */
bool blst_p1_affine_is_inf(const blst_p1_affine* a);
/** Returns a pointer to immutable static storage holding the affine G1 generator. */
const blst_p1_affine* blst_p1_affine_generator(void);
/** Batch-converts `npoints` subgroup G1 projective inputs into affine outputs.
 * This advanced helper does not promise non-subgroup behavior. When `npoints >
 * 0`, `dst` and `points[0]` must be non-null; later null descriptor entries
 * mean contiguous continuation from the preceding point. Input and output
 * storage must not overlap. Empty input performs no access. */
void blst_p1s_to_affine(blst_p1_affine dst[], const blst_p1* const points[], size_t npoints);
/** Variable-time addition of public affine G1 subgroup points. Writes the complete sum to `ret`; empty or invalid-head input writes infinity. When `npoints > 0`, `points[0]` must be non-null and later null entries mean contiguous continuation from the previous non-null element. Exact output/input overlap is unsupported. A bounded internal allocation accelerates batches of at least 16 points; allocation failure transparently uses complete serial addition. */
void blst_p1s_add(blst_p1* ret, const blst_p1_affine* const points[], size_t npoints);
/** Returns the table size needed for G1 wbits precomputation, or 0 when `wbits` is outside the supported range `[2, 14]` or the size computation would overflow. */
size_t blst_p1s_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints);
/** Precomputes a reusable G1 table of `npoints * 2^(wbits-1)` affine entries. `wbits` must lie in `[2, 14]`; for nonempty input, `table` and `points[0]` must be non-null. Later null descriptor entries mean contiguous continuation. A packed input array may occupy the final `npoints` entries of `table`; other overlap is unsupported. The routine first reuses free table storage, then may make one nonthrowing allocation capped at 144 KiB, and falls back to a fixed 8 KiB stack workspace if allocation is unavailable. */
void blst_p1s_mult_wbits_precompute(blst_p1_affine table[], size_t wbits, const blst_p1_affine* const points[], size_t npoints);
/** Returns the caller-owned scratch size needed by G1 fixed-table MSM. The result is zero for empty input and is capped at the implementation's batch capacity. */
size_t blst_p1s_mult_wbits_scratch_sizeof(size_t npoints);
/** Evaluates variable-time subgroup G1 MSM from a matching wbits table. `wbits` must lie in `[2, 14]`, `nbits` in `[1, 256]`, and nonempty input requires non-null table, scalar head, and queried scratch. Scalars use `ceil(nbits/8)` little-endian bytes; later null descriptors mean contiguous continuation. Empty input or zero bits writes infinity. Invalid parameters also fail closed to infinity. */
void blst_p1s_mult_wbits(blst_p1* ret, const blst_p1_affine table[], size_t wbits, size_t npoints, const blst_byte* const scalars[], size_t nbits, blst_limb_t* scratch);
/** Returns scratch bytes for the automatically selected G1 Pippenger window. */
size_t blst_p1s_mult_pippenger_scratch_sizeof(size_t npoints);
/** Returns G1 Pippenger scratch bytes for an explicit window in `[1, 24]`, or 0 when the window or size is invalid. */
size_t blst_p1s_mult_pippenger_window_scratch_sizeof(size_t window);
/** Evaluates variable-time subgroup G1 MSM with Pippenger. Points must already be prime-order G1 subgroup elements; scalars are interpreted modulo r and contain at most 256 meaningful little-endian bits. Descriptor arrays may contain exactly `npoints` entries: later null entries mean contiguous continuation and no sentinel is read. The one-point fallback uses `mult_vartime` and needs no scratch; multiple points require queried scratch. Empty input or zero bits writes infinity. Release does not validate membership. */
void blst_p1s_mult_pippenger(blst_p1* ret, const blst_p1_affine* const points[], size_t npoints, const blst_byte* const scalars[], size_t nbits, blst_limb_t* scratch);
/** Returns total scratch bytes for a contiguous G1 MSM from native Montgomery-form `blst_fr` scalars, including canonical scalar preparation and Pippenger buckets. Zero is valid for an empty MSM and otherwise reports an unrepresentable size. */
size_t blst_p1s_mult_pippenger_fr_scratch_sizeof(size_t npoints);
/** Evaluates a variable-time contiguous G1 MSM from public, reduced
 * Montgomery-form `blst_fr` scalars and prevalidated subgroup points. For
 * nonempty input, `ret`, both input arrays, and limb-aligned scratch must be
 * mutually disjoint; scratch must contain at least the queried byte count.
 * Empty input writes infinity without reading input or scratch. Any other
 * `false` result leaves `ret` unchanged. No Release membership/range check or
 * allocation is performed. */
bool blst_p1s_mult_pippenger_fr(blst_p1* ret, const blst_p1_affine points[], const blst_fr scalars[], size_t npoints, blst_limb_t* scratch, size_t scratch_size);
/** Returns caller-owned scratch bytes for the selected contiguous G1 MSM implementation. This is a capacity query: storage for `max_points` supports every input length through that capacity. */
size_t blst_p1s_mult_fr_scratch_sizeof(size_t max_points);
/** Evaluates the same public-scalar G1 MSM using the selected scalar-backend
 * policy. The disjointness, alignment, failure-atomicity, and empty-input
 * contract is the same as `blst_p1s_mult_pippenger_fr`. Identity points are
 * accepted and select the Pippenger fallback. */
bool blst_p1s_mult_fr(blst_p1* ret, const blst_p1_affine points[], const blst_fr scalars[], size_t npoints, blst_limb_t* scratch, size_t scratch_size);
/** Returns total scratch bytes for the explicit G1 batch-affine MSM, including canonical scalar preparation. Supported batch sizes are 128 and 384; `window` must lie in `[1, 24]`, and `npoints` must fit `uint32_t`. Zero denotes either empty input or invalid/unrepresentable parameters. */
size_t blst_p1s_mult_batch_affine_fr_scratch_sizeof(size_t npoints, size_t window, size_t batch_size);
/** Evaluates a contiguous G1 MSM with the explicit variable-time batch-affine algorithm. Inputs are public, prevalidated, non-identity subgroup points and valid native Montgomery-form scalars. `ret`, points, scalars, and scratch must be mutually disjoint; scratch must be limb-aligned and contain at least the queried byte count. The routine allocates nothing and performs no Release point validation. Empty input succeeds and writes infinity; any other failure leaves `ret` unchanged. */
bool blst_p1s_mult_batch_affine_fr(blst_p1* ret, const blst_p1_affine points[], const blst_fr scalars[], size_t npoints, blst_limb_t* scratch, size_t scratch_size, size_t window, size_t batch_size);
/** Same contract as `blst_p1s_mult_pippenger`, using an explicit window in `[1, 24]`; scratch must match `blst_p1s_mult_pippenger_window_scratch_sizeof(window)`. */
void blst_p1s_mult_pippenger_window(blst_p1* ret, const blst_p1_affine* const points[], size_t npoints, const blst_byte* const scalars[], size_t nbits, blst_limb_t* scratch, size_t window);
/** Evaluates the signed Pippenger digit beginning at the exact scalar bit offset `bit0`. Inputs are public subgroup data; `nbits` is in `[1, 256]`, `window` in `[1, 24]`, and `bit0 < nbits`. Nonempty input requires queried scratch and valid descriptor heads. Later null entries mean contiguous continuation, and exact descriptor arrays need no sentinel. Invalid input writes infinity. */
void blst_p1s_tile_pippenger(blst_p1* ret, const blst_p1_affine* const points[], size_t npoints, const blst_byte* const scalars[], size_t nbits, blst_limb_t* scratch, size_t bit0, size_t window);
/** Serializes an on-curve G1 projective point into 96 bytes. Subgroup membership
 * is not required or validated. */
void blst_p1_serialize(blst_byte out[96], const blst_p1* in);
/** Compresses an on-curve G1 projective point into 48 bytes. Subgroup membership
 * is not required or validated. */
void blst_p1_compress(blst_byte out[48], const blst_p1* in);
/** Serializes an on-curve affine G1 point into 96 bytes. Subgroup membership is
 * not required or validated. */
void blst_p1_affine_serialize(blst_byte out[96], const blst_p1_affine* in);
/** Compresses an on-curve affine G1 point into 48 bytes. Subgroup membership is
 * not required or validated. */
void blst_p1_affine_compress(blst_byte out[48], const blst_p1_affine* in);
/** Decodes a 48-byte compressed G1 encoding. Success establishes canonical
 * encoding and on-curve status, but not subgroup membership; infinity is
 * accepted. Some known non-subgroup encodings may be rejected, but callers
 * must still use `blst_p1_affine_in_g1` or an approved batch validator before
 * subgroup-only operations. On failure `out` is unchanged. */
BLST_ERROR blst_p1_uncompress(blst_p1_affine* out, const blst_byte in[48]);
/** Decodes either a 48-byte compressed or 96-byte uncompressed G1 encoding,
 * selected by the encoding flags. The caller must provide the corresponding
 * readable size. Success establishes canonical encoding and on-curve status,
 * but not subgroup membership; infinity is accepted. Some known non-subgroup
 * encodings may be rejected, but callers must still use
 * `blst_p1_affine_in_g1` or an approved batch validator before subgroup-only
 * operations. On failure `out` is unchanged. */
BLST_ERROR blst_p1_deserialize(blst_p1_affine* out, const blst_byte in[96]);

/**
 * G2 projective and affine helpers.
 *
 * Parameter interpretation matches the G1 section, with 192-byte uncompressed
 * and 96-byte compressed encodings. The unsuffixed single-point multiply is
 * constant-time for fixed public `nbits`; `_vartime` and all
 * multi-scalar/Pippenger operations require public scalars.
 * Basic single-point arithmetic, conversion, observation, and encoding accept
 * on-curve points; subgroup membership is not required unless stated. Batch
 * and other advanced operations require previously validated subgroup points.
 */
/** Adds on-curve G2 projective points using a fast incomplete formula, with infinity
 * supported. For finite inputs, "opposite" means the same affine x coordinate
 * and negated affine y coordinate, so their mathematical sum is infinity.
 * Equal or opposite finite inputs are outside this function's
 * contract; use `blst_p2_add_or_double` when those cases are possible.
 * Subgroup membership is not required or validated. */
void blst_p2_add(blst_p2* out, const blst_p2* a, const blst_p2* b);
/** Subtracts on-curve G2 projective point `b` from `a`. Subgroup membership is
 * not required or validated; `out` remains on-curve. */
void blst_p2_sub(blst_p2* out, const blst_p2* a, const blst_p2* b);
/** Complete addition of on-curve G2 points: doubles equal inputs and returns
 * infinity for opposites. Subgroup membership is not required or validated. */
void blst_p2_add_or_double(blst_p2* out, const blst_p2* a, const blst_p2* b);
/** Adds on-curve projective and affine G2 points using a fast incomplete formula, with
 * infinity supported. Finite inputs must be distinct and not opposites; use
 * `blst_p2_add_or_double_affine` when either case is possible. Subgroup
 * membership is not required or validated. */
void blst_p2_add_affine(blst_p2* out, const blst_p2* a, const blst_p2_affine* b);
/** Subtracts on-curve affine G2 point `b` from on-curve projective `a`.
 * Subgroup membership is not required or validated. */
void blst_p2_sub_affine(blst_p2* out, const blst_p2* a, const blst_p2_affine* b);
/** Complete mixed addition of on-curve G2 points, including equal and opposite
 * finite inputs. Subgroup membership is not required or validated. */
void blst_p2_add_or_double_affine(blst_p2* out, const blst_p2* a, const blst_p2_affine* b);
/** Doubles on-curve G2 projective point `a`. Subgroup membership is not
 * required or validated; `out` remains on-curve. */
void blst_p2_double(blst_p2* out, const blst_p2* a);
/** Multiplies subgroup point `p` by `scalar mod r`. For canonical scalars and fixed public `nbits`, the operation is constant-time. Noncanonical scalars are accepted, but timing may reveal that reduction was required. Release does not validate `p`. Exact `out == p` aliasing is supported; partial overlap is not. */
void blst_p2_mult(blst_p2* out, const blst_p2* p, const blst_byte* scalar, size_t nbits);
/** Optimized variable-time G2 multiplication for public scalar data only. `p` must already be an on-curve prime-order subgroup point; scalars are interpreted modulo r and Release does not validate `p`. */
void blst_p2_mult_vartime(blst_p2* out, const blst_p2* p, const blst_byte* scalar, size_t nbits);
/** Conditionally negates an on-curve point in constant time. Subgroup
 * membership is not required or validated. */
void blst_p2_cneg(blst_p2* p, bool cbit);
/** Conditionally negates an on-curve affine point in constant time. Subgroup
 * membership is not required or validated. */
void blst_p2_affine_cneg(blst_p2_affine* p, bool cbit);
/** Converts an on-curve projective point into affine form. Subgroup membership
 * is preserved but neither required nor established. */
void blst_p2_to_affine(blst_p2_affine* out, const blst_p2* in);
/** Converts an on-curve affine point into projective form. Subgroup membership
 * is preserved but neither required nor established. */
void blst_p2_from_affine(blst_p2* out, const blst_p2_affine* in);
/** Normalizes an on-curve Jacobian point. Subgroup membership is preserved but
 * neither required nor established. */
void blst_p2_from_jacobian(blst_p2* out, const blst_p2* in);
/** Returns `true` when `p` lies on the G2 curve equation. */
bool blst_p2_on_curve(const blst_p2* p);
/** Returns `true` when `p` is on-curve and lies in the prime-order G2 subgroup. */
bool blst_p2_in_g2(const blst_p2* p);
/** Compares two on-curve projective points. Subgroup membership is not required;
 * this is not a validator. */
bool blst_p2_is_equal(const blst_p2* a, const blst_p2* b);
/** Tests whether an on-curve point is infinity. Subgroup membership is not
 * required; this is not a validator. */
bool blst_p2_is_inf(const blst_p2* a);
/** Returns a pointer to immutable static storage holding the G2 generator. */
const blst_p2* blst_p2_generator(void);
/** Legacy upstream `unchecked_mult`: multiplies an on-curve `p` by the exact scalar in constant time without subgroup endomorphisms. `p` need not be in G2; Release does not validate curve membership. */
void blst_p2_unchecked_mult(blst_p2* out, const blst_p2* p, const blst_byte* scalar, size_t nbits);
/** Returns `true` when affine G2 point `p` lies on the curve. */
bool blst_p2_affine_on_curve(const blst_p2_affine* p);
/** Returns `true` when affine `p` is on-curve and lies in the prime-order G2 subgroup. */
bool blst_p2_affine_in_g2(const blst_p2_affine* p);
/** Compares two on-curve affine points. Subgroup membership is not required;
 * this is not a validator. */
bool blst_p2_affine_is_equal(const blst_p2_affine* a, const blst_p2_affine* b);
/** Tests whether an on-curve affine point is infinity. Subgroup membership is
 * not required; this is not a validator. */
bool blst_p2_affine_is_inf(const blst_p2_affine* a);
/** Returns a pointer to immutable static storage holding the affine G2 generator. */
const blst_p2_affine* blst_p2_affine_generator(void);
/** Batch-converts `npoints` subgroup G2 projective inputs into affine outputs.
 * This advanced helper does not promise non-subgroup behavior. When `npoints >
 * 0`, `dst` and `points[0]` must be non-null; later null descriptor entries
 * mean contiguous continuation from the preceding point. Input and output
 * storage must not overlap. Empty input performs no access. */
void blst_p2s_to_affine(blst_p2_affine dst[], const blst_p2* const points[], size_t npoints);
/** Variable-time addition of public affine G2 subgroup points. Writes the complete sum to `ret`; empty or invalid-head input writes infinity. When `npoints > 0`, `points[0]` must be non-null and later null entries mean contiguous continuation from the previous non-null element. Exact output/input overlap is unsupported. A bounded internal allocation accelerates batches of at least 16 points; allocation failure transparently uses complete serial addition. */
void blst_p2s_add(blst_p2* ret, const blst_p2_affine* const points[], size_t npoints);
/** Returns the table size needed for G2 wbits precomputation, or 0 when `wbits` is outside the supported range `[2, 14]` or the size computation would overflow. */
size_t blst_p2s_mult_wbits_precompute_sizeof(size_t wbits, size_t npoints);
/** Precomputes a reusable G2 table of `npoints * 2^(wbits-1)` affine entries. `wbits` must lie in `[2, 14]`; for nonempty input, `table` and `points[0]` must be non-null. Later null descriptor entries mean contiguous continuation. A packed input array may occupy the final `npoints` entries of `table`; other overlap is unsupported. The routine first reuses free table storage, then may make one nonthrowing allocation capped at 144 KiB, and falls back to a fixed 8 KiB stack workspace if allocation is unavailable. */
void blst_p2s_mult_wbits_precompute(blst_p2_affine table[], size_t wbits, const blst_p2_affine* const points[], size_t npoints);
/** Returns the caller-owned scratch size needed by G2 fixed-table MSM. The result is zero for empty input and is capped at the implementation's batch capacity. */
size_t blst_p2s_mult_wbits_scratch_sizeof(size_t npoints);
/** Evaluates variable-time subgroup G2 MSM from a matching wbits table. `wbits` must lie in `[2, 14]`, `nbits` in `[1, 256]`, and nonempty input requires non-null table, scalar head, and queried scratch. Scalars use `ceil(nbits/8)` little-endian bytes; later null descriptors mean contiguous continuation. Empty input or zero bits writes infinity. Invalid parameters also fail closed to infinity. */
void blst_p2s_mult_wbits(blst_p2* ret, const blst_p2_affine table[], size_t wbits, size_t npoints, const blst_byte* const scalars[], size_t nbits, blst_limb_t* scratch);
/** Returns scratch bytes for the automatically selected G2 Pippenger window. */
size_t blst_p2s_mult_pippenger_scratch_sizeof(size_t npoints);
/** Returns G2 Pippenger scratch bytes for an explicit window in `[1, 24]`, or 0 when the window or size is invalid. */
size_t blst_p2s_mult_pippenger_window_scratch_sizeof(size_t window);
/** Evaluates variable-time subgroup G2 MSM with Pippenger. Points must already be prime-order G2 subgroup elements; scalars are interpreted modulo r and contain at most 256 meaningful little-endian bits. Descriptor arrays may contain exactly `npoints` entries: later null entries mean contiguous continuation and no sentinel is read. The one-point fallback uses `mult_vartime` and needs no scratch; multiple points require queried scratch. Empty input or zero bits writes infinity. Release does not validate membership. */
void blst_p2s_mult_pippenger(blst_p2* ret, const blst_p2_affine* const points[], size_t npoints, const blst_byte* const scalars[], size_t nbits, blst_limb_t* scratch);
/** Same contract as `blst_p2s_mult_pippenger`, using an explicit window in `[1, 24]`; scratch must match `blst_p2s_mult_pippenger_window_scratch_sizeof(window)`. */
void blst_p2s_mult_pippenger_window(blst_p2* ret, const blst_p2_affine* const points[], size_t npoints, const blst_byte* const scalars[], size_t nbits, blst_limb_t* scratch, size_t window);
/** Returns total scratch bytes for a contiguous G2 MSM from native Montgomery-form `blst_fr` scalars, including canonical scalar preparation and Pippenger buckets. Zero is valid for an empty MSM and otherwise reports an unrepresentable size. */
size_t blst_p2s_mult_pippenger_fr_scratch_sizeof(size_t npoints);
/** Evaluates a variable-time contiguous G2 MSM from public, reduced
 * Montgomery-form `blst_fr` scalars and prevalidated subgroup points. For
 * nonempty input, `ret`, both input arrays, and limb-aligned scratch must be
 * mutually disjoint; scratch must contain at least the queried byte count.
 * Empty input writes infinity without reading input or scratch. Any other
 * `false` result leaves `ret` unchanged. No Release membership/range check or
 * allocation is performed. */
bool blst_p2s_mult_pippenger_fr(blst_p2* ret, const blst_p2_affine points[], const blst_fr scalars[], size_t npoints, blst_limb_t* scratch, size_t scratch_size);
/** Returns caller-owned scratch bytes for the selected contiguous G2 MSM implementation. This is a capacity query: storage for `max_points` supports every input length through that capacity. */
size_t blst_p2s_mult_fr_scratch_sizeof(size_t max_points);
/** Evaluates the same public-scalar G2 MSM using the selected scalar-backend
 * policy. The disjointness, alignment, failure-atomicity, and empty-input
 * contract is the same as `blst_p2s_mult_pippenger_fr`. Identity points are
 * accepted and select the Pippenger fallback. */
bool blst_p2s_mult_fr(blst_p2* ret, const blst_p2_affine points[], const blst_fr scalars[], size_t npoints, blst_limb_t* scratch, size_t scratch_size);
/** Returns total scratch bytes for the explicit G2 batch-affine MSM, including canonical scalar preparation. Supported batch sizes are 128 and 384; `window` must lie in `[1, 24]`, and `npoints` must fit `uint32_t`. Zero denotes either empty input or invalid/unrepresentable parameters. */
size_t blst_p2s_mult_batch_affine_fr_scratch_sizeof(size_t npoints, size_t window, size_t batch_size);
/** Evaluates a contiguous G2 MSM with the explicit variable-time batch-affine algorithm. Inputs are public, prevalidated, non-identity subgroup points and valid native Montgomery-form scalars. `ret`, points, scalars, and scratch must be mutually disjoint; scratch must be limb-aligned and contain at least the queried byte count. The routine allocates nothing and performs no Release point validation. Empty input succeeds and writes infinity; any other failure leaves `ret` unchanged. */
bool blst_p2s_mult_batch_affine_fr(blst_p2* ret, const blst_p2_affine points[], const blst_fr scalars[], size_t npoints, blst_limb_t* scratch, size_t scratch_size, size_t window, size_t batch_size);
/** Evaluates the signed Pippenger digit beginning at the exact scalar bit offset `bit0`. Inputs are public subgroup data; `nbits` is in `[1, 256]`, `window` in `[1, 24]`, and `bit0 < nbits`. Nonempty input requires queried scratch and valid descriptor heads. Later null entries mean contiguous continuation, and exact descriptor arrays need no sentinel. Invalid input writes infinity. */
void blst_p2s_tile_pippenger(blst_p2* ret, const blst_p2_affine* const points[], size_t npoints, const blst_byte* const scalars[], size_t nbits, blst_limb_t* scratch, size_t bit0, size_t window);
/** Serializes an on-curve G2 projective point into 192 bytes. Subgroup
 * membership is not required or validated. */
void blst_p2_serialize(blst_byte out[192], const blst_p2* in);
/** Compresses an on-curve G2 projective point into 96 bytes. Subgroup
 * membership is not required or validated. */
void blst_p2_compress(blst_byte out[96], const blst_p2* in);
/** Serializes an on-curve affine G2 point into 192 bytes. Subgroup membership
 * is not required or validated. */
void blst_p2_affine_serialize(blst_byte out[192], const blst_p2_affine* in);
/** Compresses an on-curve affine G2 point into 96 bytes. Subgroup membership is
 * not required or validated. */
void blst_p2_affine_compress(blst_byte out[96], const blst_p2_affine* in);
/** Decodes a 96-byte compressed G2 encoding. Success establishes canonical
 * encoding and on-curve status, but not subgroup membership; infinity is
 * accepted. Some known non-subgroup encodings may be rejected, but callers
 * must still use `blst_p2_affine_in_g2` or an approved batch validator before
 * subgroup-only operations. On failure `out` is unchanged. */
BLST_ERROR blst_p2_uncompress(blst_p2_affine* out, const blst_byte in[96]);
/** Decodes either a 96-byte compressed or 192-byte uncompressed G2 encoding,
 * selected by the encoding flags. The caller must provide the corresponding
 * readable size. Success establishes canonical encoding and on-curve status,
 * but not subgroup membership; infinity is accepted. Some known non-subgroup
 * encodings may be rejected, but callers must still use
 * `blst_p2_affine_in_g2` or an approved batch validator before subgroup-only
 * operations. On failure `out` is unchanged. */
BLST_ERROR blst_p2_deserialize(blst_p2_affine* out, const blst_byte in[192]);

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
 * - A context has one orientation: PK-in-G1 with signatures in G2, or PK-in-G2
 *   with signatures in G1. Mixing orientations returns
 *   `BLST_AGGR_TYPE_MISMATCH` without changing the context.
 * - The ordinary aggregate entry points fully check a non-infinite signature,
 *   reject an infinite PK, and assume the PK was previously subgroup-checked.
 *   The `chk_n` variants let callers independently request both subgroup
 *   checks; a false flag transfers that precondition to the caller.
 * - `PK` and `signature` are independently nullable. A null PK contributes
 *   only the signature; a null signature contributes only the PK/message
 *   pairing. Supplying neither selects the orientation but adds no value.
 * - `msg`, `aug`, and `DST` may be null exactly when their supplied length is
 *   zero. Aggregate calls consume `msg` and `aug` before returning; the
 *   context retains only the borrowed `DST` described above.
 * - Call `blst_pairing_commit()` before final verification or merging. More
 *   items may be aggregated after a commit and committed again. Both inputs
 *   to `blst_pairing_merge()` must already have no pending pairs.
 * - Every `BLST_ERROR` aggregate operation is failure-atomic: a non-success
 *   result leaves the context unchanged.
 * - Pairing contexts are not internally synchronized. Parallel callers use
 *   one context per worker, commit each, and merge them afterward.
 * - These operations allocate no memory and do not throw across the C ABI.
 *
 * Normal flow for PKs in G1 and signatures in G2:
 *
 *     init(ctx, mode, DST)
 *       -> aggregate_pk_in_g1(ctx, pk[0], aggregate_sig, msg[0], ...)
 *       -> aggregate_pk_in_g1(ctx, pk[1], NULL,          msg[1], ...)
 *       -> commit(ctx)
 *       -> finalverify(ctx, NULL)
 *
 * To aggregate signatures separately, pass null signatures while accumulating
 * PK/message pairs, obtain a signature Miller contribution with the matching
 * G1/G2 helper, then pass it to final verification.
 * Raw pair aggregation and the Miller-result getter form the expert path for
 * constructing and exporting custom Miller products.
 */
/**
 * Computes a Miller loop from affine `(Q, P)` inputs and writes a Miller-loop
 * value, not a GT element, into `ret`. Non-null inputs must already be valid
 * prime-order subgroup points; no curve or subgroup validation is performed.
 * A null input selects the corresponding canonical generator. `ret` must be
 * non-null. Behavior for non-subgroup inputs is not promised.
 */
void blst_miller_loop(blst_fp12* ret, const blst_p2_affine* Q, const blst_p1_affine* P);
/**
 * Multiplies `n` Miller loops together from paired affine input arrays.
 *
 * Pointer-array contract:
 * - when `n == 0`, `ret` receives the Fp12 multiplicative identity
 * - when `n > 0`, `Qs[0]` and `Ps[0]` must be non-null
 * - a later null entry means "continue from the next contiguous element after
 *   the previous non-null entry"
 * - malformed nonempty array heads fail closed by writing the all-zero Fp12
 *   sentinel; every non-null point must be a previously validated subgroup
 *   point and no curve or subgroup validation is performed
 */
void blst_miller_loop_n(blst_fp12* ret, const blst_p2_affine* const Qs[], const blst_p1_affine* const Ps[], size_t n);
/** Precomputes the 68 line-function coefficients needed for repeated Miller
 * loops with previously validated subgroup point `Q`; when `Q` is null, the G2
 * generator is used. G2 infinity produces an identity table. No point
 * validation is performed. */
void blst_precompute_lines(blst_fp6 Qlines[68], const blst_p2_affine* Q);
/** Runs a Miller loop using a valid precomputed `Qlines` table and previously
 * validated subgroup point `P`, producing a Miller-loop value rather than GT;
 * when `P` is null, the G1 generator is used. G1 infinity returns the Miller
 * identity. No point validation is performed. */
void blst_miller_loop_lines(blst_fp12* ret, const blst_fp6 Qlines[68], const blst_p1_affine* P);
/** Returns the exact storage size, in bytes, required for one `blst_pairing` object. */
size_t blst_pairing_sizeof(void);
/** Returns the required alignment for `blst_pairing` storage. */
size_t blst_pairing_alignof(void);
/**
 * Initializes caller-owned, correctly sized and aligned pairing storage.
 * `hash_or_encode == true` selects hash-to-curve; false selects encode-to-curve.
 * Reinitialization discards all prior logical state.
 */
void blst_pairing_init(blst_pairing* new_ctx, bool hash_or_encode, const blst_byte* DST BLST_CXX_DEFAULT, size_t DST_len BLST_CXX_DEFAULT);
/** Returns the borrowed DST pointer currently stored in `ctx`. */
const blst_byte* blst_pairing_get_dst(const blst_pairing* ctx);
/** Flushes the pending batch into the context's Miller product; idempotent when no pairs are pending. */
void blst_pairing_commit(blst_pairing* ctx);
/**
 * Commits pending pairs and writes the accumulated Miller-loop product to
 * caller-owned `blst_fp12` storage. This is not a GT element. An empty context
 * produces the Fp12 multiplicative identity. `out` must not overlap `ctx`.
 */
void blst_pairing_get_miller_result(blst_fp12* out, blst_pairing* ctx);
/**
 * Queues one trusted affine `(G2 q, G1 p)` pair without curve or subgroup
 * validation. Both inputs must already be subgroup points; non-subgroup
 * behavior is not promised. Either infinity contributes the Miller identity.
 * This does not establish a PK/signature orientation.
 */
void blst_pairing_raw_aggregate(blst_pairing* ctx, const blst_p2_affine* q, const blst_p1_affine* p);
/**
 * Adds a PK-in-G2 verification item. A nonnull signature is fully subgroup
 * checked; a nonnull PK must already be a subgroup point and cannot be
 * infinity. Failure leaves `ctx` unchanged.
 */
BLST_ERROR blst_pairing_aggregate_pk_in_g2(blst_pairing* ctx, const blst_p2_affine* PK, const blst_p1_affine* signature, const blst_byte* msg, size_t msg_len, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/** Adds a PK-in-G2 item with independently selected PK and signature subgroup checks; false flags require caller-established subgroup membership. */
BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g2(blst_pairing* ctx, const blst_p2_affine* PK, bool pk_grpchk, const blst_p1_affine* signature, bool sig_grpchk, const blst_byte* msg, size_t msg_len, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/** Aggregates one PK-in-G2 item with a public, variable-time weight of at most 256 little-endian bits; malformed scalar arguments return `BLST_BAD_SCALAR` without mutation. */
BLST_ERROR blst_pairing_mul_n_aggregate_pk_in_g2(blst_pairing* ctx, const blst_p2_affine* PK, const blst_p1_affine* sig, const blst_byte* scalar, size_t nbits, const blst_byte* msg, size_t msg_len, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/** Adds a weighted PK-in-G2 item with explicit subgroup-check policy; the public scalar is reduced modulo r and malformed shape leaves `ctx` unchanged. */
BLST_ERROR blst_pairing_chk_n_mul_n_aggr_pk_in_g2(blst_pairing* ctx, const blst_p2_affine* PK, bool pk_grpchk, const blst_p1_affine* sig, bool sig_grpchk, const blst_byte* scalar, size_t nbits, const blst_byte* msg, size_t msg_len, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/**
 * Adds a PK-in-G1 verification item. A nonnull signature is fully subgroup
 * checked; a nonnull PK must already be a subgroup point and cannot be
 * infinity. Failure leaves `ctx` unchanged.
 */
BLST_ERROR blst_pairing_aggregate_pk_in_g1(blst_pairing* ctx, const blst_p1_affine* PK, const blst_p2_affine* signature, const blst_byte* msg, size_t msg_len, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/** Adds a PK-in-G1 item with independently selected PK and signature subgroup checks; false flags require caller-established subgroup membership. */
BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g1(blst_pairing* ctx, const blst_p1_affine* PK, bool pk_grpchk, const blst_p2_affine* signature, bool sig_grpchk, const blst_byte* msg, size_t msg_len, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/** Aggregates one PK-in-G1 item with a public, variable-time weight of at most 256 little-endian bits; malformed scalar arguments return `BLST_BAD_SCALAR` without mutation. */
BLST_ERROR blst_pairing_mul_n_aggregate_pk_in_g1(blst_pairing* ctx, const blst_p1_affine* PK, const blst_p2_affine* sig, const blst_byte* scalar, size_t nbits, const blst_byte* msg, size_t msg_len, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/** Adds a weighted PK-in-G1 item with explicit subgroup-check policy; the public scalar is reduced modulo r and malformed shape leaves `ctx` unchanged. */
BLST_ERROR blst_pairing_chk_n_mul_n_aggr_pk_in_g1(blst_pairing* ctx, const blst_p1_affine* PK, bool pk_grpchk, const blst_p2_affine* sig, bool sig_grpchk, const blst_byte* scalar, size_t nbits, const blst_byte* msg, size_t msg_len, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/** Merges committed `ctx1` into committed `ctx`; orientations must match when both are defined, `ctx1` is unchanged, and failure leaves `ctx` unchanged. */
BLST_ERROR blst_pairing_merge(blst_pairing* ctx, const blst_pairing* ctx1);
/**
 * Verifies a committed context. A nonnull `gtsig` is a trusted external Miller
 * contribution in `blst_fp12` storage. When it is null, the internally
 * aggregated signature is used, or the Miller identity if every supplied
 * signature was infinity. An empty PK/message product returns false.
 */
bool blst_pairing_finalverify(const blst_pairing* ctx, const blst_fp12* gtsig BLST_CXX_DEFAULT);
/** Decodes and subgroup-checks one wire-format G1 contribution, adding it to
 * optional `in`; exact `out == in` aliasing is supported, infinity preserves
 * `in`, and failure leaves `out` unchanged. `zwire` must be non-null and make
 * 48 bytes readable when its first byte has bit 7 set, otherwise 96 bytes.
 * Callers with uncertain flags should provide a 96-byte readable buffer. */
BLST_ERROR blst_aggregate_in_g1(blst_p1* out, const blst_p1* in, const blst_byte* zwire);
/** Decodes and subgroup-checks one wire-format G2 contribution, adding it to
 * optional `in`; exact `out == in` aliasing is supported, infinity preserves
 * `in`, and failure leaves `out` unchanged. `zwire` must be non-null and make
 * 96 bytes readable when its first byte has bit 7 set, otherwise 192 bytes.
 * Callers with uncertain flags should provide a 192-byte readable buffer. */
BLST_ERROR blst_aggregate_in_g2(blst_p2* out, const blst_p2* in, const blst_byte* zwire);
/** Computes the pre-final Miller contribution for a previously validated affine G1 subgroup signature. */
void blst_aggregated_in_g1(blst_fp12* out, const blst_p1_affine* signature);
/** Computes the pre-final Miller contribution for a previously validated affine G2 subgroup signature. */
void blst_aggregated_in_g2(blst_fp12* out, const blst_p2_affine* signature);
/** One-shot PK-in-G1 verification that rejects an infinite PK, fully checks both subgroup operands, and returns `BLST_VERIFY_FAIL` only when the final pairing equation fails. */
BLST_ERROR blst_core_verify_pk_in_g1(const blst_p1_affine* pk, const blst_p2_affine* signature, bool hash_or_encode, const blst_byte* msg, size_t msg_len, const blst_byte* DST BLST_CXX_DEFAULT, size_t DST_len BLST_CXX_DEFAULT, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/** One-shot PK-in-G2 verification that rejects an infinite PK, fully checks both subgroup operands, and returns `BLST_VERIFY_FAIL` only when the final pairing equation fails. */
BLST_ERROR blst_core_verify_pk_in_g2(const blst_p2_affine* pk, const blst_p1_affine* signature, bool hash_or_encode, const blst_byte* msg, size_t msg_len, const blst_byte* DST BLST_CXX_DEFAULT, size_t DST_len BLST_CXX_DEFAULT, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);

/**
 * Key derivation, hashing, mapping, and utility helpers.
 *
 * Invariants:
 * - For key generation, `out_SK` names writable scalar storage. `IKM` names
 *   `IKM_len` readable bytes when `IKM_len >= 32`; a shorter IKM is rejected
 *   by clearing `out_SK` and need not be readable. A null `info` retains
 *   upstream compatibility by ignoring `info_len`; otherwise it names exactly
 *   `info_len` readable bytes.
 * - Key generation is deterministic, allocation-free, and nonthrowing across
 *   the C ABI. It performs no Release validation beyond its documented
 *   length/salt cases. The v4-family rejection step repeats only if the
 *   derived scalar is zero.
 * - `blst_uniq_sizeof()` returns the exact caller allocation size for the
 *   opaque `blst_uniq` object, or zero if the capacity is invalid.
 * - `blst_expand_message_xmd` implements the XMD expand step used by the hash
 *   to curve suite.
 */
/**
 * Derives a canonical nonzero secret scalar using upstream BLST's default v4
 * scheme and default salt. The common key-generation contract above applies.
 */
void blst_keygen(blst_scalar* out_SK, const blst_byte* IKM, size_t IKM_len, const blst_byte* info BLST_CXX_DEFAULT, size_t info_len BLST_CXX_DEFAULT);
/**
 * Retains upstream BLST's legacy v3 derivation. A successful result is reduced
 * modulo r but, unlike v4 and later, the legacy algorithm does not retry the
 * negligible zero result. The common key-generation contract applies.
 */
void blst_keygen_v3(blst_scalar* out_SK, const blst_byte* IKM, size_t IKM_len, const blst_byte* info BLST_CXX_DEFAULT, size_t info_len BLST_CXX_DEFAULT);
/**
 * Derives a canonical nonzero secret scalar using the v4.5 interface. A
 * non-null `salt` names `salt_len` readable bytes, including a valid empty
 * salt. For upstream compatibility, null `salt` selects the default salt and
 * ignores `salt_len`. The common key-generation contract applies.
 */
void blst_keygen_v4_5(blst_scalar* out_SK, const blst_byte* IKM, size_t IKM_len, const blst_byte* salt, size_t salt_len, const blst_byte* info BLST_CXX_DEFAULT, size_t info_len BLST_CXX_DEFAULT);
/**
 * Derives a canonical nonzero secret scalar using the v5 interface. `salt`
 * must be non-null and name `salt_len` readable bytes; a non-null empty salt
 * is valid. A null salt clears `out_SK`. The common key-generation contract
 * applies.
 */
void blst_keygen_v5(blst_scalar* out_SK, const blst_byte* IKM, size_t IKM_len, const blst_byte* salt, size_t salt_len, const blst_byte* info BLST_CXX_DEFAULT, size_t info_len BLST_CXX_DEFAULT);
/**
 * Derives a canonical nonzero EIP-2333 master secret key using the default-v4
 * keygen operation with empty key info. The common IKM and output contracts
 * apply.
 */
void blst_derive_master_eip2333(blst_scalar* out_SK, const blst_byte* IKM, size_t IKM_len);
/**
 * Derives a canonical nonzero EIP-2333 child key. `SK` must name a canonical
 * nonzero parent secret scalar; no Release validation is performed. `out_SK`
 * and `SK` may be exactly the same object.
 */
void blst_derive_child_eip2333(blst_scalar* out_SK, const blst_scalar* SK, uint32_t child_index);
/**
 * Encodes `msg` to a deterministic G1 subgroup point with the RFC 9380
 * XMD:SHA-256 SSWU nonuniform map. `out` must be non-null. Each byte pointer
 * must reference its stated readable length; a pointer may be null when its
 * length is zero. `aug` is prepended to `msg`. The caller selects the protocol
 * ciphersuite through `DST`. No Release input validation is performed.
 * Allocation-free, nonthrowing, and safe across the C ABI for valid inputs.
 */
void blst_encode_to_g1(blst_p1* out, const blst_byte* msg, size_t msg_len, const blst_byte* DST BLST_CXX_DEFAULT, size_t DST_len BLST_CXX_DEFAULT, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/**
 * Hashes `msg` to a G1 subgroup point with the RFC 9380 XMD:SHA-256 SSWU
 * random-oracle map. Pointer, validation, allocation, and exception contracts
 * are the same as `blst_encode_to_g1`.
 */
void blst_hash_to_g1(blst_p1* out, const blst_byte* msg, size_t msg_len, const blst_byte* DST BLST_CXX_DEFAULT, size_t DST_len BLST_CXX_DEFAULT, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/**
 * Encodes `msg` to a deterministic G2 subgroup point with the RFC 9380
 * XMD:SHA-256 SSWU nonuniform map. Pointer, validation, allocation, and
 * exception contracts are the same as `blst_encode_to_g1`.
 */
void blst_encode_to_g2(blst_p2* out, const blst_byte* msg, size_t msg_len, const blst_byte* DST BLST_CXX_DEFAULT, size_t DST_len BLST_CXX_DEFAULT, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/**
 * Hashes `msg` to a G2 subgroup point with the RFC 9380 XMD:SHA-256 SSWU
 * random-oracle map. Pointer, validation, allocation, and exception contracts
 * are the same as `blst_encode_to_g1`.
 */
void blst_hash_to_g2(blst_p2* out, const blst_byte* msg, size_t msg_len, const blst_byte* DST BLST_CXX_DEFAULT, size_t DST_len BLST_CXX_DEFAULT, const blst_byte* aug BLST_CXX_DEFAULT, size_t aug_len BLST_CXX_DEFAULT);
/**
 * Maps reduced Montgomery-form field elements into a G1 subgroup point without
 * hashing. `out` and `u` must be non-null and disjoint. When `v` is null, maps
 * only `u` (the nonuniform encode form); otherwise `v` must name a second
 * reduced Montgomery element disjoint from `out`, and both mapped points are
 * added before cofactor clearing. No Release field validation is performed.
 * Allocation-free and nonthrowing for valid inputs.
 */
void blst_map_to_g1(blst_p1* out, const blst_fp* u, const blst_fp* v BLST_CXX_DEFAULT);
/**
 * Maps reduced Montgomery-form Fp2 elements into a G2 subgroup point without
 * hashing. Nullability, disjoint-storage, validation, allocation, and exception
 * contracts are the same as `blst_map_to_g1`.
 */
void blst_map_to_g2(blst_p2* out, const blst_fp2* u, const blst_fp2* v BLST_CXX_DEFAULT);
/** Derives the projective G1 public key for a canonical nonzero `SK` produced
 * by keygen or accepted by `blst_sk_check`. No Release validation occurs;
 * runtime is independent of a valid secret key's contents. */
void blst_sk_to_pk_in_g1(blst_p1* out_pk, const blst_scalar* SK);
/** Signs a previously validated G2 subgroup hash point with a canonical
 * nonzero `SK`. No Release point/key validation occurs; runtime is independent
 * of a valid secret key's contents. */
void blst_sign_pk_in_g1(blst_p2* out_sig, const blst_p2* hash, const blst_scalar* SK);
/** G2 counterpart of `blst_sk_to_pk_in_g1`. */
void blst_sk_to_pk_in_g2(blst_p2* out_pk, const blst_scalar* SK);
/** G1 counterpart of `blst_sign_pk_in_g1`. */
void blst_sign_pk_in_g2(blst_p1* out_sig, const blst_p1* hash, const blst_scalar* SK);
/**
 * Derives the affine G1 public key from `SK` and optionally writes a 96-byte
 * pre-decorated `X || Y` buffer. The buffer is not directly a conforming
 * encoding: set bit 7 and retain the first 48 bytes for compressed form, or
 * clear bit 5 for uncompressed form. Either output may be null. `SK` must be
 * canonical and nonzero; no Release validation occurs.
 */
void blst_sk_to_pk2_in_g1(blst_byte out[96], blst_p1_affine* out_pk, const blst_scalar* SK);
/**
 * Signs a G2 hash point and optionally writes the affine result and a 192-byte
 * pre-decorated `X || Y` buffer. Set bit 7 and retain the first 96 bytes for
 * compressed form, or clear bit 5 for uncompressed form. Either output may be
 * null; `hash` must already be in G2 and `SK` canonical and nonzero.
 */
void blst_sign_pk2_in_g1(blst_byte out[192], blst_p2_affine* out_sig, const blst_p2* hash, const blst_scalar* SK);
/** G2 counterpart of `blst_sk_to_pk2_in_g1`; `out` is 192 bytes and its
 * compressed x half is 96 bytes. */
void blst_sk_to_pk2_in_g2(blst_byte out[192], blst_p2_affine* out_pk, const blst_scalar* SK);
/** G1 counterpart of `blst_sign_pk2_in_g1`; `out` is 96 bytes and its
 * compressed x half is 48 bytes. */
void blst_sign_pk2_in_g2(blst_byte out[96], blst_p1_affine* out_sig, const blst_p1* hash, const blst_scalar* SK);
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
BLST_UNIQ_RESULT blst_uniq_insert(blst_uniq* tree, const blst_byte* msg, size_t len);
/**
 * Expands `msg` using RFC 9380 XMD:SHA-256 and writes exactly `out_len`
 * bytes. `out_len` must not exceed 8160 (`255 * 32`); an excessive length
 * clears all `out_len` bytes. `out` must reference `out_len` writable bytes
 * when the length is nonzero and may be null when it is zero. `msg` and `DST`
 * must reference their stated readable lengths and may be null only when the
 * corresponding length is zero; an invalid nonzero pointer/length pair clears
 * a valid output buffer. DSTs longer than 255 bytes are shortened as required
 * by RFC 9380. Input and output regions must not overlap. The operation is
 * allocation-free and nonthrowing; zero output length is a no-op.
 */
void blst_expand_message_xmd(blst_byte* out, size_t out_len, const blst_byte* msg, size_t msg_len, const blst_byte* DST, size_t DST_len);
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

#undef BLST_CXX_DEFAULT

#ifdef __cplusplus
}
#endif

#endif
