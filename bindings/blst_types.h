/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BLST_TYPES_H
#define BLST_TYPES_H

#include <stddef.h>
#include <stdint.h>
#include "blst/common/error_codes.h"

#if !defined(__cplusplus)
# if !defined(__STDC_VERSION__) || __STDC_VERSION__ < 199901L
#  error "CatBLST requires C99 or later"
# elif __STDC_VERSION__ < 202311L
#  include <stdbool.h>
# endif
#endif

#if defined(__cplusplus) && defined(__has_cpp_attribute)
# if __has_cpp_attribute(nodiscard)
#  define BLST_NODISCARD [[nodiscard]]
# endif
#elif !defined(__cplusplus) && defined(__has_c_attribute)
# if __has_c_attribute(nodiscard)
#  define BLST_NODISCARD [[nodiscard]]
# endif
#endif
#if !defined(BLST_NODISCARD)
# if defined(__GNUC__) || defined(__clang__)
#  define BLST_NODISCARD __attribute__((warn_unused_result))
# else
#  define BLST_NODISCARD
# endif
#endif

/** Library-owned byte type used by the public C API. */
typedef uint8_t blst_byte;
/** Library-owned 64-bit arithmetic limb type used by public scratch APIs. */
typedef uint64_t blst_limb_t;

#if defined(_MSC_VER)
# define BLST_AUX_ALIGNAS(n) __declspec(align(n))
#elif defined(__GNUC__) || defined(__clang__)
# define BLST_AUX_ALIGNAS(n) __attribute__((aligned(n)))
#else
# define BLST_AUX_ALIGNAS(n)
#endif

/**
 * Scalar storage encoded as a 32-byte little-endian integer.
 *
 * The type can hold any 256-bit bit pattern; the type alone does not establish
 * reduction modulo r or nonzero secret-key validity. Use
 * `blst_scalar_fr_check` for `[0,r)` and `blst_sk_check` for `[1,r)`. This is
 * the non-Montgomery representation used by secret-key and scalar byte APIs.
 * It represents the same mathematical field as `blst_fr`, but is not directly
 * usable by Fr arithmetic. Convert with `blst_fr_from_scalar` before
 * arithmetic and with `blst_scalar_from_fr` when returning to scalar storage.
 *
 * Invariants:
 * - The byte layout is ABI-stable and matches upstream BLST.
 * - Callers may store and move this value as plain POD bytes.
 * - A `blst_scalar` must never be passed to an API expecting Montgomery-form
 *   `blst_fr`, even though both representations occupy 32 bytes.
 */
typedef struct { blst_byte b[256 / 8]; } blst_scalar;

#ifdef __cplusplus

/*
 * The native C++ implementation owns the arithmetic and curve types. The C
 * names below are exact aliases, so C++ callers cross the C ABI without a
 * wrapper object, conversion, or copy.
 */
# include "blst/scalar/fields/Fr.h"
# include "blst/scalar/fields/Fp12.h"
# include "blst/scalar/curve/point.h"
# include <type_traits>

typedef blst::scalar::Fr blst_fr;
typedef blst::scalar::Fp blst_fp;
typedef blst::scalar::Fp2 blst_fp2;
typedef blst::scalar::Fp6 blst_fp6;
typedef blst::scalar::Fp12 blst_fp12;
typedef blst::scalar::Point<blst::scalar::Fp> blst_p1;
typedef blst::scalar::AffinePoint<blst::scalar::Fp> blst_p1_affine;
typedef blst::scalar::Point<blst::scalar::Fp2> blst_p2;
typedef blst::scalar::AffinePoint<blst::scalar::Fp2> blst_p2_affine;

#else

/** Scalar field element in Montgomery form. */
typedef struct { blst_limb_t l[256 / 8 / sizeof(blst_limb_t)]; } blst_fr;
/** Base field element in Montgomery form. */
typedef struct { blst_limb_t l[384 / 8 / sizeof(blst_limb_t)]; } blst_fp;
/** Quadratic extension field element in Montgomery form. */
typedef struct { blst_fp fp[2]; } blst_fp2;
/** Sextic extension field element in Montgomery form. */
typedef struct { blst_fp2 fp2[3]; } blst_fp6;
/** Twelfth-degree extension field element in Montgomery form. */
typedef struct { blst_fp6 fp6[2]; } blst_fp12;
/** G1 Jacobian point over Montgomery Fp coordinates. */
typedef struct { blst_fp x, y, z; } blst_p1;
/** G1 affine point over Montgomery Fp coordinates. */
typedef struct { blst_fp x, y; } blst_p1_affine;
/** G2 Jacobian point over Montgomery Fp2 coordinates. */
typedef struct { blst_fp2 x, y, z; } blst_p2;
/** G2 affine point over Montgomery Fp2 coordinates. */
typedef struct { blst_fp2 x, y; } blst_p2_affine;

#endif

#ifdef __cplusplus
static_assert(std::is_same<blst_fr, blst::scalar::Fr>::value);
static_assert(std::is_same<blst_fp, blst::scalar::Fp>::value);
static_assert(std::is_same<blst_fp2, blst::scalar::Fp2>::value);
static_assert(std::is_same<blst_fp6, blst::scalar::Fp6>::value);
static_assert(std::is_same<blst_fp12, blst::scalar::Fp12>::value);
static_assert(std::is_same<blst_p1, blst::scalar::Point<blst::scalar::Fp>>::value);
static_assert(std::is_same<blst_p1_affine, blst::scalar::AffinePoint<blst::scalar::Fp>>::value);
static_assert(std::is_same<blst_p2, blst::scalar::Point<blst::scalar::Fp2>>::value);
static_assert(std::is_same<blst_p2_affine, blst::scalar::AffinePoint<blst::scalar::Fp2>>::value);

static_assert(std::is_standard_layout<blst_fr>::value && std::is_trivially_copyable<blst_fr>::value);
static_assert(std::is_standard_layout<blst_fp>::value && std::is_trivially_copyable<blst_fp>::value);
static_assert(std::is_standard_layout<blst_fp2>::value && std::is_trivially_copyable<blst_fp2>::value);
static_assert(std::is_standard_layout<blst_fp6>::value && std::is_trivially_copyable<blst_fp6>::value);
static_assert(std::is_standard_layout<blst_fp12>::value && std::is_trivially_copyable<blst_fp12>::value);
static_assert(std::is_standard_layout<blst_p1>::value && std::is_trivially_copyable<blst_p1>::value);
static_assert(std::is_standard_layout<blst_p1_affine>::value && std::is_trivially_copyable<blst_p1_affine>::value);
static_assert(std::is_standard_layout<blst_p2>::value && std::is_trivially_copyable<blst_p2>::value);
static_assert(std::is_standard_layout<blst_p2_affine>::value && std::is_trivially_copyable<blst_p2_affine>::value);

static_assert(sizeof(blst_fr) == 32 && alignof(blst_fr) == alignof(blst_limb_t));
static_assert(sizeof(blst_fp) == 48 && alignof(blst_fp) == alignof(blst_limb_t));
static_assert(sizeof(blst_fp2) == 96 && alignof(blst_fp2) == alignof(blst_limb_t));
static_assert(sizeof(blst_fp6) == 288 && alignof(blst_fp6) == alignof(blst_limb_t));
static_assert(sizeof(blst_fp12) == 576 && alignof(blst_fp12) == alignof(blst_limb_t));
static_assert(sizeof(blst_p1) == 144 && alignof(blst_p1) == alignof(blst_limb_t));
static_assert(sizeof(blst_p1_affine) == 96 && alignof(blst_p1_affine) == alignof(blst_limb_t));
static_assert(sizeof(blst_p2) == 288 && alignof(blst_p2) == alignof(blst_limb_t));
static_assert(sizeof(blst_p2_affine) == 192 && alignof(blst_p2_affine) == alignof(blst_limb_t));

static_assert(offsetof(blst_fr, l) == 0);
static_assert(offsetof(blst_fp, l) == 0);
static_assert(offsetof(blst_fp2, fp) == 0);
static_assert(offsetof(blst_fp6, fp2) == 0);
static_assert(offsetof(blst_fp12, fp6) == 0);
static_assert(offsetof(blst_p1, x) == 0);
static_assert(offsetof(blst_p1_affine, x) == 0);
static_assert(offsetof(blst_p2, x) == 0);
static_assert(offsetof(blst_p2_affine, x) == 0);
#endif

#endif /* BLST_TYPES_H */
