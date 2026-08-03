/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BLST_TYPES_H
#define BLST_TYPES_H

#include <stddef.h>
#include <stdint.h>

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

/**
 * Stable error codes returned by the plain C API.
 *
 * Invariants:
 * - Numeric values are part of the ABI.
 * - The same codes are used by the scalar, AVX, and middle-product headers.
 */
typedef enum {
    BLST_SUCCESS = 0,
    BLST_BAD_ENCODING = 1,
    BLST_POINT_NOT_ON_CURVE = 2,
    BLST_POINT_NOT_IN_GROUP = 3,
    BLST_AGGR_TYPE_MISMATCH = 4,
    BLST_VERIFY_FAIL = 5,
    BLST_PK_IS_INFINITY = 6,
    BLST_BAD_SCALAR = 7,
} BLST_ERROR;

/** Library-owned byte type used by the public C API. */
typedef uint8_t blst_byte;
/** Library-owned 64-bit arithmetic limb type used by public scratch APIs. */
typedef uint64_t blst_limb_t;

#ifdef __cplusplus
# include <cstddef>
# include <type_traits>
namespace blst {
using byte = ::blst_byte;
using limb_t = ::blst_limb_t;
}
# if defined(__CUDACC__)
#  define BLST_TYPES_HD __host__ __device__
# else
#  define BLST_TYPES_HD
# endif
#endif

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
typedef struct { blst_byte b[256 / 8]; } blst_scalar;

/**
 * Scalar field element in Montgomery form.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - All arithmetic functions interpret values in Montgomery form.
 */
#ifdef __cplusplus
namespace blst {
struct UFr;
}

namespace blst::scalar {

/**
 * @brief Scalar field Fr = Z_r, where r is the BLS12-381 group order.
 *
 * Values are stored in Montgomery form unless explicitly documented otherwise.
 * If `a'` is the canonical integer in `[0, r)`, the stored value is
 * `a = a' * R mod r`, where `R = 2^256`. Unless stated otherwise, operations
 * are modulo r, return reduced results, preserve Montgomery form, and are
 * intended to be constant-time with respect to secret data.
 */
struct blst_fr
{
    using Fr = blst_fr;

    blst_limb_t l[256 / 8 / sizeof(blst_limb_t)];

    /// Access limb @p idx in little-endian limb order.
    BLST_TYPES_HD blst_limb_t& operator[](size_t idx) { return l[idx]; }
    BLST_TYPES_HD const blst_limb_t& operator[](size_t idx) const { return l[idx]; }

    /// @brief Additive identity in Montgomery form.
    BLST_TYPES_HD static const Fr& zero();
    /// @brief Multiplicative identity in Montgomery form (`R mod r`).
    BLST_TYPES_HD static const Fr& one();
    /// @brief Raw conventional modulus r stored in Fr-sized storage.
    /// @note This value is not a reduced field element and is not Montgomery encoded.
    BLST_TYPES_HD static constexpr const Fr& modulus();
    /// @brief Return `-r^{-1} mod 2^64`, the Montgomery reduction constant.
    BLST_TYPES_HD static constexpr blst_limb_t montgomery_constant();
    /// @brief `ret = a + b mod r`.
    BLST_TYPES_HD static void add(Fr& ret, const Fr& a, const Fr& b);
    /// @brief `ret = a - b mod r`.
    BLST_TYPES_HD static void sub(Fr& ret, const Fr& a, const Fr& b);
    /// @brief `ret = a * b mod r` with Montgomery inputs and output.
    BLST_TYPES_HD static void mul(Fr& ret, const Fr& a, const Fr& b);
    /// @brief `ret = a^x mod r`.
    /// @param x Little-endian exponent bytes.
    /// @param x_bits Number of meaningful bits in @p x; zero returns one.
    BLST_TYPES_HD static void exp(Fr& ret, const Fr& a, const blst_byte* x, size_t x_bits);
    /// @brief `ret = a^2 mod r`.
    BLST_TYPES_HD static void sqr(Fr& ret, const Fr& a);
    /// @brief Convert from Montgomery form to the canonical representative.
    BLST_TYPES_HD static void from_mont(Fr& ret, const Fr& a);
    /// @brief Convert a canonical representative to Montgomery form.
    BLST_TYPES_HD static void to_mont(Fr& ret, const Fr& a);
    /// @brief Montgomery-reduce the double-width value @p a modulo r.
    BLST_TYPES_HD static void redc(Fr& ret, const ::blst::UFr& a);
    /// @brief Constant-time conditional negation: `ret = flag ? -a : a`.
    BLST_TYPES_HD static void cneg(Fr& ret, const Fr& a, blst_limb_t flag);
    /// @brief `ret = -a mod r`.
    BLST_TYPES_HD static void neg(Fr& ret, const Fr& a);
    /// @brief `ret = a * 2^count mod r`.
    BLST_TYPES_HD static void lshift(Fr& ret, const Fr& a, size_t count);
    /// @brief `ret = a / 2^count mod r`.
    BLST_TYPES_HD static void rshift(Fr& ret, const Fr& a, size_t count);
    /// @brief `ret = 3 * a mod r`.
    BLST_TYPES_HD static void mul_by_3(Fr& ret, const Fr& a);
    /// @brief Compute `out = 1 / inp` with the selected constant-time inverse.
    BLST_TYPES_HD static void reciprocal(Fr& out, const Fr& inp);
    /// @brief Constant-time equality of the represented field elements.
    BLST_TYPES_HD static blst_limb_t equal(const Fr& a, const Fr& b);
    /// @brief Constant-time Cooley-Tukey butterfly.
    BLST_TYPES_HD static void ct_bfly(Fr& x0, Fr& x1, const Fr& twiddle);
    /// @brief Constant-time Gentleman-Sande butterfly.
    BLST_TYPES_HD static void gs_bfly(Fr& x0, Fr& x1, const Fr& twiddle);
    /// @brief Return one iff @p a is strictly less than r.
    BLST_TYPES_HD static blst_limb_t is_redc(const Fr& a);
    /// @brief Return one iff @p a is reduced modulo r and nonzero.
    BLST_TYPES_HD static blst_limb_t check_mod(const Fr& a);
    /// @brief Decode a 32-byte little-endian scalar and convert it to Montgomery form.
    BLST_TYPES_HD static void from_scalar(Fr& ret, const blst_byte a[32]);
    /// @brief Convert from Montgomery form to a 32-byte little-endian scalar.
    BLST_TYPES_HD static void to_scaler(blst_byte ret[32], const Fr& a);
    /// @brief Convert four little-endian 64-bit words to Montgomery form.
    BLST_TYPES_HD static void from_uint64(Fr& ret, const uint64_t a[4]);
    /// @brief Convert a single unsigned integer to Montgomery form.
    BLST_TYPES_HD static void from_uint64(Fr& ret, uint64_t a);
    /// @brief Return a Montgomery field element constructed from @p a.
    BLST_TYPES_HD static Fr from_uint64(uint64_t a);
    /// @brief Parse a complete hexadecimal field value without modifying @p ret on failure.
    BLST_NODISCARD BLST_TYPES_HD static bool from_hexascii(Fr& ret, const char* hex);
    /// @brief Write a null-terminated hexadecimal value; @p out must hold at least 65 bytes.
    BLST_TYPES_HD static void to_hexascii(char* out, const Fr& value);
    /// @brief Number of 64-bit limbs in the representation.
    BLST_TYPES_HD static constexpr size_t num_limbs() { return 4; }

    BLST_TYPES_HD bool operator==(const Fr& other) const { return equal(*this, other) != 0; }
    BLST_TYPES_HD bool operator!=(const Fr& other) const { return !(*this == other); }
};

using Fr = blst_fr;

} // namespace blst::scalar

typedef blst::scalar::blst_fr blst_fr;

#else
typedef struct { blst_limb_t l[256 / 8 / sizeof(blst_limb_t)]; } blst_fr;
#endif

/**
 * Base field element in Montgomery form.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - All arithmetic functions interpret values in Montgomery form.
 */
#ifdef __cplusplus
namespace blst {
struct UFp;
}

namespace blst::scalar {

/**
 * @brief BLS12-381 base field Fp.
 *
 * Values are stored in Montgomery form unless explicitly documented otherwise.
 * Here `R = 2^384 mod p`; `from_mont` returns the conventional representative
 * in `[0, p)`. Unless stated otherwise, operations are modulo p, return reduced
 * Montgomery values, and are intended to be constant-time with respect to
 * secret data. `reciprocal` uses the Lehmer almost-inverse and can be
 * variable-time; use `flt_reciprocal` where strict constant time is required.
 */
struct blst_fp
{
    using Fp = blst_fp;

    blst_limb_t l[384 / 8 / sizeof(blst_limb_t)];

    /// Access limb @p idx in little-endian limb order.
    BLST_TYPES_HD blst_limb_t& operator[](size_t idx) { return l[idx]; }
    BLST_TYPES_HD const blst_limb_t& operator[](size_t idx) const { return l[idx]; }

    /** @name Canonical constants */
    ///@{
    /// @brief Additive identity in Montgomery form.
    BLST_TYPES_HD static const Fp& zero();
    /// @brief Multiplicative identity in Montgomery form (`R mod p`).
    BLST_TYPES_HD static const Fp& one();
    /// @brief Raw conventional modulus p stored in Fp-sized storage.
    /// @note This value is not a reduced field element and is not Montgomery encoded.
    BLST_TYPES_HD static constexpr const Fp& modulus();
    /// @brief `-p mod 2^384`, stored in Fp-sized storage.
    BLST_TYPES_HD static constexpr const Fp& neg_modulus();
    ///@}

    /// @brief Return `-p^{-1} mod 2^64`, the Montgomery reduction constant.
    BLST_TYPES_HD static constexpr blst_limb_t montgomery_constant();

    /** @name Basic arithmetic */
    ///@{
    /// @brief `ret = a + b mod p`.
    BLST_TYPES_HD static void add(Fp& ret, const Fp& a, const Fp& b);
    /// @brief `ret = a - b mod p`.
    BLST_TYPES_HD static void sub(Fp& ret, const Fp& a, const Fp& b);
    /// @brief `ret = a * b mod p`, preserving Montgomery form.
    BLST_TYPES_HD static void mul(Fp& ret, const Fp& a, const Fp& b);
    /// @brief `ret = a^x mod p`.
    /// @param x Little-endian exponent bytes.
    /// @param x_bits Number of meaningful bits in @p x; zero returns one.
    BLST_TYPES_HD static void exp(Fp& ret, const Fp& a, const blst_byte* x, size_t x_bits);
    /// @brief `ret = a^2 mod p`.
    BLST_TYPES_HD static void sqr(Fp& ret, const Fp& a);
    ///@}

    /** @name Montgomery transforms */
    ///@{
    /// @brief Convert from Montgomery form to the canonical representative.
    BLST_TYPES_HD static void from_mont(Fp& ret, const Fp& a);
    /// @brief Convert a canonical representative to Montgomery form.
    /// @note If @p a is already Montgomery encoded, this multiplies it by R again.
    BLST_TYPES_HD static void to_mont(Fp& ret, const Fp& a);
    /// @brief Montgomery-reduce the double-width value @p a modulo p.
    BLST_TYPES_HD static void redc(Fp& ret, const ::blst::UFp& a);
    ///@}

    /** @name Negation and conditional operations */
    ///@{
    /// @brief Constant-time conditional negation: `ret = flag ? -a : a`.
    BLST_TYPES_HD static void cneg(Fp& ret, const Fp& a, blst_limb_t flag);
    /// @brief `ret = -a mod p`.
    BLST_TYPES_HD static void neg(Fp& ret, const Fp& a);
    ///@}

    /** @name Shifts and small multipliers */
    ///@{
    /// @brief `ret = a * 2^count mod p`.
    BLST_TYPES_HD static void lshift(Fp& ret, const Fp& a, size_t count);
    /// @brief `ret = a / 2^count mod p`.
    BLST_TYPES_HD static void rshift(Fp& ret, const Fp& a, size_t count);
    /// @brief `ret = a / 2 mod p`.
    BLST_TYPES_HD static void div_by_2(Fp& ret, const Fp& a);
    /// @brief `ret = 3 * a mod p`.
    BLST_TYPES_HD static void mul_by_3(Fp& ret, const Fp& a);
    /// @brief `ret = 8 * a mod p`.
    BLST_TYPES_HD static void mul_by_8(Fp& ret, const Fp& a);
    ///@}

    /** @name Inversion and square roots */
    ///@{
    /// @brief Compute `out = inp^{-1}` with the Lehmer almost-inverse.
    /// @note This operation is generally variable-time.
    BLST_TYPES_HD static void reciprocal(Fp& out, const Fp& inp);
    /// @brief Compute `out = inp^{-1}` using the fixed Fermat addition chain.
    BLST_TYPES_HD static void flt_reciprocal(Fp& out, const Fp& inp);
    /// @brief Compute `t0 = inp^((p-3)/4)` for square-root helpers.
    BLST_TYPES_HD static void recip_sqrt_3mod4(Fp& t0, const Fp& inp);
    /// @brief Compute `out = inp^{-1/2}` when it exists.
    /// @return One on success; zero for zero or a nonsquare.
    BLST_TYPES_HD static blst_limb_t recip_sqrt(Fp& out, const Fp& inp);
    /// @brief Write a square-root candidate for @p inp.
    /// @return One when the candidate is a root; zero for a nonsquare.
    /// @note A candidate is written on success and failure; no sign normalization is applied.
    BLST_TYPES_HD static blst_limb_t sqrt(Fp& out, const Fp& inp);
    /// @brief Return one iff @p inp is a quadratic residue, including zero.
    BLST_TYPES_HD static blst_limb_t is_sqrt(const Fp& inp);
    ///@}

    /** @name Exponentiation helpers */
    ///@{
    /// @brief `out = a^(2^count) * b`.
    /// @note Count zero performs one multiplication; exact `out == b` aliasing is supported.
    BLST_TYPES_HD static void sqr_n_mul(Fp& out, const Fp& a, size_t count, const Fp& b);
    ///@}

    /** @name Comparisons and predicates */
    ///@{
    /// @brief Constant-time equality of the represented field elements.
    BLST_TYPES_HD static blst_limb_t equal(const Fp& a, const Fp& b);
    /// @brief Return `(parity) | (sgn0 << 1)` for a conventional representative.
    BLST_TYPES_HD static blst_limb_t sgn0_pty(const Fp& a);
    /// @brief Return `(parity) | (sgn0 << 1)` directly from Montgomery form.
    BLST_TYPES_HD static blst_limb_t sgn0_pty_mont(const Fp& a);
    /// @brief Return one iff @p a is strictly less than p.
    BLST_TYPES_HD static blst_limb_t is_redc(const Fp& a);
    ///@}

    /// @brief Number of 64-bit limbs in the representation.
    BLST_TYPES_HD static constexpr size_t num_limb() { return 6; }
};

using Fp = blst_fp;

} // namespace blst::scalar

typedef blst::scalar::blst_fp blst_fp;

#else
typedef struct { blst_limb_t l[384 / 8 / sizeof(blst_limb_t)]; } blst_fp;
#endif

/**
 * Quadratic extension field element in Montgomery form.
 */
#ifdef __cplusplus
namespace blst::scalar {

/**
 * @brief Quadratic extension field Fp2 used by the BLS12-381 twist.
 *
 * Elements are `a0 + u*a1`, stored as `(a0, a1)`, with `u^2 = -1`.
 * Components are Montgomery-encoded Fp values. Unless stated otherwise,
 * operations accept and return Montgomery values and are intended to be
 * constant-time. `reciprocal` uses the generally variable-time Lehmer
 * almost-inverse; use `flt_reciprocal` where strict constant time is required.
 *
 * The extension-field sign follows RFC 9380: use `sgn0(a1)` when `a1 != 0`,
 * otherwise `sgn0(a0)`.
 */
struct blst_fp2
{
    using Fp = blst::scalar::Fp;
    using Fp2 = blst_fp2;

    Fp fp[2];

    /// Access component @p idx: zero is real and one is imaginary.
    BLST_TYPES_HD Fp& operator[](size_t idx) { return fp[idx]; }
    BLST_TYPES_HD const Fp& operator[](size_t idx) const { return fp[idx]; }

    /** @name Basic arithmetic */
    ///@{
    /// @brief `ret = a + b` in Fp2.
    BLST_TYPES_HD static void add(Fp2& ret, const Fp2& a, const Fp2& b);
    /// @brief `ret = a - b` in Fp2.
    BLST_TYPES_HD static void sub(Fp2& ret, const Fp2& a, const Fp2& b);
    /// @brief `ret = a * b` in Fp2 with `u^2 = -1`.
    BLST_TYPES_HD static void mul(Fp2& ret, const Fp2& a, const Fp2& b);
    /// @brief `ret = a^2` in Fp2.
    BLST_TYPES_HD static void sqr(Fp2& ret, const Fp2& a);
    /// @brief `out = inp^pow` using a little-endian exponent.
    /// @param pow_bits Number of meaningful bits in @p pow.
    BLST_TYPES_HD static void exp(Fp2& out, const Fp2& inp, const blst_byte* pow, size_t pow_bits);
    ///@}

    /** @name Negation and conditional operations */
    ///@{
    /// @brief Constant-time conditional negation: `ret = flag ? -a : a`.
    BLST_TYPES_HD static void cneg(Fp2& ret, const Fp2& a, blst_limb_t flag);
    /// @brief `ret = -a` in Fp2.
    BLST_TYPES_HD static void neg(Fp2& ret, const Fp2& a);
    ///@}

    /** @name Montgomery transforms */
    ///@{
    /// @brief Convert both components from Montgomery form.
    BLST_TYPES_HD static void from_mont(Fp2& ret, const Fp2& a);
    /// @brief Convert both components to Montgomery form.
    BLST_TYPES_HD static void to_mont(Fp2& ret, const Fp2& a);
    ///@}

    /** @name Inversion and square roots */
    ///@{
    /// @brief Compute `out = inp^{-1}` using the Lehmer almost-inverse.
    /// @note This operation is generally variable-time.
    BLST_TYPES_HD static void reciprocal(Fp2& out, const Fp2& inp);
    /// @brief Compute `out = inp^{-1}` using a fixed Fermat addition chain.
    BLST_TYPES_HD static void flt_reciprocal(Fp2& out, const Fp2& inp);
    /// @brief Compute the canonically aligned square root when it exists.
    /// @return One on success and zero for a nonsquare.
    BLST_TYPES_HD static blst_limb_t sqrt(Fp2& out, const Fp2& inp);
    ///@}

    /** @name Small multipliers and shifts */
    ///@{
    /// @brief `ret = 3 * a` in Fp2.
    BLST_TYPES_HD static void mul_by_3(Fp2& ret, const Fp2& a);
    /// @brief `ret = 8 * a` in Fp2.
    BLST_TYPES_HD static void mul_by_8(Fp2& ret, const Fp2& a);
    /// @brief `ret = a * 2^count` in Fp2.
    BLST_TYPES_HD static void lshift(Fp2& ret, const Fp2& a, size_t count);
    ///@}

    /** @name Conjugation and special multipliers */
    ///@{
    /// @brief `ret = conjugate(a) = (a0, -a1)`.
    BLST_TYPES_HD static void conjugate(Fp2& ret, const Fp2& a);
    /// @brief `ret = (u + 1) * a`, using `u^2 = -1`.
    BLST_TYPES_HD static void mul_by_u_plus_1(Fp2& ret, const Fp2& a);
    ///@}

    /** @name Quadratic-residue helpers */
    ///@{
    /// @brief Return one iff @p inp is a quadratic residue, including zero.
    BLST_TYPES_HD static int is_square(const Fp2& inp);
    /// @brief Compute `out = inp^{-1/2}` when it exists.
    /// @return One on success and zero otherwise.
    BLST_TYPES_HD static blst_limb_t recip_sqrt(Fp2& out, const Fp2& inp);
    /**
     * @brief Align a square-root candidate to the RFC 9380 sign convention.
     *
     * Selects between @p sqrt and its negation relative to @p inp. The @p ret
     * parameter carries an implementation-specific candidate/intermediate.
     *
     * @return One on success and zero when alignment is not applicable.
     */
    BLST_TYPES_HD static blst_limb_t sqrt_align(Fp2& out, const Fp2& ret,
        const Fp2& sqrt, const Fp2& inp);
    ///@}

    /** @name Sign and parity */
    ///@{
    /// @brief Return `(parity) | (sgn0 << 1)` directly from Montgomery form.
    BLST_TYPES_HD static blst_limb_t sgn0_pty_mont(const Fp2& a);
    /// @brief Return `(parity) | (sgn0 << 1)` for conventional components.
    BLST_TYPES_HD static blst_limb_t sgn0_pty(const Fp2& a);
    ///@}

    /** @name Constants */
    ///@{
    /// @brief Additive identity `(0, 0)`.
    BLST_TYPES_HD static const Fp2& zero();
    /// @brief Multiplicative identity `(1, 0)` in Montgomery form.
    BLST_TYPES_HD static const Fp2& one();
    ///@}
};

using Fp2 = blst_fp2;

} // namespace blst::scalar

typedef blst::scalar::blst_fp2 blst_fp2;

#else
typedef struct { blst_fp fp[2]; } blst_fp2;
#endif

/**
 * Sextic extension field element in Montgomery form.
 */
#ifdef __cplusplus
namespace blst::scalar {

/**
 * @brief Sextic extension field Fp6 represented as `c0 + c1*v + c2*v^2`.
 *
 * Each coefficient is an Fp2 value in Montgomery form. This type primarily
 * provides the ABI-stable storage used by Fp12 and pairing operations.
 */
struct blst_fp6
{
    using Fp2 = blst::scalar::Fp2;

    Fp2 fp2[3];

    /// Access coefficient @p idx in increasing powers of v.
    BLST_TYPES_HD Fp2& operator[](size_t idx) { return fp2[idx]; }
    BLST_TYPES_HD const Fp2& operator[](size_t idx) const { return fp2[idx]; }
};

using Fp6 = blst_fp6;

} // namespace blst::scalar

typedef blst::scalar::blst_fp6 blst_fp6;

#else
typedef struct { blst_fp2 fp2[3]; } blst_fp6;
#endif

/**
 * Twelfth-degree extension field element in Montgomery form.
 */
#ifdef __cplusplus
namespace blst::scalar {

/// @brief Precomputed constant-time GLV-4 exponentiation schedule.
struct GLV4_Exp
{
    uint8_t sign_mask;     ///< Basis elements that require conjugation.
    uint8_t sched_len;     ///< Number of valid schedule entries.
    uint8_t schedule[70];  ///< Compact square/multiply bytecode.
};

/**
 * @brief Precomputed constant-time GLV-2 joint-sparse-form schedule.
 *
 * Each index stores two signed digits; two digit pairs are packed per byte.
 */
struct GLV_JSF_Exp
{
    uint8_t sched_len;     ///< Number of valid digit pairs.
    uint8_t schedule[71];  ///< Packed signed-digit schedule.

    /// Set the @p idx digit pair to (`d1`, `d2`).
    BLST_TYPES_HD void set(size_t idx, int8_t d1, int8_t d2)
    {
        auto enc = ((d1 + 1) | ((d2 + 1) << 2));
        enc <<= ((idx % 2) * 4);
        enc |= schedule[idx / 2] & (~(0xF << ((idx % 2) * 4)));
        schedule[idx / 2] = static_cast<uint8_t>(enc);
    }

    struct Pair { int8_t d1; int8_t d2; };

    /// Return the @p idx signed-digit pair.
    BLST_TYPES_HD Pair get(size_t idx) const
    {
        auto enc = (schedule[idx / 2] >> ((idx % 2) * 4)) & 0xF;
        return { static_cast<int8_t>((enc & 0x3) - 1),
                 static_cast<int8_t>(((enc >> 2) & 0x3) - 1) };
    }
};

/**
 * @brief Twelfth-degree extension field Fp12 used by pairings and G_T.
 *
 * Elements are represented as two Fp6 coefficients in Montgomery form.
 * General arithmetic is valid for all Fp12 values; methods explicitly named
 * `cyclotomic_*` require their input to lie in the cyclotomic subgroup.
 */
struct blst_fp12
{
    using Fp6 = blst::scalar::Fp6;
    using Fp12 = blst_fp12;
    using GLV4_Exp = blst::scalar::GLV4_Exp;
    using GLV_JSF_Exp = blst::scalar::GLV_JSF_Exp;

    Fp6 fp6[2];

    /// Access Fp6 coefficient @p idx.
    BLST_TYPES_HD Fp6& operator[](size_t idx) { return fp6[idx]; }
    BLST_TYPES_HD const Fp6& operator[](size_t idx) const { return fp6[idx]; }

    /// @brief `ret = a + b` in Fp12.
    BLST_TYPES_HD static void add(Fp12& ret, const Fp12& a, const Fp12& b);
    /// @brief `ret = a - b` in Fp12.
    BLST_TYPES_HD static void sub(Fp12& ret, const Fp12& a, const Fp12& b);
    /// @brief `ret = a^2` using the scalar implementation.
    BLST_TYPES_HD static void sqr(Fp12& ret, const Fp12& a);
    /// @brief `ret = a^2` using scalar code generation explicitly.
    BLST_TYPES_HD static void sqr_scalar(Fp12& ret, const Fp12& a);
    /// @brief `ret = a * b` using the scalar implementation.
    BLST_TYPES_HD static void mul(Fp12& ret, const Fp12& a, const Fp12& b);
    /// @brief `ret = a * b` using scalar code generation explicitly.
    BLST_TYPES_HD static void mul_scalar(Fp12& ret, const Fp12& a, const Fp12& b);
    /// @brief Multiply by @p b or its conjugate.
    /// @param conj_b When true, compute `ret = a * conjugate(b)`.
    BLST_TYPES_HD static void mul(Fp12& ret, const Fp12& a, const Fp12& b, bool conj_b);
    /// @brief Standard exponentiation `ret = a^exp`.
    /// @param exp Little-endian exponent bytes.
    /// @param bits Number of meaningful exponent bits.
    BLST_TYPES_HD static void exp(Fp12& ret, const Fp12& a, const blst_byte* exp, size_t bits);
    /// @brief Multiplicative identity of Fp12.
    BLST_TYPES_HD static const Fp12& one();
    /// @brief Cyclotomic squaring using the scalar implementation.
    /// @pre @p a lies in the cyclotomic subgroup.
    BLST_TYPES_HD static void cyclotomic_sqr(Fp12& ret, const Fp12& a);
    /// @brief Cyclotomic squaring using scalar code generation explicitly.
    /// @pre @p a lies in the cyclotomic subgroup.
    BLST_TYPES_HD static void cyclotomic_sqr_scalar(Fp12& ret, const Fp12& a);
    /// @brief Exponentiate with cyclotomic squaring.
    /// @pre @p a lies in the cyclotomic subgroup.
    BLST_TYPES_HD static void cyclotomic_exp(Fp12& ret, const Fp12& a, const blst_byte* exp, size_t bits);
    /// @brief Basic cyclotomic exponentiation implementation.
    /// @pre @p a lies in the cyclotomic subgroup.
    BLST_TYPES_HD static void cyclotomic_exp_basic(Fp12& ret, const Fp12& a, const blst_byte* exp, size_t bits);
    /// @brief Split a 256-bit scalar into two approximately 128-bit GLV components.
    /// @details The decomposition satisfies `k = k1 + k2*lambda` modulo r.
    BLST_TYPES_HD static void glv_split_2(blst_limb_t k1_out[4], blst_limb_t k2_out[4], const blst_limb_t k[4]);
    /// @brief Recode two GLV components into joint sparse form.
    BLST_TYPES_HD static GLV_JSF_Exp glv_recode_jsf(const blst_limb_t* k1, const blst_limb_t* k2);
    /// @brief GLV-2 cyclotomic exponentiation from a 32-byte scalar.
    /// @pre @p base lies in the cyclotomic subgroup.
    BLST_TYPES_HD static void cyclotomic_exp_glv(Fp12& out, const Fp12& base, const blst_byte* scalar_bytes);
    /// @brief GLV-2 cyclotomic exponentiation with pre-split scalars.
    /// @pre @p base lies in the cyclotomic subgroup.
    BLST_TYPES_HD static void cyclotomic_exp_glv(Fp12& out, const Fp12& base, const blst_limb_t* k1, const blst_limb_t* k2);
    /// @brief GLV-2 cyclotomic exponentiation with on-the-fly JSF recoding.
    /// @pre @p base lies in the cyclotomic subgroup.
    BLST_TYPES_HD static void cyclotomic_exp_glv_jsf(Fp12& out, const Fp12& base, const blst_limb_t* k1, const blst_limb_t* k2);
    /// @brief GLV-2 cyclotomic exponentiation with a precomputed JSF schedule.
    /// @pre @p base lies in the cyclotomic subgroup.
    BLST_TYPES_HD static void cyclotomic_exp_glv_jsf(Fp12& out, const Fp12& base, const GLV_JSF_Exp& sched);
    /// @brief Four-dimensional G_T scalar decomposition.
    /// @details `k = k0 + k1*lambda1 + k2*lambda2 + k3*lambda3` modulo r.
    BLST_TYPES_HD static void gt_split_scalar_4d(int64_t k_out[4], const blst_limb_t k_in[4]);
    /// @brief Split a scalar into four signed GLV components.
    BLST_TYPES_HD static void glv_split_4(int64_t k_out[4], const blst_limb_t k_in[4]);
    /// @brief Compile a scalar into a constant-time GLV-4 execution schedule.
    BLST_TYPES_HD static GLV4_Exp glv4_compile(const blst_limb_t* k_in);
    /// @brief GLV-4 cyclotomic exponentiation with a precompiled schedule.
    /// @pre @p base lies in the cyclotomic subgroup.
    BLST_TYPES_HD static void cyclotomic_exp_glv4(Fp12& out, const Fp12& base, const GLV4_Exp& props);
    /// @brief GLV-4 cyclotomic exponentiation with a runtime scalar split.
    /// @pre @p base lies in the cyclotomic subgroup.
    BLST_TYPES_HD static void cyclotomic_exp_glv4(Fp12& out, const Fp12& base, const blst_limb_t* scalar);
    /// @brief Conjugate @p a in place by negating its w coefficient.
    BLST_TYPES_HD static void conjugate(Fp12& a);
    /// @brief `ret = a^{-1}` in Fp12.
    BLST_TYPES_HD static void inverse(Fp12& ret, const Fp12& a);
    /// @brief Apply the Frobenius map: `ret = a^(p^n)`.
    BLST_TYPES_HD static void frobenius_map(Fp12& ret, const Fp12& a, size_t n);
    /// @brief Multiply by a sparse `xy00z0` element using the scalar implementation.
    BLST_TYPES_HD static void mul_by_xy00z0(Fp12& ret, const Fp12& a, const Fp6& xy00z0);
    /// @brief Multiply by a sparse `xy00z0` element using scalar code generation.
    BLST_TYPES_HD static void mul_by_xy00z0_scalar(Fp12& ret, const Fp12& a, const Fp6& xy00z0);
    /// @brief Constant-time equality predicate.
    BLST_TYPES_HD static int equal(const Fp12& a, const Fp12& b);
    /// @brief Constant-time predicate for the multiplicative identity.
    BLST_TYPES_HD static int is_one(const Fp12& a);
    /// @brief Decode 576 big-endian bytes and convert coefficients to Montgomery form.
    BLST_TYPES_HD static void from_bendian(Fp12& ret, const blst_byte serialized[12 * 48]);
    /// @brief Convert from Montgomery form and encode 576 big-endian bytes.
    BLST_TYPES_HD static void to_bendian(blst_byte serialized[12 * 48], const Fp12& a);
    /// @brief Return one iff @p f is in the cyclotomic subgroup G_phi_12.
    BLST_TYPES_HD static int in_group(const Fp12& f);
    /// @brief Faster cyclotomic-subgroup predicate using Frobenius relations.
    BLST_TYPES_HD static blst_limb_t is_cyclotomic(const Fp12& f);
    /// @brief Compute `ret = (ret * a)^(2^n)`.
    BLST_TYPES_HD static void mul_n_sqr(Fp12& ret, const Fp12& a, size_t n);
    /// @brief Final-exponentiation helper: `ret = a^((z-1)/2)`.
    BLST_TYPES_HD static void raise_to_z_div_by_2(Fp12& ret, const Fp12& a);
    /// @brief Final-exponentiation helper: `a = b^z`.
    BLST_TYPES_HD static void raise_to_z(Fp12& a, const Fp12& b);
    /// @brief Perform pairing final exponentiation into G_T.
    BLST_TYPES_HD static void final_exp(Fp12& ret, const Fp12& f);
};

using Fp12 = blst_fp12;

} // namespace blst::scalar

typedef blst::scalar::blst_fp12 blst_fp12;

#else
typedef struct { blst_fp6 fp6[2]; } blst_fp12;
#endif

/**
 * G1 projective point over Montgomery Fp coordinates.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - The point may be the point at infinity.
 */
#ifdef __cplusplus
namespace blst::scalar {

template<typename Field> struct AffinePoint;
template<typename Field> struct Pointxyzz;

/**
 * @tparam Field Either Fp for G1 or Fp2 for G2.
 * @brief Jacobian point `(X:Y:Z)` over @p Field.
 *
 * Infinity is represented by `Z == 0`. Affine inputs use the all-zero encoding
 * for infinity. Field coordinates are stored in Montgomery form. Unless a
 * method is explicitly documented otherwise, point operations are intended to
 * be constant-time with respect to secret data.
 */
template<typename Field>
struct Point
{
    using Field_type = Field;

    constexpr Point() = default;
    constexpr Point(const Point&) = default;
    constexpr Point& operator=(const Point&) = default;
    /// Construct from Jacobian coordinates already in Montgomery form.
    constexpr Point(const Field& x_in, const Field& y_in, const Field& z_in)
        : x(x_in), y(y_in), z(z_in) {}

    Field x, y, z;

    /**
     * @brief Unified constant-time add-or-double of two Jacobian points.
     * @param a4 Optional curve-a term; BLS12-381 uses `a = 0`, so pass null.
     * @note Handles equality, opposites, and infinity.
     */
    BLST_TYPES_HD static void dadd(Point& out, const Point& p1, const Point& p2,
        const Field* a4 = nullptr);
    /// @brief Unified constant-time add-or-double with an affine addend.
    /// @note Treats the all-zero affine value as infinity.
    BLST_TYPES_HD static void dadd(Point& out, const Point& p1, const AffinePoint<Field>& p2);
    /// @brief Add two distinct, nonopposite Jacobian points, handling infinity.
    /// @note Uses the EFD add-2007-bl formula.
    BLST_TYPES_HD static void add(Point& out, const Point& p1, const Point& p2);
    /// @brief Mixed Jacobian-plus-affine addition, handling infinity.
    /// @note Uses the EFD madd-2007-bl formula and requires distinct, nonopposite inputs.
    BLST_TYPES_HD static void add(Point& out, const Point& p1, const AffinePoint<Field>& p2);
    /// @brief Point doubling using the EFD dbl-2009-l formula.
    BLST_TYPES_HD static void doubling(Point& out, const Point& p1);
    /// @brief Repeated doubling: `out = 2^n * p`.
    BLST_TYPES_HD static void add_n_dbl(Point& out, const Point& p, size_t n);
    /**
     * @brief Variable-time scalar multiplication with algorithm selection.
     * @param scalar Little-endian scalar bytes.
     * @param nbits Number of meaningful scalar bits.
     *
     * Short scalars use a width-four window. Valid scalars up to 256 bits use
     * GLS/GLV acceleration; other inputs use the width-five implementation.
     */
    BLST_TYPES_HD static void mult(Point& out, const Point& a, const blst_byte* scalar, size_t nbits);
    /// @brief Width-five windowed multiplication for internal or large-scalar use.
    BLST_TYPES_HD static void mult_w5(Point& out, const Point& point, const blst_byte* scalar, size_t nbits);
    /// @brief Constant-time conditional y-coordinate negation when @p cbit is one.
    BLST_TYPES_HD static void cneg(Point& p, blst_limb_t cbit);
    /// @brief Negate the y coordinate.
    BLST_TYPES_HD static void neg(Point& out, const Point& p);
    /// @brief Convert Jacobian to affine; infinity maps to the all-zero encoding.
    BLST_TYPES_HD static void to_affine(AffinePoint<Field>& out, const Point& in);
    /// @brief Convert affine to Jacobian; finite points receive `Z = 1`.
    BLST_TYPES_HD static void from_affine(Point& out, const AffinePoint<Field>& a);
    /// @brief Normalize a Jacobian point while preserving infinity.
    BLST_TYPES_HD static void from_jacobian(Point& out, const Point& in);
    /// @brief Constant-time conditional swap when @p cbit is one.
    BLST_TYPES_HD static void cswap(Point& a, Point& b, blst_limb_t cbit);
    /// @brief Constant-time conditional copy `a = b` when @p cbit is one.
    BLST_TYPES_HD static void ccopy(Point& a, const Point& b, blst_limb_t cbit);
    /// @brief Constant-time equality of affine representatives.
    BLST_TYPES_HD static blst_limb_t equal(const Point& p1, const Point& p2);
    /// @brief Convert Jacobian to the variable-time XYZZ representation.
    BLST_TYPES_HD static void to_xyzz(Pointxyzz<Field>& out, const Point& in);
    /// @brief Group generator in Jacobian coordinates.
    BLST_TYPES_HD static const Point& generator();
    /// @brief Return the point at infinity (`Z == 0`).
    BLST_TYPES_HD static Point infinity();
    /// @brief Return the group identity, equivalent to infinity.
    BLST_TYPES_HD static Point identity() { return infinity(); }
    /// @brief Return one iff @p p is the group identity.
    BLST_TYPES_HD static blst_limb_t is_identity(const Point& p) { return is_inf(p); }
    /// @brief Return curve coefficient B for the selected group.
    BLST_TYPES_HD static const Field& get_B();
    /// @brief Multiply @p in by curve coefficient B.
    BLST_TYPES_HD static void mul_by_b(Field& out, const Field& in);
    /// @brief Multiply @p in by `4B`, used by doubling formulas.
    BLST_TYPES_HD static void mul_by_4b(Field& out, const Field& in);
    /// @brief Return one iff @p p is infinity.
    BLST_TYPES_HD static int is_inf(const Point& p);
    /// @brief Return one iff @p p is on the curve; infinity is accepted.
    BLST_TYPES_HD static blst_limb_t on_curve(const Point& p);
    /// @brief Return one iff @p p lies in the prime-order subgroup.
    BLST_TYPES_HD static blst_limb_t in_subgroup(const Point& p);
    /** @name Serialization */
    ///@{
    /// @brief Write uncompressed big-endian affine coordinates and return sgn0 parity.
    BLST_TYPES_HD static blst_limb_t serialize_BE(blst_byte* out, const Point& in);
    /// @brief Write uncompressed IETF form, using `0x40 || zeros` for infinity.
    BLST_TYPES_HD static void serialize(blst_byte* out, const Point& in);
    /// @brief Write compressed big-endian x without header flags and return y parity.
    BLST_TYPES_HD static blst_limb_t compress_BE(blst_byte* out, const Point& in);
    /// @brief Write IETF compressed form with compression, sign, and infinity flags.
    BLST_TYPES_HD static void compress(blst_byte* out, const Point& in);
    ///@}
    /**
     * @brief Endomorphism-accelerated multiplication for a valid Fr scalar.
     *
     * G1 uses the beta endomorphism; G2 uses the psi endomorphism. This is the
     * GLS/GLV path used by signing and public-key derivation.
     */
    BLST_TYPES_HD static void mult_glsv(Point& out, const Point& in, const blst_byte SK[32]);
    /// @brief Compute a BLS signature point `out = SK * in` and normalize it.
    BLST_TYPES_HD static void sign(Point& out, const Point& in, const blst_byte SK[32]);
    /**
     * @brief Signature convenience API.
     * @param out Optional compressed output bytes.
     * @param sig Optional affine output.
     * @param hash Hash-to-curve point to sign.
     * @param SK Little-endian 32-byte secret scalar.
     */
    BLST_TYPES_HD static void sign2(blst_byte* out, AffinePoint<Field>* sig,
        const Point& hash, const blst_byte SK[32]);
    /// @brief Derive public key `out = SK * generator()`.
    BLST_TYPES_HD static void to_pk(Point& out, const blst_byte SK[32]);
    /// @brief Optionally emit compressed and affine public-key outputs.
    BLST_TYPES_HD static void to_pk2(blst_byte* out, AffinePoint<Field>* pk, const blst_byte SK[32]);
};

/**
 * @tparam Field Either Fp for G1 or Fp2 for G2.
 * @brief Affine point `(X,Y)` over @p Field.
 *
 * Infinity is encoded as all-zero x and y coordinates. Field coordinates are
 * in Montgomery form. G1 encodings are 96 bytes uncompressed and 48 bytes
 * compressed; G2 encodings are 192 and 96 bytes respectively.
 */
template<typename Field>
struct AffinePoint
{
    Field x, y;

    /** @name Uncompressed serialization */
    ///@{
    /// @brief Write `X || Y` in big-endian form and return sgn0 parity.
    BLST_TYPES_HD static blst_limb_t serialize_BE(blst_byte* out, const AffinePoint& in);
    /// @brief Write uncompressed IETF form, including the infinity encoding.
    BLST_TYPES_HD static void serialize(blst_byte* out, const AffinePoint& in);
    /// @brief Parse raw big-endian `X || Y`, validating range and curve membership.
    BLST_TYPES_HD static BLST_ERROR deserialize_BE(AffinePoint& out, const blst_byte* in);
    /// @brief Parse compressed, uncompressed, or infinity-prefixed IETF input.
    BLST_TYPES_HD static BLST_ERROR deserialize(AffinePoint& out, const blst_byte* in);
    ///@}
    /** @name Compressed serialization */
    ///@{
    /// @brief Write big-endian x without header flags and return y parity.
    BLST_TYPES_HD static blst_limb_t compress_BE(blst_byte* out, const AffinePoint& in);
    /// @brief Write IETF compressed form with compression, sign, and infinity flags.
    BLST_TYPES_HD static void compress(blst_byte* out, const AffinePoint& in);
    /// @brief Parse raw big-endian x, reconstruct y, and return sgn0 parity.
    BLST_TYPES_HD static blst_limb_t uncompress_BE(AffinePoint& out, const blst_byte* in);
    /// @brief Parse and validate IETF compressed form.
    BLST_TYPES_HD static BLST_ERROR uncompress(AffinePoint& out, const blst_byte* in);
    ///@}
    /// @brief Group generator in affine coordinates.
    BLST_TYPES_HD static const AffinePoint& generator();
    /// @brief Constant-time equality on raw limbs.
    BLST_TYPES_HD static int equal(const AffinePoint& p1, const AffinePoint& p2);
    /// @brief Return one iff @p p satisfies the curve equation or is infinity.
    BLST_TYPES_HD static blst_limb_t on_curve(const AffinePoint& p);
    /// @brief Return one iff both coordinates are all zero.
    BLST_TYPES_HD static blst_limb_t is_inf(const AffinePoint& p);
    /// @brief Return one iff @p p lies in the prime-order subgroup.
    BLST_TYPES_HD static blst_limb_t in_subgroup(const AffinePoint& p);
    /// @brief Constant-time conditional y-coordinate negation.
    BLST_TYPES_HD static void cneg(AffinePoint& point, blst_limb_t c);
};

} // namespace blst::scalar

typedef blst::scalar::Point<blst::scalar::Fp> blst_p1;
#else
typedef struct { blst_fp x, y, z; } blst_p1;
#endif

/**
 * G1 affine point over Montgomery Fp coordinates.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - The point may be the point at infinity.
 */
#ifdef __cplusplus
typedef blst::scalar::AffinePoint<blst::scalar::Fp> blst_p1_affine;
#else
typedef struct { blst_fp x, y; } blst_p1_affine;
#endif

/**
 * G2 projective point over Montgomery Fp2 coordinates.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - The point may be the point at infinity.
 */
#ifdef __cplusplus
typedef blst::scalar::Point<blst::scalar::Fp2> blst_p2;
#else
typedef struct { blst_fp2 x, y, z; } blst_p2;
#endif

/**
 * G2 affine point over Montgomery Fp2 coordinates.
 *
 * Invariants:
 * - The layout is ABI-stable.
 * - The point may be the point at infinity.
 */
#ifdef __cplusplus
typedef blst::scalar::AffinePoint<blst::scalar::Fp2> blst_p2_affine;
#else
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

# undef BLST_TYPES_HD
#endif

#endif /* BLST_TYPES_H */
