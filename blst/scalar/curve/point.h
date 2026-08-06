/*
 * Copyright Supranational LLC
 * Copyright Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_SCALAR_CURVE_POINT_H
#define BLST_SCALAR_CURVE_POINT_H


 /**
  * @file
  * @brief BLS12-381 elliptic-curve point types and operations over Fp and Fp2.
  *
  * # Curves and fields
  * - **G1** is a short-Weierstrass curve over the base field @c Fp with model
  *   ` y^2 = x^3 + B ` (parameter B returned by get_BLST_B_E1()).
  * - **G2** is the twist over @c Fp2 with model
  *   ` y^2 = x^3 + B ` (parameter B returned by get_BLST_B_E2()).
  *
  * # Coordinates & representations
  * - @c Point<Field> stores Jacobian coordinates (X:Y:Z) with the point at
  *   infinity encoded as @c Z == 0. @c Field is either @c Fp (G1) or @c Fp2 (G2).
  * - @c AffinePoint<Field> stores affine (X,Y). Infinity is encoded as the
  *   all-zero struct (both X and Y are zero limbs).
  * - @c Pointxyzz<Field> is an XYZZ representation used for fast, **non-constant-time**
  *   mixed additions in multi-scalar routines (see class docs).
  *
  * # Field element format
  * All field elements in memory are in **Montgomery form**. Conversion in/out of
  * Montgomery form is handled internally where needed (e.g., serialization).
  *
  * # (Un)compressed encodings (big-endian)
  * - @c serialize_BE / @c deserialize_BE : uncompressed big-endian coordinates.
  *   - G1 (Fp): 96 bytes = X(48) || Y(48)
  *   - G2 (Fp2): 192 bytes = X.im(48) || X.re(48) || Y.im(48) || Y.re(48)
  *   - No leading flags are set; callers may inspect contents. Infinity is
  *     represented by @c serialize() (not @c serialize_BE) using a flag.
  * - @c serialize / @c deserialize : uncompressed with **infinity bit** support.
  *   - If infinity, first byte has 0x40 and the rest are zeros.
  * - @c compress_BE : compressed big-endian X with sgn0 parity in return value.
  * - @c compress / @c uncompress : IETF-style compressed form
  *   - First byte sets 0x80 (compressed). If infinity, it sets 0xC0 and zeroes.
  *   - The y-sign/parity bit is encoded using sgn0; see code for bit layout.
  *
  * # Constant-time notes
  * Arithmetic explicitly documented as constant-time avoids secret-dependent
  * control flow and memory access under its stated preconditions. Conversion,
  * serialization, decoding, membership predicates, and @c Pointxyzz helpers may
  * be variable-time and must operate only on public points or encodings.
  *
  * # CUDA / host constants
  * get_BLST_* accessors return host or device constants depending on compilation target.
  * This allows the same header to be used on CPU or CUDA kernels.
  */


#include "blst/config.h"
#include "blst/common/vect.h"
#include "blst/common/bytes.h"
#include "blst/common/errors.h"
#include <stdio.h>
#include <type_traits>
#include "blst/fields/Fr.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/scalar/curve/ec_mult.h"

namespace blst::scalar
{
	template<typename Field>
	struct AffinePoint;

	template<typename Field>
	struct Pointxyzz;

	/**
	 * @tparam Field Either Fp for G1 or Fp2 for G2.
	 * @brief Jacobian point `(X:Y:Z)` over @p Field.
	 *
	 * Infinity is represented by `Z == 0`. Affine inputs use the all-zero encoding
	 * for infinity. Field coordinates are stored in Montgomery form. Timing
	 * properties are stated per operation; no blanket constant-time guarantee
	 * applies to conversion, encoding, or validation methods. The type stores a
	 * point representation; construction or decoding does not establish subgroup
	 * membership. Unless a method states a subgroup precondition, basic
	 * single-point arithmetic, conversion, observation, and encoding require only
	 * an on-curve input and preserve (but do not establish) subgroup membership.
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
		 * @pre Both inputs are on-curve; subgroup membership is not required.
		 * @param a4 Optional curve-a term; BLS12-381 uses `a = 0`, so pass null.
		 * @note Handles equality, opposites, and infinity.
		 */
		BLST_HD static void dadd(Point& out, const Point& p1, const Point& p2,
			const Field* a4 = nullptr);
		/// @brief Unified constant-time add-or-double with an affine addend.
		/// @pre Both inputs are on-curve; subgroup membership is not required.
		/// @note Treats the all-zero affine value as infinity.
		BLST_HD static void dadd(Point& out, const Point& p1, const AffinePoint<Field>& p2);
		/// @brief Add two Jacobian points using the fast incomplete formula.
		/// @pre Both inputs are on-curve; subgroup membership is not required.
		/// @pre If both inputs are finite, they are distinct and `p2 != -p1`.
		/// @note Infinity is handled. Equal points require @c dadd for doubling;
		/// opposite points (same affine x, negated affine y) require @c dadd to
		/// produce infinity. Uses the EFD add-2007-bl formula.
		BLST_HD static void add(Point& out, const Point& p1, const Point& p2);
		/// @brief Mixed Jacobian-plus-affine addition using the fast incomplete formula.
		/// @pre Both inputs are on-curve; subgroup membership is not required.
		/// @pre If both inputs are finite, they are distinct and `p2 != -p1`.
		/// @note Infinity is handled. Use @c dadd when the finite inputs may be
		/// equal or opposite. Uses the EFD madd-2007-bl formula.
		BLST_HD static void add(Point& out, const Point& p1, const AffinePoint<Field>& p2);
		/// @brief Double an on-curve point using the EFD dbl-2009-l formula.
		/// @note Subgroup membership is not required or validated.
		BLST_HD static void doubling(Point& out, const Point& p1);
		/// @brief Accumulate @p p into @p out, then double the result @p n times:
		/// `out = 2^n * (out + p)`.
		BLST_HD static void add_n_dbl(Point& out, const Point& p, size_t n);
		/**
		 * @brief Constant-time scalar multiplication.
		 * @param scalar Little-endian scalar bytes.
		 * @param nbits Public number of meaningful scalar bits, at most 256.
		 *
		 * For canonical scalars in `[0,r)`, runtime and memory access are independent
		 * of scalar contents for a fixed @p nbits. Noncanonical values are accepted
		 * modulo r, but timing may reveal that reduction was required. The input
		 * must already be an on-curve prime-order subgroup point; Release builds do
		 * not validate that precondition. @p scalar may be null only for zero bits.
		 * Exact @p out / @p a aliasing is supported; partial overlap is not.
		 */
		BLST_HD static void mult(Point& out, const Point& a, const byte* scalar, size_t nbits);
		/**
		 * @brief Optimized variable-time scalar multiplication.
		 * @warning Scalar range and window values may affect control flow and memory
		 * access. Use only with public scalar data. The input must already be an
		 * on-curve prime-order subgroup point; Release builds do not validate it.
		 * @p nbits must not exceed 256, and @p scalar may be null only for zero bits.
		 * Exact @p out / @p a aliasing is supported; partial overlap is not.
		 */
		BLST_HD static void mult_vartime(Point& out, const Point& a, const byte* scalar, size_t nbits);
		/**
		 * @brief Constant-time generic multiplication without subgroup endomorphisms.
		 * @warning The input must be on curve, but need not be in the prime-order
		 * subgroup. Release builds do not validate that precondition. @p nbits must
		 * not exceed 256, and @p scalar may be null only for zero bits.
		 * The exact unsigned scalar is applied. Exact @p out / @p a aliasing is
		 * supported; partial overlap is not.
		 */
		BLST_HD static void unchecked_mult(Point& out, const Point& a, const byte* scalar, size_t nbits);
		/// @brief Constant-time conditional y-coordinate negation when @p cbit is one.
		/// @pre @p p is on-curve; subgroup membership is not required.
		BLST_HD static void cneg(Point& p, limb_t cbit);
		/// @brief Negate the y coordinate of an on-curve point.
		/// @note Subgroup membership is not required or validated.
		BLST_HD static void neg(Point& out, const Point& p);
		/// @brief Convert an on-curve Jacobian point to affine; infinity maps to zero.
		/// @note Subgroup membership is preserved but neither required nor established.
		/// @note Variable-time with respect to the point representation.
		BLST_HD static void to_affine(AffinePoint<Field>& out, const Point& in);
		/// @brief Convert an on-curve affine point to Jacobian; finite points get `Z = 1`.
		/// @note Subgroup membership is preserved but neither required nor established.
		BLST_HD static void from_affine(Point& out, const AffinePoint<Field>& a);
		/// @brief Normalize an on-curve Jacobian point while preserving infinity.
		/// @note Subgroup membership is preserved but neither required nor established.
		/// @note Variable-time with respect to the point representation.
		BLST_HD static void from_jacobian(Point& out, const Point& in);
		/// @brief Constant-time conditional swap when @p cbit is one.
		/// @pre @p a and @p b are distinct objects and @p cbit is zero or one.
		BLST_HD static void cswap(Point& a, Point& b, limb_t cbit);
		/// @brief Constant-time conditional copy `a = b` when @p cbit is one.
		/// @pre @p a and @p b are distinct objects and @p cbit is zero or one.
		BLST_HD static void ccopy(Point& a, const Point& b, limb_t cbit);
		/// @brief Constant-time equality of two on-curve point representatives.
		/// @note This is not a curve or subgroup validator.
		BLST_HD static limb_t equal(const Point& p1, const Point& p2);
		/// @brief Convert Jacobian to the variable-time XYZZ representation.
		BLST_HD static void to_xyzz(Pointxyzz<Field>& out, const Point& in);
		/// @brief Prime-order subgroup generator in Jacobian coordinates.
		BLST_HD static const Point& generator();
		/// @brief Return the point at infinity (`Z == 0`).
		BLST_HD static Point infinity();
		/// @brief Return the group identity, equivalent to infinity.
		BLST_HD static Point identity() { return infinity(); }
		/// @brief Return one iff @p p is the group identity.
		BLST_HD static limb_t is_identity(const Point& p) { return is_inf(p); }
		/// @brief Return curve coefficient B for the selected group.
		BLST_HD static const Field& get_B();
		/// @brief Multiply @p in by curve coefficient B.
		BLST_HD static void mul_by_b(Field& out, const Field& in);
		/// @brief Multiply @p in by `4B`, used by doubling formulas.
		BLST_HD static void mul_by_4b(Field& out, const Field& in);
		/// @brief Return one iff on-curve @p p is infinity; not a validator.
		BLST_HD static int is_inf(const Point& p);
		/// @brief Return one iff @p p is on the curve; infinity is accepted.
		BLST_HD static limb_t on_curve(const Point& p);
		/// @brief Return one iff @p p is on-curve and lies in the prime-order subgroup.
		BLST_HD static limb_t in_subgroup(const Point& p);
		/** @name Serialization */
		///@{
		/// @brief Write an on-curve point as uncompressed big-endian affine coordinates.
		/// @note Variable-time; subgroup membership is not required or validated.
		BLST_HD static limb_t serialize_BE(byte* out, const Point& in);
		/// @brief Write an on-curve point in uncompressed IETF form.
		/// @note Subgroup membership is not required or validated.
		BLST_HD static void serialize(byte* out, const Point& in);
		/// @brief Write an on-curve point as compressed big-endian x and return y parity.
		/// @note Variable-time; subgroup membership is not required or validated.
		BLST_HD static limb_t compress_BE(byte* out, const Point& in);
		/// @brief Write an on-curve point in compressed IETF form.
		/// @note Subgroup membership is not required or validated.
		BLST_HD static void compress(byte* out, const Point& in);
		///@}
		/**
		 * @brief Endomorphism-accelerated multiplication for a valid Fr scalar.
		 *
		 * G1 uses the beta endomorphism; G2 uses the psi endomorphism. This is the
		 * GLS/GLV path used by signing and public-key derivation.
		 */
		BLST_HD static void mult_glsv(Point& out, const Point& in, const byte SK[32]);
		/// @brief Variable-time GLV/GLS multiplication for public canonical scalars.
		BLST_HD static void mult_glsv_vartime(Point& out, const Point& in, const byte SK[32]);
		/// @brief Selected fixed-window GLV/GLS constant-time implementation.
		template<int Window>
		BLST_HD static void mult_glsv_window(Point& out, const Point& in, const byte SK[32]);
		/// @brief Selected odd-table wNAF GLV/GLS implementation for public scalars.
		template<int Width>
		BLST_HD static void mult_glsv_wnaf(Point& out, const Point& in, const byte SK[32]);
		/// @brief Compute a BLS signature point `out = SK * in` and normalize it.
		BLST_HD static void sign(Point& out, const Point& in, const byte SK[32]);
		/**
		 * @brief Signature convenience API.
		 * @param out Optional compressed output bytes.
		 * @param sig Optional affine output.
		 * @param hash Hash-to-curve point to sign.
		 * @param SK Little-endian 32-byte secret scalar.
		 */
		BLST_HD static void sign2(byte* out, AffinePoint<Field>* sig,
			const Point& hash, const byte SK[32]);
		/// @brief Derive public key `out = SK * generator()`.
		BLST_HD static void to_pk(Point& out, const byte SK[32]);
		/// @brief Optionally emit compressed and affine public-key outputs.
		BLST_HD static void to_pk2(byte* out, AffinePoint<Field>* pk, const byte SK[32]);
	};

	/**
	 * @tparam Field Either Fp for G1 or Fp2 for G2.
	 * @brief Affine point `(X,Y)` over @p Field.
	 *
	 * Infinity is encoded as all-zero x and y coordinates. Field coordinates are
	 * in Montgomery form. G1 encodings are 96 bytes uncompressed and 48 bytes
	 * compressed; G2 encodings are 192 and 96 bytes respectively. This storage
	 * type does not encode subgroup-validation state.
	 */
	template<typename Field>
	struct AffinePoint
	{
		Field x, y;

		/** @name Uncompressed serialization */
		///@{
		/// @brief Write an on-curve point as `X || Y`; subgroup membership is not required.
		BLST_HD static limb_t serialize_BE(byte* out, const AffinePoint& in);
		/// @brief Write an on-curve point in uncompressed IETF form.
		BLST_HD static void serialize(byte* out, const AffinePoint& in);
		/// @brief Parse raw big-endian `X || Y`, validating range and curve membership.
		/// @post Success establishes on-curve status only, not subgroup membership.
		/// @note The top three x bits are ignored. On failure @p out is unchanged.
		BLST_HD static BLST_ERROR deserialize_BE(AffinePoint& out, const byte* in);
		/// @brief Parse compressed, uncompressed, or infinity-prefixed IETF input.
		/// @post Success establishes canonical encoding and on-curve status only;
		/// subgroup membership must be checked separately before subgroup operations.
		/// @note Infinity is accepted. On failure @p out is unchanged.
		BLST_HD static BLST_ERROR deserialize(AffinePoint& out, const byte* in);
		///@}
		/** @name Compressed serialization */
		///@{
		/// @brief Write on-curve x in big-endian form; subgroup membership is not required.
		BLST_HD static limb_t compress_BE(byte* out, const AffinePoint& in);
		/// @brief Write an on-curve point in compressed IETF form.
		BLST_HD static void compress(byte* out, const AffinePoint& in);
		/// @brief Parse raw big-endian x and reconstruct an on-curve point.
		/// @post Success does not establish subgroup membership.
		/// @note The top three x bits are ignored. On failure @p out is unchanged;
		/// errors are returned as the negated @c BLST_ERROR value in @c limb_t.
		BLST_HD static limb_t uncompress_BE(AffinePoint& out, const byte* in);
		/// @brief Parse and validate IETF compressed form.
		/// @post Success establishes canonical encoding and on-curve status only;
		/// subgroup membership must be checked separately before subgroup operations.
		/// @note Infinity is accepted. On failure @p out is unchanged.
		BLST_HD static BLST_ERROR uncompress(AffinePoint& out, const byte* in);
		///@}
		/// @brief Prime-order subgroup generator in affine coordinates.
		BLST_HD static const AffinePoint& generator();
		/// @brief Constant-time equality of on-curve values; not a validator.
		BLST_HD static int equal(const AffinePoint& p1, const AffinePoint& p2);
		/// @brief Return one iff @p p satisfies the curve equation or is infinity.
		BLST_HD static limb_t on_curve(const AffinePoint& p);
		/// @brief Return one iff an on-curve value is infinity; not a validator.
		BLST_HD static limb_t is_inf(const AffinePoint& p);
		/// @brief Return one iff @p p is on-curve and lies in the prime-order subgroup.
		BLST_HD static limb_t in_subgroup(const AffinePoint& p);
		/// @brief Constant-time conditional negation of an on-curve point.
		/// @note Subgroup membership is not required or validated.
		BLST_HD static void cneg(AffinePoint& point, limb_t c);
	};

	/**
	 * @tparam Field  Either @c Fp (G1) or @c Fp2 (G2).
	 * @brief XYZZ representation for **variable-time** mixed additions.
	 *
	 * Layout: (X, Y, ZZ, ZZZ) where ZZ=Z^2 and ZZZ=Z^3. This allows quick reuse of
	 * Z-powers when chaining many additions. **Not constant-time**; use only in contexts
	 * where timing does not leak secrets (e.g., fixed-base MSM, precomputation).
	 */
	template<typename Field>
	struct Pointxyzz {
		Field x, y, zzz, zz;

		/// Unified add-or-double (variable time). Handles infinities.
		static BLST_HD void dadd(Pointxyzz<Field>& p3, const Pointxyzz<Field>& p1,
			const Pointxyzz<Field>& p2);

		/// Mixed add/sub with affine (variable time). @p subtract flips the sign of @p p2.y.
		static BLST_HD void dadd(Pointxyzz<Field>& p3, const Pointxyzz<Field>& p1,
			const AffinePoint<Field>& p2,
			bool_t subtract);

		/// Convert XYZZ -> Jacobian (rebuild Z from ZZ).
		static BLST_HD void to_Jacobian(Point<Field>& out, const Pointxyzz<Field>& in);

		/// Convert a finite affine point to XYZZ with `ZZ = ZZZ = 1`.
		/// @pre @p in is not the all-zero affine infinity encoding.
		static BLST_HD void from_affine(Pointxyzz<Field>& out, const AffinePoint<Field>& in);

		/// Compare the raw XYZZ representation; this is not mathematical point equality.
		static BLST_HD bool equal(const Pointxyzz<Field>& p1,
			const Pointxyzz<Field>& p2)
		{
			return memcmp(&p1, &p2, sizeof(Pointxyzz<Field>)) == 0;
		}

	};

	/// Aliases for convenience: G1 over Fp, G2 over Fp2.
	using G1Point = Point<Fp>;
	using G1AffinePoint = AffinePoint<Fp>;
	using G1PointZZ = Pointxyzz<Fp>;

	using G2Point = Point<Fp2>;
	using G2AffinePoint = AffinePoint<Fp2>;
	using G2PointZZ = Pointxyzz<Fp2>;


	BLST_HD BLST_INLINE void print(const char* name, const G2Point& p)
	{

		printf("%s:\n", name);
		printf("  X.Re: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.x[0][i]);
		}
		printf("\n  X.Im: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.x[1][i]);
		}
		printf("\n  Y.Re: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.y[0][i]);
		}
		printf("\n  Y.Im: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.y[1][i]);
		}
		printf("\n  Z.Re: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.z[0][i]);
		}
		printf("\n  Z.Im: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.z[1][i]);
		}
		printf("\n");
	}

	BLST_HD BLST_INLINE void print(const char* name, const G1Point& p)
	{

		printf("%s:\n", name);
		printf("  X: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.x[i]);
		}
		printf("\n  Y: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.y[i]);
		}
		printf("\n  Z: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.z[i]);
		}
		printf("\n");
	}

	BLST_HD BLST_INLINE void print(const char* name, const G2AffinePoint& p)
	{

		printf("%s:\n", name);
		printf("  X.Re: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.x[0][i]);
		}
		printf("\n  X.Im: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.x[1][i]);
		}
		printf("\n  Y.Re: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.y[0][i]);
		}
		printf("\n  Y.Im: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.y[1][i]);
		}
		printf("\n");
	}

	BLST_HD BLST_INLINE void print(const char* name, const G1AffinePoint& p)
	{

		printf("%s:\n", name);
		printf("  X: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.x[i]);
		}
		printf("\n  Y: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.y[i]);
		}
		printf("\n");
	}


	// the A' in  y^2 = x ^ 3 + A'*x + B', isogenous one
	BLST_HD const Fp& get_BLST_Aprime_E1(void);

	// the B' in  y^2 = x ^ 3 + A'*x + B', isogenous one
	BLST_HD const Fp& get_BLST_Bprime_E1(void);


	// the A' in  y^2 = x ^ 3 + A'*x + B', isogenous one
	BLST_HD const Fp2& get_BLST_Aprime_E2(void);

	// the B' in  y^2 = x ^ 3 + A'*x + B', isogenous one
	BLST_HD const Fp2& get_BLST_Bprime_E2(void);

	// return the curve parameter B for the G1 curve: y^2 = x^3 + B
	// Will return host or device memory depending on compilation target.
	BLST_HD const Fp& get_BLST_B_E1(void);

	// Return the generator point of G1 in Jacobian coordinates, in Montgomery form.
	// Will return host or device memory depending on compilation target.
	BLST_HD const Point<Fp>& get_BLST_G1(void);

	// Return the negative G1 generator in Jacobian coordinates, in Montgomery form.
	// Will return host or device memory depending on compilation target.
	BLST_HD const G1Point& get_BLST_NEG_G1(void);

	// return the curve parameter B for the G2 curve: y^2 = x^3 + B
	// Will return host or device memory depending on compilation target.
	BLST_HD const Fp2& get_BLST_B_E2(void);

	// return the curve parameter B for : y^2 = x^3 + B
	// Will return host or device memory depending on compilation target.
	template<typename Field>
	BLST_HD const Field& get_B();

	// Return the generator point of G2 in Jacobian coordinates, in Montgomery form.
	// Will return host or device memory depending on compilation target.
	BLST_HD const G2Point& get_BLST_G2(void);

	// Return the negative G2 generator in Jacobian coordinates, in Montgomery form.
	BLST_HD const G2Point& get_BLST_NEG_G2(void);


	/* such that beta^3 - 1 = 0  */
	/* -1/2 * (1 + sqrt(-3)) = ((P-2)^(P-2)) * (1 + (P-3)^((P+1)/4)) */
	/* (0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4
		  897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac << 384) % P */
	BLST_HD const Fp& get_BLST_beta(void);

	BLST_HD void psi(G2Point& out, const G2Point& in);


	template<typename Group>
	BLST_HD void times_minus_z(Group& out, const Group& in);

}
#endif
