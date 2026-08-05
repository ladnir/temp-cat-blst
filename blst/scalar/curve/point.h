/*
 * Copyright Supranational LLC
 * Copyright Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLST_SCALAR_POINT_H__
#define __BLST_SCALAR_POINT_H__


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
  * Unless explicitly stated, operations are written to be branchless/constant-time
  * w.r.t. secret data. The @c Pointxyzz helpers are explicitly marked
  * **NOT CONSTANT TIME**; use only where timing does not leak secrets (e.g., fixed-base MSMs).
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
		BLST_HD static void dadd(Point& out, const Point& p1, const Point& p2,
			const Field* a4 = nullptr);
		/// @brief Unified constant-time add-or-double with an affine addend.
		/// @note Treats the all-zero affine value as infinity.
		BLST_HD static void dadd(Point& out, const Point& p1, const AffinePoint<Field>& p2);
		/// @brief Add two distinct, nonopposite Jacobian points, handling infinity.
		/// @note Uses the EFD add-2007-bl formula.
		BLST_HD static void add(Point& out, const Point& p1, const Point& p2);
		/// @brief Mixed Jacobian-plus-affine addition, handling infinity.
		/// @note Uses the EFD madd-2007-bl formula and requires distinct, nonopposite inputs.
		BLST_HD static void add(Point& out, const Point& p1, const AffinePoint<Field>& p2);
		/// @brief Point doubling using the EFD dbl-2009-l formula.
		BLST_HD static void doubling(Point& out, const Point& p1);
		/// @brief Repeated doubling: `out = 2^n * p`.
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
		 */
		BLST_HD static void mult(Point& out, const Point& a, const byte* scalar, size_t nbits);
		/**
		 * @brief Optimized variable-time scalar multiplication.
		 * @warning Scalar range and window values may affect control flow and memory
		 * access. Use only with public scalar data. The input must already be an
		 * on-curve prime-order subgroup point; Release builds do not validate it.
		 * @p nbits must not exceed 256, and @p scalar may be null only for zero bits.
		 */
		BLST_HD static void mult_vartime(Point& out, const Point& a, const byte* scalar, size_t nbits);
		/**
		 * @brief Constant-time generic multiplication without subgroup endomorphisms.
		 * @warning The input must be on curve, but need not be in the prime-order
		 * subgroup. Release builds do not validate that precondition. @p nbits must
		 * not exceed 256, and @p scalar may be null only for zero bits.
		 */
		BLST_HD static void unchecked_mult(Point& out, const Point& a, const byte* scalar, size_t nbits);
		/// @brief Width-five windowed multiplication for internal or large-scalar use.
		BLST_HD static void mult_w5(Point& out, const Point& point, const byte* scalar, size_t nbits);
		/// @brief Constant-time conditional y-coordinate negation when @p cbit is one.
		BLST_HD static void cneg(Point& p, limb_t cbit);
		/// @brief Negate the y coordinate.
		BLST_HD static void neg(Point& out, const Point& p);
		/// @brief Convert Jacobian to affine; infinity maps to the all-zero encoding.
		BLST_HD static void to_affine(AffinePoint<Field>& out, const Point& in);
		/// @brief Convert affine to Jacobian; finite points receive `Z = 1`.
		BLST_HD static void from_affine(Point& out, const AffinePoint<Field>& a);
		/// @brief Normalize a Jacobian point while preserving infinity.
		BLST_HD static void from_jacobian(Point& out, const Point& in);
		/// @brief Constant-time conditional swap when @p cbit is one.
		BLST_HD static void cswap(Point& a, Point& b, limb_t cbit);
		/// @brief Constant-time conditional copy `a = b` when @p cbit is one.
		BLST_HD static void ccopy(Point& a, const Point& b, limb_t cbit);
		/// @brief Constant-time equality of affine representatives.
		BLST_HD static limb_t equal(const Point& p1, const Point& p2);
		/// @brief Convert Jacobian to the variable-time XYZZ representation.
		BLST_HD static void to_xyzz(Pointxyzz<Field>& out, const Point& in);
		/// @brief Group generator in Jacobian coordinates.
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
		/// @brief Return one iff @p p is infinity.
		BLST_HD static int is_inf(const Point& p);
		/// @brief Return one iff @p p is on the curve; infinity is accepted.
		BLST_HD static limb_t on_curve(const Point& p);
		/// @brief Return one iff @p p is on-curve and lies in the prime-order subgroup.
		BLST_HD static limb_t in_subgroup(const Point& p);
		/** @name Serialization */
		///@{
		/// @brief Write uncompressed big-endian affine coordinates and return sgn0 parity.
		BLST_HD static limb_t serialize_BE(byte* out, const Point& in);
		/// @brief Write uncompressed IETF form, using `0x40 || zeros` for infinity.
		BLST_HD static void serialize(byte* out, const Point& in);
		/// @brief Write compressed big-endian x without header flags and return y parity.
		BLST_HD static limb_t compress_BE(byte* out, const Point& in);
		/// @brief Write IETF compressed form with compression, sign, and infinity flags.
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
	 * compressed; G2 encodings are 192 and 96 bytes respectively.
	 */
	template<typename Field>
	struct AffinePoint
	{
		Field x, y;

		/** @name Uncompressed serialization */
		///@{
		/// @brief Write `X || Y` in big-endian form and return sgn0 parity.
		BLST_HD static limb_t serialize_BE(byte* out, const AffinePoint& in);
		/// @brief Write uncompressed IETF form, including the infinity encoding.
		BLST_HD static void serialize(byte* out, const AffinePoint& in);
		/// @brief Parse raw big-endian `X || Y`, validating range and curve membership.
		BLST_HD static BLST_ERROR deserialize_BE(AffinePoint& out, const byte* in);
		/// @brief Parse compressed, uncompressed, or infinity-prefixed IETF input.
		BLST_HD static BLST_ERROR deserialize(AffinePoint& out, const byte* in);
		///@}
		/** @name Compressed serialization */
		///@{
		/// @brief Write big-endian x without header flags and return y parity.
		BLST_HD static limb_t compress_BE(byte* out, const AffinePoint& in);
		/// @brief Write IETF compressed form with compression, sign, and infinity flags.
		BLST_HD static void compress(byte* out, const AffinePoint& in);
		/// @brief Parse raw big-endian x, reconstruct y, and return sgn0 parity.
		BLST_HD static limb_t uncompress_BE(AffinePoint& out, const byte* in);
		/// @brief Parse and validate IETF compressed form.
		BLST_HD static BLST_ERROR uncompress(AffinePoint& out, const byte* in);
		///@}
		/// @brief Group generator in affine coordinates.
		BLST_HD static const AffinePoint& generator();
		/// @brief Constant-time equality on raw limbs.
		BLST_HD static int equal(const AffinePoint& p1, const AffinePoint& p2);
		/// @brief Return one iff @p p satisfies the curve equation or is infinity.
		BLST_HD static limb_t on_curve(const AffinePoint& p);
		/// @brief Return one iff both coordinates are all zero.
		BLST_HD static limb_t is_inf(const AffinePoint& p);
		/// @brief Return one iff @p p is on-curve and lies in the prime-order subgroup.
		BLST_HD static limb_t in_subgroup(const AffinePoint& p);
		/// @brief Constant-time conditional y-coordinate negation.
		BLST_HD static void cneg(AffinePoint& point, limb_t c);
	};

	/**
	 * @tparam Field  Either @c Fp (G1) or @c Fp2 (G2).
	 * @brief XYZZ representation for **variable-time** mixed additions.
	 *
	 * Layout: (X, Y, ZZ, ZZZ) where ZZ=Z^2 and ZZZ=Z^3. This allows quick reuse of
	 * Z-powers when chaining many additions. **Not constant-time**—use only in contexts
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

		static BLST_HD void from_affine(Pointxyzz<Field>& out, const AffinePoint<Field>& in);

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

	// return the generator point of G1 in affine coordinates, in montgomery form.
	// Will return host or device memory depending on compilation target.
	BLST_HD const Point<Fp>& get_BLST_G1(void);

	// return the negative of the generator point of G1 in affine coordinates, in montgomery form.
	// Will return host or device memory depending on compilation target.
	BLST_HD const G1Point& get_BLST_NEG_G1(void);

	// return the curve parameter B for the G2 curve: y^2 = x^3 + B
	// Will return host or device memory depending on compilation target.
	BLST_HD const Fp2& get_BLST_B_E2(void);

	// return the curve parameter B for : y^2 = x^3 + B
	// Will return host or device memory depending on compilation target.
	template<typename Field>
	BLST_HD const Field& get_B();

	// return the generator point of G2 in affine coordinates, in montgomery form.
	// Will return host or device memory depending on compilation target.
	BLST_HD const G2Point& get_BLST_G2(void);

	// return the negative of the generator point of G2 in affine coordinates, in montgomery form.
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
