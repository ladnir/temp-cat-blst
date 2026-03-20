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
#include <type_traits>
#include "blst/fields/Fr.h"
#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/curve/ec_mult.h"

namespace blst::scalar
{
	template<typename Field>
	struct AffinePoint;


	template<typename Field>
	struct Pointxyzz;



	/**
	 * @tparam Field  Either @c Fp (G1) or @c Fp2 (G2).
	 * @brief Jacobian point (X:Y:Z) over @p Field.
	 *
	 * ## Semantics
	 * - Infinity is represented by @c Z == 0.
	 * - Methods that accept affine inputs treat the all-zero affine value as infinity.
	 *
	 * ## Arithmetic
	 * - @c add / @c add: standard addition (handles infinity).
	 * - @c doubling: point doubling.
	 * - @c dadd / @c dadd: unified add-or-double (handles P==Q and P==±Q) in constant time.
	 * - @c mult / @c mult_glsv: scalar multiplication; @c mult chooses a windowed scheme
	 *   or GLS/GLV decomposition depending on scalar size and group.
	 *
	 * ## Encoding
	 * See file-level notes on serialize/deserialize/compress/uncompress.
	 *
	 * ## Security
	 * Functions are intended to be constant-time unless marked otherwise.
	 * @c cswap/@c ccopy/@c cneg are constant-time conditionals.
	 */
	template<typename Field>
	struct Point
	{
		using Field_type = Field;

		constexpr Point() = default;
		constexpr Point(const Point&) = default;
		constexpr Point& operator=(const Point&) = default;

		/// Construct from Jacobian coordinates. @note Inputs are assumed in Montgomery form.
		constexpr Point(const Field& x, const Field& y, const Field& z) : X(x), Y(y), Z(z) {}

		Field X, Y, Z;


		/**
		 * @brief Unified add-or-double: ` out = p1 (+|=) p2 `.
		 * @param[out] out  Result in Jacobian coords.
		 * @param[in]  p1   Left input in Jacobian coords.
		 * @param[in]  p2   Right input in Jacobian coords.
		 * @param[in]  a4   Optional curve `a` term (for `y^2=x^3+a x + b`). For BLS12-381 `a=0`, pass nullptr.
		 * @note Handles P==Q and P==±Q and points at infinity in constant time.
		 */
		static BLST_HD void dadd(Point& out, const Point& p1, const Point& p2, const Field* a4 = nullptr);


		/**
		 * @brief Unified add-or-double with affine addend.
		 * @details Treats @p p2==0 (all-zero affine) as infinity.
		 * @note Handles P==Q and P==±Q and points at infinity in constant time.
		 */
		static BLST_HD void dadd(Point<Field>& out, const Point<Field>& p1,
			const AffinePoint<Field>& p2);

		/**
		 * @brief Add @p p1 and @p p2 (handles infinity).
		 * @note Uses EFD add-2007-bl (Jacobian).
		 * @note Requires @p p1 and @p p2 to be distinct and non-opposite for
		 */
		static BLST_HD void add(Point& out, const Point& p1, const Point& p2);


		/**
		 * @brief Mixed add (Jacobian + affine) (handles infinity).
		 * @note Uses EFD madd-2007-bl.
		 * @note Requires @p p1 and @p p2 to be distinct and non-opposite for
		 */
		static BLST_HD void add(Point& out, const Point& p1,
			const AffinePoint<Field>& p2);

		/**
		 * @brief Doubling ` out = 2 * p1 `.
		 * @note Uses EFD dbl-2009-l.
		 */
		static BLST_HD void doubling(Point& out, const Point& p1);


		/** @brief Repeated doubling: ` out = 2^n * p `.
		 */
		BLST_HD static void add_n_dbl(Point& out, const Point& p, size_t n);


		/**
		 * @brief Variable-time scalar multiplication with auto-selection.
		 * @param[out] out     Result.
		 * @param[in]  a       Base point (Jacobian).
		 * @param[in]  scalar  Little-endian scalar bytes.
		 * @param[in]  nbits   Number of meaningful bits in @p scalar.
		 * @details
		 *  - For short scalars, uses windowed w=4.
		 *  - For up to 256 bits and valid Fr range, uses GLS/GLV acceleration (@c mult_glsv).
		 *  - Otherwise falls back to windowed w=5.
		 */
		BLST_HD static void mult(
			Point& out,
			const Point& a,
			const byte* scalar,
			size_t nbits);

		/// Windowed mul (w=5). Primarily used internally or for very large scalars.
		static BLST_HD void mult_w5(Point& out, const Point& point,
			const byte* scalar, size_t nbits);

		/// Conditional negate Y in constant time if @p cbit is 1.
		static BLST_HD void cneg(Point& p, limb_t cbit);

		/// negate Y in constant time.
		static BLST_HD void neg(Point& out, const Point& p);

		/// Convert Jacobian -> affine. Infinity maps to the all-zero affine encoding.
		static BLST_HD void to_affine(AffinePoint<Field>& out, const Point& in);

		/// Convert affine -> Jacobian. Infinity maps to Z=0; otherwise Z=1.
		static BLST_HD void from_affine(Point& out, const AffinePoint<Field>& a);

		/// Normalize Jacobian to affine in @p out (computes X/Z^2, Y/Z^3). Infinity preserved.
		static BLST_HD void from_jacobian(Point& out, const Point& in);

		/// Constant-time conditional swap: if @p cbit==1, swap a and b.
		static BLST_HD void cswap(Point& BLST_restrict a,
			Point& BLST_restrict b, bool_t cbit);

		/// Constant-time conditional copy: if @p cbit==1, a <- b.
		static BLST_HD void ccopy(Point& BLST_restrict a,
			const Point& BLST_restrict b, bool_t cbit);

		/// Constant-time Jacobian equality check (compares affine representatives).
		static BLST_HD limb_t equal(const Point& p1, const Point& p2);

		/// Convert Jacobian -> XYZZ (non-normalized Z-variants) for fast mixed adds.
		static BLST_HD BLST_INLINE void to_xyzz(Pointxyzz<Field>& out, const Point<Field>& in);

		/// Group generator in Jacobian coordinates (Montgomery field limbs).
		static BLST_HD const Point& generator();

		/// returns the point at infinity (Z == 0).
		static BLST_HD Point infinity();

		// returns the identity point (same as infinity).
		static BLST_HD Point identity() { return infinity(); }

		static BLST_HD bool_t is_identity(const Point& p) { return is_inf(p); }

		static BLST_HD const Field& get_B();

		/// Multiply by curve B for G1/G2 (specialized). For BLS12-381, B = 4 (G1) or (4+4i) (G2).
		BLST_HD static void mul_by_b(Field& out, const Field& in);

		/// Multiply by 4B (used in doubling paths).
		BLST_HD static void mul_by_4b(Field& out, const Field& in);

		/// @return 1 if @p p is infinity (Z == 0), else 0.
		BLST_HD static int is_inf(const Point& p);

		/// @return 1 if on the curve (including infinity), else 0.
		BLST_HD static bool_t on_curve(const Point& p);

		/// @return 1 if in the subgroup, else 0.
		BLST_HD static bool_t in_subgroup(const Point& p);


		/** @name Serialization (big-endian) */
		///@{
		/// Uncompressed big-endian of affine representative. Returns sgn0 parity (see compress()).
		BLST_HD static limb_t serialize_BE(unsigned char* out, const Point& in);
		/// Uncompressed with infinity bit support. If infinity, sets 0x40 and zeroes.
		BLST_HD static void serialize(unsigned char* out, const Point& in);
		/// Compressed big-endian X and returns sgn0 parity (no flags set in @p out).
		BLST_HD static limb_t compress_BE(unsigned char* out, const Point& in);
		/// IETF-style compressed: sets 0x80; encodes y-parity; if infinity sets 0xC0 and zeroes.
		BLST_HD static void compress(unsigned char* out, const Point& in);
		///@}

		/**
		 * @brief Endomorphism-accelerated multiplication (GLV/GLS).
		 * @details
		 * - G1 uses the cube root of unity ` \beta ` (get_BLST_beta()) to split scalars.
		 * - G2 uses the Frobenius-like map ` \psi `.
		 * @note Intended for scalars within Fr (checked by @c Fr::check_mod).
		 *
		 * @note Galbraith-Lin-Scott, ~67% faster than G2Point_mul_w5
		 * @note Gallant-Lambert-Vanstone, ~45% faster than mult_w<G1Point,5> 
		 */
		BLST_HD static void mult_glsv(Point& out, const Point& in, const pow256 SK);

		/**
		 * @brief BLS-style signature mapping: ` out = SK * in ` in Jacobian form.
		 * @note Internally calls @c mult_glsv and normalizes to affine to avoid bias in Z.
		 */
		BLST_HD static void sign(Point& out, const Point& in, const pow256 SK);

		/**
		 * @brief Convenience signature API.
		 * @param[out] out Optional compressed bytes; first byte pre-decorated with sgn0 & infinity bits.
		 * @param[out] sig Optional affine signature.
		 * @param[in]  hash Point to sign (typically hash-to-curve result).
		 * @param[in]  SK   Little-endian 32-byte secret scalar.
		 */
		BLST_HD static void sign2(
			unsigned char* out /*optional output*/,
			AffinePoint<Field>* sig /*optional output*/,
			const Point& hash,
			const pow256 SK);

		/// Derive public key ` PK = SK * G ` into @p out.
		BLST_HD static void to_pk(Point& out, const pow256 SK);

		/// Variant that can also emit compressed bytes and/or affine PK (both optional).
		BLST_HD static void to_pk2(unsigned char* out, AffinePoint<Field>* PK, const pow256 SK);


	};
	static_assert(std::is_trivially_constructible_v<Point<Fp>>, "Point<Fp> must be an aggregate type");


	/**
	 * @tparam Field  Either @c Fp (G1) or @c Fp2 (G2).
	 * @brief Affine point (X,Y) over @p Field.
	 *
	 * - Infinity is encoded as all zero limbs for both X and Y.
	 * - Used for I/O and mixed additions. Most arithmetic is on @c Point (Jacobian).
	 *
	 * ## Encoding sizes (big-endian)
	 * - G1: 96 bytes uncompressed, 48 bytes compressed
	 * - G2: 192 bytes uncompressed, 96 bytes compressed
	 */
	template<typename Field>
	struct AffinePoint
	{
		Field X, Y;



		/** @name Uncompressed (big-endian) */
		///@{
		/// Writes X||Y in big-endian Montgomery form; returns sgn0 parity.
		BLST_HD static limb_t serialize_BE(
			unsigned char* out,
			const AffinePoint& in);
		/// Like serialize_BE but writes infinity as 0x40 || zeros (IETF style).
		BLST_HD static void serialize(unsigned char* out,
			const AffinePoint& in);
		/// Parses big-endian X||Y. Validates < p and on-curve. Outputs Montgomery.
		BLST_HD static BLST_ERROR deserialize_BE(AffinePoint& out,
			const unsigned char* in);
		/// Auto-detects compressed/uncompressed/infinity prefixes and dispatches parse.
		BLST_HD static BLST_ERROR deserialize(AffinePoint& out,
			const unsigned char* in);
		///@}


		/** @name Compressed (IETF style) */
		///@{
		/// Writes big-endian X and returns sgn0 parity (no header bits set).
		BLST_HD static limb_t compress_BE(unsigned char* out, const AffinePoint& in);
		/// Writes IETF compressed: sets 0x80; writes parity; or 0xC0 for infinity.
		BLST_HD static void compress(unsigned char* out, const AffinePoint& in);
		/// Parses big-endian X, reconstructs Y via curve equation; returns sgn0 parity.
		BLST_HD static limb_t uncompress_BE(AffinePoint& out,
			const unsigned char* in);
		/// Parses IETF compressed form with header bits; validates on-curve.
		BLST_HD static BLST_ERROR uncompress(AffinePoint& out,
			const unsigned char* in);
		///@}

		/// Group generator (affine). Backed by host/device constants.
		BLST_HD static const AffinePoint& generator();

		/// Constant-time equality test on raw limbs.
		BLST_HD static int equal(const AffinePoint& p1,
			const AffinePoint& p2);

		/// True if affine point satisfies curve equation (or is infinity).
		BLST_HD static bool_t on_curve(const AffinePoint& p);

		/// True if (X,Y) are all-zero limbs (infinity).
		BLST_HD static bool_t is_inf(const AffinePoint& p);

		BLST_HD static bool_t in_subgroup(const AffinePoint<Field>& p);

		BLST_HD static void cneg(AffinePoint& inOut, bool_t c);

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
		Field X, Y, ZZZ, ZZ;

		/// Unified add-or-double (variable time). Handles infinities.
		static BLST_HD void dadd(Pointxyzz<Field>& p3, const Pointxyzz<Field>& p1,
			const Pointxyzz<Field>& p2);

		/// Mixed add/sub with affine (variable time). @p subtract flips the sign of @p p2.Y.
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
			printf("%016llx ", (unsigned long long)p.X[0][i]);
		}
		printf("\n  X.Im: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.X[1][i]);
		}
		printf("\n  Y.Re: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.Y[0][i]);
		}
		printf("\n  Y.Im: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.Y[1][i]);
		}
		printf("\n  Z.Re: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.Z[0][i]);
		}
		printf("\n  Z.Im: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.Z[1][i]);
		}
		printf("\n");
	}

	BLST_HD BLST_INLINE void print(const char* name, const G1Point& p)
	{

		printf("%s:\n", name);
		printf("  X: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.X[i]);
		}
		printf("\n  Y: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.Y[i]);
		}
		printf("\n  Z: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.Z[i]);
		}
		printf("\n");
	}

	BLST_HD BLST_INLINE void print(const char* name, const G2AffinePoint& p)
	{

		printf("%s:\n", name);
		printf("  X.Re: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.X[0][i]);
		}
		printf("\n  X.Im: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.X[1][i]);
		}
		printf("\n  Y.Re: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.Y[0][i]);
		}
		printf("\n  Y.Im: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.Y[1][i]);
		}
		printf("\n");
	}

	BLST_HD BLST_INLINE void print(const char* name, const G1AffinePoint& p)
	{

		printf("%s:\n", name);
		printf("  X: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.X[i]);
		}
		printf("\n  Y: ");
		for (int i = NLIMBS(384) - 1; i >= 0; i--) {
			printf("%016llx ", (unsigned long long)p.Y[i]);
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
