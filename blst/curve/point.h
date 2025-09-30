/*
 * Copyright Supranational LLC
 * Copyright Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __BLS12_381_ASM_POINT_H__
#define __BLS12_381_ASM_POINT_H__


 /**
  * @file
  * @brief BLS12-381 elliptic-curve point types and operations over Fp and Fp2.
  *
  * # Curves and fields
  * - **G1** is a short-Weierstrass curve over the base field @c Fp with model
  *   @f$ y^2 = x^3 + B @f$ (parameter B returned by get_BLST_B_E1()).
  * - **G2** is the twist over @c Fp2 with model
  *   @f$ y^2 = x^3 + B @f$ (parameter B returned by get_BLST_B_E2()).
  *
  * # Coordinates & representations
  * - @c Point<Field> stores Jacobian coordinates (X:Y:Z) with the point at
  *   infinity encoded as @c Z == 0. @c Field is either @c Fp (G1) or @c Fp2 (G2).
  * - @c Point_affine<Field> stores affine (X,Y). Infinity is encoded as the
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
#include "ec_mult.h"

namespace blst
{
	template<typename Field>
	struct Point_affine;


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
	 * - @c add / @c add_affine: standard addition (handles infinity).
	 * - @c doubling: point doubling.
	 * - @c dadd / @c dadd_affine: unified add-or-double (handles P==Q and P==±Q) in constant time.
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
		 * @brief Unified add-or-double: @f$ out = p1 (+|=) p2 @f$.
		 * @param[out] out  Result in Jacobian coords.
		 * @param[in]  p1   Left input in Jacobian coords.
		 * @param[in]  p2   Right input in Jacobian coords.
		 * @param[in]  a4   Optional curve @f$a@f$ term (for @f$y^2=x^3+a x + b@f$). For BLS12-381 @f$a=0@f$, pass nullptr.
		 * @note Handles P==Q and P==±Q and points at infinity in constant time.
		 */
		static BLST_HD void dadd(Point& out, const Point& p1, const Point& p2, const Field* a4 = nullptr);


		/**
		 * @brief Unified add-or-double with affine addend.
		 * @details Treats @p p2==0 (all-zero affine) as infinity.
		 * @note Handles P==Q and P==±Q and points at infinity in constant time.
		 */
		static BLST_HD void dadd_affine(Point<Field>& out, const Point<Field>& p1,
			const Point_affine<Field>& p2);

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
		static BLST_HD void add_affine(Point& out, const Point& p1,
			const Point_affine<Field>& p2);

		/**
		 * @brief Doubling @f$ out = 2 * p1 @f$.
		 * @note Uses EFD dbl-2009-l.
		 */
		static BLST_HD void doubling(Point& out, const Point& p1);


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

		/// Convert Jacobian -> affine. Infinity maps to the all-zero affine encoding.
		static BLST_HD void to_affine(Point_affine<Field>& out, const Point& in);

		/// Convert affine -> Jacobian. Infinity maps to Z=0; otherwise Z=1.
		static BLST_HD void from_affine(Point& out, const Point_affine<Field>& a);

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

		/// Multiply by curve B for G1/G2 (specialized). For BLS12-381, B = 4 (G1) or (4+4i) (G2).
		BLST_HD static void mul_by_b(Field& out, const Field& in);

		/// Multiply by 4B (used in doubling paths).
		BLST_HD static void mul_by_4b(Field& out, const Field& in);

		/// @return 1 if @p p is infinity (Z == 0), else 0.
		BLST_HD static int is_inf(const Point& p);

		/// @return 1 if on the curve (including infinity), else 0.
		BLST_HD static bool_t on_curve(const Point& p);

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
		 * - G1 uses the cube root of unity @f$ \beta @f$ (get_BLST_beta()) to split scalars.
		 * - G2 uses the Frobenius-like map @f$ \psi @f$.
		 * @note Intended for scalars within Fr (checked by @c Fr::check_mod).
		 *
		 * @note Galbraith-Lin-Scott, ~67% faster than POINTonE2_mul_w5
		 * @note Gallant-Lambert-Vanstone, ~45% faster than mult_w<POINTonE1,5> 
		 */
		BLST_HD static void mult_glsv(Point& out, const Point& in, const pow256 SK);

		/**
		 * @brief BLS-style signature mapping: @f$ out = SK * in @f$ in Jacobian form.
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
			Point_affine<Field>* sig /*optional output*/,
			const Point& hash,
			const pow256 SK);

		/// Derive public key @f$ PK = SK * G @f$ into @p out.
		BLST_HD static void to_pk(Point& out, const pow256 SK);

		/// Variant that can also emit compressed bytes and/or affine PK (both optional).
		BLST_HD static void to_pk2(unsigned char* out, Point_affine<Field>* PK, const pow256 SK);


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
	struct Point_affine
	{
		Field X, Y;



		/** @name Uncompressed (big-endian) */
		///@{
		/// Writes X||Y in big-endian Montgomery form; returns sgn0 parity.
		BLST_HD static limb_t serialize_BE(
			unsigned char* out,
			const Point_affine& in);
		/// Like serialize_BE but writes infinity as 0x40 || zeros (IETF style).
		BLST_HD static void serialize(unsigned char* out,
			const Point_affine& in);
		/// Parses big-endian X||Y. Validates < p and on-curve. Outputs Montgomery.
		BLST_HD static BLST_ERROR deserialize_BE(Point_affine& out,
			const unsigned char* in);
		/// Auto-detects compressed/uncompressed/infinity prefixes and dispatches parse.
		BLST_HD static BLST_ERROR deserialize(Point_affine& out,
			const unsigned char* in);
		///@}


		/** @name Compressed (IETF style) */
		///@{
		/// Writes big-endian X and returns sgn0 parity (no header bits set).
		BLST_HD static limb_t compress_BE(unsigned char* out, const Point_affine& in);
		/// Writes IETF compressed: sets 0x80; writes parity; or 0xC0 for infinity.
		BLST_HD static void compress(unsigned char* out, const Point_affine& in);
		/// Parses big-endian X, reconstructs Y via curve equation; returns sgn0 parity.
		BLST_HD static limb_t uncompress_BE(Point_affine& out,
			const unsigned char* in);
		/// Parses IETF compressed form with header bits; validates on-curve.
		BLST_HD static BLST_ERROR uncompress(Point_affine& out,
			const unsigned char* in);
		///@}

		/// Group generator (affine). Backed by host/device constants.
		BLST_HD static const Point_affine& generator();

		/// Constant-time equality test on raw limbs.
		BLST_HD static int equal(const Point_affine& p1,
			const Point_affine& p2);

		/// True if affine point satisfies curve equation (or is infinity).
		BLST_HD static bool_t on_curve(const Point_affine& p);

		/// True if (X,Y) are all-zero limbs (infinity).
		BLST_HD static bool_t is_inf(const Point_affine& p);

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
		static BLST_HD void dadd_affine(Pointxyzz<Field>& p3, const Pointxyzz<Field>& p1,
			const Point_affine<Field>& p2,
			bool_t subtract);

		/// Convert XYZZ -> Jacobian (rebuild Z from ZZ).
		static BLST_HD void to_Jacobian(Point<Field>& out, const Pointxyzz<Field>& in);

	};

	/// Aliases for convenience: G1 over Fp, G2 over Fp2.
	using POINTonE1 = Point<Fp>;
	using POINTonE1_affine = Point_affine<Fp>;
	using POINTonE2 = Point<Fp2>;
	using POINTonE2_affine = Point_affine<Fp2>;


	BLST_HD BLST_INLINE void print(const char* name, const POINTonE2& p)
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

	BLST_HD BLST_INLINE void print(const char* name, const POINTonE1& p)
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

	BLST_HD BLST_INLINE void print(const char* name, const POINTonE2_affine& p)
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

	BLST_HD BLST_INLINE void print(const char* name, const POINTonE1_affine& p)
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


	static_assert(std::is_trivially_constructible_v<POINTonE1>, "POINTonE1 must be an aggregate type");

#ifdef __GNUC__
# pragma GCC diagnostic ignored "-Wunused-function"
#endif



#define BLST_B_E1_LIMBS \
    TO_LIMB_T(0xaa270000000cfff3), TO_LIMB_T(0x53cc0032fc34000a),\
    TO_LIMB_T(0x478fe97a6b0a807f), TO_LIMB_T(0xb1d37ebee6ba24d7),\
    TO_LIMB_T(0x8ec9733bbf78ab2f), TO_LIMB_T(0x09d645513d83de7e)

#define BLST_G1_LIMBS                                                     \
 /* generator point [in Montgomery] */                                    \
    /* (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905  */            \
    /*    a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb << 384) % P */ \
Fp{                                                                         \
    TO_LIMB_T(0x5cb38790fd530c16), TO_LIMB_T(0x7817fc679976fff5),         \
        TO_LIMB_T(0x154f95c7143ba1c1), TO_LIMB_T(0xf0ae6acdf3d0e747),     \
        TO_LIMB_T(0xedce6ecc21dbf440), TO_LIMB_T(0x120177419e0bfb75)      \
},                                                                        \
/* (0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af6             */     \
/*    00db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1 << 384) % P */     \
Fp{ TO_LIMB_T(0xbaac93d50ce72271), TO_LIMB_T(0x8c22631a7918fd8e),           \
  TO_LIMB_T(0xdd595f13570725ce), TO_LIMB_T(0x51ac582950405194),           \
  TO_LIMB_T(0x0e1c8c3fad0059c0), TO_LIMB_T(0x0bbc3efc5008a26a) },         \
Fp{ ONE_MONT_P_LIMBS }


#define BLST_NEG_G1_LIMBS                                                    \
/* negative generator [in Montgomery] */                                     \
    /* (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905*/                 \
     /*    a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb << 384) % P */   \
Fp{                                                                            \
    TO_LIMB_T(0x5cb38790fd530c16), TO_LIMB_T(0x7817fc679976fff5),            \
        TO_LIMB_T(0x154f95c7143ba1c1), TO_LIMB_T(0xf0ae6acdf3d0e747),        \
        TO_LIMB_T(0xedce6ecc21dbf440), TO_LIMB_T(0x120177419e0bfb75)         \
},                                                                           \
/* (0x114d1d6855d545a8aa7d76c8cf2e21f267816aef1db507c9*/                     \
/*    6655b9d5caac42364e6f38ba0ecb751bad54dcd6b939c2ca << 384) % P */        \
Fp{ TO_LIMB_T(0xff526c2af318883a), TO_LIMB_T(0x92899ce4383b0270),              \
  TO_LIMB_T(0x89d7738d9fa9d055), TO_LIMB_T(0x12caf35ba344c12a),              \
  TO_LIMB_T(0x3cff1b76964b5317), TO_LIMB_T(0x0e44d2ede9774430) },            \
Fp{ ONE_MONT_P_LIMBS }




	/*
	 * y^2 = x^3 + B
	 */
#define BLST_B_E2_LIMBS                                               \
  { TO_LIMB_T(0xaa270000000cfff3), TO_LIMB_T(0x53cc0032fc34000a),     \
    TO_LIMB_T(0x478fe97a6b0a807f), TO_LIMB_T(0xb1d37ebee6ba24d7),     \
    TO_LIMB_T(0x8ec9733bbf78ab2f), TO_LIMB_T(0x09d645513d83de7e) },   \
  { TO_LIMB_T(0xaa270000000cfff3), TO_LIMB_T(0x53cc0032fc34000a),     \
    TO_LIMB_T(0x478fe97a6b0a807f), TO_LIMB_T(0xb1d37ebee6ba24d7),     \
    TO_LIMB_T(0x8ec9733bbf78ab2f), TO_LIMB_T(0x09d645513d83de7e) }

	 /* generator point [in Montgomery] */
#define BLST_G2_LIMBS                                                      \
Fp2{ /* (0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02            */     \
/*      b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8 << 384) % P */    \
  Fp{ TO_LIMB_T(0xf5f28fa202940a10), TO_LIMB_T(0xb3f5fb2687b4961a),          \
    TO_LIMB_T(0xa1a893b53e2ae580), TO_LIMB_T(0x9894999d1a3caee9),          \
    TO_LIMB_T(0x6f67b7631863366b), TO_LIMB_T(0x058191924350bcd7) },        \
    /* (0x13e02b6052719f607dacd3a088274f65596bd0d09920b61a            */   \
    /*  b5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e << 384) % P */    \
   Fp{ TO_LIMB_T(0xa5a9c0759e23f606), TO_LIMB_T(0xaaa0c59dbccd60c3),        \
      TO_LIMB_T(0x3bb17e18e2867806), TO_LIMB_T(0x1b1ab6cc8541b367),        \
      TO_LIMB_T(0xc2b6ed0ef2158547), TO_LIMB_T(0x11922a097360edf3) }       \
  },                                                                       \
  Fp2{ /* (0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a7            */   \
    /* 6d429a695160d12c923ac9cc3baca289e193548608b82801 << 384) % P */     \
    Fp{ TO_LIMB_T(0x4c730af860494c4a), TO_LIMB_T(0x597cfa1f5e369c5a),        \
      TO_LIMB_T(0xe7e6856caa0a635a), TO_LIMB_T(0xbbefb5e96e0d495f),        \
      TO_LIMB_T(0x07d3a975f0ef25a2), TO_LIMB_T(0x0083fd8e7e80dae5) },      \
      /* (0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af          */   \
    /*    267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be << 384) % P */  \
      Fp{ TO_LIMB_T(0xadc0fc92df64b05d), TO_LIMB_T(0x18aa270a2b1461dc),      \
        TO_LIMB_T(0x86adac6a3be4eba0), TO_LIMB_T(0x79495c4ec93da33a),      \
        TO_LIMB_T(0xe7175850a43ccaed), TO_LIMB_T(0x0b2bc2a163de1bf2) },    \
    },                                                                     \
    Fp2{ Fp{ ONE_MONT_P_LIMBS }, Fp{ 0 } }

/* negative generator [in Montgomery] */
#define BLST_NEG_G2_LIMBS                                                   \
Fp2{ /* (0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02             */     \
  /*    b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8 << 384) % P  */    \
  Fp{ TO_LIMB_T(0xf5f28fa202940a10), TO_LIMB_T(0xb3f5fb2687b4961a),           \
    TO_LIMB_T(0xa1a893b53e2ae580), TO_LIMB_T(0x9894999d1a3caee9),           \
    TO_LIMB_T(0x6f67b7631863366b), TO_LIMB_T(0x058191924350bcd7) },         \
    /* (0x13e02b6052719f607dacd3a088274f65596bd0d09920b61a             */   \
    /*    b5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e << 384) % P */   \
    Fp{ TO_LIMB_T(0xa5a9c0759e23f606), TO_LIMB_T(0xaaa0c59dbccd60c3),         \
      TO_LIMB_T(0x3bb17e18e2867806), TO_LIMB_T(0x1b1ab6cc8541b367),         \
      TO_LIMB_T(0xc2b6ed0ef2158547), TO_LIMB_T(0x11922a097360edf3) }        \
  },                                                                        \
  Fp2{ /* (0x0d1b3cc2c7027888be51d9ef691d77bcb679afda66c73f17             */   \
    /*    f9ee3837a55024f78c71363275a75d75d86bab79f74782aa << 384) % P */   \
    Fp{ TO_LIMB_T(0x6d8bf5079fb65e61), TO_LIMB_T(0xc52f05df531d63a5),         \
      TO_LIMB_T(0x7f4a4d344ca692c9), TO_LIMB_T(0xa887959b8577c95f),         \
      TO_LIMB_T(0x4347fe40525c8734), TO_LIMB_T(0x197d145bbaff0bb5) },       \
      /* (0x13fa4d4a0ad8b1ce186ed5061789213d993923066dddaf10             */ \
      /*    40bc3ff59f825c78df74f2d75467e25e0f55f8a00fa030ed << 384) % P */ \
      Fp{ TO_LIMB_T(0x0c3e036d209afa4e), TO_LIMB_T(0x0601d8f4863f9e23),       \
        TO_LIMB_T(0xe0832636bacc0a84), TO_LIMB_T(0xeb2def362a476f84),       \
        TO_LIMB_T(0x64044f659f0ee1e9), TO_LIMB_T(0x0ed54f48d5a1caa7) }      \
    },                                                                      \
    Fp2{ Fp{ ONE_MONT_P_LIMBS }, Fp{ 0 } }                                         


/* such that beta^3 - 1 = 0  */
#define BLST_beta_LIMBS                                                    \
    /* -1/2 * (1 + sqrt(-3)) = ((P-2)^(P-2)) * (1 + (P-3)^((P+1)/4)) */    \
    /* (0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4*/               \
    /*    897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac << 384) % P */  \
    TO_LIMB_T(0xcd03c9e48671f071), TO_LIMB_T(0x5dab22461fcda5d2),          \
    TO_LIMB_T(0x587042afd3851b95), TO_LIMB_T(0x8eb60ebe01bacb9e),          \
    TO_LIMB_T(0x03f97d6e83d050d2), TO_LIMB_T(0x18f0206554638741)


/* y^2 = x^3 + B,   (4 << 384) % P */
	static const Fp BLST_B_E1_HOST = { BLST_B_E1_LIMBS };
	static const POINTonE1 BLST_G1_HOST{
		Fp{{TO_LIMB_T(0x5cb38790fd530c16), TO_LIMB_T(0x7817fc679976fff5),
						TO_LIMB_T(0x154f95c7143ba1c1), TO_LIMB_T(0xf0ae6acdf3d0e747),
						TO_LIMB_T(0xedce6ecc21dbf440), TO_LIMB_T(0x120177419e0bfb75) }},
		Fp{{ TO_LIMB_T(0xbaac93d50ce72271), TO_LIMB_T(0x8c22631a7918fd8e),
						TO_LIMB_T(0xdd595f13570725ce), TO_LIMB_T(0x51ac582950405194),
						TO_LIMB_T(0x0e1c8c3fad0059c0), TO_LIMB_T(0x0bbc3efc5008a26a) }},
		Fp{{ ONE_MONT_P_LIMBS }}
	};
	static const POINTonE1 BLST_NEG_G1_HOST = { BLST_NEG_G1_LIMBS };
	static const Fp BLST_beta_HOST = { BLST_beta_LIMBS };


	/*
	 * y^2 = x^3 + B, 4 + 4*i
	 */
	const Fp2 BLST_B_E2_HOST = { { BLST_B_E2_LIMBS } };
	const POINTonE2 BLST_G2_HOST{ BLST_G2_LIMBS };
	const POINTonE2 BLST_NEG_G2_HOST = { BLST_NEG_G2_LIMBS };


#if defined(BLST_CUDA) && defined(__CUDACC__)

	static __constant__ Fp BLST_B_E1_DEV = { BLST_B_E1_LIMBS };
	static __constant__ POINTonE1 BLST_G1_DEV{
		{ {TO_LIMB_T(0x5cb38790fd530c16), TO_LIMB_T(0x7817fc679976fff5),
						TO_LIMB_T(0x154f95c7143ba1c1), TO_LIMB_T(0xf0ae6acdf3d0e747),
						TO_LIMB_T(0xedce6ecc21dbf440), TO_LIMB_T(0x120177419e0bfb75) }},
		{ { TO_LIMB_T(0xbaac93d50ce72271), TO_LIMB_T(0x8c22631a7918fd8e),
						TO_LIMB_T(0xdd595f13570725ce), TO_LIMB_T(0x51ac582950405194),
						TO_LIMB_T(0x0e1c8c3fad0059c0), TO_LIMB_T(0x0bbc3efc5008a26a) }},
		{{ ONE_MONT_P_LIMBS }}
	};
	static __constant__ POINTonE1 BLST_NEG_G1_DEV = { BLST_NEG_G1_LIMBS };
	static __constant__ Fp2 BLST_B_E2_DEV = { { BLST_B_E2_LIMBS } };
	static __constant__ POINTonE2 BLST_G2_DEV = { BLST_G2_LIMBS };
	static __constant__ POINTonE2 BLST_NEG_G2_DEV = { BLST_NEG_G2_LIMBS };
	static __constant__ Fp BLST_beta_DEV = { BLST_beta_LIMBS };

#endif // __CUDA_ARCH__


	// return the curve parameter B for the G1 curve: y^2 = x^3 + B
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE const Fp& get_BLST_B_E1(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_B_E1_DEV;
#else
		return BLST_B_E1_HOST;
#endif
	}

	// return the generator point of G1 in affine coordinates, in montgomery form.
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE const Point<Fp>& get_BLST_G1(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_G1_DEV;
#else
		return BLST_G1_HOST;
#endif
	}

	// return the negative of the generator point of G1 in affine coordinates, in montgomery form.
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE const POINTonE1& get_BLST_NEG_G1(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_NEG_G1_DEV;
#else
		return BLST_NEG_G1_HOST;
#endif
	}

	// return the curve parameter B for the G2 curve: y^2 = x^3 + B
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE const Fp2& get_BLST_B_E2(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_B_E2_DEV;
#else
		return BLST_B_E2_HOST;
#endif
	}

	// return the generator point of G2 in affine coordinates, in montgomery form.
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE const POINTonE2& get_BLST_G2(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_G2_DEV;
#else
		return BLST_G2_HOST;
#endif
	}

	// return the negative of the generator point of G2 in affine coordinates, in montgomery form.
	BLST_HD BLST_FORCEINLINE const POINTonE2& get_BLST_NEG_G2(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_NEG_G2_DEV;
#else
		return BLST_NEG_G2_HOST;
#endif
	}


	/* such that beta^3 - 1 = 0  */
	/* -1/2 * (1 + sqrt(-3)) = ((P-2)^(P-2)) * (1 + (P-3)^((P+1)/4)) */
	/* (0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4
		  897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac << 384) % P */
	BLST_HD BLST_FORCEINLINE const Fp& get_BLST_beta(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_beta_DEV;
#else
		return BLST_beta_HOST;
#endif
	}

	template<>
	BLST_HD BLST_FORCEINLINE const Point<Fp>& Point<Fp>::generator()
	{
		return get_BLST_G1();
	}

	template<>
	BLST_HD BLST_FORCEINLINE const Point<Fp2>& Point<Fp2>::generator()
	{
		return get_BLST_G2();
	}

	template<>
	BLST_HD BLST_FORCEINLINE const Point_affine<Fp>& Point_affine<Fp>::generator()
	{
		return (const Point_affine<Fp>&)get_BLST_G1();
	}

	template<>
	BLST_HD BLST_FORCEINLINE const Point_affine<Fp2>& Point_affine<Fp2>::generator()
	{
		return (const Point_affine<Fp2>&)get_BLST_G2();
	}
	//return (const POINTonE1_affine&)get_BLST_G1();


	template<typename Field>
	BLST_HD BLST_INLINE  limb_t Point<Field>::equal(const Point<Field>& p1, const Point<Field>& p2)
	{
		Field Z1Z1, Z2Z2;
		Point_affine<Field> a1, a2;
		bool_t is_inf1 = vec_is_zero(&p1.Z, sizeof(p1.Z));
		bool_t is_inf2 = vec_is_zero(&p2.Z, sizeof(p2.Z));

		Field::sqr(Z1Z1, p1.Z);           /* Z1Z1 = Z1^2 */
		Field::sqr(Z2Z2, p2.Z);           /* Z2Z2 = Z2^2 */

		Field::mul(a1.X, p1.X, Z2Z2);     /* U1 = X1*Z2Z2 */
		Field::mul(a2.X, p2.X, Z1Z1);     /* U2 = X2*Z1Z1 */

		Field::mul(a1.Y, p1.Y, p2.Z);    /* Y1*Z2 */
		Field::mul(a2.Y, p2.Y, p1.Z);    /* Y2*Z1 */

		Field::mul(a1.Y, a1.Y, Z2Z2);      /* S1 = Y1*Z2*Z2Z2 */
		Field::mul(a2.Y, a2.Y, Z1Z1);      /* S2 = Y2*Z1*Z1Z1 */

		return vec_is_equal(&a1, &a2, sizeof(a1)) & (is_inf1 ^ is_inf2 ^ 1);
	}


	template<typename Field>
	BLST_HD BLST_INLINE int Point_affine<Field>::equal(const Point_affine<Field>& a, const Point_affine<Field>& b)
	{
		return (int)vec_is_equal(&a, &b, sizeof(a));

	}


	// conditional swap: if cbit is 1, swap a and b
	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::cswap(Point& BLST_restrict a,
		Point& BLST_restrict b, bool_t cbit)
	{
		BLST_ASSERT(&a != &b);
		BLST_ASSUME(&a != &b);
		vec_cswap(&a, &b, sizeof(Point), cbit);
	}

	template<typename Field>
	BLST_HD BLST_INLINE void  Point<Field>::ccopy(Point& BLST_restrict a,
		const Point& BLST_restrict b, bool_t cbit)
	{
		BLST_ASSERT(&a != &b);
		BLST_ASSUME(&a != &b);
		vec_select(&a, &b, &a, sizeof(Point), cbit);
	}



	template<>
	BLST_HD BLST_FORCEINLINE void Point<Fp>::mul_by_b(Fp& out, const Fp& in)
	{
		Fp::lshift(out, in, 2);
	}

	template<>
	BLST_HD BLST_FORCEINLINE void Point<Fp2>::mul_by_b(Fp2& out, const Fp2& in)
	{
		Fp::sub(out[0], in[0], in[1]);
		Fp::add(out[1], in[0], in[1]);
		Fp::lshift(out[0], out[0], 2);
		Fp::lshift(out[1], out[1], 2);
	}

	template<>
	BLST_HD BLST_FORCEINLINE void Point<Fp>::mul_by_4b(Fp& out, const Fp& in)
	{
		Fp::lshift(out, in, 4);
	}

	template<>
	BLST_HD BLST_INLINE void Point<Fp2>::mul_by_4b(Fp2& out, const Fp2& in)
	{
		Fp::sub(out[0], in[0], in[1]);
		Fp::add(out[1], in[0], in[1]);
		Fp::lshift(out[0], out[0], 4);
		Fp::lshift(out[1], out[1], 4);
	}



	template<typename Field>
	BLST_HD BLST_FORCEINLINE void Point<Field>::cneg(Point& p, bool_t cbit)
	{
		Field::cneg(p.Y, p.Y, cbit);
	}

	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::from_jacobian(Point<Field>& out, const Point<Field>& in)
	{
		Field Z, ZZ;
		limb_t inf = vec_is_zero(in.Z.data, sizeof(in.Z));

		Field::reciprocal(Z, in.Z);                           /* 1/Z */

		Field::sqr(ZZ, Z);
		Field::mul(out.X, in.X, ZZ);                         /* X = X/Z^2 */

		Field::mul(ZZ, ZZ, Z);
		Field::mul(out.Y, in.Y, ZZ);                         /* Y = Y/Z^3 */

		vec_select(&out.Z, &in.Z, &Point::generator().Z,
			sizeof(Point::generator().Z), inf);    /* Z = inf ? 0 : 1 */
	}

	BLST_HD BLST_INLINE void blst_p1_from_affine(POINTonE1& out, const POINTonE1_affine& a)
	{
		vec_copy(&out, &a, sizeof(a));
		vec_select(&out.Z, &a.X, &get_BLST_Rx().p, sizeof(out.Z),
			vec_is_zero(&a, sizeof(a)));
	}
	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::to_affine(Point_affine<Field>& out, const Point& in)
	{
		Point p;
		auto ptr = &in;
		if (!vec_is_equal(ptr->Z.data, &Field::one(), sizeof(in.Z))) {
			Point::from_jacobian(p, *ptr);
			ptr = &p;
		}
		vec_copy(&out, ptr, sizeof(out));
	}



	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::from_affine(Point& out, const Point_affine<Field>& a)
	{
		vec_copy(&out, &a, sizeof(a));
		vec_select(&out.Z, &a.X, &Field::one(), sizeof(out.Z),
			vec_is_zero(&a, sizeof(a)));
	}



	template<>
	BLST_HD BLST_INLINE bool_t POINTonE1_affine::on_curve(const POINTonE1_affine& p)
	{
		Fp XXX, YY;

		Fp::sqr(XXX, p.X);
		Fp::mul(XXX, XXX, p.X);                             /* X^3 */
		Fp::add(XXX, XXX, get_BLST_B_E1());                             /* X^3 + B */

		Fp::sqr(YY, p.Y);                                   /* Y^2 */

		return vec_is_equal(XXX, YY, sizeof(XXX));
	}


	template<>
	BLST_HD BLST_INLINE  bool_t POINTonE2_affine::on_curve(const POINTonE2_affine& p)
	{
		Fp2 XXX, YY;

		Fp2::sqr(XXX, p.X);
		Fp2::mul(XXX, XXX, p.X);                            /* X^3 */
		Fp2::add(XXX, XXX, get_BLST_B_E2());                /* X^3 + B */

		Fp2::sqr(YY, p.Y);                                  /* Y^2 */

		return vec_is_equal(&XXX, &YY, sizeof(XXX)) | is_inf(p);
	}


	template<typename Field>
	BLST_HD BLST_FORCEINLINE int Point<Field>::is_inf(const Point& p)
	{
		return (int)vec_is_zero(&p.Z, sizeof(p.Z));
	}

	template<typename Field>
	BLST_HD BLST_INLINE bool_t Point_affine<Field>::is_inf(const Point_affine& p)
	{
		return vec_is_zero(&p, sizeof(p));
	}


	//BLST_HD BLST_INLINE int blst_p2_affine_on_curve(const POINTonE2_affine& p)
	//{
	//    return (int)(POINTonE2_affine_on_curve(p) | vec_is_zero(&p, sizeof(p)));
	//}
	//BLST_HD BLST_FORCEINLINE int POINTonE1_affine::on_curve(const POINTonE1_affine& p)
	//{
	//    return (int)(POINTonE1_affine_on_curve(p) | vec_is_zero(&p, sizeof(p)));
	//}


	//BLST_HD BLST_INLINE void blst_p2_from_affine(POINTonE2& out, const POINTonE2_affine& a)
	//{
	//    vec_copy(&out, &a, sizeof(a));
	//    vec_select(&out.Z, &a.X, &get_BLST_Rx().p2, sizeof(out.Z),
	//        vec_is_zero(&a, sizeof(a)));
	//}

	template<typename Field>
	BLST_HD BLST_INLINE bool_t Point<Field>::on_curve(const Point& p)
	{
		Field XXX, YY, BZ6;
		limb_t inf = vec_is_zero(&p.Z, sizeof(p.Z));

		Field::sqr(BZ6, p.Z);
		Field::mul(BZ6, BZ6, p.Z);
		Field::sqr(XXX, BZ6);                                 /* Z^6 */
		Point::mul_by_b(BZ6, XXX);                      /* B*Z^6 */

		Field::sqr(XXX, p.X);
		Field::mul(XXX, XXX, p.X);                            /* X^3 */
		Field::add(XXX, XXX, BZ6);                            /* X^3 + B*Z^6 */

		Field::sqr(YY, p.Y);                                  /* Y^2 */

		return vec_is_equal(&XXX, &YY, sizeof(XXX)) | inf;
	}




	//BLST_HD BLST_INLINE limb_t POINTonE1_affine_Serialize_BE(unsigned char out[96],
	//    const POINTonE1_affine& in)
	//{
	//    Fp temp;

	//    Fp::from_mont(temp, in.X);
	//    be_bytes_from_limbs(out, temp, sizeof(temp));

	//    Fp::from_mont(temp, in.Y);
	//    be_bytes_from_limbs(out + 48, temp, sizeof(temp));

	//    return Fp::sgn0_pty(temp);
	//}



	//BLST_HD BLST_INLINE limb_t POINTonE2_affine_Serialize_BE(unsigned char out[192],
	//    const POINTonE2_affine& in)
	//{
	//    Fp2 temp;

	//    Fp::from_mont(temp[1], in.X[1]);
	//    be_bytes_from_limbs(out, temp[1], sizeof(temp[1]));
	//    Fp::from_mont(temp[0], in.X[0]);
	//    be_bytes_from_limbs(out + 48, temp[0], sizeof(temp[0]));

	//    Fp::from_mont(temp[1], in.Y[1]);
	//    be_bytes_from_limbs(out + 96, temp[1], sizeof(temp[1]));
	//    Fp::from_mont(temp[0], in.Y[0]);
	//    be_bytes_from_limbs(out + 144, temp[0], sizeof(temp[0]));

	//    return Fp2::sgn0_pty(temp);
	//}

	template<typename Field>
	BLST_HD BLST_INLINE limb_t Point_affine<Field>::serialize_BE(unsigned char* out,
		const Point_affine& in)
	{
		if constexpr (std::is_same_v<Field, Fp>) {
			Fp temp;

			Fp::from_mont(temp, in.X);
			be_bytes_from_limbs(out, temp, sizeof(temp));

			Fp::from_mont(temp, in.Y);
			be_bytes_from_limbs(out + 48, temp, sizeof(temp));

			return Fp::sgn0_pty(temp);
		}
		else {
			Fp2 temp;

			Fp::from_mont(temp[1], in.X[1]);
			be_bytes_from_limbs(out, temp[1], sizeof(temp[1]));
			Fp::from_mont(temp[0], in.X[0]);
			be_bytes_from_limbs(out + 48, temp[0], sizeof(temp[0]));

			Fp::from_mont(temp[1], in.Y[1]);
			be_bytes_from_limbs(out + 96, temp[1], sizeof(temp[1]));
			Fp::from_mont(temp[0], in.Y[0]);
			be_bytes_from_limbs(out + 144, temp[0], sizeof(temp[0]));

			return Fp2::sgn0_pty(temp);
		}
	}

	template<typename Field>
	BLST_HD BLST_INLINE void Point_affine<Field>::serialize(unsigned char* out,
		const Point_affine& in)
	{
		if (vec_is_zero(in.X, sizeof(in))) {
			bytes_zero(out, sizeof(in));
			out[0] = 0x40;    /* infinity bit */
		}
		else {
			(void)serialize_BE(out, in);
		}
	}



	template<typename Field>
	BLST_HD BLST_INLINE limb_t Point<Field>::serialize_BE(unsigned char* out,
		const Point& in)
	{
		Point p;
		auto ptr = &in;
		if (!vec_is_equal(&ptr->Z, &Field::one(), sizeof(ptr->Z))) {
			Point::from_jacobian(p, *ptr);
			ptr = &p;
		}

		return Point_affine<Field>::serialize_BE(out, (const Point_affine<Field>&) * ptr);
	}

	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::serialize(unsigned char* out, const Point& in)
	{
		if (vec_is_zero(&in.Z, sizeof(in.Z))) {
			bytes_zero(out, sizeof(Field) * 2);
			out[0] = 0x40;    /* infinity bit */
		}
		else {
			(void)Point::serialize_BE(out, in);
		}
	}



	template<typename Field>
	BLST_HD BLST_INLINE limb_t Point_affine<Field>::compress_BE(unsigned char* out, const Point_affine& in)
	{
		if constexpr (std::is_same_v<Field, Fp2>) {
			Fp temp;

			Fp::from_mont(temp, in.X[1]);
			be_bytes_from_limbs(out, temp, sizeof(temp));
			Fp::from_mont(temp, in.X[0]);
			be_bytes_from_limbs(out + 48, temp, sizeof(temp));

			return Fp2::sgn0_pty_mont(in.Y);

		}
		else
		{
			Fp temp;

			Fp::from_mont(temp, in.X);
			be_bytes_from_limbs(out, temp, sizeof(temp));

			return Fp::sgn0_pty_mont(in.Y);
		}
	}

	template<typename Field>
	BLST_HD BLST_INLINE void Point_affine<Field>::compress(unsigned char* out, const Point_affine& in)
	{
		if (is_inf(in)) {
			bytes_zero(out, sizeof(Field));
			out[0] = 0xc0;    /* compressed and infinity bits */
		}
		else {
			limb_t sign = Point_affine::compress_BE(out, in);
			out[0] |= (unsigned char)(0x80 | ((sign & 2) << 4));
		}
	}



	template<typename Field>
	BLST_HD BLST_INLINE limb_t Point<Field>::compress_BE(unsigned char* out,
		const Point& in)
	{
		Point p;
		auto ptr = &in;

		if (!vec_is_equal(&ptr->Z, &Field::one(), sizeof(in.Z))) {
			Point::from_jacobian(p, *ptr);
			ptr = &p;
		}

		return Point_affine<Field>::compress_BE(out, (const Point_affine<Field>&) * ptr);
	}

	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::compress(unsigned char* out, const Point& in)
	{
		if (vec_is_zero(&in.Z, sizeof(in.Z))) {
			bytes_zero(out, sizeof(Field));
			out[0] = 0xc0;    /* compressed and infinity bits */
		}
		else {
			limb_t sign = Point::compress_BE(out, in);
			out[0] |= (unsigned char)(0x80 | ((sign & 2) << 4));
		}
	}


	template<typename Field>
	BLST_HD BLST_INLINE limb_t Point_affine<Field>::uncompress_BE(Point_affine& out,
		const unsigned char* in)
	{
		if constexpr (std::is_same_v<Field, Fp2>) {

			POINTonE2_affine ret;
			Fp temp;

			limbs_from_be_bytes(ret.X[1], in, sizeof(ret.X[1]));
			limbs_from_be_bytes(ret.X[0], in + 48, sizeof(ret.X[0]));

			/* clear top 3 bits in case caller was conveying some information there */
			ret.X[1][sizeof(ret.X[1]) / sizeof(limb_t) - 1] &= ((limb_t)0 - 1) >> 3;
			Fp::add(temp, ret.X[1], Fp::zero());  /* less than modulus? */
			if (!vec_is_equal(temp, ret.X[1], sizeof(temp)))
				return (limb_t)0 - BLST_BAD_ENCODING;

			Fp::add(temp, ret.X[0], Fp::zero());  /* less than modulus? */
			if (!vec_is_equal(temp, ret.X[0], sizeof(temp)))
				return (limb_t)0 - BLST_BAD_ENCODING;

			Fp::mul(ret.X[0], ret.X[0], get_BLST_RR());
			Fp::mul(ret.X[1], ret.X[1], get_BLST_RR());

			Fp2::sqr(ret.Y, ret.X);
			Fp2::mul(ret.Y, ret.Y, ret.X);
			Fp2::add(ret.Y, ret.Y, get_BLST_B_E2());                        /* X^3 + B */
			if (!Fp2::sqrt(ret.Y, ret.Y))
				return (limb_t)0 - BLST_POINT_NOT_ON_CURVE;

			vec_copy(&out, &ret, sizeof(ret));

			return Fp2::sgn0_pty_mont(out.Y);
		}
		else {
			Point_affine ret;
			Fp temp;

			limbs_from_be_bytes(ret.X, in, sizeof(ret.X));
			/* clear top 3 bits in case caller was conveying some information there */
			ret.X[sizeof(ret.X) / sizeof(limb_t) - 1] &= ((limb_t)0 - 1) >> 3;
			Fp::add(temp, ret.X, Fp::zero());  /* less than modulus? */
			if (!vec_is_equal(temp, ret.X, sizeof(temp)))
				return (limb_t)0 - BLST_BAD_ENCODING;
			Fp::mul(ret.X, ret.X, get_BLST_RR());

			Fp::sqr(ret.Y, ret.X);
			Fp::mul(ret.Y, ret.Y, ret.X);
			Fp::add(ret.Y, ret.Y, get_BLST_B_E1());                         /* X^3 + B */
			if (!Fp::sqrt(ret.Y, ret.Y))
				return (limb_t)0 - BLST_POINT_NOT_ON_CURVE;

			vec_copy(&out, &ret, sizeof(ret));

			return Fp::sgn0_pty_mont(out.Y);
		}
	}

	template<typename Field>
	BLST_HD BLST_INLINE BLST_ERROR Point_affine<Field>::uncompress(Point_affine& out,
		const unsigned char* in)
	{
		unsigned char in0 = in[0];
		limb_t sgn0_pty;

		if ((in0 & 0x80) == 0)      /* compressed bit */
			return BLST_BAD_ENCODING;

		if (in0 & 0x40) {           /* infinity bit */
			if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in + 1, sizeof(Field) - 1)) {
				vec_zero(&out, sizeof(out));
				return BLST_SUCCESS;
			}
			else {
				return BLST_BAD_ENCODING;
			}
		}

		sgn0_pty = Point_affine::uncompress_BE(out, in);

		if (sgn0_pty > 3)
			return (BLST_ERROR)(0 - sgn0_pty); /* POINT_NOT_ON_CURVE */

		sgn0_pty >>= 1; /* skip over parity bit */
		sgn0_pty ^= (in0 & 0x20) >> 5;
		Field::cneg(out.Y, out.Y, sgn0_pty);

		/* (0,±2) is not in group, but application might want to ignore? */
		return vec_is_zero(&out.X, sizeof(out.X)) ? BLST_POINT_NOT_IN_GROUP
			: BLST_SUCCESS;
	}


	//BLST_HD BLST_INLINE BLST_ERROR POINTonE2_Uncompress_Z(POINTonE2_affine& out,
	//    const unsigned char in[96])
	//{
	//    unsigned char in0 = in[0];
	//    limb_t sgn0_pty;

	//    if ((in0 & 0x80) == 0)      /* compressed bit */
	//        return BLST_BAD_ENCODING;

	//    if (in0 & 0x40) {           /* infinity bit */
	//        if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in + 1, 95)) {
	//            vec_zero(&out, sizeof(out));
	//            return BLST_SUCCESS;
	//        }
	//        else {
	//            return BLST_BAD_ENCODING;
	//        }
	//    }

	//    sgn0_pty = POINTonE2_Uncompress_BE(out, in);

	//    if (sgn0_pty > 3)
	//        return (BLST_ERROR)(0 - sgn0_pty); /* POINT_NOT_ON_CURVE */

	//    sgn0_pty >>= 1; /* skip over parity bit */
	//    sgn0_pty ^= (in0 & 0x20) >> 5;
	//    Fp2::cneg(out.Y, out.Y, sgn0_pty);

	//    return BLST_SUCCESS;
	//}

	template<typename Field>
	BLST_HD BLST_INLINE BLST_ERROR Point_affine<Field>::deserialize_BE(
		Point_affine& out,
		const unsigned char* in)
	{
		if constexpr (std::is_same_v<Field, Fp2>) {

			POINTonE2_affine ret;
			Fp temp;

			limbs_from_be_bytes(ret.X[1], in, sizeof(ret.X[1]));
			limbs_from_be_bytes(ret.X[0], in + 48, sizeof(ret.X[0]));
			limbs_from_be_bytes(ret.Y[1], in + 96, sizeof(ret.Y[1]));
			limbs_from_be_bytes(ret.Y[0], in + 144, sizeof(ret.Y[0]));

			/* clear top 3 bits in case caller was conveying some information there */
			ret.X[1][sizeof(ret.X[1]) / sizeof(limb_t) - 1] &= ((limb_t)0 - 1) >> 3;
			Fp::add(temp, ret.X[1], Fp::zero());  /* less than modulus? */
			if (!vec_is_equal(temp, ret.X[1], sizeof(temp)))
				return BLST_BAD_ENCODING;

			Fp::add(temp, ret.X[0], Fp::zero());  /* less than modulus? */
			if (!vec_is_equal(temp, ret.X[0], sizeof(temp)))
				return BLST_BAD_ENCODING;

			Fp::add(temp, ret.Y[1], Fp::zero());  /* less than modulus? */
			if (!vec_is_equal(temp, ret.Y[1], sizeof(temp)))
				return BLST_BAD_ENCODING;

			Fp::add(temp, ret.Y[0], Fp::zero());  /* less than modulus? */
			if (!vec_is_equal(temp, ret.Y[0], sizeof(temp)))
				return BLST_BAD_ENCODING;

			Fp::mul(ret.X[0], ret.X[0], get_BLST_RR());
			Fp::mul(ret.X[1], ret.X[1], get_BLST_RR());
			Fp::mul(ret.Y[0], ret.Y[0], get_BLST_RR());
			Fp::mul(ret.Y[1], ret.Y[1], get_BLST_RR());

			if (!POINTonE2_affine::on_curve(ret))
				return BLST_POINT_NOT_ON_CURVE;

			vec_copy(&out, &ret, sizeof(ret));

			return BLST_SUCCESS;
		}
		else
		{

			POINTonE1_affine ret;
			Fp temp;

			limbs_from_be_bytes(ret.X, in, sizeof(ret.X));
			limbs_from_be_bytes(ret.Y, in + 48, sizeof(ret.Y));

			/* clear top 3 bits in case caller was conveying some information there */
			ret.X[sizeof(ret.X) / sizeof(limb_t) - 1] &= ((limb_t)0 - 1) >> 3;
			Fp::add(temp, ret.X, Fp::zero());  /* less than modulus? */
			if (!vec_is_equal(temp, ret.X, sizeof(temp)))
				return BLST_BAD_ENCODING;

			Fp::add(temp, ret.Y, Fp::zero());  /* less than modulus? */
			if (!vec_is_equal(temp, ret.Y, sizeof(temp)))
				return BLST_BAD_ENCODING;

			Fp::mul(ret.X, ret.X, get_BLST_RR());
			Fp::mul(ret.Y, ret.Y, get_BLST_RR());

			if (!POINTonE1_affine::on_curve(ret))
				return BLST_POINT_NOT_ON_CURVE;

			vec_copy(&out, &ret, sizeof(ret));

			/* (0,±2) is not in group, but application might want to ignore? */
			return vec_is_zero(out.X, sizeof(out.X)) ? BLST_POINT_NOT_IN_GROUP
				: BLST_SUCCESS;
		}
	}

	template<typename Field>
	BLST_HD BLST_INLINE BLST_ERROR Point_affine<Field>::deserialize(Point_affine& out,
		const unsigned char* in)
	{
		unsigned char in0 = in[0];

		if ((in0 & 0xe0) == 0)
			return Point_affine::deserialize_BE(out, in);

		if (in0 & 0x80)             /* compressed bit */
			return Point_affine::uncompress(out, in);

		if (in0 & 0x40) {           /* infinity bit */
			if (byte_is_zero(in0 & 0x3f) & bytes_are_zero(in + 1, sizeof(Point_affine) - 1)) {
				vec_zero(&out, sizeof(out));
				return BLST_SUCCESS;
			}
		}

		return BLST_BAD_ENCODING;
	}




	/* Gallant-Lambert-Vanstone, ~45% faster than mult_w<POINTonE1,5> */
	template<>
	BLST_HD BLST_INLINE void POINTonE1::mult_glsv(POINTonE1& out, const POINTonE1& in,
		const pow256 SK)
	{
		//union { Fr l; pow256 s; } val;
		Fr l;
		pow256 s;
		/* SK/z^2 [in constant time] */

		limbs_from_le_bytes(l, SK, 32);
		div_by_zz(l);
		le_bytes_from_limbs(s, l, 32);

		{
			const byte* scalars[2] = { s + 16, s };
			POINTonE1 table[2][1 << (5 - 1)];   /* 4.5KB */
			size_t i;

			precompute_w<POINTonE1, 5>(table[0], in);
			for (i = 0; i < 1 << (5 - 1); i++) {
				Fp::mul(table[1][i].X, table[0][i].X, get_BLST_beta());
				Fp::cneg(table[1][i].Y, table[0][i].Y, 1);
				vec_copy(table[1][i].Z, table[0][i].Z, sizeof(table[1][i].Z));
			}

			s_mult_w<POINTonE1, 5>(out, NULL, 2, scalars, 128, table);
			POINTonE1::cneg(out, 1);
			Fp::mul(out.Z, out.Z, get_BLST_beta());
			Fp::mul(out.Z, out.Z, get_BLST_beta());
		}

		vec_zero(l, sizeof(l));   /* scrub the copy of SK */
		vec_zero(s, sizeof(s));   /* scrub the copy of SK */
	}


	BLST_HD BLST_INLINE void psi(POINTonE2& out, const POINTonE2& in)
	{
		const Fp2 frobenius_x = { /* 1/(1 + i)^((P-1)/3) */
		  Fp{ 0 },
		  Fp{ /* (0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4
				  897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad << 384) % P */
			TO_LIMB_T(0x890dc9e4867545c3), TO_LIMB_T(0x2af322533285a5d5),
			TO_LIMB_T(0x50880866309b7e2c), TO_LIMB_T(0xa20d1b8c7e881024),
			TO_LIMB_T(0x14e4f04fe2db9068), TO_LIMB_T(0x14e56d3f1564853a) }
		};
		const Fp2 frobenius_y = { /* 1/(1 + i)^((P-1)/2) */
		  Fp{ /* (0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60
				  ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2 << 384) % P */
			TO_LIMB_T(0x3e2f585da55c9ad1), TO_LIMB_T(0x4294213d86c18183),
			TO_LIMB_T(0x382844c88b623732), TO_LIMB_T(0x92ad2afd19103e18),
			TO_LIMB_T(0x1d794e4fac7cf0b9), TO_LIMB_T(0x0bd592fc7d825ec8) },
		  Fp{ /* (0x06af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e
				  77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09 << 384) % P */
			TO_LIMB_T(0x7bcfa7a25aa30fda), TO_LIMB_T(0xdc17dec12a927e7c),
			TO_LIMB_T(0x2f088dd86b4ebef1), TO_LIMB_T(0xd1ca2087da74d4a7),
			TO_LIMB_T(0x2da2596696cebc1d), TO_LIMB_T(0x0e2b7eedbbfd87d2) },
		};

		vec_copy(&out, &in, sizeof(out));
		Fp::cneg(out.X[1], out.X[1], 1);   Fp2::mul(out.X, out.X, frobenius_x);
		Fp::cneg(out.Y[1], out.Y[1], 1);   Fp2::mul(out.Y, out.Y, frobenius_y);
		Fp::cneg(out.Z[1], out.Z[1], 1);
	}

	/* Galbraith-Lin-Scott, ~67% faster than POINTonE2_mul_w5 */
	template<>
	BLST_HD BLST_INLINE void POINTonE2::mult_glsv(POINTonE2& out, const POINTonE2& in,
		const pow256 SK)
	{
		Fr l;
		pow256 s;

		/* break down SK to "digits" with |z| as radix [in constant time] */
		limbs_from_le_bytes(l, SK, 32);
		div_by_zz(l);
		div_by_z(l);
		div_by_z(l + NLIMBS(256) / 2);


		le_bytes_from_limbs(s, l, 32);

		{
			const byte* scalars[2] = { s, NULL };
			POINTonE2 table[4][1 << (5 - 1)];   /* 18KB */
			size_t i;

			//POINTonE2_precompute_w5(table[0], in);
			precompute_w<POINTonE2, 5>(table[0], in);

			for (i = 0; i < 1 << (5 - 1); i++) {
				psi(table[1][i], table[0][i]);
				psi(table[2][i], table[1][i]);
				psi(table[3][i], table[2][i]);
				POINTonE2::cneg(table[1][i], 1); /* account for z being negative */
				POINTonE2::cneg(table[3][i], 1);
			}

			s_mult_w<POINTonE2, 5>(out, NULL, 4, scalars, 64, table);
		}

		// copy back the scalar to scrub the stack
		for (u64 i = 0; i < NLIMBS(256); i++) {
			l[i] = s[i];
		}

		vec_zero(l, sizeof(l));   /* scrub the copy of SK */
	}

	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::sign(Point& out, const Point& in, const pow256 SK)
	{
		Field Z, ZZ;
		limb_t inf;

		mult_glsv(out, in, SK);

		/* convert to affine to remove possible bias in out.Z */
		inf = vec_is_zero(&out.Z, sizeof(out.Z));
#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
		Field::flt_reciprocal(Z, out.Z);                       /* 1/Z   */
#else
		Field::reciprocal(Z, out.Z);                           /* 1/Z   */
#endif

		Field::sqr(ZZ, Z);
		Field::mul(out.X, out.X, ZZ);                         /* X = X/Z^2 */

		Field::mul(ZZ, ZZ, Z);
		Field::mul(out.Y, out.Y, ZZ);                         /* Y = Y/Z^3 */

		vec_select(&out.Z, &out.Z, &Point<Field>::generator().Z, sizeof(Point<Field>::generator().Z),
			inf);                            /* Z = inf ? 0 : 1 */
	}


	template<typename Field>
	BLST_HD void Point<Field>::sign2(
		unsigned char* out /*optional output*/,
		Point_affine<Field>* sig /*optional output*/,
		const Point& hash,
		const pow256 SK)
	{
		Point P;

		Point::sign(P, hash, SK);
		if (sig != NULL)
			vec_copy(sig, &P, sizeof(sig));
		if (out != NULL) {
			limb_t sgn0_pty = Point::serialize_BE(out, P);
			out[0] |= (sgn0_pty & 2) << 4;      /* pre-decorate */
			out[0] |= vec_is_zero(&P.Z, sizeof(P.Z)) << 6;
		}
	}

	//
	//
	//	template<>
	//	BLST_HD BLST_INLINE void POINTonE2::sign(POINTonE2& out, const POINTonE2& in, const pow256 SK)
	//	{
	//		Fp2 Z, ZZ;
	//		limb_t inf;
	//
	//		mult_glsv(out, in, SK);
	//
	//		/* convert to affine to remove possible bias in out.Z */
	//		inf = vec_is_zero(&out.Z, sizeof(out.Z));
	//#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
	//		Fp2::flt_reciprocal(Z, out.Z);                      /* 1/Z   */
	//#else
	//		Fp2::reciprocal(Z, out.Z);                          /* 1/Z   */
	//#endif
	//
	//		Fp2::sqr(ZZ, Z);
	//		Fp2::mul(out.X, out.X, ZZ);                        /* X = X/Z^2 */
	//
	//		Fp2::mul(ZZ, ZZ, Z);
	//		Fp2::mul(out.Y, out.Y, ZZ);                        /* Y = Y/Z^3 */
	//
	//		vec_select(&out.Z, &out.Z, &get_BLST_G2().Z, sizeof(get_BLST_G2().Z),
	//			inf);                            /* Z = inf ? 0 : 1 */
	//	}

	template<typename Field>
	BLST_HD BLST_FORCEINLINE void Point<Field>::to_pk(Point& out, const pow256 SK)
	{
		sign(out, Point::generator(), SK);
	}


	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::to_pk2(
		unsigned char* out,
		Point_affine<Field>* PK,
		const pow256 SK)
	{
		Point P;

		Point::sign(P, Point::generator(), SK);
		if (PK != NULL)
			vec_copy(PK, &P, sizeof(*PK));
		if (out != NULL) {
			limb_t sgn0_pty = Point::serialize_BE(out, P);
			out[0] |= (sgn0_pty & 2) << 4;      /* pre-decorate */
			out[0] |= vec_is_zero(&P.Z, sizeof(P.Z)) << 6;
		}
	}



	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::mult(
		Point& out,
		const Point& a,
		const byte* scalar,
		size_t nbits)
	{
		if (nbits < 176) {
			if (nbits)
				mult_w<Point, 4>(out, a, scalar, nbits);
			else
				vec_zero(&out, sizeof(out));
		}
		else if (nbits <= 256) {
			union { Fr l; pow256 s; } val;
			size_t i, j, top, mask = (size_t)0 - 1;

			/* this is not about constant-time-ness, but branch optimization */
			for (top = (nbits + 7) / 8, i = 0, j = 0; i < sizeof(val.s);) {
				val.s[i++] = scalar[j] & mask;
				mask = 0 - ((i - top) >> (8 * sizeof(top) - 1));
				j += 1 & mask;
			}

			if (Fr::check_mod(val.l))  /* z^4 is the formal limit */
				Point::mult_glsv(out, a, val.s);
			else    /* should never be the case, added for formal completeness */
				mult_w<Point, 5>(out, a, scalar, nbits);

			vec_zero(val.l, sizeof(val));
		}
		else {    /* should never be the case, added for formal completeness */
			mult_w<Point, 5>(out, a, scalar, nbits);
		}
	}



	/*
	 * Addition that can handle doubling [as well as points at infinity,
	 * which are encoded as Z==0] in constant time. It naturally comes at
	 * cost, but this subroutine should be called only when independent
	 * points are processed, which is considered reasonable compromise.
	 * For example, ptype##s_mult_w5 calls it, but since *major* gain is
	 * result of pure doublings being effectively divided by amount of
	 * points, slightly slower addition can be tolerated. But what is the
	 * additional cost more specifically? Best addition result is 11M+5S,
	 * while this routine takes 13M+5S (+1M+1S if a4!=0), as per
	 *
	 * -------------+-------------
	 * addition     | doubling
	 * -------------+-------------
	 * U1 = X1*Z2^2 | U1 = X1
	 * U2 = X2*Z1^2 |
	 * S1 = Y1*Z2^3 | S1 = Y1
	 * S2 = Y2*Z1^3 |
	 * zz = Z1*Z2   | zz = Z1
	 * H = U2-U1    | H' = 2*Y1
	 * R = S2-S1    | R' = 3*X1^2[+a*Z1^4]
	 * sx = U1+U2   | sx = X1+X1
	 * -------------+-------------
	 * H!=0 || R!=0 | H==0 && R==0
	 *
	 *      X3 = R^2-H^2*sx
	 *      Y3 = R*(H^2*U1-X3)-H^3*S1
	 *      Z3 = H*zz
	 *
	 * As for R!=0 condition in context of H==0, a.k.a. P-P. The result is
	 * infinity by virtue of Z3 = (U2-U1)*zz = H*zz = 0*zz == 0.
	 */
	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::dadd(Point& out, const Point& p1, const Point& p2, const Field* a4)
	{
		Point p3; /* starts as (U1, S1, zz) from addition side */
		struct { Field H, R, sx; } add, dbl;
		bool_t p1inf, p2inf, is_dbl;

		Field::add(dbl.sx, p1.X, p1.X);  /* sx = X1+X1 */
		Field::sqr(dbl.R, p1.X);          /* X1^2 */
		Field::mul_by_3(dbl.R, dbl.R);     /* R = 3*X1^2 */
		Field::add(dbl.H, p1.Y, p1.Y);   /* H = 2*Y1 */

		p2inf = vec_is_zero(&p2.Z, sizeof(p2.Z));
		Field::sqr(p3.X, p2.Z);           /* Z2^2 */
		Field::mul(p3.Z, p1.Z, p2.Z);    /* Z1*Z2 */
		p1inf = vec_is_zero(&p1.Z, sizeof(p1.Z));
		Field::sqr(add.H, p1.Z);          /* Z1^2 */

		if (a4 != NULL) {
			Field::sqr(p3.Y, add.H);       /* Z1^4, [borrow p3.Y] */
			Field::mul(p3.Y, p3.Y, *a4);
			Field::add(dbl.R, dbl.R, p3.Y);/* R = 3*X1^2+a*Z1^4 */
		}

		Field::mul(p3.Y, p1.Y, p2.Z);
		Field::mul(p3.Y, p3.Y, p3.X);      /* S1 = Y1*Z2^3 */
		Field::mul(add.R, p2.Y, p1.Z);
		Field::mul(add.R, add.R, add.H);   /* S2 = Y2*Z1^3 */
		Field::sub(add.R, add.R, p3.Y);    /* R = S2-S1 */

		Field::mul(p3.X, p3.X, p1.X);     /* U1 = X1*Z2^2 */
		Field::mul(add.H, add.H, p2.X);   /* U2 = X2*Z1^2 */

		Field::add(add.sx, add.H, p3.X);   /* sx = U1+U2 */
		Field::sub(add.H, add.H, p3.X);    /* H = U2-U1 */

		/* make the choice between addition and doubling */
		is_dbl = vec_is_zero(&add.H, 2 * sizeof(add.H));
		vec_select(&p3, &p1, &p3, sizeof(p3), is_dbl);
		vec_select(&add, &dbl, &add, sizeof(add), is_dbl);
		/* |p3| and |add| hold all inputs now, |p3| will hold output */

		Field::mul(p3.Z, p3.Z, add.H);     /* Z3 = H*Z1*Z2 */

		Field::sqr(dbl.H, add.H);          /* H^2 */
		Field::mul(dbl.R, dbl.H, add.H);   /* H^3 */
		Field::mul(dbl.R, dbl.R, p3.Y);    /* H^3*S1 */
		Field::mul(p3.Y, dbl.H, p3.X);     /* H^2*U1 */

		Field::mul(dbl.H, dbl.H, add.sx);  /* H^2*sx */
		Field::sqr(p3.X, add.R);           /* R^2 */
		Field::sub(p3.X, p3.X, dbl.H);     /* X3 = R^2-H^2*sx */

		Field::sub(p3.Y, p3.Y, p3.X);      /* H^2*U1-X3 */
		Field::mul(p3.Y, p3.Y, add.R);     /* R*(H^2*U1-X3) */
		Field::sub(p3.Y, p3.Y, dbl.R);     /* Y3 = R*(H^2*U1-X3)-H^3*S1 */

		vec_select(&p3, &p1, &p3, sizeof(Point), p2inf);
		vec_select(&out, &p2, &p3, sizeof(Point), p1inf);
	}




	/*
	 * Addition with affine point that can handle doubling [as well as
	 * points at infinity, with |p1| being encoded as Z==0 and |p2| as
	 * X,Y==0] in constant time. But at what additional cost? Best
	 * addition result is 7M+4S, while this routine takes 8M+5S, as per
	 *
	 * -------------+-------------
	 * addition     | doubling
	 * -------------+-------------
	 * U1 = X1      | U1 = X2
	 * U2 = X2*Z1^2 |
	 * S1 = Y1      | S1 = Y2
	 * S2 = Y2*Z1^3 |
	 * H = U2-X1    | H' = 2*Y2
	 * R = S2-Y1    | R' = 3*X2^2[+a]
	 * sx = X1+U2   | sx = X2+X2
	 * zz = H*Z1    | zz = H'
	 * -------------+-------------
	 * H!=0 || R!=0 | H==0 && R==0
	 *
	 *      X3 = R^2-H^2*sx
	 *      Y3 = R*(H^2*U1-X3)-H^3*S1
	 *      Z3 = zz
	 *
	 * As for R!=0 condition in context of H==0, a.k.a. P-P. The result is
	 * infinity by virtue of Z3 = (U2-U1)*zz = H*zz = 0*zz == 0.
	 */
	 //#define POINT_DADD_AFFINE_IMPL_A0(ptype, bits, field, one) 
	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::dadd_affine(Point& out, const Point& p1,
		const Point_affine<Field>& p2)
	{
		Point p3; /* starts as (,, H*Z1) from addition side */
		struct { Field H, R, sx; } add, dbl;
		bool_t p1inf, p2inf, is_dbl;

		p2inf = vec_is_zero(&p2.X, 2 * sizeof(p2.X));
		Field::add(dbl.sx, p2.X, p2.X);  /* sx = X2+X2 */
		Field::sqr(dbl.R, p2.X);          /* X2^2 */
		Field::mul_by_3(dbl.R, dbl.R);     /* R = 3*X2^2 */
		Field::add(dbl.H, p2.Y, p2.Y);   /* H = 2*Y2 */

		p1inf = vec_is_zero(&p1.Z, sizeof(p1.Z));
		Field::sqr(add.H, p1.Z);          /* Z1^2 */
		Field::mul(add.R, add.H, p1.Z);   /* Z1^3 */
		Field::mul(add.R, add.R, p2.Y);   /* S2 = Y2*Z1^3 */
		Field::sub(add.R, add.R, p1.Y);   /* R = S2-Y1 */

		Field::mul(add.H, add.H, p2.X);   /* U2 = X2*Z1^2 */

		Field::add(add.sx, add.H, p1.X);  /* sx = X1+U2 */
		Field::sub(add.H, add.H, p1.X);   /* H = U2-X1 */

		Field::mul(p3.Z, add.H, p1.Z);    /* Z3 = H*Z1 */

		/* make the choice between addition and doubling */
		is_dbl = vec_is_zero(&add.H, 2 * sizeof(add.H));
		vec_select(&p3.X, &p2, &p1, 2 * sizeof(p3.X), is_dbl);
		vec_select(&p3.Z, &dbl.H, &p3.Z, sizeof(p3.Z), is_dbl);
		vec_select(&add, &dbl, &add, sizeof(add), is_dbl);
		/* |p3| and |add| hold all inputs now, |p3| will hold output */

		Field::sqr(dbl.H, add.H);          /* H^2 */
		Field::mul(dbl.R, dbl.H, add.H);   /* H^3 */
		Field::mul(dbl.R, dbl.R, p3.Y);    /* H^3*S1 */
		Field::mul(p3.Y, dbl.H, p3.X);     /* H^2*U1 */

		Field::mul(dbl.H, dbl.H, add.sx);  /* H^2*sx */
		Field::sqr(p3.X, add.R);           /* R^2 */
		Field::sub(p3.X, p3.X, dbl.H);     /* X3 = R^2-H^2*sx */

		Field::sub(p3.Y, p3.Y, p3.X);      /* H^2*U1-X3 */
		Field::mul(p3.Y, p3.Y, add.R);     /* R*(H^2*U1-X3) */
		Field::sub(p3.Y, p3.Y, dbl.R);     /* Y3 = R*(H^2*U1-X3)-H^3*S1 */

		vec_select(&p3.X, &p2, &p3.X, 2 * sizeof(p3.X), p1inf);
		vec_select(&p3.Z, &Field::one(), &p3.Z, sizeof(p3.Z), p1inf);
		vec_select(&out, &p1, &p3, sizeof(p3), p2inf);
	}



	/*
	 * https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
	 * with twist to handle either input at infinity, which are encoded as Z==0.
	 */
	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::add(Point& out, const Point& p1, const Point& p2)
	{

		//#if (defined(DEBUG) || !defined(NDEBUG)) && !defined(__CUDA_ARCH__)
		//		// Check if we're trying to add a point to itself in debug builds
		//		if (!vec_is_zero(a.Z, sizeof(a.Z)) && !vec_is_zero(b.Z, sizeof(b.Z))) { // If neither point is infinity
		//			// Convert both points to affine to check if they're the same
		//			Fp Z_a, ZZ_a, Z_b, ZZ_b;
		//
		//			// Convert a to affine
		//			Fp::reciprocal(Z_a, a.Z);            // 1/Z_a
		//			Fp::sqr(ZZ_a, Z_a);                   // 1/Z_a^2
		//			Fp a_X_affine, a_Y_affine;
		//			Fp::mul(a_X_affine, a.X, ZZ_a);      // X_a/Z_a^2
		//
		//			// Convert b to affine
		//			Fp::reciprocal(Z_b, b.Z);            // 1/Z_b
		//			Fp::sqr(ZZ_b, Z_b);                   // 1/Z_b^2
		//			Fp b_X_affine, b_Y_affine;
		//			Fp::mul(b_X_affine, b.X, ZZ_b);      // X_b/Z_b^2
		//
		//			// Compare X coordinates
		//			if (vec_is_equal(a_X_affine, b_X_affine, sizeof(Fp))) {
		//				// Also check Y coordinates
		//				Fp::mul(ZZ_a, ZZ_a, Z_a);         // 1/Z_a^3
		//				Fp::mul(a_Y_affine, a.Y, ZZ_a);  // Y_a/Z_a^3
		//
		//				Fp::mul(ZZ_b, ZZ_b, Z_b);         // 1/Z_b^3
		//				Fp::mul(b_Y_affine, b.Y, ZZ_b);  // Y_b/Z_b^3
		//
		//				if (vec_is_equal(a_Y_affine, b_Y_affine, sizeof(Fp))) {
		//					printf("ERROR: blst_p1_add called with the same point. "
		//						"Use blst_p1_add_or_double instead.\n");
		//					exit(1);
		//				}
		//			}
		//		}
		//#endif

		Point p3;
		Field Z1Z1, Z2Z2, U1, S1, H, I, J;
		bool_t p1inf, p2inf;

		p1inf = vec_is_zero(&p1.Z, sizeof(p1.Z));
		Field::sqr(Z1Z1, p1.Z);           /* Z1Z1 = Z1^2 */

		Field::mul(p3.Z, Z1Z1, p1.Z);     /* Z1*Z1Z1 */
		Field::mul(p3.Z, p3.Z, p2.Y);     /* S2 = Y2*Z1*Z1Z1 */

		p2inf = vec_is_zero(&p2.Z, sizeof(p2.Z));
		Field::sqr(Z2Z2, p2.Z);           /* Z2Z2 = Z2^2 */

		Field::mul(S1, Z2Z2, p2.Z);       /* Z2*Z2Z2 */
		Field::mul(S1, S1, p1.Y);         /* S1 = Y1*Z2*Z2Z2 */

		Field::sub(p3.Z, p3.Z, S1);        /* S2-S1 */
		Field::add(p3.Z, p3.Z, p3.Z);      /* r = 2*(S2-S1) */

		Field::mul(U1, p1.X, Z2Z2);       /* U1 = X1*Z2Z2 */
		Field::mul(H, p2.X, Z1Z1);       /* U2 = X2*Z1Z1 */

		Field::sub(H, H, U1);              /* H = U2-U1 */

		Field::add(I, H, H);               /* 2*H */
		Field::sqr(I, I);                  /* I = (2*H)^2 */

		Field::mul(J, H, I);               /* J = H*I */
		Field::mul(S1, S1, J);             /* S1*J */

		Field::mul(p3.Y, U1, I);           /* V = U1*I */

		Field::sqr(p3.X, p3.Z);            /* r^2 */
		Field::sub(p3.X, p3.X, J);         /* r^2-J */
		Field::sub(p3.X, p3.X, p3.Y);
		Field::sub(p3.X, p3.X, p3.Y);      /* X3 = r^2-J-2*V */

		Field::sub(p3.Y, p3.Y, p3.X);      /* V-X3 */
		Field::mul(p3.Y, p3.Y, p3.Z);      /* r*(V-X3) */
		Field::sub(p3.Y, p3.Y, S1);
		Field::sub(p3.Y, p3.Y, S1);        /* Y3 = r*(V-X3)-2*S1*J */

		Field::add(p3.Z, p1.Z, p2.Z);    /* Z1+Z2 */
		Field::sqr(p3.Z, p3.Z);            /* (Z1+Z2)^2 */
		Field::sub(p3.Z, p3.Z, Z1Z1);      /* (Z1+Z2)^2-Z1Z1 */
		Field::sub(p3.Z, p3.Z, Z2Z2);      /* (Z1+Z2)^2-Z1Z1-Z2Z2 */
		Field::mul(p3.Z, p3.Z, H);         /* Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2)*H */

		vec_select(&p3, &p1, &p3, sizeof(Point), p2inf);
		vec_select(&out, &p2, &p3, sizeof(Point), p1inf);
	}

	/*
	 * https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
	 * with twist to handle either input at infinity, with |p1| encoded as Z==0,
	 * and |p2| as X==Y==0.
	 */
	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::add_affine(
		Point<Field>& out,
		const Point<Field>& p1,
		const Point_affine<Field>& p2)
	{
		//#if (defined(DEBUG) || !defined(NDEBUG)) && !defined(__CUDA_ARCH__)
		//		// Check if we're trying to add a point to itself in debug builds
		//		if (!vec_is_zero(a.Z, sizeof(a.Z))) {  // If a is not infinity
		//			Fp Z, ZZ;
		//			Fp::reciprocal(Z, a.Z);              // 1/Z
		//			Fp::sqr(ZZ, Z);                       // 1/Z^2
		//
		//			// Compare X coordinates: a.X/Z^2 == b.X?
		//			Fp a_X_affine;
		//			Fp::mul(a_X_affine, a.X, ZZ);        // X/Z^2
		//
		//			if (vec_is_equal(a_X_affine, b.X, sizeof(Fp))) {
		//				// Check Y coordinates as well
		//				Fp a_Y_affine;
		//				Fp::mul(ZZ, ZZ, Z);               // 1/Z^3
		//				Fp::mul(a_Y_affine, a.Y, ZZ);    // Y/Z^3
		//
		//				if (vec_is_equal(a_Y_affine, b.Y, sizeof(Fp))) {
		//					printf("ERROR: blst_p1_add_affine called with the same point. "
		//						"Use blst_p1_add_or_double_affine instead.\n");
		//					exit(1);
		//					// We could assert here, but let's just warn and proceed
		//					// assert(0 && "blst_p1_add_affine called with same point");
		//				}
		//			}
		//		}
		//#endif

		Point<Field> p3;
		Field Z1Z1, H, HH, I, J;
		bool_t p1inf, p2inf;

		p1inf = vec_is_zero(&p1.Z, sizeof(p1.Z));

		Field::sqr(Z1Z1, p1.Z);           /* Z1Z1 = Z1^2 */


		Field::mul(p3.Z, Z1Z1, p1.Z);     /* Z1*Z1Z1 */
		Field::mul(p3.Z, p3.Z, p2.Y);     /* S2 = Y2*Z1*Z1Z1 */

		p2inf = vec_is_zero(&p2.X, 2 * sizeof(p2.X));

		Field::mul(H, p2.X, Z1Z1);        /* U2 = X2*Z1Z1 */
		Field::sub(H, H, p1.X);           /* H = U2-X1 */

		Field::sqr(HH, H);                 /* HH = H^2 */
		Field::add(I, HH, HH);
		Field::add(I, I, I);               /* I = 4*HH */

		Field::mul(p3.Y, p1.X, I);        /* V = X1*I */
		Field::mul(J, H, I);               /* J = H*I */
		Field::mul(I, J, p1.Y);           /* Y1*J */

		Field::sub(p3.Z, p3.Z, p1.Y);     /* S2-Y1 */
		Field::add(p3.Z, p3.Z, p3.Z);      /* r = 2*(S2-Y1) */

		Field::sqr(p3.X, p3.Z);            /* r^2 */
		Field::sub(p3.X, p3.X, J);         /* r^2-J */
		Field::sub(p3.X, p3.X, p3.Y);
		Field::sub(p3.X, p3.X, p3.Y);      /* X3 = r^2-J-2*V */

		Field::sub(p3.Y, p3.Y, p3.X);      /* V-X3 */
		Field::mul(p3.Y, p3.Y, p3.Z);      /* r*(V-X3) */
		Field::sub(p3.Y, p3.Y, I);
		Field::sub(p3.Y, p3.Y, I);         /* Y3 = r*(V-X3)-2*Y1*J */

		Field::add(p3.Z, p1.Z, H);        /* Z1+H */
		Field::sqr(p3.Z, p3.Z);            /* (Z1+H)^2 */
		Field::sub(p3.Z, p3.Z, Z1Z1);      /* (Z1+H)^2-Z1Z1 */
		Field::sub(p3.Z, p3.Z, HH);        /* Z3 = (Z1+H)^2-Z1Z1-HH */

		vec_select(&p3.Z, &Fp::one(), &p3.Z, sizeof(p3.Z), p1inf);
		vec_select(&p3.X, &p2, &p3.X, 2 * sizeof(p3.X), p1inf);
		vec_select(&out, &p1, &p3, sizeof(Point), p2inf);
	}


	/*
	 * https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
	 */
	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::doubling(Point& p3, const Point& p1)
	{
		Field A, B, C;

		Field::sqr(A, p1.X);              /* A = X1^2 */
		Field::sqr(B, p1.Y);              /* B = Y1^2 */
		Field::sqr(C, B);                  /* C = B^2 */

		Field::add(B, B, p1.X);           /* X1+B */
		Field::sqr(B, B);                  /* (X1+B)^2 */
		Field::sub(B, B, A);               /* (X1+B)^2-A */
		Field::sub(B, B, C);               /* (X1+B)^2-A-C */
		Field::add(B, B, B);               /* D = 2*((X1+B)^2-A-C) */

		Field::mul_by_3(A, A);             /* E = 3*A */

		Field::sqr(p3.X, A);              /* F = E^2 */
		Field::sub(p3.X, p3.X, B);
		Field::sub(p3.X, p3.X, B);       /* X3 = F-2*D */

		Field::add(p3.Z, p1.Z, p1.Z);   /* 2*Z1 */
		Field::mul(p3.Z, p3.Z, p1.Y);   /* Z3 = 2*Z1*Y1 */

		Field::mul_by_8(C, C);             /* 8*C */
		Field::sub(p3.Y, B, p3.X);       /* D-X3 */
		Field::mul(p3.Y, p3.Y, A);       /* E*(D-X3) */
		Field::sub(p3.Y, p3.Y, C);       /* Y3 = E*(D-X3)-8*C */
	}

	/******************* !!!!! NOT CONSTANT TIME !!!!! *******************/

	/*
	 * http://hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#addition-add-2008-s
	 * http://hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#doubling-dbl-2008-s-1
	 * with twist to handle either input at infinity. Addition costs 12M+2S,
	 * while conditional doubling - 4M+6M+3S.
	 */
	template<typename Field>
	BLST_HD BLST_INLINE void Pointxyzz<Field>::dadd(Pointxyzz<Field>& p3, const Pointxyzz<Field>& p1,
		const Pointxyzz<Field>& p2)
	{
		Field U, S, P, R;

		if (vec_is_zero(p2.ZZZ, 2 * sizeof(p2.ZZZ))) {
			vec_copy(p3, p1, sizeof(*p3));
			return;
		}
		else if (vec_is_zero(p1.ZZZ, 2 * sizeof(p1.ZZZ))) {
			vec_copy(p3, p2, sizeof(*p3));
			return;
		}

		Field::mul(U, p1.X, p2.ZZ);              /* U1 = X1*ZZ2 */
		Field::mul(S, p1.Y, p2.ZZZ);             /* S1 = Y1*ZZZ2 */
		Field::mul(P, p2.X, p1.ZZ);              /* U2 = X2*ZZ1 */
		Field::mul(R, p2.Y, p1.ZZZ);             /* S2 = Y2*ZZZ1 */
		Field::sub(P, P, U);                       /* P = U2-U1 */
		Field::sub(R, R, S);                       /* R = S2-S1 */

		if (!vec_is_zero(P, sizeof(P))) {           /* X1!=X2 */
			Field PP, PPP, Q;                   /* add |p1| and |p2| */

			Field::sqr(PP, P);                     /* PP = P^2 */
			Field::mul(PPP, PP, P);                /* PPP = P*PP */
			Field::mul(Q, U, PP);                  /* Q = U1*PP */
			Field::sqr(p3.X, R);                  /* R^2 */
			Field::add(P, Q, Q);
			Field::sub(p3.X, p3.X, PPP);         /* R^2-PPP */
			Field::sub(p3.X, p3.X, P);           /* X3 = R^2-PPP-2*Q */
			Field::sub(Q, Q, p3.X);
			Field::mul(Q, Q, R);                   /* R*(Q-X3) */
			Field::mul(p3.Y, S, PPP);             /* S1*PPP */
			Field::sub(p3.Y, Q, p3.Y);           /* Y3 = R*(Q-X3)-S1*PPP */
			Field::mul(p3.ZZ, p1.ZZ, p2.ZZ);    /* ZZ1*ZZ2 */
			Field::mul(p3.ZZZ, p1.ZZZ, p2.ZZZ); /* ZZZ1*ZZZ2 */
			Field::mul(p3.ZZ, p3.ZZ, PP);        /* ZZ3 = ZZ1*ZZ2*PP */
			Field::mul(p3.ZZZ, p3.ZZZ, PPP);     /* ZZZ3 = ZZZ1*ZZZ2*PPP */
		}
		else if (vec_is_zero(R, sizeof(R))) {     /* X1==X2 && Y1==Y2 */
			Field V, W, M;                      /* double |p1| */

			Field::add(U, p1.Y, p1.Y);           /* U = 2*Y1 */
			Field::sqr(V, U);                      /* V = U^2 */
			Field::mul(W, V, U);                   /* W = U*V */
			Field::mul(S, p1.X, V);               /* S = X1*V */
			Field::sqr(M, p1.X);
			Field::mul_by_3(M, M);                 /* M = 3*X1^2[+a*ZZ1^2] */
			Field::sqr(p3.X, M);
			Field::add(U, S, S);                   /* 2*S */
			Field::sub(p3.X, p3.X, U);           /* X3 = M^2-2*S */
			Field::mul(p3.Y, W, p1.Y);           /* W*Y1 */
			Field::sub(S, S, p3.X);
			Field::mul(S, S, M);                   /* M*(S-X3) */
			Field::sub(p3.Y, S, p3.Y);           /* Y3 = M*(S-X3)-W*Y1 */
			Field::mul(p3.ZZ, p1.ZZ, V);         /* ZZ3 = V*ZZ1 */
			Field::mul(p3.ZZZ, p1.ZZZ, W);       /* ZZ3 = W*ZZZ1 */
		}
		else {                                    /* X1==X2 && Y1==-Y2 */
			vec_zero(p3.ZZZ, 2 * sizeof(p3.ZZZ));   /* set |p3| to infinity */
		}
	}

	/*
	 * http://hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#addition-madd-2008-s
	 * http://hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#doubling-mdbl-2008-s-1
	 * with twists to handle even subtractions and either input at infinity.
	 * Addition costs 8M+2S, while conditional doubling - 2M+4M+3S.
	 */
	template<typename Field>
	BLST_HD BLST_INLINE void Pointxyzz<Field>::dadd_affine(Pointxyzz<Field>& p3, const Pointxyzz<Field>& p1,
		const Point_affine<Field>& p2,
		bool_t subtract)
	{
		Field P, R;

		if (vec_is_zero(p2, sizeof(*p2))) {
			vec_copy(p3, p1, sizeof(*p3));
			return;
		}
		else if (vec_is_zero(p1.ZZZ, 2 * sizeof(p1.ZZZ))) {
			vec_copy(&p3.X, &p2.X, 2 * sizeof(p3.X));
			Field::cneg(p3.ZZZ, Field::one(), subtract);
			vec_copy(&p3.ZZ, &Field::one(), sizeof(p3.ZZ));
			return;
		}

		Field::mul(P, p2.X, p1.ZZ);              /* U2 = X2*ZZ1 */
		Field::mul(R, p2.Y, p1.ZZZ);             /* S2 = Y2*ZZZ1 */
		Field::cneg(R, R, subtract);
		Field::sub(P, P, p1.X);                   /* P = U2-X1 */
		Field::sub(R, R, p1.Y);                   /* R = S2-Y1 */

		if (!vec_is_zero(P, sizeof(P))) {           /* X1!=X2 */
			Field PP, PPP, Q;                   /* add |p2| to |p1| */

			Field::sqr(PP, P);                     /* PP = P^2 */
			Field::mul(PPP, PP, P);                /* PPP = P*PP */
			Field::mul(Q, p1.X, PP);              /* Q = X1*PP */
			Field::sqr(p3.X, R);                  /* R^2 */
			Field::add(P, Q, Q);
			Field::sub(p3.X, p3.X, PPP);         /* R^2-PPP */
			Field::sub(p3.X, p3.X, P);           /* X3 = R^2-PPP-2*Q */
			Field::sub(Q, Q, p3.X);
			Field::mul(Q, Q, R);                   /* R*(Q-X3) */
			Field::mul(p3.Y, p1.Y, PPP);         /* Y1*PPP */
			Field::sub(p3.Y, Q, p3.Y);           /* Y3 = R*(Q-X3)-Y1*PPP */
			Field::mul(p3.ZZ, p1.ZZ, PP);        /* ZZ3 = ZZ1*PP */
			Field::mul(p3.ZZZ, p1.ZZZ, PPP);     /* ZZZ3 = ZZZ1*PPP */
		}
		else if (vec_is_zero(R, sizeof(R))) {     /* X1==X2 && Y1==Y2 */
			Field U, S, M;                      /* double |p2| */

			Field::add(U, p2.Y, p2.Y);           /* U = 2*Y1 */
			Field::sqr(p3.ZZ, U);                 /* [ZZ3 =] V = U^2 */
			Field::mul(p3.ZZZ, p3.ZZ, U);        /* [ZZZ3 =] W = U*V */
			Field::mul(S, p2.X, p3.ZZ);          /* S = X1*V */
			Field::sqr(M, p2.X);
			Field::mul_by_3(M, M);                 /* M = 3*X1^2[+a] */
			Field::sqr(p3.X, M);
			Field::add(U, S, S);                   /* 2*S */
			Field::sub(p3.X, p3.X, U);           /* X3 = M^2-2*S */
			Field::mul(p3.Y, p3.ZZZ, p2.Y);     /* W*Y1 */
			Field::sub(S, S, p3.X);
			Field::mul(S, S, M);                   /* M*(S-X3) */
			Field::sub(p3.Y, S, p3.Y);           /* Y3 = M*(S-X3)-W*Y1 */
			Field::cneg(p3.ZZZ, p3.ZZZ, subtract);
		}
		else {                                    /* X1==X2 && Y1==-Y2 */
			vec_zero(p3.ZZZ, 2 * sizeof(p3.ZZZ));   /* set |p3| to infinity */
		}
	}

	template<typename Field>
	BLST_HD BLST_INLINE void Pointxyzz<Field>::to_Jacobian(Point<Field>& out, const Pointxyzz<Field>& in)
	{
		Field::mul(out.X, in.X, in.ZZ);
		Field::mul(out.Y, in.Y, in.ZZZ);
		vec_copy(out.Z, in.ZZ, sizeof(out.Z));
	}

	template<typename Field>
	BLST_HD BLST_INLINE void Point<Field>::to_xyzz(Pointxyzz<Field>& out, const Point<Field>& in)
	{
		vec_copy(out.X, in.X, 2 * sizeof(out.X));
		Field::sqr(out.ZZ, in.Z);
		Field::mul(out.ZZZ, out.ZZ, in.Z);
	}



}
#endif
