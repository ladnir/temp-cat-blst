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
	using G1Point = ::blst_p1;
	using G1AffinePoint = ::blst_p1_affine;
	using G1PointZZ = Pointxyzz<Fp>;

	using G2Point = ::blst_p2;
	using G2AffinePoint = ::blst_p2_affine;
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
