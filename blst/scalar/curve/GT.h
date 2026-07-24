#pragma once

#include "blst/scalar/fields/Fp12.h"
#include "blst/scalar/curve/point.h"

namespace blst::scalar {

	/**
	 * @brief GT group wrapper around Fp12 multiplicative group
	 * @details GT is the target group in pairing-based cryptography where:
	 * - Group operation is addition (maps to Fp12 multiplication)
	 * - Identity is Fp12::one()
	 * - Doubling maps to cyclotomic squaring
	 * - Scalar multiplication maps to exponentiation
	 */
	struct GT
	{
		using Field_type = Fp12;

		constexpr GT() = default;
		constexpr GT(const GT&) = default;
		constexpr GT& operator=(const GT&) = default;

		/// Construct from Fp12 element
		constexpr GT(const Fp12& f) : data(f) {}

		Fp12 data;

		BLST_HD auto operator[](size_t idx) -> decltype(data[idx]) {
			return data[idx];
		}
		BLST_HD auto operator[](size_t idx) const -> decltype(data[idx]) {
			return data[idx];
		}

		/**
		 * @brief Group addition: out = a + b (maps to Fp12 multiplication)
		 * @param[out] out  Result
		 * @param[in]  a    Left operand
		 * @param[in]  b    Right operand
		 */
		static BLST_HD void add(GT& out, const GT& a, const GT& b) {
			Fp12::mul(out.data, a.data, b.data);
		}


		/**
		 * @brief Group double or addition: out = a + b (maps to Fp12 multiplication).
		 * Forwards to add, provided for API symmetry.
		 * @param[out] out  Result
		 * @param[in]  a    Left operand
		 * @param[in]  b    Right operand
		 */
		static BLST_HD void dadd(GT& out, const GT& a, const GT& b) {
			add(out, a, b);
		}

		///**
		// * @brief Group subtraction: out = a - b (maps to Fp12 division)
		// * @param[out] out  Result
		// * @param[in]  a    Left operand
		// * @param[in]  b    Right operand
		// */
		//static BLST_HD void sub(GT& out, const GT& a, const GT& b) {
		//	Fp12 b_inv;
		//	Fp12::inverse(b_inv, b.data);
		//	Fp12::mul(out.data, a.data, b_inv);
		//}

		/**
		 * @brief Group subtraction: out = a - b (maps to Fp12 division)
		 * @param[out] out  Result
		 * @param[in]  a    Left operand
		 * @param[in]  b    Right operand
		 */
		static BLST_HD void sub(GT& out, const GT& a, const GT& b) {
			Fp12 b_inv = b.data;
			Fp12::conjugate(b_inv); // since in GT, inversion = conjugation
			Fp12::mul(out.data, a.data, b_inv);
		}

		/**
		 * @brief Group doubling: out = 2 * a (maps to cyclotomic squaring)
		 * @param[out] out  Result
		 * @param[in]  a    Input element
		 */
		static BLST_HD void doubling(GT& out, const GT& a) {
			Fp12::cyclotomic_sqr(out.data, a.data);
			//Fp12::sqr(out.data, a.data);
		}

		/**
		 * @brief Repeated doubling: out = 2^n * (out + a)
		 * @param[out] out  Result
		 * @param[in]  a    Input element
		 * @param[in]  n    Number of doublings
		 */
		static BLST_HD void add_n_dbl(GT& out, const GT& a, size_t n) {

			add(out, out, a);
			for (size_t i = 0; i < n; i++) {
				doubling(out, out);
			}
		}

		/**
		 * @brief Scalar multiplication: out = scalar * a (maps to Fp12 exponentiation)
		 * @param[out] out     Result
		 * @param[in]  a       Base element
		 * @param[in]  scalar  Little-endian scalar bytes
		 * @param[in]  nbits   Number of meaningful bits in scalar
		 */
		static BLST_HD void mult(GT& out, const GT& a, const byte* scalar, size_t nbits) {
			if (nbits  == 256) {
				// For 256-bit scalars, we can use the optimized GLV-4 path.
				// This is a common case for BLS signatures where scalars are 256 bits.
				// The public byte-string API does not promise limb alignment.
				alignas(limb_t) byte aligned_scalar[32];
				vec_copy(aligned_scalar, scalar, sizeof(aligned_scalar));
				Fp12::cyclotomic_exp_glv4(
					out.data,
					a.data,
					reinterpret_cast<const limb_t*>(aligned_scalar));
			} else {
				// For non-256-bit scalars, fall back to the basic cyclotomic exponentiation.
				Fp12::cyclotomic_exp(out.data, a.data, scalar, nbits);
			}
		}

		/**
		 * @brief Conditional negate in constant time if cbit is 1
		 * @details In GT, negation maps to inversion in Fp12
		 */
		static BLST_HD void cneg(GT& p, limb_t cbit) {
			if (cbit) {
				Fp12::conjugate(p.data);
			}
		}

		/**
		 * @brief Negate: out = -p (maps to Fp12 inversion)
		 */
		static BLST_HD void neg(GT& out, const GT& p) {
			out.data = p.data;
			Fp12::conjugate(out.data);
		}

		/**
		 * @brief Constant-time conditional swap: if cbit==1, swap a and b
		 */
		static BLST_HD void cswap(GT& BLST_restrict a, GT& BLST_restrict b, bool_t cbit) {
			// Implement constant-time swap for Fp12 data
			// This would need to be implemented at the Fp12 level
			vec_cswap(&a, &b, sizeof(GT), cbit);
		}

		/**
		 * @brief Constant-time conditional copy: if cbit==1, a <- b
		 */
		static BLST_HD void ccopy(GT& BLST_restrict a, const GT& BLST_restrict b, bool_t cbit) {
			vec_select(&a, &b, &a, sizeof(GT), cbit);
		}

		/**
		 * @brief Equality check
		 */
		static BLST_HD limb_t equal(const GT& a, const GT& b) {
			return Fp12::equal(a.data, b.data);
		}

		/**
		 * @brief Canonical GT generator pairing(G1::generator(), G2::generator()).
		 * @details Computed once on first host use and returned from immutable
		 * static storage.
		 */
		static const GT& generator();

		/**
		 * @brief Returns the identity element (Fp12::one())
		 */
		static BLST_HD GT identity() {
			return GT(Fp12::one());
		}

		/**
		 * @brief Check if element is the identity
		 */
		static BLST_HD int is_identity(const GT& p) {
			return Fp12::is_one(p.data);
		}

		/**
		 * @brief Check if element is in the correct subgroup
		 */
		static BLST_HD bool_t in_group(const GT& p) {
			return Fp12::in_group(p.data);
		}

		/**
		 * @brief Check if element is cyclotomic
		 */
		static BLST_HD bool_t is_cyclotomic(const GT& p) {
			return Fp12::is_cyclotomic(p.data);
		}

		/** @name Serialization (big-endian) */
		///@{
		/**
		 * @brief Serialize to big-endian bytes
		 */
		static BLST_HD void serialize(unsigned char serialized[12 * 48], const GT& p) {
			Fp12::to_bendian(serialized, p.data);
		}

		/**
		 * @brief Deserialize from big-endian bytes
		 */
		static BLST_HD void deserialize(GT& out, const unsigned char serialized[12 * 48]) {
			Fp12::from_bendian(out.data, serialized);
		}
		///@}

		/**
		 * @brief Apply final multiplication aka exponentiation for pairing results.
		 */
		static BLST_HD void final_mult(GT& out, const GT& f) {
			Fp12::final_exp(out.data, f.data);
		}

		/**
		 * @brief Raise to z (BLS12-381 specific)
		 */
		static BLST_HD void mult_by_z(GT& out, const GT& a) {
			Fp12::raise_to_z(out.data, a.data);
		}

		/**
		 * @brief Raise to z/2 (BLS12-381 specific)
		 */
		static BLST_HD void mult_by_z_div_by_2(GT& out, const GT& a) {
			Fp12::raise_to_z_div_by_2(out.data, a.data);
		}

		/**
		 * @brief Apply Frobenius map
		 * @param[out] out  Result
		 * @param[in]  a    Input element
		 * @param[in]  n    Power of Frobenius (1, 2, or 3)
		 */
		static BLST_HD void frobenius_map(GT& out, const GT& a, size_t n) {
			Fp12::frobenius_map(out.data, a.data, n);
		}

		/**
		 * @brief Conjugate (Frobenius with n=6)
		 */
			static BLST_HD void conjugate(GT& a) {
				Fp12::conjugate(a.data);
			}


			static BLST_HD void miller_loop(GT& ret, const G2AffinePoint* Q,
				const G1AffinePoint* P);

		static BLST_HD void miller_loop_n(
			GT& ret,
			const G2AffinePoint Q[],
			const G1AffinePoint P[],
			size_t n);

		static BLST_HD void miller_loop_n(
			GT& out,
			const G2AffinePoint* const Qs[],
			const G1AffinePoint* const Ps[],
			size_t n);

		static BLST_HD void precompute_lines(Fp6 Qlines[68], const G2AffinePoint& Q);

		static BLST_HD void miller_loop_lines(GT& ret, const Fp6 Qlines[68],
			const G1AffinePoint& P);

		static BLST_HD void miller_loop_lines_n(
			GT& out,
			const Fp6* const Qlines[],
			const G1AffinePoint* const Ps[],
			size_t n);

		static BLST_HD void miller_loop_lines_n_preprocessed(
			GT& out,
			const Fp6* const Qlines[],
			const G1AffinePoint Px2[],
			size_t n);

		//     /**
			  //* @brief add n doublings: out = 2^n * (a + out)
		//      */
		//     static BLST_HD void add_n_dbl(GT& out, const GT& a, size_t n) {
			//         Fp12::mul_n_sqr(out.data, a.data, n);
			//     }
		};

	} // namespace blst
