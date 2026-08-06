/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_SCALAR_FP12_H
#define BLST_SCALAR_FP12_H


#include "blst/common/defines.h"
#include "blst/backend/backend.h"
#include "blst/scalar/fields/Fp2.h"



namespace blst::scalar
{

	struct Fp4 {
		Fp2 data[2];
		BLST_HD Fp2& operator[](size_t idx) {
			return data[idx];
		}
		BLST_HD const Fp2& operator[](size_t idx)const {
			return data[idx];
		}
	};

	/**
	 * @brief Sextic extension field Fp6 represented as `c0 + c1*v + c2*v^2`.
	 *
	 * Each coefficient is an Fp2 value in Montgomery form. This type primarily
	 * provides the ABI-stable storage used by Fp12 and pairing operations.
	 */
	struct Fp6
	{
		Fp2 fp2[3];

		/// Access coefficient @p idx in increasing powers of v.
		BLST_HD Fp2& operator[](size_t idx) { return fp2[idx]; }
		BLST_HD const Fp2& operator[](size_t idx) const { return fp2[idx]; }
	};

	/**
	 * @brief Precomputed variable-time GLV-4 exponentiation schedule.
	 *
	 * This is trusted bytecode produced by `Fp12::glv4_compile`. Passing a
	 * schedule after modifying its fields violates the executor precondition.
	 */
	struct GLV4_Exp
	{
		uint8_t sign_mask;     ///< Basis elements that require conjugation.
		uint8_t sched_len;     ///< Number of valid schedule entries.
		uint8_t schedule[70];  ///< Compact square/multiply bytecode.
	};

	/**
	 * @brief Precomputed variable-time GLV-2 joint-sparse-form schedule.
	 *
	 * Each index stores two signed digits; two digit pairs are packed per byte.
	 * This is trusted bytecode produced by `Fp12::glv_recode_jsf`. Passing a
	 * schedule after modifying its fields violates the executor precondition.
	 */
	struct GLV_JSF_Exp
	{
		uint8_t sched_len;     ///< Number of valid digit pairs.
		uint8_t schedule[71];  ///< Packed signed-digit schedule.

		/// Set the @p idx digit pair to (`d1`, `d2`).
		BLST_HD void set(size_t idx, int8_t d1, int8_t d2)
		{
			unsigned enc = static_cast<unsigned>(d1 + 1) |
				(static_cast<unsigned>(d2 + 1) << 2);
			const unsigned shift = static_cast<unsigned>((idx % 2) * 4);
			enc <<= shift;
			enc |= schedule[idx / 2] & ~(0xFU << shift);
			schedule[idx / 2] = static_cast<uint8_t>(enc);
		}

		struct Pair { int8_t d1; int8_t d2; };

		/// Return the @p idx signed-digit pair.
		BLST_HD Pair get(size_t idx) const
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
	 * Ordinary arithmetic supports exact output aliasing with an input; partial
	 * overlap is not supported. Inputs must contain reduced Montgomery
	 * coefficients. Debug assertions do not constitute Release validation. Host
	 * operations allocate nothing and are nonthrowing for valid inputs.
	 */
	struct Fp12
	{
		using GLV4_Exp = scalar::GLV4_Exp;
		using GLV_JSF_Exp = scalar::GLV_JSF_Exp;

		Fp6 fp6[2];

		/// Access Fp6 coefficient @p idx.
		BLST_HD Fp6& operator[](size_t idx) { return fp6[idx]; }
		BLST_HD const Fp6& operator[](size_t idx) const { return fp6[idx]; }

		/// @brief `ret = a + b` in Fp12.
		BLST_HD static void add(Fp12& ret, const Fp12& a, const Fp12& b);
		/// @brief `ret = a - b` in Fp12.
		BLST_HD static void sub(Fp12& ret, const Fp12& a, const Fp12& b);
		/// @brief `ret = a^2` using the scalar implementation.
		BLST_HD static void sqr(Fp12& ret, const Fp12& a);
		/// @brief `ret = a^2` using scalar code generation explicitly.
		BLST_HD static void sqr_scalar(Fp12& ret, const Fp12& a);
		/// @brief `ret = a * b` using the scalar implementation.
		BLST_HD static void mul(Fp12& ret, const Fp12& a, const Fp12& b);
		/// @brief `ret = a * b` using scalar code generation explicitly.
		BLST_HD static void mul_scalar(Fp12& ret, const Fp12& a, const Fp12& b);
		/// @brief Multiply by @p b or its conjugate.
		/// @param conj_b When true, compute `ret = a * conjugate(b)`.
		BLST_HD static void mul(Fp12& ret, const Fp12& a, const Fp12& b, bool conj_b);
		/// @brief Standard exponentiation `ret = a^exp`.
		/// @param exp Little-endian exponent bytes.
		/// @param bits Number of meaningful exponent bits; zero produces one and permits a null exponent.
		BLST_HD static void exp(Fp12& ret, const Fp12& a, const byte* exp, size_t bits);
		/// @brief Multiplicative identity of Fp12.
		BLST_HD static const Fp12& one();
		/// @brief Cyclotomic squaring using the scalar implementation.
		/// @pre @p a lies in the cyclotomic subgroup.
		BLST_HD static void cyclotomic_sqr(Fp12& ret, const Fp12& a);
		/// @brief Cyclotomic squaring using scalar code generation explicitly.
		/// @pre @p a lies in the cyclotomic subgroup.
		BLST_HD static void cyclotomic_sqr_scalar(Fp12& ret, const Fp12& a);
		/// @brief Constant-time exponentiation of a G_T element.
		/// @pre @p a lies in the order-r target group G_T.
		/// @details Runtime and memory access are independent of exponent contents
		/// for a fixed public @p bits.
		BLST_HD static void cyclotomic_exp(Fp12& ret, const Fp12& a, const byte* exp, size_t bits);
		/// @brief Optimized variable-time exponentiation of a G_T element.
		/// @pre @p a lies in the order-r target group G_T.
		/// @warning Use only with public exponent data.
		BLST_HD static void cyclotomic_exp_vartime(Fp12& ret, const Fp12& a, const byte* exp, size_t bits);
		/// @brief Basic variable-time exponentiation of a G_T element using cyclotomic squaring.
		/// @pre @p a lies in the order-r target group G_T.
		/// @warning Use only with public exponent data.
		BLST_HD static void cyclotomic_exp_basic(Fp12& ret, const Fp12& a, const byte* exp, size_t bits);
		/// @brief Split a 256-bit scalar into two approximately 128-bit GLV components.
		/// @details The decomposition satisfies `k = k1 + k2*lambda` modulo r.
		BLST_HD static void glv_split_2(limb_t k1_out[4], limb_t k2_out[4], const limb_t k[4]);
		/// @brief Recode two GLV components into joint sparse form.
		BLST_HD static GLV_JSF_Exp glv_recode_jsf(const limb_t* k1, const limb_t* k2);
		/// @brief GLV-2 G_T exponentiation from a 32-byte scalar.
		/// @pre @p base lies in the order-r target group G_T.
		BLST_HD static void cyclotomic_exp_glv(Fp12& out, const Fp12& base, const byte* scalar_bytes);
		/// @brief GLV-2 G_T exponentiation with pre-split scalars.
		/// @pre @p base lies in the order-r target group G_T.
		BLST_HD static void cyclotomic_exp_glv(Fp12& out, const Fp12& base, const limb_t* k1, const limb_t* k2);
		/// @brief GLV-2 G_T exponentiation with on-the-fly JSF recoding.
		/// @pre @p base lies in the order-r target group G_T.
		BLST_HD static void cyclotomic_exp_glv_jsf(Fp12& out, const Fp12& base, const limb_t* k1, const limb_t* k2);
		/// @brief GLV-2 G_T exponentiation with a trusted precomputed JSF schedule.
		/// @pre @p base lies in G_T and @p sched is an unmodified result of `glv_recode_jsf`.
		BLST_HD static void cyclotomic_exp_glv_jsf(Fp12& out, const Fp12& base, const GLV_JSF_Exp& sched);
		/// @brief Four-dimensional G_T scalar decomposition.
		/// @details `k = k0 + k1*lambda1 + k2*lambda2 + k3*lambda3` modulo r.
		BLST_HD static void gt_split_scalar_4d(int64_t k_out[4], const limb_t k_in[4]);
		/// @brief Constant-time four-dimensional G_T scalar decomposition.
		BLST_HD static void gt_split_scalar_4d_ct(uint64_t k_out[4], const limb_t k_in[4]);
		/// @brief Compile a scalar into a variable-time GLV-4 execution schedule.
		BLST_HD static GLV4_Exp glv4_compile(const limb_t* k_in);
		/// @brief GLV-4 G_T exponentiation with a trusted precompiled schedule.
		/// @pre @p base lies in G_T and @p props is an unmodified result of `glv4_compile`.
		BLST_HD static void cyclotomic_exp_glv4(Fp12& out, const Fp12& base, const GLV4_Exp& props);
		/// @brief GLV-4 G_T exponentiation with a runtime scalar split.
		/// @details The full 256-bit unsigned scalar is accepted; decomposition
		/// preserves its value modulo the order of G_T.
		/// @pre @p base lies in the order-r target group G_T.
		BLST_HD static void cyclotomic_exp_glv4(Fp12& out, const Fp12& base, const limb_t* scalar);
		/// @brief Constant-time GLV-4 G_T exponentiation over a fixed 64-round schedule.
		/// @pre @p base lies in the order-r target group G_T.
		BLST_HD static void cyclotomic_exp_glv4_ct(Fp12& out, const Fp12& base, const limb_t* scalar);
		/// @brief Conjugate @p a in place by negating its w coefficient.
		BLST_HD static void conjugate(Fp12& a);
		/// @brief `ret = a^{-1}` in Fp12.
		/// @pre @p a is nonzero; zero has no defined field inverse.
		BLST_HD static void inverse(Fp12& ret, const Fp12& a);
		/// @brief Apply the Frobenius map: `ret = a^(p^n)`.
		/// @pre `1 <= n && n <= 3`.
		BLST_HD static void frobenius_map(Fp12& ret, const Fp12& a, size_t n);
		/// @brief Multiply by a sparse `xy00z0` element using the scalar implementation.
		BLST_HD static void mul_by_xy00z0(Fp12& ret, const Fp12& a, const Fp6& xy00z0);
		/// @brief Multiply by a sparse `xy00z0` element using scalar code generation.
		BLST_HD static void mul_by_xy00z0_scalar(Fp12& ret, const Fp12& a, const Fp6& xy00z0);
		/// @brief Constant-time equality predicate.
		BLST_HD static int equal(const Fp12& a, const Fp12& b);
		/// @brief Constant-time predicate for the multiplicative identity.
		BLST_HD static int is_one(const Fp12& a);
		/**
		 * @brief Decode 576 big-endian bytes with canonical Fp coefficients.
		 * @return `true` on success; `false` if any 48-byte coefficient is at
		 * least p.
		 * @pre @p serialized points to 576 readable bytes.
		 * @post Success writes the decoded Montgomery-form Fp12 value. Failure
		 * leaves @p ret unchanged. No Fp12 or target-group membership check is
		 * performed. Input bytes may overlap @p ret because decoding completes
		 * into temporary storage before publication.
		 */
		[[nodiscard]] BLST_HD static bool from_bendian(
			Fp12& ret, const byte serialized[12 * 48]);
		/// @brief Convert from Montgomery form and encode 576 big-endian bytes.
		BLST_HD static void to_bendian(byte serialized[12 * 48], const Fp12& a);
		/// @brief Return one iff @p f is in the cyclotomic subgroup G_phi_12.
		BLST_HD static int in_group(const Fp12& f);
		/// @brief Faster cyclotomic-subgroup predicate using Frobenius relations.
		BLST_HD static limb_t is_cyclotomic(const Fp12& f);
		/// @brief Compute `ret = (ret * a)^(2^n)`.
		BLST_HD static void mul_n_sqr(Fp12& ret, const Fp12& a, size_t n);
		/// @brief Final-exponentiation helper: `ret = a^((z-1)/2)`.
		BLST_HD static void raise_to_z_div_by_2(Fp12& ret, const Fp12& a);
		/// @brief Final-exponentiation helper: `a = b^z`.
		BLST_HD static void raise_to_z(Fp12& a, const Fp12& b);
		/// @brief Perform pairing final exponentiation into G_T.
		/// @pre @p f is a nonzero general Fp12 or Miller-loop value.
		/// @post @p ret lies in G_T.
		BLST_HD static void final_exp(Fp12& ret, const Fp12& f);
	};

	#define BLST_FP12_ONE_INITIALIZER                                      \
		Fp12{                                                              \
			Fp6{ Fp2{ Fp{ ONE_MONT_P_LIMBS }, Fp{} }, Fp2{}, Fp2{} },     \
			Fp6{}                                                           \
		}

	constexpr Fp12 BLST_Rx_HOST = BLST_FP12_ONE_INITIALIZER;

#if defined(BLST_CUDA) && defined(__CUDACC__)
	static __constant__ Fp12 BLST_Rx_DEV = BLST_FP12_ONE_INITIALIZER;
#endif

	/// Return the Fp12 multiplicative identity in host or device storage.
	BLST_HD BLST_FORCEINLINE const Fp12& get_BLST_Rx()
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_Rx_DEV;
#else
		return BLST_Rx_HOST;
#endif
	}

	#undef BLST_FP12_ONE_INITIALIZER

}

#endif /* BLST_SCALAR_FP12_H */



