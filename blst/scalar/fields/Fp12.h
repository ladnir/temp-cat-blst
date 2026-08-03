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

	struct Fp6 {
		Fp2 data[3];

		BLST_HD BLST_FORCEINLINE auto& operator[](size_t idx) {
			return data[idx];
		}
		BLST_HD BLST_FORCEINLINE auto& operator[](size_t idx) const
		{
			return data[idx];
		}
	};

	/// @brief Precomputed schedule for constant-time GLV-4 exponentiation.
	/// Contains a compact bytecode sequence of operations (squares/multiplies).
	struct GLV4_Exp {
		uint8_t sign_mask;     ///< Bitmask indicating which basis elements need conjugation.
		uint8_t sched_len;     ///< Length of the valid schedule.
		uint8_t schedule[70];  ///< Bytecode for the exponentiation loop.
	};

	/// @brief Precomputed schedule for constant-time GLV-2 JSF exponentiation.
	/// Contains a compact bytecode sequence of operations (squares/multiplies).
	/// The ecnoding has two signed digits per index. Then we pack two such
	/// pairs into a byte for storage efficiency.
	struct GLV_JSF_Exp {
		uint8_t sched_len;     ///< Length of the valid schedule.
		uint8_t schedule[71];  ///< Bytecode for the exponentiation loop.

		// sets the idx-th digit pair to (d1, d2)
		BLST_HD BLST_FORCEINLINE void set(size_t idx, int8_t d1, int8_t d2) {
			auto enc = ((d1 + 1) | ((d2 + 1) << 2));
			enc <<= ((idx % 2) * 4);
			enc |= schedule[idx / 2] & (~(0xF << ((idx % 2) * 4)));
			schedule[idx / 2] = static_cast<uint8_t>(enc);
		}

		struct Pair {
			int8_t d1;
			int8_t d2;
		};

		// gets the idx-th digit pair (d1, d2)
		BLST_HD BLST_FORCEINLINE Pair get(size_t idx) const {
			auto enc = (schedule[idx / 2] >> ((idx % 2) * 4)) & 0xF;
			auto d1 = (int8_t)((enc & 0x3) - 1);
			auto d2 = (int8_t)(((enc >> 2) & 0x3) - 1);
			return { d1, d2 };
		}
	};

	struct Fp12
	{
		Fp6 data[2];

		BLST_HD BLST_FORCEINLINE auto& operator[](size_t idx) {
			return data[idx];
		}

		BLST_HD BLST_FORCEINLINE auto& operator[](size_t idx) const {
			return data[idx];
		}

		/// @brief ret = a + b in Fp12.
		BLST_HD static void add(Fp12& ret, const Fp12& a, const Fp12& b);

		/// @brief ret = a - b in Fp12.
		BLST_HD static void sub(Fp12& ret, const Fp12& a, const Fp12& b);

		/// @brief ret = a^2 in Fp12.
		BLST_HD static void sqr(Fp12& ret, const Fp12& a);

        /// @brief ret = a^2 in Fp12. Performed with scalar codegen (no AVX).
		BLST_HD static void sqr_scalar(Fp12& ret, const Fp12& a);

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
        /// @brief ret = a^2 in Fp12. Performed with AVX codegen.
        static void sqr_avx(Fp12& ret, const Fp12& a);
#endif

		/// @brief ret = a * b in Fp12.
		BLST_HD static void mul(Fp12& ret, const Fp12& a, const Fp12& b);

        /// @brief ret = a * b in Fp12, performed with scalar codegen (no AVX).
		BLST_HD static void mul_scalar(Fp12& ret, const Fp12& a, const Fp12& b);

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
        /// @brief ret = a * b in Fp12, performed with AVX codegen.
        static void mul_avx(Fp12& ret, const Fp12& a, const Fp12& b);
#endif
		
		/// @brief ret = a * b in Fp12.
		/// @param conj_b If true, multiplies by conjugate of b (a * conj(b)).
		BLST_HD static void mul(Fp12& ret, const Fp12& a, const Fp12& b, bool conj_b);

		/// @brief ret = a^exp in Fp12 (Standard Exponentiation).
		/// @param exp Little-endian exponent bytes.
		/// @param bits Number of bits in the exponent.
		BLST_HD static void exp(Fp12& ret, const Fp12& a, const byte* exp, size_t bits);


		/// @brief Returns the multiplicative identity (1) of Fp12.
		BLST_HD static const Fp12& one();

		// --- Cyclotomic Subgroup Operations ---

		/// @brief ret = a^2 in Fp12, optimized for the Cyclotomic subgroup.
		/// Faster than standard sqr() but requires input `a` to be in G_phi_12.
		BLST_HD static void cyclotomic_sqr(Fp12& ret, const Fp12& a);

		/// @brief ret = a^2 in Fp12, optimized for the Cyclotomic subgroup.
		/// Faster than standard sqr() but requires input `a` to be in G_phi_12.
		/// Performed with scalar codegen (no AVX).
		BLST_HD static void cyclotomic_sqr_scalar(Fp12& ret, const Fp12& a);
		
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
        /// @brief ret = a^2 in Fp12, optimized for the Cyclotomic subgroup.
        /// Faster than standard sqr() but requires input `a` to be in G_phi_12.
        /// Performed with AVX codegen.
        static void cyclotomic_sqr_avx(Fp12& ret, const Fp12& a);
#endif

		/// @brief ret = a^exp in Fp12, using Cyclotomic squaring.
		/// @note Input `a` MUST be in the cyclotomic subgroup.
		BLST_HD static void cyclotomic_exp(Fp12& ret, const Fp12& a, const byte* exp, size_t bits);

		/// @brief ret = a^exp in Fp12, using Cyclotomic squaring.
		/// @note Input `a` MUST be in the cyclotomic subgroup.
		BLST_HD static void cyclotomic_exp_basic(Fp12& ret, const Fp12& a, const byte* exp, size_t bits);

		// --- GLV-2 Operations (Standard Endomorphism) ---

		/// @brief 2-Way GLV Scalar Splitter.
		/// Decomposes scalar k into k1, k2 such that k = k1 + k2*lambda.
		/// @param k1_out Output limb array for k1.
		/// @param k2_out Output limb array for k2.
		/// @param k Input scalar (256-bit).
		BLST_HD static void glv_split_2(limb_t k1_out[4], limb_t k2_out[4], const limb_t k[4]);

		/// @brief Recodes GLV-2 scalars into Joint Sparse Form (JSF).
		/// @param jsf Output buffer for signed digits.
		/// @param k1 First split scalar.
		/// @param k2 Second split scalar (multiplies lambda).
		/// @return Length of the JSF sequence.
		BLST_HD static GLV_JSF_Exp glv_recode_jsf(const limb_t* k1, const limb_t* k2);

		/// @brief ret = a^scalar in Fp12, using GLV-2 optimization.
		/// Automatically splits the scalar internally. Must be 256 bits.
		/// base MUST be in the cyclotomic subgroup.
		BLST_HD static void cyclotomic_exp_glv(Fp12& out, const Fp12& base, const byte* scalar_bytes);

		/// @brief ret = a^k1 * phi(a)^k2 in Fp12.
		/// Lower-level GLV-2 function taking pre-split scalars.
		/// k1, k2 are expected to be approx 128 bits each.
		/// base MUST be in the cyclotomic subgroup.
		BLST_HD static void cyclotomic_exp_glv(Fp12& out, const Fp12& base, const limb_t* k1, const limb_t* k2);

		/// @brief ret = a^k1 * phi(a)^k2 in Fp12, using on-the-fly JSF recoding.
		/// Fastest 2-way GLV path with runtime recoding.
		/// MUST be in the cyclotomic subgroup.
		BLST_HD static void cyclotomic_exp_glv_jsf(Fp12& out, const Fp12& base, const limb_t* k1, const limb_t* k2);


		/// @brief ret = a^k1 * phi(a)^k2 in Fp12, using precomputed JSF form.
		/// Fastest 2-way GLV path if JSF is precomputed.
		/// base MUST be in the cyclotomic subgroup.
		BLST_HD static void cyclotomic_exp_glv_jsf(Fp12& out, const Fp12& base, const GLV_JSF_Exp& sched);

		// 4-dimensional GLV split:
		//   k = k0 + k1*lambda1 + k2*lambda2 + k3*lambda3
		// where lambda1, lambda2, lambda3 are the eigenvalues of the endomorphism.
		BLST_HD static void gt_split_scalar_4d(int64_t k_out[4], const limb_t k_in[4]);

		// --- GLV-4 Operations (4-Dimensional Decomposition) ---

		/// @brief 4-Way GLV Scalar Splitter.
		/// Decomposes scalar k into k0..k3 using full endomorphism basis.
		/// k = k0 + k1*lam1 + k2*lam2 + k3*lam3.
		BLST_HD static void glv_split_4(int64_t k_out[4], const limb_t k_in[4]);

		/// @brief Compiles a scalar into a constant-time execution schedule.
		/// Ideal for fixed/constant scalars reused across many bases.
		/// @param k_in Input scalar limbs.
		/// @return Compiled properties containing schedule and sign masks.
		BLST_HD static GLV4_Exp glv4_compile(const limb_t* k_in);

		/// @brief ret = a^scalar in Fp12, using GLV-4 optimization via precompiled schedule.
		/// Fastest available exponentiation path for constant scalars.
		/// base MUST be in the cyclotomic subgroup.
		BLST_HD static void cyclotomic_exp_glv4(Fp12& out, const Fp12& base, const GLV4_Exp& props);

		/// @brief ret = a^scalar in Fp12, using GLV-4 optimization (runtime split).
		/// Splits the scalar on the fly. ~3x faster than standard exp.
		/// base MUST be in the cyclotomic subgroup.
		BLST_HD static void cyclotomic_exp_glv4(Fp12& out, const Fp12& base, const limb_t* scalar);


		// --- Core Field Operations ---

		/// @brief ret = conjugate(a). In-place.
		/// In Fp12(w), conjugation negates the w component coefficients.
		BLST_HD static void conjugate(Fp12& a);

		/// @brief ret = 1/a in Fp12.
		BLST_HD static void inverse(Fp12& ret, const Fp12& a);

		/// @brief ret = a^(p^n) in Fp12.
		/// Frobenius endomorphism map. n must be <= 3 (usually).
		BLST_HD static void frobenius_map(Fp12& ret, const Fp12& a, size_t n);

		/// @brief Fast multiplication by sparse element (xy00z0).
		/// Used in pairing final exponentiation and Miller loop.
		BLST_HD static void mul_by_xy00z0(Fp12& ret, const Fp12& a, const Fp6& xy00z0);


		/// @brief Fast multiplication by sparse element (xy00z0).
		/// Used in pairing final exponentiation and Miller loop.
		/// Performed with scalar codegen (no AVX).
		BLST_HD static void mul_by_xy00z0_scalar(Fp12& ret, const Fp12& a, const Fp6& xy00z0);

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
        /// @brief Fast multiplication by sparse element (xy00z0).
        /// Used in pairing final exponentiation and Miller loop.
        /// Performed with AVX codegen.
        static void mul_by_xy00z0_avx(Fp12& ret, const Fp12& a, const Fp6& xy00z0);
#endif


		// --- Utilities ---

		/// @brief Returns 1 if a == b, 0 otherwise. Constant-time.
		BLST_HD static int equal(const Fp12& a, const Fp12& b);

		/// @brief Returns 1 if a == 1, 0 otherwise. Constant-time.
		BLST_HD static int is_one(const Fp12& a);

		/// @brief Convert Big-Endian bytes to Fp12 element (Montgomery form).
		BLST_HD static void from_bendian(Fp12& ret, const unsigned char serialized[12 * 48]);

		/// @brief Convert Fp12 element (Montgomery form) to Big-Endian bytes.
		BLST_HD static void to_bendian(unsigned char serialized[12 * 48], const Fp12& a);

		/// @brief Checks if element is in the Cyclotomic subgroup G_phi_12.
		/// Checks condition: a^(p^4 - p^2 + 1) == 1.
		BLST_HD static int in_group(const Fp12& f);

		/// @brief Faster check: is a in Cyclotomic subgroup?
		/// Validates if frob(a, 2) * a == frob(a, 4).
		BLST_HD static bool_t is_cyclotomic(const Fp12& f);


		// --- Internal / Advanced Helpers ---

		/// @brief ret = (ret * a)^(2^n) in Fp12.
		BLST_HD static void mul_n_sqr(Fp12& ret, const Fp12& a, size_t n);

		/// @brief Helper for Final Exp: ret = a^((z-1)/2).
		BLST_HD static void raise_to_z_div_by_2(Fp12& ret, const Fp12& a);

		/// @brief Helper for Final Exp: a = b^z.
		BLST_HD static void raise_to_z(Fp12& a, const Fp12& b);

		/// @brief Performs the Final Exponentiation step of pairings.
		/// Maps extension field element to unique representative in G_T.
		BLST_HD static void final_exp(Fp12& ret, const Fp12& f);

	};



}

#endif /* BLST_SCALAR_FP12_H */



