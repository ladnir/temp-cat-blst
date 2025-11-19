/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_FP_H
#define BLST_FP_H

#include "blst/common/vect.h"
 //#include "blst/fields/consts.h"
#include "blst/config.h"
#include "blst/common/defines.h"
#include "blst/backend/backend.h"
#include "blst/fields/recip-addchain.h"
#include <type_traits>

# include "sqrt-addchain.h"
//
//#ifdef __CUDA_ARCH__
//#include "blst/backend/cuda/add_mod_cuda.h"
//#include "blst/backend/cuda/mul_mont_cuda.h"
//#endif

namespace blst
{
	// The finite field Fp for BLS12-381 base field elements.
	// All values stored in memory are in Montgomery representation unless noted.
	//
	// Notation:
	//   - R = 2^384 mod p (Montgomery radix).
	//   - from_mont(x) returns the conventional integer representative in [0, p).
	//   - Unless specified, all operations are computed modulo p and return a
	//     reduced result in [0, p).
	//
	// Constant-time notes:
	//   - Unless stated otherwise, routines are intended to be constant-time w.r.t.
	//     secret data. The Lehmer almost-inverse (reciprocal) is typically variable-time;
	//     use flt_reciprocal() when strict constant-time is required.
	struct Fp
	{
		vec_384 data;

		//BLST_HD operator limb_t* () { return data; }
		//BLST_HD operator const limb_t* () const { return data; }
		BLST_HD BLST_FORCEINLINE auto&& operator[](u64 idx) { return data[idx]; }
		BLST_HD BLST_FORCEINLINE auto&& operator[](u64 idx) const { return data[idx]; }

		/** @name Canonical constants */
		///@{
		/// @brief The additive identity 0 (Montgomery form of 0).
		BLST_HD static const Fp& zero();

		/// @brief The multiplicative identity 1 (Montgomery form of 1, i.e., R mod p).
		BLST_HD static const Fp& one();

		/// @brief The modulus p as an Fp element (slight abuse; provided for convenience).
		/// @note Returned in Montgomery form (i.e., p * R mod p == 0).
		BLST_HD static constexpr const Fp& modulus();


		/// @brief The negative modulus p as an Fp element. i.e. -p mod 2^(384).
		BLST_HD static constexpr const Fp& neg_modulus();
		///@}



		BLST_HD static constexpr limb_t montgomery_constant();

		/** @name Basic arithmetic (mod p, Montgomery inputs/outputs) */
		///@{
		/// @brief ret = a + b mod p.
		BLST_HD static void add(Fp& ret, const Fp& a, const Fp& b);

		/// @brief ret = a - b mod p.
		BLST_HD static void sub(Fp& ret, const Fp& a, const Fp& b);

		/// @brief ret = a * b mod p.
		/// @note Inputs are Montgomery; product is reduced and remains Montgomery.
		BLST_HD static void mul(Fp& ret, const Fp& a, const Fp& b);

		/// @brief ret = a^x mod p.
		/// @param x Little-endian exponent bytes.
		/// @param x_bits Number of meaningful bits in @p x.
		/// @note Input/output are Montgomery.
		BLST_HD static void exp(Fp& ret, const Fp& a, const byte* x, size_t x_bits);

		/// @brief ret = a^2 mod p (Montgomery).
		BLST_HD static void sqr(Fp& ret, const Fp& a);
		///@}

		/** @name Montgomery transforms */
		///@{
		/// @brief Convert from Montgomery: ret = from_mont(a) in [0, p).
		/// @details If a = a' * R (mod p), returns a'.
		BLST_HD static void from_mont(Fp& ret, const Fp& a);

		/// @brief Convert to Montgomery: ret = a * R (mod p).
		/// @details If input a is already Montgomery, this multiplies by R again.
		BLST_HD static void to_mont(Fp& ret, const Fp& a);

		/// @brief Montgomery reduction: ret = a (mod p), where a is a double-precision UFp.
		BLST_HD static void redc(Fp& ret, const struct UFp& a);
		///@}

		/** @name Negation and conditional operations */
		///@{
		/// @brief Conditional negate: ret = flag ? (-a) : a  (constant-time).
		BLST_HD static void cneg(Fp& ret, const Fp& a, bool_t flag);

		/// @brief ret = -a mod p.
		BLST_HD static void neg(Fp& ret, const Fp& a);
		///@}

		/** @name Bit shifts and small multipliers
		 *  @note lshift/rshift view the canonical integer in [0,p) conceptually.
		 *  Implementations keep results reduced mod p unless documented.
		 */
		 ///@{
		 /// @brief ret = a * 2^count  (mod p).
		BLST_HD static void lshift(Fp& ret, const Fp& a, size_t count);

		/// @brief ret = floor( int(a) / 2^count ) (NOT mod p).
		/// @warning This is a raw right shift of the canonical representative,
		///          not a modular multiply by 2^{-count}.
		BLST_HD static void rshift(Fp& ret, const Fp& a, size_t count);

		/// @brief ret = a / 2  (mod p), i.e., multiply by (p+1)/2.
		BLST_HD static void div_by_2(Fp& ret, const Fp& a);

		/// @brief ret = a * 3  (mod p).
		BLST_HD static void mul_by_3(Fp& ret, const Fp& a);

		/// @brief ret = a * 8  (mod p).
		BLST_HD static void mul_by_8(Fp& ret, const Fp& a);
		///@}

		/** @name Inversion & square roots */
		///@{
		/// @brief out = a^{-1} mod p via Lehmer almost-inverse.
		BLST_HD static void reciprocal(Fp& out, const Fp& inp);

		/// @brief out = a^{-1} mod p via Fermat/addition chain .
		BLST_HD static void flt_reciprocal(Fp& out, const Fp& inp);

		/// @brief t0 = a^{(p-3)/4} (helper for sqrt/recip_sqrt when p ≡ 3 (mod 4)).
		BLST_HD static void recip_sqrt_3mod4(Fp& t0, const Fp& inp);

		/// @brief out = a^{-1/2} = sqrt(1/a) if it exists.
		/// @return 1 on success, 0 if a is non-square or zero.
		BLST_HD static bool_t recip_sqrt(Fp& out, const Fp& inp);

		/// @brief out = sqrt(a) if it exists (principal root per sgn0).
		/// @return 1 on success, 0 if a is non-square.
		BLST_HD static bool_t sqrt(Fp& out, const Fp& inp);

		/// @brief Return 1 iff a is a quadratic residue mod p (including 0), else 0.
		BLST_HD static bool_t is_sqrt(const Fp& inp);
		///@}

		/** @name Exponentiation helpers */
		///@{
		/// @brief out = a^{2^count} * b  (i.e., repeated squaring then multiply).
		BLST_HD static void sqr_n_mul(Fp& out, const Fp& a, size_t count, const Fp& b);
		///@}

		/** @name Comparisons & predicates */
		///@{
		/// @brief Constant-time equality on canonical representatives: from_mont(a) == from_mont(b).
		BLST_HD static bool_t equal(const Fp& a, const Fp& b);

		/// @brief Returns two bits encoding parity and sign (sgn0) of the canonical representative.
		/// @details Layout: bit0 = parity (a mod 2); bit1 = sgn0(a) per RFC 9380 / hash-to-curve.
		///          The returned value is (parity) | (sgn0 << 1).
		BLST_HD static limb_t sgn0_pty(const Fp& a);

		/// @brief Same as sgn0_pty(from_mont(a)) but operates on the Montgomery value directly.
		/// @details Returns (parity) | (sgn0 << 1) for the canonical integer represented by a.
		BLST_HD static limb_t sgn0_pty_mont(const Fp& a);

		/// @brief Return 1 iff a is strictly less than p (i.e., reduced), else 0.
		/// @note For Montgomery inputs, valid field elements still lie in [0, p),
		///       so this checks canonical range, not representation.
		BLST_HD static limb_t is_redc(const Fp& a);
		///@}

		/// @brief Number of machine limbs used to store an element (384 bits).
		BLST_HD static constexpr size_t num_limb() { return NLIMBS(384); }
	};

	static_assert(std::is_default_constructible_v<Fp>);


	// -1/P, this is the modular inverse of P modulo 2^64
	//   used in montgomery reduction as the "magic" number
#define BLST_p0  (limb_t)0x89f3fffcfffcfffd

#define ONE_MONT_P_LIMBS TO_LIMB_T(0x760900000002fffd), \
                   TO_LIMB_T(0xebf4000bc40c0002), \
                   TO_LIMB_T(0x5f48985753c758ba), \
                   TO_LIMB_T(0x77ce585370525745), \
                   TO_LIMB_T(0x5c071a97a256ec6d), \
                   TO_LIMB_T(0x15f65ec3fa80e493)


 /* (z-1)^2 * (z^4 - z^2 + 1)/3 + z */
#define BLST_P_LIMBS  \
    TO_LIMB_T(0xb9feffffffffaaab), TO_LIMB_T(0x1eabfffeb153ffff),\
    TO_LIMB_T(0x6730d2a0f6b0f624), TO_LIMB_T(0x64774b84f38512bf),\
    TO_LIMB_T(0x4b1ba7b6434bacd7), TO_LIMB_T(0x1a0111ea397fe69a)

#define BLST_NEG_P_LIMBS  \
    TO_LIMB_T(~0xb9feffffffffaaab), TO_LIMB_T(~0x1eabfffeb153ffff),\
    TO_LIMB_T(~0x6730d2a0f6b0f624), TO_LIMB_T(~0x64774b84f38512bf),\
    TO_LIMB_T(~0x4b1ba7b6434bacd7), TO_LIMB_T(~0x1a0111ea397fe69a+1)

#define BLST_RR_LIMBS \
    TO_LIMB_T(0xf4df1f341c341746), TO_LIMB_T(0x0a76e6a609d104f1),\
    TO_LIMB_T(0x8de5476c4c95b6d5), TO_LIMB_T(0x67eb88a9939d83c0),\
    TO_LIMB_T(0x9a793e85b519952d), TO_LIMB_T(0x11988fe592cae3aa)


#define BLST_RRRR_LIMBS \
    TO_LIMB_T(0xed48ac6bd94ca1e0), TO_LIMB_T(0x315f831e03a7adf8),\
    TO_LIMB_T(0x9a53352a615e29dd), TO_LIMB_T(0x34c04e5e921e1761),\
    TO_LIMB_T(0x2512d43565724728), TO_LIMB_T(0x0aa6346091755d4d)


	// one-in-Montgomery, (1<<384) % P
	// this is a union because each one is a prefix of the next
	// everything except the least significant Fp is zero
	typedef union { Fp p12[12]; Fp p2[2]; Fp p; } radix384;

	constexpr radix384 BLST_Rx_HOST = { /* (1<<384)%P, "radix", one-in-Montgomery */
	  { { ONE_MONT_P_LIMBS },
		{ 0 } }
	};

#if defined(BLST_CUDA) && defined(__CUDACC__)

	static __constant__ Fp BLST_P_DEV = { BLST_P_LIMBS };
	static __constant__ Fp BLST_NEG_P_DEV = { BLST_NEG_P_LIMBS };
	static __constant__ Fp BLST_RR_DEV = { BLST_RR_LIMBS };
	static __constant__ Fp BLST_RRRR_DEV = { BLST_RRRR_LIMBS };
	static __constant__ radix384 BLST_Rx_DEV = { {
		{ ONE_MONT_P_LIMBS },
		{ 0 }
	} };


#endif // __CUDA_ARCH__

	// returns R, the one-in-Montgomery form. 
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE const radix384& get_BLST_Rx(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_Rx_DEV;
#else
		return BLST_Rx_HOST;
#endif
	}

	constexpr Fp BLST_P_HOST = { BLST_P_LIMBS };
	constexpr Fp BLST_NEG_P_HOST = { BLST_NEG_P_LIMBS };
	constexpr Fp BLST_RR_HOST = { BLST_RR_LIMBS };   /* (1<<768)%P, "radix"^2, to-Montgomery   */
	constexpr Fp BLST_RRRR_HOST = { BLST_RRRR_LIMBS };


	// returns the prime modulus for BLS12-381. 
	// This is the field over which G1 is defined.
	// The order of the group G1 is a different prime, which can be obtained
	// via get_BLST_r().
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE constexpr const Fp& get_BLST_P(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_P_DEV;
#else
		return BLST_P_HOST;
#endif
	}


	// returns the negative prime modulus for BLS12-381. 
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE constexpr const Fp& get_BLST_negP(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_NEG_P_DEV;
#else
		return BLST_NEG_P_HOST;
#endif
	}



	// returns R^2, the to-Montgomery form value. Perform montgomery multiplication
	// with this value to convert a normal integer into montgomery form.
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE const Fp& get_BLST_RR(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_RR_DEV;
#else
		return BLST_RR_HOST;
#endif
	}

	// returns R^4, used in hash_to_field to convert 256-bit input into montgomery form.
	BLST_HD BLST_FORCEINLINE const Fp& get_BLST_RRRR(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_RRRR_DEV;
#else
		return BLST_RRRR_HOST;
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::add(Fp& ret, const Fp& a, const Fp& b)
	{
		BLST_ASSERT(Fp::is_redc(a));
		BLST_ASSERT(Fp::is_redc(b));
#ifdef __CUDA_ARCH__
		ptx_add_mod_u64<Fp>(ret.data, a.data, b.data);
#elif defined(BLST_DEFAULT_ASM) 
		add_mod_384(&ret, &a, &b, &get_BLST_P());
#else
		cc_add_mod_n<NLIMBS(384)>(ret.data, a.data, b.data, get_BLST_P().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::sub(Fp& ret, const Fp& a, const Fp& b)
	{
		BLST_ASSERT(Fp::is_redc(a));
		BLST_ASSERT(Fp::is_redc(b));
#if defined(__CUDA_ARCH__)
		ptx_sub_mod_u64<Fp>(ret.data, a.data, b.data);
#elif defined(BLST_DEFAULT_ASM) 
		sub_mod_384(&ret, &a, &b, &get_BLST_P());
#else
		cc_sub_mod_n<NLIMBS(384)>(ret.data, a.data, b.data, get_BLST_P().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::mul(Fp& ret, const Fp& a, const Fp& b)
	{
#ifdef __CUDA_ARCH__
		ptx_mul_mont_u64<Fp>(ret.data, a.data, b.data);
#elif defined(BLST_DEFAULT_ASM) 
		mulx_mont_384(&ret, &a, &b, &get_BLST_P(), BLST_p0);
#else
		cc_mul_mont_n<NLIMBS(384)>(ret.data, a.data, b.data, get_BLST_P().data, BLST_p0);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::sqr(Fp& ret, const Fp& a)
	{
		BLST_ASSERT(Fp::is_redc(a));

#if defined(__CUDA_ARCH__)
		ptx_sqr_mont_u64<Fp>(ret.data, a.data);
#elif defined(BLST_DEFAULT_ASM) 
		sqrx_mont_384(&ret, &a, &get_BLST_P(), BLST_p0);
#else
		cc_mul_mont_n<NLIMBS(384)>(ret.data, a.data, a.data, get_BLST_P().data, BLST_p0);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::from_mont(Fp& ret, const Fp& a)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		fromx_mont_384(&ret, &a, &get_BLST_P(), BLST_p0);
#else
		cc_from_mont_n<NLIMBS(384)>(ret.data, a.data, get_BLST_P().data, BLST_p0);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::to_mont(Fp& ret, const Fp& a)
	{
		mul(ret, a, get_BLST_RR());
	}

	BLST_HD BLST_FORCEINLINE void Fp::redc(Fp& ret, const struct UFp& a)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
#ifdef BLST_ADX
		redcx_mont_384(&ret, &a, &get_BLST_P(), BLST_p0);
#else
		redc_mont_384(&ret, &a, &get_BLST_P(), BLST_p0);
#endif
#else
		cc_redc_mont_n<NLIMBS(384)>(ret.data, a.data, get_BLST_P().data, BLST_p0);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::cneg(Fp& ret, const Fp& a, bool_t flag)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		cneg_mod_384(&ret, &a, flag, &get_BLST_P());
#else
		cc_cneg_mod_n<NLIMBS(384)>(ret.data, a.data, flag, get_BLST_P().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::neg(Fp& ret, const Fp& a)
	{
		cneg(ret, a, 1);
	}

	BLST_HD BLST_FORCEINLINE void Fp::lshift(Fp& ret, const Fp& a, size_t count)
	{
		assert(count);
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		lshift_mod_384(ret.data, a.data, count, get_BLST_P().data);
#else
		cc_lshift_mod_n<NLIMBS(384)>(ret.data, a.data, count, get_BLST_P().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::rshift(Fp& ret, const Fp& a, size_t count)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		rshift_mod_384(&ret, &a, count, &get_BLST_P());
#else
		cc_rshift_mod_n<NLIMBS(384)>(ret.data, a.data, count, get_BLST_P().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::div_by_2(Fp& ret, const Fp& a)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		div_by_2_mod_384(&ret, &a, &get_BLST_P());
#else
		cc_rshift_mod_n<NLIMBS(384)>(ret.data, a.data, 1, get_BLST_P().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::mul_by_3(Fp& ret, const Fp& a)
	{
		BLST_ASSERT(Fp::is_redc(a));

#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		mul_by_3_mod_384(&ret, &a, &get_BLST_P());
#else
		cc_mul_by_3_mod_n<NLIMBS(384)>(ret.data, a.data, get_BLST_P().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::mul_by_8(Fp& ret, const Fp& a)
	{
		BLST_ASSERT(Fp::is_redc(a));

#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		mul_by_8_mod_384(&ret, &a, &get_BLST_P());
#else
		cc_lshift_mod_n<NLIMBS(384)>(ret.data, a.data, 3, get_BLST_P().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::sqr_n_mul(Fp& ret, const Fp& a, size_t count, const Fp& b)
	{
		BLST_ASSERT(Fp::is_redc(a));

#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		sqrx_n_mul_mont_383(&ret, &a, count, &get_BLST_P(), BLST_p0, &b);
#else
		cc_sqr_n_mul_mont_384(ret.data, a.data, count, get_BLST_P().data, BLST_p0, b.data);
#endif
	}


	BLST_HD BLST_FORCEINLINE limb_t Fp::sgn0_pty_mont(const Fp& a)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return sgn0x_pty_mont_384(&a, &modulus(), BLST_p0);
#else
		Fp tmp;
		cc_from_mont_n<NLIMBS(384)>(tmp.data, a.data, modulus().data, BLST_p0);
		return cc_sgn0_pty_mod_n<NLIMBS(384)>(tmp.data, modulus().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE limb_t Fp::sgn0_pty(const Fp& a)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return sgn0_pty_mod_384(a.data, modulus().data);
#else
		return cc_sgn0_pty_mod_n<NLIMBS(384)>(a.data, modulus().data);
#endif
	}


	BLST_HD BLST_FORCEINLINE void ct_inverse_mod_v384(UFp& ret, const Fp& inp, const Fp& mod,
		const Fp& modx)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		ctx_inverse_mod_384(&ret, &inp, &mod, &modx);
#else
		//cc_ct_inverse_mod_384(ret, inp, mod, modx);
		cc_ct_inverse_mod_n<NLIMBS(384)>(ret.data, inp.data, mod.data, modx.data);
#endif
	}


	BLST_HD BLST_FORCEINLINE bool_t Fp::is_sqrt(const Fp& inp)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return ct_is_square_mod_384(&inp, &modulus());
#else
		return cc_ct_is_sqr_mod_n<NLIMBS(384)>(inp.data, modulus().data);
#endif
	}

	/*
	 * |out| = |inp|^|pow|, small footprint, public exponent
	 */
	BLST_HD BLST_INLINE void Fp::exp(Fp& out, const Fp& inp, const byte* pow,
		size_t pow_bits)
	{
		Fp ret;

		// requies MSB of pow is set
		BLST_ASSERT(pow_bits && ((pow[(pow_bits - 1) >> 3] >> ((pow_bits - 1) & 7)) & 1));

		vec_copy(ret.data, inp.data, sizeof(ret));  /* ret = inp^1 */
		--pow_bits; /* most significant bit is set, skip over */
		while (pow_bits--) {
			Fp::sqr(ret, ret);
			if (is_bit_set(pow, pow_bits))
				Fp::mul(ret, ret, inp);
		}
		vec_copy(out.data, ret.data, sizeof(ret));  /* out = ret */
	}


	BLST_HD BLST_INLINE void Fp::recip_sqrt_3mod4(Fp& out, const Fp& inp)
	{
#if defined(__OPTIMIZE_SIZE__) || defined(__CUDA_ARCH__)
		const byte BLS_12_381_P_minus_3_div_4[] = {
			TO_BYTES(0xee7fbfffffffeaaa), TO_BYTES(0x07aaffffac54ffff),
			TO_BYTES(0xd9cc34a83dac3d89), TO_BYTES(0xd91dd2e13ce144af),
			TO_BYTES(0x92c6e9ed90d2eb35), TO_BYTES(0x0680447a8e5ff9a6)
		};
		Fp::exp(out, inp, BLS_12_381_P_minus_3_div_4, 379);
#else
		add_chain_recip_sqrt<Fp>(out, inp);
#endif
	}

	BLST_HD BLST_FORCEINLINE bool_t Fp::recip_sqrt(Fp& out, const Fp& inp)
	{
		Fp t0, t1;
		bool_t ret;

		recip_sqrt_3mod4(t0, inp);

		Fp::mul(t1, t0, inp);
		Fp::sqr(t1, t1);
		ret = vec_is_equal(t1.data, inp.data, sizeof(t1));
		vec_copy(out.data, t0.data, sizeof(t0));

		return ret;
	}

	BLST_HD BLST_FORCEINLINE bool_t Fp::sqrt(Fp& out, const Fp& inp)
	{
		Fp t0, t1;
		bool_t ret;

		recip_sqrt_3mod4(t0, inp);

		Fp::mul(t0, t0, inp);
		Fp::sqr(t1, t0);
		ret = vec_is_equal(t1.data, inp.data, sizeof(t1));
		vec_copy(out.data, t0.data, sizeof(t0));

		return ret;
	}


	BLST_HD BLST_FORCEINLINE bool_t Fp::equal(const Fp& a, const Fp& b)
	{
		Fp aa, bb;
		from_mont(aa, a);
		from_mont(bb, b);
		return vec_is_equal(aa.data, bb.data, sizeof(Fp));
	}




	BLST_HD BLST_FORCEINLINE void Fp::flt_reciprocal(Fp& out, const Fp& inp)
	{
#if defined(__OPTIMIZE_SIZE__) || defined(__CUDA_ARCH__)
	/*
	 * 608 multiplications for scalar inversion modulo BLS12-381 prime, 32%
	 * more than corresponding optimal addition-chain, plus mispredicted
	 * branch penalties on top of that... The addition chain below was
	 * measured to be >50% faster.
	 */
		const byte BLS12_381_P_minus_2[] = {
			TO_BYTES(0xb9feffffffffaaa9), TO_BYTES(0x1eabfffeb153ffff),
			TO_BYTES(0x6730d2a0f6b0f624), TO_BYTES(0x64774b84f38512bf),
			TO_BYTES(0x4b1ba7b6434bacd7), TO_BYTES(0x1a0111ea397fe69a)
		};

		Fp::exp(out, inp, BLS12_381_P_minus_2, 381);
#else
		// use a unrolled addition chain
		RECIPROCAL_MOD_BLS12_381_P<Fp>(out, inp);
#endif
	}


	BLST_HD BLST_FORCEINLINE void Fp::reciprocal(Fp& out, const Fp& inp) {

		const Fp Px8 = {    /* left-aligned value of the modulus */
			TO_LIMB_T(0xcff7fffffffd5558), TO_LIMB_T(0xf55ffff58a9ffffd),
			TO_LIMB_T(0x39869507b587b120), TO_LIMB_T(0x23ba5c279c2895fb),
			TO_LIMB_T(0x58dd3db21a5d66bb), TO_LIMB_T(0xd0088f51cbff34d2)
		};
		union { UFp x; Fp r[2]; } temp;

		// returns temp.x * inp = 1 (mod mod). not return a Montgomery value; 
		// that’s why we then does redc + multiply by RR to move into Mont.
		ct_inverse_mod_v384(temp.x, inp, get_BLST_P(), Px8);

		// map back to montgomery
		Fp::redc(temp.r[0], temp.x);
		Fp::mul(temp.r[0], temp.r[0], get_BLST_RR());

#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
		/* sign goes straight to flt_reciprocal */
		Fp::mul(temp.r[1], temp.r[0], inp);
		if (vec_is_equal(&temp.r[1], &get_BLST_Rx().p, sizeof(Fp)) |
			vec_is_zero(&temp.r[1], sizeof(Fp)))
			vec_copy(&out, &temp.r[0], sizeof(Fp));
		else
		{
			// fallback, should be rare
			Fp::flt_reciprocal(out, inp);
		}
#else
		vec_copy(&out, &temp.r[0], sizeof(Fp));
#endif
	}

	BLST_HD BLST_FORCEINLINE limb_t Fp::is_redc(const Fp& a)
	{
		return cc_is_redc_mod_n<NLIMBS(384)>((const byte*)a.data, modulus().data);
	}


	BLST_HD BLST_FORCEINLINE const Fp& Fp::zero()
	{
		return get_BLST_Rx().p2[1];
	}
	BLST_HD BLST_FORCEINLINE const Fp& Fp::one()
	{
		return get_BLST_Rx().p;
	}
	BLST_HD BLST_FORCEINLINE constexpr const Fp& Fp::modulus()
	{
		return get_BLST_P();
	}


	BLST_HD BLST_FORCEINLINE constexpr const Fp& Fp::neg_modulus()
	{
		return get_BLST_negP();
	}



	BLST_HD BLST_FORCEINLINE constexpr limb_t Fp::montgomery_constant()
	{
		return BLST_p0;
	}

}

#endif /* BLST_FP_H */
