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
#include "blst/fields/exp.h"


namespace blst
{


	struct Fp
	{
		vec_384 data;

		BLST_HD operator limb* () { return data; }
		BLST_HD operator const limb* () const { return data; }

		BLST_HD static const Fp& zero();
		BLST_HD static const Fp& one();
		BLST_HD static const Fp& modulus();

		BLST_HD static void add(Fp& ret, const Fp& a, const Fp& b);

		BLST_HD static void sub(Fp& ret, const Fp& a, const Fp& b);

		BLST_HD static void mul(Fp& ret, const Fp& a, const Fp& b);

		BLST_HD static void sqr(Fp& ret, const Fp& a);

		BLST_HD static void from_mont(Fp& ret, const Fp& a);

		BLST_HD static void to_mont(Fp& ret, const Fp& a);

		BLST_HD static void redc(Fp& ret, const struct UFp& a);

		BLST_HD static void cneg(Fp& ret, const Fp& a, bool_t flag);

		BLST_HD static void neg(Fp& ret, const Fp& a);

		BLST_HD static void lshift(Fp& ret, const Fp& a, size_t count);

		BLST_HD static void rshift(Fp& ret, const Fp& a, size_t count);

		BLST_HD static void div_by_2(Fp& ret, const Fp& a);

		BLST_HD static void mul_by_3(Fp& ret, const Fp& a);

		BLST_HD static void mul_by_8(Fp& ret, const Fp& a);

		BLST_HD static void reciprocal(Fp& out, const Fp& inp);

		BLST_HD static void flt_reciprocal(Fp& out, const Fp& inp);
		
		BLST_HD static bool_t recip_sqrt(Fp& out, const Fp& inp);
		
		BLST_HD static bool_t sqrt(Fp& out, const Fp& inp);
		
		BLST_HD static void sqr_n_mul(Fp& out, const Fp& a, size_t count, const Fp& b);


	};

	static_assert(std::is_default_constructible_v<Fp>);


	// -1/P, this is the modular inverse of P modulo 2^64
	//   used in montgomery reduction as the "magic" number
#define BLST_p0  (limb)0x89f3fffcfffcfffd

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

	extern const radix384 BLST_Rx_HOST; /* (1<<384)%P, "radix", one-in-Montgomery */

#if defined(BLST_CUDA) && defined(__CUDACC__)

	static __constant__ Fp BLST_P_DEV = { BLST_P_LIMBS };
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


//#define ZERO_384 (get_BLST_Rx().p2[1])



	extern const Fp BLST_P_HOST;
	extern const Fp BLST_RR_HOST;   /* (1<<768)%P, "radix"^2, to-Montgomery   */
	extern const Fp BLST_RRRR_HOST; /* (1<<1536)%P, "radix"^4, hash_to_field */


	// returns the prime modulus for BLS12-381. 
	// This is the field over which G1 is defined.
	// The order of the group G1 is a different prime, which can be obtained
	// via get_BLST_r().
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE const Fp& get_BLST_P(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_P_DEV;
#else
		return BLST_P_HOST;
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









	//BLST_HD BLST_INLINE void cc_mul_mont_384(
	//	Fp& ret, const Fp& a, const Fp& b, const Fp& p, limb n0)
	//{
	//	//cc_mul_mont_n<NLIMBS(384)>(ret.data, a.data, b.data, p.data, n0);
	//}
	//BLST_HD BLST_INLINE void cc_sqr_mont_384(
	//	Fp& ret, const Fp& a, const Fp& p, limb n0)
	//{
	//	//cc_mul_mont_n<NLIMBS(384)>(ret.data, a.data, a.data, p.data, n0);
	//}

	//BLST_HD BLST_INLINE void cc_add_mod_384(Fp& ret, const Fp& a,
	//	const Fp& b, const Fp& p)
	//{
	//	//cc_add_mod_n<NLIMBS(384)>(ret.data, a.data, b.data, p.data);
	//}


	//BLST_HD BLST_INLINE void cc_sub_mod_384(Fp& ret, const Fp& a,
	//	const Fp& b, const Fp& p)
	//{
	//	cc_sub_mod_n<NLIMBS(384)>(ret.data, a.data, b.data, p.data);
	//}

	//BLST_HD BLST_INLINE void cc_mul_by_3_mod_384(Fp& ret, const Fp& a,
	//	const Fp& p)
	//{
	//	//cc_mul_by_3_mod_n<NLIMBS(384)>(ret.data, a.data, p.data);
	//}


	BLST_HD BLST_INLINE void cc_lshift_mod_384(Fp& ret, const Fp& a, size_t count,
		const Fp& p)
	{
		cc_lshift_mod_n<NLIMBS(384)>(ret.data, a.data, count, p.data);
	}


	//BLST_HD BLST_INLINE void cc_cneg_mod_384(Fp& ret, const Fp& a, bool_t flag,
	//	const Fp& p)
	//{
	//	//cc_cneg_mod_n<NLIMBS(384)>(ret.data, a.data, flag, p.data);
	//}

	//BLST_HD BLST_INLINE void cc_from_mont_384(Fp& ret, const Fp& a,
	//	const Fp& p, limb n0)
	//{
	//	//cc_from_mont_n<NLIMBS(384)>(ret.data, a.data, p.data, n0);
	//}


	//BLST_HD BLST_INLINE void cc_redc_mont_384(Fp& ret, const UFp& a,
	//	const Fp& p, limb n0)
	//{
	//	//cc_redc_mont_n<NLIMBS(384)>(ret.data, a.data, p.data, n0);
	//}


	//BLST_HD BLST_INLINE void cc_rshift_mod_384(Fp& ret, const Fp& a, size_t count,
	//	const Fp& p)
	//{
	//	//cc_rshift_mod_n(ret.data, a.data, count, p.data, NLIMBS(384));
	//}

	//BLST_HD BLST_INLINE void cc_div_by_2_mod_384(Fp& ret, const Fp& a,
	//	const Fp& p)
	//{
	//	//cc_rshift_mod_n(ret.data, a.data, 1, p.data, NLIMBS(384));
	//}


	/* @brief Compute sign and parity of a in Montgomery form mod p. see cc_sgn0_pty_mod_n */
	BLST_HD BLST_INLINE limb cc_sgn0_pty_mod_384(const Fp& a, const Fp& p)
	{
		return cc_sgn0_pty_mod_n<NLIMBS(384)>(a.data, p.data);
	}

	/* @brief Compute sign and parity of a in Montgomery form mod p. see cc_sgn0_pty_mod_n
	 */
	BLST_HD BLST_INLINE limb cc_sgn0_pty_mont_384(const Fp& a, const Fp& p, limb n0)
	{
		Fp tmp;

		cc_from_mont_n<NLIMBS(384)>(tmp.data, a.data, p.data, n0);

		return cc_sgn0_pty_mod_n<NLIMBS(384)>(tmp.data, p.data);
	}


	/**
	 * @brief Performs a left shift modulo operation on a 384-bit integer.
	 *
	 * Computes `(a << n) mod mod` iteratively, ensuring the result remains within the range `[0, mod)`.
	 *
	 * @param[out] ret The output 384-bit integer.
	 * @param[in]  a   The input 384-bit integer.
	 * @param[in]  n   The number of left shifts to perform.
	 * @param[in]  mod The modulus.
	 *
	 * @note Assumes `n > 0` and `a` and `mod` are valid 384-bit integers.
	 */
	BLST_HD BLST_INLINE void lshift_mod_v384(Fp& ret, const Fp& a, size_t n,
		const Fp& mod)
	{
		auto ptr = &a;
		while (n--)
		{
			//cc_add_mod_384(ret, *ptr, *ptr, mod);
			Fp::add(ret, *ptr, *ptr);
			ptr = &ret;
		}
	}

	/* @brief Multiply by 8 modulo mod */
	BLST_HD BLST_INLINE void cc_mul_by_8_mod_384(Fp& ret, const Fp& a, const Fp& mod)
	{
		lshift_mod_v384(ret, a, 3, mod);
		//Fp::lshift(ret, a, 3);
	}



	//BLST_HD BLST_INLINE void cc_ct_inverse_mod_384(UFp& ret, const Fp& inp,
	//	const Fp& mod, const Fp& modx)
	//{
	//	//cc_ct_inverse_mod_n<NLIMBS(384)>(ret, inp, mod, modx);
	//}


	//BLST_HD BLST_INLINE bool_t cc_ct_is_square_mod_384(const Fp& inp,
	//	const Fp& mod)
	//{
	//	//return cc_ct_is_sqr_mod_n<NLIMBS(384)>(inp, mod);
	//}
























	BLST_HD BLST_FORCEINLINE void Fp::add(Fp& ret, const Fp& a, const Fp& b) {
		//add_mod_v384(ret, a, b, get_BLST_P());

#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		add_mod_384(&ret, &a, &b, &get_BLST_P());
#else
		//cc_add_mod_384(ret, a, b, get_BLST_P());
		cc_add_mod_n<NLIMBS(384)>(ret.data, a.data, b.data, get_BLST_P().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::sub(Fp& ret, const Fp& a, const Fp& b) {
		//sub_mod_v384(ret, a, b, get_BLST_P());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		sub_mod_384(&ret, &a, &b, &get_BLST_P());
#else
		cc_sub_mod_n<NLIMBS(384)>(ret, a, b, get_BLST_P());
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::mul(Fp& ret, const Fp& a, const Fp& b) {
		//mul_mont_v384(ret, a, b, get_BLST_P(), BLST_p0);
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		mulx_mont_384(&ret, &a, &b, &get_BLST_P(), BLST_p0);
#else
		//cc_mul_mont_384(ret, a, b, get_BLST_P(), BLST_p0);
		cc_mul_mont_n<NLIMBS(384)>(ret.data, a.data, b.data, get_BLST_P().data, BLST_p0);

#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::sqr(Fp& ret, const Fp& a) {
		//sqr_mont_v384(ret, a, get_BLST_P(), BLST_p0);
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		sqrx_mont_384(&ret, &a, &get_BLST_P(), BLST_p0);
#else
		//cc_sqr_mont_384(ret, a, get_BLST_P(), BLST_p0);
		cc_mul_mont_n<NLIMBS(384)>(ret.data, a.data, a.data, get_BLST_P().data, BLST_p0);

#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::from_mont(Fp& ret, const Fp& a) {
		//from_mont_v384(ret, a, get_BLST_P(), BLST_p0);
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		fromx_mont_384(&ret, &a, &get_BLST_P(), BLST_p0);
#else
		//cc_from_mont_384(ret, a, get_BLST_P(), BLST_p0);
		cc_from_mont_n<NLIMBS(384)>(ret.data, a.data, get_BLST_P().data, BLST_p0);

#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::to_mont(Fp& ret, const Fp& a) {
		//mul_mont_v384(ret, a, get_BLST_RR(), get_BLST_P(), BLST_p0);
		mul(ret, a, get_BLST_RR());
	}

	BLST_HD BLST_FORCEINLINE void Fp::redc(Fp& ret, const struct UFp& a) {
		//redc_mont_v384(ret, a, get_BLST_P(), BLST_p0);
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
#ifdef BLST_ADX
		redcx_mont_384(&ret, &a, &get_BLST_P(), BLST_p0);
#else
		redc_mont_384(&ret, &a, &get_BLST_P(), BLST_p0);
#endif
#else
		//cc_redc_mont_384(ret, a, get_BLST_P(), BLST_p0);
		cc_redc_mont_n<NLIMBS(384)>(ret.data, a.data, get_BLST_P().data, BLST_p0);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::cneg(Fp& ret, const Fp& a, bool_t flag) {
		//cneg_mod_v384(ret, a, flag, get_BLST_P());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		cneg_mod_384(&ret, &a, flag, &get_BLST_P());
#else
		//cc_cneg_mod_384(ret, a, flag, get_BLST_P());
		cc_cneg_mod_n<NLIMBS(384)>(ret.data, a.data, flag, get_BLST_P().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::neg(Fp& ret, const Fp& a) {
		cneg(ret, a, 1);
	}

	BLST_HD BLST_FORCEINLINE void Fp::lshift(Fp& ret, const Fp& a, size_t n) {
		lshift_mod_v384(ret, a, n, get_BLST_P());
	}

	BLST_HD BLST_FORCEINLINE void Fp::rshift(Fp& ret, const Fp& a, size_t count) {
		//rshift_mod_v384(ret, a, count, get_BLST_P());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		rshift_mod_384(&ret, &a, count, &get_BLST_P());
#else
		//cc_rshift_mod_384(ret, a, count, get_BLST_P());
		cc_rshift_mod_n<NLIMBS(384)>(ret.data, a.data, count, get_BLST_P().data);

#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::div_by_2(Fp& ret, const Fp& a) {
		//div_by_2_mod_v384(ret, a, get_BLST_P());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		div_by_2_mod_384(&ret, &a, &get_BLST_P());
#else
		//cc_div_by_2_mod_384(ret, a, p);
		cc_rshift_mod_n<NLIMBS(384)>(ret.data, a.data, 1, get_BLST_P().data);

#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::mul_by_3(Fp& ret, const Fp& a) {
		//mul_by_3_mod_v384(ret, a, get_BLST_P());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		mul_by_3_mod_384(&ret, &a, &get_BLST_P());
#else
		//cc_mul_by_3_mod_384(ret, a, p);
		cc_mul_by_3_mod_n<NLIMBS(384)>(ret.data, a.data, get_BLST_P().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::mul_by_8(Fp& ret, const Fp& a) {
		//mul_by_8_mod_v384(ret, a, get_BLST_P());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		mul_by_8_mod_384(&ret, &a, &get_BLST_P());
#else
		cc_mul_by_8_mod_384(ret, a, get_BLST_P());
		
#endif
	}



	//// 
	//BLST_HD BLST_FORCEINLINE void mul_mont_v384(Fp& ret, const Fp& a, const Fp& b,
	//	const Fp& p, limb n0)
	//{

	//}

	//BLST_HD BLST_FORCEINLINE void sqr_mont_v384(Fp& ret, const Fp& a, const Fp& p, limb n0)
	//{

	//}

	//BLST_HD BLST_FORCEINLINE void sqr_n_mul_mont_v384(Fp& ret, const Fp& a, size_t count,
	//	const Fp& p, limb n0, const Fp& b)
	//{
	//#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	//	sqrx_n_mul_mont_384(ret, a, count, get_BLST_P(), BLST_p0, b);
	//#else
	//	cc_sqr_n_mul_mont_384(ret, a, count, get_BLST_P(), BLST_p0, b);
	//#endif
	//}

	/**
	 * @brief Squares a 384-bit integer multiple times and multiplies the result in Montgomery form.
	 *
	 * This function computes `(a^(2^count) * b * R^(-1)) mod p` in Montgomery form, where `R` is the Montgomery constant.
	 * The input `a` is squared `count` times, and the result is multiplied by `b` modulo `p`.
	 *
	 * @param[out] ret   The output 384-bit integer in Montgomery form.
	 * @param[in]  a     The input 384-bit integer in Montgomery form.
	 * @param[in]  count The number of squaring operations to perform.
	 * @param[in]  p     The modulus.
	 * @param[in]  n0    The Montgomery constant.
	 * @param[in]  b     The multiplier in Montgomery form.
	 *
	 * @note Assumes `count > 0` and `a`, `b`, and `p` are valid 384-bit integers.
	 */
	BLST_HD BLST_INLINE void cc_sqr_n_mul_mont_384(Fp& ret, const Fp& a, size_t count,
		const Fp& p, limb n0, const Fp& b)
	{
		BLST_ASSUME(count != 0);
		auto ptr = &a;

		while (count--) {
			cc_mul_mont_nonred_n<NLIMBS(384)>(ret, *ptr, *ptr, p, n0);
			ptr = &ret;
		}
		cc_mul_mont_n<NLIMBS(384)>(ret, ret, b, p, n0);
	}


	/*
	 * "383"-bit variant omits full reductions at the ends of squarings,
	 * which results in up to ~15% improvement. [One can improve further
	 * by omitting full reductions even after multiplications and
	 * performing final reduction at the very end of the chain.]
	 *
	 * see sqr_n_Fp::mul
	 */
	BLST_HD BLST_FORCEINLINE void sqr_n_mul_mont_v383(Fp& ret, const Fp& a, size_t count,
		const Fp& p, limb n0, const Fp& b)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		sqrx_n_mul_mont_383(&ret, &a, count, &get_BLST_P(), BLST_p0, &b);
#else
		cc_sqr_n_mul_mont_384(ret, a, count, get_BLST_P(), BLST_p0, b);
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::sqr_n_mul(Fp& ret, const Fp& a, size_t count, const Fp& b)
	{
		sqr_n_mul_mont_v383(ret, a, count, get_BLST_P(), BLST_p0, b);
	}

	//void mul_384(UFp& ret, const Fp& a, const Fp& b)
	//{
	//#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
	//	mulx_384(ret, a, b);
	//#else
	//	cc_mul_384(ret, a, b);
	//#endif
	//}
	//void sqr_384(UFp& ret, const Fp& a);

	//BLST_HD BLST_FORCEINLINE void redc_mont_v384(Fp& ret, const UFp& a, const Fp& p, limb n0)
	//{

	//}

	//BLST_HD BLST_FORCEINLINE void from_mont_v384(Fp& ret, const Fp& a, const Fp& p, limb n0)
	//{

	//}


	BLST_HD BLST_FORCEINLINE limb sgn0_pty_mont_v384(const Fp& a, const Fp& p, limb n0)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return sgn0x_pty_mont_384(&a, &p, n0);
#else
		return cc_sgn0_pty_mont_384(a, p, n0);
#endif
	}

	BLST_HD BLST_FORCEINLINE limb sgn0_pty_mod_v384(const Fp& a, const Fp& p)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return sgn0_pty_mod_384(&a, &p);
#else
		return cc_sgn0_pty_mod_384(a, p);
#endif
	}


	//BLST_HD BLST_FORCEINLINE void add_mod_v384(Fp& ret, const Fp& a, const Fp& b, const Fp& p)
	//{
	//}

	//BLST_HD BLST_FORCEINLINE void sub_mod_v384(Fp& ret, const Fp& a, const Fp& b, const Fp& p)
	//{

	//}
	//BLST_HD BLST_FORCEINLINE void mul_by_8_mod_v384(Fp& ret, const Fp& a, const Fp& p)
	//{


	//}
	//BLST_HD BLST_FORCEINLINE void mul_by_3_mod_v384(Fp& ret, const Fp& a, const Fp& p)
	//{

	//}

	//BLST_HD BLST_FORCEINLINE void cneg_mod_v384(Fp& ret, const Fp& a, bool_t flag, const Fp& p)
	//{

	//}

	//void lshift_mod_v384(Fp& ret, const Fp& a, size_t count, const Fp& p);

	//BLST_HD BLST_FORCEINLINE void rshift_mod_v384(Fp& ret, const Fp& a, size_t count, const Fp& p)
	//{

	//}
	//BLST_HD BLST_FORCEINLINE void div_by_2_mod_v384(Fp& ret, const Fp& a, const Fp& p)
	//{

	//}


	BLST_HD BLST_FORCEINLINE void ct_inverse_mod_v384(UFp& ret, const Fp& inp, const Fp& mod,
		const Fp& modx)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		ctx_inverse_mod_384(&ret, &inp, &mod, &modx);
#else
		//cc_ct_inverse_mod_384(ret, inp, mod, modx);
		cc_ct_inverse_mod_n<NLIMBS(384)>(ret, inp, mod, modx);
#endif
	}


	BLST_HD BLST_FORCEINLINE bool_t ct_is_square_mod_v384(const Fp& inp, const Fp& mod)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return ct_is_square_mod_384(&inp, &mod);
#else
		//return cc_ct_is_square_mod_384(inp, mod);
		return cc_ct_is_sqr_mod_n<NLIMBS(384)>(inp, mod);

#endif
	}



//
//
//	//BLST_EXTERN_C_BEGIN
//
//	/*
//	 * BLS12-381-specific Fp shortcuts to assembly.
//	 */
//	BLST_HD BLST_FORCEINLINE void Fp::add(Fp& ret, const Fp& a, const Fp& b)
//	{
//		add_mod_v384(ret, a, b, get_BLST_P());
//	}
//
//	BLST_HD BLST_FORCEINLINE void Fp::sub(Fp& ret, const Fp& a, const Fp& b)
//	{
//		sub_mod_v384(ret, a, b, get_BLST_P());
//	}
//
//	BLST_HD BLST_FORCEINLINE void mul_by_3_fp(Fp& ret, const Fp& a)
//	{
//		mul_by_3_mod_v384(ret, a, get_BLST_P());
//	}
//
//	BLST_HD BLST_FORCEINLINE void mul_by_8_fp(Fp& ret, const Fp& a)
//	{
//		mul_by_8_mod_v384(ret, a, get_BLST_P());
//	}
//
//	BLST_HD BLST_FORCEINLINE void lshift_fp(Fp& ret, const Fp& a, size_t count)
//	{
//		lshift_mod_v384(ret, a, count, get_BLST_P());
//	}
//
//	BLST_HD BLST_FORCEINLINE void rshift_fp(Fp& ret, const Fp& a, size_t count)
//	{
//		rshift_mod_v384(ret, a, count, get_BLST_P());
//	}
//
//	BLST_HD BLST_FORCEINLINE void Fp::div_by_2(Fp& ret, const Fp& a)
//	{
//		div_by_2_mod_v384(ret, a, get_BLST_P());
//	}
//
//	BLST_HD BLST_FORCEINLINE void Fp::mul(Fp& ret, const Fp& a, const Fp& b)
//	{
//		mul_mont_v384(ret, a, b, get_BLST_P(), BLST_p0);
//	}
//
//	BLST_HD BLST_FORCEINLINE void Fp::sqr(Fp& ret, const Fp& a)
//	{
//		sqr_mont_v384(ret, a, get_BLST_P(), BLST_p0);
//	}
//
//	BLST_HD BLST_FORCEINLINE void cFp::neg(Fp& ret, const Fp& a, bool_t flag)
//	{
//		cneg_mod_v384(ret, a, flag, get_BLST_P());
//	}
//
//	BLST_HD BLST_FORCEINLINE void from_fp(Fp& ret, const Fp& a)
//	{
//		from_mont_v384(ret, a, get_BLST_P(), BLST_p0);
//	}
//
//	BLST_HD BLST_FORCEINLINE void to_fp(Fp& ret, const Fp& a)
//	{
//		Fp::mul(ret, a, get_BLST_RR());
//	}
//
//	BLST_HD BLST_FORCEINLINE void redc_fp(Fp& ret, const UFp& a)
//	{
//		redc_mont_v384(ret, a, get_BLST_P(), BLST_p0);
//	}
//
//#define vec_load_global vec_copy
//
//	BLST_HD void reciprocal_fp(Fp& out, const Fp& inp);
//	BLST_HD void flt_reciprocal_fp(Fp& out, const Fp& inp);
//	BLST_HD bool_t recip_sqrt_fp(Fp& out, const Fp& inp);
//	BLST_HD bool_t sqrt_fp(Fp& out, const Fp& inp);
//
//#define Fp::neg(r,a) cFp::neg((r),(a),1)
//	//#define Fp::neg2(r,a) Fp2::cneg((r),(a),1)














	//BLST_EXTERN_C_BEGIN
	/*
	 * |out| = |inp|^|pow|, small footprint, public exponent
	 */
	BLST_HD BLST_INLINE void exp_mont_384(Fp& out, const Fp& inp, const byte* pow,
		size_t pow_bits/*, const Fp& p, limb n0*/)
	{
		Fp ret;

		vec_copy(ret, inp, sizeof(ret));  /* ret = inp^1 */
		--pow_bits; /* most significant bit is set, skip over */
		while (pow_bits--) {
			Fp::sqr(ret, ret/*, p, n0*/);
			if (is_bit_set(pow, pow_bits))
				Fp::mul(ret, ret, inp/*, p, n0*/);
		}
		vec_copy(out, ret, sizeof(ret));  /* out = ret */
	}




#if defined(__OPTIMIZE_SIZE__)  || defined(__CUDA_ARCH__)
	BLST_HD BLST_FORCEINLINE void recip_sqrt_fp_3mod4(Fp& out, const Fp& inp)
	{
		const byte BLS_12_381_P_minus_3_div_4[] = {
			TO_BYTES(0xee7fbfffffffeaaa), TO_BYTES(0x07aaffffac54ffff),
			TO_BYTES(0xd9cc34a83dac3d89), TO_BYTES(0xd91dd2e13ce144af),
			TO_BYTES(0x92c6e9ed90d2eb35), TO_BYTES(0x0680447a8e5ff9a6)
		};

		exp_mont_384(out, inp, BLS_12_381_P_minus_3_div_4, 379);
	}
#else
# if 1
	/*
	 * "383"-bit variant omits full reductions at the ends of squarings,
	 * which results in up to ~15% improvement. [One can improve further
	 * by omitting full reductions even after multiplications and
	 * performing final reduction at the very end of the chain.]
	 */
	BLST_HD BLST_FORCEINLINE void sqr_n_mul_fp(Fp& out, const Fp& a, size_t count,
		const Fp& b)
	{
		sqr_n_mul_mont_v383(out, a, count, get_BLST_P(), BLST_p0, b);
	}
# else
	BLST_HD void sqr_n_Fp::mul(Fp& out, const Fp& a, size_t count,
		const Fp& b)
	{
		while (count--) {
			Fp::sqr(out, a);
			a = out;
		}
		Fp::mul(out, out, b);
	}
# endif

# define sqr(ret,a)		Fp::sqr(ret,a)
# define mul(ret,a,b)		Fp::mul(ret,a,b)
# define sqr_n_mul(ret,a,n,b)	Fp::sqr_n_mul(ret,a,n,b)

# include "sqrt-addchain.h"
	BLST_HD BLST_INLINE void recip_sqrt_fp_3mod4(Fp& out, const Fp& inp)
	{
		RECIP_SQRT_MOD_BLS12_381_P(out, inp, Fp);
	}


# undef RECIP_SQRT_MOD_BLS12_381_P
# undef sqr_n_mul
# undef sqr
# undef mul
#endif

	BLST_HD BLST_INLINE bool_t recip_sqrt_fp(Fp& out, const Fp& inp)
	{
		Fp t0, t1;
		bool_t ret;

		recip_sqrt_fp_3mod4(t0, inp);

		Fp::mul(t1, t0, inp);
		Fp::sqr(t1, t1);
		ret = vec_is_equal(t1, inp, sizeof(t1));
		vec_copy(out, t0, sizeof(t0));

		return ret;
	}

	BLST_HD BLST_INLINE bool_t sqrt_fp(Fp& out, const Fp& inp)
	{
		Fp t0, t1;
		bool_t ret;

		recip_sqrt_fp_3mod4(t0, inp);

		Fp::mul(t0, t0, inp);
		Fp::sqr(t1, t0);
		ret = vec_is_equal(t1, inp, sizeof(t1));
		vec_copy(out, t0, sizeof(t0));

		return ret;
	}

	BLST_HD BLST_FORCEINLINE int blst_fp_sqrt(Fp& out, const Fp& inp)
	{
		return (int)sqrt_fp(out, inp);
	}

	BLST_HD BLST_FORCEINLINE int blst_fp_is_square(const Fp& inp)
	{
		return (int)ct_is_square_mod_v384(inp, get_BLST_P());
	}
	BLST_HD BLST_FORCEINLINE bool_t Fp::recip_sqrt(Fp& out, const Fp& inp)
	{
		return recip_sqrt_fp(out, inp);

	}
	BLST_HD BLST_FORCEINLINE bool_t Fp::sqrt(Fp& out, const Fp& inp)
	{
		return sqrt_fp(out, inp);
	}



#if defined(__OPTIMIZE_SIZE__) || defined(__CUDA_ARCH__)
	/*
	 * 608 multiplications for scalar inversion modulo BLS12-381 prime, 32%
	 * more than corresponding optimal addition-chain, plus mispredicted
	 * branch penalties on top of that... The addition chain below was
	 * measured to be >50% faster.
	 */
	BLST_HD BLST_INLINE void flt_reciprocal_fp(Fp& out, const Fp& inp)
	{
		const byte BLS12_381_P_minus_2[] = {
			TO_BYTES(0xb9feffffffffaaa9), TO_BYTES(0x1eabfffeb153ffff),
			TO_BYTES(0x6730d2a0f6b0f624), TO_BYTES(0x64774b84f38512bf),
			TO_BYTES(0x4b1ba7b6434bacd7), TO_BYTES(0x1a0111ea397fe69a)
		};

		exp_mont_384(out, inp, BLS12_381_P_minus_2, 381);
	}
#else
#	define sqr(ret,a)		Fp::sqr(ret,a)
#	define mul(ret,a,b)		Fp::mul(ret,a,b)
#	define sqr_n_mul(ret,a,n,b)	Fp::sqr_n_mul(ret,a,n,b)

	BLST_HD BLST_INLINE void flt_reciprocal_fp(Fp& out, const Fp& inp)
	{
		RECIPROCAL_MOD_BLS12_381_P(out, inp, Fp);
	}
#	undef RECIPROCAL_MOD_BLS12_381_P
#	undef sqr_n_mul
#	undef mul
#	undef sqr
#endif

	BLST_HD BLST_FORCEINLINE void Fp::flt_reciprocal(Fp& out, const Fp& inp)
	{
		flt_reciprocal_fp(out, inp);
	}


	BLST_HD BLST_INLINE void reciprocal_fp(Fp& out, const Fp& inp)
	{
		const Fp Px8 = {    /* left-aligned value of the modulus */
			TO_LIMB_T(0xcff7fffffffd5558), TO_LIMB_T(0xf55ffff58a9ffffd),
			TO_LIMB_T(0x39869507b587b120), TO_LIMB_T(0x23ba5c279c2895fb),
			TO_LIMB_T(0x58dd3db21a5d66bb), TO_LIMB_T(0xd0088f51cbff34d2)
		};
		union { UFp x; Fp r[2]; } temp;

		ct_inverse_mod_v384(temp.x, inp, get_BLST_P(), Px8);
		Fp::redc(temp.r[0], temp.x);
		Fp::mul(temp.r[0], temp.r[0], get_BLST_RR());

#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
		/* sign goes straight to flt_reciprocal */
		Fp::mul(temp.r[1], temp.r[0], inp);
		if (vec_is_equal(&temp.r[1], &get_BLST_Rx().p, sizeof(Fp)) |
			vec_is_zero(&temp.r[1], sizeof(Fp)))
			vec_copy(&out, &temp.r[0], sizeof(Fp));
		else
			flt_reciprocal_fp(out, inp);
#else
		vec_copy(&out, &temp.r[0], sizeof(Fp));
#endif
	}

	BLST_HD BLST_FORCEINLINE void Fp::reciprocal(Fp& out, const Fp& inp) {
		reciprocal_fp(out, inp);
	}


	BLST_HD BLST_FORCEINLINE void blst_fp_inverse(Fp& out, const Fp& inp)
	{
		reciprocal_fp(out, inp);
	}

	BLST_HD BLST_FORCEINLINE void blst_fp_eucl_inverse(Fp& ret, const Fp& a)
	{
		reciprocal_fp(ret, a);
	}




	BLST_HD BLST_FORCEINLINE const Fp& Fp::zero()
	{
		return get_BLST_Rx().p2[1];
	}
	BLST_HD BLST_FORCEINLINE const Fp& Fp::one()
	{
		return get_BLST_Rx().p;
	}
	BLST_HD BLST_FORCEINLINE const Fp& Fp::modulus()
	{
		return get_BLST_P();
	}

	//BLST_EXTERN_C_END
}

#endif /* BLST_FP_H */
