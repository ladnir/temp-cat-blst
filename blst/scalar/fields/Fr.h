/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BLST_SCALAR_FR_H
#define BLST_SCALAR_FR_H

#include "blst/common/vect.h"
#include "bindings/blst_types.h"

#include "blst/backend/backend.h"
#include <bit>
#include <cstring>
#include <string.h>

namespace blst::scalar
{

	// -r^{-1} mod 2^64, used in Montgomery reduction.
#define BLST_r0 (limb_t)0xfffffffeffffffff  /* -1/r */

	/*
	 * Scalar field Fr = Z_r (r is the BLS12-381 group order).
	 *
	 * Representation:
	 *   - Values are stored in Montgomery form unless noted.
	 *     If a' is the canonical integer in [0, r), we store a = a' * R mod r,
	 *     where R = 2^256.
	 *
	 * Conventions:
	 *   - Unless stated, all ops are modulo r and return reduced results in [0, r)
	 *     and stay in Montgomery form.
	 *   - Arithmetic is intended to be constant-time w.r.t. secrets.
	 */
#define BLST_r_LIMBS \
    TO_LIMB_T(0xffffffff00000001), TO_LIMB_T(0x53bda402fffe5bfe),\
    TO_LIMB_T(0x3339d80809a1d805), TO_LIMB_T(0x73eda753299d7d48)

#define BLST_rRR_LIMBS \
    TO_LIMB_T(0xc999e990f3f29c6d), TO_LIMB_T(0x2b6cedcb87925c23),\
    TO_LIMB_T(0x05d314967254398f), TO_LIMB_T(0x0748d9d99f59ff11)

#define BLST_GT_LAMBDA_LIMBS \
    TO_LIMB_T(0x2dfefffeffff0001), TO_LIMB_T(0x53bda402fffe5bfe), TO_LIMB_T(0x3339d80809a1d805), TO_LIMB_T(0x73eda753299d7d48)

#define BLST_GT_LAMBDA2_LIMBS \
    TO_LIMB_T(0x0000000100000000), TO_LIMB_T(0xac45a4010001a402), TO_LIMB_T(0), TO_LIMB_T(0)





	/* z^4 - z^2 + 1, group order */
	constexpr Fr BLST_r_HOST = { BLST_r_LIMBS };

	/* (1<<512)%r, "radix"^2, to-Montgomery */
	constexpr Fr BLST_rRR_HOST = { BLST_rRR_LIMBS };

	constexpr Fr BLST_r_zero_HOST = { 0,0,0,0 };
	constexpr Fr BLST_r_one_HOST = { 0x1fffffffe, 0x5884b7fa00034802, 0x998c4fefecbc4ff5, 0x1824b159acc5056f };
	constexpr Fr BLST_GT_LAMBDA_HOST = { BLST_GT_LAMBDA_LIMBS };
	constexpr Fr BLST_GT_LAMBDA2_HOST = { BLST_GT_LAMBDA2_LIMBS };


#if defined(BLST_CUDA) && defined(__CUDACC__)

	static __constant__ Fr BLST_r_DEV = { BLST_r_LIMBS };
	static __constant__ Fr BLST_rRR_DEV = { BLST_rRR_LIMBS };
	static __constant__ Fr BLST_r_zero_DEV = { 0,0,0,0 };
	static __constant__ Fr BLST_r_one_DEV = { 0x1fffffffe, 0x5884b7fa00034802, 0x998c4fefecbc4ff5, 0x1824b159acc5056f };

	static __constant__ Fr BLST_GT_LAMBDA_DEV = { BLST_GT_LAMBDA_LIMBS };
	static __constant__ Fr BLST_GT_LAMBDA2_DEV = { BLST_GT_LAMBDA2_LIMBS };


#endif // __CUDA_ARCH__

	// returns the order of the group G1. This prime define the scalar field.
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE constexpr const Fr& get_BLST_r(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_r_DEV;
#else
		return BLST_r_HOST;
#endif
	}

	// returns R^2 mod r, the to-Montgomery form value for the group order r.
	// Perform montgomery multiplication with this value to convert a normal
	// integer into montgomery form.
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE const Fr& get_BLST_rRR(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_rRR_DEV;
#else
		return BLST_rRR_HOST;
#endif
	}

	// Returns lambda = p mod r (as an Fr, canonical limbs; NOT Montgomery).
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE constexpr const Fr& get_BLST_gt_lambda(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_GT_LAMBDA_DEV;
#else
		return BLST_GT_LAMBDA_HOST;
#endif
	}

	// Returns lambda^2 = p^2 mod r (as an Fr, canonical limbs; NOT Montgomery).
	// Will return host or device memory depending on compilation target.
	BLST_HD BLST_FORCEINLINE constexpr const Fr& get_BLST_gt_lambda2(void)
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_GT_LAMBDA2_DEV;
#else
		return BLST_GT_LAMBDA2_HOST;
#endif
	}


	BLST_HD BLST_FORCEINLINE limb_t Fr::check_mod(const Fr& a)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return check_mod_256(a.l, &modulus());
#else
		return cc_check_mod_n<NLIMBS(256)>((u8*)a.l, modulus().l);
#endif
	}

	BLST_HD BLST_FORCEINLINE limb_t add_n_check_mod_v256(pow256 ret, const pow256 a, const pow256 b,
		const Fr& p)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return add_n_check_mod_256(ret, a, b, &p);
#else
		return cc_add_n_check_mod_n<NLIMBS(256)>(ret, a, b, p.l);
#endif
	}

	// performs modular subtraction with zero-checking - 
	// it computes (a - b) mod p and returns whether the result is non-zero.
	BLST_HD BLST_FORCEINLINE limb_t sub_n_check_mod_v256(pow256 ret, const pow256 a, const pow256 b,
		const Fr& p)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return sub_n_check_mod_256(ret, a, b, &p);
#else
		return cc_sub_n_check_mod_n<NLIMBS(256)>(ret, a, b, p.l);
#endif
	}

	BLST_HD BLST_INLINE void cc_ct_inverse_mod_256(UFr& ret, const Fr& inp,
		const Fr& mod, const Fr& modx)
	{
		cc_ct_inverse_mod_n<NLIMBS(256)>(ret.data, inp.l, mod.l, modx.l);
	}

	BLST_HD BLST_FORCEINLINE void ct_inverse_mod_v256(UFr& ret, const Fr& inp, const Fr& mod,
		const Fr& modx)
	{
		blst::ct_inverse_mod_v256(ret.data, inp.l, mod.l, modx.l);
	}

	BLST_HD BLST_INLINE const Fr& Fr::zero()
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_r_zero_DEV;
#else
		return BLST_r_zero_HOST;
#endif
	}


	BLST_HD BLST_INLINE const Fr& Fr::one()
	{
#if defined(BLST_CUDA) && defined(__CUDA_ARCH__)
		return BLST_r_one_DEV;
#else
		return BLST_r_one_HOST;
#endif
	}


	BLST_HD BLST_INLINE void Fr::ct_bfly(Fr& x0, Fr& x1, const Fr& twiddle)
	{
		Fr x2;

		mul(x2, x1, twiddle);
		sub(x1, x0, x2);
		add(x0, x0, x2);
	}
	BLST_HD BLST_INLINE void Fr::gs_bfly(Fr& x0, Fr& x1, const Fr& twiddle)
	{
		Fr x2;
		Fr x3;

		sub(x2, x0, x1);
		add(x3, x0, x1);
		mul(x1, x2, twiddle);
		x0 = x3;
	}

	BLST_HD BLST_INLINE void Fr::from_scalar(Fr& ret, const u8 a[256 / 8])
	{
		if (static_cast<const void*>(ret.l) == static_cast<const void*>(a)) {
			Fr::to_mont(ret, ret);
		}
		else {
			Fr out;
			memcpy(&out, a, sizeof(out));
			Fr::to_mont(ret, out);
			vec_zero(&out, sizeof(out));
		}
	}

	BLST_HD BLST_INLINE void Fr::to_scaler(u8 ret[256 / 8], const Fr& a)
	{
		Fr out;
		Fr::from_mont(out, a);
		memcpy(ret, &out, sizeof(out));
		vec_zero(&out, sizeof(out));
	}


	BLST_HD BLST_INLINE void Fr::from_uint64(Fr& ret, const u64 a[256 / 64])
	{
		Fr A;
		memcpy(&A, a, sizeof(Fr));
		Fr::to_mont(ret, A);
	}

	BLST_HD BLST_INLINE void Fr::from_uint64(Fr& ret, const u64 a)
	{
		u64 aa[256 / 64];
		memset(&aa, 0, sizeof aa);
		aa[0] = a;
		from_uint64(ret, aa);
	}

	BLST_HD BLST_INLINE Fr Fr::from_uint64(const u64 a)
	{
		u64 aa[256 / 64];
		memset(&aa, 0, sizeof aa);
		aa[0] = a;
		Fr ret;
		from_uint64(ret, aa);
		return ret;
	}

	[[nodiscard]] BLST_HD BLST_INLINE bool Fr::from_hexascii(Fr& ret, const char* hex)
	{
		Fr parsed = {};
		if (!limbs_from_hexascii(parsed.l, sizeof(Fr), hex))
			return false;
		Fr::to_mont(parsed, parsed);
		ret = parsed;
		return true;
	}


	// write v to out in hex. out must at least be 65 bytes. will be 0 terminated.
	BLST_HD BLST_INLINE void Fr::to_hexascii(char* out, const Fr& v)
	{
		Fr tmp;
		Fr::from_mont(tmp, v);
		limbs_to_hexascii(out, sizeof(Fr), tmp.l);
	}

	BLST_HD BLST_FORCEINLINE int sk_check(const pow256 a)
	{
		Fr aa;
		limbs_from_le_bytes(aa.l, a, sizeof(pow256));
		const int ret = (int)Fr::check_mod(aa);
		vec_zero(&aa, sizeof(aa));
		return ret;
	}

	BLST_HD BLST_FORCEINLINE int sk_add_n_check(pow256 ret, const pow256 a, const pow256 b)
	{
		return (int)add_n_check_mod_v256(ret, a, b, get_BLST_r());
	}

	BLST_HD BLST_FORCEINLINE int sk_sub_n_check(pow256 ret, const pow256 a, const pow256 b)
	{
		return (int)sub_n_check_mod_v256(ret, a, b, get_BLST_r());
	}

	BLST_HD BLST_INLINE int sk_mul_n_check(pow256 ret, const pow256 a, const pow256 b)
	{
		Fr aa;
		Fr bb;
		Fr product;

		memcpy(&aa, a, sizeof(aa));
		memcpy(&bb, b, sizeof(bb));

		Fr::to_mont(product, aa);
		Fr::mul(product, product, bb);
		memcpy(ret, &product, sizeof(product));

		const bool_t is_zero = vec_is_zero(product.l, sizeof(product));
		vec_zero(&aa, sizeof(aa));
		vec_zero(&bb, sizeof(bb));
		vec_zero(&product, sizeof(product));

		return (int)(is_zero ^ 1);
	}

	BLST_HD BLST_INLINE void sk_inverse(pow256 ret, const pow256 a)
	{
		Fr in;
		Fr out;

		memcpy(&in, a, sizeof(in));
		Fr::to_mont(out, in);
		Fr::reciprocal(out, out);
		Fr::from_mont(out, out);
		memcpy(ret, &out, sizeof(out));

		vec_zero(&in, sizeof(in));
		vec_zero(&out, sizeof(out));
	}



	// the modulus (returned as a Fr which is a slight abuse)
	BLST_HD BLST_FORCEINLINE constexpr const Fr& Fr::modulus()
	{
		return get_BLST_r();
	}

	BLST_HD BLST_FORCEINLINE constexpr limb_t Fr::montgomery_constant()
	{
		return BLST_r0;
	}

	// ret = a + b
	BLST_HD BLST_FORCEINLINE  void Fr::add(Fr& ret, const Fr& a, const Fr& b)
	{
#if defined(__CUDA_ARCH__)
		ptx_add_mod_u64<Fr>(ret.l, a.l, b.l);
#else
		add_mod_v256(ret.l, a.l, b.l, modulus().l);
#endif
	}

	// ret = a - b
	BLST_HD BLST_FORCEINLINE  void Fr::sub(Fr& ret, const Fr& a, const Fr& b)
	{
#if defined(__CUDA_ARCH__)
		ptx_sub_mod_u64<Fr>(ret.l, a.l, b.l);
#else
		sub_mod_v256(ret.l, a.l, b.l, modulus().l);
#endif
	}

	// ret  = a * b
	BLST_HD BLST_FORCEINLINE  void Fr::mul(Fr& ret, const Fr& a, const Fr& b)
	{
#if defined(__CUDA_ARCH__)
		ptx_mul_mont_u64<Fr>(ret.l, a.l, b.l);
#else
		mul_mont_sparse_v256(ret.l, a.l, b.l, modulus().l, BLST_r0);
#endif
	}


	/*| out | = | inp | ^ | pow | , small footprint, public exponent
	*/
	BLST_HD BLST_INLINE void Fr::exp(Fr& out, const Fr& inp, const byte* pow,
		size_t pow_bits)
	{
		Fr ret;
		if (pow_bits == 0) {
			vec_copy(out.l, Fr::one().l, sizeof(out));
			return;
		}

		// Requires the most significant declared bit to be set.
		BLST_ASSERT((pow[(pow_bits - 1) >> 3] >> ((pow_bits - 1) & 7)) & 1);

		vec_copy(ret.l, inp.l, sizeof(ret));  /* ret = inp^1 */
		--pow_bits; /* most significant bit is set, skip over */
		while (pow_bits--) {
			Fr::sqr(ret, ret);
			if (is_bit_set(pow, pow_bits))
				Fr::mul(ret, ret, inp);
		}
		vec_copy(out.l, ret.l, sizeof(ret));  /* out = ret */
	}


	// ret = a^2
	BLST_HD BLST_FORCEINLINE  void Fr::sqr(Fr& ret, const Fr& a)
	{
#if defined(__CUDA_ARCH__)
		ptx_sqr_mont_u64<Fr>(ret.l, a.l);
#else
		sqr_mont_sparse_v256(ret.l, a.l, get_BLST_r().l, BLST_r0);
#endif
	}

	// Let a = a' R mod r be the Montgomery representation of a'.
	// ret = a'
	BLST_HD BLST_FORCEINLINE  void Fr::from_mont(Fr& ret, const Fr& a)
	{
		from_mont_v256(ret.l, a.l, get_BLST_r().l, BLST_r0);
	}

	// ret = a * R mod r, where R = 2^256.
	BLST_HD BLST_FORCEINLINE  void Fr::to_mont(Fr& ret, const Fr& a)
	{
		Fr::mul(ret, get_BLST_rRR(), a);
	}

	// ret = a mod r
	BLST_HD BLST_FORCEINLINE  void Fr::redc(Fr& ret, const struct UFr& a)
	{
		redc_mont_v256(ret.l, a.data, modulus().l, BLST_r0);
	}

	// ret = flag ? -a : a
	BLST_HD BLST_FORCEINLINE  void Fr::cneg(Fr& ret, const Fr& a, bool_t flag)
	{
		cneg_mod_v256(ret.l, a.l, flag, modulus().l);
	}

	// ret = -a
	BLST_HD BLST_FORCEINLINE  void Fr::neg(Fr& ret, const Fr& a)
	{
		Fr::cneg(ret, a, 1);
	}

	// ret = a * 2^count
	BLST_HD BLST_FORCEINLINE  void Fr::lshift(Fr& ret, const Fr& a, size_t count)
	{
		lshift_mod_v256(ret.l, a.l, count, get_BLST_r().l);
	}

	// ret = a / 2^count mod r
	BLST_HD BLST_FORCEINLINE  void Fr::rshift(Fr& ret, const Fr& a, size_t count)
	{
		rshift_mod_v256(ret.l, a.l, count, get_BLST_r().l);
	}

	// out = a * 3
	BLST_HD BLST_FORCEINLINE  void Fr::mul_by_3(Fr& ret, const Fr& a)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		mul_by_3_mod_256(&ret, &a, &get_BLST_r());
#else
		cc_mul_by_3_mod_n<NLIMBS(256)>(ret.l, a.l, get_BLST_r().l);
#endif
	}

	// out = 1 / inp, computed via Lehmer almost-inverse (faster)
	BLST_HD BLST_FORCEINLINE  void Fr::reciprocal(Fr& out, const Fr& inp)
	{
		const Fr rx2 = { /* left-aligned value of the modulus */
			TO_LIMB_T(0xfffffffe00000002), TO_LIMB_T(0xa77b4805fffcb7fd),
			TO_LIMB_T(0x6673b0101343b00a), TO_LIMB_T(0xe7db4ea6533afa90),
		};
		UFr temp;

		ct_inverse_mod_v256(temp, inp, get_BLST_r(), rx2);
		Fr::redc(out, temp);
		Fr::to_mont(out, out);
	}

	// return from_mont(a) == from_mont(b). 
	BLST_HD BLST_FORCEINLINE  bool_t Fr::equal(const Fr& a, const Fr& b)
	{
		return vec_is_equal(a.l, b.l, sizeof(Fr));
	}

	/**
	 * Checks if `a` is reduced modulo a given modulus.
	 *
	 * Returns 1 if `a` is less than `p`, otherwise returns 0.
	 */
	BLST_HD BLST_FORCEINLINE  limb_t Fr::is_redc(const Fr& a)
	{
		return cc_is_redc_mod_n<num_limbs()>((u8*)a.l, modulus().l);
	}


}

#endif
