

#ifndef BLST_FR_H
#define BLST_FR_H

#include "blst/common/vect.h"

#include "blst/backend/backend_asm.h"
#include "blst/backend/backend_portable.h"
#include "blst/backend/backend_cuda.h"
#include <bit>
//BLST_EXTERN_C_BEGIN

#include <iostream>

namespace blst
{

	// -1/r, this is the modular inverse of r modulo 2^256
	// used in montgomery reduction as the "magic" number
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
	struct Fr
	{
		vec_256 data;

		BLST_HD auto&& operator[](u64 idx) { return data[idx]; }
		BLST_HD auto&& operator[](u64 idx) const { return data[idx]; }


		BLST_HD static const Fr& zero();
		BLST_HD static const Fr& one();


		// the modulus (returned as a Fr which is a slight abuse)
		BLST_HD static constexpr const Fr& modulus();

		// returns the montomergy constant, -1/r, this is the modular 
		// inverse of r modulo 2^256 used in montgomery reduction as 
		// the "magic" number
		BLST_HD static constexpr limb_t montgomery_constant();

		// ret = a + b
		BLST_HD static void add(Fr& ret, const Fr& a, const Fr& b);

		// ret = a - b
		BLST_HD static void sub(Fr& ret, const Fr& a, const Fr& b);

		// ret  = a * b
		BLST_HD static void mul(Fr& ret, const Fr& a, const Fr& b);

		// ret  = a^x
		BLST_HD static void exp(Fr& ret, const Fr& a, const byte* x, size_t x_bits);

		// ret = a^2
		BLST_HD static void sqr(Fr& ret, const Fr& a);

		// let a = a' R mod p be the montgomery represtation of a'.
		// ret = a'
		BLST_HD static void from_mont(Fr& ret, const Fr& a);

		// ret = a * R, where R is the montgomery constant 2^384
		BLST_HD static void to_mont(Fr& ret, const Fr& a);

		// ret = a mod p
		BLST_HD static void redc(Fr& ret, const struct UFr& a);

		// ret = flag ? -a : a
		BLST_HD static void cneg(Fr& ret, const Fr& a, bool_t flag);

		// ret = -a
		BLST_HD static void neg(Fr& ret, const Fr& a);

		// ret = a * 2^count
		BLST_HD static void lshift(Fr& ret, const Fr& a, size_t count);

		// ret = floor(int(a) / 2^count). not mod p.
		BLST_HD static void rshift(Fr& ret, const Fr& a, size_t count);

		// out = a / 2 
		//BLST_HD static void div_by_2(Fr& ret, const Fr& a);

		// out = a * 3
		BLST_HD static void mul_by_3(Fr& ret, const Fr& a);

		// out = a * 8
		//BLST_HD static void mul_by_8(Fr& ret, const Fr& a);

		// out = 1 / inp, computed via Lehmer almost-inverse (faster)
		BLST_HD static void reciprocal(Fr& out, const Fr& inp);

		// return from_mont(a) == from_mont(b). 
		BLST_HD static bool_t equal(const Fr& a, const Fr& b);

		BLST_HD static  void ct_bfly(Fr& x0, Fr& x1, const Fr& twiddle);

		BLST_HD static void gs_bfly(Fr& x0, Fr& x1, const Fr& twiddle);

		/**
		 * Checks if `a` is reduced modulo a the modulus.
		 *
		 * Returns 1 if `a` is less than `p`, otherwise returns 0.
		 */
		BLST_HD static limb_t is_redc(const Fr& a);


		// Checks if `a` is reduced modulo a the modulus and non-zero.
		BLST_HD static limb_t check_mod(const Fr& a);

		BLST_HD static void from_scalar(Fr& ret, const u8 a[256 / 8]);

		BLST_HD static void to_scaler(u8 ret[256 / 8], const Fr& a);

		BLST_HD static void from_uint64(Fr& ret, const u64 a[256 / 64]);

		BLST_HD static void from_uint64(Fr& ret, const u64 a);

		BLST_HD static Fr from_uint64(const u64 a);

		BLST_HD static void from_hexascii(Fr& ret, const char* hex);

		// write v to out in hex. out must at least be 65 bytes. will be 0 terminated.
		BLST_HD static void to_hexascii(char* out, const Fr& v);


		BLST_HD static constexpr size_t num_limbs() { return NLIMBS(256); }


		BLST_HD BLST_FORCEINLINE bool operator==(const Fr& other) const { return equal(*this, other); }
		BLST_HD BLST_FORCEINLINE bool operator!=(const Fr& other) const { return !equal(*this, other); }
	};

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
		return check_mod_256(a.data, &modulus());
#else
		return cc_check_mod_n<NLIMBS(256)>((u8*)a.data, modulus().data);
#endif
	}

	BLST_HD BLST_FORCEINLINE limb_t add_n_check_mod_v256(pow256 ret, const pow256 a, const pow256 b,
		const Fr& p)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		return add_n_check_mod_256(ret, a, b, &p);
#else
		//return cc_add_n_check_mod_256(ret, a, b, p);
		return cc_add_n_check_mod_n<NLIMBS(256)>(ret, a, b, p.data);
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
		//return cc_sub_n_check_mod_256(ret, a, b, p);
		return cc_sub_n_check_mod_n<NLIMBS(256)>(ret, a, b, p.data);
#endif
	}

	BLST_HD BLST_INLINE void cc_ct_inverse_mod_256(UFr& ret, const Fr& inp,
		const Fr& mod, const Fr& modx)
	{
		cc_ct_inverse_mod_n<NLIMBS(256)>(ret.data, inp.data, mod.data, modx.data);
	}



	BLST_HD BLST_FORCEINLINE void ct_inverse_mod_v256(UFr& ret, const Fr& inp, const Fr& mod,
		const Fr& modx)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		ct_inverse_mod_256(&ret, &inp, &mod, &modx);
#else	
		cc_ct_inverse_mod_256(ret, inp, mod, modx);
#endif
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
	BLST_HD BLST_INLINE void print(Fr& x0)
	{
		for (size_t i = 0; i < Fr::num_limbs(); i++) {
			printf("%016llx ", (unsigned long long)x0[i]);
		}
		printf("\n");
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
		const union {
			long one;
			char little;
		} is_endian = { 1 };

		if ((void*)ret.data == (void*)a && is_endian.little) {
			Fr::to_mont((Fr&)ret, (const Fr&)a);
		}
		else {
			Fr out;
			limbs_from_le_bytes(out.data, a, 32);
			Fr::to_mont(ret, out);
			vec_zero(out.data, sizeof(out));
		}
	}

	BLST_HD BLST_INLINE void Fr::to_scaler(u8 ret[256 / 8], const Fr& a)
	{
		const union {
			long one;
			char little;
		} is_endian = { 1 };

		if ((void*)ret == (void*)a.data && is_endian.little) {
			Fr::from_mont((Fr&)ret, a);
		}
		else {
			Fr out;
			Fr::from_mont(out, a);
			le_bytes_from_limbs(ret, out.data, 32);
			vec_zero(&out, sizeof(out));
		}
	}


	BLST_HD BLST_INLINE void Fr::from_uint64(Fr& ret, const u64 a[256 / 64])
	{
		const union {
			long one;
			char little;
		} is_endian = { 1 };

		if (sizeof(limb_t) == 4 && !is_endian.little) {
			int i;
			for (i = 0; i < 4; i++) {
				unsigned long long L = a[i];
				ret[2 * i] = (limb_t)L;
				ret[2 * i + 1] = (limb_t)(L >> 32);
			}
			a = (u64*)&ret;
		}
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

	BLST_HD BLST_INLINE void Fr::from_hexascii(Fr& ret, const char* hex)
	{
		limbs_from_hexascii(ret.data, sizeof(Fr), hex);
		Fr::to_mont(ret, ret);
	}


	// write v to out in hex. out must at least be 65 bytes. will be 0 terminated.
	BLST_HD BLST_INLINE void Fr::to_hexascii(char* out, const Fr& v)
	{
		Fr tmp;
		Fr::from_mont(tmp, v);
		limbs_to_hexascii(out, sizeof(Fr), tmp.data);
	}

	BLST_HD BLST_FORCEINLINE int sk_check(const pow256 a)
	{
		Fr aa;
		limbs_from_le_bytes(aa.data, a, sizeof(pow256));
		return (int)Fr::check_mod(aa);
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
		Fr t[2];
		const union {
			long one;
			char little;
		} is_endian = { 1 };
		bool_t is_zero;

		if (((size_t)a | (size_t)b) % sizeof(limb_t) != 0 || !is_endian.little) {
			limbs_from_le_bytes(t[0].data, a, sizeof(pow256));
			limbs_from_le_bytes(t[1].data, b, sizeof(pow256));
			a = (const byte*)&t[0];
			b = (const byte*)&t[1];
		}
		Fr::to_mont(t[0], (const Fr&)a);
		Fr::mul(t[0], t[0], (const Fr&)b);
		le_bytes_from_limbs(ret, t[0].data, sizeof(pow256));
		is_zero = vec_is_zero(t[0].data, sizeof(Fr&));
		vec_zero(t, sizeof(t));

		return (int)(is_zero ^ 1);
	}

	BLST_HD BLST_INLINE void sk_inverse(pow256 ret, const pow256 a)
	{
		const union {
			long one;
			char little;
		} is_endian = { 1 };

		if (((size_t)a | (size_t)ret) % sizeof(limb_t) == 0 && is_endian.little) {
			auto& out = (Fr&)ret;
			//Fr::mul(out, (const Fr&)a, get_BLST_rRR());
			Fr::to_mont(out, (const Fr&)a);
			Fr::reciprocal(out, out);
			Fr::from_mont(out, out);
		}
		else {
			Fr out;
			limbs_from_le_bytes(out.data, a, 32);
			//Fr::mul(out, out, get_BLST_rRR());
			Fr::to_mont(out, out);
			Fr::reciprocal(out, out);
			Fr::from_mont(out, out);
			le_bytes_from_limbs(ret, out.data, 32);
			vec_zero(out.data, sizeof(out));
		}
	}


	//BLST_HD BLST_FORCEINLINE const Fr& Fr::one() {
	//	static Fr one_in_montgomery;
	//	static bool initialized = false;

	//	if (!initialized) {
	//		// Compute 1 in Montgomery form: (1 * R) mod r
	//		Fr one_normal;
	//		one_normal.data[0] = 1; // Set the first limb_t to 1 (normal integer 1)
	//		for (size_t i = 1; i < NLIMBS(256); ++i) {
	//			one_normal.data[i] = 0; // Zero out the rest of the limbs
	//		}

	//		// Convert to Montgomery form: (1 * R) mod r
	//		Fr::to_mont(one_in_montgomery, one_normal);

	//		initialized = true;
	//	}

	//	return one_in_montgomery;
	//}


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
		ptx_add_mod_u64<Fr>(ret.data, a.data, b.data);
#elif defined(BLST_DEFAULT_ASM)
		add_mod_256(&ret, &a, &b, &modulus());// linked by asm, 
#else
		cc_add_mod_n<NLIMBS(256)>(ret.data, a.data, b.data, modulus().data);
#endif
	}

	// ret = a - b
	BLST_HD BLST_FORCEINLINE  void Fr::sub(Fr& ret, const Fr& a, const Fr& b)
	{
#if defined(__CUDA_ARCH__)
		ptx_sub_mod_u64<Fr>(ret.data, a.data, b.data);
#elif defined(BLST_DEFAULT_ASM)
		sub_mod_256(&ret, &a, &b, &modulus());// linked by asm
#else
		cc_sub_mod_n<NLIMBS(256)>(ret.data, a.data, b.data, modulus().data);
#endif
	}

	// ret  = a * b
	BLST_HD BLST_FORCEINLINE  void Fr::mul(Fr& ret, const Fr& a, const Fr& b)
	{
#if defined(__CUDA_ARCH__)
		ptx_mul_mont_u64<Fr>(ret.data, a.data, b.data);
#elif defined(BLST_DEFAULT_ASM) 
		mulx_mont_sparse_256(&ret, &a, &b, &modulus(), BLST_r0);// linked by asm
#else
		cc_mul_mont_n<NLIMBS(256)>(ret.data, a.data, b.data, modulus().data, BLST_r0);
#endif
	}


	/*| out | = | inp | ^ | pow | , small footprint, public exponent
	*/
	BLST_HD BLST_INLINE void Fr::exp(Fr& out, const Fr& inp, const byte* pow,
		size_t pow_bits)
	{
		Fr ret;

		// requies MSB of pow is set
		BLST_ASSERT(pow_bits && ((pow[(pow_bits - 1) >> 3] >> ((pow_bits - 1) & 7)) & 1));

		vec_copy(ret.data, inp.data, sizeof(ret));  /* ret = inp^1 */
		--pow_bits; /* most significant bit is set, skip over */
		while (pow_bits--) {
			Fr::sqr(ret, ret);
			if (is_bit_set(pow, pow_bits))
				Fr::mul(ret, ret, inp);
		}
		vec_copy(out.data, ret.data, sizeof(ret));  /* out = ret */
	}


	// ret = a^2
	BLST_HD BLST_FORCEINLINE  void Fr::sqr(Fr& ret, const Fr& a)
	{
#if defined(__CUDA_ARCH__)
		ptx_sqr_mont_u64<Fr>(ret.data, a.data);
#elif defined(BLST_DEFAULT_ASM)
		sqrx_mont_sparse_256(&ret, &a, &get_BLST_r(), BLST_r0);// linked by asm
#else
		cc_mul_mont_n<NLIMBS(256)>(ret.data, a.data, a.data, get_BLST_r().data, BLST_r0);
#endif
	}

	// let a = a' R mod p be the montgomery represtation of a'.
	// ret = a'
	BLST_HD BLST_FORCEINLINE  void Fr::from_mont(Fr& ret, const Fr& a)
	{
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		fromx_mont_256(&ret, &a, &get_BLST_r(), BLST_r0);// linked by asm
#else
		cc_from_mont_n<NLIMBS(256)>(ret.data, a.data, get_BLST_r().data, BLST_r0);
#endif
	}

	// ret = a * R, where R is the montgomery constant 2^384
	BLST_HD BLST_FORCEINLINE  void Fr::to_mont(Fr& ret, const Fr& a)
	{
		Fr::mul(ret, a, get_BLST_rRR());
	}

	// ret = a mod p
	BLST_HD BLST_FORCEINLINE  void Fr::redc(Fr& ret, const struct UFr& a)
	{
		//redc_mont_v256(ret, a, get_BLST_r(), BLST_r0);
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		redcx_mont_256(&ret, &a, &modulus(), BLST_r0);// linked by asm
#else
		//cc_redc_mont_256(ret, a, p, n0);
		cc_redc_mont_n<NLIMBS(256)>(ret.data, a.data, modulus().data, BLST_r0);
#endif
	}

	// ret = flag ? -a : a
	BLST_HD BLST_FORCEINLINE  void Fr::cneg(Fr& ret, const Fr& a, bool_t flag)
	{
		//cneg_mod_v256(ret, a, flag, get_BLST_r());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		cneg_mod_256(&ret, &a, flag, &modulus());
#else
		//cc_cneg_mod_256(ret, a, flag, p);
		cc_cneg_mod_n<NLIMBS(256)>(ret.data, a.data, flag, modulus().data);
#endif // BLST_DEFAULT_ASM

	}

	// ret = -a
	BLST_HD BLST_FORCEINLINE  void Fr::neg(Fr& ret, const Fr& a)
	{
		Fr::cneg(ret, a, 1);
	}

	// ret = a * 2^count
	BLST_HD BLST_FORCEINLINE  void Fr::lshift(Fr& ret, const Fr& a, size_t count)
	{
		//lshift_mod_v256(ret, a, count, get_BLST_r());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		lshift_mod_256(&ret, &a, count, &get_BLST_r());
#else
		//cc_lshift_mod_256(ret, a, count, p);
		cc_lshift_mod_n<NLIMBS(256)>(ret.data, a.data, count, get_BLST_r().data);
#endif
	}

	// ret = floor(int(a) / 2^count). not mod p.
	BLST_HD BLST_FORCEINLINE  void Fr::rshift(Fr& ret, const Fr& a, size_t count)
	{
		//rshift_mod_v256(ret, a, count, get_BLST_r());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		rshift_mod_256(&ret, &a, count, &get_BLST_r());
#else
		//cc_rshift_mod_256(ret, a, count, p);
		cc_rshift_mod_n<NLIMBS(256)>(ret.data, a.data, count, get_BLST_r().data);
#endif
	}

	// out = a / 2 
	//BLST_HD BLST_FORCEINLINE  void Fr::div_by_2(Fr& ret, const Fr& a);

	// out = a * 3
	BLST_HD BLST_FORCEINLINE  void Fr::mul_by_3(Fr& ret, const Fr& a)
	{
		//mul_by_3_mod_v256(ret, a, get_BLST_r());
#if defined(BLST_DEFAULT_ASM) && !defined(__CUDA_ARCH__)
		mul_by_3_mod_256(&ret, &a, &get_BLST_r());
#else
		//cc_mul_by_3_mod_256(ret, a, p);
		cc_mul_by_3_mod_n<NLIMBS(256)>(ret.data, a.data, get_BLST_r().data);
#endif
	}

	// out = a * 8
	//BLST_HD BLST_FORCEINLINE  void Fr::mul_by_8(Fr& ret, const Fr& a);

	// out = 1 / inp, computed via Lehmer almost-inverse (faster)
	BLST_HD BLST_FORCEINLINE  void Fr::reciprocal(Fr& out, const Fr& inp)
	{
		//reciprocal_fr(out, inp);

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
		return vec_is_equal(a.data, b.data, sizeof(Fr));
	}

	/**
	 * Checks if `a` is reduced modulo a given modulus.
	 *
	 * Returns 1 if `a` is less than `p`, otherwise returns 0.
	 */
	BLST_HD BLST_FORCEINLINE  limb_t Fr::is_redc(const Fr& a)
	{
		return cc_is_redc_mod_n<num_limbs()>((u8*)a.data, modulus().data);
	}


	//BLST_EXTERN_C_END
}

#endif